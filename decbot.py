'''
This program automates our DECam observing for DECaLS.

It takes a set of three observing plans, one for each "pass" --
good/fair/poor conditions on different tilings.  As new images appear
on disk, it performs a quick-reduction to estimate the conditions and
uses that to update the exposure time of subsequent observations so
that we hit our depth goals.

We have API access to the DECam observing queue, so the bot tries to
keep some number (1 or 2) of exposures in the queue at all times.  It
polls the queue, and when there is room, it submits the next planned
exposure.
'''

from __future__ import print_function
import sys
import time
import json
import datetime
import stat
import os
from collections import OrderedDict

import numpy as np
import ephem
import fitsio
from astrometry.util.fits import fits_table

from astrometry.util.starutil_numpy import  ra2hmsstring as  ra2hms
from astrometry.util.starutil_numpy import dec2dmsstring as dec2dms

from measure_raw import measure_raw
from obsbot import (
    exposure_factor, get_tile_from_name, get_airmass,
    Obsbot, datenow, unixtime_to_ephem_date,
    ephem_date_to_mjd, choose_pass, get_forced_pass)

def main(cmdlineargs=None, get_decbot=False):
    import optparse
    parser = optparse.OptionParser(usage='%prog [<pass1.json> <pass2.json> <pass3.json>]')

    from camera_decam import (ephem_observer, default_extension, nominal_cal,
                              tile_path)
    
    parser.add_option('--rawdata', help='Directory to monitor for new images: default $DECAM_DATA if set, else "rawdata"', default=None)
    parser.add_option('--ext', help='Extension to read for computing observing conditions, default %default.  Can give comma-separated list.', default=default_extension)
    parser.add_option('--tiles', default=tile_path,
                      help='Observation status file, default %default')

    parser.add_option('--pass', dest='passnum', type=int, default=2,
                      help='Set default pass number (1/2/3), default 2')
    parser.add_option('--exptime', type=int, default=None,
                      help='Set default exposure time, default whatever is in the JSON files, usually 80 sec')

    parser.add_option('--no-adjust', dest='adjust', default=True,
                      action='store_false',
                      help='Do not adjust exposure time to compensate for previous tiles')

    parser.add_option('--nqueued', type=int, default=2,
                      help='Set maximum number of exposures in the queue; default %default')
    
    parser.add_option('--no-cut-past', dest='cut_before_now',
                      default=True, action='store_false',
                      help='Do not cut tiles that were supposed to be observed in the past')

    parser.add_option('--no-queue', dest='do_queue', default=True,
                      action='store_false',
                      help='Do not actually queue exposures.')
    parser.add_option('--remote-server', default=None,
                      help='Hostname of CommandServer for queue control')
    parser.add_option('--remote-port', default=None, type=int,
                      help='Port number of CommandServer for queue control')

    parser.add_option('--threads', type=int, default=1,
                      help='Run multi-threaded when processing multiple extensions?')
    
    parser.add_option('--verbose', default=False, action='store_true',
                      help='Turn on (even more) verbose logging')
    
    if cmdlineargs is None:
        opt,args = parser.parse_args()
    else:
        opt,args = parser.parse_args(cmdlineargs)

    if len(args) == 0:
        # Default JSON filenames
        args = ['pass1.json', 'pass2.json', 'pass3.json']

    if len(args) != 3:
        parser.print_help()
        sys.exit(-1)

    if not (opt.passnum in [1,2,3]):
        parser.print_help()
        sys.exit(-1)
        
    json1fn,json2fn,json3fn = args

    J1 = json.loads(open(json1fn,'rb').read())
    J2 = json.loads(open(json2fn,'rb').read())
    J3 = json.loads(open(json3fn,'rb').read())
    print('Read', len(J1), 'pass 1 and', len(J2), 'pass 2 and', len(J3), 'pass 3 exposures')
    if opt.cut_before_now:
        # Drop exposures that are before *now*, in all three plans.
        now = ephem.now()
        print('Now:', str(now))
        J1 = [j for j in J1 if ephem.Date(str(j['approx_datetime'])) > now]
        J2 = [j for j in J2 if ephem.Date(str(j['approx_datetime'])) > now]
        J3 = [j for j in J3 if ephem.Date(str(j['approx_datetime'])) > now]
        for i,J in enumerate([J1,J2,J3]):
            print('Pass %i: keeping %i tiles after now' % (i+1, len(J)))
            if len(J):
                print('First tile: %s' % J[0]['approx_datetime'])

    if len(J1) + len(J2) + len(J3) == 0:
        print('No tiles!')
        return

    # Annotate with 'planpass' field
    for i,J in enumerate([J1,J2,J3]):
        for j in J:
            j['planpass'] = i+1
    
    obs = ephem_observer()
    
    print('Reading tiles table', opt.tiles)
    tiles = fits_table(opt.tiles)
    
    if opt.rawdata is None:
        opt.rawdata = os.environ.get('DECAM_DATA', 'rawdata')

    if opt.do_queue:
        kw = dict()
        if opt.remote_server is not None:
            kw['cs_host'] = opt.remote_server
        if opt.remote_port is not None:
            kw['cs_port'] = opt.remote_port
        from RemoteClient import RemoteClient
        print('Creating RemoteClient connection with args:', kw)
        rc = RemoteClient(**kw)
    else:
        rc = None

    # Try loading copilot's database
    try:
        import obsdb
        from camera_decam import database_filename
        obsdb.django_setup(database_filename=database_filename)
        copilot_db = obsdb.MeasuredCCD.objects
    except:
        copilot_db = None
    
    decbot = Decbot(J1, J2, J3, opt, nominal_cal, obs, tiles, rc,
                    copilot_db=copilot_db, nqueued=opt.nqueued)
    if get_decbot:
        return decbot
    decbot.queue_initial_exposures()
    decbot.run()


class Decbot(Obsbot):
    '''
    How the control flow works and when stuff happens:

    Startup:
    - Decbot.queue_initial_exposures() -> Decbot.heartbeat()
    - Obsbot.run():
        while True:
            sleep(5)
            Obsbot.run_one()
                - look for a new file, try to open it
                - Decbot.process_file()
                    - measure image, set self.latest_measurement
                    - Decbot.update_plans()
            Decbot.heartbeat()
                - poll queue.  If not empty:
                    - Decbot.queue_exposure()
                        - save in queued_tiles
                - else if time after time of first exposure in any pass:
                    - Decbot.update_plans()

    Decbot.update_plans():
        - choose next pass
        - for each pass:
            - drop tiles before now
            - drop any tiles already observed or planned
            - update exposure times
        - Decbot.write_plans()


    Data model:

    - J1, J2, J3 are the plans for passes 1,2,3 (the name "J" comes
      from the JSON files these are read from).  We drop exposures
      from these lists using the keep_good_tiles() function.
      Exposures (aka tiles in the code) can be dropped if they were
      planned for before now, or if we have already queued that tile
      in this run of decbot.py, or if we have seen that object name in
      a file on disk, or if it has appeared previously in this pass's
      plan.

    - nextpass -- is which pass we currently think we should be
      running.  Therefore the next exposure to be queued will be
      ([J1,J2,J3][nextpass-1])[0]

    - latest_measurement -- our most recent measurement of the
      conditions.  Updated when we see a new file on disk.

    '''

    def __init__(self, J1, J2, J3, opt, nom, obs, tiles, rc,
                 nqueued=2,
                 copilot_db=None,):
        super(Decbot, self).__init__(
            opt.rawdata, backlog=False, only_process_newest=True,
            ignore_missing_dir=True, verbose=opt.verbose)
        self.timeout = None
        self.nqueued = nqueued
        self.J1 = J1
        self.J2 = J2
        self.J3 = J3
        self.opt = opt
        self.nom = nom
        self.obs = obs
        # tiles = the decam-tiles_obstatus.fits table
        self.tiles = tiles
        # rc = API access to the queue -- SISPI RemoteClient
        self.rc = rc
        self.copilot_db = copilot_db
        self.queued_tiles = []
        self.adjust_previous = opt.adjust
        self.nextpass = opt.passnum
        self.latest_measurement = None
        
        # Read existing files,recording their tile names as vetoes.
        self.observed_tiles = {}

        # Only check files timestamped since sunset; starting up was taking
        # too long since all files from a run go into one monster directory.
        sun = ephem.Sun()
        sunset = obs.previous_setting(sun)
        fns = []
        for fn in self.oldfiles:
            st = os.stat(fn)
            dd = unixtime_to_ephem_date(st.st_mtime)
            if dd > sunset:
                fns.append(fn)
                print('Checking %i of %i existing files with timestamps after sunset' %
              (len(fns), len(self.oldfiles)))

        if copilot_db is not None:
            # Further filter existing files...
            try:
                # Grab all copilot database entries from tonight
                sunset_mjd = ephem_date_to_mjd(sunset)
                ccds = copilot_db.filter(mjd_obs__gte=sunset_mjd)
                ccds = ccds.values_list('filename', 'object')
                # Build a map from filename (no directory path) -> obj name
                dbobjs = {}
                for fn,obj in ccds:
                    if len(obj) == 0:
                        continue
                    fn = str(fn).strip()
                    fn = os.path.basename(fn)
                    dbobjs[fn] = str(obj)
                print('Checking copilot database -- found', len(dbobjs),
                      'from tonight')
                # Go through existing files looking for database entries.
                keepfns = []
                for path in fns:
                    # check just the filename part, not the full path
                    fn = os.path.basename(path)
                    if fn in dbobjs:
                        obj = dbobjs[fn]
                        print('Found filename', fn, 'in copilot database: tile', obj)
                        self.observed_tiles[obj] = path
                    else:
                        keepfns.append(path)
                # Check all the ones not in the db
                fns = keepfns
            except:
                import traceback
                traceback.print_exc()
            
        self.new_observed_tiles(fns)

        # Set up initial planned_tiles
        self.nextpass = opt.passnum
        self.update_plans(exptime=opt.exptime)

    def choose_next_pass(self):
        M = self.latest_measurement
        if M is not None:
            trans  = M['transparency']
            seeing = M['seeing']
            skybright = M['skybright']
            # eg, nominal = 20, sky = 19, brighter is 1 mag brighter than nom.
            nomsky = self.nom.sky(M['band'])
            self.nextpass = choose_pass(trans, seeing, skybright, nomsky)
        return self.nextpass

    def try_open_file(self, path):
        ext = self.opt.ext
        # multiple extensions?
        exts = []
        if ext is not None:
            exts = ext.split(',')
        F = fitsio.FITS(path)
        for ext in exts:
            info = F[ext].get_info()
            self.debug('Checking file', path, ': ext', ext, ':', info)

    def queue_initial_exposures(self):
        # Queue exposures to start
        for i in range(self.nqueued):
            self.heartbeat()

    def saw_new_files(self, fns):
        # Update the veto list of tiles we have taken tonight.
        self.new_observed_tiles(fns)

    def new_observed_tiles(self, fns):
        ''' Reads the given list of files, looking for "OBJECT" keywords;
        adds them to observed_tiles.'''
        for fn in fns:
            try:
                phdr = fitsio.read_header(fn)
                obj = phdr['OBJECT']
                obj = str(obj).strip()
                print('Old file', fn, 'observed object', obj)
                self.observed_tiles[obj] = fn
            except:
                import traceback
                print('Failed to read header of file', fn)
                traceback.print_exc()

    def queue_if_ready(self):
        if self.rc is None:
            return False
        # Poll the queue to see if we can queue another one.
        nq = self.rc.get_n_queued()
        if nq >= self.nqueued:
            self.log('%i exposures in the queue, waiting until fewer than %i.'%
                     (nq, self.nqueued), uniq=True)
            return False
        self.log('%i exposures in the queue, time to queue one.' % (nq))
        e = self.queue_exposure()
        return (e is not None)

    def heartbeat(self):
        if self.queue_if_ready():
            return

        # Is "now" after the next tile in any pass?  If so, replan!
        replan = False
        for J in [self.J1, self.J2, self.J3]:
            if len(J):
                # Is the first planned tile scheduled for before now?
                jj = [J[0]]
                jj = self.tiles_after_now(jj)
                if len(jj) == 0:
                    replan = True
        if replan:
            self.update_plans()

    def queue_exposure(self):
        self.choose_next_pass()
        J = self.get_upcoming()
        if len(J) == 0:
            print('Time to queue an exposure, but none are left in the plans!')
            return None
        j = J.pop(0)

        if self.rc is None:
            print('Not actually queuing exposure (--no-queue):', j)
        else:
            print('Queuing exposure:', j)
            self.rc.addexposure(filter=j['filter'], ra=j['RA'], dec=j['dec'],
                                object=j['object'], exptime=j['expTime'],
                                verbose=self.verbose)
        self.queued_tiles.append(j)
        return j
    
    def filter_new_files(self, fns):
        return [fn for fn in fns if
                fn.endswith('.fits.fz') or fn.endswith('.fits')]

    def check_header(self, fn):
        # Read primary FITS header
        phdr = fitsio.read_header(fn)
        obstype = phdr.get('OBSTYPE','').strip()
        #print('Obstype:', obstype)
        if obstype in ['zero', 'focus', 'dome flat']:
            print('Skipping obstype =', obstype)
            return False
        elif obstype == '':
            print('Empty OBSTYPE in header:', fn)
            return False
    
        exptime = phdr.get('EXPTIME')
        if exptime == 0:
            print('Exposure time EXPTIME in header =', exptime)
            return False
    
        filt = str(phdr['FILTER'])
        filt = filt.strip()
        filt = filt.split()[0]
        ## DECam?
        if filt == 'solid':
            print('Solid (block) filter.')
            return False

        obj = phdr.get('OBJECT', '')
        print('Object:', obj, 'exptime', exptime, 'filter', filt)
        return True

    def measure_extensions(self, fn, ext, kwa):
        # If we're checking multiple extensions, build argument lists for each.
        args = []
        if ext is not None:
            exts = ext.split(',')
            for ext in exts:
                thiskwa = kwa.copy()
                thiskwa.update(ext=ext)
                args.append((fn, thiskwa))
        else:
            args.append((fn, kwa))
        # Measure extensions in parallel.
        from astrometry.util.multiproc import multiproc
        mp = multiproc(self.opt.threads)
        MM = mp.map(bounce_measure_raw, args)
        mp.close()
        del mp
        return MM

    def check_measurements(self, MM, fn):
        # Reasonableness checks
        keep = []
        for M in MM:
            ok = (M is not None) and (M['nmatched'] >= 20) and (M.get('zp',None) is not None)
            if ok:
                keep.append(M)
        if len(keep) == 0:
            print('Failed checks in our measurement of', fn,
                  '-- not updating anything')
            # FIXME -- we could fall back to pass 3 here.
            return []
        return keep

    def average_measurements(self, MM):
        # Average the measurements -- but start by copying one of
        # the measurements.
        M = MM[0].copy()
        # FIXME -- means?
        for key,nice in [('skybright','Sky'),
                         ('transparency','Transparency'),
                         ('seeing','Seeing')]:
            print('Measurements of %-12s:' % nice,
                  ', '.join(['%6.3f' % mi[key] for mi in MM]))
            M[key] = np.mean([mi[key] for mi in MM])
        return M

    def process_file(self, fn):
        ok = self.check_header(fn)
        if not ok:
            return False
        # Measure the new image
        kwa = dict(verbose=self.verbose, ps=None)
        MM = self.measure_extensions(fn, self.opt.ext, kwa)

        MM = self.check_measurements(MM, fn)
        if len(MM) == 0:
            return False
        
        if len(MM) == 1:
            M = MM[0]
        else:
            M = self.average_measurements(MM)

        trans  = M['transparency']
        seeing = M['seeing']
        skybright = M['skybright']
        
        # eg, nominal = 20, sky = 19, brighter is 1 mag brighter than nom.
        nomsky = self.nom.sky(M['band'])
        brighter = nomsky - skybright
    
        print('Transparency: %6.02f' % trans)
        print('Seeing      : %6.02f' % seeing)
        print('Sky         : %6.02f' % skybright)
        print('Nominal sky : %6.02f' % nomsky)
        print('Sky over nom: %6.02f   (positive means brighter than nom)' %
              brighter)

        if self.copilot_db is not None and (M['band'] in ['g','r']):
            self.recent_gr(M)
        
        self.latest_measurement = M
        self.update_plans()

    def recent_gr(self, M):
        '''
        Add to the given measurement dictionary *M* estimates of
        g and r sky and seeing based on the recent past.
        Updates *M* in-place.
        '''
        from copilot import recent_gr_sky_color, recent_gr_seeing

        gr, ndiff, ng, nr = recent_gr_sky_color()
        if gr is not None:
            M['grsky'] = gr

        gr = recent_gr_seeing()
        if gr is not None:
            M['grsee'] = gr
        
    def get_upcoming(self):
        return [self.J1,self.J2,self.J3][self.nextpass-1]

    def update_plans(self, exptime=None):
        # Sets self.nextpass
        self.choose_next_pass()

        self.J1 = self.keep_good_tiles(self.J1)
        self.J2 = self.keep_good_tiles(self.J2)
        self.J3 = self.keep_good_tiles(self.J3)

        # Update plans (exposure times) for all three passes
        for j,J in enumerate([self.J1, self.J2, self.J3]):
            passnum = j+1
            self.update_plans_for_pass(passnum, J, exptime)
        self.write_plans()

    def update_plans_for_pass(self, passnum, J, exptime):
        if len(J) == 0:
            print('Could not find a JSON observation in pass', passnum,
                  'with approx_datetime after now =', str(ephem.now()))
            return

        M = self.latest_measurement
        if M is not None:
            # Update the exposure times in plan J based on measured conditions.
            print('Updating exposure times for pass', passnum)
            # Keep track of expected time of observations
            # FIXME -- should add margin for the images currently in the queue.
            self.obs.date = ephem.now()
            for ii,jplan in enumerate(J):
                exptime = self.exptime_for_tile(jplan, debug=(ii < 3))
                jplan['expTime'] = exptime
                self.obs.date += (exptime + self.nom.overhead) / 86400.
        elif exptime is not None: 
            for ii,jplan in enumerate(J):
                jplan['expTime'] = exptime

        self.obs.date = ephem.now()
        for ii,jplan in enumerate(J):
            s = (('Plan tile %s (pass %i), band %s, RA,Dec (%.3f,%.3f), ' +
                  'exptime %i.') %
                  (jplan['object'], jplan['planpass'], jplan['filter'],
                   jplan['RA'], jplan['dec'], jplan['expTime']))
            if ii <= 3:
                airmass = self.airmass_for_tile(jplan)
                s += '  Airmass if observed now: %.2f' % airmass
                self.log(s)
            else:
                self.debug(s)
            
    def airmass_for_tile(self, jplan):
        '''
        Note, uses self.obs
        '''
        rastr  = ra2hms (jplan['RA' ])
        decstr = dec2dms(jplan['dec'])
        ephemstr = str('%s,f,%s,%s,20' % ('X', rastr, decstr))
        etile = ephem.readdb(ephemstr)
        etile.compute(self.obs)
        airmass = get_airmass(float(etile.alt))
        return airmass
    
    def exptime_for_tile(self, jplan, debug=False):
        '''
        Note, uses self.obs(.date) to compute airmass.
        '''
        tilename = str(jplan['object'])
        # Find this tile in the tiles table.
        tile = get_tile_from_name(tilename, self.tiles)
        ebv = tile.ebv_med
        band = str(jplan['filter'])[0]
        airmass = self.airmass_for_tile(jplan)
        
        M = self.latest_measurement
        trans  = M['transparency']
        seeing = M['seeing']
        msky  = M['skybright']
        grsky = M.get('grsky', None)
        grsee = M.get('grsee', None)
        mband = M['band']
        mairmass = M['airmass']
        assert(mairmass >= 1.0 and mairmass < 4.0)

        sky = self.predict_sky(mband, msky, band, grsky, debug)
        seeing = self.predict_seeing(band,seeing,mairmass,airmass, grsee, debug)

        fid = self.nom.fiducial_exptime(band)
        expfactor = exposure_factor(fid, self.nom, airmass,ebv,seeing,sky,trans)

        if self.adjust_previous:
            from astrometry.util.ttime import Time
            t0 = Time()
            adjfactor = self.adjust_for_previous(tile, band, fid)
            print('Adjust-for-previous:', Time()-t0)
            expfactor *= adjfactor

        exptime = expfactor * fid.exptime

        ### HACK -- safety factor!
        exptime *= 1.1

        exptime = int(np.ceil(exptime))
        exptime = np.clip(exptime, fid.exptime_min, fid.exptime_max)
        if band == 'z':
            # Compute cap on exposure time to avoid saturation /
            # loss of dynamic range.
            t_sat = self.nom.saturation_time(band, sky)
            if exptime > t_sat:
                exptime = t_sat
                # Don't print this a gajillion times
                self.log('Reduced exposure time to avoid z-band ' +
                         'saturation: %.1f s' % exptime, uniq=True)
        exptime = int(exptime)
        return exptime

    def predict_sky(self, mband, msky, band, grsky, debug):
        '''
        mband: measured band
        msky: measured sky in that band
        band: band that we want to predict the sky level in
        grsky: estimate of g-r sky color.
        '''
        if mband == band:
            sky = msky
        else:
            sky = None
            if ((grsky is not None) and
                ((mband == 'g' and band == 'r') or
                 (mband == 'r' and band == 'g'))):
                if debug:
                    print('g-r color:', grsky)
                    print('measured sky in', mband, '=', msky)
                if band == 'r':
                    sky = msky - grsky
                else:
                    sky = msky + grsky
                    if debug:
                        print('predicted sky in', band, '=', sky)
                
            if sky is None:
                # Guess that the sky is as much brighter than canonical
                # in the next band as it is in this one!
                nomsky = self.nom.sky(mband)
                sky = ((msky - nomsky) + self.nom.sky(band))
        return sky

    def predict_seeing(self, band, seeing, mairmass, airmass, grsee, debug):
        '''
        band: band we want to predict the seeing in.
        seeing: measured seeing
        mairmass: airmass of measured images
        airmass: airmass of the image we want to predict for.
        grsee: recent seeing estimates for g,r bands.
        '''
        oldsee = seeing
        if (grsee is not None) and (band in 'gr'):
            g_see,r_see,G,R = grsee
            oldair = mairmass
            if band == 'r':
                seeing = r_see
                mairmass = np.median(R.airmass)
            else:
                seeing = g_see
                mairmass = np.median(G.airmass)
            if debug:
                print('Using g,r seeing estimate', seeing, 'rather than most',
                      'recent measurement', oldsee)
                print('Updating airmass from', oldair, 'to', mairmass,
                      'used for seeing estimate')

        seeing_wrt_airmass = self.nom.seeing_wrt_airmass(band)
        seeing += (airmass - mairmass) * seeing_wrt_airmass
        if debug:
            print('Updated seeing prediction from', oldsee, 'to', seeing,
                  'for airmass %.2f to %.2f' % (mairmass, airmass))
        return seeing

    def tiles_after_now(self, J):
        now = ephem.now()
        keep = []
        for i,j in enumerate(J):
            tstart = ephem.Date(str(j['approx_datetime']))
            if tstart > now:
                #print('Found tile %s which starts at %s' %
                #      (j['object'], str(tstart)))
                keep.append(j)
        return keep
    
    def write_plans(self):
        # Write upcoming plan to a JSON file
        fn = 'decbot-plan.json'
        tmpfn = fn + '.tmp'
        upcoming = self.get_upcoming()
        jstr = json.dumps(upcoming, sort_keys=True,
                          indent=4, separators=(',', ': '))
        f = open(tmpfn, 'w')
        f.write(jstr + '\n')
        f.close()
        os.rename(tmpfn, fn)
        self.debug('Wrote', fn)

        fn = 'decbot-plan-5.json'
        tmpfn = fn + '.tmp'
        jstr = json.dumps(upcoming[:5], sort_keys=True,
                          indent=4, separators=(',', ': '))
        f = open(tmpfn, 'w')
        f.write(jstr + '\n')
        f.close()
        os.rename(tmpfn, fn)
        self.debug('Wrote', fn)

        # Write a FITS table of the exposures we think we've queued,
        # the ones we have planned, and the future tiles in passes 1,2,3.
        P = ([(j,'Q') for j in self.queued_tiles] +
             [(j,'P') for j in upcoming])

        for i,J in enumerate([self.J1,self.J2,self.J3]):
            passnum = i+1
            P.append((j,'%i' % passnum))
             
        J = [j for j,t in P]
        T = fits_table()
        T.type = np.array([t for j,t in P])
        T.tilename = np.array([str(j['object']) for j in J])
        T.filter = np.array([str(j['filter'])[0] for j in J])
        T.exptime = np.array([j['expTime'] for j in J])
        T.ra  = np.array([j['RA'] for j in J])
        T.dec = np.array([j['dec'] for j in J])
        T.planpass = np.array([j['planpass'] for j in J])
        fn = 'decbot-plan.fits'
        tmpfn = fn + '.tmp'
        T.writeto(tmpfn)
        os.rename(tmpfn, fn)
        self.debug('Wrote', fn)

    def keep_good_tiles(self, J):
        keep = []
        now = ephem.now()
        for j in J:
            tstart = ephem.Date(str(j['approx_datetime']))
            if tstart < now:
                continue
            tilename = str(j['object'])
            # Was this tile seen in a file on disk? (not incl. backlog)
            if tilename in self.observed_tiles:
                continue
            if object_name_in_list(j, self.queued_tiles):
                continue
            # Our plan files should already have this property..!
            if object_name_in_list(j, keep):
                continue
            keep.append(j)
        return keep

def object_name_in_list(j, Jlist):
    tilename = str(j['object'])
    for tile in Jlist:
        if str(tile['object']) == tilename:
            return True
    return False
        

def bounce_measure_raw(args):
    (fn, kwargs) = args
    try:
        return measure_raw(fn, **kwargs)
    except:
        print('Failed to measure image:', fn, kwargs)
        import traceback
        traceback.print_exc()
    return None

if __name__ == '__main__':
    main()
    
