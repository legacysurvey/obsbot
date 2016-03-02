'''
        # - queue new exposure from planned_tiles when we think the
        #   current exposure is nearly done; update seqnum
        # - update planned_tiles as new images come in
        #   - write out a JSON plan with the upcoming tiles, as a backup.
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
    NewFileWatcher, datenow, unixtime_to_ephem_date,
    ephem_date_to_mjd, choose_pass)

def main(cmdlineargs=None, get_decbot=False):
    import optparse
    parser = optparse.OptionParser(usage='%prog <pass1.json> <pass2.json> <pass3.json>')

    from camera_decam import (ephem_observer, default_extension, nominal_cal,
                              tile_path)
    
    parser.add_option('--rawdata', help='Directory to monitor for new images: default $DECAM_DATA if set, else "rawdata"', default=None)
    parser.add_option('--ext', help='Extension to read for computing observing conditions, default %default', default=default_extension)
    parser.add_option('--tiles', default=tile_path,
                      help='Observation status file, default %default')

    parser.add_option('--pass', dest='passnum', type=int, default=2,
                      help='Set default pass number (1/2/3), default 2')
    parser.add_option('--exptime', type=int, default=None,
                      help='Set default exposure time, default whatever is in the JSON files, usually 80 sec')

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
    
    if cmdlineargs is None:
        opt,args = parser.parse_args()
    else:
        opt,args = parser.parse_args(cmdlineargs)
    
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
        from camera import database_filename
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


class Decbot(NewFileWatcher):
    def __init__(self, J1, J2, J3, opt, nom, obs, tiles, rc,
                 nqueued=2,
                 copilot_db=None,):
        super(Decbot, self).__init__(opt.rawdata, backlog=False,
                                     only_process_newest=True)
        self.timeout = None
        self.nqueued = nqueued
        self.J1 = J1
        self.J2 = J2
        self.J3 = J3
        self.opt = opt
        self.nom = nom
        self.obs = obs
        self.tiles = tiles
        self.rc = rc
        self.copilot_db = copilot_db
        self.seqnum = 0
        self.planned_tiles = OrderedDict()
        self.upcoming = []

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
                print('Found', len(dbobjs), 'entries in the copilot database')
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
        J = [J1,J2,J3][opt.passnum - 1]

        J = self.tiles_after_now(J)
        self.plan_tiles(J, Nahead=len(J), exptime=opt.exptime)

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
    
    def heartbeat(self):
        if self.rc is None:
            return
        # Poll the queue to see if we can queue another one.
        nq = self.rc.get_n_queued()
        now = datenow().replace(microsecond=0).isoformat()
        if nq >= self.nqueued:
            print('%s: %i exposures in the queue, waiting until fewer than %i.' % (now, nq, self.nqueued))
            return
        else:
            print('%s: %i exposures in the queue, time to queue one.' % (now, nq))

        e = self.queue_exposure()
        if e is None:
            return
        self.write_plans()
        
    def queue_exposure(self):
        if not self.seqnum in self.planned_tiles:
            print('No more tiles in the plan (seqnum = %i)!' % self.seqnum)
            return None

        j = self.planned_tiles[self.seqnum]
        self.seqnum += 1

        if self.rc is None:
            print('Not actually queuing exposure (--no-queue):', j)
        else:
            print('Queuing exposure:', j)
            self.rc.addexposure(
                filter=j['filter'],
                ra=j['RA'],
                dec=j['dec'],
                object=j['object'],
                exptime=j['expTime'],
                )
        
        return j
    
    def filter_new_files(self, fns):
        return [fn for fn in fns if
                fn.endswith('.fits.fz') or fn.endswith('.fits')]

    def process_file(self, fn):
        ext = self.opt.ext
        print('%s: found new image %s' % (str(ephem.now()), fn))

        # Read primary FITS header
        phdr = fitsio.read_header(fn)
        expnum = phdr.get('EXPNUM', 0)
    
        obstype = phdr.get('OBSTYPE','').strip()
        print('Obstype:', obstype)
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
        print('Object:', obj)
        
        # Measure the new image
        kwa = {}
        if ext is not None:
            kwa.update(ext=ext)
        ps = None
        M = measure_raw(fn, ps=ps, **kwa)
    
        # Sanity checks
        ok = (M['nmatched'] >= 20) and (M.get('zp',None) is not None)
        if not ok:
            print('Failed sanity checks in our measurement of', fn, '-- not updating anything')
            # FIXME -- we could fall back to pass 3 here.
            return False
        
        # Choose the pass
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

        nextpass = choose_pass(trans, seeing, skybright, nomsky)
    
        # Choose the next tile from the right JSON tile list
        J = [self.J1,self.J2,self.J3][nextpass-1]

        J = self.tiles_after_now(J)
        if len(J) == 0:
            print('Could not find a JSON observation in pass', nextpass,
                  'with approx_datetime after now =', str(ephem.now()))
            return False
        
        # Update the exposure times in plan J based on measured conditions.
        print('Updating exposure times for pass', nextpass)
        for jplan in J:
            tilename = str(jplan['object'])
            # Find this tile in the tiles table.
            tile = get_tile_from_name(tilename, self.tiles)
            ebv = tile.ebv_med
            nextband = str(jplan['filter'])[0]
            #print('Selected tile:', tilename, nextband)
            rastr  = ra2hms (jplan['RA' ])
            decstr = dec2dms(jplan['dec'])
            ephemstr = str('%s,f,%s,%s,20' % (tilename, rastr, decstr))
            etile = ephem.readdb(ephemstr)
            etile.compute(self.obs)
            airmass = get_airmass(float(etile.alt))
            #print('Airmass of planned tile:', airmass)
    
            if M['band'] == nextband:
                nextsky = skybright
            else:
                # Guess that the sky is as much brighter than canonical
                # in the next band as it is in this one!
                nextsky = ((skybright - nomsky) + self.nom.sky(nextband))

            fid = self.nom.fiducial_exptime(nextband)
            expfactor = exposure_factor(fid, self.nom,
                                        airmass, ebv, seeing, nextsky, trans)
            #print('Tile', tilename)
            #print('Exposure factor:', expfactor)
            exptime = expfactor * fid.exptime
    
            ### HACK -- safety factor!
            #print('Exposure time:', exptime)
            exptime *= 1.1
            exptime = int(np.ceil(exptime))
            #print('Exposure time with safety factor:', exptime)

            exptime = np.clip(exptime, fid.exptime_min, fid.exptime_max)
            #print('Clipped exptime', exptime)
            if nextband == 'z':
                # Compute cap on exposure time to avoid saturation /
                # loss of dynamic range.
                t_sat = self.nom.saturation_time(nextband, nextsky)
                if exptime > t_sat:
                    exptime = t_sat
                    print('Reduced exposure time to avoid z-band saturation: %.1f' % exptime)
            exptime = int(exptime)
    
            #print('Changing exptime from', jplan['expTime'], 'to', exptime)
            jplan['expTime'] = exptime
            
        self.plan_tiles(J)
        return True

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
    
    def plan_tiles(self, J, Nahead=10, exptime=None):
        '''
        Nahead: int: How many exposures ahead should we plan?
        '''

        # Set observing conditions for computing exposure time
        now = ephem.now()
        self.obs.date = now
        
        self.upcoming = []

        iahead = 0
        for ii,jplan in enumerate(J):
            if iahead >= Nahead:
                break
            tilename = str(jplan['object'])
            nextseq = self.seqnum + iahead

            print('Considering planning tile %s for exp %i' %
                  (tilename, nextseq))

            if tilename in self.observed_tiles:
                oldfn = self.observed_tiles[tilename]
                print('Tile %s was observed in file %s' % (tilename, oldfn))
                continue
            
            # Check all planned tiles before this one for a duplicate tile.
            dup = False
            for s in range(nextseq-1, 0, -1):
                t = self.planned_tiles[s]
                if t['object'] == tilename:
                    dup = True
                    print('Wanted to plan tile %s for exp %i '
                          % (tilename, nextseq),
                          'but it was already planned for exp %i' % s)
                    break
            if dup:
                continue
    
            iahead += 1

            if exptime is not None:
                jplan['expTime'] = exptime

            print('%s: updating exposure %i to tile %s' %
                  (str(ephem.now()), nextseq, tilename))
            self.planned_tiles[nextseq] = jplan
            self.upcoming.append(jplan)
            self.obs.date += (jplan['expTime'] + self.nom.overhead) / 86400.

        self.write_plans()            

    def write_plans(self):
        # Write upcoming plan to a JSON file
        fn = 'decbot-plan.json'
        tmpfn = fn + '.tmp'
        jstr = json.dumps(self.upcoming, sort_keys=True,
                          indent=4, separators=(',', ': '))
        f = open(tmpfn, 'w')
        f.write(jstr + '\n')
        f.close()
        os.rename(tmpfn, fn)
        print('Wrote', fn)

        fn = 'decbot-plan-5.json'
        tmpfn = fn + '.tmp'
        jstr = json.dumps(self.upcoming[:5], sort_keys=True,
                          indent=4, separators=(',', ': '))
        f = open(tmpfn, 'w')
        f.write(jstr + '\n')
        f.close()
        os.rename(tmpfn, fn)
        print('Wrote', fn)

        # Write a FITS table of the exposures we think we've queued,
        # the ones we have planned, and the future tiles in passes 1,2,3.
        P = ([(self.planned_tiles[s],'Q') for s in range(self.seqnum)] +
             [(j,'P') for j in self.upcoming])
             
        # Skip ones scheduled for before now
        now = ephem.now()
        for i,J in enumerate([self.J1,self.J2,self.J3]):
            passnum = i+1
            for j in J:
                tstart = ephem.Date(str(j['approx_datetime']))
                if tstart < now:
                    continue
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
        print('Wrote', fn)
    

if __name__ == '__main__':
    main()
    
