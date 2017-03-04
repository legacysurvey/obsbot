'''
TO-DO:

-- FIXMEs -- filter changing -- here we set filter ONCE PER NIGHT!
'''

from __future__ import print_function
import sys
import time
import json
import datetime
import stat
import os
from collections import OrderedDict
from glob import glob

import matplotlib
matplotlib.use('Agg')

import numpy as np
import ephem
import fitsio
from astrometry.util.fits import fits_table

from measure_raw import measure_raw
from jnox import (jnox_preamble, jnox_filter, jnox_moveto, jnox_exposure,
                  jnox_readout, jnox_end_exposure, jnox_cmds_for_json,
                  ra2hms, dec2dms)
from obsbot import (exposure_factor, get_tile_from_name, get_airmass,
                    Obsbot, choose_pass, mjdnow)

def main(cmdlineargs=None, get_mosbot=False):
    import optparse
    parser = optparse.OptionParser(usage='%prog <pass1.json> <pass2.json> <pass3.json>')

    from camera_mosaic import (ephem_observer, default_extension, nominal_cal,
                               tile_path)

    parser.add_option('--no-filter', dest='set_filter', action='store_false', default=True,
                      help='Do not set the filter at the beginning of tonight.sh')

    parser.add_option('--pass', dest='passnum', type=int, default=2,
                      help='Set default pass number (1/2/3), default 2')
    parser.add_option('--exptime', type=int, default=None,
                      help='Set default exposure time, default whatever is in the JSON files, usually 80 sec')

    parser.add_option('--no-cut-past', dest='cut_before_now',
                      default=True, action='store_false',
                      help='Do not cut tiles that were supposed to be observed in the past')

    parser.add_option('--rawdata', help='Directory to monitor for new images: default $MOS3_DATA if set, else "rawdata"', default=None)
    parser.add_option('--script', dest='scriptfn', help='Write top-level shell script, default is %default', default='/mosaic3/exec/mosbot/tonight.sh')
    parser.add_option('--no-write-script', dest='write_script', default=True, action='store_false')
    parser.add_option('--ext', help='Extension to read for computing observing conditions, default %default', default=default_extension)
    parser.add_option('--tiles',
                      default=tile_path,
                      help='Observation status file, default %default')

    parser.add_option('--adjust', action='store_true',
                      help='Adjust exposure times based on previous passes?')

    parser.add_option('--no-db', dest='db', default=True, action='store_false',
                      help='Do not update the computed-exposure-time database')
    
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

    #ADM first delete any existing .sh scripts to guard against 
    #ADM the user executing an old script. Also delete forcepass
    #ADM and nopass files
    scriptdir = os.path.dirname(opt.scriptfn)
    if len(scriptdir) and os.path.exists(scriptdir):
        for file in (glob(os.path.join(scriptdir, '*.sh')) +
                     glob(os.path.join(scriptdir, 'forcepass?')) +
                     glob(os.path.join(scriptdir, 'nopass1'))):
            os.remove(file)

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

    obs = ephem_observer()
    
    print('Reading tiles table', opt.tiles)
    tiles = fits_table(opt.tiles)
    ## Cut tiles to only passes 1,2,3
    tiles.cut(tiles.get('pass') <= 3)
    
    if opt.rawdata is None:
        opt.rawdata = os.environ.get('MOS3_DATA', 'rawdata')
    
    mosbot = Mosbot(J1, J2, J3, opt, nominal_cal, obs, tiles)
    if get_mosbot:
        return mosbot
    mosbot.run()


class Mosbot(Obsbot):
    def __init__(self, J1, J2, J3, opt, nom, obs, tiles):
        super(Mosbot, self).__init__(opt.rawdata, backlog=False,
                                     only_process_newest=True)
        # How many exposures ahead should we write?
        self.Nahead = 10
    
        self.timeout = None
        self.J1 = J1
        self.J2 = J2
        self.J3 = J3
        self.opt = opt
        self.nom = nom
        self.obs = obs
        self.tiles = tiles

        self.scriptdir = os.path.dirname(opt.scriptfn)
        # Create scriptdir (one path component only), if necessary
        if len(self.scriptdir) and not os.path.exists(self.scriptdir):
            os.mkdir(self.scriptdir)
    
        self.seqnumfn = 'seqnum.txt'
        self.seqnumpath = os.path.join(self.scriptdir, self.seqnumfn)
    
        self.expscriptpattern  = 'expose-%i.sh'
        self.slewscriptpattern = 'slewread-%i.sh'

        self.planned_tiles = OrderedDict()

        if opt.write_script:
            # Default to Pass 2!
            J = [J1,J2,J3][opt.passnum - 1]
            self.write_initial_script(J, opt.passnum, opt.exptime,
                                      opt.scriptfn, self.seqnumfn)
            self.n_exposures = len(J)
        else:
            self.n_exposures = 0

        self.bot_runtime = mjdnow()
            
    def filter_new_files(self, fns):
        return [fn for fn in fns if
                fn.endswith('.fits.fz') or fn.endswith('.fits')]

    def write_initial_script(self, J, passnum, exptime, scriptfn, seqnumfn):
        quitfile = 'quit'

        # 775
        chmod = (stat.S_IXUSR | stat.S_IWUSR | stat.S_IRUSR |
                 stat.S_IXGRP | stat.S_IWGRP | stat.S_IRGRP |
                 stat.S_IXOTH |                stat.S_IROTH)
        
        def log(s):
            datecmd = 'date -u +"%Y-%m-%d %H:%M:%S"'
            #return 'echo "mosbot $(%s): %s"' % (datecmd, s)
            return 'echo "$(%s): %s" >> mosbot.log' % (datecmd, s)

        # Delete 'quit' file if it exists.
        quitfn = os.path.join(self.scriptdir, quitfile)
        if os.path.exists(quitfn):
            print('Removing file', quitfn)
            os.remove(quitfn)
        
        # Write "read.sh" for quitting gracefully without slew
        path = os.path.join(self.scriptdir, 'read.sh')
        f = open(path, 'w')
        f.write(jnox_readout() + '\n' + jnox_end_exposure())
        f.close()
        os.chmod(path, chmod)

        # Start writing the "tonight.sh" script.
        script = []
        # tonight.sh initial NOCS commands
        script.append(jnox_preamble(scriptfn))

        script.append(log('###############'))
        script.append(log('Starting script'))
        script.append(log('###############'))

        if self.opt.set_filter:
            # tonight.sh setting the filter once at the start of the night!
            j = J[0]
            f = str(j['filter'])
            script.append(jnox_filter(f))
        
        # Write top-level script and shell scripts for default plan.
        for i,j in enumerate(J):

            # The sequence number start at 1
            seq = i + 1

            # tonight.sh: start exposure
            script.append('\n### Exposure %i ###\n' % seq)
            script.append('echo "%i" > %s' % (seq, seqnumfn))
            
            if seq > 1:
                script.append('# Check for file "%s"; read out and quit if it exists' % quitfile)
                script.append('if [ -f %s ]; then\n  . read.sh; rm %s; exit 0;\nfi' %
                              (quitfile,quitfile))

                # tonight.sh: slew to field while reading out previous exposure
                # Write slewread-##.sh script
                fn = self.slewscriptpattern % seq
                path = os.path.join(self.scriptdir, fn)
                f = open(path, 'w')
                f.write(slewscript_for_json(j))
                f.close()
                os.chmod(path, chmod)
                print('Wrote', path)
                script.append('. %s' % fn)
                # Another check for 'quit' file; no need to run "read.sh" here.
                script.append('if [ -f %s ]; then\n  rm %s; exit 0;\nfi' %
                              (quitfile,quitfile))
                

            script.append(log('Starting exposure %i' % seq))
            expfn = self.expscriptpattern % seq
            script.append(log('Tile: $(grep "Tile:" %s)' % expfn))
            
            tilename = str(j['object'])
            self.planned_tiles[seq] = tilename
            
            if exptime is not None:
                j['expTime'] = exptime
            
            # Write expose-##.sh default exposure script
            fn = self.expscriptpattern % seq
            path = os.path.join(self.scriptdir, fn)
            f = open(path, 'w')
            s = ('# Exp %i: Tile: %s, Pass: %i; default\n' %
                 (seq, tilename, passnum))

            ra  = ra2hms (j['RA' ])
            dec = dec2dms(j['dec'])
            status = ('Exp %i: Tile %s, Pass %i, RA %s, Dec %s' %
                      (seq, tilename, passnum, ra, dec))
            
            s += expscript_for_json(j, status=status)
            f.write(s)
            f.close()
            os.chmod(path, chmod)
            print('Wrote', path)
            script.append('. %s' % fn)
    
        # Read out the last one!
        script.append('. read.sh')
    
        script = '\n'.join(script) + '\n'
        
        f = open(scriptfn, 'w')
        f.write(script)
        f.close()
        os.chmod(scriptfn, chmod)
        print('Wrote', scriptfn)

    def process_file(self, fn):
        print('%s: found new image %s' % (str(ephem.now()), fn))

        if not self.check_header(fn):
            return False

        # Measure the new image
        kwa = {}
        ext = self.opt.ext
        if ext is not None:
            kwa.update(ext=ext)
        ps = None
        M = measure_raw(fn, ps=ps, **kwa)

        # Basic checks
        ok = (M['nmatched'] >= 20) and (M.get('zp',None) is not None)
        if not ok:
            print('Failed basic checks in our measurement of', fn,
                  '-- not updating anything')
            # FIXME -- we could fall back to pass 3 here.
            return False

        return self.update_for_image(M)

    def check_header(self, fn):
        # Read primary FITS header
        phdr = fitsio.read_header(fn)
    
        obstype = phdr.get('OBSTYPE','').strip()
        print('obstype:', obstype)
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
        if filt == 'solid':
            print('Solid (block) filter.')
            return False

        return True

    def update_for_image(self, M):
        # filename patterns for the exposure and slew scripts
        expscriptpat  = os.path.join(self.scriptdir, self.expscriptpattern)
        slewscriptpat = os.path.join(self.scriptdir, self.slewscriptpattern)
        
        # Choose the pass, based on...
        trans  = M['transparency']
        seeing = M['seeing']
        skybright = M['skybright']
        # eg, nominal = 20, sky = 19, brighter is 1 mag brighter than nom.
        meas_band = M['band']
        nomsky = self.nom.sky(meas_band)
        brighter = nomsky - skybright
    
        print('Transparency: %6.02f' % trans)
        print('Seeing      : %6.02f' % seeing)
        print('Sky         : %6.02f' % skybright)
        print('Nominal sky : %6.02f' % nomsky)
        print('Sky over nom: %6.02f   (positive means brighter than nominal)' %
              brighter)

        nextpass = choose_pass(trans, seeing, skybright, nomsky,
                               forcedir=self.scriptdir)

        print()
        print('Selected pass:', nextpass)
        print()

        # Choose the next tile from the right JSON tile list Jp
        J = [self.J1,self.J2,self.J3][nextpass-1]

        now = ephem.now()

        # Read the current sequence number
        print('%s: reading sequence number from %s' %
              (str(ephem.now()), self.seqnumpath))
        f = open(self.seqnumpath, 'r')
        s = f.read()
        f.close()
        seqnum = int(s)
        print('%s: sequence number: %i' % (str(now), seqnum))
        
        # 'iplan': the tile index we will use for exposure # 'seqnum'
        iplan = None
        if self.opt.cut_before_now:
            # The usual case
            for i,j in enumerate(J):
                tstart = ephem.Date(str(j['approx_datetime']))
                if tstart <= now:
                    continue
                print('Found tile', j['object'], 'which starts at', str(tstart))
                iplan = i
                break
            if iplan is None:
                print('Could not find a JSON observation in pass', nextpass,
                      'with approx_datetime after now =', str(now),
                      '-- latest one', str(tstart))
                return False
        else:
            # For when we want to observe every tile in the plan:
            # 'seqnum' is the exposure currently running;
            # seqnum is 1-indexed, so that's the index we want for iplan.
            iplan = seqnum

        # Set observing conditions for computing exposure time
        self.obs.date = now

        # Planned exposures:
        P = fits_table()
        P.type = []
        P.tilename = []
        P.filter = []
        P.exptime = []
        P.ra = []
        P.dec = []
        P.passnumber = []
        
        iahead = 0
        for ii,jplan in enumerate(J[iplan:]):
            if iahead >= self.Nahead:
                break
            tilename = str(jplan['object'])
            nextseq = seqnum + 1 + iahead

            if self.n_exposures > 0 and nextseq > self.n_exposures:
                print('Planned tile is beyond the exposure number in ',
                      'tonight.sh -- RESTART MOSBOT')
                return False
            
            print('Considering planning tile %s for exp %i'%(tilename,nextseq))

            # Check all planned tiles before this one for a duplicate tile.
            dup = False
            for s in range(nextseq-1, 0, -1):
                if tilename == self.planned_tiles[s]:
                    dup = True
                    print('Wanted to plan tile %s (pass %i element %i) for exp %i'
                          % (tilename, nextpass, iplan+ii, nextseq),
                          'but it was already planned for exp %i' % s)
                    break
            if dup:
                continue

            self.planned_tiles[nextseq] = tilename
            iahead += 1

            largeslew = False
            if iahead == 1:
                from astrometry.util.starutil_numpy import degrees_between
                # Check for large slew
                prevname = self.planned_tiles[nextseq - 1]
                Jall = self.J1 + self.J2 + self.J3
                I, = np.nonzero([str(j['object']) == prevname for j in Jall])
                if len(I) == 0:
                    print('Could not find previous tile "%s"' % prevname)
                else:
                    jprev = Jall[I[0]]
                    slew = degrees_between(jprev['RA'], jprev['dec'],
                                           jplan['RA'], jplan['dec'])
                    print('Slew: %.1f degrees' % slew)
                    largeslew = (slew > 10)
                if largeslew:
                    print()
                    print()
                    print('Large slew from current exposure to next: ' +
                          'RA,Dec %.1f, %.1f to %.1f, %.1f ==> %.1f degrees' %
                          (jprev['RA'], jprev['dec'], jplan['RA'],
                           jplan['dec'], slew))
                    print()
                    print()
                    os.system('tput bel; sleep 0.2; ' * 5)
                
            # Find this tile in the tiles table.
            tile = get_tile_from_name(tilename, self.tiles)
            ebv = tile.ebv_med
            nextband = str(jplan['filter'])[0]

            print('Selected tile:', tile.tileid, nextband)
            rastr  = ra2hms (jplan['RA' ])
            decstr = dec2dms(jplan['dec'])
            ephemstr = str('%s,f,%s,%s,20' % (tilename, rastr, decstr))
            etile = ephem.readdb(ephemstr)
            etile.compute(self.obs)
            airmass = get_airmass(float(etile.alt))
            print('Airmass of planned tile:', airmass)

            if M['band'] == nextband:
                nextsky = skybright
            else:
                # Guess that the sky is as much brighter than canonical
                # in the next band as it is in this one!
                nextsky = ((skybright - nomsky) + self.nom.sky(nextband))

            fid = self.nom.fiducial_exptime(nextband)
            expfactor = exposure_factor(fid, self.nom,
                                        airmass, ebv, seeing, nextsky, trans)
            print('Exposure factor:', expfactor)

            expfactor_orig = expfactor

            # Adjust for previous exposures?
            if self.opt.adjust:
                debug=True
                adjfactor,others = self.adjust_for_previous(
                    tile, nextband, fid, debug=debug, get_others=True)
                # Don't adjust exposure times down, only up.
                adjfactor = max(adjfactor, 1.0)
                expfactor *= adjfactor
            else:
                adjfactor = 0.
                others = []

            exptime = expfactor * fid.exptime

            ### HACK -- safety factor!
            print('Exposure time:', exptime)
            exptime *= 1.1
            print('Exposure time with safety factor:', exptime)
            exptime_unclipped = exptime

            exptime = np.clip(exptime, fid.exptime_min, fid.exptime_max)
            print('Clipped exptime', exptime)

            exptime_satclipped = 0.
            if nextband == 'z':
                # Compute cap on exposure time to avoid saturation /
                # loss of dynamic range.
                t_sat = self.nom.saturation_time(nextband, nextsky)
                if exptime > t_sat:
                    exptime_satclipped = t_sat
                    exptime = t_sat
                    print('Reduced exposure time to avoid z-band saturation:',
                    '%.1f' % exptime)
            exptime = int(np.ceil(exptime))

            print('Changing exptime from', jplan['expTime'], 'to', exptime)
            jplan['expTime'] = exptime

            print('Predict tile will be observed at', str(self.obs.date),
                  'vs approx_datetime', jplan['approx_datetime'])

            # Update the computed exposure-time database.
            if self.opt.db:
                try:
                    # NOTE, these kwargs MUST match the names in models.py
                    self.update_exptime_db(
                        nextseq, others, tileid=tile.tileid,
                        passnumber=nextpass, band=nextband, airmass=airmass,
                        ebv=ebv, meas_band=meas_band,
                        zeropoint=M['zp'], transparency=trans,
                        seeing=seeing, sky=skybright, expfactor=expfactor_orig,
                        adjfactor=adjfactor,
                        exptime_unclipped=exptime_unclipped,
                        exptime_satclipped=exptime_satclipped,
                        exptime=exptime)
                except:
                    print('Failed to update computed-exptime database.')
                    import traceback
                    traceback.print_exc()
                    # carry on

            status = ('Exp %i: Tile %s, Pass %i, RA %s, Dec %s' %
                      (nextseq, tilename, nextpass, rastr, decstr))

            print('%s: updating exposure %i to tile %s' %
                  (str(ephem.now()), nextseq, tilename))

            expscriptfn = expscriptpat % (nextseq)
            exptmpfn = expscriptfn + '.tmp'
            f = open(exptmpfn, 'w')
            f.write(('# Exp %i Tile: %s, set at Seq %i, %s\n' %
                     (nextseq, tilename, seqnum, str(ephem.now()))) +
                    expscript_for_json(jplan, status=status))
            f.close()

            slewscriptfn = slewscriptpat % (nextseq)
            slewtmpfn = slewscriptfn + '.tmp'
            f = open(slewtmpfn, 'w')
            f.write(slewscript_for_json(jplan))
            f.close()

            os.rename(exptmpfn, expscriptfn)
            print('Wrote', expscriptfn)
            os.rename(slewtmpfn, slewscriptfn)
            print('Wrote', slewscriptfn)

            self.obs.date += (exptime + self.nom.overhead) / 86400.

            print('%s: updated exposure %i to tile %s' %
                  (str(ephem.now()), nextseq, tilename))

            P.tilename.append(tilename)
            P.filter.append(nextband)
            P.exptime.append(exptime)
            P.ra.append(jplan['RA'])
            P.dec.append(jplan['dec'])
            P.passnumber.append(nextpass)
            P.type.append('P')

        if iahead < self.Nahead:
            # We ran out of tiles.
            # Overwrite the default planned exposures with blanks to avoid
            # running the defaults (eg, at end of night).
            seqstart = seqnum + 1 + iahead
            seqend = self.n_exposures - 1
            print('Overwriting remaining exposure scripts (%i to %i inclusive) with blanks.' % (seqstart, seqend))
            for nextseq in range(seqstart, seqend+1):
                print('Blanking out exposure %i' % nextseq)
                fn = self.expscriptpattern % nextseq
                path = os.path.join(self.scriptdir, fn)
                f = open(path, 'w')
                f.write('\n')
                f.close()
                fn = self.slewscriptpattern % nextseq
                path = os.path.join(self.scriptdir, fn)
                f = open(path, 'w')
                f.write('\n')
                f.close()

        for i,J in enumerate([self.J1,self.J2,self.J3]):
            passnum = i+1
            for j in J:
                tstart = ephem.Date(str(j['approx_datetime']))
                if tstart < now:
                    continue
                P.tilename.append(str(j['object']))
                filt = str(j['filter'])[0]
                P.filter.append(filt)
                P.exptime.append(j['expTime'])
                P.ra.append(j['RA'])
                P.dec.append(j['dec'])
                P.passnumber.append(passnum)
                P.type.append('%i' % passnum)

        P.to_np_arrays()
        fn = 'mosbot-plan.fits'
        tmpfn = fn + '.tmp'
        P.writeto(tmpfn)
        os.rename(tmpfn, fn)
        print('Wrote', fn)

        return True

    def update_exptime_db(self, seqnum, others, **kwargs):
        from obsdb.models import ComputedExptime, OtherPasses
        c,created = ComputedExptime.objects.get_or_create(
            starttime=self.bot_runtime, seqnum=seqnum)
        for k,v in kwargs.items():
            setattr(c, k, v)
        c.save()

        for o in others:
            op,created = OtherPasses.objects.get_or_create(
                exposure=c, tileid=o.tileid)
            op.passnumber = o.passnum
            op.depth      = o.depth
            op.save()

    # For adjust_for_previous: find nearby tiles
    def other_passes(self, tile, tiles):
        '''
        Given tile number *tile*, return the obstatus rows for the other passes
        on this tile center.

        Returns: *otherpasses*, table object
        '''
        # For MzLS, the pass1, pass2, pass3 tiles are offset from each
        # other by quite a lot.  Therefore, we want to look at not
        # only the nearest tile, but all that overlap significantly
        # (ie, not just a little overlap at the edges).
        dd = (4096 - 400) * 2 * 0.262 / 3600.
        #print('dd', dd)
        cosdec = np.cos(np.deg2rad(tile.dec))
        dra = dd / cosdec
        #print('dra', dra)

        deltara = np.abs(tile.ra  - tiles.ra )
        I = np.flatnonzero((np.minimum(deltara, 360.-deltara) < dra) *
                           (np.abs(tile.dec - tiles.dec) <  dd))
        tt = tiles[I]
        # Omit 'tile' itself...
        tt.cut(tt.tileid != tile.tileid)

        print('Tile', tile.tileid, 'pass', tile.get('pass'), 'RA,Dec', tile.ra, tile.dec)
        print('Found overlapping tiles:')
        for t in tt:
            print(' tile', t.tileid, 'pass', t.get('pass'), 'RA,Dec', t.ra, t.dec, '(delta', t.ra - tile.ra, t.dec - tile.dec, ')')

        return tt
        

def expscript_for_json(j, status=None):

    ss = ('# Exposure scheduled for: %s UT\n' % j['approx_datetime'])

    # Expected start time of this script, in "date -u +%s" units -- seconds since unixepoch = 1970-01-01 00:00:00
    unixepoch = 25567.5
    tj = ephem.Date(str(j['approx_datetime']))

    ss = ['dt=$(($(date -u +%%s) - %i))' % (int((tj - unixepoch) * 86400))]
    #ss.append('echo DT: $dt\n')
    ss.append('if [ $dt -gt 3600 ]; then\n' +
              '  echo; echo "WARNING: This exposure is happening more than an hour late!";\n' +
              '  echo "Scheduled time: %s UT";\n' % j['approx_datetime'] +
              '  echo; tput bel; sleep 0.25; tput bel; sleep 0.25; tput bel;\n' +
              'fi\n' +
              'if [ $dt -lt -3600 ]; then\n' +
              '  echo; echo "WARNING: This exposure is happening more than an hour early!";\n' +
              '  echo "Scheduled time: %s UT";\n' % j['approx_datetime'] +
              '  echo; tput bel; sleep 0.25; tput bel; sleep 0.25; tput bel;\n' +
              'fi\n')
    #ADM force exit if observation is more than 4 hours late
    ss.append('if [ $dt -gt 14400 ]; then\n' +
              '  echo; echo "ERROR: This exposure is happening more than 4 hours late!";\n' +
              '  echo "Scheduled time: %s UT";\n' % j['approx_datetime'] +
              '  echo "YOU SHOULD RERUN mosbot.py TO CREATE A NEW tonight.sh FILE!";\n' +
              '  touch quit\n' +
              '  exit 1\n' +
              'fi\n')
    
    ra  = j['RA']
    dec = j['dec']
    ss.append(jnox_moveto(ra, dec))

    # exposure time
    et = (j['expTime'])
    # filter
    filter_name = str(j['filter'])
    tilename = str(j['object'])
    if status is None:
        ra  = ra2hms(ra)
        dec = dec2dms(dec)
        status = "Tile %s, RA %s, Dec %s" % (tilename, ra, dec)
    ss.append(jnox_exposure(et, filter_name, tilename, status))
    return '\n'.join(ss) + '\n'

def slewscript_for_json(j):
    ra  = j['RA']
    dec = j['dec']
    ss = [jnox_moveto(ra, dec)]
    # second call concludes readout
    ss.append(jnox_readout())
    # endobs
    ss.append(jnox_end_exposure())
    return '\n'.join(ss) + '\n'



if __name__ == '__main__':
    import obsdb
    from camera_mosaic import database_filename
    obsdb.django_setup(database_filename=database_filename)

    main()
    
