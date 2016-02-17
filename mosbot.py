'''
TO-DO:

-- t_sat_max should be recomputed based on skybright!!
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

import matplotlib
matplotlib.use('Agg')

import ephem

from astrometry.util.plotutils import *
from astrometry.util.fits import *

from mosaicstrategy import (
    ExposureFactor, getParserAndGlobals, setupGlobals,
    GetAirmass, StartAndEndTimes, s_to_days, readTilesTable, GetNightlyStrategy,
    WriteJSON)
from measure_raw import measure_raw, get_default_extension

from jnox import *
from copilot import get_tile_from_name

def main(cmdlineargs=None, get_mosbot=False):
    import optparse
    parser = optparse.OptionParser(usage='%prog <pass1.json> <pass2.json> <pass3.json>')
    
    parser.add_option('--rawdata', help='Directory to monitor for new images: default $MOS3_DATA if set, else "rawdata"', default=None)
    parser.add_option('--script', dest='scriptfn', help='Write top-level shell script, default is %default', default='/mosaic3/exec/mosbot/tonight.sh')
    parser.add_option('--no-write-script', dest='write_script', default=True, action='store_false')
    parser.add_option('--ext', help='Extension to read for computing observing conditions, default %default', default='im4')
    parser.add_option('--tiles',
                      default='obstatus/mosaic-tiles_obstatus.fits',
                      help='Observation status file, default %default')

    parser.add_option('--pass', dest='passnum', type=int, default=2,
                      help='Set default pass number (1/2/3), default 2')
    parser.add_option('--exptime', type=int, default=None,
                      help='Set default exposure time, default whatever is in the JSON files, usually 80 sec')

    parser.add_option('--no-cut-past', dest='cut_before_now',
                      default=True, action='store_false',
                      help='Do not cut tiles that were supposed to be observed in the past')
    
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

    # nightlystrategy setup
    parser,gvs = getParserAndGlobals()
    nsopt,nsargs = parser.parse_args('--date 2015-01-01 --pass 1 --portion 1'.split())
    obs = setupGlobals(nsopt, gvs)

    print('Reading tiles table', opt.tiles)
    tiles = fits_table(opt.tiles)
    
    if opt.rawdata is None:
        opt.rawdata = os.environ.get('MOS3_DATA', 'rawdata')
    
    mosbot = Mosbot(J1, J2, J3, opt, gvs, obs, tiles)
    if get_mosbot:
        return mosbot
    mosbot.run()


class Mosbot(object):
    def __init__(self, J1, J2, J3, opt, gvs, obs, tiles):
        self.J1 = J1
        self.J2 = J2
        self.J3 = J3
        self.gvs = gvs
        self.obs = obs
        self.opt = opt
        self.tiles = tiles
        
        self.imagedir = opt.rawdata

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
            self.write_initial_script(J, opt.passnum, opt.exptime, opt.scriptfn, self.seqnumfn)

    def write_initial_script(self, J, passnum, exptime, scriptfn, seqnumfn):
        quitfile = 'quit'
        nowscriptpattern  = 'now-%i.sh'

        # 775
        chmod = (stat.S_IXUSR | stat.S_IWUSR | stat.S_IRUSR |
                 stat.S_IXGRP | stat.S_IWGRP | stat.S_IRGRP |
                 stat.S_IXOTH |                stat.S_IROTH)
        
        def log(s):
            datecmd = 'date -u +"%Y-%m-%d %H:%M:%S"'
            #return 'echo "mosbot $(%s): %s"' % (datecmd, s)
            return 'echo "$(%s): %s" >> mosbot.log' % (datecmd, s)
    
        # Write "read.sh" for quitting gracefully without slew
        path = os.path.join(self.scriptdir, 'read.sh')
        f = open(path, 'w')
        f.write(jnox_readout() + '\n' + jnox_end_exposure())
        f.close()
        os.chmod(path, chmod)

        script = []
        script.append(jnox_preamble(scriptfn))

        script.append(log('###############'))
        script.append(log('Starting script'))
        script.append(log('###############'))

        j = J[0]
        f = str(j['filter'])
        script.append(jnox_filter(f))
        
        # Write top-level script and shell scripts for default plan.
        for i,j in enumerate(J):

            # We make the sequence numbers start at 1, just to make things
            # difficult for ourselves.
            seq = i + 1
            
            if seq > 1:

                script.append('# Check for file "%s"; read out and quit if it exists' % quitfile)
                script.append('if [ -f %s ]; then\n  . read.sh; rm %s; exit 0;\nfi' %
                              (quitfile,quitfile))
                
                # Write slewread-##.sh script
                fn = self.slewscriptpattern % seq
                path = os.path.join(self.scriptdir, fn)
                f = open(path, 'w')
                f.write(slewscript_for_json(j))
                f.close()
                os.chmod(path, chmod)
                print('Wrote', path)
                script.append('. %s' % fn)

            script.append('\n### Exposure %i ###\n' % seq)
            script.append('echo "%i" > %s' % (seq, seqnumfn))
            script.append(log('Starting exposure %i' % seq))
            expfn = self.expscriptpattern % seq
            script.append(log('Tile: $(grep "Tile:" %s)' % expfn))
            
            tilename = str(j['object'])
            self.planned_tiles[seq] = tilename
            
            # Write now-##.sh
            fn = nowscriptpattern % seq
            path = os.path.join(self.scriptdir, fn)
            f = open(path, 'w')
            f.write('echo "Do stuff before exposure %i?"' % seq)
            f.close()
            os.chmod(path, chmod)
            print('Wrote', path)
            script.append('. %s' % fn)

            ra  = ra2hms (j['RA' ])
            dec = dec2dms(j['dec'])
            status = ('Exp %i: Tile %s, Pass %i, RA %s, Dec %s' %
                      (seq, tilename, passnum, ra, dec))

            if exptime is not None:
                j['expTime'] = exptime
            
            # Write expose-##.sh
            fn = self.expscriptpattern % seq
            path = os.path.join(self.scriptdir, fn)
            f = open(path, 'w')
            s = ('# Exp %i: Tile: %s, Pass: %i; default\n' %
                 (seq, tilename, passnum))
            
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


    ### Borrowed from copilot.py ...
    def get_new_images(self):
        images = set(os.listdir(self.imagedir))
        newimgs = images - self.oldimages
        newimgs = list(newimgs)
        newimgs = [fn for fn in newimgs if
                   fn.endswith('.fits.fz') or fn.endswith('.fits')]
        return newimgs

    def get_newest_image(self, newimgs):
        if len(newimgs) == 0:
            return None
        # Take the one with the latest timestamp.
        latest = None
        newestimg = None
        for fn in newimgs:
            st = os.stat(os.path.join(self.imagedir, fn))
            t = st.st_mtime
            if latest is None or t > latest:
                newestimg = fn
                latest = t
        return newestimg

        
    def run(self):
        print('Watching directory "%s" for new files' % self.imagedir)
        self.oldimages = set(os.listdir(self.imagedir))

        while True:
            time.sleep(5.)
            newimgs = self.get_new_images()
            if len(newimgs) == 0:
                continue
            newestfn = self.get_newest_image(newimgs)

            try:
                ok = self.found_new_image(newestfn)
                if ok:
                    self.oldimages.update(newimgs)
                    
            except IOError:
                print('Failed to read FITS image:', fn)
                import traceback
                traceback.print_exc()
                continue
            
    def found_new_image(self, fn):
        ext = self.opt.ext

        expscriptpat  = os.path.join(self.scriptdir, self.expscriptpattern)
        slewscriptpat = os.path.join(self.scriptdir, self.slewscriptpattern)
        
        print('%s: found new image %s' % (str(ephem.now()), fn))

        nopass1path = os.path.join(self.scriptdir, 'nopass1')
        dopass1 = not os.path.exists(nopass1path)
        if not dopass1:
            print('Not considering Pass 1 because file exists:', nopass1path)
        nopass2path = os.path.join(self.scriptdir, 'nopass2')
        dopass2 = not os.path.exists(nopass2path)
        if not dopass2:
            print('Not considering Pass 2 because file exists:', nopass2path)
        
        # Read primary FITS header
        phdr = fitsio.read_header(fn)
        expnum = phdr.get('EXPNUM', 0)
    
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
        nomsky = self.gvs.sb_dict[M['band']]
        brighter = nomsky - skybright
    
        print('Transparency: %6.02f' % trans)
        print('Seeing      : %6.02f' % seeing)
        print('Sky         : %6.02f' % skybright)
        print('Nominal sky : %6.02f' % nomsky)
        print('Sky over nom: %6.02f   (positive means brighter than nom)' %
              brighter)
    
        transcut = 0.9
        seeingcut = 1.25
        brightcut = 0.25
        
        transok = trans > transcut
        seeingok = seeing < seeingcut
        brightok = brighter < brightcut
    
        pass1ok = transok and seeingok and brightok
        pass2ok = transok or  seeingok
        
        nextpass = 3
        if pass1ok and dopass1:
            nextpass = 1
    
        elif pass2ok and dopass2:
            nextpass = 2
    
        print('Transparency cut: %s       (%6.2f vs %6.2f)' %
              (('pass' if transok else 'fail'), trans, transcut))
        print('Seeing       cut: %s       (%6.2f vs %6.2f)' %
              (('pass' if seeingok else 'fail'), seeing, seeingcut))
        print('Brightness   cut: %s       (%6.2f vs %6.2f)' %
              (('pass' if brightok else 'fail'), skybright, nomsky+brightcut))
        print('Pass 1 = transparency AND seeing AND brightness: %s' % pass1ok)
        if pass1ok and not dopass1:
            print('Pass 1 forbidden by observer!')
        print('Pass 2 = transparency OR  seeing               : %s' % pass2ok)
        if pass2ok and not dopass2:
            print('Pass 2 forbidden by observer!')
        print('Selected pass:', nextpass)
    
        # Choose the next tile from the right JSON tile list Jp
        J = [self.J1,self.J2,self.J3][nextpass-1]
        now = ephem.now()
        iplan = None
        print('UTC now is', str(now))
        for i,j in enumerate(J):
            tstart = ephem.Date(str(j['approx_datetime']))
            if tstart > now:
                print('Found tile', j['object'], 'which starts at', str(tstart))
                iplan = i
                break
        if iplan is None:
            print('Could not find a JSON observation in pass', nextpass,
                  'with approx_datetime after now =', str(now),
                  '-- latest one', str(tstart))
            return False
    
        # Read the current sequence number
        print('%s: reading sequence number' % (str(ephem.now())))
        f = open(self.seqnumpath, 'r')
        s = f.read()
        f.close()
        seqnum = int(s)
        print('%s: sequence number: %i' % (str(ephem.now()), seqnum))
        
        # How many exposures ahead should we write?
        Nahead = 10
    
        # Set observing conditions for computing exposure time
        self.gvs.transparency = trans
        self.obs.date = now
    
        P = fits_table()
        P.type = []
        P.tilename = []
        #P.tileid = []
        P.filter = []
        P.exptime = []
        P.ra = []
        P.dec = []
        P.passnumber = []
        
        iahead = 0
        for ii,jplan in enumerate(J[iplan:]):
            if iahead >= Nahead:
                break
            tilename = str(jplan['object'])
            nextseq = seqnum + 1 + iahead
    
            print('Considering planning tile %s for exp %i' % (tilename, nextseq))
            
            # Check all planned tiles before this one for a duplicate tile.
            dup = False
            for s in range(nextseq-1, 0, -1):
                t = self.planned_tiles[s]
                if t == tilename:
                    dup = True
                    print('Wanted to plan tile %s (pass %i element %i) for exp %i'
                          % (tilename, nextpass, iplan+ii, nextseq),
                          'but it was already planned for exp %i' % s)
                    break
            if dup:
                continue
    
            self.planned_tiles[nextseq] = tilename
            iahead += 1
    
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
            airmass = GetAirmass(float(etile.alt))
            print('Airmass of planned tile:', airmass)
    
            if M['band'] == nextband:
                nextsky = skybright
            else:
                # Guess that the sky is as much brighter than canonical
                # in the next band as it is in this one!
                nextsky = ((skybright - nomsky) + self.gvs.sb_dict[nextband])
            
            fakezp = -99
            expfactor = ExposureFactor(nextband, airmass, ebv, seeing,
                                       fakezp, nextsky, self.gvs)
            # print('Exposure factor:', expfactor)
            exptime = expfactor * self.gvs.base_exptimes[nextband]
    
            ### HACK -- safety factor!
            print('Exposure time:', exptime)
            exptime *= 1.1
            exptime = int(np.ceil(exptime))
            print('Exposure time with safety factor:', exptime)
    
            exptime = np.clip(exptime, self.gvs.floor_exptimes[nextband],
                              self.gvs.ceil_exptimes[nextband])
            print('Clipped exptime', exptime)
            if nextband == 'z' and exptime > self.gvs.t_sat_max:
                exptime = self.gvs.t_sat_max
                print('Reduced exposure time to avoid z-band saturation:',
                      exptime)
            exptime = int(exptime)
    
            print('Changing exptime from', jplan['expTime'], 'to', exptime)
            jplan['expTime'] = exptime
    
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
    
            self.obs.date += (exptime + self.gvs.overheads) / 86400.
    
            print('%s: updated exposure %i to tile %s' %
                  (str(ephem.now()), nextseq, tilename))
    
            P.tilename.append(tilename)
            #P.tileid.append(tile.tileid)
            P.filter.append(nextband)
            P.exptime.append(exptime)
            P.ra.append(jplan['RA'])
            P.dec.append(jplan['dec'])
            P.passnumber.append(nextpass)
            P.type.append('P')
    
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
        for c in P.columns():
            print('col:', c, 'type', P.get(c).dtype)
    
        fn = 'mosbot-plan.fits'
        tmpfn = fn + '.tmp'
        P.writeto(tmpfn)
        os.rename(tmpfn, fn)
        print('Wrote', fn)
        
        return True
    

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
    main()
    
