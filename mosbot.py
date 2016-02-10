'''
TO-DO:

-- t_sat_max should be recomputed based on skybright!!
-- FIXMEs -- filter changing -- here we set filter ONCE PER NIGHT!

-- setting the filter was disabled for 2016-02-05 broken filter changer!

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
#from copilot import get_tile_from_name

def main():
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
    
    opt,args = parser.parse_args()
    
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
    
    # Use jnox to write one shell script per exposure
    # plus one top-level script that runs them all.
    
    scriptdir = os.path.dirname(opt.scriptfn)
    # Create scriptdir (one path component only), if necessary
    if len(scriptdir) and not os.path.exists(scriptdir):
        os.mkdir(scriptdir)
    
    seqnumfn = 'seqnum.txt'
    seqnumpath = os.path.join(scriptdir, seqnumfn)
    
    quitfile = 'quit'

    nowscriptpattern  = 'now-%i.sh'
    expscriptpattern  = 'expose-%i.sh'
    slewscriptpattern = 'slewread-%i.sh'

    # 775
    chmod = (stat.S_IXUSR | stat.S_IWUSR | stat.S_IRUSR |
             stat.S_IXGRP | stat.S_IWGRP | stat.S_IRGRP |
             stat.S_IXOTH |                stat.S_IROTH)
    
    # Default to Pass 2!
    J = [J1,J2,J3][opt.passnum - 1]

    planned_tiles = OrderedDict()

    def log(s):
        datecmd = 'date -u +"%Y-%m-%d %H:%M:%S"'
        #return 'echo "mosbot $(%s): %s"' % (datecmd, s)
        return 'echo "$(%s): %s" >> mosbot.log' % (datecmd, s)

    
    if opt.write_script:
        # Write "read.sh" for quitting gracefully without slew
        path = os.path.join(scriptdir, 'read.sh')
        f = open(path, 'w')
        f.write(jnox_readout() + '\n' + jnox_end_exposure())
        f.close()
        os.chmod(path, chmod)

        script = []
        script.append(jnox_preamble(opt.scriptfn))

        script.append(log('###############'))
        script.append(log('Starting script'))
        script.append(log('###############'))

        ## FIXME -- 2016-02-05 the filter changer is broken!
        ##j = J[0]
        ##f = j['filter']
        ##script.append(jnox_filter(f))
        
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
                fn = slewscriptpattern % seq
                path = os.path.join(scriptdir, fn)
                f = open(path, 'w')
                f.write(slewscript_for_json(j))
                f.close()
                os.chmod(path, chmod)
                print('Wrote', path)
                script.append('. %s' % fn)

            script.append('\n### Exposure %i ###\n' % seq)
            script.append('echo "%i" > %s' % (seq, seqnumfn))
            script.append(log('Starting exposure %i' % seq))
            expfn = expscriptpattern % seq
            script.append(log('Tile: $(grep "Tile:" %s)' % expfn))
            
            tilename = j['object']
            planned_tiles[seq] = tilename
            
            # Write now-##.sh
            fn = nowscriptpattern % seq
            path = os.path.join(scriptdir, fn)
            f = open(path, 'w')
            f.write('echo "Do stuff before exposure %i?"' % seq)
            f.close()
            os.chmod(path, chmod)
            print('Wrote', path)
            script.append('. %s' % fn)

            #tile = get_tile_from_name(tilename, tiles)
            #passnum = tile.get('pass') if tile is not None else 0

            ra  = ra2hms (j['RA' ])
            dec = dec2dms(j['dec'])
            status = ('Exp %i: Tile %s, Pass %i, RA %s, Dec %s' %
                      (seq, tilename, opt.passnum, ra, dec))

            if opt.exptime is not None:
                j['expTime'] = opt.exptime
            
            # Write expose-##.sh
            fn = expscriptpattern % seq
            path = os.path.join(scriptdir, fn)
            f = open(path, 'w')
            f.write(('# Exp %i: Tile: %s, Pass: %i; default\n' %
                     (seq, tilename, opt.passnum)) +
                     expscript_for_json(j, status=status))
            f.close()
            os.chmod(path, chmod)
            print('Wrote', path)
            script.append('. %s' % fn)
    
        # Read out the last one!
        script.append('. read.sh')
    
        script = '\n'.join(script) + '\n'
        
        f = open(opt.scriptfn, 'w')
        f.write(script)
        f.close()
        os.chmod(opt.scriptfn, chmod)
        print('Wrote', opt.scriptfn)
    
    imagedir = opt.rawdata
    if imagedir is None:
        imagedir = os.environ.get('MOS3_DATA', 'rawdata')

    print('Watching directory "%s" for new files' % imagedir)
    rawext = opt.ext

    lastimages = set(os.listdir(imagedir))
    
    # Now we wait for the first image to appear (while we're taking the
    # second exposure) and we use that to re-plan the third image.
    
    while True:
        first = True
        # Wait for a new image to appear
        while True:
            print
            if not first:
                time.sleep(5.)
            first = False

            print('Waiting for new images to appear in', imagedir, '...')
            
            images = set(os.listdir(imagedir))
            newimgs = images - lastimages
            newimgs = list(newimgs)
            newimgs = [fn for fn in newimgs if
                       fn.endswith('.fits.fz') or fn.endswith('.fits')]
            #print('Found new images:', newimgs)
            if len(newimgs) == 0:
                continue
    
            # Take the one with the latest timestamp.
            latest = None
            newestimg = None
            for fn in newimgs:
                st = os.stat(os.path.join(imagedir, fn))
                t = st.st_mtime
                if latest is None or t > latest:
                    newestimg = fn
                    latest = t
    
            fn = os.path.join(imagedir, newestimg)
            print('Found new file:', fn)
            try:
                print('Trying to open image:', fn, 'extension:', rawext)
                fitsio.read(fn, ext=rawext)
            except:
                print('Failed to open', fn, '-- maybe not fully written yet.')
                continue
            break
    
        try:
            ok = found_new_image(fn, rawext, opt, obs, gvs, seqnumpath,
                                 J1, J2, J3,
                                 os.path.join(scriptdir, expscriptpattern),
                                 os.path.join(scriptdir, slewscriptpattern),
                                 tiles, planned_tiles)
            if not ok:
                print('%s: found_new_image returned False' %
                      (str(ephem.now())))
                # Remove just this one image
                lastimages.add(newestimg)
            else:
                lastimages = lastimages.union(newimgs)
        except IOError:
            print('Failed to read FITS image:', fn, 'extension', rawext)
            import traceback
            traceback.print_exc()
            continue



def found_new_image(fn, ext, opt, obs, gvs, seqnumpath, J1, J2, J3,
                    expscriptpat, slewscriptpat, tiles, planned_tiles):

    print('%s: found new image %s' % (str(ephem.now()), fn))

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

    filt = phdr['FILTER']
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
    nomsky = gvs.sb_dict[M['band']]
    brighter = nomsky - skybright

    print('Transparency: %6.02f' % trans)
    print('Seeing      : %6.02f' % seeing)
    print('Sky         : %6.02f' % skybright)
    print('Nominal sky : %6.02f' % nomsky)
    print('Sky over nom: %6.02f   (positive means brighter than nom)' %
          brighter)
    
    nextpass = 3
    if trans > 0.90 and seeing < 1.25 and brighter < 0.25:
        nextpass = 1

    elif (trans > 0.90) or (seeing < 1.25):
        nextpass = 2

    print('Selected pass:', nextpass)

    # Choose the next tile from the right JSON tile list Jp
    J = [J1,J2,J3][nextpass-1]
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
    f = open(seqnumpath, 'r')
    s = f.read()
    f.close()
    seqnum = int(s)
    print('%s: sequence number: %i' % (str(ephem.now()), seqnum))
    
    # How many exposures ahead should we write?
    Nahead = 3

    # Set observing conditions for computing exposure time
    gvs.transparency = trans
    obs.date = now

    #for iahead,jplan in enumerate(J[iplan: iplan+Nahead]):

    iahead = 0
    for ii,jplan in enumerate(J[iplan:]):
        if iahead >= Nahead:
            break
        objname = jplan['object']
        nextseq = seqnum + 1 + iahead

        print('Considering planning tile %s for exp %i' % (objname, nextseq))
        
        # Check all planned tiles before this one for a duplicate tile.
        dup = False
        for s in range(nextseq-1, 0, -1):
            t = planned_tiles[s]
            if t == objname:
                dup = True
                print('Wanted to plan tile %s (pass %i element %i) for exp %i'
                      % (objname, nextpass, iplan+ii, nextseq),
                      'but it was already planned for exp %i' % s)
                break
        if dup:
            continue

        planned_tiles[nextseq] = objname
        iahead += 1

        # Parse objname like 'MzLS_5623_z'
        parts = objname.split('_')
        assert(len(parts) == 3)
        survey = parts[0]
        nextband = parts[2]
        tilenum = int(parts[1])
        print('Selected tile:', tilenum, nextband)
        rastr  = ra2hms (jplan['RA' ])
        decstr = dec2dms(jplan['dec'])
        ephemstr = str('%s,f,%s,%s,20' % (objname, rastr, decstr))
        etile = ephem.readdb(ephemstr)
        etile.compute(obs)
        airmass = GetAirmass(float(etile.alt))
        print('Airmass of planned tile:', airmass)

        # Find this tile in the tiles table.
        tile = tiles[tilenum-1]
        if tile.tileid != tilenum:
            I = np.flatnonzero(tiles.tileid == tilenum)
            assert(len(I) == 1)
            tile = tiles[I[0]]
        ebv = tile.ebv_med

        if M['band'] == nextband:
            nextsky = skybright
        else:
            # Guess that the sky is as much brighter than canonical
            # in the next band as it is in this one!
            nextsky = ((skybright - nomsky) + gvs.sb_dict[nextband])
        
        fakezp = -99
        expfactor = ExposureFactor(nextband, airmass, ebv, seeing,
                                   fakezp, nextsky, gvs)
        # print('Exposure factor:', expfactor)
        exptime = expfactor * gvs.base_exptimes[nextband]

        ### HACK -- safety factor!
        print('Exposure time:', exptime)
        exptime *= 1.1
        exptime = int(np.ceil(exptime))
        print('Exposure time with safety factor:', exptime)

        exptime = np.clip(exptime, gvs.floor_exptimes[nextband],
                          gvs.ceil_exptimes[nextband])
        print('Clipped exptime', exptime)
        if nextband == 'z' and exptime > gvs.t_sat_max:
            exptime = gvs.t_sat_max
            print('Reduced exposure time to avoid z-band saturation:',
                  exptime)
        exptime = int(exptime)

        print('Changing exptime from', jplan['expTime'], 'to', exptime)
        jplan['expTime'] = exptime

        tilename = jplan['object']
        ra  =  ra2hms(jplan['RA'])
        dec = dec2dms(jplan['dec'])
        status = ('Exp %i: Tile %s, Pass %i, RA %s, Dec %s' %
                  (nextseq, tilename, nextpass, ra, dec))

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

        obs.date += (exptime + gvs.overheads) / 86400.

        print('%s: updated exposure %i to tile %s' %
              (str(ephem.now()), nextseq, tilename))

    return True
    

def expscript_for_json(j, status=None):
    ra  = j['RA']
    dec = j['dec']
    ss = [jnox_moveto(ra, dec)]

    # exposure time
    et = (j['expTime'])
    # filter
    filter_name = j['filter']
    objname = j['object']
    if status is None:
        ra  = ra2hms(ra)
        dec = dec2dms(dec)
        status = "Tile %s, RA %s, Dec %s" % (objname, ra, dec)
    ss.append(jnox_exposure(et, filter_name, objname, status))
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
    
