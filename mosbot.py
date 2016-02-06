from __future__ import print_function
import sys
import time
import json
import datetime

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

def main():
    
    #ps = PlotSequence('bot')
    
    parser,gvs = getParserAndGlobals()
    
    #parser.add_option('--plan', help='Use the given plan file (JSON) rather than computing a new plan.')

    parser.add_option('--rawdata', help='Directory to monitor for new images: default %default', default='rawdata')
    
    parser.add_option('--script', dest='scriptfn', help='Write top-level shell script, default is %default', default='tonight.sh')
    
    parser.add_option('--ext', help='Extension to read for computing observing conditions', default=None)
    
    #parser.add_option('--check-all', action='store_true', help='Check existing files for the image we are waiting for.')
    
    opt,args = parser.parse_args()
    
    if len(args) != 3:
        print('Usage: %s <pass1.json> <pass2.json> <pass3.json>' % sys.argv[0])
        sys.exit(-1)
    
        
    if opt.date is None:
        # Figure out the date.  Note that we have to figure out the date
        # at the START of the night.
        now = datetime.datetime.now()
        # Let's make a new day start at 9am, so subtract 9 hours from now
        nightstart = now - datetime.timedelta(0, 9 * 3600)
        d = nightstart.date()
        opt.date = '%04i-%02i-%02i' % (d.year, d.month, d.day)
        print('Setting date to', opt.date)
    
    # If start time was not specified, set it to NOW.
    default_time = '00:00:00'
    if opt.start_time == default_time and opt.portion is None:
        # Note that nightlystrategy.py takes UTC dates.
        now = datetime.datetime.utcnow()
        d = now.date()
        opt.start_date = '%04i-%02i-%02i' % (d.year, d.month, d.day)
        t = now.time()
        opt.start_time = t.strftime('%H:%M:%S')
        print('Setting start to %s %s UTC' % (opt.start_date, opt.start_time))
    
    if opt.portion is None:
        opt.portion = 1
    if opt.passnumber is None:
        opt.passnumber = '1'
        
    obs = setupGlobals(opt, gvs)
    
    #
    # if opt.plan is None:
    #     # the filename written by nightlystrategy.py
    #     jsonfn = 'mosaic_%s_plan.json' % opt.date
    # 
    #     # Generate nightly strategy.
    #     tiles,survey_centers = readTilesTable(opt.tiles, gvs)
    # 
    #     if 'G_DONE' not in survey_centers.keys():
    #         survey_centers['G_DONE'] = survey_centers['Z_DONE']
    #     if 'R_DONE' not in survey_centers.keys():
    #         survey_centers['R_DONE'] = survey_centers['Z_DONE']
    # 
    #     plan = GetNightlyStrategy(obs, opt.date, opt.portion, survey_centers,
    #                               opt.passnumber, gvs)
    #     WriteJSON(plan, jsonfn)
    # 
    # else:
    #     jsonfn = opt.plan
    
    json1fn,json2fn,json3fn = args
    
    J1 = json.loads(open(json1fn,'rb').read())
    J2 = json.loads(open(json2fn,'rb').read())
    J3 = json.loads(open(json3fn,'rb').read())
    
    # Use jnox to write one shell script per exposure named by OBJECT (tile) name,
    # plus one top-level script that runs them all.
    
    # FIXME
    #if opt.scriptfn is None:
    #    opt.scriptfn = jsonfn.replace('.json', '.sh')
    
    scriptdir = os.path.dirname(opt.scriptfn)
    
    script = [jnox_preamble(opt.scriptfn)]
    last_filter = None
    
    focus_elapsed = 0
    ifocus = 1
    
    # FIXMEs -- filter changing -- here we set filter ONCE PER NIGHT!
    
    seqnumfn = 'seqnum.txt'
    seqnumpath = os.path.join(scriptdir, seqnumfn)
    
    j = J1[0]
    f = j['filter']

    ## FIXME -- 2016-02-05 the filter changer is broken!
    ##script.append(jnox_filter(f))
    
    expscriptpattern = 'expose-%i.sh'
    slewscriptpattern = 'slewread-%i.sh'

    # Also write a "read.sh" for convenience after ctrl-C
    f = open(os.path.join(scriptdir, 'read.sh'), 'w')
    # ????
    f.write(jnox_readout() + '\n' + jnox_end_exposure())
    f.close()

    # Default to Pass 2!
    J = J2

    # Write top-level script and shell scripts for default plan.
    for i,j in enumerate(J):
    
        if i > 0:
            # e#slewread.sh script
            fn = slewscriptpattern % (i+1)
            path = os.path.join(scriptdir, fn)
            f = open(path, 'w')
            f.write(slewscript_for_json(j))
            f.close()
            print('Wrote', path)

            script.append('. %s' % fn)

        # e#exp.sh
        fn = expscriptpattern % (i+1)
        path = os.path.join(scriptdir, fn)
        f = open(path, 'w')
        f.write(expscript_for_json(j))
        f.close()
        print('Wrote', path)
    
        script.append('echo "%i" > %s' % (i+1, seqnumfn))
        script.append('. %s' % fn)

        focus_elapsed += j['expTime'] + gvs.overheads
        if focus_elapsed > 3600.:
            focus_elapsed = 0.
            focusfn = 'focus-%i.sh' % ifocus
            ifocus += 1
    
            script.append('. %s' % focusfn)
            path = os.path.join(scriptdir, focusfn)
    
            ##### FIXME -- need to get current FOCUS value!!
    
            focus_start = -8000
            foc = jnox_focus(5., focus_start)
            f = open(path, 'w')
            f.write(foc)
            f.close()
            print('Wrote', path)

    # Read out the last one!
    script.append('. read.sh')

    script = '\n'.join(script) + '\n'
    
    f = open(opt.scriptfn, 'w')
    f.write(script)
    f.close()
    print('Wrote', opt.scriptfn)
    
    print('Reading tiles table', opt.tiles)
    tiles = fits_table(opt.tiles)
    
    imagedir = opt.rawdata
    
    rawext = opt.ext
    #if opt.extnum is not None:
    #    rawext = opt.extnum
    
    # if opt.check_all:
    #     lastimages = set()
    # else:
    #     lastimages = set(os.listdir(imagedir))
    
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
    
            images = set(os.listdir(imagedir))
            newimgs = images - lastimages
            newimgs = list(newimgs)
            newimgs = [fn for fn in newimgs if
                       fn.endswith('.fits.fz') or fn.endswith('.fits')]
            print('Found new images:', newimgs)
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
            found_new_image(fn, rawext, opt, obs, gvs, seqnumpath, J1, J2, J3,
                            os.path.join(scriptdir, expscriptpattern),
                            os.path.join(scriptdir, slewscriptpattern),
                            tiles)
            lastimages.add(newestimg)
        except IOError:
            print('Failed to read FITS image:', fn, 'extension', rawext)
            import traceback
            traceback.print_exc()
            continue



def found_new_image(fn, ext, opt, obs, gvs, seqnumpath, J1, J2, J3,
                    expscriptpat, slewscriptpat, tiles):
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

    # Read the current sequence number
    f = open(seqnumpath, 'r')
    s = f.read()
    f.close()
    seqnum = int(s)
    print('Current sequence number', seqnum)
    
    # Measure the new image
    kwa = {}
    if ext is not None:
        kwa.update(ext=ext)
    ps = None
    M = measure_raw(fn, ps=ps, **kwa)


    # Sanity checks
    ok = (M['nmatched'] >= 20) and (M.get('zp',None) is not None)
    if not ok:
        # FIXME -- we could fall back to pass 3 here.
        return False
    
    # Choose the pass for seq+2
    trans = M['transparency']
    seeing = M['seeing']

    # eg, nominal = 20, sky = 19, brighter is 1 mag brighter than nom.
    nomsky = gvs.sb_dict[M['band']]
    brighter = nomsky - M['skybright']

    print('Transparency: %6.02f' % trans)
    print('Seeing      : %6.02f' % seeing)
    print('Sky         : %6.02f' % M['skybright'])
    print('Nominal sky : %6.02f' % nomsky)
    print('Sky over nom: %6.02f   (positive means brighter than nom)' %
          brighter)
    
    nextpass = 3
    if trans > 0.93 and seeing < 1.25 and brighter < 0.25:
        nextpass = 1

    elif (trans > 0.93) or (seeing < 1.25):
        nextpass = 2

    print('Selected pass:', nextpass)

    # Choose the next tile from the right JSON tile list Jp
    J = [J1,J2,J3][nextpass-1]
    now = ephem.now()
    jplan = None
    print('UTC now is', str(now))
    for jplan in J:
        tstart = ephem.Date(str(jplan['approx_datetime']))
        if tstart > now:
            print('Tile', jplan['object'], 'starts at', str(tstart))
            break
    if jplan is None:
        print('Could not find a JSON observation with approx_datetime after now =', str(now))
        return False

    # Compute exposure time for this tile.
    gvs.transparency = M['transparency']
    objname = jplan['object']
    # Parse objname like 'MzLS_5623_z'
    parts = objname.split('_')
    assert(len(parts) == 3)
    nextband = parts[2]
    assert(nextband in 'grz')
    tilenum = int(parts[1])
    print('Selected tile:', tilenum, nextband)
    # Set current date
    obs.date = now
    print('Observer lat,long', obs.lat, obs.lon)
    print('Date', obs.date)
    rastr  = ra2hms (jplan['RA' ])
    decstr = dec2dms(jplan['dec'])
    ephemstr = str('%s,f,%s,%s,20' % (objname, rastr, decstr))
    planned_tile = ephem.readdb(ephemstr)
    planned_tile.compute(obs)
    airmass = GetAirmass(float(planned_tile.alt))
    print('Airmass of planned tile:', airmass)
    
    # Find this tile in the tiles table.
    tile = tiles[tilenum-1]
    if tile.tileid != tilenum:
        I = np.flatnonzero(tiles.tileid == tilenum)
        print('Matched tiles:', I)
        assert(len(I) == 1)
        tile = tiles[I[0]]
    ebv = tile.ebv_med
    print('E(b-v) of planned tile:', ebv)

    if M['band'] == nextband:
        skybright = M['skybright']
    else:
        # Guess that the sky is as much brighter than canonical
        # in the next band as it is in this one!
        skybright = ((M['skybright'] - gvs.sb_dict[M['band']]) +
                     gvs.sb_dict[nextband])
    
    fakezp = -99
    expfactor = ExposureFactor(nextband, airmass, ebv, M['seeing'], fakezp,
                               skybright, gvs)
    print('Exposure factor:', expfactor)

    band = nextband
    exptime = expfactor * gvs.base_exptimes[band]


    ### HACK -- safety factor!
    print('Exposure time:', exptime)
    exptime *= 1.1
    exptime = int(np.ceil(exptime))
    print('Exposure time with safety factor:', exptime)

    exptime = np.clip(exptime, gvs.floor_exptimes[band], gvs.ceil_exptimes[band])
    print('Clipped exptime', exptime)
    if band == 'z' and exptime > gvs.t_sat_max:
        exptime = gvs.t_sat_max
        print('Reduced exposure time to avoid z-band saturation:', exptime)
    exptime = int(exptime)

    print('Changing exptime from', jplan['expTime'], 'to', exptime)

    # UPDATE EXPTIME!
    jplan['expTime'] = exptime

    expscriptfn = expscriptpat % (seqnum + 2)
    exptmpfn = expscriptfn + '.tmp'
    f = open(exptmpfn, 'w')
    f.write(expscript_for_json(jplan))
    f.close()

    slewscriptfn = slewscriptpat % (seqnum + 2)
    slewtmpfn = slewscriptfn + '.tmp'
    f = open(slewtmpfn, 'w')
    f.write(slewscript_for_json(jplan))
    f.close()

    cmd = 'cp %s %s-orig' % (expscriptfn, expscriptfn)
    print(cmd)
    os.system(cmd)
    cmd = 'cp %s %s-orig' % (slewscriptfn, slewscriptfn)
    print(cmd)
    os.system(cmd)

    os.rename(exptmpfn, expscriptfn)
    print('Wrote', expscriptfn)
    os.rename(slewtmpfn, slewscriptfn)
    print('Wrote', slewscriptfn)
    return True
    

def expscript_for_json(j):
    ra  = j['RA']
    dec = j['dec']
    ss = [jnox_moveto(ra, dec)]

    # exposure time
    et = (j['expTime'])
    # filter
    filter_name = j['filter']
    ra  = ra2hms(ra)
    dec = dec2dms(dec)
    objname = j['object']
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
    
