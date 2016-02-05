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
from measure_raw import measure_raw_decam

from jnox import *

ps = PlotSequence('bot')

parser,gvs = getParserAndGlobals()

parser.add_option('--plan', help='Use the given plan file (JSON) rather than computing a new plan.')
parser.add_option('--ext', help='Extension to read for computing observing conditions: default %default', default='N4')

opt,args = parser.parse_args()

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
    
obs = setupGlobals(opt, gvs)

#
if opt.plan is None:
    # the filename written by nightlystrategy.py
    jsonfn = 'mosaic_%s_plan.json' % opt.date

    # Generate nightly strategy.
    tiles,survey_centers = readTilesTable(opt.tiles, gvs)

    if 'G_DONE' not in survey_centers.keys():
        survey_centers['G_DONE'] = survey_centers['Z_DONE']
    if 'R_DONE' not in survey_centers.keys():
        survey_centers['R_DONE'] = survey_centers['Z_DONE']

    plan = GetNightlyStrategy(obs, opt.date, opt.portion, survey_centers,
                              opt.passnumber, gvs)
    WriteJSON(plan, jsonfn)

else:
    jsonfn = opt.plan
    
J = json.loads(open(jsonfn,'rb').read())

# Use jnox to write one shell script per exposure named by OBJECT (tile) name,
# plus one top-level script that runs them all.

# FIXME
scriptfn = jsonfn.replace('.json', '.sh')

script = [jnox_preamble(scriptfn)]
last_filter = None

fns = []
for i,j in enumerate(J):
    obj = j['object']
    fn = 'nocs-%s.sh' % obj
    fns.append(fn)

    next_obs = None
    if i < len(J)-1:
        next_obs = J[i+1]

    f = j['filter']
    if f != last_filter:
        output.write(jnox_filter(f))
        last_filter = f

    f = open(fn, 'w')
    f.write(jnox_cmds_for_json(j, i, len(J), next_obs=next_obs))
    f.close()
    print('Wrote', fn)
    script.append('. %s' % fn)

script = '\n'.join(script)

f = open(scriptfn, 'w')
f.write(script)
f.close()
print('Wrote', scriptfn)

tiles = fits_table(opt.tiles)

imagedir = 'rawdata'
lastimages = set(os.listdir(imagedir))

# Now we wait for the first image to appear (while we're taking the
# second exposure) and we use that to re-plan the third image.

for iwait,(jwait,jplan,planfn) in enumerate(zip(J, J[2:], fns[2:])):
    iplan = iwait + 2
    
    objwait = jwait['object']
    print('Waiting for tile %s to appear to re-plan tile %s' % (jwait['object'], jplan['object']))

    # Wait for the "jwait" image to appear.
    first = True
    while True:
        if not first:
            time.sleep(5.)
        first = False

        images = set(os.listdir(imagedir))
        newimgs = images - lastimages
        newimgs = list(newimgs)
        newimgs = [fn for fn in newimgs
                   if fn.endswith('.fits.fz') or fn.endswith('.fits')]
        print('Waiting for tile %s' % objwait)
        print('New images:', newimgs)
        if len(newimgs) == 0:
            continue

        newimg = newimgs[0]
        fn = os.path.join(imagedir, newimg)
        print('Found new file:', fn)
        try:
            hdr = fitsio.read_header(fn)
            print('Read header from', fn)
        except:
            print('Failed to open', fn, '-- maybe not fully written yet.')
            import traceback
            traceback.print_exc()
            continue

        lastimages.add(newimg)
        
        obj = hdr.get('OBJECT', None)
        print('Read OBJECT card:', obj)
        if obj == objwait:
            print('We found the image we wanted!')
            break

        if len(newimgs) > 1:
            # don't wait, try the remaining new ones right away
            first = True

    fn = os.path.join(imagedir, newimg)
    for i in range(10):
        try:
            ext = opt.ext
            if ext is None:
                ext = get_default_extension(fn)
            fitsio.read(fn, ext=ext)
        except:
            print('Failed to open', fn, '-- maybe not fully written yet.')
            import traceback
            traceback.print_exc()
            time.sleep(3)
            continue

    kwa = {}
    if opt.ext is not None:
        kwa.update(ext=opt.ext)

    M = measure_raw_decam(fn, ps=ps, **kwa)

    print('Measurements:', M)

    gvs.transparency = M['transparency']

    objname = jplan['object']
    # Parse objname like 'MzLS_5623_z'
    parts = objname.split('_')
    assert(len(parts) == 3)
    nextband = parts[2]
    assert(nextband in 'grz')
    tilenum = int(parts[1])

    print('Planned tile:', tilenum, nextband)

    # Set current date
    obs.date = ephem.now()

    present_tile = ephem.readdb(str(objname+','+'f'+','+('%.6f' % jplan['RA'])+','+
                                    ('%.6f' % jplan['dec'])+','+'20'))
    present_tile.compute(obs)
    airmass = GetAirmass(float(present_tile.alt))
    print('Airmass of planned tile:', airmass)

    # Find this tile in the tiles table.
    tile = tiles[tilenum-1]
    if tile.tileid != tilenum:
        I = np.flatnonzero(tile.tileid == tilename)
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
    print('Exptime (un-clipped)', exptime)
    exptime = np.clip(exptime, gvs.floor_exptimes[band], gvs.ceil_exptimes[band])
    print('Clipped exptime', exptime)
    
    if band == 'z' and exptime > gvs.t_sat_max:
        exptime = gvs.t_sat_max
        print('Reduced exposure time to avoid z-band saturation:', exptime)

    print

    print('Changing exptime in', planfn, 'from', jplan['expTime'],
          'to', exptime)
    # UPDATE EXPTIME!
    jplan['expTime'] = exptime
    
    next_obs = None
    if iplan < len(J)-1:
        next_obs = J[iplan+1]

    tmpfn = planfn + '.tmp'
    f = open(tmpfn, 'w')
    f.write(jnox_cmds_for_json(jplan, iplan, len(J), next_obs=next_obs))
    f.close()
    cmd = 'cp %s %s-orig', planfn, planfn
    print(cmd)
    os.system(cmd)
    os.rename(tmpfn, planfn)
    print('Wrote', planfn)

