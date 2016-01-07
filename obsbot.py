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

from RemoteClient import RemoteClient

from nightlystrategy import (
    ExposureFactor, getParserAndGlobals, setupGlobals,
    GetAirmass, StartAndEndTimes, s_to_days, readTilesTable, GetNightlyStrategy,
    WriteJSON)
from RemoteClient import RemoteClient

from measure_raw_decam import measure_raw_decam


if False:
    M = measure_raw_decam('DECam_00488199.fits.fz')
    sys.exit(0)

if False:
    rc = RemoteClient()
    pid = rc.execute('get_propid')
    print('Got propid:', pid)

ps = PlotSequence('raw')



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
    jsonfn = 'decals_%s_plan.json' % opt.date

    # Generate nightly strategy.
    tiles,survey_centers = readTilesTable(opt.tiles, gvs)
    plan = GetNightlyStrategy(obs, opt.date, opt.portion, survey_centers,
                              opt.passnumber, gvs)
    WriteJSON(plan, jsonfn)

else:
    jsonfn = opt.plan
    
J = json.loads(open(jsonfn,'rb').read())

tiles = fits_table(opt.tiles)

imagedir = 'rawdata'
lastimages = set(os.listdir(imagedir))

for j in J:
    print('Planned tile:', j)

    # Wait for a new image to appear
    while True:
        images = set(os.listdir(imagedir))
        newimgs = images - lastimages
        newimgs = list(newimgs)
        newimgs = [fn for fn in newimgs if fn.endswith('.fits.fz')]
        print('New images:', newimgs)
        if len(newimgs) == 0:
            time.sleep(5.)
            continue
        lastimages = images
        break

    fn = os.path.join(imagedir, newimgs[0])
    for i in range(10):
        try:
            fitsio.read(fn, ext=opt.ext)
        except:
            print('Failed to open', fn, '-- maybe not fully written yet.')
            import traceback
            traceback.print_exc()
            time.sleep(3)
            continue

    M = measure_raw_decam(fn, ext=opt.ext, ps=ps)
    
    #M = measure_raw_decam('DECam_00488199.fits.fz')
    #M = {'seeing': 1.4890481099577366, 'airmass': 1.34,
    #'skybright': 18.383479116033314, 'transparency': 0.94488537276869045,
    #'band': 'z', 'zp': 26.442847814941093}

    print('Measurements:', M)

    gvs.transparency = M['transparency']

    objname = j['object']
    # Parse objname like 'DECaLS_5623_z'
    parts = objname.split('_')
    assert(len(parts) == 3)
    assert(parts[0] == 'DECaLS')
    nextband = parts[2]
    assert(nextband in 'grz')
    tilenum = int(parts[1])

    print('Planned tile:', tilenum, nextband)

    # Set current date
    obs.date = ephem.now()

    present_tile = ephem.readdb(str(objname+','+'f'+','+('%.6f' % j['RA'])+','+
                                    ('%.6f' % j['dec'])+','+'20'))
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

        
# rc = RemoteClient()
# pid = rc.execute('get_propid')
# print('Got propid:', pid)

