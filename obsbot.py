from __future__ import print_function
import time
import json

import matplotlib
matplotlib.use('Agg')

import ephem

from astrometry.util.plotutils import *
from astrometry.util.fits import *

from RemoteClient import RemoteClient

from nightlystrategy import (ExposureFactor, getParserAndGlobals, setupGlobals,
                             GetAirmass, StartAndEndTimes, s_to_days)
from RemoteClient import RemoteClient

from measure_raw_decam import measure_raw_decam


M = measure_raw_decam('DECam_00488199.fits.fz')
import sys
sys.exit(0)

if False:
    rc = RemoteClient()
    pid = rc.execute('get_propid')
    print('Got propid:', pid)



ps = PlotSequence('raw')

jsonfn = 'decals_2015-11-27_plan.json'
J = json.loads(open(jsonfn,'rb').read())

tiles = fits_table('decam-tiles_obstatus.fits')

ext = 'N4'
date = '2015-11-17'
passnum = 1
parser,gvs = getParserAndGlobals()
opt,args = parser.parse_args(('--date %s --pass %i --portion 1.0' %
                              (date, passnum)).split())
obs = setupGlobals(opt, gvs)

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
            fitsio.read(fn, ext=ext)
        except:
            print('Failed to open', fn, '-- maybe not fully written yet.')
            import traceback
            traceback.print_exc()
            time.sleep(3)
            continue

    M = measure_raw_decam(fn, ext=ext, ps=ps)
    
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

