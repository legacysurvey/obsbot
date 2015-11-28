from __future__ import print_function

from astrometry.util.plotutils import *

from RemoteClient import RemoteClient

from nightlystrategy import (ExposureFactor, getParserAndGlobals, setupGlobals,
                             GetAirmass, StartAndEndTimes, s_to_days)
from RemoteClient import RemoteClient

from measure_raw_decam import measure_raw_decam

if False:
    rc = RemoteClient()
    pid = rc.execute('get_propid')
    print('Got propid:', pid)



ps = PlotSequence('raw')

jsonfn = 'decals_2015-11-27_plan.json'
J = json.loads(open(jsonfn,'rb').read())

tiles = fits_table('decam-tiles_obstatus.fits')

date = '2015-11-17'
passnum = 1

j = J[0]
print('Planned tile:', j)

#M = measure_raw_decam('DECam_00488199.fits.fz')
M = {'seeing': 1.4890481099577366, 'airmass': 1.34,
     'skybright': 18.383479116033314, 'transparency': 0.94488537276869045,
     'band': 'z', 'zp': 26.442847814941093}

print('Measurements:', M)

#class Duck(object):
#    pass

parser,gvs = getParserAndGlobals()
# opt = Duck()
# opt.date = date
# opt.passnumber = passnum
# opt.portion = 1.0
# opt.start_date = '0000-00-00'
# opt.start_time = 

opt,args = parser.parse_args(('--date %s --pass %i --portion 1.0' %
                              (date, passnum)).split())

obs = setupGlobals(opt, gvs)

gvs.transparency = M['transparency']

objname = j['object']
# Parse objname like 'DECaLS_5623_z'
parts = objname.split('_')
assert(len(parts) == 3)
assert(parts[0] == 'DECaLS')
nextband = parts[2]
assert(nextband in 'grz')
tilenum = int(parts[1])

print('Next planned tile:', tilenum, nextband)

# Set current date
part = 1.0
(sn, en, lon, sn_18, en_18) = StartAndEndTimes(obs, part, gvs)
time_elapsed = 0.
obs.date = sn + time_elapsed*s_to_days

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


rc = RemoteClient()
pid = rc.execute('get_propid')
print('Got propid:', pid)

