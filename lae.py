from __future__ import print_function
import sys
import os
import json
import numpy as np
from collections import Counter
from astrometry.util.fits import fits_table
from astrometry.util.starutil_numpy import mjdtodate

fn = 'obs.fits'
if os.path.exists(fn):
    print('Reading existing copilot database', fn)
    print('To re-generate it, run')
    print('    python copilot.py --fits ', fn)
    print()
else:
    os.system('python copilot.py --no-show --fits %s' % fn)

if not os.path.exists('test_mLya.json'):
    print('Please create a symlink to the test-mLya.json file in the current directory (products/mosaic3/py/test_mLya.json)')
    sys.exit(-1)

T = fits_table(fn)
print(len(T), 'total observations')
I = np.flatnonzero((T.band == 'D51    ') * (T.obstype == 'object   ') * (T.exptime > 60))
T.cut(I)
print('Cut to', len(T), 'D51 object exposures > 60 seconds')
print()

J = json.loads(open('test_mLya.json').read())
print(len(J), 'mLya tiles defined')
obj_to_seqnum = dict([(j['object'], j['seqnum']) for j in J])

T.seqnum = np.array([obj_to_seqnum[o.strip()] for o in T.object])

print('Date (UTC)         Expnum Seeing Exptime Transparency Object   Seqnum')
for t in T:
    print('%-20s %i %6.2f %6.1f %6.2f %-16s %i' % (str(mjdtodate(t.mjd_obs))[:19], t.expnum, t.seeing, t.exptime, t.transparency, t.object.strip(), t.seqnum))

print()
T.good = ((T.seeing < 2.0) * (T.transparency > 0.4))
print(sum(T.good), 'exposures considered "good"')
tally = Counter(T.seqnum[T.good])

print()
print('Number of good exposures per seqnum:')
for seqnum in range(1, 58):
    print('%2i  %2i' % (seqnum, tally[seqnum]))





