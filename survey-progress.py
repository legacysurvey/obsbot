from __future__ import print_function
from glob import glob
import fitsio
from astrometry.util.fits import fits_table
import numpy as np

fns = glob('/project/projectdirs/cosmo/staging/mosaicz/MZLS_Raw/*/*_ori*.fits*')
fns.sort()
print(len(fns), 'files')

keys = ['EXPNUM', 'OBJECT', 'OBSTYPE', 'EXPTIME', 'RA', 'DEC', 'DATE-OBS',
        'TIME-OBS', 'MJD-OBS', 'AIRMASS', 'INSTRUME', 'FILTER']

defaults = [0, '', '', 0., '', '', '', '', 0., 0., '', '']

vals = dict([(k, []) for k in keys])
vals['RAW_FILENAME'] = []

for i,fn in enumerate(fns):
    print('Reading', i, 'of', len(fns), 'file', fn)
    hdr = fitsio.read_header(fn)
    for k,d in zip(keys, defaults):
        vals[k].append(hdr.get(k, d))
    vals['RAW_FILENAME'].append(fn)

    if i > 0 and i % 1000 == 0:
        T = fits_table()
        for k in keys:
            tk = k.lower()
            tk = tk.replace('-','_')
            T.set(tk, np.array(vals[k]))
        T.writeto('mzls-raw-part.fits')


T = fits_table()
for k in keys:
    tk = k.lower()
    tk = tk.replace('-','_')
    T.set(tk, np.array(vals[k]))
T.writeto('mzls-raw.fits')
