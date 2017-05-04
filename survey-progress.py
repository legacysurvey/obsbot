from __future__ import print_function
from glob import glob
import fitsio
from astrometry.util.fits import fits_table
import numpy as np
import os
import pylab as plt
from astrometry.util.plotutils import *

outfn = 'mzls-raw.fits'
if not os.path.exists(outfn):
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



T = fits_table(outfn)

expnum = np.zeros(len(T), int)
for i,e in enumerate(T.expnum):
    e = e.strip()
    if len(e) == 0:
        continue
    expnum[i] = int(e)
T.expnum = expnum
    
print(len(T), 'exposures')

T.tileid = np.zeros(len(T), int)
for i,t in enumerate(T.object):
    t = t.strip()
    t = t.split('_')
    if len(t) != 3:
        #print(t)
        continue
    if t[0] != 'MzLS':
        print(t)
        continue
    assert(t[0] == 'MzLS')
    assert(t[2] == 'z')
    t = int(t[1])
    T.tileid[i] = t

tiles = fits_table('obstatus/mosaic-tiles_obstatus.fits')
print('', len(tiles), 'tiles')
tiles.rename('pass', 'passnum')
#print('passes:', np.unique(tiles.passnum))
tiles.cut(tiles.passnum <= 3)
print('', len(tiles), 'tiles with pass <= 3')
tiles.cut(np.flatnonzero(tiles.in_desi))
print('', len(tiles), 'in DESI')
# tiles.cut(tiles.dec >= 30)
# print(len(tiles), 'with Dec >= 30')
# tiles.cut(tiles.dec <= 80)
# print(len(tiles), 'with Dec <= 80')

tileid_to_index = dict(zip(tiles.tileid, np.arange(len(tiles))))

T.tileindex = np.array([tileid_to_index.get(t, -1) for t in T.tileid])

from collections import Counter
counts = Counter(T.tileindex[T.tileindex >= 0])

T.tilenobs = np.zeros(len(T), np.int16)
n = T.tilenobs
for i,t in enumerate(T.tileindex):
    n[i] = counts[t]

T.tilera   = np.zeros(len(T), tiles.ra.dtype)
T.tiledec  = np.zeros(len(T), tiles.dec.dtype)
T.tilepass = np.zeros(len(T), tiles.passnum.dtype)
I = (T.tileindex >= 0)
T.tilera[I]  = tiles.ra [T.tileindex[I]]
T.tiledec[I] = tiles.dec[T.tileindex[I]]
T.tilepass[I] = tiles.passnum[T.tileindex[I]]
T.writeto('mzls-raw-2.fits')

#nexp,xe,ye = np.histogram2d(T.tilera, T.tiledec)
#nexp = nexp.T

T.cut(T.expnum > 0)
print(len(T), 'with valid EXPNUM')
T.cut(np.array([f.startswith('zd ') for f in T.filter]))
print(len(T), 'with z filter')
T.cut(T.exptime >= 30.)
print(len(T), 'with exptime >= 30 sec')
T.writeto('ok.fits')
T.cut(T.tileid > 0)
print(len(T), 'with valid tileid in OBJECT.')

# Special fields: Stripe 82, COSMOS, etc.
#print('Out of bounds tiles:', T.tileid[tileindex == -1], 'RA', T.ra[tileindex == -1], 'Dec', T.dec[tileindex == -1])
T.cut(T.tileindex > -1)
print(len(T), 'exposures with matching valid tiles')

T.writeto('mzls-raw-ok.fits')

ps = PlotSequence('progress')

print('RA', T.tilera.min(), T.tilera.max())
print('Dec', T.tiledec.min(), T.tiledec.max())

ra0,ra1 = 85,305
dec0,dec1 = 30,85

ddec = 0.587

plt.figure(1, figsize=(12,6))
plt.subplots_adjust(left=0.1, right=0.99, bottom=0.1, top=0.95)

for p in [1,2,3]:
    I = np.flatnonzero(T.tilepass == p)

    plt.clf()
    #plt.
    plothist(T.tilera[I], T.tiledec[I], #(220, 55), range=((85, 305),(30, 85)),
            (np.arange(ra0, ra1, ddec*2), np.arange(dec0, dec1, ddec)),
             imshowargs=dict(cmap='viridis', vmin=0, vmax=5), dohot=False)
    plt.title('MzLS pass %i' % p)
    plt.ylabel('Dec (deg)')
    plt.xticks(np.arange(90, ra1+30, 30))
    plt.xlim(ra1,ra0)
    plt.ylim(dec0, dec1)
    plt.xlabel('RA (deg)')
    ps.savefig()

plt.figure(2, figsize=(12,4))
plt.subplots_adjust(left=0.05, right=0.99, bottom=0.1, top=0.95)
# plt.clf()

for p in [1,2,3]:
    I = np.flatnonzero(T.tilepass == p)

    ti = np.flatnonzero(tiles.passnum == p)

    #plt.subplot(3,1,p)
    plt.clf()
    plt.plot(tiles.ra[ti], tiles.dec[ti], 'k.', mec='none', alpha=0.05)
    plt.scatter(T.tilera[I], T.tiledec[I],
                alpha=0.1)
    plt.title('MzLS pass %i' % p)
    plt.ylabel('Dec (deg)')
    plt.xticks(np.arange(90, ra1+30, 30))
    plt.axis([ra1,ra0,dec0,dec1])
    plt.xlabel('RA (deg)')
    plt.axis('scaled')
    plt.axis([ra1,ra0,dec0,dec1])
    ps.savefig()
#ps.savefig()


    
