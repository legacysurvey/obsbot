from __future__ import print_function
import fitsio
from astrometry.blind.plotstuff import Plotstuff
from astrometry.util.util import Sip, anwcs, anwcs_new_sip, wcs_pv2sip_hdr

def skyplot():
    plot = Plotstuff(size=(800,800), rdw=(103.1, 37.45, 0.8), outformat='png')
    plot.color = 'verydarkblue'
    plot.plot('fill')
    
    for ext in range(1, 17):
        fn = 'mos3.68488.fits'
        hdr = fitsio.read_header(fn, ext=ext)
        wcs = wcs_pv2sip_hdr(hdr)
        plot.color = 'red'
        plot.outline.wcs = anwcs_new_sip(wcs)
        plot.plot('outline')
    
        plot.color = 'white'
        plot.apply_settings()
        rc,dc = wcs.radec_center()
        plot.text_radec(rc, dc, hdr['EXTNAME'])
        
    plot.color = 'white'
    for ext in range(1, 17):
        fn = 'an2/mos3.68488.ext%02i.wcs' % ext
        plot.outline.wcs = anwcs(fn)
        plot.plot('outline')
    
    plot.rgb = (0.2,0.2,0.2)
    plot.plot_grid(0.1, 0.1)
    plot.color = 'gray'
    plot.plot_grid(0.5, 0.5, 0.5, 0.5)
    plot.write('plot.png')


from astrometry.util.fits import fits_table, merge_tables
import numpy as np
from glob import glob
from measure_raw import measure_raw_mosaic3
import os

T = fits_table('Almanac_2016-03-20.fits')
T.about()

#print('Extensions', np.unique(T.extname))
T.cut(T.extname == 'im16')
print(len(T), 'im16')
#print('Filenames', T.filename)
#print(' '.join(T.filename))
#print('Expnums:', ', '.join(['%i' % e for e in T.expnum]))


cofn = 'copilot.fits'
if not os.path.exists(cofn):
    expnum_map = {}
    fns = glob('/project/projectdirs/cosmo/staging/mosaicz/MZLS_Raw/20160320/*ori.fits.fz')
    for fn in fns:
        hdr = fitsio.read_header(fn)
        expnum = hdr['EXPNUM']
        if not expnum in T.expnum:
            continue
        expnum_map[expnum] = fn
        print('File', fn, 'is expnum', expnum)
        
    MM = []
    for i,t in enumerate(T):
        #extstring = t.extname
        fn = expnum_map[t.expnum]
    
        for ext in range(1, 16+1):
            extstring = 'im%i' % ext
            meas = measure_raw_mosaic3(fn, ext=extstring, n_fwhm=1)
            #print('Measurement:', meas.keys())
            M = fits_table()
            for k in ['airmass', 'extension', 'pixscale', 'nmatched', 'ra_ccd', 'dec_ccd',
                      'band', 'zp', 'rawsky', 'ndetected', 'skybright', 'dy', 'transparency',
                      'seeing', 'dx', 'exptime', 'zp_skysub', 'zp_med', 'zp_med_skysub']:
                M.set(k, np.array([meas[k]]))
            M.filename = np.array([fn])
            M.extname = np.array([extstring])
            phdr = M['primhdr']
            M.expnum = np.array([phdr['EXPNUM']])
            MM.append(M)
    M = merge_tables(MM)

    M.writeto(cofn)

C = fits_table(cofn)

if not 'expnum' in C.columns():
    fns = glob('/project/projectdirs/cosmo/staging/mosaicz/MZLS_Raw/20160320/*ori.fits.fz')
    expnum_map = {}
    for fn in fns:
        hdr = fitsio.read_header(fn)
        expnum = hdr['EXPNUM']
        if not expnum in T.expnum:
            continue
        expnum_map[fn] = expnum
        print('File', fn, 'is expnum', expnum)
    C.expnum = np.array([expnum_map[fn.strip()] for fn in C.filename])
    C.writeto('copilot2.fits')
    
TC = fits_table()
for c in C.columns():
    TC.set(c, [])
        
# First: mosstat vs copilot for im16.
for i in range(len(T)):
    J = np.flatnonzero((C.expnum == T.expnum[i]) * (C.extension == 'im16'))
    assert(len(J) == 1)
    j = J[0]
    for c in C.columns():
        TC.get(c).append(C.get(c)[j])
TC.to_np_arrays()
T.add_columns_from(TC)

import pylab as plt

plt.clf()
plt.subplot(1,2,1)
plt.plot(T.dx, T.ra_offset, 'b.')
plt.xlabel('dx')
plt.ylabel('RA offset')
    
plt.subplot(1,2,2)
plt.plot(T.dy, T.dec_offset, 'b.')
plt.xlabel('dy')
plt.ylabel('Dec offset')
plt.savefig('diff.png')
