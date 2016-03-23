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

T = fits_table('Almanac_2016-03-20.fits')
T.about()

print('Extensions', np.unique(T.extname))
T.cut(T.extname == 'im16')
print(len(T), 'im16')
print('Filenames', T.filename)
print(' '.join(T.filename))
print('Expnums:', ', '.join(['%i' % e for e in T.expnum]))



fns = glob('/project/projectdirs/cosmo/staging/mosaicz/MZLS_Raw/20160320/*ori.fits.fz')
for fn in fns:
    hdr = fitsio.read_header(fn)
    expnum = hdr['EXPNUM']
    if not expnum in T.expnum:
        continue

    MM = []
    for ext in range(1, 16+1):
        extstring = 'im%i' % ext
        meas = measure_raw_mosaic3(fn, ext=extstring, n_fwhm=1)
        print('Measurement:', meas.keys())
        M = fits_table()
        for k in []:
            M.set(k, np.array([meas[k]]))
        MM.append(M)
    M = merge_tables(MM)
    
