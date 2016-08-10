from __future__ import print_function

import os
from glob import glob

import matplotlib
matplotlib.use('Agg')
import pylab as plt

import numpy as np
import fitsio
from astrometry.util.fits import fits_table, merge_tables
from astrometry.util.util import wcs_pv2sip_hdr
from astrometry.util.plotutils import PlotSequence

if __name__ == '__main__':

    import astrometry

    ps = PlotSequence('pretty')

    catdir = os.path.join(os.path.dirname(astrometry.__file__),
                          'catalogs')
    fn = os.path.join(catdir, 'ngc2000.fits')
    NGC = fits_table(fn)
    NGC.name = np.array(['NGC %i' % ngcnum for ngcnum in NGC.ngcnum])
    #NGC.delete_column('ngcnum')

    IC = fits_table(os.path.join(catdir, 'ic2000.fits'))
    IC.name = np.array(['IC %i' % icnum for icnum in IC.icnum])
    #IC.delete_column('icnum')

    UGC = fits_table(os.path.join(catdir, 'ugc.fits'))
    UGC.name = np.array(['UGC %i' % ugcnum for ugcnum in UGC.ugcnum])
    #UGC.delete_column('ugcnum')

    UZC = fits_table(os.path.join(catdir, 'uzc2000.fits'))
    UZC.name = np.array(['UZC %s' % zname for zname in UZC.zname])
    
    ABELL = fits_table(os.path.join(catdir, 'abell-all.fits'))
    ABELL.name = np.array(['Abell %i' % aco for aco in ABELL.aco])

    cat = merge_tables([NGC, IC, UGC, UZC, ABELL],
                      columns='minimal')
    
    dirnm = os.environ['DECAM_DATA']
    fns = glob(dirnm + '/DECam_*.fits.fz')
    fns.sort()
    fns = list(reversed(fns))
    for ifn,fn in enumerate(fns):
        print('File', fn)
        F = fitsio.FITS(fn)
        primhdr = F[0].read_header()
        expnum = primhdr['EXPNUM']
        for ff in F[1:]:
            hdr = ff.read_header()
            extname = hdr['EXTNAME'].strip()
            if extname in ['S30', 'FN1','FN2', 'FN3', 'FN4', 'FS1', 'FS2', 'FS3', 'FS4']:
                continue

            wcs = wcs_pv2sip_hdr(hdr)
            ok,x,y = wcs.radec2pixelxy(cat.ra, cat.dec)

            sz = 200

            I = np.flatnonzero(ok * (x > sz) * (y > sz) *
                               (x < wcs.imagew-sz) * (y < wcs.imageh-sz))
            if len(I) == 0:
                continue
            print(len(I), 'catalog entries within WCS:', ', '.join(cat.name[I]))

            for i in I:
                xi,yi = int(x[i]), int(y[i])

                x0,x1 = xi-sz, xi+sz
                y0,y1 = yi-sz, yi+sz
                img = ff[y0:y1+1, x0:x1+1]
                mn,mx = np.percentile(img.ravel(), [25, 99])

                plt.clf()
                plt.imshow(img, interpolation='nearest', origin='lower',
                           vmin=mn, vmax=mx, cmap='gray', extent=[x0,x1,y0,y1])
                plt.title('%s in expnum %i ccd %s' % (cat.name[i], expnum, extname))
                ps.savefig()

            
        if ifn >= 5:
            break
    
