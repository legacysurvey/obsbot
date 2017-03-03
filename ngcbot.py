from __future__ import print_function

import sys
import os
import re

import numpy as np
import pylab as plt

import fitsio

from measure_raw import measure_raw, get_default_extension, camera_name, DECamMeasurer
from obsbot import (exposure_factor, get_tile_from_name, NewFileWatcher,
                    mjdnow, datenow)

from astrometry.libkd.spherematch import *
from astrometry.util.plotutils import *
from astrometry.util.fits import *

ps = PlotSequence('ngcbot')

def main():
    from camera_decam import nominal_cal, data_env_var

    import optparse
    parser = optparse.OptionParser(usage='%prog')
    
    parser.add_option('--rawdata', help='Directory to monitor for new images: default $%s if set, else "rawdata"' % data_env_var, default=None)

    parser.add_option('--no-show', dest='show', default=True, action='store_false',
                      help='Do not show plot window, just save it.')
    
    opt,args = parser.parse_args()

    imagedir = opt.rawdata
    if imagedir is None:
        imagedir = os.environ.get(data_env_var, 'rawdata')

    nom = nominal_cal
    bot = NgcBot(imagedir, nom, opt)
    bot.run()
    return 0

class NgcBot(NewFileWatcher):
    def __init__(self, imagedir, nom, opt):
        super(NgcBot, self).__init__(imagedir, backlog=True)
        self.nom = nom
        self.opt = opt

        TT = []

        import astrometry
        fn = os.path.join(os.path.dirname(astrometry.__file__),
                          'catalogs', 'ngc2000.fits')
        print('Reading', fn)
        T = fits_table(fn)
        T.name = np.array(['NGC %i' % n for n in T.ngcnum])
        TT.append(T)
        
        fn = os.path.join(os.path.dirname(astrometry.__file__),
                          'catalogs', 'ic2000.fits')
        print('Reading', fn)
        T = fits_table(fn)
        T.name = np.array(['IC %i' % n for n in T.icnum])
        TT.append(T)
        self.cat = merge_tables(TT, columns='fillzero')

        
    def try_open_file(self, path):
        print('Trying to open file: %s' % path)
        fitsio.FITS(path)

    def process_file(self, path):
        print('Reading', path)
        F = fitsio.FITS(path)
        primhdr = F[0].read_header()
        print(len(F), 'extensions')

        band = primhdr['FILTER'][0]

        rr,dd = [],[]
        exts = np.arange(1, len(F))
        for i in exts:
            ext = i
            meas = DECamMeasurer(path, ext, self.nom)
            meas.primhdr = primhdr
            hdr = F[ext].read_header()
            wcs = meas.get_wcs(hdr)
            #print('WCS:', wcs)
            rc,dc = wcs.radec_center()
            radius = wcs.radius()
            #print('RA,Dec center', rc, dc, 'radius', radius)
            rr.append(rc)
            dd.append(dc)

        rr = np.array(rr)
        dd = np.array(dd)
        I,J,d = match_radec(rr, dd, self.cat.ra, self.cat.dec, radius)

        print('Matched', len(I), 'NGC objects')
        print(self.cat.name[J])

        cutouts = []
        
        for i,j in zip(I,J):
            ext = exts[i]
            r,d = self.cat.ra[j], self.cat.dec[j]
            #hdr = F[ext].read_header()
            zpt = self.nom.zeropoint(band, ext)
            print('Nominal zeropoint:', zpt)
            
            raw,hdr = meas.read_raw(F, ext)
            wcs = meas.get_wcs(hdr)

            median,hi = np.percentile(raw.ravel(), [50, 99])
            ok,x,y = wcs.radec2pixelxy(r, d)
            x = x - 1
            y = y - 1
            print('x,y', x,y)

            pixrad = self.cat.radius[j] * 3600. / wcs.pixel_scale()
            print('radius:', self.cat.radius[j], 'pixel radius:', pixrad)
            pixrad = max(pixrad, 100)
            r = pixrad
            tt = '%s in exp %i ext %s (%i)' % (self.cat.name[j], primhdr['EXPNUM'], hdr['EXTNAME'].strip(), ext)

            # plt.clf()
            # plt.imshow(raw, interpolation='nearest', origin='lower',
            #            vmin=median, vmax=hi, cmap='hot')
            # ax = plt.axis()
            # plt.plot([x-r, x-r, x+r, x+r, x-r], [y-r, y+r, y+r, y-r, y-r],
            #          'r-')
            # plt.title(tt)
            # ps.savefig()

            H,W = wcs.shape
            xl,xh = int(np.clip(x-r, 0, W-1)), int(np.clip(x+r, 0, W-1))
            yl,yh = int(np.clip(y-r, 0, H-1)), int(np.clip(y+r, 0, H-1))
            if xl == xh or yl == yh:
                continue
            subimg = raw[yl:yh, xl:xh]
            #print('Subimage size:', subimg.shape)
            lo,hi = np.percentile(subimg.ravel(), [50, 99.5])

            # plt.clf()
            # plt.imshow(subimg, interpolation='nearest', origin='lower',
            #            vmin=lo, vmax=hi, cmap='hot')
            # plt.title(tt)
            # plt.colorbar()
            # ps.savefig()

            lo,q3,hiX = np.percentile(subimg.ravel(), [25, 75, 99.9])
            mid = median

            def nlmap(x):
                Q = 10.
                return np.arcsinh((x-mid) / (q3-mid) * Q) / np.sqrt(Q)
            
            # plt.clf()
            # plt.imshow(nlmap(subimg), interpolation='nearest', origin='lower',
            #            vmin=nlmap(lo),
            #            cmap='hot')
            # # vmin=nlmap(lo),
            # # vmax=nlmap(mid + 100.*(q3-mid)),
            # # vmax=nlmap(hi),
            # plt.title(tt)
            # plt.colorbar()
            # ps.savefig()

            tt = '%s in %s' % (self.cat.name[j], hdr['EXTNAME'].strip())
            cutouts.append((tt, nlmap(subimg), nlmap(lo), None))

        plt.clf()
        NC = int(np.ceil(np.sqrt(len(cutouts)) * 1.3))
        NR = int(np.ceil(len(cutouts) / float(NC)))
        for i,(tt,img,lo,hi) in enumerate(cutouts):
            plt.subplot(NR, NC, i+1)
            plt.imshow(img, interpolation='nearest', origin='lower',
                       vmin=lo, vmax=hi, cmap='hot')
            plt.title(tt)
            plt.xticks([])
            plt.yticks([])
        plt.suptitle('Exposure %i' % (primhdr['EXPNUM']))

        if opt.show:
            plt.draw()
            plt.show(block=False)
            plt.pause(0.001)
        
        plt.savefig('cutouts.png')
        
if __name__ == '__main__':
    main()
            
