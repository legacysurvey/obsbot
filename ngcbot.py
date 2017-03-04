from __future__ import print_function

import sys
import os
import re

import requests
import tempfile

from PIL import Image
from io import BytesIO

import numpy as np
import pylab as plt

import fitsio

from measure_raw import measure_raw, get_default_extension, camera_name, DECamMeasurer
from obsbot import (exposure_factor, get_tile_from_name, NewFileWatcher,
                    mjdnow, datenow)

from astrometry.libkd.spherematch import *
from astrometry.util.plotutils import *
from astrometry.util.fits import *
from astrometry.util.util import *
from astrometry.util.resample import *

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

    if len(args):
        for a in args:
            bot.process_file(a)

    bot.run()
    return 0

class NgcBot(NewFileWatcher):
    def __init__(self, imagedir, nom, opt):
        super(NgcBot, self).__init__(imagedir, backlog=False)
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

        # import obsdb
        # from camera_decam import database_filename
        # obsdb.django_setup(database_filename=database_filename)
        # self.copilot_db = obsdb.MeasuredCCD.objects

        
    def try_open_file(self, path):
        print('Trying to open file: %s' % path)
        F = fitsio.FITS(path)
        ## HACK
        if len(F) < 60:
            print('Got', len(F), 'extensions')
            raise RuntimeError('fewer extensions than expected')

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
            obj = self.cat[j]
            r,d = obj.ra, obj.dec
            hdr = F[ext].read_header()
            extname = hdr['EXTNAME'].strip()
            if 'F' in extname:
                # Skip focus chips
                continue

            zpt = self.nom.zeropoint(band, ext)
            print('Nominal zeropoint:', zpt)
            
            raw,hdr = meas.read_raw(F, ext)
            wcs = meas.get_wcs(hdr)

            median,hi = np.percentile(raw.ravel(), [50, 99])
            ok,x,y = wcs.radec2pixelxy(r, d)
            x = x - 1
            y = y - 1
            # print('x,y', x,y)

            pixrad = obj.radius * 3600. / wcs.pixel_scale()
            pixrad = max(pixrad, 100)
            print('radius:', obj.radius, 'pixel radius:', pixrad)
            r = pixrad
            tt = '%s in exp %i ext %s (%i)' % (obj.name, primhdr['EXPNUM'], extname, ext)
            print(tt)
            
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
            sh,sw = subimg.shape
            if sh < 25 or sw < 25:
                continue
            print('Subimage size:', subimg.shape)

            subwcs = wcs.get_subimage(xl, yl, xh-xl, yh-yl)
            print('Subwcs:', subwcs.shape)

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

            #
            
            urlpat = 'http://legacysurvey.org/viewer/%s-cutout/?ra=%.4f&dec=%.4f&pixscale=%.3f&width=%i&height=%i&layer=%s'

            scale = 1.
            if max(sh,sw) > 1024:
                scale = 4.
            elif max(sh,sw) > 512:
                scale = 2.

            jpegs = []
            for layer in ['decals-dr3', 'sdssco']:
                rh,rw = int(np.ceil(sh/scale)),int(np.ceil(sw/scale))

                url = urlpat % ('jpeg', obj.ra, obj.dec, subwcs.pixel_scale() * scale, rw, rh, layer)
                print('URL:', url)
                r = requests.get(url)

                #jpg = Image.open(BytesIO(r.content))

                ftmp = tempfile.NamedTemporaryFile()
                ftmp.write(r.content)
                ftmp.flush()
                jpg = plt.imread(ftmp.name)
                ftmp.close()

                jpg = np.flipud(jpg)

                print('Got jpg', jpg.shape)
                #jpegs.append((layer,jpg))

                # url = urlpat % ('fits', obj.ra, obj.dec, subwcs.pixel_scale() * scale, rw, rh, layer)
                # print('URL:', url)
                # r = requests.get(url)
                # ftmp = tempfile.NamedTemporaryFile()
                # ftmp.write(r.content)
                # ftmp.flush()
                # fits,hdr = fitsio.read(ftmp.name, header=True)
                # ftmp.close()
                # print('Got:', fits.shape)
                # thiswcs = Tan(hdr)
                # print('Got WCS:', thiswcs)
                # print('Bands:', hdr['BANDS'])

                pixsc = subwcs.pixel_scale() * scale / 3600.
                thiswcs = Tan(*[float(x) for x in
                                [obj.ra, obj.dec, 0.5 + rw/2., 0.5 + rh/2.,
                                 -pixsc, 0., 0., pixsc, rw, rh]])

                # print('Thiswcs shape:', thiswcs.shape)
                
                try:
                    Yo,Xo,Yi,Xi,rims = resample_with_wcs(subwcs, thiswcs)
                except:
                    continue

                # print('subwcs:', subwcs)
                # print('thiswcs:', thiswcs)
                # 
                # print('subwcs RA,Dec bounds:', subwcs.radec_bounds())
                # print('thiswcs RA,Dec bounds:', thiswcs.radec_bounds())
                # 
                # print('Yo', Yo.min(), Yo.max())
                # print('Xo', Xo.min(), Xo.max())
                # print('Yi', Yi.min(), Yi.max())
                # print('Xi', Xi.min(), Xi.max())

                resamp = np.zeros((sh,sw,3), dtype=jpg.dtype)
                #print('Resamp:', resamp.shape, resamp.dtype)

                # print('resamp shape', resamp.shape)
                # print('vs subwcs.shape', subwcs.shape)
                # print('jpg shape', jpg.shape)
                # print('vs thiswcs.shape',thiswcs.shape)

                for k in range(3):
                    resamp[Yo,Xo,k] = jpg[Yi,Xi,k]

                jpegs.append((layer, jpg))

                jpegs.append((layer, resamp))
                
            # Has the copilot measured this exposure yet?

            meas.ext = ext
            M = meas.run(n_fwhm=1, verbose=False)
            
            #print('Measured:', M)

            dx = int(np.round(M['dx']))
            dy = int(np.round(M['dy']))
            print('DX,Dy', dx,dy)

            aff = M['affine']
            x = (xl+xh)/2
            y = (yl+yh)/2
            print('x,y', x,y)
            print('Affine correction terms:', aff)
            adx = x - aff[0]
            ady = y - aff[1]

            corrx = aff[2] + aff[3] * adx + aff[4] * ady - adx
            corry = aff[5] + aff[6] * adx + aff[7] * ady - ady

            print('Affine correction', corrx, corry)

            #subimg2 = raw[(yl+dy):(yh+dy), (xl+dx):(xh+dx)]
            subimg2 = raw[(yl-dy):(yh-dy), (xl-dx):(xh-dx)]
            print('subimg2:', subimg2.shape)

            # jpegs.append(('dx', subimg2))

            tt = '%s in %s' % (obj.name, extname)
            cutouts.append((tt, nlmap(subimg), nlmap(lo), None, nlmap(subimg2), jpegs))


        if len(cutouts) == 0:
            return
            
        # plt.clf()
        # NC = int(np.ceil(np.sqrt(len(cutouts)) * 1.3))
        # NR = int(np.ceil(len(cutouts) / float(NC)))
        # for i,(tt,img,lo,hi,jpegs) in enumerate(cutouts):
        #     plt.subplot(NR, NC, i+1)
        #     plt.imshow(img, interpolation='nearest', origin='lower',
        #                vmin=lo, vmax=hi, cmap='hot')
        #     plt.title(tt)
        #     plt.xticks([])
        #     plt.yticks([])
        # plt.suptitle('Exposure %i' % (primhdr['EXPNUM']))
        # 
        # if self.opt.show:
        #     plt.draw()
        #     plt.show(block=False)
        #     plt.pause(0.001)
        # 
        # plt.savefig('cutouts.png')

        plt.clf()
        NC = len(cutouts)
        NR = 6
        for i,(tt,img,lo,hi,img2, jpegs) in enumerate(cutouts):
            plt.subplot(NR, NC, i+1)
            plt.imshow(img, interpolation='nearest', origin='lower',
                       vmin=lo, vmax=hi, cmap='hot')
            plt.title(tt)
            plt.xticks([])
            plt.yticks([])

            plt.subplot(NR, NC, NC + i+1)
            plt.imshow(img2, interpolation='nearest', origin='lower',
                       vmin=lo, vmax=hi, cmap='hot')
            plt.title(tt)
            plt.xticks([])
            plt.yticks([])

            for j,(layer,pix) in enumerate(jpegs):
                plt.subplot(NR, NC, i+1 + (j+1+1)*NC)
                #print('pix:', pix.shape, pix.dtype)
                plt.imshow(pix, interpolation='nearest', origin='lower')
                #plt.title(layer)
                plt.xticks([])
                plt.yticks([])
        plt.suptitle('Exposure %i: %s band' % (primhdr['EXPNUM'], band))
        
        if self.opt.show:
            plt.draw()
            plt.show(block=False)
            plt.pause(0.001)
        
        plt.savefig('cutouts.png')

        
if __name__ == '__main__':
    main()
            
