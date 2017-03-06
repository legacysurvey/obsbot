from __future__ import print_function

from collections import Counter

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

from legacypipe.survey import get_rgb

#rgbkwargs = dict(mnmx=(-1,100.), arcsinh=1.)
#rgbkwargs = dict(mnmx=(-1,1000.), arcsinh=1.)
rgbkwargs = dict(mnmx=(0,500.), arcsinh=1.)

ps = PlotSequence('ngcbot')

def main():
    from camera_decam import nominal_cal, data_env_var

    import optparse
    parser = optparse.OptionParser(usage='%prog')
    
    parser.add_option('--rawdata', help='Directory to monitor for new images: default $%s if set, else "rawdata"' % data_env_var, default=None)

    parser.add_option('--no-show', dest='show', default=True, action='store_false',
                      help='Do not show plot window, just save it.')

    parser.add_option('--tweet', default=False, action='store_true',
                      help='Send a tweet for each galaxy?')
    parser.add_option('--only', default=False, action='store_true',
                      help='Only run the command-line files, then quit')
    
    opt,args = parser.parse_args()

    imagedir = opt.rawdata
    if imagedir is None:
        imagedir = os.environ.get(data_env_var, 'rawdata')

    nom = nominal_cal
    bot = NgcBot(imagedir, nom, opt)

    if len(args):
        for a in args:
            bot.process_file(a)
        if opt.only:
            return 0

    bot.run()
    return 0


def sdss_rgb(rimgs, bands, scales=None,
             m = 0.02, clip=True):
    import numpy as np
    rgbscales = {'u': 1.5, #1.0,
                 'g': 2.5,
                 'r': 1.5,
                 'i': 1.0,
                 'z': 0.4, #0.3
                 }
    if scales is not None:
        rgbscales.update(scales)
        
    b,g,r = [rimg * rgbscales[b] for rimg,b in zip(rimgs, bands)]
    r = np.maximum(0, r + m)
    g = np.maximum(0, g + m)
    b = np.maximum(0, b + m)
    I = (r+g+b)/3.
    Q = 20
    fI = np.arcsinh(Q * I) / np.sqrt(Q)
    I += (I == 0.) * 1e-6
    R = fI * r / I
    G = fI * g / I
    B = fI * b / I
    # maxrgb = reduce(np.maximum, [R,G,B])
    # J = (maxrgb > 1.)
    # R[J] = R[J]/maxrgb[J]
    # G[J] = G[J]/maxrgb[J]
    # B[J] = B[J]/maxrgb[J]
    rgb = np.dstack((R,G,B))
    if clip:
        rgb = np.clip(rgb, 0, 1)
    return rgb
                



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
        T.delete_column('ngcnum')
        TT.append(T)
        
        fn = os.path.join(os.path.dirname(astrometry.__file__),
                          'catalogs', 'ic2000.fits')
        print('Reading', fn)
        T = fits_table(fn)
        T.delete_column('icnum')
        TT.append(T)
        self.cat = merge_tables(TT, columns='fillzero')

        print('Total of', len(T), 'NGC/IC objects')
        omit = [
            'OC', # open cluster
            'Ast', # asterism
            '***', # triple star
            'D*', # double star
            '*', # star
            '-', # called nonexistent in RNGC
            'PD' # plate defect
            ]
        # '?', # uncertain type or may not exist
        # '', # unidentified or type unknown
        print(Counter(T.classification))
        
        T.cut([t.strip() not in omit for t in T.classification])
        print('Cut to', len(T), 'NGC/IC objects')
        
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
        newband = band
        
        rr,dd = [],[]
        exts = np.arange(1, len(F))
        for i in exts:
            ext = i
            hdr = F[ext].read_header()
            extname = hdr['EXTNAME'].strip()
            if 'F' in extname:
                # Skip focus chips
                continue
            meas = DECamMeasurer(path, ext, self.nom)
            meas.primhdr = primhdr
            wcs = meas.get_wcs(hdr)
            #print('WCS:', wcs)
            rc,dc = wcs.radec_center()
            radius = wcs.radius()
            #print('RA,Dec center', rc, dc, 'radius', radius)
            rr.append(rc)
            dd.append(dc)

        # we use the last extension's 'radius' here...
        rr = np.array(rr)
        dd = np.array(dd)
        I,J,d = match_radec(rr, dd, self.cat.ra, self.cat.dec, radius)

        print('Matched', len(I), 'NGC objects')
        print(self.cat.name[J])

        tweets = []
        
        for i,j in zip(I,J):
            ext = exts[i]
            obj = self.cat[j]
            r,d = obj.ra, obj.dec
            hdr = F[ext].read_header()
            extname = hdr['EXTNAME'].strip()
            expnum = primhdr['EXPNUM']
            if 'F' in extname:
                # Skip focus chips
                continue

            wcs = meas.get_wcs(hdr)
            #print('Original WCS:', wcs)
            #print('WCS shape:', wcs.shape)
            ok,x,y = wcs.radec2pixelxy(obj.ra, obj.dec)
            #print('X,Y in original WCS:', x,y)
            x = x - 1
            y = y - 1

            pixrad = 1.4 * obj.radius * 3600. / wcs.pixel_scale()
            pixrad = max(pixrad, 100)
            print('radius:', obj.radius, 'pixel radius:', pixrad)
            r = pixrad
            tt = '%s in exp %i ext %s (%i)' % (obj.name, expnum, extname, ext)
            print(tt)

            H,W = wcs.shape
            xl,xh = int(np.clip(x-r, 0, W-1)), int(np.clip(x+r, 0, W-1))
            yl,yh = int(np.clip(y-r, 0, H-1)), int(np.clip(y+r, 0, H-1))
            if xl == xh or yl == yh:
                continue
            sh,sw = yh-yl, xh-xl
            if sh < 25 or sw < 25:
                continue
            #print('subimage:', xl, yl, xh-xl, yh-yl)
            #print('x range', xl, xh, 'y range', yl,yh)

            meas.ext = ext
            meas.edge_trim = 20
            M = meas.run(n_fwhm=1, verbose=False, get_image=True)
            #print('Measured:', M.keys())

            raw = M['image']
            #print('Raw image:', raw.shape)
            # Now repeat the cutout check with the trimmed image
            wcs = M['wcs']
            #print('WCS:', wcs)
            #print('WCS:', wcs.shape)
            # Trim WCS to trimmed raw image shape
            trim_x0, trim_y0 = M['trim_x0'], M['trim_y0']
            #print('Trim:', trim_x0, trim_y0)
            H,W = raw.shape
            wcs = wcs.get_subimage(trim_x0, trim_y0, W, H)
            #print('Trimmed WCS:', wcs.shape)
            #print('Trimmed WCS:', wcs)
            ok,x,y = wcs.radec2pixelxy(obj.ra, obj.dec)
            #print('X,Y in trimmed WCS:', x,y)
            x = x - 1
            y = y - 1
            xl,xh = int(np.clip(x-r, 0, W-1)), int(np.clip(x+r, 0, W-1))
            yl,yh = int(np.clip(y-r, 0, H-1)), int(np.clip(y+r, 0, H-1))
            #print('Trimmed X', xl, xh, 'Y', yl, yh)
            if xl == xh or yl == yh:
                continue
            subimg = raw[yl:yh, xl:xh]
            sh,sw = subimg.shape
            #print('Subimage shape', subimg.shape)
            if sh < 25 or sw < 25:
                continue
            subwcs = wcs.get_subimage(xl, yl, sw, sh)

            #dx = int(np.round(M['dx']))
            #dy = int(np.round(M['dy']))
            #print('DX,DY', dx,dy)
            dx = M['dx']
            dy = M['dy']
            
            aff = M['affine']
            x = (xl+xh)/2
            y = (yl+yh)/2
            #print('x,y', x,y)
            #print('Affine correction terms:', aff)
            adx = x - aff[0]
            ady = y - aff[1]
            corrx = aff[2] + aff[3] * adx + aff[4] * ady - adx
            corry = aff[5] + aff[6] * adx + aff[7] * ady - ady
            print('Affine correction', corrx, corry)

            ### Shift the 'subwcs' to account for astrometric offset
            cx,cy = subwcs.get_crpix()
            #subwcs.set_crpix((cx - dx, cy - dy))
            subwcs.set_crpix((cx - dx - corrx, cy - dy - corry))

            urlpat = 'http://legacysurvey.org/viewer/%s-cutout/?ra=%.4f&dec=%.4f&pixscale=%.3f&width=%i&height=%i&layer=%s'

            scale = 1.
            if max(sh,sw) > 1024:
                scale = 4.
            elif max(sh,sw) > 512:
                scale = 2.

            fitsimgs = []

            for layer in ['decals-dr3', 'sdssco']:
                rh,rw = int(np.ceil(sh/scale)),int(np.ceil(sw/scale))

                mx = max(rh, rw)
                rw = rh = mx
                
                # url = urlpat % ('jpeg', obj.ra, obj.dec, subwcs.pixel_scale() * scale, rw, rh, layer)
                # print('URL:', url)
                # r = requests.get(url)
                # 
                # ftmp = tempfile.NamedTemporaryFile()
                # ftmp.write(r.content)
                # ftmp.flush()
                # jpg = plt.imread(ftmp.name)
                # ftmp.close()
                # 
                # jpg = np.flipud(jpg)
                # #print('Got jpg', jpg.shape)
                # 
                # pixsc = subwcs.pixel_scale() * scale / 3600.
                # thiswcs = Tan(*[float(x) for x in
                #                 [obj.ra, obj.dec, 0.5 + rw/2., 0.5 + rh/2.,
                #                  -pixsc, 0., 0., pixsc, rw, rh]])
                # 
                # # print('Thiswcs shape:', thiswcs.shape)
                # 
                # try:
                #     Yo,Xo,Yi,Xi,rims = resample_with_wcs(subwcs, thiswcs)
                # except:
                #     continue
                # 
                # print('subwcs:', subwcs)
                # print('thiswcs:', thiswcs)
                # # 
                # # 
                # print('Yo', Yo.min(), Yo.max())
                # print('Xo', Xo.min(), Xo.max())
                # print('Yi', Yi.min(), Yi.max())
                # print('Xi', Xi.min(), Xi.max())
                # 
                # resamp = np.zeros((sh,sw,3), dtype=jpg.dtype)
                # 
                # print('resamp shape', resamp.shape)
                # print('vs subwcs.shape', subwcs.shape)
                # print('jpg shape', jpg.shape)
                # print('vs thiswcs.shape',thiswcs.shape)
                # 
                # for k in range(3):
                #     resamp[Yo,Xo,k] = jpg[Yi,Xi,k]

                url = urlpat % ('fits', obj.ra, obj.dec, subwcs.pixel_scale() * scale, rw, rh, layer)
                print('URL:', url)
                r = requests.get(url)

                # ftmp = tempfile.NamedTemporaryFile()
                # ftmp.write(r.content)
                # ftmp.flush()
                # fits,hdr = fitsio.read(ftmp.name, header=True)
                # ftmp.close()
                f,tmpfn = tempfile.mkstemp(suffix='.fits')
                os.write(f, r.content)
                os.close(f)
                fits,hdr = fitsio.read(tmpfn, header=True)
                print('Wrote FITS to', tmpfn)

                print('Got:', fits.shape)
                thiswcs = Tan(hdr)
                print('Bands:', hdr['BANDS'])
                print('WCS from FITS:', thiswcs)

                ### HACK -- surface brightness correction...
                if layer == 'sdssco':
                    s = (subwcs.pixel_scale() / 0.396)
                    fits *= s**2
                
                N,ww,hh = fits.shape
                imgs = [fits[n,:,:] for n in range(N)]
                bands = hdr['BANDS'].strip()

                # Resample the new image to this layer's WCS
                # FIXME -- we do this multiple times into all-assumed-the-same WCS.
                resamp3 = np.zeros((rh,rw), dtype=subimg.dtype)
                try:
                    #Yo,Xo,Yi,Xi,rims = resample_with_wcs(thiswcs, subwcs)
                    Yo,Xo,Yi,Xi,rims = resample_with_wcs(thiswcs, subwcs, [subimg])
                except:
                    continue
                #resamp3[Yo,Xo] = subimg[Yi,Xi]
                resamp3[Yo,Xo] = rims[0]

                if not np.all(fits == 0):
                    fitsimgs.append((layer, bands, imgs))

            if len(fitsimgs) == 0:
                return
    
            print()
            print('New image is', newband)
    
            plt.clf()
            plt.subplots_adjust(left=0.03, right=0.97, bottom=0.03)
            NC = 1 + len(fitsimgs)
            NR = 2
    
            zp = M['zp']
            zpscale = 10.**((zp - 22.5)/2.5)
            exptime = primhdr['EXPTIME']
            newimg = resamp3 / (zpscale * exptime)

            newimg_orig = newimg.copy()

            def my_rgb(imgs, bands, **kwargs):
                #return get_rgb(imgs, bands, **rgbkwargs)
                #return sdss_rgb(imgs, bands, scales=dict(g=6.0, r=3.4, i=2.5, z=2.2), m=0.03, **kwargs)
                return sdss_rgb(imgs, bands, scales=dict(g=6.0, r=3.4, i=2.5, z=2.2), m=-0.02, clip=False, **kwargs)
    
            #return sdss_rgb(rimgs, bands, 
            def grayscale(img, band):
                #oldrgb = my_rgb([img], [band])
                #oldrgb = my_rgb([img,img,img], [band,band,band])
                rgb = my_rgb([img,img,img], [band,band,band])
                #clip=False)
                index = 'zrg'.index(newband)
                gray = rgb[:,:,index]
                return gray
    
            targs = dict(fontsize=8)
            
            newgray = grayscale(newimg, newband)
    
            hi = np.percentile(newgray, 99.9)
            grayargs = dict(interpolation='nearest', origin='lower',
                            vmin=0., vmax=hi, cmap='gray')
            plt.subplot(NR, NC, 1)
            plt.imshow(newgray, **grayargs)
            plt.xticks([]); plt.yticks([])
            plt.title('New image (%s)' % newband, **targs)

            zeroimg = np.zeros_like(newimg)
            newimgs = [zeroimg, zeroimg, zeroimg]
            # sdss_rgb reverses the order, so do like grz.
            newbands = ['g','r','z']
            newindex = dict(g=0, r=1, i=2, z=2)
    
            j = newindex[newband]
            newimgs [j] = newimg
            newbands[j] = newband

            print('Setting band', newband, '(index %i)'%j, 'to new image')
            
            rgbs = []

            bestdata = None
            bestn = 0
            
            for i,(layer, bands, imgs) in enumerate(fitsimgs):

                origlayer = layer
                layer = {'decals-dr3': 'DECaLS DR3',
                         'sdssco': 'SDSS'}.get(layer, layer)
    
                # plt.subplot(NR, NC, 2+i)
                # plt.title(layer, **targs)
                # # empty white plot if not replaced
                # plt.imshow(np.ones_like(newimg), interpolation='nearest', origin='lower', vmin=0, vmax=1, cmap='hot')
                # plt.xticks([]); plt.yticks([])

                nicebands = []

                goodbands = []
                goodimgs = []
                
                for band,img in zip(bands, imgs):
                    if np.all(img == 0.):
                        print('Band', band, 'of', layer, 'is all zero')
                        nicebands.append('-')
                        continue
                    goodbands.append(band)
                    goodimgs.append(img)
                    nicebands.append(band)
                    print('  ', layer, band)
                    j = newindex[band]
                    if newimgs[j] is zeroimg:
                        print('    Setting band', band, '(index %i)'%j, 'to', layer)
                        newimgs[j] = img
                        newbands[j] = band
                    elif band == newbands[j]:
                        # patch empty regions if same band
                        print('    Patching index %i from' % j, layer, 'band', band)
                        Z = (newimgs[j] == 0)
                        newimgs[j][Z] = img[Z]

                    #if band == newband:
                    # z -> i
                    if newindex[band] == newindex[newband]:
                        # grab out grayscale
                        oldgray = grayscale(img, band)
                        plt.subplot(NR, NC, 2+i)
                        plt.imshow(oldgray, **grayargs)
                        plt.xticks([]); plt.yticks([])
                        plt.title('%s (%s)' % (layer, band), **targs)

                if len(goodbands) == 1:
                    #rgb = grayscale(goodimgs[0], goodbands[0])
                    img,band = goodimgs[0], goodbands[0]
                    rgb = my_rgb([img,img,img], [band,band,band])
                else:
                    rgb = my_rgb(imgs, bands)

                if len(goodbands) > bestn:
                    bestn = len(goodbands)
                    bestdata = origlayer

                nicebands = ''.join(nicebands)
                print('bands for', layer, ':', bands, ', actually', nicebands)
                rgbs.append((2+i+NC, rgb, '%s (%s)' % (layer, nicebands)))
    
            # list to string
            newbands = ''.join(newbands)
            print('Newbands:', newbands)
            
            plt.subplot(NR, NC, 1 + NC)
            rgb = my_rgb(newimgs, newbands)
            lo = 0.
            hi = np.percentile(rgb.ravel(), 99.9)

            rgb = np.clip((rgb - lo) / (hi - lo), 0., 1.)
            plt.imshow(rgb, interpolation='nearest', origin='lower')
            plt.xticks([]); plt.yticks([])
            plt.title('New+Old (%s)' % newbands, **targs)
    
            for sp, rgb, tt in rgbs:
                plt.subplot(NR, NC, sp)
                rgb = np.clip((rgb - lo) / (hi - lo), 0., 1.)
                plt.imshow(rgb, interpolation='nearest', origin='lower')
                plt.xticks([]); plt.yticks([])
                plt.title(tt, **targs)

            plt.suptitle('%s in DECam %i-%s: %s band' %
                         (obj.name, expnum, extname, newband))
    
            if self.opt.show:
                plt.draw()
                plt.show(block=False)
                plt.pause(0.001)

            plotfn = 'ngcbot-%i-%s-%s.png' % (expnum, extname, obj.name.replace(' ','_'))
            plt.savefig(plotfn)
            plt.savefig('ngcbot-latest.png')
            #plotfn = ps.pattern % (ps.format % ps.ploti, ps.suffixes[0])
            #ps.savefig()

            if self.opt.tweet:
                import urllib
                url = 'http://legacysurvey.org/viewer/?ra=%.3f&dec=%.3f' % (obj.ra, obj.dec)
                url2 = ('http://ned.ipac.caltech.edu/cgi-bin/objsearch?' + 
                         urllib.urlencode(dict(objname=obj.name, corr_z=1, list_limit=5, img_stamp='YES')))

                if bestdata is not None:
                    url += '&layer=%s' % bestdata

                dateobs = primhdr['DATE-OBS'].strip()
                print('Dateobs:', dateobs)
                # 2017-03-06T00:15:40.101482
                dateobs = dateobs[:19]
                dateobs = dateobs.replace('T', ' ')

                info = ''

                typenames = {
                    'Gx': 'galaxy',
                    'Gb': 'globular cluster',
                    'Nb': 'nebula',
                    'Pl': 'planetary nebula',
                    'C+N': 'cluster+nebulosity',
                    'Kt': 'knot in galaxy',
                }
                info = typenames.get(obj.classification.strip(), '')
                if len(info):
                    info = '(' + info + ') '
                
                txt = ('DECam observed %s %sin image %i-%s: %s band' %
                       (obj.name.strip(), info, expnum, extname, newband)
                       + ' at %s UT' % dateobs
                       + '\n' + url + '\n' + url2)
                print('Tweet text:', txt)

                nh,nw = newimg.shape
                coverage = np.sum(newimg_orig != 0) / float(nw*nh)
                print('Fraction', coverage, 'of new image is covered')

                if coverage > 0.75:
                    tweets.append((txt, plotfn))

        # Send one NGC object per exposure, chosen randomly.
        if self.opt.tweet and len(tweets):
            i = np.random.randint(0, len(tweets))
            txt,plotfn = tweets[i]
            send_tweet(txt, plotfn)

def send_tweet(txt, imgfn):
    from twython import Twython
    try:
        import secrets as s
    except:
        print('You need a "secrets.py" file with Twitter credentials')
        return

    twitter = Twython(s.APP_KEY, s.APP_SECRET,
                      s.OAUTH_TOKEN, s.OAUTH_TOKEN_SECRET)

    upload_response = twitter.upload_media(media=open(imgfn, 'rb'))
    media_id = upload_response['media_id']
    response = twitter.update_status(status=txt, media_ids=media_id)
    print(response)

if __name__ == '__main__':
    main()
            
