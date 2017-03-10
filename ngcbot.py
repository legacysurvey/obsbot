from __future__ import print_function

from collections import Counter
import sys
import os
import tempfile
import requests
import numpy as np
import fitsio
#import pylab as plt
plt = None

from astrometry.libkd.spherematch import match_radec
from astrometry.util.fits import fits_table, merge_tables
from astrometry.util.util import Tan
from astrometry.util.resample import resample_with_wcs

# Mosaic or Decam is selected via a symlink:
#   camera.py -> camera_decam.py    OR
#   camera.py -> camera_mosaic.py
from camera import nominal_cal, nice_camera_name, data_env_var, min_n_exts

from measure_raw import get_measurer_class_for_file
from obsbot import NewFileWatcher

if nice_camera_name == 'DECam':
    legacy_survey_layers = ['decals-dr3', 'sdssco']
else:
    legacy_survey_layers = ['mobo-dr4', 'sdssco']

nicelayernames = { 'decals-dr3': 'DECaLS DR3',
                   'mobo-dr4': 'MzLS+BASS DR4',
                   'sdssco': 'SDSS'}

ngc_typenames = {
    'Gx': 'galaxy',
    'Gb': 'globular cluster',
    'Nb': 'nebula',
    'Pl': 'planetary nebula',
    'C+N': 'cluster+nebulosity',
    'Kt': 'knot in galaxy',
}

ps = None

def main():
    import optparse
    parser = optparse.OptionParser(usage='%prog')

    parser.add_option('--rawdata', default=None,
                      help=('Directory to monitor for new images: default ' +
                            '$%s if set, else "rawdata"' % data_env_var))

    parser.add_option('--no-show', dest='show', default=True,
                      action='store_false',
                      help='Do not show plot window, just save it.')

    parser.add_option('--plotdir', help='Save plots in this directory',
                      default='ngcbot-plots')

    parser.add_option('--tweet', default=False, action='store_true',
                      help='Send a tweet for each galaxy?')
    parser.add_option('--only', default=False, action='store_true',
                      help='Only run the command-line files, then quit')

    opt,args = parser.parse_args()

    global plt
    if not opt.show:
        import matplotlib
        matplotlib.use('Agg')
    import pylab
    plt = pylab
    from astrometry.util.plotutils import PlotSequence
    global ps
    ps = PlotSequence('ngcbot')

    imagedir = opt.rawdata
    if imagedir is None:
        imagedir = os.environ.get(data_env_var, 'rawdata')

    if opt.plotdir and not os.path.exists(opt.plotdir):
        os.makedirs(opt.plotdir)

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

        # Read catalogs
        TT = []
        import astrometry
        catdir = os.path.join(os.path.dirname(astrometry.__file__),
                              'catalogs')
        fn = os.path.join(catdir, 'ngc2000.fits')
        print('Reading', fn)
        T = fits_table(fn)
        T.delete_column('ngcnum')
        TT.append(T)

        fn = os.path.join(catdir, 'ic2000.fits')
        print('Reading', fn)
        T = fits_table(fn)
        T.delete_column('icnum')
        TT.append(T)
        self.cat = merge_tables(TT, columns='fillzero')
        del TT

        print('Total of', len(self.cat), 'NGC/IC objects')
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
        print(Counter(self.cat.classification))

        self.cat.cut(np.flatnonzero(np.array([t.strip() not in omit for t in self.cat.classification])))
        print('Cut to', len(self.cat), 'NGC/IC objects')

        print('Remaining classifications:', np.unique(self.cat.classification))

    def try_open_file(self, path):
        print('Trying to open file: %s' % path)
        F = fitsio.FITS(path)
        ## HACK
        if len(F) < min_n_exts:
            print('Got', len(F), 'extensions')
            raise RuntimeError('fewer extensions than expected')

    def process_file(self, path):
        print('Reading', path)
        F = fitsio.FITS(path)
        primhdr = F[0].read_header()
        print(len(F), 'extensions')
        # measure_raw . DECamMeasurer or Mosaic3Measurer
        meas_class = get_measurer_class_for_file(path)
        if meas_class is None:
            print('Failed to identify camera in', path)
            return
        meas = meas_class(path, 0, self.nom)

        # We read the WCS headers from all extensions and then spherematch
        # against the catalog.
        rr,dd = [],[]
        exts = np.arange(1, len(F))
        keep_exts = []
        for ext in exts:
            hdr = F[ext].read_header()
            extname = hdr['EXTNAME'].strip()
            # HACK -- skip DECam focus chips.
            if 'F' in extname:
                continue
            meas.ext = ext
            meas.primhdr = primhdr
            wcs = meas.get_wcs(hdr)
            #print('WCS:', wcs)
            rc,dc = wcs.radec_center()
            radius = wcs.radius()
            #print('RA,Dec center', rc, dc, 'radius', radius)
            rr.append(rc)
            dd.append(dc)
            keep_exts.append(ext)
        exts = keep_exts
        # we use the last extension's 'radius' here...
        rr = np.array(rr)
        dd = np.array(dd)
        I,J,d = match_radec(rr, dd, self.cat.ra, self.cat.dec, radius)

        # plt.clf()
        # angles = np.linspace(0, 2.*np.pi, 40)
        # for r,d in zip(rr, dd):
        #     plt.plot(r + radius * np.sin(angles) / np.cos(np.deg2rad(d)), d + radius * np.cos(angles), 'b-')
        # plt.plot(rr, dd, 'bo')
        # ax = plt.axis()
        # plt.plot(self.cat.ra, self.cat.dec, 'r.')
        # plt.axis(ax)
        # ps.savefig()

        print('Matched', len(I), 'NGC objects')
        if len(I) == 0:
            return
        for j in J:
            info = ngc_typenames.get(self.cat.classification[j].strip(), '')
            print('  ', self.cat.name[j], info)

        # Potential tweet texts and plot filenames
        tweets = []

        for i,j in zip(I,J):
            ext = exts[i]
            obj = self.cat[j]
            obj.name = obj.name.strip()
            hdr = F[ext].read_header()
            extname = hdr['EXTNAME'].strip()
            expnum = primhdr['EXPNUM']
            wcs = meas.get_wcs(hdr)
            ok,x,y = wcs.radec2pixelxy(obj.ra, obj.dec)
            x = x - 1
            y = y - 1

            # Choose cutout area
            pixrad = 1.4 * obj.radius * 3600. / wcs.pixel_scale()
            pixrad = max(pixrad, 100)
            # print('radius:', obj.radius, 'pixel radius:', pixrad)
            r = pixrad
            tt = '%s in exp %i ext %s (%i)' % (obj.name, expnum, extname, ext)
            print(tt)

            # Find the cutout region... does it actually overlap the chip?
            H,W = wcs.shape
            xl,xh = int(np.clip(x-r, 0, W-1)), int(np.clip(x+r, 0, W-1))
            yl,yh = int(np.clip(y-r, 0, H-1)), int(np.clip(y+r, 0, H-1))
            if xl == xh or yl == yh:
                print('no actual overlap with image')
                continue
            sh,sw = yh-yl, xh-xl
            if sh < 25 or sw < 25:
                print('tiny overlap', sw, 'x', sh)
                continue

            # Measure the image!
            meas.ext = ext
            meas.edge_trim = 20
            try:
                M = meas.run(n_fwhm=1, verbose=False, get_image=True)
            except KeyboardInterrupt:
                sys.exit(0)
            except:
                import traceback
                print('Failed to measure file', path, 'ext', ext, ':')
                traceback.print_exc()
                continue
            #print('Measured:', M.keys())
            raw = M['image']

            # Now repeat the cutout check with the trimmed image
            wcs = M['wcs']
            # Trim WCS to trimmed raw image shape
            trim_x0, trim_y0 = M['trim_x0'], M['trim_y0']
            H,W = raw.shape
            wcs = wcs.get_subimage(trim_x0, trim_y0, W, H)
            ok,x,y = wcs.radec2pixelxy(obj.ra, obj.dec)
            x = x - 1
            y = y - 1
            xl,xh = int(np.clip(x-r, 0, W-1)), int(np.clip(x+r, 0, W-1))
            yl,yh = int(np.clip(y-r, 0, H-1)), int(np.clip(y+r, 0, H-1))
            if xl == xh or yl == yh:
                print('no actual overlap with image')
                continue
            subimg = raw[yl:yh, xl:xh]
            sh,sw = subimg.shape
            if sh < 25 or sw < 25:
                print('tiny overlap', sw, 'x', sh)
                continue
            subwcs = wcs.get_subimage(xl, yl, sw, sh)

            # Astrometric shifts
            dx = M['dx']
            dy = M['dy']
            aff = M['affine']
            x = (xl+xh)/2
            y = (yl+yh)/2
            #print('Affine correction terms:', aff)
            adx = x - aff[0]
            ady = y - aff[1]
            corrx = aff[2] + aff[3] * adx + aff[4] * ady - adx
            corry = aff[5] + aff[6] * adx + aff[7] * ady - ady
            print('Affine correction', corrx, corry)
            # Shift the 'subwcs' to account for astrometric offset
            cx,cy = subwcs.get_crpix()
            subwcs.set_crpix((cx - dx - corrx, cy - dy - corry))

            # Now grab existing data from the LegacySurvey site

            # What size of image are we going to request?
            scale = 1.
            if max(sh,sw) > 1024:
                scale = 4.
            elif max(sh,sw) > 512:
                scale = 2.
            rh,rw = int(np.ceil(sh/scale)),int(np.ceil(sw/scale))
            # make it square
            mx = max(rh, rw)
            rw = rh = mx

            fitsimgs = []
            # We'll resample the new image into the existing-image WCS.
            newimg = None

            for layer in legacy_survey_layers:

                url = ('http://legacysurvey.org/viewer-dev/fits-cutout/?ra=%.4f&dec=%.4f&pixscale=%.3f&width=%i&height=%i&layer=%s' %
                       (obj.ra, obj.dec, subwcs.pixel_scale() * scale, rw, rh, layer))
                print('URL:', url)
                r = requests.get(url)
                ftmp = tempfile.NamedTemporaryFile()
                ftmp.write(r.content)
                ftmp.flush()
                #fits,hdr = fitsio.read(ftmp.name, header=True)
                fitsfile = fitsio.FITS(ftmp.name)
                ftmp.close()
                print('FITS file:', len(fitsfile), 'extensions')
                hdr = fitsfile[0].read_header()
                fits = fitsfile[0].read()
                #print('fits:', fits)
                if fits is None:
                    print('no coverage in layer', layer)
                    continue

                # If you need to keep a copy (debugging...)
                # f,tmpfn = tempfile.mkstemp(suffix='.fits')
                # os.write(f, r.content)
                # os.close(f)
                # fits,hdr = fitsio.read(tmpfn, header=True)
                # print('Wrote FITS to', tmpfn)
                if np.all(fits == 0):
                    continue

                ### HACK -- surface brightness correction...
                if layer == 'sdssco':
                    s = (subwcs.pixel_scale() / 0.396)
                    fits *= s**2

                N,ww,hh = fits.shape
                # pull out the image planes
                imgs = [fits[n,:,:] for n in range(N)]
                bands = hdr['BANDS'].strip()
                fitsimgs.append((layer, bands, imgs))

                # Resample the new image to this layer's WCS
                if newimg is not None:
                    continue

                thiswcs = Tan(hdr)
                newimg = np.zeros((rh,rw), dtype=subimg.dtype)
                try:
                    #Yo,Xo,Yi,Xi,rims = resample_with_wcs(thiswcs, subwcs)
                    # Laczos
                    Yo,Xo,Yi,Xi,rims = resample_with_wcs(thiswcs, subwcs, [subimg])
                except:
                    continue
                #newimg[Yo,Xo] = subimg[Yi,Xi]
                newimg[Yo,Xo] = rims[0]

            if len(fitsimgs) == 0:
                return

            #print()
            newband = primhdr['FILTER'][0]
            #print('New image is', newband)

            plt.clf()
            plt.subplots_adjust(left=0.03, right=0.97, bottom=0.03)
            NC = 1 + len(fitsimgs)
            NR = 2

            zp = M['zp']
            zpscale = 10.**((zp - 22.5)/2.5)
            exptime = primhdr['EXPTIME']
            newimg /= (zpscale * exptime)

            nh,nw = newimg.shape
            coverage = np.sum(newimg != 0) / float(nw*nh)
            print('Fraction', coverage, 'of new image has data')

            def my_rgb(imgs, bands, **kwargs):
                return sdss_rgb(imgs, bands,
                                scales=dict(g=6.0, r=3.4, i=2.5, z=2.2),
                                m=-0.02, clip=False, **kwargs)

            def grayscale(img, band):
                rgb = my_rgb([img,img,img], [band,band,band])
                index = 'zrg'.index(newband)
                gray = rgb[:,:,index]
                return gray

            # plot title args
            targs = dict(fontsize=8)

            # New Image plot
            newgray = grayscale(newimg, newband)
            hi = np.percentile(newgray, 99.9)
            grayargs = dict(interpolation='nearest', origin='lower',
                            vmin=0., vmax=hi, cmap='gray')
            plt.subplot(NR, NC, 1)
            plt.imshow(newgray, **grayargs)
            plt.xticks([]); plt.yticks([])
            plt.title('New image (%s)' % newband, **targs)

            # Select images & bands for the New RGB image.
            zeroimg = np.zeros_like(newimg)
            newimgs = [zeroimg, zeroimg, zeroimg]
            # sdss_rgb reverses the order, so do like grz.
            newbands = ['g','r','z']
            newindex = dict(g=0, r=1, i=2, z=2)

            # Start with the new image plane of the New RGB image
            j = newindex[newband]
            newimgs [j] = newimg
            newbands[j] = newband
            #print('Setting band', newband, '(index %i)'%j, 'to new image')

            # We wait to choose the scaling until after the "New RGB" image
            # has been built, so we store the RGB images along the way...
            rgbs = []
            # In the link to the legacysurvey viewer, we select the best
            # imaging layer.  The DECaLS layer is listed first, so if it
            # has three bands, it wins.
            bestdata = None
            bestn = 0

            for i,(layer, bands, imgs) in enumerate(fitsimgs):
                nicelayer = nicelayernames.get(layer, layer)
                # For DECaLS missing bands: "--z"
                nicebands = ''

                goodbands = []
                goodimgs = []

                for band,img in zip(bands, imgs):
                    if np.all(img == 0.):
                        print('Band', band, 'of', nicelayer, 'is all zero')
                        nicebands += '-'
                        continue
                    goodbands.append(band)
                    goodimgs.append(img)
                    nicebands += band
                    #print('  ', nicelayer, band)
                    j = newindex[band]
                    if newimgs[j] is zeroimg:
                        #print('    Setting band', band, '(index %i)'%j, 'to', nicelayer)
                        newimgs[j] = img
                        newbands[j] = band
                    elif band == newbands[j]:
                        # patch empty regions if same band
                        #print('    Patching index %i from' % j, nicelayer, 'band', band)
                        Z = (newimgs[j] == 0)
                        newimgs[j][Z] = img[Z]

                    # z -> i
                    if newindex[band] == newindex[newband]:
                        # Old image grayscale
                        oldgray = grayscale(img, band)
                        plt.subplot(NR, NC, 2+i)
                        plt.imshow(oldgray, **grayargs)
                        plt.xticks([]); plt.yticks([])
                        plt.title('%s (%s)' % (nicelayer, band), **targs)

                if len(goodbands) == 1:
                    #rgb = grayscale(goodimgs[0], goodbands[0])
                    img,band = goodimgs[0], goodbands[0]
                    rgb = my_rgb([img,img,img], [band,band,band])
                else:
                    rgb = my_rgb(imgs, bands)

                if len(goodbands) > bestn:
                    bestn = len(goodbands)
                    bestdata = layer

                # print('bands for', nicelayer, ':', bands, ', actually', nicebands)
                rgbs.append((2+i+NC, rgb, '%s (%s)' % (nicelayer, nicebands)))

            # list to string
            newbands = ''.join(newbands)
            #print('Newbands:', newbands)

            # New RGB
            plt.subplot(NR, NC, 1 + NC)
            rgb = my_rgb(newimgs, newbands)
            lo = 0.
            hi = np.percentile(rgb.ravel(), 99.9)
            rgb = np.clip((rgb - lo) / (hi - lo), 0., 1.)
            plt.imshow(rgb, interpolation='nearest', origin='lower')
            plt.xticks([]); plt.yticks([])
            plt.title('New+Old (%s)' % newbands, **targs)

            # Old RGBs
            for sp, rgb, tt in rgbs:
                plt.subplot(NR, NC, sp)
                rgb = np.clip((rgb - lo) / (hi - lo), 0., 1.)
                plt.imshow(rgb, interpolation='nearest', origin='lower')
                plt.xticks([]); plt.yticks([])
                plt.title(tt, **targs)

            # NGC classification
            info = ''
            info = ngc_typenames.get(obj.classification.strip(), '')
            if len(info):
                info = '(' + info + ') '

            plt.suptitle('%s %sin %s %i-%s: %s band' %
                         (obj.name, info, nice_camera_name, expnum, extname, newband))

            # Save / show plot
            if self.opt.show:
                plt.draw()
                plt.show(block=False)
                plt.pause(0.001)
            plotfn = ('ngcbot-%i-%s-%s.png' %
                      (expnum, extname, obj.name.replace(' ','_')))
            if self.opt.plotdir:
                plotfn = os.path.join(self.opt.plotdir, plotfn)
            plt.savefig(plotfn)
            print('Saved', plotfn)
            #plt.savefig('ngcbot-latest.png')

            # Compose tweet text
            if self.opt.tweet:
                import urllib
                url = 'http://legacysurvey.org/viewer/?ra=%.3f&dec=%.3f' % (obj.ra, obj.dec)
                if bestdata is not None:
                    url += '&layer=%s' % bestdata
                url2 = ('http://ned.ipac.caltech.edu/cgi-bin/objsearch?' + 
                         urllib.urlencode(dict(objname=obj.name, corr_z=1, list_limit=5, img_stamp='YES')))

                dateobs = primhdr['DATE-OBS'].strip()
                # 2017-03-06T00:15:40.101482
                dateobs = dateobs[:19]
                dateobs = dateobs.replace('T', ' ')

                txt = ('%s observed %s %sin image %i-%s: %s band' %
                       (nice_camera_name, obj.name, info, expnum, extname, newband)
                       + ' at %s UT' % dateobs
                       + '\n' + url + '\n' + url2)
                print('Tweet text:', txt)

                if coverage > 0.75:
                    tweets.append((txt, plotfn))

        # Tweet one NGC object per exposure, chosen randomly.
        if len(tweets):
            i = np.random.randint(0, len(tweets))
            txt,plotfn = tweets[i]
            if self.opt.tweet:
                send_tweet(txt, plotfn)
            # Copy one plot to ngcbot-latest.png
            os.rename(plotfn, 'ngcbot-latest.png')

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
    #print(response)
    print('Tweet posted!')

if __name__ == '__main__':
    main()
