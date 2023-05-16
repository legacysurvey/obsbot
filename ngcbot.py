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
    #legacy_survey_layers = ['decals-dr3', 'sdssco']
    legacy_survey_layers = ['ls-dr10-grz', 'sdss']
else:
    legacy_survey_layers = ['mzls+bass-dr4', 'sdssco']

nicelayernames = { 'decals-dr3': 'DECaLS DR3',
                   'mzls+bass-dr4': 'MzLS+BASS DR4',
                   'sdssco': 'SDSS'}
nicelayernames.update({'ls-dr10': 'Legacy Surveys DR10',
                       'sdss': 'SDSS'})
# ngc_typenames = {
#     'Gx': 'galaxy',
#     'Gb': 'globular cluster',
#     'Nb': 'nebula',
#     'Pl': 'planetary nebula',
#     'C+N': 'cluster+nebulosity',
#     'Kt': 'knot in galaxy',
# }

ngc_typenames = {
    'G': 'galaxy',
    'GPair': 'galaxy pair',
    'GTrpl': 'galaxy triple',
    'GGroup': 'galaxy group',
    'Other': 'other',
    'GCl': 'globular cluster',
    'PN': 'planetary nebula',
    'Neb': 'nebula',
    'HII': 'HII region',
    'Cl+N': 'cluster+nebulosity',
    'RfN': 'reflection nebula',
    'SNR': 'SNR?',
    'EmN': 'emission nebula',
    'Nova': 'nova',
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
    parser.add_option('--coverage', type=float, help='Set coverage fraction threshold for sending tweets')
    parser.add_option('--only', default=False, action='store_true',
                      help='Only run the command-line files, then quit')
    parser.add_option('--flats', default=False, action='store_true',
                      help='Look for and use flats in the image directory?')
    parser.add_option('--ps', default=False, action='store_true',
                      help='Debugging plots?')
    
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
        #fn = os.path.join(catdir, 'ngc2.fits') #'openngc-ngc.fits') #'ngc2000.fits')
        fn = 'openngc-ngc-custom.fits'
        print('Reading', fn)
        T = fits_table(fn)
        T.delete_column('ngcnum')
        TT.append(T)

        #fn = os.path.join(catdir, 'ic2.fits') #'openngc-ic.fits') #'ic2000.fits')
        fn = 'openngc-ic-custom.fits'
        print('Reading', fn)
        T = fits_table(fn)
        T.delete_column('icnum')
        TT.append(T)
        self.cat = merge_tables(TT, columns='fillzero')
        del TT

        self.cat.rename('type', 'classification')

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

        #Counter({'G': 9742, 'OCl': 649, '*': 546, 'Other': 433, '**': 245, 'GPair': 211, 'GCl': 205, 'PN': 129, 'Neb': 93, 'HII': 81, 'Cl+N': 65, '*Ass': 60, 'RfN': 37, 'GTrpl': 23, 'NonEx': 14, 'GGroup': 14, 'SNR': 11, 'EmN': 8, 'Nova': 3})
        omit = [
            'OCl',
            '*',
            '**',
            '*Ass',
            'NonEx',
        ]

        print(Counter(self.cat.classification))

        self.cat.cut(np.flatnonzero(np.array([t.strip() not in omit for t in self.cat.classification])))
        print('Cut to', len(self.cat), 'NGC/IC objects')

        print('Remaining classifications:', np.unique(self.cat.classification))

        #self.spec = fits_table('specObj-dr12-trim-2.fits')

        self.cached_flats = {}
        self.read_flats = opt.flats
        
    def try_open_file(self, path):
        print('Trying to open file: %s' % path)
        F = fitsio.FITS(path)
        ## HACK
        if len(F) < min_n_exts:
            print('Got', len(F), 'extensions')
            raise RuntimeError('fewer extensions than expected')

    def filter_new_files(self, fns):
        # MzLS: ignore .dat files
        return [fn for fn in fns if
                fn.endswith('.fits.fz') or fn.endswith('.fits')]

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

        if self.match_ngc(rr, dd, radius, exts, F, primhdr, meas):
            return

        print('Not doing SDSS spectro jazz.')
        return

        # No match with NGC catalog -- look at SDSS spectro objects.
        I,J,d = match_radec(rr, dd, self.spec.ra, self.spec.dec, radius)
        print(len(I), 'matches to spectra')
        if len(I) == 0:
            return False

        #meas_exts = set()
        #K = []
        measures = {}

        specobjs = []
        
        for k,(i,j) in enumerate(zip(I,J)):
            ext = exts[i]
            obj = self.spec[j]
            obj.name = obj.label.strip()
            hdr = F[ext].read_header()
            extname = hdr['EXTNAME'].strip()
            expnum = primhdr['EXPNUM']
            wcs = meas.get_wcs(hdr)
            ok,x,y = wcs.radec2pixelxy(obj.ra, obj.dec)
            x = x - 1
            y = y - 1
            # Choose cutout area
            pixrad = 1.4 * 50
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

            #meas_exts.add(ext)
            #K.append(k)
            if ext in measures:
                M = measures[ext]
            else:
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
                measures[ext] = M

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

            specobjs.append((tt, subwcs, subimg))
            
            # # What size of image are we going to request?
            # scale = 1.
            # 
            # rh,rw = int(np.ceil(sh/scale)),int(np.ceil(sw/scale))
            # # make it square
            # mx = max(rh, rw)
            # rw = rh = mx
            # 
            # fitsimgs = []
            # # We'll resample the new image into the existing-image WCS.
            # newimg = None
            # 
            # # HACK
            # imgs = []
            # layer = 'sdssco'
            # bands = 'gri'
            # for band in bands:
            #     fn = 'sdssco-1679p492-%s.fits' % band
            #     fitsfile = fitsio.FITS(fn)
            #     wcs = Tan(fn, 0)
            #     hh,ww = wcs.shape
            #     ok,x,y = wcs.radec2pixelxy(obj.ra, obj.dec)
            #     x = x - 1
            #     y = y - 1
            #     print('x,y', x,y)
            #     sz = pixrad * 0.262/0.396
            #     xl,xh = int(np.clip(x-sz, 0, ww-1)), int(np.clip(x+sz, 0, ww-1))
            #     yl,yh = int(np.clip(y-sz, 0, hh-1)), int(np.clip(y+sz, 0, hh-1))
            #     if xl == xh or yl == yh:
            #         print('no overlap with SDSS image')
            #         continue
            #     
            #     img = fitsfile[1][yl:yh+1, xl:xh+1]
            #     s = (subwcs.pixel_scale() / 0.396)
            #     img *= s**2
            #     imgs.append(img)
            #     thiswcs = wcs
            #     
            # if len(imgs) == 0:
            #     continue
            # 
            # fitsimgs.append((layer, bands, imgs))
            # 
            # # Resample the new image to this layer's WCS
            # newimg = np.zeros((rh,rw), dtype=subimg.dtype)
            # try:
            #     # Laczos
            #     Yo,Xo,Yi,Xi,rims = resample_with_wcs(thiswcs, subwcs, [subimg])
            # except:
            #     continue
            # newimg[Yo,Xo] = rims[0]
            # 
            # if len(fitsimgs) == 0:
            #     # No overlap with existing surveys
            #     continue
            # 
            # newband = primhdr['FILTER'][0]
            
            
        # I = I[K]
        # J = J[K]
        # print('Cut to', len(I), 'matches in', len(meas_exts), 'extensions')
        # if len(I) == 0:
        #     return
        # measures = {}
        # for ext in meas_ext:
        #     measures[ext] = 
        print(len(specobjs), 'objects')

        def my_rgb(imgs, bands, **kwargs):
            return sdss_rgb(imgs, bands,
                            scales=dict(g=6.0, r=3.4, i=2.5, z=2.2),
                            m=-0.02, clip=False, **kwargs)

        def grayscale(img, band):
            rgb = my_rgb([img,img,img], [band,band,band])
            index = 'zrg'.index(newband)
            gray = rgb[:,:,index]
            return gray

        newband = primhdr['FILTER'][0]

        plt.figure(2)
        plt.clf()
        plt.subplots_adjust(left=0.03, right=0.97, bottom=0.03,
                            hspace=0, wspace=0)
        
        N = len(specobjs)
        C = int(np.ceil(np.sqrt(N * 1.2)))
        R = int(np.ceil(N / float(C)))
        print('Rows,cols, N', R,C, N)

        gray = grayscale(raw, newband)
        hi = np.percentile(gray, 99.9)
        grayargs = dict(interpolation='nearest', origin='lower',
                        vmin=0., vmax=hi, cmap='gray')

        print('Plotting:')
        plt.clf()
        for i,(tt, subwcs,subimg) in enumerate(specobjs):
            plt.subplot(R, C, i+1)

            print('  ', tt)
            
            newimg = subimg
            # New Image plot
            newgray = grayscale(newimg, newband)
            #hi = np.percentile(newgray, 99.9)
            plt.imshow(newgray, **grayargs)
            plt.xticks([]); plt.yticks([])
        ps.savefig()

    def get_flat(self, band, ext, meas):
        f = self.cached_flats.get((band,ext), None)
        if f is not None:
            return f
        fns = self.get_file_list()
        fns = self.filter_new_files(fns)
        fns.sort()
        print(len(fns), 'files to search for flats')
        flats = []
        for fn in fns:
            try:
                F = fitsio.FITS(fn)
                hdr = F[0].read_header()
                obstype = hdr['OBSTYPE']
                fband = meas.get_band(hdr)
                print('File', fn, 'obstype', obstype, 'filter', fband)
                if obstype != 'dome flat':
                    continue
                if fband != band:
                    continue
                flat,fhdr = meas.read_raw(F, ext)
                print('Got flat:', flat.shape)
                flats.append(flat)
            except:
                import traceback
                traceback.print_exc()
                pass
        if len(flats) == 0:
            return None
        flats = np.dstack(flats)
        print('Flats:', flats.shape)
        f = np.median(flats, axis=2)
        print('f:', f.shape)
        self.cached_flats[(band,ext)] = f
        return f

    def match_ngc(self, rr, dd, radius, exts, F, primhdr, meas):
        '''
        rr, dd: np arrays: RA,Dec centers of chips/amps
        radius: scalar, radius of chips/amps
        exts: list (same len as rr,dd) of extension indices
        F: fitsio.FITS object
        meas: measurement class
        '''
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
            return False
        for j in J:
            info = ngc_typenames.get(self.cat.classification[j].strip(), '')
            print('  ', self.cat.name[j], info)

        # Potential tweet texts and plot filenames
        tweets = []
        goodplots = []
        
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

            # HACK
            #pixrad *= 3

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
            flat = None
            if self.read_flats:
                band = meas.get_band(hdr)
                flat = self.get_flat(band, ext, meas)
                if flat is not None:
                    print('flat: range', flat.min(), flat.max(), 'median', np.median(flat))
            meas.ext = ext
            meas.edge_trim = 20
            debugps = None
            if self.opt.ps:
                debugps = ps
            try:
                M = meas.run(n_fwhm=1, verbose=False, get_image=True, flat=flat, ps=debugps)
            except KeyboardInterrupt:
                sys.exit(0)
            except:
                import traceback
                print('Failed to measure file', meas.fn, 'ext', ext, ':')
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

            if False:
                # Flat image for the same CCD region.
                flatfn = '/tmp/mzls/mos3.127506.fits'
                # meas_class = get_measurer_class_for_file(flatfn)
                # if meas_class is None:
                #     print('Failed to identify camera in', flatfn)
                #     return
                # flatmeas = meas_class(flatfn, 0, self.nom)
                # print('Flat meas:', flatmeas)
                # flatmeas.ext = ext
                # flathdr = fitsio.read_header(flatfn)
                # flatmeas.primhdr = flathdr
                # flatmeas.edge_trim = 20
                # FM = flatmeas.run(n_fwhm=1, verbose=False, get_image=True)
                # flatraw = FM['image']
                F = fitsio.FITS(flatfn)
                flatraw,flathdr = meas.read_raw(F, ext)
                flatsub = flatraw[yl:yh, xl:xh]
    
                zerofn = '/tmp/mzls/mos3.127522.fits'
                F = fitsio.FITS(zerofn)
                zeroraw,zerohdr = meas.read_raw(F, ext)
                zerosub = zeroraw[yl:yh, xl:xh]
    
                rawfn = meas.fn
                F = fitsio.FITS(rawfn)
                rawimg,rawhdr = meas.read_raw(F, ext)
                rawsub = rawimg[yl:yh, xl:xh]
            
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
                try:
                    fitsfile = fitsio.FITS(ftmp.name)
                except:
                    print('no coverage in layer', layer)
                    continue
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
                if layer in ['sdssco', 'sdss']:
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
                # No overlap with existing surveys
                continue


            #print()
            newband = primhdr['FILTER'][0]
            #print('New image is', newband)

            if False:
                mn,mx = np.percentile(flatsub, [25,99])
                plt.clf()
                #plt.imshow(newflat, interpolation='nearest', origin='lower',
                plt.imshow(flatsub, interpolation='nearest', origin='lower',
                           vmin=mn, vmax=mx)
                #plt.colorbar()
                plt.colorbar()
                plt.savefig('ngcbot-flat.png')
    
                med = np.median(newflat.ravel())
                #mn,mx = np.percentile(newimg, [25,99])
                #mn,mx = np.percentile(newimg, [25,90])
                #mn,mx = np.percentile(rawsub, [25,95])
                mn,mx = np.percentile(rawsub, [50,95])
                plt.clf()
                # plt.subplot(1,2,1)
                #plt.imshow(newimg, interpolation='nearest', origin='lower',
                # plt.imshow(subimg, interpolation='nearest', origin='lower',
                plt.imshow(rawsub, interpolation='nearest', origin='lower',
                           vmin=mn, vmax=mx)
                plt.colorbar()
                plt.savefig('ngcbot-unflattened.png')
                plt.clf()
                #plt.subplot(1,2,2)
                #plt.imshow(newimg / (newflat / med), interpolation='nearest', origin='lower',
                #plt.imshow(subimg / (flatsub / med), interpolation='nearest', origin='lower',
                plt.imshow(rawsub / (flatsub / med), interpolation='nearest', origin='lower',
                           vmin=mn, vmax=mx)
                plt.colorbar()
                plt.savefig('ngcbot-flattened.png')
    
                mn,mx = np.percentile(zerosub, [5,95])
                plt.clf()
                plt.imshow(zerosub, interpolation='nearest', origin='lower',
                           vmin=mn, vmax=mx)
                plt.colorbar()
                plt.savefig('ngcbot-zero.png')

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


            # DEBUG
            if self.opt.ps:
                ocx,ocy = subwcs.get_crpix()

                subwcs.set_crpix((cx, cy))
                newimgA = np.zeros((rh,rw), dtype=subimg.dtype)
                Yo,Xo,Yi,Xi,rims = resample_with_wcs(thiswcs, subwcs, [subimg])
                newimgA[Yo,Xo] = rims[0]

                subwcs.set_crpix((cx-dx, cy-dy))
                newimgB = np.zeros((rh,rw), dtype=subimg.dtype)
                Yo,Xo,Yi,Xi,rims = resample_with_wcs(thiswcs, subwcs, [subimg])
                newimgB[Yo,Xo] = rims[0]

                subwcs.set_crpix((cx-dx-corrx, cy-dy-corry))
                newimgC = np.zeros((rh,rw), dtype=subimg.dtype)
                Yo,Xo,Yi,Xi,rims = resample_with_wcs(thiswcs, subwcs, [subimg])
                newimgC[Yo,Xo] = rims[0]

                #newgray = grayscale(newimg, newband)
                #hi = np.percentile(newgray, 99.9)
                grayargs = dict(interpolation='nearest', origin='lower',
                                cmap='gray', vmin=0.)
                #vmin=0., vmax=hi, cmap='gray')

                (layer, bands, imgs) = fitsimgs[0]
                nicelayer = nicelayernames.get(layer, layer)
                for band,img in zip(bands, imgs):
                    newbands = ['g','r','z']
                    newindex = dict(g=0, r=1, i=2, z=2)
                    if newindex[band] == newindex[newband]:
                        oldgray = grayscale(img, band)
                        plt.clf()
                        plt.imshow(oldgray, **grayargs)
                        plt.title(nicelayer)
                        ps.savefig()
                
                plt.clf()
                plt.imshow(grayscale(newimgA, newband), **grayargs)
                #plt.imshow(newimgA, **grayargs)
                plt.title('No astromety correction')
                ps.savefig()

                plt.clf()
                plt.imshow(grayscale(newimgB, newband), **grayargs)
                #plt.imshow(newimgB, **grayargs)
                plt.title('Astromety shift only correction')
                ps.savefig()

                plt.clf()
                plt.imshow(grayscale(newimgC, newband), **grayargs)
                #plt.imshow(newimgC, **grayargs)
                plt.title('Astromety shift+affine correction')
                ps.savefig()

                subwcs.set_crpix((ocx,ocy))
            
            plt.clf()
            plt.subplots_adjust(left=0.03, right=0.97, bottom=0.03)
            NC = 1 + len(fitsimgs)
            NR = 2
            
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

            good_coverage = 0.75
            if self.opt.coverage is not None:
                good_coverage = self.opt.coverage
            if coverage > good_coverage:
                goodplots.append(plotfn)
            
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

                if coverage > good_coverage:
                    tweets.append((txt, plotfn))
                else:
                    print('Coverage', coverage, 'less than target', good_coverage, 'so not tweeting')

        # Select a random good-looking plot
        if len(goodplots):
            irandom = np.random.randint(0, len(goodplots))
            # Copy that plot to ngcbot-latest.png
            plotfn = goodplots[irandom]
            import shutil
            latest = 'ngcbot-latest.png'
            print('Copying', plotfn, 'to', latest)
            shutil.copy(plotfn, latest)
        # Tweet one NGC object per exposure, chosen randomly.
        if len(tweets):
            assert(len(tweets) == len(goodplots))
            txt,plotfn = tweets[irandom]
            if self.opt.tweet:
                send_tweet(txt, plotfn)
        return True

def send_tweet(txt, imgfn):
    from twython import Twython
    try:
        import abotsecrets as s
    except:
        print('You need a "botsecrets.py" file with Twitter credentials')
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
