from __future__ import print_function
import os
import tempfile
import json

if __name__ == '__main__':
    import matplotlib
    matplotlib.use('Agg')
import numpy as np

import fitsio

from legacypipe.ps1cat import ps1cat, ps1_to_decam, ps1_to_90prime
from legacypipe.gaiacat import GaiaCatalog

# Color terms -- no MOSAIC specific ones yet:
ps1_to_mosaic = ps1_to_decam

class RawMeasurer(object):
    def __init__(self, fn, ext, nom, aprad=3.5, skyrad_inner=7.,
                 skyrad_outer=10.,
                 minstar=5, pixscale=0.262, maxshift=240.):
        '''
        aprad: float
        Aperture photometry radius in arcsec

        skyrad_{inner,outer}: floats
        Sky annulus radius in arcsec

        nom: nominal calibration, eg NominalCalibration object
        '''

        self.fn = fn
        self.ext = ext
        self.nom = nom
        
        self.aprad = aprad
        self.skyrad = (skyrad_inner, skyrad_outer)
        self.minstar = minstar
        self.pixscale = pixscale
        self.maxshift = maxshift

        self.edge_trim = 100
        self.nominal_fwhm = 5.0
        # Detection threshold
        self.det_thresh = 20

        self.debug = True

        self.camera = 'camera'
        
    def remove_sky_gradients(self, img):
        from scipy.ndimage.filters import median_filter
        # Ugly removal of sky gradients by subtracting median in first x and then y
        H,W = img.shape
        meds = np.array([np.median(img[:,i]) for i in range(W)])
        meds = median_filter(meds, size=5)
        img -= meds[np.newaxis,:]
        meds = np.array([np.median(img[i,:]) for i in range(H)])
        meds = median_filter(meds, size=5)
        img -= meds[:,np.newaxis]

    def trim_edges(self, img):
        # Trim off some extra pixels -- image edges are often bright...
        # central part of the image
        trim = self.edge_trim
        if trim == 0:
            return img,trim,trim
        cimg = img[trim:-trim, trim:-trim]
        return cimg, trim, trim

    def read_raw(self, F, ext):
        img = F[ext].read()
        hdr = F[ext].read_header()
        return img,hdr

    def get_band(self, primhdr):
        band = primhdr['FILTER']
        band = band.split()[0]
        # HACK PTF: R -> r
        #band = band.lower()
        return band

    def match_ps1_stars(self, px, py, fullx, fully, radius, stars):
        from astrometry.libkd.spherematch import match_xy
        #print('Matching', len(px), 'PS1 and', len(fullx), 'detected stars with radius', radius)
        I,J,d = match_xy(px, py, fullx, fully, radius)
        #print(len(I), 'matches')
        dx = px[I] - fullx[J]
        dy = py[I] - fully[J]
        return I,J,dx,dy

    def detection_map(self, img, sig1, psfsig, ps):
        from scipy.ndimage.filters import gaussian_filter
        psfnorm = 1./(2. * np.sqrt(np.pi) * psfsig)
        detsn = gaussian_filter(img / sig1, psfsig) / psfnorm
        # zero out the edges -- larger margin here?
        detsn[0 ,:] = 0
        detsn[:, 0] = 0
        detsn[-1,:] = 0
        detsn[:,-1] = 0
        return detsn

    def detect_sources(self, detsn, thresh, ps):
        from scipy.ndimage.measurements import label, find_objects
        # HACK -- Just keep the brightest pixel in each blob!
        peaks = (detsn > thresh)
        blobs,nblobs = label(peaks)
        slices = find_objects(blobs)
        return slices

    def zeropoint_for_exposure(self, band, ext=None, exptime=None, primhdr=None):
        try:
            zp0 = self.nom.zeropoint(band, ext=self.ext)
        except KeyError:
            print('Unknown band "%s"; no zeropoint available.' % band)
            zp0 = None
        return zp0
        
    def get_reference_stars(self, wcs, band):
        # Read in the PS1 catalog, and keep those within 0.25 deg of CCD center
        # and those with main sequence colors
        pscat = ps1cat(ccdwcs=wcs)
        try:
            stars = pscat.get_stars()
        except:
            print('Failed to find PS1 stars -- maybe this image is outside the PS1 footprint.')
            import traceback
            traceback.print_exc()

            from astrometry.util.starutil_numpy import radectolb
            rc,dc = wcs.radec_center()
            print('RA,Dec center:', rc, dc)
            lc,bc = radectolb(rc, dc)
            print('Galactic l,b:', lc[0], bc[0])
            return None
        #print('Got PS1 stars:', len(stars))

        # we add the color term later
        try:
            ps1band = self.get_ps1_band(band)
            stars.mag = stars.median[:, ps1band]
        except KeyError:
            print('Unknown band for PS1 mags:', band)
            ps1band = None
            stars.mag = np.zeros(len(stars), np.float32)
        return stars

    def cut_reference_catalog(self, stars):
        # Now cut to just *stars* with good colors
        stars.gicolor = stars.median[:,0] - stars.median[:,2]
        keep = (stars.gicolor > 0.4) * (stars.gicolor < 2.7)
        return keep

    def get_color_term(self, stars, band):
        #print('PS1 stars mags:', stars.median)
        try:
            colorterm = self.colorterm_ps1_to_observed(stars.median, band)
        except KeyError:
            print('Color term not found for band "%s"; assuming zero.' % band)
            colorterm = 0.
        return colorterm

    def get_exptime(self, primhdr):
        return primhdr['EXPTIME']

    def run(self, ps=None, primext=0, focus=False, momentsize=5,
            n_fwhm=100, verbose=True, get_image=False, flat=None):
        import pylab as plt
        from astrometry.util.plotutils import dimshow, plothist
        import photutils
        import tractor

        if verbose:
            printmsg = print
        else:
            def printmsg(*args, **kwargs):
                pass

        fn = self.fn
        ext = self.ext
        pixsc = self.pixscale

        F = fitsio.FITS(fn)
        primhdr = F[primext].read_header()
        self.primhdr = primhdr
        img,hdr = self.read_raw(F, ext)
        self.hdr = hdr

        # pre sky-sub
        mn,mx = np.percentile(img.ravel(), [25,99])
        self.imgkwa = dict(vmin=mn, vmax=mx, cmap='gray')
        
        if self.debug and ps is not None:
            plt.clf()
            dimshow(img, **self.imgkwa)
            plt.title('Raw image')
            ps.savefig()
    
            M = 200
            plt.clf()
            plt.subplot(2,2,1)
            dimshow(img[-M:, :M], ticks=False, **self.imgkwa)
            plt.subplot(2,2,2)
            dimshow(img[-M:, -M:], ticks=False, **self.imgkwa)
            plt.subplot(2,2,3)
            dimshow(img[:M, :M], ticks=False, **self.imgkwa)
            plt.subplot(2,2,4)
            dimshow(img[:M, -M:], ticks=False, **self.imgkwa)
            plt.suptitle('Raw image corners')
            ps.savefig()


        if flat is not None:
            I = (flat != 0)
            med = np.median(flat[I]).astype(float)
            img[I] = img[I] / (flat[I] / med)
            #img /= (flat / np.median(flat).astype(float))

        img,trim_x0,trim_y0 = self.trim_edges(img)

        fullH,fullW = img.shape

        if self.debug and ps is not None:
            plt.clf()
            dimshow(img, **self.imgkwa)
            plt.title('Trimmed image')
            ps.savefig()
    
            M = 200
            plt.clf()
            plt.subplot(2,2,1)
            dimshow(img[-M:, :M], ticks=False, **self.imgkwa)
            plt.subplot(2,2,2)
            dimshow(img[-M:, -M:], ticks=False, **self.imgkwa)
            plt.subplot(2,2,3)
            dimshow(img[:M, :M], ticks=False, **self.imgkwa)
            plt.subplot(2,2,4)
            dimshow(img[:M, -M:], ticks=False, **self.imgkwa)
            plt.suptitle('Trimmed corners')
            ps.savefig()
            
        band = self.get_band(primhdr)
        exptime = self.get_exptime(primhdr)
        airmass = primhdr['AIRMASS']
        printmsg('Band', band, 'Exptime', exptime, 'Airmass', airmass)
        # airmass can be 'NaN'
        airmass = float(airmass)
        if not np.isfinite(airmass):
            printmsg('Bad airmass:', airmass, '-- setting to 1.0')
            airmass = 1.0

        zp0 = self.zeropoint_for_exposure(band, ext=self.ext, exptime=exptime, primhdr=primhdr)
        printmsg('Nominal zeropoint:', zp0)

        try:
            sky0 = self.nom.sky(band)
            printmsg('Nominal sky:', sky0)
        except KeyError:
            print('Unknown band "%s"; no nominal sky available.' % band)
            sky0 = None

        #print('Nominal calibration:', self.nom)
        #print('Fid:', self.nom.fiducial_exptime(band))
        try:
            kx = self.nom.fiducial_exptime(band).k_co
        except:
            print('Unknown band "%s"; no k_co available.' % band)
            kx = None

        # Find the sky value and noise level
        sky,sig1 = self.get_sky_and_sigma(img)
        sky1 = np.median(sky)

        if zp0 is not None:
            skybr = -2.5 * np.log10(sky1/pixsc/pixsc/exptime) + zp0
            printmsg('Sky brightness: %8.3f mag/arcsec^2' % skybr)
            printmsg('Fiducial:       %8.3f mag/arcsec^2' % sky0)
        else:
            skybr = None

        # Read WCS header and compute boresight
        wcs = self.get_wcs(hdr)
        ra_ccd,dec_ccd = wcs.pixelxy2radec((fullW+1)/2., (fullH+1)/2.)

        camera = primhdr.get('INSTRUME','').strip().lower()
        # -> "decam" / "mosaic3"
        meas = dict(band=band, airmass=airmass, exptime=exptime,
                    skybright=skybr, rawsky=sky1,
                    pixscale=pixsc, primhdr=primhdr,
                    hdr=hdr, wcs=wcs, ra_ccd=ra_ccd, dec_ccd=dec_ccd,
                    extension=ext, camera=camera)

        if get_image:
            meas.update(image=img, trim_x0=trim_x0, trim_y0=trim_y0)
        
        if skybr is not None and not np.isfinite(skybr):
            print('Measured negative sky brightness:', sky1, 'counts')
            return meas
        
        img -= sky

        self.remove_sky_gradients(img)

        # Post sky-sub
        mn,mx = np.percentile(img.ravel(), [25,98])
        self.imgkwa = dict(vmin=mn, vmax=mx, cmap='gray')
        
        if ps is not None:
            plt.clf()
            dimshow(img, **self.imgkwa)
            plt.title('Sky-sub image: %s-%s' % (os.path.basename(fn).replace('.fits','').replace('.fz',''), ext))
            plt.colorbar()
            ps.savefig()

        # Detect stars
        psfsig = self.nominal_fwhm / 2.35
        detsn = self.detection_map(img, sig1, psfsig, ps)
    
        slices = self.detect_sources(detsn, self.det_thresh, ps)
        printmsg(len(slices), 'sources detected')
        if len(slices) < 20:
            slices = self.detect_sources(detsn, 10., ps)
            printmsg(len(slices), 'sources detected')
        ndetected = len(slices)
        meas.update(ndetected=ndetected)
        if ndetected == 0:
            print('NO SOURCES DETECTED')
            return meas

        # Measure moments, etc around detected stars.
        # These measurements don't seem to be very good!
        xx,yy = [],[]
        fx,fy = [],[]
        mx2,my2,mxy = [],[],[]
        wmx2,wmy2,wmxy = [],[],[]
        # "Peak" region to centroid
        P = momentsize
        H,W = img.shape

        for i,slc in enumerate(slices):
            y0 = slc[0].start
            x0 = slc[1].start
            subimg = detsn[slc]
            imax = np.argmax(subimg)
            y,x = np.unravel_index(imax, subimg.shape)
            if (x0+x) < P or (x0+x) > W-1-P or (y0+y) < P or (y0+y) > H-1-P:
                #print('Skipping edge peak', x0+x, y0+y)
                continue
            xx.append(x0 + x)
            yy.append(y0 + y)
            pkarea = detsn[y0+y-P: y0+y+P+1, x0+x-P: x0+x+P+1]

            from scipy.ndimage.measurements import center_of_mass
            cy,cx = center_of_mass(pkarea)
            #print('Center of mass', cx,cy)
            fx.append(x0+x-P+cx)
            fy.append(y0+y-P+cy)
            #print('x,y', x0+x, y0+y, 'vs centroid', x0+x-P+cx, y0+y-P+cy)

            ### HACK -- measure source ellipticity
            # go back to the image (not detection map)
            #subimg = img[slc]
            subimg = img[y0+y-P: y0+y+P+1, x0+x-P: x0+x+P+1].copy()
            subimg /= subimg.sum()
            ph,pw = subimg.shape
            px,py = np.meshgrid(np.arange(pw), np.arange(ph))
            mx2.append(np.sum(subimg * (px - cx)**2))
            my2.append(np.sum(subimg * (py - cy)**2))
            mxy.append(np.sum(subimg * (px - cx)*(py - cy)))
            # Gaussian windowed version
            s = 1.
            wimg = subimg * np.exp(-0.5 * ((px - cx)**2 + (py - cy)**2) / s**2)
            wimg /= np.sum(wimg)
            wmx2.append(np.sum(wimg * (px - cx)**2))
            wmy2.append(np.sum(wimg * (py - cy)**2))
            wmxy.append(np.sum(wimg * (px - cx)*(py - cy)))

        mx2 = np.array(mx2)
        my2 = np.array(my2)
        mxy = np.array(mxy)
        wmx2 = np.array(wmx2)
        wmy2 = np.array(wmy2)
        wmxy = np.array(wmxy)

        # semi-major/minor axes and position angle
        theta = np.rad2deg(np.arctan2(2 * mxy, mx2 - my2) / 2.)
        theta = np.abs(theta) * np.sign(mxy)
        s = np.sqrt(((mx2 - my2)/2.)**2 + mxy**2)
        a = np.sqrt((mx2 + my2) / 2. + s)
        b = np.sqrt((mx2 + my2) / 2. - s)
        ell = 1. - b/a

        wtheta = np.rad2deg(np.arctan2(2 * wmxy, wmx2 - wmy2) / 2.)
        wtheta = np.abs(wtheta) * np.sign(wmxy)
        ws = np.sqrt(((wmx2 - wmy2)/2.)**2 + wmxy**2)
        wa = np.sqrt((wmx2 + wmy2) / 2. + ws)
        wb = np.sqrt((wmx2 + wmy2) / 2. - ws)
        well = 1. - wb/wa
            
        fx = np.array(fx)
        fy = np.array(fy)
        xx = np.array(xx)
        yy = np.array(yy)
            
        if ps is not None:
    
            plt.clf()
            dimshow(detsn, vmin=-3, vmax=50, cmap='gray')
            ax = plt.axis()
            plt.plot(fx, fy, 'go', mec='g', mfc='none', ms=10)
            plt.colorbar()
            plt.title('Detected sources')
            plt.axis(ax)
            ps.savefig()

        # Aperture photometry
        apxy = np.vstack((fx, fy)).T
        ap = []
        aprad_pix = self.aprad / pixsc
        aper = photutils.CircularAperture(apxy, aprad_pix)
        p = photutils.aperture_photometry(img, aper)
        apflux = p.field('aperture_sum')
    
        # Manual aperture photometry to get clipped means in sky annulus
        sky_inner_r, sky_outer_r = [r / pixsc for r in self.skyrad]
        sky = []
        for xi,yi in zip(fx,fy):
            ix = int(np.round(xi))
            iy = int(np.round(yi))
            skyR = int(np.ceil(sky_outer_r))
            xlo = max(0, ix-skyR)
            xhi = min(W, ix+skyR+1)
            ylo = max(0, iy-skyR)
            yhi = min(H, iy+skyR+1)
            xx,yy = np.meshgrid(np.arange(xlo,xhi), np.arange(ylo,yhi))
            r2 = (xx - xi)**2 + (yy - yi)**2
            inannulus = ((r2 >= sky_inner_r**2) * (r2 < sky_outer_r**2))
            skypix = img[ylo:yhi, xlo:xhi][inannulus]
            #print('ylo,yhi, xlo,xhi', ylo,yhi, xlo,xhi, 'img subshape', img[ylo:yhi, xlo:xhi].shape, 'inann shape', inannulus.shape)
            s,nil = sensible_sigmaclip(skypix)
            sky.append(s)
        sky = np.array(sky)

        apflux2 = apflux - sky * (np.pi * aprad_pix**2)
        good = (apflux2>0)*(apflux>0)
        apflux = apflux[good]
        apflux2 = apflux2[good]
        fx = fx[good]
        fy = fy[good]


        stars = self.get_reference_stars(wcs, band)
        if stars is None:
            return meas

        wcs2,measargs = self.update_astrometry(stars, wcs, fx, fy, trim_x0, trim_y0,
                                               pixsc, img, hdr, fullW, fullH, ps, printmsg)
        meas.update(measargs)

        keep = self.cut_reference_catalog(stars)
        stars.cut(keep)
        if len(stars) == 0:
            print('No overlap or too few stars in PS1')
            return meas
        ok,px,py = wcs2.radec2pixelxy(stars.ra, stars.dec)
        px -= 1
        py -= 1

        # Re-match with smaller search radius
        radius2 = 3. / pixsc
        I,J,dx,dy = self.match_ps1_stars(px, py, fx, fy,
                                         radius2, stars)
        printmsg('Keeping %i reference stars; matched %i' %
                 (len(stars), len(I)))
        nmatched = len(I)
        meas.update(nmatched=nmatched)

        if focus or zp0 is None:
            meas.update(img=img, hdr=hdr, primhdr=primhdr,
                        fx=fx, fy=fy,
                        sig1=sig1, stars=stars,
                        moments=(mx2,my2,mxy,theta,a,b,ell),
                        wmoments=(wmx2,wmy2,wmxy,wtheta,wa,wb,well),
                        apflux=apflux, apflux2=apflux2)
            return meas

        #print('Mean astrometric shift (arcsec): delta-ra=', -np.mean(dy)*0.263, 'delta-dec=', np.mean(dx)*0.263)
            
        # Compute photometric offset compared to PS1
        # as the PS1 minus observed mags

        colorterm = self.get_color_term(stars, band)
        #print('Color term:', colorterm)
        stars.mag += colorterm
        ps1mag = stars.mag[I]
        #print('PS1 mag:', ps1mag)

        if False and ps is not None:
            plt.clf()
            plt.semilogy(ps1mag, apflux2[J], 'b.')
            plt.xlabel('PS1 mag')
            plt.ylabel('DECam ap flux (with sky sub)')
            ps.savefig()
        
            plt.clf()
            plt.semilogy(ps1mag, apflux[J], 'b.')
            plt.xlabel('PS1 mag')
            plt.ylabel('DECam ap flux (no sky sub)')
            ps.savefig()

        # "apmag2" is the annular sky-subtracted flux
        apmag2 = -2.5 * np.log10(apflux2) + zp0 + 2.5 * np.log10(exptime)
        apmag  = -2.5 * np.log10(apflux ) + zp0 + 2.5 * np.log10(exptime)
    
        if ps is not None:
            plt.clf()
            plt.plot(ps1mag, apmag[J], 'b.', label='No sky sub')
            plt.plot(ps1mag, apmag2[J], 'r.', label='Sky sub')
            # ax = plt.axis()
            # mn = min(ax[0], ax[2])
            # mx = max(ax[1], ax[3])
            # plt.plot([mn,mx], [mn,mx], 'k-', alpha=0.1)
            # plt.axis(ax)
            plt.xlabel('PS1 mag')
            plt.ylabel('DECam ap mag')
            plt.legend(loc='upper left')
            plt.title('Zeropoint')
            ps.savefig()
    
        dm = ps1mag - apmag[J]
        dmag,dsig = sensible_sigmaclip(dm, nsigma=2.5)
        printmsg('Mag offset: %8.3f' % dmag)
        printmsg('Scatter:    %8.3f' % dsig)

        if not np.isfinite(dmag) or not np.isfinite(dsig):
            print('FAILED TO GET ZEROPOINT!')
            meas.update(zp=None)
            return meas

        from scipy.stats import sigmaclip
        goodpix,lo,hi = sigmaclip(dm, low=3, high=3)
        dmagmed = np.median(goodpix)
        #print(len(goodpix), 'stars used for zeropoint median')
        #print('Using median zeropoint:')
        zp_med = zp0 + dmagmed
        trans_med = 10.**(-0.4 * (zp0 - zp_med - kx * (airmass - 1.)))
        #print('Zeropoint %6.3f' % zp_med)
        #print('Transparency: %.3f' % trans_med)
        
        dm = ps1mag - apmag2[J]
        dmag2,dsig2 = sensible_sigmaclip(dm, nsigma=2.5)
        #print('Sky-sub mag offset', dmag2)
        #print('Scatter', dsig2)

        # Median, sky-sub
        goodpix,lo,hi = sigmaclip(dm, low=3, high=3)
        dmagmedsky = np.median(goodpix)

        if ps is not None:
            plt.clf()
            plt.plot(ps1mag, apmag[J] + dmag - ps1mag, 'b.', label='No sky sub')
            plt.plot(ps1mag, apmag2[J]+ dmag2- ps1mag, 'r.', label='Sky sub')
            plt.xlabel('PS1 mag')
            plt.ylabel('DECam ap mag - PS1 mag')
            plt.legend(loc='upper left')
            plt.ylim(-0.25, 0.25)
            plt.axhline(0, color='k', alpha=0.25)
            plt.title('Zeropoint')
            ps.savefig()

        zp_mean = zp0 + dmag

        zp_obs = zp0 + dmagmedsky
        transparency = 10.**(-0.4 * (zp0 - zp_obs - kx * (airmass - 1.)))
        meas.update(zp=zp_obs, transparency=transparency)

        printmsg('Zeropoint %6.3f' % zp_obs)
        printmsg('Fiducial  %6.3f' % zp0)
        printmsg('Transparency: %.3f' % transparency)

        # Also return the zeropoints using the sky-subtracted mags,
        # and using median rather than clipped mean.
        meas.update(zp_mean = zp_mean,
                    zp_skysub = zp0 + dmag2,
                    zp_med = zp_med,
                    zp_med_skysub = zp0 + dmagmedsky)

        # print('Using sky-subtracted values:')
        # zp_sky = zp0 + dmag2
        # trans_sky = 10.**(-0.4 * (zp0 - zp_sky - kx * (airmass - 1.)))
        # print('Zeropoint %6.3f' % zp_sky)
        # print('Transparency: %.3f' % trans_sky)
        
        fwhms = []
        psf_r = 15
        if n_fwhm not in [0, None]:
            Jf = J[:n_fwhm]
            
        for i,(xi,yi,fluxi) in enumerate(zip(fx[Jf],fy[Jf],apflux[Jf])):
            #print('Fitting source', i, 'of', len(Jf))
            ix = int(np.round(xi))
            iy = int(np.round(yi))
            xlo = max(0, ix-psf_r)
            xhi = min(W, ix+psf_r+1)
            ylo = max(0, iy-psf_r)
            yhi = min(H, iy+psf_r+1)
            xx,yy = np.meshgrid(np.arange(xlo,xhi), np.arange(ylo,yhi))
            r2 = (xx - xi)**2 + (yy - yi)**2
            keep = (r2 < psf_r**2)
            pix = img[ylo:yhi, xlo:xhi].copy()
            ie = np.zeros_like(pix)
            ie[keep] = 1. / sig1
            #print('fitting source at', ix,iy)
            #print('number of active pixels:', np.sum(ie > 0), 'shape', ie.shape)

            # plt.clf()
            # plt.imshow(pix, interpolation='nearest', origin='lower',
            #            vmin=-2.*sig1, vmax=5.*sig1)
            # plt.savefig('star-%03i.png' % i)

            psf = tractor.NCircularGaussianPSF([4.], [1.])
            psf.radius = psf_r
            tim = tractor.Image(data=pix, inverr=ie, psf=psf)
            src = tractor.PointSource(tractor.PixPos(xi-xlo, yi-ylo),
                                      tractor.Flux(fluxi))
            tr = tractor.Tractor([tim],[src])
    
            #print('Posterior before prior:', tr.getLogProb())
            src.pos.addGaussianPrior('x', 0., 1.)
            #print('Posterior after prior:', tr.getLogProb())
            
            doplot = (i < 5) * (ps is not None)
            if doplot:
                mod0 = tr.getModelImage(0)
    
            tim.freezeAllBut('psf')
            psf.freezeAllBut('sigmas')
    
            # print('Optimizing params:')
            # tr.printThawedParams()
    
            #print('Parameter step sizes:', tr.getStepSizes())
            optargs = dict(priors=False, shared_params=False)
            for step in range(50):
                dlnp,x,alpha = tr.optimize(**optargs)
                #print('dlnp', dlnp)
                #print('src', src)
                #print('psf', psf)
                if dlnp == 0:
                    break
            # Now fit only the PSF size
            tr.freezeParam('catalog')
            # print('Optimizing params:')
            # tr.printThawedParams()
    
            for step in range(50):
                dlnp,x,alpha = tr.optimize(**optargs)
                # print('dlnp', dlnp)
                # print('src', src)
                # print('psf', psf)
                if dlnp == 0:
                    break

            fwhms.append(psf.sigmas[0] * 2.35 * pixsc)
                
            if doplot:
                mod1 = tr.getModelImage(0)
                chi1 = tr.getChiImage(0)
            
                plt.clf()
                plt.subplot(2,2,1)
                plt.title('Image')
                dimshow(pix, **self.imgkwa)
                plt.subplot(2,2,2)
                plt.title('Initial model')
                dimshow(mod0, **self.imgkwa)
                plt.subplot(2,2,3)
                plt.title('Final model')
                dimshow(mod1, **self.imgkwa)
                plt.subplot(2,2,4)
                plt.title('Final chi')
                dimshow(chi1, vmin=-10, vmax=10)
                plt.suptitle('PSF fit')
                ps.savefig()
    
        fwhms = np.array(fwhms)
        fwhm = np.median(fwhms)
        printmsg('Median FWHM : %.3f arcsec' % np.median(fwhms))
        meas.update(seeing=fwhm)
    
        if False and ps is not None:
            lo,hi = np.percentile(fwhms, [5,95])
            lo -= 0.1
            hi += 0.1
            plt.clf()
            plt.hist(fwhms, 25, range=(lo,hi), histtype='step', color='b')
            plt.xlabel('FWHM (arcsec)')
            ps.savefig()
    
        if ps is not None:
            plt.clf()
            for i,(xi,yi) in enumerate(list(zip(fx[J],fy[J]))[:50]):
                ix = int(np.round(xi))
                iy = int(np.round(yi))
                xlo = max(0, ix-psf_r)
                xhi = min(W, ix+psf_r+1)
                ylo = max(0, iy-psf_r)
                yhi = min(H, iy+psf_r+1)
                pix = img[ylo:yhi, xlo:xhi]
        
                slc = pix[iy-ylo, :].copy()
                slc /= np.sum(slc)
                p1 = plt.plot(slc, 'b-', alpha=0.2)
                slc = pix[:, ix-xlo].copy()
                slc /= np.sum(slc)
                p2 = plt.plot(slc, 'r-', alpha=0.2)
                ph,pw = pix.shape
                cx,cy = pw/2, ph/2
                if i == 0:
                    xx = np.linspace(0, pw, 300)
                    dx = xx[1]-xx[0]
                    sig = fwhm / pixsc / 2.35
                    yy = np.exp(-0.5 * (xx-cx)**2 / sig**2) # * np.sum(pix)
                    yy /= (np.sum(yy) * dx)
                    p3 = plt.plot(xx, yy, 'k-', zorder=20)
            #plt.ylim(-0.2, 1.0)
            plt.legend([p1[0], p2[0], p3[0]],
                       ['image slice (y)', 'image slice (x)', 'fit'])
            plt.title('PSF fit')
            ps.savefig()

        return meas

    def update_astrometry(self, *args):
        return self.histogram_astrometric_shift(*args)

    def histogram_astrometric_shift(self, stars, wcs, fx, fy, trim_x0, trim_y0, pixsc, img, hdr,
                                    fullW, fullH, ps, printmsg):
        '''
        stars: reference stars
        wcs: header WCS
        fx,fy: detected stars
        '''
        import pylab as plt
        from astrometry.util.plotutils import dimshow, plothist
        ok,px,py = wcs.radec2pixelxy(stars.ra, stars.dec)
        px -= 1
        py -= 1
        
        if ps is not None:
            #kwa = dict(vmin=-3*sig1, vmax=50*sig1, cmap='gray')
            # Add to the 'detected sources' plot
            # mn,mx = np.percentile(img.ravel(), [50,99])
            # kwa = dict(vmin=mn, vmax=mx, cmap='gray')
            # plt.clf()
            # dimshow(img, **kwa)
            ax = plt.axis()
            #plt.plot(fx, fy, 'go', mec='g', mfc='none', ms=10)
            K = np.argsort(stars.mag)
            plt.plot(px[K[:10]]-trim_x0, py[K[:10]]-trim_y0, 'o', mec='m', mfc='none',ms=12,mew=2)
            plt.plot(px[K[10:]]-trim_x0, py[K[10:]]-trim_y0, 'o', mec='m', mfc='none', ms=8)
            plt.axis(ax)
            plt.title('PS1 stars')
            #plt.colorbar()
            ps.savefig()

        # we trimmed the image before running detection; re-add that margin
        fullx = fx + trim_x0
        fully = fy + trim_y0

        # Match PS1 to our detections, find offset
        radius = self.maxshift / pixsc

        I,J,dx,dy = self.match_ps1_stars(px, py, fullx, fully, radius, stars)
        #print(len(I), 'spatial matches with large radius', self.maxshift,
        #      'arcsec,', radius, 'pix')

        # Broad (5-pixel) bins are okay because we're going to refine them later... right?
        bins = 2*int(np.ceil(radius / 5.))
        print('Histogramming with', bins, 'bins')
        histo,xe,ye = np.histogram2d(dx, dy, bins=bins,
                                     range=((-radius,radius),(-radius,radius)))
        # smooth histogram before finding peak -- fuzzy matching
        from scipy.ndimage.filters import gaussian_filter
        histo = gaussian_filter(histo, 1.)
        histo = histo.T
        mx = np.argmax(histo)
        my,mx = np.unravel_index(mx, histo.shape)
        shiftx = (xe[mx] + xe[mx+1])/2.
        shifty = (ye[my] + ye[my+1])/2.
        
        if ps is not None:
            plt.clf()
            #plothist(dx, dy, range=((-radius,radius),(-radius,radius)), nbins=bins)
            plt.imshow(histo, interpolation='nearest', origin='lower', extent=[-radius, radius, -radius, radius],
                       cmap='hot')
            plt.xlabel('dx (pixels)')
            plt.ylabel('dy (pixels)')
            plt.title('Offsets to PS1 stars')
            ax = plt.axis()
            plt.axhline(0, color='b')
            plt.axvline(0, color='b')
            plt.plot(shiftx, shifty, 'o', mec='m', mfc='none', ms=15, mew=3)
            plt.axis(ax)
            ps.savefig()
    
        # Refine with smaller search radius
        radius2 = 3. / pixsc
        I,J,dx,dy = self.match_ps1_stars(px, py, fullx+shiftx, fully+shifty,
                                         radius2, stars)
        #print(len(J), 'matches to PS1 with small radius', 3, 'arcsec')
        shiftx2 = np.median(dx)
        shifty2 = np.median(dy)
        #print('Stage-1 shift', shiftx, shifty)
        #print('Stage-2 shift', shiftx2, shifty2)
        sx = shiftx + shiftx2
        sy = shifty + shifty2
        printmsg('Astrometric shift (%.0f, %.0f) pixels' % (sx,sy))

        # Find affine transformation
        cx = fullW/2
        cy = fullH/2

        A = np.zeros((len(J), 3))
        A[:,0] = 1.
        A[:,1] = (fullx[J] + sx) - cx
        A[:,2] = (fully[J] + sy) - cy

        R = np.linalg.lstsq(A, px[I] - cx, rcond=None)
        resx = R[0]
        #print('Affine transformation for X:', resx)
        R = np.linalg.lstsq(A, py[I] - cy, rcond=None)
        resy = R[0]
        #print('Affine transformation for Y:', resy)

        measargs = dict(affine = [cx, cy] + list(resx) + list(resy),
                        px=px-trim_x0-sx, py=py-trim_y0-sy,
                        dx=sx, dy=sy)

        subh,subw = img.shape
        wcs2 = wcs.get_subimage(sx, sy, subw, subh)

        #wcs2 = wcs.get_subimage(-trim_x0 - sx, -trim_y0 - sy, fullW, fullH)


        if True:
            r0,d0 = stars.ra[I[0]], stars.dec[I[0]]
            #print('Pan-STARRS RA,Dec:', r0,d0)
            ok,px0,py0 = wcs.radec2pixelxy(r0, d0)
            #print('Header WCS -> x,y:', px0-1, py0-1)
            mx0,my0 = fullx[J[0]], fully[J[0]]
            #print('VS Matched star x,y:', mx0, my0)
        
        if self.debug and ps is not None:
            plt.clf()
            plothist(dx, dy, range=((-radius2,radius2),(-radius2,radius2)))
            plt.xlabel('dx (pixels)')
            plt.ylabel('dy (pixels)')
            plt.title('Offsets to PS1 stars')
            ax = plt.axis()
            plt.axhline(0, color='b')
            plt.axvline(0, color='b')
            plt.plot(shiftx2, shifty2, 'o', mec='m', mfc='none', ms=15, mew=3)
            plt.axis(ax)
            ps.savefig()

        if ps is not None:
            afx = (fullx[J] + sx) - cx
            afy = (fully[J] + sy) - cy
            affx = cx + (resx[0] + resx[1] * afx + resx[2] * afy)
            affy = cy + (resy[0] + resy[1] * afx + resy[2] * afy)
            affdx = px[I] - affx
            affdy = py[I] - affy
            plt.clf()
            plothist(affdx, affdy, range=((-radius2,radius2),(-radius2,radius2)))
            plt.xlabel('dx (pixels)')
            plt.ylabel('dy (pixels)')
            plt.title('Offsets to PS1 stars after Affine correction')
            ax = plt.axis()
            plt.axhline(0, color='b')
            plt.axvline(0, color='b')
            #plt.plot(shiftx, shifty, 'o', mec='m', mfc='none', ms=15, mew=3)
            plt.axis(ax)
            ps.savefig()
        
        if ps is not None:
            mn,mx = np.percentile(img.ravel(), [50,99])
            kwa2 = dict(vmin=mn, vmax=mx, cmap='gray')
            plt.clf()
            dimshow(img, **kwa2)
            ax = plt.axis()
            plt.plot(fx[J], fy[J], 'go', mec='g', mfc='none', ms=10, mew=2)
            plt.plot(px[I]-sx-trim_x0, py[I]-sy-trim_y0, 'm+', ms=10, mew=2)
            plt.axis(ax)
            plt.title('Matched PS1 stars')
            plt.colorbar()
            ps.savefig()

            plt.clf()
            dimshow(img, **kwa2)
            ax = plt.axis()
            plt.plot(fx[J], fy[J], 'go', mec='g', mfc='none', ms=10, mew=2)
            K = np.argsort(stars.mag)
            plt.plot(px[K[:10]]-sx-trim_x0, py[K[:10]]-sy-trim_y0, 'o', mec='m', mfc='none',ms=12,mew=2)
            plt.plot(px[K[10:]]-sx-trim_x0, py[K[10:]]-sy-trim_y0, 'o', mec='m', mfc='none', ms=8,mew=2)
            plt.axis(ax)
            plt.title('All PS1 stars')
            plt.colorbar()
            ps.savefig()

            plt.clf()
            dimshow(img, **kwa2)
            ax = plt.axis()
            ok,px2,py2 = wcs2.radec2pixelxy(stars.ra, stars.dec)
            px2 -= 1
            py2 -= 1
            plt.plot(fx, fy, 'go', mec='g', mfc='none', ms=10, mew=2)
            plt.plot(px2, py2, 'mo', mec='m', mfc='none', ms=8, mew=2)
            plt.axis(ax)
            plt.title('All PS1 stars (with New WCS')
            plt.colorbar()
            ps.savefig()


        return wcs2, measargs


class DECamMeasurer(RawMeasurer):
    def __init__(self, *args, **kwargs):
        if not 'pixscale' in kwargs:
            import decam
            kwargs.update(pixscale = decam.decam_nominal_pixscale)
        super(DECamMeasurer, self).__init__(*args, **kwargs)
        self.camera = 'decam'

    def read_raw(self, F, ext):
        return read_raw_decam(F, ext)

    def get_sky_and_sigma(self, img):
        sky,sig1 = sensible_sigmaclip(img[1500:2500, 500:1000])
        return sky,sig1

    def get_wcs(self, hdr):
        from astrometry.util.util import wcs_pv2sip_hdr
        from astrometry.util.starutil_numpy import hmsstring2ra, dmsstring2dec
        # In DECam images, the CRVAL values are set to TELRA and
        # TELDEC, but for computing the RA,Dec offsets, we want to use
        # RA and DEC (commanded position, not what the telescope
        # pointing model says).
        # HACK -- convert TPV WCS header to SIP.
        oldcrval1 = hdr['CRVAL1']
        oldcrval2 = hdr['CRVAL2']
        hdr['CRVAL1'] = hmsstring2ra (self.primhdr['RA' ])
        hdr['CRVAL2'] = dmsstring2dec(self.primhdr['DEC'])
        # print(('Updated WCS CRVAL from TELRA,TELDEC = (%.4f,%.4f) to ' +
        #       'RA,Dec = (%.4f,%.4f)') %
        #       (oldcrval1, oldcrval2, hdr['CRVAL1'], hdr['CRVAL2']))
        wcs = wcs_pv2sip_hdr(hdr)
        hdr['CRVAL1'] = oldcrval1
        hdr['CRVAL2'] = oldcrval2
        #print('Converted WCS to', wcs)
        return wcs

    def colorterm_ps1_to_observed(self, ps1stars, band):
        return ps1_to_decam(ps1stars, band)

    def get_ps1_band(self, band):
        return ps1cat.ps1band[band]

class DECamCPMeasurer(DECamMeasurer):
    def read_raw(self, F, ext):
        img = F[ext].read()
        hdr = F[ext].read_header()
        img = img.astype(np.float32)
        return img,hdr

class Mosaic3Measurer(RawMeasurer):
    def __init__(self, *args, **kwargs):
        if not 'pixscale' in kwargs:
            import mosaic
            kwargs.update(pixscale = mosaic.mosaic_nominal_pixscale)
        super(Mosaic3Measurer, self).__init__(*args, **kwargs)
        self.camera = 'mosaic3'

    def get_band(self, primhdr):
        band = super(Mosaic3Measurer,self).get_band(primhdr)
        # "zd" -> "z"
        #return band[0]
        return band

    def read_raw(self, F, ext):
        '''
        F: fitsio FITS object
        '''
        img = F[ext].read()
        hdr = F[ext].read_header()
        img = img.astype(np.float32)
        #print('img qts', np.percentile(img.ravel(), [0,25,50,75,100]))
        # Subtract median overscan and multiply by gains 
        dataA = parse_section(hdr['DATASEC'], slices=True)
        biasA = parse_section(hdr['BIASSEC'], slices=True)
        gainA = hdr['GAIN']
        b = np.median(img[biasA])
        #print('subtracting bias', b)
        img[dataA] = (img[dataA] - b) * gainA
    
        # Trim the image
        trimA = parse_section(hdr['TRIMSEC'], slices=True)
        # zero out all but the trim section
        trimg = img[trimA].copy()
        img[:,:] = 0
        img[trimA] = trimg
        return img,hdr

    def get_sky_and_sigma(self, img):
        # Spline sky model to handle (?) ghost / pupil?
        from tractor.splinesky import SplineSky

        splinesky = SplineSky.BlantonMethod(img, None, 256)
        skyimg = np.zeros_like(img)
        splinesky.addTo(skyimg)
        
        mnsky,sig1 = sensible_sigmaclip(img - skyimg)
        return skyimg,sig1

    def remove_sky_gradients(self, img):
        pass

    def get_wcs(self, hdr):
        # Older images have ZPX, newer TPV.
        if hdr['CTYPE1'] == 'RA---TPV':
            from astrometry.util.util import wcs_pv2sip_hdr
            wcs = wcs_pv2sip_hdr(hdr)
        else:
            from astrometry.util.util import Tan
            hdr['CTYPE1'] = 'RA---TAN'
            hdr['CTYPE2'] = 'DEC--TAN'
            wcs = Tan(hdr)
        return wcs

    def colorterm_ps1_to_observed(self, ps1stars, band):
        bandmap = dict(zd='z', D51='g')
        band = bandmap.get(band, band)
        return ps1_to_mosaic(ps1stars, band)

    def get_ps1_band(self, band):
        bandmap = dict(zd='z', D51='g')
        band = bandmap.get(band, band)
        return ps1cat.ps1band[band]

class PointingCamMeasurer(RawMeasurer):

    # full image size: 3296 x 2472

    def __init__(self, *args, ps=None, verbose=None, **kwargs):
        super(PointingCamMeasurer, self).__init__(*args, **kwargs)
        self.camera = 'pointing'

        self.subW = 3296 // 2
        self.subH = 2472 // 2

    def get_sky_and_sigma(self, img):
        #sky,sig1 = sensible_sigmaclip(img[1000:2000, 500:1500])
        # First quadrant
        sky,sig1 = sensible_sigmaclip(img[100:self.subH-100, 100:self.subW-100])
        print('Sky, sig1:', sky, sig1)
        return sky,sig1

    def get_wcs(self, hdr):
        from astrometry.util.util import Tan
        wcs = Tan(hdr)
        wcs = wcs.get_subimage(0, 0, self.subW, self.subH)
        return wcs

    def read_raw(self, F, ext):
        img = F[ext].read()
        hdr = F[ext].read_header()
        img = img.astype(np.float32)
        img = img[:self.subH, :self.subW]
        return img, hdr

    def get_reference_stars(self, wcs, band):
        # Use a cut version of the Gaia catalog.
        cat = GaiaCatalog()
        stars = cat.get_catalog_in_wcs(wcs)
        print('Got', len(stars), 'Gaia stars')
        if not 'topring_mag' in stars.get_columns():
            raise RuntimeError('Need a specially cut Gaia catalog for the pointing cam')
        stars.mag = stars.topring_mag
        stars.cut(stars.mag < 14)
        print('Cut at 14th mag:', len(stars))
        return stars

    def cut_reference_catalog(self, stars):
        return np.ones(len(stars), bool)

    def get_color_term(self, stars, band):
        return 0.

    def get_exptime(self, primhdr):
        exptime = primhdr.get('EXPTIME', 0.)
        if exptime == 0.:
            exptime = primhdr.get('EXPOSURE', 0.) / 1000.
        return exptime

class DesiCiMeasurer(RawMeasurer):
    def __init__(self, *args, ps=None, verbose=None, **kwargs):
        super(DesiCiMeasurer, self).__init__(*args, **kwargs)
        self.camera = 'desici'
        self.pixscale = 0.130
        self.edge_trim = 0

    def get_band(self, primhdr):
        # HACK
        return 'r'

    def get_sky_and_sigma(self, img):
        sky,sig1 = sensible_sigmaclip(img[100:-100, 100:-100])
        print('Sky, sig1:', sky, sig1)
        return sky,sig1

    def update_astrometry(self, stars, wcs, fx, fy, trim_x0, trim_y0,
                          pixsc, img, hdr, fullW, fullH, ps, printmsg):
        import tempfile
        from astrometry.util.util import healpix_rangesearch_radec
        # Run Astrometry.net on cleaned-up image.
        tempimg = tempfile.NamedTemporaryFile(suffix='.fits')
        configfile = tempfile.NamedTemporaryFile(suffix='.cfg')
        fn = tempimg.name
        fitsio.write(fn, img, header=hdr, clobber=True)

        ra = hdr['CRVAL1']
        dec = hdr['CRVAL2']
        radius = 5.
        nside = 2
        healpixes = healpix_rangesearch_radec(ra, dec, radius, nside)

        configfn = configfile.name
        f = open(configfn, 'w')
        index_dir = os.environ.get('AN_INDEX_DIR', 'index-files')
        # export AN_INDEX_DIR=/global/project/projectdirs/cosmo/work/users/dstn/index-5000
        f.write('add_path %s\n' % index_dir +
                'inparallel\n' +
                '\n'.join(['index index-5001-%02i' % hp for hp in healpixes]) + '\n' +
                '\n'.join(['index index-5002-%02i' % hp for hp in healpixes]) + '\n')
        f.close()

        cmd = 'solve-field --config %s --continue' % configfn
        cmd += ' --ra %f --dec %f --radius 5' % (ra, dec)
        if self.ext == 'CIC':
            pixsc = 0.133
        else:
            pixsc = 0.123
            cmd += ' --xscale 1.09'
        cmd += ' --scale-low %f --scale-high %f --scale-units app' % (pixsc * 0.95, pixsc * 1.05)
        cmd += ' --objs 50 --crpix-center'
        cmd += ' --parity pos'
        cmd += ' --no-verify --new-fits none --solved none --match none --rdls none --corr none'
        #cmd += ' --plot-scale 0.25'
        cmd += ' --no-plots'
        cmd += ' --downsample 4'
        #cmd += ' -v'
        cmd += ' ' + fn
        print(cmd)
        os.system(cmd)

        from astrometry.util.util import Sip, Tan
        wcsfn = fn.replace('.fits', '.wcs')
        if os.path.exists(wcsfn):
            wcs2 = Sip(wcsfn)

            imh,imw = img.shape
            cx,cy = (imw + 1.) / 2., (imh + 1.) / 2.
            r,d = wcs2.pixelxy2radec(cx, cy)
            ok,x1,y1 = wcs.radec2pixelxy(r, d)
            dx = x1 - cx - trim_x0
            dy = y1 - cy - trim_y0
            # print('Computing dx,dy from initial WCS:')
            # print(wcs)
            # print('to')
            # print(wcs2)
            # print('Final WCS image center:', cx,cy)
            # print('-> RA,Dec', r,d)
            # print('-> initial WCS coords', x1,y1)
            # print('-> dx,dy', dx, dy)

            r1,d1 = wcs.radec_center()
            r2,d2 = wcs2.radec_center()
            dd = d2 - d1
            dr = (r2 - r1)*np.cos(np.deg2rad(d2))
            #print('dr,dd', dr*3600, dd*3600, 'arcsec')

            measargs = dict(dx=dx, dy=dy, dra=dr, ddec=dd)
            return wcs2, measargs
        print('Solving failed!  Trying histogram...')
        return self.histogram_astrometric_shift(stars, wcs, fx, fy, trim_x0, trim_y0, pixsc,
                                                img, hdr, fullW, fullH, ps, printmsg)

    def get_wcs(self, hdr):
        from astrometry.util.util import Tan
        # Needs extremely recent Astrometry.net
        wcs = Tan(hdr)
        return wcs

    def get_bias(self, ext):
        # From Aaron:
        # https://github.com/ameisner/ci_reduce/blob/master/py/ci_reduce/common.py#L315
        # result in in ADU
        bias_med_dict = {'CIE' : 991.0, 
                         'CIN' : 1013.0, 
                         'CIC' : 1020.0, 
                         'CIS' : 984.0, 
                         'CIW' : 990.0}
        return bias_med_dict.get(ext, 0.)

    def read_raw(self, F, ext):
        img = F[ext].read()
        hdr = F[ext].read_header()
        img = img.astype(np.float32)

        img -= self.get_bias(ext)

        # Estimate per-pixel noise via Blanton's 5-pixel MAD
        from scipy.ndimage.measurements import label, find_objects
        slice1 = (slice(0,-5,10),slice(0,-5,10))
        slice2 = (slice(5,None,10),slice(5,None,10))
        mad = np.median(np.abs(img[slice1] - img[slice2]).ravel())
        sig1 = 1.4826 * mad / np.sqrt(2.)
        #print('Computed sig1 by Blanton method:', sig1)
        med = np.median(img.ravel())
        #print('Median:', med)
        thresh = med + 5. * sig1
        blobs,nblobs = label(img > thresh)
        #print('Found', nblobs, 'blobs of pixels above threshold')
        nz = 0
        for sy,sx in find_objects(blobs):
            y0 = sy.start
            y1 = sy.stop
            if y1 > y0+1:
                continue
            x0 = sx.start
            x1 = sx.stop
            if x1 > x0+1:
                continue
            nz += 1
            img[y0,x0] = med
        print('Zeroed out', nz, 'pixels')
        return img, hdr

    def get_ps1_band(self, band):
        return ps1cat.ps1band[band]

    def colorterm_ps1_to_observed(self, ps1stars, band):
        return 0.

class BokMeasurer(RawMeasurer):
    def __init__(self, *args, **kwargs):
        if not 'pixscale' in kwargs:
            import mosaic
            kwargs.update(pixscale = mosaic.mosaic_nominal_pixscale)
        super(BokMeasurer, self).__init__(*args, **kwargs)
        self.camera = '90prime'

    def get_band(self, primhdr):
        band = super(BokMeasurer,self).get_band(primhdr)
        # bokr -> r
        band = band.strip()
        band = band.replace('bok', '')
        print('Band', band)
        return band[0]

    def read_raw(self, F, ext):
        '''
        F: fitsio FITS object
        '''
        img = F[ext].read()
        hdr = F[ext].read_header()
        img = img.astype(np.float32)
        #print('img qts', np.percentile(img.ravel(), [0,25,50,75,100]))

        primhdr = F[0].read_header()

        # Subtract median overscan and multiply by gains 
        data = parse_section(hdr['DATASEC'], slices=True)
        bias = parse_section(hdr['BIASSEC'], slices=True)
        gain = primhdr['GAIN%i' % int(ext.replace('IM','').replace('im',''),10)]
        b = np.median(img[bias])
        #print('subtracting bias', b)
        img[data] = (img[data] - b) * gain
    
        # Trim the image
        trim = parse_section(hdr['TRIMSEC'], slices=True)
        # zero out all but the trim section
        trimg = img[trim].copy()
        img[:,:] = 0
        img[trim] = trimg
        return img,hdr

    def get_wcs(self, hdr):
        from astrometry.util.util import Tan
        wcs = Tan(hdr)
        return wcs

    def get_sky_and_sigma(self, img):
        # Spline sky model to handle (?) ghost / pupil?
        from tractor.splinesky import SplineSky
        splinesky = SplineSky.BlantonMethod(img, None, 256)
        skyimg = np.zeros_like(img)
        splinesky.addTo(skyimg)
        mnsky,sig1 = sensible_sigmaclip(img - skyimg)
        return skyimg,sig1

    def remove_sky_gradients(self, img):
        pass

    def colorterm_ps1_to_observed(self, ps1stars, band):
        return ps1_to_90prime(ps1stars, band)

    def get_ps1_band(self, band):
        return ps1cat.ps1band[band]


class BokCPMeasurer(BokMeasurer):
    def read_raw(self, F, ext):
        img = F[ext].read()
        hdr = F[ext].read_header()
        img = img.astype(np.float32)
        return img,hdr

    def get_wcs(self, hdr):
        from astrometry.util.util import wcs_pv2sip_hdr
        wcs = wcs_pv2sip_hdr(hdr)
        return wcs

    def zeropoint_for_exposure(self, band, ext=None, exptime=None, **kwargs):
        zp0 = super(BokCPMeasurer, self).zeropoint_for_exposure(band, ext=ext, exptime=exptime, **kwargs)
        if zp0 is None:
            return zp0
        zp0 -= 2.5 * np.log10(exptime)
        return zp0

def measure_raw_decam(fn, ext='N4', nom=None, ps=None, measargs={}, **kwargs):
    '''
    Reads the given file *fn*, extension *ext*', and measures a number of
    quantities that are useful for depth estimates:

    Returns:
    --------
    dict containing:

    band : string, eg 'g'

    airmass : float

    seeing : float, FWHM in arcsec.  From Gaussian fit to PS1 stars.

    zp : float.  From comparison to PS1, using aperture photometry.

    skybright : float, magnitudes/arcsec^2, based on a sky estimate and the
        canonical zeropoint.

    transparency : float, fraction.  Atmospheric transparency, based
        on computed zeropoint versus canonical zeropoint.


    Details:
    --------

    What this function does:
        
    - read raw image, apply bias and gain to the two amps & trim

    - estimate sky level via sigma-clip of central pixels
    (y in [1500,255), x in [500,1000); this is like decstat.pro)

    -> sky brightness, assuming canonical zeropoint

    - hack: remove sky gradients (or flat / gain / bias problems) by
    subtracting first a column-wise and then row-wise median from the
    image.

    - hack: remove a 100-pixel margin

    - assume a FWHM of 5 pixels, smooth by a gaussian and detect sources
    above 20 sigma.  Within each 'blob' of pixels above threshold, take
    the max pixel as the peak, and then compute the centroid in +- 5
    pixel box around the peak.  Keep sources where the centroid is
    within 1 pixel of the peak.

    - aperture-photometer peaks with a radius of 3.5"

    - I also tried computing a sky annulus, but it doesn't seem to work
    better, so it's not used

    - Read PS1 stars, apply color cut for stars, apply color
      correction to DECam

      - Match PS1 stars to our peaks

      - Naively find dx,dy offset (via histogram peak)

      - Sigma-clip to find mag offset PS1 to our apmags.

      -> zeropoint & transparency (vs canonical zeropoint)

      - For each matched star, cut out a +- 15 pixel box and use the
      Tractor to fit a single-Gaussian PSF

      - Take the median of the fit FWHMs -> FWHM estimate

    '''
    if nom is None:
        import decam
        nom = decam.DecamNominalCalibration()
    meas = DECamMeasurer(fn, ext, nom, **measargs)
    results = meas.run(ps, **kwargs)
    return results

def measure_raw_mosaic3(fn, ext='im4', nom=None, ps=None,
                        measargs={}, **kwargs):
    if nom is None:
        import mosaic
        nom = mosaic.MosaicNominalCalibration()
    meas = Mosaic3Measurer(fn, ext, nom, **measargs)
    results = meas.run(ps, **kwargs)
    return results

def camera_name(primhdr):
    '''
    Returns 'mosaic3' or 'decam' or '90prime'
    '''
    return primhdr.get('INSTRUME','').strip().lower()

def get_measurer_class_for_file(fn):
    primhdr = fitsio.read_header(fn)
    cam = camera_name(primhdr)
    #print('Camera:', cam)
    if cam == 'mosaic3':
        return Mosaic3Measurer
    elif cam == '90prime':
        if 'PLVER' in primhdr:
            # CP-processed Bok image
            return BokCPMeasurer
        return BokMeasurer
    elif cam == 'decam':
        if 'PLVER' in primhdr:
            # CP-processed DECam image
            return DECamCPMeasurer
        return DECamMeasurer

def measure_raw(fn, primext=0, **kwargs):
    primhdr = fitsio.read_header(fn, ext=primext)
    cam = camera_name(primhdr)

    if cam == 'mosaic3':
        return measure_raw_mosaic3(fn, **kwargs)
    elif cam == '90prime':
        nom = kwargs.pop('nom', None)
        if nom is None:
            import bok
            nom = bok.BokNominalCalibration()

        if 'PLVER' in primhdr:
            # CP-processed Bok image
            measclass = BokCPMeasurer
            ext = kwargs.pop('ext', 'ccd1')
        else:
            measclass = BokMeasurer
            ext = kwargs.pop('ext', 'IM4')

        meas = measclass(fn, ext, nom)
        ps = kwargs.pop('ps', None)
        results = meas.run(ps, **kwargs)
        return results
    elif cam == 'decam':
        if 'PLVER' in primhdr:
            # CP-processed DECam image
            import decam
            nom = decam.DecamNominalCalibration()
            ext = kwargs.pop('ext')
            meas = DECamCPMeasurer(fn, ext, nom)
            results = meas.run(**kwargs)
            return results
        return measure_raw_decam(fn, **kwargs)

    elif cam == 'pointing':
        #return measure_raw_pointing(fn, **kwargs)
        nom = kwargs.pop('nom', None)
        if nom is None:
            import camera_pointing
            nom = camera_pointing.PointingCamNominalCalibration()
        ext = kwargs.pop('ext', 0)
        meas = PointingCamMeasurer(fn, ext, nom, **kwargs)
        results = meas.run(**kwargs)
        return results

    elif cam == 'desi':
        nom = kwargs.pop('nom', None)
        if nom is None:
            import camera_desici
            nom = camera_desici.DesiCiNominalCalibration()
        ext = kwargs.pop('ext', 0)
        meas = DesiCiMeasurer(fn, ext, nom, **kwargs)
        kwargs.update(primext=primext)
        results = meas.run(**kwargs)
        return results
    
    return None

def get_default_extension(fn):
    primhdr = fitsio.read_header(fn)
    cam = camera_name(primhdr)

    ### HACK -- must match default of measure_raw_* above
    if cam == 'mosaic3':
        return 'im4'
    elif cam == 'decam':
        return 'N4'

def sensible_sigmaclip(arr, nsigma = 4.):
    from scipy.stats import sigmaclip
    goodpix,lo,hi = sigmaclip(arr, low=nsigma, high=nsigma)
    # sigmaclip returns unclipped pixels, lo,hi, where lo,hi are
    # mean(goodpix) +- nsigma * sigma
    meanval = np.mean(goodpix)
    sigma = (meanval - lo) / nsigma
    return meanval, sigma

def parse_section(s, slices=False):
    '''
    parse '[57:2104,51:4146]' into integers; also subtract 1.
    '''
    s = s.replace('[','').replace(']','').replace(',',' ').replace(':',' ')
    #print('String', s)
    i = [int(si)-1 for si in s.split()]
    assert(len(i) == 4)
    if not slices:
        return i
    slc = slice(i[2], i[3]+1), slice(i[0], i[1]+1)
    #print('Slice', slc)
    return slc

def read_raw_decam(F, ext):
    '''
    F: fitsio FITS object
    '''
    img = F[ext].read()
    hdr = F[ext].read_header()
    #print('Image type', img.dtype, img.shape)

    img = img.astype(np.float32)
    #print('Converted image to', img.dtype, img.shape)

    if False:
        mn,mx = np.percentile(img.ravel(), [25,95])
        kwa = dict(vmin=mn, vmax=mx)

        plt.clf()
        dimshow(img, **kwa)
        plt.title('Raw image')
        ps.savefig()

        M = 200
        plt.clf()
        plt.subplot(2,2,1)
        dimshow(img[-M:, :M], ticks=False, **kwa)
        plt.subplot(2,2,2)
        dimshow(img[-M:, -M:], ticks=False, **kwa)
        plt.subplot(2,2,3)
        dimshow(img[:M, :M], ticks=False, **kwa)
        plt.subplot(2,2,4)
        dimshow(img[:M, -M:], ticks=False, **kwa)
        plt.suptitle('Raw corners')
        ps.savefig()
    
    if 'DESBIAS' in hdr:
        assert(False)

    # DECam RAW image

    # Raw image size 2160 x 4146

    # Subtract median overscan and multiply by gains 
    # print('DATASECA', hdr['DATASECA'])
    # print('BIASSECA', hdr['BIASSECA'])
    # print('DATASECB', hdr['DATASECB'])
    # print('BIASSECB', hdr['BIASSECB'])
    dataA = parse_section(hdr['DATASECA'], slices=True)
    biasA = parse_section(hdr['BIASSECA'], slices=True)
    dataB = parse_section(hdr['DATASECB'], slices=True)
    biasB = parse_section(hdr['BIASSECB'], slices=True)
    gainA = hdr['GAINA']
    gainB = hdr['GAINB']
    # print('DataA', dataA)
    # print('BiasA', biasA)
    # print('DataB', dataB)
    # print('BiasB', biasB)

    if False:
        plt.clf()
        plt.plot([np.median(img[i,:] )     for i in range(100)], 'b-')
        plt.plot([np.median(img[-(i+1),:]) for i in range(100)], 'c-')
        plt.plot([np.median(img[:,i] )     for i in range(100)], 'r-')
        plt.plot([np.median(img[:,-(i+1)]) for i in range(100)], 'm-')
        plt.title('Img')
        ps.savefig()
    
        plt.clf()
        plt.plot([np.median(img[dataA][i,:] )     for i in range(100)], 'b-')
        plt.plot([np.median(img[dataA][-(i+1),:]) for i in range(100)], 'c-')
        plt.plot([np.median(img[dataA][:,i] )     for i in range(100)], 'r-')
        plt.plot([np.median(img[dataA][:,-(i+1)]) for i in range(100)], 'm-')
        plt.title('Img DataA')
        ps.savefig()
    
        plt.clf()
        plt.plot([np.median(img[dataB][i,:] )     for i in range(100)], 'b-')
        plt.plot([np.median(img[dataB][-(i+1),:]) for i in range(100)], 'c-')
        plt.plot([np.median(img[dataB][:,i] )     for i in range(100)], 'r-')
        plt.plot([np.median(img[dataB][:,-(i+1)]) for i in range(100)], 'm-')
        plt.title('Img DataB')
        ps.savefig()

    
    img[dataA] = (img[dataA] - np.median(img[biasA])) * gainA
    img[dataB] = (img[dataB] - np.median(img[biasB])) * gainB
    
    # Trim the image -- could just take the min/max of TRIMSECA/TRIMSECB...
    trimA = parse_section(hdr['TRIMSECA'], slices=True)
    trimB = parse_section(hdr['TRIMSECB'], slices=True)
    # copy the TRIM A,B sections into a new image...
    trimg = np.zeros_like(img)
    trimg[trimA] = img[trimA]
    trimg[trimB] = img[trimB]
    # ... and then cut that new image
    trim = parse_section(hdr['TRIMSEC'], slices=True)
    img = trimg[trim]
    #print('Trimmed image:', img.dtype, img.shape)

    return img,hdr
    



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script to make measurements on raw DECam images to estimate sky brightness, PSF size, and zeropoint / transparency for exposure-time scaling.')
    parser.add_argument('--', '-o', dest='outfn', default='raw.fits',
                        help='Output file name')
    parser.add_argument('--ext', help='Extension to read for computing observing conditions: default "N4" for DECam, "im4" for Mosaic3', default=None)
    #    parser.add_argument('--ext', default=[], action='append',
    #    help='FITS image extension to read: default "N4"; may be repeated')
    parser.add_argument('--aprad', help='Aperture photometry radius in arcsec',
                        type=float, default=None)
    
    parser.add_argument('--plots', help='Make plots with this filename prefix',
                        default=None)
    parser.add_argument('images', metavar='DECam-filename.fits.fz', nargs='+',
                        help='DECam raw images to process')

    args = parser.parse_args()
    fns = args.images
    #exts = args.ext
    #if len(exts) == 0:
    #    exts = ['N4']

    ext = args.ext
    
    ps = None
    if args.plots is not None:
        from astrometry.util.plotutils import PlotSequence
        ps = PlotSequence(args.plots)

    measargs = dict()
    if args.aprad is not None:
        measargs.update(aprad=args.aprad)
        
    vals = {}
    for fn in fns:
        #for ext in exts:
        print()
        print('Measuring', fn, 'ext', ext)
        d = measure_raw(fn, ext=ext,ps=ps, measargs=measargs)
        d.update(filename=fn, ext=ext)
        for k,v in d.items():
            if not k in vals:
                vals[k] = [v]
            else:
                vals[k].append(v)

    from astrometry.util.fits import fits_table

    T = fits_table()
    for k,v in vals.items():
        if k == 'wcs':
            continue
        T.set(k, v)
    T.to_np_arrays()
    for k in T.columns():
        v = T.get(k)
        if v.dtype == np.float64:
            T.set(k, v.astype(np.float32))
    T.writeto(args.outfn)
    print('Wrote', args.outfn)
