from __future__ import print_function
if __name__ == '__main__':
    import matplotlib
    matplotlib.use('Agg')
from astrometry.util.plotutils import *
from measure_raw import *

import numpy as np
import pylab as plt

class Mosaic3FocusMeas(Mosaic3Measurer):

    def get_focus_shifts(self, hdr):
        '''
        Returns an array of pixel-space shifts for this focus
        sequence.  These are uniform except for the last one, which is
        doubled.
        '''
        nsteps = hdr['FOCNEXPO']
        steppix = hdr['FOCSHIFT']
        steps = np.arange(nsteps)
        steps[-1] += 1
        return steps * steppix

    def get_focus_values(self, hdr):
        '''
        Returns an array of focus-space (microns) steps for this focus sequence.
        It is aligned with the array from get_focus_shifts.
        '''
        nsteps = hdr['FOCNEXPO']
        step = hdr['FOCSTEP']
        f0 = hdr['FOCSTART']
        return f0 + step * np.arange(nsteps)
    
    def detection_map(self, img, sig1, psfsig, ps):
        from scipy.ndimage.filters import gaussian_filter
        # Compute detection map
        psfnorm = 1./(2. * np.sqrt(np.pi) * psfsig)
        detsn = gaussian_filter(img / sig1, psfsig) / psfnorm
        
        # Take the *minimum* of the detmap and shifted versions of itself.
        minimg = detsn.copy()
        shifts = self.get_focus_shifts(self.hdr)
        for dy in -shifts:
            if dy == 0:
                continue
            elif dy > 0:
                minimg[dy:,:] = np.minimum(minimg[dy:,:], detsn[:-dy,:])
                minimg[:dy,:] = 0
            elif dy < 0:
                minimg[:dy,:] = np.minimum(minimg[:dy,:], detsn[-dy:,:])
                minimg[dy:,:] = 0

        if ps is not None:
            kwa = dict(vmin=-3, vmax=20, cmap='gray')
            plt.clf()
            dimshow(minimg, **kwa)
            plt.title('Min image (detection S/N)')
            plt.colorbar()
            ps.savefig()

        detsn = minimg
        # zero out the edges -- larger margin here?
        detsn[0 ,:] = 0
        detsn[:, 0] = 0
        detsn[-1,:] = 0
        detsn[:,-1] = 0
        return detsn

    def trim_edges(self, img):
        # Trim off some edge pixels.
        trim = self.edge_trim

        # In the focus images, the bottom ~300 pixels are funky when shift<0
        # Top 300 pixels when shift>0
        hdr = self.hdr
        steppix = hdr['FOCSHIFT']
        if steppix < 0:
            bottom = 300
            cimg = img[bottom:-trim, trim:-trim]
            return cimg, trim, bottom
        else:
            top = 300
            cimg = img[trim:-top, trim:-trim]
            return cimg, trim, trim

    def run(self, ps=None, plotfn=None):
        from astrometry.libkd.spherematch import match_xy
        from legacyanalysis.ps1cat import ps1cat
        import photutils
        meas = super(Mosaic3FocusMeas, self).run(ps=ps, focus=True)
        # momentsize=10)

        px = meas['px']
        py = meas['py']
        fx = meas['fx']
        fy = meas['fy']
        img= meas['img']
        sig1=meas['sig1']
        stars=meas['stars']
        band=meas['band']
        apflux = meas['apflux']
        
        # (mx2 ,my2 ,mxy ,mtheta,ma,mb,mell) = meas['moments' ]
        # (wmx2,wmy2,wmxy,wtheta,wa,wb,well) = meas['wmoments']
        
        radius2 = 3. / self.pixscale
        I,J,d = match_xy(px, py, fx, fy, radius2)

        # Choose the brightest PS1 stars not near an edge
        ps1band = ps1cat.ps1band.get(band, 2)
        ps1mag = stars.median[I, ps1band]
        hdr = self.hdr
        shifts = self.get_focus_shifts(hdr)
        S = 20
        H,W = img.shape
        miny = fy[J] + min(shifts)
        maxy = fy[J] + max(shifts)
        edge = reduce(np.logical_or, [fx[J] < S, fx[J]+S >= W,
                                      miny  < S, maxy +S >= H])
        K = np.argsort(ps1mag + 1000*edge)

        nstars = min(10, len(K))
        if ps is not None:
            plt.clf()
            plt.subplots_adjust(hspace=0, wspace=0)
            sp = 1
            mn,mx = np.percentile(img.ravel(), [50,99])
            kwa = dict(vmin=mn, vmax=mx, cmap='gray')
            for shifty in shifts:
                for k in K[:nstars]:
                    plt.subplot(len(shifts), nstars, sp)
                    sp += 1
                    i = I[k]
                    j = J[k]
                    x = fx[j]
                    y = fy[j] + shifty
                    dimshow(img[y-S:y+S+1, x-S:x+S+1], ticks=False, **kwa)
                    #if shifty == 0.:
                    #    print('Plotting focus sweep star at', x,y)
            plt.suptitle('Focus sweep stars')
            ps.savefig()

        meas.update(I=I, J=J, K=K)

        shifts = self.get_focus_shifts(hdr)
        focus = self.get_focus_values(hdr)
    
        allfocus = []
        allcxx,allcyy,allcxy = [],[],[]
        nstars = 21
        for ik,k in enumerate(K[:nstars]):
            j = J[k]
            xi = fx[j]
            yi = fy[j]
            fluxi = apflux[j]
            #print('Fitting star at', xi,yi)
            for shifty,foc in zip(shifts, focus):
                j = J[k]
                xi = fx[j]
                yi = fy[j] + shifty
                p = self.fit_general_gaussian(img, sig1, xi, yi, fluxi,
                                              ps=ps if ik ==0 else None)
                #print('PSF variance:', p)
                allfocus.append(foc)
                allcxx.append(p[0])
                allcyy.append(p[1])
                allcxy.append(p[2])
    
        allfocus  = np.array(allfocus)
        allcxx = np.array(allcxx)
        allcyy = np.array(allcyy)
        allcxy = np.array(allcxy)

        meas.update(shifts=shifts, focus=focus,
                    allfocus=allfocus, allcxx=allcxx, allcyy=allcyy,
                    allcxy=allcxy)

        FF = []

        # Fit XX and YY covariances as a quadratic function of focus.
        names = ('PSF X fwhm (arcsec)', 'PSF Y fwhm (arcsec)')
        fitvals = []
        for name,Y in zip(names, [allcxx, allcyy]):
            X,I = np.unique(allfocus, return_inverse=True)
            Ymn,Ysig = np.zeros(len(X)), np.zeros(len(X))
            for i in range(len(X)):
                J = np.flatnonzero(I == i)
                yi = Y[J]
                y1,ymn,y2 = np.percentile(yi, [16,50,84])
                ysig = (y1 - y2) / 2
                Ymn [i] = ymn
                Ysig[i] = ysig

            xmed = np.median(X)
            dx = (np.max(X) - np.min(X)) / 2.
            # Rescale xx in [-1,1]
            xx = (X - xmed) / dx
            A = np.zeros((len(xx),3))
            A[:,0] = 1.
            A[:,1] = xx
            A[:,2] = xx**2
            wt = 1. / Ysig
            A *= wt[:,np.newaxis]
            b = Ymn * wt
            R = np.linalg.lstsq(A, b)
            s = R[0]
            #print('Least-squares quadratic:', s)
            mn,mx = allfocus.min(),allfocus.max()
            d = mx-mn
            mn -= 0.1 * d
            mx += 0.1 * d
            xx = np.linspace(mn, mx, 500)
            rx = (xx - xmed)/dx
            qq = s[0] + rx*s[1] + rx**2*s[2]
            # minimum
            fbest = -s[1] / (2. * s[2])
            fbest = fbest * dx + xmed
            fitvals.append((X, Ymn, Ysig, s, xx, qq, fbest))
            print('Focus position from %s: %.1f' % (name, fbest))
            FF.append(fbest)

        fmean = np.mean(FF)
        print('Mean focus: %.1f' % fmean)

        if plotfn is None and ps is None:
            return meas
        
        if plotfn is not None:
            plt.clf()

        seeings = []

        def pixvar2seeing(var):
            return np.sqrt(var) * 2.35 * self.pixscale
        
        for i,(name,Y,(X,Ymn,Ysig, s, xx, qq, fbest)) in enumerate(zip(
                names, (allcxx, allcyy), fitvals)):

            if plotfn is not None:
                plt.subplot(3,1, 1+i)
            else:
                plt.clf()

            if i in [0,1]:
                seeings.append(s[0])
                plt.errorbar(X, pixvar2seeing(Ymn), yerr=pixvar2seeing(Ysig),
                             fmt='o', color='g')
                plt.plot(allfocus, pixvar2seeing(Y), 'b.', alpha=0.25)
            ax = plt.axis()
            if i in [0,1]:
                plt.plot(xx, pixvar2seeing(qq), 'b-', alpha=0.5)
            plt.axis(ax)
            plt.ylim(0.5, 2.5)
            plt.xlabel('Focus shift (um)')
            plt.ylabel(name)
            if plotfn is not None:
                plt.text(ax[0] + 0.9*(ax[1]-ax[0]), 20,
                         'Focus: %.0f' % fbest,
                         ha='right', va='top', fontsize=12,
                         bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
            else:
                plt.title('Focus: %.1f' % fbest)
            plt.axvline(fbest, color='b')
            if plotfn is None:
                ps.savefig()

        if plotfn is not None:
            plt.subplot(3,1, 3)
        else:
            plt.clf()
        
        Y = allcxy
        plt.plot(allfocus, Y, 'b.', alpha=0.25)
        plt.ylim(-2,2)
        plt.xlabel('Focus shift (um)')
        plt.ylabel('PSF CXY')
        plt.axhline(0, color='k', alpha=0.25)
        for f in FF:
            plt.axvline(f, color='b', alpha=0.5)
        plt.axvline(fmean, color='b')

        # In pixel variance...
        meanseeing = np.mean(seeings)
        meanseeing = pixvar2seeing(meanseeing)
        
        if plotfn is not None:
            plt.subplot(3,1,1)
            plt.xlabel('')
            plt.subplot(3,1,2)
            plt.xlabel('')
            plt.suptitle('Focus (expnum = %i): Focus %.0f, Seeing %.1f' %
                         (self.primhdr.get('EXPNUM', 0), fmean, meanseeing))
            plt.savefig(plotfn)
        else:
            ps.savefig()
        
        return meas
            

    def fit_general_gaussian(self, img, sig1, xi, yi, fluxi, psf_r=15, ps=None):
        import tractor
        H,W = img.shape
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

        psf = tractor.NCircularGaussianPSF([4.], [1.])
        tim = tractor.Image(data=pix, inverr=ie, psf=psf)
        src = tractor.PointSource(tractor.PixPos(xi-xlo, yi-ylo),
                                  tractor.Flux(fluxi))
        tr = tractor.Tractor([tim],[src])

        src.pos.addGaussianPrior('x', 0., 1.)
    
        doplot = (ps is not None)
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
            if dlnp == 0:
                break

        # Now fit only the PSF size
        tr.freezeParam('catalog')
        # print('Optimizing params:')
        # tr.printThawedParams()

        for step in range(50):
            dlnp,x,alpha = tr.optimize(**optargs)
            if dlnp == 0:
                break

        # fwhms.append(psf.sigmas[0] * 2.35 * self.pixscale)
        
        if doplot:
            mod1 = tr.getModelImage(0)
            chi1 = tr.getChiImage(0)

        # Now switch to a non-isotropic PSF
        s = psf.sigmas[0]
        #print('Isotropic fit sigma', s)
        s = np.clip(s, 1., 5.)
        tim.psf = tractor.GaussianMixturePSF(1., 0., 0., s**2, s**2, 0.)
        
        #print('Optimizing params:')
        #tr.printThawedParams()

        try:
            for step in range(50):
                dlnp,x,alpha = tr.optimize(**optargs)
                #print('PSF:', tim.psf)
                if dlnp == 0:
                    break
        except:
            import traceback
            print('Error during fitting PSF in a focus frame; not to worry')
            traceback.print_exc()
            print('(The above was just an error during fitting one star in a focus frame; not to worry.)')
                
        # Don't need to re-fit source params because PSF ampl and mean
        # can fit for flux and position.

        if doplot:
            mod2 = tr.getModelImage(0)
            chi2 = tr.getChiImage(0)
            kwa = dict(vmin=-3*sig1, vmax=50*sig1, cmap='gray')

            plt.clf()
            plt.subplot(2,3,1)
            plt.title('Image')
            dimshow(pix, ticks=False, **kwa)
            plt.subplot(2,3,2)
            plt.title('Initial model')
            dimshow(mod0, ticks=False, **kwa)
            plt.subplot(2,3,3)
            plt.title('Isotropic model')
            dimshow(mod1, ticks=False, **kwa)
            plt.subplot(2,3,4)
            plt.title('Final model')
            dimshow(mod2, ticks=False, **kwa)
            plt.subplot(2,3,5)
            plt.title('Isotropic chi')
            dimshow(chi1, vmin=-10, vmax=10, ticks=False)
            plt.subplot(2,3,6)
            plt.title('Final chi')
            dimshow(chi2, vmin=-10, vmax=10, ticks=False)
            plt.suptitle('PSF fit')
            ps.savefig()
    
        return tim.psf.getParams()[-3:]
        
        
if __name__ == '__main__':
    import sys
    from astrometry.util.file import *
    import optparse

    parser = optparse.OptionParser(usage='%prog <mosaic image file.fits, eg mos3_62914.fits.gz>')
    parser.add_option('--ext', help='Extension to read: default %default',
                      default="im4")
    parser.add_option('--plot-prefix', default='focus',
                      help='Filename prefix for plots')
    opt,args = parser.parse_args()

    if len(args) == 0:
        parser.print_help()
        sys.exit(-1)
    
    ps = PlotSequence(opt.plot_prefix)

    from camera_mosaic import nominal_cal
    nom = nominal_cal

    for fn in args:
        meas = Mosaic3FocusMeas(fn, opt.ext, nom)
        meas.run(ps)
        
