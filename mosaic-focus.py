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
        #print('Focus shifts:', shifts)
        # if ps is not None:
        #     mn,mx = np.percentile(img.ravel(), [25,98])
        #     kwa = dict(vmin=mn, vmax=mx)
        #     plt.clf()
        #     dimshow(img, **kwa)
        #     ps.savefig()
            
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

            # plt.clf()
            # plt.hist(minimg.ravel(), range=(-10,30), bins=100)
            # plt.title('Min image')
            # ps.savefig()

            
        detsn = minimg
        # zero out the edges -- larger margin here?
        detsn[0 ,:] = 0
        detsn[:, 0] = 0
        detsn[-1,:] = 0
        detsn[:,-1] = 0
        return detsn

    def trim_edges(self, img):
        # Trim off some edge pixels.
        # In the focus images, the bottom ~300 pixels are funky.
        trim = self.edge_trim
        bottom = 300
        cimg = img[bottom:-trim, trim:-trim]
        return cimg, trim, bottom

    def run(self, ps=None):
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
        
        (mx2 ,my2 ,mxy ,mtheta,ma,mb,mell) = meas['moments' ]
        (wmx2,wmy2,wmxy,wtheta,wa,wb,well) = meas['wmoments']
        
        radius2 = 3. / self.pixscale
        I,J,d = match_xy(px, py, fx, fy, radius2)

        ps1band = ps1cat.ps1band[band]
        ps1mag = stars.median[I, ps1band]
        
        # if ps is not None:
        #     kwa = dict(vmin=-3*sig1, vmax=50*sig1, cmap='gray')
        #     plt.clf()
        #     dimshow(img, **kwa)
        #     ax = plt.axis()
        #     plt.plot(px, py, 'o', mec='m', mfc='none', ms=6)
        #     plt.axis(ax)
        #     plt.title('All PS1 stars')
        #     plt.colorbar()
        #     ps.savefig()

        # Choose the brightest PS1 stars not near an edge
        hdr = self.hdr
        shifts = self.get_focus_shifts(hdr)
        S = 20
        H,W = img.shape
        miny = fy[J] + min(shifts)
        maxy = fy[J] + max(shifts)
        edge = reduce(np.logical_or, [fx[J] < S, fx[J]+S >= W,
                                      miny  < S, maxy +S >= H])
        K = np.argsort(ps1mag + 1000*edge)

        nstars = 10
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

        # We no longer detect & measure each of the shifted sources,
        # so this doesn't work.
        # 
        # xl = min(shifts) - 0.5 * np.abs(steppix)
        # xh = max(shifts) + 0.5 * np.abs(steppix)
        # Jx = [np.random.normal(size=len(Ji)) *(0.1 * np.abs(steppix))
        # for Ji in JJ]
        # for YY,TT in [([mx2, my2, mxy], ['mx2', 'my2', 'mxy']),
        #               ([wmx2,wmy2,wmxy],['wx2', 'wy2', 'wxy']),
        #               ([ma,mb, mell,mtheta],['ma', 'mb', 'mell', 'mtheta']),
        #               ([wa,wb, well,wtheta],['wa', 'wb', 'well', 'wtheta']),
        #               ]:
        #     plt.clf()
        #     for i,(Y,name) in enumerate(zip(YY,TT)):
        #         plt.subplot(len(YY),1,1+i)
        #         for Ji,jx,shifty in zip(JJ, Jx, shifts):
        #             plt.plot(shifty+jx, Y[Ji], 'b.', alpha=0.25)
        #             # HACK
        #             if not 'theta' in name:
        #                 plt.plot(shifty, np.median(Y[Ji]), 'ro')
        #         plt.ylim(*np.percentile(np.hstack([Y[Ji] for Ji in JJ]), [5,95]))
        #         plt.ylabel(name)
        #         plt.xlim(xl,xh)
        #         plt.suptitle('Focus sweep measurements')
        #     ps.savefig()

        meas.update(I=I, J=J, K=K)

        shifts = self.get_focus_shifts(hdr)
        focus = self.get_focus_values(hdr)
    
        allfocus = []
        #allshifts = []
        allcxx,allcyy,allcxy = [],[],[]
        #nstars = 40
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
                #allshifts.append(shifty)
                allfocus.append(foc)
                allcxx.append(p[0])
                allcyy.append(p[1])
                allcxy.append(p[2])
    
        allfocus  = np.array(allfocus)
        #allshifts = np.array(allshifts)
        allcxx = np.array(allcxx)
        allcyy = np.array(allcyy)
        allcxy = np.array(allcxy)

        meas.update(shifts=shifts, focus=focus,
                    allfocus=allfocus, allcxx=allcxx, allcyy=allcyy,
                    allcxy=allcxy)

        FF = []

        for Y,name in [(allcxx, 'PSF CXX'), (allcyy, 'PSF CYY')]:
            plt.clf()
    
            X,I = np.unique(allfocus, return_inverse=True)
            Ymn,Ysig = np.zeros(len(X)), np.zeros(len(X))
            for i in range(len(X)):
                J = np.flatnonzero(I == i)
                yi = Y[J]
                #y1,ymn,y2 = np.percentile(yi, [25,50,75])
                # 2 * 0.6745 = inter-quartile range for a Gaussian
                #ysig = (y1 - y2) / (2. * 0.6745)
                y1,ymn,y2 = np.percentile(yi, [16,50,84])
                ysig = (y1 - y2) / 2
                Ymn [i] = ymn
                Ysig[i] = ysig
            plt.plot(allfocus, Y, 'b.', alpha=0.25)
            plt.errorbar(X, Ymn, yerr=Ysig, fmt='o', color='g')
    
            xmn = np.median(X)
            dx = (np.max(X) - np.min(X)) / 2.
            # Rescale xx in [-1,1]
            xx = (X - xmn) / dx
            #print('xx', xx)
            A = np.zeros((len(xx),3))
            A[:,0] = 1.
            A[:,1] = xx
            A[:,2] = xx**2
            wt = 1. / Ysig
            A *= wt[:,np.newaxis]
            b = Ymn * wt
            R = np.linalg.lstsq(A, b)
            s = R[0]
            print('Least-squares quadratic:', s)
    
            ax = plt.axis()
            xx = np.linspace(ax[0], ax[1], 500)
            rx = (xx - xmn)/dx
            plt.plot(xx, s[0] + rx*s[1] + rx**2*s[2], 'b-', alpha=0.5)
            plt.axis(ax)
            plt.ylim(0, 16)
            plt.xlabel('Focus shift (um)')
            plt.ylabel(name)
            fbest = -s[1] / (2. * s[2])
            fbest = fbest * dx + xmn
            print('Focus position from %s: %.1f' % (name, fbest))
            plt.title('Focus: %.1f' % fbest)
            plt.axvline(fbest, color='b')
            ps.savefig()
            FF.append(fbest)

        fmean = np.mean(FF)
        print('Mean focus: %.1f' % fmean)
        
        plt.clf()
        Y = allcxy
        #yl,yh = np.percentile(Y, [5, 95])
        plt.plot(allfocus, Y, 'b.', alpha=0.25)
        #plt.ylim(yl,yh)
        plt.ylim(-2,2)
        plt.xlabel('Focus shift (um)')
        plt.ylabel('PSF CXY')
        plt.axhline(0, color='k', alpha=0.25)
        for f in FF:
            plt.axvline(f, color='b', alpha=0.5)
        plt.axvline(fmean, color='b')
        ps.savefig()
        
        return meas
            

    def fit_general_gaussian(self, img, sig1, xi, yi, fluxi, psf_r=15, ps=None):
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
            traceback.print_exc()
            pass
                
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

    for fn in args:
        meas = Mosaic3FocusMeas(fn, opt.ext)
        meas.run(ps)
        
