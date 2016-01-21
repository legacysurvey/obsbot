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
    
    # def match_ps1_stars(self, px, py, fullx, fully, radius, stars):
    #     # check hdr for focus sequence offsets
    #     shifts = self.get_focus_shifts(self.hdr)
    #     print('Focus shifts:', shifts)
    # 
    #     # Duplicate "px","py" for each offset.
    #     # We *subtract* shift because we want to infer the un-shifted location
    #     allpx = np.hstack([px] * len(shifts))
    #     allpy = np.hstack([py - s for s in shifts])
    #     
    #     I,J,d = match_xy(allpx, allpy, fullx, fully, radius)
    #     dx = allpx[I] - fullx[J]
    #     dy = allpy[I] - fully[J]
    #     # de-duplicate
    #     I = I % len(px)
    #     return I,J,dx,dy

    def detection_map(self, img, sig1, psfsig, ps):
        # Shift-and-add the image
        shifts = self.get_focus_shifts(self.hdr)
        print('Focus shifts:', shifts)

        if ps is not None:
            mn,mx = np.percentile(img.ravel(), [25,98])
            kwa = dict(vmin=mn, vmax=mx)
            plt.clf()
            dimshow(img, **kwa)
            ps.savefig()
            
        # shadd = np.zeros_like(img)
        # for dy in -shifts:
        #     if dy == 0:
        #         shadd += img
        #     elif dy > 0:
        #         shadd[dy:,:] += img[:-dy,:]
        #     elif dy < 0:
        #         shadd[:dy,:] += img[-dy:,:]
        # 
        # if ps is not None:
        #     plt.clf()
        #     dimshow(shadd, **kwa)
        #     plt.title('Shift-n-add')
        #     ps.savefig()
        #         
        # psfnorm = 1./(2. * np.sqrt(np.pi) * psfsig)
        # detsn = gaussian_filter(shadd / sig1, psfsig) / psfnorm

        psfnorm = 1./(2. * np.sqrt(np.pi) * psfsig)
        detsn = gaussian_filter(img / sig1, psfsig) / psfnorm
        
        # Take the *minimum* of the image and shifted versions of itself
        minimg = detsn.copy()
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
            plt.clf()
            #dimshow(shadd, **kwa)
            #plt.title('Shift-n-add')
            dimshow(minimg, **kwa)
            plt.title('Min image')
            ps.savefig()

        detsn = minimg
            
        # zero out the edges -- larger margin here?
        detsn[0 ,:] = 0
        detsn[:, 0] = 0
        detsn[-1,:] = 0
        detsn[:,-1] = 0
        return detsn

    # def detect_sources(self, detsn, thresh, ps):
    #     slices = super(Mosaic3FocusMeas, self).detect_sources(detsn, thresh, ps)
    # 
    #     #shifts = self.get_focus_shifts(self.hdr)
    # 
    #     # Order by peak flux.  For each source, suppress other sources at += shifts
    #     pky,pkx,flux = [],[],[]
    #     for slc in slices:
    #         y0 = slc[0].start
    #         x0 = slc[1].start
    #         subimg = detsn[slc]
    #         imax = np.argmax(subimg)
    #         y,x = np.unravel_index(imax, subimg.shape)
    #         pky.append(y+y0)
    #         pkx.append(x+x0)
    #         flux.append(subimg.flat[imax])
    #     I = np.argsort(-np.array(flux))
    #     pkx = np.array(pkx)
    #     pky = np.array(pky)
    # 
    #     suppressed = np.zeros(detsn.shape, bool)
    #     keepi = []
    #     S = 2
    # 
    #     hdr = self.hdr
    #     nsteps = hdr['FOCNEXPO']
    #     steppix = np.abs(hdr['FOCSHIFT'])
    # 
    #     H,W = detsn.shape
    #     for i,x,y in zip(I, pkx[I], pky[I]):
    #         if suppressed[y,x]:
    #             continue
    #         keepi.append(i)
    #         for j in range(-nsteps, nsteps+1):
    #             ys = y + steppix * j
    #             suppressed[np.clip(ys-S,0,H) : np.clip(ys+S+1,0,H),
    #                        max(0,x-S) : min(W,x+S+1)] = True
    #     return [slices[i] for i in keepi]
    
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

        (mx2 ,my2 ,mxy ,mtheta,ma,mb,mell) = meas['moments' ]
        (wmx2,wmy2,wmxy,wtheta,wa,wb,well) = meas['wmoments']
        
        radius2 = 3. / self.pixscale
        # Mosaic focus: un-duplicated
        I,J,d = match_xy(px, py, fx, fy, radius2)
        print('De-duplicated:', len(I), 'matches')
        
        if ps is not None:
            kwa = dict(vmin=-3*sig1, vmax=50*sig1, cmap='gray')
            plt.clf()
            dimshow(img, **kwa)
            ax = plt.axis()
            plt.plot(fx[J], fy[J], 'go', mec='g', mfc='none', ms=10)
            plt.plot(px[I], py[I], 'm+', ms=10)
            plt.axis(ax)
            plt.title('Matched PS1 stars')
            plt.colorbar()
            ps.savefig()

            plt.clf()
            dimshow(img, **kwa)
            ax = plt.axis()
            plt.plot(px, py, 'o', mec='m', mfc='none', ms=6)
            plt.axis(ax)
            plt.title('All PS1 stars')
            plt.colorbar()
            ps.savefig()

        hdr = self.hdr
        nsteps = hdr['FOCNEXPO']
        steppix = hdr['FOCSHIFT']
        shifts = np.append(0, (2+np.arange(nsteps-1))*steppix)

        ps1band = ps1cat.ps1band[band]
        ps1mag = stars.median[I, ps1band]

        S = 20
        H,W = img.shape

        miny = fy[J] - max(shifts)
        maxy = fy[J] - min(shifts)
        edge = reduce(np.logical_or, [fx[J] < S, fx[J]+S >= W,
                                      miny  < S, maxy +S >= H])
        K = np.argsort(ps1mag + 1000*edge)
        
        nstars = 10
        plt.clf()
        plt.subplots_adjust(hspace=0, wspace=0)
        sp = 1
        S = 20
        for shifty in shifts:
            for k in K[:nstars]:
                plt.subplot(len(shifts), nstars, sp)
                sp += 1
                i = I[k]
                j = J[k]
                x = fx[j]
                y = fy[j] - shifty
                dimshow(img[y-S:y+S+1, x-S:x+S+1], ticks=False, **kwa)
                if shifty == 0.:
                    print('Plotting focus sweep star at', x,y)

        plt.suptitle('Focus sweep stars')
        ps.savefig()
            
        JJ = []
        for shifty in shifts:
            print('shifty:', shifty)
            I2,J2,d = match_xy(px, py - shifty, fx, fy, radius2)
            print('Matched', len(I2), 'PS1 stars')
            JJ.append(J2)

        xl = min(shifts) - 0.5 * np.abs(steppix)
        xh = max(shifts) + 0.5 * np.abs(steppix)
        Jx = [np.random.normal(size=len(Ji)) *(0.1 * np.abs(steppix)) for Ji in JJ]

        for YY,TT in [([mx2, my2, mxy], ['mx2', 'my2', 'mxy']),
                      ([wmx2,wmy2,wmxy],['wx2', 'wy2', 'wxy']),
                      ([ma,mb, mell,mtheta],['ma', 'mb', 'mell', 'mtheta']),
                      ([wa,wb, well,wtheta],['wa', 'wb', 'well', 'wtheta']),
                      ]:
            plt.clf()
            for i,(Y,name) in enumerate(zip(YY,TT)):
                plt.subplot(len(YY),1,1+i)
                for Ji,jx,shifty in zip(JJ, Jx, shifts):
                    plt.plot(shifty+jx, Y[Ji], 'b.', alpha=0.25)
                    # HACK
                    if not 'theta' in name:
                        plt.plot(shifty, np.median(Y[Ji]), 'ro')
                plt.ylim(*np.percentile(np.hstack([Y[Ji] for Ji in JJ]), [5,95]))
                plt.ylabel(name)
                plt.xlim(xl,xh)
                plt.suptitle('Focus sweep measurements')
            ps.savefig()

        meas.update(I=I, J=J, K=K, shifts=shifts)

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
        print('Isotropic fit sigma', s)
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
                
        # Now also fit the source params
        # tr.thawParam('catalog')
        # print('Optimizing params:')
        # tr.printThawedParams()
        # 
        # try:
        #     for step in range(50):
        #         dlnp,x,alpha = tr.optimize(**optargs)
        #         if dlnp == 0:
        #             break
        # except:
        #     import traceback
        #     traceback.print_exc()
        #     pass

        #print('Fit PSF:', tim.psf)

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
    from astrometry.util.file import *
    
    ps = PlotSequence('focus')

    meas = Mosaic3FocusMeas('mos3_62914.fits.gz', 'im4')

    R = meas.run(ps=ps)
    import sys
    sys.exit(0)
    
    pfn = 'meas.pickle'
    if os.path.exists(pfn):
        print('Reading', pfn)
        R = unpickle_from_file(pfn)
    else:
        R = meas.run(ps=ps)
        pickle_to_file(R, pfn)
        print('Wrote', pfn)

    K = R['K']
    J = R['J']
    fx = R['fx']
    fy = R['fy']
    shifts = R['shifts']
    img = R['img']
    sig1 = R['sig1']
    apflux = R['apflux']

    hdr = R['hdr']

    # these are value in microns
    focstart = hdr['FOCSTART']
    focstep  = hdr['FOCSTEP']

    # assume...
    focvals = focstart + np.arange(len(shifts)) * focstep
    
    allfocus = []
    allshifts = []
    allfocus = []
    allcxx,allcyy,allcxy = [],[],[]
    nstars = 40
    for ik,k in enumerate(K[:nstars]):
        j = J[k]
        xi = fx[j]
        yi = fy[j]
        fluxi = apflux[j]
        print('Fitting star at', xi,yi)
        for shifty,foc in zip(shifts, focvals):
            j = J[k]
            xi = fx[j]
            yi = fy[j] - shifty

            p = meas.fit_general_gaussian(img, sig1, xi, yi, fluxi, ps=ps if ik ==0 else None)

            print('PSF variance:', p)
            allshifts.append(shifty)
            allfocus.append(foc)
            allcxx.append(p[0])
            allcyy.append(p[1])
            allcxy.append(p[2])

    allfocus  = np.array(allfocus)
    allshifts = np.array(allshifts)
    allcxx = np.array(allcxx)
    allcyy = np.array(allcyy)
    allcxy = np.array(allcxy)

    plt.clf()
    Y = allcxx
    #yl,yh = np.percentile(Y, [5, 95])
    plt.plot(allfocus, Y, 'b.', alpha=0.25)
    #plt.ylim(yl,yh)
    plt.ylim(0, 16)
    plt.xlabel('Focus shift (um)')
    plt.ylabel('PSF CXX')
    ps.savefig()

    plt.clf()
    Y = allcyy
    #yl,yh = np.percentile(Y, [5, 95])
    plt.plot(allfocus, Y, 'b.', alpha=0.25)
    #plt.ylim(yl,yh)
    plt.ylim(0, 16)
    plt.xlabel('Focus shift (um)')
    plt.ylabel('PSF CYY')
    ps.savefig()

    plt.clf()
    Y = allcxy
    #yl,yh = np.percentile(Y, [5, 95])
    plt.plot(allfocus, Y, 'b.', alpha=0.25)
    #plt.ylim(yl,yh)
    plt.ylim(-2,2)
    plt.xlabel('Focus shift (um)')
    plt.ylabel('PSF CXY')
    ps.savefig()
    
