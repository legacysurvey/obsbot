from __future__ import print_function
if __name__ == '__main__':
    import matplotlib
    matplotlib.use('Agg')
from astrometry.util.plotutils import *
from measure_raw import *

import numpy as np
import pylab as plt

class Mosaic3FocusMeas(Mosaic3Measurer):
    def match_ps1_stars(self, px, py, fullx, fully, radius, stars):
        # check hdr for focus sequence offsets
        hdr = self.hdr
        nsteps = hdr['FOCNEXPO']
        steppix = hdr['FOCSHIFT']
        print('Focus:', nsteps, 'steps of', steppix, 'pixels')
        # Duplicate "px","py" for each offset?
        allpx = np.hstack([px] + [px                  for i in np.arange(nsteps-1)])
        allpy = np.hstack([py] + [py -(2 + i)*steppix for i in np.arange(nsteps-1)])
        I,J,d = match_xy(allpx, allpy, fullx, fully, radius)
        dx = allpx[I] - fullx[J]
        dy = allpy[I] - fully[J]
        # de-duplicate
        I = I % len(px)
        return I,J,dx,dy

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
        plt.suptitle('Focus sweep stars')
        ps.savefig()
            
        JJ = []
        for shifty in shifts:
            print('shifty:', shifty)
            I2,J2,d = match_xy(px, py - shifty, fx, fy, radius2)
            print('Matched', len(I2), 'PS1 stars')
            JJ.append(J2)

        xl,xh = min(shifts) - 0.5 * np.abs(steppix), max(shifts) + 0.5 * np.abs(steppix)
        Jx = [np.random.normal(size=len(J)) *(0.1 * np.abs(steppix)) for J in JJ]

        for YY,TT in [([mx2, my2, mxy], ['mx2', 'my2', 'mxy']),
                      ([wmx2,wmy2,wmxy],['wx2', 'wy2', 'wxy']),
                      ([ma,mb, mell,mtheta],['ma', 'mb', 'mell', 'mtheta']),
                      ([wa,wb, well,wtheta],['wa', 'wb', 'well', 'wtheta']),
                      ]:
            plt.clf()
            for i,(Y,name) in enumerate(zip(YY,TT)):
                plt.subplot(len(YY),1,1+i)
                for J,jx,shifty in zip(JJ, Jx, shifts):
                    plt.plot(shifty+jx, Y[J], 'b.', alpha=0.25)
                    # HACK
                    if not 'theta' in name:
                        plt.plot(shifty, np.median(Y[J]), 'ro')
                plt.ylim(*np.percentile(np.hstack([Y[J] for J in JJ]), [5,95]))
                plt.ylabel(name)
                plt.xlim(xl,xh)
                plt.suptitle('Focus sweep measurements')
            ps.savefig()

        return meas
        
        
if __name__ == '__main__':
    ps = PlotSequence('focus')
    meas = Mosaic3FocusMeas('mos3_62914.fits.gz', 'im4')
    meas.run(ps=ps)
        
        
