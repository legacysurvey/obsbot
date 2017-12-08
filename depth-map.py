from __future__ import print_function
from astrometry.blind.plotstuff import *
from astrometry.util.util import *
from astrometry.util.fits import *
import pylab as plt
import matplotlib
from collections import Counter

def mosaic_wcs(ra, dec, pixbin=1.):
    # This is pretty close to the outline of the four Mosaic chips.
    W = H = (4096 * 2 + 100) / pixbin
    cd = pixbin * 0.262 / 3600.
    tan = Tan(ra, dec, W/2., H/2., cd, 0., 0., cd,
              float(W), float(H))
    return tan

# def plot_exposure(plot, ra, dec):
#     wcs = mosaic_wcs(ra, dec)
#     plot.outline.wcs = anwcs_new_tan(wcs)
#     plot.plot('outline')

def from_ccds():
    from legacypipe.survey import LegacySurveyData

    W,H = 1300,800
    plot = Plotstuff(size=(W,H), outformat='png')
    zoom = 2.5
    plot.wcs = anwcs_create_hammer_aitoff(195., 60., zoom, W, H, True)
    
    survey = LegacySurveyData()
    T = survey.get_annotated_ccds()

    print('Tile passes from CCDs table:')
    print(Counter(T.tilepass).most_common())
    
    plt.figure(2)
    
    plt.figure(1, figsize=(13,8))
    plt.subplots_adjust(left=0.01, right=0.99, bottom=0.01, top=0.95)

    for passnum in [1,2,3, 0]:
        T.depth = T.galdepth - T.decam_extinction[:,4]
        if passnum == 0:
            I = np.flatnonzero((T.tilepass >= 1) *
                               (T.tilepass <= 3) *
                               (T.depth > 10))
            tt = 'All Passes'
        else:
            I = np.flatnonzero((T.tilepass == passnum) *
                               (T.depth > 10))
            tt = 'Pass %i' % passnum
    
        plot.color = 'black'
        plot.plot('fill')
        plot.color = 'white'
        plot.outline.fill = True
        plot.outline.stepsize = 2000
    
        targetdepth = 22.5
        maxfrac = 0.4
    
        targetsig = 10.**((targetdepth - 22.5) / -2.5)
        targetiv = 1./targetsig**2
    
        print(len(I), 'tiles for pass', passnum)
        plot.op = CAIRO_OPERATOR_ADD
        for j,i in enumerate(I):
            if j % 1000 == 0:
                print('Tile', j)
            depth = T.depth[i]
            detsig = 10.**((depth - 22.5) / -2.5)
            depthfrac = 1./detsig**2 / targetiv
    
            mwcs = survey.get_approx_wcs(T[i])
            plot.outline.wcs = anwcs_new_tan(mwcs)
            plot.alpha = np.clip(maxfrac * depthfrac, 0., 1.)
            plot.apply_settings()
            plot.plot('outline')

        img = plot.get_image_as_numpy(flip=True)
        print('Image', img.dtype, img.shape)
        img = img[:,:,0]
        # It's a uint8
        img = img * 1./255.
        print('Image', img.dtype, img.shape)
        img /= maxfrac
        # Now it's in factor of detiv of target depth.
        # back to sig -> depth
        img = -2.5 * (np.log10(1./np.sqrt(img * targetiv)) - 9.)
        img[np.logical_not(np.isfinite(img))] = 0.

        cm = matplotlib.cm.viridis
        #cm = matplotlib.cm.jet
        cm = cmap_discretize(cm, 10)
        
        plt.figure(2)

        plt.clf()
        lo,hi = 22,23
        hh = img.ravel()
        hh = hh[hh != 0]
        plt.hist(np.clip(hh, lo, hi), 50, range=(lo,hi))
        plt.xlim(lo, hi)
        plt.title('Depths: %s' % tt)
        plt.yticks([])
        plt.savefig('depth-p%i-10.png' % passnum)

        plt.clf()
        plt.hist(np.clip(T.depth[I], lo, hi), 50, range=(lo,hi))
        plt.xlim(lo, hi)
        plt.title('Depths: %s' % tt)
        plt.savefig('depth-p%i-11.png' % passnum)

        plt.figure(1)
        plt.clf()
        img[img == 0] = np.nan
        plt.imshow(img, interpolation='nearest', origin='lower',
                   vmin=lo, vmax=hi, cmap=cm)
        plt.xticks([]); plt.yticks([])
        ax = plt.axis()
        H,W = img.shape
        for d in range(30, 90, 10):
            rr = np.arange(0, 360, 1)
            dd = np.zeros_like(rr) + d
            ok,xx,yy = plot.wcs.radec2pixelxy(rr, dd)
            plt.plot(xx, H-yy, 'k-', alpha=0.1)
        for r in range(0, 360, 30):
            dd = np.arange(20, 90, 1.)
            rr = np.zeros_like(dd) + r
            ok,xx,yy = plot.wcs.radec2pixelxy(rr, dd)
            plt.plot(xx, H-yy, 'k-', alpha=0.1)
        plt.axis(ax)
        plt.colorbar()
        plt.title(tt)
        plt.savefig('depth-p%i-12.png' % passnum)

        plot.op = CAIRO_OPERATOR_OVER
    
def main():
    W,H = 1300,800
    plot = Plotstuff(size=(W,H), outformat='png')
    zoom = 2.5
    plot.wcs = anwcs_create_hammer_aitoff(195., 60., zoom, W, H, True)
    #plot.wcs = anwcs_create_mollweide(195., 60., zoom, W, H, True)
    #print('WCS:', anwcs_print_stdout(plot.wcs))
    
    T = fits_table('mosaic-obstatus-depth.fits')
    
    plt.figure(2)
    
    plt.figure(1, figsize=(13,8))
    plt.subplots_adjust(left=0.01, right=0.99, bottom=0.01, top=0.95)
    
    for passnum in [1,2,3, 0]:
        if passnum == 0:
            I = np.flatnonzero((T.get('pass') >= 1) *
                               (T.get('pass') <= 3) *
                               (T.z_done == 1) *
                               (T.z_depth > 1) * (T.z_depth < 30))
            tt = 'MzLS z-band depth (AB, 5-sigma galaxy profile)' #'All Passes'
        else:
            I = np.flatnonzero((T.get('pass') == passnum) *
                               (T.z_done == 1) *
                               (T.z_depth > 1) * (T.z_depth < 30))
            tt = 'Pass %i' % passnum
    
        plot.color = 'black'
        plot.plot('fill')
        plot.color = 'white'
        plot.outline.fill = True
        plot.outline.stepsize = 2000
    
        targetdepth = 22.5
        maxfrac = 0.4
    
        targetsig = 10.**((targetdepth - 22.5) / -2.5)
        targetiv = 1./targetsig**2
    
        print(len(I), 'tiles for pass', passnum)
        plot.op = CAIRO_OPERATOR_ADD
        for j,i in enumerate(I):
            if j % 1000 == 0:
                print('Tile', j)
            depth = T.z_depth[i]
            detsig = 10.**((depth - 22.5) / -2.5)
            depthfrac = 1./detsig**2 / targetiv
    
            #print('depth', depth, 'frac', depthfrac)
    
            mwcs = mosaic_wcs(T.ra[i], T.dec[i])
            plot.outline.wcs = anwcs_new_tan(mwcs)
            plot.alpha = np.clip(maxfrac * depthfrac, 0., 1.)
            plot.apply_settings()
            plot.plot('outline')
    
        img = plot.get_image_as_numpy(flip=True)
        print('Image', img.dtype, img.shape)
        img = img[:,:,0]
        # It's a uint8
        img = img * 1./255.
        print('Image', img.dtype, img.shape)
        # 
        img /= maxfrac
    
        # plt.clf()
        # plt.hist(img.ravel(), 120, range=(0.1,2))
        # plt.savefig('depth-p%i-4.png' % passnum)
    
        # Now it's in factor of detiv of target depth.
        # back to sig -> depth
        img = -2.5 * (np.log10(1./np.sqrt(img * targetiv)) - 9.)
        img[np.logical_not(np.isfinite(img))] = 0.

        cm = matplotlib.cm.viridis
        #cm = matplotlib.cm.jet
        cm = cmap_discretize(cm, 10)
        
        plt.figure(2)
        plt.clf()
        #lo,hi = 21.5,23.5
        lo,hi = 22,23
        hh = img.ravel()
        hh = hh[hh != 0]
        plt.hist(np.clip(hh, lo, hi), 50, range=(lo,hi))
        plt.xlim(lo, hi)
        plt.title('Depths: %s' % tt)
        plt.yticks([])
        plt.savefig('depth-p%i-4.png' % passnum)

        plt.clf()
        plt.hist(np.clip(T.z_depth[I], lo, hi), 50, range=(lo,hi))
        plt.xlim(lo, hi)
        plt.title('Depths: %s' % tt)
        plt.savefig('depth-p%i-5.png' % passnum)
    
        plt.figure(1)
        plt.clf()
        img[img == 0] = np.nan
        plt.imshow(img, interpolation='nearest', origin='lower',
                   vmin=lo, vmax=hi, cmap=cm)
        plt.xticks([]); plt.yticks([])
        ax = plt.axis()
        H,W = img.shape
        for d in range(30, 90, 10):
            rr = np.arange(0, 360, 1)
            dd = np.zeros_like(rr) + d
            ok,xx,yy = plot.wcs.radec2pixelxy(rr, dd)
            plt.plot(xx, H-yy, 'k-', alpha=0.1)
            ok,xx,yy = plot.wcs.radec2pixelxy(85, d)
            plt.text(xx, H-yy, '%+i' % d,
                     bbox=dict(edgecolor='none', facecolor='white'),
                     ha='center', va='center')
        for r in range(0, 360, 30):
            dd = np.arange(20, 90, 1.)
            rr = np.zeros_like(dd) + r
            ok,xx,yy = plot.wcs.radec2pixelxy(rr, dd)
            plt.plot(xx, H-yy, 'k-', alpha=0.1)
            ok,xx,yy = plot.wcs.radec2pixelxy(r, 28)
            plt.text(xx, H-yy, '%i' % r,
                     bbox=dict(edgecolor='none', facecolor='white'),
                     ha='center', va='center')
            
        plt.axis(ax)
        plt.colorbar()
        plt.title(tt)
        plt.savefig('depth-p%i-3.png' % passnum)

        #if passnum == 3:
        if True:
            R = fits_table('retirable-p3.fits')
            tileids = set(R.tileid)
            I = np.array([i for i,tid in enumerate(T.tileid) if tid in tileids])
            imH,imW = img.shape
            for i in I:
                mwcs = mosaic_wcs(T.ra[i], T.dec[i])
                H,W = mwcs.shape
                rr,dd = mwcs.pixelxy2radec(np.array([1,1,W,W,1]),
                                           np.array([1,H,H,1,1]))
                ok,xx,yy = plot.wcs.radec2pixelxy(rr, dd)
                plt.plot(xx, imH - yy, 'r-')
                
        plt.axis(ax)
        plt.title('Retirable: %s' % tt)
        plt.savefig('depth-p%i-6.png' % passnum)

        
        plot.op = CAIRO_OPERATOR_OVER
    
        # plot.write('depth-p%i-1.png' % passnum)
        # 
        # plot.color = 'gray'
        # plot.plot_grid(30, 10, 30, 10)
        # plot.write('depth-p%i-2.png' % passnum)
        # 
        # rgb = cm(np.clip((img - lo) / (hi - lo), 0., 1.), bytes=True)
        # #rgb = (rgb * 255.0).astype(np.uint8)
        # I,J = np.nonzero(np.logical_not(np.isfinite(img)))
        # print(len(I), 'inf pixels')
        # #for i in range(3):
        # #    rgb[I,J,i] = 0
        # rgb[I,J,:3] = 255
        # rgb = np.flipud(rgb)
        # 
        # print('RGB:', rgb.dtype, rgb.shape, rgb.min(), rgb.max())
        # plot.set_image_from_numpy(rgb)
        # plot.color = 'black'
        # plot.bg_rgba = (0,0,0,0)
        # plot.plot_grid(30, 10, 30, 10)
        # plot.write('depth-p%i-6.png' % passnum)
        
        

# From http://scipy-cookbook.readthedocs.io/items/Matplotlib_ColormapTransformations.html
def cmap_discretize(cmap, N):
    from matplotlib.cm import get_cmap
    from numpy import concatenate, linspace
    """Return a discrete colormap from the continuous colormap cmap.

        cmap: colormap instance, eg. cm.jet.
        N: number of colors.

    Example
        x = resize(arange(100), (5,100))
        djet = cmap_discretize(cm.jet, 5)
        imshow(x, cmap=djet)
    """
    if type(cmap) == str:
        cmap = get_cmap(cmap)
    colors_i = concatenate((linspace(0, 1., N), (0.,0.,0.,0.)))
    colors_rgba = cmap(colors_i)
    indices = linspace(0, 1., N+1)
    cdict = {}
    for ki,key in enumerate(('red','green','blue')):
        cdict[key] = [ (indices[i], colors_rgba[i-1,ki], colors_rgba[i,ki]) for i in xrange(N+1) ]
    # Return colormap object.
    return matplotlib.colors.LinearSegmentedColormap(cmap.name + "_%d"%N, cdict, 1024)



if __name__ == '__main__':
    from_ccds()
    main()
