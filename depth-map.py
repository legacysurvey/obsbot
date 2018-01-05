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

def from_ccds():
    from legacypipe.survey import LegacySurveyData

    # W,H = 1300,800
    # plot = Plotstuff(size=(W,H), outformat='png')
    # zoom = 2.5
    # plot.wcs = anwcs_create_hammer_aitoff(195., 60., zoom, W, H, True)

    survey = LegacySurveyData()
    T = survey.get_annotated_ccds()
    print('Read', len(T), 'CCDs from annotated table')
    T.cut(T.ccd_cuts == 0)
    print('Cut to', len(T), 'with ccd_cuts == 0')
    print('Tile passes from CCDs table:')
    print(Counter(T.tilepass).most_common())

    T.depth = T.galdepth - T.decam_extinction[:,4]
    
    imgs = []
    for passnum in [1,2,3, 9]:#, 0]:
        if passnum == 0:
            I = np.flatnonzero((T.tilepass >= 1) *
                               (T.tilepass <= 3) *
                               (T.depth > 10))
            tt = 'All Passes'
        elif passnum == 9:
            I = []
            tt = 'All Passes'
        else:
            I = np.flatnonzero((T.tilepass == passnum) *
                               (T.depth > 10))
            tt = 'Pass %i' % passnum
        print(len(I), 'tiles for pass', passnum)

        depthmap = depth_map_for_ccds(survey, ccds[I])
        iv = 1./(10.**((depthmap - 22.5) / -2.5))**2

        if passnum in [1,2,3]:
            imgs.append(iv)
        if passnum == 9:
            img = imgs[0]
            for im in imgs[1:]:
                img += im

        fitsio.write('depth-p%i.fits' % passnum, img, header=hdr, clobber=True)
        print('Unique depths:', np.unique(img.ravel()))

        cm = matplotlib.cm.viridis
        #cm = matplotlib.cm.jet
        cm = cmap_discretize(cm, 10)
        
        plt.figure(2)

        plt.clf()
        lo,hi = 21.5,23.5
        hh = img.ravel()
        hh = hh[hh != 0]
        plt.hist(np.clip(hh, lo, hi), 50, range=(lo,hi))
        plt.xlim(lo, hi)
        plt.title('Depths: %s' % tt)
        plt.yticks([])
        plt.ylabel('sky area')
        plt.savefig('depth-p%i-10.png' % passnum)

        plt.clf()
        plt.hist(np.clip(T.depth[I], lo, hi), 50, range=(lo,hi))
        plt.xlim(lo, hi)
        plt.title('Depths: %s' % tt)
        plt.ylabel('Number of CCDs')
        plt.savefig('depth-p%i-11.png' % passnum)

        tiles = fits_table('obstatus/mosaic-tiles_obstatus.fits')
        #if passnum in [1,2,3]:
        #tiles.cut(tiles.get('pass') == passnum)
        tiles.cut(tiles.get('pass') == 1)
        tiles.cut(tiles.in_desi == 1)
        ok,tiles.x,tiles.y = plot.wcs.radec2pixelxy(tiles.ra, tiles.dec)
        tiles.x -= 1
        tiles.y -= 1
        tiles.cut(ok)
        
        plt.figure(1)
        plt.clf()
        img[img == 0] = np.nan

        plt.plot(tiles.x, H-tiles.y, 'k.', alpha=0.1)
        
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


        #if passnum == 0:
        if True:
            goal = 22.5
            drange = 1.
        
            plt.clf()
            img[img == 0] = np.nan
            plt.imshow(img, interpolation='nearest', origin='lower',
                       vmin=goal-drange, vmax=goal+drange, cmap='RdBu')
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
            cb = plt.colorbar()
            cb.set_label('Depth')
            plt.title(tt)
            plt.savefig('depth-p%i-13.png' % passnum)

def main():
    W,H = 1300,800
    plot = Plotstuff(size=(W,H), outformat='png')
    zoom = 2.5
    plot.wcs = anwcs_create_hammer_aitoff(195., 60., zoom, W, H, True)
    #plot.wcs = anwcs_create_mollweide(195., 60., zoom, W, H, True)
    #print('WCS:', anwcs_print_stdout(plot.wcs))

    # generated by update-obstatus.py
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
        #maxfrac = 0.4
        # This sets the maximum (saturation) value of the image map;
        # it's in invvar fluxes, so 0.16 means 1 mag deeper than "targetdepth"
        maxfrac = 0.16

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
        lo,hi = 21.5,23.5
        #lo,hi = 22,23
        hh = img.ravel()
        hh = hh[hh != 0]
        plt.hist(np.clip(hh, lo, hi), 50, range=(lo,hi))
        plt.xlim(lo, hi)
        plt.title('Depths: %s' % tt)
        plt.yticks([])
        plt.ylabel('Sky Area')
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


        goal = 22.5
        drange = 1.
        
        plt.clf()
        img[img == 0] = np.nan
        plt.imshow(img, interpolation='nearest', origin='lower',
                   vmin=goal-drange, vmax=goal+drange, cmap='RdBu')
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
        cb = plt.colorbar()
        cb.set_label('Depth')
        plt.title(tt)
        plt.savefig('depth-p%i-2.png' % passnum)


        
        # #if passnum == 3:
        # if True:
        #     R = fits_table('retirable-p3.fits')
        #     tileids = set(R.tileid)
        #     I = np.array([i for i,tid in enumerate(T.tileid) if tid in tileids])
        #     imH,imW = img.shape
        #     for i in I:
        #         mwcs = mosaic_wcs(T.ra[i], T.dec[i])
        #         H,W = mwcs.shape
        #         rr,dd = mwcs.pixelxy2radec(np.array([1,1,W,W,1]),
        #                                    np.array([1,H,H,1,1]))
        #         ok,xx,yy = plot.wcs.radec2pixelxy(rr, dd)
        #         plt.plot(xx, imH - yy, 'r-')
        #         
        # plt.axis(ax)
        # plt.title('Retirable: %s' % tt)
        # plt.savefig('depth-p%i-6.png' % passnum)

        
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

def tiles_todo():
    '''
    Makes maps like depth-p?.fits, but for tiles that are still marked to-do
    in the obstatus file.  We give them a nominal depth (22.2 is ~ the median)
    '''
    depth = 22.2

    plt.figure(1, figsize=(13,8))
    plt.subplots_adjust(left=0.01, right=0.99, bottom=0.01, top=0.95)
    
    tiles = fits_table('obstatus/mosaic-tiles_obstatus.fits')
    tiles.cut(tiles.get('pass') <= 3)
    tiles.cut(tiles.in_desi == 1)
    tiles.cut(tiles.dec > 30)
    print(len(tiles), 'above Dec=30 in passes', np.unique(tiles.get('pass')))

    #print('Z_DATEs:', np.unique(tiles.z_date))
    #for d in np.unique(tiles.z_date):
    #    print('date', d, 'after 2017-12-08?', d > '2017-12-08')

    tiles.recent = np.array([d > '2017-12-08' for d in tiles.z_date])
    
    print('To-do:', np.sum(tiles.z_done == 0))
    print('Done recently:', np.sum(tiles.recent))
    
    tiles.cut(np.logical_or(tiles.z_done == 0, tiles.recent))
    print(len(tiles), 'not done (or done recently)')
    print('Per pass:', Counter(tiles.get('pass')))

    ## Copied from from_ccds
    W,H = 4800,3200
    plot = Plotstuff(size=(W,H), outformat='png')
    zoom = 2.7
    plot.wcs = anwcs_create_hammer_aitoff(195., 65., zoom, W, H, True)
    
    imgs = []
    for passnum in [1,2,3, 9]:
        plot.color = 'black'
        plot.plot('fill')
        plot.color = 'white'
        plot.outline.fill = True
        plot.outline.stepsize = 2000
        plot.op = CAIRO_OPERATOR_ADD

        targetdepth = 22.5
        # This gives a 1-mag margin on the target depth
        maxfrac = 0.16

        targetsig = 10.**((targetdepth - 22.5) / -2.5)
        targetiv = 1./targetsig**2

        I = np.flatnonzero(tiles.get('pass') == passnum)
        for j,i in enumerate(I):
            detsig = 10.**((depth - 22.5) / -2.5)
            depthfrac = 1./detsig**2 / targetiv
            mwcs = mosaic_wcs(tiles.ra[i], tiles.dec[i])
            plot.outline.wcs = anwcs_new_tan(mwcs)
            plot.alpha = np.clip(maxfrac * depthfrac, 0., 1.)
            plot.apply_settings()
            plot.plot('outline')
        img = plot.get_image_as_numpy(flip=True)
        img = img[:,:,0]
        # It's a uint8; scale and convert to float
        img = img * 1./255.
        img /= maxfrac
        # Now it's in factor of detiv of target depth.
        imgs.append(img)

        if passnum == 9:
            img = imgs[0]
            for im in imgs[1:]:
                img += im
        
        # back to sig -> depth
        img = -2.5 * (np.log10(1./np.sqrt(img * targetiv)) - 9.)
        img[np.logical_not(np.isfinite(img))] = 0.

        fitsio.write('depth-todo-p%i.fits' % passnum, img, clobber=True)

        lo,hi = 21.5,23.5
        cm = matplotlib.cm.viridis
        cm = cmap_discretize(cm, 10)
        plt.clf()
        img[img == 0] = np.nan
        #plt.plot(tiles.x, H-tiles.y, 'k.', alpha=0.1)
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
        plt.title('To-do: pass %i' % passnum)
        plt.savefig('depth-p%i-14.png' % passnum)

        plot.op = CAIRO_OPERATOR_OVER
        
    
def needed_tiles():
    from astrometry.util.miscutils import point_in_poly

    targetdepth = 22.5

    fn = 'depth-p9.fits'
    #wcs = anwcs(fn)
    wcs = anwcs('plot.wcs')
    depth = fitsio.read(fn)
    print('Depth:', depth.shape, depth.dtype)

    ima = dict(interpolation='nearest', origin='lower',
               vmin=22.0, vmax=23.0, cmap='RdBu')
    plt.clf()
    plt.imshow(depth, **ima)
    plt.savefig('depth-0.png')

    todo = fitsio.read('depth-todo-p9.fits')

    iv1 = 1./(10.**((depth - 22.5) / -2.5))**2
    iv2 = 1./(10.**((todo  - 22.5) / -2.5))**2
    depth = -2.5 * (np.log10(1./np.sqrt(iv1 + iv2)) - 9.)
    print('inf:', np.sum(np.logical_not(np.isfinite(depth))))
    
    plt.clf()
    plt.imshow(depth, **ima)
    plt.savefig('depth-1.png')

    tiles = fits_table('obstatus/mosaic-tiles_obstatus.fits')
    tiles.cut(tiles.get('pass') == 1)
    tiles.cut(tiles.in_desi == 1)

    ### Drop ones that are already on the to-do list.
    tiles.cut(tiles.z_done == 1)
    
    tiles.cut(np.lexsort((tiles.ra, tiles.dec)))
    
    #ok,tiles.x,tiles.y = plot.wcs.radec2pixelxy(tiles.ra, tiles.dec)
    #tiles.x -= 1
    #tiles.y -= 1
    #tiles.cut(ok)
    print(len(tiles), 'tiles in pass 1')

    tiles.cut(tiles.dec > 30)
    print(len(tiles), 'with Dec > 30')

    print('Min Dec:', min(tiles.dec))

    #tiles.cut(tiles.dec < 70)
    #print(len(tiles), 'with Dec < 70')
    
    # from mosaic_wcs
    tilesize = (4096 * 2 + 100) * 0.262 / 3600.

    dlo = tiles.dec - tilesize/2.
    dhi = tiles.dec + tilesize/2.
    cosdec = np.cos(np.deg2rad(tiles.dec))
    rlo = tiles.ra - tilesize/2./cosdec
    rhi = tiles.ra + tilesize/2./cosdec
    ok1,tiles.x1,tiles.y1 = wcs.radec2pixelxy(rlo,dlo)
    ok2,tiles.x2,tiles.y2 = wcs.radec2pixelxy(rlo,dhi)
    ok3,tiles.x3,tiles.y3 = wcs.radec2pixelxy(rhi,dhi)
    ok4,tiles.x4,tiles.y4 = wcs.radec2pixelxy(rhi,dlo)

    ishallow = []

    shallowpcts = []

    t1,t2,t3 = targetdepth, targetdepth - 0.3, targetdepth - 0.6
    
    for itile,t in enumerate(tiles):
        # -1: FITS to numpy coords
        x0 = int(np.floor(min(t.x1, t.x2, t.x3, t.x4))) - 1
        y0 = int(np.floor(min(t.y1, t.y2, t.y3, t.y4))) - 1
        x1 = int(np.ceil( max(t.x1, t.x2, t.x3, t.x4))) - 1
        y1 = int(np.ceil( max(t.y1, t.y2, t.y3, t.y4))) - 1
        print('tile', t.ra, t.dec)
        #print('x0,x1, y0,y1', x0,x1,y0,y1)
        print('  x0,y0', x0,y0)
        print('  w,h', 1+x1-x0, 1+y1-y0)
        tiledepth = depth[y0:y1+1, x0:x1+1].copy()
        print('  tiledepth', tiledepth.min(), tiledepth.max())

        if tiledepth.min() > targetdepth:
            print('  Rectangular region is already to depth')
            continue
        
        # polygon
        xx,yy = np.meshgrid(np.arange(x0, x1+1), np.arange(y0, y1+1))
        poly = np.array([[t.x1,t.y1],[t.x2,t.y2],[t.x3,t.y3],[t.x4,t.y4]])
        inpoly = point_in_poly(xx, yy, poly-1)

        print(' ', np.sum(inpoly), 'pixels')
        tmin = tiledepth[inpoly].min()
        print('  Minimum depth:', tmin)
        if tmin > targetdepth:
            print('  Tile shaped region is already to depth')
            continue

        [d1,d2,d3] = np.percentile(tiledepth[inpoly], [10, 5, 2])
        print('  Depth at completeness 90/95/98:', d1, d2, d3)
        print('  Hitting target?',
              'yes' if d1 > t1 else 'no',
              'yes' if d2 > t2 else 'no',
              'yes' if d3 > t3 else 'no')

        if (d1 > t1) and (d2 > t2) and (d3 > t3):
            print('  Reached target!')
            continue

        ishallow.append(itile)
        shallowpcts.append((d1,d2,d3))
        
        # tiledepth[inpoly == False] = np.nan
        # plt.clf()
        # #plt.imshow(depth, interpolation='nearest', origin='lower')
        # plt.imshow(tiledepth, interpolation='nearest', origin='lower')#, cmap='RdBu')
        # #plt.axis([x0,x1+1,y0,y1+1])
        # plt.savefig('depth-tile%i.png' % itile)

    ishallow = np.array(ishallow)
    shallow = tiles[ishallow]
    shallow.depth_90 = np.array([d[0] for d in shallowpcts])
    shallow.depth_95 = np.array([d[1] for d in shallowpcts])
    shallow.depth_98 = np.array([d[2] for d in shallowpcts])
    for c in 'x1 y1 x2 y2 x3 y3 x4 y4'.split():
        shallow.delete_column(c)
    shallow.writeto('shallow.fits')
    
    print(len(ishallow), 'tiles deemed shallow')

    #print('shallow z_done:', Counter(shallow.z_done))
    
    plt.clf()
    plt.imshow(depth, **ima)
    sh = tiles[ishallow]
    plt.plot(np.vstack((sh.x1, sh.x2, sh.x3, sh.x4, sh.x1)),
             np.vstack((sh.y1, sh.y2, sh.y3, sh.y4, sh.y1)),
             'r-')
    plt.savefig('depth-2.png')

    plt.clf()
    for d1,d2,d3 in shallowpcts:
        plt.plot([d1,d2,d3], 'b-', alpha=0.02)
    for i,target in enumerate([t1,t2,t3]):
        plt.plot([i-0.25,i+0.25], [target,target], 'r-')
    plt.xticks([0,1,2], ['90%', '95%', '98%'])
    plt.xlabel('Coverage fraction')
    plt.ylabel('Depth of shallow tiles')
    plt.ylim(20, 23)
    plt.savefig('depth-3.png')
    

    plt.clf()
    plt.imshow(depth, **ima)
    sh = tiles[ishallow]
    sh.xm = (sh.x1 + sh.x2 + sh.x3 + sh.x4) / 4
    sh.ym = (sh.y1 + sh.y2 + sh.y3 + sh.y4) / 4
    plt.plot(np.vstack((sh.x1, sh.x2, sh.x3, sh.x4, sh.x1)),
             np.vstack((sh.y1, sh.y2, sh.y3, sh.y4, sh.y1)),
             'r-')
    plt.colorbar()
    for j,i in enumerate(np.random.permutation(len(sh))[:20]):
        plt.axis([sh.xm[i] - 50, sh.xm[i] + 50,
                  sh.ym[i] - 50, sh.ym[i] + 50])
        plt.savefig('depth-3-%02i.png' % j)

        
    
def depth_map_for_ccds(survey, ccds):
    W,H = 4800,3200
    plot = Plotstuff(size=(W,H), outformat='png')
    zoom = 2.7
    plot.wcs = anwcs_create_hammer_aitoff(195., 65., zoom, W, H, True)
    
    plot.color = 'black'
    plot.plot('fill')
    plot.color = 'white'
    plot.outline.fill = True
    plot.outline.stepsize = 2000
    plot.op = CAIRO_OPERATOR_ADD
    
    targetdepth = 22.5
    # This gives a 1-mag margin on the target depth
    maxfrac = 0.16
    
    targetsig = 10.**((targetdepth - 22.5) / -2.5)
    targetiv = 1./targetsig**2
    
    for j,ccd in enumerate(ccds):
        if j and j % 1000 == 0:
            print('Tile', j)
        depth = ccd.depth
        mwcs = survey.get_approx_wcs(ccd)

        detsig = 10.**((depth - 22.5) / -2.5)
        depthfrac = 1./detsig**2 / targetiv

        plot.outline.wcs = anwcs_new_tan(mwcs)
        plot.alpha = np.clip(maxfrac * depthfrac, 0., 1.)
        plot.apply_settings()
        plot.plot('outline')

    img = plot.get_image_as_numpy(flip=True)
    img = img[:,:,0]
    # It's a uint8; scale and convert to float
    img = img * 1./255.
    #print('Image', img.dtype, img.shape)
    img /= maxfrac
    # Now it's in factor of detiv of target depth.
    # back to sig -> depth
    img = -2.5 * (np.log10(1./np.sqrt(img * targetiv)) - 9.)
    img[np.logical_not(np.isfinite(img))] = 0.
    return img


def depth_map_for_tiles(tiles):
    '''
    Makes a depth map given a table of tile locations and depths

    tiles.depth
    tiles.ra
    tiles.dec
    '''
    ## Copied from from_ccds
    W,H = 4800,3200
    plot = Plotstuff(size=(W,H), outformat='png')
    zoom = 2.7
    plot.wcs = anwcs_create_hammer_aitoff(195., 65., zoom, W, H, True)

    plot.color = 'black'
    plot.plot('fill')
    plot.color = 'white'
    plot.outline.fill = True
    plot.outline.stepsize = 2000
    plot.op = CAIRO_OPERATOR_ADD

    targetdepth = 22.5
    # This gives a 1-mag margin on the target depth
    maxfrac = 0.16
    targetsig = 10.**((targetdepth - 22.5) / -2.5)
    targetiv = 1./targetsig**2

    detsig = 10.**((tiles.depth - 22.5) / -2.5)
    depthfrac = 1./detsig**2 / targetiv
    for r,d,df in zip(tiles.ra, tiles.dec, depthfrac):
        plot.outline.wcs = anwcs_new_tan(mosaic_wcs(r, d))
        plot.alpha = np.clip(maxfrac * df, 0., 1.)
        plot.apply_settings()
        plot.plot('outline')
    img = plot.get_image_as_numpy(flip=True)
    img = img[:,:,0]
    # It's a uint8; scale and convert to float
    img = img * 1./255.
    img /= maxfrac

    # back to sig -> depth
    img = -2.5 * (np.log10(1./np.sqrt(img * targetiv)) - 9.)
    img[np.logical_not(np.isfinite(img))] = 0.
    return img


def measure_map_at_tiles(tiles, wcs, depth):
    from astrometry.util.miscutils import point_in_poly
    # from mosaic_wcs
    tilesize = (4096 * 2 + 100) * 0.262 / 3600.

    dlo = tiles.dec - tilesize/2.
    dhi = tiles.dec + tilesize/2.
    cosdec = np.cos(np.deg2rad(tiles.dec))
    rlo = tiles.ra - tilesize/2./cosdec
    rhi = tiles.ra + tilesize/2./cosdec
    ok1,tiles.x1,tiles.y1 = wcs.radec2pixelxy(rlo,dlo)
    ok2,tiles.x2,tiles.y2 = wcs.radec2pixelxy(rlo,dhi)
    ok3,tiles.x3,tiles.y3 = wcs.radec2pixelxy(rhi,dhi)
    ok4,tiles.x4,tiles.y4 = wcs.radec2pixelxy(rhi,dlo)

    pcts = []

    for itile,t in enumerate(tiles):
        # -1: FITS to numpy coords
        x0 = int(np.floor(min(t.x1, t.x2, t.x3, t.x4))) - 1
        y0 = int(np.floor(min(t.y1, t.y2, t.y3, t.y4))) - 1
        x1 = int(np.ceil( max(t.x1, t.x2, t.x3, t.x4))) - 1
        y1 = int(np.ceil( max(t.y1, t.y2, t.y3, t.y4))) - 1
        #print('tile', t.ra, t.dec)
        #print('x0,x1, y0,y1', x0,x1,y0,y1)
        #print('  x0,y0', x0,y0)
        #print('  w,h', 1+x1-x0, 1+y1-y0)
        tiledepth = depth[y0:y1+1, x0:x1+1].copy()
        #print('  tiledepth', tiledepth.min(), tiledepth.max())
        # polygon
        xx,yy = np.meshgrid(np.arange(x0, x1+1), np.arange(y0, y1+1))
        poly = np.array([[t.x1,t.y1],[t.x2,t.y2],[t.x3,t.y3],[t.x4,t.y4]])
        inpoly = point_in_poly(xx, yy, poly-1)
        [d1,d2,d3] = np.percentile(tiledepth[inpoly], [10, 5, 2])
        #print('  Depth at completeness 90/95/98:', d1, d2, d3)
        pcts.append((d1,d2,d3))

    tiles.depth_90 = np.array([d[0] for d in pcts])
    tiles.depth_95 = np.array([d[1] for d in pcts])
    tiles.depth_98 = np.array([d[2] for d in pcts])
    #for c in 'x1 y1 x2 y2 x3 y3 x4 y4'.split():
    #    tiles.delete_column(c)

def djs_update():
    # Evaluate David's proposed tile file update
    from astrometry.util.miscutils import point_in_poly
    from astrometry.util.plotutils import PlotSequence
    from astrometry.libkd.spherematch import match_radec
    from legacypipe.survey import LegacySurveyData
    
    ps = PlotSequence('depth-djs')
    
    targetdepth = 22.5

    tiles = fits_table('mzls-update/mosaic-tiles_obstatus-update.fits')
    tiles.cut(tiles.in_desi == 1)
    tiles.cut(tiles.get('pass') <= 3)
    tiles.cut(tiles.dec > 30)
    print(len(tiles), 'in footprint above Dec 30, pass 1/2/3')
    tiles.recent = np.array([d > '2017-12-16' for d in tiles.z_date])
    print(sum(tiles.z_done == 0), 'tiles to-do')

    # Find tiles at the edge of the footprint.
    tiles.edge = np.zeros(len(tiles), bool)
    for passnum in [1,2,3]:
        Ipass = np.flatnonzero(tiles.get('pass') == passnum)
        tt = tiles[Ipass]
        radius = (4096*2 + 100) * np.sqrt(2) * 0.262 / 3600.
        I,J,d,count = match_radec(tt.ra, tt.dec, tt.ra, tt.dec, radius, notself=True,
                                  nearest=True, count=True)
        # plt.clf()
        # plt.hist(count)
        # plt.title('Number of nearby tiles for pass %i tiles' % passnum)
        # plt.savefig('depth-djs-8-%i.png' % passnum)
        thecounts = np.zeros(len(tt), int)
        thecounts[I] = count
        I = (thecounts < 6)
        #plt.clf()
        #plt.plot(tt.ra, tt.dec, 'k.')
        #plt.plot(tt.ra[I], tt.dec[I], 'r.')
        #plt.savefig('depth-djs-9-%i.png' % passnum)
        tiles.edge[Ipass[I]] = True
        
    fn = 'depth-djs-todo.fits'
    fn2 = 'depth-djs-recent.fits'
    if not os.path.exists(fn):
        # give them nominal depth
        tiles.depth = np.zeros(len(tiles)) + 22.2

        todo_depth = depth_map_for_tiles(tiles[tiles.z_done == 0])
        fitsio.write(fn, todo_depth)

        recent_depth = depth_map_for_tiles(tiles[tiles.recent *
                                                 (tiles.z_done != 0)])
        fitsio.write(fn2, recent_depth)
        #tiles.cut(np.logical_or(tiles.z_done == 0, tiles.recent))
        #print(len(tiles), 'with z_done = 0, or taken recently')
    else:
        todo_depth = fitsio.read(fn)
        recent_depth = fitsio.read(fn2)

    fn = 'depth-p9.fits'
    wcs = anwcs('plot.wcs')
    depth = fitsio.read(fn)
    print('Depth:', depth.shape, depth.dtype)

    fn = 'depth-dr6plus2.fits'
    if not os.path.exists(fn):
        survey = LegacySurveyData()
        ccds = fits_table(os.path.join(survey.survey_dir, 'ccds-annotated-mosaic-dr6plus2.fits.gz'))
        print(len(ccds), 'CCDs')
        ccds.depth = ccds.galdepth - ccds.decam_extinction[:,4]
        depth2 = depth_map_for_ccds(survey, ccds)
        fitsio.write(fn, depth2)
    else:
        depth2 = fitsio.read(fn)

    ima = dict(interpolation='nearest', origin='lower',
               vmin=22.0, vmax=23.0, cmap='RdBu')

    plt.figure(1, figsize=(13,8))
    plt.subplots_adjust(left=0.01, right=0.99, bottom=0.01, top=0.95)
    
    plt.clf()
    plt.imshow(depth, **ima)
    plt.colorbar()
    plt.title('Existing DR6plus tile depth')
    ps.savefig()

    plt.clf()
    plt.imshow(depth2, **ima)
    plt.colorbar()
    plt.title('Existing DR6plus2 tile depth')
    ps.savefig()

    plt.clf()
    plt.imshow(recent_depth, **ima)
    plt.colorbar()
    plt.title('Nominal depth for recent tiles')
    ps.savefig()

    plt.clf()
    plt.imshow(todo_depth, **ima)
    plt.colorbar()
    plt.title('Nominal depth for to-do tiles')
    ps.savefig()

    iv1 = 1./(10.**((depth - 22.5) / -2.5))**2
    iv2 = 1./(10.**((todo_depth  - 22.5) / -2.5))**2
    iv3 = 1./(10.**((recent_depth  - 22.5) / -2.5))**2
    iv4 = 1./(10.**((depth2  - 22.5) / -2.5))**2
    depth = -2.5 * (np.log10(1./np.sqrt(iv1 + iv2 + iv3 + iv4)) - 9.)

    plt.clf()
    plt.imshow(depth, **ima)
    plt.colorbar()
    plt.title('Projected total depth')
    ps.savefig()

    # Currently how deep are the to-do tiles?
    currdepth = -2.5 * (np.log10(1./np.sqrt(iv1 + iv3 + iv4)) - 9.)
    todo_tiles = tiles[tiles.z_done == 0]
    measure_map_at_tiles(todo_tiles, wcs, currdepth)

    I = np.flatnonzero((todo_tiles.depth_90 > targetdepth) *
                       (todo_tiles.depth_95 > targetdepth-0.3) *
                       (todo_tiles.depth_98 > targetdepth-0.6))
    print(len(I), 'to-do tiles already meet depth requirements')

    plt.clf()
    plt.imshow(currdepth, **ima)
    plt.colorbar()
    plt.title('Current depth + nominal depth for recent tiles')
    ps.savefig()
    sh = todo_tiles[I]
    plt.plot(np.vstack((sh.x1, sh.x2, sh.x3, sh.x4, sh.x1)),
             np.vstack((sh.y1, sh.y2, sh.y3, sh.y4, sh.y1)),
             'r-')
    plt.title('%i to-do tiles already pass depth requirements' % len(sh))
    ps.savefig()
    
    
    #tiles.cut(tiles.dec > 30)
    tiles.cut(np.lexsort((tiles.ra, tiles.dec)))

    #fn = 'djs-tile-depths.fits'
    #if not os.path.exists(fn):
    measure_map_at_tiles(tiles, wcs, depth)
    #tiles.writeto(fn)
    #else:
    #    tiles = fits_table(fn)

    tiles.shallow = np.logical_or(
        tiles.depth_90 < targetdepth,
        np.logical_or(tiles.depth_95 < targetdepth-0.3,
                      tiles.depth_98 < targetdepth-0.6))
    Ishallow = np.flatnonzero(tiles.shallow)
    print(len(Ishallow), 'shallow tiles (Dec > 30)')
    print(sum(tiles.shallow * (tiles.edge == False)), 'non-edge shallow tiles (Dec > 30)')

    print(sum(tiles.shallow * (tiles.z_done==0)), 'of these are to-do')
    
    plt.clf()
    plt.imshow(depth, **ima)
    sh = tiles[Ishallow]
    plt.plot(np.vstack((sh.x1, sh.x2, sh.x3, sh.x4, sh.x1)),
             np.vstack((sh.y1, sh.y2, sh.y3, sh.y4, sh.y1)),
             'r-')
    plt.colorbar()
    plt.title('Shallow tiles after projected depth: %i' % len(sh))
    ps.savefig()

    plt.clf()
    plt.imshow(depth, **ima)
    sh = tiles[tiles.shallow * (tiles.edge == False)]
    plt.plot(np.vstack((sh.x1, sh.x2, sh.x3, sh.x4, sh.x1)),
             np.vstack((sh.y1, sh.y2, sh.y3, sh.y4, sh.y1)),
             'r-')
    plt.colorbar()
    plt.title('Non-edge shallow tiles after projected depth: %i' % len(sh))
    ps.savefig()
    
    for passnum in [1,2,3]:
        plt.clf()
        plt.imshow(depth, **ima)
        sh = tiles[tiles.shallow * (tiles.edge == False)]
        sh.cut(sh.get('pass') == passnum)
        print(len(sh), 'non-edge shallow in pass', passnum)

        print('Of these,', sum(sh.bstarv < 20), 'have a bright star')

        plt.plot(np.vstack((sh.x1, sh.x2, sh.x3, sh.x4, sh.x1)),
                 np.vstack((sh.y1, sh.y2, sh.y3, sh.y4, sh.y1)),
                 'r-')
        plt.colorbar()
        plt.title('Pass %i shallow tiles after projected depth: %i' %
                  (passnum, len(sh)))
        ps.savefig()
        
    plt.subplots_adjust(left=0.1, right=0.99, bottom=0.1, top=0.95)

    lo,hi = 21.75, 23.25
    ha = dict(range=(lo, hi), bins=50, histtype='step')
    plt.clf()
    plt.hist(np.clip(tiles.depth_90, lo, hi), color='b',
             label='90% coverage', **ha)
    plt.hist(np.clip(tiles.depth_95, lo, hi), color='g',
             label='95% coverage', **ha)
    plt.hist(np.clip(tiles.depth_98, lo, hi), color='r',
             label='98% coverage', **ha)
    plt.axvline(targetdepth, color='b')
    plt.axvline(targetdepth-0.3, color='g')
    plt.axvline(targetdepth-0.6, color='r')
    plt.xlim(lo,hi)
    plt.legend()
    plt.title('Projected depths of all tiles')
    ps.savefig()

    plt.clf()
    plt.hist(np.clip(todo_tiles.depth_90, lo, hi), color='b',
             label='90% coverage', **ha)
    plt.hist(np.clip(todo_tiles.depth_95, lo, hi), color='g',
             label='95% coverage', **ha)
    plt.hist(np.clip(todo_tiles.depth_98, lo, hi), color='r',
             label='98% coverage', **ha)
    plt.axvline(targetdepth, color='b')
    plt.axvline(targetdepth-0.3, color='g')
    plt.axvline(targetdepth-0.6, color='r')
    plt.xlim(lo,hi)
    plt.legend()
    plt.title('Current depths of to-do tiles')
    ps.savefig()

    # LEGACY_SURVEY_DIR=~/legacypipe/py/dr6
    survey = LegacySurveyData()
    ccds = survey.get_annotated_ccds()
    ccds.depth = ccds.galdepth - ccds.decam_extinction[:,4]

    ok1,ccds.x1,ccds.y1 = wcs.radec2pixelxy(ccds.ra0, ccds.dec0)
    ok2,ccds.x2,ccds.y2 = wcs.radec2pixelxy(ccds.ra1, ccds.dec1)
    ok3,ccds.x3,ccds.y3 = wcs.radec2pixelxy(ccds.ra2, ccds.dec2)
    ok4,ccds.x4,ccds.y4 = wcs.radec2pixelxy(ccds.ra3, ccds.dec3)
    ccds.xmean = (ccds.x1 + ccds.x2 + ccds.x3 + ccds.x4) / 4.
    ccds.ymean = (ccds.y1 + ccds.y2 + ccds.y3 + ccds.y4) / 4.

    #bright = tiles[tiles.bstarv < 20.]
    
    tiles.xmean = (tiles.x1 + tiles.x2 + tiles.x3 + tiles.x4) / 4.
    tiles.ymean = (tiles.y1 + tiles.y2 + tiles.y3 + tiles.y4) / 4.
    sh = tiles[Ishallow]
    sh.cut(sh.get('pass') == 1)
    print(len(sh), 'shallow in pass 1')
    for ii,t in enumerate(sh[:20]):
        plt.clf()
        plt.subplot(1,2,1)
        plt.imshow(depth, **ima)
        margin = 30
        plt.plot(np.vstack((t.x1, t.x2, t.x3, t.x4, t.x1)),
                 np.vstack((t.y1, t.y2, t.y3, t.y4, t.y1)),
                 'r-')

        I = np.flatnonzero((np.abs(tiles.xmean - t.xmean) < margin) *
                           (np.abs(tiles.ymean - t.ymean) < margin))
        for i in I:
            plt.text(tiles.xmean[i], tiles.ymean[i], '%.2f' % tiles.z_depth[i],
                     color='k', ha='center', va='center',
                     bbox=dict(edgecolor='none', facecolor='white'))

        plt.axis([t.xmean-margin, t.xmean+margin,
                  t.ymean-margin, t.ymean+margin])

        bstr = ''
        if t.bstarv < 20:
            bstr = ', bright star: %.1f' % t.bstarv
        plt.suptitle('Tile %i, RA,Dec %.3f,%.3f, depths %.2f / %.2f / %.2f%s' %
                  (t.tileid, t.ra, t.dec, t.depth_90, t.depth_95, t.depth_98, bstr))


        plt.subplot(1,2,2)
        plt.imshow(depth, **ima)
        plt.plot(np.vstack((t.x1, t.x2, t.x3, t.x4, t.x1)),
                 np.vstack((t.y1, t.y2, t.y3, t.y4, t.y1)),
                 'r-')
        plt.text(t.xmean, t.ymean, '%.2f' % t.z_depth,
                 color='r', ha='center', va='center',
                 bbox=dict(edgecolor='none', facecolor='white'))
        I = np.flatnonzero((np.abs(ccds.xmean - t.xmean) < margin) *
                           (np.abs(ccds.ymean - t.ymean) < margin))
        for i in I:
            plt.text(ccds.xmean[i], ccds.ymean[i], '%.2f' % ccds.depth[i],
                     color='k', ha='center', va='center',
                     bbox=dict(edgecolor='none', facecolor='white'))
        plt.plot(np.vstack((ccds.x1[I], ccds.x2[I], ccds.x3[I], ccds.x4[I], ccds.x1[I])),
                 np.vstack((ccds.y1[I], ccds.y2[I], ccds.y3[I], ccds.y4[I], ccds.y1[I])),
                 'k-')

        plt.axis([t.xmean-margin, t.xmean+margin,
                  t.ymean-margin, t.ymean+margin])
        
        ps.savefig()
    
if __name__ == '__main__':
    #from_ccds()
    #tiles_todo()
    #needed_tiles()
    #main()
    djs_update()
