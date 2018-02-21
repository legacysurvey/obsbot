from __future__ import print_function
from astrometry.blind.plotstuff import *
from astrometry.util.util import *
from astrometry.util.fits import *
import pylab as plt
import matplotlib
from collections import Counter
from scipy.ndimage.filters import minimum_filter, median_filter
from glob import glob
from astrometry.libkd.spherematch import match_radec

# Last date included in the CCDs files
recent_date = '2018-02-11'

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

    from camera_mosaic import ephem_observer
    from obsbot import mjd_to_ephem_date, get_airmass
    from jnox import ra2hms, dec2dms
    import ephem
    print('Recomputing airmasses assuming MOSAIC camera...')
    T.airmass = np.zeros(len(T), np.float32)
    obs = ephem_observer()
    for i,(mjd,ra,dec) in enumerate(zip(T.mjd_obs, T.ra, T.dec)):
        obs.date = mjd_to_ephem_date(mjd)
        rastr  = ra2hms(ra)
        decstr = dec2dms(dec)
        ephemstr = str('%s,f,%s,%s,20' % ('name', rastr, decstr))
        etile = ephem.readdb(ephemstr)
        etile.compute(obs)
        T.airmass[i] = get_airmass(float(etile.alt))
    print('done computing airmasses')
    
    T.depth = T.galdepth - T.decam_extinction[:,4]

    T.seeing = T.fwhm * T.pixscale_mean
    kx = 0.10
    zp0 = 26.4
    T.transparency = 10.**(-0.4 * (zp0 - T.ccdzpt - kx*(T.airmass-1.)))

    plt.figure(1, figsize=(13,8))
    plt.subplots_adjust(left=0.01, right=0.99, bottom=0.01, top=0.95)
    
    plt.clf()
    plt.scatter(T.ra, T.dec, c=T.transparency, s=3)
    plt.colorbar()
    plt.savefig('transp.png')

    for passnum in [1,2,3]:
        plt.clf()
        I = np.flatnonzero(T.tilepass == passnum)
        plt.scatter(T.ra[I], T.dec[I], c=T.transparency[I], s=3)
        plt.colorbar()
        plt.savefig('transp-p%i.png' % passnum)

        plt.clf()
        plt.scatter(T.ra[I], T.dec[I], c=T.airmass[I], s=3)
        plt.colorbar()
        plt.savefig('airmass-p%i.png' % passnum)
    
    imgs = []
    seeings = []
    transps = []

    for passnum in [1,2,3, 0,  9]:
        if passnum == 9:
            I = []
            tt = 'All Passes'
        else:
            I = np.flatnonzero((T.tilepass == passnum) *
                               (T.depth > 10))
            tt = 'Pass %i' % passnum
        print(len(I), 'tiles for pass', passnum)

        cm = matplotlib.cm.viridis
        #cm = matplotlib.cm.jet
        cm = cmap_discretize(cm, 10)

        if passnum != 9:
            seeing,transp,wcs = best_seeing_map_for_ccds(survey, T[I])
            seeings.append(seeing)
            transps.append(transp)
        else:
            seeing = seeings[0]
            transp = transps[0]
            for s,t in zip(seeings, transps):
                # Handle zeros...
                sI = (seeing == 0)
                seeing[sI] = s[sI]
                sI = (s != 0)
                seeing[sI] = np.minimum(seeing[sI], s[sI])
                transp = np.maximum(transp, t)
        
        lo,hi = 0.5, 2.0
        if passnum == 9:
            plt.clf()
            plt.imshow(median_filter(seeing, size=3),
                       interpolation='nearest', origin='lower',
                       vmin=lo, vmax=hi, cmap=cm)
            plt.xticks([]); plt.yticks([])
            plt.title('Best seeing (unfiltered): %s' % tt)
            plt.colorbar()
            plt.savefig('depth-p%i-17.png' % passnum)
        
        plt.clf()
        plt.imshow(median_filter(seeing, size=3),
                   interpolation='nearest', origin='lower',
                   vmin=lo, vmax=hi, cmap=cm)
        plt.xticks([]); plt.yticks([])
        plt.title('Best seeing: %s' % tt)
        plt.colorbar()
        plt.savefig('depth-p%i-15.png' % passnum)

        plt.clf()
        lo,hi = 0.5, 1.1
        plt.imshow(transp, interpolation='nearest', origin='lower',
                   vmin=lo, vmax=hi, cmap=cm)
        plt.xticks([]); plt.yticks([])
        plt.title('Best transparency: %s' % tt)
        plt.colorbar()
        plt.savefig('depth-p%i-16.png' % passnum)
        
        if passnum != 9:
            depthmap,wcs = depth_map_for_ccds(survey, T[I])
            iv = 1./(10.**((depthmap - 22.5) / -2.5))**2
            imgs.append(iv)
        else:
            iv = imgs[0]
            for im in imgs[1:]:
                iv += im
            depthmap = -2.5 * (np.log10(1./np.sqrt(iv)) - 9.)

        #wcs_unflipped = anwcs_create_hammer_aitoff(195., 65., zoom, W, H, False)
        #anwcs_write(wcs_unflipped, 'plot.wcs')
        #hdr = fitsio.read_header('plot.wcs')
        fitsio.write('depth-p%i.fits' % passnum, depthmap, clobber=True)
        #header=hdr,
        print('Unique depths:', np.unique(depthmap.ravel()))

        fitsio.write('seeing-p%i.fits' % passnum, seeing, clobber=True)
        fitsio.write('transp-p%i.fits' % passnum, transp, clobber=True)
        
        H,W = wcs.shape
        
        plt.figure(2)

        plt.clf()
        lo,hi = 21.5,23.5
        hh = depthmap.ravel()
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
        ok,tiles.x,tiles.y = wcs.radec2pixelxy(tiles.ra, tiles.dec)
        tiles.x -= 1
        tiles.y -= 1
        tiles.cut(ok)
        
        plt.figure(1)
        plt.clf()
        depthmap[depthmap == 0] = np.nan

        plt.plot(tiles.x, H-tiles.y, 'k.', alpha=0.1)
        
        plt.imshow(depthmap, interpolation='nearest', origin='lower',
                   vmin=lo, vmax=hi, cmap=cm)
        plt.xticks([]); plt.yticks([])
        ax = plt.axis()
        H,W = depthmap.shape
        for d in range(30, 90, 10):
            rr = np.arange(0, 360, 1)
            dd = np.zeros_like(rr) + d
            ok,xx,yy = wcs.radec2pixelxy(rr, dd)
            plt.plot(xx, H-yy, 'k-', alpha=0.1)
        for r in range(0, 360, 30):
            dd = np.arange(20, 90, 1.)
            rr = np.zeros_like(dd) + r
            ok,xx,yy = wcs.radec2pixelxy(rr, dd)
            plt.plot(xx, H-yy, 'k-', alpha=0.1)
        plt.axis(ax)
        plt.colorbar()
        plt.title('Depth: %s' % tt)
        plt.savefig('depth-p%i-12.png' % passnum)


        #if passnum == 0:
        if True:
            goal = 22.5
            drange = 1.
        
            plt.clf()
            depthmap[depthmap == 0] = np.nan
            plt.imshow(depthmap, interpolation='nearest', origin='lower',
                       vmin=goal-drange, vmax=goal+drange, cmap='RdBu')
            plt.xticks([]); plt.yticks([])
            ax = plt.axis()
            H,W = depthmap.shape
            for d in range(30, 90, 10):
                rr = np.arange(0, 360, 1)
                dd = np.zeros_like(rr) + d
                ok,xx,yy = wcs.radec2pixelxy(rr, dd)
                plt.plot(xx, H-yy, 'k-', alpha=0.1)
            for r in range(0, 360, 30):
                dd = np.arange(20, 90, 1.)
                rr = np.zeros_like(dd) + r
                ok,xx,yy = wcs.radec2pixelxy(rr, dd)
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

    tiles.recent = np.array([d > recent_date for d in tiles.z_date])

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
    args = (195., 65., zoom, W, H, True)
    wcs = anwcs_create_hammer_aitoff(*args)
    plot.wcs = anwcs_create_hammer_aitoff(*args)
    
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
    return img, wcs


def best_seeing_map_for_ccds(survey, ccds):
    W,H = 4800,3200
    plot = Plotstuff(size=(W,H), outformat='png')
    zoom = 2.7
    args = (195., 65., zoom, W, H, True)
    wcs = anwcs_create_hammer_aitoff(*args)
    plot.wcs = anwcs_create_hammer_aitoff(*args)
    
    plot.color = 'black'
    plot.plot('fill')
    plot.color = 'white'
    plot.outline.fill = True
    plot.outline.stepsize = 2000
    plot.op = CAIRO_OPERATOR_LIGHTEN

    plot.alpha = 1.0

    for j,ccd in enumerate(ccds):
        if j and j % 1000 == 0:
            print('Tile', j)

        mwcs = survey.get_approx_wcs(ccd)
        sval = 0.7 / ccd.seeing
        tval = ccd.transparency
        
        plot.outline.wcs = anwcs_new_tan(mwcs)
        plot.rgb = (float(np.clip(sval, 0., 1.)),
                    float(np.clip(tval, 0., 1.)),
                    0.)
        plot.apply_settings()
        plot.plot('outline')

    img = plot.get_image_as_numpy(flip=True)
    # It's a uint8; scale and convert to float
    img = img * 1./255.

    simg = img[:,:,0]
    simg = 0.7 / simg
    simg[np.logical_not(np.isfinite(simg))] = 0.

    timg = img[:,:,1]
    
    return simg, timg, wcs


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


def measure_map_at_tiles(tiles, wcs, depth, pcts=[10,5,2]):
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

    depths = []
    for itile,t in enumerate(tiles):
        # -1: FITS to numpy coords
        x0 = int(np.floor(min(t.x1, t.x2, t.x3, t.x4))) - 1
        y0 = int(np.floor(min(t.y1, t.y2, t.y3, t.y4))) - 1
        x1 = int(np.ceil( max(t.x1, t.x2, t.x3, t.x4))) - 1
        y1 = int(np.ceil( max(t.y1, t.y2, t.y3, t.y4))) - 1
        tiledepth = depth[y0:y1+1, x0:x1+1]
        xx,yy = np.meshgrid(np.arange(x0, x1+1), np.arange(y0, y1+1))
        poly = np.array([[t.x1,t.y1],[t.x2,t.y2],[t.x3,t.y3],[t.x4,t.y4]])
        inpoly = point_in_poly(xx, yy, poly-1)
        p = np.percentile(tiledepth[inpoly], pcts)
        depths.append(p)
    return np.array(depths)

def djs_update():
    # Evaluate David's proposed tile file update
    from astrometry.util.miscutils import point_in_poly
    from astrometry.util.plotutils import PlotSequence
    from astrometry.libkd.spherematch import match_radec
    from legacypipe.survey import LegacySurveyData
    
    ps = PlotSequence('depth-djs')

    targetdepth = 22.5

    #tiles = fits_table('mzls-update/mosaic-tiles_obstatus-update.fits')
    tiles = fits_table('obstatus/mosaic-tiles_obstatus.fits')
    tiles.recent = np.array([d > recent_date for d in tiles.z_date])
    #tiles = fits_table('tiles-todo-dstn.fits')
    #tiles.recent = np.logical_or(tiles.recent, tiles.maybe_good)
    #tiles = fits_table('obstatus/mosaic-tiles_remaining.fits')

    tiles.cut(tiles.in_desi == 1)
    tiles.cut(tiles.get('pass') <= 3)
    tiles.cut(tiles.dec > 30)
    print(len(tiles), 'in footprint above Dec 30, pass 1/2/3')
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

    #fn = 'depth-djs-todo.fits'
    #fn2 = 'depth-djs-recent.fits'
    #if not os.path.exists(fn):
    if True:
        t = tiles[tiles.z_done == 0]
        # give them nominal depth
        t.depth = np.zeros(len(tiles)) + 22.2
        todo_depth = depth_map_for_tiles(t)
        #fitsio.write(fn, todo_depth, clobber=True)

        t = tiles[tiles.recent * (tiles.z_done == 1)]
        t.depth = np.zeros(len(tiles)) + 22.2
        recent_depth = depth_map_for_tiles(t)
        #fitsio.write(fn2, recent_depth, clobber=True)
    else:
        todo_depth = fitsio.read(fn)
        recent_depth = fitsio.read(fn2)

    fn = 'depth-p9.fits'
    wcs = anwcs('plot.wcs')
    depth = fitsio.read(fn)
    print('Depth:', depth.shape, depth.dtype)

    brightstars = fits_table('BrightStarCatalog_Vmaglt10.fits')
    ok,brightstars.x,brightstars.y = wcs.radec2pixelxy(brightstars.ra,
                                                       brightstars.dec)
    print('ok:', np.unique(ok))
    brightstars.cut(ok)
    brightstars.x -= 1
    brightstars.y -= 1
    print(len(brightstars), 'bright stars')

    #ima = dict(interpolation='nearest', origin='lower',
    #           vmin=22.0, vmax=23.0, cmap='RdBu')

    tilestyle = dict(color='1', linestyle='-')
    
    #cm = matplotlib.cm.jet
    cm = matplotlib.cm.RdBu
    cm = cmap_discretize(cm, 10)
    cm.set_bad('0.7')
    ima = dict(interpolation='nearest', origin='lower',
               vmin=22, vmax=23, cmap=cm)
    #vmin=22.0, vmax=22.5, cmap=cm)
    def plot_depth(depth):
        depth = depth.copy()
        depth[depth < 10] = np.nan
        plt.imshow(depth, **ima)
        plt.colorbar()
        plt.xticks([]); plt.yticks([])
        ax = plt.axis()
        H,W = depth.shape
        for d in range(30, 90, 10):
            rr = np.arange(0, 360, 1)
            dd = np.zeros_like(rr) + d
            ok,xx,yy = wcs.radec2pixelxy(rr, dd)
            #plt.plot(xx, H-yy, 'k-', alpha=0.1)
            plt.plot(xx, yy, 'k-', alpha=0.1)
        for r in range(0, 360, 30):
            dd = np.arange(20, 90, 1.)
            rr = np.zeros_like(dd) + r
            ok,xx,yy = wcs.radec2pixelxy(rr, dd)
            #plt.plot(xx, H-yy, 'k-', alpha=0.1)
            plt.plot(xx, yy, 'k-', alpha=0.1)
        plt.axis(ax)
        
    plt.figure(1, figsize=(13,8))
    plt.subplots_adjust(left=0.01, right=0.99, bottom=0.01, top=0.95)
    
    plt.clf()
    plot_depth(depth)
    plt.title('Existing DR6plusX tile depth')
    ps.savefig()

    plt.clf()
    plot_depth(recent_depth)
    plt.title('Nominal depth for recent tiles')
    ps.savefig()

    plt.clf()
    plot_depth(todo_depth)
    plt.title('Nominal depth for to-do tiles')
    ps.savefig()

    iv1 = 1./(10.**((depth - 22.5) / -2.5))**2
    iv2 = 1./(10.**((todo_depth  - 22.5) / -2.5))**2
    iv3 = 1./(10.**((recent_depth  - 22.5) / -2.5))**2
    depth = -2.5 * (np.log10(1./np.sqrt(iv1 + iv2 + iv3)) - 9.)

    plt.clf()
    plot_depth(depth)
    plt.title('Projected total depth')
    ps.savefig()

    # Currently how deep are the to-do tiles?
    currdepth = -2.5 * (np.log10(1./np.sqrt(iv1 + iv3)) - 9.)
    todo_tiles = tiles[tiles.z_done == 0]
    dd = measure_map_at_tiles(todo_tiles, wcs, currdepth)
    if len(dd) == 0:
        dd = np.zeros((0,3), np.float32)
    todo_tiles.depth_90 = dd[:,0]
    todo_tiles.depth_95 = dd[:,1]
    todo_tiles.depth_98 = dd[:,2]
    I = np.flatnonzero((todo_tiles.depth_90 > targetdepth) *
                       (todo_tiles.depth_95 > targetdepth-0.3) *
                       (todo_tiles.depth_98 > targetdepth-0.6))
    print(len(I), 'to-do tiles already meet depth requirements')

    tt = fits_table()
    tt.tileid = todo_tiles.tileid[I]
    tt.writeto('todo-done-tiles.fits')
    
    plt.clf()
    plot_depth(currdepth)
    plt.title('Current depth + nominal depth for recent tiles')
    ps.savefig()
    sh = todo_tiles[I]
    plt.plot(np.vstack((sh.x1, sh.x2, sh.x3, sh.x4, sh.x1)),
             np.vstack((sh.y1, sh.y2, sh.y3, sh.y4, sh.y1)),
             **tilestyle)
    plt.title('%i to-do tiles already pass depth requirements' % len(sh))
    ps.savefig()

    I = np.flatnonzero(np.logical_not(
        (todo_tiles.depth_90 > targetdepth) *
        (todo_tiles.depth_95 > targetdepth-0.3) *
        (todo_tiles.depth_98 > targetdepth-0.6)))
    print(len(I), 'to-do tiles do not already meet depth requirements')
    
    plt.clf()
    plot_depth(currdepth)
    sh = todo_tiles[I]
    plt.plot(np.vstack((sh.x1, sh.x2, sh.x3, sh.x4, sh.x1)),
             np.vstack((sh.y1, sh.y2, sh.y3, sh.y4, sh.y1)),
             **tilestyle)
    plt.title('%i to-do tiles do not already pass depth requirements' % len(sh))
    ps.savefig()

    
    #tiles.cut(tiles.dec > 30)
    tiles.cut(np.lexsort((tiles.ra, tiles.dec)))

    print('Measuring depths...')
    dd = measure_map_at_tiles(tiles, wcs, depth)
    tiles.depth_90 = dd[:,0]
    tiles.depth_95 = dd[:,1]
    tiles.depth_98 = dd[:,2]
    tiles.writeto('tile-depths.fits')

    tiles.shallow = np.logical_or(
        tiles.depth_90 < targetdepth,
        np.logical_or(tiles.depth_95 < targetdepth-0.3,
                      tiles.depth_98 < targetdepth-0.6))
    Ishallow = np.flatnonzero(tiles.shallow)
    print(len(Ishallow), 'shallow tiles (Dec > 30)')
    print(sum(tiles.shallow * (tiles.edge == False)), 'non-edge shallow tiles (Dec > 30)')

    print(sum(tiles.shallow * (tiles.z_done==0)), 'of these are to-do')
    
    plt.clf()
    plot_depth(depth)
    sh = tiles[Ishallow]
    plt.plot(np.vstack((sh.x1, sh.x2, sh.x3, sh.x4, sh.x1)),
             np.vstack((sh.y1, sh.y2, sh.y3, sh.y4, sh.y1)),
             **tilestyle)
    plt.title('Shallow tiles after projected depth: %i' % len(sh))
    ps.savefig()

    plt.clf()
    plot_depth(depth)
    sh = tiles[tiles.shallow * (tiles.edge == False)]
    plt.plot(np.vstack((sh.x1, sh.x2, sh.x3, sh.x4, sh.x1)),
             np.vstack((sh.y1, sh.y2, sh.y3, sh.y4, sh.y1)),
             **tilestyle)
    plt.title('Non-edge shallow tiles after projected depth: %i' % len(sh))
    ps.savefig()

    bright = tiles[tiles.bstarv < 20]

    #sh = tiles[tiles.shallow * (tiles.edge == False)]
    #print(len(sh), 'shallow tiles')
    
    for passnum in [1,2,3]:
        plt.clf()
        plot_depth(depth)
        sh = tiles[tiles.shallow * (tiles.edge == False)]
        sh.cut(sh.get('pass') == passnum)
        print(len(sh), 'non-edge shallow in pass', passnum)

        print('Of these,', sum(sh.bstarv < 20), 'have a bright star')

        I,J,d = match_radec(sh.ra, sh.dec, bright.ra, bright.dec, 1.0,
                            nearest=True)
        print('Of these,', len(I), 'have a tile containing a bright star within 1 degree')
        
        plt.plot(np.vstack((sh.x1, sh.x2, sh.x3, sh.x4, sh.x1)),
                 np.vstack((sh.y1, sh.y2, sh.y3, sh.y4, sh.y1)),
                 **tilestyle)
        #ax = plt.axis()
        #plt.plot(brightstars.x, brightstars.y, 'c.')
        #plt.axis(ax)
        plt.title('Pass %i shallow tiles after projected depth: %i' %
                  (passnum, len(sh)))
        ps.savefig()

        for c in 'x1 y1 x2 y2 x3 y3 x4 y4'.split():
            sh.delete_column(c)
        sh.writeto('shallow-p%i.fits' % passnum)
            
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

    sh = tiles[tiles.shallow * (tiles.edge == False) * (tiles.bstarv > 30.) * (tiles.z_done == 1)]
    print(len(sh), 'tiles predicted to be under-depth (and not to-do, and without bright star)')

    # sort by pass
    sh = sh[np.argsort(sh.get('pass'))]

    #for ii,t in enumerate(sh):

        # Search for nearby tiles
        # Try re-taking them in order of depth and keep the first one that works
        
        #iv3 = 1./(10.**((recent_depth  - 22.5) / -2.5))**2
        

    depth_p1 = fitsio.read('depth-p1.fits')
    depth_p2 = fitsio.read('depth-p2.fits')
    depth_p3 = fitsio.read('depth-p3.fits')
    depth_p = { 1: depth_p1,
                2: depth_p2,
                3: depth_p3 }

    cm = 'bone'
    
    ima = dict(interpolation='nearest', origin='lower',
               vmin=22.0, vmax=23, cmap=cm)

    ima2 = dict(interpolation='nearest', origin='lower',
                vmin=21.5, vmax=22.5, cmap=cm)

    tilestyle = dict(color='r', linestyle='-', linewidth=3)
    
    plt.subplots_adjust(left=0.01, right=0.99, bottom=0.01, top=0.9)


    remain = fits_table('obstatus/mosaic-tiles_remaining.fits')
    print('Remaining tiles:', remain)
    rids = set(remain.tileid)
    print('rids:', len(rids))
    I, = np.nonzero([tid in rids for tid in tiles.tileid])
    plot_tiles = tiles[I]
    print(len(plot_tiles), 'remaining tiles to plot')

    plot_tiles.writeto('remaining-updated.fits')

    return
    
    # sh.cut(sh.get('pass') == 1)
    # print(len(sh), 'shallow in pass 1')
    #for ii,t in enumerate(sh[:20]):
    for ii,t in enumerate(plot_tiles):
        margin = 30
        bstr = ''
        if t.bstarv < 20:
            bstr = ', bright star: %.1f' % t.bstarv
        # plt.clf()
        # plt.subplot(1,2,1)
        # plt.imshow(depth, **ima)
        # plt.colorbar()
        # plt.plot(np.vstack((t.x1, t.x2, t.x3, t.x4, t.x1)),
        #          np.vstack((t.y1, t.y2, t.y3, t.y4, t.y1)),
        #          **tilestyle)
        # 
        # I = np.flatnonzero((np.abs(tiles.xmean - t.xmean) < margin) *
        #                    (np.abs(tiles.ymean - t.ymean) < margin))
        # for i in I:
        #     plt.text(tiles.xmean[i], tiles.ymean[i], '%.2f' % tiles.z_depth[i],
        #              color='k', ha='center', va='center',
        #              bbox=dict(edgecolor='none', facecolor='white'))
        # 
        # plt.axis([t.xmean-margin, t.xmean+margin,
        #           t.ymean-margin, t.ymean+margin])
        # 
        # plt.suptitle('Tile %i, RA,Dec %.3f,%.3f, depths %.2f / %.2f / %.2f%s, Z_done %s, Z_expnum %i' %
        #           (t.tileid, t.ra, t.dec, t.depth_90, t.depth_95, t.depth_98, bstr, t.z_done == 1, t.z_expnum))
        # 
        # 
        # plt.subplot(1,2,2)
        # plt.imshow(depth, **ima)
        # plt.colorbar()
        # plt.plot(np.vstack((t.x1, t.x2, t.x3, t.x4, t.x1)),
        #          np.vstack((t.y1, t.y2, t.y3, t.y4, t.y1)),
        #          **tilestyle)
        # plt.text(t.xmean, t.ymean, '%.2f' % t.z_depth,
        #          color='r', ha='center', va='center',
        #          bbox=dict(edgecolor='none', facecolor='white'))
        # I = np.flatnonzero((np.abs(ccds.xmean - t.xmean) < margin) *
        #                    (np.abs(ccds.ymean - t.ymean) < margin))
        # for i in I:
        #     plt.text(ccds.xmean[i], ccds.ymean[i], '%.2f' % ccds.depth[i],
        #              color='k', ha='center', va='center',
        #              bbox=dict(edgecolor='none', facecolor='white'))
        # plt.plot(np.vstack((ccds.x1[I], ccds.x2[I], ccds.x3[I], ccds.x4[I], ccds.x1[I])),
        #          np.vstack((ccds.y1[I], ccds.y2[I], ccds.y3[I], ccds.y4[I], ccds.y1[I])),
        #          'k-')
        # 
        # plt.axis([t.xmean-margin, t.xmean+margin,
        #           t.ymean-margin, t.ymean+margin])
        # 
        # ps.savefig()

        plt.clf()
        plt.subplot(2,2,1)
        plt.imshow(depth, **ima)
        plt.xticks([]); plt.yticks([])
        plt.colorbar()
        plt.plot(np.vstack((t.x1, t.x2, t.x3, t.x4, t.x1)),
                 np.vstack((t.y1, t.y2, t.y3, t.y4, t.y1)),
                 **tilestyle)
        plt.plot(brightstars.x, brightstars.y, 'ro')
        plt.axis([t.xmean-margin, t.xmean+margin,
                  t.ymean-margin, t.ymean+margin])
        plt.title('Total depth')

        for p in [1,2,3]:
            plt.subplot(2,2,p+1)
            plt.imshow(depth_p[p], **ima2)
            plt.colorbar()
            plt.xticks([]); plt.yticks([])
            margin = 30
            plt.plot(np.vstack((t.x1, t.x2, t.x3, t.x4, t.x1)),
                     np.vstack((t.y1, t.y2, t.y3, t.y4, t.y1)),
                     **tilestyle)

            I = np.flatnonzero((np.abs(ccds.xmean - t.xmean) < margin) *
                               (np.abs(ccds.ymean - t.ymean) < margin) *
                               (ccds.tilepass == p))
            expnums = np.unique(ccds.expnum[I])
            for expnum in expnums:
                J = I[ccds.expnum[I] == expnum]

                plt.text(np.mean(ccds.xmean[J]), np.mean(ccds.ymean[J]),
                         '%.2f' % np.median(ccds.depth[J]),
                         color='k', ha='center', va='center',
                         bbox=dict(edgecolor='none', facecolor='white'))
            
            plt.axis([t.xmean-margin, t.xmean+margin,
                      t.ymean-margin, t.ymean+margin])
            plt.title('Pass %i depth' % p)

        plt.suptitle('Tile %i, Pass %i, RA,Dec %.3f,%.3f, depths %.2f / %.2f / %.2f%s, Z_done %s, Z_expnum %i' %
                  (t.tileid, t.get('pass'), t.ra, t.dec, t.depth_90, t.depth_95, t.depth_98, bstr, t.z_done == 1, t.z_expnum))
        ps.savefig()


        #### Search for nearby tiles, look at their depths (90th
        # pct?), starting with the shallowest (up to some cut?), ask
        # if observing them would satisfy this tile.







        
        
def when_missing():
    from legacyzpts.psfzpt_cuts import *

    bad_expid = read_bad_expid()
    z0 = 26.20
    dz = (-0.6, 0.6)

    fns = glob(os.path.join(os.environ['LEGACY_SURVEY_DIR'],
                            'psfzpts-pre-cuts-mosaic*.fits'))
    fns.sort()
    print('Reading', fns)
    P = merge_tables([fits_table(fn) for fn in fns])
    psf_zeropoint_cuts(P, ['z'], 0.262, z0+dz[0], z0+dz[1], bad_expid, 'mosaic')
    print(len(P), 'PSF zeropoints')
    print(np.sum(P.ccd_cuts == 0), 'PSF zeropoints pass cuts')

    fns = glob(os.path.join(os.environ['LEGACY_SURVEY_DIR'],
                            'ccds-annotated-mosaic*.fits.gz'))
    fns.sort()
    print('Reading', fns)
    A = merge_tables([fits_table(fn) for fn in fns])
    
    passnum = 3
    depthmap = fitsio.read('depth-p%i.fits' % passnum)

    tiles = fits_table('obstatus/mosaic-tiles_obstatus.fits')
    tiles.cut(tiles.in_desi == 1)
    tiles.cut(tiles.get('pass') <= 3)
    tiles.cut(tiles.dec > 30)
    print(len(tiles), 'in footprint above Dec 30, pass 1/2/3')

    tiles.cut(tiles.get('pass') == passnum)
    print(len(tiles), 'in pass', passnum)
    
    print('Z_DATEs:', np.unique(tiles.z_date))

    plt.figure(1, figsize=(13,8))
    plt.subplots_adjust(left=0.01, right=0.99, bottom=0.01, top=0.95)

    depthmap[depthmap == 0] = np.nan

    wcs = anwcs_open('plot.wcs', 0)
    
    for date in np.unique(tiles.z_date):
        #if date < '2017-12-01':
        if date != '2017-12-22':
            continue
        print('Date', date)

        I = np.flatnonzero(tiles.z_date == date)
        print(len(I), 'tiles marked done on', date)

        #cm = matplotlib.cm.viridis
        #cm = cmap_discretize(cm, 10)
        #cm = 'Greys'
        cm = 'bone'
        
        plt.clf()
        lo,hi = 21.5,23.5
        plt.imshow(depthmap, interpolation='nearest', origin='lower',
                   vmin=lo, vmax=hi, cmap=cm)
        plt.xticks([]); plt.yticks([])
        plt.colorbar()

        p1 = p2 = p3 = None

        for i in I:
            mwcs = mosaic_wcs(tiles.ra[i], tiles.dec[i])
            H,W = mwcs.shape
            xx = np.array([1,1,W,W,1])
            yy = np.array([1,H,H,1,1])
            rr,dd = mwcs.pixelxy2radec(xx, yy)
            ok,xx,yy = wcs.radec2pixelxy(rr, dd)

            expnum = tiles.z_expnum[i]
            J = np.flatnonzero(P.expnum == expnum)
            #print('expnum', expnum)
            #print(len(J), 'entries in PSF zeropoints file')
            if len(J) > 0:
                if np.all(P.ccd_cuts[J] == 0):
                    p1 = plt.plot(xx, yy, 'g-', label='CCD cuts = good')
                    #print('Expnum', expnum, 'CCD cuts:', P.ccd_cuts[J])
                    K = np.flatnonzero(A.expnum == expnum)
                    #print(len(K), 'found in annotated-CCDs file.')
                    if len(K):
                        print('Expnum', expnum, 'depths:', A.galdepth[K])
                        #print('Depths', A.galdepth[K])
                    else:
                        print('Expnum', expnum, 'CCD cuts:', np.unique(P.ccd_cuts[J]), 'not found in annotated-CCDs')
                else:
                    print('Expnum', expnum, 'CCD cuts:', P.ccd_cuts[J])
                    p2 = plt.plot(xx, yy, 'r-', label='CCD cuts = bad')
            else:
                print('Expnum', expnum, 'not found in PSF Zeropoints file')
                p3 = plt.plot(xx, yy, 'c-', label='No CCD cuts available')
        plt.title('Pass %i: observed on %s' % (passnum, date))
        leg = []
        for p in [p1,p2,p3]:
            if p is not None:
                leg.append(p[0])
        plt.legend(leg)
        plt.savefig('observed-%s.png' % date)
        
def update_tiles():
    from legacyzpts.psfzpt_cuts import read_bad_expid, psf_zeropoint_cuts, psf_cuts_to_string

    tiles = fits_table('obstatus/mosaic-tiles_obstatus.fits')

    brightstars = fits_table('BrightStarCatalog_Vmaglt10.fits')

    from legacypipe.survey import LegacySurveyData
    survey = LegacySurveyData()

    ann = survey.get_annotated_ccds()
    print(len(ann), 'annotated CCDs')
    ann.depth = ann.galdepth - ann.decam_extinction[:,4]

    # What about the annotated-CCDs entries with tileid = 0?
    #ann[ann.tileid == 0].writeto('annotated-tileid-0.fits')
    # These turn out to be from different propids, etc, and not on our
    # tile centers.
    
    # max degrees between CCD center and its tile center
    max_offset = 0.24

    from astrometry.util.starutil_numpy import distsq2deg,radectoxyz
    tileidmap = dict([(tid,i) for i,tid in enumerate(tiles.tileid)])
    tileindex = np.array([tileidmap.get(tid,-1) for tid in ann.tileid])
    K = np.flatnonzero(tileindex >= 0)
    # sigh, my degrees_between function doesn't broadcast correctly, so DIY
    xyz1 = radectoxyz(ann.ra[K], ann.dec[K])
    xyz2 = radectoxyz(tiles.ra[tileindex[K]], tiles.dec[tileindex[K]])
    d2 = np.sum((xyz1 - xyz2)**2, axis=1)
    d = distsq2deg(d2)
    ann.astrom_ok = np.zeros(len(ann), bool)
    ann.astrom_ok[K] = (d < max_offset)
    print('Cut by bad astrom:', np.sum(d >= max_offset))
    print('Good astrom:', Counter(ann.astrom_ok))
    print('Bad astrom:', len(np.unique(ann.expnum[ann.astrom_ok == False])), 'exposures')
    print(np.sum(ann.astrom_ok), 'with good astrom')

    ann.cut(ann.astrom_ok)
    
    ann.ccdnum = np.zeros(len(ann), np.int16)
    for i,ccdname in enumerate(ann.ccdname):
        ccdnum = int(ccdname.replace('CCD','').replace('ccd',''), 10)
        assert(ccdnum >= 1)
        assert(ccdnum <= 4)
        ann.ccdnum[i] = ccdnum
    
    ccds = merge_tables([fits_table(fn) for fn in glob(os.path.join(survey.survey_dir, 'psfzpts-pre-*.fits*'))], columns='fillzero')
    ccds.filter = np.array([f.strip() for f in ccds.filter])

    ccds.ccdnum = np.zeros(len(ccds), np.int16)
    for i,ccdname in enumerate(ccds.ccdname):
        try:
            ccdnum = int(ccdname.replace('CCD','').replace('ccd',''), 10)
        except:
            #print('Did not parse CCDNAME', ccdname)
            continue
        assert(ccdnum >= 1)
        assert(ccdnum <= 4)
        ccds.ccdnum[i] = ccdnum
    
    print(len(ccds), 'CCDs')
    print('Filters:', Counter(ccds.filter).most_common())
    ccds.cut(ccds.filter == 'z')
    print(len(ccds), 'z-band CCDs')
    
    bad_expid = read_bad_expid()
    z0 = 26.20
    dz = (-0.6, 0.6)
    psf_zeropoint_cuts(ccds, 0.262, dict(z=z0+dz[0]), dict(z=z0+dz[1]), bad_expid, 'mosaic')

    print(len(ccds), 'pre-cuts CCDs')

    # Update depths.
    tiles.best_ccd_depths = np.zeros((len(tiles), 4), np.float32)
    tiles.n_exp = np.zeros(len(tiles), np.int16)

    print('tileids:', Counter(ann.tileid).most_common(10))
    
    tileidmap = dict()
    for i,tid in enumerate(ann.tileid):
        if tid == 0:
            continue
        if tid not in tileidmap:
            tileidmap[tid] = [i]
        else:
            tileidmap[tid].append(i)
    print('Built tile id map for', len(tileidmap), 'tile ids')

    II = np.flatnonzero((tiles.in_desi == 1) * (tiles.get('pass') <= 3) * (tiles.dec > 30))
    print('Checking', len(II), 'tiles')

    tiles.expnums_ann = np.zeros((len(tiles), 10), np.int32)
    
    unmatched = []
    for i in II:
        ai = tileidmap.get(tiles.tileid[i])
        if ai is None:
            unmatched.append(i)
            #print('No matches for tileid', tiles.tileid[i])
            continue
        ai = np.array(ai)
        tiles.n_exp[i] = len(np.unique(ann.expnum[ai]))
        ee = np.unique(ann.expnum[ai])
        tiles.expnums_ann[i,:len(ee)] = ee[:10]
        for ii in ai:
            tiles.best_ccd_depths[i, ann.ccdnum[ii]-1] = max(tiles.best_ccd_depths[i, ann.ccdnum[ii]-1],
                                                             ann.depth[ii])
    print(len(unmatched), 'tiles without matches in annotated-CCDs file')

    tiles.n_good_ccds = np.sum(tiles.best_ccd_depths > 0, axis=1)
    tiles.depth_median = np.median(tiles.best_ccd_depths, axis=1)
    tiles.depth_max = np.max(tiles.best_ccd_depths, axis=1)

    ccdidmap = dict()
    for i,obj in enumerate(ccds.object):
        words = obj.split('_')
        if len(words) != 3:
            continue
        tid = int(words[1], 10)
        if not tid in ccdidmap:
            ccdidmap[tid] = [i]
        else:
            ccdidmap[tid].append(i)

    still_unmatched = []

    flags = {}
    expnums = {}
    
    for i in unmatched:
        # Look for matches in the PSF-zeropoints files.
        tid = tiles.tileid[i]
        if not tid in ccdidmap:
            print('Tile id', tid, 'z_date', tiles.z_date[i], 'expnum', tiles.z_expnum[i], 'not in CCDs files')
            still_unmatched.append(i)
            continue
        #print('Tile id', tid, ': CCD cuts:', [psf_cuts_to_string(ccds.ccd_cuts[i]) for i in ccdidmap[tid]])

        ci = ccdidmap[tid]
        # byccdnum = {1:[],2:[],3:[],4:[]}
        # for ii in ci:
        #     byccdnum[ccds.ccdnum[ii]].append(psf_cuts_to_string(ccds.ccd_cuts[ii], join='|'))
        # print('Tile id', tid, 'cuts by CCDNUM:', byccdnum)

        byexpnum = {}
        for ii in ci:
            expnum = ccds.expnum[ii]
            if expnum in byexpnum:
                byexpnum[expnum] |= ccds.ccd_cuts[ii]
            else:
                byexpnum[expnum] = ccds.ccd_cuts[ii]
        bstar = ''
        if tiles.bstarv[i] < 10:
            bstar = 'Bright star: %.1f' % (tiles.bstarv[i])
        print('Tile id', tid, 'z_date', tiles.z_date[i], 'expnums', sorted(byexpnum.keys()),
              'with flags by exposure:', [psf_cuts_to_string(v, join='|') for v in byexpnum.values()], bstar)
        flags[i] = byexpnum.values()
        expnums[i] = byexpnum.keys()

    max_expnums = max([len(v) for v in flags.values()])
    tiles.expnums_zpts = np.zeros((len(tiles), max_expnums), np.int32)
    tiles.flags_zpts = np.zeros((len(tiles), max_expnums), np.int32)
    tiles.maybe_good = np.zeros(len(tiles), bool)
    for i in unmatched:
        if not i in flags:
            continue
        ee = expnums[i]
        ff = flags[i]
        #print('tile', tiles.tileid[i], 'expnums', ee, 'flags', ff)
        tiles.expnums_zpts[i,:len(ee)] = np.array(ee)
        tiles.flags_zpts[i,:len(ff)] = np.array(ff)

        from legacyzpts.psfzpt_cuts import CCD_CUT_BITS
        
        maybe = False
        for e,f in zip(ee,ff):
            # Images marked with only "not_third_pix" might still be okay.
            if f & ~CCD_CUT_BITS['not_third_pix'] == 0:
                maybe = True
        tiles.maybe_good[i] = maybe

    print(np.sum(tiles.maybe_good), 'tiles with matches in the zeropoints files might be okay (third-pixel)')
    
    print(len(still_unmatched), 'still umatched')

    allbits = 0
    for v in flags.values():
        for vv in v:
            allbits |= vv
    print('All flags set:', psf_cuts_to_string(allbits))

    tiles.recent = np.zeros(len(tiles), bool)
    for i in II:
        if tiles.z_date[i] > recent_date:
            #print('Tile date', tiles.z_date[i], 'is recent')
            if tiles.z_expnum[i] in bad_expid:
                pass
                #print('But in bad_expid')
            else:
                tiles.recent[i] = True

    still = []
    for i in still_unmatched:
        if not tiles.recent[i]:
            still.append(i)
    still_unmatched = still
                
    print('After dropping recent tiles, still', len(still_unmatched), 'unmatched')
    
    I = np.argsort(np.array([tiles.z_expnum[i] for i in still_unmatched]))
    still_unmatched = [still_unmatched[i] for i in I]
    
    print('Unique Z_DATEs:', np.unique([tiles.z_date[i] for i in still_unmatched]))
    print('Z_EXPNUMs:', [tiles.z_expnum[i] for i in still_unmatched])

    for i in still_unmatched:
        expnum = tiles.z_expnum[i]
        date = tiles.z_date[i]
        # print()
        # print('Tile', tiles.tileid[i])
        # print('Expnum', expnum)
        # print('Date', date)
        I = np.flatnonzero(ann.expnum == expnum)
        #print(len(I), 'matches to annotated')
        assert(len(I) == 0)
        if expnum < 100000:
            I = np.flatnonzero(ann.expnum == expnum + 300000)
            #print(len(I), 'matches to annotated 3+')
            assert(len(I) == 0)
        if expnum > 100000:
            I = np.flatnonzero(ann.expnum == expnum + 3000000)
            assert(len(I) == 0)
            #print(len(I), 'matches to annotated 3+')
        I = np.flatnonzero(ccds.expnum == expnum)
        assert(len(I) == 0)
        #print(len(I), 'matches to ccds')
        if expnum < 100000:
            I = np.flatnonzero(ccds.expnum == expnum + 300000)
            assert(len(I) == 0)
            #print(len(I), 'matches to ccds 3+')
        if expnum > 100000:
            I = np.flatnonzero(ccds.expnum == expnum + 3000000)
            assert(len(I) == 0)
            #print(len(I), 'matches to ccds 3+')

    II = np.flatnonzero((tiles.in_desi == 1) * (tiles.get('pass') <= 3) * (tiles.dec > 30))
    print('Checking', len(II), 'tiles')

    print('In annotated-CCDs:', Counter((tiles.n_exp[II] > 0)))
    print('In annotated-CCDs, or third-pixel:',
          Counter((tiles.n_exp[II] > 0) | (tiles.maybe_good[II])))
    print('In annotated-CCDs, or recent:',
          Counter((tiles.n_exp[II] > 0) | (tiles.recent[II])))
    good = ((tiles.n_exp[II] > 0) |
            (tiles.maybe_good[II]) |
            (tiles.recent[II]))
    print('In annotated-CCDs, or third-pixel, or recent:', Counter(good))

    # Find tiles near at edge of the footprint
    tiles.edge = np.zeros(len(tiles), bool)
    for passnum in [1,2,3]:
        Ipass = np.flatnonzero((tiles.get('pass') == passnum) *
                               (tiles.in_desi == 1) *
                               (tiles.dec > 30))
        tt = tiles[Ipass]
        radius = (4096*2 + 100) * np.sqrt(2) * 0.262 / 3600.
        I,J,d,count = match_radec(tt.ra, tt.dec, tt.ra, tt.dec, radius, notself=True,
                                  nearest=True, count=True)
        thecounts = np.zeros(len(tt), int)
        thecounts[I] = count
        I = (thecounts < 6)
        tiles.edge[Ipass[I]] = True

    print('Edge tiles of interest:', Counter(tiles.edge[II]))

    # Measure depth percentiles for tiles in depth maps
    wcs = anwcs('plot.wcs')
    depth = fitsio.read('depth-p9.fits')
    seeing = fitsio.read('seeing-p9.fits')
    transp = fitsio.read('transp-p9.fits')

    # todo = tiles[II[good == False]]
    # print(len(todo), 'to-do tiles')
    # recent = tiles[II[tiles.recent[II]]]
    # print(len(recent), 'recent tiles')
    # # Assume nominal depth for to-do and recent tiles.
    # todo.depth = np.zeros(len(todo)) + 22.2
    # todo_depth = depth_map_for_tiles(todo)
    # recent.depth = np.zeros(len(recent)) + 22.2
    # recent_depth = depth_map_for_tiles(recent)
    # 
    # # Combine depths
    # iv1 = 1./(10.**((depth - 22.5) / -2.5))**2
    # iv2 = 1./(10.**((recent_depth - 22.5) / -2.5))**2
    # iv3 = 1./(10.**((todo_depth - 22.5) / -2.5))**2
    # depth = -2.5 * (np.log10(1./np.sqrt(iv1 + iv2 + iv3)) - 9.)
    
    print('Measuring depth map...')
    dd = measure_map_at_tiles(tiles[II], wcs, depth)
    tiles.depth_90 = np.zeros(len(tiles), np.float32)
    tiles.depth_95 = np.zeros(len(tiles), np.float32)
    tiles.depth_98 = np.zeros(len(tiles), np.float32)
    tiles.depth_90[II] = dd[:,0]
    tiles.depth_95[II] = dd[:,1]
    tiles.depth_98[II] = dd[:,2]

    # Measure best seeing per tile area
    print('Measuring seeing map...')
    dd = measure_map_at_tiles(tiles[II], wcs, seeing, pcts=[50,90])
    tiles.seeing_50 = np.zeros(len(tiles), np.float32)
    tiles.seeing_90 = np.zeros(len(tiles), np.float32)
    tiles.seeing_50[II] = dd[:,0]
    tiles.seeing_90[II] = dd[:,1]

    # Measure best transparency per tile area
    print('Measuring transparency map...')
    dd = measure_map_at_tiles(tiles[II], wcs, transp, pcts=[50,90])
    tiles.transp_50 = np.zeros(len(tiles), np.float32)
    tiles.transp_90 = np.zeros(len(tiles), np.float32)
    tiles.transp_50[II] = dd[:,0]
    tiles.transp_90[II] = dd[:,1]

    tiles.writeto('tiles-updated.fits.gz')
    
    #tiles.z_done[II[np.logical_not(good)]] = 0
    #tiles.writeto('tiles-todo-dstn.fits')

    redo = np.logical_or(good == False,
                         (tiles.edge[II] == False) *
                         np.logical_or(tiles.seeing_50[II] > 1.4,
                                       np.logical_or(tiles.transp_50[II] < 0.75,
                                                     tiles.depth_90[II] < 22.4)))
    tiles[II[redo]].writeto('tiles-todo-dstn.fits')

    


def plot_tiles():
    from legacypipe.survey import LegacySurveyData
    from astrometry.util.plotutils import PlotSequence

    plt.figure(1, figsize=(13,8))
    fig = plt.figure(1)

    tiles = fits_table('tiles-updated.fits.gz')
    tiles.cut(tiles.in_desi == 1)
    tiles.cut(tiles.get('pass') <= 3)
    tiles.cut(tiles.dec > 30)
    print(len(tiles), 'in footprint above Dec 30, pass 1/2/3')

    tiles[tiles.seeing_50 > 1.4].writeto('tiles-todo-seeing.fits')
    tiles[tiles.transp_50 < 0.75].writeto('tiles-todo-transparency.fits')

    ps = PlotSequence('remain')
    
    wcs = anwcs('plot.wcs')
    
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
    tiles.xmean = (tiles.x1 + tiles.x2 + tiles.x3 + tiles.x4) / 4.
    tiles.ymean = (tiles.y1 + tiles.y2 + tiles.y3 + tiles.y4) / 4.
    tiles.cut(np.lexsort((tiles.ra, tiles.dec)))
    
    fn = 'depth-p9.fits'
    depth = fitsio.read(fn)

    depth_p1 = fitsio.read('depth-p1.fits')
    depth_p2 = fitsio.read('depth-p2.fits')
    depth_p3 = fitsio.read('depth-p3.fits')
    depth_p = { 1: depth_p1,
                2: depth_p2,
                3: depth_p3 }
    
    brightstars = fits_table('BrightStarCatalog_Vmaglt10.fits')
    ok,brightstars.x,brightstars.y = wcs.radec2pixelxy(brightstars.ra,
                                                       brightstars.dec)

    survey = LegacySurveyData()
    ccds = survey.get_annotated_ccds()
    ccds.depth = ccds.galdepth - ccds.decam_extinction[:,4]
    ok1,ccds.x1,ccds.y1 = wcs.radec2pixelxy(ccds.ra0, ccds.dec0)
    ok2,ccds.x2,ccds.y2 = wcs.radec2pixelxy(ccds.ra1, ccds.dec1)
    ok3,ccds.x3,ccds.y3 = wcs.radec2pixelxy(ccds.ra2, ccds.dec2)
    ok4,ccds.x4,ccds.y4 = wcs.radec2pixelxy(ccds.ra3, ccds.dec3)
    ccds.xmean = (ccds.x1 + ccds.x2 + ccds.x3 + ccds.x4) / 4.
    ccds.ymean = (ccds.y1 + ccds.y2 + ccds.y3 + ccds.y4) / 4.
    
    remain = fits_table('obstatus/mosaic-tiles_remaining.fits')
    print('Remaining tiles:', len(remain))
    rids = set(remain.tileid)
    I, = np.nonzero([tid in rids for tid in tiles.tileid])
    remain = tiles[I]
    print(len(remain), 'remaining tiles matched')
    remain.cut(np.argsort(remain.depth_90))
    remain.cut(remain.depth_90 < 22.5)
    print(len(remain), 'are below 90th-pct depth before including recent data')

    # Assume we get nominal depth for recently-observed tiles...
    recent = fits_table('obstatus/mosaic-tiles_obstatus.fits')
    recent.cut([d > recent_date for d in recent.z_date])
    print(len(recent), 'tiles observed since', recent_date)
    recent.depth = np.zeros(len(recent)) + 22.2
    recent_depth = depth_map_for_tiles(recent)

    print(len(set(remain.tileid).intersection(set(recent.tileid))), 'tiles in common between "remain" and "recent"')

    # 
    iv1 = 1./(10.**((depth - 22.5) / -2.5))**2
    iv3 = 1./(10.**((recent_depth - 22.5) / -2.5))**2
    assert(np.all(iv1 >= 0))
    assert(np.all(iv3 >= 0))
    # Update remaining tiles with depth from recent measurements...
    recdepth = -2.5 * (np.log10(1./np.sqrt(iv1 + iv3)) - 9.)
    dd = measure_map_at_tiles(remain, wcs, recdepth)
    print('Assuming nominal depth for recent tiles,',
          np.sum(dd[:,0] < 22.5), 'remaining tiles are shallow')
    # Plug in that new depth, & cut
    remain.depth_90 = dd[:,0]
    remain.cut(np.argsort(remain.depth_90))
    remain.cut(remain.depth_90 < 22.5)
    print('Cut to', len(remain), 'remaining tiles')

    print(np.sum(remain.bstarv < 10), 'remaining tiles contain bright stars')
    
    dtodo = tiles[(tiles.in_desi == 1) * (tiles.get('pass') <= 3) *
                  (tiles.dec > 30) * (tiles.edge == False) *
                  (tiles.depth_90 < 22.5)]
    print(len(dtodo), 'non-edge tiles are under-depth')
    dtodo.cut(dtodo.bstarv > 10.)
    print(len(dtodo), 'tiles are under-depth and do not contain bright stars')
    #dtodo = fits_table('tiles-todo-dstn.fits')
    #print(len(dtodo), 'dstn todo tiles')
    #dtodo.cut(dtodo.depth_90 < 22.5)
    print(len(dtodo), 'dstn todo tiles under depth')
    dtodo.cut(np.argsort(dtodo.depth_90))

    print(len(set(remain.tileid).intersection(set(dtodo.tileid))), 'tiles are in common, dtodo & remain')

    print(len(set(dtodo.tileid).intersection(set(recent.tileid))), 'tiles in common between "dtodo" and "recent"')
    
    # Assume we get nominal depth for David's tiles...
    rr = remain.copy()
    rr.depth = np.zeros(len(rr)) + 22.2
    remain_depth = depth_map_for_tiles(rr)
    iv2 = 1./(10.**((remain_depth - 22.5) / -2.5))**2
    assert(np.all(iv2 >= 0))
    rdepth = -2.5 * (np.log10(1./np.sqrt(iv1 + iv2 + iv3)) - 9.)
    dd = measure_map_at_tiles(dtodo, wcs, rdepth)
    print(np.sum(dd[:,0] < 22.5), 'dtodo tiles are shallow, assuming nominal depth for remaining+recent tiles')
    print(np.sum(dd[:,0] < 22.45), 'dtodo tiles are < 22.45, assuming nominal depth for remaining+recent tiles')
    print(np.sum(dd[:,0] < 22.4), 'dtodo tiles are < 22.4, assuming nominal depth for remaining+recent tiles')

    # Plug in that new depth
    dtodo.depth_90 = dd[:,0]
    dtodo.cut(np.argsort(dtodo.depth_90))
    dtodo.cut(dtodo.depth_90 < 22.5)
    
    print(len(set(dtodo.tileid).intersection(set(recent.tileid))), 'tiles in common between "dtodo" (after cutting on updated depth) and "recent"')

    print(len(set(dtodo.tileid).intersection(set(remain.tileid))), 'tiles in common between "dtodo" (after cutting on updated depth) and "remain"')

    print(len(set(remain.tileid).intersection(set(dtodo.tileid))), 'tiles are in common, dtodo & remain')

    cm = 'bone'
    
    ima = dict(interpolation='nearest', origin='lower',
               vmin=22.0, vmax=23, cmap=cm)
    ima2 = dict(interpolation='nearest', origin='lower',
                vmin=21.5, vmax=22.5, cmap=cm)
    tilestyle = dict(color='r', linestyle='-', linewidth=3)

    plt.subplots_adjust(left=0.01, right=0.99, bottom=0.01, top=0.9)
    import matplotlib.gridspec as gridspec

    plt.clf()
    plt.imshow(recdepth, **ima)
    plt.xticks([]); plt.yticks([])
    plt.colorbar()
    plt.title('Current+recent depth')
    ps.savefig()

    plt.clf()
    plt.imshow(recdepth, **ima)
    plt.xticks([]); plt.yticks([])
    plt.colorbar()
    ax = plt.axis()
    plt.plot(remain.x1, remain.y1, 'r.')
    plt.axis(ax)
    plt.title('Current+recent depth; %i remaining shallow tiles' % len(remain))
    ps.savefig()

    plt.clf()
    plt.imshow(rdepth, **ima)
    plt.xticks([]); plt.yticks([])
    plt.colorbar()
    ax = plt.axis()
    plt.plot(dtodo.x1, dtodo.y1, 'r.')
    plt.axis(ax)
    plt.title('Current+recent+remaining depth; %i dstn shallow tiles' % len(dtodo))
    ps.savefig()
    
    ### HACK
    plot_tiles = dtodo
    #plot_tiles = dtodo[np.argsort(dd[:,0])[:10]]
    #plot_tiles = remain

    iv_total = iv1 + iv2 + iv3

    dotiles = []
    
    #for ii,t in enumerate(plot_tiles):
    while True:

        ii = np.argmin(dtodo.depth_90)
        t = dtodo[ii]
        print('Worst tile:', t.tileid, 'at depth', t.depth_90)
        
        if t.depth_90 > 22.5:
            print('Tile', t.tileid, 'is already at depth')
            break
            #continue

        margin = 30
        bstr = ''
        if t.bstarv < 20:
            bstr = ', bright star: %.1f' % t.bstarv

        if False:
            plt.clf()
            spec = gridspec.GridSpec(ncols=3, nrows=2)
            fig.add_subplot(spec[0, :2])
            plt.imshow(depth, **ima)
            plt.xticks([]); plt.yticks([])
            plt.plot(t.x1, t.y1, 'r.')
    
            #plt.subplot(2,3,1)
            fig.add_subplot(spec[0, 2])
            plt.imshow(depth, **ima)
            plt.xticks([]); plt.yticks([])
            plt.colorbar()
            plt.plot(np.vstack((t.x1, t.x2, t.x3, t.x4, t.x1)),
                     np.vstack((t.y1, t.y2, t.y3, t.y4, t.y1)),
                     **tilestyle)
            plt.plot(brightstars.x, brightstars.y, 'ro')
            plt.axis([t.xmean-margin, t.xmean+margin,
                      t.ymean-margin, t.ymean+margin])
            plt.title('Total depth')
    
            for p in [1,2,3]:
    
                fig.add_subplot(spec[1, p-1])
                #plt.subplot(2,2,p+1)
                plt.imshow(depth_p[p], **ima2)
                plt.colorbar()
                plt.xticks([]); plt.yticks([])
                margin = 30
                plt.plot(np.vstack((t.x1, t.x2, t.x3, t.x4, t.x1)),
                         np.vstack((t.y1, t.y2, t.y3, t.y4, t.y1)),
                         **tilestyle)
    
                I = np.flatnonzero((np.abs(ccds.xmean - t.xmean) < margin) *
                                   (np.abs(ccds.ymean - t.ymean) < margin) *
                                   (ccds.tilepass == p))
                expnums = np.unique(ccds.expnum[I])
                for expnum in expnums:
                    J = I[ccds.expnum[I] == expnum]
    
                    plt.text(np.mean(ccds.xmean[J]), np.mean(ccds.ymean[J]),
                             '%.2f' % np.median(ccds.depth[J]),
                             color='k', ha='center', va='center',
                             bbox=dict(edgecolor='none', facecolor='white'))
                
                plt.axis([t.xmean-margin, t.xmean+margin,
                          t.ymean-margin, t.ymean+margin])
                plt.title('Pass %i depth' % p)
    
            cstr = ''
            if t.tileid in set(remain.tileid):
                cstr = ', in Remaining-tiles'
            plt.suptitle('Tile %i, Pass %i, RA,Dec %.3f,%.3f, depths %.2f / %.2f / %.2f, Z_expnum %i%s%s' %
                      (t.tileid, t.get('pass'), t.ra, t.dec, t.depth_90, t.depth_95, t.depth_98, t.z_expnum, bstr, cstr))
            ps.savefig()


        print('Added tile', t.tileid)
        thistile = dtodo[np.array([ii])]
        thistile.depth = np.zeros(len(thistile)) + 22.2
        tdepth = depth_map_for_tiles(thistile)
        iv = 1./(10.**((tdepth - 22.5) / -2.5))**2
        iv_total += iv
        rdepth = -2.5 * (np.log10(1./np.sqrt(iv_total)) - 9.)
        dd = measure_map_at_tiles(dtodo, wcs, rdepth)
        print('Before:', np.sum(dtodo.depth_90 < 22.5), 'tiles are under-depth; this tile', t.depth_90)
        print('After :', np.sum(dd[:,0] < 22.5), 'tiles are under-depth; this tile', dd[ii,0])
        dtodo.depth_90 = dd[:,0]

        depth = rdepth
        
        dotiles.append(ii)

    inorder = dtodo[np.array(dotiles)]
    inorder.writeto('tiles-dstn-inorder.fits')

    plt.clf()
    plt.imshow(depth, **ima)
    plt.xticks([]); plt.yticks([])
    plt.colorbar()
    ax = plt.axis()
    plt.plot(inorder.x1, inorder.y1, 'r.')
    plt.axis(ax)
    plt.title('%i dstn shallow tiles in order' % len(inorder))
    ps.savefig()


    # If I *don't* assume that all David's *remain* tiles get taken, what
    # is the ordered list?
    # Have to re-cut the 'dtodo' list, because we previously cut after
    # assuming we got the "remain" tiles.
    dtodo = tiles[(tiles.in_desi == 1) * (tiles.get('pass') <= 3) *
                  (tiles.dec > 30) * (tiles.edge == False) *
                  (tiles.depth_90 < 22.5) * (tiles.bstarv > 10.)]
    
    iv_total = iv1 + iv3

    rdepth = -2.5 * (np.log10(1./np.sqrt(iv_total)) - 9.)
    dd = measure_map_at_tiles(dtodo, wcs, rdepth)
    dtodo.depth_90 = dd[:,0]
    print('Without assuming the remaining tiles are taken:', np.sum(dtodo.depth_90 < 22.5), 'tiles are under-depth')
    
    dotiles = []
    while True:
        ii = np.argmin(dtodo.depth_90)
        t = dtodo[ii]
        print('Worst tile:', t.tileid, 'at depth', t.depth_90)
        
        if t.depth_90 > 22.5:
            print('Tile', t.tileid, 'is already at depth')
            break

        thistile = dtodo[np.array([ii])]
        thistile.depth = np.zeros(len(thistile)) + 22.2
        tdepth = depth_map_for_tiles(thistile)
        iv = 1./(10.**((tdepth - 22.5) / -2.5))**2
        iv_total += iv
        rdepth = -2.5 * (np.log10(1./np.sqrt(iv_total)) - 9.)
        dd = measure_map_at_tiles(dtodo, wcs, rdepth)
        print('Before:', np.sum(dtodo.depth_90 < 22.5), 'tiles are under-depth; this tile', t.depth_90)
        print('After :', np.sum(dd[:,0] < 22.5), 'tiles are under-depth; this tile', dd[ii,0])
        dtodo.depth_90 = dd[:,0]

        dotiles.append(ii)

    inorder = dtodo[np.array(dotiles)]
    inorder.writeto('tiles-dstn-inorder-solo.fits')
    print(len(inorder), 'tiles without assuming "remaining" tiles get taken')
    print(len(set(inorder.tileid).intersection(set(remain.tileid))), 'tiles in common with "remain"')

    plt.clf()
    plt.imshow(rdepth, **ima)
    plt.xticks([]); plt.yticks([])
    plt.colorbar()
    ax = plt.axis()
    plt.plot(inorder.x1, inorder.y1, 'r.')
    plt.axis(ax)
    plt.title('%i dstn shallow tiles in order' % len(inorder))
    ps.savefig()

    
if __name__ == '__main__':
    from_ccds()
    #update_tiles()
    djs_update()
    plot_tiles()
    #when_missing()
    #tiles_todo()
    #needed_tiles()
    #main()

    # tiles = fits_table('mosaic-tiles_remaining.fits')
    # ps = PlotSequence('remaining')
    # plot_tiles(tiles, ps)
    
