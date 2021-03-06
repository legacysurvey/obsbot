from __future__ import print_function
import matplotlib
matplotlib.use('Agg')
from astrometry.blind.plotstuff import *
from astrometry.util.util import *
from astrometry.util.fits import *
import pylab as plt
from collections import Counter
from scipy.ndimage.filters import minimum_filter, median_filter
from glob import glob
from astrometry.libkd.spherematch import match_radec
from astrometry.util.starutil_numpy import radectolb

# Last date included in the CCDs files
recent_date = '2018-02-11'

target_depths = dict(g=24.0, r=23.4, z=22.5)

compress = '[compress R; qz -1e-4]'

def mosaic_wcs(ra, dec, pixbin=1.):
    # This is pretty close to the outline of the four Mosaic chips.
    W = H = (4096 * 2 + 100) / pixbin
    cd = pixbin * 0.262 / 3600.
    tan = Tan(ra, dec, W/2., H/2., cd, 0., 0., cd,
              float(W), float(H))
    return tan

def draw_grid(wcs, ra_gridlines=None, dec_gridlines=None,
              ra_labels=None, dec_labels=None,
              ra_gridlines_decs=None, dec_gridlines_ras=None,
              ra_labels_dec=None, dec_labels_ra=None):
    if ra_gridlines_decs is None:
        ra_gridlines_decs = dec_gridlines
    if dec_gridlines_ras is None:
        ra_gridlines_decs = ra_gridlines
    ax = plt.axis()
    for d in dec_gridlines:
        rr = dec_gridlines_ras
        dd = np.zeros_like(rr) + d
        ok,xx,yy = wcs.radec2pixelxy(rr, dd)
        plt.plot(xx, yy, 'k-', alpha=0.1)
    for r in ra_gridlines:
        dd = ra_gridlines_decs
        rr = np.zeros_like(dd) + r
        ok,xx,yy = wcs.radec2pixelxy(rr, dd)
        plt.plot(xx, yy, 'k-', alpha=0.1)
    if ra_labels_dec is None:
        ra_labels_dec = dec_gridlines[0]
    if dec_labels_ra is None:
        dec_labels_ra = ra_gridlines[0]
    ok,xx,yy = wcs.radec2pixelxy(ra_labels, ra_labels_dec)
    plt.xticks(xx-0.5, ra_labels)
    ok,xx,yy = wcs.radec2pixelxy(dec_labels_ra, dec_labels)
    plt.yticks((yy-0.5), dec_labels)
    plt.axis(ax)

def get_grid_kwargs(mosaic, ngc):
    drawkwargs = {}
    if mosaic:
        dec_gridlines = list(range(30, 90, 10))
        dec_gridlines_ras = np.arange(0, 360, 1)
        ra_gridlines = range(0, 360, 30)
        ra_gridlines_decs = np.arange(20, 90, 1.)
    else:
        if ngc:
            dec_gridlines = list(range(-10, 31, 10))
            dec_gridlines_ras = np.arange(100, 280, 1)
            ra_gridlines = range(120, 271, 30)
            ra_gridlines_decs = np.arange(-10, 36, 1.)
            ra_labels = ra_gridlines
            dec_labels = dec_gridlines
        else:
            dec_gridlines = list(range(-20, 40, 10))
            dec_gridlines_ras = np.arange(0, 360, 1)
            #dec_gridlines_ras = np.arange(-80, 280, 1)
            ra_gridlines = range(0, 360, 30)
            ra_gridlines_decs = np.arange(-20, 41, 1.)
            ra_labels = ra_gridlines
            dec_labels = dec_gridlines
        drawkwargs.update(
            ra_gridlines_decs=ra_gridlines_decs,
            dec_gridlines_ras=dec_gridlines_ras,
            )
    
    drawkwargs.update(ra_gridlines=ra_gridlines, dec_gridlines=dec_gridlines,
                      ra_labels=ra_labels, dec_labels=dec_labels)
    return drawkwargs


def from_ccds(mp, mosaic=False, ngc=True,
              prefix_pat=None):
    from legacypipe.survey import LegacySurveyData

    drawkwargs = get_grid_kwargs(mosaic, ngc)

    if mosaic:
        obsfn = 'obstatus/mosaic-tiles_obstatus.fits'
    else:
        obsfn = 'obstatus/decam-tiles_obstatus.fits'

    survey = LegacySurveyData()
    T = survey.get_annotated_ccds()
    print('Read', len(T), 'CCDs from annotated table')
    T.cut(T.ccd_cuts == 0)
    print('Cut to', len(T), 'with ccd_cuts == 0')

    if ngc:
        T.l, T.b = radectolb(T.ra, T.dec)
        T.cut(T.b > 0)
        print('Cut to', len(T), 'in NGC')
    
    print('Tile passes from CCDs table:')
    print(Counter(T.tilepass).most_common())

    if mosaic:
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

    from camera import nominal_cal

    bands = ['g','r','z']
    
    print('Computing depth & transparency')
    nom = nominal_cal

    T.band = np.array([f.strip() for f in T.filter])
    
    # ugrizY
    band_index = dict(g=1, r=2, z=4)
    iband = np.array([band_index[b] for b in T.band])
    T.depth = T.galdepth - T.decam_extinction[np.arange(len(T)), iband]

    T.seeing = T.fwhm * T.pixscale_mean

    if False:
        #kx = 0.10
        #zp0 = 26.4
        #T.transparency = 10.**(-0.4 * (zp0 - T.ccdzpt - kx*(T.airmass-1.)))
        kx = dict([(b, nom.fiducial_exptime(b).k_co) for b in bands])
        zp0 = dict([(b, nom.zeropoint(b)) for b in bands])
        T.transparency = np.array([10.**(-0.4 * (zp0[b] - zpt - kx[b]*(airmass-1.)))
                                   for b,zpt,airmass in zip(T.band, T.ccdzpt, T.airmass)])
    
    # plt.clf()
    # plt.scatter(T.ra, T.dec, c=T.transparency, s=3)
    # plt.colorbar()
    # plt.savefig('transp.png')
    # 
    # for band in bands:
    #     for passnum in [1,2,3]:
    #         plt.clf()
    #         I = np.flatnonzero((T.tilepass == passnum) * (T.band == band))
    #         plt.scatter(T.ra[I], T.dec[I], c=T.transparency[I], s=3)
    #         plt.colorbar()
    #         plt.savefig('transp-%s-p%i.png' % (band, passnum))
    # 
    #         plt.clf()
    #         plt.scatter(T.ra[I], T.dec[I], c=T.airmass[I], s=3)
    #         plt.colorbar()
    #         plt.savefig('airmass-%s-p%i.png' % (band, passnum))

    print('Producing depth maps...')
    if prefix_pat is None:
        if ngc:
            prefix_pat = 'depth-ngc-%s-p%i'
        else:
            prefix_pat = 'maps/depth-%s-p%i'

    args = []
    for band in bands:
        target = target_depths[band]
        for passnum in [0,1,2,3]:
            prefix = prefix_pat % (band, passnum)
            fn = '%s.fits.fz' % prefix
            if os.path.exists(fn):
                print('Depth map already exists:', fn)
            else:
                I = np.flatnonzero((T.tilepass == passnum) *
                                   (T.depth > 10) *
                                   (T.band == band))
                print(len(I), 'CCDs for', band, 'pass', passnum)
                args.append((fn, survey, T[I], mosaic, ngc, target))
    print('mp.map depths for', len(args), 'maps...')
    mp.map(write_depth_map, args)
    print('Created maps.  Making plots...')

    assert((not mosaic) and (not ngc))
    #wcsfn = 'cea.wcs'
    wcsfn = 'cea-flip.wcs'
    wcs = anwcs(wcsfn)

    for band in bands:
        totaldepth = 0.
        seeings = []
        transps = []

        target = target_depths[band]

        for passnum in [1,2,3, 0,  9]:
            prefix = prefix_pat % (band, passnum)
                
            fn = '%s.fits.fz' % prefix
            if passnum != 9:
                tt = '%s band, Pass %i' % (band, passnum)
                print('Reading', fn)
                depthmap = fitsio.read(fn)
                iv = 1./(10.**((depthmap - 22.5) / -2.5))**2
                totaldepth = totaldepth + iv

                d = -2.5 * (np.log10(1./np.sqrt(totaldepth)) - 9.)
                plt.clf()
                plt.subplot(1,2,1)
                plt.hist(depthmap.ravel(), 50, log=True)
                plt.hist(d.ravel(), 50, log=True, range=(0,26), histtype='step', color='r')
                plt.subplot(1,2,2)
                plt.hist(iv.ravel(), 50, log=True, histtype='step', color='r')
                plt.hist(totaldepth.ravel(), 50, log=True, histtype='step', color='k')
                plt.savefig('depth-hist-%s-%i.png' % (band, passnum))

            else:
                depthmap = -2.5 * (np.log10(1./np.sqrt(totaldepth)) - 9.)
                tt = '%s band, All Passes' % band
                if os.path.exists(fn):
                    print('Already exists:', fn)
                else:
                    hdr = fitsio.read_header(wcsfn)
                    fitsio.write(fn + compress, depthmap, clobber=True, header=hdr)
                    print('Wrote', fn)

            fn = prefix + '-1.png'
            depth_map_plot(depthmap, wcs, fn, obsfn, target, drawkwargs,
                           tt='Depth: %s' % tt, redblue=False)
            print('Wrote', fn)

            fn = prefix + '-2.png'
            depth_map_plot(depthmap, wcs, fn, obsfn, target, drawkwargs,
                           tt='Depth: %s' % tt, redblue=True)
            print('Wrote', fn)


            # if passnum != 9:
            #     seeing,transp,wcs = best_seeing_map_for_ccds(survey, T[I])
            #     seeings.append(seeing)
            #     transps.append(transp)
            # else:
            #     seeing = seeings[0]
            #     transp = transps[0]
            #     for s,t in zip(seeings, transps):
            #         # Handle zeros...
            #         sI = (seeing == 0)
            #         seeing[sI] = s[sI]
            #         sI = (s != 0)
            #         seeing[sI] = np.minimum(seeing[sI], s[sI])
            #         transp = np.maximum(transp, t)
            # 
            # lo,hi = 0.5, 2.0
            # if passnum == 9:
            #     plt.clf()
            #     plt.imshow(median_filter(seeing, size=3),
            #                interpolation='nearest', origin='lower',
            #                vmin=lo, vmax=hi, cmap=cm)
            #     plt.xticks([]); plt.yticks([])
            #     plt.title('Best seeing (unfiltered): %s' % tt)
            #     plt.colorbar()
            #     plt.savefig('depth-%s-p%i-17.png' % (band, passnum))
            # 
            # plt.clf()
            # plt.imshow(median_filter(seeing, size=3),
            #            interpolation='nearest', origin='lower',
            #            vmin=lo, vmax=hi, cmap=cm)
            # plt.xticks([]); plt.yticks([])
            # plt.title('Best seeing: %s' % tt)
            # plt.colorbar()
            # plt.savefig('depth-%s-p%i-15.png' % (band, passnum))
            # 
            # plt.clf()
            # lo,hi = 0.5, 1.1
            # plt.imshow(transp, interpolation='nearest', origin='lower',
            #            vmin=lo, vmax=hi, cmap=cm)
            # plt.xticks([]); plt.yticks([])
            # plt.title('Best transparency: %s' % tt)
            # plt.colorbar()
            # plt.savefig('depth-%s-p%i-16.png' % (band, passnum))
            #fitsio.write('seeing-%s-p%i.fits' % (band, passnum), seeing, clobber=True)
            #fitsio.write('transp-%s-p%i.fits' % (band, passnum), transp, clobber=True)
            
    
            # plt.clf()
            # lo,hi = target-1, target+1
            # hh = depthmap.ravel()
            # hh = hh[hh != 0]
            # plt.hist(np.clip(hh, lo, hi), 50, range=(lo,hi))
            # plt.xlim(lo, hi)
            # plt.title('Depths: %s' % tt)
            # plt.yticks([])
            # plt.ylabel('sky area')
            # plt.savefig('depth-%s-p%i-10.png' % (band, passnum))
            # 
            # if len(I):
            #     plt.clf()
            #     plt.hist(np.clip(T.depth[I], lo+0.001, hi-0.001), 50,
            #              range=(lo,hi))
            #     plt.xlim(lo, hi)
            #     plt.title('Depths: %s' % tt)
            #     plt.ylabel('Number of CCDs')
            #     plt.savefig('depth-%s-p%i-11.png' % (band, passnum))


def depth_map_plot(depthmap, wcs, outfn, tilefn, target, drawkwargs,
                   tt=None, redblue=False):
    if redblue:
        cm = matplotlib.cm.RdBu
    else:
        cm = matplotlib.cm.viridis
    cm = cmap_discretize(cm, 10)
    cm.set_bad('0.9')

    tiles = fits_table(tilefn)
    tiles.cut(tiles.get('pass') == 1)
    tiles.cut(tiles.in_desi == 1)
    ok,tiles.x,tiles.y = wcs.radec2pixelxy(tiles.ra, tiles.dec)
    tiles.x -= 1
    tiles.y -= 1
    tiles.cut(ok)

    depthmap[depthmap < 1] = np.nan
    depthmap[np.logical_not(np.isfinite(depthmap))] = np.nan    

    H,W = wcs.shape
    imkwa = dict(interpolation='nearest', origin='lower',
                 extent=[0,W,0,H])
    lo,hi = target-1, target+1
            
    plt.figure(1)
    plt.clf()
    plt.plot(tiles.x, tiles.y, 'k.', alpha=0.1)
    plt.imshow(depthmap[::4,::4],
               vmin=lo, vmax=hi, cmap=cm, **imkwa)
    plt.xticks([]); plt.yticks([])
    draw_grid(wcs, **drawkwargs)
    cb = plt.colorbar(orientation='horizontal')
    cb.set_label('Depth')
    if tt is not None:
        plt.title(tt)
    plt.savefig(outfn)

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
        cdict[key] = [ (indices[i], colors_rgba[i-1,ki], colors_rgba[i,ki]) for i in range(N+1) ]
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

        
def write_depth_map(X):

    fn = X[0]
    args = X[1:]
    depthmap,wcs = depth_map_for_ccds(*args)
    print('Writing', fn)
    fitsio.write(fn + compress, depthmap, clobber=True)
    print('Wrote', fn)

def depth_map_for_ccds(survey, ccds, mosaic, ngc, targetdepth):
    if mosaic:
        W,H = 4800,3200
        plot = Plotstuff(size=(W,H), outformat='png')
        zoom = 2.7
        args = (195., 65., zoom, W, H, True)
        wcs = anwcs_create_hammer_aitoff(*args)
        plot.wcs = anwcs_create_hammer_aitoff(*args)
    else:
        if ngc:
            W,H = 8000,2000
            plot = Plotstuff(size=(W,H), outformat='png')

            # refra,refdec = 190., 12.
            # cd = 180. / W
            # wcs = anwcs_new_tan(Tan(refra, refdec, W/2.+0.5, H/2.+0.5,

            # zoom = 2.5
            # args = (190., 12., zoom, W, H, True)
            # wcs = anwcs_create_hammer_aitoff(*args)
            # plot.wcs = anwcs_create_hammer_aitoff(*args)

            pixscale = 180. / W
            ra,dec = 190., 12.
            refx = W/2. + 0.5
            refy = (H/2. + 0.5 - dec / -pixscale)
            args = (ra, 0., refx, refy, pixscale, W, H, True)
            wcs = anwcs_create_cea_wcs(*args)
            plot.wcs = anwcs_create_cea_wcs(*args)

            anwcs_write(wcs, 'ngc-cea.wcs')

        else:
            # W,H = 15000,2500
            # plot = Plotstuff(size=(W,H), outformat='png')
            # zoom = 1.2
            # args = (100., 5., zoom, W, H, True)
            # wcs = anwcs_create_hammer_aitoff(*args)
            # plot.wcs = anwcs_create_hammer_aitoff(*args)
            # W,H = 8000,2000
            # plot = Plotstuff(size=(W,H), outformat='png')

            # refra,refdec = 190., 12.
            # cd = 180. / W
            # wcs = anwcs_new_tan(Tan(refra, refdec, W/2.+0.5, H/2.+0.5,
            # zoom = 2.5
            # args = (190., 12., zoom, W, H, True)
            # wcs = anwcs_create_hammer_aitoff(*args)
            # plot.wcs = anwcs_create_hammer_aitoff(*args)

            #W,H = 16000,2500
            W,H = 32000,5000
            #W,H = 64000,10000
            plot = Plotstuff(size=(W,H), outformat='png')
            pixscale = 340. / W
            ra,dec = 110., 8.
            refx = W/2. + 0.5
            refy = (H/2. + 0.5 - dec / -pixscale)
            args = (ra, 0., refx, refy, pixscale, W, H, True)
            wcs = anwcs_create_cea_wcs(*args)
            plot.wcs = anwcs_create_cea_wcs(*args)

            wcsfn = 'cea.wcs'
            if os.path.exists(wcsfn):
                print('WCS file exists:', wcsfn)
            else:
                anwcs_write(wcs, wcsfn)

            # ?
            refy2 = (H/2. + 0.5 - dec / pixscale)
            args2 = (ra, 0., refx, refy2, pixscale, W, H, False)
            wcs2 = anwcs_create_cea_wcs(*args2)
            wcsfn = 'cea-flip.wcs'
            if os.path.exists(wcsfn):
                print('WCS file exists:', wcsfn)
            else:
                anwcs_write(wcs2, wcsfn)
            
    plot.color = 'black'
    plot.plot('fill')
    plot.color = 'white'
    plot.outline.fill = True
    plot.outline.stepsize = 2000
    plot.op = CAIRO_OPERATOR_ADD
    
    # This gives a 1-mag margin on the target depth
    maxfrac = 0.16
    
    targetsig = 10.**((targetdepth - 22.5) / -2.5)
    targetiv = 1./targetsig**2
    
    for j,ccd in enumerate(ccds):
        if j and j % 1000 == 0:
            print('CCD', j, '/', len(ccds))
        depth = ccd.depth
        mwcs = survey.get_approx_wcs(ccd)

        detsig = 10.**((depth - 22.5) / -2.5)
        depthfrac = 1./detsig**2 / targetiv

        plot.outline.wcs = anwcs_new_tan(mwcs)
        plot.alpha = np.clip(maxfrac * depthfrac, 0., 1.)
        plot.apply_settings()
        plot.plot('outline')

    # print('Creating numpy image...')
    # img1 = plot.get_image_as_numpy(flip=True)
    # img1 = img1[:,:,0]
    # print('img1 range', img1.min(), img1.max())

    print('Getting numpy view...')
    img = plot.get_image_as_numpy_view()
    print('Img:', img.shape, 'dtype', img.dtype)
    print('Img range:', img.min(), img.max())

    print('Converting to depth...')
    fimg = np.empty((H,W), np.float32)
    with np.errstate(invalid='ignore', divide='ignore'):
        fimg[:,:] = -2.5 * (np.log10(1./np.sqrt(
            img[::-1,:,0] * (1./255. * targetiv/maxfrac))) - 9.)
    # note the "::-1" vertical flip here

    # print('Cutting numpy image...')
    # #img = img[:,:,3]
    # img = img[:,:,0]
    # # Format is argb32
    # #print('Cutting numpy image...')
    # #img = img[:,:,0]
    # print('Scaling image...')
    # # It's a uint8; scale and convert to float
    # img = img * 1./255.
    # #print('Image', img.dtype, img.shape)
    # img /= maxfrac
    # # Now it's in factor of detiv of target depth.
    # # back to sig -> depth
    # with np.errstate(invalid='ignore', divide='ignore'):
    #     img = -2.5 * (np.log10(1./np.sqrt(img * targetiv)) - 9.)
    fimg[np.logical_not(np.isfinite(fimg))] = 0.
    # FIXME -- flip?
    print('Depth range:', fimg.min(), fimg.max())
    return fimg, wcs


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


def measure_maps_at_tiles(tiles, wcs, maps, pcts=[10,5,2], get_polygon=None):
    from astrometry.util.miscutils import point_in_poly

    if get_polygon is None:
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

    measures = np.zeros((len(tiles), len(maps), len(pcts)), np.float32)
    for itile,t in enumerate(tiles):
        if get_polygon is None:
            poly = np.array([[t.x1,t.y1],[t.x2,t.y2],[t.x3,t.y3],[t.x4,t.y4]])
        else:
            poly = get_polygon(tiles[itile], wcs)
        # -1: FITS to numpy coords
        x0 = int(np.floor(min(poly[:,0]))) - 1
        y0 = int(np.floor(min(poly[:,1]))) - 1
        x1 = int(np.ceil (max(poly[:,0]))) - 1
        y1 = int(np.ceil (max(poly[:,1]))) - 1
        # x0 = int(np.floor(min(t.x1, t.x2, t.x3, t.x4))) - 1
        # y0 = int(np.floor(min(t.y1, t.y2, t.y3, t.y4))) - 1
        # x1 = int(np.ceil( max(t.x1, t.x2, t.x3, t.x4))) - 1
        # y1 = int(np.ceil( max(t.y1, t.y2, t.y3, t.y4))) - 1
        xx,yy = np.meshgrid(np.arange(x0, x1+1), np.arange(y0, y1+1))
        inpoly = point_in_poly(xx, yy, poly-1)

        for imap,themap in enumerate(maps):
            zoom = themap[y0:y1+1, x0:x1+1]
            p = np.percentile(zoom[inpoly], pcts)
            measures[itile, imap, :] = p

        # if itile < 5:
        #     plt.clf()
        #     plt.imshow(tiledepth * inpoly,
        #                interpolation='nearest', origin='lower')
        #     plt.plot(poly[:,0]-x0-1, poly[:,1]-y0-1, 'k-')
        #     plt.savefig('tile-poly-%02i.png' % itile)

        if itile % 1000 == 0:
            print(itile, 'of', len(tiles))
        
    return measures

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
    #from legacyzpts.psfzpt_cuts import *
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

def get_decam_outline_function(close=False):
    H,W = 4094,2046
    F=fitsio.FITS('decam.wcs.gz')
    H,W = 4094,2046
    refra,refdec = 180.,0.
    wcsmap = dict()
    for i in range(1, len(F)):
        hdr = F[i].read_header()
        stepsize = 1024
        wcs = wcs_pv2sip_hdr(hdr, stepsize=stepsize, H=H)
        #print(wcs.get_crval())
        wcs.set_crval((refra, refdec))
        wcs.set_height(H)
        ccdname = hdr['EXTNAME'].strip()
        wcsmap[ccdname] = wcs
    outline = [
        ('N29', W, 0),
        ('N29', 0, 0),
        ('N25', W, 0),
        ('N25', 0, 0),
        ('N20', W, 0),
        ('N20', 0, 0),
        ('N14', W, 0),
        ('N8' , 0, 0),
        ('N1' , W, 0),
        ('S1' , 0, 0),
        ('S8' , W, 0),
        ('S14', 0, 0),
        ('S20', W, 0),
        ('S20', 0, 0),
        ('S25', W, 0),
        ('S25', 0, 0),
        ('S29', W, 0),
        ('S29', 0, 0),
        ('S31', 0, H),
        ('S31', W, H),
        ('S28', 0, H),
        ('S28', W, H),
        ('S24', 0, H),
        ('S24', W, H),
        ('S19', 0, H),
        ('S13', W, H),
        ('S7' , 0, H),
        ('N7' , W, H),
        ('N13', 0, H),
        ('N19', W, H),
        ('N24', 0, H),
        ('N24', W, H),
        ('N28', 0, H),
        ('N28', W, H),
        ('N31', 0, H),
        ('N31', W, H),
    ]
    if close:
        outline.append(outline[0])
    ra,dec = [],[]
    for ccd,x,y in outline:
        r,d = wcsmap[ccd].pixelxy2radec(x, y)
        ra.append(r)
        dec.append(d)
    dra = np.array(ra) - refra
    ddec = np.array(dec) - refdec
    # get_polygon uses "dra" and "ddec" from above.
    def get_polygon(tile, wcs):
        cosdec = np.cos(np.deg2rad(tile.dec))
        ra,dec = tile.ra + dra / cosdec, tile.dec + ddec
        ok,x,y = wcs.radec2pixelxy(ra, dec)
        H,W = wcs.shape
        x = np.clip(x, 1, W)
        y = np.clip(y, 1, H)
        return np.vstack((x,y)).T
    return get_polygon

def update_decam_tiles(done_pattern='maps/depth-%s-p9.fits.fz',
                       todo_pattern='maps/depth-todo-%s.fits.fz',
                       outfn='tiles-depths.fits'):
    tiles = fits_table('obstatus/decam-tiles_obstatus.fits')
    print(len(tiles), 'tiles')

    Idesi = np.flatnonzero(tiles.in_desi == 1)
    print(len(Idesi), 'tiles in DESI')

    bands = ['g','r','z']
    #bands = ['g']

    get_polygon = get_decam_outline_function()
    wcs = anwcs('cea-flip.wcs')
    
    depths = []
    for band in bands:
        fn = done_pattern % band
        print('Reading', fn)
        depthmap = fitsio.read(fn)
        assert(depthmap.shape == wcs.shape)
        depths.append(depthmap)

    # Adding projected depth of to-do tiles
    #todos = []
    #totals = []
    for iband,band in enumerate(bands):
        fn = todo_pattern % band
        print('Reading', fn)
        depthmap = fitsio.read(fn)
        assert(depthmap.shape == wcs.shape)

        #todos.append(depthmap)

        iv1 = 1./(10.**((depthmap - 22.5) / -2.5))**2
        # Grab the already-done map for this band
        iv2 = 1./(10.**((depths[iband] - 22.5) / -2.5))**2
        depthmap = -2.5 * (np.log10(1./np.sqrt(iv1 + iv2)) - 9.)

        #totals.append(depthmap)
        
        depths.append(depthmap)


    Ides = np.flatnonzero((tiles.in_desi == 1) * (tiles.in_des == 1))
    print(len(Ides), 'tiles in DESI & DES')


    plt.figure(1, figsize=(13,8))
    plt.subplots_adjust(left=0.03, right=0.99, bottom=0.01, top=0.95)

    # for band in bands:
    #     for passnum in [0,1,2,3]:
    #         fn = 'maps/depth-%s-p%i.fits.fz' % (band, passnum)
    #         print('Reading', fn)
    #         depth = fitsio.read(fn)
    #         H,W = depth.shape
    #         plt.clf()
    #         plt.imshow(depth[::4,::4], interpolation='nearest', origin='lower',
    #                    extent=[0,W,0,H], cmap='viridis', vmin=23, vmax=26)
    #         plt.colorbar(orientation='horizontal')
    #         plt.savefig('depth-meas-p%i.png' % passnum)
    # 
    # # g done + todo
    # #depth = depths[-3]
    # depth = totals[0]
    # H,W = depth.shape
    # plt.clf()
    # plt.imshow(depth[::4,::4], interpolation='nearest', origin='lower',
    #            extent=[0,W,0,H], cmap='viridis', vmin=23, vmax=26)
    # poly = get_polygon(tiles[Ides[0]], wcs)
    # plt.plot(poly[:,0], poly[:,1], 'k-')
    # plt.colorbar(orientation='horizontal')
    # plt.savefig('depth-meas-total.png')
    # 
    # # g todo
    # depth = todos[0]
    # H,W = depth.shape
    # plt.clf()
    # plt.imshow(depth[::4,::4], interpolation='nearest', origin='lower',
    #            extent=[0,W,0,H], cmap='viridis', vmin=23, vmax=26)
    # poly = get_polygon(tiles[Ides[0]], wcs)
    # plt.plot(poly[:,0], poly[:,1], 'k-')
    # plt.colorbar(orientation='horizontal')
    # plt.savefig('depth-meas-todo.png')
    # 
    # # g done
    # depth = depths[0]
    # H,W = depth.shape
    # plt.clf()
    # plt.imshow(depth[::4,::4], interpolation='nearest', origin='lower',
    #            extent=[0,W,0,H], cmap='viridis', vmin=23, vmax=26)
    # poly = get_polygon(tiles[Ides[0]], wcs)
    # plt.plot(poly[:,0], poly[:,1], 'k-')
    # plt.colorbar(orientation='horizontal')
    # plt.savefig('depth-meas-done.png')
    # 
    # x = np.mean(poly[:,0])
    # y = np.mean(poly[:,1])
    # plt.axis([x-200, x+200, y-200, y+200])
    # plt.savefig('depth-meas-done-zoom.png')

    dd = measure_maps_at_tiles(tiles[Idesi], wcs, depths,
                               get_polygon=get_polygon)
    for iband,band in enumerate(bands):
        depth_90 = np.zeros(len(tiles), np.float32)
        depth_95 = np.zeros(len(tiles), np.float32)
        depth_98 = np.zeros(len(tiles), np.float32)
        depth_90[Idesi] = dd[:,iband,0]
        depth_95[Idesi] = dd[:,iband,1]
        depth_98[Idesi] = dd[:,iband,2]
        tiles.set('depth_%s_90' % (band), depth_90.copy())
        tiles.set('depth_%s_95' % (band), depth_95.copy())
        tiles.set('depth_%s_98' % (band), depth_98.copy())

        # Projected depth
        depth_90[Idesi] = dd[:,len(bands)+iband,0]
        depth_95[Idesi] = dd[:,len(bands)+iband,1]
        depth_98[Idesi] = dd[:,len(bands)+iband,2]
        tiles.set('proj_depth_%s_90' % (band), depth_90)
        tiles.set('proj_depth_%s_95' % (band), depth_95)
        tiles.set('proj_depth_%s_98' % (band), depth_98)
        
    tiles.writeto(outfn)

def decam_todo_maps(mp, obsfn=None, pattern='maps/depth-todo-%s.fits.fz',
                    use_depth=False):
    from legacypipe.survey import LegacySurveyData
    bands = ['g','r','z']
    if obsfn is None:
        obsfn = 'obstatus/decam-tiles_obstatus.fits'
    tiles = fits_table(obsfn)
    #tiles.cut(tiles.get('pass') == passnum)
    #tiles.cut(tiles.get('pass') == 1)
    tiles.cut(tiles.in_desi == 1)
    tiles.cut(tiles.in_des == 0)

    F=fitsio.FITS('decam.wcs.gz')
    H,W = 4094,2046
    nccd = len(F)-1
    C = fits_table()
    C.crpix1 = np.zeros(nccd, np.float32)
    C.crpix2 = np.zeros(nccd, np.float32)
    C.cd1_1 = np.zeros(nccd, np.float32)
    C.cd1_2 = np.zeros(nccd, np.float32)
    C.cd2_1 = np.zeros(nccd, np.float32)
    C.cd2_2 = np.zeros(nccd, np.float32)
    for i in range(nccd):
        hdr = F[i+1].read_header()
        C.crpix1[i] = hdr['CRPIX1']
        C.crpix2[i] = hdr['CRPIX2']
        C.cd1_1[i] = hdr['CD1_1']
        C.cd1_2[i] = hdr['CD1_2']
        C.cd2_1[i] = hdr['CD2_1']
        C.cd2_2[i] = hdr['CD2_2']
    print('Assuming', nccd, 'CCDs')

    survey = LegacySurveyData()

    args = []
    for band in bands:
        I = np.flatnonzero(tiles.get('%s_done' % band) == 0)
        print(len(I), 'to do in', band)
        ccds = fits_table()

        c = np.zeros((len(I),nccd), np.float32)
        for col in ['crpix1','crpix2','cd1_1','cd1_2','cd2_1','cd2_2']:
            c[:,:] = C.get(col)[np.newaxis,:]
            ccds.set(col, c.ravel().copy())
        c[:,:] = tiles.ra[I,np.newaxis]
        ccds.crval1 = c.ravel().copy()
        c[:,:] = tiles.dec[I,np.newaxis]
        ccds.crval2 = c.ravel().copy()
        ccds.width  = np.zeros(len(ccds), np.int32) + 2046
        ccds.height = np.zeros(len(ccds), np.int32) + 4094

        if use_depth:
            ccds.depth = tiles.get('%s_depth' % band)[I]
        else:
            target = target_depths[band]
            depth = target - 0.3
            ccds.depth = np.zeros(len(ccds), np.float32) + depth

        fn = pattern % band
        if os.path.exists(fn):
            print('Depth map already exists:', fn)
        else:
            args.append((fn, survey, ccds, False, False, target))
    mp.map(write_depth_map, args)


    # Plots

    drawkwargs = get_grid_kwargs(False, False)
    wcsfn = 'cea-flip.wcs'
    wcs = anwcs(wcsfn)

    for band in bands:
        fn = pattern % band
        depthmap = fitsio.read(fn)
        target = target_depths[band]
        fn = fn.replace('.fits.fz', '-1.png')
        depth_map_plot(depthmap, wcs, fn, obsfn, target, drawkwargs,
                       tt='Depth: %s' % band, redblue=False)
        print('Wrote', fn)
        # fn = fn.replace('-1.png', '-2.png')
        # depth_map_plot(depthmap, wcs, fn, obsfn, target, drawkwargs,
        #                tt='Depth: %s' % band, redblue=True)
        # print('Wrote', fn)



def des_zooms(tiles_fn='tiles-depths.fits',
              done_pattern='maps/depth-%s-p9.fits.fz',
              todo_pattern='maps/depth-todo-%s.fits.fz'):

    bands = ['g','r','z']

    T = fits_table(tiles_fn)
    T.cut(T.in_desi == 1)
    print(len(T), 'tiles in DESI')
    T.cut(np.lexsort((T.ra, -T.dec)))
    
    def wrap_ra(ra):
        return ra + (ra > 180)*-360
    
    T.ra_wrap = wrap_ra(T.ra)
    T.cut((T.ra_wrap > -60) * (T.ra_wrap < 60) * (T.dec > -20) * (T.dec < 10))

    expfn = 'p0-exposures.fits'
    if not os.path.exists(expfn):
        from legacypipe.survey import LegacySurveyData
        survey = LegacySurveyData()
        ccds = survey.get_annotated_ccds()
        ccds.cut(ccds.tilepass == 0)
        #ccds.cut(np.array([f.strip() == band for f in ccds.filter]))
        print(len(ccds), 'band', band, 'pass 0')
        #bores = set(zip(ccds.ra_bore, ccds.dec_bore))
        #print(len(bores), 'unique boresights')
        e,I = np.unique(ccds.expnum, return_index=True)
        exps = ccds[I]
        print(len(exps), 'unique exposures')
        e2 = fits_table()
        e2.ra = exps.ra_bore
        e2.dec = exps.dec_bore
        e2.filter = exps.filter
        e2.writeto(expfn)

    plt.figure(1, figsize=(13,8))
    plt.subplots_adjust(left=0.03, right=0.99, bottom=0.05, top=0.95)

    get_polygon = get_decam_outline_function(close=True)
    size = 200

    shallow_set = set()

    for band in bands:
        target = target_depths[band]

        exps = fits_table(expfn)
        exps.cut(np.array([f.strip() == band for f in exps.filter]))
        exps.ra_wrap = wrap_ra(exps.ra)
        exps.cut((exps.ra_wrap > -60) * (exps.ra_wrap < 60) *
                 (exps.dec > -20) * (exps.dec < 10))
    
        fn = todo_pattern % band
        print('Reading', fn)
        todo = fitsio.FITS(fn)[1]
        fn = done_pattern % band
        print('Reading', fn)
        depthmap = fitsio.FITS(fn)[1]

        wcs = anwcs('cea-flip.wcs')
        H,W = wcs.shape
    
        # todo = 1./(10.**((todo - 22.5) / -2.5))**2
        # depthmap = 1./(10.**((depthmap - 22.5) / -2.5))**2
        # # Grab the already-done map for this band
        # depthmap = -2.5 * (np.log10(1./np.sqrt(todo + depthmap)) - 9.)
        # del todo
    

        Ides = np.flatnonzero(T.in_des == 1)

        Inotdes = np.flatnonzero(T.in_des == 0)
    
        Idone = np.flatnonzero(T.get('%s_done' % band) == 1)
        Itodo = np.flatnonzero((T.get('%s_done' % band) == 0) * (T.in_des == 0))

        # Tiles to plot...
        # I = np.flatnonzero((T.in_des == 1) *
        #                    (T.get('proj_depth_%s_90' % band) < target))
        # print(len(I), 'in DES with projected depth in', band, 'too shallow')

        # I = np.flatnonzero((T.in_desi == 1) *
        #                    (T.get('proj_depth_%s_90' % band) < target))
        # print(len(I), 'in DESI with projected depth in', band, 'too shallow')

        # I = np.flatnonzero((T.in_desi == 1) *
        #                    np.logical_or((np.abs(T.dec - 5) < 0.5) * (T.ra_wrap > 0),
        #                                  (np.abs(T.dec - 2.8) < 0.5) * (T.ra_wrap < 0)))

        J,K,d = match_radec(T.ra[Ides], T.dec[Ides], T.ra[Inotdes], T.dec[Inotdes], 3.)
        keep = np.zeros(len(T), bool)
        keep[Ides[J]] = True
        #keep[Inotdes[K]] = True

        #I = np.flatnonzero((T.in_desi == 1) * (T.dec > -14))
        #I = np.flatnonzero(keep * (T.dec > -14))
        #I = I[np.argsort(T.get('proj_depth_%s_98' % band)[I])]

        print('Band', band, 'target', target)

        ok,x,y = wcs.radec2pixelxy(T.ra, T.dec)
        x = x - 1
        y = y - 1
        inbounds = (x >= size) * (y >= size) * (x < (W-size)) * (y < (H-size))

        I = np.flatnonzero(keep * (T.dec > -14) *
                           inbounds * (T.get('proj_depth_%s_98' % band) < (target+0.2)))

        print(len(I), 'shallow-ish')
        shallow_set.update(I)

        print('In-DESI: in-des tally:', Counter(T.in_des[I]))

        if True:
            wholedepth = depthmap.read()
            wholetodo = todo.read()
            depth_sum = -2.5 * (np.log10(1./np.sqrt(
                1./(10.**((wholedepth - 22.5) / -2.5))**2 +
                1./(10.**((wholetodo  - 22.5) / -2.5))**2
            )) - 9.)
            depth_sum[depth_sum < 1] = np.nan
            depth_sum[np.logical_not(np.isfinite(depth_sum))] = np.nan    
            plt.clf()
            plt.imshow(depth_sum, interpolation='nearest', origin='lower',
                       vmin=target-1, vmax=target+1, zorder=20, cmap='RdBu')
            plt.xticks([]); plt.yticks([])
            cb = plt.colorbar(orientation='horizontal')
            cb.set_label('%s depth (90%% coverage)' % band)
            
            ok,x,y = wcs.radec2pixelxy(T.ra[I], T.dec[I])
            plt.plot(x, y, 'm.', zorder=30)
    
            ok,x1,y1 = wcs.radec2pixelxy(300, -20)
            ok,x2,y2 = wcs.radec2pixelxy( 75,  10)
            plt.axis([min(x1,x2), max(x1,x2), min(y1,y2), max(y1,y2)])
    
            plt.savefig('des-%s.png' % band)
    
            # plt.clf()
            # plt.imshow(wholedepth, interpolation='nearest', origin='lower',
            #            vmin=target-1, vmax=target+1, zorder=20, cmap='RdBu')
            # plt.xticks([]); plt.yticks([])
            # cb = plt.colorbar(orientation='horizontal')
            # cb.set_label('%s depth (90%% coverage)' % band)
            # plt.savefig('des-%s-1.png' % band)
            # 
            # plt.clf()
            # plt.imshow(wholetodo, interpolation='nearest', origin='lower',
            #            vmin=target-1, vmax=target+1, zorder=20, cmap='RdBu')
            # plt.xticks([]); plt.yticks([])
            # cb = plt.colorbar(orientation='horizontal')
            # cb.set_label('%s depth (90%% coverage)' % band)
            # plt.savefig('des-%s-2.png' % band)
    
    
            continue

        # ok,x,y = wcs.radec2pixelxy(T.ra[I], T.dec[I])
        # plt.plot(x-1, y-1, 'r.')
        # plt.savefig('depth.png')
        # plt.clf()
        # plt.imshow(wholetodo, interpolation='nearest', origin='lower')
        # ok,x,y = wcs.radec2pixelxy(T.ra[I], T.dec[I])
        # plt.plot(x-1, y-1, 'r.')
        # plt.savefig('todo.png')



        k = 1
        for itile in I:
            ok,x,y = wcs.radec2pixelxy(T.ra[itile], T.dec[itile])
            x = int(x - 1)
            y = int(y - 1)
            if x < size or y < size or x >= W-size or y >= H-size:
                print('RA,Dec', T.ra[itile], T.dec[itile], 'too close to edge, skipping')
                continue
            x0 = x-size
            y0 = y-size
            x1 = x0+size*2+1
            y1 = y0+size*2+1
            depth_todo = todo[y0:y1, x0:x1]
            depth_done = depthmap[y0:y1, x0:x1]
    
            depth_sum = -2.5 * (np.log10(1./np.sqrt(
                1./(10.**((depth_todo - 22.5) / -2.5))**2 +
                1./(10.**((depth_done - 22.5) / -2.5))**2
                )) - 9.)
    
            ok,rr,dd = wcs.pixelxy2radec(np.array([x0,x0,x1,x1,x0]), np.array([y0,y1,y1,y0,y0]))
            
            plt.clf()
            #plt.imshow(depth_done, interpolation='nearest', origin='lower')
            #plt.subplot(1,2,2)
            #plt.imshow(depth_todo, interpolation='nearest', origin='lower')
    
            ax1 = plt.subplot2grid((2, 2), (0, 0))
            #plt.subplot(1,3,1)
            plt.plot(T.ra_wrap[Inotdes], T.dec[Inotdes], 'k.', alpha=0.1)
            plt.plot(T.ra_wrap[Idone], T.dec[Idone], 'k.', alpha=0.5)
            plt.plot(exps.ra_wrap, exps.dec, 'r.', alpha=0.5)
            #plt.plot(T.ra_wrap[Ides], T.dec[Ides], 'r.', alpha=0.25)
            plt.plot(T.ra_wrap[itile], T.dec[itile], 'mo', ms=20, mec='m',
                     mew=2, mfc='none')
            plt.xlabel('RA')
            plt.ylabel('Dec')
            plt.axis('scaled')
            plt.axis([60, -60, -20, 10])
    
            ax2 = plt.subplot2grid((2, 2), (1, 0))
            #plt.subplot(1,3,2)
            plt.plot(T.ra_wrap[Inotdes], T.dec[Inotdes], 'k.', alpha=0.5)
            plt.plot(exps.ra_wrap, exps.dec, 'r.', alpha=0.5)
            plt.plot(T.ra_wrap[itile], T.dec[itile], 'mo', ms=20, mec='m',
                     mfc='none')
            plt.plot(wrap_ra(rr), dd, 'k-')
            plt.xlabel('RA')
            plt.ylabel('Dec')
            plt.axis('scaled')
            plt.axis([T.ra_wrap[itile]+10, T.ra_wrap[itile]-10,
                      T.dec[itile]-5, T.dec[itile]+5])
    
            ax3 = plt.subplot2grid((2, 2), (0, 1), rowspan=2)
            #plt.subplot(1,3,3)
            plt.imshow(depth_sum, interpolation='nearest', origin='lower',
                       vmin=target-1, vmax=target+1, zorder=20, cmap='RdBu')
            plt.xticks([]); plt.yticks([])
            cb = plt.colorbar(orientation='horizontal')
            cb.set_label('%s depth (90%% coverage)' % band)
            ax = plt.axis()
            poly = get_polygon(T[itile], wcs)
            plt.plot(poly[:,0]-x0, poly[:,1]-y0, 'k-', zorder=30)
    
            I,J,d = match_radec(T.ra[itile], T.dec[itile], T.ra[Inotdes], T.dec[Inotdes], 5.)
            Inear = Inotdes[J]
            for i in Inear:
                poly = get_polygon(T[i], wcs)
                # Gray outlined nearby tiles?
                #plt.plot(poly[:,0]-x0, poly[:,1]-y0, '-', color='0.5', lw=2, zorder=25)
            
            # I,J,d = match_radec(T.ra[itile], T.dec[itile], T.ra[Idone], T.dec[Idone], 5.)
            # Inear = Idone[J]
            # for i in Inear:
            #     poly = get_polygon(T[i], wcs)
            #     plt.plot(poly[:,0]-x0, poly[:,1]-y0, '-', color='0.5', lw=2, zorder=25)
            # 
            # I,J,d = match_radec(T.ra[itile], T.dec[itile], T.ra[Itodo], T.dec[Itodo], 5.)
            # Inear = Itodo[J]
            # for i in Inear:
            #     poly = get_polygon(T[i], wcs)
            #     plt.plot(poly[:,0]-x0, poly[:,1]-y0, '-', color='b', lw=2, zorder=25)
            #     
            # I,J,d = match_radec(T.ra[itile], T.dec[itile], exps.ra_bore, exps.dec_bore, 5.)
            # # class duck(object):
            # #     pass
            # for i in J:
            #     #t = duck()
            #     #t.ra = exps[i].ra_bore
            #     #t.dec = exps[i].dec_bore
            #     poly = get_polygon(exps[i], wcs)
            #     plt.plot(poly[:,0]-x0, poly[:,1]-y0, '-', color='r', lw=2, zorder=25, alpha=0.25)
    
            plt.axis(ax)
    
            plt.suptitle('%s band tile (pass %i) at %.3f,%.3f: proj. depths %.2f / %.2f / %.2f' %
                         (band, T.get('pass')[itile], T.ra[itile], T.dec[itile],
                          T.get('proj_depth_%s_90' % band)[itile],
                          T.get('proj_depth_%s_95' % band)[itile],
                          T.get('proj_depth_%s_98' % band)[itile]))
            
            fn = 'des-%s-%02i.png' % (band, k)
            plt.savefig(fn)
            print('Wrote', fn)
            k += 1
            if k > 100:
                break



    print(len(shallow_set), '"shallow"')
    I = np.array(list(shallow_set))
    sh = T[I]
    sh.cut(np.argsort(sh.proj_depth_g_98))
    sh.writeto('tiles-des-tiling.fits')


def zsee_fig():

    def aitoff_labels(wcs):
        # image bounds
        declo = -25
        ok,x,y = wcs.radec2pixelxy(270., declo)
        ylo = int(y-1)
        ok,x,y = wcs.radec2pixelxy(270., +90.)
        yhi = int(y-1)
        ok,x,y = wcs.radec2pixelxy([90.001, 89.999], [0.,0.])
        xlo = min(x)
        xhi = max(x)
        print('x range', xlo,xhi)    
        
        # RA labels
        lt,lv = [],[]
        for ra in [89.99, 60,30,0,330,300,270,240,210,180,150,120,90]:
            decs = np.linspace(declo-10, declo+10, 25)
            #ok,x,y = wcs.radec2pixelxy(ra, declo)
            ok,x,y = wcs.radec2pixelxy(ra, decs)
            i = np.argmin(np.abs(y - ylo))
            x = x[i]
            lt.append('%i' % np.round(ra))
            lv.append(x)
        plt.xticks(lv, lt)
        plt.yticks([])
    
        for ra in [89.999, 60,30,0,330,300,270,240,210,180,150,120,90.001]:
            decs = np.linspace(declo-5, +90., 100)
            ok,x,y = wcs.radec2pixelxy(ra + np.zeros_like(decs), decs)
            plt.plot(x-1, y-1, 'k-', alpha=0.3)
    
        for dec in [-15, 0, 15, 30, 45, 60, 75]:
            ras = np.append([89.999], np.linspace(89.999, -269.999, 100))
            ok,x,y = wcs.radec2pixelxy(ras, np.zeros_like(ras)+dec)
            #print('ras', ras)
            #print('x', x)
            plt.plot(x-1, y-1, 'k-', alpha=0.3)
    
            # i = np.argmin(x)
            # xl = x[i]-1
            # yl = y[i]-1
            # plt.text(xl-10, yl, '%i' % dec, va='center', ha='right')
    
        for dec in [-15, 0, 15, 30, 45, 60, 75]:
            ra = 300
            ok,x,y = wcs.radec2pixelxy(ra, dec)
            if dec < 30:
                ha = 'left'
                dx = +5
            else:
                ha = 'right'
                dx = -5
            plt.text(x+dx, y-5, '$%+i$' % dec, va='top', ha=ha)
    
            
        plt.axis([xlo-3,xhi+3,ylo-3,yhi+1])
        plt.xlabel('RA (deg)')




    matplotlib.rc('text', usetex=True)
    matplotlib.rc('font', family='serif')
    plt.figure(1, figsize=(13,4))
    plt.subplots_adjust(left=0.03, right=0.99, bottom=0.08, top=0.95)

    
    
    #W,H = 2000,1000
    W,H = 1000,500
    zoom = 1.0
    wcs = anwcs_create_hammer_aitoff(270., 0., zoom, W, H, False)


    # Aaron's photometric offset plots
    cm = cmap_discretize(matplotlib.cm.RdBu, 21)
    #cm.set_bad('0.9')
    plt.subplots_adjust(left=0.03, right=0.99, bottom=0.11, top=0.94)

    # background = np.zeros((H,W,4), np.uint8)
    # xx,yy = np.meshgrid(np.arange(W), np.arange(H))
    # ok,rr,dd = wcs.pixelxy2radec(xx+1, yy+1)
    # background[yy[ok], xx[ok], :3] = np.uint8(0.95*255)
    # background[yy[ok], xx[ok], 3] = 255

    # Astrometry
    if False:
        Tra = fits_table('dr7_dr6_qa_combined-dra.table.fits')
        Tdec = fits_table('dr7_dr6_qa_combined-ddec.table.fits')
        Tboth = Tra.copy()
        Tboth.dmag = np.hypot(Tra.dmag, Tdec.dmag)
        for T,name in [(Tra, 'RA'), (Tdec, 'Dec'), (Tboth, 'Total')]:
            ok,x,y = wcs.radec2pixelxy(T.ra, T.dec)
            T.ix = np.round(x-1).astype(int)
            T.iy = np.round(y-1).astype(int)
            T.cut((T.ix >= 0) * (T.ix < W) * (T.iy >= 0) * (T.iy < H))
            print('Cut to', len(T), 'in-bounds')
            T.dmag[T.dmag == 0] = np.nan
            plt.clf()
            mx = 10.
            if name == 'Total':
                mn = 0.
                thiscm = cmap_discretize(matplotlib.cm.Blues, 11)
            else:
                mn = -mx
                thiscm = cm
            plt.scatter(T.ix, T.iy, c=T.dmag * 1000., s=0.5, vmin=mn, vmax=mx, cmap=thiscm)
            aitoff_labels(wcs)
            cb = plt.colorbar()
            cb.set_label('%s offsets (milli-arcsec)' % name)
            plt.title('DR6 + DR7 vs Gaia DR1')
            plt.savefig('astrometric-offset-%s.pdf' % name.lower())
        return

    # Photometry
    if False:
        for band in 'grz':
            T = fits_table('dr7_dr6_qa_combined-%s.table.fits' % band)
            ok,x,y = wcs.radec2pixelxy(T.ra, T.dec)
            T.ix = np.round(x-1).astype(int)
            T.iy = np.round(y-1).astype(int)
            T.cut((T.ix >= 0) * (T.ix < W) * (T.iy >= 0) * (T.iy < H))
            print('Cut to', len(T), 'in-bounds')
    
            T.dmag[T.dmag == 0] = np.nan
                    
            plt.clf()
            #plt.imshow(background, interpolation='nearest', origin='lower')
    
            mx = 10.
            if band in 'rz':
                mx = 8.
            
            plt.scatter(T.ix, T.iy, c=T.dmag * 1000., s=0.5, vmin=-mx, vmax=mx, cmap=cm)
            aitoff_labels(wcs)
            cb = plt.colorbar()
            cb.set_label('%s offsets (milli-mag)' % band)
            plt.title('DR6 + DR7 vs Pan-STARRS PS1')
            plt.savefig('photometric-offset-%s.pdf' % band)



    plt.subplots_adjust(left=0.03, right=0.99, bottom=0.08, top=0.95)

    # Galdepth
    if False:
        dr6 = fits_table('/data1/dr6/survey-bricks-dr6.fits.gz',
                         columns=['ra','dec','galdepth_g','galdepth_r',
                                  'galdepth_z'])
        print(len(dr6), 'DR6 bricks')
        # dr5 = fits_table('/data1/dr5/survey-bricks-dr5.fits.gz',
        #                  columns=['ra','dec','galdepth_g', 'galdepth_r',
        #                           'galdepth_z'])
        # print(len(dr5), 'DR5 bricks')
        dr7 = fits_table('/data1/dr7/survey-bricks-dr7.fits.gz',
                         columns=['ra','dec','galdepth_g', 'galdepth_r',
                                  'galdepth_z'])
        print(len(dr7), 'DR7 bricks')
    
        ok,x,y = wcs.radec2pixelxy(dr6.ra, dr6.dec)
        dr6.ix = np.round(x-1).astype(int)
        dr6.iy = np.round(y-1).astype(int)
        ok,x,y = wcs.radec2pixelxy(dr7.ra, dr7.dec)
        dr7.ix = np.round(x-1).astype(int)
        dr7.iy = np.round(y-1).astype(int)
    
        dr6.cut((dr6.ix >= 0) * (dr6.ix < W) * (dr6.iy >= 0) * (dr6.iy < H))
        print('Cut to', len(dr6), 'in-bounds')
        dr7.cut((dr7.ix >= 0) * (dr7.ix < W) * (dr7.iy >= 0) * (dr7.iy < H))
        print('Cut to', len(dr7), 'in-bounds')
        
        gdepth = np.zeros((H,W), np.float32)
        rdepth = np.zeros((H,W), np.float32)
        zdepth = np.zeros((H,W), np.float32)
        for x,y,g,r,z in zip(dr6.ix, dr6.iy, dr6.galdepth_g, dr6.galdepth_r, dr6.galdepth_z):
            gdepth[y,x] = max(gdepth[y,x], g)
            rdepth[y,x] = max(rdepth[y,x], r)
            zdepth[y,x] = max(zdepth[y,x], z)
        for x,y,g,r,z in zip(dr7.ix, dr7.iy, dr7.galdepth_g, dr7.galdepth_r, dr7.galdepth_z):
            gdepth[y,x] = max(gdepth[y,x], g)
            rdepth[y,x] = max(rdepth[y,x], r)
            zdepth[y,x] = max(zdepth[y,x], z)
    
        depths = dict(g=gdepth, r=rdepth, z=zdepth)

        dr7 = fits_table('/data1/dr7/ccds-annotated-dr7.fits.gz',
                         columns=['ra','dec','filter','ccd_cuts','galdepth'])
        dr7.cut(dr7.ccd_cuts == 0)
        ok,x,y = wcs.radec2pixelxy(dr7.ra, dr7.dec)
        dr7.ix = np.round(x-1).astype(int)
        dr7.iy = np.round(y-1).astype(int)
        for band in 'grz':
            depth = depths[band]
            I = np.flatnonzero(dr7.filter == band)
            for x,y,d in zip(dr7.ix[I],dr7.iy[I],dr7.galdepth[I]):
                depth[y,x] = max(depth[y,x], d)
    
        gdepth[gdepth == 0] = np.nan
        rdepth[rdepth == 0] = np.nan
        zdepth[zdepth == 0] = np.nan
    
        cm = matplotlib.cm.viridis
        cm.set_bad('0.9')
        depths = dict(g=gdepth, r=rdepth, z=zdepth)
        for band in 'grz':
    
            d = depths[band]
            m = np.median(d[np.isfinite(d)])
            print('Median galdepth in', band, 'is', m)
    
            plt.clf()
            nom = dict(g=24, r=23.4, z=22.5)[band]
            mn,mx = nom-1, nom+1
            plt.imshow(depths[band], interpolation='nearest', origin='lower',
                       vmin=mn, vmax=mx, cmap=cm)
            aitoff_labels(wcs)
            cb = plt.colorbar()
            cb.set_label('%s depth (mag)' % band)
            plt.title('DR6 + DR7')
            plt.savefig('depth-%s.pdf' % band)

    # Best seeing
    if True:

        # W,H = 1000,500
        # zoom = 1.0
        # wcs = anwcs_create_hammer_aitoff(270., 0., zoom, W, H, False)

        from astrometry.util.util import Tan

        dr6g = fits_table('ccds-annotated-dr6-g_AD.fits.gz')
        dr6r = fits_table('ccds-annotated-dr6-r_AD.fits.gz')
        dr6z = fits_table('ccds-annotated-dr6-z_AD.fits.gz')
        dr7 = fits_table('ccds-annotated-dr7_AD.fits.gz')
        dr6 = merge_tables([dr6g,dr6r,dr6z])

        dr7.fwhm = dr7.fwhm_cp
        dr6.fwhm = dr6.fwhm_cp

        print(len(dr7), 'DR7')
        dr7.cut(dr7.fwhm > 2.0)
        print('Cut to', len(dr7), 'on FWHM>2')

        print(len(dr6), 'DR6')
        dr6.cut(dr6.fwhm > 2.0)
        print('Cut to', len(dr6), 'on FWHM>2')
        
        if False:
            # dr6 = fits_table('/data1/dr6/survey-ccds-dr6plus.kd.fits',
            #                  columns=['image_filename', 'ccdname',
            #                           'filter','ra','dec','fwhm',
            #                           'width','height','crval1','crval2',
            #                           'crpix1','crpix2','cd1_1','cd1_2',
            #                           'cd2_1','cd2_2'])
            # print(len(dr6), 'DR6')
            dr7 = fits_table('/data1/dr7/survey-ccds-dr7.kd.fits',
                             columns=['image_filename', 'ccdname',
                                      'filter','ra','dec','fwhm','ccd_cuts',
                                      'width','height','crval1','crval2',
                                      'crpix1','crpix2','cd1_1','cd1_2',
                                      'cd2_1','cd2_2'])
    
            dr7 = dr7[:1000]
            print(len(dr7), 'DR7')
            
            #90prime/CP20160117/ksb_160118_130951_ooi_r_v1.fits.fz
            
    
            bok = fits_table('seeing/DR8_bok_header_info.fits')
            print(len(bok), 'Bok')
            bok.cut(bok.cenra1 != 0.0)
            print(len(bok), 'ok')
    
            # Cut on date -- DR6 includes Bok data from before June 25, 2017.
            # and MzLS data before December 9, 2017.
    
            keep = np.array([i for i,fn in enumerate(bok.filename) if (fn.strip().split('/')[-1].split('_')[1] < '170625')])
            bok.cut(keep)
            print(len(bok), 'kept based on date')
            
            # def fn_norm(fn):
            #     return '_'.join(fn.strip().split('/')[-1].split('_')[:3])
            # 
            # fnmap = dict([((fn_norm(fn), ccdname.strip().upper()), i)
            #               for i,(fn,ccdname)
            #               in enumerate(zip(see.filename, see.ccdname))])
            # 
            # I = np.array([fnmap.get((fn_norm(fn), ccdname.strip().upper()), -1)
            #               for fn,ccdname in zip(dr6.image_filename, dr6.ccdname)])
            # print('Of', len(dr6), 'DR6 CCDs and', len(see), 'Bok, matched', np.sum(I>-1), 'CCDs')
            # 
            # return
    
    
            # tablist /data1/dr6/survey-ccds-dr6plus.kd.fits"[expnum==75110186][col expnum;ccdname;crpix1]"
            wcskeys = {
                'crpix1': dict(CCD1=4103.19, CCD2=-87.1133, CCD3=4124.77,
                               CCD4=-68.2412),
                'crpix2': dict(CCD1=4290.66, CCD2=4283.26, CCD3=-185.110,
                               CCD4=-198.428),
                'cd1_1': dict(CCD1=-6.25650E-07, CCD2=-1.16688E-06,
                              CCD3=-1.06214E-06, CCD4=-8.99536E-07),
                'cd1_2': dict(CCD1=-0.000125677, CCD2=-0.000125671,
                              CCD3=-0.000125729, CCD4=-0.000125735),
                'cd2_1': dict(CCD1=-0.000125699, CCD2=-0.000125717,
                              CCD3=-0.000125698, CCD4=-0.000125759),
                'cd2_2': dict(CCD1=1.12564E-06, CCD2=6.04536E-07,
                              CCD3=5.06435E-07, CCD4=1.37220E-06),
                              }
            width=4032
            height=4096
    
            for k in ['crpix1', 'crpix2', 'cd1_1', 'cd1_2', 'cd2_1', 'cd2_2']:
                bok.set(k, np.zeros(len(bok)))
                arr = bok.get(k)
                vals = wcskeys[k]
                for i,ccdname in enumerate(bok.ccdname):
                    arr[i] = vals[ccdname.strip().upper()]
            bok.width  = np.zeros(len(bok)) + 4032
            bok.height = np.zeros(len(bok)) + 4096
            bok.ra  = bok.cenra1
            bok.dec = bok.cendec1
            bok.rename('cenra1', 'crval1')
            bok.rename('cendec1', 'crval2')
            bok.rename('seeingp1', 'fwhm')
            bok.filter = np.array([fn.strip().split('/')[-1].split('_')[-2]
                                   for fn in bok.filename])
            print('Bok filter', Counter(bok.filter))
    
            dr6x = fits_table('/data1/dr6/survey-ccds-dr6plus.kd.fits',
                              columns=['expnum', 'ccdname',
                                       'filter','ra','dec','fwhm',
                                       'width','height','crval1','crval2',
                                       'crpix1','crpix2','cd1_1','cd1_2',
                                       'cd2_1','cd2_2'])
            dr6x.cut(dr6x.expnum == 107621)
            assert(len(dr6x) == 4)
    
            wcskeys = {}
            for k in ['crpix1', 'crpix2', 'cd1_1', 'cd1_2', 'cd2_1', 'cd2_2',
                      'width', 'height']:        
                vals = dict()
                for ccdname,val in zip(dr6x.ccdname, dr6x.get(k)):
                    vals[ccdname.strip().upper()] = val
                wcskeys[k] = vals
            print('WCS keys:', wcskeys)
    
            mos = fits_table('seeing/DR8_mosaic_header_info.fits')
            print(len(mos), 'Mos')
            mos.cut(mos.cenra1 != 0.0)
            print(len(mos), 'ok')
    
            keep = np.array([i for i,fn in enumerate(mos.filename) if (fn.strip().split('/')[-1].split('_')[1] < '171209')])
            mos.cut(keep)
            print(len(mos), 'kept based on date')
    
            for k in ['crpix1', 'crpix2', 'cd1_1', 'cd1_2', 'cd2_1', 'cd2_2',
                      'width', 'height']:
                mos.set(k, np.zeros(len(mos)))
                arr = mos.get(k)
                vals = wcskeys[k]
                for i,ccdname in enumerate(mos.ccdname):
                    arr[i] = vals[ccdname.strip().upper()]
            mos.ra  = mos.cenra1
            mos.dec = mos.cendec1
            mos.rename('cenra1', 'crval1')
            mos.rename('cendec1', 'crval2')
            mos.rename('seeingp1', 'fwhm')
            mos.filter = np.array([fn.strip().split('/')[-1].split('_')[-2][0]
                                   for fn in mos.filename])
    
            print('Mos filter:', Counter(mos.filter))
    
            dr7x = fits_table('/data1/dr7/survey-ccds-dr7.kd.fits',
                              columns=['expnum', 'ccdname',
                                       'filter','ra','dec','fwhm',
                                       'width','height','crval1','crval2',
                                       'crpix1','crpix2','cd1_1','cd1_2',
                                       'cd2_1','cd2_2'])
            #dr7x.cut(dr7x.expnum == 425090)
            #dr7x.cut(dr7x.expnum == 508199)
            dr7x.cut(dr7x.expnum == 663609)
            wcskeys = {}
            for k in ['crpix1', 'crpix2', 'cd1_1', 'cd1_2', 'cd2_1', 'cd2_2',
                      'width', 'height']:        
                vals = dict()
                for ccdname,val in zip(dr7x.ccdname, dr7x.get(k)):
                    vals[ccdname.strip().upper()] = val
                wcskeys[k] = vals
            #print('WCS keys:', wcskeys)
    
            dec = fits_table('seeing/DR7_header_info.fits')
            print(len(dec), 'Dec')
            dec.cut(dec.cenra1 != 0.0)
            print(len(dec), 'ok')
            for k in ['crpix1', 'crpix2', 'cd1_1', 'cd1_2', 'cd2_1', 'cd2_2',
                      'width', 'height']:
                dec.set(k, np.zeros(len(dec)))
                arr = dec.get(k)
                vals = wcskeys[k]
                for i,ccdname in enumerate(dec.ccdname):
                    arr[i] = vals[ccdname.strip().upper()]
            dec.ra  = dec.cenra1
            dec.dec = dec.cendec1
            dec.rename('cenra1', 'crval1')
            dec.rename('cendec1', 'crval2')
            dec.filter = np.array([fn.strip().split('/')[-1].split('_')[-2][0]
                                   for fn in dec.filename])
    
            dec.ccd_cuts = np.zeros(len(dec), np.int32)
        
            print('Dec filter:', Counter(dec.filter))

            for filt in ['g','r']:
                plt.clf()
                plt.hist(bok.fwhm[bok.filter == filt] * 0.454, range=(0,5), bins=50)
                plt.title('Bok %s' % filt)
                plt.savefig('see-bok-%s.png' % filt)
    
            plt.clf()
            plt.hist(mos.fwhm * 0.262, range=(0,5), bins=50)
            plt.title('Mos z')
            plt.savefig('see-mos-z.png')
    
            for filt in ['g','r','z']:
                plt.clf()
                plt.hist(dec.fwhm[dec.filter == filt] * 0.262, range=(0,5), bins=50)
                plt.title('DECam %s' % filt)
                plt.savefig('see-dec-%s.png' % filt)
    
            dr6 = merge_tables([bok, mos])
            dr7 = dec

        
        dr6_all = dr6.copy()
        dr7_all = dr7.copy()
        
        for band in 'grz':
            dr6 = dr6_all[(dr6_all.filter == band)]
            print(len(dr6), 'DR6', band)
            dr7 = dr7_all[(dr7_all.filter == band)]
            print(len(dr7), 'DR7', band)
            dr7.cut((dr7.ccd_cuts == 0))
            print(len(dr7), 'DR7 good')

            ok,x,y = wcs.radec2pixelxy(dr6.ra, dr6.dec)
            dr6.ix = np.round(x-1).astype(int)
            dr6.iy = np.round(y-1).astype(int)
            ok,x,y = wcs.radec2pixelxy(dr7.ra, dr7.dec)
            dr7.ix = np.round(x-1).astype(int)
            dr7.iy = np.round(y-1).astype(int)

            if band == 'z':
                # MzLS
                pixscale6 = 0.262
            else:
                # BASS
                pixscale6 = 0.454

            dr6.cut((dr6.ix >= 0) * (dr6.ix < W) * (dr6.iy >= 0) * (dr6.iy < H))
            print('Cut to', len(dr6), 'in-bounds')
            dr6.see = dr6.fwhm * pixscale6

            dr7.cut((dr7.ix >= 0) * (dr7.ix < W) * (dr7.iy >= 0) * (dr7.iy < H))
            print('Cut to', len(dr7), 'in-bounds')
            dr7.see = dr7.fwhm * 0.262

            print('Median seeing in DR6', band, ':', np.median(dr6.see))
            print('Median seeing in DR7', band, ':', np.median(dr7.see))
            
            dr6.x1 = np.zeros(len(dr6), int)
            dr6.x2 = np.zeros(len(dr6), int)
            dr6.y1 = np.zeros(len(dr6), int)
            dr6.y2 = np.zeros(len(dr6), int)
            for i,ccd in enumerate(dr6):
                ccdwcs = Tan(*[float(x) for x in
                    [ccd.crval1, ccd.crval2, ccd.crpix1, ccd.crpix2,
                     ccd.cd1_1,  ccd.cd1_2,  ccd.cd2_1, ccd.cd2_2,
                     ccd.width, ccd.height]])
                rr,dd = ccdwcs.pixelxy2radec(
                    [1, 1, ccd.width, ccd.width, 1],
                    [1, ccd.height, ccd.height, 1, 1])
                ok,xx,yy = wcs.radec2pixelxy(rr, dd)
                ix = np.round(xx).astype(int)
                iy = np.round(yy).astype(int)
                dr6.x1[i] = min(ix)-1
                dr6.x2[i] = max(ix)
                dr6.y1[i] = min(iy)-1
                dr6.y2[i] = max(iy)

            # Cut small number of wrap-arounds
            dr6.cut((dr6.x2 - dr6.x1) < 100)
            print('Cut to', len(dr6), 'non-wrapping')

            
            bestseeing = np.empty((H,W), np.float32)
            bestseeing[:,:] = 10.
            #for x,y,see in zip(dr6.ix, dr6.iy, dr6.see):
            #    bestseeing[y,x] = min(bestseeing[y,x], see)

            for x1,x2,y1,y2,see in zip(dr6.x1, dr6.x2, dr6.y1, dr6.y2, dr6.see):
                bestseeing[y1:y2,x1:x2] = np.minimum(bestseeing[y1:y2,x1:x2], see)

            for x,y,see in zip(dr7.ix, dr7.iy, dr7.see):
                bestseeing[y,x] = min(bestseeing[y,x], see)
            bestseeing[bestseeing == 10.] = np.nan

            plt.clf()
            plt.hist(bestseeing[np.isfinite(bestseeing)], range=(0,3), bins=60)
            plt.title('Best seeing: %s' % band)
            plt.savefig('see-best-%s.png' % band)

            
            cm = matplotlib.cm.hot
            cm.set_bad('0.9')
            plt.clf()
            plt.imshow(bestseeing, interpolation='nearest', origin='lower',
                       vmin=0.6, vmax=2.0, cmap=cm)
            aitoff_labels(wcs)
            cb = plt.colorbar()
            cb.set_label('Best %s-band seeing (arcsec)' % band)
            plt.title('DR6 + DR7')
            plt.savefig('best-seeing-%s.pdf' % band)

            
if __name__ == '__main__':

    print('zsee_fig')
    zsee_fig()
    import sys
    sys.exit(0)
    
    def aitoff_labels(wcs):
        # image bounds
        declo = -25
        dechi = +35

        ok,x,y = wcs.radec2pixelxy(0., declo)
        ylo = int(y-1)
        ok,x,y = wcs.radec2pixelxy(0., dechi)
        yhi = int(y-1)
        ok,x,y = wcs.radec2pixelxy([80, 300.], [0.,0.])
        xlo = min(x)
        xhi = max(x)
        print('x range', xlo,xhi)    
        
        # RA labels
        lt,lv = [],[]
        for ra in [89.99, 60,30,0,330,300,270]:
            decs = np.linspace(declo-10, declo+10, 25)
            ok,x,y = wcs.radec2pixelxy(ra, decs)
            i = np.argmin(np.abs(y - ylo))
            x = x[i]
            lt.append('%i' % np.round(ra))
            lv.append(x)
        plt.xticks(lv, lt)
        plt.yticks([])
    
        for ra in [89.999, 60,30,0,330,300,270]:
            decs = np.linspace(declo-5, +90., 100)
            ok,x,y = wcs.radec2pixelxy(ra + np.zeros_like(decs), decs)
            plt.plot(x-1, y-1, 'k-', alpha=0.3)
    
        for dec in [-15, 0, 15, 30]:
            ras = np.append([89.999], np.linspace(89.999, -60, 100))
            ok,x,y = wcs.radec2pixelxy(ras, np.zeros_like(ras)+dec)
            plt.plot(x-1, y-1, 'k-', alpha=0.3)
    
        for dec in [-15, 0, 15, 30]:
            ra = 300
            ok,x,y = wcs.radec2pixelxy(ra, dec)
            ha = 'left'
            dx = +5
            plt.text(x+dx, y-5, '$%+i$' % dec, va='top', ha=ha)
        plt.axis([xlo-3,xhi+3,ylo-3,yhi+1])
        plt.xlabel('RA (deg)')

    
    matplotlib.rc('text', usetex=True)
    matplotlib.rc('font', family='serif')
    plt.figure(1, figsize=(13,4))
    plt.subplots_adjust(left=0.03, right=0.99, bottom=0.08, top=0.95)
    
    W,H = 2000,1000
    zoom = 2.5
    wcs = anwcs_create_hammer_aitoff(10., 5., zoom, W, H, False)
    from astrometry.util.util import Tan


    bestall = fitsio.read('bestseeing.fits')

    bestall[bestall == 10.] = np.nan
    cm = matplotlib.cm.jet
    cm = cmap_discretize(cm, 14)
    cm.set_bad('0.9')
    plt.clf()
    plt.imshow(bestall, interpolation='nearest', origin='lower',
               vmin=0.8, vmax=1.5, cmap=cm)
    #ax = plt.axis()
    aitoff_labels(wcs)
    #plt.axis(ax)
    cb = plt.colorbar()
    cb.set_label('Best seeing (arcsec)')
    plt.title('DR7 + tile file')
    plt.savefig('bsee.pdf')

    import sys
    sys.exit(0)


    
    dr7 = fits_table('/data1/dr7/survey-ccds-dr7.kd.fits',
                     columns=['filter','ra','dec','fwhm','ccd_cuts',
                              'width','height','crval1','crval2',
                              'crpix1','crpix2','cd1_1','cd1_2',
                              'cd2_1','cd2_2', 'expnum'])
    print(len(dr7), 'DR7')
    dr7_all = dr7.copy()

    expnums = set(dr7.expnum)
    print(len(expnums), 'DR7 expnums')
    
    db = fits_table('db.fits')
    print('Got', len(db), 'copilot database exposures')
    print('Expnum range', db.expnum.min(), db.expnum.max())
    db.cut(db.expnum > 0)
    print(len(db), 'good expnums')
    E = np.unique(db.expnum)
    print(len(E), 'unique expnums')
    print('Expnum range', db.expnum.min(), db.expnum.max())
    print('Cameras:', Counter(db.camera))
    print('Obstype:', Counter(db.obstype))
    db.cut(db.obstype == 'object   ')
    print(len(db), 'obstype OBJECT')
    E,I = np.unique(db.expnum, return_index=True)
    db.cut(I)
    print(len(db), 'unique expnum')

    I = np.array([e not in expnums for e in db.expnum])
    db.cut(I)
    print(Counter(db.band))
    db.cut(np.array([b.strip() in ['g','r','z'] for b in db.band]))
    print(Counter(db.band))
    print(len(db), 'not in DR7')

    db.cut(db.exptime > 30.)
    print(len(db), 'with exptimes > 30 s')
    
    #print('Exposure times:', Counter(db.exptime))
    print('Exposure times:')
    c = Counter(db.exptime)
    for t in sorted(list(c.keys())):
        print('  ', t, c[t])
    
    db_all = db.copy()
    
    # plt.clf()
    # plt.plot(db.mjd_obs, 'b.')
    # plt.savefig('mjd.png')

    F=fitsio.FITS('decam.wcs.gz')
    #H,W = 4094,2046
    nccd = len(F)-1
    C = fits_table()
    C.crpix1 = np.zeros(nccd, np.float32)
    C.crpix2 = np.zeros(nccd, np.float32)
    C.cd1_1 = np.zeros(nccd, np.float32)
    C.cd1_2 = np.zeros(nccd, np.float32)
    C.cd2_1 = np.zeros(nccd, np.float32)
    C.cd2_2 = np.zeros(nccd, np.float32)
    for i in range(nccd):
        hdr = F[i+1].read_header()
        C.crpix1[i] = hdr['CRPIX1']
        C.crpix2[i] = hdr['CRPIX2']
        C.cd1_1[i] = hdr['CD1_1']
        C.cd1_2[i] = hdr['CD1_2']
        C.cd2_1[i] = hdr['CD2_1']
        C.cd2_2[i] = hdr['CD2_2']
    C.width = np.zeros(len(C), np.int32) + 2046
    C.height = np.zeros(len(C), np.int32)+ 4094

    print('Assuming', nccd, 'CCDs')

    bestall = np.empty((H,W), np.float32)
    bestall[:,:] = 10.
    
    for band in 'grz':
        dr7 = dr7_all[(dr7_all.filter == band)]
        print(len(dr7), 'DR7', band)
        dr7.cut((dr7.ccd_cuts == 0))
        print(len(dr7), 'DR7 good')
        ok,x,y = wcs.radec2pixelxy(dr7.ra, dr7.dec)
        dr7.ix = np.round(x-1).astype(int)
        dr7.iy = np.round(y-1).astype(int)
        dr7.cut((dr7.ix >= 0) * (dr7.ix < W) * (dr7.iy >= 0) * (dr7.iy < H))
        print('Cut to', len(dr7), 'in-bounds')
        dr7.see = dr7.fwhm * 0.262

        dr7.x1 = np.zeros(len(dr7), int)
        dr7.x2 = np.zeros(len(dr7), int)
        dr7.y1 = np.zeros(len(dr7), int)
        dr7.y2 = np.zeros(len(dr7), int)
        for i,ccd in enumerate(dr7):
            ccdwcs = Tan(*[float(x) for x in
                [ccd.crval1, ccd.crval2, ccd.crpix1, ccd.crpix2,
                 ccd.cd1_1,  ccd.cd1_2,  ccd.cd2_1, ccd.cd2_2,
                 ccd.width, ccd.height]])
            rr,dd = ccdwcs.pixelxy2radec(
                [1, 1, ccd.width, ccd.width, 1],
                [1, ccd.height, ccd.height, 1, 1])
            ok,xx,yy = wcs.radec2pixelxy(rr, dd)
            ix = np.round(xx).astype(int)
            iy = np.round(yy).astype(int)
            dr7.x1[i] = min(ix)-1
            dr7.x2[i] = max(ix)
            dr7.y1[i] = min(iy)-1
            dr7.y2[i] = max(iy)

        bestseeing = np.empty((H,W), np.float32)
        bestseeing[:,:] = 10.
        #for x,y,see in zip(dr7.ix, dr7.iy, dr7.see):
        #    bestseeing[y,x] = min(bestseeing[y,x], see)

        for x1,x2,y1,y2,see in zip(dr7.x1, dr7.x2, dr7.y1, dr7.y2, dr7.see):
            bestseeing[y1:y2,x1:x2] = np.minimum(bestseeing[y1:y2,x1:x2], see)


        best1 = bestseeing.copy()
        best1[best1 == 10.] = np.nan
        cm = matplotlib.cm.hot
        cm.set_bad('0.9')
        plt.clf()
        plt.imshow(best1, interpolation='nearest', origin='lower',
                   vmin=0.8, vmax=2.0, cmap=cm)
        #aitoff_labels(wcs)
        cb = plt.colorbar()
        cb.set_label('Best %s-band seeing (arcsec)' % band)
        plt.title('DR7')
        plt.savefig('bsee1-%s.pdf' % band)

        db = db_all[np.array([b.strip() == band for b in db_all.band])]
        print(len(db), 'DB', band)
        ok,x,y = wcs.radec2pixelxy(db.rabore, db.decbore)
        db.ix = np.round(x-1).astype(int)
        db.iy = np.round(y-1).astype(int)
        db.cut((db.ix >= 0) * (db.ix < W) * (db.iy >= 0) * (db.iy < H))
        print('Cut to', len(db), 'in-bounds')

        for i,dbi in enumerate(db):
            for j,ccd in enumerate(C):
                ccdwcs = Tan(*[float(x) for x in
                               [dbi.rabore, dbi.decbore, ccd.crpix1, ccd.crpix2,
                                ccd.cd1_1,  ccd.cd1_2,  ccd.cd2_1, ccd.cd2_2,
                                ccd.width, ccd.height]])
                rr,dd = ccdwcs.pixelxy2radec(
                    [1, 1, ccd.width, ccd.width, 1],
                    [1, ccd.height, ccd.height, 1, 1])
                ok,xx,yy = wcs.radec2pixelxy(rr, dd)
                ix = np.round(xx).astype(int)
                iy = np.round(yy).astype(int)
                x1 = min(ix)-1
                x2 = max(ix)
                y1 = min(iy)-1
                y2 = max(iy)
                bestseeing[y1:y2,x1:x2] = np.minimum(bestseeing[y1:y2,x1:x2], dbi.seeing)

        bestseeing[bestseeing == 10.] = np.nan

        cm = matplotlib.cm.hot
        cm.set_bad('0.9')
        plt.clf()
        plt.imshow(bestseeing, interpolation='nearest', origin='lower',
                   vmin=0.8, vmax=2.0, cmap=cm)
        #aitoff_labels(wcs)
        cb = plt.colorbar()
        cb.set_label('Best %s-band seeing (arcsec)' % band)
        plt.title('DR7 + tile file')
        plt.savefig('bsee-%s.pdf' % band)

        bestall = np.minimum(bestall, bestseeing)

    fitsio.write('bestseeing.fits', bestall, clobber=True)

    bestall[bestall == 10.] = np.nan
    cm = matplotlib.cm.hot
    cm.set_bad('0.9')
    plt.clf()
    plt.imshow(bestall, interpolation='nearest', origin='lower',
               vmin=0.8, vmax=2.0, cmap=cm)
    cb = plt.colorbar()
    cb.set_label('Best seeing (arcsec)')
    plt.title('DR7 + tile file')
    plt.savefig('bsee.pdf')

    import sys
    sys.exit(0)

    print('zsee_fig')
    zsee_fig()
    import sys
    sys.exit(0)
    
    fn ='maps/depth-fake-g.fits.fz'
    if not os.path.exists(fn):
        d1 = fitsio.read('maps/depth-decals-all-g.fits.fz')
        d2 = fitsio.read('maps/depth-des-fake-g.fits.fz')
        dsum = -2.5 * (np.log10(1./np.sqrt(
            1./(10.**((d1 - 22.5) / -2.5))**2 +
            1./(10.**((d2 - 22.5) / -2.5))**2
            )) - 9.)
        del d1,d2
        wcsfn = 'cea-flip.wcs'
        hdr = fitsio.read_header(wcsfn)
        # HACK
        H,W = 2600, 11001
        hdr['IMAGEW'] = W
        hdr['IMAGEH'] = H
        fitsio.write('/tmp/depth-fake-g.fits.fz', dsum, clobber=True, header=hdr)
        del dsum
        cmd = 'imcopy /tmp/depth-fake-g.fits.fz"[20000:31000, 1:2600]" %s' % fn
        print(cmd)
        os.system(cmd)

    depthmap = fitsio.read(fn)
    wcs = anwcs(fn)

    #plt.figure(1, figsize=(13,8))
    plt.figure(1, figsize=(13,6))
    plt.subplots_adjust(left=0.03, right=0.99, bottom=0.05, top=0.95)
    plt.clf()
    plt.imshow(depthmap, interpolation='nearest', origin='lower',
               vmin=23, vmax=25, cmap='bone')
    cbar = plt.colorbar(orientation='horizontal')
    plt.savefig('depth-fake-g.png')

    xy = [
    # Bottom of the barrel
        (10131, 1538),
        (9816,  1545),
        (9516,  1545),
        (9244,  1537),
        (8944,  1537),
        (8628,  1537),
        (8348,  1557),
        (8054,  1555),
        (7736,  1555),
        (7427,  1555),
        (6503,  1511),
    # Top of the barrel
        (10009, 1968),
        (9713,  1960),
        (9415,  1958),
        (9119,  1958),
        (8825,  1970),
        (8517,  1968),
        (8207,  1966),
        (7907,  1958),
        (7015,  1972),
        (6721,  1982),
    # Top of the tank
        (4753, 2246),
        (4449, 2246),
        (4143, 2246),
        (3843, 2246),
        (3581, 2246),
        (3261, 2246),
        (2967, 2246),
        (2647, 2246),
        (2355, 2246),
    # Corner
        (2055, 908),
        ]

    plt.plot([x for x,y in xy], [y for x,y in xy], 'r.')
    ax = plt.axis()

    xt,xl = [],[]
    for r in [300,310,320,330,340,350,0,10,20,30,40,50,60]:
        ok,x,y = wcs.radec2pixelxy(r, 0)
        xt.append(x)
        xl.append('%i' % r)
    plt.xticks(xt, xl)
    plt.xlabel('RA')
    yt,yl = [],[]
    for d in [-15,-10,-5,0,5,10,15]:
        ok,x,y = wcs.radec2pixelxy(0, d)
        yt.append(y)
        yl.append('%i' % d)
    plt.yticks(yt, yl)
    plt.ylabel('Dec')

    plt.axis(ax)

    plt.savefig('depth-fake-g-1.png')


    F=fitsio.FITS('decam.wcs.gz')
    H,W = 4094,2046
    nccd = len(F)-1
    C = fits_table()
    C.crpix1 = np.zeros(nccd, np.float32)
    C.crpix2 = np.zeros(nccd, np.float32)
    C.cd1_1 = np.zeros(nccd, np.float32)
    C.cd1_2 = np.zeros(nccd, np.float32)
    C.cd2_1 = np.zeros(nccd, np.float32)
    C.cd2_2 = np.zeros(nccd, np.float32)
    for i in range(nccd):
        hdr = F[i+1].read_header()
        C.crpix1[i] = hdr['CRPIX1']
        C.crpix2[i] = hdr['CRPIX2']
        C.cd1_1[i] = hdr['CD1_1']
        C.cd1_2[i] = hdr['CD1_2']
        C.cd2_1[i] = hdr['CD2_1']
        C.cd2_2[i] = hdr['CD2_2']
    print('Assuming', nccd, 'CCDs')

    H,W = depthmap.shape

    band = 'g'
    target = target_depths[band]
    depth = target - 0.3

    detsig = 10.**((depth - 22.5) / -2.5)
    #depthfrac = 1./detsig**2 / targetiv
    newiv = 1./detsig**2

    print('newiv', newiv)

    #targetsig = 10.**((targetdepth - 22.5) / -2.5)
    #targetiv = 1./targetsig**2

    get_polygon = get_decam_outline_function(close=True)
    
    plot = Plotstuff(size=(W,H), outformat='png')
    plotwcs = anwcs(fn)
    plot.wcs = plotwcs

    #print('Plot WCS:', plot.wcs)
    #anwcs_print_stdout(plot.wcs)

    plot.outline.fill = True
    plot.outline.stepsize = 2000
    #plot.op = CAIRO_OPERATOR_ADD

    ccds = fits_table()
    c = np.zeros(nccd, np.float32)
    for col in ['crpix1','crpix2','cd1_1','cd1_2','cd2_1','cd2_2']:
        ccds.set(col, C.get(col))
    ccds.crval1 = np.zeros(len(ccds))
    ccds.crval2 = np.zeros(len(ccds))
    #ccds.depth = np.zeros(len(ccds), np.float32) + depth
    ccds.width  = np.zeros(len(ccds), np.int32) + 2046
    ccds.height = np.zeros(len(ccds), np.int32) + 4094

    depth_iv = 1./(10.**((depthmap - 22.5) / -2.5))**2
    depth_iv[depthmap < 10] = 0.

    from legacypipe.survey import LegacySurveyData
    survey = LegacySurveyData()

    tan_wcses = [survey.get_approx_wcs(ccd) for ccd in ccds]

    k = 1

    class duck(object):
        pass

    quack = duck()

    tiles = fits_table()
    tiles.ra = []
    tiles.dec = []
    
    for xx,yy in xy:

        #x0 = xx - 200
        #y0 = yy - 200
        #subimg = slice(yy - 200, yy + 200), slice(xx - 200, xx + 200)
        subimg = slice(yy - 150, yy + 150), slice(xx - 150, xx + 150)

        # ok,ra,dec = wcs.pixelxy2radec(xx, yy)
        # quack.ra = ra
        # quack.dec = dec
        # poly = get_polygon(quack, wcs)
        # print('WCS shape', wcs.shape)
        print('x,y', xx,yy)
        
        DY = range(-25, 25+1)#, 5)
        DX = range(-25, 25+1)#, 5)

        metric = np.zeros((len(DY), len(DX)))
        metric99 = np.zeros((len(DY), len(DX)))

        #for dx in range(-10, 10+1):
        for ix,dx in enumerate(DX):
            for iy,dy in enumerate(DY):
                #for dy in range(-10, 10+1):
                plot.color = 'black'
                plot.plot('fill')
                plot.color = 'white'
                plot.alpha = 1.

                ok,ra,dec = wcs.pixelxy2radec(xx+dx, yy+dy)

                ccds.crval1[:] = ra
                ccds.crval2[:] = dec

                for twcs in tan_wcses:
                    twcs.set_crval(ra, dec)
                    plot.outline.wcs = anwcs_new_tan(Tan(twcs))
                    plot.plot('outline')

                img = plot.get_image_as_numpy_view()

                #iv = img[::-1,:,0] * (1./255. * newiv)
                #iv = img[:,:,0] * (1./255. * newiv)
                #print('iv', iv.dtype, iv.min(), iv.max())

                iv = img[:,:,0][subimg] * (1./255. * newiv)

                #print('Median existing depth iv:', np.median(depth_iv[subimg]))
                #print('newiv', newiv)

                ivsum = depth_iv[subimg] + iv
                dsum = -2.5 * (np.log10(1./np.sqrt(
                    ivsum
                    )) - 9.)
                dsum[ivsum == 0] = 0.

                # min?  1%ile?
                metric[iy,ix] = dsum.min()

                metric99[iy,ix] = np.percentile(dsum.ravel(), 1)

                # plt.clf()
                # plt.subplot(2,1,1)
                # plt.imshow(dsum, interpolation='nearest', origin='lower',
                #            vmin=23, vmax=25, cmap='bone')
                # #plt.plot(poly[:,0] - x0, poly[:,1] - y0, 'k-')
                # cbar = plt.colorbar()
                # plt.subplot(2,1,2)
                # plt.hist(dsum.ravel(), range=(0,25), bins=100)
                # ofn = 'depth-fake-g-%02i.png' % (k)
                # plt.savefig(ofn)
                # print('Wrote', ofn)
                # k += 1

        # plt.clf()
        # plt.imshow(metric, interpolation='nearest', origin='lower')
        # plt.colorbar()
        # plt.title('min depth')
        # ofn = 'depth-fake-g-%02i.png' % (k)
        # plt.savefig(ofn)
        # print('Wrote', ofn)
        # k += 1

        ibest = np.argmax(metric)
        iy,ix = np.unravel_index(ibest, metric.shape)
        dy = DY[iy]
        dx = DX[ix]
        print('Best for min:', dx,dy,  metric99[iy,ix], metric[iy,ix])

        ibest = np.argmax(metric99)
        iy,ix = np.unravel_index(ibest, metric.shape)
        dy = DY[iy]
        dx = DX[ix]
        print('Best for 99%:', dx,dy, metric99[iy,ix], metric[iy,ix])

        ibest = np.argmax(metric99 + metric)
        iy,ix = np.unravel_index(ibest, metric.shape)
        dy = DY[iy]
        dx = DX[ix]
        print('Best for sum:', dx,dy, metric99[iy,ix], metric[iy,ix])

        plot.color = 'black'
        plot.plot('fill')
        plot.color = 'white'
        plot.alpha = 1.

        ok,ra,dec = wcs.pixelxy2radec(xx+dx, yy+dy)

        ccds.crval1[:] = ra
        ccds.crval2[:] = dec

        for twcs in tan_wcses:
            twcs.set_crval(ra, dec)
            plot.outline.wcs = anwcs_new_tan(Tan(twcs))
            plot.plot('outline')

        img = plot.get_image_as_numpy_view()
        iv = img[:,:,0][subimg] * (1./255. * newiv)
        ivsum = depth_iv[subimg] + iv
        dsum = -2.5 * (np.log10(1./np.sqrt(
            ivsum
            )) - 9.)
        dsum[ivsum == 0] = 0.

        imin = np.argmin(dsum)
        iy,ix = np.unravel_index(imin, dsum.shape)

        quack.ra = ra
        quack.dec = dec
        poly = get_polygon(quack, wcs)

        dpre = -2.5 * (np.log10(1./np.sqrt(
            depth_iv[subimg]
            )) - 9.)
        dpre[depth_iv[subimg] == 0] = 0.

        tiles.ra.append(ra)
        tiles.dec.append(dec)

        x0 = xx - 150
        y0 = yy - 150

        plt.clf()
        plt.subplot(1,2,1)
        plt.imshow(dpre, interpolation='nearest', origin='lower',
                   vmin=23, vmax=25, cmap='bone')
        plt.title('99%%/100%% coverage: %.2f / %.2f' %
                  (np.percentile(dpre, 1), dpre.min()))
        plt.subplot(1,2,2)
        plt.imshow(dsum, interpolation='nearest', origin='lower',
                   vmin=23, vmax=25, cmap='bone')
        #plt.title('99%% coverage: %.2f' % np.percentile(dsum, 99))
        plt.title('99%%/100%% coverage: %.2f / %.2f' %
                  (np.percentile(dsum, 1), dsum.min()))
        #plt.title('min depth: %.2f' % dsum.min())
        plt.plot(poly[:,0] - x0, poly[:,1] - y0, 'k-')
        plt.plot(ix, iy, 'o', mec='m', mfc='none')
        #cbar = plt.colorbar()
        # plt.subplot(2,1,2)
        # plt.hist(dsum.ravel(), range=(0,25), bins=100)
        plt.suptitle('RA,Dec %.3f, %.3f' % (ra, dec))
        ofn = 'depth-fake-g-%02i.png' % (k)
        plt.savefig(ofn)
        print('Wrote', ofn)
        k += 1

    tiles.writeto('tiles-des-manual.fits')
    
    sys.exit(0)
    
    if False:
        W,H = 32000,5000
        pixscale = 340. / W
        ra,dec = 110., 8.
        refx = W/2. + 0.5
        refy = (H/2. + 0.5 - dec / -pixscale)
        args = (ra, 0., refx, refy, pixscale, W, H, True)
        wcs = anwcs_create_cea_wcs(*args)
        wcsfn = 'cea.wcs'
        if os.path.exists(wcsfn):
            print('WCS file exists:', wcsfn)
        else:
            anwcs_write(wcs, wcsfn)

        # ?
        refy2 = (H/2. + 0.5 - dec / pixscale)
        args2 = (ra, 0., refx, refy2, pixscale, W, H, False)
        wcs2 = anwcs_create_cea_wcs(*args2)
        wcsfn = 'cea-flip.wcs'
        if os.path.exists(wcsfn):
            print('WCS file exists:', wcsfn)
        else:
            anwcs_write(wcs2, wcsfn)


    from astrometry.util.multiproc import multiproc

    threads = 12
    #threads = 3

    mp = multiproc(threads)


    print('Setting up plots')
    plt.figure(1, figsize=(13,8))
    #plt.subplots_adjust(left=0.03, right=0.99, bottom=0.01, top=0.95)
    plt.subplots_adjust(left=0.03, right=0.99, bottom=0.05, top=0.95)
    

    #from_ccds(mp, ngc=False, prefix_pat='maps/depth-dr7a-%s-p%i')

    if False:
        tiles = fits_table('obstatus/decam-tiles_obstatus.fits')
        tiles.cut(tiles.in_desi == 1)
        tiles.cut(tiles.in_des == 0)
        # The tilings are the same for the different bands... could just run once and scale...
        tiles.g_done[:] = 0
        tiles.r_done[:] = 0
        tiles.z_done[:] = 0
        fn = 'decals-all-tiles.fits'
        tiles.writeto(fn)
        pattern='maps/depth-decals-all-%s.fits.fz'
        decam_todo_maps(mp, obsfn=fn, pattern=pattern)

    if False:
        # DES maps -- non-deep-fields
        '''
        From querying the NOAO archive with

        SELECT reference, dtpropid, surveyid, release_date, start_date, date_obs, dtpi, ra, dec, telescope, instrument, filter, exposure, obstype, obsmode, proctype, prodtype, seeing, depth, dtacqnam, reference AS archive_file, filesize, md5sum FROM voi.siap WHERE exposure >= 30 AND release_date <= '2018-01-23' AND (proctype = 'InstCal') AND (prodtype IS NULL OR prodtype <> 'png') AND ((telescope = 'ct4m' AND instrument = 'decam'))
        and (prodtype='image')
        and (filter like 'z%')
        and (dec > -16)
        ORDER BY date_obs ASC LIMIT 50000

        split by band to avoid the 50,000 query limit.
        '''
        G = fits_table('decam-g.fits')
        Gdes = G[(G.dtpropid == '2012B-0001') * (G.exposure == 90.)]
        print(len(Gdes), 'DES g')
        R = fits_table('decam-r.fits')
        Rdes = R[(R.dtpropid == '2012B-0001') * (R.exposure == 90.)]
        print(len(Rdes), 'DES r')
        Z = fits_table('decam-z.fits')
        Zdes = Z[(Z.dtpropid == '2012B-0001') * (Z.exposure == 90.)]
        print(len(Zdes), 'DES z')

        des = merge_tables([Gdes, Rdes, Zdes])
        tiles = fits_table()
        tiles.ra = des.ra
        tiles.dec = des.dec
        tiles.filter = np.array([f[0] for f in des.filter])
        tiles.in_desi = np.ones(len(tiles), int)
        tiles.in_des  = np.zeros(len(tiles), int)
        tiles.set('pass', np.ones(len(tiles), int))

        # We're calling the "to-do" code, so mark the other bands as 'done'
        tiles.g_done = np.array([0 if f.startswith('g') else 1 for f in des.filter])
        tiles.r_done = np.array([0 if f.startswith('r') else 1 for f in des.filter])
        tiles.z_done = np.array([0 if f.startswith('z') else 1 for f in des.filter])

        print('g_done', Counter(tiles.g_done).most_common())
        print('r_done', Counter(tiles.r_done).most_common())
        print('z_done', Counter(tiles.z_done).most_common())

        # Assume 90-sec DES exposures are to our nominal depth???

        fn = 'des-tiles.fits'
        tiles.writeto(fn)
        pattern='maps/depth-des-fake-%s.fits.fz'
        decam_todo_maps(mp, obsfn=fn, pattern=pattern)

    if False:
        update_decam_tiles(done_pattern='maps/depth-decals-all-%s.fits.fz',
                           todo_pattern='maps/depth-des-fake-%s.fits.fz',
                           outfn='tiles-depths-fake.fits')


    fn = 'p0-exposures.fits'
    T = fits_table('des-tiles.fits')
    T.writeto(fn)
    
    #des_zooms(done_pattern='maps/depth-decals-all-%s.fits.fz',
    #          todo_pattern='maps/depth-des-fake-%s.fits.fz')

    des_zooms(tiles_fn = 'tiles-depths-fake.fits',
              done_pattern='maps/depth-decals-all-%s.fits.fz',
              todo_pattern='maps/depth-des-fake-%s.fits.fz')

    #update_decam_tiles()
    #des_zooms()
    #djs_update()
    #plot_tiles()
    #when_missing()
    #tiles_todo()
    #needed_tiles()
    #main()

    # tiles = fits_table('mosaic-tiles_remaining.fits')
    # ps = PlotSequence('remaining')
    # plot_tiles(tiles, ps)
    
