from __future__ import print_function
from collections import Counter
import sys

import matplotlib
import pylab as plt

from legacypipe.common import LegacySurveyData

from astrometry.util.fits import *
from astrometry.libkd.spherematch import match_radec
from astrometry.util.starutil_numpy import degrees_between

def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--mzls', action='store_true',
                        help='Set MzLS (default: DECaLS)')
    opt = parser.parse_args()

    if opt.mzls:
        from mosaic import MosaicNominalCalibration
        nom = MosaicNominalCalibration()

        obstatus_fn = 'obstatus/mosaic-tiles_obstatus.fits'
        out_fn = 'mosaic-obstatus-depth.fits'
        
        bands = 'z'

        declo,dechi = 30,80

    else:
        from decam import DecamNominalCalibration
        nom = DecamNominalCalibration()

        obstatus_fn = 'obstatus/decam-tiles_obstatus.fits'
        out_fn = 'decam-obstatus-depth.fits'

        bands = 'grz'

        declo,dechi = -20,35
        
    survey = LegacySurveyData()
    ccds = survey.get_annotated_ccds()
    print(len(ccds), 'CCDs')

    O = fits_table(obstatus_fn)
    print(len(O), 'tiles')

    # Map tile IDs back to index in the obstatus file.
    tileid_to_index = np.empty(max(O.tileid)+1, int)
    tileid_to_index[:] = -1
    tileid_to_index[O.tileid] = np.arange(len(O))

    assert(len(np.unique(O.tileid)) == len(O))
    
    I = tileid_to_index[O.tileid]
    assert(np.all(I == np.arange(len(O))))

    print('Pass numbers:', np.unique(O.get('pass')))
    
    goodtiles = np.flatnonzero(O.in_desi * (O.dec > 30) * (O.get('pass') <= 3))
    print(len(goodtiles), 'tiles of interest')
    
    # Look at whether exposures from other programs are near our tile centers.
    # Basically nope.
    # plt.clf()
    # e,K = np.unique(ccds.expnum, return_index=True)
    # I,J,d = match_radec(O.ra, O.dec, ccds.ra_bore[K], ccds.dec_bore[K],
    #                     1./60., nearest=True)
    # KK = K[np.flatnonzero(ccds.tileid[K] > 0)]
    # I,J,d2 = match_radec(O.ra, O.dec, ccds.ra_bore[KK], ccds.dec_bore[KK],
    #                      1./60., nearest=True)
    # ha = dict(range=(0., 60.), bins=60, histtype='step')
    # plt.hist(d * 3600., color='b', **ha)
    # plt.hist(d2 * 3600., color='r', **ha)
    # plt.xlabel('Distance from tile to nearest DECam boresight (arcsec)')
    # plt.savefig('dists.png')

    notileids = ccds[ccds.tileid <= 0]
    print(len(notileids), 'CCDs have no tileid')
    I,J,d = match_radec(notileids.ra_bore, notileids.dec_bore, O.ra, O.dec,
                        0.5, nearest=True)
                        
    plt.clf()
    plt.hist(d, bins=50)
    plt.xlabel('Distance to nearest tile center (deg)')
    plt.savefig('tiledist.png')

    plt.clf()
    plt.hist(d*3600, bins=50, range=(0,30))
    plt.xlabel('Distance to nearest tile center (arcsec)')
    plt.savefig('tiledist2.png')
    
    ccds.cut(ccds.tileid > 0)
    print(len(ccds), 'CCDs with tileid')

    expnums,I = np.unique(ccds.expnum, return_index=True)
    print(len(expnums), 'unique exposures (with tileids)')
    # Compute the mean depth per exposure
    E = ccds[I]
    for expnum in expnums:
        I = np.flatnonzero(ccds.expnum == expnum)
        j = np.flatnonzero(E.expnum == expnum)
        assert(len(j) == 1)
        j = j[0]
        E.galdepth[j] = np.mean(ccds.galdepth[I])

        E.photometric[j] = np.all(ccds.photometric[I])
        if len(np.unique(ccds.photometric[I])) == 2:
            print('Exposure', expnum, 'has photometric and non- CCDs')
    
    # plt.clf()
    # plt.plot(O.ra, O.dec, 'k.')
    # plt.axis([360,0,-25,35])
    # plt.title('All tiles')
    # plt.savefig('tiles-all.png')
    # 
    # print('in_desi:', np.unique(O.in_desi))
    # plt.clf()
    # J = np.flatnonzero(O.in_desi == 1)
    # plt.plot(O.ra[J], O.dec[J], 'k.')
    # plt.axis([360,0,-25,35])
    # plt.title('In DESI')
    # plt.savefig('tiles-desi.png')
    # 
    # print('in_des:', np.unique(O.in_des))
    # plt.clf()
    # J = np.flatnonzero(O.in_des == 1)
    # plt.plot(O.ra[J], O.dec[J], 'k.')
    # plt.axis([360,0,-25,35])
    # plt.title('IN DES')
    # plt.savefig('tiles-des.png')
    
    #print('Number of exposures of each tile:')
    #print(Counter(E.tileid).most_common())
    print('Number of exposures of tiles:')
    for band in bands:
        I = np.flatnonzero(E.filter == band)
        c = Counter(E.tileid[I])
        c2 = Counter([v for k,v in c.most_common()])
        print('  ', band, 'band:', c2.most_common())
    
    # Detection inverse-variance is the quantity that adds when there are
    # multiple exposures.
    #   detsig1 = ccds.sig1 / ccds.galnorm_mean
    #   depth = 5. * detsig1
    #   # that's flux in nanomaggies -- convert to mag
    #   ccds.galdepth = -2.5 * (np.log10(depth) - 9)

    with np.errstate(divide='ignore', over='ignore'):
        # actually 5*detsig1...
        detsig = 10.**((E.galdepth - 22.5) / -2.5)
        E.detiv  = 1. / detsig**2

    for band in bands:
        # "I" indexes into exposures E.
        I = np.flatnonzero((E.filter == band) * E.photometric *
                           np.isfinite(E.detiv))
        print(len(I), 'photometric exposures in', band)

        # "iv" is for O.
        iv = np.zeros(len(O), np.float32)
        # "J" indexes into obstatus tiles O.
        J = tileid_to_index[E.tileid[I]]
        assert(np.all((J >= 0) * (J < len(O))))
        assert(np.all(O.tileid[J] == E.tileid[I]))

        print('tileid range', E.tileid[I].min(), E.tileid[I].max())

        d = np.array([degrees_between(*a) for a in
                      zip(E.ra_bore[I], E.dec_bore[I], O.ra[J], O.dec[J])])
        print('Degrees between tiles & exposures:', d)
        
        np.add.at(iv, J, E.detiv[I])
        print('galdepth range:', E.galdepth[I].min(), E.galdepth[I].max())
        print('detiv range:', E.detiv[I].min(), E.detiv[I].max())
        print('index range:', J.min(), J.max())
        nexp = np.zeros(len(O), int)
        np.add.at(nexp, J, 1)
        print('tile exposure counts:', Counter(nexp))

        # convert iv back to galdepth in mags
        with np.errstate(divide='ignore'):
            galdepth = -2.5 * (np.log10(np.sqrt(1. / iv)) - 9)
        # extinction correction...
        #print('galdepth deciles:', np.percentile(galdepth, [0,10,20,30,40,50,60,70,80,90,100]))

        fid = nom.fiducial_exptime(band)
        extinction = O.ebv_med * fid.A_co
        #print('Extinction range:', extinction.min(), extinction.max())
        galdepth -= extinction

        galdepth[iv == 0] = 0.

        print('galdepth deciles:', np.percentile(galdepth, [0,10,20,30,40,50,60,70,80,90,100]))
        
        # Flag tiles that have *only* non-photometric exposures with depth = 1.
        I = np.flatnonzero((E.filter == band) * (E.photometric == False))
        print(len(I), 'exposures are non-photometric in', band, 'band')
        only_nonphot = np.flatnonzero(
            galdepth[tileid_to_index[E.tileid[I]]] == 0.)
        print(len(only_nonphot), 'tiles have only non-photometric exposures')
        galdepth[only_nonphot] = 1.
        print('Marking', len(only_nonphot),'non-photometric exposures in', band)

        O.set('%s_depth' % band, galdepth)

        print('Depth deciles:', np.percentile(O.get('%s_depth' % band), [0,10,20,30,40,50,60,70,80,90,100]))

        from astrometry.util.plotutils import antigray
        rlo,rhi = 0,360
        dlo,dhi = declo,dechi
        rr,dd = np.meshgrid(np.linspace(rlo,rhi,720), np.linspace(dlo,dhi,360))
        JJ,II,d = match_radec(rr.ravel(), dd.ravel(), O.ra, O.dec, 1.5,
                              nearest=True)
        indesi = np.zeros(rr.shape, bool)
        indesi.flat[JJ] = ((O.in_desi[II] == 1) * (O.in_des[II] == 0))

        J = np.flatnonzero((O.in_desi == 1) * (O.in_des == 0) *
                           (O.dec > dlo) * (O.dec < dhi))
        print('Median E(B-V) in DECaLS area:', np.median(O.ebv_med[J]))
        print('Median extinction in DECaLS area, %s band:' % band, np.median(extinction[J]))

        I = np.flatnonzero((O.get('%s_expnum' % band) > 0) * (O.get('%s_depth' % band) == 0))
        print('Found', len(I), 'tiles with', band, 'EXPNUM but no DEPTH; setting to DEPTH=30')
        O.get('%s_depth' % band)[I] = 30.

        # print('Exposure numbers:', O.get('%s_expnum' % band)[I])
        # print('Dates:', O.get('%s_date' % band)[I])
        
        plt.figure(figsize=(10,6))
        plt.subplots_adjust(left=0.1, right=0.9)
        
        for passnum in [1,2,3]:
            print('Pass', passnum)
            plt.clf()

            J = np.flatnonzero((O.in_desi == 1) * (O.in_des == 0) *
                               (O.dec > dlo) * (O.dec < dhi) *
                               (O.get('pass') == passnum))
            #plt.plot(O.ra[J], O.dec[J], 'k.', alpha=0.5)

            # Plot the gray background showing the in_desi footprint
            plt.imshow(indesi, extent=[rlo,rhi,dlo,dhi], vmin=0, vmax=4,
                       cmap=antigray, aspect='auto', interpolation='nearest',
                       origin='lower')
            depth = O.get('%s_depth' % band)
            #J = np.flatnonzero((O.get('pass') == passnum) * (depth > 0))
            
            J = np.flatnonzero((O.get('pass') == passnum) * (depth > 1) * (depth < 30))
            # print('Depths:', depth[J])
            print('Depth deciles:', np.percentile(depth[J], [0,10,20,30,40,50,60,70,80,90,100]))

            if len(J) == 0:
                sys.exit(0)
            
            target = fid.single_exposure_depth
            print('Target depth:', target)

            cmap = cmap_discretize('RdBu', 11)
            dm = 0.275
            
            plt.scatter(O.ra[J], O.dec[J], c=depth[J] - target, linewidths=0,
                        cmap=cmap, vmin=-dm, vmax=+dm, zorder=-10, s=1)
            plt.colorbar(ticks=np.arange(-0.25, 0.251, 0.05))

            hh,ww = rr.shape
            rgba = np.zeros((hh,ww,4), np.float32)
            JJ,II,d = match_radec(rr.ravel(), dd.ravel(), O.ra[J], O.dec[J], 1.,
                                  nearest=True)
            Jy,Jx = np.unravel_index(JJ, rr.shape)
            rgba[Jy,Jx,:] = cmap((np.clip(depth[J[II]] - target, -dm, dm) - (-dm)) / (dm - (-dm)))
            plt.imshow(rgba, extent=[rlo,rhi,dlo,dhi], aspect='auto', interpolation='nearest',
                       origin='lower')

            I = np.flatnonzero((depth == 0) * (O.get('%s_done' % band) == 1) *
                               (O.get('pass') == passnum))
            plt.plot(O.ra[I], O.dec[I], 'g.')
            
            plt.title('Band %s, Pass %i' % (band, passnum))
            plt.xlabel('RA (deg)')
            plt.ylabel('Dec (deg)')
            plt.axis([rhi,rlo,dlo,dhi])
            plt.savefig('depth-%s-%i.png' % (band, passnum))

        plt.clf()
        for passnum in [1,2,3]:
            depth = O.get('%s_depth' % band)
            J = np.flatnonzero((O.get('pass') == passnum) *
                               (depth > 1) * (depth < 30))
            depth = depth[J]

            print('Pass', passnum)
            print('Fiducial single-exposure-depth:', fid.single_exposure_depth)
            print(sum(depth < fid.single_exposure_depth - 0.25), 'of', len(depth), 'tiles are more than 0.25 mag shallow')

            odepth = O.get('%s_depth' % band)
            K = np.flatnonzero((O.get('%s_done' % band) == 0) * (O.get('pass') == passnum) *
                               (odepth > 1) * (odepth < 30))
            print(sum(odepth[K] < fid.single_exposure_depth - 0.25), 'of', len(odepth[K]), 'DONE=0 tiles are more than 0.25 mag shallow')
            K = np.flatnonzero((O.get('%s_done' % band) == 1) * (O.get('pass') == passnum) *
                               (odepth > 1) * (odepth < 30))
            print(sum(odepth[K] < fid.single_exposure_depth - 0.25), 'of', len(odepth[K]), 'DONE=1 tiles are more than 0.25 mag shallow')
            
            mlo,mhi = 21,24
            plt.hist(np.clip(depth, mlo,mhi), bins=100, range=(mlo,mhi),
                     histtype='step', color=' bgr'[passnum],
                     label='Pass %i' % passnum)
        plt.axvline(fid.single_exposure_depth, color='k')
        plt.axvline(fid.single_exposure_depth - 0.25, color='k', linestyle='--')
        plt.xlabel('Depth (mag)')
        plt.legend(loc='upper left')
        plt.title('Depth: %s' % band)
        plt.savefig('depth-%s.png' % band)
        
    O.writeto(out_fn)

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
    main()

