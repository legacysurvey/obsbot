from __future__ import print_function
from collections import Counter

import matplotlib
import pylab as plt

from legacypipe.common import LegacySurveyData

from astrometry.util.fits import *
from astrometry.libkd.spherematch import match_radec

from decam import DecamNominalCalibration

def main():
    survey = LegacySurveyData()
    ccds = survey.get_annotated_ccds()
    print(len(ccds), 'CCDs')

    O = fits_table('obstatus/decam-tiles_obstatus.fits')
    print(len(O), 'tiles')
    # "tileid" = row number (1-indexed)
    assert(np.all(O.tileid == np.arange(1, len(O)+1)))

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

    ccds.cut(ccds.tileid > 0)
    print(len(ccds), 'with tileid')

    expnums,I = np.unique(ccds.expnum, return_index=True)
    print(len(expnums), 'unique exposures')
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
    for band in 'grz':
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

    # actually 5*detsig1...
    with np.errstate(divide='ignore', over='ignore'):
        detsig = 10.**((E.galdepth - 22.5) / -2.5)
        detiv  = 1. / detsig**2
    
    nom = DecamNominalCalibration()
    
    for band in 'grz':
        I = np.flatnonzero((E.filter == band) * E.photometric)
        print(len(I), 'photometric exposures in', band)
        iv = np.zeros(len(O), np.float32)
        np.add.at(iv, E.tileid[I] - 1, detiv[I])
        print('galdepth range:', E.galdepth[I].min(), E.galdepth[I].max())
        # convert iv back to galdepth in mags
        with np.errstate(divide='ignore'):
            galdepth = -2.5 * (np.log10(np.sqrt(1. / iv)) - 9)
        # extinction correction...

        fid = nom.fiducial_exptime(band)
        extinction = O.ebv_med * fid.A_co
        galdepth -= extinction

        galdepth[iv == 0] = 0.

        # Flag non-photometric exposures with depth = 1.
        I = np.flatnonzero((E.filter == band) * (E.photometric == False))
        assert(np.all(galdepth[E.tileid[I] - 1]) == 0.)
        galdepth[E.tileid[I] - 1] = 1.
        print('Marking', len(I), 'non-photometric exposures in', band)

        O.set('%s_depth' % band, galdepth)

        from astrometry.util.plotutils import antigray
        rlo,rhi = 0,360
        dlo,dhi = -20,35
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
            plt.clf()

            J = np.flatnonzero((O.in_desi == 1) * (O.in_des == 0) *
                               (O.dec > dlo) * (O.dec < dhi) *
                               (O.get('pass') == passnum))
            #plt.plot(O.ra[J], O.dec[J], 'k.', alpha=0.5)
            plt.imshow(indesi, extent=[rlo,rhi,dlo,dhi], vmin=0, vmax=4,
                       cmap=antigray, aspect='auto', interpolation='nearest',
                       origin='lower')
            depth = O.get('%s_depth' % band)
            J = np.flatnonzero((O.get('pass') == passnum) * (depth > 0))

            target = fid.single_exposure_depth
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
            J = np.flatnonzero((O.get('pass') == passnum) * (depth > 0))
            depth = depth[J]
            mlo,mhi = 21,24
            plt.hist(np.clip(depth, mlo,mhi), bins=100, range=(mlo,mhi),
                     histtype='step', color=' bgr'[passnum],
                     label='Pass %i' % passnum)
        plt.axvline(fid.single_exposure_depth, color='k')
        plt.xlabel('Depth (mag)')
        plt.legend(loc='upper left')
        plt.title('DECaLS depth: %s' % band)
        plt.savefig('depth-%s.png' % band)

        
    O.writeto('decam-obstatus-depth.fits')

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

