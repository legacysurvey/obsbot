from __future__ import print_function
from collections import Counter

import matplotlib
import pylab as plt

from legacypipe.common import LegacySurveyData

from astrometry.util.fits import *

from decam import DecamNominalCalibration

def main():
    survey = LegacySurveyData()
    ccds = survey.get_annotated_ccds()
    print(len(ccds), 'CCDs')
    ccds.cut(ccds.tileid > 0)
    print(len(ccds), 'with tileid')
    ccds.cut(ccds.photometric)
    print(len(ccds), 'photometric')

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
    
    O = fits_table('obstatus/decam-tiles_obstatus.fits')
    print(len(O), 'tiles')
    # maxdec = 34.
    # O.cut((O.in_desi == 1) * (O.in_des == 0) *
    #       (O.dec > -20) * (O.dec < maxdec) * (O.get('pass') <= 3))
    # print(len(O), 'tiles

    # "tileid" = row number (1-indexed)
    assert(np.all(O.tileid == np.arange(1, len(O)+1)))

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
    detsig = 10.**((E.galdepth - 22.5) / -2.5)
    detiv  = 1. / detsig**2
    
    nom = DecamNominalCalibration()
    
    for band in 'grz':
        I = np.flatnonzero(E.filter == band)
        print(len(I), 'exposures in', band)
        iv = np.zeros(len(O), np.float32)
        np.add.at(iv, E.tileid[I] - 1, detiv[I])
        print('galdepth range:', E.galdepth[I].min(), E.galdepth[I].max())
        # convert iv back to galdepth in mags
        with np.errstate(divide='ignore'):
            galdepth = -2.5 * (np.log10(np.sqrt(1. / iv)) - 9)
        # extinction correction...

        fid = nom.fiducial_exptime(band)
        extinction = O.ebv_med * fid.A_co
        print('Median extinction:', np.median(extinction))
        galdepth -= extinction

        galdepth[iv == 0] = 0.
        
        O.set('%s_depth' % band, galdepth)

        from astrometry.libkd.spherematch import match_radec
        from astrometry.util.plotutils import antigray
        rlo,rhi = 0,360
        dlo,dhi = -20,35
        rr,dd = np.meshgrid(np.linspace(rlo,rhi,720), np.linspace(dlo,dhi,360))
        JJ,II,d = match_radec(rr.ravel(), dd.ravel(), O.ra, O.dec, 1.5,
                              nearest=True)
        indesi = np.zeros(rr.shape, bool)
        indesi.flat[JJ] = ((O.in_desi[II] == 1) * (O.in_des[II] == 0))

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

