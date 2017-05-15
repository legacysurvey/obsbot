from __future__ import print_function
from astrometry.util.fits import fits_table
import pylab as plt
import numpy as np
import fitsio
from astrometry.util.util import Tan, wcs_pv2sip_hdr
from astrometry.libkd.spherematch import match_radec
from astrometry.util.resample import resample_with_wcs
from astrometry.util.plotutils import PlotSequence
#from astrometry.util.starutil_numpy import *

'''
A script to look into our current depth vs fill factor, to ask, "could
we retire any planned tiles", by, eg, making depth vs fill factor plots.

'''

if __name__ == '__main__':
    ps = PlotSequence('covfill')

    T = fits_table('obstatus/mosaic-tiles_obstatus.fits')
    T.cut(T.in_desi == 1)
    T.cut(T.get('pass') <= 3)
    
    # The tiles we'll examine
    P3 = T[T.get('pass') == 3]
    print(len(P3), 'pass 3, in DESI')
    todo = P3[P3.z_done == 0]
    print(len(todo), 'pass 3 tiles to do')

    # Tiles with measured depths
    T.cut((T.z_depth > 15) * (T.z_depth < 30))
    print(len(T), 'tiles with measured depths')
    # Passes other than 3... they ~ only barely overlap, so don't
    # contribute significant depth.
    T.cut(T.get('pass') < 3)
    
    udecs = np.unique(P3.dec)
    print(len(udecs), 'unique Dec values in pass 3')

    # Grab an arbitrary weight-map image and use that as a proxy!
    wtfn = 'k4m_170501_112501_oow_zd_v1.fits.fz'
    F = fitsio.FITS(wtfn)
    tilewcs = []
    #hdr = F[0].read_header()
    #tilera,tiledec = hdr['CENTRA'], hdr['CENTDEC']
    #wcsoffsets = []
    for i in range(1, len(F)):
        hdr = F[i].read_header()
        wcs = wcs_pv2sip_hdr(hdr)
        tilewcs.append(wcs)
        #cra,cdec = wcs.get_crval()
        #wcsoffsets.append(((cra - tilera) * np.cos(tiledec),
        #                   cdec - tiledec))
        #print('WCS header', i, 'has offset', wcsoffsets[-1], 'deg')
        #print('WCS bounds', wcs.radec_bounds())
        ### They make the CRVAL be the boresight for all chips... perfect!
    
    for tile in todo:
        print()
        print('Tile', tile.tileid, 'at', tile.ra, tile.dec)

        i = np.nonzero(tile.dec == udecs)[0][0]
        if i == 0 or i == len(udecs)-1:
            print('Endpoint Dec; skipping for now')
            continue
        print('  Decs:', udecs[i-1], tile.dec, udecs[i+1])

        declo = (udecs[i-1] + tile.dec) / 2.
        dechi = (udecs[i+1] + tile.dec) / 2.
        
        row = P3[P3.dec == tile.dec]
        print(' ', len(row), 'tiles in this Dec row')
        ras = np.sort(row.ra)
        i = np.nonzero(tile.ra == ras)[0][0]
        if i == 0 or i == len(ras)-1:
            print('  Endpoint RA; skipping for now')
            continue
        print('  RAs:', ras[i-1], tile.ra, ras[i+1])

        ralo = (ras[i-1] + tile.ra) / 2.
        rahi = (ras[i+1] + tile.ra) / 2.

        pixscale = 1.
        #pixscale = 0.262

        H = int(np.ceil((dechi - declo) / (pixscale/3600.)))
        W = int(np.ceil((rahi - ralo) * np.cos(np.deg2rad(tile.dec)) / (pixscale/3600.)))
        print('  Dec height', dechi-declo, 'RA width', rahi-ralo, '-> pix', W, 'x', H)

        cd = pixscale/3600.
        thiswcs = Tan(tile.ra, tile.dec, (W+1)/2., (H+1)/2.,
                      -cd, 0., 0., cd, float(W), float(H))

        # Find surrounding tiles
        radius = np.hypot(W, H) * pixscale / 3600.
        
        I,J,d = match_radec(T.ra, T.dec, tile.ra, tile.dec, radius)
        print(' ', len(I), 'tiles with measured depths nearby')

        if len(I) == 0:
            continue
        
        # Now we need to take a tile boresight position and map it to
        # boxes in RA,Dec of the good portions of the CCDs.
        
        depth = np.zeros((H,W), np.float32)
        nexp = np.zeros((H,W), np.uint8)

        for ii in I:
            for twcs in tilewcs:
                twcs.set_crval((T.ra[ii], T.dec[ii]))
                zdepth = T.z_depth[ii]

                try:
                    Yo,Xo,Yi,Xi,rims = resample_with_wcs(thiswcs, twcs)
                except:
                    continue

                dflux = 10.**((zdepth - 22.5)/-2.5)
                div = 1./dflux**2
                depth[Yo,Xo] += div
                nexp[Yo,Xo] += 1

        # Convert depth map from depth-iv back to mag.
        # flux
        with np.errstate(divide='ignore'):
            dflux = np.sqrt(1./depth)
            dflux[depth == 0] = 0.
            depth = -2.5 * (np.log10(dflux) - 9.)
        depth[dflux == 0] = 0.

        if depth.max() == 0:
            print('  Actually no overlap')
            continue

        # Extinction correction for this tile...
        ext_z = 1.211 * tile.ebv_med
        print('  Applying extinction correction', ext_z, 'mag')
        depth[depth != 0] -= ext_z
        
        #pcts = [0,10,20,30,40,50,60,70,80,90,100]
        pcts = np.arange(0, 101)
        depths = np.percentile(depth, pcts)

        target = 22.5
        req_pcts = [0, 2, 2, 5, 5, 10, 10, 100]
        req_depths = [0, 0, target-0.6, target-0.6,
                      target-0.3, target-0.3, target, target]
        
        maglo, maghi = 21,23
        
        plt.clf()
        #plt.subplot(1,2,1)
        plt.subplot2grid((2,2), (0,0))
        plt.imshow(depth, interpolation='nearest', origin='lower',
                   vmin=maglo, vmax=maghi)
        plt.colorbar(ticks=[np.arange(maglo, maghi+0.01, 0.5)])
        plt.xticks([]); plt.yticks([])
        #plt.subplot(1,2,2)
        plt.subplot2grid((2,2), (1,0))
        plt.imshow(nexp, interpolation='nearest', origin='lower',
                   vmin=0, vmax=4)
        plt.colorbar(ticks=[0,1,2,3,4])
        plt.xticks([]); plt.yticks([])

        ax = plt.subplot2grid((2,2), (0,1), rowspan=2)
        plt.plot(req_pcts, np.clip(req_depths, maglo, maghi), 'k-',
                 lw=2, alpha=0.5)
        plt.plot(pcts, np.clip(depths, maglo, maghi), 'b-')
        ax.yaxis.tick_right()
        plt.ylim(maglo, maghi)
        plt.xlim(0, 100)
        plt.xlabel('Coverage fraction')
        plt.ylabel('Existing depth')
        plt.suptitle('Tile %i' % tile.tileid)
        ps.savefig()

        if ps.ploti == 100:
            break
        
