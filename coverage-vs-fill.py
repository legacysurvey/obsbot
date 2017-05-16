from __future__ import print_function
from astrometry.util.fits import fits_table
import pylab as plt
import numpy as np
import fitsio
#import time
from astrometry.util.util import Tan, wcs_pv2sip_hdr
from astrometry.libkd.spherematch import match_radec
from astrometry.util.resample import resample_with_wcs, OverlapError
from astrometry.util.plotutils import PlotSequence
#from astrometry.util.starutil_numpy import *

from legacypipe.survey import LegacySurveyData
from collections import Counter

'''
A script to look into our current depth vs fill factor, to ask, "could
we retire any planned tiles", by, eg, making depth vs fill factor plots.

'''

if __name__ == '__main__':
    ps = PlotSequence('covfill')

    #    export LEGACY_SURVEY_DIR=~/legacypipe-dir-mzls/
    survey = LegacySurveyData()
    ccds = survey.get_annotated_ccds()
    print('Annotated CCDs:', len(ccds))
    print('Unique exposures:', len(np.unique(ccds.expnum)))

    ccds.cut(ccds.exptime > 60)
    print('Exptime > 60 sec:', len(ccds))
    
    nccds = Counter(ccds.expnum)
    for k,v in nccds.most_common():
        if v <= 4:
            break
        print('Expnum', k, 'appears', v, 'times')
    print('Tile pass numbers:', Counter(ccds.tilepass).most_common())

    # Fix parsing of OBJECT field to tileid...
    from obsbot import get_tile_id_from_name
    tileids = []
    for o in ccds.object:
        tid = get_tile_id_from_name(o.strip())
        if tid is None:
            tid = 0
        tileids.append(tid)
    tileids = np.array(tileids)
    print(len(np.unique(tileids)), 'unique tile ids in annotated file, from OBJECT')
    print(len(np.unique(ccds.tileid)), 'unique tile ids in ann file from TILEID')
    D = np.flatnonzero(tileids != ccds.tileid)
    print(len(D), 'different tileids')
    print('From OBJECT:', tileids[D])
    print('From TILEID:', ccds.tileid[D])
    ccds.tileid = tileids

    T = fits_table('obstatus/mosaic-tiles_obstatus.fits')
    
    # Update ccds.tilepass from ccds.tileid
    tileidtopass = dict(zip(T.tileid, T.get('pass')))
    tileidtoebv = dict(zip(T.tileid, T.ebv_med))
    ccds.tilepass = np.array([tileidtopass.get(tileid, 0)
                              for tileid in ccds.tileid])
    ccds.tileebv = np.array([tileidtoebv.get(tileid, 0)
                             for tileid in ccds.tileid])
    print('Tile pass numbers after update:', Counter(ccds.tilepass).most_common())
    e,I = np.unique(ccds.expnum, return_index=True)
    exps = ccds[I]
    #print('Object names,exptimes for tilepass==0:', zip(exps.object[exps.tilepass == 0], exps.exptime[exps.tilepass == 0]))

    # Average the depth per exposure
    for j,expnum in enumerate(exps.expnum):
        I = np.flatnonzero(ccds.expnum == expnum)
        if len(I) != 4:
            print('Exposure', expnum, 'has only', len(I), 'CCD entries')
            continue
        # Don't include zeros in computing average depths!
        Igood = I[(ccds.galdepth[I] > 0) * (ccds.ccdzpt[I] < 30)]
        if len(Igood) > 0:
            exps.galdepth[j] = np.mean(ccds.galdepth[Igood])
        else:
            exps.galdepth[j] = 0.
    
    # CCDs-table-based mapping from tileid to depth.
    I = np.flatnonzero((exps.tilepass > 0) * (exps.galdepth > 0) *
                       (exps.tileid > 0))
    tileid_to_depth = dict(zip(exps.tileid[I], exps.galdepth[I]))

    
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

    retirable = []

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
    
    for itile,tile in enumerate(todo):
        print()
        print('Tile', itile+1, 'of', len(todo), ':', tile.tileid, 'at', tile.ra, tile.dec)

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

        pixscale = 2.
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

        matched_exposures = set()
        #print('  tileids', T.tileid[I])
        #print('  expnums', T.z_expnum[I])

        for ii in I:
            for twcs in tilewcs:
                twcs.set_crval((T.ra[ii], T.dec[ii]))

                # Get depth from CCDs file, if available.
                zdepth = tileid_to_depth.get(T.tileid[ii], 0.)
                if zdepth == 0.:
                    zdepth = T.z_depth[ii]

                matched_exposures.add(T.z_expnum[ii])
                    
                try:
                    Yo,Xo,Yi,Xi,rims = resample_with_wcs(thiswcs, twcs)
                except OverlapError:
                    continue

                dflux = 10.**((zdepth - 22.5)/-2.5)
                div = 1./dflux**2
                depth[Yo,Xo] += div
                nexp[Yo,Xo] += 1


        # Now also look for entries in the CCDs (exposures) table not previously found.
        I,J,d = match_radec(exps.ra_bore, exps.dec_bore, tile.ra, tile.dec, radius)
        print(' ', len(I), 'exposures from CCDs file nearby')
        if len(I):
            I = np.array([i for i,expnum,gd in zip(I, exps.expnum[I], exps.galdepth[I])
                          if (not expnum in matched_exposures) and gd > 0])
            print(' ', len(I), 'exposures that were not in tile file')
        # Drop exposures from this pass, except for previous exposures of this tile!
        if len(I):
            I = I[np.logical_or(exps.tilepass[I] != 3,
                                exps.tileid[I] == tile.tileid)]
            print(' ', len(I), 'exposures not in pass 3')
        if len(I):
            # print('  objects:', [o.strip() for o in exps.object[I]])
            # print('  tileids:', exps.tileid[I])
            # print('  expnums:', exps.expnum[I])
            # print('  passes:', exps.tilepass[I])

            for ii in I:
                for twcs in tilewcs:
                    twcs.set_crval((exps.ra_bore[ii], exps.dec_bore[ii]))

                    # Get depth from CCDs file, if available.
                    zdepth = exps.galdepth[ii]

                    try:
                        Yo,Xo,Yi,Xi,rims = resample_with_wcs(thiswcs, twcs)
                    except OverlapError:
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

        print('  Depths at 2, 5, and 10th percentile vs target:',
              '%.2f' % (depths[2]  - (target - 0.6)),
              '%.2f' % (depths[5]  - (target - 0.3)),
              '%.2f' % (depths[10] -  target))

        if not ((depths[2]  > target - 0.6) and
                (depths[5]  > target - 0.3) and
                (depths[10] > target)):
            continue

        retirable.append((tile.tileid, depths))

        #if len(retirable) == 10:
        #    break
        if ps.ploti >= 100:
            continue
        
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

        #if ps.ploti == 100:
        #    break
        
    print('Tiles that could be retired:')
    print('# Tileid 0th-percentile-extcorr-depth 1st-pctile 2nd-pctile ...')
    for tileid, depths in retirable:
        print(tileid, ' '.join(['%.3f' % d for d in depths]))
