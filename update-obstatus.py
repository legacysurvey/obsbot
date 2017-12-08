from __future__ import print_function
from collections import Counter
import sys

import matplotlib
import pylab as plt
import numpy as np

from legacypipe.survey import LegacySurveyData

from astrometry.util.fits import *
from astrometry.libkd.spherematch import match_radec
from astrometry.util.starutil_numpy import degrees_between

def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--mzls', action='store_true',
                        help='Set MzLS (default: DECaLS)')
    parser.add_argument('--ann', help='Set annotated-CCDs file')
    opt = parser.parse_args()

    if opt.mzls:
        from mosaic import MosaicNominalCalibration
        from camera_mosaic import database_filename, camera_name
        nom = MosaicNominalCalibration()

        obstatus_fn = 'obstatus/mosaic-tiles_obstatus.fits'
        out_fn = 'mosaic-obstatus-depth.fits'

        bands = 'z'

        declo,dechi = -5,90

        f = open('obstatus/bad_expid.txt')
        bad_expids = set()
        for line in f:
            line = line.strip()
            if len(line) == 0:
                continue
            if line[0] == '#':
                continue
            words = line.split()
            try:
                expnum = int(words[0])
            except:
                print('Skipping line:', line)
                continue
            bad_expids.add(expnum)
        print('Read', len(bad_expids), 'bad exposure numbers')

    else:
        from decam import DecamNominalCalibration
        from camera_decam import database_filename, camera_name
        nom = DecamNominalCalibration()

        obstatus_fn = 'obstatus/decam-tiles_obstatus.fits'
        out_fn = 'decam-obstatus-depth.fits'

        bands = 'grz'

        declo,dechi = -20,35

    # Convert copilot db to fits.
    import obsdb
    from copilot import db_to_fits
    obsdb.django_setup(database_filename=database_filename)
    ccds = obsdb.MeasuredCCD.objects.all()
    copilot = db_to_fits(ccds)

    all_copilot = copilot.copy()

    copilot.writeto('copilot.fits')
    
    print(len(copilot), 'measured CCDs in copilot database')
    copilot.cut(np.array([c.strip() == camera_name for c in copilot.camera]))
    print(len(copilot), 'copilot CCDs with camera = "%s"' % camera_name)
    copilot.cut(copilot.expnum > 0)
    print(len(copilot), 'measured CCDs in copilot database with EXPNUM')

    print('Copilot expfactor extremes:', np.percentile(copilot.expfactor[copilot.expfactor != 0], [1,99]))

    survey = LegacySurveyData()
    if opt.ann:
        ccds = fits_table(opt.ann)
    else:
        print('Reading annotated CCDs files...')
        ccds = survey.get_annotated_ccds()
    print(len(ccds), 'CCDs')

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

    O = fits_table(obstatus_fn)
    print(len(O), 'tiles')

    if opt.mzls:
        from camera_mosaic import fix_expnums
        # Fix MzLS exposure numbers with wrong leading "3".
        fix_expnums(ccds.expnum)
        # Also fix leading "3" in expnums in OBSTATUS file
        fix_expnums(O.z_expnum)
        # And copilot database
        fix_expnums(copilot.expnum)

        print('Z_EXPNUM range:', O.z_expnum.min(), 'min >0:', O.z_expnum[O.z_expnum > 0].min(), O.z_expnum.max())

    print('Pass numbers:', np.unique(O.get('pass')))
    
    if opt.mzls:
        goodtiles = (O.in_desi * (O.dec > 30) * (O.get('pass') <= 3))
        print(sum(goodtiles), 'tiles of interest')
    else:
        goodtiles = (O.in_desi * (O.get('pass') <= 3))
        print(sum(goodtiles), 'tiles in the footprint')

    #O.cut(goodtiles)
    #print('Cut to', len(O), 'tiles of interest')

    # *after* fixing tileids
    allccds = ccds.copy()

    # Map tile IDs back to index in the obstatus file.
    tileid_to_index = np.empty(max(O.tileid)+1, int)
    tileid_to_index[:] = -1
    tileid_to_index[O.tileid] = np.arange(len(O))

    assert(len(np.unique(O.tileid)) == len(O))
    
    I = tileid_to_index[O.tileid]
    assert(np.all(I == np.arange(len(O))))

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

    ccds.photometric = (ccds.ccd_cuts == 0)

    # Compute the mean depth per exposure
    E = ccds[I]
    for expnum in expnums:
        I = np.flatnonzero(ccds.expnum == expnum)
        j = np.flatnonzero(E.expnum == expnum)
        assert(len(j) == 1)
        j = j[0]
        E.photometric[j] = np.all(ccds.photometric[I])
        #E.photometric[j] = np.all(ccds.ccd_cuts[I] == 0)
        if len(np.unique(ccds.photometric[I])) == 2:
            print('Exposure', expnum, 'has photometric and non-photometric CCDs')
            non = I[ccds.photometric[I] == False]
            phot = I[ccds.photometric[I]]

            if opt.mzls and len(phot) == 3:
                print('Accepting an exposure with 3 good CCDs')
                E.photometric[j] = True
                # And remove this exposure from the bad_expid list.
                if expnum in bad_expids:
                    bad_expids.remove(expnum)
                    print('Removing exposure', expnum, 'from bad_expid file')

                continue

            for ii in non:
                print('    http://legacysurvey.org/viewer-dev/?ra=%.3f&dec=%.3f&zoom=11&ccds3&bad=%i-%s' % (ccds.ra_center[ii], ccds.dec_center[ii], expnum, ccds.ccdname[ii]))
                print('    http://legacysurvey.org/viewer-dev/ccd/decals-dr5/decam-%s-%s-%s/' % (ccds.expnum[ii], ccds.ccdname[ii], ccds.filter[ii]))
            print('  image:', ccds.image_filename[I][0])
            print('  boresight:', ccds.ra_bore[I][0], ccds.dec_bore[I][0])
            #print('  ccdnames:', ccds.ccdname[I])
            print('  photometric:', len(phot), ', non-photometric:', len(non))
            print('  median phot depth:', np.median(ccds.galdepth[phot]))
            #print('  depth:', ccds.galdepth[I])
            print('  non-photometric CCDs:', ccds.ccdname[non])
            print('    depths:', ccds.galdepth[non])
            print('    ccdnmatch', ccds.ccdnmatch[non], 'vs', ccds.ccdnmatch[phot])
            print('    ccdtransp:', ccds.ccdtransp[non], 'vs', ccds.ccdtransp[phot])
            print('    ccd zpt vs frame zpt:', ccds.ccdzpt[non] - ccds.zpt[non])
            dp = ccds.ccdzpt[phot] - ccds.zpt[phot]
            print('      phot ccds zpt vs frame: range', dp.min(), dp.max(), 'mean', dp.mean())


            whitelist = [ 346662, 346664, 346665, # S3/S29 striping
                          346754, # one bad chip, wispy
                          346967, 347304, # M5 globular
                          347664, # zpt scatter
                          347744, # weird eye-shaped ghost; but lots of cov.
                          347755, 347768, 347769, 347782, # shallow, zpt scatter -- wispy pattern on focal plane
                          347918, 347920, # straddling transparency cut
                          347934, 347936, 347941, 347945, 347947, # zpt scatter
                          392377, 392380, 393173, # bright star
                          393671, # bright star
                          393672, 393673, # scatter
                          425339, 425340, # strangely low ccdnmatch
                          426225, 430808, # globular cluster
                          431640, # bright star
                          431644, # globular
                          432154, # one amp high bias
                          432159, # transp. on boundary
                          432179, # one amp high bias, +
                          432747, 432748, 432751, # scatter
                          433305, 433306, # bright star
                          497062, 497064, 497065, # low ccdnmatch
                          509516, 509517, # bright star
                          511247, 511263, # low ccdnmatch
                          511513, 511514, # bright star
                          512303, # bright star
                          520560, # bright star
                          521782, # scatter
                          522212, # bright star
                          535138, # bright stars, satellite hits?
                          535141, 535142, 535143, 535149, # low ccdnmatch
                          535210, # bright star
                          535695, # globular
                          536065, # globular
                          536385, # strangely zero ccdnmatch
                          547761, # nice galaxy
                          548257, # bright star
                          553779, # scatter 
                          553795, # shallow
                          554284, # marginal zpt
                          563659, # zpt scatter
                          563850, # mild striping
                          563852, # ??
                          583118, # bright stars?
                          592621, # marginal zpt, scatter
                          592859, # some pattern noise
                          605068, # ??
                          625710, # strangely low ccdnmatch
                          631005, # bright star
                          634493, # globular
                          634786, # strangely low ccdnmatch
                          634877, # globular
                          635535, # bright star, glob
                          635962, # low ccdnmatch
                          635973, # bias level?
                          636018, # bias level?
                          637688, # bright star
                          ]
            blacklist = [
                425328, # 2.7" seeing
                488244, 488256, 488260, 488261, 488263, 488268, # weird striping
                488270, # weird striping
                496913, 496917, 496918, 496919, 496920, 496921, 496922, # 3" seeing
                496923, 496925, 496926, 496927, 496928, 496930, # 3" seeing
                509162, 509163, 509166, 509172, 509176, 509182, 509202, # 3" seeing
                535471, # 4" seeing!
                535498, # 3" seeing
                548218, # double PSF -- telescope moved?
                563835, # striping
                563842, # striping
                ]

            if expnum in whitelist:
                print('** Exposure', expnum, 'in whitelist -- marking as photometric')
                E.photometric[j] = True
            
        # Don't include zeros in computing average depths!
        Igood = I[(ccds.galdepth[I] > 0) * (ccds.ccdzpt[I] < 30)]
        if len(Igood) > 0:
            E.galdepth[j] = np.mean(ccds.galdepth[Igood])
        else:
            E.galdepth[j] = 0.
    del expnums


    keep = np.array([not(expnum in bad_expids)
                     for expnum in ccds.expnum])
    ccds.cut(keep)
    print(len(ccds), 'CCDs NOT in the bad_expids file')

    keep = np.array([not(expnum in bad_expids)
                     for expnum in copilot.expnum])
    copilot.cut(keep)
    print(len(copilot), 'copilot exposures NOT in the bad_expids file')
    
    keep = np.array([not(expnum in bad_expids)
                     for expnum in E.expnum])
    E.cut(keep)
    print(len(E), 'CCD Exposures NOT in the bad_expids file')

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
    print()
    print()
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
        E.detiv[E.galdepth == 0] = 0.
        print('Smallest detivs:', E.detiv[np.argsort(E.detiv)[:10]])
        print('w/ galdepths:', E.galdepth[np.argsort(E.detiv)[:10]])
        print('Smallest positive detivs:', E.detiv[np.argsort(E.detiv + 1e12*(E.detiv == 0))[:10]])
        print('w/ galdepths:', E.galdepth[np.argsort(E.detiv + 1e12*(E.detiv == 0))[:10]])
        
    for band in bands:
        print()
        print('------------------')
        print(band, 'band.')

        # "I" indexes into exposures E.
        I = np.flatnonzero((E.filter == band) * E.photometric *
                           np.isfinite(E.detiv))
        print(len(I), 'photometric exposures in', band)

        # "iv" is parallel to O; will be converted to "galdepth".
        iv = np.zeros(len(O), np.float32)
        # "J" indexes into obstatus tiles O.
        J = tileid_to_index[E.tileid[I]]
        assert(np.all((J >= 0) * (J < len(O))))
        assert(np.all(O.tileid[J] == E.tileid[I]))

        #print('tileid range', E.tileid[I].min(), E.tileid[I].max())
        # d = np.array([degrees_between(*a) for a in
        #               zip(E.ra_bore[I], E.dec_bore[I], O.ra[J], O.dec[J])])
        # print('Degrees between tiles & exposures:', d)
        
        np.add.at(iv, J, E.detiv[I])
        print('galdepth range:', E.galdepth[I].min(), E.galdepth[I].max())
        print('detiv range:', E.detiv[I].min(), E.detiv[I].max())
        #print('index range:', J.min(), J.max())
        nexp = np.zeros(len(O), int)
        np.add.at(nexp, J, 1)
        print('tile exposure counts:', Counter(nexp))

        # convert iv back to galdepth in mags
        with np.errstate(divide='ignore'):
            galdepth = -2.5 * (np.log10(np.sqrt(1. / iv)) - 9)
            galdepth[iv == 0] = 0.

        # Shallowest before extinction correction
        #I = np.argsort(iv + 1e6*(iv == 0))
        I = np.argsort(galdepth + 50.*(galdepth == 0))
        print('Shallowest depth estimates from annotated CCDs file, before extinction:')
        for i in I[:10]:
            print('  ', galdepth[i], 'iv', iv[i], 'tile', O.tileid[i], 'expnum', O.get('%s_expnum' % band)[i])
            e = O.get('%s_expnum' % band)[i]
            j = np.flatnonzero(E.expnum == e)
            print('    galdepth', E.galdepth[j])
            
        fid = nom.fiducial_exptime(band)
        extinction = O.ebv_med * fid.A_co
        #print('Extinction range:', extinction.min(), extinction.max())
        galdepth -= extinction

        galdepth[iv == 0] = 0.

        # Shallowest galdepth > 0
        I = np.argsort(galdepth + 50.*(galdepth == 0))
        print('Shallowest depth estimates from annotated CCDs file:')
        for i in I[:10]:
            print('  ', galdepth[i], 'tile', O.tileid[i], 'expnum', O.get('%s_expnum' % band)[i])
        
        #print('galdepth deciles:', np.percentile(galdepth, [0,10,20,30,40,50,60,70,80,90,100]))

        # Z_DONE, Z_EXPNUM but no Z_DEPTH
        missing_depth = np.flatnonzero((O.get('%s_expnum' % band) > 0) * (O.get('%s_done' % band) == 1) * (galdepth == 0))
        print('Found', len(missing_depth), 'tiles with', band, 'DONE and EXPNUM but no DEPTH; setting to DEPTH=30')
        print('  eg, EXPNUMs', O.get('%s_expnum' % band)[missing_depth[:10]], 'DATE', O.get('%s_date' % band)[missing_depth[:10]])
        # Don't actually update 'galdepth[missing_depth]' until after this next check...

        # Flag tiles that have *only* non-photometric exposures with depth = 1.
        I = np.flatnonzero((E.filter == band) * np.logical_not(E.photometric))
        print(len(I), 'exposures are non-photometric in', band, 'band')
        J = tileid_to_index[E.tileid[I]]
        only_nonphot = J[galdepth[J] == 0.]
        print(len(only_nonphot), 'tiles have only non-photometric exposures')
        print('Marking', len(only_nonphot), 'non-photometric tiles in', band, 'with depth=1')

        orig_galdepth = galdepth.copy()
        
        galdepth[missing_depth] = 30.
        galdepth[only_nonphot] = 1.

        J = tileid_to_index[E.tileid[I]]
        nonphot = (galdepth[J] == 1.)
        print('Non-photometric galdepths:', E.galdepth[I])
        print('Non-photometric galdepths:', E.galdepth[I[nonphot]])
        plt.clf()
        phot = np.flatnonzero((E.filter == band) * E.photometric)
        plt.hist(E.galdepth[phot], range=(18, 26), bins=50, histtype='step', color='b', label='Photometric')
        plt.hist(E.galdepth[I[nonphot]], range=(18, 26), bins=50, histtype='step', color='r', label='Non-phot')
        plt.legend()
        plt.savefig('nonphot-%s.png' % band)
        
        # expnum_to_copilot = np.empty(expnums.max()+1, int)
        # expnum_to_copilot[:] = -1
        # expnum_to_copilot[copilot.expnum] = np.arange(len(copilot))
        expnum_to_copilot = dict([(e,i) for i,e in enumerate(copilot.expnum)])

        if False:
            # Let's check the accuracy of the copilot's depth estimates...
            target_exptime = copilot.expfactor * fid.exptime
            # What fraction of the target exposure time did we take?
            depth_factor = copilot.exptime / target_exptime
            nomdepth = fid.single_exposure_depth
            depth = nomdepth + 2.5 * np.log10(np.sqrt(depth_factor))
            #print('Copilot predicted depths:', depth)
            IC = np.array([expnum_to_copilot.get(e, -1) for e in allccds.expnum])
            K = np.flatnonzero(IC >= 0)
            ext = np.array([e['ugrizY'.index(f)] for e,f in zip(allccds.decam_extinction, allccds.filter)])
            dd = allccds.galdepth - ext
            print('Making scatterplot...', len(K), 'points')
            plt.clf()
            #plt.plot(dd[K], depth[IC[K]], 'b.', alpha=0.2, mec='none')
            plt.scatter(dd[K], depth[IC[K]],
                        c=np.clip(copilot.expfactor[IC[K]], 0, 2), s=10, alpha=0.2, edgecolors='none')
            plt.colorbar()
            plt.xlabel('Pipeline depth')
            plt.ylabel('Copilot depth proxy')
            plt.plot([20,25],[20,25], 'k-', alpha=0.25)
            plt.plot([20,25],[20+0.1,25+0.1], 'k--', alpha=0.25)
            plt.plot([20,25],[20-0.1,25-0.1], 'k--', alpha=0.25)
            plt.axis([20.5, 23, 21, 23.5])
            plt.title('Copilot vs Pipeline depth estimates.  (color = exp.factor)')
            plt.savefig('depth-copilot-%s.png' % band)
            print('Made scatterplot')

        plt.clf()
        ha = dict(bins=60, range=(0,30), log=True, histtype='step')
        plt.hist(O.get('%s_depth' % band), color='k', label='Before', **ha)
        plt.hist(orig_galdepth, color='b', label='Annotated CCDs', **ha)

        # Do we have measurements for any of these missing tiles in the copilot db?
        for code in [30, 0]:
            Igal = np.flatnonzero((O.get('%s_expnum' % band) > 0) *
                                  (O.get('%s_done' % band) == 1) *
                                  (galdepth == code))
            expnums = O.get('%s_expnum' % band)[Igal]
            print(len(expnums), 'still marked DONE, with EXPNUM, but missing DEPTH, with code =', code)

            Ihuh = np.flatnonzero((O.get('%s_done' % band) == 1) *
                                  (galdepth == code))
            print(len(Ihuh), 'tiles marked DONE, without EXPNUM, and DEPTH =', code)

            if len(Ihuh):
                print('Tile ids:', O.tileid[Ihuh])
                for t in O.tileid[Ihuh]:
                    I = np.flatnonzero(E.tileid == t)
                    print('  tile', t, ': exposure numbers:', E.expnum[I])
                    print('  with depths', E.galdepth[I])
                    i = tileid_to_index[t]
                    if i >= 0:
                        print('    depth', galdepth[i])
                    else:
                        print('    no depth')

                # Within an arcmin?
                I,J,d = match_radec(O.ra[Ihuh], O.dec[Ihuh], all_copilot.rabore, all_copilot.decbore, 1./60., nearest=True)
                print('For', len(Ihuh), 'weird tiles,')
                print(len(I), 'matches within an arcmin in the copilot db')
                print('Smallest distances:', d[np.argsort(d)[:10]])

                I,J,d = match_radec(O.ra[Ihuh], O.dec[Ihuh], allccds.ra_bore, allccds.dec_bore, 1./60., nearest=True)
                print(len(I), 'matches within an arcmin in the CCDs table')
                print('Smallest distances:', d[np.argsort(d)[:10]])

                O[Ihuh].writeto('weird-%s-%i.fits' % (band,code))

            
            if len(expnums) == 0:
                continue
            IC = np.array([expnum_to_copilot.get(e, -1) for e in expnums])
            K = np.flatnonzero(IC >= 0)
            expnums = expnums[K]
            # these are the indices into O / galdepth
            Igal = Igal[K]
            co = copilot[IC[K]]
            print(len(expnums), 'matched to copilot database')

            target_exptime = co.expfactor * fid.exptime
            # What fraction of the target exposure time did we take?
            depth_factor = co.exptime / target_exptime
            nomdepth = fid.single_exposure_depth
            print('Nominal single-exposure depth:', nomdepth)
            co.depth = nomdepth + 2.5 * np.log10(np.sqrt(depth_factor))
            print('Copilot predicted depths:', co.depth)
            J = np.flatnonzero(np.isfinite(co.depth))
            co = co[J]
            # indices into O
            Igal = Igal[J]
            print(len(Igal), 'good copilot depth estimates')
        
            pcts = [0,1,5,25,50,75,95,99,100]
            print('Copilot depth percentiles:', np.percentile(co.depth, pcts))

            print('Shallowest exposures:')
            I = np.argsort(co.depth)
            for i in I[:10]:
                print('  Expnum', co.expnum[i], 'depth', co.depth[i], 'exptime', co.exptime[i])
            co[I].writeto('depths.fits')

            from astrometry.util.starutil_numpy import mjdtodate
            print('Copilot-matched entries:')
            I = np.argsort(co.expnum)
            for i in I:
                print('  EXPNUM', co.expnum[i], 'date', mjdtodate(co.mjd_obs[i]), '  copilot name', co.filename[i])
                e = co.expnum[i]
                I = np.flatnonzero(allccds.expnum == e-1)
                fn1 = None
                fn2 = None
                if len(I):
                    print('    CCDs file contains', len(I), 'entries for expnum', e-1)
                    print('      filename', allccds.image_filename[I[0]])
                    fn1 = allccds.image_filename[I[0]]
                else:
                    print('    No CCDs file entries for expnum', e-1)
                I = np.flatnonzero(allccds.expnum == e+1)
                if len(I):
                    print('    CCDs file contains', len(I), 'entries for expnum', e+1)
                    print('      filename', allccds.image_filename[I[0]])
                    fn2 = allccds.image_filename[I[0]]
                else:
                    print('    No CCDs file entries for expnum', e+1)

                if fn1 is not None and fn2 is not None:
                    full1 = os.path.join(survey.get_image_dir(), fn1)
                    #print('Full path 1:', full1)
                    if os.path.exists(full1):
                        print('exists')
                    full2 = os.path.join(survey.get_image_dir(), fn2)
                    #print('Full path 2:', full2)
                    if os.path.exists(full2):
                        print('exists')
                    if os.path.exists(full1) and os.path.exists(full2):
                        dir1 = os.path.dirname(full1)
                        dir2 = os.path.dirname(full2)
                        if dir1 == dir2:
                            #print('dir:', dir1)
                            fns = os.listdir(dir1)
                            fns.sort()
                            fns = [fn for fn in fns if ('oki' in fn or 'ooi' in fn)]
                            base1 = os.path.basename(full1)
                            base2 = os.path.basename(full2)
                            i1 = fns.index(base1)
                            i2 = fns.index(base2)
                            print('Files found at list elements', i1, i2)
                            #print(fns[i1:i2+1])
                            for fn in fns[i1:i2+1]:
                                print('EXPNUM', e, 'range', os.path.join(os.path.dirname(fn1), fn))
                            if i1 + 4 == i2:
                                print('EXPNUM', e, 'expected', os.path.join(os.path.dirname(fn1), fns[i1+2]))
                            if i1 + 2 == i2:
                                print('EXPNUM', e, 'expected', os.path.join(os.path.dirname(fn1), fns[i1+1]))
                            

            #print('Before:', galdepth[Igal])
            galdepth[Igal] = co.depth
            #print('After:', galdepth[Igal])

            Igal = np.flatnonzero((O.get('%s_expnum' % band) > 0) *
                                  (O.get('%s_done' % band) == 1) *
                                  (galdepth == code))
            expnums = O.get('%s_expnum' % band)[Igal]
            print(len(expnums), 'still marked DONE, with EXPNUM, but missing DEPTH with code =', code, 'after copilot patching')
            print('Exposure numbers:', expnums)
            print('Exposure dates:', O.get('%s_date' % band)[Igal])


        
        O.set('%s_depth' % band, galdepth)

        plt.hist(O.get('%s_depth' % band), color='r', label='After', **ha)
        plt.savefig('depth-hist-%s.png' % band)

        #print('Depth deciles: [', ', '.join(['%.3f' % f for f in np.percentile(O.get('%s_depth' % band), [0,10,20,30,40,50,60,70,80,90,100])]) + ']')

        rlo,rhi = 0,360
        dlo,dhi = declo,dechi
        J = np.flatnonzero((O.in_desi == 1) * (O.in_des == 0) *
                           (O.dec > dlo) * (O.dec < dhi))
        print('Median E(B-V) in DECaLS area:', np.median(O.ebv_med[J]))
        print('Median extinction in DECaLS area, %s band:' % band, np.median(extinction[J]))

        I2 = np.flatnonzero((O.get('%s_expnum' % band) > 0) *
                            (O.get('%s_depth' % band) == 30) *
                            (O.get('%s_done' % band) == 1) *
                            (O.in_desi == 1))
        print(len(I2), 'with EXPNUM and DONE and IN_DESI, but no DEPTH')
        # Sort by expnum
        I2 = I2[np.argsort(O.get('%s_expnum' % band)[I2])]
        print('Exposure numbers:', sorted(O.get('%s_expnum' % band)[I2]))
        print('Dates:', sorted(O.get('%s_date' % band)[I2]))
        print('Dates:', np.unique(O.get('%s_date' % band)[I2]))

        for i in I2:
            print('  date', O.get('%s_date' % band)[i], 'expnum', O.get('%s_expnum' % band)[i])
        
        # for i2,o in zip(I2, O[I2]):
        #     print()
        #     e = o.get('%s_expnum' % band)
        #     print('  Expnum', e, 'orig galdepth', orig_galdepth[i2])
        #     date = o.get('%s_date' % band)
        #     print('  Date', date)
        #     print('  Pass', o.get('pass'), 'Tile', o.tileid)
        #     jj = np.flatnonzero(allccds.expnum == e)
        #     print('  In DESI:', o.in_desi, 'In DES:', o.in_des)
        #     print('  ', len(jj), 'matching CCDs')
        #     if len(jj) == 0:
        #         continue
        #     print('  CCDs OBJECT', [ob.strip() for ob in allccds.object[jj]])
        #     print('  CCDs Tileid', allccds.tileid[jj])
        #     print('  CCDs galdepth', allccds.galdepth[jj])
        #     print('  CCDs photometric', allccds.photometric[jj])
        # 
        #     ii = np.flatnonzero(E.expnum == e)
        #     print('  ', len(ii), 'Exposures matching')
        #     if len(ii):
        #         ee = E[ii[0]]
        #         print('  exposure tileid', ee.tileid)
        #         print('  index', tileid_to_index[ee.tileid])
        #         print('  vs i2=', i2)
        #         print('  only_nonphot', only_nonphot[i2], 'missing_depth', missing_depth[i2])
        #     
        #     kk = np.flatnonzero(allccds.tileid == o.tileid)
        #     kk = np.array(sorted(set(kk) - set(jj)))
        #     print('  ', len(kk), 'other CCDs of this tile')
        # #print('Dates:', O.get('%s_date' % band)[I])

        from astrometry.util.plotutils import antigray
        rr,dd = np.meshgrid(np.linspace(rlo,rhi,720), np.linspace(dlo,dhi,360))
        JJ,II,d = match_radec(rr.ravel(), dd.ravel(), O.ra, O.dec, 1.5,
                              nearest=True)
        indesi = np.zeros(rr.shape, bool)
        indesi.flat[JJ] = ((O.in_desi[II] == 1) * (O.in_des[II] == 0))
        
        plt.figure(figsize=(14,6))
        plt.subplots_adjust(left=0.1, right=0.99)
        
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
            pct = np.percentile(depth[J], [0,10,20,30,40,50,60,70,80,90,100])
            #print('Depth deciles:', np.percentile(depth[J], [0,10,20,30,40,50,60,70,80,90,100]))
            print('Depth deciles: [', ', '.join(['%.3f' % f for f in pct]) + ']')

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

        print('Fiducial single-exposure-depth:', fid.single_exposure_depth)

        for passnum in [1,2,3]:
            depth = O.get('%s_depth' % band)
            J = np.flatnonzero((O.get('pass') == passnum) *
                               (depth > 1) * (depth < 30))
            depth = depth[J]

            print('Pass', passnum)
            print(sum(depth < fid.single_exposure_depth - 0.25), 'of', len(depth), 'tiles are more than 0.25 mag shallow')

            odepth = O.get('%s_depth' % band)
            K = np.flatnonzero((O.get('%s_done' % band) == 0) * (O.get('pass') == passnum) *
                               (odepth > 1) * (odepth < 30))
            print(sum(odepth[K] < fid.single_exposure_depth - 0.25), 'of', len(odepth[K]), 'DONE=0 tiles are more than 0.25 mag shallow')

            for k in K:
                print('  EXPNUM', O.get('%s_expnum' % band)[k], 'DATE', O.get('%s_date' % band)[k], 'DEPTH', O.get('%s_depth' % band)[k])

            K = np.flatnonzero((O.get('%s_done' % band) == 1) * (O.get('pass') == passnum) *
                               (odepth > 1) * (odepth < 30))
            print(sum(odepth[K] < fid.single_exposure_depth - 0.25), 'of', len(odepth[K]), 'DONE=1 tiles are more than 0.25 mag shallow')

            K = np.flatnonzero((O.get('%s_done' % band) == 1) * (O.get('pass') == passnum) *
                               (odepth == 1))
            print(len(K), 'DONE=1 tiles have DEPTH=1 (non-photometric)')

            K = np.flatnonzero((O.get('%s_done' % band) == 1) * (O.get('pass') == passnum) *
                               (odepth == 30))
            print(len(K), 'DONE=1 tiles have DEPTH=30 (unknown depth)')

            K = np.flatnonzero((O.get('%s_done' % band) == 1) * (O.get('pass') == passnum) *
                               (odepth == 0))
            print(len(K), 'DONE=1 tiles have DEPTH=0')
            
            K = np.flatnonzero((O.get('%s_done' % band) == 0) * (O.get('pass') == passnum) *
                               (odepth != 0))
            print(len(K), 'tiles have DONE=0 but DEPTH != 0')

            
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

        for passnum in [1,2,3]:
            depth = O.get('%s_depth' % band)

            roi = ((O.in_desi == 1) * (O.in_des == 0) *
                   (O.dec > dlo) * (O.dec < dhi) *
                   (O.get('pass') == passnum))

            J = np.flatnonzero(roi)
            
            done = np.flatnonzero(roi *
                                  (O.get('%s_done' % band) == 1))
            
            redo = np.flatnonzero(roi *
                                  np.logical_or(
                                      (depth > 1) * (depth < 30) *
                                      (depth < fid.single_exposure_depth - 0.25),
                                      depth == 1))

            print('Band %s, pass %i: total tiles %i, done %i, redo %i, keep %i' %
                  (band, passnum, len(J), len(done), len(redo), len(done)-len(redo)))

            A = np.flatnonzero(roi * 
                               (depth > 1) * (depth < 30) *
                               (depth > fid.single_exposure_depth - 0.25))
            B = np.flatnonzero(roi * 
                               (depth > 1) * (depth < 30) *
                               (depth <= fid.single_exposure_depth - 0.25))

            C = np.flatnonzero(roi * (depth == 1))

            D = np.flatnonzero(roi * (depth == 30))
            
            print('Band %s, pass %i: total tiles: %i, A: %i, B: %i, C: %i, D: %i' %
                  (band, passnum, len(J), len(A), len(B), len(C), len(D)))
            
            plt.clf()
            plt.plot(O.ra[J], O.dec[J], 'ko', alpha=0.1)
            plt.plot(O.ra[done], O.dec[done], 'k.')
            plt.plot(O.ra[redo], O.dec[redo], 'r.')
            plt.axis([360, 0, -20, 38])
            plt.title('Tiles to redo: %s band, pass %i: %i of %i' % (band, passnum, len(redo), len(done)))
            plt.savefig('redo-%s-%i.png' % (band, passnum))

            if band == 'z':

                redo = np.flatnonzero((O.get('pass') == passnum) *
                                      np.logical_or(
                                          (depth > 1) * (depth < 30) *
                                          (depth < fid.single_exposure_depth - 0.5),
                                          depth == 1))

                A = np.flatnonzero(roi * 
                                   (depth > 1) * (depth < 30) *
                                   (depth > fid.single_exposure_depth - 0.5))
                B = np.flatnonzero(roi * 
                                   (depth > 1) * (depth < 30) *
                                   (depth <= fid.single_exposure_depth - 0.5))

                print('Band %s, pass %i: total tiles: %i, A: %i, B: %i, C: %i, D: %i (shallow = 0.5 mag less than target)' %
                      (band, passnum, len(J), len(A), len(B), len(C), len(D)))
                
                plt.clf()
                plt.plot(O.ra[J], O.dec[J], 'ko', alpha=0.1)
                plt.plot(O.ra[done], O.dec[done], 'k.')
                plt.plot(O.ra[redo], O.dec[redo], 'r.')
                plt.axis([360, 0, -20, 38])
                plt.title('Tiles to redo (> 0.5 mag shallow): %s band, pass %i: %i of %i' % (band, passnum, len(redo), len(done)))
                plt.savefig('redo2-%s-%i.png' % (band, passnum))


        print('Passes 1-3 combined:')
        depth = O.get('%s_depth' % band)
        J = np.flatnonzero((depth > 1) * (depth < 30))
        depth = depth[J]
        print(sum(depth < fid.single_exposure_depth - 0.25), 'of', len(depth), 'tiles are more than 0.25 mag shallow')

        odepth = O.get('%s_depth' % band)
        K = np.flatnonzero((O.get('%s_done' % band) == 0) * (odepth > 1) * (odepth < 30))
        print(sum(odepth[K] < fid.single_exposure_depth - 0.25), 'of', len(odepth[K]), 'DONE=0 tiles are more than 0.25 mag shallow')

        for k in K:
            print('  EXPNUM', O.get('%s_expnum' % band)[k], 'DATE', O.get('%s_date' % band)[k], 'DEPTH', O.get('%s_depth' % band)[k])

        K = np.flatnonzero((O.get('%s_done' % band) == 1) * (odepth > 1) * (odepth < 30))
        print(sum(odepth[K] < fid.single_exposure_depth - 0.25), 'of', len(odepth[K]), 'DONE=1 tiles are more than 0.25 mag shallow')

        K = np.flatnonzero((O.get('%s_done' % band) == 1) * (odepth > 1) * (odepth < 30) * goodtiles)
        print(sum(odepth[K] < fid.single_exposure_depth - 0.25), 'of', len(odepth[K]), 'interesting DONE=1 tiles are more than 0.25 mag shallow')

        K = np.flatnonzero((O.get('%s_done' % band) == 1) * (odepth == 1))
        print(len(K), 'DONE=1 tiles have DEPTH=1 (non-photometric)')

        K = np.flatnonzero((O.get('%s_done' % band) == 1) * (odepth == 1) * goodtiles)
        print(len(K), 'interesting DONE=1 tiles have DEPTH=1 (non-photometric)')

        K = np.flatnonzero((O.get('%s_done' % band) == 1) * (odepth == 30))
        print(len(K), 'DONE=1 tiles have DEPTH=30 (unknown depth)')

        K = np.flatnonzero((O.get('%s_done' % band) == 1) * (odepth == 30) * goodtiles)
        print(len(K), 'interesting DONE=1 tiles have DEPTH=30 (unknown depth)')

        K = np.flatnonzero((O.get('%s_done' % band) == 1) * (odepth == 0))
        print(len(K), 'DONE=1 tiles have DEPTH=0')

        K = np.flatnonzero((O.get('%s_done' % band) == 1) * (odepth == 0) * goodtiles)
        print(len(K), 'interesting DONE=1 tiles have DEPTH=0')
        
        K = np.flatnonzero((O.get('%s_done' % band) == 0) *
                           (odepth != 0))
        print(len(K), 'tiles have DONE=0 but DEPTH != 0')

                
            
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
        cdict[key] = [ (indices[i], colors_rgba[i-1,ki], colors_rgba[i,ki]) for i in range(N+1) ]
    # Return colormap object.
    return matplotlib.colors.LinearSegmentedColormap(cmap.name + "_%d"%N, cdict, 1024)


if __name__ == '__main__':
    main()

