from __future__ import print_function

from astrometry.util.fits import *

class Duck(object):
    pass

duck = Duck()
duck.rawdata = 'rawdata'
duck.passnum = 1
duck.exptime = 100
duck.verbose = True
duck.adjust = True


# def test_adj(nom, obs, tiles, bot):
#     band = 'z'
#     fid = nom.fiducial_exptime(band)
#     for i in range(100):
#         print()
#         print('Tile', i)
#         bot.adjust_for_previous(tiles[i], band, fid, debug=True)


def plot_exposure(plot, ra, dec, wcses):
    for wcs in wcses:
        wcs.set_crval((ra, dec))
        plot.outline.wcs = anwcs_new_sip(wcs)
        plot.plot('outline')
    

if __name__ == '__main__':

    if False:
        from decbot import Decbot
        import camera_decam
        # import nominal_cal, ephem_observer
        
        nom = camera_decam.nominal_cal
        obs = camera_decam.ephem_observer()
        tiles = fits_table('obstatus/decam-tiles_obstatus.fits')
        bot = Decbot([],[],[], duck, nom, obs, tiles, None)
    
        band = 'z'
        fid = nom.fiducial_exptime(band)
        for i in range(100):
            print()
            print('Tile', i)
            bot.adjust_for_previous(tiles[i], band, fid, debug=True)


    duck.write_script = False
    duck.scriptfn = './test.sh'
        
    from mosbot import Mosbot
    import camera_mosaic
    
    nom = camera_mosaic.nominal_cal
    obs = camera_mosaic.ephem_observer()
    tiles = fits_table('obstatus/mosaic-tiles_obstatus.fits')
    tiles.cut(tiles.get('pass') <= 3)
    tiles.cut(tiles.dec >= 30)
    bot = Mosbot([],[],[], duck, nom, obs, tiles)

    from astrometry.blind.plotstuff import Plotstuff
    from astrometry.util.util import Sip, anwcs, anwcs_new_sip, wcs_pv2sip_hdr
    from astrometry.util.fits import fits_table
    from astrometry.libkd.spherematch import match_radec
    from astrometry.util.plotutils import *
    import fitsio
    import numpy as np

    fn = 'k4m_160203_061211_ooi_zd_v2.fits.fz'
    wcses = []
    for ext in range(1, 4+1):
        hdr = fitsio.read_header(fn, ext=ext)
        wcs = wcs_pv2sip_hdr(hdr)
        wcses.append(wcs)

    ps = PlotSequence('tile', format='%03i')

    # i0 = 1000
    # ii = range(i0, i0+10)

    np.random.seed(42)
    ii = np.random.randint(0, len(tiles), size=20)
    
    for i in ii:
        print()
        print('Tile', i)

        tile = tiles[i]
        ra,dec = tile.ra, tile.dec
        
        PW,PH = 800,800
        plot = Plotstuff(size=(PW, PH), rdw=(ra, dec, 2), outformat='png')

        others = bot.other_passes(tiles[i], tiles)
        others.rename('pass', 'passnum')

        for passnum in [1, 2, 3]:

            plot.color = 'verydarkblue'
            plot.plot('fill')

            plot.outline.fill = False

            K = np.flatnonzero(tiles.get('pass') == passnum)
            I,J,d = match_radec(np.array([tile.ra]), np.array([tile.dec]), tiles.ra[K], tiles.dec[K], 1.)
            Tnear = tiles[K[J]]
            plot.color = 'gray'
            for r,d in zip(Tnear.ra, Tnear.dec):
                plot_exposure(plot, r, d, wcses)

            plot.color = 'red'
            plot_exposure(plot, tile.ra, tile.dec, wcses)

            
            plot.color = 'white'
            plot.alpha = 0.25
            plot.outline.fill = True
            plot.apply_settings()

            I = np.flatnonzero(others.passnum == passnum)
            for ird,(r,d) in enumerate(zip(others.ra[I], others.dec[I])):

                for wcs in wcses:
                    wcs.set_crval((r, d))
                    plot.outline.wcs = anwcs_new_sip(wcs)
                    plot.plot('outline')

            plot.write(ps.getnext())

    sys.exit(0)
            
    band = 'z'
    fid = nom.fiducial_exptime(band)
    for i in range(100):
        print()
        print('Tile', i)
        bot.adjust_for_previous(tiles[i], band, fid, debug=True)
