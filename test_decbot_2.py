from __future__ import print_function

import os
from django.test import TestCase
import obsdb
import numpy as np

from astrometry.util.fits import fits_table

class Duck(object):
    quack = True

import ephem
real_ephem_now = ephem.now
ephem_offset = 0.
def fake_ephem_now():
    now = real_ephem_now()
    # print('Real ephem.now:', now, '-- adding offset', ephem_offset)
    return now + ephem_offset
ephem.now = fake_ephem_now


def set_fake_mjd(target_mjd):
    import obsbot
    # How many ways do we get the time?
    global ephem_offset
    ephem_offset = 0.
    nownow = ephem.now()
    ephem_offset = obsbot.mjd_to_ephem_date(target_mjd) - nownow

    obsbot.mjdnow_offset = 0
    nownow = obsbot.mjdnow()
    obsbot.mjdnow_offset = target_mjd - nownow


class TestDecbot2(TestCase):

    def setUp(self):
        # fitscopy obs.fits"[mjd_obs > 57603.0 && mjd_obs < 57603.1]" \
        #      obs-test.fits
        self.data_dir = os.path.join(os.path.dirname(__file__), 'testdata')

        T = fits_table(os.path.join(self.data_dir, 'obs-test.fits'))
        for t in T:
            kwargs = {}
            for c in ['camera', 'filename', 'extension', 'obstype',
                      'object', 'band']:
                kwargs[c] = t.get(c).strip()
            for c in ['expnum', 'exptime', 'mjd_obs', 'airmass',
                      'racenter', 'deccenter', 'rabore', 'decbore',
                      'tileid', 'passnumber', 'tileebv', 'ebv', 'zeropoint',
                      'transparency', 'seeing', 'sky', 'expfactor', 'dx', 'dy',
                      'nmatched', 'md5sum', 'bad_pixcnt', 'readtime',
                      'affine_dx', 'affine_dxx', 'affine_dxy', 'affine_x0',
                      'affine_dy', 'affine_dyx', 'affine_dyy', 'affine_y0',
                      ]:
                kwargs[c] = t.get(c)
            obsdb.MeasuredCCD.objects.create(**kwargs)
            
    def test_recent(self):
        from obsbot import mjdnow
        from copilot import (get_recent_ccds, get_recent_exposures,
                             recent_gr_seeing)

        import obsbot

        # Fake the current time...
        target_mjd = 57603.1
        set_fake_mjd(target_mjd)

        now = mjdnow()
        self.assertLess(np.abs(now - target_mjd), 0.001)

        ccds = get_recent_ccds(recent = 0.1 * 24 * 60)
        self.assertEqual(len(ccds), 73)

        ccds = get_recent_ccds(recent = 30.)
        self.assertEqual(len(ccds), 15)

        exps = get_recent_exposures(recent = 30.)
        self.assertEqual(len(exps), 15)
        
        exps = get_recent_exposures(recent = 30., bands=['g','r'])
        self.assertEqual(len(exps), 15)

        gexps = get_recent_exposures(recent = 30., bands=['g'])
        self.assertEqual(len(gexps), 6)

        rexps = get_recent_exposures(recent = 30., bands=['r'])
        self.assertEqual(len(rexps), 9)

        xexps = get_recent_exposures(recent = 30., bands=[])
        self.assertEqual(xexps, None)
        
        gsee, rsee, G, R = recent_gr_seeing()
        print('gsee', gsee)
        print('rsee', rsee)
        print('G', G)
        print('R', R)

        self.assertLess(np.abs(gsee - 1.360), 0.001)
        self.assertLess(np.abs(rsee - 1.374), 0.001)
        self.assertEqual(len(G), 5)
        self.assertEqual(len(R), 5)
        
    def test_decbot_recent(self):
        from decbot import Decbot
        import json
        import obsdb
        import camera_decam as camera
        import obsbot
        from obsbot import mjdnow
        import ephem

        # Fake the current time...
        set_fake_mjd(57603.1)

        J1 = json.loads(open(os.path.join(self.data_dir, '2016-08-02-p1.json'))
                        .read())
        J2 = json.loads(open(os.path.join(self.data_dir, '2016-08-02-p2.json'))
                        .read())
        J3 = json.loads(open(os.path.join(self.data_dir, '2016-08-02-p3.json'))
                        .read())

        print('Read', len(J1), len(J2), len(J3), 'JSON plans')
        
        # Annotate with 'planpass' field
        for i,J in enumerate([J1,J2,J3]):
            for j in J:
                j['planpass'] = i+1

        opt = Duck()
        opt.cut_before_now = True
        opt.rawdata = 'no-such-directory'
        opt.ignore_mising_dir = True
        opt.verbose = False
        opt.adjust = False
        opt.passnum = 2
        opt.exptime = 100
        opt.start_double = False
        
        nom = camera.nominal_cal
        obs = camera.ephem_observer()
        # ???
        tiles = fits_table(camera.tile_path)
        rc = None
        copilot_db = obsdb.MeasuredCCD.objects
        decbot = Decbot(J1, J2, J3, opt, nom, obs, tiles, rc,
                        copilot_db=copilot_db)

        decbot.observed_tiles.clear()

        T = fits_table(os.path.join(self.data_dir, 'obs-test.fits'))
        t = T[np.argmax(T.mjd_obs)]

        M = dict(band=t.band.strip(),
                 transparency=t.transparency,
                 seeing=t.seeing,
                 skybright=t.sky,
                 airmass=t.airmass,
                 )
        decbot.latest_measurement = M
        decbot.recent_gr(M)

        print('M', M)

        self.assertEqual(M['band'], 'r')
        gsee,rsee,G,R = M['grsee']
        self.assertLess(np.abs(gsee - 1.360), 0.001)
        self.assertLess(np.abs(rsee - 1.374), 0.001)
        self.assertLess(np.abs(M['grsky'] - 0.826), 0.001)

        decbot.update_plans()

        self.assertEqual(len(decbot.J1), 416)
        self.assertEqual(len(decbot.J2), 416)
        self.assertEqual(len(decbot.J3), 416)
        self.assertEqual(decbot.nextpass, 2)

        jnext = decbot.get_upcoming()[0]
        print('Next exposure:', jnext)
        self.assertEqual(jnext['expTime'], 76.)
        self.assertEqual(jnext['object'], 'DECaLS_18778_r')

        jnext = decbot.get_upcoming()[1]
        print('Next+1 exposure:', jnext)
        self.assertEqual(jnext['expTime'], 192.)
        self.assertEqual(jnext['object'], 'DECaLS_18778_g')
        
