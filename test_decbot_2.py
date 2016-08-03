from __future__ import print_function

import os
from django.test import TestCase
import obsdb
import numpy as np

from astrometry.util.fits import fits_table

class Duck(object):
    quack = True

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
        obsbot.mjdnow_offset = 0
        nownow = mjdnow()
        obsbot.mjdnow_offset = 57603.1 - nownow

        now = mjdnow()
        self.assertLess(np.abs(now - 57603.1), 0.001)

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
        
        # Fake the current time...
        obsbot.mjdnow_offset = 0
        nownow = mjdnow()
        obsbot.mjdnow_offset = 57603.1 - nownow
        
        J1 = json.loads(open(os.path.join(self.data_dir, '2016-08-02-p1.json'))
                        .read())
        J2 = json.loads(open(os.path.join(self.data_dir, '2016-08-02-p2.json'))
                        .read())
        J3 = json.loads(open(os.path.join(self.data_dir, '2016-08-02-p3.json'))
                        .read())
        opt = Duck()
        nom = camera.nominal_cal
        obs = camera.ephem_observer
        # ???
        tiles = fits_table(camera.tile_path)
        rc = None
        copilot_db = obsdb.MeasureCCD.objects
        decbot = Decbot(J1, J2, J3, opt, nom, obs, tiles, rc,
                        copilot_db=copilot_db)

        M = dict()
        #decbot.
        decbot.recent_gr(M)

        print('M', M)
