import unittest
from django.test import TestCase
import os

class TestCopilot(TestCase):

    def setUp(self):
        # os.setenv('PS1CAT_DIR',
        self.testdatadir = os.path.join(os.path.dirname(__file__),
                                        'testdata')

    def test_cmdline_args(self):
        from copilot import main, skip_existing_files
        
        from obsdb.models import MeasuredCCD
        fn = os.path.join(self.testdatadir, 'mos3.68488.im4.fits.fz')

        args = ['--n-fwhm', '10']
        
        # Check that 'skip_existing_files' works correctly
        rawext = 'im4'
        good = skip_existing_files([fn], rawext)
        self.assertEqual(len(good), 1)

        # Database starts empty
        self.assertEqual(MeasuredCCD.objects.all().count(), 0)
        
        main(cmdlineargs = [fn, '--no-db'] + args)

        # With --no-db, database stays empty
        self.assertEqual(MeasuredCCD.objects.all().count(), 0)

        # Database empty, so not skipped
        good = skip_existing_files([fn], rawext)
        self.assertEqual(len(good), 1)

        main(cmdlineargs = [fn] + args)

        # Now database has one entry.
        self.assertEqual(MeasuredCCD.objects.all().count(), 1)

        m = MeasuredCCD.objects.first()
        self.assertTrue(abs(m.zeropoint    - 26.35) < 0.02)
        self.assertTrue(abs(m.sky          - 19.39) < 0.02)
        self.assertTrue(abs(m.transparency - 0.953) < 0.02)
        self.assertTrue(abs(m.seeing       - 0.991) < 0.02)
        self.assertTrue(abs(m.expfactor    - 0.422) < 0.01)
        
        # Now file should be skipped
        good = skip_existing_files([fn], rawext)
        self.assertEqual(len(good), 0)
        
        main(cmdlineargs = [fn, '--skip'] + args)

        self.assertEqual(MeasuredCCD.objects.all().count(), 1)

        
    def test_new_image(self):
        pass


if __name__ == '__main__':
    unittest.main()
    
