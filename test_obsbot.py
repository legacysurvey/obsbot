import unittest
from django.test import TestCase
import os

class TestCopilot(TestCase):

    def setUp(self):
        # os.setenv('PS1CAT_DIR',
        self.testdatadir = os.path.join(os.path.dirname(__file__),
                                        'testdata')

    def test_cmdline_args(self):
        from copilot import main
        from obsdb.models import MeasuredCCD
        fn = os.path.join(self.testdatadir, 'mos3.68488.im4.fits.fz')
        
        self.assertEqual(MeasuredCCD.objects.all().count(), 0)
        
        main(cmdlineargs = [fn, '--no-db'])

        self.assertEqual(MeasuredCCD.objects.all().count(), 0)

        main(cmdlineargs = [fn])

        self.assertEqual(MeasuredCCD.objects.all().count(), 1)

        main(cmdlineargs = [fn, '--skip'])

        self.assertEqual(MeasuredCCD.objects.all().count(), 1)

        
    def test_new_image(self):
        pass


if __name__ == '__main__':
    unittest.main()
    
