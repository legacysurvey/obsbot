import unittest
#from django.test import TestCase
import os

class TestMosbot(unittest.TestCase):

    def setUp(self):
        self.testdatadir = os.path.join(os.path.dirname(__file__),
                                        'testdata')

    def test_new_file(self):
        from mosbot import main
        import tempfile
        tempdir = tempfile.mkdtemp()
        args = ['--script', os.path.join(tempdir, 'tonight.sh')]
        ## FIXME -- will need to update these to have current dates; also put in git
        args += ['pass1.json', 'pass2.json', 'pass3.json']

        mosbot = main(cmdlineargs=args, get_mosbot=True)

        # write sequence number
        f = open(mosbot.seqnumpath, 'w')
        f.write('10\n')
        f.close()
        
        fn = os.path.join(self.testdatadir, 'mos3.68488.im4.fits.fz')
        mosbot.found_new_image(fn)
            

if __name__ == '__main__':
    unittest.main()
