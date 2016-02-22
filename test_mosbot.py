from __future__ import print_function
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
        fn1 = os.path.join(self.testdatadir, 'pass1.json')
        fn2 = os.path.join(self.testdatadir, 'pass2.json')
        fn3 = os.path.join(self.testdatadir, 'pass3.json')

        tmpfn1 = fn1 + '.tmp'
        tmpfn2 = fn2 + '.tmp'
        tmpfn3 = fn3 + '.tmp'
        
        for fn,tmpfn in ((fn1,tmpfn1),(fn2,tmpfn2),(fn3,tmpfn3)):
            import json
            import ephem
            
            J = json.loads(open(fn, 'r').read())
            t0 = ephem.Date(str(J[0]['approx_datetime']))
            print('First exposure:', t0)
            now = ephem.now()
            print('Now:', now)
            for j in J:
                tnew = now + ephem.Date(str(j['approx_datetime'])) - t0
                tnew = str(ephem.Date(tnew))
                j['approx_datetime'] = tnew
                print('Updated datetime to', tnew)

            f = open(tmpfn, 'w')
            json.dump(J, f, sort_keys=True, indent=4, separators=(',', ': '))
            f.close()
                
        args += [tmpfn1, tmpfn2, tmpfn3]

        mosbot = main(cmdlineargs=args, get_mosbot=True)

        # write sequence number
        f = open(mosbot.seqnumpath, 'w')
        f.write('10\n')
        f.close()
        
        fn = os.path.join(self.testdatadir, 'mos3.68488.im4.fits.fz')
        mosbot.found_new_image(fn)
            

if __name__ == '__main__':
    unittest.main()
