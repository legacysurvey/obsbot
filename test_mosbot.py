from __future__ import print_function
import unittest
#from django.test import TestCase
import os

class TestMosbot(unittest.TestCase):

    def setUp(self):
        self.testdatadir = os.path.join(os.path.dirname(__file__),
                                        'testdata')
        import tempfile

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
                
        self.jsonfiles = [tmpfn1, tmpfn2, tmpfn3]

        
    def test_new_file(self):
        from mosbot import main
        import tempfile

        tempdir = tempfile.mkdtemp()
        args = ['--script', os.path.join(tempdir, 'tonight.sh')]
        args += self.jsonfiles
        mosbot = main(cmdlineargs=args, get_mosbot=True)

        # write sequence number
        f = open(mosbot.seqnumpath, 'w')
        f.write('10\n')
        f.close()
        
        fn = os.path.join(self.testdatadir, 'mos3.68488.im4.fits.fz')
        mosbot.process_file(fn)


    def test_new_file_blank(self):
        from mosbot import main
        import tempfile

        tempdir = tempfile.mkdtemp()
        args = ['--script', os.path.join(tempdir, 'tonight.sh')]
        args += self.jsonfiles
        mosbot = main(cmdlineargs=args, get_mosbot=True)
        mosbot.J1 = mosbot.J1[:10]
        mosbot.J2 = mosbot.J2[:10]
        mosbot.J3 = mosbot.J3[:10]
        mosbot.Nahead = 15
        
        # write sequence number
        f = open(mosbot.seqnumpath, 'w')
        f.write('5\n')
        f.close()
        
        fn = os.path.join(self.testdatadir, 'mos3.68488.im4.fits.fz')
        #mosbot.process_file(fn)
        mosbot.update_for_image({
            'airmass' :1.019, 'extension': 'im4', 'pixscale': 0.26,
            'camera': 'mosaic3', 'ra_ccd': 102.97109434817698,
            'dec_ccd': 37.53083859605576, 'band': 'z',
            'zp': 26.348159194305978, 'rawsky': 3522.2537, 'ndetected': 571,
            'skybright': 19.372964228396995, 'dy': -50.904994833603581,
            'transparency': 0.94911683705518535, 'seeing': 0.9776841626480679,
            'dx': 9.0965935687085278, 'nmatched': 199})


    def test_large_slew(self):
        from mosbot import main
        import tempfile

        tempdir = tempfile.mkdtemp()
        args = ['--script', os.path.join(tempdir, 'tonight.sh')]
        args += self.jsonfiles
        mosbot = main(cmdlineargs=args, get_mosbot=True)
        mosbot.J1[0]['RA'] = 0
        mosbot.J1[0]['dec'] = -90
        
        # write sequence number
        f = open(mosbot.seqnumpath, 'w')
        f.write('5\n')
        f.close()
        
        fn = os.path.join(self.testdatadir, 'mos3.68488.im4.fits.fz')
        #mosbot.process_file(fn)
        mosbot.update_for_image({
            'airmass': 1.019,
            'skybright': 19.372964228396995,
            'transparency': 0.94911683705518535,
            'seeing': 0.9776841626480679,
            'band': 'z'})

        
    def XXXtest_run(self):
        from mosbot import main

        import tempfile
        tempdir = tempfile.mkdtemp()
        args = ['--script', os.path.join(tempdir, 'tonight.sh')]
        args += self.jsonfiles

        mosbot = main(cmdlineargs=args, get_mosbot=True)

        # write sequence number
        f = open(mosbot.seqnumpath, 'w')
        f.write('10\n')
        f.close()

        mosbot.run()

        
if __name__ == '__main__':
    unittest.main()
