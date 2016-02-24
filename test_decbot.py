from __future__ import print_function
import unittest
import os

class TestDecbot(unittest.TestCase):

    def setUp(self):
        self.testdatadir = os.path.join(os.path.dirname(__file__),
                                        'testdata')

    def test_new_file(self):
        from decbot import main
        import tempfile

        fn1 = os.path.join(self.testdatadir, 'decals-pass1.json')
        fn2 = os.path.join(self.testdatadir, 'decals-pass2.json')
        fn3 = os.path.join(self.testdatadir, 'decals-pass3.json')

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
                
        args = [tmpfn1, tmpfn2, tmpfn3]

        args += ['--remote-server', 'localhost']
        args += ['--remote-port',   '7767']
        
        decbot = main(cmdlineargs=args, get_decbot=True)

        fn = os.path.join(self.testdatadir, 'decam-00488199-n4.fits.fz')
        decbot.process_file(fn)

        decbot.run()
            

if __name__ == '__main__':
    unittest.main()
