from __future__ import print_function
import unittest
#from django.test import TestCase
import os

class TestMosbot(unittest.TestCase):

    def test_maybe_ha_flip(self):
        from mosbot import maybe_ha_flip
        self.assertEqual(maybe_ha_flip(86, 86, 3), False)
        self.assertEqual(maybe_ha_flip(86, 88, 3), True)
        self.assertEqual(maybe_ha_flip(88, 88, 3), True)
        self.assertEqual(maybe_ha_flip(88, 86, 3), True)
        self.assertEqual(maybe_ha_flip(90, 90, 3), True)
        self.assertEqual(maybe_ha_flip(90, 94, 3), True)
        self.assertEqual(maybe_ha_flip(94, 94, 3), False)
        self.assertEqual(maybe_ha_flip(94, 90, 3), True)
        self.assertEqual(maybe_ha_flip(90, 86, 3), True)

        self.assertEqual(maybe_ha_flip(-94, 94, 3), True)
        
        self.assertEqual(maybe_ha_flip(-86, -86, 3), False)
        self.assertEqual(maybe_ha_flip(-86, -88, 3), True)
        self.assertEqual(maybe_ha_flip(-88, -88, 3), True)
        self.assertEqual(maybe_ha_flip(-88, -86, 3), True)
        self.assertEqual(maybe_ha_flip(-90, -90, 3), True)
        self.assertEqual(maybe_ha_flip(-90, -94, 3), True)
        self.assertEqual(maybe_ha_flip(-94, -94, 3), False)
        self.assertEqual(maybe_ha_flip(-94, -90, 3), True)
        self.assertEqual(maybe_ha_flip(-90, -86, 3), True)
            
    
    def set_json_files(self, start_time):
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
            for j in J:
                tnew = start_time + ephem.Date(str(j['approx_datetime'])) - t0
                tnew = str(ephem.Date(tnew))
                j['approx_datetime'] = tnew
                print('Updated datetime to', tnew)

            f = open(tmpfn, 'w')
            json.dump(J, f, sort_keys=True, indent=4, separators=(',', ': '))
            f.close()
        self.jsonfiles = [tmpfn1, tmpfn2, tmpfn3]

    def setUp(self):
        self.testdatadir = os.path.join(os.path.dirname(__file__),
                                        'testdata')
        import ephem
        now = ephem.now()
        print('Now:', now)
        self.set_json_files(now)
        
    def test_new_file(self):
        from mosbot import main
        import tempfile

        tempdir = tempfile.mkdtemp()
        args = ['--script', os.path.join(tempdir, 'tonight.sh')]
        args.append('--adjust')
        #args.append('--no-db')
        args += self.jsonfiles
        mosbot = main(cmdlineargs=args, get_mosbot=True)

        # write sequence number
        f = open(mosbot.seqnumpath, 'w')
        f.write('10\n')
        f.close()
        
        fn = os.path.join(self.testdatadir, 'mos3.68488.im4.fits.fz')
        mosbot.process_file(fn)

        # Assert that entries have been added to the computed-exposure database
        import obsdb
        from obsdb.models import ComputedExptime, OtherPasses
        c = ComputedExptime.objects.all()
        for cc in c:
            print(cc)
        self.assertEqual(len(c), 19)

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
            'zp': 26.348159194305978,
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

    def test_no_cut_past(self):
        from mosbot import main
        import tempfile
        import ephem

        now = ephem.now()
        print('Now:', now)
        # Make the plans start an hour ago.
        now -= (1./24.)
        self.set_json_files(now)

        tempdir = tempfile.mkdtemp()
        args = ['--script', os.path.join(tempdir, 'tonight.sh'),
                '--no-cut-past', '--pass', '1']
        args += self.jsonfiles
        mosbot = main(cmdlineargs=args, get_mosbot=True)

        # Touch forcepass1
        f = open(os.path.join(tempdir, 'forcepass1'), 'w')
        f.close()
        
        tiles = [650124, 650123, 652527, 654902, 652528,
                 654903, 657247, 650125, 657246, 657248]
        
        # BEFORE: check tile numbers
        for i,tile in enumerate(tiles):
            fn = os.path.join(tempdir, 'expose-%i.sh' % (i+1))
            txt = open(fn).read()
            tilename = 'MzLS_%i_z' % tile
            print('Checking for tile', tilename, 'in file', fn)
            assert(tilename in txt)
        
        # write sequence number
        f = open(mosbot.seqnumpath, 'w')
        f.write('5\n')
        f.close()

        mosbot.update_for_image({
            'airmass': 1.019,
            'skybright': 19.372964228396995,
            'transparency': 0.94911683705518535,
            'zp': 26.348159194305978,
            'seeing': 0.9776841626480679,
            'band': 'z'})

        # AFTER: check that tile numbers remain the same!
        for i,tile in enumerate(tiles):
            fn = os.path.join(tempdir, 'expose-%i.sh' % (i+1))
            txt = open(fn).read()
            tilename = 'MzLS_%i_z' % tile
            print('Checking for tile', tilename, 'in file', fn)
            assert(tilename in txt)
        
        
if __name__ == '__main__':
    import obsdb
    from camera_mosaic import database_filename
    obsdb.django_setup(database_filename=database_filename)

    unittest.main()
