import unittest
from django.test import TestCase
import os

class TestCopilot(TestCase):

    def setUp(self):
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
        from copilot import main, Copilot
        import tempfile
        import time
        
        # Test reading a new image from $MOS3_DATA
        dirname = tempfile.mkdtemp()
        os.environ['MOS3_DATA'] = dirname

        args = ['--n-fwhm', '10']
        copilot = main(cmdlineargs=args, get_copilot=True)

        def fake_process_file(self, fn):
            self.processed.append(fn)

        def fake_plot_recent(self):
            self.plotted = True
            
        Copilot.process_file = fake_process_file
        Copilot.plot_recent  = fake_plot_recent
        copilot.processed = []
        copilot.plotted = False
        
        # No new files, no backlog, no timeout -- nothing should happen.
        copilot.run_one()

        self.assertEqual(len(copilot.processed), 0)
        self.assertEqual(copilot.plotted, False)

        copilot.processed = []
        copilot.plotted = False

        # Create symlink
        fn = 'mos3.68488.im4.fits.fz'
        path = os.path.join(self.testdatadir, fn)
        sym = os.path.join(dirname, fn)
        os.symlink(path, sym)

        copilot.run_one()
        
        self.assertEqual(len(copilot.processed), 1)
        self.assertEqual(copilot.plotted, True)

        copilot.processed = []
        copilot.plotted = False

        # Create another copy of the symlink, wait one second, then
        # create empty file -- Copilot should read the empty one
        # first, but fail to read "maxFail" times, then ignore it and
        # go on to read the symlink.

        # copilot remembers the first image processed above.
        self.assertEqual(len(copilot.oldfiles), 1)

        sym2 = os.path.join(dirname, 'x' + fn)
        os.symlink(path, sym2)
        time.sleep(1.)
        
        f,tmpfn = tempfile.mkstemp(dir=dirname, suffix='.fits')
        os.close(f)

        for i in range(copilot.maxFail):
            copilot.run_one()
            self.assertEqual(len(copilot.processed), 0)
            self.assertEqual(copilot.plotted, False)

        # AFTER the maxFail'th time, the bad file gets added to
        # 'oldimages', and (still) doesn't get processed.
        copilot.run_one()
        self.assertEqual(len(copilot.oldfiles), 2)
        self.assertEqual(len(copilot.processed), 0)
        self.assertEqual(copilot.plotted, False)

        # The time after THAT, the symlink should get processed.
        copilot.run_one()
        self.assertEqual(len(copilot.oldfiles), 3)
        self.assertEqual(len(copilot.processed), 1)
        self.assertEqual(copilot.plotted, True)
        self.assertEqual(copilot.processed[0], sym2)
        
        # Clean up
        os.unlink(sym2)
        os.unlink(sym)
        os.unlink(tmpfn)
        os.rmdir(dirname)

if __name__ == '__main__':
    unittest.main()
    
