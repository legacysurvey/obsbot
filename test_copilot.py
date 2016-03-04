import unittest
from django.test import TestCase
import os
import datetime

class TestCopilot(TestCase):

    def setUp(self):
        self.testdatadir = os.path.join(os.path.dirname(__file__),
                                        'testdata')

        self.mos_args = ['--ext', 'im4', '--tiles', 'obstatus/mosaic-tiles_obstatus.fits', '--no-show']
        
    def test_cmdline_args(self):
        from copilot import main, skip_existing_files
        
        from obsdb.models import MeasuredCCD
        fn = os.path.join(self.testdatadir, 'mos3.68488.im4.fits.fz')

        args = self.mos_args + ['--n-fwhm', '10']
        
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

        # add --ext im4 so that at DECam site (where default extension
        # is different) this still works... is this exposing dumbness
        # of site-specific extension defaults?
        args = self.mos_args + ['--n-fwhm', '10']
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


    def test_focus(self):
        from copilot import main, Copilot
        import tempfile
        import time
        
        # Test reading a new image from $MOS3_DATA
        dirname = tempfile.mkdtemp()
        os.environ['MOS3_DATA'] = dirname

        args = self.mos_args + ['--n-fwhm', '10']
        copilot = main(cmdlineargs=args, get_copilot=True)

        # Create symlink to focus frame
        fn = 'mos3.63127.im4.fits.fz'
        path = os.path.join(self.testdatadir, fn)
        sym = os.path.join(dirname, fn)
        os.symlink(path, sym)

        copilot.run_one()
        
        # Clean up
        os.unlink(sym)
        os.rmdir(dirname)

    def test_longtime(self):
        from copilot import main, Copilot
        from astrometry.util.fits import fits_table
        from astrometry.util.starutil_numpy import mjdtodate
        from obsdb.models import MeasuredCCD
        from obsbot import mjdnow
        import tempfile
        
        dirname = tempfile.mkdtemp()
        os.environ['MOS3_DATA'] = dirname
        
        args = self.mos_args + ['--n-fwhm', '10']
        args += ['--plot-filename', 'longtime.png']
        copilot = main(cmdlineargs=args, get_copilot=True)
        #copilot.exp

        T = fits_table(os.path.join(self.testdatadir, 'mosaic-db.fits'))

        # Move all the exposures' times so that the latest one is now-longtime
        mjdoffset = (mjdnow() - (copilot.longtime + 5.)/86400.) - max(T.mjd_obs)
        T.mjd_obs += mjdoffset
        
        for t in T:
            m,created = MeasuredCCD.objects.get_or_create(
                filename=t.filename, extension=t.extension)
            for c in T.get_columns():
                setattr(m, c, t.get(c))
            m.object = m.object.strip()
            m.obstype = m.obstype.strip()
            m.save()

        copilot.lastNewFile = mjdtodate(max(T.mjd_obs))
        copilot.plot_recent()

        # Change the mjd_obs times so that it's *not* longtime.
        copilot.opt.plot_filename = 'longtime2.png'
        
        # T.mjd_obs += 10./86400.
        # for t in T:
        #     m,created = MeasuredCCD.objects.get_or_create(
        #         filename=t.filename, extension=t.extension)
        #     m.mjd_obs = t.mjd_obs
        #     m.save()
        copilot.lastNewFile += datetime.timedelta(0, 10.)
        copilot.plot_recent()

        print('Longtime:', copilot.longtime)
        
        
if __name__ == '__main__':
    unittest.main()
    
