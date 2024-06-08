from __future__ import print_function

import os
import datetime
from collections import Counter
import time
import sys
import numpy as np

py2 = (sys.version_info[0] == 2)

class NominalExptime(object):
    def update(self, **kwargs):
        for k,v in kwargs.items():
            setattr(self, k, v)

class NominalCalibration(object):
    '''
    Overridden or instantiated by Mosaic / DECam nominal calibrations.

    Attributes (must) include:

    - pixscale -- in arcsec/pixel
    - overhead -- in seconds
    '''
    def zeropoint(self, band, ext=None):
        pass

    def sky(self, band):
        pass

    def cdmatrix(self, ext):
        pass

    def seeing_wrt_airmass(self, band):
        # From email from Arjun, 2016-08-03 "Scaling for g-band exposure times"
        return dict(g = 0.422,
                    r = 0.345,
                    z = 0.194,
                    zd = 0.194,
                    D51 = 0.422,
                    M411 = 0.422,
                    M464 = 0.422,
                    )[band]

    def fiducial_exptime(self, band):
        '''
        Returns an object with attributes:

        - skybright
        - k_co
        - A_co
        - seeing
        - exptime, exptime_min, exptime_max
        '''
        if not band in ['u',
                        'g',
                        'r',
                        'i',
                        'z','zd',
                        'D51',
                        'N419',
                        'N501',
                        'N540',
                        'N673',
                        'N708',
                        'M411',
                        'M464',
                ]:
            return None
        fid = NominalExptime()

        # 2-coverage targets (90% fill), 5-sigma extinction-corrected
        # canonical galaxy detection.
        target_depths = dict(u=24.0, g=24.0, r=23.4, i=23.0, z=22.5,
                             zd=22.5, D51=24.0,
                             N419=24.5,
                             N501=24.0,
                             N540=23.4,
                             N673=23.4,
                             N708=23.4,
                             # 25.0 after 6 passes (this is 25.0 - 2.5*np.log10(np.sqrt(3)), because of the extra sqrt(2) below!)
                             M411=24.40,
                             M464=24.40,
                             )

        target_depth = target_depths[band]
        # -> 1-coverage depth (- ~0.37 mag)
        target_depth -= 2.5*np.log10(np.sqrt(2.))

        fid.update(single_exposure_depth = target_depth)

        if band == 'u':
            fid.update(
                exptime     =  70.,
                exptime_max = 200.,
                exptime_min =  56.,
                )

        elif band == 'g':
            fid.update(
                exptime     =  70.,
                exptime_max = 200.,
                exptime_min =  56.,
                )

        elif band == 'r':
            fid.update(
                exptime     =  50.,
                exptime_max = 175.,
                exptime_min =  40.,
                )

        elif band == 'i':
            fid.update(
                exptime     =  50.,
                exptime_max = 175.,
                exptime_min =  40.,
                )

        elif band in ['z', 'zd']:
            fid.update(
                exptime     = 100.,
                exptime_max = 250.,
                exptime_min =  80.,
                )

        elif band == 'D51':
            # we're not updating exposure times, but hey
            fid.update(
                exptime     = 600.,
                exptime_max = 600.,
                exptime_min = 200.,
                )

        elif band in ['N419', 'N501', 'N540', 'N673', 'N708']:
            # we're not updating exposure times in ODIN or MERIAN, but hey
            fid.update(
                exptime     = 900.,
                exptime_max = 900.,
                exptime_min = 300.,
                )

        elif band in ['M411', 'M464']:
            # IBIS
            fid.update(
                exptime     = 300.,
                exptime_max = 500.,
                exptime_min = 200.,
                )
        else:
            raise ValueError('Unknown band "%s"' % band)

        fid.update(
            skybright = self.sky(band),
            seeing = 1.25,
            )

        # Camera-specific update:
        fid = self._fiducial_exptime(fid, band)

        return fid

    def _fiducial_exptime(self, fid, band):
        return fid

    def saturation_time(self, band, skybright):
        skyflux = 10. ** ((skybright - self.zeropoint(band)) / -2.5)
        skyflux *= self.pixscale**2
        # print('Predicted sky flux per pixel per second: %.1f electrons' %skyflux)
        # Divide by (nominal) gain to get from electrons back to ADU
        skyflux /= self.gain
        # 30k for DECam, 20k for Mosaic3
        t_sat = self.saturation_adu / skyflux
        return t_sat

def Neff(seeing, pixscale):
    r_half = 0.45 #arcsec
    # magic 2.35: convert seeing FWHM into sigmas in arcsec.
    return (4. * np.pi * (seeing / 2.35)**2 +
            8.91 * r_half**2 +
            pixscale**2/12.)

# From Anna Patej's nightlystrategy / mosaicstrategy
def exposure_factor(fid, cal,
                    airmass, ebv, seeing, skybright, transparency):
    '''
    Computes a factor by which the exposure time should be scaled
    relative to nominal.

    *fid*: fiducial exposure time properties
    *cal*: nominal camera calibration properties
    *airmass*: airmass, float > 1
    *ebv*: extinction E(B-V) mags
    *seeing*: FWHM in arcsec
    *skybright*: sky brightness
    *transparency*: sky transparency

    Returns:
    scaling: exposure time scale factor,
      scaling = T_new/T_fiducial
    '''

    ps = cal.pixscale

    # Nightlystrategy.py has:
    # pfact = 1.15
    # Neff_fid = ((4.0*np.pi*sig_fid**2)**(1.0/pfact)+(8.91*r_half**2)**(1.0/pfact))**pfact

    neff_fid = Neff(fid.seeing, ps)
    neff     = Neff(seeing, ps)

    # print('exposure_factor:')
    # print('Transparency:', transparency)
    # print('  -> factor', 1./transparency**2)
    # print('airmass:', airmass)
    # print('  -> factor', 10.**(0.8 * fid.k_co * (airmass - 1.)))
    # print('ebv:', ebv)
    # print('  -> factor', 10.**(0.8 * fid.A_co * ebv))
    # print('seeing:', seeing)
    # # print('neff:', neff, 'fid', neff_fid)
    # print('  -> factor', (neff / neff_fid))
    # print('sky:', skybright)
    # print('  -> factor', 10.**(-0.4 * (skybright - fid.skybright)))

    scaling = (1./transparency**2 *
               10.**(0.8 * fid.k_co * (airmass - 1.)) *
               10.**(0.8 * fid.A_co * ebv) *
               (neff / neff_fid) *
               10.**(-0.4 * (skybright - fid.skybright)))
    return scaling

# From Anna Patej's nightlystrategy / mosaicstrategy
def get_airmass(alt):
    if (alt < 0.07):
        alt = 0.07
    secz = 1.0/np.sin(alt)
    seczm1 = secz-1.0
    airm = secz-0.0018167*seczm1-0.002875*seczm1**2-0.0008083*seczm1**3
    return airm

def read_tiles_file(fn):
    from astrometry.util.fits import fits_table
    if fn.endswith('.fits'):
        tiles = fits_table(fn)
    else:
        # ECSV... FIXME... read with astropy, convert to astrometry fits_table
        from astropy.table import Table
        tiles = fits_table()
        t = Table.read(fn)
        colnames = list(t.columns)
        for c in colnames:
            col = t[c]
            tiles.set(c.lower(), col)
        del t
    return tiles

def get_tile_id_from_name(name):
    # Parse objname like 'MzLS_5623_z' / 'IBIS_wide_M411_100100'
    parts = str(name).split('_')
    if len(parts) == 3:
        # OLD DECaLS / MzLS
        band = parts[2]
        if not band in 'grz':
            return None
        tilestr = parts[1]

    elif len(parts) == 4:
        # IBIS
        tilestr = parts[-1]

    try:
        tileid = int(tilestr)
    except:
        return None
    return tileid

def get_tile_from_name(name, tiles):
    tileid = get_tile_id_from_name(name)
    if tileid is None:
        return None
    # Find this tile in the tiles table.
    I = np.flatnonzero(tiles.tileid == tileid)
    if len(I) != 1:
        print('Warning: tile file contains', len(I), 'matches for tile name', name, '-> tile ID', tileid, '-- expected only 1 match.')
        return None
    tile = tiles[I[0]]
    return tile

def datenow():
    return datetime.datetime.utcnow()

# For testing purposes:
mjdnow_offset = 0.
def mjdnow():
    from astrometry.util.starutil_numpy import datetomjd
    return datetomjd(datenow()) + mjdnow_offset

def unixtime_to_ephem_date(unixtime):
    import ephem
    # Unix time is seconds since 1970/1/1 00:00:00, UT.
    # 'unixepoch' here is that date in ephem.Date format.
    unixepoch = 25567.5
    # ephem.Date counts in days, hence the 86400.
    return ephem.Date(unixepoch + unixtime / 86400.)

def ephem_date_to_mjd(edate):
    from astrometry.util.starutil_numpy import datetomjd
    return datetomjd(edate.datetime())

def mjd_to_ephem_date(mjd):
    import ephem
    # MAGIC ephem.Date(datetime.datetime(1858, 11, 17, 0, 0, 0))
    return ephem.Date(mjd -15019.5)

class Logger(object):
    def __init__(self, verbose=False, timestamp=True):
        self.verbose = verbose
        self.timestamp = timestamp
        self.last_printed = None

    def log(self, *args, **kwargs):
        '''
        Keyword args:
        uniq: if True, do not print a log message if it repeats the last printed
        log message.
        kwargs: passed to print().
        '''
        from io import StringIO, BytesIO
        uniq = kwargs.pop('uniq', False)
        if py2:
            f = BytesIO()
        else:
            f = StringIO()
        pkw = kwargs.copy()
        pkw.update(file=f)
        print(*args, **pkw)
        #print(*args, file=f, **kwargs)
        s = f.getvalue()
        if uniq and s == self.last_printed:
            return
        self.last_printed = s
        if self.timestamp:
            import ephem
            now = str(ephem.now())
            print('%s: %s' % (now, s), end='')
        else:
            print(s, end='')

    def debug(self, *args, **kwargs):
        if self.verbose:
            self.log(*args, **kwargs)

class NewFileWatcher(Logger):
    def __init__(self, dir, backlog=True, only_process_newest=False,
                 ignore_missing_dir=False,
                 verbose=False, timestamp=True):
        super(NewFileWatcher, self).__init__(verbose=verbose,
                                             timestamp=timestamp)
        self.dir = dir
        self.only_process_newest = only_process_newest
        self.ignore_missing_dir = ignore_missing_dir

        # How many times to re-try processing a new file
        self.maxFail = 10

        self.sleeptime = 5.

        # How often to call the timeout function -- this is the time
        # since the last new file was seen, OR since the last time the
        # timeout was called.
        self.timeout = 60.

        # Get current file list
        files = self.get_file_list()

        if backlog:
            # (note to self, need explicit backlog because we skip existing
            # for the backlogged files, unlike new ones.)
            self.backlog = self.filter_backlog(files)
            # ... then reset oldfiles to the current file list.
            self.oldfiles = set(files)
        else:
            self.backlog = set()
            self.oldfiles = set(self.filter_new_files(files))

        # initialize timeout counter
        self.lastTimeout = datenow()
        self.lastNewFile = datenow()

        # Keep track of how many times we've failed to process a file...
        self.failCounter = Counter()

    def filter_backlog(self, backlog):
        return self.filter_new_files(backlog)

    def filter_new_files(self, fns):
        return fns

    def timed_out(self, dt):
        pass

    def processed_file(self, path):
        pass

    def get_file_list(self):
        if self.ignore_missing_dir and not os.path.exists(self.dir):
            self.log('Directory', self.dir, 'does not exist -- waiting for it',
                     uniq=True)
            return []
        files = set()
        # Note: does not follow directory symlinks.
        for (dirpath, dirnames, filenames) in os.walk(self.dir):
            for fn in filenames:
                files.add(os.path.join(dirpath, fn))
        return list(files)

    def get_new_files(self):
        files = set(self.get_file_list())
        newfiles = list(files - self.oldfiles)
        newfiles = self.filter_new_files(newfiles)
        return newfiles

    def get_newest_file(self, newfiles=None):
        if newfiles is None:
            newfiles = self.get_new_files()
        if len(newfiles) == 0:
            return None
        # Take the one with the latest timestamp.
        latest = None
        newestfile = None
        for fn in newfiles:
            try:
                st = os.stat(fn)
            except OSError as e:
                print('Failed to stat filename', fn, ':', e)
                continue
            t = st.st_mtime
            if latest is None or t > latest:
                newestfile = fn
                latest = t
        return newestfile

    def try_open_file(self, path):
        pass

    def heartbeat(self):
        pass

    def seen_files(self, fns):
        '''
        The given list of filenames has been seen, ie, will not appear
        as new files.  This can include files from the backlog as they are
        processed.
        '''
        pass

    def saw_new_files(self, fns):
        '''We found new files in the directory we're monitoring.
        Files from the backlog don't get this function called.'''
        pass

    def run_one(self):
        fns = self.get_new_files()
        if len(fns):
            self.saw_new_files(fns)
        fn = self.get_newest_file(newfiles=fns)
        if fn is None:
            if self.timeout is None:
                return False
            # Check timeout
            now = datenow()
            dt = (now - self.lastTimeout).total_seconds()
            if dt > self.timeout:
                self.timed_out(dt)
                self.lastTimeout = datenow()
            if len(self.backlog) == 0:
                return False
            fn = self.backlog.pop()
            print('Popping one file from the backlog: %s -- %i remain' %
                  (fn, len(self.backlog)))

        if self.failCounter[fn] >= self.maxFail:
            print('Failed to process file: %s, %i times.  Ignoring.' %
                  (fn, self.maxFail))
            self.oldfiles.add(fn)
            return False

        self.log('Found new file:', fn)
        try:
            self.try_open_file(fn)
        except:
            self.log('Failed to open %s: maybe not fully written yet.' % fn)
            if self.verbose:
                import traceback
                traceback.print_exc()
            self.failCounter.update([fn])
            return False

        try:
            self.process_file(fn)
            if self.only_process_newest:
                self.oldfiles.update(fns)
                self.seen_files(fns)
            else:
                self.oldfiles.add(fn)
                self.seen_files([fn])
            self.processed_file(fn)
            self.lastNewFile = self.lastTimeout = datenow()
            return True

        except IOError as e:
            self.log('Failed to process file: %s (%s)' % (fn, str(e)))
            if self.verbose:
                import traceback
                traceback.print_exc()
            self.failCounter.update([fn])
            return False

    def run(self):
        print('Checking directory for new files:', self.dir)
        sleep = False
        while True:
            print
            if sleep:
                time.sleep(self.sleeptime)
            gotone = self.run_one()
            sleep = not gotone
            self.heartbeat()
