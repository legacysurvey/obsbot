from __future__ import print_function

import os
import datetime
from collections import Counter
import time

import numpy as np

def choose_pass(trans, seeing, skybright, nomsky,
                forcedir=''):

    brighter = nomsky - skybright

    transcut = 0.9
    seeingcut = 1.25
    brightcut = 0.25

    transcut2 = 0.7
    seeingcut2 = 2.0
    
    transok = trans > transcut
    seeingok = seeing < seeingcut
    brightok = brighter < brightcut

    transfair = trans > transcut2
    seeingfair = seeing < seeingcut2

    trans_txt = 'good' if transok else ('fair' if transfair else 'poor')
    seeing_txt = 'good' if seeingok else ('fair' if seeingfair else 'poor')

    pass1ok = transok and seeingok and brightok
    pass2ok = (transok and seeingfair) or (seeingok and transfair)
    
    print('Transparency: %s       (%6.2f vs %6.2f = good, %6.2f = fair)' %
          (trans_txt, trans, transcut, transcut2))
    print('Seeing      : %s       (%6.2f vs %6.2f = good, %6.2f = fair)' %
          (seeing_txt, seeing, seeingcut, seeingcut2))
    print('Brightness  : %s       (%6.2f vs %6.2f = pass)' %
          (('pass' if brightok else 'fail'), skybright, nomsky+brightcut))
    print('Pass 1 = transparency AND seeing AND brightness: %s' % pass1ok)
    print('Pass 2 = (transparency good and seeing fair) OR (seeing good and transparency fair): %s' % pass2ok)

    p,fn = get_forced_pass(forcedir)
    if p is not None:
        print('Forcing pass %i because file exists: %s' % (p, fn))
        return p
    
    path = os.path.join(forcedir, 'nopass1')
    if os.path.exists(path):
        print('File %s exists; not allowing pass 1' % path)
        pass1ok = False
        
    if pass1ok:
        return 1
    if pass2ok:
        return 2
    return 3

def get_forced_pass(forcedir=''):
    for p in [1,2,3]:
        path = os.path.join(forcedir, 'forcepass%i' % p)
        #print('Checking for file "%s"' % path)
        if os.path.exists(path):
            #print('Forcing pass %i because file exists: %s' % (p, path))
            return p, path
    return None, None
    

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
                    z = 0.194)[band]
    
    def fiducial_exptime(self, band):
        '''
        Returns an object with attributes:

        - skybright
        - k_co
        - A_co
        - seeing
        - exptime, exptime_min, exptime_max

        
        '''
        if not band in 'grz':
            return None
        fid = NominalExptime()

        # 2-coverage targets (90% fill), 5-sigma extinction-corrected
        # canonical galaxy detection.
        target_depths = dict(g=24.0, r=23.4, z=22.5)

        target_depth = target_depths[band]
        # -> 1-coverage depth (- ~0.37 mag)
        target_depth -= 2.5*np.log10(np.sqrt(2.))

        fid.update(single_exposure_depth = target_depth)
        
        if band == 'g':
            fid.update(
                exptime     =  70.,
                exptime_max = 175.,
                exptime_min =  56.,
                )

        elif band == 'r':
            fid.update(
                exptime     =  50.,
                exptime_max = 125.,
                exptime_min =  40.,
                )

        elif band == 'z':
            fid.update(
                exptime     = 100.,
                exptime_max = 250.,
                exptime_min =  80.,
                )
        else:
            raise ValueError('Unknown band "%s"' % band)

        fid.update(
            skybright = self.sky(band),
            seeing = 1.3,
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

    r_half = 0.45 #arcsec
    ps = cal.pixscale

    def Neff(seeing):
        # magic 2.35: convert seeing FWHM into sigmas in arcsec.
        return (4. * np.pi * (seeing / 2.35)**2 +
                8.91 * r_half**2 +
                ps**2/12.)

    # Nightlystrategy.py has:
    # pfact = 1.15
    # Neff_fid = ((4.0*np.pi*sig_fid**2)**(1.0/pfact)+(8.91*r_half**2)**(1.0/pfact))**pfact
    
    neff_fid = Neff(fid.seeing)
    neff     = Neff(seeing)

    # print('exposure_factor:')
    # print('Transparency:', transparency)
    # print('airmass:', airmass)
    # print('ebv:', ebv)
    # print('seeing:', seeing)
    # print('sky:', skybright)
    # print('neff:', neff, 'fid', neff_fid)
    
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


def get_tile_from_name(name, tiles):
    # Parse objname like 'MzLS_5623_z'
    parts = name.split('_')
    ok = (len(parts) == 3)
    if ok:
        band = parts[2]
        ok = ok and (band in 'grz')
    if not ok:
        return None
    try:
        tileid = int(parts[1])
    except:
        return None
    # Find this tile in the tiles table.
    I = np.flatnonzero(tiles.tileid == tileid)
    assert(len(I) == 1)
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
        import StringIO
        uniq = kwargs.pop('uniq', False)
        f = StringIO.StringIO()
        print(*args, file=f, **kwargs)
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
        files = set(os.listdir(self.dir))
        return [os.path.join(self.dir, fn) for fn in files]
            
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
    



# Code shared between mosbot.py and decbot.py
class Obsbot(NewFileWatcher):
    def adjust_for_previous(self, tile, band, fid, debug=False,
                            get_others=False):
        '''
        Adjust the exposure time we should take for this image based
        on data we've already taken.
        '''
        # Find the other passes for this tile, and if we've taken an
        # exposure, think about adjusting the exposure time.  If the
        # depth in the other exposure is more than THRESHOLD too
        # shallow, ignore it on the assumption that we'll re-take it.
        # If there is a previous exposure for THIS tile, reduce our
        # exposure time accordingly.

        # Find other passes
        others = self.other_passes(tile, self.tiles)
        others.depth = others.get('%s_depth' % band)

        target = fid.single_exposure_depth
        if debug:
            print('Adjusting exposure for tile', tile.tileid, 'pass',
                  tile.get('pass'))

        # depth = 0 means no obs;
        # depth = 1 means non-photometric observation was taken.
        # depth = 30 means image was taken but we don't know how deep it is yet

        I = np.flatnonzero((others.depth > 1) * (others.depth < 30))
        if len(I) == 0:
            if debug:
                print('No other passes have measured depths')
            if get_others:
                return 1.0,others
            return 1.0
        if debug:
            print('Other tile passes:', others.get('pass')[I])
            print('Other tile depths:', others.get('%s_depth' % band)[I])
            print('Target depth:', target)

        thisfactor = 1.0

        threshold = 0.25

        shortfall = target - others.depth
        I = np.flatnonzero((others.depth > 1) *
                           (shortfall > 0) * (shortfall < threshold))
        if len(I) > 0:
            others.cut(I)
            if debug:
                print('Other tiles with acceptable shortfalls:', shortfall[I])
            # exposure time factors required
            factors = (10.**(-shortfall[I] / 2.5))**2
            if debug:
                print('Exposure time factors:', factors)
            extra = np.sum(1 - factors)
            if debug:
                print('Total extra fraction required:', extra)
            # Split this extra required exposure time between the remaining
            # passes...
            # Magic number 3 = three passes
            nremain = max(1, 3 - len(I))
            if debug:
                print('Extra time to be taken in this image:', extra / nremain)
            thisfactor += extra / nremain
        else:
            if debug:
                print('All other tiles reached depth or need to be retaken.')
            
        depth = tile.get('%s_depth' % band)
        if depth > 1:
            # If this tile has had previous exposure(s), subtract that.
            shortfall = target - depth
            factor = (10.**(-shortfall / 2.5))**2
            if debug:
                print('This tile had previous depth:', depth)
                print('Fraction of nominal exposure time:', factor)
            thisfactor -= factor

        if debug:
            print('Exposure time factor based on previous exposures:',
                  thisfactor)

        if get_others:
            return thisfactor,others
        return thisfactor

    def other_passes(self, tile, tiles):
        '''
        Given tile number *tile*, return the obstatus rows for the other passes
        on this tile center.

        Returns: *otherpasses*, table object
        '''
        from astrometry.libkd.spherematch import match_radec
        # Could also use the fact that the passes are offset from each other
        # by a fixed index (15872 for decam)...
        # Max separation is about 0.6 degrees for DECam...
        #### FIXME for Mosaic this is much smaller... and the three passes
        #### for a tile are not necessarily relatively close to each other.
        I,J,d = match_radec(tile.ra, tile.dec, tiles.ra, tiles.dec, 1.)
        # Omit 'tile' itself...
        K = np.flatnonzero(tiles.tileid[J] != tile.tileid)
        J = J[K]
        return tiles[J]
        
