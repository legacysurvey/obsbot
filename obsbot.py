from __future__ import print_function

import os
import datetime
from collections import Counter

import numpy as np

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

    def fiducial_exptime(self, band):
        '''
        Returns an object with attributes:

        - skybright
        - k_co
        - A_co
        - seeing
        - exptime, exptime_min, exptime_max

        
        '''
        fid = NominalExptime()

        if band == 'g':
            fid.update(
                exptime     =  50.,
                exptime_max = 125.,
                exptime_min =  40.,
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
        # print('Predicted sky flux per pixel per second: %.1f' %skyflux)
        t_sat = 30000. / skyflux
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

def mjdnow():
    from astrometry.util.starutil_numpy import datetomjd
    return datetomjd(datenow())

class NewFileWatcher(object):
    def __init__(self, dir, backlog=True):
        self.dir = dir

        # How many times to re-try processing a new file
        self.maxFail = 10

        self.sleeptime = 5.

        # How often to call the timeout function -- this is the time
        # since the last new file was seen, OR since the last time the
        # timeout was called.
        self.timeout = 60.

        # Get current file list
        files = set(os.listdir(self.dir))

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
    
    def get_new_files(self):
        files = set(os.listdir(self.dir))
        newfiles = list(files - self.oldfiles)
        newfiles = self.filter_new_files(newfiles)
        return newfiles

    def get_newest_file(self):
        newfiles = self.get_new_files()
        if len(newfiles) == 0:
            return None
        # Take the one with the latest timestamp.
        latest = None
        newestfile = None
        for fn in newfiles:
            st = os.stat(os.path.join(self.dir, fn))
            t = st.st_mtime
            if latest is None or t > latest:
                newestfile = fn
                latest = t
        return newestfile

    def try_open_file(self, path):
        pass
    
    def run_one(self):
        fn = self.get_newest_file()
        if fn is None:
            if self.timeout is None:
                return False
            # Check timeout
            now = datenow()
            dt = (now - self.lastTimeout).total_seconds()
            if dt < self.timeout:
                return False
            self.timed_out(dt)
            self.lastTimout = datenow()
            if len(self.backlog) == 0:
                return False
            fn = self.backlog.pop()

        if self.failCounter[fn] >= self.maxFail:
            print('Failed to process file: %s, %i times.  Ignoring.' %
                  (fn, self.maxFail))
            self.oldfiles.add(fn)
            return False
            
        path = os.path.join(self.dir, fn)
        print('Found new file:', path)
        try:
            self.try_open_file(path)
        except:
            print('Failed to open %s: maybe not fully written yet.' % path)
            self.failCounter.update([fn])
            return False

        try:
            self.process_file(path)
            self.oldfiles.add(fn)
            self.processed_file(path)
            self.lastNewFile = self.lastTimeout = datenow()
            return True
            
        except IOError:
            print('Failed to read file: %s' % path)
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
    
    
