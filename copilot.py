#! /usr/bin/env python2.7
'''

This script is meant to be run during DECaLS observing.  It waits for
new images to appear, measures their sky brightness, seeing, and
transparency, and advises whether & how to replan.

'''
from __future__ import print_function
try:
    from collections import OrderedDict
except:
    print('Failed to import OrderedDict.  You are probably using python 2.6.  Please re-run with python2.7')
    sys.exit(-1)

import sys
import os
import re
import time
import json
import datetime
from collections import Counter
from glob import glob

import numpy as np

import fitsio
import ephem

from astrometry.util.starutil_numpy import hmsstring2ra, dmsstring2dec, mjdtodate, datetomjd

from nightlystrategy import ExposureFactor, getParserAndGlobals, setupGlobals

from measure_raw import measure_raw, get_nominal_cal, get_default_extension, camera_name

from tractor.sfd import SFDMap

def datenow():
    return datetime.datetime.utcnow()

def mjdnow():
    return datetomjd(datenow())

def db_to_fits(mm):
    from astrometry.util.fits import fits_table
    T = fits_table()
    for field in ['filename', 'extension', 'expnum', 'exptime', 'mjd_obs',
                  'airmass', 'racenter', 'deccenter', 'rabore', 'decbore',
                  'band', 'ebv', 'zeropoint', 'transparency', 'seeing',
                  'sky', 'expfactor', 'camera', 'dx', 'dy', 'nmatched',
                  'md5sum', 'bad_pixcnt', 'readtime',
                  'obstype',
                  'object', 'tileid', 'passnumber', 'tileebv']:
        g = getattr(mm[0], field)
        if isinstance(g, basestring):
            T.set(field, np.array([str(getattr(m, field)) for m in mm]))
        else:
            T.set(field, np.array([getattr(m, field) for m in mm]))
    return T

def get_twilight(camera, date):
    '''
    Returns tuple,
    
    (sunset, -12 eve, -18 eve, -18 morn, -12 morn, sunrise)

    for the given camera ("mosaic3" or "decam"), following the night
    whose sunset starts BEFORE the given date.

    date: an ephem.Date object (in UTC)
    '''
    ## From nightlystrategy / mosaicstrategy
    R_earth = 6378.1e3 # in meters
    if camera == 'mosaic3':
        obs = ephem.Observer()
        obs.lon = '-111.6003'
        obs.lat = '31.9634'
        obs.elev = 2120.0 # meters
        #print('Assuming KPNO')

    elif cam == 'decam':
        obs = ephem.Observer()
        obs.lon = '-70.806525'
        obs.lat = '-30.169661'
        obs.elev = 2207.0 # meters
        #print('Assuming CTIO')

    else:
        raise RuntimeError('Unknown camera "%s"' % camera)

    obs.temp = 10.0 # deg celsius; average temp for August
    obs.pressure = 780.0 # mbar
    obs.horizon = -np.sqrt(2.0*obs.elev/R_earth)

    obs.date = date
    sun = ephem.Sun()
    sunset = obs.previous_setting(sun)
    obs.date = sunset
    sunrise = obs.next_rising(sun)

    obs.date = sunset
    obs.horizon = -ephem.degrees('18:00:00.0')
    eve18 = obs.next_setting(sun)
    morn18 = obs.next_rising(sun)

    obs.horizon = -ephem.degrees('12:00:00.0')
    eve12 = obs.next_setting(sun)
    morn12 = obs.next_rising(sun)

    return (sunset, eve12, eve18, morn18, morn12, sunrise)

def plot_measurements(mm, plotfn, gvs, mjds=[], mjdrange=None, allobs=None,
                      markmjds=[], show_plot=True):
    import pylab as plt
    T = db_to_fits(mm)
    print(len(T), 'exposures')
    
    T.mjd_end = T.mjd_obs + T.exptime / 86400.

    #Tall = T
    Tnonobject = T[T.obstype != 'object']
    print(len(Tnonobject), 'exposures are not OBJECTs')
    T = T[T.obstype == 'object']
    print(len(T), 'OBJECT exposures')

    ccmap = dict(g='g', r='r', z='m')

    #bands = 'grz'
    bands = np.unique(T.band)
    #print('Unique bands:', bands)
    
    TT = []
    for band in bands:
        G = T[T.band == band]
        TT.append(G)

    plt.clf()
    plt.subplots_adjust(hspace=0.1, top=0.98, right=0.95, left=0.1,
                        bottom=0.07)

    def limitstyle(band):
        return dict(mec='k', mfc=ccmap[band], ms=8, mew=1)

    # Check for bad things that can happen
    bads = []

    # bad_pixcnt
    I = np.flatnonzero(T.bad_pixcnt)
    for i in I:
        bads.append((i, 'pixcnt'))

    # low nmatched
    I = np.flatnonzero((T.nmatched >= 0) * (T.nmatched < 10))
    for i in I:
        bads.append((i, 'nmatched'))

    if allobs is not None:

        # Duplicate md5sum
        for i in range(len(T)):
            if T.md5sum[i] == '':
                continue
            a = allobs.filter(md5sum=T.md5sum[i]).exclude(
                filename=T.filename[i], extension=T.extension[i])

            if a.count():
                # Only report this one as bad if one of the repeats was earlier
                for ai in a:
                    if ai.mjd_obs < T.mjd_obs[i]:
                        bads.append((i, 'md5sum'))
                        break
                
        # Duplicate readtime
        for i in range(len(T)):
            if T.readtime[i] == 0:
                continue
            a = allobs.filter(readtime=T.readtime[i]).exclude(
                filename=T.filename[i], extension=T.extension[i])
            if a.count():

                # Only report this one as bad if one of the repeats was earlier
                # ... but not too much earlier (within a day)
                for ai in a:
                    if ai.mjd_obs < T.mjd_obs[i] and (ai.mjd_obs + 1) > T.mjd_obs[i]:
                        
                        bads.append((i, 'readtime'))
                        break

    # Group together bad things for the same image.
    bd = {}
    for i,reason in bads:
        if i in bd:
            bd[i] = bd[i] + ', ' + reason
        else:
            bd[i] = reason
    bads = bd.items()

    ilatest = np.argmax(T.mjd_obs)
    latest = T[ilatest]

    SP = 5
    mx = 2.
    plt.subplot(SP,1,1)
    for band,Tb in zip(bands, TT):
        plt.plot(Tb.mjd_obs, Tb.seeing, 'o', color=ccmap[band])

        I = np.flatnonzero(Tb.seeing > mx)
        if len(I):
            plt.plot(Tb.mjd_obs[I], [mx]*len(I), '^', **limitstyle(band))
    yl,yh = plt.ylim()
    plt.axhline(1.3, color='k', alpha=0.5)
    plt.axhline(1.2, color='k', alpha=0.1)
    plt.axhline(1.0, color='k', alpha=0.1)
    plt.axhline(0.8, color='k', alpha=0.1)

    plt.text(latest.mjd_obs, yl+0.01*(yh-yl),
             '%.2f' % latest.seeing, ha='center')

    y = yl + 0.01*(yh-yl)
    plt.plot(np.vstack((T.mjd_obs, T.mjd_end)),
             np.vstack((y, y)), '-', lw=3, alpha=0.5, color=ccmap[band],
             solid_joinstyle='bevel')

    plt.ylim(yl,min(yh,mx))
    plt.ylabel('Seeing (arcsec)')

    ax = plt.axis()
    for i,reason in bads:
        plt.axvline(T.mjd_obs[i], color='r', lw=3, alpha=0.3)
        plt.text(T.mjd_obs[i], ax[3], reason,
                 rotation=90, va='top')
    plt.axis(ax)

    plt.subplot(SP,1,2)
    mx = 20
    for band,Tb in zip(bands, TT):
        nom = get_nominal_cal(Tb.camera[0], band)
        zp0,sky0,kx0 = nom
        plt.plot(Tb.mjd_obs, Tb.sky, 'o', color=ccmap[band])
        plt.axhline(sky0, color=ccmap[band], alpha=0.5)

        I = np.flatnonzero(Tb.sky > mx)
        if len(I):
            plt.plot(Tb.mjd_obs[I], [mx]*len(I), '^', **limitstyle(band))
    yl,yh = plt.ylim()
    yh = min(yh,mx)

    plt.text(latest.mjd_obs, yl+0.01*(yh-yl),
             '%.2f' % latest.sky, ha='center')

    for t in T:
        if t.passnumber > 0:
            plt.text(t.mjd_obs, min(mx, t.sky) - 0.03*(yh-yl),
                     #yl+0.1*(yh-yl),
                     '%i' % t.passnumber, ha='center', va='top')
    
    plt.ylim(yl,yh)
    plt.ylabel('Sky (mag)')

    plt.subplot(SP,1,3)
    mx = 1.2
    mn = 0.5
    for band,Tb in zip(bands, TT):
        plt.plot(Tb.mjd_obs, Tb.transparency, 'o', color=ccmap[band])
        I = np.flatnonzero(Tb.transparency > mx)
        if len(I):
            plt.plot(Tb.mjd_obs[I], [mx]*len(I), '^', **limitstyle(band))
        I = np.flatnonzero(Tb.transparency < mn)
        if len(I):
            plt.plot(Tb.mjd_obs[I], [mn]*len(I), 'v', **limitstyle(band))

    plt.axhline(1.0, color='k', alpha=0.5)
    plt.axhline(0.9, color='k', ls='--', alpha=0.5)
    plt.ylabel('Transparency')
    yl,yh = plt.ylim()
    yl,yh = min(0.89, max(mn, yl)), min(mx, max(yh, 1.01))

    plt.text(latest.mjd_obs, yl+0.01*(yh-yl),
             '%.2f' % latest.transparency, ha='center')

    plt.ylim(yl, yh)
    plt.subplot(SP,1,4)
    mx = 300
    for band,Tb in zip(bands, TT):
        basetime = gvs.base_exptimes[band]
        lo,hi = gvs.floor_exptimes[band], gvs.ceil_exptimes[band]
        exptime = basetime * Tb.expfactor
        clipped = np.clip(exptime, lo, hi)
        if band == 'z':
            clipped = np.minimum(clipped, gvs.t_sat_max)
        Tb.clipped_exptime = clipped
        Tb.depth_factor = Tb.exptime / clipped
        I = np.flatnonzero(exptime > clipped)
        if len(I):
            plt.plot(Tb.mjd_obs[I], exptime[I], 'v', **limitstyle(band))

        I = np.flatnonzero(exptime > mx)
        if len(I):
            plt.plot(Tb.mjd_obs[I], [mx]*len(I), '^', **limitstyle(band))

        I = np.flatnonzero(exptime < clipped)
        if len(I):
            plt.plot(Tb.mjd_obs[I], exptime[I], '^', **limitstyle(band))
        plt.plot(Tb.mjd_obs, clipped, 'o', color=ccmap[band])
        plt.plot(Tb.mjd_obs, Tb.exptime, 'o', mec='k', mfc='none')
    yl,yh = plt.ylim()
    for band,Tb in zip(bands, TT):

        dt = dict(g=-0.5,r=+0.5,z=0)[band]

        basetime = gvs.base_exptimes[band]
        plt.axhline(basetime+dt, color=ccmap[band], alpha=0.2)
        lo,hi = gvs.floor_exptimes[band], gvs.ceil_exptimes[band]
        plt.axhline(lo+dt, color=ccmap[band], ls='--', alpha=0.5)
        plt.axhline(hi+dt, color=ccmap[band], ls='--', alpha=0.5)
        if band == 'z':
            plt.axhline(gvs.t_sat_max, color=ccmap[band], ls='--', alpha=0.5)
    plt.ylim(yl,min(mx, yh))
    plt.ylabel('Target exposure time (s)')

    plt.subplot(SP,1,5)
    mn,mx = 0.6, 1.4
    for band,Tb in zip(bands, TT):
        plt.plot(Tb.mjd_obs, Tb.depth_factor, 'o', color=ccmap[band])
        # lower and upper limits
        I = np.flatnonzero(Tb.depth_factor < mn)
        if len(I):
            plt.plot(Tb.mjd_obs[I], [mn]*len(I), 'v', **limitstyle(band))
        I = np.flatnonzero(Tb.depth_factor > mx)
        if len(I):
            plt.plot(Tb.mjd_obs[I], [mx]*len(I), '^', **limitstyle(band))
    plt.axhline(1.0, color='k', alpha=0.5)
    plt.axhline(0.9, color='k', ls='--', alpha=0.5)
    plt.axhline(1.1, color='k', ls='--', alpha=0.5)

    F = Tnonobject[Tnonobject.obstype == 'focus']
    print(len(F), 'focus frames')
    if len(F):
        plt.plot(F.mjd_obs, 0.9 + np.zeros(len(F)), 'ko')
        for f in F:
            plt.text(f.mjd_obs, 0.88, 'F', ha='center', va='top')
    
    if len(T) > 50:
        ii = [np.argmin(T.expnum + (T.expnum == 0)*1000000),
              np.argmax(T.expnum)]
    else:
        ii = range(len(T))
    for i in ii:
        if T.expnum[i] == 0:
            continue
        plt.text(T.mjd_obs[i], mx, '%i ' % T.expnum[i],
                 rotation=90, va='top', ha='center')
        if len(T) <= 50:
            # Mark focus frames too
            for i in range(len(F)):
                plt.text(F.mjd_obs[i], mx, '%i ' % F.expnum[i],
                         rotation=90, va='top', ha='center')

    plt.ylim(mn,mx)
    plt.ylabel('Depth factor')
    
    plt.xlabel('MJD')
    
    if mjdrange is not None:
        for sp in range(SP):
            plt.subplot(SP,1,sp+1)
            plt.xlim(mjdrange)

    plt.subplot(SP,1,SP)
    xl,xh = plt.xlim()
    if xh - xl < 2./24.:
        # range of less than 2 hours: label every ten minutes
        tx = []
        tt = []
        dl = mjdtodate(xl)
        d = datetime.datetime(dl.year, dl.month, dl.day, dl.hour)
        while True:
            mjd = datetomjd(d)
            if mjd > xh:
                break
            if mjd > xl:
                tx.append(mjd)
                tt.append(d.strftime('%H:%M'))
            d += datetime.timedelta(0, 600)
        plt.xticks(tx, tt)
        plt.xlabel(dl.strftime('Time UTC starting %Y-%m-%d'))

    elif xh - xl < 1:
        # range of less than a day: label the hours
        tx = []
        tt = []
        dl = mjdtodate(xl)
        d = datetime.datetime(dl.year, dl.month, dl.day, dl.hour)
        while True:
            mjd = datetomjd(d)
            if mjd > xh:
                break
            if mjd > xl:
                tx.append(mjd)
                tt.append(d.strftime('%H:%M'))
            d += datetime.timedelta(0, 3600)
        plt.xticks(tx, tt)
        plt.xlabel(dl.strftime('Time UTC starting %Y-%m-%d'))

    else:
        tx,tt = plt.xticks()

    # Set consistent tick marks but no labels on top plots
    tt = ['']*len(tx)
    for sp in range(1, SP):
        plt.subplot(SP,1,sp)
        plt.xticks(tx,tt)
        plt.xlim(xl,xh)

    plt.xlim(xl,xh)

    if len(markmjds):
        for sp in range(1, SP+1):
            plt.subplot(SP,1,sp)
            ax = plt.axis()
            for mjd,c in markmjds:
                plt.axvline(mjd, color=c, alpha=0.5, lw=2)
            plt.axis(ax)

    if show_plot:
        plt.draw()
        plt.show(block=False)
        plt.pause(0.001)
    plt.savefig(plotfn)
    print('Saved', plotfn)
    if show_plot:
        plt.draw()
        plt.show(block=False)
        plt.pause(0.001)
    
def ephemdate_to_mjd(edate):
    # pyephem.Date is days since noon UT on the last day of 1899.
    # MJD is days since midnight UT on 1858/11/17
    # This constant offset in days was computed via:
    #   mjdnow = datetomjd(datetime.datetime.utcnow())
    #   enow = ephem.now()
    #   mjdnow - enow ==> 15019.499915068824
    mjd = float(edate) + 15019.5
    return mjd

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

def set_tile_fields(ccd, hdr, tiles):
    obj = hdr['OBJECT']
    print('Object name', obj)
    ccd.object = obj
    tile = get_tile_from_name(obj, tiles)
    if tile is not None:
        ccd.tileid = tile.tileid
        ccd.passnumber = tile.get('pass')
        ccd.tileebv = tile.ebv_med
    print('Tile id', ccd.tileid, 'pass', ccd.passnumber)

# SFD map isn't picklable, use global instead
gSFD = None

def process_image(fn, ext, gvs, sfd, opt, obs, tiles):
    portion = opt.portion
    db = opt.db
    print('Reading', fn)

    if sfd is None:
        sfd = gSFD

    # Read primary FITS header
    phdr = fitsio.read_header(fn)

    obstype = phdr.get('OBSTYPE','').strip()
    print('obstype:', obstype)
    exptime = phdr.get('EXPTIME')
    expnum = phdr.get('EXPNUM', 0)

    filt = phdr['FILTER']
    filt = filt.strip()
    filt = filt.split()[0]

    airmass = phdr['AIRMASS']
    ra  = hmsstring2ra (phdr['RA'])
    dec = dmsstring2dec(phdr['DEC'])
    
    # Write QA plots to files named by the exposure number
    print('Exposure number:', expnum)

    skip = False
    if obstype in ['zero', 'focus', 'dome flat', '']:
        print('Skipping obstype =', obstype)
        skip = True
    if exptime == 0:
        print('Exposure time EXPTIME in header =', exptime)
        skip = True
    if expnum == '':
        print('No expnum in header')
        skip = True
    if filt == 'solid':
        print('Solid (block) filter.')
        skip = True

    if skip and not db:
        return None

    if db:
        import obsdb
        if ext is None:
            ext = get_default_extension(fn)
        m,created = obsdb.MeasuredCCD.objects.get_or_create(
            filename=fn, extension=ext)
        m.obstype = obstype
        m.camera  = camera_name(phdr)
        m.expnum  = expnum
        m.exptime = exptime
        m.mjd_obs = phdr['MJD-OBS']
        m.airmass = airmass
        m.rabore  = ra
        m.decbore = dec
        m.band = phdr['FILTER'][0]
        m.bad_pixcnt = ('PIXCNT1' in phdr)
        m.readtime = phdr.get('READTIME', 0.)

    if opt.focus and obstype == 'focus' and m.camera == 'mosaic3':
        from mosaic_focus import Mosaic3FocusMeas
        show_plot = opt.show
        if show_plot:
            import pylab as plt
            plt.figure(2, figsize=(8,10))
        if ext is None:
            ext = get_default_extension(fn)
        meas = Mosaic3FocusMeas(fn, ext)
        focusfn = 'focus.png'
        meas.run(ps=None, plotfn=focusfn)
        print('Wrote', focusfn)
        if show_plot:
            plt.draw()
            plt.show(block=False)
            plt.pause(0.001)
            plt.figure(1)
        
    if skip:
        m.save()
        return None
        
    if opt.doplots:
        from astrometry.util.plotutils import PlotSequence
        ps = PlotSequence('qa-%i' % expnum)
        ps.printfn = False
    else:
        ps = None

    # Measure the new image
    kwa = {}
    if ext is not None:
        kwa.update(ext=ext)
    if opt.n_fwhm is not None:
        kwa.update(n_fwhm=opt.n_fwhm)
    M = measure_raw(fn, ps=ps, **kwa)

    if opt.doplots:
        # Gather all the QAplots into a single pdf and clean them up.
        qafile = 'qa-%i.pdf' % expnum
        pnglist = sorted(glob('qa-%i-??.png' % expnum))
        cmd = 'convert {} {}'.format(' '.join(pnglist), qafile)
        print('Writing out {}'.format(qafile))
        #print(cmd)
        os.system(cmd)
        if not opt.keep_plots:
            [os.remove(png) for png in pnglist]

    # (example results for testig)
    #M = {'seeing': 1.4890481099577366, 'airmass': 1.34,
    #'skybright': 18.383479116033314, 'transparency': 0.94488537276869045,
    #'band': 'z', 'zp': 26.442847814941093}

    #print('Measurements:', M)

    trans = M.get('transparency', 0)
    band = M['band']
    # Look up E(B-V) in SFD map
    ebv = sfd.ebv(ra, dec)[0]
    print('E(B-V): %.3f' % ebv)

    if trans > 0:

        gvs.transparency = trans

        fakezp = -99
        expfactor = ExposureFactor(band, airmass, ebv, M['seeing'], fakezp,
                                   M['skybright'], gvs)
        print('Exposure factor:              %6.3f' % expfactor)
        t_exptime = expfactor * gvs.base_exptimes[band]
        print('Target exposure time:         %6.1f' % t_exptime)
        t_exptime = np.clip(t_exptime, gvs.floor_exptimes[band],
                          gvs.ceil_exptimes[band])
        print('Clipped exposure time:        %6.1f' % t_exptime)
    
        if band == 'z' and t_exptime > gvs.t_sat_max:
            t_exptime = gvs.t_sat_max
            print('Reduced exposure time to avoid z-band saturation: %6.1f', t_exptime)

        print

        print('Actual exposure time taken:   %6.1f' % exptime)
    
        print('Depth (exposure time) factor: %6.3f' % (exptime / t_exptime))
        
        # If you were going to re-plan, you would run with these args:
        plandict = dict(seeing=M['seeing'], transparency=trans)
        # Assume the sky is as much brighter than canonical in each band... unlikely
        dsky = M['skybright'] - gvs.sb_dict[M['band']]
        for b in 'grz':
            plandict['sb'+b] = gvs.sb_dict[b] + dsky
        # Note that nightlystrategy.py takes UTC dates.
        start = datenow()
        # Start the strategy 5 minutes from now.
        start += datetime.timedelta(0, 5*60)
        d = start.date()
        plandict['startdate'] = '%04i-%02i-%02i' % (d.year, d.month, d.day)
        t = start.time()
        plandict['starttime'] = t.strftime('%H:%M:%S')
        # Make an hour-long plan
        end = start + datetime.timedelta(0, 3600)
        d = end.date()
        plandict['enddate'] = '%04i-%02i-%02i' % (d.year, d.month, d.day)
        t = end.time()
        plandict['endtime'] = t.strftime('%H:%M:%S')
    
        # Set "--date" to be the UTC date at previous sunset.
        # (nightlystrategy will ask for the next setting of the sun below
        # -18-degrees from that date to define the sn_18).  We could
        # probably also get away with subtracting, like, 12 hours from
        # now()...
        sun = ephem.Sun()
        obs.date = datenow()
        # not the proper horizon, but this doesn't matter -- just need it to
        # be before -18-degree twilight.
        obs.horizon = 0.
        sunset = obs.previous_setting(sun)
        # pyephem's Date.tuple() splits a date into y,m,d,h,m,s
        d = sunset.tuple()
        #print('Date at sunset, UTC:', d)
        year,month,day = d[:3]
        plandict['date'] = '%04i-%02i-%02i' % (year, month, day)
    
        # Decide the pass.
        goodseeing = plandict['seeing'] < 1.3
        photometric = plandict['transparency'] > 0.9
    
        if goodseeing and photometric:
            passnum = 1
        elif goodseeing or photometric:
            passnum = 2
        else:
            passnum = 3
        plandict['pass'] = passnum
    
        plandict['portion'] = portion
        
        print('Replan command:')
        print()
        print('python2.7 nightlystrategy.py --seeg %(seeing).3f --seer %(seeing).3f --seez %(seeing).3f --sbg %(sbg).3f --sbr %(sbr).3f --sbz %(sbz).3f --transparency %(transparency).3f --start-date %(startdate)s --start-time %(starttime)s --end-date %(enddate)s --end-time %(endtime)s --date %(date)s --portion %(portion)f --pass %(pass)i' % plandict) 
        print()
    else:
        plandict = None
        expfactor = 0.

    rtn = (M, plandict, expnum)
    if not db:
        return rtn

    m.racenter  = M['ra_ccd']
    m.deccenter = M['dec_ccd']
    m.ebv  = ebv
    zp = M.get('zp', 0.)
    if zp is None:
        zp = 0.
    m.zeropoint = zp
    m.transparency = trans
    m.seeing = M.get('seeing', 0.)
    m.sky = M['skybright']
    m.expfactor = expfactor
    m.dx = M.get('dx', 0)
    m.dy = M.get('dy', 0)
    m.nmatched = M.get('nmatched',0)

    img = fitsio.read(fn, ext=1)
    cheaphash = np.sum(img)
    # cheaphash becomes an int64.
    m.md5sum = cheaphash

    set_tile_fields(m, phdr, tiles)

    m.save()

    return rtn

def bounce_process_image(X):
    process_image(*X)

def plot_recent(opt, gvs, markmjds=[], **kwargs):
    import obsdb

    if opt.mjdend is None:
        # Now
        mjd_end = mjdnow()
    else:
        mjd_end = opt.mjdend
    
    if opt.mjdstart is None:
        # an hour ago
        mjd_start = mjd_end - 3600. / (24*3600.)
    else:
        mjd_start = opt.mjdstart
        
    # mjd_start <= mjd_obs <= mjd_end
    mm = obsdb.MeasuredCCD.objects.filter(mjd_obs__gte=mjd_start,
                                          mjd_obs__lte=mjd_end)

    if not len(mm):
        print('No measurements in MJD range', mjd_start, mjd_end)
        return

    camera = mm[0].camera
    allobs = obsdb.MeasuredCCD.objects.filter(camera=camera)

    plotfn = opt.plot_filename

    (sunset, eve12, eve18, morn18, morn12, sunrise) = get_twilight(
        camera, ephem.Date(mjdtodate(mjd_end)))

    markmjds.append((ephemdate_to_mjd(eve18),'b'))
    #print('Evening twi18:', eve18, markmjds[-1])
    markmjds.append((ephemdate_to_mjd(morn18),'b'))
    #print('Morning twi18:', morn18, markmjds[-1])
    markmjds.append((ephemdate_to_mjd(eve12),'g'))
    #print('Evening twi12:', eve12, markmjds[-1])
    markmjds.append((ephemdate_to_mjd(morn12),'g'))
    #print('Morning twi12:', morn12, markmjds[-1])
    
    plot_measurements(mm, plotfn, gvs, allobs=allobs,
                      mjdrange=(mjd_start, mjd_end), markmjds=markmjds,
                      **kwargs)

    # from astrometry.util.fits import fits_table
    # tiles = fits_table(opt.tiles)
    # plt.clf()
    # I = (tiles.in_desi == 1) * (tiles.z_done == 0)
    # plt.plot(tiles.ra[I], tiles.dec[I], 'k.', alpha=0.01)
    # I = (tiles.in_desi == 1) * (tiles.z_done > 0)
    # plt.plot(tiles.ra[I], tiles.dec[I], 'k.', alpha=0.5)
    # plt.plot([m.rabore for m in mm], [m.decbore for m in mm], 'm.')
    # plt.xlabel('RA (deg)')
    # plt.ylabel('Dec (deg)')
    # plt.axis([360,0,-20,90])
    # plt.savefig('radec.png')

def skip_existing_files(imgfns, rawext):
    import obsdb
    fns = []
    for fn in imgfns:
        skipext = rawext
        if skipext is None:
            skipext = get_default_extension(fn)
        mm = obsdb.MeasuredCCD.objects.filter(filename=fn, extension=skipext)
        if mm.count():
            print('Found image', fn, 'in database.  Skipping.')
            continue
        fns.append(fn)
    return fns

    
def main(cmdlineargs=None, get_copilot=False):
    global gSFD
    import optparse
    parser = optparse.OptionParser(usage='%prog')

    plotfn_default = 'recent.png'
    
    parser.add_option('--ext', help='Extension to read for computing observing conditions: default "N4" for DECam, "im4" for Mosaic3', default=None)
    parser.add_option('--extnum', type=int, help='Integer extension to read')
    parser.add_option('--rawdata', help='Directory to monitor for new images: default $MOS3_DATA if set, else "rawdata"', default=None)
    parser.add_option('--portion', help='Portion of the night: default %default', type=float, default='1.0')

    parser.add_option('--n-fwhm', default=None, type=int, help='Number of stars on which to measure FWHM')
    
    parser.add_option('--no-db', dest='db', default=True, action='store_false',
                      help='Do not append results to database')

    parser.add_option('--no-focus', dest='focus', default=True,
                      action='store_false', help='Do not analyze focus frames')
    
    parser.add_option('--fits', help='Write database to given FITS table')
    parser.add_option('--plot', action='store_true',
                      help='Plot recent data and quit')
    parser.add_option('--plot-filename', default=None,
                      help='Save plot to given file, default %s' % plotfn_default)

    parser.add_option('--nightplot', '--night', action='store_true',
                      help="Plot tonight's data and quit")

    parser.add_option('--qa-plots', dest='doplots', default=False,
                      action='store_true', help='Create QA plots')

    parser.add_option('--keep-plots', action='store_true',
                      help='Do not remove PNG-format plots (normally merged into PDF)')
    
    parser.add_option('--mjdstart', type=float, default=None,
                      help='MJD (UTC) at which to start plot')

    now = mjdnow()
    parser.add_option('--mjdend', type=float, default=None,
                      help='MJD (UTC) at which to end plot (default: now, which is %.3f)' % now)

    parser.add_option('--skip', action='store_true',
                      help='Skip images that already exist in the database')

    parser.add_option('--threads', type=int, default=None,
                      help='Run multi-threaded when processing list of files on command-line')

    parser.add_option('--fix-db', action='store_true')

    parser.add_option('--tiles', default='obstatus/mosaic-tiles_obstatus.fits',
                      help='Tiles table, default %default')

    parser.add_option('--no-show', dest='show', default=True, action='store_false',
                      help='Do not show plot window, just save it.')

    if cmdlineargs is None:
        opt,args = parser.parse_args()
    else:
        opt,args = parser.parse_args(cmdlineargs)
        
    if not opt.show:
        import matplotlib
        matplotlib.use('Agg')

    imagedir = opt.rawdata
    if imagedir is None:
        imagedir = os.environ.get('MOS3_DATA', 'rawdata')

    rawext = opt.ext
    if opt.extnum is not None:
        rawext = opt.extnum

    from astrometry.util.fits import fits_table
    tiles = fits_table(opt.tiles)

    from django.conf import settings
    import obsdb

    import pylab as plt
    plt.figure(figsize=(8,10))

    markmjds = []

    if opt.nightplot:
        opt.plot = True

        if opt.plot_filename is None:
            opt.plot_filename = 'night.png'

        # Are we at Tololo or Kitt Peak?  Look for latest image.
        o = obsdb.MeasuredCCD.objects.all().order_by('-mjd_obs')
        cam = o[0].camera
        print('Camera:', cam)

        if opt.mjdstart is not None:
            sdate = ephem.Date(mjdtodate(opt.mjdend))
        else:
            sdate = ephem.Date(datenow())
        
        (sunset, eve12, eve18, morn18, morn12, sunrise) = get_twilight(
            cam, sdate)
        if opt.mjdstart is None:
            opt.mjdstart = ephemdate_to_mjd(sunset)
            print('Set mjd start to sunset:', sunset, opt.mjdstart)
        if opt.mjdend is None:
            opt.mjdend = ephemdate_to_mjd(sunrise)
            print('Set mjd end to sunrise', sunrise, opt.mjdend)

        markmjds.append((ephemdate_to_mjd(eve18),'b'))
        print('Evening twi18:', eve18, markmjds[-1])
        markmjds.append((ephemdate_to_mjd(morn18),'b'))
        print('Morning twi18:', morn18, markmjds[-1])
        markmjds.append((ephemdate_to_mjd(eve12),'g'))
        print('Evening twi12:', eve12, markmjds[-1])
        markmjds.append((ephemdate_to_mjd(morn12),'g'))
        print('Morning twi12:', morn12, markmjds[-1])
            
        
    if opt.plot_filename is None:
        opt.plot_filename = plotfn_default

    if opt.fits:
        ccds = obsdb.MeasuredCCD.objects.all()
        print(ccds.count(), 'measured CCDs')
        T = db_to_fits(ccds)
        T.writeto(opt.fits)
        print('Wrote', opt.fits)
        return 0

    if opt.fix_db:

        from astrometry.util.fits import fits_table
        tiles = fits_table('obstatus/mosaic-tiles_obstatus.fits')

        now = mjdnow()
        
        #ccds = obsdb.MeasuredCCD.objects.all()
        #ccds = obsdb.MeasuredCCD.objects.all().filter(mjd_obs__gt=now - 0.25)
        ccds = obsdb.MeasuredCCD.objects.all().filter(mjd_obs__gt=57434)
        
        print(ccds.count(), 'measured CCDs')
        for ccd in ccds:
            try:
                hdr = fitsio.read_header(ccd.filename, ext=0)
                # band = hdr['FILTER']
                # band = band.split()[0]
                # ccd.band = band

                set_tile_fields(ccd, hdr, tiles)

                ccd.save()
                print('Fixed', ccd.filename)
            except:
                import traceback
                traceback.print_exc()

        return 0
            
    # Get nightlystrategy data structures; use fake command-line args.
    # these don't matter at all, since we only use the ExposureFactor() function
    parser,gvs = getParserAndGlobals()
    nsopt,nsargs = parser.parse_args('--date 2015-01-01 --pass 1 --portion 1'.split())
    obs = setupGlobals(nsopt, gvs)

    if opt.plot:
        plot_recent(opt, gvs, markmjds=markmjds, show_plot=False)
        return 0
        
    print('Loading SFD maps...')
    sfd = SFDMap()
    
    if len(args) > 0:
        mp = None
        if opt.threads > 1:
            gSFD = sfd
            from astrometry.util.multiproc import multiproc
            mp = multiproc(opt.threads)

        if opt.skip:
            fns = skip_existing_files(args, rawext)
        else:
            fns = args
            
        if mp is None:
            for fn in fns:
                process_image(fn, rawext, gvs, sfd, opt, obs, tiles)
        else:
            sfd = None
            mp.map(bounce_process_image,
                   [(fn, rawext, gvs, sfd, opt, obs, tiles) for fn in fns])
        plot_recent(opt, gvs, markmjds=markmjds, show_plot=False)
        return 0
    

    copilot = Copilot(imagedir, rawext, opt, gvs, sfd, obs, tiles)

    # for testability
    if get_copilot:
        return copilot

    copilot.run()
    return 0


class Copilot(object):
    def __init__(self, imagedir, rawext,
                 opt, gvs, sfd, obs, tiles):
        self.imagedir = imagedir
        self.rawext = rawext

        # How many times to re-try processing a new image file
        self.maxFail = 10

        self.sleeptime = 5.
        # How often to re-plot
        self.plotTimeout = 60.

        # How long before we mark a line on the plot because we
        # haven't seen an image.
        self.longtime = 300.
        
        # Objects passed through to process_image, plot_recent.
        self.opt = opt
        self.gvs = gvs
        self.sfd = sfd
        self.obs = obs
        self.tiles = tiles
        
        # Set oldimages to the empty set so that we process
        # backlogged images.
        self.oldimages = set()
        backlog = self.get_new_images()
        # (note to self, need explicit backlog because we skip existing
        # for the backlogged files, unlike new ones.)
        self.backlog = skip_existing_files(
            [os.path.join(imagedir, fn) for fn in backlog], rawext)

        # ... then reset oldimages to the current file list.
        self.oldimages = set(os.listdir(imagedir))

        # initialize timers for plot_if_time_elapsed()
        self.lastPlot = self.lastNewImage = datenow()

        # Keep track of how many times we've failed to process a file...
        self.failCounter = Counter()
        
    def get_new_images(self):
        images = set(os.listdir(self.imagedir))
        newimgs = images - self.oldimages
        newimgs = list(newimgs)
        newimgs = [fn for fn in newimgs if
                   fn.endswith('.fits.fz') or fn.endswith('.fits')]
        return newimgs

    def get_newest_image(self):
        newimgs = self.get_new_images()
        if len(newimgs) == 0:
            return None
        # Take the one with the latest timestamp.
        latest = None
        newestimg = None
        for fn in newimgs:
            st = os.stat(os.path.join(self.imagedir, fn))
            t = st.st_mtime
            if latest is None or t > latest:
                newestimg = fn
                latest = t
        return newestimg

    def plot_if_time_elapsed(self):
        now = datenow()
        dtp = (now - self.lastPlot).total_seconds()
        if dtp < self.plotTimeout:
            return
        dt  = (now - self.lastNewImage).total_seconds()
        print('No new images seen for', dt, 'seconds.')
        markmjds = []
        if dt > self.longtime:
            edate = (self.lastNewImage +
                     datetime.timedelta(0, self.longtime))
            markmjds.append((datetomjd(edate),'r'))
        self.timeout_plot(markmjds=markmjds)

    def timeout_plot(self, **kwargs):
        self.plot_recent(**kwargs)
        
    def plot_recent(self, markmjds=[]):
        plot_recent(self.opt, self.gvs, markmjds=markmjds,
                    show_plot=self.opt.show)
        self.lastPlot = datenow()
            
    def process_image(self, path):
        return process_image(path, self.rawext, self.gvs, self.sfd,
                             self.opt, self.obs, self.tiles)
        
    def run_one(self):
        fn = self.get_newest_image()
        if fn is None:
            self.plot_if_time_elapsed()
            if len(self.backlog) == 0:
                return False
            fn = self.backlog.pop()

        if self.failCounter[fn] >= self.maxFail:
            print('Failed to read file: %s, ext: %s, %i times.' %
                  (fn, self.rawext, self.maxFail) + '  Ignoring.')
            self.oldimages.add(fn)
            return False
            
        path = os.path.join(self.imagedir, fn)
        print('Found new file:', path)
        try:
            print('Trying to open image: %s, ext: %s' % 
                  (path, self.rawext))
            fitsio.read(path, ext=self.rawext)
        except:
            print('Failed to open %s: maybe not fully written yet.'
                  % path)
            self.failCounter.update([fn])
            return False

        try:
            self.process_image(path)
            self.oldimages.add(fn)
            self.lastNewImage = datenow()
        except IOError:
            print('Failed to read FITS image: %s, ext %s' %
                  (path, self.rawext))
            import traceback
            traceback.print_exc()
            self.failCounter.update([fn])
            return False

        self.plot_recent()
        return True

    def run(self):
        print('Checking directory for new files:', self.imagedir)
        sleep = False
        while True:
            print
            if sleep:
                time.sleep(self.sleeptime)
            gotone = self.run_one()
            sleep = not gotone

if __name__ == '__main__':
    import obsdb
    obsdb.django_setup()
    main()
