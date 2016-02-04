#! /usr/bin/env python2.7
'''

This script is meant to be run during DECaLS observing.  It waits for
new images to appear, measures their sky brightness, seeing, and
transparency, and advises whether & how to replan.

'''
from __future__ import print_function
import sys
import os
import re

try:
    from collections import OrderedDict
except:
    print('Failed to import OrderedDict.  You are probably using python 2.6.  Please re-run with python2.7')
    sys.exit(-1)

import time
import json
import datetime
from glob import glob
import optparse

import numpy as np

import matplotlib
matplotlib.use('Agg')
import pylab as plt

import fitsio
import ephem

from astrometry.util.plotutils import PlotSequence
from astrometry.util.starutil_numpy import hmsstring2ra, dmsstring2dec, mjdtodate, datetomjd

from nightlystrategy import ExposureFactor, getParserAndGlobals, setupGlobals

from measure_raw import measure_raw, get_nominal_cal

from tractor.sfd import SFDMap


def db_to_fits(mm):
    from astrometry.util.fits import fits_table
    T = fits_table()
    for field in ['filename', 'extension', 'expnum', 'exptime', 'mjd_obs',
                  'airmass', 'racenter', 'deccenter', 'rabore', 'decbore',
                  'band', 'ebv', 'zeropoint', 'transparency', 'seeing',
                  'sky', 'expfactor', 'camera']:
        g = getattr(mm[0], field)
        if isinstance(g, basestring):
            T.set(field, np.array([str(getattr(m, field)) for m in mm]))
        else:
            T.set(field, np.array([getattr(m, field) for m in mm]))
    return T

def plot_measurements(mm, ps, gvs, mjds=[], mjdrange=None):
    T = db_to_fits(mm)

    ccmap = dict(g='g', r='r', z='m')

    #bands = 'grz'
    bands = np.unique(T.band)
    print('Unique bands:', bands)
    
    TT = []
    for band in bands:
        G = T[T.band == band]
        TT.append(G)

    plt.clf()
    plt.subplots_adjust(hspace=0.1, top=0.98, right=0.95, left=0.1)
    
    SP = 5
    plt.subplot(SP,1,1)
    for band,Tb in zip(bands, TT):
        plt.plot(Tb.mjd_obs, Tb.seeing, '.', color=ccmap[band])
    yl,yh = plt.ylim()
    plt.axhline(1.3, color='k', alpha=0.5)
    plt.axhline(1.2, color='k', alpha=0.1)
    plt.axhline(1.0, color='k', alpha=0.1)
    plt.axhline(0.8, color='k', alpha=0.1)
    plt.ylim(yl,yh)
    plt.ylabel('Seeing (arcsec)')
    #yl,yh = plt.ylim()
    #plt.ylim(0, yh)

    plt.subplot(SP,1,2)
    for band,Tb in zip(bands, TT):
        nom = get_nominal_cal(Tb.camera[0], band)
        zp0,sky0,kx0 = nom
        plt.plot(Tb.mjd_obs, Tb.sky, '.', color=ccmap[band])
        plt.axhline(sky0, color=ccmap[band], alpha=0.5)
    plt.ylabel('Sky (mag)')

    plt.subplot(SP,1,3)
    mx = 1.2
    for band,Tb in zip(bands, TT):
        plt.plot(Tb.mjd_obs, Tb.transparency, '.', color=ccmap[band])
        I = np.flatnonzero(Tb.transparency > mx)
        if len(I):
            plt.plot(Tb.mjd_obs[I], [mx]*len(I), '^', mec=ccmap[band],
                     mfc='none')
    plt.axhline(1.0, color='k', alpha=0.5)
    plt.axhline(0.9, color='k', ls='--', alpha=0.5)
    plt.ylabel('Transparency')
    yl,yh = plt.ylim()
    plt.ylim(min(0.89, yl), min(mx, max(yh, 1.01)))
    
    # plt.subplot(SP,1,4)
    # for band,Tb in zip(bands, TT):
    #     plt.plot(Tb.mjd_obs, Tb.expfactor, '.', color=ccmap[band])
    # plt.axhline(1.0, color='k', alpha=0.5)
    # plt.axhline(0.9, color='k', ls='--', alpha=0.5)
    # plt.axhline(1.1, color='k', ls='--', alpha=0.5)
    # plt.ylabel('Target exposure factor')

    plt.subplot(SP,1,4)
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
            plt.plot(Tb.mjd_obs[I], exptime[I], 'v', mec=ccmap[band],
                     mfc='none', alpha=0.5)
        I = np.flatnonzero(exptime < clipped)
        if len(I):
            plt.plot(Tb.mjd_obs[I], exptime[I], '^', mec=ccmap[band],
                     mfc='none', alpha=0.5)
        plt.plot(Tb.mjd_obs, clipped, '.', color=ccmap[band])
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
    plt.ylim(yl,yh)
    plt.ylabel('Target exposure time (s)')

    plt.subplot(SP,1,5)
    mn,mx = 0.6, 1.4
    for band,Tb in zip(bands, TT):
        plt.plot(Tb.mjd_obs, Tb.depth_factor, '.', color=ccmap[band])
        # lower and upper limits
        I = np.flatnonzero(Tb.depth_factor < mn)
        if len(I):
            plt.plot(Tb.mjd_obs[I], [mn]*len(I), 'v', mec=ccmap[band], mfc='none')
        I = np.flatnonzero(Tb.depth_factor > mx)
        if len(I):
            plt.plot(Tb.mjd_obs[I], [mx]*len(I), '^', mec=ccmap[band], mfc='none')
    plt.axhline(1.0, color='k', alpha=0.5)
    plt.axhline(0.9, color='k', ls='--', alpha=0.5)
    plt.axhline(1.1, color='k', ls='--', alpha=0.5)
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

    ps.savefig()


    
def ephemdate_to_mjd(edate):
    # pyephem.Date is days since noon UT on the last day of 1899.
    # MJD is days since midnight UT on 1858/11/17
    # This constant offset in days was computed via:
    #   mjdnow = datetomjd(datetime.datetime.utcnow())
    #   enow = ephem.now()
    #   mjdnow - enow ==> 15019.499915068824
    mjd = float(edate) + 15019.5
    return mjd
    
# SFD map isn't picklable, use global instead
gSFD = None

def process_image(fn, ext, gvs, sfd, opt, obs):
    portion = opt.portion
    db = opt.db
    print('Reading', fn)

    if sfd is None:
        sfd = gSFD

    # Read primary FITS header
    phdr = fitsio.read_header(fn)
    expnum = phdr.get('EXPNUM', 0)

    obstype = phdr['OBSTYPE'].strip()
    print('obstype:', obstype)
    if obstype == 'zero':
        print('Zero:', fn)
        return None, None, expnum
    elif obstype == 'focus':
        print('Focus:', fn)
        return None, None, expnum

    if expnum == '':
        print('No expnum in header')
        # HACK -- Mosaic3 assume filename like mos3.64083.fits
        basefn = os.path.basename(fn)
        print('basefn', basefn)
        m = re.match('mos3?\D*(\d*)\.fits.*', basefn)
        if m is None:
            return None, None, expnum
        expnum = m.group(1)
        print('expnum string guessed from filename:', expnum)
        expnum = int(expnum, 10)
        print('Parsed expnum', expnum)

    filt = phdr['FILTER']
    filt = filt.strip()
    filt = filt.split()[0]
    if filt == 'solid':
        print('Solid (block) filter.')
        return None, None, expnum

    # Write QA plots to files named by the exposure number
    print('Exposure number:', expnum)
    ps = PlotSequence('qa-%i' % expnum)

    #ps.printfn = False

    # Measure the new image
    kwa = {}
    if ext is not None:
        kwa.update(ext=ext)
    M = measure_raw(fn, ps=ps, **kwa)

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

    gvs.transparency = M['transparency']
    band = M['band']

    ra  = hmsstring2ra (phdr['RA'])
    dec = dmsstring2dec(phdr['DEC'])
    airmass = phdr['AIRMASS']
    actual_exptime = phdr['EXPTIME']
    
    # Look up E(B-V) in SFD map
    ebv = sfd.ebv(ra, dec)[0]
    print('E(B-V): %.3f' % ebv)

    fakezp = -99
    expfactor = ExposureFactor(band, airmass, ebv, M['seeing'], fakezp,
                               M['skybright'], gvs)
    print('Exposure factor:              %6.3f' % expfactor)
    exptime = expfactor * gvs.base_exptimes[band]
    print('Target exposure time:         %6.1f' % exptime)
    exptime = np.clip(exptime, gvs.floor_exptimes[band], gvs.ceil_exptimes[band])
    print('Clipped exposure time:        %6.1f' % exptime)
    
    if band == 'z' and exptime > gvs.t_sat_max:
        exptime = gvs.t_sat_max
        print('Reduced exposure time to avoid z-band saturation: %6.1f', exptime)

    print

    print('Actual exposure time taken:   %6.1f' % actual_exptime)

    print('Depth (exposure time) factor: %6.3f' % (actual_exptime / exptime))
    
    # If you were going to re-plan, you would run with these args:
    plandict = dict(seeing=M['seeing'], transparency=M['transparency'])
    # Assume the sky is as much brighter than canonical in each band... unlikely
    dsky = M['skybright'] - gvs.sb_dict[M['band']]
    for b in 'grz':
        plandict['sb'+b] = gvs.sb_dict[b] + dsky
    # Note that nightlystrategy.py takes UTC dates.
    start = datetime.datetime.utcnow()
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
    obs.date = datetime.datetime.utcnow()
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

    rtn = (M, plandict, expnum)
    if not db:
        return rtn

    import obsdb
    m,created = obsdb.MeasuredCCD.objects.get_or_create(
        filename=fn, extension=M['extension'])

    m.camera = M['camera']
    m.expnum = expnum
    m.exptime = actual_exptime
    m.mjd_obs = phdr['MJD-OBS']
    m.airmass = airmass
    m.racenter  = M['ra_ccd']
    m.deccenter = M['dec_ccd']
    m.rabore  = ra
    m.decbore = dec
    m.band = band
    m.ebv  = ebv
    m.zeropoint = M['zp']
    m.transparency = M['transparency']
    m.seeing = M['seeing']
    m.sky = M['skybright']
    m.expfactor = expfactor
    m.dx = M['dx']
    m.dy = M['dy']

    m.save()

    return rtn

def bounce_process_image(X):
    process_image(*X)

def plot_recent(opt, gvs):
    ps = PlotSequence('recent')
    import obsdb

    if opt.mjdend is None:
        # Now
        mjd_end = datetomjd(datetime.datetime.utcnow())
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
    plot_measurements(mm, ps, gvs, #mjds=[(mjd_now,'Now')],
                      mjdrange=(mjd_start, mjd_end))

    
def main():
    global gSFD
    parser = optparse.OptionParser(usage='%prog')
    
    parser.add_option('--ext', help='Extension to read for computing observing conditions: default "N4" for DECam, "im4" for Mosaic3', default=None)
    parser.add_option('--extnum', type=int, help='Integer extension to read')
    parser.add_option('--rawdata', help='Directory to monitor for new images: default %default', default='rawdata')
    parser.add_option('--portion', help='Portion of the night: default %default', type=float, default='1.0')

    parser.add_option('--no-db', dest='db', default=True, action='store_false',
                      help='Do not append results to database')

    parser.add_option('--fits', help='Write database to given FITS table')
    parser.add_option('--plot', action='store_true',
                      help='Plot recent data and quit')

    parser.add_option('--keep-plots', action='store_true',
                      help='Do not remove PNG-format plots (normally merged into PDF)')
    
    parser.add_option('--mjdstart', type=float, default=None,
                      help='MJD (UTC) at which to start plot')

    mjdnow = datetomjd(datetime.datetime.utcnow())
    parser.add_option('--mjdend', type=float, default=None,
                      help='MJD (UTC) at which to end plot (default: now, which is %.3f)' % mjdnow)

    parser.add_option('--skip', action='store_true',
                      help='Skip images that already exist in the database')

    parser.add_option('--threads', type=int, default=None,
                      help='Run multi-threaded when processing list of files on command-line')

    parser.add_option('--fix-db', action='store_true')
    
    opt,args = parser.parse_args()

    imagedir = opt.rawdata
    rawext = opt.ext
    if opt.extnum is not None:
        rawext = opt.extnum

    from django.conf import settings
    import obsdb
    obsdb.django_setup()

    plt.figure(figsize=(10,10))
    
    if opt.fits:
        ccds = obsdb.MeasuredCCD.objects.all()
        print(ccds.count(), 'measured CCDs')
        T = db_to_fits(ccds)
        T.writeto(opt.fits)
        print('Wrote', opt.fits)
        sys.exit(0)

    if opt.fix_db:
        ccds = obsdb.MeasuredCCD.objects.all()
        print(ccds.count(), 'measured CCDs')
        for ccd in ccds:
            try:
                hdr = fitsio.read_header(ccd.filename, ext=0)
                band = hdr['FILTER']
                band = band.split()[0]
                ccd.band = band
                ccd.save()
                print('Fixed', ccd.filename)
            except:
                import traceback
                traceback.print_exc()

        sys.exit(0)
            
    # ps = PlotSequence('recent')
    # mm = obsdb.MeasuredCCD.objects.all()
    # #plot_measurements(mm, ps, mjdrange=(57324, 57324.5))
    # plot_measurements(mm, ps, mjdrange=(57324.1, 57324.15))
    # sys.exit(0)
    
    # Get nightlystrategy data structures; use fake command-line args.
    # these don't matter at all, since we only use the ExposureFactor() function
    parser,gvs = getParserAndGlobals()
    nsopt,nsargs = parser.parse_args('--date 2015-01-01 --pass 1 --portion 1'.split())
    obs = setupGlobals(nsopt, gvs)

    if opt.plot:
        plot_recent(opt, gvs)
        sys.exit(0)
        
    print('Loading SFD maps...')
    sfd = SFDMap()
    
    if len(args) > 0:
        mp = None
        if opt.threads > 1:
            gSFD = sfd
            from astrometry.util.multiproc import multiproc
            mp = multiproc(opt.threads)

        fns = []
        for fn in args:
            if opt.skip:
                mm = obsdb.MeasuredCCD.objects.filter(filename=fn, extension=rawext)
                if mm.count():
                    print('Found image', fn, 'in database.  Skipping.')
                    continue
            fns.append(fn)
            
        if mp is None:
            for fn in fns:
                process_image(fn, rawext, gvs, sfd, opt, obs)
        else:
            sfd = None
            mp.map(bounce_process_image,
                   [(fn, rawext, gvs, sfd, opt, obs) for fn in fns])
        plot_recent(opt, gvs)
        sys.exit(0)
    
    
    lastimages = set(os.listdir(imagedir))
    
    print('Checking directory for new files:', imagedir)
    while True:
        first = True
        # Wait for a new image to appear
        while True:
            print
            if not first:
                time.sleep(5.)
            first = False
    
            images = set(os.listdir(imagedir))
            newimgs = images - lastimages
            newimgs = list(newimgs)
            newimgs = [fn for fn in newimgs if
                       fn.endswith('.fits.fz') or fn.endswith('.fits')]
            #print('Found new images:', newimgs)
            if len(newimgs) == 0:
                continue
    
            # Take the one with the latest timestamp.
            latest = None
            newestimg = None
            for fn in newimgs:
                st = os.stat(os.path.join(imagedir, fn))
                t = st.st_mtime
                if latest is None or t > latest:
                    newestimg = fn
                    latest = t
    
            fn = os.path.join(imagedir, newestimg)
            print('Found new file:', fn)
            try:
                print('Trying to open image:', fn, 'extension:', rawext)
                fitsio.read(fn, ext=rawext)
            except:
                print('Failed to open', fn, '-- maybe not fully written yet.')
                #import traceback
                #traceback.print_exc()
                continue
            
            images = images - set(newimgs)
            images.add(newestimg)
            lastimages = images
            break

        try:
            (M, plandict, expnum) = process_image(
                fn, rawext, gvs, sfd, opt, obs)
        except IOError:
            print('Failed to read FITS image:', fn, 'extension', rawext)
            import traceback
            traceback.print_exc()
            continue

        plot_recent(opt, gvs)
        
if __name__ == '__main__':
    main()
