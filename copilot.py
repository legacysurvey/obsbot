'''

This script is meant to be run during DECaLS observing.  It waits for
new images to appear, measures their sky brightness, seeing, and
transparency, and advises whether & how to replan.

'''
from __future__ import print_function
import sys
import os

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

from measure_raw_decam import measure_raw_decam, nominal_cal

from tractor.sfd import SFDMap

def plot_measurements(mm, ps, mjds=[], mjdrange=None):
    from astrometry.util.fits import fits_table
    T = fits_table()
    for field in ['filename', 'extension', 'expnum', 'exptime', 'mjd_obs',
                  'airmass', 'racenter', 'deccenter', 'rabore', 'decbore',
                  'band', 'ebv', 'zeropoint', 'transparency', 'seeing',
                  'sky', 'expfactor']:
        g = getattr(mm[0], field)
        if isinstance(g, basestring):
            T.set(field, np.array([str(getattr(m, field)) for m in mm]))
        else:
            T.set(field, np.array([getattr(m, field) for m in mm]))
    T.writeto('obsdb/obsdb.fits')

    ccmap = dict(g='g', r='r', z='m')

    #bands = 'grz'
    bands = np.unique(T.band)
    
    TT = []
    for band in bands:
        G = T[T.band == band]
        TT.append(G)

    plt.clf()
    plt.subplots_adjust(hspace=0.05)
    
    SP = 4
    plt.subplot(SP,1,1)
    for band,Tb in zip(bands, TT):
        plt.plot(Tb.mjd_obs, Tb.seeing, 'x', color=ccmap[band])
    plt.axhline(1.3, color='k', alpha=0.5)
    plt.ylabel('Seeing (arcsec)')
    #yl,yh = plt.ylim()
    #plt.ylim(0, yh)

    plt.subplot(SP,1,2)
    for band,Tb in zip(bands, TT):
        plt.plot(Tb.mjd_obs, Tb.sky, 'x', color=ccmap[band])
        plt.axhline(nominal_cal[band][1], color=ccmap[band], alpha=0.5)
    plt.ylabel('Sky (mag)')

    plt.subplot(SP,1,3)
    for band,Tb in zip(bands, TT):
        plt.plot(Tb.mjd_obs, Tb.transparency, 'x', color=ccmap[band])
    plt.axhline(1.0, color='k', alpha=0.5)
    plt.axhline(0.9, color='k', ls='--', alpha=0.5)
    plt.ylabel('Transparency')
    yl,yh = plt.ylim()
    plt.ylim(min(0.89, yl), max(yh, 1.01))
    
    plt.subplot(SP,1,4)
    for band,Tb in zip(bands, TT):
        plt.plot(Tb.mjd_obs, Tb.expfactor, 'x', color=ccmap[band])
    plt.axhline(1.0, color='k', alpha=0.5)
    plt.axhline(0.9, color='k', ls='--', alpha=0.5)
    plt.axhline(1.1, color='k', ls='--', alpha=0.5)
    plt.ylabel('Exposure factor')

    plt.xlabel('MJD')
    
    if mjdrange is not None:
        for sp in range(SP):
            plt.subplot(SP,1,sp+1)
            plt.xlim(mjdrange)

    for sp in range(1, SP):
        plt.subplot(SP,1,sp+1)
        plt.xticks([])

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
    
def read_normal(F, ext):
    img = F[ext].read()
    hdr = F[ext].read_header()
    return img,hdr

def process_image(fn, ext, gvs, sfd, opt, obs, db=True):
    portion = opt.portion

    # Read primary FITS header
    phdr = fitsio.read_header(fn)
    # Write QA plots to files named by the exposure number
    expnum = phdr.get('EXPNUM', 0)
    ps = PlotSequence('qa-%i' % expnum)
    ps.printfn = False
    # Measure the new image
    M = measure_raw_decam(fn, ext=ext, ps=ps)
    #M = measure_raw_decam(fn, ext=ext, ps=ps, read_raw=read_normal)

    # Gather all the QAplots into a single pdf and clean them up.
    qafile = 'qa-%i.pdf' % expnum
    pnglist = sorted(glob('qa-%i-??.png' % expnum))
    cmd = 'convert {} {}'.format(' '.join(pnglist), qafile)
    print('Writing out {}'.format(qafile))
    #print(cmd)
    os.system(cmd)
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
    print('E(B-V):', ebv)

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
    for band in 'grz':
        plandict['sb'+band] = gvs.sb_dict[band] + dsky
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
    print('Date at sunset, UTC:', d)
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
        filename=fn, extension=ext)

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
    
    m.save()

    return rtn

    
if __name__ == '__main__':
    
    parser = optparse.OptionParser(usage='%prog')
    
    parser.add_option('--ext', help='Extension to read for computing observing conditions: default %default', default='N4')
    parser.add_option('--extnum', type=int, help='Integer extension to read')
    parser.add_option('--rawdata', help='Directory to monitor for new images: default %default', default='rawdata')
    parser.add_option('--portion', help='Portion of the night: default %default', type=float, default='1.0')

    opt,args = parser.parse_args()

    imagedir = opt.rawdata
    rawext = opt.ext
    if opt.extnum is not None:
        rawext = opt.extnum

    from django.conf import settings
    import obsdb
    obsdb.django_setup()
    ccds = obsdb.MeasuredCCD.objects.all()
    print(ccds.count(), 'measured CCDs')

    ps = PlotSequence('recent')
    mm = obsdb.MeasuredCCD.objects.all()
    #plot_measurements(mm, ps, mjdrange=(57324, 57324.5))
    plot_measurements(mm, ps, mjdrange=(57324.1, 57324.15))
    sys.exit(0)
    
    # Get nightlystrategy data structures; use fake command-line args.
    # these don't matter at all, since we only use the ExposureFactor() function
    parser,gvs = getParserAndGlobals()
    nsopt,nsargs = parser.parse_args('--date 2015-01-01 --pass 1 --portion 1'.split())
    obs = setupGlobals(nsopt, gvs)

    print('Loading SFD maps...')
    sfd = SFDMap()
    
    if len(args) > 0:
        for fn in args:
            print('Reading', fn)
            process_image(fn, rawext, gvs, sfd, opt, obs)
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
            newimgs = [fn for fn in newimgs if fn.endswith('.fits.fz')]
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

        (M, plandict, expnum) = process_image(
            fn, rawext, gvs, sfd, opt, obs)

        ps = PlotSequence('recent')
        #plot_recent(ps)

        mjd_now = datetomjd(datetime.datetime.utcnow())
        # an hour ago
        mjd_start = mjd_now - 3600. / (24*3600.)

        # mjd_obs >= mjd_start
        mm = obsdb.MeasuredCCD.objects.filter(mjd_obs_gte=mjd_start)

        ps = PlotSequence('recent')
        plot_measurements(mm, ps, mjds=[(mjd_now,'Now')],
                          mjdrange=(mjd_start, mjd_now))
        
        
        
