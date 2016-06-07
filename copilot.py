#! /usr/bin/env python2.7
'''

This script is meant to be run during DECaLS/MzLS observing.  It waits
for new images to appear, measures their sky brightness, seeing, and
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

import numpy as np

import fitsio
import ephem

from astrometry.util.starutil_numpy import hmsstring2ra, dmsstring2dec, mjdtodate, datetomjd

from measure_raw import measure_raw, get_default_extension, camera_name

from obsbot import (exposure_factor, get_tile_from_name, NewFileWatcher,
                    mjdnow, datenow)

from tractor.sfd import SFDMap

def db_to_fits(mm):
    from astrometry.util.fits import fits_table
    T = fits_table()
    for field in ['filename', 'extension', 'expnum', 'exptime', 'mjd_obs',
                  'airmass', 'racenter', 'deccenter', 'rabore', 'decbore',
                  'band', 'ebv', 'zeropoint', 'transparency', 'seeing',
                  'sky', 'expfactor', 'camera', 'dx', 'dy', 'nmatched',
                  'md5sum', 'bad_pixcnt', 'readtime',
                  'obstype',
                  'object', 'tileid', 'passnumber', 'tileebv',
                  'affine_x0', 'affine_y0',
                  'affine_dx', 'affine_dxx', 'affine_dxy',
                  'affine_dy', 'affine_dyx', 'affine_dyy',]:
        g = getattr(mm[0], field)
        if isinstance(g, basestring):
            T.set(field, np.array([str(getattr(m, field)) for m in mm]))
        else:
            T.set(field, np.array([getattr(m, field) for m in mm]))
    return T

class Duck(object):
    pass

def get_twilight(camera, date):
    '''
    Returns an object with attributes:
    
    * sunset
    * "eve10": -10 degree evening
    * "eve12": -12 degree evening
    * "eve18": -18 degree evening
    * "morn18": -18 degree morning
    * "morn12": -12 degree morning
    * "morn10": -10 degree morning
    * sunrise

    for the given camera ("mosaic3" or "decam"), following the night
    whose sunset starts BEFORE the given date.

    date: an ephem.Date object (in UTC)
    '''
    t = Duck()
    ## From nightlystrategy / mosaicstrategy
    R_earth = 6378.1e3 # in meters
    if camera == 'mosaic3':
        obs = ephem.Observer()
        obs.lon = '-111.6003'
        obs.lat = '31.9634'
        obs.elev = 2120.0 # meters
        #print('Assuming KPNO')

    elif camera == 'decam':
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
    t.sunset = obs.previous_setting(sun)
    obs.date = t.sunset
    t.sunrise = obs.next_rising(sun)

    obs.date = t.sunset
    obs.horizon = -ephem.degrees('18:00:00.0')
    t.eve18 = obs.next_setting(sun)
    t.morn18 = obs.next_rising(sun)

    obs.horizon = -ephem.degrees('12:00:00.0')
    t.eve12 = obs.next_setting(sun)
    t.morn12 = obs.next_rising(sun)

    obs.horizon = -ephem.degrees('10:00:00.0')
    t.eve10 = obs.next_setting(sun)
    t.morn10 = obs.next_rising(sun)
    
    return t

def plot_measurements(mm, plotfn, nom, mjds=[], mjdrange=None, allobs=None,
                      markmjds=[], show_plot=True, nightly=False):
    import pylab as plt
    T = db_to_fits(mm)
    T.band = np.core.defchararray.replace(T.band, 'zd', 'z')
    print(len(T), 'exposures')

    T.mjd_end = T.mjd_obs + T.exptime / 86400.

    Tnonobject = T[T.obstype != 'object']
    print(len(Tnonobject), 'exposures are not OBJECTs')
    print('Obs types:', np.unique(T.obstype))
    T = T[T.obstype == 'object']
    print(len(T), 'OBJECT exposures')

    if len(T) == 0:
        return
    
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
        #return dict(mec='k', mfc=ccmap[band], ms=8, mew=1)
        return dict(mec='k', mfc='none', ms=7, mew=1)

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

    bbox = dict(facecolor='white', alpha=0.8, edgecolor='none')
    
    SP = 5
    mx = 2.05
    plt.subplot(SP,1,1)
    for band,Tb in zip(bands, TT):
        I = np.flatnonzero(Tb.seeing > 0)
        if len(I):
            plt.plot(Tb.mjd_obs[I], Tb.seeing[I], 'o', color=ccmap[band])

        I = np.flatnonzero(Tb.seeing > mx)
        if len(I):
            plt.plot(Tb.mjd_obs[I], [mx]*len(I), '^', **limitstyle(band))
    yl,yh = plt.ylim()
    plt.axhline(2.0, color='k', alpha=0.5)
    plt.axhline(1.3, color='k', alpha=0.5)
    plt.axhline(1.2, color='k', alpha=0.1)
    plt.axhline(1.0, color='k', alpha=0.1)
    plt.axhline(0.8, color='k', alpha=0.1)

    if nightly:
        plt.text(latest.mjd_obs, yl+0.03*(yh-yl),
                 'Median: %.2f' % np.median(T.seeing), ha='right', bbox=bbox)
    else:
        plt.text(latest.mjd_obs, yl+0.03*(yh-yl),
                 '%.2f' % latest.seeing, ha='center', bbox=bbox)

    #xl,xh = plt.xlim()
    #plt.text((xl+xh)/2., 
    
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

    T.dsky = np.zeros(len(T), np.float32)
    mx = 2.0
    minsky = -0.15
    nomskies = []
    medskies = []
    for band,Tb in zip(bands, TT):
        sky0 = nom.sky(band)
        T.dsky[T.band == band] = Tb.sky - sky0
        I = np.flatnonzero(Tb.sky > 0)
        if len(I):
            plt.plot(Tb.mjd_obs[I], Tb.sky[I] - sky0, 'o', color=ccmap[band])
            minsky = min(minsky, min(Tb.sky[I] - sky0))
            nomskies.append((band, sky0))
            medskies.append((band, np.median(Tb.sky[I])))
        I = np.flatnonzero((Tb.sky - sky0) > mx)
        if len(I):
            plt.plot(Tb.mjd_obs[I], [mx]*len(I), '^', **limitstyle(band))

    yl,yh = plt.ylim()
    yh = min(yh, mx)
    yl = minsky - 0.03 * (yh-yl)
            
    txt = ', '.join(['%s=%.2f' % (band,sky0) for band,sky0 in nomskies])
    xl,xh = plt.xlim()
    plt.text((xl+xh)/2., min(0., (yl + 0.95*(yh-yl))), txt,
             va='bottom', bbox=bbox)
    
    plt.axhline(0, color='k', alpha=0.5)

    if nightly:
        txt = 'Median: ' + ', '.join(['%s=%.2f' % (band,sky)
                                      for band,sky in medskies])
        plt.text(latest.mjd_obs, 0, txt, ha='right', va='top', bbox=bbox)
    else:
        latest = T[ilatest]
        plt.text(latest.mjd_obs, latest.dsky - 0.05*(yh-yl),
                 '%.2f' % latest.sky, ha='center', va='top',
                 bbox=bbox)
    
    # Plot strings of pass 1,2,3
    I = np.argsort(T.mjd_obs)
    TJ = T[I]
    TJ.cut(TJ.passnumber > 0)
    i = 0
    while i < len(TJ):
        t = TJ[i]
        p0 = t.passnumber
        j = i
        while j < len(TJ) and TJ.passnumber[j] == p0:
            j += 1
        print('Exposures from [%i,%i) have pass %i' % (i, j, p0))
        tend = TJ[j-1]
        
        y = yl + 0.1 * (yh-yl)
        if j > i+1:
            plt.plot([t.mjd_obs, tend.mjd_obs], [y, y], 'b-', lw=2, alpha=0.5)
            # add error bar caps on the endpoints
            for mjd in [t.mjd_obs, tend.mjd_obs]:
                plt.plot([mjd, mjd], [y - 0.03*(yh-yl), y + 0.03*(yh-yl)],
                         'b-', lw=2, alpha=0.5)
        plt.text((t.mjd_obs + tend.mjd_obs)/2., y, '%i' % p0,
                 ha='center', va='top')
        i = j
        
    plt.axhline(-0.25, color='k', alpha=0.25)

    plt.ylim(yl,yh)
    plt.ylabel('Sky - nominal (mag)')

    plt.subplot(SP,1,3)
    mx = 1.2
    mn = 0.5
    for band,Tb in zip(bands, TT):
        I = np.flatnonzero(Tb.transparency > 0)
        if len(I):
            plt.plot(Tb.mjd_obs[I], Tb.transparency[I], 'o', color=ccmap[band])
        I = np.flatnonzero(Tb.transparency > mx)
        if len(I):
            plt.plot(Tb.mjd_obs[I], [mx]*len(I), '^', **limitstyle(band))
        I = np.flatnonzero((Tb.transparency < mn) * (Tb.transparency > 0))
        if len(I):
            plt.plot(Tb.mjd_obs[I], [mn]*len(I), 'v', **limitstyle(band))

    plt.axhline(1.0, color='k', alpha=0.5)
    plt.axhline(0.9, color='k', ls='-', alpha=0.25)
    plt.ylabel('Transparency')
    yl,yh = plt.ylim()
    plt.axhline(0.7, color='k', ls='-', alpha=0.25)
    yl,yh = min(0.89, max(mn, yl)), min(mx, max(yh, 1.01))

    if nightly:
        plt.text(latest.mjd_obs, yl+0.03*(yh-yl),
                'Median: %.2f' % np.median(T.transparency),
                ha='right', bbox=bbox)
    else:
        plt.text(latest.mjd_obs, yl+0.03*(yh-yl),
                 '%.2f' % latest.transparency, ha='center', bbox=bbox)

    plt.ylim(yl, yh)
    plt.subplot(SP,1,4)
    mx = 300
    for band,Tb in zip(bands, TT):
        fid = nom.fiducial_exptime(band)
        basetime = fid.exptime
        lo,hi = fid.exptime_min, fid.exptime_max
        exptime = basetime * Tb.expfactor
        clipped = np.clip(exptime, lo, hi)
        if band == 'z':
            t_sat = nom.saturation_time(band, Tb.sky)
            bad = (Tb.sky == 0)
            clipped = np.minimum(clipped, t_sat + bad*1000000)
        Tb.clipped_exptime = clipped
        #Tb.depth_factor = Tb.exptime / clipped
        Tb.depth_factor = Tb.exptime / exptime
        I = np.flatnonzero(exptime > clipped)
        if len(I):
            plt.plot(Tb.mjd_obs[I], exptime[I], 'v', **limitstyle(band))

        I = np.flatnonzero(exptime > mx)
        if len(I):
            plt.plot(Tb.mjd_obs[I], [mx]*len(I), '^', **limitstyle(band))

        I = np.flatnonzero((exptime < clipped) * (exptime > 0))
        if len(I):
            plt.plot(Tb.mjd_obs[I], exptime[I], '^', **limitstyle(band))

        plt.plot(Tb.mjd_obs, clipped, 'o', mec='k', mfc='none', ms=9)

        I = np.flatnonzero(Tb.exptime > 0)
        if len(I):
            plt.plot(Tb.mjd_obs[I], Tb.exptime[I], 'o', color=ccmap[band])

            if not nightly:
                yl,yh = plt.ylim()
                for i in I:
                    plt.text(Tb.mjd_obs[i], Tb.exptime[i] + 0.04*(yh-yl),
                             '%.2f' % (Tb.depth_factor[i]),
                             rotation=90, ha='center', va='bottom')
                    # '%.0f %%' % (100. * Tb.depth_factor[i]),
                yh = max(yh, max(Tb.exptime[I] + 0.3*(yh-yl)))
                plt.ylim(yl,yh)
            
    yl,yh = plt.ylim()
    for band,Tb in zip(bands, TT):
        dt = dict(g=-0.5,r=+0.5,z=0)[band]

        fid = nom.fiducial_exptime(band)
        basetime = fid.exptime
        lo,hi = fid.exptime_min, fid.exptime_max

        plt.axhline(basetime+dt, color=ccmap[band], alpha=0.2)
        plt.axhline(lo+dt, color=ccmap[band], ls='--', alpha=0.5)
        plt.axhline(hi+dt, color=ccmap[band], ls='--', alpha=0.5)
        if band == 'z':
            I = np.flatnonzero(Tb.sky > 0)
            if len(I):
                plt.plot(Tb.mjd_obs[I], t_sat[I], color=ccmap[band],
                         ls='-', alpha=0.5)

    plt.ylim(yl,min(mx, yh))
    plt.ylabel('Exposure time (s)')

    plt.subplot(SP,1,5)

    I = np.argsort(T.mjd_obs)
    Tx = T[I]
    
    CDs = dict([(ext, nom.cdmatrix(ext)) for ext in np.unique(Tx.extension)])
    CD = np.array([CDs[ext] for ext in Tx.extension])
    dra  = (CD[:,0] * Tx.dx + CD[:,1] * Tx.dy) * 3600.
    ddec = (CD[:,2] * Tx.dx + CD[:,3] * Tx.dy) * 3600.

    refdra,refddec = None,None
    print('Camera', Tx.camera[0])
    if Tx.camera[0].strip() == 'mosaic3':
        # Convert into offsets that match Mosstat ie, offsets in im16,
        # plus magical offset of mosstat-copilot on im16.
        if not np.all(Tx.affine_x0 == 0):
            from camera_mosaic import dradec_to_ref_chip
            refdra,refddec = dradec_to_ref_chip(Tx)

    if refdra is not None:
        # imX plotted lightly
        plt.plot(Tx.mjd_obs, dra,  'bo', alpha=0.2)
        plt.plot(Tx.mjd_obs, ddec, 'go', alpha=0.2)
        plt.plot(Tx.mjd_obs, dra,  'b-', alpha=0.1)
        plt.plot(Tx.mjd_obs, ddec, 'g-', alpha=0.1)

        # Predicted im16 plotted heavy
        I = np.flatnonzero(Tx.affine_x0)
        plt.plot(Tx.mjd_obs[I], refdra[I],  'bo')
        plt.plot(Tx.mjd_obs[I], refddec[I], 'go')
        pr = plt.plot(Tx.mjd_obs[I], refdra[I],  'b-', alpha=0.5)
        pd = plt.plot(Tx.mjd_obs[I], refddec[I], 'g-', alpha=0.5)

        if not nightly:
            mjd = Tx.mjd_obs[I[-1]]
            r,d = refdra[I[-1]], refddec[I[-1]]
            yl,yh = plt.ylim()
            plt.text(mjd, yl+0.03*(yh-yl), '(%.1f, %.1f)' % (r,d),
                     ha='center', bbox=bbox)
        
    else:
        plt.plot(Tx.mjd_obs, dra,  'bo')
        plt.plot(Tx.mjd_obs, ddec, 'go')
        pr = plt.plot(Tx.mjd_obs, dra,  'b-', alpha=0.5)
        pd = plt.plot(Tx.mjd_obs, ddec, 'g-', alpha=0.5)
    #plt.legend((pr[0], pd[0]), ('RA', 'Dec'))
    yl,yh = plt.ylim()

    mx = np.percentile(np.abs(np.append(dra, ddec)), 95)
    mx *= 1.2
    
    plt.axhline(0.0, color='k', alpha=0.5)
    plt.axhline(+10., color='k', alpha=0.2)
    plt.axhline(-10., color='k', alpha=0.2)
    plt.ylim(-mx,mx)
    
    # i = np.argmax(T.mjd_obs)
    # plt.text(T.mjd_obs[i], dra[i], '  RA', ha='left', va='center')
    # plt.text(T.mjd_obs[i], ddec[i], '  Dec', ha='left', va='center')
    
    F = Tnonobject[Tnonobject.obstype == 'focus']
    Z = Tnonobject[Tnonobject.obstype == 'zero']
    print(len(F), 'focus frames,', len(Z), 'zero frames')
    if len(F):
        plt.plot(F.mjd_obs, np.zeros(len(F)), 'ro')
        for f in F:
            plt.text(f.mjd_obs, 0., 'F', ha='center', va='top')
    if len(Z):
        plt.plot(Z.mjd_obs, np.zeros(len(Z)), 'mo')
        for f in Z:
            plt.text(f.mjd_obs, 0., 'Z', ha='center', va='top')

    if len(Tx) > 50:
        ii = [np.argmin(Tx.expnum + (Tx.expnum == 0)*1000000),
              np.argmax(Tx.expnum)]
    else:
        ii = range(len(Tx))

    txty = mx * 0.8
    for i in ii:
        if Tx.expnum[i] == 0:
            continue
        plt.text(Tx.mjd_obs[i], txty, '%i ' % Tx.expnum[i],
                 rotation=90, va='center', ha='center',
                 bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
        if len(Tx) <= 50:
            # Mark focus frames too
            for i in range(len(F)):
                plt.text(F.mjd_obs[i], txty, '%i ' % F.expnum[i],
                         rotation=90, va='top', ha='center')
            for i in range(len(Z)):
                plt.text(Z.mjd_obs[i], txty, '%i ' % Z.expnum[i],
                         rotation=90, va='top', ha='center')

    plt.ylabel('dRA (blu), dDec (grn)')
    
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
            for mark in markmjds:
                if len(mark) == 2:
                    mjd,c = mark
                if len(mark) == 3:
                    # So far, this is only used for MISSING IMAGE text.
                    mjd,c,txt = mark
                    plt.text(mjd, (ax[2]+ax[3])/2, txt, rotation=90,
                             va='center', ha='right')
                plt.axvline(mjd, color=c, alpha=0.5, lw=2)

            plt.axis(ax)

    #Tcount = T[(T.passnumber > 0) * (T.bad_pixcnt == False) *
    #           (T.nmatched > 10)]

    #number of images with depth factor < 0.3 + number of images with seeing > 2"
    
    #for band in np.unique(Tcount.band):
    
    for band,Tb in zip(bands, TT):
        for passnum in [1,2,3]:

            Tcount = Tb[(Tb.passnumber == passnum) * (Tb.bad_pixcnt == False) *
                        (Tb.nmatched > 10)]
            N = len(Tcount)
            print('\nBand %s, pass %i: total of %i tiles' % (band, passnum, N))
            if N > 0:
                depth_thresh = 0.3
                seeing_thresh = 2.0
                shallow = (Tcount.depth_factor < depth_thresh)
                blurry = (Tcount.seeing > seeing_thresh)
                if np.sum(shallow):
                    print('  %i have low depth_factor < %g' % (np.sum(shallow), depth_thresh))
                if np.sum(blurry):
                    print('  %i have large seeing > %g' % (np.sum(blurry), seeing_thresh))
                Ngood = np.sum(np.logical_not(shallow) * np.logical_not(blurry))
                print('Band %s, pass %i: total of %i good tiles' % (band, passnum, Ngood))
    
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

def process_image(fn, ext, nom, sfd, opt, obs, tiles):
    db = opt.db
    print('Reading', fn)

    if sfd is None:
        sfd = gSFD

    # Read primary FITS header
    phdr = fitsio.read_header(fn)

    obstype = phdr.get('OBSTYPE','').strip()
    print('obstype:', obstype)
    exptime = phdr.get('EXPTIME', 0)
    expnum = phdr.get('EXPNUM', 0)

    filt = phdr.get('FILTER', None)
    if filt is not None:
        filt = filt.strip()
        filt = filt.split()[0]
    if filt is None:
        filt = ''

    airmass = phdr.get('AIRMASS', 0.)
    ra  = hmsstring2ra (phdr.get('RA', '0'))
    dec = dmsstring2dec(phdr.get('DEC', '0'))

    # DECam rawdata/DECam_00521494.fits.fz -- "NaN"
    if isinstance(airmass, str):
        airmass = 0.
    
    # Write QA plots to files named by the exposure number
    print('Exposure number:', expnum)

    skip = False
    if obstype in ['zero', 'focus', 'dome flat', '']:
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

    mfocus = False
    camera  = camera_name(phdr)
    if opt.focus and obstype == 'focus' and camera == 'mosaic3':
        skip = False
        mfocus = True

    if skip and not db:
        print('Skipping obstype =', obstype)
        return None

    if db:
        import obsdb
        if ext is None:
            ext = get_default_extension(fn)
        m,created = obsdb.MeasuredCCD.objects.get_or_create(
            filename=fn, extension=ext)

        if skip:
            # Also try searching by expnum and ext.
            m2 = obsdb.MeasuredCCD.objects.filter(
                expnum=expnum, extension=ext)
            if m2.count() > 0:
                print('Expnum and extension already exists in db.')
                return None
        
        m.obstype = obstype
        m.camera  = camera_name(phdr)
        m.expnum  = expnum
        m.exptime = exptime
        m.mjd_obs = phdr.get('MJD-OBS', 0.)
        m.airmass = airmass
        m.rabore  = ra
        m.decbore = dec
        m.band = filt
        m.bad_pixcnt = ('PIXCNT1' in phdr)
        m.readtime = phdr.get('READTIME', 0.)

    if mfocus:
        from mosaic_focus import Mosaic3FocusMeas
        show_plot = opt.show
        if show_plot:
            import pylab as plt
            plt.figure(2, figsize=(8,10))
            plt.subplots_adjust(top=0.95)
        if ext is None:
            ext = get_default_extension(fn)
        meas = Mosaic3FocusMeas(fn, ext, nom)
        focusfn = 'focus.png'
        M = meas.run(ps=None, plotfn=focusfn)
        if not 'focus' in M:
            if M['nmatched'] < 10:
                print('FAILED TO MATCH ENOUGH STARS IN FOCUS FRAME -- please '
                      'try another image extension, eg:')
                print('  python mosaic_focus.py --ext im16 %s' % fn)
        else:
            print('Wrote', focusfn)
            if show_plot:
                plt.draw()
                plt.show(block=False)
                plt.pause(0.001)
                plt.figure(1)
        skip = True
        
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
        from glob import glob
        # Gather all the QAplots into a single pdf and clean them up.
        qafile = 'qa-%i.pdf' % expnum
        pnglist = sorted(glob('qa-%i-??.png' % expnum))
        cmd = 'convert {} {}'.format(' '.join(pnglist), qafile)
        print('Writing out {}'.format(qafile))
        #print(cmd)
        os.system(cmd)
        if not opt.keep_plots:
            [os.remove(png) for png in pnglist]

    if not np.isfinite(M['skybright']):
        print('Negative sky measured:', M['rawsky'], '.  Bad pixcnt:',
              m.bad_pixcnt)
        m.save()
        return None

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

        fid = nom.fiducial_exptime(band)

        expfactor = exposure_factor(fid, nom,
                                    airmass, ebv, M['seeing'], M['skybright'],
                                    trans)
        print('Exposure factor:              %6.3f' % expfactor)
        t_exptime = expfactor * fid.exptime
        print('Target exposure time:         %6.1f' % t_exptime)
        t_exptime = np.clip(t_exptime, fid.exptime_min, fid.exptime_max)
        print('Clipped exposure time:        %6.1f' % t_exptime)
    
        if band == 'z':
            t_sat = nom.saturation_time(band, M['skybright'])
            if t_exptime > t_sat:
                t_exptime = t_sat
                print('Reduced exposure time to avoid z-band saturation: %.1f' % t_exptime)

        print

        print('Actual exposure time taken:   %6.1f' % exptime)
    
        print('Depth (exposure time) factor: %6.3f' % (exptime / t_exptime))
        
        # If you were going to re-plan, you would run with these args:
        plandict = dict(seeing=M['seeing'], transparency=trans)
        # Assume the sky is as much brighter than canonical in each band... unlikely
        dsky = M['skybright'] - nom.sky(M['band'])
        for b in 'grz':
            plandict['sb'+b] = nom.sky(b) + dsky
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
    
        ## ??
        plandict['portion'] = 1.0
        
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

    aff = M.get('affine', None)
    if aff is not None:
        m.affine_x0 = aff[0]
        m.affine_y0 = aff[1]
        m.affine_dx = aff[2]
        m.affine_dxx = aff[3]
        m.affine_dxy = aff[4]
        m.affine_dy = aff[5]
        m.affine_dyx = aff[6]
        m.affine_dyy = aff[7]
    
    img = fitsio.read(fn, ext=1)
    cheaphash = np.sum(img)
    # cheaphash becomes an int64.
    m.md5sum = cheaphash

    set_tile_fields(m, phdr, tiles)

    m.save()

    return rtn

def bounce_process_image(X):
    process_image(*X)

    
def mark_twilight(camera, date):
    twi = get_twilight(camera, date)
    mark = []
    mark.append((ephemdate_to_mjd(twi.eve18), 'b'))
    #print('Evening twi18:', eve18, markmjds[-1])
    mark.append((ephemdate_to_mjd(twi.morn18), 'b'))
    #print('Morning twi18:', morn18, markmjds[-1])
    gb = (0., 0.6, 0.6)
    mark.append((ephemdate_to_mjd(twi.eve12), gb))
    #print('Evening twi12:', eve12, markmjds[-1])
    mark.append((ephemdate_to_mjd(twi.morn12), gb))
    #print('Morning twi12:', morn12, markmjds[-1])
    orange = (1., 0.6, 0)
    mark.append((ephemdate_to_mjd(twi.eve10), 'g'))
    mark.append((ephemdate_to_mjd(twi.morn10),'g'))
    return mark
    
def plot_recent(opt, nom, tiles=None, markmjds=[],
                botplanfn=None, **kwargs):
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
        print('No measurements in MJD range', mjd_start, mjd_end, '=> date range',
              mjdtodate(mjd_start), mjdtodate(mjd_end), 'UTC')
        return

    camera = mm[0].camera
    allobs = obsdb.MeasuredCCD.objects.filter(camera=camera)

    plotfn = opt.plot_filename

    markmjds.extend(mark_twilight(camera, ephem.Date(mjdtodate(mjd_end))))
    
    plot_measurements(mm, plotfn, nom, allobs=allobs,
                      mjdrange=(mjd_start, mjd_end), markmjds=markmjds,
                      **kwargs)

    if botplanfn is None or (not os.path.exists(botplanfn)) or tiles is None:
        return

    import pylab as plt
    from astrometry.util.fits import fits_table
    P = fits_table(botplanfn)
    
    mlast = mm.order_by('mjd_obs').last()

    mrecent = mm.order_by('-mjd_obs')[:10]

    ccmap = dict(g='g', r='r', z='m')
    lp,lt = [],[]
    
    plt.clf()
    I = (tiles.in_desi == 1) * (tiles.z_done == 0)
    plt.plot(tiles.ra[I], tiles.dec[I], 'k.', alpha=0.05)
    I = (tiles.in_desi == 1) * (tiles.z_done > 0)
    plt.plot(tiles.ra[I], tiles.dec[I], 'k.', alpha=0.5)

    plt.plot([m.rabore for m in mm], [m.decbore for m in mm], 'k-',
             lw=2, alpha=0.5)
    pr = plt.scatter([m.rabore for m in mm], [m.decbore for m in mm],
                     color=[ccmap.get(m.band,'k') for m in mm], marker='o',
                     s=20)
    # for k,v in ccmap.items():
    #     mmb = [m for m in mm if m.band == k]
    #     if len(mmb) == 0:
    #         continue
    #     pr = plt.plot([m.rabore for m in mmb], [m.decbore for m in mmb], '.',
    #                   color=v)
    lp.append(pr)
    lt.append('Recent')
        
    P.color = np.array([ccmap.get(f,'k') for f in P.filter])
    I = np.flatnonzero(P.type == '1')
    I = I[:10]
    p1 = plt.scatter(P.ra[I], P.dec[I], c=P.color[I], marker='^', alpha=0.5,
                     s=60)
    plt.plot(P.ra[I], P.dec[I], 'k-', alpha=0.1)
    lp.append(p1)
    lt.append('Upcoming P1')
    
    I = np.flatnonzero(P.type == '2')
    I = I[:10]
    p2 = plt.scatter(P.ra[I], P.dec[I], c=P.color[I], marker='s', alpha=0.5,
                     s=60)
    plt.plot(P.ra[I], P.dec[I], 'k-', alpha=0.1)
    lp.append(p2)
    lt.append('Upcoming P2')

    I = np.flatnonzero(P.type == '3')
    I = I[:10]
    p3 = plt.scatter(P.ra[I], P.dec[I], c=P.color[I], marker='p', alpha=0.5,
                     s=60)
    plt.plot(P.ra[I], P.dec[I], 'k-', alpha=0.1)
    lp.append(p3)
    lt.append('Upcoming P3')

    pl = plt.plot(mlast.rabore, mlast.decbore, 'o',
                  color=ccmap.get(mlast.band,'k'), ms=10)
    lp.append(pl[0])
    lt.append('Last exposure')

    I = np.flatnonzero(P.type == 'P')
    plt.plot(P.ra[I], P.dec[I], 'k-', lw=3, alpha=0.5)
    pplan = plt.scatter(P.ra[I], P.dec[I], c=P.color[I], marker='*',
                        s=100)
    lp.append(pplan)
    lt.append('Planned')
    
    plt.xlabel('RA (deg)')
    plt.ylabel('Dec (deg)')
    #plt.axis([360,0,-20,90])

    plt.figlegend(lp, lt, 'upper right')
    
    ralo = min(P.ra.min(), min([m.rabore for m in mrecent]))
    rahi = max(P.ra.max(), max([m.rabore for m in mrecent]))
    declo = min(P.dec.min(), min([m.decbore for m in mrecent]))
    dechi = max(P.dec.max(), max([m.decbore for m in mrecent]))
    dr = rahi - ralo
    dd = dechi - declo
    
    plt.axis([rahi+0.1*dr, ralo-0.1*dr, declo-0.1*dd, dechi+0.1*dd])

    fn = 'radec.png'
    plt.savefig(fn)
    print('Wrote', fn)
    
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

    # Mosaic or Decam?
    from camera import (nominal_cal, ephem_observer, default_extension,
                        tile_path, camera_name, data_env_var, bot_name)
    nom = nominal_cal
    obs = ephem_observer()
    
    plotfn_default = 'recent.png'
    
    parser.add_option('--ext', default=default_extension,
                      help='Extension to read for computing observing conditions: default %default')
    parser.add_option('--extnum', type=int, help='Integer extension to read')
    parser.add_option('--rawdata', help='Directory to monitor for new images: default $%s if set, else "rawdata"' % data_env_var, default=None)

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
    parser.add_option('--ago', type=int, help='Plot N nights ago; with --night')
    
    parser.add_option('--qa-plots', dest='doplots', default=False,
                      action='store_true', help='Create QA plots')

    parser.add_option('--keep-plots', action='store_true',
                      help='Do not remove PNG-format plots (normally merged into PDF)')
    
    parser.add_option('--mjdstart', type=float, default=None,
                      help='MJD (UTC) at which to start plot')

    now = mjdnow()
    parser.add_option('--mjdend', type=float, default=None,
                      help='MJD (UTC) at which to end plot (default: now, which is %.3f)' % now)

    parser.add_option('--datestart', default=None,
                      help='UTC date at which to start plot; eg "2016/02/03 12:35"')
    parser.add_option('--dateend', default=None,
                      help='UTC date at which to end plot')
    
    parser.add_option('--skip', action='store_true',
                      help='Skip images that already exist in the database')

    parser.add_option('--threads', type=int, default=None,
                      help='Run multi-threaded when processing list of files on command-line')

    parser.add_option('--fix-db', action='store_true')

    parser.add_option('--tiles', default=tile_path,
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
        imagedir = os.environ.get(data_env_var, 'rawdata')

    rawext = opt.ext
    if opt.extnum is not None:
        rawext = opt.extnum
    assert(rawext is not None)
        
    from astrometry.util.fits import fits_table
    tiles = fits_table(opt.tiles)

    from django.conf import settings
    import obsdb

    import pylab as plt
    plt.figure(figsize=(8,10))

    if opt.datestart is not None:
        opt.mjdstart = ephemdate_to_mjd(ephem.Date(opt.datestart))
    if opt.dateend is not None:
        opt.mjdend = ephemdate_to_mjd(ephem.Date(opt.dateend))
    
    markmjds = []

    if opt.nightplot:
        opt.plot = True

        if opt.plot_filename is None:
            opt.plot_filename = 'night.png'

        if opt.mjdstart is not None:
            sdate = ephem.Date(mjdtodate(opt.mjdend))
        else:
            sdate = ephem.Date(datenow())

        if opt.ago:
            sdate -= opt.ago
            
        twi = get_twilight(camera_name, sdate)
        if opt.mjdstart is None:
            opt.mjdstart = ephemdate_to_mjd(
                twi.eve10 - np.abs(twi.eve12-twi.eve10))
            print('Set MJD start to', opt.mjdstart)
        if opt.mjdend is None:
            opt.mjdend = ephemdate_to_mjd(
                twi.morn10 + np.abs(twi.morn10-twi.morn12))
            print('Set MJD end to', opt.mjdend)
        markmjds.extend(mark_twilight(camera_name, sdate))
        
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

    botplanfn = '%s-plan.fits' % bot_name
    
    if opt.plot:
        plot_recent(opt, nom, tiles=tiles, markmjds=markmjds, show_plot=False,
                    nightly=opt.nightplot, botplanfn=botplanfn)
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
                process_image(fn, rawext, nom, sfd, opt, obs, tiles)
        else:
            sfd = None
            mp.map(bounce_process_image,
                   [(fn, rawext, nom, sfd, opt, obs, tiles) for fn in fns])
        plot_recent(opt, nom, tiles=tiles, markmjds=markmjds, show_plot=False,
                    botplanfn=botplanfn)
        return 0
    

    copilot = Copilot(imagedir, rawext, opt, nom, sfd, obs, tiles, botplanfn)

    # for testability
    if get_copilot:
        return copilot

    copilot.run()
    return 0


class Copilot(NewFileWatcher):
    def __init__(self, imagedir, rawext,
                 opt, nom, sfd, obs, tiles, botplanfn):
        self.rawext = rawext
        super(Copilot, self).__init__(imagedir, backlog=True)

        # How long before we mark a line on the plot because we
        # haven't seen an image.
        self.longtime = 300.
        
        # Objects passed through to process_file, plot_recent.
        self.opt = opt
        self.nom = nom
        self.sfd = sfd
        self.obs = obs
        self.tiles = tiles
        self.botplanfn = botplanfn
        
        # Record of exposure times we've seen
        self.exptimes = []
        
    def filter_backlog(self, backlog):
        backlog = self.filter_new_files(backlog)
        return skip_existing_files(backlog, self.rawext)

    def filter_new_files(self, fns):
        return [fn for fn in fns if
                fn.endswith('.fits.fz') or fn.endswith('.fits')]

    def try_open_file(self, path):
        print('Trying to open file: %s, ext: %s' % (path, self.rawext))
        fitsio.read(path, ext=self.rawext)

    def timed_out(self, dt):
        self.plot_recent()

    def process_file(self, path):
        R = process_image(path, self.rawext, self.nom, self.sfd,
                          self.opt, self.obs, self.tiles)
        if R is None:
            return
        (M,p,enum) = R
        exptime = M['exptime']
        self.exptimes.append(exptime)
        # Update the length of time we think we should wait for new exposures
        # based on previous few, plus some margin
        self.longtime = np.max(self.exptimes[-5:]) + 60.
        
    def processed_file(self, path):
        self.lastNewFile = datenow()
        self.plot_recent()

    def plot_recent(self, markmjds=None):
        now = datenow()
        dt  = (now - self.lastNewFile).total_seconds()
        print('No new images seen for', dt, 'seconds.')
        if markmjds is None:
            markmjds = []
        if dt > self.longtime:
            edate = self.lastNewFile + datetime.timedelta(0, self.longtime)
            markmjds.append((datetomjd(edate), 'r', 'MISSING IMAGE!'))
        
        plot_recent(self.opt, self.nom, tiles=self.tiles, markmjds=markmjds,
                    show_plot=self.opt.show, botplanfn=self.botplanfn)

if __name__ == '__main__':
    import obsdb
    from camera import database_filename
    obsdb.django_setup(database_filename=database_filename)
    main()
