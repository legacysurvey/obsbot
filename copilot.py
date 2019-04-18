#! /usr/bin/env python3.6
'''

This script is meant to be run during DECaLS/MzLS observing.  It waits
for new images to appear, measures their sky brightness, seeing, and
transparency, and makes plots of the conditions.

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
    '''Converts the obsdb database entries
    (obsdb/{mosaic3,decam}.sqlite3) into FITS format.'''
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
        if str(g) == g:
            T.set(field, np.array([str(getattr(m, field)) for m in mm]))
        else:
            T.set(field, np.array([getattr(m, field) for m in mm]))
    return T

def recent_gr_seeing(recent=30., exps=None):
    '''
    Computes estimates of seeing in g and r bands based on recent measurements.

    *recent*: how far back from now to look, in minutes

    Returns:
    None if no recent exposures are found.
    (gsee, rsee, G, R) -- seeing estimates for g,r bands; exposures used for g,r estimates.  These may include exposures from the other band!
    '''
    if exps is None:
        exps = get_recent_exposures(recent, bands='gr')
    if exps is None:
        return None

    r_exps = exps[np.flatnonzero(exps.band == 'r')]
    g_exps = exps[np.flatnonzero(exps.band == 'g')]

    r_avg = g_avg = None
    if len(r_exps) >= 5:
        r_avg = np.median(r_exps.seeing)
    if len(g_exps) >= 5:
        g_avg = np.median(g_exps.seeing)

    if r_avg is not None and g_avg is not None:
        see_ratio = r_avg / g_avg
        print('Computed recent r/g seeing ratio', see_ratio)
    else:
        see_ratio = 1.0

    recent_see   = exps.seeing[-5:]
    recent_bands = exps.band  [-5:]

    G = exps[-5:]
    g_see = recent_see.copy()
    g_see[recent_bands == 'r'] /= see_ratio
    g = np.median(g_see)
    if g_avg is not None:
        # g = max(g, g_avg)
        if g_avg > g:
            g = g_avg
            G = g_exps

    R = exps[-5:]
    r_see = recent_see.copy()
    r_see[recent_bands == 'g'] *= see_ratio
    r = np.median(r_see)
    if r_avg is not None:
        if r_avg > r:
            r = r_avg
            R = r_exps

    return g,r,G,R

def recent_gr_sky_color(recent=30., pairs=5.):
    '''
    Estimates g-r sky color based on recent measurements.

    *recent*: how far back from now to look, in minutes
    *pairs*: compare pairs of g,r exposures within this many minutes of each other.
    '''
    exps = get_recent_exposures(recent, bands='gr')
    if exps is None:
        return None,0,0,0
    gexps = exps[np.flatnonzero(exps.band == 'g')]
    rexps = exps[np.flatnonzero(exps.band == 'r')]
    print(len(gexps), 'g-band exposures')
    print(len(rexps), 'r-band exposures')

    # Find pairs of g,r exposures within 5 minutes of each other?
    # Or just difference in medians?
    diffs = []
    for gexp in gexps:
        I = np.flatnonzero(np.abs(gexp.mjd_obs - rexps.mjd_obs) < pairs/(60*24))
        #print('g', gexp.expnum, 'has', len(I), 'r-band exposures w/in',
        #      pairs, 'minutes')
        if len(I):
            diffs.append(gexp.sky - rexps.sky[I])
    if len(diffs) == 0:
        return (None, 0, len(gexps), len(rexps))
    diffs = np.hstack(diffs)
    #print('All differences:', diffs)
    diff = np.median(diffs)
    #print('Median g-r diff:', diff)
    return (diff, len(diffs), len(gexps), len(rexps))

def get_recent_ccds(recent, bands=None):
    import obsdb
    ccds = obsdb.MeasuredCCD.objects.filter(
        mjd_obs__gte=mjdnow() - recent/(60*24))
    if bands is not None:
        ccds = ccds.filter(band__in=[b for b in bands])
    print('Found', len(ccds), 'exposures in copilot db, within', recent,
          'minutes', ('' if bands is None else ('bands ' + str(bands))))
    if len(ccds) == 0:
        return None
    ccds = db_to_fits(ccds)
    return ccds

def get_recent_exposures(recent, bands=None):
    ccds = get_recent_ccds(recent, bands=bands)
    if ccds is None:
        return None
    # Find unique exposure numbers
    expnums,I = np.unique(ccds.expnum, return_index=True)
    print(len(expnums), 'unique expnums')
    # Copy the CCD measurements into the exposure measurements
    exps = ccds[I]
    # drop columns??
    # exps.delete_column('extension')

    # Average some measurements based on all extensions per exposure
    for i,expnum in enumerate(expnums):
        I = np.flatnonzero(ccds.expnum == expnum)
        exps.sky[i] = np.mean(ccds.sky[I])
        exps.seeing[i] = np.mean(ccds.seeing[I])
    # Sort by date
    exps.cut(np.argsort(exps.mjd_obs))
    return exps

    

class Duck(object):
    pass

def get_twilight(obs, date):
    '''
    *obs*: an ephem.Observer object
    *date*: an ephem.Date object (in UTC)

    Returns an object with attributes:
    
    * sunset
    * "eve10": -10 degree evening
    * "eve12": -12 degree evening
    * "eve18": -18 degree evening
    * "morn18": -18 degree morning
    * "morn12": -12 degree morning
    * "morn10": -10 degree morning
    * sunrise
    '''
    t = Duck()

    saved_vals = (obs.date, obs.horizon)
    
    obs.date = date
    sun = ephem.Sun()
    t.sunset = obs.previous_setting(sun)
    obs.date = t.sunset
    t.sunrise = obs.next_rising(sun)

    obs.date = t.sunset
    obs.horizon = -ephem.degrees('18:00:00.0')
    t.eve18 = obs.next_setting(sun)
    t.morn18 = obs.next_rising(sun)

    obs.horizon = -ephem.degrees('15:00:00.0')
    t.eve15 = obs.next_setting(sun)
    t.morn15 = obs.next_rising(sun)

    obs.horizon = -ephem.degrees('12:00:00.0')
    t.eve12 = obs.next_setting(sun)
    t.morn12 = obs.next_rising(sun)

    obs.horizon = -ephem.degrees('10:00:00.0')
    t.eve10 = obs.next_setting(sun)
    t.morn10 = obs.next_rising(sun)

    obs.date, obs.horizon = saved_vals

    return t

def plot_measurements(mm, plotfn, nom, mjds=[], mjdrange=None, allobs=None,
                      markmjds=[], show_plot=True, nightly=False,
                      label_nmatched=True, max_seeing=2.5, target_exptime=True,
                      nominal_sky=False):
    '''
    Plots our measurements of the conditions, as in the recent.png and
    night.png plots.
    '''
    import pylab as plt
    T = db_to_fits(mm)
    print('plot_measurements, nightly', nightly, 'target_exptime', target_exptime)
    print(len(T), 'exposures')
    # Replace filter "zd" -> "z", "rd" -> "r".
    T.band = np.array([dict(zd='z', rd='r').get(b, b) for b in T.band])

    T.mjd_end = T.mjd_obs + T.exptime / 86400.

    T.isobject = np.logical_or(T.obstype == 'object', T.obstype == 'science')
    
    Tnonobject = T[np.logical_not(T.isobject)]
    print(len(Tnonobject), 'exposures are not OBJECTs')
    print('Obs types:', np.unique(T.obstype))
    T = T[T.isobject]
    print(len(T), 'OBJECT exposures')

    if len(T) == 0:
        return
    
    ccmap = dict(g='g', r='r', z='m')
    xcolor = '0.5'

    bands = np.unique(T.band)
    print('Unique bands:', bands)
    
    TT = []
    for band in bands:
        TT.append(T[T.band == band])

    plt.clf()
    plt.subplots_adjust(hspace=0.1, top=0.98, right=0.95, left=0.1,
                        bottom=0.07)

    def limitstyle(band):
        return dict(mec='k', mfc='none', ms=7, mew=1)

    # Check for bad things that can happen
    bads = []

    # bad_pixcnt
    I = np.flatnonzero(T.bad_pixcnt)
    for i in I:
        bads.append((i, 'pixcnt'))

    # low nmatched
    if label_nmatched:
        I = np.flatnonzero((T.nmatched >= 0) * (T.nmatched < 10))
        for i in I:
            print('Exposure', T.expnum[i], ': nmatched', T.nmatched[i])
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
    # which ones will we use to set the scale?
    I = np.flatnonzero((T.seeing > 0) * (T.exptime > 30))
    if len(I):
        mn,mx = T.seeing[I].min(), T.seeing[I].max()
    else:
        mn,mx = 0.7, 2.5
    mx = min(mx, max_seeing)
    yl,yh = mn - 0.15*(mx-mn), mx + 0.05*(mx-mn)
    #print('mn,mx', mn,mx, 'yl,yh', yl,yh)
    
    ## Seeing
    plt.subplot(SP,1,1)
    for band,Tb in zip(bands, TT):
        I = np.flatnonzero((Tb.seeing > 0) * (Tb.exptime > 30))
        if len(I):
            plt.plot(Tb.mjd_obs[I], Tb.seeing[I], 'o',
                     color=ccmap.get(band, xcolor), mec='k')
        I = np.flatnonzero(Tb.seeing > mx)
        if len(I):
            plt.plot(Tb.mjd_obs[I], [mx]*len(I), '^', **limitstyle(band))
    plt.axhline(2.0, color='k', alpha=0.5)
    plt.axhline(1.3, color='k', alpha=0.5)
    plt.axhline(1.2, color='k', alpha=0.1)
    plt.axhline(1.0, color='k', alpha=0.1)
    plt.axhline(0.8, color='k', alpha=0.1)

    if nightly:
        I = np.flatnonzero(T.seeing > 0)
        if len(I):
            plt.text(latest.mjd_obs, yl+0.03*(yh-yl),
                     'Median: %.2f' % np.median(T.seeing[I]), ha='right', bbox=bbox)
    else:
        plt.text(latest.mjd_obs, yl+0.03*(yh-yl),
                 '%.2f' % latest.seeing, ha='center', bbox=bbox)

    y = yl + 0.01*(yh-yl)
    plt.plot(np.vstack((T.mjd_obs, T.mjd_end)),
             np.vstack((y, y)), '-', lw=3, alpha=0.5,
             color=ccmap.get(band, xcolor),
             solid_joinstyle='bevel')

    plt.ylim(yl,yh)
    plt.ylabel('Seeing (arcsec)')

    ax = plt.axis()
    for i,reason in bads:
        plt.axvline(T.mjd_obs[i], color='r', lw=3, alpha=0.3)
        plt.text(T.mjd_obs[i], ax[3], reason,
                 rotation=90, va='top')
    plt.axis(ax)

    ## Sky background
    plt.subplot(SP,1,2)

    T.dsky = np.zeros(len(T), np.float32)
    minsky = -0.15
    nomskies = []
    medskies = []

    keepbands = []
    keepTT = []
    for band,Tb in zip(bands, TT):
        if nominal_sky:
            try:
                sky0 = nom.sky(band)
            except KeyError:
                # unknown filter
                print('Unknown filter for sky:', band)
                continue
        else:
            sky0 = 0.
        T.dsky[T.band == band] = Tb.sky - sky0
        keepbands.append(band)
        keepTT.append(Tb)
    TT = keepTT
    bands = keepbands

    # set the range based on:
    # -- do we need to add a "darker than 15-deg twi" cut?
    I = np.flatnonzero((T.seeing > 0) * (T.exptime > 30))
    if len(I):
        mn,mx = T.dsky[I].min(), T.dsky[I].max()
    else:
        mn,mx = -2, 1
    mn = max(mn, -2.0)
    yl,yh = mn - 0.15*(mx-mn), mx + 0.05*(mx-mn)

    for band,Tb in zip(bands, TT):
        if nominal_sky:
            sky0 = nom.sky(band)
        else:
            sky0 = 0.
        I = np.flatnonzero(Tb.sky > 0)
        if len(I):
            plt.plot(Tb.mjd_obs[I], Tb.sky[I] - sky0, 'o', mec='k',
                     color=ccmap.get(band, xcolor))
            minsky = min(minsky, min(Tb.sky[I] - sky0))
            nomskies.append((band, sky0))
            medskies.append((band, np.median(Tb.sky[I])))
        I = np.flatnonzero((Tb.sky - sky0) > mx)
        if len(I):
            plt.plot(Tb.mjd_obs[I], [mx]*len(I), '^', **limitstyle(band))
        I = np.flatnonzero((Tb.sky - sky0) < mn)
        if len(I):
            plt.plot(Tb.mjd_obs[I], [mn]*len(I), '^', **limitstyle(band))

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
        # print('Exposures from [%i,%i) have pass %i' % (i, j, p0))
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
    if nominal_sky:
        plt.ylabel('Sky - nominal (mag)')
    else:
        plt.ylabel('Sky (mag/sq.arcsec)')

    ## Transparency
    plt.subplot(SP,1,3)
    mx = 1.2
    mn = 0.3
    for band,Tb in zip(bands, TT):
        I = np.flatnonzero(Tb.transparency > 0)
        if len(I):
            plt.plot(Tb.mjd_obs[I], Tb.transparency[I], 'o', mec='k',
                     color=ccmap.get(band, xcolor))
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
        I = np.flatnonzero(T.transparency > 0)
        if len(I):
            plt.text(latest.mjd_obs, yl+0.03*(yh-yl),
                     'Median: %.2f' % np.median(T.transparency[I]),
                ha='right', bbox=bbox)
    else:
        plt.text(latest.mjd_obs, yl+0.03*(yh-yl),
                 '%.2f' % latest.transparency, ha='center', bbox=bbox)
    plt.ylim(yl, yh)

    ## Exposure time plot
    plt.subplot(SP,1,4)
    for band,Tb in zip(bands, TT):
        fid = nom.fiducial_exptime(band)
        if fid is None:
            # Band not 'g','r', or 'z'
            print('Unknown band', band)
            continue
        basetime = fid.exptime
        lo,hi = fid.exptime_min, fid.exptime_max
        # Exposure time we should have taken
        exptime = basetime * Tb.expfactor
        clipped = np.clip(exptime, lo, hi)
        if band == 'z':
            t_sat = nom.saturation_time(band, Tb.sky)
            bad = (Tb.sky == 0)
            clipped = np.minimum(clipped, t_sat + bad*1000000)
        Tb.clipped_exptime = clipped
        #Tb.depth_factor = Tb.exptime / clipped
        Tb.depth_factor = Tb.exptime / exptime

        if target_exptime:
            I = np.flatnonzero((exptime < clipped) * (exptime > 0))
            if len(I):
                plt.plot(Tb.mjd_obs[I], exptime[I], '^', **limitstyle(band))
    
            plt.plot(Tb.mjd_obs, clipped, 'o', mec='k', mfc='none', ms=9)

        # Actual exposure times taken, marked with filled colored circles.
        #I = np.flatnonzero(Tb.exptime > 30)
        I = np.flatnonzero(Tb.exptime > 0)
        if len(I):
            plt.plot(Tb.mjd_obs[I], Tb.exptime[I], 'o', mec='k',
                     color=ccmap.get(band, xcolor))

    yl,yh = plt.ylim()
    for band,Tb in zip(bands, TT):
        fid = nom.fiducial_exptime(band)
        if fid is None:
            continue
        basetime = fid.exptime
        lo,hi = fid.exptime_min, fid.exptime_max

        dt = dict(g=-0.5,r=+0.5).get(band, 0.)

        if target_exptime:
            exptime = basetime * Tb.expfactor
            clipped = np.clip(exptime, lo, hi)
            I = np.flatnonzero(exptime > clipped)
            if len(I):
                plt.plot(Tb.mjd_obs[I], exptime[I], 'v', **limitstyle(band))

        if False:
            I = np.flatnonzero(exptime > mx)
            if len(I):
                plt.plot(Tb.mjd_obs[I], [mx]*len(I), '^', **limitstyle(band))

        if target_exptime:
            plt.axhline(basetime+dt, color=ccmap.get(band, xcolor), alpha=0.2)
            plt.axhline(lo+dt, color=ccmap.get(band, xcolor), ls='--', alpha=0.5)
            plt.axhline(hi+dt, color=ccmap.get(band, xcolor), ls='--', alpha=0.5)
        if band == 'z':
            I = np.flatnonzero(Tb.sky > 0)
            if len(I):
                plt.plot(Tb.mjd_obs[I], t_sat[I], color=ccmap.get(band, xcolor),
                         ls='-', alpha=0.5)

        if (not nightly) and target_exptime:
            I = np.flatnonzero(Tb.exptime > 0)
            if len(I):
                for i in I:
                    plt.text(Tb.mjd_obs[i], Tb.exptime[i] + 0.04*(yh-yl),
                             '%.2f' % (Tb.depth_factor[i]),
                             rotation=90, ha='center', va='bottom')
                yh = max(yh, max(Tb.exptime[I] + 0.3*(yh-yl)))

    if not nightly:
        plt.text(latest.mjd_obs, yl+0.03*(yh-yl),
                 '%i s' % int(latest.exptime), ha='center', bbox=bbox)

    plt.ylim(yl,yh)
    plt.ylabel('Exposure time (s)')

    plt.subplot(SP,1,5)

    I = np.argsort(T.mjd_obs)
    Tx = T[I]
    
    CDs = dict([(ext, nom.cdmatrix(ext)) for ext in np.unique(Tx.extension)])
    CD = np.array([CDs[ext] for ext in Tx.extension])
    dra  = (CD[:,0] * Tx.dx + CD[:,1] * Tx.dy) * 3600.
    ddec = (CD[:,2] * Tx.dx + CD[:,3] * Tx.dy) * 3600.

    mx = np.percentile(np.abs(np.append(dra, ddec)), 95)
    # make the range at least +- 10 arcsec.
    mx = max(mx, 10)
    mx *= 1.2
    yl,yh = -mx,mx

    refdra,refddec = None,None
    #print('Camera', Tx.camera[0])
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

            plt.text(mjd, yl+0.03*(yh-yl), '(%.1f, %.1f)' % (r,d),
                     ha='center', bbox=bbox)
        
    else:
        plt.plot(Tx.mjd_obs, dra,  'bo')
        plt.plot(Tx.mjd_obs, ddec, 'go')
        pr = plt.plot(Tx.mjd_obs, dra,  'b-', alpha=0.5)
        pd = plt.plot(Tx.mjd_obs, ddec, 'g-', alpha=0.5)
    #plt.legend((pr[0], pd[0]), ('RA', 'Dec'))

    plt.axhline(0.0, color='k', alpha=0.5)
    plt.axhline(+10., color='k', alpha=0.2)
    plt.axhline(-10., color='k', alpha=0.2)

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
                 rotation=90, va='center', ha='center', fontsize=10,
                 bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
        if len(Tx) <= 50:
            # Mark focus frames too
            for i in range(len(F)):
                plt.text(F.mjd_obs[i], txty, '%i ' % F.expnum[i],
                         rotation=90, va='top', ha='center')
            for i in range(len(Z)):
                plt.text(Z.mjd_obs[i], txty, '%i ' % Z.expnum[i],
                         rotation=90, va='top', ha='center')

    if refdra is None and not nightly:
        # If refdra is available, only label that.
        plt.text(Tx.mjd_obs[-1], yl+0.03*(yh-yl),
                 '(%.1f, %.1f)' % (dra[-1], ddec[-1]), ha='center', bbox=bbox)

    plt.ylim(yl, yh)
    plt.ylabel('dRA (blu), dDec (grn) (arcsec)')
    
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
                if len(mark) in [3,4]:
                    # For twilight markers: only label the first subplot.
                    if len(mark) == 4:
                        mjd,c,txt,subplots = mark
                        if not sp in subplots:
                            txt = None
                    else:
                        # So far, this is only used for MISSING IMAGE text.
                        mjd,c,txt = mark
                    if txt is not None:
                        plt.text(mjd, (ax[2]+ax[3])/2, txt, rotation=90,
                                 va='center', ha='right')
                plt.axvline(mjd, color=c, alpha=0.5, lw=2)

            plt.axis(ax)

    #Tcount = T[(T.passnumber > 0) * (T.bad_pixcnt == False) *
    #           (T.nmatched > 10)]

    #number of images with depth factor < 0.3 + number of images with seeing > 2"
    
    #for band in np.unique(Tcount.band):
    
    depth_thresh = 0.3
    seeing_thresh = 2.0
    Tbad = []
    for band,Tb in zip(bands, TT):
        for passnum in [1,2,3]:

            Tcount = Tb[(Tb.passnumber == passnum) * (Tb.bad_pixcnt == False) *
                        (Tb.nmatched > 10)]
            N = len(Tcount)
            print('\nBand %s, pass %i: total of %i tiles' % (band, passnum, N))
            if N > 0:
                shallow = (Tcount.depth_factor < depth_thresh)
                blurry = (Tcount.seeing > seeing_thresh)
                if np.sum(shallow):
                    print('  %i have low depth_factor < %g' % (np.sum(shallow), depth_thresh))
                if np.sum(blurry):
                    print('  %i have large seeing > %g' % (np.sum(blurry), seeing_thresh))
                Ngood = np.sum(np.logical_not(shallow) * np.logical_not(blurry))
                print('Band %s, pass %i: total of %i good tiles' % (band, passnum, Ngood))
                Tbad.append(Tcount[np.logical_or(shallow, blurry)])

    if len(Tbad):
        from astrometry.util.fits import merge_tables
        Tbad = merge_tables(Tbad)
        print('Possible bad_expid.txt entries:')
        for t in Tbad:
            bad = []
            if t.depth_factor < depth_thresh:
                bad.append('expfactor=%.2f' % t.depth_factor)
            if t.seeing > seeing_thresh:
                bad.append('seeing=%.2f' % t.seeing)
            print('%i %s MzLS_%i_' % (t.expnum, ','.join(bad), t.tileid))
        print()

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
    obj = hdr.get('OBJECT', '')
    #print('Object name', obj)
    ccd.object = obj
    tile = get_tile_from_name(obj, tiles)
    if tile is not None:
        ccd.tileid = tile.tileid
        ccd.passnumber = tile.get('pass')
        ccd.tileebv = tile.ebv_med
    #print('Tile id', ccd.tileid, 'pass', ccd.passnumber)

# SFD map isn't picklable, use global instead
gSFD = None

def get_expnum(phdr):
    expnum = phdr.get('EXPNUM', 0)
    instrument = phdr.get('INSTRUME')
    if instrument is None:
        return expnum
    instrument = instrument.strip()
    # Bok
    if instrument == '90prime':
        date = phdr.get('DATE-OBS').strip()
        #2017-06-15T10:42:01.301'
        yr = date[2:4]
        month = date[5:7]
        day = date[8:10]
        hr = date[11:13]
        minute = date[14:16]
        sec = date[17:19]
        expnum = int(yr + month + day + hr + minute + sec, 10)
        print('Date', date, '-> faked expnum', expnum)
    # DESI CI
    if instrument == 'DESI':
        expnum = phdr.get('EXPID', 0)
    return expnum

def get_filter(phdr):
    filt = phdr.get('FILTER', None)
    if filt is None:
        instrument = phdr.get('INSTRUME')
        if instrument is None:
            return None
        instrument = instrument.strip()
        # DESI CI
        if instrument == 'DESI':
            filt = 'r'
    if filt is not None:
        filt = filt.strip()
        filt = filt.split()[0]
    return filt

def process_image(fn, ext, nom, sfd, opt, obs, tiles):
    if obs is None:
        from camera import ephem_observer
        obs = ephem_observer()
    db = opt.db
    print('Reading', fn)

    if sfd is None:
        sfd = gSFD

    # Read primary FITS header
    phdr = fitsio.read_header(fn, ext=opt.primext)

    obstype = phdr.get('OBSTYPE','').strip()
    if len(obstype) == 0:
        obstype = phdr.get('FLAVOR','').strip()
    print('obstype:', obstype)
    exptime = phdr.get('EXPTIME', 0)
    # pointing cam
    if exptime == 0:
        exptime = phdr.get('EXPOSURE', 0) / 1000.

    expnum = get_expnum(phdr)

    filt = get_filter(phdr)
    if filt is None:
        filt = ''
    print('filter:', filt)

    airmass = phdr.get('AIRMASS', 0.)
    ra  = phdr.get('RA', '0')
    dec = phdr.get('DEC', '0')
    if not (isinstance(ra, float) and isinstance(dec, float)):
        ra  = hmsstring2ra (ra)
        dec = dmsstring2dec(dec)

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
        #m,created = obsdb.MeasuredCCD.objects.get_or_create(
        #    filename=fn, extension=ext)

        mlist = obsdb.MeasuredCCD.objects.filter(
            filename=fn, extension=ext)
        # Arbitrarily take first object if more than one found
        if mlist.count() > 0:
            m = mlist[0]
        else:
            m = obsdb.MeasuredCCD(filename=fn, extension=ext)
            m.save()

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
            plt.figure(3, figsize=(8,10))
            plt.subplots_adjust(top=0.95)
        if ext is None:
            ext = get_default_extension(fn)
        meas = Mosaic3FocusMeas(fn, ext, nom)
        focusfn = 'focus.png'
        M = meas.run(ps=None, plotfn=focusfn)
        if not 'focus' in M:
            if ('nmatched' not in M) or (M['nmatched'] < 10):
                print('FAILED TO MATCH ENOUGH STARS IN FOCUS FRAME -- please '
                      'try another image extension, eg:')
                print('  python mosaic_focus.py --ext im16 %s' % fn)
        else:
            print('Wrote', focusfn)
            if show_plot:
                plt.draw()
                plt.show(block=False)
                plt.pause(0.001)
        skip = True
        
    if skip:
        if db:
            m.save()
        return None

    if opt.doplots:
        from astrometry.util.plotutils import PlotSequence
        ps = PlotSequence('qa-%i' % expnum)
        ps.printfn = False
    else:
        ps = None

    # Measure the new image
    kwa = dict(verbose=opt.verbose, primext=opt.primext)
    if ext is not None:
        kwa.update(ext=ext)
    if opt.n_fwhm is not None:
        kwa.update(n_fwhm=opt.n_fwhm)
    if opt.maxshift is not None:
        kwa.update(measargs=dict(maxshift=opt.maxshift))
    M = measure_raw(fn, ps=ps, **kwa)
    #print(M)
    
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

    skybright = M['skybright']
    if skybright is None:
        skybright = 0.
    if not np.isfinite(skybright):
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

    dx = M.get('dx',0.)
    dy = M.get('dy',0.)
    CD = nom.cdmatrix(ext)
    dra  = (CD[0] * dx + CD[1] * dy) * 3600.
    ddec = (CD[2] * dx + CD[3] * dy) * 3600.

    print('E(B-V):          %.3f' % ebv)
    print('Airmass:         %.3f' % airmass)
    print('Sky brightness: %.3f' % skybright)
    zp = M.get('zp',0.)
    if zp and zp > 0:
        print('Zeropoint:      %.3f' % M.get('zp', 0.))
        print('Transparency:    %.3f' % trans)
    print('Seeing:          %.2f' % M.get('seeing', 0.))
    print('Astrometric offset: (%.2f, %.2f) arcsec' % (dra, ddec))

    if trans > 0:

        fid = nom.fiducial_exptime(band)

        expfactor = exposure_factor(fid, nom,
                                    airmass, ebv, M['seeing'], skybright,
                                    trans)
        #print('Exposure factor:              %6.3f' % expfactor)
        t_exptime = expfactor * fid.exptime
        #print('Target exposure time:         %6.1f' % t_exptime)
        t_unclipped = t_exptime
        t_exptime = np.clip(t_exptime, fid.exptime_min, fid.exptime_max)
        #print('Clipped exposure time:        %6.1f' % t_exptime)

        if band == 'z':
            t_sat = nom.saturation_time(band, skybright)
            if t_exptime > t_sat:
                t_exptime = t_sat
                print('Reduced exposure time to avoid z-band saturation: %.1f' % t_exptime)

        #print
        #print('Actual exposure time taken:   %6.1f' % exptime)
        #print('Depth (exposure time) factor: %6.3f' % (exptime / t_exptime))
        #print('Depth factor (on un-clipped exposure time): %6.3f' % (exptime / t_unclipped))
    else:
        expfactor = 0.

    M.update(expnum=expnum)
    if not db:
        return M

    m.racenter  = M['ra_ccd']
    m.deccenter = M['dec_ccd']
    m.ebv  = ebv
    zp = M.get('zp', 0.)
    if zp is None:
        zp = 0.
    m.zeropoint = zp
    m.transparency = trans
    m.seeing = M.get('seeing', 0.)
    m.sky = skybright
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

    #print('Meas:', M)
    
    img = fitsio.read(fn, ext=M['extension'])
    cheaphash = np.sum(img)
    # cheaphash becomes an int64.
    m.md5sum = cheaphash

    set_tile_fields(m, phdr, tiles)

    m.save()
    return M

def bounce_process_image(X):
    try:
        process_image(*X)
    except:
        print('Failed to process image:', X)
        import traceback
        traceback.print_exc()
    
def mark_twilight(obs, date):
    twi = get_twilight(obs, date)
    mark = []
    mark.append((ephemdate_to_mjd(twi.eve18), 'b', '18', [1]))
    #print('Evening twi18:', eve18, markmjds[-1])
    mark.append((ephemdate_to_mjd(twi.morn18), 'b', '18', [1]))
    #print('Morning twi18:', morn18, markmjds[-1])
    gb = (0., 0.6, 0.6)
    mark.append((ephemdate_to_mjd(twi.eve15), gb, '15', [1]))
    mark.append((ephemdate_to_mjd(twi.morn15), gb, '15', [1]))

    mark.append((ephemdate_to_mjd(twi.eve12), gb, '12', [1]))
    mark.append((ephemdate_to_mjd(twi.morn12), gb, '12', [1]))
    #orange = (1., 0.6, 0)
    mark.append((ephemdate_to_mjd(twi.eve10), 'g', '10', [1]))
    mark.append((ephemdate_to_mjd(twi.morn10),'g', '10', [1]))
    return mark
    
def plot_recent(opt, obs, nom, tiles=None, markmjds=[],
                botplanfn=None, nightly=False, **kwargs):
    ''' Creates the recent.png plot of recent measurements of the conditions.
    '''
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

    if (botplanfn is not None and os.path.exists(botplanfn) and
        tiles is not None):
        import pylab as plt
        plt.figure(2)
        radec_plot(botplanfn, mm, tiles, nightly, mjd_start)
    elif nightly:
        radec_plot(None, mm, tiles, nightly, mjd_start)

    camera = mm[0].camera
    allobs = obsdb.MeasuredCCD.objects.filter(camera=camera)

    plotfn = opt.plot_filename

    markmjds.extend(mark_twilight(obs, ephem.Date(mjdtodate(mjd_end))))

    import pylab as plt
    plt.figure(1)
    plot_measurements(mm, plotfn, nom, allobs=allobs,
                      mjdrange=(mjd_start, mjd_end), markmjds=markmjds,
                      nightly=nightly, **kwargs)


def radec_plot(botplanfn, mm, tiles, nightly, mjdstart):
    import pylab as plt
    from astrometry.util.fits import fits_table
    if botplanfn is None:
        P = None
    else:
        P = fits_table(botplanfn)
        if len(P) == 0:
            P = None
        
    msorted = list(mm.order_by('mjd_obs'))
    mlast   = msorted[-1]
    mrecent = msorted[-10:]

    ccmap = dict(g='g', r='r', z='m', zd='m', rd='r')
    lp,lt = [],[]

    plt.clf()

    if tiles is not None:
        I = (tiles.in_desi == 1) * (tiles.z_done == 0)
        plt.plot(tiles.ra[I], tiles.dec[I], 'k.', alpha=0.05)
        I = (tiles.in_desi == 1) * (tiles.z_done > 0)
        plt.plot(tiles.ra[I], tiles.dec[I], 'k.', alpha=0.25)

    plt.plot([m.rabore for m in msorted], [m.decbore for m in msorted], 'k-',
             lw=2, alpha=0.1)
    pr = plt.scatter([m.rabore for m in msorted], [m.decbore for m in msorted],
                     color=[ccmap.get(m.band[:1],'k') for m in msorted],
                     marker='o', s=20)
    lp.append(pr)
    lt.append('Recent')

    rd = []
    if nightly:
        rd.append((np.array([m.rabore for m in msorted]),
                   np.array([m.decbore for m in msorted])))

    if not nightly and P is not None:
        # Plot the planned exposures per pass.
        P.color = np.array([ccmap.get(f[:1],'k') for f in P.filter])
        I = np.flatnonzero(P.type == '1')
        I = I[:10]
        p1 = plt.scatter(P.ra[I], P.dec[I], c=P.color[I], marker='^', alpha=0.5,
                         s=60)
        rd.append((P.ra[I], P.dec[I]))
        plt.plot(P.ra[I], P.dec[I], 'k-', alpha=0.1)
        lp.append(p1)
        lt.append('Upcoming P1')
        
        I = np.flatnonzero(P.type == '2')
        I = I[:10]
        p2 = plt.scatter(P.ra[I], P.dec[I], c=P.color[I], marker='s', alpha=0.5,
                         s=60)
        rd.append((P.ra[I], P.dec[I]))
        plt.plot(P.ra[I], P.dec[I], 'k-', alpha=0.1)
        lp.append(p2)
        lt.append('Upcoming P2')
    
        I = np.flatnonzero(P.type == '3')
        I = I[:10]
        p3 = plt.scatter(P.ra[I], P.dec[I], c=P.color[I], marker='p', alpha=0.5,
                         s=60)
        rd.append((P.ra[I], P.dec[I]))
        plt.plot(P.ra[I], P.dec[I], 'k-', alpha=0.1)
        lp.append(p3)
        lt.append('Upcoming P3')

    pl = plt.plot(mlast.rabore, mlast.decbore, 'o',
                  color=ccmap.get(mlast.band,'k'), ms=10)
    rd.append(([mlast.rabore], [mlast.decbore]))
    lp.append(pl[0])
    lt.append('Last exposure')

    if not nightly and P is not None:
        I = np.flatnonzero(P.type == 'P')
        # Bold the first few ones
        II = I[:6]
        plt.plot(P.ra[II], P.dec[II], 'k-', lw=3, alpha=0.5)
        II = I[:5]
        pplan = plt.scatter(P.ra[II], P.dec[II], c=P.color[II], marker='*',
                            s=100, alpha=1.0)
        # Faint the rest
        II = I[5:]
        plt.plot(P.ra[II], P.dec[II], 'k-', lw=3, alpha=0.2)
        pplan = plt.scatter(P.ra[II], P.dec[II], c=P.color[II], marker='*',
                            s=100, alpha=0.4)

        rd.append((P.ra[I], P.dec[I]))
        lp.append(pplan)
        lt.append('Planned')

        # Bold line from "most recent" to "first planned"
        if len(I) > 0:
            p = P[I[0]]
            plt.plot([mlast.rabore, p.ra], [mlast.decbore, p.dec], '-', lw=3, alpha=0.5, color=ccmap.get(p.filter,'k'))
    
    plt.xlabel('RA (deg)')
    plt.ylabel('Dec (deg)')

    plt.figlegend(lp, lt, 'upper right', prop={'size': 10}, ncol=2, numpoints=1,
                  bbox_to_anchor=(0.9, 0.97), framealpha=0.5, frameon=True)

    rr = np.hstack([r for r,d in rd])
    dd = np.hstack([d for r,d in rd])
    
    ralo  = rr.min()
    rahi  = rr.max()
    declo = dd.min()
    dechi = dd.max()

    dr = rahi - ralo
    dd = dechi - declo
    
    plt.axis([rahi+0.1*dr, ralo-0.1*dr, declo-0.1*dd, dechi+0.1*dd])

    if nightly:
        ## HACK -- subtract 0.5 for UTC to local calendar date at start of night.
        date = mjdtodate(mjdstart - 0.5)
        plt.title('%i-%02i-%02i' % (date.year, date.month, date.day))
    
    fn = 'radec.png'
    plt.savefig(fn)
    print('Wrote', fn)
    
def skip_existing_files(imgfns, rawext, by_expnum=False, primext=0):
    import obsdb
    fns = []
    for fn in imgfns:
        skipext = rawext
        if skipext is None:
            skipext = get_default_extension(fn)
        if by_expnum:
            expnum = 0
            # check if we've seen this file before (in the database)
            m = obsdb.MeasuredCCD.objects.filter(filename = fn)
            if len(m):
                # find unique expnums, take last in case of zeros
                expnum = np.unique([mi.expnum for mi in m])[-1]
                print('file', fn, 'found in database -- expnum', expnum)
            if expnum == 0:
                from read_header import read_primary_header
                hdr = fitsio.read_header(fn, ext=primext)
                expnum = get_expnum(hdr)
                #hdr = read_primary_header(fn)
                #expnum = get_expnum(hdr)
                print('file', fn, '-> expnum', expnum)
            mm = obsdb.MeasuredCCD.objects.filter(expnum=expnum, extension=skipext)
        else:
            mm = obsdb.MeasuredCCD.objects.filter(filename=fn, extension=skipext)
        if mm.count():
            print('Found image', fn, 'in database.  Skipping.')
            continue
        fns.append(fn)
    return fns

def coverage_plots(opt, camera_name, nice_camera_name):
    import pylab as plt
    import matplotlib as mpl
    
    plt.figure(figsize=(12,4))
    plt.subplots_adjust(left=0.05, right=0.99, bottom=0.05, top=0.98)
    
    ccds = obsdb.MeasuredCCD.objects.all()
    print(ccds.count(), 'measured CCDs')
    T = db_to_fits(ccds)
    print('Camera(s):', np.unique(T.camera))
    print('Camera name for this installation:', camera_name)
    T.cut(T.camera == camera_name)
    print('Cut to', len(T), 'with camera', camera_name)
    T.cut(T.tileid > 0)
    print(len(T), 'with known tile id')
    T.cut(T.expnum > 0)
    print(len(T), 'with known exposure number')
    T.cut(T.exptime > 1)
    print(len(T), 'with exposure time > 1')
    filters = np.unique(T.band)
    print('Filters:', filters)
    # Cut to just first letter of filter name...
    T.band = np.array([f[0] for f in T.band])
    filters = np.unique(T.band)
    print('Converted to:', filters)

    ras = np.linspace(0, 360, 360)
    decs = np.linspace(-90, 90, 180)
    if camera_name == 'mosaic3':
        decs = np.linspace(30, 90, 60)
    elif camera_name == 'decam':
        decs = np.linspace(-20, 35, 30)
        ras = np.linspace(0, 360, 180)

    Igood = np.flatnonzero(T.exptime >= 30.)
    print('RA,Dec range of "good" exposures:')
    print('RA',  T.rabore [Igood].min(),  T.rabore[Igood].max())
    print('Dec', T.decbore[Igood].min(), T.decbore[Igood].max())

    sets = [(T, 'all', 'Total')]

    Ti = T[(T.mjd_obs >= opt.mjdstart) * (T.mjd_obs <= opt.mjdend)]
    d = mjdtodate(opt.mjdstart - 0.5)
    tt = '%04i-%02i-%02i' % (d.year, d.month, d.day)
    nn = '%04i%02i%02i' % (d.year, d.month, d.day)
    print('Found', len(Ti), 'for night of', tt)    
    sets.append((Ti, nn, tt))

    for T,shortname,name in sets:
        for band in filters:
            Tb = T[T.band == band]
            print(len(Tb), 'in', band, 'band')
            passes = np.unique(Tb.passnumber)
            print('Passes:', passes)
            for p in passes:
                Tp = Tb[Tb.passnumber == p]
                print(len(Tp), 'in band', band, 'and pass', p)
    
                print('RA', Tp.rabore.min(), Tp.rabore.max())
                print('Dec', Tp.decbore.min(), Tp.decbore.max())
    
                # Select one entry per exposure number.
                expnums,I = np.unique(Tp.expnum, return_index=True)
                print(len(expnums), 'unique exposure numbers')
                Tp.cut(I)
    
                plt.clf()
                plt.hist2d(Tp.rabore, Tp.decbore, bins=(ras, decs),
                            weights=Tp.exptime,
                            norm=mpl.colors.PowerNorm(0.5))
                plt.colorbar()
                plt.xticks([0,60,120,180,240,300,360])
                plt.xlim(360, 0)
                plt.xlabel('RA (deg)')
                plt.ylabel('Dec (deg)')
                plt.axis('scaled')
                fn = 'cov-%s-%s-p%i.png' % (shortname, band, p)
                plt.title('Exposure time: %s, %s, band %s, pass %i' %
                          (name, nice_camera_name, band, p))
                plt.savefig(fn)
                print('Wrote', fn)
    
def main(cmdlineargs=None, get_copilot=False):
    global gSFD
    import optparse
    parser = optparse.OptionParser(usage='%prog')

    # Mosaic or Decam?
    from camera import (nominal_cal, ephem_observer, default_extension,
                        tile_path, camera_name, data_env_var, bot_name,
                        nice_camera_name)
    try:
        from camera import default_primary_extension
    except:
        default_primary_extension = 0

    try:
        from camera import copilot_plot_args
    except:
        copilot_plot_args = {}

    nom = nominal_cal
    obs = ephem_observer()
    
    plotfn_default = 'recent.png'

    parser.add_option('--primext', default=default_primary_extension, type=int,
                      help='Extension to read for "primary" header')
    parser.add_option('--ext', default=default_extension,
                      help='Extension to read for computing observing conditions: default %default')
    parser.add_option('--extnum', type=int, help='Integer extension to read')
    parser.add_option('--rawdata', help='Directory to monitor for new images: default $%s if set, else "rawdata"' % data_env_var, default=None)

    parser.add_option('--n-fwhm', default=None, type=int, help='Number of stars on which to measure FWHM')
    parser.add_option('--maxshift', default=None, type=int, help='Maximum search radius for astrometric offset vs Pan-STARRS, in arcsec, default 240.')
    
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
                      help="Plot tonight's data and quit", default=False)
    parser.add_option('--ago', type=int, help='Plot N nights ago; with --night')
    
    parser.add_option('--qa-plots', dest='doplots', default=False,
                      action='store_true', help='Create QA plots')

    parser.add_option('--cov', dest='covplots', default=False,
                      action='store_true', help='Create coverage plots')

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
                      help='Skip image filenames that already exist in the database')
    parser.add_option('--skip-expnum', action='store_true',
                      help='Skip images (by exposure number) that already exist in the database')

    parser.add_option('--threads', type=int, default=None,
                      help='Run multi-threaded when processing list of files on command-line')

    parser.add_option('--fix-db', action='store_true')

    parser.add_option('--tiles', default=tile_path,
                      help='Tiles table, default %default')

    parser.add_option('--no-show', dest='show', default=True, action='store_false',
                      help='Do not show plot window, just save it.')

    parser.add_option('--verbose', default=False, action='store_true',
                      help='Turn on (even more) verbose logging')

    if cmdlineargs is None:
        opt,args = parser.parse_args()
    else:
        opt,args = parser.parse_args(cmdlineargs)

    try:
        import pylab
    except:
        opt.show = False
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
    if opt.tiles is None:
        tiles = None
    else:
        tiles = fits_table(opt.tiles)

    from django.conf import settings
    import obsdb

    import pylab as plt

    # Mosaic focus plot
    # figure 3
    # RA,Dec plot
    plt.figure(num=2, figsize=(10,6))
    # Copilot plot
    plt.figure(num=1, figsize=(8,10))

    if opt.datestart is not None:
        opt.mjdstart = ephemdate_to_mjd(ephem.Date(opt.datestart))
    if opt.dateend is not None:
        opt.mjdend = ephemdate_to_mjd(ephem.Date(opt.dateend))
    
    markmjds = []

    if opt.nightplot or opt.covplots:
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
            
        twi = get_twilight(obs, sdate)
        if opt.mjdstart is None:
            opt.mjdstart = ephemdate_to_mjd(
                twi.eve10 - np.abs(twi.eve12-twi.eve10))
            print('Set MJD start to', opt.mjdstart)
        if opt.mjdend is None:
            opt.mjdend = ephemdate_to_mjd(
                twi.morn10 + np.abs(twi.morn10-twi.morn12))
            print('Set MJD end to', opt.mjdend)
        markmjds.extend(mark_twilight(obs, sdate))

    if opt.covplots:
        coverage_plots(opt, camera_name, nice_camera_name)
        return 0
        
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
        
        ccds = obsdb.MeasuredCCD.objects.all()
        #ccds = obsdb.MeasuredCCD.objects.all().filter(mjd_obs__gt=now - 0.25)
        #ccds = obsdb.MeasuredCCD.objects.all().filter(mjd_obs__gt=57434)
        
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
        plot_recent(opt, obs, nom,
                    tiles=tiles, markmjds=markmjds, show_plot=False,
                    nightly=opt.nightplot, botplanfn=botplanfn, **copilot_plot_args)
        return 0
        
    print('Loading SFD maps...')
    sfd = SFDMap()
    
    if len(args) > 0:
        mp = None
        if (opt.threads or 0) > 1:
            gSFD = sfd
            from astrometry.util.multiproc import multiproc
            mp = multiproc(opt.threads)

        if opt.skip:
            fns = skip_existing_files(args, rawext, primext=opt.primext)
        elif opt.skip_expnum:
            fns = skip_existing_files(args, rawext, by_expnum=True, primext=opt.primext)
        else:
            fns = args
            
        if mp is None:
            for fn in fns:
                process_image(fn, rawext, nom, sfd, opt, obs, tiles)
        else:
            sfd = None
            realobs = obs
            obs = None
            mp.map(bounce_process_image,
                   [(fn, rawext, nom, sfd, opt, obs, tiles) for fn in fns])
            obs = realobs
        plot_recent(opt, obs, nom, tiles=tiles, markmjds=markmjds,
                    show_plot=False, botplanfn=botplanfn)
        return 0

    copilot = Copilot(imagedir, rawext, opt, nom, sfd, obs, tiles, botplanfn,
                      copilot_plot_args)

    # for testability
    if get_copilot:
        return copilot

    copilot.run()
    return 0


class Copilot(NewFileWatcher):
    def __init__(self, imagedir, rawext,
                 opt, nom, sfd, obs, tiles, botplanfn, plot_kwargs):
        self.rawext = rawext
        self.plot_kwargs = plot_kwargs
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
        #print('Filter_new_files:', fns)
        return [fn for fn in fns if
                fn.endswith('.fits.fz') or fn.endswith('.fits')]

    def try_open_file(self, path):
        print('Trying to open file: %s, ext: %s' % (path, self.rawext))
        fitsio.read(path, ext=self.rawext)

    def timed_out(self, dt):
        self.plot_recent()

    def process_file(self, path):
        M = process_image(path, self.rawext, self.nom, self.sfd,
                          self.opt, self.obs, self.tiles)
        if M is None:
            return
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

        plot_recent(self.opt, self.obs, self.nom, tiles=self.tiles,
                    markmjds=markmjds, show_plot=self.opt.show,
                    botplanfn=self.botplanfn, **self.plot_kwargs)

if __name__ == '__main__':
    import obsdb
    from camera import database_filename
    obsdb.django_setup(database_filename=database_filename)
    main()
