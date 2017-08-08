from __future__ import print_function

from astrometry.util.fits import *

from astrometry.blind.plotstuff import Plotstuff
from astrometry.util.util import Sip, anwcs, anwcs_new_sip, wcs_pv2sip_hdr, anwcs_new_tan, Tan
from astrometry.util.fits import fits_table
from astrometry.libkd.spherematch import match_radec
from astrometry.util.plotutils import *
import fitsio
import numpy as np

from mosbot import Mosbot
import camera_mosaic

import obsdb
obsdb.django_setup()
from obsdb.models import ComputedExptime, OtherPasses


band = 'z'
nom = camera_mosaic.nominal_cal
fid = nom.fiducial_exptime(band)
target = fid.single_exposure_depth

threshold = 0.25

def depth_to_factor(depth):
    shortfall = target - depth
    factor = (10.**(-shortfall / 2.5))**2
    return factor


def db_to_fits(adj):
    T = fits_table()
    for field in ['starttime', 'seqnum', 'tileid', 'passnumber', 'band',
                  'airmass', 'ebv', 'meas_band', 'zeropoint',
                  'transparency', 'seeing', 'sky', 'expfactor',
                  'adjfactor', 'exptime_unclipped', 'exptime_clipped',
                  'exptime_satclipped', 'exptime']:
        g = getattr(adj[0], field)
        #if isinstance(g, basestring):
        if str(g) == g:
            T.set(field, np.array([str(getattr(m, field)) for m in adj]))
        else:
            T.set(field, np.array([getattr(m, field) for m in adj]))
    return T

def mosaic_wcs(ra, dec, pixbin=1.):
    # This is pretty close to the outline of the four Mosaic chips.
    W = H = (4096 * 2 + 100) / pixbin
    cd = pixbin * 0.262 / 3600.
    tan = Tan(ra, dec, W/2., H/2., cd, 0., 0., cd,
              float(W), float(H))
    return tan

def plot_exposure(plot, ra, dec):
    wcs = mosaic_wcs(ra, dec)
    plot.outline.wcs = anwcs_new_tan(wcs)
    plot.plot('outline')


mjd0 = 57789.0

db = ComputedExptime.objects.filter(starttime__gt=mjd0, starttime__lt=mjd0+0.8)
print(db.count(), 'entries from', mjd0)
adj = db.filter(adjfactor__gt=1.)
print(adj.count(), 'entries with adjfactor > 1')

#adj = db_to_fits(adj)

ps = PlotSequence('adjust', format='%03i')

tiles = fits_table('obstatus/mosaic-tiles_obstatus.fits')
tiles.cut(tiles.get('pass') <= 3)
tiles.cut(tiles.dec >= 30)

for a in adj[:10]:
    tileid = a.tileid
    print()
    print('tile id', tileid)
    i = np.flatnonzero(tiles.tileid == tileid)
    tile = tiles[i[0]]
    ra,dec = tile.ra, tile.dec
    tilepass = tile.get('pass')
    
    others = a.otherpasses_set.all()
    print(others.count(), 'other passes involved')
    others = others.filter(depth__gt=1, depth__lt=30)
    print(others.count(), 'with measured depths')

    o = fits_table()
    o.tileid  = np.array([t.tileid     for t in others])
    o.passnum = np.array([t.passnumber for t in others])
    o.depth   = np.array([t.depth      for t in others])
    others = o
    #print('tileids:', [t.tileid for t in others])
    tid,I = np.unique(others.tileid, return_index=True)
    others.cut(I)
    print('Cut to', len(others), 'tiles based on unique tileids')

    others.factor = [depth_to_factor(t.depth) for t in others]
    # look up ra,dec
    others.ra  = np.zeros(len(others))
    others.dec = np.zeros(len(others))
    for ii,t in enumerate(others):
        i = np.flatnonzero(tiles.tileid == t.tileid)
        i = i[0]
        others.ra [ii] = tiles.ra [i]
        others.dec[ii] = tiles.dec[i]

    print('tileids:', others.tileid)
    print('passes:', others.passnum)
    #print('depths:', [t.depth for t in others])
    print('factors:', ', '.join(['%.02f'%f for f in others.factor]))

    if False:
        PW,PH = 800,800
        plot = Plotstuff(size=(PW, PH), rdw=(ra, dec, 2), outformat='png')
        plot.color = 'verydarkblue'
        plot.plot('fill')
        
        plot.outline.fill = False
        plot.color = 'red'
        plot_exposure(plot, tile.ra, tile.dec)
    
        plot.color = 'white'
        plot.outline.fill = True
        for t in others:
            plot.alpha = 0.25 * t.factor
            plot.apply_settings()
            plot_exposure(plot, t.ra, t.dec)
    
        plot.write(ps.getnext())

    pixbin = 8
    mywcs = mosaic_wcs(tile.ra, tile.dec, pixbin=pixbin)
    H,W = mywcs.shape

    haspass = dict([(p, np.zeros((H,W), bool)) for p in [1,2,3]])
    covs = dict([(p, np.zeros((H,W), np.float32)) for p in [1,2,3]])
    cov = np.zeros((H,W), np.float32)

    for t in others:
        if t.factor == 0:
            continue
        ok,x,y = mywcs.radec2pixelxy(t.ra, t.dec)
        xlo = np.clip(int(x - W/2), 0, W)
        xhi = np.clip(int(x + W/2), 0, W)
        ylo = np.clip(int(y - H/2), 0, H)
        yhi = np.clip(int(y + H/2), 0, H)
        if xlo == xhi or ylo == yhi:
            continue
        haspass[t.passnum][ylo:yhi, xlo:xhi] = True
        covs[t.passnum][ylo:yhi, xlo:xhi] += t.factor
        cov [ylo:yhi, xlo:xhi] += t.factor
    
    # Previous exposure for this tile
    depth = tile.get('%s_depth' % band)
    shortfall = target - depth
    if depth == 30:
        oldfactor = 1.
    elif shortfall > threshold:
        oldfactor = 0.
    else:
        oldfactor = (10.**(-shortfall / 2.5))**2

    if oldfactor > 0:
        p = tilepass
        haspass[p][:,:] = True
        covs[p] += oldfactor
        cov  += oldfactor

    ncov = haspass[1]*1 + haspass[2]*1 + haspass[3]*1

    def makeplots(tt):

        plt.subplots_adjust(hspace=0.1)

        plt.clf()
        for p in [1,2,3]:
            plt.subplot(2,2,p)
            plt.imshow(covs[p], interpolation='nearest', origin='lower',
                       vmin=0, vmax=2, cmap='RdBu')
            plt.colorbar(ticks=[0,1,2])
            plt.title('Pass %i' % p)
            plt.xticks([]); plt.yticks([])
        plt.subplot(2,2,4)
        plt.imshow(cov - ncov, interpolation='nearest', origin='lower',
                   vmin=-1, vmax=1, cmap='RdBu')
        plt.xticks([]); plt.yticks([])
        plt.title('Total Depth - N passes')
        plt.colorbar(ticks=[-1,0,1])
        plt.suptitle('Depth factor: ' + tt)
        ps.savefig()

        hmax = 5
        ha = dict(range=(0,hmax), bins=40, histtype='step')
        cmap = { 1:'r', 2:'g', 3:'b', 4:'m', 5:'c' }

        plt.subplots_adjust(hspace=0)
        
        plt.clf()
        for p in [1,2,3]:
            plt.subplot(4,1,p)
            plt.hist(covs[p].ravel(), color=cmap.get(p, 'k'),
                     label='Pass %i' % (p), **ha)
            plt.legend()
            plt.xticks([])
            plt.yticks([])
            plt.axvline(1, color='k', lw=2, alpha=0.5)
        plt.subplot(4,1,4)
        plt.hist(np.clip(cov.ravel(), 0, hmax), color='m',
                 label='Total depth', **ha)
        plt.axvline(1, color='k', lw=2, alpha=0.5)
        plt.axvline(2, color='k', lw=2, alpha=0.5)
        plt.axvline(3, color='k', lw=2, alpha=0.5)
        plt.yticks([])
        plt.xticks(np.arange(hmax+1))
        plt.legend()
        plt.suptitle(tt)
        ps.savefig()

        
    makeplots('Before (pass %i)' % tilepass)
    
    factor = a.adjfactor
    p = tilepass
    haspass[p][:,:] = True
    cov += factor
    covs[p] += factor
    ncov = haspass[1]*1 + haspass[2]*1 + haspass[3]*1

    makeplots('After (factor = %.2f, pass = %i)' % (factor, tilepass))
