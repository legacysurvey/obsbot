from __future__ import print_function
import os
import tempfile
import fitsio
import numpy as np
import matplotlib
matplotlib.use('Agg')
import pylab as plt

from scipy.stats import sigmaclip
from scipy.ndimage.filters import gaussian_filter
from scipy.ndimage.measurements import label, find_objects, center_of_mass
from scipy.ndimage.filters import median_filter

from astrometry.util.plotutils import *
from astrometry.util.util import wcs_pv2sip_hdr
from astrometry.libkd.spherematch import match_xy

from legacyanalysis.ps1cat import ps1cat, ps1_to_decam

import photutils

def sensible_sigmaclip(arr, nsigma = 4.):
    goodpix,lo,hi = sigmaclip(arr, low=nsigma, high=nsigma)
    # sigmaclip returns unclipped pixels, lo,hi, where lo,hi are
    # mean(goodpix) +- nsigma * sigma
    meanval = np.mean(goodpix)
    sigma = (meanval - lo) / nsigma
    return meanval, sigma

def parse_section(s, slices=False):
    '''
    parse '[57:2104,51:4146]' into integers; also subtract 1.
    '''
    s = s.replace('[','').replace(']','').replace(',',' ').replace(':',' ')
    #print('String', s)
    i = [int(si)-1 for si in s.split()]
    assert(len(i) == 4)
    if not slices:
        return i
    slc = slice(i[2], i[3]+1), slice(i[0], i[1]+1)
    #print('Slice', slc)
    return slc

def read_raw_decam(F, ext):
    '''
    F: fitsio FITS object
    '''
    img = F[ext].read()
    hdr = F[ext].read_header()
    print('Image type', img.dtype, img.shape)

    img = img.astype(np.float32)
    print('Converted image to', img.dtype, img.shape)

    if False:
        mn,mx = np.percentile(img.ravel(), [25,95])
        kwa = dict(vmin=mn, vmax=mx)

        plt.clf()
        dimshow(img, **kwa)
        plt.title('Raw image')
        ps.savefig()

        M = 200
        plt.clf()
        plt.subplot(2,2,1)
        dimshow(img[-M:, :M], ticks=False, **kwa)
        plt.subplot(2,2,2)
        dimshow(img[-M:, -M:], ticks=False, **kwa)
        plt.subplot(2,2,3)
        dimshow(img[:M, :M], ticks=False, **kwa)
        plt.subplot(2,2,4)
        dimshow(img[:M, -M:], ticks=False, **kwa)
        plt.suptitle('Raw corners')
        ps.savefig()
    
    if 'DESBIAS' in hdr:
        assert(False)
    # DECam RAW image

    # Raw image size 2160 x 4146

    # Subtract median overscan and multiply by gains 
    # print('DATASECA', hdr['DATASECA'])
    # print('BIASSECA', hdr['BIASSECA'])
    # print('DATASECB', hdr['DATASECB'])
    # print('BIASSECB', hdr['BIASSECB'])
    dataA = parse_section(hdr['DATASECA'], slices=True)
    biasA = parse_section(hdr['BIASSECA'], slices=True)
    dataB = parse_section(hdr['DATASECB'], slices=True)
    biasB = parse_section(hdr['BIASSECB'], slices=True)
    gainA = hdr['GAINA']
    gainB = hdr['GAINB']
    # print('DataA', dataA)
    # print('BiasA', biasA)
    # print('DataB', dataB)
    # print('BiasB', biasB)

    if False:
        plt.clf()
        plt.plot([np.median(img[i,:] )     for i in range(100)], 'b-')
        plt.plot([np.median(img[-(i+1),:]) for i in range(100)], 'c-')
        plt.plot([np.median(img[:,i] )     for i in range(100)], 'r-')
        plt.plot([np.median(img[:,-(i+1)]) for i in range(100)], 'm-')
        plt.title('Img')
        ps.savefig()
    
        plt.clf()
        plt.plot([np.median(img[dataA][i,:] )     for i in range(100)], 'b-')
        plt.plot([np.median(img[dataA][-(i+1),:]) for i in range(100)], 'c-')
        plt.plot([np.median(img[dataA][:,i] )     for i in range(100)], 'r-')
        plt.plot([np.median(img[dataA][:,-(i+1)]) for i in range(100)], 'm-')
        plt.title('Img DataA')
        ps.savefig()
    
        plt.clf()
        plt.plot([np.median(img[dataB][i,:] )     for i in range(100)], 'b-')
        plt.plot([np.median(img[dataB][-(i+1),:]) for i in range(100)], 'c-')
        plt.plot([np.median(img[dataB][:,i] )     for i in range(100)], 'r-')
        plt.plot([np.median(img[dataB][:,-(i+1)]) for i in range(100)], 'm-')
        plt.title('Img DataB')
        ps.savefig()

    
    img[dataA] = (img[dataA] - np.median(img[biasA])) * gainA
    img[dataB] = (img[dataB] - np.median(img[biasB])) * gainB
    
    # Trim the image -- could just take the min/max of TRIMSECA/TRIMSECB...
    trimA = parse_section(hdr['TRIMSECA'], slices=True)
    trimB = parse_section(hdr['TRIMSECB'], slices=True)
    # copy the TRIM A,B sections into a new image...
    trimg = np.zeros_like(img)
    trimg[trimA] = img[trimA]
    trimg[trimB] = img[trimB]
    # ... and then cut that new image
    trim = parse_section(hdr['TRIMSEC'], slices=True)
    img = trimg[trim]
    print('Trimmed image:', img.dtype, img.shape)

    return img,hdr
    

def measure_raw_decam(fn, ext='N4'):

    ps = PlotSequence('raw')
    
    # aperture phot radii in arcsec
    aprad = 3.5
    skyrad = [7., 10.]

    minstar = 5
    pixsc = 0.262
    maxshift = 30 # arcsec to PS1
    
    F = fitsio.FITS(fn)
    primhdr = F[0].read_header()

    img,hdr = read_raw_decam(F, ext)

    M = 200
    mn,mx = np.percentile(img.ravel(), [25,98])
    kwa = dict(vmin=mn, vmax=mx)
    
    plt.clf()
    dimshow(img, **kwa)
    plt.title('Trimmed image')
    ps.savefig()

    plt.clf()
    plt.subplot(2,2,1)
    dimshow(img[-M:, :M], ticks=False, **kwa)
    plt.subplot(2,2,2)
    dimshow(img[-M:, -M:], ticks=False, **kwa)
    plt.subplot(2,2,3)
    dimshow(img[:M, :M], ticks=False, **kwa)
    plt.subplot(2,2,4)
    dimshow(img[:M, -M:], ticks=False, **kwa)
    plt.suptitle('Trimmed corners')
    ps.savefig()

    expnum = primhdr['EXPNUM']
    print('Expnum', expnum)
    obstype = primhdr['OBSTYPE']
    print('Obstype', obstype)
    band = primhdr['FILTER']
    band = band.split()[0]
    print('Band', band)
    exptime = primhdr['EXPTIME']
    airmass = primhdr['AIRMASS']
    print('Exptime', exptime)
    print('Airmass', airmass)

    nominal_cal = dict(
        g = (26.610,
             22.04,
             0.17,),
        r = (26.818,
             20.91,
             0.10,),
        z = (26.484,
             18.46,
             0.06,),
        )
    
    zp0, sky0, kx = nominal_cal[band]


    # Find the sky value and noise level
    sky,sig1 = sensible_sigmaclip(img[1500:2500, 500:1000])

    skybr = -2.5 * np.log10(sky/pixsc/pixsc/exptime) + zp0
    print('Sky brightness: %8.2f mag/arcsec^2' % skybr)
    print('Fiducial:       %8.2f mag/arcsec^2' % sky0)

    img -= sky
    
    # Ugly removal of sky gradients by subtracting median in first x and then y
    H,W = img.shape
    meds = np.array([np.median(img[:,i]) for i in range(W)])
    meds = median_filter(meds, size=5)
    img -= meds[np.newaxis,:]
    meds = np.array([np.median(img[i,:]) for i in range(H)])
    meds = median_filter(meds, size=5)
    img -= meds[:,np.newaxis]

    mn,mx = np.percentile(img.ravel(), [25,98])
    kwa = dict(vmin=mn, vmax=mx)
    
    plt.clf()
    dimshow(img, **kwa)
    plt.title('Sky-sub image')
    ps.savefig()
    
    # Trim off some extra pixels -- image edges are often bright...
    trim = 100
    # central part of the image
    cimg = img[trim:-trim, trim:-trim]
    
    # Detect & aperture-photometer stars
    fwhm = 5.0
    psfsig = fwhm / 2.35
    psfnorm = 1./(2. * np.sqrt(np.pi) * psfsig)
    detsn = gaussian_filter(cimg / sig1, psfsig) / psfnorm

    plt.clf()
    dimshow(detsn, vmin=-3, vmax=50, cmap='hot')
    plt.colorbar()
    plt.title('Detection S/N')
    ps.savefig()

    # zero out the edges -- larger margin here?
    detsn[0 ,:] = 0
    detsn[:, 0] = 0
    detsn[-1,:] = 0
    detsn[:,-1] = 0

    # Detection threshold
    det_thresh = 20
    peaks = (detsn > det_thresh)

    # "Peak" region to centroid
    P = 5

    # HACK -- Just keep the brightest pixel in each blob!
    blobs,nblobs = label(peaks)
    slices = find_objects(blobs)
    xx,yy = [],[]
    fx,fy = [],[]
    for i,slc in enumerate(slices):
        y0 = slc[0].start
        x0 = slc[1].start
        subimg = detsn[slc]
        imax = np.argmax(subimg)
        y,x = np.unravel_index(imax, subimg.shape)
        if (x0+x) < P or (x0+x) > W-1-P or (y0+y) < P or (y0+y) > H-1-P:
            print('Skipping edge peak', x0+x, y0+y)
            continue
        xx.append(x0 + x)
        yy.append(y0 + y)
        pkarea = detsn[y0+y-P: y0+y+P+1, x0+x-P: x0+x+P+1]
        cy,cx = center_of_mass(pkarea)
        #print('Center of mass', cx,cy)
        fx.append(x0+x-P+cx)
        fy.append(y0+y-P+cy)
        #print('x,y', x0+x, y0+y, 'vs centroid', x0+x-P+cx, y0+y-P+cy)

    fx = np.array(fx)
    fy = np.array(fy)
    xx = np.array(xx)
    yy = np.array(yy)
        
    plt.clf()
    dimshow(detsn, vmin=-3, vmax=50, cmap='hot')
    ax = plt.axis()
    plt.plot(xx, yy, 'go', mec='g', mfc='none', ms=10)
    plt.colorbar()
    plt.title('Detected sources')
    plt.axis(ax)
    ps.savefig()

    plt.clf()
    dimshow(detsn, vmin=-3, vmax=50, cmap='hot')
    ax = plt.axis()
    plt.plot(fx, fy, 'go', mec='g', mfc='none', ms=10)
    plt.colorbar()
    plt.title('Detected sources (2)')
    plt.axis(ax)
    ps.savefig()

    keep = (np.hypot(fx - xx, fy - yy) < 1)
    print(sum(keep), 'of', len(keep), 'stars have centroids within 1 of peaks')
    print('mean dx', np.mean(fx-xx), 'dy', np.mean(fy-yy))
    
    # Cut down to stars whose centroids are within 1 pixel of their peaks...
    assert(float(sum(keep)) / len(keep) > 0.9)

    fx = fx[keep]
    fy = fy[keep]

    # we trimmed the image before running detection; re-add that margin
    fx += trim
    fy += trim
    
    #aps = np.append(aprad, skyrad).astype(np.float32)
    #aps /= pixsc
    apxy = np.vstack((fx, fy)).T
    ap = []
    #for rad in aps:
    aprad_pix = aprad / pixsc
    aper = photutils.CircularAperture(apxy, aprad_pix)
    p = photutils.aperture_photometry(img, aper)
    apflux = p.field('aperture_sum')
    #apflux,sky1,sky2 = ap

    # Manual aperture photometry to get clipped means in sky annulus
    sky_inner_r, sky_outer_r = [r / pixsc for r in skyrad]
    sky = []
    for xi,yi in zip(fx,fy):
        ix = int(np.round(xi))
        iy = int(np.round(yi))
        skyR = sky_outer_r
        xlo = max(0, ix-skyR)
        xhi = min(W, ix+skyR+1)
        ylo = max(0, iy-skyR)
        yhi = min(H, iy+skyR+1)
        xx,yy = np.meshgrid(np.arange(xlo,xhi), np.arange(ylo,yhi))
        r2 = (xx - xi)**2 + (yy - yi)**2
        inannulus = ((r2 >= sky_inner_r**2) * (r2 < sky_outer_r**2))
        skypix = img[ylo:yhi, xlo:xhi][inannulus]
        s,nil = sensible_sigmaclip(skypix)
        sky.append(s)
    sky = np.array(sky)

    apflux2 = apflux - sky * (np.pi * aprad_pix**2)

    # HACK -- convert TPV WCS header to SIP.
    wcs = wcs_pv2sip_hdr(hdr)
    print('Converted WCS to', wcs)

    ra_ccd,dec_ccd = wcs.pixelxy2radec((W+1)/2., (H+1)/2.)
    
    # Read in the PS1 catalog, and keep those within 0.25 deg of CCD center
    # and those with main sequence colors
    pscat = ps1cat(ccdwcs=wcs)
    stars = pscat.get_stars()
    print('Got PS1 stars:', stars)

    stars.gicolor = stars.median[:,0] - stars.median[:,2]
    keep = (stars.gicolor > 0.4) * (stars.gicolor < 2.7)
    stars.cut(keep)
    if len(stars) == 0:
        print('No overlap or too few stars in PS1')
        return None

    ok,px,py = wcs.radec2pixelxy(stars.ra, stars.dec)
    px -= 1
    py -= 1

    kwa = dict(vmin=-3*sig1, vmax=50*sig1, cmap='hot')
    plt.clf()
    dimshow(img, **kwa)
    ax = plt.axis()
    plt.plot(fx, fy, 'go', mec='g', mfc='none', ms=10)
    plt.plot(px, py, 'm.')
    plt.axis(ax)
    plt.title('PS1 stars')
    ps.savefig()
    
    # Match PS1 to our detections, find offset
    radius = maxshift / pixsc
    I,J,d = match_xy(px, py, fx, fy, radius)
    print(len(I), 'matches')

    plt.clf()
    plothist(px[I] - fx[J], py[I] - fy[J])
    plt.xlabel('dx (pixels)')
    plt.ylabel('dy (pixels)')
    plt.title('DECam to PS1 matches')
    ps.savefig()

    dx = px[I] - fx[J]
    dy = py[I] - fy[J]

    # Hmm, let's try being dumb and see if that works...
    shiftx = np.median(dx)
    shifty = np.median(dy)

    ax = plt.axis()
    plt.plot(shiftx, shifty, 'm.')
    plt.axis(ax)
    ps.savefig()

    # Refine with smaller search radius
    radius2 = 3. / pixsc
    I,J,d = match_xy(px, py, fx+shiftx, fy+shifty, radius2)
    
    dx = px[I] - fx[J]
    dy = py[I] - fy[J]
    shiftx = np.median(dx)
    shifty = np.median(dy)

    plt.clf()
    plothist(dx, dy)
    plt.xlabel('dx (pixels)')
    plt.ylabel('dy (pixels)')
    plt.title('DECam to PS1 matches')
    ax = plt.axis()
    plt.plot(shiftx, shifty, 'm.')
    plt.axis(ax)
    ps.savefig()

    # Compute photometric offset compared to PS1
    # as the PS1 minus DECam-observed mags
    colorterm = ps1_to_decam(stars.median[I,:], band)

    ps1band = ps1cat.ps1band[band]
    ps1mag = stars.median[I, ps1band] + colorterm
    
    plt.clf()
    plt.semilogy(ps1mag, apflux2[J], 'b.')
    plt.xlabel('PS1 mag')
    plt.ylabel('DECam ap flux')
    ps.savefig()

    plt.clf()
    plt.semilogy(ps1mag, apflux[J], 'b.')
    plt.xlabel('PS1 mag')
    plt.ylabel('DECam ap flux 1')
    ps.savefig()

    
    apmag = -2.5 * np.log10(apflux2) + zp0 + 2.5 * np.log10(exptime)

    plt.clf()
    plt.plot(ps1mag, apmag[J], 'b.')
    plt.xlabel('PS1 mag')
    plt.ylabel('DECam ap mag')
    ps.savefig()

    dmag = ps1mag - apmag[J]
    dm,dsig = sensible_sigmaclip(dmag, nsigma=2.5)
    print('Mag offset', dm)
    print('Scatter', dsig)
    
    # dmagA = dmag[iA]
    # dmagB = dmag[iB]
    # djs_iterstat, dmag, median=mag_offset, sigma=mag_rms, sigrej=2.5
    # zpt_observed = zpt0 + mag_offset
    # djs_iterstat, dmagA, median=mag_offsetA, sigma=mag_rmsA, sigrej=2.5
    # djs_iterstat, dmagB, median=mag_offsetB, sigma=mag_rmsB, sigrej=2.5
    # zpt_observedA = zpt0 + mag_offsetA
    # zpt_observedB = zpt0 + mag_offsetB
    # transparency = 10.d0^(-0.4*(zpt0-zpt_observed - kx*(1.0-airmass)))


    
    
measure_raw_decam('DECam_00488199.fits.fz')
    
