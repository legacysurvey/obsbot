from __future__ import print_function
import fitsio
from astrometry.blind.plotstuff import Plotstuff
from astrometry.util.util import Sip, anwcs, anwcs_new_sip, wcs_pv2sip_hdr

def skyplot():
    plot = Plotstuff(size=(800,800), rdw=(103.1, 37.45, 0.8), outformat='png')
    plot.color = 'verydarkblue'
    plot.plot('fill')
    
    for ext in range(1, 17):
        fn = 'mos3.68488.fits'
        hdr = fitsio.read_header(fn, ext=ext)
        wcs = wcs_pv2sip_hdr(hdr)
        plot.color = 'red'
        plot.outline.wcs = anwcs_new_sip(wcs)
        plot.plot('outline')
    
        plot.color = 'white'
        plot.apply_settings()
        rc,dc = wcs.radec_center()
        plot.text_radec(rc, dc, hdr['EXTNAME'])
        
    plot.color = 'white'
    for ext in range(1, 17):
        fn = 'an2/mos3.68488.ext%02i.wcs' % ext
        plot.outline.wcs = anwcs(fn)
        plot.plot('outline')
    
    plot.rgb = (0.2,0.2,0.2)
    plot.plot_grid(0.1, 0.1)
    plot.color = 'gray'
    plot.plot_grid(0.5, 0.5, 0.5, 0.5)
    plot.write('plot.png')


from astrometry.util.fits import fits_table, merge_tables
import numpy as np
from glob import glob
from measure_raw import measure_raw_mosaic3
import os
import pylab as plt


A = fits_table('Almanac_2016-03-31.fits')
A.cut(A.ra_offset < 99)
T = fits_table('db.fits')
T.cut(np.array([e in A.expnum for e in T.expnum]))
print(len(A), len(T))
from camera_mosaic import dradec_to_ref_chip
refdra,refddec = dradec_to_ref_chip(T)

plt.clf()
plt.subplot(2,1,1)
p1 = plt.plot(A.expnum, A.ra_offset, 'b.')
p2 = plt.plot(T.expnum, refdra, 'r.')
plt.figlegend([p1[0],p2[0]], ('Mosstat', 'Copilot'), 'upper right')
plt.ylabel('dRA (arcsec)')
plt.subplot(2,1,2)
plt.plot(A.expnum, A.dec_offset, 'b.')
plt.plot(T.expnum, refddec, 'r.')
plt.ylabel('dDec (arcsec)')
plt.xlabel('Expnum')
plt.suptitle('Copilot and Mosstat, 2016-03-31')
plt.savefig('astrom.png')

plt.clf()
p1 = plt.plot(A.expnum, A.ra_offset - refdra, 'b.')
p2 = plt.plot(A.expnum, A.dec_offset - refddec, 'r.')
plt.figlegend([p1[0],p2[0]], ('RA', 'Dec'), 'upper right')
plt.ylabel('Mosstat - Copilot (arcsec)')
plt.xlabel('Expnum')
plt.suptitle('Copilot and Mosstat, 2016-03-31')
plt.savefig('astrom2.png')

A = fits_table('Almanac_2016-03-20.fits')
A.about()

#print('Extensions', np.unique(A.extname))
A.cut(A.extname == 'im16')
print(len(A), 'im16')
#print('Filenames', A.filename)
#print(' '.join(A.filename))
#print('Expnums:', ', '.join(['%i' % e for e in A.expnum]))

#A = A[:10]

cofn = 'copilot2.fits'
#cofn = 'copilot3.fits'
if not os.path.exists(cofn):
    expnum_map = {}
    fns = glob('/project/projectdirs/cosmo/staging/mosaicz/MZLS_Raw/20160320/*ori.fits.fz')
    for fn in fns:
        hdr = fitsio.read_header(fn)
        expnum = hdr['EXPNUM']
        if not expnum in A.expnum:
            continue
        expnum_map[expnum] = fn
        print('File', fn, 'is expnum', expnum)
        
    MM = []
    for i,t in enumerate(A):
        #extstring = t.extname
        fn = expnum_map[t.expnum]
    
        for ext in range(1, 16+1):
            extstring = 'im%i' % ext
            meas = measure_raw_mosaic3(fn, ext=extstring, n_fwhm=1)
            #print('Measurement:', meas.keys())
            M = fits_table()
            for k in ['airmass', 'extension', 'pixscale', 'nmatched', 'ra_ccd', 'dec_ccd',
                      'band', 'zp', 'rawsky', 'ndetected', 'skybright', 'dy', 'transparency',
                      'seeing', 'dx', 'exptime', 'zp_skysub', 'zp_med', 'zp_med_skysub', 'affine']:
                M.set(k, np.array([meas[k]]))
            M.filename = np.array([fn])
            M.extname = np.array([extstring])
            phdr = meas['primhdr']
            M.expnum = np.array([phdr['EXPNUM']])
            MM.append(M)
    M = merge_tables(MM)

    M.writeto(cofn)

C = fits_table(cofn)

from camera_mosaic import nominal_cal
nom = nominal_cal
C.extension = np.array([ext.strip() for ext in C.extension])
CDs = dict([(ext, nom.cdmatrix(ext)) for ext in np.unique(C.extension)])
C.cd = np.array([CDs[ext] for ext in C.extension])
C.dra  = (C.cd[:,0] * C.dx + C.cd[:,1] * C.dy) * 3600.
C.ddec = (C.cd[:,2] * C.dx + C.cd[:,3] * C.dy) * 3600.



if not 'expnum' in C.columns():
    fns = glob('/project/projectdirs/cosmo/staging/mosaicz/MZLS_Raw/20160320/*ori.fits.fz')
    expnum_map = {}
    for fn in fns:
        hdr = fitsio.read_header(fn)
        expnum = hdr['EXPNUM']
        if not expnum in A.expnum:
            continue
        expnum_map[fn] = expnum
        print('File', fn, 'is expnum', expnum)
    C.expnum = np.array([expnum_map[fn.strip()] for fn in C.filename])
    C.writeto('copilot2.fits')
    
AC = fits_table()
for c in C.columns():
    AC.set(c, [])
        
# First: mosstat vs copilot for im16.
for i in range(len(A)):
    J = np.flatnonzero((C.expnum == A.expnum[i]) * (C.extension == 'im16'))
    assert(len(J) == 1)
    j = J[0]
    for c in C.columns():
        AC.get(c).append(C.get(c)[j])
AC.to_np_arrays()
A.add_columns_from(AC)

A.about()

plt.clf()
p1 = plt.plot(A.zpt, A.zp, 'b.')
p2 = plt.plot(A.zpt, A.zp_med, 'g.')
p3 = plt.plot(A.zpt, A.zp_skysub, 'r.')
p4 = plt.plot(A.zpt, A.zp_med_skysub, 'm.')
plt.legend([p1[0],p2[0],p3[0],p4[0]], ['Current', 'Median', 'Sky', 'Med,Sky'],
           'upper left')
plt.xlabel('Mosstat Zeropoint')
plt.ylabel('Copilot Zeropoint')
plt.savefig('zpt.png')

plt.clf()
p1 = plt.plot(A.zpt, 0.0 + A.zp, 'b.')
p2 = plt.plot(A.zpt, 0.2 + A.zp_med, 'g.')
p3 = plt.plot(A.zpt, 0.4 + A.zp_skysub, 'r.')
p4 = plt.plot(A.zpt, 0.6 + A.zp_med_skysub, 'm.')
plt.legend([p1[0],p2[0],p3[0],p4[0]], ['Current', 'Median', 'Sky', 'Med,Sky'],
           'upper left')
plt.xlabel('Mosstat Zeropoint')
plt.ylabel('Copilot Zeropoint')
plt.savefig('zpt1.png')


plt.clf()
p1 = plt.plot(A.zpt, A.zp - A.zpt, 'b.')
p2 = plt.plot(A.zpt, A.zp_med - A.zpt, 'g.')
p3 = plt.plot(A.zpt, A.zp_skysub - A.zpt, 'r.')
p4 = plt.plot(A.zpt, A.zp_med_skysub - A.zpt, 'm.')
plt.legend([p1[0],p2[0],p3[0],p4[0]], ['Current', 'Median', 'Sky', 'Med,Sky'],
           'lower left')
plt.xlabel('Mosstat Zeropoint')
plt.ylabel('Copilot Zeropoint - Mosstat Zeropoint')
plt.savefig('zpt2.png')


plt.clf()
p1 = plt.plot(A.zpt, 0.6 + A.zp - A.zpt, 'b.')
p2 = plt.plot(A.zpt, 0.4 + A.zp_med - A.zpt, 'g.')
p3 = plt.plot(A.zpt, 0.2 + A.zp_skysub - A.zpt, 'r.')
p4 = plt.plot(A.zpt, 0.0 + A.zp_med_skysub - A.zpt, 'm.')
plt.axhline(0.0, color='k', alpha=0.1)
plt.axhline(0.2, color='k', alpha=0.1)
plt.axhline(0.4, color='k', alpha=0.1)
plt.axhline(0.6, color='k', alpha=0.1)
plt.legend([p1[0],p2[0],p3[0],p4[0]], [
    'Current (%.0f +- %.0f)' % (1000. * np.mean(A.zp - A.zpt),
                                1000. * np.std(A.zp - A.zpt)),
    'Median  (%.0f +- %.0f)' % (1000. * np.mean(A.zp_med - A.zpt),
                                1000. * np.std(A.zp_med - A.zpt)),
    'Sky (%.0f +- %.0f)' % (1000. * np.mean(A.zp_skysub - A.zpt),
                            1000. * np.std(A.zp_skysub - A.zpt)),
    'Med,Sky (%.0f +- %.0f)' % (1000. * np.mean(A.zp_med_skysub - A.zpt),
                                1000. * np.std(A.zp_med_skysub - A.zpt)),],
           'lower left')
plt.xlabel('Mosstat Zeropoint')
plt.ylabel('Copilot Zeropoint - Mosstat Zeropoint')
plt.title('Zeropoints -- 2016-03-20')
plt.savefig('zpt3.png')


#sys.exit(0)

plt.clf()
# plt.subplot(1,2,1)
# plt.plot(A.dx, A.ra_offset, 'b.')
# plt.xlabel('dx')
# plt.ylabel('RA offset')
# plt.subplot(1,2,2)
# plt.plot(A.dy, A.dec_offset, 'b.')
# plt.xlabel('dy')
# plt.ylabel('Dec offset')
p1 = plt.plot(A.ra_offset,  A.dy * 0.262, 'b.')
p2 = plt.plot(A.dec_offset, A.dx * 0.262, 'r.')
plt.legend([p1[0],p2[0]], ['dy,dRA','dx,dDec'], 'upper left')
plt.xlabel('dDec | dRA')
plt.ylabel('dx | dy')
plt.savefig('diff.png')


dra = np.median(A.ra_offset - A.dra)
ddec = np.median(A.dec_offset - A.ddec)
print('Shift in im16 dRA,dDec, arcsec: %.2f, %.2f' % (dra, ddec))
off_dra = dra
off_ddec = ddec

plt.clf()
p1 = plt.plot(A.ra_offset,  A.dra, 'b.')
p2 = plt.plot(A.dec_offset, A.ddec, 'r.')
ax = plt.axis()
xx = np.array([-20,20])
plt.plot(xx, xx - dra, 'b-', alpha=0.2)
plt.plot(xx, xx - ddec, 'r-', alpha=0.2)
plt.axis(ax)
plt.legend([p1[0],p2[0]], ['dRA','dDec'], 'upper left')
plt.xlabel('mosstat')
plt.ylabel('copilot')
plt.title('Mosaic3 im16')
plt.savefig('diffrd.png')



plt.clf()
p1 = plt.plot(A.expnum, A.dx * 0.262, 'b.')
p2 = plt.plot(A.expnum, A.dy * 0.262, 'r.')
p3 = plt.plot(A.expnum, A.ra_offset, 'g.')
p4 = plt.plot(A.expnum, A.dec_offset, 'm.')
plt.legend([p1[0],p2[0],p3[0],p4[0]], ['dx','dy','dRA','dDec'], 'lower right')
plt.xlabel('expnum')
plt.savefig('difft.png')


# Now, look at each copilot extension vs im16.
from astrometry.util.plotutils import PlotSequence
ps = PlotSequence('astrom')


C.affdx = C.affine[:,2]
C.affdy = C.affine[:,5]

refext = 'im16'
Cref = C[C.extension == refext]
print(len(Cref), 'im16 exposures')
ref_expnum = dict([(expnum,i) for i,expnum in enumerate(Cref.expnum)])


plt.clf()
p1 = plt.plot(C.expnum, C.affine[:,3] - 1, 'b.')
p2 = plt.plot(C.expnum, C.affine[:,4], 'r.')
plt.plot(Cref.expnum, Cref.affine[:,4], 'k-')
p3 = plt.plot(C.expnum, C.affine[:,6], 'c.')
p4 = plt.plot(C.expnum, C.affine[:,7] - 1, 'm.')
plt.xlabel('Exposure number')
plt.ylabel('Affine transformation element')
plt.legend((p1[0],p2[0],p3[0],p4[0]),
           ('dX/dx - 1', 'dX/dy', 'dY/dx', 'dY/dy-1'), 'upper right')
plt.title('Mosaic astrometry -- 2016-03-20')
ps.savefig()


C.nomzp = np.array([nom.zeropoint(b, ext=ext.strip())
                    for b,ext in zip(C.band, C.extension)])

CICR = []
#for ext in [4, 16]:
for ext in range(1, 16):
    Ci = C[C.extension == 'im%i' % ext]
    print('Extension', ext, ':', len(Ci), 'exposures')
    Cr = Cref[np.array([ref_expnum[expnum] for expnum in Ci.expnum])]

    CICR.append((ext,Ci,Cr))

# Photometric zeropoint offsets
for ext,Ci,Cr in CICR:
    plt.clf()
    plt.hist(Ci.zp_med_skysub - Cr.zp_med_skysub, bins=20)
    plt.xlabel('Zeropoint of ext %i vs ref im16' % ext)
    ps.savefig()

# Photometric zeropoint offsets
for ext,Ci,Cr in CICR:
    plt.clf()
    plt.hist(Ci.zp_med_skysub - Ci.nomzp, bins=20)
    plt.xlabel('Zeropoint of ext %i vs nominal' % ext)
    ps.savefig()

sys.exit(0)    


for ext,Ci,Cr in CICR:
    from camera_mosaic import dradec_to_ref_chip
    Ci.affine_x0  = Ci.affine[:,0]
    Ci.affine_y0  = Ci.affine[:,1]
    Ci.affine_dx  = Ci.affine[:,2]
    Ci.affine_dxx = Ci.affine[:,3]
    Ci.affine_dxy = Ci.affine[:,4]
    Ci.affine_dy  = Ci.affine[:,5]
    Ci.affine_dyx = Ci.affine[:,6]
    Ci.affine_dyy = Ci.affine[:,7]
    
    cdra,cddec = dradec_to_ref_chip(Ci, refext=refext)

    amap = dict([(expnum, i) for i,expnum in enumerate(A.expnum)])
    ai = np.array([amap[e] for e in Ci.expnum])
    adra  = A.ra_offset[ai]
    addec = A.dec_offset[ai]

    print('Ci expnum:', Ci.expnum)
    print('A expnum:', A.expnum[ai])

    
    plt.clf()
    plt.plot(adra,  cdra - adra, 'b.')
    plt.plot(addec, cddec- addec, 'r.')
    plt.axhline(0, color='k', alpha=0.1)
    plt.legend([p1[0],p2[0]], ['dRA','dDec'], 'upper left')
    plt.xlabel('Mosstat %s offset (arcsec)' % refext)
    plt.ylabel('Copilot im%i - Mosstat %s' % (ext, refext))
    plt.title('Mosaic3: Copilot offsets w/ rotation (arcsec)')
    ps.savefig()

    
    # Take the offset between the reference chip CRPIX and this chip's
    # CRPIX and push that through the affine rotation matrix.
    
    # Assume the affine rotation elements are due to a whole-camera
    # rigid rotation.  Push this through the difference in CRPIX
    # values (ie, distance from boresight) through that rotation
    # matrix to convert an offset in imX to an offset in im16.
    (refcrx, refcry) = nom.crpix(refext)
    (crx,    cry)    = nom.crpix('im%i' % ext)
    dcrx = refcrx - crx
    dcry = refcry - cry
    dcx = Ci.affine[:, 3] * dcrx + Ci.affine[:,4] * dcry - dcrx
    dcy = Ci.affine[:, 6] * dcrx + Ci.affine[:,7] * dcry - dcry
    print('Scatter:', np.std(Cr.dx - (Ci.dx - dcx)), np.std(Cr.dy - (Ci.dy - dcy)))
    
    drai = np.median(Cr.dra - Ci.dra)
    ddeci = np.median(Cr.ddec - Ci.ddec)

    ax = [-20,20,-20,20]
    
    plt.clf()
    plt.plot(Cr.dra, Ci.dra, 'b.')
    plt.plot(Cr.ddec, Ci.ddec, 'r.')
    #ax = plt.axis()
    xx = np.array([-20,20])
    p1 = plt.plot(xx, xx - drai, 'b-', alpha=0.2)
    p2 = plt.plot(xx, xx - ddeci, 'r-', alpha=0.2)
    plt.axis(ax)
    plt.legend([p1[0],p2[0]], ['dRA','dDec'], 'upper left')
    plt.xlabel(refext)
    plt.ylabel('im%i' % ext)
    plt.title('Mosaic3: Copilot offsets (arcsec)')
    ps.savefig()
    
    
    # plt.clf()
    # plt.plot(Cr.dx, Ci.dx, 'b.')
    # plt.plot(Cr.dy, Ci.dy, 'r.')
    # 
    # plt.plot(Cr.dx, Ci.dx - dcx, 'c.')
    # plt.plot(Cr.dy, Ci.dy - dcy, 'm.')
    # 
    # dxi = np.median(Cr.dx - Ci.dx)
    # dyi = np.median(Cr.dy - Ci.dy)
    # 
    # ax = plt.axis()
    # xx = np.array([-100,100])
    # p1 = plt.plot(xx, xx - dxi, 'b-', alpha=0.2)
    # p2 = plt.plot(xx, xx - dyi, 'r-', alpha=0.2)
    # plt.axis(ax)
    # plt.legend([p1[0],p2[0]], ['dx','dy'], 'upper left')
    # plt.xlabel(refext)
    # plt.ylabel('im%i' % ext)
    # plt.title('Mosaic3: Copilot offsets (pixels)')
    # ps.savefig()
    
    cdx = Ci.dx - dcx
    cdy = Ci.dy - dcy
    
    cdra  = (Ci.cd[:,0] * cdx + Ci.cd[:,1] * cdy) * 3600.
    cddec = (Ci.cd[:,2] * cdx + Ci.cd[:,3] * cdy) * 3600.
    
    plt.clf()
    plt.plot(Cr.dra,  cdra, 'b.')
    plt.plot(Cr.ddec, cddec, 'r.')
    #ax = plt.axis()
    xx = np.array([-20,20])
    p1 = plt.plot(xx, xx - np.median(Cr.dra - cdra), 'b-', alpha=0.2)
    p2 = plt.plot(xx, xx - np.median(Cr.ddec - cddec), 'r-', alpha=0.2)
    plt.axis(ax)
    plt.legend([p1[0],p2[0]], ['dRA','dDec'], 'upper left')
    plt.xlabel('Copilot %s' % refext)
    plt.ylabel('Copilot im%i' % ext)
    plt.title('Mosaic3: Copilot offsets w/ rotation (arcsec)')
    ps.savefig()


    plt.clf()
    plt.plot(adra,  cdra  + off_dra, 'b.')
    plt.plot(addec, cddec + off_ddec, 'r.')
    plt.axis(ax)
    xx = np.array([-20,20])
    plt.plot(xx, xx, 'k-', alpha=0.2)
    plt.xlabel('Mosstat %s' % refext)
    plt.ylabel('Copilot im%i -> Mosstat %s' % (ext, refext))
    plt.title('Mosaic3: Copilot offsets w/ rotation (arcsec)')
    ps.savefig()

    
    plt.clf()
    plt.plot(Cr.dra,  cdra - Cr.dra, 'b.')
    plt.plot(Cr.ddec, cddec- Cr.ddec, 'r.')
    plt.axhline(0, color='k', alpha=0.1)
    plt.legend([p1[0],p2[0]], ['dRA','dDec'], 'upper left')
    plt.xlabel(refext)
    plt.ylabel('im%i - %s' % (ext, refext))
    plt.title('Mosaic3: Copilot offsets w/ rotation (arcsec)')
    ps.savefig()
    
    # plt.clf()
    # plt.plot(Cr.affdx, Ci.affdx, 'b.')
    # plt.plot(Cr.affdy, Ci.affdy, 'r.')
    # ax = plt.axis()
    # xx = np.array([-100,100])
    # p1 = plt.plot(xx, xx - np.median(Cr.affdx - Ci.affdx), 'b-', alpha=0.2)
    # p2 = plt.plot(xx, xx - np.median(Cr.affdy - Ci.affdy), 'r-', alpha=0.2)
    # plt.axis(ax)
    # plt.legend([p1[0],p2[0]], ['dx','dy'], 'upper left')
    # plt.xlabel(refext)
    # plt.ylabel('im%i' % ext)
    # plt.title('Mosaic3: Copilot offsets (affine; pixels)')
    # ps.savefig()

    
    
sys.exit(0)    
    
for ext,Ci,Cr in CICR:
    dxi = np.median(Cr.dx - Ci.dx)
    dyi = np.median(Cr.dy - Ci.dy)

    A = np.zeros((len(Ci), 3))
    A[:,0] = 1.
    A[:,1] = Ci.dx
    A[:,2] = Ci.dy

    R = np.linalg.lstsq(A, Cr.dx)
    resx = R[0]
    print('Im16 dx = (1,dx,dy) *', resx)

    R = np.linalg.lstsq(A, Cr.dy)
    resy = R[0]
    print('Im16 dy = (1,dx,dy) *', resy)

    fitx = resx[0] + Ci.dx * resx[1] + Ci.dy * resx[2]
    fity = resy[0] + Ci.dx * resy[1] + Ci.dy * resy[2]
    
    plt.clf()
    p1 = plt.plot(Cr.dx, Ci.dx, 'b.')
    p2 = plt.plot(Cr.dy, Ci.dy, 'r.')
    plt.plot(Cr.dx, fitx, 'c.')
    plt.plot(Cr.dy, fity, 'm.')
    ax = plt.axis()
    xx = np.array([-100,100])
    plt.plot(xx, xx - dxi, 'b-', alpha=0.2)
    plt.plot(xx, xx - dyi, 'r-', alpha=0.2)
    plt.legend([p1[0],p2[0]], ['dx','dy'], 'upper left')
    plt.axis(ax)
    plt.xlabel(refext)
    plt.ylabel('im%i' % ext)
    plt.title('Mosaic3: Copilot offsets')
    ps.savefig()
    
