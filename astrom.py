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

T = fits_table('Almanac_2016-03-20.fits')
T.about()

#print('Extensions', np.unique(T.extname))
T.cut(T.extname == 'im16')
print(len(T), 'im16')
#print('Filenames', T.filename)
#print(' '.join(T.filename))
#print('Expnums:', ', '.join(['%i' % e for e in T.expnum]))


cofn = 'copilot2.fits'
if not os.path.exists(cofn):
    expnum_map = {}
    fns = glob('/project/projectdirs/cosmo/staging/mosaicz/MZLS_Raw/20160320/*ori.fits.fz')
    for fn in fns:
        hdr = fitsio.read_header(fn)
        expnum = hdr['EXPNUM']
        if not expnum in T.expnum:
            continue
        expnum_map[expnum] = fn
        print('File', fn, 'is expnum', expnum)
        
    MM = []
    for i,t in enumerate(T[:2]):
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
            phdr = M['primhdr']
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
### Are these the right way around?
C.dra  = (C.cd[:,0] * C.dx + C.cd[:,1] * C.dy) * 3600.
C.ddec = (C.cd[:,2] * C.dx + C.cd[:,3] * C.dy) * 3600.



if not 'expnum' in C.columns():
    fns = glob('/project/projectdirs/cosmo/staging/mosaicz/MZLS_Raw/20160320/*ori.fits.fz')
    expnum_map = {}
    for fn in fns:
        hdr = fitsio.read_header(fn)
        expnum = hdr['EXPNUM']
        if not expnum in T.expnum:
            continue
        expnum_map[fn] = expnum
        print('File', fn, 'is expnum', expnum)
    C.expnum = np.array([expnum_map[fn.strip()] for fn in C.filename])
    C.writeto('copilot2.fits')
    
TC = fits_table()
for c in C.columns():
    TC.set(c, [])
        
# First: mosstat vs copilot for im16.
for i in range(len(T)):
    J = np.flatnonzero((C.expnum == T.expnum[i]) * (C.extension == 'im16'))
    assert(len(J) == 1)
    j = J[0]
    for c in C.columns():
        TC.get(c).append(C.get(c)[j])
TC.to_np_arrays()
T.add_columns_from(TC)

import pylab as plt

plt.clf()
# plt.subplot(1,2,1)
# plt.plot(T.dx, T.ra_offset, 'b.')
# plt.xlabel('dx')
# plt.ylabel('RA offset')
# plt.subplot(1,2,2)
# plt.plot(T.dy, T.dec_offset, 'b.')
# plt.xlabel('dy')
# plt.ylabel('Dec offset')
p1 = plt.plot(T.ra_offset,  T.dy * 0.262, 'b.')
p2 = plt.plot(T.dec_offset, T.dx * 0.262, 'r.')
plt.legend([p1[0],p2[0]], ['dy,dRA','dx,dDec'], 'upper left')
plt.xlabel('dDec | dRA')
plt.ylabel('dx | dy')
plt.savefig('diff.png')


dra = np.median(T.ra_offset - T.dra)
ddec = np.median(T.dec_offset - T.ddec)
print('Shift in im16 dRA,dDec, arcsec: %.2f, %.2f' % (dra, ddec))

plt.clf()
p1 = plt.plot(T.ra_offset,  T.dra, 'b.')
p2 = plt.plot(T.dec_offset, T.ddec, 'r.')
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
p1 = plt.plot(T.expnum, T.dx * 0.262, 'b.')
p2 = plt.plot(T.expnum, T.dy * 0.262, 'r.')
p3 = plt.plot(T.expnum, T.ra_offset, 'g.')
p4 = plt.plot(T.expnum, T.dec_offset, 'm.')
plt.legend([p1[0],p2[0],p3[0],p4[0]], ['dx','dy','dRA','dDec'], 'lower right')
plt.xlabel('expnum')
plt.savefig('difft.png')



# Now, look at each copilot extension vs im16.
from astrometry.util.plotutils import PlotSequence
ps = PlotSequence('astrom')

Cref = C[C.extension == 'im16']
print(len(Cref), 'im16 exposures')
ref_expnum = dict([(expnum,i) for i,expnum in enumerate(Cref.expnum)])

CICR = []
for ext in range(1, 16):
    Ci = C[C.extension == 'im%i' % ext]
    print('Extension', ext, ':', len(Ci), 'exposures')
    Cr = Cref[np.array([ref_expnum[expnum] for expnum in Ci.expnum])]

    CICR.append((ext,Ci,Cr))
    
    drai = np.median(Cr.dra - Ci.dra)
    ddeci = np.median(Cr.ddec - Ci.ddec)
    
    plt.clf()
    plt.plot(Cr.dra, Ci.dra, 'b.')
    plt.plot(Cr.ddec, Ci.ddec, 'r.')
    ax = plt.axis()
    xx = np.array([-20,20])
    p1 = plt.plot(xx, xx - drai, 'b-', alpha=0.2)
    p2 = plt.plot(xx, xx - ddeci, 'r-', alpha=0.2)
    plt.axis(ax)
    plt.legend([p1[0],p2[0]], ['dRA','dDec'], 'upper left')
    plt.xlabel('im16')
    plt.ylabel('im%i' % ext)
    plt.title('Mosaic3: Copilot offsets')
    ps.savefig()
    
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
    plt.xlabel('im16')
    plt.ylabel('im%i' % ext)
    plt.title('Mosaic3: Copilot offsets')
    ps.savefig()
    
