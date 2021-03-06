from __future__ import print_function
import os
import fitsio
from astrometry.blind.plotstuff import *
from astrometry.util.util import *
from astrometry.util.fits import fits_table
from astrometry.libkd.spherematch import match_radec
from astrometry.util.starutil import hmsstring2ra, dmsstring2dec
# from camera_mosaic import ephem_observer
# import ephem
# 
# obs = ephem_observer()

tycho = fits_table('/data1/tycho2.fits.gz')
tycho.cut(tycho.mag < 8)
print(len(tycho), 'Tycho-2 stars < mag 8')

flipcd = True

#for expnum in (list(range(136934, 136942+1)) +
#               list(range(136949, 136951+1))):
#for expnum in [137031, 137027, 137036, 137039, 137040]:

#for expnum in list(range(140881, 140923+1)):
for expnum in list(range(140919, 140923+1)):

    #outfn = 'img-%i.png' % expnum
    #if os.path.exists(outfn):
    #    continue

    origwcs = []
    truewcs = []

    #fn = '/mosaic3/data3/observer/20170820/mos3_%i.fits' % expnum
    #fn = '/mosaic3/data3/observer/20170828/mos3_%i.fits' % expnum
    fn = '/mosaic3/data3/observer/20170925/mos3_%i.fits' % expnum

    hdr = fitsio.read_header(fn)
    # date = hdr['DATE-OBS'] # UTC
    # date = ephem.Date(date.replace('-','/').replace('T',' '))
    # obs.date = date

    ra = hdr['RA']
    dec = hdr['DEC']
    ra = hmsstring2ra(ra)
    dec = dmsstring2dec(dec)

    print()
    print('-----------------------------')
    print('Exposure number:', expnum)
    print('RA,Dec', ra, dec)
    print('Altitude:', hdr['ALT'])

    for ext in range(1, 17):

        print('Reading header from', fn, 'ext', ext)
        hdr = fitsio.read_header(fn, ext=ext)

        if flipcd:
            for key in ['CD1_1', 'CD1_2', 'CD2_1', 'CD2_2']:
                val = hdr[key]
                hdr[key] = -val

        #print('Header:', hdr)
        owcs = wcs_pv2sip_hdr(hdr)

        hdr2 = fitsio.FITSHDR()
        owcs.add_to_header(hdr2)
        fitsio.write('/tmp/orig-%i-%i.wcs.fits' % (expnum, ext), np.zeros((1,1)), header=hdr2, clobber=True)

        origwcs.append(owcs)
        rc,dc = owcs.get_crval()

        print('CRVAL:', rc,dc)
        hp = radecdegtohealpix(rc, dc, 2)
        print('Healpix:', hp)

        wcsfn = '/tmp/mos3-%i-%i.wcs' % (expnum, ext)
        if not os.path.exists(wcsfn):
            cmd = 'solve-field %s --config /data1/astrometry-index/cfg --dir /tmp --out mos3-%i-%i --crpix-center --ext %i --continue --no-plots -N none --no-tweak --objs 40 --scale-low 0.25 --scale-high 0.27 --scale-units app --ra %.3f --dec %.3f --radius 5.' % (fn, expnum, ext, ext, rc, dc)
            print(cmd)
            os.system(cmd)
        if os.path.exists(wcsfn):
            print('Reading WCS from', wcsfn)
            wcs = Tan(wcsfn, 0)
            truewcs.append((ext,wcs))

            rc1,dc1 = owcs.radec_center()
            rc2,dc2 = wcs.radec_center()

            dra = (rc2 - rc1) * np.cos(np.deg2rad(dc1))
            ddec = (dc2 - dc1)
            print('RA,Dec offset: (%.3f, %.3f) deg' % (dra, ddec))
            print('RA,Dec offset: (%.1f, %.1f) arcsec' % (dra*3600., ddec*3600.))
            ok,x2,y2 = wcs.radec2pixelxy(rc2,dc2)
            ok,x1,y1 = wcs.radec2pixelxy(rc1,dc1)
            #print('x1,y1', x1,y1)
            #print('x2,y2', x2,y2)
            print('X,Y offset: (%.1f, %.1f)' % ((x2-x1), (y2-y1)))

    I,J,d = match_radec(tycho.ra, tycho.dec, rc, dc, 1.)
    print(len(I), 'Tycho-2 stars nearby')

    plot = Plotstuff(outformat='png', size=(800,800), rdw=(rc,dc,1))

    plot.color = 'verydarkblue'
    plot.plot('fill')
    plot.color = 'green'
    plot.lw = 5
    for wcs in origwcs:
        plot.outline.wcs = anwcs_new_sip(wcs)
        plot.plot('outline')
        plot.color = 'white'
        plot.lw = 1

    plot.color = 'orange'
    plot.lw = 5
    for e,wcs in truewcs:
        plot.outline.wcs = anwcs_new_tan(wcs)
        plot.plot('outline')
        plot.color = 'yellow'
        plot.lw = 1
    plot.color = 'gray'
    plot.plot_grid(1., 1., 1., 1.)

    plot.color = 'white'
    for r,d,m in zip(tycho.ra[I], tycho.dec[I], tycho.mag[I]):
        plot.set_markersize(15. * np.clip(np.exp(6 - m), 0.05, 1.))
        plot.apply_settings()
        plot.marker_radec(r,d)
    plot.fill()

    outfn = 'wcs-%i.png' % expnum
    plot.write(outfn)
    print('Wrote', outfn)

    if True:
        for ext,wcs in truewcs:
            img = fitsio.read(fn, ext=ext)
            lo,hi = np.percentile(img.ravel(), [25, 98])
            rgb = (np.clip((img - lo) / (hi - lo), 0., 1.) * 255.).astype(np.uint8)
            #rgb = rgb[np.newaxis,:,:].repeat(3, axis=0)
            rgb = rgb[:,:,np.newaxis].repeat(3, axis=2)
            plot.image.set_image_from_numpy(rgb)
            plot.image.wcs = anwcs_new_tan(wcs)
            plot.plot('image')

        outfn = 'img-%i.png' % expnum
        plot.write(outfn)
        print('Wrote', outfn)
