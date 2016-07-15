from __future__ import print_function
import os
import sys

from astrometry.util.fits import *
from astrometry.util.plotutils import *
from astrometry.util.starutil_numpy import mjdtodate

import pylab as plt
import numpy as np

from legacypipe.common import *

from decam import DecamNominalCalibration
import tractor

#target_mjd = 57445.5
#mjd_diff = 0.5
#mjd_diff = 7

ps = PlotSequence('bot')

nom = DecamNominalCalibration()

botfn = 'bot-matched.fits'
if not os.path.exists(botfn):
    obsfn = 'obsbot.fits'
    if not os.path.exists(obsfn):
        cmd = 'python copilot.py --fits %s' % obsfn
        print('Running:', cmd)
        os.system(cmd)
        
    bot = fits_table(obsfn)

    date = np.array(map(mjdtodate, bot.mjd_obs))

    # plt.clf()
    # for i,band in enumerate('grz'):
    #     I = np.flatnonzero(np.array([b.strip() == band for b in bot.band]))
    #     plt.plot(date[I], np.zeros(len(I))+i, '.', color=dict(z='m').get(band,band))
    # plt.xlabel('Date')
    # plt.ylabel('band')
    # xt,xl = plt.xticks()
    # print('Ticks', xt)
    # print('Labels', xl)
    # plt.xticks(xt, ['']*len(xt), rotation=90)
    # ml,mh = np.min(bot.mjd_obs[bot.mjd_obs > 0])-10, np.max(bot.mjd_obs)+10
    # xl,xh = mjdtodate(ml), mjdtodate(mh)
    # plt.xlim(xl,xh)
    # plt.twiny()
    # plt.xlim(ml,mh)
    # plt.xlabel('MJD')
    # plt.ylim(-0.5, 2.5)
    # ps.savefig()

    bot.cut(bot.mjd_obs != 0)
    
    #bot.cut(np.abs(bot.mjd_obs - target_mjd) < mjd_diff)
    #print(len(bot), 'images for MJD')
    bot.cut(np.argsort(bot.mjd_obs))
    for x in zip(bot.mjd_obs, bot.filename, bot.expnum):
        print(x)
    
    botccds = set(zip(bot.expnum, bot.extension))
    print(len(botccds), 'unique CCDs')
    
    survey = LegacySurveyData()
    #ccds = survey.get_ccds_readonly()
    ccds = survey.get_annotated_ccds()
    
    # HACK
    #ccds.cut(np.abs(ccds.mjd_obs - target_mjd) < mjd_diff)
    #print('Cut to', len(ccds), 'CCDs near target MJD')

    ccds.cut(ccds.ccdzpt < 99.)
    
    ccdmap = dict([((expnum,ext),i) for i,(expnum,ext) in enumerate(zip(ccds.expnum, ccds.ccdname))])

    #ccds.index = np.arange(len(ccds))
    imatched = []
    matched = []
    #for expnum,ccdname in botccds:
    for i,(expnum,ccdname,f) in enumerate(
            zip(bot.expnum, bot.extension, bot.band)):
        print('expnum', expnum, 'ccdname', ccdname, 'filter', f)
        try:
            iccd = ccdmap[(expnum, ccdname)]
        except KeyError:
            continue
        thisccd = ccds[np.array([iccd])]
        #thisccd = ccds[(ccds.expnum == expnum) * (ccds.ccdname == ccdname)]
        #print('    found', len(thisccd))
        #assert(len(thisccd) <= 1)
        matched.append(thisccd)
        #if len(thisccd) == 1:
        imatched.append(i)
    matched = merge_tables(matched)
    print('Found', len(matched), 'CCDs in survey CCDs file')
    imatched = np.array(imatched)
    print('Matched to', len(imatched), 'from bot file')
    bot.cut(imatched)
    bot.add_columns_from(matched)
    bot.writeto(botfn)

else:
    bot = fits_table(botfn)

print('Filters:', np.unique(bot.band))
print('Exptimes:', np.unique(bot.exptime))
print('Transparencies:', np.unique(bot.transparency)[:10])

bot.cut(bot.exptime > 30.)
# HACK
bot.cut(bot.transparency > 0.5)
print('CUT to', len(bot), 'images')

allbot = bot

# FROM OBSBOT:
r_half = 0.45 #arcsec
pixsc = 0.262 # pixscale
def Neff(seeing):
    # magic 2.35: convert seeing FWHM into sigmas in arcsec.
    return (4. * np.pi * (seeing / 2.35)**2 +
            8.91 * r_half**2 +
            pixsc**2/12.)

def Neff2(seeing):
    # magic 2.35: convert seeing FWHM into sigmas in arcsec.
    pfact = 1.15
    return ((4. * np.pi * (seeing / 2.35)**2)**(1./pfact) +
            (8.91 * r_half**2)**(1./pfact)) ** pfact

# Compute galaxy norm for a 1.3" Gaussian PSF, given FWHM seeing in arcsec
# and pixel scale in arcsec/pixel.
def get_galnorm(seeing, pixsc):
    S = 32
    W,H = S*2+1, S*2+1
    psf_sigma = seeing / 2.35 / pixsc
    tim = tractor.Image(data=np.zeros((H,W)), inverr=np.ones((H,W)),
                        psf=tractor.NCircularGaussianPSF([psf_sigma], [1.]),
                        wcs=tractor.NullWCS(pixscale=pixsc))
    gal = SimpleGalaxy(tractor.PixPos(S,S), tractor.Flux(1.))
    mm = Patch(0, 0, np.ones((H,W), bool))
    galmod = gal.getModelPatch(tim, modelMask=mm).patch
    galmod = np.maximum(0, galmod)
    galmod /= galmod.sum()
    return np.sqrt(np.sum(galmod**2))

pixsc = nom.pixscale
fidseeing = 1.3


if True:
    # Does averaging multiple CCDs per exposure reduce the scatter?
    for band in np.unique(bot.band):
        band = band.strip()
        bot = allbot[np.array([b.strip() == band for b in allbot.band])]
        print(len(bot), 'in', band, 'band')

        from collections import Counter
        c = Counter(bot.expnum)
        keep_expnum = set()
        for k,v in c.most_common():
            print('Exposure number', k, 'appears', v, 'times')
            if v <= 1:
                break
            keep_expnum.add(k)

        bot.cut(np.nonzero([expnum in keep_expnum
                            for expnum in bot.expnum])[0])
        print('Cut to', len(bot), 'with > 1 CCD per exposure')
        
        bad = np.flatnonzero((bot.photometric == False))
        print(len(bad), 'rows are non-photometric')
        
        tt = 'Obsbot vs Pipeline depth: %s band' % band
    
        fid = nom.fiducial_exptime(band)
        extinction = bot.ebv * fid.A_co
        plt.clf()
        notbad = np.flatnonzero((bot.photometric == True))
        equivtime = (bot.exptime / bot.expfactor)
        extdepth = bot.galdepth - extinction

        # 2-coverage target (90% fill)
        target_depth = dict(g=24.0, r=23.4, z=22.5)[band]
        # -> 1-coverage depth (- ~0.37 mag)
        target_depth -= 2.5*np.log10(np.sqrt(2.))

        p1 = plt.plot(equivtime[notbad], extdepth[notbad], 'b.', alpha=0.5)
        (xl,xh,yl,yh) = plt.axis()
        xl = 10.

        meantime = []
        meandepth = []
        for expnum in keep_expnum:
            I = np.flatnonzero(bot.expnum == expnum)
            print(len(I), 'CCDs for expnum', expnum)
            print('Exp factors:', bot.expfactor[I])
            meanf = np.mean(bot.expfactor[I])
            i = I[0]
            meantime.append(bot.exptime[i] / meanf)
            meandepth.append(np.mean(bot.galdepth[I] - extinction[I]))

        # plot fit line (x = f(y), just to confuse you...)
        # 1 mag = factor of 6.3 in exposure time
        slope = (10.**(1./2.5))**2
        diff = np.median(np.log(equivtime[notbad]) / np.log(slope) -
                         extdepth[notbad])
        yy = np.linspace(20, 25, 100)
        xx = slope**(yy + diff)
        plt.plot(xx, yy, 'k-', alpha=0.3)
        p2 = plt.plot(xx, yy+0.05, 'k--', alpha=0.3)
        plt.plot(xx, yy-0.05, 'k--', alpha=0.3)

        plt.ylabel('Pipeline galdepth (unextincted) (mag)')
        plt.xlabel('Bot-predicted equivalent exposure time (s)')

        plt.axhline(target_depth, color='b', alpha=0.3)
        plt.axvline(fid.exptime, color='b', alpha=0.3)
    
        plt.xscale('log')
        xt = [10, 20, 50, 100, 200, 300, 400, 500]
        plt.xticks(xt, ['%i'%t for t in xt])
        plt.axis([xl,xh,yl,yh])
        plt.title(tt)
        ps.savefig()

        plt.clf()
        plt.plot(meantime, meandepth, 'bo')

        plt.plot(xx, yy, 'k-', alpha=0.3)
        p2 = plt.plot(xx, yy+0.05, 'k--', alpha=0.3)
        plt.plot(xx, yy-0.05, 'k--', alpha=0.3)
        plt.ylabel('Pipeline galdepth (unextincted) (mag)')
        plt.xlabel('Bot-predicted equivalent exposure time (s)')
        plt.axhline(target_depth, color='b', alpha=0.3)
        plt.axvline(fid.exptime, color='b', alpha=0.3)
        plt.xscale('log')
        xt = [10, 20, 50, 100, 200, 300, 400, 500]
        plt.xticks(xt, ['%i'%t for t in xt])
        plt.axis([xl,xh,yl,yh])
        plt.title(tt)
        ps.savefig()

        # for expnum,mnt,mnd in zip(keep_expnum, meantime, meandepth):
        #     I = np.flatnonzero(bot.expnum == expnum)
        #     plt.plot([mnt + np.zeros(len(I)), equivtime[I]],
        #              [mnd + np.zeros(len(I)), extdepth[I]], 'b.-')
        #     
        # plt.axis([xl,xh,yl,yh])
        # ps.savefig()
        
        
    sys.exit(0)
    
bot = allbot




gauss_galnorm_nom = get_galnorm(fidseeing, pixsc)
# correct for ~17% difference between Gaussian galnorm and real PSF
galnorm_nom = gauss_galnorm_nom / 1.17

galneff_nom = Neff(fidseeing)

see = np.arange(0.7, 2.01, 0.1)
gauss_galnorms = np.array([get_galnorm(s, pixsc) for s in see])
galneffs = Neff(see)
galneff2_nom = Neff2(fidseeing)
galneffs2 = Neff2(see)

plt.clf()
p1 = plt.plot(see, (1. / gauss_galnorms**2) / (1. / gauss_galnorm_nom**2),
              'bo-')
p2 = plt.plot(see, galneffs / galneff_nom, 'g--')
p2 = plt.plot(see, galneffs2 / galneff2_nom, 'go-')
plt.legend((p1[0],p2[0]),('Galnorm (Gaussian PSF)', 'Gal Neff'), loc='upper left')
plt.xlabel('Seeing (arcsec)')
plt.ylabel('Exposure factor')
ps.savefig()




for band in np.unique(bot.band):
    band = band.strip()

    bot = allbot[np.array([b.strip() == band for b in allbot.band])]
    print(len(bot), 'in', band, 'band')

    bad = np.flatnonzero((bot.photometric == False))
    print(len(bad), 'rows are non-photometric')
    
    tt = 'Obsbot vs Pipeline depth: %s band' % band

    fid = nom.fiducial_exptime(band)
    extinction = bot.ebv * fid.A_co
    zp0 = nom.zeropoint(band)

    psfnorm_seeing = bot.pixscale_mean * 2.35 * (1. / (bot.psfnorm_mean * 2. * np.sqrt(np.pi)))
    
    factor = np.median(bot.seeing / psfnorm_seeing)
    xx = np.array([0, 5])
    
    plt.clf()
    plt.plot(psfnorm_seeing, bot.seeing, 'b.')
    plt.plot(psfnorm_seeing[bad], bot.seeing[bad], 'r.')
    plt.xlabel('Pipeline PSF norm -> seeing')
    plt.ylabel('Bot seeing')
    ax = plt.axis()
    plt.plot(xx, xx*factor, 'k-', alpha=0.3)
    p2 = plt.plot(xx, xx*factor*0.95, 'k--', alpha=0.3)
    plt.plot(xx, xx*factor*1.05, 'k--', alpha=0.3)
    plt.axis(ax)
    plt.title(tt)
    plt.legend([p2[0]], ['+- 5%'], loc='lower right')
    ps.savefig()

    # galnorm_seeing = bot.pixscale_mean * 2.35 * 1. / (bot.galnorm_mean * 2. * np.sqrt(np.pi))
    # 
    # A = np.zeros((len(bot),2))
    # A[:,0] = 1.
    # A[:,1] = galnorm_seeing
    # b = np.linalg.lstsq(A, bot.seeing)[0]
    # print('Lstsq:', b)
    # offset = b[0]
    # slope = b[1]
    # xx = np.array([0, 5])
    # 
    # plt.clf()
    # #plt.plot(bot.galnorm_mean, bot.seeing, 'b.')
    # plt.plot(galnorm_seeing, bot.seeing, 'b.')
    # plt.xlabel('Pipeline galaxy norm -> seeing')
    # plt.ylabel('Bot seeing')
    # ax = plt.axis()
    # p = plt.plot(xx, offset + xx*slope, 'k-', alpha=0.3)
    # plt.plot(xx, offset + xx*slope * 0.9, 'k--', alpha=0.3)
    # plt.plot(xx, offset + xx*slope * 1.1, 'k--', alpha=0.3)
    # plt.legend([p[0]], ['offset %0.2f, slope %0.3f' % (offset, slope)],
    #            loc='lower right')
    # plt.axis(ax)
    # plt.title(tt)
    # ps.savefig()


    galneff = Neff(bot.seeing)
    plotx = 1. / (bot.galnorm_mean)**2 * pixsc**2

    factor = np.median(galneff / plotx)
    xx = np.array([0, 1000])
    
    plt.clf()
    plt.plot(plotx, galneff, 'b.')
    plt.plot(plotx[bad], galneff[bad], 'r.')
    plt.xlabel('Pipeline Neff = 1 / galaxy norm^2 (arcsec^2)')
    plt.ylabel('Bot galaxy Neff (arcsec^2)')
    plt.title(tt)
    ax = plt.axis()
    p1 = plt.plot(xx, xx*factor, 'k-', alpha=0.3)
    plt.plot(xx, xx*factor*0.9, 'k--', alpha=0.3)
    p2 = plt.plot(xx, xx*factor*1.1, 'k--', alpha=0.3)
    plt.plot(xx, xx, 'r-', alpha=0.3)
    plt.axhline(galneff_nom, color='b', alpha=0.3)
    plt.axis(ax)
    plt.legend([p1[0], p2[0]], ['slope %0.3f' % (factor), '+- 10%'],
               loc='lower right')
    ps.savefig()

    # What does the Neff imply about depth factor?

    galnorm_factor = (1. / bot.galnorm_mean**2) / (1. / galnorm_nom**2)
    print('Nominal galnorm: Gaussian', gauss_galnorm_nom, 'corrected', galnorm_nom)
    print('Median galnorm:', np.median(bot.galnorm_mean))

    galneff_factor = galneff / galneff_nom
    print('Nominal Neff:', galneff_nom)
    print('Median Neff:', np.median(galneff))


    diff = np.median(bot.galdepth - bot.gaussgaldepth)
    print('Median different between galdepth and Gaussgaldepth:', diff)
    # plt.clf()
    # plt.plot(bot.gaussgaldepth, bot.galdepth, 'b.')
    # plt.xlabel('Gaussian galdepth')
    # plt.ylabel('Galdepth')
    # ax = plt.axis()
    # plt.plot(xx, xx, 'r-', alpha=0.3)
    # plt.plot(xx, xx+diff, 'k-', alpha=0.3)
    # plt.axis(ax)
    # plt.title(tt)
    # ps.savefig()

    # Duh, Gaussian Galnorm is just a factor ~ 17% larger than Galnorm.
    # galnorm_x = 1. / 10.**((bot.galdepth / -2.5) + 9)
    # gaussgalnorm_x = 1. / 10.**((bot.gaussgaldepth / -2.5) + 9)
    # factor = 10. ** (0.17 / 2.5)
    # plt.clf()
    # plt.plot(galnorm_x, gaussgalnorm_x, 'b.')
    # plt.xlabel('galnorm')
    # plt.ylabel('Gauss galnorm')
    # plt.title(tt)
    # ax = plt.axis()
    # plt.plot(xx, xx, 'r-', alpha=0.3)
    # plt.plot(xx, xx*factor, 'b-', alpha=0.3)
    # plt.axis(ax)
    # ps.savefig()
    
    factor = np.median(galneff_factor / galnorm_factor)
    xx = np.array([0, 10])

    scatter = np.std(galneff_factor / (galnorm_factor * factor))
    print('PSF Scatter:', scatter)
    
    plt.clf()
    p1 = plt.plot(galnorm_factor, galneff_factor, 'b.')
    plt.plot(galnorm_factor[bad], galneff_factor[bad], 'r.')
    plt.xlabel('Pipeline exposure factor from galnorm')
    plt.ylabel('Bot exposure factor from Neff')
    plt.title(tt)
    ax = plt.axis()
    p = plt.plot(xx, xx*factor, 'k-', alpha=0.3)
    p2 = plt.plot(xx, xx*factor*0.9, 'k--', alpha=0.3)
    plt.plot(xx, xx*factor*1.1, 'k--', alpha=0.3)
    plt.plot(xx, xx, 'r-', alpha=0.3)
    plt.legend([p[0],p2[0],p1[0]], ['slope %0.3f' % (factor), '+- 10%',
                                    'scatter %.1f %%' % (100.*scatter)],
               loc='lower right')
    plt.axis(ax)
    ps.savefig()

    
    skyflux = 10.**((bot.sky - zp0) / -2.5) * bot.exptime * pixsc**2

    xx = np.array([0, 30000])
    A = np.zeros((len(skyflux),2))
    A[:,0] = 1.
    A[:,1] = bot.avsky
    b = np.linalg.lstsq(A, skyflux)[0]
    print('Lstsq:', b)
    offset = b[0]
    slope = b[1]
    
    plt.clf()
    plt.plot(bot.avsky, skyflux, 'b.')
    plt.plot(bot.avsky[bad], skyflux[bad], 'r.')
    plt.xlabel('CP avsky')
    plt.ylabel('Bot sky flux')
    ax = plt.axis()
    p = plt.plot(xx, offset + xx*slope, 'k-', alpha=0.3)
    p2 = plt.plot(xx, offset + xx*slope*0.9, 'k--', alpha=0.3)
    plt.plot(xx, offset + xx*slope*1.1, 'k--', alpha=0.3)
    plt.axis(ax)
    plt.title(tt)
    plt.legend([p[0],p2[0]], ['offset %0.2f, slope %0.3f' % (offset, slope),
                              '+- 10%'],
               loc='lower right')
    ps.savefig()
    
    # Convert skyflux into a sig1 estimate
    # in Poisson process, mean = variance; total sky counts are distributed
    # like that.
    # ignore gain
    skyvar = skyflux
    skysig1 = np.sqrt(skyvar)
    # use (bot) zeropoint to scale to nanomaggy units
    zpscale = NanoMaggies.zeropointToScale(bot.zeropoint)
    skysig1 /= zpscale
    # extra factor of time from the zeropoint...
    skysig1 /= bot.exptime

    factor = np.median(skysig1 / bot.sig1)
    xx = np.array([0, 1])
    
    plt.clf()
    plt.plot(bot.sig1, skysig1, 'b.')
    plt.plot(bot.sig1[bad], skysig1[bad], 'r.')
    plt.xlabel('Pipeline sig1')
    plt.ylabel('Bot sig1 = f(sky, zpt)')
    ax = plt.axis()
    p = plt.plot(xx, xx*factor, 'k-', alpha=0.3)
    p2 = plt.plot(xx, (xx*factor)*0.95, 'k--', alpha=0.3)
    plt.plot(xx, (xx*factor)*1.05, 'k--', alpha=0.3)
    plt.plot(xx, xx, 'r-', alpha=0.3)
    plt.axis(ax)
    plt.legend([p[0],p2[0]], ['slope %0.3f' % (factor), '+- 5%'],
               loc='lower right')
    plt.title(tt)
    ps.savefig()

    # # sig1 includes effect of zeropoint...
    # # Larger noise, longer exptime.
    zpscale = NanoMaggies.zeropointToScale(bot.zeropoint)
    expfactor_sig1 = (bot.sig1 * zpscale)**2
    expfactor_sig1 *= bot.exptime / fid.exptime
    expfactor_sig1 /= np.median(expfactor_sig1)

    expfactor_sky = 10.**(-0.4 * (bot.sky - fid.skybright))

    factor = np.median(expfactor_sky / expfactor_sig1)
    xx = np.array([0, 10])

    scatter = np.std(expfactor_sky / (expfactor_sig1 * factor))
    print('Sky Scatter:', scatter)
    
    plt.clf()
    p1 = plt.plot(expfactor_sig1, expfactor_sky, 'b.')
    plt.plot(expfactor_sig1[bad], expfactor_sky[bad], 'r.')
    plt.xlabel('Pipeline exposure factor from sig1 (and zpt; arb. scale)')
    plt.ylabel('Bot exposure factor from sky')
    plt.title(tt)
    ax = plt.axis()
    p = plt.plot(xx, xx*factor, 'k-', alpha=0.3)
    p2 = plt.plot(xx, (xx*factor)*0.9, 'k--', alpha=0.3)
    plt.plot(xx, (xx*factor)*1.1, 'k--', alpha=0.3)
    plt.legend([p2[0], p1[0]], ['+- 10%', 'scatter %.1f %%' % (100.*scatter)],
               loc='lower right')
    plt.axis(ax)
    ps.savefig()

    
    diff = np.median(bot.zeropoint - bot.ccdzpt)
    xx = np.array([20,30])

    plt.clf()
    plt.plot(bot.ccdzpt, bot.zeropoint, 'b.')
    plt.plot(bot.ccdzpt[bad], bot.zeropoint[bad], 'r.')
    plt.xlabel('Pipeline zeropoint')
    plt.ylabel('Bot zeropoint')
    ax = plt.axis()
    plt.plot(xx, xx + diff, 'k-', alpha=0.3)
    p2 = plt.plot(xx, xx + diff+0.05, 'k--', alpha=0.3)
    plt.plot(xx, xx + diff-0.05, 'k--', alpha=0.3)
    plt.legend([p2[0]], ['+- 0.05 mag'], loc='lower right')
    plt.axis(ax)
    plt.title(tt)
    ps.savefig()

    bot_trans = 10.**(-0.4 * (zp0 - bot.zeropoint))
    expfactor_bot_zpt = 1./bot_trans**2
    p_trans = 10.**(-0.4 * (np.median(bot.ccdzpt) - bot.ccdzpt))
    expfactor_p_zpt = 1./p_trans**2

    factor = np.median(expfactor_bot_zpt / expfactor_p_zpt)
    xx = np.array([0, 10])

    scatter = np.std(expfactor_bot_zpt / (expfactor_p_zpt * factor))
    print('Zpt Scatter:', scatter)
    
    plt.clf()
    p1 = plt.plot(expfactor_p_zpt, expfactor_bot_zpt, 'b.')
    plt.plot(expfactor_p_zpt[bad], expfactor_bot_zpt[bad], 'r.')
    plt.xlabel('Pipeline exposure factor from zpt (arb. scale)')
    plt.ylabel('Bot exposure factor from zeropoint')
    plt.title(tt)
    ax = plt.axis()
    p = plt.plot(xx, xx*factor, 'k-', alpha=0.3)
    p2 = plt.plot(xx, (xx*factor)*0.9, 'k--', alpha=0.3)
    plt.plot(xx, (xx*factor)*1.1, 'k--', alpha=0.3)
    plt.legend([p2[0], p1[0]], ['+- 10%', 'scatter %.1f %%' % (100.*scatter)],
               loc='lower right')
    plt.axis(ax)
    ps.savefig()
    
    
    # unextinct both so that the plot range tells us achieved depth
    # in the quantity we care about
    fid = nom.fiducial_exptime(band)
    extinction = bot.ebv * fid.A_co

    # 2-coverage target (90% fill)
    target_depth = dict(g=24.0, r=23.4, z=22.5)[band]
    # -> 1-coverage depth (- ~0.37 mag)
    target_depth -= 2.5*np.log10(np.sqrt(2.))
    
    # Using pipeline galdepth estimate, compute expfactor.

    depthfactor = 10.**(-0.4 * (bot.galdepth - target_depth))
    # airmass factor and transparency factor are already included
    # in the zeropoint, hence in depthfactor.
    expfactor = (depthfactor**2 *
                 10.**(0.8 * fid.A_co * bot.ebv))
    # this is on top of the exposure time of this image / nominal
    expfactor *= bot.exptime / fid.exptime

    factor = np.median(bot.expfactor / expfactor)
    xx = np.array([0, 10])

    scatter = np.std(bot.expfactor / (expfactor * factor))
    print('Overall Scatter:', scatter)
    
    plt.clf()
    p1 = plt.plot(expfactor, bot.expfactor, 'b.')
    plt.plot(expfactor[bad], bot.expfactor[bad], 'r.')
    plt.xlabel('Pipeline expfactor from galdepth')
    plt.ylabel('Bot expfactor')
    ax = plt.axis()
    p = plt.plot(xx, xx*factor, 'k-', alpha=0.3)
    p2 = plt.plot(xx, (xx*factor)*0.9, 'k--', alpha=0.3)
    plt.plot(xx, (xx*factor)*1.1, 'k--', alpha=0.3)
    plt.plot(xx, xx, 'r-', alpha=0.3)
    plt.legend([p[0],p2[0],p1[0]], ['slope %0.3f' % (factor), '+- 10%',
                                    'scatter %.1f %%' % (100.*scatter)],
               loc='lower right')
    plt.axis(ax)
    plt.title(tt)
    ps.savefig()

    galnorm = 1. / np.sqrt(galneff)
    galsig1 = skysig1 / galnorm
    galdepth = -2.5 * (np.log10(5. * galsig1) - 9)

    diff = np.median(galdepth - bot.galdepth)
    xx = np.array([20, 25])
    plt.clf()
    plt.plot(bot.galdepth - extinction, galdepth - extinction, 'b.')
    plt.plot((bot.galdepth - extinction)[bad], (galdepth - extinction)[bad],
             'r.')
    plt.xlabel('Pipeline galdepth (unextincted)')
    plt.ylabel('Bot galdepth (unextincted)')
    plt.axvline(target_depth, color='b', alpha=0.3)
    ax = plt.axis()
    plt.plot(xx, xx+diff, 'k-', alpha=0.3)
    p2 = plt.plot(xx, xx+diff+0.05, 'k--', alpha=0.3)
    plt.plot(xx, xx+diff-0.05, 'k--', alpha=0.3)
    plt.axis(ax)
    plt.legend([p2[0]], ['+- 0.05 mag'], loc='lower right')
    plt.title(tt)
    ps.savefig()

    equivtime = (bot.exptime / bot.expfactor)

    extdepth = bot.galdepth - extinction
    
    plt.clf()

    notbad = np.flatnonzero((bot.photometric == True))

    p1 = plt.plot(equivtime[notbad], extdepth[notbad], 'b.', alpha=0.5)

    yl,yh = plt.ylim()
    xl,xh = np.min(equivtime[notbad]) * 0.9, np.max(equivtime[notbad]) * 1.1
    
    I = np.flatnonzero(reduce(np.logical_or,
                              [bot.exptime == fid.exptime,
                               bot.exptime == fid.exptime_min,
                               bot.exptime == fid.exptime_max]))
    plt.plot(equivtime[I], extdepth[I], 'g.')

    plt.plot(equivtime[bad], extdepth[bad], 'r.')
    plt.axis(ax)
    
    # plot fit line (x = f(y), just to confuse you...)
        # 1 mag = factor of 6.3 in exposure time
    slope = (10.**(1./2.5))**2
    diff = np.median(np.log(equivtime[notbad]) / np.log(slope) - extdepth[notbad])
    yy = np.linspace(20, 25, 100)
    xx = slope**(yy + diff)
    plt.plot(xx, yy, 'k-', alpha=0.3)
    p2 = plt.plot(xx, yy+0.05, 'k--', alpha=0.3)
    plt.plot(xx, yy-0.05, 'k--', alpha=0.3)

    plt.ylabel('Pipeline galdepth (unextincted) (mag)')
    plt.xlabel('Bot-predicted equivalent exposure time (s)')

    plt.axhline(target_depth, color='b', alpha=0.3)
    plt.axvline(fid.exptime, color='b', alpha=0.3)

    # max x limits symmetric around fiducial
    mid = fid.exptime
    maxrange = max(mid / xl, xh / mid)
    xl = mid / maxrange
    xh = mid * maxrange

    # y too
    mid = target_depth
    maxrange = max(mid - yl, yh - mid)
    yl = mid - maxrange
    yh = mid + maxrange*0.8

    # 
    targettime = slope**(target_depth + diff)
    print('Exposure time to hit nominal depth:', targettime)
    plt.axvline(targettime, color='k', alpha=0.1)
    plt.text(targettime, (yl+mid)/2, '%.1f s' % targettime, color='k', ha='left', va='center')
    
    plt.xscale('log')
    xt = [10, 20, 50, 100, 200, 300, 400, 500]
    plt.xticks(xt, ['%i'%t for t in xt])
    plt.axis([xl,xh,yl,yh])
    
    plt.twiny()
    bins = np.arange(20., 25., 0.05)
    n,b,p = plt.hist(extdepth[notbad], bins=bins, orientation='horizontal',
                     histtype='step', color='b')
    #plt.xlim(0, np.max(n)*4)
    plt.xlim(np.max(n)*4, 0)
    plt.ylim(yl,yh)
    plt.xticks([])
    
    plt.legend([p2[0]], ['+- 0.05 mag'], loc='lower right')
    plt.title(tt)
    ps.savefig()
