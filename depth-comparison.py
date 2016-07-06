from __future__ import print_function

from astrometry.util.fits import *
from astrometry.util.plotutils import *

import pylab as plt
import numpy as np

from legacypipe.common import *

target_mjd = 57445.5
#mjd_diff = 0.5
mjd_diff = 7


botfn = 'bot-matched.fits'
if not os.path.exists(botfn):
    bot = fits_table('obsbot.fits')
    bot.cut(np.abs(bot.mjd_obs - target_mjd) < mjd_diff)
    print(len(bot), 'images for MJD')
    bot.cut(np.argsort(bot.mjd_obs))
    for x in zip(bot.mjd_obs, bot.filename, bot.expnum):
        print(x)
    
    botccds = set(zip(bot.expnum, bot.extension))
    print(len(botccds), 'unique CCDs')
    
    survey = LegacySurveyData()
    #ccds = survey.get_ccds_readonly()
    ccds = survey.get_annotated_ccds()
    
    # HACK
    ccds.cut(np.abs(ccds.mjd_obs - target_mjd) < mjd_diff)
    print('Cut to', len(ccds), 'CCDs near target MJD')
    
    #ccds.index = np.arange(len(ccds))
    imatched = []
    matched = []
    #for expnum,ccdname in botccds:
    for i,(expnum,ccdname,f) in enumerate(
            zip(bot.expnum, bot.extension, bot.band)):
        print('expnum', expnum, 'ccdname', ccdname, 'filter', f)
        thisccd = ccds[(ccds.expnum == expnum) * (ccds.ccdname == ccdname)]
        print('    found', len(thisccd))
        assert(len(thisccd) <= 1)
        matched.append(thisccd)
        if len(thisccd) == 1:
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

ps = PlotSequence('bot')

bot.cut(bot.exptime > 30.)

# HACK
bot.cut(bot.transparency > 0.5)
print('CUT to', len(bot), 'images')

allbot = bot

for band in np.unique(bot.band):
    band = band.strip()

    bot = allbot[np.array([b.strip() == band for b in allbot.band])]
    print(len(bot), 'in', band, 'band')
    
    tt = 'Obsbot depth vs Pipeline depth: 2016-02-25 (%s)' % band

    plt.clf()
    plt.plot(bot.galdepth, bot.expfactor, 'b.')
    plt.xlabel('galdepth')
    plt.ylabel('bot expfactor')
    plt.title(tt)
    ps.savefig()
    
    # plt.clf()
    # plt.plot(bot.gaussgaldepth, bot.expfactor, 'b.')
    # plt.xlabel('gaussian galdepth')
    # plt.ylabel('bot expfactor')
    # ps.savefig()
    
    expfactor_depth = 1. / np.sqrt(bot.expfactor)
    
    diff = np.median(expfactor_depth - bot.galdepth)
    xx = np.array([20, 25])
    
    plt.clf()
    plt.plot(bot.galdepth, expfactor_depth, 'b.')
    plt.xlabel('Pipeline galdepth')
    plt.ylabel('Bot expfactor -> depth')
    ax = plt.axis()
    plt.plot(xx, xx+diff, 'k-', alpha=0.3)
    plt.plot(xx, xx+diff+0.1, 'k--', alpha=0.3)
    plt.plot(xx, xx+diff-0.1, 'k--', alpha=0.3)
    plt.axis(ax)
    plt.title(tt)
    ps.savefig()
    
    diff = np.median(expfactor_depth - bot.gaussgaldepth)
    
    # plt.clf()
    # plt.plot(bot.gaussgaldepth, expfactor_depth, 'b.')
    # plt.xlabel('gaussian galdepth')
    # plt.ylabel('bot expfactor -> depth')
    # ax = plt.axis()
    # plt.plot(xx, xx+diff, 'k-', alpha=0.3)
    # plt.axis(ax)
    # ps.savefig()
    
    psfnorm_seeing = bot.pixscale_mean * 2.35 * 1. / (bot.psfnorm_mean * 2. * np.sqrt(np.pi))
    
    factor = np.median(bot.seeing / psfnorm_seeing)
    xx = np.array([0, 5])
    
    plt.clf()
    plt.plot(psfnorm_seeing, bot.seeing, 'b.')
    plt.xlabel('Pipeline PSF norm -> seeing')
    plt.ylabel('Bot seeing')
    ax = plt.axis()
    plt.plot(xx, xx*factor, 'k-', alpha=0.3)
    plt.plot(xx, xx*factor*0.9, 'k--', alpha=0.3)
    plt.plot(xx, xx*factor*1.1, 'k--', alpha=0.3)
    plt.axis(ax)
    plt.title(tt)
    ps.savefig()

    galnorm_seeing = bot.pixscale_mean * 2.35 * 1. / (bot.galnorm_mean * 2. * np.sqrt(np.pi))

    #factor = np.median(bot.seeing / galnorm_seeing)
    A = np.zeros((len(bot),2))
    A[:,0] = 1.
    A[:,1] = galnorm_seeing
    b = np.linalg.lstsq(A, bot.seeing)[0]
    print('Lstsq:', b)
    offset = b[0]
    slope = b[1]
    xx = np.array([0, 5])
    
    plt.clf()
    #plt.plot(bot.galnorm_mean, bot.seeing, 'b.')
    plt.plot(galnorm_seeing, bot.seeing, 'b.')
    plt.xlabel('Pipeline galaxy norm -> seeing')
    plt.ylabel('Bot seeing')
    ax = plt.axis()
    p = plt.plot(xx, offset + xx*slope, 'k-', alpha=0.3)
    plt.plot(xx, offset + xx*slope * 0.9, 'k--', alpha=0.3)
    plt.plot(xx, offset + xx*slope * 1.1, 'k--', alpha=0.3)
    plt.legend([p[0]], ['offset %0.2f, slope %0.3f' % (offset, slope)],
               loc='lower right')
    plt.axis(ax)
    plt.title(tt)
    ps.savefig()


    # FROM OBSBOT:
    r_half = 0.45 #arcsec
    pixsc = 0.262 # pixscale
    def Neff(seeing):
        # magic 2.35: convert seeing FWHM into sigmas in arcsec.
        return (4. * np.pi * (seeing / 2.35)**2 +
                8.91 * r_half**2 +
                pixsc**2/12.)

    galneff = Neff(bot.seeing)

    plotx = 1. / (bot.galnorm_mean)**2
    
    factor = np.median(galneff / plotx)
    xx = np.array([0, 1000])
    
    plt.clf()
    plt.plot(plotx, galneff, 'b.')
    plt.xlabel('Pipeline 1 / galaxy norm^2')
    plt.ylabel('Bot galaxy Neff')
    plt.title(tt)
    ax = plt.axis()
    plt.plot(xx, xx*factor, 'k-', alpha=0.3)
    plt.plot(xx, xx*factor*0.9, 'k--', alpha=0.3)
    plt.plot(xx, xx*factor*1.1, 'k--', alpha=0.3)
    plt.axis(ax)
    ps.savefig()
    
    
    skyflux = 10.**(bot.sky / -2.5) * bot.exptime
    
    #factor = np.median(skyflux / bot.avsky)
    xx = np.array([0, 30000])
    A = np.zeros((len(skyflux),2))
    A[:,0] = 1.
    A[:,1] = bot.avsky
    b = np.linalg.lstsq(A, skyflux)[0]
    print('Lstsq:', b)
    offset = b[0]
    slope = b[1]
    
    plt.clf()
    #plt.plot(bot.avsky, bot.sky, 'b.')
    plt.plot(bot.avsky, skyflux, 'b.')
    #plt.scatter(bot.avsky, skyflux, c=bot.exptime)
    plt.xlabel('CP avsky')
    plt.ylabel('Bot sky flux')
    ax = plt.axis()
    plt.plot(xx, offset + xx*slope, 'k-', alpha=0.3)
    plt.plot(xx, offset + xx*slope*0.8, 'k--', alpha=0.3)
    plt.plot(xx, offset + xx*slope*1.2, 'k--', alpha=0.3)
    plt.axis(ax)
    plt.title(tt)
    ps.savefig()
    
    # # in Poisson process, mean = variance
    # sig1sky = bot.sig1**2
    # 
    # xx = np.array([0, 1])
    # A = np.zeros((len(skyflux),2))
    # A[:,0] = 1.
    # A[:,1] = sig1sky
    # b = np.linalg.lstsq(A, skyflux)[0]
    # print('Lstsq:', b)
    # offset = b[0]
    # slope = b[1]
    # 

    # sig1x = bot.sig1 * 10.**(bot.ccdzpt / 2.5)
    # plt.clf()
    # #plt.plot(sig1sky, skyflux, 'b.')
    # #plt.plot(bot.sig1, skyflux, 'b.')
    # plt.plot(sig1x, skyflux, 'b.')
    # #plt.xlabel('Pipeline sig1 -> sky')
    # plt.xlabel('Pipeline sig1 x zpt')
    # plt.ylabel('Bot sky flux')
    # # ax = plt.axis()
    # # plt.plot(xx, offset + xx*slope, 'k-', alpha=0.3)
    # # plt.plot(xx, offset + xx*slope*0.8, 'k--', alpha=0.3)
    # # plt.plot(xx, offset + xx*slope*1.2, 'k--', alpha=0.3)
    # # plt.axis(ax)
    # plt.title(tt)
    # ps.savefig()
    
    
    diff = np.median(bot.zeropoint - bot.ccdzpt)
    xx = np.array([20,30])
    
    plt.clf()
    plt.plot(bot.ccdzpt, bot.zeropoint, 'b.')
    plt.xlabel('Pipeline zeropoint')
    plt.ylabel('Bot zeropoint')
    ax = plt.axis()
    plt.plot(xx, xx + diff, 'k-', alpha=0.3)
    plt.plot(xx, xx + diff+0.1, 'k--', alpha=0.3)
    plt.plot(xx, xx + diff-0.1, 'k--', alpha=0.3)
    plt.axis(ax)
    plt.title(tt)
    ps.savefig()
