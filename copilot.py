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

import fitsio
import ephem

from astrometry.util.plotutils import PlotSequence
from astrometry.util.starutil_numpy import hmsstring2ra, dmsstring2dec

from nightlystrategy import ExposureFactor, getParserAndGlobals, setupGlobals

from measure_raw_decam import measure_raw_decam

from tractor.sfd import SFDMap

def process_image(fn, ext, gvs, sfd, opt):
    portion = opt.portion

    # Read primary FITS header
    phdr = fitsio.read_header(fn)
    # Write QA plots to files named by the exposure number
    expnum = phdr['EXPNUM']
    ps = PlotSequence('qa-%i' % expnum)
    # Measure the new image
    M = measure_raw_decam(fn, ext=ext, ps=ps)

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
    print('Exposure factor:', expfactor)
    exptime = expfactor * gvs.base_exptimes[band]
    print('Exptime (un-clipped)', exptime)
    exptime = np.clip(exptime, gvs.floor_exptimes[band], gvs.ceil_exptimes[band])
    print('Clipped exptime', exptime)
    
    if band == 'z' and exptime > gvs.t_sat_max:
        exptime = gvs.t_sat_max
        print('Reduced exposure time to avoid z-band saturation:', exptime)

    print

    print('Actual exposure time taken:', actual_exptime)

    print('Depth fraction: %6.3f' % (actual_exptime / exptime))
    
    # If you were going to re-plan, you would run with these args:
    plandict = dict(seeing=M['seeing'], transparency=M['transparency'])
    # Assume the sky is as much brighter than canonical in each band... unlikely
    dsky = M['skybright'] - gvs.sb_dict[M['band']]
    for band in 'grz':
        plandict['sb'+band] = gvs.sb_dict[band] + dsky
    # Note that nightlystrategy.py takes UTC dates.
    now = datetime.datetime.utcnow()
    d = now.date()
    plandict['startdate'] = '%04i-%02i-%02i' % (d.year, d.month, d.day)
    t = now.time()
    plandict['starttime'] = t.strftime('%H:%M:%S')
    # Make an hour-long plan
    end = now + datetime.timedelta(0, 3600)
    d = end.date()
    plandict['enddate'] = '%04i-%02i-%02i' % (d.year, d.month, d.day)
    t = end.time()
    plandict['endtime'] = t.strftime('%H:%M:%S')

    # Decide the pass.
    goodseeing = plandict['seeing']<1.3
    photometric = plandict['transparency']>0.95

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
    print('python2.7 nightlystrategy.py --start-twi 10 --end-twi 10 --seeg %(seeing).3f --seer %(seeing).3f --seez %(seeing).3f --sbg %(sbg).3f --sbr %(sbr).3f --sbz %(sbz).3f --transparency %(transparency).3f --start-date %(startdate)s --start-time %(starttime)s --end-date %(enddate)s --end-time %(endtime)s --date %(startdate)s --portion %(portion)f --pass %(pass)i' % plandict) 
    print()

    return M, plandict, expnum
    
    
if __name__ == '__main__':
    
    parser = optparse.OptionParser(usage='%prog')
    
    parser.add_option('--ext', help='Extension to read for computing observing conditions: default %default', default='N4')
    parser.add_option('--rawdata', help='Directory to monitor for new images: default %default', default='rawdata')
    parser.add_option('--portion', help='Portion of the night: default %default', type=float, default='1.0')

    opt,args = parser.parse_args()

    imagedir = opt.rawdata
    rawext = opt.ext

    # Get nightlystrategy data structures; use fake command-line args.
    # these don't matter at all, since we only use ExposureFactor
    parser,gvs = getParserAndGlobals()
    nsopt,nsargs = parser.parse_args('--date 2015-01-01 --pass 1 --portion 1'.split())
    obs = setupGlobals(nsopt, gvs)

    print('Loading SFD maps...')
    sfd = SFDMap()
    
    if len(args) > 0:
        for fn in args:
            print('Reading', fn)
            process_image(fn, rawext, gvs, sfd, opt)
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

        process_image(fn, rawext, gvs, sfd, opt)

