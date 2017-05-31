from __future__ import print_function
from astrometry.util.fits import fits_table
import pylab as plt
import numpy as np
import fitsio
#import time
from astrometry.util.util import Tan, wcs_pv2sip_hdr
from astrometry.libkd.spherematch import match_radec
from astrometry.util.resample import resample_with_wcs, OverlapError
from astrometry.util.plotutils import PlotSequence
#from astrometry.util.starutil_numpy import *

from legacypipe.survey import LegacySurveyData
from collections import Counter

'''
A script to look into our current depth vs fill factor, to ask, "could
we retire any planned tiles", by, eg, making depth vs fill factor plots.

'''

if __name__ == '__main__':
    ps = PlotSequence('covfill2')

    T = fits_table('all-depths-1.fits')
    print('Depths for', len(T), 'tiles')
    depths = T.depths

    print('Depths:', depths.shape)
    
    target = 22.5

    I = np.flatnonzero((depths[:,2] >= target - 0.6) *
                       (depths[:,5] >= target - 0.3) *
                       (depths[:,10] >= target))
    print(len(I), 'tiles are retirable')

    I = np.flatnonzero((depths[:,5] >= target - 0.3) *
                       (depths[:,10] >= target))
    print(len(I), 'tiles meet 5% and 10% targets')

    I = np.flatnonzero((depths[:,10] >= target))
    print(len(I), 'tiles meet 10% target')

    I = np.flatnonzero((depths[:,3] >= target - 0.6) *
                       (depths[:,5] >= target - 0.3) *
                       (depths[:,10] >= target))
    print(len(I), 'tiles meet 2->3%, 5%, and 10% targets')

    I = np.flatnonzero((depths[:,2] >= target - 0.6) *
                       (depths[:,6] >= target - 0.3) *
                       (depths[:,10] >= target))
    print(len(I), 'tiles meet 2%, 5->6%, and 10% targets')

    I = np.flatnonzero((depths[:,2] >= target - 0.6) *
                       (depths[:,5] >= target - 0.3) *
                       (depths[:,11] >= target))
    print(len(I), 'tiles meet 2%, 5%, and 10->11% targets')

    I = np.flatnonzero((depths[:,2] >= target - 0.6 - 0.01) *
                       (depths[:,5] >= target - 0.3 - 0.01) *
                       (depths[:,10] >= target - 0.01))
    print(len(I), 'tiles meet 2%, 5%, and 10-% targets - 0.01 mag')

    I = np.flatnonzero((depths[:,2] >= target - 0.6 - 0.02) *
                       (depths[:,5] >= target - 0.3 - 0.02) *
                       (depths[:,10] >= target - 0.02))
    print(len(I), 'tiles meet 2%, 5%, and 10-% targets - 0.02 mag')

