import os
import numpy as np
from scipy.ndimage.measurements import *
import fitsio
from collections import Counter

from astrometry.util.util import healpix_rangesearch_radec

def ci_astrom(fn):
    for ext in ['CIC', 'CIN', 'CIS', 'CIE', 'CIW']:
        img,hdr = fitsio.read(fn, ext=ext, header=True)
        print('img', img.dtype, img.shape)
        img = img.astype(np.float32)
        
        # Estimate per-pixel noise via Blanton's 5-pixel MAD
        slice1 = (slice(0,-5,10),slice(0,-5,10))
        slice2 = (slice(5,None,10),slice(5,None,10))
        mad = np.median(np.abs(img[slice1] - img[slice2]).ravel())
        sig1 = 1.4826 * mad / np.sqrt(2.)
        print('Computed sig1 by Blanton method:', sig1)
        med = np.median(img.ravel())
        print('Median:', med)

        thresh = med + 5. * sig1
        blobs,nblobs = label(img > thresh)
        print('Found', nblobs, 'blobs of pixels above threshold')
        # How many pixels in each blob?
        #npix = Counter(blobs.ravel())
        #print('Min label:', blobs.min())
        # zero_out = np.zeros(img.shape, bool)
        # for k,n in npix.items():
        #     if n == 1:
        #         remap[k] = 0

        nz = 0
        for sy,sx in find_objects(blobs):
            #print('slice', sy, sx)
            y0 = sy.start
            y1 = sy.stop
            if y1 > y0+1:
                continue
            x0 = sx.start
            x1 = sx.stop
            if x1 > x0+1:
                continue
            #print('Single pixel:', x0, y0)
            nz += 1
            img[y0,x0] = med
        print('Zeroed out', nz, 'pixels')

        #img = remap[img]

        filtfn = '/tmp/%s.fits' % ext
        fitsio.write(filtfn, img, header=hdr)

        ra = hdr['CRVAL1']
        dec = hdr['CRVAL2']
        radius = 5.
        nside = 2
        healpixes = healpix_rangesearch_radec(ra, dec, radius, nside)

        configfn = '/tmp/cfg'
        f = open(configfn, 'w')
        f.write('add_path /global/project/projectdirs/cosmo/work/users/dstn/index-5000\n' +
                'inparallel\n' +
                '\n'.join(['index index-5001-%02i' % hp for hp in healpixes]) + '\n' +
                '\n'.join(['index index-5002-%02i' % hp for hp in healpixes]) + '\n')
        f.close()
        #configfn = '/global/project/projectdirs/cosmo/work/users/dstn/index-5000/cfg'

        cmd = 'solve-field --config %s --scale-low 5 --scale-high 8 --scale-units amw --continue' % configfn
        cmd += ' --ra %f --dec %f --radius 5' % (ra, dec)
        if ext != 'CIC':
            cmd += ' --xscale 1.1'
        cmd += ' --objs 100'
        cmd += ' --plot-scale 0.25'
        cmd += ' --downsample 4'
        #cmd += ' -v'
        cmd += ' ' + filtfn
        print(cmd)
        os.system(cmd)


ci_astrom('/global/project/projectdirs/desi/spectro/data/20190414/00007000/ci-00007000.fits.fz')

