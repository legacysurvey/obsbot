from obsbot.measure_raw import RawMeasurer

class DesiCiMeasurer(RawMeasurer):
    def __init__(self, *args, **kwargs):
        ps = kwargs.pop('ps', None)
        verbose = kwargs.pop('verbose', None)
        super(DesiCiMeasurer, self).__init__(*args, **kwargs)
        self.camera = 'desici'
        self.pixscale = 0.130
        self.edge_trim = 0

    def get_band(self, primhdr):
        # HACK
        return 'r'

    def get_sky_and_sigma(self, img):
        sky,sig1 = sensible_sigmaclip(img[100:-100, 100:-100])
        print('Sky, sig1:', sky, sig1)
        sky1 = np.median(sky)
        return sky,sky1,sig1

    def update_astrometry(self, stars, wcs, fx, fy, trim_x0, trim_y0,
                          pixsc, img, hdr, fullW, fullH, ps, printmsg):
        import tempfile
        from astrometry.util.util import healpix_rangesearch_radec
        # Run Astrometry.net on cleaned-up image.
        tempimg = tempfile.NamedTemporaryFile(suffix='.fits')
        configfile = tempfile.NamedTemporaryFile(suffix='.cfg')
        fn = tempimg.name
        fitsio.write(fn, img, header=hdr, clobber=True)

        ra = hdr['CRVAL1']
        dec = hdr['CRVAL2']
        radius = 5.
        nside = 2
        healpixes = healpix_rangesearch_radec(ra, dec, radius, nside)

        configfn = configfile.name
        f = open(configfn, 'w')
        index_dir = os.environ.get('AN_INDEX_DIR', 'index-files')
        # export AN_INDEX_DIR=/global/project/projectdirs/cosmo/work/users/dstn/index-5000
        f.write('add_path %s\n' % index_dir +
                'inparallel\n' +
                '\n'.join(['index index-5001-%02i' % hp for hp in healpixes]) + '\n' +
                '\n'.join(['index index-5002-%02i' % hp for hp in healpixes]) + '\n')
        f.close()

        cmd = 'solve-field --config %s --continue' % configfn
        cmd += ' --ra %f --dec %f --radius 5' % (ra, dec)
        if self.ext == 'CIC':
            pixsc = 0.133
        else:
            pixsc = 0.123
            cmd += ' --xscale 1.09'
        cmd += ' --scale-low %f --scale-high %f --scale-units app' % (pixsc * 0.95, pixsc * 1.05)
        cmd += ' --objs 50 --crpix-center'
        cmd += ' --parity pos'
        cmd += ' --no-verify --new-fits none --solved none --match none --rdls none --corr none'
        #cmd += ' --plot-scale 0.25'
        cmd += ' --no-plots'
        cmd += ' --downsample 4'
        #cmd += ' -v'
        cmd += ' ' + fn
        print(cmd)
        os.system(cmd)

        from astrometry.util.util import Sip, Tan
        wcsfn = fn.replace('.fits', '.wcs')
        if os.path.exists(wcsfn):
            wcs2 = Sip(wcsfn)

            imh,imw = img.shape
            cx,cy = (imw + 1.) / 2., (imh + 1.) / 2.
            r,d = wcs2.pixelxy2radec(cx, cy)
            ok,x1,y1 = wcs.radec2pixelxy(r, d)
            dx = x1 - cx - trim_x0
            dy = y1 - cy - trim_y0
            # print('Computing dx,dy from initial WCS:')
            # print(wcs)
            # print('to')
            # print(wcs2)
            # print('Final WCS image center:', cx,cy)
            # print('-> RA,Dec', r,d)
            # print('-> initial WCS coords', x1,y1)
            # print('-> dx,dy', dx, dy)

            r1,d1 = wcs.radec_center()
            r2,d2 = wcs2.radec_center()
            dd = d2 - d1
            dr = (r2 - r1)*np.cos(np.deg2rad(d2))
            #print('dr,dd', dr*3600, dd*3600, 'arcsec')
            print('dx,dy (%.1f, %.1f), dr,dd (%.1f, %.1f)' % (dx, dy, dr*3600., dd*3600.))

            measargs = dict(dx=dx, dy=dy, dra=dr, ddec=dd)
            return wcs2, measargs
        print('Solving failed!  Trying histogram...')
        return self.histogram_astrometric_shift(stars, wcs, fx, fy, trim_x0, trim_y0, pixsc,
                                                img, hdr, fullW, fullH, ps, printmsg)

    def get_wcs(self, hdr):
        from astrometry.util.util import Tan
        # Needs extremely recent Astrometry.net
        wcs = Tan(hdr)
        return wcs

    def get_bias(self, ext):
        # From Aaron:
        # https://github.com/ameisner/ci_reduce/blob/master/py/ci_reduce/common.py#L315
        # result in in ADU
        bias_med_dict = {'CIE' : 991.0, 
                         'CIN' : 1013.0, 
                         'CIC' : 1020.0, 
                         'CIS' : 984.0, 
                         'CIW' : 990.0}
        return bias_med_dict.get(ext, 0.)

    def read_raw(self, F, ext):
        img = F[ext].read()
        hdr = F[ext].read_header()
        img = img.astype(np.float32)

        img -= self.get_bias(ext)

        # Estimate per-pixel noise via Blanton's 5-pixel MAD
        from scipy.ndimage.measurements import label, find_objects
        slice1 = (slice(0,-5,10),slice(0,-5,10))
        slice2 = (slice(5,None,10),slice(5,None,10))
        mad = np.median(np.abs(img[slice1] - img[slice2]).ravel())
        sig1 = 1.4826 * mad / np.sqrt(2.)
        #print('Computed sig1 by Blanton method:', sig1)
        med = np.median(img.ravel())
        #print('Median:', med)
        thresh = med + 5. * sig1
        blobs,nblobs = label(img > thresh)
        #print('Found', nblobs, 'blobs of pixels above threshold')
        nz = 0
        for sy,sx in find_objects(blobs):
            y0 = sy.start
            y1 = sy.stop
            if y1 > y0+1:
                continue
            x0 = sx.start
            x1 = sx.stop
            if x1 > x0+1:
                continue
            nz += 1
            img[y0,x0] = med
        print('Zeroed out', nz, 'pixels')
        return img, hdr

    def get_ps1_band(self, band):
        return ps1cat.ps1band[band]

    def colorterm_ps1_to_observed(self, ps1stars, band):
        return 0.

