from obsbot.measure_raw import RawMeasurer

class PointingCamMeasurer(RawMeasurer):

    # full image size: 3296 x 2472

    def __init__(self, *args, **kwargs):
        ps = kwargs.pop('ps', None)
        verbose = kwargs.pop('verbose', None)
        super(PointingCamMeasurer, self).__init__(*args, **kwargs)
        self.camera = 'pointing'

        self.subW = 3296 // 2
        self.subH = 2472 // 2

    def get_sky_and_sigma(self, img):
        #sky,sig1 = sensible_sigmaclip(img[1000:2000, 500:1500])
        # First quadrant
        sky,sig1 = sensible_sigmaclip(img[100:self.subH-100, 100:self.subW-100])
        print('Sky, sig1:', sky, sig1)
        sky1 = np.median(sky)
        return sky,sky1,sig1

    def get_wcs(self, hdr):
        from astrometry.util.util import Tan
        wcs = Tan(hdr)
        wcs = wcs.get_subimage(0, 0, self.subW, self.subH)
        return wcs

    def read_raw(self, F, ext):
        img = F[ext].read()
        hdr = F[ext].read_header()
        img = img.astype(np.float32)
        img = img[:self.subH, :self.subW]
        return img, hdr

    def get_reference_stars(self, wcs, band):
        # Use a cut version of the Gaia catalog.
        cat = GaiaCatalog()
        stars = cat.get_catalog_in_wcs(wcs)
        print('Got', len(stars), 'Gaia stars')
        if not 'topring_mag' in stars.get_columns():
            raise RuntimeError('Need a specially cut Gaia catalog for the pointing cam')
        stars.mag = stars.topring_mag
        stars.cut(stars.mag < 14)
        print('Cut at 14th mag:', len(stars))
        return stars

    def cut_reference_catalog(self, stars):
        return np.ones(len(stars), bool)

    def get_color_term(self, stars, band):
        return 0.

    def get_exptime(self, primhdr):
        exptime = primhdr.get('EXPTIME', 0.)
        if exptime == 0.:
            exptime = primhdr.get('EXPOSURE', 0.) / 1000.
        return exptime

