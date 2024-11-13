from osbot.measure_raw import RawMeasurer

class BokMeasurer(RawMeasurer):
    def __init__(self, *args, **kwargs):
        if not 'pixscale' in kwargs:
            import mosaic
            kwargs.update(pixscale = mosaic.mosaic_nominal_pixscale)
        super(BokMeasurer, self).__init__(*args, **kwargs)
        self.camera = '90prime'

    def get_band(self, primhdr):
        band = super(BokMeasurer,self).get_band(primhdr)
        # bokr -> r
        band = band.strip()
        band = band.replace('bok', '')
        print('Band', band)
        return band[0]

    def read_raw(self, F, ext):
        '''
        F: fitsio FITS object
        '''
        img = F[ext].read()
        hdr = F[ext].read_header()
        img = img.astype(np.float32)
        #print('img qts', np.percentile(img.ravel(), [0,25,50,75,100]))

        primhdr = F[0].read_header()

        # Subtract median overscan and multiply by gains 
        data = parse_section(hdr['DATASEC'], slices=True)
        bias = parse_section(hdr['BIASSEC'], slices=True)
        gain = primhdr['GAIN%i' % int(ext.replace('IM','').replace('im',''),10)]
        b = np.median(img[bias])
        #print('subtracting bias', b)
        img[data] = (img[data] - b) * gain
    
        # Trim the image
        trim = parse_section(hdr['TRIMSEC'], slices=True)
        # zero out all but the trim section
        trimg = img[trim].copy()
        img[:,:] = 0
        img[trim] = trimg
        return img,hdr

    def get_wcs(self, hdr):
        from astrometry.util.util import Tan
        wcs = Tan(hdr)
        return wcs

    def get_sky_and_sigma(self, img):
        # Spline sky model to handle (?) ghost / pupil?
        from tractor.splinesky import SplineSky
        splinesky = SplineSky.BlantonMethod(img, None, 256)
        skyimg = np.zeros_like(img)
        splinesky.addTo(skyimg)
        mnsky,sig1 = sensible_sigmaclip(img - skyimg)
        sky1 = np.median(skyimg)
        return skyimg,sky1,sig1

    def remove_sky_gradients(self, img):
        pass

    def colorterm_ps1_to_observed(self, ps1stars, band):
        return ps1_to_90prime(ps1stars, band)

    def get_ps1_band(self, band):
        return ps1cat.ps1band[band]


class BokCPMeasurer(BokMeasurer):
    def read_raw(self, F, ext):
        img = F[ext].read()
        hdr = F[ext].read_header()
        img = img.astype(np.float32)
        return img,hdr

    def get_wcs(self, hdr):
        from astrometry.util.util import wcs_pv2sip_hdr
        wcs = wcs_pv2sip_hdr(hdr)
        return wcs

    def zeropoint_for_exposure(self, band, ext=None, exptime=None, **kwargs):
        zp0 = super(BokCPMeasurer, self).zeropoint_for_exposure(band, ext=ext, exptime=exptime, **kwargs)
        if zp0 is None:
            return zp0
        zp0 -= 2.5 * np.log10(exptime)
        return zp0
