from obsbot.measure_raw import RawMeasurer

# Color terms -- no MOSAIC specific ones yet:
ps1_to_mosaic = ps1_to_decam

def measure_raw_mosaic3(fn, ext='im4', nom=None, ps=None,
                        measargs={}, **kwargs):
    if nom is None:
        import mosaic
        nom = mosaic.MosaicNominalCalibration()
    meas = Mosaic3Measurer(fn, ext, nom, **measargs)
    results = meas.run(ps, **kwargs)
    return results

class Mosaic3Measurer(RawMeasurer):
    def __init__(self, *args, **kwargs):
        if not 'pixscale' in kwargs:
            import mosaic
            kwargs.update(pixscale = mosaic.mosaic_nominal_pixscale)
        super(Mosaic3Measurer, self).__init__(*args, **kwargs)
        self.camera = 'mosaic3'

    def get_band(self, primhdr):
        band = super(Mosaic3Measurer,self).get_band(primhdr)
        # "zd" -> "z"
        #return band[0]
        return band

    def read_raw(self, F, ext):
        '''
        F: fitsio FITS object
        '''
        img = F[ext].read()
        hdr = F[ext].read_header()
        img = img.astype(np.float32)
        #print('img qts', np.percentile(img.ravel(), [0,25,50,75,100]))
        # Subtract median overscan and multiply by gains 
        dataA = parse_section(hdr['DATASEC'], slices=True)
        biasA = parse_section(hdr['BIASSEC'], slices=True)
        gainA = hdr['GAIN']
        b = np.median(img[biasA])
        #print('subtracting bias', b)
        img[dataA] = (img[dataA] - b) * gainA
    
        # Trim the image
        trimA = parse_section(hdr['TRIMSEC'], slices=True)
        # zero out all but the trim section
        trimg = img[trimA].copy()
        img[:,:] = 0
        img[trimA] = trimg
        return img,hdr

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

    def get_wcs(self, hdr):
        # Older images have ZPX, newer TPV.
        if hdr['CTYPE1'] == 'RA---TPV':
            from astrometry.util.util import wcs_pv2sip_hdr
            wcs = wcs_pv2sip_hdr(hdr)
        else:
            from astrometry.util.util import Tan
            hdr['CTYPE1'] = 'RA---TAN'
            hdr['CTYPE2'] = 'DEC--TAN'
            wcs = Tan(hdr)
        return wcs

    def colorterm_ps1_to_observed(self, ps1stars, band):
        bandmap = dict(zd='z', D51='g')
        band = bandmap.get(band, band)
        return ps1_to_mosaic(ps1stars, band)

    def get_ps1_band(self, band):
        bandmap = dict(zd='z', D51='g')
        band = bandmap.get(band, band)
        return ps1cat.ps1band[band]
