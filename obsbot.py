from __future__ import print_function

import numpy as np

class NominalCalibration(object):
    '''
    Overridden or instantiated by Mosaic / DECam nominal calibrations.

    Attributes (must) include:

    - pixscale -- in arcsec/pixel
    
    '''

    def zeropoint(self, band, ext=None):
        pass

    def sky(self, band):
        pass

    def fiducial_exptime(self, band):
        '''
        Returns an object with attributes:

        - skybright
        - k_co
        - A_co
        
        '''
        pass

    

def exposure_factor(self, fid, cal,
                    airmass, ebv, seeing, skybright, transparency):
        '''
        Computes a factor by which the exposure time should be scaled
        relative to nominal.

        *band*: string, 'g','r' or 'z'
        *airmass*: airmass, float > 1
        *ebv*: extinction E(B-V) mags
        *seeing*: FWHM in arcsec
        *skybright*: sky brightness
        *transparency*: sky transparency
        '''

        r_half = 0.45 #arcsec
        ps = cal.pixscale

        def Neff(self, seeing):
            # magic 2.35: convert seeing FWHM into sigmas in arcsec.
            return (4. * np.pi * (seeing / 2.35)**2 +
                    8.91 * r_half**2 +
                    ps**2/12.)

        neff_fid = Neff(fid.seeing)
        neff     = Neff(seeing)

        # scaling = T_new/T_fid
    
        scaling = (1./transparency**2 *
                   10.**(0.8 * fid.k_co * (airmass - 1.)) *
                   10.**(0.8 * fid.A_co * ebv) *
                   (neff / neff_fid) *
                   10.**(-0.4 * (skybright - fid.skybright)))
        return scaling

