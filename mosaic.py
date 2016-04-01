from __future__ import print_function
from obsbot import NominalCalibration

mosaic_nominal_pixscale = 0.260

class MosaicNominalCalibration(NominalCalibration):
    '''
    '''
    def __init__(self):
        self.pixscale = mosaic_nominal_pixscale
        self.overhead = 30
        
        self.zp0 = dict(
            g = 26.930,
            r = 27.014,
            z = 26.518,
            )
        self.zp0.update({
            ('z', 'im4' ): 26.406,
            ('z', 'im7' ): 26.609,
            ('z', 'im11'): 26.556,
            ('z', 'im16'): 26.499,
            })

        self.sky0 = dict(
            g = 22.04,
            r = 20.91,
            z = 18.46,
            )

    def crpix(self, ext):
        # listhead mos3.68488.fits | grep '\(CRPIX\|EXTNAME\)'
        crpixes = dict(
            im1  = (-2161,  4209),
            im2  = (-  99,  4209),
            im3  = (-2161,  2191),
            im4  = (-  99,  2191),
            im5  = ( 2153,  4211),
            im6  = ( 4215,  4211),
            im7  = ( 2153,  2193),
            im8  = ( 4215,  2193),
            im9  = (-  89, -2106),
            im10 = (-2151, -2106),
            im11 = (-  89, -  88),
            im12 = (-2151, -  88),
            im13 = ( 4222, -2099),
            im14 = ( 2160, -2099),
            im15 = ( 4222, -  81),
            im16 = ( 2160, -  81))
        return crpixes[ext]
        
    def zeropoint(self, band, ext=None):
        if ext is not None:
            try:
                return self.zp0[(band, ext)]
            except KeyError:
                pass
        return self.zp0[band]

    def sky(self, band):
        return self.sky0[band]

    def cdmatrix(self, ext):
        # listhead mos3.68488.fits | grep 'EXTNAME\|CD._.'
        return dict(im1 =( 2.5e-07, 7.2e-05, 7.2e-05,-9.1e-08),
                    im2 =( 2.5e-07, 7.2e-05, 7.2e-05,-9.1e-08),
                    im3 =( 2.5e-07, 7.2e-05, 7.2e-05,-9.1e-08),
                    im4 =( 2.5e-07, 7.2e-05, 7.2e-05,-9.1e-08),
                    im5 =( 3.4e-08, 7.2e-05,  7.e-05,-1.9e-07),
                    im6 =( 3.4e-08, 7.2e-05,  7.e-05,-1.9e-07),
                    im7 =( 3.4e-08, 7.2e-05,  7.e-05,-1.9e-07),
                    im8 =( 3.4e-08, 7.2e-05,  7.e-05,-1.9e-07),
                    im9 =(-3.0e-08, 7.2e-05, 7.2e-05,-1.3e-07),
                    im10=(-3.0e-08, 7.2e-05, 7.2e-05,-1.3e-07),
                    im11=(-3.0e-08, 7.2e-05, 7.2e-05,-1.3e-07),
                    im12=(-3.0e-08, 7.2e-05, 7.2e-05,-1.3e-07),
                    im13=( 2.0e-07, 7.2e-05, 7.2e-05,-4.1e-08),
                    im14=( 2.0e-07, 7.2e-05, 7.2e-05,-4.1e-08),
                    im15=( 2.0e-07, 7.2e-05, 7.2e-05,-4.1e-08),
                    im16=( 2.0e-07, 7.2e-05, 7.2e-05,-4.1e-08))[ext]

    def _fiducial_exptime(self, fid, band):
        if band == 'g':
            fid.update(
                k_co = 0.211,
                A_co = 3.303,
                )

        elif band == 'r':
            fid.update(
                k_co = 0.109,
                A_co = 2.285,
                )

        elif band == 'z':
            fid.update(
                k_co = 0.10,
                A_co = 1.263,
                )
        else:
            raise ValueError('Unknown band "%s"' % band)
        return fid
