from __future__ import print_function
from obsbot import NominalCalibration

mosaic_nominal_pixscale = 0.262

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
        
    def zeropoint(self, band, ext=None):
        if ext is not None:
            try:
                return self.zp0[(band, ext)]
            except KeyError:
                pass
        return self.zp0[band]

    def sky(self, band):
        return self.sky0[band]

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
                k_co = 0.06,
                A_co = 1.263,
                )
        else:
            raise ValueError('Unknown band "%s"' % band)
        return fid
