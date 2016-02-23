from __future__ import print_function
from obsbot import NominalCalibration

class DecamNominalCalibration(NominalCalibration):
    '''
    '''
    def __init__(self):
        self.pixscale = 0.262
        self.overhead = 30

        self.zp0 = dict(
            g = 26.610,
            r = 26.818,
            z = 26.484,
            )

        self.sky0 = dict(
            g = 22.04,
            r = 20.91,
            z = 18.46,
            )
        
    def zeropoint(self, band, ext=None):
        return self.zp0[band]

    def sky(self, band):
        return self.sky0[band]

    def _fiducial_exptime(self, fid, band):
        if band == 'g':
            fid.update(
                k_co = 0.17,
                A_co = 3.303,
                )

        elif band == 'r':
            fid.update(
                k_co = 0.1,
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

