from __future__ import print_function
from osbbot import NominalCalibration, ExposureTimeCalculator

class DecamExptime(object):
    def update(self, **kwargs):
        for k,v in kwargs.items():
            setattr(self, k, v)

class DecamNominalCalibration(object):
    '''
    '''
    def __init__(self):
        self.pixscale = 0.262

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

    def fiducial_exptime(self, band):
        fid = DecamExptime()

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

        fid.update(
            skybright = self.sky(band),
            seeing = 1.3,
            )
        
        return fid
