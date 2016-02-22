from __future__ import print_function
from osbbot import NominalCalibration, ExposureTimeCalculator

class MosaicExptime(object):
    def update(self, **kwargs):
        for k,v in kwargs.items():
            setattr(self, k, v)

class MosaicNominalCalibration(object):
    '''
    '''
    def __init__(self):
        self.pixscale = 0.26

        self.zp0 = dict(
            g = 26.930,
            r = 27.014,
            z = 26.455,
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
        '''
        Returns an object with attributes:

        - skybright
        - k_co
        - A_co
        - seeing
        
        '''
        fid = MosaicExptime()

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

        fid.update(
            skybright = self.sky(band),
            seeing = 1.3,
            )
        
        return fid
