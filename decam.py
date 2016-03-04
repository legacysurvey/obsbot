from __future__ import print_function
from obsbot import NominalCalibration

decam_nominal_pixscale = 0.262

class DecamNominalCalibration(NominalCalibration):
    '''
    '''
    def __init__(self):
        self.pixscale = decam_nominal_pixscale
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

    def cdmatrix(self, ext):
        # listhead DECam_00488199.fits.fz | grep 'EXTNAME\|CD._.'
        science = (-1.824E-07,    7.2853E-05, -7.284990E-05,   -1.8218E-07)
        focus   = (       0.0, -7.5E-05, 7.5E-05,       0.0)
        if ext in ['FS1', 'FS2', 'FS3', 'FS4', 'FN1', 'FN2', 'FN3', 'FN4']:
            return focus
        return science

    def _fiducial_exptime(self, fid, band):
        if band == 'g':
            fid.update(
                k_co = 0.178,
                A_co = 3.303,
                )

        elif band == 'r':
            fid.update(
                k_co = 0.094,
                A_co = 2.285,
                )

        elif band == 'z':
            fid.update(
                k_co = 0.060,
                A_co = 1.263,
                )
        else:
            raise ValueError('Unknown band "%s"' % band)

        return fid

