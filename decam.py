from __future__ import print_function
from obsbot import NominalCalibration

decam_nominal_pixscale = 0.262

class DecamNominalCalibration(NominalCalibration):
    '''
    '''
    def __init__(self):
        self.pixscale = decam_nominal_pixscale
        self.overhead = 30
        self.gain = 4.0
        self.saturation_adu = 15000

        self.zp0 = dict(
            u = 22.46,
            g = 26.610,
            r = 26.818,
            #i = 26.86, adjusted based on kentool clouds measure
            i = 26.95,
            z = 26.484,
            #N419 = 23.259,
            N419 = 22.68,
            ## arbitrary set from exp 961101,102
            N501 = 23.812,
            N540 = 25.32,
            #            N540 = 24.92,
            N673 = 24.151,
            N708 = 24.92,
            # made up, = g
            M411 = 26.610,
            M464 = 26.610,
            )

        self.sky0 = dict(
            u = 22.04,
            g = 22.04,
            r = 20.91,
            # i just set as the average of r,z
            i = 19.69,
            z = 18.46,
            ##
            N419 = 19.128,
            N501 = 19.128,
            N540 = 19.614,
            N673 = 19.614,
            N708 = 19.614,
            # made up, = g
            M411 = 22.04,
            M464 = 22.04,
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
        if band == 'u':
            fid.update(
                # COPIED FROM 'g'!
                k_co = 0.17,
                A_co = 3.214,
                )
        elif band == 'g':
            fid.update(
                k_co = 0.17,
                A_co = 3.214,
                )

        elif band == 'r':
            fid.update(
                k_co = 0.10,
                A_co = 2.165,
                )

        elif band == 'i':
            # just set to the average of r,z
            fid.update(
                k_co = 0.080,
                A_co = 1.878,
                )

        elif band == 'z':
            fid.update(
                k_co = 0.060,
                A_co = 1.592,
                )

        elif band == 'N419':
            fid.update(
                k_co = 0.17,
                A_co = 3.214,
                )

        elif band == 'N501':
            fid.update(
                k_co = 0.17,
                A_co = 3.214,
                )

        elif band == 'N540':
            fid.update(
                k_co = 0.17,
                A_co = 3.214,
                )

        elif band == 'N673':
            fid.update(
                k_co = 0.10,
                A_co = 2.165,
                )

        elif band == 'N708':
            fid.update(
                k_co = 0.10,
                A_co = 2.165,
                )

        elif band == 'M411':
            # made up, just = g band
            fid.update(
                k_co = 0.17,
                A_co = 3.214,
                )
        elif band == 'M464':
            # made up, just = g band
            fid.update(
                k_co = 0.17,
                A_co = 3.214,
                )

        else:
            raise ValueError('Unknown band "%s"' % band)

        return fid
