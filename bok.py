from __future__ import print_function
from obsbot import NominalCalibration

bok_nominal_pixscale = 0.455

class BokNominalCalibration(NominalCalibration):
    '''
    '''
    def __init__(self):
        self.pixscale = bok_nominal_pixscale
        self.overhead = 30
        self.gain = 1.4
        self.saturation_adu = 30000

        self.zp0 = dict(
            # Arjun's values for CP-processed data
            g = 25.800,
            r = 25.580,
            )

        #Ian?s (based on raw):
        #bok_zpt0 = {'g':25.55,'bokr':25.38}

        # DECam
        # g = 26.610,
        # r = 26.818,
        # z = 26.484,

        # And these are the mosaic values...
        # 'g': begin
        # zpt0 = 26.93
        # sky0 = 22.04
        # kx = 0.17 ; atmospheric extinction coefficient
        # end
        # 'r': begin
        # zpt0 = 27.01
        # sky0 = 20.91
        # kx = 0.10 ; atmospheric extinction coefficient
        # end

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
        a = (0.00012638, 0., 0., 0.00012638)
        b = (-0.00012638, 0., 0., 0.00012638)
        c = (0.00012638, 0., 0., -0.00012638)
        d = (-0.00012638, 0., 0., -0.00012638)
        cds = dict(IM4=a,
                   IM3=b,
                   IM2=c,
                   IM1=d,
                   IM8=a,
                   IM7=b,
                   IM6=c,
                   IM5=d,
                   IM9=d,
                   IM10=c,
                   IM11=b,
                   IM12=a,
                   IM13=d,
                   IM14=c,
                   IM15=b,
                   IM16=a,
                   )
        return cds[ext]

    def _fiducial_exptime(self, fid, band):

        if band == 'g':
            fid.update(
                k_co = 0.17,
                A_co = 3.214,
                )

        elif band == 'r':
            fid.update(
                k_co = 0.10,
                A_co = 2.165,
                )

        elif band == 'z':
            fid.update(
                k_co = 0.060,
                A_co = 1.592,
                )
        else:
            raise ValueError('Unknown band "%s"' % band)

        return fid

