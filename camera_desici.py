from obsbot import NominalCalibration
from obsbot import NominalExptime

camera_name = 'desici'

nice_camera_name = 'DESI-CI'
# minimum number of extensions in a valid raw FITS file from this cam
min_n_exts = 1
default_extension = 'CIC'
default_primary_extension = 1

copilot_plot_args = dict(label_nmatched=False, max_seeing=4., target_exptime=False,
                         nominal_sky=False)

database_filename = camera_name + '.sqlite3'

tile_path = None
data_env_var = 'MPC_DATA'
bot_name = 'cibot'

def ephem_observer():
    import ephem
    import numpy as np
    # pyephem
    mayall = ephem.Observer()
    mayall.lon = '-111.6003'
    mayall.lat = '31.9634'
    mayall.elev = 2120.0 # meters
    #mayall.temp = 10.0 # deg celsius; average temp for August
    #mayall.pressure = 780.0 # mbar
    R_earth = 6378.1e3 # in meters
    mayall.horizon = -np.sqrt(2.0*mayall.elev/R_earth)
    return mayall

desici_nominal_pixscale = 0.128

class DesiCiNominalCalibration(NominalCalibration):
    def __init__(self):
        self.pixscale = desici_nominal_pixscale
        self.saturation_adu = 20000
        self.zp0 = dict(
            r = 26.0)
        self.sky0 = dict(
            r = 20.0)

    def zeropoint(self, band, ext=None):
        if ext is not None:
            try:
                return self.zp0[(band, ext)]
            except KeyError:
                pass
        return self.zp0[band]

    def sky(self, band):
        return self.sky0[band]

    def fiducial_exptime(self, band):
        if band == 'r':
            fid = NominalExptime()
            fid.update(
                seeing = 1.3,
                skybright = self.sky(band),
                exptime = 10.,
                exptime_min = 1.,
                exptime_max = 100.,
                k_co = 0.109,
                A_co = 2.165,
                )
            return fid

    def cdmatrix(self, ext):
        # listhead DECam_00488199.fits.fz | grep 'EXTNAME\|CD._.'
        cds = dict(CIS = (-3.555E-05, 0., 0., 3.285E-05),
                   CIE = (0., -3.285E-05, -3.555E-05, 0.),
                   CIN = (3.555E-05, 0., 0., -3.285E-05),
                   CIW = (0., 3.285E-05, 3.555E-05, 0.),
                   CIC = (3.69E-05, 0., 0., -3.69E-05))
        return cds[ext]


nominal_cal = DesiCiNominalCalibration()

