from obsbot import NominalCalibration

camera_name = 'pointingcam'

nice_camera_name = 'MayallPointing'
# minimum number of extensions in a valid raw FITS file from this cam
min_n_exts = 1
default_extension = 0

database_filename = camera_name + '.sqlite3'

tile_path = None
data_env_var = 'MPC_DATA'
bot_name = 'pointbot'

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


pointing_nominal_pixscale = 8.7

class PointingCamNominalCalibration(NominalCalibration):
    def __init__(self):
        self.pixscale = pointing_nominal_pixscale
        self.saturation_adu = 20000
        self.zp0 = dict(
            r = 25.0)
        self.sky0 = dict(
            r = 20.91)

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
        if band == 'r':
            fid.update(
                k_co = 0.109,
                A_co = 2.165,
                )

nominal_cal = PointingCamNominalCalibration()

