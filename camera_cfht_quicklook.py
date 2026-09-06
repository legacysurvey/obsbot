from obsbot import NominalCalibration, NominalExptime

camera_name = 'cfht-quicklook'

nice_camera_name = 'CFHT Quicklook'
# minimum number of extensions in a valid raw FITS file from this cam
min_n_exts = 1

bot_name = 'cfhtbot'

data_env_var = 'CFHT_DATA'

database_filename = camera_name + '.sqlite3'

class CFHTQuicklookNominalCalibration(NominalCalibration):
    def __init__(self):
        self.pixscale = 0.185
        self.overhead = 40
        self.gain = 2.5
        self.saturation_adu = 15000
        self.zp0 = dict(
            M4376= 24.78,
        )
        self.sky0 = dict(
            M4376 = 22.00,
        )

    def zeropoint(self, band, ext=None):
        return self.zp0[band]

    def sky(self, band):
        return self.sky0[band]

    def cdmatrix(self, ext):
        # approx...
        science = (5.194E-05, 0., 0., -5.194E-05)
        return science

    def fiducial_exptime(self, band):
        #def _fiducial_exptime(self, fid, band)
        fid = NominalExptime()
        print('fid: band', band)
        fid.update(seeing=1.0)
        if band == 'M4376':
            fid.update(
                k_co = 0.273,
                A_co = 4.103,
            )
        else:
            raise ValueError('Unknown band "%s"' % band)
        fid.update(skybright = self.sky(band))
        return fid

nominal_cal = CFHTQuicklookNominalCalibration()

default_extension = 'ccd00'

def ephem_observer():
    import ephem
    import numpy as np
    # Pyephem set-up:
    # coords from airmass.org: LAT 19 49 30.96, LON -155 28 07.67 ALT 4215
    cfht = ephem.Observer()
    cfht.lon = '-155.46879722'
    cfht.lat = '19.82526666'
    cfht.elev = 4215.0 # meters
    R_earth = 6378.1e3 # in meters
    cfht.horizon = -np.sqrt(2.0*cfht.elev/R_earth)
    return cfht

tile_path = 'obstatus/cfht-tiles.ecsv'
