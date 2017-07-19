from decam import *

camera_name = 'decam'

nice_camera_name = 'DECam'
# arbitrary minimum number of extensions in a valid raw FITS file from this cam
min_n_exts = 60

bot_name = 'decbot'

data_env_var = 'DECAM_DATA'

database_filename = camera_name + '.sqlite3'

nominal_cal = DecamNominalCalibration()

default_extension = 'N4'

def ephem_observer():
    import ephem
    import numpy as np
    # Pyephem set-up for DECam:
    decam = ephem.Observer()
    decam.lon = '-70.806525'
    decam.lat = '-30.169661'
    decam.elev = 2207.0 # meters
    #decam.temp = 10.0 # deg celsius; average temp for August
    #decam.pressure = 780.0 # mbar
    R_earth = 6378.1e3 # in meters
    decam.horizon = -np.sqrt(2.0*decam.elev/R_earth)
    return decam

tile_path = 'obstatus/decam-tiles_obstatus.fits'
