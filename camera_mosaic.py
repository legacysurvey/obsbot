from mosaic import *

camera_name = 'mosaic3'

nominal_cal = MosaicNominalCalibration()

default_extension = 'im4'

def ephem_observer():
    import ephem
    # Pyephem set-up for mosaic:
    mosaic = ephem.Observer()
    mosaic.lon = '-111.6003'
    mosaic.lat = '31.9634'
    mosaic.elev = 2120.0 # meters
    #mosaic.temp = 10.0 # deg celsius; average temp for August
    #mosaic.pressure = 780.0 # mbar
    #mosaic.horizon = -np.sqrt(2.0*mosaic.elev/R_earth)
    return mosaic
