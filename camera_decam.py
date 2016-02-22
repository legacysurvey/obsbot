from decam import *

camera_name = 'decam'

nominal_cal = DecamNominalCalibration()

default_extension = 'N4'

def ephem_observer():
    import ephem
    # Pyephem set-up for DECam:
    decam = ephem.Observer()
    decam.lon = '-70.806525'
    decam.lat = '-30.169661'
    decam.elev = 2207.0 # meters
    #decam.temp = 10.0 # deg celsius; average temp for August
    #decam.pressure = 780.0 # mbar
    #decam.horizon = -np.sqrt(2.0*decam.elev/R_earth)
    return decam
