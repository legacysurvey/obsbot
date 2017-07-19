from mosaic import *

camera_name = 'mosaic3'

nice_camera_name = 'Mosaic3'
# minimum number of extensions in a valid raw FITS file from this cam
min_n_exts = 16

bot_name = 'mosbot'

data_env_var = 'MOS3_DATA'

database_filename = camera_name + '.sqlite3'

nominal_cal = MosaicNominalCalibration()

default_extension = 'im4'

def fix_expnums(expnum):
    '''expnum: numpy array, fixed in-place.
    '''
    import numpy as np
    # Fix 6-digit 3xxxxx exposure numbers
    I = np.flatnonzero((expnum >= 300000) * (expnum < 400000))
    if len(I):
        #print('Set range of EXPNUMs', expnum[I].min(), expnum[I].max())
        expnum[I] -= 300000
        #print('  to', expnum[I].min(), expnum[I].max())

    # Also, in the mosaic obstatus file, there are some in the range
    # 260,223 to 260,355
    I = np.flatnonzero((expnum >= 260000) * (expnum < 270000))
    if len(I):
        expnum[I] -= 200000
        
    # Fix 7-digit 3xxxxxx exposure numbers
    I = np.flatnonzero((expnum >= 3000000))
    if len(I):
        #print('Set range of EXPNUMs', expnum[I].min(), expnum[I].max())
        expnum[I] -= 3000000
        #print('  to', expnum[I].min(), expnum[I].max())
    return expnum
        
def dradec_to_ref_chip(T, refext = 'im16'):
    import numpy as np
    nom = nominal_cal
    crx, cry = T.affine_x0, T.affine_y0
    (refcrx, refcry) = nom.crpix(refext)
    # Distance between im16 and this chip
    dcrx = refcrx - crx
    dcry = refcry - cry
    # Apply rotation
    dcx = T.affine_dxx * dcrx + T.affine_dxy * dcry - dcrx
    dcy = T.affine_dyx * dcrx + T.affine_dyy * dcry - dcry
    # Predicted pixel shift in im16
    cdx = T.dx - dcx
    cdy = T.dy - dcy
    # Convert to dRA, dDec in im16
    CDs = dict([(ext, nom.cdmatrix(ext)) for ext in np.unique(T.extension)])
    CD = np.array([CDs[ext] for ext in T.extension])
    xdra  = (CD[:,0] * cdx + CD[:,1] * cdy) * 3600.
    xddec = (CD[:,2] * cdx + CD[:,3] * cdy) * 3600.
    # Apply remaining mosstat - copilot offset
    offsets = dict(
        im16 = (-5.4, 3.4),
        )
    offra,offdec = offsets.get(refext, (0,0))
    refdra  = xdra  + offra
    refddec = xddec + offdec
    return refdra, refddec

def ephem_observer():
    import ephem
    import numpy as np
    # Pyephem set-up for mosaic:
    mosaic = ephem.Observer()
    mosaic.lon = '-111.6003'
    mosaic.lat = '31.9634'
    mosaic.elev = 2120.0 # meters
    #mosaic.temp = 10.0 # deg celsius; average temp for August
    #mosaic.pressure = 780.0 # mbar
    R_earth = 6378.1e3 # in meters
    mosaic.horizon = -np.sqrt(2.0*mosaic.elev/R_earth)
    return mosaic

tile_path = 'obstatus/mosaic-tiles_obstatus.fits'
