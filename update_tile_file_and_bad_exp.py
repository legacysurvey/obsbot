'''

This script reads the db log files, identifies bad exposures, and updates the effective exposure times in the tiles file. 

'''

import numpy as np
import glob
from astropy.table import Table, join
import matplotlib.pyplot as plt
import os
from pathlib import Path


def identify_bad_exposures(log_file, bad_exp_file):
    '''
    Identifies the bad exposures to be repeated from the wide survey.
    '''

    # Conditions the exposures need to satisfy: seeing < 2.0, Teff > 100, and survey speed > 20%
    below_seeing_bool = log_file['seeing'] < 2.0
    above_t_eff_bool = log_file['efftime'] > 100.
    above_survey_speed_bool = log_file['efftime']/log_file['exptime'] > 0.20

    out_of_spec_bool = ~(below_seeing_bool * above_t_eff_bool * above_survey_speed_bool)

    return out_of_spec_bool


def sort_bad_exp_file(filename):
    '''
    Sorts the bad exposure file on expnum
    '''    
    if not os.path.exists(filename) or os.path.getsize(filename) == 0:
        return  # nothing to do

    # Read table (use same format settings as before)
    table = Table.read(filename, format="ascii.commented_header", guess=False)

    if len(table) == 0:
        return

    # Sort by expnum
    table.sort("expnum")

    # Overwrite file with sorted content
    table.write(filename, format="ascii.commented_header", overwrite=True)


def update_bad_expid_file(bad_exp_file, bad_obs_table):
    '''
    Updates the bad exposure ID file
    '''

    expnum = bad_obs_table['expnum']
    seeing = bad_obs_table['seeing']
    efftime = bad_obs_table['efftime']
    survey_speed = efftime / bad_obs_table['exptime']
    comment = np.full(len(bad_obs_table), 'Automatically flagged')
    
    # The information that will be written to bad exposure file
    data = np.column_stack((expnum, seeing, efftime, survey_speed, comment))

    # Read existing expnums in bad_exp_file
    existing_expnums = set()
    if os.path.exists(bad_exp_file):
        with open(bad_exp_file, "r") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                existing_expnums.add(int(line.split()[0]))

    # Remove rows that are already in bad_exp_file
    mask = ~np.isin(expnum, list(existing_expnums))
    new_data = data[mask]

    if len(new_data) == 0: # Nothing new to add
        sort_bad_exp_file(bad_exp_file) # Make sure sorting is in order
        return  
    
    # Write bad exposures to file
    with open(bad_exp_file, "a") as f:
        for i, row in enumerate(new_data):
            f.write(f"{int(new_data[i][0])} "
                    f"{float(new_data[i][1]):.3f} "
                    f"{float(new_data[i][2]):.1f} "
                    f"{float(new_data[i][3]):.3f} "
                    f"\"{new_data[i][4]}\"\n")

    # Sort by expnum
    sort_bad_exp_file(bad_exp_file)


def update_tile_and_bad_exp_file(logs_dir_path, tile_file, bad_exp_file):
    '''
    Reads db-*.ecsv log files and updates the efftime_tot in tile file. 
    Identifies bad exposures and writes bad exposure file.
    '''

    print('Updating tile & bad exposure file')
    db_log_files = glob.glob(logs_dir_path+'db-*.ecsv')

    fmt = 'ascii.ecsv'
    tiles = Table.read(tile_file, format=fmt)

    # Check if bad exposure file already exists, if not create one
    if not Path(bad_exp_file).exists():
        with open(bad_exp_file, "w") as f:
            f.write("# expnum seeing efftime survey_speed comment\n")
            
    updated_efftime_array = np.full(len(tiles), 0.0)
    done_array = np.full(len(tiles), False)
    count = 0 # To keep count of number of bad exposures

    for i, file_name in enumerate(db_log_files): # Loop over all nights

        log_file = Table.read(file_name, format=fmt)

        # Select wide & deep survey exposures 
        ibis_wide_survey_idx = np.where(np.char.startswith(log_file['object'], "IBIS_wide_"))[0]
        ibis_deep_survey_idx = np.where(np.char.startswith(log_file['object'], "IBIS_deep_"))[0]
        survey_idx = np.concatenate((ibis_wide_survey_idx, ibis_deep_survey_idx), axis=None)

        log_file = log_file[survey_idx]
        
        out_of_spec_boolean = identify_bad_exposures(log_file, bad_exp_file)

        # Add the out of spec exposures to the bad exposure file, if they are not already in there
        update_bad_expid_file(bad_exp_file, log_file[out_of_spec_boolean]) 

        # Read the file to ensure also bad exposures added by hand are removed from tile file
        bad_exps =  Table.read(bad_exp_file, format="ascii.basic", guess=False, names=["expnum", "seeing", "efftime", "survey_speed", "comment"])

        # Set all eff exposure time in tile file to zero & update with new efftime in tile file 
        for i, row in enumerate(log_file):
                
            if log_file['expnum'][i] in bad_exps['expnum']: # Skip if in bad exposure file
                count += 1
                continue

            idx = np.where(tiles['OBJECT'] == log_file[i]['object'])[0]
            updated_efftime_array[idx] += log_file[i]['efftime']

            if len(idx) == 0:
                print('No object ID available')
                print(log_file[i])

            if len(idx) > 1:
                print('Multiple object IDs available')

            # Set to done if exposure time 0.75x"EFFTIME_GOAL" is reached, set "DONE" to 1
            if updated_efftime_array[idx] >= 0.75*tiles['EFFTIME_GOAL'][idx]:   
                done_array[idx] = True

    print('Num tiles "DONE" = ', np.sum(done_array), ' out of ', len(done_array), ' in the survey')
    print('Total number of bad exposures = ', count)

    # Write the new efftime array to the EFFTIME_TOT column
    tiles['EFFTIME_TOT'] = np.round(updated_efftime_array, decimals=2)
    tiles.write(tile_file, overwrite=True)

    print('Wrote new tile file. Done!')


def main(ibis_tile_file, logs_dir_path, bad_exp_filename):
    '''
    ibis_tile_file (str) - path to IBIS tile file (ibis-observing/obstatus/ibis-tiles.ecsv) 
    logs_dir_path (str) - path to logs directory with db-* files (ibis-observing/logs/)
    bad_exp_filename (str) - path to bad exposure file (ibis-observing/obstatus/bad-exp-file.txt)
    '''

    update_tile_and_bad_exp_file(logs_dir_path, ibis_tile_file, bad_exp_filename)


