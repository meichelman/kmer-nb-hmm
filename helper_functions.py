import numpy as np
import os

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Functions for handling observertions/bed files
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
def load_obs_mut(obs_file, mutrates_file, window_size):
    with open(obs_file) as infile:
        num_windows = sum(1 for _ in infile)
            
    obs_arr = np.zeros(num_windows, dtype=np.int16)
    
    with open(obs_file) as infile:
        for idx, line in enumerate(infile):
            count = line.split('\t')[3]
            obs_arr[idx] = int(count)

    mutrates_arr = np.zeros(num_windows, dtype=float)

    with open(mutrates_file) as infile:
        for line in infile:
            if line.startswith('contig'):
                continue
            contig, start, end, mutrate = line.strip().split('\t')
            start, end, mutrate = int(start), int(end), float(mutrate)

            # Convert mutation-rate windows to observation windows
            obs_idx_start = start // window_size
            obs_idx_end = (end + window_size - 1) // window_size 
            obs_idx_end = min(obs_idx_end, num_windows)

            mutrates_arr[obs_idx_start:obs_idx_end] = mutrate

    for idx, (obs, m) in enumerate(zip(obs_arr, mutrates_arr)):
        if obs > 0 and m == 0:
            print(f"Warning: observation={obs} but mutrate=0 at index={idx}")

    return obs_arr, mutrates_arr

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Various helper functions
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------

def Make_folder_if_not_exists(path):
    '''
    Check if path exists - otherwise make it;
    '''
    path = os.path.dirname(path)
    if path != '':
        if not os.path.exists(path):
            os.makedirs(path)
