import os
import glob
import numpy as np
from netCDF4 import Dataset
from datetime import datetime, timedelta
from multiprocessing import Pool, Manager
import re

def read_value_from_config(key_to_find, file_path):
    """Reads a floating-point value from a config file given a specific key."""
    with open(file_path, 'r') as file:
        for line in file:
            if key_to_find in line:
                parts = line.split('=')
                if len(parts) > 1:
                    value_part = parts[1].split(',')[0]
                    match = re.search(r'[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?', value_part)
                    if match:
                        return float(match.group(0))
                    else:
                        print("No valid number found after the key.")
                        return None
        print(f"Key '{key_to_find}' not found in the file.")
        return None

# ---------------------------------------------------------------------------
# Global parameter setup (original logic preserved)
# ---------------------------------------------------------------------------
# config_path = '../run.conf'
# DX = read_value_from_config('DX', config_path)
# DY = read_value_from_config('DY', config_path)
# DZ = read_value_from_config('DZ', config_path)

# Time step interval in milliseconds (original).
TIME_STEP_INTERVAL = 1000  # or 100 for coalescence-event tracing

# Global structured dtype for particle data.
particle_dtype = np.dtype([
    ('sd_r',           np.float64),
    ('sd_n',           np.int64),
    ('sd_x',           np.float64),
    ('sd_y',           np.float64),
    ('sd_z',           np.float64),
    ('pre_sdid',       np.int32),
    ('pre_dmid',       np.int32),
    ('if_coal',        np.int8),
    ('other_pre_dmid', np.int32),
    ('other_pre_sdid', np.int32),
    ('num_col',        np.int32),
    ('other_sd_r',     np.float64),
    ('other_sd_n',     np.int64)
])

def generate_time_steps(start_time, end_time, step_ms=TIME_STEP_INTERVAL):
    """
    Generates a list of time strings (HHMMSS.mmm) from start_time down to end_time
    in step_ms millisecond increments.
    """
    time_format = "%H%M%S.%f"
    current_time = datetime.strptime(start_time, time_format)
    end_time = datetime.strptime(end_time, time_format)
    time_delta = timedelta(milliseconds=step_ms)

    time_steps = []
    if current_time <= end_time:
        while current_time <= end_time:
            time_steps.append(current_time.strftime(time_format)[:-3])
            current_time += time_delta
    else:
        while current_time >= end_time:
            time_steps.append(current_time.strftime(time_format)[:-3])
            current_time -= time_delta

    return time_steps

def generate_coal_filename(time, domain_id, input_directory):
    domain_id_str = f"{domain_id:06d}"
    fixed_name = os.path.join(input_directory, f"SD_coal_output_NetCDF_00000101-{time}.pe{domain_id_str}")
    if os.path.isfile(fixed_name):
        return fixed_name
    matched = sorted(glob.glob(os.path.join(input_directory, f"SD_coal_output_NetCDF_*-{time}.pe{domain_id_str}")))
    if matched:
        return matched[-1]
    return fixed_name

def generate_filename(time, domain_id, input_directory):
    domain_id_str = f"{domain_id:06d}"
    selected_file = os.path.join(input_directory, f"SD_selected_NetCDF_00000101-{time}.pe{domain_id_str}")
    if not os.path.isfile(selected_file):
        selected_candidates = sorted(glob.glob(os.path.join(input_directory, f"SD_selected_NetCDF_*-{time}.pe{domain_id_str}")))
        if selected_candidates:
            selected_file = selected_candidates[-1]
    if os.path.isfile(selected_file):
        return selected_file
    all_file = os.path.join(input_directory, f"SD_all_NetCDF_00000101-{time}.pe{domain_id_str}")
    if os.path.isfile(all_file):
        return all_file
    all_candidates = sorted(glob.glob(os.path.join(input_directory, f"SD_all_NetCDF_*-{time}.pe{domain_id_str}")))
    if all_candidates:
        return all_candidates[-1]
    return all_file

# ---------------------------------------------------------------------------
# Below are the core algorithm functions. We apply caching to reduce I/O overhead.
# ---------------------------------------------------------------------------

def open_netcdf_cached(filename, cache):
    """
    Opens a NetCDF file and stores the Dataset object in a cache dict.
    If the file is already opened, returns the cached object instead.
    The caller is responsible for closing all cached files once done.
    """
    if filename not in cache:
        cache[filename] = Dataset(filename, 'r')
    return cache[filename]

def read_tracking_var(nc, names, prefixes=None):
    for name in names:
        if name in nc.variables:
            return nc.variables[name][:]
    if prefixes is not None:
        all_names = list(nc.variables.keys())
        for prefix in prefixes:
            matched = sorted([name for name in all_names if name.startswith(prefix)])
            if matched:
                return nc.variables[matched[-1]][:]
    return None

def get_particle_data(filename, sd_id, cache=None):
    """
    Retrieves particle information for a specific sd_id from the specified
    NetCDF file. If a cache dict is provided, repeated open/close is avoided.
    """
    try:
        # Use the cached NetCDF dataset if available
        if cache is not None:
            if filename not in cache:
                # Open the file once and only read required fields
                nc = Dataset(filename, 'r')
                cache[filename] = {
                    'sd_r': nc.variables['sd_r'][:],
                    'sd_n': nc.variables['sd_n'][:],
                    'sd_x': nc.variables['sd_x'][:],
                    'sd_y': nc.variables['sd_y'][:],
                    'sd_z': nc.variables['sd_z'][:],
                    'pre_sdid': read_tracking_var(nc, ['pre_sdid', 'sd_id'], ['pre_sdid_', 'sd_id_']),
                    'pre_dmid': read_tracking_var(nc, ['pre_dmid', 'dm_id'], ['pre_dmid_', 'dm_id_']),
                    'if_coal': read_tracking_var(nc, ['if_coal'], ['if_coal_']),
                }
                nc.close()
            all_data = cache[filename]
        else:
            nc = Dataset(filename, 'r')
            all_data = {
                'sd_r': nc.variables['sd_r'][:],
                'sd_n': nc.variables['sd_n'][:],
                'sd_x': nc.variables['sd_x'][:],
                'sd_y': nc.variables['sd_y'][:],
                'sd_z': nc.variables['sd_z'][:],
                'pre_sdid': read_tracking_var(nc, ['pre_sdid', 'sd_id'], ['pre_sdid_', 'sd_id_']),
                'pre_dmid': read_tracking_var(nc, ['pre_dmid', 'dm_id'], ['pre_dmid_', 'dm_id_']),
                'if_coal': read_tracking_var(nc, ['if_coal'], ['if_coal_']),
            }
            nc.close()

        idx = sd_id - 1  # Convert to 0-based index
        data = {key: all_data[key][idx] for key in all_data}

        return data
    except Exception as e:
        print(f"Error processing {filename}: {e}")
        return None

def check_coalescence(coal_file, pre_dmid, pre_sdid, current_file, cache=None, missing_cache=None):
    """
    Reads the coalescence-event NetCDF file to obtain coalescence info.
    If file doesn't exist, or if no matching record is found, returns (None, None, None).
    """
    # --- 1) Quick file-existence check --------------------------------------
    if coal_file in missing_cache:
        return None, None, None  # skip quickly
    if not os.path.isfile(coal_file):
        missing_cache.append(coal_file)  # Add to the list of missing files
        # Print a short message (or none) to keep logs minimal:
        print(f"[Info] Coalescence file not found: {coal_file} (skipping)")
        return None, None, None

    try:
        # --- 2) Use cached Dataset if provided --------------------------------
        if cache is not None:
            nc = open_netcdf_cached(coal_file, cache)
        else:
            nc = Dataset(coal_file, 'r')
        # ----------------------------------------------------------------------

        # Retrieve arrays
        pre_sdid1 = read_tracking_var(nc, ['pre_sdid1', 'sd_id1'], ['pre_sdid1_', 'sd_id1_'])
        pre_dmid1 = read_tracking_var(nc, ['pre_dmid1', 'dm_id1'], ['pre_dmid1_', 'dm_id1_'])
        pre_sdid2 = read_tracking_var(nc, ['pre_sdid2', 'sd_id2'], ['pre_sdid2_', 'sd_id2_'])
        pre_dmid2 = read_tracking_var(nc, ['pre_dmid2', 'dm_id2'], ['pre_dmid2_', 'dm_id2_'])
        num_col = nc.variables['num_col'][:]

        # Identify records matching (pre_dmid, pre_sdid)
        mask1 = (pre_dmid1 == pre_dmid) & (pre_sdid1 == pre_sdid)
        mask2 = (pre_dmid2 == pre_dmid) & (pre_sdid2 == pre_sdid)

        if np.any(mask1):
            index = np.where(mask1)[0][0]
            result_dmid, result_sdid, result_numcol = pre_dmid2[index], pre_sdid2[index], num_col[index]
        elif np.any(mask2):
            index = np.where(mask2)[0][0]
            result_dmid, result_sdid, result_numcol = pre_dmid1[index], pre_sdid1[index], num_col[index]
        else:
            # Print a minimal info or skip printing entirely:
            print(f"[Info] No matching coalescence for (pre_dmid={pre_dmid}, pre_sdid={pre_sdid}) in {coal_file}.")
            if cache is None:
                nc.close()
            return None, None, None

        if cache is None:
            nc.close()

        return result_dmid, result_sdid, result_numcol

    except Exception as e:
        # Print or log a minimal message to avoid clutter
        print(f"[Warning] Failed to open or read coalescence file: {coal_file} - {e}")
        return None, None, None

def process_file_for_particles(filename):
    """
    Worker function to process one NetCDF file.
    Returns a list of tuples: (filename, idx, domain_id) for valid particles.
    A particle is valid if its sd_z > 0.
    """
    results = []
    try:
        parts = filename.split(".pe")
        if len(parts) < 2:
            return results
        domain_id = int(parts[1])
        with Dataset(filename, 'r') as nc:
            sd_z = nc.variables['sd_z'][:]
            valid_idx = np.where(sd_z > 0)[0]
            for idx in valid_idx:
                results.append((filename, int(idx), domain_id))
    except Exception as e:
        print(f"[Error] Processing file {filename} failed: {e}")
    return results

def initialize_particles(input_directory, time_steps):
    """
    Initializes particles at the initial time step by scanning NetCDF files.
    
    Returns:
      - part_array: a structured array of shape (num_times, num_particles) with initial particle data.
      - global_mapping: a dictionary mapping (filename, idx) to global particle id.
      - particle_infos: a list of tuples (filename, idx, domain_id).
    
    Notes:
      - The valid indices from np.where are 0-based.
      - The fields pre_sdid and pre_dmid read from the file are 1-based.
      - Uses parallel processing to scan files.
    """
    initial_time = time_steps[0]
    selected_pattern = os.path.join(input_directory, f"SD_selected_NetCDF_00000101-{initial_time}.pe*")
    files = glob.glob(selected_pattern)
    if len(files) == 0:
        file_pattern = os.path.join(input_directory, f"SD_all_NetCDF_00000101-{initial_time}.pe*")
        files = glob.glob(file_pattern)
    
    with Pool() as pool:
        results = pool.map(process_file_for_particles, files)
    
    particle_infos = [item for sublist in results for item in sublist]
    
    if len(particle_infos) == 0:
        print("No particle found satisfying the conditions.")
        exit()
    
    num_particles = len(particle_infos)
    num_times = len(time_steps)
    
    # Build global mapping as a dictionary.
    global_mapping = {}
    for i, (fname, idx, _) in enumerate(particle_infos):
        global_mapping[(fname, idx)] = i

    # Allocate structured array for initial particle data.
    part_array = np.empty((num_times, num_particles), dtype=particle_dtype)
    part_array['sd_r'][0, :]     = -1.0
    part_array['sd_n'][0, :]     = -1
    part_array['sd_x'][0, :]     = -1.0
    part_array['sd_y'][0, :]     = -1.0
    part_array['sd_z'][0, :]     = -1.0
    part_array['pre_sdid'][0, :] = -1
    part_array['pre_dmid'][0, :] = -1
    part_array['if_coal'][0, :]  = -1
    part_array['other_pre_dmid'][0, :] = -1
    part_array['other_pre_sdid'][0, :] = -1
    part_array['num_col'][0, :]        = -1
    part_array['other_sd_r'][0, :]     = -1.0
    part_array['other_sd_n'][0, :]     = -1

    # Fill initial data from each file.
    for i, (fname, idx, domain_id) in enumerate(particle_infos):
        try:
            with Dataset(fname, 'r') as nc:
                part_array['sd_r'][0, i] = nc.variables['sd_r'][idx]
                part_array['sd_n'][0, i] = nc.variables['sd_n'][idx]
                part_array['sd_x'][0, i] = nc.variables['sd_x'][idx]
                part_array['sd_y'][0, i] = nc.variables['sd_y'][idx]
                part_array['sd_z'][0, i] = nc.variables['sd_z'][idx]
                pre_sdid_arr = read_tracking_var(nc, ['pre_sdid', 'sd_id'], ['pre_sdid_', 'sd_id_'])
                pre_dmid_arr = read_tracking_var(nc, ['pre_dmid', 'dm_id'], ['pre_dmid_', 'dm_id_'])
                if_coal_arr = read_tracking_var(nc, ['if_coal'], ['if_coal_'])
                part_array['pre_sdid'][0, i] = pre_sdid_arr[idx]
                part_array['pre_dmid'][0, i] = pre_dmid_arr[idx]
                part_array['if_coal'][0, i]  = if_coal_arr[idx]
        except Exception as e:
            print(f"[Error] Reading initial data for {fname}, idx={idx}: {e}")
    
    return part_array, global_mapping, particle_infos

def process_single_particle_block(particle_block, input_directory, time_steps, num_times,
                                  output_dir, block_idx, missing_cache, global_mapping):
    """
    Processes a subset (block) of particles in parallel, writing each time step's data
    to a separate NetCDF file. Incorporates:
      (1) Domain-based grouping to reduce I/O on the same NetCDF file,
      (2) Dictionary-based coalescence lookups,
      (4) Memory management improvements (clearing unused data).

    Parameters
    ----------
    particle_block : list of (filename, idx, domain_id)
        The subset of particles (by file+index) that this block is responsible for.
    input_directory : str
        Directory containing the NetCDF files.
    time_steps : list of str
        The time-step strings (e.g. ["001000.000", "000900.000", ...]) in descending order.
    num_times : int
        The number of time steps to process.
    output_dir : str
        Where to write temporary output files for this block.
    block_idx : int
        The block index (for naming output files).
    missing_cache : multiprocessing.Manager().list
        A shared list for tracking missing coalescence files or other files to avoid repeated logs.
    global_mapping : dict
        Maps (filename, idx) -> global particle ID for referencing across blocks.

    Returns
    -------
    None
        (Writes one NetCDF file per time step: temp_block_{block_idx}_time_{time_steps[t]}.nc)
    """
    # A cache for SD_all_NetCDF data: {filename: {var: array, ...}}
    netcdf_cache = {}
    # A cache for coalescence data: {(time_str, domain_id): dict_of_{(dmid, sdid) -> (odmid, osdid, num_col)}}
    coal_cache = {}

    # We allocate memory for each time step as a dictionary of lists.
    # Key structure matches your original script.
    particles_by_time = [
        {
            'particle_index': [],
            'sd_x': [], 'sd_y': [], 'sd_z': [],
            'sd_r': [], 'sd_n': [],
            'pre_sdid': [], 'pre_dmid': [],
            'if_coal': [], 'domain_id': [],
            'other_pre_dmid': [], 'other_pre_sdid': [],
            'num_col': [], 'other_sd_r': [], 'other_sd_n': []
        }
        for _ in range(num_times)
    ]

    # ---------------------
    # 1) Fill t=0
    # ---------------------
    for (filename, idx, domain_id) in particle_block:
        pindex = global_mapping[(filename, idx)]
        try:
            # Open or re-use netCDF file from cache
            if filename not in netcdf_cache:
                from netCDF4 import Dataset
                nc = Dataset(filename, 'r')
                netcdf_cache[filename] = {
                    'sd_r': nc.variables['sd_r'][:],
                    'sd_n': nc.variables['sd_n'][:],
                    'sd_x': nc.variables['sd_x'][:],
                    'sd_y': nc.variables['sd_y'][:],
                    'sd_z': nc.variables['sd_z'][:],
                    'pre_sdid': read_tracking_var(nc, ['pre_sdid', 'sd_id'], ['pre_sdid_', 'sd_id_']),
                    'pre_dmid': read_tracking_var(nc, ['pre_dmid', 'dm_id'], ['pre_dmid_', 'dm_id_']),
                    'if_coal': read_tracking_var(nc, ['if_coal'], ['if_coal_']),
                }
                nc.close()
            arr = netcdf_cache[filename]

            # Read data for idx
            sd_x = arr['sd_x'][idx]
            sd_y = arr['sd_y'][idx]
            sd_z = arr['sd_z'][idx]
            sd_r = arr['sd_r'][idx]
            sd_nval = arr['sd_n'][idx]
            psdid = arr['pre_sdid'][idx]
            pdmid = arr['pre_dmid'][idx]
            if_coal_val = arr['if_coal'][idx]

            # Store in particles_by_time[0]
            particles_by_time[0]['particle_index'].append(pindex)
            particles_by_time[0]['sd_x'].append(sd_x)
            particles_by_time[0]['sd_y'].append(sd_y)
            particles_by_time[0]['sd_z'].append(sd_z)
            particles_by_time[0]['sd_r'].append(sd_r)
            particles_by_time[0]['sd_n'].append(sd_nval)
            particles_by_time[0]['pre_sdid'].append(psdid)
            particles_by_time[0]['pre_dmid'].append(pdmid)
            particles_by_time[0]['if_coal'].append(if_coal_val)
            particles_by_time[0]['domain_id'].append(int(filename.split(".pe")[1]))

            # Initialize "other_" coalescence fields as -1
            particles_by_time[0]['other_pre_dmid'].append(-1)
            particles_by_time[0]['other_pre_sdid'].append(-1)
            particles_by_time[0]['num_col'].append(-1)
            particles_by_time[0]['other_sd_r'].append(-1.0)
            particles_by_time[0]['other_sd_n'].append(-1)

        except Exception as e:
            print(f"[Error] Reading t=0 data from {filename}, idx={idx}, block={block_idx}: {e}")

    # ---------------------
    # 2) Fill t=1..(num_times-1)
    # ---------------------
    from collections import defaultdict
    from netCDF4 import Dataset

    for t in range(1, num_times):
        current_t_str = time_steps[t]
        prev_t_data = particles_by_time[t - 1]

        # 2.1) Group the particles by domain_id from (t-1),
        #      so we only open each domain's NetCDF file once.
        domain_groups = defaultdict(list)

        num_prev = len(prev_t_data['particle_index'])
        for i in range(num_prev):
            prev_dmid = prev_t_data['pre_dmid'][i]
            prev_sdid = prev_t_data['pre_sdid'][i]
            pindex    = prev_t_data['particle_index'][i]
            if prev_dmid > 0:  # Only group valid domain_id
                domain_groups[prev_dmid].append((i, pindex, prev_sdid))

        # We'll build an internal mapping from pindex -> row index in the new t,
        # so we can easily fill coalescence info after reading the main droplet data.
        pindex_map_t = {}

        # 2.2) For each domain group, read the corresponding file once
        for dmid, info_list in domain_groups.items():
            # Generate the SD_all_NetCDF filename
            current_file = generate_filename(current_t_str, dmid, input_directory)
            # Check existence
            if not os.path.isfile(current_file):
                print(f"[Warning] Missing file for block={block_idx}, time={current_t_str}, dmid={dmid}")
                # Fill placeholders for each particle in this domain group
                for (i_in_prev, pidx, psdid) in info_list:
                    row_idx_new = len(particles_by_time[t]['particle_index'])
                    pindex_map_t[pidx] = row_idx_new
                    # Create a blank entry
                    particles_by_time[t]['particle_index'].append(pidx)
                    # Fill -1 or -1.0
                    for key in particles_by_time[t]:
                        if key == 'particle_index':
                            continue
                        if key.startswith('sd_'):
                            particles_by_time[t][key].append(-1.0)
                        else:
                            particles_by_time[t][key].append(-1)
                continue

            # Load from cache or read arrays
            if current_file not in netcdf_cache:
                try:
                    nc_data = Dataset(current_file, 'r')
                    netcdf_cache[current_file] = {
                        'sd_r': nc_data.variables['sd_r'][:],
                        'sd_n': nc_data.variables['sd_n'][:],
                        'sd_x': nc_data.variables['sd_x'][:],
                        'sd_y': nc_data.variables['sd_y'][:],
                        'sd_z': nc_data.variables['sd_z'][:],
                        'pre_sdid': nc_data.variables['pre_sdid'][:],
                        'pre_dmid': nc_data.variables['pre_dmid'][:],
                        'if_coal': nc_data.variables['if_coal'][:],
                    }
                    nc_data.close()
                except Exception as e:
                    print(f"[Error] Reading {current_file}: {e}")
                    # Fill placeholders
                    for (i_in_prev, pidx, psdid) in info_list:
                        row_idx_new = len(particles_by_time[t]['particle_index'])
                        pindex_map_t[pidx] = row_idx_new
                        for key in particles_by_time[t]:
                            if key == 'particle_index':
                                continue
                            if key.startswith('sd_'):
                                particles_by_time[t][key].append(-1.0)
                            else:
                                particles_by_time[t][key].append(-1)
                    continue

            arr_d = netcdf_cache[current_file]

            # For each particle in this domain group, read index
            for (i_in_prev, pidx, psdid) in info_list:
                row_idx_new = len(particles_by_time[t]['particle_index'])
                pindex_map_t[pidx] = row_idx_new

                if psdid < 1:
                    # No valid sdid => fill placeholders
                    particles_by_time[t]['particle_index'].append(pidx)
                    for key in particles_by_time[t]:
                        if key == 'particle_index':
                            continue
                        if key.startswith('sd_'):
                            particles_by_time[t][key].append(-1.0)
                        else:
                            particles_by_time[t][key].append(-1)
                    continue

                idx_0b = psdid - 1
                sd_r_val   = arr_d['sd_r'][idx_0b]
                sd_n_val   = arr_d['sd_n'][idx_0b]
                sd_x_val   = arr_d['sd_x'][idx_0b]
                sd_y_val   = arr_d['sd_y'][idx_0b]
                sd_z_val   = arr_d['sd_z'][idx_0b]
                psdid_val  = arr_d['pre_sdid'][idx_0b]
                pdmid_val  = arr_d['pre_dmid'][idx_0b]
                if_coalval = arr_d['if_coal'][idx_0b]

                particles_by_time[t]['particle_index'].append(pidx)
                particles_by_time[t]['sd_x'].append(sd_x_val)
                particles_by_time[t]['sd_y'].append(sd_y_val)
                particles_by_time[t]['sd_z'].append(sd_z_val)
                particles_by_time[t]['sd_r'].append(sd_r_val)
                particles_by_time[t]['sd_n'].append(sd_n_val)
                particles_by_time[t]['pre_sdid'].append(psdid_val)
                particles_by_time[t]['pre_dmid'].append(pdmid_val)
                particles_by_time[t]['if_coal'].append(if_coalval)
                particles_by_time[t]['domain_id'].append(dmid)

                # Initialize coalescence placeholders
                particles_by_time[t]['other_pre_dmid'].append(-1)
                particles_by_time[t]['other_pre_sdid'].append(-1)
                particles_by_time[t]['num_col'].append(-1)
                particles_by_time[t]['other_sd_r'].append(-1.0)
                particles_by_time[t]['other_sd_n'].append(-1)

        # 2.3) Now handle coalescence for particles that had if_coal == 1 at (t-1)
        coals_needed = []
        for i in range(num_prev):
            if prev_t_data['if_coal'][i] == 1:
                pidx = prev_t_data['particle_index'][i]
                pdmid = prev_t_data['pre_dmid'][i]
                psdid = prev_t_data['pre_sdid'][i]
                coals_needed.append((pidx, pdmid, psdid))

        # Group these by domain_id
        from collections import defaultdict
        coal_groups = defaultdict(list)
        for (pidx, pdmid, psdid) in coals_needed:
            if pdmid > 0 and psdid > 0:
                coal_groups[pdmid].append((pidx, psdid))

        # For each domain, read the coalescence file once
        for dmid_coal, c_list in coal_groups.items():
            cfile = generate_coal_filename(current_t_str, dmid_coal, input_directory)
            if not os.path.isfile(cfile):
                if cfile not in missing_cache:
                    print(f"[Info] Coalescence file not found: {cfile}")
                    missing_cache.append(cfile)
                continue

            # If not cached, read & build dictionary
            if (current_t_str, dmid_coal) not in coal_cache:
                try:
                    with Dataset(cfile, 'r') as nc_c:
                        p1 = read_tracking_var(nc_c, ['pre_sdid1', 'sd_id1'], ['pre_sdid1_', 'sd_id1_'])
                        d1 = read_tracking_var(nc_c, ['pre_dmid1', 'dm_id1'], ['pre_dmid1_', 'dm_id1_'])
                        p2 = read_tracking_var(nc_c, ['pre_sdid2', 'sd_id2'], ['pre_sdid2_', 'sd_id2_'])
                        d2 = read_tracking_var(nc_c, ['pre_dmid2', 'dm_id2'], ['pre_dmid2_', 'dm_id2_'])
                        nc_ = nc_c.variables['num_col'][:]

                    c_dict = {}
                    for irow in range(len(p1)):
                        key1 = (d1[irow], p1[irow])
                        c_dict[key1] = (d2[irow], p2[irow], nc_[irow])
                        key2 = (d2[irow], p2[irow])
                        c_dict[key2] = (d1[irow], p1[irow], nc_[irow])

                    coal_cache[(current_t_str, dmid_coal)] = c_dict
                except Exception as e:
                    print(f"[Warning] Reading coalescence {cfile} failed: {e}")
                    continue

            c_dict = coal_cache.get((current_t_str, dmid_coal), {})
            # For each (pidx, psdid) in c_list, find the row in t
            for (pidx, psdid) in c_list:
                # row at time t
                row_t = pindex_map_t.get(pidx, None)
                if row_t is None:
                    continue
                # check if we have an entry
                if (dmid_coal, psdid) not in c_dict:
                    # no matching coalescence in file
                    continue

                (odmid, osdid, numc) = c_dict[(dmid_coal, psdid)]
                particles_by_time[t]['other_pre_dmid'][row_t] = odmid
                particles_by_time[t]['other_pre_sdid'][row_t] = osdid
                particles_by_time[t]['num_col'][row_t]        = numc

                # Optionally fetch "other" droplet data from SD_all_NetCDF...
                if odmid > 0 and osdid > 0:
                    other_file = generate_filename(current_t_str, odmid, input_directory)
                    if os.path.isfile(other_file):
                        if other_file not in netcdf_cache:
                            try:
                                ofnc = Dataset(other_file, 'r')
                                netcdf_cache[other_file] = {
                                    'sd_r': ofnc.variables['sd_r'][:],
                                    'sd_n': ofnc.variables['sd_n'][:],
                                    'sd_x': ofnc.variables['sd_x'][:],
                                    'sd_y': ofnc.variables['sd_y'][:],
                                    'sd_z': ofnc.variables['sd_z'][:],
                                    'pre_sdid': read_tracking_var(ofnc, ['pre_sdid', 'sd_id'], ['pre_sdid_', 'sd_id_']),
                                    'pre_dmid': read_tracking_var(ofnc, ['pre_dmid', 'dm_id'], ['pre_dmid_', 'dm_id_']),
                                    'if_coal': read_tracking_var(ofnc, ['if_coal'], ['if_coal_'])
                                }
                                ofnc.close()
                            except Exception as e:
                                print(f"[Error] {other_file} reading for coalescence: {e}")

                        arr_o = netcdf_cache.get(other_file, None)
                        if arr_o:
                            oidx_0b = osdid - 1
                            if 0 <= oidx_0b < len(arr_o['sd_r']):
                                # Check who is bigger, etc. (example logic)
                                main_r = particles_by_time[t]['sd_r'][row_t]
                                other_r = arr_o['sd_r'][oidx_0b]
                                if main_r < other_r:
                                    # Swap: the bigger droplet becomes the "main"
                                    particles_by_time[t]['other_sd_r'][row_t] = main_r
                                    particles_by_time[t]['other_sd_n'][row_t] = particles_by_time[t]['sd_n'][row_t]

                                    # Overwrite main droplet with bigger
                                    particles_by_time[t]['sd_r'][row_t]  = other_r
                                    particles_by_time[t]['sd_n'][row_t]  = arr_o['sd_n'][oidx_0b]
                                    particles_by_time[t]['if_coal'][row_t] = arr_o['if_coal'][oidx_0b]
                                    particles_by_time[t]['pre_sdid'][row_t] = arr_o['pre_sdid'][oidx_0b]
                                    particles_by_time[t]['pre_dmid'][row_t] = arr_o['pre_dmid'][oidx_0b]
                                    particles_by_time[t]['sd_x'][row_t] = arr_o['sd_x'][oidx_0b]
                                    particles_by_time[t]['sd_y'][row_t] = arr_o['sd_y'][oidx_0b]
                                    particles_by_time[t]['sd_z'][row_t] = arr_o['sd_z'][oidx_0b]
                                    particles_by_time[t]['domain_id'][row_t] = odmid
                                else:
                                    # Keep main, record the "other" droplet
                                    particles_by_time[t]['other_sd_r'][row_t] = other_r
                                    particles_by_time[t]['other_sd_n'][row_t] = arr_o['sd_n'][oidx_0b]

        # --------------------------------------------
        # 2.4) Write the results
        # --------------------------------------------
        temp_file = os.path.join(output_dir, f"temp_block_{block_idx}_time_{current_t_str}.nc")
        with Dataset(temp_file, 'w', format='NETCDF4') as nc_out:
            nparts = len(particles_by_time[t]['sd_x'])
            nc_out.createDimension('particle', nparts)

            if nparts > 0:
                # Write the "particle_index" variable
                var_idx = nc_out.createVariable('particle_index', np.int64, ('particle',))
                var_idx[:] = particles_by_time[t]['particle_index']

            # Define a small helper for picking dtype
            def pick_dtype(key):
                if key in ['sd_n', 'other_sd_n']:
                    return np.int64
                elif key in ['if_coal']:
                    return np.int8
                elif key in ['sd_x', 'sd_y', 'sd_z', 'sd_r', 'other_sd_r']:
                    return np.float32
                else:
                    return np.int32

            for key, arr in particles_by_time[t].items():
                if key == 'particle_index':
                    continue  # already written
                if len(arr) == 0:
                    continue  # skip if empty
                dtype_ = pick_dtype(key)
                var = nc_out.createVariable(key, dtype_, ('particle',), zlib=True, complevel=4)
                var[:] = arr

        # 2.5) Clean cache
        if t > 1:
            # we no longer need time step t-2
            particles_by_time[t-2].clear()

    if num_times >= 2:
        particles_by_time[num_times-2].clear()
    particles_by_time[num_times-1].clear()

    netcdf_cache.clear()
    coal_cache.clear()

def merge_temp_files(output_file, output_dir, time_steps, num_blocks, block_sizes):
    """
    Merges all temporary files into a final NetCDF file.
    
    The temporary file names are constructed proactively based on known time_steps
    and block indices. For each time step, the temporary file for block b is assumed to be:
        temp_block_{b}_time_{time_str}.nc

    The final file will have dimensions (time, total_particles), where total_particles is
    the sum of block_sizes. In addition, the function sorts the particles by their global id
    (particle_index) for each time step, and writes the global particle_index variable in the output.
    
    Parameters:
      output_file: Path to the final output NetCDF file.
      output_dir: Directory where temporary files are stored.
      time_steps: List of time strings.
      num_blocks: Total number of blocks.
      block_sizes: List of particle counts for each block.
    """
    num_times = len(time_steps)
    total_particles = sum(block_sizes)

    # Allocate final arrays for each variable.
    final_sd_r = np.full((num_times, total_particles), -1.0, dtype=np.float64)
    final_sd_n = np.full((num_times, total_particles), -1, dtype=np.int64)
    final_sd_x = np.full((num_times, total_particles), -1.0, dtype=np.float64)
    final_sd_y = np.full((num_times, total_particles), -1.0, dtype=np.float64)
    final_sd_z = np.full((num_times, total_particles), -1.0, dtype=np.float64)
    final_pre_sdid = np.full((num_times, total_particles), -1, dtype=np.int32)
    final_pre_dmid = np.full((num_times, total_particles), -1, dtype=np.int32)
    final_if_coal = np.full((num_times, total_particles), -1, dtype=np.int8)
    final_other_pre_dmid = np.full((num_times, total_particles), -1, dtype=np.int32)
    final_other_pre_sdid = np.full((num_times, total_particles), -1, dtype=np.int32)
    final_num_col = np.full((num_times, total_particles), -1, dtype=np.int32)
    final_other_sd_r = np.full((num_times, total_particles), -1.0, dtype=np.float64)
    final_other_sd_n = np.full((num_times, total_particles), -1, dtype=np.int64)
    # Allocate an array for particle_index (global id)
    final_particle_index = np.full((num_times, total_particles), -1, dtype=np.int64)

    # Compute column offsets for each block.
    offsets = [0] * num_blocks
    offsets[0] = 0
    for b in range(1, num_blocks):
        offsets[b] = offsets[b-1] + block_sizes[b-1]

    # For each time step, for each block, read the corresponding temporary file.
    for t_idx, t_str in enumerate(reversed(time_steps)):
        for b in range(num_blocks):
            temp_file = os.path.join(output_dir, f"temp_block_{b}_time_{t_str}.nc")
            if not os.path.exists(temp_file):
                print(f"[Warning] Temporary file not found: {temp_file}")
                continue
            with Dataset(temp_file, 'r') as nc:
                nparts = len(nc.dimensions['particle'])
                if nparts == 0:
                    continue
                # Read variables from the temporary file.
                sd_r_arr = nc.variables['sd_r'][:] if 'sd_r' in nc.variables else None
                sd_n_arr = nc.variables['sd_n'][:] if 'sd_n' in nc.variables else None
                sd_x_arr = nc.variables['sd_x'][:] if 'sd_x' in nc.variables else None
                sd_y_arr = nc.variables['sd_y'][:] if 'sd_y' in nc.variables else None
                sd_z_arr = nc.variables['sd_z'][:] if 'sd_z' in nc.variables else None
                pre_sdid_arr = read_tracking_var(nc, ['pre_sdid', 'sd_id'], ['pre_sdid_', 'sd_id_'])
                pre_dmid_arr = read_tracking_var(nc, ['pre_dmid', 'dm_id'], ['pre_dmid_', 'dm_id_'])
                if_coal_arr = read_tracking_var(nc, ['if_coal'], ['if_coal_'])
                other_pre_dmid_arr = nc.variables['other_pre_dmid'][:] if 'other_pre_dmid' in nc.variables else None
                other_pre_sdid_arr = nc.variables['other_pre_sdid'][:] if 'other_pre_sdid' in nc.variables else None
                num_col_arr = nc.variables['num_col'][:] if 'num_col' in nc.variables else None
                other_sd_r_arr = nc.variables['other_sd_r'][:] if 'other_sd_r' in nc.variables else None
                other_sd_n_arr = nc.variables['other_sd_n'][:] if 'other_sd_n' in nc.variables else None
                pindex_arr = nc.variables['particle_index'][:] if 'particle_index' in nc.variables else None

                start_col = offsets[b]
                end_col = start_col + nparts
                if sd_r_arr is not None:
                    final_sd_r[t_idx, start_col:end_col] = sd_r_arr
                if sd_n_arr is not None:
                    final_sd_n[t_idx, start_col:end_col] = sd_n_arr
                if sd_x_arr is not None:
                    final_sd_x[t_idx, start_col:end_col] = sd_x_arr
                if sd_y_arr is not None:
                    final_sd_y[t_idx, start_col:end_col] = sd_y_arr
                if sd_z_arr is not None:
                    final_sd_z[t_idx, start_col:end_col] = sd_z_arr
                if pre_sdid_arr is not None:
                    final_pre_sdid[t_idx, start_col:end_col] = pre_sdid_arr
                if pre_dmid_arr is not None:
                    final_pre_dmid[t_idx, start_col:end_col] = pre_dmid_arr
                if if_coal_arr is not None:
                    final_if_coal[t_idx, start_col:end_col] = if_coal_arr
                if other_pre_dmid_arr is not None:
                    final_other_pre_dmid[t_idx, start_col:end_col] = other_pre_dmid_arr
                if other_pre_sdid_arr is not None:
                    final_other_pre_sdid[t_idx, start_col:end_col] = other_pre_sdid_arr
                if num_col_arr is not None:
                    final_num_col[t_idx, start_col:end_col] = num_col_arr
                if other_sd_r_arr is not None:
                    final_other_sd_r[t_idx, start_col:end_col] = other_sd_r_arr
                if other_sd_n_arr is not None:
                    final_other_sd_n[t_idx, start_col:end_col] = other_sd_n_arr
                # Read particle_index from temporary file.
                if pindex_arr is not None:
                    final_particle_index[t_idx, start_col:end_col] = pindex_arr

        # Now, for each timestep, sort the data by particle_index.
        sort_idx = np.argsort(final_particle_index[t_idx, :])
        final_sd_r[t_idx, :] = final_sd_r[t_idx, :][sort_idx]
        final_sd_n[t_idx, :] = final_sd_n[t_idx, :][sort_idx]
        final_sd_x[t_idx, :] = final_sd_x[t_idx, :][sort_idx]
        final_sd_y[t_idx, :] = final_sd_y[t_idx, :][sort_idx]
        final_sd_z[t_idx, :] = final_sd_z[t_idx, :][sort_idx]
        final_pre_sdid[t_idx, :] = final_pre_sdid[t_idx, :][sort_idx]
        final_pre_dmid[t_idx, :] = final_pre_dmid[t_idx, :][sort_idx]
        final_if_coal[t_idx, :] = final_if_coal[t_idx, :][sort_idx]
        final_other_pre_dmid[t_idx, :] = final_other_pre_dmid[t_idx, :][sort_idx]
        final_other_pre_sdid[t_idx, :] = final_other_pre_sdid[t_idx, :][sort_idx]
        final_num_col[t_idx, :] = final_num_col[t_idx, :][sort_idx]
        final_other_sd_r[t_idx, :] = final_other_sd_r[t_idx, :][sort_idx]
        final_other_sd_n[t_idx, :] = final_other_sd_n[t_idx, :][sort_idx]
        final_particle_index[t_idx, :] = final_particle_index[t_idx, :][sort_idx]

        # Write final data for this time step into a NetCDF file.
        temp_file = os.path.join(output_dir, f"temp_block_merged_time_{t_str}.nc")
        with Dataset(temp_file, 'w', format='NETCDF4') as nc_out:
            nc_out.createDimension('particle', total_particles)  # Number of particles
            nc_out.createDimension('time', 1)  # Only 1 time step for each file

            # Create time variable and store the current time step
            time_var = nc_out.createVariable('time', str, ('time',), zlib=True, complevel=4)
            time_var[0] = t_str

            # Store the data for each variable
            var_names = ['particle_index', 'sd_r', 'sd_n', 'sd_x', 'sd_y', 'sd_z',
                         'pre_sdid', 'pre_dmid', 'if_coal', 'other_pre_dmid', 'other_pre_sdid',
                         'num_col', 'other_sd_r', 'other_sd_n']
            var_dtypes = [np.int64, np.float64, np.int64, np.float64, np.float64, np.float64,
                          np.int32, np.int32, np.int8, np.int32, np.int32, np.int32,
                          np.float64, np.int64]
            
            for var_name, dtype in zip(var_names, var_dtypes):
                var = nc_out.createVariable(var_name, dtype, ('particle',), zlib=True, complevel=4)
                var[:] = eval(f"final_{var_name}[t_idx, :]")

    print(f"[Info] Successfully merged and written data for all time steps into {output_file}.")

def split_particles(particle_infos, num_blocks):
    """
    Splits the list of particle_infos into num_blocks sublists for parallel processing.
    Returns:
      - blocks: A list of sublists.
      - block_sizes: A list of particle counts for each block.
    """
    num_sd = len(particle_infos)
    base_size = num_sd // num_blocks
    remainder = num_sd % num_blocks

    blocks = []
    block_sizes = []
    start = 0
    for i in range(num_blocks):
        size = base_size + 1 if i < remainder else base_size
        end = start + size
        blocks.append(particle_infos[start:end])
        block_sizes.append(size)
        start = end

    return blocks, block_sizes

def main(input_directory, output_dir, output_file):
    """
    Main driver function that:
      1. Generates time steps,
      2. Initializes particles (using parallel file scanning),
      3. Splits them into blocks,
      4. Processes each block in parallel,
      5. Merges the results into a single NetCDF file.
    """
    start_time = "011000.000"
    end_time = "010000.000"
    time_steps = generate_time_steps(start_time, end_time)

    # Initialize particles; returns part_array, mapping array, and particle_infos.
    part_array, map_arr, particle_infos = initialize_particles(input_directory, time_steps)
    
    num_times = len(time_steps)
    num_sd = len(particle_infos)
    num_blocks = 4  # Adjust as needed

    if num_blocks > num_sd:
        raise ValueError("num_blocks must not exceed the total number of particles.")

    # Split particle_infos into blocks and also get block sizes.
    particle_blocks, block_sizes = split_particles(particle_infos, num_blocks)
    
    # Use a Manager for shared missing_cache.
    with Manager() as manager:
        missing_cache = manager.list()  # Shared list for missing files
        
        # Process each block in parallel.
        with Pool(num_blocks) as pool:
            args = []
            for block_idx, block in enumerate(particle_blocks):
                # Each block is a list of (filename, idx). Pass global mapping (map_arr)
                args.append((block, input_directory, time_steps, num_times, output_dir, block_idx, missing_cache, map_arr))
            pool.starmap(process_single_particle_block, args)
        
        # After processing, build a list of all temporary file names.
        temp_files = []
        for b in range(num_blocks):
            for t in range(num_times):
                temp_files.append(os.path.join(output_dir, f"temp_block_{b}_time_{time_steps[t]}.nc"))
        
        # Merge temporary files into final output.
        merge_temp_files(output_file, output_dir, time_steps, num_blocks, block_sizes)
        
if __name__ == "__main__":
    input_directory = "../"
    output_dir = "./temp/"
    output_file = "./SD_011000.000.nc"
    os.makedirs(output_dir, exist_ok=True)
    main(input_directory, output_dir, output_file)
