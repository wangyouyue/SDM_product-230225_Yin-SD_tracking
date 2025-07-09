import os
import glob
import numpy as np
from netCDF4 import Dataset
from multiprocessing import Pool, cpu_count
import sys
from decimal import Decimal, getcontext

# Set floating-point precision to avoid precision issues
getcontext().prec = 20

num_processes = 4  # Set the number of processes

def time_str_to_seconds(time_str):
    hours = int(time_str[0:2])
    minutes = int(time_str[2:4])
    seconds = Decimal(time_str[4:])
    total_seconds = Decimal(hours) * Decimal(3600) + Decimal(minutes) * Decimal(60) + seconds
    return total_seconds

def seconds_to_time_str(seconds):
    hours = int(seconds // Decimal(3600))
    minutes = int((seconds % Decimal(3600)) // Decimal(60))
    seconds = seconds % Decimal(60)

    if seconds >= Decimal('59.999999999999'):
        print(f"Warning: Seconds reached {seconds}, correcting to 0 and incrementing minutes.")
        seconds = Decimal('0.000')
        minutes += 1

    if minutes >= 60:
        print(f"Warning: Minutes reached {minutes}, correcting to 0 and incrementing hours.")
        hours += 1
        minutes -= 60

    time_str = f"{hours:02}{minutes:02}{seconds:06.3f}"
    print(f"Generated time string: {time_str} (from {seconds})")  # Debug output
    return time_str

def combine_1d_variables(file_paths, variable_name, time_str):
    all_data = []

    # Check if file_paths is empty
    if not file_paths:
        print(f"Error: No files found for variable '{variable_name}' at '{time_str}'.")
        sys.exit(1)  # Exit the program after outputting error message

    for file_path in file_paths:
        with Dataset(file_path, 'r') as dataset:
            if variable_name not in dataset.variables:
                print(f"Error: Variable '{variable_name}' not found in file {file_path}.")
                sys.exit(1)  # Exit the program after outputting error message

            variable_data = dataset.variables[variable_name][:]
            all_data.append(variable_data)

    if not all_data:
        print(f"Error: No data found for variable '{variable_name}' in files: {file_paths}")
        sys.exit(1)  # Exit the program after outputting error message

    combined_data = np.concatenate(all_data)
    return combined_data

def process_and_sort_variables(args):
    input_directory, time_str, variable_names = args
    file_paths = sorted(glob.glob(os.path.join(input_directory, f"SD_selected_NetCDF_00000101-{time_str}.pe*")))
    combined_vars = {}
    for var_name in variable_names:
        combined_vars[var_name] = combine_1d_variables(file_paths, var_name, time_str)
    sorted_indices = np.lexsort((combined_vars['sd_id'], combined_vars['dm_id']))
    for var_name in variable_names:
        combined_vars[var_name] = combined_vars[var_name][sorted_indices]
    return time_str, combined_vars

def find_and_update_missing_pairs(initial_pairs, current_pairs, initial_index_map, historical_missing_pairs):
    current_set = set(current_pairs)
    missing_pairs = [pair for pair in initial_pairs if pair not in current_set and pair not in historical_missing_pairs]

    # Update historically missing pairs
    historical_missing_pairs.update(missing_pairs)

    # Get indices of missing pairs
    missing_pairs_with_index = [(pair, initial_index_map[pair]) for pair in missing_pairs]

    return missing_pairs_with_index

def fill_missing_data(initial_pairs, missing_pairs_with_index, combined_vars, initial_vars, initial_index_map):
    filled_vars = {}
    for var_name in combined_vars:
        dtype = combined_vars[var_name].dtype
        if var_name in ['dm_id', 'sd_id', 'if_coal']:
            fill_value = -9999
        elif var_name == 'sd_n':
            fill_value = np.iinfo(np.int64).min
        else:
            fill_value = np.nan
        filled_vars[var_name] = np.full_like(initial_vars[var_name], fill_value, dtype=dtype)

    # Fill in existing data
    for var_name in combined_vars:
        for idx, pair in enumerate(zip(combined_vars['dm_id'], combined_vars['sd_id'])):
            original_idx = initial_index_map[pair]
            filled_vars[var_name][original_idx] = combined_vars[var_name][idx]

    # Complete missing data
    for pair, idx in missing_pairs_with_index:
        filled_vars['dm_id'][idx] = pair[0]
        filled_vars['sd_id'][idx] = pair[1]

    return filled_vars

def create_or_append_to_file(output_file, variable_names, combined_vars, time_step, time_seconds):
    if time_step == 0:
        with Dataset(output_file, 'w') as output_dataset:
            n_particles = len(combined_vars[variable_names[0]])
            output_dataset.createDimension('particles', n_particles)
            output_dataset.createDimension('time', None)  # Unlimited time dimension

            # Create time variable
            time_var = output_dataset.createVariable('time', 'f8', ('time',), zlib=True)
            time_var.units = 'seconds since start'

            for var_name in variable_names:
                if var_name in ['dm_id', 'sd_id', 'if_coal']:
                    dtype = 'i4'
                elif var_name == 'sd_n':
                    dtype = 'i8'
                else:
                    dtype = 'f8'
                output_dataset.createVariable(var_name, dtype, ('time', 'particles'), zlib=True)

    with Dataset(output_file, 'a') as output_dataset:
        # Append time variable
        time_var = output_dataset.variables['time']
        time_var[time_step] = time_seconds

        for var_name in variable_names:
            var = output_dataset.variables[var_name]
            var[time_step, :] = combined_vars[var_name]

def main():
    TIME_STEP_INTERVAL = Decimal('60.0')  # Time step (in seconds)
    time_str = "000000.000"
    end_time_str = "001000.000"
    input_directory = "../"  # Input directory
    variable_names = ['sd_x', 'sd_y', 'sd_z', 'sd_r', 'sd_n', 'dm_id', 'sd_id', 'if_coal']
    output_file = 'combined_sorted_output.nc'
    current_time_str = time_str
    current_time_seconds = time_str_to_seconds(current_time_str)
    end_time_seconds = time_str_to_seconds(end_time_str)
    time_step = 0

    # Get variables for the initial time and sort them
    initial_vars = process_and_sort_variables((input_directory, current_time_str, variable_names))[1]
    initial_dm_id = initial_vars['dm_id']
    initial_sd_id = initial_vars['sd_id']

    # Save (dm_id, sd_id) pairs for the initial time and their indices
    initial_pairs = list(zip(initial_dm_id, initial_sd_id))
    initial_index_map = {pair: idx for idx, pair in enumerate(initial_pairs)}

    # Record historically missing pairs
    historical_missing_pairs = set()

    # Write data for the initial time to the file
    create_or_append_to_file(output_file, variable_names, initial_vars, time_step, current_time_seconds)
    time_step += 1

    current_time_seconds = round(current_time_seconds + TIME_STEP_INTERVAL, 3)
    current_time_str = seconds_to_time_str(current_time_seconds)

    # Create task list
    tasks = []
    while current_time_seconds <= end_time_seconds:
        print(f"Current time string: {current_time_str}, current seconds: {current_time_seconds}")  # Debug output
        tasks.append((input_directory, current_time_str, variable_names))
        current_time_seconds += TIME_STEP_INTERVAL
        current_time_str = seconds_to_time_str(current_time_seconds)

    # Process each time step in parallel
    with Pool(min(num_processes, cpu_count())) as pool:
        results = pool.map(process_and_sort_variables, tasks)

    # Sort results in time order
    results.sort(key=lambda x: time_str_to_seconds(x[0]))

    # Process results and write to file
    for time_str, sorted_vars in results:
        current_pairs = list(zip(sorted_vars['dm_id'], sorted_vars['sd_id']))
        current_time_seconds = time_str_to_seconds(time_str)

        # Find and update missing pairs
        if len(current_pairs) < len(initial_pairs):
            missing_pairs_with_index = find_and_update_missing_pairs(initial_pairs, current_pairs, initial_index_map, historical_missing_pairs)
            sorted_vars = fill_missing_data(initial_pairs, missing_pairs_with_index, sorted_vars, initial_vars, initial_index_map)

        create_or_append_to_file(output_file, variable_names, sorted_vars, time_step, current_time_seconds)
        time_step += 1

    print(f"Process finished, data has been saved in {output_file}")

if __name__ == "__main__":
    main()