import os
import glob
import numpy as np
from netCDF4 import Dataset
from datetime import datetime, timedelta
from multiprocessing import Pool

DX, DY, DZ = 50, 50, 5  # Grid dimensions
TIME_STEP_INTERVAL = 60000  # Time step interval in milliseconds

def initialize_netcdf(output_file, time_steps, num_sd):
    with Dataset(output_file, 'w', format='NETCDF4') as nc:
        nc.createDimension('time', len(time_steps))
        nc.createDimension('num_sd', num_sd)

        time_var = nc.createVariable('time', str, ('time',))
        time_var[:] = time_steps

        for key in ['x_idx', 'y_idx', 'z_idx', 'sd_r', 'sd_n', 'pre_sdid', 'pre_dmid', 'if_coal', 'other_pre_dmid', 'other_pre_sdid', 'num_col', 'other_sd_r', 'other_sd_n']:
            dtype = 'f4'
            if key in ['sd_n', 'other_sd_n']:
                dtype = 'i8'
            elif key in ['pre_sdid', 'pre_dmid', 'if_coal', 'other_pre_dmid', 'other_pre_sdid', 'num_col']:
                dtype = 'i4'
            elif key in ['if_coal']:
                dtype = 'i2'
            nc.createVariable(key, dtype, ('time', 'num_sd'))

def generate_time_steps(start_time, end_time, step_ms=100):
    time_format = "%H%M%S.%f"
    current_time = datetime.strptime(start_time, time_format)
    end_time = datetime.strptime(end_time, time_format)
    time_delta = timedelta(milliseconds=step_ms)

    time_steps = []
    while current_time >= end_time:
        time_steps.append(current_time.strftime(time_format)[:-3])  # Remove microseconds
        current_time -= time_delta

    return time_steps

def get_previous_time(current_time):
    try:
        hours = int(current_time[:2])
        minutes = int(current_time[2:4])
        seconds = int(current_time[4:6])
        millis = int(current_time[7:])

        total_millis = hours * 3600000 + minutes * 60000 + seconds * 1000 + millis - TIME_STEP_INTERVAL
        new_hours = (total_millis // 3600000) % 24
        new_minutes = (total_millis % 3600000) // 60000
        new_seconds = (total_millis % 60000) // 1000
        new_millis = total_millis % 1000

        return f"{new_hours:02}{new_minutes:02}{new_seconds:02}.{new_millis:03}"
    except Exception as e:
        print(f"Error calculating previous time from {current_time}: {e}")
        return "000000.000"

def generate_coal_filename(time, domain_id, input_directory):
    domain_id_str = f"{domain_id:06d}"
    return os.path.join(input_directory, f"SD_coal_output_NetCDF_00000101-{time}.pe{domain_id_str}")

def generate_filename(time, domain_id, input_directory):
    domain_id_str = f"{domain_id:06d}"
    return os.path.join(input_directory, f"SD_all_NetCDF_00000101-{time}.pe{domain_id_str}")

def get_particle_data(filename, sd_id):
    try:
        with Dataset(filename, 'r') as nc:
            sd_r = nc.variables['sd_r'][sd_id - 1]
            sd_n = nc.variables['sd_n'][sd_id - 1]
            sd_x = nc.variables['sd_x'][sd_id - 1]
            sd_y = nc.variables['sd_y'][sd_id - 1]
            sd_z = nc.variables['sd_z'][sd_id - 1]
            pre_sdid = nc.variables['pre_sdid'][sd_id - 1]
            pre_dmid = nc.variables['pre_dmid'][sd_id - 1]
            if_coal = nc.variables['if_coal'][sd_id - 1]
            sd_x = np.floor(sd_x / DX) * DX + DX / 2
            sd_y = np.floor(sd_y / DY) * DY + DY / 2
            sd_z = np.floor(sd_z / DZ) * DZ + DZ / 2
            data = {
                'sd_r': sd_r,
                'sd_n': sd_n,
                'sd_x': sd_x,
                'sd_y': sd_y,
                'sd_z': sd_z,
                'pre_sdid': pre_sdid,
                'pre_dmid': pre_dmid,
                'if_coal': if_coal
            }
            return data
    except Exception as e:
        print(f"Error processing {filename}: {e}")
    return None

def check_coalescence(coal_file, pre_dmid, pre_sdid, current_file):
    try:
        with Dataset(coal_file, 'r') as nc:
            pre_sdid1 = nc.variables['pre_sdid1'][:]
            pre_dmid1 = nc.variables['pre_dmid1'][:]
            pre_sdid2 = nc.variables['pre_sdid2'][:]
            pre_dmid2 = nc.variables['pre_dmid2'][:]
            num_col = nc.variables['num_col'][:]

            mask1 = (pre_dmid1 == pre_dmid) & (pre_sdid1 == pre_sdid)
            mask2 = (pre_dmid2 == pre_dmid) & (pre_sdid2 == pre_sdid)

            if np.any(mask1):
                index = np.where(mask1)[0][0]
                return pre_dmid2[index], pre_sdid2[index], num_col[index]
            elif np.any(mask2):
                index = np.where(mask2)[0][0]
                return pre_dmid1[index], pre_sdid1[index], num_col[index]
            else:
                print(f"No matching coalescence found in {coal_file} for (pre_dmid, pre_sdid) = ({pre_dmid}, {pre_sdid}) in file {current_file}")
                exit()
    except Exception as e:
        print(f"Error processing {coal_file} for (pre_dmid, pre_sdid) = ({pre_dmid}, {pre_sdid}): {e}")
    return None, None, None

def initialize_particles(input_directory, time_steps):
    initial_time_step = time_steps[0]
    file_pattern = os.path.join(input_directory, f"SD_all_NetCDF_00000101-{initial_time_step}.pe*")
    files = glob.glob(file_pattern)

    particle_infos = []
    domain_ids = []

    for filename in files:
        domain_id = int(filename.split(".pe")[1])
        with Dataset(filename, 'r') as nc:
            # sd_r = nc.variables['sd_r'][:]
            sd_z = nc.variables['sd_z'][:]
            mask = sd_z > 0
            indices = np.where(mask)[0]
            for idx in indices:
                particle_infos.append((filename, idx))
                domain_ids.append(domain_id)

    if not particle_infos:
        print("No particle found satisfying the conditions.")
        exit()

    # particle_infos = particle_infos[:2000]  # Select only 2000 particles for testing
    # domain_ids = domain_ids[:2000]  # Select only 2000 domain_ids for testing

    num_times = len(time_steps)
    num_particles = len(particle_infos)
    particles = {
        'x_idx': np.full((num_times, num_particles), -1.0),
        'y_idx': np.full((num_times, num_particles), -1.0),
        'z_idx': np.full((num_times, num_particles), -1.0),
        'sd_r': np.full((num_times, num_particles), -1.0),
        'sd_n': np.full((num_times, num_particles), -1, dtype=np.int64),
        'pre_sdid': np.full((num_times, num_particles), -1, dtype=np.int32),
        'pre_dmid': np.full((num_times, num_particles), -1, dtype=np.int32),
        'if_coal': np.full((num_times, num_particles), -1, dtype=np.int8),
        'other_pre_dmid': np.full((num_times, num_particles), -1, dtype=np.int32),
        'other_pre_sdid': np.full((num_times, num_particles), -1, dtype=np.int32),
        'num_col': np.full((num_times, num_particles), -1, dtype=np.int32),
        'other_sd_r': np.full((num_times, num_particles), -1.0),
        'other_sd_n': np.full((num_times, num_particles), -1, dtype=np.int64),
        'domain_id': np.full((num_times, num_particles), -1, dtype=np.int32)
    }

    for particle_index, ((filename, idx), domain_id) in enumerate(zip(particle_infos, domain_ids)):
        try:
            with Dataset(filename, 'r') as nc:
                sd_x = nc.variables['sd_x'][idx]
                sd_y = nc.variables['sd_y'][idx]
                sd_z_filtered = nc.variables['sd_z'][idx]
                sd_r_filtered = nc.variables['sd_r'][idx]
                sd_n = nc.variables['sd_n'][idx]
                pre_sdid = nc.variables['pre_sdid'][idx]
                pre_dmid = nc.variables['pre_dmid'][idx]
                if_coal = nc.variables['if_coal'][idx]

                x_idx = np.floor(sd_x / DX) * DX + DX / 2
                y_idx = np.floor(sd_y / DY) * DY + DY / 2
                z_idx = np.floor(sd_z_filtered / DZ) * DZ + DZ / 2

                particles['x_idx'][0, particle_index] = x_idx
                particles['y_idx'][0, particle_index] = y_idx
                particles['z_idx'][0, particle_index] = z_idx
                particles['sd_r'][0, particle_index] = sd_r_filtered
                particles['sd_n'][0, particle_index] = sd_n
                particles['pre_sdid'][0, particle_index] = pre_sdid
                particles['pre_dmid'][0, particle_index] = pre_dmid
                particles['if_coal'][0, particle_index] = if_coal
                particles['domain_id'][0, particle_index] = domain_id
        except Exception as e:
            print(f"Error processing {filename}: {e}")

    return particles, particle_infos

def process_single_particle_block(particle_block, input_directory, time_steps, num_times, output_dir, block_idx):
    for particle_index, (filename, idx) in enumerate(particle_block):
        particles = {
            'x_idx': np.full(num_times, -1.0),
            'y_idx': np.full(num_times, -1.0),
            'z_idx': np.full(num_times, -1.0),
            'sd_r': np.full(num_times, -1.0),
            'sd_n': np.full(num_times, -1, dtype=np.int64),
            'pre_sdid': np.full(num_times, -1, dtype=np.int32),
            'pre_dmid': np.full(num_times, -1, dtype=np.int32),
            'if_coal': np.full(num_times, -1, dtype=np.int8),
            'other_pre_dmid': np.full(num_times, -1, dtype=np.int32),
            'other_pre_sdid': np.full(num_times, -1, dtype=np.int32),
            'num_col': np.full(num_times, -1, dtype=np.int32),
            'other_sd_r': np.full(num_times, -1.0),
            'other_sd_n': np.full(num_times, -1, dtype=np.int64),
            'domain_id': np.full(num_times, -1, dtype=np.int32)
        }

        try:
            with Dataset(filename, 'r') as nc:
                sd_x = nc.variables['sd_x'][idx]
                sd_y = nc.variables['sd_y'][idx]
                sd_z_filtered = nc.variables['sd_z'][idx]
                sd_r_filtered = nc.variables['sd_r'][idx]
                sd_n = nc.variables['sd_n'][idx]
                pre_sdid = nc.variables['pre_sdid'][idx]
                pre_dmid = nc.variables['pre_dmid'][idx]
                if_coal = nc.variables['if_coal'][idx]

                x_idx = np.floor(sd_x / DX) * DX + DX / 2
                y_idx = np.floor(sd_y / DY) * DY + DY / 2
                z_idx = np.floor(sd_z_filtered / DZ) * DZ + DZ / 2

                particles['x_idx'][0] = x_idx
                particles['y_idx'][0] = y_idx
                particles['z_idx'][0] = z_idx
                particles['sd_r'][0] = sd_r_filtered
                particles['sd_n'][0] = sd_n
                particles['pre_sdid'][0] = pre_sdid
                particles['pre_dmid'][0] = pre_dmid
                particles['if_coal'][0] = if_coal
                particles['domain_id'][0] = int(filename.split(".pe")[1])
        except Exception as e:
            print(f"Error processing {filename} in process {block_idx}: {e}")

        for t in range(1, num_times):
            current_time_step = time_steps[t]
            previous_sd_id = particles['pre_sdid'][t-1]
            previous_dmid = particles['pre_dmid'][t-1]
            previous_domain_id = particles['domain_id'][t-1]
            current_file = generate_filename(current_time_step, previous_dmid, input_directory)
            particle_data = get_particle_data(current_file, previous_sd_id)

            if particle_data:
                particles['x_idx'][t] = particle_data['sd_x']
                particles['y_idx'][t] = particle_data['sd_y']
                particles['z_idx'][t] = particle_data['sd_z']
                particles['sd_r'][t] = particle_data['sd_r']
                particles['sd_n'][t] = particle_data['sd_n']
                particles['pre_sdid'][t] = particle_data['pre_sdid']
                particles['pre_dmid'][t] = particle_data['pre_dmid']
                particles['if_coal'][t] = particle_data['if_coal']
                particles['domain_id'][t] = previous_dmid

                if particles['if_coal'][t-1] == 1:
                    coal_time = current_time_step
                    coal_file = generate_coal_filename(coal_time, previous_domain_id, input_directory)
                    other_pre_dmid, other_pre_sdid, num_col = check_coalescence(coal_file, previous_dmid, previous_sd_id, current_file)

                    if other_pre_dmid is not None:
                        other_file = generate_filename(coal_time, other_pre_dmid, input_directory)
                        other_particle_data = get_particle_data(other_file, other_pre_sdid)

                        if other_particle_data:
                            if (particle_data['sd_r'] < other_particle_data['sd_r'] & np.max(particle_data['sd_r'], other_particle_data['sd_r']) < particles['sd_r'][t-1]):
                                particles['other_pre_dmid'][t] = particle_data['pre_dmid']
                                particles['other_pre_sdid'][t] = particle_data['pre_sdid']
                                particles['pre_dmid'][t], particles['pre_sdid'][t] = other_particle_data['pre_dmid'], other_particle_data['pre_sdid']

                                particles['other_sd_r'][t] = particle_data['sd_r']
                                particles['other_sd_n'][t] = particle_data['sd_n']
                                particles['sd_r'][t], particles['sd_n'][t] = other_particle_data['sd_r'], other_particle_data['sd_n']
                                particles['if_coal'][t] = other_particle_data['if_coal']
                                particles['x_idx'][t] = other_particle_data['sd_x']
                                particles['y_idx'][t] = other_particle_data['sd_y']
                                particles['z_idx'][t] = other_particle_data['sd_z']
                                particles['domain_id'][t] = other_pre_dmid
                            else:
                                particles['other_pre_dmid'][t] = other_particle_data['pre_dmid']
                                particles['other_pre_sdid'][t] = other_particle_data['pre_sdid']
                                particles['other_sd_r'][t] = other_particle_data['sd_r']
                                particles['other_sd_n'][t] = other_particle_data['sd_n']
                            particles['num_col'][t] = num_col
                        else:
                            print(f"No matching coalescence found in {coal_file} for (pre_dmid, pre_sdid) = ({previous_dmid}, {previous_sd_id}) in file {current_file} in process {block_idx}")
                            particles['other_pre_dmid'][t] = -1
                            particles['other_pre_sdid'][t] = -1
                            particles['num_col'][t] = -1
                            particles['other_sd_r'][t] = -1.0
                            particles['other_sd_n'][t] = -1
                    else:
                        print(f"No matching coalescence found in {coal_file} for (pre_dmid, pre_sdid) = ({previous_dmid}, {previous_sd_id}) in file {current_file} in process {block_idx}")
                        particles['other_pre_dmid'][t] = -1
                        particles['other_pre_sdid'][t] = -1
                        particles['num_col'][t] = -1
                        particles['other_sd_r'][t] = -1.0
                        particles['other_sd_n'][t] = -1
                else:
                    particles['other_pre_dmid'][t] = -1
                    particles['other_pre_sdid'][t] = -1
                    particles['num_col'][t] = -1
                    particles['other_sd_r'][t] = -1.0
                    particles['other_sd_n'][t] = -1
            else:
                print(f"Particle with pre_dmid={previous_dmid} and pre_sdid={previous_sd_id} not found in {current_file} in process {block_idx}")

        temp_file = os.path.join(output_dir, f"temp_block_{block_idx}_{idx}.nc")
        with Dataset(temp_file, 'w', format='NETCDF4') as nc:
            nc.createDimension('time', num_times)
            for key in particles:
                if key != 'domain_id':  # do not save domain_id
                    dtype = particles[key].dtype
                    var = nc.createVariable(key, dtype, ('time',))
                    var[:] = particles[key]

def write_to_netcdf(output_file, combined_data, time_steps):
    num_sd = combined_data['sd_r'].shape[1]  # get number of sds

    with Dataset(output_file, 'w', format='NETCDF4') as nc:
        nc.createDimension('time', len(time_steps))
        nc.createDimension('num_sd', num_sd)

        time_var = nc.createVariable('time', str, ('time',))
        for i, time in enumerate(time_steps):
            time_var[i] = time

        for key in combined_data:
            dtype = combined_data[key].dtype
            var = nc.createVariable(key, dtype, ('time', 'num_sd'))
            var[:, :] = combined_data[key]

def merge_temp_files(output_file, temp_files, time_steps, num_sd):
    num_times = len(time_steps)

    combined_data = {
        'x_idx': np.full((num_times, num_sd), -1.0),
        'y_idx': np.full((num_times, num_sd), -1.0),
        'z_idx': np.full((num_times, num_sd), -1.0),
        'sd_r': np.full((num_times, num_sd), -1.0),
        'sd_n': np.full((num_times, num_sd), -1, dtype=np.int64),
        'pre_sdid': np.full((num_times, num_sd), -1, dtype=np.int32),
        'pre_dmid': np.full((num_times, num_sd), -1, dtype=np.int32),
        'if_coal': np.full((num_times, num_sd), -1, dtype=np.int8),
        'other_pre_dmid': np.full((num_times, num_sd), -1, dtype=np.int32),
        'other_pre_sdid': np.full((num_times, num_sd), -1, dtype=np.int32),
        'num_col': np.full((num_times, num_sd), -1, dtype=np.int32),
        'other_sd_r': np.full((num_times, num_sd), -1.0),
        'other_sd_n': np.full((num_times, num_sd), -1, dtype=np.int64),
    }

    particle_idx = 0
    for temp_file in temp_files:
        with Dataset(temp_file, 'r') as nc:
            for key in combined_data:
                combined_data[key][:, particle_idx] = nc.variables[key][:]
        particle_idx += 1

    write_to_netcdf(output_file, combined_data, time_steps)

    # for temp_file in temp_files:
    #     os.remove(temp_file)

def split_particles(particle_infos, num_blocks):
    num_sd = len(particle_infos)
    base_size = num_sd // num_blocks  # size of every block
    remainder = num_sd % num_blocks   # number of the remaining sds

    blocks = []
    start = 0

    for i in range(num_blocks):
        # for the first {remainder} blocks, size of every block is {base_size} + 1
        # for other blocks, size is {base_size}
        size = base_size + 1 if i < remainder else base_size
        end = start + size
        blocks.append(particle_infos[start:end])
        start = end

    return blocks

def main(input_directory, output_dir, output_file):
    start_time = "001000.000"
    end_time = "000000.000"
    time_steps = generate_time_steps(start_time, end_time)

    particles, particle_infos = initialize_particles(input_directory, time_steps)

    num_times = len(time_steps)
    num_sd = len(particle_infos)
    num_blocks = 40  # Number of blocks for parallel processing

    if num_blocks > num_sd:
       raise ValueError("Number of blocks (num_blocks) cannot be greater than the number of particles (num_sd)")

    particle_blocks = split_particles(particle_infos, num_blocks)

    with Pool(num_blocks) as pool:
        args = [(particle_block, input_directory, time_steps, num_times, output_dir, block_idx) for block_idx, particle_block in enumerate(particle_blocks)]
        pool.starmap(process_single_particle_block, args)

    # Collect all temporary files generated
    temp_files = []
    for block_idx in range(num_blocks):
        for particle_index in range(len(particle_blocks[block_idx])):
            temp_files.append(os.path.join(output_dir, f"temp_block_{block_idx}_{particle_blocks[block_idx][particle_index][1]}.nc"))

    merge_temp_files(output_file, temp_files, time_steps, num_sd)

if __name__ == "__main__":
    input_directory = "../"  # Replace with the path to your NetCDF files
    output_dir = "./temp/"  # Directory to store temporary files
    output_file = "./SD_001000.000.nc"  # Output NetCDF file
    os.makedirs(output_dir, exist_ok=True)
    main(input_directory, output_dir, output_file)
