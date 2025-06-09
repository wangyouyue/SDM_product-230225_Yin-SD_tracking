import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import random
from matplotlib import cm
from matplotlib.colors import LogNorm

# Configuration parameters
TIME_STEP = 0.1  # Time step in seconds
INITIAL_RADIUS_THRESHOLD = 1.0  # Initial radius threshold in microns
MAX_RADIUS_THRESHOLD = 5.0  # Maximum radius threshold in microns

def split_trajectory_at_boundary(x_trajectory, z_trajectory, boundary=400):
    """
    Splits the trajectory into multiple segments if the particle crosses the periodic boundary.
    """
    segments = []
    start_idx = 0

    # Loop through the x_trajectory to detect jumps greater than the boundary
    for i in range(1, len(x_trajectory)):
        # Check if the particle crosses the boundary
        if abs(x_trajectory[i] - x_trajectory[i - 1]) > boundary / 2:
            # If a jump is detected, split the trajectory at this point
            segments.append((x_trajectory[start_idx:i], z_trajectory[start_idx:i]))
            start_idx = i

    # Add the last segment (from the last split point to the end)
    segments.append((x_trajectory[start_idx:], z_trajectory[start_idx:]))

    return segments

def filter_valid_data(data, fill_value=-9999):
    """
    Filter out fill values from the data.
    """
    return data[data != fill_value]

def main():
    # 1) Open the NetCDF file
    file_name = "particle_trajectories.nc"  # Modify if the filename is different
    
    try:
        ds = nc.Dataset(file_name, "r")
    except FileNotFoundError:
        print(f"Error: File '{file_name}' not found. Please check the file path.")
        return
    except Exception as e:
        print(f"Error opening file: {e}")
        return

    # 2) Read the relevant variables
    try:
        # Note: Data dimensions are [particle_index, time]
        x_data = ds.variables["x_coord"][:]  # x position in meters
        z_data = ds.variables["z_coord"][:]  # z position in meters
        r_data = ds.variables["radius"][:] * 1e6  # Convert radius from meters to microns
        particle_index = ds.variables["sd_id_in_file_at_step"][:, 0]  # Read particle_index at time=0
    except KeyError as e:
        print(f"Error: Variable {e} not found in the NetCDF file.")
        ds.close()
        return

    # 3) Get the number of particles and time steps
    num_sd, time_len = x_data.shape  # First dimension is particles, second dimension is time
    print(f"Data dimensions: {num_sd} particles, {time_len} time steps")
    
    # 4) Use all time steps
    test_time_steps = time_len
    print(f"Using all {test_time_steps} time steps")
    print(f"Radius criteria: initial < {INITIAL_RADIUS_THRESHOLD} μm AND maximum > {MAX_RADIUS_THRESHOLD} μm")

    # 5) Find particles that meet the radius criteria
    valid_particles = []
    for i in range(num_sd):
        # Check if particle has valid position data
        valid_x = filter_valid_data(x_data[i, :])  # x coordinates for particle i across all time steps
        valid_z = filter_valid_data(z_data[i, :])  # z coordinates for particle i across all time steps
        valid_r = filter_valid_data(r_data[i, :])  # radius for particle i across all time steps
        
        # Require at least 2 valid points to draw a trajectory
        if len(valid_x) >= 2 and len(valid_z) >= 2 and len(valid_r) >= 2:
            # Check radius criteria
            initial_radius = valid_r[0]  # First valid radius
            max_radius = np.max(valid_r)  # Maximum radius during the trajectory
            
            if initial_radius < INITIAL_RADIUS_THRESHOLD and max_radius > MAX_RADIUS_THRESHOLD:
                valid_particles.append(i)
    
    print(f"Found {len(valid_particles)} particles meeting radius criteria")
    
    if len(valid_particles) == 0:
        print(f"No particles meeting the radius criteria found!")
        print(f"Criteria: initial radius < {INITIAL_RADIUS_THRESHOLD} μm AND maximum radius > {MAX_RADIUS_THRESHOLD} μm")
        ds.close()
        return

    # 6) Randomly select 10 particles from valid ones
    num_pick = min(10, len(valid_particles))
    picked_particles = random.sample(valid_particles, num_pick)
    print(f"Randomly selected {num_pick} particles for trajectory plotting")
    
    # Print selected particle IDs and their radius info
    selected_ids = [particle_index[p] for p in picked_particles]
    print(f"Selected particle IDs: {selected_ids}")
    
    for p_idx in picked_particles:
        valid_r = filter_valid_data(r_data[p_idx, :])
        if len(valid_r) > 0:
            print(f"Particle {particle_index[p_idx]}: initial radius = {valid_r[0]:.3f} μm, max radius = {np.max(valid_r):.3f} μm")

    # 7) Prepare color mapping based on radius range
    all_valid_radii = []
    for p_idx in picked_particles:
        valid_radii = filter_valid_data(r_data[p_idx, :])  # radius for particle p_idx across all time steps
        if len(valid_radii) > 0:
            all_valid_radii.extend(valid_radii)
    
    if len(all_valid_radii) > 0:
        all_valid_radii = np.array(all_valid_radii)
        min_radius = np.min(all_valid_radii)
        max_radius = np.max(all_valid_radii)
        print(f"Radius range for color mapping: {min_radius:.6f} - {max_radius:.6f} μm")
    else:
        print("Warning: No valid radius data found for color mapping.")
        min_radius = 1e-5
        max_radius = 1e-1

    # Use logarithmic normalization for the color mapping
    norm = LogNorm(vmin=max(min_radius, 1e-6), vmax=max_radius)
    cmap = cm.viridis

    # 8) Create the plot with larger fonts
    plt.rcParams.update({'font.size': 16})  # Increase global font size
    
    # Use consistent colors for each particle across both subplots
    colors = plt.cm.tab10(np.linspace(0, 1, num_pick))  # Different colors for each particle
    
    # 9) Create subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 14))
    
    # Top plot: Trajectories
    ax1.set_xlabel("X Position (m)", fontsize=18)
    ax1.set_ylabel("Z Position (m)", fontsize=18)
    ax1.set_title(f"Trajectories of {num_pick} Randomly Selected Particles ({test_time_steps} time steps)", fontsize=20)
    ax1.set_xlim(0, 400)
    
    # Set z-axis limit based on actual data
    all_z_data = z_data[z_data != -9999]
    if len(all_z_data) > 0:
        ax1.set_ylim(0, np.max(all_z_data) * 1.05)
    else:
        ax1.set_ylim(0, 1000)
    
    ax1.grid(True, alpha=0.3)
    ax1.tick_params(axis='both', which='major', labelsize=16)
    
    # Store legend handles for particles
    particle_handles = []
    particle_labels = []
    
    # Plot trajectories on top subplot
    for idx, p_idx in enumerate(picked_particles):
        x_trajectory = x_data[p_idx, :]  # x coordinates for particle p_idx across all time steps
        z_trajectory = z_data[p_idx, :]  # z coordinates for particle p_idx across all time steps
        
        valid_mask = (x_trajectory != -9999) & (z_trajectory != -9999)
        x_trajectory = x_trajectory[valid_mask]
        z_trajectory = z_trajectory[valid_mask]
        
        if len(x_trajectory) < 2:
            continue

        x_trajectory = np.mod(x_trajectory, 400)
        segments = split_trajectory_at_boundary(x_trajectory, z_trajectory)
        
        particle_color = colors[idx]
        particle_label = f"Particle {particle_index[p_idx]}"
        
        # Plot trajectory segments
        for seg_idx, segment in enumerate(segments):
            segment_x, segment_z = segment
            if len(segment_x) < 2:
                continue
            
            # Only add to legend for the first segment
            if seg_idx == 0:
                line, = ax1.plot(segment_x, segment_z, color=particle_color, linewidth=1, 
                        alpha=0.8, marker='o', markersize=1.5)
                particle_handles.append(line)
                particle_labels.append(particle_label)
            else:
                ax1.plot(segment_x, segment_z, color=particle_color, linewidth=1, 
                        alpha=0.8, marker='o', markersize=1.5)

    # Add start and end point markers (only once for legend)
    start_handle = None
    end_handle = None
    
    for idx, p_idx in enumerate(picked_particles):
        x_trajectory = x_data[p_idx, :]
        z_trajectory = z_data[p_idx, :]
        
        valid_mask = (x_trajectory != -9999) & (z_trajectory != -9999)
        x_trajectory = x_trajectory[valid_mask]
        z_trajectory = z_trajectory[valid_mask]
        
        if len(x_trajectory) > 0:
            x_trajectory = np.mod(x_trajectory, 400)
            
            # Add start and end points
            if start_handle is None:
                start_handle = ax1.scatter(x_trajectory[0], z_trajectory[0], color='blue', marker='o', 
                           s=80, edgecolors='white', linewidth=2, zorder=10)
            else:
                ax1.scatter(x_trajectory[0], z_trajectory[0], color='blue', marker='o', 
                           s=80, edgecolors='white', linewidth=2, zorder=10)
                           
            if end_handle is None:
                end_handle = ax1.scatter(x_trajectory[-1], z_trajectory[-1], color='red', marker='*', 
                           s=100, edgecolors='white', linewidth=1, zorder=10)
            else:
                ax1.scatter(x_trajectory[-1], z_trajectory[-1], color='red', marker='*', 
                           s=100, edgecolors='white', linewidth=1, zorder=10)
    
    # Add start/end point legend inside the first subplot
    if start_handle is not None and end_handle is not None:
        ax1.legend([start_handle, end_handle], ['Start Point', 'End Point'],
                  loc='upper right', fontsize=12,
                  title="Trajectory Markers", title_fontsize=14,
                  framealpha=0.9, fancybox=True, shadow=True)
    
    # Bottom plot: Radius evolution over time (converted to minutes)
    ax2.set_xlabel("Time (minutes)", fontsize=18)
    ax2.set_ylabel("Radius (μm)", fontsize=18)
    ax2.set_title("Radius Evolution Over Time", fontsize=20)
    ax2.grid(True, alpha=0.3)
    ax2.tick_params(axis='both', which='major', labelsize=16)
    
    time_steps = np.arange(test_time_steps)
    time_minutes = time_steps * TIME_STEP / 60.0  # Convert to minutes
    
    for idx, p_idx in enumerate(picked_particles):
        r_trajectory = r_data[p_idx, :]  # radius for particle p_idx across all time steps
        valid_mask = r_trajectory != -9999
        
        if np.sum(valid_mask) < 2:
            continue
            
        valid_time = time_minutes[valid_mask]
        valid_radius = r_trajectory[valid_mask]
        
        particle_color = colors[idx]  # Use same color as in top plot
        
        # Plot without label (will use shared legend)
        ax2.plot(valid_time, valid_radius, color=particle_color, linewidth=1, 
                alpha=0.8, marker='o', markersize=2)
    
    ax2.set_yscale('log')
    # Limit x-axis range to all time range in minutes
    ax2.set_xlim(0, (test_time_steps - 1) * TIME_STEP / 60.0)

    # Create shared legend for particles at the bottom of the figure
    fig.legend(particle_handles, particle_labels, 
              loc='lower center', bbox_to_anchor=(0.5, -0.02), 
              ncol=min(5, len(particle_labels)), fontsize=14,
              title="Particles", title_fontsize=16)

    # 10) Save and show the plot
    output_filename = "random_particle_trajectories.png"
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.15)  # Add space for bottom legend
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"\nPlot saved as '{output_filename}'")
    
    plt.show()
    
    # Close the dataset
    ds.close()
    print("Analysis completed successfully.")

if __name__ == "__main__":
    main()

