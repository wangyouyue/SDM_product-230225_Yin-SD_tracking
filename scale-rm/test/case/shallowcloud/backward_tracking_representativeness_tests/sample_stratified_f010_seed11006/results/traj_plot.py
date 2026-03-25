#!/usr/bin/env python3
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import random
from matplotlib import cm
from matplotlib.colors import LogNorm

def split_trajectory_at_boundary(x_trajectory, z_trajectory, boundary=6000):
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

def main():
    # 1) Open the NetCDF file
    file_name = "SD_011000.000.nc"  # Modify if the filename is different
    ds = nc.Dataset(file_name, "r")

    # 2) Read the relevant variables: 'sd_x', 'sd_z', 'sd_r', and 'particle_index'
    #    Shapes: [time, num_sd]
    x_data = ds.variables["sd_x"][:]  # x position in meters
    z_data = ds.variables["sd_z"][:]  # z position in meters
    r_data = ds.variables["sd_r"][:] * 1e6  # Convert radius from meters to microns (now in microns)
    particle_index = ds.variables["particle_index"][0, :]  # Read particle_index at time=0

    # 3) Get the time and number of super-droplets
    time_len, num_sd = x_data.shape

    # 4) Select particles that initially have radius < 5 microns but have reached 15 microns at any time
    selected_particles = []
    for i in range(num_sd):
        initial_radius = r_data[0, i]  # Radius in microns
        max_radius = np.max(r_data[:, i])  # Maximum radius during the simulation

        if initial_radius < 5.0 and max_radius >= 15.0:
            selected_particles.append(i)

    # 5) Randomly pick 10 particles from the selected ones (if there are fewer than 10, pick all)
    num_pick = min(10, len(selected_particles))  # In case there are fewer than 10 selected particles
    picked_particles = random.sample(selected_particles, num_pick)

    # 6) Plot the trajectories of the selected particles
    plt.figure(figsize=(10, 6))

    # Find the minimum and maximum radius values for the selected particles for color mapping
    selected_radii = r_data[:, selected_particles].flatten()  # Get all radii for the selected particles
    min_radius = np.min(selected_radii)
    max_radius = np.max(selected_radii)

    # Use logarithmic normalization for the color mapping (color by radius of selected particles)
    norm = LogNorm(vmin=min_radius, vmax=max_radius)
    cmap = cm.viridis  # Use 'viridis' color map (you can change this)

    for p_idx in picked_particles:
        # Periodic boundary correction for x position
        x_trajectory = x_data[:, p_idx]
        z_trajectory = z_data[:, p_idx]

        # Apply periodic boundary for x (wrap around at 6000 meters)
        x_trajectory = np.mod(x_trajectory, 6000)

        # Split the trajectory if there is a jump across the periodic boundary
        segments = split_trajectory_at_boundary(x_trajectory, z_trajectory)

        # Plot each segment separately
        for segment in segments:
            # Extract the radius values for the current segment
            segment_x, segment_z = segment
            start_idx = x_trajectory.tolist().index(segment_x[0])
            end_idx = start_idx + len(segment_x)  # Find corresponding radius values for the segment
            radii_segment = r_data[start_idx:end_idx, p_idx]  # Extract the radius values for this segment

            # Color by the radius (use the color map)
            scatter = plt.scatter(segment_x, segment_z, c=radii_segment, cmap=cmap, norm=norm, marker='o', s=0.5, label=f"Particle {particle_index[p_idx]}" if start_idx == 0 else "")

        # Mark the start with a smaller 'o' and the end with a smaller 'x'
        plt.scatter(x_trajectory[0], z_trajectory[0], color='b', marker='o', s=10)  # Start point (larger marker)
        plt.scatter(x_trajectory[-1], z_trajectory[-1], color='r', marker='x', s=10)  # End point (larger marker)

    # 7) Customize the plot
    plt.xlabel("X Position (m)", fontsize=14)  # Increase font size for x-axis label
    plt.ylabel("Z Position (m)", fontsize=14)  # Increase font size for y-axis label
    plt.title("Trajectories of Selected Particles", fontsize=16)  # Increase font size for title
    plt.xlim(0, 6000)  # x-coordinate range: 0 to 6000 meters
    plt.ylim(0, np.max(z_data))  # y-coordinate range: 0 to the max value of z_data
    plt.legend(loc='upper right')

    # 8) Add a color bar for the radius (in microns)
    cbar = plt.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=plt.gca())
    cbar.set_label("Radius (microns)", fontsize=12)

    # 9) Save the plot as a PNG file with 300 DPI
    plt.savefig("particle_trajectories.png", dpi=300)

    # 10) Show the plot (optional, you can comment out if not needed)
    plt.show()

if __name__ == "__main__":
    main()
