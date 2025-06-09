PROGRAM ParticleTracer
  ! Module area: contains NetCDF interface, custom types, global parameters
  USE netcdf
  IMPLICIT NONE

  ! ----------------------------------------------------------------------
  ! Parameter definition area (user configurable)
  ! ----------------------------------------------------------------------
  CHARACTER(LEN=256) :: path_to_data = "./"  ! Path to data files
  CHARACTER(LEN=64)  :: sd_all_prefix = "SD_all_NetCDF_"
  CHARACTER(LEN=64)  :: sd_coal_prefix = "SD_coal_output_NetCDF_"
  CHARACTER(LEN=256) :: output_nc_filename = "particle_trajectories.nc"
  CHARACTER(LEN=8)   :: date_string = "00000101" ! Date part of the filename

  CHARACTER(LEN=10)  :: time_start_str = "000000.000" ! Tracking start time (earliest time)
  CHARACTER(LEN=10)  :: time_end_str   = "001000.000" ! MODIFIED: Tracking end time for testing
  INTEGER            :: time_interval_milliseconds = 100 
  REAL(KIND=8)       :: filter_radius_gt = 1.0E-15_8 ! sd_r filter threshold (m)
  REAL(KIND=8)       :: filter_height_gt = 0.0_8     ! sd_z filter threshold (m)

  INTEGER            :: PRC_NUM_X = 4 ! Example value, number of domains in X direction
  INTEGER            :: PRC_NUM_Y = 1 ! Example value, number of domains in Y direction
  
  INTEGER, PARAMETER :: TRAJECTORY_BATCH_SIZE = 500 

  INTEGER, PARAMETER :: FILENAME_MAX_LEN = 256
  INTEGER, PARAMETER :: TIME_STR_LEN = 10
  INTEGER, PARAMETER :: DOMAIN_ID_STR_LEN = 6
  INTEGER, PARAMETER :: SAFE_DIM_NAME_LEN = 256 

  ! NetCDF fill values
  REAL(KIND=8), PARAMETER :: FILL_DOUBLE = -9999.0_8
  REAL(KIND=4), PARAMETER :: FILL_FLOAT  = -9999.0_4
  INTEGER, PARAMETER      :: FILL_INT    = -9999
  INTEGER(KIND=8), PARAMETER :: FILL_INT64 = -9999_8
  INTEGER(KIND=2), PARAMETER :: FILL_SHORT = -9999_2
  INTEGER(KIND=1), PARAMETER :: FILL_BYTE = -127

  ! ----------------------------------------------------------------------
  ! Type definitions
  ! ----------------------------------------------------------------------
  TYPE TimeRecord
     INTEGER :: hours, minutes, seconds, milliseconds
     INTEGER :: total_milliseconds_since_midnight 
  END TYPE TimeRecord

  TYPE ParticleStepData
     LOGICAL :: is_valid = .FALSE.
     CHARACTER(LEN=TIME_STR_LEN) :: timestamp_str
     INTEGER :: domain_id = FILL_INT
     INTEGER :: sd_id_in_file = FILL_INT 
     REAL(KIND=8) :: x_coord = FILL_DOUBLE
     REAL(KIND=8) :: y_coord = FILL_DOUBLE
     REAL(KIND=8) :: z_coord = FILL_DOUBLE
     REAL(KIND=8) :: radius = FILL_DOUBLE
     INTEGER(KIND=8) :: multiplicity = FILL_INT64
     INTEGER(KIND=2) :: if_coal_flag = FILL_SHORT 
     INTEGER :: parent_sdid_for_trace = FILL_INT 
     INTEGER :: parent_dmid_for_trace = FILL_INT
     INTEGER :: coalescence_num_col = FILL_INT
     INTEGER :: other_parent_prev_dmid = FILL_INT
     INTEGER :: other_parent_prev_sdid = FILL_INT 
     REAL(KIND=8) :: other_parent_radius = FILL_DOUBLE
     INTEGER(KIND=8) :: other_parent_multiplicity = FILL_INT64
  END TYPE ParticleStepData

  TYPE Trajectory
     INTEGER :: trajectory_id_global = FILL_INT 
     CHARACTER(LEN=TIME_STR_LEN) :: initial_selection_time_str
     INTEGER :: initial_selection_domain_id = FILL_INT
     INTEGER :: initial_selection_sd_id = FILL_INT 
     INTEGER :: actual_length = 0
     TYPE(ParticleStepData), ALLOCATABLE :: steps(:) 
  END TYPE Trajectory

  TYPE SeedParticle
      CHARACTER(LEN=TIME_STR_LEN) :: time_str
      INTEGER :: domain_id
      INTEGER :: sd_id 
  END TYPE SeedParticle

  TYPE ActiveTrajectoryState
      INTEGER :: original_seed_list_idx     
      INTEGER :: global_trajectory_array_idx 
      LOGICAL :: is_active = .TRUE.
      TYPE(TimeRecord) :: current_time_rec ! Time of the particle *being processed* in this step
      INTEGER :: current_domain_id = FILL_INT
      INTEGER :: current_sd_id_0based = FILL_INT      
      INTEGER :: next_parent_dmid_to_find = FILL_INT
      INTEGER :: next_parent_sdid_1based_to_find = FILL_INT 
      INTEGER :: steps_recorded_this_trajectory = 0
  END TYPE ActiveTrajectoryState

  TYPE CachedParticleInfo 
      REAL(KIND=8) :: x_coord, y_coord, z_coord, radius
      INTEGER(KIND=8) :: multiplicity
      INTEGER(KIND=2) :: if_coal_flag
      INTEGER :: pre_sdid_1based, pre_dmid
  END TYPE CachedParticleInfo

  TYPE CachedSDAllFile
      CHARACTER(LEN=FILENAME_MAX_LEN) :: filename
      TYPE(TimeRecord) :: file_time_rec
      INTEGER :: domain_id = FILL_INT
      LOGICAL :: is_loaded = .FALSE.
      INTEGER :: num_particles = 0
      TYPE(CachedParticleInfo), ALLOCATABLE :: particles(:)
  END TYPE CachedSDAllFile

  TYPE CachedCoalPairInfo 
      INTEGER :: pre_sdid1_1based, pre_dmid1
      INTEGER :: pre_sdid2_1based, pre_dmid2
      INTEGER :: num_col
  END TYPE CachedCoalPairInfo

  TYPE CachedSDCoalFile
      CHARACTER(LEN=FILENAME_MAX_LEN) :: filename
      TYPE(TimeRecord) :: file_time_rec 
      INTEGER :: domain_id = FILL_INT    
      LOGICAL :: is_loaded = .FALSE.
      INTEGER :: num_pairs = 0
      TYPE(CachedCoalPairInfo), ALLOCATABLE :: coal_pairs(:)
  END TYPE CachedSDCoalFile

  ! ----------------------------------------------------------------------
  ! Global variables
  ! ----------------------------------------------------------------------
  TYPE(TimeRecord) :: time_start_rec, time_end_rec
  INTEGER :: total_domains
  INTEGER :: num_time_steps_possible 

  TYPE(SeedParticle), ALLOCATABLE :: seed_particles_list(:) 
  INTEGER :: num_seed_particles_found = 0

  TYPE(Trajectory), ALLOCATABLE :: all_trajectories(:) 
  
  INTEGER :: ncid_out, err
  INTEGER :: dimid_traj, dimid_time_step, dimid_char_time
  INTEGER :: varid_init_sel_time, varid_init_sel_dom, varid_init_sel_sdid
  INTEGER :: varid_traj_len
  INTEGER :: varid_is_valid_step, varid_ts_at_step, varid_dom_at_step, varid_sdid_at_step
  INTEGER :: varid_x, varid_y, varid_z, varid_r, varid_n, varid_if_coal
  INTEGER :: varid_parent_sdid, varid_parent_dmid
  INTEGER :: varid_coal_num_col, varid_other_dmid, varid_other_sdid
  INTEGER :: varid_other_r, varid_other_n

  ! ----------------------------------------------------------------------
  ! Main program logic starts
  ! ----------------------------------------------------------------------
  PRINT *, "Starting Particle Tracer Program..."

  total_domains = PRC_NUM_X * PRC_NUM_Y
  CALL parse_time_string(time_start_str, time_start_rec)
  CALL parse_time_string(time_end_str, time_end_rec)

  IF (time_end_rec%total_milliseconds_since_midnight < time_start_rec%total_milliseconds_since_midnight .OR. &
      time_interval_milliseconds <= 0) THEN
      PRINT *, "Error: Invalid time configuration."
      STOP "Invalid time configuration"
  END IF
  num_time_steps_possible = (time_end_rec%total_milliseconds_since_midnight - &
                             time_start_rec%total_milliseconds_since_midnight) / time_interval_milliseconds + 1
  IF (num_time_steps_possible <= 0) THEN
      PRINT *, "Error: num_time_steps_possible calculated to be non-positive: ", num_time_steps_possible
      STOP "Invalid num_time_steps_possible"
  END IF

  PRINT *, "Total domains to process: ", total_domains
  PRINT *, "Max possible trajectory length (num_time_steps_possible): ", num_time_steps_possible

  PRINT *, "Step 1: Identifying seed particles at time ", time_end_str, " ..."
  CALL find_seed_particles(time_end_rec)
  PRINT *, "Found ", num_seed_particles_found, " seed particles."

  IF (num_seed_particles_found == 0) THEN
    PRINT *, "No seed particles found. Exiting."
    STOP
  END IF

  ALLOCATE(all_trajectories(num_seed_particles_found))

  PRINT *, "Step 2: Initializing output NetCDF file: ", TRIM(output_nc_filename), " ..."
  CALL initialize_output_netcdf()

  PRINT *, "Step 3: Tracing trajectories in batches..."
  BLOCK
      INTEGER :: num_batches, i_batch, batch_start_idx, batch_end_idx, current_batch_size
      INTEGER :: trajectory_idx_offset
      TYPE(SeedParticle), ALLOCATABLE :: current_batch_seeds(:)

      num_batches = (num_seed_particles_found + TRAJECTORY_BATCH_SIZE - 1) / TRAJECTORY_BATCH_SIZE

      DO i_batch = 1, num_batches
          batch_start_idx = (i_batch - 1) * TRAJECTORY_BATCH_SIZE + 1
          batch_end_idx = MIN(i_batch * TRAJECTORY_BATCH_SIZE, num_seed_particles_found)
          current_batch_size = batch_end_idx - batch_start_idx + 1
          trajectory_idx_offset = batch_start_idx -1 
          
          PRINT *, "Processing batch ", i_batch, " of ", num_batches, &
                   " (seeds ", batch_start_idx, " to ", batch_end_idx, ")"

          ALLOCATE(current_batch_seeds(current_batch_size))
          current_batch_seeds = seed_particles_list(batch_start_idx : batch_end_idx)
          
          CALL trace_trajectories_in_batch(current_batch_seeds, all_trajectories, trajectory_idx_offset, current_batch_size)
          
          DEALLOCATE(current_batch_seeds)
      END DO
  END BLOCK
  
  PRINT *, "All trajectories processed. Writing to NetCDF file..."
  BLOCK
      INTEGER :: traj_write_idx
      DO traj_write_idx = 1, num_seed_particles_found
          CALL write_trajectory_to_netcdf(all_trajectories(traj_write_idx))
          IF (ALLOCATED(all_trajectories(traj_write_idx)%steps)) THEN
              DEALLOCATE(all_trajectories(traj_write_idx)%steps) 
          END IF
      END DO
  END BLOCK
  PRINT *, "All trajectories written."


  err = nf90_close(ncid_out)
  IF (err /= nf90_noerr) THEN
      CALL handle_nc_error(err, "nf90_close output file")
  END IF

  IF (ALLOCATED(seed_particles_list)) THEN
      DEALLOCATE(seed_particles_list)
  END IF
  IF (ALLOCATED(all_trajectories)) THEN
      DEALLOCATE(all_trajectories)
  END IF

  PRINT *, "Particle Tracer Program Finished Successfully."

CONTAINS

  SUBROUTINE parse_time_string(time_str, time_rec)
    CHARACTER(LEN=*), INTENT(IN) :: time_str
    TYPE(TimeRecord), INTENT(OUT) :: time_rec
    CHARACTER(LEN=2) :: hh_str, mm_str, ss_str
    CHARACTER(LEN=3) :: sss_str
    IF (LEN_TRIM(time_str) /= TIME_STR_LEN .OR. time_str(7:7) /= '.') THEN
        PRINT *, "Error: Invalid time string format: ", TRIM(time_str)
        STOP "Invalid time string"
    END IF
    hh_str = time_str(1:2)
    mm_str = time_str(3:4)
    ss_str = time_str(5:6)
    sss_str = time_str(8:10)
    READ(hh_str, '(I2)') time_rec%hours
    READ(mm_str, '(I2)') time_rec%minutes
    READ(ss_str, '(I2)') time_rec%seconds
    READ(sss_str, '(I3)') time_rec%milliseconds
    time_rec%total_milliseconds_since_midnight = ((time_rec%hours*60 + time_rec%minutes)*60 + time_rec%seconds)*1000 + time_rec%milliseconds
  END SUBROUTINE parse_time_string

  SUBROUTINE format_time_string(time_rec, time_str)
    TYPE(TimeRecord), INTENT(IN) :: time_rec
    CHARACTER(LEN=*), INTENT(OUT) :: time_str
    WRITE(time_str, '(I2.2, I2.2, I2.2, ".", I3.3)') time_rec%hours, time_rec%minutes, time_rec%seconds, time_rec%milliseconds
  END SUBROUTINE format_time_string

  SUBROUTINE decrement_time(time_rec_io, interval_ms)
    TYPE(TimeRecord), INTENT(INOUT) :: time_rec_io
    INTEGER, INTENT(IN) :: interval_ms
    INTEGER :: total_ms_current, total_ms_new
    total_ms_current = time_rec_io%total_milliseconds_since_midnight
    total_ms_new = total_ms_current - interval_ms
    IF (total_ms_new < 0) THEN
        PRINT *, "Warning: Decrementing time resulted in negative total_milliseconds. Clamping to 0."
        total_ms_new = 0
    END IF
    time_rec_io%total_milliseconds_since_midnight = total_ms_new
    time_rec_io%milliseconds = MOD(total_ms_new, 1000)
    total_ms_new = total_ms_new / 1000 
    time_rec_io%seconds = MOD(total_ms_new, 60)
    total_ms_new = total_ms_new / 60 
    time_rec_io%minutes = MOD(total_ms_new, 60)
    total_ms_new = total_ms_new / 60 
    time_rec_io%hours = total_ms_new
  END SUBROUTINE decrement_time

  SUBROUTINE generate_nc_filename(time_rec, domain_id, file_prefix, filename_out)
    TYPE(TimeRecord), INTENT(IN) :: time_rec
    INTEGER, INTENT(IN) :: domain_id
    CHARACTER(LEN=*), INTENT(IN) :: file_prefix
    CHARACTER(LEN=*), INTENT(OUT) :: filename_out
    CHARACTER(LEN=TIME_STR_LEN) :: time_s
    CHARACTER(LEN=DOMAIN_ID_STR_LEN) :: domain_s
    CALL format_time_string(time_rec, time_s)
    WRITE(domain_s, '(I6.6)') domain_id
    filename_out = TRIM(path_to_data) // TRIM(file_prefix) // TRIM(date_string) // "-" // TRIM(time_s) // ".pe" // TRIM(domain_s)
  END SUBROUTINE generate_nc_filename

  SUBROUTINE handle_nc_error(status, context)
    USE netcdf
    INTEGER, INTENT(IN) :: status
    CHARACTER(LEN=*), INTENT(IN) :: context
    PRINT *, "NetCDF Error in ", TRIM(context), ": ", nf90_strerror(status)
    STOP "NetCDF Error"
  END SUBROUTINE handle_nc_error

  SUBROUTINE find_seed_particles(target_time_rec)
    USE netcdf
    TYPE(TimeRecord), INTENT(IN) :: target_time_rec
    INTEGER :: d_idx, p_idx_1based, ncid_in, varid_r, varid_z, dimid_sd_num, sd_num_val
    CHARACTER(LEN=FILENAME_MAX_LEN) :: current_filename
    TYPE(SeedParticle), ALLOCATABLE :: temp_seeds(:) 
    INTEGER :: temp_seed_count, initial_alloc_size
    CHARACTER(LEN=SAFE_DIM_NAME_LEN) :: local_dim_name_buffer 
    INTEGER :: err_dim_inq
    LOGICAL :: condition1, condition2
    initial_alloc_size = 1000
    ALLOCATE(temp_seeds(initial_alloc_size))
    temp_seed_count = 0
    DO d_idx = 0, total_domains - 1
        CALL generate_nc_filename(target_time_rec, d_idx, sd_all_prefix, current_filename)
        err = nf90_open(TRIM(current_filename), NF90_NOWRITE, ncid_in)
        IF (err == NF90_ENOTNC) THEN 
            PRINT *, "    Warning: File is not a NetCDF file (ENOTNC), skipping: ", TRIM(current_filename)
            CYCLE
        ELSE IF (err /= nf90_noerr) THEN
            PRINT *, "    Warning: File open failed (err=", err, "), possibly not found or other issue, skipping: ", TRIM(current_filename)
            CYCLE
        END IF
        err = nf90_inq_dimid(ncid_in, "sd_num", dimid_sd_num)
        IF (err /= nf90_noerr) THEN 
            CALL handle_nc_error(err, "nf90_inq_dimid sd_num in " // TRIM(current_filename))
        END IF
        local_dim_name_buffer = " " 
        err_dim_inq = nf90_inquire_dimension(ncid_in, dimid_sd_num, local_dim_name_buffer, sd_num_val) 
        IF (err_dim_inq /= nf90_noerr) THEN
            PRINT *, "    Error inquiring dimension sd_num in ", TRIM(current_filename)
            CALL handle_nc_error(err_dim_inq, "nf90_inquire_dimension for sd_num in " // TRIM(current_filename))
            err = nf90_close(ncid_in)
            IF (err /= nf90_noerr) THEN 
                CALL handle_nc_error(err, "nf90_close after dim_inq_err for " // TRIM(current_filename))
            END IF
            CYCLE
        END IF
        IF (sd_num_val == 0) THEN
            err = nf90_close(ncid_in)
            IF (err /= nf90_noerr) THEN 
                CALL handle_nc_error(err, "nf90_close empty file for " // TRIM(current_filename))
            END IF
            CYCLE
        END IF
        BLOCK
            REAL(KIND=8) :: sd_r_data(sd_num_val), sd_z_data(sd_num_val)
            err = nf90_inq_varid(ncid_in, "sd_r", varid_r)
            IF (err /= nf90_noerr) THEN 
                CALL handle_nc_error(err, "inq sd_r " // TRIM(current_filename))
            END IF
            err = nf90_get_var(ncid_in, varid_r, sd_r_data)
            IF (err /= nf90_noerr) THEN 
                CALL handle_nc_error(err, "get sd_r " // TRIM(current_filename))
            END IF
            err = nf90_inq_varid(ncid_in, "sd_z", varid_z)
            IF (err /= nf90_noerr) THEN 
                CALL handle_nc_error(err, "inq sd_z " // TRIM(current_filename))
            END IF
            err = nf90_get_var(ncid_in, varid_z, sd_z_data)
            IF (err /= nf90_noerr) THEN 
                CALL handle_nc_error(err, "get sd_z " // TRIM(current_filename))
            END IF
            err = nf90_close(ncid_in)
            IF (err /= nf90_noerr) THEN 
                CALL handle_nc_error(err, "close in find_seed for " // TRIM(current_filename))
            END IF
            DO p_idx_1based = 1, sd_num_val
                condition1 = (sd_r_data(p_idx_1based) > filter_radius_gt)
                condition2 = (sd_z_data(p_idx_1based) > filter_height_gt)
                IF (condition1 .AND. condition2) THEN
                    temp_seed_count = temp_seed_count + 1
                    IF (temp_seed_count > SIZE(temp_seeds)) THEN
                        BLOCK
                            TYPE(SeedParticle), ALLOCATABLE :: new_temp_seeds(:)
                            ALLOCATE(new_temp_seeds(SIZE(temp_seeds) + initial_alloc_size))
                            new_temp_seeds(1:SIZE(temp_seeds)) = temp_seeds
                            DEALLOCATE(temp_seeds)
                            CALL MOVE_ALLOC(new_temp_seeds, temp_seeds)
                        END BLOCK
                    END IF
                    CALL format_time_string(target_time_rec, temp_seeds(temp_seed_count)%time_str)
                    temp_seeds(temp_seed_count)%domain_id = d_idx
                    temp_seeds(temp_seed_count)%sd_id = p_idx_1based - 1
                END IF
            END DO
        END BLOCK 
    END DO 
    IF (temp_seed_count > 0) THEN
        ALLOCATE(seed_particles_list(temp_seed_count))
        seed_particles_list = temp_seeds(1:temp_seed_count)
    END IF
    num_seed_particles_found = temp_seed_count
    IF (ALLOCATED(temp_seeds)) DEALLOCATE(temp_seeds)
  END SUBROUTINE find_seed_particles

  SUBROUTINE trace_trajectories_in_batch(seeds_this_batch, all_trajectories_global, &
                                         trajectory_idx_offset_global, actual_batch_size)
    USE netcdf
    TYPE(SeedParticle), INTENT(IN) :: seeds_this_batch(:)
    TYPE(Trajectory), INTENT(INOUT) :: all_trajectories_global(:)
    INTEGER, INTENT(IN) :: trajectory_idx_offset_global 
    INTEGER, INTENT(IN) :: actual_batch_size

    TYPE(ActiveTrajectoryState) :: active_states(actual_batch_size)
    INTEGER :: i_traj_batch, num_active_in_this_batch
    TYPE(TimeRecord) :: current_batch_time, previous_batch_time
    
    TYPE(CachedSDAllFile), ALLOCATABLE :: sd_all_cache(:)
    TYPE(CachedSDCoalFile), ALLOCATABLE :: sd_coal_cache(:)
    INTEGER :: max_cache_entries_sd_all, max_cache_entries_sd_coal 
    INTEGER :: current_cache_count_sd_all, current_cache_count_sd_coal

    CHARACTER(LEN=FILENAME_MAX_LEN) :: temp_fname
    INTEGER :: ncid_temp, varid_temp, dimid_temp, temp_num_particles, temp_num_pairs
    INTEGER :: temp_err_dim_inq, temp_close_status
    CHARACTER(LEN=SAFE_DIM_NAME_LEN) :: local_dim_name_buffer
    LOGICAL :: file_found_in_cache
    INTEGER :: cache_file_idx
    INTEGER :: k_particle, k_coal_pair
    
    REAL(KIND=8), DIMENSION(1) :: temp_r8_val_arr
    INTEGER(KIND=8), DIMENSION(1) :: temp_i8_val_arr
    INTEGER(KIND=2), DIMENSION(1) :: temp_i2_val_arr
    INTEGER, DIMENSION(1) :: temp_int_val_arr 

    DO i_traj_batch = 1, actual_batch_size
        active_states(i_traj_batch)%original_seed_list_idx = trajectory_idx_offset_global + i_traj_batch 
        active_states(i_traj_batch)%global_trajectory_array_idx = trajectory_idx_offset_global + i_traj_batch
        active_states(i_traj_batch)%is_active = .TRUE.
        CALL parse_time_string(seeds_this_batch(i_traj_batch)%time_str, active_states(i_traj_batch)%current_time_rec)
        active_states(i_traj_batch)%current_domain_id = seeds_this_batch(i_traj_batch)%domain_id
        active_states(i_traj_batch)%current_sd_id_0based = seeds_this_batch(i_traj_batch)%sd_id
        active_states(i_traj_batch)%steps_recorded_this_trajectory = 0
        
        BLOCK
            INTEGER :: g_idx, step_num_for_seed
            CHARACTER(LEN=FILENAME_MAX_LEN) :: seed_particle_filename
            INTEGER :: ncid_seed_file
            REAL(KIND=8) :: seed_x, seed_y, seed_z, seed_r
            INTEGER(KIND=8) :: seed_n
            INTEGER(KIND=2) :: seed_if_coal
            INTEGER :: seed_pre_sdid, seed_pre_dmid
            
            g_idx = active_states(i_traj_batch)%global_trajectory_array_idx
            
            all_trajectories_global(g_idx)%trajectory_id_global = g_idx 
            all_trajectories_global(g_idx)%initial_selection_time_str = seeds_this_batch(i_traj_batch)%time_str
            all_trajectories_global(g_idx)%initial_selection_domain_id = seeds_this_batch(i_traj_batch)%domain_id
            all_trajectories_global(g_idx)%initial_selection_sd_id = seeds_this_batch(i_traj_batch)%sd_id
            ALLOCATE(all_trajectories_global(g_idx)%steps(num_time_steps_possible))
            all_trajectories_global(g_idx)%actual_length = 0 
            
            CALL generate_nc_filename(active_states(i_traj_batch)%current_time_rec, &
                                      active_states(i_traj_batch)%current_domain_id, &
                                      sd_all_prefix, seed_particle_filename)
            err = nf90_open(TRIM(seed_particle_filename), NF90_NOWRITE, ncid_seed_file)
            IF (err /= nf90_noerr) THEN
                PRINT *, "Error opening initial seed file for batch: ", TRIM(seed_particle_filename)
                active_states(i_traj_batch)%is_active = .FALSE.
                CYCLE 
            END IF

            err = nf90_inq_varid(ncid_seed_file, "sd_x", varid_temp)
            IF (err /= nf90_noerr) THEN
                CALL handle_nc_error(err, "inq sd_x seed " // TRIM(seed_particle_filename))
            END IF
            err = nf90_get_var(ncid_seed_file, varid_temp, temp_r8_val_arr, start=(/active_states(i_traj_batch)%current_sd_id_0based + 1/), count=(/1/))
            IF (err /= nf90_noerr) THEN
                CALL handle_nc_error(err, "get sd_x seed " // TRIM(seed_particle_filename))
            END IF
            seed_x = temp_r8_val_arr(1)

            err = nf90_inq_varid(ncid_seed_file, "sd_y", varid_temp)
            IF (err /= nf90_noerr) THEN
                CALL handle_nc_error(err, "inq sd_y seed " // TRIM(seed_particle_filename))
            END IF
            err = nf90_get_var(ncid_seed_file, varid_temp, temp_r8_val_arr, start=(/active_states(i_traj_batch)%current_sd_id_0based + 1/), count=(/1/))
            IF (err /= nf90_noerr) THEN
                CALL handle_nc_error(err, "get sd_y seed " // TRIM(seed_particle_filename))
            END IF
            seed_y = temp_r8_val_arr(1)

            err = nf90_inq_varid(ncid_seed_file, "sd_z", varid_temp)
            IF (err /= nf90_noerr) THEN
                CALL handle_nc_error(err, "inq sd_z seed " // TRIM(seed_particle_filename))
            END IF
            err = nf90_get_var(ncid_seed_file, varid_temp, temp_r8_val_arr, start=(/active_states(i_traj_batch)%current_sd_id_0based + 1/), count=(/1/))
            IF (err /= nf90_noerr) THEN
                CALL handle_nc_error(err, "get sd_z seed " // TRIM(seed_particle_filename))
            END IF
            seed_z = temp_r8_val_arr(1)

            err = nf90_inq_varid(ncid_seed_file, "sd_r", varid_temp)
            IF (err /= nf90_noerr) THEN
                CALL handle_nc_error(err, "inq sd_r seed " // TRIM(seed_particle_filename))
            END IF
            err = nf90_get_var(ncid_seed_file, varid_temp, temp_r8_val_arr, start=(/active_states(i_traj_batch)%current_sd_id_0based + 1/), count=(/1/))
            IF (err /= nf90_noerr) THEN
                CALL handle_nc_error(err, "get sd_r seed " // TRIM(seed_particle_filename))
            END IF
            seed_r = temp_r8_val_arr(1)

            err = nf90_inq_varid(ncid_seed_file, "sd_n", varid_temp)
            IF (err /= nf90_noerr) THEN
                CALL handle_nc_error(err, "inq sd_n seed " // TRIM(seed_particle_filename))
            END IF
            err = nf90_get_var(ncid_seed_file, varid_temp, temp_i8_val_arr, start=(/active_states(i_traj_batch)%current_sd_id_0based + 1/), count=(/1/))
            IF (err /= nf90_noerr) THEN
                CALL handle_nc_error(err, "get sd_n seed " // TRIM(seed_particle_filename))
            END IF
            seed_n = temp_i8_val_arr(1)

            err = nf90_inq_varid(ncid_seed_file, "if_coal", varid_temp)
            IF (err /= nf90_noerr) THEN
                CALL handle_nc_error(err, "inq if_coal seed " // TRIM(seed_particle_filename))
            END IF
            err = nf90_get_var(ncid_seed_file, varid_temp, temp_i2_val_arr, start=(/active_states(i_traj_batch)%current_sd_id_0based + 1/), count=(/1/))
            IF (err /= nf90_noerr) THEN
                CALL handle_nc_error(err, "get if_coal seed " // TRIM(seed_particle_filename))
            END IF
            seed_if_coal = temp_i2_val_arr(1)

            err = nf90_inq_varid(ncid_seed_file, "pre_sdid", varid_temp)
            IF (err /= nf90_noerr) THEN 
                CALL handle_nc_error(err, "inq pre_sdid seed " // TRIM(seed_particle_filename))
            END IF
            err = nf90_get_var(ncid_seed_file, varid_temp, temp_int_val_arr, & 
                               start=(/active_states(i_traj_batch)%current_sd_id_0based + 1/), count=(/1/))
            IF (err /= nf90_noerr) THEN 
                CALL handle_nc_error(err, "get pre_sdid seed " // TRIM(seed_particle_filename))
            END IF
            seed_pre_sdid = temp_int_val_arr(1)

            err = nf90_inq_varid(ncid_seed_file, "pre_dmid", varid_temp)
            IF (err /= nf90_noerr) THEN 
                CALL handle_nc_error(err, "inq pre_dmid seed " // TRIM(seed_particle_filename))
            END IF
            err = nf90_get_var(ncid_seed_file, varid_temp, temp_int_val_arr, & 
                               start=(/active_states(i_traj_batch)%current_sd_id_0based + 1/), count=(/1/))
            IF (err /= nf90_noerr) THEN 
                CALL handle_nc_error(err, "get pre_dmid seed " // TRIM(seed_particle_filename))
            END IF
            seed_pre_dmid = temp_int_val_arr(1)
            
            err = nf90_close(ncid_seed_file)
            IF (err /= nf90_noerr) THEN 
                CALL handle_nc_error(err, "close seed file " // TRIM(seed_particle_filename))
            END IF

            all_trajectories_global(g_idx)%actual_length = 1
            step_num_for_seed = 1
            all_trajectories_global(g_idx)%steps(step_num_for_seed)%is_valid = .TRUE.
            CALL format_time_string(active_states(i_traj_batch)%current_time_rec, all_trajectories_global(g_idx)%steps(step_num_for_seed)%timestamp_str)
            all_trajectories_global(g_idx)%steps(step_num_for_seed)%domain_id = active_states(i_traj_batch)%current_domain_id
            all_trajectories_global(g_idx)%steps(step_num_for_seed)%sd_id_in_file = active_states(i_traj_batch)%current_sd_id_0based
            all_trajectories_global(g_idx)%steps(step_num_for_seed)%x_coord = seed_x
            all_trajectories_global(g_idx)%steps(step_num_for_seed)%y_coord = seed_y 
            all_trajectories_global(g_idx)%steps(step_num_for_seed)%z_coord = seed_z 
            all_trajectories_global(g_idx)%steps(step_num_for_seed)%radius = seed_r   
            all_trajectories_global(g_idx)%steps(step_num_for_seed)%multiplicity = seed_n 
            all_trajectories_global(g_idx)%steps(step_num_for_seed)%if_coal_flag = seed_if_coal 
            all_trajectories_global(g_idx)%steps(step_num_for_seed)%parent_dmid_for_trace = seed_pre_dmid
            all_trajectories_global(g_idx)%steps(step_num_for_seed)%parent_sdid_for_trace = seed_pre_sdid
            
            active_states(i_traj_batch)%steps_recorded_this_trajectory = 1
            active_states(i_traj_batch)%next_parent_dmid_to_find = seed_pre_dmid
            active_states(i_traj_batch)%next_parent_sdid_1based_to_find = seed_pre_sdid

            IF (active_states(i_traj_batch)%next_parent_sdid_1based_to_find <= 0) THEN
                active_states(i_traj_batch)%is_active = .FALSE. 
            END IF
        END BLOCK
    END DO
    
    num_active_in_this_batch = COUNT(active_states(:)%is_active)
    IF (num_active_in_this_batch == 0) THEN
        PRINT *, "No active trajectories to trace in this batch after initial seed processing."
        RETURN
    END IF

    current_batch_time = active_states(1)%current_time_rec 

    max_cache_entries_sd_all = total_domains 
    max_cache_entries_sd_coal = total_domains
    ALLOCATE(sd_all_cache(max_cache_entries_sd_all))
    ALLOCATE(sd_coal_cache(max_cache_entries_sd_coal))

BATCH_TIME_STEP_LOOP : DO WHILE (num_active_in_this_batch > 0)
        
        CALL decrement_time(current_batch_time, time_interval_milliseconds)
        previous_batch_time = current_batch_time 

        IF (previous_batch_time%total_milliseconds_since_midnight < time_start_rec%total_milliseconds_since_midnight) THEN
            EXIT BATCH_TIME_STEP_LOOP
        END IF

        current_cache_count_sd_all = 0
        current_cache_count_sd_coal = 0
        
        DO i_traj_batch = 1, SIZE(sd_all_cache)
            IF (ALLOCATED(sd_all_cache(i_traj_batch)%particles)) THEN
                 DEALLOCATE(sd_all_cache(i_traj_batch)%particles)
            END IF
            sd_all_cache(i_traj_batch)%is_loaded = .FALSE.
        END DO
        DO i_traj_batch = 1, SIZE(sd_coal_cache)
            IF (ALLOCATED(sd_coal_cache(i_traj_batch)%coal_pairs)) THEN
                DEALLOCATE(sd_coal_cache(i_traj_batch)%coal_pairs)
            END IF
            sd_coal_cache(i_traj_batch)%is_loaded = .FALSE.
        END DO

        DO i_traj_batch = 1, actual_batch_size
            IF (.NOT. active_states(i_traj_batch)%is_active) CYCLE
            IF (active_states(i_traj_batch)%next_parent_sdid_1based_to_find <= 0) THEN
                 active_states(i_traj_batch)%is_active = .FALSE.
                 num_active_in_this_batch = num_active_in_this_batch -1
                 CYCLE
            END IF

            CALL generate_nc_filename(previous_batch_time, active_states(i_traj_batch)%next_parent_dmid_to_find, sd_all_prefix, temp_fname)
            CALL load_sd_all_file_to_cache(TRIM(temp_fname), previous_batch_time, active_states(i_traj_batch)%next_parent_dmid_to_find, &
                                           sd_all_cache, current_cache_count_sd_all, max_cache_entries_sd_all)
        END DO

        DO i_traj_batch = 1, actual_batch_size
            IF (.NOT. active_states(i_traj_batch)%is_active) CYCLE
            
            BLOCK
                INTEGER :: g_idx, step_num
                TYPE(CachedParticleInfo) :: parent_particle_data
                LOGICAL :: parent_found_in_cache
                INTEGER :: parent_domain, parent_sdid_1based
                
                g_idx = active_states(i_traj_batch)%global_trajectory_array_idx
                parent_domain = active_states(i_traj_batch)%next_parent_dmid_to_find
                parent_sdid_1based = active_states(i_traj_batch)%next_parent_sdid_1based_to_find

                parent_found_in_cache = .FALSE.
                CALL find_particle_in_sd_all_cache(previous_batch_time, parent_domain, parent_sdid_1based, &
                                                   sd_all_cache, current_cache_count_sd_all, &
                                                   parent_particle_data, parent_found_in_cache)
                
                IF (parent_found_in_cache) THEN
                    all_trajectories_global(g_idx)%actual_length = all_trajectories_global(g_idx)%actual_length + 1
                    step_num = all_trajectories_global(g_idx)%actual_length
                    IF (step_num > num_time_steps_possible) THEN
                        PRINT *, "Error: Exceeded max trajectory steps for traj ", g_idx
                        active_states(i_traj_batch)%is_active = .FALSE.
                        num_active_in_this_batch = num_active_in_this_batch - 1
                        CYCLE
                    END IF

                    CALL format_time_string(previous_batch_time, all_trajectories_global(g_idx)%steps(step_num)%timestamp_str)
                    all_trajectories_global(g_idx)%steps(step_num)%domain_id = parent_domain
                    all_trajectories_global(g_idx)%steps(step_num)%sd_id_in_file = parent_sdid_1based - 1 
                    all_trajectories_global(g_idx)%steps(step_num)%x_coord = parent_particle_data%x_coord
                    all_trajectories_global(g_idx)%steps(step_num)%y_coord = parent_particle_data%y_coord
                    all_trajectories_global(g_idx)%steps(step_num)%z_coord = parent_particle_data%z_coord
                    all_trajectories_global(g_idx)%steps(step_num)%radius = parent_particle_data%radius
                    all_trajectories_global(g_idx)%steps(step_num)%multiplicity = parent_particle_data%multiplicity
                    all_trajectories_global(g_idx)%steps(step_num)%if_coal_flag = parent_particle_data%if_coal_flag
                    all_trajectories_global(g_idx)%steps(step_num)%parent_dmid_for_trace = parent_particle_data%pre_dmid
                    all_trajectories_global(g_idx)%steps(step_num)%parent_sdid_for_trace = parent_particle_data%pre_sdid_1based
                    all_trajectories_global(g_idx)%steps(step_num)%is_valid = .TRUE.

                    active_states(i_traj_batch)%current_time_rec = previous_batch_time
                    active_states(i_traj_batch)%current_domain_id = parent_domain
                    active_states(i_traj_batch)%current_sd_id_0based = parent_sdid_1based - 1
                    active_states(i_traj_batch)%next_parent_dmid_to_find = parent_particle_data%pre_dmid
                    active_states(i_traj_batch)%next_parent_sdid_1based_to_find = parent_particle_data%pre_sdid_1based
                    active_states(i_traj_batch)%steps_recorded_this_trajectory = step_num

                    IF (parent_particle_data%if_coal_flag == 1) THEN
                        ! TODO: Implement SD_coal caching and lookup here
                    END IF

                    IF (parent_particle_data%pre_sdid_1based <= 0) THEN
                        active_states(i_traj_batch)%is_active = .FALSE.
                        num_active_in_this_batch = num_active_in_this_batch - 1
                    END IF
                ELSE
                    PRINT *, "    Trace Warning: Parent not found in cache for trajectory ", g_idx, &
                             " at time ", previous_batch_time%total_milliseconds_since_midnight, &
                             " domain ", parent_domain, " sdid ", parent_sdid_1based
                    active_states(i_traj_batch)%is_active = .FALSE.
                    num_active_in_this_batch = num_active_in_this_batch - 1
                END IF
            END BLOCK
        END DO 

        IF (num_active_in_this_batch == 0) THEN
            EXIT BATCH_TIME_STEP_LOOP
        END IF

    END DO BATCH_TIME_STEP_LOOP

    IF (ALLOCATED(sd_all_cache)) THEN
        DO i_traj_batch = 1, SIZE(sd_all_cache)
            IF (ALLOCATED(sd_all_cache(i_traj_batch)%particles)) DEALLOCATE(sd_all_cache(i_traj_batch)%particles)
        END DO
        DEALLOCATE(sd_all_cache)
    END IF
    IF (ALLOCATED(sd_coal_cache)) THEN
         DO i_traj_batch = 1, SIZE(sd_coal_cache)
            IF (ALLOCATED(sd_coal_cache(i_traj_batch)%coal_pairs)) DEALLOCATE(sd_coal_cache(i_traj_batch)%coal_pairs)
        END DO
        DEALLOCATE(sd_coal_cache)
    END IF

  END SUBROUTINE trace_trajectories_in_batch

  SUBROUTINE load_sd_all_file_to_cache(filename_to_load, file_time, file_domain, &
                                       cache_array, current_cache_count, max_cache_size)
    USE netcdf
    CHARACTER(LEN=*), INTENT(IN) :: filename_to_load
    TYPE(TimeRecord), INTENT(IN) :: file_time
    INTEGER, INTENT(IN) :: file_domain
    TYPE(CachedSDAllFile), INTENT(INOUT), ALLOCATABLE :: cache_array(:)
    INTEGER, INTENT(INOUT) :: current_cache_count
    INTEGER, INTENT(IN) :: max_cache_size
    
    INTEGER :: i, ncid, num_p, varid_tmp, dimid_tmp
    CHARACTER(LEN=SAFE_DIM_NAME_LEN) :: local_dim_name_buffer
    INTEGER :: file_cache_slot

    DO i = 1, current_cache_count
        IF (cache_array(i)%is_loaded .AND. TRIM(cache_array(i)%filename) == TRIM(filename_to_load)) THEN
            RETURN 
        END IF
    END DO

    IF (current_cache_count >= max_cache_size) THEN
        PRINT *, "Error: SD_all cache full. Increase max_cache_entries_sd_all."
        STOP "SD_all cache overflow"
    END IF

    current_cache_count = current_cache_count + 1
    file_cache_slot = current_cache_count
    
    cache_array(file_cache_slot)%filename = TRIM(filename_to_load)
    cache_array(file_cache_slot)%file_time_rec = file_time
    cache_array(file_cache_slot)%domain_id = file_domain
    cache_array(file_cache_slot)%is_loaded = .FALSE. 

    err = nf90_open(TRIM(filename_to_load), NF90_NOWRITE, ncid)
    IF (err /= nf90_noerr) THEN
        PRINT *, "Cache Load Error: Could not open ", TRIM(filename_to_load), nf90_strerror(err)
        current_cache_count = current_cache_count - 1 
        RETURN
    END IF

    err = nf90_inq_dimid(ncid, "sd_num", dimid_tmp)
    IF (err /= nf90_noerr) THEN 
        CALL handle_nc_error(err, "cache_load inq_dimid " // TRIM(filename_to_load))
    END IF
    local_dim_name_buffer = " "
    err = nf90_inquire_dimension(ncid, dimid_tmp, local_dim_name_buffer, num_p)
    IF (err /= nf90_noerr) THEN 
        CALL handle_nc_error(err, "cache_load inquire_dim " // TRIM(filename_to_load))
    END IF
    
    cache_array(file_cache_slot)%num_particles = num_p
    IF (num_p > 0) THEN
        ALLOCATE(cache_array(file_cache_slot)%particles(num_p))
        BLOCK 
            REAL(KIND=8) :: temp_x(num_p), temp_y(num_p), temp_z(num_p), temp_r(num_p)
            INTEGER(KIND=8) :: temp_n(num_p)
            INTEGER(KIND=2) :: temp_if_coal(num_p)
            INTEGER :: temp_pre_sdid(num_p), temp_pre_dmid(num_p)

            err = nf90_inq_varid(ncid, "sd_x", varid_tmp)
            IF(err/=nf90_noerr) THEN 
                CALL handle_nc_error(err,"cache_sd_x_vid")
            END IF
            err = nf90_get_var(ncid, varid_tmp, temp_x)
            IF(err/=nf90_noerr) THEN 
                CALL handle_nc_error(err,"cache_sd_x_get")
            END IF
            err = nf90_inq_varid(ncid, "sd_y", varid_tmp)
            IF(err/=nf90_noerr) THEN 
                CALL handle_nc_error(err,"cache_sd_y_vid")
            END IF
            err = nf90_get_var(ncid, varid_tmp, temp_y)
            IF(err/=nf90_noerr) THEN 
                CALL handle_nc_error(err,"cache_sd_y_get")
            END IF
            err = nf90_inq_varid(ncid, "sd_z", varid_tmp)
            IF(err/=nf90_noerr) THEN 
                CALL handle_nc_error(err,"cache_sd_z_vid")
            END IF
            err = nf90_get_var(ncid, varid_tmp, temp_z)
            IF(err/=nf90_noerr) THEN 
                CALL handle_nc_error(err,"cache_sd_z_get")
            END IF
            err = nf90_inq_varid(ncid, "sd_r", varid_tmp)
            IF(err/=nf90_noerr) THEN 
                CALL handle_nc_error(err,"cache_sd_r_vid")
            END IF
            err = nf90_get_var(ncid, varid_tmp, temp_r)
            IF(err/=nf90_noerr) THEN 
                CALL handle_nc_error(err,"cache_sd_r_get")
            END IF
            err = nf90_inq_varid(ncid, "sd_n", varid_tmp)
            IF(err/=nf90_noerr) THEN 
                CALL handle_nc_error(err,"cache_sd_n_vid")
            END IF
            err = nf90_get_var(ncid, varid_tmp, temp_n)
            IF(err/=nf90_noerr) THEN 
                CALL handle_nc_error(err,"cache_sd_n_get")
            END IF
            err = nf90_inq_varid(ncid, "if_coal", varid_tmp)
            IF(err/=nf90_noerr) THEN 
                CALL handle_nc_error(err,"cache_if_coal_vid")
            END IF
            err = nf90_get_var(ncid, varid_tmp, temp_if_coal)
            IF(err/=nf90_noerr) THEN 
                CALL handle_nc_error(err,"cache_if_coal_get")
            END IF
            err = nf90_inq_varid(ncid, "pre_sdid", varid_tmp)
            IF(err/=nf90_noerr) THEN 
                CALL handle_nc_error(err,"cache_pre_sdid_vid")
            END IF
            err = nf90_get_var(ncid, varid_tmp, temp_pre_sdid)
            IF(err/=nf90_noerr) THEN 
                CALL handle_nc_error(err,"cache_pre_sdid_get")
            END IF
            err = nf90_inq_varid(ncid, "pre_dmid", varid_tmp)
            IF(err/=nf90_noerr) THEN 
                CALL handle_nc_error(err,"cache_pre_dmid_vid")
            END IF
            err = nf90_get_var(ncid, varid_tmp, temp_pre_dmid)
            IF(err/=nf90_noerr) THEN 
                CALL handle_nc_error(err,"cache_pre_dmid_get")
            END IF

            DO i = 1, num_p
                cache_array(file_cache_slot)%particles(i)%x_coord = temp_x(i)
                cache_array(file_cache_slot)%particles(i)%y_coord = temp_y(i)
                cache_array(file_cache_slot)%particles(i)%z_coord = temp_z(i)
                cache_array(file_cache_slot)%particles(i)%radius = temp_r(i)
                cache_array(file_cache_slot)%particles(i)%multiplicity = temp_n(i)
                cache_array(file_cache_slot)%particles(i)%if_coal_flag = temp_if_coal(i)
                cache_array(file_cache_slot)%particles(i)%pre_sdid_1based = temp_pre_sdid(i)
                cache_array(file_cache_slot)%particles(i)%pre_dmid = temp_pre_dmid(i)
            END DO
        END BLOCK
    END IF
    cache_array(file_cache_slot)%is_loaded = .TRUE.
    err = nf90_close(ncid)
    IF (err /= nf90_noerr) THEN 
        CALL handle_nc_error(err, "cache_load close " // TRIM(filename_to_load))
    END IF

  END SUBROUTINE load_sd_all_file_to_cache

  SUBROUTINE find_particle_in_sd_all_cache(target_time, target_domain, target_sdid_1based, &
                                           cache_array, current_cache_size, &
                                           particle_data_out, found)
    TYPE(TimeRecord), INTENT(IN) :: target_time
    INTEGER, INTENT(IN) :: target_domain, target_sdid_1based
    TYPE(CachedSDAllFile), INTENT(IN) :: cache_array(:)
    INTEGER, INTENT(IN) :: current_cache_size
    TYPE(CachedParticleInfo), INTENT(OUT) :: particle_data_out
    LOGICAL, INTENT(OUT) :: found
    
    INTEGER :: i_cache, target_sdid_0based
    CHARACTER(LEN=FILENAME_MAX_LEN) :: target_filename
    
    found = .FALSE.
    IF (target_sdid_1based <= 0) THEN
        RETURN
    END IF

    CALL generate_nc_filename(target_time, target_domain, sd_all_prefix, target_filename)
    target_sdid_0based = target_sdid_1based - 1

    DO i_cache = 1, current_cache_size
        IF (cache_array(i_cache)%is_loaded .AND. TRIM(cache_array(i_cache)%filename) == TRIM(target_filename)) THEN
            IF (target_sdid_0based >= 0 .AND. target_sdid_0based < cache_array(i_cache)%num_particles) THEN
                particle_data_out = cache_array(i_cache)%particles(target_sdid_0based + 1)
                found = .TRUE.
                RETURN
            ELSE
                PRINT *, "Warning: Particle sdid ", target_sdid_1based, " out of bounds for cached file ", TRIM(target_filename)
                RETURN
            END IF
        END IF
    END DO
  END SUBROUTINE find_particle_in_sd_all_cache
  
  SUBROUTINE initialize_output_netcdf()
    USE netcdf
    INTEGER :: cmode 
    cmode = NF90_CLOBBER + NF90_NETCDF4 
    err = nf90_create(TRIM(output_nc_filename), cmode, ncid_out)
    IF (err /= nf90_noerr) THEN 
        CALL handle_nc_error(err, "nf90_create output file (NetCDF-4)")
    END IF
    err = nf90_def_dim(ncid_out, "trajectory_id", NF90_UNLIMITED, dimid_traj)
    IF (err /= nf90_noerr) THEN 
        CALL handle_nc_error(err, "nf90_def_dim trajectory_id")
    END IF
    err = nf90_def_dim(ncid_out, "time_step_index", num_time_steps_possible, dimid_time_step)
    IF (err /= nf90_noerr) THEN 
        CALL handle_nc_error(err, "nf90_def_dim time_step_index")
    END IF
    err = nf90_def_dim(ncid_out, "char_dim_time", TIME_STR_LEN, dimid_char_time)
    IF (err /= nf90_noerr) THEN 
        CALL handle_nc_error(err, "nf90_def_dim char_dim_time")
    END IF
    err = nf90_put_att(ncid_out, NF90_GLOBAL, "description", "Traced particle data")
    IF (err /= nf90_noerr) THEN 
        CALL handle_nc_error(err, "put_att description")
    END IF
    ! err = nf90_def_var(ncid_out, "initial_selection_time", NF90_CHAR, (/dimid_char_time, dimid_traj /), varid_init_sel_time)
    ! IF (err /= nf90_noerr) THEN 
    !     CALL handle_nc_error(err, "def_var initial_selection_time")
    ! END IF
    ! err = nf90_def_var(ncid_out, "initial_selection_domain_id", NF90_INT, (/dimid_traj/), varid_init_sel_dom)
    ! IF (err /= nf90_noerr) THEN 
    !     CALL handle_nc_error(err, "def_var initial_selection_domain_id")
    ! END IF
    ! err = nf90_def_var(ncid_out, "initial_selection_sd_id", NF90_INT, (/dimid_traj/), varid_init_sel_sdid)
    ! IF (err /= nf90_noerr) THEN 
    !     CALL handle_nc_error(err, "def_var initial_selection_sd_id")
    ! END IF
    ! err = nf90_def_var(ncid_out, "trajectory_actual_length", NF90_INT, (/dimid_traj/), varid_traj_len)
    ! IF (err /= nf90_noerr) THEN 
    !     CALL handle_nc_error(err, "def_var trajectory_actual_length")
    ! END IF
    ! err = nf90_def_var(ncid_out, "is_valid_step", NF90_BYTE, (/dimid_time_step, dimid_traj/), varid_is_valid_step)
    ! IF (err /= nf90_noerr) THEN 
    !     CALL handle_nc_error(err, "def_var is_valid_step")
    ! END IF
    ! err = nf90_put_att(ncid_out, varid_is_valid_step, "_FillValue", FILL_BYTE)
    ! IF (err /= nf90_noerr) THEN 
    !     CALL handle_nc_error(err, "put_att is_valid_step _FillValue")
    ! END IF
    ! err = nf90_def_var(ncid_out, "timestamp_at_step", NF90_CHAR, (/dimid_char_time, dimid_time_step, dimid_traj/), varid_ts_at_step) 
    ! IF (err /= nf90_noerr) THEN 
    !     CALL handle_nc_error(err, "def_var timestamp_at_step")
    ! END IF
    err = nf90_def_var(ncid_out, "domain_id_at_step", NF90_INT, (/dimid_time_step, dimid_traj/), varid_dom_at_step) 
    IF (err /= nf90_noerr) THEN 
        CALL handle_nc_error(err, "def_var domain_id_at_step")
    END IF
    err = nf90_put_att(ncid_out, varid_dom_at_step, "_FillValue", FILL_INT)
    IF (err /= nf90_noerr) THEN 
        CALL handle_nc_error(err, "put_att domain_id_at_step _FillValue")
    END IF
    err = nf90_def_var(ncid_out, "sd_id_in_file_at_step", NF90_INT, (/dimid_time_step, dimid_traj/), varid_sdid_at_step) 
    IF (err /= nf90_noerr) THEN 
        CALL handle_nc_error(err, "def_var sd_id_in_file_at_step")
    END IF
    err = nf90_put_att(ncid_out, varid_sdid_at_step, "_FillValue", FILL_INT)
    IF (err /= nf90_noerr) THEN 
        CALL handle_nc_error(err, "put_att sd_id_in_file_at_step _FillValue")
    END IF
    err = nf90_def_var(ncid_out, "x_coord", NF90_DOUBLE, (/dimid_time_step, dimid_traj/), varid_x) 
    IF (err /= nf90_noerr) THEN 
        CALL handle_nc_error(err, "def_var x_coord")
    END IF
    err = nf90_put_att(ncid_out, varid_x, "_FillValue", FILL_DOUBLE)
    IF (err /= nf90_noerr) THEN 
        CALL handle_nc_error(err,"FA_x_fill")
    END IF
    err = nf90_def_var(ncid_out, "y_coord", NF90_DOUBLE, (/dimid_time_step, dimid_traj/), varid_y) 
    IF (err /= nf90_noerr) THEN 
        CALL handle_nc_error(err, "def_var y_coord")
    END IF
    err = nf90_put_att(ncid_out, varid_y, "_FillValue", FILL_DOUBLE)
    IF (err /= nf90_noerr) THEN 
        CALL handle_nc_error(err,"FA_y_fill")
    END IF
    err = nf90_def_var(ncid_out, "z_coord", NF90_DOUBLE, (/dimid_time_step, dimid_traj/), varid_z) 
    IF (err /= nf90_noerr) THEN 
        CALL handle_nc_error(err, "def_var z_coord")
    END IF
    err = nf90_put_att(ncid_out, varid_z, "_FillValue", FILL_DOUBLE)
    IF (err /= nf90_noerr) THEN 
        CALL handle_nc_error(err,"FA_z_fill")
    END IF
    err = nf90_def_var(ncid_out, "radius", NF90_DOUBLE, (/dimid_time_step, dimid_traj/), varid_r) 
    IF (err /= nf90_noerr) THEN 
        CALL handle_nc_error(err, "def_var radius")
    END IF
    err = nf90_put_att(ncid_out, varid_r, "_FillValue", FILL_DOUBLE)
    IF (err /= nf90_noerr) THEN 
        CALL handle_nc_error(err,"FA_r_fill")
    END IF
    err = nf90_def_var(ncid_out, "multiplicity", NF90_INT64, (/dimid_time_step, dimid_traj/), varid_n) 
    IF (err /= nf90_noerr) THEN 
        CALL handle_nc_error(err, "def_var multiplicity")
    END IF
    err = nf90_put_att(ncid_out, varid_n, "_FillValue", FILL_INT64)
    IF (err /= nf90_noerr) THEN 
        CALL handle_nc_error(err, "put_att n _FillValue")
    END IF
    err = nf90_def_var(ncid_out, "if_coal_flag", NF90_SHORT, (/dimid_time_step, dimid_traj/), varid_if_coal) 
    IF (err /= nf90_noerr) THEN 
        CALL handle_nc_error(err, "def_var if_coal_flag")
    END IF
    err = nf90_put_att(ncid_out, varid_if_coal, "_FillValue", FILL_SHORT)
    IF (err /= nf90_noerr) THEN 
        CALL handle_nc_error(err, "put_att if_coal _FillValue")
    END IF
    err = nf90_def_var(ncid_out, "parent_sdid_for_trace", NF90_INT, (/dimid_time_step, dimid_traj/), varid_parent_sdid) 
    IF (err /= nf90_noerr) THEN 
        CALL handle_nc_error(err, "def_var parent_sdid_for_trace")
    END IF
    err = nf90_put_att(ncid_out, varid_parent_sdid, "_FillValue", FILL_INT)
    IF (err /= nf90_noerr) THEN 
        CALL handle_nc_error(err, "put_att parent_sdid _FillValue")
    END IF
    err = nf90_def_var(ncid_out, "parent_dmid_for_trace", NF90_INT, (/dimid_time_step, dimid_traj/), varid_parent_dmid) 
    IF (err /= nf90_noerr) THEN 
        CALL handle_nc_error(err, "def_var parent_dmid_for_trace")
    END IF
    err = nf90_put_att(ncid_out, varid_parent_dmid, "_FillValue", FILL_INT)
    IF (err /= nf90_noerr) THEN 
        CALL handle_nc_error(err, "put_att parent_dmid _FillValue")
    END IF
    err = nf90_def_var(ncid_out, "coalescence_num_col", NF90_INT, (/dimid_time_step, dimid_traj/), varid_coal_num_col) 
    IF (err /= nf90_noerr) THEN 
        CALL handle_nc_error(err, "def_var coalescence_num_col")
    END IF
    err = nf90_put_att(ncid_out, varid_coal_num_col, "_FillValue", FILL_INT)
    IF (err /= nf90_noerr) THEN 
        CALL handle_nc_error(err, "put_att coal_num_col _FillValue")
    END IF
    err = nf90_def_var(ncid_out, "other_parent_prev_dmid", NF90_INT, (/dimid_time_step, dimid_traj/), varid_other_dmid) 
    IF (err /= nf90_noerr) THEN 
        CALL handle_nc_error(err, "def_var other_parent_prev_dmid")
    END IF
    err = nf90_put_att(ncid_out, varid_other_dmid, "_FillValue", FILL_INT)
    IF (err /= nf90_noerr) THEN 
        CALL handle_nc_error(err, "put_att other_dmid _FillValue")
    END IF
    err = nf90_def_var(ncid_out, "other_parent_prev_sdid", NF90_INT, (/dimid_time_step, dimid_traj/), varid_other_sdid) 
    IF (err /= nf90_noerr) THEN 
        CALL handle_nc_error(err, "def_var other_parent_prev_sdid")
    END IF
    err = nf90_put_att(ncid_out, varid_other_sdid, "_FillValue", FILL_INT)
    IF (err /= nf90_noerr) THEN 
        CALL handle_nc_error(err, "put_att other_sdid _FillValue")
    END IF
    err = nf90_def_var(ncid_out, "other_parent_radius", NF90_DOUBLE, (/dimid_time_step, dimid_traj/), varid_other_r) 
    IF (err /= nf90_noerr) THEN 
        CALL handle_nc_error(err, "def_var other_parent_radius")
    END IF
    err = nf90_put_att(ncid_out, varid_other_r, "_FillValue", FILL_DOUBLE)
    IF (err /= nf90_noerr) THEN 
        CALL handle_nc_error(err,"FA_other_r_fill")
    END IF
    err = nf90_def_var(ncid_out, "other_parent_multiplicity", NF90_INT64, (/dimid_time_step, dimid_traj/), varid_other_n) 
    IF (err /= nf90_noerr) THEN 
        CALL handle_nc_error(err, "def_var other_parent_multiplicity")
    END IF
    err = nf90_put_att(ncid_out, varid_other_n, "_FillValue", FILL_INT64)
    IF (err /= nf90_noerr) THEN 
        CALL handle_nc_error(err, "put_att other_n _FillValue")
    END IF
    err = nf90_enddef(ncid_out)
    IF (err /= nf90_noerr) THEN 
        CALL handle_nc_error(err, "nf90_enddef output file")
    END IF
  END SUBROUTINE initialize_output_netcdf

  SUBROUTINE write_trajectory_to_netcdf(traj_data)
    USE netcdf
    TYPE(Trajectory), INTENT(IN) :: traj_data
    INTEGER :: i
    INTEGER(KIND=1) :: is_valid_byte_arr(num_time_steps_possible)
    CHARACTER(LEN=TIME_STR_LEN) :: ts_char_arr(num_time_steps_possible)
    INTEGER :: dom_id_arr(num_time_steps_possible)
    INTEGER :: sd_id_0based_arr(num_time_steps_possible)
    REAL(KIND=8) :: x_arr(num_time_steps_possible), y_arr(num_time_steps_possible), z_arr(num_time_steps_possible)
    REAL(KIND=8) :: r_arr(num_time_steps_possible)
    INTEGER(KIND=8) :: n_arr(num_time_steps_possible)
    INTEGER(KIND=2) :: if_coal_arr(num_time_steps_possible)
    INTEGER :: parent_sdid_1based_arr(num_time_steps_possible)
    INTEGER :: parent_dmid_arr(num_time_steps_possible)
    INTEGER :: coal_num_arr(num_time_steps_possible)
    INTEGER :: other_dmid_arr(num_time_steps_possible)
    INTEGER :: other_sdid_1based_arr(num_time_steps_possible)
    REAL(KIND=8) :: other_r_arr(num_time_steps_possible)
    INTEGER(KIND=8) :: other_n_arr(num_time_steps_possible)
    INTEGER :: current_record_idx
    current_record_idx = traj_data%trajectory_id_global

    ! err = nf90_put_var(ncid_out,varid_init_sel_time,TRIM(traj_data%initial_selection_time_str),start=(/1,current_record_idx/),count=(/TIME_STR_LEN,1/))
    ! IF(err/=nf90_noerr)THEN 
    !     CALL handle_nc_error(err,"put_var initial_selection_time")
    ! END IF
    ! err = nf90_put_var(ncid_out,varid_init_sel_dom,(/traj_data%initial_selection_domain_id/),start=(/current_record_idx/),count=(/1/))
    ! IF(err/=nf90_noerr)THEN 
    !     CALL handle_nc_error(err,"put_var initial_selection_domain_id")
    ! END IF
    ! err = nf90_put_var(ncid_out,varid_init_sel_sdid,(/traj_data%initial_selection_sd_id/),start=(/current_record_idx/),count=(/1/))
    ! IF(err/=nf90_noerr)THEN 
    !     CALL handle_nc_error(err,"put_var initial_selection_sd_id")
    ! END IF
    ! err = nf90_put_var(ncid_out,varid_traj_len,(/traj_data%actual_length/),start=(/current_record_idx/),count=(/1/))
    ! IF(err/=nf90_noerr)THEN 
    !     CALL handle_nc_error(err,"put_var trajectory_actual_length")
    ! END IF

    DO i=1,num_time_steps_possible
        IF(i<=traj_data%actual_length.AND.traj_data%steps(i)%is_valid)THEN
            is_valid_byte_arr(i)=1_1;ts_char_arr(i)=traj_data%steps(i)%timestamp_str
            dom_id_arr(i)=traj_data%steps(i)%domain_id;sd_id_0based_arr(i)=traj_data%steps(i)%sd_id_in_file
            x_arr(i)=traj_data%steps(i)%x_coord;y_arr(i)=traj_data%steps(i)%y_coord;z_arr(i)=traj_data%steps(i)%z_coord
            r_arr(i)=traj_data%steps(i)%radius;n_arr(i)=traj_data%steps(i)%multiplicity;if_coal_arr(i)=traj_data%steps(i)%if_coal_flag
            parent_sdid_1based_arr(i)=traj_data%steps(i)%parent_sdid_for_trace;parent_dmid_arr(i)=traj_data%steps(i)%parent_dmid_for_trace
            coal_num_arr(i)=traj_data%steps(i)%coalescence_num_col;other_dmid_arr(i)=traj_data%steps(i)%other_parent_prev_dmid
            other_sdid_1based_arr(i)=traj_data%steps(i)%other_parent_prev_sdid;other_r_arr(i)=traj_data%steps(i)%other_parent_radius
            other_n_arr(i)=traj_data%steps(i)%other_parent_multiplicity
        ELSE
            is_valid_byte_arr(i)=0_1;ts_char_arr(i)=REPEAT(' ',TIME_STR_LEN)
            dom_id_arr(i)=FILL_INT;sd_id_0based_arr(i)=FILL_INT
            x_arr(i)=FILL_DOUBLE;y_arr(i)=FILL_DOUBLE;z_arr(i)=FILL_DOUBLE
            r_arr(i)=FILL_DOUBLE;n_arr(i)=FILL_INT64;if_coal_arr(i)=FILL_SHORT
            parent_sdid_1based_arr(i)=FILL_INT;parent_dmid_arr(i)=FILL_INT
            coal_num_arr(i)=FILL_INT;other_dmid_arr(i)=FILL_INT;other_sdid_1based_arr(i)=FILL_INT
            other_r_arr(i)=FILL_DOUBLE;other_n_arr(i)=FILL_INT64
        END IF
    END DO

    ! err=nf90_put_var(ncid_out,varid_is_valid_step,is_valid_byte_arr,start=(/1,current_record_idx/),count=(/num_time_steps_possible,1/))
    ! IF(err/=nf90_noerr)THEN 
    !     CALL handle_nc_error(err,"PVS")
    ! END IF
    ! DO i=1,num_time_steps_possible
    !     err=nf90_put_var(ncid_out,varid_ts_at_step,TRIM(ts_char_arr(i)),start=(/1,i,current_record_idx/),count=(/TIME_STR_LEN,1,1/))
    !     IF(err/=nf90_noerr)THEN 
    !         CALL handle_nc_error(err,"PVT_slice_" //achar(i+iachar('0'))) 
    !     END IF
    ! END DO
    err=nf90_put_var(ncid_out,varid_dom_at_step,dom_id_arr,start=(/1,current_record_idx/),count=(/num_time_steps_possible,1/))
    IF(err/=nf90_noerr)THEN 
        CALL handle_nc_error(err,"PVD")
    END IF
    err=nf90_put_var(ncid_out,varid_sdid_at_step,sd_id_0based_arr,start=(/1,current_record_idx/),count=(/num_time_steps_possible,1/))
    IF(err/=nf90_noerr)THEN 
        CALL handle_nc_error(err,"PVSID")
    END IF
    err=nf90_put_var(ncid_out,varid_x,x_arr,start=(/1,current_record_idx/),count=(/num_time_steps_possible,1/))
    IF(err/=nf90_noerr)THEN 
        CALL handle_nc_error(err,"PVX")
    END IF
    err=nf90_put_var(ncid_out,varid_y,y_arr,start=(/1,current_record_idx/),count=(/num_time_steps_possible,1/))
    IF(err/=nf90_noerr)THEN 
        CALL handle_nc_error(err,"PVY")
    END IF
    err=nf90_put_var(ncid_out,varid_z,z_arr,start=(/1,current_record_idx/),count=(/num_time_steps_possible,1/))
    IF(err/=nf90_noerr)THEN 
        CALL handle_nc_error(err,"PVZ")
    END IF
    err=nf90_put_var(ncid_out,varid_r,r_arr,start=(/1,current_record_idx/),count=(/num_time_steps_possible,1/))
    IF(err/=nf90_noerr)THEN 
        CALL handle_nc_error(err,"PVR")
    END IF
    err=nf90_put_var(ncid_out,varid_n,n_arr,start=(/1,current_record_idx/),count=(/num_time_steps_possible,1/))
    IF(err/=nf90_noerr)THEN 
        CALL handle_nc_error(err,"PVN")
    END IF
    err=nf90_put_var(ncid_out,varid_if_coal,if_coal_arr,start=(/1,current_record_idx/),count=(/num_time_steps_possible,1/))
    IF(err/=nf90_noerr)THEN 
        CALL handle_nc_error(err,"PVIC")
    END IF
    err=nf90_put_var(ncid_out,varid_parent_sdid,parent_sdid_1based_arr,start=(/1,current_record_idx/),count=(/num_time_steps_possible,1/))
    IF(err/=nf90_noerr)THEN 
        CALL handle_nc_error(err,"PVPSID")
    END IF
    err=nf90_put_var(ncid_out,varid_parent_dmid,parent_dmid_arr,start=(/1,current_record_idx/),count=(/num_time_steps_possible,1/))
    IF(err/=nf90_noerr)THEN 
        CALL handle_nc_error(err,"PVPDMID")
    END IF
    err=nf90_put_var(ncid_out,varid_coal_num_col,coal_num_arr,start=(/1,current_record_idx/),count=(/num_time_steps_possible,1/))
    IF(err/=nf90_noerr)THEN 
        CALL handle_nc_error(err,"PVCNC")
    END IF
    err=nf90_put_var(ncid_out,varid_other_dmid,other_dmid_arr,start=(/1,current_record_idx/),count=(/num_time_steps_possible,1/))
    IF(err/=nf90_noerr)THEN 
        CALL handle_nc_error(err,"PVODM")
    END IF
    err=nf90_put_var(ncid_out,varid_other_sdid,other_sdid_1based_arr,start=(/1,current_record_idx/),count=(/num_time_steps_possible,1/))
    IF(err/=nf90_noerr)THEN 
        CALL handle_nc_error(err,"PVOSDID")
    END IF
    err=nf90_put_var(ncid_out,varid_other_r,other_r_arr,start=(/1,current_record_idx/),count=(/num_time_steps_possible,1/))
    IF(err/=nf90_noerr)THEN 
        CALL handle_nc_error(err,"PVOR")
    END IF
    err=nf90_put_var(ncid_out,varid_other_n,other_n_arr,start=(/1,current_record_idx/),count=(/num_time_steps_possible,1/))
    IF(err/=nf90_noerr)THEN 
        CALL handle_nc_error(err,"PVON")
    END IF
  END SUBROUTINE write_trajectory_to_netcdf

END PROGRAM ParticleTracer

