# Introduction to the Super-Droplet Method (SDM) in SCALE-SDM

## What is SDM?
The **Super-Droplet Method (SDM)**, originally introduced by **Shima et al. (2009)**, represents a significant advancement in Lagrangian cloud microphysics simulations. Unlike traditional Eulerian microphysics schemes, SDM employs “super-droplets,” which are computational particles representing a large number of real droplets, aerosols, or precipitation particles with similar properties. This method allows for detailed simulation of key cloud microphysical processes such as condensation, evaporation, and coalescence, enabling an improved understanding of cloud formation and precipitation dynamics.

## Integration with SCALE-SDM
The **SCALE (Scalable Computing for Advanced Library and Environment)** framework, developed with co-design by computational and computer science researchers, provides a robust, scalable platform for high-resolution atmospheric simulations (Nishizawa et al., 2015; Sato et al., 2015). SCALE-SDM integrates the SDM within this framework, offering a powerful tool for simulating cloud microphysics and dynamics. For more information on the SCALE platform, please visit the [SCALE's official website](http://scale.aics.riken.jp/).

This repository integrates SDM into SCALE version 5.2.6, leveraging both SDM’s microphysical precision and SCALE’s computational scalability to facilitate simulations of cloud systems that capture the complex interactions between microphysics and atmospheric dynamics. Further details about this version of SCALE can be explored in the [SCALE version 5.2.6 archives](https://scale.riken.jp/archives/5.2.6/).

# Super-Droplet Forward Tracking Method

## Overview
The **forward tracking method** allows for the tracking of a selected subset of super-droplets (SDs) from the initial stage of the simulation, focusing on SDs that meet certain size criteria within specific altitude ranges. This method is particularly suited for high-resolution simulations with a large number of SDs, as it enables efficient data storage and extraction of relevant microphysical statistics while maintaining computational feasibility.

## Methodology
1. **Random Selection of SDs**
At the beginning of the simulation, a random selection of SDs with radii larger than a specified threshold is made within each domain across a given altitude range. The number of SDs selected at each altitude level is proportional to the number of SDs at that level that meet the radius criterion, relative to the total number of eligible SDs within the specified altitude range. The actual number of selected SDs is always less than or equal to the pre-set total number of SDs to be selected.

2. **Initialization of Selected SD IDs**
- **dm_id:** The domain ID of the selected SD at the initial simulation time.
- **sd_id:** The sequential identifier of the selected SD within the domain at the initial time.
  
These IDs remain unchanged throughout the simulation.

3. **Output Results**
Currently, output is only supported in netCDF format for selected SDs. The output filename begins with `SD_selected_NetCDF_*` and is controlled by setting `sdm_dmpvar=100` in [run.conf](https://github.com/wangyouyue/SDM_product-230225_Yin-SD_tracking/blob/SDM_selected_SD/scale-rm/test/case/shallowcloud/SD_tracking/run.conf#L149). The output variables include:
- **dm_id** and **sd_id** of the selected SDs.
- **if_coal:** A flag indicating whether collision-coalescence occurred during the previous microphysical time step (0 = no, 1 = yes).

4. **Collision-Coalescence Algorithm**
The collision-coalescence algorithm is similar to that used in backward tracking (refer to [the backward tracking branch](https://github.com/wangyouyue/SDM_product-230225_Yin-SD_tracking). The IDs of the SDs involved in a collision-coalescence event remain unchanged before and after the event. Only collision-coalescence events involving the selected SDs are recorded. If one of the colliding SDs is not selected, its `dm_id1/dm_id2` and `sd_id1/sd_id2` are set to a default value of `-990`. Additionally, the file records:
- **sd_r1/sd_r2:** The radii of the two SDs before the collision-coalescence.
- **sd_n1/sd_n2:** The multiplicities of the two SDs before the collision-coalescence.
- **num_col:** The number of coalescence between the pair of SDs.

These collision event files are **output at each microphysical time step**. If no collision-coalescence occurs or if no selected SD is involved in a collision-coalescence event, no `SD_coal_output_NetCDF*` file is generated.

5. **Random Perturbations in SD Motion**
Random perturbations are applied to the displacement of each SD during its motion. The magnitude of the displacement perturbation is proportional to the square root of the time step and inversely proportional to the square root of the SD radius. The displacement formula is:

$$\Delta x = (\text{random number} - 0.5) \times \sqrt{\frac{\Delta t}{r}} \times \text{sdm\_noise\_amp}$$

where the **random number** is between 0 and 1 (the same for all three directions), $\Delta x$ is the displacement, $\Delta t$ is the time step, $r$ is the SD radius, and **sdm_noise_amp** is the user-specified noise amplitude.

6. **Parameter Configuration in run.conf**
For parameter configuration in the `run.conf` file (see [run.conf](https://github.com/wangyouyue/SDM_product-230225_Yin-SD_tracking/blob/SDM_selected_SD/scale-rm/test/case/shallowcloud/SD_tracking/run.conf#L160-L165)):

- **num_selected:** Number of super-droplets per domain to select.
- **height_min and height_max:** Minimum and maximum heights for selection [m].
- **radius_min:** Minimum radius for selection [m].
- **coal_output:** Control flag to output collision-coalescence events (0: off, 1: on).
- **sdm_noise_amp:** Amplitude of random noise [m^1.5 * s^-0.5], default value is 0.D0.

## Advantages and Disadvantages of the Forward Tracking Method
The forward tracking method allows for the random sampling of SDs from specific altitude ranges at the beginning of the simulation, and tracks their future trajectories. This approach is ideal for high-resolution simulations with a large number of SDs, as it reduces storage usage while providing microphysical statistics (e.g., mean radius and relative dispersion at each altitude) with a confidence level close to the entire population.

However, one potential drawback is that the spatial distribution of particles may become uneven over long simulation durations, which could compromise the accuracy of the statistical information. Additionally, the computational resources required for post-processing SD information are higher than for the backward tracking method, and not all SD information is retained.

# Installation and Usage
For more details, please see [the official user guide of SACLE](https://scale.riken.jp/archives/scale_users_guide_En.v5.2.6.pdf).

## Prerequisites
- **Compilers:** Fortran and C compilers are required.
- **Libraries:** MPI, NetCDF4, and HDF5 libraries must be installed.

## Setting Up SCALE-SDM
1. **Environment Preparation:**
`$ export SCALE_SYS=Linux64-intel-impi`

2. **Clean:**
   Please clean up past compilation files before compiling.
  ```
  $ cd scale-rm/test/case/shallowcloud/dycoms2_rf02_sdm_hokudai/
  $ make allclean
  $ make allclean SCALE_ENABLE_SDM=T SCALE_DISABLE_LOCALBIN=T SCALE_DYCOMS2_RF02_SDM=T
  ```

3. **Compilation:**
   Navigate to the desired test case directory and compile using the provided Makefile. For SDM-specific cases, ensure the appropriate flags are set.
  ```
  $ cd scale-rm/test/case/shallowcloud/dycoms2_rf02_sdm_hokudai/
  $ make SCALE_ENABLE_SDM=T SCALE_DISABLE_LOCALBIN=T SCALE_DYCOMS2_RF02_SDM=T
  $ ln -fsv  `grep ^TOPDIR Makefile | sed s/\)//g | awk '{print $NF}'`/bin/scale-rm* .
  ```

4. **Running Simulations:**
   Modify and submit the job script according to your system's job scheduler. For detailed execution steps and analysis programs for SDM, refer to the provided scripts and the SCALE user guide.
   - **Run a batch job on Hokudai supercomputer:**
     Modify the job script according to your system. The job scheduler on Hokudai supercomputer is [PJM](https://www.hucc.hokudai.ac.jp/en_supercomputer/basic/en_job_execution/). Run `$ pjsub hokudai_run.sh` to submit the job.

5. **Running Analysis Program:**
  ```
  $ cd scale-rm/test/case/shallowcloud/dycoms2_rf02_sdm_hokudai/
  $ pjsub --step --sparam "sn=1" ncl.sh
  $ pjsub --step --sparam "jid=JOB_ID, sn=2, sd=ec!=0:after:1" merge.sh
  ```

## Running the 2D Test on Supercomputer at University of Hyogo
1. **Switching the Branch:**
  `$ git checkout SDM_seleted_SD`
  
2. **Clean:**
   Please clean up past compilation files before compiling.
  ```
  $ cd scale-rm/test/case/shallowcloud/dycoms2_rf02_sdm_2D_forward
  $ module purge
  $ module load intel/2022.3.1 mpt hdf5/1.14.3 netcdf-c/4.9.2 netcdf-fortran/4.6.1
  $ make allclean
  $ make allclean SCALE_ENABLE_SDM=T SCALE_DISABLE_LOCALBIN=T SCALE_DYCOMS2_RF02_SDM=T
  ```

3. **Compilation:**
  ```
  $ cd scale-rm/test/case/shallowcloud/dycoms2_rf02_sdm_2D_forward
  $ make SCALE_ENABLE_SDM=T SCALE_DISABLE_LOCALBIN=T SCALE_DYCOMS2_RF02_SDM=T
  $ ln -fsv  `grep ^TOPDIR Makefile | sed s/\)//g | awk '{print $NF}'`/bin/scale-rm* .
  ```

4. **Running Simulations:**
  `$ qsub UoH_run.pbs`

5. **Running Analysis Program:**
  Before submitting the job, ensure that the output interval of super-droplets (`TIME_STEP_INTERVAL` in seconds) in the [Python script](https://github.com/wangyouyue/SDM_product-230225_Yin-SD_tracking/blob/SDM_selected_SD/scale-rm/test/case/shallowcloud/dycoms2_rf02_sdm_2D_forward/results/sd_output.py) `sd_output.py` match the settings in `init.conf` and `run.conf`. Additionally, note that `sdm_dmpitvl` (the time interval for binary output of selected droplets) in `run.conf` is also specified in seconds. You can adjust the processing duration in the Python script by modifying `time_str` and `end_time_str` (in the format “HHMMSS.sss”). Since this is a forward tracking process, `time_str` should be earlier than `end_time_str`. Finally, `num_processes` represents the number of parallel processes, so ensure it is consistent with the job script `run_py.pbs`.
  ```
  $ cd scale-rm/test/case/shallowcloud/dycoms2_rf02_sdm_2D_forward/results
  $ qsub run_py.pbs
  ```

  ## Running the 2D Simulation on SQUID Supercomputer at Osaka University

1. **Compilation using the Intel compilers and Intel MPI:**
   
  - Set the SCALE system environment variable for your session:
    ` $ export SCALE_SYS=Linux64-intel-impi`
  - Purge existing modules:
    ` $ module purge` 
  - Load the required modules:
    ` $ module load intel/2022.3.1 mpt hdf5/1.14.3 netcdf-c/4.9.2 netcdf-fortran/4.6.1`
  - After completing these steps, your environment should be ready. You can then compile SCALE-SDM:
    ```
    $ cd scale-rm/test/case/shallowcloud/dycoms2_rf02_sdm_2D_forward/results
    $ make allclean
    $ make allclean SCALE_ENABLE_SDM=T SCALE_DISABLE_LOCALBIN=T SCALE_DYCOMS2_RF02_SDM=T
    $ make SCALE_ENABLE_SDM=T SCALE_DISABLE_LOCALBIN=T SCALE_DYCOMS2_RF02_SDM=T
    $ ln -fsv  `grep ^TOPDIR Makefile | sed s/\)//g | awk '{print $NF}'`/bin/scale-rm* .
    ```
  
  2. **Submit the job:**
    `$ qsub squid_run.sh`

## Support and Community
Questions, issues, and discussions about SCALE-SDM can be directed here. Contributions and feedback are highly encouraged to enhance the model's capabilities and user experience. Please feel free to contact me: yinchongzhi@gmail.com. :grin:

## Reference
*Nishizawa, S., Yashiro, H., Sato, Y., Miyamoto, Y., and Tomita, H.: Influence of grid aspect ratio on planetary boundary layer turbulence in large-eddy simulations, Geoscientific Model Development, 8, 3393-3419, [https://doi.org/10.5194/gmd-8-33932015](https://doi.org/10.5194/gmd-8-33932015), 2015.*

*Sato, Y., Nishizawa, S., Yashiro, H., Miyamoto, Y., Kajikawa, Y., and Tomita, H.: Impacts of cloud microphysics on trade wind cumulus: which cloud microphysics processes contribute to the diversity in a large eddy simulation?, Progress in Earth and Planetary Science, 2, 1-16, [https://doi.org/10.1186/s40645-015-0053-6](https://doi.org/10.1186/s40645-015-0053-6), 2015.*

*Shima, S.-i., Kusano, K., Kawano, A., Sugiyama, T., and Kawahara, S.: The super-droplet method for the numerical simulation of clouds and precipitation: A particle-based and probabilistic microphysics model coupled with a non-hydrostatic model, Quarterly Journal of the Royal Meteorological Society, 135, 1307-1320, [https://doi.org/10.1002/qj.441](https://doi.org/10.1002/qj.441), 2009.*
