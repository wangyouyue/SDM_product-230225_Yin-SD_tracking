# Introduction to the Super-Droplet Method (SDM) in SCALE-SDM

## What is SDM?
The **Super-Droplet Method (SDM)**, originally introduced by **Shima et al. (2009)**, represents a significant advancement in Lagrangian cloud microphysics simulations. Unlike traditional Eulerian microphysics schemes, SDM employs “super-droplets,” which are computational particles representing a large number of real droplets, aerosols, or precipitation particles with similar properties. This method allows for detailed simulation of key cloud microphysical processes such as condensation, evaporation, and coalescence, enabling an improved understanding of cloud formation and precipitation dynamics.

## Integration with SCALE-SDM
The **SCALE (Scalable Computing for Advanced Library and Environment)** framework, developed with co-design by computational and computer science researchers, provides a robust, scalable platform for high-resolution atmospheric simulations (Nishizawa et al., 2015; Sato et al., 2015). SCALE-SDM integrates the SDM within this framework, offering a powerful tool for simulating cloud microphysics and dynamics. For more information on the SCALE platform, please visit the [SCALE's official website](http://scale.aics.riken.jp/).

This repository integrates SDM into SCALE version 5.2.6, leveraging both SDM’s microphysical precision and SCALE’s computational scalability to facilitate simulations of cloud systems that capture the complex interactions between microphysics and atmospheric dynamics. Further details about this version of SCALE can be explored in the [SCALE version 5.2.6 archives](https://scale.riken.jp/archives/5.2.6/).

## SD Tracking in SCALE-SDM
To enhance the analytical capabilities of the SCALE-SDM model, we have implemented a **Super-Droplet (SD)** backward tracking algorithm, which enables the tracing of super-droplet lifecycles and interactions throughout the simulation. This feature is essential for investigating the intricate microphysical processes governing cloud and precipitation formation. Specifically, the tracking mechanism records the “previous” state (`pre_id` and `pre_dmid`) of each SD at user-defined output intervals, thus facilitating detailed analysis of SD trajectories and their interactions over time.

# Super-Droplet Backward Tracking Method

## Overview
The SD tracking method is designed to trace the lifecycle and interactions of super-droplets within cloud microphysics simulations. By capturing the evolution of SDs, including their growth, transport, and interaction, this method provides a unique perspective on the underlying microphysical processes in clouds.

## Methodology
At each output interval, two key variables are recorded for each SD:

- **pre_dmid:** The domain identifier, indicating the simulation domain of the SD at the previous output interval.
- **pre_id:** The unique identifier of the SD in the simulation domain at the previous output interval.
  
This data allows backward tracing of SDs’ paths through the simulation, providing insights into the microphysical events and processes each SD undergoes. By analyzing these trajectories, researchers can better understand the microphysical history of cloud droplets, such as the processes of condensation, evaporation, and collision-coalescence.

## Advantages
1. **Efficiency in post-processing:** The method efficiently tracks super-droplets without the need for continuous tracking. By recording only at specific intervals, it reduces the computational and storage overhead compared to other methods.
2. **Simplicity in implementation:** The tracking method integrates seamlessly with existing cloud microphysics frameworks, making it accessible for researchers working with large-scale atmospheric simulations.
  
## Limitations and Proposed Solutions
1. **Limited forward tracing capability:** Since the current implementation is primarily designed for backward tracing, forward tracking of SDs can become cumbersome.
  - **Solution:** A forward tracking algorithm has been developed as a separate branch, enabling the tracking of SDs’ future paths based on their current state. For more details on forward tracking, please visit [this branch](https://github.com/wangyouyue/SDM_product-230225_Yin-SD_tracking/tree/SDM_selected_SD).
2. **Potential missed microphysical events:** Due to the intervals between outputs, important microphysical events may not be captured in real time.
  - **Solution A:** Increase the output frequency to capture more detailed interactions, though this will increase the data volume.
  - **Solution B (currently implemented):** Record key microphysical events, such as coalescence, in real time, ensuring that significant events are not missed. For further details, see [Coalescence Event Tracking and Output](#coalescence-event-tracking-and-output).
  - **Solution C:** Employ event-driven outputs, where specific changes in SD properties trigger additional data output.

## **Output Methods and Compatibility**
TThe SD tracking method outputs simulation results in several formats:
- **sdm_outasci:** Outputs all SDs in ASCII format.
- **sdm_outnetcdf:** Outputs all SDs in netCDF format.
- **sdm_copy_selected_sd:** Outputs selected SDs in netCDF format (currently incompatible with the SD tracking feature).

   Users must choose between `sdm_outasci` and `sdm_outnetcdf`, as both cannot be used simultaneously. To address this, users can set the `sdm_dmpvar` option in the SDM settings section (`&PARAM_ATMOS_PHY_MP_SDM`) of the run.conf namelist file to either `010` or `001`.

## Conclusion
The **Super-Droplet tracking method** provides an essential tool for gaining detailed insights into the dynamics and interactions of cloud droplets in Lagrangian cloud microphysics simulations. By allowing researchers to trace individual SDs and record key microphysical events, this method facilitates the investigation of processes such as condensation, evaporation, and collision-coalescence, which are critical to understanding cloud development and precipitation formation. While the method is computationally efficient and straightforward to implement, challenges remain in capturing fine-scale microphysical processes and enabling efficient forward tracking. The proposed solutions, such as event-driven outputs and real-time coalescence event recording, aim to mitigate these limitations, providing a more comprehensive and detailed framework for SD tracking.

## Future Work
Future advancements in the SD tracking framework could involve the development of more sophisticated event detection algorithms that trigger outputs based on significant changes in the SD states, such as rapid growth or collision events. This could enhance the resolution and accuracy of the recorded microphysical processes without a substantial increase in data volume. Additionally, improving the real-time analysis capabilities of the model would allow for on-the-fly tracking and analysis of droplet interactions. Another potential avenue for improvement is optimizing the trade-off between detailed tracking and computational efficiency, particularly for large-scale, long-duration cloud simulations, ensuring that the model remains scalable while providing accurate microphysical insights.

# Coalescence Event Tracking and Output

## Coalescence Event Tracking
To extend the analytical capabilities of the **SCALE-SDM model**, we have integrated a robust feature for recording **coalescence events** between super-droplets (SDs). This enhancement is pivotal for understanding the evolution and interactions of cloud particles, particularly during collision-coalescence processes that play a fundamental role in precipitation formation. By explicitly recording the microphysical details of coalescence, researchers can gain deeper insights into the dynamics that lead to cloud droplet growth and the initiation of rainfall.

## Methodology
The tracking of coalescence events follows a structured approach designed to capture the key interactions between SDs:

1. **Event Identification:** At each simulation time step, potential coalescence events between pairs of SDs are identified based on their calculated collision-coalescence probability. This probability is derived from the physical properties of the SDs, such as size and relative velocity, as well as the environmental conditions within the simulation domain.
2. **Event Recording:** For each coalescence event, the following information is recorded:
- **pre_sdid1 and pre_sdid2:** The unique identifiers of the colliding SDs before the coalescence event.
- **pre_dmid1 and pre_dmid2:** The domain identifiers of the colliding SDs, indicating their location within the simulation domain.
- **sd_r1/sd_r2:** The radii of the two SDs before the collision-coalescence.
- **sd_n1/sd_n2:** The multiplicities of the two SDs before the collision-coalescence.
- **num_col:** The number of coalescence events occurring between the pair of SDs within a coalescence time step.
3. **Dynamic Data Management:** Given the stochastic nature of coalescence events, the model employs dynamic memory allocation to manage the varying number of coalescence events that occur during each simulation time step. This approach ensures computational efficiency and prevents unnecessary memory allocation when coalescence events are infrequent.
4. **Output Generation:** The coalescence data is output using the `sdm_coal_outnetcdf` subroutine, which writes detailed information about each event to **NetCDF** files (only NetCDF format is supported currently). This output format includes not only the identifiers and domain information of the colliding SDs but also the frequency of their interactions, enabling post-simulation analysis of the collision dynamics in a standardized and accessible format.

## Integration with SDM Introduction
The coalescence event tracking feature is seamlessly integrated into the broader SDM framework, augmenting its ability to simulate cloud microphysical processes. This addition provides researchers with a powerful tool to investigate the growth mechanisms of cloud droplets, particularly through the lens of collision-coalescence, which is crucial for understanding precipitation initiation. The tracking and output of coalescence events complement the SDM’s existing capabilities, enabling a more detailed examination of cloud dynamics and the processes that lead to precipitation.

## Output File Naming Convention
The output files generated by the `sdm_coal_outnetcdf` subroutine adhere to a structured naming convention that facilitates easy identification and organization. Each output file name includes:

- **SD_coal_output:** Indicates that the file contains coalescence event data.
- **Timestamp (YYYYMMDD-HHMMSS.sss):** Reflects the exact time in the simulation when the coalescence data was recorded.
- **Process Identifier:** Ensures that each file is uniquely associated with a specific simulation instance, providing a clear reference for subsequent analysis.

This systematic naming convention allows researchers to efficiently manage and analyze the large volumes of data generated by high-resolution cloud microphysics simulations.

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

5. **Running analysis program:**
  ```
  $ cd scale-rm/test/case/shallowcloud/dycoms2_rf02_sdm_hokudai/
  $ pjsub --step --sparam "sn=1" ncl.sh
  $ pjsub --step --sparam "jid=JOB_ID, sn=2, sd=ec!=0:after:1" merge.sh
  ```
  
## Running the 2D Test on Supercomputer at University of Hyogo
1. **Clean:**
  ```
  $ cd scale-rm/test/case/shallowcloud/dycoms2_rf02_sdm_2D_backward
  $ module purge
  $ module load intel/2022.3.1 mpt hdf5/1.14.3 netcdf-c/4.9.2 netcdf-fortran/4.6.1
  $ make allclean
  $ make allclean SCALE_ENABLE_SDM=T SCALE_DISABLE_LOCALBIN=T SCALE_DYCOMS2_RF02_SDM=T
  ```

2. **Compilation:**
  ```
  $ cd scale-rm/test/case/shallowcloud/dycoms2_rf02_sdm_2D_backward
  $ make SCALE_ENABLE_SDM=T SCALE_DISABLE_LOCALBIN=T SCALE_DYCOMS2_RF02_SDM=T
  $ ln -fsv  `grep ^TOPDIR Makefile | sed s/\)//g | awk '{print $NF}'`/bin/scale-rm* .
  ```

3. **Running Simulations:**
  `$ qsub UoH_run.sh`

4. **Running analysis program:**
  Before submitting the job, ensure that the grid resolution (`DX`, `DY`, and `DZ`) and the output interval of super-droplets (`TIME_STEP_INTERVAL` in milliseconds) in the [Python script](https://github.com/wangyouyue/SDM_product-230225_Yin-SD_tracking/tree/SDM_SD_tracking/scale-rm/test/case/shallowcloud/dycoms2_rf02_sdm_2D_backward/results/sd_output.py) `sd_output.py` match the settings in `init.conf` and `run.conf`. Additionally, note that `sdm_dmpitvb` (the time interval for binary output of all droplets) in `run.conf` is specified in seconds.

  You can adjust the processing duration in the Python script by modifying `start_time` and `end_time` (in the format “HHMMSS.sss”). Since this is a backward tracking process, `end_time` should be earlier than `start_time`. Finally, `num_blocks` represents the number of parallel processes, so ensure it is consistent with the job script `run_py.pbs`.
  ```
  $ cd scale-rm/test/case/shallowcloud/dycoms2_rf02_sdm_2D_backward/results
  $ qsub run_py.pbs
  ```

## Support and Community
Questions, issues, and discussions about SCALE-SDM can be directed here. Contributions and feedback are highly encouraged to enhance the model's capabilities and user experience. Please feel free to contact me: yinchongzhi@gmail.com. :grin:

## Reference
*Nishizawa, S., Yashiro, H., Sato, Y., Miyamoto, Y., and Tomita, H.: Influence of grid aspect ratio on planetary boundary layer turbulence in large-eddy simulations, Geoscientific Model Development, 8, 3393-3419, [https://doi.org/10.5194/gmd-8-33932015](https://doi.org/10.5194/gmd-8-33932015), 2015.*

*Sato, Y., Nishizawa, S., Yashiro, H., Miyamoto, Y., Kajikawa, Y., and Tomita, H.: Impacts of cloud microphysics on trade wind cumulus: which cloud microphysics processes contribute to the diversity in a large eddy simulation?, Progress in Earth and Planetary Science, 2, 1-16, [https://doi.org/10.1186/s40645-015-0053-6](https://doi.org/10.1186/s40645-015-0053-6), 2015.*

*Shima, S.-i., Kusano, K., Kawano, A., Sugiyama, T., and Kawahara, S.: The super-droplet method for the numerical simulation of clouds and precipitation: A particle-based and probabilistic microphysics model coupled with a non-hydrostatic model, Quarterly Journal of the Royal Meteorological Society, 135, 1307-1320, [https://doi.org/10.1002/qj.441](https://doi.org/10.1002/qj.441), 2009.*
