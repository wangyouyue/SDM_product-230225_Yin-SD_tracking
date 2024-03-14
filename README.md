# Introduction to the Super Droplet Method (SDM) in SCALE-SDM

## What is SDM?
The Super Droplet Method (SDM), introduced by Shima et al. in 2009, offers a novel Lagrangian approach to cloud microphysics simulations. Unlike traditional schemes, SDM tracks "super droplets," each representing a collection of real particles (aerosols, cloud droplets and precipitation particles) with similar attributes. This method enables detailed simulations of microphysical processes such as condensation, evaporation, and coalescence, providing insights into cloud formation and precipitation dynamics.

## Integration with SCALE-SDM
The SCALE (Scalable Computing for Advanced Library and Environment) library is developed with co-design by researchers of computational science and computer science, ensuring a robust and scalable platform for high-resolution atmospheric simulations (Nishizawa et al., 2015; Sato et al., 2015).  For more information about SCALE and its features, visit the  [SCALE's official website](http://scale.aics.riken.jp/). The specific implementation of SDM into SCALE version 5.2.6 can be explored in the [SCALE version 5.2.6 archives](https://scale.riken.jp/archives/5.2.6/)

We have incorporated the SDM into the SCALE framework, version 5.2.6. This integration, known as SCALE-SDM, leverages the strengths of both SDM and SCALE to offer a powerful tool for studying cloud microphysics and dynamics. The SCALE-SDM model facilitates high-fidelity simulations of cloud systems, capturing the complex interplay of microphysical processes and atmospheric dynamics.

## SD Tracking in SCALE-SDM
To further enhance the capabilities of SCALE-SDM, we have introduced a Super Droplet (SD) tracking feature. This addition allows researchers to trace the lifecycle and interactions of super droplets within the simulation, providing deeper insights into the microphysical processes driving cloud and precipitation development. The tracking feature records the "previous" state (pre_id and pre_dmid) of each SD at specified output intervals, enabling detailed analysis of SD trajectories and interactions over time.

# Super Droplet Tracking Method

## Overview
The Super Droplet tracking method is designed to monitor the lifecycle and interactions of super droplets within cloud microphysics simulations. By capturing the dynamics of SDs, including their formation, growth, and interactions, this method provides a unique window into the microphysical processes at play within clouds.

## Methodology
At each output interval, the method records two key pieces of information for each SD:

- **pre_id:** The identifier of the SD at the previous output interval.
- **pre_dmid:** The domain identifier, indicating the simulation domain of the SD at the previous output interval.
  
This information allows for backward tracing of each SD's path through the simulation, providing insights into the microphysical processes and interactions that each droplet undergoes.

## Advantages
1. **Simplicity:** The method is straightforward to implement and integrate into existing cloud microphysics simulation frameworks.
2. **Efficiency:** By only recording data at output intervals, the method minimizes the computational overhead and storage requirements compared to continuous tracking.
  
## Disadvantages and Solutions
1. **Backward Tracing Complexity:** Forward tracing of SDs becomes cumbersome, as the method is inherently designed for backward tracing.
  - **Solution:** Implement real-time tracking of each SD's trajectory, recording interactions and movements at every time step, albeit at a higher computational and storage cost.
2. **Missed Microphysical Processes:** The output interval may skip over significant microphysical events due to its length being typically much greater than the time steps of microphysical processes.
  - **Solution A:** Increase the output frequency to capture more detailed microphysical interactions, balancing against the increased data volume.
  - **Solution B (what we currently use):** Record key events (e.g., coalescence) as they occur, in addition to the regular output, to ensure significant changes are not missed. For more information, see [Coalescence Event Tracking and Output](#coalescence-event-tracking-and-output).
  - **Solution C:** Use event-driven outputs that trigger based on specific changes in SD states, offering a compromise between data volume and detail.

## **Output Methods and Compatibility**
The model supports three output formats for SD results:
- **sdm_outasci:** Outputs all SDs in ASCII format.
- **sdm_outnetcdf:** Outputs all SDs in netCDF format.
- **sdm_copy_selected_sd:** Outputs selected SDs in netCDF format. However, this method is currently not compatible with the SD tracking feature.

   Additionally, users must choose between `sdm_outasci` and `sdm_outnetcdf` as these cannot be used simultaneously. To address this, users can set the `sdm_dmpvar` option in the SDM settings section (&PARAM_ATMOS_PHY_MP_SDM) of the run.conf namelist file to either `010` or `001`.

## Conclusion
The SD tracking method provides a valuable tool for understanding the dynamics and interactions of super droplets in cloud microphysics simulations. While it offers simplicity and efficiency, the method also faces challenges in capturing detailed microphysical processes and facilitating forward tracing. The proposed solutions aim to address these challenges, enabling more comprehensive and detailed tracking of SDs.

## Future Work
Future enhancements could include developing more sophisticated event detection mechanisms, improving real-time analysis capabilities, and optimizing the balance between detail and computational efficiency.

# Coalescence Event Tracking and Output

## Coalescence Event Tracking
To extend the capabilities of SCALE-SDM further, we have integrated a feature for recording coalescence events between super droplets (SDs). This development is pivotal for understanding the evolution and dynamics of cloud particles through collision-coalescence processes, which are fundamental to precipitation formation.

## Methodology
The tracking of coalescence events involves several key steps:

1. **Event Identification:** At each time step, potential coalescence events between pairs of SDs are identified based on their collision-coalescence probability, which is calculated considering their physical properties and environmental conditions.
2. **Event Recording:** For each identified coalescence event, we record:
- **pre_sdid1 and pre_sdid2:** The identifiers of the colliding SDs before the coalescence event.
- **pre_dmid1 and pre_dmid2:** The domain identifiers of the colliding SDs, indicating their location within the simulation domain.
- **num_col:** The number of coalescence events occurring between the pair of SDs within a coalescence time step.
3. **Dynamic Data Management:** Given the unpredictable nature of coalescence events, dynamic allocation of arrays is employed to efficiently manage the varying number of coalescence events.
4. **Output Generation:** Utilizing the `sdm_coal_outnetcdf` subroutine, detailed information about each coalescence event is outputted in NetCDF format for further analysis. This output includes not only the identifiers and domain information of the colliding SDs but also the frequency of their interactions.

## Integration with SDM Introduction
This coalescence event tracking feature enriches the SDM framework by providing a comprehensive toolset for dissecting the microphysical interactions within clouds. Researchers can now delve deeper into the mechanisms of cloud particle formation and growth, gaining insights into the complexities of precipitation processes.

By offering this additional layer of analysis, SCALE-SDM enhances its utility as a state-of-the-art tool for cloud microphysics research, facilitating a closer examination of the intricate dance of droplets that leads to the myriad forms of precipitation observed in nature.

## Output File Naming Convention
The output files generated by the `sdm_coal_outnetcdf` subroutine follow a naming convention that incorporates the type of output (`SD_coal_output`), a timestamp, and the process identifier. This ensures that each output file is uniquely identifiable and related to its specific simulation context.

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

5. **Run analysis program:**
  ```
  $ cd scale-rm/test/case/shallowcloud/dycoms2_rf02_sdm_hokudai/
  $ pjsub --step --sparam "sn=1" ncl.sh
  $ pjsub --step --sparam "jid=JOB_ID, sn=2, sd=ec!=0:after:1" merge.sh
  ```

## Support and Community
Questions, issues, and discussions about SCALE-SDM can be directed here. Contributions and feedback are highly encouraged to enhance the model's capabilities and user experience. Please feel free to contact me: yinchongzhi@gmail.com. :grin:

## Reference
*Nishizawa, S., Yashiro, H., Sato, Y., Miyamoto, Y., and Tomita, H.: Influence of grid aspect ratio on planetary boundary layer turbulence in large-eddy simulations, Geoscientific Model Development, 8, 3393-3419, [https://doi.org/10.5194/gmd-8-33932015](https://doi.org/10.5194/gmd-8-33932015), 2015.*

*Sato, Y., Nishizawa, S., Yashiro, H., Miyamoto, Y., Kajikawa, Y., and Tomita, H.: Impacts of cloud microphysics on trade wind cumulus: which cloud microphysics processes contribute to the diversity in a large eddy simulation?, Progress in Earth and Planetary Science, 2, 1-16, [https://doi.org/10.1186/s40645-015-0053-6](https://doi.org/10.1186/s40645-015-0053-6), 2015.*

*Shima, S.-i., Kusano, K., Kawano, A., Sugiyama, T., and Kawahara, S.: The super-droplet method for the numerical simulation of clouds and precipitation: A particle-based and probabilistic microphysics model coupled with a non-hydrostatic model, Quarterly Journal of the Royal Meteorological Society, 135, 1307-1320, [https://doi.org/10.1002/qj.441](https://doi.org/10.1002/qj.441), 2009.*
