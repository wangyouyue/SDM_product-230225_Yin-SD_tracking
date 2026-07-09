!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Cloud Microphysics / SDM
!!
!! @par Description
!!          Input and Output of the SDM variables
!!
!! - Reference
!!  - Shima et al., 2009:
!!    The super-droplet method for the numerical simulation of clouds and precipitation:
!!    A particle-based and probabilistic microphysics model coupled with a non-hydrostatic model.
!!    Quart. J. Roy. Meteorol. Soc., 135: 1307-1320
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2014-06-27 (S.Shima) [new] sdm_outasci is added
!! @li      2016-07-11 (S.Shima) [mod] modified to support sdice output
!! @li      2018-06-25 (S.Shima) [add] netcdf output
!! @li      2018-06-30 (S.Shima) [add] rime mass and number of monomers as SD attributes
!! @li      2020-07-16 (S.Shima) [add] netcdf compression and shuffle options
!! @li      2020-07-17 (S.Shima) [add] sdm_outnetcdf_hist
!! @li      2020-07-23 (S.Shima) [mod] sdm_outnetcdf and sdm_outnetcdf_hist for sdm_dmpvar == 1?? and sdm_dmpvar == 2?? 
!! @li      2020-10-30 (S.Shima) [fix] sdm_outnetcdf and sdm_outnetcdf_hist for sdm_dmpvar == 1?? and sdm_dmpvar == 2?? 
!!
!<
!-------------------------------------------------------------------------------
module m_sdm_io

  implicit none
  private
  public :: sdm_outasci,sdm_outnetcdf,sdm_outnetcdf_hist,sdm_coal_outnetcdf
  public :: sdm_event_collision_outnetcdf
  public :: sdm_event_singleproc_outnetcdf,sdm_event_diag_outnetcdf
  public :: sdm_assign_tracking_subset,sdm_interest_id_outnetcdf
  public :: sdm_lifecycle_outnetcdf
  public :: sdm_tracking_lifecycle_write_domain_entry
  public :: sdm_tracking_lifecycle_write_global_halo_entry
  public :: sdm_tracking_lifecycle_write_seeding_entry

contains
  subroutine sdm_outasci(otime,sd_num,sd_numasl,sd_n,sd_liqice,sd_x,sd_y,sd_z,sd_r,sd_asl,sd_vz,sdi,sdn_dmpnskip)
    use scale_precision
    use scale_stdio
    use scale_time
    use scale_process, only: &
         mype => PRC_myrank, &
         PRC_MPIstop
    use m_sdm_common, only: &
         i2, sdm_cold, STAT_LIQ, STAT_ICE, STAT_MIX, sdicedef, &
         INVALID_i4, coalescence_output_enable

    implicit none

    real(DP), intent(in) :: otime
    integer, intent(in) :: sd_num    ! number of super-droplets
    integer, intent(in) :: sd_numasl ! number of chemical species contained in super droplets
    integer(DP), intent(in) :: sd_n(1:sd_num) ! multiplicity of super-droplets
    integer(kind=i2), intent(in) :: sd_liqice(1:sd_num)
                       ! status of super-droplets (liquid/ice)
                       ! 01 = all liquid, 10 = all ice
                       ! 11 = mixture of ice and liquid
    real(RP), intent(in) :: sd_x(1:sd_num) ! x-coordinate of super-droplets
    real(RP), intent(in) :: sd_y(1:sd_num) ! y-coordinate of super-droplets
    real(RP), intent(in) :: sd_z(1:sd_num) ! z-coordinate of super-droplets
    real(RP), intent(in) :: sd_r(1:sd_num) ! equivalent radius of super-droplets
    real(RP), intent(in) :: sd_asl(1:sd_num,1:sd_numasl) ! aerosol mass of super-droplets
    real(RP), intent(in) :: sd_vz(1:sd_num) ! terminal velocity of super-droplets
    type(sdicedef), intent(in) :: sdi   ! ice phase super-droplets
    integer, intent(in) :: sdn_dmpnskip ! Base skip to store super droplets in text format

    character(len=17) :: fmt2="(A, '.', A, I*.*)"
    character(len=17) :: fmt3="(3A)"
    character(len=H_LONG) :: ftmp
    character(len=H_LONG) :: basename_sd_out
    character(len=19) :: basename_time
    integer :: fid_sdm_o
    integer :: n, m, ierr
    character(len=80) :: fmt     ! output formate
    character(len=5)  :: cstat   ! status character
    
    !--- output Super Droplets in ASCII format
    call TIME_gettimelabel(basename_time)
    fid_sdm_o = IO_get_available_fid()

    write(fmt2(14:14),'(I1)') 6
    write(fmt2(16:16),'(I1)') 6
    write(ftmp,fmt3) 'SD_output', '_ASCII_', trim(basename_time)
    write(basename_sd_out,fmt2) trim(ftmp), 'pe',mype

    open (fid_sdm_o, file = trim(basename_sd_out), & !action = "write", &
          access = "sequential", status = "replace", form = "formatted", &
          iostat = ierr)

    if( ierr /= 0 ) then
      write(*,*) "sdm_ascii_out", "Write error"
      call PRC_MPIstop
    endif 

    if( .not.sdm_cold ) then

       write(fid_sdm_o,'(3a,i2.2,2a)') '# x[m],y[m],z[m],vz[m],',       &
            &                'radius(droplet)[m],',                        &
            &                'mass_of_aerosol_in_droplet(1:',sd_numasl,')[g],',&
            &                'multiplicity[-],status[-],index'

       write(fmt,'( "(", i2.2, "e16.8,i20,a5,i10)" )')(sd_numasl+5)

       do m=1,sd_num,sdn_dmpnskip

          cstat = '  LIQ'
          
          write(fid_sdm_o,trim(fmt)) sd_x(m),           &
               &                   sd_y(m),           &
               &                   sd_z(m), sd_vz(m), sd_r(m),        &
               &                   (sd_asl(m,n),n=1,sd_numasl),       &
               &                   sd_n(m), cstat, m
       end do

    else
!       write(fid_sdm_o,'(3a,i2.2,6a)') '# x[m],y[m],z[m],vz[m],',       &
       write(fid_sdm_o,'(3a,i2.2,5a)') '# x[m],y[m],z[m],vz[m],',       &
            &            'radius(droplet)[m],',                            &
            &            'mass_of_aerosol_in_droplet/ice(1:',sd_numasl,')[g],',&
            &            'radius_eq(ice)[m],radius_pol(ice)[m],',              &
            &            'density(droplet/ice)[kg/m3],',                       &
!            &            'temperature[K],',                                    &
            &            'freezing_temp.(ice)[deg],',                          &
            &            'multiplicity[-],status[-],index,rime_mass[kg],num_of_monomers[-]'

!       write(fmt,'( "(", i2.2, "e16.8,i20,a5,i10)" )')(sd_numasl+10)
       write(fmt,'( "(", i2.2, "e16.8,i20,a5,i10,e16.8,i10)" )')(sd_numasl+9)

       do m=1,sd_num,sdn_dmpnskip

          if( sd_liqice(m)==STAT_LIQ ) then
             cstat = '  LIQ'
          else if( sd_liqice(m)==STAT_ICE ) then
             cstat = '  ICE'
          else if( sd_liqice(m)==STAT_MIX ) then
             cstat = ' MELT'
          end if

          write(fid_sdm_o,trim(fmt)) sd_x(m),           &
               &                   sd_y(m),           &
               &                   sd_z(m), sd_vz(m), sd_r(m),        &
               &                   (sd_asl(m,n),n=1,sd_numasl),       &
               &                   sdi%re(m), sdi%rp(m), sdi%rho(m),  &
!               &                   sdi%t(m), sdi%tf(m),               &
               &                   sdi%tf(m),               &
               &                   sd_n(m), cstat, m, sdi%mrime(m), sdi%nmono(m)
       end do


    end if

    close(fid_sdm_o)
    if( IO_L ) write(IO_FID_LOG,*) '*** Closed output file (ASCII) of Super Droplet'

    return

  end subroutine sdm_outasci
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sdm_assign_tracking_subset(sd_num, sd_z, sd_r, sd_id, dm_id, if_coal, sd_liqice, sdi)
    use scale_precision
    use rng_uniform_mt, only: &
         rng_init, rng_generate
    use scale_stdio, only: &
         IO_L, IO_FID_LOG
    use scale_process, only: &
         mype => PRC_myrank
    use m_sdm_common, only: &
         i2, INVALID_i4, tracking_mode, tracking_selection_mode, tracking_fraction, max_tracked_sds, &
         tracking_height_min, tracking_height_max, tracking_radius_min, tracking_radius_max, tracking_nz_bin, tracking_nr_bin, &
         tracking_min_per_bin, tracking_fallback_to_random, tracking_sample_initialized, tracking_sampling_seed, rng_tracking_s2c, &
         sdicedef
    use m_sdm_idutil, only: &
         sdm_tracking_effective_radius
    use m_sdm_tracking_cold, only: &
         sdm_tracking_valid_id_pair

    implicit none

    integer, intent(in) :: sd_num
    real(RP), intent(in) :: sd_z(1:sd_num)
    real(RP), intent(in) :: sd_r(1:sd_num)
    integer, intent(inout) :: sd_id(1:sd_num)
    integer, intent(inout) :: dm_id(1:sd_num)
    integer(kind=i2), intent(inout) :: if_coal(1:sd_num)
    integer(kind=i2), intent(in), optional :: sd_liqice(1:sd_num)
    type(sdicedef), intent(in), optional :: sdi

    integer :: n, ibin, iz, ir
    integer :: tracked_cnt, candidate_cnt, target_cnt
    integer :: nbin, sum_quota, extra_needed, reduce_needed, idx_best
    integer :: needed_in_bin, non_empty_bins, min_per_bin_eff
    integer :: upper_exceed_cnt, near_upper_cnt
    real(RP) :: rand_tracking
    real(RP) :: z_span, r_min_eff, r_max_eff, log_r_span
    real(RP) :: best_frac, max_radius_in_height, near_upper_threshold
    logical :: do_track, use_stratified, stratified_selected, use_radius_upper_bound
    integer, allocatable :: bin_cnt(:), bin_quota(:), bin_selected(:), bin_remaining(:), bin_min(:)
    real(RP), allocatable :: bin_frac(:)
    real(RP), allocatable :: sd_track_radius(:)

    allocate(sd_track_radius(sd_num))
    do n=1,sd_num
      if( present(sd_liqice) .and. present(sdi) ) then
        sd_track_radius(n) = sdm_tracking_effective_radius(sd_r(n), sd_liqice(n), sdi%re(n), sdi%rp(n))
      else
        sd_track_radius(n) = sd_r(n)
      end if
    end do

    if( tracking_mode /= 2 .or. tracking_fraction <= 0.0_RP ) then
      do n=1,sd_num
        sd_id(n) = INVALID_i4
        dm_id(n) = INVALID_i4
        if_coal(n)  = 0_i2
      end do
      tracking_sample_initialized = .false.
      return
    end if

    tracked_cnt = 0
    if( .not. tracking_sample_initialized ) then
      call rng_init( rng_tracking_s2c, mype + tracking_sampling_seed )
      if( tracking_fraction >= 1.0_RP .and. max_tracked_sds <= 0 ) then
        do n=1,sd_num
          if( sdm_tracking_valid_id_pair(sd_id(n), dm_id(n)) ) then
            sd_id(n) = n
            dm_id(n) = mype
            if_coal(n) = 0_i2
          else
            sd_id(n) = INVALID_i4
            dm_id(n) = INVALID_i4
            if_coal(n) = 0_i2
          end if
        end do
        tracking_sample_initialized = .true.
        return
      end if

      use_stratified = trim(adjustl(tracking_selection_mode)) == 'stratified' .or. &
           trim(adjustl(tracking_selection_mode)) == 'STRATIFIED' .or. &
           trim(adjustl(tracking_selection_mode)) == 'Stratified'
      use_radius_upper_bound = tracking_radius_max > tracking_radius_min
      stratified_selected = .false.

      if( use_stratified ) then
        if( tracking_nz_bin > 0 .and. tracking_nr_bin > 0 .and. tracking_height_max > tracking_height_min ) then
          nbin = tracking_nz_bin * tracking_nr_bin
          allocate(bin_cnt(nbin), bin_quota(nbin), bin_selected(nbin), bin_remaining(nbin), bin_min(nbin), bin_frac(nbin))

          candidate_cnt = 0
          r_min_eff = max(tracking_radius_min, 1.0E-12_RP)
          r_max_eff = r_min_eff
          upper_exceed_cnt = 0
          near_upper_cnt = 0
          max_radius_in_height = 0.0_RP
          near_upper_threshold = tracking_radius_max * 0.98_RP
          do n=1,sd_num
            if( sd_z(n) >= tracking_height_min .and. sd_z(n) <= tracking_height_max ) then
              if( sd_track_radius(n) > max_radius_in_height ) max_radius_in_height = sd_track_radius(n)
              if( use_radius_upper_bound ) then
                if( sd_track_radius(n) > tracking_radius_max ) then
                  upper_exceed_cnt = upper_exceed_cnt + 1
                else if( sd_track_radius(n) >= near_upper_threshold ) then
                  near_upper_cnt = near_upper_cnt + 1
                end if
              end if
            end if
            if( sd_z(n) >= tracking_height_min .and. sd_z(n) <= tracking_height_max .and. &
                sd_track_radius(n) >= tracking_radius_min .and. &
                ( .not. use_radius_upper_bound .or. sd_track_radius(n) <= tracking_radius_max ) ) then
              candidate_cnt = candidate_cnt + 1
              if( sd_track_radius(n) > r_max_eff ) r_max_eff = sd_track_radius(n)
            end if
          end do
          if( use_radius_upper_bound .and. IO_L ) then
            if( upper_exceed_cnt > 0 ) then
              write(IO_FID_LOG,*) '*** WARNING (sdm_assign_tracking_subset): stratified candidates exceed tracking_radius_max.'
              write(IO_FID_LOG,*) '    upper_exceed_cnt=', upper_exceed_cnt, ' tracking_radius_max[m]=', tracking_radius_max, &
                   ' max_radius_in_height[m]=', max_radius_in_height
            else if( near_upper_cnt > 0 ) then
              write(IO_FID_LOG,*) '*** WARNING (sdm_assign_tracking_subset): candidate radii are close to tracking_radius_max.'
              write(IO_FID_LOG,*) '    near_upper_cnt=', near_upper_cnt, ' near_upper_threshold[m]=', near_upper_threshold, &
                   ' tracking_radius_max[m]=', tracking_radius_max
            end if
          end if

          target_cnt = int( real(candidate_cnt,kind=RP) * tracking_fraction + 0.5_RP )
          if( candidate_cnt > 0 .and. tracking_fraction > 0.0_RP .and. target_cnt == 0 ) target_cnt = 1
          if( target_cnt > candidate_cnt ) target_cnt = candidate_cnt
          if( max_tracked_sds > 0 ) target_cnt = min(target_cnt, max_tracked_sds)

          if( candidate_cnt > 0 .and. target_cnt > 0 ) then
            z_span = tracking_height_max - tracking_height_min
            if( r_max_eff > r_min_eff ) then
              log_r_span = log(r_max_eff / r_min_eff)
            else
              log_r_span = 0.0_RP
            end if

            bin_cnt(:) = 0
            bin_quota(:) = 0
            bin_selected(:) = 0
            bin_remaining(:) = 0
            bin_min(:) = 0
            bin_frac(:) = 0.0_RP

            do n=1,sd_num
              if( sd_z(n) >= tracking_height_min .and. sd_z(n) <= tracking_height_max .and. &
                  sd_track_radius(n) >= tracking_radius_min .and. &
                  ( .not. use_radius_upper_bound .or. sd_track_radius(n) <= tracking_radius_max ) ) then
                iz = int( (sd_z(n)-tracking_height_min) / z_span * real(tracking_nz_bin,kind=RP) ) + 1
                iz = min( tracking_nz_bin, max(1,iz) )
                if( tracking_nr_bin == 1 ) then
                  ir = 1
                else
                  if( log_r_span > 0.0_RP .and. sd_track_radius(n) > r_min_eff ) then
                    ir = int( log(sd_track_radius(n)/r_min_eff) / log_r_span * real(tracking_nr_bin,kind=RP) ) + 1
                  else
                    ir = 1
                  end if
                  ir = min( tracking_nr_bin, max(1,ir) )
                end if
                ibin = (iz-1)*tracking_nr_bin + ir
                bin_cnt(ibin) = bin_cnt(ibin) + 1
              end if
            end do

            non_empty_bins = count(bin_cnt > 0)
            if( non_empty_bins > 0 ) then
              min_per_bin_eff = max(0, tracking_min_per_bin)
              min_per_bin_eff = min(min_per_bin_eff, target_cnt / non_empty_bins)
            else
              min_per_bin_eff = 0
            end if

            do ibin=1,nbin
              if( bin_cnt(ibin) > 0 ) then
                bin_quota(ibin) = int( real(target_cnt,kind=RP) * real(bin_cnt(ibin),kind=RP) / real(candidate_cnt,kind=RP) )
                bin_quota(ibin) = min(bin_quota(ibin), bin_cnt(ibin))
                if( min_per_bin_eff > 0 ) then
                  bin_min(ibin) = min(min_per_bin_eff, bin_cnt(ibin))
                  if( bin_quota(ibin) < bin_min(ibin) ) bin_quota(ibin) = bin_min(ibin)
                end if
              end if
            end do

            sum_quota = sum(bin_quota)
            if( sum_quota < target_cnt ) then
              extra_needed = target_cnt - sum_quota
              do while( extra_needed > 0 )
                idx_best = 0
                best_frac = -1.0_RP
                do ibin=1,nbin
                  if( bin_quota(ibin) < bin_cnt(ibin) ) then
                    if( bin_cnt(ibin) > 0 ) then
                      bin_frac(ibin) = real(bin_quota(ibin),kind=RP) / real(bin_cnt(ibin),kind=RP)
                    else
                      bin_frac(ibin) = 1.0_RP
                    end if
                    if( idx_best == 0 .or. bin_frac(ibin) < best_frac ) then
                      idx_best = ibin
                      best_frac = bin_frac(ibin)
                    end if
                  end if
                end do
                if( idx_best == 0 ) exit
                bin_quota(idx_best) = bin_quota(idx_best) + 1
                extra_needed = extra_needed - 1
              end do
            else if( sum_quota > target_cnt ) then
              reduce_needed = sum_quota - target_cnt
              do while( reduce_needed > 0 )
                idx_best = 0
                best_frac = 2.0_RP
                do ibin=1,nbin
                  if( bin_quota(ibin) > bin_min(ibin) ) then
                    if( bin_cnt(ibin) > 0 ) then
                      bin_frac(ibin) = real(bin_quota(ibin),kind=RP) / real(bin_cnt(ibin),kind=RP)
                    else
                      bin_frac(ibin) = 0.0_RP
                    end if
                    if( idx_best == 0 .or. bin_frac(ibin) > best_frac ) then
                      idx_best = ibin
                      best_frac = bin_frac(ibin)
                    end if
                  end if
                end do
                if( idx_best == 0 ) exit
                bin_quota(idx_best) = bin_quota(idx_best) - 1
                reduce_needed = reduce_needed - 1
              end do
            end if

            bin_remaining(:) = bin_cnt(:)
            do n=1,sd_num
              do_track = .false.
              if( sd_z(n) >= tracking_height_min .and. sd_z(n) <= tracking_height_max .and. &
                  sd_track_radius(n) >= tracking_radius_min .and. &
                  ( .not. use_radius_upper_bound .or. sd_track_radius(n) <= tracking_radius_max ) ) then
                iz = int( (sd_z(n)-tracking_height_min) / z_span * real(tracking_nz_bin,kind=RP) ) + 1
                iz = min( tracking_nz_bin, max(1,iz) )
                if( tracking_nr_bin == 1 ) then
                  ir = 1
                else
                  if( log_r_span > 0.0_RP .and. sd_track_radius(n) > r_min_eff ) then
                    ir = int( log(sd_track_radius(n)/r_min_eff) / log_r_span * real(tracking_nr_bin,kind=RP) ) + 1
                  else
                    ir = 1
                  end if
                  ir = min( tracking_nr_bin, max(1,ir) )
                end if
                ibin = (iz-1)*tracking_nr_bin + ir
                needed_in_bin = bin_quota(ibin) - bin_selected(ibin)
                if( needed_in_bin > 0 .and. bin_remaining(ibin) > 0 ) then
                  rand_tracking = real(rng_generate(rng_tracking_s2c), kind=RP)
                  if( rand_tracking <= real(needed_in_bin,kind=RP) / real(bin_remaining(ibin),kind=RP) ) then
                    do_track = .true.
                    bin_selected(ibin) = bin_selected(ibin) + 1
                  end if
                end if
                if( bin_remaining(ibin) > 0 ) bin_remaining(ibin) = bin_remaining(ibin) - 1
              end if
              if( do_track ) then
                tracked_cnt = tracked_cnt + 1
                sd_id(n) = n
                dm_id(n) = mype
              else
                sd_id(n) = INVALID_i4
                dm_id(n) = INVALID_i4
              end if
              if_coal(n) = 0_i2
            end do
            tracking_sample_initialized = .true.
            stratified_selected = .true.
          end if

          deallocate(bin_cnt, bin_quota, bin_selected, bin_remaining, bin_min, bin_frac)
        end if
      end if

      if( .not. stratified_selected ) then
        if( use_stratified .and. .not. tracking_fallback_to_random ) then
          do n=1,sd_num
            sd_id(n) = INVALID_i4
            dm_id(n) = INVALID_i4
            if_coal(n) = 0_i2
          end do
          tracking_sample_initialized = .true.
        else
          do n=1,sd_num
            do_track = .true.
            if( tracking_fraction < 1.0_RP ) then
              rand_tracking = real(rng_generate(rng_tracking_s2c), kind=RP)
              if( rand_tracking > tracking_fraction ) do_track = .false.
            end if
            if( max_tracked_sds > 0 ) then
              if( tracked_cnt >= max_tracked_sds ) do_track = .false.
            end if
            if( do_track ) then
              tracked_cnt = tracked_cnt + 1
              sd_id(n) = n
              dm_id(n) = mype
            else
              sd_id(n) = INVALID_i4
              dm_id(n) = INVALID_i4
            end if
            if_coal(n) = 0_i2
          end do
          tracking_sample_initialized = .true.
        end if
      end if
    else
      do n=1,sd_num
        do_track = sdm_tracking_valid_id_pair(sd_id(n), dm_id(n))
        if( do_track ) then
          tracked_cnt = tracked_cnt + 1
          sd_id(n) = n
          dm_id(n) = mype
        else
          sd_id(n) = INVALID_i4
          dm_id(n) = INVALID_i4
        end if
        if_coal(n) = 0_i2
      end do
    end if

    return
  end subroutine sdm_assign_tracking_subset
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sdm_outnetcdf(otime,sd_num,sd_numasl,sd_n,sd_liqice,sd_x,sd_y,sd_z,sd_r,sd_asl,sd_vz,sdi,sd_id,dm_id,if_coal, &
       &                   sdn_dmpnskip,filetag,sd_event_mask,sd_event_sig_mask,sd_diag_mask,sd_phase_change_flag,                       &
       &                   sd_spatial_visit_flag,                                                                       &
       &                   sd_liq_radius_max_interval,sd_ice_rvol_max_interval,sd_mixed_rvol_max_interval,              &
       &                   sd_rime_mass_max_interval,sd_rime_frac_max_interval,sd_nmono_max_interval,                   &
       &                   sd_aspect_ratio_max_interval)
    use netcdf
    use scale_precision
    use scale_stdio
    use scale_time
    use scale_process, only: &
         mype => PRC_myrank
    use m_sdm_common, only: &
         i2, sdm_cold, STAT_LIQ, STAT_ICE, STAT_MIX, sdicedef, &
         INVALID_i4, tracking_mode, forward_tracking_enable, backward_tracking_enable, coalescence_output_enable

    implicit none

    real(DP), intent(in) :: otime
    integer, intent(in) :: sd_num    ! number of super-droplets
    integer, intent(in) :: sd_numasl ! number of chemical species contained in super droplets
    integer(DP), intent(in) :: sd_n(1:sd_num) ! multiplicity of super-droplets
    integer(kind=i2), intent(inout) :: if_coal(1:sd_num)
                       ! flag of coalescence
                       ! 0 = Super Droplet hasn't undergone coalescence during the previous output interval
                       ! 1 = Super Droplet has undergone coalescence during the previous output interval
    integer(kind=i2), intent(in) :: sd_liqice(1:sd_num)
                       ! status of super-droplets (liquid/ice)
                       ! 01 = all liquid, 10 = all ice
                       ! 11 = mixture of ice and liquid
    real(RP), intent(in) :: sd_x(1:sd_num) ! x-coordinate of super-droplets
    real(RP), intent(in) :: sd_y(1:sd_num) ! y-coordinate of super-droplets
    real(RP), intent(in) :: sd_z(1:sd_num) ! z-coordinate of super-droplets
    real(RP), intent(in) :: sd_r(1:sd_num) ! equivalent radius of super-droplets
    real(RP), intent(in) :: sd_asl(1:sd_num,1:sd_numasl) ! aerosol mass of super-droplets
    real(RP), intent(in) :: sd_vz(1:sd_num) ! terminal velocity of super-droplets
    type(sdicedef), intent(in) :: sdi   ! ice phase super-droplets
    integer, intent(in) :: sd_id(1:sd_num) ! save index of super-droplets
    integer, intent(in) :: dm_id(1:sd_num) ! domain index of super-droplets
    integer, intent(in) :: sdn_dmpnskip ! Base skip to store super droplets in text format
    character(len=*),intent(in),optional :: filetag ! user defined text string tag to be added to the filenames
    integer, intent(in), optional :: sd_event_mask(1:sd_num)
    integer, intent(in), optional :: sd_event_sig_mask(1:sd_num)
    integer, intent(in), optional :: sd_diag_mask(1:sd_num)
    integer, intent(in), optional :: sd_phase_change_flag(1:sd_num)
    integer, intent(in), optional :: sd_spatial_visit_flag(1:sd_num)
    real(RP), intent(in), optional :: sd_liq_radius_max_interval(1:sd_num)
    real(RP), intent(in), optional :: sd_ice_rvol_max_interval(1:sd_num)
    real(RP), intent(in), optional :: sd_mixed_rvol_max_interval(1:sd_num)
    real(RP), intent(in), optional :: sd_rime_mass_max_interval(1:sd_num)
    real(RP), intent(in), optional :: sd_rime_frac_max_interval(1:sd_num)
    real(RP), intent(in), optional :: sd_nmono_max_interval(1:sd_num)
    real(RP), intent(in), optional :: sd_aspect_ratio_max_interval(1:sd_num)

    character(len=17) :: fmt2="(A, '.', A, I*.*)"
    character(len=H_LONG) :: ftmp
    character(len=H_LONG) :: basename_sd_out
    character(len=19) :: basename_time
    integer :: fid_sdm_o
    integer :: n, m, ierr
    character(len=80) :: fmt     ! output formate
    character(len=5)  :: cstat   ! status character
    integer :: nf90_real_precision
    integer :: ncid, sd_num_id, sd_numasl_id
    integer :: sd_x_id, sd_y_id, sd_z_id, sd_vz_id, sd_r_id, sd_asl_id, sd_n_id, sd_liqice_id
    integer :: sd_id_id, domain_id, if_coal_id
    integer :: sd_event_mask_id, sd_event_sig_mask_id, sd_diag_mask_id, sd_phase_change_flag_id, sd_spatial_visit_flag_id
    integer :: sd_liq_radius_max_interval_id, sd_ice_rvol_max_interval_id, sd_mixed_rvol_max_interval_id
    integer :: sd_rime_mass_max_interval_id, sd_rime_frac_max_interval_id
    integer :: sd_nmono_max_interval_id, sd_aspect_ratio_max_interval_id
    integer :: sdi_re_id, sdi_rp_id, sdi_rho_id, sdi_tf_id, sdi_mrime_id, sdi_nmono_id
    logical :: write_forward_tracking, write_backward_tracking, write_tracking, write_coal, write_cold_tracking
    character(len=80) :: tracking_id_label

    integer,parameter :: nc_deflate_level = 1      ! NetCDF compression level {1,..,9}
    integer,parameter :: nc_deflate       = 1 ! turn on NetCDF compresion
    integer,parameter :: nc_shuffle       = 1 ! turn on NetCDF shuffle filter
    character(len=100) :: ftag ! =filetag or ''(default)

    !--- output Super Droplets in NetCDF format
    call TIME_gettimelabel(basename_time)
    fid_sdm_o = IO_get_available_fid()

    if(present(filetag)) then
       write(ftag,'(2A)') 'SD_', trim(filetag)
    else
       ftag = 'SD_output'
    end if

    write(fmt2(14:14),'(I1)') 6
    write(fmt2(16:16),'(I1)') 6
    write(ftmp,'(3A)') trim(ftag), '_NetCDF_', trim(basename_time)
    write(basename_sd_out,fmt2) trim(ftmp), 'pe',mype

    ! Check the presision
    if(RP == SP)then
       nf90_real_precision = NF90_FLOAT
    else
       nf90_real_precision = NF90_DOUBLE
    end if

    ! Open the output file
    call check_netcdf( nf90_create(trim(basename_sd_out),NF90_NETCDF4, ncid) )

    ! Definition of dimensions and variables
    !!! sd_num
    call check_netcdf( nf90_def_dim(ncid, "sd_num", sd_num, sd_num_id) )
    !!! sd_numasl
    call check_netcdf( nf90_def_dim(ncid, "sd_numasl", sd_numasl, sd_numasl_id) )

    !!! sd_x
    call check_netcdf( nf90_def_var(ncid, "sd_x", nf90_real_precision, sd_num_id, sd_x_id) )
    call check_netcdf( nf90_def_var_deflate(ncid, sd_x_id, &
         & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
    call check_netcdf( nf90_put_att(ncid, sd_x_id, 'long_name', 'x-coordinate') )
    call check_netcdf( nf90_put_att(ncid, sd_x_id, 'units', 'm') )
    !!! sd_y
    call check_netcdf( nf90_def_var(ncid, "sd_y", nf90_real_precision, sd_num_id, sd_y_id) )
    call check_netcdf( nf90_def_var_deflate(ncid, sd_y_id, &
         & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
    call check_netcdf( nf90_put_att(ncid, sd_y_id, 'long_name', 'y-coordinate') )
    call check_netcdf( nf90_put_att(ncid, sd_y_id, 'units', 'm') )
    !!! sd_z
    call check_netcdf( nf90_def_var(ncid, "sd_z", nf90_real_precision, sd_num_id, sd_z_id) )
    call check_netcdf( nf90_def_var_deflate(ncid, sd_z_id, &
         & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
    call check_netcdf( nf90_put_att(ncid, sd_z_id, 'long_name', 'z-coordinate') )
    call check_netcdf( nf90_put_att(ncid, sd_z_id, 'units', 'm') )
    !!! sd_vz
    call check_netcdf( nf90_def_var(ncid, "sd_vz", nf90_real_precision, sd_num_id, sd_vz_id) )
    call check_netcdf( nf90_def_var_deflate(ncid, sd_vz_id, &
         & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
    call check_netcdf( nf90_put_att(ncid, sd_vz_id, 'long_name', 'terminal velocity') )
    call check_netcdf( nf90_put_att(ncid, sd_vz_id, 'units', 'm/s') )
    !!! sd_r
    call check_netcdf( nf90_def_var(ncid, "sd_r", nf90_real_precision, sd_num_id, sd_r_id) )
    call check_netcdf( nf90_def_var_deflate(ncid, sd_r_id, &
         & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
    call check_netcdf( nf90_put_att(ncid, sd_r_id, 'long_name', 'equivalent radius of liquid droplets') )
    call check_netcdf( nf90_put_att(ncid, sd_r_id, 'units', 'm') )
    !!! sd_asl
    call check_netcdf( nf90_def_var(ncid, "sd_asl", nf90_real_precision, (/sd_num_id, sd_numasl_id/), sd_asl_id) )
    call check_netcdf( nf90_def_var_deflate(ncid, sd_asl_id, &
         & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
    call check_netcdf( nf90_put_att(ncid, sd_asl_id, 'long_name', 'aerosol mass') )
    call check_netcdf( nf90_put_att(ncid, sd_asl_id, 'units', 'g') )

    !!! sd_n
    call check_netcdf( nf90_def_var(ncid, "sd_n", NF90_INT64, sd_num_id, sd_n_id) )
    call check_netcdf( nf90_def_var_deflate(ncid, sd_n_id, &
         & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
    call check_netcdf( nf90_put_att(ncid, sd_n_id, 'long_name', 'multiplicity') )
    call check_netcdf( nf90_put_att(ncid, sd_n_id, 'units', '') )
    !!! sd_liqice
    call check_netcdf( nf90_def_var(ncid, "sd_liqice", NF90_SHORT, sd_num_id, sd_liqice_id) )
    call check_netcdf( nf90_def_var_deflate(ncid, sd_liqice_id, &
         & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
    call check_netcdf( nf90_put_att(ncid, sd_liqice_id, 'long_name', 'status of droplets: 01=liquid, 10=ice, 11=mixture') )
    call check_netcdf( nf90_put_att(ncid, sd_liqice_id, 'units', '') )

    write_forward_tracking = tracking_mode == 1
    write_backward_tracking = tracking_mode == 2
    write_tracking = write_forward_tracking .or. write_backward_tracking
    write_coal = coalescence_output_enable .and. (.not. sdm_cold)
    write_cold_tracking = sdm_cold .and. present(sd_event_mask) .and. present(sd_event_sig_mask) .and. &
         &                 present(sd_diag_mask) .and. &
         &                 present(sd_phase_change_flag) .and. present(sd_spatial_visit_flag) .and. &
         &                 present(sd_liq_radius_max_interval) .and. &
         &                 present(sd_ice_rvol_max_interval) .and. present(sd_mixed_rvol_max_interval) .and. &
         &                 present(sd_rime_mass_max_interval) .and. present(sd_rime_frac_max_interval) .and. &
         &                 present(sd_nmono_max_interval) .and. present(sd_aspect_ratio_max_interval)
    if( write_backward_tracking ) then
      tracking_id_label = 'pre_sdid+pre_dmid'
      call check_netcdf( nf90_def_var(ncid, "pre_sdid", NF90_INT, sd_num_id, sd_id_id) )
      call check_netcdf( nf90_def_var_deflate(ncid, sd_id_id, &
           & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
      call check_netcdf( nf90_put_att(ncid, sd_id_id, 'long_name', trim(tracking_id_label)) )
      call check_netcdf( nf90_put_att(ncid, sd_id_id, 'units', '') )
      call check_netcdf( nf90_def_var(ncid, "pre_dmid", NF90_INT, sd_num_id, domain_id) )
      call check_netcdf( nf90_def_var_deflate(ncid, domain_id, &
           & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
      call check_netcdf( nf90_put_att(ncid, domain_id, 'long_name', trim(tracking_id_label)) )
      call check_netcdf( nf90_put_att(ncid, domain_id, 'units', '') )
    else if( write_forward_tracking ) then
      call check_netcdf( nf90_def_var(ncid, "sd_id", NF90_INT, sd_num_id, sd_id_id) )
      call check_netcdf( nf90_def_var_deflate(ncid, sd_id_id, &
           & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
      call check_netcdf( nf90_put_att(ncid, sd_id_id, 'long_name', 'SD ID') )
      call check_netcdf( nf90_put_att(ncid, sd_id_id, 'units', '') )
      call check_netcdf( nf90_def_var(ncid, "dm_id", NF90_INT, sd_num_id, domain_id) )
      call check_netcdf( nf90_def_var_deflate(ncid, domain_id, &
           & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
      call check_netcdf( nf90_put_att(ncid, domain_id, 'long_name', 'domain ID') )
      call check_netcdf( nf90_put_att(ncid, domain_id, 'units', '') )
    end if
    if( write_coal ) then
      call check_netcdf( nf90_def_var(ncid, "if_coal", NF90_SHORT, sd_num_id, if_coal_id) )
      call check_netcdf( nf90_def_var_deflate(ncid, if_coal_id, &
           & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
      call check_netcdf( nf90_put_att(ncid, if_coal_id, 'long_name', 'coalescence flag: 0=has not occurred, 1=occurred') )
      call check_netcdf( nf90_put_att(ncid, if_coal_id, 'units', '') )
    end if
    if( write_cold_tracking ) then
      call check_netcdf( nf90_def_var(ncid, "sd_event_mask", NF90_INT, sd_num_id, sd_event_mask_id) )
      call check_netcdf( nf90_def_var_deflate(ncid, sd_event_mask_id, &
           & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
      call check_netcdf( nf90_put_att(ncid, sd_event_mask_id, 'long_name', 'cold SD event bitmask over output interval') )
      call check_netcdf( nf90_put_att(ncid, sd_event_mask_id, 'units', '') )
      call check_netcdf( nf90_def_var(ncid, "sd_event_sig_mask", NF90_INT, sd_num_id, sd_event_sig_mask_id) )
      call check_netcdf( nf90_def_var_deflate(ncid, sd_event_sig_mask_id, &
           & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
      call check_netcdf( nf90_put_att(ncid, sd_event_sig_mask_id, 'long_name', &
           'cold SD significant event subset bitmask over output interval') )
      call check_netcdf( nf90_put_att(ncid, sd_event_sig_mask_id, 'units', '') )
      call check_netcdf( nf90_def_var(ncid, "sd_diag_mask", NF90_INT, sd_num_id, sd_diag_mask_id) )
      call check_netcdf( nf90_def_var_deflate(ncid, sd_diag_mask_id, &
           & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
      call check_netcdf( nf90_put_att(ncid, sd_diag_mask_id, 'long_name', 'cold SD diagnostic bitmask evaluated at output time') )
      call check_netcdf( nf90_put_att(ncid, sd_diag_mask_id, 'units', '') )
      call check_netcdf( nf90_def_var(ncid, "sd_phase_change_flag", NF90_INT, sd_num_id, sd_phase_change_flag_id) )
      call check_netcdf( nf90_def_var_deflate(ncid, sd_phase_change_flag_id, &
           & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
      call check_netcdf( nf90_put_att(ncid, sd_phase_change_flag_id, 'long_name', 'cold phase-state category-change flag over output interval') )
      call check_netcdf( nf90_put_att(ncid, sd_phase_change_flag_id, 'units', '') )
      call check_netcdf( nf90_def_var(ncid, "sd_spatial_visit_flag", NF90_INT, sd_num_id, sd_spatial_visit_flag_id) )
      call check_netcdf( nf90_def_var_deflate(ncid, sd_spatial_visit_flag_id, &
           & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
      call check_netcdf( nf90_put_att(ncid, sd_spatial_visit_flag_id, 'long_name', &
           'cold spatial-visit flag over output interval') )
      call check_netcdf( nf90_put_att(ncid, sd_spatial_visit_flag_id, 'units', '') )
      call check_netcdf( nf90_def_var(ncid, "sd_liq_radius_max_interval", nf90_real_precision, sd_num_id, sd_liq_radius_max_interval_id) )
      call check_netcdf( nf90_def_var_deflate(ncid, sd_liq_radius_max_interval_id, &
           & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
      call check_netcdf( nf90_put_att(ncid, sd_liq_radius_max_interval_id, 'long_name', 'maximum liquid radius over output interval') )
      call check_netcdf( nf90_put_att(ncid, sd_liq_radius_max_interval_id, 'units', 'm') )
      call check_netcdf( nf90_def_var(ncid, "sd_ice_rvol_max_interval", nf90_real_precision, sd_num_id, sd_ice_rvol_max_interval_id) )
      call check_netcdf( nf90_def_var_deflate(ncid, sd_ice_rvol_max_interval_id, &
           & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
      call check_netcdf( nf90_put_att(ncid, sd_ice_rvol_max_interval_id, 'long_name', 'maximum ice volume-equivalent radius over output interval') )
      call check_netcdf( nf90_put_att(ncid, sd_ice_rvol_max_interval_id, 'units', 'm') )
      call check_netcdf( nf90_def_var(ncid, "sd_mixed_rvol_max_interval", nf90_real_precision, sd_num_id, sd_mixed_rvol_max_interval_id) )
      call check_netcdf( nf90_def_var_deflate(ncid, sd_mixed_rvol_max_interval_id, &
           & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
      call check_netcdf( nf90_put_att(ncid, sd_mixed_rvol_max_interval_id, 'long_name', 'maximum mixed-phase volume-equivalent radius over output interval') )
      call check_netcdf( nf90_put_att(ncid, sd_mixed_rvol_max_interval_id, 'units', 'm') )
      call check_netcdf( nf90_def_var(ncid, "sd_rime_mass_max_interval", nf90_real_precision, sd_num_id, sd_rime_mass_max_interval_id) )
      call check_netcdf( nf90_def_var_deflate(ncid, sd_rime_mass_max_interval_id, &
           & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
      call check_netcdf( nf90_put_att(ncid, sd_rime_mass_max_interval_id, 'long_name', 'maximum rime mass over output interval') )
      call check_netcdf( nf90_put_att(ncid, sd_rime_mass_max_interval_id, 'units', 'kg') )
      call check_netcdf( nf90_def_var(ncid, "sd_rime_frac_max_interval", nf90_real_precision, sd_num_id, sd_rime_frac_max_interval_id) )
      call check_netcdf( nf90_def_var_deflate(ncid, sd_rime_frac_max_interval_id, &
           & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
      call check_netcdf( nf90_put_att(ncid, sd_rime_frac_max_interval_id, 'long_name', 'maximum rime mass fraction over output interval') )
      call check_netcdf( nf90_put_att(ncid, sd_rime_frac_max_interval_id, 'units', '') )
      call check_netcdf( nf90_def_var(ncid, "sd_nmono_max_interval", nf90_real_precision, sd_num_id, sd_nmono_max_interval_id) )
      call check_netcdf( nf90_def_var_deflate(ncid, sd_nmono_max_interval_id, &
           & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
      call check_netcdf( nf90_put_att(ncid, sd_nmono_max_interval_id, 'long_name', 'maximum monomer count over output interval') )
      call check_netcdf( nf90_put_att(ncid, sd_nmono_max_interval_id, 'units', '') )
      call check_netcdf( nf90_def_var(ncid, "sd_aspect_ratio_max_interval", nf90_real_precision, sd_num_id, sd_aspect_ratio_max_interval_id) )
      call check_netcdf( nf90_def_var_deflate(ncid, sd_aspect_ratio_max_interval_id, &
           & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
      call check_netcdf( nf90_put_att(ncid, sd_aspect_ratio_max_interval_id, 'long_name', 'maximum ice aspect ratio over output interval') )
      call check_netcdf( nf90_put_att(ncid, sd_aspect_ratio_max_interval_id, 'units', '') )
    end if

    if( sdm_cold ) then
       !!! sdi%re
       call check_netcdf( nf90_def_var(ncid, "sdi_re", nf90_real_precision, sd_num_id, sdi_re_id) )
       call check_netcdf( nf90_def_var_deflate(ncid, sdi_re_id, &
            & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
       call check_netcdf( nf90_put_att(ncid, sdi_re_id, 'long_name', 'equatorial radius of ice crystals') )
       call check_netcdf( nf90_put_att(ncid, sdi_re_id, 'units', 'm') )
       !!! sdi%rp
       call check_netcdf( nf90_def_var(ncid, "sdi_rp", nf90_real_precision, sd_num_id, sdi_rp_id) )
       call check_netcdf( nf90_def_var_deflate(ncid, sdi_rp_id, &
            & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
       call check_netcdf( nf90_put_att(ncid, sdi_rp_id, 'long_name', 'polar radius of ice crystals') )
       call check_netcdf( nf90_put_att(ncid, sdi_rp_id, 'units', 'm') )
       !!! sdi%rho
       call check_netcdf( nf90_def_var(ncid, "sdi_rho", nf90_real_precision, sd_num_id, sdi_rho_id) )
       call check_netcdf( nf90_def_var_deflate(ncid, sdi_rho_id, &
            & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
       call check_netcdf( nf90_put_att(ncid, sdi_rho_id, 'long_name', 'density of ice crystals') )
       call check_netcdf( nf90_put_att(ncid, sdi_rho_id, 'units', 'kg/m3') )
       !!! sdi%tf
       call check_netcdf( nf90_def_var(ncid, "sdi_tf", nf90_real_precision, sd_num_id, sdi_tf_id) )
       call check_netcdf( nf90_def_var_deflate(ncid, sdi_tf_id, &
            & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
       call check_netcdf( nf90_put_att(ncid, sdi_tf_id, 'long_name', 'freezing temperature of particles') )
       call check_netcdf( nf90_put_att(ncid, sdi_tf_id, 'units', 'degC') )
       !!! sdi%mrime
       call check_netcdf( nf90_def_var(ncid, "sdi_mrime", nf90_real_precision, sd_num_id, sdi_mrime_id) )
       call check_netcdf( nf90_def_var_deflate(ncid, sdi_mrime_id, &
            & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
       call check_netcdf( nf90_put_att(ncid, sdi_mrime_id, 'long_name', 'rime mass') )
       call check_netcdf( nf90_put_att(ncid, sdi_mrime_id, 'units', 'kg') )
       !!! sdi%nmono
       call check_netcdf( nf90_def_var(ncid, "sdi_nmono", NF90_INT, sd_num_id, sdi_nmono_id) )
       call check_netcdf( nf90_def_var_deflate(ncid, sdi_nmono_id, &
            & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
       call check_netcdf( nf90_put_att(ncid, sdi_nmono_id, 'long_name', 'number of monomers (primary ice crystals)') )
       call check_netcdf( nf90_put_att(ncid, sdi_nmono_id, 'units', '') )
    end if

    !!! End of definition
    call check_netcdf( nf90_enddef(ncid) )

    ! Save data
    !!! sd_x
    call check_netcdf( nf90_put_var(ncid, sd_x_id, sd_x) )
    !!! sd_y
    call check_netcdf( nf90_put_var(ncid, sd_y_id, sd_y) )
    !!! sd_z
    call check_netcdf( nf90_put_var(ncid, sd_z_id, sd_z) )
    !!! sd_vz
    call check_netcdf( nf90_put_var(ncid, sd_vz_id, sd_vz) )
    !!! sd_r
    call check_netcdf( nf90_put_var(ncid, sd_r_id, sd_r) )
    !!! sd_asl
    call check_netcdf( nf90_put_var(ncid, sd_asl_id, sd_asl) )

    !!! sd_n
    call check_netcdf( nf90_put_var(ncid, sd_n_id, sd_n) )
    !!! sd_liqice
    call check_netcdf( nf90_put_var(ncid, sd_liqice_id, sd_liqice) )

    if( write_tracking ) then
      call check_netcdf( nf90_put_var(ncid, sd_id_id, sd_id) )
      call check_netcdf( nf90_put_var(ncid, domain_id, dm_id) )
    end if
    if( write_coal ) then
      call check_netcdf( nf90_put_var(ncid, if_coal_id, if_coal) )
    end if
    if( write_cold_tracking ) then
      call check_netcdf( nf90_put_var(ncid, sd_event_mask_id, sd_event_mask) )
      call check_netcdf( nf90_put_var(ncid, sd_event_sig_mask_id, sd_event_sig_mask) )
      call check_netcdf( nf90_put_var(ncid, sd_diag_mask_id, sd_diag_mask) )
      call check_netcdf( nf90_put_var(ncid, sd_phase_change_flag_id, sd_phase_change_flag) )
      call check_netcdf( nf90_put_var(ncid, sd_spatial_visit_flag_id, sd_spatial_visit_flag) )
      call check_netcdf( nf90_put_var(ncid, sd_liq_radius_max_interval_id, sd_liq_radius_max_interval) )
      call check_netcdf( nf90_put_var(ncid, sd_ice_rvol_max_interval_id, sd_ice_rvol_max_interval) )
      call check_netcdf( nf90_put_var(ncid, sd_mixed_rvol_max_interval_id, sd_mixed_rvol_max_interval) )
      call check_netcdf( nf90_put_var(ncid, sd_rime_mass_max_interval_id, sd_rime_mass_max_interval) )
      call check_netcdf( nf90_put_var(ncid, sd_rime_frac_max_interval_id, sd_rime_frac_max_interval) )
      call check_netcdf( nf90_put_var(ncid, sd_nmono_max_interval_id, sd_nmono_max_interval) )
      call check_netcdf( nf90_put_var(ncid, sd_aspect_ratio_max_interval_id, sd_aspect_ratio_max_interval) )
    end if

    if( sdm_cold ) then
       !!! sdi_re
       call check_netcdf( nf90_put_var(ncid, sdi_re_id, sdi%re(1:sd_num)) )
       !!! sdi_rp
       call check_netcdf( nf90_put_var(ncid, sdi_rp_id, sdi%rp(1:sd_num)) )
       !!! sdi_rho
       call check_netcdf( nf90_put_var(ncid, sdi_rho_id, sdi%rho(1:sd_num)) )
       !!! sdi_tf
       call check_netcdf( nf90_put_var(ncid, sdi_tf_id, sdi%tf(1:sd_num)) )
       !!! sdi_mrime
       call check_netcdf( nf90_put_var(ncid, sdi_mrime_id, sdi%mrime(1:sd_num)) )
       !!! sdi_nmono
       call check_netcdf( nf90_put_var(ncid, sdi_nmono_id, sdi%nmono(1:sd_num)) )
    end if

    ! Close the output file
    call check_netcdf( nf90_close(ncid) )

    if( IO_L ) write(IO_FID_LOG,*) '*** Closed output file (NetCDF) of Super Droplet'

    if( tracking_mode == 1 .and. present(filetag) ) then
      if( trim(filetag) == 'selected' ) then
        call sdm_write_tracking_id_file(basename_sd_out, sd_num, sd_id, dm_id)
      end if
    end if

    return

  end subroutine sdm_outnetcdf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sdm_outnetcdf_hist(otime,sd_num,sd_numasl,sd_n,sd_liqice,sd_x,sd_y,sd_z,sd_r,sd_asl,sd_vz,sdi,sd_id,dm_id,if_coal, &
       &                        sdn_dmpnskip,filetag,sd_event_mask,sd_event_sig_mask,sd_diag_mask,sd_phase_change_flag,                       &
       &                        sd_spatial_visit_flag,                                                                       &
       &                        sd_liq_radius_max_interval,sd_ice_rvol_max_interval,sd_mixed_rvol_max_interval,              &
       &                        sd_rime_mass_max_interval,sd_rime_frac_max_interval,sd_nmono_max_interval,                   &
       &                        sd_aspect_ratio_max_interval)
    use netcdf
    use scale_precision
    use scale_stdio
    use scale_time
    use scale_process, only: &
         mype => PRC_myrank, &
         PRC_MPIstop
    use m_sdm_common, only: &
         i2, sdm_cold, STAT_LIQ, STAT_ICE, STAT_MIX, sdicedef, &
         INVALID_i4, tracking_mode, coalescence_output_enable

    implicit none

    real(DP), intent(in) :: otime
    integer, intent(in) :: sd_num    ! number of super-droplets
    integer, intent(in) :: sd_numasl ! number of chemical species contained in super droplets
    integer(DP), intent(in) :: sd_n(1:sd_num) ! multiplicity of super-droplets
    integer(kind=i2), intent(in) :: sd_liqice(1:sd_num)
                       ! status of super-droplets (liquid/ice)
                       ! 01 = all liquid, 10 = all ice
                       ! 11 = mixture of ice and liquid
    real(RP), intent(in) :: sd_x(1:sd_num) ! x-coordinate of super-droplets
    real(RP), intent(in) :: sd_y(1:sd_num) ! y-coordinate of super-droplets
    real(RP), intent(in) :: sd_z(1:sd_num) ! z-coordinate of super-droplets
    real(RP), intent(in) :: sd_r(1:sd_num) ! equivalent radius of super-droplets
    real(RP), intent(in) :: sd_asl(1:sd_num,1:sd_numasl) ! aerosol mass of super-droplets
    real(RP), intent(in) :: sd_vz(1:sd_num) ! terminal velocity of super-droplets
    type(sdicedef), intent(in) :: sdi   ! ice phase super-droplets
    integer, intent(inout) :: sd_id(1:sd_num)
    integer, intent(inout) :: dm_id(1:sd_num)
    integer(kind=i2), intent(inout) :: if_coal(1:sd_num)
    integer, intent(in) :: sdn_dmpnskip ! Base skip to store super droplets in text format
    character(len=*),intent(in),optional :: filetag ! user defined text string tag to be added to the filenames
    integer, intent(in), optional :: sd_event_mask(1:sd_num)
    integer, intent(in), optional :: sd_event_sig_mask(1:sd_num)
    integer, intent(in), optional :: sd_diag_mask(1:sd_num)
    integer, intent(in), optional :: sd_phase_change_flag(1:sd_num)
    integer, intent(in), optional :: sd_spatial_visit_flag(1:sd_num)
    real(RP), intent(in), optional :: sd_liq_radius_max_interval(1:sd_num)
    real(RP), intent(in), optional :: sd_ice_rvol_max_interval(1:sd_num)
    real(RP), intent(in), optional :: sd_mixed_rvol_max_interval(1:sd_num)
    real(RP), intent(in), optional :: sd_rime_mass_max_interval(1:sd_num)
    real(RP), intent(in), optional :: sd_rime_frac_max_interval(1:sd_num)
    real(RP), intent(in), optional :: sd_nmono_max_interval(1:sd_num)
    real(RP), intent(in), optional :: sd_aspect_ratio_max_interval(1:sd_num)

    character(len=17) :: fmt2="(A, '.', A, I*.*)"
    character(len=17) :: fmt3="(3A)"
    character(len=H_LONG) :: ftmp
    character(len=H_LONG) :: basename_sd_out
    character(len=19) :: basename_time
    integer :: fid_sdm_o
    integer :: n, m, ierr
    character(len=80) :: fmt     ! output formate
    character(len=5)  :: cstat   ! status character
    integer :: nf90_real_precision
    integer :: ncid, sd_num_id, sd_numasl_id
    integer :: sd_x_id, sd_y_id, sd_z_id, sd_vz_id, sd_r_id, sd_asl_id, sd_n_id, sd_liqice_id
    integer :: sd_id_id, domain_id, if_coal_id
    integer :: sd_event_mask_id, sd_event_sig_mask_id, sd_diag_mask_id, sd_phase_change_flag_id, sd_spatial_visit_flag_id
    integer :: sd_liq_radius_max_interval_id, sd_ice_rvol_max_interval_id, sd_mixed_rvol_max_interval_id
    integer :: sd_rime_mass_max_interval_id, sd_rime_frac_max_interval_id
    integer :: sd_nmono_max_interval_id, sd_aspect_ratio_max_interval_id
    integer :: sdi_re_id, sdi_rp_id, sdi_rho_id, sdi_tf_id, sdi_mrime_id, sdi_nmono_id
    integer :: sd_dmp_time_idx_id, sd_dmp_time_id

    integer,parameter :: nc_deflate_level = 1      ! NetCDF compression level {1,..,9}
    integer,parameter :: nc_deflate       = 1 ! turn on NetCDF compresion
    integer,parameter :: nc_shuffle       = 1 ! turn on NetCDF shuffle filter

    logical :: newfile
    character(len=100) :: var_name
    character(len=100) :: ftag ! =filetag or ''(default)
    integer, parameter :: max_filenum = 2
    integer, save :: filenum = 0
    character(len=100), save :: ftag_list(1:max_filenum)
    integer, save :: time_count(1:max_filenum)
    integer :: nf, fileid
    logical :: write_forward_tracking, write_backward_tracking, write_tracking, write_coal, write_cold_tracking
    character(len=80) :: tracking_id_label

    !--- output Super Droplets in NetCDF format
    call TIME_gettimelabel(basename_time)
    fid_sdm_o = IO_get_available_fid()

    if(present(filetag)) then
       write(ftag,'(2A)') 'SD_', trim(filetag)
    else
       ftag = 'SD_output'
    end if

    newfile = .true.
    if(filenum /= 0) then
       do nf=1,filenum
          if(trim(ftag_list(nf)) == trim(ftag)) then
             fileid=nf
             newfile = .false.
             exit
          end if
       end do
    end if

    if(newfile) then
       if(filenum == max_filenum) then
          write(*,*) "sdm_outnetcdf_hist: ", "Exceeds the maximum file tag numbers allowed.", &
               &     "Consider increasing 'max_filenum' in sdm_io.f90"
          fileid = 1 ! Provide a safe fallback before aborting to avoid undefined behavior
          call PRC_MPIstop
       else
          filenum = filenum + 1
          fileid = filenum
          write(ftag_list(fileid),'(A)') trim(ftag)
       end if
    end if

    write(fmt2(14:14),'(I1)') 6
    write(fmt2(16:16),'(I1)') 6
    write(ftmp,fmt3) trim(ftag), '_history'
    write(basename_sd_out,fmt2) trim(ftmp), 'pe',mype

    ! Check the presision
    if(RP == SP)then
       nf90_real_precision = NF90_FLOAT
    else
       nf90_real_precision = NF90_DOUBLE
    end if

    ! Create or Open the output file
    if(newfile) then
       time_count(fileid) = 1
       call check_netcdf( nf90_create(trim(basename_sd_out),NF90_NETCDF4, ncid) ) 
    else
       time_count(fileid) = time_count(fileid) + 1
       call check_netcdf( nf90_open(trim(basename_sd_out),NF90_WRITE, ncid) )
       call check_netcdf( nf90_redef(ncid) )
    end if

    ! Definition of dimensions and variables
    !!! sd_dmp_time_idx
    if( NF90_NOERR .ne. nf90_inq_dimid(ncid, "sd_dmp_time_idx", sd_dmp_time_idx_id) ) then
       call check_netcdf( nf90_def_dim(ncid, "sd_dmp_time_idx", NF90_UNLIMITED, sd_dmp_time_idx_id) )
    end if
    !!! sd_numasl
    if( NF90_NOERR .ne. nf90_inq_dimid(ncid, "sd_numasl", sd_numasl_id) ) then
       call check_netcdf( nf90_def_dim(ncid, "sd_numasl", sd_numasl, sd_numasl_id) )
    end if
    !!! sd_num
    write(var_name,fmt='(A,I0.4)') "sd_num_", time_count(fileid)
    var_name = trim(var_name)
    call check_netcdf( nf90_def_dim(ncid, var_name, sd_num, sd_num_id) )

    !!! sd_dmp_time
    if( NF90_NOERR .ne. nf90_inq_varid(ncid, "sd_dmp_time", sd_dmp_time_id) ) then
       call check_netcdf( nf90_def_var(ncid, "sd_dmp_time", NF90_DOUBLE, sd_dmp_time_idx_id, sd_dmp_time_id) )
       call check_netcdf( nf90_def_var_deflate(ncid, sd_dmp_time_id, &
            & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
       call check_netcdf( nf90_put_att(ncid, sd_dmp_time_id, 'long_name', 'time') )
       call check_netcdf( nf90_put_att(ncid, sd_dmp_time_id, 'units', 's') )
    end if

    !!! sd_x
    write(var_name,fmt='(A,I0.4)') "sd_x_", time_count(fileid)
    var_name = trim(var_name)
    call check_netcdf( nf90_def_var(ncid, var_name, nf90_real_precision, sd_num_id, sd_x_id) )
    call check_netcdf( nf90_def_var_deflate(ncid, sd_x_id, &
         & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
    call check_netcdf( nf90_put_att(ncid, sd_x_id, 'long_name', 'x-coordinate') )
    call check_netcdf( nf90_put_att(ncid, sd_x_id, 'units', 'm') )
    !!! sd_y
    write(var_name,fmt='(A,I0.4)') "sd_y_", time_count(fileid)
    var_name = trim(var_name)
    call check_netcdf( nf90_def_var(ncid, var_name, nf90_real_precision, sd_num_id, sd_y_id) )
    call check_netcdf( nf90_def_var_deflate(ncid, sd_y_id, &
         & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
    call check_netcdf( nf90_put_att(ncid, sd_y_id, 'long_name', 'y-coordinate') )
    call check_netcdf( nf90_put_att(ncid, sd_y_id, 'units', 'm') )
    !!! sd_z
    write(var_name,fmt='(A,I0.4)') "sd_z_", time_count(fileid)
    var_name = trim(var_name)
    call check_netcdf( nf90_def_var(ncid, var_name, nf90_real_precision, sd_num_id, sd_z_id) )
    call check_netcdf( nf90_def_var_deflate(ncid, sd_z_id, &
         & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
    call check_netcdf( nf90_put_att(ncid, sd_z_id, 'long_name', 'z-coordinate') )
    call check_netcdf( nf90_put_att(ncid, sd_z_id, 'units', 'm') )
    !!! sd_vz
    write(var_name,fmt='(A,I0.4)') "sd_vz_", time_count(fileid)
    var_name = trim(var_name)
    call check_netcdf( nf90_def_var(ncid, var_name, nf90_real_precision, sd_num_id, sd_vz_id) )
    call check_netcdf( nf90_def_var_deflate(ncid, sd_vz_id, &
         & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
    call check_netcdf( nf90_put_att(ncid, sd_vz_id, 'long_name', 'terminal velocity') )
    call check_netcdf( nf90_put_att(ncid, sd_vz_id, 'units', 'm/s') )
    !!! sd_r
    write(var_name,fmt='(A,I0.4)') "sd_r_", time_count(fileid)
    var_name = trim(var_name)
    call check_netcdf( nf90_def_var(ncid, var_name, nf90_real_precision, sd_num_id, sd_r_id) )
    call check_netcdf( nf90_def_var_deflate(ncid, sd_r_id, &
         & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
    call check_netcdf( nf90_put_att(ncid, sd_r_id, 'long_name', 'equivalent radius of liquid droplets') )
    call check_netcdf( nf90_put_att(ncid, sd_r_id, 'units', 'm') )
    !!! sd_asl
    write(var_name,fmt='(A,I0.4)') "sd_asl_", time_count(fileid)
    var_name = trim(var_name)
    call check_netcdf( nf90_def_var(ncid, var_name, nf90_real_precision, (/sd_num_id, sd_numasl_id/), sd_asl_id) )
    call check_netcdf( nf90_def_var_deflate(ncid, sd_asl_id, &
         & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
    call check_netcdf( nf90_put_att(ncid, sd_asl_id, 'long_name', 'aerosol mass') )
    call check_netcdf( nf90_put_att(ncid, sd_asl_id, 'units', 'g') )

    !!! sd_n
    write(var_name,fmt='(A,I0.4)') "sd_n_", time_count(fileid)
    var_name = trim(var_name)
    call check_netcdf( nf90_def_var(ncid, var_name, NF90_INT64, sd_num_id, sd_n_id) )
    call check_netcdf( nf90_def_var_deflate(ncid, sd_n_id, &
         & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
    call check_netcdf( nf90_put_att(ncid, sd_n_id, 'long_name', 'multiplicity') )
    call check_netcdf( nf90_put_att(ncid, sd_n_id, 'units', '') )
    !!! sd_liqice
    write(var_name,fmt='(A,I0.4)') "sd_liqice_", time_count(fileid)
    var_name = trim(var_name)
    call check_netcdf( nf90_def_var(ncid, var_name, NF90_SHORT, sd_num_id, sd_liqice_id) )
    call check_netcdf( nf90_def_var_deflate(ncid, sd_liqice_id, &
         & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
    call check_netcdf( nf90_put_att(ncid, sd_liqice_id, 'long_name', 'status of droplets: 01=liquid, 10=ice, 11=mixture') )
    call check_netcdf( nf90_put_att(ncid, sd_liqice_id, 'units', '') )

    write_forward_tracking = tracking_mode == 1
    write_backward_tracking = tracking_mode == 2
    write_tracking = write_forward_tracking .or. write_backward_tracking
    write_coal = coalescence_output_enable .and. (.not. sdm_cold)
    write_cold_tracking = sdm_cold .and. present(sd_event_mask) .and. present(sd_event_sig_mask) .and. &
         &                 present(sd_diag_mask) .and. &
         &                 present(sd_phase_change_flag) .and. present(sd_spatial_visit_flag) .and. &
         &                 present(sd_liq_radius_max_interval) .and. &
         &                 present(sd_ice_rvol_max_interval) .and. present(sd_mixed_rvol_max_interval) .and. &
         &                 present(sd_rime_mass_max_interval) .and. present(sd_rime_frac_max_interval) .and. &
         &                 present(sd_nmono_max_interval) .and. present(sd_aspect_ratio_max_interval)

    if( write_backward_tracking ) then
      tracking_id_label = 'pre_sdid+pre_dmid'
      write(var_name,fmt='(A,I0.4)') "pre_sdid_", time_count(fileid)
      var_name = trim(var_name)
      call check_netcdf( nf90_def_var(ncid, var_name, NF90_INT, sd_num_id, sd_id_id) )
      call check_netcdf( nf90_def_var_deflate(ncid, sd_id_id, &
           & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
      call check_netcdf( nf90_put_att(ncid, sd_id_id, 'long_name', trim(tracking_id_label)) )
      call check_netcdf( nf90_put_att(ncid, sd_id_id, 'units', '') )
      write(var_name,fmt='(A,I0.4)') "pre_dmid_", time_count(fileid)
      var_name = trim(var_name)
      call check_netcdf( nf90_def_var(ncid, var_name, NF90_INT, sd_num_id, domain_id) )
      call check_netcdf( nf90_def_var_deflate(ncid, domain_id, &
           & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
      call check_netcdf( nf90_put_att(ncid, domain_id, 'long_name', trim(tracking_id_label)) )
      call check_netcdf( nf90_put_att(ncid, domain_id, 'units', '') )
    else if( write_forward_tracking ) then
      write(var_name,fmt='(A,I0.4)') "sd_id_", time_count(fileid)
      var_name = trim(var_name)
      call check_netcdf( nf90_def_var(ncid, var_name, NF90_INT, sd_num_id, sd_id_id) )
      call check_netcdf( nf90_def_var_deflate(ncid, sd_id_id, &
           & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
      call check_netcdf( nf90_put_att(ncid, sd_id_id, 'long_name', 'SD ID') )
      call check_netcdf( nf90_put_att(ncid, sd_id_id, 'units', '') )
      write(var_name,fmt='(A,I0.4)') "dm_id_", time_count(fileid)
      var_name = trim(var_name)
      call check_netcdf( nf90_def_var(ncid, var_name, NF90_INT, sd_num_id, domain_id) )
      call check_netcdf( nf90_def_var_deflate(ncid, domain_id, &
           & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
      call check_netcdf( nf90_put_att(ncid, domain_id, 'long_name', 'domain ID') )
      call check_netcdf( nf90_put_att(ncid, domain_id, 'units', '') )
    end if
    if( write_coal ) then
      write(var_name,fmt='(A,I0.4)') "if_coal_", time_count(fileid)
      var_name = trim(var_name)
      call check_netcdf( nf90_def_var(ncid, var_name, NF90_SHORT, sd_num_id, if_coal_id) )
      call check_netcdf( nf90_def_var_deflate(ncid, if_coal_id, &
           & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
      call check_netcdf( nf90_put_att(ncid, if_coal_id, 'long_name', 'coalescence flag: 0=has not occurred, 1=occurred') )
      call check_netcdf( nf90_put_att(ncid, if_coal_id, 'units', '') )
    end if
    if( write_cold_tracking ) then
      write(var_name,fmt='(A,I0.4)') "sd_event_mask_", time_count(fileid)
      var_name = trim(var_name)
      call check_netcdf( nf90_def_var(ncid, var_name, NF90_INT, sd_num_id, sd_event_mask_id) )
      call check_netcdf( nf90_def_var_deflate(ncid, sd_event_mask_id, &
           & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
      call check_netcdf( nf90_put_att(ncid, sd_event_mask_id, 'long_name', 'cold SD event bitmask over output interval') )
      call check_netcdf( nf90_put_att(ncid, sd_event_mask_id, 'units', '') )
      write(var_name,fmt='(A,I0.4)') "sd_event_sig_mask_", time_count(fileid)
      var_name = trim(var_name)
      call check_netcdf( nf90_def_var(ncid, var_name, NF90_INT, sd_num_id, sd_event_sig_mask_id) )
      call check_netcdf( nf90_def_var_deflate(ncid, sd_event_sig_mask_id, &
           & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
      call check_netcdf( nf90_put_att(ncid, sd_event_sig_mask_id, 'long_name', &
           'cold SD significant event subset bitmask over output interval') )
      call check_netcdf( nf90_put_att(ncid, sd_event_sig_mask_id, 'units', '') )
      write(var_name,fmt='(A,I0.4)') "sd_diag_mask_", time_count(fileid)
      var_name = trim(var_name)
      call check_netcdf( nf90_def_var(ncid, var_name, NF90_INT, sd_num_id, sd_diag_mask_id) )
      call check_netcdf( nf90_def_var_deflate(ncid, sd_diag_mask_id, &
           & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
      call check_netcdf( nf90_put_att(ncid, sd_diag_mask_id, 'long_name', 'cold SD diagnostic bitmask evaluated at output time') )
      call check_netcdf( nf90_put_att(ncid, sd_diag_mask_id, 'units', '') )
      write(var_name,fmt='(A,I0.4)') "sd_phase_change_flag_", time_count(fileid)
      var_name = trim(var_name)
      call check_netcdf( nf90_def_var(ncid, var_name, NF90_INT, sd_num_id, sd_phase_change_flag_id) )
      call check_netcdf( nf90_def_var_deflate(ncid, sd_phase_change_flag_id, &
           & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
      call check_netcdf( nf90_put_att(ncid, sd_phase_change_flag_id, 'long_name', 'cold phase-state category-change flag over output interval') )
      call check_netcdf( nf90_put_att(ncid, sd_phase_change_flag_id, 'units', '') )
      write(var_name,fmt='(A,I0.4)') "sd_spatial_visit_flag_", time_count(fileid)
      var_name = trim(var_name)
      call check_netcdf( nf90_def_var(ncid, var_name, NF90_INT, sd_num_id, sd_spatial_visit_flag_id) )
      call check_netcdf( nf90_def_var_deflate(ncid, sd_spatial_visit_flag_id, &
           & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
      call check_netcdf( nf90_put_att(ncid, sd_spatial_visit_flag_id, 'long_name', &
           'cold spatial-visit flag over output interval') )
      call check_netcdf( nf90_put_att(ncid, sd_spatial_visit_flag_id, 'units', '') )
      write(var_name,fmt='(A,I0.4)') "sd_liq_radius_max_interval_", time_count(fileid)
      var_name = trim(var_name)
      call check_netcdf( nf90_def_var(ncid, var_name, nf90_real_precision, sd_num_id, sd_liq_radius_max_interval_id) )
      call check_netcdf( nf90_def_var_deflate(ncid, sd_liq_radius_max_interval_id, &
           & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
      call check_netcdf( nf90_put_att(ncid, sd_liq_radius_max_interval_id, 'long_name', 'maximum liquid radius over output interval') )
      call check_netcdf( nf90_put_att(ncid, sd_liq_radius_max_interval_id, 'units', 'm') )
      write(var_name,fmt='(A,I0.4)') "sd_ice_rvol_max_interval_", time_count(fileid)
      var_name = trim(var_name)
      call check_netcdf( nf90_def_var(ncid, var_name, nf90_real_precision, sd_num_id, sd_ice_rvol_max_interval_id) )
      call check_netcdf( nf90_def_var_deflate(ncid, sd_ice_rvol_max_interval_id, &
           & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
      call check_netcdf( nf90_put_att(ncid, sd_ice_rvol_max_interval_id, 'long_name', 'maximum ice volume-equivalent radius over output interval') )
      call check_netcdf( nf90_put_att(ncid, sd_ice_rvol_max_interval_id, 'units', 'm') )
      write(var_name,fmt='(A,I0.4)') "sd_mixed_rvol_max_interval_", time_count(fileid)
      var_name = trim(var_name)
      call check_netcdf( nf90_def_var(ncid, var_name, nf90_real_precision, sd_num_id, sd_mixed_rvol_max_interval_id) )
      call check_netcdf( nf90_def_var_deflate(ncid, sd_mixed_rvol_max_interval_id, &
           & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
      call check_netcdf( nf90_put_att(ncid, sd_mixed_rvol_max_interval_id, 'long_name', 'maximum mixed-phase volume-equivalent radius over output interval') )
      call check_netcdf( nf90_put_att(ncid, sd_mixed_rvol_max_interval_id, 'units', 'm') )
      write(var_name,fmt='(A,I0.4)') "sd_rime_mass_max_interval_", time_count(fileid)
      var_name = trim(var_name)
      call check_netcdf( nf90_def_var(ncid, var_name, nf90_real_precision, sd_num_id, sd_rime_mass_max_interval_id) )
      call check_netcdf( nf90_def_var_deflate(ncid, sd_rime_mass_max_interval_id, &
           & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
      call check_netcdf( nf90_put_att(ncid, sd_rime_mass_max_interval_id, 'long_name', 'maximum rime mass over output interval') )
      call check_netcdf( nf90_put_att(ncid, sd_rime_mass_max_interval_id, 'units', 'kg') )
      write(var_name,fmt='(A,I0.4)') "sd_rime_frac_max_interval_", time_count(fileid)
      var_name = trim(var_name)
      call check_netcdf( nf90_def_var(ncid, var_name, nf90_real_precision, sd_num_id, sd_rime_frac_max_interval_id) )
      call check_netcdf( nf90_def_var_deflate(ncid, sd_rime_frac_max_interval_id, &
           & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
      call check_netcdf( nf90_put_att(ncid, sd_rime_frac_max_interval_id, 'long_name', 'maximum rime mass fraction over output interval') )
      call check_netcdf( nf90_put_att(ncid, sd_rime_frac_max_interval_id, 'units', '') )
      write(var_name,fmt='(A,I0.4)') "sd_nmono_max_interval_", time_count(fileid)
      var_name = trim(var_name)
      call check_netcdf( nf90_def_var(ncid, var_name, nf90_real_precision, sd_num_id, sd_nmono_max_interval_id) )
      call check_netcdf( nf90_def_var_deflate(ncid, sd_nmono_max_interval_id, &
           & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
      call check_netcdf( nf90_put_att(ncid, sd_nmono_max_interval_id, 'long_name', 'maximum monomer count over output interval') )
      call check_netcdf( nf90_put_att(ncid, sd_nmono_max_interval_id, 'units', '') )
      write(var_name,fmt='(A,I0.4)') "sd_aspect_ratio_max_interval_", time_count(fileid)
      var_name = trim(var_name)
      call check_netcdf( nf90_def_var(ncid, var_name, nf90_real_precision, sd_num_id, sd_aspect_ratio_max_interval_id) )
      call check_netcdf( nf90_def_var_deflate(ncid, sd_aspect_ratio_max_interval_id, &
           & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
      call check_netcdf( nf90_put_att(ncid, sd_aspect_ratio_max_interval_id, 'long_name', 'maximum ice aspect ratio over output interval') )
      call check_netcdf( nf90_put_att(ncid, sd_aspect_ratio_max_interval_id, 'units', '') )
    end if

    if( sdm_cold ) then
       !!! sdi%re
       write(var_name,fmt='(A,I0.4)') "sdi_re_", time_count(fileid)
       var_name = trim(var_name)
       call check_netcdf( nf90_def_var(ncid, var_name, nf90_real_precision, sd_num_id, sdi_re_id) )
       call check_netcdf( nf90_def_var_deflate(ncid, sdi_re_id, &
            & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
       call check_netcdf( nf90_put_att(ncid, sdi_re_id, 'long_name', 'equatorial radius of ice crystals') )
       call check_netcdf( nf90_put_att(ncid, sdi_re_id, 'units', 'm') )
       !!! sdi%rp
       write(var_name,fmt='(A,I0.4)') "sdi_rp_", time_count(fileid)
       var_name = trim(var_name)
       call check_netcdf( nf90_def_var(ncid, var_name, nf90_real_precision, sd_num_id, sdi_rp_id) )
       call check_netcdf( nf90_def_var_deflate(ncid, sdi_rp_id, &
            & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
       call check_netcdf( nf90_put_att(ncid, sdi_rp_id, 'long_name', 'polar radius of ice crystals') )
       call check_netcdf( nf90_put_att(ncid, sdi_rp_id, 'units', 'm') )
       !!! sdi%rho
       write(var_name,fmt='(A,I0.4)') "sdi_rho_", time_count(fileid)
       var_name = trim(var_name)
       call check_netcdf( nf90_def_var(ncid, var_name, nf90_real_precision, sd_num_id, sdi_rho_id) )
       call check_netcdf( nf90_def_var_deflate(ncid, sdi_rho_id, &
            & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
       call check_netcdf( nf90_put_att(ncid, sdi_rho_id, 'long_name', 'density of ice crystals') )
       call check_netcdf( nf90_put_att(ncid, sdi_rho_id, 'units', 'kg/m3') )
       !!! sdi%tf
       write(var_name,fmt='(A,I0.4)') "sdi_tf_", time_count(fileid)
       var_name = trim(var_name)
       call check_netcdf( nf90_def_var(ncid, var_name, nf90_real_precision, sd_num_id, sdi_tf_id) )
       call check_netcdf( nf90_def_var_deflate(ncid, sdi_tf_id, &
            & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
       call check_netcdf( nf90_put_att(ncid, sdi_tf_id, 'long_name', 'freezing temperature of particles') )
       call check_netcdf( nf90_put_att(ncid, sdi_tf_id, 'units', 'degC') )
       !!! sdi%mrime
       write(var_name,fmt='(A,I0.4)') "sdi_mrime_", time_count(fileid)
       var_name = trim(var_name)
       call check_netcdf( nf90_def_var(ncid, var_name, nf90_real_precision, sd_num_id, sdi_mrime_id) )
       call check_netcdf( nf90_def_var_deflate(ncid, sdi_mrime_id, &
            & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
       call check_netcdf( nf90_put_att(ncid, sdi_mrime_id, 'long_name', 'rime mass') )
       call check_netcdf( nf90_put_att(ncid, sdi_mrime_id, 'units', 'kg') )
       !!! sdi%nmono
       write(var_name,fmt='(A,I0.4)') "sdi_nmono_", time_count(fileid)
       var_name = trim(var_name)
       call check_netcdf( nf90_def_var(ncid, var_name, NF90_INT, sd_num_id, sdi_nmono_id) )
       call check_netcdf( nf90_def_var_deflate(ncid, sdi_nmono_id, &
            & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
       call check_netcdf( nf90_put_att(ncid, sdi_nmono_id, 'long_name', 'number of monomers (primary ice crystals)') )
       call check_netcdf( nf90_put_att(ncid, sdi_nmono_id, 'units', '') )
    end if

    !!! End of definition
    call check_netcdf( nf90_enddef(ncid) )

    ! Save data
    !!! sd_dmp_time
    call check_netcdf( nf90_put_var(ncid, sd_dmp_time_id, otime, start=(/time_count(fileid)/)) )

    !!! sd_x
    call check_netcdf( nf90_put_var(ncid, sd_x_id, sd_x) )
    !!! sd_y
    call check_netcdf( nf90_put_var(ncid, sd_y_id, sd_y) )
    !!! sd_z
    call check_netcdf( nf90_put_var(ncid, sd_z_id, sd_z) )
    !!! sd_vz
    call check_netcdf( nf90_put_var(ncid, sd_vz_id, sd_vz) )
    !!! sd_r
    call check_netcdf( nf90_put_var(ncid, sd_r_id, sd_r) )
    !!! sd_asl
    call check_netcdf( nf90_put_var(ncid, sd_asl_id, sd_asl) )

    !!! sd_n
    call check_netcdf( nf90_put_var(ncid, sd_n_id, sd_n) )
    !!! sd_liqice
    call check_netcdf( nf90_put_var(ncid, sd_liqice_id, sd_liqice) )

    if( write_tracking ) then
      call check_netcdf( nf90_put_var(ncid, sd_id_id, sd_id) )
      call check_netcdf( nf90_put_var(ncid, domain_id, dm_id) )
    end if
    if( write_coal ) then
      call check_netcdf( nf90_put_var(ncid, if_coal_id, if_coal) )
    end if
    if( write_cold_tracking ) then
      call check_netcdf( nf90_put_var(ncid, sd_event_mask_id, sd_event_mask) )
      call check_netcdf( nf90_put_var(ncid, sd_event_sig_mask_id, sd_event_sig_mask) )
      call check_netcdf( nf90_put_var(ncid, sd_diag_mask_id, sd_diag_mask) )
      call check_netcdf( nf90_put_var(ncid, sd_phase_change_flag_id, sd_phase_change_flag) )
      call check_netcdf( nf90_put_var(ncid, sd_spatial_visit_flag_id, sd_spatial_visit_flag) )
      call check_netcdf( nf90_put_var(ncid, sd_liq_radius_max_interval_id, sd_liq_radius_max_interval) )
      call check_netcdf( nf90_put_var(ncid, sd_ice_rvol_max_interval_id, sd_ice_rvol_max_interval) )
      call check_netcdf( nf90_put_var(ncid, sd_mixed_rvol_max_interval_id, sd_mixed_rvol_max_interval) )
      call check_netcdf( nf90_put_var(ncid, sd_rime_mass_max_interval_id, sd_rime_mass_max_interval) )
      call check_netcdf( nf90_put_var(ncid, sd_rime_frac_max_interval_id, sd_rime_frac_max_interval) )
      call check_netcdf( nf90_put_var(ncid, sd_nmono_max_interval_id, sd_nmono_max_interval) )
      call check_netcdf( nf90_put_var(ncid, sd_aspect_ratio_max_interval_id, sd_aspect_ratio_max_interval) )
    end if

    if( sdm_cold ) then
       !!! sdi_re
       call check_netcdf( nf90_put_var(ncid, sdi_re_id, sdi%re(1:sd_num)) )
       !!! sdi_rp
       call check_netcdf( nf90_put_var(ncid, sdi_rp_id, sdi%rp(1:sd_num)) )
       !!! sdi_rho
       call check_netcdf( nf90_put_var(ncid, sdi_rho_id, sdi%rho(1:sd_num)) )
       !!! sdi_tf
       call check_netcdf( nf90_put_var(ncid, sdi_tf_id, sdi%tf(1:sd_num)) )
       !!! sdi_mrime
       call check_netcdf( nf90_put_var(ncid, sdi_mrime_id, sdi%mrime(1:sd_num)) )
       !!! sdi_nmono
       call check_netcdf( nf90_put_var(ncid, sdi_nmono_id, sdi%nmono(1:sd_num)) )
    end if

    ! Close the output file
    call check_netcdf( nf90_close(ncid) )

    if( IO_L ) write(IO_FID_LOG,*) '*** Closed output file (NetCDF_HIST) of Super Droplet'

    if( tracking_mode == 1 .and. present(filetag) ) then
      if( trim(filetag) == 'selected' ) then
        call sdm_write_tracking_id_file(basename_sd_out, sd_num, sd_id, dm_id)
      end if
    end if

    if( tracking_mode == 2 ) then
      if( sdm_cold ) then
        call sdm_assign_tracking_subset(sd_num, sd_z, sd_r, sd_id, dm_id, if_coal, sd_liqice, sdi)
      else
        call sdm_assign_tracking_subset(sd_num, sd_z, sd_r, sd_id, dm_id, if_coal)
      end if
    end if

    return

  end subroutine sdm_outnetcdf_hist
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sdm_write_tracking_id_file(basename_sd_out, sd_num, sd_id, dm_id)
    use scale_stdio
    use scale_process, only: &
         mype => PRC_myrank
    use scale_process, only: &
         PRC_MPIstop
    use m_sdm_common, only: &
         INVALID_i4
    use m_sdm_tracking_cold, only: &
         sdm_tracking_valid_id_pair

    implicit none

    character(len=*), intent(in) :: basename_sd_out
    integer, intent(in) :: sd_num
    integer, intent(in) :: sd_id(1:sd_num)
    integer, intent(in) :: dm_id(1:sd_num)

    integer :: fid, ierr, n
    character(len=H_LONG) :: tracking_id_filename

    tracking_id_filename = trim(basename_sd_out) // '.ids'
    fid = IO_get_available_fid()

    open(fid, file=trim(tracking_id_filename), access='sequential', status='replace', form='formatted', iostat=ierr)
    if( ierr /= 0 ) then
      write(*,*) 'sdm_write_tracking_id_file', 'Write error', trim(tracking_id_filename), mype
      call PRC_MPIstop
    end if

    do n = 1, sd_num
      if( .not. sdm_tracking_valid_id_pair(sd_id(n), dm_id(n)) ) cycle
      write(fid,'(I0,1X,I0)') dm_id(n), sd_id(n)
    end do

    close(fid)

    if( IO_L ) write(IO_FID_LOG,*) '*** Wrote tracking ID list:', trim(tracking_id_filename)

    return
  end subroutine sdm_write_tracking_id_file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sdm_append_tracking_id_pairs(tracking_id_filename, pair_cnt, dm_id, sd_id)
    use scale_stdio
    use scale_process, only: &
         mype => PRC_myrank
    use scale_process, only: &
         PRC_MPIstop

    implicit none

    character(len=*), intent(in) :: tracking_id_filename
    integer, intent(in) :: pair_cnt
    integer, intent(in) :: dm_id(:)
    integer, intent(in) :: sd_id(:)

    integer :: fid, ierr, n

    if( pair_cnt <= 0 ) return

    fid = IO_get_available_fid()
    open(fid, file=trim(tracking_id_filename), access='sequential', status='unknown', position='append', &
         form='formatted', iostat=ierr)
    if( ierr /= 0 ) then
      write(*,*) 'sdm_append_tracking_id_pairs', 'Write error', trim(tracking_id_filename), mype
      call PRC_MPIstop
    end if

    do n = 1, pair_cnt
      write(fid,'(I0,1X,I0)') dm_id(n), sd_id(n)
    end do

    close(fid)

    if( IO_L ) write(IO_FID_LOG,*) '*** Appended tracking ID list:', trim(tracking_id_filename), pair_cnt

    return
  end subroutine sdm_append_tracking_id_pairs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sdm_interest_id_outnetcdf(otime, sd_num, sd_r, sd_id, dm_id, if_coal, records_written, &
       &                               sd_liqice, sdi, sd_event_mask)
    use scale_precision
    use scale_stdio
    use scale_process, only: &
         mype => PRC_myrank, &
         PRC_nprocs
    use scale_process, only: &
         PRC_MPIstop
    use scale_rm_process, only: &
         PRC_NUM_X, PRC_NUM_Y
    use m_sdm_common, only: &
         i2, INVALID_i4, tracking_id_output_basename, sdicedef, STAT_ICE, STAT_MIX, &
         tracking_interest_radius_enable, tracking_interest_radius_threshold, &
         tracking_interest_ice_radius_enable, tracking_interest_ice_radius_threshold, &
         tracking_interest_ice_phase_enable, tracking_interest_rime_mass_enable, &
         tracking_interest_rime_mass_threshold, tracking_interest_coalescence_enable
            use m_sdm_idutil, only: &
                 sdm_tracking_effective_radius
            use m_sdm_tracking_cold, only: &
                 EVENT_LIQ_LIQ_COAL, EVENT_RIMING, EVENT_AGGREGATION, &
                 sdm_tracking_valid_id_pair

    implicit none

    real(DP), intent(in) :: otime
    integer, intent(in) :: sd_num
    real(RP), intent(in) :: sd_r(1:sd_num)
    integer, intent(in) :: sd_id(1:sd_num)
    integer, intent(in) :: dm_id(1:sd_num)
    integer(kind=i2), intent(in) :: if_coal(1:sd_num)
    integer, intent(out), optional :: records_written
            integer(kind=i2), intent(in), optional :: sd_liqice(1:sd_num)
            type(sdicedef), intent(in), optional :: sdi
            integer, intent(in), optional :: sd_event_mask(1:sd_num)

    character(len=H_LONG) :: tracking_id_filename
    character(len=H_LONG) :: tracking_id_text_filename
    integer :: n, fid, ierr
    integer :: candidate_count
            logical :: selected
            logical :: has_cold_tracking_fields
            logical :: has_cold_event_mask
            real(RP) :: tracking_radius
            integer :: collision_event_mask

    if( present(records_written) ) records_written = 0

    if( len_trim(tracking_id_output_basename) == 0 ) return
    if( .not. tracking_interest_radius_enable .and. .not. tracking_interest_coalescence_enable .and. &
         .not. tracking_interest_ice_radius_enable .and. .not. tracking_interest_ice_phase_enable .and. &
         .not. tracking_interest_rime_mass_enable ) return

            has_cold_tracking_fields = present(sd_liqice) .and. present(sdi)
            has_cold_event_mask = present(sd_event_mask)
            collision_event_mask = EVENT_LIQ_LIQ_COAL + EVENT_RIMING + EVENT_AGGREGATION

    if( len_trim(tracking_id_output_basename) >= 3 ) then
      if( tracking_id_output_basename(len_trim(tracking_id_output_basename)-2:len_trim(tracking_id_output_basename)) == '.nc' ) then
        tracking_id_filename = trim(adjustl(tracking_id_output_basename))
      else
        write(tracking_id_filename,'(A,".pe",I6.6,".nc")') trim(adjustl(tracking_id_output_basename)), mype
      end if
    else
      write(tracking_id_filename,'(A,".pe",I6.6,".nc")') trim(adjustl(tracking_id_output_basename)), mype
    end if

    if( len_trim(tracking_id_filename) >= 3 ) then
      if( tracking_id_filename(len_trim(tracking_id_filename)-2:len_trim(tracking_id_filename)) == '.nc' ) then
        tracking_id_text_filename = trim(tracking_id_filename(1:len_trim(tracking_id_filename)-3)) // '.ids'
      else
        tracking_id_text_filename = trim(tracking_id_filename) // '.ids'
      end if
    else
      tracking_id_text_filename = trim(tracking_id_filename) // '.ids'
    end if

    candidate_count = 0
    do n = 1, sd_num
      if( .not. sdm_tracking_valid_id_pair(sd_id(n), dm_id(n)) ) cycle
      selected = .false.
      if( tracking_interest_radius_enable ) then
        if( sd_r(n) >= tracking_interest_radius_threshold ) selected = .true.
      end if
      if( has_cold_tracking_fields ) then
        if( tracking_interest_ice_radius_enable .and. &
             ( sd_liqice(n) == STAT_ICE .or. sd_liqice(n) == STAT_MIX ) ) then
          tracking_radius = sdm_tracking_effective_radius(sd_r(n), sd_liqice(n), sdi%re(n), sdi%rp(n))
          if( tracking_radius >= tracking_interest_ice_radius_threshold ) selected = .true.
        end if
        if( tracking_interest_ice_phase_enable ) then
          if( sd_liqice(n) == STAT_ICE .or. sd_liqice(n) == STAT_MIX ) selected = .true.
        end if
        if( tracking_interest_rime_mass_enable ) then
          if( sdi%mrime(n) >= tracking_interest_rime_mass_threshold ) selected = .true.
        end if
      end if
              if( tracking_interest_coalescence_enable ) then
                if( has_cold_event_mask ) then
                  if( iand(sd_event_mask(n), collision_event_mask) /= 0 ) selected = .true.
                else
                  if( if_coal(n) > 0_i2 ) selected = .true.
                end if
              end if
      if( selected ) candidate_count = candidate_count + 1
    end do

    if( candidate_count <= 0 ) return

    fid = IO_get_available_fid()
    ! TPHT FW emits a raw per-rank ID stream here and relies on offline merge/dedup later.
    open(fid, file=trim(tracking_id_text_filename), access='sequential', status='unknown', position='append', &
         form='formatted', iostat=ierr)
    if( ierr /= 0 ) then
      write(*,*) 'sdm_interest_id_outnetcdf', 'Write error', trim(tracking_id_text_filename), mype
      call PRC_MPIstop
    end if

    write(fid,'(A,1X,I0,1X,I0,1X,I0)') '# TPHT_META', PRC_NUM_X, PRC_NUM_Y, PRC_nprocs

    do n = 1, sd_num
      if( .not. sdm_tracking_valid_id_pair(sd_id(n), dm_id(n)) ) cycle
      selected = .false.
      if( tracking_interest_radius_enable ) then
        if( sd_r(n) >= tracking_interest_radius_threshold ) selected = .true.
      end if
      if( has_cold_tracking_fields ) then
        if( tracking_interest_ice_radius_enable .and. &
             ( sd_liqice(n) == STAT_ICE .or. sd_liqice(n) == STAT_MIX ) ) then
          tracking_radius = sdm_tracking_effective_radius(sd_r(n), sd_liqice(n), sdi%re(n), sdi%rp(n))
          if( tracking_radius >= tracking_interest_ice_radius_threshold ) selected = .true.
        end if
        if( tracking_interest_ice_phase_enable ) then
          if( sd_liqice(n) == STAT_ICE .or. sd_liqice(n) == STAT_MIX ) selected = .true.
        end if
        if( tracking_interest_rime_mass_enable ) then
          if( sdi%mrime(n) >= tracking_interest_rime_mass_threshold ) selected = .true.
        end if
      end if
              if( tracking_interest_coalescence_enable ) then
                if( has_cold_event_mask ) then
                  if( iand(sd_event_mask(n), collision_event_mask) /= 0 ) selected = .true.
                else
                  if( if_coal(n) > 0_i2 ) selected = .true.
                end if
              end if
      if( .not. selected ) cycle

      write(fid,'(I0,1X,I0)') dm_id(n), sd_id(n)
      if( present(records_written) ) records_written = records_written + 1
    end do

    close(fid)

    return
          end subroutine sdm_interest_id_outnetcdf
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sdm_event_singleproc_outnetcdf(otime, sd_num, sd_id, dm_id, sd_x, sd_y, sd_z, &
       &                                    sd_r, sd_liqice, sdi, phase_pre,               &
       &                                    liq_mass_pre, ice_mass_pre, hydro_radius_pre,  &
       &                                    ice_re_pre, ice_rp_pre, ice_rho_pre,           &
       &                                    track_meltfreeze_events, track_liq_vapor_events, &
       &                                    track_ice_vapor_events, track_kohler_activation_events, &
       &                                    kohler_rcrit_pre, kohler_rcrit_post, &
       &                                    kohler_margin_pre, kohler_margin_post, &
       &                                    kohler_active_pre, kohler_active_post, &
       &                                    aerosol_total_mass_pre, aerosol_total_mass_post, &
       &                                    aerosol_kohler_solute_pre, aerosol_kohler_solute_post, &
       &                                    air_temperature, air_pressure, water_vapor_mixing_ratio)
    use netcdf
    use scale_precision
    use scale_stdio
    use scale_process, only: &
         mype => PRC_myrank
    use m_sdm_common, only: &
         i2, sdm_cold, sdm_dmpitvl, tracking_mode, TRACK_COLD_OUTPUT_ICE_GEOMETRY, &
         TRACK_COLD_OUTPUT_AEROSOL_CONTEXT, TRACK_COLD_OUTPUT_THERMO_CONTEXT, &
         TRACK_COLD_OUTPUT_KOHLER_CONTEXT, sdicedef, &
         tracking_evt_freezing_enable, tracking_evt_melting_enable, &
         tracking_evt_deposition_enable, tracking_evt_sublimation_enable, &
         tracking_evt_condensation_enable, tracking_evt_evaporation_enable, &
         tracking_evt_activation_enable, tracking_evt_deactivation_enable, &
         tracking_sig_deposition_enable, tracking_sig_deposition_relmass_threshold, &
         tracking_sig_sublimation_enable, tracking_sig_sublimation_relmass_threshold, &
         tracking_sig_condensation_enable, tracking_sig_condensation_relmass_threshold, &
         tracking_sig_evaporation_enable, tracking_sig_evaporation_relmass_threshold, &
         tracking_sig_freezing_enable, tracking_sig_freezing_relmass_threshold, &
         tracking_sig_melting_enable, tracking_sig_melting_relmass_threshold, &
         tracking_sig_activation_enable, tracking_sig_activation_radius_threshold, &
         tracking_sig_deactivation_enable, tracking_sig_deactivation_radius_threshold
    use m_sdm_tracking_cold, only: &
         PHASE_LIQUID, PHASE_ICE, &
         TRIG_PROC_FREEZING, TRIG_PROC_MELTING, TRIG_PROC_DEPOSITION, &
         TRIG_PROC_SUBLIMATION, TRIG_PROC_CONDENSATION, TRIG_PROC_EVAPORATION, &
         TRIG_PROC_ACTIVATION, TRIG_PROC_DEACTIVATION, &
         TRIG_LEVEL_OCCURRENCE, TRIG_LEVEL_SIGNIFICANT, TARGET_BY_EVENT, &
         sdm_cold_phase_state, &
         sdm_cold_liq_mass, sdm_cold_ice_mass, sdm_cold_hydro_radius

    implicit none

    real(DP), intent(in) :: otime
    integer, intent(in) :: sd_num
    integer, intent(in) :: sd_id(1:sd_num)
    integer, intent(in) :: dm_id(1:sd_num)
    real(RP), intent(in) :: sd_x(1:sd_num)
    real(RP), intent(in) :: sd_y(1:sd_num)
    real(RP), intent(in) :: sd_z(1:sd_num)
    real(RP), intent(in) :: sd_r(1:sd_num)
    integer(i2), intent(in) :: sd_liqice(1:sd_num)
    type(sdicedef), intent(in) :: sdi
    integer, intent(in) :: phase_pre(1:sd_num)
    real(RP), intent(in) :: liq_mass_pre(1:sd_num)
    real(RP), intent(in) :: ice_mass_pre(1:sd_num)
    real(RP), intent(in) :: hydro_radius_pre(1:sd_num)
    real(RP), intent(in) :: ice_re_pre(1:sd_num)
    real(RP), intent(in) :: ice_rp_pre(1:sd_num)
    real(RP), intent(in) :: ice_rho_pre(1:sd_num)
    logical, intent(in) :: track_meltfreeze_events
    logical, intent(in) :: track_liq_vapor_events
    logical, intent(in) :: track_ice_vapor_events
    logical, intent(in) :: track_kohler_activation_events
    real(RP), intent(in), optional :: kohler_rcrit_pre(1:sd_num)
    real(RP), intent(in), optional :: kohler_rcrit_post(1:sd_num)
    real(RP), intent(in), optional :: kohler_margin_pre(1:sd_num)
    real(RP), intent(in), optional :: kohler_margin_post(1:sd_num)
    integer, intent(in), optional :: kohler_active_pre(1:sd_num)
    integer, intent(in), optional :: kohler_active_post(1:sd_num)
    real(RP), intent(in), optional :: aerosol_total_mass_pre(1:sd_num)
    real(RP), intent(in), optional :: aerosol_total_mass_post(1:sd_num)
    real(RP), intent(in), optional :: aerosol_kohler_solute_pre(1:sd_num)
    real(RP), intent(in), optional :: aerosol_kohler_solute_post(1:sd_num)
    real(RP), intent(in), optional :: air_temperature(1:sd_num)
    real(RP), intent(in), optional :: air_pressure(1:sd_num)
    real(RP), intent(in), optional :: water_vapor_mixing_ratio(1:sd_num)

    character(len=17) :: fmt2="(A, '.', A, I*.*)"
    character(len=H_LONG) :: ftmp
    character(len=H_LONG) :: basename_sd_out
    character(len=19) :: basename_time
    integer :: nf90_real_precision
    integer :: ncid, event_dim_id, time_id, trigger_code_id, trigger_level_id, target_reason_mask_id
    integer :: phase_state_pre_id, phase_state_post_id, sd_id_id, dm_id_id
    integer :: x_id, y_id, z_id
    integer :: hydro_radius_pre_id, hydro_radius_post_id
    integer :: liq_mass_pre_id, liq_mass_post_id
    integer :: ice_mass_pre_id, ice_mass_post_id
    integer :: ice_re_pre_id, ice_re_post_id
    integer :: ice_rp_pre_id, ice_rp_post_id
    integer :: ice_rho_pre_id, ice_rho_post_id
    integer :: kohler_rcrit_pre_id, kohler_rcrit_post_id
    integer :: kohler_margin_pre_id, kohler_margin_post_id
    integer :: activated_state_pre_id, activated_state_post_id
    integer :: aerosol_total_mass_pre_id, aerosol_total_mass_post_id
    integer :: aerosol_kohler_solute_pre_id, aerosol_kohler_solute_post_id
    integer :: air_temperature_id, air_pressure_id, water_vapor_mixing_ratio_id
    integer :: event_offset, start1(1), count1(1)
    integer :: n, rec_count, rec, phase_post
    integer :: trigger_code, trigger_level
    real(RP) :: liq_mass_post, ice_mass_post
    real(RP) :: delta_liq_mass, delta_ice_mass, rel_mass_change
    real(DP) :: otime_bucket, otime_bucket_daysec
    integer :: otime_bucket_hh, otime_bucket_mm, otime_bucket_ss, otime_bucket_ms
    logical :: file_exists, write_tracking_ids, write_extended, write_kohler_context, has_kohler_context
    logical :: write_aerosol_context, has_aerosol_context
    logical :: write_thermo_context, has_thermo_context
    character(len=16) :: id_name, dm_name
    real(DP), allocatable :: time_out(:)
    integer, allocatable :: trigger_code_out(:), trigger_level_out(:), target_reason_mask_out(:)
    integer, allocatable :: phase_state_pre_out(:), phase_state_post_out(:), sd_id_out(:), dm_id_out(:)
    real(RP), allocatable :: x_out(:), y_out(:), z_out(:)
    real(RP), allocatable :: hydro_radius_pre_out(:), hydro_radius_post_out(:)
    real(RP), allocatable :: liq_mass_pre_out(:), liq_mass_post_out(:)
    real(RP), allocatable :: ice_mass_pre_out(:), ice_mass_post_out(:)
    real(RP), allocatable :: ice_re_pre_out(:), ice_re_post_out(:)
    real(RP), allocatable :: ice_rp_pre_out(:), ice_rp_post_out(:)
    real(RP), allocatable :: ice_rho_pre_out(:), ice_rho_post_out(:)
    real(RP), allocatable :: kohler_rcrit_pre_out(:), kohler_rcrit_post_out(:)
    real(RP), allocatable :: kohler_margin_pre_out(:), kohler_margin_post_out(:)
    integer, allocatable :: activated_state_pre_out(:), activated_state_post_out(:)
    real(RP), allocatable :: aerosol_total_mass_pre_out(:), aerosol_total_mass_post_out(:)
    real(RP), allocatable :: aerosol_kohler_solute_pre_out(:), aerosol_kohler_solute_post_out(:)
    real(RP), allocatable :: air_temperature_out(:), air_pressure_out(:), water_vapor_mixing_ratio_out(:)

    integer, parameter :: nc_deflate_level = 1
    integer, parameter :: nc_deflate       = 1
    integer, parameter :: nc_shuffle       = 1

    if( .not. sdm_cold ) return

    if(RP == SP)then
       nf90_real_precision = NF90_FLOAT
    else
       nf90_real_precision = NF90_DOUBLE
    end if

    write_kohler_context = TRACK_COLD_OUTPUT_KOHLER_CONTEXT
    has_kohler_context = present(kohler_rcrit_pre) .and. present(kohler_rcrit_post) .and. &
         present(kohler_margin_pre) .and. present(kohler_margin_post) .and. &
         present(kohler_active_pre) .and. present(kohler_active_post)
    write_aerosol_context = TRACK_COLD_OUTPUT_AEROSOL_CONTEXT
    has_aerosol_context = present(aerosol_total_mass_pre) .and. present(aerosol_total_mass_post) .and. &
         present(aerosol_kohler_solute_pre) .and. present(aerosol_kohler_solute_post)
    write_thermo_context = TRACK_COLD_OUTPUT_THERMO_CONTEXT
    has_thermo_context = present(air_temperature) .and. present(air_pressure) .and. &
         present(water_vapor_mixing_ratio)

    rec_count = 0
    do n = 1, sd_num
       phase_post = sdm_cold_phase_state(sd_liqice(n), sd_r(n), sdi%re(n), sdi%rp(n))
       liq_mass_post = sdm_cold_liq_mass(sd_r(n))
       ice_mass_post = sdm_cold_ice_mass(sdi%re(n), sdi%rp(n), sdi%rho(n))
       delta_liq_mass = liq_mass_post - liq_mass_pre(n)
       delta_ice_mass = ice_mass_post - ice_mass_pre(n)
       if( track_meltfreeze_events ) then
          if( singleproc_writes_record(max(delta_ice_mass, -delta_liq_mass, 0.0_RP) > 0.0_RP, &
               tracking_evt_freezing_enable, tracking_sig_freezing_enable, &
               max(delta_ice_mass, -delta_liq_mass, 0.0_RP) / max(liq_mass_pre(n), tiny(1.0_RP)) >= &
               tracking_sig_freezing_relmass_threshold) ) rec_count = rec_count + 1
          if( singleproc_writes_record(max(delta_liq_mass, -delta_ice_mass, 0.0_RP) > 0.0_RP, &
               tracking_evt_melting_enable, tracking_sig_melting_enable, &
               max(delta_liq_mass, -delta_ice_mass, 0.0_RP) / max(ice_mass_pre(n), tiny(1.0_RP)) >= &
               tracking_sig_melting_relmass_threshold) ) rec_count = rec_count + 1
       end if
       if( track_liq_vapor_events ) then
          rel_mass_change = abs(delta_liq_mass) / max(liq_mass_pre(n), tiny(1.0_RP))
          if( singleproc_writes_record(delta_liq_mass > 0.0_RP, tracking_evt_condensation_enable, &
               tracking_sig_condensation_enable, rel_mass_change >= tracking_sig_condensation_relmass_threshold) ) &
               rec_count = rec_count + 1
	       if( singleproc_writes_record(delta_liq_mass < 0.0_RP, tracking_evt_evaporation_enable, &
	            tracking_sig_evaporation_enable, rel_mass_change >= tracking_sig_evaporation_relmass_threshold) ) &
	            rec_count = rec_count + 1
       end if
       if( track_kohler_activation_events .and. has_kohler_context ) then
          if( singleproc_writes_record(kohler_active_pre(n) == 0 .and. kohler_active_post(n) == 1, &
               tracking_evt_activation_enable, tracking_sig_activation_enable, &
               kohler_margin_post(n) >= tracking_sig_activation_radius_threshold) ) rec_count = rec_count + 1
          if( singleproc_writes_record(kohler_active_pre(n) == 1 .and. kohler_active_post(n) == 0, &
               tracking_evt_deactivation_enable, tracking_sig_deactivation_enable, &
               -kohler_margin_post(n) >= tracking_sig_deactivation_radius_threshold) ) rec_count = rec_count + 1
       end if
       ! Significant dep/sub records are thresholded process hits, not budgets.
       if( track_ice_vapor_events ) then
          rel_mass_change = abs(delta_ice_mass) / max(ice_mass_pre(n), tiny(1.0_RP))
          if( singleproc_writes_record(delta_ice_mass > 0.0_RP, tracking_evt_deposition_enable, &
               tracking_sig_deposition_enable, rel_mass_change >= tracking_sig_deposition_relmass_threshold) ) &
               rec_count = rec_count + 1
          if( singleproc_writes_record(delta_ice_mass < 0.0_RP, tracking_evt_sublimation_enable, &
               tracking_sig_sublimation_enable, rel_mass_change >= tracking_sig_sublimation_relmass_threshold) ) &
               rec_count = rec_count + 1
       end if
    end do

    if( rec_count <= 0 ) return

    allocate(time_out(rec_count))
    allocate(trigger_code_out(rec_count))
    allocate(trigger_level_out(rec_count))
    allocate(target_reason_mask_out(rec_count))
    allocate(phase_state_pre_out(rec_count))
    allocate(phase_state_post_out(rec_count))
    allocate(x_out(rec_count), y_out(rec_count), z_out(rec_count))
    allocate(hydro_radius_pre_out(rec_count), hydro_radius_post_out(rec_count))
    allocate(liq_mass_pre_out(rec_count), liq_mass_post_out(rec_count))
    allocate(ice_mass_pre_out(rec_count), ice_mass_post_out(rec_count))

    write_tracking_ids = tracking_mode /= 0
    write_extended = TRACK_COLD_OUTPUT_ICE_GEOMETRY
    if( write_tracking_ids ) then
       allocate(sd_id_out(rec_count))
       allocate(dm_id_out(rec_count))
    end if
    if( write_extended ) then
       allocate(ice_re_pre_out(rec_count), ice_re_post_out(rec_count))
       allocate(ice_rp_pre_out(rec_count), ice_rp_post_out(rec_count))
       allocate(ice_rho_pre_out(rec_count), ice_rho_post_out(rec_count))
    end if
    if( write_kohler_context ) then
       allocate(kohler_rcrit_pre_out(rec_count), kohler_rcrit_post_out(rec_count))
       allocate(kohler_margin_pre_out(rec_count), kohler_margin_post_out(rec_count))
       allocate(activated_state_pre_out(rec_count), activated_state_post_out(rec_count))
    end if
    if( write_aerosol_context ) then
       allocate(aerosol_total_mass_pre_out(rec_count), aerosol_total_mass_post_out(rec_count))
       allocate(aerosol_kohler_solute_pre_out(rec_count), aerosol_kohler_solute_post_out(rec_count))
    end if
    if( write_thermo_context ) then
       allocate(air_temperature_out(rec_count), air_pressure_out(rec_count), water_vapor_mixing_ratio_out(rec_count))
    end if

    rec = 0
    do n = 1, sd_num
       phase_post = sdm_cold_phase_state(sd_liqice(n), sd_r(n), sdi%re(n), sdi%rp(n))
       liq_mass_post = sdm_cold_liq_mass(sd_r(n))
       ice_mass_post = sdm_cold_ice_mass(sdi%re(n), sdi%rp(n), sdi%rho(n))
       delta_liq_mass = liq_mass_post - liq_mass_pre(n)
       delta_ice_mass = ice_mass_post - ice_mass_pre(n)
       if( track_meltfreeze_events ) then
          rel_mass_change = max(delta_ice_mass, -delta_liq_mass, 0.0_RP) / max(liq_mass_pre(n), tiny(1.0_RP))
          call maybe_put_singleproc_record(n, TRIG_PROC_FREEZING, &
               max(delta_ice_mass, -delta_liq_mass, 0.0_RP) > 0.0_RP, &
               tracking_evt_freezing_enable, tracking_sig_freezing_enable, &
               rel_mass_change >= tracking_sig_freezing_relmass_threshold, phase_post)
          rel_mass_change = max(delta_liq_mass, -delta_ice_mass, 0.0_RP) / max(ice_mass_pre(n), tiny(1.0_RP))
          call maybe_put_singleproc_record(n, TRIG_PROC_MELTING, &
               max(delta_liq_mass, -delta_ice_mass, 0.0_RP) > 0.0_RP, &
               tracking_evt_melting_enable, tracking_sig_melting_enable, &
               rel_mass_change >= tracking_sig_melting_relmass_threshold, phase_post)
       end if
       if( track_liq_vapor_events ) then
          rel_mass_change = abs(delta_liq_mass) / max(liq_mass_pre(n), tiny(1.0_RP))
          call maybe_put_singleproc_record(n, TRIG_PROC_CONDENSATION, delta_liq_mass > 0.0_RP, &
               tracking_evt_condensation_enable, tracking_sig_condensation_enable, &
               rel_mass_change >= tracking_sig_condensation_relmass_threshold, phase_post)
	       call maybe_put_singleproc_record(n, TRIG_PROC_EVAPORATION, delta_liq_mass < 0.0_RP, &
	            tracking_evt_evaporation_enable, tracking_sig_evaporation_enable, &
	            rel_mass_change >= tracking_sig_evaporation_relmass_threshold, phase_post)
	    end if
       if( track_kohler_activation_events .and. has_kohler_context ) then
          call maybe_put_singleproc_record(n, TRIG_PROC_ACTIVATION, &
               kohler_active_pre(n) == 0 .and. kohler_active_post(n) == 1, &
               tracking_evt_activation_enable, tracking_sig_activation_enable, &
               kohler_margin_post(n) >= tracking_sig_activation_radius_threshold, phase_post)
          call maybe_put_singleproc_record(n, TRIG_PROC_DEACTIVATION, &
               kohler_active_pre(n) == 1 .and. kohler_active_post(n) == 0, &
               tracking_evt_deactivation_enable, tracking_sig_deactivation_enable, &
               -kohler_margin_post(n) >= tracking_sig_deactivation_radius_threshold, phase_post)
       end if
       if( track_ice_vapor_events ) then
          rel_mass_change = abs(delta_ice_mass) / max(ice_mass_pre(n), tiny(1.0_RP))
          call maybe_put_singleproc_record(n, TRIG_PROC_DEPOSITION, delta_ice_mass > 0.0_RP, &
               tracking_evt_deposition_enable, tracking_sig_deposition_enable, &
               rel_mass_change >= tracking_sig_deposition_relmass_threshold, phase_post)
          call maybe_put_singleproc_record(n, TRIG_PROC_SUBLIMATION, delta_ice_mass < 0.0_RP, &
               tracking_evt_sublimation_enable, tracking_sig_sublimation_enable, &
               rel_mass_change >= tracking_sig_sublimation_relmass_threshold, phase_post)
       end if
    end do

    if( tracking_mode == 2 ) then
       id_name = 'pre_sdid'
       dm_name = 'pre_dmid'
    else if( tracking_mode == 1 ) then
       id_name = 'sd_id'
       dm_name = 'dm_id'
    else
       id_name = ''
       dm_name = ''
    end if

    if( sdm_dmpitvl > 0.0_RP ) then
       otime_bucket = sdm_dmpitvl * real( int( otime / sdm_dmpitvl + 1.0e-10_RP ), kind=DP )
    else
       otime_bucket = otime
    end if
    otime_bucket_daysec = mod( otime_bucket, 86400.0_DP )
    otime_bucket_hh = int( otime_bucket_daysec / 3600.0_DP )
    otime_bucket_mm = int( mod(otime_bucket_daysec, 3600.0_DP) / 60.0_DP )
    otime_bucket_ss = int( mod(otime_bucket_daysec, 60.0_DP) )
    otime_bucket_ms = nint( (otime_bucket_daysec - real(int(otime_bucket_daysec),kind=DP)) * 1000.0_DP )
    if( otime_bucket_ms >= 1000 ) then
       otime_bucket_ms = otime_bucket_ms - 1000
       otime_bucket_ss = otime_bucket_ss + 1
    end if
    write(basename_time,'(I8.8,A1,I2.2,I2.2,I2.2,A1,I3.3)') 101, '-', &
         otime_bucket_hh, otime_bucket_mm, otime_bucket_ss, '.', otime_bucket_ms
    write(fmt2(14:14),'(I1)') 6
    write(fmt2(16:16),'(I1)') 6
    write(ftmp,'(3A)') 'SD_event_singleproc', '_NetCDF_', trim(basename_time)
    write(basename_sd_out,fmt2) trim(ftmp), 'pe',mype

    inquire(file=trim(basename_sd_out), exist=file_exists)
    if( file_exists ) then
       call check_netcdf( nf90_open(trim(basename_sd_out), NF90_WRITE, ncid) )
       call check_netcdf( nf90_inq_dimid(ncid, "event", event_dim_id) )
       call check_netcdf( nf90_inquire_dimension(ncid, event_dim_id, len=event_offset) )
       call check_netcdf( nf90_inq_varid(ncid, "time", time_id) )
       call check_netcdf( nf90_inq_varid(ncid, "trigger_code", trigger_code_id) )
       call check_netcdf( nf90_inq_varid(ncid, "trigger_level", trigger_level_id) )
       call check_netcdf( nf90_inq_varid(ncid, "target_reason_mask", target_reason_mask_id) )
       call check_netcdf( nf90_inq_varid(ncid, "phase_state_pre", phase_state_pre_id) )
       call check_netcdf( nf90_inq_varid(ncid, "phase_state_post", phase_state_post_id) )
       if( write_tracking_ids ) then
          call check_netcdf( nf90_inq_varid(ncid, trim(id_name), sd_id_id) )
          call check_netcdf( nf90_inq_varid(ncid, trim(dm_name), dm_id_id) )
       end if
       call check_netcdf( nf90_inq_varid(ncid, "x", x_id) )
       call check_netcdf( nf90_inq_varid(ncid, "y", y_id) )
       call check_netcdf( nf90_inq_varid(ncid, "z", z_id) )
       call check_netcdf( nf90_inq_varid(ncid, "hydro_radius_pre", hydro_radius_pre_id) )
       call check_netcdf( nf90_inq_varid(ncid, "hydro_radius_post", hydro_radius_post_id) )
       call check_netcdf( nf90_inq_varid(ncid, "liq_mass_pre", liq_mass_pre_id) )
       call check_netcdf( nf90_inq_varid(ncid, "liq_mass_post", liq_mass_post_id) )
       call check_netcdf( nf90_inq_varid(ncid, "ice_mass_pre", ice_mass_pre_id) )
       call check_netcdf( nf90_inq_varid(ncid, "ice_mass_post", ice_mass_post_id) )
	       if( write_extended ) then
	          call check_netcdf( nf90_inq_varid(ncid, "ice_re_pre", ice_re_pre_id) )
	          call check_netcdf( nf90_inq_varid(ncid, "ice_re_post", ice_re_post_id) )
	          call check_netcdf( nf90_inq_varid(ncid, "ice_rp_pre", ice_rp_pre_id) )
	          call check_netcdf( nf90_inq_varid(ncid, "ice_rp_post", ice_rp_post_id) )
	          call check_netcdf( nf90_inq_varid(ncid, "ice_rho_pre", ice_rho_pre_id) )
	          call check_netcdf( nf90_inq_varid(ncid, "ice_rho_post", ice_rho_post_id) )
	       end if
       if( write_kohler_context ) then
          call check_netcdf( nf90_inq_varid(ncid, "kohler_rcrit_pre", kohler_rcrit_pre_id) )
          call check_netcdf( nf90_inq_varid(ncid, "kohler_rcrit_post", kohler_rcrit_post_id) )
          call check_netcdf( nf90_inq_varid(ncid, "kohler_margin_pre", kohler_margin_pre_id) )
          call check_netcdf( nf90_inq_varid(ncid, "kohler_margin_post", kohler_margin_post_id) )
          call check_netcdf( nf90_inq_varid(ncid, "activated_state_pre", activated_state_pre_id) )
          call check_netcdf( nf90_inq_varid(ncid, "activated_state_post", activated_state_post_id) )
       end if
       if( write_aerosol_context ) then
          call check_netcdf( nf90_inq_varid(ncid, "aerosol_total_mass_pre", aerosol_total_mass_pre_id) )
          call check_netcdf( nf90_inq_varid(ncid, "aerosol_total_mass_post", aerosol_total_mass_post_id) )
          call check_netcdf( nf90_inq_varid(ncid, "aerosol_kohler_solute_pre", aerosol_kohler_solute_pre_id) )
          call check_netcdf( nf90_inq_varid(ncid, "aerosol_kohler_solute_post", aerosol_kohler_solute_post_id) )
       end if
       if( write_thermo_context ) then
          call check_netcdf( nf90_inq_varid(ncid, "air_temperature", air_temperature_id) )
          call check_netcdf( nf90_inq_varid(ncid, "air_pressure", air_pressure_id) )
          call check_netcdf( nf90_inq_varid(ncid, "water_vapor_mixing_ratio", water_vapor_mixing_ratio_id) )
       end if
    else
       call check_netcdf( nf90_create(trim(basename_sd_out), NF90_NETCDF4, ncid) )
       call check_netcdf( nf90_def_dim(ncid, "event", NF90_UNLIMITED, event_dim_id) )
       call check_netcdf( nf90_def_var(ncid, "time", NF90_DOUBLE, event_dim_id, time_id) )
       call check_netcdf( nf90_def_var(ncid, "trigger_code", NF90_INT, event_dim_id, trigger_code_id) )
       call check_netcdf( nf90_def_var(ncid, "trigger_level", NF90_INT, event_dim_id, trigger_level_id) )
       call check_netcdf( nf90_def_var(ncid, "target_reason_mask", NF90_INT, event_dim_id, target_reason_mask_id) )
       call check_netcdf( nf90_def_var(ncid, "phase_state_pre", NF90_INT, event_dim_id, phase_state_pre_id) )
       call check_netcdf( nf90_def_var(ncid, "phase_state_post", NF90_INT, event_dim_id, phase_state_post_id) )
       if( write_tracking_ids ) then
          call check_netcdf( nf90_def_var(ncid, trim(id_name), NF90_INT, event_dim_id, sd_id_id) )
          call check_netcdf( nf90_def_var(ncid, trim(dm_name), NF90_INT, event_dim_id, dm_id_id) )
       end if
       call check_netcdf( nf90_def_var(ncid, "x", nf90_real_precision, event_dim_id, x_id) )
       call check_netcdf( nf90_def_var(ncid, "y", nf90_real_precision, event_dim_id, y_id) )
       call check_netcdf( nf90_def_var(ncid, "z", nf90_real_precision, event_dim_id, z_id) )
       call check_netcdf( nf90_def_var(ncid, "hydro_radius_pre", nf90_real_precision, event_dim_id, hydro_radius_pre_id) )
       call check_netcdf( nf90_def_var(ncid, "hydro_radius_post", nf90_real_precision, event_dim_id, hydro_radius_post_id) )
       call check_netcdf( nf90_def_var(ncid, "liq_mass_pre", nf90_real_precision, event_dim_id, liq_mass_pre_id) )
       call check_netcdf( nf90_def_var(ncid, "liq_mass_post", nf90_real_precision, event_dim_id, liq_mass_post_id) )
       call check_netcdf( nf90_def_var(ncid, "ice_mass_pre", nf90_real_precision, event_dim_id, ice_mass_pre_id) )
       call check_netcdf( nf90_def_var(ncid, "ice_mass_post", nf90_real_precision, event_dim_id, ice_mass_post_id) )
	       if( write_extended ) then
	          call check_netcdf( nf90_def_var(ncid, "ice_re_pre", nf90_real_precision, event_dim_id, ice_re_pre_id) )
	          call check_netcdf( nf90_def_var(ncid, "ice_re_post", nf90_real_precision, event_dim_id, ice_re_post_id) )
	          call check_netcdf( nf90_def_var(ncid, "ice_rp_pre", nf90_real_precision, event_dim_id, ice_rp_pre_id) )
	          call check_netcdf( nf90_def_var(ncid, "ice_rp_post", nf90_real_precision, event_dim_id, ice_rp_post_id) )
	          call check_netcdf( nf90_def_var(ncid, "ice_rho_pre", nf90_real_precision, event_dim_id, ice_rho_pre_id) )
	          call check_netcdf( nf90_def_var(ncid, "ice_rho_post", nf90_real_precision, event_dim_id, ice_rho_post_id) )
	       end if
       if( write_kohler_context ) then
          call check_netcdf( nf90_def_var(ncid, "kohler_rcrit_pre", nf90_real_precision, event_dim_id, kohler_rcrit_pre_id) )
          call check_netcdf( nf90_def_var(ncid, "kohler_rcrit_post", nf90_real_precision, event_dim_id, kohler_rcrit_post_id) )
          call check_netcdf( nf90_def_var(ncid, "kohler_margin_pre", nf90_real_precision, event_dim_id, kohler_margin_pre_id) )
          call check_netcdf( nf90_def_var(ncid, "kohler_margin_post", nf90_real_precision, event_dim_id, kohler_margin_post_id) )
          call check_netcdf( nf90_def_var(ncid, "activated_state_pre", NF90_INT, event_dim_id, activated_state_pre_id) )
          call check_netcdf( nf90_def_var(ncid, "activated_state_post", NF90_INT, event_dim_id, activated_state_post_id) )
       end if
       if( write_aerosol_context ) then
          call check_netcdf( nf90_def_var(ncid, "aerosol_total_mass_pre", nf90_real_precision, event_dim_id, &
               aerosol_total_mass_pre_id) )
          call check_netcdf( nf90_def_var(ncid, "aerosol_total_mass_post", nf90_real_precision, event_dim_id, &
               aerosol_total_mass_post_id) )
          call check_netcdf( nf90_def_var(ncid, "aerosol_kohler_solute_pre", nf90_real_precision, event_dim_id, &
               aerosol_kohler_solute_pre_id) )
          call check_netcdf( nf90_def_var(ncid, "aerosol_kohler_solute_post", nf90_real_precision, event_dim_id, &
               aerosol_kohler_solute_post_id) )
       end if
       if( write_thermo_context ) then
          call check_netcdf( nf90_def_var(ncid, "air_temperature", nf90_real_precision, event_dim_id, air_temperature_id) )
          call check_netcdf( nf90_def_var(ncid, "air_pressure", nf90_real_precision, event_dim_id, air_pressure_id) )
          call check_netcdf( nf90_def_var(ncid, "water_vapor_mixing_ratio", nf90_real_precision, event_dim_id, &
               water_vapor_mixing_ratio_id) )
       end if
       call check_netcdf( nf90_put_att(ncid, trigger_code_id, 'long_name', 'cold-SDM single-process trigger code') )
       call check_netcdf( nf90_put_att(ncid, trigger_level_id, 'long_name', &
            'cold-SDM process trigger level: 1=occurrence, 2=significant') )
       call check_netcdf( nf90_put_att(ncid, target_reason_mask_id, 'long_name', 'target reason bitmask') )
       call check_netcdf( nf90_put_att(ncid, phase_state_pre_id, 'long_name', &
            'pre-event phase state for phase-aware radius interpretation') )
       call check_netcdf( nf90_put_att(ncid, phase_state_post_id, 'long_name', &
            'post-event phase state for phase-aware radius interpretation') )
       call check_netcdf( nf90_put_att(ncid, x_id, 'units', 'm') )
       call check_netcdf( nf90_put_att(ncid, y_id, 'units', 'm') )
       call check_netcdf( nf90_put_att(ncid, z_id, 'units', 'm') )
       call check_netcdf( nf90_put_att(ncid, hydro_radius_pre_id, 'units', 'm') )
       call check_netcdf( nf90_put_att(ncid, hydro_radius_post_id, 'units', 'm') )
       call check_netcdf( nf90_put_att(ncid, liq_mass_pre_id, 'units', 'kg') )
       call check_netcdf( nf90_put_att(ncid, liq_mass_post_id, 'units', 'kg') )
       call check_netcdf( nf90_put_att(ncid, ice_mass_pre_id, 'units', 'kg') )
       call check_netcdf( nf90_put_att(ncid, ice_mass_post_id, 'units', 'kg') )
	       if( write_extended ) then
	          call check_netcdf( nf90_put_att(ncid, ice_re_pre_id, 'long_name', &
	               'ice equatorial radius before event') )
          call check_netcdf( nf90_put_att(ncid, ice_re_post_id, 'long_name', &
               'ice equatorial radius after event') )
          call check_netcdf( nf90_put_att(ncid, ice_rp_pre_id, 'long_name', &
               'ice polar radius before event') )
          call check_netcdf( nf90_put_att(ncid, ice_rp_post_id, 'long_name', &
               'ice polar radius after event') )
          call check_netcdf( nf90_put_att(ncid, ice_rho_pre_id, 'long_name', &
               'ice crystal density before event') )
          call check_netcdf( nf90_put_att(ncid, ice_rho_post_id, 'long_name', &
               'ice crystal density after event') )
          call check_netcdf( nf90_put_att(ncid, ice_re_pre_id, 'units', 'm') )
          call check_netcdf( nf90_put_att(ncid, ice_re_post_id, 'units', 'm') )
          call check_netcdf( nf90_put_att(ncid, ice_rp_pre_id, 'units', 'm') )
          call check_netcdf( nf90_put_att(ncid, ice_rp_post_id, 'units', 'm') )
          call check_netcdf( nf90_put_att(ncid, ice_rho_pre_id, 'units', 'kg/m3') )
	          call check_netcdf( nf90_put_att(ncid, ice_rho_post_id, 'units', 'kg/m3') )
	       end if
       if( write_kohler_context ) then
          call check_netcdf( nf90_put_att(ncid, kohler_rcrit_pre_id, 'long_name', &
               'Kohler critical radius before derived activation evaluation') )
          call check_netcdf( nf90_put_att(ncid, kohler_rcrit_post_id, 'long_name', &
               'Kohler critical radius after derived activation evaluation') )
          call check_netcdf( nf90_put_att(ncid, kohler_margin_pre_id, 'long_name', &
               'sd_r - Kohler critical radius before event') )
          call check_netcdf( nf90_put_att(ncid, kohler_margin_post_id, 'long_name', &
               'sd_r - Kohler critical radius after event') )
          call check_netcdf( nf90_put_att(ncid, activated_state_pre_id, 'long_name', &
               'derived Kohler activated state before event: 0=false, 1=true') )
          call check_netcdf( nf90_put_att(ncid, activated_state_post_id, 'long_name', &
               'derived Kohler activated state after event: 0=false, 1=true') )
          call check_netcdf( nf90_put_att(ncid, kohler_rcrit_pre_id, 'units', 'm') )
          call check_netcdf( nf90_put_att(ncid, kohler_rcrit_post_id, 'units', 'm') )
          call check_netcdf( nf90_put_att(ncid, kohler_margin_pre_id, 'units', 'm') )
          call check_netcdf( nf90_put_att(ncid, kohler_margin_post_id, 'units', 'm') )
       end if
       if( write_aerosol_context ) then
          call check_netcdf( nf90_put_att(ncid, aerosol_total_mass_pre_id, 'long_name', &
               'total aerosol mass before single-process event') )
          call check_netcdf( nf90_put_att(ncid, aerosol_total_mass_post_id, 'long_name', &
               'total aerosol mass after single-process event') )
          call check_netcdf( nf90_put_att(ncid, aerosol_kohler_solute_pre_id, 'long_name', &
               'Kohler solute sum before single-process event: sum(sd_asl*ion/molecular_weight)') )
          call check_netcdf( nf90_put_att(ncid, aerosol_kohler_solute_post_id, 'long_name', &
               'Kohler solute sum after single-process event: sum(sd_asl*ion/molecular_weight)') )
          call check_netcdf( nf90_put_att(ncid, aerosol_total_mass_pre_id, 'units', 'kg') )
          call check_netcdf( nf90_put_att(ncid, aerosol_total_mass_post_id, 'units', 'kg') )
          call check_netcdf( nf90_put_att(ncid, aerosol_kohler_solute_pre_id, 'units', 'model_units') )
          call check_netcdf( nf90_put_att(ncid, aerosol_kohler_solute_post_id, 'units', 'model_units') )
       end if
       if( write_thermo_context ) then
          call check_netcdf( nf90_put_att(ncid, air_temperature_id, 'long_name', &
               'local air temperature used by the single-process wrapper') )
          call check_netcdf( nf90_put_att(ncid, air_pressure_id, 'long_name', &
               'local air pressure used by the single-process wrapper') )
          call check_netcdf( nf90_put_att(ncid, water_vapor_mixing_ratio_id, 'long_name', &
               'local water-vapor mixing ratio used by the single-process wrapper') )
          call check_netcdf( nf90_put_att(ncid, air_temperature_id, 'units', 'K') )
          call check_netcdf( nf90_put_att(ncid, air_pressure_id, 'units', 'Pa') )
          call check_netcdf( nf90_put_att(ncid, water_vapor_mixing_ratio_id, 'units', 'kg/kg') )
       end if
       call check_netcdf( nf90_enddef(ncid) )
       event_offset = 0
    end if

    start1(1) = event_offset + 1
    count1(1) = rec_count
    call check_netcdf( nf90_put_var(ncid, time_id, time_out, start=start1, count=count1) )
    call check_netcdf( nf90_put_var(ncid, trigger_code_id, trigger_code_out, start=start1, count=count1) )
    call check_netcdf( nf90_put_var(ncid, trigger_level_id, trigger_level_out, start=start1, count=count1) )
    call check_netcdf( nf90_put_var(ncid, target_reason_mask_id, target_reason_mask_out, start=start1, count=count1) )
    call check_netcdf( nf90_put_var(ncid, phase_state_pre_id, phase_state_pre_out, start=start1, count=count1) )
    call check_netcdf( nf90_put_var(ncid, phase_state_post_id, phase_state_post_out, start=start1, count=count1) )
    if( write_tracking_ids ) then
       call check_netcdf( nf90_put_var(ncid, sd_id_id, sd_id_out, start=start1, count=count1) )
       call check_netcdf( nf90_put_var(ncid, dm_id_id, dm_id_out, start=start1, count=count1) )
    end if
    call check_netcdf( nf90_put_var(ncid, x_id, x_out, start=start1, count=count1) )
    call check_netcdf( nf90_put_var(ncid, y_id, y_out, start=start1, count=count1) )
    call check_netcdf( nf90_put_var(ncid, z_id, z_out, start=start1, count=count1) )
    call check_netcdf( nf90_put_var(ncid, hydro_radius_pre_id, hydro_radius_pre_out, start=start1, count=count1) )
    call check_netcdf( nf90_put_var(ncid, hydro_radius_post_id, hydro_radius_post_out, start=start1, count=count1) )
    call check_netcdf( nf90_put_var(ncid, liq_mass_pre_id, liq_mass_pre_out, start=start1, count=count1) )
    call check_netcdf( nf90_put_var(ncid, liq_mass_post_id, liq_mass_post_out, start=start1, count=count1) )
    call check_netcdf( nf90_put_var(ncid, ice_mass_pre_id, ice_mass_pre_out, start=start1, count=count1) )
    call check_netcdf( nf90_put_var(ncid, ice_mass_post_id, ice_mass_post_out, start=start1, count=count1) )
	    if( write_extended ) then
	       call check_netcdf( nf90_put_var(ncid, ice_re_pre_id, ice_re_pre_out, start=start1, count=count1) )
	       call check_netcdf( nf90_put_var(ncid, ice_re_post_id, ice_re_post_out, start=start1, count=count1) )
	       call check_netcdf( nf90_put_var(ncid, ice_rp_pre_id, ice_rp_pre_out, start=start1, count=count1) )
	       call check_netcdf( nf90_put_var(ncid, ice_rp_post_id, ice_rp_post_out, start=start1, count=count1) )
	       call check_netcdf( nf90_put_var(ncid, ice_rho_pre_id, ice_rho_pre_out, start=start1, count=count1) )
	       call check_netcdf( nf90_put_var(ncid, ice_rho_post_id, ice_rho_post_out, start=start1, count=count1) )
	    end if
    if( write_kohler_context ) then
       call check_netcdf( nf90_put_var(ncid, kohler_rcrit_pre_id, kohler_rcrit_pre_out, start=start1, count=count1) )
       call check_netcdf( nf90_put_var(ncid, kohler_rcrit_post_id, kohler_rcrit_post_out, start=start1, count=count1) )
       call check_netcdf( nf90_put_var(ncid, kohler_margin_pre_id, kohler_margin_pre_out, start=start1, count=count1) )
       call check_netcdf( nf90_put_var(ncid, kohler_margin_post_id, kohler_margin_post_out, start=start1, count=count1) )
       call check_netcdf( nf90_put_var(ncid, activated_state_pre_id, activated_state_pre_out, start=start1, count=count1) )
       call check_netcdf( nf90_put_var(ncid, activated_state_post_id, activated_state_post_out, start=start1, count=count1) )
    end if
    if( write_aerosol_context ) then
       call check_netcdf( nf90_put_var(ncid, aerosol_total_mass_pre_id, aerosol_total_mass_pre_out, start=start1, count=count1) )
       call check_netcdf( nf90_put_var(ncid, aerosol_total_mass_post_id, aerosol_total_mass_post_out, start=start1, count=count1) )
       call check_netcdf( nf90_put_var(ncid, aerosol_kohler_solute_pre_id, aerosol_kohler_solute_pre_out, &
            start=start1, count=count1) )
       call check_netcdf( nf90_put_var(ncid, aerosol_kohler_solute_post_id, aerosol_kohler_solute_post_out, &
            start=start1, count=count1) )
    end if
    if( write_thermo_context ) then
       call check_netcdf( nf90_put_var(ncid, air_temperature_id, air_temperature_out, start=start1, count=count1) )
       call check_netcdf( nf90_put_var(ncid, air_pressure_id, air_pressure_out, start=start1, count=count1) )
       call check_netcdf( nf90_put_var(ncid, water_vapor_mixing_ratio_id, water_vapor_mixing_ratio_out, &
            start=start1, count=count1) )
    end if
	    call check_netcdf( nf90_close(ncid) )

    if( allocated(sd_id_out) ) deallocate(sd_id_out)
    if( allocated(dm_id_out) ) deallocate(dm_id_out)
    deallocate(time_out, trigger_code_out, trigger_level_out, target_reason_mask_out, phase_state_pre_out, phase_state_post_out)
    deallocate(x_out, y_out, z_out)
    deallocate(hydro_radius_pre_out, hydro_radius_post_out)
    deallocate(liq_mass_pre_out, liq_mass_post_out, ice_mass_pre_out, ice_mass_post_out)
	    if( allocated(ice_re_pre_out) ) deallocate(ice_re_pre_out, ice_re_post_out)
	    if( allocated(ice_rp_pre_out) ) deallocate(ice_rp_pre_out, ice_rp_post_out)
	    if( allocated(ice_rho_pre_out) ) deallocate(ice_rho_pre_out, ice_rho_post_out)
    if( allocated(kohler_rcrit_pre_out) ) deallocate(kohler_rcrit_pre_out, kohler_rcrit_post_out)
    if( allocated(kohler_margin_pre_out) ) deallocate(kohler_margin_pre_out, kohler_margin_post_out)
    if( allocated(activated_state_pre_out) ) deallocate(activated_state_pre_out, activated_state_post_out)
    if( allocated(aerosol_total_mass_pre_out) ) deallocate(aerosol_total_mass_pre_out, aerosol_total_mass_post_out)
    if( allocated(aerosol_kohler_solute_pre_out) ) deallocate(aerosol_kohler_solute_pre_out, aerosol_kohler_solute_post_out)
    if( allocated(air_temperature_out) ) deallocate(air_temperature_out, air_pressure_out, water_vapor_mixing_ratio_out)

  contains
    logical function singleproc_writes_record(process_active, occurrence_enabled, significant_enabled, significant_hit) result(writes_record)
      logical, intent(in) :: process_active
      logical, intent(in) :: occurrence_enabled
      logical, intent(in) :: significant_enabled
      logical, intent(in) :: significant_hit

      writes_record = process_active .and. ( (significant_enabled .and. significant_hit) .or. occurrence_enabled )
    end function singleproc_writes_record

    subroutine maybe_put_singleproc_record(n, process_code, process_active, occurrence_enabled, &
         significant_enabled, significant_hit, phase_post)
      integer, intent(in) :: n, process_code, phase_post
      logical, intent(in) :: process_active
      logical, intent(in) :: occurrence_enabled
      logical, intent(in) :: significant_enabled
      logical, intent(in) :: significant_hit

      if( .not. process_active ) return
      if( significant_enabled .and. significant_hit ) then
         trigger_level = TRIG_LEVEL_SIGNIFICANT
      else if( occurrence_enabled ) then
         trigger_level = TRIG_LEVEL_OCCURRENCE
      else
         return
      end if
      rec = rec + 1
      trigger_code = process_code
      call put_singleproc_record(n, rec, trigger_code, trigger_level, phase_post)
    end subroutine maybe_put_singleproc_record

    subroutine put_singleproc_record(n, rec, trigger_code, trigger_level, phase_post)
      integer, intent(in) :: n, rec, trigger_code, trigger_level, phase_post

      time_out(rec) = otime
      trigger_code_out(rec) = trigger_code
      trigger_level_out(rec) = trigger_level
      ! event_multiplicity is collision-family only; single-process hits do not write it.
      target_reason_mask_out(rec) = TARGET_BY_EVENT
      phase_state_pre_out(rec) = phase_pre(n)
      phase_state_post_out(rec) = phase_post
      x_out(rec) = sd_x(n)
      y_out(rec) = sd_y(n)
      z_out(rec) = sd_z(n)
      hydro_radius_pre_out(rec) = hydro_radius_pre(n)
      hydro_radius_post_out(rec) = sdm_cold_hydro_radius( &
           phase_post, sd_r(n), sdi%re(n), sdi%rp(n))
      liq_mass_pre_out(rec) = liq_mass_pre(n)
      liq_mass_post_out(rec) = sdm_cold_liq_mass(sd_r(n))
      ice_mass_pre_out(rec) = ice_mass_pre(n)
      ice_mass_post_out(rec) = sdm_cold_ice_mass(sdi%re(n), sdi%rp(n), sdi%rho(n))
	      if( write_extended ) then
	         ice_re_pre_out(rec) = ice_re_pre(n)
	         ice_re_post_out(rec) = sdi%re(n)
	         ice_rp_pre_out(rec) = ice_rp_pre(n)
	         ice_rp_post_out(rec) = sdi%rp(n)
	         ice_rho_pre_out(rec) = ice_rho_pre(n)
	         ice_rho_post_out(rec) = sdi%rho(n)
	      end if
      if( write_kohler_context ) then
         if( has_kohler_context ) then
            kohler_rcrit_pre_out(rec) = kohler_rcrit_pre(n)
            kohler_rcrit_post_out(rec) = kohler_rcrit_post(n)
            kohler_margin_pre_out(rec) = kohler_margin_pre(n)
            kohler_margin_post_out(rec) = kohler_margin_post(n)
            activated_state_pre_out(rec) = kohler_active_pre(n)
            activated_state_post_out(rec) = kohler_active_post(n)
         else
            kohler_rcrit_pre_out(rec) = -huge(1.0_RP)
            kohler_rcrit_post_out(rec) = -huge(1.0_RP)
            kohler_margin_pre_out(rec) = -huge(1.0_RP)
            kohler_margin_post_out(rec) = -huge(1.0_RP)
            activated_state_pre_out(rec) = 0
            activated_state_post_out(rec) = 0
         end if
      end if
      if( write_aerosol_context ) then
         if( has_aerosol_context ) then
            aerosol_total_mass_pre_out(rec) = aerosol_total_mass_pre(n)
            aerosol_total_mass_post_out(rec) = aerosol_total_mass_post(n)
            aerosol_kohler_solute_pre_out(rec) = aerosol_kohler_solute_pre(n)
            aerosol_kohler_solute_post_out(rec) = aerosol_kohler_solute_post(n)
         else
            aerosol_total_mass_pre_out(rec) = -huge(1.0_RP)
            aerosol_total_mass_post_out(rec) = -huge(1.0_RP)
            aerosol_kohler_solute_pre_out(rec) = -huge(1.0_RP)
            aerosol_kohler_solute_post_out(rec) = -huge(1.0_RP)
         end if
      end if
      if( write_thermo_context ) then
         if( has_thermo_context ) then
            air_temperature_out(rec) = air_temperature(n)
            air_pressure_out(rec) = air_pressure(n)
            water_vapor_mixing_ratio_out(rec) = water_vapor_mixing_ratio(n)
         else
            air_temperature_out(rec) = -huge(1.0_RP)
            air_pressure_out(rec) = -huge(1.0_RP)
            water_vapor_mixing_ratio_out(rec) = -huge(1.0_RP)
         end if
      end if
	      if( write_tracking_ids ) then
	         sd_id_out(rec) = sd_id(n)
	         dm_id_out(rec) = dm_id(n)
      end if
    end subroutine put_singleproc_record
  end subroutine sdm_event_singleproc_outnetcdf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sdm_event_diag_outnetcdf(otime, sd_num, sd_id, dm_id, sd_x, sd_y, sd_z, &
       &                              sd_event_mask, sd_event_sig_mask, sd_diag_mask, sd_phase_change_flag, &
       &                              sd_liq_radius_max_interval, sd_ice_rvol_max_interval, &
       &                              sd_mixed_rvol_max_interval, sd_rime_mass_max_interval, &
       &                              sd_rime_frac_max_interval, sd_nmono_max_interval, &
       &                              sd_aspect_ratio_max_interval)
    use netcdf
    use scale_precision
    use scale_stdio
    use scale_process, only: &
         mype => PRC_myrank
    use m_sdm_common, only: &
         sdm_cold, sdm_dmpitvl, tracking_mode
    use m_sdm_tracking_cold, only: &
         DIAG_LIQ_RADIUS_LARGE, DIAG_ICE_RVOL_LARGE, DIAG_MIXED_RVOL_LARGE, &
         DIAG_RIME_MASS_LARGE, DIAG_RIME_FRAC_LARGE, DIAG_NMONO_LARGE, &
         DIAG_ASPECT_RATIO_MATCH, TARGET_BY_DIAG, sdm_cold_diag_trigger_from_bit

    implicit none

    real(DP), intent(in) :: otime
    integer, intent(in) :: sd_num
    integer, intent(in) :: sd_id(1:sd_num)
    integer, intent(in) :: dm_id(1:sd_num)
    real(RP), intent(in) :: sd_x(1:sd_num)
    real(RP), intent(in) :: sd_y(1:sd_num)
    real(RP), intent(in) :: sd_z(1:sd_num)
    integer, intent(in) :: sd_event_mask(1:sd_num)
    integer, intent(in) :: sd_event_sig_mask(1:sd_num)
    integer, intent(in) :: sd_diag_mask(1:sd_num)
    integer, intent(in) :: sd_phase_change_flag(1:sd_num)
    real(RP), intent(in) :: sd_liq_radius_max_interval(1:sd_num)
    real(RP), intent(in) :: sd_ice_rvol_max_interval(1:sd_num)
    real(RP), intent(in) :: sd_mixed_rvol_max_interval(1:sd_num)
    real(RP), intent(in) :: sd_rime_mass_max_interval(1:sd_num)
    real(RP), intent(in) :: sd_rime_frac_max_interval(1:sd_num)
    real(RP), intent(in) :: sd_nmono_max_interval(1:sd_num)
    real(RP), intent(in) :: sd_aspect_ratio_max_interval(1:sd_num)

    integer, parameter :: ndiag = 7
    integer, parameter :: diag_bits(ndiag) = (/ DIAG_LIQ_RADIUS_LARGE, &
         DIAG_ICE_RVOL_LARGE, DIAG_MIXED_RVOL_LARGE, DIAG_RIME_MASS_LARGE, &
         DIAG_RIME_FRAC_LARGE, DIAG_NMONO_LARGE, DIAG_ASPECT_RATIO_MATCH /)
    integer, parameter :: nc_deflate_level = 1
    integer, parameter :: nc_deflate       = 1
    integer, parameter :: nc_shuffle       = 1

    character(len=17) :: fmt2="(A, '.', A, I*.*)"
    character(len=H_LONG) :: ftmp
    character(len=H_LONG) :: basename_sd_out
    character(len=19) :: basename_time
    integer :: nf90_real_precision
    integer :: ncid, event_dim_id, time_id, trigger_code_id, target_reason_mask_id
    integer :: sd_id_id, dm_id_id, x_id, y_id, z_id
    integer :: sd_event_mask_id, sd_event_sig_mask_id, sd_diag_mask_id, sd_phase_change_flag_id
    integer :: liq_radius_id, ice_rvol_id, mixed_rvol_id, rime_mass_id
    integer :: rime_frac_id, nmono_id, aspect_ratio_id
    integer :: event_offset, start1(1), count1(1)
    integer :: n, b, rec_count, rec
    real(DP) :: otime_bucket, otime_bucket_daysec
    integer :: otime_bucket_hh, otime_bucket_mm, otime_bucket_ss, otime_bucket_ms
    logical :: file_exists, write_tracking_ids
    character(len=16) :: id_name, dm_name
    real(DP), allocatable :: time_out(:)
    integer, allocatable :: trigger_code_out(:), target_reason_mask_out(:)
    integer, allocatable :: sd_id_out(:), dm_id_out(:)
    integer, allocatable :: sd_event_mask_out(:), sd_event_sig_mask_out(:), sd_diag_mask_out(:)
    integer, allocatable :: sd_phase_change_flag_out(:)
    real(RP), allocatable :: x_out(:), y_out(:), z_out(:)
    real(RP), allocatable :: liq_radius_out(:), ice_rvol_out(:), mixed_rvol_out(:)
    real(RP), allocatable :: rime_mass_out(:), rime_frac_out(:), nmono_out(:), aspect_ratio_out(:)

    if( .not. sdm_cold ) return

    if(RP == SP)then
       nf90_real_precision = NF90_FLOAT
    else
       nf90_real_precision = NF90_DOUBLE
    end if

    rec_count = 0
    do n = 1, sd_num
       do b = 1, ndiag
          if( iand(sd_diag_mask(n), diag_bits(b)) /= 0 ) rec_count = rec_count + 1
       end do
    end do
    if( rec_count <= 0 ) return

    allocate(time_out(rec_count), trigger_code_out(rec_count), target_reason_mask_out(rec_count))
    allocate(x_out(rec_count), y_out(rec_count), z_out(rec_count))
    allocate(sd_event_mask_out(rec_count), sd_event_sig_mask_out(rec_count), sd_diag_mask_out(rec_count))
    allocate(sd_phase_change_flag_out(rec_count))
    allocate(liq_radius_out(rec_count), ice_rvol_out(rec_count), mixed_rvol_out(rec_count))
    allocate(rime_mass_out(rec_count), rime_frac_out(rec_count), nmono_out(rec_count), aspect_ratio_out(rec_count))

    write_tracking_ids = tracking_mode /= 0
    if( write_tracking_ids ) then
       allocate(sd_id_out(rec_count))
       allocate(dm_id_out(rec_count))
    end if

    rec = 0
    do n = 1, sd_num
       do b = 1, ndiag
          if( iand(sd_diag_mask(n), diag_bits(b)) == 0 ) cycle
          rec = rec + 1
          time_out(rec) = otime
          trigger_code_out(rec) = sdm_cold_diag_trigger_from_bit(diag_bits(b))
          target_reason_mask_out(rec) = TARGET_BY_DIAG
          ! Diagnostic coordinates are SD output-time positions, not maxima-time positions.
          x_out(rec) = sd_x(n)
          y_out(rec) = sd_y(n)
          z_out(rec) = sd_z(n)
          sd_event_mask_out(rec) = sd_event_mask(n)
          sd_event_sig_mask_out(rec) = sd_event_sig_mask(n)
          sd_diag_mask_out(rec) = sd_diag_mask(n)
          sd_phase_change_flag_out(rec) = sd_phase_change_flag(n)
          liq_radius_out(rec) = sd_liq_radius_max_interval(n)
          ice_rvol_out(rec) = sd_ice_rvol_max_interval(n)
          mixed_rvol_out(rec) = sd_mixed_rvol_max_interval(n)
          rime_mass_out(rec) = sd_rime_mass_max_interval(n)
          rime_frac_out(rec) = sd_rime_frac_max_interval(n)
          nmono_out(rec) = sd_nmono_max_interval(n)
          aspect_ratio_out(rec) = sd_aspect_ratio_max_interval(n)
          if( write_tracking_ids ) then
             sd_id_out(rec) = sd_id(n)
             dm_id_out(rec) = dm_id(n)
          end if
       end do
    end do

    if( tracking_mode == 2 ) then
       id_name = 'pre_sdid'
       dm_name = 'pre_dmid'
    else if( tracking_mode == 1 ) then
       id_name = 'sd_id'
       dm_name = 'dm_id'
    else
       id_name = ''
       dm_name = ''
    end if

    if( sdm_dmpitvl > 0.0_RP ) then
       otime_bucket = sdm_dmpitvl * real( int( otime / sdm_dmpitvl + 1.0e-10_RP ), kind=DP )
    else
       otime_bucket = otime
    end if
    otime_bucket_daysec = mod( otime_bucket, 86400.0_DP )
    otime_bucket_hh = int( otime_bucket_daysec / 3600.0_DP )
    otime_bucket_mm = int( mod(otime_bucket_daysec, 3600.0_DP) / 60.0_DP )
    otime_bucket_ss = int( mod(otime_bucket_daysec, 60.0_DP) )
    otime_bucket_ms = nint( (otime_bucket_daysec - real(int(otime_bucket_daysec),kind=DP)) * 1000.0_DP )
    if( otime_bucket_ms >= 1000 ) then
       otime_bucket_ms = otime_bucket_ms - 1000
       otime_bucket_ss = otime_bucket_ss + 1
    end if
    write(basename_time,'(I8.8,A1,I2.2,I2.2,I2.2,A1,I3.3)') 101, '-', &
         otime_bucket_hh, otime_bucket_mm, otime_bucket_ss, '.', otime_bucket_ms
    write(fmt2(14:14),'(I1)') 6
    write(fmt2(16:16),'(I1)') 6
    write(ftmp,'(3A)') 'SD_event_diag', '_NetCDF_', trim(basename_time)
    write(basename_sd_out,fmt2) trim(ftmp), 'pe',mype

    inquire(file=trim(basename_sd_out), exist=file_exists)
    if( file_exists ) then
       call check_netcdf( nf90_open(trim(basename_sd_out), NF90_WRITE, ncid) )
       call check_netcdf( nf90_inq_dimid(ncid, "event", event_dim_id) )
       call check_netcdf( nf90_inquire_dimension(ncid, event_dim_id, len=event_offset) )
       call check_netcdf( nf90_inq_varid(ncid, "time", time_id) )
       call check_netcdf( nf90_inq_varid(ncid, "trigger_code", trigger_code_id) )
       call check_netcdf( nf90_inq_varid(ncid, "target_reason_mask", target_reason_mask_id) )
       if( write_tracking_ids ) then
          call check_netcdf( nf90_inq_varid(ncid, trim(id_name), sd_id_id) )
          call check_netcdf( nf90_inq_varid(ncid, trim(dm_name), dm_id_id) )
       end if
       call check_netcdf( nf90_inq_varid(ncid, "x", x_id) )
       call check_netcdf( nf90_inq_varid(ncid, "y", y_id) )
       call check_netcdf( nf90_inq_varid(ncid, "z", z_id) )
       call check_netcdf( nf90_inq_varid(ncid, "sd_event_mask", sd_event_mask_id) )
       call check_netcdf( nf90_inq_varid(ncid, "sd_event_sig_mask", sd_event_sig_mask_id) )
       call check_netcdf( nf90_inq_varid(ncid, "sd_diag_mask", sd_diag_mask_id) )
       call check_netcdf( nf90_inq_varid(ncid, "sd_phase_change_flag", sd_phase_change_flag_id) )
       call check_netcdf( nf90_inq_varid(ncid, "sd_liq_radius_max_interval", liq_radius_id) )
       call check_netcdf( nf90_inq_varid(ncid, "sd_ice_rvol_max_interval", ice_rvol_id) )
       call check_netcdf( nf90_inq_varid(ncid, "sd_mixed_rvol_max_interval", mixed_rvol_id) )
       call check_netcdf( nf90_inq_varid(ncid, "sd_rime_mass_max_interval", rime_mass_id) )
       call check_netcdf( nf90_inq_varid(ncid, "sd_rime_frac_max_interval", rime_frac_id) )
       call check_netcdf( nf90_inq_varid(ncid, "sd_nmono_max_interval", nmono_id) )
       call check_netcdf( nf90_inq_varid(ncid, "sd_aspect_ratio_max_interval", aspect_ratio_id) )
    else
       call check_netcdf( nf90_create(trim(basename_sd_out), NF90_NETCDF4, ncid) )
       call check_netcdf( nf90_def_dim(ncid, "event", NF90_UNLIMITED, event_dim_id) )
       call check_netcdf( nf90_def_var(ncid, "time", NF90_DOUBLE, event_dim_id, time_id) )
       call check_netcdf( nf90_def_var(ncid, "trigger_code", NF90_INT, event_dim_id, trigger_code_id) )
       call check_netcdf( nf90_def_var(ncid, "target_reason_mask", NF90_INT, event_dim_id, target_reason_mask_id) )
       if( write_tracking_ids ) then
          call check_netcdf( nf90_def_var(ncid, trim(id_name), NF90_INT, event_dim_id, sd_id_id) )
          call check_netcdf( nf90_def_var(ncid, trim(dm_name), NF90_INT, event_dim_id, dm_id_id) )
       end if
       call check_netcdf( nf90_def_var(ncid, "x", nf90_real_precision, event_dim_id, x_id) )
       call check_netcdf( nf90_def_var(ncid, "y", nf90_real_precision, event_dim_id, y_id) )
       call check_netcdf( nf90_def_var(ncid, "z", nf90_real_precision, event_dim_id, z_id) )
       call check_netcdf( nf90_def_var(ncid, "sd_event_mask", NF90_INT, event_dim_id, sd_event_mask_id) )
       call check_netcdf( nf90_def_var(ncid, "sd_event_sig_mask", NF90_INT, event_dim_id, sd_event_sig_mask_id) )
       call check_netcdf( nf90_def_var(ncid, "sd_diag_mask", NF90_INT, event_dim_id, sd_diag_mask_id) )
       call check_netcdf( nf90_def_var(ncid, "sd_phase_change_flag", NF90_INT, event_dim_id, sd_phase_change_flag_id) )
       call check_netcdf( nf90_def_var(ncid, "sd_liq_radius_max_interval", nf90_real_precision, event_dim_id, liq_radius_id) )
       call check_netcdf( nf90_def_var(ncid, "sd_ice_rvol_max_interval", nf90_real_precision, event_dim_id, ice_rvol_id) )
       call check_netcdf( nf90_def_var(ncid, "sd_mixed_rvol_max_interval", nf90_real_precision, event_dim_id, mixed_rvol_id) )
       call check_netcdf( nf90_def_var(ncid, "sd_rime_mass_max_interval", nf90_real_precision, event_dim_id, rime_mass_id) )
       call check_netcdf( nf90_def_var(ncid, "sd_rime_frac_max_interval", nf90_real_precision, event_dim_id, rime_frac_id) )
       call check_netcdf( nf90_def_var(ncid, "sd_nmono_max_interval", nf90_real_precision, event_dim_id, nmono_id) )
       call check_netcdf( nf90_def_var(ncid, "sd_aspect_ratio_max_interval", nf90_real_precision, event_dim_id, aspect_ratio_id) )
       call check_netcdf( nf90_put_att(ncid, trigger_code_id, 'long_name', 'cold-SDM diagnostic trigger code') )
       call check_netcdf( nf90_put_att(ncid, target_reason_mask_id, 'long_name', 'target reason bitmask') )
       call check_netcdf( nf90_put_att(ncid, x_id, 'units', 'm') )
       call check_netcdf( nf90_put_att(ncid, y_id, 'units', 'm') )
       call check_netcdf( nf90_put_att(ncid, z_id, 'units', 'm') )
       call check_netcdf( nf90_enddef(ncid) )
       event_offset = 0
    end if

    start1(1) = event_offset + 1
    count1(1) = rec_count
    call check_netcdf( nf90_put_var(ncid, time_id, time_out, start=start1, count=count1) )
    call check_netcdf( nf90_put_var(ncid, trigger_code_id, trigger_code_out, start=start1, count=count1) )
    call check_netcdf( nf90_put_var(ncid, target_reason_mask_id, target_reason_mask_out, start=start1, count=count1) )
    if( write_tracking_ids ) then
       call check_netcdf( nf90_put_var(ncid, sd_id_id, sd_id_out, start=start1, count=count1) )
       call check_netcdf( nf90_put_var(ncid, dm_id_id, dm_id_out, start=start1, count=count1) )
    end if
    call check_netcdf( nf90_put_var(ncid, x_id, x_out, start=start1, count=count1) )
    call check_netcdf( nf90_put_var(ncid, y_id, y_out, start=start1, count=count1) )
    call check_netcdf( nf90_put_var(ncid, z_id, z_out, start=start1, count=count1) )
    call check_netcdf( nf90_put_var(ncid, sd_event_mask_id, sd_event_mask_out, start=start1, count=count1) )
    call check_netcdf( nf90_put_var(ncid, sd_event_sig_mask_id, sd_event_sig_mask_out, start=start1, count=count1) )
    call check_netcdf( nf90_put_var(ncid, sd_diag_mask_id, sd_diag_mask_out, start=start1, count=count1) )
    call check_netcdf( nf90_put_var(ncid, sd_phase_change_flag_id, sd_phase_change_flag_out, start=start1, count=count1) )
    call check_netcdf( nf90_put_var(ncid, liq_radius_id, liq_radius_out, start=start1, count=count1) )
    call check_netcdf( nf90_put_var(ncid, ice_rvol_id, ice_rvol_out, start=start1, count=count1) )
    call check_netcdf( nf90_put_var(ncid, mixed_rvol_id, mixed_rvol_out, start=start1, count=count1) )
    call check_netcdf( nf90_put_var(ncid, rime_mass_id, rime_mass_out, start=start1, count=count1) )
    call check_netcdf( nf90_put_var(ncid, rime_frac_id, rime_frac_out, start=start1, count=count1) )
    call check_netcdf( nf90_put_var(ncid, nmono_id, nmono_out, start=start1, count=count1) )
    call check_netcdf( nf90_put_var(ncid, aspect_ratio_id, aspect_ratio_out, start=start1, count=count1) )
    call check_netcdf( nf90_close(ncid) )

    if( allocated(sd_id_out) ) deallocate(sd_id_out)
    if( allocated(dm_id_out) ) deallocate(dm_id_out)
    deallocate(time_out, trigger_code_out, target_reason_mask_out)
    deallocate(x_out, y_out, z_out)
    deallocate(sd_event_mask_out, sd_event_sig_mask_out, sd_diag_mask_out, sd_phase_change_flag_out)
    deallocate(liq_radius_out, ice_rvol_out, mixed_rvol_out)
    deallocate(rime_mass_out, rime_frac_out, nmono_out, aspect_ratio_out)

  end subroutine sdm_event_diag_outnetcdf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sdm_lifecycle_outnetcdf(otime, lifecycle_code, n_parent, n_child, &
       parent_sd_id, parent_dm_id, child_sd_id, child_dm_id, x, y, z)
    use netcdf
    use scale_precision
    use scale_stdio
    use scale_process, only: &
         mype => PRC_myrank
    use m_sdm_common, only: &
         sdm_cold, sdm_dmpitvl, tracking_mode, TRACK_ID_INVALID

    implicit none

    real(DP), intent(in) :: otime
    integer, intent(in) :: lifecycle_code
    integer, intent(in) :: n_parent
    integer, intent(in) :: n_child
    integer, intent(in) :: parent_sd_id(:)
    integer, intent(in) :: parent_dm_id(:)
    integer, intent(in) :: child_sd_id(:)
    integer, intent(in) :: child_dm_id(:)
    real(RP), intent(in) :: x
    real(RP), intent(in) :: y
    real(RP), intent(in) :: z

    character(len=17) :: fmt2="(A, '.', A, I*.*)"
    character(len=H_LONG) :: ftmp
    character(len=H_LONG) :: basename_sd_out
    character(len=19) :: basename_time
    integer :: nf90_real_precision
    integer :: ncid, event_dim_id, event_offset
    integer :: time_id, lifecycle_code_id, n_parent_id, n_child_id
    integer :: parent_sd_id1_id, parent_dm_id1_id, parent_sd_id2_id, parent_dm_id2_id
    integer :: child_sd_id1_id, child_dm_id1_id, child_sd_id2_id, child_dm_id2_id
    integer :: child_sd_id3_id, child_dm_id3_id, child_sd_id4_id, child_dm_id4_id
    integer :: x_id, y_id, z_id
    integer :: start1(1), count1(1)
    integer :: otime_bucket_hh, otime_bucket_mm, otime_bucket_ss, otime_bucket_ms
    integer :: parent_sd(2), parent_dm(2), child_sd(4), child_dm(4)
    integer :: n, ncopy
    real(DP) :: otime_bucket, otime_bucket_daysec
    logical :: file_exists

    if( .not. sdm_cold ) return
    if( tracking_mode == 0 ) return

    if(RP == SP)then
       nf90_real_precision = NF90_FLOAT
    else
       nf90_real_precision = NF90_DOUBLE
    end if

    parent_sd(:) = TRACK_ID_INVALID
    parent_dm(:) = TRACK_ID_INVALID
    child_sd(:) = TRACK_ID_INVALID
    child_dm(:) = TRACK_ID_INVALID
    ncopy = min(2, max(0, n_parent), size(parent_sd_id), size(parent_dm_id))
    do n = 1, ncopy
       parent_sd(n) = parent_sd_id(n)
       parent_dm(n) = parent_dm_id(n)
    end do
    ncopy = min(4, max(0, n_child), size(child_sd_id), size(child_dm_id))
    do n = 1, ncopy
       child_sd(n) = child_sd_id(n)
       child_dm(n) = child_dm_id(n)
    end do

    if( sdm_dmpitvl > 0.0_RP ) then
       otime_bucket = sdm_dmpitvl * real( int( otime / sdm_dmpitvl + 1.0e-10_RP ), kind=DP )
    else
       otime_bucket = otime
    end if
    otime_bucket_daysec = mod( otime_bucket, 86400.0_DP )
    otime_bucket_hh = int( otime_bucket_daysec / 3600.0_DP )
    otime_bucket_mm = int( mod(otime_bucket_daysec, 3600.0_DP) / 60.0_DP )
    otime_bucket_ss = int( mod(otime_bucket_daysec, 60.0_DP) )
    otime_bucket_ms = nint( (otime_bucket_daysec - real(int(otime_bucket_daysec),kind=DP)) * 1000.0_DP )
    if( otime_bucket_ms >= 1000 ) then
       otime_bucket_ms = otime_bucket_ms - 1000
       otime_bucket_ss = otime_bucket_ss + 1
       if( otime_bucket_ss >= 60 ) then
          otime_bucket_ss = otime_bucket_ss - 60
          otime_bucket_mm = otime_bucket_mm + 1
          if( otime_bucket_mm >= 60 ) then
             otime_bucket_mm = otime_bucket_mm - 60
             otime_bucket_hh = otime_bucket_hh + 1
          end if
       end if
    end if
    write(basename_time,'(I8.8,A1,I2.2,I2.2,I2.2,A1,I3.3)') 101, '-', &
         otime_bucket_hh, otime_bucket_mm, otime_bucket_ss, '.', otime_bucket_ms
    write(fmt2(14:14),'(I1)') 6
    write(fmt2(16:16),'(I1)') 6
    write(ftmp,'(3A)') 'SD_lifecycle', '_NetCDF_', trim(basename_time)
    write(basename_sd_out,fmt2) trim(ftmp), 'pe',mype

    inquire(file=trim(basename_sd_out), exist=file_exists)
    if( file_exists ) then
       call check_netcdf( nf90_open(trim(basename_sd_out), NF90_WRITE, ncid) )
       call check_netcdf( nf90_inq_dimid(ncid, "event", event_dim_id) )
       call check_netcdf( nf90_inquire_dimension(ncid, event_dim_id, len=event_offset) )
       call inq_var("time", time_id)
       call inq_var("lifecycle_code", lifecycle_code_id)
       call inq_var("n_parent", n_parent_id)
       call inq_var("n_child", n_child_id)
       call inq_var("parent_sd_id1", parent_sd_id1_id)
       call inq_var("parent_dm_id1", parent_dm_id1_id)
       call inq_var("parent_sd_id2", parent_sd_id2_id)
       call inq_var("parent_dm_id2", parent_dm_id2_id)
       call inq_var("child_sd_id1", child_sd_id1_id)
       call inq_var("child_dm_id1", child_dm_id1_id)
       call inq_var("child_sd_id2", child_sd_id2_id)
       call inq_var("child_dm_id2", child_dm_id2_id)
       call inq_var("child_sd_id3", child_sd_id3_id)
       call inq_var("child_dm_id3", child_dm_id3_id)
       call inq_var("child_sd_id4", child_sd_id4_id)
       call inq_var("child_dm_id4", child_dm_id4_id)
       call inq_var("x", x_id)
       call inq_var("y", y_id)
       call inq_var("z", z_id)
    else
       call check_netcdf( nf90_create(trim(basename_sd_out), NF90_NETCDF4, ncid) )
       call check_netcdf( nf90_def_dim(ncid, "event", NF90_UNLIMITED, event_dim_id) )
       call def_time("time", "lifecycle record time", "s", time_id)
       call def_int("lifecycle_code", "SD lifecycle code", "", lifecycle_code_id)
       call def_int("n_parent", "number of parent IDs in this fixed-width record", "", n_parent_id)
       call def_int("n_child", "number of child IDs in this fixed-width record", "", n_child_id)
       call def_int("parent_sd_id1", "parent-1 SD ID", "", parent_sd_id1_id)
       call def_int("parent_dm_id1", "parent-1 domain ID", "", parent_dm_id1_id)
       call def_int("parent_sd_id2", "parent-2 SD ID", "", parent_sd_id2_id)
       call def_int("parent_dm_id2", "parent-2 domain ID", "", parent_dm_id2_id)
       call def_int("child_sd_id1", "child-1 SD ID", "", child_sd_id1_id)
       call def_int("child_dm_id1", "child-1 domain ID", "", child_dm_id1_id)
       call def_int("child_sd_id2", "child-2 SD ID", "", child_sd_id2_id)
       call def_int("child_dm_id2", "child-2 domain ID", "", child_dm_id2_id)
       call def_int("child_sd_id3", "child-3 SD ID", "", child_sd_id3_id)
       call def_int("child_dm_id3", "child-3 domain ID", "", child_dm_id3_id)
       call def_int("child_sd_id4", "child-4 SD ID", "", child_sd_id4_id)
       call def_int("child_dm_id4", "child-4 domain ID", "", child_dm_id4_id)
       call def_real("x", "lifecycle record x", "m", x_id)
       call def_real("y", "lifecycle record y", "m", y_id)
       call def_real("z", "lifecycle record z", "m", z_id)
       call check_netcdf( nf90_enddef(ncid) )
       event_offset = 0
    end if

    start1(1) = event_offset + 1
    count1(1) = 1
    call put_real_dp(time_id, otime)
    call put_int_scalar(lifecycle_code_id, lifecycle_code)
    call put_int_scalar(n_parent_id, n_parent)
    call put_int_scalar(n_child_id, n_child)
    call put_int_scalar(parent_sd_id1_id, parent_sd(1))
    call put_int_scalar(parent_dm_id1_id, parent_dm(1))
    call put_int_scalar(parent_sd_id2_id, parent_sd(2))
    call put_int_scalar(parent_dm_id2_id, parent_dm(2))
    call put_int_scalar(child_sd_id1_id, child_sd(1))
    call put_int_scalar(child_dm_id1_id, child_dm(1))
    call put_int_scalar(child_sd_id2_id, child_sd(2))
    call put_int_scalar(child_dm_id2_id, child_dm(2))
    call put_int_scalar(child_sd_id3_id, child_sd(3))
    call put_int_scalar(child_dm_id3_id, child_dm(3))
    call put_int_scalar(child_sd_id4_id, child_sd(4))
    call put_int_scalar(child_dm_id4_id, child_dm(4))
    call put_real_rp(x_id, x)
    call put_real_rp(y_id, y)
    call put_real_rp(z_id, z)
    call check_netcdf( nf90_close(ncid) )

  contains
    subroutine def_int(name, long_name, units, varid)
      character(len=*), intent(in) :: name, long_name, units
      integer, intent(out) :: varid
      call check_netcdf( nf90_def_var(ncid, trim(name), NF90_INT, event_dim_id, varid) )
      call check_netcdf( nf90_put_att(ncid, varid, 'long_name', trim(long_name)) )
      call check_netcdf( nf90_put_att(ncid, varid, 'units', trim(units)) )
    end subroutine def_int

    subroutine def_real(name, long_name, units, varid)
      character(len=*), intent(in) :: name, long_name, units
      integer, intent(out) :: varid
      call check_netcdf( nf90_def_var(ncid, trim(name), nf90_real_precision, event_dim_id, varid) )
      call check_netcdf( nf90_put_att(ncid, varid, 'long_name', trim(long_name)) )
      call check_netcdf( nf90_put_att(ncid, varid, 'units', trim(units)) )
    end subroutine def_real

    subroutine def_time(name, long_name, units, varid)
      character(len=*), intent(in) :: name, long_name, units
      integer, intent(out) :: varid
      call check_netcdf( nf90_def_var(ncid, trim(name), NF90_DOUBLE, event_dim_id, varid) )
      call check_netcdf( nf90_put_att(ncid, varid, 'long_name', trim(long_name)) )
      call check_netcdf( nf90_put_att(ncid, varid, 'units', trim(units)) )
    end subroutine def_time

    subroutine inq_var(name, varid)
      character(len=*), intent(in) :: name
      integer, intent(out) :: varid
      call check_netcdf( nf90_inq_varid(ncid, trim(name), varid) )
    end subroutine inq_var

    subroutine put_int_scalar(varid, value)
      integer, intent(in) :: varid
      integer, intent(in) :: value
      call check_netcdf( nf90_put_var(ncid, varid, (/value/), start=start1, count=count1) )
    end subroutine put_int_scalar

    subroutine put_real_rp(varid, value)
      integer, intent(in) :: varid
      real(RP), intent(in) :: value
      call check_netcdf( nf90_put_var(ncid, varid, (/value/), start=start1, count=count1) )
    end subroutine put_real_rp

    subroutine put_real_dp(varid, value)
      integer, intent(in) :: varid
      real(DP), intent(in) :: value
      call check_netcdf( nf90_put_var(ncid, varid, (/value/), start=start1, count=count1) )
    end subroutine put_real_dp
  end subroutine sdm_lifecycle_outnetcdf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sdm_tracking_lifecycle_write_domain_entry(otime, child_sd_id, child_dm_id, x, y, z)
    use scale_precision
    use m_sdm_common, only: &
         TRACK_ID_INVALID
    use m_sdm_tracking_cold, only: &
         LIFE_DOMAIN_ENTRY

    real(DP), intent(in) :: otime
    integer, intent(in) :: child_sd_id
    integer, intent(in) :: child_dm_id
    real(RP), intent(in) :: x
    real(RP), intent(in) :: y
    real(RP), intent(in) :: z
    integer :: parent_sd(1), parent_dm(1), child_sd(1), child_dm(1)

    parent_sd(1) = TRACK_ID_INVALID
    parent_dm(1) = TRACK_ID_INVALID
    child_sd(1) = child_sd_id
    child_dm(1) = child_dm_id
    call sdm_lifecycle_outnetcdf(otime, LIFE_DOMAIN_ENTRY, 0, 1, &
         parent_sd, parent_dm, child_sd, child_dm, x, y, z)
  end subroutine sdm_tracking_lifecycle_write_domain_entry
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sdm_tracking_lifecycle_write_global_halo_entry(otime, child_sd_id, child_dm_id, x, y, z)
    use scale_precision
    use m_sdm_common, only: &
         TRACK_ID_INVALID
    use m_sdm_tracking_cold, only: &
         LIFE_GLOBAL_HALO_ENTRY

    real(DP), intent(in) :: otime
    integer, intent(in) :: child_sd_id
    integer, intent(in) :: child_dm_id
    real(RP), intent(in) :: x
    real(RP), intent(in) :: y
    real(RP), intent(in) :: z
    integer :: parent_sd(1), parent_dm(1), child_sd(1), child_dm(1)

    parent_sd(1) = TRACK_ID_INVALID
    parent_dm(1) = TRACK_ID_INVALID
    child_sd(1) = child_sd_id
    child_dm(1) = child_dm_id
    call sdm_lifecycle_outnetcdf(otime, LIFE_GLOBAL_HALO_ENTRY, 0, 1, &
         parent_sd, parent_dm, child_sd, child_dm, x, y, z)
  end subroutine sdm_tracking_lifecycle_write_global_halo_entry
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sdm_tracking_lifecycle_write_seeding_entry(otime, child_sd_id, child_dm_id, x, y, z)
    use scale_precision
    use m_sdm_common, only: &
         TRACK_ID_INVALID
    use m_sdm_tracking_cold, only: &
         LIFE_SEEDING_ENTRY

    real(DP), intent(in) :: otime
    integer, intent(in) :: child_sd_id
    integer, intent(in) :: child_dm_id
    real(RP), intent(in) :: x
    real(RP), intent(in) :: y
    real(RP), intent(in) :: z
    integer :: parent_sd(1), parent_dm(1), child_sd(1), child_dm(1)

    parent_sd(1) = TRACK_ID_INVALID
    parent_dm(1) = TRACK_ID_INVALID
    child_sd(1) = child_sd_id
    child_dm(1) = child_dm_id
    call sdm_lifecycle_outnetcdf(otime, LIFE_SEEDING_ENTRY, 0, 1, &
         parent_sd, parent_dm, child_sd, child_dm, x, y, z)
  end subroutine sdm_tracking_lifecycle_write_seeding_entry
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine check_netcdf(status)
    use netcdf
    use scale_process, only: &
         PRC_MPIstop

    integer, intent (in) :: status
    
    if(status /= nf90_noerr) then 
       write(*,*) "sdm_netcdf_out: Write error ",status
      call PRC_MPIstop
    end if
  end subroutine check_netcdf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sdm_event_collision_outnetcdf(otime, num_pair, sd_id1, sd_id2, dm_id1, dm_id2, event_multiplicity, &
                                trigger_code, trigger_level, phase_state1_pre, phase_state2_pre, &
                                phase_state1_post, phase_state2_post, event_x, event_y, event_z, &
                                sd_n1_pre, sd_n2_pre, sd_n1_post, sd_n2_post, &
                                hydro_radius1_pre, hydro_radius2_pre, hydro_radius1_post, hydro_radius2_post, &
                                hydro_mass1_pre, hydro_mass2_pre, hydro_mass1_post, hydro_mass2_post, &
                                x1, y1, z1, x2, y2, z2, &
                                ice_re1_pre, ice_rp1_pre, ice_rho1_pre, ice_re1_post, ice_rp1_post, ice_rho1_post, &
                                ice_re2_pre, ice_rp2_pre, ice_rho2_pre, ice_re2_post, ice_rp2_post, ice_rho2_post, &
                                rime_mass1_pre, rime_mass1_post, rime_mass2_pre, rime_mass2_post, &
                                rime_frac1_pre, rime_frac1_post, rime_frac2_pre, rime_frac2_post, &
                                aspect_ratio1_pre, aspect_ratio1_post, aspect_ratio2_pre, aspect_ratio2_post)
    use netcdf
    use scale_precision
    use scale_stdio
    use scale_process, only: &
         mype => PRC_myrank
    use m_sdm_common, only: &
         sdm_cold, sdm_dmpitvl, tracking_mode, &
         TRACK_COLD_OUTPUT_ICE_GEOMETRY, TRACK_COLD_OUTPUT_RIME_MORPHOLOGY
    use m_sdm_tracking_cold, only: &
         TARGET_BY_EVENT

    implicit none

    real(DP), intent(in) :: otime
    integer, intent(in) :: num_pair
    integer, allocatable, intent(in) :: sd_id1(:)
    integer, allocatable, intent(in) :: sd_id2(:)
    integer, allocatable, intent(in) :: dm_id1(:)
    integer, allocatable, intent(in) :: dm_id2(:)
    integer, allocatable, intent(in) :: event_multiplicity(:)
    integer, allocatable, intent(in) :: trigger_code(:)
    integer, allocatable, intent(in) :: trigger_level(:)
    integer, allocatable, intent(in) :: phase_state1_pre(:)
    integer, allocatable, intent(in) :: phase_state2_pre(:)
    integer, allocatable, intent(in) :: phase_state1_post(:)
    integer, allocatable, intent(in) :: phase_state2_post(:)
    real(RP), allocatable, intent(in) :: event_x(:)
    real(RP), allocatable, intent(in) :: event_y(:)
    real(RP), allocatable, intent(in) :: event_z(:)
    integer(DP), allocatable, intent(in) :: sd_n1_pre(:)
    integer(DP), allocatable, intent(in) :: sd_n2_pre(:)
    integer(DP), allocatable, intent(in) :: sd_n1_post(:)
    integer(DP), allocatable, intent(in) :: sd_n2_post(:)
    real(RP), allocatable, intent(in) :: hydro_radius1_pre(:)
    real(RP), allocatable, intent(in) :: hydro_radius2_pre(:)
    real(RP), allocatable, intent(in) :: hydro_radius1_post(:)
    real(RP), allocatable, intent(in) :: hydro_radius2_post(:)
    real(RP), allocatable, intent(in) :: hydro_mass1_pre(:)
    real(RP), allocatable, intent(in) :: hydro_mass2_pre(:)
    real(RP), allocatable, intent(in) :: hydro_mass1_post(:)
    real(RP), allocatable, intent(in) :: hydro_mass2_post(:)
    real(RP), allocatable, intent(in) :: x1(:)
    real(RP), allocatable, intent(in) :: y1(:)
    real(RP), allocatable, intent(in) :: z1(:)
    real(RP), allocatable, intent(in) :: x2(:)
    real(RP), allocatable, intent(in) :: y2(:)
    real(RP), allocatable, intent(in) :: z2(:)
    real(RP), allocatable, intent(in) :: ice_re1_pre(:)
    real(RP), allocatable, intent(in) :: ice_rp1_pre(:)
    real(RP), allocatable, intent(in) :: ice_rho1_pre(:)
    real(RP), allocatable, intent(in) :: ice_re1_post(:)
    real(RP), allocatable, intent(in) :: ice_rp1_post(:)
    real(RP), allocatable, intent(in) :: ice_rho1_post(:)
    real(RP), allocatable, intent(in) :: ice_re2_pre(:)
    real(RP), allocatable, intent(in) :: ice_rp2_pre(:)
    real(RP), allocatable, intent(in) :: ice_rho2_pre(:)
    real(RP), allocatable, intent(in) :: ice_re2_post(:)
    real(RP), allocatable, intent(in) :: ice_rp2_post(:)
    real(RP), allocatable, intent(in) :: ice_rho2_post(:)
    real(RP), allocatable, intent(in) :: rime_mass1_pre(:)
    real(RP), allocatable, intent(in) :: rime_mass1_post(:)
    real(RP), allocatable, intent(in) :: rime_mass2_pre(:)
    real(RP), allocatable, intent(in) :: rime_mass2_post(:)
    real(RP), allocatable, intent(in) :: rime_frac1_pre(:)
    real(RP), allocatable, intent(in) :: rime_frac1_post(:)
    real(RP), allocatable, intent(in) :: rime_frac2_pre(:)
    real(RP), allocatable, intent(in) :: rime_frac2_post(:)
    real(RP), allocatable, intent(in) :: aspect_ratio1_pre(:)
    real(RP), allocatable, intent(in) :: aspect_ratio1_post(:)
    real(RP), allocatable, intent(in) :: aspect_ratio2_pre(:)
    real(RP), allocatable, intent(in) :: aspect_ratio2_post(:)

    character(len=17) :: fmt2="(A, '.', A, I*.*)"
    character(len=H_LONG) :: ftmp
    character(len=H_LONG) :: basename_sd_out
    character(len=19) :: basename_time
    character(len=16) :: id1_name, id2_name, dm1_name, dm2_name
    integer :: nf90_real_precision
    integer :: ncid, event_dim_id, event_offset
    integer :: time_id, trigger_code_id, trigger_level_id, target_reason_mask_id, event_multiplicity_id
    integer :: sd_id1_id, sd_id2_id, dm_id1_id, dm_id2_id
    integer :: phase_state1_pre_id, phase_state2_pre_id
    integer :: phase_state1_post_id, phase_state2_post_id
    integer :: event_x_id, event_y_id, event_z_id
    integer :: sd_n1_pre_id, sd_n2_pre_id, sd_n1_post_id, sd_n2_post_id
    integer :: hydro_radius1_pre_id, hydro_radius2_pre_id
    integer :: hydro_radius1_post_id, hydro_radius2_post_id
    integer :: hydro_mass1_pre_id, hydro_mass2_pre_id
    integer :: hydro_mass1_post_id, hydro_mass2_post_id
    integer :: x1_id, y1_id, z1_id, x2_id, y2_id, z2_id
    integer :: ice_re1_pre_id, ice_rp1_pre_id, ice_rho1_pre_id
    integer :: ice_re1_post_id, ice_rp1_post_id, ice_rho1_post_id
    integer :: ice_re2_pre_id, ice_rp2_pre_id, ice_rho2_pre_id
    integer :: ice_re2_post_id, ice_rp2_post_id, ice_rho2_post_id
    integer :: rime_mass1_pre_id, rime_mass1_post_id
    integer :: rime_mass2_pre_id, rime_mass2_post_id
    integer :: rime_frac1_pre_id, rime_frac1_post_id
    integer :: rime_frac2_pre_id, rime_frac2_post_id
    integer :: aspect_ratio1_pre_id, aspect_ratio1_post_id
    integer :: aspect_ratio2_pre_id, aspect_ratio2_post_id
    integer :: start1(1), count1(1)
    integer :: otime_bucket_hh, otime_bucket_mm, otime_bucket_ss, otime_bucket_ms
    real(DP) :: otime_bucket, otime_bucket_daysec
    logical :: file_exists, write_ids, write_ice_geometry, write_rime_morphology
    real(DP), allocatable :: time_out(:)
    integer, allocatable :: target_reason_mask(:)
    integer, parameter :: nc_deflate_level = 1
    integer, parameter :: nc_deflate       = 1
    integer, parameter :: nc_shuffle       = 1

    if( .not. sdm_cold ) return
    if( num_pair <= 0 ) return

    if(RP == SP)then
       nf90_real_precision = NF90_FLOAT
    else
       nf90_real_precision = NF90_DOUBLE
    end if

    write_ids = allocated(sd_id1) .and. allocated(sd_id2) .and. allocated(dm_id1) .and. allocated(dm_id2)
    write_ice_geometry = TRACK_COLD_OUTPUT_ICE_GEOMETRY
    write_rime_morphology = TRACK_COLD_OUTPUT_RIME_MORPHOLOGY
    if( tracking_mode == 2 ) then
       id1_name = 'pre_sdid1'
       id2_name = 'pre_sdid2'
       dm1_name = 'pre_dmid1'
       dm2_name = 'pre_dmid2'
    else
       id1_name = 'sd_id1'
       id2_name = 'sd_id2'
       dm1_name = 'dm_id1'
       dm2_name = 'dm_id2'
    end if

    if( sdm_dmpitvl > 0.0_RP ) then
       otime_bucket = sdm_dmpitvl * real( int( otime / sdm_dmpitvl + 1.0e-10_RP ), kind=DP )
    else
       otime_bucket = otime
    end if
    otime_bucket_daysec = mod( otime_bucket, 86400.0_DP )
    otime_bucket_hh = int( otime_bucket_daysec / 3600.0_DP )
    otime_bucket_mm = int( mod(otime_bucket_daysec, 3600.0_DP) / 60.0_DP )
    otime_bucket_ss = int( mod(otime_bucket_daysec, 60.0_DP) )
    otime_bucket_ms = nint( (otime_bucket_daysec - real(int(otime_bucket_daysec),kind=DP)) * 1000.0_DP )
    if( otime_bucket_ms >= 1000 ) then
       otime_bucket_ms = otime_bucket_ms - 1000
       otime_bucket_ss = otime_bucket_ss + 1
       if( otime_bucket_ss >= 60 ) then
          otime_bucket_ss = otime_bucket_ss - 60
          otime_bucket_mm = otime_bucket_mm + 1
          if( otime_bucket_mm >= 60 ) then
             otime_bucket_mm = otime_bucket_mm - 60
             otime_bucket_hh = otime_bucket_hh + 1
          end if
       end if
    end if
    write(basename_time,'(I8.8,A1,I2.2,I2.2,I2.2,A1,I3.3)') 101, '-', &
         otime_bucket_hh, otime_bucket_mm, otime_bucket_ss, '.', otime_bucket_ms
    write(fmt2(14:14),'(I1)') 6
    write(fmt2(16:16),'(I1)') 6
    write(ftmp,'(3A)') 'SD_event_collision', '_NetCDF_', trim(basename_time)
    write(basename_sd_out,fmt2) trim(ftmp), 'pe',mype

    inquire(file=trim(basename_sd_out), exist=file_exists)
    if( file_exists ) then
       call check_netcdf( nf90_open(trim(basename_sd_out), NF90_WRITE, ncid) )
       call check_netcdf( nf90_inq_dimid(ncid, "event", event_dim_id) )
       call check_netcdf( nf90_inquire_dimension(ncid, event_dim_id, len=event_offset) )
       call inq_int("trigger_code", trigger_code_id)
       call inq_int("trigger_level", trigger_level_id)
       call inq_int("target_reason_mask", target_reason_mask_id)
       call inq_int("event_multiplicity", event_multiplicity_id)
       call inq_int("phase_state1_pre", phase_state1_pre_id)
       call inq_int("phase_state2_pre", phase_state2_pre_id)
       call inq_int("phase_state1_post", phase_state1_post_id)
       call inq_int("phase_state2_post", phase_state2_post_id)
       call inq_real("time", time_id)
       call inq_real("x", event_x_id)
       call inq_real("y", event_y_id)
       call inq_real("z", event_z_id)
       call inq_i8("sd_n1_pre", sd_n1_pre_id)
       call inq_i8("sd_n2_pre", sd_n2_pre_id)
       call inq_i8("sd_n1_post", sd_n1_post_id)
       call inq_i8("sd_n2_post", sd_n2_post_id)
       call inq_real("hydro_radius1_pre", hydro_radius1_pre_id)
       call inq_real("hydro_radius2_pre", hydro_radius2_pre_id)
       call inq_real("hydro_radius1_post", hydro_radius1_post_id)
       call inq_real("hydro_radius2_post", hydro_radius2_post_id)
       call inq_real("hydro_mass1_pre", hydro_mass1_pre_id)
       call inq_real("hydro_mass2_pre", hydro_mass2_pre_id)
       call inq_real("hydro_mass1_post", hydro_mass1_post_id)
       call inq_real("hydro_mass2_post", hydro_mass2_post_id)
       call inq_real("x1", x1_id)
       call inq_real("y1", y1_id)
       call inq_real("z1", z1_id)
       call inq_real("x2", x2_id)
       call inq_real("y2", y2_id)
       call inq_real("z2", z2_id)
       if( write_ids ) then
          call inq_int(trim(id1_name), sd_id1_id)
          call inq_int(trim(dm1_name), dm_id1_id)
          call inq_int(trim(id2_name), sd_id2_id)
          call inq_int(trim(dm2_name), dm_id2_id)
       end if
       if( write_ice_geometry ) call inq_ice_geometry()
       if( write_rime_morphology ) call inq_rime_morphology()
    else
       call check_netcdf( nf90_create(trim(basename_sd_out), NF90_NETCDF4, ncid) )
       call check_netcdf( nf90_def_dim(ncid, "event", NF90_UNLIMITED, event_dim_id) )
       call def_time(time_id)
       call def_int("trigger_code", "cold-SDM collision process trigger code", "", trigger_code_id)
       call def_int("trigger_level", "cold-SDM process trigger level: 1=occurrence, 2=significant", "", &
            trigger_level_id)
       call def_int("target_reason_mask", "target reason bitmask", "", target_reason_mask_id)
       if( write_ids ) then
          call def_int(trim(id1_name), "participant-1 SD ID", "", sd_id1_id)
          call def_int(trim(dm1_name), "participant-1 domain ID", "", dm_id1_id)
          call def_int(trim(id2_name), "participant-2 SD ID", "", sd_id2_id)
          call def_int(trim(dm2_name), "participant-2 domain ID", "", dm_id2_id)
       end if
       call def_int("event_multiplicity", "collision-family event multiplicity", "", event_multiplicity_id)
       call def_int("phase_state1_pre", "participant-1 pre-event phase state", "", phase_state1_pre_id)
       call def_int("phase_state2_pre", "participant-2 pre-event phase state", "", phase_state2_pre_id)
       call def_int("phase_state1_post", "participant-1 post-event phase state", "", phase_state1_post_id)
       call def_int("phase_state2_post", "participant-2 post-event phase state", "", phase_state2_post_id)
       ! x/y/z is the pair midpoint; x1/y1/z1 and x2/y2/z2 are participant positions.
       call def_real("x", "collision pair pre-event midpoint x", "m", event_x_id)
       call def_real("y", "collision pair pre-event midpoint y", "m", event_y_id)
       call def_real("z", "collision pair pre-event midpoint z", "m", event_z_id)
       call def_i8("sd_n1_pre", "participant-1 multiplicity before event", "", sd_n1_pre_id)
       call def_i8("sd_n2_pre", "participant-2 multiplicity before event", "", sd_n2_pre_id)
       call def_i8("sd_n1_post", "participant-1 multiplicity after event", "", sd_n1_post_id)
       call def_i8("sd_n2_post", "participant-2 multiplicity after event", "", sd_n2_post_id)
       call def_real("hydro_radius1_pre", "participant-1 phase-aware volume-equivalent radius before event", "m", &
            hydro_radius1_pre_id)
       call def_real("hydro_radius2_pre", "participant-2 phase-aware volume-equivalent radius before event", "m", &
            hydro_radius2_pre_id)
       call def_real("hydro_radius1_post", "participant-1 phase-aware volume-equivalent radius after event", "m", &
            hydro_radius1_post_id)
       call def_real("hydro_radius2_post", "participant-2 phase-aware volume-equivalent radius after event", "m", &
            hydro_radius2_post_id)
       call def_real("hydro_mass1_pre", "participant-1 hydrometeor mass before event", "kg", hydro_mass1_pre_id)
       call def_real("hydro_mass2_pre", "participant-2 hydrometeor mass before event", "kg", hydro_mass2_pre_id)
       call def_real("hydro_mass1_post", "participant-1 hydrometeor mass after event", "kg", hydro_mass1_post_id)
       call def_real("hydro_mass2_post", "participant-2 hydrometeor mass after event", "kg", hydro_mass2_post_id)
       call def_real("x1", "participant-1 pre-event x", "m", x1_id)
       call def_real("y1", "participant-1 pre-event y", "m", y1_id)
       call def_real("z1", "participant-1 pre-event z", "m", z1_id)
       call def_real("x2", "participant-2 pre-event x", "m", x2_id)
       call def_real("y2", "participant-2 pre-event y", "m", y2_id)
       call def_real("z2", "participant-2 pre-event z", "m", z2_id)
       if( write_ice_geometry ) call def_ice_geometry()
       if( write_rime_morphology ) call def_rime_morphology()
       call check_netcdf( nf90_enddef(ncid) )
       event_offset = 0
    end if

    allocate(time_out(num_pair))
    allocate(target_reason_mask(num_pair))
    time_out(:) = otime
    target_reason_mask(:) = TARGET_BY_EVENT
    start1(1) = event_offset + 1
    count1(1) = num_pair

    call put_time(time_id, time_out)
    call put_int(trigger_code_id, trigger_code)
    call put_int(trigger_level_id, trigger_level)
    call put_int(target_reason_mask_id, target_reason_mask)
    if( write_ids ) then
       call put_int(sd_id1_id, sd_id1)
       call put_int(dm_id1_id, dm_id1)
       call put_int(sd_id2_id, sd_id2)
       call put_int(dm_id2_id, dm_id2)
    end if
    call put_int(event_multiplicity_id, event_multiplicity)
    call put_int(phase_state1_pre_id, phase_state1_pre)
    call put_int(phase_state2_pre_id, phase_state2_pre)
    call put_int(phase_state1_post_id, phase_state1_post)
    call put_int(phase_state2_post_id, phase_state2_post)
    call put_real(event_x_id, event_x)
    call put_real(event_y_id, event_y)
    call put_real(event_z_id, event_z)
    call put_i8(sd_n1_pre_id, sd_n1_pre)
    call put_i8(sd_n2_pre_id, sd_n2_pre)
    call put_i8(sd_n1_post_id, sd_n1_post)
    call put_i8(sd_n2_post_id, sd_n2_post)
    call put_real(hydro_radius1_pre_id, hydro_radius1_pre)
    call put_real(hydro_radius2_pre_id, hydro_radius2_pre)
    call put_real(hydro_radius1_post_id, hydro_radius1_post)
    call put_real(hydro_radius2_post_id, hydro_radius2_post)
    call put_real(hydro_mass1_pre_id, hydro_mass1_pre)
    call put_real(hydro_mass2_pre_id, hydro_mass2_pre)
    call put_real(hydro_mass1_post_id, hydro_mass1_post)
    call put_real(hydro_mass2_post_id, hydro_mass2_post)
    call put_real(x1_id, x1)
    call put_real(y1_id, y1)
    call put_real(z1_id, z1)
    call put_real(x2_id, x2)
    call put_real(y2_id, y2)
    call put_real(z2_id, z2)
    if( write_ice_geometry ) call put_ice_geometry()
    if( write_rime_morphology ) call put_rime_morphology()

    deallocate(time_out, target_reason_mask)
    call check_netcdf( nf90_close(ncid) )

  contains
    subroutine def_int(name, long_name, units, varid)
      character(len=*), intent(in) :: name, long_name, units
      integer, intent(out) :: varid
      call check_netcdf( nf90_def_var(ncid, trim(name), NF90_INT, event_dim_id, varid) )
      call check_netcdf( nf90_def_var_deflate(ncid, varid, &
           shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
      call check_netcdf( nf90_put_att(ncid, varid, 'long_name', trim(long_name)) )
      call check_netcdf( nf90_put_att(ncid, varid, 'units', trim(units)) )
    end subroutine def_int

    subroutine def_i8(name, long_name, units, varid)
      character(len=*), intent(in) :: name, long_name, units
      integer, intent(out) :: varid
      call check_netcdf( nf90_def_var(ncid, trim(name), NF90_INT64, event_dim_id, varid) )
      call check_netcdf( nf90_def_var_deflate(ncid, varid, &
           shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
      call check_netcdf( nf90_put_att(ncid, varid, 'long_name', trim(long_name)) )
      call check_netcdf( nf90_put_att(ncid, varid, 'units', trim(units)) )
    end subroutine def_i8

    subroutine def_real(name, long_name, units, varid)
      character(len=*), intent(in) :: name, long_name, units
      integer, intent(out) :: varid
      call check_netcdf( nf90_def_var(ncid, trim(name), nf90_real_precision, event_dim_id, varid) )
      call check_netcdf( nf90_def_var_deflate(ncid, varid, &
           shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
      call check_netcdf( nf90_put_att(ncid, varid, 'long_name', trim(long_name)) )
      call check_netcdf( nf90_put_att(ncid, varid, 'units', trim(units)) )
    end subroutine def_real

    subroutine def_time(varid)
      integer, intent(out) :: varid
      call check_netcdf( nf90_def_var(ncid, "time", NF90_DOUBLE, event_dim_id, varid) )
      call check_netcdf( nf90_def_var_deflate(ncid, varid, &
           shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
      call check_netcdf( nf90_put_att(ncid, varid, 'long_name', 'event time') )
      call check_netcdf( nf90_put_att(ncid, varid, 'units', 's') )
    end subroutine def_time

    subroutine inq_int(name, varid)
      character(len=*), intent(in) :: name
      integer, intent(out) :: varid
      call check_netcdf( nf90_inq_varid(ncid, trim(name), varid) )
    end subroutine inq_int

    subroutine inq_i8(name, varid)
      character(len=*), intent(in) :: name
      integer, intent(out) :: varid
      call check_netcdf( nf90_inq_varid(ncid, trim(name), varid) )
    end subroutine inq_i8

    subroutine inq_real(name, varid)
      character(len=*), intent(in) :: name
      integer, intent(out) :: varid
      call check_netcdf( nf90_inq_varid(ncid, trim(name), varid) )
    end subroutine inq_real

    subroutine put_int(varid, values)
      integer, intent(in) :: varid
      integer, intent(in) :: values(:)
      call check_netcdf( nf90_put_var(ncid, varid, values, start=start1, count=count1) )
    end subroutine put_int

    subroutine put_i8(varid, values)
      integer, intent(in) :: varid
      integer(DP), intent(in) :: values(:)
      call check_netcdf( nf90_put_var(ncid, varid, values, start=start1, count=count1) )
    end subroutine put_i8

    subroutine put_real(varid, values)
      integer, intent(in) :: varid
      real(RP), intent(in) :: values(:)
      call check_netcdf( nf90_put_var(ncid, varid, values, start=start1, count=count1) )
    end subroutine put_real

    subroutine put_time(varid, values)
      integer, intent(in) :: varid
      real(DP), intent(in) :: values(:)
      call check_netcdf( nf90_put_var(ncid, varid, values, start=start1, count=count1) )
    end subroutine put_time

    subroutine def_ice_geometry()
      call def_real("ice_re1_pre", "participant-1 ice equatorial radius before event", "m", ice_re1_pre_id)
      call def_real("ice_rp1_pre", "participant-1 ice polar radius before event", "m", ice_rp1_pre_id)
      call def_real("ice_rho1_pre", "participant-1 ice density before event", "kg/m3", ice_rho1_pre_id)
      call def_real("ice_re1_post", "participant-1 ice equatorial radius after event", "m", ice_re1_post_id)
      call def_real("ice_rp1_post", "participant-1 ice polar radius after event", "m", ice_rp1_post_id)
      call def_real("ice_rho1_post", "participant-1 ice density after event", "kg/m3", ice_rho1_post_id)
      call def_real("ice_re2_pre", "participant-2 ice equatorial radius before event", "m", ice_re2_pre_id)
      call def_real("ice_rp2_pre", "participant-2 ice polar radius before event", "m", ice_rp2_pre_id)
      call def_real("ice_rho2_pre", "participant-2 ice density before event", "kg/m3", ice_rho2_pre_id)
      call def_real("ice_re2_post", "participant-2 ice equatorial radius after event", "m", ice_re2_post_id)
      call def_real("ice_rp2_post", "participant-2 ice polar radius after event", "m", ice_rp2_post_id)
      call def_real("ice_rho2_post", "participant-2 ice density after event", "kg/m3", ice_rho2_post_id)
    end subroutine def_ice_geometry

    subroutine def_rime_morphology()
      call def_real("rime_mass1_pre", "participant-1 rime mass before event", "kg", rime_mass1_pre_id)
      call def_real("rime_mass1_post", "participant-1 rime mass after event", "kg", rime_mass1_post_id)
      call def_real("rime_mass2_pre", "participant-2 rime mass before event", "kg", rime_mass2_pre_id)
      call def_real("rime_mass2_post", "participant-2 rime mass after event", "kg", rime_mass2_post_id)
      call def_real("rime_frac1_pre", "participant-1 rime mass fraction before event", "", rime_frac1_pre_id)
      call def_real("rime_frac1_post", "participant-1 rime mass fraction after event", "", rime_frac1_post_id)
      call def_real("rime_frac2_pre", "participant-2 rime mass fraction before event", "", rime_frac2_pre_id)
      call def_real("rime_frac2_post", "participant-2 rime mass fraction after event", "", rime_frac2_post_id)
      call def_real("aspect_ratio1_pre", "participant-1 ice aspect ratio before event", "", aspect_ratio1_pre_id)
      call def_real("aspect_ratio1_post", "participant-1 ice aspect ratio after event", "", aspect_ratio1_post_id)
      call def_real("aspect_ratio2_pre", "participant-2 ice aspect ratio before event", "", aspect_ratio2_pre_id)
      call def_real("aspect_ratio2_post", "participant-2 ice aspect ratio after event", "", aspect_ratio2_post_id)
    end subroutine def_rime_morphology

    subroutine inq_ice_geometry()
      call inq_real("ice_re1_pre", ice_re1_pre_id)
      call inq_real("ice_rp1_pre", ice_rp1_pre_id)
      call inq_real("ice_rho1_pre", ice_rho1_pre_id)
      call inq_real("ice_re1_post", ice_re1_post_id)
      call inq_real("ice_rp1_post", ice_rp1_post_id)
      call inq_real("ice_rho1_post", ice_rho1_post_id)
      call inq_real("ice_re2_pre", ice_re2_pre_id)
      call inq_real("ice_rp2_pre", ice_rp2_pre_id)
      call inq_real("ice_rho2_pre", ice_rho2_pre_id)
      call inq_real("ice_re2_post", ice_re2_post_id)
      call inq_real("ice_rp2_post", ice_rp2_post_id)
      call inq_real("ice_rho2_post", ice_rho2_post_id)
    end subroutine inq_ice_geometry

    subroutine inq_rime_morphology()
      call inq_real("rime_mass1_pre", rime_mass1_pre_id)
      call inq_real("rime_mass1_post", rime_mass1_post_id)
      call inq_real("rime_mass2_pre", rime_mass2_pre_id)
      call inq_real("rime_mass2_post", rime_mass2_post_id)
      call inq_real("rime_frac1_pre", rime_frac1_pre_id)
      call inq_real("rime_frac1_post", rime_frac1_post_id)
      call inq_real("rime_frac2_pre", rime_frac2_pre_id)
      call inq_real("rime_frac2_post", rime_frac2_post_id)
      call inq_real("aspect_ratio1_pre", aspect_ratio1_pre_id)
      call inq_real("aspect_ratio1_post", aspect_ratio1_post_id)
      call inq_real("aspect_ratio2_pre", aspect_ratio2_pre_id)
      call inq_real("aspect_ratio2_post", aspect_ratio2_post_id)
    end subroutine inq_rime_morphology

    subroutine put_ice_geometry()
      call put_real(ice_re1_pre_id, ice_re1_pre)
      call put_real(ice_rp1_pre_id, ice_rp1_pre)
      call put_real(ice_rho1_pre_id, ice_rho1_pre)
      call put_real(ice_re1_post_id, ice_re1_post)
      call put_real(ice_rp1_post_id, ice_rp1_post)
      call put_real(ice_rho1_post_id, ice_rho1_post)
      call put_real(ice_re2_pre_id, ice_re2_pre)
      call put_real(ice_rp2_pre_id, ice_rp2_pre)
      call put_real(ice_rho2_pre_id, ice_rho2_pre)
      call put_real(ice_re2_post_id, ice_re2_post)
      call put_real(ice_rp2_post_id, ice_rp2_post)
      call put_real(ice_rho2_post_id, ice_rho2_post)
    end subroutine put_ice_geometry

    subroutine put_rime_morphology()
      call put_real(rime_mass1_pre_id, rime_mass1_pre)
      call put_real(rime_mass1_post_id, rime_mass1_post)
      call put_real(rime_mass2_pre_id, rime_mass2_pre)
      call put_real(rime_mass2_post_id, rime_mass2_post)
      call put_real(rime_frac1_pre_id, rime_frac1_pre)
      call put_real(rime_frac1_post_id, rime_frac1_post)
      call put_real(rime_frac2_pre_id, rime_frac2_pre)
      call put_real(rime_frac2_post_id, rime_frac2_post)
      call put_real(aspect_ratio1_pre_id, aspect_ratio1_pre)
      call put_real(aspect_ratio1_post_id, aspect_ratio1_post)
      call put_real(aspect_ratio2_pre_id, aspect_ratio2_pre)
      call put_real(aspect_ratio2_post_id, aspect_ratio2_post)
    end subroutine put_rime_morphology
  end subroutine sdm_event_collision_outnetcdf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sdm_coal_outnetcdf(otime, num_pair, sd_id1, sd_id2, dm_id1, dm_id2, num_col, sdr1_out, sdr2_out, sdn1_out, sdn2_out, &
       &                        event_x, event_y, event_z, sdn1_post, sdn2_post, &
       &                        hydro_radius1_pre, hydro_radius2_pre, hydro_radius1_post, hydro_radius2_post, &
       &                        hydro_mass1_pre, hydro_mass2_pre, hydro_mass1_post, hydro_mass2_post)
    use netcdf
    use scale_precision
    use scale_stdio
    use scale_process, only: &
         mype => PRC_myrank
    use m_sdm_common, only: &
         sdm_dmpitvl, tracking_mode

    implicit none

    real(DP), intent(in) :: otime
    integer, intent(in) :: num_pair
    integer, allocatable, intent(in) :: sd_id1(:)
    integer, allocatable, intent(in) :: dm_id1(:)
    integer, allocatable, intent(in) :: sd_id2(:)
    integer, allocatable, intent(in) :: dm_id2(:)
    integer, allocatable, intent(in) :: num_col(:)
    real(RP), allocatable, intent(in) :: sdr1_out(:)
    real(RP), allocatable, intent(in) :: sdr2_out(:)
    integer(DP), allocatable, intent(in) :: sdn1_out(:)
    integer(DP), allocatable, intent(in) :: sdn2_out(:)
    real(RP), allocatable, intent(in) :: event_x(:)
    real(RP), allocatable, intent(in) :: event_y(:)
    real(RP), allocatable, intent(in) :: event_z(:)
    integer(DP), allocatable, intent(in) :: sdn1_post(:)
    integer(DP), allocatable, intent(in) :: sdn2_post(:)
    real(RP), allocatable, intent(in) :: hydro_radius1_pre(:)
    real(RP), allocatable, intent(in) :: hydro_radius2_pre(:)
    real(RP), allocatable, intent(in) :: hydro_radius1_post(:)
    real(RP), allocatable, intent(in) :: hydro_radius2_post(:)
    real(RP), allocatable, intent(in) :: hydro_mass1_pre(:)
    real(RP), allocatable, intent(in) :: hydro_mass2_pre(:)
    real(RP), allocatable, intent(in) :: hydro_mass1_post(:)
    real(RP), allocatable, intent(in) :: hydro_mass2_post(:)

    character(len=17) :: fmt2="(A, '.', A, I*.*)"
    character(len=H_LONG) :: ftmp
    character(len=H_LONG) :: basename_sd_out
    character(len=19) :: basename_time
    integer :: nf90_real_precision
    integer :: ncid, event_dim_id, sdr1_out_id, sdr2_out_id, sdn1_out_id, sdn2_out_id
    integer :: sd_id1_id, sd_id2_id, dm_id1_id, dm_id2_id, event_multiplicity_id, event_time_id
    integer :: event_x_id, event_y_id, event_z_id
    integer :: sd_n1_pre_id, sd_n2_pre_id, sd_n1_post_id, sd_n2_post_id
    integer :: hydro_radius1_pre_id, hydro_radius2_pre_id
    integer :: hydro_radius1_post_id, hydro_radius2_post_id
    integer :: hydro_mass1_pre_id, hydro_mass2_pre_id
    integer :: hydro_mass1_post_id, hydro_mass2_post_id
    integer :: event_offset
    integer :: start1(1), count1(1)
    logical :: file_exists
    real(DP), allocatable :: event_time(:)
    real(DP) :: otime_bucket
    real(DP) :: otime_bucket_daysec
    integer :: otime_bucket_hh
    integer :: otime_bucket_mm
    integer :: otime_bucket_ss
    integer :: otime_bucket_ms
    character(len=16) :: id1_name, id2_name, dm1_name, dm2_name
    logical :: write_tracking_ids

    integer,parameter :: nc_deflate_level = 1
    integer,parameter :: nc_deflate       = 1
    integer,parameter :: nc_shuffle       = 1

    if( num_pair <= 0 ) return

    write_tracking_ids = tracking_mode /= 0
    if( tracking_mode == 2 ) then
       id1_name = 'pre_sdid1'
       id2_name = 'pre_sdid2'
       dm1_name = 'pre_dmid1'
       dm2_name = 'pre_dmid2'
    else if( tracking_mode == 1 ) then
       id1_name = 'sd_id1'
       id2_name = 'sd_id2'
       dm1_name = 'dm_id1'
       dm2_name = 'dm_id2'
    else
       id1_name = ''
       id2_name = ''
       dm1_name = ''
       dm2_name = ''
    end if

    write(fmt2(14:14),'(I1)') 6
    write(fmt2(16:16),'(I1)') 6

    if(RP == SP)then
       nf90_real_precision = NF90_FLOAT
    else
       nf90_real_precision = NF90_DOUBLE
    end if

    if( sdm_dmpitvl > 0.0_RP ) then
       otime_bucket = sdm_dmpitvl * real( int( otime / sdm_dmpitvl + 1.0e-10_RP ), kind=DP )
    else
       otime_bucket = otime
    end if

    otime_bucket_daysec = mod( otime_bucket, 86400.0_DP )
    otime_bucket_hh = int( otime_bucket_daysec / 3600.0_DP )
    otime_bucket_mm = int( mod(otime_bucket_daysec, 3600.0_DP) / 60.0_DP )
    otime_bucket_ss = int( mod(otime_bucket_daysec,   60.0_DP) )
    otime_bucket_ms = nint( (otime_bucket_daysec - real(int(otime_bucket_daysec),kind=DP)) * 1000.0_DP )

    if( otime_bucket_ms >= 1000 ) then
       otime_bucket_ms = otime_bucket_ms - 1000
       otime_bucket_ss = otime_bucket_ss + 1
       if( otime_bucket_ss >= 60 ) then
          otime_bucket_ss = otime_bucket_ss - 60
          otime_bucket_mm = otime_bucket_mm + 1
          if( otime_bucket_mm >= 60 ) then
             otime_bucket_mm = otime_bucket_mm - 60
             otime_bucket_hh = otime_bucket_hh + 1
          end if
       end if
    end if

    write(basename_time,'(I8.8,A1,I2.2,I2.2,I2.2,A1,I3.3)') 101, '-', &
         otime_bucket_hh, otime_bucket_mm, otime_bucket_ss, '.', otime_bucket_ms
    write(ftmp,'(3A)') 'SD_coal_output', '_NetCDF_', trim(basename_time)
    write(basename_sd_out,fmt2) trim(ftmp), 'pe',mype

    inquire(file=trim(basename_sd_out), exist=file_exists)

    if( file_exists ) then
       call check_netcdf( nf90_open(trim(basename_sd_out), NF90_WRITE, ncid) )
       call check_netcdf( nf90_inq_dimid(ncid, "event", event_dim_id) )
       call check_netcdf( nf90_inquire_dimension(ncid, event_dim_id, len=event_offset) )
       call check_netcdf( nf90_inq_varid(ncid, "event_time", event_time_id) )
       if( write_tracking_ids ) then
          call check_netcdf( nf90_inq_varid(ncid, trim(id1_name), sd_id1_id) )
          call check_netcdf( nf90_inq_varid(ncid, trim(dm1_name), dm_id1_id) )
       end if
       call check_netcdf( nf90_inq_varid(ncid, "sd_r1", sdr1_out_id) )
       call check_netcdf( nf90_inq_varid(ncid, "sd_n1", sdn1_out_id) )
       if( write_tracking_ids ) then
          call check_netcdf( nf90_inq_varid(ncid, trim(id2_name), sd_id2_id) )
          call check_netcdf( nf90_inq_varid(ncid, trim(dm2_name), dm_id2_id) )
       end if
       call check_netcdf( nf90_inq_varid(ncid, "sd_r2", sdr2_out_id) )
       call check_netcdf( nf90_inq_varid(ncid, "sd_n2", sdn2_out_id) )
       call check_netcdf( nf90_inq_varid(ncid, "event_multiplicity", event_multiplicity_id) )
       call check_netcdf( nf90_inq_varid(ncid, "x", event_x_id) )
       call check_netcdf( nf90_inq_varid(ncid, "y", event_y_id) )
       call check_netcdf( nf90_inq_varid(ncid, "z", event_z_id) )
       call check_netcdf( nf90_inq_varid(ncid, "sd_n1_pre", sd_n1_pre_id) )
       call check_netcdf( nf90_inq_varid(ncid, "sd_n2_pre", sd_n2_pre_id) )
       call check_netcdf( nf90_inq_varid(ncid, "sd_n1_post", sd_n1_post_id) )
       call check_netcdf( nf90_inq_varid(ncid, "sd_n2_post", sd_n2_post_id) )
       call check_netcdf( nf90_inq_varid(ncid, "hydro_radius1_pre", hydro_radius1_pre_id) )
       call check_netcdf( nf90_inq_varid(ncid, "hydro_radius2_pre", hydro_radius2_pre_id) )
       call check_netcdf( nf90_inq_varid(ncid, "hydro_radius1_post", hydro_radius1_post_id) )
       call check_netcdf( nf90_inq_varid(ncid, "hydro_radius2_post", hydro_radius2_post_id) )
       call check_netcdf( nf90_inq_varid(ncid, "hydro_mass1_pre", hydro_mass1_pre_id) )
       call check_netcdf( nf90_inq_varid(ncid, "hydro_mass2_pre", hydro_mass2_pre_id) )
       call check_netcdf( nf90_inq_varid(ncid, "hydro_mass1_post", hydro_mass1_post_id) )
       call check_netcdf( nf90_inq_varid(ncid, "hydro_mass2_post", hydro_mass2_post_id) )
    else
       call check_netcdf( nf90_create(trim(basename_sd_out),NF90_NETCDF4, ncid) )
       call check_netcdf( nf90_def_dim(ncid, "event", NF90_UNLIMITED, event_dim_id) )
       call check_netcdf( nf90_def_var(ncid, "event_time", NF90_DOUBLE, event_dim_id, event_time_id) )
       call check_netcdf( nf90_def_var_deflate(ncid, event_time_id, &
            & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
       call check_netcdf( nf90_put_att(ncid, event_time_id, 'long_name', 'event time') )
       call check_netcdf( nf90_put_att(ncid, event_time_id, 'units', 's') )
       if( write_tracking_ids ) then
          call check_netcdf( nf90_def_var(ncid, trim(id1_name), NF90_INT, event_dim_id,sd_id1_id) )
          call check_netcdf( nf90_def_var_deflate(ncid, sd_id1_id, &
               & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
          call check_netcdf( nf90_put_att(ncid, sd_id1_id, 'long_name', 'index of SDs with large multiplicity') )
          call check_netcdf( nf90_put_att(ncid, sd_id1_id, 'units', '') )
          call check_netcdf( nf90_def_var(ncid, trim(dm1_name), NF90_INT, event_dim_id, dm_id1_id) )
          call check_netcdf( nf90_def_var_deflate(ncid, dm_id1_id, &
               & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
          call check_netcdf( nf90_put_att(ncid, dm_id1_id, 'long_name', 'domain ID of SDs with large multiplicity') )
          call check_netcdf( nf90_put_att(ncid, dm_id1_id, 'units', '') )
       end if
       call check_netcdf( nf90_def_var(ncid, "sd_r1", nf90_real_precision, event_dim_id, sdr1_out_id) )
       call check_netcdf( nf90_def_var_deflate(ncid, sdr1_out_id, &
            & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
       call check_netcdf( nf90_put_att(ncid, sdr1_out_id, 'long_name', 'equivalent radius of SDs with large multiplicity') )
       call check_netcdf( nf90_put_att(ncid, sdr1_out_id, 'units', 'm') )
       call check_netcdf( nf90_def_var(ncid, "sd_n1", NF90_INT64, event_dim_id, sdn1_out_id) )
       call check_netcdf( nf90_def_var_deflate(ncid, sdn1_out_id, &
            & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
       call check_netcdf( nf90_put_att(ncid, sdn1_out_id, 'long_name', 'multiplicity of SDs with large multiplicity') )
       call check_netcdf( nf90_put_att(ncid, sdn1_out_id, 'units', '') )
       if( write_tracking_ids ) then
          call check_netcdf( nf90_def_var(ncid, trim(id2_name), NF90_INT, event_dim_id,sd_id2_id) )
          call check_netcdf( nf90_def_var_deflate(ncid, sd_id2_id, &
               & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
          call check_netcdf( nf90_put_att(ncid, sd_id2_id, 'long_name', 'index of SDs with small multiplicity') )
          call check_netcdf( nf90_put_att(ncid, sd_id2_id, 'units', '') )
          call check_netcdf( nf90_def_var(ncid, trim(dm2_name), NF90_INT, event_dim_id, dm_id2_id) )
          call check_netcdf( nf90_def_var_deflate(ncid, dm_id2_id, &
               & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
          call check_netcdf( nf90_put_att(ncid, dm_id2_id, 'long_name', 'domain ID of SDs with small multiplicity') )
          call check_netcdf( nf90_put_att(ncid, dm_id2_id, 'units', '') )
       end if
       call check_netcdf( nf90_def_var(ncid, "sd_r2", nf90_real_precision, event_dim_id, sdr2_out_id) )
       call check_netcdf( nf90_def_var_deflate(ncid, sdr2_out_id, &
            & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
       call check_netcdf( nf90_put_att(ncid, sdr2_out_id, 'long_name', 'equivalent radius of SDs with small multiplicity') )
       call check_netcdf( nf90_put_att(ncid, sdr2_out_id, 'units', 'm') )
       call check_netcdf( nf90_def_var(ncid, "sd_n2", NF90_INT64, event_dim_id, sdn2_out_id) )
       call check_netcdf( nf90_def_var_deflate(ncid, sdn2_out_id, &
            & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
       call check_netcdf( nf90_put_att(ncid, sdn2_out_id, 'long_name', 'multiplicity of SDs with small multiplicity') )
       call check_netcdf( nf90_put_att(ncid, sdn2_out_id, 'units', '') )
       call check_netcdf( nf90_def_var(ncid, "event_multiplicity", NF90_INT, event_dim_id, event_multiplicity_id) )
       call check_netcdf( nf90_def_var_deflate(ncid, event_multiplicity_id, &
            & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
       call check_netcdf( nf90_put_att(ncid, event_multiplicity_id, 'long_name', 'collision-family event multiplicity') )
       call check_netcdf( nf90_put_att(ncid, event_multiplicity_id, 'units', '') )
       call check_netcdf( nf90_def_var(ncid, "x", nf90_real_precision, event_dim_id, event_x_id) )
       call check_netcdf( nf90_def_var_deflate(ncid, event_x_id, &
            & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
       call check_netcdf( nf90_put_att(ncid, event_x_id, 'long_name', 'pre-event pair midpoint x coordinate') )
       call check_netcdf( nf90_put_att(ncid, event_x_id, 'units', 'm') )
       call check_netcdf( nf90_def_var(ncid, "y", nf90_real_precision, event_dim_id, event_y_id) )
       call check_netcdf( nf90_def_var_deflate(ncid, event_y_id, &
            & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
       call check_netcdf( nf90_put_att(ncid, event_y_id, 'long_name', 'pre-event pair midpoint y coordinate') )
       call check_netcdf( nf90_put_att(ncid, event_y_id, 'units', 'm') )
       call check_netcdf( nf90_def_var(ncid, "z", nf90_real_precision, event_dim_id, event_z_id) )
       call check_netcdf( nf90_def_var_deflate(ncid, event_z_id, &
            & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
       call check_netcdf( nf90_put_att(ncid, event_z_id, 'long_name', 'pre-event pair midpoint z coordinate') )
       call check_netcdf( nf90_put_att(ncid, event_z_id, 'units', 'm') )
       call check_netcdf( nf90_def_var(ncid, "sd_n1_pre", NF90_INT64, event_dim_id, sd_n1_pre_id) )
       call check_netcdf( nf90_def_var_deflate(ncid, sd_n1_pre_id, &
            & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
       call check_netcdf( nf90_put_att(ncid, sd_n1_pre_id, 'long_name', 'participant-1 multiplicity before event') )
       call check_netcdf( nf90_put_att(ncid, sd_n1_pre_id, 'units', '') )
       call check_netcdf( nf90_def_var(ncid, "sd_n2_pre", NF90_INT64, event_dim_id, sd_n2_pre_id) )
       call check_netcdf( nf90_def_var_deflate(ncid, sd_n2_pre_id, &
            & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
       call check_netcdf( nf90_put_att(ncid, sd_n2_pre_id, 'long_name', 'participant-2 multiplicity before event') )
       call check_netcdf( nf90_put_att(ncid, sd_n2_pre_id, 'units', '') )
       call check_netcdf( nf90_def_var(ncid, "sd_n1_post", NF90_INT64, event_dim_id, sd_n1_post_id) )
       call check_netcdf( nf90_def_var_deflate(ncid, sd_n1_post_id, &
            & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
       call check_netcdf( nf90_put_att(ncid, sd_n1_post_id, 'long_name', 'participant-1 multiplicity after event') )
       call check_netcdf( nf90_put_att(ncid, sd_n1_post_id, 'units', '') )
       call check_netcdf( nf90_def_var(ncid, "sd_n2_post", NF90_INT64, event_dim_id, sd_n2_post_id) )
       call check_netcdf( nf90_def_var_deflate(ncid, sd_n2_post_id, &
            & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
       call check_netcdf( nf90_put_att(ncid, sd_n2_post_id, 'long_name', 'participant-2 multiplicity after event') )
       call check_netcdf( nf90_put_att(ncid, sd_n2_post_id, 'units', '') )
       call check_netcdf( nf90_def_var(ncid, "hydro_radius1_pre", nf90_real_precision, event_dim_id, hydro_radius1_pre_id) )
       call check_netcdf( nf90_def_var_deflate(ncid, hydro_radius1_pre_id, &
            & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
       call check_netcdf( nf90_put_att(ncid, hydro_radius1_pre_id, 'long_name', 'participant-1 liquid radius before event') )
       call check_netcdf( nf90_put_att(ncid, hydro_radius1_pre_id, 'units', 'm') )
       call check_netcdf( nf90_def_var(ncid, "hydro_radius2_pre", nf90_real_precision, event_dim_id, hydro_radius2_pre_id) )
       call check_netcdf( nf90_def_var_deflate(ncid, hydro_radius2_pre_id, &
            & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
       call check_netcdf( nf90_put_att(ncid, hydro_radius2_pre_id, 'long_name', 'participant-2 liquid radius before event') )
       call check_netcdf( nf90_put_att(ncid, hydro_radius2_pre_id, 'units', 'm') )
       call check_netcdf( nf90_def_var(ncid, "hydro_radius1_post", nf90_real_precision, event_dim_id, hydro_radius1_post_id) )
       call check_netcdf( nf90_def_var_deflate(ncid, hydro_radius1_post_id, &
            & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
       call check_netcdf( nf90_put_att(ncid, hydro_radius1_post_id, 'long_name', 'participant-1 liquid radius after event') )
       call check_netcdf( nf90_put_att(ncid, hydro_radius1_post_id, 'units', 'm') )
       call check_netcdf( nf90_def_var(ncid, "hydro_radius2_post", nf90_real_precision, event_dim_id, hydro_radius2_post_id) )
       call check_netcdf( nf90_def_var_deflate(ncid, hydro_radius2_post_id, &
            & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
       call check_netcdf( nf90_put_att(ncid, hydro_radius2_post_id, 'long_name', 'participant-2 liquid radius after event') )
       call check_netcdf( nf90_put_att(ncid, hydro_radius2_post_id, 'units', 'm') )
       call check_netcdf( nf90_def_var(ncid, "hydro_mass1_pre", nf90_real_precision, event_dim_id, hydro_mass1_pre_id) )
       call check_netcdf( nf90_def_var_deflate(ncid, hydro_mass1_pre_id, &
            & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
       call check_netcdf( nf90_put_att(ncid, hydro_mass1_pre_id, 'long_name', 'participant-1 liquid-water mass before event') )
       call check_netcdf( nf90_put_att(ncid, hydro_mass1_pre_id, 'units', 'kg') )
       call check_netcdf( nf90_def_var(ncid, "hydro_mass2_pre", nf90_real_precision, event_dim_id, hydro_mass2_pre_id) )
       call check_netcdf( nf90_def_var_deflate(ncid, hydro_mass2_pre_id, &
            & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
       call check_netcdf( nf90_put_att(ncid, hydro_mass2_pre_id, 'long_name', 'participant-2 liquid-water mass before event') )
       call check_netcdf( nf90_put_att(ncid, hydro_mass2_pre_id, 'units', 'kg') )
       call check_netcdf( nf90_def_var(ncid, "hydro_mass1_post", nf90_real_precision, event_dim_id, hydro_mass1_post_id) )
       call check_netcdf( nf90_def_var_deflate(ncid, hydro_mass1_post_id, &
            & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
       call check_netcdf( nf90_put_att(ncid, hydro_mass1_post_id, 'long_name', 'participant-1 liquid-water mass after event') )
       call check_netcdf( nf90_put_att(ncid, hydro_mass1_post_id, 'units', 'kg') )
       call check_netcdf( nf90_def_var(ncid, "hydro_mass2_post", nf90_real_precision, event_dim_id, hydro_mass2_post_id) )
       call check_netcdf( nf90_def_var_deflate(ncid, hydro_mass2_post_id, &
            & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
       call check_netcdf( nf90_put_att(ncid, hydro_mass2_post_id, 'long_name', 'participant-2 liquid-water mass after event') )
       call check_netcdf( nf90_put_att(ncid, hydro_mass2_post_id, 'units', 'kg') )
       call check_netcdf( nf90_enddef(ncid) )
       event_offset = 0
    end if

    allocate(event_time(num_pair))
    event_time(:) = otime
    start1(1) = event_offset + 1
    count1(1) = num_pair
    call check_netcdf( nf90_put_var(ncid, event_time_id, event_time, start=start1, count=count1) )
    if( write_tracking_ids ) then
       call check_netcdf( nf90_put_var(ncid, sd_id1_id, sd_id1, start=start1, count=count1) )
       call check_netcdf( nf90_put_var(ncid, dm_id1_id, dm_id1, start=start1, count=count1) )
    end if
    call check_netcdf( nf90_put_var(ncid, sdr1_out_id, sdr1_out, start=start1, count=count1) )
    call check_netcdf( nf90_put_var(ncid, sdn1_out_id, sdn1_out, start=start1, count=count1) )
    if( write_tracking_ids ) then
       call check_netcdf( nf90_put_var(ncid, sd_id2_id, sd_id2, start=start1, count=count1) )
       call check_netcdf( nf90_put_var(ncid, dm_id2_id, dm_id2, start=start1, count=count1) )
    end if
    call check_netcdf( nf90_put_var(ncid, sdr2_out_id, sdr2_out, start=start1, count=count1) )
    call check_netcdf( nf90_put_var(ncid, sdn2_out_id, sdn2_out, start=start1, count=count1) )
    call check_netcdf( nf90_put_var(ncid, event_multiplicity_id, num_col, start=start1, count=count1) )
    call check_netcdf( nf90_put_var(ncid, event_x_id, event_x, start=start1, count=count1) )
    call check_netcdf( nf90_put_var(ncid, event_y_id, event_y, start=start1, count=count1) )
    call check_netcdf( nf90_put_var(ncid, event_z_id, event_z, start=start1, count=count1) )
    call check_netcdf( nf90_put_var(ncid, sd_n1_pre_id, sdn1_out, start=start1, count=count1) )
    call check_netcdf( nf90_put_var(ncid, sd_n2_pre_id, sdn2_out, start=start1, count=count1) )
    call check_netcdf( nf90_put_var(ncid, sd_n1_post_id, sdn1_post, start=start1, count=count1) )
    call check_netcdf( nf90_put_var(ncid, sd_n2_post_id, sdn2_post, start=start1, count=count1) )
    call check_netcdf( nf90_put_var(ncid, hydro_radius1_pre_id, hydro_radius1_pre, start=start1, count=count1) )
    call check_netcdf( nf90_put_var(ncid, hydro_radius2_pre_id, hydro_radius2_pre, start=start1, count=count1) )
    call check_netcdf( nf90_put_var(ncid, hydro_radius1_post_id, hydro_radius1_post, start=start1, count=count1) )
    call check_netcdf( nf90_put_var(ncid, hydro_radius2_post_id, hydro_radius2_post, start=start1, count=count1) )
    call check_netcdf( nf90_put_var(ncid, hydro_mass1_pre_id, hydro_mass1_pre, start=start1, count=count1) )
    call check_netcdf( nf90_put_var(ncid, hydro_mass2_pre_id, hydro_mass2_pre, start=start1, count=count1) )
    call check_netcdf( nf90_put_var(ncid, hydro_mass1_post_id, hydro_mass1_post, start=start1, count=count1) )
    call check_netcdf( nf90_put_var(ncid, hydro_mass2_post_id, hydro_mass2_post, start=start1, count=count1) )
    deallocate(event_time)

    call check_netcdf( nf90_close(ncid) )

    if( IO_L ) write(IO_FID_LOG,*) '*** Closed output file (NetCDF) of Super Droplet'

    return

  end subroutine sdm_coal_outnetcdf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module m_sdm_io
