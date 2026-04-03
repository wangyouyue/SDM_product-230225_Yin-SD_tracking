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
  public :: sdm_outasci,sdm_outnetcdf,sdm_outnetcdf_hist,sdm_coal_outnetcdf,sdm_assign_tracking_subset

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
         INVALID_i4, forward_tracking_enable, backward_tracking_enable, coalescence_output_enable

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
  subroutine sdm_assign_tracking_subset(sd_num, sd_z, sd_r, sd_id, dm_id, if_coal)
    use scale_precision
    use scale_stdio, only: &
         IO_L, IO_FID_LOG
    use scale_process, only: &
         mype => PRC_myrank
    use m_sdm_common, only: &
         i2, INVALID_i4, backward_tracking_enable, tracking_selection_mode, tracking_fraction, max_tracked_sds, &
         tracking_height_min, tracking_height_max, tracking_radius_min, tracking_radius_max, tracking_nz_bin, tracking_nr_bin, &
         tracking_min_per_bin, tracking_fallback_to_random, tracking_sample_initialized

    implicit none

    integer, intent(in) :: sd_num
    real(RP), intent(in) :: sd_z(1:sd_num)
    real(RP), intent(in) :: sd_r(1:sd_num)
    integer, intent(inout) :: sd_id(1:sd_num)
    integer, intent(inout) :: dm_id(1:sd_num)
    integer(kind=i2), intent(inout) :: if_coal(1:sd_num)

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

    if( (.not. backward_tracking_enable) .or. (tracking_fraction <= 0.0_RP) ) then
      do n=1,sd_num
        sd_id(n) = INVALID_i4
        dm_id(n) = INVALID_i4
        if_coal(n)  = 0_i2
      end do
      tracking_sample_initialized = .false.
      return
    end if

    if( tracking_fraction >= 1.0_RP .and. max_tracked_sds <= 0 ) then
      do n=1,sd_num
        sd_id(n) = n
        dm_id(n) = mype
        if_coal(n) = 0_i2
      end do
      tracking_sample_initialized = .true.
      return
    end if

    tracked_cnt = 0
    if( .not. tracking_sample_initialized ) then
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
              if( sd_r(n) > max_radius_in_height ) max_radius_in_height = sd_r(n)
              if( use_radius_upper_bound ) then
                if( sd_r(n) > tracking_radius_max ) then
                  upper_exceed_cnt = upper_exceed_cnt + 1
                else if( sd_r(n) >= near_upper_threshold ) then
                  near_upper_cnt = near_upper_cnt + 1
                end if
              end if
            end if
            if( sd_z(n) >= tracking_height_min .and. sd_z(n) <= tracking_height_max .and. &
                sd_r(n) >= tracking_radius_min .and. &
                ( .not. use_radius_upper_bound .or. sd_r(n) <= tracking_radius_max ) ) then
              candidate_cnt = candidate_cnt + 1
              if( sd_r(n) > r_max_eff ) r_max_eff = sd_r(n)
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
                  sd_r(n) >= tracking_radius_min .and. &
                  ( .not. use_radius_upper_bound .or. sd_r(n) <= tracking_radius_max ) ) then
                iz = int( (sd_z(n)-tracking_height_min) / z_span * real(tracking_nz_bin,kind=RP) ) + 1
                iz = min( tracking_nz_bin, max(1,iz) )
                if( tracking_nr_bin == 1 ) then
                  ir = 1
                else
                  if( log_r_span > 0.0_RP .and. sd_r(n) > r_min_eff ) then
                    ir = int( log(sd_r(n)/r_min_eff) / log_r_span * real(tracking_nr_bin,kind=RP) ) + 1
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
                  sd_r(n) >= tracking_radius_min .and. &
                  ( .not. use_radius_upper_bound .or. sd_r(n) <= tracking_radius_max ) ) then
                iz = int( (sd_z(n)-tracking_height_min) / z_span * real(tracking_nz_bin,kind=RP) ) + 1
                iz = min( tracking_nz_bin, max(1,iz) )
                if( tracking_nr_bin == 1 ) then
                  ir = 1
                else
                  if( log_r_span > 0.0_RP .and. sd_r(n) > r_min_eff ) then
                    ir = int( log(sd_r(n)/r_min_eff) / log_r_span * real(tracking_nr_bin,kind=RP) ) + 1
                  else
                    ir = 1
                  end if
                  ir = min( tracking_nr_bin, max(1,ir) )
                end if
                ibin = (iz-1)*tracking_nr_bin + ir
                needed_in_bin = bin_quota(ibin) - bin_selected(ibin)
                if( needed_in_bin > 0 .and. bin_remaining(ibin) > 0 ) then
                  call random_number(rand_tracking)
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
              call random_number(rand_tracking)
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
        do_track = ( sd_id(n) > INVALID_i4 .and. dm_id(n) > INVALID_i4 )
        if( do_track .and. max_tracked_sds > 0 ) then
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
    end if

    return
  end subroutine sdm_assign_tracking_subset
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sdm_outnetcdf(otime,sd_num,sd_numasl,sd_n,sd_liqice,sd_x,sd_y,sd_z,sd_r,sd_asl,sd_vz,sdi,sd_id,dm_id,if_coal,sdn_dmpnskip,filetag)
    use netcdf
    use scale_precision
    use scale_stdio
    use scale_time
    use scale_process, only: &
         mype => PRC_myrank
    use m_sdm_common, only: &
         i2, sdm_cold, STAT_LIQ, STAT_ICE, STAT_MIX, sdicedef, &
         INVALID_i4, backward_tracking_enable, coalescence_output_enable

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
    integer :: sd_x_id, sd_y_id, sd_z_id, sd_vz_id, sd_r_id, sd_asl_id, sd_n_id, sd_liqice_id,sd_id_id, domain_id, if_coal_id
    integer :: sdi_re_id, sdi_rp_id, sdi_rho_id, sdi_tf_id, sdi_mrime_id, sdi_nmono_id
    logical :: write_forward_tracking, write_backward_tracking, write_tracking, write_coal
    character(len=80) :: tracking_id_label

    integer,parameter :: nc_deflate_level = 1      ! NetCDF compression level {1,..,9}
    integer,parameter :: nc_deflate       = .true. ! turn on NetCDF compresion
    integer,parameter :: nc_shuffle       = .true. ! turn on NetCDF shuffle filter
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

    write_forward_tracking = forward_tracking_enable
    write_backward_tracking = backward_tracking_enable
    write_tracking = write_forward_tracking .or. write_backward_tracking
    write_coal = coalescence_output_enable
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

    if( write_coal ) then
      if_coal(1:sd_num) = 0
    end if

    return

  end subroutine sdm_outnetcdf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sdm_outnetcdf_hist(otime,sd_num,sd_numasl,sd_n,sd_liqice,sd_x,sd_y,sd_z,sd_r,sd_asl,sd_vz,sdi,sd_id,dm_id,if_coal,sdn_dmpnskip,filetag)
    use netcdf
    use scale_precision
    use scale_stdio
    use scale_time
    use scale_process, only: &
         mype => PRC_myrank, &
         PRC_MPIstop
    use m_sdm_common, only: &
         i2, sdm_cold, STAT_LIQ, STAT_ICE, STAT_MIX, sdicedef, &
         INVALID_i4, forward_tracking_enable, backward_tracking_enable, coalescence_output_enable

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
    integer :: sd_x_id, sd_y_id, sd_z_id, sd_vz_id, sd_r_id, sd_asl_id, sd_n_id, sd_liqice_id, sd_id_id, domain_id, if_coal_id
    integer :: sdi_re_id, sdi_rp_id, sdi_rho_id, sdi_tf_id, sdi_mrime_id, sdi_nmono_id
    integer :: sd_dmp_time_idx_id, sd_dmp_time_id

    integer,parameter :: nc_deflate_level = 1      ! NetCDF compression level {1,..,9}
    integer,parameter :: nc_deflate       = .true. ! turn on NetCDF compresion
    integer,parameter :: nc_shuffle       = .true. ! turn on NetCDF shuffle filter

    logical :: newfile
    character(len=100) :: var_name
    character(len=100) :: ftag ! =filetag or ''(default)
    integer, parameter :: max_filenum = 2
    integer, save :: filenum = 0
    character(len=100), save :: ftag_list(1:max_filenum)
    integer, save :: time_count(1:max_filenum)
    integer :: nf, fileid
    logical :: write_forward_tracking, write_backward_tracking, write_tracking, write_coal
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

    write_forward_tracking = forward_tracking_enable
    write_backward_tracking = backward_tracking_enable
    write_tracking = write_forward_tracking .or. write_backward_tracking
    write_coal = coalescence_output_enable

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

    if( backward_tracking_enable ) then
      call sdm_assign_tracking_subset(sd_num, sd_z, sd_r, sd_id, dm_id, if_coal)
    end if

    return

  end subroutine sdm_outnetcdf_hist
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
  subroutine sdm_coal_outnetcdf(otime, num_pair, sd_id1, sd_id2, dm_id1, dm_id2, num_col, sdr1_out, sdr2_out, sdn1_out, sdn2_out)
    use netcdf
    use scale_precision
    use scale_stdio
    use scale_process, only: &
         mype => PRC_myrank
    use m_sdm_common, only: &
         sdm_cold, sdm_dmpitvl, forward_tracking_enable, backward_tracking_enable

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

    character(len=17) :: fmt2="(A, '.', A, I*.*)"
    character(len=H_LONG) :: ftmp
    character(len=H_LONG) :: basename_sd_out
    character(len=19) :: basename_time
    integer :: nf90_real_precision
    integer :: ncid, event_dim_id, sdr1_out_id, sdr2_out_id, sdn1_out_id, sdn2_out_id
    integer :: sd_id1_id, sd_id2_id, dm_id1_id, dm_id2_id, num_col_id, event_time_id
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
    integer,parameter :: nc_deflate       = .true.
    integer,parameter :: nc_shuffle       = .true.
    character(len=100) :: ftag

    if( num_pair <= 0 ) return

    ftag = 'SD_coal_output'
    write_tracking_ids = forward_tracking_enable .or. backward_tracking_enable
    if( backward_tracking_enable ) then
      id1_name = 'pre_sdid1'
      id2_name = 'pre_sdid2'
      dm1_name = 'pre_dmid1'
      dm2_name = 'pre_dmid2'
    else if( forward_tracking_enable ) then
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

    ! NOTE:
    ! Unified collision-output strategy for both forward and backward tracking:
    ! append events into one file per SDM-output interval and MPI rank.
    if( sdm_dmpitvl > 0.0_RP ) then
       otime_bucket = sdm_dmpitvl * real( int( otime / sdm_dmpitvl + 1.0e-10_RP ), kind=DP )
    else
       otime_bucket = otime
    end if

    ! Keep naming style close to existing SD_coal_output_NetCDF_00000101-hhmmss.mmm.
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
    write(ftmp,'(3A)') trim(ftag), '_NetCDF_', trim(basename_time)
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
       call check_netcdf( nf90_inq_varid(ncid, "num_col", num_col_id) )
    else
       call check_netcdf( nf90_create(trim(basename_sd_out),NF90_NETCDF4, ncid) )
       call check_netcdf( nf90_def_dim(ncid, "event", NF90_UNLIMITED, event_dim_id) )
       call check_netcdf( nf90_def_var(ncid, "event_time", NF90_DOUBLE, event_dim_id, event_time_id) )
       call check_netcdf( nf90_def_var_deflate(ncid, event_time_id, &
            & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
       call check_netcdf( nf90_put_att(ncid, event_time_id, 'long_name', 'coalescence event time') )
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
       call check_netcdf( nf90_def_var(ncid, "num_col", NF90_INT, event_dim_id, num_col_id) )
       call check_netcdf( nf90_def_var_deflate(ncid, num_col_id, &
            & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
       call check_netcdf( nf90_put_att(ncid, num_col_id, 'long_name', 'coalescence count of SD pairs') )
       call check_netcdf( nf90_put_att(ncid, num_col_id, 'units', '') )
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
    call check_netcdf( nf90_put_var(ncid, num_col_id, num_col, start=start1, count=count1) )
    deallocate(event_time)

    call check_netcdf( nf90_close(ncid) )

    if( IO_L ) write(IO_FID_LOG,*) '*** Closed output file (NetCDF) of Super Droplet'

    return

  end subroutine sdm_coal_outnetcdf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module m_sdm_io
