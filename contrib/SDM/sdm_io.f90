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
!! @li      2023-03-01 (C.Yin)   [add] output the index of Super-droplets
!!
!<
!-------------------------------------------------------------------------------
module m_sdm_io

  implicit none
  private
  public :: sdm_outasci,sdm_outnetcdf,sdm_outnetcdf_hist,sdm_coal_outnetcdf

contains
  subroutine sdm_outasci(otime,sd_num,sd_numasl,sd_n,sd_liqice,sd_x,sd_y,sd_z,sd_r,sd_asl,sd_vz,sdi,pre_sdid,pre_dmid,sdn_dmpnskip)
    use scale_precision
    use scale_stdio
    use scale_time
    use scale_process, only: &
         mype => PRC_myrank, &
         PRC_MPIstop
    use m_sdm_common, only: &
         i2, sdm_cold, STAT_LIQ, STAT_ICE, STAT_MIX, sdicedef

    implicit none

    real(DP), intent(in) :: otime
    integer, intent(in) :: sd_num    ! number of super-droplets
    integer, intent(in) :: sd_numasl ! number of chemical species contained in super droplets
    integer(DP), intent(in) :: sd_n(1:sd_num) ! multiplicity of super-droplets
    integer, intent(inout) :: pre_sdid(1:sd_num) ! save index of super-droplets
    integer, intent(inout) :: pre_dmid(1:sd_num) ! domain index of super-droplets
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
            &                'multiplicity[-],status[-],index,previous index,previous domain'

       write(fmt,'( "(", i2.2, "e16.8,i20,a5,i10)" )')(sd_numasl+5)

       do m=1,sd_num,sdn_dmpnskip

          cstat = '  LIQ'
          
          write(fid_sdm_o,trim(fmt)) sd_x(m),           &
               &                   sd_y(m),           &
               &                   sd_z(m), sd_vz(m), sd_r(m),        &
               &                   (sd_asl(m,n),n=1,sd_numasl),       &
               &                   sd_n(m), cstat, m, pre_sdid(m), pre_dmid(m)
          pre_sdid(m) = m
          pre_dmid(m) = mype
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
            &            'multiplicity[-],status[-],index,rime_mass[kg],num_of_monomers[-],', &
            &            'previous index,previous domain'

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
               &                   sd_n(m), cstat, m, sdi%mrime(m), sdi%nmono(m), &
               &                   pre_sdid(m), pre_dmid(m)
               pre_sdid(m) = m
               pre_dmid(m) = mype
       end do


    end if

    close(fid_sdm_o)
    if( IO_L ) write(IO_FID_LOG,*) '*** Closed output file (ASCII) of Super Droplet'

    return

  end subroutine sdm_outasci
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sdm_outnetcdf(otime,sd_num,sd_numasl,sd_n,sd_liqice,sd_x,sd_y,sd_z,sd_r,sd_asl,sd_vz,sdi,pre_sdid,pre_dmid,sdn_dmpnskip,filetag)
    use netcdf
    use scale_precision
    use scale_stdio
    use scale_time
    use scale_process, only: &
         mype => PRC_myrank
    use m_sdm_common, only: &
         i2, sdm_cold, STAT_LIQ, STAT_ICE, STAT_MIX, sdicedef

    implicit none

    real(DP), intent(in) :: otime
    integer, intent(in) :: sd_num    ! number of super-droplets
    integer, intent(in) :: sd_numasl ! number of chemical species contained in super droplets
    integer(DP), intent(in) :: sd_n(1:sd_num) ! multiplicity of super-droplets
    integer, intent(inout) :: pre_sdid(1:sd_num) ! save index of super-droplets
    integer, intent(inout) :: pre_dmid(1:sd_num) ! domain index of super-droplets
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
    integer :: sd_x_id, sd_y_id, sd_z_id, sd_vz_id, sd_r_id, sd_asl_id, sd_n_id, sd_liqice_id, sd_id_id, domain_id
    integer :: sdi_re_id, sdi_rp_id, sdi_rho_id, sdi_tf_id, sdi_mrime_id, sdi_nmono_id

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
    !!! pre_sdid
    call check_netcdf( nf90_def_var(ncid, "pre_sdid", NF90_INT, sd_num_id, sd_id_id) )
    call check_netcdf( nf90_def_var_deflate(ncid, sd_id_id, &
         & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
    call check_netcdf( nf90_put_att(ncid, sd_id_id, 'long_name', 'previous index') )
    call check_netcdf( nf90_put_att(ncid, sd_id_id, 'units', '') )
    !!! pre_dmid
    call check_netcdf( nf90_def_var(ncid, "pre_dmid", NF90_INT, sd_num_id, domain_id) )
    call check_netcdf( nf90_def_var_deflate(ncid, domain_id, &
         & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
    call check_netcdf( nf90_put_att(ncid, domain_id, 'long_name', 'previous domain') )
    call check_netcdf( nf90_put_att(ncid, domain_id, 'units', '') )

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
    !!! pre_sdid
    call check_netcdf( nf90_put_var(ncid, sd_id_id, pre_sdid) )
    !!! pre_dmid
    call check_netcdf( nf90_put_var(ncid, domain_id, pre_dmid) )

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

    do n=1,sd_num
       pre_sdid(n) = n
       pre_dmid(n) = mype
    end do

    return

  end subroutine sdm_outnetcdf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sdm_outnetcdf_hist(otime,sd_num,sd_numasl,sd_n,sd_liqice,sd_x,sd_y,sd_z,sd_r,sd_asl,sd_vz,sdi,pre_sdid,pre_dmid,sdn_dmpnskip,filetag)
    use netcdf
    use scale_precision
    use scale_stdio
    use scale_time
    use scale_process, only: &
         mype => PRC_myrank, &
         PRC_MPIstop
    use m_sdm_common, only: &
         i2, sdm_cold, STAT_LIQ, STAT_ICE, STAT_MIX, sdicedef

    implicit none

    real(DP), intent(in) :: otime
    integer, intent(in) :: sd_num    ! number of super-droplets
    integer, intent(in) :: sd_numasl ! number of chemical species contained in super droplets
    integer(DP), intent(in) :: sd_n(1:sd_num) ! multiplicity of super-droplets
    integer, intent(inout) :: pre_sdid(1:sd_num) ! save index of super-droplets
    integer, intent(inout) :: pre_dmid(1:sd_num) ! domain index of super-droplets
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
    integer :: sd_x_id, sd_y_id, sd_z_id, sd_vz_id, sd_r_id, sd_asl_id, sd_n_id, sd_liqice_id, sd_id_id, domain_id
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
    !!! pre_sdid
    call check_netcdf( nf90_def_var(ncid, "pre_sdid", NF90_INT, sd_num_id, sd_id_id) )
    call check_netcdf( nf90_def_var_deflate(ncid, sd_id_id, &
         & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
    call check_netcdf( nf90_put_att(ncid, sd_id_id, 'long_name', 'previous index') )
    call check_netcdf( nf90_put_att(ncid, sd_id_id, 'units', '') )
    !!! pre_dmid
    call check_netcdf( nf90_def_var(ncid, "pre_dmid", NF90_INT, sd_num_id, domain_id) )
    call check_netcdf( nf90_def_var_deflate(ncid, domain_id, &
         & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
    call check_netcdf( nf90_put_att(ncid, domain_id, 'long_name', 'previous domain') )
    call check_netcdf( nf90_put_att(ncid, domain_id, 'units', '') )

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
    !!! pre_sdid
    call check_netcdf( nf90_put_var(ncid, sd_id_id, pre_sdid) )
    !!! pre_dmid
    call check_netcdf( nf90_put_var(ncid, domain_id, pre_dmid) )

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

    do n=1,sd_num
       pre_sdid(n) = n
       pre_dmid(n) = mype
    end do

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
  subroutine sdm_coal_outnetcdf(otime, num_pair, pre_sdid1, pre_sdid2, pre_dmid1, pre_dmid2, num_col)
    use netcdf
    use scale_precision
    use scale_stdio
    use scale_time
    use scale_process, only: &
         mype => PRC_myrank
    use m_sdm_common, only: &
         sdm_cold

    implicit none

    real(DP), intent(in) :: otime
    integer, intent(in) :: num_pair    ! number of super-droplet pairs
    integer, allocatable, intent(in) :: pre_sdid1(:) ! save index of super-droplets
    integer, allocatable, intent(in) :: pre_dmid1(:) ! domain index of super-droplets
    integer, allocatable, intent(in) :: pre_sdid2(:) ! save index of super-droplets
    integer, allocatable, intent(in) :: pre_dmid2(:) ! domain index of super-droplets
    integer, allocatable, intent(in) :: num_col(:)   ! number of coalesence of pairs of SDs

    character(len=17) :: fmt2="(A, '.', A, I*.*)"
    character(len=H_LONG) :: ftmp
    character(len=H_LONG) :: basename_sd_out
    character(len=19) :: basename_time
    integer :: fid_sdm_o
    integer :: ierr
    integer :: nf90_real_precision
    integer :: ncid, num_pair_id
    integer :: pre_sdid1_id, pre_sdid2_id, pre_dmid1_id, pre_dmid2_id, num_col_id

    integer,parameter :: nc_deflate_level = 1      ! NetCDF compression level {1,..,9}
    integer,parameter :: nc_deflate       = .true. ! turn on NetCDF compresion
    integer,parameter :: nc_shuffle       = .true. ! turn on NetCDF shuffle filter
    character(len=100) :: ftag ! =filetag or ''(default)

    !--- output Super Droplets in NetCDF format
    call TIME_gettimelabel(basename_time)
    fid_sdm_o = IO_get_available_fid()

    ftag = 'SD_coal_output'

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
    !!! num_pair
    call check_netcdf( nf90_def_dim(ncid, "num_pair", num_pair, num_pair_id) )

    !!! pre_sdid1
    call check_netcdf( nf90_def_var(ncid, "pre_sdid1", NF90_INT, num_pair_id,pre_sdid1_id) )
    call check_netcdf( nf90_def_var_deflate(ncid, pre_sdid1_id, &
         & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
    call check_netcdf( nf90_put_att(ncid, pre_sdid1_id, 'long_name', 'previous index of large SD') )
    call check_netcdf( nf90_put_att(ncid, pre_sdid1_id, 'units', '') )
    !!! pre_dmid1
    call check_netcdf( nf90_def_var(ncid, "pre_dmid1", NF90_INT, num_pair_id, pre_dmid1_id) )
    call check_netcdf( nf90_def_var_deflate(ncid, pre_dmid1_id, &
         & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
    call check_netcdf( nf90_put_att(ncid, pre_dmid1_id, 'long_name', 'previous domain of large SD') )
    call check_netcdf( nf90_put_att(ncid, pre_dmid1_id, 'units', '') )
    !!! pre_sdid2
    call check_netcdf( nf90_def_var(ncid, "pre_sdid2", NF90_INT, num_pair_id,pre_sdid2_id) )
    call check_netcdf( nf90_def_var_deflate(ncid, pre_sdid2_id, &
         & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
    call check_netcdf( nf90_put_att(ncid, pre_sdid2_id, 'long_name', 'previous index of small SD') )
    call check_netcdf( nf90_put_att(ncid, pre_sdid2_id, 'units', '') )
    !!! pre_dmid2
    call check_netcdf( nf90_def_var(ncid, "pre_dmid2", NF90_INT, num_pair_id, pre_dmid2_id) )
    call check_netcdf( nf90_def_var_deflate(ncid, pre_dmid2_id, &
         & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
    call check_netcdf( nf90_put_att(ncid, pre_dmid2_id, 'long_name', 'previous domain of small SD') )
    call check_netcdf( nf90_put_att(ncid, pre_dmid2_id, 'units', '') )
    !!! num_col
    call check_netcdf( nf90_def_var(ncid, "num_col", NF90_INT, num_pair_id, num_col_id) )
    call check_netcdf( nf90_def_var_deflate(ncid, num_col_id, &
         & shuffle=nc_shuffle, deflate=nc_deflate, deflate_level=nc_deflate_level) )
    call check_netcdf( nf90_put_att(ncid, num_col_id, 'long_name', 'coalescence count of SD pairs') )
    call check_netcdf( nf90_put_att(ncid, num_col_id, 'units', '') )

    !if( sdm_cold ) then
        !!! to be continued
    !end if

    !!! End of definition
    call check_netcdf( nf90_enddef(ncid) )

    ! Save data
    !!! pre_sdid1
    call check_netcdf( nf90_put_var(ncid, pre_sdid1_id, pre_sdid1) )
    !!! pre_dmid1
    call check_netcdf( nf90_put_var(ncid, pre_dmid1_id, pre_dmid1) )
    !!! pre_sdid2
    call check_netcdf( nf90_put_var(ncid, pre_sdid2_id, pre_sdid2) )
    !!! pre_dmid2
    call check_netcdf( nf90_put_var(ncid, pre_dmid2_id, pre_dmid2) )
    !!! num_col
    call check_netcdf( nf90_put_var(ncid, num_col_id, num_col) )

    !if( sdm_cold ) then
       !!! to be continued
    !end if

    ! Close the output file
    call check_netcdf( nf90_close(ncid) )

    if( IO_L ) write(IO_FID_LOG,*) '*** Closed output file (NetCDF) of Super Droplet'

    return

  end subroutine sdm_coal_outnetcdf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module m_sdm_io
