!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Cloud Microphysics / SDM
!!
!! @par Description
!!          Utility subroutines dealing super-droplet IDs
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
!! @li      2014-07-12 (S.Shima) [new] Separated from scale_atmos_phy_mp_sdm.F90
!! @li      2020-07-23 (S.Shima) [add] sdm_copy_selected_sd for sdm_dmpvar == 1?? and sdm_dmpvar == 2??
!! @li      2020-07-24 (S.Shima) [add] 'all', 'large', and 'selected' options to sdm_copy_selected_sd 
!! @li      2020-07-30 (S.Shima) [fix] sdm_copy_selected_sd 
!! @li      2020-10-30 (S.Shima) [fix] sdm_copy_selected_sd 
!!
!<
!-------------------------------------------------------------------------------
module sdm_sorting_module
  use scale_precision
  implicit none
contains
  recursive subroutine quicksort(array, left, right, sd_rand, layer_indices)
    integer, intent(inout) :: array(:)
    integer, intent(in) :: left, right
    real(RP), intent(in) :: sd_rand(:)
    integer, intent(in) :: layer_indices(:)

    integer :: pivot_index

    if (left < right) then
      pivot_index = partition(array, left, right, sd_rand, layer_indices)
      call quicksort(array, left, pivot_index - 1, sd_rand, layer_indices)
      call quicksort(array, pivot_index + 1, right, sd_rand, layer_indices)
    end if
  end subroutine quicksort
!---------------------------------------------------------------------------------------------------------------------------------
  integer function partition(array, left, right, sd_rand, layer_indices)
    implicit none
    integer, intent(inout) :: array(:)
    integer, intent(in) :: left, right
    real(RP), intent(in) :: sd_rand(:)
    integer, intent(in) :: layer_indices(:)

    integer :: i, j, temp
    real(RP) :: pivot

    pivot = sd_rand(layer_indices(array(right)))
    i = left - 1

    do j = left, right - 1
      if (sd_rand(layer_indices(array(j))) <= pivot) then
        i = i + 1
        temp = array(i)
        array(i) = array(j)
        array(j) = temp
      end if
    end do

    temp = array(i + 1)
    array(i + 1) = array(right)
    array(right) = temp

    partition = i + 1
  end function partition
!---------------------------------------------------------------------------------------------------------------------------------
  subroutine sort_index_array_by_sd_rand(index_array, sd_rand, layer_indices, count)
    implicit none
    integer, intent(inout) :: index_array(:)
    real(RP), intent(in) :: sd_rand(:)
    integer, intent(in) :: layer_indices(:)
    integer, intent(in) :: count

    call quicksort(index_array, 1, count, sd_rand, layer_indices)
  end subroutine sort_index_array_by_sd_rand
end module sdm_sorting_module
!---------------------------------------------------------------------------------------------------------------------------------

module m_sdm_idutil
  use scale_precision

  implicit none
  private
  public :: sdm_sort,sdm_getperm,sdm_copy_selected_sd,sdm_select_stratified_random_particles, &
       sdm_select_particles_from_id_file

contains
  subroutine sdm_getperm(freq_max,ni_sdm,nj_sdm,nk_sdm,sd_num,    &
                         sort_tag0,fsort_tag,fsort_id,            &
                         sd_rand,sd_perm)

    ! Input variables
    integer, intent(in) :: freq_max ! maximum number of SD in each sd-grid
    integer, intent(in) :: ni_sdm   ! SDM model dimension in x direction
    integer, intent(in) :: nj_sdm   ! SDM model dimension in y direction
    integer, intent(in) :: nk_sdm   ! SDM model dimension in z direction
    integer, intent(in) :: sd_num   ! Number of super-droplets
    integer, intent(in) :: sort_tag0(1:ni_sdm*nj_sdm*nk_sdm+2)
    ! = sort_tag(n) - 1
    ! sort_tag(m) : accumulated number of
    !               super-droplets in each
    !               SDM-grid
    integer, intent(in) :: fsort_tag(0:freq_max+1)
    ! accumulated number with respect to
    ! the number contaiend super-droplets
    integer, intent(in) :: fsort_id(1:ni_sdm*nj_sdm*nk_sdm+2)
    ! super-droplets sorted by the number
    ! contaiend super-droplets
    ! Input and output variables
    real(RP), intent(inout) :: sd_rand(1:sd_num) ! random numbers
    ! Output variables
    integer, intent(out) :: sd_perm(1:sd_num)    ! random permutations
    ! Work variables
    integer :: m, n, t, ss, tt         ! index
    !---------------------------------------------------------------

    ! Initialize for calculating random permutations
    
    do t=fsort_tag(2),fsort_tag(freq_max+1)-1

       m = fsort_id(t)

       tt = sort_tag0(m) + 1
       sd_perm(tt) = 1

    end do

    ! Calculate random permutations

    do n=2,freq_max

       do t=fsort_tag(n),fsort_tag(freq_max+1)-1

          m = fsort_id(t)

          tt = sort_tag0(m) + n
          !ORG        ss = sort_tag0(m) + ( int(sd_rand(tt)*n) + 1 )
          ss = sort_tag0(m) + ( int(sd_rand(tt)*real(n,kind=RP)) + 1 )

          !### swap data ###!

          sd_perm(tt) = sd_perm(ss)
          sd_perm(ss) = n

       end do

    end do

    return
  end subroutine sdm_getperm
  !----------------------------------------------------------------------------
  subroutine sdm_sort(ni_sdm,nj_sdm,nk_sdm,                       &
                      sd_num,sd_n,sd_ri,sd_rj,sd_rk,                &
                      sort_id,sort_key,sort_freq,sort_tag,jdgtype)
    use scale_grid_index, only: &
         IS,JS,KS
    use m_sdm_common, only: &
         VALID2INVALID,knum_sdm
    use gadg_algorithm, only: &
         gadg_count_sort
    ! Input variables
    integer, intent(in) :: ni_sdm  ! SDM model dimension in x direction
    integer, intent(in) :: nj_sdm  ! SDM model dimension in y direction
    integer, intent(in) :: nk_sdm  ! SDM model dimension in z direction
    integer, intent(in) :: sd_num  ! number of super-droplets
    integer(DP), intent(in) :: sd_n(1:sd_num)   ! multiplicity of super-droplets
    real(RP), intent(in) :: sd_ri(1:sd_num)  ! index[i/real] of super-droplets
    real(RP), intent(in) :: sd_rj(1:sd_num)  ! index[j/real] of super-droplets
    real(RP), intent(in) :: sd_rk(1:sd_num)  ! index[k/real] of super-droplets
    character(len=5), intent(in) :: jdgtype   ! flag for sorting
    ! Output variables
    integer, intent(out) :: sort_id(1:sd_num)   ! id that super-droplets sorted by sd-grids
    integer, intent(out) :: sort_key(1:sd_num)  ! sort key
    integer, intent(inout) :: sort_freq(1:ni_sdm*nj_sdm*nk_sdm+1)  ! number of super-droplets in each sd-grid
    integer, intent(out) :: sort_tag(1:ni_sdm*nj_sdm*nk_sdm+2)     ! accumulated number of super-droplets in each sd-grid
    integer :: adr                ! grid id
    integer :: max_key            ! total grid number
    integer :: i, j, k, n         ! index
    integer :: ix, jy
    !--------------------------------------------------------------------

    ! Initialize

    max_key = ni_sdm * nj_sdm * knum_sdm + 1

    ! Sorting [step-1] -- initialize IKEY
    if( jdgtype .eq. 'valid' ) then

       do n=1,sd_num

          if( sd_rk(n)>VALID2INVALID ) then

             ! i \in {0,...,ni_sdm-1}, j \in {0,...,nj_sdm-1}, k \in {0,...,knum_sdm-1}
             i = floor(sd_ri(n)) - (IS-1) 
             j = floor(sd_rj(n)) - (JS-1)
             k = floor(sd_rk(n)) - (KS-1)

             adr = ni_sdm*nj_sdm*k + ni_sdm*j + i + 1

          else

             adr = max_key    !! invalid super-droplets

          end if

          sort_key(n) = adr

       end do

    else if( jdgtype .eq. 'multi' ) then

       do n=1,sd_num

          if( sd_rk(n)>VALID2INVALID .and. sd_n(n)>1 ) then

             i = floor(sd_ri(n)) - (IS-1) 
             j = floor(sd_rj(n)) - (JS-1)
             k = floor(sd_rk(n)) - (KS-1)

             adr = ni_sdm*nj_sdm*k + ni_sdm*j + i + 1

          else

             adr = max_key    !! invalid super-droplets

          end if

          sort_key(n) = adr

       end do

    end if

    
    ! Sorting [step-2] -- counting sort
    call gadg_count_sort( sort_key, 1, max_key,                       &
         sort_freq, sort_tag, sort_id )

    sort_tag(max_key+1) = sort_tag(max_key) + sort_freq(max_key)

    return
  end subroutine sdm_sort
  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine sdm_copy_selected_sd(sd_num,    sd_numasl,    sd_n,    sd_x,    sd_y,    sd_ri,    sd_rj,    sd_rk,     &
       &                          sd_liqice,    sd_asl,    sd_r,    sdi,     sd_id,   dm_id,    if_coal,             &
       &                          sd_num_tmp,sd_numasl_tmp,sd_n_tmp,sd_x_tmp,sd_y_tmp,sd_ri_tmp,sd_rj_tmp,sd_rk_tmp, &
       &                          sd_liqice_tmp,sd_asl_tmp,sd_r_tmp,sdi_tmp,sd_id_tmp,dm_id_tmp,if_coal_tmp,         &
       &                          TEMP0,ilist,sdtype)
    use scale_process, only: &
         & PRC_MPIstop
    use scale_grid_index, only: &
         & IA,JA,KA
    use m_sdm_common, only: &
         & i2, sdicedef, sdm_cold, num_threads, VALID2INVALID, STAT_LIQ, STAT_ICE, &
         & sdm_aslset, mass_amsul, ion_amsul, mass_nacl, ion_nacl, CurveF, ASL_FF, &
         & INVALID_i4
    use m_sdm_coordtrans, only: &
         & sdm_x2ri, sdm_y2rj

    integer,  intent(in)  :: sd_num      ! number of super-droplets
    integer,  intent(in)  :: sd_numasl   ! number of kind of chemical material contained as water-soluble aerosol in super droplets
    integer(DP), intent(in) :: sd_n(1:sd_num) ! multiplicity of super-droplets
    real(RP), intent(in)  :: sd_x(1:sd_num) ! x-coordinate of super-droplets
    real(RP), intent(in)  :: sd_y(1:sd_num) ! y-coordinate of super-droplets
    real(RP), intent(inout) :: sd_ri(1:sd_num)   ! index[i/real] of super-droplets
    real(RP), intent(inout) :: sd_rj(1:sd_num)   ! index[j/real] of super-droplets
    real(RP), intent(in)  :: sd_rk(1:sd_num) ! face index-k(real) of super-droplets
    integer(i2), intent(in) :: if_coal(1:sd_num)
                       ! flag of coalescence
                       ! 0 = Super Droplet hasn't undergone coalescence during the previous output interval
                       ! 1 = Super Droplet has undergone coalescence during the previous output interval
    integer(i2), intent(in) :: sd_liqice(1:sd_num)
                       ! status of super-droplets (liquid/ice)
                       ! 01 = all liquid, 10 = all ice
                       ! 11 = mixture of ice and liquid
    real(RP), intent(in) :: sd_asl(1:sd_num,1:sd_numasl) ! aerosol mass of super-droplets
    real(RP), intent(in) :: sd_r(1:sd_num) ! equivalent radius of super-droplets
    type(sdicedef), intent(in) :: sdi   ! ice phase super-droplets
    integer, intent(in) :: sd_id(1:sd_num)
    integer, intent(in) :: dm_id(1:sd_num)

    integer,  intent(out) :: sd_num_tmp  ! number of super-droplets
    integer,  intent(out)  :: sd_numasl_tmp   ! number of kind of chemical material contained as water-soluble aerosol in super droplets
    integer(DP), intent(out) :: sd_n_tmp(1:sd_num) ! multiplicity of super-droplets
    real(RP), intent(out)  :: sd_x_tmp(1:sd_num) ! x-coordinate of super-droplets
    real(RP), intent(out)  :: sd_y_tmp(1:sd_num) ! x-coordinate of super-droplets
    real(RP), intent(out) :: sd_ri_tmp(1:sd_num)   ! index[i/real] of super-droplets
    real(RP), intent(out) :: sd_rj_tmp(1:sd_num)   ! index[j/real] of super-droplets
    real(RP), intent(out)  :: sd_rk_tmp(1:sd_num) ! face index-k(real) of super-droplets
    integer(i2), intent(out) :: sd_liqice_tmp(1:sd_num)
    real(RP), intent(out) :: sd_asl_tmp(1:sd_num,1:sd_numasl) ! aerosol mass of super-droplets
    real(RP), intent(out) :: sd_r_tmp(1:sd_num) ! equivalent radius of super-droplets
    type(sdicedef), intent(inout) :: sdi_tmp   ! ice phase super-droplets
    integer, intent(out) :: sd_id_tmp(1:sd_num)
    integer, intent(out) :: dm_id_tmp(1:sd_num)
    integer(i2), intent(out) :: if_coal_tmp(1:sd_num)

    real(RP), intent(in)  :: TEMP0(KA,IA,JA)       ! temperature [K]

    integer, intent(out) :: ilist(1:sd_num)   ! buffer for index list
    character(len=*), intent(in) :: sdtype

    real(RP) :: sd_thld_radi ! threshold radius [m]
    integer :: n,k,cnt,m,i,j,t,s
    integer :: i_threads
    real(RP):: sd_aslmw(1:22) ! Molecular mass of chemical material contained as water-soluble aerosol in super droplets (default+20)
    real(RP):: sd_aslion(1:22) ! Degree of ion dissociation of chemical material contained as water-soluble aerosol in super droplets (default+20)
    integer :: idx_nasl(1:22)  ! index for vactorization
    real(RP) :: dmask(1:22)  ! mask for vactorization
    real(RP):: coef_a, coef_b ! Coefficients of Kohler curve
    real(RP) :: t_sd      ! temperature of the grid contained the SD
    real(RP) :: ivt_sd    ! 1.d0 / t_sd
    real(RP) :: dtmp      ! temporary for interation  

    call sdm_x2ri(sd_num,sd_x,sd_ri,sd_rk)
    call sdm_y2rj(sd_num,sd_y,sd_rj,sd_rk)
    
    !### Copy the same aerosol chemical componets
    sd_numasl_tmp = sd_numasl

    !### Setup aerosol related parameters
    if( abs(mod(sdm_aslset,10))==1 ) then

       !### numasl=1 @ init+rest : (NH4)2SO4 ###!

       sd_aslmw(1)  = mass_amsul
       sd_aslion(1) = ion_amsul

    else if( abs(mod(sdm_aslset,10))==2 ) then

       if( abs(sdm_aslset)==2 ) then

          !### numasl=1 @ init : NaCl ###!

          sd_aslmw(1)  = mass_nacl
          sd_aslion(1) = ion_nacl

       else if( abs(sdm_aslset)==12 ) then

          !### numasl=2 @ init : NaCl, rest : (NH4)2SO4 ###!

          sd_aslmw(1) = mass_amsul
          sd_aslmw(2) = mass_nacl
          sd_aslion(1) = ion_amsul
          sd_aslion(2) = ion_nacl

       end if

    else if( abs(mod(sdm_aslset,10))==3 ) then

       !### numasl>=2 @ init+rest : (NH4)2SO4, NaCl, ... ###!

       sd_aslmw(1) = mass_amsul
       sd_aslmw(2) = mass_nacl

       sd_aslion(1) = ion_amsul
       sd_aslion(2) = ion_nacl

       !! Must be a Bug. This cannot be simply commented out
       do n=1,20
       !            call getrname( id_sdm_aslmw  + (n-1), sd_aslmw(n+2)  )
       !            call getrname( id_sdm_aslion + (n-1), sd_aslion(n+2) )
       end do

    else if( abs(mod(sdm_aslset,10))==5 ) then

       !### numasl=1 @ init+rest : (NH4)HSO4 ###!

       sd_aslmw(1)  = mass_amsul
       sd_aslion(1) = ion_amsul

    end if

    do n=1,22

       if( n<=sd_numasl ) then
          idx_nasl(n) = n
          dmask(n) = 1.0_RP
       else
          idx_nasl(n) = sd_numasl
          dmask(n) = 0.0_RP
       end if

    end do


    ! Get index list of the selected SDs
    cnt = 0
    if(sdtype == 'all') then
       do n=1,sd_num
          cnt = cnt + 1
          ilist(cnt) = n
       end do

    else if (sdtype == 'large') then
       sd_thld_radi = 1.0e-6_RP ! threshold radius [m]
       do n=1,sd_num
          if( sd_rk(n)<VALID2INVALID ) cycle

          if( sd_liqice(n) == STAT_LIQ ) then
             if( sd_r(n)>sd_thld_radi ) then
                cnt = cnt + 1
                ilist(cnt) = n
             end if

          else if( sdm_cold .and. (sd_liqice(n) == STAT_ICE)) then
             if( (sdi%re(n)>sd_thld_radi) .or. (sdi%rp(n)>sd_thld_radi) ) then
                cnt = cnt + 1
                ilist(cnt) = n
             end if
          end if
       end do

    else if (sdtype == 'activated') then
       do n=1,sd_num
          if( sd_rk(n)<VALID2INVALID ) cycle

          if( sd_liqice(n) == STAT_LIQ ) then
             !! calculate the coefficient a of Kohler curve
             i = floor(sd_ri(n))+1
             j = floor(sd_rj(n))+1
             k = floor(sd_rk(n))+1

             t_sd  = TEMP0(k,i,j)

             ivt_sd = 1.0_RP / t_sd
             coef_a  = CurveF * ivt_sd

             !! calculate the coefficient b of Kohler curve
             coef_b = 0.0_RP

!OCL UNROLL('full'),NOSWP  
             do t=1,22

                s = idx_nasl(t)

                dtmp = sd_asl(n,s) * (real(sd_aslion(s),kind=RP)            &
                     / real(sd_aslmw(s),kind=RP))
                coef_b = coef_b + dmask(t) * dtmp
                
             end do

             coef_b = coef_b * ASL_FF

             !! calculate critical radius
             sd_thld_radi = sqrt(3.0_RP*coef_b/coef_a)

             if( sd_r(n)>sd_thld_radi ) then
                cnt = cnt + 1
                ilist(cnt) = n
             end if

          else if( sdm_cold .and. (sd_liqice(n) == STAT_ICE)) then
             sd_thld_radi = 1.0e-6_RP ! threshold radius [m]
             if( (sdi%re(n)>sd_thld_radi) .or. (sdi%rp(n)>sd_thld_radi) ) then
                cnt = cnt + 1
                ilist(cnt) = n
             end if
          end if
       end do

    else if (sdtype == 'selected') then
       do n=1,sd_num
          if( dm_id(n)<0 .or. sd_rk(n)<VALID2INVALID ) cycle

          cnt = cnt + 1
          ilist(cnt) = n

       end do

    else
       ! stop if unsupported sdtype option is specified
       write(*,*) "sdm_copy_selected_sd: Unsupported sdtype option is specified"
       call PRC_MPIstop

    end if
    sd_num_tmp = cnt

    ! Copy data of the selcted SDs
    if(sd_num_tmp /= 0) then
       do i_threads=1,num_threads
       do m=i_threads,sd_num_tmp,num_threads
          n = ilist(m)
       
          sd_n_tmp(m)      = sd_n(n)
          sd_x_tmp(m)      = sd_x(n)
          sd_y_tmp(m)      = sd_y(n)
          sd_ri_tmp(m)     = sd_ri(n)
          sd_rj_tmp(m)     = sd_rj(n)
          sd_rk_tmp(m)     = sd_rk(n)
          sd_liqice_tmp(m) = sd_liqice(n)
          sd_r_tmp(m)      = sd_r(n)
          sd_id_tmp(m)      = sd_id(n)
          dm_id_tmp(m)      = dm_id(n)
          if_coal_tmp(m)      = if_coal(n)

       end do
       end do

       do k=1,sd_numasl_tmp
          do i_threads=1,num_threads
          do m=i_threads,sd_num_tmp,num_threads
             n = ilist(m)

             sd_asl_tmp(m,k) = sd_asl(n,k)

          end do
          end do
       end do

       if( sdm_cold ) then
          do i_threads=1,num_threads
          do m=i_threads,sd_num_tmp,num_threads
             n = ilist(m)

             sdi_tmp%re(m) = sdi%re(n)
             sdi_tmp%rp(m) = sdi%rp(n)
             sdi_tmp%rho(m) = sdi%rho(n)
             sdi_tmp%tf(m) = sdi%tf(n)
             sdi_tmp%mrime(m) = sdi%mrime(n)
             sdi_tmp%nmono(m) = sdi%nmono(n)

          end do
          end do
       end if
    end if

  end subroutine sdm_copy_selected_sd
!---------------------------------------------------------------------------------------------------------------------------------
  subroutine sdm_select_particles_from_id_file(sd_num, tracking_id_input_basename, tracking_sample_initialized, sd_id, dm_id, if_coal, status_rdm)
    use netcdf
    use scale_stdio, only: &
         H_LONG, IO_L, IO_FID_LOG, IO_get_available_fid
    use scale_process, only: &
         mype => PRC_myrank
    use m_sdm_common, only: &
         i2, INVALID_i4

    implicit none

    integer, intent(in) :: sd_num
    character(len=*), intent(in) :: tracking_id_input_basename
    logical, intent(inout) :: tracking_sample_initialized
    integer, intent(inout) :: sd_id(1:sd_num)
    integer, intent(inout) :: dm_id(1:sd_num)
    integer(i2), intent(inout) :: if_coal(1:sd_num)
    integer, intent(out) :: status_rdm

    character(len=H_LONG) :: input_filename
    character(len=512) :: linebuf
    integer :: fid, ierr, n, m, pair_cnt, ios, ncid, dimid, varid
    integer, allocatable :: target_sd_id(:), target_dm_id(:)
    logical :: do_track, file_exists, is_netcdf

    status_rdm = 0
    input_filename = trim(adjustl(tracking_id_input_basename))
    is_netcdf = .false.

    if( len_trim(input_filename) == 0 ) then
      status_rdm = 1
      return
    end if

    if( len_trim(input_filename) >= 3 ) then
      if( input_filename(len_trim(input_filename)-2:len_trim(input_filename)) == '.nc' ) then
        is_netcdf = .true.
      end if
    end if

    if( .not. is_netcdf ) then
      if( len_trim(input_filename) >= 4 ) then
        if( input_filename(len_trim(input_filename)-3:len_trim(input_filename)) /= '.ids' ) then
          write(input_filename,'(A,".pe",I6.6,".nc")') trim(adjustl(tracking_id_input_basename)), mype
          inquire(file=trim(input_filename), exist=file_exists)
          if( file_exists ) then
            is_netcdf = .true.
          else
            write(input_filename,'(A,".pe",I6.6,".ids")') trim(adjustl(tracking_id_input_basename)), mype
          end if
        end if
      else
        write(input_filename,'(A,".pe",I6.6,".nc")') trim(adjustl(tracking_id_input_basename)), mype
        inquire(file=trim(input_filename), exist=file_exists)
        if( file_exists ) then
          is_netcdf = .true.
        else
          write(input_filename,'(A,".pe",I6.6,".ids")') trim(adjustl(tracking_id_input_basename)), mype
        end if
      end if
    end if

    if( is_netcdf ) then
      ierr = nf90_open(trim(input_filename), NF90_NOWRITE, ncid)
      if( ierr /= nf90_noerr ) then
        if( IO_L ) write(IO_FID_LOG,*) '*** WARNING (sdm_select_particles_from_id_file): failed to open tracking ID NetCDF file:', trim(input_filename)
        status_rdm = 2
        return
      end if

      ierr = nf90_inq_dimid(ncid, 'tracked_id', dimid)
      if( ierr /= nf90_noerr ) then
        ios = nf90_close(ncid)
        if( IO_L ) write(IO_FID_LOG,*) '*** WARNING (sdm_select_particles_from_id_file): missing tracked_id dimension:', trim(input_filename)
        status_rdm = 3
        return
      end if

      ierr = nf90_inquire_dimension(ncid, dimid, len=pair_cnt)
      if( ierr /= nf90_noerr ) then
        ios = nf90_close(ncid)
        if( IO_L ) write(IO_FID_LOG,*) '*** WARNING (sdm_select_particles_from_id_file): failed to inspect tracked_id dimension:', trim(input_filename)
        status_rdm = 3
        return
      end if

      if( pair_cnt <= 0 ) then
        ios = nf90_close(ncid)
        do n = 1, sd_num
          sd_id(n) = INVALID_i4
          dm_id(n) = INVALID_i4
          if_coal(n) = 0_i2
        end do
        tracking_sample_initialized = .true.
        return
      end if

      allocate(target_dm_id(pair_cnt), target_sd_id(pair_cnt))

      ierr = nf90_inq_varid(ncid, 'dm_id', varid)
      if( ierr == nf90_noerr ) ierr = nf90_get_var(ncid, varid, target_dm_id)
      if( ierr == nf90_noerr ) ierr = nf90_inq_varid(ncid, 'sd_id', varid)
      if( ierr == nf90_noerr ) ierr = nf90_get_var(ncid, varid, target_sd_id)
      ios = nf90_close(ncid)
      if( ierr /= nf90_noerr ) then
        if( IO_L ) write(IO_FID_LOG,*) '*** WARNING (sdm_select_particles_from_id_file): failed to read tracking IDs from NetCDF:', trim(input_filename)
        deallocate(target_dm_id, target_sd_id)
        status_rdm = 3
        return
      end if
    else
      fid = IO_get_available_fid()
      open(unit=fid, file=trim(input_filename), form='formatted', status='old', action='read', iostat=ierr)
      if( ierr /= 0 ) then
        if( IO_L ) write(IO_FID_LOG,*) '*** WARNING (sdm_select_particles_from_id_file): failed to open tracking ID file:', trim(input_filename)
        status_rdm = 2
        return
      end if

      pair_cnt = 0
      do
        read(fid,'(A)',iostat=ierr) linebuf
        if( ierr /= 0 ) exit
        if( len_trim(linebuf) == 0 ) cycle
        if( linebuf(1:1) == '#' ) cycle
        read(linebuf,*,iostat=ios) m, n
        if( ios /= 0 ) cycle
        pair_cnt = pair_cnt + 1
      end do

      rewind(fid)

      if( pair_cnt <= 0 ) then
        close(fid)
        do n = 1, sd_num
          sd_id(n) = INVALID_i4
          dm_id(n) = INVALID_i4
          if_coal(n) = 0_i2
        end do
        tracking_sample_initialized = .true.
        return
      end if

      allocate(target_dm_id(pair_cnt), target_sd_id(pair_cnt))
      pair_cnt = 0
      do
        read(fid,'(A)',iostat=ierr) linebuf
        if( ierr /= 0 ) exit
        if( len_trim(linebuf) == 0 ) cycle
        if( linebuf(1:1) == '#' ) cycle
        read(linebuf,*,iostat=ios) m, n
        if( ios /= 0 ) cycle
        pair_cnt = pair_cnt + 1
        target_dm_id(pair_cnt) = m
        target_sd_id(pair_cnt) = n
      end do
      close(fid)
    end if

    do n = 1, sd_num
      do_track = .false.
      if( sd_id(n) > INVALID_i4 .and. dm_id(n) > INVALID_i4 ) then
        do m = 1, pair_cnt
          if( dm_id(n) == target_dm_id(m) .and. sd_id(n) == target_sd_id(m) ) then
            do_track = .true.
            exit
          end if
        end do
      end if
      if( .not. do_track ) then
        sd_id(n) = INVALID_i4
        dm_id(n) = INVALID_i4
      end if
      if_coal(n) = 0_i2
    end do

    tracking_sample_initialized = .true.
    deallocate(target_dm_id, target_sd_id)

    return
  end subroutine sdm_select_particles_from_id_file
!---------------------------------------------------------------------------------------------------------------------------------
  subroutine sdm_select_stratified_random_particles(sd_num, sd_rk, sd_r,         &
                                                    tracking_selection_mode,      &
                                                    tracking_fraction,             &
                                                    max_tracked_sds,              &
                                                    tracking_height_min,          &
                                                    tracking_height_max,          &
                                                    tracking_radius_min,          &
                                                    tracking_radius_max,          &
                                                    tracking_nz_bin,              &
                                                    tracking_nr_bin,              &
                                                    tracking_min_per_bin,         &
                                                    tracking_fallback_to_random,  &
                                                    tracking_sample_initialized,  &
                                                    dm_id, sd_id, if_coal, status_rdm)
    use scale_precision
    use scale_grid, only: DZ
    use m_sdm_common, only: VALID2INVALID, INVALID_i4, i2
    use scale_process, only: mype => PRC_myrank

    implicit none

    integer, intent(in) :: sd_num
    real(RP), intent(in) :: sd_rk(1:sd_num)
    real(RP), intent(in) :: sd_r(1:sd_num)
    character(len=*), intent(in) :: tracking_selection_mode
    real(RP), intent(in) :: tracking_fraction
    integer, intent(in) :: max_tracked_sds
    real(RP), intent(in) :: tracking_height_min
    real(RP), intent(in) :: tracking_height_max
    real(RP), intent(in) :: tracking_radius_min
    real(RP), intent(in) :: tracking_radius_max
    integer, intent(in) :: tracking_nz_bin
    integer, intent(in) :: tracking_nr_bin
    integer, intent(in) :: tracking_min_per_bin
    logical, intent(in) :: tracking_fallback_to_random
    logical, intent(inout) :: tracking_sample_initialized
    integer, intent(inout) :: sd_id(1:sd_num)
    integer, intent(inout) :: dm_id(1:sd_num)
    integer(i2), intent(inout) :: if_coal(1:sd_num)
    integer, intent(out) :: status_rdm

    integer :: n, iz, ir, ibin
    integer :: tracked_cnt, candidate_cnt, target_cnt
    integer :: nbin, sum_quota, extra_needed, reduce_needed
    integer :: needed_in_bin, non_empty_bins, min_per_bin_eff
    real(RP) :: rand_tracking
    real(RP) :: z_span, r_min_eff, r_max_eff, log_r_span, z_height
    logical :: do_track, use_stratified, stratified_selected, use_radius_upper_bound
    integer, allocatable :: bin_cnt(:), bin_quota(:), bin_selected(:), bin_remaining(:), bin_min(:)
    real(RP), allocatable :: bin_frac(:)

    status_rdm = 0

    if( tracking_fraction <= 0.0_RP ) then
      do n = 1, sd_num
        sd_id(n) = INVALID_i4
        dm_id(n) = INVALID_i4
        if_coal(n) = 0_i2
      end do
      tracking_sample_initialized = .false.
      return
    end if

    tracked_cnt = 0
    if( .not. tracking_sample_initialized ) then
      use_stratified = trim(adjustl(tracking_selection_mode)) == 'stratified' .or. &
                       trim(adjustl(tracking_selection_mode)) == 'STRATIFIED' .or. &
                       trim(adjustl(tracking_selection_mode)) == 'Stratified'
      use_radius_upper_bound = tracking_radius_max > tracking_radius_min
      stratified_selected = .false.

      if( use_stratified .and. tracking_nz_bin > 0 .and. tracking_nr_bin > 0 .and. &
          tracking_height_max > tracking_height_min ) then
        nbin = tracking_nz_bin * tracking_nr_bin
        allocate(bin_cnt(nbin), bin_quota(nbin), bin_selected(nbin), bin_remaining(nbin), bin_min(nbin), bin_frac(nbin))

        candidate_cnt = 0
        r_min_eff = max(tracking_radius_min, 1.0E-12_RP)
        r_max_eff = r_min_eff
        do n = 1, sd_num
          if( sd_rk(n) <= VALID2INVALID ) cycle
          z_height = sd_rk(n) * DZ
          if( z_height >= tracking_height_min .and. z_height <= tracking_height_max .and. &
              sd_r(n) >= tracking_radius_min .and. &
              ( .not. use_radius_upper_bound .or. sd_r(n) <= tracking_radius_max ) ) then
            candidate_cnt = candidate_cnt + 1
            if( sd_r(n) > r_max_eff ) r_max_eff = sd_r(n)
          end if
        end do

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

          do n = 1, sd_num
            if( sd_rk(n) <= VALID2INVALID ) cycle
            z_height = sd_rk(n) * DZ
            if( z_height >= tracking_height_min .and. z_height <= tracking_height_max .and. &
                sd_r(n) >= tracking_radius_min .and. &
                ( .not. use_radius_upper_bound .or. sd_r(n) <= tracking_radius_max ) ) then
              iz = int( (z_height-tracking_height_min) / z_span * real(tracking_nz_bin,kind=RP) ) + 1
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

          do ibin = 1, nbin
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
              do ibin = 1, nbin
                if( extra_needed <= 0 ) exit
                if( bin_quota(ibin) < bin_cnt(ibin) ) then
                  bin_quota(ibin) = bin_quota(ibin) + 1
                  extra_needed = extra_needed - 1
                end if
              end do
              if( all(bin_quota >= bin_cnt) ) exit
            end do
          else if( sum_quota > target_cnt ) then
            reduce_needed = sum_quota - target_cnt
            do while( reduce_needed > 0 )
              do ibin = nbin, 1, -1
                if( reduce_needed <= 0 ) exit
                if( bin_quota(ibin) > bin_min(ibin) ) then
                  bin_quota(ibin) = bin_quota(ibin) - 1
                  reduce_needed = reduce_needed - 1
                end if
              end do
              if( all(bin_quota <= bin_min) ) exit
            end do
          end if

          bin_remaining(:) = bin_cnt(:)
          do n = 1, sd_num
            do_track = .false.
            if( sd_rk(n) > VALID2INVALID ) then
              z_height = sd_rk(n) * DZ
              if( z_height >= tracking_height_min .and. z_height <= tracking_height_max .and. &
                  sd_r(n) >= tracking_radius_min .and. &
                  ( .not. use_radius_upper_bound .or. sd_r(n) <= tracking_radius_max ) ) then
                iz = int( (z_height-tracking_height_min) / z_span * real(tracking_nz_bin,kind=RP) ) + 1
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

      if( .not. stratified_selected ) then
        if( use_stratified .and. .not. tracking_fallback_to_random ) then
          do n = 1, sd_num
            sd_id(n) = INVALID_i4
            dm_id(n) = INVALID_i4
            if_coal(n) = 0_i2
          end do
          tracking_sample_initialized = .true.
        else
          do n = 1, sd_num
            do_track = ( sd_rk(n) > VALID2INVALID )
            if( do_track .and. tracking_fraction < 1.0_RP ) then
              call random_number(rand_tracking)
              if( rand_tracking > tracking_fraction ) do_track = .false.
            end if
            if( do_track .and. max_tracked_sds > 0 .and. tracked_cnt >= max_tracked_sds ) do_track = .false.
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
      do n = 1, sd_num
        if( sd_id(n) > INVALID_i4 .and. dm_id(n) > INVALID_i4 ) then
          tracked_cnt = tracked_cnt + 1
        end if
        if_coal(n) = 0_i2
      end do
    end if

    return
  end subroutine sdm_select_stratified_random_particles

end module m_sdm_idutil
