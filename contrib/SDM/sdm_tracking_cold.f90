!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Cloud Microphysics / SDM cold tracking helpers
!!
!! @par Description
!!          Constants and small utilities for cold-SDM event tracking.
!<
!-------------------------------------------------------------------------------
module m_sdm_tracking_cold
  use scale_precision
  use scale_const, only: &
       ONE_PI => CONST_PI, &
       dens_w_mks => CONST_DWATR
  use m_sdm_common, only: &
       i2, sdicedef, STAT_LIQ, STAT_ICE, STAT_MIX, VALID2INVALID, &
       CurveF, ASL_FF, &
       mass_amsul, ion_amsul, mass_nacl, ion_nacl, &
       TRACK_ID_INVALID, TRACK_SD_ID_DYNAMIC_START, &
       tracking_diag_liq_radius_enable, tracking_diag_liq_radius_threshold, &
       tracking_diag_ice_rvol_enable, tracking_diag_ice_rvol_threshold, &
       tracking_diag_mixed_rvol_enable, tracking_diag_mixed_rvol_threshold, &
       tracking_diag_rime_mass_enable, tracking_diag_rime_mass_threshold, &
       tracking_diag_rime_frac_enable, tracking_diag_rime_frac_threshold, &
       tracking_diag_nmono_enable, tracking_diag_nmono_threshold, &
       tracking_diag_aspect_ratio_enable, tracking_diag_aspect_ratio_threshold, &
       tracking_spatial_visit_enable, tracking_spatial_active_rank, &
       tracking_spatial_x_min, tracking_spatial_x_max, &
       tracking_spatial_y_min, tracking_spatial_y_max, &
       tracking_spatial_z_min, tracking_spatial_z_max

  implicit none
  private

  public :: EVENT_COALESCENCE, EVENT_LIQ_LIQ_COAL, EVENT_RIMING, EVENT_AGGREGATION
  public :: EVENT_FREEZING, EVENT_MELTING, EVENT_DEPOSITION, EVENT_SUBLIMATION
  public :: EVENT_CONDENSATION, EVENT_EVAPORATION, EVENT_ACTIVATION, EVENT_DEACTIVATION
  public :: EVENT_SIG_DEPOSITION, EVENT_SIG_SUBLIMATION
  public :: EVENT_SIG_CONDENSATION, EVENT_SIG_EVAPORATION
  public :: EVENT_SIG_FREEZING, EVENT_SIG_MELTING
  public :: EVENT_SIG_RIMING, EVENT_SIG_AGGREGATION, EVENT_SIG_LIQ_LIQ_COAL
  public :: DIAG_LIQ_RADIUS_LARGE, DIAG_ICE_RVOL_LARGE, DIAG_MIXED_RVOL_LARGE
  public :: DIAG_RIME_MASS_LARGE, DIAG_RIME_FRAC_LARGE, DIAG_NMONO_LARGE
  public :: DIAG_ASPECT_RATIO_MATCH
  public :: TRIG_PROC_COALESCENCE, TRIG_PROC_RIMING, TRIG_PROC_AGGREGATION
  public :: TRIG_PROC_FREEZING, TRIG_PROC_MELTING, TRIG_PROC_DEPOSITION, TRIG_PROC_SUBLIMATION
  public :: TRIG_PROC_CONDENSATION, TRIG_PROC_EVAPORATION, TRIG_PROC_ACTIVATION, TRIG_PROC_DEACTIVATION
  public :: TRIG_LEVEL_OCCURRENCE, TRIG_LEVEL_SIGNIFICANT
  public :: TRIG_EVT_LIQ_LIQ_COAL, TRIG_EVT_RIMING, TRIG_EVT_AGGREGATION
  public :: TRIG_EVT_FREEZING, TRIG_EVT_MELTING, TRIG_EVT_SIG_DEPOSITION, TRIG_EVT_SIG_SUBLIMATION
  public :: TRIG_EVT_SIG_CONDENSATION, TRIG_EVT_SIG_EVAPORATION
  public :: TRIG_EVT_SIG_FREEZING, TRIG_EVT_SIG_MELTING
  public :: TRIG_EVT_SIG_RIMING, TRIG_EVT_SIG_AGGREGATION, TRIG_EVT_SIG_LIQ_LIQ_COAL
  public :: TRIG_EVT_DEPOSITION, TRIG_EVT_SUBLIMATION, TRIG_EVT_CONDENSATION, TRIG_EVT_EVAPORATION
  public :: TRIG_RESERVED_HOMOGENEOUS_FREEZING, TRIG_RESERVED_IMMERSION_FREEZING
  public :: TRIG_RESERVED_DEPOSITION_FREEZING, TRIG_RESERVED_CONDENSATION_FREEZING
  public :: TRIG_RESERVED_RIME_SPLINTERING, TRIG_RESERVED_ICE_COLLISIONAL_BREAKUP
  public :: TRIG_RESERVED_SPONTANEOUS_ICE_BREAKUP
  public :: TRIG_RESERVED_DELIQUESCENCE, TRIG_RESERVED_EFFLORESCENCE
  public :: TRIG_DIAG_LIQ_RADIUS_LARGE, TRIG_DIAG_ICE_RVOL_LARGE
  public :: TRIG_DIAG_MIXED_RVOL_LARGE, TRIG_DIAG_RIME_MASS_LARGE
  public :: TRIG_DIAG_RIME_FRAC_LARGE, TRIG_DIAG_NMONO_LARGE
  public :: TRIG_DIAG_ASPECT_RATIO_MATCH
  public :: TARGET_BY_EVENT, TARGET_BY_DIAG, TARGET_BY_TRANSITION
  public :: PHASE_DRY_AEROSOL, PHASE_LIQUID, PHASE_ICE, PHASE_MIXED, PHASE_NONE
  public :: sdm_cold_phase_from_liqice, sdm_cold_phase_name
  public :: sdm_cold_phase_state, sdm_cold_collision_trigger
  public :: sdm_cold_event_bit_from_trigger, sdm_cold_diag_trigger_from_bit
  public :: sdm_cold_liq_mass, sdm_cold_ice_mass, sdm_cold_hydro_mass
  public :: sdm_cold_ice_rvol, sdm_cold_mixed_rvol, sdm_cold_hydro_radius
  public :: sdm_cold_rime_fraction, sdm_cold_aspect_ratio
  public :: sdm_kohler_solute_params
  public :: sdm_kohler_critical_radius
  public :: sdm_kohler_activated_state
  public :: sdm_cold_diag_mask_from_values
  public :: sdm_cold_tracking_update_interval
  public :: sdm_cold_tracking_update_spatial_visit
  public :: sdm_cold_tracking_update_spatial_visit_segment
  public :: sdm_cold_segment_intersects_box
  public :: sdm_cold_tracking_reset_interval
  public :: sdm_cold_tracking_evaluate_diag_masks
  public :: TRACK_ID_INVALID, TRACK_SD_ID_DYNAMIC_START
  public :: sdm_tracking_valid_sd_id, sdm_tracking_invalid_sd_id
  public :: sdm_tracking_valid_dm_id, sdm_tracking_invalid_dm_id
  public :: sdm_tracking_valid_id_pair
  public :: sdm_tracking_count_valid_id_pairs
  public :: sdm_tracking_assign_dynamic_id, sdm_tracking_find_slot_by_id
  public :: sdm_tracking_reset_invalid_slot
  public :: LIFE_SDADD_SPLIT, LIFE_SDREMOVE_INVALIDATION
  public :: LIFE_DOMAIN_ENTRY, LIFE_GLOBAL_HALO_ENTRY, LIFE_SEEDING_ENTRY, LIFE_ASLFORM_NEW_SD

  ! sd_event_mask and sd_event_sig_mask share this process bit table.
  ! sd_event_mask marks a tracked occurrence or significant hit; sd_event_sig_mask
  ! marks the significant subset.  Bits are interval summaries, not counters.
  integer, parameter :: EVENT_COALESCENCE   = 1
  integer, parameter :: EVENT_LIQ_LIQ_COAL  = EVENT_COALESCENCE ! deprecated alias
  integer, parameter :: EVENT_RIMING        = 2
  integer, parameter :: EVENT_AGGREGATION   = 4
  integer, parameter :: EVENT_FREEZING      = 8
  integer, parameter :: EVENT_MELTING       = 16
  integer, parameter :: EVENT_DEPOSITION    = 32
  integer, parameter :: EVENT_SUBLIMATION   = 64
  integer, parameter :: EVENT_CONDENSATION  = 128
  integer, parameter :: EVENT_EVAPORATION   = 256
  integer, parameter :: EVENT_ACTIVATION    = 512
  integer, parameter :: EVENT_DEACTIVATION  = 1024

  ! Deprecated v1.2 significant aliases.  Do not use these for current output.
  integer, parameter :: EVENT_SIG_DEPOSITION   = EVENT_DEPOSITION
  integer, parameter :: EVENT_SIG_SUBLIMATION  = EVENT_SUBLIMATION
  integer, parameter :: EVENT_SIG_CONDENSATION = EVENT_CONDENSATION
  integer, parameter :: EVENT_SIG_EVAPORATION  = EVENT_EVAPORATION
  integer, parameter :: EVENT_SIG_FREEZING     = EVENT_FREEZING
  integer, parameter :: EVENT_SIG_MELTING      = EVENT_MELTING
  integer, parameter :: EVENT_SIG_RIMING       = EVENT_RIMING
  integer, parameter :: EVENT_SIG_AGGREGATION  = EVENT_AGGREGATION
  integer, parameter :: EVENT_SIG_LIQ_LIQ_COAL = EVENT_COALESCENCE

  integer, parameter :: DIAG_LIQ_RADIUS_LARGE    = 1
  integer, parameter :: DIAG_ICE_RVOL_LARGE      = 2
  integer, parameter :: DIAG_MIXED_RVOL_LARGE    = 4
  integer, parameter :: DIAG_RIME_MASS_LARGE     = 8
  integer, parameter :: DIAG_RIME_FRAC_LARGE     = 16
  integer, parameter :: DIAG_NMONO_LARGE         = 32
  integer, parameter :: DIAG_ASPECT_RATIO_MATCH  = 64

  ! For process event files, trigger_code is process identity and trigger_level
  ! distinguishes occurrence from significant.  Diagnostic trigger_code remains
  ! 101-107 and has no trigger_level.
  integer, parameter :: TRIG_PROC_COALESCENCE   = 1
  integer, parameter :: TRIG_PROC_RIMING        = 2
  integer, parameter :: TRIG_PROC_AGGREGATION   = 3
  integer, parameter :: TRIG_PROC_FREEZING      = 4
  integer, parameter :: TRIG_PROC_MELTING       = 5
  integer, parameter :: TRIG_PROC_DEPOSITION    = 6
  integer, parameter :: TRIG_PROC_SUBLIMATION   = 7
  integer, parameter :: TRIG_PROC_CONDENSATION  = 8
  integer, parameter :: TRIG_PROC_EVAPORATION   = 9
  integer, parameter :: TRIG_PROC_ACTIVATION    = 10
  integer, parameter :: TRIG_PROC_DEACTIVATION  = 11

  integer, parameter :: TRIG_LEVEL_OCCURRENCE  = 1
  integer, parameter :: TRIG_LEVEL_SIGNIFICANT = 2

  ! Deprecated compatibility aliases.  Significant-vs-occurrence is no longer
  ! encoded in trigger_code for active v1.3 output.
  integer, parameter :: TRIG_EVT_LIQ_LIQ_COAL     = TRIG_PROC_COALESCENCE
  integer, parameter :: TRIG_EVT_RIMING           = TRIG_PROC_RIMING
  integer, parameter :: TRIG_EVT_AGGREGATION      = TRIG_PROC_AGGREGATION
  integer, parameter :: TRIG_EVT_FREEZING         = TRIG_PROC_FREEZING
  integer, parameter :: TRIG_EVT_MELTING          = TRIG_PROC_MELTING
  integer, parameter :: TRIG_EVT_SIG_DEPOSITION   = TRIG_PROC_DEPOSITION
  integer, parameter :: TRIG_EVT_SIG_SUBLIMATION  = TRIG_PROC_SUBLIMATION
  integer, parameter :: TRIG_EVT_SIG_CONDENSATION = TRIG_PROC_CONDENSATION
  integer, parameter :: TRIG_EVT_SIG_EVAPORATION  = TRIG_PROC_EVAPORATION
  integer, parameter :: TRIG_EVT_SIG_FREEZING     = TRIG_PROC_FREEZING
  integer, parameter :: TRIG_EVT_SIG_MELTING      = TRIG_PROC_MELTING
  integer, parameter :: TRIG_EVT_SIG_RIMING       = TRIG_PROC_RIMING
  integer, parameter :: TRIG_EVT_SIG_AGGREGATION  = TRIG_PROC_AGGREGATION
  integer, parameter :: TRIG_EVT_SIG_LIQ_LIQ_COAL = TRIG_PROC_COALESCENCE
  integer, parameter :: TRIG_EVT_DEPOSITION       = TRIG_PROC_DEPOSITION
  integer, parameter :: TRIG_EVT_SUBLIMATION      = TRIG_PROC_SUBLIMATION
  integer, parameter :: TRIG_EVT_CONDENSATION     = TRIG_PROC_CONDENSATION
  integer, parameter :: TRIG_EVT_EVAPORATION      = TRIG_PROC_EVAPORATION

  ! Reserved trigger ranges:
  ! 20-29 primary ice nucleation source candidates
  ! 30-39 secondary ice / breakup / fragmentation candidates
  ! 40-49 aerosol-water phase candidates
  ! Reserved constants are inactive unless explicitly wired to a process hook.
  integer, parameter :: TRIG_RESERVED_HOMOGENEOUS_FREEZING    = 20
  integer, parameter :: TRIG_RESERVED_IMMERSION_FREEZING      = 21
  integer, parameter :: TRIG_RESERVED_DEPOSITION_FREEZING     = 22
  integer, parameter :: TRIG_RESERVED_CONDENSATION_FREEZING   = 23
  integer, parameter :: TRIG_RESERVED_RIME_SPLINTERING        = 30
  integer, parameter :: TRIG_RESERVED_ICE_COLLISIONAL_BREAKUP = 31
  integer, parameter :: TRIG_RESERVED_SPONTANEOUS_ICE_BREAKUP = 32
  integer, parameter :: TRIG_RESERVED_DELIQUESCENCE           = 40
  integer, parameter :: TRIG_RESERVED_EFFLORESCENCE           = 41

  integer, parameter :: TRIG_DIAG_LIQ_RADIUS_LARGE   = 101
  integer, parameter :: TRIG_DIAG_ICE_RVOL_LARGE     = 102
  integer, parameter :: TRIG_DIAG_MIXED_RVOL_LARGE   = 103
  integer, parameter :: TRIG_DIAG_RIME_MASS_LARGE    = 104
  integer, parameter :: TRIG_DIAG_RIME_FRAC_LARGE    = 105
  integer, parameter :: TRIG_DIAG_NMONO_LARGE        = 106
  integer, parameter :: TRIG_DIAG_ASPECT_RATIO_MATCH = 107

  integer, parameter :: TARGET_BY_EVENT      = 1
  integer, parameter :: TARGET_BY_DIAG       = 2
  integer, parameter :: TARGET_BY_TRANSITION = 4

  ! phase_state_pre/post uses the model sd_liqice convention where possible.
  ! 0 is dry/aerosol-only; 99 is the missing or not-applicable sentinel.
  integer, parameter :: PHASE_DRY_AEROSOL = 0
  integer, parameter :: PHASE_LIQUID      = 1
  integer, parameter :: PHASE_ICE         = 10
  integer, parameter :: PHASE_MIXED       = 11
  integer, parameter :: PHASE_NONE        = 99

  integer, parameter :: LIFE_SDADD_SPLIT = 1
  integer, parameter :: LIFE_SDREMOVE_INVALIDATION = 2
  integer, parameter :: LIFE_DOMAIN_ENTRY = 3
  integer, parameter :: LIFE_GLOBAL_HALO_ENTRY = 4
  integer, parameter :: LIFE_SEEDING_ENTRY = 5
  integer, parameter :: LIFE_ASLFORM_NEW_SD = 6

contains
  subroutine sdm_kohler_solute_params(sdm_aslset, user_aslmw, user_aslion, sd_aslmw, sd_aslion)
    integer, intent(in) :: sdm_aslset
    real(RP), intent(in) :: user_aslmw(:)
    real(RP), intent(in) :: user_aslion(:)
    real(RP), intent(out) :: sd_aslmw(:)
    real(RP), intent(out) :: sd_aslion(:)

    integer :: n, nuser, nout

    sd_aslmw(:) = 0.0_RP
    sd_aslion(:) = 0.0_RP

    select case(abs(mod(sdm_aslset,10)))
    case(1,5)
       if( size(sd_aslmw) >= 1 .and. size(sd_aslion) >= 1 ) then
          sd_aslmw(1)  = mass_amsul
          sd_aslion(1) = ion_amsul
       end if
    case(2)
       if( abs(sdm_aslset) == 2 ) then
          if( size(sd_aslmw) >= 1 .and. size(sd_aslion) >= 1 ) then
             sd_aslmw(1)  = mass_nacl
             sd_aslion(1) = ion_nacl
          end if
       else if( abs(sdm_aslset) == 12 ) then
          if( size(sd_aslmw) >= 2 .and. size(sd_aslion) >= 2 ) then
             sd_aslmw(1)  = mass_amsul
             sd_aslmw(2)  = mass_nacl
             sd_aslion(1) = ion_amsul
             sd_aslion(2) = ion_nacl
          end if
       end if
    case(3)
       if( size(sd_aslmw) >= 2 .and. size(sd_aslion) >= 2 ) then
          sd_aslmw(1)  = mass_amsul
          sd_aslmw(2)  = mass_nacl
          sd_aslion(1) = ion_amsul
          sd_aslion(2) = ion_nacl
       end if
       nuser = min(size(user_aslmw), size(user_aslion))
       nout = min(size(sd_aslmw), size(sd_aslion))
       do n = 1, min(nuser, max(0, nout - 2))
          if( user_aslmw(n) > 0.0_RP ) then
             sd_aslmw(n+2) = user_aslmw(n)
             sd_aslion(n+2) = user_aslion(n)
          end if
       end do
    end select
  end subroutine sdm_kohler_solute_params

  real(RP) function sdm_kohler_critical_radius(sd_numasl, sd_asl, sd_aslmw, sd_aslion, temperature) result(rcrit)
    integer, intent(in) :: sd_numasl
    real(RP), intent(in) :: sd_asl(:)
    real(RP), intent(in) :: sd_aslmw(:)
    real(RP), intent(in) :: sd_aslion(:)
    real(RP), intent(in) :: temperature

    integer :: s, ncomp
    real(RP) :: coef_a, coef_b

    rcrit = huge(1.0_RP)
    if( temperature <= 0.0_RP ) return

    coef_a = CurveF / temperature
    if( coef_a <= 0.0_RP ) return

    coef_b = 0.0_RP
    ncomp = min(sd_numasl, size(sd_asl), size(sd_aslmw), size(sd_aslion))
    do s = 1, ncomp
       if( sd_aslmw(s) > 0.0_RP ) then
          coef_b = coef_b + sd_asl(s) * (sd_aslion(s) / sd_aslmw(s))
       end if
    end do
    coef_b = coef_b * ASL_FF

    if( coef_b <= 0.0_RP ) return

    ! This is the selected-output sdtype='activated' critical radius.
    ! It is intentionally distinct from sdm_condevp's internal Rc=sqrt(eq_b/eq_a).
    rcrit = sqrt(3.0_RP * coef_b / coef_a)
  end function sdm_kohler_critical_radius

  logical function sdm_kohler_activated_state(valid, sd_liqice, sd_r, rcrit, eps) result(is_activated)
    logical, intent(in) :: valid
    integer(i2), intent(in) :: sd_liqice
    real(RP), intent(in) :: sd_r
    real(RP), intent(in) :: rcrit
    real(RP), intent(in) :: eps

    is_activated = valid .and. sd_liqice == STAT_LIQ .and. &
         sd_r > rcrit * (1.0_RP + eps)
  end function sdm_kohler_activated_state

  logical function sdm_tracking_valid_sd_id(id) result(is_valid)
    integer, intent(in) :: id

    is_valid = (id >= 0) .or. (id <= TRACK_SD_ID_DYNAMIC_START)
  end function sdm_tracking_valid_sd_id

  logical function sdm_tracking_invalid_sd_id(id) result(is_invalid)
    integer, intent(in) :: id

    is_invalid = (id == TRACK_ID_INVALID)
  end function sdm_tracking_invalid_sd_id

  logical function sdm_tracking_valid_dm_id(id) result(is_valid)
    integer, intent(in) :: id

    is_valid = (id >= 0)
  end function sdm_tracking_valid_dm_id

  logical function sdm_tracking_invalid_dm_id(id) result(is_invalid)
    integer, intent(in) :: id

    is_invalid = (id == TRACK_ID_INVALID)
  end function sdm_tracking_invalid_dm_id

  logical function sdm_tracking_valid_id_pair(sd_id, dm_id) result(is_valid)
    integer, intent(in) :: sd_id
    integer, intent(in) :: dm_id

    is_valid = sdm_tracking_valid_sd_id(sd_id) .and. sdm_tracking_valid_dm_id(dm_id)
  end function sdm_tracking_valid_id_pair

  integer function sdm_tracking_count_valid_id_pairs(sd_num, sd_id, dm_id) result(valid_count)
    integer, intent(in) :: sd_num
    integer, intent(in) :: sd_id(1:sd_num)
    integer, intent(in) :: dm_id(1:sd_num)
    integer :: n

    valid_count = 0
    do n = 1, sd_num
       if( sdm_tracking_valid_id_pair(sd_id(n), dm_id(n)) ) valid_count = valid_count + 1
    end do
  end function sdm_tracking_count_valid_id_pairs

  subroutine sdm_tracking_assign_dynamic_id(n, sd_id, dm_id, next_dynamic_sd_id, mype)
    integer, intent(in) :: n
    integer, intent(inout) :: sd_id(:)
    integer, intent(inout) :: dm_id(:)
    integer, intent(inout) :: next_dynamic_sd_id
    integer, intent(in) :: mype

    sd_id(n) = next_dynamic_sd_id
    dm_id(n) = mype
    next_dynamic_sd_id = next_dynamic_sd_id - 1
  end subroutine sdm_tracking_assign_dynamic_id

  subroutine sdm_tracking_find_slot_by_id(target_sd_id, target_dm_id, sd_id, dm_id, sd_num, found, slot)
    integer, intent(in) :: target_sd_id
    integer, intent(in) :: target_dm_id
    integer, intent(in) :: sd_num
    integer, intent(in) :: sd_id(1:sd_num)
    integer, intent(in) :: dm_id(1:sd_num)
    logical, intent(out) :: found
    integer, intent(out) :: slot
    integer :: n

    found = .false.
    slot = 0
    if( .not. sdm_tracking_valid_id_pair(target_sd_id, target_dm_id) ) return

    if( target_sd_id >= 1 .and. target_sd_id <= sd_num ) then
       if( sd_id(target_sd_id) == target_sd_id .and. dm_id(target_sd_id) == target_dm_id ) then
          found = .true.
          slot = target_sd_id
          return
       end if
    end if

    do n = 1, sd_num
       if( .not. sdm_tracking_valid_id_pair(sd_id(n), dm_id(n)) ) cycle
       if( sd_id(n) == target_sd_id .and. dm_id(n) == target_dm_id ) then
          found = .true.
          slot = n
          return
       end if
    end do
  end subroutine sdm_tracking_find_slot_by_id

  subroutine sdm_tracking_reset_invalid_slot(n, sd_id, dm_id, &
       sd_event_mask, sd_event_sig_mask, sd_diag_mask, sd_phase_change_flag, sd_spatial_visit_flag, &
       sd_liq_radius_max_interval, sd_ice_rvol_max_interval, &
       sd_mixed_rvol_max_interval, sd_rime_mass_max_interval, &
       sd_rime_frac_max_interval, sd_nmono_max_interval, sd_aspect_ratio_max_interval)
    integer, intent(in) :: n
    integer, intent(inout) :: sd_id(:)
    integer, intent(inout) :: dm_id(:)
    integer, intent(inout) :: sd_event_mask(:)
    integer, intent(inout) :: sd_event_sig_mask(:)
    integer, intent(inout) :: sd_diag_mask(:)
    integer, intent(inout) :: sd_phase_change_flag(:)
    integer, intent(inout) :: sd_spatial_visit_flag(:)
    real(RP), intent(inout) :: sd_liq_radius_max_interval(:)
    real(RP), intent(inout) :: sd_ice_rvol_max_interval(:)
    real(RP), intent(inout) :: sd_mixed_rvol_max_interval(:)
    real(RP), intent(inout) :: sd_rime_mass_max_interval(:)
    real(RP), intent(inout) :: sd_rime_frac_max_interval(:)
    real(RP), intent(inout) :: sd_nmono_max_interval(:)
    real(RP), intent(inout) :: sd_aspect_ratio_max_interval(:)

    ! Only approved lifecycle/reuse paths should call this after event/MPI capture.
    sd_id(n) = TRACK_ID_INVALID
    dm_id(n) = TRACK_ID_INVALID
    sd_event_mask(n) = 0
    sd_event_sig_mask(n) = 0
    sd_diag_mask(n) = 0
    sd_phase_change_flag(n) = 0
    sd_spatial_visit_flag(n) = 0
    sd_liq_radius_max_interval(n) = 0.0_RP
    sd_ice_rvol_max_interval(n) = 0.0_RP
    sd_mixed_rvol_max_interval(n) = 0.0_RP
    sd_rime_mass_max_interval(n) = 0.0_RP
    sd_rime_frac_max_interval(n) = 0.0_RP
    sd_nmono_max_interval(n) = 0.0_RP
    sd_aspect_ratio_max_interval(n) = 0.0_RP
  end subroutine sdm_tracking_reset_invalid_slot

  integer function sdm_cold_phase_from_liqice(sd_liqice_value) result(phase_state)
    integer(i2), intent(in) :: sd_liqice_value

    select case(sd_liqice_value)
    case(0_i2)
       phase_state = PHASE_DRY_AEROSOL
    case(STAT_LIQ)
       phase_state = PHASE_LIQUID
    case(STAT_ICE)
       phase_state = PHASE_ICE
    case(STAT_MIX)
       phase_state = PHASE_MIXED
    case default
       phase_state = PHASE_NONE
    end select
  end function sdm_cold_phase_from_liqice

  character(len=20) function sdm_cold_phase_name(phase_state) result(phase_name)
    integer, intent(in) :: phase_state

    select case(phase_state)
    case(PHASE_DRY_AEROSOL)
       phase_name = 'PHASE_DRY_AEROSOL'
    case(PHASE_LIQUID)
       phase_name = 'PHASE_LIQUID'
    case(PHASE_ICE)
       phase_name = 'PHASE_ICE'
    case(PHASE_MIXED)
       phase_name = 'PHASE_MIXED'
    case default
       phase_name = 'PHASE_NONE'
    end select
  end function sdm_cold_phase_name

  integer function sdm_cold_phase_state(sd_liqice, sd_r, re, rp_ice) result(phase_state)
    integer(i2), intent(in) :: sd_liqice
    real(RP), intent(in) :: sd_r
    real(RP), intent(in) :: re
    real(RP), intent(in) :: rp_ice

    phase_state = sdm_cold_phase_from_liqice(sd_liqice)
  end function sdm_cold_phase_state

  integer function sdm_cold_collision_trigger(sd_liqice1, sd_liqice2) result(trigger_code)
    integer(i2), intent(in) :: sd_liqice1
    integer(i2), intent(in) :: sd_liqice2

    if( (sd_liqice1 == STAT_LIQ) .and. (sd_liqice2 == STAT_LIQ) ) then
       trigger_code = TRIG_EVT_LIQ_LIQ_COAL
    else if( ((sd_liqice1 == STAT_ICE) .and. (sd_liqice2 == STAT_LIQ)) .or. &
         &   ((sd_liqice1 == STAT_LIQ) .and. (sd_liqice2 == STAT_ICE)) ) then
       trigger_code = TRIG_EVT_RIMING
    else if( (sd_liqice1 == STAT_ICE) .and. (sd_liqice2 == STAT_ICE) ) then
       trigger_code = TRIG_EVT_AGGREGATION
    else
       trigger_code = 0
    end if
  end function sdm_cold_collision_trigger

  integer function sdm_cold_event_bit_from_trigger(trigger_code) result(event_bit)
    integer, intent(in) :: trigger_code

    select case(trigger_code)
    case(TRIG_EVT_LIQ_LIQ_COAL);    event_bit = EVENT_LIQ_LIQ_COAL
    case(TRIG_EVT_RIMING);          event_bit = EVENT_RIMING
    case(TRIG_EVT_AGGREGATION);     event_bit = EVENT_AGGREGATION
    case(TRIG_EVT_FREEZING);        event_bit = EVENT_FREEZING
    case(TRIG_EVT_MELTING);         event_bit = EVENT_MELTING
    case(TRIG_PROC_DEPOSITION);     event_bit = EVENT_DEPOSITION
    case(TRIG_PROC_SUBLIMATION);    event_bit = EVENT_SUBLIMATION
    case(TRIG_PROC_CONDENSATION);   event_bit = EVENT_CONDENSATION
    case(TRIG_PROC_EVAPORATION);    event_bit = EVENT_EVAPORATION
    case(TRIG_PROC_ACTIVATION);     event_bit = EVENT_ACTIVATION
    case(TRIG_PROC_DEACTIVATION);   event_bit = EVENT_DEACTIVATION
    case default;                   event_bit = 0
    end select
  end function sdm_cold_event_bit_from_trigger

  integer function sdm_cold_diag_trigger_from_bit(diag_bit) result(trigger_code)
    integer, intent(in) :: diag_bit

    select case(diag_bit)
    case(DIAG_LIQ_RADIUS_LARGE);    trigger_code = TRIG_DIAG_LIQ_RADIUS_LARGE
    case(DIAG_ICE_RVOL_LARGE);      trigger_code = TRIG_DIAG_ICE_RVOL_LARGE
    case(DIAG_MIXED_RVOL_LARGE);    trigger_code = TRIG_DIAG_MIXED_RVOL_LARGE
    case(DIAG_RIME_MASS_LARGE);     trigger_code = TRIG_DIAG_RIME_MASS_LARGE
    case(DIAG_RIME_FRAC_LARGE);     trigger_code = TRIG_DIAG_RIME_FRAC_LARGE
    case(DIAG_NMONO_LARGE);         trigger_code = TRIG_DIAG_NMONO_LARGE
    case(DIAG_ASPECT_RATIO_MATCH);  trigger_code = TRIG_DIAG_ASPECT_RATIO_MATCH
    case default;                   trigger_code = 0
    end select
  end function sdm_cold_diag_trigger_from_bit

  real(RP) function sdm_cold_liq_mass(sd_r) result(liq_mass)
    real(RP), intent(in) :: sd_r

    if( sd_r > 0.0_RP ) then
       liq_mass = (4.0_RP/3.0_RP) * ONE_PI * sd_r**3 * dens_w_mks
    else
       liq_mass = 0.0_RP
    end if
  end function sdm_cold_liq_mass

  real(RP) function sdm_cold_ice_mass(re, rp_ice, rho) result(ice_mass)
    real(RP), intent(in) :: re
    real(RP), intent(in) :: rp_ice
    real(RP), intent(in) :: rho

    if( re > 0.0_RP .and. rp_ice > 0.0_RP .and. rho > 0.0_RP ) then
       ice_mass = (4.0_RP/3.0_RP) * ONE_PI * re**2 * rp_ice * rho
    else
       ice_mass = 0.0_RP
    end if
  end function sdm_cold_ice_mass

  real(RP) function sdm_cold_hydro_mass(phase_state, sd_r, re, rp_ice, rho) result(hydro_mass)
    integer, intent(in) :: phase_state
    real(RP), intent(in) :: sd_r
    real(RP), intent(in) :: re
    real(RP), intent(in) :: rp_ice
    real(RP), intent(in) :: rho

    select case(phase_state)
    case(PHASE_LIQUID)
       hydro_mass = sdm_cold_liq_mass(sd_r)
    case(PHASE_ICE)
       hydro_mass = sdm_cold_ice_mass(re, rp_ice, rho)
    case(PHASE_MIXED)
       hydro_mass = sdm_cold_liq_mass(sd_r) + sdm_cold_ice_mass(re, rp_ice, rho)
    case default
       hydro_mass = 0.0_RP
    end select
  end function sdm_cold_hydro_mass

  real(RP) function sdm_cold_ice_rvol(re, rp_ice) result(radius)
    real(RP), intent(in) :: re
    real(RP), intent(in) :: rp_ice

    if( re > 0.0_RP .and. rp_ice > 0.0_RP ) then
       radius = (re**2 * rp_ice)**(1.0_RP/3.0_RP)
    else
       radius = 0.0_RP
    end if
  end function sdm_cold_ice_rvol

  real(RP) function sdm_cold_mixed_rvol(sd_r, re, rp_ice) result(radius)
    real(RP), intent(in) :: sd_r
    real(RP), intent(in) :: re
    real(RP), intent(in) :: rp_ice

    radius = (max(sd_r, 0.0_RP)**3 + max(re, 0.0_RP)**2 * max(rp_ice, 0.0_RP))**(1.0_RP/3.0_RP)
  end function sdm_cold_mixed_rvol

  real(RP) function sdm_cold_hydro_radius(phase_state, sd_r, re, rp_ice) result(radius)
    integer, intent(in) :: phase_state
    real(RP), intent(in) :: sd_r
    real(RP), intent(in) :: re
    real(RP), intent(in) :: rp_ice

    ! Use a phase-aware volume-equivalent radius for cold hydrometeors.
    select case(phase_state)
    case(PHASE_LIQUID)
       radius = max(sd_r, 0.0_RP)
    case(PHASE_ICE)
       radius = sdm_cold_ice_rvol(re, rp_ice)
    case(PHASE_MIXED)
       radius = sdm_cold_mixed_rvol(sd_r, re, rp_ice)
    case default
       radius = 0.0_RP
    end select
  end function sdm_cold_hydro_radius

  real(RP) function sdm_cold_rime_fraction(mrime, re, rp_ice, rho) result(frac)
    real(RP), intent(in) :: mrime
    real(RP), intent(in) :: re
    real(RP), intent(in) :: rp_ice
    real(RP), intent(in) :: rho
    real(RP) :: ice_mass

    ice_mass = sdm_cold_ice_mass(re, rp_ice, rho)
    if( ice_mass > 0.0_RP ) then
       frac = max(mrime, 0.0_RP) / ice_mass
    else
       frac = 0.0_RP
    end if
  end function sdm_cold_rime_fraction

  real(RP) function sdm_cold_aspect_ratio(re, rp_ice) result(aspect_ratio)
    real(RP), intent(in) :: re
    real(RP), intent(in) :: rp_ice

    if( re > 0.0_RP .and. rp_ice > 0.0_RP ) then
       aspect_ratio = max(re, rp_ice) / min(re, rp_ice)
    else
       aspect_ratio = 0.0_RP
    end if
  end function sdm_cold_aspect_ratio

  integer function sdm_cold_diag_mask_from_values( &
       liq_radius_max, ice_rvol_max, mixed_rvol_max, rime_mass_max, &
       rime_frac_max, nmono_max, aspect_ratio_max, &
       liq_radius_enable, liq_radius_threshold, &
       ice_rvol_enable, ice_rvol_threshold, &
       mixed_rvol_enable, mixed_rvol_threshold, &
       rime_mass_enable, rime_mass_threshold, &
       rime_frac_enable, rime_frac_threshold, &
       nmono_enable, nmono_threshold, &
       aspect_ratio_enable, aspect_ratio_threshold) result(diag_mask)
    real(RP), intent(in) :: liq_radius_max
    real(RP), intent(in) :: ice_rvol_max
    real(RP), intent(in) :: mixed_rvol_max
    real(RP), intent(in) :: rime_mass_max
    real(RP), intent(in) :: rime_frac_max
    real(RP), intent(in) :: nmono_max
    real(RP), intent(in) :: aspect_ratio_max
    logical, intent(in) :: liq_radius_enable
    real(RP), intent(in) :: liq_radius_threshold
    logical, intent(in) :: ice_rvol_enable
    real(RP), intent(in) :: ice_rvol_threshold
    logical, intent(in) :: mixed_rvol_enable
    real(RP), intent(in) :: mixed_rvol_threshold
    logical, intent(in) :: rime_mass_enable
    real(RP), intent(in) :: rime_mass_threshold
    logical, intent(in) :: rime_frac_enable
    real(RP), intent(in) :: rime_frac_threshold
    logical, intent(in) :: nmono_enable
    real(RP), intent(in) :: nmono_threshold
    logical, intent(in) :: aspect_ratio_enable
    real(RP), intent(in) :: aspect_ratio_threshold

    diag_mask = 0
    if( liq_radius_enable .and. liq_radius_max >= liq_radius_threshold ) &
         diag_mask = diag_mask + DIAG_LIQ_RADIUS_LARGE
    if( ice_rvol_enable .and. ice_rvol_max >= ice_rvol_threshold ) &
         diag_mask = diag_mask + DIAG_ICE_RVOL_LARGE
    if( mixed_rvol_enable .and. mixed_rvol_max >= mixed_rvol_threshold ) &
         diag_mask = diag_mask + DIAG_MIXED_RVOL_LARGE
    if( rime_mass_enable .and. rime_mass_max >= rime_mass_threshold ) &
         diag_mask = diag_mask + DIAG_RIME_MASS_LARGE
    if( rime_frac_enable .and. rime_frac_max >= rime_frac_threshold ) &
         diag_mask = diag_mask + DIAG_RIME_FRAC_LARGE
    if( nmono_enable .and. nmono_max >= nmono_threshold ) &
         diag_mask = diag_mask + DIAG_NMONO_LARGE
    if( aspect_ratio_enable .and. aspect_ratio_max >= aspect_ratio_threshold ) &
         diag_mask = diag_mask + DIAG_ASPECT_RATIO_MATCH
  end function sdm_cold_diag_mask_from_values

  subroutine sdm_cold_tracking_update_interval( &
       sd_num, sd_liqice, sd_r, sdi, &
       liq_radius_max, ice_rvol_max, mixed_rvol_max, &
       rime_mass_max, rime_frac_max, nmono_max, aspect_ratio_max)
    integer, intent(in) :: sd_num
    integer(i2), intent(in) :: sd_liqice(1:sd_num)
    real(RP), intent(in) :: sd_r(1:sd_num)
    type(sdicedef), intent(in) :: sdi
    real(RP), intent(inout) :: liq_radius_max(1:sd_num)
    real(RP), intent(inout) :: ice_rvol_max(1:sd_num)
    real(RP), intent(inout) :: mixed_rvol_max(1:sd_num)
    real(RP), intent(inout) :: rime_mass_max(1:sd_num)
    real(RP), intent(inout) :: rime_frac_max(1:sd_num)
    real(RP), intent(inout) :: nmono_max(1:sd_num)
    real(RP), intent(inout) :: aspect_ratio_max(1:sd_num)
    integer :: n
    integer :: phase_state

    ! Maintain maxima until the next SD output/reset/restart boundary.
    do n = 1, sd_num
       phase_state = sdm_cold_phase_state(sd_liqice(n), sd_r(n), sdi%re(n), sdi%rp(n))
       select case(phase_state)
       case(PHASE_LIQUID)
          liq_radius_max(n) = max(liq_radius_max(n), max(sd_r(n), 0.0_RP))
       case(PHASE_ICE)
          ice_rvol_max(n) = max(ice_rvol_max(n), sdm_cold_ice_rvol(sdi%re(n), sdi%rp(n)))
       case(PHASE_MIXED)
          mixed_rvol_max(n) = max(mixed_rvol_max(n), &
               sdm_cold_mixed_rvol(sd_r(n), sdi%re(n), sdi%rp(n)))
       end select

       if( phase_state == PHASE_ICE .or. phase_state == PHASE_MIXED ) then
          rime_mass_max(n) = max(rime_mass_max(n), max(sdi%mrime(n), 0.0_RP))
          rime_frac_max(n) = max(rime_frac_max(n), &
               sdm_cold_rime_fraction(sdi%mrime(n), sdi%re(n), sdi%rp(n), sdi%rho(n)))
          nmono_max(n) = max(nmono_max(n), real(max(sdi%nmono(n), 0), kind=RP))
          aspect_ratio_max(n) = max(aspect_ratio_max(n), &
               sdm_cold_aspect_ratio(sdi%re(n), sdi%rp(n)))
       end if
    end do
  end subroutine sdm_cold_tracking_update_interval

  subroutine sdm_cold_tracking_reset_interval( &
       sd_num, sd_liqice, sd_r, sdi, event_mask, event_sig_mask, diag_mask, phase_change_flag, &
       spatial_visit_flag, &
       liq_radius_max, ice_rvol_max, mixed_rvol_max, &
       rime_mass_max, rime_frac_max, nmono_max, aspect_ratio_max)
    integer, intent(in) :: sd_num
    integer(i2), intent(in) :: sd_liqice(1:sd_num)
    real(RP), intent(in) :: sd_r(1:sd_num)
    type(sdicedef), intent(in) :: sdi
    integer, intent(inout) :: event_mask(1:sd_num)
    integer, intent(inout) :: event_sig_mask(1:sd_num)
    integer, intent(inout) :: diag_mask(1:sd_num)
    integer, intent(inout) :: phase_change_flag(1:sd_num)
    integer, intent(inout) :: spatial_visit_flag(1:sd_num)
    real(RP), intent(inout) :: liq_radius_max(1:sd_num)
    real(RP), intent(inout) :: ice_rvol_max(1:sd_num)
    real(RP), intent(inout) :: mixed_rvol_max(1:sd_num)
    real(RP), intent(inout) :: rime_mass_max(1:sd_num)
    real(RP), intent(inout) :: rime_frac_max(1:sd_num)
    real(RP), intent(inout) :: nmono_max(1:sd_num)
    real(RP), intent(inout) :: aspect_ratio_max(1:sd_num)

    ! Interval masks are occurrence proxies, not event counters.
    event_mask(:) = 0
    event_sig_mask(:) = 0
    diag_mask(:) = 0
    phase_change_flag(:) = 0
    spatial_visit_flag(:) = 0
    liq_radius_max(:) = 0.0_RP
    ice_rvol_max(:) = 0.0_RP
    mixed_rvol_max(:) = 0.0_RP
    rime_mass_max(:) = 0.0_RP
    rime_frac_max(:) = 0.0_RP
    nmono_max(:) = 0.0_RP
    aspect_ratio_max(:) = 0.0_RP
    call sdm_cold_tracking_update_interval( &
         sd_num, sd_liqice, sd_r, sdi, &
         liq_radius_max, ice_rvol_max, mixed_rvol_max, &
         rime_mass_max, rime_frac_max, nmono_max, aspect_ratio_max)
  end subroutine sdm_cold_tracking_reset_interval

  subroutine sdm_cold_tracking_update_spatial_visit(sd_num, sd_x, sd_y, sd_z, sd_rk, spatial_visit_flag)
    integer, intent(in) :: sd_num
    real(RP), intent(in) :: sd_x(1:sd_num)
    real(RP), intent(in) :: sd_y(1:sd_num)
    real(RP), intent(in) :: sd_z(1:sd_num)
    real(RP), intent(in) :: sd_rk(1:sd_num)
    integer, intent(inout) :: spatial_visit_flag(1:sd_num)
    integer :: n

    if( .not. tracking_spatial_visit_enable ) return
    if( .not. tracking_spatial_active_rank ) return

    do n = 1, sd_num
       if( sd_rk(n) < VALID2INVALID ) cycle
       if( sd_x(n) >= tracking_spatial_x_min .and. sd_x(n) <= tracking_spatial_x_max .and. &
           sd_y(n) >= tracking_spatial_y_min .and. sd_y(n) <= tracking_spatial_y_max .and. &
           sd_z(n) >= tracking_spatial_z_min .and. sd_z(n) <= tracking_spatial_z_max ) then
          spatial_visit_flag(n) = 1
       end if
    end do
  end subroutine sdm_cold_tracking_update_spatial_visit

  subroutine sdm_cold_tracking_update_spatial_visit_segment( &
       sd_num, sd_x_prev, sd_y_prev, sd_z_prev, sd_x, sd_y, sd_z, sd_rk, spatial_visit_flag)
    integer, intent(in) :: sd_num
    real(RP), intent(in) :: sd_x_prev(1:sd_num)
    real(RP), intent(in) :: sd_y_prev(1:sd_num)
    real(RP), intent(in) :: sd_z_prev(1:sd_num)
    real(RP), intent(in) :: sd_x(1:sd_num)
    real(RP), intent(in) :: sd_y(1:sd_num)
    real(RP), intent(in) :: sd_z(1:sd_num)
    real(RP), intent(in) :: sd_rk(1:sd_num)
    integer, intent(inout) :: spatial_visit_flag(1:sd_num)
    integer :: n

    if( .not. tracking_spatial_visit_enable ) return
    if( .not. tracking_spatial_active_rank ) return

    do n = 1, sd_num
       if( sd_rk(n) < VALID2INVALID ) cycle
       if( sdm_cold_segment_intersects_box( &
            sd_x_prev(n), sd_y_prev(n), sd_z_prev(n), sd_x(n), sd_y(n), sd_z(n), &
            tracking_spatial_x_min, tracking_spatial_x_max, &
            tracking_spatial_y_min, tracking_spatial_y_max, &
            tracking_spatial_z_min, tracking_spatial_z_max) ) then
          spatial_visit_flag(n) = 1
       end if
    end do
  end subroutine sdm_cold_tracking_update_spatial_visit_segment

  logical function sdm_cold_segment_intersects_box( &
       x0, y0, z0, x1, y1, z1, xmin, xmax, ymin, ymax, zmin, zmax) result(intersects)
    real(RP), intent(in) :: x0, y0, z0, x1, y1, z1
    real(RP), intent(in) :: xmin, xmax, ymin, ymax, zmin, zmax
    real(RP) :: tmin, tmax
    logical :: ok

    tmin = 0.0_RP
    tmax = 1.0_RP
    call update_axis(x0, x1, xmin, xmax, tmin, tmax, ok)
    if( .not. ok ) then
       intersects = .false.
       return
    end if
    call update_axis(y0, y1, ymin, ymax, tmin, tmax, ok)
    if( .not. ok ) then
       intersects = .false.
       return
    end if
    call update_axis(z0, z1, zmin, zmax, tmin, tmax, ok)
    intersects = ok .and. tmax >= tmin

  contains
    subroutine update_axis(p0, p1, bmin, bmax, tmin_io, tmax_io, ok_axis)
      real(RP), intent(in) :: p0, p1, bmin, bmax
      real(RP), intent(inout) :: tmin_io, tmax_io
      logical, intent(out) :: ok_axis
      real(RP) :: dp, t1, t2, tlo, thi

      dp = p1 - p0
      if( abs(dp) <= tiny(1.0_RP) ) then
         ok_axis = p0 >= bmin .and. p0 <= bmax
         return
      end if

      t1 = (bmin - p0) / dp
      t2 = (bmax - p0) / dp
      tlo = min(t1, t2)
      thi = max(t1, t2)
      tmin_io = max(tmin_io, tlo)
      tmax_io = min(tmax_io, thi)
      ok_axis = tmax_io >= tmin_io
    end subroutine update_axis
  end function sdm_cold_segment_intersects_box

  subroutine sdm_cold_tracking_evaluate_diag_masks( &
       sd_num, diag_mask, &
       liq_radius_max, ice_rvol_max, mixed_rvol_max, &
       rime_mass_max, rime_frac_max, nmono_max, aspect_ratio_max)
    integer, intent(in) :: sd_num
    integer, intent(inout) :: diag_mask(1:sd_num)
    real(RP), intent(in) :: liq_radius_max(1:sd_num)
    real(RP), intent(in) :: ice_rvol_max(1:sd_num)
    real(RP), intent(in) :: mixed_rvol_max(1:sd_num)
    real(RP), intent(in) :: rime_mass_max(1:sd_num)
    real(RP), intent(in) :: rime_frac_max(1:sd_num)
    real(RP), intent(in) :: nmono_max(1:sd_num)
    real(RP), intent(in) :: aspect_ratio_max(1:sd_num)
    integer :: n

    do n = 1, sd_num
       ! sd_diag_mask is evaluated from interval maxima at SD output time.
       diag_mask(n) = sdm_cold_diag_mask_from_values( &
            liq_radius_max(n), ice_rvol_max(n), mixed_rvol_max(n), &
            rime_mass_max(n), rime_frac_max(n), nmono_max(n), aspect_ratio_max(n), &
            tracking_diag_liq_radius_enable, tracking_diag_liq_radius_threshold, &
            tracking_diag_ice_rvol_enable, tracking_diag_ice_rvol_threshold, &
            tracking_diag_mixed_rvol_enable, tracking_diag_mixed_rvol_threshold, &
            tracking_diag_rime_mass_enable, tracking_diag_rime_mass_threshold, &
            tracking_diag_rime_frac_enable, tracking_diag_rime_frac_threshold, &
            tracking_diag_nmono_enable, tracking_diag_nmono_threshold, &
            tracking_diag_aspect_ratio_enable, tracking_diag_aspect_ratio_threshold)
    end do
  end subroutine sdm_cold_tracking_evaluate_diag_masks
end module m_sdm_tracking_cold
