!==============================================================================
! Earth System Modeling Framework
! Copyright (c) 2002-2024, University Corporation for Atmospheric Research,
! Massachusetts Institute of Technology, Geophysical Fluid Dynamics
! Laboratory, University of Michigan, National Centers for Environmental
! Prediction, Los Alamos National Laboratory, Argonne National Laboratory,
! NASA Goddard Space Flight Center.
! Licensed under the University of Illinois-NCSA License.
!==============================================================================

module genio_mod_config

  !-----------------------------------------------------------------------------
  ! Generic IO Config
  !-----------------------------------------------------------------------------

  use ESMF
  use genio_mod_struct
  use genio_mod_util

  implicit none

  private

  public io_comp_read_output

  contains

  !-----------------------------------------------------------------------------
  ! Generic IO Config Subroutines
  !-----------------------------------------------------------------------------

  subroutine io_comp_read_output(geniostate, rc)
    ! arguments
    type(genio_state), pointer, intent(inout) :: geniostate
    integer, intent(out)                      :: rc
    ! local variables
    logical            :: isPresent
    integer            :: stat
    logical            :: check
    type(ESMF_HConfig) :: outcfg
    character(:), allocatable :: cfgval
    character(:), allocatable :: badKey

    rc = ESMF_SUCCESS

    if (.not. associated(geniostate)) then
      call ESMF_LogSetError(ESMF_RC_PTR_NOTALLOC, &
        msg='GENIO: geniostate has not been associated', &
        line=__LINE__, file=__FILE__, rcToReturn=rc)
      return
    endif

    ! read output configuration
    if (geniostate%cfgPresent) then
      isPresent = ESMF_HConfigIsDefined(geniostate%cfg, &
        keyString="output", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    else
      isPresent = .false.
    endif
    if (isPresent) then
      ! access output
      outcfg = ESMF_HConfigCreateAt(geniostate%cfg, &
        keyString="output", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      check = ESMF_HConfigValidateMapKeys(outcfg, &
        vocabulary=["write_final"], badKey=badKey, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      if (.not. check) then
        call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
          msg=trim(geniostate%cname)//": unknown output option key - "//badKey, &
        line=__LINE__,file=__FILE__, rcToReturn=rc)
        return
      endif
      ! options
      geniostate%write_final = genio_hconfig2logical(outcfg, "write_final", &
        defaultValue=.true., rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call ESMF_HConfigDestroy(outcfg, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    endif ! outcfg

  endsubroutine io_comp_read_output

  !-----------------------------------------------------------------------------

endmodule genio_mod_config
