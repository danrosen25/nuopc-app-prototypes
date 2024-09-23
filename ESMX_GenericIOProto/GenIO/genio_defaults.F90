!==============================================================================
! Earth System Modeling Framework
! Copyright (c) 2002-2024, University Corporation for Atmospheric Research,
! Massachusetts Institute of Technology, Geophysical Fluid Dynamics
! Laboratory, University of Michigan, National Centers for Environmental
! Prediction, Los Alamos National Laboratory, Argonne National Laboratory,
! NASA Goddard Space Flight Center.
! Licensed under the University of Illinois-NCSA License.
!==============================================================================

module genio_mod_defaults

  !-----------------------------------------------------------------------------
  ! Generic IO Defaults
  !-----------------------------------------------------------------------------

  use ESMF
  use NUOPC
  use genio_mod_struct
  use genio_mod_params
  use genio_mod_util

  implicit none

  private

  public genio_dflts_create

  contains

  !-----------------------------------------------------------------------------
  ! Generic IO Defaults Subroutines
  !-----------------------------------------------------------------------------

  function genio_dflts_create(dfltcfg, rc)
    ! return value
    type(genio_dflts) :: genio_dflts_create
    ! arguments
    type(ESMF_HConfig), intent(in) :: dfltcfg
    integer, intent(out)                     :: rc
    ! local variables
    logical                   :: check
    character(:), allocatable :: badKey

    rc = ESMF_SUCCESS

    check = ESMF_HConfigValidateMapKeys(dfltcfg, &
      vocabulary=["fillValue ", &
                  "nx        ", &
                  "ny        ", &
                  "nz        ", &
                  "minx      ", &
                  "maxx      ", &
                  "miny      ", &
                  "maxy      ", &
                  "coordSys  ", &
                  "outputType"  &
                 ], badKey=badKey, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    if (.not. check) then
      call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
        msg="GENIO: unknown option in defaults "//badKey, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)
      return
    endif
    genio_dflts_create%filv = genio_hconfig2r8(dfltcfg, "fillValue", &
      defaultValue=GENIO_FILV, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    genio_dflts_create%nx   = genio_hconfig2i4(dfltcfg, "nx", &
      defaultValue=GENIO_DFLT_NX, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    genio_dflts_create%ny   = genio_hconfig2i4(dfltcfg, "ny", &
      defaultValue=GENIO_DFLT_NY, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    genio_dflts_create%nz   = genio_hconfig2i4(dfltcfg, "nz", &
      defaultValue=GENIO_DFLT_NZ, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    genio_dflts_create%minx = genio_hconfig2r8(dfltcfg, "minx", &
      defaultValue=GENIO_DFLT_MINX, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    genio_dflts_create%maxx = genio_hconfig2r8(dfltcfg, "maxx", &
      defaultValue=GENIO_DFLT_MAXX, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    genio_dflts_create%miny = genio_hconfig2r8(dfltcfg, "miny", &
      defaultValue=GENIO_DFLT_MINY, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    genio_dflts_create%maxy = genio_hconfig2r8(dfltcfg, "maxy", &
      defaultValue=GENIO_DFLT_MAXY, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    genio_dflts_create%csys = genio_hconfig2csys(dfltcfg, "coordSys", &
      defaultValue=GENIO_DFLT_CSYS, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    genio_dflts_create%otyp = genio_hconfig2otyp(dfltcfg, "outputType", &
      defaultValue=GENIO_DFLT_OTYP, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    genio_dflts_create%ffmt = genio_hconfig2ffmt(dfltcfg, "fileformat", &
      defaultValue=GENIO_DFLT_FFMT, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

  endfunction genio_dflts_create

endmodule genio_mod_defaults
