!==============================================================================
! Earth System Modeling Framework
! Copyright (c) 2002-2024, University Corporation for Atmospheric Research,
! Massachusetts Institute of Technology, Geophysical Fluid Dynamics
! Laboratory, University of Michigan, National Centers for Environmental
! Prediction, Los Alamos National Laboratory, Argonne National Laboratory,
! NASA Goddard Space Flight Center.
! Licensed under the University of Illinois-NCSA License.
!==============================================================================

module genio_mod_util

  !-----------------------------------------------------------------------------
  ! Generic IO Utilities
  !-----------------------------------------------------------------------------

  use ESMF
  use genio_mod_param

  implicit none

  private

  public genio_hconfig2i4
  public genio_hconfig2r8
  public genio_hconfig2str
  public genio_hconfig2logical
  public genio_hconfig2priority

  contains

  !-----------------------------------------------------------------------------
  ! Generic IO HConfig Subroutines
  !-----------------------------------------------------------------------------

  function genio_hconfig2i4(hconfig, key, defaultValue, rc)
    ! return value
    integer(ESMF_KIND_I4) :: genio_hconfig2i4
    ! arguments
    type(ESMF_HConfig), intent(in)              :: hconfig
    character(*), intent(in)                    :: key
    integer(ESMF_KIND_I4), intent(in), optional :: defaultValue
    integer, intent(out)                        :: rc
    ! local variables
    logical :: isPresent
    logical :: check

    rc = ESMF_SUCCESS
    genio_hconfig2i4 = 0

    isPresent = ESMF_HConfigIsDefined(hconfig, keyString=key, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    if (isPresent) then
      genio_hconfig2i4 = ESMF_HConfigAsI4(hconfig, keyString=key, &
        asOkay=check, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      if (.not.check) then
        call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
          msg="GENIO: Value cannot be converted to I4 - "//trim(key), &
          line=__LINE__, file=__FILE__, rcToReturn=rc)
        return
      endif
    elseif (present(defaultValue)) then
      genio_hconfig2i4 = defaultValue
    else
      call ESMF_LogSetError(ESMF_RC_NOT_FOUND, &
        msg="GENIO: Key not found - "//trim(key), &
        line=__LINE__, file=__FILE__, rcToReturn=rc)
      return
    endif
  endfunction genio_hconfig2i4

  !-----------------------------------------------------------------------------

  function genio_hconfig2r8(hconfig, key, defaultValue, rc)
    ! return value
    real(ESMF_KIND_R8) :: genio_hconfig2r8
    ! arguments
    type(ESMF_HConfig), intent(in)           :: hconfig
    character(*), intent(in)                 :: key
    real(ESMF_KIND_R8), intent(in), optional :: defaultValue
    integer, intent(out)                     :: rc
    ! local variables
    logical :: isPresent
    logical :: check

    rc = ESMF_SUCCESS
    genio_hconfig2r8 = 0.0_ESMF_KIND_R8

    isPresent = ESMF_HConfigIsDefined(hconfig, keyString=key, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    if (isPresent) then
      genio_hconfig2r8 = ESMF_HConfigAsR8(hconfig, keyString=key, &
        asOkay=check, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      if (.not.check) then
        call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
          msg="GENIO: Value cannot be converted to R8 - "//trim(key), &
          line=__LINE__, file=__FILE__, rcToReturn=rc)
        return
      endif
    elseif (present(defaultValue)) then
      genio_hconfig2r8 = defaultValue
    else
      call ESMF_LogSetError(ESMF_RC_NOT_FOUND, &
        msg="GENIO: Key not found - "//trim(key), &
        line=__LINE__, file=__FILE__, rcToReturn=rc)
      return
    endif
  endfunction genio_hconfig2r8

  !-----------------------------------------------------------------------------

  function genio_hconfig2str(hconfig, key, defaultValue, rc)
    ! return value
    character(:), allocatable :: genio_hconfig2str
    ! arguments
    type(ESMF_HConfig), intent(in)     :: hconfig
    character(*), intent(in)           :: key
    character(*), intent(in), optional :: defaultValue
    integer, intent(out)               :: rc
    ! local variables
    logical :: isPresent
    logical :: check

    rc = ESMF_SUCCESS
    genio_hconfig2str = ' '

    isPresent = ESMF_HConfigIsDefined(hconfig, keyString=key, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    if (isPresent) then
      genio_hconfig2str = ESMF_HConfigAsString(hconfig, keyString=key, &
        asOkay=check, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      if (.not.check) then
        call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
          msg="GENIO: Value cannot be converted to String - "//trim(key), &
          line=__LINE__, file=__FILE__, rcToReturn=rc)
        return
      endif
    elseif (present(defaultValue)) then
      genio_hconfig2str = defaultValue
    else
      call ESMF_LogSetError(ESMF_RC_NOT_FOUND, &
        msg="GENIO: Key not found - "//trim(key), &
        line=__LINE__, file=__FILE__, rcToReturn=rc)
      return
    endif
  endfunction genio_hconfig2str

  !-----------------------------------------------------------------------------

  function genio_hconfig2logical(hconfig, key, defaultValue, rc)
    ! return value
    logical :: genio_hconfig2logical
    ! arguments
    type(ESMF_HConfig), intent(in) :: hconfig
    character(*), intent(in)       :: key
    logical, intent(in), optional  :: defaultValue
    integer, intent(out)           :: rc
    ! local variables
    logical :: isPresent
    logical :: check

    rc = ESMF_SUCCESS
    genio_hconfig2logical = .false.

    isPresent = ESMF_HConfigIsDefined(hconfig, keyString=key, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    if (isPresent) then
      genio_hconfig2logical = ESMF_HConfigAsLogical(hconfig, keyString=key, &
        asOkay=check, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      if (.not.check) then
        call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
          msg="GENIO: Value cannot be converted to Logical - "//trim(key), &
          line=__LINE__, file=__FILE__, rcToReturn=rc)
        return
      endif
    elseif (present(defaultValue)) then
      genio_hconfig2logical = defaultValue
    else
      call ESMF_LogSetError(ESMF_RC_NOT_FOUND, &
        msg="GENIO: Key not found - "//trim(key), &
        line=__LINE__, file=__FILE__, rcToReturn=rc)
      return
    endif
  endfunction genio_hconfig2logical

  !-----------------------------------------------------------------------------

  function genio_hconfig2priority(hconfig, defaultValue, rc)
    ! return value
    integer :: genio_hconfig2priority
    ! arguments
    type(ESMF_HConfig), intent(in) :: hconfig
    integer, intent(in), optional  :: defaultValue
    integer, intent(out)           :: rc
    ! local variables
    character(len=:), allocatable :: strValue
    logical :: isPresent
    logical :: check

    rc = ESMF_SUCCESS

    genio_hconfig2priority = p_optional

    isPresent = ESMF_HConfigIsDefined(hconfig, keyString="priority", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    if (isPresent) then
      strValue = ESMF_HConfigAsString(hconfig, keyString="priority", &
        asOkay=check, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      if (.not.check) then
        call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
          msg="GENIO: Field priority cannot be converted to String", &
          line=__LINE__, file=__FILE__, rcToReturn=rc)
        return
      endif
      strValue = ESMF_UtilStringLowerCase(strValue, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      select case (strValue)
        case ("disabled","off")
          genio_hconfig2priority = p_disabled
        case ("required","on")
          genio_hconfig2priority = p_required
        case ("optional")
          genio_hconfig2priority = p_optional
        case default 
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
            msg="GENIO: Field priority is invalid - "//trim(strValue), &
            line=__LINE__, file=__FILE__, rcToReturn=rc)
        return
      endselect
      deallocate(strValue)
    elseif (present(defaultValue)) then
      genio_hconfig2priority = p_optional
    else
      call ESMF_LogSetError(ESMF_RC_NOT_FOUND, &
        msg="GENIO: Field priroty not found", &
        line=__LINE__, file=__FILE__, rcToReturn=rc)
      return
    endif

  endfunction genio_hconfig2priority

  !-----------------------------------------------------------------------------

endmodule genio_mod_util
