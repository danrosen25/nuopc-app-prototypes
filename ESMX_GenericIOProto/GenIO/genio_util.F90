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
  use genio_mod_struct
  use genio_mod_params

  implicit none

  private

  public genio_hconfig2i4
  public genio_hconfig2r8
  public genio_hconfig2str
  public genio_hconfig2logical
  public genio_hconfig2timeint
  public genio_hconfig2csys
  public genio_hconfig2ftyp

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
        msg="GENIO: logical key not found - "//trim(key), &
        line=__LINE__, file=__FILE__, rcToReturn=rc)
      return
    endif
  endfunction genio_hconfig2logical

  !-----------------------------------------------------------------------------

  function genio_hconfig2timeint(hconfig, key, defaultValue, rc)
    ! return value
    type(ESMF_TimeInterval) :: genio_hconfig2timeint
    ! arguments
    type(ESMF_HConfig), intent(in)                 :: hconfig
    character(*), intent(in)                       :: key
    type(ESMF_TimeInterval), intent(in), optional  :: defaultValue
    integer, intent(out)                           :: rc
    ! local variables
    character(len=:), allocatable :: strValue
    logical                       :: isPresent
    logical                       :: check
    integer(ESMF_KIND_I4)         :: yy
    integer(ESMF_KIND_I4)         :: mm
    integer(ESMF_KIND_I4)         :: d
    integer(ESMF_KIND_I4)         :: h
    integer(ESMF_KIND_I4)         :: m
    integer(ESMF_KIND_I4)         :: s
    character                     :: x1
    character                     :: x2
    character                     :: x3
    character                     :: x4
    character                     :: x5

    rc = ESMF_SUCCESS

    call ESMF_TimeIntervalSet(genio_hconfig2timeint, s=0, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    isPresent = ESMF_HConfigIsDefined(hconfig, keyString=key, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    if (isPresent) then
      strValue = ESMF_HConfigAsString(hconfig, keyString=key, asOkay=check, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      if (.not.check) then
        call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
          msg="GENIO: TimeInterval cannot be converted to String", &
          line=__LINE__, file=__FILE__, rcToReturn=rc)
        return
      endif
      read(strValue, '(I4,5(A1,I2))', err=8) yy, x1, mm, x2, d, x3, h, x4, m, x5, s
      deallocate(strValue)
      call ESMF_TimeIntervalSet(genio_hconfig2timeint, &
        yy=yy, mm=mm, d=d, h=h, m=m, s=s, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    elseif (present(defaultValue)) then
      genio_hconfig2timeint = defaultValue
    else
      call ESMF_LogSetError(ESMF_RC_NOT_FOUND, &
        msg="GENIO: timeint key not found - "//trim(key), &
        line=__LINE__, file=__FILE__, rcToReturn=rc)
      return
    endif

4   return
8   call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
      msg="GENIO: TimeInterval format ISO 8601 - "//strValue, &
      line=__LINE__, file=__FILE__, rcToReturn=rc)

  endfunction genio_hconfig2timeint

  !-----------------------------------------------------------------------------

  function genio_hconfig2csys(hconfig, key, defaultValue, rc)
    ! return value
    type(ESMF_CoordSys_Flag) :: genio_hconfig2csys
    ! arguments
    type(ESMF_HConfig), intent(in)                 :: hconfig
    character(*), intent(in)                       :: key
    type(ESMF_CoordSys_Flag), intent(in), optional :: defaultValue
    integer, intent(out)                           :: rc
    ! local variables
    character(len=:), allocatable :: strValue
    logical :: isPresent
    logical :: check

    rc = ESMF_SUCCESS

    genio_hconfig2csys = ESMF_COORDSYS_SPH_DEG

    isPresent = ESMF_HConfigIsDefined(hconfig, keyString=key, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    if (isPresent) then
      strValue = ESMF_HConfigAsString(hconfig, keyString=key, asOkay=check, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      if (.not.check) then
        call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
          msg="GENIO: Field coordSys cannot be converted to String", &
          line=__LINE__, file=__FILE__, rcToReturn=rc)
        return
      endif
      strValue = ESMF_UtilStringUpperCase(strValue, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      select case (strValue)
        case ("ESMF_COORDSYS_CART")
          genio_hconfig2csys = ESMF_COORDSYS_CART
        case ("ESMF_COORDSYS_SPH_DEG")
          genio_hconfig2csys = ESMF_COORDSYS_SPH_DEG
        case ("ESMF_COORDSYS_SPH_RAD")
          genio_hconfig2csys = ESMF_COORDSYS_SPH_RAD
        case default
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
            msg="GENIO: Field coordSys is invalid - "//trim(strValue), &
            line=__LINE__, file=__FILE__, rcToReturn=rc)
        return
      endselect
      deallocate(strValue)
    elseif (present(defaultValue)) then
      genio_hconfig2csys = defaultValue
    else
      call ESMF_LogSetError(ESMF_RC_NOT_FOUND, &
        msg="GENIO: csys key not found - "//trim(key), &
        line=__LINE__, file=__FILE__, rcToReturn=rc)
      return
    endif

  endfunction genio_hconfig2csys

  !-----------------------------------------------------------------------------

  function genio_hconfig2ftyp(hconfig, key, defaultValue, rc)
    ! return value
    integer :: genio_hconfig2ftyp
    ! arguments
    type(ESMF_HConfig), intent(in) :: hconfig
    character(*), intent(in)       :: key
    integer, intent(in), optional  :: defaultValue
    integer, intent(out)           :: rc
    ! local variables
    character(len=:), allocatable :: strValue
    logical                       :: isPresent
    logical                       :: check

    rc = ESMF_SUCCESS

    genio_hconfig2ftyp = GENIO_ERROR

    isPresent = ESMF_HConfigIsDefined(hconfig, keyString=key, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    if (isPresent) then
      strValue = ESMF_HConfigAsString(hconfig, keyString=key, asOkay=check, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      if (.not.check) then
        call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
          msg="GENIO: Field ftyp cannot be converted to String", &
          line=__LINE__, file=__FILE__, rcToReturn=rc)
        return
      endif
      strValue = ESMF_UtilStringUpperCase(strValue, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      select case (strValue)
        case ("MEAN","GENIO_MEAN")
          genio_hconfig2ftyp = GENIO_MEAN
        case ("INST","GENIO_INST","INSTANTANEOUS")
          genio_hconfig2ftyp = GENIO_INST
        case ("ACCM","GENIO_ACCM","ACCUMULATE")
          genio_hconfig2ftyp = GENIO_ACCM
        case default
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
            msg="GENIO: Field output type is invalid - "//trim(strValue), &
            line=__LINE__, file=__FILE__, rcToReturn=rc)
        return
      endselect
      deallocate(strValue)
    elseif (present(defaultValue)) then
      genio_hconfig2ftyp = defaultValue
    else
      call ESMF_LogSetError(ESMF_RC_NOT_FOUND, &
        msg="GENIO: ftyp key not found - "//trim(key), &
        line=__LINE__, file=__FILE__, rcToReturn=rc)
      return
    endif

  endfunction genio_hconfig2ftyp

  !-----------------------------------------------------------------------------

endmodule genio_mod_util
