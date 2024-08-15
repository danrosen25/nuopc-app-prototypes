!==============================================================================
! Earth System Modeling Framework
! Copyright (c) 2002-2024, University Corporation for Atmospheric Research,
! Massachusetts Institute of Technology, Geophysical Fluid Dynamics
! Laboratory, University of Michigan, National Centers for Environmental
! Prediction, Los Alamos National Laboratory, Argonne National Laboratory,
! NASA Goddard Space Flight Center.
! Licensed under the University of Illinois-NCSA License.
!==============================================================================

module genio

  !-----------------------------------------------------------------------------
  ! Generic IO Model
  !-----------------------------------------------------------------------------

  use ESMF
  use NUOPC
  use NUOPC_Model, &
    modelSS => SetServices
  use genio_mod_struct
  use genio_mod_params
  use genio_mod_util
  use genio_mod_defaults
  use genio_mod_geom
  use genio_mod_field

  implicit none

  private

  public SetServices, SetVM

  contains

  !-----------------------------------------------------------------------------
  ! Generic IO Model Specialization
  !-----------------------------------------------------------------------------

  subroutine SetServices(genio, rc)
    ! arguments
    type(ESMF_GridComp)  :: genio
    integer, intent(out) :: rc
    ! local variables
    integer                    :: stat
    type(geniostate_wrap)      :: is
    type(genio_state), pointer :: geniostate
    character(len=64)          :: value
    type(ESMF_Config)          :: config
    type(ESMF_HConfig)         :: hconfig
    character(80)              :: compLabel
    character(:), allocatable  :: badKey
    logical                    :: isFlag
    logical                    :: check
    logical                    :: compcfg_p
    logical                    :: dfltcfg_p
    logical                    :: outpcfg_p
    logical                    :: geomcfg_p
    logical                    :: flstcfg_p
    type(ESMF_HConfig)         :: compcfg

    rc = ESMF_SUCCESS

    ! derive generic model phases
    call NUOPC_CompDerive(genio, modelSS, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! allocate memory for this internal state and set it in the component
    allocate(is%ptr, stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg='GENIO: Memory allocation failed.', &
      line=__LINE__, &
      file=__FILE__, &
      rcToReturn=rc)) return
    call ESMF_GridCompSetInternalState(genio, is, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    geniostate => is%ptr
    if (.not. associated(geniostate)) then
      call ESMF_LogSetError(ESMF_RC_PTR_NOTALLOC, &
        msg='GENIO: geniostate has not been associated', &
        line=__LINE__, file=__FILE__, rcToReturn=rc)
      return
    endif

    ! query component for information
    call NUOPC_CompGet(genio, name=geniostate%cname, &
      verbosity=geniostate%verbosity, diagnostic=geniostate%diagnostic, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! specialize model
    call NUOPC_CompSpecialize(genio, specLabel=label_Advertise, &
      specRoutine=Advertise, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call NUOPC_CompSpecialize(genio, specLabel=label_ModifyAdvertised, &
      specRoutine=ModifyAdvertised, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call NUOPC_CompSpecialize(genio, specLabel=label_RealizeProvided, &
      specRoutine=RealizeProvided, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call NUOPC_CompSpecialize(genio, specLabel=label_RealizeAccepted, &
      specRoutine=RealizeAccepted, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call NUOPC_CompSpecialize(genio, specLabel=label_DataInitialize, &
      specRoutine=DataInitialize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call NUOPC_CompSpecialize(genio, specLabel=label_CheckImport, &
       specRoutine=CheckImport, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call NUOPC_CompSpecialize(genio, specLabel=label_Advance, &
      specRoutine=ModelAdvance, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call NUOPC_CompSpecialize(genio, specLabel=label_Finalize, &
      specRoutine=ModelFinalize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! setup component variables
    geniostate%outid = 0
    call ESMF_GridCompGet(genio, vm=geniostate%vm, &
      localPet=geniostate%myid, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call ESMF_TimeIntervalSet(geniostate%zerotime, s=0, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    call ESMF_GridCompGet(genio, name=compLabel, configIsPresent=isFlag, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return  ! bail out
    if (isFlag) then
      ! Config present, assert it is in the ESMX YAML format
      call ESMF_GridCompGet(genio, config=config, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,    file=__FILE__)) return  ! bail out
      call ESMF_ConfigGet(config, hconfig=hconfig, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return  ! bail out
      compcfg_p = ESMF_HConfigIsDefined(hconfig, keyString=compLabel, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    else
      compcfg_p = .false.
    endif

    ! create component configuration
    if (compcfg_p) then
      compcfg = ESMF_HConfigCreateAt(hconfig, keyString=compLabel, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return  ! bail out
      check = ESMF_HConfigValidateMapKeys(compcfg, &
        vocabulary=["model         ", &  ! ESMX_Driver handled option
                    "petList       ", &  ! ESMX_Driver handled option
                    "ompNumThreads ", &  ! ESMX_Driver handled option
                    "attributes    ", &  ! ESMX_Driver handled option
                    "defaultOptions", &  ! GenIO handled option
                    "outputOptions ", &  ! GenIO handled option
                    "fieldOptions  ", &  ! GenIO handled option
                    "geomOptions   "  &  ! GenIO handled option
                   ], badKey=badKey, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return  ! bail out
      if (.not.check) then
        call ESMF_LogSetError(ESMF_RC_ARG_WRONG, &
          msg="An invalid key was found in config under "//trim(compLabel)// &
            " (maybe a typo?): "//badKey, &
          line=__LINE__, file=__FILE__, rcToReturn=rc)
        return
      endif
      dfltcfg_p = ESMF_HConfigIsDefined(compcfg, &
        keyString="defaultOptions", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      outpcfg_p = ESMF_HConfigIsDefined(compcfg, &
        keyString="outputOptions", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      geomcfg_p = ESMF_HConfigIsDefined(compcfg, &
        keyString="geom", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      flstcfg_p = ESMF_HConfigIsDefined(compcfg, &
        keyString="fieldOptions", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    else
      dfltcfg_p = .false.
      outpcfg_p = .false.
      geomcfg_p = .false.
      flstcfg_p = .false.
    endif

    if (dfltcfg_p) then
      geniostate%dfltcfg = ESMF_HConfigCreateAt(compcfg, &
        keyString="defaultOptions", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    else
      geniostate%dfltcfg = ESMF_HConfigCreate(content="{}", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    endif
    if (outpcfg_p) then
      geniostate%outpcfg = ESMF_HConfigCreateAt(compcfg, &
        keyString="outputOptions", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    else
      geniostate%outpcfg = ESMF_HConfigCreate(content="{}", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    endif
    if (geomcfg_p) then
      geniostate%geomcfg = ESMF_HConfigCreateAt(compcfg, &
        keyString="geom", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    else
      geniostate%geomcfg = ESMF_HConfigCreate(content="{}", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    endif
    if (flstcfg_p) then
      geniostate%flstcfg = ESMF_HConfigCreateAt(compcfg, &
        keyString="fieldOptions", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    else
      geniostate%flstcfg = ESMF_HConfigCreate(content="{}", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    endif

    ! destory component configuration
    if (compcfg_p) then
      call ESMF_HConfigDestroy(compcfg, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    endif

  endsubroutine SetServices

  !-----------------------------------------------------------------------------

  subroutine Advertise(genio, rc)
    ! arguments
    type(ESMF_GridComp)  :: genio
    integer, intent(out) :: rc
    ! local variables
    integer                    :: stat
    type(geniostate_wrap)      :: is
    type(genio_state), pointer :: geniostate
    type(ESMF_State)           :: importState

    rc = ESMF_SUCCESS

    ! query component for internal state
    nullify(is%ptr)
    call ESMF_GridCompGetInternalState(genio, is, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    geniostate => is%ptr
    if (.not. associated(geniostate)) then
      call ESMF_LogSetError(ESMF_RC_PTR_NOTALLOC, &
        msg='GENIO: geniostate has not been associated', &
        line=__LINE__, file=__FILE__, rcToReturn=rc)
      return
    endif

    ! query component for information
    call NUOPC_CompGet(genio, name=geniostate%cname, &
      verbosity=geniostate%verbosity, diagnostic=geniostate%diagnostic, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! query for importState
    call NUOPC_ModelGet(genio, importState=importState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! accept all import fields
    call NUOPC_SetAttribute(importState, "FieldTransferPolicy", "transferAll", &
      rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

  endsubroutine Advertise

  !-----------------------------------------------------------------------------

  subroutine ModifyAdvertised(genio, rc)
    type(ESMF_GridComp)  :: genio
    integer, intent(out) :: rc

    ! local variables
    integer                    :: stat
    type(geniostate_wrap)      :: is
    type(genio_state), pointer :: geniostate
    type(ESMF_State)           :: importState
    logical                    :: dfltcfg_p
    logical                    :: outpcfg_p
    logical                    :: geomcfg_p
    logical                    :: flstcfg_p

    rc = ESMF_SUCCESS

    ! query component for internal state
    nullify(is%ptr)
    call ESMF_GridCompGetInternalState(genio, is, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    geniostate => is%ptr
    if (.not. associated(geniostate)) then
      call ESMF_LogSetError(ESMF_RC_PTR_NOTALLOC, &
        msg='GENIO: geniostate has not been associated', &
        line=__LINE__, file=__FILE__, rcToReturn=rc)
      return
    endif

    ! query component for information
    call NUOPC_CompGet(genio, name=geniostate%cname, &
      verbosity=geniostate%verbosity, diagnostic=geniostate%diagnostic, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! query for importState
    call NUOPC_ModelGet(genio, importState=importState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! create defaults
    geniostate%dflts = genio_dflts_create(geniostate%dfltcfg, &
      rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    ! create source data state
    geniostate%srclst = genio_fldlst_create(importstate, geniostate%outpcfg, &
      geniostate%dflts, name=trim(geniostate%cname)//"SrcList", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    ! create geom object
    geniostate%geom = genio_geom_create(geniostate%geomcfg, &
      geniostate%dflts, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

  end subroutine ModifyAdvertised

  !-----------------------------------------------------------------------------

  subroutine RealizeProvided(genio, rc)
    ! arguments
    type(ESMF_GridComp)  :: genio
    integer, intent(out) :: rc
    ! local variables
    integer                    :: stat
    type(geniostate_wrap)      :: is
    type(genio_state), pointer :: geniostate
    type(ESMF_State)           :: importState
    integer                    :: i

    rc = ESMF_SUCCESS

    ! query component for internal state
    nullify(is%ptr)
    call ESMF_GridCompGetInternalState(genio, is, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    geniostate => is%ptr
    if (.not. associated(geniostate)) then
      call ESMF_LogSetError(ESMF_RC_PTR_NOTALLOC, &
        msg='GENIO: geniostate has not been associated', &
        line=__LINE__, file=__FILE__, rcToReturn=rc)
      return
    endif

    ! query component for information
    call NUOPC_CompGet(genio, name=geniostate%cname, &
      verbosity=geniostate%verbosity, diagnostic=geniostate%diagnostic, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! query for importState
    call NUOPC_ModelGet(genio, importState=importState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! write grid to NetCDF file
    if (btest(geniostate%diagnostic,16)) then
      call genio_geom_write(geniostate%geom, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    endif

  endsubroutine RealizeProvided

  !-----------------------------------------------------------------------------

  subroutine RealizeAccepted(genio, rc)
    ! arguments
    type(ESMF_GridComp)  :: genio
    integer, intent(out) :: rc
    ! local variables
    integer                    :: stat
    type(geniostate_wrap)      :: is
    type(genio_state), pointer :: geniostate
    type(ESMF_State)           :: importState
    type(ESMF_TimeInterval)    :: dfltTimeInt

    rc = ESMF_SUCCESS

    ! query component for internal state
    nullify(is%ptr)
    call ESMF_GridCompGetInternalState(genio, is, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    geniostate => is%ptr
    if (.not. associated(geniostate)) then
      call ESMF_LogSetError(ESMF_RC_PTR_NOTALLOC, &
        msg='GENIO: geniostate has not been associated', &
        line=__LINE__, file=__FILE__, rcToReturn=rc)
      return
    endif

    ! query component for information
    call NUOPC_CompGet(genio, name=geniostate%cname, &
      verbosity=geniostate%verbosity, diagnostic=geniostate%diagnostic, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! query for importState
    call NUOPC_ModelGet(genio, importState=importState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! realize import fields
    call genio_fldlst_realize(importState, geniostate%srclst, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! default time interval
    call ESMF_TimeIntervalSet(dfltTimeInt, s=0, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! create output state
    geniostate%outlst = genio_fldlst_create(geniostate%srclst, &
      geniostate%outpcfg, geniostate%flstcfg, geniostate%dflts, &
      name=trim(geniostate%cname)//"OutLst", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! reset source data
    call genio_fldlst_reset(geniostate%srclst, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

  endsubroutine RealizeAccepted

  !-----------------------------------------------------------------------------

  subroutine DataInitialize(genio, rc)
    ! arguments
    type(ESMF_GridComp)  :: genio
    integer, intent(out) :: rc
    ! local variables
    integer                    :: stat
    type(geniostate_wrap)      :: is
    type(genio_state), pointer :: geniostate
    type(ESMF_Clock)           :: modelClock
    logical                    :: allDataSatisfied

    rc = ESMF_SUCCESS

    ! query component for internal state
    nullify(is%ptr)
    call ESMF_GridCompGetInternalState(genio, is, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    geniostate => is%ptr
    if (.not. associated(geniostate)) then
      call ESMF_LogSetError(ESMF_RC_PTR_NOTALLOC, &
        msg='GENIO: geniostate has not been associated', &
        line=__LINE__, file=__FILE__, rcToReturn=rc)
      return
    endif

    ! query component for information
    call NUOPC_CompGet(genio, name=geniostate%cname, &
      verbosity=geniostate%verbosity, diagnostic=geniostate%diagnostic, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! query component for clock
    call NUOPC_ModelGet(genio, modelClock=modelClock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    call genio_fldlst_reset(geniostate%outlst, modelClock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    allDataSatisfied = genio_fldlst_alldata(geniostate%outlst, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    if (allDataSatisfied) then
      call NUOPC_CompAttributeSet(genio, &
        name="InitializeDataComplete", value="true", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    else
      call NUOPC_CompAttributeSet(genio, &
        name="InitializeDataComplete", value="false", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    endif

  endsubroutine DataInitialize

  !-----------------------------------------------------------------------------

  subroutine CheckImport(genio, rc)
    ! arguments
    type(ESMF_GridComp) :: genio
    integer,intent(out) :: rc
    ! local variables
    integer                    :: stat
    type(geniostate_wrap)      :: is
    type(genio_state), pointer :: geniostate
    type(ESMF_Clock)           :: modelClock
    type(ESMF_Time)            :: currTime
    type(ESMF_Time)            :: modelCurrTime
    type(ESMF_State)           :: importState
    logical                    :: allCurrTime

    rc = ESMF_SUCCESS

    ! query component for internal State
    nullify(is%ptr)
    call ESMF_GridCompGetInternalState(genio, is, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    geniostate => is%ptr
    if (.not. associated(geniostate)) then
      call ESMF_LogSetError(ESMF_RC_PTR_NOTALLOC, &
        msg='GENIO: geniostate has not been associated', &
        line=__LINE__, file=__FILE__, rcToReturn=rc)
      return
    endif

    ! query component for information
    call NUOPC_CompGet(genio, name=geniostate%cname, &
      verbosity=geniostate%verbosity, diagnostic=geniostate%diagnostic, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! query the component for its clock and import state
    call NUOPC_ModelGet(genio, modelClock=modelClock, &
      importState=importState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! get the stop time out of the clock
    call ESMF_ClockGet(modelClock, currTime=modelCurrTime, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    allCurrTime = NUOPC_IsAtTime(importState, modelCurrTime, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    if (.NOT.allCurrTime) then
      call ESMF_LogWrite(trim(geniostate%cname)//": "// &
        "NUOPC INCOMPATIBILITY DETECTED: Import Fields not at current time", &
        ESMF_LOGMSG_WARNING)
    endif
  endsubroutine CheckImport

  !-----------------------------------------------------------------------------

  subroutine ModelAdvance(genio, rc)
    ! arguments
    type(ESMF_GridComp)  :: genio
    integer, intent(out) :: rc
    ! local variables
    integer                    :: stat
    type(geniostate_wrap)      :: is
    type(genio_state), pointer :: geniostate
    type(ESMF_Clock)           :: modelClock
    type(ESMF_Time)            :: currTime
    type(ESMF_TimeInterval)    :: accStep
    logical                    :: write_out

    rc = ESMF_SUCCESS

    ! query component for internal state
    nullify(is%ptr)
    call ESMF_GridCompGetInternalState(genio, is, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    geniostate => is%ptr
    if (.not. associated(geniostate)) then
      call ESMF_LogSetError(ESMF_RC_PTR_NOTALLOC, &
        msg='GENIO: geniostate has not been associated', &
        line=__LINE__, file=__FILE__, rcToReturn=rc)
      return
    endif

    ! query component for information
    call NUOPC_CompGet(genio, name=geniostate%cname, &
      verbosity=geniostate%verbosity, diagnostic=geniostate%diagnostic, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! query component for modelClock
    call NUOPC_ModelGet(genio, modelClock=modelClock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call ESMF_ClockGet(modelClock, currTime=currTime, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! write source state
    if (geniostate%srclst%write) then
      call genio_fldlst_write(geniostate%srclst, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    endif
    if (geniostate%srclst%diagn) then
      call genio_fldlst_print(geniostate%srclst, myid=geniostate%myid, &
        outid=geniostate%outid, vm=geniostate%vm, cname=geniostate%cname, &
        label="ModelAdvance", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    endif

    ! source fields -> output state
    call genio_fldlst_addsrc(geniostate%outlst, modelClock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    if (geniostate%outlst%diagn) then
      call genio_fldlst_print(geniostate%outlst, myid=geniostate%myid, &
        outid=geniostate%outid, vm=geniostate%vm, cname=geniostate%cname, &
        label="ModelAdvance", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    endif

    ! write output state
    write_out = genio_fldlst_alarm(geniostate%outlst, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    if (write_out) then
      call genio_fldlst_write(geniostate%outlst, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call genio_fldlst_reset(geniostate%outlst, modelClock, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    endif

  endsubroutine ModelAdvance

  !-----------------------------------------------------------------------------

  subroutine ModelFinalize(genio, rc)
    ! arguments
    type(ESMF_GridComp)  :: genio
    integer, intent(out) :: rc
    ! local variables
    integer                    :: stat
    type(geniostate_wrap)      :: is
    type(genio_state), pointer :: geniostate
    type(ESMF_Clock)           :: modelClock
    type(ESMF_Time)            :: currTime
    integer                    :: i

    rc = ESMF_SUCCESS

    ! query component for internal state
    nullify(is%ptr)
    call ESMF_GridCompGetInternalState(genio, is, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    geniostate => is%ptr
    if (.not. associated(geniostate)) then
      call ESMF_LogSetError(ESMF_RC_PTR_NOTALLOC, &
        msg='GENIO: geniostate has not been associated', &
        line=__LINE__, file=__FILE__, rcToReturn=rc)
      return
    endif

    ! query component for information
    call NUOPC_CompGet(genio, name=geniostate%cname, &
      verbosity=geniostate%verbosity, diagnostic=geniostate%diagnostic, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! query component for modelClock
    call NUOPC_ModelGet(genio, modelClock=modelClock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call ESMF_ClockGet(modelClock, currTime=currTime, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! write source state
    if (geniostate%srclst%write) then
      call genio_fldlst_write(geniostate%srclst, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    endif
    if (geniostate%srclst%diagn) then
      call genio_fldlst_print(geniostate%srclst, myid=geniostate%myid, &
        outid=geniostate%outid, vm=geniostate%vm, cname=geniostate%cname, &
        label="ModelFinalize", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    endif

    ! source fields -> output state
    call genio_fldlst_addsrc(geniostate%outlst, modelClock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    if (geniostate%outlst%diagn) then
      call genio_fldlst_print(geniostate%outlst, myid=geniostate%myid, &
        outid=geniostate%outid, vm=geniostate%vm, cname=geniostate%cname, &
        label="ModelFinalize", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    endif

    ! write output state
    if (geniostate%outlst%wfinl) then
      call genio_fldlst_write(geniostate%outlst, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    endif

    call genio_fldlst_destroy(geniostate%outlst, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call genio_fldlst_destroy(geniostate%srclst, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call genio_geom_destroy(geniostate%geom, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    call ESMF_HConfigDestroy(geniostate%geomcfg, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call ESMF_HConfigDestroy(geniostate%outpcfg, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call ESMF_HConfigDestroy(geniostate%dfltcfg, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call ESMF_HConfigDestroy(geniostate%flstcfg, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    deallocate(is%ptr, stat=stat)
    if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
      msg='GENIO: Memory deallocation failed.', &
      line=__LINE__, &
      file=__FILE__, &
      rcToReturn=rc)) return
  endsubroutine ModelFinalize

endmodule genio
