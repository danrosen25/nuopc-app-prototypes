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
  use genio_mod_util
  use genio_mod_config
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
    type(ESMF_HConfig)         :: hconfig, hconfigNode
    character(80)              :: compLabel
    character(:), allocatable  :: badKey
    logical                    :: isFlag

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

    ! query component for vm and local pet
    geniostate%outid = 0
    call ESMF_GridCompGet(genio, vm=geniostate%vm, &
      localPet=geniostate%myid, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! validate config
    call ESMF_GridCompGet(genio, name=compLabel, configIsPresent=isFlag, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return  ! bail out
    if (isFlag) then
      geniostate%cfgPresent = .true.
      ! Config present, assert it is in the ESMX YAML format
      call ESMF_GridCompGet(genio, config=config, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,    file=__FILE__)) return  ! bail out
      call ESMF_ConfigGet(config, hconfig=hconfig, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return  ! bail out
      geniostate%cfg = ESMF_HConfigCreateAt(hconfig, keyString=compLabel, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return  ! bail out
      ! component responsibility to validate ESMX handled options here, and
      ! potentially locally handled options
      isFlag = ESMF_HConfigValidateMapKeys(geniostate%cfg, &
        vocabulary=["model        ", &  ! ESMX_Driver handled option
                    "petList      ", &  ! ESMX_Driver handled option
                    "ompNumThreads", &  ! ESMX_Driver handled option
                    "attributes   ", &  ! ESMX_Driver handled option
                    "output       ", &  ! GenIO handled option
                    "geom         ", &  ! GenIO handled option
                    "fieldOptions "  &  ! GenIO handled option
                   ], badKey=badKey, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return  ! bail out
      if (.not.isFlag) then
        call ESMF_LogSetError(ESMF_RC_ARG_WRONG, &
          msg="An invalid key was found in config under "//trim(compLabel)// &
            " (maybe a typo?): "//badKey, &
          line=__LINE__, file=__FILE__, rcToReturn=rc)
        return
      endif
    else
      geniostate%cfgPresent = .false.
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
    logical                    :: outpcfg_p
    logical                    :: geomcfg_p
    logical                    :: flstcfg_p
    type(ESMF_HConfig)         :: outpcfg
    type(ESMF_HConfig)         :: geomcfg
    type(ESMF_HConfig)         :: flstcfg
    character(256)             :: errmsg

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

    ! check configuration information
    if (geniostate%cfgPresent) then
      outpcfg_p = ESMF_HConfigIsDefined(geniostate%cfg, &
        keyString="output", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      geomcfg_p = ESMF_HConfigIsDefined(geniostate%cfg, &
        keyString="geom", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      flstcfg_p = ESMF_HConfigIsDefined(geniostate%cfg, &
        keyString="fieldOptions", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    else
      outpcfg_p = .false.
      geomcfg_p = .false.
      flstcfg_p = .false.
    endif

    ! geometry configuration
    if (geomcfg_p) then
      allocate(geniostate%geolst(1), stat=stat)
      if (ESMF_LogFoundAllocError(statusToCheck=stat, &
        msg='GENIO: Memory allocation failed.', &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)) return
      geomcfg = ESMF_HConfigCreateAt(geniostate%cfg, &
        keyString="geom", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      geniostate%geolst(1) = genio_geom_create(geomcfg, errmsg, rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call ESMF_HConfigDestroy(geomcfg, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    else
      allocate(geniostate%geolst(0), stat=stat)
      if (ESMF_LogFoundAllocError(statusToCheck=stat, &
        msg='GENIO: Memory allocation failed.', &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)) return
    endif
    ! field configuration
    if (flstcfg_p) then
      flstcfg = ESMF_HConfigCreateAt(geniostate%cfg, &
        keyString="fieldOptions", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      geniostate%srclst = genio_fldlst_create(importstate, flstcfg, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call ESMF_HConfigDestroy(flstcfg, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    else
      geniostate%srclst = genio_fldlst_create(importstate, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    endif

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
      do i=1, size(geniostate%geolst)
        call genio_geom_write(geniostate%geolst(i), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
      enddo
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

    ! query component for import and export states
    call NUOPC_ModelGet(genio, importState=importState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    call genio_fldlst_reset(geniostate%srclst, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    call NUOPC_CompAttributeSet(genio, &
      name="InitializeDataComplete", value="true", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

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
    type(ESMF_State)           :: importState
    character(len=160)         :: timeString

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

    ! query component for import and export states
    call NUOPC_ModelGet(genio, modelClock=modelClock, &
      importState=importState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call ESMF_ClockGet(modelClock, currTime=currTime, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call ESMF_TimeGet(currTime, timeString=timeString, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! write field state
    call genio_fldlst_write(geniostate%srclst, &
      fileName=trim(geniostate%cname)//"_"//trim(timeString)//".nc", &
      rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    if (geniostate%verbosity .gt. 0) then
      call genio_fldlst_diagnostics(geniostate%srclst, myid=geniostate%myid, &
        outid=geniostate%outid, vm=geniostate%vm, cname=geniostate%cname, &
        label="ModelAdvance", timeString=timeString, rc=rc)
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
    type(ESMF_State)           :: importState
    character(len=160)         :: timeString
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

    ! query component for import state
    call NUOPC_ModelGet(genio, modelClock=modelClock, &
      importState=importState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call ESMF_ClockGet(modelClock, currTime=currTime, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call ESMF_TimeGet(currTime, timeString=timeString, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    if (geniostate%verbosity .gt. 0) then
      call genio_fldlst_diagnostics(geniostate%srclst, myid=geniostate%myid, &
        outid=geniostate%outid, vm=geniostate%vm, cname=geniostate%cname, &
        label="ModelFinalize", timeString=timeString, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    endif ! diagnostic

    call genio_fldlst_write(geniostate%srclst, &
      fileName=trim(geniostate%cname)//"_ModelFinalize.nc", &
      rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    call genio_fldlst_remove(importState, geniostate%srclst, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    call genio_fldlst_destroy(geniostate%srclst, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    do i=1, size(geniostate%geolst)
      call genio_geom_destroy(geniostate%geolst(i), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    enddo

    call ESMF_HConfigDestroy(geniostate%cfg, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    geniostate%cfgPresent = .false.

    deallocate(is%ptr, stat=stat)
    if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
      msg='GENIO: Memory deallocation failed.', &
      line=__LINE__, &
      file=__FILE__, &
      rcToReturn=rc)) return
  endsubroutine ModelFinalize

endmodule genio
