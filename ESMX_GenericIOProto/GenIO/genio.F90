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
    modelSS    => SetServices
  use genio_mod_param
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
                    "fieldsOptions"  &  ! GenIO handled option
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
    type(genio_field), pointer :: iofield => null()

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
    type(ESMF_Field)           :: field

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

    ! get configuration information
    if (geniostate%cfgPresent) then
      call io_comp_read_output(geniostate, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call io_comp_read_geom(geniostate, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call io_comp_add_fields(importstate, geniostate, rc=rc)
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
    type(genio_field), pointer :: iofield

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

    ! create grid
    geniostate%grid = ESMF_GridCreateNoPeriDimUfrm( &
      name=trim(geniostate%cname)//"_grid", &
      minIndex=(/1, 1/), maxIndex=(/geniostate%nx, geniostate%ny/), &
      minCornerCoord=(/geniostate%minx,geniostate%miny/), &
      maxCornerCoord=(/geniostate%maxx,geniostate%maxy/), &
      staggerLocList=(/ESMF_STAGGERLOC_CENTER, ESMF_STAGGERLOC_CORNER/), &
      coordSys=geniostate%coordSys, &
      rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! write grid to NetCDF file
    if (btest(geniostate%diagnostic,16)) then
      call io_comp_grid_diag(geniostate, trim(geniostate%cname)//"_grid.nc", rc=rc)
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
    type(genio_field), pointer :: iofield

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
    iofield => geniostate%imp_flds_head
    do while (associated(iofield))
      call io_comp_realize_field(geniostate, iofield, importState, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      iofield => iofield%nfld
    enddo

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
    type(genio_field), pointer :: iofield
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

    ! query component for import and export states
    call NUOPC_ModelGet(genio, importState=importState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! reset import fields
    iofield => geniostate%imp_flds_head
    do while (associated(iofield))
      if (iofield%rlze) then
        call ESMF_FieldFill(iofield%efld, dataFillScheme="const", &
          const1=iofield%dflt, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
      endif
      iofield => iofield%nfld
    enddo

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
    type(genio_field), pointer :: iofield
    type(ESMF_Clock)           :: modelClock
    type(ESMF_State)           :: importState
    character(len=160)         :: clockString

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
    call ESMF_ClockPrint(modelClock, options="currTime", &
      unit=clockString, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    if (geniostate%verbosity .gt. 0) then
      ! write to standard out
      if (geniostate%myid .eq. geniostate%outid) then
        write(*,'(A,X,A)') trim(geniostate%cname)//": Model Advance",trim(clockString)
      endif

      ! sum import data from all PETs
      iofield => geniostate%imp_flds_head
      if (geniostate%myid .eq. geniostate%outid) then
        write(*,'(A)') trim(geniostate%cname)//": Import Fields"
        write(*,'(A,X,A25,X,A9,3(X,A9))') &
          trim(geniostate%cname)//":", "FIELD", &
          "COUNT", "MEAN", &
          "MIN", "MAX"
      endif
      do while (associated(iofield))
        call io_comp_field_diagnostics(geniostate, iofield, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        if (geniostate%myid .eq. geniostate%outid) then
          write(*,'(A,X,A25,X,I9,3(X,E9.2))') &
            trim(geniostate%cname)//":", trim(iofield%stdn), &
            int(iofield%gsum(2)), iofield%gavg, &
            iofield%gmin(1), iofield%gmax(1)
        endif
        iofield => iofield%nfld
      enddo
    endif ! diagnostic

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
    type(genio_field), pointer :: iofield
    type(ESMF_Clock)           :: modelClock
    type(ESMF_State)           :: importState
    character(len=160)         :: clockString
    integer                    :: fc
    type(ESMF_Field), pointer  :: fl(:)
    type(ESMF_FieldBundle)     :: fb
    character(ESMF_MAXSTR)     :: fieldName

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
    call ESMF_ClockPrint(modelClock, options="currTime", &
      unit=clockString, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    if (geniostate%verbosity .gt. 0) then 
      ! write to standard out
      if (geniostate%myid .eq. geniostate%outid) then
        write(*,'(A,X,A)') trim(geniostate%cname)//": Model Finalize",trim(clockString)
      endif
          
      ! sum import data from all PETs
      iofield => geniostate%imp_flds_head
      if (geniostate%myid .eq. geniostate%outid) then
        write(*,'(A)') trim(geniostate%cname)//": Import Fields"
        write(*,'(A,X,A25,X,A9,3(X,A9))') &
          trim(geniostate%cname)//":", "FIELD", &
          "COUNT", "MEAN", &
          "MIN", "MAX"
      endif
      do while (associated(iofield))
        call io_comp_field_diagnostics(geniostate, iofield, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        if (geniostate%myid .eq. geniostate%outid) then
          write(*,'(A,X,A25,X,I9,3(X,E9.2))') &
            trim(geniostate%cname)//":", trim(iofield%stdn), &
            int(iofield%gsum(2)), iofield%gavg, &
            iofield%gmin(1), iofield%gmax(1)
        endif
        iofield => iofield%nfld                 
      enddo
    endif ! diagnostic

    ! write final import and export states
    if (geniostate%write_final) then
      call NUOPC_GetStateMemberCount(importState, fieldCount=fc, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      if (fc .gt. 0) then
        nullify(fl)
        call NUOPC_GetStateMemberLists(importState, fieldList=fl, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        fb = ESMF_FieldBundleCreate(fieldList=fl, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call ESMF_FieldBundleWrite(fb, &
          fileName=trim(geniostate%cname)//"_final_import.nc", &
          overwrite=.true., rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call ESMF_FieldBundleDestroy(fb, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        deallocate(fl)
      endif
    endif

    ! remove import fields from importState and destroy
    do while (associated(geniostate%imp_flds_head))
      iofield => geniostate%imp_flds_head
      geniostate%imp_flds_head => iofield%nfld
      call ESMF_FieldGet(iofield%efld, name=fieldName, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call ESMF_StateRemove(importState, (/fieldName/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call ESMF_FieldDestroy(iofield%efld, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      deallocate(iofield, stat=stat)
      if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
        msg=trim(geniostate%cname)//': Memory deallocation failed.', &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)) return
      nullify(iofield)
    enddo
    geniostate%imp_flds_tail => null()

    ! destroy grid
    call ESMF_GridDestroy(geniostate%grid, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! destroy config
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
