!==============================================================================
! Earth System Modeling Framework
! Copyright (c) 2002-2024, University Corporation for Atmospheric Research,
! Massachusetts Institute of Technology, Geophysical Fluid Dynamics
! Laboratory, University of Michigan, National Centers for Environmental
! Prediction, Los Alamos National Laboratory, Argonne National Laboratory,
! NASA Goddard Space Flight Center.
! Licensed under the University of Illinois-NCSA License.
!==============================================================================

module OCN

  !-----------------------------------------------------------------------------
  ! OCN Component.
  !-----------------------------------------------------------------------------

  use ESMF
  use NUOPC
  use NUOPC_Model, &
    modelSS     => SetServices

  implicit none

  private

  public SetServices

  real(ESMF_KIND_R8),parameter :: zeroValue    =      0.0_ESMF_KIND_R8
  real(ESMF_KIND_R8),parameter :: missingValue = 999999.0_ESMF_KIND_R8
  real(ESMF_KIND_R8),parameter :: islandLoc(4) = (/ 289.5_ESMF_KIND_R8, &
                                                    295.5_ESMF_KIND_R8, &
                                                      9.5_ESMF_KIND_R8, &
                                                     25.5_ESMF_KIND_R8 /)
  type(ESMF_Time) :: prevCurrTime

  !-----------------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------------

  subroutine SetServices(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    rc = ESMF_SUCCESS

    ! derive from NUOPC_Model
    call NUOPC_CompDerive(model, modelSS, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! specialize model
    call NUOPC_CompSpecialize(model, specLabel=label_Advertise, &
      specRoutine=Advertise, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompSpecialize(model, specLabel=label_RealizeProvided, &
      specRoutine=Realize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompSpecialize(model, specLabel=label_SetClock, &
      specRoutine=SetClock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompSpecialize(model, specLabel=label_CheckImport, &
      specRoutine=CheckImport, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompSpecialize(model, specLabel=label_Advance, &
      specRoutine=Advance, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompSpecialize(model, specLabel=label_Finalize, &
      specRoutine=Finalize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine Advertise(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_State)        :: importState, exportState

    rc = ESMF_SUCCESS

    ! query for importState and exportState
    call NUOPC_ModelGet(model, importState=importState, &
      exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! Disabling the following macro, e.g. renaming to WITHIMPORTFIELDS_disable,
    ! will result in a model component that does not advertise any importable
    ! Fields. Use this if you want to drive the model independently.
#define WITHIMPORTFIELDS
#ifdef WITHIMPORTFIELDS
    ! importable field: air_pressure_at_sea_level
    call NUOPC_Advertise(importState, &
      StandardName="air_pressure_at_sea_level", name="pmsl", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! importable field: surface_net_downward_shortwave_flux
    call NUOPC_Advertise(importState, &
      StandardName="surface_net_downward_shortwave_flux", name="rsns", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
#endif

    ! exportable field: sea_surface_temperature
    call NUOPC_Advertise(exportState, &
      StandardName="sea_surface_temperature", name="sst", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine Realize(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_State)              :: importState, exportState
    type(ESMF_TimeInterval)        :: stabilityTimeStep
    type(ESMF_Field)               :: field
    type(ESMF_Grid)                :: gridIn
    type(ESMF_Grid)                :: gridOut
    integer                        :: tlb(2), tub(2)
    real(ESMF_KIND_R8), pointer    :: lon_fptr(:)
    real(ESMF_KIND_R8), pointer    :: lat_fptr(:)
    integer(ESMF_KIND_I4), pointer :: msk_fptr(:,:)
    integer                        :: i,j

    rc = ESMF_SUCCESS

    ! query for importState and exportState
    call NUOPC_ModelGet(model, importState=importState, &
      exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! create a Grid object for Fields
    gridIn = ESMF_GridCreateNoPeriDimUfrm(maxIndex=(/15, 16/), &
      minCornerCoord=(/260._ESMF_KIND_R8,  0._ESMF_KIND_R8/), &
      maxCornerCoord=(/355._ESMF_KIND_R8, 45._ESMF_KIND_R8/), &
      coordSys=ESMF_COORDSYS_SPH_DEG, staggerLocList=(/ESMF_STAGGERLOC_CENTER/), &
      rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! get grid coordinates
    call ESMF_GridGetCoord(gridIn, coordDim=1, &
      staggerLoc=ESMF_STAGGERLOC_CENTER, farrayPtr=lon_fptr, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_GridGetCoord(gridIn, coordDim=2, &
      staggerLoc=ESMF_STAGGERLOC_CENTER, farrayPtr=lat_fptr, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! add mask and island
    call ESMF_GridAddItem(gridIn, itemflag=ESMF_GRIDITEM_MASK, &
      staggerLoc=ESMF_STAGGERLOC_CENTER, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_GridGetItem(gridIn, itemflag=ESMF_GRIDITEM_MASK, &
      staggerLoc=ESMF_STAGGERLOC_CENTER, &
      totalLBound=tlb, totalUBound=tub, farrayPtr=msk_fptr, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    msk_fptr = 0
    do j=tlb(2), tub(2)
    do i=tlb(1), tub(1)
      if ((lon_fptr(i).ge.islandLoc(1)) .AND. &
          (lon_fptr(i).le.islandLoc(2)) .AND. &
          (lat_fptr(j).ge.islandLoc(3)) .AND. &
          (lat_fptr(j).le.islandLoc(4)) ) then
        msk_fptr(i,j) = 1
      endif
    enddo
    enddo

    gridOut = gridIn ! for now out same as in

#ifdef WITHIMPORTFIELDS
    ! importable field: air_pressure_at_sea_level
    field = ESMF_FieldCreate(name="pmsl", grid=gridIn, &
      typekind=ESMF_TYPEKIND_R8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_Realize(importState, field=field, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! fill import field pmsl with missingValue
    call ESMF_FieldFill(field, dataFillScheme="const", &
      const1=missingValue, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! importable field: surface_net_downward_shortwave_flux
    field = ESMF_FieldCreate(name="rsns", grid=gridIn, &
      typekind=ESMF_TYPEKIND_R8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_Realize(importState, field=field, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! fill import field rsns with missingValue
    call ESMF_FieldFill(field, dataFillScheme="const", &
      const1=missingValue, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
#endif

    ! exportable field: sea_surface_temperature
    field = ESMF_FieldCreate(name="sst", grid=gridOut, &
      typekind=ESMF_TYPEKIND_R8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_Realize(exportState, field=field, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! fill export field sst with 270
    call ESMF_FieldFill(field, dataFillScheme="const", &
      const1=270.0_ESMF_KIND_R8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine SetClock(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)              :: clock
    type(ESMF_TimeInterval)       :: stabilityTimeStep

    rc = ESMF_SUCCESS

    ! query the Component for its clock, importState and exportState
    call NUOPC_ModelGet(model, modelClock=clock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! initialize internal clock
    ! here: parent Clock and stability timeStep determine actual model timeStep
    !TODO: stabilityTimeStep should be read in from configuation
    !TODO: or computed from internal Grid information
    call ESMF_TimeIntervalSet(stabilityTimeStep, m=15, rc=rc) ! 15 minute steps
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompSetClock(model, clock, stabilityTimeStep, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_ClockGet(clock, startTime=prevCurrTime, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine CheckImport(model, rc)
    type(ESMF_GridComp)   :: model
    integer, intent(out)  :: rc

    ! local variables
    type(ESMF_Clock)                :: clock
    type(ESMF_Time)                 :: currTime
    type(ESMF_Time)                 :: invalidTime
    type(ESMF_State)                :: importState
    logical                         :: timeCheck
    integer                         :: i

    rc = ESMF_SUCCESS

    ! query for clock and importState
    call ESMF_GridCompGet(model, clock=clock, importState=importState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! get the current time out of the clock
    call ESMF_ClockGet(clock, currTime=currTime, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! set up invalid time (by convention)
    call ESMF_TimeSet(invalidTime, yy=99999999, mm=01, dd=01, &
      h=00, m=00, s=00, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    timeCheck = NUOPC_IsAtTime(importState, invalidTime, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    if (timeCheck) then
      ! The field is at invalidTime
      call ESMF_LogWrite("OCN: detected import field at invalidTime", &
        ESMF_LOGMSG_INFO, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    else
      ! The field is NOT at invalidTime
      ! check that Fields in the importState show correct timestamp
      timeCheck = NUOPC_IsAtTime(importState, currTime, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      if (.not.timeCheck) then
        call ESMF_LogWrite("OCN: NUOPC INCOMPATIBILITY DETECTED: "// &
          "Import Field not at current time", &
          ESMF_LOGMSG_WARNING)
      endif

    endif

  end subroutine


  !-----------------------------------------------------------------------------

  subroutine Advance(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    ! local variables
    integer                     :: localPet
    type(ESMF_Clock)            :: clock
    type(ESMF_State)            :: importState, exportState
    type(ESMF_Time)             :: currTime
    type(ESMF_TimeInterval)     :: timeStep
    character(len=160)          :: msgString

#define NUOPC_TRACE__OFF
#ifdef NUOPC_TRACE
    call ESMF_TraceRegionEnter("OCN:Advance")
#endif

    rc = ESMF_SUCCESS

    ! query the Component for its localPet
    call ESMF_GridCompGet(model, localPet=localPet, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! query for clock, importState and exportState
    call NUOPC_ModelGet(model, modelClock=clock, importState=importState, &
      exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_ClockGet(clock, currTime=currTime, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    if (currTime .lt. prevCurrTime) then
      msgString="OCN: Clock rewind caught"
      call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      if (localPet .eq. 0) print *, trim(msgString)
    endif

    call ESMF_ClockPrint(clock, options="currTime", &
      preString="OCN: current time = ", unit=msgString, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    if (localPet .eq. 0) print *, trim(msgString)

    call ESMF_ClockGet(clock, currTime=prevCurrTime, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

#ifdef NUOPC_TRACE
    call ESMF_TraceRegionExit("OCN:Advance")
#endif
  end subroutine

  !-----------------------------------------------------------------------------

  subroutine Finalize(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc
    ! local variables
    type(ESMF_VM)                           :: vm
    character(ESMF_MAXSTR)                  :: name
    character(len=8)                        :: value
    logical                                 :: zeroValues
    logical                                 :: missingValues
    integer                                 :: localPet
    type(ESMF_State)                        :: importState
    integer                                 :: itemCount, i
    integer                                 :: fieldCount
    character(len=80), allocatable          :: itemNameList(:)
    type(ESMF_StateItem_Flag), allocatable  :: itemTypeList(:)
    type(ESMF_Field)                        :: field
    real(ESMF_KIND_R8), pointer             :: farrayPtr(:,:)
    integer                                 :: lclZero(1)
    integer                                 :: gblZero(1)
    integer                                 :: lclMissing(1)
    integer                                 :: gblMissing(1)

    rc = ESMF_SUCCESS

    ! query the Component for its vm
    call ESMF_GridCompGet(model, vm=vm, name=name, localPet=localPet, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! get diagnostic settings
    call ESMF_AttributeGet(model, name="zeroValues", &
      value=value, defaultValue="false", &
      convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    value = ESMF_UtilStringLowerCase(value, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    select case (value)
      case ('true','t','yes')
        zeroValues=.true.
      case default
        zeroValues=.false.
    endselect

    call ESMF_AttributeGet(model, name="missingValues", &
      value=value, defaultValue="false", &
      convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    value = ESMF_UtilStringLowerCase(value, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    select case (value)
      case ('true','t','yes')
        missingValues=.true.
      case default
        missingValues=.false.
    endselect

    ! query the Component for its importState
    call NUOPC_ModelGet(model, importState=importState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! fill import state with zeros
    call ESMF_StateGet(importState, nestedFlag=.true., &
      itemCount=itemCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    allocate(itemNameList(itemCount), itemTypeList(itemCount))
    call ESMF_StateGet(importState, nestedFlag=.true., &
      itemNameList=itemNameList, itemTypeList=itemTypeList, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    do i=1, itemCount
      if (itemTypeList(i)==ESMF_STATEITEM_FIELD) then
        call ESMF_StateGet(importState, field=field, itemName=itemNameList(i), &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        call ESMF_FieldGet(field, farrayPtr=farrayPtr, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        lclZero(1) = COUNT(farrayPtr(:,:).eq.zeroValue)
        lclMissing(1) = COUNT(farrayPtr(:,:).eq.missingValue)
        call ESMF_VMReduce(vm, lclZero, gblZero, &
          reduceflag=ESMF_REDUCE_SUM, count=1, rootPet=0, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        call ESMF_VMReduce(vm, lclMissing, gblMissing, &
          reduceflag=ESMF_REDUCE_SUM, count=1, rootPet=0, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        if (localPet.eq.0) then
          select case (itemNameList(i))
            case ('pmsl','rsns')
              write (*,"(A,A,A,A,I0)") trim(name),": ", &
                trim(itemNameList(i)), &
                " zero values = ", gblZero(1)
              write (*,"(A,A,A,A,I0)") trim(name),": ", &
                trim(itemNameList(i)), &
                " missing values = ", gblMissing(1)
              if ((gblZero(1).gt.0).neqv.zeroValues) then
                write (*,"(A,A,A,L1,A)") "ERROR: ",trim(name), &
                  " zero values must be ", zeroValues, &
                  ", see [runconfig] file "
                call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
                  msg="Zero values does not match config settings", &
                  line=__LINE__, &
                  file=__FILE__, &
                  rcToReturn=rc)
                return ! bail out
              elseif ((gblMissing(1).gt.0).neqv.missingValues) then
                write (*,"(A,A,A,L1,A)") "ERROR: ",trim(name), &
                  " missing values must be ", missingValues, &
                  ", see [runconfig] file "
                call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
                  msg="Missing values does not match config settings", &
                  line=__LINE__, &
                  file=__FILE__, &
                  rcToReturn=rc)
                return ! bail out
              endif
            case default
              write (*,"(A)") "Field is unknown "//trim(itemNameList(i))
              call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
                msg="Field is unknown "//trim(itemNameList(i)), &
                line=__LINE__, &
                file=__FILE__, &
                rcToReturn=rc)
              return ! bail out
          endselect
        endif
      endif
    enddo
    deallocate(itemNameList, itemTypeList)

    call ESMF_LogWrite("OCN: Field check passed.", ESMF_LOGMSG_INFO, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  end subroutine

  !-----------------------------------------------------------------------------

end module
