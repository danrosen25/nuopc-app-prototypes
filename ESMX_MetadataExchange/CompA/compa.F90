!==============================================================================
! Earth System Modeling Framework
! Copyright (c) 2002-2025, University Corporation for Atmospheric Research,
! Massachusetts Institute of Technology, Geophysical Fluid Dynamics
! Laboratory, University of Michigan, National Centers for Environmental
! Prediction, Los Alamos National Laboratory, Argonne National Laboratory,
! NASA Goddard Space Flight Center.
! Licensed under the University of Illinois-NCSA License.
!==============================================================================

#define RC_CHECK(rc) \
if (ESMF_LogFoundError(rcToCheck=rc, \
msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) \
return

module COMPA

  !-----------------------------------------------------------------------------
  ! COMPA: The Absolute Worst Atmosphere Stand-in
  !-----------------------------------------------------------------------------

  use ESMF
  use NUOPC
  use NUOPC_Model, &
    modelSS    => SetServices

  implicit none

  private

  public SetVM, SetServices

  !-----------------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------------

  subroutine SetServices(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_HConfig)        :: hconfig, hconfigNode
    character(80)             :: compLabel
    character(:), allocatable :: badKey
    logical                   :: isFlag

    rc = ESMF_SUCCESS

    ! derive from NUOPC_Model
    call NUOPC_CompDerive(model, modelSS, rc=rc)
    RC_CHECK(rc)

    ! specialize model
    call NUOPC_CompSpecialize(model, specLabel=label_Advertise, &
      specRoutine=Advertise, rc=rc)
    RC_CHECK(rc)
    call NUOPC_CompSpecialize(model, specLabel=label_RealizeProvided, &
      specRoutine=Realize, rc=rc)
    RC_CHECK(rc)
    call NUOPC_CompSpecialize(model, specLabel=label_Advance, &
      specRoutine=Advance, rc=rc)
    RC_CHECK(rc)

    ! validate config
    call ESMF_GridCompGet(model, name=compLabel, hconfigIsPresent=isFlag, rc=rc)
    RC_CHECK(rc)
    if (isFlag) then
      ! Hconfig object present
      call ESMF_GridCompGet(model, hconfig=hconfig, rc=rc)
      RC_CHECK(rc)
      hconfigNode = ESMF_HConfigCreateAt(hconfig, keyString=compLabel, rc=rc)
      RC_CHECK(rc)
      ! component responsibility to validate ESMX handled options here, and
      ! potentially locally handled options
      isFlag = ESMF_HConfigValidateMapKeys(hconfigNode, &
        vocabulary=["model        ", &  ! ESMX handled option
                    "petList      ", &  ! ESMX handled option
                    "ompNumThreads", &  ! ESMX handled option
                    "stdout       ", &  ! ESMX handled option
                    "stderr       ", &  ! ESMX handled option
                    "attributes   "  &  ! ESMX handled option
                   ], badKey=badKey, rc=rc)
      RC_CHECK(rc)
      if (.not.isFlag) then
        call ESMF_LogSetError(ESMF_RC_ARG_WRONG, &
          msg="An invalid key was found in config under "//trim(compLabel)// &
            " (maybe a typo?): "//badKey, &
          line=__LINE__, &
          file=__FILE__, rcToReturn=rc)
        return
      endif
    endif

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
    RC_CHECK(rc)

    ! importable field: sea_surface_temperature
    call NUOPC_Advertise(importState, &
      StandardName="sea_surface_temperature", name="sst", rc=rc)
    RC_CHECK(rc)

    ! importable field: sea_surface_salinity
    call NUOPC_Advertise(importState, &
      StandardName="sea_surface_salinity", name="sss", rc=rc)
    RC_CHECK(rc)

    ! exportable field: air_pressure_at_sea_level
    call NUOPC_Advertise(exportState, &
      StandardName="air_pressure_at_sea_level", name="pmsl", rc=rc)
    RC_CHECK(rc)

    ! exportable field: surface_net_downward_shortwave_flux
    call NUOPC_Advertise(exportState, &
      StandardName="surface_net_downward_shortwave_flux", name="rsns", rc=rc)
    RC_CHECK(rc)

    ! exportable field: precipitation_flux
    call NUOPC_Advertise(exportState, &
      StandardName="precipitation_flux", name="precip", rc=rc)
    RC_CHECK(rc)

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine Realize(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    ! local variables
    character(80)           :: compLabel
    type(ESMF_State)        :: importState, exportState
    type(ESMF_Field)        :: field
    type(ESMF_Grid)         :: gridIn, gridOut
    type(ESMF_Mesh)         :: meshIn, meshOut
    type(ESMF_LocStream)    :: locsIn, locsOut

    integer, parameter              :: totalNumPoints=100
    integer(ESMF_KIND_I4), pointer  :: mask(:)
    real(ESMF_KIND_R8), pointer     :: lon(:), lat(:)
    real(ESMF_KIND_R8), pointer     :: fptr(:)
    integer                         :: clb(1), cub(1), i
    type(ESMF_VM)                   :: vm
    type(ESMF_Info)                 :: einfo

    rc = ESMF_SUCCESS

    call NUOPC_CompGet(model, name=compLabel, rc=rc)
    RC_CHECK(rc)

    ! query for importState and exportState
    call NUOPC_ModelGet(model, importState=importState, &
      exportState=exportState, rc=rc)
    RC_CHECK(rc)

    ! create Grid objects for Fields
    gridIn = ESMF_GridCreateNoPeriDimUfrm(maxIndex=(/10, 100/), &
      minCornerCoord=(/10._ESMF_KIND_R8, 20._ESMF_KIND_R8/), &
      maxCornerCoord=(/100._ESMF_KIND_R8, 200._ESMF_KIND_R8/), &
      coordSys=ESMF_COORDSYS_CART, &
      staggerLocList=(/ESMF_STAGGERLOC_CENTER, ESMF_STAGGERLOC_CORNER/), &
      rc=rc)
    RC_CHECK(rc)
    gridOut = gridIn ! for now out same as in

    ! create Mesh objects for Fields
    meshIn = ESMF_MeshCreate(grid=gridIn, rc=rc)
    RC_CHECK(rc)
    meshOut = ESMF_MeshCreate(grid=gridOut, rc=rc)
    RC_CHECK(rc)

    ! create LocStream objects for Fields
    locsIn=ESMF_LocStreamCreate(name="Equatorial Measurements", &
        maxIndex=totalNumPoints, coordSys=ESMF_COORDSYS_SPH_DEG, &
        indexFlag=ESMF_INDEX_GLOBAL, rc=rc)
    RC_CHECK(rc)
    ! Add key data (internally allocating memory).
    call ESMF_LocStreamAddKey(locsIn,                 &
         keyName="ESMF:Lat",           &
         KeyTypeKind=ESMF_TYPEKIND_R8, &
         keyUnits="Degrees",           &
         keyLongName="Latitude", rc=rc)
    RC_CHECK(rc)
    call ESMF_LocStreamAddKey(locsIn,                 &
         keyName="ESMF:Lon",           &
         KeyTypeKind=ESMF_TYPEKIND_R8, &
         keyUnits="Degrees",           &
         keyLongName="Longitude", rc=rc)
    RC_CHECK(rc)
    call ESMF_LocStreamAddKey(locsIn,                 &
         keyName="ESMF:Mask",           &
         KeyTypeKind=ESMF_TYPEKIND_I4, &
         keyUnits="none",           &
         keyLongName="mask values", rc=rc)
    RC_CHECK(rc)
    ! Get coordinate memory
    call ESMF_LocStreamGetKey(locsIn,                 &
         localDE=0,                    &
         keyName="ESMF:Lat",           &
         farray=lat,                   &
         rc=rc)
    RC_CHECK(rc)
    call ESMF_LocStreamGetKey(locsIn,                 &
         localDE=0,                    &
         keyName="ESMF:Lon",           &
         farray=lon,                   &
         rc=rc)
    RC_CHECK(rc)
    ! Get mask memory
    call ESMF_LocStreamGetKey(locsIn,                 &
         localDE=0,                    &
         keyName="ESMF:Mask",           &
         farray=mask,                   &
         rc=rc)
    RC_CHECK(rc)
    locsOut = locsIn ! for now out same as in

    ! importable field on Grid: sea_surface_temperature
    field = ESMF_FieldCreate(name="sst", grid=gridIn, &
      typekind=ESMF_TYPEKIND_R8, rc=rc)
    RC_CHECK(rc)
    call NUOPC_Realize(importState, field=field, rc=rc)
    RC_CHECK(rc)

    ! importable field on Mesh: sea_surface_salinity
    field = ESMF_FieldCreate(name="sss", mesh=meshIn, &
      typekind=ESMF_TYPEKIND_R8, rc=rc)
    RC_CHECK(rc)
    call NUOPC_Realize(importState, field=field, rc=rc)
    RC_CHECK(rc)

    ! exportable field on Grid: air_pressure_at_sea_level
    field = ESMF_FieldCreate(name="pmsl", grid=gridOut, &
      typekind=ESMF_TYPEKIND_R8, rc=rc)
    RC_CHECK(rc)
    call NUOPC_Realize(exportState, field=field, rc=rc)
    RC_CHECK(rc)

    ! exportable field on Grid: surface_net_downward_shortwave_flux
    field = ESMF_FieldCreate(name="rsns", grid=gridOut, &
      typekind=ESMF_TYPEKIND_R8, rc=rc)
    RC_CHECK(rc)
    call NUOPC_Realize(exportState, field=field, rc=rc)
    RC_CHECK(rc)

    ! exportable field on Mesh: precipitation_flux
    field = ESMF_FieldCreate(name="precip", mesh=meshOut, &
      typekind=ESMF_TYPEKIND_R8, rc=rc)
    RC_CHECK(rc)
    call NUOPC_Realize(exportState, field=field, rc=rc)
    RC_CHECK(rc)


    call NUOPC_ModelTagAdd(model, tag="ATM_INI", rc=rc)
    RC_CHECK(rc)

!    ! set transfer info on export state
!    call ESMF_InfoGetFromHost(exportState, info=einfo, rc=rc)
!    RC_CHECK(rc)
!    call ESMF_InfoSet(einfo, "/NUOPC/Transfer/SOURCE", compLabel, rc=rc)
!    RC_CHECK(rc)
!    call ESMF_InfoSet(einfo, "/NUOPC/Transfer/TYPE", "ATMOSPHERE", rc=rc)
!    RC_CHECK(rc)

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine Advance(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    ! local variables
    character(80)               :: compLabel
    type(ESMF_Clock)            :: clock
    type(ESMF_State)            :: importState, exportState
    type(ESMF_VM)               :: vm
    integer                     :: currentSsiPe, localPet
    character(len=160)          :: msgString
    character(len=160)          :: currTimeString
    type(ESMF_Info)             :: einfo

    rc = ESMF_SUCCESS

    call NUOPC_CompGet(model, name=compLabel, rc=rc)
    RC_CHECK(rc)

    call NUOPC_ModelGet(model, modelClock=clock, &
      importState=importState, exportState=exportState, rc=rc)
    RC_CHECK(rc)

    call ESMF_GridCompGet(model, vm=vm, localPet=localPet, rc=rc)
    RC_CHECK(rc)

    call ESMF_VMLog(vm, prefix="COMPA Advance(): ", &
      logMsgFlag=ESMF_LOGMSG_INFO, rc=rc)
    RC_CHECK(rc)

    call NUOPC_ModelTagAdd(model, tag="ATM_ADV", rc=rc)
    RC_CHECK(rc)

!    ! set transfer info on export state
!    call ESMF_ClockPrint(clock, options="currTime", &
!      unit=currTimeString, rc=rc)
!    RC_CHECK(rc)
!    call ESMF_InfoGetFromHost(exportState, info=einfo, rc=rc)
!    RC_CHECK(rc)
!    call ESMF_InfoSet(einfo, "/NUOPC/Transfer/TIME", currTimeString, rc=rc)
!    RC_CHECK(rc)

!    ! print info from importState
!    call ESMF_InfoGetFromHost(importState, info=info, rc=rc)
!    RC_CHECK(rc)
!    call ESMF_InfoPrint(info, &
!      preString="COMPA_INFO:"//trim(currTimeString), rc=rc)
!    RC_CHECK(rc)

  end subroutine

  !-----------------------------------------------------------------------------

end module

#ifdef SHARED_OBJECT

! External access to SetVM
subroutine SetVM(comp, rc)
  use ESMF
  use COMPA, only: SetVMModule => SetVM
  type(ESMF_GridComp) :: comp
  integer, intent(out) :: rc
  call SetVMModule(comp, rc)
  RC_CHECK(rc)
end subroutine

! External access to SetServices
subroutine SetServices(comp, rc)
  use ESMF
  use COMPA, only: SetServicesModule => SetServices
  type(ESMF_GridComp) :: comp
  integer, intent(out) :: rc
  call SetServicesModule(comp, rc)
  RC_CHECK(rc)
end subroutine

#endif
