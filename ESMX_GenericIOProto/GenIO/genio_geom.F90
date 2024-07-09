!==============================================================================
! Earth System Modeling Framework
! Copyright (c) 2002-2024, University Corporation for Atmospheric Research,
! Massachusetts Institute of Technology, Geophysical Fluid Dynamics
! Laboratory, University of Michigan, National Centers for Environmental
! Prediction, Los Alamos National Laboratory, Argonne National Laboratory,
! NASA Goddard Space Flight Center.
! Licensed under the University of Illinois-NCSA License.
!==============================================================================

module genio_mod_geom

  !-----------------------------------------------------------------------------
  ! Generic IO Geom
  !-----------------------------------------------------------------------------

  use ESMF
  use genio_mod_param
  use genio_mod_struct
  use genio_mod_util

  implicit none

  private

  public io_comp_read_geom
  public io_comp_grid_diag

  contains

  !-----------------------------------------------------------------------------
  ! Generic IO Geom Subroutines
  !-----------------------------------------------------------------------------

  subroutine io_comp_read_geom(geniostate, rc)
    ! arguments
    type(genio_state), pointer, intent(inout) :: geniostate
    integer, intent(out)                      :: rc
    ! local variables
    logical            :: isPresent
    integer            :: stat
    logical            :: check
    type(ESMF_HConfig) :: geomcfg
    character(len=64)  :: cfgval
    character(:), allocatable :: badKey

    rc = ESMF_SUCCESS

    if (.not. associated(geniostate)) then
      call ESMF_LogSetError(ESMF_RC_PTR_NOTALLOC, &
        msg='GENIO: geniostate has not been associated', &
        line=__LINE__, file=__FILE__, rcToReturn=rc)
      return
    endif

    ! read geom configuration
    if (geniostate%cfgPresent) then
      isPresent = ESMF_HConfigIsDefined(geniostate%cfg, &
        keyString="geom", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    else
      isPresent = .false.
    endif
    if (isPresent) then
      ! access geom
      geomcfg = ESMF_HConfigCreateAt(geniostate%cfg, &
        keyString="geom", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      check = ESMF_HConfigValidateMapKeys(geomcfg, &
        vocabulary=["nx      ", &
                    "ny      ", &
                    "nz      ", &
                    "coordSys", &
                    "minx    ", &
                    "maxx    ", &
                    "miny    ", &
                    "maxy    "  &
                   ], badKey=badKey, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      if (.not. check) then
        call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
          msg=trim(geniostate%cname)//": unknown geom option key - "//badKey, &
        line=__LINE__,file=__FILE__, rcToReturn=rc)
        return
      endif
      ! dimensions
      geniostate%nx = genio_hconfig2i4(geomcfg, "nx", &
        defaultValue=64, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      geniostate%ny = genio_hconfig2i4(geomcfg, "ny", &
        defaultValue=32, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      geniostate%nz = genio_hconfig2i4(geomcfg, "nz", &
        defaultValue=4, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      ! coordinate system
      cfgval = genio_hconfig2str(geomcfg, "coordSys", &
        defaultValue="ESMF_COORDSYS_SPH_DEG", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      cfgval = ESMF_UtilStringUpperCase(cfgval, rc=rc)
      select case (cfgval)
        case ("ESMF_COORDSYS_CART")
          geniostate%coordSys = ESMF_COORDSYS_CART
        case ("ESMF_COORDSYS_SPH_DEG")
          geniostate%coordSys = ESMF_COORDSYS_SPH_DEG
        case ("ESMF_COORDSYS_SPH_RAD")
          geniostate%coordSys = ESMF_COORDSYS_SPH_RAD
        case default
          call ESMF_LogSetError(ESMF_RC_ARG_VALUE, &
            msg=trim(geniostate%cname)//': invalid value - coordSys', &
            line=__LINE__, file=__FILE__, rcToReturn=rc)
          return
      endselect
      ! coordinates
      geniostate%minx = genio_hconfig2r8(geomcfg, "minx", &
        defaultValue=-126.0_ESMF_KIND_R8, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      geniostate%maxx = genio_hconfig2r8(geomcfg, "maxx", &
        defaultValue=-64.0_ESMF_KIND_R8, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      geniostate%miny = genio_hconfig2r8(geomcfg, "miny", &
        defaultValue=22.0_ESMF_KIND_R8, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      geniostate%maxy = genio_hconfig2r8(geomcfg, "maxy", &
        defaultValue=50.0_ESMF_KIND_R8, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call ESMF_HConfigDestroy(geomcfg, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    endif ! geomcfg

  endsubroutine io_comp_read_geom

  !-----------------------------------------------------------------------------

  subroutine io_comp_grid_diag(geniostate, fileName, overwrite, status, timeslice, &
  iofmt, relaxedflag, rc)
    ! arguments
    type(genio_state), pointer, intent(in)           :: geniostate
    character(len=*), intent(in), optional           :: fileName
    logical, intent(in), optional                    :: overwrite
    type(ESMF_FileStatus_Flag), intent(in), optional :: status
    integer, intent(in), optional                    :: timeslice
    type(ESMF_IOFmt_Flag), intent(in), optional      :: iofmt
    logical, intent(in), optional                    :: relaxedflag
    integer, intent(out)                             :: rc
    ! local variables
    logical                 :: ioCapable
    logical                 :: doItFlag
    character(len=64)       :: lfileName
    character(len=64)       :: gridName
    type(ESMF_Array)        :: array
    type(ESMF_ArrayBundle)  :: arraybundle
    logical                 :: isPresent
    integer                 :: dimCount
    integer                 :: dimIndex
    integer,allocatable     :: coordDimCount(:)
    integer                 :: coordDimMax
    integer                 :: stat
    logical                 :: hasCorners

    rc = ESMF_SUCCESS

    if (.not. associated(geniostate)) then
      call ESMF_LogSetError(ESMF_RC_PTR_NOTALLOC, &
        msg='GENIO: geniostate has not been associated', &
        line=__LINE__, file=__FILE__, rcToReturn=rc)
      return
    endif

    ioCapable = (ESMF_IO_PIO_PRESENT .and. &
      (ESMF_IO_NETCDF_PRESENT .or. ESMF_IO_PNETCDF_PRESENT))

    doItFlag = .true. ! default
    if (present(relaxedFlag)) then
      doItFlag = .not.relaxedflag .or. (relaxedflag.and.ioCapable)
    endif

    if (doItFlag) then

      if (present(fileName)) then
        lfileName = trim(fileName)
      else
        call ESMF_GridGet(geniostate%grid, name=gridName, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        lfileName = trim(gridName)//".nc"
      endif

      arraybundle = ESMF_ArrayBundleCreate(rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return

      ! -- centers --

      call ESMF_GridGetCoord(geniostate%grid, staggerLoc=ESMF_STAGGERLOC_CENTER, &
        isPresent=isPresent, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      if (isPresent) then
        call ESMF_GridGetCoord(geniostate%grid, coordDim=1, &
          staggerLoc=ESMF_STAGGERLOC_CENTER, array=array, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call ESMF_ArraySet(array, name="lon_center", rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call ESMF_ArrayBundleAdd(arraybundle,(/array/),rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call ESMF_GridGetCoord(geniostate%grid, coordDim=2, &
          staggerLoc=ESMF_STAGGERLOC_CENTER, array=array, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call ESMF_ArraySet(array, name="lat_center", rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call ESMF_ArrayBundleAdd(arraybundle,(/array/),rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
      endif

      ! -- corners --

      call ESMF_GridGetCoord(geniostate%grid, staggerLoc=ESMF_STAGGERLOC_CORNER, &
        isPresent=hasCorners, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      if (hasCorners) then
        call ESMF_GridGetCoord(geniostate%grid, coordDim=1, &
          staggerLoc=ESMF_STAGGERLOC_CORNER, array=array, rc=rc)
        if (.not. ESMF_LogFoundError(rcToCheck=rc)) then
          call ESMF_ArraySet(array, name="lon_corner", rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
          call ESMF_ArrayBundleAdd(arraybundle,(/array/),rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
        endif
        call ESMF_GridGetCoord(geniostate%grid, coordDim=2, &
          staggerLoc=ESMF_STAGGERLOC_CORNER, array=array, rc=rc)
        if (.not. ESMF_LogFoundError(rcToCheck=rc)) then
          call ESMF_ArraySet(array, name="lat_corner", rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
          call ESMF_ArrayBundleAdd(arraybundle,(/array/),rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
        endif
      endif

      ! -- mask --

      call ESMF_GridGetItem(geniostate%grid, itemflag=ESMF_GRIDITEM_MASK, &
        staggerLoc=ESMF_STAGGERLOC_CENTER, isPresent=isPresent, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      if (isPresent) then
        call ESMF_GridGetItem(geniostate%grid, staggerLoc=ESMF_STAGGERLOC_CENTER, &
          itemflag=ESMF_GRIDITEM_MASK, array=array, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call ESMF_ArraySet(array, name="mask", rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call ESMF_ArrayBundleAdd(arraybundle,(/array/),rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
      endif

      ! -- area --

      call ESMF_GridGetItem(geniostate%grid, itemflag=ESMF_GRIDITEM_AREA, &
        staggerLoc=ESMF_STAGGERLOC_CENTER, isPresent=isPresent, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      if (isPresent) then
        call ESMF_GridGetItem(geniostate%grid, staggerLoc=ESMF_STAGGERLOC_CENTER, &
          itemflag=ESMF_GRIDITEM_AREA, array=array, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call ESMF_ArraySet(array, name="area", rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call ESMF_ArrayBundleAdd(arraybundle,(/array/),rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
      endif

      call ESMF_ArrayBundleWrite(arraybundle, &
        fileName=trim(lfileName),rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return

      call ESMF_ArrayBundleDestroy(arraybundle,rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    endif
  endsubroutine io_comp_grid_diag

  !-----------------------------------------------------------------------------

endmodule genio_mod_geom
