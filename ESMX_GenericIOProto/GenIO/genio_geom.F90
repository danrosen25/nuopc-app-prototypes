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
  use genio_mod_struct
  use genio_mod_params
  use genio_mod_util

  implicit none

  private

  public genio_geolst_create
  public genio_geolst_destroy
  public genio_geolst_write

  contains

  !-----------------------------------------------------------------------------
  ! Generic IO Geom Subroutines
  !-----------------------------------------------------------------------------

  function genio_geolst_create(geolcfg, dflts, rc)
    ! return value
    type(genio_geom), allocatable :: genio_geolst_create(:)
    ! arguments
    type(ESMF_HConfig), intent(in) :: geolcfg
    type(genio_dflts), intent(in)  :: dflts
    integer, intent(out)           :: rc
    ! local variables
    integer            :: i
    integer            :: gcnt
    integer            :: stat
    type(ESMF_HConfig) :: geomcfg
    type(ESMF_Grid)    :: grid
    type(ESMF_Mesh)    :: mesh

    rc = ESMF_SUCCESS

    gcnt = ESMF_HConfigGetSize(geolcfg, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    allocate(genio_geolst_create(gcnt), stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg='GENIO: Memory allocation failed.', &
      line=__LINE__, file=__FILE__, rcToReturn=rc)) return

    do i=1, gcnt
        geomcfg = ESMF_HConfigCreateAt(geolcfg, index=i, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        genio_geolst_create(i)%name = genio_hconfig2str(geomcfg, &
          key="name", rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        genio_geolst_create(i)%gtyp = genio_hconfig2gtyp(geomcfg, &
          key="geomtype", rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        if (genio_geolst_create(i)%gtyp .eq. ESMF_GEOMTYPE_GRID) then
          grid = genio_grid_create(genio_geolst_create(i)%name, &
            gridcfg=geomcfg, dflts=dflts, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
          genio_geolst_create(i)%geom = ESMF_GeomCreate(grid, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
        elseif (genio_geolst_create(i)%gtyp .eq. ESMF_GEOMTYPE_MESH) then
          mesh = genio_mesh_create(genio_geolst_create(i)%name, &
            meshcfg=geomcfg, dflts=dflts, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
          genio_geolst_create(i)%geom = ESMF_GeomCreate(mesh, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
        else
          call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
            msg="GENIO: geom type not supported", &
            line=__LINE__, file=__FILE__, rcToReturn=rc)
          return
        endif
        call ESMF_HConfigDestroy(geomcfg, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
    enddo

  endfunction genio_geolst_create

  !-----------------------------------------------------------------------------

  function genio_grid_create(name, gridcfg, dflts, rc)
    ! return value
    type(ESMF_Grid) :: genio_grid_create
    ! arguments
    character(*), intent(in)       :: name
    type(ESMF_HConfig), intent(in) :: gridcfg
    type(genio_dflts), intent(in)  :: dflts
    integer, intent(out)           :: rc
    ! local variables
    integer                         :: i
    integer                         :: stat
    logical                         :: check
    character(:), allocatable       :: badKey
    type(ESMF_CoordSys_Flag)        :: coordSys
    integer                         :: dimCount
    integer, allocatable            :: maxIndex(:)
    real(ESMF_KIND_R8), allocatable :: minCornerCoord(:)
    real(ESMF_KIND_R8), allocatable :: maxCornerCoord(:)

    rc = ESMF_SUCCESS

    genio_grid_create = ESMF_GridEmptyCreate(rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! check configuration
    check = ESMF_HConfigValidateMapKeys(gridcfg, &
      vocabulary=["name          ", &
                  "geomtype      ", &
                  "coordSys      ", &
                  "maxIndex      ", &
                  "minCornerCoord", &
                  "maxCornerCoord"  &
                 ], badKey=badKey, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    if (.not. check) then
      call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
        msg="GENIO: unknown option in defaults "//badKey, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)
      return
    endif

    ! coordSys
    coordSys = genio_hconfig2csys(gridcfg, key="coordSys", &
      defaultValue=dflts%csys, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    ! maxIndex
    dimCount = ESMF_HConfigGetSize(gridcfg, keyString="maxIndex", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    allocate(maxIndex(dimCount), stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg='GENIO: Memory allocation failed.', &
      line=__LINE__, file=__FILE__, rcToReturn=rc)) return
    maxIndex = ESMF_HConfigAsI4Seq(gridcfg, keyString="maxIndex", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    ! minCornerCoord
    dimCount = ESMF_HConfigGetSize(gridcfg, keyString="minCornerCoord", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    allocate(minCornerCoord(dimCount), stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg='GENIO: Memory allocation failed.', &
      line=__LINE__, file=__FILE__, rcToReturn=rc)) return
    minCornerCoord = ESMF_HConfigAsR8Seq(gridcfg, keyString="minCornerCoord", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    ! maxCornerCoord
    dimCount = ESMF_HConfigGetSize(gridcfg, keyString="maxCornerCoord", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    allocate(maxCornerCoord(dimCount), stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg='GENIO: Memory allocation failed.', &
      line=__LINE__, file=__FILE__, rcToReturn=rc)) return
    maxCornerCoord = ESMF_HConfigAsR8Seq(gridcfg, keyString="maxCornerCoord", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    ! ESMF_GridCreate
    genio_grid_create = ESMF_GridCreateNoPeriDimUfrm(maxIndex=maxIndex, &
      minCornerCoord=minCornerCoord, maxCornerCoord=maxCornerCoord, &
      coordSys=coordSys, name=name, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    deallocate(maxIndex, stat=stat)
    deallocate(minCornerCoord, stat=stat)
    deallocate(maxCornerCoord, stat=stat)

  endfunction genio_grid_create

  !-----------------------------------------------------------------------------

  function genio_mesh_create(name, meshcfg, dflts, rc)
    ! return value
    type(ESMF_Mesh) :: genio_mesh_create
    ! arguments
    character(*), intent(in)       :: name
    type(ESMF_HConfig), intent(in) :: meshcfg
    type(genio_dflts), intent(in)  :: dflts
    integer, intent(out)           :: rc
    ! local variables
    integer                         :: i
    integer                         :: stat
    logical                         :: check
    character(:), allocatable       :: badKey
    type(ESMF_CoordSys_Flag)        :: coordSys
    character(len=ESMF_MAXSTR)      :: filename
    type(ESMF_FileFormat_Flag)      :: fileformat

    rc = ESMF_SUCCESS

    genio_mesh_create = ESMF_MeshEmptyCreate(rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! check configuration
    check = ESMF_HConfigValidateMapKeys(meshcfg, &
      vocabulary=["name          ", &
                  "geomtype      ", &
                  "coordSys      ", &
                  "filename      ", &
                  "fileformat    "  &
                 ], badKey=badKey, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    if (.not. check) then
      call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
        msg="GENIO: unknown option in defaults "//badKey, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)
      return
    endif

    ! coordSys
    coordSys = genio_hconfig2csys(meshcfg, key="coordSys", &
      defaultValue=dflts%csys, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    ! filename
    filename = ESMF_HConfigAsString(meshcfg, keyString="filename", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    ! fileformat
    fileformat = genio_hconfig2ffmt(meshcfg, key="fileformat", &
      defaultValue=dflts%ffmt, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! ESMF_MeshCreate
    genio_mesh_create = ESMF_MeshCreate(filename=filename, &
      fileformat=fileformat, coordSys=coordSys, name=name, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

  endfunction genio_mesh_create

  !-----------------------------------------------------------------------------

  subroutine genio_geolst_destroy(geolst, rc)
    ! arguments
    type(genio_geom), allocatable, intent(inout) :: geolst(:)
    integer, intent(out)                         :: rc
    ! local variables
    integer :: i
    integer :: stat

    rc = ESMF_SUCCESS

    do i=1, size(geolst)
      call ESMF_GeomDestroy(geolst(i)%geom, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    enddo
    deallocate(geolst, stat=stat)

  endsubroutine genio_geolst_destroy

  !-----------------------------------------------------------------------------

  subroutine genio_geolst_write(geolst, overwrite, status, timeslice, iofmt, &
  relaxedflag, rc)
    ! arguments
    type(genio_geom), allocatable, intent(in)        :: geolst(:)
    logical, intent(in), optional                    :: overwrite
    type(ESMF_FileStatus_Flag), intent(in), optional :: status
    integer, intent(in), optional                    :: timeslice
    type(ESMF_IOFmt_Flag), intent(in), optional      :: iofmt
    logical, intent(in), optional                    :: relaxedflag
    integer, intent(out)                             :: rc
    ! local variables
    integer         :: i
    type(ESMF_Grid) :: grid

    rc = ESMF_SUCCESS

    do i=1, size(geolst)
      if (geolst(i)%gtyp .eq. ESMF_GEOMTYPE_GRID) then
        call ESMF_GeomGet(geolst(i)%geom, grid=grid, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call genio_grid_write(grid, overwrite=overwrite, status=status, &
          timeslice=timeslice, iofmt=iofmt, relaxedflag=relaxedflag, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
      endif
    enddo

  endsubroutine genio_geolst_write

  !-----------------------------------------------------------------------------

  subroutine genio_grid_write(grid, overwrite, status, timeslice, iofmt, &
  relaxedflag, rc)
    ! arguments
    type(ESMF_Grid), intent(in)                      :: grid
    logical, intent(in), optional                    :: overwrite
    type(ESMF_FileStatus_Flag), intent(in), optional :: status
    integer, intent(in), optional                    :: timeslice
    type(ESMF_IOFmt_Flag), intent(in), optional      :: iofmt
    logical, intent(in), optional                    :: relaxedflag
    integer, intent(out)                             :: rc
    ! local variables
    logical                 :: ioCapable
    logical                 :: doItFlag
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

    ioCapable = (ESMF_IO_PIO_PRESENT .and. &
      (ESMF_IO_NETCDF_PRESENT .or. ESMF_IO_PNETCDF_PRESENT))

    doItFlag = .true. ! default
    if (present(relaxedFlag)) then
      doItFlag = .not.relaxedflag .or. (relaxedflag.and.ioCapable)
    endif

    if (doItFlag) then

      call ESMF_GridGet(grid, name=gridName, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return

      arraybundle = ESMF_ArrayBundleCreate(rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return

      ! -- centers --

      call ESMF_GridGetCoord(grid, staggerLoc=ESMF_STAGGERLOC_CENTER, &
        isPresent=isPresent, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      if (isPresent) then
        call ESMF_GridGetCoord(grid, coordDim=1, &
          staggerLoc=ESMF_STAGGERLOC_CENTER, array=array, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call ESMF_ArraySet(array, name="lon_center", rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call ESMF_ArrayBundleAdd(arraybundle,(/array/),rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call ESMF_GridGetCoord(grid, coordDim=2, &
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

      call ESMF_GridGetCoord(grid, staggerLoc=ESMF_STAGGERLOC_CORNER, &
        isPresent=hasCorners, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      if (hasCorners) then
        call ESMF_GridGetCoord(grid, coordDim=1, &
          staggerLoc=ESMF_STAGGERLOC_CORNER, array=array, rc=rc)
        if (.not. ESMF_LogFoundError(rcToCheck=rc)) then
          call ESMF_ArraySet(array, name="lon_corner", rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
          call ESMF_ArrayBundleAdd(arraybundle,(/array/),rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
        endif
        call ESMF_GridGetCoord(grid, coordDim=2, &
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

      call ESMF_GridGetItem(grid, itemflag=ESMF_GRIDITEM_MASK, &
        staggerLoc=ESMF_STAGGERLOC_CENTER, isPresent=isPresent, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      if (isPresent) then
        call ESMF_GridGetItem(grid, staggerLoc=ESMF_STAGGERLOC_CENTER, &
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

      call ESMF_GridGetItem(grid, itemflag=ESMF_GRIDITEM_AREA, &
        staggerLoc=ESMF_STAGGERLOC_CENTER, isPresent=isPresent, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      if (isPresent) then
        call ESMF_GridGetItem(grid, staggerLoc=ESMF_STAGGERLOC_CENTER, &
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
        fileName=trim(gridName)//".nc",rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return

      call ESMF_ArrayBundleDestroy(arraybundle,rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    endif
  endsubroutine genio_grid_write

  !-----------------------------------------------------------------------------

endmodule genio_mod_geom
