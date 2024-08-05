!==============================================================================
! Earth System Modeling Framework
! Copyright (c) 2002-2024, University Corporation for Atmospheric Research,
! Massachusetts Institute of Technology, Geophysical Fluid Dynamics
! Laboratory, University of Michigan, National Centers for Environmental
! Prediction, Los Alamos National Laboratory, Argonne National Laboratory,
! NASA Goddard Space Flight Center.
! Licensed under the University of Illinois-NCSA License.
!==============================================================================

module genio_mod_struct

  !-----------------------------------------------------------------------------
  ! Generic IO Data Structures
  !-----------------------------------------------------------------------------

  use ESMF

  implicit none

  private

  ! parameters
  public GENIO_FILV
  public GENIO_FCTRL_DISABLED
  public GENIO_FCTRL_OPTIONAL
  public GENIO_FCTRL_REQUIRED
  public GENIO_DFLT_NX
  public GENIO_DFLT_NY
  public GENIO_DFLT_NZ
  public GENIO_DFLT_MINX
  public GENIO_DFLT_MAXX
  public GENIO_DFLT_MINY
  public GENIO_DFLT_MAXY
  public GENIO_DFLT_CSYS
  ! data structures
  public genio_fld
  public genio_srclst
  public genio_outlst
  public genio_geom
  public genio_state
  public geniostate_wrap

  ! parameters
  real(ESMF_KIND_R8), parameter :: &
    GENIO_FILV = -1.0E34_ESMF_KIND_R8
  integer, parameter ::        &
    GENIO_FCTRL_DISABLED = -1, &
    GENIO_FCTRL_OPTIONAL =  0, &
    GENIO_FCTRL_REQUIRED =  1
  integer, parameter :: &
    GENIO_DFLT_NX = 64, &
    GENIO_DFLT_NY = 32, &
    GENIO_DFLT_NZ =  4
  real(ESMF_KIND_R8), parameter ::           &
    GENIO_DFLT_MINX = -126.000_ESMF_KIND_R8, &
    GENIO_DFLT_MAXX =  -64.000_ESMF_KIND_R8, &
    GENIO_DFLT_MINY =   22.000_ESMF_KIND_R8, &
    GENIO_DFLT_MAXY =   50.000_ESMF_KIND_R8
  type(ESMF_CoordSys_Flag), parameter :: &
    GENIO_DFLT_CSYS = ESMF_COORDSYS_SPH_DEG

  type genio_fld
    character(len=64)           :: name        = "uninitialized"
    type(ESMF_TypeKind_Flag)    :: ftyp        = ESMF_TYPEKIND_R8
    real(ESMF_KIND_R8)          :: dflt        = GENIO_FILV
    integer                     :: fdim        = 2
    type(genio_fld), pointer    :: sfld        => null()
    type(ESMF_Field)            :: efld
    real(ESMF_KIND_R8)          :: lmin(1)     = GENIO_FILV
    real(ESMF_KIND_R8)          :: gmin(1)     = GENIO_FILV
    real(ESMF_KIND_R8)          :: lmax(1)     = GENIO_FILV
    real(ESMF_KIND_R8)          :: gmax(1)     = GENIO_FILV
    real(ESMF_KIND_R8)          :: lsum(2)     = GENIO_FILV
    real(ESMF_KIND_R8)          :: gsum(2)     = GENIO_FILV
    real(ESMF_KIND_R8)          :: gavg        = GENIO_FILV
    logical, allocatable        :: msk2D(:,:)
    logical, allocatable        :: msk3D(:,:,:)
    real(ESMF_KIND_R8), pointer :: ptrR82D(:,:)   => null()
    real(ESMF_KIND_R8), pointer :: ptrR83D(:,:,:) => null()
  endtype genio_fld

  type genio_srclst
    character(len=64)            :: name = "SourceFields"
    type(genio_fld), allocatable :: fld(:)
  endtype genio_srclst

  type genio_outlst
    character(len=64)            :: name = "OutputFields"
    type(ESMF_TimeInterval)      :: outfreq
    type(ESMF_TimeInterval)      :: outtimr
    type(ESMF_TimeInterval)      :: accintv
    type(genio_fld), allocatable :: fld(:)
  endtype genio_outlst

  type genio_geom
    character(len=64)        :: name = "uninitialized"
    integer                  :: nx = GENIO_DFLT_NX
    integer                  :: ny = GENIO_DFLT_NY
    integer                  :: nz = GENIO_DFLT_NZ
    real(ESMF_KIND_R8)       :: minx = GENIO_DFLT_MINX
    real(ESMF_KIND_R8)       :: maxx = GENIO_DFLT_MAXX
    real(ESMF_KIND_R8)       :: miny = GENIO_DFLT_MINY
    real(ESMF_KIND_R8)       :: maxy = GENIO_DFLT_MAXY
    type(ESMF_CoordSys_Flag) :: csys = ESMF_COORDSYS_SPH_DEG
    type(ESMF_GeomType_Flag) :: gtyp = ESMF_GEOMTYPE_GRID
    type(ESMF_Grid)          :: grid
    type(ESMF_Mesh)          :: mesh
  endtype genio_geom

  type genio_state
    ! component information
    character(32)                   :: cname      = "GENIO"
    logical                         :: cfgPresent = .false.
    type(ESMF_HConfig)              :: cfg
    integer                         :: verbosity  =  0
    integer                         :: diagnostic =  0
    integer                         :: myid       = -1
    integer                         :: outid      =  0
    type(ESMF_VM)                   :: vm
    type(genio_geom), allocatable   :: geolst(:)
    type(genio_srclst)              :: srclst
    type(genio_outlst), allocatable :: outlst(:)
  endtype genio_state

  type geniostate_wrap
    type(genio_state), pointer :: ptr
  endtype geniostate_wrap

endmodule genio_mod_struct
