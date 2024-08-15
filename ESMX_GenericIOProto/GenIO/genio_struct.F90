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
  use genio_mod_params

  implicit none

  private

  ! data structures
  public genio_fld
  public genio_srclst
  public genio_outlst
  public genio_geom
  public genio_state
  public geniostate_wrap

  type genio_fld
    character(len=64)           :: name        = "uninitialized"
    character(len=32)           :: unts        = "uninitialized"
    logical                     :: lavg        = .true.
    type(ESMF_Time)             :: prvt
    logical                     :: strt        = .false.
    real(ESMF_KIND_R8)          :: timw        = 0
    integer                     :: ftyp        = GENIO_ERROR
    integer                     :: otyp        = GENIO_INST
    character(len=6)            :: ostr        = "(----)"
    real(ESMF_KIND_R8)          :: dflt        = GENIO_FILV
    type(ESMF_Field)            :: srcf
    type(ESMF_Field)            :: prvf
    type(ESMF_Field)            :: sumf
    type(ESMF_Field)            :: outf
    character(len=19)           :: tstr
    real(ESMF_KIND_R8)          :: lmin(1)     = GENIO_FILV
    real(ESMF_KIND_R8)          :: gmin(1)     = GENIO_FILV
    real(ESMF_KIND_R8)          :: lmax(1)     = GENIO_FILV
    real(ESMF_KIND_R8)          :: gmax(1)     = GENIO_FILV
    real(ESMF_KIND_R8)          :: lsum(2)     = GENIO_FILV
    real(ESMF_KIND_R8)          :: gsum(2)     = GENIO_FILV
    real(ESMF_KIND_R8)          :: gavg        = GENIO_FILV
    real(ESMF_KIND_R8)          :: eval
  endtype genio_fld

  type genio_srclst
    character(len=64)            :: name  = "uninitialized"
    logical                      :: write = .false.
    logical                      :: diagn = .false.
    integer                      :: slice = 0
    type(genio_fld), allocatable :: fld(:)
  endtype genio_srclst

  type genio_outlst
    character(len=64)            :: name  = "uninitialized"
    logical                      :: wfinl = .false.
    logical                      :: diagn = .false.
    integer                      :: slice = 0
    type(ESMF_TimeInterval)      :: outfrq
    type(ESMF_TimeInterval)      :: accwnd
    type(ESMF_Time)              :: prvtme
    type(genio_fld), allocatable :: fld(:)
    type(ESMF_TimeInterval)      :: zrofrq
    real(ESMF_KIND_R8)           :: eval
  endtype genio_outlst

  type genio_geom
    character(len=64)        :: name = "uninitialized"
    integer                  :: nx   = GENIO_DFLT_NX
    integer                  :: ny   = GENIO_DFLT_NY
    integer                  :: nz   = GENIO_DFLT_NZ
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
    character(32)           :: cname      = "GENIO"
    logical                 :: cfgPresent = .false.
    type(ESMF_HConfig)      :: cfg
    integer                 :: verbosity  =  0
    integer                 :: diagnostic =  0
    integer                 :: myid       = -1
    integer                 :: outid      =  0
    type(ESMF_TimeInterval) :: zerotime
    type(ESMF_VM)           :: vm
    type(genio_geom)        :: geom
    type(genio_srclst)      :: srclst
    type(genio_outlst)      :: outlst
  endtype genio_state

  type geniostate_wrap
    type(genio_state), pointer :: ptr
  endtype geniostate_wrap

endmodule genio_mod_struct
