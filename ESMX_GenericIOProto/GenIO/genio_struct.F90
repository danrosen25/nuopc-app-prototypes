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
  use genio_mod_param

  implicit none

  private

  public genio_field
  public genio_state
  public geniostate_wrap

  ! derived types
  type genio_field
    character(len=64)           :: stdn        = "dummy"
    real(ESMF_KIND_R8)          :: dflt        = filv
    integer                     :: prio        = p_optional
    logical                     :: rlze        = .false.
    integer                     :: fdim        = 2
    real(ESMF_KIND_R8)          :: lmin(1)     = filv
    real(ESMF_KIND_R8)          :: gmin(1)     = filv
    real(ESMF_KIND_R8)          :: lmax(1)     = filv
    real(ESMF_KIND_R8)          :: gmax(1)     = filv
    real(ESMF_KIND_R8)          :: lsum(2)     = filv
    real(ESMF_KIND_R8)          :: gsum(2)     = filv
    real(ESMF_KIND_R8)          :: gavg        = filv
    type(ESMF_Field), pointer   :: efld        => null()
    real(ESMF_KIND_R8), pointer :: ptr2(:,:)   => null()
    real(ESMF_KIND_R8), pointer :: ptr3(:,:,:) => null()
    type(genio_field), pointer  :: nfld        => null()
  endtype genio_field

  type genio_state
    ! component information
    character(32)      :: cname       = "GENIO"
    logical            :: cfgPresent  = .false.
    type(ESMF_HConfig) :: cfg
    integer            :: verbosity   =  0
    integer            :: diagnostic  =  0
    logical            :: write_final = .true.
    integer            :: myid        = -1
    integer            :: outid       =  0
    type(ESMF_VM)      :: vm
    ! grid information
    integer                  :: nx = 64
    integer                  :: ny = 32
    integer                  :: nz = 4
    real(ESMF_KIND_R8)       :: minx = -126.000_ESMF_KIND_R8
    real(ESMF_KIND_R8)       :: maxx =  -64.000_ESMF_KIND_R8
    real(ESMF_KIND_R8)       :: miny =   22.000_ESMF_KIND_R8
    real(ESMF_KIND_R8)       :: maxy =   50.000_ESMF_KIND_R8
    type(ESMF_CoordSys_Flag) :: coordSys = ESMF_COORDSYS_SPH_DEG
    type(ESMF_Grid)          :: grid
    ! field information
    type(genio_field), pointer :: imp_flds_head => null()
    type(genio_field), pointer :: imp_flds_tail => null()
  endtype genio_state

  type geniostate_wrap
    type(genio_state), pointer :: ptr
  endtype geniostate_wrap

endmodule genio_mod_struct
