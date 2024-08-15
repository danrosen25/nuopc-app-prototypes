!==============================================================================
! Earth System Modeling Framework
! Copyright (c) 2002-2024, University Corporation for Atmospheric Research,
! Massachusetts Institute of Technology, Geophysical Fluid Dynamics
! Laboratory, University of Michigan, National Centers for Environmental
! Prediction, Los Alamos National Laboratory, Argonne National Laboratory,
! NASA Goddard Space Flight Center.
! Licensed under the University of Illinois-NCSA License.
!==============================================================================

module genio_mod_params

  !-----------------------------------------------------------------------------
  ! Generic IO Parameters
  !-----------------------------------------------------------------------------

  use ESMF

  implicit none

  private

  ! parameters
  public GENIO_ERROR
  public GENIO_FILV
  public GENIO_I4D1
  public GENIO_I4D2
  public GENIO_I4D3
  public GENIO_R8D1
  public GENIO_R8D2
  public GENIO_R8D3
  public GENIO_SKIP
  public GENIO_INST
  public GENIO_MEAN
  public GENIO_ACCM
  public GENIO_DFLT_NX
  public GENIO_DFLT_NY
  public GENIO_DFLT_NZ
  public GENIO_DFLT_MINX
  public GENIO_DFLT_MAXX
  public GENIO_DFLT_MINY
  public GENIO_DFLT_MAXY
  public GENIO_DFLT_CSYS
  public GENIO_DFLT_OTYP

  ! parameters
  integer, parameter :: &
    GENIO_ERROR = -999999
  real(ESMF_KIND_R8), parameter :: &
    GENIO_FILV = -1.0E34_ESMF_KIND_R8
  integer, parameter :: &
    GENIO_I4D1 =  1, &
    GENIO_I4D2 =  2, &
    GENIO_I4D3 =  3, &
    GENIO_R8D1 =  4, &
    GENIO_R8D2 =  5, &
    GENIO_R8D3 =  6
  integer, parameter :: &
    GENIO_SKIP =  7, &
    GENIO_INST =  8, &
    GENIO_MEAN =  9, &
    GENIO_ACCM = 10
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
  integer, parameter :: &
    GENIO_DFLT_OTYP = GENIO_INST

endmodule genio_mod_params
