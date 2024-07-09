!==============================================================================
! Earth System Modeling Framework
! Copyright (c) 2002-2024, University Corporation for Atmospheric Research,
! Massachusetts Institute of Technology, Geophysical Fluid Dynamics
! Laboratory, University of Michigan, National Centers for Environmental
! Prediction, Los Alamos National Laboratory, Argonne National Laboratory,
! NASA Goddard Space Flight Center.
! Licensed under the University of Illinois-NCSA License.
!==============================================================================

module genio_mod_param

  !-----------------------------------------------------------------------------
  ! Generic IO Parameters
  !-----------------------------------------------------------------------------

  use ESMF

  implicit none

  private

  public filv
  public p_disabled
  public p_required
  public p_optional

  ! parameters
  real(ESMF_KIND_R8), parameter :: filv = -1.0E34_ESMF_KIND_R8
  integer, parameter :: p_disabled = -1
  integer, parameter :: p_required =  1
  integer, parameter :: p_optional =  0

endmodule genio_mod_param
