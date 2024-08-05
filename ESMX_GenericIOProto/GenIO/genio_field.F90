!==============================================================================
! Earth System Modeling Framework
! Copyright (c) 2002-2024, University Corporation for Atmospheric Research,
! Massachusetts Institute of Technology, Geophysical Fluid Dynamics
! Laboratory, University of Michigan, National Centers for Environmental
! Prediction, Los Alamos National Laboratory, Argonne National Laboratory,
! NASA Goddard Space Flight Center.
! Licensed under the University of Illinois-NCSA License.
!==============================================================================

module genio_mod_field

  !-----------------------------------------------------------------------------
  ! Generic IO Field
  !-----------------------------------------------------------------------------

  use ESMF
  use NUOPC
  use genio_mod_struct
  use genio_mod_util

  implicit none

  private

  public genio_fldlst_create
  public genio_fldlst_destroy
  public genio_fldlst_realize
  public genio_fldlst_remove
  public genio_fldlst_zero
  public genio_fldlst_reset
  public genio_fldlst_diagnostics
  public genio_fldlst_write

  interface genio_fldlst_create
    module procedure genio_srclst_create
  end interface

  interface genio_fldlst_destroy
    module procedure genio_srclst_destroy
  end interface

  interface genio_fldlst_realize
    module procedure genio_srclst_realize
  end interface

  interface genio_fldlst_remove
    module procedure genio_srclst_remove
  end interface

  interface genio_fldlst_zero
    module procedure genio_srclst_zero
  end interface

  interface genio_fldlst_reset
    module procedure genio_srclst_reset
  end interface

  interface genio_fldlst_diagnostics
    module procedure genio_srclst_diagnostics
  end interface

  interface genio_fldlst_write
    module procedure genio_srclst_write
  end interface

  contains

  !-----------------------------------------------------------------------------
  ! Generic IO Field Subroutines
  !-----------------------------------------------------------------------------

  function genio_fldcfg_disabled(flstcfg, fldnm, rc)
    ! return value
    logical :: genio_fldcfg_disabled
    ! arguments
    type(ESMF_HConfig), intent(in) :: flstcfg
    character(*), intent(in)       :: fldnm
    integer, intent(out)           :: rc
    ! local variables
    logical                        :: fldOpts
    type(ESMF_HConfig)             :: fieldcfg
    integer                        :: control
    character(:), allocatable      :: badKey

    rc = ESMF_SUCCESS
    genio_fldcfg_disabled = .false.

    ! access fieldcfg
    fldOpts = ESMF_HConfigIsDefined(flstcfg, &
      keyString=fldnm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    if (fldOpts) then
      fieldcfg = ESMF_HConfigCreateAt(flstcfg, &
        keyString=fldnm, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      control = genio_hconfig2control(fieldcfg, key="control", &
        defaultValue=GENIO_FCTRL_OPTIONAL, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      if (control .eq. GENIO_FCTRL_DISABLED) then
        genio_fldcfg_disabled = .true.
      else
        genio_fldcfg_disabled = .false.
      endif
      call ESMF_HConfigDestroy(fieldcfg, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    else
      genio_fldcfg_disabled = .false.
    endif
  endfunction genio_fldcfg_disabled

  !-----------------------------------------------------------------------------

  subroutine genio_state_check(state, flstcfg, rc)
    ! arguments
    type(ESMF_State), intent(in)   :: state
    type(ESMF_HConfig), intent(in) :: flstcfg
    integer, intent(out)           :: rc
    ! local variables
    type(ESMF_HConfig)                     :: fieldcfg
    type(ESMF_HConfigIter)                 :: flistcur
    type(ESMF_HConfigIter)                 :: flistbeg
    type(ESMF_HConfigIter)                 :: flistend
    character(:), allocatable              :: fname
    integer                                :: control
    type(ESMF_StateItem_Flag)              :: itemType

    rc = ESMF_SUCCESS

    flistbeg = ESMF_HConfigIterBegin(flstcfg, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    flistend = ESMF_HConfigIterEnd(flstcfg, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    flistcur = flistbeg
    do while (ESMF_HConfigIterLoop(flistcur, flistbeg, flistend, rc=rc))
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      fname = ESMF_HConfigAsStringMapKey(flistcur, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      fieldcfg = ESMF_HConfigCreateAt(flstcfg, keyString=fname, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      control = genio_hconfig2control(fieldcfg, key="control", &
        defaultValue=GENIO_FCTRL_OPTIONAL, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      if (control .eq. GENIO_FCTRL_REQUIRED) then
        call ESMF_StateGet(state, itemName=fname, itemType=itemType, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        if (itemType .ne. ESMF_STATEITEM_FIELD) then
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
            msg="GENIO: required field is not available - "//trim(fname), &
            line=__LINE__, file=__FILE__, rcToReturn=rc)
          return
        endif
      endif
      call ESMF_HConfigDestroy(fieldcfg, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      deallocate(fname)
    enddo
  endsubroutine genio_state_check

  !-----------------------------------------------------------------------------

  function genio_srclst_create(state, flstcfg, rc)
    ! return value
    type(genio_srclst) :: genio_srclst_create
    ! arguments
    type(ESMF_State), intent(inout)          :: state
    type(ESMF_HConfig), intent(in), optional :: flstcfg
    integer, intent(out)                     :: rc
    ! local variables
    integer                                :: stat
    integer                                :: itemCount
    integer                                :: fieldCount
    character(len=64),allocatable          :: itemNameList(:)
    type(ESMF_StateItem_Flag), allocatable :: itemTypeList(:)
    integer                                :: i
    integer                                :: n
    type(ESMF_Field)                       :: field
    logical                                :: disabled

    rc = ESMF_SUCCESS

    ! get items in state
    call ESMF_StateGet(state, itemCount=itemCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    allocate(itemNameList(itemCount), stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg='GENIO: Memory allocation failed.', &
      line=__LINE__, file=__FILE__, rcToReturn=rc)) return
    allocate(itemTypeList(itemCount), stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg='GENIO: Memory allocation failed.', &
      line=__LINE__, file=__FILE__, rcToReturn=rc)) return
    call ESMF_StateGet(state, itemNameList=itemNameList, &
      itemTypeList=itemTypeList, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! count fields in state (if not disabled)
    fieldCount = 0
    if (present(flstcfg)) then
      call genio_state_check(state, flstcfg, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      do n=1, itemCount
        if (itemTypeList(n) == ESMF_STATEITEM_FIELD) then
          disabled = genio_fldcfg_disabled(flstcfg, itemNameList(n), rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
          if (.not. disabled) fieldCount = fieldCount + 1
        endif
      enddo
    else
      do n=1, itemCount
        if (itemTypeList(n) == ESMF_STATEITEM_FIELD) then
          fieldCount = fieldCount + 1
        endif
      enddo
    endif

    ! create srclst from enabled fields
    allocate(genio_srclst_create%fld(fieldCount), stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg='GENIO: Memory allocation failed.', &
      line=__LINE__, file=__FILE__, rcToReturn=rc)) return
    i = 0
    do n=1, itemCount
      if (itemTypeList(n) == ESMF_STATEITEM_FIELD) then
        if (present(flstcfg)) then
          disabled = genio_fldcfg_disabled(flstcfg, itemNameList(n), rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
        else
          disabled = .false.
        endif
        if (.not. disabled) then
          i = i + 1
          genio_srclst_create%fld(i)%name = itemNameList(n)
          genio_srclst_create%fld(i)%fdim = 2
          call ESMF_StateGet(state, itemName=itemNameList(n), &
            field=genio_srclst_create%fld(i)%efld, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
        else
          call ESMF_StateRemove(state, itemNameList=(/itemNameList(n)/), &
            rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
        endif
      endif
    enddo

  endfunction genio_srclst_create

  !-----------------------------------------------------------------------------

  subroutine genio_srclst_destroy(srclst, rc)
    ! arguments
    type(genio_srclst), intent(inout) :: srclst
    integer, intent(out)              :: rc
    ! local variables
    integer :: stat
    integer :: i

    rc = ESMF_SUCCESS

    do i=1, size(srclst%fld)
      call ESMF_FieldDestroy(srclst%fld(i)%efld, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    enddo

    deallocate(srclst%fld, stat=stat)
    if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
      msg='GENIO: Memory deallocation failed.', &
      line=__LINE__, file=__FILE__, rcToReturn=rc)) return

  endsubroutine genio_srclst_destroy

  !-----------------------------------------------------------------------------

  subroutine genio_srclst_realize(state, srclst, rc)
    ! arguments
    type(ESMF_State), intent(inout)   :: state
    type(genio_srclst), intent(inout) :: srclst
    integer, intent(out)              :: rc
    ! local variables
    integer :: stat
    integer :: i

    rc = ESMF_SUCCESS

    do i=1, size(srclst%fld)
      call NUOPC_Realize(state, fieldName=srclst%fld(i)%name, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call ESMF_FieldGet(srclst%fld(i)%efld, &
        dimCount=srclst%fld(i)%fdim, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      if(srclst%fld(i)%fdim .eq. 3) then
        call ESMF_FieldGet(srclst%fld(i)%efld, &
          farrayPtr=srclst%fld(i)%ptrR82D, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
      elseif(srclst%fld(i)%fdim .eq. 2) then
        call ESMF_FieldGet(srclst%fld(i)%efld, &
          farrayPtr=srclst%fld(i)%ptrR82D, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
      else
        call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
          msg="GENIO: field dimension - "//trim(srclst%fld(i)%name), &
          line=__LINE__, file=__FILE__, rcToReturn=rc)
        return
      endif
      call ESMF_FieldFill(srclst%fld(i)%efld, dataFillScheme="const", &
        const1=0.0_ESMF_KIND_R8, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    enddo ! srclst

  endsubroutine genio_srclst_realize

  !-----------------------------------------------------------------------------

  subroutine genio_srclst_remove(state, srclst, rc)
    ! arguments
    type(ESMF_State), intent(inout)   :: state
    type(genio_srclst), intent(inout) :: srclst
    integer, intent(out)              :: rc
    ! local variables
    integer :: stat
    integer :: i

    rc = ESMF_SUCCESS

    do i=1, size(srclst%fld)
      call ESMF_StateRemove(state, (/srclst%fld(i)%name/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    enddo ! srclst

  endsubroutine genio_srclst_remove

  !-----------------------------------------------------------------------------

  subroutine genio_srclst_zero(srclst, rc)
    ! arguments
    type(genio_srclst), intent(inout) :: srclst
    integer, intent(out)              :: rc
    ! local variables
    integer :: i

    rc = ESMF_SUCCESS

    do i=1, size(srclst%fld)
      call ESMF_FieldFill(srclst%fld(i)%efld, dataFillScheme="const", &
        const1=0.0_ESMF_KIND_R8, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    enddo ! srclst

  endsubroutine genio_srclst_zero

  !-----------------------------------------------------------------------------

  subroutine genio_srclst_reset(srclst, rc)
    ! arguments
    type(genio_srclst), intent(inout) :: srclst
    integer, intent(out)              :: rc
    ! local variables
    integer :: i

    rc = ESMF_SUCCESS

    do i=1, size(srclst%fld)
      call ESMF_FieldFill(srclst%fld(i)%efld, dataFillScheme="const", &
        const1=srclst%fld(i)%dflt, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    enddo ! srclst

  endsubroutine genio_srclst_reset

  !-----------------------------------------------------------------------------

  subroutine genio_srclst_diagnostics(srclst, myid, outid, vm, cname, label, &
  timeString, rc)
    ! arguments
    type(genio_srclst), intent(inout) :: srclst
    integer, intent(in)               :: myid
    integer, intent(in)               :: outid
    type(ESMF_VM), intent(in)         :: vm
    character(*), intent(in)          :: cname
    character(*), intent(in)          :: label
    character(*), intent(in)          :: timeString
    integer, intent(out)              :: rc
    ! local variables
    real(ESMF_KIND_R8) :: eval
    integer            :: i

    rc = ESMF_SUCCESS
    do i = 1, size(srclst%fld)
      call genio_field_stats(srclst%fld(i))
    enddo

    if (myid .eq. outid) then
      write(*,'(A,X,A)') trim(cname)//": "//trim(label), &
        trim(timeString)
      write(*,'(A)') trim(cname)//": "//trim(srclst%name)
      write(*,'(A,X,A25,X,A9,3(X,A9))') &
        trim(cname)//":", "FIELD", "COUNT", "MEAN", "MIN", "MAX"
      do i = 1, size(srclst%fld)
        write(*,'(A,X,A25,X,I9,3(X,E9.2))') &
          trim(cname)//":", trim(srclst%fld(i)%name), &
          int(srclst%fld(i)%gsum(2)), srclst%fld(i)%gavg, &
          srclst%fld(i)%gmin(1), srclst%fld(i)%gmax(1)
      enddo
    endif

    contains

    subroutine genio_field_stats(fld)
      ! arguments
      type(genio_fld), intent(inout) :: fld

      if(fld%fdim .eq. 3) then
        eval = epsilon(fld%ptrR83D)
        fld%lsum(1)=sum(fld%ptrR83D,abs(fld%ptrR83D-GENIO_FILV)>eval)
        fld%lsum(2)=count(abs(fld%ptrR83D-GENIO_FILV)>eval)
        fld%lmin(1)=minval(fld%ptrR83D,abs(fld%ptrR83D-GENIO_FILV)>eval)
        fld%lmax(1)=maxval(fld%ptrR83D,abs(fld%ptrR83D-GENIO_FILV)>eval)
      elseif(fld%fdim .eq. 2) then
        eval = epsilon(fld%ptrR82D)
        fld%lsum(1)=sum(fld%ptrR82D,abs(fld%ptrR82D-GENIO_FILV)>eval)
        fld%lsum(2)=count(abs(fld%ptrR82D-GENIO_FILV)>eval)
        fld%lmin(1)=minval(fld%ptrR82D,abs(fld%ptrR82D-GENIO_FILV)>eval)
        fld%lmax(1)=maxval(fld%ptrR82D,abs(fld%ptrR82D-GENIO_FILV)>eval)
      else
        call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
          msg=trim(cname)//": field dimension - "// &
              trim(fld%name), &
          line=__LINE__, file=__FILE__, rcToReturn=rc)
        return
      endif
      call ESMF_VMReduce(vm, sendData=fld%lsum, &
        recvData=fld%gsum, count=2, &
        reduceflag=ESMF_REDUCE_SUM, rootPet=outid, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call ESMF_VMReduce(vm=vm, sendData=fld%lmin, &
        recvData=fld%gmin, count=1, &
        reduceflag=ESMF_REDUCE_MIN, rootPet=outid, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call ESMF_VMReduce(vm=vm, sendData=fld%lmax, &
        recvData=fld%gmax, count=1, &
        reduceflag=ESMF_REDUCE_MAX, rootPet=outid, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      if (myid .eq. outid) then
        ! calculate average
        if(fld%gsum(2) .lt. 1) then
          fld%gavg = 0.0_ESMF_KIND_R8
        else
          fld%gavg = fld%gsum(1) / fld%gsum(2)
        endif
      else
        fld%gsum = (/GENIO_FILV, 0.0_ESMF_KIND_R8/)
        fld%gmin = GENIO_FILV
        fld%gmax = GENIO_FILV
        fld%gavg = 0.0_ESMF_KIND_R8
      endif

    endsubroutine genio_field_stats

  endsubroutine genio_srclst_diagnostics

  !-----------------------------------------------------------------------------

  subroutine genio_outlst_write(outlst, fileName, rc)
    ! arguments
    type(genio_outlst), intent(in) :: outlst
    character(*), intent(in)       :: fileName
    integer, intent(out)           :: rc
    ! local variables
    type(ESMF_FieldBundle) :: fb
    integer                :: i

    rc = ESMF_SUCCESS

    fb = ESMF_FieldBundleCreate(rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    do i = 1, size(outlst%fld)
      call ESMF_FieldBundleAdd(fb, fieldList=(/outlst%fld(i)%efld/), rc=rc)
    enddo
    call ESMF_FieldBundleWrite(fb, &
      fileName=filename, &
      overwrite=.true., rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call ESMF_FieldBundleDestroy(fb, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

  endsubroutine genio_outlst_write

  !-----------------------------------------------------------------------------

  subroutine genio_srclst_write(srclst, fileName, rc)
    ! arguments
    type(genio_srclst), intent(in) :: srclst
    character(*), intent(in)       :: fileName
    integer, intent(out)           :: rc
    ! local variables
    type(ESMF_FieldBundle) :: fb
    integer                :: i

    rc = ESMF_SUCCESS

    fb = ESMF_FieldBundleCreate(rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    do i = 1, size(srclst%fld)
      call ESMF_FieldBundleAdd(fb, fieldList=(/srclst%fld(i)%efld/), rc=rc)
    enddo
    call ESMF_FieldBundleWrite(fb, &
      fileName=filename, &
      overwrite=.true., rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call ESMF_FieldBundleDestroy(fb, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

  endsubroutine genio_srclst_write

endmodule genio_mod_field
