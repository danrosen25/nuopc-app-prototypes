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
  use genio_mod_param
  use genio_mod_struct
  use genio_mod_util

  implicit none

  private

  public io_comp_add_fields
  public io_comp_realize_field
  public io_comp_field_diagnostics

  contains

  !-----------------------------------------------------------------------------
  ! Generic IO Field Subroutines
  !-----------------------------------------------------------------------------

  subroutine io_comp_add_fields(state, geniostate, rc)
    ! arguments
    type(ESMF_State), intent(in)              :: state
    type(genio_state), pointer, intent(inout) :: geniostate
    integer, intent(out)                      :: rc
    ! local variables
    integer                                :: stat
    integer                                :: n
    integer                                :: itemCount
    character(len=64),allocatable          :: itemNameList(:)
    type(ESMF_StateItem_Flag), allocatable :: itemTypeList(:)
    type(ESMF_Field)                       :: field
    logical                                :: fldOptsList
    logical                                :: fldOpts
    character(len=64)                      :: fldName
    logical                                :: check
    type(ESMF_HConfig)                     :: flistcfg
    type(ESMF_HConfig)                     :: fieldcfg
    type(ESMF_HConfigIter)                 :: flistcur
    type(ESMF_HConfigIter)                 :: flistbeg
    type(ESMF_HConfigIter)                 :: flistend
    type(genio_field), pointer             :: iofield
    character(:), allocatable              :: badKey

    rc = ESMF_SUCCESS

    if (.not. associated(geniostate)) then
      call ESMF_LogSetError(ESMF_RC_PTR_NOTALLOC, &
        msg='GENIO: geniostate has not been associated', &
        line=__LINE__, file=__FILE__, rcToReturn=rc)
      return
    endif

    ! read field options
    if (geniostate%cfgPresent) then
      fldOptsList = ESMF_HConfigIsDefined(geniostate%cfg, &
        keyString="fieldOptions", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    else
      fldOptsList = .false.
    endif
    if (fldOptsList) then
      flistcfg = ESMF_HConfigCreateAt(geniostate%cfg, &
        keyString="fieldOptions", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    endif

    call ESMF_StateGet(state, itemCount=itemCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    allocate(itemNameList(itemCount), stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg='GENIO: Memory allocation failed.', &
      line=__LINE__, &
      file=__FILE__, &
      rcToReturn=rc)) return
    allocate(itemTypeList(itemCount), stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg='GENIO: Memory allocation failed.', &
      line=__LINE__, &
      file=__FILE__, &
      rcToReturn=rc)) return

    call ESMF_StateGet(state, itemNameList=itemNameList, &
      itemTypeList=itemTypeList, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    do n=1, itemCount
      if (itemTypeList(n) == ESMF_STATEITEM_FIELD) then
        nullify(iofield)
        allocate(iofield, stat=stat)
        if (ESMF_LogFoundAllocError(statusToCheck=stat, &
          msg=trim(geniostate%cname)//': Memory allocation failed.', &
          line=__LINE__, &
          file=__FILE__, &
          rcToReturn=rc)) return
        ! field name
        iofield%stdn = itemNameList(n)
        ! configurable options
        iofield%prio = p_optional
        ! basic initialization
        iofield%fdim = 2
        iofield%lsum = (/filv, 0.0_ESMF_KIND_R8/)
        iofield%lmin = filv
        iofield%lmax = filv
        iofield%gsum = (/filv, 0.0_ESMF_KIND_R8/)
        iofield%gmin = filv
        iofield%gmax = filv
        iofield%gavg = filv
        iofield%dflt = filv
        iofield%nfld => null()
        if (.not. associated(geniostate%imp_flds_head)) then
          geniostate%imp_flds_head => iofield
          geniostate%imp_flds_tail => iofield
        else
          geniostate%imp_flds_tail%nfld => iofield
          geniostate%imp_flds_tail => iofield
        endif
        ! read field options
        if (fldOptsList) then
          ! access fieldcfg
          fldOpts = ESMF_HConfigIsDefined(flistcfg, &
            keyString=iofield%stdn, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
          if (fldOpts) then
            fieldcfg = ESMF_HConfigCreateAt(flistcfg, &
              keyString=iofield%stdn, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, file=__FILE__)) return
            check = ESMF_HConfigValidateMapKeys(fieldcfg, &
              vocabulary=["priority" &
                         ], badKey=badKey, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, file=__FILE__)) return
            if (.not. check) then
              call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
                msg=trim(geniostate%cname)//": ("//iofield%stdn//")"// &
                    " unknown fieldOption key - "//badKey, &
              line=__LINE__,file=__FILE__, rcToReturn=rc)
              return
            endif
            iofield%prio = genio_hconfig2priority(fieldcfg, &
              defaultValue=p_optional, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, file=__FILE__)) return
            call ESMF_HConfigDestroy(fieldcfg, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, file=__FILE__)) return
          endif ! fldOpts
        endif ! fldOptsList
      endif ! ESMF_STATEITEM_FIELD
    enddo ! itemCount

    if (fldOptsList) then  
      call ESMF_HConfigDestroy(flistcfg, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    endif ! fldOptsList

  endsubroutine io_comp_add_fields

  !-----------------------------------------------------------------------------

  subroutine io_comp_realize_field(geniostate, iofield, state, rc)
    ! arguments
    type(genio_state), pointer, intent(inout) :: geniostate
    type(genio_field), pointer, intent(inout) :: iofield
    type(ESMF_State), intent(inout)  :: state
    integer, intent(out)             :: rc
    ! local variables
    integer :: stat

    rc = ESMF_SUCCESS

    if (.not. associated(geniostate)) then
      call ESMF_LogSetError(ESMF_RC_PTR_NOTALLOC, &
        msg='GENIO: geniostate has not been associated', &
        line=__LINE__, file=__FILE__, rcToReturn=rc)
      return
    endif

    if (.not. associated(iofield)) then
      call ESMF_LogSetError(ESMF_RC_MEM_ALLOCATE, &
        msg=trim(geniostate%cname)//": iofield error", &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)
      return
    endif

    call NUOPC_Realize(state, fieldName=iofield%stdn, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    if (associated(iofield%efld)) then
      call ESMF_LogSetError(ESMF_RC_MEM_ALLOCATE, &
        msg=trim(geniostate%cname)//": ESMF_Field error - "//trim(iofield%stdn), &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)
      return
    endif
    allocate(iofield%efld, stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg=trim(geniostate%cname)//': Memory allocation failed.', &
      line=__LINE__, &
      file=__FILE__, &
      rcToReturn=rc)) return
    call ESMF_StateGet(state, itemName=iofield%stdn, field=iofield%efld, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call ESMF_FieldGet(iofield%efld, dimCount=iofield%fdim, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    if(iofield%fdim .eq. 3) then
      call ESMF_FieldGet(iofield%efld, farrayPtr=iofield%ptr3, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    elseif(iofield%fdim .eq. 2) then
      call ESMF_FieldGet(iofield%efld, farrayPtr=iofield%ptr2, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    else
      call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
        msg=trim(geniostate%cname)//": field dimension - "//trim(iofield%stdn), &
        line=__LINE__, file=__FILE__, rcToReturn=rc)
      return
    endif
    call ESMF_FieldFill(iofield%efld, dataFillScheme="const", &
      const1=0.0_ESMF_KIND_R8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    iofield%rlze = .true.

  endsubroutine io_comp_realize_field

  !-----------------------------------------------------------------------------

  subroutine io_comp_field_diagnostics(geniostate, iofield, rc)
    ! arguments
    type(genio_state), pointer, intent(in)    :: geniostate
    type(genio_field), pointer, intent(inout) :: iofield
    integer, intent(out)                      :: rc
    ! local variables
    real(ESMF_KIND_R8) :: eval

    rc = ESMF_SUCCESS

    if (.not. associated(geniostate)) then
      call ESMF_LogSetError(ESMF_RC_PTR_NOTALLOC, &
        msg='GENIO: geniostate has not been associated', &
        line=__LINE__, file=__FILE__, rcToReturn=rc)
      return
    endif

    if (iofield%rlze) then
      if(iofield%fdim .eq. 3) then
        eval = epsilon(iofield%ptr3)
        iofield%lsum(1)=sum(iofield%ptr3,abs(iofield%ptr3-filv)>eval)
        iofield%lsum(2)=count(abs(iofield%ptr3-filv)>eval)
        iofield%lmin(1)=minval(iofield%ptr3,abs(iofield%ptr3-filv)>eval)
        iofield%lmax(1)=maxval(iofield%ptr3,abs(iofield%ptr3-filv)>eval)
      elseif(iofield%fdim .eq. 2) then
        eval = epsilon(iofield%ptr2)
        iofield%lsum(1)=sum(iofield%ptr2,abs(iofield%ptr2-filv)>eval)
        iofield%lsum(2)=count(abs(iofield%ptr2-filv)>eval)
        iofield%lmin(1)=minval(iofield%ptr2,abs(iofield%ptr2-filv)>eval)
        iofield%lmax(1)=maxval(iofield%ptr2,abs(iofield%ptr2-filv)>eval)
      else
        call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
          msg=trim(geniostate%cname)//": field dimension - "//trim(iofield%stdn), &
          line=__LINE__, file=__FILE__, rcToReturn=rc)
        return
      endif
      call ESMF_VMReduce(vm=geniostate%vm, sendData=iofield%lsum, &
        recvData=iofield%gsum, count=2, &
        reduceflag=ESMF_REDUCE_SUM, rootPet=geniostate%outid, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call ESMF_VMReduce(vm=geniostate%vm, sendData=iofield%lmin, &
        recvData=iofield%gmin, count=1, &
        reduceflag=ESMF_REDUCE_MIN, rootPet=geniostate%outid, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call ESMF_VMReduce(vm=geniostate%vm, sendData=iofield%lmax, &
        recvData=iofield%gmax, count=1, &
        reduceflag=ESMF_REDUCE_MAX, rootPet=geniostate%outid, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      if (geniostate%myid .eq. geniostate%outid) then
        ! calculate average
        if(iofield%gsum(2) .lt. 1) then
          iofield%gavg = 0.0_ESMF_KIND_R8
        else
          iofield%gavg = iofield%gsum(1) / iofield%gsum(2)
        endif
      endif
    else
      iofield%gsum = (/filv, 0.0_ESMF_KIND_R8/)
      iofield%gmin = filv
      iofield%gmax = filv
      iofield%gavg = 0.0_ESMF_KIND_R8
    endif
  endsubroutine io_comp_field_diagnostics

endmodule genio_mod_field
