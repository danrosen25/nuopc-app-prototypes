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
  use genio_mod_params
  use genio_mod_util

  implicit none

  private

  public genio_fldlst_create
  public genio_fldlst_destroy
  public genio_fldlst_realize
  public genio_fldlst_zero
  public genio_fldlst_reset
  public genio_fldlst_addsrc
  public genio_fldlst_print
  public genio_fldlst_alldata
  public genio_fldlst_alarm
  public genio_fldlst_write

  interface genio_fldlst_create
    module procedure genio_srclst_create
    module procedure genio_outlst_create
  end interface

  interface genio_fldlst_destroy
    module procedure genio_srclst_destroy
    module procedure genio_outlst_destroy
  end interface

  interface genio_fldlst_realize
    module procedure genio_srclst_realize
  end interface

  interface genio_fldlst_zero
    module procedure genio_srclst_zero
  end interface

  interface genio_fldlst_reset
    module procedure genio_srclst_reset
    module procedure genio_outlst_reset
  end interface

  interface genio_fldlst_addsrc
    module procedure genio_outlst_addsrc
  end interface

  interface genio_fldlst_print
    module procedure genio_srclst_print
    module procedure genio_outlst_print
  end interface

  interface genio_fldlst_alldata
    module procedure genio_outlst_alldata
  end interface

  interface genio_fldlst_alarm
    module procedure genio_outlst_alarm
  end interface

  interface genio_fldlst_write
    module procedure genio_srclst_write
    module procedure genio_outlst_write
  end interface

  contains

  !-----------------------------------------------------------------------------
  ! Generic IO Field Subroutines
  !-----------------------------------------------------------------------------

  function genio_fldcfg_enabled(fsetcfg, fldnm, rc)
    ! return value
    logical :: genio_fldcfg_enabled
    ! arguments
    type(ESMF_HConfig), intent(in) :: fsetcfg
    character(*), intent(in)       :: fldnm
    integer, intent(out)           :: rc
    ! local variables
    integer                   :: fsetCount
    character(:), allocatable :: fname
    integer                   :: i

    rc = ESMF_SUCCESS
    genio_fldcfg_enabled = .false.

    fsetCount = ESMF_HConfigGetSize(fsetcfg, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    do i=1, fsetCount
      fname = ESMF_HConfigAsString(fsetcfg, index=i, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      if (fname .eq. fldnm) then
        genio_fldcfg_enabled = .true.
        deallocate(fname)
        return
      endif
      deallocate(fname)
    enddo
  endfunction genio_fldcfg_enabled

  !-----------------------------------------------------------------------------

  function genio_srclst_create(state, outpcfg, dflts, name, rc)
    ! return value
    type(genio_srclst) :: genio_srclst_create
    ! arguments
    type(ESMF_State), intent(inout) :: state
    type(ESMF_HConfig), intent(in)  :: outpcfg
    type(genio_dflts), intent(in)   :: dflts
    character(*), intent(in)        :: name
    integer, intent(out)            :: rc
    ! local variables
    logical                                :: fsetcfg_p
    type(ESMF_HConfig)                     :: fsetcfg
    integer                                :: stat
    integer                                :: itemCount
    integer                                :: fieldCount
    integer                                :: fsetCount
    character(len=64),allocatable          :: itemNameList(:)
    type(ESMF_StateItem_Flag), allocatable :: itemTypeList(:)
    character(:), allocatable              :: fname
    type(ESMF_StateItem_Flag)              :: itemType
    integer                                :: i
    integer                                :: n
    logical                                :: enabled

    rc = ESMF_SUCCESS

    genio_srclst_create%name = name
    genio_srclst_create%slice = 0

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

    genio_srclst_create%write = genio_hconfig2logical(outpcfg, &
      "writeSrc", defaultValue=.false., rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    genio_srclst_create%diagn = genio_hconfig2logical(outpcfg, &
      "diagnosticSrc", defaultValue=.false., rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    fsetcfg_p = ESMF_HConfigIsDefined(outpcfg, keyString="fieldset", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    fieldCount = 0
    if (fsetcfg_p) then
      ! count fields in fieldset
      fsetcfg = ESMF_HConfigCreateAt(outpcfg, keyString="fieldset", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      fsetCount = ESMF_HConfigGetSize(fsetcfg, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      do i=1, fsetCount
        fname = ESMF_HConfigAsString(fsetcfg, index=i, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call ESMF_StateGet(state, itemName=fname, itemType=itemType, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        if (itemType .ne. ESMF_STATEITEM_FIELD) then
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
            msg="GENIO: required field is not available - "//trim(fname), &
            line=__LINE__, file=__FILE__, rcToReturn=rc)
          deallocate(fname)
          return
        else
          fieldCount = fieldCount + 1
        endif
        deallocate(fname)
      enddo
    else
      ! count all fields in state
      do i=1, itemCount
        if (itemTypeList(i) == ESMF_STATEITEM_FIELD) then
          fieldCount = fieldCount + 1
        endif
      enddo
    endif

    ! create srclst from enabled fields
    allocate(genio_srclst_create%fld(fieldCount), stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg='GENIO: Memory allocation failed.', &
      line=__LINE__, file=__FILE__, rcToReturn=rc)) return
    n = 0
    do i=1, itemCount
      if (itemTypeList(i) == ESMF_STATEITEM_FIELD) then
        if (fsetcfg_p) then
          enabled = genio_fldcfg_enabled(fsetcfg, itemNameList(i), rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
        else
          enabled = .true.
        endif
        if (enabled) then
          n = n + 1
          genio_srclst_create%fld(n)%name = itemNameList(i)
          call ESMF_StateGet(state, itemName=itemNameList(i), &
            field=genio_srclst_create%fld(n)%srcf, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
        else
          call ESMF_StateRemove(state, itemNameList=(/itemNameList(i)/), &
            rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
        endif ! enabled
      endif ! ESMF_STATEITEM_FIELD
    enddo

    if (fsetcfg_p) then
      call ESMF_HConfigDestroy(fsetcfg, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    endif

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
      call ESMF_FieldDestroy(srclst%fld(i)%srcf, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    enddo

    deallocate(srclst%fld, stat=stat)
    if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
      msg='GENIO: Memory deallocation failed.', &
      line=__LINE__, file=__FILE__, rcToReturn=rc)) return

    srclst%name = "uninitialized"
    srclst%slice = 0

  endsubroutine genio_srclst_destroy

  !-----------------------------------------------------------------------------

  subroutine genio_srclst_realize(state, srclst, rc)
    ! arguments
    type(ESMF_State), intent(inout)   :: state
    type(genio_srclst), intent(inout) :: srclst
    integer, intent(out)              :: rc
    ! local variables
    integer                  :: stat
    integer                  :: i
    integer                  :: dimcount
    type(ESMF_TypeKind_Flag) :: typekind
    real(ESMF_KIND_R8)       :: r8


    rc = ESMF_SUCCESS

    do i=1, size(srclst%fld)
      call NUOPC_Realize(state, fieldName=srclst%fld(i)%name, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call ESMF_FieldGet(srclst%fld(i)%srcf, &
        dimcount=dimcount, typekind=typekind, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      if (typekind .eq. ESMF_TYPEKIND_I4) then
        if (dimcount .eq. 1) then
          srclst%fld(i)%ftyp = GENIO_I4D1
        elseif (dimcount .eq. 2) then
          srclst%fld(i)%ftyp = GENIO_I4D2
        elseif (dimcount .eq. 3) then
          srclst%fld(i)%ftyp = GENIO_I4D3
        else
          call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
            msg="GENIO: field dimcount - "//trim(srclst%fld(i)%name), &
            line=__LINE__, file=__FILE__, rcToReturn=rc)
          return
        endif
        srclst%fld(i)%eval = 0
      elseif (typekind .eq. ESMF_TYPEKIND_R8) then
        if (dimcount .eq. 1) then
          srclst%fld(i)%ftyp = GENIO_R8D1
        elseif (dimcount .eq. 2) then
          srclst%fld(i)%ftyp = GENIO_R8D2
        elseif (dimcount .eq. 3) then
          srclst%fld(i)%ftyp = GENIO_R8D3
        else
          call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
            msg="GENIO: field dimcount - "//trim(srclst%fld(i)%name), &
            line=__LINE__, file=__FILE__, rcToReturn=rc)
          return
        endif
        srclst%fld(i)%eval = epsilon(r8)
      else
        call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
          msg="GENIO: field typekind - "//trim(srclst%fld(i)%name), &
          line=__LINE__, file=__FILE__, rcToReturn=rc)
        return
      endif
      call ESMF_FieldFill(srclst%fld(i)%srcf, dataFillScheme="const", &
        const1=0.0_ESMF_KIND_R8, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    enddo ! srclst

  endsubroutine genio_srclst_realize

  !-----------------------------------------------------------------------------

  subroutine genio_srclst_zero(srclst, rc)
    ! arguments
    type(genio_srclst), intent(inout) :: srclst
    integer, intent(out)              :: rc
    ! local variables
    integer :: i

    rc = ESMF_SUCCESS

    do i=1, size(srclst%fld)
      call ESMF_FieldFill(srclst%fld(i)%srcf, dataFillScheme="const", &
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
      call ESMF_FieldFill(srclst%fld(i)%srcf, dataFillScheme="const", &
        const1=srclst%fld(i)%dflt, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    enddo ! srclst

  endsubroutine genio_srclst_reset

  !-----------------------------------------------------------------------------

  subroutine genio_srclst_print(srclst, myid, outid, vm, cname, label, rc)
    ! arguments
    type(genio_srclst), intent(inout) :: srclst
    integer, intent(in)               :: myid
    integer, intent(in)               :: outid
    type(ESMF_VM), intent(in)         :: vm
    character(*), intent(in)          :: cname
    character(*), intent(in)          :: label
    integer, intent(out)              :: rc
    ! local variables
    integer            :: i

    rc = ESMF_SUCCESS

    do i = 1, size(srclst%fld)
      call genio_srcfld_stats(srclst%fld(i), myid, outid, vm, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    enddo

    if (myid .eq. outid) then
      write(*,'(A,A,A,A,A)') trim(cname)//": ### ", &
        trim(label)," - ",trim(srclst%name)," ###"
      write(*,'(A,5X,A9,3(X,A9))') trim(cname)//":", &
        "COUNT", "MEAN", "MIN", "MAX"
      do i = 1, size(srclst%fld)
        write(*,'(A,3X,A,A1,A)') trim(cname)//":", &
          trim(srclst%fld(i)%name), "@", srclst%fld(i)%tstr
        write(*,'(A,5X,I9,3(X,E9.2))') trim(cname)//":", &
          int(srclst%fld(i)%gsum(2)), srclst%fld(i)%gavg, &
          srclst%fld(i)%gmin(1), srclst%fld(i)%gmax(1)
      enddo
      write(*,'(A)') trim(cname)//":"
    endif

  endsubroutine genio_srclst_print

  !-----------------------------------------------------------------------------

  subroutine genio_srclst_write(srclst, rc)
    ! arguments
    type(genio_srclst), intent(inout) :: srclst
    integer, intent(out)              :: rc
    ! local variables
    type(ESMF_FieldBundle)      :: fb
    integer                     :: i

    rc = ESMF_SUCCESS

    fb = ESMF_FieldBundleCreate(rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    do i = 1, size(srclst%fld)
      call ESMF_FieldBundleAdd(fb, fieldList=(/srclst%fld(i)%srcf/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    enddo
    srclst%slice = srclst%slice + 1
    call ESMF_FieldBundleWrite(fb, &
      fileName=trim(srclst%name)//".nc", &
      overwrite=.true., timeslice=srclst%slice, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call ESMF_FieldBundleDestroy(fb, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

  endsubroutine genio_srclst_write

  !-----------------------------------------------------------------------------

  subroutine genio_srcfld_stats(fld, myid, outid, vm, rc)
    ! arguments
    type(genio_fld), intent(inout) :: fld
    integer, intent(in)            :: myid
    integer, intent(in)            :: outid
    type(ESMF_VM), intent(in)      :: vm
    integer, intent(out)           :: rc
    ! local variables
    type(ESMF_Time)             :: fieldtime
    integer                     :: localDeCount
    integer                     :: i
    real(ESMF_KIND_R8), pointer :: srcR8D1(:)     => null()
    real(ESMF_KIND_R8), pointer :: srcR8D2(:,:)   => null()
    real(ESMF_KIND_R8), pointer :: srcR8D3(:,:,:) => null()

    rc = ESMF_SUCCESS

    call NUOPC_GetTimestamp(fld%srcf, time=fieldtime, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call ESMF_TimeGet(fieldtime, timeString=fld%tstr, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call ESMF_FieldGet(fld%srcf, localDeCount=localDeCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    if (localDeCount .lt. 1) then
      fld%lsum(1)=0.0_ESMF_KIND_R8
      fld%lsum(2)=0
      fld%lmin(1)=0.0_ESMF_KIND_R8
      fld%lmax(1)=0.0_ESMF_KIND_R8
    else
      if(fld%ftyp .eq. GENIO_R8D1) then
        call ESMF_FieldGet(fld%srcf, localDE=0, farrayPtr=srcR8D1, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        fld%lsum(1)=sum(srcR8D1,abs(srcR8D1-GENIO_FILV)>fld%eval)
        fld%lsum(2)=count(abs(srcR8D1-GENIO_FILV)>fld%eval)
        fld%lmin(1)=minval(srcR8D1,abs(srcR8D1-GENIO_FILV)>fld%eval)
        fld%lmax(1)=maxval(srcR8D1,abs(srcR8D1-GENIO_FILV)>fld%eval)
        do i=1, localDeCount-1
          call ESMF_FieldGet(fld%srcf, localDE=i, farrayPtr=srcR8D1, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
          fld%lsum(1)=fld%lsum(1) + &
                      sum(srcR8D1,abs(srcR8D1-GENIO_FILV)>fld%eval)
          fld%lsum(2)=fld%lsum(2) + &
                      count(abs(srcR8D1-GENIO_FILV)>fld%eval)
          fld%lmin(1)=min(fld%lmin(1), &
                          minval(srcR8D1,abs(srcR8D1-GENIO_FILV)>fld%eval))
          fld%lmax(1)=max(fld%lmax(1), &
                          maxval(srcR8D1,abs(srcR8D1-GENIO_FILV)>fld%eval))
        enddo
      elseif(fld%ftyp .eq. GENIO_R8D2) then
        call ESMF_FieldGet(fld%srcf, localDE=0, farrayPtr=srcR8D2, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        fld%lsum(1)=sum(srcR8D2,abs(srcR8D2-GENIO_FILV)>fld%eval)
        fld%lsum(2)=count(abs(srcR8D2-GENIO_FILV)>fld%eval)
        fld%lmin(1)=minval(srcR8D2,abs(srcR8D2-GENIO_FILV)>fld%eval)
        fld%lmax(1)=maxval(srcR8D2,abs(srcR8D2-GENIO_FILV)>fld%eval)
        do i=1, localDeCount-1
          call ESMF_FieldGet(fld%srcf, localDE=i, farrayPtr=srcR8D2, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
          fld%lsum(1)=fld%lsum(1) + &
                      sum(srcR8D2,abs(srcR8D2-GENIO_FILV)>fld%eval)
          fld%lsum(2)=fld%lsum(2) + &
                      count(abs(srcR8D2-GENIO_FILV)>fld%eval)
          fld%lmin(1)=min(fld%lmin(1), &
                          minval(srcR8D2,abs(srcR8D2-GENIO_FILV)>fld%eval))
          fld%lmax(1)=max(fld%lmax(1), &
                          maxval(srcR8D2,abs(srcR8D2-GENIO_FILV)>fld%eval))
        enddo
      elseif(fld%ftyp .eq. GENIO_R8D3) then
        call ESMF_FieldGet(fld%srcf, localDE=0, farrayPtr=srcR8D3, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        fld%lsum(1)=sum(srcR8D3,abs(srcR8D3-GENIO_FILV)>fld%eval)
        fld%lsum(2)=count(abs(srcR8D3-GENIO_FILV)>fld%eval)
        fld%lmin(1)=minval(srcR8D3,abs(srcR8D3-GENIO_FILV)>fld%eval)
        fld%lmax(1)=maxval(srcR8D3,abs(srcR8D3-GENIO_FILV)>fld%eval)
        do i=1, localDeCount-1
          call ESMF_FieldGet(fld%srcf, localDE=i, farrayPtr=srcR8D3, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
          fld%lsum(1)=fld%lsum(1) + &
                      sum(srcR8D3,abs(srcR8D3-GENIO_FILV)>fld%eval)
          fld%lsum(2)=fld%lsum(2) + &
                      count(abs(srcR8D3-GENIO_FILV)>fld%eval)
          fld%lmin(1)=min(fld%lmin(1), &
                          minval(srcR8D3,abs(srcR8D3-GENIO_FILV)>fld%eval))
          fld%lmax(1)=max(fld%lmax(1), &
                          maxval(srcR8D3,abs(srcR8D3-GENIO_FILV)>fld%eval))
        enddo
      else
        call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
          msg="GENIO: field typekind-dim "//trim(fld%name), &
          line=__LINE__, file=__FILE__, rcToReturn=rc)
        return
      endif
    endif ! localDeCount

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

  endsubroutine genio_srcfld_stats

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------

  function genio_outlst_create(srclst, outpcfg, fldscfg, dflts, name, rc)
    ! return value
    type(genio_outlst) :: genio_outlst_create
    ! arguments
    type(genio_srclst), intent(in)      :: srclst
    type(ESMF_HConfig), intent(in)      :: outpcfg
    type(ESMF_HConfig), intent(in)      :: fldscfg
    type(genio_dflts), intent(in)       :: dflts
    character(*), intent(in)            :: name
    integer, intent(out)                :: rc
    ! local variables
    type(ESMF_TimeInterval)   :: zerotimeint
    integer                   :: stat
    integer                   :: i
    type(ESMF_HConfig)        :: fldsetcfg
    character(:), allocatable :: fname
    type(ESMF_Geom)           :: geom
    real(ESMF_KIND_R8)        :: r8
    logical                    :: crrfcfg_p
    logical                    :: otypcfg_p
    type(ESMF_HConfig)         :: crrfcfg

    rc = ESMF_SUCCESS

    genio_outlst_create%name = name
    genio_outlst_create%slice = 0
    genio_outlst_create%wfinl = genio_hconfig2logical(outpcfg, &
      "writeFinal", defaultValue=.false., rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    ! diagnostics
    genio_outlst_create%diagn = genio_hconfig2logical(outpcfg, &
      "diagnosticOut", defaultValue=.false., rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    ! default time interval
    call ESMF_TimeIntervalSet(zerotimeint, s=0, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    genio_outlst_create%zrofrq = zerotimeint
    genio_outlst_create%accwnd = zerotimeint
    genio_outlst_create%outfrq = genio_hconfig2timeint(outpcfg, &
      "frequency", defaultValue=zerotimeint, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    genio_outlst_create%eval = epsilon(r8)

    allocate(genio_outlst_create%fld(size(srclst%fld)), stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg='GENIO: Memory allocation failed.', &
      line=__LINE__, file=__FILE__, rcToReturn=rc)) return
    do i = 1, size(genio_outlst_create%fld)
      genio_outlst_create%fld(i)%name = srclst%fld(i)%name
      genio_outlst_create%fld(i)%ftyp = srclst%fld(i)%ftyp
      genio_outlst_create%fld(i)%srcf = srclst%fld(i)%srcf
      genio_outlst_create%fld(i)%eval = srclst%fld(i)%eval
      genio_outlst_create%fld(i)%prvf = ESMF_FieldCreate(srclst%fld(i)%srcf, &
        datacopyflag=ESMF_DATACOPY_VALUE, name=srclst%fld(i)%name, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      genio_outlst_create%fld(i)%outf = ESMF_FieldCreate(srclst%fld(i)%srcf, &
        datacopyflag=ESMF_DATACOPY_VALUE, name=srclst%fld(i)%name, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call ESMF_FieldGet(srclst%fld(i)%srcf, geom=geom, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      genio_outlst_create%fld(i)%sumf = ESMF_FieldCreate(geom, &
        typekind=ESMF_TYPEKIND_R8, name=srclst%fld(i)%name, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call ESMF_FieldFill(genio_outlst_create%fld(i)%sumf, &
        dataFillScheme="const", const1=0.0_ESMF_KIND_R8, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      ! field configuration
      crrfcfg_p = ESMF_HConfigIsDefined(fldscfg, &
        keyString=srclst%fld(i)%name, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      if (crrfcfg_p) then
        crrfcfg = ESMF_HConfigCreateAt(fldscfg, &
          keyString=srclst%fld(i)%name, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        genio_outlst_create%fld(i)%otyp = genio_hconfig2otyp(crrfcfg, &
          "outputType", defaultValue=dflts%otyp, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call ESMF_HConfigDestroy(crrfcfg, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
      else
        genio_outlst_create%fld(i)%otyp = dflts%otyp
      endif
      if (genio_outlst_create%fld(i)%otyp .eq. GENIO_INST) then
        genio_outlst_create%fld(i)%ostr = "(inst)"
      elseif (genio_outlst_create%fld(i)%otyp .eq. GENIO_ACCM) then
        genio_outlst_create%fld(i)%ostr = "(accm)"
      elseif (genio_outlst_create%fld(i)%otyp .eq. GENIO_MEAN) then
        genio_outlst_create%fld(i)%ostr = "(mean)"
      else
        genio_outlst_create%fld(i)%ostr = "(----)"
      endif
      genio_outlst_create%fld(i)%timw = 0
    enddo

  endfunction genio_outlst_create

  !-----------------------------------------------------------------------------

  subroutine genio_outlst_destroy(outlst, rc)
    ! arguments
    type(genio_outlst), intent(inout) :: outlst
    integer, intent(out)              :: rc
    ! local variables
    integer                   :: stat
    integer                   :: i

    rc = ESMF_SUCCESS

    do i = 1, size(outlst%fld)
      call ESMF_FieldDestroy(outlst%fld(i)%prvf, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call ESMF_FieldDestroy(outlst%fld(i)%outf, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call ESMF_FieldDestroy(outlst%fld(i)%sumf, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    enddo
    deallocate(outlst%fld, stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg='GENIO: Memory allocation failed.', &
      line=__LINE__, file=__FILE__, rcToReturn=rc)) return

    outlst%outfrq = outlst%zrofrq
    outlst%name = "uninitialized"
    outlst%slice = 0
  endsubroutine genio_outlst_destroy

  !-----------------------------------------------------------------------------

  subroutine genio_outlst_reset(outlst, clock, rc)
    ! arguments
    type(genio_outlst), intent(inout) :: outlst
    type(ESMF_Clock), intent(in)      :: clock
    integer, intent(out)              :: rc
    ! local variables
    integer                 :: i
    type(ESMF_Time)         :: currtime
    type(ESMF_Time)         :: stoptime
    logical                 :: iscurr
    logical                 :: isstop
    integer                 :: localDeCount
    integer                 :: j

    rc = ESMF_SUCCESS

    call ESMF_ClockGet(clock, currtime=currtime, stoptime=stoptime, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    ! reset field accumulation
    do i=1, size(outlst%fld)
      iscurr = NUOPC_IsAtTime(outlst%fld(i)%srcf, currtime, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      isstop = NUOPC_IsAtTime(outlst%fld(i)%srcf, stoptime, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      if (iscurr .or. isstop) then
        outlst%fld(i)%strt = .true.
      else
        outlst%fld(i)%strt = .false.
      endif
      call ESMF_FieldCopy(outlst%fld(i)%prvf, outlst%fld(i)%srcf, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call ESMF_FieldFill(outlst%fld(i)%outf, dataFillScheme="const", &
        const1=0.0_ESMF_KIND_R8, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call ESMF_FieldFill(outlst%fld(i)%sumf, dataFillScheme="const", &
        const1=0.0_ESMF_KIND_R8, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      outlst%fld(i)%timw = 0.0_ESMF_KIND_R8
      iscurr = NUOPC_IsAtTime(outlst%fld(i)%srcf, currtime, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      isstop = NUOPC_IsAtTime(outlst%fld(i)%srcf, stoptime, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      if (iscurr .or. isstop) then
        outlst%fld(i)%strt = .true.
        call NUOPC_GetTimestamp(outlst%fld(i)%srcf, &
          time=outlst%fld(i)%prvt, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
      else
        outlst%fld(i)%strt = .false.
        outlst%fld(i)%prvt = currtime
      endif
    enddo
    ! reset list clocks
    outlst%accwnd = outlst%zrofrq
    outlst%prvtme = currtime

  endsubroutine genio_outlst_reset

  !-----------------------------------------------------------------------------

  subroutine genio_outlst_addsrc(outlst, clock, rc)
    ! arguments
    type(genio_outlst), intent(inout) :: outlst
    type(ESMF_Clock), intent(in)      :: clock
    integer, intent(out)              :: rc
    ! local variables
    integer                     :: i
    logical                     :: tsUpdate
    integer                     :: atype
    type(ESMF_TimeInterval)     :: timestep
    type(ESMF_Time)             :: fieldtime
    type(ESMF_TimeInterval)     :: modelts
    type(ESMF_TimeInterval)     :: fieldts
    real(ESMF_KIND_R8)          :: ts_s

    rc = ESMF_SUCCESS

    call ESMF_ClockGet(clock, timestep=timestep, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    do i=1, size(outlst%fld)
      call NUOPC_GetTimestamp(outlst%fld(i)%srcf, &
        time=fieldtime, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      fieldts = fieldtime - outlst%fld(i)%prvt
      call ESMF_TimeIntervalGet(fieldts, s_r8=ts_s, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      ts_s = abs(ts_s)
      if (ts_s .le. outlst%eval) then
        continue ! no timestamp update
      elseif (outlst%fld(i)%otyp .eq. GENIO_INST) then
        continue ! instantaneous output only
      elseif (outlst%fld(i)%otyp .eq. GENIO_ACCM) then
        call genio_outfld_accmsumf(outlst%fld(i), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
      elseif (outlst%fld(i)%otyp .eq. GENIO_MEAN) then
        if (outlst%fld(i)%strt) then
          call genio_outfld_meansumf(outlst%fld(i), ts_s, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
        endif
      else
        call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
          msg="GENIO: field output type "//trim(outlst%fld(i)%name), &
          line=__LINE__, file=__FILE__, rcToReturn=rc)
        return
      endif
      call ESMF_FieldCopy(outlst%fld(i)%prvf, outlst%fld(i)%srcf, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      outlst%fld(i)%strt = .true.
      outlst%fld(i)%prvt = fieldtime
      outlst%fld(i)%timw = outlst%fld(i)%timw + ts_s
    enddo
    ! set output window
    outlst%accwnd = outlst%accwnd + timestep

  endsubroutine genio_outlst_addsrc

  !-----------------------------------------------------------------------------

  function genio_outlst_alldata(outlst, rc)
    ! return value
    logical :: genio_outlst_alldata
    ! arguments
    type(genio_outlst), intent(in) :: outlst
    integer, intent(out)           :: rc
    ! local variables
    integer :: i

    rc = ESMF_SUCCESS

    genio_outlst_alldata = .true.
    do i=1, size(outlst%fld)
      if (.not. outlst%fld(i)%strt) then
        genio_outlst_alldata = .false.
        return
      endif
    enddo

  endfunction genio_outlst_alldata

  !-----------------------------------------------------------------------------

  function genio_outlst_alarm(outlst, rc)
    ! return value
    logical :: genio_outlst_alarm
    ! arguments
    type(genio_outlst), intent(in) :: outlst
    integer, intent(out)           :: rc

    rc = ESMF_SUCCESS

    if (outlst%outfrq .eq. outlst%zrofrq) then
      genio_outlst_alarm = .true.
    elseif (outlst%accwnd .ge. outlst%outfrq) then
      genio_outlst_alarm = .true.
    else
      genio_outlst_alarm = .false.
    endif

  endfunction genio_outlst_alarm

  !-----------------------------------------------------------------------------

  subroutine genio_outlst_print(outlst, myid, outid, vm, cname, label, &
  rc)
    ! arguments
    type(genio_outlst), intent(inout) :: outlst
    integer, intent(in)               :: myid
    integer, intent(in)               :: outid
    type(ESMF_VM), intent(in)         :: vm
    character(*), intent(in)          :: cname
    character(*), intent(in)          :: label
    integer, intent(out)              :: rc
    ! local variables
    integer            :: i

    rc = ESMF_SUCCESS
    do i = 1, size(outlst%fld)
      if (outlst%fld(i)%otyp .eq. GENIO_INST) then
        call ESMF_FieldCopy(outlst%fld(i)%outf, outlst%fld(i)%prvf, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
      elseif (outlst%fld(i)%otyp .eq. GENIO_ACCM) then
        call ESMF_FieldCopy(outlst%fld(i)%outf, outlst%fld(i)%sumf, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
      elseif (outlst%fld(i)%otyp .eq. GENIO_MEAN) then
        if (outlst%fld(i)%timw .gt. outlst%eval) then
          call genio_outfld_meanoutf(outlst%fld(i), rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
        else
          call ESMF_FieldCopy(outlst%fld(i)%outf, outlst%fld(i)%prvf, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
        endif
      endif
      call genio_outfld_stats(outlst%fld(i), myid, outid, vm, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    enddo

    if (myid .eq. outid) then
      write(*,'(A,A,A,A,A)') trim(cname)//": ### ", &
        trim(label)," - ",trim(outlst%name)," ###"
      write(*,'(A,5X,A11,X,A9,3(X,A9))') trim(cname)//":", &
        "TIMEW", "COUNT", "MEAN", "MIN", "MAX"
      do i = 1, size(outlst%fld)
        write(*,'(A,3X,A,A1,A,X,A)') trim(cname)//":", &
          trim(outlst%fld(i)%name), "@", outlst%fld(i)%tstr, outlst%fld(i)%ostr
        write(*,'(A,5X,F11.1,X,I9,3(X,E9.2))') trim(cname)//":", &
          outlst%fld(i)%timw, &
          int(outlst%fld(i)%gsum(2)), outlst%fld(i)%gavg, &
          outlst%fld(i)%gmin(1), outlst%fld(i)%gmax(1)
      enddo
      write(*,'(A)') trim(cname)//":"
    endif

  endsubroutine genio_outlst_print

  !-----------------------------------------------------------------------------

  subroutine genio_outlst_write(outlst, rc)
    ! arguments
    type(genio_outlst), intent(inout) :: outlst
    integer, intent(out)              :: rc
    ! local variables
    type(ESMF_FieldBundle)      :: fb
    integer                     :: i

    rc = ESMF_SUCCESS

    fb = ESMF_FieldBundleCreate(rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    do i = 1, size(outlst%fld)
      if (outlst%fld(i)%otyp .eq. GENIO_INST) then
        call ESMF_FieldCopy(outlst%fld(i)%outf, outlst%fld(i)%prvf, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
      elseif (outlst%fld(i)%otyp .eq. GENIO_ACCM) then
        call ESMF_FieldCopy(outlst%fld(i)%outf, outlst%fld(i)%sumf, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
      elseif (outlst%fld(i)%otyp .eq. GENIO_MEAN) then
        if (outlst%fld(i)%timw .gt. outlst%eval) then
          call genio_outfld_meanoutf(outlst%fld(i), rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
        else
          call ESMF_FieldCopy(outlst%fld(i)%outf, outlst%fld(i)%prvf, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
        endif
      endif
      call ESMF_FieldBundleAdd(fb, fieldList=(/outlst%fld(i)%outf/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    enddo
    outlst%slice = outlst%slice + 1
    call ESMF_FieldBundleWrite(fb, &
      fileName=trim(outlst%name)//".nc", &
      overwrite=.true., timeslice=outlst%slice, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call ESMF_FieldBundleDestroy(fb, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

  endsubroutine genio_outlst_write

  !---------------------------------------------------------------------------

  subroutine genio_outfld_stats(fld, myid, outid, vm, rc)
    ! arguments
    type(genio_fld), intent(inout) :: fld
    integer, intent(in)            :: myid
    integer, intent(in)            :: outid
    type(ESMF_VM), intent(in)      :: vm
    integer, intent(out)           :: rc
    ! local variables
    integer                     :: localDeCount
    integer                     :: i
    real(ESMF_KIND_R8), pointer :: outR8D1(:)     => null()
    real(ESMF_KIND_R8), pointer :: outR8D2(:,:)   => null()
    real(ESMF_KIND_R8), pointer :: outR8D3(:,:,:) => null()

    rc = ESMF_SUCCESS

    call ESMF_TimeGet(fld%prvt, timeString=fld%tstr, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call ESMF_FieldGet(fld%outf, localDeCount=localDeCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    if (localDeCount .lt. 1) then
      fld%lsum(1)=0.0_ESMF_KIND_R8
      fld%lsum(2)=0
      fld%lmin(1)=0.0_ESMF_KIND_R8
      fld%lmax(1)=0.0_ESMF_KIND_R8
    else
      if(fld%ftyp .eq. GENIO_R8D1) then
        call ESMF_FieldGet(fld%outf, localDE=0, farrayPtr=outR8D1, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        fld%lsum(1)=sum(outR8D1,abs(outR8D1-GENIO_FILV)>fld%eval)
        fld%lsum(2)=count(abs(outR8D1-GENIO_FILV)>fld%eval)
        fld%lmin(1)=minval(outR8D1,abs(outR8D1-GENIO_FILV)>fld%eval)
        fld%lmax(1)=maxval(outR8D1,abs(outR8D1-GENIO_FILV)>fld%eval)
        do i=1, localDeCount-1
          call ESMF_FieldGet(fld%outf, localDE=i, farrayPtr=outR8D1, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
          fld%lsum(1)=fld%lsum(1) + &
                      sum(outR8D1,abs(outR8D1-GENIO_FILV)>fld%eval)
          fld%lsum(2)=fld%lsum(2) + &
                      count(abs(outR8D1-GENIO_FILV)>fld%eval)
          fld%lmin(1)=min(fld%lmin(1), &
                          minval(outR8D1,abs(outR8D1-GENIO_FILV)>fld%eval))
          fld%lmax(1)=max(fld%lmax(1), &
                          maxval(outR8D1,abs(outR8D1-GENIO_FILV)>fld%eval))
        enddo
      elseif(fld%ftyp .eq. GENIO_R8D2) then
        call ESMF_FieldGet(fld%outf, localDE=0, farrayPtr=outR8D2, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        fld%lsum(1)=sum(outR8D2,abs(outR8D2-GENIO_FILV)>fld%eval)
        fld%lsum(2)=count(abs(outR8D2-GENIO_FILV)>fld%eval)
        fld%lmin(1)=minval(outR8D2,abs(outR8D2-GENIO_FILV)>fld%eval)
        fld%lmax(1)=maxval(outR8D2,abs(outR8D2-GENIO_FILV)>fld%eval)
        do i=1, localDeCount-1
          call ESMF_FieldGet(fld%outf, localDE=i, farrayPtr=outR8D2, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
          fld%lsum(1)=fld%lsum(1) + &
                      sum(outR8D2,abs(outR8D2-GENIO_FILV)>fld%eval)
          fld%lsum(2)=fld%lsum(2) + &
                      count(abs(outR8D2-GENIO_FILV)>fld%eval)
          fld%lmin(1)=min(fld%lmin(1), &
                          minval(outR8D2,abs(outR8D2-GENIO_FILV)>fld%eval))
          fld%lmax(1)=max(fld%lmax(1), &
                          maxval(outR8D2,abs(outR8D2-GENIO_FILV)>fld%eval))
        enddo
      elseif(fld%ftyp .eq. GENIO_R8D3) then
        call ESMF_FieldGet(fld%outf, localDE=0, farrayPtr=outR8D3, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        fld%lsum(1)=sum(outR8D3,abs(outR8D3-GENIO_FILV)>fld%eval)
        fld%lsum(2)=count(abs(outR8D3-GENIO_FILV)>fld%eval)
        fld%lmin(1)=minval(outR8D3,abs(outR8D3-GENIO_FILV)>fld%eval)
        fld%lmax(1)=maxval(outR8D3,abs(outR8D3-GENIO_FILV)>fld%eval)
        do i=1, localDeCount-1
          call ESMF_FieldGet(fld%outf, localDE=i, farrayPtr=outR8D3, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
          fld%lsum(1)=fld%lsum(1) + &
                      sum(outR8D3,abs(outR8D3-GENIO_FILV)>fld%eval)
          fld%lsum(2)=fld%lsum(2) + &
                      count(abs(outR8D3-GENIO_FILV)>fld%eval)
          fld%lmin(1)=min(fld%lmin(1), &
                          minval(outR8D3,abs(outR8D3-GENIO_FILV)>fld%eval))
          fld%lmax(1)=max(fld%lmax(1), &
                          maxval(outR8D3,abs(outR8D3-GENIO_FILV)>fld%eval))
        enddo
      else
        call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
          msg="GENIO: field typekind-dim "//trim(fld%name), &
          line=__LINE__, file=__FILE__, rcToReturn=rc)
        return
      endif
    endif ! localDeCount

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

  endsubroutine genio_outfld_stats

  !---------------------------------------------------------------------------

  subroutine genio_outfld_accmsumf(fld, rc)
    ! arguments
    type(genio_fld), intent(inout)      :: fld
    integer, intent(out)                :: rc
    ! local variables
    real(ESMF_KIND_R8), pointer :: srcR8D1(:)
    real(ESMF_KIND_R8), pointer :: srcR8D2(:,:)
    real(ESMF_KIND_R8), pointer :: srcR8D3(:,:,:)
    real(ESMF_KIND_R8), pointer :: sumR8D1(:)
    real(ESMF_KIND_R8), pointer :: sumR8D2(:,:)
    real(ESMF_KIND_R8), pointer :: sumR8D3(:,:,:)
    integer                     :: localDeCount
    integer                     :: i

    call ESMF_FieldGet(fld%srcf, localDeCount=localDeCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    do i=0, localDeCount-1
      if (fld%ftyp .eq. GENIO_R8D1) then
        call ESMF_FieldGet(fld%srcf, localDe=i, farrayPtr=srcR8D1, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call ESMF_FieldGet(fld%sumf, localDe=i, farrayPtr=sumR8D1, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        sumR8D1 = sumR8D1 + srcR8D1
      elseif (fld%ftyp .eq. GENIO_R8D2) then
        call ESMF_FieldGet(fld%srcf, localDe=i, farrayPtr=srcR8D2, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call ESMF_FieldGet(fld%sumf, localDe=i, farrayPtr=sumR8D2, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        sumR8D2 = sumR8D2 + srcR8D2
      elseif (fld%ftyp .eq. GENIO_R8D3) then
        call ESMF_FieldGet(fld%srcf, localDe=i, farrayPtr=srcR8D3, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call ESMF_FieldGet(fld%sumf, localDe=i, farrayPtr=sumR8D3, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        sumR8D3 = sumR8D3 + srcR8D3
      else
        call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
          msg="GENIO: field output type "//trim(fld%name), &
          line=__LINE__, file=__FILE__, rcToReturn=rc)
        return
      endif
    enddo
  endsubroutine genio_outfld_accmsumf

  !---------------------------------------------------------------------------

  subroutine genio_outfld_meansumf(fld, ts_s, rc)
    ! arguments
    type(genio_fld), intent(inout) :: fld
    real(ESMF_KIND_R8), intent(in) :: ts_s
    integer, intent(out)           :: rc
    ! local variables
    real(ESMF_KIND_R8), pointer :: srcR8D1(:)
    real(ESMF_KIND_R8), pointer :: srcR8D2(:,:)
    real(ESMF_KIND_R8), pointer :: srcR8D3(:,:,:)
    real(ESMF_KIND_R8), pointer :: prvR8D1(:)
    real(ESMF_KIND_R8), pointer :: prvR8D2(:,:)
    real(ESMF_KIND_R8), pointer :: prvR8D3(:,:,:)
    real(ESMF_KIND_R8), pointer :: sumR8D1(:)
    real(ESMF_KIND_R8), pointer :: sumR8D2(:,:)
    real(ESMF_KIND_R8), pointer :: sumR8D3(:,:,:)
    integer                     :: localDeCount
    integer                     :: i

    call ESMF_FieldGet(fld%srcf, localDeCount=localDeCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    do i=0, localDeCount-1
      if (fld%ftyp .eq. GENIO_R8D1) then
        call ESMF_FieldGet(fld%srcf, localDe=i, farrayPtr=srcR8D1, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call ESMF_FieldGet(fld%prvf, localDe=i, farrayPtr=prvR8D1, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call ESMF_FieldGet(fld%sumf, localDe=i, farrayPtr=sumR8D1, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        sumR8D1 = sumR8D1 + &
                  (real((srcR8D1 + prvR8D1), ESMF_KIND_R8) / 2 * ts_s)
      elseif (fld%ftyp .eq. GENIO_R8D2) then
        call ESMF_FieldGet(fld%srcf, localDe=i, farrayPtr=srcR8D2, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call ESMF_FieldGet(fld%prvf, localDe=i, farrayPtr=prvR8D2, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call ESMF_FieldGet(fld%sumf, localDe=i, farrayPtr=sumR8D2, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        sumR8D2 = sumR8D2 + &
                  (real((srcR8D2 + prvR8D2), ESMF_KIND_R8) / 2 * ts_s)
      elseif (fld%ftyp .eq. GENIO_R8D3) then
        call ESMF_FieldGet(fld%srcf, localDe=i, farrayPtr=srcR8D3, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call ESMF_FieldGet(fld%prvf, localDe=i, farrayPtr=prvR8D3, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call ESMF_FieldGet(fld%sumf, localDe=i, farrayPtr=sumR8D3, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        sumR8D3 = sumR8D3 + &
                  (real((srcR8D3 + prvR8D3), ESMF_KIND_R8) / 2 * ts_s)
      else
        call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
          msg="GENIO: field typekind-dim "//trim(fld%name), &
          line=__LINE__, file=__FILE__, rcToReturn=rc)
        return
      endif
    enddo
  endsubroutine genio_outfld_meansumf

  !---------------------------------------------------------------------------

  subroutine genio_outfld_meanoutf(fld, rc)
    ! arguments
    type(genio_fld), intent(inout) :: fld
    integer, intent(out)           :: rc
    ! local variables
    real(ESMF_KIND_R8), pointer :: outR8D1(:)
    real(ESMF_KIND_R8), pointer :: outR8D2(:,:)
    real(ESMF_KIND_R8), pointer :: outR8D3(:,:,:)
    real(ESMF_KIND_R8), pointer :: sumR8D1(:)
    real(ESMF_KIND_R8), pointer :: sumR8D2(:,:)
    real(ESMF_KIND_R8), pointer :: sumR8D3(:,:,:)
    integer                     :: localDeCount
    integer                     :: i

    call ESMF_FieldGet(fld%sumf, localDeCount=localDeCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    do i=0, localDeCount-1
      if (fld%ftyp .eq. GENIO_R8D1) then
        call ESMF_FieldGet(fld%outf, localDe=i, farrayPtr=outR8D1, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call ESMF_FieldGet(fld%sumf, localDe=i, farrayPtr=sumR8D1, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        outR8D1 = sumR8D1 / fld%timw
      elseif (fld%ftyp .eq. GENIO_R8D2) then
        call ESMF_FieldGet(fld%outf, localDe=i, farrayPtr=outR8D2, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call ESMF_FieldGet(fld%sumf, localDe=i, farrayPtr=sumR8D2, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        outR8D2 = sumR8D2 / fld%timw
      elseif (fld%ftyp .eq. GENIO_R8D3) then
        call ESMF_FieldGet(fld%outf, localDe=i, farrayPtr=outR8D3, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call ESMF_FieldGet(fld%sumf, localDe=i, farrayPtr=sumR8D3, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        outR8D3 = sumR8D3 / fld%timw
      else
        call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
          msg="GENIO: field typekind-dim "//trim(fld%name), &
          line=__LINE__, file=__FILE__, rcToReturn=rc)
        return
      endif
    enddo
  endsubroutine genio_outfld_meanoutf

endmodule genio_mod_field
