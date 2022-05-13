module Ensemble_driver

  !-----------------------------------------------------------------------------
  ! Code that creates the ensemble driver layer above the esm driver instance.
  ! The ensmeble driver is configured to run a single clock cycle in nuopc with time step
  ! length of stop_time - start_time.  It's purpose is to instantiate NINST copies of the
  ! esm driver and its components layed out concurently across mpi tasks.
  !-----------------------------------------------------------------------------

  use shr_kind_mod  , only : cl=>shr_kind_cl, cs=>shr_kind_cs
  use shr_log_mod   , only : shrlogunit=> shr_log_unit
  use shr_file_mod  , only : shr_file_setLogUnit
  use esm_utils_mod , only : mastertask, logunit, chkerr

  implicit none
  private

  public  :: SetServices
  private :: SetModelServices
!  private :: InitializeP1
  private :: InitializeIO

  integer,  allocatable  :: asyncio_petlist(:)
  character(*),parameter :: u_FILE_u = &
       __FILE__

!================================================================================
contains
!================================================================================

  subroutine SetServices(ensemble_driver, rc)

    use NUOPC        , only : NUOPC_CompDerive, NUOPC_CompSpecialize
    use NUOPC_Driver , only : driver_routine_SS             => SetServices
    use NUOPC_Driver , only : ensemble_label_SetModelServices => label_SetModelServices
    use NUOPC_Driver , only : ensemble_label_SetRunSequence => label_SetRunSequence
    use NUOPC_Driver , only : ensemble_label_ModifyCplLists => label_ModifyCplLists
    use ESMF         , only : ESMF_GridComp, ESMF_GridCompSet
    use ESMF         , only : ESMF_Config, ESMF_ConfigCreate, ESMF_ConfigLoadFile
    use ESMF         , only : ESMF_SUCCESS, ESMF_LogWrite, ESMF_LOGMSG_INFO
    use ESMF         , only : ESMF_GridCompSetEntryPoint, ESMF_Method_Initialize
    
    type(ESMF_GridComp)  :: ensemble_driver
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Config) :: config
    character(len=*), parameter :: subname = u_FILE_u//":SetServices)"
    !---------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)

    ! NUOPC_Driver registers the generic methods
    call NUOPC_CompDerive(ensemble_driver, driver_routine_SS, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! attach specializing method(s)
    call NUOPC_CompSpecialize(ensemble_driver, specLabel=ensemble_label_SetModelServices, &
         specRoutine=SetModelServices, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

!    call ESMF_GridCompSetEntryPoint(ensemble_driver, ESMF_METHOD_INITIALIZE, userRoutine=InitializeP1, phase=1, rc=rc)

!    call ESMF_GridCompSetEntryPoint(ensemble_driver, ESMF_METHOD_INITIALIZE, userRoutine=InitializeIO, phase=0, rc=rc)
!    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSpecialize(ensemble_driver, specLabel=ensemble_label_ModifyCplLists, &
!    call NUOPC_CompSpecialize(ensemble_driver, specLabel=ensemble_label_SetRunSequence, &
         specRoutine=InitializeIO, rc=rc)

    ! Create, open and set the config
    config = ESMF_ConfigCreate(rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ConfigLoadFile(config, "nuopc.runconfig", rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_GridCompSet(ensemble_driver, config=config, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

  end subroutine SetServices
#ifdef USETHIS
  !================================================================================
  recursive subroutine InitializeP1(driver, importState, exportState, clock, rc)
    use ESMF
    use NUOPC
    type(ESMF_GridComp)   :: driver
    type(ESMF_State)      :: importState, exportState
    type(ESMF_Clock)      :: clock
    integer, intent(out)  :: rc
    
    ! local variables
    character(*), parameter   :: rName="InitializeP1"
    character(ESMF_MAXSTR)    :: name
    integer                   :: verbosity, profiling
!    type(type_InternalState)  :: is
    logical                   :: isSet
    character(len=*), parameter :: subname = u_FILE_u//":InitializeP1)"
    rc = ESMF_SUCCESS
    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)

    ! query the component for info
    call NUOPC_CompGet(driver, name=name, verbosity=verbosity, &
      profiling=profiling, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=subname)) return  ! bail out

    ! handle profiling
    if (btest(profiling,9)) then
      call ESMF_TraceRegionEnter("Leading Barrier", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=subname, rcToReturn=rc)) &
        return  ! bail out
      call ESMF_VMBarrier(rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=subname, rcToReturn=rc)) &
        return  ! bail out
      call ESMF_TraceRegionExit("Leading Barrier", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=subname, rcToReturn=rc)) &
        return  ! bail out
    endif
    if (btest(profiling,0)) then
      call ESMF_TraceRegionEnter(rName, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=subname, rcToReturn=rc)) &
        return  ! bail out
    endif

    ! intro
    call NUOPC_LogIntro(name, rName, verbosity, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=subname)) return  ! bail out

    ! check if HierarchyProtocol attribute was set
    call NUOPC_CompAttributeGet(driver, name="HierarchyProtocol", isSet=isSet, &
      rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=subname)) return  ! bail out

    if (.not.isSet) then
      ! turn hierarchy support to connect outside NUOPC
      call NUOPC_CompAttributeSet(driver, &
        name="HierarchyProtocol", value="ConnectProvidedFields", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=subname)) return  ! bail out
    endif
    call ESMF_LogWrite(trim(subname)//"v02p01: called", ESMF_LOGMSG_INFO)

    ! call the actual initialize routines
    !    call InitializeIPDv02p1(driver, importState, exportState, clock, rc=rc)
    call NUOPC_CompSearchPhaseMap(driver, methodflag=ESMF_METHOD_INITIALIZE, & 
         phaseLabel=label_ExternalAdvertise, phaseIndex=phase, rc=rc)
    call ESMF_GridCompInitialize(driver, phase=phase, clock=clock, & 
         importState=importState, exportState=exportState, userRc=urc, rc=rc)
    
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=subname)) return  ! bail out
    call ESMF_LogWrite(trim(subname)//"v02p03: called", ESMF_LOGMSG_INFO)
    !    call InitializeIPDv02p3(driver, importState, exportState, clock, rc=rc)
    call NUOPC_CompSearchPhaseMap(driver, methodflag=ESMF_METHOD_INITIALIZE, & 
         phaseLabel=label_ExternalRealize, phaseIndex=phase, rc=rc) ! check rc 
    
    call ESMF_GridCompInitialize(driver, phase=phase, clock=clock, & 
         importState=importState, exportState=exportState, userRc=urc, rc=rc) ! check rc and urc        
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=subname)) return  ! bail out
    call ESMF_LogWrite(trim(subname)//"v02p05: called", ESMF_LOGMSG_INFO)
!    call InitializeIPDv02p5(driver, importState, exportState, clock, rc=rc)
    call NUOPC_CompSearchPhaseMap(driver, methodflag=ESMF_METHOD_INITIALIZE, & 
         phaseLabel=label_ExternalDataInit, phaseIndex=phase, rc=rc) ! check rc 

    call ESMF_GridCompInitialize(driver, phase=phase, clock=clock, & 
         importState=importState, exportState=exportState, userRc=urc, rc=rc) ! check rc and urc    
  
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=subname)) return  ! bail out
      
    ! query Component for the internal State
!    nullify(is%wrap)
!#ifdef ESMF_NO_F2018ASSUMEDTYPE
!    call ESMF_UserCompGetInternalState(driver, label_InternalState, is, rc)
!#else
!    call ESMF_UserCompGetInternalState(driver, label_InternalState, is, rc=rc)
!#endif
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, file=subname, rcToReturn=rc)) &
!      return  ! bail out

    ! check for dead-lock condition in hierarchical data dependency resolution
!    if (.not.is%wrap%dataDepAllComplete) then
!      ! this indicates a dead-lock condition
!      call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, &
!        msg="Initialize data-dependency resolution loop "// &
!        "has entered a dead-lock situation.", &
!        line=__LINE__, file=subname, rcToReturn=rc)
!      return  ! bail out
!    endif

    ! extro
    call NUOPC_LogExtro(name, rName, verbosity, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=subname)) return  ! bail out

    ! handle profiling
    if (btest(profiling,0)) then
      call ESMF_TraceRegionExit(rName, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=subname, rcToReturn=rc)) &
        return  ! bail out
    endif
    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

  end subroutine
#endif

  !================================================================================

  subroutine SetModelServices(ensemble_driver, rc)

    use ESMF          , only : ESMF_GridComp, ESMF_VM, ESMF_Config, ESMF_Clock, ESMF_VMGet
    use ESMF          , only : ESMF_GridCompGet, ESMF_VMGet, ESMF_ConfigGetAttribute
    use ESMF          , only : ESMF_ConfigGetLen, ESMF_RC_NOT_VALID, ESMF_LogFoundAllocError
    use ESMF          , only : ESMF_LogSetError, ESMF_LogWrite, ESMF_LOGMSG_INFO
    use ESMF          , only : ESMF_GridCompSet, ESMF_SUCCESS, ESMF_METHOD_INITIALIZE, ESMF_RC_ARG_BAD
    use ESMF          , only : ESMF_CalendarSetDefault
    use ESMF          , only : ESMF_CALKIND_NOLEAP, ESMF_CALKIND_GREGORIAN
    use NUOPC         , only : NUOPC_CompAttributeGet, NUOPC_CompAttributeSet, NUOPC_CompAttributeAdd
    use NUOPC_Driver  , only : NUOPC_DriverAddComp
    use esm           , only : ESMSetServices => SetServices, ReadAttributes
    use esm_time_mod  , only : esm_time_clockInit
    ! input/output variables
    type(ESMF_GridComp)    :: ensemble_driver
    integer, intent(out)   :: rc

    ! local variables
    type(ESMF_VM)          :: vm
    type(ESMF_GridComp)    :: driver, gridcomptmp
    type(ESMF_Config)      :: config
    integer                :: n, n1, stat
    integer, pointer       :: petList(:)
    character(len=20)      :: model, prefix
    integer                :: petCount, i
    integer                :: localPet
    logical                :: is_set
    character(len=512)     :: diro
    character(len=512)     :: logfile
    integer                :: global_comm
    logical                :: read_restart
    character(len=CS)      :: read_restart_string
    integer                :: inst
    integer                :: number_of_members
    integer                :: ntasks_per_member
    integer                :: asyncio_ntasks  ! input from config
    integer                :: asyncio_stride  ! input from config
    integer                :: asyncio_tasks   ! local counter
    logical                :: inComp
    character(CL)          :: start_type     ! Type of startup
    character(len=7)       :: drvrinst
    character(len=5)       :: inst_suffix
    character(len=CL)      :: msgstr
    character(len=CL)      :: cvalue
    character(len=CL)      :: calendar
    character(len=*) , parameter :: start_type_start = "startup"
    character(len=*) , parameter :: start_type_cont  = "continue"
    character(len=*) , parameter :: start_type_brnch = "branch"
    character(len=*) , parameter :: subname = u_FILE_u//":SetModelServices)"
    !-------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)

    call ESMF_GridCompGet(ensemble_driver, config=config, vm=vm, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    !-------------------------------------------
    ! Initialize clocks
    !-------------------------------------------

    call ReadAttributes(ensemble_driver, config, "ALLCOMP_attributes::", rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ReadAttributes(ensemble_driver, config, "CLOCK_attributes::", rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompAttributeGet(ensemble_driver, 'calendar', calendar, rc=rc) 
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (calendar == 'NO_LEAP') then
       call ESMF_CalendarSetDefault(ESMF_CALKIND_NOLEAP, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    else if (calendar == 'GREGORIAN') then
       call ESMF_CalendarSetDefault(ESMF_CALKIND_GREGORIAN, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    else
       write (msgstr, *) "Only NO_LEAP and GREGORIAN calendars currently supported"
       call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
       return  ! bail out
    end if

    ! Check valid values of start type
    call NUOPC_CompAttributeGet(ensemble_driver, name="start_type", value=start_type, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if ((trim(start_type) /= start_type_start) .and.  &
        (trim(start_type) /= start_type_cont ) .and.  &
        (trim(start_type) /= start_type_brnch)) then
       write (msgstr, *) subname//': start_type invalid = '//trim(start_type)
       call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
       return
    end if

    if (trim(start_type) == trim(start_type_cont) .or. trim(start_type) == trim(start_type_brnch)) then
       read_restart = .true.
    else
       read_restart = .false.
    endif
    write(read_restart_string,*) read_restart

    ! Add read_restart to ensemble_driver attributes
    call NUOPC_CompAttributeAdd(ensemble_driver, attrList=(/'read_restart'/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompAttributeSet(ensemble_driver, name='read_restart', value=trim(read_restart_string), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    !-------------------------------------------
    ! Extract the config object from the ensemble_driver
    !-------------------------------------------

    call ReadAttributes(ensemble_driver, config, "PELAYOUT_attributes::", rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ReadAttributes(ensemble_driver, config, "DRIVER_attributes::", rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    !-------------------------------------------
    ! Determine number of ensemble members and the number of tasks per member
    !-------------------------------------------

    call NUOPC_CompAttributeGet(ensemble_driver, name="ninst", value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) number_of_members

    call NUOPC_CompAttributeGet(ensemble_driver, name="asyncio_ntasks", value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) asyncio_ntasks

    call NUOPC_CompAttributeGet(ensemble_driver, name="asyncio_stride", value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) asyncio_stride

    allocate(asyncio_petlist(asyncio_ntasks))
    do i=1,asyncio_ntasks
       asyncio_petlist(i) = asyncio_stride*i
    enddo

    call ESMF_VMGet(vm, localPet=localPet, PetCount=PetCount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ntasks_per_member = (PetCount - asyncio_ntasks)/number_of_members
    if(ntasks_per_member*number_of_members + asyncio_ntasks .ne. PetCount) then
       write (msgstr,'(a,i5,a,i5,a,i3,a)') &
            "PetCount - asyncio_ntasks (",PetCount," - ", asyncio_ntasks, ") must be evenly divisable by number of members (",number_of_members,")"
       call ESMF_LogSetError(ESMF_RC_ARG_BAD, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
       return
    endif

    !-------------------------------------------
    ! Loop over number of ensemble members
    !-------------------------------------------

    allocate(petList(ntasks_per_member))
    asyncio_tasks = 0
    do inst=1,number_of_members

       ! Determine pet list for driver instance
       petList(1) = (inst-1) * ntasks_per_member
       do n=2,ntasks_per_member
          petList(n) = petList(n-1) + 1
          if(asyncio_ntasks > 0) then
             if(petlist(n) == asyncio_petlist(asyncio_tasks+1)) then
                asyncio_tasks = asyncio_tasks + 1
                petList(n) = petList(n) + 1
             endif
          endif
       enddo

       ! Add driver instance to ensemble driver
       write(drvrinst,'(a,i4.4)') "ESM",inst
       call NUOPC_DriverAddComp(ensemble_driver, drvrinst, ESMSetServices, petList=petList, comp=gridcomptmp, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       inComp = .false.
       do n=1,ntasks_per_member
          if(localpet == petlist(n)) then
             inComp = .true.
             exit
          endif
       enddo
       if (inComp) then

          driver = gridcomptmp

          if(number_of_members > 1) then
             call NUOPC_CompAttributeAdd(driver, attrList=(/'inst_suffix'/), rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             write(inst_suffix,'(a,i4.4)') '_',inst
             call NUOPC_CompAttributeSet(driver, name='inst_suffix', value=inst_suffix, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
          else
             inst_suffix = ''
          endif

          ! Set the driver instance attributes
          call NUOPC_CompAttributeAdd(driver, attrList=(/'read_restart'/), rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call NUOPC_CompAttributeSet(driver, name='read_restart', value=trim(read_restart_string), rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          call ReadAttributes(driver, config, "CLOCK_attributes::", rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          call ReadAttributes(driver, config, "DRIVER_attributes::", rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          call ReadAttributes(driver, config, "DRV_modelio::", rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          ! Set the driver log to the driver task 0 
          if (mod(localPet, ntasks_per_member) == 0) then
             call NUOPC_CompAttributeGet(driver, name="diro", value=diro, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             call NUOPC_CompAttributeGet(driver, name="logfile", value=logfile, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             open (newunit=logunit,file=trim(diro)//"/"//trim(logfile))
             mastertask = .true.
          else
             logUnit = shrlogunit
             mastertask = .false.
          endif
          call shr_file_setLogUnit (logunit)
       endif

       ! Create a clock for each driver instance
       call esm_time_clockInit(ensemble_driver, driver, logunit, mastertask, inComp, rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

    enddo

    deallocate(petList)
    
    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

  end subroutine SetModelServices
  
!  recursive subroutine InitializeIO(ensemble_driver, importState, exportState, clock, rc)
   subroutine InitializeIO(ensemble_driver, rc)
    use ESMF, only: ESMF_GridComp, ESMF_LOGMSG_INFO, ESMF_LogWrite
    use ESMF, only: ESMF_SUCCESS, ESMF_VM, ESMF_GridCompGet, ESMF_VMGet
    use ESMF, only: ESMF_CONFIG, ESMF_GridCompIsPetLocal, ESMF_State, ESMF_Clock
    use NUOPC, only: NUOPC_CompAttributeGet
    use NUOPC_DRIVER, only: NUOPC_DriverGetComp
    use shr_pio_mod   , only: shr_pio_init, shr_pio_component_init

    type(ESMF_GridComp) :: ensemble_driver
!    type(ESMF_State)    :: importState, exportState
!    type(ESMF_Clock)    :: clock
    integer, intent(out) :: rc

    type(ESMF_VM) :: ensemble_vm
    type(ESMF_CONFIG) :: config
    type(ESMF_GridComp), pointer :: dcomp(:)
    integer :: Global_Comm
    integer :: number_of_members
    logical :: asyncio_task
    integer :: i
    integer :: iam
    character(len=CL) :: cvalue
    character(len=*), parameter :: subname=u_FILE_u//"InitializeIO"

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)

    call ESMF_GridCompGet(ensemble_driver, vm=ensemble_vm, config=config, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompAttributeGet(ensemble_driver, name="ninst", value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue, *) number_of_members
    allocate(dcomp(number_of_members))
    
    call ESMF_VMGet(ensemble_vm, localpet=iam, mpiCommunicator=Global_Comm, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    
    ! Initialize PIO
    ! This reads in the pio parameters that are independent of component    
    call shr_pio_init(ensemble_driver, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    ! Read in component dependent PIO parameters and initialize
    ! IO systems
    nullify(dcomp)
    call NUOPC_DriverGetComp(ensemble_driver, complist=dcomp, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    asyncio_task = .false.
    do i=1,size(asyncio_petlist)
       if(iam == asyncio_petlist(i)) then
          asyncio_task = .true.
          exit
       end if
    enddo
    do i=1,size(dcomp)
       if (ESMF_GridCompIsPetLocal(dcomp(i), rc=rc) .or. asyncio_task) then
          call shr_pio_component_init(dcomp(i), Global_Comm, asyncio_petlist, rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       endif
    enddo
    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
  end subroutine InitializeIO
  
end module Ensemble_driver
