!> Module for the coupling of AeroDyn and InflowWind
!! Also contains routine to couple with a "Fake ElastoDyn" (minimalistic structural solver) 
module AeroDyn_Inflow
   use NWTC_Library
   use AeroDyn_Inflow_Types
   use AeroDyn_Types
   use AeroDyn, only: AD_Init, AD_ReInit, AD_CalcOutput, AD_UpdateStates, AD_End
   use AeroDyn_IO, only: AD_SetVTKSurface
   use InflowWind, only: InflowWind_Init, InflowWind_CalcOutput, InflowWind_End

   implicit none

   private

   type(ProgDesc), parameter  :: ADI_Ver = ProgDesc( 'ADI', '', '' )

   public   :: ADI_Init 
   public   :: ADI_ReInit 
   public   :: ADI_End
   public   :: ADI_CalcOutput
   public   :: ADI_UpdateStates

   ! Convenient routines for driver
   public   :: concatOutputHeaders
   public   :: Init_MeshMap_For_ADI
   public   :: Set_Inputs_For_ADI

   real(ReKi), parameter :: myNaN = -99.9_ReKi
contains

!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is called at the start of the simulation to perform initialization steps.
!! The parameters are set here and not changed during the simulation.
!! The initial states and initial guess for the input are defined.
subroutine ADI_Init(InitInp, u, p, x, xd, z, OtherState, y, m, Interval, InitOut, errStat, errMsg)
   type(ADI_InitInputType),         intent(inout)  :: InitInp        !< Input data for initialization routine  (inout so we can use MOVE_ALLOC)
   type(ADI_InputType),             intent(  out)  :: u              !< An initial guess for the input; input mesh must be defined
   type(ADI_ParameterType),         intent(  out)  :: p              !< Parameters
   type(ADI_ContinuousStateType),   intent(  out)  :: x              !< Initial continuous states
   type(ADI_DiscreteStateType),     intent(  out)  :: xd             !< Initial discrete states
   type(ADI_ConstraintStateType),   intent(  out)  :: z              !< Initial guess of the constraint states
   type(ADI_OtherStateType),        intent(  out)  :: OtherState     !< Initial other states
   type(ADI_OutputType),            intent(  out)  :: y              !< Initial system outputs (outputs are not calculated;
   type(ADI_MiscVarType),           intent(  out)  :: m              !< Initial misc/optimization variables
   real(DbKi),                      intent(inout)  :: interval       !< Coupling interval in seconds
   type(ADI_InitOutputType),        intent(inout)  :: InitOut        !< Output for initialization routine. NOTE: inout to allow for reinit?
   integer(IntKi),                  intent(  out)  :: errStat        !< Error status of the operation
   character(*),                    intent(  out)  :: errMsg         !< Error message if errStat /= ErrID_None
   ! Local variables
   type(InflowWind_InitOutputType) :: InitOut_IW  ! Output data from initialization
   type(AD_InitOutputType) :: InitOut_AD  ! Output data from initialization
   integer(IntKi)          :: errStat2       ! temporary error status of the operation
   character(errMsgLen)    :: errMsg2        ! temporary error message

   ! Initialize variables for this routine
   errStat = ErrID_None
   errMsg  = ""

   ! Initialize the NWTC Subroutine Library
   call NWTC_Init( EchoLibVer=.FALSE. )

   ! Display the module information
   call DispNVD( ADI_Ver )

   ! Clear writeoutputs
   if (allocated(InitOut%WriteOutputHdr)) deallocate(InitOut%WriteOutputHdr)
   if (allocated(InitOut%WriteOutputUnt)) deallocate(InitOut%WriteOutputUnt)

   ! Set parameters
   p%dt         = interval
   p%storeHHVel = InitInp%storeHHVel
   p%WrVTK      = InitInp%WrVTK
   p%MHK        = InitInp%AD%MHK
   p%WtrDpth    = InitInp%AD%WtrDpth

   ! --- Initialize Inflow Wind 
   call ADI_InitInflowWind(InitInp%RootName, InitInp%IW_InitInp, u%AD, OtherState%AD, m%IW, Interval, InitOut_IW, errStat2, errMsg2); if (Failed()) return
   ! Concatenate AD outputs to IW outputs
   call concatOutputHeaders(InitOut%WriteOutputHdr, InitOut%WriteOutputUnt, InitOut_IW%WriteOutputHdr, InitOut_IW%WriteOutputUnt, errStat2, errMsg2); if(Failed()) return

   ! --- Initialize AeroDyn
   ! Link InflowWind's FlowField to AeroDyn's FlowField
   InitInp%AD%FlowField => InitOut_IW%FlowField

   call AD_Init(InitInp%AD, u%AD, p%AD, x%AD, xd%AD, z%AD, OtherState%AD, y%AD, m%AD, Interval, InitOut_AD, errStat2, errMsg2); if (Failed()) return
   InitOut%Ver = InitOut_AD%ver
   ! Add writeoutput units and headers to driver, same for all cases and rotors!
   !TODO: this header is too short if we add more rotors.  Should also add a rotor identifier
   call concatOutputHeaders(InitOut%WriteOutputHdr, InitOut%WriteOutputUnt, InitOut_AD%rotors(1)%WriteOutputHdr, InitOut_AD%rotors(1)%WriteOutputUnt, errStat2, errMsg2); if(Failed()) return

   ! --- Initialize grouped outputs
   !TODO: assumes one rotor
   p%NumOuts = p%AD%rotors(1)%NumOuts + p%AD%rotors(1)%BldNd_TotNumOuts + m%IW%p%NumOuts
   call AllocAry(y%WriteOutput, p%NumOuts, 'WriteOutput', errStat2, errMsg2); if (Failed()) return

   ! --- Initialize outputs
   call AllocAry(y%IW_WriteOutput, size(m%IW%y%WriteOutput),'IW_WriteOutput', errStat2, errMsg2); if(Failed()) return
   y%IW_WriteOutput = myNaN
   if (p%storeHHVel) then
      call AllocAry(y%HHVel, 3, size(InitInp%AD%rotors), 'HHVel', errStat2, errMsg2); if(Failed()) return
      y%HHVel= myNaN
   else
      call AllocAry(y%HHVel, 0,          0             , 'HHVel', errStat2, errMsg2); if(Failed()) return
   endif

   ! --- Initialize VTK
   if (p%WrVTK>0) then
      call AD_SetVTKSurface(InitOut_AD, u%AD, m%VTK_Surfaces, errStat2, errMsg2); if(Failed()) return
   endif

    call cleanup()

contains

   subroutine cleanup()
      call AD_DestroyInitInput (InitInp%AD, errStat2, errMsg2)   
      call AD_DestroyInitOutput(InitOut_AD, errStat2, errMsg2)      
      call InflowWind_DestroyInitOutput(InitOut_IW, errStat2, errMsg2)
   end subroutine cleanup

   logical function Failed()
        call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'ADI_Init') 
        Failed =  errStat >= AbortErrLev
         if (Failed) call cleanup()
   end function Failed

end subroutine ADI_Init
!----------------------------------------------------------------------------------------------------------------------------------
!> ReInit
subroutine ADI_ReInit(p, x, xd, z, OtherState, m, Interval, errStat, errMsg)
   type(ADI_ParameterType),         intent(in   )  :: p              !< Parameters
   type(ADI_ContinuousStateType),   intent(inout)  :: x              !< Initial continuous states
   type(ADI_DiscreteStateType),     intent(inout)  :: xd             !< Initial discrete states
   type(ADI_ConstraintStateType),   intent(inout)  :: z              !< Initial guess of the constraint states
   type(ADI_OtherStateType),        intent(inout)  :: OtherState     !< Initial other states
   type(ADI_MiscVarType),           intent(inout)  :: m              !< Initial misc/optimization variables
   real(DbKi),                      intent(inout)  :: interval       !< Coupling interval in seconds
   integer(IntKi),                  intent(  out)  :: errStat        !< Error status of the operation
   character(*),                    intent(  out)  :: errMsg         !< Error message if errStat /= ErrID_None
   ! Local variables
   integer(IntKi)          :: errStat2       ! temporary error status of the operation
   character(errMsgLen)    :: errMsg2        ! temporary error message
   errStat = ErrID_None
   errMsg  = ""

   ! Reinitialize AeroDyn without reopening input file
   call AD_ReInit(p%AD, x%AD, xd%AD, z%AD, OtherState%AD, m%AD, Interval, errStat2, errMsg2); if(Failed()) return
   ! Set parameters
   !p%dt = interval ! dt shouldn't change
contains

   logical function Failed()
        call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'ADI_ReInit') 
        Failed =  errStat >= AbortErrLev
   end function Failed

end subroutine ADI_ReInit
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is called at the end of the simulation.
subroutine ADI_End( u, p, x, xd, z, OtherState, y, m, errStat, errMsg )
   type(ADI_InputType),             intent(inout)  :: u(:)        !< System inputs NOTE: used to be allocatable
   type(ADI_ParameterType),         intent(inout)  :: p           !< Parameters
   type(ADI_ContinuousStateType),   intent(inout)  :: x           !< Continuous states
   type(ADI_DiscreteStateType),     intent(inout)  :: xd          !< Discrete states
   type(ADI_ConstraintStateType),   intent(inout)  :: z           !< Constraint states
   type(ADI_OtherStateType),        intent(inout)  :: OtherState  !< Other states
   type(ADI_OutputType),            intent(inout)  :: y           !< System outputs
   type(ADI_MiscVarType),           intent(inout)  :: m           !< Misc/optimization variables
   integer(IntKi),                  intent(  out)  :: errStat     !< Error status of the operation
   character(*),                    intent(  out)  :: errMsg      !< Error message if errStat /= ErrID_None
   integer(IntKi) :: i
   errStat = ErrID_None
   errMsg  = ""

   ! End modules
   call AD_End(u(1)%AD, p%AD, x%AD, xd%AD, z%AD, OtherState%AD, y%AD, m%AD, ErrStat, ErrMsg)
   call InflowWind_End(m%IW%u, m%IW%p, m%IW%x, m%IW%xd, m%IW%z, m%IW%OtherSt, m%IW%y, m%IW%m, ErrStat, ErrMsg)

   ! Destroy the input data:
   !if (allocated(u)) then
      do i=1,size(u)
         call ADI_DestroyInput( u(i), errStat, errMsg )
      enddo
   !endif

   ! Destroy the parameter data:
   call ADI_DestroyParam( p, errStat, errMsg )

   ! Destroy the state data:
   call ADI_DestroyContState(   x,           errStat, errMsg )
   call ADI_DestroyDiscState(   xd,          errStat, errMsg )
   call ADI_DestroyConstrState( z,           errStat, errMsg )
   call ADI_DestroyOtherState(  OtherState,  errStat, errMsg )
   call ADI_DestroyMisc(        m,           errStat, errMsg )

   ! Destroy the output data:
   call ADI_DestroyOutput( y, errStat, errMsg )

end subroutine ADI_End
!----------------------------------------------------------------------------------------------------------------------------------
!> Loose coupling routine for solving for constraint states, integrating continuous states, and updating discrete and other states.
!! Continuous, constraint, discrete, and other states are updated for t + Interval
subroutine ADI_UpdateStates( t, n, u, utimes, p, x, xd, z, OtherState, m, errStat, errMsg)
   real(DbKi),                      intent(in   )  :: t           !< Current simulation time in seconds
   integer(IntKi),                  intent(in   )  :: n           !< Current simulation time step n = 0,1,...
   type(ADI_InputType),             intent(inout)  :: u(:)        !< Inputs at utimes (out only for mesh record-keeping in ExtrapInterp routine)
   real(DbKi),                      intent(in   )  :: utimes(:)   !< Times associated with u(:), in seconds
   type(ADI_ParameterType),         intent(in   )  :: p           !< Parameters
   type(ADI_ContinuousStateType),   intent(inout)  :: x           !< Input: Continuous states at t; Output: at t+DTaero
   type(ADI_DiscreteStateType),     intent(inout)  :: xd          !< Input: Discrete states at t;   Output: at t+DTaero
   type(ADI_ConstraintStateType),   intent(inout)  :: z           !< Input: Constraint states at t; Output: at t+DTaero
   type(ADI_OtherStateType),        intent(inout)  :: OtherState  !< Input: Other states at t;      Output: at t+DTaero
   type(ADI_MiscVarType),           intent(inout)  :: m           !< Misc/optimization variables
   integer(IntKi),                  intent(  out)  :: errStat     !< Error status of the operation
   character(*),                    intent(  out)  :: errMsg      !< Error message if errStat /= ErrID_None
   ! local variables
   integer :: it ! Index on times
   type(AD_InputType)  ::  u_AD(size(utimes))
   integer(IntKi)                :: errStat2                                                           ! temporary Error status
   character(errMsgLen)          :: errMsg2                                                            ! temporary Error message
   !type(ADI_InputType)           :: uInterp     ! Interpolated/Extrapolated input
   errStat = ErrID_None
   errMsg  = ""

   ! Compute InflowWind inputs for each time
   do it=1,size(utimes)
      call AD_CopyInput(u(it)%AD,u_AD(it),MESH_NEWCOPY,ErrStat2,ErrMsg2); if(Failed()) return
   enddo

   ! Get state variables at next step: INPUT at step nt - 1, OUTPUT at step nt
   call AD_UpdateStates(t, n, u_AD(:), utimes(:), p%AD, x%AD, xd%AD, z%AD, OtherState%AD, m%AD, errStat2, errMsg2); if(Failed()) return

   call CleanUp()

contains

   subroutine CleanUp()
      !call ADI_DestroyConstrState(z_guess, errStat2, errMsg2); if(Failed()) return
      do it=1,size(utimes)
         call AD_DestroyInput(u_AD(it), errStat2, errMsg2)  ! ignore errors here
      enddo
   end subroutine

   logical function Failed()
      call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'ADI_UpdateStates') 
      Failed =  errStat >= AbortErrLev
      if (Failed) call CleanUp()
   end function Failed

end subroutine ADI_UpdateStates
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine for computing outputs, used in both loose and tight coupling.
subroutine ADI_CalcOutput(t, u, p, x, xd, z, OtherState, y, m, errStat, errMsg)
   use IfW_FlowField, only: IfW_FlowField_GetVelAcc
   real(DbKi),                      intent(in   )  :: t           !< Current simulation time in seconds
   type(ADI_InputType),             intent(inout)  :: u           !< Inputs at Time t  ! NOTE: set as in-out since "Inflow" needs to be set
   type(ADI_ParameterType),         intent(in   )  :: p           !< Parameters
   type(ADI_ContinuousStateType),   intent(in   )  :: x           !< Continuous states at t
   type(ADI_DiscreteStateType),     intent(in   )  :: xd          !< Discrete states at t
   type(ADI_ConstraintStateType),   intent(in   )  :: z           !< Constraint states at t
   type(ADI_OtherStateType),        intent(in   )  :: OtherState  !< Other states at t
   type(ADI_OutputType),            intent(inout)  :: y           !< Outputs computed at t (Input only so that mesh con-
                                                                  !!   nectivity information does not have to be recalculated)
   type(ADI_MiscVarType),           intent(inout)  :: m           !< Misc/optimization variables
   integer(IntKi),                  intent(  out)  :: errStat     !< Error status of the operation
   character(*),                    intent(  out)  :: errMsg      !< Error message if errStat /= ErrID_None

   integer(IntKi)                :: errStat2
   character(errMsgLen)          :: errMsg2
   integer(IntKi)                :: node
   character(*), parameter       :: RoutineName = 'ADI_CalcOutput'
   integer :: iWT
   errStat = ErrID_None
   errMsg  = ""

   !----------------------------------------------------------------------------
   ! Calculate InflowWind outputs if module was initialized
   !----------------------------------------------------------------------------

   if (m%IW%CompInflow == 1) then
      call InflowWind_CalcOutput(t, m%IW%u, m%IW%p, m%IW%x, m%IW%xd, m%IW%z,  &
                                 m%IW%OtherSt, m%IW%y, m%IW%m, errStat2, errMsg2)
      if(Failed()) return

      ! Copy InflowWind outputs to ADI outputs
      y%IW_WriteOutput(:) = m%IW%y%WriteOutput(:)
   end if

   !----------------------------------------------------------------------------
   ! Calculate aerodyn output
   !----------------------------------------------------------------------------

   ! Calculate outputs at t
   call AD_CalcOutput(t, u%AD, p%AD, x%AD, xd%AD, z%AD, OtherState%AD, y%AD, m%AD, errStat2, errMsg2)
   if(Failed()) return

   ! --- Outputs for driver
   y%PLExp = m%IW%PLExp

   ! --- Set outputs
   !TODO: this assumes one rotor!!!
   associate(AD_NumOuts => p%AD%rotors(1)%NumOuts + p%AD%rotors(1)%BldNd_TotNumOuts, &
             IW_NumOuts => m%IW%p%NumOuts)
      y%WriteOutput(1:IW_NumOuts) = y%IW_WriteOutput(1:IW_NumOuts)
      y%WriteOutput(IW_NumOuts+1:p%NumOuts) = y%AD%rotors(1)%WriteOutput(1:AD_NumOuts)
   end associate

   !----------------------------------------------------------------------------
   ! Store hub height velocity calculated in CalcOutput
   !----------------------------------------------------------------------------

   if (p%storeHHVel) then
      do iWT = 1, size(u%AD%rotors)
         y%HHVel(:,iWT) = m%AD%Inflow(1)%RotInflow(iWT)%InflowOnHub(:,1)
      end do
   endif

contains

   subroutine CleanUp()
      !call ADI_DestroyConstrState(z_guess, errStat2, errMsg2); if(Failed()) return
   end subroutine

   logical function Failed()
      call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'ADI_CalcOutput') 
      Failed =  errStat >= AbortErrLev
      if (Failed) call CleanUp()
   end function Failed
end subroutine ADI_CalcOutput
!----------------------------------------------------------------------------------------------------------------------------------
!>
subroutine ADI_InitInflowWind(Root, i_IW, u_AD, o_AD, IW, dt, InitOutData, errStat, errMsg)
   use InflowWind, only: InflowWind_Init
   use InflowWind_IO, only: IfW_SteadyWind_Init
   use IfW_FlowField, only: IfW_UniformField_CalcAccel
   character(*),                 intent(in   ) :: Root          ! Rootname for input files
   type(ADI_IW_InputData),       intent(in   ) :: i_IW          ! Inflow Wind "pseudo init input" data
   type(AD_InputType),           intent(in   ) :: u_AD          ! AeroDyn data 
   type(AD_OtherStateType),      intent(in   ) :: o_AD          ! AeroDyn data 
   type(ADI_InflowWindData),     intent(inout) :: IW            ! InflowWind data 
   real(DbKi),                   intent(inout) :: dt            ! interval
   type(InflowWind_InitOutputType), intent(out) :: InitOutData  ! Output data from initialization
   integer(IntKi)              , intent(  out) :: errStat       ! Status of error message
   character(*)                , intent(  out) :: errMsg        ! Error message if errStat /= ErrID_None

   integer(IntKi)                   :: errStat2             ! local status of error message
   character(errMsgLen)             :: errMsg2              ! local error message if errStat /= ErrID_None
   type(InflowWind_InitInputType)   :: InitInData           ! Input data for initialization
   Type(Steady_InitInputType)       :: Steady_InitInput
   Type(WindFileDat)                :: FileDat
   errStat = ErrID_None
   errMsg  = ''

   ! --- Init InflowWind
   if (i_IW%CompInflow == 0) then
      ! Initialze only the flow field with steady wind
      allocate(InitOutData%WriteOutputHdr(0))
      allocate(InitOutData%WriteOutputUnt(0))
      allocate(IW%y%WriteOutput(0))
      Steady_InitInput%HWindSpeed = i_IW%HWindSpeed
      Steady_InitInput%RefHt = i_IW%RefHt
      Steady_InitInput%PLExp = i_IW%PLExp
      allocate(IW%p%FlowField)
      IW%p%FlowField%PropagationDir = 0.0_ReKi
      IW%p%FlowField%VFlowAngle = 0.0_ReKi
      IW%p%FlowField%RotateWindBox = .false.
      IW%p%FlowField%FieldType = Uniform_FieldType
      IW%p%FlowField%RefPosition = [0.0_ReKi, 0.0_ReKi, i_IW%RefHt]
      InitOutData%FlowField => IW%p%FlowField
      call IfW_SteadyWind_Init(Steady_InitInput, 0, IW%p%FlowField%Uniform, &
                               FileDat, errStat2, errMsg2)
      if(Failed()) return
      if (i_IW%MHK == MHK_FixedBottom .or. i_IW%MHK == MHK_FLoating) then
         call IfW_UniformField_CalcAccel(IW%p%FlowField%Uniform, errStat2, errMsg2)
         if(Failed()) return
         IW%p%FlowField%AccFieldValid = .true.
      end if
   else
      ! Initialze InflowWind module
      InitInData%InputFileName    = i_IW%InputFile
      InitInData%Linearize        = i_IW%Linearize
      InitInData%FilePassingMethod= i_IW%FilePassingMethod
      InitInData%NumWindPoints = 1
      if (i_IW%FilePassingMethod == 1_IntKi) then     ! passing input file as an FileInfoType structure
         call NWTC_Library_Copyfileinfotype( i_IW%PassedFileInfo, InitInData%PassedFileInfo, MESH_NEWCOPY, errStat2, errMsg2 ); if (Failed()) return
      elseif (i_IW%FilePassingMethod == 2_IntKi) then ! passing input file as an IfW_InputFile structure
         call InflowWind_CopyInputFile( i_IW%PassedFileData, InitInData%PassedFileData, MESH_NEWCOPY, errStat2, errMsg2 ); if (Failed()) return
      endif
      InitInData%RootName         = trim(Root)//'.IfW'
      InitInData%MHK              = i_IW%MHK
      ! OLAF might be used in AD, in which case we need to allow out of bounds for some calcs. To do that
      ! the average values for the entire wind profile must be calculated and stored (we don't know if OLAF
      ! is used until after AD_Init below).
      InitInData%BoxExceedAllow = .true.
      CALL InflowWind_Init( InitInData, IW%u, IW%p, &
                     IW%x, IW%xd, IW%z, IW%OtherSt, &
                     IW%y, IW%m, dt,  InitOutData, errStat2, errMsg2 )
      if(Failed()) return
   endif

   ! --- Store main init input data (data that don't use InfloWind directly)
   IW%CompInflow = i_IW%CompInflow
   IW%HWindSpeed = i_IW%HWindSpeed
   IW%RefHt      = i_IW%RefHt
   IW%PLExp      = i_IW%PLExp

   call cleanup()

contains
   subroutine cleanup()
      call InflowWind_DestroyInitInput( InitInData, errStat2, errMsg2 )   
   end subroutine cleanup

   logical function Failed()
      CALL SetErrStat( errStat2, errMsg2, errStat, errMsg, 'ADI_InitInflowWind' )
      Failed = errStat >= AbortErrLev
      if (Failed) call cleanup()
   end function Failed
end subroutine ADI_InitInflowWind
!----------------------------------------------------------------------------------------------------------------------------------
!> Concatenate new output channels info to the extisting ones in the driver
subroutine concatOutputHeaders(WriteOutputHdr0, WriteOutputUnt0, WriteOutputHdr, WriteOutputUnt, errStat, errMsg)
   character(ChanLen), dimension(:), allocatable, intent(inout) ::  WriteOutputHdr0 !< Channel headers
   character(ChanLen), dimension(:), allocatable, intent(inout) ::  WriteOutputUnt0 !< Channel units
   character(ChanLen), dimension(:), allocatable, intent(inout) ::  WriteOutputHdr !< Channel headers
   character(ChanLen), dimension(:), allocatable, intent(inout) ::  WriteOutputUnt !< Channel units
   integer(IntKi)              , intent(  out) :: errStat       !< Status of error message
   character(*)                , intent(  out) :: errMsg        !< Error message if errStat /= ErrID_None
   ! Locals
   character(ChanLen), allocatable :: TmpHdr(:)
   character(ChanLen), allocatable :: TmpUnt(:)
   integer :: nOld, nAdd
   errStat = ErrID_None
   errMsg  = ''
   !print*,'>>> Concat',allocated(WriteOutputHdr0), allocated(WriteOutputUnt0), allocated(WriteOutputHdr), allocated(WriteOutputUnt)
   if (.not.allocated(WriteOutputHdr)) return
   if (.not.allocated(WriteOutputHdr0)) then
      call move_alloc(WriteOutputHdr, WriteOutputHdr0)
      call move_alloc(WriteOutputUnt, WriteOutputUnt0)   
   else
      nOld = size(WriteOutputHdr0)
      nAdd = size(WriteOutputHdr)

      call move_alloc(WriteOutputHdr0, TmpHdr)
      call move_alloc(WriteOutputUnt0, TmpUnt)   

      allocate(WriteOutputHdr0(nOld+nAdd))
      allocate(WriteOutputUnt0(nOld+nAdd))
      WriteOutputHdr0(1:nOld) = TmpHdr
      WriteOutputUnt0(1:nOld) = TmpUnt
      WriteOutputHdr0(nOld+1:nOld+nAdd) = WriteOutputHdr
      WriteOutputUnt0(nOld+1:nOld+nAdd) = WriteOutputUnt
      deallocate(TmpHdr)
      deallocate(TmpUnt)
   endif
end subroutine concatOutputHeaders
! --------------------------------------------------------------------------------}
! --- ROUTINES RELEVANT FOR COUPLING WITH "FED": Fake ElastoDyn 
! --------------------------------------------------------------------------------{
!> Initialize the mesh mappings between the structure and aerodyn
!! Also adjust the tower mesh so that is is aligned with the tower base and tower top
!! Similar to FAST_Solver.f90, InitModuleMappings
subroutine Init_MeshMap_For_ADI(FED, p, uAD, errStat, errMsg)
   type(FED_Data), target,       intent(inout) :: FED       ! Elastic wind turbine data (Fake ElastoDyn)
   type(ADI_ParameterType),      intent(   in) :: p             ! Parameters
   type(AD_InputType),           intent(inout) :: uAD           ! AeroDyn input data 
   integer(IntKi)              , intent(  out) :: errStat       ! Status of error message
   character(*)                , intent(  out) :: errMsg        ! Error message if errStat /= ErrID_None
   ! locals
   real(ReKi)            :: pos(3), Pbase(3), Ptop(3), DeltaP(3)
   real(R8Ki)            :: orientation(3,3)
   real(ReKi)            :: twrHeightAD , twrHeight
   real(ReKi)            :: zBar ! dimensionsless tower height
   integer(IntKi)        :: iWT, iB, i
   integer(IntKi)        :: errStat2      ! local status of error message
   character(ErrMsgLen)  :: errMsg2       ! local error message if errStat /= ErrID_None
   type(RotFED), pointer :: y_ED ! Alias to shorten notation
   errStat = ErrID_None
   errMsg  = ''

   ! --- Create Mappings from structure to AeroDyn
   do iWT=1,size(FED%WT)
      y_ED => FED%WT(iWT)
      ! hub 2 hubAD
      call MeshMapCreate(y_ED%HubPtMotion, uAD%rotors(iWT)%hubMotion, y_ED%ED_P_2_AD_P_H, errStat2, errMsg2); if(Failed())return

      ! nac 2 nacAD
      call MeshMapCreate(y_ED%NacelleMotion, uAD%rotors(iWT)%nacelleMotion, y_ED%ED_P_2_AD_P_N, errStat2, errMsg2); if(Failed())return

      ! nac 2 tfinAD
      if (uAD%rotors(iWT)%TFinMotion%Committed) then
         call MeshMapCreate(y_ED%NacelleMotion, uAD%rotors(iWT)%TFinMotion, y_ED%ED_P_2_AD_P_TF, errStat2, errMsg2); if(Failed())return
      endif

      ! bldroot 2 bldroot AD
      allocate(y_ED%ED_P_2_AD_P_R(y_ED%numBlades))
      do iB = 1, y_ED%numBlades
         call MeshMapCreate(y_ED%BladeRootMotion(iB), uAD%rotors(iWT)%BladeRootMotion(iB), y_ED%ED_P_2_AD_P_R(iB), errStat2, errMsg2); if(Failed())return
      enddo

      if (y_ED%rigidBlades) then
         ! TODO Only for Rigid
         ! AD bld root 2 AD blade line
         allocate(y_ED%AD_P_2_AD_L_B(y_ED%numBlades))
         do iB = 1, y_ED%numBlades
            call MeshMapCreate(uAD%rotors(iWT)%BladeRootMotion(iB), uAD%rotors(iWT)%BladeMotion(iB), y_ED%AD_P_2_AD_L_B(iB), errStat2, errMsg2); if(Failed())return
         enddo
      else
         print*,'>>> Init_MeshMap_For_ADI, TODO coupling with elastic blades'
         STOP
      endif

      if (uAD%rotors(iWT)%TowerMotion%nNodes>0) then
         if (y_ED%hasTower) then
            twrHeightAD=uAD%rotors(iWT)%TowerMotion%Position(3,uAD%rotors(iWT)%TowerMotion%nNodes)-uAD%rotors(iWT)%TowerMotion%Position(3,1)
            ! Check tower height
            if ( p%MHK==MHK_Floating ) then
               if (twrHeightAD>0) then
                  errStat=ErrID_Fatal
                  errMsg='First AeroDyn tower height should be larger than last AD tower height for a floating MHK turbine'
               endif
            else
               if (twrHeightAD<0) then
                  errStat=ErrID_Fatal
                  errMsg='First AeroDyn tower height should be smaller than last AD tower height'
               endif
            endif

            twrHeightAD=uAD%rotors(iWT)%TowerMotion%Position(3,uAD%rotors(iWT)%TowerMotion%nNodes) ! NOTE: assuming start a z=0
            if ( p%MHK==MHK_FixedBottom ) then
               twrHeightAD = twrHeightAD + p%WtrDpth
            elseif ( p%MHK==MHK_Floating ) then
               twrHeightAD = abs(twrHeightAD)
            endif

            twrHeight=TwoNorm(y_ED%NacelleMotion%Position(:,1) - y_ED%TwrPtMesh%Position(:,1)  )
            ! KEEP ME, in summary file
            !print*,'Tower Height',twrHeight, twrHeightAD
            if (abs(twrHeightAD-twrHeight)> twrHeight*0.1) then
               errStat=ErrID_Fatal
               errMsg='More than 10% difference between AeroDyn tower length ('//trim(num2lstr(twrHeightAD))//&
                  'm), and the distance from tower base to nacelle ('//trim(num2lstr(twrHeight))//'m) for turbine '//trim(num2lstr(iWT))
            endif

            ! Adjust tower position (AeroDyn return values assuming (0,0,0) for tower base
            Pbase = y_ED%TwrPtMesh%Position(:,1)
            Ptop = y_ED%NacelleMotion%Position(:,1)
            if ( p%MHK==MHK_Floating ) then
               DeltaP = Pbase-Ptop
            else
               DeltaP = Ptop-Pbase
            endif
            do i = 1, uAD%rotors(iWT)%TowerMotion%nNodes
               if ( p%MHK==MHK_FixedBottom ) then
                  zBar = (uAD%rotors(iWT)%TowerMotion%Position(3,i) + p%WtrDpth) / twrHeight
               else
                  zBar = uAD%rotors(iWT)%TowerMotion%Position(3,i)/twrHeight
               endif
               uAD%rotors(iWT)%TowerMotion%Position(:,i)= Pbase+ zBar * DeltaP
               uAD%rotors(iWT)%TowerMotion%RefOrientation(:,:,i)= y_ED%TwrPtMesh%RefOrientation(:,:,1)
            enddo
            ! Create AD tower base point mesh
            pos         = y_ED%TwrPtMesh%Position(:,1)
            orientation = y_ED%TwrPtMesh%RefOrientation(:,:,1)
            call Eye(orientation, errStat2, errMsg2)
            call CreateInputPointMesh(y_ED%TwrPtMeshAD, pos, orientation, errStat2, errMsg2, hasMotion=.True., hasLoads=.False.); if(Failed())return

            ! TowerBase to AD tower base
            call MeshMapCreate(y_ED%TwrPtMesh, y_ED%TwrPtMeshAD, y_ED%ED_P_2_AD_P_T, errStat2, errMsg2); if(Failed()) return

            ! AD TowerBase to AD tower line
            call MeshMapCreate(y_ED%TwrPtMeshAD, uAD%rotors(iWT)%TowerMotion, y_ED%AD_P_2_AD_L_T, errStat2, errMsg2); if(Failed()) return
         endif ! hasTower
      else
         ! Do Nothing for now
      endif

   enddo ! Loop on WT/rotors

contains

   logical function Failed()
      call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'Init_MeshMap_For_ADI')
      Failed = errStat >= AbortErrLev
   end function Failed
end subroutine Init_MeshMap_For_ADI
!----------------------------------------------------------------------------------------------------------------------------------
!> Set aerodyn inputs based on FED meshes
!  - set AD input meshes and inflow
subroutine Set_Inputs_For_ADI(u_ADI, FED, errStat, errMsg)
   type(ADI_InputType),          intent(inout) :: u_ADI     !< AeroDyn/InflowWind Data inputs
   type(FED_Data), target,       intent(inout) :: FED       !< Elastic wind turbine data (Fake ElastoDyn)
   integer(IntKi)              , intent(  out) :: errStat   !< Status of error message
   character(*)                , intent(  out) :: errMsg    !< Error message if errStat /= ErrID_None
   ! local variables
   integer(intKi)          :: iWT ! loop counter for rotors
   integer(intKi)          :: iB ! loop counter for blades
   integer(IntKi)          :: errStat2      ! local status of error message
   character(ErrMsgLen)    :: errMsg2       ! local error message if errStat /= ErrID_None
   type(RotFED), pointer :: y_ED ! Alias to shorten notation
   errStat = ErrID_None
   errMsg  = ""

   ! --- Transfer motion from "ED" to AeroDyn
   do iWT=1,size(FED%WT)
      y_ED => FED%WT(iWT)
      ! Hub 2 Hub AD 
      call Transfer_Point_to_Point(y_ED%HubPtMotion, u_ADI%AD%rotors(iWT)%hubMotion, y_ED%ED_P_2_AD_P_H, errStat2, errMsg2); if(Failed()) return

      ! Nac 2 Nac AD 
      call Transfer_Point_to_Point(y_ED%NacelleMotion, u_ADI%AD%rotors(iWT)%nacelleMotion, y_ED%ED_P_2_AD_P_N, errStat2, errMsg2); if(Failed()) return

      ! Nac 2 TailFin AD (Transfer ElastoDyn CM motion (taken as Nacelle) to AeroDyn tailfin ref point motion
      if (u_ADI%AD%rotors(iWT)%TFinMotion%Committed) then
         call Transfer_Point_to_Point( y_ED%NacelleMotion, u_ADI%AD%rotors(IWT)%TFinMotion, y_ED%ED_P_2_AD_P_TF, errStat2, errMsg2 ); if(Failed()) return
      end if

      ! Blade root to blade root AD
      do iB = 1,y_ED%numBlades
         call Transfer_Point_to_Point(y_ED%BladeRootMotion(iB), u_ADI%AD%rotors(iWT)%BladeRootMotion(iB), y_ED%ED_P_2_AD_P_R(iB), errStat2, errMsg2); if(Failed()) return
      enddo
            
      ! Blade root AD to blade line AD
      if (y_ED%rigidBlades) then
         do iB = 1,y_ED%numBlades
            call Transfer_Point_to_Line2(u_ADI%AD%rotors(iWT)%BladeRootMotion(iB), u_ADI%AD%rotors(iWT)%BladeMotion(iB), y_ED%AD_P_2_AD_L_B(iB), errStat2, errMsg2); if(Failed()) return
         enddo
      else
         print*,'>>> Set_Inputs_For_ADI: TODO Elastic Blades'
         STOP
      endif

      ! Tower motion
      if (y_ED%hasTower) then
         if (u_ADI%AD%rotors(iWT)%TowerMotion%nNodes>0) then
            call Transfer_Point_to_Point(y_ED%TwrPtMesh,  y_ED%TwrPtMeshAD, y_ED%ED_P_2_AD_P_T, errStat2, errMsg2); if(Failed()) return
            call Transfer_Point_to_Line2(y_ED%TwrPtMeshAD, u_ADI%AD%rotors(iWT)%TowerMotion, y_ED%AD_P_2_AD_L_T, errStat2, errMsg2); if(Failed()) return
         endif
      endif
   enddo ! iWT, rotors

contains
   logical function Failed()
      call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'Set_Inputs_For_ADI')
      Failed = errStat >= AbortErrLev
   end function Failed
end subroutine Set_Inputs_For_ADI


end module AeroDyn_Inflow
