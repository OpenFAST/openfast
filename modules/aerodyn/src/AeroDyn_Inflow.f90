!>
!!
!!
module AeroDyn_Inflow
   use NWTC_Library
   use AeroDyn_Inflow_Types
   use AeroDyn_Types
   use AeroDyn, only: AD_Init, AD_ReInit, AD_CalcOutput, AD_UpdateStates, AD_NumWindPoints
   use AeroDyn_IO, only: AD_SetVTKSurface
   use InflowWind, only: InflowWind_Init, InflowWind_CalcOutput

   implicit none

   private

   type(ProgDesc), parameter  :: ADI_Ver = ProgDesc( 'ADI', '', '' )

   public   :: ADI_Init 
   public   :: ADI_ReInit 
   public   :: ADI_End
   public   :: ADI_CalcOutput
   public   :: ADI_UpdateStates

   ! Convenient routines for driver
   public   :: ADI_ADIW_Solve
   public   :: concatOutputHeaders

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
   integer(IntKi)          :: UnEcho         ! Unit number for the echo file
   character(len=1054) :: DirName
   integer :: iW

   ! Initialize variables for this routine
   errStat = ErrID_None
   errMsg  = ""

   ! Initialize the NWTC Subroutine Library
   call NWTC_Init( EchoLibVer=.FALSE. )

   ! Display the module information
   call DispNVD( ADI_Ver )

   ! Set parameters
   p%dt         = interval
   p%storeHHVel = InitInp%storeHHVel
   p%WrVTK      = InitInp%WrVTK

   ! --- Initialize AeroDyn
   if (allocated(InitOut%WriteOutputHdr)) deallocate(InitOut%WriteOutputHdr)
   if (allocated(InitOut%WriteOutputUnt)) deallocate(InitOut%WriteOutputUnt)

   call AD_Init(InitInp%AD, u%AD, p%AD, x%AD, xd%AD, z%AD, OtherState%AD, y%AD, m%AD, Interval, InitOut_AD, errStat2, errMsg2); if (Failed()) return
   InitOut%Ver = InitOut_AD%ver
   ! Add writeoutput units and headers to driver, same for all cases and rotors!
   call concatOutputHeaders(InitOut%WriteOutputHdr, InitOut%WriteOutputUnt, InitOut_AD%rotors(1)%WriteOutputHdr, InitOut_AD%rotors(1)%WriteOutputUnt, errStat2, errMsg2); if(Failed()) return

   ! --- Initialize Inflow Wind 
   call ADI_InitInflowWind(InitInp%RootName, InitInp%IW_InitInp, u%AD, OtherState%AD, m%IW, Interval, InitOut_IW, errStat2, errMsg2); if (Failed()) return
   ! Concatenate AD outputs to IW outputs
   call concatOutputHeaders(InitOut%WriteOutputHdr, InitOut%WriteOutputUnt, InitOut_IW%WriteOutputHdr, InitOut_IW%WriteOutputUnt, errStat2, errMsg2); if(Failed()) return

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
   real(DbKi) :: t
   errStat = ErrID_None
   errMsg  = ""

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
      call ADI_ADIW_Solve(utimes(it), u(it)%AD, OtherState%AD, m%IW%u, m%IW, p%storeHHVel, errStat2, errMsg2)
      u_AD(it) = u(it)%AD
   enddo

   ! Get state variables at next step: INPUT at step nt - 1, OUTPUT at step nt
   call AD_UpdateStates(t, n, u_AD(:), utimes(:), p%AD, x%AD, xd%AD, z%AD, OtherState%AD, m%AD, errStat2, errMsg2); if(Failed()) return

contains

   subroutine CleanUp()
      !call ADI_DestroyConstrState(z_guess, errStat2, errMsg2); if(Failed()) return
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
   ! Local variables
   integer(IntKi)                :: errStat2
   character(errMsgLen)          :: errMsg2
   character(*), parameter       :: RoutineName = 'ADI_CalcOutput'
   integer :: iWT
   errStat = ErrID_None
   errMsg  = ""

   ! --- CalcOutputs for IW (Sets u_AD%rotors(:)%InflowOnBlade, etc,  and m%IW%y)
   ! TODO TODO TODO Uncomment
   !call ADI_ADIW_Solve(t, u%AD, OtherState%AD, m%IW%u, m%IW, p%storeHHVel, errStat2, errMsg2)
   y%IW_WriteOutput(:) = m%IW%y%WriteOutput(:)


   ! Calculate outputs at t
   call AD_CalcOutput(t, u%AD, p%AD, x%AD, xd%AD, z%AD, OtherState%AD, y%AD, m%AD, errStat2, errMsg2); if(Failed()) return

   ! --- Outputs for driver
   ! Hub Height velocity outputs
   if (p%storeHHVel) then
      do iWT = 1, size(p%AD%rotors)
         y%HHVel(1, iWT) = m%IW%y%VelocityUVW(1, iWT)
         y%HHVel(2, iWT) = m%IW%y%VelocityUVW(2, iWT)
         y%HHVel(3, iWT) = m%IW%y%VelocityUVW(3, iWT)
      enddo
   endif
   y%PLExp = m%IW%PLExp

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
   character(len=*),             intent(in   ) :: Root          ! Rootname for input files
   type(ADI_IW_InputData),       intent(in   ) :: i_IW          ! Inflow Wind "pseudo init input" data
   type(AD_InputType),           intent(in   ) :: u_AD          ! AeroDyn data 
   type(AD_OtherStateType),      intent(in   ) :: o_AD          ! AeroDyn data 
   type(ADI_InflowWindData),     intent(inout) :: IW            ! InflowWind data 
   real(DbKi),                   intent(inout) :: dt            ! interval
   type(InflowWind_InitOutputType), intent(out) :: InitOutData  ! Output data from initialization
   integer(IntKi)              , intent(  out) :: errStat       ! Status of error message
   character(*)                , intent(  out) :: errMsg        ! Error message if errStat /= ErrID_None
   ! locals
   real(reKi)                      :: theta(3)
   integer(IntKi)                  :: j, k, nOut_AD, nOut_IW, nOut_Dvr
   integer(IntKi)                  :: iWT
   integer(IntKi)                  :: errStat2      ! local status of error message
   character(errMsgLen)            :: errMsg2       ! local error message if errStat /= ErrID_None
   type(InflowWind_InitInputType)  :: InitInData     ! Input data for initialization
   errStat = ErrID_None
   errMsg  = ''

   ! --- Count number of points required by AeroDyn
   InitInData%NumWindPoints = AD_NumWindPoints(u_AD, o_AD)
   ! Adding Hub windspeed for each turbine
   InitInData%NumWindPoints = InitInData%NumWindPoints + size(u_AD%rotors)

   ! --- Init InflowWind
   if (i_IW%CompInflow==0) then
      ! Fake "InflowWind" init
      allocate(InitOutData%WriteOutputHdr(0))
      allocate(InitOutData%WriteOutputUnt(0))
      allocate(IW%y%WriteOutput(0))
      call AllocAry(IW%u%PositionXYZ, 3, InitInData%NumWindPoints, 'PositionXYZ', errStat2, errMsg2); if (Failed()) return
      call AllocAry(IW%y%VelocityUVW, 3, InitInData%NumWindPoints, 'VelocityUVW', errStat2, errMsg2); if (Failed()) return
      IW%u%PositionXYZ = myNaN
      IW%y%VelocityUVW = myNaN
   else
      ! Module init
      InitInData%InputFileName    = i_IW%InputFile
      InitInData%Linearize        = .false.
      InitInData%UseInputFile     = .true.
      InitInData%RootName         = Root
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
!----------------------------------------------------------------------------------------------------------------------------------
!> Solve for the wind speed at the location necessary for AeroDyn
subroutine ADI_ADIW_Solve(t, u_AD, o_AD, u_IfW, IW, hubHeightFirst, errStat, errMsg)
   real(DbKi),                   intent(in   ) :: t             ! Time of evaluation
   type(AD_InputType),           intent(inout) :: u_AD          ! AeroDyn data 
   type(AD_OtherStateType),      intent(in   ) :: o_AD          ! AeroDyn data 
   type(InflowWind_InputType),   intent(inout) :: u_IfW         ! InflowWind data 
   type(ADI_InflowWindData),     intent(inout) :: IW            ! InflowWind data 
   logical,                      intent(in   ) :: hubHeightFirst ! Hub Height velocity is packed at beginning
   integer(IntKi)              , intent(  out) :: errStat       ! Status of error message
   character(*)                , intent(  out) :: errMsg        ! Error message if errStat /= ErrID_None
   integer(IntKi)       :: errStat2      ! Status of error message
   character(errMsgLen) :: errMsg2       ! Error message if errStat /= ErrID_None
   errStat = ErrID_None
   errMsg  = ''

   ! Set u_ifW%PositionXYZ
   call ADI_Set_IW_Inputs(u_AD, o_AD, u_IfW, hubHeightFirst, errStat2, errMsg2); if(Failed()) return
   ! Compute IW%y%VelocityUVW
   call ADI_CalcOutput_IW(t, u_IfW, IW, errStat2, errMsg2); if(Failed()) return
   ! Set u_AD%..%InflowOnBlade, u_AD%..%InflowOnTower, etc
   call ADI_AD_InputSolve_IfW(u_AD, IW%y, hubHeightFirst, errStat2, errMsg2); if(Failed()) return

contains
   logical function Failed()
      call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'ADI_ADIW_Solve')
      Failed = errStat >= AbortErrLev
   end function Failed
end subroutine ADI_ADIW_Solve
!----------------------------------------------------------------------------------------------------------------------------------
!> Set inputs for inflow wind
!! NOTE: order is important and should match AD_NumWindPoints
!! Similar to FAST_Solver, IfW_InputSolve
subroutine ADI_Set_IW_Inputs(u_AD, o_AD, u_IfW, hubHeightFirst, errStat, errMsg)
   type(AD_InputType),           intent(in   ) :: u_AD          ! AeroDyn data 
   type(AD_OtherStateType),      intent(in   ) :: o_AD          ! AeroDyn data 
   type(InflowWind_InputType),   intent(inout) :: u_IfW         ! InflowWind data 
   logical,                      intent(in   ) :: hubHeightFirst ! Hub Height velocity is packed at beginning
   integer(IntKi)              , intent(  out) :: errStat       ! Status of error message
   character(*)                , intent(  out) :: errMsg        ! Error message if errStat /= ErrID_None
   integer :: K, J, Node, iWT
   errStat = ErrID_None
   errMsg  = ''
   Node=0

   if (hubHeightFirst) then
      ! Hub Height point for each turbine
      do iWT=1,size(u_AD%rotors)
         Node = Node + 1
         u_IfW%PositionXYZ(:,Node) = u_AD%rotors(iWT)%hubMotion%Position(:,1) + u_AD%rotors(iWT)%hubMotion%TranslationDisp(:,1)
      enddo
   endif

   do iWT=1,size(u_AD%rotors)
      ! Blade
      do K = 1,size(u_AD%rotors(iWT)%BladeMotion)
         do J = 1,u_AD%rotors(iWT)%BladeMotion(k)%Nnodes
            Node = Node + 1
            u_IfW%PositionXYZ(:,Node) = u_AD%rotors(iWT)%BladeMotion(k)%TranslationDisp(:,j) + u_AD%rotors(iWT)%BladeMotion(k)%Position(:,j)
         end do !J = 1,p%BldNodes ! Loop through the blade nodes / elements
      end do !K = 1,p%NumBl         
      ! Tower
      do J=1,u_AD%rotors(iWT)%TowerMotion%nnodes
         Node = Node + 1
         u_IfW%PositionXYZ(:,Node) = u_AD%rotors(iWT)%TowerMotion%TranslationDisp(:,J) + u_AD%rotors(iWT)%TowerMotion%Position(:,J)
      end do      
      ! Nacelle
      if (u_AD%rotors(iWT)%NacelleMotion%Committed) then
         Node = Node + 1
         u_IfW%PositionXYZ(:,Node) = u_AD%rotors(iWT)%NacelleMotion%TranslationDisp(:,1) + u_AD%rotors(iWT)%NacelleMotion%Position(:,1)
      end if
      ! Hub

   enddo ! iWT
   ! vortex points from FVW in AD15
   if (allocated(o_AD%WakeLocationPoints)) then
      do J=1,size(o_AD%WakeLocationPoints,dim=2)
         Node = Node + 1
         u_IfW%PositionXYZ(:,Node) = o_AD%WakeLocationPoints(:,J)
         ! rewrite the history of this so that extrapolation doesn't make a mess of things
!          do k=2,size(IW%u)
!             if (allocated(IW%u(k)%PositionXYZ))   IW%u(k)%PositionXYZ(:,Node) = IW%u(1)%PositionXYZ(:,Node)
!          end do
      enddo !j, wake points
   end if
end subroutine ADI_Set_IW_Inputs
!----------------------------------------------------------------------------------------------------------------------------------
!> Calculate Wind at desired points
!! NOTE: order is important and should match AD_NumWindPoints
!! Similar to FAST_Solver, IfW_InputSolve
subroutine ADI_CalcOutput_IW(t, u_IfW, IW, errStat, errMsg)
   real(DbKi),                   intent(in   ) :: t             ! Time of evaluation
   type(InflowWind_InputType),   intent(inout) :: u_IfW         ! InflowWind data 
   type(ADI_InflowWindData),     intent(inout) :: IW            ! InflowWind data 
   integer(IntKi)              , intent(  out) :: errStat       ! Status of error message
   character(*)                , intent(  out) :: errMsg        ! Error message if errStat /= ErrID_None
   integer              :: j
   real(ReKi)           :: z
   integer(IntKi)       :: errStat2      ! Status of error message
   character(errMsgLen) :: errMsg2       ! Error message if errStat /= ErrID_None
   errStat = ErrID_None
   errMsg  = ''
   if (IW%CompInflow==1) then
      call InflowWind_CalcOutput(t, u_IfW, IW%p, IW%x, IW%xd, IW%z, IW%OtherSt, IW%y, IW%m, errStat2, errMsg2)
      call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'ADI_CalcOutput_IW') 
   else
      do j=1,size(u_IfW%PositionXYZ,2)
         z = u_IfW%PositionXYZ(3,j)
         IW%y%VelocityUVW(1,j) = IW%HWindSpeed*(z/IW%RefHt)**IW%PLExp
         IW%y%VelocityUVW(2,j) = 0.0_ReKi !V
         IW%y%VelocityUVW(3,j) = 0.0_ReKi !W      
      end do 
   endif
end subroutine ADI_CalcOutput_IW
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine sets the wind claculated by InflowWind to the AeroDyn arrays
!! See similar routine in FAST_Solver
subroutine ADI_AD_InputSolve_IfW(u_AD, y_IfW, hubHeightFirst, errStat, errMsg)
   ! Passed variables
   TYPE(AD_InputType),          INTENT(INOUT)   :: u_AD        !< The inputs to AeroDyn
   TYPE(InflowWind_OutputType), INTENT(IN)      :: y_IfW       !< The outputs from InflowWind
   logical,                      intent(in   ) :: hubHeightFirst ! Hub Height velocity is packed at beginning
   INTEGER(IntKi)                               :: errStat     !< Error status of the operation
   CHARACTER(*)                                 :: errMsg      !< Error message if errStat /= ErrID_None
   ! Local variables:
   INTEGER(IntKi)                               :: J           ! Loops through nodes / elements.
   INTEGER(IntKi)                               :: K           ! Loops through blades.
   INTEGER(IntKi)                               :: NumBl
   INTEGER(IntKi)                               :: NNodes
   INTEGER(IntKi)                               :: node
   INTEGER(IntKi)                               :: iWT
   errStat = ErrID_None
   errMsg  = ""
   node = 1
   ! Order important!
   if (hubHeightFirst) then
      do iWT=1,size(u_AD%rotors)
         node = node + 1 ! Hub velocities for each rotor
      enddo
   endif

   do iWT=1,size(u_AD%rotors)
      NumBl  = size(u_AD%rotors(iWT)%InflowOnBlade,3)
      Nnodes = size(u_AD%rotors(iWT)%InflowOnBlade,2)
      ! Blades
      do k=1,NumBl
         do j=1,Nnodes
            u_AD%rotors(iWT)%InflowOnBlade(:,j,k) = y_IfW%VelocityUVW(:,node)
            node = node + 1
         end do
      end do
      ! Tower
      if ( allocated(u_AD%rotors(iWT)%InflowOnTower) ) then
         Nnodes = size(u_AD%rotors(iWT)%InflowOnTower,2)
         do j=1,Nnodes
            u_AD%rotors(iWT)%InflowOnTower(:,j) = y_IfW%VelocityUVW(:,node)
            node = node + 1
         end do      
      end if
      ! Nacelle
      if (u_AD%rotors(iWT)%NacelleMotion%NNodes > 0) then
         u_AD%rotors(iWT)%InflowOnNacelle(:) = y_IfW%VelocityUVW(:,node)
         node = node + 1
      else
         u_AD%rotors(iWT)%InflowOnNacelle = 0.0_ReKi
      end if
      ! Hub 
!      if (u_AD%HubMotion%NNodes > 0) then
!         u_AD%InflowOnHub(:) = y_IfW%VelocityUVW(:,node)
!         node = node + 1
!      else
!         u_AD%InflowOnHub = 0.0_ReKi
!      end if
   enddo ! rotors
   ! OLAF points
   if ( allocated(u_AD%InflowWakeVel) ) then
      Nnodes = size(u_AD%InflowWakeVel,DIM=2)
      do j=1,Nnodes
         u_AD%InflowWakeVel(:,j) = y_IfW%VelocityUVW(:,node)
         node = node + 1
      end do !j, wake points
   end if
end subroutine ADI_AD_InputSolve_IfW

end module AeroDyn_Inflow
