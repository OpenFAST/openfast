!**********************************************************************************************************************************
!> LinDyn, module for a second order linear dynamical system with mass, stiffness and damping matrix
!!   
!!  The state is q = [x; xdot], of shape nq = 2*nx
!!  The input is F_ext of shape nx
!!  The equation of motion is:
!!   
!!     qdot = [xdot ] = [  0              I    ] [ x   ]  + [ 0    ] F_ext
!!            [xddot]   [-M^{-1} K   -M^{-1} C ] [ xdot]  + [M^{-1}]
!!  
!! ..................................................................................................................................
!! ## Licensing
!! Copyright (C) 2012-2013, 2015-2016  National Renewable Energy Laboratory
!!
!!    This file is part of LinDyn.
!!
!! Licensed under the Apache License, Version 2.0 (the "License");
!! you may not use this file except in compliance with the License.
!! You may obtain a copy of the License at
!!
!!     http://www.apache.org/licenses/LICENSE-2.0
!!
!! Unless required by applicable law or agreed to in writing, software
!! distributed under the License is distributed on an "AS IS" BASIS,
!! WITHout WARRANTIES OR CONDITIONS OF ANY KinD, either express or implied.
!! See the License for the specific language governing permissions and
!! limitations under the License.
!**********************************************************************************************************************************
module LinDyn

   use LinDyn_Types
   use NWTC_Library
   USE NWTC_LAPACK
! 
   implicit none

   type(ProgDesc), parameter  :: LD_Ver = ProgDesc( 'LinDyn', '', '' )

   private

   public :: LD_Init                           !  Initialization routine
   public :: LD_InitInputData                  ! Set default values and allocations for init
   public :: LD_End                            !  Ending routine (includes clean up)
   public :: LD_UpdateStates                   !  Loose coupling routine for solving for constraint states, integrating continuous states, and updating discrete states
   public :: LD_CalcOutput                     !  Routine for computing outputs
   public :: LD_CalcContStateDeriv             !  Tight coupling routine for computing derivatives of continuous states
   public :: LD_JacobianPInput                 !  Jacobians of (y, x, xd, z) with respect to the inputs (u)
   public :: LD_JacobianPContState             !  Jacobians of (y, x, xd, z) with respect to the continuous (x)
   public :: LD_GetOP                          !  Routine to get the operating-point values for linearization (from data structures to arrays)
! 
contains

subroutine LD_Init(InitInp, u, p, x, xd, z, OtherState, y, m, InitOut, errStat, errMsg)
   type(LD_InitInputType),       intent(in   ) :: InitInp     !< Input data for initialization routine
   type(LD_InputType),           intent(out)   :: u           !< An initial guess for the input; input mesh must be defined
   type(LD_ParameterType),       intent(out)   :: p           !< Parameters
   type(LD_ContinuousStateType), intent(out)   :: x           !< Initial continuous states
   type(LD_DiscreteStateType),   intent(out)   :: xd          !< Initial discrete states
   type(LD_ConstraintStateType), intent(out)   :: z           !< Initial guess of the constraint states
   type(LD_OtherStateType),      intent(out)   :: OtherState  !< Initial other states (logical, etc)
   type(LD_OutputType),          intent(out)   :: y           !< Initial system outputs (outputs are not calculated;
   type(LD_MiscVarType),         intent(out)   :: m           !< Misc variables for optimization (not copied in glue code)
   type(LD_InitOutputType),      intent(out)   :: InitOut     !< Output for initialization routine
   integer(IntKi),               intent(out)   :: errStat     !< Error status of the operation
   character(*),                 intent(out)   :: errMsg      !< Error message if errStat /= ErrID_None
   integer(IntKi)  :: errStat2    ! Status of error message
   character(1024) :: errMsg2     ! Error message if ErrStat /= ErrID_None
   ! Misc Init
   errStat = ErrID_None
   errMsg  = ""
   call NWTC_Init( )      ! Initialize the NWTC Subroutine Library
   call DispNVD( LD_Ver ) ! Display the module information

   ! --- Setting Params from InitInp
   p%nx = size(InitInp%MM,1)
   p%nq = 2*p%nx
   call AllocAry(p%MM        , p%nx, p%nx, 'MM', errStat2, errMsg2); if(Failed()) return
   call AllocAry(p%CC        , p%nx, p%nx, 'CC', errStat2, errMsg2); if(Failed()) return
   call AllocAry(p%KK        , p%nx, p%nx, 'KK', errStat2, errMsg2); if(Failed()) return
   call AllocAry(p%activeDOFs, p%nx      , 'activeDOFs', errStat2, errMsg2); if(Failed()) return
   p%dt         = InitInp%dt
   p%IntMethod  = InitInp%IntMethod
   p%MM         = InitInp%MM
   p%CC         = InitInp%CC
   p%KK         = InitInp%KK
   p%activeDOFs = InitInp%activeDOFs
   ! Prescribed motion
   if (len_trim(InitInp%PrescribedMotionFile)>0) then
      if( count(p%activeDOFs)/=0) then
         errStat2 = errID_Fatal
         errMsg2  = 'Currently, prescribed motion is only allowed if all degrees of freedom are turned off'
         if(Failed()) return
      endif
      call WrScr('    Using prescribed motion.')
      call ReadDelimFile(InitInp%PrescribedMotionFile, (p%nx*3+1), p%PrescribedValues, errStat2, errMsg2); if(Failed()) return
   else
      if (allocated(p%PrescribedValues)) deallocate(p%PrescribedValues)
   endif
   call StateMatrices(p%MM, p%CC, p%KK, p%AA, p%BB, errStat2, errMsg2); if(Failed()) return

   ! --- Misc
   call allocAry(m%qPrescribed, 3*p%nx, 'qPrescribed', errStat2, errMsg2); if(Failed()) return
   m%qPrescribed = 0.0_ReKi ! NOTE: will be updated by LD_SetInitialConditions

   ! --- Allocate States
   call AllocAry( x%q    , p%nq, 'DOFs', errStat, errMsg); if(Failed()) return
   call LD_SetInitialConditions(x, InitInp%x0, InitInp%xd0, p, OtherState, m, errStat, errMsg); if(Failed()) return
   if ( ( p%IntMethod .eq. 2) .OR. ( p%IntMethod .eq. 3)) then !Multi-step methods
       allocate( OtherState%xdot(4), STAT=errStat2); errMsg2='Error allocating OtherState%xdot'
       if(Failed()) return
   endif

   ! --- Guess inputs
   call AllocAry(u%Fext, p%nx, 'Fext', errStat2, errMsg2); if(Failed()) return
   u%Fext=0.0_ReKi

   ! --- Outputs &  Write Outputs
   call Init_Outputs(p, m, y, InitInp, InitOut, errStat, errMsg); if(Failed()) return
   InitOut%Ver = LD_Ver

   ! --- Linearization
   if (InitInp%Linearize) then
      call Init_Lin(p, InitOut, errStat, errMsg); if(Failed()) return
   endif
! 
!    ! --- Summary file 
!    if (InputFileData%SumPrint) then
!    TODO use yaml  
!    print*,''
!    print*,'M',p%MM(1,:)
!    print*,'M',p%MM(2,:)
!    print*,'M',p%MM(3,:)
!    print*,''
!    print*,'C',p%CC(1,:)
!    print*,'C',p%CC(2,:)
!    print*,'C',p%CC(3,:)
!    print*,''
!    print*,'K',p%KK(1,:)
!    print*,'K',p%KK(2,:)
!    print*,'K',p%KK(3,:)
!    print*,''
! 
!    print*,''
!    print*,'A',p%AA(1,:)
!    print*,'A',p%AA(2,:)
!    print*,'A',p%AA(3,:)
!    print*,'A',p%AA(4,:)
!    print*,'A',p%AA(5,:)
!    print*,'A',p%AA(6,:)
!    print*,''
!    print*,'B',p%BB(1,:)
!    print*,'B',p%BB(2,:)
!    print*,'B',p%BB(3,:)
!    print*,'B',p%BB(4,:)
!    print*,'B',p%BB(5,:)
!    print*,'B',p%BB(6,:)
!    endif
! 
contains
   logical function Failed()
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'LD_Init' )
      Failed = ErrStat >= AbortErrLev
      if (Failed) call CleanUp()
   end function Failed
   subroutine CleanUp()
   end subroutine CleanUp
end subroutine LD_Init
!----------------------------------------------------------------------------------------------------------------------------------
subroutine LD_SetInitialConditions(x, x0, xd0, p, OtherState, m, errStat, errMsg)
   type(LD_ContinuousStateType), intent(inout) :: x          !< Initial continuous states
   real(ReKi),                   intent(in)    :: x0(:)      !< Values of the positions at t=0
   real(ReKi),                   intent(in)    :: xd0(:)     !< Velocity values at t=0
   type(LD_ParameterType),       intent(in   ) :: p          !< Parameters
   type(LD_OtherStateType),      intent(inout) :: OtherState !< Other states
   type(LD_MiscVarType),         intent(inout) :: m          !< Misc variables for optimization (not copied in glue code)
   integer(IntKi),               intent(out)   :: errStat    !< Error status of the operation
   character(*),                 intent(out)   :: errMsg     !< Error message if errStat /= ErrID_None
   integer :: nx
   nx = int(size(x%q)/2)
   errStat = ErrID_Fatal
   if (size(x0)/=size(xd0)) then
      errMsg ='Shape of x0 and xd0 should match when setting intial conditions'; return
   endif
   if (size(x0)/=nx) then
      errMsg ='Shape of x0 should match nx when setting intial conditions'; return
   endif
   errMsg  = ''
   errStat = ErrID_None

   if (allocated(p%PrescribedValues)) then
      call interpTimeValue(p%PrescribedValues, 0.0_DbKi, OtherState%iMotionInterpLast, m%qPrescribed(:))
      ! TODO the code below will need to be updated if a subset of the DOFs are active
      x%q(1:p%nq) = m%qPrescribed(1:p%nq)
   else
      x%q(   1:nx)   = x0
      x%q(nx+1:2*nx) = xd0
   endif
end subroutine LD_SetInitialConditions
!----------------------------------------------------------------------------------------------------------------------------------
!> Allocate init input data for module based on number of degrees of freedom
subroutine LD_InitInputData(nx, InitInp, errStat, errMsg)
   integer(IntKi),               intent(in ) :: nx        !< Number of degrees of freedom
   type(LD_InitInputType),       intent(out) :: InitInp     !< Input data for initialization routine
   integer(IntKi),               intent(out) :: errStat   !< Error status of the operation
   character(*),                 intent(out) :: errMsg    !< Error message if errStat /= ErrID_None
   integer(IntKi)  :: iDOF
   integer(IntKi)  :: errStat2    ! Status of error message
   character(1024) :: errMsg2     ! Error message if ErrStat /= ErrID_None
   ! Initialize errStat
   errStat   = ErrID_None           ! no error has occurred
   errMsg    = ""
   call AllocAry(InitInp%MM        , nx, nx, 'MM' , errStat2, errMsg2); if(Failed()) return
   call AllocAry(InitInp%CC        , nx, nx, 'CC' , errStat2, errMsg2); if(Failed()) return
   call AllocAry(InitInp%KK        , nx, nx, 'KK' , errStat2, errMsg2); if(Failed()) return
   call AllocAry(InitInp%x0        , nx    , 'x0' , errStat2, errMsg2); if(Failed()) return
   call AllocAry(InitInp%xd0       , nx    , 'xd0', errStat2, errMsg2); if(Failed()) return
   call AllocAry(InitInp%activeDOFs, nx    , 'activeDOFs', errStat2, errMsg2); if(Failed()) return
   call AllocAry(InitInp%DOFsNames , nx    , 'DOFsNames' , errStat2, errMsg2); if(Failed()) return
   call AllocAry(InitInp%DOFsUnits , nx    , 'DOFsUnits' , errStat2, errMsg2); if(Failed()) return
   InitInp%MM         = 0.0_ReKi
   InitInp%CC         = 0.0_ReKi
   InitInp%KK         = 0.0_ReKi
   InitInp%x0         = 0.0_ReKi
   InitInp%xd0        = 0.0_ReKi
   InitInp%activeDOFs = .True.
   ! Default DOFs Names and Units
   do iDOF=1,nx
      InitInp%DOFsNames(iDOF)='x'//trim(num2lstr(iDOF))
      InitInp%DOFsUnits(iDOF)='-'
   enddo

contains
   logical function Failed()
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'LD_Init' )
      Failed = ErrStat >= AbortErrLev
   end function Failed
end subroutine LD_InitInputData
!----------------------------------------------------------------------------------------------------------------------------------
!> Compute A and B state matrices for a linear mechanical system
!! NOTE: Generic function (no derived types), keep it that way
!! A = [   0        I       ]     B = [ 0      ]
!!     [-M^{-1}K   -M^{-1}C ]       = [ M^{-1} ]
subroutine StateMatrices(MM, CC, KK, AA, BB, errStat, errMsg)
   real(ReKi),              intent(in ) :: MM(:,:)
   real(ReKi),              intent(in ) :: CC(:,:)
   real(ReKi),              intent(in ) :: KK(:,:)
   real(ReKi), allocatable, intent(out) :: AA(:,:)
   real(ReKi), allocatable, intent(out) :: BB(:,:)
   integer(IntKi),          intent(out) :: errStat   !< Error status of the operation
   character(*),            intent(out) :: errMsg    !< Error message if errStat /= ErrID_None
   integer(IntKi)                          :: errStat2    ! Status of error message
   character(1024)                         :: errMsg2     ! Error message if ErrStat /= ErrID_None
   integer                                 :: nx, nq, i
   
   real(ReKi), dimension(:,:), allocatable :: MLU     ! LU factorization of M matrix
   real(ReKi), dimension(:,:), allocatable :: MinvX   ! Tmp array to store either: M^{-1} C, M^{-1} K , or M^{-1}
   real(ReKi), dimension(:) , allocatable  :: WORK    ! LAPACK variable
   integer, allocatable                    :: IPIV(:) ! LAPACK variable
   integer                                 :: LWORK   ! LAPACK variable
   ! Initialize errStat
   errStat = ErrID_None
   errMsg  = ""

   ! --- Init A and B matrix
   nx = size(MM,1)
   nq = 2*nx
   call AllocAry(AA, nq, nq, 'AA', errStat2, errMsg2); if(Failed()) return
   call AllocAry(BB, nq, nx, 'BB', errStat2, errMsg2); if(Failed()) return
   AA(:,:) = 0.0_ReKi
   BB(:,:) = 0.0_ReKi
   do i=1,nx ; AA(i,i+nx)=1; enddo ! Identity matrix for upper right block

   ! --- Compute misc inverse of M and put in A and B matrices
   call AllocAry(IPIV        , nx    , 'IPIV' , errStat2, errMsg2); if(Failed()) return
   call AllocAry(MinvX       , nx, nx, 'MinvX', errStat2, errMsg2); if(Failed()) return
   call AllocAry(MLU         , nx, nx, 'MLU'  , errStat2, errMsg2); if(Failed()) return

   ! LU Factorization of M
   MLU = MM ! temp copy
   call LAPACK_getrf(nx, nx, MLU, IPIV, errStat2, errMsg2); if(Failed()) return

   ! M^-1 C
   MinvX = CC
   call LAPACK_getrs('n', nx, MLU, IPIV, MinvX, errStat2, errMsg2); if(Failed()) return
   AA(nx+1:nq,nx+1:nq) = -MinvX

   ! M^-1 K
   MinvX = KK
   call LAPACK_getrs('n', nx, MLU, IPIV, MinvX, errStat2, errMsg2); if(Failed()) return
   AA(nx+1:nq, 1:nx) = -MinvX

   ! Inverse of M
   MinvX = MLU
   LWORK=nx*nx ! Somehow LWORK = -1 does not work
   call AllocAry(WORK, LWORk, 'WORK', errStat2, errMsg2); if(Failed()) return
   call LAPACK_getri(nx, MinvX, IPIV, WORK, LWORK, errStat2, errMsg2); if(Failed()) return
   BB(nx+1:nq, : ) = MinvX

   call CleanUp()
   
contains
   logical function Failed()
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'StateMatrices' )
      Failed = ErrStat >= AbortErrLev
      if (Failed) call CleanUp()
   end function Failed
   subroutine CleanUp()
      if (allocated(MLU))   deallocate(MLU)
      if (allocated(IPIV))  deallocate(IPIV)
      if (allocated(WORK )) deallocate(WORK)
      if (allocated(MinvX)) deallocate(MinvX)
   end subroutine CleanUp
end subroutine StateMatrices
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is called at the end of the simulation.
subroutine LD_End( u, p, x, xd, z, OtherState, y, m, errStat, errMsg )
   type(LD_InputType),           intent(inout) :: u          !< System inputs
   type(LD_ParameterType),       intent(inout) :: p          !< Parameters
   type(LD_ContinuousStateType), intent(inout) :: x          !< Continuous states
   type(LD_DiscreteStateType),   intent(inout) :: xd         !< Discrete states
   type(LD_ConstraintStateType), intent(inout) :: z          !< Constraint states
   type(LD_OtherStateType),      intent(inout) :: OtherState !< Other states
   type(LD_OutputType),          intent(inout) :: y          !< System outputs
   type(LD_MiscVarType),         intent(inout) :: m          !< Misc variables for optimization (not copied in glue code)
   integer(IntKi),               intent(out)   :: errStat    !< Error status of the operation
   character(*),                 intent(out)   :: errMsg     !< Error message if errStat /= ErrID_None
   ! Initialize errStat
   errStat   = ErrID_None           ! no error has occurred
   errMsg    = ""
   call LD_DestroyInput      (u         ,errStat,errMsg)
   call LD_DestroyParam      (p         ,errStat,errMsg)
   call LD_DestroyContState  (x         ,errStat,errMsg)
   call LD_DestroyDiscState  (xd        ,errStat,errMsg)
   call LD_DestroyConstrState(z         ,errStat,errMsg)
   call LD_DestroyOtherState (OtherState,errStat,errMsg)
   call LD_DestroyOutput     (y         ,errStat,errMsg)
   call LD_DestroyMisc       (m         ,errStat,errMsg)
end subroutine LD_End
!----------------------------------------------------------------------------------------------------------------------------------
!> Fourth-order Adams-Bashforth Method (RK4) for numerically integration (see ElastoDyn.f9)
subroutine LD_AB4( t, n, u, utimes, p, x, xd, z, OtherState, m, errStat, errMsg )
   real(DbKi),                   intent(in   ) :: t          !< Current simulation time in seconds
   integer(IntKi),               intent(in   ) :: n          !< time step number
   type(LD_InputType),           intent(inout) :: u(:)       !< Inputs at t
   real(DbKi),                   intent(in   ) :: utimes(:)  !< times of input
   type(LD_ParameterType),       intent(in   ) :: p          !< Parameters
   type(LD_ContinuousStateType), intent(inout) :: x          !< Continuous states at t on input at t + dt on output
   type(LD_DiscreteStateType),   intent(in   ) :: xd         !< Discrete states at t
   type(LD_ConstraintStateType), intent(in   ) :: z          !< Constraint states at t (possibly a guess)
   type(LD_OtherStateType),      intent(inout) :: OtherState !< Other states at t on input at t + dt on output
   type(LD_MiscVarType),         intent(inout) :: m          !< Misc/optimization variables
   integer(IntKi),               intent(out)   :: errStat    !< Error status of the operation
   character(*),                 intent(out)   :: errMsg     !< Error message if errStat /= ErrID_None
   ! local variables
   type(LD_ContinuousStateType) :: xdot       ! Continuous state derivs at t
   type(LD_InputType)           :: u_interp
   ! Initialize errStat
   errStat = ErrID_None
   errMsg  = "" 
   
   ! need xdot at t
   call LD_CopyInput(u(1), u_interp, MESH_NEWCOPY, errStat, errMsg  )  ! we need to allocate input arrays/meshes before calling ExtrapInterp...
   call LD_Input_ExtrapInterp(u, utimes, u_interp, t, errStat, errMsg)
   call LD_CalcContStateDeriv( t, u_interp, p, x, xd, z, OtherState, m, xdot, errStat, errMsg ) ! initializes xdot
   call LD_DestroyInput( u_interp, errStat, errMsg)   ! we don't need this local copy anymore
   if (n .le. 2) then
      OtherState%n = n
      call LD_CopyContState(xdot, OtherState%xdot(3-n), MESH_UPDATECOPY, errStat, errMsg )
      call LD_RK4(t, n, u, utimes, p, x, xd, z, OtherState, m, errStat, errMsg )
   else
      if (OtherState%n .lt. n) then
         OtherState%n = n
         call LD_CopyContState(OtherState%xdot(3), OtherState%xdot(4), MESH_UPDATECOPY, errStat, errMsg )
         call LD_CopyContState(OtherState%xdot(2), OtherState%xdot(3), MESH_UPDATECOPY, errStat, errMsg )
         call LD_CopyContState(OtherState%xdot(1), OtherState%xdot(2), MESH_UPDATECOPY, errStat, errMsg )
      elseif (OtherState%n .gt. n) then
         errStat = ErrID_Fatal
         errMsg = ' Backing up in time is not supported with a multistep method '
         return
      endif
      call LD_CopyContState( xdot, OtherState%xdot ( 1 ), MESH_UPDATECOPY, errStat, errMsg )
      !OtherState%xdot ( 1 )     = xdot  ! make sure this is most up to date
      x%q = x%q + (p%dt / 24._ReKi) * (55._ReKi*OtherState%xdot(1)%q - 59._ReKi*OtherState%xdot(2)%q + 37._ReKi*OtherState%xdot(3)%q - 9._ReKi * OtherState%xdot(4)%q)
   endif
   call LD_DestroyContState(xdot, errStat, errMsg)
   call LD_DestroyInput(u_interp, errStat, errMsg)
end subroutine LD_AB4
!----------------------------------------------------------------------------------------------------------------------------------
!> Fourth-order Adams-Bashforth-Moulton Method (RK4) for numerically integrating (see ElastoDyn.f90)
subroutine LD_ABM4( t, n, u, utimes, p, x, xd, z, OtherState, m, errStat, errMsg )
   real(DbKi),                   intent(in   ) :: t          !< Current simulation time in seconds
   integer(IntKi),               intent(in   ) :: n          !< time step number
   type(LD_InputType),           intent(inout) :: u(:)       !< Inputs at t
   real(DbKi),                   intent(in   ) :: utimes(:)  !< times of input
   type(LD_ParameterType),       intent(in   ) :: p          !< Parameters
   type(LD_ContinuousStateType), intent(inout) :: x          !< Continuous states at t on input at t + dt on output ! TODO TODO TODO in
   type(LD_DiscreteStateType),   intent(in   ) :: xd         !< Discrete states at t
   type(LD_ConstraintStateType), intent(in   ) :: z          !< Constraint states at t (possibly a guess)
   type(LD_OtherStateType),      intent(inout) :: OtherState !< Other states at t on input at t + dt on output
   type(LD_MiscVarType),         intent(inout) :: m          !< Misc/optimization variables
   integer(IntKi),               intent(out)   :: errStat    !< Error status of the operation
   character(*),                 intent(out)   :: errMsg     !< Error message if errStat /= ErrID_None
   ! local variables
   type(LD_InputType)            :: u_interp        ! Continuous states at t
   type(LD_ContinuousStateType)  :: x_pred          ! Continuous states at t
   type(LD_ContinuousStateType)  :: xdot_pred       ! Continuous states at t
   ! Initialize errStat
   errStat = ErrID_None
   errMsg  = "" 
   call LD_CopyContState(x, x_pred, MESH_NEWCOPY, errStat, errMsg) !initialize x_pred      
   call LD_AB4( t, n, u, utimes, p, x_pred, xd, z, OtherState, m, errStat, errMsg )
   if (n .gt. 2) then
      call LD_CopyInput( u(1), u_interp, MESH_NEWCOPY, errStat, errMsg) ! make copy so that arrays/meshes get initialized/allocated for ExtrapInterp
      call LD_Input_ExtrapInterp(u, utimes, u_interp, t + p%dt, errStat, errMsg)
      call LD_CalcContStateDeriv(t + p%dt, u_interp, p, x_pred, xd, z, OtherState, m, xdot_pred, errStat, errMsg ) ! initializes xdot_pred
      call LD_DestroyInput( u_interp, errStat, errMsg) ! local copy no longer needed
   
      x%q    = x%q    + (p%dt / 24.) * ( 9. * xdot_pred%q +  19. * OtherState%xdot(1)%q - 5. * OtherState%xdot(2)%q  + 1. * OtherState%xdot(3)%q )
      call LD_DestroyContState( xdot_pred, errStat, errMsg) ! local copy no longer needed
   else
      x%q    = x_pred%q
   endif
   call LD_DestroyContState( x_pred, errStat, errMsg) ! local copy no longer needed
end subroutine LD_ABM4
!----------------------------------------------------------------------------------------------------------------------------------
!> Fourth-order Runge-Kutta Method (RK4) for numerically integration (see ElastoDyn.f90)
subroutine LD_RK4( t, n, u, utimes, p, x, xd, z, OtherState, m, errStat, errMsg )
   real(DbKi),                   intent(in   ) :: t          !< Current simulation time in seconds
   integer(IntKi),               intent(in   ) :: n          !< time step number
   type(LD_InputType),           intent(inout) :: u(:)       !< Inputs at t
   real(DbKi),                   intent(in   ) :: utimes(:)  !< times of input
   type(LD_ParameterType),       intent(in   ) :: p          !< Parameters
   type(LD_ContinuousStateType), intent(inout) :: x          !< Continuous states at t on input at t + dt on output
   type(LD_DiscreteStateType),   intent(in   ) :: xd         !< Discrete states at t
   type(LD_ConstraintStateType), intent(in   ) :: z          !< Constraint states at t (possibly a guess)
   type(LD_OtherStateType),      intent(inout) :: OtherState !< Other states at t on input at t + dt on output
   type(LD_MiscVarType),         intent(inout) :: m          !< Misc/optimization variables
   integer(IntKi),               intent(out)   :: errStat    !< Error status of the operation
   character(*),                 intent(out)   :: errMsg     !< Error message if errStat /= ErrID_None
   ! local variables
   type(LD_ContinuousStateType)                 :: xdot        ! time derivatives of continuous states      
   type(LD_ContinuousStateType)                 :: k1          ! RK4 constant; see above
   type(LD_ContinuousStateType)                 :: k2          ! RK4 constant; see above 
   type(LD_ContinuousStateType)                 :: k3          ! RK4 constant; see above 
   type(LD_ContinuousStateType)                 :: k4          ! RK4 constant; see above 
   type(LD_ContinuousStateType)                 :: x_tmp       ! Holds temporary modification to x
   type(LD_InputType)                           :: u_interp    ! interpolated value of inputs 
   ! Initialize errStat
   errStat = ErrID_None
   errMsg  = "" 
   
   ! Initialize interim vars
   call LD_CopyContState( x, k1,       MESH_NEWCOPY, errStat, errMsg )
   call LD_CopyContState( x, k2,       MESH_NEWCOPY, errStat, errMsg )
   call LD_CopyContState( x, k3,       MESH_NEWCOPY, errStat, errMsg )
   call LD_CopyContState( x, k4,       MESH_NEWCOPY, errStat, errMsg )
   call LD_CopyContState( x, x_tmp,    MESH_NEWCOPY, errStat, errMsg )
   
   ! interpolate u to find u_interp = u(t)
   call LD_CopyInput(u(1), u_interp, MESH_NEWCOPY, errStat, errMsg  )  ! we need to allocate input arrays/meshes before calling ExtrapInterp...     
   call LD_Input_ExtrapInterp( u, utimes, u_interp, t, errStat, errMsg )
   
   ! find xdot at t
   call LD_CalcContStateDeriv( t, u_interp, p, x, xd, z, OtherState, m, xdot, errStat, errMsg ) !initializes xdot
   
   k1%q    = p%dt * xdot%q
   x_tmp%q = x%q    + 0.5_ReKi * k1%q
   
   ! interpolate u to find u_interp = u(t + dt/2)
   call LD_Input_ExtrapInterp(u, utimes, u_interp, t+0.5_ReKi*p%dt, errStat, errMsg)
   
   ! find xdot at t + dt/2
   call LD_CalcContStateDeriv( t + 0.5_ReKi*p%dt, u_interp, p, x_tmp, xd, z, OtherState, m, xdot, errStat, errMsg )
   
   k2%q    = p%dt * xdot%q
   x_tmp%q = x%q    + 0.5_ReKi * k2%q
   
   ! find xdot at t + dt/2
   call LD_CalcContStateDeriv( t + 0.5_ReKi*p%dt, u_interp, p, x_tmp, xd, z, OtherState, m, xdot, errStat, errMsg )
   
   k3%q    = p%dt * xdot%q
   x_tmp%q = x%q    + k3%q
   
   ! interpolate u to find u_interp = u(t + dt)
   call LD_Input_ExtrapInterp(u, utimes, u_interp, t + p%dt, errStat, errMsg)
   
   ! find xdot at t + dt
   call LD_CalcContStateDeriv( t + p%dt, u_interp, p, x_tmp, xd, z, OtherState, m, xdot, errStat, errMsg )
   k4%q   = p%dt * xdot%q
   x%q    = x%q    +  ( k1%q    + 2._ReKi * k2%q    + 2._ReKi * k3%q    + k4%q    ) / 6._ReKi
   call CleanUp()
contains      
   subroutine CleanUp()
      integer(IntKi)             :: errStat3    ! The error identifier (errStat)
      character(1024)            :: errMsg3     ! The error message (errMsg)
      call LD_DestroyContState( xdot,     errStat3, errMsg3 )
      call LD_DestroyContState( k1,       errStat3, errMsg3 )
      call LD_DestroyContState( k2,       errStat3, errMsg3 )
      call LD_DestroyContState( k3,       errStat3, errMsg3 )
      call LD_DestroyContState( k4,       errStat3, errMsg3 )
      call LD_DestroyContState( x_tmp,    errStat3, errMsg3 )
      call LD_DestroyInput(     u_interp, errStat3, errMsg3 )
   end subroutine CleanUp            
end subroutine LD_RK4
!----------------------------------------------------------------------------------------------------------------------------------
!> Loose coupling routine for solving states at t+dt
subroutine LD_UpdateStates( t, n, Inputs, InputTimes, p, x, xd, z, OtherState, m, errStat, errMsg )
   real(DbKi),                   intent(in   ) :: t             !< Current simulation time in seconds
   integer(IntKi),               intent(in   ) :: n             !< Current step of the simulation: t = n*dt
   type(LD_InputType),           intent(inout) :: Inputs(:)     !< Inputs at InputTimes (output from this routine only
   real(DbKi),                   intent(in   ) :: InputTimes(:) !< Times in seconds associated with Inputs
   type(LD_ParameterType),       intent(in   ) :: p             !< Parameters
   type(LD_ContinuousStateType), intent(inout) :: x             !< Input: Continuous states at t;  Output: at t+dt
   type(LD_DiscreteStateType),   intent(inout) :: xd            !< Input: Discrete states at t;    Output: at t+dt
   type(LD_ConstraintStateType), intent(inout) :: z             !< Input: Constraint states at t;  Output: at t+dt
   type(LD_OtherStateType),      intent(inout) :: OtherState    !< Other states: Other states at t;Output: at t+dt
   type(LD_MiscVarType),         intent(inout) :: m             !< Misc variables for optimization (not copied in glue code)
   integer(IntKi),               intent(out)   :: errStat       !< Error status of the operation
   character(*),                 intent(out)   :: errMsg        !< Error message if errStat /= ErrID_None
   ! Initialize variables
   errStat   = ErrID_None           ! no error has occurred
   errMsg    = ""
   if (allocated(p%PrescribedValues)) then
      call interpTimeValue(p%PrescribedValues, t+p%dt, OtherState%iMotionInterpLast, m%qPrescribed(:))
      x%q(1:p%nq) = m%qPrescribed(1:p%nq)
   endif
   if ( p%nq == 0) return 
   if (p%IntMethod .eq. 1) then 
      call LD_RK4( t, n, Inputs, InputTimes, p, x, xd, z, OtherState, m, errStat, errMsg )
   elseif (p%IntMethod .eq. 2) then
      call LD_AB4( t, n, Inputs, InputTimes, p, x, xd, z, OtherState, m, errStat, errMsg )
   elseif (p%IntMethod .eq. 3) then
      call LD_ABM4( t, n, Inputs, InputTimes, p, x, xd, z, OtherState, m, errStat, errMsg )
   else  
      call SeterrStat(ErrID_Fatal,'Invalid time integration method:'//Num2LStr(p%IntMethod),errStat,errMsg,'LD_UpdateState') 
   end if
end subroutine LD_UpdateStates
!----------------------------------------------------------------------------------------------------------------------------------
!> This is a routine for computing outputs, used in both loose and tight coupling.
subroutine LD_CalcOutput( t, u, p, x, xd, z, OtherState, y, m, errStat, errMsg )
   real(DbKi),                   intent(in   ) :: t          !< Current simulation time in seconds
   type(LD_InputType),           intent(in   ) :: u          !< Inputs at t
   type(LD_ParameterType),       intent(in   ) :: p          !< Parameters
   type(LD_ContinuousStateType), intent(in   ) :: x          !< Continuous states at t
   type(LD_DiscreteStateType),   intent(in   ) :: xd         !< Discrete states at t
   type(LD_ConstraintStateType), intent(in   ) :: z          !< Constraint states at t
   type(LD_OtherStateType),      intent(in   ) :: OtherState !< Other states at t
   type(LD_MiscVarType),         intent(inout) :: m          !< Misc variables for optimization (not copied in glue code)
   type(LD_OutputType),          intent(inout) :: y          !< Outputs computed at t (Input only so that mesh con-
   integer(IntKi),               intent(out)   :: errStat    !< Error status of the operation
   character(*),                 intent(out)   :: errMsg     !< Error message if errStat /= ErrID_None
   type(LD_ContinuousStateType) :: dxdt                                                 !< 
   integer(IntKi)               :: errStat2    ! Status of error message
   character(1024)              :: errMsg2     ! Error message if ErrStat /= ErrID_None
   errStat   = ErrID_None           ! no error has occurred
   errMsg    = ""

   ! --- Compute accelerations
   if (allocated(p%PrescribedValues)) then
      y%xdd(1:p%nx) = m%qPrescribed(p%nq+1:p%nq+p%nx)
   else
      call LD_CalcContStateDeriv(t, u, p, x, xd, z, OtherState, m, dxdt, errStat2, errMsg2)
      y%xdd(1:p%nx) = dxdt%q(p%nx+1:p%nq)
   endif

   !--- Computing outputs:  y = Cx + Du  (optional)

   ! --- Write Outputs
   y%WriteOutput(1:2*p%nx) = x%q(1:p%nq)            ! Positions and velocities
   y%WriteOutput(2*p%nx+1:3*p%nx) = y%xdd(1:p%nx)   ! Accelerations
   y%WriteOutput(3*p%nx+1:4*p%nx) = u%Fext(1:p%nx)  ! Forces

contains
   logical function Failed()
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'LD_CalcOutput' )
      Failed =  errStat >= AbortErrLev
   end function Failed
end subroutine LD_CalcOutput
!----------------------------------------------------------------------------------------------------------------------------------
!> Tight coupling routine for computing derivatives of continuous states.
subroutine LD_CalcContStateDeriv( t, u, p, x, xd, z, OtherState, m, dxdt, errStat, errMsg )
   real(DbKi),                   intent(in   ) :: t          !< Current simulation time in seconds
   type(LD_InputType),           intent(in   ) :: u          !< Inputs at t
   type(LD_ParameterType),       intent(in   ) :: p          !< Parameters
   type(LD_ContinuousStateType), intent(in   ) :: x          !< Continuous states at t
   type(LD_DiscreteStateType),   intent(in   ) :: xd         !< Discrete states at t
   type(LD_ConstraintStateType), intent(in   ) :: z          !< Constraint states at t
   type(LD_OtherStateType),      intent(in   ) :: OtherState !< Other states at t
   type(LD_MiscVarType),         intent(inout) :: m          !< Misc variables for optimization (not copied in glue code)
   type(LD_ContinuousStateType), intent(out)   :: dxdt       !< Continuous state derivatives at t
   integer(IntKi),               intent(out)   :: errStat    !< Error status of the operation
   character(*),                 intent(out)   :: errMsg     !< Error message if errStat /= ErrID_None
   ! Local variables
   integer(IntKi)  :: iDOF
   integer(IntKi)  :: errStat2    ! Status of error message
   character(1024) :: errMsg2     ! Error message if ErrStat /= ErrID_None
   ! Initialize variables
   errStat   = ErrID_None           ! no error has occurred
   errMsg    = ""
   ! Allocation of output dxdt (since intent(out))
   call AllocAry(dxdt%q, p%nq, 'dxdt%q', errStat2, errMsg2); if(Failed()) return
   if ( p%nq == 0 ) return

   ! --- Computation of dq
   ! >>> MATMUL IMPLEMENTATION 
   dxdt%q = matmul(p%AA,x%q) + matmul(p%BB,u%Fext)
   ! >>> BLAS IMPLEMENTATION 
   !           COPY( N   , X                    , inCX, Y      , inCY)
   !call LAPACK_COPY(p%nCB, x%qmdot              , 1  , dxdt%qm    , 1  ) ! qmdot=qmdot
   !!           GEMV(TRS, M    , N     , alpha    , A  , LDA , X  ,inCX, Beta   ,  Y         , IncY)
   !call LAPACK_GEMV('n', p%nq, p%nq  ,  1.0_ReKi, p%AA, p%nq, x%q   , 1  , 1.0_ReKi, dxdt%qmdot, 1   ) !        - K22 x2
   !call LAPACK_GEMV('n', p%nq, p%nx  ,  1.0_ReKi, p%BB, p%nq, u%Fext, 1  , 1.0_ReKi, dxdt%qmdot, 1   ) !        - M21 \ddot{x1}
   ! --- Desactivating Constant DOFs
   do iDOF = 1,p%nx
      if (.not. p%activeDOFs(iDOF)) then
         dxdt%q(iDOF     ) = 0.0_ReKi
         dxdt%q(iDOF+p%nx) = 0.0_ReKi
      endif
   enddo

contains
   logical function Failed()
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, 'LD_CalcContStateDeriv' )
      Failed =  errStat >= AbortErrLev
   end function Failed
end subroutine LD_CalcContStateDeriv
!----------------------------------------------------------------------------------------------------------------------------------
!> Setup outputs
subroutine Init_Outputs(p, m, y, InitInp, InitOut, errStat, errMsg)
   ! character(ChanLen),   intent(in)    :: OutList(:) !< list of user-requested outputs
   type(LD_ParameterType), intent(inout) :: p       !< module parameters
   type(LD_MiscVarType),   intent(inout) :: m       !< module misc
   type(LD_OutputType),    intent(inout) :: y       !< module outputs
   type(LD_InitInputType), intent(in   ) :: InitInp !< module init inputs
   type(LD_InitOutputType),intent(inout) :: InitOut !< module init outputs
   integer(intki),         intent(out)   :: errStat !< error status code
   character(*),           intent(out)   :: errMsg  !< error message, if an error occurred
   integer         :: errStat2    ! temporary (local) error status
   character(1024) :: errMsg2     ! Error message if ErrStat /= ErrID_None
   integer         :: i, iOut
   errStat = ErrID_None
   errMsg = ""

   ! --- Regular outputs
   call AllocAry(y%xdd, p%nx, 'qd', errStat2, errMsg2); if(Failed()) return
   y%xdd = 0.0_ReKi

   ! --- Write Outputs
   p%NumOuts = (p%nx) * (1 + 1 + 1 + 1) ! Pos, Vel, Acc, Force

   !call AllocAry(m%AllOuts, p%NumOuts, "LinDyn AllOut", errStat,errMsg ); if(Failed()) return; m%AllOuts(:) = 0.0_ReKi
   call AllocAry(y%WriteOutput,         p%NumOuts,'WriteOutput',   errStat,errMsg); if(Failed()) return
   call AllocAry(InitOut%WriteOutputHdr,p%NumOuts,'WriteOutputHdr',errStat,errMsg); if(Failed()) return
   call AllocAry(InitOut%WriteOutputUnt,p%NumOuts,'WriteOutputUnt',errStat,errMsg); if(Failed()) return
   y%WriteOutput(1:p%NumOuts) = 0.0

   ! Sanity checks
   if (.not. allocated(InitInp%DOFsNames)) then
      errStat2 = errID_Fatal; errMsg2='DOFs Names not allocated'; if(Failed()) return
   else
      if(size(InitInp%DOFsNames)/=p%nx) then
         errStat2 = errID_Fatal; errMsg2='Shape of DOFs Names incorrect'; if(Failed()) return
      endif
      if (.not.allocated(InitInp%DOFsUnits)) then
         errStat2 = errID_Fatal; errMsg2='DOFs Units should be allocated if Names are provided'; if(Failed()) return
      endif
      if(size(InitInp%DOFsUnits)/=p%nx) then
         errStat2 = errID_Fatal; errMsg2='Shape of DOFs Units incorrect'; if(Failed()) return
      endif
   endif

   iOut = 0 ! Cumulative counter
   call SetWriteOutputsForDOFs(''  ) ! Positions
   call SetWriteOutputsForDOFs('d' ) ! Velocities
   call SetWriteOutputsForDOFs('dd') ! Accelerations
   call SetWriteOutputsForDOFs('f' ) ! Forces

   ! If using OutParam instead
   !InitOut%WriteOutputHdr(1:p%NumOuts) = p%OutParam(1:p%NumOuts)%Name
   !InitOut%WriteOutputUnt(1:p%NumOuts) = p%OutParam(1:p%NumOuts)%Units     
   ! Debug output to screen
   !do i = 1,p%NumOuts
   !   print*,i, InitOut%WriteOutputHdr(i), InitOut%WriteOutputUnt(i)
   !enddo
   
contains
   subroutine SetWriteOutputsForDOFs(sPrefix)
      character(len=*) :: sPrefix
      do i = 1, p%nx
         iOut = iOut+1
         InitOut%WriteOutputHdr(iOut) = trim(InitInp%prefix)//trim(sPrefix)//trim(InitInp%DOFsNames(i))
         ! Units 
         if (sPrefix == '')   InitOut%WriteOutputUnt(iOut) ='('//trim(InitInp%DOFsUnits(i))//')'
         if (sPrefix == 'd')  InitOut%WriteOutputUnt(iOut) ='('//trim(InitInp%DOFsUnits(i))//'/s)'
         if (sPrefix == 'dd') InitOut%WriteOutputUnt(iOut) ='('//trim(InitInp%DOFsUnits(i))//'/s^2)'
         if (sPrefix == 'f') then 
            if     (InitInp%DOFsUnits(i)=='m')   then; InitOut%WriteOutputUnt(iOut) ='(N)' ;
            elseif (InitInp%DOFsUnits(i)=='rad') then; InitOut%WriteOutputUnt(iOut) ='(Nm)' ;
            else;    InitOut%WriteOutputUnt(iOut) ='(-)'
            endif
         endif
      enddo
   endsubroutine

   logical function Failed()
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, 'Init_Outputs' )
      Failed =  errStat >= AbortErrLev
   end function Failed
end subroutine Init_Outputs
!----------------------------------------------------------------------------------------------------------------------------------
!> Setup Linearization data
subroutine Init_Lin(p, InitOut, errStat, errMsg)
   type(LD_ParameterType), intent(in   ) :: p       !< module parameters
   type(LD_InitOutputType),intent(inout) :: InitOut !< module init outputs
   integer(intki),         intent(out)   :: errStat !< error status code
   character(*),           intent(out)   :: errMsg  !< error message, if an error occurred
   integer                 :: errStat2                                        ! temporary (local) error status
   character(1024)         :: errMsg2     ! Error message if ErrStat /= ErrID_None
   integer                 :: i, nu
   errStat = ErrID_None
   errMsg = ""
   nu = p%nx

!   LinNames_y                      {:}           -         -        "Names of the outputs used in linearization" -
!   LinNames_x                      {:}           -         -        "Names of the continuous states used in linearization" -
!   LinNames_u                      {:}           -         -        "Names of the inputs used in linearization" -
!   RotFrame_y                      {:}           -         -        "Flag that tells FAST/MBC3 if the outputs used in linearization are in the rotating frame"	-
!   RotFrame_x                      {:}           -         -        "Flag that tells FAST/MBC3 if the continuous states used in linearization are in the rotating frame"	-
!   RotFrame_u                      {:}           -         -        "Flag that tells FAST/MBC3 if the inputs used in linearization are in the rotating frame"	-
!   IsLoad_u                        {:}           -         -        "Flag that tells FAST if the inputs used in linearization are loads (for preconditioning matrix)" -
!   DerivOrder_x                    {:}           -         -        "Integer that tells FAST/MBC3 the maximum derivative order of continuous states used in linearization" -
   !Appropriate Jacobian row/column names and rotating-frame flags here:   
   call AllocAry(InitOut%LinNames_y  , p%NumOuts , 'LinNames_y', errStat, errMsg); if(Failed()) return
   call AllocAry(InitOut%RotFrame_y  , p%NumOuts , 'RotFrame_y', errStat, errMsg); if(Failed()) return
   call AllocAry(InitOut%LinNames_x  , p%nq      , 'LinNames_x', errStat, errMsg); if(Failed()) return
   call AllocAry(InitOut%RotFrame_x  , p%nq      , 'RotFrame_x', errStat, errMsg); if(Failed()) return
   call AllocAry(InitOut%DerivOrder_x, p%nq      , 'DerivOrd_x', errStat, errMsg); if(Failed()) return
   call AllocAry(InitOut%LinNames_u  , nu        , 'LinNames_u', errStat, errMsg); if(Failed()) return
   call AllocAry(InitOut%RotFrame_u  , nu        , 'RotFrame_u', errStat, errMsg); if(Failed()) return
   call AllocAry(InitOut%IsLoad_u    , nu        , 'IsLoad_u'  , errStat, errMsg); if(Failed()) return
   InitOut%DerivOrder_x(:)=2
   ! LinNames_y
   do i=1, p%NumOuts
      InitOut%LinNames_y(i) = trim(InitOut%WriteOutputHdr(i))//', '//trim(InitOut%WriteOutputUnt(i))
      print*,'y',i, trim(InitOut%LinNames_y(i))
   enddo
   ! LinNames_u
   do i=1, p%nx
      InitOut%LinNames_u(i) = trim(InitOut%WriteOutputHdr(3*p%nx+ i))//', '//trim(InitOut%WriteOutputUnt(3*p%nx+i))
      print*,'u',i, trim(InitOut%LinNames_u(i))
   enddo
   ! LinNames_x
   do I=1,p%nq; 
       InitOut%LinNames_x(I) = trim(InitOut%WriteOutputHdr(i))//', '//trim(InitOut%WriteOutputUnt(i))
      print*,'x',i, trim(InitOut%LinNames_x(i))
   enddo
   InitOut%RotFrame_x = .false.
   InitOut%RotFrame_y = .false.
   InitOut%RotFrame_u = .false.
   InitOut%IsLoad_u   = .true. 
   !   
contains
   logical function Failed()
      if (errStat >= AbortErrLev) errMsg = 'LD_JacobianLin:'//trim(errMsg)
      Failed =  errStat >= AbortErrLev
   end function Failed
end subroutine Init_Lin

!----------------------------------------------------------------------------------------------------------------------------------
!> Jacobians with respect to inputs (u)
subroutine LD_JacobianPInput( t, u, p, x, xd, z, OtherState, y, m, errStat, errMsg, dYdu, dXdu, dXddu, dZdu)
   real(DbKi),                        intent(in   ) :: t          !< Time in seconds at operating point
   type(LD_InputType),                intent(in   ) :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
   type(LD_ParameterType),            intent(in   ) :: p          !< Parameters
   type(LD_ContinuousStateType),      intent(in   ) :: x          !< Continuous states at operating point
   type(LD_DiscreteStateType),        intent(in   ) :: xd         !< Discrete states at operating point
   type(LD_ConstraintStateType),      intent(in   ) :: z          !< Constraint states at operating point
   type(LD_OtherStateType),           intent(in   ) :: OtherState !< Other states at operating point
   type(LD_OutputType),               intent(in   ) :: y          !< Output (change to inout if a mesh copy is required); 
   type(LD_MiscVarType),              intent(inout) :: m          !< Misc/optimization variables
   integer(IntKi),                    intent(  out) :: errStat    !< Error status of the operation
   character(*),                      intent(  out) :: errMsg     !< Error message if errStat /= ErrID_None
   real(R8Ki), allocatable, optional, intent(inout) :: dYdu(:,:)  !< Jacobians of output functions (Y) with respect to (u)
   real(R8Ki), allocatable, optional, intent(inout) :: dXdu(:,:)  !< Jacobians of continuous state functions (X) with respect to (u)
   real(R8Ki), allocatable, optional, intent(inout) :: dXddu(:,:) !< Jacobians of discrete state functions (Xd) with respect to (u)
   real(R8Ki), allocatable, optional, intent(inout) :: dZdu(:,:)  !< Jacobians of constraint state functions (Z) with respect to (u)
   integer(IntKi) :: i, nu  ! Loop index
   ! Initialize errStat
   errStat = ErrID_None
   errMsg  = ''
   nu = p%nx
   if (present(dYdu)) then
      if (.not. allocated(dYdu)) then
          call AllocAry(dYdu, p%NumOuts, nu, 'dYdu', errStat, errMsg); if(Failed()) return
          dYdu(:,:) = 0.0_ReKi
      end if
      !dYdu(1        :  p%nx, :)  = 0.0_ReKi        ! Positions 
      dYdu(  p%nx+1 : 3*p%nx, :)  = p%BB            ! Velocities and accelerations
      do i=1, p%nx ; dYdu(3*p%nx+i, i)  = 1.0_ReKi;  enddo ! Forces (which are inputs)
   end if
   if (present(dXdu)) then
      if (.not. allocated(dXdu)) then
          call AllocAry(dXdu, p%nq, nu, 'dXdu', errStat, errMsg); if(Failed()) return
          dXdu(:,:) = 0.0_ReKi
      end if
      dXdu = p%BB
   end if
   if (present(dXddu)) then
   end if
   if (present(dZdu)) then
   end if
contains
    logical function Failed()
        if (errStat >= AbortErrLev) errMsg = 'LD_JacobianPInput:'//trim(errMsg)
        Failed =  errStat >= AbortErrLev
    end function Failed
end subroutine LD_JacobianPInput
!----------------------------------------------------------------------------------------------------------------------------------
!> Jacobians with respect to continuous states (x)
subroutine LD_JacobianPContState( t, u, p, x, xd, z, OtherState, y, m, errStat, errMsg, dYdx, dXdx, dXddx, dZdx )
   real(DbKi),                        intent(in   ) :: t          !< Time in seconds at operating point
   type(LD_InputType),                intent(in   ) :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
   type(LD_ParameterType),            intent(in   ) :: p          !< Parameters
   type(LD_ContinuousStateType),      intent(in   ) :: x          !< Continuous states at operating point
   type(LD_DiscreteStateType),        intent(in   ) :: xd         !< Discrete states at operating point
   type(LD_ConstraintStateType),      intent(in   ) :: z          !< Constraint states at operating point
   type(LD_OtherStateType),           intent(in   ) :: OtherState !< Other states at operating point
   type(LD_OutputType),               intent(in   ) :: y          !< Output (change to inout if a mesh copy is required);
   type(LD_MiscVarType),              intent(inout) :: m          !< Misc/optimization variables
   integer(IntKi),                    intent(out)   :: errStat    !< Error status of the operation
   character(*),                      intent(out)   :: errMsg     !< Error message if errStat /= ErrID_None
   real(R8Ki), allocatable, optional, intent(inout) :: dYdx(:,:)  !< Jacobians of output functions (Y) with respect to (x)
   real(R8Ki), allocatable, optional, intent(inout) :: dXdx(:,:)  !< Jacobians of continuous state functions (X) with respect to (x)
   real(R8Ki), allocatable, optional, intent(inout) :: dXddx(:,:) !< Jacobians of discrete state functions (Xd) with respect to (x)
   real(R8Ki), allocatable, optional, intent(inout) :: dZdx(:,:)  !< Jacobians of constraint state functions (Z) with respect to (x)
   integer(IntKi) :: i    ! Loop index
   ! Initialize errStat
   errStat = ErrID_None
   errMsg  = ''
   if (present(dYdx)) then
      ! allocate and set dYdx
      if (.not. allocated(dYdx)) then
          call AllocAry(dYdx, p%NumOuts, p%nq, 'dYdx', errStat, errMsg); if(Failed()) return
          dYdx(:,:) = 0.0_ReKi
      end if
      do i=1,p%nx;  dYdx(i,i) = 1.0_ReKi; enddo ! Position
      dYdx(p%nx+1:3*p%nx,:  ) = p%AA            ! Velocity and acceleration
      !dYdx(3*p%nx+1:,:) = 0                    ! Forces
   end if
   if (present(dXdx)) then
      ! allocate and set dXdx
      if (.not. allocated(dXdx)) then
          call AllocAry(dXdx, p%nq, p%nq, 'dXdx', errStat, errMsg); if(Failed()) return
          dXdx(:,:) = 0.0_ReKi
      end if
      dXdx = p%AA
   end if
   if (present(dXddx)) then
   end if
   if (present(dZdx)) then
   end if
contains
   logical function Failed()
      if (errStat >= AbortErrLev) errMsg = 'LD_JacobianPContState:'//trim(errMsg)
      Failed =  errStat >= AbortErrLev
   end function Failed
end subroutine LD_JacobianPContState
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine to pack the data structures representing the operating points into arrays for linearization.
subroutine LD_GetOP( t, u, p, x, xd, z, OtherState, y, m, errStat, errMsg, u_op, y_op, x_op, dx_op, xd_op, z_op )
   real(DbKi),                        intent(in   ) :: t          !< Time in seconds at operating point
   type(LD_InputType),                intent(in   ) :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
   type(LD_ParameterType),            intent(in   ) :: p          !< Parameters
   type(LD_ContinuousStateType),      intent(in   ) :: x          !< Continuous states at operating point
   type(LD_DiscreteStateType),        intent(in   ) :: xd         !< Discrete states at operating point
   type(LD_ConstraintStateType),      intent(in   ) :: z          !< Constraint states at operating point
   type(LD_OtherStateType),           intent(in   ) :: OtherState !< Other states at operating point
   type(LD_OutputType),               intent(in   ) :: y          !< Output at operating point
   type(LD_MiscVarType),              intent(inout) :: m          !< Misc/optimization variables
   integer(IntKi),                    intent(out)   :: errStat    !< Error status of the operation
   character(*),                      intent(out)   :: errMsg     !< Error message if errStat /= ErrID_None
   real(ReKi), allocatable, optional, intent(inout) :: u_op(:)    !< values of linearized inputs
   real(ReKi), allocatable, optional, intent(inout) :: y_op(:)    !< values of linearized outputs
   real(ReKi), allocatable, optional, intent(inout) :: x_op(:)    !< values of linearized continuous states
   real(ReKi), allocatable, optional, intent(inout) :: dx_op(:)   !< values of first time derivatives of linearized continuous states
   real(ReKi), allocatable, optional, intent(inout) :: xd_op(:)   !< values of linearized discrete states
   real(ReKi), allocatable, optional, intent(inout) :: z_op(:)    !< values of linearized constraint states
   integer(IntKi) :: i, nu
   type(LD_ContinuousStateType) :: dx          !< derivative of continuous states at operating point
   ! Initialize errStat
   errStat = ErrID_None
   errMsg  = ''
   nu = p%nx

   if ( present( u_op ) ) then
       if (.not. allocated(u_op)) then
           call AllocAry(u_op, nu, 'u_op', errStat, errMsg); if(Failed())return
       endif
       u_op(:) = u%Fext
   end if

   if ( present( y_op ) ) then
       if (.not. allocated(y_op)) then
           call AllocAry(y_op, p%NumOuts, 'y_op', errStat, errMsg); if(Failed())return
       endif
       ! Update the output mesh
       do i=1,p%NumOuts         
           y_op(i) = y%WriteOutput(i)
       end do      
   end if

   if ( present( x_op ) ) then
       if (.not. allocated(x_op)) then
           call AllocAry(x_op, p%nq, 'x_op', errStat, errMsg); if (Failed())return
       endif
       x_op = x%q
   end if

   if ( present( dx_op ) ) then
       if (.not. allocated(dx_op)) then
           call AllocAry(dx_op, p%nq, 'dx_op', errStat, errMsg); if (Failed())return
       endif
       call LD_CalcContStateDeriv(t, u, p, x, xd, z, OtherState, m, dx, errStat, errMsg); if(Failed()) return
       dx_op = dx%q
   end if

   if ( present( xd_op ) ) then
   end if
   
   if ( present( z_op ) ) then
   end if

contains
    logical function Failed()
      if (errStat >= AbortErrLev) errMsg = 'LD_GetOP:'//trim(errMsg)
      Failed =  errStat >= AbortErrLev
    end function Failed
end subroutine LD_GetOP

end module LinDyn
!**********************************************************************************************************************************
