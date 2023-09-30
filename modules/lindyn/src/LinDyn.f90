!**********************************************************************************************************************************
!> LinDyn, module for linear dynamical system with mass, stiffness and damping matrix
! ..................................................................................................................................
!! ## LICENSinG
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
! 
!    private

!    public :: LD_Init                           !  Initialization routine
!    public :: LD_End                            !  Ending routine (includes clean up)
!    public :: LD_UpdateStates                   !  Loose coupling routine for solving for constraint states, integrating continuous states, and updating discrete states
!    public :: LD_CalcOutput                     !  Routine for computing outputs
!    public :: LD_CalcConstrStateResidual        !  Tight coupling routine for returning the constraint state residual
!    public :: LD_CalcContStateDeriv             !  Tight coupling routine for computing derivatives of continuous states
!    public :: LD_UpdateDiscState                !  Tight coupling routine for updating discrete states
!    public :: LD_JacobianPInput                 !  Jacobians of (y, x, xd, z) with respect to the inputs (u)
!    public :: LD_JacobianPContState             !  Jacobians of (y, x, xd, z) with respect to the continuous (x)
!    public :: LD_JacobianPDiscState             !  Jacobians of (y, x, xd, z) with respect to the discrete states (xd)
!    public :: LD_JacobianPConstrState           !  Jacobians of (y, x, xd, z) with respect to the constraint states (z)
!    public :: LD_GetOP                          !  Routine to get the operating-point values for linearization (from data structures to arrays)
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

   print*,''
   print*,'M',p%MM(1,:)
   print*,'M',p%MM(2,:)
   print*,'M',p%MM(3,:)
   print*,''
   print*,'C',p%CC(1,:)
   print*,'C',p%CC(2,:)
   print*,'C',p%CC(3,:)
   print*,''
   print*,'K',p%KK(1,:)
   print*,'K',p%KK(2,:)
   print*,'K',p%KK(3,:)
   print*,''

   call StateMatrices(p%MM, p%CC, p%KK, p%AA, p%BB, errStat2, errMsg2); if(Failed()) return

   print*,''
   print*,'A',p%AA(1,:)
   print*,'A',p%AA(2,:)
   print*,'A',p%AA(3,:)
   print*,'A',p%AA(4,:)
   print*,'A',p%AA(5,:)
   print*,'A',p%AA(6,:)
   print*,''
   print*,'B',p%BB(1,:)
   print*,'B',p%BB(2,:)
   print*,'B',p%BB(3,:)
   print*,'B',p%BB(4,:)
   print*,'B',p%BB(5,:)
   print*,'B',p%BB(6,:)

   ! --- Allocate STates
   call AllocAry( x%q    , p%nq,'DOFs' , errStat,errMsg); if(Failed()) return
   x%q(     1:p%nx) = InitInp%x0
   x%q(p%nx+1:p%nq) = InitInp%xd0

   ! allocate OtherState%xdot if using multi-step method; initialize n
   if ( ( p%IntMethod .eq. 2) .OR. ( p%IntMethod .eq. 3)) THEN
       allocate( OtherState%xdot(4), STAT=errStat2); errMsg2='Error allocating OtherState%xdot'
       if(Failed()) return
   endif

   ! --- Initialize Misc Variables:

!    ! Define initial guess (set up mesh first) for the system inputs here:
!    call Init_meshes(u, y, InitInp, errStat, errMsg); if(Failed()) return
   ! --- Guess inputs
   call AllocAry(u%Fext, p%nx, 'Fext', errStat2, errMsg2); if(Failed()) return
   u%Fext=0.0_ReKi

   ! --- Outputs
   call AllocAry(y%qd, p%nx, 'qd', errStat2, errMsg2); if(Failed()) return
   y%qd = 0.0_ReKi
   y%qd(1:p%nx) = InitInp%xd0

   ! --- Write Outputs
   p%NumOuts = 0
!    ! Setting p%OutParam from OutList
!    call SetOutParam(InputFileData%OutList, InputFileData%NumOuts, p, errStat, errMsg); if(Failed()) return
!   
   call AllocAry( m%AllOuts, p%NumOuts, "LinDyn AllOut", errStat,errMsg ); if(Failed()) return
   m%AllOuts(:) = 0.0_ReKi
   call AllocAry( y%WriteOutput,        p%NumOuts,'WriteOutput',   errStat,errMsg); if(Failed()) return
   call AllocAry(InitOut%WriteOutputHdr,p%NumOuts,'WriteOutputHdr',errStat,errMsg); if(Failed()) return
   call AllocAry(InitOut%WriteOutputUnt,p%NumOuts,'WriteOutputUnt',errStat,errMsg); if(Failed()) return
   y%WriteOutput(1:p%NumOuts) = 0.0
   !InitOut%WriteOutputHdr(1:p%NumOuts) = p%OutParam(1:p%NumOuts)%Name
   !InitOut%WriteOutputUnt(1:p%NumOuts) = p%OutParam(1:p%NumOuts)%Units     
   InitOut%Ver = LD_Ver
!       
!    if (InitInp%Linearize) then
!       ! TODO The linearization features are in place but waiting for glue-code changes, and testing.
!       call SeterrStat( ErrID_Fatal, 'LinDyn linearization analysis is currently not supported by the glue code.', errStat, errMsg, 'LD_Init');
!       if(Failed())return
!       !Appropriate Jacobian row/column names and rotating-frame flags here:   
!       call AllocAry(InitOut%LinNames_y, 6+p%NumOuts , 'LinNames_y', errStat, errMsg); if(Failed()) return
!       call AllocAry(InitOut%RotFrame_y, 6+p%NumOuts , 'RotFrame_y', errStat, errMsg); if(Failed()) return
!       call AllocAry(InitOut%LinNames_x, 2*p%nCB     , 'LinNames_x', errStat, errMsg); if(Failed()) return
!       call AllocAry(InitOut%RotFrame_x, 2*p%nCB     , 'RotFrame_x', errStat, errMsg); if(Failed()) return
!       call AllocAry(InitOut%DerivOrder_x, 2*p%nCB   , 'DerivOrd_x', errStat, errMsg); if(Failed()) return
!       call AllocAry(InitOut%LinNames_u, N_inPUTS    , 'LinNames_u', errStat, errMsg); if(Failed()) return
!       call AllocAry(InitOut%RotFrame_u, N_inPUTS    , 'RotFrame_u', errStat, errMsg); if(Failed()) return
!       call AllocAry(InitOut%IsLoad_u  , N_inPUTS    , 'IsLoad_u'  , errStat, errMsg); if(Failed()) return
!       InitOut%DerivOrder_x(:)=2
!       ! LinNames_y
!       do I=1,3; 
!           InitOut%LinNames_y(I)   = 'Interface node '//XYZ(I)//' force, N' 
!           InitOut%LinNames_y(I+3) = 'Interface node '//XYZ(I)//' moment, Nm' 
!       enddo
!       do i=1,p%NumOuts
!           InitOut%LinNames_y(N_outPUTS+i) = trim(p%OutParam(i)%Name)//', '//p%OutParam(i)%Units
!       end do
!       ! LinNames_u
!       do I=1,3;
!           InitOut%LinNames_u(I+ 0) = 'Interface node '//XYZ(I)//' translation displacement, m'
!           InitOut%LinNames_u(I+ 3) = 'Interface node '//XYZ(I)//' rotation, rad'
!           InitOut%LinNames_u(I+ 6) = 'Interface node '//XYZ(I)//' translation velocity, m/s'
!           InitOut%LinNames_u(I+ 9) = 'Interface node '//XYZ(I)//' rotation velocity, rad/s'
!           InitOut%LinNames_u(I+12) = 'Interface node '//XYZ(I)//' translation acceleration, m/s^2'
!           InitOut%LinNames_u(I+15) = 'Interface node '//XYZ(I)//' rotation acceleration, rad/s^2'
!       enddo
!       ! LinNames_x
!       do I=1,p%nCB; 
!           InitOut%LinNames_x(I)       = 'Mode '//trim(Num2LStr(p%ActiveCBDOF(I)))//' displacement, -';
!           InitOut%LinNames_x(I+p%nCB) = 'Mode '//trim(Num2LStr(p%ActiveCBDOF(I)))//' velocity, -';
!       enddo
!       InitOut%RotFrame_x = .false. ! note that meshes are in the global, not rotating frame
!       InitOut%RotFrame_y = .false. ! note that meshes are in the global, not rotating frame
!       InitOut%RotFrame_u = .false. ! note that meshes are in the global, not rotating frame
!       InitOut%IsLoad_u   = .false. ! the inputs are not loads but kinematics
!    end if
! 
!    ! --- Summary file 
!    if (InputFileData%SumPrint) then
!        call LD_PrintSum(x, p, m, InitInp%RootName, errStat, errMsg); if(Failed()) return
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
!> Allocate init input data for module based on number of degrees of freedom
subroutine LD_InitInputData(nx, InitInp, errStat, errMsg)
   integer(IntKi),               intent(in ) :: nx        !< Number of degrees of freedom
   type(LD_InitInputType),       intent(out) :: InitInp     !< Input data for initialization routine
   integer(IntKi),               intent(out) :: errStat   !< Error status of the operation
   character(*),                 intent(out) :: errMsg    !< Error message if errStat /= ErrID_None
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
   InitInp%MM         = 0.0_ReKi
   InitInp%CC         = 0.0_ReKi
   InitInp%KK         = 0.0_ReKi
   InitInp%x0         = 0.0_ReKi
   InitInp%xd0        = 0.0_ReKi
   InitInp%activeDOFs = .True.
contains
   logical function Failed()
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'LD_Init' )
      Failed = ErrStat >= AbortErrLev
   end function Failed
end subroutine LD_InitInputData
!----------------------------------------------------------------------------------------------------------------------------------
!> Compute A and B state matrices for a linear mechanical system
!! NOTE: Generic function (no derived types), keep it that way
!! A = [   0        I       ]     B = [0      ]
!!     [-M^{-1}K   -M^{-1}C ]       = [-M^{-1}]
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
   allocate(WORK(LWORk))
   call LAPACK_getri(nx, MinvX, IPIV, WORK, LWORK, errStat2, errMsg2); if(Failed()) return
   BB(nx+1:nq, : ) = -MinvX

contains
   logical function Failed()
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'LD_Init' )
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
! subroutine Init_meshes(u, y, InitInp, errStat, errMsg)
!    type(LD_InputType),           intent(inout)  :: u           !< System inputs
!    type(LD_OutputType),          intent(inout)  :: y           !< System outputs
!    type(LD_InitInputType),       intent(in   )  :: InitInp     !< Input data for initialization routine
!    integer(IntKi),                    intent(  out)  :: errStat     !< Error status of the operation
!    character(*),                      intent(  out)  :: errMsg      !< Error message if errStat /= ErrID_None
!    ! Create the input and output meshes associated with platform loads
!    call MeshCreate(  BlankMesh         = u%PtfmMesh       , &
!                      IOS               = COMPONENT_inPUT  , &
!                      Nnodes            = 1                , &
!                      errStat           = errStat          , &
!                      ErrMess           = errMsg           , &
!                      TranslationDisp   = .TRUE.           , &
!                      Orientation       = .TRUE.           , &
!                      TranslationVel    = .TRUE.           , &
!                      RotationVel       = .TRUE.           , &
!                      TranslationAcc    = .TRUE.           , &
!                      RotationAcc       = .TRUE.)
!    if(Failed()) return
!       
!    ! Create the node on the mesh, the node is located at the PlatformRefzt, to match ElastoDyn
!    call MeshPositionNode (u%PtfmMesh, 1, (/0.0_ReKi, 0.0_ReKi, InitInp%PtfmRefzt/), errStat, errMsg ); if(Failed()) return
!    ! Create the mesh element
!    call MeshConstructElement (  u%PtfmMesh, ELEMENT_POinT, errStat, errMsg, 1 ); if(Failed()) return
!    call MeshCommit ( u%PtfmMesh, errStat, errMsg ); if(Failed()) return
!    ! the output mesh is a sibling of the input:
!    call MeshCopy( SrcMesh=u%PtfmMesh, DestMesh=y%PtfmMesh, CtrlCode=MESH_SIBLinG, IOS=COMPONENT_outPUT, &
!                   errStat=errStat, ErrMess=errMsg, Force=.TRUE., Moment=.TRUE. )
!    if(Failed()) return
! CONTAinS
!     logical function Failed()
!         call SeterrStatSimple(errStat, errMsg, 'Init_meshes')
!         Failed =  errStat >= AbortErrLev
!     end function Failed
! end subroutine Init_meshes
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
   ! Local variables
   type(LD_ContinuousStateType) :: dxdt !< 
   integer(IntKi)  :: errStat2    ! Status of error message
   character(1024) :: errMsg2     ! Error message if ErrStat /= ErrID_None
!    integer(IntKi)                                  :: I                 !< Generic counters
!    real(ReKi), dimension(6)                        :: Fc                !< Output coupling force
   ! Initialize variables
   errStat   = ErrID_None           ! no error has occurred
   errMsg    = ""
! 
   ! --- Compute accelerations
   call LD_CalcContStateDeriv(t, u, p, x, xd, z, OtherState, m, dxdt, errStat2, errMsg2)
   y%qd = dxdt%q

   !--- Computing output:  y = Cx + Du + Fy  
! 
!    ! Update the output mesh
!    do i=1,3
!       y%PtfmMesh%Force(I,1)  = Fc(I)
!       y%PtfmMesh%Moment(I,1) = Fc(I+3)
!    enddo
! 
!    ! --- All Outputs
!    m%AllOuts(ID_PtfFx) = y%PtfmMesh%Force (1,1)
!    m%AllOuts(ID_PtfFy) = y%PtfmMesh%Force (2,1)
!    m%AllOuts(ID_PtfFz) = y%PtfmMesh%Force (3,1)
!    m%AllOuts(ID_PtfMx) = y%PtfmMesh%Moment(1,1)
!    m%AllOuts(ID_PtfMy) = y%PtfmMesh%Moment(2,1)
!    m%AllOuts(ID_PtfMz) = y%PtfmMesh%Moment(3,1)
!    m%AllOuts(ID_InpFx) = m%F_at_t(1)
!    m%AllOuts(ID_InpFy) = m%F_at_t(2)
!    m%AllOuts(ID_InpFz) = m%F_at_t(3)
!    m%AllOuts(ID_InpMx) = m%F_at_t(4)
!    m%AllOuts(ID_InpMy) = m%F_at_t(5)
!    m%AllOuts(ID_InpMz) = m%F_at_t(6)
!    !y%WriteOutput(ID_WaveElev) = .. ! TODO
!    do i=1,p%nCB
!       m%AllOuts(ID_QStart + 0*p%nCBFull -1 + p%ActiveCBDOF(I)) = x%qm   (I)    ! CBQ  - DOF Positions
!       m%AllOuts(ID_QStart + 1*p%nCBFull -1 + p%ActiveCBDOF(I)) = x%qmdot(I)    ! CBQD - DOF Velocities
!       m%AllOuts(ID_QStart + 2*p%nCBFull -1 + p%ActiveCBDOF(I)) = m%F_at_t(6+I) ! CBF  - DOF Forces
!    enddo
!    ! --- Selected output channels only
!    do I = 1,p%NumOuts
!       if (p%OutParam(I)%Indx>0) then
!           y%WriteOutput(I) = p%OutParam(I)%SignM * m%AllOuts( p%OutParam(I)%Indx )
!       else
!           y%WriteOutput(I) = -9.9999e20
!       endif
!    enddo    
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
   integer(IntKi)  :: errStat2    ! Status of error message
   character(1024) :: errMsg2     ! Error message if ErrStat /= ErrID_None
!    integer(IntKi)                                    :: I
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
!
contains
   logical function Failed()
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'LD_CalcContStateDeriv' )
      Failed =  errStat >= AbortErrLev
   end function Failed
end subroutine LD_CalcContStateDeriv
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!>
! subroutine LD_JacobianPInput( t, u, p, x, xd, z, OtherState, y, m, errStat, errMsg, dYdu, dXdu, dXddu, dZdu)
!    real(DbKi),                         intent(in   ) :: t          !< Time in seconds at operating point
!    type(LD_InputType),            intent(in   ) :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
!    type(LD_ParameterType),        intent(in   ) :: p          !< Parameters
!    type(LD_ContinuousStateType),  intent(in   ) :: x          !< Continuous states at operating point
!    type(LD_DiscreteStateType),    intent(in   ) :: xd         !< Discrete states at operating point
!    type(LD_ConstraintStateType),  intent(in   ) :: z          !< Constraint states at operating point
!    type(LD_OtherStateType),       intent(in   ) :: OtherState !< Other states at operating point
!    type(LD_OutputType),           intent(in   ) :: y          !< Output (change to inout if a mesh copy is required); 
!                                                                    !!   Output fields are not used by this routine, but type is   
!                                                                    !!   available here so that mesh parameter information (i.e.,  
!                                                                    !!   connectivity) does not have to be recalculated for dYdu.
!    type(LD_MiscVarType),          intent(inout) :: m          !< Misc/optimization variables
!    integer(IntKi),                     intent(  out) :: errStat    !< Error status of the operation
!    character(*),                       intent(  out) :: errMsg     !< Error message if errStat /= ErrID_None
!    real(R8Ki), allocatable, OPTIONAL,  intent(inout) :: dYdu(:,:)  !< Partial derivatives of output functions (Y) with respect 
!                                                                    !!   to the inputs (u) [intent in to avoid deallocation]
!    real(R8Ki), allocatable, OPTIONAL,  intent(inout) :: dXdu(:,:)  !< Partial derivatives of continuous state functions (X) with 
!                                                                    !!   respect to the inputs (u) [intent in to avoid deallocation]
!    real(R8Ki), allocatable, OPTIONAL,  intent(inout) :: dXddu(:,:) !< Partial derivatives of discrete state functions (Xd) with 
!                                                                    !!   respect to the inputs (u) [intent in to avoid deallocation]
!    real(R8Ki), allocatable, OPTIONAL,  intent(inout) :: dZdu(:,:)  !< Partial derivatives of constraint state functions (Z) with 
!                                                                    !!   respect to the inputs (u) [intent in to avoid deallocation]
!    integer(IntKi) :: i,j  ! Loop index
!    integer(IntKi) :: idx  ! Index of output channel in AllOuts
!    ! Initialize errStat
!    errStat = ErrID_None
!    errMsg  = ''
!    if (present(dYdu)) then
!       ! allocate and set dYdu
!       if (.not. allocated(dYdu)) then
!           call AllocAry(dYdu, N_outPUTS+p%NumOuts, N_inPUTS, 'dYdu', errStat, errMsg); if(Failed()) return
!           do i=1,size(dYdu,1); do j=1,size(dYdu,2); dYdu(i,j)=0.0_ReKi; enddo;enddo
!       end if
!       dYdu(1:6,1:N_inPUTS) = p%DMat(1:6,1:N_inPUTS)
!       !dYdu is zero except if WriteOutput is the interface loads 
!       do i = 1,p%NumOuts
!           idx  = p%OutParam(i)%Indx
!           if     (idx==ID_PtfFx) then; dYdu(6+i,1:N_inPUTS) = p%DMat(1,1:N_inPUTS)
!           elseif (idx==ID_PtfFy) then; dYdu(6+i,1:N_inPUTS) = p%DMat(2,1:N_inPUTS)
!           elseif (idx==ID_PtfFx) then; dYdu(6+i,1:N_inPUTS) = p%DMat(3,1:N_inPUTS)
!           elseif (idx==ID_PtfMz) then; dYdu(6+i,1:N_inPUTS) = p%DMat(4,1:N_inPUTS)
!           elseif (idx==ID_PtfMy) then; dYdu(6+i,1:N_inPUTS) = p%DMat(5,1:N_inPUTS)
!           elseif (idx==ID_PtfMz) then; dYdu(6+i,1:N_inPUTS) = p%DMat(6,1:N_inPUTS)
!           else                       ; dYdu(6+i,1:N_inPUTS) = 0.0_ReKi
!           endif 
!       end do
!   end if
!    if (present(dXdu)) then
!       ! allocate and set dXdu
!       if (.not. allocated(dXdu)) then
!           call AllocAry(dXdu, 2*p%nCB, N_inPUTS, 'dXdu', errStat, errMsg); if(Failed()) return
!           do i=1,size(dXdu,1); do j=1,size(dXdu,2); dXdu(i,j)=0.0_ReKi; enddo;enddo
!       end if
!       dXdu(1:2*p%nCB,1:N_inPUTS) = p%BMat(1:2*p%nCB,1:N_inPUTS)
!    end if
!    if (present(dXddu)) then
!    end if
!    if (present(dZdu)) then
!    end if
! CONTAinS
!     logical function Failed()
!         call SeterrStatSimple(errStat, errMsg, 'LD_JacobianPInput')
!         Failed =  errStat >= AbortErrLev
!     end function Failed
! end subroutine LD_JacobianPInput
! !----------------------------------------------------------------------------------------------------------------------------------
! !> Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) functions
! !! with respect to the continuous states (x). The partial derivatives dY/dx, dX/dx, dXd/dx, and DZ/dx are returned.
! subroutine LD_JacobianPContState( t, u, p, x, xd, z, OtherState, y, m, errStat, errMsg, dYdx, dXdx, dXddx, dZdx )
! !..................................................................................................................................
!    real(DbKi),                         intent(in   ) :: t          !< Time in seconds at operating point
!    type(LD_InputType),            intent(in   ) :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
!    type(LD_ParameterType),        intent(in   ) :: p          !< Parameters
!    type(LD_ContinuousStateType),  intent(in   ) :: x          !< Continuous states at operating point
!    type(LD_DiscreteStateType),    intent(in   ) :: xd         !< Discrete states at operating point
!    type(LD_ConstraintStateType),  intent(in   ) :: z          !< Constraint states at operating point
!    type(LD_OtherStateType),       intent(in   ) :: OtherState !< Other states at operating point
!    type(LD_OutputType),           intent(in   ) :: y          !< Output (change to inout if a mesh copy is required); 
!                                                                    !!   Output fields are not used by this routine, but type is   
!                                                                    !!   available here so that mesh parameter information (i.e.,  
!                                                                    !!   connectivity) does not have to be recalculated for dYdx.
!    type(LD_MiscVarType),          intent(inout) :: m          !< Misc/optimization variables
!    integer(IntKi),                     intent(  out) :: errStat    !< Error status of the operation
!    character(*),                       intent(  out) :: errMsg     !< Error message if errStat /= ErrID_None
!    real(R8Ki), allocatable, OPTIONAL,  intent(inout) :: dYdx(:,:)  !< Partial derivatives of output functions
!                                                                    !!   (Y) with respect to the continuous
!                                                                    !!   states (x) [intent in to avoid deallocation]
!    real(R8Ki), allocatable, OPTIONAL,  intent(inout) :: dXdx(:,:)  !< Partial derivatives of continuous state
!                                                                    !!   functions (X) with respect to
!                                                                    !!   the continuous states (x) [intent in to avoid deallocation]
!    real(R8Ki), allocatable, OPTIONAL,  intent(inout) :: dXddx(:,:) !< Partial derivatives of discrete state
!                                                                    !!   functions (Xd) with respect to
!                                                                    !!   the continuous states (x) [intent in to avoid deallocation]
!    real(R8Ki), allocatable, OPTIONAL,  intent(inout) :: dZdx(:,:)  !< Partial derivatives of constraint state
!                                                                    !!   functions (Z) with respect to
!                                                                    !!   the continuous states (x) [intent in to avoid deallocation]
!    integer(IntKi) :: i,j    ! Loop index
!    integer(IntKi) :: idx  ! Index of output channel in AllOuts
!    integer(IntKi) :: iDOF ! Mode number
!    ! Initialize errStat
!    errStat = ErrID_None
!    errMsg  = ''
!    if (present(dYdx)) then
!       ! allocate and set dYdx
!       if (.not. allocated(dYdx)) then
!           call AllocAry(dYdx, N_outPUTS+p%NumOuts, 2*p%nCB, 'dYdx', errStat, errMsg); if(Failed()) return
!           do i=1,size(dYdx,1); do j=1,size(dYdx,2); dYdx(i,j)=0.0_ReKi; enddo;enddo
!       end if
!       dYdx(1:6,1:2*p%nCB) = p%CMat(1:6, 1:2*p%nCB)
!       ! WriteOutputs
!       do i = 1,p%NumOuts
!           idx  = p%OutParam(i)%Indx
!           iDOF = mod(idx-ID_QSTART, p%nCB)+1
!            ! if output is an interface load dYdx is a row of the Cmatrix
!            if     (idx==ID_PtfFx) then; dYdx(6+i,1:2*p%nCB) = p%CMat(1,1:2*p%nCB)
!            elseif (idx==ID_PtfFy) then; dYdx(6+i,1:2*p%nCB) = p%CMat(2,1:2*p%nCB)
!            elseif (idx==ID_PtfFx) then; dYdx(6+i,1:2*p%nCB) = p%CMat(3,1:2*p%nCB)
!            elseif (idx==ID_PtfMx) then; dYdx(6+i,1:2*p%nCB) = p%CMat(4,1:2*p%nCB)
!            elseif (idx==ID_PtfMy) then; dYdx(6+i,1:2*p%nCB) = p%CMat(5,1:2*p%nCB)
!            elseif (idx==ID_PtfMz) then; dYdx(6+i,1:2*p%nCB) = p%CMat(6,1:2*p%nCB)
!            ! Below we look at the index, we assumed an order for the outputs
!            ! where after the index ID_Qstart, the AllOutputs are: Q,QDot and Qf
!            ! An alternative coulbe to look at the name of the DOF instead:
!            ! e.g. if (index(p%OutParam,'CBQ_')>0) then ... (see SetOutParam) 
!            else if ((idx-ID_QStart>=  0    ) .and. (idx-ID_QStart<p%nCB) ) then
!                ! Output is a DOF position, dYdx has a 1 at the proper location
!                dYdx(6+i,1:2*p%nCB   ) = 0.0_ReKi
!                dYdx(6+i,        iDOF) = 1.0_ReKi ! TODO TODO TODO ALLDOF_2_DOF
!            else if ((idx-ID_QStart>=  p%nCB) .and. (idx-ID_QStart<2*p%nCB) ) then
!                ! Output is a DOF velocity, dYdx has a 1 at the proper location
!                dYdx(6+i,1:2*p%nCB   ) = 0.0_ReKi
!                dYdx(6+i,p%nCB + iDOF) = 1.0_ReKi ! TODO TODO TODO ALLDOF_2_DOF
!            else ! e.g. WaveElevation or CB Forces
!                dYdx(6+i,1:2*p%nCB  ) = 0.0_ReKi
!            endif 
!       end do
!    end if
!    if (present(dXdx)) then
!       ! allocate and set dXdx
!       if (.not. allocated(dXdx)) then
!           call AllocAry(dXdx, 2*p%nCB, 2*p%nCB, 'dXdx', errStat, errMsg); if(Failed()) return
!           do i=1,size(dXdx,1); do j=1,size(dXdx,2); dXdx(i,j)=0.0_ReKi; enddo;enddo
!       end if
!       dXdx(1:2*p%nCB,1:2*p%nCB) = p%AMat(1:2*p%nCB,1:2*p%nCB)
!    end if
!    if (present(dXddx)) then
!    end if
!    if (present(dZdx)) then
!    end if
! CONTAinS
!     logical function Failed()
!         call SeterrStatSimple(errStat, errMsg, 'LD_JacobianPInput')
!         Failed =  errStat >= AbortErrLev
!     end function Failed
! end subroutine LD_JacobianPContState
! !----------------------------------------------------------------------------------------------------------------------------------
! !> Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) functions
! !! with respect to the discrete states (xd). The partial derivatives dY/dxd, dX/dxd, dXd/dxd, and DZ/dxd are returned.
! subroutine LD_JacobianPDiscState( t, u, p, x, xd, z, OtherState, y, m, errStat, errMsg, dYdxd, dXdxd, dXddxd, dZdxd )
! !..................................................................................................................................
! 
!    real(DbKi),                                intent(in   )           :: t          !< Time in seconds at operating point
!    type(LD_InputType),                   intent(in   )           :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
!    type(LD_ParameterType),               intent(in   )           :: p          !< Parameters
!    type(LD_ContinuousStateType),         intent(in   )           :: x          !< Continuous states at operating point
!    type(LD_DiscreteStateType),           intent(in   )           :: xd         !< Discrete states at operating point
!    type(LD_ConstraintStateType),         intent(in   )           :: z          !< Constraint states at operating point
!    type(LD_OtherStateType),              intent(in   )           :: OtherState !< Other states at operating point
!    type(LD_OutputType),                  intent(in   )           :: y          !< Output (change to inout if a mesh copy is required); 
!                                                                                     !!   Output fields are not used by this routine, but type is   
!                                                                                     !!   available here so that mesh parameter information (i.e.,  
!                                                                                     !!   connectivity) does not have to be recalculated for dYdxd.
!    type(LD_MiscVarType),                 intent(inout)           :: m          !< Misc/optimization variables
!    integer(IntKi),                            intent(  out)           :: errStat    !< Error status of the operation
!    character(*),                              intent(  out)           :: errMsg     !< Error message if errStat /= ErrID_None
!    real(R8Ki), allocatable, OPTIONAL,         intent(inout)           :: dYdxd(:,:) !< Partial derivatives of output functions
!                                                                                     !!  (Y) with respect to the discrete
!                                                                                     !!  states (xd) [intent in to avoid deallocation]
!    real(R8Ki), allocatable, OPTIONAL,         intent(inout)           :: dXdxd(:,:) !< Partial derivatives of continuous state
!                                                                                     !!   functions (X) with respect to the
!                                                                                     !!   discrete states (xd) [intent in to avoid deallocation]
!    real(R8Ki), allocatable, OPTIONAL,         intent(inout)           :: dXddxd(:,:)!< Partial derivatives of discrete state
!                                                                                     !!   functions (Xd) with respect to the
!                                                                                     !!   discrete states (xd) [intent in to avoid deallocation]
!    real(R8Ki), allocatable, OPTIONAL,         intent(inout)           :: dZdxd(:,:) !< Partial derivatives of constraint state
!                                                                                     !!   functions (Z) with respect to the
!    ! Initialize errStat
!    errStat = ErrID_None
!    errMsg  = ''
!    if (present(dYdxd)) then
!    end if
!    if (present(dXdxd)) then
!    end if
!    if (present(dXddxd)) then
!    end if
!    if (present(dZdxd)) then
!    end if
! end subroutine LD_JacobianPDiscState
! !----------------------------------------------------------------------------------------------------------------------------------
! !> Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) functions
! !! with respect to the constraint states (z). The partial derivatives dY/dz, dX/dz, dXd/dz, and DZ/dz are returned.
! subroutine LD_JacobianPConstrState( t, u, p, x, xd, z, OtherState, y, m, errStat, errMsg, dYdz, dXdz, dXddz, dZdz )
! !..................................................................................................................................
!    real(DbKi),                                intent(in   )           :: t          !< Time in seconds at operating point
!    type(LD_InputType),                   intent(in   )           :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
!    type(LD_ParameterType),               intent(in   )           :: p          !< Parameters
!    type(LD_ContinuousStateType),         intent(in   )           :: x          !< Continuous states at operating point
!    type(LD_DiscreteStateType),           intent(in   )           :: xd         !< Discrete states at operating point
!    type(LD_ConstraintStateType),         intent(in   )           :: z          !< Constraint states at operating point
!    type(LD_OtherStateType),              intent(in   )           :: OtherState !< Other states at operating point
!    type(LD_OutputType),                  intent(in   )           :: y          !< Output (change to inout if a mesh copy is required); 
!                                                                                     !!   Output fields are not used by this routine, but type is   
!                                                                                     !!   available here so that mesh parameter information (i.e.,  
!                                                                                     !!   connectivity) does not have to be recalculated for dYdz.
!    type(LD_MiscVarType),                 intent(inout)           :: m          !< Misc/optimization variables
!    integer(IntKi),                            intent(  out)           :: errStat    !< Error status of the operation
!    character(*),                              intent(  out)           :: errMsg     !< Error message if errStat /= ErrID_None
!    real(R8Ki), allocatable, OPTIONAL,         intent(inout)           :: dYdz(:,:)  !< Partial derivatives of output
!                                                                                     !!  functions (Y) with respect to the
!                                                                                     !!  constraint states (z) [intent in to avoid deallocation]
!    real(R8Ki), allocatable, OPTIONAL,         intent(inout)           :: dXdz(:,:)  !< Partial derivatives of continuous
!                                                                                     !!  state functions (X) with respect to
!                                                                                     !!  the constraint states (z) [intent in to avoid deallocation]
!    real(R8Ki), allocatable, OPTIONAL,         intent(inout)           :: dXddz(:,:) !< Partial derivatives of discrete state
!                                                                                     !!  functions (Xd) with respect to the
!                                                                                     !!  constraint states (z) [intent in to avoid deallocation]
!    real(R8Ki), allocatable, OPTIONAL,         intent(inout)           :: dZdz(:,:)  !< Partial derivatives of constraint
!    ! Initialize errStat
!    errStat = ErrID_None
!    errMsg  = ''
!    if (present(dYdz)) then
!    end if
!    if (present(dXdz)) then
!    end if
!    if (present(dXddz)) then
!    end if
!    if (present(dZdz)) then
!    end if
! end subroutine LD_JacobianPConstrState
! !----------------------------------------------------------------------------------------------------------------------------------
! !> Routine to pack the data structures representing the operating points into arrays for linearization.
! subroutine LD_GetOP( t, u, p, x, xd, z, OtherState, y, m, errStat, errMsg, u_op, y_op, x_op, dx_op, xd_op, z_op )
!    real(DbKi),                           intent(in   )           :: t          !< Time in seconds at operating point
!    type(LD_InputType),              intent(in   )           :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
!    type(LD_ParameterType),          intent(in   )           :: p          !< Parameters
!    type(LD_ContinuousStateType),    intent(in   )           :: x          !< Continuous states at operating point
!    type(LD_DiscreteStateType),      intent(in   )           :: xd         !< Discrete states at operating point
!    type(LD_ConstraintStateType),    intent(in   )           :: z          !< Constraint states at operating point
!    type(LD_OtherStateType),         intent(in   )           :: OtherState !< Other states at operating point
!    type(LD_OutputType),             intent(in   )           :: y          !< Output at operating point
!    type(LD_MiscVarType),            intent(inout)           :: m          !< Misc/optimization variables
!    integer(IntKi),                       intent(  out)           :: errStat    !< Error status of the operation
!    character(*),                         intent(  out)           :: errMsg     !< Error message if errStat /= ErrID_None
!    real(ReKi), allocatable, OPTIONAL,    intent(inout)           :: u_op(:)    !< values of linearized inputs
!    real(ReKi), allocatable, OPTIONAL,    intent(inout)           :: y_op(:)    !< values of linearized outputs
!    real(ReKi), allocatable, OPTIONAL,    intent(inout)           :: x_op(:)    !< values of linearized continuous states
!    real(ReKi), allocatable, OPTIONAL,    intent(inout)           :: dx_op(:)   !< values of first time derivatives of linearized continuous states
!    real(ReKi), allocatable, OPTIONAL,    intent(inout)           :: xd_op(:)   !< values of linearized discrete states
!    real(ReKi), allocatable, OPTIONAL,    intent(inout)           :: z_op(:)    !< values of linearized constraint states
!    integer(IntKi)                    :: I
!    type(LD_ContinuousStateType) :: dx          !< derivative of continuous states at operating point
!    ! Initialize errStat
!    errStat = ErrID_None
!    errMsg  = ''
! 
!    if ( present( u_op ) ) then
!        if (.not. allocated(u_op)) then
!            call AllocAry(u_op, N_inPUTS, 'u_op', errStat, errMsg); if(Failed())return
!        endif
!        u_op(1:3)   = u%PtfmMesh%TranslationDisp(:,1)
!        u_op(4:6)   = GetSmllRotAngs(u%PtfmMesh%Orientation(:,:,1), errStat, errMsg); if(Failed())return
!        u_op(7:9  ) = u%PtfmMesh%TranslationVel(:,1)
!        u_op(10:12) = u%PtfmMesh%RotationVel   (:,1)
!        u_op(13:15) = u%PtfmMesh%TranslationAcc(:,1)
!        u_op(16:18) = u%PtfmMesh%RotationAcc   (:,1)
!    end if
! 
!    if ( present( y_op ) ) then
!        if (.not. allocated(y_op)) then
!            call AllocAry(y_op, N_outPUTS+p%NumOuts, 'y_op', errStat, errMsg); if(Failed())return
!        endif
!        ! Update the output mesh
!        y_op(1:3)=y%PtfmMesh%Force(1:3,1)
!        y_op(4:6)=y%PtfmMesh%Moment(1:3,1)
!        do i=1,p%NumOuts         
!            y_op(i+N_outPUTS) = y%WriteOutput(i)
!        end do      
!    end if
! 
!    if ( present( x_op ) ) then
!        if (.not. allocated(x_op)) then
!            call AllocAry(x_op, 2*p%nCB, 'x_op', errStat, errMsg); if (Failed())return
!        endif
!        x_op(1:p%nCB)         = x%qm(1:p%nCB)
!        x_op(p%nCB+1:2*p%nCB) = x%qmdot(1:p%nCB)
!    end if
! 
!    if ( present( dx_op ) ) then
!        if (.not. allocated(dx_op)) then
!            call AllocAry(dx_op, 2*p%nCB, 'dx_op', errStat, errMsg); if (Failed())return
!        endif
!        call LD_CalcContStateDeriv(t, u, p, x, xd, z, OtherState, m, dx, errStat, errMsg); if(Failed()) return
!        dx_op(1:p%nCB)         = dx%qm(1:p%nCB)
!        dx_op(p%nCB+1:2*p%nCB) = dx%qmdot(1:p%nCB)
!    end if
! 
!    if ( present( xd_op ) ) then
!    end if
!    
!    if ( present( z_op ) ) then
!    end if
! 
! contains
!     logical function Failed()
!         call SeterrStatSimple(errStat, errMsg, 'LD_GetOP')
!         Failed =  errStat >= AbortErrLev
!     end function Failed
! end subroutine LD_GetOP
! !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end module LinDyn
!**********************************************************************************************************************************
