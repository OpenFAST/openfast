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

subroutine LD_Init(InitInp, u, p, x, xd, z, OtherState, y, m, dt_gluecode, InitOut, errStat, errMsg)
   type(LD_InitInputType),       intent(in   ) :: InitInp     !< Input data for initialization routine
   type(LD_InputType),           intent(out)   :: u           !< An initial guess for the input; input mesh must be defined
   type(LD_ParameterType),       intent(out)   :: p           !< Parameters
   type(LD_ContinuousStateType), intent(out)   :: x           !< Initial continuous states
   type(LD_DiscreteStateType),   intent(out)   :: xd          !< Initial discrete states
   type(LD_ConstraintStateType), intent(out)   :: z           !< Initial guess of the constraint states
   type(LD_OtherStateType),      intent(out)   :: OtherState  !< Initial other states (logical, etc)
   type(LD_OutputType),          intent(out)   :: y           !< Initial system outputs (outputs are not calculated;
   type(LD_MiscVarType),         intent(out)   :: m           !< Misc variables for optimization (not copied in glue code)
   real(DbKi),                   intent(inout) :: dt_gluecode !< Coupling interval in seconds: the rate that
   type(LD_InitOutputType),      intent(out)   :: InitOut     !< Output for initialization routine
   integer(IntKi),               intent(out)   :: errStat     !< Error status of the operation
   character(*),                 intent(out)   :: errMsg      !< Error message if errStat /= ErrID_None
   integer(IntKi)  :: errStat2    ! Status of error message
   character(1024) :: errMsg2     ! Error message if ErrStat /= ErrID_None
   integer(IntKi)  :: i, n ! Loop counter
   ! Misc Init
   errStat = ErrID_None
   errMsg  = ""
   call NWTC_Init( )      ! Initialize the NWTC Subroutine Library
   call DispNVD( LD_Ver ) ! Display the module information

   ! --- Setting Params from InitInp
   n = size(p%MM,1)
   call AllocAry(p%MM        , n, n, 'MM', errStat2, errMsg2); if(Failed()) return
   call AllocAry(p%CC        , n, n, 'CC', errStat2, errMsg2); if(Failed()) return
   call AllocAry(p%KK        , n, n, 'KK', errStat2, errMsg2); if(Failed()) return
   call AllocAry(p%activeDOFs, n   , 'activeDOFs', errStat2, errMsg2); if(Failed()) return
   p%dt         = InitInp%dt
   p%IntMethod  = InitInp%IntMethod
   p%MM         = InitInp%MM
   p%CC         = InitInp%CC
   p%KK         = InitInp%KK
   p%activeDOFs = InitInp%activeDOFs

!    INTERFACE LAPACK_getri 
!    SUBROUTINE LAPACK_DGETRI( N, A, IPIV, WORK, LWORK, ErrStat, ErrMsg )


!    ! Setting p%OutParam from OutList
!    call SetOutParam(InputFileData%OutList, InputFileData%NumOuts, p, errStat, errMsg); if(Failed()) return
!    ! Set the constant state matrices A,B,C,D
!    call SetStateMatrices(p, errStat, errMsg)
!   
!    ! --- Allocate and init continuous states
!    call AllocAry( x%qm    , p%nCB,'CB DOF positions' , errStat,errMsg); if(Failed()) return
!    call AllocAry( x%qmdot , p%nCB,'CB DOF velocities', errStat,errMsg); if(Failed()) return
!    if (allocated(InputFileData%InitPosList)) then
!        if (size(InputFileData%InitPosList)/=p%nCB) then
!            call SeterrStat(ErrID_Fatal, 'The number of elements of `InitPosList` ('//trim(Num2LStr(size(InputFileData%InitPosList)))//') does not match the number of CB modes: '//trim(Num2LStr(p%nCB)), errStat, errMsg, 'LD_Init'); 
!            return
!        endif
!        do I=1,p%nCB;
!            x%qm(I)=InputFileData%InitPosList(I);
!        end do
!    else
!        do I=1,p%nCB; x%qm   (I)=0; end do
!    endif
!    if (allocated(InputFileData%InitVelList)) then
!        if (size(InputFileData%InitVelList)/=p%nCB) then
!            call SeterrStat(ErrID_Fatal, 'The number of elements of `InitVelList` ('//trim(Num2LStr(size(InputFileData%InitVelList)))//') does not match the number of CB modes: '//trim(Num2LStr(p%nCB)), errStat, errMsg, 'LD_Init'); 
!            return
!        endif
!        do I=1,p%nCB;
!            x%qmdot(I)=InputFileData%InitVelList(I);
!        enddo
!    else
!        do I=1,p%nCB; x%qmdot(I)=0; end do
!    endif
! 
!    ! Other states
!    xd%DummyDiscState          = 0.0_ReKi
!    z%DummyConstrState         = 0.0_ReKi
!    ! allocate OtherState%xdot if using multi-step method; initialize n
!    if ( ( p%IntMethod .eq. 2) .OR. ( p%IntMethod .eq. 3)) THEN
!        allocate( OtherState%xdot(4), STAT=errStat )
!        errMsg='Error allocating OtherState%xdot'
!        if(Failed()) return
!    endif
! 
!    ! Initialize Misc Variables:
!    !m%EquilStart = InputFileData%EquilStart
!    m%EquilStart = .False. ! Feature not yet implemented
! 
!    m%Indx = 1 ! used to optimize interpolation of loads in time
!    call AllocAry( m%F_at_t, p%nTot,'Loads at t', errStat,errMsg); if(Failed()) return
!    do I=1,p%nTot; m%F_at_t(I)=0; end do
!    call AllocAry( m%xFlat, 2*p%nCB,'xFlat', errStat,errMsg); if(Failed()) return
!    do I=1,2*p%nCB; m%xFlat(I)=0; end do
!    do I=1,N_inPUTS; m%uFlat(I)=0; end do
!    
!    ! Define initial guess (set up mesh first) for the system inputs here:
!    call Init_meshes(u, y, InitInp, errStat, errMsg); if(Failed()) return
! 
!    ! --- Outputs
!    call AllocAry( m%AllOuts, ID_QStart+3*p%nCBFull-1, "LinDyn AllOut", errStat,errMsg ); if(Failed()) return
!    m%AllOuts(1:ID_QStart+3*p%nCBFull-1) = 0.0
!    call AllocAry( y%WriteOutput,        p%NumOuts,'WriteOutput',   errStat,errMsg); if(Failed()) return
!    call AllocAry(InitOut%WriteOutputHdr,p%NumOuts,'WriteOutputHdr',errStat,errMsg); if(Failed()) return
!    call AllocAry(InitOut%WriteOutputUnt,p%NumOuts,'WriteOutputUnt',errStat,errMsg); if(Failed()) return
!    y%WriteOutput(1:p%NumOuts) = 0.0
!    InitOut%WriteOutputHdr(1:p%NumOuts) = p%OutParam(1:p%NumOuts)%Name
!    InitOut%WriteOutputUnt(1:p%NumOuts) = p%OutParam(1:p%NumOuts)%Units     
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
!       ! 
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
        call SeterrStatSimple(errStat, errMsg, 'LD_Init')
        Failed =  errStat >= AbortErrLev
    end function Failed
end subroutine LD_Init
! 
! 
!----------------------------------------------------------------------------------------------------------------------------------
! subroutine SetStateMatrices( p, errStat, errMsg)
! subroutine SetStateMatrices( p, errStat, errMsg)
! !..................................................................................................................................
!    type(LD_ParameterType), intent(inout) :: p                                   !< All the parameter matrices stored in this input file
!    integer(IntKi),              intent(out)   :: errStat                             !< Error status                              
!    character(*),                intent(out)   :: errMsg                              !< Error message
!    ! Local variables:
!    integer(IntKi)                          :: I                                         ! loop counter
!    integer(IntKi)                          :: nX                                        ! Number of states
!    integer(IntKi)                          :: nU                                        ! Number of inputs
!    integer(IntKi)                          :: nY                                        ! Number of ouputs
!    integer(IntKi)                          :: n1                                        ! Number of interface DOF
!    integer(IntKi)                          :: n2                                        ! Number of CB DOF
!    real(ReKi), dimension(:,:), allocatable :: I22
!    ! Init 
!    nX = 2*p%nCB
!    nU = 3*6
!    nY = 6
!    n1 = 6
!    n2 = p%nCB
!    if (allocated(p%AMat)) deallocate(p%AMat)
!    if (allocated(p%BMat)) deallocate(p%BMat)
!    if (allocated(p%CMat)) deallocate(p%CMat)
!    if (allocated(p%DMat)) deallocate(p%DMat)
!    if (allocated(p%M11))  deallocate(p%M11)
!    if (allocated(p%M12))  deallocate(p%M12)
!    if (allocated(p%M22))  deallocate(p%M22)
!    if (allocated(p%M21))  deallocate(p%M21)
!    if (allocated(p%C11))  deallocate(p%C11)
!    if (allocated(p%C12))  deallocate(p%C12)
!    if (allocated(p%C22))  deallocate(p%C22)
!    if (allocated(p%C21))  deallocate(p%C21)
!    if (allocated(p%K11))  deallocate(p%C11)
!    if (allocated(p%K22))  deallocate(p%C22)
!    ! Allocation
!    call allocAry(p%AMat, nX, nX, 'p%AMat', errStat, errMsg); if(Failed()) return ; p%AMat(1:nX,1:nX) =0
!    call allocAry(p%BMat, nX, nU, 'p%BMat', errStat, errMsg); if(Failed()) return ; p%BMat(1:nX,1:nU) =0
!    call allocAry(p%FX  , nX,     'p%FX'  , errStat, errMsg); if(Failed()) return ; p%Fx  (1:nX)      =0
!    call allocAry(p%CMat, nY, nX, 'p%CMat', errStat, errMsg); if(Failed()) return ; p%CMat(1:nY,1:nX) =0
!    call allocAry(p%DMat, nY, nU, 'p%DMat', errStat, errMsg); if(Failed()) return ; p%DMat(1:nY,1:nU) =0
!    call allocAry(p%FY  , nY,     'p%FY'  , errStat, errMsg); if(Failed()) return ; p%FY  (1:nY)      =0
!    call allocAry(p%M11 , n1, n1, 'p%M11' , errStat, errMsg); if(Failed()) return ; p%M11 (1:n1,1:n1) =0
!    call allocAry(p%K11 , n1, n1, 'p%K11' , errStat, errMsg); if(Failed()) return ; p%K11 (1:n1,1:n1) =0
!    call allocAry(p%C11 , n1, n1, 'p%C11' , errStat, errMsg); if(Failed()) return ; p%C11 (1:n1,1:n1) =0
!    call allocAry(p%M22 , n2, n2, 'p%M22' , errStat, errMsg); if(Failed()) return ; p%M22 (1:n2,1:n2) =0
!    call allocAry(p%K22 , n2, n2, 'p%K22' , errStat, errMsg); if(Failed()) return ; p%K22 (1:n2,1:n2) =0
!    call allocAry(p%C22 , n2, n2, 'p%C22' , errStat, errMsg); if(Failed()) return ; p%C22 (1:n2,1:n2) =0
!    call allocAry(p%M12 , n1, n2, 'p%M12' , errStat, errMsg); if(Failed()) return ; p%M12 (1:n1,1:n2) =0
!    call allocAry(p%C12 , n1, n2, 'p%C12' , errStat, errMsg); if(Failed()) return ; p%C12 (1:n1,1:n2) =0
!    call allocAry(p%M21 , n2, n1, 'p%M21' , errStat, errMsg); if(Failed()) return ; p%M21 (1:n2,1:n1) =0
!    call allocAry(p%C21 , n2, n1, 'p%C21' , errStat, errMsg); if(Failed()) return ; p%C21 (1:n2,1:n1) =0
!    call allocAry(  I22 , n2, n2, '  I22' , errStat, errMsg); if(Failed()) return ;   I22 (1:n2,1:n2) =0
!    do I=1,n2 ; I22(I,I)=1; enddo ! Identity matrix
!    ! Submatrices
!    p%M11(1:n1,1:n1) = p%Mass(1:n1      ,1:n1      )
!    p%C11(1:n1,1:n1) = p%Damp(1:n1      ,1:n1      )
!    p%K11(1:n1,1:n1) = p%Stff(1:n1      ,1:n1      )
!    p%M12(1:n1,1:n2) = p%Mass(1:n1      ,n1+1:n1+n2)
!    p%C12(1:n1,1:n2) = p%Damp(1:n1      ,n1+1:n1+n2)
!    p%M21(1:n2,1:n1) = p%Mass(n1+1:n1+n2,1:n1      )
!    p%C21(1:n2,1:n1) = p%Damp(n1+1:n1+n2,1:n1      )
!    p%M22(1:n2,1:n2) = p%Mass(n1+1:n1+n2,n1+1:n1+n2)
!    p%C22(1:n2,1:n2) = p%Damp(n1+1:n1+n2,n1+1:n1+n2)
!    p%K22(1:n2,1:n2) = p%Stff(n1+1:n1+n2,n1+1:n1+n2)
!    ! A matrix
!    p%AMat(1:n2   ,n2+1:nX) = I22   (1:n2,1:n2)
!    p%AMat(n2+1:nX,1:n2   ) = -p%K22(1:n2,1:n2)
!    p%AMat(n2+1:nX,n2+1:nX) = -p%C22(1:n2,1:n2)
!    ! B matrix
!    p%BMat(n2+1:nX,7 :12  ) = -p%C21(1:n2,1:6)
!    p%BMat(n2+1:nX,13:18  ) = -p%M21(1:n2,1:6)
!    ! C matrix
!    p%CMat(1:nY,1:n2   ) = matmul(p%M12,p%K22)
!    p%CMat(1:nY,n2+1:nX) = matmul(p%M12,p%C22) - p%C12
!    ! D matrix
!    p%DMat(1:nY,1:6   ) = -p%K11
!    p%DMat(1:nY,7:12  ) = -p%C11 + matmul(p%M12,p%C21)
!    p%DMat(1:nY,13:18 ) = -p%M11 + matmul(p%M12,p%M21)
! CONTAinS
!     logical function Failed()
!         call SeterrStatSimple(errStat, errMsg, 'LD_SetStateMatrices')
!         Failed =  errStat >= AbortErrLev
!     end function Failed
! end subroutine SetStateMatrices
! !----------------------------------------------------------------------------------------------------------------------------------
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
! !----------------------------------------------------------------------------------------------------------------------------------
! !> This routine is called at the end of the simulation.
! subroutine LD_End( u, p, x, xd, z, OtherState, y, m, errStat, errMsg )
! !..................................................................................................................................
!    type(LD_InputType),           intent(inout)  :: u           !< System inputs
!    type(LD_ParameterType),       intent(inout)  :: p           !< Parameters
!    type(LD_ContinuousStateType), intent(inout)  :: x           !< Continuous states
!    type(LD_DiscreteStateType),   intent(inout)  :: xd          !< Discrete states
!    type(LD_ConstraintStateType), intent(inout)  :: z           !< Constraint states
!    type(LD_OtherStateType),      intent(inout)  :: OtherState  !< Other states
!    type(LD_OutputType),          intent(inout)  :: y           !< System outputs
!    type(LD_MiscVarType),         intent(inout)  :: m           !< Misc variables for optimization (not copied in glue code)
!    integer(IntKi),                    intent(  out)  :: errStat     !< Error status of the operation
!    character(*),                      intent(  out)  :: errMsg      !< Error message if errStat /= ErrID_None
!    ! Place any last minute operations or calculations here:
!    ! Close files here (but because of checkpoint-restart capability, it is not recommended to have files open during the simulation):
!    ! Destroy the input data:
!    call LD_DestroyInput( u, errStat, errMsg ); if(Failed()) return
!    ! Destroy the parameter data:
!    call LD_DestroyParam( p, errStat, errMsg ); if(Failed()) return
!    ! Destroy the state data:
!    call LD_DestroyContState(   x,          errStat,errMsg); if(Failed()) return
!    call LD_DestroyDiscState(   xd,         errStat,errMsg); if(Failed()) return
!    call LD_DestroyConstrState( z,          errStat,errMsg); if(Failed()) return
!    call LD_DestroyOtherState(  OtherState, errStat,errMsg); if(Failed()) return
!    ! Destroy the output data:
!    call LD_DestroyOutput( y, errStat, errMsg ); if(Failed()) return
!    ! Destroy the misc data:
!    call LD_DestroyMisc( m, errStat, errMsg ); if(Failed()) return
! CONTAinS
!     logical function Failed()
!         call SeterrStatSimple(errStat, errMsg, 'LD_End')
!         Failed =  errStat >= AbortErrLev
!     end function Failed
! end subroutine LD_End
! 
! 
! !----------------------------------------------------------------------------------------------------------------------------------
! !> This subroutine implements the fourth-order Adams-Bashforth Method (RK4) for numerically integrating ordinary differential 
! !! equations:
! !!  Let f(t, x) = xdot denote the time (t) derivative of the continuous states (x). 
! !!      x(t+dt) = x(t)  + (dt / 24.) * ( 55.*f(t,x) - 59.*f(t-dt,x) + 37.*f(t-2.*dt,x) - 9.*f(t-3.*dt,x) )
! !!  See, e.g.,
! !!      http://en.wikipedia.org/wiki/Linear_multistep_method
! !!      K. E. Atkinson, "An Introduction to Numerical Analysis", 1989, John Wiley & Sons, Inc, Second Edition.
! subroutine LD_AB4( t, n, u, utimes, p, x, xd, z, OtherState, m, errStat, errMsg )
! !..................................................................................................................................
!    real(DbKi),                        intent(in   )  :: t           !< Current simulation time in seconds
!    integer(IntKi),                    intent(in   )  :: n           !< time step number
!    type(LD_InputType),           intent(inout)  :: u(:)        !< Inputs at t
!    real(DbKi),                        intent(in   )  :: utimes(:)   !< times of input
!    type(LD_ParameterType),       intent(in   )  :: p           !< Parameters
!    type(LD_ContinuousStateType), intent(inout)  :: x           !< Continuous states at t on input at t + dt on output
!    type(LD_DiscreteStateType),   intent(in   )  :: xd          !< Discrete states at t
!    type(LD_ConstraintStateType), intent(in   )  :: z           !< Constraint states at t (possibly a guess)
!    type(LD_OtherStateType),      intent(inout)  :: OtherState  !< Other states at t on input at t + dt on output
!    type(LD_MiscVarType),         intent(inout)  :: m           !< Misc/optimization variables
!    integer(IntKi),                    intent(  out)  :: errStat     !< Error status of the operation
!    character(*),                      intent(  out)  :: errMsg      !< Error message if errStat /= ErrID_None
!    ! local variables
!    type(LD_ContinuousStateType) :: xdot       ! Continuous state derivs at t
!    type(LD_InputType)           :: u_interp
!    ! Initialize errStat
!    errStat = ErrID_None
!    errMsg  = "" 
!    
!    ! need xdot at t
!    call LD_CopyInput(u(1), u_interp, MESH_NEWCOPY, errStat, errMsg  )  ! we need to allocate input arrays/meshes before calling ExtrapInterp...
!    call LD_Input_ExtrapInterp(u, utimes, u_interp, t, errStat, errMsg)
!    call LD_CalcContStateDeriv( t, u_interp, p, x, xd, z, OtherState, m, xdot, errStat, errMsg ) ! initializes xdot
!    call LD_DestroyInput( u_interp, errStat, errMsg)   ! we don't need this local copy anymore
!    if (n .le. 2) then
!       OtherState%n = n
!       call LD_CopyContState(xdot, OtherState%xdot(3-n), MESH_UPDATECOPY, errStat, errMsg )
!       call LD_RK4(t, n, u, utimes, p, x, xd, z, OtherState, m, errStat, errMsg )
!    else
!       if (OtherState%n .lt. n) then
!          OtherState%n = n
!          call LD_CopyContState(OtherState%xdot(3), OtherState%xdot(4), MESH_UPDATECOPY, errStat, errMsg )
!          call LD_CopyContState(OtherState%xdot(2), OtherState%xdot(3), MESH_UPDATECOPY, errStat, errMsg )
!          call LD_CopyContState(OtherState%xdot(1), OtherState%xdot(2), MESH_UPDATECOPY, errStat, errMsg )
!       elseif (OtherState%n .gt. n) then
!          errStat = ErrID_Fatal
!          errMsg = ' Backing up in time is not supported with a multistep method '
!          RETURN
!       endif
!       call LD_CopyContState( xdot, OtherState%xdot ( 1 ), MESH_UPDATECOPY, errStat, errMsg )
!       !OtherState%xdot ( 1 )     = xdot  ! make sure this is most up to date
!       x%qm    = x%qm    + (p%EP_DeltaT / 24.) * ( 55.*OtherState%xdot(1)%qm - 59.*OtherState%xdot(2)%qm    + 37.*OtherState%xdot(3)%qm  &
!                                     - 9. * OtherState%xdot(4)%qm )
!       x%qmdot = x%qmdot + (p%EP_DeltaT / 24.) * ( 55.*OtherState%xdot(1)%qmdot - 59.*OtherState%xdot(2)%qmdot  &
!                                        + 37.*OtherState%xdot(3)%qmdot  - 9.*OtherState%xdot(4)%qmdot )
!    endif
!    call LD_DestroyContState(xdot, errStat, errMsg)
!    call LD_DestroyInput(u_interp, errStat, errMsg)
!    
! end subroutine LD_AB4
! !----------------------------------------------------------------------------------------------------------------------------------
! !> This subroutine implements the fourth-order Adams-Bashforth-Moulton Method (RK4) for numerically integrating ordinary 
! !! differential equations:
! !!   Let f(t, x) = xdot denote the time (t) derivative of the continuous states (x). 
! !!   Adams-Bashforth Predictor:
! !!      x^p(t+dt) = x(t)  + (dt / 24.) * ( 55.*f(t,x) - 59.*f(t-dt,x) + 37.*f(t-2.*dt,x) - 9.*f(t-3.*dt,x) )
! !!   Adams-Moulton Corrector:
! !!      x(t+dt) = x(t)  + (dt / 24.) * ( 9.*f(t+dt,x^p) + 19.*f(t,x) - 5.*f(t-dt,x) + 1.*f(t-2.*dt,x) )
! !!  See, e.g.,
! !!      http://en.wikipedia.org/wiki/Linear_multistep_method
! !!      K. E. Atkinson, "An Introduction to Numerical Analysis", 1989, John Wiley & Sons, Inc, Second Edition.
! subroutine LD_ABM4( t, n, u, utimes, p, x, xd, z, OtherState, m, errStat, errMsg )
! !..................................................................................................................................
!    real(DbKi),                          intent(in   ) :: t           !< Current simulation time in seconds
!    integer(IntKi),                      intent(in   ) :: n           !< time step number
!    type(LD_InputType),             intent(inout) :: u(:)        !< Inputs at t
!    real(DbKi),                          intent(in   ) :: utimes(:)   !< times of input
!    type(LD_ParameterType),         intent(in   ) :: p           !< Parameters
!    type(LD_ContinuousStateType),   intent(inout) :: x           !< Continuous states at t on input at t + dt on output ! TODO TODO TODO in
!    type(LD_DiscreteStateType),     intent(in   ) :: xd          !< Discrete states at t
!    type(LD_ConstraintStateType),   intent(in   ) :: z           !< Constraint states at t (possibly a guess)
!    type(LD_OtherStateType),        intent(inout) :: OtherState  !< Other states at t on input at t + dt on output
!    type(LD_MiscVarType),           intent(inout) :: m           !< Misc/optimization variables
!    integer(IntKi),                      intent(  out) :: errStat     !< Error status of the operation
!    character(*),                        intent(  out) :: errMsg      !< Error message if errStat /= ErrID_None
!    ! local variables
!    type(LD_InputType)            :: u_interp        ! Continuous states at t
!    type(LD_ContinuousStateType)  :: x_pred          ! Continuous states at t
!    type(LD_ContinuousStateType)  :: xdot_pred       ! Continuous states at t
!    
!    ! Initialize errStat
!    errStat = ErrID_None
!    errMsg  = "" 
!    
!    call LD_CopyContState(x, x_pred, MESH_NEWCOPY, errStat, errMsg) !initialize x_pred      
!    call LD_AB4( t, n, u, utimes, p, x_pred, xd, z, OtherState, m, errStat, errMsg )
!    if (n .gt. 2) then
!       call LD_CopyInput( u(1), u_interp, MESH_NEWCOPY, errStat, errMsg) ! make copy so that arrays/meshes get initialized/allocated for ExtrapInterp
!       call LD_Input_ExtrapInterp(u, utimes, u_interp, t + p%EP_DeltaT, errStat, errMsg)
!       call LD_CalcContStateDeriv(t + p%EP_DeltaT, u_interp, p, x_pred, xd, z, OtherState, m, xdot_pred, errStat, errMsg ) ! initializes xdot_pred
!       call LD_DestroyInput( u_interp, errStat, errMsg) ! local copy no longer needed
!    
!       x%qm    = x%qm    + (p%EP_DeltaT / 24.) * ( 9. * xdot_pred%qm +  19. * OtherState%xdot(1)%qm - 5. * OtherState%xdot(2)%qm &
!                                        + 1. * OtherState%xdot(3)%qm )
!    
!       x%qmdot = x%qmdot + (p%EP_DeltaT / 24.) * ( 9. * xdot_pred%qmdot + 19. * OtherState%xdot(1)%qmdot - 5. * OtherState%xdot(2)%qmdot &
!                                        + 1. * OtherState%xdot(3)%qmdot )
!       call LD_DestroyContState( xdot_pred, errStat, errMsg) ! local copy no longer needed
!    else
!       x%qm    = x_pred%qm
!       x%qmdot = x_pred%qmdot
!    endif
!    call LD_DestroyContState( x_pred, errStat, errMsg) ! local copy no longer needed
! end subroutine LD_ABM4
! 
! !----------------------------------------------------------------------------------------------------------------------------------
! !> This subroutine implements the fourth-order Runge-Kutta Method (RK4) for numerically integrating ordinary differential equations:
! !!   Let f(t, x) = xdot denote the time (t) derivative of the continuous states (x). 
! !!   Define constants k1, k2, k3, and k4 as 
! !!        k1 = dt * f(t        , x_t        )
! !!        k2 = dt * f(t + dt/2 , x_t + k1/2 )
! !!        k3 = dt * f(t + dt/2 , x_t + k2/2 ), and
! !!        k4 = dt * f(t + dt   , x_t + k3   ).
! !!   Then the continuous states at t = t + dt are
! !!        x_(t+dt) = x_t + k1/6 + k2/3 + k3/3 + k4/6 + O(dt^5)
! !! For details, see:
! !!   Press, W. H.; Flannery, B. P.; Teukolsky, S. A.; and Vetterling, W. T. "Runge-Kutta Method" and "Adaptive Step Size Control for 
! !!   Runge-Kutta." sections 16.1 and 16.2 in Numerical Recipes in FORTRAN: The Art of Scientific Computing, 2nd ed. Cambridge, England: 
! !!   Cambridge University Press, pp. 704-716, 1992.
! subroutine LD_RK4( t, n, u, utimes, p, x, xd, z, OtherState, m, errStat, errMsg )
! !..................................................................................................................................
!    real(DbKi),                     intent(in   )      :: t           !< Current simulation time in seconds
!    integer(IntKi),                 intent(in   )      :: n           !< time step number
!    type(LD_InputType),             intent(inout) :: u(:)        !< Inputs at t
!    real(DbKi),                     intent(in   )      :: utimes(:)   !< times of input
!    type(LD_ParameterType),         intent(in   ) :: p           !< Parameters
!    type(LD_ContinuousStateType),   intent(inout) :: x           !< Continuous states at t on input at t + dt on output
!    type(LD_DiscreteStateType),     intent(in   ) :: xd          !< Discrete states at t
!    type(LD_ConstraintStateType),   intent(in   ) :: z           !< Constraint states at t (possibly a guess)
!    type(LD_OtherStateType),        intent(inout) :: OtherState  !< Other states at t on input at t + dt on output
!    type(LD_MiscVarType),           intent(inout) :: m           !< Misc/optimization variables
!    integer(IntKi),                 intent(  out)      :: errStat     !< Error status of the operation
!    character(*),                   intent(  out)      :: errMsg      !< Error message if errStat /= ErrID_None
!    ! local variables
!    type(LD_ContinuousStateType)                 :: xdot        ! time derivatives of continuous states      
!    type(LD_ContinuousStateType)                 :: k1          ! RK4 constant; see above
!    type(LD_ContinuousStateType)                 :: k2          ! RK4 constant; see above 
!    type(LD_ContinuousStateType)                 :: k3          ! RK4 constant; see above 
!    type(LD_ContinuousStateType)                 :: k4          ! RK4 constant; see above 
!    type(LD_ContinuousStateType)                 :: x_tmp       ! Holds temporary modification to x
!    type(LD_InputType)                           :: u_interp    ! interpolated value of inputs 
!    ! Initialize errStat
!    errStat = ErrID_None
!    errMsg  = "" 
!    
!    ! Initialize interim vars
!    !bjj: the state type contains allocatable arrays, so we must first allocate space:
!    call LD_CopyContState( x, k1,       MESH_NEWCOPY, errStat, errMsg )
!    call LD_CopyContState( x, k2,       MESH_NEWCOPY, errStat, errMsg )
!    call LD_CopyContState( x, k3,       MESH_NEWCOPY, errStat, errMsg )
!    call LD_CopyContState( x, k4,       MESH_NEWCOPY, errStat, errMsg )
!    call LD_CopyContState( x, x_tmp,    MESH_NEWCOPY, errStat, errMsg )
!    
!    ! interpolate u to find u_interp = u(t)
!    call LD_CopyInput(u(1), u_interp, MESH_NEWCOPY, errStat, errMsg  )  ! we need to allocate input arrays/meshes before calling ExtrapInterp...     
!    call LD_Input_ExtrapInterp( u, utimes, u_interp, t, errStat, errMsg )
!    
!    ! find xdot at t
!    call LD_CalcContStateDeriv( t, u_interp, p, x, xd, z, OtherState, m, xdot, errStat, errMsg ) !initializes xdot
!    
!    k1%qm    = p%EP_DeltaT * xdot%qm
!    k1%qmdot = p%EP_DeltaT * xdot%qmdot
!    x_tmp%qm    = x%qm    + 0.5 * k1%qm
!    x_tmp%qmdot = x%qmdot + 0.5 * k1%qmdot
!    
!    ! interpolate u to find u_interp = u(t + dt/2)
!    call LD_Input_ExtrapInterp(u, utimes, u_interp, t+0.5*p%EP_DeltaT, errStat, errMsg)
!    
!    ! find xdot at t + dt/2
!    call LD_CalcContStateDeriv( t + 0.5*p%EP_DeltaT, u_interp, p, x_tmp, xd, z, OtherState, m, xdot, errStat, errMsg )
!    
!    k2%qm    = p%EP_DeltaT * xdot%qm
!    k2%qmdot = p%EP_DeltaT * xdot%qmdot
!    x_tmp%qm    = x%qm    + 0.5 * k2%qm
!    x_tmp%qmdot = x%qmdot + 0.5 * k2%qmdot
!    
!    ! find xdot at t + dt/2
!    call LD_CalcContStateDeriv( t + 0.5*p%EP_DeltaT, u_interp, p, x_tmp, xd, z, OtherState, m, xdot, errStat, errMsg )
!    
!    k3%qm    = p%EP_DeltaT * xdot%qm
!    k3%qmdot = p%EP_DeltaT * xdot%qmdot
!    x_tmp%qm    = x%qm    + k3%qm
!    x_tmp%qmdot = x%qmdot + k3%qmdot
!    
!    ! interpolate u to find u_interp = u(t + dt)
!    call LD_Input_ExtrapInterp(u, utimes, u_interp, t + p%EP_DeltaT, errStat, errMsg)
!    
!    ! find xdot at t + dt
!    call LD_CalcContStateDeriv( t + p%EP_DeltaT, u_interp, p, x_tmp, xd, z, OtherState, m, xdot, errStat, errMsg )
!    
!    k4%qm    = p%EP_DeltaT * xdot%qm
!    k4%qmdot = p%EP_DeltaT * xdot%qmdot
!    x%qm    = x%qm    +  ( k1%qm    + 2. * k2%qm    + 2. * k3%qm    + k4%qm    ) / 6.      
!    x%qmdot = x%qmdot +  ( k1%qmdot + 2. * k2%qmdot + 2. * k3%qmdot + k4%qmdot ) / 6.      
!    call ExitThisRoutine()
! CONTAinS      
!    !...............................................................................................................................
!    subroutine ExitThisRoutine()
!       ! This subroutine destroys all the local variables
!       integer(IntKi)             :: errStat3    ! The error identifier (errStat)
!       character(1024)            :: errMsg3     ! The error message (errMsg)
!       call LD_DestroyContState( xdot,     errStat3, errMsg3 )
!       call LD_DestroyContState( k1,       errStat3, errMsg3 )
!       call LD_DestroyContState( k2,       errStat3, errMsg3 )
!       call LD_DestroyContState( k3,       errStat3, errMsg3 )
!       call LD_DestroyContState( k4,       errStat3, errMsg3 )
!       call LD_DestroyContState( x_tmp,    errStat3, errMsg3 )
!       call LD_DestroyInput(     u_interp, errStat3, errMsg3 )
!    end subroutine ExitThisRoutine            
!       
! end subroutine LD_RK4
! 
! 
! !----------------------------------------------------------------------------------------------------------------------------------
! !> This is a loose coupling routine for solving constraint states, integrating continuous states, and updating discrete and other
! !! states. Continuous, constraint, discrete, and other states are updated to values at t + Interval.
! subroutine LD_UpdateStates( t, n, Inputs, InputTimes, p, x, xd, z, OtherState, m, errStat, errMsg )
! !..................................................................................................................................
!    real(DbKi),                         intent(in   ) :: t               !< Current simulation time in seconds
!    integer(IntKi),                     intent(in   ) :: n               !< Current step of the simulation: t = n*Interval
!    type(LD_InputType),            intent(inout) :: Inputs(:)       !< Inputs at InputTimes (output from this routine only
!                                                                         !!  because of record keeping in routines that copy meshes)
!    real(DbKi),                         intent(in   ) :: InputTimes(:)   !< Times in seconds associated with Inputs
!    type(LD_ParameterType),        intent(in   ) :: p               !< Parameters
!    type(LD_ContinuousStateType),  intent(inout) :: x               !< Input: Continuous states at t;
!                                                                         !!   Output: Continuous states at t + Interval
!    type(LD_DiscreteStateType),    intent(inout) :: xd              !< Input: Discrete states at t;
!                                                                         !!   Output: Discrete states at t + Interval
!    type(LD_ConstraintStateType),  intent(inout) :: z               !< Input: Constraint states at t;
!                                                                         !!   Output: Constraint states at t + Interval
!    type(LD_OtherStateType),       intent(inout) :: OtherState      !< Other states: Other states at t;
!                                                                         !!   Output: Other states at t + Interval
!    type(LD_MiscVarType),          intent(inout) :: m               !<  Misc variables for optimization (not copied in glue code)
!    integer(IntKi),                     intent(  out) :: errStat         !< Error status of the operation
!    character(*),                       intent(  out) :: errMsg          !< Error message if errStat /= ErrID_None
!    ! Initialize variables
!    errStat   = ErrID_None           ! no error has occurred
!    errMsg    = ""
!    if ( p%nCB == 0) return ! no modes = no states
!    if (p%IntMethod .eq. 1) then 
!       call LD_RK4( t, n, Inputs, InputTimes, p, x, xd, z, OtherState, m, errStat, errMsg )
!    elseif (p%IntMethod .eq. 2) then
!       call LD_AB4( t, n, Inputs, InputTimes, p, x, xd, z, OtherState, m, errStat, errMsg )
!    elseif (p%IntMethod .eq. 3) then
!       call LD_ABM4( t, n, Inputs, InputTimes, p, x, xd, z, OtherState, m, errStat, errMsg )
!    else  
!       call SeterrStat(ErrID_Fatal,'Invalid time integration method:'//Num2LStr(p%IntMethod),errStat,errMsg,'LD_UpdateState') 
!    end IF
! end subroutine LD_UpdateStates
! !----------------------------------------------------------------------------------------------------------------------------------
! !> This is a routine for computing outputs, used in both loose and tight coupling.
! subroutine LD_CalcOutput( t, u, p, x, xd, z, OtherState, y, m, errStat, errMsg )
!    real(DbKi),                        intent(in   )  :: t           !< Current simulation time in seconds
!    type(LD_InputType),           intent(in   )  :: u           !< Inputs at t
!    type(LD_ParameterType),       intent(in   )  :: p           !< Parameters
!    type(LD_ContinuousStateType), intent(in   )  :: x           !< Continuous states at t
!    type(LD_DiscreteStateType),   intent(in   )  :: xd          !< Discrete states at t
!    type(LD_ConstraintStateType), intent(in   )  :: z           !< Constraint states at t
!    type(LD_OtherStateType),      intent(in   )  :: OtherState  !< Other states at t
!    type(LD_MiscVarType),         intent(inout)  :: m           !< Misc variables for optimization (not copied in glue code)
!    type(LD_OutputType),          intent(inout)  :: y           !< Outputs computed at t (Input only so that mesh con-
!                                                                      !!   nectivity information does not have to be recalculated)
!    integer(IntKi),                    intent(  out)  :: errStat     !< Error status of the operation
!    character(*),                      intent(  out)  :: errMsg      !< Error message if errStat /= ErrID_None
!    ! Local variables
!    integer(IntKi)                                  :: I                 !< Generic counters
!    real(ReKi), dimension(6)                        :: Fc                !< Output coupling force
!    ! Compute the loads `fr1 fr2` at t (fr1 without added mass) by time interpolation of the inputs loads p%Forces
!    call InterpStpMat(real(t,ReKi), p%times, p%Forces, m%Indx, p%nTimeSteps, m%F_at_t)
! 
!    ! --- Flatening vectors and using linear state formulation y=Cx+Du+Fy
!    ! u flat (x1, \dot{x1}, \ddot{x1})
!    m%uFlat(1:3)   = u%PtfmMesh%TranslationDisp(:,1)
!    m%uFlat(4:6)   = GetSmllRotAngs(u%PtfmMesh%Orientation(:,:,1), errStat, errMsg); call SeterrStatSimple(errStat, errMsg, 'LD_CalcOutput')
!    m%uFlat(7:9  ) = u%PtfmMesh%TranslationVel(:,1)
!    m%uFlat(10:12) = u%PtfmMesh%RotationVel   (:,1)
!    m%uFlat(13:15) = u%PtfmMesh%TranslationAcc(:,1)
!    m%uFlat(16:18) = u%PtfmMesh%RotationAcc   (:,1)
! 
!    !--- Computing output:  y = Cx + Du + Fy  
!    ! 
!    if (p%nCB>0) then
!        ! x flat
!        m%xFlat(      1:p%nCB  ) = x%qm   (1:p%nCB)
!        m%xFlat(p%nCB+1:2*p%nCB) = x%qmdot(1:p%nCB)
! 
!        ! >>> MATMUL implementation
!        !Fc = matmul(p%CMat, m%xFlat) + matmul(p%DMat, m%uFlat) + m%F_at_t(1:6) - matmul(p%M12, m%F_at_t(6+1:6+p%nCB))
! 
!        ! >>> LAPACK implementation
!        Fc(1:6) = m%F_at_t(1:6) ! Fc = F1r + ...
!        !           GEMV(TRS, M  , N      , alpha    , A     , LDA, X                    ,inCX, Beta  ,  Y, IncY)
!        call LAPACK_GEMV('n', 6  , 2*p%nCB,  1.0_ReKi, p%CMat, 6  , m%xFlat              , 1, 1.0_ReKi, Fc, 1   ) ! = C*x + (F1r)
!        call LAPACK_GEMV('n', 6  ,   18   ,  1.0_ReKi, p%DMat, 6  , m%uFlat              , 1, 1.0_ReKi, Fc, 1   ) ! + D*u
!        call LAPACK_GEMV('n', 6  , p%nCB  , -1.0_ReKi, p%M12 , 6  , m%F_at_t(6+1:6+p%nCB), 1, 1.0_ReKi, Fc, 1   ) ! - M12*F2r
!    else
!        Fc =                           matmul(p%DMat, m%uFlat) + m%F_at_t(1:6) 
!    endif
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
! end subroutine LD_CalcOutput
! !----------------------------------------------------------------------------------------------------------------------------------
! 
! 
! !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! !> This is a tight coupling routine for computing derivatives of continuous states.
! subroutine LD_CalcContStateDeriv( t, u, p, x, xd, z, OtherState, m, dxdt, errStat, errMsg )
! !..................................................................................................................................
!    real(DbKi),                        intent(in   )  :: t           !< Current simulation time in seconds
!    type(LD_InputType),           intent(in   )  :: u           !< Inputs at t
!    type(LD_ParameterType),       intent(in   )  :: p           !< Parameters
!    type(LD_ContinuousStateType), intent(in   )  :: x           !< Continuous states at t
!    type(LD_DiscreteStateType),   intent(in   )  :: xd          !< Discrete states at t
!    type(LD_ConstraintStateType), intent(in   )  :: z           !< Constraint states at t
!    type(LD_OtherStateType),      intent(in   )  :: OtherState  !< Other states at t
!    type(LD_MiscVarType),         intent(inout)  :: m           !< Misc variables for optimization (not copied in glue code)
!    type(LD_ContinuousStateType), intent(  out)  :: dxdt        !< Continuous state derivatives at t
!    integer(IntKi),                    intent(  out)  :: errStat     !< Error status of the operation
!    character(*),                      intent(  out)  :: errMsg      !< Error message if errStat /= ErrID_None
!    ! Local variables
!    integer(IntKi)                                    :: I
!    ! Allocation of output dxdt (since intent(out))
!    call AllocAry(dxdt%qm,    p%nCB, 'dxdt%qm',    errStat, errMsg); if(Failed()) return
!    call AllocAry(dxdt%qmdot, p%nCB, 'dxdt%qmdot', errStat, errMsg); if(Failed()) return
!    if ( p%nCB == 0 ) return
!    do I=1,p%nCB; dxdt%qm   (I)=0; enddo
!    do I=1,p%nCB; dxdt%qmdot(I)=0; enddo
! 
!    ! Compute the loads `fr1 fr2` at t (fr1 without added mass) by time interpolation of the inputs loads p%F
!    call InterpStpMat(real(t,ReKi), p%times, p%Forces, m%Indx, p%nTimeSteps, m%F_at_t)
! 
!    ! u flat (x1, \dot{x1}, \ddot{x1})
!    m%uFlat(1:3)   = u%PtfmMesh%TranslationDisp(:,1)
!    m%uFlat(4:6)   = GetSmllRotAngs(u%PtfmMesh%Orientation(:,:,1), errStat, errMsg); if(Failed()) return
!    m%uFlat(7:9  ) = u%PtfmMesh%TranslationVel(:,1)
!    m%uFlat(10:12) = u%PtfmMesh%RotationVel   (:,1)
!    m%uFlat(13:15) = u%PtfmMesh%TranslationAcc(:,1)
!    m%uFlat(16:18) = u%PtfmMesh%RotationAcc   (:,1)
! 
!    ! --- Computation of qm and qmdot
!    ! >>> Latex formulae:
!    ! \ddot{x2} = -K22 x2 - C22 \dot{x2}  - C21 \dot{x1} - M21 \ddot{x1} + fr2
!    ! >>> MATMUL IMPLEMENTATION 
!    !dxdt%qm= x%qmdot
!    !dxdt%qmdot = - matmul(p%K22,x%qm) - matmul(p%C22,x%qmdot) &
!    !             - matmul(p%C21,m%uFlat(7:12)) - matmul(p%M21, m%uFlat(13:18)) + m%F_at_t(6+1:6+p%nCB)
!    ! >>> BLAS IMPLEMENTATION 
!    !           COPY( N   , X                    , inCX, Y      , inCY)
!    call LAPACK_COPY(p%nCB, x%qmdot              , 1  , dxdt%qm    , 1  ) ! qmdot=qmdot
!    call LAPACK_COPY(p%nCB, m%F_at_t(6+1:6+p%nCB), 1  , dxdt%qmdot , 1  )                                          ! qmddot = fr2
!    !           GEMV(TRS, M    , N     , alpha    , A    , LDA  , X              ,inCX, Beta   ,  Y        , IncY)
!    call LAPACK_GEMV('n', p%nCB, p%nCB , -1.0_ReKi, p%K22, p%nCB, x%qm          , 1  , 1.0_ReKi, dxdt%qmdot, 1   ) !        - K22 x2
!    call LAPACK_GEMV('n', p%nCB, 6     , -1.0_ReKi, p%C21, p%nCB, m%uFlat(7:12) , 1  , 1.0_ReKi, dxdt%qmdot, 1   ) !        - C21 \dot{x1}
!    call LAPACK_GEMV('n', p%nCB, p%nCB , -1.0_ReKi, p%C22, p%nCB, x%qmdot       , 1  , 1.0_ReKi, dxdt%qmdot, 1   ) !        - C22 \dot{x2}
!    call LAPACK_GEMV('n', p%nCB, 6     , -1.0_ReKi, p%M21, p%nCB, m%uFlat(13:18), 1  , 1.0_ReKi, dxdt%qmdot, 1   ) !        - M21 \ddot{x1}
! 
! CONTAinS
!     logical function Failed()
!         call SeterrStatSimple(errStat, errMsg, 'LD_CalcContStateDeriv')
!         Failed =  errStat >= AbortErrLev
!     end function Failed
! end subroutine LD_CalcContStateDeriv
! !----------------------------------------------------------------------------------------------------------------------------------
! !> This is a tight coupling routine for updating discrete states.
! subroutine LD_UpdateDiscState( t, n, u, p, x, xd, z, OtherState, m, errStat, errMsg )
! !..................................................................................................................................
!    real(DbKi),                        intent(in   )  :: t           !< Current simulation time in seconds
!    integer(IntKi),                    intent(in   )  :: n           !< Current step of the simulation: t = n*Interval
!    type(LD_InputType),           intent(in   )  :: u           !< Inputs at t
!    type(LD_ParameterType),       intent(in   )  :: p           !< Parameters
!    type(LD_ContinuousStateType), intent(in   )  :: x           !< Continuous states at t
!    type(LD_DiscreteStateType),   intent(inout)  :: xd          !< Input: Discrete states at t, Output: Discrete states at t + Interval
!    type(LD_ConstraintStateType), intent(in   )  :: z           !< Constraint states at t
!    type(LD_OtherStateType),      intent(in   )  :: OtherState  !< Other states at t
!    type(LD_MiscVarType),         intent(inout)  :: m           !< Misc variables for optimization (not copied in glue code)
!    integer(IntKi),                    intent(  out)  :: errStat     !< Error status of the operation
!    character(*),                      intent(  out)  :: errMsg      !< Error message if errStat /= ErrID_None
!    ! Initialize errStat
!    errStat = ErrID_None
!    errMsg  = ""
!    ! Update discrete states here:
!    xd%DummyDiscState = 0.0_Reki
! end subroutine LD_UpdateDiscState
! !----------------------------------------------------------------------------------------------------------------------------------
! !> This is a tight coupling routine for solving for the residual of the constraint state functions.
! subroutine LD_CalcConstrStateResidual( t, u, p, x, xd, z, OtherState, m, Z_residual, errStat, errMsg )
! !..................................................................................................................................
!    real(DbKi),                        intent(in   )  :: t           !< Current simulation time in seconds
!    type(LD_InputType),           intent(in   )  :: u           !< Inputs at t
!    type(LD_ParameterType),       intent(in   )  :: p           !< Parameters
!    type(LD_ContinuousStateType), intent(in   )  :: x           !< Continuous states at t
!    type(LD_DiscreteStateType),   intent(in   )  :: xd          !< Discrete states at t
!    type(LD_ConstraintStateType), intent(in   )  :: z           !< Constraint states at t (possibly a guess)
!    type(LD_OtherStateType),      intent(in   )  :: OtherState  !< Other states at t
!    type(LD_MiscVarType),         intent(inout)  :: m           !< Misc variables for optimization (not copied in glue code)
!    type(LD_ConstraintStateType), intent(  out)  :: Z_residual  !< Residual of the constraint state functions using
!                                                                     !!     the input values described above
!    integer(IntKi),                    intent(  out)  :: errStat     !< Error status of the operation
!    character(*),                      intent(  out)  :: errMsg      !< Error message if errStat /= ErrID_None
!    ! Initialize errStat
!    errStat = ErrID_None
!    errMsg  = ""
!    ! Solve for the residual of the constraint state functions here:
!    Z_residual%DummyConstrState = 0.0_ReKi
! 
! end subroutine LD_CalcConstrStateResidual
! !----------------------------------------------------------------------------------------------------------------------------------
! 
! !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ! ###### The following four routines are Jacobian routines for linearization capabilities #######
! ! If the module does not implement them, set errStat = ErrID_Fatal in LD_Init() when InitInp%Linearize is .true.
! !----------------------------------------------------------------------------------------------------------------------------------
! !> Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) functions
! !! with respect to the inputs (u). The partial derivatives dY/du, dX/du, dXd/du, and DZ/du are returned.
! 
! subroutine LD_JacobianPInput( t, u, p, x, xd, z, OtherState, y, m, errStat, errMsg, dYdu, dXdu, dXddu, dZdu)
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
