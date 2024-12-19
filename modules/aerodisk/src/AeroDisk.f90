!**********************************************************************************************************************************
!> ## AeroDisk
!! The AeroDisk module solves a quasi-steady actuator disk representation of the rotor to calculate the 3 forces and 3 moments of
!! the rotor dependent on the tip-speed ratio (TSR), rotor speed (RotSpeed), relative wind velocity vector (VRel), and the rotor-
!! collective blade-pitch (BlPitch).
!!
! ..................................................................................................................................
!! ## LICENSING
!! Copyright (C) 2024  National Renewable Energy Laboratory
!!
!!    This file is part of AeroDisk.
!!
!! Licensed under the Apache License, Version 2.0 (the "License");
!! you may not use this file except in compliance with the License.
!! You may obtain a copy of the License at
!!
!!     http://www.apache.org/licenses/LICENSE-2.0
!!
!! Unless required by applicable law or agreed to in writing, software
!! distributed under the License is distributed on an "AS IS" BASIS,
!! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!! See the License for the specific language governing permissions and
!! limitations under the License.
!**********************************************************************************************************************************
MODULE AeroDisk

   USE AeroDisk_Types
   USE AeroDisk_IO
   USE NWTC_Library
   use IfW_FlowField, only: IfW_FlowField_GetVelAcc


   implicit none
   private
   type(ProgDesc), parameter :: ADsk_Ver = ProgDesc( 'AeroDisk', '', '' )

   public :: ADsk_Init
   public :: ADsk_End
   public :: ADsk_UpdateStates
   public :: ADsk_CalcOutput
   public :: ADsk_CalcContStateDeriv

   ! Linearization is not supported by this module, so the following routines are omitted
   !public :: ADsk_CalcConstrStateResidual
   !public :: ADsk_UpdateDiscState
   !public :: ADsk_JacobianPInput
   !public :: ADsk_JacobianPContState
   !public :: ADsk_JacobianPDiscState
   !public :: ADsk_JacobianPConstrState
   !public :: ADsk_GetOP

CONTAINS


!----------------------------------------------------------------------------------------------------------------------------------
!> Initialize the AeroDisk module:
!!    - load settings (passed or from file)
!!    - setup meshes
!!    - initialize outputs and other data storage
SUBROUTINE ADsk_Init( InitInp, u, p, x, xd, z, OtherState, y, m, Interval, InitOut, ErrStat, ErrMsg )
   type(ADsk_InitInputType),        intent(in   )  :: InitInp     !< Input data for initialization routine
   type(ADsk_InputType),            intent(  out)  :: u           !< An initial guess for the input; input mesh must be defined
   type(ADsk_ParameterType),        intent(  out)  :: p           !< Parameters
   type(ADsk_ContinuousStateType),  intent(  out)  :: x           !< Initial continuous states
   type(ADsk_DiscreteStateType),    intent(  out)  :: xd          !< Initial discrete states
   type(ADsk_ConstraintStateType),  intent(  out)  :: z           !< Initial guess of the constraint states
   type(ADsk_OtherStateType),       intent(  out)  :: OtherState  !< Initial other states (logical, etc)
   type(ADsk_OutputType),           intent(  out)  :: y           !< Initial system outputs (outputs are not calculated)
   type(ADsk_MiscVarType),          intent(  out)  :: m           !< Misc variables for optimization (not copied in glue code)
   real(DbKi),                      intent(inout)  :: Interval    !< Coupling interval in seconds
   type(ADsk_InitOutputType),       intent(  out)  :: InitOut     !< Output for initialization routine
   integer(IntKi),                  intent(  out)  :: ErrStat     !< Error status of the operation
   character(*),                    intent(  out)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

   ! local variables
   type(ADsk_InputFile)                            :: InputFileData  !< Data from input file as a string array
   type(FileInfoType)                              :: FileInfo_In !< The derived type for holding the full input file for parsing -- we may pass this in the future
   integer(IntKi)                                  :: UnEc        ! unit number for the echo file (-1 for not in use)
   integer(IntKi)                                  :: ErrStat2    ! local error status
   character(ErrMsgLen)                            :: ErrMsg2     ! local error message
   character(*), parameter                         :: RoutineName = 'ADsk_Init'

      ! Initialize variables
   ErrStat = ErrID_None
   ErrMsg  = ""

   ! Initialize the NWTC Subroutine Library
   call NWTC_Init( )

   ! Display the module information
   call DispNVD( ADsk_Ver )

   ! set rootname
   p%RootName = trim(InitInp%RootName)//".ADsk"

   ! Get primary input file
   if ( InitInp%UseInputFile ) then
      CALL ProcessComFile( InitInp%InputFile, FileInfo_In, ErrStat2, ErrMsg2 )
   else
      CALL NWTC_Library_CopyFileInfoType( InitInp%PassedFileData, FileInfo_In, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
   endif
   if (Failed()) return

   ! For diagnostic purposes, the following can be used to display the contents
   ! of the FileInfo_In data structure.
   !call Print_FileInfo_Struct( CU, FileInfo_In ) ! CU is the screen -- different number on different systems.

   ! Parse all ADsk-related input and populate the InputFileData structure
   call ADsk_ParsePrimaryFileData( InitInp, p%RootName, Interval, FileInfo_In, InputFileData, UnEc, ErrStat2, ErrMsg2 )
   if (Failed()) return;

   ! Verify all the necessary initialization and input file data
   CALL ADskInput_ValidateInput( InitInp, InputFileData, ErrStat2, ErrMsg2 )
   if (Failed()) return;

   ! Set parameters
   CALL ADskInput_SetParameters( InitInp, Interval, InputFileData, p, ErrStat2, ErrMsg2 )
   if (Failed()) return;

   ! For diagnostic purposes.  If we add a summary file, use this to write table
   !call WriteAeroTab(p%AeroTable,Cu)


   ! Set pointer to FlowField data
   if (associated(InitInp%FlowField)) then
      p%FlowField => InitInp%FlowField
      call SetDiskAvgPoints(ErrStat2,ErrMsg2);  if (Failed()) return
   else
      ErrStat2 = ErrID_Fatal
      ErrMsg2  = "No flow field data available.  AeroDisk cannot continue."
      if (Failed()) return
   endif


   ! Set inputs
   call Init_U(ErrStat2,ErrMsg2);   if (Failed())  return

   ! Set outputs
   call Init_Y(ErrStat2,ErrMsg2);   if (Failed())  return

   ! Set InitOutputs
   call Init_InitY(ErrStat2,ErrMsg2);  if (Failed())  return

   ! Set some other stuff that the framework requires
   call Init_OtherStuff(ErrStat2,ErrMsg2);  if (Failed())  return

contains
   logical function Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        Failed =  ErrStat >= AbortErrLev
        if (Failed) call CleanUp()
   end function Failed

   subroutine Cleanup()
      if (UnEc > 0_IntKi)  close (UnEc)
   end subroutine Cleanup

   !> Setup points for disk average velocity calculations
   subroutine SetDiskAvgPoints(ErrStat3,ErrMsg3)
      integer(IntKi),   intent(  out)  :: ErrStat3
      character(*),     intent(  out)  :: ErrMsg3
      integer(IntKi) :: i
      real(ReKi)     :: R,theta
      ! positions relative to hub
      call AllocAry(p%DiskWindPosRel,3,ADsk_NumPtsDiskAvg,'ADsk_NumPtsDiskAvg',ErrStat3,ErrMsg3);  if (errStat3 >= AbortErrLev) return
      ! absolute point positions for call to GetWindVelAcc
      call AllocAry(m%DiskWindPosAbs,3,ADsk_NumPtsDiskAvg,'ADsk_NumPtsDiskAvg',ErrStat3,ErrMsg3);  if (errStat3 >= AbortErrLev) return
      ! wind velocity at all requested points
      call AllocAry(m%DiskWindVel   ,3,ADsk_NumPtsDiskAvg,'DiskWindVel'       ,ErrStat3,ErrMsg3);  if (errStat3 >= AbortErrLev) return
      ! Calculate relative points on disk (do this once up front to save computational time).
      ! NOTE: this is in the XY plane, and will be multiplied by the hub orientation vector
      R = real(p%RotorRad,ReKi) * 0.7_reKi !70% radius
      do i=1,ADsk_NumPtsDiskAvg
         theta = pi +(i-1)*TwoPi/ADsk_NumPtsDiskAvg
         p%DiskWindPosRel(1,i) = 0.0_ReKi          ! Hub X (perpindicular to rotor plane)
         p%DiskWindPosRel(2,i) = R*cos(theta)      ! Hub Y
         p%DiskWindPosRel(3,i) = R*sin(theta)      ! Hub Z (in vertical plane when azimuth=0)
      end do
   end subroutine SetDiskAvgPoints

   !> Initialize the inputs in u
   subroutine Init_U(ErrStat3,ErrMsg3)
      integer(IntKi),   intent(  out)  :: ErrStat3
      character(*),     intent(  out)  :: ErrMsg3
      ! HubMotion mesh
      call MeshCreate ( BlankMesh  = u%HubMotion  &
                     ,IOS       = COMPONENT_INPUT &
                     ,Nnodes    = 1               &
                     ,ErrStat   = ErrStat3        &
                     ,ErrMess   = ErrMsg3         &
                     ,Orientation     = .true.    &
                     ,TranslationDisp = .true.    &
                     ,RotationVel     = .true.    &
                     )
         if (errStat3 >= AbortErrLev) return
      call MeshPositionNode(u%HubMotion, 1, InitInp%HubPosition, errStat3, errMsg3, InitInp%HubOrientation);   if (errStat3 >= AbortErrLev) return
      call MeshConstructElement( u%HubMotion, ELEMENT_POINT, errStat3, errMsg3, p1=1 );   if (errStat3 >= AbortErrLev) return
      call MeshCommit(u%HubMotion, errStat3, errMsg3 );   if (errStat3 >= AbortErrLev) return
      u%HubMotion%Orientation     = u%HubMotion%RefOrientation
      u%HubMotion%TranslationDisp = 0.0_R8Ki
      u%HubMotion%RotationVel     = 0.0_ReKi
      return
   end subroutine Init_U

   !> Initialize the outputs in Y
   subroutine Init_Y(ErrStat3,ErrMSg3)
      integer(IntKi),   intent(  out)  :: ErrStat3
      character(*),     intent(  out)  :: ErrMsg3
      ! Set output loads mesh
      call MeshCopy (  SrcMesh  = u%HubMotion        &
                     , DestMesh = y%AeroLoads        &
                     , CtrlCode = MESH_SIBLING       &
                     , IOS      = COMPONENT_OUTPUT   &
                     , force    = .TRUE.             &
                     , moment   = .TRUE.             &
                     , ErrStat  = ErrStat3           &
                     , ErrMess  = ErrMsg3            )
         if (ErrStat3 >= AbortErrLev) return

      ! Initialize all outputs to zero (will be set by CalcOutput)
      y%YawErr    = 0.0_ReKi
      y%PsiSkew   = 0.0_ReKi
      y%ChiSkew   = 0.0_ReKi
      y%VRel      = 0.0_ReKi
      y%Ct        = 0.0_ReKi
      y%Cq        = 0.0_ReKi
      call AllocAry(y%WriteOutput,p%NumOuts,'WriteOutput',Errstat3,ErrMsg3);  if (ErrStat3 >= AbortErrLev) return
      y%WriteOutput = 0.0_ReKi
   end subroutine Init_Y

   !> Initialize other stuff that the framework requires, but isn't used here
   subroutine Init_OtherStuff(ErrStat3,ErRMsg3)
      integer(IntKi),   intent(  out)  :: ErrStat3
      character(*),     intent(  out)  :: ErrMsg3
      ErrStat3 = ErrID_None
      ErrMsg3  = ""
      x%DummyContState           = 0.0_ReKi
      xd%DummyDiscreteState      = 0.0_ReKi
      z%DummyConstrState         = 0.0_ReKi
      OtherState%DummyOtherState = 0_IntKi
      m%idx_last                 = 1_IntKi      ! Aerotable lookup indice
      if (allocated(m%AllOuts)) deallocate(m%AllOuts)
      allocate(m%AllOuts(0:MaxOutPts),STAT=ErrStat3)
      if (ErrStat3 /= 0) then
         ErrStat3 = ErrID_Fatal
         ErrMsg3  = "Cannot allocate m%AllOuts"
         return
      endif
      m%AllOuts = 0.0_SiKi
   end subroutine Init_OtherStuff

   !> Initialize the InitOutput
   subroutine Init_InitY(ErrStat3,ErrMsg3)
      integer(IntKi),   intent(  out)  :: ErrStat3
      character(*),     intent(  out)  :: ErrMsg3
      integer(IntKi)                   :: i
      call AllocAry(InitOut%WriteOutputHdr,p%NumOuts,'WriteOutputHdr',Errstat3,ErrMsg3);  if (ErrStat3 >= AbortErrLev) return
      call AllocAry(InitOut%WriteOutputUnt,p%NumOuts,'WriteOutputUnt',Errstat3,ErrMsg3);  if (ErrStat3 >= AbortErrLev) return
      do i=1,p%NumOuts
         InitOut%WriteOutputHdr(i) = p%OutParam(i)%Name
         InitOut%WriteOutputUnt(i) = p%OutParam(i)%Units
      end do
      ! Version
      InitOut%Ver     = ADsk_Ver
      InitOut%AirDens = p%AirDens
   end subroutine Init_InitY
END SUBROUTINE ADsk_Init


!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is called at the end of the simulation.
SUBROUTINE ADsk_End( u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )
   type(ADsk_InputType),            intent(inout)  :: u           !< System inputs
   type(ADsk_ParameterType),        intent(inout)  :: p           !< Parameters
   type(ADsk_ContinuousStateType),  intent(inout)  :: x           !< Continuous states
   type(ADsk_DiscreteStateType),    intent(inout)  :: xd          !< Discrete states
   type(ADsk_ConstraintStateType),  intent(inout)  :: z           !< Constraint states
   type(ADsk_OtherStateType),       intent(inout)  :: OtherState  !< Other states
   type(ADsk_OutputType),           intent(inout)  :: y           !< System outputs
   type(ADsk_MiscVarType),          intent(inout)  :: m           !< Misc variables for optimization (not copied in glue code)
   integer(IntKi),                  intent(  out)  :: ErrStat     !< Error status of the operation
   character(*),                    intent(  out)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

   ! local variables
   integer(IntKi)                                  :: ErrStat2    ! local error status
   character(ErrMsgLen)                            :: ErrMsg2     ! local error message
   character(*), parameter                         :: RoutineName = 'ADsk_End'

      ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ""

      !! Place any last minute operations or calculations here:

      !! Close files here (but because of checkpoint-restart capability, it is not recommended to have files open during the simulation):

   ! Destroy the input data:
   call ADsk_DestroyInput( u, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

   ! Destroy the parameter data:
   call ADsk_DestroyParam( p, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

   ! Destroy the state data:
   call ADsk_DestroyContState(   x,          ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call ADsk_DestroyDiscState(   xd,         ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call ADsk_DestroyConstrState( z,          ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call ADsk_DestroyOtherState(  OtherState, ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

   ! Destroy the output data:
   call ADsk_DestroyOutput( y, ErrStat2, ErrMsg2 ); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

   ! Destroy the misc data:
   call ADsk_DestroyMisc( m, ErrStat2, ErrMsg2 ); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
END SUBROUTINE ADsk_End


!----------------------------------------------------------------------------------------------------------------------------------
!> This is a routine for computing outputs, used in both loose and tight coupling.
!! It calculates the forces and moments on the rotor disk given an orientation, motion, blade pitch, rotor speed, and wind velocity.
!! The calculations are based on an interpolation into input table values.  Since the table is stored as single kind and this method
!! makes several simplifying assumptions that introduuce error, all calculations here are performed in single precision.
SUBROUTINE ADsk_CalcOutput( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, NeedWriteOutput )
   real(DbKi),                      intent(in   )  :: t           !< Current simulation time in seconds
   type(ADsk_InputType),            intent(in   )  :: u           !< Inputs at t
   type(ADsk_ParameterType),        intent(in   )  :: p           !< Parameters
   type(ADsk_ContinuousStateType),  intent(in   )  :: x           !< Continuous states at t
   type(ADsk_DiscreteStateType),    intent(in   )  :: xd          !< Discrete states at t
   type(ADsk_ConstraintStateType),  intent(in   )  :: z           !< Constraint states at t
   type(ADsk_OtherStateType),       intent(in   )  :: OtherState  !< Other states at t
   type(ADsk_MiscVarType),          intent(inout)  :: m           !< Misc variables for optimization (not copied in glue code)
   type(ADsk_OutputType),           intent(inout)  :: y           !< Outputs computed at t (Input only for mesh)
   integer(IntKi),                  intent(  out)  :: ErrStat     !< Error status of the operation
   character(*),                    intent(  out)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None
   logical,             optional,   intent(in   )  :: NeedWriteOutput   !< Flag to determine if WriteOutput values need to be calculated in this call

   ! local variables
   real(ReKi), allocatable                         :: NoAcc(:,:)     ! Placeholder array not used when accelerations not required.
   real(SiKi)                                      :: x_hatDisk(3)   ! X direction unit vector of rotor disk (global)
   real(SiKi)                                      :: y_hatDisk(3)   ! Y direction unit vector of rotor disk (global)
   real(SiKi)                                      :: z_hatDisk(3)   ! Z direction unit vector of rotor disk (global)
   real(SiKi)                                      :: VRel_vec(3)    ! relative velocity of wind in moving rotor disk frame
   real(SiKi)                                      :: tmp1,tmp3(3)   ! temporary variables for calculations
   integer(IntKi)                                  :: i              ! generic counter
   integer(IntKi)                                  :: ErrStat2       ! local error status
   character(ErrMsgLen)                            :: ErrMsg2        ! local error message
   character(*), parameter                         :: RoutineName = 'ADsk_CalcOutput'
   logical                                         :: CalcWriteOutput

   ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ""
   m%AllOuts = 0.0_SiKi
   
   if (present(NeedWriteOutput)) then
      CalcWriteOutput = NeedWriteOutput
   else
      CalcWriteOutput = .true. ! by default, calculate WriteOutput unless told that we do not need it
   end if


   !--------------------------------
   !> Disk average wind speed
   call CalcDiskAvgVel(ErrStat2,ErrMsg2);  if (FAiled()) return

   !--------------------------------
   !> Disk vectors (in global frame)
   !!
   !! | Vector                      | Description                                                                                 | Name      |
   !! | :-------------------------- | :------------------------------------------------------------------------------------------ | :-------- |
   !! | \f$\hat{x}_\textrm{disk}\f$ | normal to disk                                                                              | x_hatDisk |
   !! | \f$\hat{y}_\textrm{disk}\f$ | perpendicular to \f$\hat{x}_\textrm{disk}\f$ in the \f$x-y\f$ plane                         | y_hatDisk |
   !! | \f$\hat{z}_\textrm{disk}\f$ | right handed coordinate axis to \f$\hat{x}_\textrm{disk}\f$ and \f$\hat{y}_\textrm{disk}\f$ | z_hatDisk |
   ! Normal to disk
   x_hatDisk(1:3) = real(u%HubMotion%Orientation(1:3,1,1), SiKi)
   ! unit vector of disk normal projected to XY plane (calculate from cross_product(Z_global, x_hatDisk) )
   y_hatDisk(1) = real(-x_hatDisk(2),SiKi)
   y_hatDisk(2) = real(x_hatDisk(1),SiKi)
   y_hatDisk(3) = 0.0_SiKi
   y_hatDisk = y_hatDisk / TwoNorm( y_hatDisk )     ! normalize
   ! unit vector normal to x_hatDisk and y_hatDisk -- cross_product(x_hatDisk, y_hatDisk)
   z_hatDisk(1) = real(-x_hatDisk(1)*x_hatDisk(3),SiKi)
   z_hatDisk(2) = real(-x_hatDisk(2)*x_hatDisk(3),SiKi)
   z_hatDisk(3) = real(x_hatDisk(1)*x_hatDisk(1) + x_hatDisk(2)*x_hatDisk(2),SiKi)


   !-------------
   ! Error checks
   !-------------
   ! Verify rotor and wind orientation are not pointed vertically
   if ((EqualRealNos(m%DiskAvgVel(1), 0.0_ReKi) .and. EqualRealNos(m%DiskAvgVel(2), 0.0_ReKi)) .and. (.not. EqualRealNos(m%DiskAvgVel(3), 0.0_ReKi))) then
      ErrStat2 = ErrID_Fatal;    ErrMsg2 = "AeroDisk cannot calculate aero loads with wind in the vertical direction. Nacelle yaw-error undefined."
      if (Failed())  return
   endif
   if (EqualRealNos(x_hatDisk(1), 0.0_SiKi) .and. EqualRealNos(x_hatDisk(2), 0.0_SiKi)) then
      ErrStat2 = ErrID_Fatal;    ErrMsg2 = "AeroDisk cannot calculate aero loads with rotor pointed in vertical direction. Nacelle yaw-error undefined."
      if (Failed())  return
   endif

   !-------------------------
   !> Calculate some constants
   !!  - \f$\vec{V}_\textrm{rel} = \vec{V}_\textrm{wind} - \vec{v}_\textrm{rotor}\f$
   !!  - \f$V_\textrm{rel} = \left\| \vec{V}_\textrm{rel} \right\|_2\f$
   !!  - \f$V_\textrm{rel,x-disk} = \vec{V}_\textrm{rel_g} \bullet \hat{x}_\textrm{disk}\f$
   !-------------------------
   ! Calculate relative wind velocity (global)
   VRel_vec(1:3) = real(m%DiskAvgVel(1:3) - u%HubMotion%TranslationVel(1:3,1), SiKi)
   ! Magnitude of relative wind velocity
   m%VRel = TwoNorm(VRel_vec)
   ! relative wind velocity along disk normal
   m%VRel_xd = abs(dot_product(VRel_vec, x_hatDisk))
   ! set output
   y%VRel = m%VRel

   !> calculate Lambda (TSR) and ChiSkew (inflow skew angle)
   !!  - \f$\lambda = \left\{  \begin{matrix}
   !!              0  &\textrm{for}& \vec{V}_\textrm{rel,x-disk} = 0 \\
   !!              \frac{\Omega R}{\vec{V}_\textrm{rel,x-disk}} &\textrm{otherwise}&
   !!              \end{matrix} \right.\f$
   !!  - \f$ \chi = \left\{  \begin{matrix}
   !!              0  &\textrm{for}& V_\textrm{rel} = 0 \\
   !!              \text{ACOS}\left(\frac{\vec{V}_\textrm{rel}\bullet\hat{x}_\textrm{disk}}{V_\textrm{rel}}\right) &\textrm{otherwise}&
   !!              \end{matrix} \right.\f$
   !!
   if (EqualRealNos(m%VRel_xd,0.0_SiKi)) then
      m%lambda = 0.0_SiKi
   else
      m%lambda = real((u%RotSpeed * p%RotorRad),SiKi) / abs(m%VRel_xd)
   endif
   if (EqualRealNos(m%VRel,0.0_SiKi)) then
      m%Chi = 0.0_SiKi
   else
      m%Chi = acos( m%VRel_xd / m%VRel )
   endif
   y%ChiSkew = m%Chi

   !--------------
   !> x,y,z vectors -- convert global coordinates to disk coordinates
   !! - \f$ A_\textrm{tmp} = \left\| \left( \vec{V}_\textrm{rel} \bullet \hat{x}_\textrm{disk} \right) \hat{x}_\textrm{disk} - \vec{V}_\textrm{rel} \right\|_2 \f$
   !! - \f$ \hat{x} = \hat{x}_\textrm{disk}\f$
   !! - \f$ \hat{y} = \left\{ \begin{matrix}
   !!                \hat{y}_\textrm{disk} = \hat{y}_\textrm{disk} &\textrm{for}& A_\textrm{tmp} = 0\\
   !!                \hat{y}_\textrm{disk} = \frac{\left(\vec{V}_\textrm{rel} \bullet \hat{x}_\textrm{disk}\right) \hat{x}_\textrm{disk} - \vec{V}_\textrm{rel}}{A_\textrm{tmp}} &\textrm{otherwise}&
   !!             \end{matrix}\right. \f$
   !! - \f$ \hat{z} = \left\{ \begin{matrix}
   !!                \hat{z}_\textrm{disk} = \hat{z}_\textrm{disk} &\textrm{for}& A_\textrm{tmp} = 0\\
   !!                \hat{z}_\textrm{disk} = \frac{\vec{V}_\textrm{rel} \times \hat{x}_\textrm{disk}}{A_\textrm{tmp}} &\textrm{otherwise}&
   !!             \end{matrix}\right.\f$
   tmp3 = dot_product(VRel_vec, x_hatDisk) * x_hatDisk - VRel_vec
   tmp1 = TwoNorm(tmp3)
   if (EqualRealNos(tmp1, 0.0_SiKi))  then
      m%x_hat = x_hatDisk
      m%y_hat = y_hatDisk
      m%z_hat = z_hatDisk
   else
      m%x_hat = x_hatDisk
      m%y_hat = tmp3 / tmp1
      m%z_hat = cross_product( VRel_vec, x_hatDisk ) / tmp1
   endif

   !---------------
   !> YawErr and Skew
   !! - YawErr:       \f$ \gamma_\textrm{YawErr} = \text{ATAN2}\left(
   !!          \vec{V}_\textrm{rel}[2] \hat{x}_\textrm{disk}[1] - \vec{V}_\textrm{rel}[1] \hat{x}_\textrm{disk}[2],
   !!          \vec{V}_\textrm{rel}[1] \hat{x}_\textrm{disk}[1] + \vec{V}_\textrm{rel}[2] \hat{x}_\textrm{disk}[2] \right) \f$
   !! - PsiSkew:    \f$ \Psi_\textrm{skew} = \text{ATAN2}\left( \hat{z} \bullet \hat{y}_\textrm{disk}, -\hat{z}\bullet \hat{z}_\textrm{disk} \right) \f$
   y%YawErr  = atan2( VRel_vec(2)*x_hatDisk(1) - VRel_vec(1)*x_hatDisk(2),   VRel_vec(1)*x_hatDisk(1) + VRel_vec(2)*x_hatDisk(2) )
   y%PsiSkew = atan2( -1.0_ReKi * real( dot_product(m%z_hat,y_hatDisk), ReKi), real( dot_product(m%z_hat, z_hatDisk), ReKi) )
  

   !-------------------------------------------
   !> Interpolate Force and Moment coefficients
   call ADskTableInterp(p%AeroTable, p%UseTSR, m%lambda, real(u%RotSpeed,SiKi), m%VRel_xd, real(u%BlPitch,SiKi), m%Chi, m%idx_last, m%C_F, m%C_M, ErrStat2, ErrMsg2)
      if (Failed()) return

   !> Apply skew if not in table
   !! - \f$ \vec{F} = \vec{F} \left( \cos(\chi) \right)^2 \f$
   !! - \f$ \vec{M} = \vec{M} \left( \cos(\chi) \right)^2 \f$
   if (p%AeroTable%N_Skew  <= 0_IntKi) then
      tmp1  = cos(m%Chi) * cos(m%Chi)
      m%C_F(1:3)  = m%C_F(1:3) * tmp1
      m%C_M(1:3)  = m%C_M(1:3) * tmp1
   endif

   !-------------------------------------------
   !> Calculate forces using force coefficients (disk coordinates)
   !! - \f$ F_x = \frac{1}{2} \rho A \left( V_\textrm{rel,x} \right)^2 * C_\textrm{F,x}\left(\text{TSR}@\lambda,\text{RtSpd}@\Omega,\text{V}_\text{rel}@V_\textrm{rel},\text{Pitch}@\theta,\text{Skew}@\chi\right) \f$
   !! - \f$ F_y = \frac{1}{2} \rho A \left( V_\textrm{rel,x} \right)^2 * C_\textrm{F,y}\left(\text{TSR}@\lambda,\text{RtSpd}@\Omega,\text{V}_\text{rel}@V_\textrm{rel},\text{Pitch}@\theta,\text{Skew}@\chi\right) \f$
   !! - \f$ F_z = \frac{1}{2} \rho A \left( V_\textrm{rel,x} \right)^2 * C_\textrm{F,z}\left(\text{TSR}@\lambda,\text{RtSpd}@\Omega,\text{V}_\text{rel}@V_\textrm{rel},\text{Pitch}@\theta,\text{Skew}@\chi\right) \f$
   !! - \f$ M_x = \frac{1}{2} \rho A \left( V_\textrm{rel,x} \right)^2 * C_\textrm{M,x}\left(\text{TSR}@\lambda,\text{RtSpd}@\Omega,\text{V}_\text{rel}@V_\textrm{rel},\text{Pitch}@\theta,\text{Skew}@\chi\right) \f$
   !! - \f$ M_y = \frac{1}{2} \rho A \left( V_\textrm{rel,x} \right)^2 * C_\textrm{M,y}\left(\text{TSR}@\lambda,\text{RtSpd}@\Omega,\text{V}_\text{rel}@V_\textrm{rel},\text{Pitch}@\theta,\text{Skew}@\chi\right) \f$
   !! - \f$ M_z = \frac{1}{2} \rho A \left( V_\textrm{rel,x} \right)^2 * C_\textrm{M,z}\left(\text{TSR}@\lambda,\text{RtSpd}@\Omega,\text{V}_\text{rel}@V_\textrm{rel},\text{Pitch}@\theta,\text{Skew}@\chi\right) \f$
   tmp1 = real(p%halfRhoA,SiKi) * m%VRel_xd * m%VRel_xd
   m%Force(1:3)  = tmp1 * m%C_F(1:3)
   m%Moment(1:3) = tmp1 * real(p%RotorRad,SiKi) * m%C_M(1:3)



   !-------------
   !> Set outputs
   !! - Thrust force:       \f$ C_t = C_\textrm{F,x} \f$
   !! - Torque coefficient: \f$ C_t = C_\textrm{M,x} \f$
   !! - \f$ \vec{F} = F_\textrm{x,disk} \hat{x} + F_\textrm{y,disk} \hat{y} + F_\textrm{z,disk} \hat{z} \f$ 
   !! - \f$ \vec{M} = M_\textrm{x,disk} \hat{x} + M_\textrm{y,disk} \hat{y} + M_\textrm{z,disk} \hat{z} \f$ 
   y%Ct = m%C_F(1)    ! Fx in disk reference frame
   y%Cq = m%C_M(1)    ! Mx in disk reference frame

   ! AeroLoads in global coordinates
   y%AeroLoads%Force( 1:3,1) = m%Force( 1)*m%x_hat + m%Force( 2)*m%y_hat + m%Force( 3)*m%z_hat
   y%AeroLoads%Moment(1:3,1) = m%Moment(1)*m%x_hat + m%Moment(2)*m%y_hat + m%Moment(3)*m%z_hat


   !----------------------------
   !> Set requested WriteOutputs
   call Calc_WriteOutput( u, p, y, m, ErrStat2, ErrMsg2, CalcWriteOutput )

   if (CalcWriteOutput) then
      ! Place the selected output channels into the WriteOutput(:)
      do i = 1,p%NumOuts  ! Loop through all selected output channels
         y%WriteOutput(i) = p%OutParam(i)%SignM * m%AllOuts( p%OutParam(i)%Indx )
      end do
   endif


contains
   logical function Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        Failed =  ErrStat >= AbortErrLev
        !if (Failed) call CleanUp()
   end function Failed
   subroutine CalcDiskAvgVel(ErrStat3,ErrMsg3)
      integer(IntKi),         intent(  out)  :: ErrStat3
      character(ErrMsgLen),   intent(  out)  :: ErrMsg3
      integer(IntKi) :: i
      integer(IntKi), parameter :: StartNode = 1  ! index to start returning wind info from external flow field (note: this will not work with ExtInflow)
      do i=1,ADsk_NumPtsDiskAvg
         m%DiskWindPosAbs(:,i) = real(u%HubMotion%Position(1:3,1)+u%HubMotion%TranslationDisp(1:3,1),ReKi) + matmul(real(u%HubMotion%Orientation(1:3,1:3,1),ReKi),p%DiskWindPosRel(:,i))
      end do
      call IfW_FlowField_GetVelAcc(p%FlowField, StartNode, t, m%DiskWindPosAbs, m%DiskWindVel, NoAcc, ErrStat3, ErrMsg3)
      if (ErrStat2 >= AbortErrLev) return
      ! calculate average
      m%DiskAvgVel = sum(m%DiskWindVel, dim=2) / REAL(ADsk_NumPtsDiskAvg,SiKi)
   end subroutine CalcDiskAvgVel
END SUBROUTINE ADsk_CalcOutput

subroutine ADskTableInterp(ATab, UseTSR, lambda, RotSpeed, VRel, BlPitch, Chi, idx_last, C_F, C_M, ErrStat, ErrMsg)
   type(ADsk_AeroTable),            intent(in   )  :: ATab           !< AeroTable
   logical,                         intent(in   )  :: UseTSR         !< flag to use TSR instead of RotSpeed and VRel
   real(SiKi),                      intent(in   )  :: lambda         !< TSR - tip speed ratio
   real(SiKi),                      intent(in   )  :: RotSpeed       !< Rotor Speed (rad/s)
   real(SiKi),                      intent(in   )  :: VRel           !< relative wind velocity along disk normal
   real(SiKi),                      intent(in   )  :: BlPitch        !< Blade pitch (collective)
   real(SiKi),                      intent(in   )  :: Chi            !< Inflow skew angle
   integer(IntKi),                  intent(inout)  :: idx_last(5)    !< m%idx_last -- for slight speedup in searching
   real(SiKi),                      intent(  out)  :: C_F(3)         !< Interpolated force coefficient
   real(SiKi),                      intent(  out)  :: C_M(3)         !< Interpolated moment coefficient
   integer(IntKi),                  intent(  out)  :: ErrStat        !< Error status of the operation
   character(*),                    intent(  out)  :: ErrMsg         !< Error message if ErrStat /= ErrID_None

   ! Local variables
   real(ReKi),parameter                            :: Tol = 1.0E-4   ! a tolerance for determining if two reals are the same (for interpolation)
   real(SiKi)                                      :: r_TSR          ! Location between bouding indices into TSR   dimension
   real(SiKi)                                      :: r_RtSpd        ! Location between bounding indices into RtSpd dimension
   real(SiKi)                                      :: r_VRel         ! Location between bounding indices into VRel  dimension
   real(SiKi)                                      :: r_Pitch        ! Location between bounding indices into Pitch dimension
   real(SiKi)                                      :: r_Skew         ! Location between bounding indices into Skew  dimension
   integer(IntKi)                                  :: i_TSR(2)       ! Bounding indices into TSR   dimension
   integer(IntKi)                                  :: i_RtSpd(2)     ! Bounding indices into RtSpd dimension
   integer(IntKi)                                  :: i_VRel(2)      ! Bounding indices into VRel  dimension
   integer(IntKi)                                  :: i_Pitch(2)     ! Bounding indices into Pitch dimension
   integer(IntKi)                                  :: i_Skew(2)      ! Bounding indices into Skew  dimension
   real(SiKi)                                      :: N3D( 8)        ! interpolation coefficients for 3D interp -- size 2^n
   real(SiKi)                                      :: U3D( 8)        ! interpolation values       for 3D interp  -- size 2^n
   real(SiKi)                                      :: N4D(16)        ! interpolation coefficients for 4D interp -- size 2^n
   real(SiKi)                                      :: U4D(16)        ! interpolation values       for 4D interp  -- size 2^n
   integer(IntKi)                                  :: ErrStat2       ! local error status
   character(ErrMsgLen)                            :: ErrMsg2        ! local error message
   character(*), parameter                         :: RoutineName = 'ADskTableInterp'

      ! Initialize variables
   ErrStat   = ErrID_None
   ErrMsg    = ""

   !---------------------------------------------
   ! Check that we can interpolate into the table
   if (UseTSR) then
      if ((lambda < ATab%TSR(1)-Tol) .or. (lambda > ATab%TSR(ATab%N_TSR)+Tol)) then
         ErrMsg2  = " TSR value of "//trim(Num2LStr(lambda))//" is outside bounds of aero table ["// &
               trim(Num2LStr(ATab%TSR(1)))//":"//trim(Num2LStr(ATab%TSR(ATab%N_TSR)))//"]"
         call SetErrStat(ErrID_Fatal,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      endif
   else
      if ((RotSpeed < ATab%RtSpd(1)-Tol) .or. (RotSpeed > ATab%RtSpd(ATab%N_RtSpd)+Tol)) then
         ErrMsg2  = " Rotor Speed value of "//trim(Num2LStr(RotSpeed))//" is outside bounds of aero table ["// &
               trim(Num2LStr(ATab%RtSpd(1)))//":"//trim(Num2LStr(ATab%RtSpd(ATab%N_RtSpd)))//"]"
         call SetErrStat(ErrID_Fatal,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      endif
      if ((VRel < ATab%VRel(1)-Tol) .or. (VRel > ATab%VRel(ATab%N_VRel)+Tol)) then
         ErrMsg2  = " VRel value of "//trim(Num2LStr(VRel))//" is outside bounds of aero table ["// &
               trim(Num2LStr(ATab%VRel(1)))//":"//trim(Num2LStr(ATab%VRel(ATab%N_VRel)))//"]"
         call SetErrStat(ErrID_Fatal,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      endif
   endif
   if (ATab%N_Pitch > 0_IntKi) then
      if ((BlPitch < ATab%Pitch(1)-Tol) .or. (BlPitch > ATab%Pitch(ATab%N_Pitch)+Tol)) then
         ErrMsg2  = " Blade pitch value of "//trim(Num2LStr(BlPitch))//" is outside bounds of aero table ["// &
               trim(Num2LStr(ATab%Pitch(1)))//":"//trim(Num2LStr(ATab%Pitch(ATab%N_Pitch)))//"]"
         call SetErrStat(ErrID_Fatal,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      endif
   endif
   if (ATab%N_Skew  > 0_IntKi) then
      if ((Chi < ATab%Skew(1)-Tol) .or. (Chi > ATab%Skew(ATab%N_Skew)+Tol)) then
         ErrMsg2  = " Skew value of "//trim(Num2LStr(Chi))//" is outside bounds of aero table ["// &
               trim(Num2LStr(ATab%Skew(1)))//":"//trim(Num2LStr(ATab%Skew(ATab%N_Skew)))//"]"
         call SetErrStat(ErrID_Fatal,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      endif
   endif
   if (ErrStat >= AbortErrLev)   return


   !-------------------------------
   ! Find indices for interpolation

   ! Initialize all indices to 1.  If a column is not used, both indices will be 1
   i_TSR   = 1_IntKi
   i_RtSpd = 1_IntKi
   i_VRel  = 1_IntKi
   i_Pitch = 1_IntKi
   i_Skew  = 1_IntKi
   ! Set the ratios for all indices to 0 (centered between indices)
   r_TSR   = 0.0_SiKi
   r_RtSpd = 0.0_SiKi
   r_VRel  = 0.0_SiKi
   r_Pitch = 0.0_SiKi
   r_Skew  = 0.0_SiKi

   ! Find indices to TSR or RtSpd + VRel entries
   if (UseTSR) then ! use TSR entry
      if (ATab%N_TSR > 1_IntKi) then
         ! Find lambda in the table (idx_last(1) is for TSR)
         CALL LocateStp( lambda, ATab%TSR, idx_last(1), ATab%N_TSR )
         i_TSR(1) = idx_last(1)                                         ! point before our point of interest
         i_TSR(2) = min(idx_last(1)+1,ATab%N_TSR)                       ! next point, but not out of bounds
         call CalcIsoparCoords( lambda, ATab%TSR(i_TSR(1)), ATab%TSR(i_TSR(2)), r_TSR )
      endif
   else  ! Use RtSpd and VRel instead
      if (ATab%N_RtSpd > 1_IntKi) then
         ! Find RotSpeed in the table (idx_last(2) is for RtSpd)
         CALL LocateStp( RotSpeed, ATab%RtSpd, idx_last(2), ATab%N_RtSpd )
         i_RtSpd(1) = idx_last(2)
         i_RtSpd(2) = min(idx_last(2)+1,ATab%N_RtSpd)
         call CalcIsoparCoords( RotSpeed, ATab%RtSpd(i_RtSpd(1)), ATab%RtSpd(i_RtSpd(2)), r_RtSpd )
      endif
      if (ATab%N_VRel > 1_IntKi) then
         ! Find VRel in the table (idx_last(3) is for VRel)
         CALL LocateStp( VRel, ATab%VRel, idx_last(3), ATab%N_VRel )
         i_VRel(1) = idx_last(3)
         i_VRel(2) = min(idx_last(3)+1,ATab%N_VRel)
         call CalcIsoparCoords( VRel, ATab%VRel(i_VRel(1)), ATab%VRel(i_VRel(2)), r_VRel )
      endif
   endif 

   ! Find indices to pitch 
   if (ATab%N_Pitch > 1_IntKi) then
      ! Find Pitch in the table (idx_last(4) is for Pitch)
      CALL LocateStp( BlPitch, ATab%Pitch, idx_last(4), ATab%N_Pitch )
      i_Pitch(1) = idx_last(4)
      i_Pitch(2) = min(idx_last(4)+1,ATab%N_Pitch)
      call CalcIsoparCoords( BlPitch, ATab%Pitch(i_Pitch(1)), ATab%Pitch(i_Pitch(2)), r_Pitch )
   endif

   ! Find indices to Skew 
   if (ATab%N_Skew > 1_IntKi) then
      ! Find Chi in the table (idx_last(5) is for Skew)
      CALL LocateStp( Chi, ATab%Skew, idx_last(5), ATab%N_Skew )
      i_Skew(1) = idx_last(5)
      i_Skew(2) = min(idx_last(5)+1,ATab%N_Skew)
      call CalcIsoparCoords( Chi, ATab%Skew(i_Skew(1)), ATab%Skew(i_Skew(2)), r_Skew )
   endif


   !------------------------------------------------
   ! Interpolate values from the coefficients tables
   !  For speed, the TSR and RtSpd + VRel cases are
   !  handled separately.  The table could be
   !  interpolated using a 5D interpolation, but 
   !  since it is known that some indices will not
   !  be needed (one ignored indice for RtSpd + VRel
   !  case, or two indices for TSR case), simpler
   !  interpolations can be used.

   if (UseTSR) then ! use TSR entry
      ! Coefficients -- same for all calculations
      N3D = getN3D()

      ! Force coefficients
      U3D = getU3D(ATab%C_Fx);      C_F(1) = sum(N3D * U3D)
      U3D = getU3D(ATab%C_Fy);      C_F(2) = sum(N3D * U3D)
      U3D = getU3D(ATab%C_Fz);      C_F(3) = sum(N3D * U3D)

      ! Moment coefficients
      U3D = getU3D(ATab%C_Mx);      C_M(1) = sum(N3D * U3D)
      U3D = getU3D(ATab%C_My);      C_M(2) = sum(N3D * U3D)
      U3D = getU3D(ATab%C_Mz);      C_M(3) = sum(N3D * U3D)
   else  ! Use RtSpd and VRel instead
      ! Coefficients -- same for all calculations
      N4D = getN4D()

      ! Force coefficients
      U4D = getU4D(ATab%C_Fx);      C_F(1) = sum(N4D * U4D)
      U4D = getU4D(ATab%C_Fy);      C_F(2) = sum(N4D * U4D)
      U4D = getU4D(ATab%C_Fz);      C_F(3) = sum(N4D * U4D)

      ! Moment coefficients
      U4D = getU4D(ATab%C_Mx);      C_M(1) = sum(N4D * U4D)
      U4D = getU4D(ATab%C_My);      C_M(2) = sum(N4D * U4D)
      U4D = getU4D(ATab%C_Mz);      C_M(3) = sum(N4D * U4D)
   endif

   return

contains
   logical function Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        Failed =  ErrStat >= AbortErrLev
        !if (Failed) call CleanUp()
   end function Failed
   function getN3D() result(Nr)        ! For when TSR is given
      real(SiKi)  :: Nr(8)
      Nr( 1) = ( 1.0_SiKi - r_TSR ) * ( 1.0_SiKi - r_Pitch ) * ( 1.0_SiKi - r_Skew )
      Nr( 2) = ( 1.0_SiKi - r_TSR ) * ( 1.0_SiKi - r_Pitch ) * ( 1.0_SiKi + r_Skew )
      Nr( 3) = ( 1.0_SiKi - r_TSR ) * ( 1.0_SiKi + r_Pitch ) * ( 1.0_SiKi - r_Skew )
      Nr( 4) = ( 1.0_SiKi - r_TSR ) * ( 1.0_SiKi + r_Pitch ) * ( 1.0_SiKi + r_Skew )
      Nr( 5) = ( 1.0_SiKi + r_TSR ) * ( 1.0_SiKi - r_Pitch ) * ( 1.0_SiKi - r_Skew )
      Nr( 6) = ( 1.0_SiKi + r_TSR ) * ( 1.0_SiKi - r_Pitch ) * ( 1.0_SiKi + r_Skew )
      Nr( 7) = ( 1.0_SiKi + r_TSR ) * ( 1.0_SiKi + r_Pitch ) * ( 1.0_SiKi - r_Skew )
      Nr( 8) = ( 1.0_SiKi + r_TSR ) * ( 1.0_SiKi + r_Pitch ) * ( 1.0_SiKi + r_Skew )
      Nr     = Nr/ REAL( SIZE(Nr), SiKi )  ! normalize
   end function getN3D
   function getN4D() result(Nr)        ! For when TSR is not given
      real(SiKi)  :: Nr(16)
      Nr( 1) = ( 1.0_SiKi - r_RtSpd ) * ( 1.0_SiKi - r_VRel ) * ( 1.0_SiKi - r_Pitch ) * ( 1.0_SiKi - r_Skew )
      Nr( 2) = ( 1.0_SiKi - r_RtSpd ) * ( 1.0_SiKi - r_VRel ) * ( 1.0_SiKi - r_Pitch ) * ( 1.0_SiKi + r_Skew )
      Nr( 3) = ( 1.0_SiKi - r_RtSpd ) * ( 1.0_SiKi - r_VRel ) * ( 1.0_SiKi + r_Pitch ) * ( 1.0_SiKi - r_Skew )
      Nr( 4) = ( 1.0_SiKi - r_RtSpd ) * ( 1.0_SiKi - r_VRel ) * ( 1.0_SiKi + r_Pitch ) * ( 1.0_SiKi + r_Skew )
      Nr( 5) = ( 1.0_SiKi - r_RtSpd ) * ( 1.0_SiKi + r_VRel ) * ( 1.0_SiKi - r_Pitch ) * ( 1.0_SiKi - r_Skew )
      Nr( 6) = ( 1.0_SiKi - r_RtSpd ) * ( 1.0_SiKi + r_VRel ) * ( 1.0_SiKi - r_Pitch ) * ( 1.0_SiKi + r_Skew )
      Nr( 7) = ( 1.0_SiKi - r_RtSpd ) * ( 1.0_SiKi + r_VRel ) * ( 1.0_SiKi + r_Pitch ) * ( 1.0_SiKi - r_Skew )
      Nr( 8) = ( 1.0_SiKi - r_RtSpd ) * ( 1.0_SiKi + r_VRel ) * ( 1.0_SiKi + r_Pitch ) * ( 1.0_SiKi + r_Skew )
      Nr( 9) = ( 1.0_SiKi + r_RtSpd ) * ( 1.0_SiKi - r_VRel ) * ( 1.0_SiKi - r_Pitch ) * ( 1.0_SiKi - r_Skew )
      Nr(10) = ( 1.0_SiKi + r_RtSpd ) * ( 1.0_SiKi - r_VRel ) * ( 1.0_SiKi - r_Pitch ) * ( 1.0_SiKi + r_Skew )
      Nr(11) = ( 1.0_SiKi + r_RtSpd ) * ( 1.0_SiKi - r_VRel ) * ( 1.0_SiKi + r_Pitch ) * ( 1.0_SiKi - r_Skew )
      Nr(12) = ( 1.0_SiKi + r_RtSpd ) * ( 1.0_SiKi - r_VRel ) * ( 1.0_SiKi + r_Pitch ) * ( 1.0_SiKi + r_Skew )
      Nr(13) = ( 1.0_SiKi + r_RtSpd ) * ( 1.0_SiKi + r_VRel ) * ( 1.0_SiKi - r_Pitch ) * ( 1.0_SiKi - r_Skew )
      Nr(14) = ( 1.0_SiKi + r_RtSpd ) * ( 1.0_SiKi + r_VRel ) * ( 1.0_SiKi - r_Pitch ) * ( 1.0_SiKi + r_Skew )
      Nr(15) = ( 1.0_SiKi + r_RtSpd ) * ( 1.0_SiKi + r_VRel ) * ( 1.0_SiKi + r_Pitch ) * ( 1.0_SiKi - r_Skew )
      Nr(16) = ( 1.0_SiKi + r_RtSpd ) * ( 1.0_SiKi + r_VRel ) * ( 1.0_SiKi + r_Pitch ) * ( 1.0_SiKi + r_Skew )
      Nr     = Nr / REAL( SIZE(Nr), SiKi )  ! normalize
   end function getN4D
   function getU3D(CT) result(Ur)      ! For when TSR is given (i_RtSpd(1)=i_VRel(1)=1)
      real(SiKi),  intent(in) :: CT(:,:,:,:,:)      ! Coefficient table
      real(SiKi)              :: Ur(8)
      Ur( 1) = CT( i_TSR(1), i_RtSpd(1), i_VRel(1), i_Pitch(1), i_Skew(1) )
      Ur( 2) = CT( i_TSR(1), i_RtSpd(1), i_VRel(1), i_Pitch(1), i_Skew(2) )
      Ur( 3) = CT( i_TSR(1), i_RtSpd(1), i_VRel(1), i_Pitch(2), i_Skew(1) )
      Ur( 4) = CT( i_TSR(1), i_RtSpd(1), i_VRel(1), i_Pitch(2), i_Skew(2) )
      Ur( 5) = CT( i_TSR(2), i_RtSpd(1), i_VRel(1), i_Pitch(1), i_Skew(1) )
      Ur( 6) = CT( i_TSR(2), i_RtSpd(1), i_VRel(1), i_Pitch(1), i_Skew(2) )
      Ur( 7) = CT( i_TSR(2), i_RtSpd(1), i_VRel(1), i_Pitch(2), i_Skew(1) )
      Ur( 8) = CT( i_TSR(2), i_RtSpd(1), i_VRel(1), i_Pitch(2), i_Skew(2) )
   end function getU3D
   function getU4D(CT) result(Ur)      ! For when TSR is not given (i_TSR(1)=1)
      real(SiKi),  intent(in) :: CT(:,:,:,:,:)      ! Coefficient table
      real(SiKi)              :: Ur(16)
      Ur( 1) = CT( i_TSR(1), i_RtSpd(1), i_VRel(1), i_Pitch(1), i_Skew(1) )
      Ur( 2) = CT( i_TSR(1), i_RtSpd(1), i_VRel(1), i_Pitch(1), i_Skew(2) )
      Ur( 3) = CT( i_TSR(1), i_RtSpd(1), i_VRel(1), i_Pitch(2), i_Skew(1) )
      Ur( 4) = CT( i_TSR(1), i_RtSpd(1), i_VRel(1), i_Pitch(2), i_Skew(2) )
      Ur( 5) = CT( i_TSR(1), i_RtSpd(1), i_VRel(2), i_Pitch(1), i_Skew(1) )
      Ur( 6) = CT( i_TSR(1), i_RtSpd(1), i_VRel(2), i_Pitch(1), i_Skew(2) )
      Ur( 7) = CT( i_TSR(1), i_RtSpd(1), i_VRel(2), i_Pitch(2), i_Skew(1) )
      Ur( 8) = CT( i_TSR(1), i_RtSpd(1), i_VRel(2), i_Pitch(2), i_Skew(2) )
      Ur( 9) = CT( i_TSR(1), i_RtSpd(2), i_VRel(1), i_Pitch(1), i_Skew(1) )
      Ur(10) = CT( i_TSR(1), i_RtSpd(2), i_VRel(1), i_Pitch(1), i_Skew(2) )
      Ur(11) = CT( i_TSR(1), i_RtSpd(2), i_VRel(1), i_Pitch(2), i_Skew(1) )
      Ur(12) = CT( i_TSR(1), i_RtSpd(2), i_VRel(1), i_Pitch(2), i_Skew(2) )
      Ur(13) = CT( i_TSR(1), i_RtSpd(2), i_VRel(2), i_Pitch(1), i_Skew(1) )
      Ur(14) = CT( i_TSR(1), i_RtSpd(2), i_VRel(2), i_Pitch(1), i_Skew(2) )
      Ur(15) = CT( i_TSR(1), i_RtSpd(2), i_VRel(2), i_Pitch(2), i_Skew(1) )
      Ur(16) = CT( i_TSR(1), i_RtSpd(2), i_VRel(2), i_Pitch(2), i_Skew(2) )
   end function getU4D
end subroutine ADskTableInterp
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine calculates the iosparametric coordinates, isopc, which is a value between -1 and 1 (for each dimension of a dataset)
!! indicating where InCoord falls between posLo and posHi.
!! This routine is copied from WAMIT_Interp.f90
subroutine CalcIsoparCoords( InCoord, posLo, posHi, isopc )
   real(SiKi),    intent(in   )  :: InCoord     !< 
   real(SiKi),    intent(in   )  :: posLo       !< coordinate values associated with Indx_Lo 
   real(SiKi),    intent(in   )  :: posHi       !< coordinate values associated with Indx_Hi
   real(SiKi),    intent(  out)  :: isopc       !< isoparametric coordinates 
   ! local variables
   real(SiKi)                     :: dx         ! difference between high and low coordinates in the bounding "box"
   dx = posHi - posLo 
   if (EqualRealNos(dx, 0.0_SiKi)) then
      isopc = 1.0_SiKi
   else
      isopc = ( 2.0_SiKi*InCoord - posLo - posHi ) / dx
         ! to verify that we don't extrapolate, make sure this is bound between -1 and 1 (effectively nearest neighbor)
      isopc = min( 1.0_SiKi, isopc )
      isopc = max(-1.0_SiKi, isopc )
   end if
end subroutine CalcIsoparCoords

!----------------------------------------------------------------------------------------------------------------------------------
!> This is a loose coupling routine for solving constraint states, integrating continuous states, and updating discrete and other
!! states. Continuous, constraint, discrete, and other states are updated to values at t + Interval.
!! NOTE: there are no states in AeroDisk
SUBROUTINE ADsk_UpdateStates( t, n, Inputs, InputTimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
   real(DbKi),                         intent(in   )  :: t               !< Current simulation time in seconds
   integer(IntKi),                     intent(in   )  :: n               !< Current step of the simulation: t = n*Interval
   type(ADsk_InputType),               intent(inout)  :: Inputs(:)       !< Inputs at InputTimes (output for mesh connect)
   real(DbKi),                         intent(in   )  :: InputTimes(:)   !< Times in seconds associated with Inputs
   type(ADsk_ParameterType),           intent(in   )  :: p               !< Parameters
   type(ADsk_ContinuousStateType),     intent(inout)  :: x               !< Input: Continuous states at t;
   type(ADsk_DiscreteStateType),       intent(inout)  :: xd              !< Input: Discrete states at t;
   type(ADsk_ConstraintStateType),     intent(inout)  :: z               !< Input: Constraint states at t;
   type(ADsk_OtherStateType),          intent(inout)  :: OtherState      !< Other states: Other states at t;
   type(ADsk_MiscVarType),             intent(inout)  :: m               !<  Misc variables for optimization
   integer(IntKi),                     intent(  out)  :: ErrStat         !< Error status of the operation
   character(*),                       intent(  out)  :: ErrMsg          !< Error message if ErrStat /= ErrID_None

   ! Local variables
   character(*), parameter                            :: RoutineName = 'ADsk_UpdateStates'

      ! Initialize variables
   ErrStat   = ErrID_None
   ErrMsg    = ""

      ! There are no states.
   x%DummyContState = 0.0_ReKi

end subroutine ADsk_UpdateStates


!----------------------------------------------------------------------------------------------------------------------------------
!> This is a tight coupling routine for computing derivatives of continuous states.
!! NOTE: there are no states in AeroDisk
SUBROUTINE ADsk_CalcContStateDeriv( t, u, p, x, xd, z, OtherState, m, dxdt, ErrStat, ErrMsg )
   real(DbKi),                      intent(in   )  :: t           !< Current simulation time in seconds
   type(ADsk_InputType),            intent(in   )  :: u           !< Inputs at t
   type(ADsk_ParameterType),        intent(in   )  :: p           !< Parameters
   type(ADsk_ContinuousStateType),  intent(in   )  :: x           !< Continuous states at t
   type(ADsk_DiscreteStateType),    intent(in   )  :: xd          !< Discrete states at t
   type(ADsk_ConstraintStateType),  intent(in   )  :: z           !< Constraint states at t
   type(ADsk_OtherStateType),       intent(in   )  :: OtherState  !< Other states at t
   type(ADsk_MiscVarType),          intent(inout)  :: m           !< Misc variables for optimization (not copied in glue code)
   type(ADsk_ContinuousStateType),  intent(  out)  :: dxdt        !< Continuous state derivatives at t
   INTEGER(IntKi),                  intent(  out)  :: ErrStat     !< Error status of the operation
   character(*),                    intent(  out)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

   ! local variables
   character(*), parameter                         :: RoutineName = 'ADsk_CalcContStateDeriv'

      ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ""

      ! There are no states
   dxdt%DummyContState = 0.0_ReKi

END SUBROUTINE ADsk_CalcContStateDeriv


END MODULE AeroDisk
!**********************************************************************************************************************************
