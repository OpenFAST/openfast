!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2015-2016  National Renewable Energy Laboratory
! Copyright (C) 2016-2017  Envision Energy USA, LTD
!
!    This file is part of AeroDyn.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!
!**********************************************************************************************************************************
module BEMT
    
   use NWTC_Library
   
   use BEMT_Types
   use BEMTUncoupled
   use DBEMT
   
   use UnsteadyAero
   !USE AeroDyn_Types
   use AirfoilInfo
   

   implicit none
         
   
   private
   
   type(ProgDesc), parameter  :: BEMT_Ver = ProgDesc( 'BEM', '', '' )
   character(*),   parameter  :: BEMT_Nickname = 'BEM'
      
   
   ! ..... Public Subroutines ...................................................................................................

   public :: BEMT_Init                           ! Initialization routine
   public :: BEMT_End                            ! Ending routine (includes clean up)

   public :: BEMT_UpdateStates                   ! Loose coupling routine for solving for constraint states, integrating
                                                 !   continuous states, and updating discrete states
   public :: BEMT_CalcOutput                     ! Routine for computing outputs

   public :: BEMT_CalcConstrStateResidual        ! Tight coupling routine for returning the constraint state residual
   public :: BEMT_CalcContStateDeriv             ! Tight coupling routine for computing derivatives of continuous states
   public :: BEMT_UpdateDiscState                ! Tight coupling routine for updating discrete states
   
   public :: BEMT_ReInit
   ! routines for linearization
   public :: Get_phi_perturbations
   public :: ComputeFrozenWake
   public :: CheckLinearizationInput
   public :: UpdatePhi
   public :: BEMT_InitStates
   
   contains


    
!----------------------------------------------------------------------------------------------------------------------------------   
real(ReKi) function ComputePhiWithInduction( Vx, Vy, a, aprime )
! This routine is used to compute the inflow angle, phi, from the local velocities and the induction factors.
!..................................................................................................................................
   real(ReKi),                    intent(in   )  :: Vx          ! Local velocity component along the thrust direction
   real(ReKi),                    intent(in   )  :: Vy          ! Local velocity component along the rotor plane-of-rotation direction
   real(ReKi),                    intent(in   )  :: a           ! Axial induction factor
   real(ReKi),                    intent(in   )  :: aprime      ! Tangential induction factor
   
   real(ReKi)                                    :: x
   real(ReKi)                                    :: y
      
   x = Vx*(1.0_ReKi-a)
   y = Vy*(1.0_ReKi + aprime)
   
   if ( EqualRealNos(y, 0.0_ReKi) .AND. EqualRealNos(x, 0.0_ReKi) ) then
      ComputePhiWithInduction = 0.0_ReKi
   else
      ComputePhiWithInduction  = atan2( x , y )
   end if
   
   
end function ComputePhiWithInduction
 
!----------------------------------------------------------------------------------------------------------------------------------   
subroutine BEMT_Set_UA_InitData( InitInp, interval, Init_UA_Data, errStat, errMsg )
! This routine is called from BEMT_Init.
! The parameters are set here and not changed during the simulation.
!..................................................................................................................................
   type(BEMT_InitInputType),       intent(inout)  :: InitInp     ! Input data for initialization routine
   real(DbKi),                     intent(in   )  :: interval    ! time interval  
   type(UA_InitInputType),         intent(  out)  :: Init_UA_Data           ! Parameters
   integer(IntKi),                 intent(  out)  :: errStat     ! Error status of the operation
   character(*),                   intent(  out)  :: errMsg      ! Error message if ErrStat /= ErrID_None

   integer                                        :: i,j
   integer(intKi)                                 :: ErrStat2           ! temporary Error status
   
      ! Set up initialization data
   
   Allocate(Init_UA_Data%c(InitInp%numBladeNodes,InitInp%numBlades), STAT = errStat2)
   if (ErrStat2 /= 0) then
      ErrStat = ErrID_Fatal
      ErrMsg = "BEMT_Set_UA_InitData:Error allocating Init_UA_Data%c."
      return
   else
      ErrStat = ErrID_None
      ErrMsg = ""
   end if
   
   do j = 1,InitInp%numBlades
      do i = 1,InitInp%numBladeNodes
         Init_UA_Data%c(i,j)      = InitInp%chord(i,j)
      end do
   end do
   
   call move_alloc(InitInp%UAOff_innerNode, Init_UA_Data%UAOff_innerNode)
   call move_alloc(InitInp%UAOff_outerNode, Init_UA_Data%UAOff_outerNode)
   
   Init_UA_Data%dt              = interval          
   Init_UA_Data%OutRootName     = InitInp%RootName ! was 'Debug.UA'
               
   Init_UA_Data%numBlades       = InitInp%numBlades 
   Init_UA_Data%nNodesPerBlade  = InitInp%numBladeNodes
                                  
   Init_UA_Data%UAMod           = InitInp%UAMod  
   Init_UA_Data%Flookup         = InitInp%Flookup
   Init_UA_Data%a_s             = InitInp%a_s ! m/s  
   Init_UA_Data%ShedEffect      = .true. ! This should be true when coupled to BEM
   Init_UA_Data%WrSum           = InitInp%SumPrint
   
end subroutine BEMT_Set_UA_InitData

   
!----------------------------------------------------------------------------------------------------------------------------------   
subroutine BEMT_SetParameters( InitInp, p, errStat, errMsg )
! This routine is called from BEMT_Init.
! The parameters are set here and not changed during the simulation.
!..................................................................................................................................
   type(BEMT_InitInputType),       intent(in   )  :: InitInp     ! Input data for initialization routine
   type(BEMT_ParameterType),       intent(  out)  :: p           ! Parameters
   integer(IntKi),                 intent(  out)  :: errStat     ! Error status of the operation
   character(*),                   intent(  out)  :: errMsg      ! Error message if ErrStat /= ErrID_None


      ! Local variables
   integer(IntKi)                                :: errStat2                ! temporary Error status of the operation
   character(*), parameter                       :: RoutineName = 'BEMT_SetParameters'
   integer(IntKi)                                :: i, j

   
      ! Initialize variables for this routine

   errStat = ErrID_None
   errMsg  = ""

   p%numBladeNodes  = InitInp%numBladeNodes 
   p%numBlades      = InitInp%numBlades    
   p%UA_Flag        = InitInp%UA_Flag   
   p%DBEMT_Mod      = InitInp%DBEMT_Mod

   allocate ( p%chord(p%numBladeNodes, p%numBlades), STAT = errStat2 )
   if ( errStat2 /= 0 ) then
      call SetErrStat( ErrID_Fatal, 'Error allocating memory for p%chord.', errStat, errMsg, RoutineName )
      return
   end if 
   
   allocate ( p%zHub(p%numBlades), STAT = errStat2 )
   if ( errStat2 /= 0 ) then
      call SetErrStat( ErrID_Fatal, 'Error allocating memory for p%zHub.', errStat, errMsg, RoutineName )
      return
   end if 
   
   allocate ( p%AFindx(p%numBladeNodes,p%numBlades), STAT = errStat2 )
   if ( errStat2 /= 0 ) then
      call SetErrStat( ErrID_Fatal, 'Error allocating memory for p%AFindx.', errStat, errMsg, RoutineName )
      return
   end if 
   
   allocate ( p%tipLossConst(p%numBladeNodes, p%numBlades), STAT = errStat2 )
   if ( errStat2 /= 0 ) then
      call SetErrStat( ErrID_Fatal, 'Error allocating memory for p%tipLossConst.', errStat, errMsg, RoutineName )
      return
   end if 
   
   allocate ( p%hubLossConst(p%numBladeNodes, p%numBlades), STAT = errStat2 )
   if ( errStat2 /= 0 ) then
      call SetErrStat( ErrID_Fatal, 'Error allocating memory for p%hubLossConst.', errStat, errMsg, RoutineName )
      return
   end if 
   
   allocate ( p%FixedInductions(p%numBladeNodes, p%numBlades), STAT = errStat2 )
   if ( errStat2 /= 0 ) then
      call SetErrStat( ErrID_Fatal, 'Error allocating memory for p%numBladeNodes.', errStat, errMsg, RoutineName )
      return
   end if
   
   p%AFindx = InitInp%AFindx 
   
      ! Compute the tip and hub loss constants using the distances along the blade (provided as input for now) 
   do j=1,p%numBlades
      p%zHub(j) = InitInp%zHub(j)
      do i=1,p%numBladeNodes
         p%chord(i,j)        = InitInp%chord(i,j)
         p%tipLossConst(i,j) = p%numBlades*(InitInp%zTip    (j) - InitInp%zLocal(i,j)) / (2.0*InitInp%zLocal(i,j))
         p%hubLossConst(i,j) = p%numBlades*(InitInp%zLocal(i,j) - InitInp%zHub    (j)) / (2.0*InitInp%zHub    (j))
      end do
   end do
   
   
  !p%DT               = InitInp%DT
   p%airDens          = InitInp%airDens
   p%kinVisc          = InitInp%kinVisc
   p%skewWakeMod      = InitInp%skewWakeMod
   p%yawCorrFactor    = InitInp%yawCorrFactor
   p%useTipLoss       = InitInp%useTipLoss
   p%useHubLoss       = InitInp%useHubLoss
   p%useInduction     = InitInp%useInduction
   p%useTanInd        = InitInp%useTanInd
   p%useAIDrag        = InitInp%useAIDrag
   p%useTIDrag        = InitInp%useTIDrag
   p%numReIterations  = InitInp%numReIterations
   p%maxIndIterations = InitInp%maxIndIterations
   p%aTol             = InitInp%aTol
   
   
   ! setting this condition here so we don't have to do some many EqualRealNos() checks later in the code.
   do j=1,p%numBlades
      do i=1,p%numBladeNodes
         p%FixedInductions(i,j) = ( p%useTiploss .and. EqualRealNos(p%tipLossConst(i,j),0.0_ReKi) ) .or. ( p%useHubloss .and. EqualRealNos(p%hubLossConst(i,j),0.0_ReKi) )
      end do
   end do
   
end subroutine BEMT_SetParameters

!----------------------------------------------------------------------------------------------------------------------------------
subroutine BEMT_InitContraintStates( z, p, errStat, errMsg )
! This routine is called from BEMT_Init.
! The constraint state data is allocated and set to zero.
!..................................................................................................................................

   type(BEMT_ConstraintStateType), intent(  out)  :: z           ! Input data for initialization routine
   type(BEMT_ParameterType),       intent(in   )  :: p           ! Parameters
   integer(IntKi),                 intent(  out)  :: errStat     ! Error status of the operation
   character(*),                   intent(  out)  :: errMsg      ! Error message if ErrStat /= ErrID_None


      ! Local variables
   character(*), parameter                       :: RoutineName = 'BEMT_InitContraintStates'
   integer(IntKi)                                :: errStat2                ! temporary Error status of the operation
   
      ! Initialize variables for this routine

   errStat = ErrID_None
   errMsg  = ""
      
   allocate ( z%phi( p%numBladeNodes, p%numBlades ), STAT = errStat2 )
   if ( errStat2 /= 0 ) then
      call SetErrStat( ErrID_Fatal, 'Error allocating memory for z%phi.', errStat, errMsg, RoutineName )
      return
   end if 
  
end subroutine BEMT_InitContraintStates


!----------------------------------------------------------------------------------------------------------------------------------
subroutine BEMT_InitOtherStates( OtherState, p, errStat, errMsg )
! This routine is called from BEMT_Init.
! The OtherState data is allocated and set to zero.
!..................................................................................................................................

   type(BEMT_OtherStateType),      intent(inout)  :: OtherState  ! OtherState data
   type(BEMT_ParameterType),       intent(in   )  :: p           ! Parameters
   integer(IntKi),                intent(  out)  :: errStat     ! Error status of the operation
   character(*),                  intent(  out)  :: errMsg      ! Error message if ErrStat /= ErrID_None


      ! Local variables
   !character(ErrMsgLen)                          :: errMsg2                 ! temporary Error message if ErrStat /= ErrID_None
   integer(IntKi)                                :: errStat2                ! temporary Error status of the operation
   character(*), parameter                       :: RoutineName = 'BEMT_InitOtherStates'
   
      ! Initialize variables for this routine

   errStat = ErrID_None
   errMsg  = ""
   
   
   if (p%UseInduction) then
      
      allocate ( OtherState%ValidPhi( p%numBladeNodes, p%numBlades ), STAT = errStat2 )
      if ( errStat2 /= 0 ) then
         call SetErrStat( ErrID_Fatal, 'Error allocating memory for OtherState%ValidPhi.', errStat, errMsg, RoutineName )
         return
      end if

   end if
   
   ! values of the OtherStates are initialized in BEMT_ReInit()

end subroutine BEMT_InitOtherStates

!----------------------------------------------------------------------------------------------------------------------------------
subroutine BEMT_AllocInput( u, p, errStat, errMsg )
! This routine is called from BEMT_Init.
!  
!  
!..................................................................................................................................

   type(BEMT_InputType),           intent(  out)  :: u           ! Input data
   type(BEMT_ParameterType),       intent(in   )  :: p           ! Parameters
   integer(IntKi),                 intent(  out)  :: errStat     ! Error status of the operation
   character(*),                   intent(  out)  :: errMsg      ! Error message if ErrStat /= ErrID_None


      ! Local variables
   integer(IntKi)                                :: errStat2                ! temporary Error status of the operation
   character(*), parameter                       :: RoutineName = 'BEMT_AllocInput'

      ! Initialize variables for this routine

   errStat = ErrID_None
   errMsg  = ""

   allocate ( u%theta( p%numBladeNodes, p%numBlades ), STAT = errStat2 )
   if ( errStat2 /= 0 ) then
      call SetErrStat( ErrID_Fatal, 'Error allocating memory for u%theta.', errStat, errMsg, RoutineName )
      return
   end if 
   u%theta = 0.0_ReKi
   
   allocate ( u%psi( p%numBlades ), STAT = errStat2 )
   if ( errStat2 /= 0 ) then
      call SetErrStat( ErrID_Fatal, 'Error allocating memory for u%psi.', errStat, errMsg, RoutineName )
      return
   end if 
   u%psi = 0.0_ReKi
   
   allocate ( u%Vx( p%numBladeNodes, p%numBlades ), STAT = errStat2 )
   if ( errStat2 /= 0 ) then
      call SetErrStat( ErrID_Fatal, 'Error allocating memory for u%Vx.', errStat, errMsg, RoutineName )
      return
   end if 
   u%Vx = 0.0_ReKi
   
   allocate ( u%Vy( p%numBladeNodes, p%numBlades ), STAT = errStat2 )
   if ( errStat2 /= 0 ) then
      call SetErrStat( ErrID_Fatal, 'Error allocating memory for u%Vy.', errStat, errMsg, RoutineName )
      return
   end if 
   u%Vy = 0.0_ReKi
 
   if (p%DBEMT_Mod==DBEMT_cont_tauConst) then
      allocate ( u%Vx_elast_dot( p%numBladeNodes, p%numBlades ), STAT = errStat2 )
      if ( errStat2 /= 0 ) then
         call SetErrStat( ErrID_Fatal, 'Error allocating memory for u%Vx_dot.', errStat, errMsg, RoutineName )
         return
      end if 
      u%Vx_elast_dot = 0.0_ReKi
   
      allocate ( u%Vy_elast_dot( p%numBladeNodes, p%numBlades ), STAT = errStat2 )
      if ( errStat2 /= 0 ) then
         call SetErrStat( ErrID_Fatal, 'Error allocating memory for u%Vy_dot.', errStat, errMsg, RoutineName )
         return
      end if 
      u%Vy_elast_dot = 0.0_ReKi
   end if
   
   allocate ( u%omega_z( p%numBladeNodes, p%numBlades ), STAT = errStat2 )
   if ( errStat2 /= 0 ) then
      call SetErrStat( ErrID_Fatal, 'Error allocating memory for u%omega_z.', errStat, errMsg, RoutineName )
      return
   end if 
   u%omega_z = 0.0_ReKi
   
   allocate ( u%rLocal( p%numBladeNodes, p%numBlades ), STAT = errStat2 )
   if ( errStat2 /= 0 ) then
      call SetErrStat( ErrID_Fatal, 'Error allocating memory for u%rLocal.', errStat, errMsg, RoutineName )
      return
   end if 
   u%rLocal = 0.0_ReKi
   
   allocate ( u%UserProp( p%numBladeNodes, p%numBlades ), STAT = errStat2 )
   if ( errStat2 /= 0 ) then
      call SetErrStat( ErrID_Fatal, 'Error allocating memory for u%UserProp.', errStat, errMsg, RoutineName )
      return
   end if 
   u%UserProp = 0.0_ReKi
   
   
   u%omega  = 0.0_ReKi
   
end subroutine BEMT_AllocInput


!----------------------------------------------------------------------------------------------------------------------------------
subroutine BEMT_AllocOutput( y, p, errStat, errMsg )
! This routine is called from BEMT_Init.
!  
!  
!..................................................................................................................................

   type(BEMT_OutputType),         intent(  out)  :: y           ! output data
   type(BEMT_ParameterType),      intent(in   )  :: p           ! Parameters
   integer(IntKi),                intent(  out)  :: errStat     ! Error status of the operation
   character(*),                  intent(  out)  :: errMsg      ! Error message if ErrStat /= ErrID_None


      ! Local variables
   character(ErrMsgLen  )                        :: errMsg2                 ! temporary Error message if ErrStat /= ErrID_None
   integer(IntKi)                                :: errStat2                ! temporary Error status of the operation
   character(*), parameter                       :: RoutineName = 'BEMT_AllocOutput'
   
      ! Initialize variables for this routine

   errStat = ErrID_None
   errMsg  = ""

   call allocAry( y%Vrel, p%numBladeNodes, p%numBlades, 'y%Vrel', errStat2, errMsg2); call setErrStat(errStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call allocAry( y%phi, p%numBladeNodes, p%numBlades, 'y%phi', errStat2, errMsg2); call setErrStat(errStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call allocAry( y%chi, p%numBladeNodes, p%numBlades, 'y%chi', errStat2, errMsg2); call setErrStat(errStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call allocAry( y%Re, p%numBladeNodes, p%numBlades, 'y%Re', errStat2, errMsg2); call setErrStat(errStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call allocAry( y%axInduction, p%numBladeNodes, p%numBlades, 'y%axInduction', errStat2, errMsg2); call setErrStat(errStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call allocAry( y%tanInduction, p%numBladeNodes, p%numBlades, 'y%tanInduction', errStat2, errMsg2); call setErrStat(errStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call allocAry( y%AOA, p%numBladeNodes, p%numBlades, 'y%AOA', errStat2, errMsg2); call setErrStat(errStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call allocAry( y%Cx, p%numBladeNodes, p%numBlades, 'y%Cx', errStat2, errMsg2); call setErrStat(errStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call allocAry( y%Cy, p%numBladeNodes, p%numBlades, 'y%Cy', errStat2, errMsg2); call setErrStat(errStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call allocAry( y%Cm, p%numBladeNodes, p%numBlades, 'y%Cm', errStat2, errMsg2); call setErrStat(errStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call allocAry( y%Cl, p%numBladeNodes, p%numBlades, 'y%Cl', errStat2, errMsg2); call setErrStat(errStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call allocAry( y%Cd, p%numBladeNodes, p%numBlades, 'y%Cd', errStat2, errMsg2); call setErrStat(errStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call allocAry( y%Cpmin, p%numBladeNodes, p%numBlades, 'm%Cpmin', errStat2, errMsg2); call setErrStat(errStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

   
   if (ErrStat >= AbortErrLev) RETURN
   
      ! outputs documented in AeroDyn
   y%Vrel = 0.0_ReKi
   y%phi = 0.0_ReKi   
   y%Cx = 0.0_ReKi
   y%Cy = 0.0_ReKi
   y%Cm = 0.0_ReKi

      ! others:
   y%chi = 0.0_ReKi
   y%Re = 0.0_ReKi   
   y%axInduction = 0.0_ReKi   
   y%tanInduction = 0.0_ReKi
   y%AOA = 0.0_ReKi   
   y%Cl = 0.0_ReKi
   y%Cd = 0.0_ReKi
   y%Cpmin = 0.0_ReKi
end subroutine BEMT_AllocOutput


!----------------------------------------------------------------------------------------------------------------------------------
subroutine BEMT_Init( InitInp, u, p, x, xd, z, OtherState, AFInfo, y, misc, Interval, InitOut, ErrStat, ErrMsg )
! This routine is called at the start of the simulation to perform initialization steps.
! The parameters are set here and not changed during the simulation.
! The initial states and initial guess for the input are defined.
!..................................................................................................................................

   type(BEMT_InitInputType),       intent(inout)  :: InitInp     ! Input data for initialization routine; intent out so that we can use MOVE_ALLOC
   type(BEMT_InputType),           intent(  out)  :: u           ! An initial guess for the input; input mesh must be defined
   type(BEMT_ParameterType),       intent(  out)  :: p           ! Parameters
   type(BEMT_ContinuousStateType), intent(  out)  :: x           ! Initial continuous states
   type(BEMT_DiscreteStateType),   intent(  out)  :: xd          ! Initial discrete states
   type(BEMT_ConstraintStateType), intent(  out)  :: z           ! Initial guess of the constraint states
   type(BEMT_OtherStateType),      intent(  out)  :: OtherState  ! Initial other states
   type(BEMT_MiscVarType),         intent(  out)  :: misc        ! Initial misc/optimization variables
   type(AFI_ParameterType),        intent(in   )  :: AFInfo(:)   ! The airfoil parameter data
   type(BEMT_OutputType),          intent(  out)  :: y           ! Initial system outputs (outputs are not calculated;
                                                                 !   only the output mesh is initialized)
   real(DbKi),                     intent(in   )  :: interval    ! Coupling interval in seconds: the rate that
                                                                 !   (1) BEMT_UpdateStates() is called in loose coupling &
                                                                 !   (2) BEMT_UpdateDiscState() is called in tight coupling.
                                                                 !   Input is the suggested time from the glue code;
                                                                 !   Output is the actual coupling interval that will be used
                                                                 !   by the glue code.
   type(BEMT_InitOutputType),      intent(  out)  :: InitOut     ! Output for initialization routine
   integer(IntKi),                 intent(  out)  :: errStat     ! Error status of the operation
   character(*),                   intent(  out)  :: errMsg      ! Error message if ErrStat /= ErrID_None


      ! Local variables
   character(ErrMsgLen)                           :: errMsg2     ! temporary Error message if ErrStat /= ErrID_None
   integer(IntKi)                                 :: errStat2    ! temporary Error status of the operation
   character(*), parameter                        :: RoutineName = 'BEMT_Init'
   type(UA_InitInputType)                         :: Init_UA_Data
   type(UA_InitOutputType)                        :: InitOutData_UA

   type(DBEMT_InitInputType)                      :: InitInp_DBEMT
   type(DBEMT_InitOutputType)                     :: InitOut_DBEMT

      ! Initialize variables for this routine
   errStat = ErrID_None
   errMsg  = ""


      ! Initialize the NWTC Subroutine Library
   call NWTC_Init( EchoLibVer=.FALSE. )

      ! Display the module information
!   call DispNVD( BEMT_Ver )


      !............................................................................................
      ! Define parameters here
      !............................................................................................
       
   call BEMT_SetParameters( InitInp, p, errStat, errMsg )
   if (errStat >= AbortErrLev) return
   p%DT = interval
      !............................................................................................
      ! Define states here
      !............................................................................................
   
   
      ! initialize the constraint states
   call BEMT_InitContraintStates( z, p, errStat2, errMsg2 )     ! initialize the constraint states
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
      if (errStat >= AbortErrLev) return
   
      ! initialize the continuous states in DBEMT and/or UA
      
      ! Initialize other states
   call BEMT_InitOtherStates( OtherState, p,  errStat, errMsg )    ! initialize the other states
   if (errStat >= AbortErrLev) return

   if ( p%DBEMT_Mod /= DBEMT_none ) then
      InitInp_DBEMT%DBEMT_Mod  = p%DBEMT_Mod
      InitInp_DBEMT%numBlades  = p%numBlades 
      InitInp_DBEMT%numNodes   = p%numBladeNodes
      InitInp_DBEMT%tau1_const = InitInp%tau1_const

      allocate(misc%u_DBEMT(2),stat=errStat2)
         if (errStat2 /= 0) then
            call SetErrStat(ErrID_Fatal,"Error allocating u_DBEMT",errStat,errMsg,RoutineName)
            return
         end if
      
      if (allocated(InitInp%rlocal)) then
         call MOVE_ALLOC( InitInp%rlocal, InitInp_DBEMT%rlocal )
      else   
         ! If not allocated we have a problem!  Issue an error and return
         call SetErrStat( ErrID_FATAL, " An InitInp%rlocal array has not been allocated and is required for DBEMT_Mod /= 0.", errStat, errMsg, RoutineName )
      end if
      if (errStat>=AbortErrLev) return
      
      call DBEMT_Init(InitInp_DBEMT, misc%u_DBEMT(1), p%DBEMT, x%DBEMT, OtherState%DBEMT, misc%DBEMT, interval, InitOut_DBEMT, errStat2, errMsg2)
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
         if (errStat >= AbortErrLev) return

      call MOVE_ALLOC( InitInp_DBEMT%rlocal, InitInp%rlocal )
      
      call DBEMT_CopyInput(misc%u_DBEMT(1),misc%u_DBEMT(2), MESH_NEWCOPY, ErrStat2, ErrMsg2)
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
   else
      p%DBEMT%lin_nx = 0
   end if
      
   
      ! in calcOutput, we will use the UA inputs for output calculations, so we must allocate them regardless of UA_Flag:
   allocate(misc%u_UA( p%numBladeNodes, p%numBlades, 2), stat=errStat2)
      if (errStat2 /= 0) then
         call SetErrStat(ErrID_Fatal,"Error allocating u_UA",errStat,errMsg,RoutineName)
         call cleanup()
         return
      end if
   
   if ( p%UA_Flag ) then
      call BEMT_Set_UA_InitData( InitInp, interval, Init_UA_Data, errStat2, errMsg2 )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
         if (errStat >= AbortErrLev) then
            call cleanup()
            return
         end if
      
      call UA_Init( Init_UA_Data, misc%u_UA(1,1,1), p%UA, x%UA, xd%UA, OtherState%UA, misc%y_UA, misc%UA, interval, AFInfo, p%AFIndx, InitOutData_UA, errStat2, errMsg2 )       
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
         if (errStat >= AbortErrLev) then
            call cleanup()
            return
         end if
   else
      p%UA%lin_nx = 0
   end if ! unsteady aero is used
   
   
      !............................................................................................
      ! Define initial guess for the system inputs here:
      !............................................................................................

         ! allocate all the arrays that store data in the input type:
   call BEMT_AllocInput( u, p, errStat2, errMsg2 )      
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
      if (errStat >= AbortErrLev) then
         call cleanup()
         return
      end if


   call BEMT_AllocOutput(y, p, errStat2, errMsg2) !u is sent so we can create sibling meshes
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
      if (errStat >= AbortErrLev) then
         call cleanup()
         return
      end if
   

   InitOut%Version = BEMT_Ver
   
 
   call AllocAry(misc%AxInduction, p%numBladeNodes,p%numBlades,'misc%AxInduction',  errStat2,errMsg2); call SetErrStat(errStat2,errMsg2,errStat,errMsg,RoutineName)
   call AllocAry(misc%TanInduction,p%numBladeNodes,p%numBlades,'misc%TanInduction', errStat2,errMsg2); call SetErrStat(errStat2,errMsg2,errStat,errMsg,RoutineName)
   call AllocAry(misc%Rtip,p%numBlades,'misc%Rtip', errStat2,errMsg2); call SetErrStat(errStat2,errMsg2,errStat,errMsg,RoutineName)
   call AllocAry(misc%phi,p%numBladeNodes,p%numBlades,'misc%phi', errStat2,errMsg2); call SetErrStat(errStat2,errMsg2,errStat,errMsg,RoutineName)
   call AllocAry(misc%chi,p%numBladeNodes,p%numBlades,'misc%chi', errStat2,errMsg2); call SetErrStat(errStat2,errMsg2,errStat,errMsg,RoutineName)
   call AllocAry(misc%ValidPhi,p%numBladeNodes,p%numBlades,'misc%ValidPhi', errStat2,errMsg2); call SetErrStat(errStat2,errMsg2,errStat,errMsg,RoutineName)
   
      if (errStat >= AbortErrLev) then
         call cleanup()
         return
      end if
   
      ! set initial values for states and misc vars
   call BEMT_ReInit(p,x,xd,z,OtherState,misc,ErrStat2,ErrMsg2)
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
   call Cleanup()
   
CONTAINS
   !...............................................................................................................................
   SUBROUTINE Cleanup()
   ! This subroutine cleans up local variables that may have allocatable arrays
   !...............................................................................................................................

   call UA_DestroyInitInput( Init_UA_Data, ErrStat2, ErrMsg2 )
   call UA_DestroyInitOutput( InitOutData_UA, ErrStat2, ErrMsg2 )

   END SUBROUTINE Cleanup

END SUBROUTINE BEMT_Init
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine reinitializes BEMT and UA, assuming that we will start the simulation over again, with only the inputs being different.
!! This allows us to bypass reading input files and allocating arrays.
subroutine BEMT_ReInit(p,x,xd,z,OtherState,misc,ErrStat,ErrMsg)

   type(BEMT_ParameterType),       intent(in   )  :: p           ! Parameters
   type(BEMT_ContinuousStateType), intent(inout)  :: x           ! Initial continuous states
   type(BEMT_DiscreteStateType),   intent(inout)  :: xd          ! Initial discrete states
   type(BEMT_ConstraintStateType), intent(inout)  :: z           ! Initial guess of the constraint states
   type(BEMT_OtherStateType),      intent(inout)  :: OtherState  ! Initial other states
   type(BEMT_MiscVarType),         intent(inout)  :: misc        ! Initial misc/optimization variables
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

   character(*), parameter                        :: RoutineName = 'BEMT_ReInit'

   ErrStat = ErrID_None
   ErrMsg  = ""
   
   misc%useFrozenWake = .FALSE.
   misc%FirstWarn_Skew = .true.
   misc%FirstWarn_Phi = .true.
   misc%FirstWarn_BEMoff = .true.
   misc%BEM_weight = 0.0_ReKi
   
   OtherState%DBEMT%tau1 = 0.0_ReKi !we're going to output this value, so let's initialize it
   
   if (p%UseInduction) then
      OtherState%ValidPhi = .true.
      
      if (p%DBEMT_Mod /= DBEMT_none ) then
         call DBEMT_ReInit(p%DBEMT, x%DBEMT, OtherState%DBEMT, misc%DBEMT)
      end if
   
   else
      OtherState%ValidPhi = .false.
      misc%AxInduction  = 0.0_ReKi
      misc%TanInduction = 0.0_ReKi
   end if
    
   z%phi = 0.0_ReKi
   OtherState%nodesInitialized = .false. ! z%phi hasn't been initialized properly, so make sure we compute a value for phi until we've updated them in the first call to BEMT_UpdateStates()

   
end subroutine BEMT_ReInit
!----------------------------------------------------------------------------------------------------------------------------------
subroutine BEMT_End( u, p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
! This routine is called at the end of the simulation.
!..................................................................................................................................

      TYPE(BEMT_InputType),           INTENT(INOUT)  :: u           ! System inputs
      TYPE(BEMT_ParameterType),       INTENT(INOUT)  :: p           ! Parameters
      TYPE(BEMT_ContinuousStateType), INTENT(INOUT)  :: x           ! Continuous states
      TYPE(BEMT_DiscreteStateType),   INTENT(INOUT)  :: xd          ! Discrete states
      TYPE(BEMT_ConstraintStateType), INTENT(INOUT)  :: z           ! Constraint states
      TYPE(BEMT_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other states
      TYPE(BEMT_OutputType),          INTENT(INOUT)  :: y           ! System outputs
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None



         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""


         ! Place any last minute operations or calculations here:


         ! Close files here:
      if ( p%UA_Flag ) then
         CALL UA_End(p%UA)
      end if


         ! Destroy the input data:

      CALL BEMT_DestroyInput( u, ErrStat, ErrMsg )
      

         ! Destroy the parameter data:

      CALL BEMT_DestroyParam( p, ErrStat, ErrMsg )


         ! Destroy the state data:

      CALL BEMT_DestroyContState(   x,           ErrStat, ErrMsg )
      CALL BEMT_DestroyDiscState(   xd,          ErrStat, ErrMsg )
      CALL BEMT_DestroyConstrState( z,           ErrStat, ErrMsg )
      CALL BEMT_DestroyOtherState(  OtherState,  ErrStat, ErrMsg )


         ! Destroy the output data:

      CALL BEMT_DestroyOutput( y, ErrStat, ErrMsg )

END SUBROUTINE BEMT_End


!----------------------------------------------------------------------------------------------------------------------------------
subroutine BEMT_UpdateStates( t, n, u1, u2,  p, x, xd, z, OtherState, AFInfo, m, errStat, errMsg )
! Loose coupling routine for solving for constraint states, integrating continuous states, and updating discrete states
! Continuous, constraint, discrete, and other states are updated for t + Interval
!
! NOTE:  This is a non-standard framework interface!!!!!  GJH
!..................................................................................................................................

   real(DbKi),                          intent(in   ) :: t          ! Current simulation time in seconds
   integer(IntKi),                      intent(in   ) :: n          ! Current simulation time step n = 0,1,...
   type(BEMT_InputType),                intent(in   ) :: u1,u2      ! Input at t and t+ dt 
   !real(DbKi),                         intent(in   ) :: utime      ! Times associated with u(:), in seconds
   type(BEMT_ParameterType),            intent(in   ) :: p          ! Parameters   
   type(BEMT_ContinuousStateType),      intent(inout) :: x          ! Input: Continuous states at t;
                                                                    !   Output: Continuous states at t + Interval
   type(BEMT_DiscreteStateType),        intent(inout) :: xd         ! Input: Discrete states at t;
                                                                    !   Output: Discrete states at t  + Interval
   type(BEMT_ConstraintStateType),      intent(inout) :: z          ! Input: Constraint states at t;
                                                                    !   Output: Constraint states at t + Interval
   type(BEMT_OtherStateType),           intent(inout) :: OtherState ! Input: Other states at t;
                                                                    !   Output: Other states at t + Interval
   type(BEMT_MiscVarType),              intent(inout) :: m          ! Misc/optimization variables
   type(AFI_ParameterType),             intent(in   ) :: AFInfo(:)  ! The airfoil parameter data
   integer(IntKi),                      intent(  out) :: errStat    ! Error status of the operation
   character(*),                        intent(  out) :: errMsg     ! Error message if ErrStat /= ErrID_None


      
   integer(IntKi)                                    :: i,j
   
   character(ErrMsgLen)                              :: errMsg2     ! temporary Error message if ErrStat /= ErrID_None
   integer(IntKi)                                    :: errStat2    ! temporary Error status of the operation
   character(*), parameter                           :: RoutineName = 'BEMT_UpdateStates'
   real(DbKi)                                        :: uTimes(2)
   
   ErrStat = ErrID_None
   ErrMsg = ""
   
   uTimes(1) = t
   uTimes(2) = t+p%dt
   
   !...............................................................................................................................
   ! if we haven't initialized z%phi, we want to get a better guess as to what the actual values of phi at t are:
   !...............................................................................................................................

   if (.not. OtherState%nodesInitialized) then
      call UpdatePhi( u1, p, z%phi, AFInfo, m, OtherState%ValidPhi, errStat2, errMsg2 )
      OtherState%nodesInitialized = .true.         ! otherState updated to t+dt (i.e., n+1)
   end if
   
   !...............................................................................................................................
   !  compute inputs to DBEMT at step n (also setting inductions--including DBEMT and skewed wake corrections--at time n)
   !...............................................................................................................................
   call BEMT_CalcOutput_Inductions( 1, t, .true., .true., z%phi, u1, p, x, xd, z, OtherState, AFInfo, m%axInduction, m%tanInduction, m%chi, m, errStat2, errMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (ErrStat >= AbortErrLev) return

#ifdef DEBUG_BEMT_RESIDUAL
      if (p%useInduction) call WriteDEBUGValuesToFile(t, u1, p, x, xd, z, OtherState, m, AFInfo)
#endif   
   !...............................................................................................................................
   !  compute inputs to UA at time n (also setting inductions--including DBEMT and skewed wake corrections--at time n)
   !...............................................................................................................................
   if (p%UA_Flag) then
      m%phi = z%phi
      call SetInputs_for_UA_AllNodes(u1, p, m%phi, m%axInduction, m%tanInduction, m%u_UA(:,:,1))
   end if
   
   !...............................................................................................................................
   !  update BEMT states to step n+1
   !...............................................................................................................................
   call UpdatePhi( u2, p, z%phi, AFInfo, m, OtherState%ValidPhi, errStat2, errMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (errStat >= AbortErrLev) return
  
   
   !...............................................................................................................................
   !  compute inputs to DBEMT at step n+1 (also setting inductions--WITHOUT DBEMT or skewed wake corrections--at step n+1)
   !...............................................................................................................................
   call BEMT_CalcOutput_Inductions( 2, t, .true., .false., z%phi, u2, p, x, xd, z, OtherState, AFInfo, m%axInduction, m%tanInduction, m%chi, m, errStat2, errMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (ErrStat >= AbortErrLev) return

   !...............................................................................................................................
   !  update DBEMT states to step n+1
   !...............................................................................................................................
   if (p%DBEMT_Mod /= DBEMT_none) then

      !........................
      ! update DBEMT states to t+dt
      !........................
      do j = 1,p%numBlades
         do i = 1,p%numBladeNodes
            call DBEMT_UpdateStates(i, j, t, n, m%u_DBEMT, p%DBEMT, x%DBEMT, OtherState%DBEMT, m%DBEMT, errStat2, errMsg2)
               if (ErrStat2 /= ErrID_None) then
                  call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//trim(NodeText(i,j)))
                  if (errStat >= AbortErrLev) return
               end if
         end do
      end do
      
   end if
   
   
   !...............................................................................................................................
   !  compute inputs to UA at time n+1 (also applying corrections to inductions--including DBEMT and skewed wake corrections)
   !...............................................................................................................................
   if (p%UA_Flag) then
         ! after updating DBEMT states, we can now apply the corrections we omitted on the last call to BEMT_CalcOutput_Inductions()
      if ( p%useInduction .and. .not. m%UseFrozenWake) then
            !............................................
            ! apply DBEMT correction to axInduction and tanInduction:
            !............................................
            if (p%DBEMT_Mod /= DBEMT_none) then
               call calculate_Inductions_from_DBEMT_AllNodes(2, uTimes(2), u2, p, x, OtherState, m, m%axInduction, m%tanInduction)
            end if
         
            call ApplySkewedWakeCorrection_AllNodes(p, u2, m, m%axInduction, m%chi)

            !............................................
            ! If TSR is too low, (start to) turn off induction
            !............................................
            call check_turnOffBEMT(p, u2, m%BEM_weight, m%axInduction, m%tanInduction, m%FirstWarn_BEMoff)
            
      end if
   
      m%phi = z%phi
      call SetInputs_for_UA_AllNodes(u2, p, m%phi, m%axInduction, m%tanInduction, m%u_UA(:,:,2))
   
      !...............................................................................................................................
      !  compute UA states at t+dt
      !...............................................................................................................................
      do j = 1,p%numBlades
         do i = 1,p%numBladeNodes

               ! COMPUTE: x%UA and/or xd%UA, OtherState%UA
            call UA_UpdateStates( i, j, t, n, m%u_UA(i,j,:), uTimes, p%UA, x%UA, xd%UA, OtherState%UA, AFInfo(p%AFIndx(i,j)), m%UA, errStat2, errMsg2 )
               if (ErrStat2 /= ErrID_None) then
                  call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//trim(NodeText(i,j)))
                  if (errStat >= AbortErrLev) return
               end if

         end do
      end do


   end if ! is UA used?
    
end subroutine BEMT_UpdateStates
!..................................................................................................................................
subroutine SetInputs_For_DBEMT(u_DBEMT, u, p, axInduction, tanInduction, Rtip)
   ! note that this subroutine inherits all data from BEMT_UpdateStates
   type(BEMT_InputType),                intent(in   ) :: u         ! BEMT Input
   type(BEMT_ParameterType),            intent(in   ) :: p         ! BEMT parameters
   type(DBEMT_InputType),               intent(inout) :: u_DBEMT   ! DBEMT Input
   real(ReKi),                          intent(in   ) :: axInduction(:,:)
   real(ReKi),                          intent(in   ) :: tanInduction(:,:)
   real(ReKi),                          intent(in   ) :: Rtip(:)

   integer                                            :: i, j
   
   
      ! Locate the maximum rlocal value for all blades.
   u_DBEMT%R_disk   = Rtip(1)
   do j = 2,p%numBlades
      u_DBEMT%R_disk   = max( u_DBEMT%R_disk  , Rtip(j) )
   end do

   if (p%DBEMT_Mod == DBEMT_tauVaries ) then
            
         ! We need to generate a disk-averaged axial induction for this timestep
      u_DBEMT%AxInd_disk = 0.0_ReKi
      do j = 1,p%numBlades
         do i = 1,p%numBladeNodes
            u_DBEMT%AxInd_disk = u_DBEMT%AxInd_disk + axInduction(i,j)
            
            u_DBEMT%element(i,j)%spanRatio  =   u%rlocal(i,j)/u_DBEMT%R_disk
         end do
      end do
      u_DBEMT%AxInd_disk = u_DBEMT%AxInd_disk / (p%numBladeNodes*p%numBlades)
         
      u_DBEMT%Un_disk    = u%Un_disk
   end if
   
   
   do j = 1,p%numBlades
      do i = 1,p%numBladeNodes
         u_DBEMT%element(i,j)%vind_s(1)  =  -axInduction( i,j)*u%Vx(i,j)  ! Eq. 38
         u_DBEMT%element(i,j)%vind_s(2)  =   tanInduction(i,j)*u%Vy(i,j)  ! Eq. 38
      end do
   end do
   
   if( allocated(u%Vx_elast_dot)) then ! only for DBEMT_Mod=DBEMT_cont_tauConst
      do j = 1,p%numBlades
         do i = 1,p%numBladeNodes
            u_DBEMT%element(i,j)%vind_s_dot(1)  =   axInduction( i,j)*u%Vx_elast_dot(i,j) - u%omega_z(i,j)*tanInduction(i,j)*u%Vy(i,j) ! Eq. 41
            u_DBEMT%element(i,j)%vind_s_dot(2)  =  -tanInduction(i,j)*u%Vy_elast_dot(i,j) - u%omega_z(i,j)*axInduction( i,j)*u%Vx(i,j) ! Eq. 41
         end do
      end do
   end if
      

end subroutine SetInputs_For_DBEMT
!..................................................................................................................................
subroutine UpdatePhi( u, p, phi, AFInfo, m, ValidPhi, errStat, errMsg )
!..................................................................................................................................

   type(BEMT_InputType),                intent(in   ) :: u              ! Input at t
   type(BEMT_ParameterType),            intent(in   ) :: p              ! Parameters
   real(ReKi),                          intent(inout) :: phi(:,:) 
   type(BEMT_MiscVarType),              intent(inout) :: m              ! Misc/optimization variables
   logical,                             intent(inout) :: ValidPhi(:,:)  ! if this is a valid BEM solution of phi
   type(AFI_ParameterType),             intent(in   ) :: AFInfo(:)      ! The airfoil parameter data
   integer(IntKi),                      intent(  out) :: errStat        ! Error status of the operation
   character(*),                        intent(  out) :: errMsg         ! Error message if ErrStat /= ErrID_None

      
   integer(IntKi)                                     :: i,j
   character(ErrMsgLen)                               :: errMsg2     ! temporary Error message if ErrStat /= ErrID_None
   integer(IntKi)                                     :: errStat2    ! temporary Error status of the operation
   character(*), parameter                            :: RoutineName = 'UpdatePhi'
   
   
   ErrStat = ErrID_None
   ErrMsg = ""
         
   !...............................................................................................................................
   ! take the current guess of z%phi and update it using inputs
   ! [called if we haven't initialized z%phi, we want to get a better guess as to what the actual values of phi at t are
   !  *or* if we are updating the BEMT states]
   !...............................................................................................................................

      if (p%useInduction) then
         
         do j = 1,p%numBlades ! Loop through all blades
            do i = 1,p%numBladeNodes ! Loop through the blade nodes / elements
               
               call BEMT_UnCoupledSolve(p, u, i, j, phi(i,j), AFInfo(p%AFIndx(i,j)), &
                                        ValidPhi(i,j), m%FirstWarn_Phi,errStat2, errMsg2)

                     if (errStat2 /= ErrID_None) then
                        call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//trim(NodeText(i,j)))
                        if (errStat >= AbortErrLev) return 
                     end if
            end do
         end do
      else
            ! We'll simply compute a geometrical phi based on both induction factors being 0.0
         do j = 1,p%numBlades ! Loop through all blades
            do i = 1,p%numBladeNodes ! Loop through the blade nodes / elements
               phi(i,j) = ComputePhiWithInduction(u%Vx(i,j), u%Vy(i,j),  0.0_ReKi, 0.0_ReKi)
            end do
         end do
         
      end if

      return
   
end subroutine UpdatePhi
!..................................................................................................................................
subroutine GetRTip( u, p, RTip )
!..................................................................................................................................

   type(BEMT_InputType),                intent(in   ) :: u              ! Input at t
   type(BEMT_ParameterType),            intent(in   ) :: p              ! Parameters
   real(ReKi),                          intent(  out) :: Rtip(:)

   integer(IntKi)                                    :: i,j
         
   do j = 1,p%numBlades

      ! Locate the maximum rlocal value for this time step and this blade.  This is passed to the solve as Rtip.
      Rtip(j) = u%rlocal(1,j)
      do i = 2,p%numBladeNodes
         Rtip(j) = max( Rtip(j), u%rlocal(i,j) )
      end do

   end do
      
   return

end subroutine GetRTip
!..................................................................................................................................
subroutine calculate_Inductions_from_BEMT(p,phi,u,OtherState,AFInfo,axInduction,tanInduction, ErrStat,ErrMsg)

   type(BEMT_ParameterType),        intent(in   ) :: p                  !< Parameters
   real(ReKi),                      intent(in   ) :: phi(:,:)           !< phi
   type(BEMT_InputType),            intent(in   ) :: u                  !< Inputs at t
   type(BEMT_OtherStateType),       intent(in   ) :: OtherState         !< Other/logical states at t
   type(AFI_ParameterType),         intent(in   ) :: AFInfo(:)          !< The airfoil parameter data
   real(ReKi),                      intent(inout) :: axInduction(:,:)   !< axial induction
   real(ReKi),                      intent(inout) :: tanInduction(:,:)  !< tangential induction
   integer(IntKi),                  intent(  out) :: errStat            !< Error status of the operation
   character(*),                    intent(  out) :: errMsg             !< Error message if ErrStat /= ErrID_None

   integer(IntKi)                                 :: i                  !< blade node counter
   integer(IntKi)                                 :: j                  !< blade counter

   real(ReKi)                                     :: fzero              !< residual from induction equation (not used here)
   logical                                        :: IsValidSolution    !< this indicates BEMT found a geometric solution because a BEMT solution could not be found
   integer(IntKi)                                 :: errStat2           !< Error status of the operation
   character(ErrMsgLen)                           :: errMsg2            !< Error message if ErrStat /= ErrID_None
   character(*), parameter                        :: RoutineName = 'calculate_Inductions_from_BEMT'
   
   ErrStat = ErrID_None
   ErrMsg = ""
   
   
   do j = 1,p%numBlades ! Loop through all blades
      do i = 1,p%numBladeNodes ! Loop through the blade nodes / elements

         if (OtherState%ValidPhi(i,j)) then
      
            ! Need to get the induction factors for these conditions without skewed wake correction and without UA
            ! COMPUTE: axInduction, tanInduction  
            fzero = BEMTU_InductionWithResidual(p, u, i, j, phi(i,j), AFInfo(p%AFIndx(i,j)), IsValidSolution, ErrStat2, ErrMsg2, a=axInduction(i,j), ap=tanInduction(i,j))
         
               if (ErrStat2 /= ErrID_None) then
                  call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//trim(NodeText(i,j)))
                  if (errStat >= AbortErrLev) return
               end if
               
            ! modify inductions based on max/min allowed values (note that we do this before calling DBEMT so that its disk-averaged induction input isn't dominated by very large values here
            if (.not. IsValidSolution) then
               axInduction(i,j) = 0.0_ReKi
               tanInduction(i,j) = 0.0_ReKi
            else
               call limitInductionFactors(axInduction(i,j), tanInduction(i,j))
            end if
      
         else
      
            axInduction(i,j) = 0.0_ReKi
            tanInduction(i,j) = 0.0_ReKi
      
         end if

      end do
   end do
   
end subroutine calculate_Inductions_from_BEMT
!..................................................................................................................................
subroutine calculate_Inductions_from_DBEMT(i, j, Vx, Vy, t, p, u, x, OtherState, m, axInduction, tanInduction)

   integer(IntKi),                  intent(in   ) :: i               !< blade node counter
   integer(IntKi),                  intent(in   ) :: j               !< blade counter
   real(ReKi),                      intent(in   ) :: Vx              !< Velocity
   real(ReKi),                      intent(in   ) :: Vy              !< Velocity
   real(DbKi),                      intent(in   ) :: t               !< Current simulation time in seconds
   type(DBEMT_ParameterType),       intent(in   ) :: p               !< Parameters
   type(DBEMT_InputType),           intent(in   ) :: u               !< Inputs for DBEMT module (note that there aren't any allocatable arrays or pointers); CalcOutput needs only vind_s to be set, and only at initialization
   type(DBEMT_ContinuousStateType), intent(in   ) :: x               !< Continuous states at t
   type(DBEMT_MiscVarType),         intent(inout) :: m               !< Initial misc/optimization variables
   type(DBEMT_OtherStateType),      intent(in   ) :: OtherState      !< Other/logical states at t
   real(ReKi),                      intent(inout) :: axInduction     !< axial induction
   real(ReKi),                      intent(inout) :: tanInduction    !< tangential induction

   real(ReKi)                      :: dbemt_vind(2)
   
   ! DBEMT doesn't set ErrStat/ErrMsg to anything but ErrStat=ErrID_None so these are local variables
   integer(IntKi)                  :: errStat              !< Error status of the operation
   character(ErrMsgLen)            :: errMsg               !< Error message if ErrStat /= ErrID_None

      ! Note that the outputs of DBEMT are the state variables x%vind, except on the first call when they uninitialized and use the inputs in place of the states

   call DBEMT_CalcOutput( i, j, t, u, dbemt_vind, p, x, OtherState, m, errStat, errMsg )
      
   if ( EqualRealnos(Vx, 0.0_ReKi) ) then
      axInduction   = 0.0_ReKi
   else
      axInduction   = -dbemt_vind(1)/Vx
   end if

   if ( EqualRealnos(Vy, 0.0_ReKi) ) then
      tanInduction  = 0.0_ReKi
   else
      tanInduction  = dbemt_vind(2)/Vy
   end if

   call limitInductionFactors(axInduction, tanInduction)

end subroutine calculate_Inductions_from_DBEMT
!----------------------------------------------------------------------------------------------------------------------------------
subroutine check_turnOffBEMT(p, u, Weight, axInduction, tanInduction, FirstWarn)

   type(BEMT_ParameterType),        intent(in   ) :: p                  !< Parameters
   type(BEMT_InputType),            intent(in   ) :: u                  !< Inputs at t
   real(ReKi),                      intent(  out) :: Weight             !< scaling value from BEM-to-non_BEM solution (blends inductions from BEMT and non-BEMT[zeros])
   real(ReKi),                      intent(inout) :: axInduction(:,:)   !< axial induction
   real(ReKi),                      intent(inout) :: tanInduction(:,:)  !< tangential induction
   logical,                         intent(inout) :: FirstWarn          !< Whether BEM had a warning about being shut off

   integer(IntKi)                                 :: i                  !< blade node counter
   integer(IntKi)                                 :: j                  !< blade counter

   if( u%TSR < BEMT_upperBoundTSR ) then
   
      Weight = BlendCosine( u%TSR, BEMT_lowerBoundTSR, BEMT_upperBoundTSR )
      
      if (FirstWarn) then
         if (Weight < 1.0_ReKi) then
               call WrScr( 'The BEM solution is being turned off due to low TSR.  (TSR = ' &
                  //TRIM(Num2Lstr(u%TSR))//'). This warning will not be repeated though the condition may persist. (See GeomPhi output channel.)' )
            FirstWarn = .false.
         end if
      end if
   
      do j = 1,p%numBlades ! Loop through all blades
         do i = 1,p%numBladeNodes ! Loop through the blade nodes / elements
            ! weighted average with BEM turned off (if TSR falls in certain range):
            axInduction( i,j) = axInduction( i,j)*Weight
            tanInduction(i,j) = tanInduction(i,j)*Weight
         end do
      end do
      
   else
      Weight = 1.0_ReKi
   end if

end subroutine check_turnOffBEMT
!----------------------------------------------------------------------------------------------------------------------------------
subroutine BEMT_CalcOutput( t, u, p, x, xd, z, OtherState, AFInfo, y, m, errStat, errMsg )
! Routine for computing outputs, used in both loose and tight coupling.
!..................................................................................................................................

   
   real(DbKi),                     intent(in   )  :: t           ! Current simulation time in seconds
   type(BEMT_InputType),           intent(in   )  :: u           ! Inputs at Time t
   type(BEMT_ParameterType),       intent(in   )  :: p           ! Parameters
   type(BEMT_ContinuousStateType), intent(in   )  :: x           ! Continuous states at t
   type(BEMT_DiscreteStateType),   intent(in   )  :: xd          ! Discrete states at t
   type(BEMT_ConstraintStateType), intent(in   )  :: z           ! Constraint states at t
   type(BEMT_OtherStateType),      intent(in   )  :: OtherState  ! Other states at t
   type(BEMT_MiscVarType),         intent(inout)  :: m           ! Misc/optimization variables
   type(AFI_ParameterType),        intent(in   )  :: AFInfo(:)   ! The airfoil parameter data
   type(BEMT_OutputType),          intent(inout)  :: y           ! Outputs computed at t (Input only so that mesh con-
                                                                 !   nectivity information does not have to be recalculated)
   integer(IntKi),                 intent(  out)  :: errStat     ! Error status of the operation
   character(*),                   intent(  out)  :: errMsg      ! Error message if ErrStat /= ErrID_None


      ! Local variables:

   integer(IntKi)                                 :: i                                               ! Generic index
   integer(IntKi)                                 :: j                                               ! Loops through nodes / elements
   integer(IntKi), parameter                      :: InputIndex=1      ! we will always use values at t in this routine
   
   character(ErrMsgLen)                           :: errMsg2     ! temporary Error message if ErrStat /= ErrID_None
   integer(IntKi)                                 :: errStat2    ! temporary Error status of the operation
   character(*), parameter                        :: RoutineName = 'BEMT_CalcOutput'
   
   type(AFI_OutputType)                           :: AFI_interp

         ! Initialize some output values
   errStat = ErrID_None
   errMsg  = ""

!!#ifdef DEBUG_BEMT_RESIDUAL
!!   call WriteDEBUGValuesToFile(t, u, p, x, xd, z, OtherState, m, AFInfo)
!!#endif

   y%phi = z%phi ! set this before possibly calling UpdatePhi() because phi is intent(inout) in the solve
   m%ValidPhi = OtherState%ValidPhi ! set this so that we don't overwrite OtherSTate%ValidPhi
   
   !...............................................................................................................................
   ! if we haven't initialized z%phi, we want to get a better guess as to what the actual values of phi are:
   !...............................................................................................................................
   if (.not. OtherState%nodesInitialized) then
      call UpdatePhi( u, p, y%phi, AFInfo, m, m%ValidPhi, errStat2, errMsg2 )
   end if
   
   !............................................
   ! calculate inductions using BEMT, applying the DBEMT, and/or skewed wake corrections as applicable:
   !............................................
   call BEMT_CalcOutput_Inductions( InputIndex, t, .false., .true., y%phi, u, p, x, xd, z, OtherState, AFInfo, y%axInduction, y%tanInduction, y%chi, m, errStat, errMsg )
   
   !............................................
   ! update phi if necessary (consistent with inductions) and calculate inputs to UA (EVEN if UA isn't used, because we use the inputs later):
   !............................................
   call SetInputs_for_UA_AllNodes(u, p, y%phi, y%axInduction, y%tanInduction, m%u_UA(:,:,InputIndex))
   
   !............................................
   ! unsteady aero and related outputs (cl, cd, cm, AoA, Vrel, Re)
   !............................................
   do j = 1,p%numBlades ! Loop through all blades
      do i = 1,p%numBladeNodes ! Loop through the blade nodes / elements

            ! set output values:
         y%AOA( i,j) = m%u_UA(i,j,InputIndex)%alpha
         y%Vrel(i,j) = m%u_UA(i,j,InputIndex)%U
         y%Re(  i,j) = m%u_UA(i,j,InputIndex)%Re
         
      enddo             ! I - Blade nodes / elements
   enddo          ! J - All blades
   
      ! Now depending on the option for UA get the airfoil coefs, Cl, Cd, Cm for unsteady or steady implementation
   if (p%UA_Flag ) then
   
      do j = 1,p%numBlades ! Loop through all blades
         do i = 1,p%numBladeNodes ! Loop through the blade nodes / elements

            call UA_CalcOutput(i, j, t, m%u_UA(i,j,InputIndex), p%UA, x%UA, xd%UA, OtherState%UA, AFInfo(p%AFindx(i,j)), m%y_UA, m%UA, errStat2, errMsg2 )
               if (ErrStat2 /= ErrID_None) then
                  call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//trim(NodeText(i,j)))
                  if (errStat >= AbortErrLev) return
               end if
               
            y%Cl(i,j) = m%y_UA%Cl
            y%Cd(i,j) = m%y_UA%Cd
            y%Cm(i,j) = m%y_UA%Cm
            y%Cpmin(i,j) = 0.0_ReKi !bjj: this isn't set anywhere... ???? 
         enddo             ! I - Blade nodes / elements
      enddo          ! J - All blades
   
      ! if ( mod(REAL(t,ReKi),.1) < p%dt) then
         call UA_WriteOutputToFile(t, p%UA, m%y_UA)
      ! end if
      
   else
            ! compute steady Airfoil Coefs
      do j = 1,p%numBlades ! Loop through all blades
         do i = 1,p%numBladeNodes ! Loop through the blade nodes / elements
         
            call AFI_ComputeAirfoilCoefs( y%AOA(i,j), y%Re(i,j), u%UserProp(i,j),  AFInfo(p%AFindx(i,j)), AFI_interp, errStat2, errMsg2 )
               if (ErrStat2 /= ErrID_None) then
                  call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//trim(NodeText(i,j)))
                  if (errStat >= AbortErrLev) return
               end if
            y%Cl(i,j) = AFI_interp%Cl
            y%Cd(i,j) = AFI_interp%Cd
            y%Cm(i,j) = AFI_interp%Cm
            y%Cpmin(i,j) = AFI_interp%Cpmin
         
         enddo             ! I - Blade nodes / elements
      enddo          ! J - All blades
      
   end if


   !............................................
   ! Compute Cx, Cy given Cl, Cd and phi
   !............................................
   do j = 1,p%numBlades ! Loop through all blades
      do i = 1,p%numBladeNodes ! Loop through the blade nodes / elements
            ! NOTE: For these calculations we force the useAIDrag and useTIDrag flags to .TRUE.
         call Transform_ClCd_to_CxCy( y%phi(i,j), .TRUE., .TRUE., y%Cl(i,j), y%Cd(i,j),y%Cx(i,j), y%Cy(i,j) )
            
      enddo             ! I - Blade nodes / elements
   enddo          ! J - All blades

   return

end subroutine BEMT_CalcOutput

!----------------------------------------------------------------------------------------------------------------------------------
!> Routine used in linearization to get the states initialized properly at t=0 (before perturbing things)
subroutine BEMT_InitStates(t, u, p, x, xd, z, OtherState, m, AFInfo, ErrStat, ErrMsg )
   REAL(DbKi),                     intent(in   )  :: t           ! current simulation time
   type(BEMT_InputType),           intent(in   )  :: u           ! Inputs at Time t
   type(BEMT_ParameterType),       intent(in   )  :: p           ! Parameters
   type(BEMT_ContinuousStateType), intent(inout)  :: x           ! Continuous states at t
   type(BEMT_DiscreteStateType),   intent(in   )  :: xd          ! Discrete states at t
   type(BEMT_ConstraintStateType), intent(in   )  :: z           ! Constraint states at t
   type(BEMT_OtherStateType),      intent(inout)  :: OtherState  ! Other states at t
   type(BEMT_MiscVarType),         intent(inout)  :: m           ! Misc/optimization variables
   type(AFI_ParameterType),        intent(in   )  :: AFInfo(:)   ! The airfoil parameter data
   integer(IntKi),                 intent(  out)  :: errStat     ! Error status of the operation
   character(*),                   intent(  out)  :: errMsg      ! Error message if ErrStat /= ErrID_None

   INTEGER(IntKi),                 parameter      :: InputIndex  = 1
   LOGICAL,                        parameter      :: CalculateDBEMTInputs = .true.
   LOGICAL,                        parameter      :: ApplyCorrections = .true.

   character(*), parameter                        :: RoutineName = 'BEMT_InitStates'
   
   
   ErrStat = ErrID_None
   ErrMsg = ""
   
   if (OtherState%nodesInitialized) return
   
   m%phi = z%phi
   call BEMT_CalcOutput_Inductions( InputIndex, t, CalculateDBEMTInputs, ApplyCorrections, m%phi, u, p, x, xd, z, OtherState, AFInfo, m%axInduction, m%tanInduction, m%chi, m, errStat, errMsg )

   if (p%DBEMT_Mod /= DBEMT_none) then
      call DBEMT_InitStates_AllNodes( m%u_DBEMT(InputIndex), p%DBEMT, x%DBEMT, OtherState%DBEMT )
   end if
   
   if (p%UA_Flag) then
      call SetInputs_for_UA_AllNodes(u, p, m%phi, m%axInduction, m%tanInduction, m%u_UA(:,:,InputIndex))
      
      call UA_InitStates_AllNodes( m%u_UA(:,:,InputIndex), p%UA, x%UA, OtherState%UA, AFInfo, p%AFIndx )
   end if ! is UA used?
   
   
end subroutine BEMT_InitStates
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine for computing inductions outputs, used in both loose and tight coupling.
subroutine BEMT_CalcOutput_Inductions( InputIndex, t, CalculateDBEMTInputs, ApplyCorrections, phi, u, p, x, xd, z, OtherState, AFInfo, axInduction, tanInduction, chi, m, errStat, errMsg )
!..................................................................................................................................

   REAL(DbKi),                     intent(in   )  :: t           ! current simulation time
   INTEGER(IntKi),                 intent(in   )  :: InputIndex  ! index into m%u_DEBMT (1 or 2)
   LOGICAL,                        intent(in   )  :: CalculateDBEMTInputs ! if called from UpdateStates, we'd want the DBEMT inputs calculated. Otherwise it saves time not to calculate them
   LOGICAL,                        intent(in   )  :: ApplyCorrections     ! if called from UpdateStates and we haven't updated DBEMT states, we'd want to avoid applying DBEMT and the skewed wake corrections
   REAL(ReKi),                     intent(in   )  :: phi(:,:)   
   type(BEMT_InputType),           intent(in   )  :: u           ! Inputs at Time t
   type(BEMT_ParameterType),       intent(in   )  :: p           ! Parameters
   type(BEMT_ContinuousStateType), intent(in   )  :: x           ! Continuous states at t
   type(BEMT_DiscreteStateType),   intent(in   )  :: xd          ! Discrete states at t
   type(BEMT_ConstraintStateType), intent(in   )  :: z           ! Constraint states at t
   type(BEMT_OtherStateType),      intent(in   )  :: OtherState  ! Other states at t
   type(BEMT_MiscVarType),         intent(inout)  :: m           ! Misc/optimization variables
   type(AFI_ParameterType),        intent(in   )  :: AFInfo(:)   ! The airfoil parameter data
   REAL(ReKi),                     intent(inout)  :: axInduction(:,:)
   REAL(ReKi),                     intent(inout)  :: tanInduction(:,:)
   REAL(ReKi),                     intent(inout)  :: chi(:,:)    ! value used in skewed wake correction
   integer(IntKi),                 intent(  out)  :: errStat     ! Error status of the operation
   character(*),                   intent(  out)  :: errMsg      ! Error message if ErrStat /= ErrID_None


      ! Local variables:

   integer(IntKi)                                 :: i                                               ! Generic index
   integer(IntKi)                                 :: j                                               ! Loops through nodes / elements
   
   character(ErrMsgLen)                           :: errMsg2     ! temporary Error message if ErrStat /= ErrID_None
   integer(IntKi)                                 :: errStat2    ! temporary Error status of the operation
   character(*), parameter                        :: RoutineName = 'BEMT_CalcOutput_Inductions'
    
   
         ! Initialize some output values
   errStat = ErrID_None
   errMsg  = ""

   
   chi = u%chi0    ! with no induction, chi = chi0 at all nodes (overwritten in ApplySkewedWakeCorrection)
   
   if ( p%useInduction ) then
   
      if (m%UseFrozenWake) then ! note that this is used only when DBEMT_Mod==DBEMT_none
      
         do j = 1,p%numBlades ! Loop through all blades
            do i = 1,p%numBladeNodes ! Loop through the blade nodes / elements
                  ! note that we've already checked if these velocities are zero before calling this routine in linearization
            !BJJ: THIS ISN'T NECESSARIALLY TRUE ANYMORE
               if (EqualRealNos(u%Vx(i,j),0.0_ReKi)) then
                  axInduction(i,j)  =  0.0_ReKi ! does this violate the frozen wake requirements????
               else
                  axInduction(i,j)  = -m%AxInd_op(i,j) / u%Vx(i,j)
               end if
               if (EqualRealNos(u%Vy(i,j),0.0_ReKi)) then
                  tanInduction(i,j) =  0.0_ReKi ! does this violate the frozen wake requirements????
               else
                  tanInduction(i,j) =  m%TnInd_op(i,j) / u%Vy(i,j)
               end if               
            enddo             ! I - Blade nodes / elements
         enddo          ! J - All blades
         
      else
      
         ! Locate the maximum rlocal value for this time step and this blade.  This is passed to the solve as Rtip.
         call GetRTip( u, p, m%Rtip )

         !............................................
         ! get BEMT inductions (axInduction and tanInduction):
         !............................................
         call calculate_Inductions_from_BEMT(p, phi, u, OtherState, AFInfo, axInduction, tanInduction, ErrStat2, ErrMsg2)
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
               if (errStat >= AbortErrLev) return

         !............................................
         ! apply DBEMT correction to axInduction and tanInduction:
         !............................................
         if (p%DBEMT_Mod /= DBEMT_none) then
            ! If we are using DBEMT, then we will obtain the time-filtered versions of axInduction(i,j), tanInduction(i,j)
         
            ! Note that the outputs of DBEMT are the state variables x%vind, so we don't NEED to set the inputs except on initialization step (when we output the inputs instead of the states)
            ! If calling from UpdateStates, we'd want to have a copy of the DBEMT inputs for updating its states later
            if (CalculateDBEMTInputs .or. .not. OtherState%nodesInitialized) then
               call SetInputs_For_DBEMT(m%u_DBEMT(InputIndex), u, p, axInduction, tanInduction, m%Rtip)  ! DBEMT inputs at InputIndex time (1=n; 2=n+1)
            end if

            if (ApplyCorrections) then ! avoid updating inductions when called before updating DBEMT states
               call calculate_Inductions_from_DBEMT_AllNodes(InputIndex, t, u, p, x, OtherState, m, axInduction, tanInduction)
            end if
         end if ! DBEMT

         if (ApplyCorrections) then ! avoid updating inductions when called before updating DBEMT states
            !............................................
            ! Apply skewed wake correction to the axial induction (axInduction)
            !............................................
            call ApplySkewedWakeCorrection_AllNodes(p, u, m, axInduction, chi)
            
            !............................................
            ! If TSR is too low, (start to) turn off induction
            !............................................
            call check_turnOffBEMT(p, u, m%BEM_weight, axInduction, tanInduction, m%FirstWarn_BEMoff)
         end if

      end if ! UseFrozenWake

   else ! no induction:
      axInduction  = 0.0_ReKi ! all nodes
      tanInduction = 0.0_ReKi ! all nodes
   end if ! UseInduction

   return
end subroutine BEMT_CalcOutput_Inductions
!----------------------------------------------------------------------------------------------------------------------------------
subroutine calculate_Inductions_from_DBEMT_AllNodes(InputIndex, t, u, p, x, OtherState, m, axInduction, tanInduction)
   INTEGER(IntKi),                 intent(in   )  :: InputIndex  ! index into m%u_DEBMT (1 or 2)
   REAL(DbKi),                     intent(in   )  :: t           ! current simulation time
   type(BEMT_InputType),           intent(in   )  :: u           ! Inputs at Time t
   type(BEMT_ParameterType),       intent(in   )  :: p           ! Parameters
   type(BEMT_ContinuousStateType), intent(in   )  :: x           ! Continuous states at t
   type(BEMT_OtherStateType),      intent(in   )  :: OtherState  ! Other states at t
   type(BEMT_MiscVarType),         intent(inout)  :: m           ! Misc/optimization variables
   REAL(ReKi),                     intent(inout)  :: axInduction(:,:)
   REAL(ReKi),                     intent(inout)  :: tanInduction(:,:)


      ! Local variables:

   integer(IntKi)                                 :: i                                               ! Generic index
   integer(IntKi)                                 :: j                                               ! Loops through nodes / elements
   
   
   do j = 1,p%numBlades ! Loop through all blades
      do i = 1,p%numBladeNodes ! Loop through the blade nodes / elements
         if ( .not. p%FixedInductions(i,j) ) then
            call calculate_Inductions_from_DBEMT(i, j, u%Vx(i,j), u%Vy(i,j), t, p%DBEMT, m%u_DBEMT(InputIndex), x%DBEMT, OtherState%DBEMT, m%DBEMT, axInduction(i,j), tanInduction(i,j))
         end if ! .not. p%FixedInductions (special case for tip and/or hub loss)
      enddo             ! I - Blade nodes / elements
   enddo          ! J - All blades


end subroutine calculate_Inductions_from_DBEMT_AllNodes
!----------------------------------------------------------------------------------------------------------------------------------
subroutine ApplySkewedWakeCorrection_AllNodes(p, u, m, axInduction, chi)
   type(BEMT_InputType),           intent(in   )  :: u           ! Inputs at Time t
   type(BEMT_ParameterType),       intent(in   )  :: p           ! Parameters
   type(BEMT_MiscVarType),         intent(inout)  :: m           ! Misc/optimization variables
   REAL(ReKi),                     intent(inout)  :: axInduction(:,:)
   REAL(ReKi),                     intent(inout)  :: chi(:,:)    ! value used in skewed wake correction

   integer(IntKi)                                 :: i                                               ! Generic index
   integer(IntKi)                                 :: j                                               ! Loops through nodes / elements
   
   !............................................
   ! Apply skewed wake correction to the axial induction (y%axInduction)
   !............................................
   if ( p%skewWakeMod == SkewMod_PittPeters ) then
      do j = 1,p%numBlades ! Loop through all blades
         do i = 1,p%numBladeNodes ! Loop through the blade nodes / elements
            if ( .not. p%FixedInductions(i,j) ) then
               call ApplySkewedWakeCorrection( p%yawCorrFactor, u%psi(j), u%chi0, u%rlocal(i,j)/m%Rtip(j), axInduction(i,j), chi(i,j), m%FirstWarn_Skew )
            end if ! .not. p%FixedInductions (special case for tip and/or hub loss)
         enddo    ! I - Blade nodes / elements
      enddo       ! J - All blades
   end if

end subroutine ApplySkewedWakeCorrection_AllNodes

!----------------------------------------------------------------------------------------------------------------------------------
subroutine BEMT_CalcContStateDeriv( t, u, p, x, xd, z, OtherState, m, dxdt, AFInfo, ErrStat, ErrMsg )
! Tight coupling routine for computing derivatives of continuous states
!..................................................................................................................................

   REAL(DbKi),                     INTENT(IN   )  :: t           ! Current simulation time in seconds
   TYPE(BEMT_InputType),           INTENT(IN   )  :: u           ! Inputs at t
   TYPE(BEMT_ParameterType),       INTENT(IN   )  :: p           ! Parameters
   TYPE(BEMT_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at t
   TYPE(BEMT_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at t
   TYPE(BEMT_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at t
   TYPE(BEMT_OtherStateType),      INTENT(IN   )  :: OtherState  ! Other states at t
   TYPE(BEMT_MiscVarType),         INTENT(INOUT)  :: m           ! Misc/optimization variables
   TYPE(BEMT_ContinuousStateType), INTENT(INOUT)  :: dxdt        ! Continuous state derivatives at t
   TYPE(AFI_ParameterType),        INTENT(IN   )  :: AFInfo(:)   ! The airfoil parameter data
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   ! local variables
   CHARACTER(ErrMsgLen)                           :: ErrMsg2     ! temporary Error message if ErrStat /= ErrID_None
   INTEGER(IntKi)                                 :: ErrStat2    ! temporary Error status of the operation
   CHARACTER(*), PARAMETER                        :: RoutineName = 'BEMT_CalcContStateDeriv'
   
   INTEGER(IntKi), parameter                      :: InputIndex = 1
   INTEGER(IntKi)                                 :: i, j

      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ""

   
   
   m%phi = z%phi ! set this before possibly calling InitPhi() because phi is intent(inout) in the solve
   m%ValidPhi = OtherState%ValidPhi ! set this so that we don't overwrite OtherSTate%ValidPhi
   
   !...............................................................................................................................
   ! THIS IS CALLED ONLY DURING LINEARIZATION, where we want to treat z%phi as an implied input, thus we need to recalculate phi
   ! consistent with the inputs:
   !...............................................................................................................................
   call UpdatePhi( u, p, m%phi, AFInfo, m, m%ValidPhi, errStat2, errMsg2 )
   
   !...............................................................................................................................
   !  compute inputs to DBEMT (also setting inductions needed for UA inputs--including DBEMT and skewed wake corrections)
   !...............................................................................................................................
   call BEMT_CalcOutput_Inductions( InputIndex, t, .true., .true., m%phi, u, p, x, xd, z, OtherState, AFInfo, m%axInduction, m%tanInduction, m%chi, m, errStat2, errMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (ErrStat >= AbortErrLev) return
      
   !...............................................................................................................................
   !  compute derivatives for DBEMT continuous states:
   !...............................................................................................................................
   if (p%DBEMT_Mod /= DBEMT_none) then
      if (.not. allocated(dxdt%DBEMT%element)) then
         call DBEMT_CopyContState( x%DBEMT, dxdt%DBEMT, MESH_UPDATECOPY, ErrStat2, ErrMsg2 )
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            if (ErrStat >= AbortErrLev) return
      end if
      
      do j = 1,p%numBlades ! Loop through all blades
         do i = 1,p%numBladeNodes ! Loop through the blade nodes / elements
            call DBEMT_CalcContStateDeriv( i, j, t, m%u_DBEMT(InputIndex)%element(i,j), p%DBEMT, x%DBEMT%element(i,j), OtherState%DBEMT, m%DBEMT, dxdt%DBEMT%element(i,j), ErrStat2, ErrMsg2 )
               call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
               if (ErrStat >= AbortErrLev) return
         enddo             ! I - Blade nodes / elements
      enddo          ! J - All blades
    end if

      
   
   !...............................................................................................................................
   !  compute inputs to UA
   !...............................................................................................................................
   if (p%UA_Flag) then
      call SetInputs_for_UA_AllNodes(u, p, m%phi, m%axInduction, m%tanInduction, m%u_UA(:,:,InputIndex))

   !...............................................................................................................................
   !  compute derivatives for UA continuous states:
   !...............................................................................................................................
      if (.not. allocated(dxdt%UA%element)) then
         call UA_CopyContState( x%UA, dxdt%UA, MESH_UPDATECOPY, ErrStat2, ErrMsg2 )
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            if (ErrStat >= AbortErrLev) return
      end if
      
      do j = 1,p%numBlades
         do i = 1,p%numBladeNodes
            call UA_CalcContStateDeriv( i, j, t, m%u_UA(i,j,InputIndex), p%UA, x%UA%element(i,j), OtherState%UA, AFInfo(p%AFIndx(i,j)), m%UA, dxdt%UA%element(i,j), ErrStat2, ErrMsg2 )
         end do
      end do
     
   end if
   
END SUBROUTINE BEMT_CalcContStateDeriv
!----------------------------------------------------------------------------------------------------------------------------------
subroutine BEMT_UpdateDiscState( t, n, u, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
! Tight coupling routine for updating discrete states
!..................................................................................................................................

   REAL(DbKi),                     INTENT(IN   )  :: t           ! Current simulation time in seconds
   INTEGER(IntKi),                 INTENT(IN   )  :: n           ! Current step of the simulation: t = n*Interval
   TYPE(BEMT_InputType),           INTENT(IN   )  :: u           ! Inputs at t
   TYPE(BEMT_ParameterType),       INTENT(IN   )  :: p           ! Parameters
   TYPE(BEMT_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at t
   TYPE(BEMT_DiscreteStateType),   INTENT(INOUT)  :: xd          ! Input: Discrete states at t;
                                                                 !   Output: Discrete states at t + Interval
   TYPE(BEMT_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at t
   TYPE(BEMT_OtherStateType),      INTENT(IN   )  :: OtherState  ! Other states at t
   type(BEMT_MiscVarType),         intent(inout)  :: m           ! Misc/optimization variables
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None


         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""


         ! Update discrete states here:

      ! StateData%DiscState =

END SUBROUTINE BEMT_UpdateDiscState
!----------------------------------------------------------------------------------------------------------------------------------
subroutine BEMT_CalcConstrStateResidual( Time, u, p, x, xd, z, OtherState, m, z_residual, AFInfo, ErrStat, ErrMsg )
! Tight coupling routine for solving for the residual of the constraint state equations
!..................................................................................................................................

   real(DbKi),                     intent(in   )  :: Time        ! Current simulation time in seconds
   type(BEMT_InputType),           intent(in   )  :: u           ! Inputs at Time
   type(BEMT_ParameterType),       intent(in   )  :: p           ! Parameters
   type(BEMT_ContinuousStateType), intent(in   )  :: x           ! Continuous states at Time
   type(BEMT_DiscreteStateType),   intent(in   )  :: xd          ! Discrete states at Time
   type(BEMT_ConstraintStateType), intent(in   )  :: z           ! Constraint states at Time (possibly a guess)
   type(BEMT_OtherStateType),      intent(in   )  :: OtherState  ! Other states at Time
   type(BEMT_MiscVarType),         intent(inout)  :: m           ! Misc/optimization variables
   type(BEMT_ConstraintStateType), intent(inout)  :: z_residual  ! Residual of the constraint state equations using
                                                                 !     the input values described above
   type(AFI_ParameterType),        intent(in   )  :: AFInfo(:)   ! The airfoil parameter data
   integer(IntKi),                 intent(  out)  :: ErrStat     ! Error status of the operation
   character(*),                   intent(  out)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None


      !# set epsilon
   !REAL(ReKi), PARAMETER     ::epsilon = 1e-6
   
      ! Local variables
   INTEGER    :: i,j
   character(ErrMsgLen)                           :: errMsg2     ! temporary Error message if ErrStat /= ErrID_None
   integer(IntKi)                                 :: errStat2    ! temporary Error status of the operation
   character(*), parameter                        :: RoutineName = 'BEMT_CalcConstrStateResidual'
   logical                                        :: IsValidSolution ! placeholder for flag to determine if the residual solution is invalid
  
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   
   if (p%useInduction) then 
      
      if ( m%UseFrozenWake ) then ! we are linearizing with frozen wake assumption; i.e., p%FrozenWake is true and this was called from the linearization routine
         do j = 1,p%numBlades            
            do i = 1,p%numBladeNodes
               Z_residual%phi(i,j) = sin(z%phi(i,j)) * (u%Vy(i,j) + m%TnInd_op(i,j)) - cos(z%phi(i,j)) * (u%Vx(i,j) + m%AxInd_op(i,j))
            end do
         end do
         
      else
            
         do j = 1,p%numBlades            
            do i = 1,p%numBladeNodes
         
                  ! Solve for the constraint states here:
               Z_residual%phi(i,j) = BEMTU_InductionWithResidual(p, u, i, j, z%phi(i,j), AFInfo(p%AFindx(i,j)), IsValidSolution, ErrStat2, ErrMsg2)
               call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
               if (ErrStat >= AbortErrLev) return                  
           
            end do
         end do
         
      end if ! not frozen wake
      
   else
      
      do j = 1,p%numBlades            
         do i = 1,p%numBladeNodes     
            Z_residual%Phi(i,j) = sin(z%phi(i,j)) * u%Vy(i,j) - cos(z%phi(i,j)) * u%Vx(i,j)
         end do
      end do
      

   end if
   
   
END SUBROUTINE BEMT_CalcConstrStateResidual

!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine computes the BEMT inductions that are frozen when linearizing with the FrozenWake flag.
!> It first calls BEMT_CalcOutput to compute y%tanInduction and y%axInduction at this operating point.
!SUBROUTINE computeFrozenWake( t, u, p, x, xd, z, OtherState, AFInfo, y, m, errStat, errMsg )
SUBROUTINE computeFrozenWake( u, p, y, m )
!..................................................................................................................................

   
   !real(DbKi),                     intent(in   )  :: t           ! Current simulation time in seconds
   type(BEMT_InputType),           intent(in   )  :: u           ! Inputs at Time t
   type(BEMT_ParameterType),       intent(in   )  :: p           ! Parameters
   !type(BEMT_ContinuousStateType), intent(in   )  :: x           ! Continuous states at t
   !type(BEMT_DiscreteStateType),   intent(in   )  :: xd          ! Discrete states at t
   !type(BEMT_ConstraintStateType), intent(in   )  :: z           ! Constraint states at t
   !type(BEMT_OtherStateType),      intent(in   )  :: OtherState  ! Other states at t
   type(BEMT_MiscVarType),         intent(inout)  :: m           ! Misc/optimization variables
   !type(AFI_ParameterType),        intent(in   )  :: AFInfo(:)   ! The airfoil parameter data
   type(BEMT_OutputType),          intent(inout)  :: y           ! Outputs computed at t 
   !integer(IntKi),                 intent(  out)  :: errStat     ! Error status of the operation
   !character(*),                   intent(  out)  :: errMsg      ! Error message if ErrStat /= ErrID_None


      ! local variables
   INTEGER(IntKi)                                 :: j,k  ! loop counters
   character(*), parameter                        :: RoutineName = 'computeFrozenWake'
   
      ! get a and aprime            
   !call BEMT_CalcOutput(t, u, p, x, xd, z, OtherState, AFInfo, y, m, errStat, errMsg)
      
   do k = 1,p%numBlades            
      do j = 1,p%numBladeNodes
            
         m%AxInd_op(j,k) = - u%Vx(j,k) * y%axInduction( j,k)
         m%TnInd_op(j,k) =   u%Vy(j,k) * y%tanInduction(j,k)
            
      end do
   end do
      
      
     
END SUBROUTINE computeFrozenWake
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine checks for BEMT inputs that are invalid when linearizing constraint state equations.
SUBROUTINE CheckLinearizationInput(p, u, z, m, OtherState, ErrStat, ErrMsg)

   type(BEMT_ParameterType),       intent(in   )  :: p           ! Parameters
   type(BEMT_InputType),           intent(in   )  :: u           ! Inputs at the operating point
   type(BEMT_ConstraintStateType), intent(in   )  :: z           ! Constraint states at the operating point
   type(BEMT_MiscVarType),         intent(in   )  :: m           ! Misc/optimization variables
   type(BEMT_OtherStateType),      intent(in   )  :: OtherState  ! Other state at the operating point
   integer(IntKi),                 intent(  out)  :: ErrStat     ! Error status of the operation
   character(*),                   intent(  out)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! local variables
   INTEGER(IntKi)                                 :: j,k  ! loop counters
   character(*), parameter                        :: RoutineName = 'CheckLinearizationInput'
   
   
   ErrStat = ErrID_None
   ErrMsg  = ''
   
   if (p%UseInduction) then
      
      do k = 1,p%numBlades            
         do j = 1,p%numBladeNodes
            if (.not. OtherState%ValidPhi(j,k)) then
               call SetErrStat(ErrID_Fatal,"Blade"//trim(num2lstr(k))//', node '//trim(num2lstr(j))//&
                      ": Current operating point does not contain a valid value of phi.",ErrStat,ErrMsg,RoutineName)
               return
            end if
            
         end do
      end do
         
      

      if ( m%UseFrozenWake )  then ! we are linearizing with frozen wake assumption (i.e., p%FrozenWake is true and this is called from linearization routine)
      
         do k = 1,p%numBlades            
            do j = 1,p%numBladeNodes
            
               if ( VelocityIsZero( u%Vy(j,k)+m%TnInd_op(j,k)) .and. VelocityIsZero( u%Vx(j,k)+m%AxInd_op(j,k) ) ) then
                  call SetErrStat(ErrID_Fatal,"Blade"//trim(num2lstr(k))//', node '//trim(num2lstr(j))//&
                      ": residual is undefined because u%Vy + TnInd_op = u%Vx + AxInd_op = 0.",ErrStat,ErrMsg,RoutineName)
                  return
               end if
            
            end do
         end do
      
      else
      
         do k = 1,p%numBlades            
            do j = 1,p%numBladeNodes
            
               if ( EqualRealNos(z%phi(j,k), 0.0_ReKi) ) then
                  call SetErrStat(ErrID_Fatal,"Blade"//trim(num2lstr(k))//', node '//trim(num2lstr(j))//&
                       ": residual is discontinuous or undefined because z%phi = 0.",ErrStat,ErrMsg,RoutineName)
                  return
               else if ( VelocityIsZero(u%Vy(j,k)) ) then
                  call SetErrStat(ErrID_Fatal,"Blade"//trim(num2lstr(k))//', node '//trim(num2lstr(j))//&
                       ": residual is discontinuous or undefined because u%Vy = 0.",ErrStat,ErrMsg,RoutineName)
                  return
               else if ( VelocityIsZero(u%Vx(j,k)) ) then
                  call SetErrStat(ErrID_Fatal,"Blade"//trim(num2lstr(k))//', node '//trim(num2lstr(j))//&
                       ": residual is discontinuous or undefined because u%Vx = 0.",ErrStat,ErrMsg,RoutineName)
                  return
               end if
                           
            end do
         end do
               
      end if      
      
      
            
   else ! .not. p%UseInduction:
      
      do k = 1,p%numBlades            
         do j = 1,p%numBladeNodes
            
            if ( EqualRealNos( u%Vy(j,k), 0.0_ReKi ) .and. EqualRealNos( u%Vx(j,k), 0.0_ReKi ) ) then
               call SetErrStat(ErrID_Fatal,"Blade"//trim(num2lstr(k))//', node '//trim(num2lstr(j))//&
                    ": residual is undefined because u%Vy = u%Vx = 0.",ErrStat,ErrMsg,RoutineName)
               return
            end if
            
         end do
      end do      
      

                  
   end if
                  
END SUBROUTINE CheckLinearizationInput
!----------------------------------------------------------------------------------------------------------------------------------
!> This gets the constraint-state perturbations for linearization about a given operating point. This returns two deltas (one plus 
!! and one minus) such that z_op+dz_p and z_op-dz_m are in the same solution region (or that they are in adjacent regions that have 
!! continuous solution regions [i.e, the pi/2 boundary is okay because it is continuous across the momentum/empirical regions]).
SUBROUTINE Get_phi_perturbations(p, m, z_op, dz_p, dz_m)

   type(BEMT_ParameterType),       intent(in   )  :: p           ! Parameters
   type(BEMT_MiscVarType),         intent(in   )  :: m           ! Misc/optimization variables

   REAL(ReKi), intent(in    ) :: z_op     !< value of z%phi(i,j) at the operating point
   REAL(R8Ki), intent(   out) :: dz_p     !< change in z_op in the plus direction 
   REAL(R8Ki), intent(   out) :: dz_m     !< change in z_op in the minus direction

      ! local variables
   real(R8Ki)                 :: dz         ! size of perturbation
   real(R8Ki)                 :: zp         ! z_op+dz
   real(R8Ki)                 :: zm         ! z_op-dz

   
   dz = 2.0_R8Ki*D2R_D
   
      ! we'll assume a central difference unless we are on the boundaries below  [default]
   dz_p = dz
   dz_m = dz
   
   
   if (p%UseInduction .and. .not. m%UseFrozenWake) then
   
      zp = z_op+dz
      zm = z_op-dz
   
         ! check if it goes past the pi-eps upper boundary
      if ( zp > pi_D - BEMT_epsilon2 ) then
         ! 1-sided difference
         dz_p = 0.0_R8Ki
         dz_m = dz
      
         ! next we care about the -pi/4-eps boundary:
      else if ( zm < -pi_D/4.0_R8Ki - BEMT_epsilon2) then
         ! 1-sided difference
         dz_p = dz
         dz_m = 0.0_R8Ki
      
         ! next let's check the +eps boundaries:
      else if ( z_op > 0.0_R8Ki .and. zm < BEMT_epsilon2 ) then
         ! 1-sided difference
         dz_p = dz
         dz_m = 0.0_R8Ki
      
         ! next let's check the -eps boundaries:
      else if ( z_op < 0.0_ReKi .and. zp > -BEMT_epsilon2 ) then
         ! 1-sided difference
         dz_p = 0.0_R8Ki
         dz_m = dz
      
      ! else ! we don't care about the pi/2 boundary, so let's do a central difference for everything else
      end if
   
   end if
      
   
END SUBROUTINE Get_phi_perturbations
!----------------------------------------------------------------------------------------------------------------------------------
subroutine GetSolveRegionOrdering(Vx, phiIn, test_lower, test_upper)
   real(ReKi),             intent(in   ) :: Vx
   real(ReKi),             intent(in   ) :: phiIn
   real(ReKi),             intent(  out) :: test_lower(3)
   real(ReKi),             intent(  out) :: test_upper(3)


   if (Vx > 0) then
   
      test_lower(1) = BEMT_epsilon2
      test_upper(1) = PiBy2 - BEMT_epsilon2
   
      if (phiIn < pi/4.0_ReKi  .and. phiIn > -pi/4.0_ReKi) then !bjj: added the negative for those cases where the previously calculated non-BEMT phi is in the [-pi,-pi/4] range
         test_lower(2) = -pi/4.0_ReKi
         test_upper(2) = -BEMT_epsilon2

         test_lower(3) = PiBy2 + BEMT_epsilon2
         test_upper(3) = pi - BEMT_epsilon2
      else
         test_lower(3) = -pi/4.0_ReKi
         test_upper(3) = -BEMT_epsilon2

         test_lower(2) = PiBy2 + BEMT_epsilon2
         test_upper(2) = pi - BEMT_epsilon2
      end if
      
   else
      
      test_lower(1) = -BEMT_epsilon2
      test_upper(1) = -PiBy2 + BEMT_epsilon2
   
      if (phiIn > -pi/4.0_ReKi  .and. phiIn < pi/4.0_ReKi) then !bjj: added the negative for those cases where the previously calculated non-BEMT phi is in the [-pi,-pi/4] range
         test_lower(2) = pi/4.0_ReKi
         test_upper(2) = BEMT_epsilon2

         test_lower(3) = -PiBy2 - BEMT_epsilon2
         test_upper(3) = -pi + BEMT_epsilon2
      else
         test_lower(3) = pi/4.0_ReKi
         test_upper(3) = BEMT_epsilon2

         test_lower(2) = -PiBy2 - BEMT_epsilon2
         test_upper(2) = -pi + BEMT_epsilon2
      end if

   end if
   
   
end subroutine GetSolveRegionOrdering
   
integer function FindTestRegion(p, u, iBladeNode, jBlade, phiLower, phiUpper, AFInfo, &
                        phiIn_IsValidSolution, phiIn, f_phiIn, f1, f2, errStat, errMsg) result(TestRegion)

   type(BEMT_ParameterType),intent(in  ) :: p
   type(BEMT_InputType),   intent(in   ) :: u
   integer(IntKi),         intent(in   ) :: iBladeNode         !< index for blade node
   integer(IntKi),         intent(in   ) :: jBlade             !< index for blade
   real(ReKi),             intent(inout) :: phiLower !intent "out" in case the previous solution can alter the test region bounds
   real(ReKi),             intent(inout) :: phiUpper !intent "out" in case the previous solution can alter the test region bounds
   !integer,                intent(in   ) :: numBladeNodes
   type(AFI_ParameterType),intent(in   ) :: AFInfo
   logical,                intent(in   ) :: phiIn_IsValidSolution
   real(ReKi),             intent(in   ) :: phiIn
   real(ReKi),             intent(in   ) :: f_phiIn
   real(ReKi),             intent(  out) :: f1 !< value of residual at phiLower
   real(ReKi),             intent(  out) :: f2 !< value of residual at phiUpper
   integer(IntKi),         intent(  out) :: ErrStat       ! Error status of the operation
   character(*),           intent(  out) :: ErrMsg        ! Error message if ErrStat /= ErrID_None
   
      ! Local variables  
   character(errMsgLen)                  :: errMsg2                 ! temporary Error message if ErrStat /= ErrID_None
   integer(IntKi)                        :: errStat2                ! temporary Error status of the operation
   character(*), parameter               :: RoutineName='TestRegion'
   logical                               :: IsValidSolution, IsValidSolution2 ! placeholder for flag to determine if the residual solution is invalid (we'll handle that after the brent solve) 
   

   ErrStat = ErrID_None
   ErrMsg  = ""
   
   
   f1 = BEMTU_InductionWithResidual(p, u, iBladeNode, jBlade, phiLower, AFInfo, IsValidSolution, errStat2, errMsg2)
   
   call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName ) 
   if (errStat >= AbortErrLev) return
   
   f2 = BEMTU_InductionWithResidual(p, u, iBladeNode, jBlade, phiUpper, AFInfo, IsValidSolution2, errStat2, errMsg2)
   call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName ) 
   if (errStat >= AbortErrLev) return
   
   
      ! Look for zero-crossing
   if ( EqualRealNos(f1, 0.0_ReKi) .and. EqualRealNos(f2, 0.0_ReKi) .and. IsValidSolution .and. IsValidSolution2) then
      TestRegion = 0  ! all solutions yield zero -- special case
      return
   else
      if ( abs(f1) < p%aTol ) then
         if ( abs(f2) < abs(f1) .and. IsValidSolution2 ) then
            TestRegion = 4 ! special case: upper end point is a zero (and it's smaller than the solution at the lower end point)
            return
         elseif ( IsValidSolution ) then
            TestRegion = 3 ! special case: lower end point is a zero
            return
         end if
      elseif ( abs(f2) < p%aTol .and. IsValidSolution2 ) then
            TestRegion = 4 ! special case: upper end point is a zero
            return
      end if
      
   end if
   
   if ( sign(1.0_ReKi,f1) /= sign(1.0_ReKi,f2) ) then
      TestRegion = 1
      
      if (phiIn_IsValidSolution) then
         if ( (phiLower < phiIn .and. phiIn < phiUpper ) .or. (phiUpper < phiIn .and. phiIn < phiLower ) ) then
         
            ! the previous solution was in this region
            if ( sign(1.0_ReKi,f1) /= sign(1.0_ReKi,f_phiIn) ) then     
               phiUpper = phiIn
               f2 = f_phiIn
            else
               phiLower = phiIn
               f1 = f_phiIn
            end if
            
         end if
      end if
      
      
   else
      TestRegion = 2  ! No zero
   end if
   
end function FindTestRegion
      
subroutine BEMT_UnCoupledSolve(p, u, iBladeNode, jBlade, phi, AFInfo, ValidPhi, FirstWarn, ErrStat, ErrMsg)

   use mod_root1dim
   !use fminMod
   type(BEMT_ParameterType),intent(in  ) :: p
   type(BEMT_InputType),    intent(in  ) :: u
   integer(IntKi),          intent(in  ) :: iBladeNode         !< index for blade node
   integer(IntKi),          intent(in  ) :: jBlade             !< index for blade
   real(ReKi),             intent(inout) :: phi
   TYPE(AFI_ParameterType),INTENT(IN   ) :: AFInfo
   logical,                intent(inout) :: ValidPhi
   logical,                intent(inout) :: FirstWarn
   integer(IntKi),         intent(  out) :: errStat       ! Error status of the operation
   character(*),           intent(  out) :: errMsg        ! Error message if ErrStat /= ErrID_None
   
   
   
      ! Local variables
   character(ErrMsgLen)                  :: errMsg2       ! temporary Error message if ErrStat /= ErrID_None
   integer(IntKi)                        :: errStat2      ! temporary Error status of the operation
   character(*), parameter               :: RoutineName = 'BEMT_UnCoupledSolve'
   real(ReKi), parameter                 :: MsgLimit = 0.07_ReKi  !BEMT_epsilon2*100.0_ReKi don't print a message if we're within about 4 degrees of 0 or +/- pi/2 [bjj arbitrary number]
   
   real(ReKi) :: f1, f_lower, f_upper, phiIn
   real(ReKi) :: phi_lower(3), phi_upper(3)         ! upper and lower bounds for region of phi in which we are trying to find a solution to the BEM equations
   integer    :: i, TestRegionResult
   logical    :: IsValidSolution
   
   ErrStat = ErrID_None
   ErrMsg  = ""
  
   
   if ( VelocityIsZero(u%Vx(iBladeNode,jBlade)) ) then
      phi =  0.0_ReKi
      ValidPhi = .true.
      return
   else if ( VelocityIsZero(u%Vy(iBladeNode,jBlade)) ) then
      if (u%Vx(iBladeNode,jBlade) > 0.0_ReKi) then
         phi =  PiBy2
      else
         phi = -PiBy2
      end if
      ValidPhi = .true.
      return
   end if
   
   
   
   !# ------ BEM solution method see (Ning, doi:10.1002/we.1636) ------
   
      ! See if the previous value of phi still satisfies the residual equation.
      ! (If the previous phi wasn't a valid solution to BEMT equations, skip this check and just perform the solve)
   if (ValidPhi .and. .NOT. EqualRealNos(phi, 0.0_ReKi) .and. .not. EqualRealNos(abs(phi),PiBy2) ) then  
      f1 = BEMTU_InductionWithResidual(p, u, iBladeNode, jBlade, phi, AFInfo, IsValidSolution, errStat2, errMsg2)
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName ) 
      if (errStat >= AbortErrLev) return

      if ( abs(f1) < p%aTol .and. IsValidSolution ) then
         !phiStar =  phiIn
         return
      end if
   else
      IsValidSolution = .false.
   end if
   
   
      
   ! 
   ValidPhi = .false. ! initialize to false while we try to find a new valid solution
   
   call GetSolveRegionOrdering(u%Vx(iBladeNode,jBlade), phi, phi_lower, phi_upper)
   phiIn = phi
   
   do i = 1,size(phi_upper)   ! Need to potentially test 4 regions
       TestRegionResult = FindTestRegion(p, u, iBladeNode, jBlade, phi_lower(i), phi_upper(i), AFInfo, &
                        IsValidSolution, phiIn, f1, f_lower, f_upper, errStat2, errMsg2)
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName ) 
      
      if ( TestRegionResult == 1 ) then
         !............
         ! There is a zero in the solution region [phi_lower(i), phi_upper(i)] because the endpoints have residuals with different signs (SolutionRegion=1)
         ! We use Brent's Method to find the zero-residual solution in this region
                  
         call sub_brent(p, u, iBladeNode, jBlade, phi, phi_lower(i), phi_upper(i), AFInfo, IsValidSolution, ErrStat2, ErrMsg2, f_lower, f_upper)
            call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   
         if (IsValidSolution) then ! we have a valid BEMT solution
            ValidPhi = .true.
            exit
         end if    
      elseif (TestRegionResult == 3) then
         phi = phi_lower(i) !this boundary is a solution
         ValidPhi = .true.
         exit         
      elseif (TestRegionResult == 4) then
         phi = phi_upper(i) !this boundary is a solution
         ValidPhi = .true.
         exit         
      elseif (TestRegionResult == 0) then ! Special case where both end points return 0 residual; return value closest to 0 as the solution
         if (phi_lower(i) > 0.0_ReKi) then
            phi = phi_lower(i) !this boundary is a solution
            ValidPhi = .true.
            exit
         else
            phi = phi_upper(i) !this boundary is a solution
            ValidPhi = .true.
            exit
         end if         
      end if 
      
   end do
   

   if (.not. ValidPhi) then
      phi = ComputePhiWithInduction(u%Vx(iBladeNode,jBlade), u%Vy(iBladeNode,jBlade),  0.0_ReKi, 0.0_ReKi)
      
      if (abs(phi)>MsgLimit .and. abs(abs(phi)-PiBy2) > MsgLimit ) then
         if (FirstWarn) then
            call SetErrStat( ErrID_Info, 'There is no valid value of phi for these operating conditions: Vx = '//TRIM(Num2Lstr(u%Vx(iBladeNode,jBlade)))//&
               ', Vy = '//TRIM(Num2Lstr(u%Vy(iBladeNode,jBlade)))//', rlocal = '//TRIM(Num2Lstr(u%rLocal(iBladeNode,jBlade)))//', theta = '//TRIM(Num2Lstr(u%theta(iBladeNode,jBlade)))//', geometric phi = '//TRIM(Num2Lstr(phi)) &
               //'. This warning will not be repeated though the condition may persist. (See GeomPhi output channel.)', errStat, errMsg, RoutineName )
            FirstWarn = .false.
         end if !FirstWarn
      end if
   end if
   
         
end subroutine BEMT_UnCoupledSolve
!----------------------------------------------------------------------------------------------------------------------------------
function NodeText(i,j)
   integer(IntKi), intent(in) :: i ! node number
   integer(IntKi), intent(in) :: j ! blade number
   character(25)              :: NodeText
   
   NodeText = '(node '//trim(num2lstr(i))//', blade '//trim(num2lstr(j))//')'
end function NodeText
!----------------------------------------------------------------------------------------------------------------------------------
subroutine SetInputs_for_UA(phi, theta, axInduction, tanInduction, Vx, Vy, omega, chord, kinVisc, UserProp, u_UA)
   real(ReKi),                   intent(in   ) :: UserProp           ! User property (for 2D Airfoil interp)
   real(ReKi),                   intent(in   ) :: phi
   real(ReKi),                   intent(in   ) :: theta
   real(ReKi),                   intent(in   ) :: axInduction
   real(ReKi),                   intent(in   ) :: tanInduction
   real(ReKi),                   intent(in   ) :: Vx
   real(ReKi),                   intent(in   ) :: Vy
   real(ReKi),                   intent(in   ) :: omega ! aka PitchRate
   real(ReKi),                   intent(in   ) :: chord
   real(ReKi),                   intent(in   ) :: kinVisc
   type(UA_InputType),           intent(  out) :: u_UA


      ! ....... compute inputs to UA ...........
   u_UA%alpha   =  phi - theta  ! angle of attack
   u_UA%UserProp = UserProp
            
      ! Need to compute relative velocity at the aerodynamic center, including both axial and tangential induction 
      ! COMPUTE: u_UA%U, u_UA%Re, u_UA%v_ac
   call GetRelativeVelocity( axInduction, tanInduction, Vx, Vy, u_UA%U, u_UA%v_ac )
   call GetReynoldsNumber(   axInduction, tanInduction, Vx, Vy, chord, kinVisc, u_UA%Re)

   u_UA%v_ac(1) = sin(u_UA%alpha)*u_UA%U
   u_UA%v_ac(2) = cos(u_UA%alpha)*u_UA%U
   
   u_UA%omega = omega

end subroutine SetInputs_for_UA
!----------------------------------------------------------------------------------------------------------------------------------
subroutine SetInputs_for_UA_AllNodes(u, p, phi, axInduction, tanInduction, u_UA)
   type(BEMT_InputType),         intent(in   ) :: u
   type(BEMT_ParameterType),     intent(in   ) :: p
   real(ReKi),                   intent(inout) :: phi(:,:)
   real(ReKi),                   intent(in   ) :: axInduction(:,:)
   real(ReKi),                   intent(in   ) :: tanInduction(:,:)
   type(UA_InputType),           intent(  out) :: u_UA(:,:)

   INTEGER(IntKi)                              :: i, j

   if ( p%useInduction ) then
      !............................................
      ! we're recomputing phi because it may not be consistent with the computed axInduction or tanInduction values
      !............................................
      do j = 1,p%numBlades ! Loop through all blades
         do i = 1,p%numBladeNodes ! Loop through the blade nodes / elements
            phi(i,j) = ComputePhiWithInduction( u%Vx(i,j), u%Vy(i,j),  axInduction(i,j), tanInduction(i,j) )
         enddo             ! I - Blade nodes / elements
      enddo          ! J - All blades
   end if
   
   !............................................
   ! unsteady aero inputs (AoA, Vrel, Re)
   !............................................
   do j = 1,p%numBlades ! Loop through all blades
      do i = 1,p%numBladeNodes ! Loop through the blade nodes / elements
      
            ! Compute AoA, Re, Vrel (inputs to UA) based on current values of axInduction, tanInduction:
         call SetInputs_for_UA(phi(i,j), u%theta(i,j), axInduction(i,j), tanInduction(i,j), u%Vx(i,j), u%Vy(i,j), u%omega_z(i,j), p%chord(i,j), p%kinVisc, u%UserProp(i,j), u_UA(i,j))
         
      end do
   end do

end subroutine SetInputs_for_UA_AllNodes
!----------------------------------------------------------------------------------------------------------------------------------
#ifdef DEBUG_BEMT_RESIDUAL
subroutine WriteDEBUGValuesToFile(t, u, p, x, xd, z, OtherState, m, AFInfo)
   real(DbKi),                     intent(in   )  :: t           ! Current simulation time in seconds
   type(BEMT_InputType),           intent(in   )  :: u           ! Inputs at Time t
   type(BEMT_ParameterType),       intent(in   )  :: p           ! Parameters
   type(BEMT_ContinuousStateType), intent(in   )  :: x           ! Continuous states at t
   type(BEMT_DiscreteStateType),   intent(in   )  :: xd          ! Discrete states at t
   type(BEMT_ConstraintStateType), intent(in   )  :: z           ! Constraint states at t
   type(BEMT_OtherStateType),      intent(in   )  :: OtherState  ! Other states at t
   type(BEMT_MiscVarType),         intent(in   )  :: m           ! Misc/optimization variables
   type(AFI_ParameterType),        intent(in   )  :: AFInfo(:)   ! The airfoil parameter data

   character(ErrMsgLen)                           :: errMsg       ! temporary Error message if ErrStat /= ErrID_None
   integer(IntKi)                                 :: errStat      ! temporary Error status of the operation
   integer(IntKi)                                 :: OutInt
   integer(IntKi)                                 :: UnOut
   logical                                        :: ValidPhi
   integer(IntKi)                                 :: ii
   integer(IntKi), parameter                      :: numPhi = 3001
   real(ReKi)                                     :: phi
   real(ReKi)                                     :: fzero
   real(ReKi)                                     :: delta_phi
   real(ReKi)                                     :: axInd, tnInd
   
   ! I'm using persistent data.... don't copy this code anywhere else :)
   integer,  save         :: DEBUG_BLADE
   integer,  save         :: DEBUG_BLADENODE
   integer,  save         :: DEBUG_nStep = 1
   integer,  save         :: DEBUG_FILE_UNIT

!   character(*), parameter               :: RoutineName = 'BEMT_UnCoupledSolve'
   
   DEBUG_BLADE     = 1 !size(u%Vx,2)
   DEBUG_BLADENODE = 23 !max(1, size(u%Vx,1) / 2 )
   
   if (DEBUG_nStep == 1) then
      call GetNewUnit(DEBUG_FILE_UNIT, ErrStat, ErrMsg )
      call OpenFOutFile ( DEBUG_FILE_UNIT, "CheckBEMT.inputs.out", ErrStat, ErrMsg )
      
      WRITE(DEBUG_FILE_UNIT,'(A7,200(1x,A15))') "Step","IsGeomPhi","Time", "GeomPhi", "phi" &
                                             , "chi0", "omega", "Un_disk", "TSR" &
                                             , "theta" , "psi" &
                                             , "Vx" , "Vy"  &
                                             , "omega_z"  &
                                             , "rLocal" , "UserProp"  &
                                             , "AxInd", "TanInd"
!                                            , "Vx_elast_dot" , "Vy_elast_dot" &
          
   end if
   
   ! get geometric phi for these inputs:
   phi = ComputePhiWithInduction(u%Vx(DEBUG_BLADENODE,DEBUG_BLADE), u%Vy(DEBUG_BLADENODE,DEBUG_BLADE),  0.0_ReKi, 0.0_ReKi)
   
   if (OtherState%ValidPhi(DEBUG_BLADENODE,DEBUG_BLADE)) then
      OutInt = 0
   else
      OutInt = 1
   end if
   
   WRITE(DEBUG_FILE_UNIT, '(I7,1x,I15,200(1x,ES15.6))') DEBUG_nStep, OutInt, t, phi*R2D & 
                                             , z%phi(         DEBUG_BLADENODE,DEBUG_BLADE)*R2D & ! remember this does not have any corrections 
                                             , u%chi0, u%omega, u%Un_disk, u%TSR     &
                                             , u%theta(       DEBUG_BLADENODE,DEBUG_BLADE)     &
                                             , u%psi(                         DEBUG_BLADE)     &
                                             , u%Vx(          DEBUG_BLADENODE,DEBUG_BLADE)     &
                                             , u%Vy(          DEBUG_BLADENODE,DEBUG_BLADE)     &
                                             , u%Vz(          DEBUG_BLADENODE,DEBUG_BLADE)     &
                                             , u%omega_z(     DEBUG_BLADENODE,DEBUG_BLADE)     &
                                             , u%rLocal(      DEBUG_BLADENODE,DEBUG_BLADE)     &
                                             , u%UserProp(    DEBUG_BLADENODE,DEBUG_BLADE)     &
                                             , m%axInduction( DEBUG_BLADENODE,DEBUG_BLADE)     &
                                             , m%tanInduction(DEBUG_BLADENODE,DEBUG_BLADE)
! these are not always allocated
!                                             , u%Vx_elast_dot(DEBUG_BLADENODE,DEBUG_BLADE) &
!                                             , u%Vy_elast_dot(DEBUG_BLADENODE,DEBUG_BLADE) &
         
   ! now write the residual function to a separate file:
   if ((DEBUG_nStep >= 0).AND.(DEBUG_nStep <= 450000).AND.(MOD(DEBUG_nStep,25) == 0)) then
                                             
      call GetNewUnit( UnOut, ErrStat, ErrMsg )
      call OpenFOutFile ( UnOut, "CheckBEMT.residual."//trim(num2lstr(DEBUG_nStep))//".out", ErrStat, ErrMsg )
      !WRITE(UnOut, '(2(A15,1x),A2,5(1x,A15))') 'phi', 'residual', 'OK', 'AxialInd', 'TangentialInd', 'k', 'kp'

      do ii = 1, numPhi

         ! nonlinear mapping of ii --> phi
         phi = smoothStep( real(ii,ReKi), 1.0, -pi+BEMT_epsilon2, real(numPhi,ReKi)/2.0, 0.0_ReKi ) + smoothStep( real(ii,ReKi), real(numPhi,ReKi)/2.0, 0.0_ReKi, real(numPhi,ReKi), pi-BEMT_epsilon2 )
      
         fzero = BEMTU_InductionWithResidual(p, u, DEBUG_BLADENODE, DEBUG_BLADE, phi, AFInfo(p%AFIndx(DEBUG_BLADENODE,DEBUG_BLADE)), ValidPhi, errStat, errMsg, a=axInd, ap=tnInd )
         if (ValidPhi) then
            OutInt = 0
         else
            OutInt = 1
         end if
      
         WRITE(UnOut, '(2(ES15.6,1x),I2,8(1x,ES15.6))') phi, fzero, OutInt, axInd, tnInd
      
      end do
   
      CLOSE(UnOut)
   end if
   
   DEBUG_nStep = DEBUG_nStep + 1
   
end subroutine WriteDEBUGValuesToFile
#endif
!----------------------------------------------------------------------------------------------------------------------------------
end module BEMT
    
