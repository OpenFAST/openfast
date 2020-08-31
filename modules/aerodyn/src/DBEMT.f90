!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2015-2016  National Renewable Energy Laboratory
! Copyright (C) 2017  Envision Energy USA, LTD
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
!**********************************************************************************************************************************
module DBEMT
   
   use NWTC_Library   
   use DBEMT_Types
   
   implicit none 

private


   public :: DBEMT_Init
   public :: DBEMT_UpdateStates
   public :: DBEMT_CalcOutput
   public :: DBEMT_End
   PUBLIC :: DBEMT_CalcContStateDeriv             !  Tight coupling routine for computing derivatives of continuous states
   

   public :: DBEMT_ReInit
   public :: DBEMT_InitStates_AllNodes

   contains
   
   
subroutine DBEMT_ValidateInitInp(interval, InitInp, errStat, errMsg)
   real(DbKi),                      intent(in   ) :: interval      !< Coupling interval in seconds
   type(DBEMT_InitInputType),       intent(in   ) :: InitInp       !< Input data for initialization routine
   integer(IntKi),                  intent(  out) :: errStat       !< Error status of the operation
   character(*),                    intent(  out) :: errMsg        !< Error message if ErrStat /= ErrID_None
   
   character(*), parameter                     :: RoutineName = 'DBEMT_ValidateInitInp'
   real(ReKi)                                  :: rlocalMax
   integer(IntKi)                              :: i,j
     ! Initialize variables for this routine

   errStat = ErrID_None
   errMsg  = ""
   
   if ( interval <= sqrt(epsilon(1.0_ReKi)) ) call SetErrStat( ErrID_Fatal, " The timestep size for DBEMT (interval) must be larger than sqrt(epsilon).", ErrStat, ErrMsg, RoutineName)
   if ( (InitInp%DBEMT_Mod .ne. DBEMT_tauConst) .and. (InitInp%DBEMT_Mod .ne. DBEMT_tauVaries) .and. (InitInp%DBEMT_Mod .ne. DBEMT_cont_tauConst)) then
      call SetErrStat( ErrID_Fatal, " DBEMT_Mod must be set to 1, 2, or 3.", ErrStat, ErrMsg, RoutineName)
   end if
   
   if (InitInp%numBlades < 1) call SetErrStat( ErrID_Fatal, " InitInp%numBlades must set to 1 or more.", ErrStat, ErrMsg, RoutineName)
   if (InitInp%numNodes < 2) call SetErrStat( ErrID_Fatal, " InitInp%numNodes must set to 2 or more.", ErrStat, ErrMsg, RoutineName)


   if ( InitInp%DBEMT_Mod == DBEMT_tauConst .or. InitInp%DBEMT_Mod == DBEMT_cont_tauConst ) then
   
      if (InitInp%tau1_const <= 0.0_ReKi)  call SetErrStat( ErrID_Fatal, " InitInp%tau1_const must be greater than zero.", ErrStat, ErrMsg, RoutineName)
        ! Default = 0.33_ReKi
   
      if (.not. allocated(InitInp%rlocal) ) then
         call SetErrStat( ErrID_Fatal, " InitInput%rlocal must be allocated to size (InitInp%numNodes,InitInp%numBlades).", ErrStat, ErrMsg, RoutineName)
         return
      end if
   
      do j = 1,InitInp%numBlades
         rlocalMax = 0.0_ReKi
         do i= 1,InitInp%numNodes
            if (InitInp%rlocal(i,j) <= 0.0_ReKi ) then
               call SetErrStat( ErrID_Fatal, " Blades nodes must be located at a positive radial distance (rlocal) greater than zero.", ErrStat, ErrMsg, RoutineName)
               return
            end if
            
            rlocalMax = max(rlocalMax,InitInp%rlocal(i,j))
         end do
         if ( EqualRealNos(rlocalMax, 0.0_ReKi) ) call SetErrStat( ErrID_Fatal, " Blades must have nodes located at a radial distance (rlocal) greater than zero.  ", ErrStat, ErrMsg, RoutineName)
      end do
      
   end if

end subroutine DBEMT_ValidateInitInp


!----------------------------------------------------------------------------------------------------------------------------------   
!> This routine is called at the start of the simulation to perform initialization steps.
!! The parameters are set here and not changed during the simulation.
!! The initial states and initial guess for the input are defined.
subroutine DBEMT_Init( InitInp, u, p, x, OtherState, m, Interval, InitOut, ErrStat, ErrMsg )
!..................................................................................................................................

   type(DBEMT_InitInputType),       intent(in   ) :: InitInp       !< Input data for initialization routine
   type(DBEMT_InputType),           intent(  out) :: u             !< An initial guess for the input; input mesh must be defined
   type(DBEMT_ParameterType),       intent(  out) :: p             !< Parameters
   type(DBEMT_ContinuousStateType), intent(  out) :: x             !< Initial continuous states
   type(DBEMT_OtherStateType),      intent(  out) :: OtherState    !< Initial other/logical states
   type(DBEMT_MiscVarType),         intent(  out) :: m             !< Initial misc/optimization variables
   real(DbKi),                      intent(in   ) :: interval      !< Coupling interval in seconds: the rate that
                                                                   !!   (1) DBEMT_UpdateStates() is called in loose coupling &
                                                                   !!   (2) DBEMT_UpdateDiscState() is called in tight coupling.
                                                                   !!   Input is the suggested time from the glue code;
                                                                   !!   Output is the actual coupling interval that will be used
                                                                   !!   by the glue code.
   type(DBEMT_InitOutputType),      intent(  out) :: InitOut       !< Output for initialization routine
   integer(IntKi),                  intent(  out) :: errStat       !< Error status of the operation
   character(*),                    intent(  out) :: errMsg        !< Error message if ErrStat /= ErrID_None
   
   
     ! Local variables
   integer(IntKi)                              :: i,j             ! loop counter
   real(ReKi)                                  :: rTip          ! Maximum tip radius for all blades 
   integer(IntKi)                              :: errStat2      ! temporary error status of the operation
   character(ErrMsgLen)                        :: errMsg2       ! temporary error message 
   character(*), parameter                     :: RoutineName = 'DBEMT_Init'
   
   
   
      ! Initialize variables for this routine

   errStat = ErrID_None
   errMsg  = ""
   
      ! initialize InitOut to avoid compiler warnings
   InitOut%Ver = ProgDesc( 'DBEMT', '', '' ) 
   
   
  !if (p%DBEMT_Mod == DBEMT_none) return  ! DBEMT is turned off.
   
      ! Validate the Initialization inputs
   call DBEMT_ValidateInitInp(interval, InitInp, errStat2, errMsg2)
      call SetErrStat( errStat2, errMsg2, ErrStat, ErrMsg, RoutineName)
      if (errStat >= AbortErrLev) return
      
      ! Set parameter data using the initialization inputs
   
   p%DBEMT_Mod  = InitInp%DBEMT_Mod
   
   
   p%dt = interval
   p%numBlades  = InitInp%numBlades
   p%numNodes   = InitInp%numNodes
   p%k_0ye      = 0.6_ReKi
   !>>>>> bjj: these are unused:
   !p%c5         = 1.1_ReKi
   !p%c6         = 1.0_ReKi
   !p%c7         = 1.3_ReKi
   !p%c8         = 0.39_ReKi
   !p%c9         = 0.26_ReKi
   !<<<<<<
   p%tau1_const = InitInp%tau1_const  ! Default = 0.33_ReKi
   
   allocate(u%element(p%numNodes, p%numBlades), stat=ErrStat2)
      if (ErrStat2 /= 0 ) then
         call SetErrStat( ErrID_Fatal, " Error allocating u%element.", ErrStat, ErrMsg, RoutineName)
         return
      end if

   if (p%DBEMT_Mod == DBEMT_tauConst .or. p%DBEMT_Mod == DBEMT_cont_tauConst) then
      ! DBEMT_Mod = 1 and 3 use constant tau1 and tau2
      allocate( p%spanRatio(p%numNodes, p%numBlades), STAT=ErrStat2)
      if (ErrStat2 /= 0 ) then
         call SetErrStat( ErrID_Fatal, " Error allocating p%spanRatio.  ", ErrStat, ErrMsg, RoutineName)
         return
      end if
         ! compute the largest tip radius of all blades
      rTip = 0.0_ReKi
      do j = 1,p%numBlades
         do i= 1,p%numNodes
            rTip = max(rTip, InitInp%rlocal(i,j))
         end do
      end do
      
      do j = 1,p%numBlades
         do i= 1,p%numNodes
            p%spanRatio(i,j) = InitInp%rlocal(i,j)/rTip  ! normalized radial distance from center of rotation to node
         end do
      end do
   end if
   
      ! Initialize the continuous states
   allocate( x%element(p%numNodes, p%numBlades), STAT=ErrStat2) ! This is the axial and tangential induced velocity at node i on blade j
      if (ErrStat2 /= 0 ) then
         call SetErrStat( ErrID_Fatal, " Error allocating x%element.  ", ErrStat, ErrMsg, RoutineName)
         return
      end if

   if (p%DBEMT_Mod == DBEMT_cont_tauConst) then
      allocate( OtherState%n(p%numNodes, p%numBlades), STAT=ErrStat2)
         if (ErrStat2 /= 0 ) then
            call SetErrStat( ErrID_Fatal, " Error allocating OtherState%n.", ErrStat, ErrMsg, RoutineName)
            return
         end if
         
      do i=1,size(OtherState%xdot)
         call DBEMT_CopyContState( x, OtherState%xdot(i), MESH_NEWCOPY, ErrStat2, ErrMsg2 )
            if (ErrStat2 /= 0 ) then
               call SetErrStat( ErrID_Fatal, " Error allocating OtherState%xdot.", ErrStat, ErrMsg, RoutineName)
               return
            end if
      end do

      p%lin_nx = p%numNodes*p%numBlades*4 ! vind and vind_dot
   else
      p%lin_nx = 0
   end if

   allocate( OtherState%areStatesInitialized(p%numNodes, p%numBlades), STAT=ErrStat2)
      if (ErrStat2 /= 0 ) then
         call SetErrStat( ErrID_Fatal, " Error allocating OtherState%areStatesInitialized.  ", ErrStat, ErrMsg, RoutineName)
         return
      end if
   
   call DBEMT_ReInit(p, x,OtherState,m)
      
end subroutine DBEMT_Init

!..................................................................................................................................
subroutine DBEMT_ReInit( p, x, OtherState, m )

   type(DBEMT_ParameterType),       intent(in   ) :: p             !< parameters
   type(DBEMT_ContinuousStateType), intent(inout) :: x             !< Initial continuous states
   type(DBEMT_OtherStateType),      intent(inout) :: OtherState    !< Initial other/logical states
   type(DBEMT_MiscVarType),         intent(inout) :: m             !< Initial misc/optimization variables

   integer                                        :: i
   integer                                        :: j
   integer                                        :: n
   integer(IntKi)                                 :: ErrStat2
   character(ErrMsgLen)                           :: ErrMsg2
   
      ! Initialize variables for this routine

   do j=1,size(x%element,2)
      do i=1,size(x%element,1)
         x%element(i,j)%vind     = 0.0_ReKi
         x%element(i,j)%vind_dot = 0.0_ReKi
         x%element(i,j)%vind_1   = 0.0_ReKi
      end do
   end do

   OtherState%areStatesInitialized = .false.
   
   if (p%DBEMT_Mod == DBEMT_tauConst .or. p%DBEMT_Mod == DBEMT_cont_tauConst) then
      OtherState%tau1  = p%tau1_const
   else
      OtherState%tau1  = 0.0_ReKi
   end if
   
   if (allocated(OtherState%n)) then
      OtherState%n = -1
   
      do n=1,size(OtherState%xdot)
         do j=1,size(x%element,2)
            do i=1,size(x%element,1)
               call DBEMT_CopyElementContinuousStateType( x%element(i,j), OtherState%xdot(n)%element(i,j), MESH_UPDATECOPY, ErrStat2, ErrMsg2)
            end do
         end do
      end do
   end if
   
   m%FirstWarn_tau1 = .true.

end subroutine DBEMT_ReInit
!!----------------------------------------------------------------------------------------------------------------------------------
!> routine to initialize the states based on inputs
subroutine DBEMT_InitStates_AllNodes( u, p, x, OtherState )
   type(DBEMT_InputType),           intent(in   ) :: u          !< Inputs at t
   type(DBEMT_ParameterType),       intent(in   ) :: p          !< Parameters
   type(DBEMT_ContinuousStateType), intent(inout) :: x          !< Input: Continuous states at t;
                                                                !!   Output: Continuous states at t + Interval
   type(DBEMT_OtherStateType),      intent(inout) :: OtherState !< Other/logical states at t on input; at t+dt on output
   
   integer(IntKi)                                 :: i          !< blade node counter
   integer(IntKi)                                 :: j          !< blade counter

   
   do j=1,size(x%element,2)
      do i=1,size(x%element,1)
         call DBEMT_InitStates( i, j, u, p, x, OtherState )
      end do
   end do

end subroutine DBEMT_InitStates_AllNodes
!!----------------------------------------------------------------------------------------------------------------------------------
!> routine to initialize the states based on inputs
subroutine DBEMT_InitStates( i, j, u, p, x, OtherState )
   integer(IntKi),                  intent(in   ) :: i          !< blade node counter
   integer(IntKi),                  intent(in   ) :: j          !< blade counter
   type(DBEMT_InputType),           intent(in   ) :: u          !< Inputs at t
   type(DBEMT_ParameterType),       intent(in   ) :: p          !< Parameters
   type(DBEMT_ContinuousStateType), intent(inout) :: x          !< Input: Continuous states at t;
                                                                !!   Output: Continuous states at t + Interval
   type(DBEMT_OtherStateType),      intent(inout) :: OtherState !< Other/logical states at t on input; at t+dt on output

   if ( .not. OtherState%areStatesInitialized(i,j) ) then
      x%element(i,j)%vind(1) = u%element(i,j)%vind_s(1)
      x%element(i,j)%vind(2) = u%element(i,j)%vind_s(2)
      
      if (p%DBEMT_Mod == DBEMT_cont_tauConst) then
         x%element(i,j)%vind_dot(1) = u%element(i,j)%vind_s_dot(1)
         x%element(i,j)%vind_dot(2) = u%element(i,j)%vind_s_dot(2)
      else
         x%element(i,j)%vind_1(1) = u%element(i,j)%vind_s(1)
         x%element(i,j)%vind_1(2) = u%element(i,j)%vind_s(2)
      end if
      
      OtherState%areStatesInitialized(i,j) = .true.
      return
   end if
   
end subroutine DBEMT_InitStates
!!----------------------------------------------------------------------------------------------------------------------------------
!> Loose coupling routine for solving for constraint states, integrating continuous states, and updating discrete and other states.
!! Continuous, constraint, discrete, and other states are updated for t + Interval
subroutine DBEMT_UpdateStates( i, j, t, n, u, p, x, OtherState, m, errStat, errMsg )
!..................................................................................................................................
   integer(IntKi),                  intent(in   ) :: i          !< blade node counter
   integer(IntKi),                  intent(in   ) :: j          !< blade counter
   real(DbKi),                      intent(in   ) :: t          !< Current simulation time in seconds
   integer(IntKi),                  intent(in   ) :: n          !< Current simulation time step n = 0,1,...
   type(DBEMT_InputType),           intent(in   ) :: u(2)       !< Inputs at t and t+dt
   type(DBEMT_ParameterType),       intent(in   ) :: p          !< Parameters
   type(DBEMT_ContinuousStateType), intent(inout) :: x          !< Input: Continuous states at t;
                                                                !!   Output: Continuous states at t + Interval
   type(DBEMT_MiscVarType),         intent(inout) :: m          !< Initial misc/optimization variables
   type(DBEMT_OtherStateType),      intent(inout) :: OtherState !< Other/logical states at t on input; at t+dt on output
   integer(IntKi),                  intent(  out) :: errStat    !< Error status of the operation
   character(*),                    intent(  out) :: errMsg     !< Error message if ErrStat /= ErrID_None

   ! local variables
   real(ReKi)                                     :: A, B, C0, k_tau, C0_2 ! tau1_plus1, C_tau1, C, K1
   integer(IntKi)                                 :: indx
   real(DbKi)                                     :: utimes(2)
   
   TYPE(DBEMT_ElementInputType)                   :: u_elem(2)        !< Inputs at utimes
   
   character(*), parameter                        :: RoutineName = 'DBEMT_UpdateStates'
      
   ErrStat = ErrID_None
   ErrMsg  = ""
   
  !if (p%DBEMT_Mod == DBEMT_none) return  ! DBEMT is turned off.
   
   call ComputeTau1(u(1), p, m, OtherState%tau1, errStat, errMsg) ! place this here for DBEMTau1 output reasons
      if (errStat >= AbortErrLev) return
   call ComputeTau2(i, j, u(1)%element(i,j), p, OtherState%tau1, OtherState%tau2, k_tau)

   call DBEMT_InitStates( i, j, u(1), p, x, OtherState )

   if (p%DBEMT_Mod == DBEMT_cont_tauConst) then ! continuous formulation:
      utimes(1) = t
      utimes(2) = t + p%dt
      
      u_elem(1) = u(1)%element(i,j)
      u_elem(2) = u(2)%element(i,j)
      call DBEMT_ABM4( i, j, t, n, u_elem, utimes, p, x, OtherState, m, ErrStat, ErrMsg )
   
   else ! finite difference formulation:

      do indx=1,2 !Axial and tangential components.  jmj questions if this should be done for tangential induction.
      
         B =  ( u(2)%element(i,j)%vind_s(indx) - u(1)%element(i,j)%vind_s(indx) ) / p%dt                    ! Eqn. 1.17c  ! bjj: note that interpOrder will affect this numerical derivative
         A =    u(1)%element(i,j)%vind_s(indx) + B*p%k_0ye*OtherState%tau1                                  ! Eqn. 1.17b
      
         C0 = x%element(i,j)%vind_1(indx) - A + B*OtherState%tau1                                           ! Eqn. 1.21b
      
         x%element(i,j)%vind_1(indx) = C0*exp(-p%dt/OtherState%tau1) + A + B*(p%dt-OtherState%tau1)         ! Eqn. 1.21a, but this is using p%dt instead of t
      
         C0_2 = x%element(i,j)%vind(indx) - C0/(1-k_tau) - A + B*(OtherState%tau1 + OtherState%tau2)        ! Eqn. 1.24c
         x%element(i,j)%vind(indx) = C0_2*exp(-p%dt/OtherState%tau2) + A + B*(p%dt-OtherState%tau1-OtherState%tau2) &
                               + (C0/(1-k_tau))*exp(-p%dt/OtherState%tau1)                                  ! Eqn. 1.25
      
      
         !C      = (u(2)%element(i,j)%vind_s(indx) - u(1)%element(i,j)%vind_s(indx))/p%dt ! v_ind_s_future could come from BEMT update states, but this seems to violate the framework  !Eqn. 1.27C
         !A      = u(1)%element(i,j)%vind_s(indx) + C*p%k_0ye*tau1  !
         !B      = C*(p%k_0ye*C_tau1+1)
         !K1     = (A + A*C_tau1 - B*tau1) / (1+C_tau1)
         !C0     = ( x%element(i,j)%vind_1(indx) - K1 )*tau1**(1.0/C_tau1)
         !C0_2   = ( x%element(i,j)%vind(indx) - K1 - C0 / ((1-k_tau)**(1.0/C_tau1)) + B*k_tau*tau1/((1+C_tau1)*(1+k_tau*C_tau1)) ) * tau1**(1.0/(k_tau*C_tau1))
         !x%element(i,j)%vind(indx) = K1 + B*(p%dt-k_tau*tau1)/((1+C_tau1)*(1+k_tau*C_tau1)) + C0 / ((1-k_tau)*(tau1+C_tau1*p%dt)**(1.0/C_tau1)) + C0_2 / ((tau1 + C_tau1*p%dt)**(1.0/(k_tau*C_tau1)))
      end do
      
   end if
end subroutine DBEMT_UpdateStates

!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine computes the (rotor) value of tau1 for DBEMT
!----------------------------------------------------------------------------------------------------------------------------------
subroutine ComputeTau1(u, p, m, tau1, errStat, errMsg)
   type(DBEMT_InputType),         intent(in   ) :: u          !< Inputs at u(1)
   type(DBEMT_ParameterType),     intent(in   ) :: p          !< Parameters
   type(DBEMT_MiscVarType),       intent(inout) :: m          !< Initial misc/optimization variables
   real(ReKi)           ,         intent(  out) :: tau1       !< tau1 value used in DBEMT filter
   integer(IntKi),                intent(  out) :: errStat    !< Error status of the operation
   character(*),                  intent(  out) :: errMsg     !< Error message if ErrStat /= ErrID_None

   ! local variables
   real(ReKi)                                   :: temp

   real(ReKi)                                   :: AxInd_disk
   real(ReKi), parameter                        :: max_AxInd = 0.5_ReKi

   real(ReKi)                                   :: Un_disk
   real(ReKi), parameter                        :: min_Un = 0.1_ReKi
   
   character(*), parameter                      :: RoutineName = 'ComputeTau'


   ErrStat = ErrID_None
   ErrMsg  = ""

   if ( p%DBEMT_Mod == DBEMT_tauConst .or. p%DBEMT_Mod == DBEMT_cont_tauConst) then
      tau1       = p%tau1_const
   else
      
   ! We need to extrapolate the radius and disk velocity to the i+1 timestep
   ! We will grab the i+1 version of vind,s and the disk averaged induction by using the
   ! the already updated states of the BEMT module.
   
      !bjj: I believe u(1) is at t, which seems inconsistant with this comment
   
         ! Check if input values are valid for this formulation:
      if ( u%AxInd_disk > max_AxInd ) then
         AxInd_disk = max_AxInd
         if (m%FirstWarn_tau1) then
            call setErrStat( ErrID_Severe, 'Rotor-averaged axial induction factor is greater than '//trim(num2lstr(max_AxInd)) &
                    //'; limiting time-varying tau1. This message will not be repeated though the condition may persist.', &
                    ErrStat, ErrMsg, RoutineName ) ! don't print this error more than one time
            m%FirstWarn_tau1 = .false.
         end if
      else
         AxInd_disk = u%AxInd_disk
      end if
      
      if ( u%Un_disk < min_Un ) then
         Un_disk = min_Un
         if (m%FirstWarn_tau1) then
            call setErrStat( ErrID_Severe, 'Uninduced axial relative air speed, Un, is less than '//trim(num2lstr(min_Un)) &
                     // ' m/s; limiting time-varying tau1. This message will not be repeated though the ' &
                     //'condition may persist.', ErrStat, ErrMsg, RoutineName ) ! don't print this error more than one time
            m%FirstWarn_tau1 = .false.
         end if
      else
         Un_disk = u%Un_disk
      end if
      
      temp   = (1.0-1.3*AxInd_disk)*Un_disk
      
      tau1   = 1.1*u%R_disk/temp          ! Eqn. 1.2 (note that we've eliminated possibility of temp being 0)
      tau1   = min(tau1, 100.0_ReKi)      ! put a limit on this time constant so it isn't unrealistically long (particularly at initialization)
      
   end if

end subroutine ComputeTau1
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine computes the (rotor) value of tau1, tau2, and k_tau for DBEMT
!----------------------------------------------------------------------------------------------------------------------------------
subroutine ComputeTau2(i, j, u, p, tau1, tau2, k_tau_out)
   integer(IntKi),                intent(in   ) :: i          !< blade node counter
   integer(IntKi),                intent(in   ) :: j          !< blade counter
   type(DBEMT_ElementInputType),  intent(in   ) :: u          !< element inputs at u(1)
   type(DBEMT_ParameterType),     intent(in   ) :: p          !< Parameters
   real(ReKi)           ,         intent(in   ) :: tau1       !< tau1 value used in DBEMT filter
   real(ReKi)           ,         intent(  out) :: tau2       !< tau2 value used in DBEMT filter
   real(ReKi), optional ,         intent(  out) :: k_tau_out  !< k_tau value used in DBEMT filter

   ! local variables
   real(ReKi)                                   :: spanRatio  ! local version of r / R
   real(ReKi)                                   :: k_tau      


   if ( p%DBEMT_Mod == DBEMT_tauConst .or. p%DBEMT_Mod == DBEMT_cont_tauConst ) then
      spanRatio = p%spanRatio(i,j)
   else
      spanRatio = u%spanRatio
   end if
   
   k_tau = 0.39 - 0.26*spanRatio**2                                                 ! Eqn. 1.23b
   tau2  = k_tau*tau1                                                               ! Eqn. 1.7 or Eqn 1.23a
   
   if (present(k_tau_out) ) k_tau_out = k_tau
   
end subroutine ComputeTau2

!----------------------------------------------------------------------------------------------------------------------------------
!> Routine for computing outputs, used in both loose and tight coupling.
!! This subroutine is used to compute the output channels (motions and loads) and place them in the WriteOutput() array.
!! The descriptions of the output channels are not given here. Please see the included OutListParameters.xlsx sheet for
!! for a complete description of each output parameter.
subroutine DBEMT_CalcOutput( i, j, t, u, y_vind, p, x, OtherState, m, errStat, errMsg )
! NOTE: no matter how many channels are selected for output, all of the outputs are calculated
! All of the calculated output channels are placed into the m%AllOuts(:), while the channels selected for outputs are
! placed in the y%WriteOutput(:) array.
!..................................................................................................................................
   integer(IntKi),                  intent(in   ) :: i          !< blade node counter
   integer(IntKi),                  intent(in   ) :: j          !< blade counter
   real(DbKi),                      intent(in   ) :: t          !< Current simulation time in seconds
   type(DBEMT_InputType),           intent(in   ) :: u          !< Inputs at t 
   real(ReKi),                      intent(  out) :: y_vind(2)
   !type(DBEMT_OutputType),          intent(inout) :: y          !< Inputs at utimes (out only for mesh record-keeping in ExtrapInterp routine)
   type(DBEMT_ParameterType),       intent(in   ) :: p          !< Parameters
   type(DBEMT_ContinuousStateType), intent(in   ) :: x          !< Input: Continuous states at t;
                                                                !!   Output: Continuous states at t + Interval
   type(DBEMT_MiscVarType),         intent(inout) :: m          !< Initial misc/optimization variables
   type(DBEMT_OtherStateType),      intent(in   ) :: OtherState !< Other/logical states at t

   integer(IntKi),                  intent(  out) :: errStat    !< Error status of the operation
   character(*),                    intent(  out) :: errMsg     !< Error message if ErrStat /= ErrID_None
  
   ! local variables
   
   character(*), parameter                      :: RoutineName = 'DBEMT_CalcOutput'
      
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   if ( .not. OtherState%areStatesInitialized(i,j) ) then
      y_vind = u%element(i,j)%vind_s
   else     
      y_vind = x%element(i,j)%vind ! array copy
   end if
   
      
end subroutine DBEMT_CalcOutput

!----------------------------------------------------------------------------------------------------------------------------------
!> Tight coupling routine for computing derivatives of continuous states.
SUBROUTINE DBEMT_CalcContStateDeriv( i, j, t, u, p, x, OtherState, m, dxdt, ErrStat, ErrMsg )
!..................................................................................................................................

   INTEGER(IntKi),                  INTENT(IN   )  :: i           !< blade node counter
   INTEGER(IntKi),                  INTENT(IN   )  :: j           !< blade counter
   REAL(DbKi),                      INTENT(IN   )  :: t           !< Current simulation time in seconds
   TYPE(DBEMT_ElementInputType),    INTENT(IN   )  :: u           !< Inputs at t
   TYPE(DBEMT_ParameterType),       INTENT(IN   )  :: p           !< Parameters
   TYPE(DBEMT_ElementContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at t (and i and j)
   TYPE(DBEMT_OtherStateType),      INTENT(IN   )  :: OtherState  !< Other states
   TYPE(DBEMT_MiscVarType),         INTENT(INOUT)  :: m           !< Misc variables for optimization (not copied in glue code)
   TYPE(DBEMT_ElementContinuousStateType), INTENT(OUT)  :: dxdt        !< Continuous state derivatives at t (note that since we are operating on only a portion of the continuous states, I have made a separate, smaller type to avoid excessive memory allocation/deallocation)
   INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      ! LOCAL variables
   CHARACTER(*), PARAMETER                         :: RoutineName = 'DBEMT_CalcContStateDeriv'
   
   REAL(ReKi)                                      :: tauConst
   REAL(ReKi)                                      :: tau1
   REAL(ReKi)                                      :: tau2
   
      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ""
   
         ! Compute the first time derivatives of the continuous states here:

   if (p%DBEMT_Mod /= DBEMT_cont_tauConst) then
      call SetErrStat(ErrID_Fatal,"Continuous state derivatives cannot be calculated unless DBEMT_Mod is 3.",ErrStat,ErrMsg,RoutineName)
      return
   end if

   tau1 = p%tau1_const
   !call ComputeTau1( u, p, m, tau1, errStat, errMsg)
   call ComputeTau2(i, j, u, p, tau1, tau2)

   ! Implement Equation 37 from E.Branlard 16-Dec-2019 doc:
   
   dxdt%vind = x%vind_dot
   
   tauConst = -1.0_ReKi/(tau1 * tau2)
   
   dxdt%vind_dot = tauConst * ( x%vind(:) + (tau1 + tau2)*x%vind_dot(:) &
                            - u%vind_s(:)  - p%k_0ye*tau1*u%vind_s_dot(:) )
                            
END SUBROUTINE DBEMT_CalcContStateDeriv
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine implements the fourth-order Runge-Kutta Method (RK4) for numerically integrating ordinary differential equations:
!!
!!   Let f(t, x) = xdot denote the time (t) derivative of the continuous states (x). 
!!   Define constants k1, k2, k3, and k4 as 
!!        k1 = dt * f(t        , x_t        )
!!        k2 = dt * f(t + dt/2 , x_t + k1/2 )
!!        k3 = dt * f(t + dt/2 , x_t + k2/2 ), and
!!        k4 = dt * f(t + dt   , x_t + k3   ).
!!   Then the continuous states at t = t + dt are
!!        x_(t+dt) = x_t + k1/6 + k2/3 + k3/3 + k4/6 + O(dt^5)
!!
!! For details, see:
!! Press, W. H.; Flannery, B. P.; Teukolsky, S. A.; and Vetterling, W. T. "Runge-Kutta Method" and "Adaptive Step Size Control for 
!!   Runge-Kutta." Sections 16.1 and 16.2 in Numerical Recipes in FORTRAN: The Art of Scientific Computing, 2nd ed. Cambridge, England: 
!!   Cambridge University Press, pp. 704-716, 1992.
SUBROUTINE DBEMT_RK4( i, j, t, n, u, utimes, p, x, OtherState, m, ErrStat, ErrMsg )
!..................................................................................................................................

   integer(IntKi),                  INTENT(IN   )  :: i           !< blade node counter
   integer(IntKi),                  INTENT(IN   )  :: j           !< blade counter
   REAL(DbKi),                      INTENT(IN   )  :: t           !< Current simulation time in seconds
   INTEGER(IntKi),                  INTENT(IN   )  :: n           !< time step number
   TYPE(DBEMT_ElementInputType),    INTENT(IN   )  :: u(:)        !< Inputs at utimes
   REAL(DbKi),                      INTENT(IN   )  :: utimes(:)   !< times of input
   TYPE(DBEMT_ParameterType),       INTENT(IN   )  :: p           !< Parameters
   TYPE(DBEMT_ContinuousStateType), INTENT(INOUT)  :: x           !< Continuous states at t on input at t + dt on output
   TYPE(DBEMT_OtherStateType),      INTENT(INOUT)  :: OtherState  !< Other states
   TYPE(DBEMT_MiscVarType),         INTENT(INOUT)  :: m           !< misc/optimization variables
   INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

   ! local variables
         
   TYPE(DBEMT_ElementContinuousStateType)          :: k1          ! RK4 constant; see above
   TYPE(DBEMT_ElementContinuousStateType)          :: k2          ! RK4 constant; see above 
   TYPE(DBEMT_ElementContinuousStateType)          :: k3          ! RK4 constant; see above 
   TYPE(DBEMT_ElementContinuousStateType)          :: k4          ! RK4 constant; see above 
   TYPE(DBEMT_ElementContinuousStateType)          :: x_tmp       ! Holds temporary modification to x
   TYPE(DBEMT_ElementInputType)                    :: u_interp    ! interpolated value of inputs
   
   REAL(DbKi)                                      :: TPlusHalfDt
   REAL(DbKi)                                      :: TPlusDt
   INTEGER(IntKi)                                  :: ErrStat2    ! local error status
   CHARACTER(ErrMsgLen)                            :: ErrMsg2     ! local error message (ErrMsg)
   CHARACTER(*), PARAMETER                         :: RoutineName = 'DBEMT_RK4'
      
      
   !NOTE: the error handling here assumes that we do not have any allocatable data in the inputs (u_interp) to be concerned with.
   !      Also, We assume that if there is going to be an error in DBEMT_CalcContStateDeriv, it will happen only on the first call 
   !      to the routine.
   
      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 
      
      ! interpolate u to find u_interp = u(t)
      CALL DBEMT_ElementInputType_ExtrapInterp( u, utimes, u_interp, t, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) RETURN

      x_tmp     = x%element(i,j)

      ! find xdot at t
      CALL DBEMT_CalcContStateDeriv( i, j, t, u_interp, p, x_tmp, OtherState, m, k1, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) RETURN
         

      k1%vind     = p%dt * k1%vind
      k1%vind_dot = p%dt * k1%vind_dot
  
      x_tmp%vind     = x%element(i,j)%vind     + 0.5 * k1%vind
      x_tmp%vind_dot = x%element(i,j)%vind_dot + 0.5 * k1%vind_dot

      ! interpolate u to find u_interp = u(t + dt/2)
      TPlusHalfDt = t+0.5_DbKi*p%dt
      CALL DBEMT_ElementInputType_ExtrapInterp(u, utimes, u_interp, TPlusHalfDt, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) RETURN

      ! find xdot at t + dt/2
      CALL DBEMT_CalcContStateDeriv( i, j, TPlusHalfDt, u_interp, p, x_tmp, OtherState, m, k2, ErrStat2, ErrMsg2 )

      k2%vind     = p%dt * k2%vind
      k2%vind_dot = p%dt * k2%vind_dot

      x_tmp%vind     = x%element(i,j)%vind     + 0.5 * k2%vind
      x_tmp%vind_dot = x%element(i,j)%vind_dot + 0.5 * k2%vind_dot

      ! find xdot at t + dt/2 (note x_tmp has changed)
      CALL DBEMT_CalcContStateDeriv( i, j, TPlusHalfDt, u_interp, p, x_tmp, OtherState, m, k3, ErrStat2, ErrMsg2 )

      k3%vind     = p%dt * k3%vind
      k3%vind_dot = p%dt * k3%vind_dot

      x_tmp%vind     = x%element(i,j)%vind     + k3%vind
      x_tmp%vind_dot = x%element(i,j)%vind_dot + k3%vind_dot

      ! interpolate u to find u_interp = u(t + dt)
      TPlusDt = t + p%dt
      CALL DBEMT_ElementInputType_ExtrapInterp(u, utimes, u_interp, TPlusDt, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) RETURN

      ! find xdot at t + dt
      CALL DBEMT_CalcContStateDeriv( i, j, TPlusDt, u_interp, p, x_tmp, OtherState, m, k4, ErrStat2, ErrMsg2 )

      k4%vind     = p%dt * k4%vind
      k4%vind_dot = p%dt * k4%vind_dot

      x%element(i,j)%vind     = x%element(i,j)%vind     + ( k1%vind     + 2. * k2%vind     + 2. * k3%vind     + k4%vind     ) / 6.
      x%element(i,j)%vind_dot = x%element(i,j)%vind_dot + ( k1%vind_dot + 2. * k2%vind_dot + 2. * k3%vind_dot + k4%vind_dot ) / 6.

END SUBROUTINE DBEMT_RK4
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine implements the fourth-order Adams-Bashforth Method (RK4) for numerically integrating ordinary differential 
!! equations:
!!
!!   Let f(t, x) = xdot denote the time (t) derivative of the continuous states (x). 
!!
!!   x(t+dt) = x(t)  + (dt / 24.) * ( 55.*f(t,x) - 59.*f(t-dt,x) + 37.*f(t-2.*dt,x) - 9.*f(t-3.*dt,x) )
!!
!!  See, e.g.,
!!  http://en.wikipedia.org/wiki/Linear_multistep_method
!!
!!  or
!!
!!  K. E. Atkinson, "An Introduction to Numerical Analysis", 1989, John Wiley & Sons, Inc, Second Edition.
SUBROUTINE DBEMT_AB4( i, j, t, n, u, utimes, p, x, OtherState, m, ErrStat, ErrMsg )
!..................................................................................................................................

   integer(IntKi),                  INTENT(IN   )  :: i           !< blade node counter
   integer(IntKi),                  INTENT(IN   )  :: j           !< blade counter
   REAL(DbKi),                      INTENT(IN   )  :: t           !< Current simulation time in seconds
   INTEGER(IntKi),                  INTENT(IN   )  :: n           !< time step number
   TYPE(DBEMT_ElementInputType),    INTENT(IN   )  :: u(:)        !< Inputs at utimes
   REAL(DbKi),                      INTENT(IN   )  :: utimes(:)   !< times of input
   TYPE(DBEMT_ParameterType),       INTENT(IN   )  :: p           !< Parameters
   TYPE(DBEMT_ContinuousStateType), INTENT(INOUT)  :: x           !< Continuous states at t on input at t + dt on output
   TYPE(DBEMT_OtherStateType),      INTENT(INOUT)  :: OtherState  !< Other states
   TYPE(DBEMT_MiscVarType),         INTENT(INOUT)  :: m           !< misc/optimization variables
   INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None


   ! local variables
   TYPE(DBEMT_ElementInputType)                    :: u_interp
   TYPE(DBEMT_ElementContinuousStateType)          :: x_tmp
         
   INTEGER(IntKi)                                  :: ErrStat2    ! local error status
   CHARACTER(ErrMsgLen)                            :: ErrMsg2     ! local error message (ErrMsg)
   CHARACTER(*), PARAMETER                         :: RoutineName = 'DBEMT_AB4'


      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 
      
      if (OtherState%n(i,j) < n) then

         OtherState%n(i,j) = n
            ! these types don't contain pointers or allocatable arrays, so we can copy them directly here:
         OtherState%xdot(4)%element(i,j) = OtherState%xdot(3)%element(i,j)
         OtherState%xdot(3)%element(i,j) = OtherState%xdot(2)%element(i,j)
         OtherState%xdot(2)%element(i,j) = OtherState%xdot(1)%element(i,j)

      elseif (OtherState%n(i,j) > n) then

         CALL SetErrStat(ErrID_Fatal,'Backing up in time is not supported with a multistep method.',ErrStat,ErrMsg,RoutineName)
         RETURN

      endif
      
      ! need xdot at t, get inputs at t
      CALL DBEMT_ElementInputType_ExtrapInterp(u, utimes, u_interp, t, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) RETURN
         
      x_tmp     = x%element(i,j)
      
         ! calculate OtherState%xdot( 1 )%element(i,j)
      CALL DBEMT_CalcContStateDeriv( i, j, t, u_interp, p, x_tmp, OtherState, m, OtherState%xdot( 1 )%element(i,j), ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) RETURN
         
                                                    
      if (n <= 2) then

         CALL DBEMT_RK4(i, j, t, n, u, utimes, p, x, OtherState, m, ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            IF ( ErrStat >= AbortErrLev ) RETURN

      else
         
         x%element(i,j)%vind     = x%element(i,j)%vind     + p%DT/24. * ( 55.*OtherState%xdot(1)%element(i,j)%vind     - 59.*OtherState%xdot(2)%element(i,j)%vind   &
                                                                        + 37.*OtherState%xdot(3)%element(i,j)%vind     -  9.*OtherState%xdot(4)%element(i,j)%vind )

         x%element(i,j)%vind_dot = x%element(i,j)%vind_dot + p%DT/24. * ( 55.*OtherState%xdot(1)%element(i,j)%vind_dot - 59.*OtherState%xdot(2)%element(i,j)%vind_dot  &
                                                                        + 37.*OtherState%xdot(3)%element(i,j)%vind_dot -  9.*OtherState%xdot(4)%element(i,j)%vind_dot )

      endif

      
END SUBROUTINE DBEMT_AB4
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine implements the fourth-order Adams-Bashforth-Moulton Method (RK4) for numerically integrating ordinary 
!! differential equations:
!!
!!   Let f(t, x) = xdot denote the time (t) derivative of the continuous states (x). 
!!
!!   Adams-Bashforth Predictor: \n
!!   x^p(t+dt) = x(t)  + (dt / 24.) * ( 55.*f(t,x) - 59.*f(t-dt,x) + 37.*f(t-2.*dt,x) - 9.*f(t-3.*dt,x) )
!!
!!   Adams-Moulton Corrector: \n
!!   x(t+dt) = x(t)  + (dt / 24.) * ( 9.*f(t+dt,x^p) + 19.*f(t,x) - 5.*f(t-dt,x) + 1.*f(t-2.*dt,x) )
!!
!!  See, e.g.,
!!  http://en.wikipedia.org/wiki/Linear_multistep_method
!!
!!  or
!!
!!  K. E. Atkinson, "An Introduction to Numerical Analysis", 1989, John Wiley & Sons, Inc, Second Edition.
SUBROUTINE DBEMT_ABM4( i, j, t, n, u, utimes, p, x, OtherState, m, ErrStat, ErrMsg )
!..................................................................................................................................

   integer(IntKi),                  INTENT(IN   )  :: i           !< blade node counter
   integer(IntKi),                  INTENT(IN   )  :: j           !< blade counter
   REAL(DbKi),                      INTENT(IN   )  :: t           !< Current simulation time in seconds
   INTEGER(IntKi),                  INTENT(IN   )  :: n           !< time step number
   TYPE(DBEMT_ElementInputType),    INTENT(IN   )  :: u(:)        !< Inputs at utimes
   REAL(DbKi),                      INTENT(IN   )  :: utimes(:)   !< times of input
   TYPE(DBEMT_ParameterType),       INTENT(IN   )  :: p           !< Parameters
   TYPE(DBEMT_ContinuousStateType), INTENT(INOUT)  :: x           !< Continuous states at t on input at t + dt on output
   TYPE(DBEMT_OtherStateType),      INTENT(INOUT)  :: OtherState  !< Other states
   TYPE(DBEMT_MiscVarType),         INTENT(INOUT)  :: m           !< misc/optimization variables
   INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

   ! local variables

   TYPE(DBEMT_ElementInputType)                    :: u_interp    ! Inputs at t
   TYPE(DBEMT_ElementContinuousStateType)          :: x_in        ! Continuous states at t
   TYPE(DBEMT_ElementContinuousStateType)          :: xdot_pred   ! Derivative of continuous states at t

   INTEGER(IntKi)                                  :: ErrStat2    ! local error status
   CHARACTER(ErrMsgLen)                            :: ErrMsg2     ! local error message (ErrMsg)
   CHARACTER(*), PARAMETER                         :: RoutineName = 'DBEMT_ABM4'
      
      
      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 

      ! save copy of x(t):
      x_in = x%element(i,j)
      

         ! predict: (note that we are overwritting x%element(i,j)%vind and x%element(i,j)%vind_dot here):
      CALL DBEMT_AB4( i, j, t, n, u, utimes, p, x, OtherState, m, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) RETURN

      if (n > 2_IntKi) then
      
            ! correct:
         
         CALL DBEMT_ElementInputType_ExtrapInterp(u, utimes, u_interp, t + p%dt, ErrStat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            IF ( ErrStat >= AbortErrLev ) RETURN

         CALL DBEMT_CalcContStateDeriv(i, j, t + p%dt, u_interp, p, x%element(i,j), OtherState, m, xdot_pred, ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            IF ( ErrStat >= AbortErrLev ) RETURN

         
         x%element(i,j)%vind     = x_in%vind     + p%DT/24. * ( 9. * xdot_pred%vind     + 19. * OtherState%xdot(1)%element(i,j)%vind &
                                                                                        -  5. * OtherState%xdot(2)%element(i,j)%vind &
                                                                                        +  1. * OtherState%xdot(3)%element(i,j)%vind )

         x%element(i,j)%vind_dot = x_in%vind_dot + p%DT/24. * ( 9. * xdot_pred%vind_dot + 19. * OtherState%xdot(1)%element(i,j)%vind_dot &
                                                                                        -  5. * OtherState%xdot(2)%element(i,j)%vind_dot &
                                                                                        +  1. * OtherState%xdot(3)%element(i,j)%vind_dot )
      endif
      
END SUBROUTINE DBEMT_ABM4
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is called at the end of the simulation.
subroutine DBEMT_End( u, p, x, OtherState, m, ErrStat, ErrMsg )
!..................................................................................................................................

      TYPE(DBEMT_InputType),           INTENT(INOUT)  :: u(2)           !< System inputs
      TYPE(DBEMT_ParameterType),       INTENT(INOUT)  :: p              !< Parameters
      TYPE(DBEMT_ContinuousStateType), INTENT(INOUT)  :: x              !< Continuous states
      type(DBEMT_MiscVarType),         intent(inout)  :: m              !< Initial misc/optimization variables
      type(DBEMT_OtherStateType),      intent(inout)  :: OtherState     !< Initial misc/optimization variables
      INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat        !< Error status of the operation
      CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg         !< Error message if ErrStat /= ErrID_None



         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""


         ! Place any last minute operations or calculations here:


         ! Close files here:



         ! Destroy the input data:

      CALL DBEMT_DestroyInput( u(1), ErrStat, ErrMsg )
      CALL DBEMT_DestroyInput( u(2), ErrStat, ErrMsg )


         ! Destroy the parameter data:

      CALL DBEMT_DestroyParam( p, ErrStat, ErrMsg )


         ! Destroy the state data:

      CALL DBEMT_DestroyContState(   x,           ErrStat, ErrMsg )


      CALL DBEMT_DestroyMisc(   m,           ErrStat, ErrMsg )
      CALL DBEMT_DestroyOtherState(   OtherState,           ErrStat, ErrMsg )


END SUBROUTINE DBEMT_End

end module DBEMT