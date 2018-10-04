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

   public :: DBEMT_ReInit

   contains
   
   
subroutine DBEMT_ValidateInitInp(interval, InitInp, errStat, errMsg)
   real(DbKi),                      intent(inout) :: interval      !< Coupling interval in seconds: the rate that
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
   if ( (InitInp%DBEMT_Mod .ne. DBEMT_tauConst) .and. (InitInp%DBEMT_Mod .ne. DBEMT_tauVaries) ) call SetErrStat( ErrID_Fatal, " DBEMT_Mod must be set to 1 or 2.", ErrStat, ErrMsg, RoutineName)
   
   if (InitInp%numBlades < 1) call SetErrStat( ErrID_Fatal, " InitInp%numBlades must set to 1 or more.", ErrStat, ErrMsg, RoutineName)
   if (InitInp%numNodes < 2) call SetErrStat( ErrID_Fatal, " InitInp%numNodes must set to 2 or more.", ErrStat, ErrMsg, RoutineName)
  
   if ( (InitInp%DBEMT_Mod == DBEMT_tauConst) )then
   
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
   real(DbKi),                      intent(inout) :: interval      !< Coupling interval in seconds: the rate that
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
   
      ! initialize InitOut and u to avoid compiler warnings
   InitOut%Ver = ProgDesc( 'DBEMT', 'v1.00.00', '28-Jul-2017' ) 
   u%spanRatio = 0.0_ReKi
   
   
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
   
   if (p%DBEMT_Mod == DBEMT_tauConst) then
      ! DBEMT_Mod = 1 uses constant tau1 and tau2
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
   allocate( x%vind(2,p%numNodes, p%numBlades), STAT=ErrStat2) ! This is the axial and tangential induced velocity at node i on blade j
      if (ErrStat2 /= 0 ) then
         call SetErrStat( ErrID_Fatal, " Error allocating x%vind.  ", ErrStat, ErrMsg, RoutineName)
         return
      end if

   allocate( x%vind_1(2,p%numNodes, p%numBlades), STAT=ErrStat2) ! This is the axial and tangential induced velocity at node i on blade j
      if (ErrStat2 /= 0 ) then
         call SetErrStat( ErrID_Fatal, " Error allocating x%vind_1.  ", ErrStat, ErrMsg, RoutineName)
         return
      end if
   
   allocate( OtherState%areStatesInitialized(p%numNodes, p%numBlades), STAT=ErrStat2)
      if (ErrStat2 /= 0 ) then
         call SetErrStat( ErrID_Fatal, " Error allocating OtherState%areStatesInitialized.  ", ErrStat, ErrMsg, RoutineName)
         return
      end if
   
   call DBEMT_ReInit(x,OtherState,m)
      
end subroutine DBEMT_Init

!..................................................................................................................................
subroutine DBEMT_ReInit( x, OtherState, m )

   type(DBEMT_ContinuousStateType), intent(inout) :: x             !< Initial continuous states
   type(DBEMT_OtherStateType),      intent(inout) :: OtherState    !< Initial other/logical states
   type(DBEMT_MiscVarType),         intent(inout) :: m             !< Initial misc/optimization variables

      ! Initialize variables for this routine

   x%vind = 0.0                      ! This is the axial and tangential induced velocity
   x%vind_1 = 0.0                    ! This is the axial and tangential induced velocity
   OtherState%areStatesInitialized = .false.
   OtherState%tau1 = 0.0_ReKi
   m%FirstWarn_tau1 = .true.

end subroutine DBEMT_ReInit

!!----------------------------------------------------------------------------------------------------------------------------------
!> Loose coupling routine for solving for constraint states, integrating continuous states, and updating discrete and other states.
!! Continuous, constraint, discrete, and other states are updated for t + Interval
subroutine DBEMT_UpdateStates( i, j, t, u,  p, x, OtherState, m, errStat, errMsg )
!..................................................................................................................................
   integer(IntKi),                  intent(in   ) :: i          !< blade node counter
   integer(IntKi),                  intent(in   ) :: j          !< blade counter
   real(DbKi),                      intent(in   ) :: t          !< Current simulation time in seconds
   type(DBEMT_InputType),           intent(inout) :: u(2)       !< Inputs at utimes (out only for mesh record-keeping in ExtrapInterp routine)
   type(DBEMT_ParameterType),       intent(in   ) :: p          !< Parameters
   type(DBEMT_ContinuousStateType), intent(inout) :: x          !< Input: Continuous states at t;
                                                                !!   Output: Continuous states at t + Interval
   type(DBEMT_MiscVarType),         intent(inout) :: m          !< Initial misc/optimization variables
   type(DBEMT_OtherStateType),      intent(inout) :: OtherState !< Other/logical states at t on input; at t+dt on output
   integer(IntKi),                  intent(  out) :: errStat    !< Error status of the operation
   character(*),                    intent(  out) :: errMsg     !< Error message if ErrStat /= ErrID_None

   ! local variables
   real(ReKi)                                   :: spanRatio       ! local version of r / R
   real(ReKi)                                   :: temp, tau2 , A, B, C0, k_tau, C0_2 ! tau1_plus1, C_tau1, C, K1
   real(ReKi)                                   :: Un_disk
   real(ReKi)                                   :: AxInd_disk
   integer(IntKi)                               :: indx
   
   character(*), parameter                      :: RoutineName = 'DBEMT_UpdateStates'
      
   ErrStat = ErrID_None
   ErrMsg  = ""
   
  !if (p%DBEMT_Mod == DBEMT_none) return  ! DBEMT is turned off.
   
   call ComputeTau1(u(1), p, m, OtherState%tau1, errStat, errMsg) ! place this here for DBEMTau1 output reasons
      if (errStat >= AbortErrLev) return

   if ( .not. OtherState%areStatesInitialized(i,j) ) then
      x%vind_1(1,i,j) = u(1)%vind_s(1)
      x%vind_1(2,i,j) = u(1)%vind_s(2)
      x%vind(1,i,j) = u(1)%vind_s(1)
      x%vind(2,i,j) = u(1)%vind_s(2)
      OtherState%areStatesInitialized(i,j) = .true.
      return
   end if
   
   
   if ( p%DBEMT_Mod == DBEMT_tauConst ) then
      spanRatio   = p%spanRatio(i,j)
   else
      spanRatio = u(1)%spanRatio
   end if
   

   do indx=1,2 !Axial and tangential components.  jmj questions if this should be done for tangential induction.
      
      
      ! TODO: Deal with initialization so that we avoid spikes???
      
      B =  ( u(2)%vind_s(indx) - u(1)%vind_s(indx) ) / p%dt                                              ! Eqn. 1.17c  ! bjj: note that interpOrder will affect this numerical derivative
      A = u(1)%vind_s(indx) + B*p%k_0ye*OtherState%tau1                                                  ! Eqn. 1.17b
      
      C0 = x%vind_1(indx,i,j) - A + B*OtherState%tau1                                                    ! Eqn. 1.21b
      
      x%vind_1(indx,i,j) = C0*exp(-p%dt/OtherState%tau1) + A + B*(p%dt-OtherState%tau1)                  ! Eqn. 1.21a, but this is using p%dt instead of t
      
      k_tau = 0.39 - 0.26*spanRatio**2                                                                   ! Eqn. 1.23b
      tau2 = k_tau*OtherState%tau1                                                                       ! Eqn. 1.7 or Eqn 1.23a
      C0_2 = x%vind(indx,i,j) - C0/(1-k_tau) - A + B*(OtherState%tau1 + tau2)                            ! Eqn. 1.24c 
      x%vind(indx,i,j) = C0_2*exp(-p%dt/tau2) + A + B*(p%dt-OtherState%tau1-tau2) + (C0/(1-k_tau))*exp(-p%dt/OtherState%tau1)  ! Eqn. 1.25
      
      
      !C      = (u(2)%vind_s(indx) - u(1)%vind_s(indx))/p%dt ! v_ind_s_future could come from BEMT update states, but this seems to violate the framework  !Eqn. 1.27C
      !A      = u(1)%vind_s(indx) + C*p%k_0ye*tau1  !
      !B      = C*(p%k_0ye*C_tau1+1)
      !K1     = (A + A*C_tau1 - B*tau1) / (1+C_tau1)
      !C0     = ( x%vind_1(indx,i,j) - K1 )*tau1**(1.0/C_tau1)
      !k_tau  = 0.39 - 0.26*(spanRatio)**2
      !C0_2   = ( x%vind(indx,i,j) - K1 - C0 / ((1-k_tau)**(1.0/C_tau1)) + B*k_tau*tau1/((1+C_tau1)*(1+k_tau*C_tau1)) ) * tau1**(1.0/(k_tau*C_tau1))
      !x%vind(indx,i,j) = K1 + B*(p%dt-k_tau*tau1)/((1+C_tau1)*(1+k_tau*C_tau1)) + C0 / ((1-k_tau)*(tau1+C_tau1*p%dt)**(1.0/C_tau1)) + C0_2 / ((tau1 + C_tau1*p%dt)**(1.0/(k_tau*C_tau1)))
   end do
end subroutine DBEMT_UpdateStates

!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine computes the (rotor) value of tau1 for DBEMT
!----------------------------------------------------------------------------------------------------------------------------------
subroutine ComputeTau1(u, p, m, tau1, errStat, errMsg)
   type(DBEMT_InputType),         intent(inout) :: u          !< Inputs at u(1)
   type(DBEMT_ParameterType),     intent(in   ) :: p          !< Parameters
   type(DBEMT_MiscVarType),       intent(inout) :: m          !< Initial misc/optimization variables
   real(ReKi)           ,         intent(inout) :: tau1       !< tau1 value used in DBEMT filter
   integer(IntKi),                intent(  out) :: errStat    !< Error status of the operation
   character(*),                  intent(  out) :: errMsg     !< Error message if ErrStat /= ErrID_None

   ! local variables
   real(ReKi)                                   :: temp

   real(ReKi)                                   :: AxInd_disk
   real(ReKi), parameter                        :: max_AxInd = 0.5_ReKi

   real(ReKi)                                   :: Un_disk
   real(ReKi), parameter                        :: min_Un = 0.1_ReKi
   
   character(*), parameter                      :: RoutineName = 'ComputeTau1'


   ErrStat = ErrID_None
   ErrMsg  = ""

   if ( p%DBEMT_Mod == DBEMT_tauConst ) then
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
            call setErrStat( ErrID_Severe, 'Induced axial relative air speed, Un, is less than '//trim(num2lstr(min_Un)) &
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
      y_vind = u%vind_s
   else     
      y_vind(:) = x%vind(:,i,j)
   end if
   
      
end subroutine DBEMT_CalcOutput

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