!**********************************************************************************************************************************
! The Conv_Radiation and Conv_Radiation_Types modules make up a template for creating user-defined calculations in the FAST Modularization
! Framework. Conv_Radiations_Types will be auto-generated based on a description of the variables for the module.
!
! "Conv_Radiation" should be replaced with the name of your module. Example: HydroDyn
! "Conv_Rdtn" (in Conv_Rdtn_*) should be replaced with the module name or an abbreviation of it. Example: HD
!..................................................................................................................................
! LICENSING
! Copyright (C) 2012-2015  National Renewable Energy Laboratory
!
!    This file is part of Conv_Radiation.
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
MODULE Conv_Radiation

   USE Conv_Radiation_Types
   USE NWTC_Library
   USE NWTC_FFTPACK

   IMPLICIT NONE

   PRIVATE

   REAL(DbKi), PARAMETER, PRIVATE       :: OnePlusEps  = 1.0 + EPSILON(OnePlusEps)   !< The number slighty greater than unity in the precision of DbKi.


   TYPE(ProgDesc), PARAMETER            :: Conv_Rdtn_ProgDesc = ProgDesc( 'Conv_Radiation', '', '' )


      ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: Conv_Rdtn_Init                           ! Initialization routine
   PUBLIC :: Conv_Rdtn_End                            ! Ending routine (includes clean up)

   PUBLIC :: Conv_Rdtn_UpdateStates                   ! Loose coupling routine for solving for constraint states, integrating
                                                      !   continuous states, and updating discrete states
   PUBLIC :: Conv_Rdtn_CalcOutput                     ! Routine for computing outputs

   PUBLIC :: Conv_Rdtn_CalcConstrStateResidual        ! Tight coupling routine for returning the constraint state residual
   PUBLIC :: Conv_Rdtn_CalcContStateDeriv             ! Tight coupling routine for computing derivatives of continuous states
   PUBLIC :: Conv_Rdtn_UpdateDiscState                ! Tight coupling routine for updating discrete states


CONTAINS



!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is called at the start of the simulation to perform initialization steps.
!! The parameters are set here and not changed during the simulation.
!! The initial states and initial guess for the input are defined.
SUBROUTINE Conv_Rdtn_Init( InitInp, u, p, x, xd, z, OtherState, y, m, InitOut, ErrStat, ErrMsg )
!..................................................................................................................................

      TYPE(Conv_Rdtn_InitInputType),       INTENT(IN   )  :: InitInp     !< Input data for initialization routine
      TYPE(Conv_Rdtn_InputType),           INTENT(  OUT)  :: u           !< An initial guess for the input; input mesh must be defined
      TYPE(Conv_Rdtn_ParameterType),       INTENT(  OUT)  :: p           !< Parameters
      TYPE(Conv_Rdtn_ContinuousStateType), INTENT(  OUT)  :: x           !< Initial continuous states
      TYPE(Conv_Rdtn_DiscreteStateType),   INTENT(  OUT)  :: xd          !< Initial discrete states
      TYPE(Conv_Rdtn_ConstraintStateType), INTENT(  OUT)  :: z           !< Initial guess of the constraint states
      TYPE(Conv_Rdtn_OtherStateType),      INTENT(  OUT)  :: OtherState  !< Initial other states
      TYPE(Conv_Rdtn_OutputType),          INTENT(  OUT)  :: y           !< Initial system outputs (outputs are not calculated;
                                                                         !!   only the output mesh is initialized)
      TYPE(Conv_Rdtn_MiscVarType),         INTENT(  OUT)  :: m           !< Initial misc/optimization variables
      TYPE(Conv_Rdtn_InitOutputType),      INTENT(  OUT)  :: InitOut     !< Output for initialization routine
      INTEGER(IntKi),                      INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                        INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None


         ! Local variables

      REAL(SiKi)                             :: Omega                                ! Wave frequency (rad/s)
      REAL(DbKi)                             :: Krnl_Fact                            ! Factor used to scale the magnitude of the RdtnKnrl  as required by the discrete time (co)sine transform (-)
      REAL(DbKi)                             :: RdtnTMax                             ! Analysis time for wave radiation kernel calculations (sec), may be different from Init_Data%RdtnTMax
      REAL(ReKi)                             :: RdtnDOmega                           ! Frequency step for wave radiation kernel calculations (rad/s)
      REAL(ReKi)                             :: RdtnOmegaMax                         ! Maximum frequency used in the (co)sine transform to fine the radiation impulse response functions (rad/s)
      REAL(DbKi), ALLOCATABLE                :: RdtnTime  (:)                        ! Simulation times at which the instantaneous values of the wave radiation kernel are determined (sec)
      LOGICAL                                :: RdtnFrmAM                            ! Determine the wave radiation kernel from the frequency-dependent hydrodynamic added mass matrix? (.TRUE = yes, .FALSE. = determine the wave radiation kernel from the frequency-dependent hydrodynamic damping matrix) !JASON: SHOULD YOU MAKE THIS AN INPUT???<--JASON: IT IS NOT WISE TO COMPUTE THE RADIATION KERNEL FROM THE FREQUENCY-DEPENDENT ADDED MASS MATRIX, UNLESS A CORRECTION IS APPLIED.  THIS IS DESCRIBED IN THE WAMIT USER'S GUIDE!!!!
      INTEGER                                :: NStepRdtn2                           ! ( NStepRdtn-1 )/2
      INTEGER                                :: I                                    ! Generic index
      INTEGER                                :: J                                    ! Generic index
      INTEGER                                :: K                                    ! Generic index
      INTEGER                                :: LastInd                              ! Index into the arrays saved from the last call as a starting point for this call

      TYPE(FFT_DataType)                     :: FFT_Data                             ! the instance of the FFT module we're using

         ! Error handling
      CHARACTER(1024)                        :: ErrMsg2                              ! Temporary error message for calls
      INTEGER(IntKi)                         :: ErrStat2                             ! Temporary error status for calls

         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""

      ! For now, this is the only model we have implemented
      RdtnFrmAM = .FALSE.


         !    If HighFreq is greater than
         !   RdtnOmegaMax, Abort because RdtnDT must be reduced in order to have
         !   sufficient accuracy in the computation of the radiation impulse response
         !   functions:
      p%NBody      = InitInp%NBody
      p%RdtnDT     = InitInp%RdtnDT ! this is also Interval
      RdtnOmegaMax = Pi / InitInp%RdtnDT

      IF ( InitInp%HighFreq > RdtnOmegaMax      )  THEN   ! .TRUE. if the highest frequency component (not counting infinity) in the WAMIT file is greater than RdtnOmegaMax
         ErrMsg =   ' Based on the frequency range found in "'//TRIM(InitInp%WAMITFile)//'.1",'       // &
                      ' RdtnDT must be set smaller than '//TRIM(Num2LStr( Pi/InitInp%HighFreq ))//' sec'// &
                      ' in order to accurately compute the radiation impulse response functions.'
         ErrStat = ErrID_Fatal
         RETURN
      END IF

      call AllocAry( u%Velocity, 6*p%NBody, "u%Velocity"    , ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Conv_Rdtn_Init' )
      call AllocAry( y%F_Rdtn  , 6*p%NBody, "y%F_Rdtn"      , ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Conv_Rdtn_Init' )

      u%Velocity = 0.0 !this is an initial guess;


         ! Perform some initialization computations including calculating the total
         !   number of frequency components = total number of time steps in the wave,
         !   radiation kernel, calculating the frequency step, and ALLOCATing the
         !   arrays:
         ! NOTE: RdtnDOmega = Pi/RdtnTMax since, in the (co)sine transforms:
         !          Omega = (K-1)*RdtnDOmega
         !          Time  = (J-1)*RdtnDT
         !       and therefore:
         !          Omega*Time = (K-1)*(J-1)*RdtnDOmega*RdtnDT
         !                     = (K-1)*(J-1)*Pi/(NStepRdtn-1) [see NWTC_FFTPACK]
         !       or:
         !          RdtnDOmega = Pi/((NStepRdtn-1)*RdtnDT)
         !                     = Pi/RdtnTMax

      p%NStepRdtn  = CEILING ( InitInp%RdtnTMax/p%RdtnDT )                 ! Set NStepRdtn to an odd integer

      IF ( MOD(p%NStepRdtn,2) == 0 )  p%NStepRdtn = p%NStepRdtn + 1  !   larger or equal to RdtnTMax/RdtnDT.

      NStepRdtn2   = MAX( ( p%NStepRdtn-1 )/2, 1 )                                 ! Make sure that NStepRdtn-1 is an even product of small factors (PSF) that is greater
      p%NStepRdtn  = 2*PSF ( NStepRdtn2, 9 ) + 1                                 !   or equal to RdtnTMax/RdtnDT to ensure that the (co)sine transform is efficient.

      p%NStepRdtn1 = p%NStepRdtn + 1                                       ! Save the value of NStepRdtn + 1 for future use.
      !NStepRdtn2   = ( p%NStepRdtn-1 )/2                                           ! Update the value of NStepRdtn2 based on the value needed for NStepRdtn.
      RdtnTMax     = ( p%NStepRdtn-1 )*p%RdtnDT                      ! Update the value of RdtnTMax   based on the value needed for NStepRdtn.
      RdtnDOmega   = Pi/RdtnTMax                                                 ! Compute the frequency step for wave radiation kernel calculations.

      ALLOCATE ( RdtnTime (0:p%NStepRdtn-1    ) , STAT=ErrStat )
      IF ( ErrStat /= ErrID_None )  THEN
         ErrMsg = ' Error allocating memory for the RdtnTime array.'
         ErrStat = ErrID_Fatal
         RETURN
      END IF

      ALLOCATE ( p%RdtnKrnl (0:p%NStepRdtn-1,6*p%NBody,6*p%NBody) , STAT=ErrStat )
      IF ( ErrStat /= ErrID_None )  THEN
         ErrMsg = ' Error allocating memory for the RdtnKrnl array.'
         ErrStat = ErrID_Fatal
         RETURN
      END IF

      ALLOCATE ( xd%XDHistory(0:p%NStepRdtn  ,6*p%NBody  ) , STAT=ErrStat )   ! In the numerical convolution we must have NStepRdtn1 elements within the XDHistory array, which is one more than the NStepRdtn elements that are in the RdtnKrnl array
      IF ( ErrStat /= ErrID_None )  THEN
         ErrMsg = ' Error allocating memory for the XDHistory array.'
         ErrStat = ErrID_Fatal
         RETURN
      END IF

         ! Initialize all elements of the xd%XDHistory array with the intial values of u%Velocity
      DO K = 0,p%NStepRdtn-1
         DO J = 1,6*p%NBody                 ! Loop through all DOFs
            xd%XDHistory(K,J) = u%Velocity(J)
         END DO
      END DO

      LastInd = 1
      IF ( RdtnFrmAM )  THEN  ! .TRUE. if we will determine the wave radiation kernel from the frequency-dependent hydrodynamic added mass matrix



         ! Calculate the factor needed by the discrete sine transform in the
         !   calculation of the wave radiation kernel:

         Krnl_Fact = -1.0_DbKi/p%RdtnDT ! This factor is needed by the discrete time sine transform



         ! Compute all frequency components (including zero) of the sine transform
         !   of the wave radiation kernel:

         DO I = 0,p%NStepRdtn-1 ! Loop through all frequency components (including zero) of the sine transform


         ! Calculate the array of simulation times at which the instantaneous values
         !   of the wave radiation kernel are to be determined:

            RdtnTime(I) = I*p%RdtnDT


         ! Compute the frequency of this component:

            Omega = I*RdtnDOmega


         ! Compute the upper-triangular portion (diagonal and above) of the sine
         !   transform of the wave radiation kernel:

           ! Indx = 0
            DO J = 1,6*p%NBody        ! Loop through all rows    of RdtnKrnl
               DO K = 1,6*p%NBody     ! Loop through all columns of RdtnKrnl above and including the diagonal
                  !Indx = Indx + 1
                  p%RdtnKrnl(I,J,K) = Krnl_Fact*Omega*( InterpStp( Omega, InitInp%HdroFreq(:), &
                                                                                InitInp%HdroAddMs(:       ,J,K), LastInd, InitInp%NInpFreq ) &
                                                      -                         InitInp%HdroAddMs(InitInp%NInpFreq,J,K)                      )
               END DO          ! K - All columns of RdtnKrnl above and including the diagonal
            END DO             ! J - All rows    of RdtnKrnl


         END DO                ! I - All frequency components (including zero) of the sine transform



         ! Compute the sine transforms to find the time-domain representation of
         !   the wave radiation kernel:

         CALL InitSINT ( p%NStepRdtn, FFT_Data, .TRUE., ErrStat )

         IF ( ErrStat /= ErrID_None ) THEN
            ErrMsg  = 'Error Initializating Sine Transforms'
            ErrStat = ErrID_Fatal
            RETURN
         END IF

         DO J = 1,6*p%NBody                 ! Loop through all rows    of RdtnKrnl
            DO K = 1,6*p%NBody              ! Loop through all columns of RdtnKrnl above and including the diagonal
               CALL ApplySINT( p%RdtnKrnl(:,J,K), FFT_Data, ErrStat )
               IF ( ErrStat /= ErrID_None ) RETURN
            END DO                   ! K - All columns of RdtnKrnl above and including the diagonal
         END DO                      ! J - All rows    of RdtnKrnl

         CALL ExitSINT(FFT_Data, ErrStat)
         IF ( ErrStat /= ErrID_None ) THEN
            ErrMsg  = 'Error Cleaning up Sine Transforms'
            ErrStat = ErrID_Fatal
            RETURN
         END IF


      ELSE                    ! We must be determining the wave radiation kernel from the frequency-dependent hydrodynamic damping matrix



         ! Calculate the factor needed by the discrete cosine transform in the
         !   calculation of the wave radiation kernel:

         Krnl_Fact = 1.0/p%RdtnDT  ! This factor is needed by the discrete time cosine transform



         ! Compute all frequency components (including zero) of the cosine transform
         !   of the wave radiation kernel:

         DO I = 0,p%NStepRdtn-1 ! Loop through all frequency components (including zero) of the cosine transform


         ! Calculate the array of simulation times at which the instantaneous values
         !   of the wave radiation kernel are to be determined:

            RdtnTime(I) = I*p%RdtnDT


         ! Compute the frequency of this component:

            Omega = I*RdtnDOmega


         ! Compute the upper-triangular portion (diagonal and above) of the cosine
         !   transform of the wave radiation kernel:

            !Indx = 0
            DO J = 1,6*p%NBody        ! Loop through all rows    of RdtnKrnl
               DO K = 1,6*p%NBody     ! Loop through all columns of RdtnKrnl above and including the diagonal
                  !Indx = Indx + 1
                  p%RdtnKrnl(I,J,K) = Krnl_Fact*InterpStp ( Omega, InitInp%HdroFreq(:), InitInp%HdroDmpng(:,J,K), LastInd, InitInp%NInpFreq )
               END DO          ! K - All columns of RdtnKrnl above and including the diagonal
            END DO             ! J - All rows    of RdtnKrnl


         END DO                ! I - All frequency components (including zero) of the cosine transform



         ! Compute the cosine transforms to find the time-domain representation of
         !   the wave radiation kernel:

         CALL InitCOST ( p%NStepRdtn, FFT_Data, .TRUE., ErrStat )
         IF ( ErrStat /= ErrID_None ) THEN
            ErrMsg  = 'Error Initializating Cosine Transforms'
            ErrStat = ErrID_Fatal
            RETURN
         END IF

         DO J = 1,6*p%NBody                          ! Loop through all rows    of RdtnKrnl
            DO K = 1,6*p%NBody                       ! Loop through all columns of RdtnKrnl above and including the diagonal
               CALL ApplyCOST( p%RdtnKrnl(:,J,K), FFT_Data, ErrStat )
               IF ( ErrStat /= ErrID_None ) THEN
                  ErrMsg  = 'Error applying Cosine Transform'
                  ErrStat = ErrID_Fatal
                  RETURN
               END IF
            END DO                            ! K - All columns of RdtnKrnl above and including the diagonal
         END DO                               ! J - All rows    of RdtnKrnl

         CALL ExitCOST(FFT_Data, ErrStat)
         IF ( ErrStat /= ErrID_None ) THEN
            ErrMsg  = 'Error cleaning up Cosine Transforms'
            ErrStat = ErrID_Fatal
            RETURN
         END IF


      END IF


   IF ( ALLOCATED( RdtnTime     ) ) DEALLOCATE( RdtnTime     )


   m%LastIndRdtn = 0
   OtherState%IndRdtn = 0

      ! bjj: these initializations don't matter, but I don't like seeing the compilation warning in IVF:
   x%DummyContState = 0.0
   z%DummyConstrState = 0.0
   y%F_Rdtn = 0.0
   InitOut%DummyInitOut = 0

END SUBROUTINE Conv_Rdtn_Init
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is called at the end of the simulation.
SUBROUTINE Conv_Rdtn_End( u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )
!..................................................................................................................................

      TYPE(Conv_Rdtn_InputType),           INTENT(INOUT)  :: u           !< System inputs
      TYPE(Conv_Rdtn_ParameterType),       INTENT(INOUT)  :: p           !< Parameters
      TYPE(Conv_Rdtn_ContinuousStateType), INTENT(INOUT)  :: x           !< Continuous states
      TYPE(Conv_Rdtn_DiscreteStateType),   INTENT(INOUT)  :: xd          !< Discrete states
      TYPE(Conv_Rdtn_ConstraintStateType), INTENT(INOUT)  :: z           !< Constraint states
      TYPE(Conv_Rdtn_OtherStateType),      INTENT(INOUT)  :: OtherState  !< Other/optimization states
      TYPE(Conv_Rdtn_OutputType),          INTENT(INOUT)  :: y           !< System outputs
      TYPE(Conv_Rdtn_MiscVarType),         INTENT(INOUT)  :: m           !< Initial misc/optimization variables
      INTEGER(IntKi),                      INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                        INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None



         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""


         ! Place any last minute operations or calculations here:


         ! Close files here:



         ! Destroy the input data:

      CALL Conv_Rdtn_DestroyInput( u, ErrStat, ErrMsg )


         ! Destroy the parameter data:

      CALL Conv_Rdtn_DestroyParam( p, ErrStat, ErrMsg )


         ! Destroy the state data:

      CALL Conv_Rdtn_DestroyContState(   x,           ErrStat, ErrMsg )
      CALL Conv_Rdtn_DestroyDiscState(   xd,          ErrStat, ErrMsg )
      CALL Conv_Rdtn_DestroyConstrState( z,           ErrStat, ErrMsg )
      CALL Conv_Rdtn_DestroyOtherState(  OtherState,  ErrStat, ErrMsg )

      CALL Conv_Rdtn_DestroyMisc( m, ErrStat, ErrMsg )

         ! Destroy the output data:

      CALL Conv_Rdtn_DestroyOutput( y, ErrStat, ErrMsg )




END SUBROUTINE Conv_Rdtn_End

!----------------------------------------------------------------------------------------------------------------------------------
!> Loose coupling routine for solving constraint states, integrating continuous states, and updating discrete states.
!! Continuous, constraint, and discrete states are updated to values at t + Interval.
SUBROUTINE Conv_Rdtn_UpdateStates( t, n, Inputs, InputTimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
!..................................................................................................................................

      REAL(DbKi),                           INTENT(IN   ) :: t               !< Current simulation time in seconds
      INTEGER(IntKi),                       INTENT(IN   ) :: n               !< Current step of the simulation: t = n*Interval
      TYPE(Conv_Rdtn_InputType),            INTENT(INOUT) :: Inputs(:)       !< Inputs at InputTimes
      REAL(DbKi),                           INTENT(IN   ) :: InputTimes(:)   !< Times in seconds associated with Inputs
      TYPE(Conv_Rdtn_ParameterType),        INTENT(IN   ) :: p               !< Parameters
      TYPE(Conv_Rdtn_ContinuousStateType),  INTENT(INOUT) :: x               !< Input: Continuous states at t;
                                                                             !!   Output: Continuous states at t + Interval
      TYPE(Conv_Rdtn_DiscreteStateType),    INTENT(INOUT) :: xd              !< Input: Discrete states at t;
                                                                             !!   Output: Discrete states at t + Interval
      TYPE(Conv_Rdtn_ConstraintStateType),  INTENT(INOUT) :: z               !< Input: Constraint states at t;
                                                                             !!   Output: Constraint states at t + Interval
      TYPE(Conv_Rdtn_OtherStateType),       INTENT(INOUT) :: OtherState      !< Input: Other states at t;
                                                                             !!   Output: Other states at t + Interval
      TYPE(Conv_Rdtn_MiscVarType),          INTENT(INOUT) :: m               !< Initial misc/optimization variables
      INTEGER(IntKi),                       INTENT(  OUT) :: ErrStat         !< Error status of the operation
      CHARACTER(*),                         INTENT(  OUT) :: ErrMsg          !< Error message if ErrStat /= ErrID_None

         ! Local variables

      TYPE(Conv_Rdtn_InputType)                           :: u               !< Instantaneous inputs
      INTEGER(IntKi)                                      :: ErrStat2        !< Error status of the operation (secondary error)
      CHARACTER(ErrMsgLen)                                :: ErrMsg2         !< Error message if ErrStat2 /= ErrID_None
      character(*), parameter                             :: RoutineName = 'Conv_Rdtn_UpdateStates'


         ! Initialize variables

      ErrStat   = ErrID_None           ! no error has occurred
      ErrMsg    = ""


      ! This subroutine contains an example of how the states could be updated. Developers will
      ! want to adjust the logic as necessary for their own situations.


         ! Get the inputs at time t, based on the array of values sent by the glue code:
      call Conv_Rdtn_CopyInput( Inputs(1), u, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
         call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )


      CALL Conv_Rdtn_Input_ExtrapInterp( Inputs, InputTimes, u, t, ErrStat2, ErrMsg2 )
         call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) RETURN


         ! Update discrete states:
         !   Note that xd [discrete state] is changed in Conv_Rdtn_UpdateDiscState() so xd will now contain values at t+Interval
         !   We'll first make a copy that contains xd at time t, which will be used in computing the constraint states

      CALL Conv_Rdtn_UpdateDiscState( t, n, u, p, x, xd, z, OtherState, m, ErrStat2, ErrMsg2 )
         call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )


         ! Integrate (update) continuous states (x) here:

      !x = function of dxdt and x


         ! Destroy local variables before returning

      CALL Conv_Rdtn_DestroyInput(   u, ErrStat2, ErrMsg2)


END SUBROUTINE Conv_Rdtn_UpdateStates

!----------------------------------------------------------------------------------------------------------------------------------
!> Routine for computing outputs, used in both loose and tight coupling.
SUBROUTINE Conv_Rdtn_CalcOutput( Time, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )
!..................................................................................................................................

      REAL(DbKi),                          INTENT(IN   )  :: Time        !< Current simulation time in seconds
      TYPE(Conv_Rdtn_InputType),           INTENT(IN   )  :: u           !< Inputs at Time
      TYPE(Conv_Rdtn_ParameterType),       INTENT(IN   )  :: p           !< Parameters
      TYPE(Conv_Rdtn_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at Time
      TYPE(Conv_Rdtn_DiscreteStateType),   INTENT(IN   )  :: xd          !< Discrete states at Time
      TYPE(Conv_Rdtn_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at Time
      TYPE(Conv_Rdtn_OtherStateType),      INTENT(IN   )  :: OtherState  !< Other states at Time
      TYPE(Conv_Rdtn_OutputType),          INTENT(INOUT)  :: y           !< Outputs computed at Time (Input only so that mesh con-
                                                                         !!   nectivity information does not have to be recalculated)
      TYPE(Conv_Rdtn_MiscVarType),          INTENT(INOUT) :: m           !< Initial misc/optimization variables
      INTEGER(IntKi),                       INTENT(  OUT) :: ErrStat     !< Error status of the operation
      CHARACTER(*),                         INTENT(  OUT) :: ErrMsg      !< Error message if ErrStat /= ErrID_None

!      REAL(ReKi)                           :: F_Rdtn (6)
      REAL(ReKi)                           :: F_RdtnDT (6*p%NBody)                            ! The portion of the total load contribution from wave radiation damping associated with the convolution integral proportional to ( RdtnDT - RdtnRmndr ) (N, N-m)

      INTEGER                              :: I                                       ! Generic index
      INTEGER                              :: J                                       ! Generic index
      INTEGER                              :: K                                       ! Generic index

      INTEGER(IntKi)                       :: MaxInd
         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""


      ! Perform numerical convolution to determine the load contribution from wave
      !   radiation damping:

      MaxInd = MIN(p%NStepRdtn-1,OtherState%IndRdtn)  ! Note: xd%IndRdtn index is from the previous time-step since this state was for the previous time-step

      DO I = 1,6*p%NBody                 ! Loop through all wave radiation damping forces and moments

         F_RdtnDT   (I) = 0.0
       !  F_RdtnRmndr(I) = 0.0

         DO J = 1,6*p%NBody              ! Loop through all platform DOFs
            ! Contribution from the first and last time steps are halved to make the integration 2nd-order accurate
            F_RdtnDT(I) = F_RdtnDT(I) - 0.5_SiKi * p%RdtnKrnl(MaxInd,I,J)*xd%XDHistory(0,J) &
                                      - 0.5_SiKi * p%RdtnKrnl(0,I,J)*xd%XDHistory(MaxInd,J)
            DO K = 1, MaxInd-1 ! Loop through all remaining NStepRdtn-2 time steps in the radiation Kernel (less than NStepRdtn time steps are used when ZTime < RdtnTmax)
               F_RdtnDT(I) = F_RdtnDT(I) - p%RdtnKrnl(MaxInd-K,I,J)*xd%XDHistory(K,J)
            END DO
            !DO K = MAX(0,xd%IndRdtn-p%NStepRdtn  ),xd%IndRdtn-1  ! Loop through all NStepRdtn time steps in the radiation Kernel (less than NStepRdtn time steps are used when ZTime < RdtnTmax)
            !   F_RdtnDT   (I) = F_RdtnDT   (I) - p%RdtnKrnl(xd%IndRdtn-1-K,I,J)*xd%XDHistory(MOD(K,p%NStepRdtn1),J)
            !END DO                                        ! K - All NStepRdtn time steps in the radiation Kernel (less than NStepRdtn time steps are used when ZTime < RdtnTmax)

            !DO K = MAX(0,xd%IndRdtn-p%NStepRdtn+1),xd%IndRdtn    ! Loop through all NStepRdtn time steps in the radiation Kernel (less than NStepRdtn time steps are used when ZTime < RdtnTmax)
            !   F_RdtnRmndr(I) = F_RdtnRmndr(I) - p%RdtnKrnl(xd%IndRdtn  -K,I,J)*xd%XDHistory(MOD(K,p%NStepRdtn1),J)
            !END DO                                        ! K - All NStepRdtn time steps in the radiation Kernel (less than NStepRdtn time steps are used when ZTime < RdtnTmax)

         END DO                   ! J - All platform DOFs

         !F_Rdtn     (I) = ( p%RdtnDT - xd%RdtnRmndr )*F_RdtnDT(I) + xd%RdtnRmndr*F_RdtnRmndr(I)

      END DO                      ! I - All wave radiation damping forces and moments

      y%F_Rdtn = p%RdtnDT*F_RdtnDT !F_Rdtn

END SUBROUTINE Conv_Rdtn_CalcOutput
!----------------------------------------------------------------------------------------------------------------------------------
!> Tight coupling routine for computing derivatives of continuous states.
SUBROUTINE Conv_Rdtn_CalcContStateDeriv( Time, u, p, x, xd, z, OtherState, m, dxdt, ErrStat, ErrMsg )
!..................................................................................................................................

      REAL(DbKi),                           INTENT(IN   )  :: Time        ! Current simulation time in seconds
      TYPE(Conv_Rdtn_InputType),            INTENT(IN   )  :: u           ! Inputs at Time
      TYPE(Conv_Rdtn_ParameterType),        INTENT(IN   )  :: p           ! Parameters
      TYPE(Conv_Rdtn_ContinuousStateType),  INTENT(IN   )  :: x           ! Continuous states at Time
      TYPE(Conv_Rdtn_DiscreteStateType),    INTENT(IN   )  :: xd          ! Discrete states at Time
      TYPE(Conv_Rdtn_ConstraintStateType),  INTENT(IN   )  :: z           ! Constraint states at Time
      TYPE(Conv_Rdtn_OtherStateType),       INTENT(IN   )  :: OtherState  ! Other states at Time
      TYPE(Conv_Rdtn_MiscVarType),          INTENT(INOUT)  :: m           ! Initial misc/optimization variables
      TYPE(Conv_Rdtn_ContinuousStateType),  INTENT(  OUT)  :: dxdt        ! Continuous state derivatives at Time
      INTEGER(IntKi),                       INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                         INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None


         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""


         ! Compute the first time derivatives of the continuous states here:

      dxdt%DummyContState = 0.0


END SUBROUTINE Conv_Rdtn_CalcContStateDeriv

!----------------------------------------------------------------------------------------------------------------------------------
!> Tight coupling routine for updating discrete states.
SUBROUTINE Conv_Rdtn_UpdateDiscState( Time, n, u, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
!..................................................................................................................................

      REAL(DbKi),                          INTENT(IN   )  :: Time        !< Current simulation time in seconds
      INTEGER(IntKi),                      INTENT(IN   )  :: n           !< Current step of the simulation: t = n*Interval
      TYPE(Conv_Rdtn_InputType),           INTENT(IN   )  :: u           !< Inputs at Time
      TYPE(Conv_Rdtn_ParameterType),       INTENT(IN   )  :: p           !< Parameters
      TYPE(Conv_Rdtn_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at Time
      TYPE(Conv_Rdtn_DiscreteStateType),   INTENT(INOUT)  :: xd          !< Input: Discrete states at Time;
                                                                         !! Output: Discrete states at Time + Interval
      TYPE(Conv_Rdtn_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at Time
      TYPE(Conv_Rdtn_OtherStateType),      INTENT(INOUT)  :: OtherState  !< Other states at Time (output: other states at Time + Interval)
                                                                         !! THIS (intent out) BREAKS THE FRAMEWORK (but we don't care at this level)
      TYPE(Conv_Rdtn_MiscVarType),         INTENT(INOUT)  :: m           !< Initial misc/optimization variables
      INTEGER(IntKi),                      INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                        INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

          ! Local Variables
      REAL(ReKi)                           :: RdtnRmndr                  ! Fractional amount of the p%RdtnDT timestep
      INTEGER(IntKi)                       :: J                          ! Generic index
      INTEGER(IntKi)                       :: K                          ! Generic index

         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""


         ! Find the index xd%IndRdtn, where RdtnTime(IndRdtn) is the largest value in
         !   RdtnTime(:) that is less than or equal to Time and find the amount of
         !   time remaining from this calculation:
         ! NOTE: Time is scaled by OnePlusEps to ensure that xd%IndRdtn increments in
         !       steps of 1 when RdtnDT = DT, even in the presence of numerical
         !       precision errors.

      OtherState%IndRdtn   = FLOOR ( ( Time*OnePlusEps )/p%RdtnDT )
      RdtnRmndr    = Time - ( OtherState%IndRdtn*p%RdtnDT ) ! = Time - RdtnTime(IndRdtn); however, RdtnTime(:) has a maximum index of NStepRdtn-1

         ! This subroutine can only be called at integer multiples of p%RdtnDT, if RdtnRmdr > 0, then this requirement has been violated!
      IF (RdtnRmndr > EPSILON(0.0_ReKi) ) THEN
         ErrStat = ErrID_FATAL
         ErrMsg  = "Conv_Rdtn_UpdateDiscState() must be called at integer multiples of the radiation timestep."
         RETURN
      END IF

         ! Save the new values of u%Velocity in xd%XDHistory if:
         !   (1) we are on the initialization pass where Time = 0.0,
         !   (2) we have increased in time by at least RdtnDT, or
         !   (3) the time has not changed since the last time we have increased in
         !       time by at least RdtnDT (i.e., on a call to the corrector)
         !   When saving the new values, interpolate to find all of the values
         !   between index LastIndRdtn and index IndRdtn.  Also, if the XDHistory
         !   array is full, use MOD(Index,NStepRdtn1) to replace the oldest values
         !   with the newest values:
         ! NOTE: When IndRdtn > LastIndRdtn, IndRdtn will equal           LastIndRdtn + 1 if DT <= RdtnDT;
         !       When IndRdtn > LastIndRdtn, IndRdtn will be greater than LastIndRdtn + 1 if DT >  RdtnDT.
   !BJJ: this needs a better check so that it is ALWAYS done (MATLAB/Simulink could possibly avoid this step by starting at Time>0, OR there may be some numerical issues where this is NOT EXACTLY zero)

      IF ( OtherState%IndRdtn < (p%NStepRdtn) )  THEN
         DO J = 1,6*p%NBody  ! Loop through all platform DOFs
            xd%XDHistory(OtherState%IndRdtn,J) = u%Velocity(J)  ! XDHistory was allocated as a zero-based array!
         END DO       ! J - All platform DOFs
      ELSE

         ! Shift the stored history by one index
         DO K = 0,p%NStepRdtn-2
            DO J = 1,6*p%NBody                 ! Loop through all DOFs
               xd%XDHistory(K,J) = xd%XDHistory(K+1,J)
            END DO
         END DO
         DO J = 1,6*p%NBody  ! Loop through all platform DOFs
            xd%XDHistory(p%NStepRdtn-1,J) = u%Velocity(J) ! Set the last array element to the current velocity
         END DO       ! J - All platform DOFs
      END IF

END SUBROUTINE Conv_Rdtn_UpdateDiscState
!----------------------------------------------------------------------------------------------------------------------------------
!> Tight coupling routine for solving for the residual of the constraint state equations.
SUBROUTINE Conv_Rdtn_CalcConstrStateResidual( Time, u, p, x, xd, z, OtherState, m, z_residual, ErrStat, ErrMsg )
!..................................................................................................................................

      REAL(DbKi),                          INTENT(IN   )  :: Time        !< Current simulation time in seconds
      TYPE(Conv_Rdtn_InputType),           INTENT(IN   )  :: u           !< Inputs at Time
      TYPE(Conv_Rdtn_ParameterType),       INTENT(IN   )  :: p           !< Parameters
      TYPE(Conv_Rdtn_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at Time
      TYPE(Conv_Rdtn_DiscreteStateType),   INTENT(IN   )  :: xd          !< Discrete states at Time
      TYPE(Conv_Rdtn_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at Time (possibly a guess)
      TYPE(Conv_Rdtn_OtherStateType),      INTENT(INOUT)  :: OtherState  !< Other states at Time
      TYPE(Conv_Rdtn_MiscVarType),         INTENT(INOUT)  :: m           !< Initial misc/optimization variables
      TYPE(Conv_Rdtn_ConstraintStateType), INTENT(  OUT)  :: z_residual  !< Residual of the constraint state equations using
                                                                         !!     the input values described above
      INTEGER(IntKi),                      INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                        INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None


         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""


         ! Solve for the constraint states here:

      z_residual%DummyConstrState = 0

END SUBROUTINE Conv_Rdtn_CalcConstrStateResidual
!----------------------------------------------------------------------------------------------------------------------------------

END MODULE Conv_Radiation
!**********************************************************************************************************************************
