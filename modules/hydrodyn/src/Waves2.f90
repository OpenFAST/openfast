!**********************************************************************************************************************************
!  NOTE: documentation in this file is written for use with Doxygen 1.8.6 and higher.
!
!> Waves2 module
!!
!!  This module calculates the second order wave forces on a structure. This module is used with HydroDyn in FAST.
!!
!!  This software is written in the FAST modular framework.
!!
!..................................................................................................................................
! LICENSING
! Copyright (C) 2012-2015  National Renewable Energy Laboratory
!
!    This file is part of the Waves2 sub-module of HydroDyn.
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
!
!**********************************************************************************************************************************
MODULE Waves2



         !> Known issues:
         !!

   USE Waves2_Types
!   USE WAMIT_Interp
   USE Waves2_Output
   USE NWTC_Library
   USE NWTC_FFTPACK
   USE Waves,  ONLY : WaveNumber

   IMPLICIT NONE

   PRIVATE

!   INTEGER(IntKi), PARAMETER                             :: DataFormatID = 1  !< Update this value if the data types change (used in Waves2_Pack)
   TYPE(ProgDesc), PARAMETER                             :: Waves2_ProgDesc = ProgDesc( 'Waves2', '', '' )
                                                                              !< This holds the name of the program, version info, and date.

   REAL(DbKi), PARAMETER, PRIVATE                        :: OnePlusEps  = 1.0 + EPSILON(OnePlusEps)   ! The number slighty greater than unity in the precision of DbKi.


      ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: Waves2_Init                           !< Initialization routine
   PUBLIC :: Waves2_End                            !< Ending routine (includes clean up)

   PUBLIC :: Waves2_UpdateStates                   !< Loose coupling routine for solving for constraint states, integrating
                                                   !!   continuous states, and updating discrete states
   PUBLIC :: Waves2_CalcOutput                     !< Routine for computing outputs

   PUBLIC :: Waves2_CalcConstrStateResidual        !< Tight coupling routine for returning the constraint state residual
   PUBLIC :: Waves2_CalcContStateDeriv             !< Tight coupling routine for computing derivatives of continuous states
   PUBLIC :: Waves2_UpdateDiscState                !< Tight coupling routine for updating discrete states


CONTAINS
!----------------------------------------------------------------------------------------------------------------------------------
!> @brief
!!    This routine is called at the start of the simulation to perform initialization steps.
!!    The parameters that are set here are not changed during the simulation.
!!    The initial states and initial guess for the input are defined.
SUBROUTINE Waves2_Init( InitInp, u, p, x, xd, z, OtherState, y, misc, Interval, InitOut, ErrStat, ErrMsg )
!..................................................................................................................................

      TYPE(Waves2_InitInputType),         INTENT(IN   )  :: InitInp              !< Input data for initialization routine
      TYPE(Waves2_InputType),             INTENT(  OUT)  :: u                    !< An initial guess for the input; input mesh must be defined
      TYPE(Waves2_ParameterType),         INTENT(  OUT)  :: p                    !< Parameters
      TYPE(Waves2_ContinuousStateType),   INTENT(  OUT)  :: x                    !< Initial continuous states
      TYPE(Waves2_DiscreteStateType),     INTENT(  OUT)  :: xd                   !< Initial discrete states
      TYPE(Waves2_ConstraintStateType),   INTENT(  OUT)  :: z                    !< Initial guess of the constraint states
      TYPE(Waves2_OtherStateType),        INTENT(  OUT)  :: OtherState           !< Initial other states
      TYPE(Waves2_OutputType),            INTENT(  OUT)  :: y                    !< Initial system outputs (outputs are not calculated; only the output mesh is initialized)
      TYPE(Waves2_MiscVarType),           INTENT(  OUT)  :: misc                 !< Misc/optimization variables
      REAL(DbKi),                         INTENT(INOUT)  :: Interval             !< Coupling interval in seconds: don't change it from the glue code provided value.
      TYPE(Waves2_InitOutputType),        INTENT(  OUT)  :: InitOut              !< Output for initialization routine
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat              !< Error status of the operation
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg               !< Error message if ErrStat /= ErrID_None


         ! Local Variables
      COMPLEX(SiKi)                                      :: ImagNmbr = (0.0,1.0) !< The imaginary number, \f$ \sqrt{-1.0} \f$

      INTEGER(IntKi)                                     :: I                    !< Generic counter
      INTEGER(IntKi)                                     :: J                    !< Generic counter
      INTEGER(IntKi)                                     :: n                    !< Generic counter for calculations
      INTEGER(IntKi)                                     :: m                    !< Generic counter for calculations
      INTEGER(IntKi)                                     :: mu_minus             !< Generic counter for difference kinematics calculations
      INTEGER(IntKi)                                     :: mu_plus              !< Generic counter for sum        kinematics calculations

      REAL(SiKi)                                         :: B_minus              !< The value of the \f$ B^-_{nm} \f$ transfer function for the current n,m,z.
      REAL(SiKi)                                         :: B_plus               !< The value of the \f$ B^+_{nm} \f$ transfer function for the current n,m,z.

      REAL(SiKi)                                         :: Omega_n              !< The frequency corresponding to index n
      REAL(SiKi)                                         :: Omega_m              !< The frequency corresponding to index m
      REAL(SiKi)                                         :: Omega_minus          !< The difference frequency corresponding to \f$ \omega_{\mu^-} \f$
      REAL(SiKi)                                         :: Omega_plus           !< The sum frequency corresponding to \f$ \omega_{\mu^+} \f$

      REAL(SiKi)                                         :: k_n                  !< The wavenumber corresponding to \f$ \omega_n \f$
      REAL(SiKi)                                         :: k_m                  !< The wavenumber corresponding to \f$ \omega_m \f$
      REAL(SiKi)                                         :: k_nm                 !< Value of \f$ k_{nm}^{-} \f$

      COMPLEX(SiKi)                                      :: WaveElevxyPrime0     !< The dot product of the wave vector differece and location \f$ \exp[ -i * (\vec{k_n}-\vec{k_m})\cdot\vec{x}] \f$

      COMPLEX(SiKi)                                      :: WaveElevC_n          !< The complex wave elevation for the nth frequency component
      COMPLEX(SiKi)                                      :: WaveElevC_m          !< The complex wave elevation for the mth frequency component
      COMPLEX(SiKi),                         ALLOCATABLE :: WaveElevC0Norm(:)    !< The complex wave amplitude, normalized for these equations.

         ! Velocity calculations
      REAL(SiKi)                                         :: Ux_nm_minus          !< The value of \f$ _xU^-_{nm}      = B_{nm}^- \cdot (|\vec{k_n}|\cos \theta_n - |\vec{k_m}|\cos \theta_m) \f$ used in calculating the x-component of the second order wave velocity
      REAL(SiKi)                                         :: Ux_nm_plus           !< The value of \f$ _xU^+_{nm}      = B_{nm}^+ \cdot (|\vec{k_n}|\cos \theta_n + |\vec{k_m}|\cos \theta_m) \f$ used in calculating the x-component of the second order sum-frequency wave velocity
      REAL(SiKi)                                         :: Uy_nm_minus          !< The value of \f$ _yU^-_{nm}      = B_{nm}^- \cdot (|\vec{k_n}|\sin \theta_n - |\vec{k_m}|\sin \theta_m) \f$ used in calculating the y-component of the second order wave velocity
      REAL(SiKi)                                         :: Uy_nm_plus           !< The value of \f$ _yU^+_{nm}      = B_{nm}^+ \cdot (|\vec{k_n}|\sin \theta_n + |\vec{k_m}|\sin \theta_m) \f$ used in calculating the y-component of the second order sum-frequency wave velocity
      COMPLEX(SiKi)                                      :: Uz_nm_minus          !< The value of \f$ _z{U}^-_{nm}    = (\imath) \cdot {B}_{nm}^- \cdot k^-_{nm} \cdot \tanh[k^-_{nm}(h+z)]   \f$ used in calculating the z-component of the second order wave velocity
      COMPLEX(SiKi)                                      :: Uz_nm_plus           !< The value of \f$ _z{U}^+_{nm}    = (\imath) \cdot {B}_{nm}^+ \cdot k^+_{nm} \cdot \tanh[k^+_{nm}(h+z)]   \f$ used in calculating the z-component of the second order sum-frequency wave velocity

         ! Acceleration calculations
      COMPLEX(SiKi)                                      :: Accx_nm_minus        !< The value of \f$ _xAcc^-_{nm}    = (\imath) \cdot _xU^-_{nm} \omega_{\mu^-} \f$
      COMPLEX(SiKi)                                      :: Accx_nm_plus         !< The value of \f$ _xAcc^+_{nm}    = (\imath) \cdot _xU^+_{nm} \omega_{\mu^+} \f$
      COMPLEX(SiKi)                                      :: Accy_nm_minus        !< The value of \f$ _yAcc^-_{nm}    = (\imath) \cdot _yU^-_{nm} \omega_{\mu^-} \f$
      COMPLEX(SiKi)                                      :: Accy_nm_plus         !< The value of \f$ _yAcc^+_{nm}    = (\imath) \cdot _yU^+_{nm} \omega_{\mu^+} \f$
      COMPLEX(SiKi)                                      :: Accz_nm_minus        !< The value of \f$ _z{Acc}^-_{nm}  = (\imath) \cdot _zU^-_{nm} \omega_{\mu^-} \f$
      COMPLEX(SiKi)                                      :: Accz_nm_plus         !< The value of \f$ _z{Acc}^+_{nm}  = (\imath) \cdot _zU^+_{nm} \omega_{\mu^+} \f$

         ! Pressure calculations
      REAL(SiKi)                                         :: DynP_nm_minus        !< The value of \f$ \rho_\mathrm{w} B_{nm}^- \omega_{\mu^-} \f$
      REAL(SiKi)                                         :: DynP_nm_plus         !< The value of \f$ \rho_\mathrm{w} B_{nm}^+ \omega_{\mu^+} \f$

         ! Tracking of joints for which we are doing calculations
      REAL(SiKi),                            ALLOCATABLE :: WaveKinzi0Prime(:)   !< zi-coordinates for points where the incident wave kinematics will be computed before applying stretching; these are relative to the mean see level (meters)
      INTEGER(IntKi),                        ALLOCATABLE :: WaveKinPrimeMap(:)   !< Mapping function for the wave kinematics to calculate (based on depth)
      INTEGER(IntKi)                                     :: NWaveKin0Prime       !< Number of points where the incident wave kinematics will be computed before applying stretching to the instantaneous free surface (-)


         ! Second order wave elevation calculations
      REAL(SiKi),                            ALLOCATABLE :: TmpTimeSeries(:)     !< Temporary storage for a wave elevation time series for a single point.
      REAL(SiKi),                            ALLOCATABLE :: TmpTimeSeries2(:)    !< Temporary storage for a wave elevation time series for a single point.
      COMPLEX(SiKi),                         ALLOCATABLE :: TmpFreqSeries(:)     !< Temporary storage for a wave elevation frequency series for a single point.
      COMPLEX(SiKi),                         ALLOCATABLE :: TmpFreqSeries2(:)    !< Temporary storage for a wave elevation frequency series for a single point.


         ! Calculation of 2nd order particle acceleration, velocity, and pressure terms
      COMPLEX(SiKi),                         ALLOCATABLE :: WaveVel2xCDiff(:)    !< Frequency space difference frequency particle velocity     term in the x direction
      COMPLEX(SiKi),                         ALLOCATABLE :: WaveVel2yCDiff(:)    !< Frequency space difference frequency particle velocity     term in the y direction
      COMPLEX(SiKi),                         ALLOCATABLE :: WaveVel2zCDiff(:)    !< Frequency space difference frequency particle velocity     term in the z direction
      COMPLEX(SiKi),                         ALLOCATABLE :: WaveAcc2xCDiff(:)    !< Frequency space difference frequency particle acceleration term in the x direction
      COMPLEX(SiKi),                         ALLOCATABLE :: WaveAcc2yCDiff(:)    !< Frequency space difference frequency particle acceleration term in the y direction
      COMPLEX(SiKi),                         ALLOCATABLE :: WaveAcc2zCDiff(:)    !< Frequency space difference frequency particle acceleration term in the z direction
      REAL(SiKi),                            ALLOCATABLE :: WaveVel2xDiff(:)     !< Time domain     difference frequency particle velocity     term in the x direction
      REAL(SiKi),                            ALLOCATABLE :: WaveVel2yDiff(:)     !< Time domain     difference frequency particle velocity     term in the y direction
      REAL(SiKi),                            ALLOCATABLE :: WaveVel2zDiff(:)     !< Time domain     difference frequency particle velocity     term in the z direction
      REAL(SiKi),                            ALLOCATABLE :: WaveAcc2xDiff(:)     !< Time domain     difference frequency particle acceleration term in the x direction
      REAL(SiKi),                            ALLOCATABLE :: WaveAcc2yDiff(:)     !< Time domain     difference frequency particle acceleration term in the y direction
      REAL(SiKi),                            ALLOCATABLE :: WaveAcc2zDiff(:)     !< Time domain     difference frequency particle acceleration term in the z direction
      COMPLEX(SiKi),                         ALLOCATABLE :: WaveDynP2CDiff(:)    !< Frequency space difference frequency dynamic pressure term
      REAL(SiKi),                            ALLOCATABLE :: WaveDynP2Diff(:)     !< Time domain     difference frequency dynamic pressure term

      COMPLEX(SiKi),                         ALLOCATABLE :: WaveVel2xCSumT1(:)   !< Frequency space difference frequency particle velocity     term 1 in the x direction
      COMPLEX(SiKi),                         ALLOCATABLE :: WaveVel2yCSumT1(:)   !< Frequency space difference frequency particle velocity     term 1 in the y direction
      COMPLEX(SiKi),                         ALLOCATABLE :: WaveVel2zCSumT1(:)   !< Frequency space difference frequency particle velocity     term 1 in the z direction
      COMPLEX(SiKi),                         ALLOCATABLE :: WaveAcc2xCSumT1(:)   !< Frequency space difference frequency particle acceleration term 1 in the x direction
      COMPLEX(SiKi),                         ALLOCATABLE :: WaveAcc2yCSumT1(:)   !< Frequency space difference frequency particle acceleration term 1 in the y direction
      COMPLEX(SiKi),                         ALLOCATABLE :: WaveAcc2zCSumT1(:)   !< Frequency space difference frequency particle acceleration term 1 in the z direction
      REAL(SiKi),                            ALLOCATABLE :: WaveVel2xSumT1(:)    !< Time domain     difference frequency particle velocity     term 1 in the x direction
      REAL(SiKi),                            ALLOCATABLE :: WaveVel2ySumT1(:)    !< Time domain     difference frequency particle velocity     term 1 in the y direction
      REAL(SiKi),                            ALLOCATABLE :: WaveVel2zSumT1(:)    !< Time domain     difference frequency particle velocity     term 1 in the z direction
      REAL(SiKi),                            ALLOCATABLE :: WaveAcc2xSumT1(:)    !< Time domain     difference frequency particle acceleration term 1 in the x direction
      REAL(SiKi),                            ALLOCATABLE :: WaveAcc2ySumT1(:)    !< Time domain     difference frequency particle acceleration term 1 in the y direction
      REAL(SiKi),                            ALLOCATABLE :: WaveAcc2zSumT1(:)    !< Time domain     difference frequency particle acceleration term 1 in the z direction
      COMPLEX(SiKi),                         ALLOCATABLE :: WaveDynP2CSumT1(:)   !< Frequency space difference frequency dynamic pressure term 1
      REAL(SiKi),                            ALLOCATABLE :: WaveDynP2SumT1(:)    !< Time domain     difference frequency dynamic pressure term 1

      COMPLEX(SiKi),                         ALLOCATABLE :: WaveVel2xCSumT2(:)   !< Frequency space difference frequency particle velocity     term 2 in the x direction
      COMPLEX(SiKi),                         ALLOCATABLE :: WaveVel2yCSumT2(:)   !< Frequency space difference frequency particle velocity     term 2 in the y direction
      COMPLEX(SiKi),                         ALLOCATABLE :: WaveVel2zCSumT2(:)   !< Frequency space difference frequency particle velocity     term 2 in the z direction
      COMPLEX(SiKi),                         ALLOCATABLE :: WaveAcc2xCSumT2(:)   !< Frequency space difference frequency particle acceleration term 2 in the x direction
      COMPLEX(SiKi),                         ALLOCATABLE :: WaveAcc2yCSumT2(:)   !< Frequency space difference frequency particle acceleration term 2 in the y direction
      COMPLEX(SiKi),                         ALLOCATABLE :: WaveAcc2zCSumT2(:)   !< Frequency space difference frequency particle acceleration term 2 in the z direction
      REAL(SiKi),                            ALLOCATABLE :: WaveVel2xSumT2(:)    !< Time domain     difference frequency particle velocity     term 2 in the x direction
      REAL(SiKi),                            ALLOCATABLE :: WaveVel2ySumT2(:)    !< Time domain     difference frequency particle velocity     term 2 in the y direction
      REAL(SiKi),                            ALLOCATABLE :: WaveVel2zSumT2(:)    !< Time domain     difference frequency particle velocity     term 2 in the z direction
      REAL(SiKi),                            ALLOCATABLE :: WaveAcc2xSumT2(:)    !< Time domain     difference frequency particle acceleration term 2 in the x direction
      REAL(SiKi),                            ALLOCATABLE :: WaveAcc2ySumT2(:)    !< Time domain     difference frequency particle acceleration term 2 in the y direction
      REAL(SiKi),                            ALLOCATABLE :: WaveAcc2zSumT2(:)    !< Time domain     difference frequency particle acceleration term 2 in the z direction
      COMPLEX(SiKi),                         ALLOCATABLE :: WaveDynP2CSumT2(:)   !< Frequency space difference frequency dynamic pressure term 2
      REAL(SiKi),                            ALLOCATABLE :: WaveDynP2SumT2(:)    !< Time domain     difference frequency dynamic pressure term 2

         ! Stuff for the FFT calculations
      TYPE(FFT_DataType)                                 :: FFT_Data             !< the instance of the FFT module we're using




         ! Temporary error trapping variables
      INTEGER(IntKi)                                     :: ErrStatTmp           !< Temporary variable for holding the error status  returned from a CALL statement
      CHARACTER(2048)                                    :: ErrMsgTmp            !< Temporary variable for holding the error message returned from a CALL statement


         ! Subroutine contents

         ! Initialize Error handling variables

      ErrStat     = ErrID_None
      ErrStatTmp  = ErrID_None
      ErrMsg      = ""
      ErrMsgTmp   = ""

         ! Initialize the data storage
      misc%LastIndWave = 1_IntKi

         ! Initialize the NWTC Subroutine Library and display the information about this module.

      CALL NWTC_Init( )


      !-----------------------------------------------------------------------------
      !> Before attempting to do any real calculations, we first check what was
      !! passed in through _InitInp_ to make sure it makes sense.  That routine will
      !! then copy over the relevant information that should be kept in parameters
      !! (_p_).
      !!
      !! _InitInp_ will also check the flags, existence of files, and set flags
      !! accordingly.
      !-----------------------------------------------------------------------------


      !--------------------------------------------------------------------------------
      ! Check the Min and Max frequencies for the full QTF cases
      !  -- these checks are performed based on the DiffQTFF and SumQTFF flags
      !--------------------------------------------------------------------------------

         ! 1. Check that the min / max diff frequencies make sense if using DiffQTF

      IF ( InitInp%WvDiffQTFF .eqv. .TRUE. ) THEN
         IF ( ( InitInp%WvHiCOffD < InitInp%WvLowCOffD ) .OR. ( InitInp%WvLowCOffD < 0.0 ) ) THEN
            CALL SetErrStat( ErrID_Fatal, ' Programming Error in call to Waves2_Init: '//NewLine// &
                  '           WvHiCOffD must be larger than WvLowCOffD. Both must be positive.'// &
                  '              --> This should have been checked by the calling program.', ErrStat, ErrMsg, 'Waves2_Init')
            CALL CleanUp
            RETURN
         END IF
      END IF


         ! 2. Check that the min / max diff frequencies make sense if using SumQTF

      IF ( InitInp%WvSumQTFF .eqv. .TRUE. ) THEN
         IF ( ( InitInp%WvHiCOffS < InitInp%WvLowCOffS ) .OR. ( InitInp%WvLowCOffS < 0.0 ) ) THEN
            CALL SetErrStat( ErrID_Fatal, ' Programming Error in call to Waves2_Init: '//NewLine// &
                  '           WvHiCOffS must be larger than WvLowCOffS. Both must be positive.'// &
                  '              --> This should have been checked by the calling program.', ErrStat, ErrMsg, 'Waves2_Init')
            CALL CleanUp
            RETURN
         END IF
      END IF



      !--------------------------------------------------------------------------------
      ! Check the size of arrays that were passed in containing the wave info
      !--------------------------------------------------------------------------------


         ! Check that WaveElevC0 is a 2x(NStepWave2+1) sized array (0 index start)

      IF ( SIZE( InitInp%WaveElevC0, DIM=2 ) /= (InitInp%NStepWave2 + 1) ) THEN    ! Expect a 2x(0:NStepWave2) array
         CALL SetErrStat( ErrID_Fatal, ' Programming error in call to Waves2_Init:'//NewLine// &
               '        --> Expected array for WaveElevC0 to be of size 2x'//TRIM(Num2LStr(InitInp%NStepWave2 + 1))// &
               ' (2x(NStepWave2+1)), but instead received array of size '// &
               TRIM(Num2LStr(SIZE(InitInp%WaveElevC0,1)))//'x'//TRIM(Num2LStr(SIZE(InitInp%WaveElevC0,2)))//'.', &
               ErrStat, ErrMsg, 'Waves2_Init')
         CALL CleanUp
         RETURN
      END IF


         ! Check that WaveTime is of size (NStepWave+1)

      IF ( SIZE( InitInp%WaveTime ) /= (InitInp%NStepWave + 1) ) THEN    ! Expect a 2x(0:NStepWave2) array
         CALL SetErrStat( ErrID_Fatal, ' Programming error in call to Waves2_Init:'//NewLine// &
               '        --> Expected array for WaveTime to be of size '//TRIM(Num2LStr(InitInp%NStepWave + 1))// &
               ' (NStepWave+1), but instead received array of size '// &
               TRIM(Num2LStr(SIZE(InitInp%WaveTime)))//'.', &
               ErrStat, ErrMsg, 'Waves2_Init')
         CALL CleanUp
         RETURN
      END IF


      !--------------------------------------------------------------------------------
      ! Now copy over things to parameters...
      !--------------------------------------------------------------------------------

         ! Wave information we need to keep

      p%NWaveElev    = InitInp%NWaveElev
      p%NStepWave    = InitInp%NStepWave
      p%NStepWave2   = InitInp%NStepWave2


         ! Time related information

      p%DT                    =  Interval                   ! Timestep from calling program


         ! Allocate array for the WaveTime information -- array of times to generate output for.  NOTE: can't use MOVE_ALLOC since InitInp is intent in.
      CALL AllocAry( p%WaveTime, SIZE(InitInp%WaveTime,DIM=1), 'array to hold WaveTime', ErrStatTmp, ErrMsgTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveTime.',ErrStat,ErrMsg,'Waves2_Init')
      p%WaveTime              =  InitInp%WaveTime

         ! Difference QTF
      p%WvDiffQTFF            =  InitInp%WvDiffQTFF           ! Flag for calculation

         ! Summation QTF
      p%WvSumQTFF             =  InitInp%WvSumQTFF            ! Flag for calculation


         ! Initialize the channel outputs
      p%NumOuts               =  InitInp%NumOuts
      p%NumOutAll             =  InitInp%NumOutAll

      CALL Wvs2OUT_Init( InitInp, y, p, InitOut, ErrStatTmp, ErrMsgTmp )
      CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, 'Waves2_Init')
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL CleanUp
         RETURN
      END IF



         ! The wave elevation information in frequency space -- we need to normalize this by NStepWave2
      ALLOCATE ( WaveElevC0Norm(0:InitInp%NStepWave2) , STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveElevC0Norm.',ErrStat,ErrMsg,'Waves2_Init')

      DO I=0,InitInp%NStepWave2
         WaveElevC0Norm(I) = CMPLX( InitInp%WaveElevC0(1,I), InitInp%WaveElevC0(2,I), SiKi ) / REAL(InitInp%NStepWave2,SiKi)
      ENDDO

      !--------------------------------------------------------------------------------
      ! Setup WaveKin0Prime -- points from the mesh that are passed in
      !--------------------------------------------------------------------------------

      !> @note Wave stretching will need to be incorporated here when we add it to
      !!       the waves module.

      ! Determine the number of, NWaveKin0Prime, and the zi-coordinates for,
      !   WaveKinzi0Prime(:), points where the incident wave kinematics will be
      !   computed before applying stretching to the instantaneous free surface.
      !   The locations are relative to the mean see level.  These depend on
      !   which incident wave kinematics stretching method is being used:


     ! SELECT CASE ( InitInp%WaveStMod )  ! Which model are we using to extrapolate the incident wave kinematics to the instantaneous free surface?

     ! CASE ( 0 )                 ! None=no stretching.


         ! Since we have no stretching, NWaveKin0Prime and WaveKinzi0Prime(:) are
         !   equal to the number of, and the zi-coordinates for, the points in the
         !   WaveKinzi(:) array between, and including, -WtrDpth and 0.0.

         ! Determine NWaveKin0Prime here:

         NWaveKin0Prime = 0
         DO J = 1,InitInp%NWaveKin   ! Loop through all mesh points  where the incident wave kinematics will be computed
               ! NOTE: We test to 0 instead of MSL2SWL because the locations of WaveKinzi and WtrDpth have already been adjusted using MSL2SWL
            IF (    InitInp%WaveKinzi(J) >= -InitInp%WtrDpth .AND. InitInp%WaveKinzi(J) <= 0 )  THEN
               NWaveKin0Prime = NWaveKin0Prime + 1
            END IF
         END DO                ! J - All Morison nodes where the incident wave kinematics will be computed



         ! ALLOCATE the WaveKinzi0Prime(:) array and compute its elements here:

         ALLOCATE ( WaveKinzi0Prime(NWaveKin0Prime) , STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveKinzi0Prime.',ErrStat,ErrMsg,'Waves2_Init')

         ALLOCATE ( WaveKinPrimeMap(NWaveKin0Prime) , STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveKinPrimeMap.',ErrStat,ErrMsg,'Waves2_Init')

         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF


         I = 1

         DO J = 1,InitInp%NWaveKin ! Loop through all points where the incident wave kinematics will be computed without stretching
               ! NOTE: We test to 0 instead of MSL2SWL because the locations of WaveKinzi and WtrDpth have already been adjusted using MSL2SWL
            IF (    InitInp%WaveKinzi(J) >= -InitInp%WtrDpth .AND. InitInp%WaveKinzi(J) <= 0 )  THEN

               WaveKinzi0Prime(I) =  InitInp%WaveKinzi(J)
               WaveKinPrimeMap(I) =  J
               I = I + 1

            END IF

         END DO                   ! J - All points where the incident wave kinematics will be computed without stretching



      !CASE ( 1, 2 )              ! Vertical stretching or extrapolation stretching.
      !   CALL SetErrStat(ErrID_Fatal,' Vertical and extrapolation stretching not supported in second order calculations.',ErrStat,ErrMsg,'Waves2_Init')
      !
      !
      !CASE ( 3 )                 ! Wheeler stretching.
      !   CALL SetErrStat(ErrID_Fatal,' Wheeler stretching not supported in second order calculations.',ErrStat,ErrMsg,'Waves2_Init')
      !
      !CASE DEFAULT
      !   CALL SetErrStat(ErrID_Fatal,' Stretching is not supported in the second order waves kinematics calculations.',ErrStat,ErrMsg,'Waves2_Init')
      !
      !
      !ENDSELECT


      IF ( ErrStat >= AbortErrLev ) THEN
         CALL CleanUp()
         RETURN
      END IF




      !--------------------------------------------------------------------------------
      ! Setup the output arrays
      !--------------------------------------------------------------------------------


      ALLOCATE ( p%WaveElev2 (0:InitInp%NStepWave,InitInp%NWaveElev  ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array p%WaveElev2.', ErrStat,ErrMsg,'Waves2_Init')

      ALLOCATE ( InitOut%WaveVel2D  (0:InitInp%NStepWave,InitInp%NWaveKin,3), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array InitOut%WaveVel2D.',  ErrStat,ErrMsg,'Waves2_Init')

      ALLOCATE ( InitOut%WaveAcc2D  (0:InitInp%NStepWave,InitInp%NWaveKin,3), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array InitOut%WaveAcc2D.',  ErrStat,ErrMsg,'Waves2_Init')

      ALLOCATE ( InitOut%WaveDynP2D (0:InitInp%NStepWave,InitInp%NWaveKin  ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array InitOut%WaveDynP2D.', ErrStat,ErrMsg,'Waves2_Init')

      ALLOCATE ( InitOut%WaveVel2S  (0:InitInp%NStepWave,InitInp%NWaveKin,3), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array InitOut%WaveVel2S.',  ErrStat,ErrMsg,'Waves2_Init')

      ALLOCATE ( InitOut%WaveAcc2S  (0:InitInp%NStepWave,InitInp%NWaveKin,3), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array InitOut%WaveAcc2S.',  ErrStat,ErrMsg,'Waves2_Init')

      ALLOCATE ( InitOut%WaveDynP2S (0:InitInp%NStepWave,InitInp%NWaveKin  ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array InitOut%WaveDynP2S.', ErrStat,ErrMsg,'Waves2_Init')

         ! Now check if all the allocations worked properly
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL CleanUp()
         RETURN
      END IF


         !Initialize the output arrays to zero.  We will only fill it in for the points we calculate.
      p%WaveElev2          =  0.0_SiKi
      InitOut%WaveVel2D    =  0.0_SiKi
      InitOut%WaveAcc2D    =  0.0_SiKi
      InitOut%WaveDynP2D   =  0.0_SiKi
      InitOut%WaveVel2S  =  0.0_SiKi
      InitOut%WaveAcc2S  =  0.0_SiKi
      InitOut%WaveDynP2S =  0.0_SiKi



         ! For creating animations of the sea surface, the WaveElevXY array is passed in with a series of x,y coordinates
         ! (index 1).  The second index corresponds to the number of points passed in.  A two dimensional time series
         ! is created with the first index corresponding to the timestep, and second index corresponding to the second
         ! index of the WaveElevXY array.
      IF ( ALLOCATED(InitInp%WaveElevXY)) THEN
         ALLOCATE ( InitOut%WaveElevSeries2 (0:InitInp%NStepWave, 1:SIZE(InitInp%WaveElevXY, DIM=2)) , STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) THEN
            CALL SetErrStat(ErrID_Fatal,'Cannot allocate array InitOut%WaveElevSeries2.',ErrStat,ErrMsg,'Waves2_Init')
            CALL CleanUp()
            RETURN
         END IF
      ENDIF


         ! For calculating the 2nd-order wave elevation corrections, we need a temporary array to hold the information.
      ALLOCATE ( TmpTimeSeries(0:InitInp%NStepWave), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array TmpTimeSeries.', ErrStat,ErrMsg,'Waves2_Init')
      ALLOCATE ( TmpTimeSeries2(0:InitInp%NStepWave), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array TmpTimeSeries2.', ErrStat,ErrMsg,'Waves2_Init')

      ALLOCATE ( TmpFreqSeries(0:InitInp%NStepWave2), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array TmpFreqSeries.', ErrStat,ErrMsg,'Waves2_Init')
      ALLOCATE ( TmpFreqSeries2(0:InitInp%NStepWave2), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array TmpFreqSeries2.', ErrStat,ErrMsg,'Waves2_Init')



      !--------------------------------------------------------------------------------
      ! Setup the FFT working arrays
      !--------------------------------------------------------------------------------

      CALL InitFFT ( InitInp%NStepWave, FFT_Data, .FALSE., ErrStatTmp )
      CALL SetErrStat(ErrStatTmp,'Error occured while initializing the FFT.',ErrStat,ErrMsg,'Waves2_Init')
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL CleanUp()
         RETURN
      END IF


      !--------------------------------------------------------------------------------
      !> # Difference Frequency #
      !! The second order difference frequency corrections to the velocity, acceleration,
      !! and pressure can be written in generic form as:
      !!
      !! \f$
      !!    V^{(2)-}(t) =  2  \Re \left[ \sum_{\mu^-=1}^{\frac{N}{2}-1} H^-(\omega_{\mu^-})
      !!                      \exp(\imath \omega_{\mu^-} t) \right]
      !!                =  2  \operatorname{IFFT}\left[H^-\right]     \f$
      !!
      !! Notice that in the equations that follow, there is no term analagous to the mean
      !! drift term in the WAMIT2 module.  Rather, for \f$ \omega_{\mu^-} = 0 \f$, the
      !! result is zero.  So this term is not included.
      !--------------------------------------------------------------------------------


      IF(p%WvDiffQTFF) THEN

            ! Tell our nice users what is about to happen that may take a while:
         CALL WrScr ( ' Calculating second order difference frequency wave kinematics.' )


         !--------------------------------------------------------------------------------
         ! Setup arrays for the calculations
         !--------------------------------------------------------------------------------

            ! Frequency space arrays:

         ALLOCATE ( WaveVel2xCDiff   (0:InitInp%NStepWave2), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveVel2xCDiff.',  ErrStat,ErrMsg,'Waves2_Init')
         ALLOCATE ( WaveVel2yCDiff   (0:InitInp%NStepWave2), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveVel2yCDiff.',  ErrStat,ErrMsg,'Waves2_Init')
         ALLOCATE ( WaveVel2zCDiff   (0:InitInp%NStepWave2), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveVel2zCDiff.',  ErrStat,ErrMsg,'Waves2_Init')

         ALLOCATE ( WaveAcc2xCDiff   (0:InitInp%NStepWave2), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveAcc2xCDiff.',  ErrStat,ErrMsg,'Waves2_Init')
         ALLOCATE ( WaveAcc2yCDiff   (0:InitInp%NStepWave2), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveAcc2yCDiff.',  ErrStat,ErrMsg,'Waves2_Init')
         ALLOCATE ( WaveAcc2zCDiff   (0:InitInp%NStepWave2), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveAcc2zCDiff.',  ErrStat,ErrMsg,'Waves2_Init')

         ALLOCATE ( WaveDynP2CDiff   (0:InitInp%NStepWave2), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveDynP2CDiff.',  ErrStat,ErrMsg,'Waves2_Init')

            ! Now check if all the allocations worked properly
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF


            ! Time domain arrays:
         ALLOCATE ( WaveVel2xDiff   (0:InitInp%NStepWave), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveVel2xDiff.',  ErrStat,ErrMsg,'Waves2_Init')
         ALLOCATE ( WaveVel2yDiff   (0:InitInp%NStepWave), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveVel2yDiff.',  ErrStat,ErrMsg,'Waves2_Init')
         ALLOCATE ( WaveVel2zDiff   (0:InitInp%NStepWave), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveVel2zDiff.',  ErrStat,ErrMsg,'Waves2_Init')

         ALLOCATE ( WaveAcc2xDiff   (0:InitInp%NStepWave), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveAcc2xDiff.',  ErrStat,ErrMsg,'Waves2_Init')
         ALLOCATE ( WaveAcc2yDiff   (0:InitInp%NStepWave), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveAcc2yDiff.',  ErrStat,ErrMsg,'Waves2_Init')
         ALLOCATE ( WaveAcc2zDiff   (0:InitInp%NStepWave), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveAcc2zDiff.',  ErrStat,ErrMsg,'Waves2_Init')

         ALLOCATE ( WaveDynP2Diff   (0:InitInp%NStepWave), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveDynP2Diff.',  ErrStat,ErrMsg,'Waves2_Init')

            ! Now check if all the allocations worked properly
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF



            !--------------------------------------------------------------------------------
            !> ## Calculate the surface elevation corrections ##
            !!
            !! For each (x,y) coordinate that a wave elevation is requested at (both from the
            !! (WaveElevxi,WaveElevyi) pairs, and the WaveElevXY pairs), a call is made to the
            !! subroutine waves2::waveelevtimeseriesatxy_diff to calculate the full time series for
            !! that point.  The results are added to the wave elevation results from the sum
            !! frequency calculations later in the code.
            !--------------------------------------------------------------------------------

            ! Step through the requested points
         DO I=1,InitInp%NWaveElev
            CALL WaveElevTimeSeriesAtXY_Diff(InitInp%WaveElevxi(I), InitInp%WaveElevyi(I), TmpTimeSeries, ErrStatTmp, ErrMsgTmp )
            CALL SetErrStat(ErrStatTmp,'Error occured while applying the FFT to InitOut%WaveElev.',ErrStat,ErrMsg,'Waves2_Init')
            IF ( ErrStat >= AbortErrLev ) THEN
               CALL CleanUp()
               RETURN
            END IF
            p%WaveElev2(:,I) = TmpTimeSeries(:)
         ENDDO    ! Wave elevation points requested


            ! Calculate the wave elevation at all points requested in the array WaveElevXY
         IF ( ALLOCATED(InitInp%WaveElevXY) ) THEN
            DO I = 1,SIZE(InitInp%WaveElevXY, DIM=2)
                  ! This subroutine call applies the FFT at the correct location.
               CALL WaveElevTimeSeriesAtXY_Diff( InitInp%WaveElevXY(1,I), InitInp%WaveElevXY(2,I), TmpTimeSeries, ErrStatTmp, ErrMsgTmp )
               CALL SetErrStat(ErrStatTmp,ErrMsgTmp,ErrStat,ErrMsg,'Waves2_Init')
               IF ( ErrStat >= AbortErrLev ) THEN
                  CALL CleanUp()
                  RETURN
               END IF
               InitOut%WaveElevSeries2(:,I) = TmpTimeSeries(:)
            ENDDO
         ENDIF



         !--------------------------------------------------------------------------------
         !> ## Calculate the second order velocity, acceleration, and pressure corrections for all joints below surface. ##
         !--------------------------------------------------------------------------------


            ! NWaveKin0Prime loop start
         DO I=1,NWaveKin0Prime


               ! Reset the \f$ H_{\mu^-} \f$ terms to zero before calculating.
            WaveVel2xCDiff = CMPLX(0.0_SiKi, 0.0_SiKi, SiKi)
            WaveVel2yCDiff = CMPLX(0.0_SiKi, 0.0_SiKi, SiKi)
            WaveVel2zCDiff = CMPLX(0.0_SiKi, 0.0_SiKi, SiKi)
            WaveAcc2xCDiff = CMPLX(0.0_SiKi, 0.0_SiKi, SiKi)
            WaveAcc2yCDiff = CMPLX(0.0_SiKi, 0.0_SiKi, SiKi)
            WaveAcc2zCDiff = CMPLX(0.0_SiKi, 0.0_SiKi, SiKi)
            WaveDynP2CDiff = CMPLX(0.0_SiKi, 0.0_SiKi, SiKi)


               ! \f$ \mu^- \f$ loop.  This loop is used to construct the full set of \f$ H_{\mu^-} \f$ terms used in the IFFT to find the timeseries.
               !> * \f$ \mu^- = n -m \f$
            DO mu_minus=1,InitInp%NStepWave2-1

                  ! The frequency we are dealing with
                  !> * \f$ \omega^- = \mu^- \Delta \omega \f$
               Omega_minus =  mu_minus * InitInp%WaveDOmega

               IF ( Omega_minus >= InitInp%WvLowCOffD .AND. Omega_minus <= InitInp%WvHiCOffD ) THEN

                     ! The inner \f$ m \f$ loop for calculating the \f$ H_{\mu^-} \f$ terms at each frequency.
                  DO m=1,InitInp%NStepWave2-mu_minus
                        ! Calculate the value of the n index from \f$ \mu^- = n - m \f$.  Calculate corresponding wavenumbers and frequencies.
                     n           =  mu_minus + m
                     Omega_n     =  n * InitInp%WaveDOmega
                     Omega_m     =  m * InitInp%WaveDOmega
                     k_n         =  WaveNumber( Omega_n, InitInp%Gravity, InitInp%WtrDpth )
                     k_m         =  WaveNumber( Omega_m, InitInp%Gravity, InitInp%WtrDpth )
                     k_nm        =  k_nm_minus( n, m, k_n, k_m )


                        ! Calculate the terms \f$ n,m \f$ necessary for calculations

                        !> Calculate the dot product of the wavenumbers with the (x,y) location
                        !! This is given by:
                        !!
                        !! *  \f$ \exp\left(-\imath \left[\vec{k_n} - \vec{k_m} \right] \cdot \vec{x} \right)
                        !!       = \exp \left( -\imath \left[
                        !!                \left( |\vec{k_n}| \cos \theta_n - |\vec{k_m}| cos \theta_m \right) ~ x
                        !!             +  \left( |\vec{k_n}| \sin \theta_n - |\vec{k_m}| sin \theta_m \right) ~ y \right] \right) \f$

                     WaveElevxyPrime0  = exp( - ImagNmbr &
                              *  (  ( k_n * COS( D2R_S*InitInp%WaveDirArr(n) ) - k_m * COS( D2R_S*InitInp%WaveDirArr(m) ) ) * InitInp%WaveKinxi(WaveKinPrimeMap(I))  &
                                 +  ( k_n * SIN( D2R_S*InitInp%WaveDirArr(n) ) - k_m * SIN( D2R_S*InitInp%WaveDirArr(m) ) ) * InitInp%WaveKinyi(WaveKinPrimeMap(I))  ))


                        ! Get value for \f$ B^- \f$ for the n,m index pair
                     B_minus  =  TransFuncB_minus( n, m, k_n, k_m, WaveKinzi0Prime(I) )


                        !> Calculate \f$ U^- \f$ terms for the velocity calculations (\f$B^-\f$ provided by waves2::transfuncb_minus)
                        ! NOTE: InitInp%WtrDpth + WaveKinzi0Prime(I) is the height above the ocean floor
                        !> * \f$ _x{U}_{nm}^- = B_{nm}^- \left(k_n \cos \theta_n - k_m \cos \theta_m \right) \f$
                     Ux_nm_minus = B_minus * ( k_n * COS( D2R_S*InitInp%WaveDirArr(n) ) - k_m * COS( D2R_S*InitInp%WaveDirArr(m) ) )

                        !> * \f$ _y{U}_{nm}^- = B_{nm}^- \left(k_n \sin \theta_n - k_m \sin \theta_m \right) \f$
                     Uy_nm_minus = B_minus * ( k_n * SIN( D2R_S*InitInp%WaveDirArr(n) ) - k_m * SIN( D2R_S*InitInp%WaveDirArr(m) ) )

                        !> * \f$ _z{U}_{nm}^- = \imath B_{nm}^- k_{nm} \tanh \left( k_{nm} ( h + z ) \right) \f$
                     Uz_nm_minus = ImagNmbr * B_minus * k_nm * tanh( k_nm * ( InitInp%WtrDpth + WaveKinzi0Prime(I) ) )


                        !> Acceleration calculations
                     Accx_nm_minus = ImagNmbr * Ux_nm_minus * Omega_minus     !> * \f$ _x\dot{U}_{nm}^- = \imath * _xU_{nm}^- * \omega_{\mu^-} \f$
                     Accy_nm_minus = ImagNmbr * Uy_nm_minus * Omega_minus     !> * \f$ _y\dot{U}_{nm}^- = \imath * _yU_{nm}^- * \omega_{\mu^-} \f$
                     Accz_nm_minus = ImagNmbr * Uz_nm_minus * Omega_minus     !> * \f$ _z\dot{U}_{nm}^- = \imath * _zU_{nm}^- * \omega_{\mu^-} \f$


                        !> Dynamic pressure
                        !> * \f$ P_{nm}^- = \rho_\mathrm{w} B_{nm}^- \omega_{\mu^-} \f$
                     DynP_nm_minus  = REAL(InitInp%WtrDens,SiKi) * B_minus * Omega_minus



                        !> ### Calculate the inner summation \f$ H^-(\omega_{\mu^-}) \f$ terms for the velocity, acceleration, and pressure. ###


                        ! First get the wave amplitude -- must be reconstructed from the WaveElevC0 array.  First index is the real (1) or
                        ! imaginary (2) part.  Divide by NStepWave2 to remove the built in normalization in WaveElevC0.  Note that the phase
                        ! shift associated with the (x,y) location is accounted for by the WaveElevxyPrime0 variable.
                     WaveElevC_n =  WaveElevC0Norm(n)
                     WaveElevC_m =  WaveElevC0Norm(m)

                        !> Velocity terms:
                        !!    *  \f$ H^-(\omega_{\mu^-}) =  {\sum_{m=1}^{\frac{N}{2}-\mu^{-}}}  A_n  A^*_m U_{nm}^-
                        !!                                  \exp\left(-\imath (\vec{k_n} - \vec{k_m})\cdot\vec{x}\right) \f$
                     WaveVel2xCDiff(mu_minus)   =  WaveVel2xCDiff(mu_minus) + WaveElevC_n * CONJG( WaveElevC_m ) * Ux_nm_minus * WaveElevxyPrime0
                     WaveVel2yCDiff(mu_minus)   =  WaveVel2yCDiff(mu_minus) + WaveElevC_n * CONJG( WaveElevC_m ) * Uy_nm_minus * WaveElevxyPrime0
                     WaveVel2zCDiff(mu_minus)   =  WaveVel2zCDiff(mu_minus) + WaveElevC_n * CONJG( WaveElevC_m ) * Uz_nm_minus * WaveElevxyPrime0

                        !> Acceleration terms:
                        !!    *  \f$ H^-(\omega_{\mu^-}) =  {\sum_{m=1}^{\frac{N}{2}-\mu^{-}}}  A_n  A^*_m \dot{U}_{nm}^-
                        !!                                  \exp\left(-\imath (\vec{k_n} - \vec{k_m})\cdot\vec{x}\right) \f$
                     WaveAcc2xCDiff(mu_minus)   =  WaveAcc2xCDiff(mu_minus) + WaveElevC_n * CONJG( WaveElevC_m ) * Accx_nm_minus * WaveElevxyPrime0
                     WaveAcc2yCDiff(mu_minus)   =  WaveAcc2yCDiff(mu_minus) + WaveElevC_n * CONJG( WaveElevC_m ) * Accy_nm_minus * WaveElevxyPrime0
                     WaveAcc2zCDiff(mu_minus)   =  WaveAcc2zCDiff(mu_minus) + WaveElevC_n * CONJG( WaveElevC_m ) * Accz_nm_minus * WaveElevxyPrime0

                        !> Pressure term:
                        !!    *  \f$ H^-(\omega_{\mu^-}) =  {\sum_{m=1}^{\frac{N}{2}-\mu^{-}}}  A_n  A^*_m P_{nm}^-
                        !!                                  \exp\left(-\imath (\vec{k_n} - \vec{k_m})\cdot\vec{x}\right) \f$
                     WaveDynP2CDiff(mu_minus)   =  WaveDynP2CDiff(mu_minus) + WaveElevC_n * CONJG( WaveElevC_m ) * DynP_nm_minus * WaveElevxyPrime0


                  ENDDO ! m loop


               ENDIF ! Check to see if WvLowCOffD <= mu_minus <= WvHiCOffD

            ENDDO ! mu_minus loop (diff frequency)

                  !  Divide by two for the single sided FFT given in the documentation.
            WaveVel2xCDiff =  WaveVel2xCDiff / 2.0_SiKi
            WaveVel2yCDiff =  WaveVel2yCDiff / 2.0_SiKi
            WaveVel2zCDiff =  WaveVel2zCDiff / 2.0_SiKi
            WaveAcc2xCDiff =  WaveAcc2xCDiff / 2.0_SiKi
            WaveAcc2yCDiff =  WaveAcc2yCDiff / 2.0_SiKi
            WaveAcc2zCDiff =  WaveAcc2zCDiff / 2.0_SiKi
            WaveDynP2CDiff =  WaveDynP2CDiff / 2.0_SiKi




               !> ### Apply the inverse FFT to each of the components to get the time domain result ###
               !> *   \f$ V(t) = 2 \operatorname{IFFT}\left[H^-\right] \f$
            CALL ApplyFFT_cx(  WaveVel2xDiff(:),  WaveVel2xCDiff(:), FFT_Data, ErrStatTmp )
               CALL SetErrStat(ErrStatTmp,'Error occured while applying the FFT on V_x.',ErrStat,ErrMsg,'Waves2_Init')
            CALL ApplyFFT_cx(  WaveVel2yDiff(:),  WaveVel2yCDiff(:), FFT_Data, ErrStatTmp )
               CALL SetErrStat(ErrStatTmp,'Error occured while applying the FFT on V_y.',ErrStat,ErrMsg,'Waves2_Init')
            CALL ApplyFFT_cx(  WaveVel2zDiff(:),  WaveVel2zCDiff(:), FFT_Data, ErrStatTmp )
               CALL SetErrStat(ErrStatTmp,'Error occured while applying the FFT on V_z.',ErrStat,ErrMsg,'Waves2_Init')

            CALL ApplyFFT_cx(  WaveAcc2xDiff(:),  WaveAcc2xCDiff(:), FFT_Data, ErrStatTmp )
               CALL SetErrStat(ErrStatTmp,'Error occured while applying the FFT on Acc_x.',ErrStat,ErrMsg,'Waves2_Init')
            CALL ApplyFFT_cx(  WaveAcc2yDiff(:),  WaveAcc2yCDiff(:), FFT_Data, ErrStatTmp )
               CALL SetErrStat(ErrStatTmp,'Error occured while applying the FFT on Acc_y.',ErrStat,ErrMsg,'Waves2_Init')
            CALL ApplyFFT_cx(  WaveAcc2zDiff(:),  WaveAcc2zCDiff(:), FFT_Data, ErrStatTmp )
               CALL SetErrStat(ErrStatTmp,'Error occured while applying the FFT on Acc_z.',ErrStat,ErrMsg,'Waves2_Init')

            CALL ApplyFFT_cx(  WaveDynP2Diff(:),  WaveDynP2CDiff(:), FFT_Data, ErrStatTmp )
               CALL SetErrStat(ErrStatTmp,'Error occured while applying the FFT on DynP2.',ErrStat,ErrMsg,'Waves2_Init')



            IF ( ErrStat >= AbortErrLev ) THEN
               CALL CleanUp()
               RETURN
            END IF


               ! Copy the results to the output
            InitOut%WaveVel2D(:,WaveKinPrimeMap(I),1) =  2.0_SiKi * WaveVel2xDiff(:)     ! x-component of velocity
            InitOut%WaveVel2D(:,WaveKinPrimeMap(I),2) =  2.0_SiKi * WaveVel2yDiff(:)     ! y-component of velocity
            InitOut%WaveVel2D(:,WaveKinPrimeMap(I),3) =  2.0_SiKi * WaveVel2zDiff(:)     ! z-component of velocity

            InitOut%WaveAcc2D(:,WaveKinPrimeMap(I),1) =  2.0_SiKi * WaveAcc2xDiff(:)     ! x-component of acceleration
            InitOut%WaveAcc2D(:,WaveKinPrimeMap(I),2) =  2.0_SiKi * WaveAcc2yDiff(:)     ! y-component of acceleration
            InitOut%WaveAcc2D(:,WaveKinPrimeMap(I),3) =  2.0_SiKi * WaveAcc2zDiff(:)     ! z-component of acceleration

            InitOut%WaveDynP2D(:,WaveKinPrimeMap(I))  =  2.0_SiKi * WaveDynP2Diff(:)     ! Dynamic pressure


               ! Copy the first point to the last to make it easier.
            InitOut%WaveVel2D(InitInp%NStepWave,WaveKinPrimeMap(I),1)   =  WaveVel2xDiff(0)
            InitOut%WaveVel2D(InitInp%NStepWave,WaveKinPrimeMap(I),2)   =  WaveVel2yDiff(0)
            InitOut%WaveVel2D(InitInp%NStepWave,WaveKinPrimeMap(I),3)   =  WaveVel2zDiff(0)

            InitOut%WaveAcc2D(InitInp%NStepWave,WaveKinPrimeMap(I),1)   =  WaveAcc2xDiff(0)
            InitOut%WaveAcc2D(InitInp%NStepWave,WaveKinPrimeMap(I),2)   =  WaveAcc2yDiff(0)
            InitOut%WaveAcc2D(InitInp%NStepWave,WaveKinPrimeMap(I),3)   =  WaveAcc2zDiff(0)

            InitOut%WaveDynP2D(InitInp%NStepWave,WaveKinPrimeMap(I))    =  WaveDynP2Diff(0)


         ENDDO    ! I=1,NWaveKin0Prime loop end


            ! Deallocate working arrays.
         IF (ALLOCATED(WaveVel2xCDiff))   DEALLOCATE(WaveVel2xCDiff,    STAT=ErrStatTmp)
         IF (ALLOCATED(WaveVel2yCDiff))   DEALLOCATE(WaveVel2yCDiff,    STAT=ErrStatTmp)
         IF (ALLOCATED(WaveVel2zCDiff))   DEALLOCATE(WaveVel2zCDiff,    STAT=ErrStatTmp)
         IF (ALLOCATED(WaveAcc2xCDiff))   DEALLOCATE(WaveAcc2xCDiff,    STAT=ErrStatTmp)
         IF (ALLOCATED(WaveAcc2yCDiff))   DEALLOCATE(WaveAcc2yCDiff,    STAT=ErrStatTmp)
         IF (ALLOCATED(WaveAcc2zCDiff))   DEALLOCATE(WaveAcc2zCDiff,    STAT=ErrStatTmp)
         IF (ALLOCATED(WaveDynP2CDiff))   DEALLOCATE(WaveDynP2CDiff,    STAT=ErrStatTmp)

         IF (ALLOCATED(WaveVel2xDiff))    DEALLOCATE(WaveVel2xDiff,     STAT=ErrStatTmp)
         IF (ALLOCATED(WaveVel2yDiff))    DEALLOCATE(WaveVel2yDiff,     STAT=ErrStatTmp)
         IF (ALLOCATED(WaveVel2zDiff))    DEALLOCATE(WaveVel2zDiff,     STAT=ErrStatTmp)
         IF (ALLOCATED(WaveAcc2xDiff))    DEALLOCATE(WaveAcc2xDiff,     STAT=ErrStatTmp)
         IF (ALLOCATED(WaveAcc2yDiff))    DEALLOCATE(WaveAcc2yDiff,     STAT=ErrStatTmp)
         IF (ALLOCATED(WaveAcc2zDiff))    DEALLOCATE(WaveAcc2zDiff,     STAT=ErrStatTmp)
         IF (ALLOCATED(WaveDynP2Diff))    DEALLOCATE(WaveDynP2Diff,     STAT=ErrStatTmp)

         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF


      ENDIF    ! p%WvDiffQTFF





      !--------------------------------------------------------------------------------
      !> # Sum Frequency #
      !!
      !! The sum frequency corrections to the velocity, acceleration, and pressure can
      !! be written in generic form as:
      !!
      !! \f{eqnarray*}{
      !!    V^{(2)+}(t) &=&   \Re \left[ \sum_{n=1}^{\lfloor \frac{N}{4} \rfloor} K^+(\omega_n) 
      !!                         \exp(\imath 2 \omega_n t) \right]
      !!                      +  2 \Re \left[ \sum_{\mu^+=2}^{\frac{N}{2}} H^+(\omega_{\mu^+})
      !!                         \exp(\imath \omega_{\mu^+} t) \right]\\
      !!                &=&      \operatorname{IFFT}\left[K(\omega_n)\right]
      !!                      + 2\operatorname{IFFT}\left[H(\omega_{\mu^+})\right]     \f}
      !!
      !! Notice that the first term only contains the \f$ 2 \omega_n \f$ terms; the others
      !! are zero.
      !--------------------------------------------------------------------------------


      IF(p%WvSumQTFF) THEN

            ! Tell our nice users what is about to happen that may take a while:
         CALL WrScr ( ' Calculating second order sum frequency wave kinematics.' )


         !--------------------------------------------------------------------------------
         ! Setup arrays for the calculations
         !--------------------------------------------------------------------------------

            ! Frequency space arrays:  Term 1 (n=m term)

         ALLOCATE ( WaveVel2xCSumT1    (0:InitInp%NStepWave2), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveVel2xCSumT1.',  ErrStat,ErrMsg,'Waves2_Init')
         ALLOCATE ( WaveVel2yCSumT1    (0:InitInp%NStepWave2), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveVel2yCSumT1.',  ErrStat,ErrMsg,'Waves2_Init')
         ALLOCATE ( WaveVel2zCSumT1    (0:InitInp%NStepWave2), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveVel2zCSumT1.',  ErrStat,ErrMsg,'Waves2_Init')

         ALLOCATE ( WaveAcc2xCSumT1    (0:InitInp%NStepWave2), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveAcc2xCSumT1.',  ErrStat,ErrMsg,'Waves2_Init')
         ALLOCATE ( WaveAcc2yCSumT1    (0:InitInp%NStepWave2), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveAcc2yCSumT1.',  ErrStat,ErrMsg,'Waves2_Init')
         ALLOCATE ( WaveAcc2zCSumT1    (0:InitInp%NStepWave2), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveAcc2zCSumT1.',  ErrStat,ErrMsg,'Waves2_Init')

         ALLOCATE ( WaveDynP2CSumT1    (0:InitInp%NStepWave2), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveDynP2CSumT1.',  ErrStat,ErrMsg,'Waves2_Init')

            ! Term 2 (n/=m term)
         ALLOCATE ( WaveVel2xCSumT2    (0:InitInp%NStepWave2), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveVel2xCSumT2.',  ErrStat,ErrMsg,'Waves2_Init')
         ALLOCATE ( WaveVel2yCSumT2    (0:InitInp%NStepWave2), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveVel2yCSumT2.',  ErrStat,ErrMsg,'Waves2_Init')
         ALLOCATE ( WaveVel2zCSumT2    (0:InitInp%NStepWave2), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveVel2zCSumT2.',  ErrStat,ErrMsg,'Waves2_Init')

         ALLOCATE ( WaveAcc2xCSumT2    (0:InitInp%NStepWave2), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveAcc2xCSumT2.',  ErrStat,ErrMsg,'Waves2_Init')
         ALLOCATE ( WaveAcc2yCSumT2    (0:InitInp%NStepWave2), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveAcc2yCSumT2.',  ErrStat,ErrMsg,'Waves2_Init')
         ALLOCATE ( WaveAcc2zCSumT2    (0:InitInp%NStepWave2), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveAcc2zCSumT2.',  ErrStat,ErrMsg,'Waves2_Init')

         ALLOCATE ( WaveDynP2CSumT2    (0:InitInp%NStepWave2), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveDynP2CSumT2.',  ErrStat,ErrMsg,'Waves2_Init')

            ! Now check if all the allocations worked properly
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF


            ! Time domain arrays: Term 1 (n=m term)

         ALLOCATE ( WaveVel2xSumT1   (0:InitInp%NStepWave), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveVel2xSumT1.',  ErrStat,ErrMsg,'Waves2_Init')
         ALLOCATE ( WaveVel2ySumT1   (0:InitInp%NStepWave), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveVel2ySumT1.',  ErrStat,ErrMsg,'Waves2_Init')
         ALLOCATE ( WaveVel2zSumT1   (0:InitInp%NStepWave), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveVel2zSumT1.',  ErrStat,ErrMsg,'Waves2_Init')

         ALLOCATE ( WaveAcc2xSumT1   (0:InitInp%NStepWave), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveAcc2xSumT1.',  ErrStat,ErrMsg,'Waves2_Init')
         ALLOCATE ( WaveAcc2ySumT1   (0:InitInp%NStepWave), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveAcc2ySumT1.',  ErrStat,ErrMsg,'Waves2_Init')
         ALLOCATE ( WaveAcc2zSumT1   (0:InitInp%NStepWave), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveAcc2zSumT1.',  ErrStat,ErrMsg,'Waves2_Init')

         ALLOCATE ( WaveDynP2SumT1   (0:InitInp%NStepWave), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveDynP2SumT1.',  ErrStat,ErrMsg,'Waves2_Init')

            ! Term 2 (n/=m term)
         ALLOCATE ( WaveVel2xSumT2   (0:InitInp%NStepWave), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveVel2xSumT2.',  ErrStat,ErrMsg,'Waves2_Init')
         ALLOCATE ( WaveVel2ySumT2   (0:InitInp%NStepWave), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveVel2ySumT2.',  ErrStat,ErrMsg,'Waves2_Init')
         ALLOCATE ( WaveVel2zSumT2   (0:InitInp%NStepWave), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveVel2zSumT2.',  ErrStat,ErrMsg,'Waves2_Init')

         ALLOCATE ( WaveAcc2xSumT2   (0:InitInp%NStepWave), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveAcc2xSumT2.',  ErrStat,ErrMsg,'Waves2_Init')
         ALLOCATE ( WaveAcc2ySumT2   (0:InitInp%NStepWave), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveAcc2ySumT2.',  ErrStat,ErrMsg,'Waves2_Init')
         ALLOCATE ( WaveAcc2zSumT2   (0:InitInp%NStepWave), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveAcc2zSumT2.',  ErrStat,ErrMsg,'Waves2_Init')

         ALLOCATE ( WaveDynP2SumT2   (0:InitInp%NStepWave), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveDynP2SumT2.',  ErrStat,ErrMsg,'Waves2_Init')

            ! Now check if all the allocations worked properly
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF




         !--------------------------------------------------------------------------------
         !> ## Calculate the surface elevation corrections ##
         !!
         !! For each (x,y) coordinate that a wave elevation is requested at (both from the
         !! (WaveElevxi,WaveElevyi) pairs, and the WaveElevXY pairs), a call is made to the
         !! subroutine waves2::waveelevtimeseriesatxy_sum to calculate the full time series for
         !! that point.  The results are added to the wave elevation results from the diff
         !! frequency calculations earlier in the code.
         !--------------------------------------------------------------------------------

             ! Step through the requested points
         DO I=1,InitInp%NWaveElev
            CALL WaveElevTimeSeriesAtXY_Sum(InitInp%WaveElevxi(I), InitInp%WaveElevyi(I), TmpTimeSeries, ErrStatTmp, ErrMsgTmp )
            CALL SetErrStat(ErrStatTmp,'Error occured while applying the FFT to InitOut%WaveElev.',ErrStat,ErrMsg,'Waves2_Init')
            IF ( ErrStat >= AbortErrLev ) THEN
               CALL CleanUp()
               RETURN
            END IF
               ! Add to the series since the difference is already included
            p%WaveElev2(:,I) = p%WaveElev2(:,I) + TmpTimeSeries(:)
         ENDDO    ! Wave elevation points requested


            ! Calculate the wave elevation at all points requested in the array WaveElevXY
         IF ( ALLOCATED(InitInp%WaveElevXY) ) THEN
            DO I = 1,SIZE(InitInp%WaveElevXY, DIM=2)
                  ! This subroutine call applies the FFT at the correct location.
               CALL WaveElevTimeSeriesAtXY_Sum( InitInp%WaveElevXY(1,I), InitInp%WaveElevXY(2,I), TmpTimeSeries, ErrStatTmp, ErrMsgTmp )
               CALL SetErrStat(ErrStatTmp,ErrMsgTmp,ErrStat,ErrMsg,'Waves2_Init')
               IF ( ErrStat >= AbortErrLev ) THEN
                  CALL CleanUp()
                  RETURN
               END IF
                  ! Add to the series since the difference is already included
               InitOut%WaveElevSeries2(:,I) = InitOut%WaveElevSeries2(:,I) + TmpTimeSeries(:)
            ENDDO
         ENDIF



         !--------------------------------------------------------------------------------
         !> ## Calculate the second order velocity, acceleration, and pressure corrections for all joints below surface. ##
         !--------------------------------------------------------------------------------
            ! NWaveKin0Prime loop start
         DO I=1,NWaveKin0Prime


               ! Reset the \f$ H_{\mu^+} \f$ terms to zero before calculating.
            WaveVel2xCSumT1 = CMPLX(0.0_SiKi, 0.0_SiKi, SiKi)
            WaveVel2yCSumT1 = CMPLX(0.0_SiKi, 0.0_SiKi, SiKi)
            WaveVel2zCSumT1 = CMPLX(0.0_SiKi, 0.0_SiKi, SiKi)
            WaveAcc2xCSumT1 = CMPLX(0.0_SiKi, 0.0_SiKi, SiKi)
            WaveAcc2yCSumT1 = CMPLX(0.0_SiKi, 0.0_SiKi, SiKi)
            WaveAcc2zCSumT1 = CMPLX(0.0_SiKi, 0.0_SiKi, SiKi)
            WaveDynP2CSumT1 = CMPLX(0.0_SiKi, 0.0_SiKi, SiKi)

            WaveVel2xCSumT2 = CMPLX(0.0_SiKi, 0.0_SiKi, SiKi)
            WaveVel2yCSumT2 = CMPLX(0.0_SiKi, 0.0_SiKi, SiKi)
            WaveVel2zCSumT2 = CMPLX(0.0_SiKi, 0.0_SiKi, SiKi)
            WaveAcc2xCSumT2 = CMPLX(0.0_SiKi, 0.0_SiKi, SiKi)
            WaveAcc2yCSumT2 = CMPLX(0.0_SiKi, 0.0_SiKi, SiKi)
            WaveAcc2zCSumT2 = CMPLX(0.0_SiKi, 0.0_SiKi, SiKi)
            WaveDynP2CSumT2 = CMPLX(0.0_SiKi, 0.0_SiKi, SiKi)


               !---------------
               !> ### First term ###
               !! This term is only the FFT over the diagonal elements where \f$ \omega_n = \omega_m \f$
               !! (note that the sum frequency is \f$ 2 \omega_n \f$).  The index for the sum frequency is
               !! therefore \f$ 2 n \f$.  Since we are placing the calculated value for the \f$ A_n A_n 
               !! H^+ \f$ term in the \f$ 2 \omega \f$ location, we will only run through the first
               !! half of the frequencies (the sum frequency will exceed the bounds of the frequencies
               !! used in the FFT otherwise).
               !! The IFFT will be calculated later with the IFFT of the second term.
               !---------------

               ! The limits look a little funny.  But remember we are placing the value in the 2*J location,
               ! so we cannot overun the end of the array.  The floor function is just in case NStepWave2 is
               ! an odd number
            DO n=1,FLOOR( REAL(InitInp%NStepWave2-1) / 2.0_SiKi )   ! Only

               Omega_n  =  n * InitInp%WaveDOmega

               ! The frequency we are dealing with
               !> * \f$ \omega^+ = \mu^+ \Delta \omega = 2 \omega_n \f$
               mu_plus     =  2 * n
               Omega_plus  =  2.0_SiKi * Omega_n

               IF ( Omega_plus >= InitInp%WvLowCOffS .AND. Omega_plus <= InitInp%WvHiCOffS ) THEN
                  k_n         =  WaveNumber( Omega_n, InitInp%Gravity, InitInp%WtrDpth )
                  k_nm        =  k_nm_plus( n, n, k_n, k_n )


                     ! Calculate the terms \f$ n,m \f$ necessary for calculations

                     !> Calculate the dot product of the wavenumbers with the (x,y) location
                     !! This is given by:
                     !!
                     !! *  \f$ \exp\left(-\imath 2 \vec{k_n} \cdot \vec{x} \right)
                     !!       = \exp \left( -\imath 2 \left[
                     !!                |\vec{k_n}| \cos \theta_n ~ x
                     !!             +  |\vec{k_n}| \sin \theta_n ~ y \right] \right) \f$

                  WaveElevxyPrime0  = exp( - ImagNmbr &
                           *  (  2.0_SiKi * k_n * COS( D2R_S*InitInp%WaveDirArr(n) ) * InitInp%WaveKinxi(WaveKinPrimeMap(I))  &
                              +  2.0_SiKi * k_n * SIN( D2R_S*InitInp%WaveDirArr(n) ) * InitInp%WaveKinyi(WaveKinPrimeMap(I))  ))


                     ! Get value for \f$ B+ \f$ for the n,m index pair
                  B_plus  =  TransFuncB_plus( n, n, k_n, k_n, WaveKinzi0Prime(I) )


                     !> Calculate \f$ U^+ \f$ terms for the velocity calculations (\f$B^+\f$ provided by waves2::transfuncb_plus)
                     ! NOTE: InitInp%WtrDpth + WaveKinzi0Prime(I) is the height above the ocean floor
                     !> * \f$ _x{U}_{nn}^+ = B_{nn}^+ 2 k_n \cos \theta_n \f$
                  Ux_nm_plus = B_plus * 2.0_SiKi * k_n * COS( D2R_S*InitInp%WaveDirArr(n) )

                     !> * \f$ _y{U}_{nn}^+ = B_{nn}^+ 2 k_n \sin \theta_n \f$
                  Uy_nm_plus = B_plus * 2.0_SiKi * k_n * SIN( D2R_S*InitInp%WaveDirArr(n) )

                     !> * \f$ _z{U}_{nn}^+ = \imath B_{nn}^+ k_{nn} \tanh \left( k_{nn} ( h + z ) \right) \f$
                  Uz_nm_plus = ImagNmbr * B_plus * k_nm * tanh( k_nm * ( InitInp%WtrDpth + WaveKinzi0Prime(I) ) )


                     !> Acceleration calculations
                  Accx_nm_plus = ImagNmbr * Ux_nm_plus * Omega_plus     !> * \f$ _x\dot{U}_{nn}^+ = \imath * _xU_{nn}^+ * \omega_{\mu^+} \f$
                  Accy_nm_plus = ImagNmbr * Uy_nm_plus * Omega_plus     !> * \f$ _y\dot{U}_{nn}^+ = \imath * _yU_{nn}^+ * \omega_{\mu^+} \f$
                  Accz_nm_plus = ImagNmbr * Uz_nm_plus * Omega_plus     !> * \f$ _z\dot{U}_{nn}^+ = \imath * _zU_{nn}^+ * \omega_{\mu^+} \f$


                     !> Dynamic pressure
                     !> * \f$ P_{nn}^+ = \rho_\mathrm{w} B_{nn}^+ \omega_{\mu^+} \f$
                  DynP_nm_plus  = REAL(InitInp%WtrDens, SiKi) * B_plus * Omega_plus



                  !> ### Calculate the array of \f$ K^+(\omega_n) \f$ for the first term of the velocity, acceleration, and pressure. ###

                     ! First get the wave amplitude -- must be reconstructed from the WaveElevC0 array.  First index is the real (1) or
                     ! imaginary (2) part.  Divide by NStepWave2 to remove the built in normalization in WaveElevC0.  Note that the phase
                     ! shift associated with the (x,y) location is accounted for by the WaveElevxyPrime0 variable.
                  WaveElevC_n =  WaveElevC0Norm(n)
 
                     !> Velocity terms:
                     !!    *  \f$ K^+(\omega_n) =  A_n A_n U_{nn}^+         \exp\left(-\imath 2 \vec{k_n} \cdot\vec{x}\right) \f$
                  WaveVel2xCSumT1(mu_plus)   =  WaveElevC_n * WaveElevC_n * Ux_nm_plus * WaveElevxyPrime0
                  WaveVel2yCSumT1(mu_plus)   =  WaveElevC_n * WaveElevC_n * Uy_nm_plus * WaveElevxyPrime0
                  WaveVel2zCSumT1(mu_plus)   =  WaveElevC_n * WaveElevC_n * Uz_nm_plus * WaveElevxyPrime0

                     !> Acceleration terms:
                     !!    *  \f$ K^+(\omega_n) =  A_n A_n \dot{U}_{nn}^+   \exp\left(-\imath 2 \vec{k_n} \cdot\vec{x}\right) \f$
                  WaveAcc2xCSumT1(mu_plus)   =  WaveElevC_n * WaveElevC_n * Accx_nm_plus * WaveElevxyPrime0
                  WaveAcc2yCSumT1(mu_plus)   =  WaveElevC_n * WaveElevC_n * Accy_nm_plus * WaveElevxyPrime0
                  WaveAcc2zCSumT1(mu_plus)   =  WaveElevC_n * WaveElevC_n * Accz_nm_plus * WaveElevxyPrime0

                     !> Pressure term:
                     !!    *  \f$ K^+(\omega_n) =  A_n A_n P_{nn}^+         \exp\left(-\imath 2 \vec{k_n} \cdot\vec{x}\right) \f$
                  WaveDynP2CSumT1(mu_plus)   =  WaveElevC_n * WaveElevC_n * DynP_nm_plus * WaveElevxyPrime0


               ENDIF ! Check to see if WvLowCOffS <= mu_plus <= WvHiCOffS

            ENDDO ! n loop (sum frequency)

               ! NOTE: The IFFT of the these terms is performed below.


               !---------------
               !> ### Second term ###
               !! In this term, we are are now stepping through the sum frequencies.  The inner
               !! sum essentially covers all the off diagonal terms (omega_m /= omega_n).  The limits
               !! on the outer integral that is the FFT run through the full frequency range that
               !! we are using.
               !---------------

               ! \f$ \mu^+ \f$ loop.  This loop is used to construct the full set of \f$ H_{\mu^+} \f$ terms used in the IFFT to find the timeseries.
               !> * \f$ \mu^+ = n + m \f$
            DO mu_plus=2,InitInp%NStepWave2-1

                  ! The frequency we are dealing with
                  !> * \f$ \omega^+ = \mu^+ \Delta \omega \f$
               Omega_plus =  mu_plus * InitInp%WaveDOmega

               IF ( Omega_plus >= InitInp%WvLowCOffS .AND. Omega_plus <= InitInp%WvHiCOffS ) THEN
                     ! The inner \f$ m \f$ loop for calculating the \f$ H_{\mu^+} \f$ terms at each frequency.
                  DO m=1,FLOOR( REAL(mu_plus - 1) / 2.0_SiKi )
                        ! Calculate the value of the n index from \f$ \mu^+ = n + m \f$.  Calculate corresponding wavenumbers and frequencies.
                     n           =  mu_plus - m
                     Omega_n     =  n * InitInp%WaveDOmega
                     Omega_m     =  m * InitInp%WaveDOmega
                     k_n         =  WaveNumber( Omega_n, InitInp%Gravity, InitInp%WtrDpth )
                     k_m         =  WaveNumber( Omega_m, InitInp%Gravity, InitInp%WtrDpth )
                     k_nm        =  k_nm_plus( n, m, k_n, k_m )


                        ! Calculate the terms \f$ n,m \f$ necessary for calculations

                        !> Calculate the dot product of the wavenumbers with the (x,y) location
                        !! This is given by:
                        !!
                        !! *  \f$ \exp\left(-\imath \left[\vec{k_n} + \vec{k_m} \right] \cdot \vec{x} \right)
                        !!       = \exp \left( -\imath \left[
                        !!                \left( |\vec{k_n}| \cos \theta_n + |\vec{k_m}| cos \theta_m \right) ~ x
                        !!             +  \left( |\vec{k_n}| \sin \theta_n + |\vec{k_m}| sin \theta_m \right) ~ y \right] \right) \f$

                     WaveElevxyPrime0  = exp( - ImagNmbr &
                              *  (  ( k_n * COS( D2R_S*InitInp%WaveDirArr(n) ) + k_m * COS( D2R_S*InitInp%WaveDirArr(m) ) ) * InitInp%WaveKinxi(WaveKinPrimeMap(I))  &
                                 +  ( k_n * SIN( D2R_S*InitInp%WaveDirArr(n) ) + k_m * SIN( D2R_S*InitInp%WaveDirArr(m) ) ) * InitInp%WaveKinyi(WaveKinPrimeMap(I))  ))


                        ! Get value for \f$ B+ \f$ for the n,m index pair
                     B_plus  =  TransFuncB_plus( n, m, k_n, k_m, WaveKinzi0Prime(I) )


                        !> Calculate \f$ U^+ \f$ terms for the velocity calculations (\f$B^+\f$ provided by waves2::transfuncb_plus)
                        ! NOTE: InitInp%WtrDpth + WaveKinzi0Prime(I) is the height above the ocean floor
                        !> * \f$ _x{U}_{nm}^+ = B_{nm}^+ \left(k_n \cos \theta_n + k_m \cos \theta_m \right) \f$
                     Ux_nm_plus = B_plus * ( k_n * COS( D2R_S*InitInp%WaveDirArr(n) ) + k_m * COS( D2R_S*InitInp%WaveDirArr(m) ) )

                        !> * \f$ _y{U}_{nm}^+ = B_{nm}^+ \left(k_n \sin \theta_n + k_m \sin \theta_m \right) \f$
                     Uy_nm_plus = B_plus * ( k_n * SIN( D2R_S*InitInp%WaveDirArr(n) ) + k_m * SIN( D2R_S*InitInp%WaveDirArr(m) ) )

                        !> * \f$ _z{U}_{nm}^+ = \imath B_{nm}^+ k_{nm} \tanh \left( k_{nm} ( h + z ) \right) \f$
                     Uz_nm_plus = ImagNmbr * B_plus * k_nm * tanh( k_nm * ( InitInp%WtrDpth + WaveKinzi0Prime(I) ) )


                        !> Acceleration calculations
                     Accx_nm_plus = ImagNmbr * Ux_nm_plus * Omega_plus     !> * \f$ _x\dot{U}_{nm}^+ = \imath * _xU_{nm}^+ * \omega_{\mu^+} \f$
                     Accy_nm_plus = ImagNmbr * Uy_nm_plus * Omega_plus     !> * \f$ _y\dot{U}_{nm}^+ = \imath * _yU_{nm}^+ * \omega_{\mu^+} \f$
                     Accz_nm_plus = ImagNmbr * Uz_nm_plus * Omega_plus     !> * \f$ _z\dot{U}_{nm}^+ = \imath * _zU_{nm}^+ * \omega_{\mu^+} \f$


                        !> Dynamic pressure
                        !> * \f$ P_{nm}^+ = \rho_\mathrm{w} B_{nm}^+ \omega_{\mu^+} \f$
                     DynP_nm_plus  = REAL(InitInp%WtrDens,SiKi) * B_plus * Omega_plus



                        !> ### Calculate the inner summation \f$ H^+(\omega_{\mu^+}) \f$ terms for the velocity, acceleration, and pressure. ###


                        ! First get the wave amplitude -- must be reconstructed from the WaveElevC0 array.  First index is the real (1) or
                        ! imaginary (2) part.  Divide by NStepWave2 to remove the built in normalization in WaveElevC0.  Note that the phase
                        ! shift associated with the (x,y) location is accounted for by the WaveElevxyPrime0 variable.
                     WaveElevC_n =  WaveElevC0Norm(n)
                     WaveElevC_m =  WaveElevC0Norm(m)
 

                        !> Velocity terms:
                        !!    *  \f$ H^+(\omega_{\mu^+}) =  \sum_{m=1}^{\lfloor \frac{\mu^+-1}{2}\rfloor}  A_n  A_m U_{nm}^+
                        !!                                  \exp\left(-\imath (\vec{k_n} + \vec{k_m})\cdot\vec{x}\right) \f$
                     WaveVel2xCSumT2(mu_plus)   =  WaveVel2xCSumT2(mu_plus) + WaveElevC_n * WaveElevC_m * Ux_nm_plus * WaveElevxyPrime0
                     WaveVel2yCSumT2(mu_plus)   =  WaveVel2yCSumT2(mu_plus) + WaveElevC_n * WaveElevC_m * Uy_nm_plus * WaveElevxyPrime0
                     WaveVel2zCSumT2(mu_plus)   =  WaveVel2zCSumT2(mu_plus) + WaveElevC_n * WaveElevC_m * Uz_nm_plus * WaveElevxyPrime0

                        !> Acceleration terms:
                        !!    *  \f$ H^+(\omega_{\mu^+}) =  \sum_{m=1}^{\lfloor \frac{\mu^+-1}{2}\rfloor}  A_n  A_m \dot{U}_{nm}^+
                        !!                                  \exp\left(-\imath (\vec{k_n} + \vec{k_m})\cdot\vec{x}\right) \f$
                     WaveAcc2xCSumT2(mu_plus)   =  WaveAcc2xCSumT2(mu_plus) + WaveElevC_n * WaveElevC_m * Accx_nm_plus * WaveElevxyPrime0
                     WaveAcc2yCSumT2(mu_plus)   =  WaveAcc2yCSumT2(mu_plus) + WaveElevC_n * WaveElevC_m * Accy_nm_plus * WaveElevxyPrime0
                     WaveAcc2zCSumT2(mu_plus)   =  WaveAcc2zCSumT2(mu_plus) + WaveElevC_n * WaveElevC_m * Accz_nm_plus * WaveElevxyPrime0

                        !> Pressure term:
                        !!    *  \f$ H^+(\omega_{\mu^+}) =  \sum_{m=1}^{\lfloor \frac{\mu^+-1}{2}\rfloor}  A_n  A_m P_{nm}^+
                        !!                                  \exp\left(-\imath (\vec{k_n} + \vec{k_m})\cdot\vec{x}\right) \f$
                     WaveDynP2CSumT2(mu_plus)   =  WaveDynP2CSumT2(mu_plus) + WaveElevC_n * WaveElevC_m * DynP_nm_plus * WaveElevxyPrime0

                  ENDDO ! m loop

               ENDIF ! Check to see if WvLowCOffS <= mu_plus <= WvHiCOffS

            ENDDO ! mu_plus loop (diff frequency)


                  !  Divide by two for the single sided FFT given in the documentation.
            WaveVel2xCSumT1 =  WaveVel2xCSumT1 / 2.0_SiKi
            WaveVel2yCSumT1 =  WaveVel2yCSumT1 / 2.0_SiKi
            WaveVel2zCSumT1 =  WaveVel2zCSumT1 / 2.0_SiKi
            WaveAcc2xCSumT1 =  WaveAcc2xCSumT1 / 2.0_SiKi
            WaveAcc2yCSumT1 =  WaveAcc2yCSumT1 / 2.0_SiKi
            WaveAcc2zCSumT1 =  WaveAcc2zCSumT1 / 2.0_SiKi
            WaveDynP2CSumT1 =  WaveDynP2CSumT1 / 2.0_SiKi
            WaveVel2xCSumT2 =  WaveVel2xCSumT2 / 2.0_SiKi
            WaveVel2yCSumT2 =  WaveVel2yCSumT2 / 2.0_SiKi
            WaveVel2zCSumT2 =  WaveVel2zCSumT2 / 2.0_SiKi
            WaveAcc2xCSumT2 =  WaveAcc2xCSumT2 / 2.0_SiKi
            WaveAcc2yCSumT2 =  WaveAcc2yCSumT2 / 2.0_SiKi
            WaveAcc2zCSumT2 =  WaveAcc2zCSumT2 / 2.0_SiKi
            WaveDynP2CSumT2 =  WaveDynP2CSumT2 / 2.0_SiKi




               !> ### Apply the inverse FFT to the first and second terms of each of the components to get the time domain result ###
               !> *   \f$ V^{(2)+}(t)  =  \operatorname{IFFT}\left[K^+\right]
               !!                      + 2\operatorname{IFFT}\left[H^+\right]     \f$
            CALL ApplyFFT_cx(  WaveVel2xSumT1(:),  WaveVel2xCSumT1(:), FFT_Data, ErrStatTmp )
               CALL SetErrStat(ErrStatTmp,'Error occured while applying the FFT on V_x.',ErrStat,ErrMsg,'Waves2_Init')
            CALL ApplyFFT_cx(  WaveVel2ySumT1(:),  WaveVel2yCSumT1(:), FFT_Data, ErrStatTmp )
               CALL SetErrStat(ErrStatTmp,'Error occured while applying the FFT on V_y.',ErrStat,ErrMsg,'Waves2_Init')
            CALL ApplyFFT_cx(  WaveVel2zSumT1(:),  WaveVel2zCSumT1(:), FFT_Data, ErrStatTmp )
               CALL SetErrStat(ErrStatTmp,'Error occured while applying the FFT on V_z.',ErrStat,ErrMsg,'Waves2_Init')

            CALL ApplyFFT_cx(  WaveAcc2xSumT1(:),  WaveAcc2xCSumT1(:), FFT_Data, ErrStatTmp )
               CALL SetErrStat(ErrStatTmp,'Error occured while applying the FFT on Acc_x.',ErrStat,ErrMsg,'Waves2_Init')
            CALL ApplyFFT_cx(  WaveAcc2ySumT1(:),  WaveAcc2yCSumT1(:), FFT_Data, ErrStatTmp )
               CALL SetErrStat(ErrStatTmp,'Error occured while applying the FFT on Acc_y.',ErrStat,ErrMsg,'Waves2_Init')
            CALL ApplyFFT_cx(  WaveAcc2zSumT1(:),  WaveAcc2zCSumT1(:), FFT_Data, ErrStatTmp )
               CALL SetErrStat(ErrStatTmp,'Error occured while applying the FFT on Acc_z.',ErrStat,ErrMsg,'Waves2_Init')

            CALL ApplyFFT_cx(  WaveDynP2SumT1(:),  WaveDynP2CSumT1(:), FFT_Data, ErrStatTmp )
               CALL SetErrStat(ErrStatTmp,'Error occured while applying the FFT on DynP2.',ErrStat,ErrMsg,'Waves2_Init')

            CALL ApplyFFT_cx(  WaveVel2xSumT2(:),  WaveVel2xCSumT2(:), FFT_Data, ErrStatTmp )
               CALL SetErrStat(ErrStatTmp,'Error occured while applying the FFT on V_x.',ErrStat,ErrMsg,'Waves2_Init')
            CALL ApplyFFT_cx(  WaveVel2ySumT2(:),  WaveVel2yCSumT2(:), FFT_Data, ErrStatTmp )
               CALL SetErrStat(ErrStatTmp,'Error occured while applying the FFT on V_y.',ErrStat,ErrMsg,'Waves2_Init')
            CALL ApplyFFT_cx(  WaveVel2zSumT2(:),  WaveVel2zCSumT2(:), FFT_Data, ErrStatTmp )
               CALL SetErrStat(ErrStatTmp,'Error occured while applying the FFT on V_z.',ErrStat,ErrMsg,'Waves2_Init')

            CALL ApplyFFT_cx(  WaveAcc2xSumT2(:),  WaveAcc2xCSumT2(:), FFT_Data, ErrStatTmp )
               CALL SetErrStat(ErrStatTmp,'Error occured while applying the FFT on Acc_x.',ErrStat,ErrMsg,'Waves2_Init')
            CALL ApplyFFT_cx(  WaveAcc2ySumT2(:),  WaveAcc2yCSumT2(:), FFT_Data, ErrStatTmp )
               CALL SetErrStat(ErrStatTmp,'Error occured while applying the FFT on Acc_y.',ErrStat,ErrMsg,'Waves2_Init')
            CALL ApplyFFT_cx(  WaveAcc2zSumT2(:),  WaveAcc2zCSumT2(:), FFT_Data, ErrStatTmp )
               CALL SetErrStat(ErrStatTmp,'Error occured while applying the FFT on Acc_z.',ErrStat,ErrMsg,'Waves2_Init')

            CALL ApplyFFT_cx(  WaveDynP2SumT2(:),  WaveDynP2CSumT2(:), FFT_Data, ErrStatTmp )
               CALL SetErrStat(ErrStatTmp,'Error occured while applying the FFT on DynP2.',ErrStat,ErrMsg,'Waves2_Init')

            IF ( ErrStat >= AbortErrLev ) THEN
               CALL CleanUp()
               RETURN
            END IF


               ! Add the results to the output
            InitOut%WaveVel2S(:,WaveKinPrimeMap(I),1) =  WaveVel2xSumT1(:) +  2.0_SiKi * WaveVel2xSumT2(:)     ! x-component of velocity
            InitOut%WaveVel2S(:,WaveKinPrimeMap(I),2) =  WaveVel2ySumT1(:) +  2.0_SiKi * WaveVel2ySumT2(:)     ! y-component of velocity
            InitOut%WaveVel2S(:,WaveKinPrimeMap(I),3) =  WaveVel2zSumT1(:) +  2.0_SiKi * WaveVel2zSumT2(:)     ! z-component of velocity

            InitOut%WaveAcc2S(:,WaveKinPrimeMap(I),1) =  WaveAcc2xSumT1(:) +  2.0_SiKi * WaveAcc2xSumT2(:)     ! x-component of acceleration
            InitOut%WaveAcc2S(:,WaveKinPrimeMap(I),2) =  WaveAcc2ySumT1(:) +  2.0_SiKi * WaveAcc2ySumT2(:)     ! y-component of acceleration
            InitOut%WaveAcc2S(:,WaveKinPrimeMap(I),3) =  WaveAcc2zSumT1(:) +  2.0_SiKi * WaveAcc2zSumT2(:)     ! z-component of acceleration

            InitOut%WaveDynP2S(:,WaveKinPrimeMap(I))  =  WaveDynP2SumT1(:) +  2.0_SiKi * WaveDynP2SumT2(:)     ! Dynamic pressure


               ! Copy the first point to the last to make it easier.
            InitOut%WaveVel2S(InitInp%NStepWave,WaveKinPrimeMap(I),:)     =  InitOut%WaveVel2S(0,WaveKinPrimeMap(I),:)
            InitOut%WaveAcc2S(InitInp%NStepWave,WaveKinPrimeMap(I),:)     =  InitOut%WaveAcc2S(0,WaveKinPrimeMap(I),:)
            InitOut%WaveDynP2S(InitInp%NStepWave,WaveKinPrimeMap(I))    =  InitOut%WaveDynP2S(0,WaveKinPrimeMap(I))


         ENDDO    ! I=1,NWaveKin0Prime loop end


            ! Deallocate working arrays.
         IF (ALLOCATED(WaveVel2xCSumT1))   DEALLOCATE(WaveVel2xCSumT1,    STAT=ErrStatTmp)
         IF (ALLOCATED(WaveVel2yCSumT1))   DEALLOCATE(WaveVel2yCSumT1,    STAT=ErrStatTmp)
         IF (ALLOCATED(WaveVel2zCSumT1))   DEALLOCATE(WaveVel2zCSumT1,    STAT=ErrStatTmp)
         IF (ALLOCATED(WaveAcc2xCSumT1))   DEALLOCATE(WaveAcc2xCSumT1,    STAT=ErrStatTmp)
         IF (ALLOCATED(WaveAcc2yCSumT1))   DEALLOCATE(WaveAcc2yCSumT1,    STAT=ErrStatTmp)
         IF (ALLOCATED(WaveAcc2zCSumT1))   DEALLOCATE(WaveAcc2zCSumT1,    STAT=ErrStatTmp)
         IF (ALLOCATED(WaveDynP2CSumT1))   DEALLOCATE(WaveDynP2CSumT1,    STAT=ErrStatTmp)

         IF (ALLOCATED(WaveVel2xSumT1))    DEALLOCATE(WaveVel2xSumT1,     STAT=ErrStatTmp)
         IF (ALLOCATED(WaveVel2ySumT1))    DEALLOCATE(WaveVel2ySumT1,     STAT=ErrStatTmp)
         IF (ALLOCATED(WaveVel2zSumT1))    DEALLOCATE(WaveVel2zSumT1,     STAT=ErrStatTmp)
         IF (ALLOCATED(WaveAcc2xSumT1))    DEALLOCATE(WaveAcc2xSumT1,     STAT=ErrStatTmp)
         IF (ALLOCATED(WaveAcc2ySumT1))    DEALLOCATE(WaveAcc2ySumT1,     STAT=ErrStatTmp)
         IF (ALLOCATED(WaveAcc2zSumT1))    DEALLOCATE(WaveAcc2zSumT1,     STAT=ErrStatTmp)
         IF (ALLOCATED(WaveDynP2SumT1))    DEALLOCATE(WaveDynP2SumT1,     STAT=ErrStatTmp)

         IF (ALLOCATED(WaveVel2xCSumT2))   DEALLOCATE(WaveVel2xCSumT2,    STAT=ErrStatTmp)
         IF (ALLOCATED(WaveVel2yCSumT2))   DEALLOCATE(WaveVel2yCSumT2,    STAT=ErrStatTmp)
         IF (ALLOCATED(WaveVel2zCSumT2))   DEALLOCATE(WaveVel2zCSumT2,    STAT=ErrStatTmp)
         IF (ALLOCATED(WaveAcc2xCSumT2))   DEALLOCATE(WaveAcc2xCSumT2,    STAT=ErrStatTmp)
         IF (ALLOCATED(WaveAcc2yCSumT2))   DEALLOCATE(WaveAcc2yCSumT2,    STAT=ErrStatTmp)
         IF (ALLOCATED(WaveAcc2zCSumT2))   DEALLOCATE(WaveAcc2zCSumT2,    STAT=ErrStatTmp)
         IF (ALLOCATED(WaveDynP2CSumT2))   DEALLOCATE(WaveDynP2CSumT2,    STAT=ErrStatTmp)

         IF (ALLOCATED(WaveVel2xSumT2))    DEALLOCATE(WaveVel2xSumT2,     STAT=ErrStatTmp)
         IF (ALLOCATED(WaveVel2ySumT2))    DEALLOCATE(WaveVel2ySumT2,     STAT=ErrStatTmp)
         IF (ALLOCATED(WaveVel2zSumT2))    DEALLOCATE(WaveVel2zSumT2,     STAT=ErrStatTmp)
         IF (ALLOCATED(WaveAcc2xSumT2))    DEALLOCATE(WaveAcc2xSumT2,     STAT=ErrStatTmp)
         IF (ALLOCATED(WaveAcc2ySumT2))    DEALLOCATE(WaveAcc2ySumT2,     STAT=ErrStatTmp)
         IF (ALLOCATED(WaveAcc2zSumT2))    DEALLOCATE(WaveAcc2zSumT2,     STAT=ErrStatTmp)
         IF (ALLOCATED(WaveDynP2SumT2))    DEALLOCATE(WaveDynP2SumT2,     STAT=ErrStatTmp)

         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF



      ENDIF    ! p%WvSumQTFF






         CALL  ExitFFT(FFT_Data, ErrStatTmp)
         CALL  SetErrStat(ErrStatTmp,'Error occured while cleaning up after the FFTs.', ErrStat,ErrMsg,'Waves2_Init')
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF


            ! If we calculated any second order wave elevation corrections, the array TmpTimeSeries was used.  Deallocate it.
         IF (ALLOCATED(TmpTimeSeries))    DEALLOCATE(TmpTimeSeries,     STAT=ErrStatTmp)
         IF (ALLOCATED(TmpTimeSeries2))   DEALLOCATE(TmpTimeSeries2,    STAT=ErrStatTmp)
         IF (ALLOCATED(TmpFreqSeries))    DEALLOCATE(TmpFreqSeries,     STAT=ErrStatTmp)
         IF (ALLOCATED(TmpFreqSeries2))   DEALLOCATE(TmpFreqSeries2,    STAT=ErrStatTmp)


         ! initialize dummy variables for the framework, so that compilers don't complain that the INTENT(OUT) variables have not been set:
         u%DummyInput               = 0.0_SiKi
         x%DummyContState           = 0.0_SiKi
         xd%DummyDiscState          = 0.0_SiKi
         z%DummyConstrState         = 0.0_SiKi
         OtherState%DummyOtherState = 0_IntKi

         RETURN


   CONTAINS

   
      !> This subroutine calculates the second order difference frequency correction to the wave elevation.  The transfer function
      !! for the surface elevation is \f$ L^-_{nm} \f$, and is calculated within this subroutine rather than in a separate subroutine.
      !! The calculations in this routine follow the same basic structure usind in ::Waves2_Init for the difference frequency
      !! calculations.
      !!
      !! \f$
      !!    \eta^{(2)-}(t) =  2  \Re \left[ \sum_{\mu^-=1}^{\frac{N}{2}-1} H^-(\omega_{\mu^-})
      !!                         \exp(\imath \omega_{\mu^-} t) \right]
      !!                   =  \operatorname{IFFT}\left[2 H^-\right]     \f$
      !!
      !! Notice that in the equations that follow, there is no term analagous to the mean
      !! drift term in the WAMIT2 module.  Rather, for \f$ \omega_{\mu^-} = 0 \f$, the
      !! result is zero.  So this term is not included.
      !!
      !! Also notice that the multiplier 2 is moved inside the IFFT.  This was done purely to make the programming simpler.
      SUBROUTINE WaveElevTimeSeriesAtXY_Diff(Xcoord,Ycoord, WaveElevSeriesAtXY, ErrStatLcl, ErrMsgLcl )
   
         REAL(SiKi),       INTENT(IN   )              :: Xcoord
         REAL(SiKi),       INTENT(IN   )              :: Ycoord
         REAL(SiKi),       INTENT(  OUT)              :: WaveElevSeriesAtXY(0:InitInp%NStepWave)
         INTEGER(IntKi),   INTENT(  OUT)              :: ErrStatLcl
         INTEGER(IntKi)                               :: ErrStatLcl2
         CHARACTER(*),     INTENT(  OUT)              :: ErrMsgLcl

            ! Local variables
         INTEGER(IntKi)                               :: n              !< Index to the first  frequency we are dealing with
         INTEGER(IntKi)                               :: m              !< Index to the second frequency we are dealing with
         REAL(SiKi)                                   :: k_n            !< WaveNumber for Omega_n -- note no direction associated with this
         REAL(SiKi)                                   :: k_m            !< WaveNumber for Omega_m -- note no direction associated with this
         REAL(SiKi)                                   :: L_minus        !< Resulting \f$ L^{-}_{nm} \f$ value.  Calculated in this routine.
         REAL(SiKi)                                   :: R_n            !< Effect scaling relationship of depth and wavenumber
         REAL(SiKi)                                   :: R_m            !< Effect scaling relationship of depth and wavenumber
         REAL(SiKi)                                   :: Omega_n        !< First  frequency of index n
         REAL(SiKi)                                   :: Omega_m        !< Second frequency of index m
         REAL(SiKi)                                   :: D_minus        !< Value of \f$ D^-_{nm} \f$ found by ::TransFuncD_minus

            ! Initializations
         ErrMsgLcl   = ''
         ErrStatLcl  = ErrID_None
  
            ! Note that TmpFreqSeries was allocated in the calling routine.  Probably bad programming
            ! practice, but I didn't want to have to allocate it at each point.
         TmpFreqSeries  =  CMPLX(0.0_SiKi, 0.0_SiKi, SiKi)
         WaveElevSeriesAtXY   =  0.0_SiKi

            ! \f$ \mu^- \f$ loop.  This loop is used to construct the full set of \f$ H_{\mu^-} \f$ terms used in the IFFT to find the timeseries.
            !> * \f$ \mu^- = n -m \f$
         DO mu_minus=1,InitInp%NStepWave2-1

               ! The frequency we are dealing with
               !> * \f$ \omega^- = \mu^- \Delta \omega \f$
            Omega_minus =  mu_minus * InitInp%WaveDOmega

            IF ( Omega_minus >= InitInp%WvLowCOffD .AND. Omega_minus <= InitInp%WvHiCOffD ) THEN

                  ! The inner \f$ m \f$ loop for calculating the \f$ H_{\mu^-} \f$ terms at each frequency.
               DO m=1,InitInp%NStepWave2-mu_minus
                     ! Calculate the value of the n index from \f$ \mu^- = n - m \f$.  Calculate corresponding wavenumbers and frequencies.
                  n           =  mu_minus + m
                  Omega_n     =  n * InitInp%WaveDOmega
                  Omega_m     =  m * InitInp%WaveDOmega
                  k_n         =  WaveNumber( Omega_n, InitInp%Gravity, InitInp%WtrDpth )
                  k_m         =  WaveNumber( Omega_m, InitInp%Gravity, InitInp%WtrDpth )
                  R_n         =  k_n * tanh( k_n * InitInp%WtrDpth )
                  R_m         =  k_m * tanh( k_m * InitInp%WtrDpth )
                  D_minus     =  TransFuncD_minus(n,m,k_n,k_m,R_n,R_m)

                     !> Calculate the value of 
                     !!    \f$ L^-_{nm} = \frac{1}{4} \left[ 
                     !!             \frac{D^-_{nm} - |\vec{k}_n| |\vec{k}_m| \cos(\theta_n - \theta_m) - R_n R_m}{\sqrt{R_n R_m}}
                     !!          +  (R_n+R_m) \right] \f$
                     !!
                     !!    The value of \f$ D^-_{nm} \f$ is found from by the ::TransFuncD_minus routine.

                  L_minus  =  (( D_minus - k_n * k_m * COS(D2R_S*InitInp%WaveDirArr(n) - D2R_S*InitInp%WaveDirArr(m)) - R_n * R_m )/SQRT( R_n * R_m ) + R_n + R_m) / 4.0_SiKi !4.0_SiKi


                     ! Calculate the terms \f$ n,m \f$ necessary for calculations

                     !> Calculate the dot product of the wavenumbers with the (x,y) location
                     !! This is given by:
                     !!
                     !! *  \f$ \exp\left(-\imath \left[\vec{k_n} - \vec{k_m} \right] \cdot \vec{x} \right)
                     !!       = \exp \left( -\imath \left[
                     !!                \left( |\vec{k_n}| \cos \theta_n - |\vec{k_m}| cos \theta_m \right) ~ x
                     !!             +  \left( |\vec{k_n}| \sin \theta_n - |\vec{k_m}| sin \theta_m \right) ~ y \right] \right) \f$

                  WaveElevxyPrime0  = exp( - ImagNmbr &
                           *  ( ( k_n * COS( D2R_S*InitInp%WaveDirArr(n) ) - k_m * COS( D2R_S*InitInp%WaveDirArr(m) ) ) * XCoord  &
                              + ( k_n * SIN( D2R_S*InitInp%WaveDirArr(n) ) - k_m * SIN( D2R_S*InitInp%WaveDirArr(m) ) ) * YCoord  ))


                     !> ### Calculate the inner summation \f$ H^-(\omega_{\mu^-}) \f$ terms for the velocity, acceleration, and pressure. ###

                     ! First get the wave amplitude -- must be reconstructed from the WaveElevC0 array.  First index is the real (1) or
                     ! imaginary (2) part.  Divide by NStepWave2 to remove the built in normalization in WaveElevC0.  Note that the phase
                     ! shift associated with the (x,y) location is accounted for by the WaveElevxyPrime0 variable.
                  WaveElevC_n =  WaveElevC0Norm(n)
                  WaveElevC_m =  WaveElevC0Norm(m)
 

                     !> Wave elevation term:
                     !!    *  \f$ 2 H^-(\omega_{\mu^-}) =  {\sum_{m=1}^{\frac{N}{2}-\mu^{-}}}  2 A_n  A^*_m L_{nm}^-
                     !!                                  \exp\left(-\imath (\vec{k_n} - \vec{k_m})\cdot\vec{x}\right) \f$
                  TmpFreqSeries(mu_minus)   =  TmpFreqSeries(mu_minus) + 2.0_SiKi * WaveElevC_n * CONJG( WaveElevC_m ) * L_minus * WaveElevxyPrime0


               ENDDO ! m loop

            ENDIF ! Check to see if WvLowCOffD <= mu_minus <= WvHiCOffD

         ENDDO ! mu_minus loop (diff frequency)


                  !  Divide by two for the single sided FFT given in the documentation.
            TmpFreqSeries = TmpFreqSeries / 2.0_SiKi


            !> ### Apply the inverse FFT to each of the components to get the time domain result ###
            !> *   \f$ \eta(t) = \operatorname{IFFT}\left[2 H^-\right] \f$
         CALL ApplyFFT_cx( WaveElevSeriesAtXY(:), TmpFreqSeries(:), FFT_Data, ErrStatLcl2 )
         CALL SetErrStat(ErrStatLcl2,'Error occured while applying the FFT on WaveElevSeriesAtXY.',ErrStatLcl,ErrMsgLcl,'WaveElevSeriesAtXY_Diff')
 
            ! Append first datapoint as the last as aid for repeated wave data
         WaveElevSeriesAtXY(InitInp%NStepWave) = WaveElevSeriesAtXY(0)
   

      END SUBROUTINE WaveElevTimeSeriesAtXY_Diff



      !> This subroutine calculates the second order sum frequency correction to the wave elevation.  The transfer function
      !! for the surface elevation is \f$ L^+_{nm} \f$, and is calculated within this subroutine rather than in a separate subroutine.
      !! The calculations in this routine follow the same basic structure usind in ::Waves2_Init for the sum frequency
      !! calculations.
      !!
      !! \f$
      !!    \eta^{(2)+}(t) =  \Re \left[ \sum_{n=1}^{\lfloor\frac{N}{4}\rfloor} K^+ \exp(\imath 2\omega_n t) \right]
      !!                   + 2\Re \left[ \sum_{\mu^+=2}^{\frac{N}{2}} H^+(\omega_{\mu^+})
      !!                         \exp(\imath \omega_{\mu^+} t) \right]
      !!                   =  \operatorname{IFFT}\left[K^+\right] + 2\operatorname{IFFT}\left[H^+\right]     \f$
      !!
      SUBROUTINE WaveElevTimeSeriesAtXY_Sum(Xcoord,Ycoord, WaveElevSeriesAtXY, ErrStatLcl, ErrMsgLcl )
   
         REAL(SiKi),       INTENT(IN   )              :: Xcoord
         REAL(SiKi),       INTENT(IN   )              :: Ycoord
         REAL(SiKi),       INTENT(  OUT)              :: WaveElevSeriesAtXY(0:InitInp%NStepWave)
         INTEGER(IntKi),   INTENT(  OUT)              :: ErrStatLcl
         INTEGER(IntKi)                               :: ErrStatLcl2
         CHARACTER(*),     INTENT(  OUT)              :: ErrMsgLcl

            ! Local variables
         INTEGER(IntKi)                               :: n              !< Index to the first  frequency we are dealing with
         INTEGER(IntKi)                               :: m              !< Index to the second frequency we are dealing with
         INTEGER(IntKi)                               :: Ctr            !< Generic counter
         REAL(SiKi)                                   :: k_n            !< WaveNumber for Omega_n -- note no direction associated with this
         REAL(SiKi)                                   :: k_m            !< WaveNumber for Omega_m -- note no direction associated with this
         REAL(SiKi)                                   :: L_plus         !< Resulting \f$ L^{+}_{nm} \f$ value.  Calculated in this routine.
         REAL(SiKi)                                   :: R_n            !< Effect scaling relationship of depth and wavenumber
         REAL(SiKi)                                   :: R_m            !< Effect scaling relationship of depth and wavenumber
         REAL(SiKi)                                   :: Omega_n        !< First  frequency of index n
         REAL(SiKi)                                   :: Omega_m        !< Second frequency of index m
         REAL(SiKi)                                   :: D_plus         !< Value of \f$ D^+_{nm} \f$ found by ::TransFuncD_plus

            ! Initializations
         ErrMsgLcl   = ''
         ErrStatLcl  = ErrID_None
  
            ! Note that TmpFreqSeries was allocated in the calling routine.  Probably bad programming
            ! practice, but I didn't want to have to allocate it at each point.
         TmpFreqSeries  =  CMPLX(0.0_SiKi, 0.0_SiKi, SiKi)     ! used for first term
         TmpFreqSeries2 =  CMPLX(0.0_SiKi, 0.0_SiKi, SiKi)     ! used for second term
         WaveElevSeriesAtXY   =  0.0_SiKi


            !> ## First term ##
            ! First term results are stored in TmpFreqSeries.

         DO n=1,FLOOR( REAL(InitInp%NStepWave2-1) / 2.0_SiKi )   ! Only

            Omega_n  =  n * InitInp%WaveDOmega

            ! The frequency we are dealing with
            !> * \f$ \omega^+ = \mu^+ \Delta \omega = 2 \omega_n \f$
            mu_plus     =  2 * n
            Omega_plus  =  2.0_SiKi * Omega_n

            IF ( Omega_plus >= InitInp%WvLowCOffS .AND. Omega_plus <= InitInp%WvHiCOffS ) THEN
               k_n         =  WaveNumber( Omega_n, InitInp%Gravity, InitInp%WtrDpth )
               R_n         =  k_n * tanh( k_n * InitInp%WtrDpth )
               D_plus      =  TransFuncD_plus(n,n,k_n,k_n,R_n,R_n)

                  !> Calculate the value of 
                  !!    \f$ L^+_{nn} = \frac{1}{4} \left[ 
                  !!             \frac{D^+_{nn} - |\vec{k}_n| |\vec{k}_n| \cos(\theta_n - \theta_n) + R_n R_n}{\sqrt{R_n R_n}}
                  !!          +  (R_n+R_n) \right]
                  !!       =  \frac{1}{4} \left[ \frac{ D^+_{nn} - |\vec{k}_n|^2 + R_n^2 }{ R_n } + 2 R_n \right] \f$
                  !!
                  !!    The value of \f$ D^+_{nn} \f$ is found from by the ::TransFuncD_plus routine.
               L_plus  =  (( D_plus - k_n * k_n + R_n * R_n )/R_n + 2.0_SiKi * R_n ) / 4.0_SiKi

                  !> Calculate the dot product of the wavenumbers with the (x,y) location
                  !! This is given by:
                  !!
                  !! *  \f$ \exp\left(-\imath 2 \vec{k_n} \cdot \vec{x} \right)
                  !!       = \exp \left( -\imath 2 \left[
                  !!                |\vec{k_n}| \cos \theta_n ~ x
                  !!             +  |\vec{k_n}| \sin \theta_n ~ y \right] \right) \f$

               WaveElevxyPrime0  = exp( - ImagNmbr &
                        *  (  2.0_SiKi * k_n * COS( D2R_S*InitInp%WaveDirArr(n) ) * XCoord  &
                           +  2.0_SiKi * k_n * SIN( D2R_S*InitInp%WaveDirArr(n) ) * YCoord  ))

                  ! First get the wave amplitude -- must be reconstructed from the WaveElevC0 array.  First index is the real (1) or
                  ! imaginary (2) part.  Divide by NStepWave2 to remove the built in normalization in WaveElevC0.  Note that the phase
                  ! shift associated with the (x,y) location is accounted for by the WaveElevxyPrime0 variable.
               WaveElevC_n =  WaveElevC0Norm(n)
 
                 !> ### Calculate the array of \f$ K^+(\omega_n) \f$ for the first term of the velocity, acceleration, and pressure. ###
                 !! *  \f$ K^+(\omega_n) =  A_n A_n L_{nn}^+         \exp\left(-\imath 2 \vec{k_n} \cdot\vec{x}\right) \f$
               TmpFreqSeries(mu_plus) = WaveElevC_n * WaveElevC_n * L_plus * WaveElevxyPrime0

            ENDIF ! Check to see if WvLowCOffS <= mu_plus <= WvHiCOffS

         ENDDO ! n loop (diff frequency)

            ! NOTE: The IFFT of the these terms is performed below.


            !---------------
            !> ## Second term ##
            !! In this term, we are are now stepping through the sum frequencies.  The inner
            !! sum essentially covers all the off diagonal terms (omega_m /= omega_n).  The limits
            !! on the outer integral that is the FFT run through the full frequency range that
            !! we are using.
            !---------------
            ! The frequency information will be stored in TmpFreqSeries2

            ! \f$ \mu^+ \f$ loop.  This loop is used to construct the full set of \f$ H_{\mu^+} \f$ terms used in the IFFT to find the timeseries.
            !> * \f$ \mu^+ = n + m \f$
         DO mu_plus=2,InitInp%NStepWave2-1

               ! The frequency we are dealing with
               !> * \f$ \omega^+ = \mu^+ \Delta \omega \f$
            Omega_plus =  mu_plus * InitInp%WaveDOmega

            IF ( Omega_plus >= InitInp%WvLowCOffS .AND. Omega_plus <= InitInp%WvHiCOffS ) THEN

                  ! The inner \f$ m \f$ loop for calculating the \f$ H_{\mu^+} \f$ terms at each frequency.
               DO m=1,FLOOR( REAL(mu_plus - 1) / 2.0_SiKi )
                     ! Calculate the value of the n index from \f$ \mu^+ = n + m \f$.  Calculate corresponding wavenumbers and frequencies.
                  n           =  mu_plus - m
                  Omega_n     =  n * InitInp%WaveDOmega
                  Omega_m     =  m * InitInp%WaveDOmega
                  k_n         =  WaveNumber( Omega_n, InitInp%Gravity, InitInp%WtrDpth )
                  k_m         =  WaveNumber( Omega_m, InitInp%Gravity, InitInp%WtrDpth )
                  R_n         =  k_n * tanh( k_n * InitInp%WtrDpth )
                  R_m         =  k_m * tanh( k_m * InitInp%WtrDpth )
                  D_plus      =  TransFuncD_plus(n,m,k_n,k_m,R_n,R_m)

                     !> Calculate the value of 
                     !!    \f$ L^+_{nm} = \frac{1}{4} \left[ 
                     !!             \frac{D^+_{nm} - |\vec{k}_n| |\vec{k}_m| \cos(\theta_n - \theta_m) + R_n R_m}{\sqrt{R_n R_m}}
                     !!          +  (R_n+R_m) \right] \f$
                     !!
                     !!    The value of \f$ D^-_{nm} \f$ is found from by the ::TransFuncD_plus routine.
                  L_plus  =  (( D_plus - k_n * k_m * COS(D2R_S*InitInp%WaveDirArr(n) - D2R_S*InitInp%WaveDirArr(m)) + R_n * R_m )/SQRT( R_n * R_m ) + R_n + R_m) / 4.0_SiKi

                     !> Calculate the dot product of the wavenumbers with the (x,y) location
                     !! This is given by:
                     !!
                     !! *  \f$ \exp\left(-\imath \left[\vec{k_n} + \vec{k_m} \right] \cdot \vec{x} \right)
                     !!       = \exp \left( -\imath \left[
                     !!                \left( |\vec{k_n}| \cos \theta_n + |\vec{k_m}| cos \theta_m \right) ~ x
                     !!             +  \left( |\vec{k_n}| \sin \theta_n + |\vec{k_m}| sin \theta_m \right) ~ y \right] \right) \f$

                  WaveElevxyPrime0  = exp( - ImagNmbr &
                           *  (  ( k_n * COS( D2R_S*InitInp%WaveDirArr(n) ) + k_m * COS( D2R_S*InitInp%WaveDirArr(m) ) ) * XCoord  &
                              +  ( k_n * SIN( D2R_S*InitInp%WaveDirArr(n) ) + k_m * SIN( D2R_S*InitInp%WaveDirArr(m) ) ) * YCoord  ))



                     ! First get the wave amplitude -- must be reconstructed from the WaveElevC0 array.  First index is the real (1) or
                     ! imaginary (2) part.  Divide by NStepWave2 to remove the built in normalization in WaveElevC0.  Note that the phase
                     ! shift associated with the (x,y) location is accounted for by the WaveElevxyPrime0 variable.
                  WaveElevC_n =  WaveElevC0Norm(n)
                  WaveElevC_m =  WaveElevC0Norm(m)
 
                     !> ### Calculate the inner summation \f$ H^+(\omega_{\mu^+}) \f$ term. ###
                     !! * \f$ H^+(\omega_{\mu^+}) =  \sum_{m=1}^{\lfloor \frac{\mu^+-1}{2}\rfloor}  A_n  A_m L^+_{nm}
                     !!                                  \exp\left(-\imath (\vec{k_n} + \vec{k_m})\cdot\vec{x}\right) \f$
                  TmpFreqSeries2(mu_plus) = TmpFreqSeries2(mu_plus) + WaveElevC_n * WaveElevC_m * L_plus * WaveElevxyPrime0

               ENDDO ! m loop

            ENDIF ! Check to see if WvLowCOffS <= mu_plus <= WvHiCOffS

         ENDDO ! mu_plus loop (diff frequency)


                  !  Divide by two for the single sided FFT given in the documentation.
            TmpFreqSeries  = TmpFreqSeries / 2.0_SiKi
            TmpFreqSeries2 = TmpFreqSeries2 / 2.0_SiKi

            !> ## Apply the inverse FFT to the first and second terms to get the time domain result ##
            !> *   \f$ \eta^{(2)+}(t)  =  \operatorname{IFFT}\left[K^+\right]
            !!                         + 2\operatorname{IFFT}\left[H^+\right]     \f$
         CALL ApplyFFT_cx( WaveElevSeriesAtXY(:),  TmpFreqSeries(:), FFT_Data, ErrStatLcl2 )
         CALL SetErrStat(ErrStatLcl2,'Error occured while applying the FFT on WaveElevSeriesAtXY.',ErrStatLcl,ErrMsgLcl,'WaveElevSeriesAtXY_Sum')
         CALL ApplyFFT_cx( TmpTimeSeries2(:),      TmpFreqSeries2(:), FFT_Data, ErrStatLcl2 )
         CALL SetErrStat(ErrStatLcl2,'Error occured while applying the FFT on WaveElevSeriesAtXY.',ErrStatLcl,ErrMsgLcl,'WaveElevSeriesAtXY_Sum')

            ! Add the two terms together
         DO Ctr=0,InitInp%NStepWave
            WaveElevSeriesAtXY(Ctr) =  WaveElevSeriesAtXY(Ctr)  +  2.0_SiKi * TmpTimeSeries2(Ctr)
         ENDDO
 
            ! Append first datapoint as the last as aid for repeated wave data
         WaveElevSeriesAtXY(InitInp%NStepWave) = WaveElevSeriesAtXY(0)
   

      END SUBROUTINE WaveElevTimeSeriesAtXY_Sum






      !> This function calculates the term \f$ B^-_{nm} \f$ used in calculating the veloicty, acceleration, and dynamic pressure terms.
      !! The equation is given by:
      !!
      !! \f$ B_{nm}^-(z, \omega_n, \omega_m, \theta_n, \theta_m) =\frac{g^2}{\omega_n \omega_m}
      !!             \cdot    \frac{1}{4} \frac{\cosh\left[k_{nm}^-(h+z)\right]}{\cosh\left[k_{nm}^-(h)\right]}
      !!                      \frac{D_{nm}^-}{\omega_n - \omega_m}   \f$
      !!
      FUNCTION TransFuncB_minus(n,m,k_n,k_m,z)

            ! Passed variables
         INTEGER(IntKi),               INTENT(IN   )  :: n              !< Index to the first  frequency we are dealing with
         INTEGER(IntKi),               INTENT(IN   )  :: m              !< Index to the second frequency we are dealing with
         REAL(SiKi),                   INTENT(IN   )  :: k_n            !< WaveNumber for Omega_n -- note no direction associated with this
         REAL(SiKi),                   INTENT(IN   )  :: k_m            !< WaveNumber for Omega_m -- note no direction associated with this
         REAL(SiKi),                   INTENT(IN   )  :: z              !< The depth of the point of interest from the surface of the water.
         REAL(SiKi)                                   :: TransFuncB_minus  !< Resulting \f$ B^{-}_{nm} \f$ value.

            ! Local variables
         REAL(SiKi)                                   :: R_n            !< Effect scaling relationship of depth and wavenumber
         REAL(SiKi)                                   :: R_m            !< Effect scaling relationship of depth and wavenumber
         REAL(SiKi)                                   :: k_nm           !< Value of \f$ k_{nm}^{-} \f$ found by ::k_nm_minus
         REAL(SiKi)                                   :: Omega_n        !< First  frequency of index n
         REAL(SiKi)                                   :: Omega_m        !< Second frequency of index m
         REAL(SiKi)                                   :: D_minus        !< Value of \f$ D^-_{nm} \f$ found by ::TransFuncD_minus

            ! Check that we are not trying to compute a transfer function with any zero frequencies in it.  Those are by definition zero.
         IF ( n==0 .or. m==0 ) THEN

            TransFuncB_minus  = 0.0_SiKi

         ELSEIF ( n==m ) THEN

               ! If the frequencies are the same, we get a zero in the denominator.  These should be defined as zero.
            TransFuncB_minus  = 0.0_SiKi

         ELSE

               ! Frequencies
            Omega_n     =  n * InitInp%WaveDOmega
            Omega_m     =  m * InitInp%WaveDOmega

               ! Wavenumbers
            k_nm        =  k_nm_minus( n,m,k_n,k_m )

               ! Effect of depth scaling
            R_n         =  k_n * tanh( k_n * InitInp%WtrDpth )
            R_m         =  k_m * tanh( k_m * InitInp%WtrDpth )

               ! Transfer function D_minus
            D_minus     =  TransFuncD_minus(n,m,k_n,k_m,R_n,R_m)


               ! Calculation of B_minus
            TransFuncB_minus  =  REAL(InitInp%Gravity*InitInp%Gravity,SiKi) / ( 4.0_SiKi * Omega_n * Omega_m ) &          
                                 * COSHNumOvrCOSHDen(k_nm, REAL(InitInp%WtrDpth,SiKi), z)  * D_minus / ( Omega_n - Omega_m )


         ENDIF



      END FUNCTION TransFuncB_minus




      !> This function calculates the term \f$ B^+_{nm} \f$ used in calculating the velocity, acceleration, and dynamic pressure terms.
      !! The equation is given by:
      !!
      !! \f$ B_{nm}^+(z, \omega_n, \omega_m, \theta_n, \theta_m) =\frac{g^2}{\omega_n \omega_m}
      !!             \cdot    \frac{1}{4} \frac{\cosh\left[k_{nm}^-(h+z)\right]}{\cosh\left[k_{nm}^-(h)\right]}
      !!                      \frac{D_{nm}^+}{\omega_n + \omega_m}   \f$
      !!
      FUNCTION TransFuncB_plus(n,m,k_n,k_m,z)

            ! Passed variables
         INTEGER(IntKi),               INTENT(IN   )  :: n              !< Index to the first  frequency we are dealing with
         INTEGER(IntKi),               INTENT(IN   )  :: m              !< Index to the second frequency we are dealing with
         REAL(SiKi),                   INTENT(IN   )  :: k_n            !< WaveNumber for Omega_n -- note no direction associated with this
         REAL(SiKi),                   INTENT(IN   )  :: k_m            !< WaveNumber for Omega_m -- note no direction associated with this
         REAL(SiKi),                   INTENT(IN   )  :: z              !< The depth of the point of interest from the surface of the water.
         REAL(SiKi)                                   :: TransFuncB_plus  !< Resulting \f$ B^{-}_{nm} \f$ value.

            ! Local variables
         REAL(SiKi)                                   :: R_n            !< Effect scaling relationship of depth and wavenumber
         REAL(SiKi)                                   :: R_m            !< Effect scaling relationship of depth and wavenumber
         REAL(SiKi)                                   :: k_nm           !< Value of \f$ k_{nm}^{+} \f$ found by ::k_nm_plus
         REAL(SiKi)                                   :: Omega_n        !< First  frequency of index n
         REAL(SiKi)                                   :: Omega_m        !< Second frequency of index m
         REAL(SiKi)                                   :: D_plus        !< Value of \f$ D^+_{nm} \f$ found by ::TransFuncD_plus

            ! Check that we are not trying to compute a transfer function with any zero frequencies in it.  Those are by definition zero.
         IF ( n==0 .or. m==0 ) THEN

            TransFuncB_plus  = 0.0_SiKi


         ELSE

               ! Frequencies
            Omega_n     =  n * InitInp%WaveDOmega
            Omega_m     =  m * InitInp%WaveDOmega

               ! Wavenumbers
            k_nm        =  k_nm_plus( n,m,k_n,k_m )

               ! Effect of depth scaling
            R_n         =  k_n * tanh( k_n * InitInp%WtrDpth )
            R_m         =  k_m * tanh( k_m * InitInp%WtrDpth )

               ! Transfer function D_plus
            D_plus     =  TransFuncD_plus(n,m,k_n,k_m,R_n,R_m)

               ! Calculation of B_plus
            TransFuncB_plus  =  REAL(InitInp%Gravity*InitInp%Gravity,SiKi) / ( 4.0_SiKi * Omega_n * Omega_m ) &
                                 * COSHNumOvrCOSHDen(k_nm, REAL(InitInp%WtrDpth,SiKi), z)  * D_plus / ( Omega_n + Omega_m )


         ENDIF



      END FUNCTION TransFuncB_plus






      !> This function was taken directly from the Waves.f90 file and should be identical to it.
      FUNCTION COSHNumOvrCOSHDen ( k, h, z )

            ! This FUNCTION computes the shallow water hyperbolic numerator
            ! over denominator term in the wave kinematics expressions:
            !
            !                    COSH( k*( z + h ) )/COSH( k*h )
            !
            ! given the wave number, k, water depth, h, and elevation z, as
            ! inputs.

         IMPLICIT                        NONE

            ! Passed Variables:
         REAL(SiKi)                   :: COSHNumOvrCOSHDen                               ! This function = COSH( k*( z + h ) )/COSH( k*h ) (-)
         REAL(SiKi), INTENT(IN )      :: h                                               ! Water depth ( h      >  0 ) (meters)
         REAL(SiKi), INTENT(IN )      :: k                                               ! Wave number ( k      >= 0 ) (1/m)
         REAL(SiKi), INTENT(IN )      :: z                                               ! Elevation   (-h <= z <= 0 ) (meters)


            ! Compute the hyperbolic numerator over denominator:

         IF ( k*h  > 89.4_SiKi )  THEN    ! When .TRUE., the shallow water formulation will trigger a floating point overflow error;
                                          ! however, COSH( k*( z + h ) )/COSH( k*h ) = EXP( k*z ) + EXP( -k*( z + 2*h ) ) for large k*h.
                                          ! This equals the deep water formulation, EXP( k*z ), except near z = -h, because
                                          ! h > 14.23*wavelength (since k = 2*Pi/wavelength) in this case.

            COSHNumOvrCOSHDen = EXP( k*z ) + EXP( -k*( z + 2.0_SiKi*h ) )

         ELSE                       ! 0 < k*h <= 89.4; use the shallow water formulation.

            COSHNumOvrCOSHDen = COSH( k*( z + h ) )/COSH( k*h )

         END IF

         RETURN
      END FUNCTION COSHNumOvrCOSHDen


      !> This function calculates the term \f$ D^-_{nm} \f$ used in finding the transfer functions.
      !! The equation is given by:
      !!
      !! \f$ {D}_{nm}^{-} =
      !!          \frac {   \left(\sqrt{R_n} - \sqrt{R_m}\right) \left[ \sqrt{R_m} \left( k_n^2 - R_n^2 \right) - \sqrt{R_n} \left( k_m^{2} - R_m^2 \right) \right]
      !!             ~+~  2 \left(\sqrt{R_n} - \sqrt{R_m}\right)^2
      !!                    \left[ \left|\vec{k_n}\right| \left|\vec{k_m}\right| \cos \left( \theta_n-\theta_m \right)  + R_n R_m  \right] }
      !!                {   \left(\sqrt{R_n} - \sqrt{R_m}\right)^2   -  k_{nm}^{-} \tanh \left( k_{nm}^{-} h \right)   }  \f$
      !!
      !! where \f$ k_{nm}^{-} \f$ is handled by the function ::k_nm_minus and \f$R_n\f$ is given by
      !!    \f$ R_n = \left| \overrightarrow{k}_n \right| \tanh\left(\left| \overrightarrow{k}_n \right| h \right) \f$
      !! where \f$ h \f$ is the depth from MSL (or the still water line).
      !!
      !! To calculate this, we simplify some of the common pieces:
      !!
      !! \f$ {D}_{nm}^{-}  = \frac  {     R_{nm}         \left[ \sqrt{R_m}  \left( k_n^2   - R_n^2 \right)  -  \sqrt{R_n} \left( k_m^2 -R_m^2 \right) \right]
      !!                           ~+~ 2  R_{nm}^2       \left[ \left|\vec{k_n}\right| \left|\vec{k_m}\right| \cos \left( \theta_n-\theta_m \right)  + R_n R_m  \right] }
      !!                            {     R_{nm}^2    -  k_{nm}^{-} \tanh \left( k_{nm}^- h \right)   }  \f$
      !!
      !! where \f$   R_{nm} \equiv \sqrt{R_n} - \sqrt{R_m} \f$.
      !!
      !! The denominator goes to zero when \f$ n = m \f$, as does the numerator.  So, using L'Hopital's rule to find the limit, it may be
      !! possible to prove that \f$  \stackrel{\lim}{n \to m} D_{nm}^{-} = 0 \f$.  It will take more than one derivative
      !! to check this, so due to time, this has not been verified.  For now we we simply assume this to be true and proceed to set
      !! _TransFuncD_minus_ to zero for all \f$ n = m \f$ cases. _This should not be done!_ Plotting this cross sections of this function
      !! shows that when \f$ n = m \f$, the result should converge to something non-zero, but close to zero.
      !!
      !! @note update this function when the limit has been derived.
      !!
      FUNCTION TransFuncD_minus(n,m,k_n,k_m,R_n,R_m)

            ! Passed variables
         INTEGER(IntKi),               INTENT(IN   )  :: n              !< Index to the first  frequency we are dealing with
         INTEGER(IntKi),               INTENT(IN   )  :: m              !< Index to the second frequency we are dealing with
         REAL(SiKi),                   INTENT(IN   )  :: k_n            !< WaveNumber for Omega_n -- note no direction associated with this
         REAL(SiKi),                   INTENT(IN   )  :: k_m            !< WaveNumber for Omega_m -- note no direction associated with this
         REAL(SiKi),                   INTENT(IN   )  :: R_n            !< Effect scaling relationship of depth and wavenumber
         REAL(SiKi),                   INTENT(IN   )  :: R_m            !< Effect scaling relationship of depth and wavenumber
         REAL(SiKi)                                   :: TransFuncD_minus  !< Resulting \f$ D^{-}_{nm} \f$ value.

            ! Local variables
         REAL(SiKi)                                   :: k_nm           !< Value of \f$ k_{nm}^{-} \f$
         REAL(SiKi)                                   :: SqrtRnMinusRm  !< Value of \f$ \sqrt{R_n} - \sqrt{R_m} \f$

         REAL(SiKi)                                   :: Den            !< Denominator
         REAL(SiKi)                                   :: Num1           !< Numerator first  term
         REAL(SiKi)                                   :: Num2           !< Numerator second term


            ! If n == m, D^- is set to zero.  It should be set to the limit as n -> m.
         IF ( n==m ) THEN
            TransFuncD_minus  =  0.0_SiKi
         ELSE

            k_nm  = k_nm_minus(n,m,k_n,k_m)

               ! Calculate R_nm that appears in multiple places
            SqrtRnMinusRm  = SQRT(R_n) - SQRT(R_m)

               ! Calculate the two pieces of the numerator
            Num1  = SqrtRnMinusRm * ( SQRT(R_m) * ( k_n*k_n - R_n*R_n ) - SQRT(R_n) * ( k_m*k_m - R_m*R_m ) )

            Num2  = 2*SqrtRnMinusRm*SqrtRnMinusRm*( k_n * k_m * COS( D2R_S*InitInp%WaveDirArr(n) - D2R_S*InitInp%WaveDirArr(m) ) + R_n*R_m )

               ! Calculate the denominator
            Den   = SqrtRnMinusRm*SqrtRnMinusRm - k_nm * tanh( k_nm * InitInp%WtrDpth )

            TransFuncD_minus  = (Num1+Num2) / Den

         ENDIF

         RETURN
      END FUNCTION TransFuncD_minus



      !> This function calculates the term \f$ D^+_{nm} \f$ used in finding the transfer functions.
      !! The equation is given by:
      !!
      !! \f$ {D}_{nm}^{+} =
      !!          \frac {   \left(\sqrt{R_n} + \sqrt{R_m}\right) \left[ \sqrt{R_m} \left( k_n^2 - R_n^2 \right) + \sqrt{R_n} \left( k_m^{2} - R_m^2 \right) \right]
      !!             ~+~  2 \left(\sqrt{R_n} + \sqrt{R_m}\right)^2
      !!                    \left[ \left|\vec{k_n}\right| \left|\vec{k_m}\right| \cos \left( \theta_n-\theta_m \right)  - R_n R_m  \right] }
      !!                {   \left(\sqrt{R_n} + \sqrt{R_m}\right)^2   -  k_{nm}^{+} \tanh \left( k_{nm}^{+} h \right)   }  \f$
      !!
      !! where \f$ k_{nm}^{+} \f$ is handled by the function ::k_nm_plus and \f$R_n\f$ is given by
      !!    \f$ R_n = \left| \overrightarrow{k}_n \right| \tanh\left(\left| \overrightarrow{k}_n \right| h \right) \f$
      !! where \f$ h \f$ is the depth from MSL (or the still water line).
      !!
      !! To calculate this, we simplify some of the common pieces:
      !!
      !! \f$ {D}_{nm}^{+}  = \frac  {     R_{nm}         \left[ \sqrt{R_m}  \left( k_n^2   - R_n^2 \right)  -  \sqrt{R_n} \left( k_m^2 - R_m^2 \right) \right]
      !!                           ~+~ 2  R_{nm}^2       \left[ \left|\vec{k_n}\right| \left|\vec{k_m}\right| \cos \left( \theta_n-\theta_m \right)  + R_n R_m  \right] }
      !!                            {     R_{nm}^2    -  k_{nm}^{+} \tanh \left( k_{nm}^+ h \right)   }  \f$
      !!
      !! where \f$   R_{nm} \equiv \sqrt{R_n} - \sqrt{R_m} \f$.
      !!
      FUNCTION TransFuncD_plus(n,m,k_n,k_m,R_n,R_m)

            ! Passed variables
         INTEGER(IntKi),               INTENT(IN   )  :: n              !< Index to the first  frequency we are dealing with
         INTEGER(IntKi),               INTENT(IN   )  :: m              !< Index to the second frequency we are dealing with
         REAL(SiKi),                   INTENT(IN   )  :: k_n            !< WaveNumber for Omega_n -- note no direction associated with this
         REAL(SiKi),                   INTENT(IN   )  :: k_m            !< WaveNumber for Omega_m -- note no direction associated with this
         REAL(SiKi),                   INTENT(IN   )  :: R_n            !< Effect scaling relationship of depth and wavenumber
         REAL(SiKi),                   INTENT(IN   )  :: R_m            !< Effect scaling relationship of depth and wavenumber
         REAL(SiKi)                                   :: TransFuncD_plus  !< Resulting \f$ D^{+}_{nm} \f$ value.

            ! Local variables
         REAL(SiKi)                                   :: k_nm           !< Value of \f$ k_{nm}^{+} \f$
         REAL(SiKi)                                   :: SqrtRnPlusRm  !< Value of \f$ \sqrt{R_n} + \sqrt{R_m} \f$

         REAL(SiKi)                                   :: Den            !< Denominator
         REAL(SiKi)                                   :: Num1           !< Numerator first  term
         REAL(SiKi)                                   :: Num2           !< Numerator second term



         k_nm  = k_nm_plus(n,m,k_n,k_m)

            ! Calculate R_nm that appears in multiple places
         SqrtRnPlusRm  = SQRT(R_n) + SQRT(R_m)

            ! Calculate the two pieces of the numerator
         Num1  = SqrtRnPlusRm * ( SQRT(R_m) * ( k_n*k_n - R_n*R_n ) + SQRT(R_n) * ( k_m*k_m - R_m*R_m ) )

         Num2  = 2*SqrtRnPlusRm*SqrtRnPlusRm*( k_n * k_m * COS( D2R_S*InitInp%WaveDirArr(n) - D2R_S*InitInp%WaveDirArr(m) ) - R_n*R_m )

            ! Calculate the denominator
         Den   = SqrtRnPlusRm*SqrtRnPlusRm - k_nm * tanh( k_nm * InitInp%WtrDpth )

         TransFuncD_plus  = (Num1+Num2) / Den


         RETURN
      END FUNCTION TransFuncD_plus






      !> This function calculates the amplitude of the combined WaveNumber, \f$ k^-_{nm} \f$ of the wave numbers
      !! for \f$ k_n \f$ and \f$ k_m \f$ for the difference frequency.  The equation is given by
      !! \f$ {k}_{nm}^{-} = \sqrt{{k_n}^2 +{k_m}^2 - 2{k_n}{k_m}\cos(\theta_n-\theta_m)} \f$
      !!
      !! @note \f$ \theta_n \f$ is given by _InitInp\%WaveDirArr(n)_
      FUNCTION k_nm_minus(n,m,k_n,k_m)

            ! Passed variables
         INTEGER(IntKi),               INTENT(IN   )  :: n              !< Index to the first  frequency we are dealing with
         INTEGER(IntKi),               INTENT(IN   )  :: m              !< Index to the second frequency we are dealing with
         REAL(SiKi),                   INTENT(IN   )  :: k_n            !< WaveNumber for \f$\omega_n\f$ -- note the direction is found in _InitInp\%WaveDirArr(n)_
         REAL(SiKi),                   INTENT(IN   )  :: k_m            !< WaveNumber for \f$\omega_m\f$ -- note the direction is found in _InitInp\%WaveDirArr(m)_
         REAL(SiKi)                                   :: k_nm_minus

         IF (n == m ) THEN
            k_nm_minus = 0.0_SiKi            ! This is just to eliminate any numerical error
         ELSE
               !bjj: added abs() because we were getting very small negative numbers here (which should be 0). 
            k_nm_minus = sqrt( abs( k_n * k_n + k_m * k_m - 2 * k_n * k_m * cos( D2R_S*InitInp%WaveDirArr(n) - D2R_S*InitINp%WaveDirArr(m) )  ) )
         ENDIF

      END FUNCTION k_nm_minus




      !> This function calculates the amplitude of the combined WaveNumber, \f$ k^+_{nm} \f$ of the wave numbers
      !! for \f$ k_n \f$ and \f$ k_m \f$ for the difference frequency.  The equation is given by
      !! \f$ {k}_{nm}^{+} = \sqrt{{k_n}^2 +{k_m}^2 + 2{k_n}{k_m}\cos(\theta_n-\theta_m)} \f$
      !!
      !! @note \f$ \theta_n \f$ is given by _InitInp\%WaveDirArr(n)_
      FUNCTION k_nm_plus(n,m,k_n,k_m)

            ! Passed variables
         INTEGER(IntKi),               INTENT(IN   )  :: n              !< Index to the first  frequency we are dealing with
         INTEGER(IntKi),               INTENT(IN   )  :: m              !< Index to the second frequency we are dealing with
         REAL(SiKi),                   INTENT(IN   )  :: k_n            !< WaveNumber for \f$\omega_n\f$ -- note the direction is found in _InitInp\%WaveDirArr(n)_
         REAL(SiKi),                   INTENT(IN   )  :: k_m            !< WaveNumber for \f$\omega_m\f$ -- note the direction is found in _InitInp\%WaveDirArr(m)_
         REAL(SiKi)                                   :: k_nm_plus

         IF (n == m ) THEN
            k_nm_plus = 2.0_SiKi * k_n       ! This is just to eliminate any numerical error.
         ELSE
            k_nm_plus = sqrt( k_n * k_n + k_m * k_m + 2_SiKi * k_n * k_m * cos( D2R_S*InitInp%WaveDirArr(n) - D2R_S*InitINp%WaveDirArr(m) )  )
         ENDIF

      END FUNCTION k_nm_plus









      SUBROUTINE CleanUp()

         CALL  ExitFFT(FFT_Data, ErrStatTmp)

         IF (ALLOCATED(TmpTimeSeries))    DEALLOCATE(TmpTimeSeries,     STAT=ErrStatTmp)

         IF (ALLOCATED(WaveVel2xCDiff))   DEALLOCATE(WaveVel2xCDiff,    STAT=ErrStatTmp)
         IF (ALLOCATED(WaveVel2yCDiff))   DEALLOCATE(WaveVel2yCDiff,    STAT=ErrStatTmp)
         IF (ALLOCATED(WaveVel2zCDiff))   DEALLOCATE(WaveVel2zCDiff,    STAT=ErrStatTmp)
         IF (ALLOCATED(WaveAcc2xCDiff))   DEALLOCATE(WaveAcc2xCDiff,    STAT=ErrStatTmp)
         IF (ALLOCATED(WaveAcc2yCDiff))   DEALLOCATE(WaveAcc2yCDiff,    STAT=ErrStatTmp)
         IF (ALLOCATED(WaveAcc2zCDiff))   DEALLOCATE(WaveAcc2zCDiff,    STAT=ErrStatTmp)
         IF (ALLOCATED(WaveVel2xDiff))    DEALLOCATE(WaveVel2xDiff,     STAT=ErrStatTmp)
         IF (ALLOCATED(WaveVel2yDiff))    DEALLOCATE(WaveVel2yDiff,     STAT=ErrStatTmp)
         IF (ALLOCATED(WaveVel2zDiff))    DEALLOCATE(WaveVel2zDiff,     STAT=ErrStatTmp)
         IF (ALLOCATED(WaveAcc2xDiff))    DEALLOCATE(WaveAcc2xDiff,     STAT=ErrStatTmp)
         IF (ALLOCATED(WaveAcc2yDiff))    DEALLOCATE(WaveAcc2yDiff,     STAT=ErrStatTmp)
         IF (ALLOCATED(WaveAcc2zDiff))    DEALLOCATE(WaveAcc2zDiff,     STAT=ErrStatTmp)
         IF (ALLOCATED(WaveDynP2CDiff))   DEALLOCATE(WaveDynP2CDiff,    STAT=ErrStatTmp)
         IF (ALLOCATED(WaveDynP2Diff))    DEALLOCATE(WaveDynP2Diff,     STAT=ErrStatTmp)

         IF (ALLOCATED(WaveVel2xCSumT1))  DEALLOCATE(WaveVel2xCSumT1,   STAT=ErrStatTmp)
         IF (ALLOCATED(WaveVel2yCSumT1))  DEALLOCATE(WaveVel2yCSumT1,   STAT=ErrStatTmp)
         IF (ALLOCATED(WaveVel2zCSumT1))  DEALLOCATE(WaveVel2zCSumT1,   STAT=ErrStatTmp)
         IF (ALLOCATED(WaveAcc2xCSumT1))  DEALLOCATE(WaveAcc2xCSumT1,   STAT=ErrStatTmp)
         IF (ALLOCATED(WaveAcc2yCSumT1))  DEALLOCATE(WaveAcc2yCSumT1,   STAT=ErrStatTmp)
         IF (ALLOCATED(WaveAcc2zCSumT1))  DEALLOCATE(WaveAcc2zCSumT1,   STAT=ErrStatTmp)
         IF (ALLOCATED(WaveVel2xSumT1))   DEALLOCATE(WaveVel2xSumT1,    STAT=ErrStatTmp)
         IF (ALLOCATED(WaveVel2ySumT1))   DEALLOCATE(WaveVel2ySumT1,    STAT=ErrStatTmp)
         IF (ALLOCATED(WaveVel2zSumT1))   DEALLOCATE(WaveVel2zSumT1,    STAT=ErrStatTmp)
         IF (ALLOCATED(WaveAcc2xSumT1))   DEALLOCATE(WaveAcc2xSumT1,    STAT=ErrStatTmp)
         IF (ALLOCATED(WaveAcc2ySumT1))   DEALLOCATE(WaveAcc2ySumT1,    STAT=ErrStatTmp)
         IF (ALLOCATED(WaveAcc2zSumT1))   DEALLOCATE(WaveAcc2zSumT1,    STAT=ErrStatTmp)
         IF (ALLOCATED(WaveDynP2CSumT1))  DEALLOCATE(WaveDynP2CSumT1,   STAT=ErrStatTmp)
         IF (ALLOCATED(WaveDynP2SumT1))   DEALLOCATE(WaveDynP2SumT1,    STAT=ErrStatTmp)

         IF (ALLOCATED(WaveVel2xCSumT2))  DEALLOCATE(WaveVel2xCSumT2,   STAT=ErrStatTmp)
         IF (ALLOCATED(WaveVel2yCSumT2))  DEALLOCATE(WaveVel2yCSumT2,   STAT=ErrStatTmp)
         IF (ALLOCATED(WaveVel2zCSumT2))  DEALLOCATE(WaveVel2zCSumT2,   STAT=ErrStatTmp)
         IF (ALLOCATED(WaveAcc2xCSumT2))  DEALLOCATE(WaveAcc2xCSumT2,   STAT=ErrStatTmp)
         IF (ALLOCATED(WaveAcc2yCSumT2))  DEALLOCATE(WaveAcc2yCSumT2,   STAT=ErrStatTmp)
         IF (ALLOCATED(WaveAcc2zCSumT2))  DEALLOCATE(WaveAcc2zCSumT2,   STAT=ErrStatTmp)
         IF (ALLOCATED(WaveVel2xSumT2))   DEALLOCATE(WaveVel2xSumT2,    STAT=ErrStatTmp)
         IF (ALLOCATED(WaveVel2ySumT2))   DEALLOCATE(WaveVel2ySumT2,    STAT=ErrStatTmp)
         IF (ALLOCATED(WaveVel2zSumT2))   DEALLOCATE(WaveVel2zSumT2,    STAT=ErrStatTmp)
         IF (ALLOCATED(WaveAcc2xSumT2))   DEALLOCATE(WaveAcc2xSumT2,    STAT=ErrStatTmp)
         IF (ALLOCATED(WaveAcc2ySumT2))   DEALLOCATE(WaveAcc2ySumT2,    STAT=ErrStatTmp)
         IF (ALLOCATED(WaveAcc2zSumT2))   DEALLOCATE(WaveAcc2zSumT2,    STAT=ErrStatTmp)
         IF (ALLOCATED(WaveDynP2CSumT2))  DEALLOCATE(WaveDynP2CSumT2,   STAT=ErrStatTmp)
         IF (ALLOCATED(WaveDynP2SumT2))   DEALLOCATE(WaveDynP2SumT2,    STAT=ErrStatTmp)

      END SUBROUTINE CleanUp




END SUBROUTINE Waves2_Init





!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is called at the end of the simulation.  The purpose of this routine is to destroy any data that is leftover.  If
!! we don't do this, we may leave memory tied up after the simulation ends.
!! To destroy the data, we call several routines that are generated by the FAST registry, so any issues with the destroy routines
!! should be addressed by the registry.exe which generates the Waves2_Types.f90 file.
!!
SUBROUTINE Waves2_End( u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )
!..................................................................................................................................

      TYPE(Waves2_InputType),             INTENT(INOUT)  :: u              !< System inputs
      TYPE(Waves2_ParameterType),         INTENT(INOUT)  :: p              !< Parameters
      TYPE(Waves2_ContinuousStateType),   INTENT(INOUT)  :: x              !< Continuous states
      TYPE(Waves2_DiscreteStateType),     INTENT(INOUT)  :: xd             !< Discrete states
      TYPE(Waves2_ConstraintStateType),   INTENT(INOUT)  :: z              !< Constraint states
      TYPE(Waves2_OtherStateType),        INTENT(INOUT)  :: OtherState     !< Other states
      TYPE(Waves2_OutputType),            INTENT(INOUT)  :: y              !< System outputs
      TYPE(Waves2_MiscVarType),           INTENT(INOUT)  :: m              !< Misc/optimization variables
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat        !< Error status of the operation
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg         !< Error message if ErrStat /= ErrID_None



         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""


         !> Place any last minute operations or calculations here.  For Waves2, most calculations all performed
         !! during the initialization, so there are no final calculations that need to be performed.


         ! Close files here.  The Waves2 module does not open any files, so there should be nothing to close. 


         !> Destroy the input data:

      CALL Waves2_DestroyInput( u, ErrStat, ErrMsg )


         !> Destroy the parameter data:

      CALL Waves2_DestroyParam( p, ErrStat, ErrMsg )


         !> Destroy the state data:

      CALL Waves2_DestroyContState(   x,           ErrStat, ErrMsg )
      CALL Waves2_DestroyDiscState(   xd,          ErrStat, ErrMsg )
      CALL Waves2_DestroyConstrState( z,           ErrStat, ErrMsg )
      CALL Waves2_DestroyOtherState(  OtherState,  ErrStat, ErrMsg )


         !> Destroy the output data:

      CALL Waves2_DestroyOutput( y, ErrStat, ErrMsg )


END SUBROUTINE Waves2_End



!----------------------------------------------------------------------------------------------------------------------------------
!> Loose coupling routine for solving constraint states, integrating continuous states, and updating discrete states.
!> Continuous, constraint, discrete, and other states are updated to values at t + Interval.
SUBROUTINE Waves2_UpdateStates( t, n, Inputs, InputTimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
!..................................................................................................................................

      REAL(DbKi),                         INTENT(IN   )  :: t              !< Current simulation time in seconds
      INTEGER(IntKi),                     INTENT(IN   )  :: n              !< Current step of the simulation: t = n*Interval
      TYPE(Waves2_InputType),             INTENT(IN   )  :: Inputs(:)      !< Inputs at InputTimes
      REAL(DbKi),                         INTENT(IN   )  :: InputTimes(:)  !< Times in seconds associated with Inputs
      TYPE(Waves2_ParameterType),         INTENT(IN   )  :: p              !< Parameters
      TYPE(Waves2_ContinuousStateType),   INTENT(INOUT)  :: x              !< Input: Continuous states at t;
                                                                           !!Output: Continuous states at t + Interval
      TYPE(Waves2_DiscreteStateType),     INTENT(INOUT)  :: xd             !< Input: Discrete states at t;
                                                                           !!Output: Discrete states at t + Interval
      TYPE(Waves2_ConstraintStateType),   INTENT(INOUT)  :: z              !< Input: Constraint states at t;
                                                                           !!Output: Constraint states at t + Interval
      TYPE(Waves2_OtherStateType),        INTENT(INOUT)  :: OtherState     !< Input: Other states at t;
                                                                           !!Output: Other states at t + Interval
      TYPE(Waves2_MiscVarType),           INTENT(INOUT)  :: m              !< Misc/optimization variables
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat        !< Error status of the operation
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg         !< Error message if ErrStat /= ErrID_None


         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "Warning: No States to update in Waves2 module. *Waves2_UpdateStates was called*"


END SUBROUTINE Waves2_UpdateStates



!----------------------------------------------------------------------------------------------------------------------------------
!> Routine for computing outputs, used in both loose and tight coupling.
!! The Waves2 module second order wave kinematic corrections are processed at initialization and passed to other modules (such as
!! Morrison) for processing.  As a result, there is nothing that needs to be calculated by the CalcOutput routine other than the
!! WriteOutput values at each timestep.
SUBROUTINE Waves2_CalcOutput( Time, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )
!..................................................................................................................................

      REAL(DbKi),                         INTENT(IN   )  :: Time           !< Current simulation time in seconds
      TYPE(Waves2_InputType),             INTENT(IN   )  :: u              !< Inputs at Time
      TYPE(Waves2_ParameterType),         INTENT(IN   )  :: p              !< Parameters
      TYPE(Waves2_ContinuousStateType),   INTENT(IN   )  :: x              !< Continuous states at Time
      TYPE(Waves2_DiscreteStateType),     INTENT(IN   )  :: xd             !< Discrete states at Time
      TYPE(Waves2_ConstraintStateType),   INTENT(IN   )  :: z              !< Constraint states at Time
      TYPE(Waves2_OtherStateType),        INTENT(IN   )  :: OtherState     !< Other states at Time
      TYPE(Waves2_OutputType),            INTENT(INOUT)  :: y              !< Outputs computed at Time (Input only so that mesh
                                                                           !!   connectivity information does not have to be recalculated)
      TYPE(Waves2_MiscVarType),           INTENT(INOUT)  :: m              !< Misc/optimization variables
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat        !< Error status of the operation
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg         !< Error message if ErrStat /= ErrID_None



         ! Local Variables:
      INTEGER(IntKi)                                     :: I                          ! Generic index
      REAL(SiKi)                                         :: WaveElev2Temp(p%NWaveElev)
      REAL(ReKi)                                         :: AllOuts(MaxWaves2Outputs)

 

         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""




         ! Abort if the Waves2 module did not calculate anything 

      IF ( .NOT. ALLOCATED ( p%WaveElev2 ) )  RETURN
      IF ( p%NumOuts < 1 ) RETURN


      DO I=1,p%NWaveElev
         WaveElev2Temp(I)  = InterpWrappedStpReal ( REAL(Time, SiKi), p%WaveTime(:), p%WaveElev2(:,I), &
                                                     m%LastIndWave, p%NStepWave + 1       )
      ENDDO

         ! Map the calculated results into the AllOuts Array
      CALL Wvs2Out_MapOutputs(Time, y, p%NWaveElev, WaveElev2Temp, AllOuts, ErrStat, ErrMsg)



              ! Put the output data in the OutData array
      DO I = 1,p%NumOuts
         y%WriteOutput(I) = p%OutParam(I)%SignM * AllOuts( p%OutParam(I)%Indx )
      END DO



END SUBROUTINE Waves2_CalcOutput




!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is required for the FAST framework, but is not actually needed for this module.
!! In the framework, this routine calculates the derivative of the continuous states.
!! As this routine is not necessary in the Waves2 module, it simply issues a warning and returns.
!! @note A few values will be set so that compilers are happy, but nothing of value is done.
SUBROUTINE Waves2_CalcContStateDeriv( Time, u, p, x, xd, z, OtherState, m, dxdt, ErrStat, ErrMsg )
!..................................................................................................................................

      REAL(DbKi),                         INTENT(IN   )  :: Time           !< Current simulation time in seconds
      TYPE(Waves2_InputType),             INTENT(IN   )  :: u              !< Inputs at Time
      TYPE(Waves2_ParameterType),         INTENT(IN   )  :: p              !< Parameters
      TYPE(Waves2_ContinuousStateType),   INTENT(IN   )  :: x              !< Continuous states at Time
      TYPE(Waves2_DiscreteStateType),     INTENT(IN   )  :: xd             !< Discrete states at Time
      TYPE(Waves2_ConstraintStateType),   INTENT(IN   )  :: z              !< Constraint states at Time
      TYPE(Waves2_OtherStateType),        INTENT(IN   )  :: OtherState     !< Other states at Time
      TYPE(Waves2_MiscVarType),           INTENT(INOUT)  :: m              !< Misc/optimization variables
      TYPE(Waves2_ContinuousStateType),   INTENT(  OUT)  :: dxdt           !< Continuous state derivatives at Time
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat        !< Error status of the operation
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg         !< Error message if ErrStat /= ErrID_None


         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "Warning: No States to take derivative of in Waves2 module. *Waves2::CalcContStateDeriv was called.  It "// &
                  "is not necessary in the Waves2 module, so it does nothing.*"


         ! Compute the first time derivatives of the continuous states here: None to calculate, so no code here.

         ! Dummy output value for dxdt -- this is only here to prevent the compiler from complaining.
   dxdt%DummyContState = 0.0_SiKi


END SUBROUTINE Waves2_CalcContStateDeriv




!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is required for the FAST framework, but is not actually needed for this module.
!! In the framework, this routine is used to update discrete states, by
!! So, this routine will simply issue a warning and return.
SUBROUTINE Waves2_UpdateDiscState( Time, n, u, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
!..................................................................................................................................

      REAL(DbKi),                         INTENT(IN   )  :: Time           !< Current simulation time in seconds
      INTEGER(IntKi),                     INTENT(IN   )  :: n              !< Current step of the simulation: t = n*Interval
      TYPE(Waves2_InputType),             INTENT(IN   )  :: u              !< Inputs at Time
      TYPE(Waves2_ParameterType),         INTENT(IN   )  :: p              !< Parameters
      TYPE(Waves2_ContinuousStateType),   INTENT(IN   )  :: x              !< Continuous states at Time
      TYPE(Waves2_DiscreteStateType),     INTENT(INOUT)  :: xd             !< Input: Discrete states at Time;
                                                                           !!   Output: Discrete states at Time + Interval
      TYPE(Waves2_ConstraintStateType),   INTENT(IN   )  :: z              !< Constraint states at Time
      TYPE(Waves2_OtherStateType),        INTENT(IN   )  :: OtherState     !< Other states at Time
      TYPE(Waves2_MiscVarType),           INTENT(INOUT)  :: m              !< Misc/optimization variables
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat        !< Error status of the operation
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg         !< Error message if ErrStat /= ErrID_None


         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "Warning: No Discrete States to update in Waves2 module. *Waves2::UpdateDiscState was called.  It is not "// &
                  "necessary in the Waves2 module, so it does nothing.*"

         ! Code to update the discrete states would live here, but there are no discrete states to update, hence no code.


END SUBROUTINE Waves2_UpdateDiscState




!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is required for the FAST framework, but is not actually needed for this module.
!! In the framework, this is a tight coupling routine for solving for the residual of the constraint state equations
!! So, this routine will simply issue a warning and return.
!! @note A few values will be set so that compilers are happy, but nothing of value is done.
SUBROUTINE Waves2_CalcConstrStateResidual( Time, u, p, x, xd, z, OtherState, m, z_residual, ErrStat, ErrMsg )
!..................................................................................................................................

      REAL(DbKi),                         INTENT(IN   )  :: Time           !< Current simulation time in seconds
      TYPE(Waves2_InputType),             INTENT(IN   )  :: u              !< Inputs at Time
      TYPE(Waves2_ParameterType),         INTENT(IN   )  :: p              !< Parameters
      TYPE(Waves2_ContinuousStateType),   INTENT(IN   )  :: x              !< Continuous states at Time
      TYPE(Waves2_DiscreteStateType),     INTENT(IN   )  :: xd             !< Discrete states at Time
      TYPE(Waves2_ConstraintStateType),   INTENT(IN   )  :: z              !< Constraint states at Time (possibly a guess)
      TYPE(Waves2_OtherStateType),        INTENT(IN   )  :: OtherState     !< Other states at Time
      TYPE(Waves2_MiscVarType),           INTENT(INOUT)  :: m              !< Misc/optimization variables
      TYPE(Waves2_ConstraintStateType),   INTENT(  OUT)  :: z_residual     !< Residual of the constraint state equations using
                                                                           !!  the input values described above
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat        !< Error status of the operation
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg         !< Error message if ErrStat /= ErrID_None


         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "Warning: No States in Waves2 module. *Waves2::CalcConstrStateResidual was called.  It is not needed in "//&
                  "the Waves2 module, so it does nothing useful."



         ! Solve for the constraint states here: Since there are no constraint states to solve for in Waves2, there is no code here.

      z_residual%DummyConstrState = 0.0_SiKi    ! This exists just so that we can make the compiler happy.

END SUBROUTINE Waves2_CalcConstrStateResidual



!----------------------------------------------------------------------------------------------------------------------------------

END MODULE Waves2
!**********************************************************************************************************************************
