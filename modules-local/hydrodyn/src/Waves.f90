!**********************************************************************************************************************************
! The Waves and Waves_Types modules make up a template for creating user-defined calculations in the FAST Modularization 
! Framework. Wavess_Types will be auto-generated based on a description of the variables for the module.
!
! "Waves" should be replaced with the name of your module. Example: HydroDyn
! "Waves" (in Waves_*) should be replaced with the module name or an abbreviation of it. Example: HD
!..................................................................................................................................
! LICENSING
! Copyright (C) 2013  National Renewable Energy Laboratory
!
!    This file is part of Waves.
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
! File last committed: $Date$
! (File) Revision #: $Rev$
! URL: $HeadURL$
!**********************************************************************************************************************************
MODULE Waves

   USE Waves_Types   
   USE NWTC_Library
   USE NWTC_FFTPACK
      
   IMPLICIT NONE
   
   PRIVATE

!   INTEGER(IntKi), PARAMETER            :: DataFormatID = 1   ! Update this value if the data types change (used in Waves_Pack)
   TYPE(ProgDesc), PARAMETER            :: Waves_ProgDesc = ProgDesc( 'Waves', 'v1.00.02', '22-Jun-2014' )

   
      ! ..... Public Subroutines ...................................................................................................
   PUBLIC :: WavePkShpDefault                     ! Return the default value of the peak shape parameter of the incident wave spectrum
   PUBLIC :: Waves_Init                           ! Initialization routine
   PUBLIC :: Waves_End                            ! Ending routine (includes clean up)
   
   PUBLIC :: Waves_UpdateStates                   ! Loose coupling routine for solving for constraint states, integrating 
                                                    !   continuous states, and updating discrete states
   PUBLIC :: Waves_CalcOutput                     ! Routine for computing outputs
   
   PUBLIC :: Waves_CalcConstrStateResidual        ! Tight coupling routine for returning the constraint state residual
   PUBLIC :: Waves_CalcContStateDeriv             ! Tight coupling routine for computing derivatives of continuous states
   PUBLIC :: Waves_UpdateDiscState                ! Tight coupling routine for updating discrete states
      
   !PUBLIC :: Waves_JacobianPInput                 ! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete-
   !                                                 !   (Xd), and constraint-state (Z) equations all with respect to the inputs (u)
   !PUBLIC :: Waves_JacobianPContState             ! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete-
   !                                                 !   (Xd), and constraint-state (Z) equations all with respect to the continuous 
   !                                                 !   states (x)
   !PUBLIC :: Waves_JacobianPDiscState             ! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete-
   !                                                 !   (Xd), and constraint-state (Z) equations all with respect to the discrete 
   !                                                 !   states (xd)
   !PUBLIC :: Waves_JacobianPConstrState           ! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete-
   !                                                 !   (Xd), and constraint-state (Z) equations all with respect to the constraint 
   !                                                 !   states (z)
   
   PUBLIC :: InterpWrappedStpReal
   
   PRIVATE:: WheelerStretching                    ! This FUNCTION applies the principle of Wheeler stretching to (1-Forward) find the elevation where the wave kinematics are to be applied using Wheeler stretching or (2-Backword)   
   PRIVATE:: BoxMuller
   PRIVATE:: JONSWAP
   PUBLIC :: WaveNumber
   PRIVATE:: UserWaveSpctrm   
   PRIVATE:: StillWaterWaves_Init
   PRIVATE:: VariousWaves_Init
   PRIVATE:: GH_BladedWaves_Init
  ! PRIVATE:: WhiteNoiseWaves_Init

CONTAINS

!=======================================================================

   FUNCTION WavePkShpDefault ( Hs, Tp )


      ! This FUNCTION is used to return the default value of the peak shape
      ! parameter of the incident wave spectrum, conditioned on significant
      ! wave height and peak spectral period.
      !
      ! There are several different versions of the JONSWAP spectrum
      ! formula.  This version is based on the one documented in the
      ! IEC61400-3 wind turbine design standard for offshore wind turbines.



   IMPLICIT                        NONE


      ! Passed Variables:

   REAL(ReKi), INTENT(IN )      :: Hs                                              ! Significant wave height (meters)
   REAL(ReKi), INTENT(IN )      :: Tp                                              ! Peak spectral period (sec)
   REAL(ReKi)                   :: WavePkShpDefault                                ! This function = default value of the peak shape parameter of the incident wave spectrum conditioned on significant wave height and peak spectral period (-)


      ! Local Variables:

   REAL(ReKi)                   :: TpOvrSqrtHs                                     ! = Tp/SQRT(Hs) (s/SQRT(m))



      ! Compute the default peak shape parameter of the incident wave spectrum,
      !   conditioned on significant wave height and peak spectral period:

   TpOvrSqrtHs = Tp/SQRT(Hs)

   IF (     TpOvrSqrtHs <= 3.6 )  THEN
      WavePkShpDefault = 5.0
   ELSEIF ( TpOvrSqrtHs >= 5.0 )  THEN
      WavePkShpDefault = 1.0
   ELSE
      WavePkShpDefault = EXP( 5.75 - 1.15*TpOvrSqrtHs )
   END IF



   RETURN
   END FUNCTION WavePkShpDefault
   
!=======================================================================
      FUNCTION BoxMuller ( NDAmp, Phase )


         ! This FUNCTION uses the Box-Muller method to turn two uniformly
         ! distributed randoms into two unit normal randoms, which are
         ! returned as real and imaginary components.



      IMPLICIT                             NONE


         ! Passed Variables:

      COMPLEX(ReKi)                     :: BoxMuller                                  ! This function

      REAL(ReKi), INTENT(IN ), OPTIONAL :: Phase                                      ! Optional phase to override random phase (radians)

      LOGICAL,    INTENT(IN )           :: NDAmp                                      ! Flag for normally-distributed amplitudes


         ! Local Variables:

      REAL(ReKi)                   :: C1                                              ! Intermediate variable
      REAL(ReKi)                   :: C2                                              ! Intermediate variable
      REAL(ReKi)                   :: U1                                              ! First  uniformly distributed random
      REAL(ReKi)                   :: U2                                              ! Second uniformly distributed random



         ! Compute the two uniformly distributed randoms:
         ! NOTE: The first random, U1, cannot be zero else the LOG() function
         !       below will blow up; there is no restriction on the value of the
         !       second random, U2.

      U1 = 0.0
      DO WHILE ( U1 == 0.0 )
         CALL RANDOM_NUMBER(U1)
      END DO
      CALL    RANDOM_NUMBER(U2)


         ! Compute intermediate variables:

      IF ( NDAmp )  THEN            ! Normally-distributed amplitudes
         C1 = SQRT( -2.0*LOG(U1) )
      ELSE                          ! Constant amplitudes (ignore U1); therefore, C1 = SQRT( 2.0 ) = MEAN( SQRT( -2.0*LOG(U1) ) for a uniform distribution of U1 between 0 and 1
         C1 = SQRT(  2.0         )
      END IF
      
      IF ( PRESENT( Phase ) )  THEN ! Specified phase to replace random phase (ignore U2)
         C2 = Phase
      ELSE                          ! Uniformly-distributed phase
         C2 = TwoPi*U2
      END IF


         ! Compute the unit normal randoms:

      BoxMuller = CMPLX( C1*COS(C2), C1*SIN(C2) )



      RETURN
      END FUNCTION BoxMuller
!=======================================================================
      FUNCTION JONSWAP ( Omega, Hs, Tp, Gamma )


         ! This FUNCTION computes the JOint North Sea WAve Project
         ! (JONSWAP) representation of the one-sided power spectral density
         ! or wave spectrum given the frequency, Omega, peak shape
         ! parameter, Gamma, significant wave height, Hs, and peak spectral
         ! period, Tp, as inputs.  If the value of Gamma is 1.0, the
         ! Pierson-Moskowitz wave spectrum is returned.
         !
         ! There are several different versions of the JONSWAP spectrum
         ! formula.  This version is based on the one documented in the
         ! IEC61400-3 wind turbine design standard for offshore wind
         ! turbines.




      IMPLICIT                        NONE


         ! Passed Variables:

      REAL(ReKi), INTENT(IN )      :: Gamma                                           ! Peak shape parameter (-)
      REAL(ReKi), INTENT(IN )      :: Hs                                              ! Significant wave height (meters)
      REAL(ReKi)                   :: JONSWAP                                         ! This function = JONSWAP wave spectrum, S (m^2/(rad/s))
      REAL(ReKi), INTENT(IN )      :: Omega                                           ! Wave frequency (rad/s)
      REAL(ReKi), INTENT(IN )      :: Tp                                              ! Peak spectral period (sec)


         ! Local Variables:

      REAL(ReKi)                   :: Alpha                                           ! Exponent on Gamma used in the spectral formulation (-)
      REAL(ReKi)                   :: C                                               ! Normalising factor used in the spectral formulation (-)
      REAL(ReKi)                   :: f                                               ! Wave frequency (Hz)
      REAL(ReKi)                   :: fp                                              ! Peak spectral frequency (Hz)
      REAL(ReKi)                   :: fpOvrf4                                         ! (fp/f)^4
      REAL(ReKi)                   :: Sigma                                           ! Scaling factor used in the spectral formulation (-)

       REAL(ReKi)                  :: Inv2Pi   =  0.15915494

         ! Compute the JONSWAP wave spectrum, unless Omega is zero, in which case,
         !   return zero:

      IF ( Omega == 0.0 )  THEN  ! When .TRUE., the formulation below is ill-conditioned; thus, the known value of zero is returned.


         JONSWAP  = 0.0


      ELSE                       ! Omega > 0.0; forumulate the JONSWAP spectrum.


         ! Compute the wave frequency and peak spectral frequency in Hz:

         f        = Inv2Pi*Omega
         fp       = 1/Tp
         fpOvrf4  = (fp/f)**4.0


         ! Compute the normalising factor:

         C        = 1.0 - ( 0.287*LOG(GAMMA) )


         ! Compute Alpha:

         IF ( f <= fp )  THEN
            Sigma = 0.07
         ELSE
            Sigma = 0.09
         END IF

         Alpha    = EXP( ( -0.5*( ( (f/fp) - 1.0 )/Sigma )**2.0 ) )


         ! Compute the wave spectrum:

         JONSWAP  = Inv2Pi*C*( 0.3125*Hs*Hs*fpOvrf4/f )*EXP( ( -1.25*fpOvrf4 ) )*( GAMMA**Alpha )


      END IF



      RETURN
      END FUNCTION JONSWAP
      !=======================================================================
!JASON: MOVE THIS USER-DEFINED ROUTINE (UserWaveSpctrm) TO THE UserSubs.f90 OF HydroDyn WHEN THE PLATFORM LOADING FUNCTIONALITY HAS BEEN DOCUMENTED!!!!!
      SUBROUTINE UserWaveSpctrm ( Omega, WaveDir, DirRoot, WaveS1Sdd )


         ! This is a dummy routine for holding the place of a user-specified
         ! wave spectrum.  Modify this code to create your own spectrum.



      IMPLICIT                        NONE


         ! Passed Variables:

      REAL(ReKi), INTENT(IN )      :: Omega                                           ! Wave frequency, rad/s.
      REAL(ReKi), INTENT(IN )      :: WaveDir                                         ! Incident wave propagation heading direction, degrees
      REAL(ReKi), INTENT(OUT)      :: WaveS1Sdd                                       ! One-sided power spectral density of the wave spectrum per unit time for the current frequency component and heading direction, m^2/(rad/s).

      CHARACTER(1024), INTENT(IN ) :: DirRoot                                         ! The name of the root file including the full path to the current working directory.  This may be useful if you want this routine to write a permanent record of what it does to be stored with the simulation results: the results should be stored in a file whose name (including path) is generated by appending any suitable extension to DirRoot.



      WaveS1Sdd = 0.0



      RETURN
      END SUBROUTINE UserWaveSpctrm
      !=======================================================================
      FUNCTION WaveNumber ( Omega, g, h )


         ! This FUNCTION solves the finite depth dispersion relationship:
         !
         !                   k*tanh(k*h)=(Omega^2)/g
         !
         ! for k, the wavenumber (WaveNumber) given the frequency, Omega,
         ! gravitational constant, g, and water depth, h, as inputs.  A
         ! high order initial guess is used in conjunction with a quadratic
         ! Newton's method for the solution with seven significant digits
         ! accuracy using only one iteration pass.  The method is due to
         ! Professor J.N. Newman of M.I.T. as found in routine EIGVAL of
         ! the SWIM-MOTION-LINES (SML) software package in source file
         ! Solve.f of the SWIM module.



      IMPLICIT                        NONE


         ! Passed Variables:

      REAL(ReKi), INTENT(IN )      :: g                                               ! Gravitational acceleration (m/s^2)
      REAL(ReKi), INTENT(IN )      :: h                                               ! Water depth (meters)
      REAL(ReKi), INTENT(IN )      :: Omega                                           ! Wave frequency (rad/s)
      REAL(ReKi)                   :: WaveNumber                                      ! This function = wavenumber, k (1/m)


         ! Local Variables:

      REAL(ReKi)                   :: A                                               ! A temporary variable used in the solution.
      REAL(ReKi)                   :: B                                               ! A temporary variable used in the solution.
      REAL(ReKi)                   :: C                                               ! A temporary variable used in the solution.
      REAL(ReKi)                   :: C2                                              ! A temporary variable used in the solution.
      REAL(ReKi)                   :: CC                                              ! A temporary variable used in the solution.
      REAL(ReKi)                   :: E2                                              ! A temporary variable used in the solution.
      REAL(ReKi)                   :: X0                                              ! A temporary variable used in the solution.



         ! Compute the wavenumber, unless Omega is zero, in which case, return
         !   zero:

      IF ( Omega == 0.0 )  THEN  ! When .TRUE., the formulation below is ill-conditioned; thus, the known value of zero is returned.


         WaveNumber = 0.0


      ELSE                       ! Omega > 0.0; solve for the wavenumber as usual.


         C  = Omega*Omega*h/g
         CC = C*C


         ! Find X0:

         IF ( C <= 2.0 )  THEN

            X0 = SQRT(C)*( 1.0 + C*( 0.169 + (0.031*C) ) )

         ELSE

            E2 = EXP(-2.0*C)

            X0 = C*( 1.0 + ( E2*( 2.0 - (12.0*E2) ) ) )

         END IF


         ! Find the WaveNumber:

         IF ( C <= 4.8 )  THEN

            C2 = CC - X0*X0
            A  = 1.0/( C - C2 )
            B  = A*( ( 0.5*LOG( ( X0 + C )/( X0 - C ) ) ) - X0 )

            WaveNumber = ( X0 - ( B*C2*( 1.0 + (A*B*C*X0) ) ) )/h

         ELSE

            WaveNumber = X0/h

         END IF


      END IF



      RETURN
      END FUNCTION WaveNumber
      
      !=======================================================================
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

      REAL(ReKi)                   :: COSHNumOvrCOSHDen                               ! This function = COSH( k*( z + h ) )/COSH( k*h ) (-)
      REAL(ReKi), INTENT(IN )      :: h                                               ! Water depth ( h      >  0 ) (meters)
      REAL(ReKi), INTENT(IN )      :: k                                               ! Wave number ( k      >= 0 ) (1/m)
      REAL(ReKi), INTENT(IN )      :: z                                               ! Elevation   (-h <= z <= 0 ) (meters)



         ! Compute the hyperbolic numerator over denominator:

      IF ( k*h  > 89.4_ReKi )  THEN   ! When .TRUE., the shallow water formulation will trigger a floating point overflow error; however, COSH( k*( z + h ) )/COSH( k*h ) = EXP( k*z ) + EXP( -k*( z + 2*h ) ) for large k*h.  This equals the deep water formulation, EXP( k*z ), except near z = -h, because h > 14.23*wavelength (since k = 2*Pi/wavelength) in this case.

         COSHNumOvrCOSHDen = EXP( k*z ) + EXP( -k*( z + 2.0_ReKi*h ) )

      ELSE                       ! 0 < k*h <= 89.4; use the shallow water formulation.

         COSHNumOvrCOSHDen = COSH( k*( z + h ) )/COSH( k*h )

      END IF



      RETURN
      END FUNCTION COSHNumOvrCOSHDen
!=======================================================================
      FUNCTION COSHNumOvrSINHDen ( k, h, z )


         ! This FUNCTION computes the shallow water hyperbolic numerator
         ! over denominator term in the wave kinematics expressions:
         !
         !                    COSH( k*( z + h ) )/SINH( k*h )
         !
         ! given the wave number, k, water depth, h, and elevation z, as
         ! inputs.



      IMPLICIT                        NONE


         ! Passed Variables:

      REAL(ReKi)                   :: COSHNumOvrSINHDen                               ! This function = COSH( k*( z + h ) )/SINH( k*h ) (-)
      REAL(ReKi), INTENT(IN )      :: h                                               ! Water depth ( h      >  0 ) (meters)
      REAL(ReKi), INTENT(IN )      :: k                                               ! Wave number ( k      >= 0 ) (1/m)
      REAL(ReKi), INTENT(IN )      :: z                                               ! Elevation   (-h <= z <= 0 ) (meters)



         ! Compute the hyperbolic numerator over denominator:


      IF (   k  < EPSILON(0.0_ReKi)  )  THEN  ! When .TRUE., the shallow water formulation is ill-conditioned; thus, HUGE(k) is returned to approximate the known value of infinity.

         COSHNumOvrSINHDen = HUGE( k )

      ELSEIF ( k*h  > 89.4_ReKi )  THEN  ! When .TRUE., the shallow water formulation will trigger a floating point overflow error; however, COSH( k*( z + h ) )/SINH( k*h ) = EXP( k*z ) + EXP( -k*( z + 2*h ) ) for large k*h.  This equals the deep water formulation, EXP( k*z ), except near z = -h, because h > 14.23*wavelength (since k = 2*Pi/wavelength) in this case.

         COSHNumOvrSINHDen = EXP( k*z ) + EXP( -k*( z + 2*h ) )

      ELSE                          ! 0 < k*h <= 89.4; use the shallow water formulation.

         COSHNumOvrSINHDen = COSH( k*( z + h ) )/SINH( k*h )

      END IF



      RETURN
      END FUNCTION COSHNumOvrSINHDen
!=======================================================================
      FUNCTION COTH ( X )


         ! This FUNCTION computes the hyperbolic cotangent,
         ! COSH(X)/SINH(X).


      USE                             Precision


      IMPLICIT                        NONE


         ! Passed Variables:

      REAL(ReKi)                   :: COTH                                            ! This function = COSH( X )/SINH( X ) (-)
      REAL(ReKi), INTENT(IN )      :: X                                               ! The argument (-)



         ! Compute the hyperbolic cotangent:

      IF ( X == 0.0 )  THEN   ! When .TRUE., the formulation below is ill-conditioned; thus, HUGE(X) is returned to approximate the known value of infinity.

         COTH = HUGE( X )

      ELSE                    ! X /= 0.0; use the numerically-stable computation of COTH(X) by means of TANH(X).

         COTH = 1.0/TANH( X ) ! = COSH( X )/SINH( X )

      END IF



      RETURN
      END FUNCTION COTH
    
      !=======================================================================
      FUNCTION SINHNumOvrSINHDen ( k, h, z )


         ! This FUNCTION computes the shallow water hyperbolic numerator
         ! over denominator term in the wave kinematics expressions:
         !
         !                    SINH( k*( z + h ) )/SINH( k*h )
         !
         ! given the wave number, k, water depth, h, and elevation z, as
         ! inputs.


      IMPLICIT                        NONE


         ! Passed Variables:

      REAL(ReKi)                   :: SINHNumOvrSINHDen                               ! This function = SINH( k*( z + h ) )/SINH( k*h ) (-)
      REAL(ReKi), INTENT(IN )      :: h                                               ! Water depth ( h      >  0 ) (meters)
      REAL(ReKi), INTENT(IN )      :: k                                               ! Wave number ( k      >= 0 ) (1/m)
      REAL(ReKi), INTENT(IN )      :: z                                               ! Elevation   (-h <= z <= 0 ) (meters)



         ! Compute the hyperbolic numerator over denominator:

      IF (     k   == 0.0  )  THEN  ! When .TRUE., the shallow water formulation is ill-conditioned; thus, the known value of unity is returned.

         SINHNumOvrSINHDen = 1.0

      ELSEIF ( k*h >  89.4_ReKi )  THEN  ! When .TRUE., the shallow water formulation will trigger a floating point overflow error; however, SINH( k*( z + h ) )/SINH( k*h ) = EXP( k*z ) - EXP( -k*( z + 2*h ) ) for large k*h.  This equals the deep water formulation, EXP( k*z ), except near z = -h, because h > 14.23*wavelength (since k = 2*Pi/wavelength) in this case.

         SINHNumOvrSINHDen = EXP( k*z ) - EXP( -k*( z + 2.0_ReKi*h ) )

      ELSE                          ! 0 < k*h <= 89.4; use the shallow water formulation.

         SINHNumOvrSINHDen = SINH( k*( z + h ) )/SINH( k*h )

      END IF



      RETURN
      END FUNCTION SINHNumOvrSINHDen
      
      
      




!----------------------------------------------------------------------------------------------------------------------------------
!SUBROUTINE WhiteNoiseWaves_Init ( InitInp, InitOut, ErrStat, ErrMsg )    
!! This routine initializes the waves data for WaveMod = 5 , or white noise waves option
!!----------------------------------------------------------------------------------------------------------------------------------
! 
!
!   TYPE(Waves_InitInputType),       INTENT(IN   )  :: InitInp     ! Input data for initialization routine    
!   TYPE(Waves_InitOutputType),      INTENT(INOUT)  :: InitOut     ! Initialization output data      
!   INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat     ! Error status of the operation
!   CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
!
!   
!      ! Local Variables
!      WaveS1Sdd = WaveDT*WaveHs*WaveHs/(8.0*Pi)
!END SUBROUTINE WhiteNoiseWaves_Init

!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE StillWaterWaves_Init ( InitInp, InitOut, ErrStat, ErrMsg )    
! This routine initializes the waves data for WaveMod = 0 , or still water waves option
!----------------------------------------------------------------------------------------------------------------------------------
 

   TYPE(Waves_InitInputType),       INTENT(IN   )  :: InitInp     ! Input data for initialization routine    
   TYPE(Waves_InitOutputType),      INTENT(INOUT)  :: InitOut     ! Initialization output data      
   INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   
      ! Local Variables
   INTEGER                      :: J                                               ! Generic index
   INTEGER(IntKi)                :: ErrStatTmp                    ! Temporary error status
   
      ! Initialize ErrStat
         
   ErrStat = ErrID_None         
   ErrStatTmp  =  ErrID_None
   ErrMsg  = "" 
   
   
   ! Initialize everything to zero:

      InitOut%NStepWave  = 2                ! We must have at least two elements in order to interpolate later on
      InitOut%NStepWave2 = 1

      ALLOCATE ( InitOut%WaveTime   (0:InitOut%NStepWave                    ) , STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array InitOut%WaveTime.',  ErrStat,ErrMsg,'StillWaterWaves_Init')

      ALLOCATE ( InitOut%WaveElevC0 (2, 0:InitOut%NStepWave2                ) , STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array InitOut%WaveElevC0.',ErrStat,ErrMsg,'StillWaterWaves_Init')

      ALLOCATE ( InitOut%WaveElev   (0:InitOut%NStepWave,InitInp%NWaveElev  ) , STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array InitOut%WaveElev.',  ErrStat,ErrMsg,'StillWaterWaves_Init')

      ALLOCATE ( InitOut%WaveDynP0  (0:InitOut%NStepWave,InitInp%NWaveKin0  ) , STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array InitOut%WaveDynP0.', ErrStat,ErrMsg,'StillWaterWaves_Init')

      ALLOCATE ( InitOut%WaveVel0   (0:InitOut%NStepWave,InitInp%NWaveKin0,3) , STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array InitOut%WaveVel0.',  ErrStat,ErrMsg,'StillWaterWaves_Init')

      ALLOCATE ( InitOut%WaveAcc0   (0:InitOut%NStepWave,InitInp%NWaveKin0,3) , STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array InitOut%WaveAcc0.',  ErrStat,ErrMsg,'StillWaterWaves_Init')

      IF ( ErrStat >= AbortErrLev ) RETURN 

      InitOut%WaveDOmega = 0.0
      InitOut%WaveTime   = (/ 0.0_DbKi, 1.0_DbKi, 2.0_DbKi /)   ! We must have at least two different time steps in the interpolation
      InitOut%WaveElevC0(1,0) = 0.0
      InitOut%WaveElevC0(2,0) = 0.0
      InitOut%WaveElev   = 0.0
      InitOut%WaveDynP0  = 0.0
      InitOut%WaveVel0   = 0.0
      InitOut%WaveAcc0   = 0.0
      
      ! Add the current velocities to the wave velocities:

      DO J = 1,InitInp%NWaveKin0   ! Loop through all Morison element nodes where the incident wave kinematics will be computed
         
         InitOut%WaveVel0(:,J,1) =  InitInp%CurrVxi(J)  ! xi-direction
         InitOut%WaveVel0(:,J,2) =  InitInp%CurrVyi(J)  ! yi-direction

      END DO                ! J - All points where the incident wave kinematics will be computed

END SUBROUTINE StillWaterWaves_Init


!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE VariousWaves_Init ( InitInp, InitOut, ErrStat, ErrMsg )
! Compute the wave kinematics and related information for  Plane progressive (regular) wave, JONSWAP/Pierson-Moskowitz spectrum  
! (irregular) wave, or user-defined spectrum (irregular) wave.
!----------------------------------------------------------------------------------------------------------------------------------

   TYPE(Waves_InitInputType),       INTENT(IN   )  :: InitInp     ! Input data for initialization routine
   TYPE(Waves_InitOutputType),       INTENT(INOUT)  :: InitOut     ! Output data      
   INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   
   
   ! Local Variables
   COMPLEX(ReKi), PARAMETER     :: ImagNmbr = (0.0,1.0)                            ! The imaginary number, SQRT(-1.0)
   COMPLEX(ReKi)                :: ImagOmega                                       ! = ImagNmbr*Omega (rad/s)
   REAL(ReKi), ALLOCATABLE      :: WaveElev0  (:)                                  ! Instantaneous elevation of incident waves at the platform reference point (meters)
   COMPLEX(ReKi), ALLOCATABLE   :: PWaveAccC0HxiPz0 (:)                            ! Partial derivative of WaveAccC0Hxi(:) with respect to zi at zi = 0 (1/s^2) 
   COMPLEX(ReKi), ALLOCATABLE   :: PWaveAccC0HyiPz0 (:)                            ! Partial derivative of WaveAccC0Hyi(:) with respect to zi at zi = 0 (1/s^2) 
   COMPLEX(ReKi), ALLOCATABLE   :: PWaveAccC0VPz0 (:)                              ! Partial derivative of WaveAccC0V  (:) with respect to zi at zi = 0 (1/s^2) 
   COMPLEX(ReKi), ALLOCATABLE   :: PWaveDynPC0BPz0(:)                              ! Partial derivative of WaveDynPC0B (:) with respect to zi at zi = 0 (N/m  ) 
   COMPLEX(ReKi), ALLOCATABLE   :: PWaveVelC0HxiPz0 (:)                            ! Partial derivative of WaveVelC0Hxi(:) with respect to zi at zi = 0 (1/s  ) 
   COMPLEX(ReKi), ALLOCATABLE   :: PWaveVelC0HyiPz0 (:)                            ! Partial derivative of WaveVelC0Hyi(:) with respect to zi at zi = 0 (1/s  ) 
   COMPLEX(ReKi), ALLOCATABLE   :: PWaveVelC0VPz0 (:)                              ! Partial derivative of WaveVelC0V  (:) with respect to zi at zi = 0 (1/s  ) 
   COMPLEX(ReKi), ALLOCATABLE   :: WaveAccC0Hxi(:,:)                               ! Discrete Fourier transform of the instantaneous horizontal acceleration in x-direction of incident waves before applying stretching at the zi-coordinates for points (m/s^2)
   COMPLEX(ReKi), ALLOCATABLE   :: WaveAccC0Hyi(:,:)                               ! Discrete Fourier transform of the instantaneous horizontal acceleration in y-direction of incident waves before applying stretching at the zi-coordinates for points (m/s^2)
   COMPLEX(ReKi), ALLOCATABLE   :: WaveAccC0V(:,:)                                 ! Discrete Fourier transform of the instantaneous vertical   acceleration                of incident waves before applying stretching at the zi-coordinates for points (m/s^2)
   COMPLEX(ReKi), ALLOCATABLE   :: WaveDynPC0(:,:)                                 ! Discrete Fourier transform of the instantaneous dynamic pressure                       of incident waves before applying stretching at the zi-coordinates for points (N/m^2)
   COMPLEX(ReKi), ALLOCATABLE   :: WaveElevC (:,:)                                 ! Discrete Fourier transform of the instantaneous elevation of incident waves (meters)
   COMPLEX(ReKi), ALLOCATABLE   :: WaveVelC0Hxi(:,:)                               ! Discrete Fourier transform of the instantaneous horizontal velocity                    of incident waves before applying stretching at the zi-coordinates for points (m/s)
   COMPLEX(ReKi), ALLOCATABLE   :: WaveVelC0Hyi(:,:)                               ! Discrete Fourier transform of the instantaneous horizontal velocity in x-direction     of incident waves before applying stretching at the zi-coordinates for points (m/s)
   COMPLEX(ReKi), ALLOCATABLE   :: WaveVelC0V(:,:)                                 ! Discrete Fourier transform of the instantaneous vertical   velocity in y-direction     of incident waves before applying stretching at the zi-coordinates for points (m/s)
   COMPLEX(ReKi)                :: WGNC                                            ! Discrete Fourier transform of the realization of a White Gaussian Noise (WGN) time series process with unit variance for the current frequency component (-)
   REAL(ReKi)                   :: CurrVxi                                         ! xi-component of the current velocity at the instantaneous elevation (m/s)
   REAL(ReKi)                   :: CurrVyi                                         ! yi-component of the current velocity at the instantaneous elevation (m/s)
   REAL(ReKi)                   :: CurrVxi0                                        ! xi-component of the current velocity at zi =  0.0 meters            (m/s)
   REAL(ReKi)                   :: CurrVyi0                                        ! yi-component of the current velocity at zi =  0.0 meters            (m/s)
   REAL(ReKi)                   :: CurrVxiS                                        ! xi-component of the current velocity at zi = -SmllNmbr meters       (m/s)
   REAL(ReKi)                   :: CurrVyiS                                        ! yi-component of the current velocity at zi = -SmllNmbr meters       (m/s)
   REAL(ReKi), ALLOCATABLE      :: CosWaveDir(:)                                   ! COS( WaveDirArr(I) ) -- Each wave frequency has a unique wave direction.
   REAL(ReKi), ALLOCATABLE      :: GHWaveAcc (:,:)                                 ! Instantaneous acceleration of incident waves in the xi- (1), yi- (2), and zi- (3) directions, respectively, at each of the GHNWvDpth vertical locations in GH Bladed wave data files (m/s^2)
   REAL(ReKi), ALLOCATABLE      :: GHWaveDynP(:  )                                 ! Instantaneous dynamic pressure of incident waves                                                            at each of the GHNWvDpth vertical locations in GH Bladed wave data files (N/m^2)
   REAL(ReKi)                   :: GHWaveTime                                      ! Instantaneous simulation times in GH Bladed wave data files (sec)
   REAL(ReKi), ALLOCATABLE      :: GHWaveVel (:,:)                                 ! Instantaneous velocity     of incident waves in the xi- (1), yi- (2), and zi- (3) directions, respectively, at each of the GHNWvDpth vertical locations in GH Bladed wave data files (m/s  )
   REAL(ReKi), ALLOCATABLE      :: GHWvDpth  (:)                                   ! Vertical locations in GH Bladed wave data files.
!UNUSED:   !REAL(ReKi), PARAMETER        :: n_Massel = 3.0                                  ! Factor used to the scale the peak spectral frequency in order to find the cut-off frequency based on the suggestion in: Massel, S. R., Ocean Surface Waves: Their Physics and Prediction, Advanced Series on Ocean Engineering - Vol. 11, World Scientific Publishing, Singapore - New Jersey - London - Hong Kong, 1996.  This reference recommends n_Massel > 3.0 (higher for higher-order wave kinemetics); the ">" designation is accounted for by checking if ( Omega > OmegaCutOff ).
   REAL(ReKi)                   :: Omega                                           ! Wave frequency (rad/s)
!UNUSED:   !REAL(ReKi)                   :: OmegaCutOff                                     ! Cut-off frequency or upper frequency limit of the wave spectrum beyond which the wave spectrum is zeroed (rad/s)
   REAL(ReKi)                   :: PCurrVxiPz0                                     ! Partial derivative of CurrVxi        with respect to zi at zi = 0 (1/s  )
   REAL(ReKi)                   :: PCurrVyiPz0                                     ! Partial derivative of CurrVyi        with respect to zi at zi = 0 (1/s  )
   REAL(ReKi), ALLOCATABLE      :: PWaveAcc0HxiPz0(:)                              ! Partial derivative of WaveAcc0Hxi(:) with respect to zi at zi = 0 (1/s^2)
   REAL(ReKi), ALLOCATABLE      :: PWaveAcc0HyiPz0(:)                              ! Partial derivative of WaveAcc0Hyi(:) with respect to zi at zi = 0 (1/s^2)
   REAL(ReKi), ALLOCATABLE      :: PWaveAcc0VPz0  (:)                              ! Partial derivative of WaveAcc0V  (:) with respect to zi at zi = 0 (1/s^2)
   REAL(ReKi), ALLOCATABLE      :: PWaveDynP0BPz0 (:)                              ! Partial derivative of WaveDynP0B (:) with respect to zi at zi = 0 (N/m  ) 
   REAL(ReKi), ALLOCATABLE      :: PWaveVel0HxiPz0(:)                              ! Partial derivative of WaveVel0Hxi(:) with respect to zi at zi = 0 (1/s  )
   REAL(ReKi), ALLOCATABLE      :: PWaveVel0HyiPz0(:)                              ! Partial derivative of WaveVel0Hyi(:) with respect to zi at zi = 0 (1/s  )
   REAL(ReKi), ALLOCATABLE      :: PWaveVel0VPz0  (:)                              ! Partial derivative of WaveVel0V  (:) with respect to zi at zi = 0 (1/s  )
   REAL(ReKi)                   :: Slope                                           ! Miscellanous slope used in an interpolation (-)
   REAL(ReKi), PARAMETER        :: SmllNmbr  = 9.999E-4                            ! A small number representing epsilon for taking numerical derivatives.  !bjj: how about using SQRT(EPSILON())?
   REAL(ReKi)                   :: SQRTNStepWave2                                  ! SQRT( NStepWave/2 )
   REAL(ReKi), ALLOCATABLE      :: SinWaveDir     (:)                              ! SIN( WaveDirArr(I) )
   REAL(ReKi), ALLOCATABLE      :: WaveAcc0Hxi (:,:)                               ! Instantaneous horizontal acceleration in x-direction of incident waves before applying stretching at the zi-coordinates for points (m/s^2)
   REAL(ReKi), ALLOCATABLE      :: WaveAcc0Hyi (:,:)                               ! Instantaneous horizontal acceleration in y-direction of incident waves before applying stretching at the zi-coordinates for points (m/s^2)
!UNUSED:  for future wave stretching
!   REAL(ReKi)                   :: WaveAcc0HxiExtrap                               ! Temporary value extrapolated from the WaveAcc0Hxi(:,:) array (m/s^2)
!   REAL(ReKi)                   :: WaveAcc0HxiExtrap                               ! Temporary value extrapolated from the WaveAcc0Hxi(:,:) array (m/s^2)
!   REAL(ReKi)                   :: WaveAcc0HyiInterp                               ! Temporary value interpolated from the WaveAcc0Hyi(:,:) array (m/s^2)
!   REAL(ReKi)                   :: WaveAcc0HyiInterp                               ! Temporary value interpolated from the WaveAcc0Hyi(:,:) array (m/s^2)
   REAL(ReKi), ALLOCATABLE      :: WaveAcc0V (:,:)                                 ! Instantaneous vertical   acceleration of incident waves before applying stretching at the zi-coordinates for points (m/s^2)
!UNUSED:  for future wave stretching
!   REAL(ReKi)                   :: WaveAcc0VExtrap                                 ! Temporary value extrapolated from the WaveAcc0V  (:,:) array (m/s^2)
!   REAL(ReKi)                   :: WaveAcc0VInterp                                 ! Temporary value interpolated from the WaveAcc0V  (:,:) array (m/s^2)
   REAL(ReKi), ALLOCATABLE      :: WaveDynP0B(:,:)                                 ! Instantaneous dynamic pressure        of incident waves before applying stretching at the zi-coordinates for points (N/m^2)
!UNUSED:  for future wave stretching
!   REAL(ReKi)                   :: WaveDynP0BExtrap                                ! Temporary value extrapolated from the WaveDynP0B (:,:) array (N/m^2)
!   REAL(ReKi)                   :: WaveDynP0BInterp                                ! Temporary value interpolated from the WaveDynP0B (:,:) array (N/m^2)
   REAL(ReKi)                   :: WaveElev_Max                                    ! Maximum expected value of the instantaneous elevation of incident waves (meters)
   REAL(ReKi)                   :: WaveElev_Min                                    ! Minimum expected value of the instantaneous elevation of incident waves (meters)
   COMPLEX(ReKi)                :: WaveElevxiPrime0
   REAL(ReKi), ALLOCATABLE      :: WaveKinzi0Prime(:)                              ! zi-coordinates for points where the incident wave kinematics will be computed before applying stretching; these are relative to the mean see level (meters)
   INTEGER   , ALLOCATABLE      :: WaveKinPrimeMap(:)
   REAL(ReKi), ALLOCATABLE      :: WaveKinzi0St   (:)                              ! Array of elevations found by stretching the elevations in the WaveKinzi0Prime(:) array using the instantaneous wave elevation; these are relative to the mean see level (meters)
   REAL(ReKi)                   :: WaveNmbr                                        ! Wavenumber of the current frequency component (1/meter)
   REAL(ReKi)                   :: WaveS1Sdd                                       ! One-sided power spectral density of the wave spectrum per unit time for the current frequency component (m^2/(rad/s))
   REAL(ReKi)                   :: WaveS2Sdd                                       ! Two-sided power spectral density of the wave spectrum per unit time for the current frequency component (m^2/(rad/s))
   REAL(DbKi)                   :: WaveTMax                                        ! Analysis time for incident wave calculations (sec)
   REAL(ReKi), ALLOCATABLE      :: WaveVel0Hxi    (:,:)                            ! Instantaneous xi-direction velocity   of incident waves before applying stretching at the zi-coordinates for points (m/s  )
!UNUSED:  for future wave stretching
!   REAL(ReKi)                   :: WaveVel0HxiExtrap                               ! Temporary value extrapolated from the WaveVel0Hxi(:,:) array (m/s  )
!   REAL(ReKi)                   :: WaveVel0HxiInterp                               ! Temporary value interpolated from the WaveVel0Hxi(:,:) array (m/s  )
   REAL(ReKi), ALLOCATABLE      :: WaveVel0Hyi    (:,:)                            ! Instantaneous yi-direction velocity   of incident waves before applying stretching at the zi-coordinates for points (m/s  )
!UNUSED:  for future wave stretching
!   REAL(ReKi)                   :: WaveVel0HyiExtrap                               ! Temporary value extrapolated from the WaveVel0Hyi(:,:) array (m/s  )
!   REAL(ReKi)                   :: WaveVel0HyiInterp                               ! Temporary value interpolated from the WaveVel0Hyi(:,:) array (m/s  )
   REAL(ReKi), ALLOCATABLE      :: WaveVel0V (:,:)                                 ! Instantaneous vertical     velocity   of incident waves before applying stretching at the zi-coordinates for points (m/s  )
!UNUSED:  for future wave stretching
!   REAL(ReKi)                   :: WaveVel0VExtrap                                 ! Temporary value extrapolated from the WaveVel0V  (:,:) array (m/s  )
!   REAL(ReKi)                   :: WaveVel0VInterp                                 ! Temporary value interpolated from the WaveVel0V  (:,:) array (m/s  )
   REAL(ReKi)                   :: zi_Max                                          ! Maximum elevation where the wave kinematics are to be applied using      stretching to the instantaneous free surface (meters)
   REAL(ReKi)                   :: zi_Min                                          ! Minimum elevation where the wave kinematics are to be applied using      stretching to the instantaneous free surface (meters)
   REAL(ReKi)                   :: ziPrime_Max                                     ! Maximum elevation where the wave kinematics are computed before applying stretching to the instantaneous free surface (meters)
   REAL(ReKi)                   :: ziPrime_Min                                     ! Minimum elevation where the wave kinematics are computed before applying stretching to the instantaneous free surface (meters)

   REAL(ReKi)                   :: WGNC_Fact
   INTEGER                      :: GHNStepWave                                     ! Total number of time steps in the GH Bladed wave data files.
   INTEGER                      :: GHNWvDpth                                       ! Number of vertical locations in GH Bladed wave data files.
   INTEGER                      :: I                                               ! Generic index
   INTEGER                      :: I_Orig                                          ! The index of the time step from original (input) part of data
   INTEGER                      :: I_WaveTp                                        ! The index of the frequency component nearest to WaveTp
   INTEGER                      :: J                                               ! Generic index
   INTEGER                      :: J_Min                                           ! The minimum value of index J such that WaveKinzi0(J) >= -WtrDpth
   INTEGER                      :: K                                               ! Generic index
   INTEGER                      :: LastInd                                         ! Index into the arrays saved from the last call as a starting point for this call
   INTEGER                      :: nSeeds                                          ! number of seeds required to initialize the intrinsic random number generator
   INTEGER                      :: NWaveKin0Prime                                  ! Number of points where the incident wave kinematics will be computed before applying stretching to the instantaneous free surface (-)
   INTEGER,    ALLOCATABLE      :: TmpWaveSeeds   (:)                              ! A temporary array used for portability. IVF/CVF use a random number generator initialized with 2 seeds; other platforms can use different implementations (e.g. gfortran needs 8 or 12 seeds)
   COMPLEX(ReKi)                :: tmpComplex                                      ! A temporary varible to hold the complex value of the wave elevation before storing it into a REAL array
   COMPLEX(ReKi),ALLOCATABLE    :: tmpComplexArr(:)                                ! A temporary array (0:NStepWave2-1) for FFT use. 
   TYPE(FFT_DataType)           :: FFT_Data                                        ! the instance of the FFT module we're using

      ! Variables for mult-direction waves
   INTEGER(IntKi)               :: WaveNDirMax                                     !< Maximum value we can change WaveNDir to (relative to original value passed in). Used in finding new WaveNDir value.
   INTEGER(IntKi)               :: WvSpreadNDir                                    !< Number of wave spreading directions for intermediate calculations.  Set later to be MAX(15*InitOut%WaveNDir,1000)
   INTEGER(IntKi)               :: WvSpreadFreqPerDir                              !< Number of wave frequencies per direction
   REAL(ReKi),    ALLOCATABLE   :: WvSpreadCos2SArr(:)                             !< Wave spreading function results array.  Used in equal energy wave spreading function.
   REAL(ReKi)                   :: WvSpreadCos2SConst                              !< Normalization constant for wave spreading function.
   REAL(ReKi),    ALLOCATABLE   :: WvSpreadIntegral(:)                             !< Cumulative integral of the wave spreading function.  Used in finding equal energy wave directions.
   REAL(ReKi)                   :: WvSpreadDTheta                                  !< Wave direction step size for intermediate calculations.  Used in finding equal energy wave directions.
   REAL(ReKi),    ALLOCATABLE   :: WvSpreadThetas(:)                               !< Wave direction used in calculations and interpolations
   REAL(ReKi),    ALLOCATABLE   :: WvSpreadThetaIdx(:)                             !< Indices for wave directions 
   REAL(ReKi),    ALLOCATABLE   :: WvTheta(:)                                      !< Final set of wave directions (degrees)
   REAL(ReKi)                   :: WvSpreadIntegralTmp                             !< Temporary variable for the interpolation

   
      ! Variables for error handling
   INTEGER(IntKi)               :: ErrStatTmp                                      !< Temporary error status
   CHARACTER(1024)              :: ErrMsgTmp                                       !< Temporary error message
   CHARACTER(1024)              :: ErrMsgTmp2                                      !< Another temporary error message



      ! Initialize ErrStat
         
   ErrStat = ErrID_None         
   ErrStatTmp  = ErrID_None
   ErrMsg  = "" 
   ErrMsgTmp   =  ""
   
      !  Set the WaveNDir information (number of wave directions).
      !     -> Since this must be adjusted later, put it in InitOut first and adjust that later in the code.
   InitOut%WaveNDir  = InitInp%WaveNDir
   InitOut%WaveDir   = InitInp%WaveDir                         ! We may want this value later (I had a nasty surprise when this wasn't set)
   WaveNDirMax       = CEILING(InitOut%WaveNDir*1.25_ReKi)     ! Value we allow WaveNDir to reach before aborting

   

!FIXME: what does this comment go to now?  It appears to have gotten separated from its source code.
   ! Plane progressive (regular) wave, JONSWAP/Pierson-Moskowitz spectrum (irregular) wave, or user-defined spectrum (irregular) wave.




         ! Tell our nice users what is about to happen that may take a while:

      CALL WrScr ( ' Generating incident wave kinematics and current time history.' )





      ! Determine the number of, NWaveKin0Prime, and the zi-coordinates for,
      !   WaveKinzi0Prime(:), points where the incident wave kinematics will be
      !   computed before applying stretching to the instantaneous free surface.
      !   The locations are relative to the mean see level.  Also determine J_Min,
      !   which is the minimum value of index J such that WaveKinzi0(J) >=
      !   -WtrDpth.  These depend on which incident wave kinematics stretching
      !   method is being used:

!JASON: ADD OTHER STRETCHING METHODS HERE, SUCH AS: DELTA STRETCHING (SEE ISO 19901-1) OR CHAKRABARTI STRETCHING (SEE OWTES)???
!JASON: APPLY STRETCHING TO THE DYNAMIC PRESSURE, IF YOU EVER COMPUTE THAT HERE!!!

      SELECT CASE ( InitInp%WaveStMod )  ! Which model are we using to extrapolate the incident wave kinematics to the instantaneous free surface?

      CASE ( 0 )                 ! None=no stretching.


      ! Since we have no stretching, NWaveKin0Prime and WaveKinzi0Prime(:) are
      !   equal to the number of, and the zi-coordinates for, the points in the
      !   WaveKinzi0(:) array between, and including, -WtrDpth and 0.0.
       
      ! Determine J_Min and NWaveKin0Prime here:

         J_Min          = 0
         NWaveKin0Prime = 0
         DO J = 1,InitInp%NWaveKin0   ! Loop through all mesh points  where the incident wave kinematics will be computed
               ! NOTE: We test to 0 instead of MSL2SWL because the locations of WaveKinzi0 and WtrDpth have already been adjusted using MSL2SWL
            IF (    InitInp%WaveKinzi0(J) >= -InitInp%WtrDpth .AND. InitInp%WaveKinzi0(J) <= 0 )  THEN
               NWaveKin0Prime = NWaveKin0Prime + 1
            END IF
         END DO                ! J - All Morison nodes where the incident wave kinematics will be computed



      ! ALLOCATE the WaveKinzi0Prime(:) array and compute its elements here:

         ALLOCATE ( WaveKinzi0Prime(NWaveKin0Prime) , STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveKinzi0Prime.',ErrStat,ErrMsg,'VariousWaves_Init')

         ALLOCATE ( WaveKinPrimeMap(NWaveKin0Prime) , STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveKinPrimeMap.',ErrStat,ErrMsg,'VariousWaves_Init')

         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF
         
         
         I = 1
         
         DO J = 1,InitInp%NWaveKin0 ! Loop through all points where the incident wave kinematics will be computed without stretching
               ! NOTE: We test to 0 instead of MSL2SWL because the locations of WaveKinzi0 and WtrDpth have already been adjusted using MSL2SWL
            IF (    InitInp%WaveKinzi0(J) >= -InitInp%WtrDpth .AND. InitInp%WaveKinzi0(J) <= 0 )  THEN
               
               WaveKinzi0Prime(I) =  InitInp%WaveKinzi0(J)
               WaveKinPrimeMap(I) =  J 
               I = I + 1
               
            END IF
            
         END DO                   ! J - All points where the incident wave kinematics will be computed without stretching



      CASE ( 1, 2 )              ! Vertical stretching or extrapolation stretching.


      ! Vertical stretching says that the wave kinematics above the mean sea level
      !   equal the wave kinematics at the mean sea level.  The wave kinematics
      !   below the mean sea level are left unchanged.
      !
      ! Extrapolation stretching uses a linear Taylor expansion of the wave
      !   kinematics (and their partial derivatives with respect to z) at the mean
      !   sea level to find the wave kinematics above the mean sea level.  The
      !   wave kinematics below the mean sea level are left unchanged.
      !
      ! Vertical stretching and extrapolation stretching do not effect the wave
      !   kinematics below the mean sea level; also, vertical stretching and
      !   extrapolation stretching say the wave kinematics above the mean sea
      !   level depend only on the mean sea level values.  Consequently,
      !   NWaveKin0Prime and WaveKinzi0Prime(:) are equal to the number of, and
      !   the zi-coordinates for, the points in the WaveKinzi0(:) array between,
      !   and including, -WtrDpth and 0.0; the WaveKinzi0Prime(:) array must also
      !   include 0.0 even if the WaveKinzi0(:) array does not.

  


      CASE ( 3 )                 ! Wheeler stretching.


      ! Wheeler stretching says that wave kinematics calculated using Airy theory
      !   at the mean sea level should actually be applied at the instantaneous
      !   free surface and that Airy wave kinematics computed at locations between
      !   the seabed and the mean sea level should be shifted vertically to new
      !   locations in proportion to their elevation above the seabed.
      !
      ! Thus, given a range of zi(:) where we want to know the wave kinematics
      !   after applying Wheeler stretching, the required range of ziPrime(:)
      !   where the wave kinematics need to be computed before applying
      !   stretching, is as follows:
      !
      ! ziPrime_Min <= ziPrime(:) <= ziPrime_Max
      !
      ! ziPrime_Min = MAX{ ( zi_Min - WaveElev_Max )/( 1 + WaveElev_Max/WtrDpth ), -WtrDpth }
      ! ziPrime_Max = MIN{ ( zi_Max - WaveElev_Min )/( 1 + WaveElev_Min/WtrDpth ),        0 }
      !
      ! where,
      !   zi_Max        = maximum elevation where the wave kinematics are to be
      !                   applied using stretching to the instantaneous free
      !                   surface
      !   zi_Min        = minimum elevation where the wave kinematics are to be
      !                   applied using stretching to the instantaneous free
      !                   surface
      !   ziPrime_Max   = maximum elevation where the wave kinematics are computed
      !                   before applying stretching to the instantaneous free
      !                   surface
      !   ziPrime_Min   = minimum elevation where the wave kinematics are computed
      !                   before applying stretching to the instantaneous free
      !                   surface
      !   WaveElev_Max  = maximum expected value of the instantaneous elevation of
      !                   incident waves
      !   WaveElev_Min  = minimum expected value of the instantaneous elevation of
      !                   incident waves
      !
      ! Thus, in order to account for Wheeler stretching when computing the wave
      !   kinematics at each of the NWaveKin0 points along a vertical line passing
      !   through the platform reference point [defined by the zi-coordinates
      !   relative to the mean see level as specified in the WaveKinzi0(:) array],
      !   we must first compute the wave kinematics without stretching at
      !   alternative elevations [indicated here by the NWaveKin0Prime-element
      !   array WaveKinzi0Prime(:)]:

   



      ENDSELECT




      ! Perform some initialization computations including initializing the
      !   pseudorandom number generator, calculating the total number of frequency
      !   components = total number of time steps in the incident wave,
      !   calculating the frequency step, calculating the index of the frequency
      !   component nearest to WaveTp, and ALLOCATing the arrays:
      ! NOTE: WaveDOmega = 2*Pi/WaveTMax since, in the FFT:
      !          Omega = (K-1)*WaveDOmega
      !          Time  = (J-1)*WaveDT
      !       and therefore:
      !          Omega*Time = (K-1)*(J-1)*WaveDOmega*WaveDT
      !                     = (K-1)*(J-1)*2*Pi/NStepWave [see NWTC_FFTPACK]
      !       or:
      !          WaveDOmega = 2*Pi/(NStepWave*WaveDT)
      !                     = 2*Pi/WaveTMax

      CALL RANDOM_SEED ( SIZE = nSeeds )
      
      IF ( nSeeds /= 2 ) THEN
         ErrMsgTmp   = ' The random number generator in use differs from the original code provided by NREL. This pRNG uses ' &
                                  //TRIM(Int2LStr(nSeeds))//' seeds instead of the 2 in the HydroDyn input file.'
         CALL SetErrStat(ErrID_Warn,ErrMsgTmp,ErrStat,ErrMsg,'VariousWaves_Init')
      END IF

      ALLOCATE ( TmpWaveSeeds ( nSeeds ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array TmpWaveSeeds.',ErrStat,ErrMsg,'VariousWaves_Init')
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL CleanUp()
         RETURN
      END IF   

         ! We'll just populate this with odd seeds = Seed(1) and even seeds = Seed(2)
      DO I = 1,nSeeds,2
         TmpWaveSeeds(I) = InitInp%WaveSeed(1)
      END DO
      DO I = 2,nSeeds,2
         TmpWaveSeeds(I) = InitInp%WaveSeed(2)
      END DO
                     
                  
      CALL RANDOM_SEED ( PUT=TmpWaveSeeds )
      DEALLOCATE(TmpWaveSeeds, STAT=ErrStatTmp)
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot deallocate array TmpWaveSeeds.',ErrStat,ErrMsg,'VariousWaves_Init')
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL CleanUp()
         RETURN
      END IF            
                  

         ! Set new value for NStepWave so that the FFT algorithms are efficient.  Note that if this method is changed, the method
         ! used to calculate the number of multidirectional wave directions (WaveNDir) may need to be updated.

      InitOut%NStepWave  = CEILING ( InitInp%WaveTMax/InitInp%WaveDT )                       ! Set NStepWave to an even integer
      IF ( MOD(InitOut%NStepWave,2) == 1 )  InitOut%NStepWave = InitOut%NStepWave + 1        !   larger or equal to WaveTMax/WaveDT.
      InitOut%NStepWave2 = MAX( InitOut%NStepWave/2, 1 )                                     ! Make sure that NStepWave is an even product of small factors (PSF) that is
      InitOut%NStepWave  = 2*PSF ( InitOut%NStepWave2, 9 )                                   !   greater or equal to WaveTMax/WaveDT to ensure that the FFT is efficient.

      InitOut%NStepWave2 = InitOut%NStepWave/2                                               ! Update the value of NStepWave2 based on the value needed for NStepWave.
      WaveTMax           = InitOut%NStepWave*InitInp%WaveDT                                  ! Update the value of WaveTMax   based on the value needed for NStepWave.
      InitOut%WaveTMax   = WaveTMax                                                          ! Update the value of WaveTMax in the output.  Needed by glue code later.
      SQRTNStepWave2 = SQRT( REAL( InitOut%NStepWave2, ReKi ) )                              ! Compute SQRT( NStepWave/2 ).
      InitOut%WaveDOmega = TwoPi/WaveTMax                                                    ! Compute the frequency step for incident wave calculations.
      I_WaveTp           = NINT ( TwoPi/(InitOut%WaveDOmega*InitInp%WaveTp) )                ! Compute the index of the frequency component nearest to WaveTp.
      


         ! Allocate all the arrays we need.

      ALLOCATE ( InitOut%WaveTime  (0:InitOut%NStepWave                    ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array InitOut%WaveTime.',  ErrStat,ErrMsg,'VariousWaves_Init')

      ALLOCATE ( InitOut%WaveElevC0(2, 0:InitOut%NStepWave2                ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array InitOut%WaveElevC0.',ErrStat,ErrMsg,'VariousWaves_Init')

      ALLOCATE ( tmpComplexArr(0:InitOut%NStepWave2                        ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array tmpComplexArr.',     ErrStat,ErrMsg,'VariousWaves_Init')

      ALLOCATE ( WaveElevC         (0:InitOut%NStepWave2 ,InitInp%NWaveElev), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveElevC.',         ErrStat,ErrMsg,'VariousWaves_Init')

      ALLOCATE ( WaveDynPC0        (0:InitOut%NStepWave2 ,NWaveKin0Prime   ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveDynPC0.',        ErrStat,ErrMsg,'VariousWaves_Init')

      ALLOCATE ( WaveVelC0Hxi      (0:InitOut%NStepWave2 ,NWaveKin0Prime   ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveVelC0Hxi.',      ErrStat,ErrMsg,'VariousWaves_Init')

      ALLOCATE ( WaveVelC0Hyi      (0:InitOut%NStepWave2 ,NWaveKin0Prime   ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveVelC0Hyi.',      ErrStat,ErrMsg,'VariousWaves_Init')

      ALLOCATE ( WaveVelC0V        (0:InitOut%NStepWave2 ,NWaveKin0Prime   ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveVelC0V.',        ErrStat,ErrMsg,'VariousWaves_Init')

      ALLOCATE ( WaveAccC0Hxi      (0:InitOut%NStepWave2 ,NWaveKin0Prime   ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveAccC0Hxi.',      ErrStat,ErrMsg,'VariousWaves_Init')

      ALLOCATE ( WaveAccC0Hyi      (0:InitOut%NStepWave2 ,NWaveKin0Prime   ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveAccC0Hyi.',      ErrStat,ErrMsg,'VariousWaves_Init')

      ALLOCATE ( WaveAccC0V        (0:InitOut%NStepWave2 ,NWaveKin0Prime   ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveAccC0V.',        ErrStat,ErrMsg,'VariousWaves_Init')

      ALLOCATE ( PWaveDynPC0BPz0   (0:InitOut%NStepWave2                   ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array PWaveDynPC0BPz0.',   ErrStat,ErrMsg,'VariousWaves_Init')

      ALLOCATE ( PWaveVelC0HxiPz0  (0:InitOut%NStepWave2                   ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array PWaveVelC0HxiPz0.',  ErrStat,ErrMsg,'VariousWaves_Init')

      ALLOCATE ( PWaveVelC0HyiPz0  (0:InitOut%NStepWave2                   ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array PWaveVelC0HyiPz0.',  ErrStat,ErrMsg,'VariousWaves_Init')

      ALLOCATE ( PWaveVelC0VPz0    (0:InitOut%NStepWave2                   ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array PWaveVelC0VPz0.',    ErrStat,ErrMsg,'VariousWaves_Init')

      ALLOCATE ( PWaveAccC0HxiPz0  (0:InitOut%NStepWave2                   ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array PWaveAccC0HxiPz0.',  ErrStat,ErrMsg,'VariousWaves_Init')

      ALLOCATE ( PWaveAccC0HyiPz0  (0:InitOut%NStepWave2                   ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array PWaveAccC0HyiPz0.',  ErrStat,ErrMsg,'VariousWaves_Init')

      ALLOCATE ( PWaveAccC0VPz0    (0:InitOut%NStepWave2                   ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array PWaveAccC0VPz0.',    ErrStat,ErrMsg,'VariousWaves_Init')

      ALLOCATE ( WaveElev0 (0:InitOut%NStepWave-1                          ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveElev0.',         ErrStat,ErrMsg,'VariousWaves_Init')

      ALLOCATE ( InitOut%WaveElev  (0:InitOut%NStepWave,InitInp%NWaveElev  ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array InitOut%WaveElev.',  ErrStat,ErrMsg,'VariousWaves_Init')

      ALLOCATE ( WaveDynP0B        (0:InitOut%NStepWave-1,NWaveKin0Prime   ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveDynP0B.',        ErrStat,ErrMsg,'VariousWaves_Init')

      ALLOCATE ( WaveVel0Hxi       (0:InitOut%NStepWave-1,NWaveKin0Prime   ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveVel0Hxi.',       ErrStat,ErrMsg,'VariousWaves_Init')

      ALLOCATE ( WaveVel0Hyi       (0:InitOut%NStepWave-1,NWaveKin0Prime   ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveVel0Hyi.',       ErrStat,ErrMsg,'VariousWaves_Init')

      ALLOCATE ( WaveVel0V         (0:InitOut%NStepWave-1,NWaveKin0Prime   ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveVel0V.',         ErrStat,ErrMsg,'VariousWaves_Init')

      ALLOCATE ( WaveAcc0Hxi       (0:InitOut%NStepWave-1,NWaveKin0Prime   ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveAcc0Hxi.',       ErrStat,ErrMsg,'VariousWaves_Init')

      ALLOCATE ( WaveAcc0Hyi       (0:InitOut%NStepWave-1,NWaveKin0Prime   ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveAcc0Hyi.',       ErrStat,ErrMsg,'VariousWaves_Init')

      ALLOCATE ( WaveAcc0V         (0:InitOut%NStepWave-1,NWaveKin0Prime   ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveAcc0V.',         ErrStat,ErrMsg,'VariousWaves_Init')

      ALLOCATE ( PWaveDynP0BPz0    (0:InitOut%NStepWave-1                  ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array PWaveDynP0BPz0.',    ErrStat,ErrMsg,'VariousWaves_Init')

      ALLOCATE ( PWaveVel0HxiPz0   (0:InitOut%NStepWave-1                  ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array PWaveVel0HxiPz0.',   ErrStat,ErrMsg,'VariousWaves_Init')

      ALLOCATE ( PWaveVel0HyiPz0   (0:InitOut%NStepWave-1                  ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array PWaveVel0HyiPz0.',   ErrStat,ErrMsg,'VariousWaves_Init')

      ALLOCATE ( PWaveVel0VPz0     (0:InitOut%NStepWave-1                  ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array PWaveVel0Pz0.',      ErrStat,ErrMsg,'VariousWaves_Init')

      ALLOCATE ( PWaveAcc0HxiPz0   (0:InitOut%NStepWave-1                  ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array PWaveAcc0HxiPz0.',   ErrStat,ErrMsg,'VariousWaves_Init')

      ALLOCATE ( PWaveAcc0HyiPz0   (0:InitOut%NStepWave-1                  ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array PWaveAcc0HyiPz0.',   ErrStat,ErrMsg,'VariousWaves_Init')

      ALLOCATE ( PWaveAcc0VPz0     (0:InitOut%NStepWave-1                  ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array PWaveAcc0VPz0.',     ErrStat,ErrMsg,'VariousWaves_Init')

      ALLOCATE ( InitOut%WaveDynP0 (0:InitOut%NStepWave,InitInp%NWaveKin0  ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array InitOut%WaveDynP0.', ErrStat,ErrMsg,'VariousWaves_Init')

      ALLOCATE ( InitOut%WaveVel0  (0:InitOut%NStepWave,InitInp%NWaveKin0,3), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array InitOut%WaveVel0.',  ErrStat,ErrMsg,'VariousWaves_Init')

      ALLOCATE ( InitOut%WaveAcc0  (0:InitOut%NStepWave,InitInp%NWaveKin0,3), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array InitOut%WaveAcc0.',  ErrStat,ErrMsg,'VariousWaves_Init')

         ! Wave direction associated with each frequency
      ALLOCATE ( InitOut%WaveDirArr( 0:InitOut%NStepWave2                  ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array InitOut%WaveDirArr.',ErrStat,ErrMsg,'VariousWaves_Init')

         ! Arrays for the Sin and Cos of the wave direction for each frequency.  Used in calculating wave elevation, velocity, acceleration etc.
      ALLOCATE ( CosWaveDir( 0:InitOut%NStepWave2                          ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array CosWaveDir.',        ErrStat,ErrMsg,'VariousWaves_Init')

      ALLOCATE ( SinWaveDir( 0:InitOut%NStepWave2                          ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array SinWaveDir.',        ErrStat,ErrMsg,'VariousWaves_Init')


         ! Now check if all the allocations worked properly
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL CleanUp()
         RETURN
      END IF


         ! For creating animations of the sea surface, the WaveElevXY array is passed in with a series of x,y coordinates
         ! (index 1).  The second index corresponds to the number of points passed in.  A two dimensional time series
         ! is created with the first index corresponding to the timestep, and second index corresponding to the second
         ! index of the WaveElevXY array.
      IF ( ALLOCATED(InitInp%WaveElevXY)) THEN
         ALLOCATE ( InitOut%WaveElevSeries (0:InitOut%NStepWave, 1:SIZE(InitInp%WaveElevXY, DIM=2)) , STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) THEN
            CALL SetErrStat(ErrID_Fatal,'Cannot allocate array InitOut%WaveElevSeries.',ErrStat,ErrMsg,'VariousWaves_Init')
            CALL CleanUp()
            RETURN
         END IF
      ENDIF

      
!FIXME: Is this piece still needed?  If so, why is it commented out?
     ! Calculate the factors needed by the discrete time inverse Fourier
      !   transform in the calculations of the White Gaussian Noise (WGN) and
      !   the two-sided power spectral density of the wave spectrum per unit time:

         ! This factor is needed by the discrete time inverse Fourier transform to ensure that the time series WGN
         ! process has unit variance
    !  WGNC_Fact = SQRT( Pi/(InitOut%WaveDOmega*InitInp%WaveDT) )





      !--------------------------------------------------------------------------------
      !> #  Multi Directional Waves
      !> ## Adjust WaveNDir
      !!
      !! If multi-directional waves will be used, the value of WaveNDir may need to be adjusted.  The reason is that
      !! for the equal energy approach used here, the following condition must be met:
      !!
      !!       CONDITION:  (NStepWave2) / WaveNDir    must be an integer
      !!
      !! If this is true, then an equal number of frequencies is assigned to each of the WaveNDir directions which
      !! gives the proper wave direction distribution function.  Otherwise, the energy distribution by direction
      !! will not be correct.
      !!
      !! _WaveNDir_ could not be adjusted before _NStepWave2_ was finalized above.
      !!
      !! @note    Use the value of WaveNDir stored in InitOut since InitInp cannot be changed.
      !!
      !! @note    Originally, the criteria had been that (NStepWave2 - 1) / WaveNDir is an integer.  This criteria
      !!          was relaxed by setting the direction for Omega = 0 (which has no amplitude) since it was found that
      !!          (NStepWave2 - 1) is often a prime number due to how NStepWave is calculated above to be a product
      !!          of smallish numbers.
      
      IF ( InitInp%WaveMultiDir ) THEN    ! Multi-directional waves in use

            ! Check that the number of wave directions is a positive odd number.  In theory this has been
            ! done before the Waves module was called.  We repeat it here in the event that the Waves module
            ! gets used in some other code.
            !  -> If it is less than 0, error out.
            !  -> If it is even, we will increment it by 1.
         IF ( InitOut%WaveNDir <= 0_IntKi ) THEN
            ErrMsgTmp   = 'WaveNDir must be an odd number greater than 0.'
            ErrStatTmp  = ErrID_Fatal
            CALL SetErrStat(ErrStatTmp,ErrMsgTmp,ErrStat,ErrMsg,'VariousWaves_Init')
            RETURN
         END IF

         IF ( MODULO( InitOut%WaveNDir, 2_IntKi) == 0_IntKi ) THEN
            InitOut%WaveNDir  = InitOut%WaveNDir + 1
            ErrMsgTmp = 'WaveNDir must be odd.  Changing the value to '//TRIM(Num2LStr(InitOut%WaveNDir))
            CALL SetErrStat(ErrID_Warn,ErrMsgTmp,ErrStat,ErrMsg,'VariousWaves_Init')
         END IF

            ! Now adjust WaveNDir as necessary so that (NStepWave2) / WaveNDir is integer
         IF ( .NOT. EqualRealNos(REAL((     InitOut%NStepWave2 )/     InitOut%WaveNDir), &
                                     ((REAL(InitOut%NStepWave2))/REAL(InitOut%WaveNDir)) )) THEN
            DO WHILE ( InitOut%WaveNDir <= WaveNDirMax )

               InitOut%WaveNDir = InitOut%WaveNDir + 2.0_ReKi
               IF ( EqualRealNos(REAL((     InitOut%NStepWave2 )/     InitOut%WaveNDir), &
                                     ((REAL(InitOut%NStepWave2))/REAL(InitOut%WaveNDir)) )) THEN
                  ErrMsgTmp   =  'Changed WaveNDir from '//TRIM(Num2LStr(InitInp%WaveNDir))//' to '//  &
                                 TRIM(Num2LStr(InitOut%WaveNDir))//' so that an equal number of frequencies are assigned to '// &
                                 'each direction.'
                  CALL SetErrStat(ErrID_Warn,ErrMsgTmp,ErrStat,ErrMsg,'VariousWaves_Init')
                  EXIT
               END IF
            END DO
         END IF

            ! If we exited because we hit a limit (in which case the condition is not satisfied), then we cannot continue.
            ! We warn the user that a value for WaveNDir was not found, and that they should try a different value, or try
            ! a different value for WaveTMax.  The reason for suggesting the latter is that NStepWave is derived from
            ! WaveTMax and adjusted until it is a product of smallish numbers (most likely even, but not necessarily so).
            ! So, there is a very small possibility then that NStepWave2 is a prime number, in which case we won't find a
            ! value for WaveNDir, so we suggest that the user change WaveTMax.  To make this a little easier for the user,
            ! we will report the first 5 possible values for WaveNDir between their requested value and 1/4 of NStepWave2,
            ! if there are any.
         IF ( .NOT.  EqualRealNos(REAL((     InitOut%NStepWave2 )/     InitOut%WaveNDir), &
                                      ((REAL(InitOut%NStepWave2))/REAL(InitOut%WaveNDir)) )) THEN
            ErrMsgTmp   = 'Could not find value for WaveNDir between '//TRIM(Num2LStr(InitInp%WaveNDir))//' and '// &
                           TRIM(Num2LStr(WaveNDirMax))//' such that an equal number of frequencies are assigned to each '// &
                           'direction.'
            ErrStatTmp  = ErrID_Fatal

            ! Now check for the possible values of WaveNDir so that we can tell the user about it.  The variable 'I' contains
            ! the count of the number of values of WaveNDir found.
            I = 0
            ErrMsgTmp2 = 'The next values of WaveNDir that work with the selected values for WaveTMax and WaveDT:'
            DO WHILE ( InitOut%WaveNDir <= INT(InitOut%NStepWave2/4.0) )
               IF ( EqualRealNos(REAL((     InitOut%NStepWave2 )/     InitOut%WaveNDir), &
                                     ((REAL(InitOut%NStepWave2))/REAL(InitOut%WaveNDir)) )) THEN
                  ErrMsgTmp2  = TRIM(ErrMsgTmp2)//"  "//TRIM(Num2LStr(InitOut%WaveNDir))
                  I = I + 1
               END IF

               InitOut%WaveNDir = InitOut%WaveNDir + 2.0_ReKi

               IF ( I >= 5 ) EXIT

            END DO
 
            ! If there were no additional values for WaveNDir found, I will be 0, so we rewrite the error message.
            IF ( I == 0 ) THEN
               ErrMsgTmp2  =  'There are no values for WaveNDir between '//TRIM(Num2LStr(WaveNDirMax))//' and '// &
                              TRIM(Num2LStr(INT(InitOut%NStepWave2/4.0)))//' (4 frequencies per wave direction)'// &
                              ' that will work with the selected values for WaveTMax ('//TRIM(Num2Lstr(InitOut%WaveTMax))// &
                              ') and WaveDT ('//TRIM(Num2LStr(InitInp%WaveDT))//').  Change either'// &
                              ' WaveTMax or WaveDT.'
            ELSE
               ErrMsgTmp2  = TRIM(ErrMsgTmp2)//'.'
            ENDIF
            
            ! Append the message about the possible values for WaveNDir (if any were found) and set the error status before
            ! returning to the calling program.
            ErrMsgTmp   =  TRIM(ErrMsgTmp)//NewLine//' '//TRIM(ErrMsgTmp2)

            CALL SetErrStat(ErrStatTmp,ErrMsgTmp,ErrStat,ErrMsg,'VariousWaves_Init')
            RETURN
         END IF

            ! Save the number of frequencies per direction so that we can use it later in assigning the directios.
         WvSpreadFreqPerDir   =  (InitOut%NStepWave2)/InitOut%WaveNDir


         !> ## Calculate the wave directions based on an equal energy approach.
         !!
         !! All the angles are supplied in degrees and are converted as needed.  For the cosine function,
         !! we could convert degrees to radians, but the conversion constant cancels out.
         !!
         !! |  Variable          |  Fortran Name         |  Location |  Units       |   Description                                          |
         !! | :----------------: | :-------------------: | :-------: | :--------:   | :----------------------------------------------------- |
         !! | \f$\bar\theta\f$   | _WaveDir_             | _InitInp_ |  (degrees)   |  Mean direction heading (_WaveDir_)                    |
         !! | \f$\Theta\f$       | _WaveNDir_            | _InitOut_ |  (-)         |  Number of wave directions                             |
         !! | \f$\delta\theta\f$ | _WaveDirRange_        | _InitInp_ |  (degrees)   |  Full range of spreading function                      |
         !! | \f$S\f$            | _WaveDirSpread_       | _InitInp_ |  (-)         |  The spreading coefficient                             |
         !! |                    | _WvSpreadNDir_        |  local    |  (-)         |  Number of angles discretizing the spreading function  |
         !! | \f$C\f$            | _WvSpreadCos2SConst_  |  local    |  (1/degrees) |  The normalization coefficient                         |
         !! |                    | _WvTheta_             |  local    |  (degrees)   |  The interpolated wave directions to assign to         |
         !! | \f$\theta_i\f$     | _WvSpreadThetas_      |  local    |  (degrees)   |  Array of wave directions associated with _WvSpreadIntegral_ |
         !! |                    | _D2R_                 |  global   | (rad/degree) |  Constant from library to convert degrees to radians   |
         !!
         !! The equal energy approach is used to set the wave directions such that each direction has the same
         !! number of frequencies.  To ensure that direction spreading function (Cosine^2S in this case) has
         !! the correct overal energy distribution shape, the wave directions are adjusted.  The spacing between
         !! directions is closer near the central direction than in the tails of the spreading function.  The
         !! method distributes the wave directions so that the energy integral between wave directions is kept
         !! constant.  The following steps are taken:
         !!
         !! 1. Discretize the spreading function over the range _InitInp%WaveDirRange_ into _WvSpreadNDir_.
         !!
         !! 2. Calculate the spreading function, _WvSpreadCos2SArr_, in the range.\n
         !!          \f$ D(\theta) = C \left| \cos\left(\frac{\pi (\theta-\bar\theta)}{\delta\theta}\right)\right|^{2S} \f$\n
         !!       where\n
         !!          \f$ C = \frac{\sqrt{\pi} \: \Gamma(S+1)}{\delta\theta \: \Gamma(S+1/2)} \f$,
         !!       and
         !!          \f$ \Gamma \f$ is the gamma function.
         !!
         !! 3. Calculate the integral of WvSpreadCos2SArr up to the current angle, and save it as
         !!       WvSpreadIntegral. The integral can be written as:\n
         !!          \f$P(\theta) = \int\limits^{\theta}_{\bar\theta - \delta\theta/2} D(\theta') \: \mathrm{d}\theta'\f$
         !!
         !! 4. Do a sanity check on the result of \f$P(\theta)\f$ over the range.
         !!
         !! 5. Divide the integrated area of _WvSpreadCos2SArr_ into _InitOut%WaveNDir_ directions (the final number
         !!       of wave directions that was solved for above).  To do this, simply find the _1/WaveNDir_ values
         !!       of the integral and interpolate to find the values of the _WvSpreadThetas_ that match.  These are the
         !!       new wave directions to use.  These results are stored in the array _WvTheta_.
         !!
         !! 6. Cleanup
         !!

         !> ### Code Implimentation order
         !! 1. Discretize the spreading function range and calculate the values of the wave spreading function

            ! Now that we have the value for _WaveNDir_ found above, we set the value of _WvSpreadNDir_ to be 15x as
            ! large, or 1000 (whichever is larger).  WvSpreadNDir is used only in discretization for later
            ! interpolation of actual wave directions.
         WvSpreadNDir   = MAX(15*InitOut%WaveNDir,1000)
         WvSpreadDTheta = InitInp%WaveDirRange/REAL(WvSpreadNDir,ReKi)

            ! Calculate the normalization constant for the wave spreading.
         WvSpreadCos2SConst   = sqrt(Pi)* (NWTC_GAMMA(InitInp%WaveDirSpread + 1.0_ReKi))/               &
                                (InitInp%WaveDirRange * NWTC_GAMMA(InitInp%WaveDirSpread + 0.5_ReKi))

            ! Allocate arrays to use for storing the intermediate values
         ALLOCATE( WvSpreadCos2SArr(0:WvSpreadNDir),  STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WvSpreadCos2SArr.',  ErrStat,ErrMsg,'VariousWaves_Init')

         ALLOCATE( WvSpreadIntegral(0:WvSpreadNDir),  STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WvSpreadIntegral.',  ErrStat,ErrMsg,'VariousWaves_Init')
      
         ALLOCATE( WvSpreadThetas(0:WvSpreadNDir),    STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WvSpreadThetas.',    ErrStat,ErrMsg,'VariousWaves_Init')

         ALLOCATE( WvTheta(1:InitOut%WaveNDir),       STAT=ErrStatTmp )    ! Dealocate this at very end of routine.
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WvTheta.',           ErrStat,ErrMsg,'VariousWaves_Init')

         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF

            !> 2. Calculate the spreading function as a function of angle.  Step through all _WvSpreadNDir_ steps.
         DO I=0,WvSpreadNDir
               ! The current angle as we step through the range
            WvSpreadThetas(I) =  I*WvSpreadDTheta  + InitInp%WaveDir - InitInp%WaveDirRange/(2.0_ReKi)

               ! Calculate the wave spreading for the current value of WvSpreadThetas
            WvSpreadCos2SArr(I)  =  WvSpreadCos2SConst*abs( cos(Pi*(WvSpreadThetas(I)-InitInp%WaveDir)/InitInp%WaveDirRange) ) &
                                                            **(2*InitInp%WaveDirSpread)

            !> 3. Calculate the integral of the spreading function up to the current angle and save it.
            !     Remember that the first element can't refer to one before it.
            IF (I == 0) THEN
               WvSpreadIntegral(I)  =  WvSpreadCos2SArr(I) * WvSpreadDTheta
            ELSE
               WvSpreadIntegral(I)  =  WvSpreadCos2SArr(I) * WvSpreadDTheta + WvSpreadIntegral(I-1)
            END IF
         ENDDO

      
            !> 4. Perform a quick sanity check.  The last value of the integral table should be 1.0 exactly.
            !!    We will allow for a 1% deviation.  If for some reason an error occurs, it may be due to the
            !!    GAMMA function calculation for the normalization constant, _WvSpreadCos2SConst_.
         IF ( WvSpreadIntegral(WvSpreadNDir) < 0.99_ReKi .OR. WvSpreadIntegral(WvSpreadNDir) > 1.01_ReKi ) THEN
            CALL SetErrStat(ErrID_Fatal,' Something went wrong in evaluating the multidirectional wave spreading function.  '// &
                           'Integral is '//TRIM(Num2LStr(WvSpreadIntegral(WvSpreadNDir))),ErrStat,ErrMsg,'VariousWaves_Init')
            RETURN
         END IF


            !> 5. Set the wave directions using the results from the integral.
            !  We will use the variable LastInd as a simple index for figuring out where in the array we are.  First set to 0
         LastInd  =  0_IntKi
         DO I=1,InitOut%WaveNDir
            WvSpreadIntegralTmp  = (REAL(I)-0.5_ReKi)/REAL(InitOut%WaveNDir)
            WvTheta(I)    =  InterpStp( WvSpreadIntegralTmp, WvSpreadIntegral, WvSpreadThetas, LastInd, WvSpreadNDir )
         ENDDO       ! I=1,InitOut%WaveNDir

            ! Store the minimum and maximum wave directions
         InitOut%WaveDirMin   = MINVAL(WvTheta)
         InitOut%WaveDirMax   = MAXVAL(WvTheta)

            !> 6. Done with equal energy wavedirection calculations.  Deallocate the arrays used during calculations.
         IF(ALLOCATED( WvSpreadCos2SArr ))      DEALLOCATE( WvSpreadCos2SArr, STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot deallocate array WvSpreadCos2SArr.',  ErrStat,ErrMsg,'VariousWaves_Init')
         IF(ALLOCATED( WvSpreadIntegral ))      DEALLOCATE( WvSpreadIntegral, STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot deallocate array WvSpreadIntegral.',  ErrStat,ErrMsg,'VariousWaves_Init')
         IF(ALLOCATED( WvSpreadThetas   ))      DEALLOCATE( WvSpreadThetas,   STAT=ErrStatTmp   )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot deallocate array WvSpreadThetas.',  ErrStat,ErrMsg,'VariousWaves_Init')

         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF

      ELSE     ! Multi-directional waves not used
      
         InitOut%WaveDirMin   = InitInp%WaveDir
         InitOut%WaveDirMax   = InitInp%WaveDir

      ENDIF    ! Multi-directional waves in use (InitInp%WaveMultiDir == .TRUE.)






      

      ! JASON: IMPLEMENT EQUATIONS (2.12 - 2.13) IN MY DISSERTATION SO THAT ONE CAN READ IN EXTERNAL WAVE
      !        DATA?<--BETTER YET, IMPLEMENT WaveElevC0 = DFT(WaveElev) WHERE WaveElev CAN BE READ IN AS
      !        GH BLADED WAVE DATA.  THAT IS, ADD AN OPTION TO READ IN WAVE DATA FOR FLOATERS!

      ! Compute the positive-frequency components (including zero) of the discrete
      !   Fourier transforms of the wave kinematics:

      DO I = 0,InitOut%NStepWave2  ! Loop through the positive frequency components (including zero) of the discrete Fourier transforms


      ! Compute the frequency of this component and its imaginary value:

             Omega = I*       InitOut%WaveDOmega
         ImagOmega = ImagNmbr*Omega



      ! Compute the discrete Fourier transform of the realization of a White
      !   Gaussian Noise (WGN) time series process with unit variance:
      !
      ! NOTE: For the time series process to be real with zero mean, the values at
      !       Omega == 0.0 and Omega == NStepWave2*WaveDOmega (= WaveOmegaMax)
      !       must be zero.

         IF ( ( I == 0 ) .OR. ( I == InitOut%NStepWave2 ) )  THEN   ! .TRUE. if ( Omega == 0.0 ) or ( Omega == NStepWave2*WaveDOmega (= WaveOmegaMax) )
            WGNC = (0.0,0.0)
         ELSEIF ( InitInp%WaveMod == 10 )  THEN                     ! .TRUE. for plane progressive (regular) waves with a specified phase
            WGNC = BoxMuller ( InitInp%WaveNDAmp, InitInp%WavePhase )
               ! This scaling of WGNC is used to ensure that the Box-Muller method is only providing a random phase,
               ! not a magnitude change, at the frequency of the plane progressive wave.  The SQRT(2.0) is used to
               ! ensure that the time series WGN process has unit variance (i.e. sinusoidal with amplitude SQRT(2.0)).
               !
               ! NOTE: the denominator here will never equal zero since U1 cannot equal 1.0, and thus, C1 cannot be
               !       0.0 in the Box-Muller method.
            IF (  ( I == I_WaveTp ) )  WGNC = WGNC*( SQRT(2.0)/ABS(WGNC) )
         ELSE                                               ! All other Omega
            WGNC = BoxMuller ( InitInp%WaveNDAmp                    )
               ! This scaling of WGNC is used to ensure that the Box-Muller method is only providing a random phase,
               ! not a magnitude change, at the frequency of the plane progressive wave.  The SQRT(2.0) is used to
               ! ensure that the time series WGN process has unit variance (i.e. sinusoidal with amplitude SQRT(2.0)).
               !
               ! NOTE: the denominator here will never equal zero since U1 cannot equal 1.0, and thus, C1 cannot be
               !       0.0 in the Box-Muller method.
            IF ( ( InitInp%WaveMod == 1 ) .AND. ( I == I_WaveTp ) )  WGNC = WGNC*( SQRT(2.0)/ABS(WGNC) )
         END IF


      ! Compute the one-sided power spectral density of the wave spectrum per unit
      !   time; zero-out the wave spectrum above the cut-off frequency:

         SELECT CASE ( InitInp%WaveMod ) ! Which incident wave kinematics model are we using?

         CASE ( 1, 10 )          ! Plane progressive (regular) wave; the wave spectrum is an impulse function centered on frequency component closest to WaveTp.
            IF ( I == I_WaveTp )  THEN       ! .TRUE. if we are at the Omega closest to WaveTp.
               WaveS1Sdd = 0.5*(InitInp%WaveHs/2.0)*(InitInp%WaveHs/2.0)/InitOut%WaveDOmega
            ELSE                             ! All other Omega
               WaveS1Sdd = 0.0
            END IF

         CASE ( 2 )              ! JONSWAP/Pierson-Moskowitz spectrum (irregular) wave.
               !  Zero-out the wave spectrum above the cut-off frequency.  We must cut-off the frequency in order to
               !  void nonphysical wave forces.  Waves that have wavelengths much smaller than the platform diameter
               !  (high frequency) do not contribute to the net force because regions of positive and negative
               !  velocity/acceleration are experienced by the platform at the same time and cancel out.
               !
               !  JASON: OTHER FREQUENCY CUT-OFF CONDITIONS ARE USED THROUGHOUT THE INDUSTRY.  SHOULD YOU USE ONE OF
               !         THEM INSTEAD?  SEE, FOR EXAMPLE, MY E-MAIL EXCHANGES WITH PAUL SCLAVOUNOS DATED 5/26/2006 OR:
               !         "GH Bladed Thoery Manual" OR: Trumars, Jenny M. V.; Tarp-Johansen, Niels Jacob; Krogh, Thomas;
               !         "The Effect of Wave Modelling on Offshore Wind Turbine Fatigue Loads," 2005 Copenhagen Offshore
               !         Wind Conference and Exhibition, 26-28 October 2005, Copenhagen, Denmark [CD-ROM].
            IF ( Omega < InitInp%WvLowCOff .OR. Omega > InitInp%WvHiCOff )  THEN ! .TRUE. if Omega is above or below the cut-off frequency
               WaveS1Sdd = 0.0
            ELSE                             ! All other Omega
               WaveS1Sdd = JONSWAP ( Omega, InitInp%WaveHs, InitInp%WaveTp, InitInp%WavePkShp )
            END IF
         CASE ( 3 )              ! White-noise
            IF ( Omega < InitInp%WvLowCOff .OR. Omega > InitInp%WvHiCOff )  THEN ! .TRUE. if Omega is above or below the cut-off frequency
               WaveS1Sdd = 0.0  
            ELSE  
               WaveS1Sdd =  InitInp%WaveHs * InitInp%WaveHs / ( 8.0 * (InitInp%WvHiCOff - InitInp%WvLowCOff) )
            END IF
         CASE ( 4 )              ! User-defined spectrum (irregular) wave.
            IF ( Omega < InitInp%WvLowCOff .OR. Omega > InitInp%WvHiCOff )  THEN ! .TRUE. if Omega is above or below the cut-off frequency
               WaveS1Sdd = 0.0  
            ELSE  
               CALL UserWaveSpctrm ( Omega, InitInp%WaveDir, InitInp%DirRoot, WaveS1Sdd )
            END IF

         ENDSELECT


      ! Compute the two-sided power spectral density of the wave spectrum per unit
      !   time:

         WaveS2Sdd  = 0.5*WaveS1Sdd


      ! Compute the discrete Fourier transform of the instantaneous elevation of
      !   incident waves at the WAMIT reference point:
         tmpComplex                   = SQRTNStepWave2*WGNC*SQRT( TwoPi*WaveS2Sdd/REAL(InitInp%WaveDT,ReKi) )
         InitOut%WaveElevC0     (1,I) = REAL( tmpComplex)
         InitOut%WaveElevC0     (2,I) = AIMAG(tmpComplex)

      END DO                ! I - The positive frequency components (including zero) of the discrete Fourier transforms



      !--------------------------------------------------------------------------------
      !=== Multi-Directional Waves ===
      !> ## Assign Wave directions
      !!  For the equal energy approach to the multi-directional waves, we need to use the random number generator to
      !!  select which direction each wave frequency is assigned to.  We also require that the phase and amplitudes
      !!  assigned to each frequency are the same regardless of whether or not multiple directions are used, we must
      !!  first finish assigning all the amplitudes and phases before using the random number generator again.  For this
      !!  reason, the above do loop is completed, the multiple wave directions are computed, and then we run through the
      !!  all wave frequencies again to set up the remaining pieces.  If we did not do this, we would change the seed
      !!  used by the random number generator before selecting the next amplitude and phase pair.
      !!
      !!  The wave directions are assigned in groups of _WaveNDir_ frequencies such that each frequency is assigned to
      !!  one of the _WaveNDir_ unique wave directions.  Each wave direction is used only once within each group of
      !!  frequencies.
      !!
      !!  When complete, we deallocate the _WvSpreadThetas_ array that was used to store the assigned directions.

      IF ( InitInp%WaveMultiDir .AND. InitInp%WaveNDir > 1 )   THEN     ! Multi-directional waves in use


            ! Allocate the index array for each group of frequencies.  This array is used to randomize the directions
            ! within each WaveNDir sized group of frequencies.  This is a REAL array used to hold the random numbers.
         ALLOCATE( WvSpreadThetaIdx(1:InitOut%WaveNDir), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WvSpreadThetaIdx while assigning '//  &
                                                          'wave directions.',ErrStat,ErrMsg,'VariousWaves_Init')
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF


            ! Reset the K so that we can use it to count the frequency index.
            ! It should be exactly NStepWave2 when done assigning directions. The the Omega = 0 has
            ! no amplitude, but gets a direction anyhow (to simplify the calculation of WaveNDir).
         K  = 0


            ! Work through the frequencies in groups of directions.
         DO I = 1,WvSpreadFreqPerDir

               ! Populate the array with random numbers
            CALL RANDOM_NUMBER(WvSpreadThetaIdx)

            DO J = 1, InitOut%WaveNDir

                  ! Find the index lowest value in the WvSpreadThetaIdx array.  This is the index to
                  ! use for this wave direction.
               LastInd  = MINLOC( WvSpreadThetaIdx, DIM=1 )

                  ! Assign the direction for this frequency piece to the LastInd value.
               InitOut%WaveDirArr(K)   =  WvTheta( LastInd )

                  ! Now make that element in the WvSpreadThetaIdx really big so we don't pick it again
               WvSpreadThetaIdx( LastInd )   = HUGE(1.0)

               K  = K + 1     ! Increment the frequency index

            ENDDO
         ENDDO

            ! Perform a quick sanity check.  We should have assigned all wave frequencies a direction, so K should be
            ! K = NStepWave2 (K is incrimented afterwards).
         IF ( K /= (InitOut%NStepWave2 ) )    CALL SetErrStat(ErrID_Fatal,  &
                     'Something went wrong while assigning wave directions.',ErrStat,ErrMsg,'VariousWaves_Init')
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF

            !  We are done with the indexing array, so deallocate it
         IF(ALLOCATED( WvSpreadThetaIdx ))   DEALLOCATE( WvSpreadThetaIdx )

            !  This had been allocated above (calculation section of equal energy portion of multidirectional waves).
            !  Deallocate it here after we have completed all the calculations involving it.
         IF(ALLOCATED( WvTheta ))            DEALLOCATE( WvTheta )

      ELSE     ! Not really multi-directional waves

            ! Since we do not have multi-directional waves, we must set the wave direction array to the single wave heading.
         InitOut%WaveDirArr   = InitInp%WaveDir
      ENDIF    ! Multi-directional waves in use.


         ! Set the CosWaveDir and SinWaveDir arrays
      CosWaveDir=COS(D2R*InitOut%WaveDirArr)
      SinWaveDir=SIN(D2R*InitOut%WaveDirArr)


      !--------------------------------------------------------------------------------
      !> ## Compute IFFTs
      !> Compute the discrete Fourier transform of the instantaneous elevation of
      !!   incident waves at each desired point on the still water level plane
      !!   where it can be output:

      DO I = 0,InitOut%NStepWave2  ! Loop through the positive frequency components (including zero) of the discrete Fourier transforms

      
            ! Set tmpComplex to the Ith element of the WAveElevC0 array

         tmpComplex  = CMPLX(  InitOut%WaveElevC0(1,I),   InitOut%WaveElevC0(2,I))
         

      ! Compute the frequency of this component and its imaginary value:

             Omega = I*       InitOut%WaveDOmega
         ImagOmega = ImagNmbr*Omega

      ! Compute the wavenumber:  

         WaveNmbr   = WaveNumber ( Omega, InitInp%Gravity, InitInp%WtrDpth )
         

      ! Compute the discrete Fourier transform of the incident wave kinematics
      !   before applying stretching at the zi-coordinates for the WAMIT reference point, and all
      !   points where are Morison loads will be calculated.


         DO J = 1,NWaveKin0Prime ! Loop through all points where the incident wave kinematics will be computed without stretching

            WaveElevxiPrime0 = EXP( -ImagNmbr*WaveNmbr*( InitInp%WaveKinxi0(WaveKinPrimeMap(J))*CosWaveDir(I) + &
                                                         InitInp%WaveKinyi0(WaveKinPrimeMap(J))*SinWaveDir(I) ))
            WaveDynPC0 (I,J)     = InitOut%RhoXg*tmpComplex*WaveElevxiPrime0 * COSHNumOvrCOSHDen ( WaveNmbr, InitInp%WtrDpth, WaveKinzi0Prime(J) )       
            WaveVelC0Hxi (I,J)   = CosWaveDir(I)*Omega*tmpComplex* WaveElevxiPrime0 * COSHNumOvrSINHDen ( WaveNmbr, InitInp%WtrDpth, WaveKinzi0Prime(J) )
            WaveVelC0Hyi (I,J)   = SinWaveDir(I)*Omega*tmpComplex* WaveElevxiPrime0 * COSHNumOvrSINHDen ( WaveNmbr, InitInp%WtrDpth, WaveKinzi0Prime(J) )
            WaveVelC0V (I,J)     =           ImagOmega*tmpComplex* WaveElevxiPrime0 * SINHNumOvrSINHDen ( WaveNmbr, InitInp%WtrDpth, WaveKinzi0Prime(J) )
            WaveAccC0Hxi (I,J)   = ImagOmega*        WaveVelC0Hxi (I,J)
            WaveAccC0Hyi (I,J)   = ImagOmega*        WaveVelC0Hyi (I,J)
            WaveAccC0V (I,J)     = ImagOmega*        WaveVelC0V   (I,J)

         END DO                   ! J - All points where the incident wave kinematics will be computed without stretching

         PWaveDynPC0BPz0(I )  = InitOut%RhoXg*tmpComplex*WaveNmbr*TANH ( WaveNmbr*InitInp%WtrDpth )
         PWaveVelC0HxiPz0(I ) = CosWaveDir(I)*Omega*tmpComplex*WaveNmbr
         PWaveVelC0HyiPz0(I ) = SinWaveDir(I)*Omega*tmpComplex*WaveNmbr
         PWaveVelC0VPz0(I )   =           ImagOmega*tmpComplex*WaveNmbr*COTH ( WaveNmbr*InitInp%WtrDpth )
         PWaveAccC0HxiPz0(I ) =      ImagOmega*       PWaveVelC0HxiPz0(I  )
         PWaveAccC0HyiPz0(I ) =      ImagOmega*       PWaveVelC0HyiPz0(I  )
         PWaveAccC0VPz0(I )   =      ImagOmega*       PWaveVelC0VPz0(I  )


      END DO                ! I - The positive frequency components (including zero) of the discrete Fourier transforms
      


      ! Calculate the array of simulation times at which the instantaneous
      !   elevation of, velocity of, acceleration of, and loads associated with
      !   the incident waves are to be determined:

      DO I = 0,InitOut%NStepWave ! Loop through all time steps
         InitOut%WaveTime(I) = I*REAL(InitInp%WaveDT,ReKi)
      END DO                ! I - All time steps

      DO I = 0,InitOut%NStepWave2  ! Loop through the positive frequency components (including zero) of the discrete Fourier transform
         tmpComplexArr(I)    =  CMPLX(InitOut%WaveElevC0(1,I), InitOut%WaveElevC0(2,I))
      END DO
      
      ! Compute the inverse discrete Fourier transforms to find the time-domain
      !   representations of the wave kinematics without stretcing:

      CALL InitFFT ( InitOut%NStepWave, FFT_Data, .TRUE., ErrStatTmp )
      CALL SetErrStat(ErrStatTmp,'Error occured while initializing the FFT.',ErrStat,ErrMsg,'VariousWaves_Init')
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL CleanUp()
         RETURN
      END IF
      
         ! We'll need the following for wave stretching once we implement it.
      CALL    ApplyFFT_cx (  WaveElev0    (:),  tmpComplexArr    (:  ), FFT_Data, ErrStatTmp )
      CALL SetErrStat(ErrStatTmp,'Error occured while applying the FFT to WaveElev0.',ErrStat,ErrMsg,'VariousWaves_Init')
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL CleanUp()
         RETURN
      END IF

      DO J = 1,InitInp%NWaveElev      ! Loop through all points where the incident wave elevations can be output
               ! This subroutine call applies the FFT at the correct location.
         CALL WaveElevTimeSeriesAtXY( InitInp%WaveElevxi(J), InitInp%WaveElevyi(J), InitOut%WaveElev(:,J), ErrStatTmp, ErrMsgTmp )
         CALL SetErrStat(ErrStatTmp,'Error occured while applying the FFT to InitOut%WaveElev.',ErrStat,ErrMsg,'VariousWaves_Init')
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF
      END DO                   ! J - All points where the incident wave elevations can be output


         ! Calculate the wave elevation at all points requested in the array WaveElevXY
      IF ( ALLOCATED(InitInp%WaveElevXY) ) THEN
         DO J = 1,SIZE(InitInp%WaveElevXY, DIM=2)
               ! This subroutine call applies the FFT at the correct location.
            CALL WaveElevTimeSeriesAtXY( InitInp%WaveElevXY(1,J),  InitInp%WaveElevXY(2,J),   InitOut%WaveElevSeries(0:InitOut%NStepWave,J), ErrStatTmp, ErrMsgTmp )
            CALL SetErrStat(ErrStatTmp,ErrMsgTmp,ErrStat,ErrMsg,'VariousWaves_Init')
            IF ( ErrStat >= AbortErrLev ) THEN
               CALL CleanUp()
               RETURN
            END IF
         ENDDO
      ENDIF

      
         ! User requested data points -- Do all the FFT calls first, then return if something failed.
      DO J = 1,NWaveKin0Prime ! Loop through all points where the incident wave kinematics will be computed without stretching
         CALL ApplyFFT_cx (          WaveDynP0B   (:,J),          WaveDynPC0    (:,J), FFT_Data, ErrStatTmp )
         CALL SetErrStat(ErrStatTmp,'Error occured while applying the FFT to WaveDynP0B.',       ErrStat,ErrMsg,'VariousWaves_Init')

         CALL ApplyFFT_cx (          WaveVel0Hxi  (:,J),          WaveVelC0Hxi  (:,J), FFT_Data, ErrStatTmp )
         CALL SetErrStat(ErrStatTmp,'Error occured while applying the FFT to WaveVel0Hxi.',      ErrStat,ErrMsg,'VariousWaves_Init')

         CALL ApplyFFT_cx (          WaveVel0Hyi  (:,J),          WaveVelC0Hyi  (:,J), FFT_Data, ErrStatTmp )
         CALL SetErrStat(ErrStatTmp,'Error occured while applying the FFT to WaveVel0Hyi.',      ErrStat,ErrMsg,'VariousWaves_Init')

         CALL ApplyFFT_cx (          WaveVel0V    (:,J),          WaveVelC0V    (:,J), FFT_Data, ErrStatTmp )
         CALL SetErrStat(ErrStatTmp,'Error occured while applying the FFT to WaveVel0V.',        ErrStat,ErrMsg,'VariousWaves_Init')

         CALL ApplyFFT_cx (          WaveAcc0Hxi  (:,J),          WaveAccC0Hxi  (:,J), FFT_Data, ErrStatTmp )
         CALL SetErrStat(ErrStatTmp,'Error occured while applying the FFT to WaveAcc0Hxi.',      ErrStat,ErrMsg,'VariousWaves_Init')

         CALL ApplyFFT_cx (          WaveAcc0Hyi  (:,J),          WaveAccC0Hyi  (:,J), FFT_Data, ErrStatTmp )
         CALL SetErrStat(ErrStatTmp,'Error occured while applying the FFT to WaveAcc0Hyi.',      ErrStat,ErrMsg,'VariousWaves_Init')

         CALL ApplyFFT_cx (          WaveAcc0V    (:,J),          WaveAccC0V    (:,J), FFT_Data, ErrStatTmp )
         CALL SetErrStat(ErrStatTmp,'Error occured while applying the FFT to WaveAcc0V.',        ErrStat,ErrMsg,'VariousWaves_Init')

         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF

      END DO                   ! J - All points where the incident wave kinematics will be computed without stretching

         ! FFT's of the partial derivatives
      CALL  ApplyFFT_cx (         PWaveDynP0BPz0(:  ),         PWaveDynPC0BPz0(:  ), FFT_Data, ErrStatTmp )
      CALL  SetErrStat(ErrStatTmp,'Error occured while applying the FFT to PWaveDynP0BPz0.',   ErrStat,ErrMsg,'VariousWaves_Init')

      CALL  ApplyFFT_cx (         PWaveVel0HxiPz0 (:  ),       PWaveVelC0HxiPz0( : ),FFT_Data, ErrStatTmp )
      CALL  SetErrStat(ErrStatTmp,'Error occured while applying the FFT to PWaveVel0HxiPz0.',  ErrStat,ErrMsg,'VariousWaves_Init')

      CALL  ApplyFFT_cx (         PWaveVel0HyiPz0 (:  ),       PWaveVelC0HyiPz0( : ),FFT_Data, ErrStatTmp )
      CALL  SetErrStat(ErrStatTmp,'Error occured while applying the FFT to PWaveVel0HyiPz0.',  ErrStat,ErrMsg,'VariousWaves_Init')

      CALL  ApplyFFT_cx (         PWaveVel0VPz0 (:  ),         PWaveVelC0VPz0 (:  ), FFT_Data, ErrStatTmp )
      CALL  SetErrStat(ErrStatTmp,'Error occured while applying the FFT to PWaveVel0VPz0.',    ErrStat,ErrMsg,'VariousWaves_Init')

      CALL  ApplyFFT_cx (         PWaveAcc0HxiPz0 (:  ),       PWaveAccC0HxiPz0(:  ),FFT_Data, ErrStatTmp )
      CALL  SetErrStat(ErrStatTmp,'Error occured while applying the FFT to PWaveAcc0HxiPz0.',  ErrStat,ErrMsg,'VariousWaves_Init')

      CALL  ApplyFFT_cx (         PWaveAcc0HyiPz0 (:  ),       PWaveAccC0HyiPz0(:  ),FFT_Data, ErrStatTmp )
      CALL  SetErrStat(ErrStatTmp,'Error occured while applying the FFT to PWaveAcc0HyiPz0.',  ErrStat,ErrMsg,'VariousWaves_Init')

      CALL  ApplyFFT_cx (         PWaveAcc0VPz0 (:  ),         PWaveAccC0VPz0( :  ), FFT_Data, ErrStatTmp )
      CALL  SetErrStat(ErrStatTmp,'Error occured while applying the FFT to PWaveAcc0VPz0.',    ErrStat,ErrMsg,'VariousWaves_Init')

      IF ( ErrStat >= AbortErrLev ) THEN
         CALL CleanUp()
         RETURN
      END IF

      CALL  ExitFFT(FFT_Data, ErrStatTmp)
      CALL  SetErrStat(ErrStatTmp,'Error occured while cleaning up after the FFTs.', ErrStat,ErrMsg,'VariousWaves_Init')
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL CleanUp()
         RETURN
      END IF



      ! Add the current velocities to the wave velocities:
      ! NOTE: Both the horizontal velocities and the partial derivative of the
      !       horizontal velocities with respect to zi at zi = 0 are found here.
      !
      ! NOTE:  The current module must be called prior to the waves module.  If that was not done, then we 
      !        don't have a current to add to the wave velocity.  So, check if the current velocity components
      !        exist.


      ! If there is a current, we need to add that (the current module was called prior to calling this module

      IF(ALLOCATED(InitInp%CurrVxi)) THEN

         DO J = 1,NWaveKin0Prime ! Loop through all points where the incident wave kinematics will be computed without stretching

            WaveVel0Hxi (:,J) =  WaveVel0Hxi (:,J) +  InitInp%CurrVxi(WaveKinPrimeMap(J))     ! xi-direction
            WaveVel0Hyi (:,J) =  WaveVel0Hyi (:,J) +  InitInp%CurrVyi(WaveKinPrimeMap(J))     ! yi-direction

         END DO                   ! J - All points where the incident wave kinematics will be computed without stretching

         PWaveVel0HxiPz0(:  ) =  PWaveVel0HxiPz0(:  ) + InitInp%PCurrVxiPz0  ! xi-direction
         PWaveVel0HyiPz0(:  ) =  PWaveVel0HyiPz0(:  ) + InitInp%PCurrVyiPz0  ! yi-direction

      ENDIF


      ! Apply stretching to obtain the wave kinematics, WaveDynP0, WaveVel0, and
      !   WaveAcc0, at the desired locations from the wave kinematics at
      !   alternative locations, WaveDynP0B, WaveVel0Hxi, WaveVel0Hyi, WaveVel0V,
      !   WaveAcc0Hxi, WaveAcc0Hyi, WaveAcc0V, if the elevation of the point defined by
      !   WaveKinzi0(J) lies between the seabed and the instantaneous free
      !   surface, else set WaveDynP0, WaveVel0, and WaveAcc0 to zero.  This
      !   depends on which incident wave kinematics stretching method is being
      !   used:

      SELECT CASE ( InitInp%WaveStMod )  ! Which model are we using to extrapolate the incident wave kinematics to the instantaneous free surface?

      CASE ( 0 )                 ! None=no stretching.


      ! Since we have no stretching, the wave kinematics between the seabed and
      !   the mean sea level are left unchanged; below the seabed or above the
      !   mean sea level, the wave kinematics are zero:
         
         DO I = 0,InitOut%NStepWave-1       ! Loop through all time steps
            K = 1
            DO J = 1,InitInp%NWaveKin0      ! Loop through all points where the incident wave kinematics will be computed
                  ! NOTE: We test to 0 instead of MSL2SWL because the locations of WaveKinzi0 and WtrDpth have already been adjusted using MSL2SWL
               IF (   ( InitInp%WaveKinzi0(J) < -InitInp%WtrDpth ) .OR. ( InitInp%WaveKinzi0(J) > 0.0          ) )  THEN   ! .TRUE. if the elevation of the point defined by WaveKinzi0(J) lies below the seabed or above mean sea level (exclusive)

                  InitOut%WaveDynP0(I,J  )  = 0.0
                  InitOut%WaveVel0 (I,J,:)  = 0.0
                  InitOut%WaveAcc0 (I,J,:)  = 0.0

               ELSE                                                                                 ! The elevation of the point defined by WaveKinzi0(J) must lie between the seabed and the mean sea level (inclusive)

                  InitOut%WaveDynP0(I,J  )  = WaveDynP0B (I,K     )
                  InitOut%WaveVel0 (I,J,1)  = WaveVel0Hxi(I,K     )
                  InitOut%WaveVel0 (I,J,2)  = WaveVel0Hyi(I,K     )
                  InitOut%WaveVel0 (I,J,3)  = WaveVel0V  (I,K     )
                  InitOut%WaveAcc0 (I,J,1)  = WaveAcc0Hxi(I,K     )
                  InitOut%WaveAcc0 (I,J,2)  = WaveAcc0Hyi(I,K     )
                  InitOut%WaveAcc0 (I,J,3)  = WaveAcc0V  (I,K     )
                  K = K + 1
               END IF

            END DO                   ! J - All points where the incident wave kinematics will be computed

         END DO                      ! I - All time steps



      CASE ( 1 )                 ! Vertical stretching.


      ! Vertical stretching says that the wave kinematics above the mean sea level
      !   equal the wave kinematics at the mean sea level.  The wave kinematics
      !   below the mean sea level are left unchanged:

   



      CASE ( 2 )                 ! Extrapolation stretching.


      ! Extrapolation stretching uses a linear Taylor expansion of the wave
      !   kinematics (and their partial derivatives with respect to z) at the mean
      !   sea level to find the wave kinematics above the mean sea level.  The
      !   wave kinematics below the mean sea level are left unchanged:

        



      CASE ( 3 )                 ! Wheeler stretching.


      ! Wheeler stretching says that wave kinematics calculated using Airy theory
      !   at the mean sea level should actually be applied at the instantaneous
      !   free surface and that Airy wave kinematics computed at locations between
      !   the seabed and the mean sea level should be shifted vertically to new
      !   locations in proportion to their elevation above the seabed.
      !
      ! Computing the wave kinematics with Wheeler stretching requires that first
      !   say that the wave kinematics we computed at the elevations defined by
      !   the WaveKinzi0Prime(:) array are actual applied at the elevations found
      !   by stretching the elevations in the WaveKinzi0Prime(:) array using the
      !   instantaneous wave elevation--these new elevations are stored in the
      !   WaveKinzi0St(:) array.  Next, we interpolate the wave kinematics
      !   computed without stretching to the desired elevations (defined in the
      !   WaveKinzi0(:) array) using the WaveKinzi0St(:) array:

      


      ENDSELECT

      ! Set the ending timestep to the same as the first timestep
      InitOut%WaveDynP0(InitOut%NStepWave,:  )  = InitOut%WaveDynP0(0,:  )
      InitOut%WaveVel0 (InitOut%NStepWave,:,:)  = InitOut%WaveVel0 (0,:,:)
      InitOut%WaveAcc0 (InitOut%NStepWave,:,:)  = InitOut%WaveAcc0 (0,:,:)




   CALL CleanUp ( )
   

CONTAINS


   SUBROUTINE WaveElevTimeSeriesAtXY(Xcoord,Ycoord, WaveElevSeriesAtXY, ErrStatLcl, ErrMsgLcl )

      REAL(ReKi),       INTENT(IN   )                 :: Xcoord
      REAL(ReKi),       INTENT(IN   )                 :: Ycoord
      REAL(ReKi),       INTENT(  OUT)                 :: WaveElevSeriesAtXY(0:InitOut%NStepWave)
      INTEGER(IntKi),   INTENT(  OUT)                 :: ErrStatLcl
      INTEGER(IntKi)                                  :: ErrStatLcl2
      CHARACTER(1024),  INTENT(  OUT)                 :: ErrMsgLcl

      ! This is probably poor programming practice, but we will use I, Omega, WaveNmbr, and tmpComplexArr from the calling routine here.

      ErrStatLcl  = ErrID_None

         ! Zero out the temporary array.
      tmpComplexArr  = CMPLX(0.0_ReKi,0.0_ReKi)

         ! Loop through the positive frequency components (including zero).  Skip the last point since that is zero by definition.
      DO I = 0,InitOut%NStepWave2-1

         Omega             = I*       InitOut%WaveDOmega
         WaveNmbr          = WaveNumber ( Omega, InitInp%Gravity, InitInp%WtrDpth )
         tmpComplexArr(I)  =  CMPLX(  InitOut%WaveElevC0(1,I),   InitOut%WaveElevC0(2,I))   *          &
                  EXP( -ImagNmbr*WaveNmbr*(  Xcoord*CosWaveDir(I)+    &
                                             Ycoord*SinWaveDir(I) )   )
      ENDDO

      CALL ApplyFFT_cx (   WaveElevSeriesAtXY(0:InitOut%NStepWave-1),   tmpComplexArr, FFT_Data,   ErrStatLcl2  )
      CALL SetErrStat(ErrStatLcl2,'Error occured while applying the FFT to InitOut%WaveElevSeries.',ErrStatLcl,ErrMsgLcl,'WaveElevTimeSeriesAtXY')

         ! Append first datpoint as the last as aid for repeated wave data
      WaveElevSeriesAtXY(InitOut%NStepWave) = WaveElevSeriesAtXY(0)

   END SUBROUTINE WaveElevTimeSeriesAtXY



   SUBROUTINE CleanUp( )

      IF (ALLOCATED( WaveKinPrimeMap ))   DEALLOCATE( WaveKinPrimeMap,  STAT=ErrStatTmp)
      IF (ALLOCATED( WaveKinzi0Prime ))   DEALLOCATE( WaveKinzi0Prime,  STAT=ErrStatTmp)
      IF (ALLOCATED( WvSpreadCos2SArr ))  DEALLOCATE( WvSpreadCos2SArr, STAT=ErrStatTmp)
      IF (ALLOCATED( WvSpreadIntegral ))  DEALLOCATE( WvSpreadIntegral, STAT=ErrStatTmp)
      IF (ALLOCATED( WvSpreadThetas ))    DEALLOCATE( WvSpreadThetas,   STAT=ErrStatTmp)
      IF (ALLOCATED( TmpWaveSeeds ))      DEALLOCATE( TmpWaveSeeds,     STAT=ErrStatTmp)
      IF (ALLOCATED( GHWaveAcc ))         DEALLOCATE( GHWaveAcc,        STAT=ErrStatTmp)
      IF (ALLOCATED( GHWaveDynP ))        DEALLOCATE( GHWaveDynP,       STAT=ErrStatTmp)
      IF (ALLOCATED( GHWaveVel ))         DEALLOCATE( GHWaveVel,        STAT=ErrStatTmp)
      IF (ALLOCATED( GHWvDpth ))          DEALLOCATE( GHWvDpth,         STAT=ErrStatTmp)
      IF (ALLOCATED( PWaveAcc0HxiPz0 ))   DEALLOCATE( PWaveAcc0HxiPz0,  STAT=ErrStatTmp)
      IF (ALLOCATED( PWaveAcc0HyiPz0 ))   DEALLOCATE( PWaveAcc0HyiPz0,  STAT=ErrStatTmp)
      IF (ALLOCATED( PWaveAcc0VPz0 ))     DEALLOCATE( PWaveAcc0VPz0,    STAT=ErrStatTmp)
      IF (ALLOCATED( PWaveAccC0HxiPz0 ))  DEALLOCATE( PWaveAccC0HxiPz0, STAT=ErrStatTmp)
      IF (ALLOCATED( PWaveAccC0HyiPz0 ))  DEALLOCATE( PWaveAccC0HyiPz0, STAT=ErrStatTmp)
      IF (ALLOCATED( PWaveAccC0VPz0 ))    DEALLOCATE( PWaveAccC0VPz0,   STAT=ErrStatTmp)
      IF (ALLOCATED( PWaveDynP0BPz0 ))    DEALLOCATE( PWaveDynP0BPz0,   STAT=ErrStatTmp)
      IF (ALLOCATED( PWaveDynPC0BPz0 ))   DEALLOCATE( PWaveDynPC0BPz0,  STAT=ErrStatTmp)
      IF (ALLOCATED( PWaveVel0HxiPz0 ))   DEALLOCATE( PWaveVel0HxiPz0,  STAT=ErrStatTmp)
      IF (ALLOCATED( PWaveVel0HyiPz0 ))   DEALLOCATE( PWaveVel0HyiPz0,  STAT=ErrStatTmp)
      IF (ALLOCATED( PWaveVel0VPz0 ))     DEALLOCATE( PWaveVel0VPz0,    STAT=ErrStatTmp)
      IF (ALLOCATED( PWaveVelC0HxiPz0 ))  DEALLOCATE( PWaveVelC0HxiPz0, STAT=ErrStatTmp)
      IF (ALLOCATED( PWaveVelC0HyiPz0 ))  DEALLOCATE( PWaveVelC0HyiPz0, STAT=ErrStatTmp)
      IF (ALLOCATED( PWaveVelC0VPz0 ))    DEALLOCATE( PWaveVelC0VPz0,   STAT=ErrStatTmp)
      IF (ALLOCATED( WaveAcc0Hxi ))       DEALLOCATE( WaveAcc0Hxi,      STAT=ErrStatTmp)
      IF (ALLOCATED( WaveAcc0Hyi ))       DEALLOCATE( WaveAcc0Hyi,      STAT=ErrStatTmp)
      IF (ALLOCATED( WaveAcc0V ))         DEALLOCATE( WaveAcc0V,        STAT=ErrStatTmp)
      IF (ALLOCATED( WaveAccC0Hxi ))      DEALLOCATE( WaveAccC0Hxi,     STAT=ErrStatTmp)
      IF (ALLOCATED( WaveAccC0Hyi ))      DEALLOCATE( WaveAccC0Hyi,     STAT=ErrStatTmp)
      IF (ALLOCATED( WaveAccC0V ))        DEALLOCATE( WaveAccC0V,       STAT=ErrStatTmp)
      IF (ALLOCATED( WaveDynP0B ))        DEALLOCATE( WaveDynP0B,       STAT=ErrStatTmp)
      IF (ALLOCATED( WaveDynPC0 ))        DEALLOCATE( WaveDynPC0,       STAT=ErrStatTmp)
      IF (ALLOCATED( WaveElev0 ))         DEALLOCATE( WaveElev0,        STAT=ErrStatTmp)
      IF (ALLOCATED( WaveElevC ))         DEALLOCATE( WaveElevC,        STAT=ErrStatTmp)
      IF (ALLOCATED( WaveVel0Hxi ))       DEALLOCATE( WaveVel0Hxi,      STAT=ErrStatTmp)
      IF (ALLOCATED( WaveVel0Hyi ))       DEALLOCATE( WaveVel0Hyi,      STAT=ErrStatTmp)
      IF (ALLOCATED( WaveVel0V ))         DEALLOCATE( WaveVel0V,        STAT=ErrStatTmp)
      IF (ALLOCATED( WaveVelC0Hxi ))      DEALLOCATE( WaveVelC0Hxi,     STAT=ErrStatTmp)
      IF (ALLOCATED( WaveVelC0Hyi ))      DEALLOCATE( WaveVelC0Hyi,     STAT=ErrStatTmp)
      IF (ALLOCATED( WaveVelC0V ))        DEALLOCATE( WaveVelC0V,       STAT=ErrStatTmp)
      IF (ALLOCATED( WvSpreadThetaIdx ))  DEALLOCATE( WvSpreadThetaIdx, STAT=ErrStatTmp)
      IF (ALLOCATED( tmpComplexArr ))     DEALLOCATE( tmpComplexArr,    STAT=ErrStatTmp)
      RETURN

   END SUBROUTINE CleanUp


END SUBROUTINE VariousWaves_Init

!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE GH_BladedWaves_Init ( InitInp, InitOut, ErrStat, ErrMsg )              
!  This routine initializes GH BLADED Waves
!----------------------------------------------------------------------------------------------------------------------------------
   TYPE(Waves_InitInputType),       INTENT(IN   )  :: InitInp     ! Input data for initialization routine
   TYPE(Waves_InitOutputType),      INTENT(INOUT)  :: InitOut     ! Initialization outputs      
   INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   ! GH Bladed wave data.
   REAL(ReKi), ALLOCATABLE      :: WaveElev0  (:)                                  ! Instantaneous elevation of incident waves at the WAMIT reference point (meters)
   
   INTEGER                       :: UnFA                          ! I/O unit number for the file needed for the GH Bladed wave data by FAST.
   INTEGER                       :: UnKi                          ! I/O unit number for the GH Bladed wave data file containing wave particle kinematics time history.
   INTEGER                       :: UnSu                          ! I/O unit number for the GH Bladed wave data file containing surface elevation time history.
   REAL(ReKi), ALLOCATABLE      :: GHWaveAcc (:,:)                                 ! Instantaneous acceleration of incident waves in the xi- (1), yi- (2), and zi- (3) directions, respectively, at each of the GHNWvDpth vertical locations in GH Bladed wave data files (m/s^2)
   REAL(ReKi), ALLOCATABLE      :: GHWaveDynP(:  )                                 ! Instantaneous dynamic pressure of incident waves                                                            at each of the GHNWvDpth vertical locations in GH Bladed wave data files (N/m^2)
   REAL(ReKi)                   :: GHWaveTime                                      ! Instantaneous simulation times in GH Bladed wave data files (sec)
   REAL(ReKi), ALLOCATABLE      :: GHWaveVel (:,:)                                 ! Instantaneous velocity     of incident waves in the xi- (1), yi- (2), and zi- (3) directions, respectively, at each of the GHNWvDpth vertical locations in GH Bladed wave data files (m/s  )
   REAL(ReKi), ALLOCATABLE      :: GHWvDpth  (:)                                   ! Vertical locations in GH Bladed wave data files.
   INTEGER                      :: GHNStepWave                                     ! Total number of time steps in the GH Bladed wave data files.
   INTEGER                      :: GHNWvDpth                                       ! Number of vertical locations in GH Bladed wave data files.
   INTEGER                      :: I                                               ! Generic index
   INTEGER                      :: I_Orig                                          ! The index of the time step from original (input) part of data
   INTEGER                      :: J                                               ! Generic index
   INTEGER                      :: K                                               ! Generic index
   INTEGER                      :: LastInd                                         ! Index into the arrays saved from the last call as a starting point for this call
   
   LOGICAL                      :: Reading                                        ! Flag to say whether or not we are still reading from the GH Bladed wave data files (files not exhausted).

      ! Temporary error handling variables
   INTEGER(IntKi)                :: ErrStatTmp                    ! Temporarary error status for procesing
   CHARACTER(1024)               :: ErrMsgTmp                     ! Temporary error message for processing
   

      ! Initialize ErrStat
         
   ErrStat = ErrID_None         
   ErrMsg  = "" 
   
   LastInd  = 1
   Reading  = .TRUE.

   ! Tell our nice users what is about to happen that may take a while:

   CALL WrScr1 ( ' Reading in wave data from GH Bladed files with root name "'//TRIM(InitInp%GHWvFile)//'".' )



   ! Perform some initialization computations including calculating the
   !   total number of time steps in the incident wave and ALLOCATing the
   !   arrays; initialize the unneeded values to zero:

   InitOut%NStepWave  = CEILING ( InitInp%WaveTMax/InitInp%WaveDT )                              ! Set NStepWave to an even integer
   IF ( MOD(InitOut%NStepWave,2) == 1 )  InitOut%NStepWave = InitOut%NStepWave + 1              !   larger or equal to WaveTMax/WaveDT.
   InitOut%NStepWave2 = InitOut%NStepWave/2

   ALLOCATE ( InitOut%WaveTime   (0:InitOut%NStepWave                    ) , STAT=ErrStatTmp )
   IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array InitOut%WaveTime.',  ErrStat,ErrMsg,'GH_BladedWaves_Init')

   ALLOCATE ( InitOut%WaveElevC0 (2,0:InitOut%NStepWave2                 ) , STAT=ErrStatTmp )
   IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array InitOut%WaveElevC0.',ErrStat,ErrMsg,'GH_BladedWaves_Init')

   ALLOCATE ( WaveElev0          (0:InitOut%NStepWave-1                  ) , STAT=ErrStatTmp )
   IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveElev0.',         ErrStat,ErrMsg,'GH_BladedWaves_Init')

   ALLOCATE ( InitOut%WaveElev   (0:InitOut%NStepWave,InitInp%NWaveElev  ) , STAT=ErrStatTmp )
   IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array InitOut%WaveElev.',  ErrStat,ErrMsg,'GH_BladedWaves_Init')

   ALLOCATE ( InitOut%WaveDynP0  (0:InitOut%NStepWave,InitInp%NWaveKin0  ) , STAT=ErrStatTmp )
   IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array InitOut%WaveDynP0.', ErrStat,ErrMsg,'GH_BladedWaves_Init')

   ALLOCATE ( InitOut%WaveVel0   (0:InitOut%NStepWave,InitInp%NWaveKin0,3) , STAT=ErrStatTmp )
   IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array InitOut%WaveVel0.',  ErrStat,ErrMsg,'GH_BladedWaves_Init')

   ALLOCATE ( InitOut%WaveAcc0   (0:InitOut%NStepWave,InitInp%NWaveKin0,3) , STAT=ErrStatTmp )
   IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array InitOut%WaveAcc0.',  ErrStat,ErrMsg,'GH_BladedWaves_Init')

      ! Now check if all the allocations worked properly
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL CleanUp()
      RETURN
   END IF

   InitOut%WaveDOmega = 0.0
   InitOut%WaveElevC0(1,0) = 0.0
   InitOut%WaveElevC0(2,0) = 0.0
   



   ! Open the file needed for the GH Bladed wave data by FAST, read in the
   !   input parameters, then close it again:

   CALL  GetNewUnit( UnFA )
   CALL  OpenFInpFile ( UnFA, TRIM(InitInp%GHWvFile)//'_FAST.txt', ErrStatTmp, ErrMsgTmp ) ! Open file.
   CALL  SetErrStat(ErrStatTmp,ErrMsgTmp,ErrStat,ErrMsg,'GH_BladedWaves_Init')
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL CleanUp()
      RETURN
   END IF



   ! GHNWvDpth - Number of vertical locations in GH Bladed wave data files.

   READ (UnFA,*)  GHNWvDpth

   IF ( GHNWvDpth <= 0 )  THEN
      CALL  SetErrStat(ErrID_Fatal,'GHNWvDpth must be greater than zero.',ErrStat,ErrMsg,'GH_BladedWaves_Init')
      CLOSE ( UnFA ) 
      RETURN
   END IF

   ! GHWvDpth - Vertical locations in GH Bladed wave data files.

   ALLOCATE ( GHWvDpth(GHNWvDpth) , STAT=ErrStatTmp )
   IF (ErrStatTmp /= 0) THEN
      CALL SetErrStat(ErrID_Fatal,'Cannot allocate array GHWvDpth.',ErrStat,ErrMsg,'GH_BladedWaves_Init')
      CALL CleanUp()
      CLOSE ( UnFA ) 
      RETURN
   END IF

   DO J = 1,GHNWvDpth   ! Loop through all vertical locations in the GH Bladed wave data files
      READ (UnFA,*)  GHWvDpth(J)
   END DO                ! J - All vertical locations in the GH Bladed wave data files

   IF ( GHWvDpth(1) /= -InitInp%WtrDpth )  THEN  !TODO: Verify this check is valid if MSL2SWL is /= 0  GJH 9/6/13
      CALL SetErrStat(ErrID_Fatal,'GHWvDpth(1) must be set to -WtrDpth when WaveMod is set to 5.',ErrStat,ErrMsg,'GH_BladedWaves_Init')
      CALL CleanUp()
      CLOSE ( UnFA ) 
      RETURN
   END IF

   CLOSE ( UnFA )                                        ! Close file.



   ! ALLOCATE arrays associated with the GH Bladed wave data:

   ALLOCATE ( GHWaveDynP(GHNWvDpth) , STAT=ErrStatTmp )
   IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array GHWaveDynP.',  ErrStat,ErrMsg,'GH_BladedWaves_Init')

   ALLOCATE ( GHWaveVel(GHNWvDpth,3) , STAT=ErrStatTmp )
   IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array GHWaveVel.',   ErrStat,ErrMsg,'GH_BladedWaves_Init')

   ALLOCATE ( GHWaveAcc(GHNWvDpth,3) , STAT=ErrStatTmp )
   IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array GHWaveAcc.',   ErrStat,ErrMsg,'GH_BladedWaves_Init')

      ! Now check if all the allocations worked properly
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL CleanUp()
      RETURN
   END IF




   ! Open the GH Bladed wave data files:

   CALL  GetNewUnit( UnKi )
   CALL  OpenFInpFile ( UnKi, TRIM(InitInp%GHWvFile)//'_kinematics.txt', ErrStatTmp, ErrMsgTmp )
   CALL  SetErrStat(ErrStatTmp,ErrMsgTmp,ErrStat,ErrMsg,'GH_BladedWaves_Init')
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL CleanUp()
      RETURN
   END IF

   CALL  GetNewUnit( UnSu )
   CALL  OpenFInpFile ( UnSu, TRIM(InitInp%GHWvFile)//'_surface.txt',ErrStatTmp, ErrMsgTmp )
   CALL  SetErrStat(ErrStatTmp,ErrMsgTmp,ErrStat,ErrMsg,'GH_BladedWaves_Init')
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL CleanUp()
      RETURN
   END IF


   ! Skip first line in the surface file:

   READ (UnSu,'()')


   ! Process data for all the time steps:

   DO I = 0,InitOut%NStepWave ! Loop through all time steps


   ! Calculate the array of simulation times at which the instantaneous
   !   elevation of, velocity of, acceleration of, and loads associated with
   !   the incident waves are to be determined:

      InitOut%WaveTime(I) = I*REAL(InitInp%WaveDT,ReKi)


      IF ( Reading )  THEN       ! .TRUE. if we are still reading from the GH Bladed wave data files.


   ! Let's read in data for this time step:

         READ (UnSu,*,IOSTAT=ErrStatTmp)  GHWaveTime, WaveElev0(I)

         IF ( ErrStatTmp == 0 )  THEN ! .TRUE. if there was no error reading in the line of data

            IF ( NINT( GHWaveTime/InitInp%WaveDT ) /= I )  THEN ! This is the same as: IF ( GHWaveTime /= WaveTime(I) ), but works better numerically
               ErrMsgTmp  = ' The input value of WaveDT is not consistent with the'// &
                              ' time step inherent in the GH Bladed wave data files.'
               CALL  SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,'GH_BladedWaves_Init')
               RETURN
            END IF

            DO J = 1,GHNWvDpth   ! Loop through all vertical locations in the GH Bladed wave data files
               READ (UnKi,*)  ( GHWaveVel(J,K), K=1,3 ), ( GHWaveAcc(J,K), K=1,3 ), GHWaveDynP(J)
            END DO                ! J - All vertical locations in the GH Bladed wave data files


   ! Let's interpolate GHWaveDynP, GHWaveVel, and GHWaveAcc to find
   !   WaveDynP0, WaveVel0, and WaveAcc0 if the elevation of the point
   !   defined by WaveKinzi0(J) lies within the range of GHWvDpth, else set
   !   WaveDynP0, WaveVel0, and WaveAcc0 to zero:

            DO J = 1,InitInp%NWaveKin0   ! Loop through all points where the incident wave kinematics will be computed
               IF ( ( InitInp%WaveKinzi0(J) < GHWvDpth(1) ) .OR. ( InitInp%WaveKinzi0(J) > GHWvDpth(GHNWvDpth) ) )  THEN ! .TRUE. if the elevation of the point defined by WaveKinzi0(J) lies outside the range of GHWvDpth
                  InitOut%WaveDynP0   (I,J  ) = 0.0
                  InitOut%WaveVel0    (I,J,:) = 0.0
                  InitOut%WaveAcc0    (I,J,:) = 0.0
               ELSE                                                                                      ! The elevation of the point defined by WaveKinzi0(J) must lie within the range of GHWvDpth; therefore, interpolate to find the incident wave kinematics at that elevation
                  InitOut%WaveDynP0   (I,J  ) = InterpStp ( InitInp%WaveKinzi0(J), GHWvDpth(:), GHWaveDynP(:),  &
                                                            LastInd,               GHNWvDpth                    )
                  DO K = 1,3     ! Loop through all xi- (1), yi- (2), and zi- (3) directions
                     InitOut%WaveVel0 (I,J,K) = InterpStp ( InitInp%WaveKinzi0(J), GHWvDpth(:), GHWaveVel(:,K), &
                                                            LastInd,               GHNWvDpth                    )
                     InitOut%WaveAcc0 (I,J,K) = InterpStp ( InitInp%WaveKinzi0(J), GHWvDpth(:), GHWaveAcc(:,K), &
                                                            LastInd,               GHNWvDpth                    )
                  END DO          ! K - All xi- (1), yi- (2), and zi- (3) directions
               END IF
            END DO                ! J - All points where the incident wave kinematics will be computed

         ELSE                    ! There must have been an error reading in the line of data

            GHNStepWave = I
            Reading     = .FALSE.

         END IF


      END IF


      IF ( .NOT. Reading )  THEN ! .TRUE. if we have finished reading from the GH Bladed wave data files.


   ! Let's reuse the input data to fill out the array:

         I_Orig = MOD( I, GHNStepWave )

         InitOut%WaveDynP0(I,:  ) = InitOut%WaveDynP0(I_Orig,:  )
         WaveElev0(I    ) = WaveElev0(I_Orig    )
         InitOut%WaveVel0 (I,:,:) = InitOut%WaveVel0 (I_Orig,:,:)
         InitOut%WaveAcc0 (I,:,:) = InitOut%WaveAcc0 (I_Orig,:,:)


      END IF


   END DO                ! I - All time steps


   ! Close the GH Bladed wave data files:

   CLOSE ( UnKi )
   CLOSE ( UnSu )



   ! Compute the incident wave elevations at each desired point on the still
   !   water level plane where it can be output; the only available point in
   !   the GH Bladed wave data is (xi=0.0,yi=0.0):

   DO J = 1,InitInp%NWaveElev   ! Loop through all points where the incident wave elevations can be output
      InitOut%WaveElev (:,J) = WaveElev0(:)
   END DO                ! J - All points where the incident wave elevations can be output


CONTAINS
   SUBROUTINE CleanUp( )

      IF (ALLOCATED( GHWaveAcc ))         DEALLOCATE( GHWaveAcc,        STAT=ErrStatTmp)
      IF (ALLOCATED( GHWaveDynP ))        DEALLOCATE( GHWaveDynP,       STAT=ErrStatTmp)
      IF (ALLOCATED( GHWaveVel ))         DEALLOCATE( GHWaveVel,        STAT=ErrStatTmp)
      IF (ALLOCATED( GHWvDpth ))          DEALLOCATE( GHWvDpth,         STAT=ErrStatTmp)
      IF (ALLOCATED( WaveElev0 ))         DEALLOCATE( WaveElev0,        STAT=ErrStatTmp)
      RETURN
   END SUBROUTINE CleanUp

END SUBROUTINE GH_BladedWaves_Init


!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Waves_Init( InitInp, u, p, x, xd, z, OtherState, y, Interval, InitOut, ErrStat, ErrMsg )
! This routine is called at the start of the simulation to perform initialization steps. 
! The parameters are set here and not changed during the simulation.
! The initial states and initial guess for the input are defined.
!..................................................................................................................................

      TYPE(Waves_InitInputType),       INTENT(IN   )  :: InitInp     ! Input data for initialization routine
      TYPE(Waves_InputType),           INTENT(  OUT)  :: u           ! An initial guess for the input; input mesh must be defined
      TYPE(Waves_ParameterType),       INTENT(  OUT)  :: p           ! Parameters      
      TYPE(Waves_ContinuousStateType), INTENT(  OUT)  :: x           ! Initial continuous states
      TYPE(Waves_DiscreteStateType),   INTENT(  OUT)  :: xd          ! Initial discrete states
      TYPE(Waves_ConstraintStateType), INTENT(  OUT)  :: z           ! Initial guess of the constraint states
      TYPE(Waves_OtherStateType),      INTENT(  OUT)  :: OtherState  ! Initial other/optimization states            
      TYPE(Waves_OutputType),          INTENT(  OUT)  :: y           ! Initial system outputs (outputs are not calculated; 
                                                                     !   only the output mesh is initialized)
      REAL(DbKi),                      INTENT(INOUT)  :: Interval    ! Coupling interval in seconds: the rate that 
                                                                     !   (1) Waves_UpdateStates() is called in loose coupling &
                                                                     !   (2) Waves_UpdateDiscState() is called in tight coupling.
                                                                     !   Input is the suggested time from the glue code; 
                                                                     !   Output is the actual coupling interval that will be used 
                                                                     !   by the glue code.
      TYPE(Waves_InitOutputType),      INTENT(  OUT)  :: InitOut     ! Output for initialization routine
      INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      
      
         ! Local Variables:
      INTEGER(IntKi)                                  :: ErrStatTmp  ! Temporary error status for processing
      CHARACTER(1024)                                 :: ErrMsgTmp   ! Temporary error message for procesing

     
!      TYPE(FFT_DataType)           :: FFT_Data                                        ! the instance of the FFT module we're using

      
   
         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrStatTmp  = ErrID_None
      ErrMsg  = ""               
      ErrMsgTmp   = ""
      
      
         ! Initialize the NWTC Subroutine Library
         
      CALL NWTC_Init(  )

  
      
         
         ! Define initialization-routine output here:
         
      !InitOut%WriteOutputHdr = (/ 'Time', 'Column2' /)
      !InitOut%WriteOutputUnt = (/ '(s)',  '(-)'     /)

      InitOut%RhoXg         = InitInp%WtrDens*InitInp%Gravity
      

         ! Set the minimum and maximum wave directions to WaveDir.  These are reset in the 
         ! subroutine calls as necessary.
      InitOut%WaveDirMin   = InitInp%WaveDir
      InitOut%WaveDirMax   = InitInp%WaveDir
      

            ! Initialize the variables associated with the incident wave:

      SELECT CASE ( InitInp%WaveMod ) ! Which incident wave kinematics model are we using?
      

      CASE ( 0 )              ! None=still water.

         CALL StillWaterWaves_Init( InitInp, InitOut, ErrStatTmp, ErrMsgTmp )
         CALL  SetErrStat(ErrStatTmp,ErrMsgTmp,ErrStat,ErrMsg,'Waves_Init')
         IF ( ErrStat >= AbortErrLev ) RETURN
            
      
      !CASE ( 3 )                 ! White-noise
      !   
      !   CALL WhiteNoiseWaves_Init( InitInp, InitOut, ErrStatTmp, ErrMsgTmp )
      !   CALL  SetErrStat(ErrStatTmp,ErrMsgTmp,ErrStat,ErrMsg,'Waves_Init')
      !   IF ( ErrStat >= AbortErrLev ) RETURN
         
      CASE ( 1, 2, 3, 4, 10 )       ! 1, 10: Plane progressive (regular) wave, 2: JONSWAP/Pierson-Moskowitz spectrum (irregular) wave, 3: white-noise, or 4: user-defined spectrum (irregular) wave.
      
         CALL VariousWaves_Init( InitInp, InitOut, ErrStatTmp, ErrMsgTmp )
         CALL  SetErrStat(ErrStatTmp,ErrMsgTmp,ErrStat,ErrMsg,'Waves_Init')
         IF ( ErrStat >= AbortErrLev ) RETURN

        
      CASE ( 5 )              ! GH Bladed wave data.   
         
         CALL GH_BladedWaves_Init( InitInp, InitOut, ErrStatTmp, ErrMsgTmp )
         CALL  SetErrStat(ErrStatTmp,ErrMsgTmp,ErrStat,ErrMsg,'Waves_Init')
         IF ( ErrStat >= AbortErrLev ) RETURN

      ENDSELECT
   
      
   u%DummyInput = 0.0
   p%DT = Interval
   p%WaveMultiDir = InitOut%WaveMultiDir     ! Flag to indicate multidirectional waves
   x%DummyContState = 0.0
   xd%DummyDiscState = 0.0
   z%DummyConstrState = 0.0
   OtherState%DummyOtherState = 0
   y%DummyOutput = 0.0
      
      
      
END SUBROUTINE Waves_Init


!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Waves_End( u, p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
! This routine is called at the end of the simulation.
!..................................................................................................................................

      TYPE(Waves_InputType),           INTENT(INOUT)  :: u           ! System inputs
      TYPE(Waves_ParameterType),       INTENT(INOUT)  :: p           ! Parameters     
      TYPE(Waves_ContinuousStateType), INTENT(INOUT)  :: x           ! Continuous states
      TYPE(Waves_DiscreteStateType),   INTENT(INOUT)  :: xd          ! Discrete states
      TYPE(Waves_ConstraintStateType), INTENT(INOUT)  :: z           ! Constraint states
      TYPE(Waves_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states            
      TYPE(Waves_OutputType),          INTENT(INOUT)  :: y           ! System outputs
      INTEGER(IntKi),                    INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                      INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

         ! Local error handling variables
      INTEGER(IntKi)                                  :: ErrStatTmp
      CHARACTER(1024)                                 :: ErrMsgTmp


         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrStatTmp  = ErrID_None
      ErrMsg  = ""               
      ErrMsgTmp   = ""
      
      
         ! Place any last minute operations or calculations here:


         ! Close files here:     
                  
                  

         ! Destroy the input data:
         
      CALL Waves_DestroyInput( u, ErrStatTmp, ErrMsgTmp )
      CALL  SetErrStat(ErrStatTmp,ErrMsgTmp,ErrStat,ErrMsg,'Waves_End')


         ! Destroy the parameter data:
         
      CALL Waves_DestroyParam( p, ErrStatTmp, ErrMsgTmp )
      CALL  SetErrStat(ErrStatTmp,ErrMsgTmp,ErrStat,ErrMsg,'Waves_End')


         ! Destroy the state data:
         
      CALL Waves_DestroyContState(   x,           ErrStatTmp, ErrMsgTmp )
      CALL  SetErrStat(ErrStatTmp,ErrMsgTmp,ErrStat,ErrMsg,'Waves_End')
      CALL Waves_DestroyDiscState(   xd,          ErrStatTmp, ErrMsgTmp )
      CALL  SetErrStat(ErrStatTmp,ErrMsgTmp,ErrStat,ErrMsg,'Waves_End')
      CALL Waves_DestroyConstrState( z,           ErrStatTmp, ErrMsgTmp )
      CALL  SetErrStat(ErrStatTmp,ErrMsgTmp,ErrStat,ErrMsg,'Waves_End')
      CALL Waves_DestroyOtherState(  OtherState,  ErrStatTmp, ErrMsgTmp )
      CALL  SetErrStat(ErrStatTmp,ErrMsgTmp,ErrStat,ErrMsg,'Waves_End')
         

         ! Destroy the output data:
         
      CALL Waves_DestroyOutput( y, ErrStatTmp, ErrMsgTmp )
      CALL  SetErrStat(ErrStatTmp,ErrMsgTmp,ErrStat,ErrMsg,'Waves_End')


      

END SUBROUTINE Waves_End
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Waves_UpdateStates( Time, u, p, x, xd, z, OtherState, ErrStat, ErrMsg )
! Loose coupling routine for solving for constraint states, integrating continuous states, and updating discrete states
! Constraint states are solved for input Time; Continuous and discrete states are updated for Time + Interval
!..................................................................................................................................
   
      REAL(DbKi),                       INTENT(IN   ) :: Time        ! Current simulation time in seconds
      TYPE(Waves_InputType),            INTENT(IN   ) :: u           ! Inputs at Time                    
      TYPE(Waves_ParameterType),        INTENT(IN   ) :: p           ! Parameters                              
      TYPE(Waves_ContinuousStateType),  INTENT(INOUT) :: x           ! Input: Continuous states at Time; 
                                                                       !   Output: Continuous states at Time + Interval
      TYPE(Waves_DiscreteStateType),    INTENT(INOUT) :: xd          ! Input: Discrete states at Time; 
                                                                      !   Output: Discrete states at Time  + Interval
      TYPE(Waves_ConstraintStateType),  INTENT(INOUT) :: z           ! Input: Initial guess of constraint states at Time;
                                                                      !   Output: Constraint states at Time
      TYPE(Waves_OtherStateType),       INTENT(INOUT) :: OtherState  ! Other/optimization states
      INTEGER(IntKi),                   INTENT(  OUT) :: ErrStat     ! Error status of the operation     
      CHARACTER(*),                     INTENT(  OUT) :: ErrMsg      ! Error message if ErrStat /= ErrID_None

         ! Local variables
         
      TYPE(Waves_ContinuousStateType)                 :: dxdt        ! Continuous state derivatives at Time
      TYPE(Waves_ConstraintStateType)                 :: z_Residual  ! Residual of the constraint state equations (Z)
         
      INTEGER(IntKi)                                  :: ErrStatTmp  ! Temporary Error status  for calls
      CHARACTER(1024)                                 :: ErrMsgTmp   ! Temporary Error message for calls
                        
         ! Initialize ErrStat
         
      ErrStat     = ErrID_None         
      ErrMsg      = ""               
      
      ErrStatTmp  = ErrID_None         
      ErrMsgTmp   = ""               
           
      
         ! Solve for the constraint states (z) here:
                           
         ! Check if the z guess is correct and update z with a new guess.
         ! Iterate until the value is within a given tolerance. 
                                    
      CALL Waves_CalcConstrStateResidual( Time, u, p, x, xd, z, OtherState, z_Residual, ErrStatTmp, ErrMsgTmp )
      CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, 'Waves_UpdateStates' )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Waves_DestroyConstrState( z_Residual, ErrStatTmp, ErrMsgTmp)
         CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, 'Waves_UpdateStates' )
         RETURN      
      END IF
         
      ! DO WHILE ( z_Residual% > tolerance )
      !
      !  z = 
      !
      !  CALL Waves_CalcConstrStateResidual( Time, u, p, x, xd, z, OtherState, z_Residual, ErrStatTmp, ErrMsgTmp )
      !  CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, 'Waves_UpdateStates' )
      !  IF ( ErrStat >= AbortErrLev ) THEN      
      !     CALL Waves_DestroyConstrState( z_Residual, ErrStatTmp, ErrMsgTmp)
      !     CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, 'Waves_UpdateStates' )
      !     RETURN      
      !  END IF
      !           
      ! END DO         
      
      
         ! Destroy z_Residual because it is not necessary for the rest of the subroutine:
            
      CALL Waves_DestroyConstrState( z_Residual, ErrStatTmp, ErrMsgTmp)
      CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, 'Waves_UpdateStates' )
      IF ( ErrStat >= AbortErrLev ) RETURN      
         
         
         
         ! Get first time derivatives of continuous states (dxdt):
      
      CALL Waves_CalcContStateDeriv( Time, u, p, x, xd, z, OtherState, dxdt, ErrStatTmp, ErrMsgTmp )
      CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, 'Waves_UpdateStates' )
      IF ( ErrStat >= AbortErrLev ) THEN      
         CALL Waves_DestroyContState( dxdt, ErrStatTmp, ErrMsgTmp)
         CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, 'Waves_UpdateStates')
         RETURN
      END IF
               
               
         ! Update discrete states:
         !   Note that xd [discrete state] is changed in Waves_UpdateDiscState(), so Waves_CalcOutput(),  
         !   Waves_CalcContStateDeriv(), and Waves_CalcConstrStates() must be called first (see above).
      
      CALL Waves_UpdateDiscState(Time, u, p, x, xd, z, OtherState, ErrStatTmp, ErrMsgTmp )   
      CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, 'Waves_UpdateStates' )
      IF ( ErrStat >= AbortErrLev ) THEN      
         CALL Waves_DestroyContState( dxdt, ErrStatTmp, ErrMsgTmp)
         CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, 'Waves_UpdateStates' )
         RETURN      
      END IF
         
         
         ! Integrate (update) continuous states (x) here:
         
      !x = function of dxdt and x


         ! Destroy dxdt because it is not necessary for the rest of the subroutine
            
      CALL Waves_DestroyContState( dxdt, ErrStatTmp, ErrMsgTmp)
      CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, 'Waves_UpdateStates' )
      IF ( ErrStat >= AbortErrLev ) RETURN      
     
   
      
END SUBROUTINE Waves_UpdateStates
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Waves_CalcOutput( Time, u, p, x, xd, z, OtherState, y, ErrStat, ErrMsg )   
! Routine for computing outputs, used in both loose and tight coupling.
!..................................................................................................................................
   
      REAL(DbKi),                      INTENT(IN   )  :: Time        ! Current simulation time in seconds
      TYPE(Waves_InputType),           INTENT(IN   )  :: u           ! Inputs at Time
      TYPE(Waves_ParameterType),       INTENT(IN   )  :: p           ! Parameters
      TYPE(Waves_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at Time
      TYPE(Waves_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at Time
      TYPE(Waves_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at Time
      TYPE(Waves_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
      TYPE(Waves_OutputType),          INTENT(INOUT)  :: y           ! Outputs computed at Time (Input only so that mesh con-
                                                                     !   nectivity information does not have to be recalculated)
      INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      
         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      
      
         ! Compute outputs here:
      y%DummyOutput    = 2.0_ReKi

     
               

END SUBROUTINE Waves_CalcOutput
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Waves_CalcContStateDeriv( Time, u, p, x, xd, z, OtherState, dxdt, ErrStat, ErrMsg )  
! Tight coupling routine for computing derivatives of continuous states
!..................................................................................................................................
   
      REAL(DbKi),                        INTENT(IN   )  :: Time        ! Current simulation time in seconds
      TYPE(Waves_InputType),           INTENT(IN   )  :: u           ! Inputs at Time                    
      TYPE(Waves_ParameterType),       INTENT(IN   )  :: p           ! Parameters                             
      TYPE(Waves_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at Time
      TYPE(Waves_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at Time
      TYPE(Waves_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at Time
      TYPE(Waves_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states                    
      TYPE(Waves_ContinuousStateType), INTENT(  OUT)  :: dxdt        ! Continuous state derivatives at Time
      INTEGER(IntKi),                    INTENT(  OUT)  :: ErrStat     ! Error status of the operation     
      CHARACTER(*),                      INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

               
         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      
      
         ! Compute the first time derivatives of the continuous states here:
      
      dxdt%DummyContState = 0.0_ReKi
         

END SUBROUTINE Waves_CalcContStateDeriv
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Waves_UpdateDiscState( Time, u, p, x, xd, z, OtherState, ErrStat, ErrMsg )   
! Tight coupling routine for updating discrete states
!..................................................................................................................................
   
      REAL(DbKi),                        INTENT(IN   )  :: Time        ! Current simulation time in seconds   
      TYPE(Waves_InputType),           INTENT(IN   )  :: u           ! Inputs at Time                       
      TYPE(Waves_ParameterType),       INTENT(IN   )  :: p           ! Parameters                                 
      TYPE(Waves_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at Time
      TYPE(Waves_DiscreteStateType),   INTENT(INOUT)  :: xd          ! Input: Discrete states at Time; 
                                                                       !   Output: Discrete states at Time + Interval
      TYPE(Waves_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at Time
      TYPE(Waves_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states           
      INTEGER(IntKi),                    INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                      INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

               
         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      
      
         ! Update discrete states here:
      
      ! StateData%DiscState = 

END SUBROUTINE Waves_UpdateDiscState
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Waves_CalcConstrStateResidual( Time, u, p, x, xd, z, OtherState, z_residual, ErrStat, ErrMsg )   
! Tight coupling routine for solving for the residual of the constraint state equations
!..................................................................................................................................
   
      REAL(DbKi),                        INTENT(IN   )  :: Time        ! Current simulation time in seconds   
      TYPE(Waves_InputType),           INTENT(IN   )  :: u           ! Inputs at Time                       
      TYPE(Waves_ParameterType),       INTENT(IN   )  :: p           ! Parameters                           
      TYPE(Waves_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at Time
      TYPE(Waves_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at Time
      TYPE(Waves_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at Time (possibly a guess)
      TYPE(Waves_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states                    
      TYPE(Waves_ConstraintStateType), INTENT(  OUT)  :: z_residual  ! Residual of the constraint state equations using  
                                                                       !     the input values described above      
      INTEGER(IntKi),                    INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                      INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

               
         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      
      
         ! Solve for the constraint states here:
      
      z_residual%DummyConstrState = 0.0_ReKi

END SUBROUTINE Waves_CalcConstrStateResidual
!!----------------------------------------------------------------------------------------------------------------------------------
!SUBROUTINE Waves_JacobianPInput( Time, u, p, x, xd, z, OtherState, dYdu, dXdu, dXddu, dZdu, ErrStat, ErrMsg )   
!! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) equations 
!! with respect to the inputs (u). The partial derivatives dY/du, dX/du, dXd/du, and DZ/du are returned.
!!..................................................................................................................................
!   
!      REAL(DbKi),                                INTENT(IN   )           :: Time       ! Current simulation time in seconds   
!      TYPE(Waves_InputType),                   INTENT(IN   )           :: u          ! Inputs at Time                       
!      TYPE(Waves_ParameterType),               INTENT(IN   )           :: p          ! Parameters                           
!      TYPE(Waves_ContinuousStateType),         INTENT(IN   )           :: x          ! Continuous states at Time
!      TYPE(Waves_DiscreteStateType),           INTENT(IN   )           :: xd         ! Discrete states at Time
!      TYPE(Waves_ConstraintStateType),         INTENT(IN   )           :: z          ! Constraint states at Time
!      TYPE(Waves_OtherStateType),              INTENT(INOUT)           :: OtherState ! Other/optimization states                    
!      !TYPE(Waves_PartialOutputPInputType),     INTENT(  OUT), OPTIONAL :: dYdu       ! Partial derivatives of output equations
!      !                                                                                 !   (Y) with respect to the inputs (u)
!      !TYPE(Waves_PartialContStatePInputType),  INTENT(  OUT), OPTIONAL :: dXdu       ! Partial derivatives of continuous state
!      !                                                                                 !   equations (X) with respect to inputs (u)
!      !TYPE(Waves_PartialDiscStatePInputType),  INTENT(  OUT), OPTIONAL :: dXddu      ! Partial derivatives of discrete state 
!      !                                                                                 !   equations (Xd) with respect to inputs (u)
!      !TYPE(Waves_PartialConstrStatePInputType),INTENT(  OUT), OPTIONAL :: dZdu       ! Partial derivatives of constraint state 
!                                                                                       !   equations (Z) with respect to inputs (u)
!      INTEGER(IntKi),                            INTENT(  OUT)           :: ErrStat    ! Error status of the operation
!      CHARACTER(*),                              INTENT(  OUT)           :: ErrMsg     ! Error message if ErrStat /= ErrID_None
!
!               
!         ! Initialize ErrStat
!         
!      ErrStat = ErrID_None         
!      ErrMsg  = ""               
!      
!      
!      IF ( PRESENT( dYdu ) ) THEN
!      
!         ! Calculate the partial derivative of the output equations (Y) with respect to the inputs (u) here:
!
!         dYdu%DummyOutput%DummyInput = 0
!
!      END IF
!      
!      IF ( PRESENT( dXdu ) ) THEN
!      
!         ! Calculate the partial derivative of the continuous state equations (X) with respect to the inputs (u) here:
!      
!         dXdu%DummyContState%DummyInput = 0
!
!      END IF
!      
!      IF ( PRESENT( dXddu ) ) THEN
!
!         ! Calculate the partial derivative of the discrete state equations (Xd) with respect to the inputs (u) here:
!
!         dXddu%DummyDiscState%DummyInput = 0
!
!      END IF
!      
!      IF ( PRESENT( dZdu ) ) THEN
!
!         ! Calculate the partial derivative of the constraint state equations (Z) with respect to the inputs (u) here:
!      
!         dZdu%DummyConstrState%DummyInput = 0
!
!      END IF
!
!
!END SUBROUTINE Waves_JacobianPInput
!!----------------------------------------------------------------------------------------------------------------------------------
!SUBROUTINE Waves_JacobianPContState( Time, u, p, x, xd, z, OtherState, dYdx, dXdx, dXddx, dZdx, ErrStat, ErrMsg )   
!! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) equations
!! with respect to the continuous states (x). The partial derivatives dY/dx, dX/dx, dXd/dx, and DZ/dx are returned.
!!..................................................................................................................................
!   
!      REAL(DbKi),                                    INTENT(IN   )           :: Time       ! Current simulation time in seconds   
!      TYPE(Waves_InputType),                       INTENT(IN   )           :: u          ! Inputs at Time                       
!      TYPE(Waves_ParameterType),                   INTENT(IN   )           :: p          ! Parameters                           
!      TYPE(Waves_ContinuousStateType),             INTENT(IN   )           :: x          ! Continuous states at Time
!      TYPE(Waves_DiscreteStateType),               INTENT(IN   )           :: xd         ! Discrete states at Time
!      TYPE(Waves_ConstraintStateType),             INTENT(IN   )           :: z          ! Constraint states at Time
!      TYPE(Waves_OtherStateType),                  INTENT(INOUT)           :: OtherState ! Other/optimization states                    
!      !TYPE(Waves_PartialOutputPContStateType),     INTENT(  OUT), OPTIONAL :: dYdx       ! Partial derivatives of output equations
!      !                                                                                     !   (Y) with respect to the continuous 
!      !                                                                                     !   states (x)
!      !TYPE(Waves_PartialContStatePContStateType),  INTENT(  OUT), OPTIONAL :: dXdx       ! Partial derivatives of continuous state
!      !                                                                                     !   equations (X) with respect to 
!      !                                                                                     !   the continuous states (x)
!      !TYPE(Waves_PartialDiscStatePContStateType),  INTENT(  OUT), OPTIONAL :: dXddx      ! Partial derivatives of discrete state 
!      !                                                                                     !   equations (Xd) with respect to 
!      !                                                                                     !   the continuous states (x)
!      !TYPE(Waves_PartialConstrStatePContStateType),INTENT(  OUT), OPTIONAL :: dZdx       ! Partial derivatives of constraint state
!                                                                                           !   equations (Z) with respect to 
!                                                                                           !   the continuous states (x)
!      INTEGER(IntKi),                                INTENT(  OUT)           :: ErrStat    ! Error status of the operation
!      CHARACTER(*),                                  INTENT(  OUT)           :: ErrMsg     ! Error message if ErrStat /= ErrID_None
!
!               
!         ! Initialize ErrStat
!         
!      ErrStat = ErrID_None         
!      ErrMsg  = ""               
!      
!      
!     
!      IF ( PRESENT( dYdx ) ) THEN
!
!         ! Calculate the partial derivative of the output equations (Y) with respect to the continuous states (x) here:
!
!         dYdx%DummyOutput%DummyContState = 0
!
!      END IF
!      
!      IF ( PRESENT( dXdx ) ) THEN
!      
!         ! Calculate the partial derivative of the continuous state equations (X) with respect to the continuous states (x) here:
!      
!         dXdx%DummyContState%DummyContState = 0
!
!      END IF
!      
!      IF ( PRESENT( dXddx ) ) THEN
!
!         ! Calculate the partial derivative of the discrete state equations (Xd) with respect to the continuous states (x) here:
!
!         dXddx%DummyDiscState%DummyContState = 0
!         
!      END IF
!      
!      IF ( PRESENT( dZdx ) ) THEN
!
!
!         ! Calculate the partial derivative of the constraint state equations (Z) with respect to the continuous states (x) here:
!      
!         dZdx%DummyConstrState%DummyContState = 0
!
!      END IF
!      
!
!   END SUBROUTINE Waves_JacobianPContState
!!----------------------------------------------------------------------------------------------------------------------------------
!SUBROUTINE Waves_JacobianPDiscState( Time, u, p, x, xd, z, OtherState, dYdxd, dXdxd, dXddxd, dZdxd, ErrStat, ErrMsg )   
!! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) equations
!! with respect to the discrete states (xd). The partial derivatives dY/dxd, dX/dxd, dXd/dxd, and DZ/dxd are returned.
!!..................................................................................................................................
!
!      REAL(DbKi),                                    INTENT(IN   )           :: Time       ! Current simulation time in seconds   
!      TYPE(Waves_InputType),                       INTENT(IN   )           :: u          ! Inputs at Time                       
!      TYPE(Waves_ParameterType),                   INTENT(IN   )           :: p          ! Parameters                           
!      TYPE(Waves_ContinuousStateType),             INTENT(IN   )           :: x          ! Continuous states at Time
!      TYPE(Waves_DiscreteStateType),               INTENT(IN   )           :: xd         ! Discrete states at Time
!      TYPE(Waves_ConstraintStateType),             INTENT(IN   )           :: z          ! Constraint states at Time
!      TYPE(Waves_OtherStateType),                  INTENT(INOUT)           :: OtherState ! Other/optimization states                    
!      !TYPE(Waves_PartialOutputPDiscStateType),     INTENT(  OUT), OPTIONAL :: dYdxd      ! Partial derivatives of output equations
!      !                                                                                     !  (Y) with respect to the discrete 
!      !                                                                                     !  states (xd)
!      !TYPE(Waves_PartialContStatePDiscStateType),  INTENT(  OUT), OPTIONAL :: dXdxd      ! Partial derivatives of continuous state
!      !                                                                                     !   equations (X) with respect to the 
!      !                                                                                     !   discrete states (xd)
!      !TYPE(Waves_PartialDiscStatePDiscStateType),  INTENT(  OUT), OPTIONAL :: dXddxd     ! Partial derivatives of discrete state 
!      !                                                                                     !   equations (Xd) with respect to the
!      !                                                                                     !   discrete states (xd)
!      !TYPE(Waves_PartialConstrStatePDiscStateType),INTENT(  OUT), OPTIONAL :: dZdxd      ! Partial derivatives of constraint state
!                                                                                           !   equations (Z) with respect to the 
!                                                                                           !   discrete states (xd)
!      INTEGER(IntKi),                                INTENT(  OUT)           :: ErrStat    ! Error status of the operation
!      CHARACTER(*),                                  INTENT(  OUT)           :: ErrMsg     ! Error message if ErrStat /= ErrID_None
!
!               
!         ! Initialize ErrStat
!         
!      ErrStat = ErrID_None         
!      ErrMsg  = ""               
!      
!      
!      IF ( PRESENT( dYdxd ) ) THEN
!      
!         ! Calculate the partial derivative of the output equations (Y) with respect to the discrete states (xd) here:
!
!         dYdxd%DummyOutput%DummyDiscState = 0
!
!      END IF
!      
!      IF ( PRESENT( dXdxd ) ) THEN
!
!         ! Calculate the partial derivative of the continuous state equations (X) with respect to the discrete states (xd) here:
!      
!         dXdxd%DummyContState%DummyDiscState = 0
!
!      END IF
!      
!      IF ( PRESENT( dXddxd ) ) THEN
!
!         ! Calculate the partial derivative of the discrete state equations (Xd) with respect to the discrete states (xd) here:
!
!         dXddxd%DummyDiscState%DummyDiscState = 0
!
!      END IF
!      
!      IF ( PRESENT( dZdxd ) ) THEN
!
!         ! Calculate the partial derivative of the constraint state equations (Z) with respect to the discrete states (xd) here:
!      
!         dZdxd%DummyConstrState%DummyDiscState = 0
!
!      END IF
!      
!
!
!END SUBROUTINE Waves_JacobianPDiscState
!!----------------------------------------------------------------------------------------------------------------------------------    
!SUBROUTINE Waves_JacobianPConstrState( Time, u, p, x, xd, z, OtherState, dYdz, dXdz, dXddz, dZdz, ErrStat, ErrMsg )   
!! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) equations
!! with respect to the constraint states (z). The partial derivatives dY/dz, dX/dz, dXd/dz, and DZ/dz are returned.
!!..................................................................................................................................
!   
!      REAL(DbKi),                                      INTENT(IN   )           :: Time       ! Current simulation time in seconds   
!      TYPE(Waves_InputType),                         INTENT(IN   )           :: u          ! Inputs at Time                       
!      TYPE(Waves_ParameterType),                     INTENT(IN   )           :: p          ! Parameters                           
!      TYPE(Waves_ContinuousStateType),               INTENT(IN   )           :: x          ! Continuous states at Time
!      TYPE(Waves_DiscreteStateType),                 INTENT(IN   )           :: xd         ! Discrete states at Time
!      TYPE(Waves_ConstraintStateType),               INTENT(IN   )           :: z          ! Constraint states at Time
!      TYPE(Waves_OtherStateType),                    INTENT(INOUT)           :: OtherState ! Other/optimization states                    
!      !TYPE(Waves_PartialOutputPConstrStateType),     INTENT(  OUT), OPTIONAL :: dYdz       ! Partial derivatives of output 
!                                                                                             !  equations (Y) with respect to the 
!                                                                                             !  constraint states (z)
!      !TYPE(Waves_PartialContStatePConstrStateType),  INTENT(  OUT), OPTIONAL :: dXdz       ! Partial derivatives of continuous
!      !                                                                                       !  state equations (X) with respect to 
!      !                                                                                       !  the constraint states (z)
!      !TYPE(Waves_PartialDiscStatePConstrStateType),  INTENT(  OUT), OPTIONAL :: dXddz      ! Partial derivatives of discrete state
!      !                                                                                       !  equations (Xd) with respect to the 
!      !                                                                                       !  constraint states (z)
!      !TYPE(Waves_PartialConstrStatePConstrStateType),INTENT(  OUT), OPTIONAL :: dZdz       ! Partial derivatives of constraint 
!                                                                                             ! state equations (Z) with respect to 
!                                                                                             !  the constraint states (z)
!      INTEGER(IntKi),                                  INTENT(  OUT)           :: ErrStat    ! Error status of the operation
!      CHARACTER(*),                                    INTENT(  OUT)           :: ErrMsg     ! Error message if ErrStat /= ErrID_None
!
!               
!         ! Initialize ErrStat
!         
!      ErrStat = ErrID_None         
!      ErrMsg  = ""               
!      
!      IF ( PRESENT( dYdz ) ) THEN
!      
!            ! Calculate the partial derivative of the output equations (Y) with respect to the constraint states (z) here:
!        
!         dYdz%DummyOutput%DummyConstrState = 0
!         
!      END IF
!      
!      IF ( PRESENT( dXdz ) ) THEN
!      
!            ! Calculate the partial derivative of the continuous state equations (X) with respect to the constraint states (z) here:
!         
!         dXdz%DummyContState%DummyConstrState = 0
!
!      END IF
!      
!      IF ( PRESENT( dXddz ) ) THEN
!
!            ! Calculate the partial derivative of the discrete state equations (Xd) with respect to the constraint states (z) here:
!
!         dXddz%DummyDiscState%DummyConstrState = 0
!
!      END IF
!      
!      IF ( PRESENT( dZdz ) ) THEN
!
!            ! Calculate the partial derivative of the constraint state equations (Z) with respect to the constraint states (z) here:
!         
!         dZdz%DummyConstrState%DummyConstrState = 0
!
!      END IF
!      
!
!END SUBROUTINE Waves_JacobianPConstrState

!----------------------------------------------------------------------------------------------------------------------------------

!=======================================================================
FUNCTION WheelerStretching ( zOrzPrime, Zeta, h, ForwardOrBackward, ErrStat, ErrMsg )


   ! This FUNCTION applies the principle of Wheeler stretching to
   ! (1-Forward) find the elevation where the wave kinematics are to
   ! be applied using Wheeler stretching or (2-Backword) find the
   ! elevation where the wave kinematics are computed before applying
   ! Wheeler stretching.  Wheeler stretching says that wave
   ! kinematics calculated using Airy theory at the mean sea level
   ! should actually be applied at the instantaneous free surface and
   ! that Airy wave kinematics computed at locations between the
   ! seabed and the mean sea level should be shifted vertically to
   ! new locations in proportion to their elevation above the seabed
   ! as follows:
   !
   ! Forward:  z(zPrime,Zeta,h) = ( 1 + Zeta/h )*zPrime + Zeta
   !
   ! or equivalently:
   !
   ! Backword: zPrime(z,Zeta,h) = ( z - Zeta )/( 1 + Zeta/h )
   !
   ! where,
   !   Zeta   = instantaneous elevation of incident waves
   !   h      = water depth
   !   z      = elevations where the wave kinematics are to be
   !            applied using Wheeler stretching
   !   zPrime = elevations where the wave kinematics are computed
   !            before applying Wheeler stretching



   IMPLICIT                        NONE


      ! Passed Variables:

   REAL(ReKi),     INTENT(IN )    :: h                                               ! Water depth (meters)
   REAL(ReKi)                     :: WheelerStretching                               ! This function = zPrime [forward] or z [backward] (meters)
   REAL(ReKi),     INTENT(IN )    :: Zeta                                            ! Instantaneous elevation of incident waves (meters)
   REAL(ReKi),     INTENT(IN )    :: zOrzPrime                                       ! Elevations where the wave kinematics are to be applied using Wheeler stretching, z, [forward] or elevations where the wave kinematics are computed before applying Wheeler stretching, zPrime, [backward] (meters)
   CHARACTER(1),   INTENT(IN )    :: ForwardOrBackWard                               ! A string holding the direction ('F'=Forward, 'B'=Backward) for applying Wheeler stretching.
   INTEGER(IntKi), INTENT(OUT)    :: ErrStat                                         ! Error status of the operation
   CHARACTER(*),   INTENT(OUT)    :: ErrMsg                                        ! Error message if ErrStat /= ErrID_None

   
      ! Initialize ErrStat
         
   ErrStat = ErrID_None         
   ErrMsg  = "" 

      ! Apply Wheeler stretching, depending on the direction:

   SELECT CASE ( ForwardOrBackWard )

   CASE ( 'F'  )  ! Forward

      WheelerStretching = ( 1.0 + Zeta/h )*zOrzPrime + Zeta


   CASE ( 'B' )   ! Backward

      WheelerStretching = ( zOrzPrime - Zeta )/( 1.0 + Zeta/h )


   CASE DEFAULT

      WheelerStretching = 0.0_ReKi
      
      ErrMsg = 'The last argument in routine WheelerStretching() must be ''F'' or ''B''.'
      ErrStat = ErrID_Fatal
      RETURN


   END SELECT



   RETURN
END FUNCTION WheelerStretching
   
END MODULE Waves
!**********************************************************************************************************************************
