!**********************************************************************************************************************************
! The Waves and Waves_Types modules make up a template for creating user-defined calculations in the FAST Modularization
! Framework. Waves_Types will be auto-generated based on a description of the variables for the module.
!..................................................................................................................................
! LICENSING
! Copyright (C) 2013-2015  National Renewable Energy Laboratory
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
MODULE Waves

   USE Waves_Types
   USE UserWaves
   USE NWTC_Library
   USE NWTC_FFTPACK
   USE NWTC_RandomNumber

   IMPLICIT NONE

   PRIVATE

   TYPE(ProgDesc), PARAMETER            :: Waves_ProgDesc = ProgDesc( 'Waves', '', '' )
   
   COMPLEX(SiKi),  PARAMETER, PUBLIC    :: ImagNmbr = (0.0_SiKi,1.0_SiKi)  ! The imaginary number, SQRT(-1.0)


      ! ..... Public Subroutines ...................................................................................................
   PUBLIC :: WavePkShpDefault                     ! Return the default value of the peak shape parameter of the incident wave spectrum
   PUBLIC :: Waves_Init                           ! Initialization routine


   PRIVATE:: WheelerStretching                    ! This FUNCTION applies the principle of Wheeler stretching to (1-Forward) find the elevation where the wave kinematics are to be applied using Wheeler stretching or (2-Backword)
   PRIVATE:: BoxMuller
   PRIVATE:: JONSWAP
   PUBLIC :: WaveNumber
   PRIVATE:: UserWaveSpctrm
   PRIVATE:: StillWaterWaves_Init
   PRIVATE:: VariousWaves_Init
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

   REAL(SiKi), INTENT(IN )      :: Hs                                              ! Significant wave height (meters)
   REAL(SiKi), INTENT(IN )      :: Tp                                              ! Peak spectral period (sec)
   REAL(SiKi)                   :: WavePkShpDefault                                ! This function = default value of the peak shape parameter of the incident wave spectrum conditioned on significant wave height and peak spectral period (-)


      ! Local Variables:

   REAL(SiKi)                   :: TpOvrSqrtHs                                     ! = Tp/SQRT(Hs) (s/SQRT(m))



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
   FUNCTION BoxMuller ( RNGType, NDAmp, Phase )

         ! This FUNCTION uses the Box-Muller method to turn two uniformly
         ! distributed randoms into two unit normal randoms, which are
         ! returned as real and imaginary components.

      IMPLICIT NONE

      COMPLEX(SiKi)                     :: BoxMuller                                  ! This function

         ! Passed Variables:

      INTEGER,    INTENT(IN)           :: RNGType
      LOGICAL,    INTENT(IN)           :: NDAmp                                       ! Flag for normally-distributed amplitudes
      REAL(SiKi), INTENT(IN), OPTIONAL :: Phase                                       ! Optional phase to override random phase (radians)

         ! Local Variables:

      REAL(SiKi)                   :: C1                                              ! Intermediate variable
      REAL(SiKi)                   :: C2                                              ! Intermediate variable
      REAL(SiKi)                   :: U1(1)                                           ! First  uniformly distributed random
      REAL(SiKi)                   :: U2(1)                                           ! Second uniformly distributed random

         ! Compute the two uniformly distributed randoms:
         ! NOTE: The first random, U1, cannot be zero else the LOG() function
         !       below will blow up; there is no restriction on the value of the
         !       second random, U2.

      U1 = 0.0
      DO WHILE ( U1(1) == 0.0 )
         CALL UniformRandomNumbers(RNGType, U1)
      END DO
      CALL UniformRandomNumbers(RNGType, U2)

         ! Compute intermediate variables:

      IF ( NDAmp )  THEN            ! Normally-distributed amplitudes
         C1 = SQRT( -2.0*LOG(U1(1)) )
      ELSE                          ! Constant amplitudes (ignore U1); therefore, C1 = SQRT( 2.0 ) = MEAN( SQRT( -2.0*LOG(U1) ) for a uniform distribution of U1 between 0 and 1
         C1 = SQRT(  2.0         )
      END IF

      IF ( PRESENT( Phase ) )  THEN ! Specified phase to replace random phase (ignore U2)
         C2 = Phase
      ELSE                          ! Uniformly-distributed phase
         C2 = TwoPi*U2(1)
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

      REAL(SiKi), INTENT(IN )      :: Gamma                                           ! Peak shape parameter (-)
      REAL(SiKi), INTENT(IN )      :: Hs                                              ! Significant wave height (meters)
      REAL(SiKi)                   :: JONSWAP                                         ! This function = JONSWAP wave spectrum, S (m^2/(rad/s))
      REAL(SiKi), INTENT(IN )      :: Omega                                           ! Wave frequency (rad/s)
      REAL(SiKi), INTENT(IN )      :: Tp                                              ! Peak spectral period (sec)


         ! Local Variables:

      REAL(SiKi)                   :: Alpha                                           ! Exponent on Gamma used in the spectral formulation (-)
      REAL(SiKi)                   :: C                                               ! Normalising factor used in the spectral formulation (-)
      REAL(SiKi)                   :: f                                               ! Wave frequency (Hz)
      REAL(SiKi)                   :: fp                                              ! Peak spectral frequency (Hz)
      REAL(SiKi)                   :: fpOvrf4                                         ! (fp/f)^4
      REAL(SiKi)                   :: Sigma                                           ! Scaling factor used in the spectral formulation (-)

       REAL(SiKi)                  :: Inv2Pi   =  0.15915494

         ! Compute the JONSWAP wave spectrum, unless Omega is zero, in which case,
         !   return zero:

      IF ( EqualRealNos(Omega, 0.0_SiKi) )  THEN  ! When .TRUE., the formulation below is ill-conditioned; thus, the known value of zero is returned.


         JONSWAP  = 0.0


      ELSE                       ! Omega > 0.0; forumulate the JONSWAP spectrum.


         ! Compute the wave frequency and peak spectral frequency in Hz:

         f        = Inv2Pi*Omega
         fp       = 1/Tp
         fpOvrf4  = (fp/f)**4


         ! Compute the normalising factor:

         C        = 1.0 - ( 0.287*LOG(GAMMA) )


         ! Compute Alpha:

         IF ( f <= fp )  THEN
            Sigma = 0.07
         ELSE
            Sigma = 0.09
         END IF

!bjj:         Alpha    = EXP( ( -0.5*( ( (f/fp) - 1.0 )/Sigma )**2 ) )
         Alpha    = EXP( ( -0.5*( ( (f*Tp) - 1.0 )/Sigma )**2 ) ) !this works even if Tp is 0 (but using f/fp doesn't)


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

      REAL(SiKi), INTENT(IN )      :: Omega                                           ! Wave frequency, rad/s.
      REAL(SiKi), INTENT(IN )      :: WaveDir                                         ! Incident wave propagation heading direction, degrees
      REAL(SiKi), INTENT(OUT)      :: WaveS1Sdd                                       ! One-sided power spectral density of the wave spectrum per unit time for the current frequency component and heading direction, m^2/(rad/s).

      CHARACTER(*), INTENT(IN )    :: DirRoot                                         ! The name of the root file including the full path to the current working directory.  This may be useful if you want this routine to write a permanent record of what it does to be stored with the simulation results: the results should be stored in a file whose name (including path) is generated by appending any suitable extension to DirRoot.



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
      REAL(SiKi), INTENT(IN )      :: Omega                                           ! Wave frequency (rad/s)
      REAL(SiKi)                   :: WaveNumber                                      ! This function = wavenumber, k (1/m)


         ! Local Variables:

      REAL(SiKi)                   :: A                                               ! A temporary variable used in the solution.
      REAL(SiKi)                   :: B                                               ! A temporary variable used in the solution.
      REAL(SiKi)                   :: C                                               ! A temporary variable used in the solution.
      REAL(SiKi)                   :: C2                                              ! A temporary variable used in the solution.
      REAL(SiKi)                   :: CC                                              ! A temporary variable used in the solution.
      REAL(SiKi)                   :: E2                                              ! A temporary variable used in the solution.
      REAL(SiKi)                   :: X0                                              ! A temporary variable used in the solution.



         ! Compute the wavenumber, unless Omega is zero, in which case, return
         !   zero:

      IF ( Omega == 0.0 )  THEN  ! When .TRUE., the formulation below is ill-conditioned; thus, the known value of zero is returned.


         WaveNumber = 0.0


      ELSE                       ! Omega > 0.0; solve for the wavenumber as usual.


         C  = Omega*Omega*REAL(h,SiKi)/REAL(g,SiKi)
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

            WaveNumber = ( X0 - ( B*C2*( 1.0 + (A*B*C*X0) ) ) )/REAL(h,SiKi)

         ELSE

            WaveNumber = X0/REAL(h,SiKi)

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

      REAL(SiKi)                   :: COSHNumOvrCOSHDen                               ! This function = COSH( k*( z + h ) )/COSH( k*h ) (-)
      REAL(ReKi), INTENT(IN )      :: h                                               ! Water depth ( h      >  0 ) (meters)
      REAL(SiKi), INTENT(IN )      :: k                                               ! Wave number ( k      >= 0 ) (1/m)
      REAL(SiKi), INTENT(IN )      :: z                                               ! Elevation   (-h <= z <= 0 ) (meters)



         ! Compute the hyperbolic numerator over denominator:

      IF ( k*h  > 89.4_SiKi )  THEN   ! When .TRUE., the shallow water formulation will trigger a floating point overflow error; however, COSH( k*( z + h ) )/COSH( k*h ) = EXP( k*z ) + EXP( -k*( z + 2*h ) ) for large k*h.  This equals the deep water formulation, EXP( k*z ), except near z = -h, because h > 14.23*wavelength (since k = 2*Pi/wavelength) in this case.

         COSHNumOvrCOSHDen = EXP( k*z ) + EXP( -k*( z + 2.0_SiKi*REAL(h,SiKi) ) )

      ELSE                       ! 0 < k*h <= 89.4; use the shallow water formulation.

         COSHNumOvrCOSHDen =REAL( COSH( k*( z + REAL(h,SiKi) ) ),R8Ki)/COSH( k*REAL(h,SiKi) )

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

      REAL(SiKi)                   :: COSHNumOvrSINHDen                               ! This function = COSH( k*( z + h ) )/SINH( k*h ) (-)
      REAL(ReKi), INTENT(IN )      :: h                                               ! Water depth ( h      >  0 ) (meters)
      REAL(SiKi), INTENT(IN )      :: k                                               ! Wave number ( k      >= 0 ) (1/m)
      REAL(SiKi), INTENT(IN )      :: z                                               ! Elevation   (-h <= z <= 0 ) (meters)



         ! Compute the hyperbolic numerator over denominator:


      IF (   k  < EPSILON(0.0_SiKi)  )  THEN  ! When .TRUE., the shallow water formulation is ill-conditioned; thus, HUGE(k) is returned to approximate the known value of infinity.

         COSHNumOvrSINHDen = HUGE( k )

      ELSEIF ( k*REAL(h,SiKi)  > 89.4_SiKi )  THEN  ! When .TRUE., the shallow water formulation will trigger a floating point overflow error; however, COSH( k*( z + h ) )/SINH( k*h ) = EXP( k*z ) + EXP( -k*( z + 2*h ) ) for large k*h.  This equals the deep water formulation, EXP( k*z ), except near z = -h, because h > 14.23*wavelength (since k = 2*Pi/wavelength) in this case.

         COSHNumOvrSINHDen = EXP( k*z ) + EXP( -k*( z + 2*REAL(h,SiKi) ) )

      ELSE                          ! 0 < k*h <= 89.4; use the shallow water formulation.

         COSHNumOvrSINHDen = COSH( k*( z + REAL(h,SiKi) ) )/SINH( k*REAL(h,SiKi) )

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

      REAL(SiKi)                   :: COTH                                            ! This function = COSH( X )/SINH( X ) (-)
      REAL(SiKi), INTENT(IN )      :: X                                               ! The argument (-)



         ! Compute the hyperbolic cotangent:

      IF ( X == 0.0_SiKi )  THEN   ! When .TRUE., the formulation below is ill-conditioned; thus, HUGE(X) is returned to approximate the known value of infinity.

         COTH = HUGE( X )

      ELSE                    ! X /= 0.0; use the numerically-stable computation of COTH(X) by means of TANH(X).

         COTH = 1.0_SiKi/TANH( X ) ! = COSH( X )/SINH( X )

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

      REAL(SiKi)                   :: SINHNumOvrSINHDen                               ! This function = SINH( k*( z + h ) )/SINH( k*h ) (-)
      REAL(ReKi), INTENT(IN )      :: h                                               ! Water depth ( h      >  0 ) (meters)
      REAL(SiKi), INTENT(IN )      :: k                                               ! Wave number ( k      >= 0 ) (1/m)
      REAL(SiKi), INTENT(IN )      :: z                                               ! Elevation   (-h <= z <= 0 ) (meters)



         ! Compute the hyperbolic numerator over denominator:

      IF (     k   == 0.0_SiKi  )  THEN  ! When .TRUE., the shallow water formulation is ill-conditioned; thus, the known value of unity is returned.

         SINHNumOvrSINHDen = 1.0

      ELSEIF ( k*REAL(h,SiKi) >  89.4_SiKi )  THEN  ! When .TRUE., the shallow water formulation will trigger a floating point overflow error; however, SINH( k*( z + h ) )/SINH( k*h ) = EXP( k*z ) - EXP( -k*( z + 2*h ) ) for large k*h.  This equals the deep water formulation, EXP( k*z ), except near z = -h, because h > 14.23*wavelength (since k = 2*Pi/wavelength) in this case.

         SINHNumOvrSINHDen = EXP( k*z ) - EXP( -k*( z + 2.0_SiKi*h ) )

      ELSE                          ! 0 < k*h <= 89.4; use the shallow water formulation.

         SINHNumOvrSINHDen = SINH( k*( z + REAL(h,SiKi) ) )/SINH( k*REAL(h,SiKi) )

      END IF



      RETURN
      END FUNCTION SINHNumOvrSINHDen


!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE StillWaterWaves_Init ( InitInp, InitOut, ErrStat, ErrMsg )
! This routine initializes the waves data for WaveMod = 0 , or still water waves option
!----------------------------------------------------------------------------------------------------------------------------------


   TYPE(Waves_InitInputType),       INTENT(IN   )  :: InitInp     ! Input data for initialization routine
   TYPE(Waves_InitOutputType),      INTENT(INOUT)  :: InitOut     ! Initialization output data
   INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   ! Local Variables
   INTEGER                      :: I, J,k, count                          ! Generic index
   INTEGER(IntKi)               :: ErrStatTmp                    ! Temporary error status
   CHARACTER(ErrMsgLen)         :: ErrMsgTmp                     ! Temporary error message
   character(*), parameter      :: RoutineName = 'StillWaterWaves_Init'

      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ""


   ! Initialize everything to zero:

      !>>>>>> COMPUTE INITOUT SCALARS InitOut%NStepWave, InitOut%NStepWave2, InitOut%WaveTMax, and InitOut%WaveDOmega for WAVEMOD = 0
      InitOut%NStepWave  = 2                ! We must have at least two elements in order to interpolate later on
      InitOut%NStepWave2 = 1
      InitOut%WaveTMax   = InitInp%WaveTMax ! bjj added this... I don't think it was set anywhere for this wavemod.
      InitOut%WaveDOmega = 0.0
      
      ! >>> Allocate and initialize (set to 0) InitOut arrays
      call Initial_InitOut_Arrays(InitOut, InitInp, 1.0_DbKi, ErrStatTmp, ErrMsgTmp);    CALL SetErrStat(ErrStatTmp,ErrMsgTmp,  ErrStat,ErrMsg,RoutineName)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      IF ( ErrStat >= AbortErrLev ) RETURN


      ! Add the current velocities to the wave velocities:
      count = 0

      !DO J = 1,InitInp%NWaveKinGrid      ! Loop through all points where the incident wave kinematics will be computed
      do k = 1, InitInp%NGrid(3)
         do j = 1, InitInp%NGrid(2)
            do i = 1, InitInp%NGrid(1)
               count = count + 1
               InitOut%WaveVel(:,i,j,k,1) =  InitInp%CurrVxi(count)  ! xi-direction
               InitOut%WaveVel(:,i,j,k,2) =  InitInp%CurrVyi(count)  ! yi-direction
            end do
         end do
      end do

     ! END DO                ! J - All points where the incident wave kinematics will be computed

END SUBROUTINE StillWaterWaves_Init


!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE VariousWaves_Init ( InitInp, InitOut, ErrStat, ErrMsg )
! Compute the wave kinematics and related information for  Plane progressive (regular) wave, JONSWAP/Pierson-Moskowitz spectrum
! (irregular) wave, or user-defined spectrum (irregular) wave.
!----------------------------------------------------------------------------------------------------------------------------------

   TYPE(Waves_InitInputType),       INTENT(IN   )  :: InitInp     ! Input data for initialization routine
   TYPE(Waves_InitOutputType),      INTENT(INOUT)  :: InitOut     ! Output data
   INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None



   ! Local Variables
   COMPLEX(SiKi)                :: ImagOmega                                       ! = ImagNmbr*Omega (rad/s)
   COMPLEX(SiKi), ALLOCATABLE   :: PWaveAccC0HxiPz0 (:,:)                            ! Partial derivative of WaveAccC0Hxi(:) with respect to zi at zi = 0 (1/s^2)
   COMPLEX(SiKi), ALLOCATABLE   :: PWaveAccC0HyiPz0 (:,:)                            ! Partial derivative of WaveAccC0Hyi(:) with respect to zi at zi = 0 (1/s^2)
   COMPLEX(SiKi), ALLOCATABLE   :: PWaveAccC0VPz0 (:,:)                              ! Partial derivative of WaveAccC0V  (:) with respect to zi at zi = 0 (1/s^2)
   COMPLEX(SiKi), ALLOCATABLE   :: PWaveDynPC0BPz0(:,:)                              ! Partial derivative of WaveDynPC0B (:) with respect to zi at zi = 0 (N/m  )
   COMPLEX(SiKi), ALLOCATABLE   :: PWaveVelC0HxiPz0 (:,:)                            ! Partial derivative of WaveVelC0Hxi(:) with respect to zi at zi = 0 (1/s  )
   COMPLEX(SiKi), ALLOCATABLE   :: PWaveVelC0HyiPz0 (:,:)                            ! Partial derivative of WaveVelC0Hyi(:) with respect to zi at zi = 0 (1/s  )
   COMPLEX(SiKi), ALLOCATABLE   :: PWaveVelC0VPz0 (:,:)                              ! Partial derivative of WaveVelC0V  (:) with respect to zi at zi = 0 (1/s  )
   COMPLEX(SiKi), ALLOCATABLE   :: WaveAccC0Hxi(:,:)                               ! Discrete Fourier transform of the instantaneous horizontal acceleration in x-direction of incident waves before applying stretching at the zi-coordinates for points (m/s^2)
   COMPLEX(SiKi), ALLOCATABLE   :: WaveAccC0Hyi(:,:)                               ! Discrete Fourier transform of the instantaneous horizontal acceleration in y-direction of incident waves before applying stretching at the zi-coordinates for points (m/s^2)
   COMPLEX(SiKi), ALLOCATABLE   :: WaveAccC0V(:,:)                                 ! Discrete Fourier transform of the instantaneous vertical   acceleration                of incident waves before applying stretching at the zi-coordinates for points (m/s^2)
   COMPLEX(SiKi), ALLOCATABLE   :: WaveDynPC0(:,:)                                 ! Discrete Fourier transform of the instantaneous dynamic pressure                       of incident waves before applying stretching at the zi-coordinates for points (N/m^2)
   COMPLEX(SiKi), ALLOCATABLE   :: WaveVelC0Hxi(:,:)                               ! Discrete Fourier transform of the instantaneous horizontal velocity                    of incident waves before applying stretching at the zi-coordinates for points (m/s)
   COMPLEX(SiKi), ALLOCATABLE   :: WaveVelC0Hyi(:,:)                               ! Discrete Fourier transform of the instantaneous horizontal velocity in x-direction     of incident waves before applying stretching at the zi-coordinates for points (m/s)
   COMPLEX(SiKi), ALLOCATABLE   :: WaveVelC0V(:,:)                                 ! Discrete Fourier transform of the instantaneous vertical   velocity in y-direction     of incident waves before applying stretching at the zi-coordinates for points (m/s)

   REAL(SiKi), ALLOCATABLE      :: CosWaveDir(:)                                   ! COS( WaveDirArr(I) ) -- Each wave frequency has a unique wave direction.
   REAL(SiKi), ALLOCATABLE      :: GHWaveAcc (:,:)                                 ! Instantaneous acceleration of incident waves in the xi- (1), yi- (2), and zi- (3) directions, respectively, at each of the GHNWvDpth vertical locations in GH Bladed wave data files (m/s^2)
   REAL(SiKi), ALLOCATABLE      :: GHWaveDynP(:  )                                 ! Instantaneous dynamic pressure of incident waves                                                            at each of the GHNWvDpth vertical locations in GH Bladed wave data files (N/m^2)

   REAL(SiKi), ALLOCATABLE      :: GHWaveVel (:,:)                                 ! Instantaneous velocity     of incident waves in the xi- (1), yi- (2), and zi- (3) directions, respectively, at each of the GHNWvDpth vertical locations in GH Bladed wave data files (m/s  )
   REAL(SiKi), ALLOCATABLE      :: GHWvDpth  (:)                                   ! Vertical locations in GH Bladed wave data files.

   REAL(SiKi), ALLOCATABLE      :: PWaveAcc0HxiPz0(:,:)                              ! Partial derivative of WaveAcc0Hxi(:) with respect to zi at zi = 0 (1/s^2)
   REAL(SiKi), ALLOCATABLE      :: PWaveAcc0HyiPz0(:,:)                              ! Partial derivative of WaveAcc0Hyi(:) with respect to zi at zi = 0 (1/s^2)
   REAL(SiKi), ALLOCATABLE      :: PWaveAcc0VPz0  (:,:)                              ! Partial derivative of WaveAcc0V  (:) with respect to zi at zi = 0 (1/s^2)
   REAL(SiKi), ALLOCATABLE      :: PWaveDynP0BPz0 (:,:)                              ! Partial derivative of WaveDynP0B (:) with respect to zi at zi = 0 (N/m  )
   REAL(SiKi), ALLOCATABLE      :: PWaveVel0HxiPz0(:,:)                              ! Partial derivative of WaveVel0Hxi(:) with respect to zi at zi = 0 (1/s  )
   REAL(SiKi), ALLOCATABLE      :: PWaveVel0HyiPz0(:,:)                              ! Partial derivative of WaveVel0Hyi(:) with respect to zi at zi = 0 (1/s  )
   REAL(SiKi), ALLOCATABLE      :: PWaveVel0VPz0  (:,:)                              ! Partial derivative of WaveVel0V  (:) with respect to zi at zi = 0 (1/s  )

   REAL(SiKi), ALLOCATABLE      :: SinWaveDir     (:)                              ! SIN( WaveDirArr(I) )
   REAL(SiKi), ALLOCATABLE      :: WaveAcc0Hxi (:,:)                               ! Instantaneous horizontal acceleration in x-direction of incident waves before applying stretching at the zi-coordinates for points (m/s^2)
   REAL(SiKi), ALLOCATABLE      :: WaveAcc0Hyi (:,:)                               ! Instantaneous horizontal acceleration in y-direction of incident waves before applying stretching at the zi-coordinates for points (m/s^2)

   REAL(SiKi), ALLOCATABLE      :: WaveAcc0V (:,:)                                 ! Instantaneous vertical   acceleration of incident waves before applying stretching at the zi-coordinates for points (m/s^2)

   REAL(SiKi), ALLOCATABLE      :: WaveDynP0B(:,:)                                 ! Instantaneous dynamic pressure        of incident waves before applying stretching at the zi-coordinates for points (N/m^2)

   COMPLEX(SiKi)                :: WaveElevxiPrime0
   REAL(SiKi), ALLOCATABLE      :: WaveKinzi0Prime(:)                              ! zi-coordinates for points where the incident wave kinematics will be computed before applying stretching; these are relative to the mean see level (meters)
   INTEGER   , ALLOCATABLE      :: WaveKinPrimeMap(:)
   REAL(SiKi)                   :: WaveNmbr                                        ! Wavenumber of the current frequency component (1/meter)
   REAL(SiKi), ALLOCATABLE      :: WaveVel0Hxi    (:,:)                            ! Instantaneous xi-direction velocity   of incident waves before applying stretching at the zi-coordinates for points (m/s  )
   REAL(SiKi), ALLOCATABLE      :: WaveVel0Hyi    (:,:)                            ! Instantaneous yi-direction velocity   of incident waves before applying stretching at the zi-coordinates for points (m/s  )
   REAL(SiKi), ALLOCATABLE      :: WaveVel0V (:,:)                                 ! Instantaneous vertical     velocity   of incident waves before applying stretching at the zi-coordinates for points (m/s  )
   INTEGER                      :: I,count                                               ! Generic index
   INTEGER                      :: J                                               ! Generic index
   INTEGER                      :: K                                               ! Generic index
   INTEGER                      :: NWaveKin0Prime                                  ! Number of points where the incident wave kinematics will be computed before applying stretching to the instantaneous free surface (-)
   integer                      :: primeCount                                      ! Counter for locations before applying stretching
   COMPLEX(SiKi)                :: tmpComplex                                      ! A temporary varible to hold the complex value of the wave elevation before storing it into a REAL array
   COMPLEX(SiKi),ALLOCATABLE    :: tmpComplexArr(:)                                ! A temporary array (0:NStepWave2-1) for FFT use.
   TYPE(FFT_DataType)           :: FFT_Data                                        ! the instance of the FFT module we're using

   REAL(SiKi), ALLOCATABLE      :: WaveS1SddArr(:)                                 !< One-sided power spectral density of the wave spectrum at all non-negative frequencies (m^2/(rad/s))
   REAL(SiKi), ALLOCATABLE      :: OmegaArr(:)                                     !< Array of all non-negative angular frequencies (rad/s)
   
   ! Variables for MacCamy-Fuchs model
   REAL(SiKi)                   :: ka
   REAL(SiKi)                   :: JPrime
   REAL(SiKi)                   :: YPrime
   REAL(SiKi)                   :: HPrime
   REAL(SiKi)                   :: MCFC
   COMPLEX(SiKi), ALLOCATABLE   :: WaveAccC0HxiMCF(:,:)                               ! Discrete Fourier transform of the instantaneous horizontal acceleration in x-direction of incident waves before applying stretching at the zi-coordinates for points (m/s^2)
   COMPLEX(SiKi), ALLOCATABLE   :: WaveAccC0HyiMCF(:,:)                               ! Discrete Fourier transform of the instantaneous horizontal acceleration in y-direction of incident waves before applying stretching at the zi-coordinates for points (m/s^2)
   COMPLEX(SiKi), ALLOCATABLE   :: WaveAccC0VMCF(:,:)                                 ! Discrete Fourier transform of the instantaneous vertical   acceleration                of incident waves before applying stretching at the zi-coordinates for points (m/s^2)
   COMPLEX(SiKi), ALLOCATABLE   :: PWaveAccC0HxiMCFPz0(:,:)                               ! Discrete Fourier transform of the instantaneous horizontal acceleration in x-direction of incident waves before applying stretching at the zi-coordinates for points (m/s^2)
   COMPLEX(SiKi), ALLOCATABLE   :: PWaveAccC0HyiMCFPz0(:,:)                               ! Discrete Fourier transform of the instantaneous horizontal acceleration in y-direction of incident waves before applying stretching at the zi-coordinates for points (m/s^2)
   COMPLEX(SiKi), ALLOCATABLE   :: PWaveAccC0VMCFPz0(:,:)                                 ! Discrete Fourier transform of the instantaneous vertical   acceleration                of incident waves before applying stretching at the zi-coordinates for points (m/s^2)
   REAL(SiKi), ALLOCATABLE   :: WaveAcc0HxiMCF(:,:)                               ! Discrete Fourier transform of the instantaneous horizontal acceleration in x-direction of incident waves before applying stretching at the zi-coordinates for points (m/s^2)
   REAL(SiKi), ALLOCATABLE   :: WaveAcc0HyiMCF(:,:)                               ! Discrete Fourier transform of the instantaneous horizontal acceleration in y-direction of incident waves before applying stretching at the zi-coordinates for points (m/s^2)
   REAL(SiKi), ALLOCATABLE   :: WaveAcc0VMCF(:,:)                                 ! Discrete Fourier transform of the instantaneous vertical   acceleration                of incident waves before applying stretching at the zi-coordinates for points (m/s^2)
   REAL(SiKi), ALLOCATABLE   :: PWaveAcc0HxiMCFPz0(:,:)                               ! Discrete Fourier transform of the instantaneous horizontal acceleration in x-direction of incident waves before applying stretching at the zi-coordinates for points (m/s^2)
   REAL(SiKi), ALLOCATABLE   :: PWaveAcc0HyiMCFPz0(:,:)                               ! Discrete Fourier transform of the instantaneous horizontal acceleration in y-direction of incident waves before applying stretching at the zi-coordinates for points (m/s^2)
   REAL(SiKi), ALLOCATABLE   :: PWaveAcc0VMCFPz0(:,:)                                 ! Discrete Fourier transform of the instantaneous vertical   acceleration                of incident waves before applying stretching at the zi-coordinates for points (m/s^2)
      

   ! Variables for error handling
   INTEGER(IntKi)               :: ErrStatTmp                                      !< Temporary error status
   CHARACTER(ErrMsgLen)         :: ErrMsgTmp                                       !< Temporary error message
   character(*), parameter      :: RoutineName = 'VariousWaves_Init'

      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ""

      ! Tell our users what is about to happen that may take a while:
      CALL WrScr ( ' Generating incident wave kinematics and current time history.' )



      ! Determine the number of, NWaveKin0Prime, and the zi-coordinates for,
      !   WaveKinzi0Prime(:), points where the incident wave kinematics will be
      !   computed before applying stretching to the instantaneous free surface.
      !   The locations are relative to the mean see level.  

         NWaveKin0Prime = 0
         DO J = 1,InitInp%NWaveKinGrid   ! Loop through all mesh points  where the incident wave kinematics will be computed
               ! NOTE: We test to 0 instead of MSL2SWL because the locations of WaveKinGridzi and WtrDpth have already been adjusted using MSL2SWL
           IF (    InitInp%WaveKinGridzi(J) >= -InitInp%WtrDpth .AND. InitInp%WaveKinGridzi(J) <= 0 )  THEN
               NWaveKin0Prime = NWaveKin0Prime + 1
           END IF
         END DO                ! J - All Morison nodes where the incident wave kinematics will be computed



      ! ALLOCATE the WaveKinzi0Prime(:) array and compute its elements here:

         ALLOCATE ( WaveKinzi0Prime(NWaveKin0Prime) , STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveKinzi0Prime.',ErrStat,ErrMsg,RoutineName)

         ALLOCATE ( WaveKinPrimeMap(NWaveKin0Prime) , STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveKinPrimeMap.',ErrStat,ErrMsg,RoutineName)

         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF


         I = 1

         DO J = 1,InitInp%NWaveKinGrid ! Loop through all points where the incident wave kinematics will be computed without stretching
               ! NOTE: We test to 0 instead of MSL2SWL because the locations of WaveKinGridzi and WtrDpth have already been adjusted using MSL2SWL
            IF (    InitInp%WaveKinGridzi(J) >= -InitInp%WtrDpth .AND. InitInp%WaveKinGridzi(J) <= 0 )  THEN

               WaveKinzi0Prime(I) =  InitInp%WaveKinGridzi(J)
               WaveKinPrimeMap(I) =  J
               I = I + 1

            END IF

         END DO                   ! J - All points where the incident wave kinematics will be computed without stretching



      ! Perform some initialization computations including calculating the total number of frequency
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




      ! Set new value for NStepWave so that the FFT algorithms are efficient.  Note that if this method is changed, the method
      ! used to calculate the number of multidirectional wave directions (WaveNDir) and the UserWaveElevations_Init subroutine
      ! will need to be updated.

      !>>>>>> COMPUTE INITOUT SCALARS InitOut%NStepWave, InitOut%NStepWave2, InitOut%WaveTMax, and InitOut%WaveDOmega for WAVEMOD = 1,2,3,4,10 (5 and 7 also call this routine, but have been set already)
      ! NOTE:  For WaveMod = 5, NStepWave and several other things were already set in the UserWaveElevations_Init routine
      !        using file information (an FFT was performed there, so the information was needed before now).
      !        Same with WaveMod = 7. With WaveMod = 7, WaveDirArr is also populated in UserWaveComponents_Init routine. 
      !        Need to make sure the wave-direction in formation is not overwritten later. 
      IF (InitInp%WaveMod /= 5 .AND. InitInp%WaveMod /= 7) THEN
         InitOut%NStepWave    = CEILING ( InitInp%WaveTMax/InitInp%WaveDT )               ! Set NStepWave to an even integer ...
         IF ( MOD(InitOut%NStepWave,2) == 1 )  InitOut%NStepWave = InitOut%NStepWave + 1  !   ... larger or equal to WaveTMax/WaveDT.
         
         InitOut%NStepWave2   = MAX( InitOut%NStepWave/2, 1 )                             ! Make sure that NStepWave is an even product of small factors (PSF) that is
         InitOut%NStepWave    = 2 * PSF( InitOut%NStepWave2, 9 )                          !   greater or equal to WaveTMax/WaveDT to ensure that the FFT is efficient.

         InitOut%NStepWave2   = InitOut%NStepWave/2                                       ! Update the value of NStepWave2 based on the value needed for NStepWave.
         InitOut%WaveTMax     = InitOut%NStepWave*InitInp%WaveDT                          ! Update the value of WaveTMax   based on the value needed for NStepWave.
         InitOut%WaveDOmega   = TwoPi/InitOut%WaveTMax                                    ! Compute the frequency step for incident wave calculations.
      
         ! >>> Allocate and initialize (set to 0) InitOut arrays
         call Initial_InitOut_Arrays(InitOut, InitInp, InitInp%WaveDT, ErrStatTmp, ErrMsgTmp);    CALL SetErrStat(ErrStatTmp,ErrMsgTmp,  ErrStat,ErrMsg,RoutineName)
      ENDIF
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      

      ! Allocate all the arrays we need.
      ALLOCATE ( tmpComplexArr(0:InitOut%NStepWave2                        ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array tmpComplexArr.',     ErrStat,ErrMsg,RoutineName)

      ALLOCATE ( WaveDynPC0        (0:InitOut%NStepWave2 ,NWaveKin0Prime   ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveDynPC0.',        ErrStat,ErrMsg,RoutineName)

      ALLOCATE ( WaveVelC0Hxi      (0:InitOut%NStepWave2 ,NWaveKin0Prime   ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveVelC0Hxi.',      ErrStat,ErrMsg,RoutineName)

      ALLOCATE ( WaveVelC0Hyi      (0:InitOut%NStepWave2 ,NWaveKin0Prime   ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveVelC0Hyi.',      ErrStat,ErrMsg,RoutineName)

      ALLOCATE ( WaveVelC0V        (0:InitOut%NStepWave2 ,NWaveKin0Prime   ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveVelC0V.',        ErrStat,ErrMsg,RoutineName)

      ALLOCATE ( WaveAccC0Hxi      (0:InitOut%NStepWave2 ,NWaveKin0Prime   ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveAccC0Hxi.',      ErrStat,ErrMsg,RoutineName)

      ALLOCATE ( WaveAccC0Hyi      (0:InitOut%NStepWave2 ,NWaveKin0Prime   ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveAccC0Hyi.',      ErrStat,ErrMsg,RoutineName)

      ALLOCATE ( WaveAccC0V        (0:InitOut%NStepWave2 ,NWaveKin0Prime   ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveAccC0V.',        ErrStat,ErrMsg,RoutineName)


      ALLOCATE ( WaveDynP0B        (0:InitOut%NStepWave-1,NWaveKin0Prime   ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveDynP0B.',        ErrStat,ErrMsg,RoutineName)

      ALLOCATE ( WaveVel0Hxi       (0:InitOut%NStepWave-1,NWaveKin0Prime   ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveVel0Hxi.',       ErrStat,ErrMsg,RoutineName)

      ALLOCATE ( WaveVel0Hyi       (0:InitOut%NStepWave-1,NWaveKin0Prime   ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveVel0Hyi.',       ErrStat,ErrMsg,RoutineName)

      ALLOCATE ( WaveVel0V         (0:InitOut%NStepWave-1,NWaveKin0Prime   ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveVel0V.',         ErrStat,ErrMsg,RoutineName)

      ALLOCATE ( WaveAcc0Hxi       (0:InitOut%NStepWave-1,NWaveKin0Prime   ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveAcc0Hxi.',       ErrStat,ErrMsg,RoutineName)

      ALLOCATE ( WaveAcc0Hyi       (0:InitOut%NStepWave-1,NWaveKin0Prime   ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveAcc0Hyi.',       ErrStat,ErrMsg,RoutineName)

      ALLOCATE ( WaveAcc0V         (0:InitOut%NStepWave-1,NWaveKin0Prime   ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveAcc0V.',         ErrStat,ErrMsg,RoutineName)
      
      
      IF (InitInp%MCFD > 0.0_SiKi) THEN ! MacCamy-Fuchs model
       
         ALLOCATE ( WaveAccC0HxiMCF(0:InitOut%NStepWave2 ,NWaveKin0Prime   ), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveAccC0HxiMCF.',      ErrStat,ErrMsg,RoutineName)

         ALLOCATE ( WaveAccC0HyiMCF(0:InitOut%NStepWave2 ,NWaveKin0Prime   ), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveAccC0HyiMCF.',      ErrStat,ErrMsg,RoutineName)

         ALLOCATE ( WaveAccC0VMCF  (0:InitOut%NStepWave2 ,NWaveKin0Prime   ), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveAccC0VMCF.',        ErrStat,ErrMsg,RoutineName)
  
         ALLOCATE ( WaveAcc0HxiMCF (0:InitOut%NStepWave-1,NWaveKin0Prime   ), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveAcc0HxiMCF.',       ErrStat,ErrMsg,RoutineName)

         ALLOCATE ( WaveAcc0HyiMCF (0:InitOut%NStepWave-1,NWaveKin0Prime   ), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveAcc0HyiMCF.',       ErrStat,ErrMsg,RoutineName)

         ALLOCATE ( WaveAcc0VMCF   (0:InitOut%NStepWave-1,NWaveKin0Prime   ), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveAcc0VMCF.',         ErrStat,ErrMsg,RoutineName)
      
         ALLOCATE ( InitOut%WaveAccMCF  (0:InitOut%NStepWave,InitInp%NGrid(1),InitInp%NGrid(2),InitInp%NGrid(3),3), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array InitOut%WaveAccMCF.',  ErrStat,ErrMsg,RoutineName)
      END IF
      
      
      IF (InitInp%WaveStMod .EQ. 2_IntKi) THEN ! Extrapolation Wave Stretching

         ALLOCATE ( PWaveDynPC0BPz0   (0:InitOut%NStepWave2 ,InitInp%NWaveElevGrid), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array PWaveDynPC0BPz0.',   ErrStat,ErrMsg,RoutineName)
      
         ALLOCATE ( PWaveVelC0HxiPz0  (0:InitOut%NStepWave2 ,InitInp%NWaveElevGrid), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array PWaveVelC0HxiPz0.',  ErrStat,ErrMsg,RoutineName)
      
         ALLOCATE ( PWaveVelC0HyiPz0  (0:InitOut%NStepWave2 ,InitInp%NWaveElevGrid), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array PWaveVelC0HyiPz0.',  ErrStat,ErrMsg,RoutineName)
      
         ALLOCATE ( PWaveVelC0VPz0    (0:InitOut%NStepWave2 ,InitInp%NWaveElevGrid), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array PWaveVelC0VPz0.',    ErrStat,ErrMsg,RoutineName)
      
         ALLOCATE ( PWaveAccC0HxiPz0  (0:InitOut%NStepWave2 ,InitInp%NWaveElevGrid), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array PWaveAccC0HxiPz0.',  ErrStat,ErrMsg,RoutineName)
      
         ALLOCATE ( PWaveAccC0HyiPz0  (0:InitOut%NStepWave2 ,InitInp%NWaveElevGrid), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array PWaveAccC0HyiPz0.',  ErrStat,ErrMsg,RoutineName)
      
         ALLOCATE ( PWaveAccC0VPz0    (0:InitOut%NStepWave2 ,InitInp%NWaveElevGrid), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array PWaveAccC0VPz0.',    ErrStat,ErrMsg,RoutineName)

         ALLOCATE ( PWaveDynP0BPz0    (0:InitOut%NStepWave-1,InitInp%NWaveElevGrid), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array PWaveDynP0BPz0.',    ErrStat,ErrMsg,RoutineName)
      
         ALLOCATE ( PWaveVel0HxiPz0   (0:InitOut%NStepWave-1,InitInp%NWaveElevGrid), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array PWaveVel0HxiPz0.',   ErrStat,ErrMsg,RoutineName)
      
         ALLOCATE ( PWaveVel0HyiPz0   (0:InitOut%NStepWave-1,InitInp%NWaveElevGrid), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array PWaveVel0HyiPz0.',   ErrStat,ErrMsg,RoutineName)
      
         ALLOCATE ( PWaveVel0VPz0     (0:InitOut%NStepWave-1,InitInp%NWaveElevGrid), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array PWaveVel0Pz0.',      ErrStat,ErrMsg,RoutineName)
      
         ALLOCATE ( PWaveAcc0HxiPz0   (0:InitOut%NStepWave-1,InitInp%NWaveElevGrid), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array PWaveAcc0HxiPz0.',   ErrStat,ErrMsg,RoutineName)
      
         ALLOCATE ( PWaveAcc0HyiPz0   (0:InitOut%NStepWave-1,InitInp%NWaveElevGrid), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array PWaveAcc0HyiPz0.',   ErrStat,ErrMsg,RoutineName)
      
         ALLOCATE ( PWaveAcc0VPz0     (0:InitOut%NStepWave-1,InitInp%NWaveElevGrid), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array PWaveAcc0VPz0.',     ErrStat,ErrMsg,RoutineName)

         ALLOCATE ( InitOut%PWaveDynP0 (0:InitOut%NStepWave,InitInp%NGrid(1),InitInp%NGrid(2)  ), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array InitOut%PWaveDynP0.', ErrStat,ErrMsg,RoutineName)
      
         ALLOCATE ( InitOut%PWaveVel0  (0:InitOut%NStepWave,InitInp%NGrid(1),InitInp%NGrid(2),3), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array InitOut%PWaveVel0.',  ErrStat,ErrMsg,RoutineName)
      
         ALLOCATE ( InitOut%PWaveAcc0  (0:InitOut%NStepWave,InitInp%NGrid(1),InitInp%NGrid(2),3), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array InitOut%PWaveAcc0.',  ErrStat,ErrMsg,RoutineName)
         
         IF (InitInp%MCFD > 0.0_ReKi) THEN ! MacCamy-Fuchs model
         
            ALLOCATE ( PWaveAccC0HxiMCFPz0  (0:InitOut%NStepWave2 ,InitInp%NWaveElevGrid), STAT=ErrStatTmp )
            IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array PWaveAccC0HxiMCFPz0.',  ErrStat,ErrMsg,RoutineName)
      
            ALLOCATE ( PWaveAccC0HyiMCFPz0  (0:InitOut%NStepWave2 ,InitInp%NWaveElevGrid), STAT=ErrStatTmp )
            IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array PWaveAccC0HyiMCFPz0.',  ErrStat,ErrMsg,RoutineName)
      
            ALLOCATE ( PWaveAccC0VMCFPz0    (0:InitOut%NStepWave2 ,InitInp%NWaveElevGrid), STAT=ErrStatTmp )
            IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array PWaveAccC0VMCFPz0.',    ErrStat,ErrMsg,RoutineName)

            ALLOCATE ( PWaveAcc0HxiMCFPz0   (0:InitOut%NStepWave-1,InitInp%NWaveElevGrid), STAT=ErrStatTmp )
            IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array PWaveAcc0HxiMCFPz0.',   ErrStat,ErrMsg,RoutineName)
      
            ALLOCATE ( PWaveAcc0HyiMCFPz0   (0:InitOut%NStepWave-1,InitInp%NWaveElevGrid), STAT=ErrStatTmp )
            IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array PWaveAcc0HyiMCFPz0.',   ErrStat,ErrMsg,RoutineName)
      
            ALLOCATE ( PWaveAcc0VMCFPz0     (0:InitOut%NStepWave-1,InitInp%NWaveElevGrid), STAT=ErrStatTmp )
            IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array PWaveAcc0VMCFPz0.',     ErrStat,ErrMsg,RoutineName)
         
            ALLOCATE ( InitOut%PWaveAccMCF0  (0:InitOut%NStepWave,InitInp%NGrid(1),InitInp%NGrid(2),3), STAT=ErrStatTmp )
            IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array InitOut%PWaveAccMCF0.',  ErrStat,ErrMsg,RoutineName)
         
         END IF

      END IF

! END TODO SECTION


      ! Arrays for the Sin and Cos of the wave direction for each frequency.  Used in calculating wave elevation, velocity, acceleration etc.
      ALLOCATE ( CosWaveDir( 0:InitOut%NStepWave2                          ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array CosWaveDir.',        ErrStat,ErrMsg,RoutineName)

      ALLOCATE ( SinWaveDir( 0:InitOut%NStepWave2                          ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array SinWaveDir.',        ErrStat,ErrMsg,RoutineName)

      ALLOCATE ( OmegaArr( 0:InitOut%NStepWave2                            ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array OmegaArr.',   ErrStat,ErrMsg,RoutineName)

      
      ! Arrays for the constrained wave
      ALLOCATE ( WaveS1SddArr( 0:InitOut%NStepWave2                        ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveS1SddArr.',   ErrStat,ErrMsg,RoutineName)

      ! Now check if all the allocations worked properly
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL CleanUp()
         RETURN
      END IF



      ! Compute the positive-frequency components (including zero) of the discrete
      !   Fourier transforms of the wave kinematics:
      DO I = 0,InitOut%NStepWave2  ! Loop through the positive frequency components (including zero) of the discrete Fourier transforms
          OmegaArr(I) = I*InitOut%WaveDOmega
      END DO

      call Get_1Spsd_and_WaveElevC0(InitInp, InitOut, OmegaArr, WaveS1SddArr)

      
      !> #  Multi Directional Waves
      call CalculateWaveDirection(InitInp, InitOut, ErrStatTmp, ErrMsgTmp)
         call SetErrStat(ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName)
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF
         
         ! Store the minimum and maximum wave directions
      InitOut%WaveDirMin   = MINVAL(InitOut%WaveDirArr)
      InitOut%WaveDirMax   = MAXVAL(InitOut%WaveDirArr)
         

      ! Set the CosWaveDir and SinWaveDir arrays
      CosWaveDir=COS(D2R*InitOut%WaveDirArr)
      SinWaveDir=SIN(D2R*InitOut%WaveDirArr)

      
      ! make sure this is called before calling ConstrainedNewWaves
      CALL InitFFT ( InitOut%NStepWave, FFT_Data, .TRUE., ErrStatTmp )
         CALL SetErrStat(ErrStatTmp,'Error occured while initializing the FFT.',ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF
      
      !--------------------------------------------------------------------------------
      !=== Constrained New Waves ===
      ! Modify the wave components to implement the constrained wave
      ! Only do this if WaveMod = 2 (JONSWAP/Pierson-Moskowitz Spectrum) and ConstWaveMod > 0
      IF ( InitInp%WaveMod == 2 .AND. InitInp%ConstWaveMod > 0) THEN
         ! adjust InitOut%WaveElevC0 for constrained wave:
         call ConstrainedNewWaves(InitInp, InitOut, OmegaArr, WaveS1SddArr, CosWaveDir, SinWaveDir, FFT_Data, ErrStatTmp, ErrMsgTmp)
            call SetErrStat(ErrStatTmp,ErrMsgTmp, ErrStat,ErrMsg,RoutineName)
            if (ErrStat >= AbortErrLev) then
               call cleanup()
               return
            end if
      ENDIF
      ! End of Constrained Wave

      !--------------------------------------------------------------------------------
      !> ## Phase shift the discrete Fourier transform of wave elevations at the WRP
      !> This changes the phasing of all wave kinematics and loads to reflect the turbine's
      !! location in the larger farm, in the case of FAST.Farm simulations, based on
      !! specified PtfmLocationX and PtfmLocationY.
      
      IF (InitInp%WaveFieldMod == 2) THEN             ! case 2: adjust wave phases based on turbine offsets from farm origin
      
         CALL WrScr ( ' Adjusting incident wave kinematics for turbine offset from array origin.' )
      
         DO I = 0,InitOut%NStepWave2  

            tmpComplex  = CMPLX(  InitOut%WaveElevC0(1,I),   InitOut%WaveElevC0(2,I))
            
            ! some redundant calculations with later, but insignificant
            WaveNmbr   = WaveNumber ( OmegaArr(I), InitInp%Gravity, InitInp%WtrDpth )
            
            ! apply the phase shift
            tmpComplex = tmpComplex * EXP( -ImagNmbr*WaveNmbr*( InitInp%PtfmLocationX*CosWaveDir(I) + InitInp%PtfmLocationY*SinWaveDir(I) ))
      
            ! put shifted complex amplitudes back into the array for use in the remainder of this module and other modules (Waves2, WAMIT, WAMIT2)
            InitOut%WaveElevC0 (1,I) = REAL( tmpComplex)
            InitOut%WaveElevC0 (2,I) = AIMAG(tmpComplex)
      
         END DO
      END IF


      !--------------------------------------------------------------------------------
      !> ## Compute IFFTs
      !> Compute the discrete Fourier transform of the instantaneous elevation of
      !!   incident waves at each desired point on the still water level plane
      !!   where it can be output:

      DO I = 0,InitOut%NStepWave2  ! Loop through the positive frequency components (including zero) of the discrete Fourier transforms


            ! Set tmpComplex to the Ith element of the WAveElevC0 array

         tmpComplex  = CMPLX(  InitOut%WaveElevC0(1,I),   InitOut%WaveElevC0(2,I))


      ! Compute the frequency of this component and its imaginary value:

         ImagOmega = ImagNmbr*OmegaArr(I)

      ! Compute the wavenumber:

         WaveNmbr   = WaveNumber ( OmegaArr(I), InitInp%Gravity, InitInp%WtrDpth )

      ! Wavenumber-dependent acceleration scaling for MacCamy-Fuchs model
      MCFC = 0.0_ReKi
      IF (InitInp%MCFD > 0.0_SiKi .AND. I>0_IntKi) THEN
         ka = 0.5_ReKi * WaveNmbr * InitInp%MCFD
         JPrime = BESSEL_JN(1,ka) / ka - BESSEL_JN(2,ka)
         YPrime = BESSEL_YN(1,ka) / ka - BESSEL_YN(2,ka)
         HPrime = SQRT(JPrime*JPrime + YPrime*YPrime)
         MCFC = 4.0_ReKi/( PI * ka * ka * HPrime )
      END IF

      ! Compute the discrete Fourier transform of the incident wave kinematics
      !   before applying stretching at the zi-coordinates for the WAMIT reference point, and all
      !   points where are Morison loads will be calculated.

         DO J = 1,NWaveKin0Prime ! Loop through all points where the incident wave kinematics will be computed without stretching

            WaveElevxiPrime0 = EXP( -ImagNmbr*WaveNmbr*( InitInp%WaveKinGridxi(WaveKinPrimeMap(J))*CosWaveDir(I) + &
                                                         InitInp%WaveKinGridyi(WaveKinPrimeMap(J))*SinWaveDir(I) ))

            WaveDynPC0 (I,J)     = InitOut%RhoXg*tmpComplex*WaveElevxiPrime0 * COSHNumOvrCOSHDen ( WaveNmbr, InitInp%WtrDpth, WaveKinzi0Prime(J) )

            WaveVelC0Hxi (I,J)   = CosWaveDir(I)*OmegaArr(I)*tmpComplex* WaveElevxiPrime0 * COSHNumOvrSINHDen ( WaveNmbr, InitInp%WtrDpth, WaveKinzi0Prime(J) )
            WaveVelC0Hyi (I,J)   = SinWaveDir(I)*OmegaArr(I)*tmpComplex* WaveElevxiPrime0 * COSHNumOvrSINHDen ( WaveNmbr, InitInp%WtrDpth, WaveKinzi0Prime(J) )

            WaveVelC0V (I,J)     = ImagOmega*tmpComplex* WaveElevxiPrime0 * SINHNumOvrSINHDen ( WaveNmbr, InitInp%WtrDpth, WaveKinzi0Prime(J) )
            WaveAccC0Hxi (I,J)   = ImagOmega*        WaveVelC0Hxi (I,J)

            WaveAccC0Hyi (I,J)   = ImagOmega*        WaveVelC0Hyi (I,J)
            WaveAccC0V (I,J)     = ImagOmega*        WaveVelC0V   (I,J)

            IF (InitInp%MCFD > 0.0_SiKi) THEN
               WaveAccC0HxiMCF(I,J) = WaveAccC0Hxi(I,J) * MCFC
               WaveAccC0HyiMCF(I,J) = WaveAccC0Hyi(I,J) * MCFC
               WaveAccC0VMCF(I,J)   = WaveAccC0V(I,J)   * MCFC
            END IF


         END DO                   ! J - All points where the incident wave kinematics will be computed without stretching

         !===================================
         IF (InitInp%WaveStMod .EQ. 2_IntKi) THEN ! Extrapolation wave stretching
            DO J = 1,InitInp%NWaveElevGrid ! Loop through all points on the SWL
               WaveElevxiPrime0 = EXP( -ImagNmbr*WaveNmbr*( InitInp%WaveKinGridxi(J)*CosWaveDir(I) + &
                                                            InitInp%WaveKinGridyi(J)*SinWaveDir(I) ))
               ! Partial derivatives at zi = 0
               PWaveDynPC0BPz0 (I,J) = InitOut%RhoXg*      tmpComplex*WaveElevxiPrime0*WaveNmbr*TANH ( WaveNmbr*InitInp%WtrDpth )
               PWaveVelC0HxiPz0(I,J) = CosWaveDir(I)*OmegaArr(I)*tmpComplex*WaveElevxiPrime0*WaveNmbr
               PWaveVelC0HyiPz0(I,J) = SinWaveDir(I)*OmegaArr(I)*tmpComplex*WaveElevxiPrime0*WaveNmbr
            
               IF (I == 0_IntKi) THEN ! Zero frequency component - Need to avoid division by zero.
                 PWaveVelC0VPz0  (I,J) =         0.0_ReKi
               ELSE
                 PWaveVelC0VPz0  (I,J) =         ImagOmega*tmpComplex*WaveElevxiPrime0*WaveNmbr/TANH ( WaveNmbr*InitInp%WtrDpth )
               END IF
            
               PWaveAccC0HxiPz0(I,J) =           ImagOmega*PWaveVelC0HxiPz0(I,J)
               PWaveAccC0HyiPz0(I,J) =           ImagOmega*PWaveVelC0HyiPz0(I,J)
               PWaveAccC0VPz0  (I,J) =           ImagOmega*PWaveVelC0VPz0  (I,J)
               
               
               IF (InitInp%MCFD > 0.0_SiKi) THEN
                  PWaveAccC0HxiMCFPz0(I,J) = PWaveAccC0HxiPz0(I,J) * MCFC
                  PWaveAccC0HyiMCFPz0(I,J) = PWaveAccC0HyiPz0(I,J) * MCFC
                  PWaveAccC0VMCFPz0(I,J)   = PWaveAccC0VPz0(I,J)   * MCFC
               END IF
               
            END DO                   ! J - All points where the incident wave kinematics will be computed without stretching
         END IF
        !===================================

      END DO                ! I - The positive frequency components (including zero) of the discrete Fourier transforms

      ! Calculate the array of simulation times at which the instantaneous
      !   elevation of, velocity of, acceleration of, and loads associated with
      !   the incident waves are to be determined:
      DO I = 0,InitOut%NStepWave ! Loop through all time steps
         InitOut%WaveTime(I) = I*REAL(InitInp%WaveDT,SiKi)
      END DO                ! I - All time steps
      
      
      DO I = 0,InitOut%NStepWave2  ! Loop through the positive frequency components (including zero) of the discrete Fourier transform
         tmpComplexArr(I)    =  CMPLX(InitOut%WaveElevC0(1,I), InitOut%WaveElevC0(2,I))
      END DO

      ! Compute the inverse discrete Fourier transforms to find the time-domain
      !   representations of the wave kinematics without stretcing:

      CALL    ApplyFFT_cx (  InitOut%WaveElev0    (0:InitOut%NStepWave-1),  tmpComplexArr    (:  ), FFT_Data, ErrStatTmp )
      CALL SetErrStat(ErrStatTmp,'Error occured while applying the FFT to WaveElev0.',ErrStat,ErrMsg,RoutineName)
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL CleanUp()
         RETURN
      END IF
!NOTE:  For all grid points
      DO k = 1,InitInp%NWaveElevGrid     ! Loop through all points where the incident wave elevations are to be computed (normally all the XY grid points)
               ! This subroutine call applies the FFT at the correct location.
         i = mod(k-1, InitInp%NGrid(1)) + 1
         j = (k-1) / InitInp%NGrid(1) + 1
            ! note that this subroutine resets tmpComplexArr
         CALL WaveElevTimeSeriesAtXY( InitInp%WaveKinGridxi(k), InitInp%WaveKinGridyi(k), InitOut%WaveElev(:,i,j), InitOut%WaveElevC(:,:,k), tmpComplexArr, ErrStatTmp, ErrMsgTmp ) ! Note this sets tmpComplexArr
         CALL SetErrStat(ErrStatTmp,'Error occured while applying the FFT to InitOut%WaveElev.',ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF
      END DO                   ! J - All points where the incident wave elevations can be output




         ! User requested data points -- Do all the FFT calls first, then return if something failed.
      DO J = 1,NWaveKin0Prime ! Loop through all points where the incident wave kinematics will be computed without stretching
         CALL ApplyFFT_cx (          WaveDynP0B   (:,J),          WaveDynPC0    (:,J), FFT_Data, ErrStatTmp )
         CALL SetErrStat(ErrStatTmp,'Error occured while applying the FFT to WaveDynP0B.',       ErrStat,ErrMsg,RoutineName)

         CALL ApplyFFT_cx (          WaveVel0Hxi  (:,J),          WaveVelC0Hxi  (:,J), FFT_Data, ErrStatTmp )
         CALL SetErrStat(ErrStatTmp,'Error occured while applying the FFT to WaveVel0Hxi.',      ErrStat,ErrMsg,RoutineName)

         CALL ApplyFFT_cx (          WaveVel0Hyi  (:,J),          WaveVelC0Hyi  (:,J), FFT_Data, ErrStatTmp )
         CALL SetErrStat(ErrStatTmp,'Error occured while applying the FFT to WaveVel0Hyi.',      ErrStat,ErrMsg,RoutineName)

         CALL ApplyFFT_cx (          WaveVel0V    (:,J),          WaveVelC0V    (:,J), FFT_Data, ErrStatTmp )
         CALL SetErrStat(ErrStatTmp,'Error occured while applying the FFT to WaveVel0V.',        ErrStat,ErrMsg,RoutineName)

         CALL ApplyFFT_cx (          WaveAcc0Hxi  (:,J),          WaveAccC0Hxi  (:,J), FFT_Data, ErrStatTmp )
         CALL SetErrStat(ErrStatTmp,'Error occured while applying the FFT to WaveAcc0Hxi.',      ErrStat,ErrMsg,RoutineName)

         CALL ApplyFFT_cx (          WaveAcc0Hyi  (:,J),          WaveAccC0Hyi  (:,J), FFT_Data, ErrStatTmp )
         CALL SetErrStat(ErrStatTmp,'Error occured while applying the FFT to WaveAcc0Hyi.',      ErrStat,ErrMsg,RoutineName)

         CALL ApplyFFT_cx (          WaveAcc0V    (:,J),          WaveAccC0V    (:,J), FFT_Data, ErrStatTmp )
         CALL SetErrStat(ErrStatTmp,'Error occured while applying the FFT to WaveAcc0V.',        ErrStat,ErrMsg,RoutineName)

         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF

      END DO                   ! J - All points where the incident wave kinematics will be computed without stretching

      IF (InitInp%MCFD > 0.0_SiKi) THEN
         DO J = 1,NWaveKin0Prime ! Loop through all points where the incident wave kinematics will be computed without stretching
            CALL ApplyFFT_cx (          WaveAcc0HxiMCF  (:,J),          WaveAccC0HxiMCF  (:,J), FFT_Data, ErrStatTmp )
            CALL SetErrStat(ErrStatTmp,'Error occured while applying the FFT to WaveAcc0HxiMCF.',      ErrStat,ErrMsg,RoutineName)

            CALL ApplyFFT_cx (          WaveAcc0HyiMCF  (:,J),          WaveAccC0HyiMCF  (:,J), FFT_Data, ErrStatTmp )
            CALL SetErrStat(ErrStatTmp,'Error occured while applying the FFT to WaveAcc0HyiMCF.',      ErrStat,ErrMsg,RoutineName)

            CALL ApplyFFT_cx (          WaveAcc0VMCF    (:,J),          WaveAccC0VMCF    (:,J), FFT_Data, ErrStatTmp )
            CALL SetErrStat(ErrStatTmp,'Error occured while applying the FFT to WaveAcc0VMCF.',        ErrStat,ErrMsg,RoutineName)
        
            IF ( ErrStat >= AbortErrLev ) THEN
               CALL CleanUp()
               RETURN
            END IF
         END DO
      END IF

      !===================================
      IF (InitInp%WaveStMod .EQ. 2_IntKi) THEN ! Extrapolation Wave Stretching
         DO J = 1,InitInp%NWaveElevGrid ! Loop through all points on the SWL where z-partial derivatives will be computed for extrapolated stretching
            ! FFT's of the partial derivatives
            CALL  ApplyFFT_cx (         PWaveDynP0BPz0(:,J  ),         PWaveDynPC0BPz0(:,J  ), FFT_Data, ErrStatTmp )
            CALL  SetErrStat(ErrStatTmp,'Error occured while applying the FFT to PWaveDynP0BPz0.',   ErrStat,ErrMsg,RoutineName)
      
            CALL  ApplyFFT_cx (         PWaveVel0HxiPz0 (:,J  ),       PWaveVelC0HxiPz0( :,J ),FFT_Data, ErrStatTmp )
            CALL  SetErrStat(ErrStatTmp,'Error occured while applying the FFT to PWaveVel0HxiPz0.',  ErrStat,ErrMsg,RoutineName)
      
            CALL  ApplyFFT_cx (         PWaveVel0HyiPz0 (:,J  ),       PWaveVelC0HyiPz0( :,J ),FFT_Data, ErrStatTmp )
            CALL  SetErrStat(ErrStatTmp,'Error occured while applying the FFT to PWaveVel0HyiPz0.',  ErrStat,ErrMsg,RoutineName)
      
            CALL  ApplyFFT_cx (         PWaveVel0VPz0 (:,J  ),         PWaveVelC0VPz0 (:,J  ), FFT_Data, ErrStatTmp )
            CALL  SetErrStat(ErrStatTmp,'Error occured while applying the FFT to PWaveVel0VPz0.',    ErrStat,ErrMsg,RoutineName)
      
            CALL  ApplyFFT_cx (         PWaveAcc0HxiPz0 (:,J  ),       PWaveAccC0HxiPz0(:,J  ),FFT_Data, ErrStatTmp )
            CALL  SetErrStat(ErrStatTmp,'Error occured while applying the FFT to PWaveAcc0HxiPz0.',  ErrStat,ErrMsg,RoutineName)
      
            CALL  ApplyFFT_cx (         PWaveAcc0HyiPz0 (:,J  ),       PWaveAccC0HyiPz0(:,J  ),FFT_Data, ErrStatTmp )
            CALL  SetErrStat(ErrStatTmp,'Error occured while applying the FFT to PWaveAcc0HyiPz0.',  ErrStat,ErrMsg,RoutineName)
      
            CALL  ApplyFFT_cx (         PWaveAcc0VPz0 (:,J  ),         PWaveAccC0VPz0( :,J  ), FFT_Data, ErrStatTmp )
            CALL  SetErrStat(ErrStatTmp,'Error occured while applying the FFT to PWaveAcc0VPz0.',    ErrStat,ErrMsg,RoutineName)
      
            IF ( ErrStat >= AbortErrLev ) THEN
               CALL CleanUp()
               RETURN
            END IF
      
         END DO                   ! J - All points where the incident wave kinematics will be computed without stretching
         
         IF (InitInp%MCFD > 0.0_SiKi) THEN ! MacCamy-Fuchs scaled acceleration field
            DO J = 1,InitInp%NWaveElevGrid
      
               CALL  ApplyFFT_cx (         PWaveAcc0HxiMCFPz0 (:,J  ),       PWaveAccC0HxiMCFPz0(:,J  ),FFT_Data, ErrStatTmp )
               CALL  SetErrStat(ErrStatTmp,'Error occured while applying the FFT to PWaveAcc0HxiMCFPz0.',  ErrStat,ErrMsg,RoutineName)
      
               CALL  ApplyFFT_cx (         PWaveAcc0HyiMCFPz0 (:,J  ),       PWaveAccC0HyiMCFPz0(:,J  ),FFT_Data, ErrStatTmp )
               CALL  SetErrStat(ErrStatTmp,'Error occured while applying the FFT to PWaveAcc0HyiMCFPz0.',  ErrStat,ErrMsg,RoutineName)
      
               CALL  ApplyFFT_cx (         PWaveAcc0VMCFPz0 (:,J  ),         PWaveAccC0VMCFPz0( :,J  ), FFT_Data, ErrStatTmp )
               CALL  SetErrStat(ErrStatTmp,'Error occured while applying the FFT to PWaveAcc0VMCFPz0.',    ErrStat,ErrMsg,RoutineName)
            
               IF ( ErrStat >= AbortErrLev ) THEN
                  CALL CleanUp()
                  RETURN
               END IF
            
            END DO
         END IF
         
      END IF
!===================================


      CALL  ExitFFT(FFT_Data, ErrStatTmp)
      CALL  SetErrStat(ErrStatTmp,'Error occured while cleaning up after the FFTs.', ErrStat,ErrMsg,RoutineName)
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

         ! Commented out - We do not extrapolate the current profile with extrapolated wave stretching
         !PWaveVel0HxiPz0(:  ) =  PWaveVel0HxiPz0(:  ) + InitInp%PCurrVxiPz0  ! xi-direction
         !PWaveVel0HyiPz0(:  ) =  PWaveVel0HyiPz0(:  ) + InitInp%PCurrVyiPz0  ! yi-direction

      ENDIF


      ! Apply stretching to obtain the wave kinematics, WaveDynP0, WaveVel0, and
      !   WaveAcc0, at the desired locations from the wave kinematics at
      !   alternative locations, WaveDynP0B, WaveVel0Hxi, WaveVel0Hyi, WaveVel0V,
      !   WaveAcc0Hxi, WaveAcc0Hyi, WaveAcc0V, if the elevation of the point defined by
      !   WaveKinGridzi(J) lies between the seabed and the instantaneous free
      !   surface, else set WaveDynP0, WaveVel0, and WaveAcc0 to zero.  This
      !   depends on which incident wave kinematics stretching method is being
      !   used:

    !  SELECT CASE ( InitInp%WaveStMod )  ! Which model are we using to extrapolate the incident wave kinematics to the instantaneous free surface?

    !  CASE ( 0 )                 ! None=no stretching.


      ! Since we have no stretching, the wave kinematics between the seabed and
      !   the mean sea level are left unchanged; below the seabed or above the
      !   mean sea level, the wave kinematics are zero:

   !   InitOut%PWaveDynP0(:,:,:,:)   = 0.0
   !   InitOut%PWaveVel0 (:,:,:,:,:) = 0.0
   !   InitOut%PWaveAcc0 (:,:,:,:,:) = 0.0

      primeCount = 1
      count = 1
      !DO J = 1,InitInp%NWaveKinGrid      ! Loop through all points where the incident wave kinematics will be computed
      do k = 1, InitInp%NGrid(3)
         do j = 1, InitInp%NGrid(2)
            do i = 1, InitInp%NGrid(1)

             !  ii = mod(count-1, InitInp%NGrid(1)) + 1
             !  jj = mod( (count-1) /InitInp%NGrid(1), InitInp%NGrid(2) ) + 1
             !  kk = (count-1) / (InitInp%NGrid(1)*InitInp%NGrid(2)) + 1

               IF (   ( InitInp%WaveKinGridzi(count) < -InitInp%WtrDpth ) .OR. ( InitInp%WaveKinGridzi(count) > 0.0 ) ) THEN
                  ! .TRUE. if the elevation of the point defined by WaveKinGridzi(J) lies below the seabed or above mean sea level (exclusive)
                  ! NOTE: We test to 0 instead of MSL2SWL because the locations of WaveKinGridzi and WtrDpth have already been adjusted using MSL2SWL

                  InitOut%WaveDynP(:,i,j,k  )  = 0.0
                  InitOut%WaveVel (:,i,j,k,:)  = 0.0
                  InitOut%WaveAcc (:,i,j,k,:)  = 0.0

               ELSE
                  ! The elevation of the point defined by WaveKinGridzi(J) must lie between the seabed and the mean sea level (inclusive)

                  InitOut%WaveDynP(0:InitOut%NStepWave-1,i,j,k  ) = WaveDynP0B( 0:InitOut%NStepWave-1,primeCount)
                  InitOut%WaveVel (0:InitOut%NStepWave-1,i,j,k,1) = WaveVel0Hxi(0:InitOut%NStepWave-1,primeCount)
                  InitOut%WaveVel (0:InitOut%NStepWave-1,i,j,k,2) = WaveVel0Hyi(0:InitOut%NStepWave-1,primeCount)
                  InitOut%WaveVel (0:InitOut%NStepWave-1,i,j,k,3) = WaveVel0V(  0:InitOut%NStepWave-1,primeCount)
                  InitOut%WaveAcc (0:InitOut%NStepWave-1,i,j,k,1) = WaveAcc0Hxi(0:InitOut%NStepWave-1,primeCount)
                  InitOut%WaveAcc (0:InitOut%NStepWave-1,i,j,k,2) = WaveAcc0Hyi(0:InitOut%NStepWave-1,primeCount)
                  InitOut%WaveAcc (0:InitOut%NStepWave-1,i,j,k,3) = WaveAcc0V(  0:InitOut%NStepWave-1,primeCount)
                  primeCount = primeCount + 1
               END IF
               count = count + 1
            end do
         end do
      end do

      ! MacCamy-Fuchs scaled fluid acceleration
      IF (InitInp%MCFD > 0.0_SiKi) THEN
         primeCount = 1
         count = 1
         do k = 1, InitInp%NGrid(3)
            do j = 1, InitInp%NGrid(2)
               do i = 1, InitInp%NGrid(1)
                  IF (   ( InitInp%WaveKinGridzi(count) < -InitInp%WtrDpth ) .OR. ( InitInp%WaveKinGridzi(count) > 0.0 ) ) THEN
                     ! .TRUE. if the elevation of the point defined by WaveKinGridzi(J) lies below the seabed or above mean sea level (exclusive)
                     ! NOTE: We test to 0 instead of MSL2SWL because the locations of WaveKinGridzi and WtrDpth have already been adjusted using MSL2SWL
                     InitOut%WaveAccMCF(:,i,j,k,:)  = 0.0
                  ELSE
                     ! The elevation of the point defined by WaveKinGridzi(J) must lie between the seabed and the mean sea level (inclusive)
                     InitOut%WaveAccMCF (0:InitOut%NStepWave-1,i,j,k,1) = WaveAcc0HxiMCF(0:InitOut%NStepWave-1,primeCount)
                     InitOut%WaveAccMCF (0:InitOut%NStepWave-1,i,j,k,2) = WaveAcc0HyiMCF(0:InitOut%NStepWave-1,primeCount)
                     InitOut%WaveAccMCF (0:InitOut%NStepWave-1,i,j,k,3) = WaveAcc0VMCF(  0:InitOut%NStepWave-1,primeCount)
                     primeCount = primeCount + 1
                  END IF
                  count = count + 1
               end do
            end do
         end do
      END IF

      IF (InitInp%WaveStMod .EQ. 2_IntKi) THEN ! Extrapolation Wave Stretching
         
         primeCount = 1
         DO j = 1, InitInp%NGrid(2)  ! Loop through all points on the SWL where partial derivatives about z were computed
            DO i = 1, InitInp%NGrid(1)
               InitOut%PWaveDynP0(0:InitOut%NStepWave-1,i,j  ) = PWaveDynP0BPz0( 0:InitOut%NStepWave-1,primeCount)
               InitOut%PWaveVel0 (0:InitOut%NStepWave-1,i,j,1) = PWaveVel0HxiPz0(0:InitOut%NStepWave-1,primeCount)
               InitOut%PWaveVel0 (0:InitOut%NStepWave-1,i,j,2) = PWaveVel0HyiPz0(0:InitOut%NStepWave-1,primeCount)
               InitOut%PWaveVel0 (0:InitOut%NStepWave-1,i,j,3) = PWaveVel0VPz0(  0:InitOut%NStepWave-1,primeCount)
               InitOut%PWaveAcc0 (0:InitOut%NStepWave-1,i,j,1) = pWaveAcc0HxiPz0(0:InitOut%NStepWave-1,primeCount)
               InitOut%PWaveAcc0 (0:InitOut%NStepWave-1,i,j,2) = pWaveAcc0HyiPz0(0:InitOut%NStepWave-1,primeCount)
               InitOut%PWaveAcc0 (0:InitOut%NStepWave-1,i,j,3) = PWaveAcc0VPz0(  0:InitOut%NStepWave-1,primeCount)
               primeCount = primeCount + 1
            END DO
         END DO
         
         IF (InitInp%MCFD > 0.0_SiKi) THEN
            primeCount = 1
            DO j = 1, InitInp%NGrid(2)  ! Loop through all points on the SWL where partial derivatives about z were computed
               DO i = 1, InitInp%NGrid(1)
                  InitOut%PWaveAccMCF0 (0:InitOut%NStepWave-1,i,j,1) = pWaveAcc0HxiMCFPz0(0:InitOut%NStepWave-1,primeCount)
                  InitOut%PWaveAccMCF0 (0:InitOut%NStepWave-1,i,j,2) = pWaveAcc0HyiMCFPz0(0:InitOut%NStepWave-1,primeCount)
                  InitOut%PWaveAccMCF0 (0:InitOut%NStepWave-1,i,j,3) = PWaveAcc0VMCFPz0(  0:InitOut%NStepWave-1,primeCount)
                  primeCount = primeCount + 1
               END DO
            END DO
         END IF

      END IF


    !  END DO                   ! J - All points where the incident wave kinematics will be computed

    !  CASE ( 1 )                 ! Vertical stretching.


      ! Vertical stretching says that the wave kinematics above the mean sea level
      !   equal the wave kinematics at the mean sea level.  The wave kinematics
      !   below the mean sea level are left unchanged:





    !  CASE ( 2 )                 ! Extrapolation stretching.


      ! Extrapolation stretching uses a linear Taylor expansion of the wave
      !   kinematics (and their partial derivatives with respect to z) at the mean
      !   sea level to find the wave kinematics above the mean sea level.  The
      !   wave kinematics below the mean sea level are left unchanged:





    !  CASE ( 3 )                 ! Wheeler stretching.


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
      !   WaveKinGridzi(:) array) using the WaveKinzi0St(:) array:




  !    ENDSELECT

      ! Set the ending timestep to the same as the first timestep
      InitOut%WaveElev0 (InitOut%NStepWave)          = InitOut%WaveElev0 (0    )
      InitOut%WaveDynP  (InitOut%NStepWave,:,:,:  )  = InitOut%WaveDynP  (0,:,:,:  )
      InitOut%WaveVel   (InitOut%NStepWave,:,:,:,:)  = InitOut%WaveVel   (0,:,:,:,:)
      InitOut%WaveAcc   (InitOut%NStepWave,:,:,:,:)  = InitOut%WaveAcc   (0,:,:,:,:)
      IF (InitInp%MCFD > 0.0_SiKi) THEN
         InitOut%WaveAccMCF (InitOut%NStepWave,:,:,:,:) = InitOut%WaveAccMCF(0,:,:,:,:)
      END IF
      
      IF (InitInp%WaveStMod .EQ. 2_IntKi) THEN ! Extrapolation Wave Stretching
         InitOut%PWaveDynP0(InitOut%NStepWave,:,:  )    = InitOut%PWaveDynP0(0,:,:  )
         InitOut%PWaveVel0 (InitOut%NStepWave,:,:,:)    = InitOut%PWaveVel0 (0,:,:,:)
         InitOut%PWaveAcc0 (InitOut%NStepWave,:,:,:)    = InitOut%PWaveAcc0 (0,:,:,:)
         IF (InitInp%MCFD > 0.0_SiKi) THEN
            InitOut%PWaveAccMCF0 (InitOut%NStepWave,:,:,:) = InitOut%PWaveAccMCF0(0,:,:,:)
         END IF
      END IF

   CALL CleanUp ( )


CONTAINS

!--------------------------------------------------------------------------------
   SUBROUTINE WaveElevTimeSeriesAtXY(Xcoord,Ycoord, WaveElevAtXY, WaveElevCAtXY, tmpComplexArr, ErrStatLcl, ErrMsgLcl )

      REAL(SiKi),       INTENT(IN   )                 :: Xcoord
      REAL(SiKi),       INTENT(IN   )                 :: Ycoord
      REAL(SiKi),       INTENT(  OUT)                 :: WaveElevAtXY(0:InitOut%NStepWave)
      real(SiKi),       INTENT(  OUT)                 :: WaveElevCAtXY(2,0:InitOut%NStepWave2)
      COMPLEX(SiKi),    INTENT(INOUT)                 :: tmpComplexArr(0:InitOut%NStepWave2)            ! A temporary array (0:NStepWave2-1) for FFT use.
      INTEGER(IntKi),   INTENT(  OUT)                 :: ErrStatLcl
      CHARACTER(*),     INTENT(  OUT)                 :: ErrMsgLcl
      
      integer                                         :: i
      REAL(SiKi)                                      :: WaveNmbr                                        ! Wavenumber of the current frequency component (1/meter)
      INTEGER(IntKi)                                  :: ErrStatLcl2

      ! note that InitOut, InitInp, FFT_Data, CosWaveDir and SinWaveDir are used here, but their values are not changed
      ErrStatLcl  = ErrID_None
      ErrMsgLcl = ""

         ! Zero out the temporary array.
      tmpComplexArr  = CMPLX(0.0_SiKi,0.0_SiKi)

         ! Loop through the positive frequency components (including zero).
      DO I = 0,InitOut%NStepWave2

         WaveNmbr          = WaveNumber ( OmegaArr(I), InitInp%Gravity, InitInp%WtrDpth )
         tmpComplexArr(I)  =  CMPLX(  InitOut%WaveElevC0(1,I),   InitOut%WaveElevC0(2,I))   *          &
                                      EXP( -ImagNmbr*WaveNmbr*(  Xcoord*CosWaveDir(I)+    &
                                                                 Ycoord*SinWaveDir(I) )   )
      ENDDO

      CALL ApplyFFT_cx (   WaveElevAtXY(0:InitOut%NStepWave-1),   tmpComplexArr, FFT_Data,   ErrStatLcl2  )
      CALL SetErrStat(ErrStatLcl2,'Error occured while applying the FFT.',ErrStatLcl,ErrMsgLcl,'WaveElevTimeSeriesAtXY')

      WaveElevCAtXY( 1,: ) = REAL(tmpComplexArr(:))
      WaveElevCAtXY( 2,: ) = AIMAG(tmpComplexArr(:))
      
         ! Append first datpoint as the last as aid for repeated wave data
      WaveElevAtXY(InitOut%NStepWave) = WaveElevAtXY(0)

   END SUBROUTINE WaveElevTimeSeriesAtXY

!--------------------------------------------------------------------------------
   SUBROUTINE CleanUp( )

      IF (ALLOCATED( WaveKinPrimeMap ))   DEALLOCATE( WaveKinPrimeMap,  STAT=ErrStatTmp)
      IF (ALLOCATED( WaveKinzi0Prime ))   DEALLOCATE( WaveKinzi0Prime,  STAT=ErrStatTmp)
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
      IF (ALLOCATED( WaveVel0Hxi ))       DEALLOCATE( WaveVel0Hxi,      STAT=ErrStatTmp)
      IF (ALLOCATED( WaveVel0Hyi ))       DEALLOCATE( WaveVel0Hyi,      STAT=ErrStatTmp)
      IF (ALLOCATED( WaveVel0V ))         DEALLOCATE( WaveVel0V,        STAT=ErrStatTmp)
      IF (ALLOCATED( WaveVelC0Hxi ))      DEALLOCATE( WaveVelC0Hxi,     STAT=ErrStatTmp)
      IF (ALLOCATED( WaveVelC0Hyi ))      DEALLOCATE( WaveVelC0Hyi,     STAT=ErrStatTmp)
      IF (ALLOCATED( WaveVelC0V ))        DEALLOCATE( WaveVelC0V,       STAT=ErrStatTmp)
      IF (ALLOCATED( tmpComplexArr ))     DEALLOCATE( tmpComplexArr,    STAT=ErrStatTmp)

      IF (ALLOCATED( WaveS1SddArr ))      DEALLOCATE( WaveS1SddArr,     STAT=ErrStatTmp)
      IF (ALLOCATED( OmegaArr ))          DEALLOCATE( OmegaArr,         STAT=ErrStatTmp)

      IF (ALLOCATED( WaveAccC0HxiMCF ))     DEALLOCATE( WaveAccC0HxiMCF,     STAT=ErrStatTmp)
      IF (ALLOCATED( WaveAccC0HyiMCF ))     DEALLOCATE( WaveAccC0HyiMCF,     STAT=ErrStatTmp)
      IF (ALLOCATED( WaveAccC0VMCF ))       DEALLOCATE( WaveAccC0VMCF,       STAT=ErrStatTmp)
      IF (ALLOCATED( WaveAcc0HxiMCF ))      DEALLOCATE( WaveAcc0HxiMCF,      STAT=ErrStatTmp)
      IF (ALLOCATED( WaveAcc0HyiMCF ))      DEALLOCATE( WaveAcc0HyiMCF,      STAT=ErrStatTmp)
      IF (ALLOCATED( WaveAcc0VMCF ))        DEALLOCATE( WaveAcc0VMCF,        STAT=ErrStatTmp)
      IF (ALLOCATED( PWaveAccC0HxiMCFPz0 )) DEALLOCATE( PWaveAccC0HxiMCFPz0, STAT=ErrStatTmp)
      IF (ALLOCATED( PWaveAccC0HyiMCFPz0 )) DEALLOCATE( PWaveAccC0HyiMCFPz0, STAT=ErrStatTmp)
      IF (ALLOCATED( PWaveAccC0VMCFPz0 ))   DEALLOCATE( PWaveAccC0VMCFPz0,   STAT=ErrStatTmp)
      IF (ALLOCATED( PWaveAcc0HxiMCFPz0 ))  DEALLOCATE( PWaveAcc0HxiMCFPz0,  STAT=ErrStatTmp)
      IF (ALLOCATED( PWaveAcc0HyiMCFPz0 ))  DEALLOCATE( PWaveAcc0HyiMCFPz0,  STAT=ErrStatTmp)
      IF (ALLOCATED( PWaveAcc0VMCFPz0 ))    DEALLOCATE( PWaveAcc0VMCFPz0,    STAT=ErrStatTmp)


      RETURN

   END SUBROUTINE CleanUp


END SUBROUTINE VariousWaves_Init




!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is called at the start of the simulation to perform initialization steps.
!! The initial states and initial guess for the input are defined.
SUBROUTINE Waves_Init( InitInp, InitOut, ErrStat, ErrMsg )
!..................................................................................................................................

      TYPE(Waves_InitInputType),       INTENT(INOUT)  :: InitInp     !< Input data for initialization routine !NOTE: We are making this INOUT because UserWaveComponents_Init changes the value of InitInp%WaveDT
      TYPE(Waves_InitOutputType),      INTENT(  OUT)  :: InitOut     !< Output for initialization routine
      INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None



         ! Local Variables:
      INTEGER(IntKi)                                  :: ErrStatTmp  ! Temporary error status for processing
      CHARACTER(ErrMsgLen)                            :: ErrMsgTmp   ! Temporary error message for procesing
!      REAL(ReKi), ALLOCATABLE                         :: tmpWaveKinzi(:)

!      TYPE(FFT_DataType)           :: FFT_Data                                        ! the instance of the FFT module we're using



         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrStatTmp  = ErrID_None
      ErrMsg  = ""
      ErrMsgTmp   = ""


         ! Initialize the pRNG
      CALL RandNum_Init(InitInp%RNG, ErrStat, ErrMsg)
      IF ( ErrStat >= AbortErrLev ) RETURN

         ! Define initialization-routine output here:

      !InitOut%WriteOutputHdr = (/ 'Time', 'Column2' /)
      !InitOut%WriteOutputUnt = (/ '(s)',  '(-)'     /)

      InitOut%RhoXg         = InitInp%WtrDens*InitInp%Gravity




            ! Initialize the variables associated with the incident wave:

      SELECT CASE ( InitInp%WaveMod ) ! Which incident wave kinematics model are we using?


      CASE ( 0 )              ! None=still water.

         CALL StillWaterWaves_Init( InitInp, InitOut, ErrStatTmp, ErrMsgTmp )
         CALL  SetErrStat(ErrStatTmp,ErrMsgTmp,ErrStat,ErrMsg,'Waves_Init')
         IF ( ErrStat >= AbortErrLev ) RETURN



      CASE ( 1, 2, 3, 4, 10 )       ! 1, 10: Plane progressive (regular) wave, 2: JONSWAP/Pierson-Moskowitz spectrum (irregular) wave, 3: white-noise, or 4: user-defined spectrum (irregular) wave.

            ! Now call the init with all the zi locations for the Morrison member nodes
         CALL VariousWaves_Init( InitInp, InitOut, ErrStatTmp, ErrMsgTmp )
            CALL  SetErrStat(ErrStatTmp,ErrMsgTmp,ErrStat,ErrMsg,'Waves_Init')
            IF ( ErrStat >= AbortErrLev ) RETURN


      CASE ( 5 )              ! User-supplied wave elevation time history; HD derives full wave kinematics from this elevation time series data.

            ! Get the wave frequency information from the file (by FFT of the elevation)
         CALL UserWaveElevations_Init( InitInp, InitOut, ErrStatTmp, ErrMsgTmp )
         CALL  SetErrStat(ErrStatTmp,ErrMsgTmp,ErrStat,ErrMsg,'Waves_Init')
         IF ( ErrStat >= AbortErrLev ) RETURN

            ! Now call VariousWaves to continue using the wave elevation and derived frequency information from the file
         CALL VariousWaves_Init( InitInp, InitOut, ErrStatTmp, ErrMsgTmp )
         CALL  SetErrStat(ErrStatTmp,ErrMsgTmp,ErrStat,ErrMsg,'Waves_Init')
         IF ( ErrStat >= AbortErrLev ) RETURN


      CASE ( 6 )              ! User-supplied wave kinematics data.

         CALL UserWaves_Init( InitInp, InitOut, ErrStatTmp, ErrMsgTmp )
         CALL  SetErrStat(ErrStatTmp,ErrMsgTmp,ErrStat,ErrMsg,'Waves_Init')
         IF ( ErrStat >= AbortErrLev ) RETURN
         
      CASE ( 7 )
         
         ! Get the wave frequency information from the file (by reading in wave frequency components)
         CALL UserWaveComponents_Init( InitInp, InitOut, ErrStatTmp, ErrMsgTmp )
         CALL  SetErrStat(ErrStatTmp,ErrMsgTmp,ErrStat,ErrMsg,'Waves_Init')
         IF ( ErrStat >= AbortErrLev ) RETURN

         ! Now call VariousWaves to continue using the wave frequency information from the file
         CALL VariousWaves_Init( InitInp, InitOut, ErrStatTmp, ErrMsgTmp )
         CALL  SetErrStat(ErrStatTmp,ErrMsgTmp,ErrStat,ErrMsg,'Waves_Init')
         IF ( ErrStat >= AbortErrLev ) RETURN 
         
      ENDSELECT

END SUBROUTINE Waves_Init



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
   REAL(SiKi)                     :: WheelerStretching                               ! This function = zPrime [forward] or z [backward] (meters)
   REAL(SiKi),     INTENT(IN )    :: Zeta                                            ! Instantaneous elevation of incident waves (meters)
   REAL(SiKi),     INTENT(IN )    :: zOrzPrime                                       ! Elevations where the wave kinematics are to be applied using Wheeler stretching, z, [forward] or elevations where the wave kinematics are computed before applying Wheeler stretching, zPrime, [backward] (meters)
   CHARACTER(1),   INTENT(IN )    :: ForwardOrBackWard                               ! A string holding the direction ('F'=Forward, 'B'=Backward) for applying Wheeler stretching.
   INTEGER(IntKi), INTENT(OUT)    :: ErrStat                                         ! Error status of the operation
   CHARACTER(*),   INTENT(OUT)    :: ErrMsg                                        ! Error message if ErrStat /= ErrID_None


      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ""

      ! Apply Wheeler stretching, depending on the direction:

   SELECT CASE ( ForwardOrBackWard )

   CASE ( 'F'  )  ! Forward

      WheelerStretching = ( 1.0 + Zeta/REAL(h,SiKi) )*zOrzPrime + Zeta


   CASE ( 'B' )   ! Backward

      WheelerStretching = ( zOrzPrime - Zeta )/( 1.0 + Zeta/REAL(h,SiKi) )


   CASE DEFAULT

      WheelerStretching = 0.0_SiKi

      ErrMsg = 'The last argument in routine WheelerStretching() must be ''F'' or ''B''.'
      ErrStat = ErrID_Fatal
      RETURN


   END SELECT



   RETURN
END FUNCTION WheelerStretching

!------------------------------------------------------------------------------------------------------------------------
SUBROUTINE CalculateWaveNDir(InitInp, InitOut, ErrStat, ErrMsg)
   TYPE(Waves_InitInputType),       INTENT(IN   )  :: InitInp     ! Input data for initialization routine
   TYPE(Waves_InitOutputType),      INTENT(INOUT)  :: InitOut     ! Output data
   INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   INTEGER(IntKi)                                  :: I           ! loop counter
   INTEGER(IntKi)                                  :: WaveNDirMax !< Maximum value we can change WaveNDir to (relative to original value passed in). Used in finding new WaveNDir value.
   
   INTEGER(IntKi)                                  :: ErrStatTmp  !< Temporary error status
   CHARACTER(ErrMsgLen)                            :: ErrMsgTmp   !< Temporary error message
   character(*), parameter                         :: RoutineName = 'CalculateWaveNDir'

   
         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""
   

      WaveNDirMax      = CEILING(InitInp%WaveNDir*1.25_SiKi)     ! Value we allow WaveNDir to reach before aborting
      InitOut%WaveNDir = InitInp%WaveNDir

         ! Check that the number of wave directions is a positive odd number.  In theory this has been
         ! done before the Waves module was called.  We repeat it here in the event that the Waves module
         ! gets used in some other code.
         !  -> If it is less than 0, error out.
         !  -> If it is even, we will increment it by 1.
      IF ( InitOut%WaveNDir <= 0_IntKi ) THEN
         CALL SetErrStat(ErrID_Fatal,'WaveNDir must be an odd number greater than 0.',ErrStat,ErrMsg,RoutineName)
         RETURN
      END IF

      IF ( MODULO( InitOut%WaveNDir, 2_IntKi) == 0_IntKi ) THEN
         InitOut%WaveNDir  = InitOut%WaveNDir + 1
         CALL SetErrStat(ErrID_Warn,'WaveNDir must be odd.  Changing the value to '//TRIM(Num2LStr(InitOut%WaveNDir)),ErrStat,ErrMsg,RoutineName)
      END IF

         
      ErrStatTmp = ErrID_None
      ErrMsgTmp = ""

      DO WHILE ( .NOT. EqualRealNos( REAL(InitOut%NStepWave2/InitOut%WaveNDir), REAL(InitOut%NStepWave2)/REAL(InitOut%WaveNDir) ))
                                  
         IF (InitOut%WaveNDir > WaveNDirMax ) THEN
            ErrMsgTmp   = 'Could not find value for WaveNDir between '//TRIM(Num2LStr(InitInp%WaveNDir))//' and '// &
                           TRIM(Num2LStr(WaveNDirMax))//' such that an equal number of frequencies are assigned to each direction.'
            ErrStatTmp  = ErrID_Fatal
            EXIT
         ELSE
            InitOut%WaveNDir = InitOut%WaveNDir + 2
            ErrMsgTmp   =  'Changed WaveNDir from '//TRIM(Num2LStr(InitInp%WaveNDir))//' to '// TRIM(Num2LStr(InitOut%WaveNDir))//  &
                           ' so that an equal number of frequencies are assigned to each direction.'
            ErrStatTmp  = ErrID_Warn
         END IF
               
      END DO
         
      CALL SetErrStat(ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName)


      IF (ErrStat == ErrID_Fatal) THEN

         ! If we exited because we hit a limit (in which case the condition is not satisfied), then we cannot continue.
         ! We warn the user that a value for WaveNDir was not found, and that they should try a different value, or try
         ! a different value for WaveTMax.  The reason for suggesting the latter is that NStepWave is derived from
         ! WaveTMax and adjusted until it is a product of smallish numbers (most likely even, but not necessarily so).
         ! So, there is a very small possibility then that NStepWave2 is a prime number, in which case we won't find a
         ! value for WaveNDir, so we suggest that the user change WaveTMax.  To make this a little easier for the user,
         ! we will report the first 5 possible values for WaveNDir between their requested value and 1/4 of NStepWave2,
         ! if there are any.
         
         
         ! Now check for the possible values of WaveNDir (up to I=5) so that we can tell the user about it.
         I = 0
         ErrMsgTmp = 'The next values of WaveNDir that work with the selected values for WaveTMax and WaveDT:'
         DO WHILE ( InitOut%WaveNDir <= INT(InitOut%NStepWave2/4.0) )
            IF ( EqualRealNos(REAL(InitOut%NStepWave2/InitOut%WaveNDir), &
                              REAL(InitOut%NStepWave2)/REAL(InitOut%WaveNDir) )) THEN
               ErrMsgTmp  = TRIM(ErrMsgTmp)//"  "//TRIM(Num2LStr(InitOut%WaveNDir))
               I = I + 1
               IF (I >= 5) EXIT ! limit the number of choices for WaveNDir that are printed
            END IF

            InitOut%WaveNDir = InitOut%WaveNDir + 2
         END DO

         ! If there were no additional values for WaveNDir found, I will be 0, so we rewrite the error message.
         IF ( I == 0 ) THEN
            ErrMsgTmp  =  'There are no values for WaveNDir between '//TRIM(Num2LStr(WaveNDirMax))//' and '// &
                           TRIM(Num2LStr(INT(InitOut%NStepWave2/4.0)))//' (4 frequencies per wave direction)'// &
                           ' that will work with the selected values for WaveTMax ('//TRIM(Num2Lstr(InitOut%WaveTMax))// &
                           ') and WaveDT ('//TRIM(Num2LStr(InitInp%WaveDT))//').  Change either WaveTMax or WaveDT.'
         ELSE
            ErrMsgTmp  = TRIM(ErrMsgTmp)//'.'
         ENDIF

         ! Append the message about the possible values for WaveNDir (if any were found) and set the error status before
         ! returning to the calling program.
         CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
         RETURN
      END IF

END SUBROUTINE CalculateWaveNDir
!------------------------------------------------------------------------------------------------------------------------
SUBROUTINE CalculateWaveDirection(InitInp, InitOut, ErrStat, ErrMsg )
! Compute the wave direction array, InitOut%WaveDirArr
!----------------------------------------------------------------------------------------------------------------------------------

   TYPE(Waves_InitInputType),       INTENT(IN   )  :: InitInp     ! Input data for initialization routine
   TYPE(Waves_InitOutputType),      INTENT(INOUT)  :: InitOut     ! Output data
   INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   ! Local Variables
   REAL(SiKi),                       ALLOCATABLE   :: WvTheta(:)                                      !< Final set of wave directions (degrees)
   REAL(SiKi),                       ALLOCATABLE   :: WvSpreadThetaIdx(:)                             !< Indices for wave directions
   INTEGER(IntKi)                                  :: WvSpreadFreqPerDir                              !< Number of wave frequencies per direction
   INTEGER                                         :: I                                               ! Generic index
   INTEGER                                         :: J                                               ! Generic index
   INTEGER                                         :: K                                               ! Generic index
   INTEGER                                         :: LastInd                                         ! Index into the arrays saved from the last call as a starting point for this call

   ! Variables for error handling
   INTEGER(IntKi)                                  :: ErrStatTmp                                      !< Temporary error status
   CHARACTER(ErrMsgLen)                            :: ErrMsgTmp                                       !< Temporary error message
   character(*), parameter                         :: RoutineName = 'CalculateWaveDirection'

   
      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ""


      IF (InitInp%WaveMod == 7) THEN ! wavemod 0 and 6 aren't called from this routine, but they fall into this case, too

         RETURN
         !InitOut%WaveDirArr set in UserWaveComponents_Init for WaveMod 7
         !InitOut%WaveDirArr = 0, set in Initial_InitOut_Arrays for WaveMod 0 and 6

      ELSEIF(.not. InitInp%WaveMultiDir .or. InitInp%WaveNDir <= 1) THEN ! we have a single wave direction
      
         InitOut%WaveDirArr = InitInp%WaveDir

      ELSE ! multi directional waves
         
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
         !!          was relaxed by setting the direction for OmegaArr(I) = 0 (which has no amplitude) since it was found that
         !!          (NStepWave2 - 1) is often a prime number due to how NStepWave is calculated above to be a product
         !!          of smallish numbers.

            ! this sets InitOut%WaveNDir:
         call CalculateWaveNDir(InitInp, InitOut, ErrStatTmp, ErrMsgTmp)
            call SetErrStat(ErrStatTmp,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
            if (ErrStat >= AbortErrLev) then
               call Cleanup()
               return
            end if
            

            ! This allocates and sets WvTheta:
         call CalculateWaveSpreading(InitInp, InitOut, WvTheta, ErrStatTmp, ErrMsgTmp)
            call SetErrStat(ErrStatTmp,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
            if (ErrStat >= AbortErrLev) then
               call Cleanup()
               return
            end if


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


            ! Allocate the index array for each group of frequencies.  This array is used to randomize the directions
            ! within each WaveNDir sized group of frequencies.  This is a REAL array used to hold the random numbers.
         ALLOCATE( WvSpreadThetaIdx(1:InitOut%WaveNDir), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) THEN
            CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WvSpreadThetaIdx while assigning wave directions.',ErrStat,ErrMsg,RoutineName)
            CALL CleanUp()
            RETURN
         END IF


            ! K should be exactly NStepWave2 when done assigning directions. The the Omega = 0 has
            ! no amplitude, but gets a direction anyhow (to simplify the calculation of WaveNDir).
         WvSpreadFreqPerDir   =  (InitOut%NStepWave2)/InitOut%WaveNDir
         K  = 0
            ! Work through the frequencies in groups of directions.
         DO I = 1,WvSpreadFreqPerDir

               ! Populate the array with random numbers
            CALL UniformRandomNumbers(InitInp%RNG%pRNG, WvSpreadThetaIdx)

            DO J = 1, InitOut%WaveNDir

                  ! Find the index lowest value in the WvSpreadThetaIdx array.  This is the index to
                  ! use for this wave direction.
               LastInd  = MINLOC( WvSpreadThetaIdx, DIM=1 )

                  ! Assign the direction for this frequency piece to the LastInd value.
               InitOut%WaveDirArr(K)   =  WvTheta( LastInd )

                  ! Now make that element in the WvSpreadThetaIdx really big so we don't pick it again
               WvSpreadThetaIdx( LastInd )   = HUGE(1.0_SiKi)

               K  = K + 1     ! Increment the frequency index

            ENDDO
         ENDDO
         
         ! Filling last value since it is not reached by the loop above
         CALL UniformRandomNumbers(InitInp%RNG%pRNG, WvSpreadThetaIdx)
         LastInd  = MINLOC( WvSpreadThetaIdx, DIM=1 )
         InitOut%WaveDirArr(K)   =  WvTheta( LastInd )

            ! Perform a quick sanity check.  We should have assigned all wave frequencies a direction, so K should be
            ! K = NStepWave2 (K is incrimented afterwards).
         IF ( K /= (InitOut%NStepWave2 ) ) THEN
            CALL SetErrStat(ErrID_Fatal, 'Something went wrong while assigning wave directions.',ErrStat,ErrMsg,RoutineName)
            CALL CleanUp()
            RETURN
         END IF

      ENDIF    ! Multi-directional waves in use.
      
      
      CALL Cleanup()
      
CONTAINS
   SUBROUTINE Cleanup()
   
      IF(ALLOCATED( WvSpreadThetaIdx ))   DEALLOCATE( WvSpreadThetaIdx )
      IF(ALLOCATED( WvTheta ))            DEALLOCATE( WvTheta )
   
   END SUBROUTINE Cleanup

END SUBROUTINE CalculateWaveDirection
!------------------------------------------------------------------------------------------------------------------------
SUBROUTINE CalculateWaveSpreading(InitInp, InitOut, WvTheta, ErrStat, ErrMsg )
! Compute the wave direction array
!----------------------------------------------------------------------------------------------------------------------------------

   TYPE(Waves_InitInputType),       INTENT(IN   )  :: InitInp     !< Input data for initialization routine
   TYPE(Waves_InitOutputType),      INTENT(INOUT)  :: InitOut     !< Output data
   REAL(SiKi),         ALLOCATABLE, INTENT(  OUT)  :: WvTheta(:)  !< Final set of wave directions (degrees)
   INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

   INTEGER(IntKi)                                  :: i                                               !< loop counter
   INTEGER(IntKi)                                  :: LastInd                                         !< last index found (for faster finding of next index)
   INTEGER(IntKi)                                  :: WvSpreadNDir                                    !< Number of wave spreading directions for intermediate calculations.  Set later to be MAX(15*InitOut%WaveNDir,1000)
   REAL(SiKi),                       ALLOCATABLE   :: WvSpreadCos2SArr(:)                             !< Wave spreading function results array.  Used in equal energy wave spreading function.
   REAL(SiKi)                                      :: WvSpreadCos2SConst                              !< Normalization constant for wave spreading function.
   REAL(SiKi),                       ALLOCATABLE   :: WvSpreadIntegral(:)                             !< Cumulative integral of the wave spreading function.  Used in finding equal energy wave directions.
   REAL(SiKi)                                      :: WvSpreadDTheta                                  !< Wave direction step size for intermediate calculations.  Used in finding equal energy wave directions.
   REAL(SiKi),                       ALLOCATABLE   :: WvSpreadThetas(:)                               !< Wave direction used in calculations and interpolations
   REAL(SiKi)                                      :: WvSpreadIntegralTmp                             !< Temporary variable for the interpolation
   
   ! Variables for error handling
   INTEGER(IntKi)                                  :: ErrStatTmp                                      !< Temporary error status
   character(*), parameter                         :: RoutineName = 'CalculateWaveSpreading'
   
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
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

         !> ### Code Implementation order
         !! 1. Discretize the spreading function range and calculate the values of the wave spreading function

            ! Now that we have the value for _WaveNDir_ found above, we set the value of _WvSpreadNDir_ to be 15x as
            ! large, or 1000 (whichever is larger).  WvSpreadNDir is used only in discretization for later
            ! interpolation of actual wave directions.
         WvSpreadNDir   = MAX(15*InitOut%WaveNDir,1000)
         WvSpreadDTheta = InitInp%WaveDirRange/REAL(WvSpreadNDir,SiKi)

            ! Calculate the normalization constant for the wave spreading.
         IF ( InitInp%WaveDirSpread < 25.0_SiKi ) THEN ! Use exact expression
            WvSpreadCos2SConst   = sqrt(Pi) * (NWTC_GAMMA(InitInp%WaveDirSpread + 1.0_SiKi))/(InitInp%WaveDirRange * NWTC_GAMMA(InitInp%WaveDirSpread + 0.5_SiKi))
         ELSE ! Use asymptotic approximation for large argument
            WvSpreadCos2SConst   = sqrt(Pi*InitInp%WaveDirSpread)*(1.0_SiKi+0.125_SiKi/InitInp%WaveDirSpread)/InitInp%WaveDirRange
         ENDIF

            ! Allocate arrays to use for storing the intermediate values
         ALLOCATE( WvSpreadCos2SArr(0:WvSpreadNDir),  STAT=ErrStatTmp ); IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WvSpreadCos2SArr.',  ErrStat,ErrMsg,RoutineName)
         ALLOCATE( WvSpreadIntegral(0:WvSpreadNDir),  STAT=ErrStatTmp ); IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WvSpreadIntegral.',  ErrStat,ErrMsg,RoutineName)
         ALLOCATE( WvSpreadThetas(  0:WvSpreadNDir),  STAT=ErrStatTmp ); IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WvSpreadThetas.',    ErrStat,ErrMsg,RoutineName)
         ALLOCATE( WvTheta(1:InitOut%WaveNDir),       STAT=ErrStatTmp ); IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WvTheta.',           ErrStat,ErrMsg,RoutineName)

         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF

            !> 2. Calculate the spreading function as a function of angle.  Step through all _WvSpreadNDir_ steps.
         DO I=0,WvSpreadNDir
               ! The current angle as we step through the range
            WvSpreadThetas(I) =  I*WvSpreadDTheta  + InitInp%WaveDir - InitInp%WaveDirRange/(2.0_SiKi)

               ! Calculate the wave spreading for the current value of WvSpreadThetas
            WvSpreadCos2SArr(I)  =  WvSpreadCos2SConst*abs( cos(Pi*(WvSpreadThetas(I)-InitInp%WaveDir)/InitInp%WaveDirRange) ) **(2*InitInp%WaveDirSpread)

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
         IF ( WvSpreadIntegral(WvSpreadNDir) < 0.99_SiKi .OR. WvSpreadIntegral(WvSpreadNDir) > 1.01_SiKi ) THEN
            CALL SetErrStat(ErrID_Fatal,' Something went wrong in evaluating the multidirectional wave spreading function.  '// &
                           'Integral is '//TRIM(Num2LStr(WvSpreadIntegral(WvSpreadNDir))),ErrStat,ErrMsg,RoutineName)
            call Cleanup()
            RETURN
         END IF


            !> 5. Set the wave directions using the results from the integral.
            !  We will use the variable LastInd as a simple index for figuring out where in the array we are.  First set to 0
         LastInd  =  0_IntKi
         DO I=1,InitOut%WaveNDir
            WvSpreadIntegralTmp  = (REAL(I)-0.5_SiKi)/REAL(InitOut%WaveNDir)
            WvTheta(I)    =  InterpStp( WvSpreadIntegralTmp, WvSpreadIntegral, WvSpreadThetas, LastInd, WvSpreadNDir )
         ENDDO       ! I=1,InitOut%WaveNDir


            !> 6. Done with equal energy wavedirection calculations.  Deallocate the arrays used during calculations.

         CALL CleanUp()


contains
   subroutine Cleanup()
      IF(ALLOCATED( WvSpreadCos2SArr ))      DEALLOCATE( WvSpreadCos2SArr, STAT=ErrStatTmp )
      IF(ALLOCATED( WvSpreadIntegral ))      DEALLOCATE( WvSpreadIntegral, STAT=ErrStatTmp )
      IF(ALLOCATED( WvSpreadThetas   ))      DEALLOCATE( WvSpreadThetas,   STAT=ErrStatTmp )
   end subroutine Cleanup
   
END SUBROUTINE CalculateWaveSpreading
!------------------------------------------------------------------------------------------------------------------------
!> sets WaveS1SddArr(:) and InitOut%WaveElevC0
SUBROUTINE Get_1Spsd_and_WaveElevC0(InitInp, InitOut, OmegaArr, WaveS1SddArr)

   TYPE(Waves_InitInputType),       INTENT(IN   )  :: InitInp                                       ! Input data for initialization routine
   TYPE(Waves_InitOutputType),      INTENT(INOUT)  :: InitOut                                       ! Output data
   REAL(SiKi),                      INTENT(IN   )  :: OmegaArr(0:InitOut%NStepWave2)                !< Array of all non-negative angular frequencies (rad/s)
   REAL(SiKi),                      INTENT(  OUT)  :: WaveS1SddArr(0:InitOut%NStepWave2)            !< One-sided power spectral density of the wave spectrum at all non-negative frequencies (m^2/(rad/s))

   COMPLEX(SiKi)                                   :: WGNC(0:InitOut%NStepWave2)                    ! Discrete Fourier transform of the realization of a White Gaussian Noise (WGN) time series process with unit variance for the current frequency component (-)
   INTEGER                                         :: I                                             ! Loop counter
   INTEGER                                         :: I_WaveTp                                      ! The index of the frequency component nearest to WaveTp
   REAL(SiKi)                                      :: SQRTNStepWave2                                ! SQRT( NStepWave/2 )
   COMPLEX(SiKi)                                   :: tmpComplex                                    ! A temporary varible to hold the complex value of the wave elevation before storing it into a REAL array
   REAL(SiKi)                                      :: WaveS2Sdd                                     ! Two-sided power spectral density of the wave spectrum per unit time for the current frequency component (m^2/(rad/s))
   
   
      IF ( InitInp%WaveMod == 5 .OR. InitInp%WaveMod == 7) THEN    ! Wave elevation or frequency component data read in
   
         DO I = 0,InitOut%NStepWave2
         
            ! Apply limits to the existing WaveElevC0 arrays if outside frequency range
            IF ( OmegaArr(I) < InitInp%WvLowCOff .OR. OmegaArr(I) > InitInp%WvHiCOff )  THEN
               InitOut%WaveElevC0(:,I) = 0.0_SiKi
            ENDIF
            
         END DO
      
         WaveS1SddArr = 0 ! unused here
         RETURN
      
      END IF
   
   
      I_WaveTp  = NINT ( TwoPi/(InitOut%WaveDOmega*InitInp%WaveTp) )        ! Compute the index of the frequency component nearest to WaveTp. Note, we don't check if it's a valid index into the arrays
   
      ! Compute the discrete Fourier transform of the realization of a White
      !   Gaussian Noise (WGN) time series process with unit variance:

      ! ---------------------------------
      ! Set White Gaussian Noise with unit variance
      !
      ! NOTE: For the time series process to be real with zero mean, the values at
      !       OmegaArr(I) == 0.0 and OmegaArr(I) == NStepWave2*WaveDOmega (= WaveOmegaMax)
      !       must be zero.
      !---------------------------------
      ! I == 1 or InitOut%NStepWave2 if ( OmegaArr(I) == 0.0 ) or ( OmegaArr(I) == NStepWave2*WaveDOmega (= WaveOmegaMax) )
      WGNC(1)                  = (0.0,0.0)
      WGNC(InitOut%NStepWave2) = (0.0,0.0)
      
      IF ( InitInp%WaveMod == 10 )  THEN                     ! .TRUE. for plane progressive (regular) waves with a specified phase
         DO I = 0,InitOut%NStepWave2-1                       ! Loop through the positive frequency components (including zero) of the discrete Fourier transforms
            IF (I==1) CYCLE
            
            WGNC(I) = BoxMuller ( InitInp%RNG%pRNG, InitInp%WaveNDAmp, InitInp%WavePhase )
         END DO
      ELSE                                               ! All other OmegaArr(I)
         DO I = 0,InitOut%NStepWave2-1  ! Loop through the positive frequency components (including zero) of the discrete Fourier transforms
            IF (I==1) CYCLE
            
            WGNC(I) = BoxMuller ( InitInp%RNG%pRNG, InitInp%WaveNDAmp )
         END DO
      END IF
      
      !------------------------------------
      ! For (WaveMod=1 plane progressive (regular); and WaveMod=10 plane progressive (regular) waves with a specified phase)
      ! adjust WGNC and set PSD at specified frequency
      !------------------------------------
      IF (InitInp%WaveMod == 10 .or. InitInp%WaveMod == 1) THEN
         WaveS1SddArr = 0.0
         
         IF (I_WaveTp < InitOut%NStepWave2 .and. (I_WaveTp > 1 .or. I_WaveTp == 0) ) THEN
             
               ! This scaling of WGNC is used to ensure that the Box-Muller method is only providing a random phase,
               ! not a magnitude change, at the frequency of the plane progressive wave.  The SQRT(2.0) is used to
               ! ensure that the time series WGN process has unit variance (i.e. sinusoidal with amplitude SQRT(2.0)).
               !
               ! NOTE: the denominator here will never equal zero since U1 cannot equal 1.0, and thus, C1 cannot be 0.0 in the Box-Muller method.
             
              WGNC(I_WaveTp)         = WGNC(I_WaveTp) * ( SQRT(2.0_SiKi) / ABS(WGNC(I_WaveTp)) )
              
               ! Plane progressive (regular) wave; the wave spectrum is an impulse function centered on frequency component closest to WaveTp.              
              WaveS1SddArr(I_WaveTp) = 0.5_SiKi * (InitInp%WaveHs/2.0_SiKi)**2 / InitOut%WaveDOmega
              
         END IF
      ELSE
      
         DO I = 0,InitOut%NStepWave2

            IF ( OmegaArr(I) < InitInp%WvLowCOff .OR. OmegaArr(I) > InitInp%WvHiCOff )  THEN ! .TRUE. if OmegaArr(I) is above or below the cut-off frequency
               !  Zero-out the wave spectrum above the cut-off frequency.  We must cut-off the frequency in order to
               !  void nonphysical wave forces.  Waves that have wavelengths much smaller than the platform diameter
               !  (high frequency) do not contribute to the net force because regions of positive and negative
               !  velocity/acceleration are experienced by the platform at the same time and cancel out.
            
               WaveS1SddArr(I) = 0.0
               
            ELSE
            
               SELECT CASE ( InitInp%WaveMod ) ! Which incident wave kinematics model are we using?
                  CASE ( 2 )              ! JONSWAP/Pierson-Moskowitz spectrum (irregular) wave.
                        WaveS1SddArr(I) = JONSWAP ( OmegaArr(I), InitInp%WaveHs, InitInp%WaveTp, InitInp%WavePkShp )
                  CASE ( 3 )              ! White-noise
                        WaveS1SddArr(I) =  InitInp%WaveHs * InitInp%WaveHs / ( 16.0 * (InitInp%WvHiCOff - InitInp%WvLowCOff) )
                  CASE ( 4 )              ! User-defined spectrum (irregular) wave.
                        CALL UserWaveSpctrm ( OmegaArr(I), InitInp%WaveDir, InitInp%DirRoot, WaveS1SddArr(I) )
               ENDSELECT
         
            END IF
            
         END DO
      
      
      END IF

      
      ! ---------------------------------
      ! Compute the one-sided power spectral density of the wave spectrum per unit
      !   time; zero-out the wave spectrum above the cut-off frequency:
      !---------------------------------
      SQRTNStepWave2 = SQRT( REAL( InitOut%NStepWave2, SiKi ) )                  ! Compute SQRT( NStepWave/2 ).
      
      DO I = 0,InitOut%NStepWave2
            ! Compute the two-sided power spectral density of the wave spectrum per unit
            !   time:

         WaveS2Sdd = 0.5_SiKi*WaveS1SddArr(I)

            ! Compute the discrete Fourier transform of the instantaneous elevation of
            !   incident waves at the WAMIT reference point:
         tmpComplex                   = SQRTNStepWave2 * WGNC(I) *SQRT( TwoPi_R4 * WaveS2Sdd / REAL(InitInp%WaveDT,SiKi) )
         InitOut%WaveElevC0     (1,I) = REAL( tmpComplex)
         InitOut%WaveElevC0     (2,I) = AIMAG(tmpComplex)

      END DO   ! I - The positive frequency components (including zero) of the discrete Fourier transforms
      
END SUBROUTINE Get_1Spsd_and_WaveElevC0
!------------------------------------------------------------------------------------------------------------------------
!> update InitOut%WaveElevC0; call InitFFT before calling this routine!
SUBROUTINE ConstrainedNewWaves(InitInp, InitOut, OmegaArr, WaveS1SddArr, CosWaveDir, SinWaveDir, FFT_Data, ErrStat, ErrMsg)

   TYPE(Waves_InitInputType),       INTENT(IN   )  :: InitInp                                       ! Input data for initialization routine
   TYPE(Waves_InitOutputType),      INTENT(INOUT)  :: InitOut                                       ! Output data
   REAL(SiKi),                      INTENT(IN   )  :: OmegaArr(0:InitOut%NStepWave2)                !< Array of all non-negative angular frequencies (rad/s)
   REAL(SiKi),                      INTENT(IN   )  :: WaveS1SddArr(0:InitOut%NStepWave2)            !< One-sided power spectral density of the wave spectrum at all non-negative frequencies (m^2/(rad/s))
   REAL(SiKi),                      INTENT(IN   )  :: CosWaveDir(0:InitOut%NStepWave2)              !< COS( WaveDirArr(I) ) -- Each wave frequency has a unique wave direction
   REAL(SiKi),                      INTENT(IN   )  :: SinWaveDir(0:InitOut%NStepWave2)              !< SIN( WaveDirArr(I) ) -- Each wave frequency has a unique wave direction
   TYPE(FFT_DataType),              INTENT(IN   )  :: FFT_Data                                      !< data for FFT computations, already initialized
   INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat                                       !< error level/status
   CHARACTER(ErrMsgLen),            INTENT(  OUT)  :: ErrMsg                                        !< error message

   
   REAL(SiKi)                                      :: WaveNmbr                                      ! Wavenumber of the current frequency component (1/meter)
   INTEGER                                         :: I                                             ! Generic index
   
   ! Variables for constrained wave
   REAL(SiKi)                                      :: WaveElevC0ReSum                               !< Sum of the wave DFT amplitudes (real part) across all frequencies (m)
   REAL(SiKi)                                      :: WaveElevC0ImOmegaSum                          !< Sum of the wave DFT amplitudes (imaginary part) times the angular frequency across all frequencies (m(rad/s))
   REAL(SiKi)                                      :: Crest                                         !< Crest elevation measured from SWL (m)
   REAL(SiKi)                                      :: CrestHeight                                   !< Crest height measured from the crest to the preceding or following trough (m)
   REAL(SiKi)                                      :: CrestHeight1                                  !< Crest height with purturbed crest elevation (m)
   REAL(SiKi)                                      :: CrestHeightError                              !< Error in crest height relative to the specified crest height (m)
   REAL(SiKi)                                      :: ConstWavePhase                                !< Phase adjustment to wave DFT amplitudes due to constrained wave (m)
   REAL(SiKi)                                      :: Trough                                        !< The trough preceding or following the crest, whichever is lower (m)
   REAL(SiKi)                                      :: m0                                            !< Zeroth spectral moment of the wave spectrum (m^2)
   REAL(SiKi)                                      :: m2                                            !< First spectral moment of the wave spectrum (m^2(rad/s)^2)
   REAL(SiKi)                                      :: CrestHeightTol = 1.0E-3                       !< Relative tolerance for the crest height when ConstWaveMod = 2
   INTEGER(IntKi)                                  :: NStepTp                                       !< Number of time steps per peak period when waveMod = 2 (-)
   INTEGER(IntKi)                                  :: Iter                                          !< Number of iterations when trying to meet the prescribed crest height (-)
   INTEGER(IntKi)                                  :: MaxCrestIter = 20                             !< Maximum number of iterations when trying to meet the prescribed crest height (-)
   
   REAL(SiKi)                                      :: tmpArr(0:InitOut%NStepWave2)                  !< A temporary array of real numbers of constrained wave (-)
   COMPLEX(SiKi)                                   :: tmpComplexArr(0:InitOut%NStepWave2)           !< A temporary array for FFT use
   
   COMPLEX(SiKi)                                   :: tmpComplex                                    ! A temporary varible to hold the complex value of the wave elevation before storing it into a REAL array
   
   INTEGER(IntKi)                                  :: ErrStatTmp                                    !< error level/status
!   CHARACTER(ErrMsgLen)                            :: ErrMsgTmp                                     !< error message
   CHARACTER(*), PARAMETER                         :: RoutineName = 'ConstrainedNewWaves'
   
   
   ErrStat = ErrID_None
   ErrMsg = ""
   
      !=== Constrained New Waves ===
      ! Modify the wave components to implement the constrained wave
   
      ! Compute the relevant sums
      m0                   = InitOut%WaveDOmega * SUM(WaveS1SddArr)
      m2                   = InitOut%WaveDOmega * SUM(WaveS1SddArr*OmegaArr*OmegaArr)
      WaveElevC0ReSum      = SUM(InitOut%WaveElevC0(1,:))/m0
      WaveElevC0ImOmegaSum = SUM(InitOut%WaveElevC0(2,:) * OmegaArr)/m2
      ! Apply the part of the modification that is independent from the crest elevation
      InitOut%WaveElevC0(1,:) = InitOut%WaveElevC0(1,:) - WaveElevC0ReSum                 * WaveS1SddArr * InitOut%WaveDOmega
      InitOut%WaveElevC0(2,:) = InitOut%WaveElevC0(2,:) - WaveElevC0ImOmegaSum * OmegaArr * WaveS1SddArr * InitOut%WaveDOmega

      Crest = 0.5_SiKi * InitInp%CrestHmax ! Set crest elevation to half of crest height
      tmpArr = InitOut%NStepWave2/m0 * InitOut%WaveDOmega * WaveS1SddArr
        
      IF (InitInp%ConstWaveMod == 1) THEN  ! Crest elevation prescribed
      
         ! Apply the remaining part of the modification proportional to crest elevation
         InitOut%WaveElevC0(1,:) = InitOut%WaveElevC0(1,:) + Crest * tmpArr
         
      ELSE IF (InitInp%ConstWaveMod == 2) THEN ! Crest height prescribed - Need to interate
      
         NStepTp = CEILING(InitInp%WaveTp/InitInp%WaveDT)

         Iter = 0
         CrestHeightError = InitInp%CrestHmax
         DO WHILE(CrestHeightError>CrestHeightTol .AND. Iter<=MaxCrestIter)
            Iter = Iter + 1

            ! Compute the crest height based on the current guess of crest elevation
            tmpComplexArr = CMPLX(  InitOut%WaveElevC0(1,:) + Crest * tmpArr, &
                                    InitOut%WaveElevC0(2,:))
            CALL ApplyFFT_cx (  InitOut%WaveElev0    (0:InitOut%NStepWave-1),  tmpComplexArr    (:  ), FFT_Data, ErrStatTmp )
            CALL SetErrStat(ErrStatTmp,'Error occured while applying the FFT to WaveElev0.',ErrStat,ErrMsg,RoutineName)
            IF ( ErrStat >= AbortErrLev ) RETURN

            ! Find the preceding or following trough, whichever is lower
            Trough = MIN(MINVAL(InitOut%WaveElev0(1:MIN(NStepTp,InitOut%NStepWave-1))), &
                         MINVAL(InitOut%WaveElev0(MAX(InitOut%NStepWave-NStepTp,0):InitOut%NStepWave-1)))
            CrestHeight = Crest-Trough
            CrestHeightError = ABS(CrestHeight - InitInp%CrestHmax)
            ! print *, CrestHeight

            If (CrestHeightError>CrestHeightTol) THEN ! If crest height tolerance is not satisfied
               ! Compute the crest height based on a slightly nudged crest elevation
               tmpComplexArr = CMPLX(  InitOut%WaveElevC0(1,:) + (Crest+CrestHeightTol) * tmpArr, &
                                       InitOut%WaveElevC0(2,:))
               CALL ApplyFFT_cx (  InitOut%WaveElev0    (0:InitOut%NStepWave-1),  tmpComplexArr    (:  ), FFT_Data, ErrStatTmp )
               CALL SetErrStat(ErrStatTmp,'Error occured while applying the FFT to WaveElev0.',ErrStat,ErrMsg,RoutineName)
               IF ( ErrStat >= AbortErrLev ) RETURN

               
               ! Find the preceding or following trough, whichever is lower
               Trough = MIN(MINVAL(InitOut%WaveElev0(1:MIN(NStepTp,InitOut%NStepWave-1))), &
                           MINVAL(InitOut%WaveElev0(MAX(InitOut%NStepWave-NStepTp,0):InitOut%NStepWave-1)))
               CrestHeight1 = Crest+CrestHeightTol-Trough
               ! Update crest elevation with Newton-Raphson Method
               Crest = Crest - (CrestHeight-InitInp%CrestHmax)*CrestHeightTol/(CrestHeight1-CrestHeight)
            ENDIF
         END DO
         
         ! Apply the remaining part of the modification based on the final crest elevation
         InitOut%WaveElevC0(1,:) = InitOut%WaveElevC0(1,:) + Crest * tmpArr
      ENDIF
      
      ! Modify the wave phase so that the crest shows up at the right place and the right time
      DO I = 1,InitOut%NStepWave2-1
         WaveNmbr   = WaveNumber ( OmegaArr(I), InitInp%Gravity, InitInp%WtrDpth )
         ConstWavePhase = WaveNmbr*(CosWaveDir(I)*InitInp%CrestXi  + &
                                    SinWaveDir(I)*InitInp%CrestYi) - &
                                    OmegaArr(I)*InitInp%CrestTime
         tmpComplex = CMPLX( InitOut%WaveElevC0(1,I) , InitOut%WaveElevC0(2,I)  )
         tmpComplex = tmpComplex * CMPLX( cos(ConstWavePhase), sin(ConstWavePhase)  )
         InitOut%WaveElevC0(1,I) = REAL(tmpComplex)
         InitOut%WaveElevC0(2,I) = AIMAG(tmpComplex)
      END DO

END SUBROUTINE ConstrainedNewWaves
!------------------------------------------------------------------------------------------------------------------------
END MODULE Waves
!**********************************************************************************************************************************
