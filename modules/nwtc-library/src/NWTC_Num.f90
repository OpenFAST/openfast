!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2013-2016  National Renewable Energy Laboratory
!
!    This file is part of the NWTC Subroutine Library.
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

!..................................................................................................................................
!> This module contains numeric-type routines with non-system-specific logic and references.
MODULE NWTC_Num
!..................................................................................................................................
   
   USE                                          NWTC_IO

   IMPLICIT NONE

!=======================================================================

      ! Global numeric-related variables.

   REAL(DbKi)                                :: D2R_D                         !< Factor to convert degrees to radians in double precision
   REAL(DbKi)                                :: Inf_D                         !< IEEE value for NaN (not-a-number) in double precision
   REAL(DbKi)                                :: Inv2Pi_D                      !< 0.5/Pi (1/(2*Pi)) in double precision
   REAL(DbKi)                                :: NaN_D                         !< IEEE value for Inf (infinity) in double precision
   REAL(DbKi)                                :: Pi_D                          !< Ratio of a circle's circumference to its diameter in double precision
   REAL(DbKi)                                :: PiBy2_D                       !< Pi/2 in double precision
   REAL(DbKi)                                :: R2D_D                         !< Factor to convert radians to degrees in double precision
   REAL(DbKi)                                :: RPM2RPS_D                     !< Factor to convert revolutions per minute to radians per second in double precision
   REAL(DbKi)                                :: RPS2RPM_D                     !< Factor to convert radians per second to revolutions per minute in double precision
   REAL(DbKi)                                :: TwoByPi_D                     !< 2/Pi in double precision
   REAL(DbKi)                                :: TwoPi_D                       !< 2*Pi in double precision


   REAL(ReKi)                                :: D2R                           !< Factor to convert degrees to radians
   REAL(ReKi)                                :: Inf                           !< IEEE value for NaN (not-a-number)
   REAL(ReKi)                                :: Inv2Pi                        !< 0.5/Pi = 1 / (2*pi)
   REAL(ReKi)                                :: NaN                           !< IEEE value for Inf (infinity)
   REAL(ReKi)                                :: Pi                            !< Ratio of a circle's circumference to its diameter
   REAL(ReKi)                                :: PiBy2                         !< Pi/2
   REAL(ReKi)                                :: R2D                           !< Factor to convert radians to degrees
   REAL(ReKi)                                :: RPM2RPS                       !< Factor to convert revolutions per minute to radians per second
   REAL(ReKi)                                :: RPS2RPM                       !< Factor to convert radians per second to revolutions per minute
   REAL(ReKi)                                :: TwoByPi                       !< 2/Pi
   REAL(ReKi)                                :: TwoPi                         !< 2*Pi

   REAL(SiKi)                                :: D2R_S                         !< Factor to convert degrees to radians in single precision
   REAL(SiKi)                                :: Inf_S                         !< IEEE value for NaN (not-a-number) in single precision
   REAL(SiKi)                                :: Inv2Pi_S                      !< 0.5/Pi (1/(2*Pi)) in single precision
   REAL(SiKi)                                :: NaN_S                         !< IEEE value for Inf (infinity) in single precision
   REAL(SiKi)                                :: Pi_S                          !< Ratio of a circle's circumference to its diameter in single precision
   REAL(SiKi)                                :: PiBy2_S                       !< Pi/2 in single precision
   REAL(SiKi)                                :: R2D_S                         !< Factor to convert radians to degrees in single precision
   REAL(SiKi)                                :: RPM2RPS_S                     !< Factor to convert revolutions per minute to radians per second in single precision
   REAL(SiKi)                                :: RPS2RPM_S                     !< Factor to convert radians per second to revolutions per minute in single precision
   REAL(SiKi)                                :: TwoByPi_S                     !< 2/Pi in single precision
   REAL(SiKi)                                :: TwoPi_S                       !< 2*Pi in single precision

   REAL(SiKi)                                :: Pi_R4                         !< Ratio of a circle's circumference to its diameter in 4-byte precision
   REAL(R8Ki)                                :: Pi_R8                         !< Ratio of a circle's circumference to its diameter in 8-byte precision

   REAL(SiKi)                                :: TwoPi_R4                      !< 2*pi in 4-byte precision
   REAL(R8Ki)                                :: TwoPi_R8                      !< 2*pi in 8-byte precision
   
   ! constants for kernel smoothing
   INTEGER, PARAMETER :: kernelType_EPANECHINIKOV = 1
   INTEGER, PARAMETER :: kernelType_QUARTIC       = 2
   INTEGER, PARAMETER :: kernelType_BIWEIGHT      = 3
   INTEGER, PARAMETER :: kernelType_TRIWEIGHT     = 4
   INTEGER, PARAMETER :: kernelType_TRICUBE       = 5
   INTEGER, PARAMETER :: kernelType_GAUSSIAN      = 6

   
      ! constants for output formats
   INTEGER, PARAMETER                        :: Output_in_Native_Units = 0
   INTEGER, PARAMETER                        :: Output_in_SI_Units     = 1
   INTEGER, PARAMETER                        :: Output_in_Engr_Units   = 2
!=======================================================================

      ! Create interfaces for generic routines that use specific routines.

      !> \copydoc nwtc_num::equalrealnos4()
   INTERFACE EqualRealNos
      MODULE PROCEDURE EqualRealNos4
      MODULE PROCEDURE EqualRealNos8
   END INTERFACE

      !> \copydoc nwtc_num::eulerconstructr4()
   INTERFACE EulerConstruct
      MODULE PROCEDURE EulerConstructR4
      MODULE PROCEDURE EulerConstructR8
   END INTERFACE

   INTERFACE EulerConstructZYX
      MODULE PROCEDURE EulerConstructZYXR8
   END INTERFACE
   
      !> \copydoc nwtc_num::eulerextractr4()
   INTERFACE EulerExtract
      MODULE PROCEDURE EulerExtractR4
      MODULE PROCEDURE EulerExtractR8
   END INTERFACE

      !> \copydoc nwtc_num::taitbryanyxzextractr4()
      !! See nwtc_num::taitbryanyxzextractr4() for details on the algorithm
   INTERFACE TaitBryanYXZExtract
      MODULE PROCEDURE TaitBryanYXZExtractR4
      MODULE PROCEDURE TaitBryanYXZExtractR8
   END INTERFACE
   
   INTERFACE TaitBryanYXZConstruct
      MODULE PROCEDURE TaitBryanYXZConstructR4
      MODULE PROCEDURE TaitBryanYXZConstructR8
   END INTERFACE

      !> \copydoc nwtc_num::outerproductr4
   INTERFACE OuterProduct
      MODULE PROCEDURE OuterProductR4
      MODULE PROCEDURE OuterProductR8
   END INTERFACE

      !> \copydoc nwtc_num::cross_productr4()
   INTERFACE Cross_Product
      MODULE PROCEDURE Cross_ProductR4
      MODULE PROCEDURE Cross_ProductR4R8
      MODULE PROCEDURE Cross_ProductR8
      MODULE PROCEDURE Cross_ProductR8R4
   END INTERFACE
   
      !> \copydoc nwtc_num::smllrottransd()
   INTERFACE SmllRotTrans
      MODULE PROCEDURE SmllRotTransD
      MODULE PROCEDURE SmllRotTransR
   END INTERFACE

      !> \copydoc nwtc_num::getsmllrotangsd()
   INTERFACE GetSmllRotAngs
      MODULE PROCEDURE GetSmllRotAngsD
      MODULE PROCEDURE GetSmllRotAngsR
   END INTERFACE
  
      !> \copydoc nwtc_num::zero2twopir4
   INTERFACE Zero2TwoPi
      MODULE PROCEDURE Zero2TwoPiR4
      MODULE PROCEDURE Zero2TwoPiR8
   END INTERFACE
   
      !> \copydoc nwtc_num::twonormr4
   INTERFACE TwoNorm
      MODULE PROCEDURE TwoNormR4
      MODULE PROCEDURE TwoNormR8
   END INTERFACE
   
      !> \copydoc nwtc_num::tracer4
   INTERFACE trace
      MODULE PROCEDURE traceR4
      MODULE PROCEDURE traceR8
   END INTERFACE
   
      !> \copydoc nwtc_num::dcm_expd
   INTERFACE DCM_exp  
      MODULE PROCEDURE DCM_expR
      MODULE PROCEDURE DCM_expD
   END INTERFACE
   
      !> \copydoc nwtc_num::dcm_logmapd
   INTERFACE DCM_logMap
      MODULE PROCEDURE DCM_logMapR
      MODULE PROCEDURE DCM_logMapD
   END INTERFACE

      !> \copydoc nwtc_num::dcm_setlogmapforinterpd
   INTERFACE DCM_SetLogMapForInterp
      MODULE PROCEDURE DCM_SetLogMapForInterpR
      MODULE PROCEDURE DCM_SetLogMapForInterpD
   END INTERFACE
   
      !> \copydoc nwtc_num::eye2
   INTERFACE Eye
      MODULE PROCEDURE Eye2   ! matrix of two dimensions
      MODULE PROCEDURE Eye2D  ! matrix of two dimensions (double precision)
      MODULE PROCEDURE Eye3   ! matrix of three dimensions
      MODULE PROCEDURE Eye3D  ! matrix of three dimensions
   END INTERFACE

      !> \copydoc nwtc_num::interpbincomp
   INTERFACE InterpBin
      MODULE PROCEDURE InterpBinComp
      MODULE PROCEDURE InterpBinReal
   END INTERFACE

      !> \copydoc nwtc_num::interpstpcomp4
   INTERFACE InterpStp
      MODULE PROCEDURE InterpStpComp4
      MODULE PROCEDURE InterpStpComp8
      MODULE PROCEDURE InterpStpReal4
      MODULE PROCEDURE InterpStpReal4_8
      MODULE PROCEDURE InterpStpReal8
   END INTERFACE

      !> \copydoc nwtc_num::interparrayr4
   INTERFACE InterpArray
      MODULE PROCEDURE InterpArrayR4
      MODULE PROCEDURE InterpArrayR8
   END INTERFACE

      !> \copydoc nwtc_num::interpwrappedstpreal4
   INTERFACE InterpWrappedStpReal
      MODULE PROCEDURE InterpWrappedStpReal4
      MODULE PROCEDURE InterpWrappedStpReal4_8
      MODULE PROCEDURE InterpWrappedStpReal8
   END INTERFACE
   
      !> \copydoc nwtc_num::locatestpr4
   INTERFACE LocateStp
      MODULE PROCEDURE LocateStpR4
      MODULE PROCEDURE LocateStpR8
   END INTERFACE

   !> \copydoc nwtc_num::skewsymmatr4
   INTERFACE SkewSymMat
      MODULE PROCEDURE SkewSymMatR4
      MODULE PROCEDURE SkewSymMatR8
   END INTERFACE
   
      !> \copydoc nwtc_num::angle_extrapinterp2_r4
   INTERFACE Angles_ExtrapInterp
      MODULE PROCEDURE Angles_ExtrapInterp1_R4
      MODULE PROCEDURE Angles_ExtrapInterp1_R8
      MODULE PROCEDURE Angles_ExtrapInterp1_R4R
      MODULE PROCEDURE Angles_ExtrapInterp1_R8R
      MODULE PROCEDURE Angles_ExtrapInterp2_R4
      MODULE PROCEDURE Angles_ExtrapInterp2_R8
      MODULE PROCEDURE Angles_ExtrapInterp2_R4R
      MODULE PROCEDURE Angles_ExtrapInterp2_R8R
   END INTERFACE

      !> \copydoc nwtc_num::addorsub2pi_r4
   INTERFACE AddOrSub2Pi
      MODULE PROCEDURE AddOrSub2Pi_R4
      MODULE PROCEDURE AddOrSub2Pi_R8
   END INTERFACE
   
      !> \copydoc nwtc_num::mpi2pi_r4
   INTERFACE MPi2Pi
      MODULE PROCEDURE MPi2Pi_R4
      MODULE PROCEDURE MPi2Pi_R8
   END INTERFACE
   
CONTAINS

!=======================================================================
!> This routine is used to convert NewAngle to an angle within Pi of
!!   OldAngle by adding or subtracting 2*Pi accordingly.  
!!   This routine is useful for converting
!!   angles returned from a call to the ATAN2() FUNCTION into angles that may
!!   exceed the -Pi to Pi limit of ATAN2().  For example, if the nacelle yaw
!!   angle was 179deg in the previous time step and the yaw angle increased
!!   by 2deg in the new time step, we want the new yaw angle returned from a
!!   call to the ATAN2() FUNCTION to be 181deg instead of -179deg.  This
!!   routine assumes that the angle change between calls is not more than
!!   Pi in absolute value.
!! Use AddOrSub2Pi (nwtc_num::addorsub2pi) instead of directly calling a specific routine in the generic interface.
   SUBROUTINE AddOrSub2Pi_R4 ( OldAngle, NewAngle )
        ! Argument declarations:

   REAL(SiKi), INTENT(IN   )    :: OldAngle                                     !< Angle from which NewAngle will be converted to within Pi of, rad.
   REAL(SiKi), INTENT(INOUT)    :: NewAngle                                     !< Angle to be converted to within 2*Pi of OldAngle, rad.


      ! Local declarations:

   REAL(SiKi)                   :: DelAngle                                     ! The difference between OldAngle and NewAngle, rad.



      ! Add or subtract 2*Pi in order to convert NewAngle two within Pi of OldAngle:

   
   DelAngle = OldAngle - NewAngle

   DO WHILE ( ABS( DelAngle ) > Pi_R4 )

      NewAngle = NewAngle + SIGN( TwoPi_R4, DelAngle )
      DelAngle = OldAngle - NewAngle

   END DO

   RETURN
   END SUBROUTINE AddOrSub2Pi_R4
!=======================================================================
!> \copydoc nwtc_num::addorsub2pi_r4
   SUBROUTINE AddOrSub2Pi_R8 ( OldAngle, NewAngle )

      ! Argument declarations:

   REAL(R8Ki), INTENT(IN   )    :: OldAngle                                     ! Angle from which NewAngle will be converted to within Pi of, rad.
   REAL(R8Ki), INTENT(INOUT)    :: NewAngle                                     ! Angle to be converted to within Pi of OldAngle, rad.


      ! Local declarations:

   REAL(R8Ki)                   :: DelAngle                                     ! The difference between OldAngle and NewAngle, rad.



      ! Add or subtract 2*Pi in order to convert NewAngle two within Pi of OldAngle:

   
   DelAngle = OldAngle - NewAngle

   DO WHILE ( ABS( DelAngle ) > Pi_R8 )

      NewAngle = NewAngle + SIGN( TwoPi_R8, DelAngle )
      DelAngle = OldAngle - NewAngle

   END DO

   RETURN
   END SUBROUTINE AddOrSub2Pi_R8
!=======================================================================
   FUNCTION BlendCosine( x, LowerBound, UpperBound ) RESULT(S)
   
      REAL(ReKi), INTENT(IN) :: x            !
      REAL(ReKi), INTENT(IN) :: LowerBound   !< if x <= LowerBound, S=0 
      REAL(ReKi), INTENT(IN) :: UpperBound   !< if x >= UpperBound, S=1
      REAL(ReKi)             :: S 
      
      if (x >= UpperBound) then
         S = 1.0_ReKi
      elseif (x <= LowerBound) then
         S = 0.0_ReKi
      elseif (LowerBound < UpperBound) then
         S = 0.5_ReKi*(1.0_ReKi - cos((x-LowerBound) / (UpperBound-LowerBound)*pi))
      else ! can only get here if LowerBound>=UpperBound>x , which should be an error
         S = 0.0_ReKi
      end if
      
   END FUNCTION BlendCosine
!=======================================================================
!> This routine sorts a list of real numbers. It uses the bubble sort algorithm,
!! which is only suitable for short lists.
   SUBROUTINE BSortReal ( RealAry, NumPts )

      ! Argument declarations:

   INTEGER, INTENT(IN)          :: NumPts                                       !< The length of the list to be sorted.

   REAL(ReKi), INTENT(INOUT)    :: RealAry(NumPts)                              !< The list of real numbers to be sorted.


      ! Local declarations:

   REAL(ReKi)                   :: Temp                                         ! Temporary variable to hold the current element.

   INTEGER                      :: I                                            ! Index into the array.

   LOGICAL                      :: Change                                       ! Flag to indicate if a change of order was made.


      ! Sort the list

   Change = .TRUE.

   DO WHILE ( Change )

      Change = .FALSE.

      DO I=2,NumPts
         IF ( RealAry(I) < RealAry(I-1) )  THEN
            Temp           = RealAry(I)
            RealAry(I)   = RealAry(I-1)
            RealAry(I-1) = Temp
            Change         = .TRUE.
         END IF
      END DO ! I

   END DO ! WHILE


   RETURN
   END SUBROUTINE BSortReal ! ( RealAry, NumPts )
!=======================================================================
!> This subroutine takes an "oldUnits" array, compares the strings
!! to a list of units that will be converted to SI, and returns two arrays
!! that give the new units and the multiplicative scaling factor to convert 
!! the old units to the new ones. The three arrays must be the same size.
   SUBROUTINE ConvertUnitsToSI(Units,ScaleFactor)
      CHARACTER(*), INTENT(INOUT)   :: Units(:)          !< in: the old units; out: the new units
      REAL(ReKi),   INTENT(  OUT)   :: ScaleFactor(:)    !< scaling factor to convert old to new units (old*SF = new)

      
      ! local variables
      INTEGER                       :: i

      DO i=1,SIZE(Units)

         SELECT CASE( TRIM(Units(i)) )  ! Note that this IS case sensitive!

         CASE ('(kN)','kN')
            Units(i)    = '(N)'
            ScaleFactor(i) = 1000.0_ReKi
         CASE ('(kN-m)','kN-m')
            Units(i)    = '(N-m)'
            ScaleFactor(i) = 1000.0_ReKi
         CASE ('(deg)','deg')
            Units(i)    = '(rad)'
            ScaleFactor(i) = D2R
         CASE ('(deg/s)','deg/s')
            Units(i)    = '(rad/s)'
            ScaleFactor(i) = D2R
         CASE ('(deg/s^2)','deg/s^2')
            Units(i)    = '(rad/s^2)'
            ScaleFactor(i) = D2R
         CASE ('(rpm)','rpm')
            Units(i)    = '(rad/s)'
            ScaleFactor(i) = RPM2RPS
         CASE ('(kW)','kW')
            Units(i)    = '(W)'
            ScaleFactor(i) = 1000.0_ReKi
         CASE DEFAULT
            ScaleFactor(i) = 1.0_ReKi
         END SELECT

      END DO

   END SUBROUTINE ConvertUnitsToSI
!=======================================================================
!> This subroutine takes an "oldUnits" array, compares the strings
!! to a list of units that will be converted to engineering units (kN and deg), and returns two arrays
!! that give the new units and the multiplicative scaling factor to convert 
!! the old units to the new ones. The three arrays must be the same size.
   SUBROUTINE ConvertUnitsToEngr(Units,ScaleFactor)
      CHARACTER(*), INTENT(INOUT)   :: Units(:)          !< in: the old units; out: the new units
      REAL(ReKi),   INTENT(  OUT)   :: ScaleFactor(:)    !< scaling factor to convert old to new units (old*SF = new)

      
      ! local variables
      INTEGER                       :: i

      DO i=1,SIZE(Units)

         SELECT CASE( TRIM(Units(i)) )  ! Note that this IS case sensitive!

         CASE ('(N)','N')
            Units(i)    = '(kN)'
            ScaleFactor(i) = 0.001_ReKi
         CASE ('(N-m)','N-m', '(Nm)', 'Nm')
            Units(i)    = '(kN-m)'
            ScaleFactor(i) = 0.001_ReKi
         CASE ('(rad)','rad')
            Units(i)    = '(deg)'
            ScaleFactor(i) = R2D
         CASE ('(rad/s)','rad/s')
            Units(i)    = '(deg/s)'
            ScaleFactor(i) = R2D
         CASE ('(rad/s^2)','rad/s^2')
            Units(i)    = '(deg/s^2)'
            ScaleFactor(i) = R2D
         CASE ('(rps)','rps')
            Units(i)    = '(rpm)'
            ScaleFactor(i) = 60.0_ReKi
         CASE ('(W)','W')
            Units(i)    = '(kW)'
            ScaleFactor(i) = 0.001_ReKi
         CASE DEFAULT
            ScaleFactor(i) = 1.0_ReKi
         END SELECT

      END DO

   END SUBROUTINE ConvertUnitsToEngr
!=======================================================================
!> This function computes the cross product of two 3-element arrays (resulting in a vector): \n
!! cross_product = Vector1 \f$\times\f$ Vector2 \n
!! Use cross_product (nwtc_num::cross_product) instead of directly calling a specific routine in the generic interface.
   FUNCTION Cross_ProductR4(Vector1, Vector2) result(CProd)

      ! Argument declarations.

   REAL(SiKi), INTENT(IN )         :: Vector1       (3)
   REAL(SiKi), INTENT(IN )         :: Vector2       (3)

      ! Function definition
   REAL(SiKi)                      :: CProd (3)        ! = Vector1 X Vector2 (resulting in a vector)


   CProd(1) = Vector1(2)*Vector2(3) - Vector1(3)*Vector2(2)
   CProd(2) = Vector1(3)*Vector2(1) - Vector1(1)*Vector2(3)
   CProd(3) = Vector1(1)*Vector2(2) - Vector1(2)*Vector2(1)


   RETURN
   END FUNCTION Cross_ProductR4
!=======================================================================
!> \copydoc nwtc_num::cross_productr4
   FUNCTION Cross_ProductR4R8(Vector1, Vector2) result(CProd)

      ! Argument declarations.

   REAL(SiKi), INTENT(IN )         :: Vector1       (3)
   REAL(R8Ki), INTENT(IN )         :: Vector2       (3)

      ! Function definition
   REAL(R8Ki)                      :: CProd (3)        ! = Vector1 X Vector2 (resulting in a vector)


   CProd(1) = Vector1(2)*Vector2(3) - Vector1(3)*Vector2(2)
   CProd(2) = Vector1(3)*Vector2(1) - Vector1(1)*Vector2(3)
   CProd(3) = Vector1(1)*Vector2(2) - Vector1(2)*Vector2(1)


   RETURN
   END FUNCTION Cross_ProductR4R8
!=======================================================================
!> \copydoc nwtc_num::cross_productr4
   FUNCTION Cross_ProductR8(Vector1, Vector2) result(CProd)

      ! Argument declarations.

   REAL(R8Ki), INTENT(IN )         :: Vector1       (3)
   REAL(R8Ki), INTENT(IN )         :: Vector2       (3)

      ! Function definition
   REAL(R8Ki)                      :: CProd (3)        ! = Vector1 X Vector2 (resulting in a vector)


   CProd(1) = Vector1(2)*Vector2(3) - Vector1(3)*Vector2(2)
   CProd(2) = Vector1(3)*Vector2(1) - Vector1(1)*Vector2(3)
   CProd(3) = Vector1(1)*Vector2(2) - Vector1(2)*Vector2(1)


   RETURN
   END FUNCTION Cross_ProductR8
!=======================================================================
!> \copydoc nwtc_num::cross_productr4
   FUNCTION Cross_ProductR8R4(Vector1, Vector2) result(CProd)

      ! Argument declarations.

   REAL(R8Ki), INTENT(IN )         :: Vector1       (3)
   REAL(SiKi), INTENT(IN )         :: Vector2       (3)

      ! Function definition
   REAL(R8Ki)                      :: CProd (3)        ! = Vector1 X Vector2 (resulting in a vector)


   CProd(1) = Vector1(2)*Vector2(3) - Vector1(3)*Vector2(2)
   CProd(2) = Vector1(3)*Vector2(1) - Vector1(1)*Vector2(3)
   CProd(3) = Vector1(1)*Vector2(2) - Vector1(2)*Vector2(1)


   RETURN
   END FUNCTION Cross_ProductR8R4
!=======================================================================
!> This routine calculates the parameters needed to compute a irregularly-spaced natural cubic spline.
!! Natural cubic splines are used in that the curvature at the end points is zero.
!! This routine does not require that the XAry be regularly spaced.
   SUBROUTINE CubicSplineInit ( AryLen, XAry, YAry, Coef, ErrStat, ErrMsg )

      ! Argument declarations:

   INTEGER, INTENT(IN)          :: AryLen                                     !< Length of the array

   REAL(ReKi), INTENT(OUT)      :: Coef  (AryLen-1,0:3)                       !< The coefficients for the cubic polynomials
   REAL(ReKi), INTENT(IN)       :: XAry  (AryLen)                             !< Input array of x values
   REAL(ReKi), INTENT(IN)       :: YAry  (AryLen)                             !< Input array of y values

   INTEGER(IntKi), INTENT(OUT)  :: ErrStat                                    !< Error status

   CHARACTER(*), INTENT(OUT)    :: ErrMsg                                     !< Error message


      ! Local declarations.

   REAL(ReKi), ALLOCATABLE      :: DelX  (:)                                  ! The distances between the randomly spaced points.
   REAL(ReKi), ALLOCATABLE      :: Slope (:)                                  ! The AryLen-1 length array of slopes between points.
   REAL(ReKi), ALLOCATABLE      :: U     (:)                                  ! An AryLen-1 length array used in the Gaussian elimination.
   REAL(ReKi), ALLOCATABLE      :: V     (:)                                  ! An AryLen-1 length array used in the Gaussian elimination.
   REAL(ReKi)                   :: ZHi                                        ! A parameter used to calculate the polynomial coefficients.
   REAL(ReKi)                   :: ZLo                                        ! A parameter used to calculate the polynomial coefficients.

   INTEGER(IntKi)               :: ErrStatLcL                                 ! Local error status.
   INTEGER                      :: I                                          ! The index into the arrays.
   CHARACTER(*), PARAMETER      :: RoutineName = 'CubicSplineInit'

   ErrStat = ErrID_None
   ErrMsg  = ""

      ! Allocate the various intermediate arrays.

   ALLOCATE ( DelX( AryLen - 1 ), STAT=ErrStatLcL )
   IF ( ErrStatLcL /= 0 )  THEN
      CALL ExitThisRoutine ( ErrID_Fatal, RoutineName//':Error allocating memory for the DelX array.' )
      RETURN
   ENDIF

   ALLOCATE ( Slope( AryLen - 1 ), STAT=ErrStatLcL )
   IF ( ErrStatLcL /= 0 )  THEN
      CALL ExitThisRoutine ( ErrID_Fatal, RoutineName//':Error allocating memory for the Slope array.' )
      RETURN
   ENDIF

   ALLOCATE ( U( AryLen - 1 ), STAT=ErrStatLcL )
   IF ( ErrStatLcL /= 0 )  THEN
      CALL ExitThisRoutine ( ErrID_Fatal, RoutineName//':Error allocating memory for the U array.' )
      RETURN
   ENDIF

   ALLOCATE ( V( AryLen - 1 ), STAT=ErrStatLcL )
   IF ( ErrStatLcL /= 0 )  THEN
      CALL ExitThisRoutine ( ErrID_Fatal, RoutineName//':Error allocating memory for the V array.' )
      RETURN
   ENDIF


      ! Compute the distance between XAry values and the slopes between points.

   DO I=1,AryLen-1
      DelX (I) =   XAry(I+1) - XAry(I)
      Slope(I) = ( YAry(I+1) - YAry(I) )/DelX(I)
   END DO ! I


      ! Use Gaussian elimination to solve the tri-diagonal matrix.

   U(1) = 2.0_ReKi*( DelX (2) + DelX (1) )
   V(1) = 6.0_ReKi*( Slope(2) - Slope(1) )

   DO I=2,AryLen-1
      U(I) = 2.0_ReKi*( DelX(I-1) + DelX(I)    ) - DelX(I-1)*DelX(I-1)/U(I-1)
      V(I) = 6.0_ReKi*( Slope(I)  - Slope(I-1) ) - DelX(I-1)*   V(I-1)/U(I-1)
   END DO ! I


      ! Determine the coefficients of the polynomials.

   Coef(:,0) = YAry(:)

   ZHi = 0.0_ReKi

   DO I=AryLen-1,1,-1
      ZLo       = ( V(I) - DelX(I)*ZHi )/U(I)
      Coef(I,1) = Slope(I) - DelX(I)*( ZHi/6.0_ReKi + ZLo/3.0_ReKi )
      Coef(I,2) = 0.5_ReKi*ZLo
      Coef(I,3) = ( ZHi - ZLo )/( 6.0_ReKi*DelX(I) )
      ZHi       = ZLo
   END DO ! I



   CALL ExitThisRoutine ( ErrID_None, 'No Problemo' )

   RETURN

   !=======================================================================
   CONTAINS
   !=======================================================================
      SUBROUTINE ExitThisRoutine ( ErrID, Msg )

         ! This subroutine cleans up the parent routine before exiting.


            ! Argument declarations.

         INTEGER(IntKi), INTENT(IN)       :: ErrID                            ! The error identifier (ErrLev)

         CHARACTER(*),   INTENT(IN)       :: Msg                              ! The error message (ErrMsg)


            ! Local declarations.

         LOGICAL                          :: IsOpen                           ! A flag that indicates if the input unit is still open.


            ! Set error status/message

         ErrStat = ErrID
         ErrMsg  = Msg


            ! Deallocate the Words array if it had been allocated.

         IF ( ALLOCATED( DelX  ) )  DEALLOCATE( DelX  )
         IF ( ALLOCATED( Slope ) )  DEALLOCATE( Slope )
         IF ( ALLOCATED( U     ) )  DEALLOCATE( U     )
         IF ( ALLOCATED( V     ) )  DEALLOCATE( V     )


         RETURN

      END SUBROUTINE ExitThisRoutine ! ( ErrID, Msg )

   END SUBROUTINE CubicSplineInit ! ( AryLen, XAry, YAry, YAry, Coef, ErrStat, ErrMsg )
!=======================================================================
!> This routine calculates the parameters needed to compute a irregularly-spaced natural cubic spline.
!! Natural cubic splines are used in that the curvature at the end points is zero.
!! This routine does not require that the XAry be regularly spaced.
!! This version of the routine works with multiple curves that share the same X values.
   SUBROUTINE CubicSplineInitM ( XAry, YAry, Coef, ErrStat, ErrMsg )

      ! Argument declarations:

   REAL(ReKi), INTENT(OUT)      :: Coef  (:,:,0:)                             !< The coefficients for the cubic polynomials
   REAL(ReKi), INTENT(IN)       :: XAry  (:)                                  !< Input array of x values
   REAL(ReKi), INTENT(IN)       :: YAry  (:,:)                                !< Input array of y values with multiple curves

   INTEGER(IntKi), INTENT(OUT)  :: ErrStat                                    !< Error status

   CHARACTER(*), INTENT(OUT)    :: ErrMsg                                     !< Error message


      ! Local declarations:

   REAL(ReKi), ALLOCATABLE      :: DelX  (:)                                  ! The distances between the randomly spaced points.
   REAL(ReKi), ALLOCATABLE      :: Slope (:,:)                                ! The NumPts-1 length array of slopes between points.
   REAL(ReKi), ALLOCATABLE      :: U     (:)                                  ! An NumPts-1 length array used in the Gaussian elimination.
   REAL(ReKi), ALLOCATABLE      :: V     (:,:)                                ! An NumPts-1 by NumCrvs length array used in the Gaussian elimination.
   REAL(ReKi), ALLOCATABLE      :: ZHi   (:)                                  ! A parameter used to calculate the polynomial coefficients.
   REAL(ReKi), ALLOCATABLE      :: ZLo   (:)                                  ! A parameter used to calculate the polynomial coefficients.

   INTEGER(IntKi)               :: ErrStatLcL                                 ! Local error status.

   INTEGER                      :: I                                          ! The index into the arrays.
   INTEGER                      :: NumCrvs                                    ! Number of curves to be interpolated.
   INTEGER                      :: NumPts                                     ! Number of points in each curve.
   CHARACTER(*), PARAMETER      :: RoutineName = 'CubicSplineInitM'


      ! How big are the arrays?

   NumPts  = SIZE( XAry )
   NumCrvs = SIZE( YAry, 2 )


      ! Allocate the various intermediate arrays.

   ALLOCATE ( ZLo( NumCrvs ), STAT=ErrStatLcL )
   IF ( ErrStatLcL /= 0 )  THEN
      CALL ExitThisRoutine ( ErrID_Fatal, RoutineName//':Error allocating memory for the ZLo array.' )
      RETURN
   ENDIF

   ALLOCATE ( ZHi( NumCrvs ), STAT=ErrStatLcL )
   IF ( ErrStatLcL /= 0 )  THEN
      CALL ExitThisRoutine ( ErrID_Fatal, RoutineName//':Error allocating memory for the ZHi array.' )
      RETURN
   ENDIF

   ALLOCATE ( DelX( NumPts - 1 ), STAT=ErrStatLcL )
   IF ( ErrStatLcL /= 0 )  THEN
      CALL ExitThisRoutine ( ErrID_Fatal, RoutineName//':Error allocating memory for the DelX array.' )
      RETURN
   ENDIF

   ALLOCATE ( Slope( NumPts-1, NumCrvs ), STAT=ErrStatLcL )
   IF ( ErrStatLcL /= 0 )  THEN
      CALL ExitThisRoutine ( ErrID_Fatal, RoutineName//':Error allocating memory for the Slope array.' )
      RETURN
   ENDIF

   ALLOCATE ( U( NumPts - 1 ), STAT=ErrStatLcL )
   IF ( ErrStatLcL /= 0 )  THEN
      CALL ExitThisRoutine ( ErrID_Fatal, RoutineName//':Error allocating memory for the U array.' )
      RETURN
   ENDIF

   ALLOCATE ( V( NumPts-1, NumCrvs ), STAT=ErrStatLcL )
   IF ( ErrStatLcL /= 0 )  THEN
      CALL ExitThisRoutine ( ErrID_Fatal, RoutineName//':Error allocating memory for the V array.' )
      RETURN
   ENDIF


      ! Compute the distance between XAry values and the slopes between points.

   DO I=1,NumPts-1
      DelX (I  ) =   XAry(I+1  ) - XAry(I  )
      
      if ( equalRealNos( DelX(I), 0.0_ReKi ) ) then
         CALL ExitThisRoutine ( ErrID_Fatal, RoutineName//':XAry must have unique values.' )
         RETURN
      ENDIF
      
      Slope(I,:) = ( YAry(I+1,:) - YAry(I,:) )/DelX(I)
   END DO ! I


      ! Use Gaussian elimination to solve the tri-diagonal matrix.

   U(1  ) = 2.0_ReKi*( DelX (2)   + DelX (1  ) )
   V(1,:) = 6.0_ReKi*( Slope(2,:) - Slope(1,:) )

   DO I=2,NumPts-1
      if ( equalRealNos( U(I-1), 0.0_ReKi ) ) then
         CALL ExitThisRoutine ( ErrID_Fatal, RoutineName//':XAry must be monotonic.' )
         RETURN
      ENDIF
            
      U(I)   = 2.0_ReKi*( DelX (I-1) + DelX (I)     ) - DelX(I-1)*DelX(I-1  )/U(I-1)
      V(I,:) = 6.0_ReKi*( Slope(I,:) - Slope(I-1,:) ) - DelX(I-1)*   V(I-1,:)/U(I-1)
   END DO ! I


      ! Determine the coefficients of the polynomials.

   Coef(:,:,0) = YAry(1:NumPts-1,:)

   ZHi(:) = 0.0_ReKi

   DO I=NumPts-1,1,-1
      ZLo(:)      = ( V(I,:) - DelX(I)*ZHi(:) )/U(I)                             ! bjj: already checked for u(I) == 0
      Coef(I,:,1) = Slope(I,:) - DelX(I)*( ZHi(:)/6.0_ReKi + ZLo(:)/3.0_ReKi )
      Coef(I,:,2) = 0.5_ReKi*ZLo(:)
      Coef(I,:,3) = ( ZHi(:) - ZLo(:) )/( 6.0_ReKi*DelX(I) )                     ! bjj: already checked for DelX(I) == 0
      ZHi(:)      = ZLo(:)
   END DO ! I
   
   CALL ExitThisRoutine ( ErrID_None, 'No Problemo' )

   RETURN

   !=======================================================================
   CONTAINS
   !=======================================================================
      SUBROUTINE ExitThisRoutine ( ErrID, Msg )

         ! This subroutine cleans up the parent routine before exiting.


            ! Argument declarations.

         INTEGER(IntKi), INTENT(IN)       :: ErrID                            ! The error identifier (ErrLev)

         CHARACTER(*),   INTENT(IN)       :: Msg                              ! The error message (ErrMsg)


            ! Local declarations.

         LOGICAL                          :: IsOpen                           ! A flag that indicates if the input unit is still open.


            ! Set error status/message

         ErrStat = ErrID
         ErrMsg  = Msg


            ! Deallocate the Words array if it had been allocated.

         IF ( ALLOCATED( DelX  ) )  DEALLOCATE( DelX  )
         IF ( ALLOCATED( Slope ) )  DEALLOCATE( Slope )
         IF ( ALLOCATED( U     ) )  DEALLOCATE( U     )
         IF ( ALLOCATED( V     ) )  DEALLOCATE( V     )


         RETURN

      END SUBROUTINE ExitThisRoutine ! ( ErrID, Msg )

   END SUBROUTINE CubicSplineInitM ! ( XAry, YAry, Coef, ErrStat, ErrMsg )
!=======================================================================
!> This routine calculates the parameters needed to compute a irregularly-spaced natural linear spline.      
!! This routine does not require that the XAry be regularly spaced.
   SUBROUTINE CubicLinSplineInitM ( XAry, YAry, Coef, ErrStat, ErrMsg )

      ! Argument declarations:

   REAL(ReKi), INTENT(OUT)      :: Coef  (:,:,0:)                             !< The coefficients for the cubic polynomials
   REAL(ReKi), INTENT(IN)       :: XAry  (:)                                  !< Input array of x values
   REAL(ReKi), INTENT(IN)       :: YAry  (:,:)                                !< Input array of y values with multiple curves

   INTEGER(IntKi), INTENT(OUT)  :: ErrStat                                    !< Error status

   CHARACTER(*), INTENT(OUT)    :: ErrMsg                                     !< Error message


      ! Local declarations.

   REAL(ReKi)                   :: DelX                                       ! The distances between the randomly spaced points.

   INTEGER                      :: I                                          ! The index into the arrays.
   INTEGER                      :: NumPts                                     ! Number of points in each curve.
   CHARACTER(*), PARAMETER      :: RoutineName = 'CubicLinSplineInitM'


      ! How big are the arrays?

   NumPts  = SIZE( XAry )
   ErrStat = ErrID_None
   ErrMsg  = ""

      ! Determine the coefficients of the polynomials.

   
   DO I=NumPts-1,1,-1
      DelX =   XAry(I+1  ) - XAry(I  )
      
      if ( equalRealNos( DelX, 0.0_ReKi ) ) then
         CALL SetErrStat ( ErrID_Fatal, 'XAry must have unique values.',ErrStat,ErrMsg,RoutineName )
         RETURN
      ENDIF
      
            
      Coef(I,:,0) = YAry(I,:)
      Coef(I,:,1) = (YAry(I+1,: ) - YAry(I,:  )) / DelX
      Coef(I,:,2) = 0.0_ReKi
      Coef(I,:,3) = 0.0_ReKi
   END DO ! I


   RETURN


   END SUBROUTINE CubicLinSplineInitM ! ( XAry, YAry, Coef, ErrStat, ErrMsg )
!=======================================================================
!> This routine interpolates a pair of arrays using cubic splines to find the function value at X.
!! One must call cubicsplineinit first to compute the coefficients of the cubics.
!! This routine does not require that the XAry be regularly spaced.
   FUNCTION CubicSplineInterp ( X, AryLen, XAry, YAry, Coef, ErrStat, ErrMsg )

      ! Function declaration.

   REAL(ReKi)                   :: CubicSplineInterp                          !  This function


      ! Argument declarations:

   INTEGER, INTENT(IN)          :: AryLen                                     !< Length of the array

   REAL(ReKi), INTENT(IN)       :: Coef  (AryLen-1,0:3)                       !< The coefficients for the cubic polynomials
   REAL(ReKi), INTENT(IN)       :: X                                          !< The value we are trying to interpolate for
   REAL(ReKi), INTENT(IN)       :: XAry (AryLen)                              !< Input array of regularly spaced x values
   REAL(ReKi), INTENT(IN)       :: YAry (AryLen)                              !< Input array of y values

   INTEGER(IntKi), INTENT(OUT)  :: ErrStat                                    !< Error status

   CHARACTER(*), INTENT(OUT)    :: ErrMsg                                     !< Error message


      ! Local declarations.

   REAL(ReKi)                   :: XOff                                       ! The distance from X to XAry(ILo).

!   INTEGER(IntKi)               :: ErrStatLcL                                 ! Local error status.
   INTEGER                      :: ILo                                        ! The index into the array for which X is just above or equal to XAry(ILo).


   ErrStat = ErrID_None
   ErrMsg  = ""

      ! See if X is within the range of XAry.  Return the end point if it is not.

   IF ( X <= XAry(1) )  THEN
      CubicSplineInterp = YAry(1)
      RETURN
   ELSEIF ( X >= XAry(AryLen) )  THEN
      CubicSplineInterp = YAry(AryLen)
      RETURN
   ENDIF ! ( X <= XAry(1) )


      ! We are somewhere inside XAry.  Find the segment that bounds X using binary search.

   CALL LocateBin( X, XAry, ILo, AryLen )

   XOff = X - XAry(ILo)

   CubicSplineInterp = Coef(ILo,0) + XOff*( Coef(ILo,1) + XOff*( Coef(ILo,2) + XOff*Coef(ILo,3) ) )


   RETURN
   END FUNCTION CubicSplineInterp ! ( X, AryLen, XAry, YAry, Coef, ErrStat, ErrMsg )
!=======================================================================
!> This routine interpolates a pair of arrays using cubic splines to find the function value at X.
!! One must call cubicsplineinit first to compute the coefficients of the cubics.
!! This routine does not require that the XAry be regularly spaced.
!! This version of the routine works with multiple curves that share the same X values.
   FUNCTION CubicSplineInterpM ( X, XAry, YAry, Coef, ErrStat, ErrMsg ) RESULT( Res )

      ! Function declaration.

   REAL(ReKi), ALLOCATABLE      :: Res(:)                                     ! The result of this function


      ! Argument declarations:

   REAL(ReKi), INTENT(IN)       :: Coef  (:,:,0:)                             !< The coefficients for the cubic polynomials
   REAL(ReKi), INTENT(IN)       :: X                                          !< The value we are trying to interpolate for
   REAL(ReKi), INTENT(IN)       :: XAry (:)                                   !< Input array of regularly spaced x values
   REAL(ReKi), INTENT(IN)       :: YAry (:,:)                                 !< Input array of y values with multiple curves

   INTEGER(IntKi), INTENT(OUT)  :: ErrStat                                    !< Error status

   CHARACTER(*),    INTENT(OUT) :: ErrMsg                                     !< Error message


      ! Local declarations.

   REAL(ReKi)                   :: XOff                                       ! The distance from X to XAry(ILo).

   INTEGER(IntKi)               :: ErrStatLcL                                 ! Local error status.
   INTEGER                      :: ILo                                        ! The index into the array for which X is just above or equal to XAry(ILo).
   INTEGER                      :: NumCrvs                                    ! Number of curves to be interpolated.
   INTEGER                      :: NumPts                                     ! Number of points in each curve.

   CHARACTER(*), PARAMETER      :: RoutineName = 'RegCubicSplineInterpM'

      ! How big are the arrays?

   NumPts  = SIZE( XAry )
   NumCrvs = SIZE( YAry, 2 )

   ALLOCATE ( Res( NumCrvs ) , STAT=ErrStatLcl )
   IF ( ErrStatLcl /= 0 )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg = RoutineName//':Error allocating memory for the function result array.'
      RETURN
   ELSE
      ErrStat = ErrID_None
      ErrMsg  = ""
   ENDIF


      ! See if X is within the range of XAry.  Return the end point if it is not.

   IF ( X <= XAry(1) )  THEN
      Res(:) = YAry(1,:)
      RETURN
   ELSEIF ( X >= XAry(NumPts) )  THEN
      Res(:) = YAry(NumPts,:)
      RETURN
   ENDIF ! ( X <= XAry(1) )


      ! We are somewhere inside XAry.  Find the segment that bounds X using binary search.

   CALL LocateBin( X, XAry, ILo, NumPts )

   XOff = X - XAry(ILo)

   Res(:) = Coef(ILo,:,0) + XOff*( Coef(ILo,:,1) + XOff*( Coef(ILo,:,2) + XOff*Coef(ILo,:,3) ) )

   RETURN

   END FUNCTION CubicSplineInterpM ! ( X, XAry, YAry, Coef, ErrStat, ErrMsg )
!=======================================================================         
!> This function returns the matrix exponential, \f$\Lambda = \exp(\lambda)\f$, of an input skew-symmetric matrix, \f$\lambda\f$.
!!
!! \f$\lambda\f$ is defined as:
!!
!! \f{equation}{  \lambda = \begin{bmatrix}
!!  0          &  \lambda_3 & -\lambda_2 \\
!!  -\lambda_3 &  0         &  \lambda_1 \\
!!   \lambda_2 & -\lambda_1 &  0          
!! 	\end{bmatrix}
!! \f}   
!! The angle of rotation for \f$\lambda\f$ is 
!! \f{equation}{ \theta = \sqrt{{\lambda_1}^2+{\lambda_2}^2+{\lambda_3}^2} \f}
!!
!! The matrix exponential is calculated as
!! \f{equation}{
!!  \Lambda = \exp(\lambda) = \left\{ \begin{matrix}
!!  I                                                                              &  \theta = 0 \\
!!  I + \frac{\sin\theta}{\theta}\lambda + \frac{1-\cos\theta}{\theta^2}\lambda^2  &  \theta > 0 
!!  \end{matrix}  \right.
!! \f}
!!
!! This routine is the inverse of DCM_logMap (nwtc_num::dcm_logmap). \n
!! Use DCM_exp (nwtc_num::dcm_exp) instead of directly calling a specific routine in the generic interface.   
   FUNCTION DCM_expD(lambda)
   
      
   REAL(DbKi), INTENT(IN)  :: lambda(3)            !< vector containing \f$\lambda_1\f$, \f$\lambda_2\f$, and \f$\lambda_3\f$, the unique components of skew-symmetric matrix \f$\lambda\f$ 
   REAL(DbKi)              :: DCM_expD(3,3)        !< the computed matrix exponential, \f$\Lambda\f$
   
      ! local variables
   REAL(DbKi)              :: stheta         ! sine of angle of rotation   
   REAL(DbKi)              :: theta          ! angle of rotation   
   REAL(DbKi)              :: theta2         ! angle of rotation squared
   REAL(DbKi)              :: tmp_Mat(3,3)
   
   INTEGER(IntKi)          :: ErrStat
   CHARACTER(30)           :: ErrMsg  
   
   
   theta = TwoNorm(lambda)                   ! Eq. 32
   theta2 = theta**2

   IF ( EqualRealNos(theta,   0.0_DbKi)   .or. &
        EqualRealNos(theta2,  0.0_DbKi) ) THEN  !
      
      CALL eye(DCM_expD, ErrStat, ErrMsg)    ! Eq. 33a
      
   ELSE   
      
         ! convert lambda to skew-symmetric matrix:
      tmp_mat(1,1) =  0.0_DbKi                                            
      tmp_mat(2,1) = -lambda(3)                                           
      tmp_mat(3,1) =  lambda(2)                                           
      tmp_mat(1,2) =              lambda(3)                               
      tmp_mat(2,2) =              0.0_DbKi                                
      tmp_mat(3,2) =             -lambda(1)                               
      tmp_mat(1,3) =                               -lambda(2)             
      tmp_mat(2,3) =                                lambda(1)             
      tmp_mat(3,3) =                                0.0_DbKi            
      
      
         ! Eq. 33b
      !DCM_exp = I + sin(theta)/theta*tmp_mat + (1-cos(theta))/theta**2)*matmul(tmp_mat,tmp_mat)
      
         ! one method:
      !CALL eye(DCM_exp, ErrStat, ErrMsg)                  
      !DCM_exp = DCM_exp + sin(theta)/theta*tmp_mat 
      !DCM_exp = DCM_exp + (1-cos(theta))/theta2 * MATMUL(tmp_mat, tmp_mat) 
      
         ! hopefully this order of calculations gives better numerical results:
      stheta = sin(theta)
      DCM_expD      = (1-cos(theta))/theta * tmp_mat      
      DCM_expD(1,1) = DCM_expD(1,1) + stheta
      DCM_expD(2,2) = DCM_expD(2,2) + stheta
      DCM_expD(3,3) = DCM_expD(3,3) + stheta
      
      DCM_expD = matmul( DCM_expD, tmp_mat )
      DCM_expD = DCM_expD / theta
      DCM_expD(1,1) = DCM_expD(1,1) + 1.0_DbKi ! add identity
      DCM_expD(2,2) = DCM_expD(2,2) + 1.0_DbKi
      DCM_expD(3,3) = DCM_expD(3,3) + 1.0_DbKi
            
   END IF

   
      
   END FUNCTION DCM_expD
!=======================================================================  
!> \copydoc nwtc_num::dcm_expd
   FUNCTION DCM_expR(lambda)
      ! This function computes a matrix exponential.
      !
      ! "'Interpolation' of DCMs", M.A. Sprague, 11 March 2014, Eq. 31-33
      
   REAL(SiKi), INTENT(IN)  :: lambda(3)      !< vector containing unique components of skew-symmetric matrix: \f$\lambda_1\f$, \f$\lambda_2\f$, and \f$\lambda_3\f$
   REAL(SiKi)              :: DCM_expR(3,3)  !< the computed matrix exponential, \f$\Lambda\f$
   
      ! local variables
   REAL(SiKi)              :: stheta         ! sine of angle of rotation   
   REAL(SiKi)              :: theta          ! angle of rotation   
   REAL(SiKi)              :: theta2         ! angle of rotation squared
   REAL(SiKi)              :: tmp_Mat(3,3)
   
   INTEGER(IntKi)          :: ErrStat
   CHARACTER(30)           :: ErrMsg  
   
   
   theta = TwoNorm(lambda)                   ! Eq. 32
   theta2 = theta**2

   IF ( EqualRealNos(theta,   0.0_SiKi)   .or. &
        EqualRealNos(theta2,  0.0_SiKi) ) THEN  !
      
      CALL eye(DCM_expR, ErrStat, ErrMsg)    ! Eq. 33a
      
   ELSE   
      
         ! convert lambda to skew-symmetric matrix:
      !tmp_mat = -SkewSymMat(lambda)
      tmp_mat(1,1) =  0.0_SiKi                                            
      tmp_mat(2,1) = -lambda(3)                                           
      tmp_mat(3,1) =  lambda(2)                                           
      tmp_mat(1,2) =              lambda(3)                               
      tmp_mat(2,2) =              0.0_SiKi                                
      tmp_mat(3,2) =             -lambda(1)                               
      tmp_mat(1,3) =                               -lambda(2)             
      tmp_mat(2,3) =                                lambda(1)             
      tmp_mat(3,3) =                                0.0_SiKi            
      
      
         ! Eq. 33b
      !DCM_exp = I + sin(theta)/theta*tmp_mat + (1-cos(theta))/theta**2)*matmul(tmp_mat,tmp_mat)
      
         ! one method:
      !CALL eye(DCM_exp, ErrStat, ErrMsg)                  
      !DCM_exp = DCM_exp + sin(theta)/theta*tmp_mat 
      !DCM_exp = DCM_exp + (1-cos(theta))/theta2 * MATMUL(tmp_mat, tmp_mat) 
      
         ! hopefully this order of calculations gives better numerical results:
      stheta = sin(theta)
      DCM_expR      = (1-cos(theta))/theta * tmp_mat      
      DCM_expR(1,1) = DCM_expR(1,1) + stheta
      DCM_expR(2,2) = DCM_expR(2,2) + stheta
      DCM_expR(3,3) = DCM_expR(3,3) + stheta
      
      DCM_expR = matmul( DCM_expR, tmp_mat )
      DCM_expR = DCM_expR / theta
      DCM_expR(1,1) = DCM_expR(1,1) + 1.0_ReKi ! add identity
      DCM_expR(2,2) = DCM_expR(2,2) + 1.0_ReKi
      DCM_expR(3,3) = DCM_expR(3,3) + 1.0_ReKi
            
   END IF

   
      
   END FUNCTION DCM_expR
!=======================================================================  
!> For any direction cosine matrix (DCM), \f$\Lambda\f$, this routine calculates the
!! logarithmic map, \f$\lambda\f$, which a skew-symmetric matrix:
!!
!! \f{equation}{
!! \lambda 
!! = \log( \Lambda )
!! = \begin{bmatrix}
!!       0          &  \lambda_3 & -\lambda_2 \\
!!      -\lambda_3  &  0         &  \lambda_1 \\
!!       \lambda_2 & -\lambda_1 &  0          
!! \end{bmatrix}
!! \f}
!! The angle of rotation for \f$\Lambda\f$ is
!! \f{equation}{
!! \theta= \begin{matrix} \cos^{-1}\left(\frac{1}{2}\left(\mathrm{trace}(\Lambda)-1\right)\right) & \theta \in \left[0,\pi\right]\end{matrix}
!! \f}
!! And the logarithmic map is
!! \f{equation}{
!!  \lambda = \left\{ \begin{matrix}
!! 0                                                             &  \theta = 0 \\
!! \frac{\theta}{2\sin\theta} \left( \Lambda - \Lambda^T\right)  &  \theta \in \left(0,\pi\right) \\
!! \pm\pi v  &  \theta = \pi 
!!  \end{matrix}  \right.
!! \f}
!! where \f$v\f$ is the skew-symmetric matrix associated with the unit-length eigenvector of \f$\Lambda\f$ associated with the eigenvalue 1.
!! However, this equation has numerical issues near \f$\theta = \pi\f$, so for \f$\theta > 3.1\f$  we instead implement
!! a separate equation to find lambda * sign(lambda(indx_max))
!! and use \f$\Lambda - \Lambda^T\f$ to choose the appropriate signs. 
!!   
!! This routine is the inverse of DCM_exp (nwtc_num::dcm_exp). \n
!! Use DCM_logMap (nwtc_num::dcm_logmap) instead of directly calling a specific routine in the generic interface. 
   SUBROUTINE DCM_logMapD(DCM, logMap, ErrStat, ErrMsg, thetaOut)
   
   REAL(DbKi),         INTENT(IN)    :: DCM(3,3)                  !< the direction cosine matrix, \f$\Lambda\f$             
   REAL(DbKi),         INTENT(  OUT) :: logMap(3)                 !< vector containing \f$\lambda_1\f$, \f$\lambda_2\f$, and \f$\lambda_3\f$, the unique components of skew-symmetric matrix \f$\lambda\f$ 
   REAL(DbKi),OPTIONAL,INTENT(  OUT) :: thetaOut                  !< the angle of rotation, \f$\theta\f$; output only for debugging
   INTEGER(IntKi),     INTENT(  OUT) :: ErrStat                   !< Error status of the operation
   CHARACTER(*),       INTENT(  OUT) :: ErrMsg                    !< Error message if ErrStat /= ErrID_None
   
      ! local variables
   REAL(DbKi)                        :: theta
   REAL(DbKi)                        :: cosTheta
   REAL(DbKi)                        :: TwoSinTheta
   REAL(DbKi)                        :: v(3)
   REAL(DbKi)                        :: divisor
   INTEGER(IntKi)                    :: indx_max
      
         ! initialization
      ErrStat = ErrID_None
      ErrMsg  = ""   
   
   
      cosTheta = 0.5_DbKi*( trace(DCM) - 1.0_DbKi )
      cosTheta = min( max(cosTheta,-1.0_DbKi), 1.0_DbKi ) !make sure it's in a valid range (to avoid cases where this is slightly outside the +/-1 range)
      theta    = ACOS( cosTheta )                                                   ! Eq. 25 ( 0<=theta<=pi )

      IF ( PRESENT( thetaOut ) ) THEN
         thetaOut = theta
      END IF      
      
      !> Note that \f$ DCM = \begin{bmatrix}
      !!     1-\frac{1-\cos\theta}{\theta^2}\left( \lambda_3^2 + \lambda_2^2\right) 
      !!   &   \frac{\sin\theta}{\theta}\lambda_3+\frac{1-\cos\theta}{\theta^2}\lambda_1\lambda_2  
      !!   &  -\frac{\sin\theta}{\theta}\lambda_2+\frac{1-\cos\theta}{\theta^2}\lambda_1\lambda_3 \\
      !!      -\frac{\sin\theta}{\theta}\lambda_3+\frac{1-\cos\theta}{\theta^2}\lambda_1\lambda_2
      !!   &  1-\frac{1-\cos\theta}{\theta^2}\left( \lambda_3^2 + \lambda_1^2\right) 
      !!   &  \frac{\sin\theta}{\theta}\lambda_1+\frac{1-\cos\theta}{\theta^2}\lambda_2\lambda_3  \\
      !!      \frac{\sin\theta}{\theta}\lambda_2+\frac{1-\cos\theta}{\theta^2}\lambda_1\lambda_3 
      !!   & -\frac{\sin\theta}{\theta}\lambda_1+\frac{1-\cos\theta}{\theta^2}\lambda_2\lambda_3 
      !!   &  1-\frac{1-\cos\theta}{\theta^2}\left( \lambda_2^2 + \lambda_1^2\right)                    \\
      !!    \end{bmatrix} \f$
      
      
      !IF ( EqualRealNos( pi_D, theta )  ) THEN
      IF ( theta > 3.1_DbKi ) THEN  ! theta/(2*sin(theta)) blows up quickly as theta approaches pi, 
         ! so I'm putting a pretty large tolerance on pi here, and using a different equation to find the solution near pi
       
         logMap(1) = 1.0_DbKi + DCM(1,1) - DCM(2,2) - DCM(3,3);
         logMap(2) = 1.0_DbKi - DCM(1,1) + DCM(2,2) - DCM(3,3);
         logMap(3) = 1.0_DbKi - DCM(1,1) - DCM(2,2) + DCM(3,3);
             
         indx_max = maxloc( abs(logMap), 1 )
             
         divisor = sqrt(abs( logMap(indx_max) *  2.0_DbKi*(1.0_DbKi - cosTheta)  )) / theta  ! 2*(1-cosTheta)/theta^2 * abs(lambda(indx_max))
         if (indx_max == 1) then
           !logMap(1) = 1.0 + DCM(1,1) - DCM(2,2) - DCM(3,3)               ! 2*(1-cosTheta)/theta^2 * lambda(1) * lambda(1)
            logMap(2) = DCM(1,2) + DCM(2,1)                                ! 2*(1-cosTheta)/theta^2 * lambda(1) * lambda(2)
            logMap(3) = DCM(1,3) + DCM(3,1)                                ! 2*(1-cosTheta)/theta^2 * lambda(1) * lambda(3)
         elseif (indx_max == 2) then
            logMap(1) = DCM(1,2) + DCM(2,1)                                ! 2*(1-cosTheta)/theta^2 * lambda(2) * lambda(1)
           !logMap(2) = 1.0 - DCM(1,1) + DCM(2,2) - DCM(3,3)               ! 2*(1-cosTheta)/theta^2 * lambda(2) * lambda(2)
            logMap(3) = DCM(2,3) + DCM(3,2)                                ! 2*(1-cosTheta)/theta^2 * lambda(2) * lambda(3)
         else
            logMap(1) = DCM(1,3) + DCM(3,1)                                ! 2*(1-cosTheta)/theta^2 * lambda(3) * lambda(1)
            logMap(2) = DCM(2,3) + DCM(3,2)                                ! 2*(1-cosTheta)/theta^2 * lambda(3) * lambda(2)
           !logMap(3) = 1.0 - DCM(1,1) - DCM(2,2) + DCM(3,3)               ! 2*(1-cosTheta)/theta^2 * lambda(3) * lambda(3)
         end if
         logMap = logMap / divisor                                         ! lambda * sign(lambda(indx_max))

         ! at this point we may have the wrong sign for logMap (though if theta==pi, it doesn't matter because we can change it in the DCM_setLogMapforInterp() routines)
         ! we'll do a little checking to see if we should change the sign:
         
         IF ( EqualRealNos( pi_D, theta )  ) RETURN
         
         v(1) = -DCM(3,2) + DCM(2,3) !-skewSym(3,2) = 2*sin(theta)/theta * lambda(1) = (small positive value with theta near pi) * lambda(1)
         v(2) =  DCM(3,1) - DCM(1,3) ! skewSym(3,1) = 2*sin(theta)/theta * lambda(2) = (small positive value with theta near pi) * lambda(2)
         v(3) = -DCM(2,1) + DCM(1,2) !-skewSym(2,1) = 2*sin(theta)/theta * lambda(3) = (small positive value with theta near pi) * lambda(3)
         
         indx_max = maxloc( abs(v), 1 )  ! find component with largest magnitude
         if ( .not. EqualRealNos( sign(1.0_DbKi,v(indx_max)), sign(1.0_DbKi,logMap(indx_max)) )) logMap = -logMap
         
      ELSE
         
         TwoSinTheta = 2.0_DbKi*sin(theta)
         
         IF ( EqualRealNos(0.0_DbKi, theta) .or. EqualRealNos( 0.0_DbKi, TwoSinTheta ) ) THEN
         
            !skewSym = DCM - TRANSPOSE(DCM)
            !
            !logMap(1) = -skewSym(3,2)
            !logMap(2) =  skewSym(3,1)
            !logMap(3) = -skewSym(2,1)
            !
            !logMap = 0.5_DbKi * logMap   ! Eq. 26b with limit as x approaches 0 of (x/sin(x)) = 1
         
         
            logMap = 0.0_DbKi                                                   ! Eq. 26a
                  
         ELSE ! 0 < theta < pi 
      
            !skewSym = DCM - TRANSPOSE(DCM)
      
            logMap(1) = -DCM(3,2) + DCM(2,3) !-skewSym(3,2)
            logMap(2) =  DCM(3,1) - DCM(1,3) ! skewSym(3,1)
            logMap(3) = -DCM(2,1) + DCM(1,2) !-skewSym(2,1)
      
            logMap    = theta / TwoSinTheta * logMap   ! Eq. 26b
         END IF
         
      END IF

      
   END SUBROUTINE DCM_logMapD
!=======================================================================
!> \copydoc nwtc_num::dcm_logmapd
   SUBROUTINE DCM_logMapR(DCM, logMap, ErrStat, ErrMsg, thetaOut)
   
      ! This function computes the logarithmic map for a direction cosine matrix.
   
   REAL(SiKi),         INTENT(IN)    :: DCM(3,3)
   REAL(SiKi),         INTENT(  OUT) :: logMap(3)
   REAL(SiKi),OPTIONAL,INTENT(  OUT) :: thetaOut
   INTEGER(IntKi),     INTENT(  OUT) :: ErrStat                   ! Error status of the operation
   CHARACTER(*),       INTENT(  OUT) :: ErrMsg                    ! Error message if ErrStat /= ErrID_None
   
      ! local variables
   REAL(SiKi)                        :: cosTheta
   REAL(SiKi)                        :: theta
   REAL(SiKi)                        :: TwoSinTheta
   REAL(SiKi)                        :: v(3)
   REAL(SiKi)                        :: divisor
   INTEGER(IntKi)                    :: indx_max
      
         ! initialization
      ErrStat = ErrID_None
      ErrMsg  = ""   
   
   
      cosTheta  = 0.5_SiKi*( trace(DCM) - 1.0_SiKi )
      cosTheta  = min( max(cosTheta,-1.0_SiKi), 1.0_SiKi ) !make sure it's in a valid range (to avoid cases where this is slightly outside the +/-1 range)
      theta     = ACOS( cosTheta )                         ! Eq. 25 ( 0<=theta<=pi )
      
      
      !IF ( EqualRealNos( pi, theta )  ) THEN
      IF ( theta > 3.1_SiKi ) THEN  ! theta/(2*sin(theta)) blows up quickly as theta approaches pi, 
         ! so I'm putting a pretty large tolerance on pi here, and using a different equation to find the solution near pi

         logMap(1) = 1.0_ReKi + DCM(1,1) - DCM(2,2) - DCM(3,3);
         logMap(2) = 1.0_ReKi - DCM(1,1) + DCM(2,2) - DCM(3,3);
         logMap(3) = 1.0_ReKi - DCM(1,1) - DCM(2,2) + DCM(3,3);
             
         indx_max = maxloc( abs(logMap), 1 )
             
         divisor = sqrt(abs( logMap(indx_max) *  2.0_SiKi*(1.0_SiKi - cosTheta)  )) / theta  ! 2*(1-cosTheta)/theta^2 * abs(lambda(indx_max))
         if (indx_max == 1) then
           !logMap(1) = 1.0 + DCM(1,1) - DCM(2,2) - DCM(3,3)               ! 2*(1-cosTheta)/theta^2 * lambda(1) * lambda(1)
            logMap(2) = DCM(1,2) + DCM(2,1)                                ! 2*(1-cosTheta)/theta^2 * lambda(1) * lambda(2)
            logMap(3) = DCM(1,3) + DCM(3,1)                                ! 2*(1-cosTheta)/theta^2 * lambda(1) * lambda(3)
         elseif (indx_max == 2) then
            logMap(1) = DCM(1,2) + DCM(2,1)                                ! 2*(1-cosTheta)/theta^2 * lambda(2) * lambda(1)
           !logMap(2) = 1.0 - DCM(1,1) + DCM(2,2) - DCM(3,3)               ! 2*(1-cosTheta)/theta^2 * lambda(2) * lambda(2)
            logMap(3) = DCM(2,3) + DCM(3,2)                                ! 2*(1-cosTheta)/theta^2 * lambda(2) * lambda(3)
         else
            logMap(1) = DCM(1,3) + DCM(3,1)                                ! 2*(1-cosTheta)/theta^2 * lambda(3) * lambda(1)
            logMap(2) = DCM(2,3) + DCM(3,2)                                ! 2*(1-cosTheta)/theta^2 * lambda(3) * lambda(2)
           !logMap(3) = 1.0 - DCM(1,1) - DCM(2,2) + DCM(3,3)               ! 2*(1-cosTheta)/theta^2 * lambda(3) * lambda(3)
         end if
         logMap = logMap / divisor                                         ! lambda * sign(lambda(indx))
      
         ! at this point we may have the wrong sign for logMap (though if theta==pi, it doesn't matter because we can change it in the DCM_setLogMapforInterp() routines)
         ! we'll do a little checking to see if we should change the sign:
         
         IF ( EqualRealNos( Pi_S, theta )  ) RETURN
         
         v(1) = -DCM(3,2) + DCM(2,3) !-skewSym(3,2)
         v(2) =  DCM(3,1) - DCM(1,3) ! skewSym(3,1)
         v(3) = -DCM(2,1) + DCM(1,2) !-skewSym(2,1)
 
         indx_max = maxloc( abs(v), 1 )  ! find component with largest magnitude
         if ( .not. EqualRealNos( sign(1.0_SiKi,v(indx_max)), sign(1.0_SiKi,logMap(indx_max)) )) logMap = -logMap
         
      ELSE
         
         TwoSinTheta = 2.0_SiKi*sin(theta)
         
         IF ( EqualRealNos(0.0_SiKi, theta) .or. EqualRealNos( 0.0_SiKi, TwoSinTheta ) ) THEN
         
            !skewSym = DCM - TRANSPOSE(DCM)
            !
            !logMap(1) = -skewSym(3,2)
            !logMap(2) =  skewSym(3,1)
            !logMap(3) = -skewSym(2,1)
            !
            !logMap = 0.5_ReKi * logMap   ! Eq. 26b with limit as x approaches 0 of (x/sin(x)) = 1
         
         
            logMap = 0.0_SiKi                                                   ! Eq. 26a
                  
         ELSE ! 0 < theta < pi 
      
            logMap(1) = -DCM(3,2) + DCM(2,3) !-skewSym(3,2)
            logMap(2) =  DCM(3,1) - DCM(1,3) ! skewSym(3,1)
            logMap(3) = -DCM(2,1) + DCM(1,2) !-skewSym(2,1)
      
            logMap    = theta / TwoSinTheta * logMap   ! Eq. 26b
         END IF
         
      END IF
      
      IF ( PRESENT( thetaOut ) ) THEN
         thetaOut = theta
      END IF      
      
   END SUBROUTINE DCM_logMapR  
!=======================================================================  
!> This routine sets the rotation parameters (logMap tensors from dcm_logmap)
!! so that they can be appropriately interpolated, based on
!! continunity of the neighborhood. The tensor input matrix has columns
!! of rotational parameters; one column for each set of values to be 
!! interpolated (i.e., for each column, i, tensor(:,i) is the returned logMap value from the routine dcm_logmap).
!!
!! This is based on the \f$2\pi\f$ periodicity of rotations: \n
!! if \f$\lambda\f$ is one solution to \f$\log(\Lambda)\f$, then so is 
!! \f$\lambda_k = \lambda \left( 1 + \frac{2k\pi}{\left\| \lambda \right\|}\right)\f$ for any integer k.  
!! 
!! Use DCM_SetLogMapForInterp (nwtc_num::dcm_setlogmapforinterp) instead of directly calling a specific routine in the generic interface. 
   SUBROUTINE DCM_SetLogMapForInterpD( tensor )
         
   REAL(DbKi),     INTENT(INOUT) :: tensor(:,:)       !< a 3xn matrix, whose columns represent individual skew-symmmetric matrices. On exit,
                                                      !! each column will be within \f$2\pi\f$ of the previous column, allowing for interpolation 
                                                      !! of the quantities.

   REAL(DbKi)                    :: diff1, diff2      ! magnitude-squared of difference between two adjacent values
   REAL(DbKi)                    :: temp(3), temp1(3) ! difference between two tensors
   REAL(DbKi)                    :: period(3)         ! the period to add to the rotational parameters
   INTEGER(IntKi)                :: nc                ! size of the tensors matrix
   INTEGER(IntKi)                :: ic                ! loop counters for each array dimension
   
   nc = size(tensor,2)
          
      ! 
   do ic=2,nc      
      
      diff1 = TwoNorm( tensor(:,ic) )
      
      if ( .NOT. EqualRealNos( diff1, 0.0_DbKi) ) then
            ! check if we're going around a 2pi boundary:
      
         period = tensor(:,ic) * ( Twopi_D/diff1 )
      
         temp1 = tensor(:,ic-1) - tensor(:,ic)
         diff1 = DOT_PRODUCT( temp1, temp1 )
                            
            ! try for k < 0
         temp = temp1 + period !k=-1; 
         diff2 = DOT_PRODUCT( temp, temp )
      
         if (diff2 < diff1) then
         
            do while (diff2 < diff1)
               tensor(:,ic) = tensor(:,ic) - period  !k=k-1
                              
               diff1 = diff2
               temp  = temp + period !k=k-1; % = tensor(:,ic-1) - tensor(:,ic)
               diff2 = DOT_PRODUCT( temp, temp )
            end do
         
         else
            ! try for k > 0
         
               ! check if the new value is too small:
            temp = temp1 - period !k=+1; 
            diff2 = DOT_PRODUCT( temp, temp )
            
            do while (diff2 < diff1)
               tensor(:,ic) = tensor(:,ic) + period  !k=k+1

               diff1 = diff2
               temp  = temp - period !k=k+1; % = tensor(:,ic-1) - tensor(:,ic)
               diff2 = DOT_PRODUCT( temp, temp )
            end do
   
         end if
      
      end if ! tensor vector isn't zero=length
            
   end do
                 
   END SUBROUTINE DCM_SetLogMapForInterpD
!=======================================================================         
!> \copydoc nwtc_num::dcm_setlogmapforinterpd
   SUBROUTINE DCM_SetLogMapForInterpR( tensor )

   ! this routine sets the rotation parameters (tensors from DCM_logMap)
   ! so that they can be appropriately interpolated, based on
   ! continunity of the neighborhood. The tensor input matrix has columns
   ! of rotational parameters; one column for each set of values to be 
   ! interpolated
   !
   ! This is based on the 2pi periodicity of rotations:
   ! if tensor is one solution to DCM_logMap( DCM ), then so is
   !  tensor*( 1 + TwoPi*k/TwoNorm(tensor) ) for any integer k
      
   
   REAL(SiKi),     INTENT(INOUT) :: tensor(:,:)

   REAL(SiKi)                    :: diff1, diff2      ! magnitude-squared of difference between two adjacent values
   REAL(SiKi)                    :: temp(3), temp1(3) ! difference between two tensors
   REAL(SiKi)                    :: period(3)         ! the period to add to the rotational parameters
   INTEGER(IntKi)                :: nc                ! size of the tensors matrix
   INTEGER(IntKi)                :: ic                ! loop counters for each array dimension
   
   nc = size(tensor,2)
          
      ! 
   do ic=2,nc      
      
      diff1 = TwoNorm( tensor(:,ic) )
      
      if ( .NOT. EqualRealNos( diff1, 0.0_SiKi) ) then
            ! check if we're going around a 2pi boundary:
      
         period = tensor(:,ic) * ( Twopi/diff1 )
      
         temp1 = tensor(:,ic-1) - tensor(:,ic)
         diff1 = DOT_PRODUCT( temp1, temp1 )
                            
            ! try for k < 0
         temp = temp1 + period !k=-1; 
         diff2 = DOT_PRODUCT( temp, temp )
      
         if (diff2 < diff1) then
         
            do while (diff2 < diff1)
               tensor(:,ic) = tensor(:,ic) - period  !k=k-1
                              
               diff1 = diff2
               temp  = temp + period !k=k-1; % = tensor(:,ic-1) - tensor(:,ic)
               diff2 = DOT_PRODUCT( temp, temp )
            end do
         
         else
            ! try for k > 0
         
               ! check if the new value is too small:
            temp = temp1 - period !k=+1; 
            diff2 = DOT_PRODUCT( temp, temp )
            
            do while (diff2 < diff1)
               tensor(:,ic) = tensor(:,ic) + period  !k=k+1

               diff1 = diff2
               temp  = temp - period !k=k+1; % = tensor(:,ic-1) - tensor(:,ic)
               diff2 = DOT_PRODUCT( temp, temp )
            end do
   
         end if
      
      end if ! tensor vector isn't zero=length
            
   end do
                 
   END SUBROUTINE DCM_SetLogMapForInterpR
!=======================================================================     
!> This function compares two real numbers and determines if they
!! are "almost" equal, i.e. within some relative tolerance (basically ignoring the last 2 significant digits)
!! (see "Safe Comparisons" suggestion from http://www.lahey.com/float.htm)
!!
!! Note that the numbers are added together in this routine, so overflow can result if comparing two "huge" numbers. \n
!! Use EqualRealNos (nwtc_num::equalrealnos) instead of directly calling a specific routine in the generic interface. 
   FUNCTION EqualRealNos4 ( ReNum1, ReNum2 )

      ! passed variables

   REAL(SiKi), INTENT(IN )         :: ReNum1                            !< the first  real number to compare
   REAL(SiKi), INTENT(IN )         :: ReNum2                            !< the second real number to compare

   LOGICAL                         :: EqualRealNos4                     !< .true. if and only if the numbers are almost equal

      ! local variables
   REAL(SiKi), PARAMETER           :: Eps = EPSILON(ReNum1)             ! machine precision
   REAL(SiKi), PARAMETER           :: Tol = 100.0_SiKi*Eps / 2.0_SiKi   ! absolute tolerance (ignore the last 2 significant digits)

   REAL(SiKi)                      :: Fraction


      ! make sure we're never trying to get more precision than Tol

   Fraction = MAX( ABS(ReNum1+ReNum2), 1.0_SiKi )



      ! determine if ReNum1 and ReNum2 are approximately equal

   IF ( ABS(ReNum1 - ReNum2) <= Fraction*Tol ) THEN  ! the relative error
      EqualRealNos4 = .TRUE.
   ELSE
      EqualRealNos4 = .FALSE.
   ENDIF


   END FUNCTION EqualRealNos4
!=======================================================================
!> \copydoc nwtc_num::equalrealnos4
   FUNCTION EqualRealNos8 ( ReNum1, ReNum2 )

      ! passed variables

   REAL(R8Ki), INTENT(IN )         :: ReNum1                            ! the first  real number to compare
   REAL(R8Ki), INTENT(IN )         :: ReNum2                            ! the second real number to compare

   LOGICAL                         :: EqualRealNos8                     !< .true. if and only if the numbers are almost equal

      ! local variables
   REAL(R8Ki), PARAMETER           :: Eps = EPSILON(ReNum1)             ! machine precision
   REAL(R8Ki), PARAMETER           :: Tol = 100.0_R8Ki*Eps / 2.0_R8Ki   ! absolute tolerance (ignore the last 2 significant digits)

   REAL(R8Ki)                      :: Fraction


      ! make sure we're never trying to get more precision than Tol

   Fraction = MAX( ABS(ReNum1+ReNum2), 1.0_R8Ki )



      ! determine if ReNum1 and ReNum2 are approximately equal

   IF ( ABS(ReNum1 - ReNum2) <= Fraction*Tol ) THEN  ! the relative error
      EqualRealNos8 = .TRUE.
   ELSE
      EqualRealNos8 = .FALSE.
   ENDIF


   END FUNCTION EqualRealNos8
!=======================================================================
!> This function creates a rotation matrix, M, from a 1-2-3 rotation
!! sequence of the 3 Euler angles, \f$\theta_x\f$, \f$\theta_y\f$, and \f$\theta_z\f$, in radians.
!! M represents a change of basis (from global to local coordinates; 
!! not a physical rotation of the body). It is the inverse of EulerExtract (nwtc_num::eulerextract).
!!
!! \f{eqnarray*}{   
!! M & = & R(\theta_z) R(\theta_y) R(\theta_x) \\
!!   & = & \begin{bmatrix}  \cos(\theta_z) & \sin(\theta_z) & 0 \\
!!                         -\sin(\theta_z) & \cos(\theta_z) & 0 \\
!!                           0      &  0      & 1 \end{bmatrix}
!!         \begin{bmatrix}  \cos(\theta_y) & 0 & -\sin(\theta_y) \\
!!                                0 & 1 & 0        \\
!!                          \sin(\theta_y) & 0 & \cos(\theta_y)  \end{bmatrix}
!!         \begin{bmatrix}   1 &  0       & 0       \\
!!                           0 &  \cos(\theta_x) & \sin(\theta_x) \\
!!                           0 & -\sin(\theta_x) & \cos(\theta_x) \end{bmatrix} \\
!!   & = & \begin{bmatrix}  
!!    \cos(\theta_y)\cos(\theta_z) &   \cos(\theta_x)\sin(\theta_z)+\sin(\theta_x)\sin(\theta_y)\cos(\theta_z) &
!!                                     \sin(\theta_x)\sin(\theta_z)-\cos(\theta_x)\sin(\theta_y)\cos(\theta_z) \\
!!    -\cos(\theta_y)\sin(\theta_z)  & \cos(\theta_x)\cos(\theta_z)-\sin(\theta_x)\sin(\theta_y)\sin(\theta_z) & 
!!                                     \sin(\theta_x)\cos(\theta_z)+\cos(\theta_x)\sin(\theta_y)\sin(\theta_z) \\
!!    \sin(\theta_y)                & -\sin(\theta_x)\cos(\theta_y) & \cos(\theta_x)\cos(\theta_y) \\
!!         \end{bmatrix}   
!! \f}
!! Use EulerConstruct (nwtc_num::eulerconstruct) instead of directly calling a specific routine in the generic interface. 
   FUNCTION EulerConstructR4(theta) result(M)
   
   
      REAL(SiKi)             :: M(3,3)    !< rotation matrix, M 
      REAL(SiKi), INTENT(IN) :: theta(3)  !< the 3 rotation angles: \f$\theta_x, \theta_y, \theta_z\f$
      
      REAL(SiKi)             :: cx        ! cos(theta_x)
      REAL(SiKi)             :: sx        ! sin(theta_x)
      REAL(SiKi)             :: cy        ! cos(theta_y)
      REAL(SiKi)             :: sy        ! sin(theta_y)
      REAL(SiKi)             :: cz        ! cos(theta_z)
      REAL(SiKi)             :: sz        ! sin(theta_z)
   

      cx = cos( theta(1) )
      sx = sin( theta(1) )
      
      cy = cos( theta(2) )
      sy = sin( theta(2) )
      
      cz = cos( theta(3) )
      sz = sin( theta(3) )
         
      M(1,1) =  cy*cz            
      M(2,1) = -cy*sz            
      M(3,1) =  sy    
      
      M(1,2) =  cx*sz+sx*sy*cz            
      M(2,2) =  cx*cz-sx*sy*sz            
      M(3,2) =       -sx*cy     
      
      M(1,3) =  sx*sz-cx*sy*cz            
      M(2,3) =  sx*cz+cx*sy*sz            
      M(3,3) =        cx*cy               
   
   END FUNCTION EulerConstructR4
!=======================================================================
!> \copydoc nwtc_num::eulerconstructr4
   FUNCTION EulerConstructR8(theta) result(M)
   
      ! this function creates a rotation matrix, M, from a 1-2-3 rotation
      ! sequence of the 3 Euler angles, theta_x, theta_y, and theta_z, in radians.
      ! M represents a change of basis (from global to local coordinates; 
      ! not a physical rotation of the body). it is the inverse of EulerExtract (nwtc_num::eulerextract).
      !
      ! M = R(theta_z) * R(theta_y) * R(theta_x)
      !   = [ cz sz 0 |   [ cy  0 -sy |   [ 1   0   0 |
      !     |-sz cz 0 | * |  0  1   0 | * | 0  cx  sx |
      !     |  0  0 1 ]   | sy  0  cy ]   | 0 -sx  cx ]
      !   = [ cy*cz   cx*sz+sx*sy*cz    sx*sz-cx*sy*cz |
      !     |-cy*sz   cx*cz-sx*sy*sz    sx*cz+cx*sy*sz |
      !     | sy           -sx*cy             cx*cy    ]
      ! where cz = cos(theta_z), sz = sin(theta_z), cy = cos(theta_y), etc.
   
      REAL(R8Ki)             :: M(3,3)    ! rotation matrix M 
      REAL(R8Ki), INTENT(IN) :: theta(3)  ! the 3 rotation angles: theta_x, theta_y, theta_z
      
      REAL(R8Ki)             :: cx        ! cos(theta_x)
      REAL(R8Ki)             :: sx        ! sin(theta_x)
      REAL(R8Ki)             :: cy        ! cos(theta_y)
      REAL(R8Ki)             :: sy        ! sin(theta_y)
      REAL(R8Ki)             :: cz        ! cos(theta_z)
      REAL(R8Ki)             :: sz        ! sin(theta_z)
   

      cx = cos( theta(1) )
      sx = sin( theta(1) )
      
      cy = cos( theta(2) )
      sy = sin( theta(2) )
      
      cz = cos( theta(3) )
      sz = sin( theta(3) )
         
      M(1,1) =  cy*cz            
      M(2,1) = -cy*sz            
      M(3,1) =  sy    
      
      M(1,2) =  cx*sz+sx*sy*cz            
      M(2,2) =  cx*cz-sx*sy*sz            
      M(3,2) =       -sx*cy     
      
      M(1,3) =  sx*sz-cx*sy*cz            
      M(2,3) =  sx*cz+cx*sy*sz            
      M(3,3) =        cx*cy               
   
   END FUNCTION EulerConstructR8
!=======================================================================
!> if M is a rotation matrix from a 1-2-3 rotation sequence, this function returns 
!! the 3 Euler angles, \f$\theta_x\f$, \f$\theta_y\f$, and \f$\theta_z\f$ (in radians), that formed 
!! the matrix. M represents a change of basis (from global to local coordinates; 
!! not a physical rotation of the body). M is the inverse of EulerConstruct (nwtc_num::eulerconstruct).
!!
!! \f{eqnarray*}{   
!! M & = & R(\theta_z) R(\theta_y) R(\theta_x) \\
!!   & = & \begin{bmatrix}  \cos(\theta_z) & \sin(\theta_z) & 0 \\
!!                         -\sin(\theta_z) & \cos(\theta_z) & 0 \\
!!                           0      &  0      & 1 \end{bmatrix}
!!         \begin{bmatrix}  \cos(\theta_y) & 0 & -\sin(\theta_y) \\
!!                                0 & 1 & 0        \\
!!                          \sin(\theta_y) & 0 & \cos(\theta_y)  \end{bmatrix}
!!         \begin{bmatrix}   1 &  0       & 0       \\
!!                           0 &  \cos(\theta_x) & \sin(\theta_x) \\
!!                           0 & -\sin(\theta_x) & \cos(\theta_x) \end{bmatrix} \\
!!   & = & \begin{bmatrix}  
!!    \cos(\theta_y)\cos(\theta_z) &   \cos(\theta_x)\sin(\theta_z)+\sin(\theta_x)\sin(\theta_y)\cos(\theta_z) &
!!                                     \sin(\theta_x)\sin(\theta_z)-\cos(\theta_x)\sin(\theta_y)\cos(\theta_z) \\
!!    -\cos(\theta_y)\sin(\theta_z)  & \cos(\theta_x)\cos(\theta_z)-\sin(\theta_x)\sin(\theta_y)\sin(\theta_z) & 
!!                                     \sin(\theta_x)\cos(\theta_z)+\cos(\theta_x)\sin(\theta_y)\sin(\theta_z) \\
!!    \sin(\theta_y)                & -\sin(\theta_x)\cos(\theta_y) & \cos(\theta_x)\cos(\theta_y) \\
!!         \end{bmatrix}   
!! \f}
!! returned angles are in the range \f$\theta_x,\theta_y, \theta_z \in \left[ \pi, -\pi \right]\f$ \n
!! Use EulerExtract (nwtc_num::eulerextract)  instead of directly calling a specific routine in the generic interface. 
   FUNCTION EulerExtractR4(M) result(theta)
   
   
      REAL(SiKi), INTENT(IN) :: M(3,3)    !< rotation matrix, M 
      REAL(SiKi)             :: theta(3)  !< the 3 rotation angles: \f$\theta_x, \theta_y, \theta_z\f$
      
      REAL(SiKi)             :: cx        ! cos(theta_x)
      REAL(SiKi)             :: sx        ! sin(theta_x)
      REAL(SiKi)             :: cy        ! cos(theta_y)
!     REAL(SiKi)             :: sy        ! sin(theta_y)
      REAL(SiKi)             :: cz        ! cos(theta_z)
      REAL(SiKi)             :: sz        ! sin(theta_z)
   
         ! use trig identity sz**2 + cz**2 = 1 to get abs(cy):
      cy = sqrt( m(1,1)**2 + m(2,1)**2 ) 
!      cy = sqrt( m(3,3)**2 + m(3,2)**2 ) 
            
      if ( EqualRealNos(cy,0.0_SiKi) ) then
      !if ( cy < 16*epsilon(0.0_ReKi) ) then
         
         theta(2) = atan2( m(3,1), cy )               ! theta_y
         
         ! cy = 0 -> sy = +/- 1
         ! M  = [  0   cx*sz+/-sx*cz    sx*sz-/+cx*cz |
         !      |  0   cx*cz-/+sx*sz    sx*cz+/-cx*sz |
         !      |+/-1        0                0       ]
         
         ! gimbal lock allows us to choose theta_z = 0
         theta(3) = 0.0_SiKi                          ! theta_z
         
         ! which reduces the matrix to 
         ! M  = [  0  +/-sx  -/+cx |
         !      |  0     cx     sx |
         !      |+/-1    0       0 ]
         
         theta(1) = atan2(  m(2,3), m(2,2) )          ! theta_x
         
      else
         ! atan2( cy*sz, cy*cz )
         theta(3) = atan2( -m(2,1), m(1,1) )          ! theta_z         
         cz       = cos( theta(3) )
         sz       = sin( theta(3) )

            ! get the appropriate sign for cy:
         if ( EqualRealNos(cz, 0.0_SiKi) ) then
            cy = sign( cy, -m(2,1)/sz )
            !cy = -m(2,1)/sz
         else
            cy = sign( cy, m(1,1)/cz )
            !cy = -m(1,1)/cz
         end if
         theta(2) = atan2( m(3,1), cy )               ! theta_y
         
        !theta(1) = atan2( -m(3,2), m(3,3) )          ! theta_x
         
         ! for numerical reasons, we're going to get theta_x using
         ! M' = (R(theta_z) * R(theta_y))^T * M = R(theta_x)
         !    = [ cy  0  sy |   [ cz -sz 0 |       [ 1   0   0 |
         !      |  0  1   0 | * | sz  cz 0 | * M = | 0  cx  sx |
         !      |-sy  0  cy ]   |  0   0 1 ]       | 0 -sx  cx ]
         !    = [ cy*cz  -cy*sz  sy |       [ 1   0   0 |
         !      |    sz      cz   0 | * M = | 0  cx  sx |
         !      |-sy*cz   sy*sz  cy ]       | 0 -sx  cx ]
         ! taking M'(2,2) and M'(2,3) , we get cx and sx:
         ! sz*m(1,2) + cz*m(2,2) = cx
         ! sz*m(1,3) + cz*m(2,3) = sx

         cz = cos( theta(3) )
         sz = sin( theta(3) )
         
         cx = sz*m(1,2) + cz*m(2,2)
         sx = sz*m(1,3) + cz*m(2,3)
         
         theta(1) = atan2( sx, cx )
         
      end if
            
      
   END FUNCTION EulerExtractR4
!=======================================================================
!> \copydoc nwtc_num::eulerextractr4 
   FUNCTION EulerExtractR8(M) result(theta)
   
      ! if M is a rotation matrix from a 1-2-3 rotation sequence, this function returns 
      ! the 3 Euler angles, theta_x, theta_y, and theta_z (in radians), that formed 
      ! the matrix. M represents a change of basis (from global to local coordinates; 
      ! not a physical rotation of the body). M is the inverse of EulerConstruct().
      !
      ! M = R(theta_z) * R(theta_y) * R(theta_x)
      !   = [ cz sz 0 |   [ cy  0 -sy |   [ 1   0   0 |
      !     |-sz cz 0 | * |  0  1   0 | * | 0  cx  sx |
      !     |  0  0 1 ]   | sy  0  cy ]   | 0 -sx  cx ]
      !   = [ cy*cz   cx*sz+sx*sy*cz    sx*sz-cx*sy*cz |
      !     |-cy*sz   cx*cz-sx*sy*sz    sx*cz+cx*sy*sz |
      !     | sy           -sx*cy             cx*cy    ]
      ! where cz = cos(theta_z), sz = sin(theta_z), cy = cos(theta_y), etc.
      ! 
      ! returned angles are in the range [-pi, pi]
   
      REAL(R8Ki), INTENT(IN) :: M(3,3)    ! rotation matrix M 
      REAL(R8Ki)             :: theta(3)  ! the 3 rotation angles: theta_x, theta_y, theta_z
      
      REAL(R8Ki)             :: cx        ! cos(theta_x)
      REAL(R8Ki)             :: sx        ! sin(theta_x)
      REAL(R8Ki)             :: cy        ! cos(theta_y)
!     REAL(R8Ki)             :: sy        ! sin(theta_y)
      REAL(R8Ki)             :: cz        ! cos(theta_z)
      REAL(R8Ki)             :: sz        ! sin(theta_z)
   
         ! use trig identity sz**2 + cz**2 = 1 to get abs(cy):
      cy = sqrt( m(1,1)**2 + m(2,1)**2 ) 
!      cy = sqrt( m(3,3)**2 + m(3,2)**2 ) 
            
      if ( EqualRealNos(cy,0.0_R8Ki) ) then
      !if ( cy < 16*epsilon(0.0_ReKi) ) then
         
         theta(2) = atan2( m(3,1), cy )               ! theta_y
         
         ! cy = 0 -> sy = +/- 1
         ! M  = [  0   cx*sz+/-sx*cz    sx*sz-/+cx*cz |
         !      |  0   cx*cz-/+sx*sz    sx*cz+/-cx*sz |
         !      |+/-1        0                0       ]
         
         ! gimbal lock allows us to choose theta_z = 0
         theta(3) = 0.0_R8Ki                          ! theta_z
         
         ! which reduces the matrix to 
         ! M  = [  0  +/-sx  -/+cx |
         !      |  0     cx     sx |
         !      |+/-1    0       0 ]
         
         theta(1) = atan2(  m(2,3), m(2,2) )          ! theta_x
         
      else
         ! atan2( cy*sz, cy*cz )
         theta(3) = atan2( -m(2,1), m(1,1) )          ! theta_z         
         cz       = cos( theta(3) )
         sz       = sin( theta(3) )

            ! get the appropriate sign for cy:
         if ( EqualRealNos(cz, 0.0_R8Ki) ) then
            cy = sign( cy, -m(2,1)/sz )
            !cy = -m(2,1)/sz
         else
            cy = sign( cy, m(1,1)/cz )
            !cy = -m(1,1)/cz
         end if
         theta(2) = atan2( m(3,1), cy )               ! theta_y
         
        !theta(1) = atan2( -m(3,2), m(3,3) )          ! theta_x
         
         ! for numerical reasons, we're going to get theta_x using
         ! M' = (R(theta_z) * R(theta_y))^T * M = R(theta_x)
         !    = [ cy  0  sy |   [ cz -sz 0 |       [ 1   0   0 |
         !      |  0  1   0 | * | sz  cz 0 | * M = | 0  cx  sx |
         !      |-sy  0  cy ]   |  0   0 1 ]       | 0 -sx  cx ]
         !    = [ cy*cz  -cy*sz  sy |       [ 1   0   0 |
         !      |    sz      cz   0 | * M = | 0  cx  sx |
         !      |-sy*cz   sy*sz  cy ]       | 0 -sx  cx ]
         ! taking M'(2,2) and M'(2,3) , we get cx and sx:
         ! sz*m(1,2) + cz*m(2,2) = cx
         ! sz*m(1,3) + cz*m(2,3) = sx

         cz = cos( theta(3) )
         sz = sin( theta(3) )
         
         cx = sz*m(1,2) + cz*m(2,2)
         sx = sz*m(1,3) + cz*m(2,3)
         
         theta(1) = atan2( sx, cx )
         
      end if
            
      
   END FUNCTION EulerExtractR8

!=======================================================================
!> 
   FUNCTION EulerConstructZYXR8(theta) result(M)
   
      ! this function creates a rotation matrix, M, from a 3-2-1 rotation
      ! sequence of the 3 Euler angles, theta_z, theta_y, and theta_x, in radians.
      ! M represents a change of basis (from global to local coordinates; 
      ! not a physical rotation of the body). 
      !
      REAL(R8Ki)             :: M(3,3)    ! rotation matrix M 
      REAL(R8Ki), INTENT(IN) :: theta(3)  ! the 3 rotation angles: theta_x, theta_y, theta_z
      
      REAL(R8Ki)             :: cx        ! cos(theta_x)
      REAL(R8Ki)             :: sx        ! sin(theta_x)
      REAL(R8Ki)             :: cy        ! cos(theta_y)
      REAL(R8Ki)             :: sy        ! sin(theta_y)
      REAL(R8Ki)             :: cz        ! cos(theta_z)
      REAL(R8Ki)             :: sz        ! sin(theta_z)
   

      cx = cos( theta(1) )
      sx = sin( theta(1) )
      
      cy = cos( theta(2) )
      sy = sin( theta(2) )
      
      cz = cos( theta(3) )
      sz = sin( theta(3) )
         
      M(1,1) =  cy*cz            
      M(2,1) =  sx*sy*cz - sz*cx
      M(3,1) =  sx*sz + sy*cx*cz
      
      M(1,2) =  sz*cy
      M(2,2) =  sx*sy*sz + cx*cz
      M(3,2) =  -sx*cz + sy*sz*cx
      
      M(1,3) =  -sy
      M(2,3) =  sx*cy
      M(3,3) =  cx*cy
   
   END FUNCTION EulerConstructZYXR8
!=======================================================================
!> This routine sets the matrices in the first two dimensions of A equal 
!! to the identity matrix (all zeros, with ones on the diagonal).
!! If the first two dimensions of A are not equal (i.e., matrix A(:,:,n)    
!! is non-square), this routine returns the pseudo-identity.  
!!
!! Use eye (nwtc_num::eye) instead of directly calling a specific routine in the generic interface. 
   SUBROUTINE Eye2( A, ErrStat, ErrMsg )


   REAL(SiKi),     INTENT(INOUT) :: A (:,:)                        !< Array to set to the identity matrix (nr,nc,n)
   INTEGER(IntKi), INTENT(OUT)   :: ErrStat                        !< Error level
   CHARACTER(*),   INTENT(OUT)   :: ErrMsg                         !< ErrMsg corresponding to ErrStat

      ! local variables
   INTEGER                       :: j                              ! loop counter
   INTEGER                       :: nr                             ! number of rows
   INTEGER                       :: nc                             ! number of columns


   nr = SIZE(A,1)
   nc = SIZE(A,2)

   IF (nr /= nc) THEN
      ErrStat = ErrID_Info
      ErrMsg  = 'NWTC Library, Eye(): Matrix is not square.'
   ELSE
      ErrStat = ErrID_None
      ErrMsg = ''
   END IF

      ! initialize to zero:
   A = 0.0_SiKi

      ! set the diagonals to one:
   DO j = 1, MIN(nr,nc) ! the diagonal of the matrix
      A(j,j) = 1.0_SiKi
   END DO

   END SUBROUTINE Eye2
!=======================================================================
!> \copydoc nwtc_num::eye2 
   SUBROUTINE Eye2D( A, ErrStat, ErrMsg )

   REAL(DbKi),     INTENT(INOUT) :: A (:,:)                        !< Array to set to the identity matrix (nr,nc,n)
   INTEGER(IntKi), INTENT(OUT)   :: ErrStat                        !< Error level
   CHARACTER(*),   INTENT(OUT)   :: ErrMsg                         !< ErrMsg corresponding to ErrStat

      ! local variables
   INTEGER                       :: j                              ! loop counter
   INTEGER                       :: nr                             ! number of rows
   INTEGER                       :: nc                             ! number of columns


   nr = SIZE(A,1)
   nc = SIZE(A,2)

   IF (nr /= nc) THEN
      ErrStat = ErrID_Info
      ErrMsg  = 'NWTC Library, Eye(): Matrix is not square.'
   ELSE
      ErrStat = ErrID_None
      ErrMsg = ''
   END IF

      ! initialize to zero:
   A = 0._DbKi

      ! set the diagonals to one:
   DO j = 1, MIN(nr,nc) ! the diagonal of the matrix
      A(j,j) = 1._DbKi
   END DO

   END SUBROUTINE Eye2D
!=======================================================================
!> \copybrief nwtc_num::eye2 
   SUBROUTINE Eye3( A, ErrStat, ErrMsg )

      ! This routine sets each of the n matries A(:,:,n) to the identity
      ! matrix (all zeros, with ones on the diagonal).
      ! Note that this also returns the "pseudo-identity" when A(:,:)
      ! is not square (i.e., nr/=nc).

   REAL(SiKi),     INTENT(INOUT) :: A (:,:,:)                      ! Array to set to the identity matrix (nr,nc,n)
   INTEGER(IntKi), INTENT(OUT)   :: ErrStat                        ! Error level
   CHARACTER(*),   INTENT(OUT)   :: ErrMsg                         ! ErrMsg corresponding to ErrStat

      ! local variables
   INTEGER                       :: i, j                           ! loop counters
   INTEGER                       :: nr                             ! number of rows
   INTEGER                       :: nc                             ! number of columns
   INTEGER                       :: n                              ! number of matricies


   nr = SIZE(A,1)
   nc = SIZE(A,2)
   n  = SIZE(A,3)

   IF (nr /= nc) THEN
      ErrStat = ErrID_Info
      ErrMsg  = 'NWTC Library, Eye(): Matrix is not square.'
   ELSE
      ErrStat = ErrID_None
      ErrMsg = ''
   END IF

      ! initialize to zero:
   A = 0.0_SiKi

      ! set the diagonals to one:
   DO i = 1, n ! loop through the matrices
      DO j = 1, MIN(nr,nc) ! the diagonal of the matrix
         A(j,j,i) = 1.0_SiKi
      END DO
   END DO

   END SUBROUTINE Eye3
!=======================================================================
!> \copybrief nwtc_num::eye2 
   SUBROUTINE Eye3D( A, ErrStat, ErrMsg )

      ! This routine sets each of the n matries A(:,:,n) to the identity
      ! matrix (all zeros, with ones on the diagonal).
      ! Note that this also returns the "pseudo-identity" when A(:,:)
      ! is not square (i.e., nr/=nc).

   REAL(DbKi),     INTENT(INOUT) :: A (:,:,:)                      !< Array to set to the identity matrix (nr,nc,n)
   INTEGER(IntKi), INTENT(OUT)   :: ErrStat                        !< Error level
   CHARACTER(*),   INTENT(OUT)   :: ErrMsg                         !< ErrMsg corresponding to ErrStat

      ! local variables
   INTEGER                       :: i, j                           ! loop counters
   INTEGER                       :: nr                             ! number of rows
   INTEGER                       :: nc                             ! number of columns
   INTEGER                       :: n                              ! number of matricies


   nr = SIZE(A,1)
   nc = SIZE(A,2)
   n  = SIZE(A,3)

   IF (nr /= nc) THEN
      ErrStat = ErrID_Info
      ErrMsg  = 'NWTC Library, Eye(): Matrix is not square.'
   ELSE
      ErrStat = ErrID_None
      ErrMsg = ''
   END IF

      ! initialize to zero:
   A = 0._ReKi

      ! set the diagonals to one:
   DO i = 1, n ! loop through the matrices
      DO j = 1, MIN(nr,nc) ! the diagonal of the matrix
         A(j,j,i) = 1._DbKi
      END DO
   END DO

   END SUBROUTINE Eye3D
!====================================================================================================
INTEGER FUNCTION FindValidChannelIndx(OutListVal, ValidParamAry, SignM_out) RESULT( Indx )

   CHARACTER(*),                INTENT(IN)  :: OutListVal
   CHARACTER(OutStrLenM1),      INTENT(IN)  :: ValidParamAry(:)
   INTEGER,           OPTIONAL, INTENT(OUT) :: SignM_out
   
   CHARACTER(ChanLen)             :: OutListTmp                                      ! A string to temporarily hold OutList(I)
   INTEGER                        :: SignM
   LOGICAL                        :: CheckOutListAgain                               ! Flag used to determine if output parameter starting with "M" is valid (or the negative of another parameter)
   
      OutListTmp          = OutListVal

      ! Reverse the sign (+/-) of the output channel if the user prefixed the
      !   channel name with a "-", "_", "m", or "M" character indicating "minus".
      CheckOutListAgain = .FALSE.

      IF      ( INDEX( "-_", OutListTmp(1:1) ) > 0 ) THEN
         SignM = -1                         ! ex, "-TipDxc1" causes the sign of TipDxc1 to be switched.
         OutListTmp          = OutListTmp(2:)
      ELSE IF ( INDEX( "mM", OutListTmp(1:1) ) > 0 ) THEN ! We'll assume this is a variable name for now, (if not, we will check later if OutListTmp(2:) is also a variable name)
         CheckOutListAgain   = .TRUE.
         SignM = 1
      ELSE
         SignM = 1
      END IF

      CALL Conv2UC( OutListTmp )    ! Convert OutListTmp to upper case


      Indx = IndexCharAry( OutListTmp(1:OutStrLenM1), ValidParamAry )


         ! If it started with an "M" (CheckOutListAgain) we didn't find the value in our list (Indx < 1)

      IF ( CheckOutListAgain .AND. Indx < 1 ) THEN    ! Let's assume that "M" really meant "minus" and then test again
         SignM         = -1                     ! ex, "MTipDxc1" causes the sign of TipDxc1 to be switched.
         OutListTmp    = OutListTmp(2:)

         Indx = IndexCharAry( OutListTmp(1:OutStrLenM1), ValidParamAry )
      END IF
      
      IF (PRESENT(SignM_out))  SignM_out = SignM
      
END FUNCTION FindValidChannelIndx
!=======================================================================
!> This routine uses the Gauss-Jordan elimination method for the
!!   solution of a given set of simultaneous linear equations.
!! NOTE: this routine works if no pivot points are zero and you
!!   don't want the eschelon or reduced eschelon form of the
!!   augmented matrix.  The form of the original augmented matrix
!!   IS preserved in this call.
!! This routine was originally in FAST.f90.
!! When AugMatIn = [ A b ], this routine returns the solution
!! vector x to the equation Ax = b.
   SUBROUTINE GaussElim( AugMatIn, NumEq, x, ErrStat, ErrMsg )


   IMPLICIT                        NONE


      ! Passed variables:

   INTEGER(IntKi), INTENT(IN )  :: NumEq                                           !< Number of equations in augmented matrix

   REAL(ReKi),     INTENT(IN )  :: AugMatIn (NumEq, NumEq+1 )                      !< Augmented matrix passed into this subroutine ( AugMatIn = [ A b ]
   REAL(ReKi),     INTENT(OUT)  :: x (NumEq)                                       !< Solution vector

   INTEGER(IntKi), INTENT(OUT)  :: ErrStat                                         !< Error level
   CHARACTER(*),   INTENT(OUT)  :: ErrMsg                                          !< ErrMsg corresponding to ErrStat


      ! Local variables:

   REAL(ReKi)                   :: AugMat   (NumEq,( NumEq + 1 ))                  ! The augmented matrix [A b]

   INTEGER(IntKi)               :: I                                               ! Steps through columns
   INTEGER(IntKi)               :: J                                               ! Steps through rows
   INTEGER(IntKi)               :: L                                               ! Steps through rows
   INTEGER(IntKi)               :: NAug                                            ! Column dimension of augmented matrix


      ! Initialize variables:

   ErrStat = ErrID_None                ! No error has occurred
   NAug    = NumEq + 1                 ! The column dimension of the augmented matrix


      ! Create the augmented matrix, AugMat = [A b] (we make a copy so we don't overwrite the existing matrix):

   AugMat = AugMatIn



      ! Perform Gauss-Jordan elimination and store the solution vector
      !   in the last column of the augmented matrix:

   DO L = 1,NumEq             ! Loop through all rows

      IF ( EqualRealNos( AugMat(L,L), 0.0_ReKi ) ) THEN
         ErrStat = ErrID_Fatal
         ErrMsg  = 'Division by zero in NWTC Library subroutine GaussElim.'
         RETURN
      END IF

      DO I = ( L + 1 ), NAug  ! Loop through all columns above current row number

         AugMat(L,I) = AugMat(L,I) / AugMat(L,L)

         DO J = 1,NumEq       ! Loop through all rows except L
            IF ( J /= L )  AugMat(J,I) = AugMat(J,I) - ( AugMat(J,L)*AugMat(L,I) )
         ENDDO                ! J - All rows except L

      ENDDO                   ! I - All columns above current row number

   ENDDO                      ! L - All rows


      ! Transfer the solution vector from AugMat() to x():

   x = AugMat(:,NAug)



   RETURN

   END SUBROUTINE GaussElim
!=======================================================================
!> Determine index of the point in Ary just below Val and the fractional distance to the next point in the array.
!! The elements of the array are assumed to be regularly spaced.
   SUBROUTINE GetOffsetReg ( Ary, NumPts, Val, Ind, Fract, ErrStat, ErrMsg )

      ! Argument declarations:

   INTEGER, INTENT(IN)          :: NumPts                                     !< Length of the array.

   REAL(ReKi), INTENT(IN)       :: Ary  (NumPts)                              !< Input array of regularly spaced values.
   REAL(ReKi), INTENT(OUT)      :: Fract                                      !< The fractional distance of Val between the surrounding array elements.
   REAL(ReKi), INTENT(IN)       :: Val                                        !< The value we hope to bound in the array.

   INTEGER(IntKi), INTENT(OUT)  :: ErrStat                                    !< Error status.
   INTEGER(IntKi), INTENT(OUT)  :: Ind                                        !< The index of the point in Ary just below Val.

   CHARACTER(*),   INTENT(OUT)  :: ErrMsg                                     !< Error message.


      ! Local declarations.

   REAL(ReKi)                   :: Del                                        ! The distances between the regularly spaced points.

!   INTEGER(IntKi)               :: ErrStatLcL                                 ! Local error status.



      ! Check the validity of the data.

   IF ( NumPts == 0 )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg = 'GetOffsetReg:The value of NumPts cannot be zero.'
      RETURN
   ELSE
      ErrStat = ErrID_None
      ErrMsg  = ""
   END IF


      ! Compute the distance between Ary values.

   Del = ( Ary(NumPts) - Ary(1) )/REAL( NumPts-1, ReKi )                    ! QUESTION: Is this more accurate than computing the distance between two adjacent points?


      ! Find the index of the array element just below Val.

   IF ( Val <= Ary(1) )  THEN
      Ind   = 1
      Fract = 0.0_ReKi
      RETURN
   ELSEIF ( Val >= Ary(NumPts) )  THEN
      Ind   = NumPts
      Fract = 0.0_ReKi
      RETURN
   ENDIF ! ( X <= XAry(1) )

   Ind   = INT( ( Val - Ary(1) )/Del ) + 1
   Fract = ( Val - Ary(Ind) )/Del

   RETURN

   END SUBROUTINE GetOffsetReg ! ( Ary, NumPts, Val, Ind, Fract, ErrStat, ErrMsg )
!=======================================================================
!   SUBROUTINE GetPermMat ( InpMat, PMat, ErrStat )
!
!      ! This subroutine computes a permutation matrix, PMat, for a given
!      ! input matrix, InpMat. It assumes that InpMat is of full rank
!      ! and for now, the matrices are 3 x 3.
!
!      ! passed variables
!
!   REAL(ReKi), INTENT(IN )         :: InpMat       (3,3)
!   REAL(ReKi), INTENT(OUT )        :: PMat         (3,3) !this could be integer, but we'll leave it real now
!   INTEGER,    INTENT(OUT )        :: ErrStat            ! a non-zero value indicates an error in the permutation matrix algorithm
!
!      ! local variables
!   INTEGER                         :: iCol               ! loop counter
!   INTEGER                         :: iRow               ! loop counter
!   INTEGER                         :: MaxCol             ! holds index of maximum value in a column
!
!   LOGICAL                         :: ChkCols     (3)    ! a check to make sure we have only one non-zero element per column
!
!      ! initialize some variables
!   PMat    = 0.0
!   ChkCols = .FALSE.
!   ErrStat = 0
!
!      ! find the pivots
!   DO iRow = 1,3
!
!      MaxCol = 1        ! initialize max index
!      DO iCol = 2,3
!         IF ( ABS(InpMat(iRow,iCol)) > ABS(InpMat(iRow,MaxCol)) ) &
!            MaxCol = iCol
!      END DO ! iCol
!
!      IF ( ChkCols(MaxCol) ) THEN   ! we can have only 1 non-zero entry per row and column, but we've just violated that!
!         CALL ProgAbort( ' Error in GetPermMat(): InpMat is not full rank.', TrapErrors = .TRUE. )
!         ErrStat = 1
!      END IF
!
!      PMat(MaxCol, iRow) = SIGN( 1.0_ReKi, InpMat(iRow,MaxCol) )  ! technically a permutation matrix would only have +1.0 (not -1.0)
!      ChkCols(MaxCol)    = .TRUE.
!
!   END DO ! iRow
!
!   RETURN
!   END SUBROUTINE GetPermMat ! ( InpMat, PMat, ErrStat )
!=======================================================================
!> This subroutine computes the angles that make up the input direction cosine matrix, DCMat,
!! assuming small angles. It is the inverse of SmllRotTrans (nwtc_num::smllrottrans). \n
!! Use GetSmllRotAngs (nwtc_num::getsmllrotangs) instead of directly calling a specific routine in the generic interface. 
!=======================================================================
   FUNCTION GetSmllRotAngsD ( DCMat, ErrStat, ErrMsg )
            
      ! passed variables

   REAL(DbKi), INTENT(IN )            :: DCMat          (3,3)  !< a direction cosine matrix
   INTEGER,    INTENT(OUT )           :: ErrStat               !< a non-zero value indicates an error in the permutation matrix algorithm
   CHARACTER(*),INTENT(OUT ),OPTIONAL :: ErrMsg                !< a non-zero value indicates an error in the permutation matrix algorithm

   REAL(DbKi)                         :: GetSmllRotAngsD ( 3 ) !< the rotational angles

      ! local variables
   REAL(DbKi)                         :: denom                 ! the denominator of the resulting matrix
   REAL(DbKi), PARAMETER              :: LrgAngle  = 0.4_DbKi  ! Threshold for when a small angle becomes large (about 23deg).  This comes from: COS(SmllAngle) ~ 1/SQRT( 1 + SmllAngle^2 ) and SIN(SmllAngle) ~ SmllAngle/SQRT( 1 + SmllAngle^2 ) results in ~5% error when SmllAngle = 0.4rad.



      ! initialize output angles (just in case there is an error that prevents them from getting set)

   GetSmllRotAngsD = 0.0_DbKi
   ErrStat         = ErrID_None
   ErrMsg          = ""

      ! calculate the small angles
   GetSmllRotAngsD(1) = DCMat(2,3) - DCMat(3,2)
   GetSmllRotAngsD(2) = DCMat(3,1) - DCMat(1,3)
   GetSmllRotAngsD(3) = DCMat(1,2) - DCMat(2,1)

   denom             = DCMat(1,1) + DCMat(2,2) + DCMat(3,3) - 1.0_DbKi

   IF ( .NOT. EqualRealNos( denom, 0.0_DbKi ) ) THEN
      GetSmllRotAngsD = GetSmllRotAngsD / denom

         ! check that the angles are, in fact, small
      IF ( ANY( ABS(GetSmllRotAngsD) > LrgAngle ) ) THEN
         ErrStat = ErrID_Severe

         IF (PRESENT(ErrMsg)) THEN
            ErrMsg = ' Angles in GetSmllRotAngs() are larger than '//TRIM(Num2LStr(LrgAngle))//' radians.'
         ELSE
            CALL ProgWarn( ' Angles in GetSmllRotAngs() are larger than '//TRIM(Num2LStr(LrgAngle))//' radians.' )
         END IF

      END IF

   ELSE
         ! check that the angles are, in fact, small (denom should be close to 2 if angles are small)
      ErrStat = ErrID_Fatal

      IF (PRESENT(ErrMsg)) THEN
         ErrMsg = ' Denominator is zero in GetSmllRotAngs().'
      ELSE
         CALL ProgAbort( ' Denominator is zero in GetSmllRotAngs().', TrapErrors = .TRUE. )
      END IF

   END IF


   END FUNCTION GetSmllRotAngsD
!=======================================================================
!> \copydoc nwtc_num::getsmllrotangsd 
   FUNCTION GetSmllRotAngsR ( DCMat, ErrStat, ErrMsg )

      ! passed variables

   REAL(SiKi), INTENT(IN )            :: DCMat          (3,3)
   INTEGER,    INTENT(OUT )           :: ErrStat               ! a non-zero value indicates an error in the permutation matrix algorithm
   CHARACTER(*),INTENT(OUT ),OPTIONAL :: ErrMsg                ! a non-zero value indicates an error in the permutation matrix algorithm

   REAL(SiKi)                         :: GetSmllRotAngsR ( 3 )

      ! local variables
   REAL(SiKi)                         :: denom                 ! the denominator of the resulting matrix
   REAL(SiKi), PARAMETER              :: LrgAngle  = 0.4_SiKi  ! Threshold for when a small angle becomes large (about 23deg).  This comes from: COS(SmllAngle) ~ 1/SQRT( 1 + SmllAngle^2 ) and SIN(SmllAngle) ~ SmllAngle/SQRT( 1 + SmllAngle^2 ) results in ~5% error when SmllAngle = 0.4rad.



      ! initialize output angles (just in case there is an error that prevents them from getting set)

   GetSmllRotAngsR = 0.0_SiKi
   ErrStat         = ErrID_None
   ErrMsg          = ""

      ! calculate the small angles
   GetSmllRotAngsR(1) = DCMat(2,3) - DCMat(3,2)
   GetSmllRotAngsR(2) = DCMat(3,1) - DCMat(1,3)
   GetSmllRotAngsR(3) = DCMat(1,2) - DCMat(2,1)

   denom             = DCMat(1,1) + DCMat(2,2) + DCMat(3,3) - 1.0_SiKi

   IF ( .NOT. EqualRealNos( denom, 0.0_SiKi ) ) THEN
      GetSmllRotAngsR = GetSmllRotAngsR / denom

         ! check that the angles are, in fact, small
      IF ( ANY( ABS(GetSmllRotAngsR) > LrgAngle ) ) THEN
         ErrStat = ErrID_Severe

         IF (PRESENT(ErrMsg)) THEN
            ErrMsg = ' Angles in GetSmllRotAngs() are larger than '//TRIM(Num2LStr(LrgAngle))//' radians.'
         ELSE
            CALL ProgWarn( ' Angles in GetSmllRotAngs() are larger than '//TRIM(Num2LStr(LrgAngle))//' radians.' )
         END IF

      END IF

   ELSE
         ! check that the angles are, in fact, small (denom should be close to 2 if angles are small)
      ErrStat = ErrID_Fatal

      IF (PRESENT(ErrMsg)) THEN
         ErrMsg = ' Denominator is zero in GetSmllRotAngs().'
      ELSE
         CALL ProgAbort( ' Denominator is zero in GetSmllRotAngs().', TrapErrors = .TRUE. )
      END IF

   END IF


   END FUNCTION GetSmllRotAngsR
!=======================================================================
!> This funtion returns the non-dimensional (-1:+1) location of the given Gauss-Legendre Quadrature point and its weight.
!! It works for NPts \f$\in \left[{1,6\right]\f$.
!! The values came from Carnahan, Brice; Luther, H.A.; Wilkes, James O.  (1969)  "Applied Numerical Methods."
   SUBROUTINE GL_Pts ( IPt, NPts, Loc, Wt, ErrStat, ErrMsg )

      ! Argument declarations.

   REAL(ReKi), INTENT(OUT)        :: Loc                                         !< The location of the specified point.
   REAL(ReKi), INTENT(OUT)        :: Wt                                          !< The weight for the specified point.

   INTEGER,     INTENT(OUT)       :: ErrStat                                     !< Error status
   CHARACTER(*),INTENT(OUT)       :: ErrMsg                                      !< Error message
   INTEGER, INTENT(IN   )         :: IPt                                         !< The quadrature point in question.
   INTEGER, INTENT(IN   )         :: NPts                                        !< The number of points used in the quadrature.


   ErrStat = ErrID_None
   ErrMsg  = ''


      ! Check to see if the number of points and the specific point are valid values.

   IF ( ( NPts < 1 ) .OR. ( NPts > 6 ) )  THEN
      ErrMsg = 'In function GL_Loc, the number of points used for Gauss-Legendre Quadrature must be between 1 and 6' &
                    //' (inclusive).  Instead, it is "'//TRIM( Int2LStr( NPts ) )//'".'
      ErrStat = ErrID_Fatal
      RETURN
   END IF

   IF ( ( Ipt < 1 ) .OR. ( Ipt > NPts ) )  THEN
      ErrMsg = 'In function GL_Loc, the point being used for Gauss-Legendre Quadrature must be between 1 and ' &
                   //TRIM( Int2LStr( NPts ) )//' (inclusive).  Instead, it is "'//TRIM( Int2LStr( Ipt ) )//'".'
      ErrStat = ErrID_Fatal
      RETURN
   END IF


      ! Set the location and weight of the point.

   SELECT CASE ( NPts )
      CASE ( 1 )                         ! Case 1 is really just rectangular integration.
         Loc = 0.0
         Wt  = 2.0
      CASE ( 2 )
         SELECT CASE ( Ipt )
            CASE ( 1 )
               Loc = -0.5773503
               Wt  =  1.0
            CASE ( 2 )
               Loc = 0.5773503
               Wt  = 1.0
          END SELECT ! Ipt
      CASE ( 3 )
         SELECT CASE ( Ipt )
            CASE ( 1 )
               Loc = -0.7745967
               Wt  =  0.5555556
            CASE ( 2 )
               Loc =  0.0
               Wt  =  0.8888889
            CASE ( 3 )
               Loc =  0.7745967
               Wt  =  0.5555556
         END SELECT ! Ipt
      CASE ( 4 )
         SELECT CASE ( Ipt )
            CASE ( 1 )
               Loc = -0.8611363
               Wt  =  0.3478548
            CASE ( 2 )
               Loc = -0.3399810
               Wt  =  0.6521452
            CASE ( 3 )
               Loc =  0.3399810
               Wt  =  0.6521452
            CASE ( 4 )
               Loc =  0.8611363
               Wt  =  0.3478548
         END SELECT ! Ipt
      CASE ( 5 )
         SELECT CASE ( Ipt )
            CASE ( 1 )
               Loc = -0.9061798
               Wt  =  0.2369269
            CASE ( 2 )
               Loc = -0.5384693
               Wt  =  0.4786287
            CASE ( 3 )
               Loc =  0.0
               Wt  =  0.5688889
            CASE ( 4 )
               Loc =  0.5384693
               Wt  =  0.4786287
            CASE ( 5 )
               Loc =  0.9061798
               Wt  =  0.2369269
         END SELECT ! Ipt
      CASE ( 6 )
         SELECT CASE ( Ipt )
            CASE ( 1 )
               Loc = -0.9324695
               Wt  =  0.1713245
            CASE ( 2 )
               Loc = -0.6612094
               Wt  =  0.3607616
            CASE ( 3 )
               Loc = -0.2386192
               Wt  =  0.4679139
            CASE ( 4 )
               Loc =  0.2386192
               Wt  =  0.4679139
            CASE ( 5 )
               Loc =  0.6612094
               Wt  =  0.3607616
            CASE ( 6 )
               Loc =  0.9324695
               Wt  =  0.1713245
         END SELECT ! Ipt
   END SELECT ! Npts

   RETURN
   END SUBROUTINE GL_Pts ! ( IPt, NPts, Loc, Wt [, ErrStat] )
!=======================================================================
!> This funtion returns an integer index such that CAry(IndexCharAry) = CVal. If
!! no element in the array matches CVal, the value -1 is returned.  The routine
!! performs a binary search on the input array to determine if CVal is an
!! element of the array; thus, CAry must be sorted and stored in increasing
!! alphebetical (ASCII) order. The routine does not check that the array is
!! sorted.  The routine assumes that CVal is type CHARACTER and CAry
!! is an array of CHARACTERS.
   FUNCTION IndexCharAry( CVal, CAry )

      ! Function declaration.


   INTEGER                      :: IndexCharAry                                   !< integer index such that CAry(IndexCharAry) = CVal

      ! Argument declarations.

   CHARACTER(*), INTENT(IN)     :: CVal                                           !< String to find
   CHARACTER(*), INTENT(IN)     :: CAry(:)                                        !< Array of strings to search



      ! Local declarations.

   INTEGER                      :: IHi                                             ! The high index into the arrays.
   INTEGER                      :: IMid                                            ! The mid-point index between IHi and ILo.
   INTEGER                      :: ILo


      ! Initialize some variables

   ILo = 1
   IHi = SIZE(CAry)

   IF (     CVal == CAry(ILo) ) THEN
      IndexCharAry = ILo
   ELSEIF ( CVal == CAry(IHi) ) THEN
      IndexCharAry = IHi
   ELSE
      IndexCharAry = -1


         ! Let's search!

      DO WHILE ( IHi-ILo > 1 )

         IMid = ( IHi + ILo )/2

         IF( CVal > CAry(IMid) ) THEN
            ILo = IMid
         ELSEIF (CVal < CAry(IMid) ) THEN
            IHi = IMid
         ELSE !Found it
            IndexCharAry = IMid
            EXIT
         END IF

      END DO

   END IF


   RETURN

   END FUNCTION IndexCharAry
!=======================================================================
!> This funtion returns a y-value that corresponds to an input x-value by interpolating into the arrays.
!! It uses a binary interpolation scheme that takes about log(AryLen) / log(2) steps to converge.
!! It returns the first or last YAry() value if XVal is outside the limits of XAry(). 
!!
!! Use InterpBin (nwtc_num::interpbin) instead of directly calling a specific routine in the generic interface. 
   FUNCTION InterpBinComp( XVal, XAry, YAry, ILo, AryLen )

      ! Function declaration.


   COMPLEX(ReKi)                :: InterpBinComp                                   !< The interpolated value of Y at XVal


      ! Argument declarations.

   INTEGER, INTENT(IN)          :: AryLen                                          !< Length of the arrays.
   INTEGER, INTENT(INOUT)       :: ILo                                             !< The low index into the arrays.

   REAL(ReKi), INTENT(IN)       :: XAry    (AryLen)                                !< Array of X values to be interpolated.
   REAL(ReKi), INTENT(IN)       :: XVal                                            !< X value to be interpolated.

   COMPLEX(ReKi), INTENT(IN)    :: YAry    (AryLen)                                !< Array of Y values to be interpolated.


      ! Local declarations.

   INTEGER                      :: IHi                                             ! The high index into the arrays.
   INTEGER                      :: IMid                                            ! The mid-point index between IHi and ILo.



      ! Let's check the limits first.

   IF ( XVal <= XAry(1) )  THEN
      InterpBinComp = YAry(1)
      ILo           = 1
      RETURN
   ELSE IF ( XVal >= XAry(AryLen) )  THEN
      InterpBinComp = YAry(AryLen)
      ILo           = AryLen - 1
      RETURN
   END IF


      ! Let's interpolate!

   ILo  = 1
   IHi  = AryLen

   DO WHILE ( IHi-ILo > 1 )

      IMid = ( IHi + ILo )/2

      IF ( XVal >= XAry(IMid) ) THEN
         ILo = IMid
      ELSE
         IHi = IMid
      END IF

   END DO

   InterpBinComp = YAry(ILo) + ( YAry(IHi) - YAry(ILo) )*( XVal - XAry(ILo) )/( XAry(IHi) - XAry(ILo) )


   RETURN
   END FUNCTION InterpBinComp ! ( XVal, XAry, YAry, ILo, AryLen )
!=======================================================================
!> \copydoc nwtc_num::interpbincomp
   FUNCTION InterpBinReal( XVal, XAry, YAry, ILo, AryLen )

      ! Function declaration.


   REAL(ReKi)                   :: InterpBinReal                                   !< The interpolated value of Y at XVal


      ! Argument declarations.

   INTEGER, INTENT(IN)          :: AryLen                                          !< Length of the arrays.
   INTEGER, INTENT(INOUT)       :: ILo                                             !< The low index into the arrays.

   REAL(ReKi), INTENT(IN)       :: XAry    (AryLen)                                !< Array of X values to be interpolated.
   REAL(ReKi), INTENT(IN)       :: XVal                                            !< X value to be interpolated.
   REAL(ReKi), INTENT(IN)       :: YAry    (AryLen)                                !< Array of Y values to be interpolated.


      ! Local declarations.

   INTEGER                      :: IHi                                             ! The high index into the arrays.
   INTEGER                      :: IMid                                            ! The mid-point index between IHi and ILo.



      ! Let's check the limits first.

   IF ( XVal <= XAry(1) )  THEN
      InterpBinReal = YAry(1)
      ILo           = 1
      RETURN
   ELSE IF ( XVal >= XAry(AryLen) )  THEN
      InterpBinReal = YAry(AryLen)
      ILo           = AryLen - 1
      RETURN
   END IF


      ! Let's interpolate!

   ILo  = 1
   IHi  = AryLen

   DO WHILE ( IHi-ILo > 1 )

      IMid = ( IHi + ILo )/2

      IF ( XVal >= XAry(IMid) ) THEN
         ILo = IMid
      ELSE
         IHi = IMid
      END IF

   END DO

   InterpBinReal = YAry(ILo) + ( YAry(IHi) - YAry(ILo) )*( XVal - XAry(ILo) )/( XAry(IHi) - XAry(ILo) )


   RETURN
   END FUNCTION InterpBinReal ! ( XVal, XAry, YAry, ILo, AryLen )
!=======================================================================
!> This funtion returns a y-value that corresponds to an input x-value by interpolating into the arrays.
!! It uses the passed index as the starting point and does a stepwise interpolation from there. This is
!! especially useful when the calling routines save the value from the last time this routine was called
!! for a given case where XVal does not change much from call to call. When there is no correlation
!! from one interpolation to another, InterpBin() (nwtc_num::interpbin) may be a better choice.
!! It returns the first or last YAry() value if XVal is outside the limits of XAry().
!!
!! Use InterpStp (nwtc_num::interpstp) instead of directly calling a specific routine in the generic interface. 
   FUNCTION InterpStpComp4( XVal, XAry, YAry, Ind, AryLen )

      ! Function declaration.

   COMPLEX(SiKi)                :: InterpStpComp4                                  ! The interpolated value of Y at XVal


      ! Argument declarations.

   INTEGER, INTENT(IN)          :: AryLen                                          !< Length of the arrays.
   INTEGER, INTENT(INOUT)       :: Ind                                             !< Initial and final index into the arrays.

   REAL(SiKi), INTENT(IN)       :: XAry    (AryLen)                                !< Array of X values to be interpolated.
   REAL(SiKi), INTENT(IN)       :: XVal                                            !< X value to be interpolated.

   COMPLEX(SiKi), INTENT(IN)    :: YAry    (AryLen)                                !< Array of Y values to be interpolated.



      ! Let's check the limits first.

   IF ( XVal <= XAry(1) )  THEN
      InterpStpComp4 = YAry(1)
      Ind            = 1
      RETURN
   ELSE IF ( XVal >= XAry(AryLen) )  THEN
      InterpStpComp4 = YAry(AryLen)
      Ind            = MAX(AryLen - 1, 1)
      RETURN
   END IF


     ! Let's interpolate!

   Ind = MAX( MIN( Ind, AryLen-1 ), 1 )

   DO

      IF ( XVal < XAry(Ind) )  THEN

         Ind = Ind - 1

      ELSE IF ( XVal >= XAry(Ind+1) )  THEN

         Ind = Ind + 1

      ELSE

         InterpStpComp4 = ( YAry(Ind+1) - YAry(Ind) )*( XVal - XAry(Ind) )/( XAry(Ind+1) - XAry(Ind) ) + YAry(Ind)
         RETURN

      END IF

   END DO


   RETURN
   END FUNCTION InterpStpComp4
!=======================================================================
!> \copydoc nwtc_num::interpstpcomp4
   FUNCTION InterpStpComp8( XVal, XAry, YAry, Ind, AryLen )

      ! Function declaration.

   COMPLEX(R8Ki)                :: InterpStpComp8                                  !< The interpolated value of Y at XVal


      ! Argument declarations.

   INTEGER, INTENT(IN)          :: AryLen                                          !< Length of the arrays.
   INTEGER, INTENT(INOUT)       :: Ind                                             !< Initial and final index into the arrays.

   REAL(R8Ki), INTENT(IN)       :: XAry    (AryLen)                                !< Array of X values to be interpolated.
   REAL(R8Ki), INTENT(IN)       :: XVal                                            !< X value to be interpolated.

   COMPLEX(R8Ki), INTENT(IN)    :: YAry    (AryLen)                                !< Array of Y values to be interpolated.



      ! Let's check the limits first.

   IF ( XVal <= XAry(1) )  THEN
      InterpStpComp8 = YAry(1)
      Ind            = 1
      RETURN
   ELSE IF ( XVal >= XAry(AryLen) )  THEN
      InterpStpComp8 = YAry(AryLen)
      Ind            = MAX(AryLen - 1, 1)
      RETURN
   END IF


     ! Let's interpolate!

   Ind = MAX( MIN( Ind, AryLen-1 ), 1 )

   DO

      IF ( XVal < XAry(Ind) )  THEN

         Ind = Ind - 1

      ELSE IF ( XVal >= XAry(Ind+1) )  THEN

         Ind = Ind + 1

      ELSE

         InterpStpComp8 = ( YAry(Ind+1) - YAry(Ind) )*( XVal - XAry(Ind) )/( XAry(Ind+1) - XAry(Ind) ) + YAry(Ind)
         RETURN

      END IF

   END DO


   RETURN
   END FUNCTION InterpStpComp8

!=======================================================================
!> \copydoc nwtc_num::interpstpcomp4
   FUNCTION InterpStpReal4( XVal, XAry, YAry, Ind, AryLen )

      ! Function declaration.

   REAL(SiKi)                   :: InterpStpReal4                                  !< The interpolated value of Y at XVal


      ! Argument declarations.

   INTEGER, INTENT(IN)          :: AryLen                                          ! Length of the arrays.
   INTEGER, INTENT(INOUT)       :: Ind                                             ! Initial and final index into the arrays.

   REAL(SiKi), INTENT(IN)       :: XAry    (AryLen)                                ! Array of X values to be interpolated.
   REAL(SiKi), INTENT(IN)       :: XVal                                            ! X value to be interpolated.
   REAL(SiKi), INTENT(IN)       :: YAry    (AryLen)                                ! Array of Y values to be interpolated.



      ! Let's check the limits first.

   IF ( XVal <= XAry(1) )  THEN
      InterpStpReal4 = YAry(1)
      Ind            = 1
      RETURN
   ELSE IF ( XVal >= XAry(AryLen) )  THEN
      InterpStpReal4 = YAry(AryLen)
      Ind            = MAX(AryLen - 1, 1)
      RETURN
   END IF


     ! Let's interpolate!

   Ind = MAX( MIN( Ind, AryLen-1 ), 1 )

   DO

      IF ( XVal < XAry(Ind) )  THEN

         Ind = Ind - 1

      ELSE IF ( XVal >= XAry(Ind+1) )  THEN

         Ind = Ind + 1

      ELSE

         InterpStpReal4 = ( YAry(Ind+1) - YAry(Ind) )*( XVal - XAry(Ind) )/( XAry(Ind+1) - XAry(Ind) ) + YAry(Ind)
         RETURN

      END IF

   END DO


   RETURN
   END FUNCTION InterpStpReal4
!=======================================================================
!> \copydoc nwtc_num::interpstpcomp4
   FUNCTION InterpStpReal4_8( XVal, XAry, YAry, Ind, AryLen )

      ! Function declaration.

   REAL(R8Ki)                   :: InterpStpReal4_8                                !< The interpolated value of Y at XVal


      ! Argument declarations.

   INTEGER, INTENT(IN)          :: AryLen                                          ! Length of the arrays.
   INTEGER, INTENT(INOUT)       :: Ind                                             ! Initial and final index into the arrays.

   REAL(SiKi), INTENT(IN)       :: XAry    (AryLen)                                ! Array of X values to be interpolated.
   REAL(SiKi), INTENT(IN)       :: XVal                                            ! X value to be interpolated.
   REAL(R8Ki), INTENT(IN)       :: YAry    (AryLen)                                ! Array of Y values to be interpolated.



      ! Let's check the limits first.

   IF ( XVal <= XAry(1) )  THEN
      InterpStpReal4_8 = YAry(1)
      Ind            = 1
      RETURN
   ELSE IF ( XVal >= XAry(AryLen) )  THEN
      InterpStpReal4_8 = YAry(AryLen)
      Ind            = MAX(AryLen - 1, 1)
      RETURN
   END IF


     ! Let's interpolate!

   Ind = MAX( MIN( Ind, AryLen-1 ), 1 )

   DO

      IF ( XVal < XAry(Ind) )  THEN

         Ind = Ind - 1

      ELSE IF ( XVal >= XAry(Ind+1) )  THEN

         Ind = Ind + 1

      ELSE

         InterpStpReal4_8 = ( YAry(Ind+1) - YAry(Ind) )*( XVal - XAry(Ind) )/( XAry(Ind+1) - XAry(Ind) ) + YAry(Ind)
         RETURN

      END IF

   END DO


   RETURN
   END FUNCTION InterpStpReal4_8 
!=======================================================================
!> \copydoc nwtc_num::interpstpcomp4
   FUNCTION InterpStpReal8( XVal, XAry, YAry, Ind, AryLen )

      ! Function declaration.

   REAL(R8Ki)                   :: InterpStpReal8                                  !< The interpolated value of Y at XVal


      ! Argument declarations.

   INTEGER, INTENT(IN)          :: AryLen                                          ! Length of the arrays.
   INTEGER, INTENT(INOUT)       :: Ind                                             ! Initial and final index into the arrays.

   REAL(R8Ki), INTENT(IN)       :: XAry    (AryLen)                                ! Array of X values to be interpolated.
   REAL(R8Ki), INTENT(IN)       :: XVal                                            ! X value to be interpolated.
   REAL(R8Ki), INTENT(IN)       :: YAry    (AryLen)                                ! Array of Y values to be interpolated.



      ! Let's check the limits first.

   IF ( XVal <= XAry(1) )  THEN
      InterpStpReal8 = YAry(1)
      Ind            = 1
      RETURN
   ELSE IF ( XVal >= XAry(AryLen) )  THEN
      InterpStpReal8 = YAry(AryLen)
      Ind            = MAX(AryLen - 1, 1)
      RETURN
   END IF


     ! Let's interpolate!

   Ind = MAX( MIN( Ind, AryLen-1 ), 1 )

   DO

      IF ( XVal < XAry(Ind) )  THEN

         Ind = Ind - 1

      ELSE IF ( XVal >= XAry(Ind+1) )  THEN

         Ind = Ind + 1

      ELSE

         InterpStpReal8 = ( YAry(Ind+1) - YAry(Ind) )*( XVal - XAry(Ind) )/( XAry(Ind+1) - XAry(Ind) ) + YAry(Ind)
         RETURN

      END IF

   END DO


   RETURN
   END FUNCTION InterpStpReal8 

!=======================================================================
!> This funtion returns a y-value array that corresponds to an input x-value by interpolating into the arrays.
!! It uses the passed index as the starting point and does a stepwise interpolation from there. This is
!! especially useful when the calling routines save the value from the last time this routine was called
!! for a given case where XVal does not change much from call to call. 
!! It returns the first or last Y() row value if XVal is outside the limits of XAry().
   SUBROUTINE InterpStpMat( XVal, XAry, Y, Ind, AryLen, yInterp )

      ! Function declaration.

   REAL(ReKi), intent(out)      :: yInterp(:)                                      !< The interpolated value(s) of Y(dim=2) at XVal


      ! Argument declarations.

   INTEGER, INTENT(IN)          :: AryLen                                          !< Length of the arrays.
   INTEGER, INTENT(INOUT)       :: Ind                                             !< Initial and final index into the arrays.

   REAL(ReKi), INTENT(IN)       :: XAry    (AryLen)                                !< Array of X values to be interpolated.
   REAL(ReKi), INTENT(IN)       :: XVal                                            !< X value to be interpolated.
   REAL(ReKi), INTENT(IN)       :: Y       (:,:)                                   !< Matrix of Y values to be interpolated; First dimension is AryLen.



      ! Let's check the limits first.

   IF ( XVal <= XAry(1) )  THEN
      yInterp = Y(1,:)
      Ind     = 1
      RETURN
   ELSE IF ( XVal >= XAry(AryLen) )  THEN
      yInterp = Y(AryLen,:)
      Ind     = MAX(AryLen - 1, 1)
      RETURN
   END IF


     ! Let's interpolate!

   Ind = MAX( MIN( Ind, AryLen-1 ), 1 )

   DO

      IF ( XVal < XAry(Ind) )  THEN

         Ind = Ind - 1

      ELSE IF ( XVal >= XAry(Ind+1) )  THEN

         Ind = Ind + 1

      ELSE

         yInterp = ( Y(Ind+1,:) - Y(Ind,:) )*( XVal - XAry(Ind) )/( XAry(Ind+1) - XAry(Ind) ) + Y(Ind,:)
         RETURN

      END IF

   END DO


   RETURN
   END SUBROUTINE InterpStpMat
!=======================================================================   
!< This routine linearly interpolates Dataset. It is
!! set for a 2-d interpolation on x and y of the input point.
!! x and y must be in increasing order. Each dimension may contain only 1 value.
!! The method is described in this paper: 
!!   http://www.colorado.edu/engineering/CAS/courses.d/AFEM.d/AFEM.Ch11.d/AFEM.Ch11.pdf
   SUBROUTINE InterpStpReal2D( InCoord, Dataset, x, y, LastIndex, InterpData )

      INTEGER, PARAMETER :: NumDimensions = 2

         ! I/O variables

      REAL(ReKi),                     INTENT(IN   ) :: InCoord(NumDimensions)                       !< Arranged as (x, y)
      REAL(ReKi),                     INTENT(IN   ) :: Dataset(:,:)                                 !< Arranged as (x, y)
      REAL(ReKi),                     INTENT(IN   ) :: x(:)                                         !< first dimension in increasing order
      REAL(ReKi),                     INTENT(IN   ) :: y(:)                                         !< second dimension in increasing order
      INTEGER(IntKi),                 INTENT(INOUT) :: LastIndex(NumDimensions)                     !< Index for the last (x, y) used
      REAL(ReKi),                     INTENT(  OUT) :: InterpData                                   !< The interpolated value of Dataset(:,:) at InCoord


         ! Local variables

      INTEGER(IntKi)                                :: Indx_Lo(NumDimensions)                       ! index associated with lower bound of dimension 1,2 where val(Indx_lo(i)) <= InCoord(i) <= val(Indx_hi(i))
      INTEGER(IntKi)                                :: Indx_Hi(NumDimensions)                       ! index associated with upper bound of dimension 1,2 where val(Indx_lo(i)) <= InCoord(i) <= val(Indx_hi(i))
      REAL(ReKi)                                    :: Pos_Lo(NumDimensions)                        ! coordinate value with lower bound of dimension 1,2
      REAL(ReKi)                                    :: Pos_Hi(NumDimensions)                        ! coordinate value with upper bound of dimension 1,2

      REAL(ReKi)                                    :: isopc(NumDimensions)                         ! isoparametric coordinates

      REAL(ReKi)                                    :: N(2**NumDimensions)                          ! size 2^n
      REAL(ReKi)                                    :: u(2**NumDimensions)                          ! size 2^n

      INTEGER(IntKi)                                :: nx, ny


         ! find the indices into the arrays representing coordinates of each dimension:
         !  (by using LocateStp, we do not require equally spaced arrays)

      nx = SIZE(x)
      ny = SIZE(y)

      CALL LocateStp( InCoord(1), x, LastIndex(1), nx )
      CALL LocateStp( InCoord(2), y, LastIndex(2), ny )

      Indx_Lo = LastIndex  ! at this point, 0 <= Indx_Lo(i) <= n(i) for all i


      ! x (indx 1)
      IF (Indx_Lo(1) == 0) THEN
         Indx_Lo(1) = 1
      ELSEIF (Indx_Lo(1) == nx ) THEN
         Indx_Lo(1) = max( nx - 1, 1 )                ! make sure it's a valid index
      END IF
      Indx_Hi(1) = min( Indx_Lo(1) + 1 , nx )         ! make sure it's a valid index

      ! y (indx 2)
      IF (Indx_Lo(2) == 0) THEN
         Indx_Lo(2) = 1
      ELSEIF (Indx_Lo(2) == ny ) THEN
         Indx_Lo(2) = max( ny - 1, 1 )                ! make sure it's a valid index
      END IF
      Indx_Hi(2) = min( Indx_Lo(2) + 1 , ny )         ! make sure it's a valid index


         ! calculate the bounding box; the positions of all dimensions:

      pos_Lo(1) = x( Indx_Lo(1) )
      pos_Hi(1) = x( Indx_Hi(1) )

      pos_Lo(2) = y( Indx_Lo(2) )
      pos_Hi(2) = y( Indx_Hi(2) )


         ! 2-D linear interpolation:

      CALL IsoparametricCoords( InCoord, pos_Lo, pos_Hi, isopc )      ! Calculate iospc

      N(1)  = ( 1.0_ReKi + isopc(1) )*( 1.0_ReKi - isopc(2) )
      N(2)  = ( 1.0_ReKi + isopc(1) )*( 1.0_ReKi + isopc(2) )
      N(3)  = ( 1.0_ReKi - isopc(1) )*( 1.0_ReKi + isopc(2) )
      N(4)  = ( 1.0_ReKi - isopc(1) )*( 1.0_ReKi - isopc(2) )
      N     = N / REAL( SIZE(N), ReKi )  ! normalize


      u(1)  = Dataset( Indx_Hi(1), Indx_Lo(2) )
      u(2)  = Dataset( Indx_Hi(1), Indx_Hi(2) )
      u(3)  = Dataset( Indx_Lo(1), Indx_Hi(2) )
      u(4)  = Dataset( Indx_Lo(1), Indx_Lo(2) )

      InterpData = SUM ( N * u )


   END SUBROUTINE InterpStpReal2D   
!=======================================================================
!< This routine linearly interpolates Dataset. It is set for a 3-d 
!! interpolation on x and y of the input point. x, y, and z must be 
!! in increasing order. Each dimension may contain only 1 value.
!! The method is described in this paper: 
!!   http://www.colorado.edu/engineering/CAS/courses.d/AFEM.d/AFEM.Ch11.d/AFEM.Ch11.pdf
   SUBROUTINE InterpStpReal3D( InCoord, Dataset, x, y, z, LastIndex, InterpData )
   ! This routine linearly interpolates Dataset. It is set for a 3-d 
   ! interpolation on x and y of the input point. x, y, and z must be 
   ! in increasing order. Each dimension may contain only 1 value.
   ! The method is described in this paper: 
   !   http://www.colorado.edu/engineering/CAS/courses.d/AFEM.d/AFEM.Ch11.d/AFEM.Ch11.pdf

      INTEGER, PARAMETER :: NumDimensions = 3

         ! I/O variables

      REAL(ReKi),                     INTENT(IN   ) :: InCoord(NumDimensions)                       !< Arranged as (x, y, z)
      REAL(ReKi),                     INTENT(IN   ) :: Dataset(:,:,:)                               !< Arranged as (x, y, z)
      REAL(ReKi),                     INTENT(IN   ) :: x(:)                                         !< first dimension in increasing order
      REAL(ReKi),                     INTENT(IN   ) :: y(:)                                         !< second dimension in increasing order
      REAL(ReKi),                     INTENT(IN   ) :: z(:)                                         !< third dimension in increasing order
      INTEGER(IntKi),                 INTENT(INOUT) :: LastIndex(NumDimensions)                     !< Index for the last (x, y, z) used
      REAL(ReKi),                     INTENT(  OUT) :: InterpData                                   !< The interpolated value of Dataset(:,:,:) at InCoord


         ! Local variables

      INTEGER(IntKi)                                :: Indx_Lo(NumDimensions)                       ! index associated with lower bound of dimension i where val(Indx_lo(i)) <= InCoord(i) <= val(Indx_hi(i))
      INTEGER(IntKi)                                :: Indx_Hi(NumDimensions)                       ! index associated with upper bound of dimension i where val(Indx_lo(i)) <= InCoord(i) <= val(Indx_hi(i))
      REAL(ReKi)                                    :: Pos_Lo(NumDimensions)                        ! coordinate value with lower bound of dimension i
      REAL(ReKi)                                    :: Pos_Hi(NumDimensions)                        ! coordinate value with upper bound of dimension i

      REAL(ReKi)                                    :: isopc(NumDimensions)                         ! isoparametric coordinates

      REAL(ReKi)                                    :: N(2**NumDimensions)                          ! size 2^NumDimensions
      REAL(ReKi)                                    :: u(2**NumDimensions)                          ! size 2^NumDimensions

      INTEGER(IntKi)                                :: nd(NumDimensions)                            ! size of each dimension
      INTEGER(IntKi)                                :: i
   

         ! find the indices into the arrays representing coordinates of each dimension:
         !  (by using LocateStp, we do not require equally spaced frequencies or points)

      nd(1) = SIZE(x)
      nd(2) = SIZE(y)
      nd(3) = SIZE(z)

      CALL LocateStp( InCoord(1), x, LastIndex(1), nd(1) )
      CALL LocateStp( InCoord(2), y, LastIndex(2), nd(2) )
      CALL LocateStp( InCoord(3), z, LastIndex(3), nd(3) )

      Indx_Lo = LastIndex  ! at this point, 0 <= Indx_Lo(i) <= n(i) for all i


      DO i=1,NumDimensions
         IF (Indx_Lo(i) == 0) THEN
            Indx_Lo(i) = 1
         ELSEIF (Indx_Lo(i) == nd(i) ) THEN
            Indx_Lo(i) = max( nd(i) - 1, 1 )                ! make sure it's a valid index
         END IF
         Indx_Hi(i) = min( Indx_Lo(i) + 1 , nd(i) )         ! make sure it's a valid index
      END DO
   
 

         ! calculate the bounding box; the positions of all dimensions:

      pos_Lo(1) = x( Indx_Lo(1) )
      pos_Hi(1) = x( Indx_Hi(1) )

      pos_Lo(2) = y( Indx_Lo(2) )
      pos_Hi(2) = y( Indx_Hi(2) )

      pos_Lo(3) = z( Indx_Lo(3) )
      pos_Hi(3) = z( Indx_Hi(3) )
   

         ! 2-D linear interpolation:

      CALL IsoparametricCoords( InCoord, pos_Lo, pos_Hi, isopc )      ! Calculate iospc

   
      N(1)  = ( 1.0_ReKi + isopc(1) )*( 1.0_ReKi - isopc(2) )*( 1.0_ReKi - isopc(3) )
      N(2)  = ( 1.0_ReKi + isopc(1) )*( 1.0_ReKi + isopc(2) )*( 1.0_ReKi - isopc(3) )
      N(3)  = ( 1.0_ReKi - isopc(1) )*( 1.0_ReKi + isopc(2) )*( 1.0_ReKi - isopc(3) )
      N(4)  = ( 1.0_ReKi - isopc(1) )*( 1.0_ReKi - isopc(2) )*( 1.0_ReKi - isopc(3) )
      N(5)  = ( 1.0_ReKi + isopc(1) )*( 1.0_ReKi - isopc(2) )*( 1.0_ReKi + isopc(3) )
      N(6)  = ( 1.0_ReKi + isopc(1) )*( 1.0_ReKi + isopc(2) )*( 1.0_ReKi + isopc(3) )
      N(7)  = ( 1.0_ReKi - isopc(1) )*( 1.0_ReKi + isopc(2) )*( 1.0_ReKi + isopc(3) )
      N(8)  = ( 1.0_ReKi - isopc(1) )*( 1.0_ReKi - isopc(2) )*( 1.0_ReKi + isopc(3) )
      N     = N / REAL( SIZE(N), ReKi )  ! normalize
      
      u(1)  = Dataset( Indx_Hi(1), Indx_Lo(2), Indx_Lo(3) )
      u(2)  = Dataset( Indx_Hi(1), Indx_Hi(2), Indx_Lo(3) )
      u(3)  = Dataset( Indx_Lo(1), Indx_Hi(2), Indx_Lo(3) )
      u(4)  = Dataset( Indx_Lo(1), Indx_Lo(2), Indx_Lo(3) )
      u(5)  = Dataset( Indx_Hi(1), Indx_Lo(2), Indx_Hi(3) )
      u(6)  = Dataset( Indx_Hi(1), Indx_Hi(2), Indx_Hi(3) )
      u(7)  = Dataset( Indx_Lo(1), Indx_Hi(2), Indx_Hi(3) )
      u(8)  = Dataset( Indx_Lo(1), Indx_Lo(2), Indx_Hi(3) )   
   
      InterpData = SUM ( N * u )     ! could use dot_product, though I'm not sure it's the came for complex numbers
      

   END SUBROUTINE InterpStpReal3D   
!=======================================================================
!> This funtion returns a y-value that corresponds to an input x-value which is wrapped back
!! into the range [0-XAry(AryLen)] by interpolating into the arrays.  
!! It is assumed that XAry is sorted in ascending order.
!! It uses the passed index as the starting point and does a stepwise interpolation from there.  This is
!! especially useful when the calling routines save the value from the last time this routine was called
!! for a given case where XVal does not change much from call to call.  When there is no correlation
!! from one interpolation to another, InterpBin() may be a better choice.
!! It returns the first or last YAry() value if XVal is outside the limits of XAry().
!!
!! Use InterpWrappedStpReal (nwtc_num::interpwrappedstpreal) instead of directly calling a specific routine in the generic interface. 
   FUNCTION InterpWrappedStpReal4( XValIn, XAry, YAry, Ind, AryLen )

      ! Function declaration.

   REAL(SiKi)                   :: InterpWrappedStpReal4                           !< The interpolated value of Y at XVal


      ! Argument declarations.

   INTEGER, INTENT(IN)          :: AryLen                                          !< Length of the arrays.
   INTEGER, INTENT(INOUT)       :: Ind                                             !< Initial and final index into the arrays.

   REAL(SiKi), INTENT(IN)       :: XAry    (AryLen)                                !< Array of X values to be interpolated.
   REAL(SiKi), INTENT(IN)       :: XValIn                                          !< X value to be interpolated.
   REAL(SiKi), INTENT(IN)       :: YAry    (AryLen)                                !< Array of Y values to be interpolated.

   REAL(SiKi)                   :: XVal                                            !< X value to be interpolated.
   
   
   
      ! Wrap XValIn into the range XAry(1) to XAry(AryLen)
   XVal = MOD(XValIn, XAry(AryLen))

      ! Set the Ind to the first index if we are at the beginning of XAry
   IF ( XVal <= XAry(2) )  THEN  
      Ind           = 1
   END IF
   
   InterpWrappedStpReal4 = InterpStp( XVal, XAry, YAry, Ind, AryLen )
   
   
   END FUNCTION InterpWrappedStpReal4
!=======================================================================
!> \copydoc nwtc_num::interpwrappedstpreal4
   FUNCTION InterpWrappedStpReal4_8( XValIn, XAry, YAry, Ind, AryLen )

      ! Function declaration.

   REAL(R8Ki)                   :: InterpWrappedStpReal4_8                         !< The interpolated value of Y at XVal


      ! Argument declarations.

   INTEGER, INTENT(IN)          :: AryLen                                          ! Length of the arrays.
   INTEGER, INTENT(INOUT)       :: Ind                                             ! Initial and final index into the arrays.

   REAL(SiKi), INTENT(IN)       :: XAry    (AryLen)                                ! Array of X values to be interpolated.
   REAL(SiKi), INTENT(IN)       :: XValIn                                          ! X value to be interpolated.
   REAL(R8Ki), INTENT(IN)       :: YAry    (AryLen)                                ! Array of Y values to be interpolated.

   REAL(SiKi)                   :: XVal                                            ! X value to be interpolated.
   
   
   
      ! Wrap XValIn into the range XAry(1) to XAry(AryLen)
   XVal = MOD(XValIn, XAry(AryLen))

      ! Set the Ind to the first index if we are at the beginning of XAry
   IF ( XVal <= XAry(2) )  THEN  
      Ind           = 1
   END IF
   
   InterpWrappedStpReal4_8 = InterpStp( XVal, XAry, YAry, Ind, AryLen )
   
   
   END FUNCTION InterpWrappedStpReal4_8 
!=======================================================================
!> \copydoc nwtc_num::interpwrappedstpreal4
   FUNCTION InterpWrappedStpReal8( XValIn, XAry, YAry, Ind, AryLen )

      ! Function declaration.

   REAL(R8Ki)                   :: InterpWrappedStpReal8                           !< The interpolated value of Y at XVal


      ! Argument declarations.

   INTEGER, INTENT(IN)          :: AryLen                                          ! Length of the arrays.
   INTEGER, INTENT(INOUT)       :: Ind                                             ! Initial and final index into the arrays.

   REAL(R8Ki), INTENT(IN)       :: XAry    (AryLen)                                ! Array of X values to be interpolated.
   REAL(R8Ki), INTENT(IN)       :: XValIn                                           ! X value to be interpolated.
   REAL(R8Ki), INTENT(IN)       :: YAry    (AryLen)                                ! Array of Y values to be interpolated.

   REAL(R8Ki)                   :: XVal                                           ! X value to be interpolated.
   
   
   
      ! Wrap XValIn into the range XAry(1) to XAry(AryLen)
   XVal = MOD(XValIn, XAry(AryLen))

      ! Set the Ind to the first index if we are at the beginning of XAry
   IF ( XVal <= XAry(2) )  THEN  
      Ind           = 1
   END IF
   
   InterpWrappedStpReal8 = InterpStp( XVal, XAry, YAry, Ind, AryLen )
   
   
   END FUNCTION InterpWrappedStpReal8
!=======================================================================
!> This subroutine calculates interpolated values for an array of input values.
!! The size of the xknown and yknown arrays must match, and the size of the
!! xnew and ynew arrays must match. Xknown must be in ascending order.
!! Values outside the range of xknown are fixed to the end points.
   SUBROUTINE InterpArrayR4( xknown, yknown, xnew, ynew )
      REAL(SiKi), INTENT(IN   ) :: xknown(:)
      REAL(SiKi), INTENT(IN   ) :: yknown(:)
      REAL(SiKi), INTENT(IN   ) :: xnew(:)
      REAL(SiKi), INTENT(  OUT) :: ynew(:)
      integer(IntKi) i,itmp,nknown
      nknown=size(xknown)
      do i=1,size(xnew)
         itmp=minloc(abs(xnew(i)-xknown),dim=1)
         if (itmp==nknown) then
            if (xknown(itmp)>xnew(i)) then
               ynew(i)=interp_lin0(xnew(i),xknown(itmp-1),xknown(itmp),yknown(itmp-1),yknown(itmp))
            else
               ! The current x is above the max of xknown
               ! extrapolation required, here fixed to upper bound
               ynew(i)=yknown(nknown)
            endif
         elseif (xknown(itmp)<xnew(i)) then
            ! normal case, x between itmp and itmp+1
            ynew(i)=interp_lin0(xnew(i),xknown(itmp),xknown(itmp+1),yknown(itmp),yknown(itmp+1))
         elseif (itmp==1) then
            ! The current x is below the min of xknown
            ynew(i)=yknown(1)
         else
            ! normal case but inverted, x between itmp-1 and itmp
            ynew(i)=interp_lin0(xnew(i),xknown(itmp-1),xknown(itmp),yknown(itmp-1),yknown(itmp))
         endif
      enddo
      CONTAINS
         function interp_lin0(x,x0,x1,f0,f1)   ! Linear interpolation function                                     
            real(SiKi) ::interp_lin0
            real(SiKi),intent(in):: x,x0,x1,f0,f1
            if (EqualRealNos(x0,x1)) then    ! to avoid division by zero
               interp_lin0=f0
            else
               interp_lin0=(x-x1)/(x0-x1)*f0+(x-x0)/(x1-x0)*f1
            endif
         end function interp_lin0
   END SUBROUTINE InterpArrayR4
!=======================================================================
!> \copydoc nwtc_num::interparrayr4
   SUBROUTINE InterpArrayR8( xknown, yknown, xnew, ynew )
      REAL(R8Ki), INTENT(IN   ) :: xknown(:)
      REAL(R8Ki), INTENT(IN   ) :: yknown(:)
      REAL(R8Ki), INTENT(IN   ) :: xnew(:)
      REAL(R8Ki), INTENT(  OUT) :: ynew(:)
      integer(IntKi) i,itmp,nknown
      nknown=size(xknown)
      do i=1,size(xnew)
         itmp=minloc(abs(xnew(i)-xknown),dim=1)
         if (itmp==nknown) then
            if (xknown(itmp)>xnew(i)) then
               ynew(i)=interp_lin0(xnew(i),xknown(itmp-1),xknown(itmp),yknown(itmp-1),yknown(itmp))
            else
               ! The current x is above the max of xknown
               ! extrapolation required, here fixed to upper bound
               ynew(i)=yknown(nknown)
            endif
         elseif (xknown(itmp)<xnew(i)) then
            ! normal case, x between itmp and itmp+1
            ynew(i)=interp_lin0(xnew(i),xknown(itmp),xknown(itmp+1),yknown(itmp),yknown(itmp+1))
         elseif (itmp==1) then
            ! The current x is below the min of xknown
            ynew(i)=yknown(1)
         else
            ! normal case but inverted, x between itmp-1 and itmp
            ynew(i)=interp_lin0(xnew(i),xknown(itmp-1),xknown(itmp),yknown(itmp-1),yknown(itmp))
         endif
      enddo
      CONTAINS
         function interp_lin0(x,x0,x1,f0,f1)   ! Linear interpolation function                                     
            real(R8Ki) ::interp_lin0
            real(R8Ki),intent(in):: x,x0,x1,f0,f1
            if (EqualRealNos(x0,x1)) then    ! to avoid division by zero
               interp_lin0=f0
            else
               interp_lin0=(x-x1)/(x0-x1)*f0+(x-x0)/(x1-x0)*f1
            endif
         end function interp_lin0
   END SUBROUTINE InterpArrayR8
!=======================================================================
!> This subroutine calculates the iosparametric coordinates, isopc, which is a value between -1 and 1 
!! (for each dimension of a dataset), indicating where InCoord falls between posLo and posHi.
!! It is used in InterpStpReal2D (nwtcnum::interpstpreal2d) and InterpStpReal3D (nwtcnum::interpstpreal3d).
   SUBROUTINE IsoparametricCoords( InCoord, posLo, posHi, isopc )
   
      REAL(ReKi),     INTENT(IN   )          :: InCoord(:)                             !< Coordinate values we're interpolating to; (size = number of interpolation dimensions)
      REAL(ReKi),     INTENT(IN   )          :: posLo(:)                               !< coordinate values associated with Indx_Lo; (size = number of interpolation dimensions)
      REAL(ReKi),     INTENT(IN   )          :: posHi(:)                               !< coordinate values associated with Indx_Hi; (size = number of interpolation dimensions)
      REAL(ReKi),     INTENT(  OUT)          :: isopc(:)                               !< isoparametric coordinates; (position within the box)

      ! local variables
      REAL(ReKi)                             :: dx                                     ! difference between high and low coordinates in the bounding "box"
      INTEGER(IntKi)                         :: i                                      ! loop counter
   
   
      do i=1,size(isopc)
      
         dx = posHi(i) - posLo(i) 
         if (EqualRealNos(dx, 0.0_ReKi)) then
            isopc(i) = 1.0_ReKi
         else
            isopc(i) = ( 2.0_ReKi*InCoord(i) - posLo(i) - posHi(i) ) / dx
               ! to verify that we don't extrapolate, make sure this is bound between -1 and 1 (effectively nearest neighbor)
            isopc(i) = min( 1.0_ReKi, isopc(i) )
            isopc(i) = max(-1.0_ReKi, isopc(i) )
         end if
      
      end do
            
   END SUBROUTINE IsoparametricCoords   
!=======================================================================   
!> This function returns a logical TRUE/FALSE value that indicates
!! if the given (2-dimensional) matrix, A, is symmetric. If A is not
!! square it returns FALSE.
   FUNCTION IsSymmetric( A )

         ! passed variables

      REAL(ReKi), INTENT(IN) :: A(:,:)                   !< a real matrix A, whose symmetry is questioned
      LOGICAL                :: IsSymmetric              !< true if A is symmetric, false if not

         ! local variables

      INTEGER(IntKi)         :: i                        ! counter for rows
      INTEGER(IntKi)         :: j                        ! counter for columns
      INTEGER(IntKi)         :: N                        ! size of A


         ! If A is non-square, it is not symmetric:

      N = SIZE(A,1)

      IF ( N /= SIZE(A,2) ) THEN
         IsSymmetric = .FALSE.
         RETURN
      END IF


         ! If A(i,j) /= A(j,i), it is not symmetric:

      IsSymmetric = .TRUE.

      DO i = 1,(N-1)          ! Loop through the 1st N-1 rows of A
         DO j = (i+1),N       ! Loop through upper triangular part of A

            IsSymmetric = EqualRealNos( A(i,j), A(j,i) )
            IF ( .NOT. IsSymmetric ) RETURN

         END DO               ! j - All columns (rows) past I
      END DO                  ! i - The 1st N-1 rows (columns) of A


   END FUNCTION IsSymmetric
!=======================================================================
!> KERNELSMOOTHING Kernel smoothing of vector data
!!
!!   fNew = kernelSmoothing( x, f, KERNELTYPE, RADIUS ) generates a smoothed
!!   version of the data f(x) in fNew.  Supported KERNELTYPE values are
!!   'EPANECHINIKOV', 'QUARTIC' or 'BIWEIGHT', 'TRIWEIGHT', 'TRICUBE' and
!!   'GAUSSIAN'.  RADIUS controls the width of the kernel relative to the
!!   vector x.
!!
!!   See also: https://en.wikipedia.org/wiki/Kernel_(statistics)#Kernel_functions_in_common_use
subroutine kernelSmoothing(x, f, kernelType, radius, fNew)

   REAL(ReKi),             INTENT(in   ) :: x(:)         !> independent axis
   REAL(ReKi),             INTENT(in   ) :: f(:)         !> function values, f(x), to be smoothed
   INTEGER,                INTENT(in   ) :: kernelType   !> what kind of smoothing function to use
   REAL(ReKi),             INTENT(in   ) :: radius       !> width of the "window", in the units of x
   REAL(ReKi),             INTENT(  out) :: fNew(:)      !> smoothed function values
   
   REAL(ReKi)                            :: k
   REAL(ReKi)                            :: k_sum
   REAL(ReKi)                            :: w
   INTEGER(IntKi)                        :: Exp1
   INTEGER(IntKi)                        :: Exp2
   REAL(ReKi)                            :: u(size(x))
   INTEGER                               :: i, j
   INTEGER                               :: n
   
   ! check that radius > 0
   ! check that size(x) = size(f)=size(fNew)
   ! check that kernelType is a valid number
   
   n = size(x)
   
   
   ! make sure that the value of u is in [-1 and 1] for these kernels:
   if (kernelType /= kernelType_GAUSSIAN) then

      select case ( kernelType )
         case (kernelType_EPANECHINIKOV)
            w = 3.0_ReKi/4.0_ReKi
            Exp1 = 2
            Exp2 = 1
         case (kernelType_QUARTIC, kernelType_BIWEIGHT)
            w = 15.0_ReKi/16.0_ReKi
            Exp1 = 2
            Exp2 = 2
         case (kernelType_TRIWEIGHT)
            w = 35.0_ReKi/32.0_ReKi
            Exp1 = 2
            Exp2 = 3
         case (kernelType_TRICUBE)
            w = 70.0_ReKi/81.0_ReKi
            Exp1 = 3
            Exp2 = 3
      end select
         
      fNew = 0.0_ReKi ! whole array operation
      do j=1,n ! for each value in f:
      
         u = (x - x(j)) / radius ! whole array operation
         do i=1,n
            u(i) = min( 1.0_ReKi, max( -1.0_ReKi, u(i) ) )
         end do
         
         k_sum   = 0.0_ReKi
         do i=1,n
            k = w*(1.0_ReKi-abs(u(i))**Exp1)**Exp2;
            k_sum = k_sum + k
            fNew(j) = fNew(j) + k*f(i)
         end do
         if (k_sum > 0.0_ReKi) then
            fNew(j) = fNew(j) / k_sum
         end if
         
      end do ! j (each output value)
      
   else ! kernelType_GAUSSIAN
      w = 1.0_ReKi/sqrt(TwoPi)
      
      fNew = 0.0_ReKi ! whole array operation
      do j=1,n ! for each value in f:
      
         u = (x - x(j)) / radius ! whole array operation
      
         k_sum   = 0.0_ReKi
         do i=1,n
            k = w*exp(-0.5*u(i)**2);
            k_sum = k_sum + k
            fNew(j) = fNew(j) + k*f(i)
         end do
         if (k_sum > 0.0_ReKi) then
            fNew(j) = fNew(j) / k_sum
         end if
         
      end do ! j (each output value)
      
   end if

end subroutine kernelSmoothing
!=======================================================================
!> This subroutine finds the lower-bound index of an input x-value located in an array.
!! On return, Ind has a value such that
!!           XAry(Ind) <= XVal < XAry(Ind+1), with the exceptions that
!!             Ind = 0 when XVal < XAry(1), and
!!          Ind = AryLen when XAry(AryLen) <= XVal.
!!
!! It uses a binary interpolation scheme that takes about log(AryLen)/log(2) steps to converge.
!! If the index doesn't change much between calls, LocateStp() (nwtc_num::locatestp) may be a better option.
   SUBROUTINE LocateBin( XVal, XAry, Ind, AryLen )

      ! Argument declarations.

   INTEGER, INTENT(IN)          :: AryLen                                          !< Length of the array.
   INTEGER, INTENT(OUT)         :: Ind                                             !< Final (low) index into the array.

   REAL(ReKi), INTENT(IN)       :: XAry    (AryLen)                                !< Array of X values to be interpolated.
   REAL(ReKi), INTENT(IN)       :: XVal                                            !< X value to be interpolated.


      ! Local declarations.

   INTEGER                      :: IHi                                             ! The high index into the arrays.
   INTEGER                      :: IMid                                            ! The mid-point index between IHi and Ind.



      ! Let's check the limits first.

   IF ( XVal < XAry(1) )  THEN
      Ind = 0
   ELSE IF ( XVal >= XAry(AryLen) )  THEN
      Ind = AryLen
   ELSE
         ! Let's interpolate!

      Ind  = 1
      IHi  = AryLen

      DO WHILE ( IHi-Ind > 1 )

         IMid = ( IHi + Ind )/2

         IF ( XVal >= XAry(IMid) ) THEN
            Ind = IMid
         ELSE
            IHi = IMid
         END IF

      END DO

   END IF

   RETURN
   END SUBROUTINE LocateBin
!=======================================================================
!> This subroutine finds the lower-bound index of an input x-value located in an array.
!! On return, Ind has a value such that
!!           XAry(Ind) <= XVal < XAry(Ind+1), with the exceptions that
!!             Ind = 0 when XVal < XAry(1), and
!!          Ind = AryLen when XAry(AryLen) <= XVal.
!!
!! It uses the passed index as the starting point and does a stepwise search from there.  This is
!! especially useful when the calling routines save the value from the last time this routine was called
!! for a given case where XVal does not change much from call to call.  When there is no correlation
!! from one interpolation to another, a binary search may be a better choice (see nwtc_num::locatebin).
!!
!! Use LocateStp (nwtc_num::locatestp) instead of directly calling a specific routine in the generic interface.    
   SUBROUTINE LocateStpR4( XVal, XAry, Ind, AryLen )

      ! Argument declarations.

   INTEGER, INTENT(IN)          :: AryLen                                          !< Length of the array.
   INTEGER, INTENT(INOUT)       :: Ind                                             !< Initial and final index into the array.

   REAL(SiKi), INTENT(IN)       :: XAry    (AryLen)                                !< Array of X values to be interpolated.
   REAL(SiKi), INTENT(IN)       :: XVal                                            !< X value to be interpolated.



      ! Let's check the limits first.

   IF ( XVal < XAry(1) )  THEN
      Ind = 0
   ELSE IF ( XVal >= XAry(AryLen) )  THEN
      Ind = AryLen
   ELSE

      Ind = MAX( MIN( Ind, AryLen-1 ), 1 )

      DO

         IF ( XVal < XAry(Ind) )  THEN

            Ind = Ind - 1

         ELSE IF ( XVal >= XAry(Ind+1) )  THEN

            Ind = Ind + 1

         ELSE

            RETURN

         END IF

      END DO


   END IF

   RETURN

   END SUBROUTINE LocateStpR4
!=======================================================================
!> \copydoc nwtc_num::locatestpr4
   SUBROUTINE LocateStpR8( XVal, XAry, Ind, AryLen )

      ! Argument declarations.

   INTEGER, INTENT(IN)          :: AryLen                                          ! Length of the array.
   INTEGER, INTENT(INOUT)       :: Ind                                             ! Initial and final index into the array.

   REAL(R8Ki), INTENT(IN)       :: XAry    (AryLen)                                ! Array of X values to be interpolated.
   REAL(R8Ki), INTENT(IN)       :: XVal                                            ! X value to be interpolated.



      ! Let's check the limits first.

   IF ( XVal < XAry(1) )  THEN
      Ind = 0
   ELSE IF ( XVal >= XAry(AryLen) )  THEN
      Ind = AryLen
   ELSE

      Ind = MAX( MIN( Ind, AryLen-1 ), 1 )

      DO

         IF ( XVal < XAry(Ind) )  THEN

            Ind = Ind - 1

         ELSE IF ( XVal >= XAry(Ind+1) )  THEN

            Ind = Ind + 1

         ELSE

            RETURN

         END IF

      END DO


   END IF

   RETURN

   END SUBROUTINE LocateStpR8
!=======================================================================
!> This routine calculates the mean value of an array.
   FUNCTION Mean ( Ary, AryLen )
      
   !NOTE: We should make AryLen an optional argument and use SIZE( Ary ) if it is not present.

      ! Function declaration.

   REAL(ReKi)                   :: Mean                                         ! The mean of the values in Ary.


      ! Argument declarations:

   INTEGER, INTENT(IN)          :: AryLen                                       !< Length of the array.

   REAL(ReKi), INTENT(IN)       :: Ary  (AryLen)                                !< Input array.


      ! Local declarations.

   REAL(DbKi)                   :: Sum                                          ! A temporary sum.

   INTEGER                      :: I                                            ! The index into the array.



   Sum = 0.0_DbKi

   DO I=1,AryLen
      Sum = Sum + Ary(I)
   END DO ! I

   Mean = Sum/AryLen


   RETURN
   END FUNCTION Mean ! ( Ary, AryLen )
!=======================================================================
!> This routine is used to convert Angle to an equivalent value
!!  between \f$-\pi\f$ and \f$pi\f$.
!!
!! Use MPi2Pi (nwtc_num::mpi2pi) instead of directly calling a specific routine in the generic interface.
   SUBROUTINE MPi2Pi_R4 ( Angle )

                 
      ! Argument declarations:

   REAL(SiKi), INTENT(INOUT)    :: Angle  !< Angle (in radians) to be converted


      ! Get the angle between 0 and 2Pi.

   Angle = MODULO( Angle, TwoPi_R4 )


      ! Get the angle between -Pi and Pi.

   IF ( Angle > Pi_R4 )  THEN
      Angle = Angle - TwoPi_R4
   END IF


   RETURN
   END SUBROUTINE MPi2Pi_R4
!=======================================================================
!> \copydoc nwtc_num::mpi2pi_r4
   SUBROUTINE MPi2Pi_R8 ( Angle )

                 
      ! Argument declarations:

   REAL(R8Ki), INTENT(INOUT)    :: Angle


      ! Get the angle between 0 and 2Pi.

   Angle = MODULO( Angle, TwoPi_R8 )


      ! Get the angle between -Pi and Pi.

   IF ( Angle > Pi_R8 )  THEN
      Angle = Angle - TwoPi_R8
   END IF


   RETURN
   END SUBROUTINE MPi2Pi_R8
!=======================================================================
!> This function takes an angle in radians and converts it to 
!! an angle in degrees in the range [-180,180]
real(reKi) function Rad2M180to180Deg(Angle) result(Alpha)

   real(ReKi),             intent(in   ) :: Angle !< input angle, radians
   
   Alpha = Angle
   
   call MPi2Pi ( Alpha ) ! change Angle into range of -pi to pi
   Alpha = Alpha*R2D     ! change Angle into degrees
   
end function Rad2M180to180Deg
!=======================================================================
   
!> This routine calculates the outer product of two vectors, 
!! \f$u = \left(u_1, u_2, \ldots, u_m\right)\f$ and 
!! \f$v = \left(v_1, v_2, \ldots ,v_n\right)\f$. The outer product is defined as
!! \f{equation}{
!!   A = u \otimes v = \begin{bmatrix}
!!   u_1 v_1 & u_1 v_2 & \dots  & u_1 v_n \\
!!   u_2 v_1 & u_2 v_2 & \dots  & u_2 v_n \\
!!    \vdots & \vdots  & \ddots & \vdots \\
!!   u_m v_1 & u_m v_2 & \dots  & u_m v_n
!!   \end{bmatrix}  
!! \f}   
!!
!! Use OuterProduct (nwtc_num::outerproduct) instead of directly calling a specific routine in the generic interface.    
   FUNCTION OuterProductR4(u,v)
   

   REAL(SiKi),  INTENT(IN) :: u(:)                                !< first vector, \f$u\f$, in the outer product
   REAL(SiKi),  INTENT(IN) :: v(:)                                !< second vector, \f$v\f$, in the outer product
   REAL(SiKi)              :: OuterProductR4(SIZE(u),SIZE(v))     !< the resultant matrix, A

   INTEGER(IntKi)::i,j,n1,n2

   n1=SIZE(u)
   n2=SIZE(v)

   DO i=1,n1
       DO j=1,n2
           OuterProductR4(i,j) = u(i) * v(j)
       ENDDO
   ENDDO

   END FUNCTION OuterProductR4   
!=======================================================================
!> \copydoc nwtc_num::outerproductr4
   FUNCTION OuterProductR8(u,v)
   
   ! this routine calculates the outer product of two vectors

   REAL(R8Ki),INTENT(IN):: u(:)
   REAL(R8Ki),INTENT(IN):: v(:)
   REAL(R8Ki)::OuterProductR8(SIZE(u),SIZE(v))

   INTEGER(IntKi)::i,j,n1,n2

   n1=SIZE(u)
   n2=SIZE(v)

   DO i=1,n1
       DO j=1,n2
           OuterProductR8(i,j) = u(i) * v(j)
       ENDDO
   ENDDO

   END FUNCTION OuterProductR8   
!=======================================================================
!> This subroutine perturbs an orientation matrix by a small angle.  Two methods
!! are used:
!!    small angle DCM:  perturb small angles extracted from DCM
!!    large angle DCM:  multiply input DCM with DCM created with small angle
!!                      perturbations
!! NOTE1: this routine originally used logarithmic mapping for small angle
!!          perturbations
!! NOTE2: all warnings from SmllRotTrans are ignored.
!! NOTE3: notice no checks are made to verify correct set of inputs given
!!          one of the follwing combinations must be provided (others truly
!!          optional):
!!             Perturbations
!!             Perturbation + AngleDim
   SUBROUTINE PerturbOrientationMatrix( Orientation, Perturbation, AngleDim, Perturbations, UseSmlAngle )
      REAL(R8Ki),           INTENT(INOUT)  :: Orientation(3,3)
      REAL(R8Ki), OPTIONAL, INTENT(IN)     :: Perturbation ! angle (radians) of the perturbation
      INTEGER,    OPTIONAL, INTENT(IN)     :: AngleDim
      REAL(R8Ki), OPTIONAL, INTENT(IN)     :: Perturbations(3) ! angles (radians) of the perturbations
      LOGICAL,    OPTIONAL, INTENT(IN)     :: UseSmlAngle
   
           ! Local variables
      REAL(R8Ki)                 :: angles(3)
      REAL(R8Ki)                 :: OrientationTmp(3,3)
      LOGICAL                    :: OutputSmlAngle
      integer(intKi)             :: ErrStat2
      character(ErrMsgLen)       :: ErrMsg2
      
      if (present(UseSmlAngle)) then
         OutputSmlAngle = UseSmlAngle
      else
         OutputSmlAngle = .false.
      end if

      if (OutputSmlAngle) then
!         CALL DCM_LogMap( Orientation, angles, ErrStat2, ErrMsg2 )
         angles =  GetSmllRotAngs ( Orientation, ErrStat2, ErrMsg2 )
         IF (PRESENT(Perturbations)) THEN
            angles = angles + Perturbations
         ELSE
            angles(AngleDim) = angles(AngleDim) + Perturbation
         END IF
!         Orientation = DCM_exp( angles )
         call SmllRotTrans( 'linearization perturbation', angles(1), angles(2), angles(3), Orientation, ErrStat=ErrStat2, ErrMsg=ErrMsg2 )
      else !Only works if AngleDim is specified
         IF (PRESENT(Perturbations)) THEN
            angles = Perturbations
         ELSE
            angles = 0.0_R8Ki
            angles(AngleDim) = Perturbation
         END IF
         call SmllRotTrans( 'linearization perturbation', angles(1), angles(2), angles(3), OrientationTmp, ErrStat=ErrStat2, ErrMsg=ErrMsg2 )
         Orientation = matmul(Orientation, OrientationTmp)
      endif
      

   END SUBROUTINE PerturbOrientationMatrix
!=======================================================================
!> This routine factors the number N into its primes. If any of those
!! prime factors is greater than the NumPrimes'th prime, a value of 1
!! is added to N and the new number is factored.  This process is 
!! repeated until no prime factors are greater than the NumPrimes'th 
!! prime.
!!
!! If subract is .true., we will subtract 1 from the value of N instead
!! of adding it.
   FUNCTION PSF ( Npsf, NumPrimes, subtract )


    IMPLICIT                 NONE

    !Passed variables
    INTEGER,         INTENT(IN) :: Npsf                   !< Initial number we're trying to factor.
    INTEGER,         INTENT(IN) :: NumPrimes              !< Number of unique primes.
    INTEGER                     :: PSF                    !< The smallest number at least as large as Npsf, that is the product of small factors when we return.
                                                          !! IF subtract is present and .TRUE., PSF is the largest number not greater than Npsf that is a  product of small factors.
    LOGICAL,OPTIONAL,INTENT(IN) :: subtract               !< if PRESENT and .TRUE., we will subtract instead of add 1 to the number when looking for the value of PSF to return.
    
    !Other variables
    INTEGER                     :: sign                   ! +1 or -1 
    INTEGER                     :: IPR                    ! A counter for the NPrime array
    INTEGER, PARAMETER          :: NFact = 9              ! The number of prime numbers (the first NFact primes)
    INTEGER                     :: NP                     ! A temp variable to determine if NPr divides NTR
    INTEGER                     :: NPr                    ! A small prime number
    INTEGER                     :: NT                     ! A temp variable to determine if NPr divides NTR: INT( NTR / NPr )
    INTEGER                     :: NTR                    ! The number we're trying to factor in each iteration
    INTEGER, PARAMETER          :: NPrime(NFact) = (/ 2, 3, 5, 7, 11, 13, 17, 19, 23 /) ! The first 9 prime numbers
                              
    LOGICAL                     :: DividesN1(NFact)       ! Does this factor divide NTR-1?



    DividesN1(:) = .FALSE.                              ! We need to check all of the primes the first time through
    
    sign = 1
    IF ( PRESENT( subtract ) ) THEN
       IF (subtract) THEN
          sign = -1
       END IF
    END IF
    
    PSF = Npsf

    DO
           ! First:  Factor NTR into its primes.

       NTR = PSF

       DO IPR=1,MIN( NumPrimes, NFact ) 

           IF ( DividesN1(IPR) ) THEN

                   ! If P divides N-1, then P cannot divide N.

               DividesN1(IPR) = .FALSE.               ! This prime number does not divide psf; We'll check it next time.

           ELSE

               NPr = NPrime(IPR)                      ! The small prime number we will try to find the the factorization of NTR

               DO
                   NT = NTR/NPr                       ! Doing some modular arithmetic to see if
                   NP = NT*NPr                        ! MOD( NTR, NPr ) == 0, i.e. if NPr divides NTR

                   IF ( NP /= NTR )  EXIT             ! There aren't any more of this prime number in the factorization

                   NTR = NT                           ! This is the new number we need to get factors for
                   DividesN1(IPR) = .TRUE.            ! This prime number divides psf, so we won't check it next time (on Npsf+1).

               ENDDO

               IF ( NTR .EQ. 1 )  RETURN              ! We've found all the prime factors, so we're finished

           ENDIF !  DividesN1

       ENDDO ! IPR

           ! Second:  There is at least one prime larger than NPrime(NumPrimes).  Add or subtract
           !          a point to NTR and factor again.

       PSF = PSF + sign*1

    ENDDO


    RETURN
    END FUNCTION PSF   
!=======================================================================  
!> This function computes the conjugate of a quaternion, q
   FUNCTION Quaternion_Conjugate(q)

      !
      ! "'Interpolation' of DCMs", M.A. Sprague, 11 March 2014, Eq. 6
   
   TYPE(Quaternion), INTENT(IN)    :: q                     !< quaternion
   
   TYPE(Quaternion)                :: Quaternion_Conjugate  !< conjugate of the quaternion
   
      
   Quaternion_Conjugate%q0 =  q%q0 
   Quaternion_Conjugate%v  = -q%v
      
   END FUNCTION Quaternion_Conjugate   
!=======================================================================  
!> This function computes the 2-norm of a quaternion, q
   FUNCTION Quaternion_Norm(q)
      
      ! "'Interpolation' of DCMs", M.A. Sprague, 11 March 2014, Eq. 5
   
   TYPE(Quaternion), INTENT(IN)    :: q                  !< quaternion
   
   REAL(ReKi)                      :: Quaternion_Norm    !< 2-norm of q
   
      
   Quaternion_Norm = sqrt( q%q0**2 + DOT_PRODUCT(q%v, q%v) )
   
   
   END FUNCTION Quaternion_Norm   
!=======================================================================  
!> This function computes the quaternion, q, raised to an arbitrary
!! real exponent, alpha.
   FUNCTION Quaternion_Power(q,alpha)

      ! "'Interpolation' of DCMs", M.A. Sprague, 11 March 2014, Eq. 7-8
   
   TYPE(Quaternion), INTENT(IN)    :: q                  !< quaternion
   REAL(ReKi)      , INTENT(IN)    :: alpha              !< exponent
   
   TYPE(Quaternion)                :: Quaternion_Power   !< q^alpha
   
   
      ! local variables
   REAL(ReKi)                      :: greek   ! the product of alpha and theta
   REAL(ReKi)                      :: n(3)
   REAL(ReKi)                      :: q_norm
   REAL(ReKi)                      :: q_norm_power
   REAL(ReKi)                      :: theta
      
   
   q_norm       = Quaternion_Norm( q )     
   theta        = acos( max(-1.0_ReKi, min(1.0_ReKi, q%q0 / q_norm)) )
   n            = q%v / TwoNorm(q%v)
   
   greek        = alpha * theta
   q_norm_power = q_norm ** alpha
   
   Quaternion_Power%q0 =  q_norm_power * cos( greek )
   Quaternion_Power%v  =  q_norm_power * sin( greek ) * n
      
   END FUNCTION Quaternion_Power   
!=======================================================================  
!> This function computes the product of two quaternions, p and q
   FUNCTION Quaternion_Product(p, q)

      ! "'Interpolation' of DCMs", M.A. Sprague, 11 March 2014, Eq. 4
   
   TYPE(Quaternion), INTENT(IN)    :: p                     !< quaternion   
   TYPE(Quaternion), INTENT(IN)    :: q                     !< quaternion 
   
   TYPE(Quaternion)                :: Quaternion_Product    !< quaternion product, p*q
   
      
   Quaternion_Product%q0 = p%q0 * q%q0 - DOT_PRODUCT(p%v, q%v)
   Quaternion_Product%v  = p%q0*q%v + q%q0*p%v + CROSS_PRODUCT( p%v, q%v ) 
   
   
   END FUNCTION Quaternion_Product   
!=======================================================================  
!> This function converts a quaternion to an equivalent direction cosine matrix.
   FUNCTION Quaternion_to_DCM(q)

      ! "'Interpolation' of DCMs", M.A. Sprague, 11 March 2014, Eq. 9-17
   
   TYPE(Quaternion), INTENT(IN)    :: q                        !< quaternion    
   
   REAL(ReKi)                      :: Quaternion_to_DCM (3,3)  ! equivalent direction cosine matrix
   
      ! local variables (products of quaternion terms)
   REAL(ReKi)                      :: q0q0, q0q1, q0q2, q0q3
   REAL(ReKi)                      :: q1q1, q1q2, q1q3 
   REAL(ReKi)                      :: q2q2, q2q3
   REAL(ReKi)                      :: q3q3
   
   q0q0 = q%q0**2
   q0q1 = q%q0      * q%v(1)
   q0q2 = q%q0      * q%v(2)
   q0q3 = q%q0      * q%v(3)
   
   q1q1 = q%v(1)**2
   q1q2 = q%v(1)    * q%v(2)
   q1q3 = q%v(1)    * q%v(3)

   q2q2 = q%v(2)**2
   q2q3 = q%v(2)    * q%v(3)
   
   q3q3 = q%v(2)**2
   
   
   Quaternion_to_DCM(1,1) =          q0q0 +          q1q1 - q2q2 - q3q3  ! Eq.  9
   Quaternion_to_DCM(1,2) = 2.0_ReKi*q1q2 + 2.0_ReKi*q0q3                ! Eq. 10
   Quaternion_to_DCM(1,3) = 2.0_ReKi*q1q3 + 2.0_ReKi*q0q2                ! Eq. 11

   Quaternion_to_DCM(2,1) = 2.0_ReKi*q1q2 - 2.0_ReKi*q0q3                ! Eq. 12
   Quaternion_to_DCM(2,2) =          q0q0 -          q1q1 + q2q2 - q3q3  ! Eq. 13
   Quaternion_to_DCM(2,3) = 2.0_ReKi*q2q3 +          q0q1                ! Eq. 14
   
   
   Quaternion_to_DCM(3,1) = 2.0_ReKi*q1q3 +          q0q2                ! Eq. 15
   Quaternion_to_DCM(3,2) = 2.0_ReKi*q2q3 -          q0q1                ! Eq. 16 
   Quaternion_to_DCM(3,3) =          q0q0 -          q1q1 - q2q2 + q3q3  ! Eq. 17
   
   
   END FUNCTION Quaternion_to_DCM   
!=======================================================================  
!> This function converts a direction cosine matrix to an equivalent quaternion.
   FUNCTION DCM_to_Quaternion(DCM)

      ! "'Interpolation' of DCMs", M.A. Sprague, 11 March 2014, Eq. 18-21
   
   REAL(ReKi)      , INTENT(IN)    :: DCM (3,3)          !< direction cosine matrix
   TYPE(Quaternion)                :: DCM_to_Quaternion  !< equivalent quaternion      
   
         
   DCM_to_Quaternion%q0   =      0.5_ReKi * sqrt( 1.0_ReKi + DCM(1,1) + DCM(2,2) + DCM(3,3) )                         ! Eq. 18
   DCM_to_Quaternion%v(1) = sign(0.5_ReKi * sqrt( 1.0_ReKi + DCM(1,1) - DCM(2,2) - DCM(3,3) ) , DCM(2,3) - DCM(3,2) ) ! Eq. 19
   DCM_to_Quaternion%v(2) = sign(0.5_ReKi * sqrt( 1.0_ReKi - DCM(1,1) + DCM(2,2) - DCM(3,3) ) , DCM(3,1) - DCM(1,3) ) ! Eq. 20
   DCM_to_Quaternion%v(3) = sign(0.5_ReKi * sqrt( 1.0_ReKi - DCM(1,1) - DCM(2,2) + DCM(3,3) ) , DCM(1,2) - DCM(2,1) ) ! Eq. 21

   
   
   END FUNCTION DCM_to_Quaternion
!=======================================================================         
!> This function computes the interpolated quaternion at time
!! t1 + s*(t2-t1) and s is in [0,1]
   FUNCTION Quaternion_Interp(q1,q2,s)

      ! "'Interpolation' of DCMs", M.A. Sprague, 11 March 2014, Eq. 23
   
   TYPE(Quaternion), INTENT(IN)    :: q1              !< quaternion at t1    
   TYPE(Quaternion), INTENT(IN)    :: q2              !< quaternion at t2
   REAL(ReKi),       INTENT(IN)    :: s               !< fraction of distance between t1 and t2: s must be in [0,1]
   
   TYPE(Quaternion)                :: Quaternion_Interp !< interpolated quaternion at s
   
      
   Quaternion_Interp = Quaternion_Conjugate(q1)
   Quaternion_Interp = Quaternion_Product(Quaternion_Interp, q2)
   Quaternion_Interp = Quaternion_Power(  Quaternion_Interp, s )
   Quaternion_Interp = Quaternion_Product(q1, Quaternion_Interp)
   
   
! bjj: this function has not been tested. I have not tested any of the quaternion routines, either. 

   END FUNCTION Quaternion_Interp
!=======================================================================
!> This routine calculates the parameters needed to compute a regularly-spaced natural cubic spline.
!! Natural cubic splines are used in that the curvature at the end points is zero.
!! It assumes the XAry values are equally spaced for speed. If you have multiple curves that share the 
!! same value, use RegCubicSplineInitM (nwtc_num::regcubicsplineinitm) instead of calling this routine multiple times.
   SUBROUTINE RegCubicSplineInit ( AryLen, XAry, YAry, DelX, Coef, ErrStat, ErrMsg )

      ! Argument declarations:

   INTEGER, INTENT(IN)          :: AryLen                                     !< Length of the array.

   REAL(ReKi), INTENT(OUT)      :: Coef  (AryLen-1,0:3)                       !< The coefficients for the cubic polynomials.
   REAL(ReKi), INTENT(OUT)      :: DelX                                       !< The distance between the equally spaced points.
   REAL(ReKi), INTENT(IN)       :: XAry  (AryLen)                             !< Input array of x values.
   REAL(ReKi), INTENT(IN)       :: YAry  (AryLen)                             !< Input array of y values.

   INTEGER(IntKi), INTENT(OUT)  :: ErrStat                                    !< Error status.

   CHARACTER(*),   INTENT(OUT)  :: ErrMsg                                     !< Error message.


      ! Local declarations.

   REAL(ReKi)                   :: DelX2                                      ! The square of the distance between points.
   REAL(ReKi)                   :: DelX4                                      ! Four times the distance between points.
   REAL(ReKi)                   :: DelX6                                      ! Six times the distance between points.
   REAL(ReKi), ALLOCATABLE      :: Slope (:)                                  ! The AryLen-1 length array of slopes between points.
   REAL(ReKi), ALLOCATABLE      :: U     (:)                                  ! An AryLen-1 length array used in the Gaussian elimination.
   REAL(ReKi), ALLOCATABLE      :: V     (:)                                  ! An AryLen-1 length array used in the Gaussian elimination.
   REAL(ReKi)                   :: ZHi                                        ! A parameter used to calculate the polynomial coefficients.
   REAL(ReKi)                   :: ZLo                                        ! A parameter used to calculate the polynomial coefficients.

   INTEGER(IntKi)               :: ErrStatLcL                                 ! Local error status.
   INTEGER                      :: I                                          ! The index into the arrays.
   CHARACTER(*), PARAMETER      :: RoutineName = 'RegCubicSplineInit'


      ! Allocate the various intermediate arrays.

   ALLOCATE ( Slope( AryLen - 1 ), STAT=ErrStatLcL )
   IF ( ErrStatLcL /= 0 )  THEN
      CALL ExitThisRoutine ( ErrID_Fatal, RoutineName//':Error allocating memory for the Slope array.' )
      RETURN
   ENDIF

   ALLOCATE ( U( AryLen - 1 ), STAT=ErrStatLcL )
   IF ( ErrStatLcL /= 0 )  THEN
      CALL ExitThisRoutine ( ErrID_Fatal, RoutineName//':Error allocating memory for the U array.' )
      RETURN
   ENDIF

   ALLOCATE ( V( AryLen - 1 ), STAT=ErrStatLcL )
   IF ( ErrStatLcL /= 0 )  THEN
      CALL ExitThisRoutine ( ErrID_Fatal, RoutineName//':Error allocating memory for the V array.' )
      RETURN
   ENDIF


      ! Compute the distance between XAry values and the slopes between points.

   DelX  = ( XAry(AryLen) - XAry(1) )/REAL( AryLen-1, ReKi )                   ! Is this more accurate than XAry(2) - XAry(1)?
   DelX2 = DelX*DelX
   DelX4 = 4_ReKI*DelX
   DelX6 = 6_ReKI*DelX

   DO I=1,AryLen-1
      Slope(I) = ( YAry(I+1) - YAry(I) )/DelX
   END DO ! I


      ! Use Gaussian elimination to solve the tri-diagonal matrix.

   U(1) = DelX4
   V(1) = 6.0_ReKi*( Slope(2) - Slope(1) )

   DO I=2,AryLen-1
      U(I) = DelX4 - DelX2/U(I-1)
      V(I) = 6.0_ReKi*( Slope(I) - Slope(I-1) ) - DelX*V(I-1)/U(I-1)
   END DO ! I


      ! Determine the coefficients of the polynomials.

   Coef(:,0) = YAry(1:AryLen-1)

   ZHi = 0.0_ReKi

   DO I=AryLen-1,1,-1
      ZLo       = ( V(I) - DelX*ZHi )/U(I)
      Coef(I,1) = Slope(I) - DelX*( ZHi/6.0_ReKi + ZLo/3.0_ReKi )
      Coef(I,2) = 0.5_ReKi*ZLo
      Coef(I,3) = ( ZHi - ZLo )/DelX6
      ZHi       = ZLo
   END DO ! I


   CALL ExitThisRoutine ( ErrID_None, 'No Problemo' )

   RETURN

   !=======================================================================
   CONTAINS
   !=======================================================================
      SUBROUTINE ExitThisRoutine ( ErrID, Msg )

         ! This subroutine cleans up the parent routine before exiting.


            ! Argument declarations.

         INTEGER(IntKi), INTENT(IN)       :: ErrID                            ! The error identifier (ErrLev)

         CHARACTER(*),   INTENT(IN)       :: Msg                              ! The error message (ErrMsg)


            ! Local declarations.

         LOGICAL                          :: IsOpen                           ! A flag that indicates if the input unit is still open.


            ! Set error status/message

         ErrStat = ErrID
         ErrMsg  = Msg


            ! Deallocate the Words array if it had been allocated.

         IF ( ALLOCATED( Slope ) )  DEALLOCATE( Slope )
         IF ( ALLOCATED( U     ) )  DEALLOCATE( U     )
         IF ( ALLOCATED( V     ) )  DEALLOCATE( V     )


         RETURN

      END SUBROUTINE ExitThisRoutine ! ( ErrID, Msg )

   END SUBROUTINE RegCubicSplineInit ! ( AryLen, XAry, YAry, DelX, Coef, ErrStat, ErrMsg )
!=======================================================================
!> This routine calculates the parameters needed to compute a regularly-spaced natural cubic spline.
!! Natural cubic splines are used in that the curvature at the end points is zero.
!! It assumes the XAry values are equally spaced for speed.
!! This version of the routine works with multiple curves that share the same X values.
   SUBROUTINE RegCubicSplineInitM ( XAry, YAry, DelX, Coef, ErrStat, ErrMsg )

      ! Argument declarations:

   REAL(ReKi), INTENT(OUT)      :: Coef  (:,:,0:)                             ! The coefficients for the cubic polynomials.
   REAL(ReKi), INTENT(OUT)      :: DelX                                       ! The distance between X values in XAry.
   REAL(ReKi), INTENT(IN)       :: XAry  (:)                                  ! Input array of regularly spaced x values.
   REAL(ReKi), INTENT(IN)       :: YAry  (:,:)                                ! Input array of y values.

   INTEGER(IntKi), INTENT(OUT)  :: ErrStat                                    ! Error status.

   CHARACTER(*),   INTENT(OUT)  :: ErrMsg                                     ! Error message.
   

      ! Local declarations.

   REAL(ReKi)                   :: DelX2                                      ! The square of the distance between points.
   REAL(ReKi)                   :: DelX4                                      ! Four times the distance between points.
   REAL(ReKi)                   :: DelX6                                      ! Six times the distance between points.
   REAL(ReKi), ALLOCATABLE      :: Slope (:,:)                                ! The NumPts-1 length array of slopes between points.
   REAL(ReKi), ALLOCATABLE      :: U     (:)                                  ! An NumPts-1 length array used in the Gaussian elimination.
   REAL(ReKi), ALLOCATABLE      :: V     (:,:)                                ! An NumPts-1 length array used in the Gaussian elimination.
   REAL(ReKi), ALLOCATABLE      :: ZHi   (:)                                  ! A parameter used to calculate the polynomial coefficients.
   REAL(ReKi), ALLOCATABLE      :: ZLo   (:)                                  ! A parameter used to calculate the polynomial coefficients.

   INTEGER(IntKi)               :: ErrStatLcL                                 ! Local error status.
   INTEGER                      :: I                                          ! The index into the arrays.
!   INTEGER                      :: IC                                         ! The curve index into the arrays.
   INTEGER                      :: NumCrvs                                    ! Number of curves to be interpolated.
   INTEGER                      :: NumPts                                     ! Number of points in each curve.

   CHARACTER(*), PARAMETER      :: RoutineName = 'RegCubicSplineInitM'

      ! How big are the arrays?

   NumPts  = SIZE( XAry )
   NumCrvs = SIZE( YAry, 2 )


      ! Allocate the various intermediate arrays.

   ALLOCATE ( ZLo( NumCrvs ), STAT=ErrStatLcL )
   IF ( ErrStatLcL /= 0 )  THEN
      CALL ExitThisRoutine ( ErrID_Fatal, RoutineName//':Error allocating memory for the ZLo array.' )
      RETURN
   ENDIF

   ALLOCATE ( ZHi( NumCrvs ), STAT=ErrStatLcL )
   IF ( ErrStatLcL /= 0 )  THEN
      CALL ExitThisRoutine ( ErrID_Fatal, RoutineName//':Error allocating memory for the ZHi array.' )
      RETURN
   ENDIF

   ALLOCATE ( Slope( NumPts-1, NumCrvs ), STAT=ErrStatLcL )
   IF ( ErrStatLcL /= 0 )  THEN
      CALL ExitThisRoutine ( ErrID_Fatal, RoutineName//':Error allocating memory for the Slope array.' )
      RETURN
   ENDIF

   ALLOCATE ( U( NumPts - 1 ), STAT=ErrStatLcL )
   IF ( ErrStatLcL /= 0 )  THEN
      CALL ExitThisRoutine ( ErrID_Fatal, RoutineName//':Error allocating memory for the U array.' )
      RETURN
   ENDIF

   ALLOCATE ( V( NumPts-1, NumCrvs ), STAT=ErrStatLcL )
   IF ( ErrStatLcL /= 0 )  THEN
      CALL ExitThisRoutine ( ErrID_Fatal, RoutineName//':Error allocating memory for the V array.' )
      RETURN
   ENDIF


      ! Compute the distance between XAry values and the slopes between points.

   DelX  = ( XAry(NumPts) - XAry(1) )/REAL( NumPts-1, ReKi )                  ! Is this more accurate than XAry(2) - XAry(1)?
   DelX2 = DelX*DelX
   DelX4 = 4_ReKI*DelX
   DelX6 = 6_ReKI*DelX

   DO I=1,NumPts-1
      Slope(I,:) = ( YAry(I+1,:) - YAry(I,:) )/DelX
   END DO ! I


      ! Use Gaussian elimination to solve the tri-diagonal matrix.

   U(1) = DelX4

   DO I=2,NumPts-1
      U(I) = DelX4 - DelX2/U(I-1)
   END DO ! I

   V(1,:) = 6.0_ReKi*( Slope(2,:) - Slope(1,:) )

   DO I=2,NumPts-1
      V(I,:) = 6.0_ReKi*( Slope(I,:) - Slope(I-1,:) ) - DelX*V(I-1,:)/U(I-1)
   END DO ! I


      ! Determine the coefficients of the polynomials.

   Coef(:,:,0) = YAry(1:NumPts-1,:)

   ZHi(:) = 0.0_ReKi

   DO I=NumPts-1,1,-1
      ZLo(:)      = ( V(I,:) - DelX*ZHi(:) )/U(I)
      Coef(I,:,1) = Slope(I,:) - DelX*( ZHi(:)/6.0_ReKi + ZLo(:)/3.0_ReKi )
      Coef(I,:,2) = 0.5_ReKi*ZLo(:)
      Coef(I,:,3) = ( ZHi(:) - ZLo(:) )/DelX6
      ZHi(:)      = ZLo(:)
   END DO ! I


   CALL ExitThisRoutine ( ErrID_None, 'No Problemo' )

   RETURN

   !=======================================================================
   CONTAINS
   !=======================================================================
      SUBROUTINE ExitThisRoutine ( ErrID, Msg )

         ! This subroutine cleans up the parent routine before exiting.


            ! Argument declarations.

         INTEGER(IntKi), INTENT(IN)       :: ErrID                            ! The error identifier (ErrLev)

         CHARACTER(*),   INTENT(IN)       :: Msg                              ! The error message (ErrMsg)


            ! Local declarations.

         LOGICAL                          :: IsOpen                           ! A flag that indicates if the input unit is still open.


            ! Set error status/message

         ErrStat = ErrID
         ErrMsg  = Msg


            ! Deallocate the Words array if it had been allocated.

         IF ( ALLOCATED( ZHi   ) )  DEALLOCATE( ZHi   )
         IF ( ALLOCATED( ZLo   ) )  DEALLOCATE( ZLo   )
         IF ( ALLOCATED( Slope ) )  DEALLOCATE( Slope )
         IF ( ALLOCATED( U     ) )  DEALLOCATE( U     )
         IF ( ALLOCATED( V     ) )  DEALLOCATE( V     )


         RETURN

      END SUBROUTINE ExitThisRoutine ! ( ErrID, Msg )

   END SUBROUTINE RegCubicSplineInitM ! ( XAry, YAry, DelX, Coef, ErrStat, ErrMsg )
!=======================================================================
!> This routine interpolates a pair of arrays using cubic splines to find the function value at X.
!! One must call RegCubicSplineInit() (nwtc_num::regcubicsplineinit) first to compute the coefficients of the cubics.
!! This routine requires that the XAry be regularly spaced, which improves performance.
   FUNCTION RegCubicSplineInterp ( X, AryLen, XAry, YAry, DelX, Coef, ErrStat, ErrMsg )

      ! Function declaration.

   REAL(ReKi)                   :: RegCubicSplineInterp                       !< This function.


      ! Argument declarations:

   INTEGER, INTENT(IN)          :: AryLen                                     !< Length of the array.

   REAL(ReKi), INTENT(IN)       :: Coef  (AryLen-1,0:3)                       !< The coefficients for the cubic polynomials.
   REAL(ReKi), INTENT(IN)       :: DelX                                       !< The distance between X values in XAry.
   REAL(ReKi), INTENT(IN)       :: X                                          !< The value we are trying to interpolate for.
   REAL(ReKi), INTENT(IN)       :: XAry (AryLen)                              !< Input array of regularly spaced x values.
   REAL(ReKi), INTENT(IN)       :: YAry (AryLen)                              !< Input array of y values.

   INTEGER(IntKi), INTENT(OUT)  :: ErrStat                                    !< Error status.

   CHARACTER(*),   INTENT(OUT)  :: ErrMsg                                     !< Error message.


      ! Local declarations.

   REAL(ReKi)                   :: XOff                                       ! The distance from X to XAry(ILo).

   INTEGER                      :: ILo                                        ! The index into the array for which X is just above or equal to XAry(ILo).

   ErrStat = ErrID_None
   ErrMsg  = ""

      ! See if X is within the range of XAry.  Return the end point if it is not.

   IF ( X <= XAry(1) )  THEN
      RegCubicSplineInterp = YAry(1)
      RETURN
   ELSEIF ( X >= XAry(AryLen) )  THEN
      RegCubicSplineInterp = YAry(AryLen)
      RETURN
   ENDIF ! ( X <= XAry(1) )


      ! We are somewhere inside XAry.  Find the segment that bounds X.

   ILo = INT( ( X - XAry(1) )/DelX ) + 1

   XOff = X - XAry(ILo)

   RegCubicSplineInterp = Coef(ILo,0) + XOff*( Coef(ILo,1) + XOff*( Coef(ILo,2) + XOff*Coef(ILo,3) ) )


   RETURN

   END FUNCTION RegCubicSplineInterp ! ( X, AryLen, XAry, YAry, DelX, Coef, ErrStat, ErrMsg )
!=======================================================================
!> This routine interpolates a pair of arrays using cubic splines to find the function value at X.
!! One must call RegCubicSplineInitM() (nwtc_num::regcubicsplineinitm) first to compute the coefficients of the cubics.
!! This routine requires that the XAry be regularly spaced, which improves performance.
!! This version of the routine works with multiple curves that share the same X values.
   FUNCTION RegCubicSplineInterpM ( X, XAry, YAry, DelX, Coef, ErrStat, ErrMsg ) RESULT( Res )

      ! Function declaration.

   REAL(ReKi), ALLOCATABLE      :: Res(:)                                     !< The result of this function.


      ! Argument declarations:

   REAL(ReKi), INTENT(IN)       :: Coef  (:,:,0:)                             !< The coefficients for the cubic polynomials.
   REAL(ReKi), INTENT(IN)       :: DelX                                       !< The distance between X values in XAry.
   REAL(ReKi), INTENT(IN)       :: X                                          !< The value we are trying to interpolate for.
   REAL(ReKi), INTENT(IN)       :: XAry (:)                                   !< Input array of regularly spaced x values.
   REAL(ReKi), INTENT(IN)       :: YAry (:,:)                                 !< Input array of y values.

   INTEGER(IntKi), INTENT(OUT)  :: ErrStat                                    !< Error status.

   CHARACTER(*),   INTENT(OUT)  :: ErrMsg                                     !< Error message.


      ! Local declarations.

   REAL(ReKi)                   :: XOff                                       ! The distance from X to XAry(ILo).

   INTEGER                      :: ErrStatLcL                                 ! Local error status.
   INTEGER                      :: ILo                                        ! The index into the array for which X is just above or equal to XAry(ILo).
   INTEGER                      :: NumCrvs                                    ! Number of curves.
   INTEGER                      :: NumPts                                     ! Number of points in each curve.


   ErrStat = ErrID_None
   ErrMsg  = ""

      ! How big are the arrays?  Use the size to allocate the result.

   NumPts  = SIZE( XAry )
   NumCrvs = SIZE( YAry, 2 )

   ALLOCATE ( Res( NumCrvs ) , STAT=ErrStatLcl )
   IF ( ErrStatLcl /= 0 )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = "RegCubicSplineInterpM:Error allocating memory for the function result."
      RETURN
   ENDIF


      ! See if X is within the range of XAry.  Return the end point if it is not.

   IF ( X <= XAry(1) )  THEN
      Res(:) = YAry(1,:)
      RETURN
   ELSEIF ( X >= XAry(NumPts) )  THEN
      Res(:) = YAry(NumPts,:)
      RETURN
   ENDIF ! ( X <= XAry(1) )


      ! We are somewhere inside XAry.  Find the segment that bounds X.

   ILo = INT( ( X - XAry(1) )/DelX ) + 1

   XOff = X - XAry(ILo)

   Res(:) = Coef(ILo,:,0) + XOff*( Coef(ILo,:,1) + XOff*( Coef(ILo,:,2) + XOff*Coef(ILo,:,3) ) )


   RETURN
   END FUNCTION RegCubicSplineInterpM ! ( X, XAry, YAry, DelX, Coef, ErrStat, ErrMsg )
!=======================================================================
!> This routine is used to integrate funciton f over the interval [a, b]. This routine
!! is useful for sufficiently smooth (e.g., analytic) integrands, integrated over
!! intervals which contain no singularities, and where the endpoints are also nonsingular.
!!
!! f is an external function. For example f(x) = 1 + x.
!!
!!   FUNCTION f(x)
!!      USE PRECISION
!!      IMPLICIT NONE
!!
!!      REAL(ReKi) f
!!      REAL(ReKi) x
!!
!!      f = 1 + x
!!
!!      RETURN
!!   END FUNCTION f
   SUBROUTINE RombergInt(f, a, b, R, err, eps, ErrStat)


      IMPLICIT NONE

         ! Argument declarations:

      REAL(ReKi), EXTERNAL              :: f               !< Integrand function name
      REAL(ReKi), INTENT(IN)            :: a               !< Lower integration limit
      REAL(ReKi), INTENT(IN)            :: b               !< Upper integration limit
      REAL(ReKi), INTENT(IN)            :: eps             !< Absolute error bound
      REAL(ReKi), INTENT(OUT)           :: R               !< The result of integration
      REAL(ReKi), INTENT(OUT)           :: err             !< Actual absolute error
      INTEGER, INTENT(OUT), OPTIONAL    :: ErrStat         !< Error status; if present, program does not abort on error

         ! Local declarations:

      INTEGER                           :: m, i, j, k, IOS
      INTEGER, PARAMETER                :: mmax = 50       ! Maximum iteration number for m
      INTEGER, PARAMETER                :: imax = 50       ! Maximum iteration number for i

      REAL(ReKi), ALLOCATABLE           :: T(:,:)
      REAL(ReKi)                        :: h               ! Step length
      REAL(ReKi)                        :: sumf

         ! Initialize T
      ALLOCATE( T( mmax, imax ), Stat=ios )
      IF (IOS /= 0) THEN
         CALL ProgAbort ( 'RombergInt: Error allocating T.', PRESENT(ErrStat) )
         IF ( PRESENT(ErrStat) ) THEN
            ErrStat = ErrID_Fatal
            RETURN
         END IF
      END IF
      
      T = 0

      T(1, 1) = 0.5*(b - a)*( f(a) + f(b) )

      k = 2
      DO m = 1, mmax-2
         h = (b-a)*(0.5)**m

         sumf = 0
         DO i = 1, 2**(m-1)
            sumf = sumf + f(a + (2*i-1)*h)
            k = k + 1
         END DO


         T( m+1, 1) = 0.5*T( m, 1 )+ h * sumf

         DO j = 1, m
            T(m-j+1, j+1) = ( 4.0**j * T(m-j+2, j) - T(m-j+1, j) )/(4.0**j - 1.0)

               ! absolute error
            err = ABS( T(m-j+1, j+1) - T( m-j+2, j ) )

               ! set k >=9 to prevent early terminations
            IF( (err .LT. eps) .and. (k >= 9) ) THEN

                  ! return the intergration result if the conditions are met
               R = T(m-j+1, j+1)

               IF( ALLOCATED(T) ) DEALLOCATE(T)

               RETURN
            END IF

         END DO

      END DO

      err = ABS( T(m-j+1, j+1) - T( m-j+2, j ) )
      R = T(m-j+1, j+1)

      IF( ALLOCATED(T) ) DEALLOCATE(T)

         ! Return error message if the maximum iteration number is reached.
      CALL ProgAbort ( ' In subroutine RombergInt, the iteration reaches the maximum number. The integration did NOT converge! ', &
                       PRESENT(ErrStat) )
      IF ( PRESENT(ErrStat) ) THEN
         ErrStat = ErrID_Fatal
         RETURN
      END IF

      RETURN
   END SUBROUTINE RombergInt
!=======================================================================
!> This routine displays a message that gives that status of the simulation and the predicted end time of day.
!! It is intended to be used with SimStatus (nwtc_num::simstatus) and SimStatus_FirstTime (nwtc_num::simstatus_firsttime).
   SUBROUTINE RunTimes( StrtTime, UsrTime1, SimStrtTime, UsrTime2, ZTime, UnSum, UsrTime_out, DescStrIn )

      IMPLICIT                        NONE

         ! Passed variables

      INTEGER   ,     INTENT(IN)          :: StrtTime (8)                              !< Start time of simulation (including initialization)
      INTEGER   ,     INTENT(IN)          :: SimStrtTime (8)                           !< Start time of simulation (after initialization)
      REAL(ReKi),     INTENT(IN)          :: UsrTime1                                  !< User CPU time for simulation initialization.
      REAL(ReKi),     INTENT(IN)          :: UsrTime2                                  !< User CPU time for simulation (without intialization)
      REAL(DbKi),     INTENT(IN)          :: ZTime                                     !< The final simulation time (not necessarially TMax)
      INTEGER(IntKi), INTENT(IN), OPTIONAL:: UnSum                                     !< optional unit number of file. If present and > 0,
      REAL(ReKi),     INTENT(OUT),OPTIONAL:: UsrTime_out                               !< User CPU time for entire run - optional value returned to calling routine

      CHARACTER(*), INTENT(IN), OPTIONAL :: DescStrIn                                 !< optional additional string to print for SimStatus
                  
         ! Local variables

      REAL(ReKi)                      :: ClckTime                                        ! Elapsed clock time for the entire run.
      REAL(ReKi)                      :: ClckTimeSim                                     ! Elapsed clock time for the simulation phase of the run.
      REAL(ReKi)                      :: Factor                                          ! Ratio of seconds to a specified time period.
      REAL(ReKi)                      :: TRatio                                          ! Ratio of simulation time to elapsed clock time.
      REAL(ReKi), PARAMETER           :: SecPerDay = 24*60*60.0_ReKi                     ! Number of seconds per day
                                      
      REAL(ReKi)                      :: UsrTime                                         ! User CPU time for entire run.
      REAL(ReKi)                      :: UsrTimeSim                                      ! User CPU time for simulation (not including initialization).
      INTEGER                         :: EndTimes (8)                                    ! An array holding the ending clock time of the simulation.
                                      
      CHARACTER( 8)                   :: TimePer
      CHARACTER(MaxWrScrLen)          :: BlankLine
      CHARACTER(10)                   :: DescStr                                        !< optional additional string to print for SimStatus


      if (present(DescStrIn)) then
         DescStr = DescStrIn
      else
         DescStr = ""
      end if

         ! Get the end times to compare with start times.

      CALL DATE_AND_TIME ( VALUES=EndTimes )
      CALL CPU_TIME ( UsrTime )
      UsrTime = MAX( 0.0_ReKi, UsrTime )  ! CPU_TIME: If a meaningful time cannot be returned, a processor-dependent negative value is returned
   

      ! Calculate the elapsed wall-clock time in seconds.

      ClckTime     = GetClockTime(StrtTime,      EndTimes)
     !ClckTimeInit = GetClockTime(StrtTime,   SimStrtTime)
      ClckTimeSim  = GetClockTime(SimStrtTime,   EndTimes)

         ! Calculate CPU times.

      UsrTime    = MAX( 0.0_ReKi, UsrTime - UsrTime1 )
      UsrTimeSim = MAX( 0.0_ReKi, UsrTime - UsrTime2 )


      IF ( .NOT. EqualRealNos( UsrTimeSim, 0.0_ReKi ) .AND. ZTime > 0.0_DbKi )  THEN

         TRatio = REAL(ZTime,ReKi) / UsrTimeSim

         IF     ( UsrTime > SecPerDay )  THEN
            Factor = 1.0_ReKi/SecPerDay
            TimePer = ' days'
         ELSEIF ( UsrTime >  3600.0_ReKi )  THEN
            Factor = 1.0_ReKi/3600.0_ReKi
            TimePer = ' hours'
         ELSEIF ( UsrTime >    60.0_ReKi )  THEN
            Factor = 1.0_ReKi/60.0_ReKi
            TimePer = ' minutes'
         ELSE
            Factor = 1.0_ReKi
            TimePer = ' seconds'
         ENDIF

         BlankLine = ""
         CALL WrOver( BlankLine )  ! BlankLine contains MaxWrScrLen spaces
         CALL WrScr ( DescStr   )
         CALL WrScr ( ' Total Real Time:       '//TRIM( Num2LStr( Factor*ClckTime      ) )//TRIM( TimePer ) )
         CALL WrScr ( ' Total CPU Time:        '//TRIM( Num2LStr( Factor*UsrTime       ) )//TRIM( TimePer ) )
   !     CALL WrScr ( ' ')
   !     CALL WrScr ( ' Simulation Real Time:  '//TRIM( Num2LStr( Factor*ClckTimeSim   ) )//TRIM( TimePer ) )
         CALL WrScr ( ' Simulation CPU Time:   '//TRIM( Num2LStr( Factor*UsrTimeSim    ) )//TRIM( TimePer ) )      
         CALL WrScr ( ' Simulated Time:        '//TRIM( Num2LStr( Factor*REAL( ZTime ) ) )//TRIM( TimePer ) )
         CALL WrScr ( ' Time Ratio (Sim/CPU):  '//TRIM( Num2LStr( TRatio ) ) )

         IF (PRESENT(UnSum)) THEN
            IF (UnSum>0) THEN
               WRITE( UnSum, '(//)' )
               WRITE( UnSum, '(A)') ' Total Real Time:       '//TRIM( Num2LStr( Factor*ClckTime      ) )//TRIM( TimePer )
               WRITE( UnSum, '(A)') ' Total CPU Time:        '//TRIM( Num2LStr( Factor*UsrTime       ) )//TRIM( TimePer )
               WRITE( UnSum, '(A)') ' Simulation CPU Time:   '//TRIM( Num2LStr( Factor*UsrTimeSim    ) )//TRIM( TimePer )
               WRITE( UnSum, '(A)') ' Simulated Time:        '//TRIM( Num2LStr( Factor*REAL( ZTime ) ) )//TRIM( TimePer )
               WRITE( UnSum, '(A)') ' Time Ratio (Sim/CPU):  '//TRIM( Num2LStr( TRatio ) )
            END IF
         END IF
            
         
         
      ENDIF

      IF (PRESENT(UsrTime_out)) UsrTime_out = UsrTime
      RETURN
   CONTAINS

      FUNCTION GetClockTime(StartClockTime, EndClockTime)
      ! return the number of seconds between StartClockTime and EndClockTime
   
         REAL                         :: GetClockTime          ! Elapsed clock time for the simulation phase of the run.
         INTEGER   , INTENT(IN)       :: StartClockTime (8)    ! Start time of simulation (after initialization)
         INTEGER   , INTENT(IN)       :: EndClockTime (8)      ! Start time of simulation (after initialization)
   
      !bjj: This calculation will be wrong at certain times (e.g. if it's near midnight on the last day of the month), but to my knowledge, no one has complained...
         GetClockTime =       0.001*( EndClockTime(8) - StartClockTime(8) ) &  ! Is the milliseconds of the second (range 0 to 999) - local time
                        +           ( EndClockTime(7) - StartClockTime(7) ) &  ! Is the seconds of the minute (range 0 to 59) - local time
                        +      60.0*( EndClockTime(6) - StartClockTime(6) ) &  ! Is the minutes of the hour (range 0 to 59) - local time
                        +    3600.0*( EndClockTime(5) - StartClockTime(5) ) &  ! Is the hour of the day (range 0 to 23) - local time
                        + SecPerDay*( EndClockTime(3) - StartClockTime(3) )    ! Is the day of the month
   
   
      END FUNCTION GetClockTime
   
   END SUBROUTINE RunTimes   
!=======================================================================   
!> this routine takes angles (in radians) and converts them to appropriate
!! ranges so they can be interpolated appropriately
!! (i.e., interpolating between pi+.1 and -pi should give pi+0.5 
!! instead of of 0.05 radians, so we return the angles pi+.1 and -pi+2pi=pi
!! we assume the interpolation occurs in the second dimension of angles
!! and it is done for each angle in the first dimension
   SUBROUTINE SetAnglesForInterp( angles )

   
      REAL(ReKi), INTENT(INOUT)     :: angles(:,:)  !

      REAL(ReKi)                    :: diff         ! difference between two adjacent angles 
      INTEGER(IntKi)                :: nr, nc       ! size of the angles matrix
      INTEGER(IntKi)                :: ir, ic       ! loop counters for each array dimension
   
      nr = size(angles,1)
      nc = size(angles,2)
   
   
         ! now let's make sure they don't cross a 2pi boundary (max |difference| can be pi):
         ! bjj: this is a dumb algorithm that should be revisited sometime
   
      do ic=2,nc            
         do ir=1,nr
            diff = angles(ir,ic-1) - angles(ir,ic)
            do while ( diff > pi )
               angles(ir,ic) = angles(ir,ic) + TwoPi
               diff = angles(ir,ic-1) - angles(ir,ic)
            end do
            do while ( diff < -pi )
               angles(ir,ic) = angles(ir,ic) - TwoPi
               diff = angles(ir,ic-1) - angles(ir,ic)
            end do                     
         end do      
      end do
   
   
   END SUBROUTINE SetAnglesForInterp
!=======================================================================
!> This routine computes numeric constants stored in the NWTC Library
   SUBROUTINE SetConstants( )

         ! Constants based upon Pi:

      Pi_D      = ACOS( -1.0_DbKi )
      D2R_D     = Pi_D/180.0_DbKi
      R2D_D     = 180.0_DbKi/Pi_D
      PiBy2_D   = Pi_D/2.0_DbKi
      RPM2RPS_D = Pi_D/30.0_DbKi
      RPS2RPM_D = 30.0_DbKi/Pi_D
      TwoByPi_D =  2.0_DbKi/Pi_D
      TwoPi_D   =  2.0_DbKi*Pi_D
      Inv2Pi_D  =  0.5_DbKi/Pi_D    ! 1.0_DbKi/TwoPi_D

      Pi      = ACOS( -1.0_ReKi )
      D2R     = Pi/180.0_ReKi
      R2D     = 180.0_ReKi/Pi
      PiBy2   = Pi/2.0_ReKi
      RPM2RPS = Pi/30.0_ReKi
      RPS2RPM = 30.0_ReKi/Pi
      TwoByPi =  2.0_ReKi/Pi
      TwoPi   =  2.0_ReKi*Pi
      Inv2Pi  =  0.5_ReKi/Pi        ! 1.0/TwoPi

      Pi_S      = ACOS( -1.0_SiKi )
      D2R_S     = Pi_S/180.0_SiKi
      R2D_S     = 180.0_SiKi/Pi_S
      PiBy2_S   = Pi_S/2.0_SiKi
      RPM2RPS_S = Pi_S/30.0_SiKi
      RPS2RPM_S = 30.0_SiKi/Pi_S
      TwoByPi_S =  2.0_SiKi/Pi_S
      TwoPi_S   =  2.0_SiKi*Pi_S
      Inv2Pi_S  =  0.5_SiKi/Pi_S    ! 1.0_SiKi/TwoPi_S
      Pi_R4   = ACOS( -1.0_SiKi )
      Pi_R8   = ACOS( -1.0_R8Ki )

      TwoPi_R4  = Pi_R4 *2.0_SiKi
      TwoPi_R8  = Pi_R8 *2.0_R8Ki
      
         ! IEEE constants:
      CALL Set_IEEE_Constants( NaN_D, Inf_D, NaN, Inf, NaN_S, Inf_S )
      

   RETURN
   END SUBROUTINE SetConstants
!=======================================================================   
!> This routine displays a message that gives that status of the simulation.
!! It is intended to be used with RunTimes (nwtc_num::runtimes) and SimStatus (nwtc_num::simstatus).
   SUBROUTINE SimStatus_FirstTime( PrevSimTime, PrevClockTime, SimStrtTime, UsrTimeSim, ZTime, TMax, DescStrIn )

      IMPLICIT                        NONE

         ! Passed variables
      REAL(DbKi), INTENT(IN   )    :: ZTime                                           !< Current simulation time (s)
      REAL(DbKi), INTENT(IN   )    :: TMax                                            !< Expected simulation time (s)
      REAL(DbKi), INTENT(  OUT)    :: PrevSimTime                                     !< Previous time message was written to screen (s > 0)
      REAL(ReKi), INTENT(  OUT)    :: PrevClockTime                                   !< Previous clock time in seconds past midnight
      INTEGER,    INTENT(  OUT)    :: SimStrtTime (8)                                 !< An array containing the elements of the start time.
      REAL(ReKi), INTENT(  OUT)    :: UsrTimeSim                                      !< User CPU time for simulation (without intialization)

      CHARACTER(*), INTENT(IN), OPTIONAL :: DescStrIn                                 !< optional additional string to print for SimStatus
      
         ! Local variables.

      REAL(ReKi)                   :: CurrClockTime                                   ! Current time in seconds past midnight.
      CHARACTER(10)                :: DescStr                                        !< optional additional string to print for SimStatus


      if (present(DescStrIn)) then
         DescStr = DescStrIn
      else
         DescStr = ""
      end if


         ! How many seconds past midnight?

      CALL DATE_AND_TIME ( Values=SimStrtTime )
      CALL CPU_TIME ( UsrTimeSim )                                                    ! Initial CPU time   
      UsrTimeSim = MAX( 0.0_ReKi, UsrTimeSim )  ! CPU_TIME: If a meaningful time cannot be returned, a processor-dependent negative value is returned

      CurrClockTime = TimeValues2Seconds( SimStrtTime )


      CALL WrScr ( trim(DescStr)//' Time: '//TRIM( Num2LStr( NINT( ZTime ) ) )//' of '//TRIM( Num2LStr( TMax ) )//' seconds.')


      ! Let's save this time as the previous time for the next call to the routine
      PrevClockTime = CurrClockTime
      PrevSimTime   = ZTime

      RETURN
   END SUBROUTINE SimStatus_FirstTime
!=======================================================================
!> This routine displays a message that gives that status of the simulation and the predicted end time of day.
!! It is intended to be used with RunTimes (nwtc_num::runtimes) and SimStatus_FirstTime (nwtc_num::simstatus_firsttime).
   SUBROUTINE SimStatus( PrevSimTime, PrevClockTime, ZTime, TMax, DescStrIn, StatInfoIn)
   

      IMPLICIT                        NONE

         ! Passed variables
      REAL(DbKi), INTENT(IN)       :: ZTime                                !< Current simulation time (s)
      REAL(DbKi), INTENT(IN)       :: TMax                                 !< Expected simulation time (s)
      REAL(DbKi), INTENT(INOUT)    :: PrevSimTime                          !< Previous time message was written to screen (s > 0)
      REAL(ReKi), INTENT(INOUT)    :: PrevClockTime                        !< Previous clock time in seconds past midnight

      CHARACTER(*), INTENT(IN), OPTIONAL :: DescStrIn                      !< optional additional string to print at start of SimStatus
      CHARACTER(*), INTENT(IN), OPTIONAL :: StatInfoIn                     !< optional additional string to print at end of SimStatus

         ! Local variables.

      REAL(ReKi)                   :: CurrClockTime                        ! Current time in seconds past midnight.
      REAL(ReKi)                   :: DeltTime                             ! The amount of time elapsed since the last call.
      REAL(ReKi)                   :: EndTime                              ! Approximate time of day when simulation will complete.
      REAL(ReKi), PARAMETER        :: InSecHr  = 1.0_ReKi/3600.0_ReKi      ! Inverse of the number of seconds in an hour
      REAL(ReKi), PARAMETER        :: InSecMn  = 1.0_ReKi/  60.0_ReKi      ! Inverse of the number of seconds in a minute
      REAL(ReKi)                   :: SimTimeLeft                          ! Approximate clock time remaining before simulation completes

      REAL(ReKi), PARAMETER        :: SecPerDay = 24*60*60.0_ReKi          ! Number of seconds per day

      INTEGER(4)                   :: EndHour                              ! The hour when the simulations is expected to complete.
      INTEGER(4)                   :: EndMin                               ! The minute when the simulations is expected to complete.
      INTEGER(4)                   :: EndSec                               ! The second when the simulations is expected to complete.
      INTEGER(4)                   :: TimeAry  (8)                         ! An array containing the elements of the start time.

      CHARACTER(MaxWrScrLen)       :: BlankLine
      CHARACTER( 8)                :: ETimeStr                             ! String containing the end time.
      CHARACTER( 10)               :: DescStr                              !< optional additional string to print for SimStatus
      CHARACTER(200)               :: StatInfo                             !< optional additional string to print for SimStatus


      IF ( ZTime <= PrevSimTime ) RETURN


      if (present(DescStrIn)) then
         DescStr = DescStrIn
      else
         DescStr = ""
      end if

      if (present(StatInfoIn)) then
         StatInfo = StatInfoIn
      else
         StatInfo = ""
      end if
      
         ! How many seconds past midnight?

      CALL DATE_AND_TIME ( Values=TimeAry )
      CurrClockTime = TimeValues2Seconds( TimeAry )

         ! Calculate elapsed clock time

      DeltTime = CurrClockTime - PrevClockTime


         ! We may have passed midnight since the last revoultion.  We will assume that (ZTime - PrevSimTime) of 
        !  simulation time doesn't take more than a day.

      IF ( CurrClockTime < PrevClockTime )  THEN
         DeltTime = DeltTime + SecPerDay
      ENDIF


         ! Estimate the end time in hours, minutes, and seconds

      SimTimeLeft = REAL( ( TMax - ZTime )*DeltTime/( ZTime - PrevSimTime ), ReKi )          ! DeltTime/( ZTime - PrevSimTime ) is the delta_ClockTime divided by the delta_SimulationTime
      EndTime  =  MOD( CurrClockTime+SimTimeLeft, SecPerDay )
      EndHour  =  INT(   EndTime*InSecHr )
      EndMin   =  INT( ( EndTime - REAL( 3600*EndHour ) )*InSecMn )
      EndSec   = NINT(   EndTime - REAL( 3600*EndHour + 60*EndMin ) ) !bjj: this NINT can make the seconds say "60"

      WRITE (ETimeStr,"(I2.2,2(':',I2.2))")  EndHour, EndMin, EndSec

      BlankLine = ""
      CALL WrOver( BlankLine )  ! BlankLine contains MaxWrScrLen spaces
      CALL WrOver ( trim(DescStr)//' Time: '//TRIM( Num2LStr( NINT( ZTime ) ) )//' of '//TRIM( Num2LStr( TMax ) )// &
                    ' seconds. '//trim(StatInfo)// &
                    ' Estimated final completion at '//ETimeStr//'.')

         ! Let's save this time as the previous time for the next call to the routine
      PrevClockTime = CurrClockTime
      PrevSimTime   = ZTime

      RETURN
   END SUBROUTINE SimStatus
!=======================================================================
!>  This routine computes the 3x3 transformation matrix, \f$TransMat\f$,
!!   to a coordinate system \f$x\f$ (with orthogonal axes \f$x_1, x_2, x_3\f$)
!!   resulting from three rotations (\f$\theta_1\f$, \f$\theta_2\f$, \f$\theta_3\f$) about the
!!   orthogonal axes (\f$X_1, X_2, X_3\f$) of coordinate system \f$X\f$.  All angles
!!   are assummed to be small, as such, the order of rotations does
!!   not matter and Euler angles do not need to be used.  This routine
!!   is used to compute the transformation matrix (\f$TransMat\f$) between
!!   undeflected (\f$X\f$) and deflected (\f$x\f$) coordinate systems.  In matrix
!!   form:
!! \f{equation}{   
!!   \left\{ \begin{matrix} x_1 \\ x_2 \\ x_3 \end{matrix} \right\} =
!!   \left[ TransMat(\theta_1, \theta_2, \theta_3) \right]
!!   \left\{ \begin{matrix} X_1 \\ X_2 \\ X_3 \end{matrix} \right\}   
!! \f}
!!
!! The transformation matrix, \f$TransMat\f$, is the closest orthonormal
!!   matrix to the nonorthonormal, but skew-symmetric, Bernoulli-Euler
!!   matrix:
!! \f{equation}{   A =
!!   \begin{bmatrix}    1      &  \theta_3 & -\theta_2 \\
!!                   -\theta_3 &  1        &  \theta_1 \\
!!                    \theta_2 & -\theta_1 &  1 \end{bmatrix}   
!! \f}
!!   In the Frobenius Norm sense, the closest orthornormal matrix is:
!!      \f$TransMat = U V^T\f$,
!!   where the columns of \f$U\f$ contain the eigenvectors of\f$ AA^T\f$ and the
!!   columns of \f$V\f$ contain the eigenvectors of \f$A^TA\f$ (note that \f$^T\f$ = transpose).
!!   This result comes directly from the Singular Value Decomposition
!!   (SVD) of \f$A = USV^T\f$ where \f$S\f$ is a diagonal matrix containing the
!!   singular values of \f$A\f$, which are sqrt( eigenvalues of \f$AA^T\f$ ) =
!!   sqrt( eigenvalues of \f$A^TA\f$ ).
!!
!! The algebraic form of the transformation matrix, as implemented
!!   below, was derived symbolically by J. Jonkman by computing \f$UV^T\f$
!!   by hand with verification in Mathematica.
!!
!! This routine is the inverse of GetSmllRotAngs (nwtc_num::getsmllrotangs). \n
!! Use SmllRotTrans (nwtc_num::smllrottrans) instead of directly calling a specific routine in the generic interface. 
   SUBROUTINE SmllRotTransD( RotationType, Theta1, Theta2, Theta3, TransMat, ErrTxt, ErrStat, ErrMsg )

      ! Passed Variables:

   REAL(R8Ki), INTENT(IN )             :: Theta1                                          !< \f$\theta_1\f$: the small rotation about \f$X_1\f$, (rad).
   REAL(R8Ki), INTENT(IN )             :: Theta2                                          !< \f$\theta_2\f$: the small rotation about \f$X_2\f$, (rad).
   REAL(R8Ki), INTENT(IN )             :: Theta3                                          !< \f$\theta_3\f$: the small rotation about \f$X_3\f$, (rad).
   REAL(R8Ki), INTENT(OUT)             :: TransMat (3,3)                                  !< The resulting transformation matrix from \f$X\f$ to \f$x\f$, (-).

   INTEGER(IntKi),INTENT(OUT)          :: ErrStat                                         !< Error status 
   CHARACTER(*), INTENT(OUT)           :: ErrMsg                                          !< Error message corresponding to ErrStat

   CHARACTER(*), INTENT(IN)            :: RotationType                                    !< The type of rotation; used to inform the user where a large rotation is occuring upon such an event.
   CHARACTER(*), INTENT(IN ), OPTIONAL :: ErrTxt                                          !< an additional message to be displayed as a warning (typically the simulation time)

      ! Local Variables:

   REAL(DbKi)                          :: ComDenom                                        ! = ( Theta1^2 + Theta2^2 + Theta3^2 )*SQRT( 1.0 + Theta1^2 + Theta2^2 + Theta3^2 )
   REAL(DbKi), PARAMETER               :: LrgAngle  = 0.4                                 ! Threshold for when a small angle becomes large (about 23deg).  This comes from: COS(SmllAngle) ~ 1/SQRT( 1 + SmllAngle^2 ) and SIN(SmllAngle) ~ SmllAngle/SQRT( 1 + SmllAngle^2 ) results in ~5% error when SmllAngle = 0.4rad.
   REAL(DbKi)                          :: Theta11                                         ! = Theta1^2
   REAL(DbKi)                          :: Theta12S                                        ! = Theta1*Theta2*[ SQRT( 1.0 + Theta1^2 + Theta2^2 + Theta3^2 ) - 1.0 ]
   REAL(DbKi)                          :: Theta13S                                        ! = Theta1*Theta3*[ SQRT( 1.0 + Theta1^2 + Theta2^2 + Theta3^2 ) - 1.0 ]
   REAL(DbKi)                          :: Theta22                                         ! = Theta2^2
   REAL(DbKi)                          :: Theta23S                                        ! = Theta2*Theta3*[ SQRT( 1.0 + Theta1^2 + Theta2^2 + Theta3^2 ) - 1.0 ]
   REAL(DbKi)                          :: Theta33                                         ! = Theta3^2
   REAL(DbKi)                          :: SqrdSum                                         ! = Theta1^2 + Theta2^2 + Theta3^2
   REAL(DbKi)                          :: SQRT1SqrdSum                                    ! = SQRT( 1.0 + Theta1^2 + Theta2^2 + Theta3^2 )

   LOGICAL,    SAVE                    :: FrstWarn  = .TRUE.                              ! When .TRUE., indicates that we're on the first warning.


   ErrStat = ErrID_None
   ErrMsg  = ''

      ! Display a warning message if at least one angle gets too large in magnitude:

   IF ( ( ( ABS(Theta1) > LrgAngle ) .OR. ( ABS(Theta2) > LrgAngle ) .OR. ( ABS(Theta3) > LrgAngle ) ) .AND. FrstWarn )  THEN

      ErrStat= ErrID_Severe
      ErrMsg = 'Small angle assumption violated in SUBROUTINE SmllRotTrans() due to a large '//TRIM(RotationType)//'. '// &
               'The solution may be inaccurate. Simulation continuing, but future warnings from SmllRotTrans() will be suppressed.'

      IF ( PRESENT(ErrTxt) ) THEN
         ErrMsg = TRIM(ErrMsg)//NewLine//' Additional debugging message from SUBROUTINE SmllRotTrans(): '//TRIM(ErrTxt)
      END IF

      !CALL ProgWarn( TRIM(ErrMsg) )

      FrstWarn = .FALSE.   ! Don't enter here again!

   ENDIF


      ! Compute some intermediate results:

   Theta11      = Theta1*Theta1
   Theta22      = Theta2*Theta2
   Theta33      = Theta3*Theta3

   SqrdSum      = Theta11 + Theta22 + Theta33
   SQRT1SqrdSum = SQRT( 1.0_DbKi + SqrdSum )
   ComDenom     = SqrdSum*SQRT1SqrdSum

   Theta12S     = Theta1*Theta2*( SQRT1SqrdSum - 1.0_DbKi )
   Theta13S     = Theta1*Theta3*( SQRT1SqrdSum - 1.0_DbKi )
   Theta23S     = Theta2*Theta3*( SQRT1SqrdSum - 1.0_DbKi )


      ! Define the transformation matrix:

   IF ( ComDenom == 0.0_DbKi )  THEN  ! All angles are zero and matrix is ill-conditioned (the matrix is derived assuming that the angles are not zero); return identity

      TransMat(1,:) = (/ 1.0_DbKi, 0.0_DbKi, 0.0_DbKi /)
      TransMat(2,:) = (/ 0.0_DbKi, 1.0_DbKi, 0.0_DbKi /)
      TransMat(3,:) = (/ 0.0_DbKi, 0.0_DbKi, 1.0_DbKi /)

   ELSE                          ! At least one angle is nonzero

      TransMat(1,1) = ( Theta11*SQRT1SqrdSum + Theta22              + Theta33              )/ComDenom
      TransMat(2,2) = ( Theta11              + Theta22*SQRT1SqrdSum + Theta33              )/ComDenom
      TransMat(3,3) = ( Theta11              + Theta22              + Theta33*SQRT1SqrdSum )/ComDenom
      TransMat(1,2) = (  Theta3*SqrdSum + Theta12S )/ComDenom
      TransMat(2,1) = ( -Theta3*SqrdSum + Theta12S )/ComDenom
      TransMat(1,3) = ( -Theta2*SqrdSum + Theta13S )/ComDenom
      TransMat(3,1) = (  Theta2*SqrdSum + Theta13S )/ComDenom
      TransMat(2,3) = (  Theta1*SqrdSum + Theta23S )/ComDenom
      TransMat(3,2) = ( -Theta1*SqrdSum + Theta23S )/ComDenom

   ENDIF


   RETURN
   END SUBROUTINE SmllRotTransD
!=======================================================================
!> \copydoc nwtc_num::smllrottransd
   SUBROUTINE SmllRotTransR( RotationType, Theta1, Theta2, Theta3, TransMat, ErrTxt, ErrStat, ErrMsg )


      ! This routine computes the 3x3 transformation matrix, TransMat,
      !   to a coordinate system x (with orthogonal axes x1, x2, x3)
      !   resulting from three rotations (Theta1, Theta2, Theta3) about the
      !   orthogonal axes (X1, X2, X3) of coordinate system X.  All angles
      !   are assummed to be small, as such, the order of rotations does
      !   not matter and Euler angles do not need to be used.  This routine
      !   is used to compute the transformation matrix (TransMat) between
      !   undeflected (X) and deflected (x) coordinate systems.  In matrix
      !   form:
      !      {x1}   [TransMat(Theta1, ] {X1}
      !      {x2} = [         Theta2, ]*{X2}
      !      {x3}   [         Theta3 )] {X3}
      !
      ! The transformation matrix, TransMat, is the closest orthonormal
      !   matrix to the nonorthonormal, but skew-symmetric, Bernoulli-Euler
      !   matrix:
      !          [   1.0    Theta3 -Theta2 ]
      !      A = [ -Theta3   1.0    Theta1 ]
      !          [  Theta2 -Theta1   1.0   ]
      !
      !   In the Frobenius Norm sense, the closest orthornormal matrix is:
      !      TransMat = U*V^T,
      !
      !   where the columns of U contain the eigenvectors of A*A^T and the
      !   columns of V contain the eigenvectors of A^T*A (^T = transpose).
      !   This result comes directly from the Singular Value Decomposition
      !   (SVD) of A = U*S*V^T where S is a diagonal matrix containing the
      !   singular values of A, which are SQRT( eigenvalues of A*A^T ) =
      !   SQRT( eigenvalues of A^T*A ).
      !
      ! The algebraic form of the transformation matrix, as implemented
      !   below, was derived symbolically by J. Jonkman by computing U*V^T
      !   by hand with verification in Mathematica.
      !
      ! This routine is the inverse of GetSmllRotAngs()

      ! Passed Variables:

   REAL(SiKi), INTENT(IN )             :: Theta1                                          ! The small rotation about X1, (rad).
   REAL(SiKi), INTENT(IN )             :: Theta2                                          ! The small rotation about X2, (rad).
   REAL(SiKi), INTENT(IN )             :: Theta3                                          ! The small rotation about X3, (rad).
   REAL(SiKi), INTENT(OUT)             :: TransMat (3,3)                                  ! The resulting transformation matrix from X to x, (-).

   INTEGER(IntKi),INTENT(OUT)          :: ErrStat
   CHARACTER(*), INTENT(OUT)           :: ErrMsg

   CHARACTER(*), INTENT(IN)            :: RotationType                                    ! The type of rotation; used to inform the user where a large rotation is occuring upon such an event.
   CHARACTER(*), INTENT(IN ), OPTIONAL :: ErrTxt                                          ! an additional message to be displayed as a warning (typically the simulation time)

      ! Local Variables:

   REAL(SiKi)                          :: ComDenom                                        ! = ( Theta1^2 + Theta2^2 + Theta3^2 )*SQRT( 1.0 + Theta1^2 + Theta2^2 + Theta3^2 )
   REAL(SiKi), PARAMETER               :: LrgAngle  = 0.4                                 ! Threshold for when a small angle becomes large (about 23deg).  This comes from: COS(SmllAngle) ~ 1/SQRT( 1 + SmllAngle^2 ) and SIN(SmllAngle) ~ SmllAngle/SQRT( 1 + SmllAngle^2 ) results in ~5% error when SmllAngle = 0.4rad.
   REAL(SiKi)                          :: Theta11                                         ! = Theta1^2
   REAL(SiKi)                          :: Theta12S                                        ! = Theta1*Theta2*[ SQRT( 1.0 + Theta1^2 + Theta2^2 + Theta3^2 ) - 1.0 ]
   REAL(SiKi)                          :: Theta13S                                        ! = Theta1*Theta3*[ SQRT( 1.0 + Theta1^2 + Theta2^2 + Theta3^2 ) - 1.0 ]
   REAL(SiKi)                          :: Theta22                                         ! = Theta2^2
   REAL(SiKi)                          :: Theta23S                                        ! = Theta2*Theta3*[ SQRT( 1.0 + Theta1^2 + Theta2^2 + Theta3^2 ) - 1.0 ]
   REAL(SiKi)                          :: Theta33                                         ! = Theta3^2
   REAL(SiKi)                          :: SqrdSum                                         ! = Theta1^2 + Theta2^2 + Theta3^2
   REAL(SiKi)                          :: SQRT1SqrdSum                                    ! = SQRT( 1.0 + Theta1^2 + Theta2^2 + Theta3^2 )

   LOGICAL,    SAVE                    :: FrstWarn  = .TRUE.                              ! When .TRUE., indicates that we're on the first warning.


   ErrStat = ErrID_None
   ErrMsg  = ''

      ! Display a warning message if at least one angle gets too large in magnitude:

   IF ( ( ( ABS(Theta1) > LrgAngle ) .OR. ( ABS(Theta2) > LrgAngle ) .OR. ( ABS(Theta3) > LrgAngle ) ) .AND. FrstWarn )  THEN

      ErrStat= ErrID_Severe
      ErrMsg = 'Small angle assumption violated in SUBROUTINE SmllRotTrans() due to a large '//TRIM(RotationType)//'. '// &
               'The solution may be inaccurate. Simulation continuing, but future warnings from SmllRotTrans() will be suppressed.'

      IF ( PRESENT(ErrTxt) ) THEN
         ErrMsg = TRIM(ErrMsg)//NewLine//' Additional debugging message from SUBROUTINE SmllRotTrans(): '//TRIM(ErrTxt)
      END IF

      !CALL ProgWarn( TRIM(ErrMsg) )

      FrstWarn = .FALSE.   ! Don't enter here again!

   ENDIF


      ! Compute some intermediate results:

   Theta11      = Theta1*Theta1
   Theta22      = Theta2*Theta2
   Theta33      = Theta3*Theta3

   SqrdSum      = Theta11 + Theta22 + Theta33
   SQRT1SqrdSum = SQRT( 1.0_ReKi + SqrdSum )
   ComDenom     = SqrdSum*SQRT1SqrdSum

   Theta12S     = Theta1*Theta2*( SQRT1SqrdSum - 1.0_Siki )
   Theta13S     = Theta1*Theta3*( SQRT1SqrdSum - 1.0_Siki )
   Theta23S     = Theta2*Theta3*( SQRT1SqrdSum - 1.0_Siki )


      ! Define the transformation matrix:

   IF ( ComDenom == 0.0_ReKi )  THEN  ! All angles are zero and matrix is ill-conditioned (the matrix is derived assuming that the angles are not zero); return identity

      TransMat(1,:) = (/ 1.0_SiKi, 0.0_SiKi, 0.0_SiKi /)
      TransMat(2,:) = (/ 0.0_SiKi, 1.0_SiKi, 0.0_SiKi /)
      TransMat(3,:) = (/ 0.0_SiKi, 0.0_SiKi, 1.0_SiKi /)

   ELSE                          ! At least one angle is nonzero

      TransMat(1,1) = ( Theta11*SQRT1SqrdSum + Theta22              + Theta33              )/ComDenom
      TransMat(2,2) = ( Theta11              + Theta22*SQRT1SqrdSum + Theta33              )/ComDenom
      TransMat(3,3) = ( Theta11              + Theta22              + Theta33*SQRT1SqrdSum )/ComDenom
      TransMat(1,2) = (  Theta3*SqrdSum + Theta12S )/ComDenom
      TransMat(2,1) = ( -Theta3*SqrdSum + Theta12S )/ComDenom
      TransMat(1,3) = ( -Theta2*SqrdSum + Theta13S )/ComDenom
      TransMat(3,1) = (  Theta2*SqrdSum + Theta13S )/ComDenom
      TransMat(2,3) = (  Theta1*SqrdSum + Theta23S )/ComDenom
      TransMat(3,2) = ( -Theta1*SqrdSum + Theta23S )/ComDenom

   ENDIF


   RETURN
   END SUBROUTINE SmllRotTransR
!=======================================================================
!> This routine takes two sorted arrays and finds the sorted union of the two.
!!
!! Note: If the same value is found in both arrays, only one is kept. However, if either
!!       array as multiple occurances of the same value, the largest multiple will be
!!       kept. Duplicates should be eliminated externally if this is not desirable.
   SUBROUTINE SortUnion ( Ary1, N1, Ary2, N2, Ary, N )

      ! Argument declarations:

   INTEGER, INTENT(OUT)         :: N                                            !< The length of the output array.
   INTEGER, INTENT(IN)          :: N1                                           !< The length of the first input array.
   INTEGER, INTENT(IN)          :: N2                                           !< The length of the second input array.

   REAL(ReKi), INTENT(OUT)      :: Ary(N1+N2)                                   !< The sorted union.
   REAL(ReKi), INTENT(IN)       :: Ary1(N1)                                     !< The first list of sorted real numbers.
   REAL(ReKi), INTENT(IN)       :: Ary2(N2)                                     !< The second list of sorted real numbers.


      ! Local declarations:

   INTEGER                      :: I1                                           ! Index into the first array.
   INTEGER                      :: I2                                           ! Index into the second array.



   I1 = 1
   I2 = 1
   N  = 1

   DO WHILE ( ( I1 <= N1 ) .AND. ( I2 <= N2 ) )

      IF ( Ary1(I1) < Ary2(I2) )  THEN
         Ary(N) = Ary1(I1)
         I1 = I1 + 1
      ELSE IF ( Ary1(I1) > Ary2(I2) )  THEN
         Ary(N) = Ary2(I2)
         I2 = I2 + 1
      ELSE
         Ary(N) = Ary1(I1)
         I1 = I1 + 1
         I2 = I2 + 1
      END IF

      N  = N  + 1

   END DO ! WHILE


     ! We've reached the end of one array, but we need to add the end
     ! of the other array if we haven't reached the end of it yet.

   IF ( I1 <= N1 ) THEN
      Ary(N:N+N1-I1) = Ary1(I1:)
      N = N+N1-I1
   ELSEIF ( I2 <= N2 ) THEN
      Ary(N:N+N2-I2) = Ary2(I2:)
      N = N+N2-I2
   ELSE
      N = N - 1
   ENDIF


   RETURN
   END SUBROUTINE SortUnion ! ( Ary1, N1, Ary2, N2, Ary, N )
!=======================================================================
!> This routine calculates the standard deviation of a population contained in Ary.
!!
!! This can be calculated as either\n
!! \f$ \sqrt{ \frac{\sum_{i=1}^N \left(x_i -\bar{x}\right)^2 }{N-1} } \f$   \n
!! or \n
!! \f$ \sqrt{ \frac{\sum_{i=1}^N \left(x_i -\bar{x}\right)^2 }{N} } \f$ if `UseN` is true \n
   FUNCTION StdDevFn ( Ary, AryLen, Mean, UseN )

      ! Function declaration.

   REAL(ReKi)                   :: StdDevFn                                     !< This function.


      ! Argument declarations:

   INTEGER, INTENT(IN)          :: AryLen                                       !< Length of the array.

   REAL(ReKi), INTENT(IN)       :: Ary  (AryLen)                                !< Input array.
   REAL(ReKi), INTENT(IN)       :: Mean                                         !< The previously calculated mean of the array.
   LOGICAL, OPTIONAL, INTENT(IN) :: UseN                                        !< Use `N` insted of `N-1` in denomenator


      ! Local declarations.

   REAL(DbKi)                   :: Sum                                          ! A temporary sum.

   INTEGER                      :: I                                            ! The index into the array.
   INTEGER                      :: Denom                                        ! Denominator

   IF(PRESENT(UseN)) THEN
      IF (UseN) THEN
         Denom = AryLen
      ELSE
         Denom = AryLen-1
      ENDIF
   ELSE
      Denom = AryLen-1
   ENDIF

   Sum = 0.0_DbKi

   DO I=1,AryLen
      Sum = Sum + ( Ary(I) - Mean )**2
   END DO ! I

   StdDevFn = SQRT( Sum/( Denom ) )


   RETURN
   END FUNCTION StdDevFn
!=======================================================================
!> This function returns the 3x3 skew-symmetric matrix for cross-product
!! calculation of vector \f$\vec{x}\f$ via matrix multiplication, defined as
!! \f{equation}{
!!   f_{_\times}\left( \vec{x} \right) = 
!!  \begin{bmatrix}
!!       0  & -x_3 &  x_2 \\
!!      x_3 &  0   & -x_1 \\
!!     -x_2 &  x_1 &  0          
!! \end{bmatrix}
!! \f}   
!> Use SkewSymMat (nwtc_num::skewsymmat) instead of directly calling a specific routine in the generic interface.
   FUNCTION SkewSymMatR4 ( x ) RESULT(M)

      ! Function arguments

   REAL(SiKi)                   :: M(3,3)                          !< skew-symmetric matrix formed from input vector \f$x\f$
   REAL(SiKi), INTENT(IN)       :: x(3)                            !< input vector \f$x\f$

   M(1,1) =    0.0_SiKi
   M(2,1) =  x(3)
   M(3,1) = -x(2)

   M(1,2) = -x(3)
   M(2,2) =    0.0_SiKi
   M(3,2) =  x(1)

   M(1,3) =  x(2)
   M(2,3) = -x(1)
   M(3,3) =    0.0_SiKi
   
   RETURN
   END FUNCTION SkewSymMatR4 
!=======================================================================
!> \copydoc nwtc_num::skewsymmatr4
   FUNCTION SkewSymMatR8 ( x ) RESULT(M)

      ! Function arguments

   REAL(R8Ki)                   :: M(3,3)                          ! skew-symmetric matrix formed from input vector \f$x\f$
   REAL(R8Ki), INTENT(IN)       :: x(3)                            ! input vector \f$x\f$

   M(1,1) =    0.0_R8Ki
   M(2,1) =  x(3)
   M(3,1) = -x(2)

   M(1,2) = -x(3)
   M(2,2) =    0.0_R8Ki
   M(3,2) =  x(1)

   M(1,3) =  x(2)
   M(2,3) = -x(1)
   M(3,3) =    0.0_R8Ki
   
   RETURN
   END FUNCTION SkewSymMatR8
!=======================================================================
!> If M is a rotation matrix from a 1-2-3 rotation sequence about Y-X-Z, this function returns 
!! the 3 sequential angles, \f$\theta_y\f$, \f$\theta_x\f$, and \f$\theta_z\f$ (in radians), that formed 
!! the matrix. M represents a change of basis (from global to local coordinates; 
!! not a physical rotation of the body; passive rotation).
!!
!! See Tait-Bryan angle \f$ Y_1 X_2 Z_3 \f$ at https://en.wikipedia.org/wiki/Euler_angles#Rotation_matrix
!! Note that what we are using here is the passive rotation, which is the transpose of what appears in the
!! wikipedia article.
!! 
!!
!! \f{eqnarray*}{   
!! M & = & R(\theta_z) R(\theta_x) R(\theta_y)
!!   & = & R(\theta_3) R(\theta_2) R(\theta_1) \\
!!   & = & \begin{bmatrix}    \cos(\theta_z)    & \sin(\theta_z)     & 0                  \\
!!                            -\sin(\theta_z)   & \cos(\theta_z)     & 0                  \\
!!                            0                 &  0                 & 1                  \end{bmatrix}
!!         \begin{bmatrix}    1                 &  0                 & 0                  \\
!!                            0                 &  \cos(\theta_x)    & \sin(\theta_x)     \\
!!                            0                 & -\sin(\theta_x)    & \cos(\theta_x)     \end{bmatrix}
!!         \begin{bmatrix}    \cos(\theta_y)    & 0                  & -\sin(\theta_y)    \\
!!                            0                 & 1                  & 0                  \\
!!                            \sin(\theta_y)    & 0                  & \cos(\theta_y)     \end{bmatrix}
!!   & = & \begin{bmatrix}    C_3               & S_3                & 0                  \\
!!                            -S_3              & C_3                & 0                  \\
!!                            0                 & 0                  & 1                  \end{bmatrix}
!!         \begin{bmatrix}    1                 & 0                  & 0                  \\
!!                            0                 & C_2                & S_2                \\
!!                            0                 & -S_2               & C_2                \end{bmatrix}
!!         \begin{bmatrix}    C_1               & 0                  & -S_1               \\
!!                            0                 & 1                  & 0                  \\
!!                            S_1               & 0                  & C_1                \end{bmatrix}  \\
!!   & = & \begin{bmatrix}  
!!             \cos(\theta_y) \cos(\theta_z) + \sin(\theta_y) \sin(\theta_x) \sin(\theta_z)     &  \cos(\theta_x) \sin(\theta_z)    &  \cos(\theta_y) \sin(\theta_x) \sin(\theta_z) - \sin(\theta_y) \cos(\theta_z)  \\
!!             \sin(\theta_y) \sin(\theta_x) \cos(\theta_z) - \cos(\theta_y) \sin(\theta_z)     &  \cos(\theta_x) \cos(\theta_z)    &  \cos(\theta_y) \sin(\theta_x) \cos(\theta_z) + \sin(\theta_y) \sin(\theta_z)  \\
!!             \sin(\theta_y) \cos(\theta_x)                                                    &  -\sin(\theta_x)                  &  \cos(\theta_y) \cos(\theta_x)                                                 \\
!!         \end{bmatrix}   
!!   & = & \begin{bmatrix}  
!!             C_1 C_3 + S_1 S_2 S_3   &  C_2 S_3     &  C_1 S_2 S_3 - S_1 C_3   \\
!!             S_1 S_2 C_3 - C_1 S_3   &  C_2 C_3     &  C_1 S_2 C_3 + S_1 S_3   \\
!!             S_1 C_2                 &  -S_2        &  C_1 C_2                 \\
!!         \end{bmatrix}   
!! \f}
!! returned angles are in the range \f$\theta_y,\theta_x, \theta_z \in \left[ \pi, -\pi \right]\f$ \n
!! Use TaitBryanYXZExtract (nwtc_num::taitbryanyxzextract)  instead of directly calling a specific routine in the generic interface. 
   FUNCTION TaitBryanYXZExtractR4(M) result(theta)
   
   
      REAL(SiKi), INTENT(IN) :: M(3,3)    !< rotation matrix, M 
      REAL(SiKi)             :: theta(3)  !< the 3 rotation angles, \f$(\theta_y, \theta_x, \theta_z)\f$, corresponding to the Tait-Bryan rotation angle corresponding to cant-toe-twist
     
      REAL(SiKi)              :: C_1       ! C_1 > cos(theta_y) 
      REAL(SiKi)              :: S_1       ! S_1 > sin(theta_y) 
      REAL(SiKi)              :: C_2       ! C_2 > cos(theta_x) 
      REAL(SiKi)              :: C_3       ! C_3 > cos(theta_z) 
      REAL(SiKi)              :: S_3       ! S_3 > sin(theta_z)

         !> # Algorithm

         !> Starting with the trig identity of \f$ \sin(\theta_3)^2 + \cos(\theta_3)^2 = S_3^2 + C_3^2 \equiv 1\f$, we can find \f$ \cos(\theta_2) \f$
         !! from matrix elements \f$M(1,2)\f$ and \f$M(2,2)\f$ by 
         !! \f{equation}{
         !!    \cos(\theta_2) = C_2 = \sqrt{ M(1,2)^2 + M(2,2)^2} = \sqrt{ C_2^2 S_3^2 + C_2^2 C_3^2 } = \sqrt{ C_2^2 ( S_3^2 + C_3^2 ) }.
         !! \f}

         ! use trig identity S_3**2 + C_3**2 = 1 to get abs( C_2 )
      C_2 = sqrt( m(1,2)**2 + m(2,2)**2 )

      if ( EqualRealNos( C_2, 0.0_SiKi ) ) then

         !> ## If \f$ \cos(\theta_2) = C_2 = 0\f$:
         !! If \f$\cos(\theta_2) = C_2 = 0\f$, then \f$ \theta_2 \f$ is \f$ \pm\pi/2 \f$ and \f$ S_2 = \pm 1\f$.  We can solve for the sign of \f$\theta_2\f$ by using 
         !! \f{equation}{
         !!    \theta_2 = \arctan{\left( \frac{-M(3,2)}{C_2} \right)} = \arctan{\left( \frac{S_2}{C_2} \right)}
         !! \f}
         !! (but using the _atan2_ function in the complex plane instead of \f$ \arctan \f$).

         theta(2) = atan2( -m(3,2), C_2 )     ! theta_2 -> theta_x


         !> Considering \f$ C_2 = 0 \f$ and \f$ S_2 = \pm 1\f$, the matrix \f$ M \f$ reduces to
         !! \f{equation}
         !!    M =  \begin{bmatrix}  
         !!             C_1 C_3 \pm S_1 S_3     &  0        &  \pm C_1 S_3 - S_1 C_3   \\
         !!             \pm S_1 C_3 - C_1 S_3   &  0        &  \pm C_1 C_3 + S_1 S_3   \\
         !!             0                       &  \mp 1    &  0                       \\
         !!         \end{bmatrix}
         !! \f}
         !!
         !! At this point we can choose \f$ \theta_3 = 0 \f$ due to gimbal lock giving \f$ \sin(\theta_3) = 0 \f$, \f$ \cos(\theta_3) = 1\f$.

         theta(3) = 0.0_SiKi                          ! theta_z = theta_3

         !> This further reduces \f$ M \f$ to 
         !! \f{equation}
         !!    M =  \begin{bmatrix}  
         !!             C_1      &  0        &  - S_1    \\
         !!             \pm S_1  &  0        &  \pm C_1  \\
         !!             0        &  \mp 1    &  0        \\
         !!         \end{bmatrix},
         !! \f}
         !! allowing us to solve for \f$ \theta_1 \f$ by \f$ \theta_1 = \arctan{\left( \frac{-M(1,3)}{M(1,1)} \right)} = \arctan{\left( \frac{S_1}{C_1} \right)}\f$.

         theta(1) = atan2( -m(1,3), m(1,1) )

      else
         !> ## Else \f$ \cos(\theta_2) = C_2 \neq 0\f$:
         !!
         !! First, start by finding \f$ \theta(1) \f$ from \f$ M(3,1) \f$ and \f$ M(3,3) \f$ using
         !! \f{equation}{
         !!    \theta_1 = \arctan{\left( \frac{M(3,1)}{M(3,3)} \right)} = \arctan{\left( \frac{S_1 C_2}{C_1 C_2} \right)}.
         !! \f}
         !! With this we calculate values for \f$S_1\f$ and \f$C_1\f$.

         theta(1) = atan2( m(3,1), m(3,3) )     ! theta_1 -> theta_y
         C_1   = cos( theta(1) )
         S_1   = sin( theta(1) )

         !> We already know \f$ \text{abs}( C_2 ) \f$, but need the sign of it.  This can be found by comparing the
         !! \f$ S_1 C_2 \f$ and \f$ C_1 C_2 \f$ terms with the \f$ C_1 \f$ and \f$ S_1 \f$ terms we just found.
         !! If \f$ C_1 = 0 \f$, then we use 
         !! \f{equation}{
         !!    C_2 = C_2 \cdot \text{sgn}{\left( \frac{M(3,1)}{S_1} \right)} = C_2 \cdot \text{sgn}{( C_2 )},
         !! \f}
         !! otherwise
         !! \f{equation}{
         !!    C_2 = C_2 \cdot \text{sgn}{\left( \frac{M(3,3)}{C_1} \right)} = C_2 \cdot \text{sgn}{( C_2 )}
         !! \f}
         !!
         if ( EqualRealNos( C_1, 0.0_SiKi ) ) then
            C_2 = sign( C_2, m(3,1) / S_1 )
         else
            C_2 = sign( C_2, m(3,3) / C_1 )
         endif

         !> Now can calculate \f$ \theta_2 \f$ from
         !! \f{equation}{
         !!    \theta_2 = \arctan{\left( \frac{-M(3,2)}{C_2} \right)} = \arctan{\left( \frac{S_2}{C_2} \right)}
         !! \f}
         theta(2) = atan2( -m(3,2), C_2 )
        

         !> For numerical reasons, we're going to get \f$ \theta_3 \f$ (\f$\theta_z\f$) using
         !! \f{eqnarray*}{
         !!    M' &=& M \cdot (R(\theta_2) \cdot R(\theta_1))^\text{T} = M \cdot R(\theta_1)^\text{T} \cdot R(\theta_2)^\text{T} & = & R(\theta_3) \\
         !!       &=& R(\theta_3) R(\theta_2) R(\theta_1) R(\theta_1)^T R(\theta_2)^T  &=& R(\theta_3) \\
         !!       &=&   M \cdot
         !!             \begin{bmatrix}
         !!                C_1   &  0     &  S_1   \\
         !!                0     &  1     &  0     \\
         !!                -S_1  &  0     &  C_1
         !!             \end{bmatrix}
         !!                \cdot
         !!             \begin{bmatrix}
         !!                1     &  0     &  0     \\
         !!                0     &  C_2   & -S_2   \\
         !!                0     &  S_2   &  C_2
         !!             \end{bmatrix}
         !!       &=&  
         !!             \begin{bmatrix}
         !!                C_3   &  S_3   &  0     \\
         !!                -S_3  &  C_3   &  0     \\
         !!                0     &  0     &  1
         !!             \end{bmatrix}  \\
         !!       &=&   M \cdot
         !!             \begin{bmatrix}
         !!                C_1   &  S_1 S_2  &  S_1 C_2  \\
         !!                0     &  C_2      &  -S_2     \\
         !!                -S_1  &  C_1 S_2  &  C_1 C_2
         !!             \end{bmatrix}
         !!       &=&   \begin{bmatrix}
         !!                C_3   &  S_3   &  0     \\
         !!                -S_3  &  C_3   &  0     \\
         !!                0     &  0     &  1
         !!             \end{bmatrix}  \\
         !! \f}
         !!
         !! From this we can find \f$ -S_3 \f$ and \f$ C_3 \f$ as
         !! \f{eqnarray}{
         !!    -S_3  &=&   M(2,1) C_1 + M(2,3) (- S_1 )  &=&   ( S_1 S_2 C_3 - C_1 S_3 ) C_1 + ( C_1 S_2 C_3 + S_1 S_3 ) ( - S_1 )  \\
         !!          &&                                  &=&   S_1 C_1 S_2 C_3 - C_1^2 S_3 - S_1^2 S_3 - S_1 C_1 S_2 C_3            \\
         !!          &&                                  &=&   -( C_1^2 + S_1^2 ) S_3                                               \\
         !!          &&                                  &=&   -S_3
         !! \f}
         !! and
         !! \f{eqnarray}{
         !!    C_3   &=&   M(1,1) C_1 + M(1,3) (- S_1 )  &=&   ( C_1 C_3 + S_1 S_2 S_3 ) C_1 + ( C_1 S_2 S_3 - S_1 C_3 ) (- S_1 )   \\
         !!          &&                                  &=&   C_1^2 C_3 + S_1 C_1 S_2 S_3 - S_1 C_1 S_2 S_3 + S_1^2 C_3            \\
         !!          &&                                  &=&   ( C_1^2 + S_1^2 ) C_3                                                \\
         !!          &&                                  &=&   C_3
         !! \f}
         !!
         !! \f$\theta_3\f$ is then found as \f$\theta_3 = \arctan{\left( \frac{S_3}{C_3} \right)}\f$.


         S_3 =    -( m(2,1) * C_1   + m(2,3) * (- S_1)   )
         C_3 =       m(1,1) * C_1   + m(1,3) * (- S_1)

         theta(3) = atan2( S_3, C_3)

      endif

      
   END FUNCTION TaitBryanYXZExtractR4

!> See nwtc_num::taitbryanyxzextractr4 for detailed explanation of algorithm
   FUNCTION TaitBryanYXZExtractR8(M) result(theta)
   
   
      REAL(R8Ki), INTENT(IN) :: M(3,3)    !< rotation matrix, M 
      REAL(R8Ki)             :: theta(3)  !< the 3 rotation angles, \f$(\theta_y, \theta_x, \theta_z)\f$, corresponding to the Tait-Bryan rotation angle corresponding to cant-toe-twist
     
      REAL(R8Ki)              :: C_1       ! C_1 > cos(theta_y) 
      REAL(R8Ki)              :: S_1       ! S_1 > sin(theta_y) 
      REAL(R8Ki)              :: C_2       ! C_2 > cos(theta_x) 
      REAL(R8Ki)              :: C_3       ! C_3 > cos(theta_z) 
      REAL(R8Ki)              :: S_3       ! S_3 > sin(theta_z)

         !> See nwtc_num::taitbryanyxzextractr4 for detailed description of how this works.

         ! use trig identity S_3**2 + C_3**2 = 1 to get abs( C_2 )
      C_2 = sqrt( m(1,2)**2 + m(2,2)**2 )

         ! If C_2 is zero, we can simplifiy some things since theta(2) is +/- pi/2
      if ( EqualRealNos( C_2, 0.0_R8Ki ) ) then

         ! find sign of theta(2) based on sin(theta_2)
         theta(2) = atan2( -m(3,2), C_2 )     ! theta_2 -> theta_x

         ! Considering C_2 = 0  and  S_2 = \pm 1, the matrix M reduces to
         !     M =  [   C_1 C_3 \pm S_1 S_3        0           \pm C_1 S_3 - S_1 C_3  |
         !          |   \pm S_1 C_3 - C_1 S_3      0           \pm C_1 C_3 + S_1 S_3  |
         !          |   0                          0  \mp 1    0                      ]
         ! 
         !  At this point we can choose \theta_3 = 0 due to gimbal lock giving sin(theta(3)) = 0, cos(theta(3)) = 1.

         theta(3) = 0.0_R8Ki                          ! theta_z = theta_3

         ! This further reduces M to 
         !     M =  [   C_1         0           - S_1    |
         !          |   \pm S_1     0           \pm C_1  |
         !          |   0           \mp 1       0        ]
         !
         !
         ! allowing us to solve for theta_1  by  theta_1 = atan2( -M(1,3), M(1,1) ) = atan2( S_1, C_1).

         theta(1) = atan2( -m(1,3), m(1,1) )

      else
         ! First, start by finding \f$ \theta(1) \f$ from \f$ M(3,1) \f$ and \f$ M(3,3) \f$ using
         !
         !    theta(1) = atan2( M(3,1), M(3,3) ) = atan2( S_1 * C_2, C_1 * C_2 ).
         ! With this we calculate values for S_1 and C_1.

         theta(1) = atan2( m(3,1), m(3,3) )     ! theta_1 -> theta_y
         C_1   = cos( theta(1) )
         S_1   = sin( theta(1) )

         !  We already know abs( C_2 ), but need the sign of it.  This can be found by comparing the
         !  S_1 * C_2 and C_1 * C_2 terms with the C_1 and S_1 terms we just found.
          
         if ( EqualRealNos( C_1, 0.0_R8Ki ) ) then
            C_2 = sign( C_2, m(3,1) / S_1 )
         else
            C_2 = sign( C_2, m(3,3) / C_1 )
         endif

         ! Now can calculate theta(2)
         theta(2) = atan2( -m(3,2), C_2 )
        

         ! For numerical reasons, we're going to get theta(3) using some matrix math and identities about M.
         !  See nwtc_num::taitbryanyxzextractr4 for complete documentation on the matrix math used here

         S_3 =    -( m(2,1) * C_1   + m(2,3) * (- S_1)   )
         C_3 =       m(1,1) * C_1   + m(1,3) * (- S_1)

         theta(3) = atan2( S_3, C_3)

      endif

      
   END FUNCTION TaitBryanYXZExtractR8


      FUNCTION TaitBryanYXZConstructR4(theta) result(M)
            ! this function creates a rotation matrix, M, from a 1-2-3 rotation
      ! sequence of the 3 TaitBryan angles, theta_x, theta_y, and theta_z, in radians.
      ! M represents a change of basis (from global to local coordinates; 
      ! not a physical rotation of the body). it is the inverse of TaitBryanYXZExtract().
      ! 
      ! M = R(theta_z) * R(theta_x) * R(theta_y)
      !   = [ cz sz 0 |  [ 1   0   0 |    [ cy  0 -sy | 
      !     |-sz cz 0 |* | 0  cx  sx |  * |  0  1   0 | 
      !     |  0  0 1 ]  | 0 -sx  cx ]    | sy  0  cy ] 
      !   = [ cy*cz+sy*sx*sz   cx*sz    cy*sx*sz-cz*sy |
      !     |cz*sy*sx-cy*sz   cx*cz    cy*cz*sx+sy*sz |
      !     |cx*sy           -sx             cx*cy    ]
      ! where cz = cos(theta_z), sz = sin(theta_z), cy = cos(theta_y), etc.
   
      REAL(SiKi)             :: M(3,3)    !< rotation matrix, M 
      REAL(SiKi), INTENT(IN) :: theta(3)  !< the 3 rotation angles: \f$\theta_x, \theta_y, \theta_z\f$
      
      REAL(SiKi)             :: cx        ! cos(theta_x)
      REAL(SiKi)             :: sx        ! sin(theta_x)
      REAL(SiKi)             :: cy        ! cos(theta_y)
      REAL(SiKi)             :: sy        ! sin(theta_y)
      REAL(SiKi)             :: cz        ! cos(theta_z)
      REAL(SiKi)             :: sz        ! sin(theta_z)
   
      cx = cos( theta(1) )
      sx = sin( theta(1) )
      
      cy = cos( theta(2) )
      sy = sin( theta(2) )
      
      cz = cos( theta(3) )
      sz = sin( theta(3) )
         
      M(1,1) =  cy*cz+sy*sx*sz            
      M(1,2) =  cx*sz            
      M(1,3) =  cy*sx*sz-cz*sy    
      
      M(2,1) =  cz*sy*sx-cy*sz            
      M(2,2) =  cx*cz            
      M(2,3) =  cy*cz*sx+sy*sz     
      
      M(3,1) =   cx*sy            
      M(3,2) =  -sx            
      M(3,3) =   cy*cx   

   END FUNCTION TaitBryanYXZConstructR4
   
   FUNCTION TaitBryanYXZConstructR8(theta) result(M)
   
      ! this function creates a rotation matrix, M, from a 1-2-3 rotation
      ! sequence of the 3 TaitBryan angles, theta_x, theta_y, and theta_z, in radians.
      ! M represents a change of basis (from global to local coordinates; 
      ! not a physical rotation of the body). it is the inverse of TaitBryanYXZExtract().
      ! 
      ! M = R(theta_z) * R(theta_x) * R(theta_y)
      !   = [ cz sz 0 |  [ 1   0   0 |    [ cy  0 -sy | 
      !     |-sz cz 0 |* | 0  cx  sx |  * |  0  1   0 | 
      !     |  0  0 1 ]  | 0 -sx  cx ]    | sy  0  cy ] 
      !   = [ cy*cz+sy*sx*sz   cx*sz    cy*sx*sz-cz*sy |
      !     |cz*sy*sx-cy*sz   cx*cz    cy*cz*sx+sy*sz |
      !     |cx*sy           -sx             cx*cy    ]
      ! where cz = cos(theta_z), sz = sin(theta_z), cy = cos(theta_y), etc.
   
      REAL(R8Ki)             :: M(3,3)    ! rotation matrix M 
      REAL(R8Ki), INTENT(IN) :: theta(3)  ! the 3 rotation angles: theta_x, theta_y, theta_z
      
      REAL(R8Ki)             :: cx        ! cos(theta_x)
      REAL(R8Ki)             :: sx        ! sin(theta_x)
      REAL(R8Ki)             :: cy        ! cos(theta_y)
      REAL(R8Ki)             :: sy        ! sin(theta_y)
      REAL(R8Ki)             :: cz        ! cos(theta_z)
      REAL(R8Ki)             :: sz        ! sin(theta_z)
   
      cx = cos( theta(1) )
      sx = sin( theta(1) )
      
      cy = cos( theta(2) )
      sy = sin( theta(2) )
      
      cz = cos( theta(3) )
      sz = sin( theta(3) )
         
      M(1,1) =  cy*cz+sy*sx*sz            
      M(1,2) =  cx*sz            
      M(1,3) =  cy*sx*sz-cz*sy    
      
      M(2,1) =  cz*sy*sx-cy*sz            
      M(2,2) =  cx*cz            
      M(2,3) =  cy*cz*sx+sy*sz     
      
      M(3,1) =   cx*sy            
      M(3,2) =  -sx            
      M(3,3) =   cy*cx               
   
   END FUNCTION TaitBryanYXZConstructR8
   

!=======================================================================
!> This routine takes an array of time values such as that returned from
!!     CALL DATE_AND_TIME ( Values=TimeAry )
!! and converts TimeAry to the number of seconds past midnight.
   FUNCTION TimeValues2Seconds( TimeAry )

      ! Passed variables:
   INTEGER, INTENT(IN)          :: TimeAry  (8)                                    ! An array containing the elements of the time
   REAL(ReKi)                   :: TimeValues2Seconds                              ! Current time in seconds past midnight


   TimeValues2Seconds = 3600*TimeAry(5) + 60*TimeAry(6) + TimeAry(7) + 0.001_ReKi*TimeAry(8)

   END FUNCTION TimeValues2Seconds
!=======================================================================
!> This function computes the trace of a matrix \f$A \in \mathbb{R}^{m,n}\f$. The 
!! trace of \f$A\f$, \f$\mathrm{Tr}\left[ A \right]\f$, is the sum of the diagonal elements of \f$A\f$:   
!! \f{equation}{   
!!   \mathrm{Tr}\left[ A \right] = \sum_{i=1}^{\min(m,n)} A(i,i)
!! \f}   
!!
!! Use trace (nwtc_num::trace) instead of directly calling a specific routine in the generic interface.    
   FUNCTION traceR4(A)
         
   REAL(SiKi), INTENT(IN)  :: A(:,:)     !< matrix A
   REAL(SiKi)              :: traceR4    !< sum of the diagonal elements of A
   
   INTEGER(IntKi)          :: n     ! rows/cols in A
   INTEGER(IntKi)          :: i     ! loop counter
   
   n = min( SIZE(A,1), SIZE(A,2) )

   traceR4 = 0.0_ReKi
   do i=1,n
      traceR4 = traceR4 + A(i,i)
   end do
   
   END FUNCTION traceR4
!=======================================================================
!> \copydoc nwtc_num::tracer4
   FUNCTION traceR8(A)
         
   REAL(R8Ki), INTENT(IN)  :: A(:,:)  !< matrix A
   REAL(R8Ki)              :: traceR8 !< sum of the diagonal elements of A
   
   INTEGER(IntKi)          :: n     ! rows/cols in A
   INTEGER(IntKi)          :: i     ! loop counter
   
   n = min( SIZE(A,1), SIZE(A,2) )

   traceR8 = 0.0_ReKi
   do i=1,n
      traceR8 = traceR8 + A(i,i)
   end do
   
   END FUNCTION traceR8

!=======================================================================
!> This function returns the \f$l_2\f$ (Euclidian) norm of a vector, 
!! \f$v = \left(v_1, v_2, \ldots ,v_n\right)\f$. The \f$l_2\f$-norm is defined as   
!! \f{equation}{   
!!  \lVert v \rVert_2 = \left( \sum_{i=1}^{n} {v_i}^2 \right)^{1/2}
!! \f} \n
!! Use TwoNorm (nwtc_num::twonorm) instead of directly calling a specific routine in the generic interface.    
   FUNCTION TwoNormR4(v)
   
      ! fortran 2008 has Norm2() built in
      
      REAL(SiKi), INTENT(IN)  :: v(:)           !< vector, v
      REAL(SiKi)              :: TwoNormR4      !< two-norm of v
      
      TwoNormR4 = SQRT( DOT_PRODUCT(v, v) )
      
      
   END FUNCTION
!=======================================================================
!> \copydoc nwtc_num::twonormr4
   FUNCTION TwoNormR8(v)
   
      ! this function returns the 2-norm of a vector v
      ! fortran 2008 has Norm2() built in
      
      REAL(R8Ki), INTENT(IN)  :: v(:)      
      REAL(R8Ki)              :: TwoNormR8      
      
      TwoNormR8 = SQRT( DOT_PRODUCT(v, v) )
      
      
   END FUNCTION

!=======================================================================  
!> This routine is used to convert Angle to an equivalent value
!!  in the range \f$[0, 2\pi)\f$. \n
!! Use Zero2TwoPi (nwtc_num::zero2twopi) instead of directly calling a specific routine in the generic interface.    
   SUBROUTINE Zero2TwoPiR4 ( Angle )
      
      ! Argument declarations:

   REAL(SiKi), INTENT(INOUT)    :: Angle     !< angle that is input and converted to equivalent in range \f$[0, 2\pi)\f$


      ! Get the angle between 0 and 2Pi.

   Angle = MODULO( Angle, TwoPi_R4 )


      ! Check numerical case where Angle == 2Pi.

   IF ( Angle == TwoPi_R4 )  THEN
      Angle = 0.0_ReKi
   END IF


   RETURN
   END SUBROUTINE Zero2TwoPiR4
!=======================================================================  
!> \copydoc nwtc_num::zero2twopir4
   SUBROUTINE Zero2TwoPiR8 ( Angle )

      ! This routine is used to convert Angle to an equivalent value
      !  in the range [0, 2*pi).
      

      ! Argument declarations:

   REAL(R8Ki), INTENT(INOUT)    :: Angle



      ! Get the angle between 0 and 2Pi.

   Angle = MODULO( Angle, TwoPi_R8 )


      ! Check numerical case where Angle == 2Pi.

   IF ( Angle == TwoPi_R8 )  THEN
      Angle = 0.0_DbKi
   END IF


   RETURN
   END SUBROUTINE Zero2TwoPiR8   
!=======================================================================
   !< This routine extrapolates or interpolates between angles
   SUBROUTINE Angles_ExtrapInterp1_R4(Angle1, Angle2, tin, Angle_out, tin_out )
       REAL(SiKi),          INTENT(IN   )  :: Angle1 !< Angle at t1 > t2
       REAL(SiKi),          INTENT(IN   )  :: Angle2 !< Angle at t2
       REAL(R8Ki),          INTENT(IN   )  :: tin(:)                    !< Times associated with the inputs
       REAL(SiKi),          INTENT(INOUT)  :: Angle_out                 !< Input at tin_out
       REAL(R8Ki),          INTENT(IN   )  :: tin_out                   !< time to be extrap/interp'd to
                                                                     
         ! local variables                                              
       INTEGER(IntKi), parameter           :: order = 1                 ! order of polynomial fit (max 2)
       REAL(R8Ki)                          :: t(SIZE(tin))              ! Times associated with the inputs
       REAL(R8Ki)                          :: t_out                     ! Time to which to be extrap/interpd
                                                                     
       REAL(SiKi)                          :: Angle2_mod
    
          ! we'll subtract a constant from the times to resolve some
          ! numerical issues when t gets large (and to simplify the equations)
       t = tin - tin(1)
       t_out = tin_out - tin(1)

      !    ! some error checking:
      !
      ! if ( size(t) .ne. order+1) then
      !    ErrStat = ErrID_Fatal
      !    ErrMsg = 'Angles_ExtrapInterp1: size(t) must equal 2.'
      !    RETURN
      ! end if
      !
      !IF ( EqualRealNos( t(1), t(2) ) ) THEN
      !   ErrStat = ErrID_Fatal
      !   ErrMsg  = 'Angles_ExtrapInterp1: t(1) must not equal t(2) to avoid a division-by-zero error.'
      !   RETURN
      !END IF

      Angle2_mod = Angle2
      call AddOrSub2Pi( Angle1, Angle2_mod )
      
      Angle_out = Angle1 + (Angle2_mod - Angle1) * t_out / t(2)
      
!     call Zero2TwoPi(Angle_out)
!      call MPi2Pi(Angle_out)

   END SUBROUTINE Angles_ExtrapInterp1_R4
!=======================================================================  
   !< This routine extrapolates or interpolates between angles
   SUBROUTINE Angles_ExtrapInterp1_R8(Angle1, Angle2, tin, Angle_out, tin_out)
       REAL(R8Ki),          INTENT(IN   )  :: Angle1 !< Angle at t1 > t2
       REAL(R8Ki),          INTENT(IN   )  :: Angle2 !< Angle at t2
       REAL(R8Ki),          INTENT(IN   )  :: tin(:)                    !< Times associated with the inputs
       REAL(R8Ki),          INTENT(INOUT)  :: Angle_out                 !< Input at tin_out
       REAL(R8Ki),          INTENT(IN   )  :: tin_out                   !< time to be extrap/interp'd to
         
         ! local variables                                              
       INTEGER(IntKi), parameter           :: order = 1                 ! order of polynomial fit (max 2)
       REAL(R8Ki)                          :: t(SIZE(tin))              ! Times associated with the inputs
       REAL(R8Ki)                          :: t_out                     ! Time to which to be extrap/interpd
                                                                     
       REAL(R8Ki)                          :: Angle2_mod
    
          ! we'll subtract a constant from the times to resolve some
          ! numerical issues when t gets large (and to simplify the equations)
       t = tin - tin(1)
       t_out = tin_out - tin(1)

      !    ! some error checking:
      !
      ! if ( size(t) .ne. order+1) then
      !    ErrStat = ErrID_Fatal
      !    ErrMsg = 'Angles_ExtrapInterp1: size(t) must equal 2.'
      !    RETURN
      ! end if
      !
      !IF ( EqualRealNos( t(1), t(2) ) ) THEN
      !   ErrStat = ErrID_Fatal
      !   ErrMsg  = 'Angles_ExtrapInterp1: t(1) must not equal t(2) to avoid a division-by-zero error.'
      !   RETURN
      !END IF

      Angle2_mod = Angle2
      call AddOrSub2Pi( Angle1, Angle2_mod )
      
      Angle_out = Angle1 + (Angle2_mod - Angle1) * t_out / t(2)
!     call Zero2TwoPi(Angle_out)
!      call MPi2Pi(Angle_out)

   END SUBROUTINE Angles_ExtrapInterp1_R8
!=======================================================================
   !< This routine extrapolates or interpolates between angles
   SUBROUTINE Angles_ExtrapInterp1_R4R(Angle1, Angle2, tin, Angle_out, tin_out )
       REAL(SiKi),          INTENT(IN   )  :: Angle1 !< Angle at t1 > t2
       REAL(SiKi),          INTENT(IN   )  :: Angle2 !< Angle at t2
       REAL(R4Ki),          INTENT(IN   )  :: tin(:)                    !< Times associated with the inputs
       REAL(SiKi),          INTENT(INOUT)  :: Angle_out                 !< Input at tin_out
       REAL(R4Ki),          INTENT(IN   )  :: tin_out                   !< time to be extrap/interp'd to
                                                                     
         ! local variables                                              
       INTEGER(IntKi), parameter           :: order = 1                 ! order of polynomial fit (max 2)
       REAL(SiKi)                          :: t(SIZE(tin))              ! Times associated with the inputs
       REAL(SiKi)                          :: t_out                     ! Time to which to be extrap/interpd
                                                                     
       REAL(SiKi)                          :: Angle2_mod
    
          ! we'll subtract a constant from the times to resolve some
          ! numerical issues when t gets large (and to simplify the equations)
       t = tin - tin(1)
       t_out = tin_out - tin(1)

      !    ! some error checking:
      !
      ! if ( size(t) .ne. order+1) then
      !    ErrStat = ErrID_Fatal
      !    ErrMsg = 'Angles_ExtrapInterp1: size(t) must equal 2.'
      !    RETURN
      ! end if
      !
      !IF ( EqualRealNos( t(1), t(2) ) ) THEN
      !   ErrStat = ErrID_Fatal
      !   ErrMsg  = 'Angles_ExtrapInterp1: t(1) must not equal t(2) to avoid a division-by-zero error.'
      !   RETURN
      !END IF

      Angle2_mod = Angle2
      call AddOrSub2Pi( Angle1, Angle2_mod )
      
      Angle_out = Angle1 + (Angle2_mod - Angle1) * t_out / t(2)
      
!     call Zero2TwoPi(Angle_out)
!      call MPi2Pi(Angle_out)

   END SUBROUTINE Angles_ExtrapInterp1_R4R
!=======================================================================  
   !< This routine extrapolates or interpolates between angles
   SUBROUTINE Angles_ExtrapInterp1_R8R(Angle1, Angle2, tin, Angle_out, tin_out)
       REAL(R8Ki),          INTENT(IN   )  :: Angle1 !< Angle at t1 > t2
       REAL(R8Ki),          INTENT(IN   )  :: Angle2 !< Angle at t2
       REAL(SiKi),          INTENT(IN   )  :: tin(:)                    !< Times associated with the inputs
       REAL(R8Ki),          INTENT(INOUT)  :: Angle_out                 !< Input at tin_out
       REAL(SiKi),          INTENT(IN   )  :: tin_out                   !< time to be extrap/interp'd to
         
         ! local variables                                              
       INTEGER(IntKi), parameter           :: order = 1                 ! order of polynomial fit (max 2)
       REAL(SiKi)                          :: t(SIZE(tin))              ! Times associated with the inputs
       REAL(SiKi)                          :: t_out                     ! Time to which to be extrap/interpd
                                                                     
       REAL(R8Ki)                          :: Angle2_mod
    
          ! we'll subtract a constant from the times to resolve some
          ! numerical issues when t gets large (and to simplify the equations)
       t = tin - tin(1)
       t_out = tin_out - tin(1)

      !    ! some error checking:
      !
      ! if ( size(t) .ne. order+1) then
      !    ErrStat = ErrID_Fatal
      !    ErrMsg = 'Angles_ExtrapInterp1: size(t) must equal 2.'
      !    RETURN
      ! end if
      !
      !IF ( EqualRealNos( t(1), t(2) ) ) THEN
      !   ErrStat = ErrID_Fatal
      !   ErrMsg  = 'Angles_ExtrapInterp1: t(1) must not equal t(2) to avoid a division-by-zero error.'
      !   RETURN
      !END IF

      Angle2_mod = Angle2
      call AddOrSub2Pi( Angle1, Angle2_mod )
      
      Angle_out = Angle1 + (Angle2_mod - Angle1) * t_out / t(2)
!     call Zero2TwoPi(Angle_out)
!      call MPi2Pi(Angle_out)

   END SUBROUTINE Angles_ExtrapInterp1_R8R
!=======================================================================  
   !< This routine extrapolates or interpolates between angles
   SUBROUTINE Angles_ExtrapInterp2_R4(Angle1, Angle2, Angle3, tin, Angle_out, tin_out )
       REAL(SiKi),          INTENT(IN   )  :: Angle1 !< Angle at t1 > t2 > t3
       REAL(SiKi),          INTENT(IN   )  :: Angle2 !< Angle at t2 > t3
       REAL(SiKi),          INTENT(IN   )  :: Angle3 !< Angle at t3
       REAL(DbKi),          INTENT(IN   )  :: tin(:)                    !< Times associated with the inputs
       REAL(SiKi),          INTENT(INOUT)  :: Angle_out                 !< Input at tin_out
       REAL(DbKi),          INTENT(IN   )  :: tin_out                   !< time to be extrap/interp'd to
                                                                     
         ! local variables                                              
       INTEGER(IntKi), parameter           :: order = 2                 ! order of polynomial fit (max 2)
       REAL(DbKi)                          :: t(SIZE(tin))              ! Times associated with the inputs
       REAL(DbKi)                          :: t_out                     ! Time to which to be extrap/interpd
                                                                     
       REAL(DbKi)                          :: scaleFactor               ! temporary for extrapolation/interpolation
       REAL(SiKi)                          :: Angle2_mod
       REAL(SiKi)                          :: Angle3_mod
    
          ! we'll subtract a constant from the times to resolve some
          ! numerical issues when t gets large (and to simplify the equations)
       t = tin - tin(1)
       t_out = tin_out - tin(1)

      !    ! some error checking:
      !
      !if ( size(t) .ne. order+1) then
      !   ErrStat = ErrID_Fatal
      !   ErrMsg = 'Angles_ExtrapInterp2: size(t) must equal 3.'
      !   RETURN
      !end if
      !
      !IF ( EqualRealNos( t(1), t(2) ) ) THEN
      !   ErrStat = ErrID_Fatal
      !   ErrMsg  = 'Angles_ExtrapInterp2: t(1) must not equal t(2) to avoid a division-by-zero error.'
      !   RETURN
      !END IF
      !IF ( EqualRealNos( t(2), t(3) ) ) THEN
      !   ErrStat = ErrID_Fatal
      !   ErrMsg  = 'Angles_ExtrapInterp2: t(2) must not equal t(3) to avoid a division-by-zero error.'
      !   RETURN
      !END IF
      !IF ( EqualRealNos( t(1), t(3) ) ) THEN
      !   ErrStat = ErrID_Fatal
      !   ErrMsg  = 'Angles_ExtrapInterp2: t(1) must not equal t(3) to avoid a division-by-zero error.'
      !   RETURN
      !END IF

      Angle2_mod = Angle2
      Angle3_mod = Angle3
      call AddOrSub2Pi( Angle1, Angle2_mod )
      call AddOrSub2Pi( Angle2_mod, Angle3_mod )
      
      scaleFactor = t_out / ( t(2) * t(3) * (t(2) - t(3)) )

      Angle_out =   Angle1 &
                     + ( t(3)**2 * (Angle1 - Angle2_mod) + t(2)**2*(-Angle1 + Angle3_mod) ) * scaleFactor &
                     + ( (t(2)-t(3))*Angle1 + t(3)*Angle2_mod - t(2)*Angle3_mod ) *scaleFactor * t_out
                     
!     call Zero2TwoPi(Angle_out)
!      call MPi2Pi(Angle_out)
      
   END SUBROUTINE Angles_ExtrapInterp2_R4
!=======================================================================  
   !< This routine extrapolates or interpolates between angles
   SUBROUTINE Angles_ExtrapInterp2_R8(Angle1, Angle2, Angle3, tin, Angle_out, tin_out)
       REAL(R8Ki),          INTENT(IN   )  :: Angle1 !< Angle at t1 > t2 > t3
       REAL(R8Ki),          INTENT(IN   )  :: Angle2 !< Angle at t2 > t3
       REAL(R8Ki),          INTENT(IN   )  :: Angle3 !< Angle at t3
       REAL(DbKi),          INTENT(IN   )  :: tin(:)                    !< Times associated with the inputs
       REAL(R8Ki),          INTENT(INOUT)  :: Angle_out                 !< Input at tin_out
       REAL(DbKi),          INTENT(IN   )  :: tin_out                   !< time to be extrap/interp'd to
                                                                     
         ! local variables                                              
       INTEGER(IntKi), parameter           :: order = 2                 ! order of polynomial fit (max 2)
       REAL(DbKi)                          :: t(SIZE(tin))              ! Times associated with the inputs
       REAL(DbKi)                          :: t_out                     ! Time to which to be extrap/interpd
                                                                     
       REAL(DbKi)                          :: scaleFactor               ! temporary for extrapolation/interpolation
       REAL(R8Ki)                          :: Angle2_mod
       REAL(R8Ki)                          :: Angle3_mod
    
          ! we'll subtract a constant from the times to resolve some
          ! numerical issues when t gets large (and to simplify the equations)
       t = tin - tin(1)
       t_out = tin_out - tin(1)

          ! some error checking:

      !if ( size(t) .ne. order+1) then
      !   ErrStat = ErrID_Fatal
      !   ErrMsg = 'Angles_ExtrapInterp2: size(t) must equal 3.'
      !   RETURN
      !end if
      !
      !IF ( EqualRealNos( t(1), t(2) ) ) THEN
      !   ErrStat = ErrID_Fatal
      !   ErrMsg  = 'Angles_ExtrapInterp2: t(1) must not equal t(2) to avoid a division-by-zero error.'
      !   RETURN
      !END IF
      !IF ( EqualRealNos( t(2), t(3) ) ) THEN
      !   ErrStat = ErrID_Fatal
      !   ErrMsg  = 'Angles_ExtrapInterp2: t(2) must not equal t(3) to avoid a division-by-zero error.'
      !   RETURN
      !END IF
      !IF ( EqualRealNos( t(1), t(3) ) ) THEN
      !   ErrStat = ErrID_Fatal
      !   ErrMsg  = 'Angles_ExtrapInterp2: t(1) must not equal t(3) to avoid a division-by-zero error.'
      !   RETURN
      !END IF

      Angle2_mod = Angle2
      Angle3_mod = Angle3
      call AddOrSub2Pi( Angle1, Angle2_mod )
      call AddOrSub2Pi( Angle2_mod, Angle3_mod )
      
      scaleFactor = t_out / ( t(2) * t(3) * (t(2) - t(3)) )

      Angle_out =   Angle1 &
                     + ( t(3)**2 * (Angle1 - Angle2_mod) + t(2)**2*(-Angle1 + Angle3_mod) ) * scaleFactor &
                     + ( (t(2)-t(3))*Angle1 + t(3)*Angle2_mod - t(2)*Angle3_mod ) *scaleFactor * t_out
!     call Zero2TwoPi(Angle_out)
!      call MPi2Pi(Angle_out)
      
   END SUBROUTINE Angles_ExtrapInterp2_R8

!=======================================================================  
   !< This routine extrapolates or interpolates between angles
   SUBROUTINE Angles_ExtrapInterp2_R4R(Angle1, Angle2, Angle3, tin, Angle_out, tin_out )
       REAL(SiKi),          INTENT(IN   )  :: Angle1 !< Angle at t1 > t2 > t3
       REAL(SiKi),          INTENT(IN   )  :: Angle2 !< Angle at t2 > t3
       REAL(SiKi),          INTENT(IN   )  :: Angle3 !< Angle at t3
       REAL(R4Ki),          INTENT(IN   )  :: tin(:)                    !< Times associated with the inputs
       REAL(SiKi),          INTENT(INOUT)  :: Angle_out                 !< Input at tin_out
       REAL(R4Ki),          INTENT(IN   )  :: tin_out                   !< time to be extrap/interp'd to
                                                                     
         ! local variables                                              
       INTEGER(IntKi), parameter           :: order = 2                 ! order of polynomial fit (max 2)
       REAL(R4Ki)                          :: t(SIZE(tin))              ! Times associated with the inputs
       REAL(R4Ki)                          :: t_out                     ! Time to which to be extrap/interpd
                                                                     
       REAL(R8Ki)                          :: scaleFactor               ! temporary for extrapolation/interpolation
       REAL(SiKi)                          :: Angle2_mod
       REAL(SiKi)                          :: Angle3_mod
    
          ! we'll subtract a constant from the times to resolve some
          ! numerical issues when t gets large (and to simplify the equations)
       t = tin - tin(1)
       t_out = tin_out - tin(1)

      !    ! some error checking:
      !
      !if ( size(t) .ne. order+1) then
      !   ErrStat = ErrID_Fatal
      !   ErrMsg = 'Angles_ExtrapInterp2: size(t) must equal 3.'
      !   RETURN
      !end if
      !
      !IF ( EqualRealNos( t(1), t(2) ) ) THEN
      !   ErrStat = ErrID_Fatal
      !   ErrMsg  = 'Angles_ExtrapInterp2: t(1) must not equal t(2) to avoid a division-by-zero error.'
      !   RETURN
      !END IF
      !IF ( EqualRealNos( t(2), t(3) ) ) THEN
      !   ErrStat = ErrID_Fatal
      !   ErrMsg  = 'Angles_ExtrapInterp2: t(2) must not equal t(3) to avoid a division-by-zero error.'
      !   RETURN
      !END IF
      !IF ( EqualRealNos( t(1), t(3) ) ) THEN
      !   ErrStat = ErrID_Fatal
      !   ErrMsg  = 'Angles_ExtrapInterp2: t(1) must not equal t(3) to avoid a division-by-zero error.'
      !   RETURN
      !END IF

      Angle2_mod = Angle2
      Angle3_mod = Angle3
      call AddOrSub2Pi( Angle1, Angle2_mod )
      call AddOrSub2Pi( Angle2_mod, Angle3_mod )
      
      scaleFactor = t_out / ( t(2) * t(3) * (t(2) - t(3)) )

      Angle_out =   Angle1 &
                     + ( t(3)**2 * (Angle1 - Angle2_mod) + t(2)**2*(-Angle1 + Angle3_mod) ) * scaleFactor &
                     + ( (t(2)-t(3))*Angle1 + t(3)*Angle2_mod - t(2)*Angle3_mod ) *scaleFactor * t_out
                     
!     call Zero2TwoPi(Angle_out)
!      call MPi2Pi(Angle_out)
      
   END SUBROUTINE Angles_ExtrapInterp2_R4R
!=======================================================================  
   !< This routine extrapolates or interpolates between angles
   SUBROUTINE Angles_ExtrapInterp2_R8R(Angle1, Angle2, Angle3, tin, Angle_out, tin_out)
       REAL(R8Ki),          INTENT(IN   )  :: Angle1 !< Angle at t1 > t2 > t3
       REAL(R8Ki),          INTENT(IN   )  :: Angle2 !< Angle at t2 > t3
       REAL(R8Ki),          INTENT(IN   )  :: Angle3 !< Angle at t3
       REAL(R4Ki),          INTENT(IN   )  :: tin(:)                    !< Times associated with the inputs
       REAL(R8Ki),          INTENT(INOUT)  :: Angle_out                 !< Input at tin_out
       REAL(R4Ki),          INTENT(IN   )  :: tin_out                   !< time to be extrap/interp'd to
                                                                     
         ! local variables                                              
       INTEGER(IntKi), parameter           :: order = 2                 ! order of polynomial fit (max 2)
       REAL(R4Ki)                          :: t(SIZE(tin))              ! Times associated with the inputs
       REAL(R4Ki)                          :: t_out                     ! Time to which to be extrap/interpd
                                                                     
       REAL(R8Ki)                          :: scaleFactor               ! temporary for extrapolation/interpolation
       REAL(R8Ki)                          :: Angle2_mod
       REAL(R8Ki)                          :: Angle3_mod
    
          ! we'll subtract a constant from the times to resolve some
          ! numerical issues when t gets large (and to simplify the equations)
       t = tin - tin(1)
       t_out = tin_out - tin(1)

          ! some error checking:

      !if ( size(t) .ne. order+1) then
      !   ErrStat = ErrID_Fatal
      !   ErrMsg = 'Angles_ExtrapInterp2: size(t) must equal 3.'
      !   RETURN
      !end if
      !
      !IF ( EqualRealNos( t(1), t(2) ) ) THEN
      !   ErrStat = ErrID_Fatal
      !   ErrMsg  = 'Angles_ExtrapInterp2: t(1) must not equal t(2) to avoid a division-by-zero error.'
      !   RETURN
      !END IF
      !IF ( EqualRealNos( t(2), t(3) ) ) THEN
      !   ErrStat = ErrID_Fatal
      !   ErrMsg  = 'Angles_ExtrapInterp2: t(2) must not equal t(3) to avoid a division-by-zero error.'
      !   RETURN
      !END IF
      !IF ( EqualRealNos( t(1), t(3) ) ) THEN
      !   ErrStat = ErrID_Fatal
      !   ErrMsg  = 'Angles_ExtrapInterp2: t(1) must not equal t(3) to avoid a division-by-zero error.'
      !   RETURN
      !END IF

      Angle2_mod = Angle2
      Angle3_mod = Angle3
      call AddOrSub2Pi( Angle1, Angle2_mod )
      call AddOrSub2Pi( Angle2_mod, Angle3_mod )
      
      scaleFactor = t_out / ( t(2) * t(3) * (t(2) - t(3)) )

      Angle_out =   Angle1 &
                     + ( t(3)**2 * (Angle1 - Angle2_mod) + t(2)**2*(-Angle1 + Angle3_mod) ) * scaleFactor &
                     + ( (t(2)-t(3))*Angle1 + t(3)*Angle2_mod - t(2)*Angle3_mod ) *scaleFactor * t_out
!     call Zero2TwoPi(Angle_out)
!      call MPi2Pi(Angle_out)
      
   END SUBROUTINE Angles_ExtrapInterp2_R8R
!=======================================================================  
END MODULE NWTC_Num
