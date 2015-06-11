!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2013-2015  National Renewable Energy Laboratory
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
!
!**********************************************************************************************************************************
! File last committed: $Date$
! (File) Revision #: $Rev$
! URL: $HeadURL$
!**********************************************************************************************************************************
MODULE NWTC_Num


   ! This module contains numeric-type routines with non-system-specific logic and references.


   ! It contains the following routines:

   !  SUBROUTINE AddOrSub2Pi           ( OldAngle, NewAngle )
   !  SUBROUTINE BSortReal             ( RealAry, NumPts )
   !  FUNCTION   CROSS_PRODUCT         ( Vector1, Vector2 )
   !  SUBROUTINE CubicSplineInit       ( AryLen, XAry, YAry, Coef )                                   ! Calculate coefficients for irregularly spaced array to use cubic splines.
   !  SUBROUTINE CubicSplineInitM      ( XAry, YAry, Coef, ErrStat, ErrMsg )                          ! Initialize cubic splines for multiple tables of irregularly space data.
   !  FUNCTION   CubicSplineInterp     ( X, AryLen, XAry, YAry,       Coef )                          ! Interpolate a irregularly spaced array using cubic splines.
   !  FUNCTION   CubicSplineInterpM    ( X, XAry, YAry, Coef, ErrStat, ErrMsg )                       ! Interpolate using cubic splines for multiple tables of irregularly space data.
   !  FUNCTION   EqualRealNos          ( ReNum1, ReNum2 )
   !  SUBROUTINE Eye                   ( A, ErrStat, ErrMsg )                                         ! sets A equal to the identity matrix (A can have 2 or 3 dimensions)
   !  FUNCTION   DCM_exp               ( lambda )                         
   !  SUBROUTINE DCM_logMap            ( DCM, logMap, ErrStat, ErrMsg )
   !  SUBROUTINE DCM_SetLogMapForInterp( tensor )
   !  SUBROUTINE GaussElim             ( AugMat, NumEq, x, ErrStat, ErrMsg )                          ! Performs Gauss-Jordan elimination to solve Ax=b for x; AugMat = [A b]
   !  SUBROUTINE GetOffsetReg          ( Ary, NumPts, Val, Ind, Fract, ErrStat, ErrMsg )              ! Determine index of the point in Ary just below Val and the fractional distance to the next point in the array.
   !  FUNCTION   GetSmllRotAngs        ( DCMat, ErrStat, ErrMsg )
   !  SUBROUTINE GL_Pts                ( IPt, NPts, Loc, Wt [, ErrStat] )
   !  FUNCTION   IndexCharAry          ( CVal, CAry )
   !  FUNCTION   InterpBin             ( XVal, XAry, YAry, ILo, AryLen )                              ! Generic interface for InterpBinComp and InterpBinReal.
   !     FUNCTION   InterpBinComp      ( XVal, XAry, YAry, ILo, AryLen )
   !     FUNCTION   InterpBinReal      ( XVal, XAry, YAry, ILo, AryLen )
   !  FUNCTION   InterpStp             ( XVal, XAry, YAry, ILo, AryLen )                              ! Generic interface for InterpStpComp and InterpStpReal.
   !     FUNCTION   InterpStpComp      ( XVal, XAry, YAry, Ind, AryLen )
   !     FUNCTION   InterpStpReal      ( XVal, XAry, YAry, Ind, AryLen )
   !  SUBROUTINE InterpStpReal2D       ( InCoord, Dataset, x, y, z, LastIndex, InterpData )
   !  SUBROUTINE InterpStpReal3D       ( InCoord, Dataset, x, y,    LastIndex, InterpData )   
   !  FUNCTION   InterpWrappedStpReal  ( XValIn, XAry, YAry, Ind, AryLen )
   !  SUBROUTINE IsoparametricCoords   ( InCoord, posLo, posHi, isopc )
   !  FUNCTION   IsSymmetric           ( A )                                                          ! Function to determine if A(:,:) is symmetric
   !  SUBROUTINE LocateBin             ( XVal, XAry, Ind, AryLen )
   !  SUBROUTINE LocateStp             ( XVal, XAry, Ind, AryLen )
   !  FUNCTION   Mean                  ( Ary, AryLen )                                                ! Function to calculate the mean value of a vector array.
   !  SUBROUTINE MPi2Pi                ( Angle )
   !  FUNCTION   PSF                   ( N, NumPrimes )                                               ! This routine factors the number N into its primes.  
   !  FUNCTION   Quaternion_Conjugate( q )
   !  FUNCTION   Quaternion_Norm( q )
   !  FUNCTION   Quaternion_Power( q, alpha )
   !  FUNCTION   Quaternion_Product( p, q )
   !  FUNCTION   Quaternion_to_DCM( q )
   !  FUNCTION   DCM_to_Quaternion( DCM )
   !  FUNCTION   Quaternion_Interp( q1, q2, s )
   !  SUBROUTINE RegCubicSplineInit    ( AryLen, XAry, YAry, DelX, Coef )                             ! Calculate coefficients for regularly spaced array to use cubic splines.
   !  SUBROUTINE RegCubicSplineInitM   ( XAry, YAry, DelX, Coef, ErrStat, ErrMsg )                    ! Interpolate using cubic splines for multiple tables of regularly space data.
   !  FUNCTION   RegCubicSplineInterp  ( X, AryLen, XAry, YAry, DelX, Coef )                          ! Interpolate a regularly spaced array using cubic splines.
   !  FUNCTION   RegCubicSplineInterpM ( X, XAry, YAry, DelX, Coef, ErrStat, ErrMsg )                 ! Initialize cubic splines for multiple tables of regularly space data.
   !  SUBROUTINE RombergInt            ( f, a, b, R, err, eps, ErrStat )
   !  SUBROUTINE SetAnglesForInterp    ( angles )                                                     ! uses 2pi periodicity of angles to set angles for interpolation (makes sure no two adjacent entries are more than pi apart)
   !  SUBROUTINE SetConstants
   !  SUBROUTINE SmllRotTrans          ( RotationType, Theta1, Theta2, Theta3, TransMat, ErrTxt )
   !  SUBROUTINE SortUnion             ( Ary1, N1, Ary2, N2, Ary, N )
   !  FUNCTION   StdDevFn              ( Ary, AryLen, Mean )                                          ! Function to calculate the standard deviation of a vector array.
   !  FUNCTION   trace                 ( A )                                                          ! computes the trace (sum of diagonal elements) of a matrix  (2-dimension array)
   !  FUNCTION   TwoNorm               ( v )                                                          ! computes the l2 norm of a vector (1-dimension array) 
   !  SUBROUTINE Zero2TwoPi            ( Angle )
   
   USE                                          NWTC_IO

   IMPLICIT NONE

!=======================================================================


      ! Global numeric-related variables.

   REAL(DbKi)                                :: D2R_D                         ! Factor to convert degrees to radians in double precision
   REAL(DbKi)                                :: Inf_D                         ! IEEE value for NaN (not-a-number) in double precision
   REAL(DbKi)                                :: Inv2Pi_D                      ! 0.5/Pi (1/(2*Pi)) in double precision
   REAL(DbKi)                                :: NaN_D                         ! IEEE value for Inf (infinity) in double precision
   REAL(DbKi)                                :: Pi_D                          ! Ratio of a circle's circumference to its diameter in double precision
   REAL(DbKi)                                :: PiBy2_D                       ! Pi/2 in double precision
   REAL(DbKi)                                :: R2D_D                         ! Factor to convert radians to degrees in double precision
   REAL(DbKi)                                :: RPM2RPS_D                     ! Factor to convert revolutions per minute to radians per second in double precision
   REAL(DbKi)                                :: RPS2RPM_D                     ! Factor to convert radians per second to revolutions per minute in double precision
   REAL(DbKi)                                :: TwoByPi_D                     ! 2/Pi in double precision
   REAL(DbKi)                                :: TwoPi_D                       ! 2*Pi in double precision


   REAL(ReKi)                                :: D2R                           ! Factor to convert degrees to radians
   REAL(ReKi)                                :: Inf                           ! IEEE value for NaN (not-a-number)
   REAL(ReKi)                                :: Inv2Pi                        ! 0.5/Pi = 1 / (2*pi)
   REAL(ReKi)                                :: NaN                           ! IEEE value for Inf (infinity)
   REAL(ReKi)                                :: Pi                            ! Ratio of a circle's circumference to its diameter
   REAL(ReKi)                                :: PiBy2                         ! Pi/2
   REAL(ReKi)                                :: R2D                           ! Factor to convert radians to degrees
   REAL(ReKi)                                :: RPM2RPS                       ! Factor to convert revolutions per minute to radians per second
   REAL(ReKi)                                :: RPS2RPM                       ! Factor to convert radians per second to revolutions per minute
   REAL(ReKi)                                :: TwoByPi                       ! 2/Pi
   REAL(ReKi)                                :: TwoPi                         ! 2*Pi

!bjj: I'm not sure why MLB added these 3 types; they don't seem to be used anywhere
   !TYPE, PUBLIC               :: CubSplineType                                ! This derived type is used to hold data for performing cubic splines.
   !   INTEGER                                :: NumPts                        ! The number of points in the XAry and YAry arrays.
   !   REAL(ReKi), ALLOCATABLE                :: Coef      (:,:)               ! The NumPts-1 length array of cubic coefficients.  The second dimension must be "0:3".
   !   REAL(ReKi), ALLOCATABLE                :: XAry      (:)                 ! The NumPts length array of x values for the interpolation.
   !   REAL(ReKi), ALLOCATABLE                :: YAry      (:)                 ! The NumPts length array of y values for the interpolation.
   !END TYPE CubSplineType
   !
   !TYPE, PUBLIC               :: RegCubSplineType                             ! This derived type is used to hold data for performing cubic splines wuth regularly-spaced data.
   !   INTEGER                                :: NumPts                        ! The number of points in the XAry and YAry arrays.
   !   REAL(ReKi), ALLOCATABLE                :: Coef      (:,:)               ! The NumPts-1 length array of cubic coefficients.  The second dimension must be "0:3".
   !   REAL(ReKi)                             :: DelX                          ! The distance between the equally spaced points in XAry.
   !   REAL(ReKi), ALLOCATABLE                :: XAry      (:)                 ! The NumPts length array of x values for the interpolation.
   !   REAL(ReKi), ALLOCATABLE                :: YAry      (:)                 ! The NumPts length array of y values for the interpolation.
   !END TYPE RegCubSplineType
   !
   !TYPE, PUBLIC               :: RegGridType                                  ! This derived type is used to hold the contents of a regular grid of data.
   !   INTEGER                                :: NumDims                       ! The number of dimensions for this grid.
   !   REAL(ReKi), ALLOCATABLE                :: Mins      (:)                 ! The set of minimums for the grid in each NumDims dimensions.
   !   REAL(ReKi), ALLOCATABLE                :: Steps     (:)                 ! The set of step sizes for the grid in each NumDims dimensions.
   !   REAL(ReKi), ALLOCATABLE                :: Grid      (:)                 ! The NumDims dimensional grid.
   !END TYPE RegGridType
!bjj: end of MLB types removed

!=======================================================================

      ! Create interface for a generic EqualRealNos that uses specific routines.

   INTERFACE EqualRealNos
      MODULE PROCEDURE EqualRealNos4
      MODULE PROCEDURE EqualRealNos8
      MODULE PROCEDURE EqualRealNos16
   END INTERFACE


      ! Create interface for a generic Eye that uses specific routines.

   INTERFACE Eye
      MODULE PROCEDURE Eye2   ! matrix of two dimensions
      MODULE PROCEDURE Eye3   ! matrix of three dimensions
   END INTERFACE


      ! Create interface for a generic InterpBin that actually uses specific routines.

   INTERFACE InterpBin
      MODULE PROCEDURE InterpBinComp
      MODULE PROCEDURE InterpBinReal
   END INTERFACE


      ! Create interface for a generic InterpStp that actually uses specific routines.

   INTERFACE InterpStp
      MODULE PROCEDURE InterpStpComp
      MODULE PROCEDURE InterpStpReal
   END INTERFACE




CONTAINS

!=======================================================================
   SUBROUTINE AddOrSub2Pi ( OldAngle, NewAngle )


      ! This routine is used to convert NewAngle to an angle within 2*Pi of
      !   OldAngle by adding or subtracting 2*Pi accordingly; it then sets
      !   OldAngle equal to NewAngle.  This routine is useful for converting
      !   angles returned from a call to the ATAN2() FUNCTION into angles that may
      !   exceed the -Pi to Pi limit of ATAN2().  For example, if the nacelle yaw
      !   angle was 179deg in the previous time step and the yaw angle increased
      !   by 2deg in the new time step, we want the new yaw angle returned from a
      !   call to the ATAN2() FUNCTION to be 181deg instead of -179deg.  This
      !   routine assumes that the angle change between calls is not more than
      !   2*Pi in absolute value.  OldAngle should be SAVEd in the calling
      !   routine.


      ! Argument declarations:

   REAL(ReKi), INTENT(INOUT)    :: OldAngle                                     ! Angle from which NewAngle will be converted to within 2*Pi of, rad.
   REAL(ReKi), INTENT(INOUT)    :: NewAngle                                     ! Angle to be converted to within 2*Pi of OldAngle, rad.


      ! Local declarations:

   REAL(ReKi)                   :: DelAngle                                     ! The difference between OldAngle and NewAngle, rad.



      ! Add or subtract 2*Pi in order to convert NewAngle two within 2*Pi of
      !   OldAngle:

   DelAngle = OldAngle - NewAngle

   DO WHILE ( ABS( DelAngle ) >= TwoPi )

      NewAngle = NewAngle + SIGN( TwoPi, DelAngle )
      DelAngle = OldAngle - NewAngle

   END DO


      ! Set OldAngle to equal NewAngle:

   OldAngle = NewAngle



   RETURN
   END SUBROUTINE AddOrSub2Pi
!=======================================================================
   SUBROUTINE BSortReal ( RealAry, NumPts )


      ! This routine sorts a list of real numbers.  It uses the bubble sort algorithm,
      ! which is only suitable for short lists.


      ! Argument declarations:

   INTEGER, INTENT(IN)          :: NumPts                                       ! The length of the list to be sorted.

   REAL(ReKi), INTENT(INOUT)    :: RealAry(NumPts)                              ! The list of real numbers to be sorted.


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
   FUNCTION Cross_Product(Vector1, Vector2)

      ! This function computes the cross product of two 3-element arrays:
      ! Cross_Product = Vector1 X Vector2 (resulting in a vector)


      ! Argument declarations.

   REAL(ReKi), INTENT(IN )         :: Vector1       (3)
   REAL(ReKi), INTENT(IN )         :: Vector2       (3)

      ! Function definition
   REAL(ReKi)                      :: Cross_Product (3)        ! = Vector1 X Vector2 (resulting in a vector)


   Cross_Product(1) = Vector1(2)*Vector2(3) - Vector1(3)*Vector2(2)
   Cross_Product(2) = Vector1(3)*Vector2(1) - Vector1(1)*Vector2(3)
   Cross_Product(3) = Vector1(1)*Vector2(2) - Vector1(2)*Vector2(1)


   RETURN
   END FUNCTION Cross_Product
!=======================================================================
   SUBROUTINE CubicSplineInit ( AryLen, XAry, YAry, Coef, ErrStat, ErrMsg )


      ! This routine calculates the parameters needed to compute a irregularly-spaced natural cubic spline.
      ! Natural cubic splines are used in that the curvature at the end points is zero.
      ! This routine does not require that the XAry be regularly spaced.


      ! Argument declarations:

   INTEGER, INTENT(IN)          :: AryLen                                     ! Length of the array.

   REAL(ReKi), INTENT(OUT)      :: Coef  (AryLen-1,0:3)                       ! The coefficients for the cubic polynomials.
   REAL(ReKi), INTENT(IN)       :: XAry  (AryLen)                             ! Input array of x values.
   REAL(ReKi), INTENT(IN)       :: YAry  (AryLen)                             ! Input array of y values.

   INTEGER(IntKi), INTENT(OUT)  :: ErrStat                                    ! Error status.

   CHARACTER(*), INTENT(OUT)    :: ErrMsg                                     ! Error message.


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
   SUBROUTINE CubicSplineInitM ( XAry, YAry, Coef, ErrStat, ErrMsg )


      ! This routine calculates the parameters needed to compute a irregularly-spaced natural cubic spline.
      ! Natural cubic splines are used in that the curvature at the end points is zero.
      ! This routine does not require that the XAry be regularly spaced.
      ! This version of the routine works with multiple curves that share the same X values.


      ! Argument declarations:

   REAL(ReKi), INTENT(OUT)      :: Coef  (:,:,0:)                             ! The coefficients for the cubic polynomials.
   REAL(ReKi), INTENT(IN)       :: XAry  (:)                                  ! Input array of x values.
   REAL(ReKi), INTENT(IN)       :: YAry  (:,:)                                ! Input array of y values with multiple curves.

   INTEGER(IntKi), INTENT(OUT)  :: ErrStat                                    ! Error status.

   CHARACTER(*), INTENT(OUT)    :: ErrMsg                                     ! Error message.


      ! Local declarations.

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
      Slope(I,:) = ( YAry(I+1,:) - YAry(I,:) )/DelX(I)
   END DO ! I


      ! Use Gaussian elimination to solve the tri-diagonal matrix.

   U(1  ) = 2.0_ReKi*( DelX (2)   + DelX (1  ) )
   V(1,:) = 6.0_ReKi*( Slope(2,:) - Slope(1,:) )

   DO I=2,NumPts-1
      U(I)   = 2.0_ReKi*( DelX (I-1) + DelX (I)     ) - DelX(I-1)*DelX(I-1  )/U(I-1)
      V(I,:) = 6.0_ReKi*( Slope(I,:) - Slope(I-1,:) ) - DelX(I-1)*   V(I-1,:)/U(I-1)
   END DO ! I


      ! Determine the coefficients of the polynomials.

   Coef(:,:,0) = YAry(1:NumPts-1,:)

   ZHi(:) = 0.0_ReKi

   DO I=NumPts-1,1,-1
      ZLo(:)      = ( V(I,:) - DelX(I)*ZHi(:) )/U(I)
      Coef(I,:,1) = Slope(I,:) - DelX(I)*( ZHi(:)/6.0_ReKi + ZLo(:)/3.0_ReKi )
      Coef(I,:,2) = 0.5_ReKi*ZLo(:)
      Coef(I,:,3) = ( ZHi(:) - ZLo(:) )/( 6.0_ReKi*DelX(I) )
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
   FUNCTION CubicSplineInterp ( X, AryLen, XAry, YAry, Coef, ErrStat, ErrMsg )


      ! This routine interpolates a pair of arrays using cubic splines to find the function value at X.
      ! One must call CubicSplineInit() first to compute the coefficients of the cubics.
      ! This routine does not require that the XAry be regularly spaced.


      ! Function declaration.

   REAL(ReKi)                   :: CubicSplineInterp                          ! This function.


      ! Argument declarations:

   INTEGER, INTENT(IN)          :: AryLen                                     ! Length of the array.

   REAL(ReKi), INTENT(IN)       :: Coef  (AryLen-1,0:3)                       ! The coefficients for the cubic polynomials.
   REAL(ReKi), INTENT(IN)       :: X                                          ! The value we are trying to interpolate for.
   REAL(ReKi), INTENT(IN)       :: XAry (AryLen)                              ! Input array of regularly spaced x values.
   REAL(ReKi), INTENT(IN)       :: YAry (AryLen)                              ! Input array of y values.

   INTEGER(IntKi), INTENT(OUT)  :: ErrStat                                    ! Error status.

   CHARACTER(*), INTENT(OUT)    :: ErrMsg                                     ! Error message.


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
   FUNCTION CubicSplineInterpM ( X, XAry, YAry, Coef, ErrStat, ErrMsg ) RESULT( Res )


      ! This routine interpolates a pair of arrays using cubic splines to find the function value at X.
      ! One must call CubicSplineInit() first to compute the coefficients of the cubics.
      ! This routine does not require that the XAry be regularly spaced.
      ! This version of the routine works with multiple curves that share the same X values.


      ! Function declaration.

   REAL(ReKi), ALLOCATABLE      :: Res(:)                                     ! The result of this function.


      ! Argument declarations:

   REAL(ReKi), INTENT(IN)       :: Coef  (:,:,0:)                             ! The coefficients for the cubic polynomials.
   REAL(ReKi), INTENT(IN)       :: X                                          ! The value we are trying to interpolate for.
   REAL(ReKi), INTENT(IN)       :: XAry (:)                                   ! Input array of regularly spaced x values.
   REAL(ReKi), INTENT(IN)       :: YAry (:,:)                                 ! Input array of y values with multiple curves.

   INTEGER(IntKi), INTENT(OUT)  :: ErrStat                                    ! Error status.

   CHARACTER(*),    INTENT(OUT) :: ErrMsg                                     ! Error message.


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
   FUNCTION DCM_exp(lambda)
   
      ! This function computes a matrix exponential.
      !
      ! "'Interpolation' of DCMs", M.A. Sprague, 11 March 2014, Eq. 31-33
      
   REAL(ReKi), INTENT(IN)  :: lambda(3)
   REAL(ReKi)              :: DCM_exp(3,3)
   
      ! local variables
   REAL(ReKi)              :: theta          ! angle of rotation   
   REAL(ReKi)              :: tmp_Mat(3,3)
   
   INTEGER(IntKi)          :: ErrStat
   CHARACTER(30)           :: ErrMsg  
   
   
   theta = TwoNorm(lambda)                   ! Eq. 32
   

   IF (EqualRealNos(theta, 0.0_ReKi) ) THEN    
      CALL eye(DCM_exp, ErrStat, ErrMsg)    ! Eq. 33a
   ELSE   
      
         ! convert lambda to skew-symmetric matrix:
      tmp_mat(1,1) =  0.0_ReKi                                            
      tmp_mat(2,1) = -lambda(3)                                           
      tmp_mat(3,1) =  lambda(2)                                           
      tmp_mat(1,2) =              lambda(3)                               
      tmp_mat(2,2) =              0.0_ReKi                                
      tmp_mat(3,2) =             -lambda(1)                               
      tmp_mat(1,3) =                               -lambda(2)             
      tmp_mat(2,3) =                                lambda(1)             
      tmp_mat(3,3) =                                0.0_ReKi            
      
         ! Eq. 33b
      CALL eye(DCM_exp, ErrStat, ErrMsg)                  
      DCM_exp = DCM_exp + sin(theta)/theta*tmp_mat + (1-cos(theta))/theta**2 * MATMUL(tmp_mat, tmp_mat) 
      
   END IF

   
      
   END FUNCTION DCM_exp
!=======================================================================  
   SUBROUTINE DCM_logMap(DCM, logMap, ErrStat, ErrMsg)

      ! This function computes the logarithmic map for a direction 
      ! cosine matrix.
      !
      ! "'Interpolation' of DCMs", M.A. Sprague, 11 March 2014, Eq. 24-30
      ! with eigenvector equations updated to account for numerics
   
   REAL(ReKi),       INTENT(IN)    :: DCM(3,3)
   REAL(ReKi),       INTENT(   OUT):: logMap(3)
   INTEGER(IntKi),   INTENT(  OUT) :: ErrStat                   ! Error status of the operation
   CHARACTER(*),     INTENT(  OUT) :: ErrMsg                    ! Error message if ErrStat /= ErrID_None
   
      ! local variables
   REAL(ReKi)                      :: temp
   REAL(ReKi)                      :: theta
   REAL(ReKi)                      :: v(3)
   REAL(ReKi)                      :: skewSym(3,3) ! an anti-symmetric matrix
      
         ! initialization
      ErrStat = ErrID_None
      ErrMsg  = ""   
   
   
      temp  = 0.5_ReKi*( trace(DCM) - 1.0_ReKi )
      temp  = min( max(temp,-1.0_ReKi), 1.0_ReKi ) !make sure it's in a valid range (to avoid cases where this is slightly outside the +/-1 range)
      theta = ACOS( temp )                                                      ! Eq. 25
   
      IF ( EqualRealNos(0.0_ReKi, theta) ) THEN
         logMap = 0.0_ReKi                                                   ! Eq. 26a
      ELSEIF ( EqualRealNos( pi, theta ) ) THEN
      
         ! calculate the eigenvector of DCM associated with eigenvalue +1:
      
         temp = -1.0_ReKi + DCM(2,2) + DCM(2,3)*DCM(3,2) + DCM(3,3) - DCM(2,2)*DCM(3,3)
         if ( .NOT. EqualRealNos(temp, 0.0_ReKi) ) then

            v(1) = 1.0_ReKi
            v(2) = -(DCM(2,1) + DCM(2,3)*DCM(3,1) - DCM(2,1)*DCM(3,3))/temp
            v(3) = -(DCM(3,1) - DCM(2,2)*DCM(3,1) + DCM(2,1)*DCM(3,2))/temp

         else 
            temp = -1.0_ReKi + DCM(1,1) + DCM(1,3)*DCM(3,1) + DCM(3,3) - DCM(1,1)*DCM(3,3)
            if ( .NOT. EqualRealNos(temp, 0.0_ReKi) ) then

               v(1) = -(DCM(1,2) + DCM(1,3)*DCM(3,2) - DCM(1,2)*DCM(3,3))/temp
               v(2) =  1.0_ReKi
               v(3) = -(DCM(3,2) + DCM(1,2)*DCM(3,1) - DCM(1,1)*DCM(3,2))/temp

            else 
               temp = -1.0_ReKi + DCM(1,1) + DCM(1,2)*DCM(2,1) + DCM(2,2) - DCM(1,1)*DCM(2,2)
               if ( .NOT. EqualRealNos(temp, 0.0_ReKi) ) then

                  v(1) = -(DCM(1,3) - DCM(1,3)*DCM(2,2) + DCM(1,2)*DCM(2,3))/temp
                  v(2) = -(DCM(1,3)*DCM(2,1) + DCM(2,3) - DCM(1,1)*DCM(2,3))/temp
                  v(3) = 1.0_ReKi
            
               else
                     ! break with error
                  ErrStat = ErrID_Fatal
                  WRITE( ErrMsg, '("DCM_logMap:invalid DCM matrix",3("'//Newline//'",4x,3(ES10.3E2,1x)))') DCM(1,:),DCM(2,:),DCM(3,:)               
                  RETURN
               end if         
            end if         
         endif
                         
            ! normalize the eigenvector:
         v = v / TwoNorm(v)                                                       ! Eq. 27                  
      
            ! calculate the skew-symmetric tensor (note we could change sign here for continuity)
         v =  pi*v                                                                ! Eq. 26c  
      
         logMap(1) = -v(1)
         logMap(2) =  v(2)
         logMap(3) = -v(3)
      
      ELSE ! 0 < theta < pi 
      
         skewSym = DCM - TRANSPOSE(DCM)
      
         logMap(1) = -skewSym(3,2)
         logMap(2) =  skewSym(3,1)
         logMap(3) = -skewSym(2,1)
      
         logMap = 0.5_ReKi * theta / sin(theta) * logMap   ! Eq. 26b
      END IF
   
   END SUBROUTINE DCM_logMap   
!=======================================================================
SUBROUTINE DCM_SetLogMapForInterp( tensor )

   ! this routine sets the rotation parameters (tensors from DCM_logMap)
   ! so that they can be appropriately interpolated, based on
   ! continunity of the neighborhood. The tensor input matrix has columns
   ! of rotational parameters; one column for each set of values to be 
   ! interpolated
   !
   ! This is based on the 2pi periodicity of rotations:
   ! if tensor is one solution to DCM_logMap( DCM ), then so is
   !  tensor*( 1 + TwoPi*k/TwoNorm(tensor) ) for any integer k
      
   
   REAL(ReKi),     INTENT(INOUT) :: tensor(:,:)

   REAL(ReKi)                    :: diff1, diff2      ! magnitude-squared of difference between two adjacent values
   REAL(ReKi)                    :: temp(3), temp1(3) ! difference between two tensors
   REAL(ReKi)                    :: period(3)         ! the period to add to the rotational parameters
   INTEGER(IntKi)                :: nc                ! size of the tensors matrix
   INTEGER(IntKi)                :: ic                ! loop counters for each array dimension
   
   nc = size(tensor,2)
          
      ! 
   do ic=2,nc      
      
      diff1 = TwoNorm( tensor(:,ic) )
      
      if ( .NOT. EqualRealNos( diff1, 0.0_ReKi) ) then
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
               temp  = temp1 + period !k=k-1; 
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
               temp  = temp1 + period !k=k-1; 
               diff2 = DOT_PRODUCT( temp, temp )
            end do
   
         end if
      
      end if ! tensor vector isn't zero=length
            
   end do
                 
END SUBROUTINE DCM_SetLogMapForInterp
!=======================================================================     
   FUNCTION EqualRealNos4 ( ReNum1, ReNum2 )

      ! This function compares 2 real numbers and determines if they
      ! are "almost" equal, i.e. within some relative tolerance
      ! ("Safe Comparisons" suggestion from http://www.lahey.com/float.htm)

      ! passed variables

   REAL(SiKi), INTENT(IN )         :: ReNum1                            ! the first  real number to compare
   REAL(SiKi), INTENT(IN )         :: ReNum2                            ! the second real number to compare

   LOGICAL                         :: EqualRealNos4                     ! the function definition -- returns .true. if the numbers are almost equal

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
   FUNCTION EqualRealNos8 ( ReNum1, ReNum2 )

      ! This function compares 2 real numbers and determines if they
      ! are "almost" equal, i.e. within some relative tolerance
      ! ("Safe Comparisons" suggestion from http://www.lahey.com/float.htm)

      ! passed variables

   REAL(R8Ki), INTENT(IN )         :: ReNum1                            ! the first  real number to compare
   REAL(R8Ki), INTENT(IN )         :: ReNum2                            ! the second real number to compare

   LOGICAL                         :: EqualRealNos8                     ! the function definition -- returns .true. if the numbers are almost equal

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
   FUNCTION EqualRealNos16 ( ReNum1, ReNum2 )

      ! This function compares 2 real numbers and determines if they
      ! are "almost" equal, i.e. within some relative tolerance
      ! ("Safe Comparisons" suggestion from http://www.lahey.com/float.htm)

      ! passed variables

   REAL(QuKi), INTENT(IN )         :: ReNum1                            ! the first  real number to compare
   REAL(QuKi), INTENT(IN )         :: ReNum2                            ! the second real number to compare

   LOGICAL                         :: EqualRealNos16                    ! the function definition -- returns .true. if the numbers are almost equal

      ! local variables
   REAL(QuKi), PARAMETER           :: Eps = EPSILON(ReNum1)             ! machine precision
   REAL(QuKi), PARAMETER           :: Tol = 100.0_QuKi*Eps / 2.0_QuKi   ! absolute tolerance (ignore the last 2 significant digits)

   REAL(QuKi)                      :: Fraction


      ! make sure we're never trying to get more precision than Tol

   Fraction = MAX( ABS(ReNum1+ReNum2), 1.0_QuKi )



      ! determine if ReNum1 and ReNum2 are approximately equal

   IF ( ABS(ReNum1 - ReNum2) <= Fraction*Tol ) THEN  ! the relative error
      EqualRealNos16 = .TRUE.
   ELSE
      EqualRealNos16 = .FALSE.
   ENDIF


   END FUNCTION EqualRealNos16
!=======================================================================
   FUNCTION EulerConstruct(theta) result(M)
   
      ! this function creates a rotation matrix, M, from a 1-2-3 rotation
      ! sequence of the 3 Euler angles, theta_x, theta_y, and theta_z.
      ! M represents a change of basis (from global to local coordinates; 
      ! not a physical rotation of the body). it is the inverse of EulerExtract().
      !
      ! M = R(theta_z) * R(theta_y) * R(theta_x)
      !   = [ cz sz 0 |   [ cy  0 -sy |   [ 1   0   0 |
      !     |-sz cz 0 | * |  0  1   0 | * | 0  cx  sx |
      !     |  0  0 1 ]   | sy  0  cy ]   | 0 -sx  cx ]
      !   = [ cy*cz   cx*sz+sx*sy*cz    sx*sz-cx*sy*cz |
      !     |-cy*sz   cx*cz-sx*sy*sz    sx*cz+cx*sy*sz |
      !     | sy           -sx*cy             cx*cy    ]
      ! where cz = cos(theta_z), sz = sin(theta_z), cy = cos(theta_y), etc.
   
      REAL(ReKi)             :: M(3,3)    ! rotation matrix M 
      REAL(ReKi), INTENT(IN) :: theta(3)  ! the 3 rotation angles: theta_x, theta_y, theta_z
      
      REAL(ReKi)             :: cx        ! cos(theta_x)
      REAL(ReKi)             :: sx        ! sin(theta_x)
      REAL(ReKi)             :: cy        ! cos(theta_y)
      REAL(ReKi)             :: sy        ! sin(theta_y)
      REAL(ReKi)             :: cz        ! cos(theta_z)
      REAL(ReKi)             :: sz        ! sin(theta_z)
   

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
   
   END FUNCTION EulerConstruct
!=======================================================================
   FUNCTION EulerExtract(M) result(theta)
   
      ! if M is a rotation matrix from a 1-2-3 rotation sequence, this function
      ! returns the 3 Euler angles, theta_x, theta_y, and theta_z, that formed 
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
   
      REAL(ReKi), INTENT(IN) :: M(3,3)    ! rotation matrix M 
      REAL(ReKi)             :: theta(3)  ! the 3 rotation angles: theta_x, theta_y, theta_z
      
      REAL(ReKi)             :: cx        ! cos(theta_x)
      REAL(ReKi)             :: sx        ! sin(theta_x)
      REAL(ReKi)             :: cy        ! cos(theta_y)
      REAL(ReKi)             :: sy        ! sin(theta_y)
      REAL(ReKi)             :: cz        ! cos(theta_z)
      REAL(ReKi)             :: sz        ! sin(theta_z)
   
         ! use trig identity sx**2 + cx**2 = 1 to get cy:
      cy = sqrt( m(1,1)**2 + m(2,1)**2 ) 
            
      if ( EqualRealNos(cy,0.0_ReKi) ) then
      !if ( cy < 16*epsilon(0.0_ReKi) ) then
         
         theta(2) = atan2( m(3,1), cy )               ! theta_y
         
         ! cy = 0 -> sy = +/- 1
         ! M  = [  0   cx*sz+/-sx*cz    sx*sz-/+cx*cz |
         !      |  0   cx*cz-/+sx*sz    sx*cz+/-cx*sz |
         !      |+/-1        0                0       ]
         
         ! gimbal lock allows us to choose theta_z = 0
         theta(3) = 0.0_ReKi                          ! theta_z
         
         ! which reduces the matrix to 
         ! M  = [  0  +/-sx  -/+cx |
         !      |  0     cx     sx |
         !      |+/-1    0       0 ]
         
         theta(1) = atan2(  m(2,3), m(2,2) )          ! theta_x
      else

         theta(3) = atan2( -m(2,1), m(1,1) )          ! theta_z         
         cz       = cos( theta(3) )
         sz       = sin( theta(3) )

            ! get the appropriate sign for cy:
         if ( EqualRealNos(cz, 0.0_ReKi) ) then
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
            
      
   END FUNCTION EulerExtract
!=======================================================================
   SUBROUTINE Eye2( A, ErrStat, ErrMsg )

      ! This routine sets the matrix A(:,:) to the identity
      ! matrix (all zeros, with ones on the diagonal)
      ! Note that this also returns the "pseudo-identity" when A(:,:)
      ! is not square (i.e., nr/=nc).

   REAL(ReKi),     INTENT(INOUT) :: A (:,:)                        ! Array to matricies to set to the identity matrix (nr,nc,n)
   INTEGER(IntKi), INTENT(OUT)   :: ErrStat                        ! Error level
   CHARACTER(*),   INTENT(OUT)   :: ErrMsg                         ! ErrMsg corresponding to ErrStat

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
   A = 0._ReKi

      ! set the diagonals to one:
   DO j = 1, MIN(nr,nc) ! the diagonal of the matrix
      A(j,j) = 1._ReKi
   END DO

   END SUBROUTINE Eye2
!=======================================================================
   SUBROUTINE Eye3( A, ErrStat, ErrMsg )

      ! This routine sets each of the n matries A(:,:,n) to the identity
      ! matrix (all zeros, with ones on the diagonal).
      ! Note that this also returns the "pseudo-identity" when A(:,:)
      ! is not square (i.e., nr/=nc).

   REAL(ReKi),     INTENT(INOUT) :: A (:,:,:)                      ! Array to matricies to set to the identity matrix (nr,nc,n)
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
   A = 0._ReKi

      ! set the diagonals to one:
   DO i = 1, n ! loop through the matrices
      DO j = 1, MIN(nr,nc) ! the diagonal of the matrix
         A(j,j,i) = 1._ReKi
      END DO
   END DO

   END SUBROUTINE Eye3
!=======================================================================
   SUBROUTINE GaussElim( AugMatIn, NumEq, x, ErrStat, ErrMsg )

      ! This routine uses the Gauss-Jordan elimination method for the
      !   solution of a given set of simultaneous linear equations.
      ! NOTE: this routine works if no pivot points are zero and you
      !   don't want the eschelon or reduced eschelon form of the
      !   augmented matrix.  The form of the original augmented matrix
      !   IS preserved in this call.
      ! This routine was originally in FAST.f90.
      ! When AugMatIn = [ A b ], this routine returns the solution
      ! vector x to the equation Ax = b.

   IMPLICIT                        NONE


      ! Passed variables:

   INTEGER(IntKi), INTENT(IN )  :: NumEq                                           ! Number of equations in augmented matrix

   REAL(ReKi),     INTENT(IN )  :: AugMatIn (NumEq, NumEq+1 )                      ! Augmented matrix passed into this subroutine ( AugMatIn = [ A b ]
   REAL(ReKi),     INTENT(OUT)  :: x (NumEq)                                       ! Solution vector

   INTEGER(IntKi), INTENT(OUT)  :: ErrStat                                         ! Error level
   CHARACTER(*),   INTENT(OUT)  :: ErrMsg                                          ! ErrMsg corresponding to ErrStat


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
   SUBROUTINE GetOffsetReg ( Ary, NumPts, Val, Ind, Fract, ErrStat, ErrMsg )


      ! Determine index of the point in Ary just below Val and the fractional distance to the next point in the array.
      ! The elements of the array are assumed to be regularly spaced.


      ! Argument declarations:

   INTEGER, INTENT(IN)          :: NumPts                                     ! Length of the array.

   REAL(ReKi), INTENT(IN)       :: Ary  (NumPts)                              ! Input array of regularly spaced values.
   REAL(ReKi), INTENT(OUT)      :: Fract                                      ! The fractional distance of Val between the surrounding array elements.
   REAL(ReKi), INTENT(IN)       :: Val                                        ! The value we hope to bound in the array.

   INTEGER(IntKi), INTENT(OUT)  :: ErrStat                                    ! Error status.
   INTEGER(IntKi), INTENT(OUT)  :: Ind                                        ! The index of the point in Ary just below Val.

   CHARACTER(*),   INTENT(OUT)  :: ErrMsg                                     ! Error message.


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
   FUNCTION GetSmllRotAngs ( DCMat, ErrStat, ErrMsg )

      ! This subroutine computes the angles that make up the input direction cosine matrix, DCMat
      ! It is the inverse of SmllRotTrans()
      
      
      ! passed variables

   REAL(ReKi), INTENT(IN )            :: DCMat          (3,3)
   INTEGER,    INTENT(OUT )           :: ErrStat               ! a non-zero value indicates an error in the permutation matrix algorithm
   CHARACTER(*),INTENT(OUT ),OPTIONAL :: ErrMsg                ! a non-zero value indicates an error in the permutation matrix algorithm

   REAL(ReKi)                         :: GetSmllRotAngs ( 3 )

      ! local variables
   REAL(ReKi)                         :: denom                 ! the denominator of the resulting matrix
   REAL(ReKi), PARAMETER              :: LrgAngle  = 0.4       ! Threshold for when a small angle becomes large (about 23deg).  This comes from: COS(SmllAngle) ~ 1/SQRT( 1 + SmllAngle^2 ) and SIN(SmllAngle) ~ SmllAngle/SQRT( 1 + SmllAngle^2 ) results in ~5% error when SmllAngle = 0.4rad.



      ! initialize output angles (just in case there is an error that prevents them from getting set)

   GetSmllRotAngs = 0.0
   ErrStat        = ErrID_None
   ErrMsg         = ""

      ! calculate the small angles
   GetSmllRotAngs(1) = DCMat(2,3) - DCMat(3,2)
   GetSmllRotAngs(2) = DCMat(3,1) - DCMat(1,3)
   GetSmllRotAngs(3) = DCMat(1,2) - DCMat(2,1)

   denom             = DCMat(1,1) + DCMat(2,2) + DCMat(3,3) - 1

   IF ( .NOT. EqualRealNos( denom, 0.0_ReKi ) ) THEN
      GetSmllRotAngs = GetSmllRotAngs / denom

         ! check that the angles are, in fact, small
      IF ( ANY( ABS(GetSmllRotAngs) > LrgAngle ) ) THEN
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


   END FUNCTION GetSmllRotAngs ! ( DCMat, PMat, ErrStat [, ErrMsg] )
!=======================================================================
   SUBROUTINE GL_Pts ( IPt, NPts, Loc, Wt, ErrStat )

      ! This funtion returns the non-dimensional (-1:+1) location of the given Gauss-Legendre Quadrature point and its weight.
      ! The values came from Carnahan, Brice; Luther, H.A.; Wilkes, James O.  (1969)  "Applied Numerical Methods."


      ! Argument declarations.

   REAL(ReKi)                     :: Loc                                         ! The location of the specified point.
   REAL(ReKi)                     :: Wt                                          ! The weight for the specified point.

   INTEGER, INTENT(OUT), OPTIONAL :: ErrStat                                     ! Error status; if present, program does not abort on error
   INTEGER, INTENT(INOUT)         :: IPt                                         ! The quadrature point in question.
   INTEGER, INTENT(INOUT)         :: NPts                                        ! The number of points used in the quadrature.


   IF ( PRESENT(ErrStat) ) ErrStat = 0


      ! Check to see if the number of points and the specific point are valid values.

   IF ( ( NPts < 1 ) .OR. ( NPts > 6 ) )  THEN
      CALL ProgAbort ( ' In function GL_Loc, the number of points used for Gauss-Legendre Quadrature must be between 1 and 6' &
                    //' (inclusive).  Instead, it is "'//TRIM( Int2LStr( NPts ) )//'".', PRESENT(ErrStat) )
      IF ( PRESENT(ErrStat) ) THEN ! this should always be true here
         ErrStat = ErrID_Fatal
         RETURN
      END IF
   END IF

   IF ( ( Ipt < 1 ) .OR. ( Ipt > NPts ) )  THEN
      CALL ProgAbort ( ' In function GL_Loc, the point being used for Gauss-Legendre Quadrature must be between 1 and ' &
                   //TRIM( Int2LStr( NPts ) )//' (inclusive).  Instead, it is "'//TRIM( Int2LStr( Ipt ) )//'".', PRESENT(ErrStat) )
      IF ( PRESENT(ErrStat) ) THEN
         ErrStat = ErrID_Fatal
         RETURN
      END IF
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
   FUNCTION IndexCharAry( CVal, CAry )


      ! This funtion returns an integer index such that CAry(IndexCharAry) = CVal. If
      ! no element in the array matches CVal, the value -1 is returned.  The routine
      ! performs a binary search on the input array to determine if CVal is an
      ! element of the array; thus, CAry must be sorted and stored in increasing
      ! alphebetical (ASCII) order. The routine does not check that the array is
      ! sorted.  The routine assumes that CVal is type CHARACTER and CAry
      ! is an array of CHARACTERS.


      ! Function declaration.


   INTEGER                      :: IndexCharAry                                   ! This function

      ! Argument declarations.

   CHARACTER(*), INTENT(IN)     :: CVal                                           ! String to find.
   CHARACTER(*), INTENT(IN)     :: CAry(:)                                        ! Array of strings to search.



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
   FUNCTION InterpBinComp( XVal, XAry, YAry, ILo, AryLen )


      ! This funtion returns a y-value that corresponds to an input x-value by interpolating into the arrays.
      ! It uses a binary interpolation scheme that takes about log(AryLen)/log(2) steps to converge.
      ! It returns the first or last YAry() value if XVal is outside the limits of XAry().
      ! This routine assumes YAry is COMPLEX.


      ! Function declaration.


   COMPLEX(ReKi)                :: InterpBinComp                                   ! This function.


      ! Argument declarations.

   INTEGER, INTENT(IN)          :: AryLen                                          ! Length of the arrays.
   INTEGER, INTENT(INOUT)       :: ILo                                             ! The low index into the arrays.

   REAL(ReKi), INTENT(IN)       :: XAry    (AryLen)                                ! Array of X values to be interpolated.
   REAL(ReKi), INTENT(IN)       :: XVal                                            ! X value to be interpolated.

   COMPLEX(ReKi), INTENT(IN)    :: YAry    (AryLen)                                ! Array of Y values to be interpolated.


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
   FUNCTION InterpBinReal( XVal, XAry, YAry, ILo, AryLen )


      ! This funtion returns a y-value that corresponds to an input x-value by interpolating into the arrays.
      ! It uses a binary interpolation scheme that takes about log(AryLen)/log(2) steps to converge.
      ! It returns the first or last YAry() value if XVal is outside the limits of XAry().
      ! This routine assumes YAry is REAL.


      ! Function declaration.


   REAL(ReKi)                   :: InterpBinReal                                   ! This function.


      ! Argument declarations.

   INTEGER, INTENT(IN)          :: AryLen                                          ! Length of the arrays.
   INTEGER, INTENT(INOUT)       :: ILo                                             ! The low index into the arrays.

   REAL(ReKi), INTENT(IN)       :: XAry    (AryLen)                                ! Array of X values to be interpolated.
   REAL(ReKi), INTENT(IN)       :: XVal                                            ! X value to be interpolated.
   REAL(ReKi), INTENT(IN)       :: YAry    (AryLen)                                ! Array of Y values to be interpolated.


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
   FUNCTION InterpStpComp( XVal, XAry, YAry, Ind, AryLen )


      ! This funtion returns a y-value that corresponds to an input x-value by interpolating into the arrays.
      ! It uses the passed index as the starting point and does a stepwise interpolation from there.  This is
      ! especially useful when the calling routines save the value from the last time this routine was called
      ! for a given case where XVal does not change much from call to call.  When there is no correlation
      ! from one interpolation to another, InterpBin() may be a better choice.
      ! It returns the first or last YAry() value if XVal is outside the limits of XAry().
      ! This routine assumes YAry is COMPLEX.


      ! Function declaration.


   COMPLEX(ReKi)                :: InterpStpComp                                   ! This function.


      ! Argument declarations.

   INTEGER, INTENT(IN)          :: AryLen                                          ! Length of the arrays.
   INTEGER, INTENT(INOUT)       :: Ind                                             ! Initial and final index into the arrays.

   REAL(ReKi), INTENT(IN)       :: XAry    (AryLen)                                ! Array of X values to be interpolated.
   REAL(ReKi), INTENT(IN)       :: XVal                                            ! X value to be interpolated.

   COMPLEX(ReKi), INTENT(IN)    :: YAry    (AryLen)                                ! Array of Y values to be interpolated.



      ! Let's check the limits first.

   IF ( XVal <= XAry(1) )  THEN
      InterpStpComp = YAry(1)
      Ind           = 1
      RETURN
   ELSE IF ( XVal >= XAry(AryLen) )  THEN
      InterpStpComp = YAry(AryLen)
      Ind           = MAX(AryLen - 1, 1)
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

         InterpStpComp = ( YAry(Ind+1) - YAry(Ind) )*( XVal - XAry(Ind) )/( XAry(Ind+1) - XAry(Ind) ) + YAry(Ind)
         RETURN

      END IF

   END DO


   RETURN
   END FUNCTION InterpStpComp ! ( XVal, XAry, YAry, Ind, AryLen )
!=======================================================================
   FUNCTION InterpStpReal( XVal, XAry, YAry, Ind, AryLen )


      ! This funtion returns a y-value that corresponds to an input x-value by interpolating into the arrays.
      ! It uses the passed index as the starting point and does a stepwise interpolation from there.  This is
      ! especially useful when the calling routines save the value from the last time this routine was called
      ! for a given case where XVal does not change much from call to call.  When there is no correlation
      ! from one interpolation to another, InterpBin() may be a better choice.
      ! It returns the first or last YAry() value if XVal is outside the limits of XAry().
      ! This routine assumes YAry is REAL.


      ! Function declaration.

   REAL(ReKi)                   :: InterpStpReal                                   ! This function.


      ! Argument declarations.

   INTEGER, INTENT(IN)          :: AryLen                                          ! Length of the arrays.
   INTEGER, INTENT(INOUT)       :: Ind                                             ! Initial and final index into the arrays.

   REAL(ReKi), INTENT(IN)       :: XAry    (AryLen)                                ! Array of X values to be interpolated.
   REAL(ReKi), INTENT(IN)       :: XVal                                            ! X value to be interpolated.
   REAL(ReKi), INTENT(IN)       :: YAry    (AryLen)                                ! Array of Y values to be interpolated.



      ! Let's check the limits first.

   IF ( XVal <= XAry(1) )  THEN
      InterpStpReal = YAry(1)
      Ind           = 1
      RETURN
   ELSE IF ( XVal >= XAry(AryLen) )  THEN
      InterpStpReal = YAry(AryLen)
      Ind           = MAX(AryLen - 1, 1)
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

         InterpStpReal = ( YAry(Ind+1) - YAry(Ind) )*( XVal - XAry(Ind) )/( XAry(Ind+1) - XAry(Ind) ) + YAry(Ind)
         RETURN

      END IF

   END DO


   RETURN
   END FUNCTION InterpStpReal ! ( XVal, XAry, YAry, Ind, AryLen )

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
   FUNCTION InterpWrappedStpReal( XValIn, XAry, YAry, Ind, AryLen )


      ! This funtion returns a y-value that corresponds to an input x-value which is wrapped back
      ! into the range [0-XAry(AryLen) by interpolating into the arrays.  
      ! It is assumed that XAry is sorted in ascending order.
      ! It uses the passed index as the starting point and does a stepwise interpolation from there.  This is
      ! especially useful when the calling routines save the value from the last time this routine was called
      ! for a given case where XVal does not change much from call to call.  When there is no correlation
      ! from one interpolation to another, InterpBin() may be a better choice.
      ! It returns the first or last YAry() value if XVal is outside the limits of XAry().
      ! This routine assumes YAry is REAL.


      ! Function declaration.

   REAL(ReKi)                   :: InterpWrappedStpReal                                   ! This function.


      ! Argument declarations.

   INTEGER, INTENT(IN)          :: AryLen                                          ! Length of the arrays.
   INTEGER, INTENT(INOUT)       :: Ind                                             ! Initial and final index into the arrays.

   REAL(ReKi), INTENT(IN)       :: XAry    (AryLen)                                ! Array of X values to be interpolated.
   REAL(ReKi), INTENT(IN)       :: XValIn                                           ! X value to be interpolated.
   REAL(ReKi), INTENT(IN)       :: YAry    (AryLen)                                ! Array of Y values to be interpolated.

   REAL(ReKi)                   :: XVal                                           ! X value to be interpolated.
   
   
   
      ! Wrap XValIn into the range XAry(1) to XAry(AryLen)
   XVal = MOD(XValIn, XAry(AryLen))

      ! Set the Ind to the first index if we are at the beginning of XAry
   IF ( XVal <= XAry(2) )  THEN  
      Ind           = 1
   END IF
   
   InterpWrappedStpReal = InterpStpReal( XVal, XAry, YAry, Ind, AryLen )
   
   
   END FUNCTION InterpWrappedStpReal ! ( XVal, XAry, YAry, Ind, AryLen )
!=======================================================================
!> This subroutine calculates the iosparametric coordinates, isopc, which is a value between -1 and 1 
!! (for each dimension of a dataset), indicating where InCoord falls between posLo and posHi.
!! It is used in InterpStpReal2D and InterpStpReal3D.
   SUBROUTINE IsoparametricCoords( InCoord, posLo, posHi, isopc )

! This subroutine calculates the iosparametric coordinates, isopc, which is a value between -1 and 1 
! (for each dimension of a dataset), indicating where InCoord falls between posLo and posHi.
! It is used in InterpStpReal2D and InterpStpReal3D.
   
   
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
   FUNCTION IsSymmetric( A )

      ! This function returns a logical TRUE/FALSE value that indicates
      ! if the given (2-dimensional) matrix, A, is symmetric. If A is not
      ! square it returns FALSE.


         ! passed variables

      REAL(ReKi), INTENT(IN) :: A(:,:)                   ! a real matrix A, whose symmetry is questioned
      LOGICAL                :: IsSymmetric              ! true if A is symmetric, false if not

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
   SUBROUTINE LocateBin( XVal, XAry, Ind, AryLen )

      ! This subroutine finds the lower-bound index of an input x-value located in an array.
      ! On return, Ind has a value such that
      !           XAry(Ind) <= XVal < XAry(Ind+1), with the exceptions that
      !             Ind = 0 when XVal < XAry(1), and
      !          Ind = AryLen when XAry(AryLen) <= XVal.
      !
      ! It uses a binary interpolation scheme that takes about log(AryLen)/log(2) steps to converge.
      ! If the index doesn't change much between calls, LocateStp() may be a better option.


      ! Argument declarations.

   INTEGER, INTENT(IN)          :: AryLen                                          ! Length of the array.
   INTEGER, INTENT(OUT)         :: Ind                                             ! Final (low) index into the array.

   REAL(ReKi), INTENT(IN)       :: XAry    (AryLen)                                ! Array of X values to be interpolated.
   REAL(ReKi), INTENT(IN)       :: XVal                                            ! X value to be interpolated.


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
   SUBROUTINE LocateStp( XVal, XAry, Ind, AryLen )

      ! This subroutine finds the lower-bound index of an input x-value located in an array.
      ! On return, Ind has a value such that
      !           XAry(Ind) <= XVal < XAry(Ind+1), with the exceptions that
      !             Ind = 0 when XVal < XAry(1), and
      !          Ind = AryLen when XAry(AryLen) <= XVal.
      !
      ! It uses the passed index as the starting point and does a stepwise search from there.  This is
      ! especially useful when the calling routines save the value from the last time this routine was called
      ! for a given case where XVal does not change much from call to call.  When there is no correlation
      ! from one interpolation to another, a binary search may be a better choice.



      ! Argument declarations.

   INTEGER, INTENT(IN)          :: AryLen                                          ! Length of the array.
   INTEGER, INTENT(INOUT)       :: Ind                                             ! Initial and final index into the array.

   REAL(ReKi), INTENT(IN)       :: XAry    (AryLen)                                ! Array of X values to be interpolated.
   REAL(ReKi), INTENT(IN)       :: XVal                                            ! X value to be interpolated.



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

   END SUBROUTINE LocateStp
!=======================================================================
   FUNCTION Mean ( Ary, AryLen )

      ! This routine calculates the mean value of an array.

   !NOTE: We should make AryLen an optional argument and use SIZE( Ary ) if it is not present.

      ! Function declaration.

   REAL(ReKi)                   :: Mean                                         ! This function.


      ! Argument declarations:

   INTEGER, INTENT(IN)          :: AryLen                                       ! Length of the array.

   REAL(ReKi), INTENT(IN)       :: Ary  (AryLen)                                ! Input array.


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
   SUBROUTINE MPi2Pi ( Angle )

      ! This routine is used to convert Angle to an equivalent value
      !  between -pi and pi.
                 
      ! Argument declarations:

   REAL(ReKi), INTENT(INOUT)    :: Angle



      ! Get the angle between 0 and 2Pi.

   Angle = MODULO( Angle, TwoPi )


      ! Get the angle between -Pi and Pi.

   IF ( Angle > Pi )  THEN
      Angle = Angle - TwoPi
   END IF


   RETURN
   END SUBROUTINE MPi2Pi
!=======================================================================
   FUNCTION OuterProduct(vec1,vec2)
   
   ! this routine calculates the outer product of two vectors

   REAL(ReKi),INTENT(IN):: vec1(:),vec2(:)
   REAL(ReKi)::OuterProduct(SIZE(vec1),SIZE(vec2))

   INTEGER(IntKi)::i,j,n1,n2

   n1=SIZE(vec1)
   n2=SIZE(vec2)

   DO i=1,n1
       DO j=1,n2
           OuterProduct(i,j) = vec1(i) * vec2(j)
       ENDDO
   ENDDO

   END FUNCTION OuterProduct   
!=======================================================================
   FUNCTION PSF ( Npsf, NumPrimes, subtract )

    ! This routine factors the number N into its primes.  If any of those
    ! prime factors is greater than the NumPrimes'th prime, a value of 1
    ! is added to N and the new number is factored.  This process is 
    ! repeated until no prime factors are greater than the NumPrimes'th 
    ! prime.
    !
    ! If subract is .true., we will subtract 1 from the value of N instead
    ! of adding it.

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
   FUNCTION Quaternion_Conjugate(q)

      ! This function computes the conjugate of a quaternion, q
      !
      ! "'Interpolation' of DCMs", M.A. Sprague, 11 March 2014, Eq. 6
   
   TYPE(Quaternion), INTENT(IN)    :: q     
   
   TYPE(Quaternion)                :: Quaternion_Conjugate
   
      
   Quaternion_Conjugate%q0 =  q%q0 
   Quaternion_Conjugate%v  = -q%v
      
   END FUNCTION Quaternion_Conjugate   
!=======================================================================  
   FUNCTION Quaternion_Norm(q)

      ! This function computes the 2-norm of a quaternion, q
      !
      ! "'Interpolation' of DCMs", M.A. Sprague, 11 March 2014, Eq. 5
   
   TYPE(Quaternion), INTENT(IN)    :: q     
   
   REAL(ReKi)                      :: Quaternion_Norm
   
      
   Quaternion_Norm = sqrt( q%q0**2 + DOT_PRODUCT(q%v, q%v) )
   
   
   END FUNCTION Quaternion_Norm   
!=======================================================================  
   FUNCTION Quaternion_Power(q,alpha)

      ! This function computes the quaternion, q, raised to an arbitrary
      ! real exponent, alpha.
      !
      ! "'Interpolation' of DCMs", M.A. Sprague, 11 March 2014, Eq. 7-8
   
   TYPE(Quaternion), INTENT(IN)    :: q     
   REAL(ReKi)      , INTENT(IN)    :: alpha
   
   TYPE(Quaternion)                :: Quaternion_Power
   
   
      ! local variables
   REAL(ReKi)                      :: greek   ! the product of alpha and theta
   REAL(ReKi)                      :: n(3)
   REAL(ReKi)                      :: q_norm
   REAL(ReKi)                      :: q_norm_power
   REAL(ReKi)                      :: theta
      
   
   q_norm       = Quaternion_Norm( q )     
   theta        = acos( q%q0 / q_norm )
   n            = q%v / TwoNorm(q%v)
   
   greek        = alpha * theta
   q_norm_power = q_norm ** alpha
   
   Quaternion_Power%q0 =  q_norm_power * cos( greek )
   Quaternion_Power%v  =  q_norm_power * sin( greek ) * n
      
   END FUNCTION Quaternion_Power   
!=======================================================================  
   FUNCTION Quaternion_Product(p, q)

      ! This function computes the product of two quaternions, p and q
      !
      ! "'Interpolation' of DCMs", M.A. Sprague, 11 March 2014, Eq. 4
   
   TYPE(Quaternion), INTENT(IN)    :: p      
   TYPE(Quaternion), INTENT(IN)    :: q     
   
   TYPE(Quaternion)                :: Quaternion_Product
   
      
   Quaternion_Product%q0 = p%q0 * q%q0 - DOT_PRODUCT(p%v, q%v)
   Quaternion_Product%v  = p%q0*q%v + q%q0*p%v + CROSS_PRODUCT( p%v, q%v ) 
   
   
   END FUNCTION Quaternion_Product   
!=======================================================================  
   FUNCTION Quaternion_to_DCM(q)

      ! This function converts a quaternion to an equivalent direction cosine matrix
      !
      ! "'Interpolation' of DCMs", M.A. Sprague, 11 March 2014, Eq. 9-17
   
   TYPE(Quaternion), INTENT(IN)    :: q     
   
   REAL(ReKi)                      :: Quaternion_to_DCM (3,3)
   
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
   FUNCTION DCM_to_Quaternion(DCM)

      ! This function converts a direction cosine matrix to an equivalent quaternion 
      !
      ! "'Interpolation' of DCMs", M.A. Sprague, 11 March 2014, Eq. 18-21
   
   REAL(ReKi)      , INTENT(IN)    :: DCM (3,3)          ! Direction cosine matrix
   TYPE(Quaternion)                :: DCM_to_Quaternion        
   
         
   DCM_to_Quaternion%q0   =      0.5_ReKi * sqrt( 1.0_ReKi + DCM(1,1) + DCM(2,2) + DCM(3,3) )                         ! Eq. 18
   DCM_to_Quaternion%v(1) = sign(0.5_ReKi * sqrt( 1.0_ReKi + DCM(1,1) - DCM(2,2) - DCM(3,3) ) , DCM(2,3) - DCM(3,2) ) ! Eq. 19
   DCM_to_Quaternion%v(2) = sign(0.5_ReKi * sqrt( 1.0_ReKi - DCM(1,1) + DCM(2,2) - DCM(3,3) ) , DCM(3,1) - DCM(1,3) ) ! Eq. 20
   DCM_to_Quaternion%v(3) = sign(0.5_ReKi * sqrt( 1.0_ReKi - DCM(1,1) - DCM(2,2) + DCM(3,3) ) , DCM(1,2) - DCM(2,1) ) ! Eq. 21

   
   
   END FUNCTION DCM_to_Quaternion
!=======================================================================         
   FUNCTION Quaternion_Interp(q1,q2,s)

      ! This function computes the interpolated quaternion at time
      ! t1 + s*(t2-t1) and s is in [0,1]
      ! 
      ! "'Interpolation' of DCMs", M.A. Sprague, 11 March 2014, Eq. 23
   
   TYPE(Quaternion), INTENT(IN)    :: q1      
   TYPE(Quaternion), INTENT(IN)    :: q2    
   REAL(ReKi),       INTENT(IN)    :: s
   
   TYPE(Quaternion)                :: Quaternion_Interp
   
      
   Quaternion_Interp = Quaternion_Conjugate(q1)
   Quaternion_Interp = Quaternion_Product(Quaternion_Interp, q2)
   Quaternion_Interp = Quaternion_Power(  Quaternion_Interp, s )
   Quaternion_Interp = Quaternion_Product(q1, Quaternion_Interp)
   
   
! bjj: this function has not been tested. I have not tested any of the quaternion routines, either. 

   END FUNCTION Quaternion_Interp
!=======================================================================
   SUBROUTINE RegCubicSplineInit ( AryLen, XAry, YAry, DelX, Coef, ErrStat, ErrMsg )


      ! This routine calculates the parameters needed to compute a regularly-spaced natural cubic spline.
      ! Natural cubic splines are used in that the curvature at the end points is zero.
      ! It assumes the XAry values are equally spaced for speed.


      ! Argument declarations:

   INTEGER, INTENT(IN)          :: AryLen                                     ! Length of the array.

   REAL(ReKi), INTENT(OUT)      :: Coef  (AryLen-1,0:3)                       ! The coefficients for the cubic polynomials.
   REAL(ReKi), INTENT(OUT)      :: DelX                                       ! The distance between the equally spaced points.
   REAL(ReKi), INTENT(IN)       :: XAry  (AryLen)                             ! Input array of x values.
   REAL(ReKi), INTENT(IN)       :: YAry  (AryLen)                             ! Input array of y values.

   INTEGER(IntKi), INTENT(OUT)  :: ErrStat                                    ! Error status.

   CHARACTER(*),   INTENT(OUT)  :: ErrMsg                                     ! Error message.


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
   SUBROUTINE RegCubicSplineInitM ( XAry, YAry, DelX, Coef, ErrStat, ErrMsg )


      ! This routine calculates the parameters needed to compute a regularly-spaced natural cubic spline.
      ! Natural cubic splines are used in that the curvature at the end points is zero.
      ! It assumes the XAry values are equally spaced for speed.
      ! This version of the routine works with multiple curves that share the same X values.


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
   FUNCTION RegCubicSplineInterp ( X, AryLen, XAry, YAry, DelX, Coef, ErrStat, ErrMsg )


      ! This routine interpolates a pair of arrays using cubic splines to find the function value at X.
      ! One must call RegCubicSplineInit() first to compute the coefficients of the cubics.
      ! This routine requires that the XAry be regularly spaced, which improves performance.


      ! Function declaration.

   REAL(ReKi)                   :: RegCubicSplineInterp                       ! This function.


      ! Argument declarations:

   INTEGER, INTENT(IN)          :: AryLen                                     ! Length of the array.

   REAL(ReKi), INTENT(IN)       :: Coef  (AryLen-1,0:3)                       ! The coefficients for the cubic polynomials.
   REAL(ReKi), INTENT(IN)       :: DelX                                       ! The distance between X values in XAry.
   REAL(ReKi), INTENT(IN)       :: X                                          ! The value we are trying to interpolate for.
   REAL(ReKi), INTENT(IN)       :: XAry (AryLen)                              ! Input array of regularly spaced x values.
   REAL(ReKi), INTENT(IN)       :: YAry (AryLen)                              ! Input array of y values.

   INTEGER(IntKi), INTENT(OUT)  :: ErrStat                                    ! Error status.

   CHARACTER(*),   INTENT(OUT)  :: ErrMsg                                     ! Error message.


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
   FUNCTION RegCubicSplineInterpM ( X, XAry, YAry, DelX, Coef, ErrStat, ErrMsg ) RESULT( Res )


      ! This routine interpolates a pair of arrays using cubic splines to find the function value at X.
      ! One must call RegCubicSplineInit() first to compute the coefficients of the cubics.
      ! This routine requires that the XAry be regularly spaced, which improves performance.
      ! This version of the routine works with multiple curves that share the same X values.


      ! Function declaration.

   REAL(ReKi), ALLOCATABLE      :: Res(:)                                     ! The result of this function.


      ! Argument declarations:

   REAL(ReKi), INTENT(IN)       :: Coef  (:,:,0:)                             ! The coefficients for the cubic polynomials.
   REAL(ReKi), INTENT(IN)       :: DelX                                       ! The distance between X values in XAry.
   REAL(ReKi), INTENT(IN)       :: X                                          ! The value we are trying to interpolate for.
   REAL(ReKi), INTENT(IN)       :: XAry (:)                                   ! Input array of regularly spaced x values.
   REAL(ReKi), INTENT(IN)       :: YAry (:,:)                                 ! Input array of y values.

   INTEGER(IntKi), INTENT(OUT)  :: ErrStat                                    ! Error status.

   CHARACTER(*),   INTENT(OUT)  :: ErrMsg                                     ! Error message.


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
   SUBROUTINE RombergInt(f, a, b, R, err, eps, ErrStat)

         ! This routine is used to integrate funciton f over the interval [a, b]. This routine
         ! is useful for sufficiently smooth (e.g., analytic) integrands, integrated over
         ! intervals which contain no singularities, and where the endpoints are also nonsingular.
         !
         ! f is an external function. For example f(x) = 1 + x.
         !
         !   FUNCTION f(x)
         !      USE PRECISION
         !      IMPLICIT NONE
         !
         !      REAL(ReKi) f
         !      REAL(ReKi) x
         !
         !      f = 1 + x
         !
         !      RETURN
         !   END FUNCTION f

      IMPLICIT NONE

         ! Argument declarations:

      REAL(ReKi), EXTERNAL              :: f               ! Integrand function name
      REAL(ReKi), INTENT(IN)            :: a               ! Lower integration limit
      REAL(ReKi), INTENT(IN)            :: b               ! Upper integration limit
      REAL(ReKi), INTENT(IN)            :: eps             ! Absolute error bound
      REAL(ReKi), INTENT(OUT)           :: R               ! The result of integration
      REAL(ReKi), INTENT(OUT)           :: err             ! Actual absolute error
      INTEGER, INTENT(OUT), OPTIONAL    :: ErrStat         ! Error status; if present, program does not abort on error

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
   SUBROUTINE RunTimes( StrtTime, UsrTime1, SimStrtTime, UsrTime2, ZTime, UsrTime_out )
   ! This routine displays a message that gives that status of the simulation and the predicted end time of day.

      IMPLICIT                        NONE

         ! Passed variables

      INTEGER   , INTENT(IN)       :: StrtTime (8)                                    ! Start time of simulation (including initialization)
      INTEGER   , INTENT(IN)       :: SimStrtTime (8)                                 ! Start time of simulation (after initialization)
      REAL      , INTENT(IN)       :: UsrTime1                                        ! User CPU time for simulation initialization.
      REAL,       INTENT(IN)       :: UsrTime2                                        ! User CPU time for simulation (without intialization)
      REAL(DbKi), INTENT(IN)       :: ZTime                                           ! The final simulation time (not necessarially TMax)
      REAL,OPTIONAL, INTENT(OUT)   :: UsrTime_out                                     ! User CPU time for entire run - optional value returned to calling routine

         ! Local variables

      REAL                         :: ClckTime                                        ! Elapsed clock time for the entire run.
      REAL                         :: ClckTimeSim                                     ! Elapsed clock time for the simulation phase of the run.
      REAL                         :: Factor                                          ! Ratio of seconds to a specified time period.
      REAL                         :: TRatio                                          ! Ratio of simulation time to elapsed clock time.
      REAL(ReKi), PARAMETER        :: SecPerDay = 24*60*60.0_ReKi                     ! Number of seconds per day

      REAL                         :: UsrTime                                         ! User CPU time for entire run.
      REAL                         :: UsrTimeSim                                      ! User CPU time for simulation (not including initialization).
      INTEGER                      :: EndTimes (8)                                    ! An array holding the ending clock time of the simulation.

      CHARACTER( 8)                :: TimePer
      CHARACTER(MaxWrScrLen)       :: BlankLine

         ! Get the end times to compare with start times.

      CALL DATE_AND_TIME ( VALUES=EndTimes )
      CALL CPU_TIME ( UsrTime )
      UsrTime = MAX( 0.0, UsrTime )  ! CPU_TIME: If a meaningful time cannot be returned, a processor-dependent negative value is returned
   

      ! Calculate the elapsed wall-clock time in seconds.

      ClckTime     = GetClockTime(StrtTime,      EndTimes)
     !ClckTimeInit = GetClockTime(StrtTime,   SimStrtTime)
      ClckTimeSim  = GetClockTime(SimStrtTime,   EndTimes)

         ! Calculate CPU times.

      UsrTime    = MAX( 0.0, UsrTime - UsrTime1 )
      UsrTimeSim = MAX( 0.0, UsrTime - UsrTime2 )


      IF ( .NOT. EqualRealNos( UsrTimeSim, 0.0 ) .AND. ZTime > 0.0_DbKi )  THEN

         TRatio = REAL(ZTime) / UsrTimeSim

         IF     ( UsrTime > SecPerDay )  THEN
            Factor = 1.0/SecPerDay
            TimePer = ' days'
         ELSEIF ( UsrTime >  3600.0 )  THEN
            Factor = 1.0/3600.0
            TimePer = ' hours'
         ELSEIF ( UsrTime >    60.0 )  THEN
            Factor = 1.0/60.0
            TimePer = ' minutes'
         ELSE
            Factor = 1.0
            TimePer = ' seconds'
         ENDIF

         BlankLine = ""
         CALL WrOver( BlankLine )  ! BlankLine contains MaxWrScrLen spaces
         CALL WrScr1( ' Total Real Time:       '//TRIM( Num2LStr( Factor*ClckTime      ) )//TRIM( TimePer ) )
         CALL WrScr ( ' Total CPU Time:        '//TRIM( Num2LStr( Factor*UsrTime       ) )//TRIM( TimePer ) )
   !     CALL WrScr ( ' ')
   !     CALL WrScr ( ' Simulation Real Time:  '//TRIM( Num2LStr( Factor*ClckTimeSim   ) )//TRIM( TimePer ) )
         CALL WrScr ( ' Simulation CPU Time:   '//TRIM( Num2LStr( Factor*UsrTimeSim    ) )//TRIM( TimePer ) )      
         CALL WrScr ( ' Simulated Time:        '//TRIM( Num2LStr( Factor*REAL( ZTime ) ) )//TRIM( TimePer ) )
         CALL WrScr ( ' Time Ratio (Sim/CPU):  '//TRIM( Num2LStr( TRatio ) ) )

      ENDIF

      IF (PRESENT(UsrTime_out)) UsrTime_out = UsrTime
      RETURN
   CONTAINS

      FUNCTION GetClockTime(StartClockTime, EndClockTime)
      ! return the number of seconds between StartClockTime and EndClockTime
   
         REAL                         :: GetClockTime          ! Elapsed clock time for the simulation phase of the run.
         INTEGER   , INTENT(IN)       :: StartClockTime (8)                                 ! Start time of simulation (after initialization)
         INTEGER   , INTENT(IN)       :: EndClockTime (8)                                 ! Start time of simulation (after initialization)
   
      !bjj: This calculation will be wrong at certain times (e.g. if it's near midnight on the last day of the month), but to my knowledge, no one has complained...
         GetClockTime =       0.001*( EndClockTime(8) - StartClockTime(8) ) &  ! Is the milliseconds of the second (range 0 to 999) - local time
                        +           ( EndClockTime(7) - StartClockTime(7) ) &  ! Is the seconds of the minute (range 0 to 59) - local time
                        +      60.0*( EndClockTime(6) - StartClockTime(6) ) &  ! Is the minutes of the hour (range 0 to 59) - local time
                        +    3600.0*( EndClockTime(5) - StartClockTime(5) ) &  ! Is the hour of the day (range 0 to 23) - local time
                        + SecPerDay*( EndClockTime(3) - StartClockTime(3) )    ! Is the day of the month
   
   
      END FUNCTION GetClockTime
   
   END SUBROUTINE RunTimes   
!=======================================================================   
   SUBROUTINE SetAnglesForInterp( angles )

      ! this routine takes angles (in radians) and converts them to appropriate
      ! ranges so they can be interpolated appropriately
      ! (i.e., interpolating between pi+.1 and -pi should give pi+0.5 
      ! instead of of 0.05 radians, so we give return the angles pi+.1 and -pi+2pi=pi
      ! we assume the interpolation occurs in the second dimension of angles
      ! and it is done for each angle in the first dimension
   
      REAL(ReKi), INTENT(INOUT)     :: angles(:,:)

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
   SUBROUTINE SetConstants( )

         ! This routine computes numeric constants stored in the NWTC Library

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


         ! IEEE constants:
      CALL Set_IEEE_Constants( NaN_D, Inf_D, NaN, Inf )
      

   RETURN
   END SUBROUTINE SetConstants
!=======================================================================   
   SUBROUTINE SimStatus_FirstTime( PrevSimTime, PrevClockTime, SimStrtTime, UsrTimeSim, ZTime, TMax )
      ! This routine displays a message that gives that status of the simulation.

      IMPLICIT                        NONE

         ! Passed variables
      REAL(DbKi), INTENT(IN   )    :: ZTime                                           ! Current simulation time (s)
      REAL(DbKi), INTENT(IN   )    :: TMax                                            ! Expected simulation time (s)
      REAL(DbKi), INTENT(  OUT)    :: PrevSimTime                                     ! Previous time message was written to screen (s > 0)
      REAL(ReKi), INTENT(  OUT)    :: PrevClockTime                                   ! Previous clock time in seconds past midnight
      INTEGER,    INTENT(  OUT)    :: SimStrtTime (8)                                 ! An array containing the elements of the start time.
      REAL,       INTENT(  OUT)    :: UsrTimeSim                                      ! User CPU time for simulation (without intialization)

         ! Local variables.

      REAL(ReKi)                   :: CurrClockTime                                   ! Current time in seconds past midnight.


         ! How many seconds past midnight?

      CALL DATE_AND_TIME ( Values=SimStrtTime )
      CALL CPU_TIME ( UsrTimeSim )                                                    ! Initial CPU time   
      UsrTimeSim = MAX( 0.0, UsrTimeSim )  ! CPU_TIME: If a meaningful time cannot be returned, a processor-dependent negative value is returned

      CurrClockTime = TimeValues2Seconds( SimStrtTime )


      CALL WrScr ( ' Timestep: '//TRIM( Num2LStr( NINT( ZTime ) ) )//' of '//TRIM( Num2LStr( TMax ) )//' seconds.')


      ! Let's save this time as the previous time for the next call to the routine
      PrevClockTime = CurrClockTime
      PrevSimTime   = ZTime

      RETURN
   END SUBROUTINE SimStatus_FirstTime
!=======================================================================
   SUBROUTINE SimStatus( PrevSimTime, PrevClockTime, ZTime, TMax )
   
   ! This routine displays a message that gives that status of the simulation and the predicted end time of day.

      IMPLICIT                        NONE

         ! Passed variables
      REAL(DbKi), INTENT(IN)       :: ZTime                                ! Current simulation time (s)
      REAL(DbKi), INTENT(IN)       :: TMax                                 ! Expected simulation time (s)
      REAL(DbKi), INTENT(INOUT)    :: PrevSimTime                          ! Previous time message was written to screen (s > 0)
      REAL(ReKi), INTENT(INOUT)    :: PrevClockTime                        ! Previous clock time in seconds past midnight


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

      CHARACTER( 8)                :: ETimeStr                             ! String containing the end time.


      IF ( ZTime <= PrevSimTime ) RETURN


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

      CALL WrOver ( ' Timestep: '//TRIM( Num2LStr( NINT( ZTime ) ) )//' of '//TRIM( Num2LStr( TMax ) )// &
                    ' seconds.  Estimated final completion at '//ETimeStr//'.'                             )

         ! Let's save this time as the previous time for the next call to the routine
      PrevClockTime = CurrClockTime
      PrevSimTime   = ZTime

      RETURN
   END SUBROUTINE SimStatus   
!=======================================================================
   SUBROUTINE SmllRotTrans( RotationType, Theta1, Theta2, Theta3, TransMat, ErrTxt, ErrStat, ErrMsg )


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

   REAL(ReKi), INTENT(IN )             :: Theta1                                          ! The small rotation about X1, (rad).
   REAL(ReKi), INTENT(IN )             :: Theta2                                          ! The small rotation about X2, (rad).
   REAL(ReKi), INTENT(IN )             :: Theta3                                          ! The small rotation about X3, (rad).
   REAL(ReKi), INTENT(OUT)             :: TransMat (3,3)                                  ! The resulting transformation matrix from X to x, (-).

   INTEGER(IntKi),INTENT(OUT)          :: ErrStat
   CHARACTER(*), INTENT(OUT)           :: ErrMsg

   CHARACTER(*), INTENT(IN)            :: RotationType                                    ! The type of rotation; used to inform the user where a large rotation is occuring upon such an event.
   CHARACTER(*), INTENT(IN ), OPTIONAL :: ErrTxt                                          ! an additional message to be displayed as a warning (typically the simulation time)

      ! Local Variables:

   REAL(ReKi)                          :: ComDenom                                        ! = ( Theta1^2 + Theta2^2 + Theta3^2 )*SQRT( 1.0 + Theta1^2 + Theta2^2 + Theta3^2 )
   REAL(ReKi), PARAMETER               :: LrgAngle  = 0.4                                 ! Threshold for when a small angle becomes large (about 23deg).  This comes from: COS(SmllAngle) ~ 1/SQRT( 1 + SmllAngle^2 ) and SIN(SmllAngle) ~ SmllAngle/SQRT( 1 + SmllAngle^2 ) results in ~5% error when SmllAngle = 0.4rad.
   REAL(ReKi)                          :: Theta11                                         ! = Theta1^2
   REAL(ReKi)                          :: Theta12S                                        ! = Theta1*Theta2*[ SQRT( 1.0 + Theta1^2 + Theta2^2 + Theta3^2 ) - 1.0 ]
   REAL(ReKi)                          :: Theta13S                                        ! = Theta1*Theta3*[ SQRT( 1.0 + Theta1^2 + Theta2^2 + Theta3^2 ) - 1.0 ]
   REAL(ReKi)                          :: Theta22                                         ! = Theta2^2
   REAL(ReKi)                          :: Theta23S                                        ! = Theta2*Theta3*[ SQRT( 1.0 + Theta1^2 + Theta2^2 + Theta3^2 ) - 1.0 ]
   REAL(ReKi)                          :: Theta33                                         ! = Theta3^2
   REAL(ReKi)                          :: SqrdSum                                         ! = Theta1^2 + Theta2^2 + Theta3^2
   REAL(ReKi)                          :: SQRT1SqrdSum                                    ! = SQRT( 1.0 + Theta1^2 + Theta2^2 + Theta3^2 )

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
   SQRT1SqrdSum = SQRT( 1.0 + SqrdSum )
   ComDenom     = SqrdSum*SQRT1SqrdSum

   Theta12S     = Theta1*Theta2*( SQRT1SqrdSum - 1.0 )
   Theta13S     = Theta1*Theta3*( SQRT1SqrdSum - 1.0 )
   Theta23S     = Theta2*Theta3*( SQRT1SqrdSum - 1.0 )


      ! Define the transformation matrix:

   IF ( ComDenom == 0.0 )  THEN  ! All angles are zero and matrix is ill-conditioned (the matrix is derived assuming that the angles are not zero); return identity

      TransMat(1,:) = (/ 1.0, 0.0, 0.0 /)
      TransMat(2,:) = (/ 0.0, 1.0, 0.0 /)
      TransMat(3,:) = (/ 0.0, 0.0, 1.0 /)

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
   END SUBROUTINE SmllRotTrans
!=======================================================================
   SUBROUTINE SortUnion ( Ary1, N1, Ary2, N2, Ary, N )


      ! This routine takes two sorted arrays and finds the sorted union of the two.

      ! Note: If the same value is found in both arrays, only one is kept.  However, if either
      !       array as multiple occurances of the same value, the largest multiple will be
      !       kept.  Duplicates should be eliminated externally if this is not desirable.


      ! Argument declarations:

   INTEGER, INTENT(OUT)         :: N                                            ! The length of the output array.
   INTEGER, INTENT(IN)          :: N1                                           ! The length of the first input array.
   INTEGER, INTENT(IN)          :: N2                                           ! The length of the second input array.

   REAL(ReKi), INTENT(OUT)      :: Ary(N1+N2)                                   ! The sorted union.
   REAL(ReKi), INTENT(IN)       :: Ary1(N1)                                     ! The first list of sorted real numbers.
   REAL(ReKi), INTENT(IN)       :: Ary2(N2)                                     ! The second list of sorted real numbers.


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
   FUNCTION StdDevFn ( Ary, AryLen, Mean )


      ! This routine calculates the standard deviation of a population contained in Ary.


      ! Function declaration.

   REAL(ReKi)                   :: StdDevFn                                     ! This function.


      ! Argument declarations:

   INTEGER, INTENT(IN)          :: AryLen                                       ! Length of the array.

   REAL(ReKi), INTENT(IN)       :: Ary  (AryLen)                                ! Input array.
   REAL(ReKi), INTENT(IN)       :: Mean                                         ! The previously calculated mean of the array.


      ! Local declarations.

   REAL(DbKi)                   :: Sum                                          ! A temporary sum.

   INTEGER                      :: I                                            ! The index into the array.



   Sum = 0.0_DbKi

   DO I=1,AryLen
      Sum = Sum + ( Ary(I) - Mean )**2
   END DO ! I

   StdDevFn = SQRT( Sum/( AryLen - 1 ) )


   RETURN
   END FUNCTION StdDevFn ! ( Ary, AryLen, Mean )
!=======================================================================
   FUNCTION TimeValues2Seconds( TimeAry )
      ! This routine takes an array of time values such as that returned from
      !     CALL DATE_AND_TIME ( Values=TimeAry )
      ! and converts TimeAry to the number of seconds past midnight.

      ! Passed variables:
   INTEGER, INTENT(IN)          :: TimeAry  (8)                                    ! An array containing the elements of the time
   REAL(ReKi)                   :: TimeValues2Seconds                              ! Current time in seconds past midnight


   TimeValues2Seconds = 3600*TimeAry(5) + 60*TimeAry(6) + TimeAry(7) + 0.001_ReKi*TimeAry(8)

   END FUNCTION TimeValues2Seconds
!=======================================================================
   FUNCTION trace(A)
   
      ! This function computes the trace of a square matrix:
      ! SUM ( A(i,i) ) for i=1, min( SIZE(A,1), SIZE(A,2) )
      
   REAL(ReKi), INTENT(IN)  :: A(:,:)
   REAL(ReKi)              :: trace
   
   INTEGER(IntKi)          :: n     ! rows/cols in A
   INTEGER(IntKi)          :: i     ! loop counter
   
   n = min( SIZE(A,1), SIZE(A,2) )

   trace = 0.0_ReKi
   do i=1,n
      trace = trace + A(i,i)
   end do
   
   END FUNCTION trace
!=======================================================================
   FUNCTION TwoNorm(v)
   
      ! this function returns the 2-norm of a vector v
      ! fortran 2008 has Norm2() built in
      
      REAL(ReKi), INTENT(IN)  :: v(:)      
      REAL(ReKi)              :: TwoNorm      
      
      TwoNorm = SQRT( DOT_PRODUCT(v, v) )
      
      
   END FUNCTION
!=======================================================================  
   SUBROUTINE Zero2TwoPi ( Angle )

      ! This routine is used to convert Angle to an equivalent value
      !  in the range [0, 2*pi).
      

      ! Argument declarations:

   REAL(ReKi), INTENT(INOUT)    :: Angle



      ! Get the angle between 0 and 2Pi.

   Angle = MODULO( Angle, TwoPi )   


      ! Check numerical case where Angle == 2Pi.

   IF ( Angle == TwoPi )  THEN
      Angle = 0.0_ReKi
   END IF


   RETURN
   END SUBROUTINE Zero2TwoPi   
!=======================================================================  
END MODULE NWTC_Num
