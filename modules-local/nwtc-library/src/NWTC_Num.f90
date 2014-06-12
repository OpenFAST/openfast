!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2013  National Renewable Energy Laboratory
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
   !  FUNCTION   IsSymmetric           ( A )                                                          ! Function to determine if A(:,:) is symmetric
   !  SUBROUTINE LocateBin             ( XVal, XAry, Ind, AryLen )
   !  SUBROUTINE LocateStp             ( XVal, XAry, Ind, AryLen )
   !  FUNCTION   Mean                  ( Ary, AryLen )                                                ! Function to calculate the mean value of a vector array.
   !  SUBROUTINE MPi2Pi                ( Angle )
   !  SUBROUTINE Zero2TwoPi            ( Angle )   
   !  SUBROUTINE RegCubicSplineInit    ( AryLen, XAry, YAry, DelX, Coef )                             ! Calculate coefficients for regularly spaced array to use cubic splines.
   !  SUBROUTINE RegCubicSplineInitM   ( XAry, YAry, DelX, Coef, ErrStat, ErrMsg )                    ! Interpolate using cubic splines for multiple tables of regularly space data.
   !  FUNCTION   RegCubicSplineInterp  ( X, AryLen, XAry, YAry, DelX, Coef )                          ! Interpolate a regularly spaced array using cubic splines.
   !  FUNCTION   RegCubicSplineInterpM ( X, XAry, YAry, DelX, Coef, ErrStat, ErrMsg )                 ! Initialize cubic splines for multiple tables of regularly space data.
   !  SUBROUTINE RombergInt            ( f, a, b, R, err, eps, ErrStat )
   !  SUBROUTINE SetConstants
   !  SUBROUTINE SmllRotTrans          ( RotationType, Theta1, Theta2, Theta3, TransMat, ErrTxt )
   !  SUBROUTINE SortUnion             ( Ary1, N1, Ary2, N2, Ary, N )
   !  FUNCTION   StdDevFn              ( Ary, AryLen, Mean )                                          ! Function to calculate the standard deviation of a vector array.
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


   TYPE, PUBLIC               :: CubSplineType                                ! This derived type is used to hold data for performing cubic splines.
      INTEGER                                :: NumPts                        ! The number of points in the XAry and YAry arrays.
      REAL(ReKi), ALLOCATABLE                :: Coef      (:,:)               ! The NumPts-1 length array of cubic coefficients.  The second dimension must be "0:3".
      REAL(ReKi), ALLOCATABLE                :: XAry      (:)                 ! The NumPts length array of x values for the interpolation.
      REAL(ReKi), ALLOCATABLE                :: YAry      (:)                 ! The NumPts length array of y values for the interpolation.
   END TYPE CubSplineType

   TYPE, PUBLIC               :: RegCubSplineType                             ! This derived type is used to hold data for performing cubic splines wuth regularly-spaced data.
      INTEGER                                :: NumPts                        ! The number of points in the XAry and YAry arrays.
      REAL(ReKi), ALLOCATABLE                :: Coef      (:,:)               ! The NumPts-1 length array of cubic coefficients.  The second dimension must be "0:3".
      REAL(ReKi)                             :: DelX                          ! The distance between the equally spaced points in XAry.
      REAL(ReKi), ALLOCATABLE                :: XAry      (:)                 ! The NumPts length array of x values for the interpolation.
      REAL(ReKi), ALLOCATABLE                :: YAry      (:)                 ! The NumPts length array of y values for the interpolation.
   END TYPE RegCubSplineType

   TYPE, PUBLIC               :: RegGridType                                  ! This derived type is used to hold the contents of a regular grid of data.
      INTEGER                                :: NumDims                       ! The number of dimensions for this grid.
      REAL(ReKi), ALLOCATABLE                :: Mins      (:)                 ! The set of minimums for the grid in each NumDims dimensions.
      REAL(ReKi), ALLOCATABLE                :: Steps     (:)                 ! The set of step sizes for the grid in each NumDims dimensions.
      REAL(ReKi), ALLOCATABLE                :: Grid      (:)                 ! The NumDims dimensional grid.
   END TYPE RegGridType


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

   CHARACTER(4096), INTENT(OUT) :: ErrMsg                                     ! Error message.


      ! Local declarations.

   REAL(ReKi), ALLOCATABLE      :: DelX  (:)                                  ! The distances between the randomly spaced points.
   REAL(ReKi), ALLOCATABLE      :: Slope (:)                                  ! The AryLen-1 length array of slopes between points.
   REAL(ReKi), ALLOCATABLE      :: U     (:)                                  ! An AryLen-1 length array used in the Gaussian elimination.
   REAL(ReKi), ALLOCATABLE      :: V     (:)                                  ! An AryLen-1 length array used in the Gaussian elimination.
   REAL(ReKi)                   :: ZHi                                        ! A parameter used to calculate the polynomial coefficients.
   REAL(ReKi)                   :: ZLo                                        ! A parameter used to calculate the polynomial coefficients.

   INTEGER(IntKi)               :: ErrStatLcL                                 ! Local error status.
   INTEGER                      :: I                                          ! The index into the arrays.



      ! Allocate the various intermediate arrays.

   ALLOCATE ( DelX( AryLen - 1 ), STAT=ErrStatLcL )
   IF ( ErrStatLcL /= 0 )  THEN
      CALL ExitThisRoutine ( ErrID_Fatal, NewLine//' >> Error allocating memory for the DelX array in CubicSplineInit.' )
      RETURN
   ENDIF

   ALLOCATE ( Slope( AryLen - 1 ), STAT=ErrStatLcL )
   IF ( ErrStatLcL /= 0 )  THEN
      CALL ExitThisRoutine ( ErrID_Fatal, NewLine//' >> Error allocating memory for the Slope array in CubicSplineInit.' )
      RETURN
   ENDIF

   ALLOCATE ( U( AryLen - 1 ), STAT=ErrStatLcL )
   IF ( ErrStatLcL /= 0 )  THEN
      CALL ExitThisRoutine ( ErrID_Fatal, NewLine//' >> Error allocating memory for the U array in CubicSplineInit.' )
      RETURN
   ENDIF

   ALLOCATE ( V( AryLen - 1 ), STAT=ErrStatLcL )
   IF ( ErrStatLcL /= 0 )  THEN
      CALL ExitThisRoutine ( ErrID_Fatal, NewLine//' >> Error allocating memory for the V array in CubicSplineInit.' )
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

   CHARACTER(4096), INTENT(OUT) :: ErrMsg                                     ! Error message.


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



      ! How big are the arrays?

   NumPts  = SIZE( XAry )
   NumCrvs = SIZE( YAry, 2 )


      ! Allocate the various intermediate arrays.

   ALLOCATE ( ZLo( NumCrvs ), STAT=ErrStatLcL )
   IF ( ErrStatLcL /= 0 )  THEN
      CALL ExitThisRoutine ( ErrID_Fatal, NewLine//' >> Error allocating memory for the ZLo array in CubicSplineInitM.' )
      RETURN
   ENDIF

   ALLOCATE ( ZHi( NumCrvs ), STAT=ErrStatLcL )
   IF ( ErrStatLcL /= 0 )  THEN
      CALL ExitThisRoutine ( ErrID_Fatal, NewLine//' >> Error allocating memory for the ZHi array in CubicSplineInitM.' )
      RETURN
   ENDIF

   ALLOCATE ( DelX( NumPts - 1 ), STAT=ErrStatLcL )
   IF ( ErrStatLcL /= 0 )  THEN
      CALL ExitThisRoutine ( ErrID_Fatal, NewLine//' >> Error allocating memory for the DelX array in CubicSplineInitM.' )
      RETURN
   ENDIF

   ALLOCATE ( Slope( NumPts-1, NumCrvs ), STAT=ErrStatLcL )
   IF ( ErrStatLcL /= 0 )  THEN
      CALL ExitThisRoutine ( ErrID_Fatal, NewLine//' >> Error allocating memory for the Slope array in CubicSplineInitM.' )
      RETURN
   ENDIF

   ALLOCATE ( U( NumPts - 1 ), STAT=ErrStatLcL )
   IF ( ErrStatLcL /= 0 )  THEN
      CALL ExitThisRoutine ( ErrID_Fatal, NewLine//' >> Error allocating memory for the U array in CubicSplineInitM.' )
      RETURN
   ENDIF

   ALLOCATE ( V( NumPts-1, NumCrvs ), STAT=ErrStatLcL )
   IF ( ErrStatLcL /= 0 )  THEN
      CALL ExitThisRoutine ( ErrID_Fatal, NewLine//' >> Error allocating memory for the V array in CubicSplineInitM.' )
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

   CHARACTER(4096), INTENT(OUT) :: ErrMsg                                     ! Error message.


      ! Local declarations.

   REAL(ReKi)                   :: XOff                                       ! The distance from X to XAry(ILo).

!   INTEGER(IntKi)               :: ErrStatLcL                                 ! Local error status.
   INTEGER                      :: ILo                                        ! The index into the array for which X is just above or equal to XAry(ILo).



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


         RETURN

      END SUBROUTINE ExitThisRoutine ! ( ErrID, Msg )

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

   CHARACTER(4096), INTENT(OUT) :: ErrMsg                                     ! Error message.


      ! Local declarations.

   REAL(ReKi)                   :: XOff                                       ! The distance from X to XAry(ILo).

   INTEGER(IntKi)               :: ErrStatLcL                                 ! Local error status.
   INTEGER                      :: ILo                                        ! The index into the array for which X is just above or equal to XAry(ILo).
   INTEGER                      :: NumCrvs                                    ! Number of curves to be interpolated.
   INTEGER                      :: NumPts                                     ! Number of points in each curve.



      ! How big are the arrays?

   NumPts  = SIZE( XAry )
   NumCrvs = SIZE( YAry, 2 )
   NumCrvs = SIZE( YAry, 2 )

   ALLOCATE ( Res( NumCrvs ) , STAT=ErrStatLcl )
   IF ( ErrStatLcl /= 0 )  THEN
      CALL ExitThisRoutine ( ErrID_Fatal, '  >> Error allocating memory for the function result array in RegCubicSplineInterpM.' )
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


      ! We are somewhere inside XAry.  Find the segment that bounds X using binary search.

   CALL LocateBin( X, XAry, ILo, NumPts )

   XOff = X - XAry(ILo)

   Res(:) = Coef(ILo,:,0) + XOff*( Coef(ILo,:,1) + XOff*( Coef(ILo,:,2) + XOff*Coef(ILo,:,3) ) )


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


         RETURN

      END SUBROUTINE ExitThisRoutine ! ( ErrID, Msg )

   END FUNCTION CubicSplineInterpM ! ( X, XAry, YAry, Coef, ErrStat, ErrMsg )
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

   CHARACTER(4096), INTENT(OUT) :: ErrMsg                                     ! Error message.


      ! Local declarations.

   REAL(ReKi)                   :: Del                                        ! The distances between the regularly spaced points.

!   INTEGER(IntKi)               :: ErrStatLcL                                 ! Local error status.



      ! Check the validity of the data.

   IF ( NumPts == 0 )  THEN
      CALL ExitThisRoutine ( ErrID_Fatal, ' >> The value of NumPts cannot be zero when calling GetOffsetReg.' )
      RETURN
   END IF

   IF ( NumPts == 0 )  THEN
      CALL ExitThisRoutine ( ErrID_Fatal, ' >> The value of NumPts cannot be zero when calling GetOffsetReg.' )
      RETURN
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


         RETURN

      END SUBROUTINE ExitThisRoutine ! ( ErrID, Msg )

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
         ErrStat = 1
         RETURN
      END IF
   END IF

   IF ( ( Ipt < 1 ) .OR. ( Ipt > NPts ) )  THEN
      CALL ProgAbort ( ' In function GL_Loc, the point being used for Gauss-Legendre Quadrature must be between 1 and ' &
                   //TRIM( Int2LStr( NPts ) )//' (inclusive).  Instead, it is "'//TRIM( Int2LStr( Ipt ) )//'".', PRESENT(ErrStat) )
      IF ( PRESENT(ErrStat) ) THEN
         ErrStat = 1
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

   CHARACTER(4096), INTENT(OUT) :: ErrMsg                                     ! Error message.


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



      ! Allocate the various intermediate arrays.

   ALLOCATE ( Slope( AryLen - 1 ), STAT=ErrStatLcL )
   IF ( ErrStatLcL /= 0 )  THEN
      CALL ExitThisRoutine ( ErrID_Fatal, NewLine//' >> Error allocating memory for the Slope array in RegCubicSplineInit.' )
      RETURN
   ENDIF

   ALLOCATE ( U( AryLen - 1 ), STAT=ErrStatLcL )
   IF ( ErrStatLcL /= 0 )  THEN
      CALL ExitThisRoutine ( ErrID_Fatal, NewLine//' >> Error allocating memory for the U array in RegCubicSplineInit.' )
      RETURN
   ENDIF

   ALLOCATE ( V( AryLen - 1 ), STAT=ErrStatLcL )
   IF ( ErrStatLcL /= 0 )  THEN
      CALL ExitThisRoutine ( ErrID_Fatal, NewLine//' >> Error allocating memory for the V array in RegCubicSplineInit.' )
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

   CHARACTER(4096), INTENT(OUT) :: ErrMsg                                     ! Error message.


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



      ! How big are the arrays?

   NumPts  = SIZE( XAry )
   NumCrvs = SIZE( YAry, 2 )


      ! Allocate the various intermediate arrays.

   ALLOCATE ( ZLo( NumCrvs ), STAT=ErrStatLcL )
   IF ( ErrStatLcL /= 0 )  THEN
      CALL ExitThisRoutine ( ErrID_Fatal, NewLine//' >> Error allocating memory for the ZLo array in CubicSplineInitM.' )
      RETURN
   ENDIF

   ALLOCATE ( ZHi( NumCrvs ), STAT=ErrStatLcL )
   IF ( ErrStatLcL /= 0 )  THEN
      CALL ExitThisRoutine ( ErrID_Fatal, NewLine//' >> Error allocating memory for the ZHi array in CubicSplineInitM.' )
      RETURN
   ENDIF

   ALLOCATE ( Slope( NumPts-1, NumCrvs ), STAT=ErrStatLcL )
   IF ( ErrStatLcL /= 0 )  THEN
      CALL ExitThisRoutine ( ErrID_Fatal, NewLine//' >> Error allocating memory for the Slope array in RegCubicSplineInitM.' )
      RETURN
   ENDIF

   ALLOCATE ( U( NumPts - 1 ), STAT=ErrStatLcL )
   IF ( ErrStatLcL /= 0 )  THEN
      CALL ExitThisRoutine ( ErrID_Fatal, NewLine//' >> Error allocating memory for the U array in RegCubicSplineInitM.' )
      RETURN
   ENDIF

   ALLOCATE ( V( NumPts-1, NumCrvs ), STAT=ErrStatLcL )
   IF ( ErrStatLcL /= 0 )  THEN
      CALL ExitThisRoutine ( ErrID_Fatal, NewLine//' >> Error allocating memory for the V array in RegCubicSplineInitM.' )
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

   CHARACTER(4096), INTENT(OUT) :: ErrMsg                                     ! Error message.


      ! Local declarations.

   REAL(ReKi)                   :: XOff                                       ! The distance from X to XAry(ILo).

   INTEGER                      :: ILo                                        ! The index into the array for which X is just above or equal to XAry(ILo).



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


         RETURN

      END SUBROUTINE ExitThisRoutine ! ( ErrID, Msg )

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

   CHARACTER(4096), INTENT(OUT) :: ErrMsg                                     ! Error message.


      ! Local declarations.

   REAL(ReKi)                   :: XOff                                       ! The distance from X to XAry(ILo).

   INTEGER                      :: ErrStatLcL                                 ! Local error status.
   INTEGER                      :: ILo                                        ! The index into the array for which X is just above or equal to XAry(ILo).
   INTEGER                      :: NumCrvs                                    ! Number of curves.
   INTEGER                      :: NumPts                                     ! Number of points in each curve.



      ! How big are the arrays?  Use the size to allocate the result.

   NumPts  = SIZE( XAry )
   NumCrvs = SIZE( YAry, 2 )

   ALLOCATE ( Res( NumCrvs ) , STAT=ErrStatLcl )
   IF ( ErrStatLcl /= 0 )  THEN
      CALL ExitThisRoutine ( ErrID_Fatal, '  >> Error allocating memory for the function result array in RegCubicSplineInterpM.' )
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


         RETURN

      END SUBROUTINE ExitThisRoutine ! ( ErrID, Msg )

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

      INTEGER                           :: m, i, j, k
      INTEGER, PARAMETER                :: mmax = 50       ! Maximum iteration number for m
      INTEGER, PARAMETER                :: imax = 50       ! Maximum iteration number for i

      REAL(ReKi), ALLOCATABLE           :: T(:,:)
      REAL(ReKi)                        :: h               ! Step length
      REAL(ReKi)                        :: sumf

         ! Initialize T
      ALLOCATE( T( mmax, imax ) )
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
         ErrStat = 1
         RETURN
      END IF

      RETURN
   END SUBROUTINE RombergInt
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
   SUBROUTINE Zero2TwoPi ( Angle )

      ! This routine is used to convert Angle to an equivalent value
      !  between 0 and 2*pi.
      

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
