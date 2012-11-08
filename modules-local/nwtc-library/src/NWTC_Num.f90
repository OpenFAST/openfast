MODULE NWTC_Num


   ! This module contains numeric-type routines with non-system-specific logic and references.


   ! It contains the following routines:

   !     SUBROUTINE AddOrSub2Pi   ( OldAngle, NewAngle )
   !     SUBROUTINE BSortReal     ( RealAry, NumPts )
   !     FUNCTION   CROSS_PRODUCT ( Vector1, Vector2 )
   !     FUNCTION   EqualRealNos  ( ReNum1, ReNum2 )
   !     SUBROUTINE GL_Pts        ( IPt, NPts, Loc, Wt [, ErrStat] )
   !     FUNCTION   IndexCharAry  ( CVal, CAry )
   !     FUNCTION   InterpBin     ( XVal, XAry, YAry, ILo, AryLen )             ! Generic interface for InterpBinComp and InterpBinReal.
   !     FUNCTION   InterpBinComp ( XVal, XAry, YAry, ILo, AryLen )
   !     FUNCTION   InterpBinReal ( XVal, XAry, YAry, ILo, AryLen )
   !     FUNCTION   InterpStp     ( XVal, XAry, YAry, ILo, AryLen )             ! Generic interface for InterpStpComp and InterpStpReal.
   !     FUNCTION   InterpStpComp ( XVal, XAry, YAry, Ind, AryLen )
   !     FUNCTION   InterpStpReal ( XVal, XAry, YAry, Ind, AryLen )
   !     SUBROUTINE LocateStp     ( XVal, XAry, Ind, AryLen )
   !     FUNCTION   Mean          ( Ary, AryLen )                               ! Function to calculate the mean value of a vector array.
   !     SUBROUTINE MPi2Pi        ( Angle )
   !     SUBROUTINE RombergInt    ( f, a, b, R, err, eps, ErrStat )
   !     SUBROUTINE SetConstants
   !     SUBROUTINE SmllRotTrans  ( RotationType, Theta1, Theta2, Theta3, TransMat, ErrTxt )
   !     SUBROUTINE SortUnion     ( Ary1, N1, Ary2, N2, Ary, N )
   !     FUNCTION   StdDevFn      ( Ary, AryLen, Mean )                         ! Function to calculate the standard deviation of a vector array.


   USE                             NWTC_IO

   IMPLICIT NONE

!=======================================================================


      ! Global numeric-related variables.

   REAL(DbKi)                   :: D2R_D                                        ! Factor to convert degrees to radians in double precision.
   REAL(DbKi)                   :: Inf_D                                        ! IEEE value for NaN (not-a-number) in double precision 
   REAL(DbKi)                   :: NaN_D                                        ! IEEE value for Inf (infinity) in double precision
   REAL(DbKi)                   :: Pi_D                                         ! Ratio of a circle's circumference to its diameter in double precision.
   REAL(DbKi)                   :: PiBy2_D                                      ! Pi/2 in double precision.
   REAL(DbKi)                   :: R2D_D                                        ! Factor to convert radians to degrees in double precision.
   REAL(DbKi)                   :: RPM2RPS_D                                    ! Factor to convert revolutions per minute to radians per second in double precision.
   REAL(DbKi)                   :: RPS2RPM_D                                    ! Factor to convert radians per second to revolutions per minute in double precision.
   REAL(DbKi)                   :: TwoByPi_D                                    ! 2/Pi in double precision.
   REAL(DbKi)                   :: TwoPi_D                                      ! 2*Pi in double precision.
   
   
   REAL(ReKi)                   :: D2R                                          ! Factor to convert degrees to radians.
   REAL(ReKi)                   :: Inf                                          ! IEEE value for NaN (not-a-number)
   REAL(ReKi)                   :: NaN                                          ! IEEE value for Inf (infinity)
   REAL(ReKi)                   :: Pi                                           ! Ratio of a circle's circumference to its diameter.
   REAL(ReKi)                   :: PiBy2                                        ! Pi/2.
   REAL(ReKi)                   :: R2D                                          ! Factor to convert radians to degrees.
   REAL(ReKi)                   :: RPM2RPS                                      ! Factor to convert revolutions per minute to radians per second.
   REAL(ReKi)                   :: RPS2RPM                                      ! Factor to convert radians per second to revolutions per minute.
   REAL(ReKi)                   :: TwoByPi                                      ! 2/Pi.
   REAL(ReKi)                   :: TwoPi                                        ! 2*Pi.

   INTEGER, ALLOCATABLE         :: IntIndx  (:,:)                               ! The array of indices holding that last index used for interpolation in IntBlade().


!=======================================================================

      ! Create interface for a generic EqualRealNos that uses specific routines.

   INTERFACE EqualRealNos
      MODULE PROCEDURE EqualRealNos4
      MODULE PROCEDURE EqualRealNos8
      MODULE PROCEDURE EqualRealNos16
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


      ! This routine sorts a list of real numbers.  It uses the buble sort algorithm,
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
   FUNCTION GetSmllRotAngs ( DCMat, ErrStat )

      ! This subroutine computes the angles that make up the input
      ! direction cosine matrix, DCMat

      ! passed variables

   REAL(ReKi), INTENT(IN )         :: DCMat          (3,3)
   INTEGER,    INTENT(OUT )        :: ErrStat               ! a non-zero value indicates an error in the permutation matrix algorithm

   REAL(ReKi)                      :: GetSmllRotAngs ( 3 )

      ! local variables
   REAL(ReKi)                      :: denom                 ! the denominator of the resulting matrix
   REAL(ReKi), PARAMETER           :: LrgAngle  = 0.4       ! Threshold for when a small angle becomes large (about 23deg).  This comes from: COS(SmllAngle) ~ 1/SQRT( 1 + SmllAngle^2 ) and SIN(SmllAngle) ~ SmllAngle/SQRT( 1 + SmllAngle^2 ) results in ~5% error when SmllAngle = 0.4rad.
   REAL(ReKi), PARAMETER           :: TOL = EPSILON(TOL)    ! tolerance for division by zero



      ! initialize output angles (just in case there is an error that prevents them from getting set)

   GetSmllRotAngs = 0.0
   ErrStat        = 0


      ! calculate the small angles
   GetSmllRotAngs(1) = DCMat(2,3) - DCMat(3,2)
   GetSmllRotAngs(2) = DCMat(3,1) - DCMat(1,3)
   GetSmllRotAngs(3) = DCMat(1,2) - DCMat(2,1)

   denom             = DCMat(1,1) + DCMat(2,2) + DCMat(3,3) - 1

   IF ( ABS(denom) > TOL ) THEN
      GetSmllRotAngs = GetSmllRotAngs / denom

               ! check that the angles are, in fact, small
      IF ( ANY( ABS(GetSmllRotAngs) > LrgAngle ) ) THEN
         CALL ProgWarn( ' Angles in GetSmllRotAngs() are larger than '//TRIM(Num2LStr(LrgAngle))//' radians.' )
         ErrStat = 1
      END IF

   ELSE
            ! check that the angles are, in fact, small (denom should be close to 2 if angles are small)
      CALL ProgAbort( ' Denominator is zero in GetSmllRotAngs().', TrapErrors = .TRUE. )
      ErrStat = -1

   END IF


   END FUNCTION GetSmllRotAngs ! ( DCMat, PMat, ErrStat )
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
      Ind           = AryLen - 1
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
      Ind           = AryLen - 1
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


      ! Function declaration.

   REAL(ReKi)                   :: Mean                                         ! This function.


      ! Argument declarations:

   INTEGER, INTENT(IN)          :: AryLen                                       ! Length of the array.

   REAL(ReKi), INTENT(IN)       :: Ary  (AryLen)                                ! Input array.


      ! Local declarations.

   INTEGER                      :: I                                            ! The index into the array.



   Mean = 0.0

   DO I=1,AryLen
      Mean = Mean + Ary(I)
   END DO ! I

   Mean = Mean/AryLen


   RETURN
   END FUNCTION Mean ! ( Ary, AryLen )
!=======================================================================
   SUBROUTINE MPi2Pi ( Angle )


      ! This routine ensures that Angle lies between -pi and pi.


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

!   USE, INTRINSIC :: ieee_arithmetic  !use this for compilers that have implemented 

      ! local variables for getting values of NaN and Inf (not necessary when using ieee_arithmetic)
   REAL(DbKi)                          :: Neg_D          ! a negative real(DbKi) number
   REAL(ReKi)                          :: Neg            ! a negative real(ReKi) number                                        

      
      ! Constants based upon Pi:

   Pi_D      = ACOS( -1.0_DbKi )
   D2R_D     = Pi_D/180.0_DbKi
   R2D_D     = 180.0_DbKi/Pi_D
   PiBy2_D   = Pi_D/2.0_DbKi
   RPM2RPS_D = Pi_D/30.0_DbKi
   RPS2RPM_D = 30.0_DbKi/Pi_D
   TwoByPi_D = 2.0_DbKi/Pi_D
   TwoPi_D   = 2.0_DbKi*Pi_D

   Pi      = ACOS( -1.0_ReKi )
   D2R     = Pi/180.0_ReKi
   R2D     = 180.0_ReKi/Pi
   PiBy2   = Pi/2.0_ReKi
   RPM2RPS = Pi/30.0_ReKi
   RPS2RPM = 30.0_ReKi/Pi
   TwoByPi =  2.0_ReKi/Pi
   TwoPi   =  2.0_ReKi*Pi
   
   
      ! IEEE constants:
      
!   NaN_D = ieee_value(0.0_DbKi, ieee_quiet_nan)
!   Inf_D = ieee_value(0.0_DbKi, ieee_positive_inf)
!
!   NaN   = ieee_value(0.0_ReKi, ieee_quiet_nan)
!   Inf   = ieee_value(0.0_DbKi, ieee_positive_inf)
   
      ! set variables to negative numbers to calculate NaNs (compilers may complain when taking sqrt of negative constants)
   Neg   = -1.0_ReKi
   Neg_D = -1.0_DbKi
   
   NaN_D = SQRT ( Neg_D )
   Inf_D = Pi_D / 0.0_DbKi

   NaN   = SQRT ( Neg )
   Inf   = Pi / 0.0_ReKi
   
   
   RETURN
   END SUBROUTINE SetConstants
!=======================================================================
   SUBROUTINE SmllRotTrans( RotationType, Theta1, Theta2, Theta3, TransMat, ErrTxt )


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

      ! The algebraic form of the transformation matrix, as implemented
      !   below, was derived symbolically by J. Jonkman by computing U*V^T
      !   by hand with verification in Mathematica.



      ! Passed Variables:

   REAL(ReKi), INTENT(IN )             :: Theta1                                          ! The small rotation about X1, (rad).
   REAL(ReKi), INTENT(IN )             :: Theta2                                          ! The small rotation about X2, (rad).
   REAL(ReKi), INTENT(IN )             :: Theta3                                          ! The small rotation about X3, (rad).
   REAL(ReKi), INTENT(OUT)             :: TransMat (3,3)                                  ! The resulting transformation matrix from X to x, (-).

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



      ! Display a warning message if at least one angle gets too large in
      !   magnitude:

   IF ( ( ( ABS(Theta1) > LrgAngle ) .OR. ( ABS(Theta2) > LrgAngle ) .OR. ( ABS(Theta3) > LrgAngle ) ) .AND. FrstWarn )  THEN

               
      CALL ProgWarn(' Small angle assumption violated in SUBROUTINE SmllRotTrans() due to'// &
                     ' a large '//TRIM(RotationType)//'. The solution may be inaccurate.'// &
                     ' Simulation continuing, but future warnings will be suppressed.')
      IF ( PRESENT(ErrTxt) ) THEN
         CALL WrScr(' Additional debugging message from SUBROUTINE SmllRotTrans(): '//TRIM(ErrTxt) )
      END IF

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

END MODULE NWTC_Num
