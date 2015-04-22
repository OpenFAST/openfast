PROGRAM Test_RegCubicSpline


      ! This program is used to test the routines to interpolate regularly spaced arrays using cubic splines.


   USE                                    :: NWTC_Library

   IMPLICIT                                  NONE


      ! Type declarations.

   INTEGER, PARAMETER                        :: AryLen     = 10               ! Index into the arrays.

   REAL(ReKi)                                :: X                             ! The value to interpolate for.
   REAL(ReKi)                                :: Y                             ! The interpolated value from   regularly spaced data.
   REAL(ReKi)                                :: Z                             ! The interpolated value from irregularly spaced data.

   INTEGER                                   :: I                             ! Index into the arrays.
   INTEGER(IntKi)                            :: ErrStat                       ! Error status.

   CHARACTER(4096)                           :: ErrMsg                        ! Error message.

   TYPE(   CubSplineType)                    ::    CubSplineData              ! Derived type to hold data for irregularly spaced cubic splines.
   TYPE(RegCubSplineType)                    :: RegCubSplineData              ! Derived type to hold data for regularly spaced cubic splines.



      ! Initialize the NWTC Library.

   CALL NWTC_Init ( 'Test_RegCubicSpline', 'v1.00.00, 17-May-2013', EchoLibVer=.FALSE. )
   CALL WrScr     ( NewLine//' Running Test_RegCubicSpline (v1.00.00, 17-May-2013)' )


      ! Set up the regularly spaced data to be interpolated.

   RegCubSplineData%NumPts = 10

   ALLOCATE ( RegCubSplineData%XAry( RegCubSplineData%NumPts ), STAT=ErrStat )
   IF ( ErrStat /= 0 )  THEN
      CALL ProgAbort ( ' >> Error allocating memory for the RegCubSplineData%XAry array in Test_RegCubicSpline.', .FALSE., 1.0, &
                       ErrStat )
   ENDIF

   ALLOCATE ( RegCubSplineData%YAry( RegCubSplineData%NumPts ), STAT=ErrStat )
   IF ( ErrStat /= 0 )  THEN
      CALL ProgAbort ( ' >> Error allocating memory for the RegCubSplineData%YAry array in Test_RegCubicSpline.', .FALSE., 1.0, &
                       ErrStat )
   ENDIF

   ALLOCATE ( RegCubSplineData%Coef( RegCubSplineData%NumPts-1, 0:3 ), STAT=ErrStat )
   IF ( ErrStat /= 0 )  THEN
      CALL ProgAbort ( ' >> Error allocating memory for the RegCubSplineData%Coef array in Test_RegCubicSpline.', .FALSE., 1.0, &
                       ErrStat )
   ENDIF


      ! Set up the irregularly spaced data to be interpolated.

   CubSplineData%NumPts = RegCubSplineData%NumPts

   ALLOCATE ( CubSplineData%XAry( CubSplineData%NumPts ), STAT=ErrStat )
   IF ( ErrStat /= 0 )  THEN
      CALL ProgAbort ( ' >> Error allocating memory for the CubSplineData%XAry array in Test_RegCubicSpline.', .FALSE., 1.0, &
                       ErrStat )
   ENDIF

   ALLOCATE ( CubSplineData%YAry( CubSplineData%NumPts ), STAT=ErrStat )
   IF ( ErrStat /= 0 )  THEN
      CALL ProgAbort ( ' >> Error allocating memory for the CubSplineData%YAry array in Test_RegCubicSpline.', .FALSE., 1.0, &
                       ErrStat )
   ENDIF

   ALLOCATE ( CubSplineData%Coef( CubSplineData%NumPts-1, 0:3 ), STAT=ErrStat )
   IF ( ErrStat /= 0 )  THEN
      CALL ProgAbort ( ' >> Error allocating memory for the CubSplineData%Coef array in Test_RegCubicSpline.', .FALSE., 1.0, &
                       ErrStat )
   ENDIF


      ! Create an array to interpolate.

   PRINT *

   DO I=1,AryLen
      RegCubSplineData%XAry(I) = REAL( 2*I - 2 )*Pi/REAL( AryLen - 1 )           ! Values from 0 to 2*Pi in steps of 2*Pi/(AryLen-1)
      RegCubSplineData%YAry(I) = SIN( RegCubSplineData%XAry(I) )
         CubSplineData%XAry(I) = RegCubSplineData%XAry(I)
         CubSplineData%YAry(I) = RegCubSplineData%YAry(I)
      PRINT '(I5,1X,F12.4)', NINT( R2D*RegCubSplineData%XAry(I) ), RegCubSplineData%YAry(I)
   ENDDO ! I


      ! Compute the coefficients of the piecewise polynomials for the regularly-spaced data.

   CALL RegCubicSplineInit ( RegCubSplineData%NumPts, RegCubSplineData%XAry, RegCubSplineData%YAry, RegCubSplineData%DelX &
                           , RegCubSplineData%Coef, ErrStat, ErrMsg )
   IF ( ErrStat == ErrID_Fatal )  THEN
      CALL ProgAbort ( ErrMsg, .FALSE., 1.0, ErrStat )
   ELSEIF ( ErrStat /= ErrID_None )  THEN
      CALL WrScr ( ErrMsg )
   ENDIF


      ! Compute the coefficients of the piecewise polynomials for the irregularly-spaced data.

   CALL CubicSplineInit ( CubSplineData%NumPts, CubSplineData%XAry, CubSplineData%YAry &
                                              , CubSplineData%Coef, ErrStat, ErrMsg )
   IF ( ErrStat == ErrID_Fatal )  THEN
      CALL ProgAbort ( ErrMsg, .FALSE., 1.0, ErrStat )
   ELSEIF ( ErrStat /= ErrID_None )  THEN
      CALL WrScr ( ErrMsg )
   ENDIF


      ! Write out the interpolated values for both sets of interpolated data.

   PRINT *

   DO I=0,360,1
      X = REAL( I )*D2R
      Y = RegCubicSplineInterp( X, RegCubSplineData%NumPts, RegCubSplineData%XAry, RegCubSplineData%YAry &
                                 , RegCubSplineData%DelX  , RegCubSplineData%Coef, ErrStat, ErrMsg )
      Z = CubicSplineInterp   ( X,    CubSplineData%NumPts,    CubSplineData%XAry,    CubSplineData%YAry &
                                                          ,    CubSplineData%Coef, ErrStat, ErrMsg )
      IF ( ErrStat == ErrID_Fatal )  THEN
         CALL ProgAbort ( ErrMsg, .FALSE., 1.0, ErrStat )
      ELSEIF ( ErrStat /= ErrID_None )  THEN
         CALL WrScr ( ErrMsg )
      ENDIF
      PRINT '(I5,2(1X,F12.4))', I, Y, Z
   ENDDO ! I


   CALL NormStop

END PROGRAM Test_RegCubicSpline