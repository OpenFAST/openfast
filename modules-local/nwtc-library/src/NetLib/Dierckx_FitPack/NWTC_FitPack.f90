!**********************************************************************************************************************************
! Copyright (C) 2014  National Renewable Energy Laboratory
!
! This code provides a wrapper for Dierckx's FitPack routines currently used at the NWTC (mainly codes in the FAST framework).
!
!**********************************************************************************************************************************
MODULE NWTC_FitPack


      ! This Fortran module consists of single-precision wrappers for some of Paul Dierckx's FitPack routines.


   USE NWTC_Library
!   USE, INTRINSIC               :: ISO_C_Binding, only: C_FLOAT, C_DOUBLE          ! this is included in NWTC_Library

      ! Notes:

         ! Your project must include the following files:
         ! From the NWTC Subroutine Library:
         !     SingPrec.f90          [from NWTC Library]
         !     Sys*.f90              [from NWTC Library]
         !     NWTC_Base.f90         [from NWTC Library]
         ! fitpack library (preferably a binary, but available in source form from http://www.netlib.org/dierckx/, too).
         ! This wrapper file:
         !     NWTC_FitPack.f90

   IMPLICIT  NONE

   INTERFACE FitPack_BispevLite                 ! Computes bicubic splines of a rectangular grid of data.
!      MODULE PROCEDURE FitPack_dBispevLite      ! Dummy routine that prints an error message if double precision is used.
      MODULE PROCEDURE FitPack_sBispevLite      ! Single-precision version.
   END INTERFACE

   INTERFACE FitPack_RegridLite                 ! Computes bicubic splines of a rectangular grid of data.
      MODULE PROCEDURE FitPack_dRegridLite      ! Dummy routine that prints an error message if double precision is used.
      MODULE PROCEDURE FitPack_sRegridLite      ! Single-precision version.
   END INTERFACE

!   INTERFACE FitPack_surfit                     ! Fits a 2D surface to a randomly located set of points.
!      MODULE PROCEDURE FitPack_dsurfit
!      MODULE PROCEDURE FitPack_ssurfit
!   END INTERFACE


CONTAINS

!=======================================================================
   SUBROUTINE FitPack_sBispevLite ( NumX, X, NumY, Y, SCsize, SplCoef, ErrStat, ErrMsg )
!      subroutine bispev(tx,nx,ty,ny,c,kx,ky,x,mx,y,my,z,wrk,lwrk,
!     * iwrk,kwrk,ier)


      ! This Fortran routine is a wrapper for Paul Dierckx's FitPack bispev.f routine.
      ! This works only with single-precision values (SiKi).
      ! It *may* assume that the length of REALs and INTEGERs are the same.
      ! Hard-coded assumptions:
      !     The calling routine passes it a rectangular grid of data stored as a vector.
      ! You can blame Marshall Buhl for this mess.

      ! This routine calls bispev and bispev uses the following routines either directly or indirectly:
      !     fpbisp and fpbspl.
      !
      !  bispev
      !     fpbisp
      !        fpbspl

      USE                                          NWTC_Base                  ! We only need the precision and error level constants

      IMPLICIT                                     NONE


         ! Argument declarations:

      INTEGER,                 INTENT(IN)     :: NumX                ! The length of the X array.
      INTEGER,                 INTENT(IN)     :: NumY                ! The length of the Y array.
      INTEGER,                 INTENT(IN)     :: SCsize              ! The length of the SplCoef array.

      REAL(SiKi),              INTENT(OUT)    :: SplCoef(SCsize)     ! The resulting spline coefficients.
      REAL(SiKi),              INTENT(IN)     :: X      (NumX)       ! X coordinates for the Elev array.
      REAL(SiKi),              INTENT(IN)     :: Y      (NumY)       ! Y coordinates for the Elev array.
                                 
      INTEGER,                 INTENT(OUT)    :: ErrStat             ! The error status to be returned to the calling program.
                                 
      CHARACTER(*),            INTENT(OUT)    :: ErrMsg              ! Error message.


         ! Local declarations.

!NOTE:  For the surface elevation, regrid() needs a rank-1 array with the Y values varying most rapidly.  It would be best to make that for AoA and use X for Control, which there will most often be only one value.

      REAL(SiKi)                              :: ResidSq                   ! The sum of squared residuals of the spline approximation
      REAL(SiKi), ALLOCATABLE                 :: ReWorkAry  (:)            ! A working array of type REAL.
!MLB: Should SmthFact be a user-specified parameter?  It almost seems like it would be dependent on airfoils.  Andrew used 0.1 for lift and 0.2 for drag because of the relative sizes of the curves.
!     We decided that for the first cut, force it to zero (no smoothing), but maybe eventually making it user-specified.
      REAL(SiKi), PARAMETER                   :: SmthFact   = 0.0_SiKi     ! The non-negative smoothing factor.  Hard-coded to 0 for no smoothing.
      REAL(SiKi), ALLOCATABLE                 :: Xknots     (:)            ! The array of knots in the X direction.
      REAL(SiKi), ALLOCATABLE                 :: Yknots     (:)            ! The array of knots in the Y direction.

      INTEGER, PARAMETER                      :: CubicSpl   = 3            ! The degree of the cubic-spline polynomials.
      INTEGER(IntKi)                          :: ErrStatLcl                ! The local version of the error status.
      INTEGER(IntKi)                          :: I                         ! A DO counter.
      INTEGER(IntKi), ALLOCATABLE             :: InWorkAry  (:)            ! A working array of type INTEGER.
      INTEGER(IntKi)                          :: LenIWrkAry                ! The size of InWorkAry.
      INTEGER(IntKi)                          :: LenRWrkAry                ! The size of ReWorkAry.
      INTEGER(IntKi)                          :: MaxXknots                 ! The maximum number of knots in the X direction.
      INTEGER(IntKi)                          :: MaxYknots                 ! The maximum number of knots in the Y direction.
      INTEGER(IntKi)                          :: NumXknots                 ! The number of knots in the X direction.
      INTEGER(IntKi)                          :: NumYknots                 ! The number of knots in the Y direction.
!mlb: testing      INTEGER(IntKi), PARAMETER               :: Option     = -1           ! The degree of the spline polynomials.
      INTEGER(IntKi), PARAMETER               :: Option     = 1            ! The degree of the spline polynomials.
      INTEGER                                 :: SplDegX                   ! The degree of the spline polynomials in the X direction.
      INTEGER                                 :: SplDegY                   ! The degree of the spline polynomials in the Y direction.



         ! Check to see if this program was compiled with non-standard REALs (not four bytes).
         ! We do this because the FitPack routines *may* require the word length of REALs and INTEGERs to be the same.
                                  
      IF ( KIND( 1.0 ) /= SiKi )  THEN
         CALL ExitThisRoutine( ErrID_SEVERE, ' The FitPack routines have not been tested with REALs that are not four bytes.' )
         RETURN
      END IF       

      SplCoef = 0.0  


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


            ! Deallocate the work arrays.

         IF ( ALLOCATED( ReWorkAry ) ) DEALLOCATE( ReWorkAry )
         IF ( ALLOCATED( InWorkAry ) ) DEALLOCATE( InWorkAry )


            ! Set error status/message

         ErrStat = ErrID
         ErrMsg  = Msg


         RETURN

      END SUBROUTINE ExitThisRoutine ! ( ErrID, Msg )

   END SUBROUTINE FitPack_sBispevLite
!=======================================================================
   SUBROUTINE FitPack_dRegridLite ( NumX, X, NumY, Y, ElSize, Elev, SCsize, SplCoef, ErrStat, ErrMsg )


      ! This Fortran routine is a wrapper for Paul Dierckx's FitPack regrid.f routine.
      ! This works only with single-precision values (SiKi).
      ! It *may* assume that the length of REALs and INTEGERs are the same.
      ! Hard-coded assumptions:
      !     The calling routine passes it a rectangular grid of data.
      !     The knots will be at the original X,Y points (Option = -1).
      !     It will generate cubic splines in both directions.
      !     No smoothing will be done.
      ! You can blame Marshall Buhl for this mess.

      ! This routine calls regrid and regrid uses the following routines either directly or indirectly:
      !     fpback, fpbspl, fpchec, fpdisc, fpgivs, fpgrre, fpknot, fpregr, and fprota.
      !
      !  regrid
      !     fpchec
      !     fpregr
      !        fpgrre
      !           fpbspl
      !           fpdisc
      !           fpgivs
      !           fprota
      !           fpback
      !        fpknot

      USE                                          NWTC_Base                  ! We only need the precision and error level constants

      IMPLICIT                                     NONE


         ! Argument declarations:

      INTEGER,          INTENT(IN)     :: ElSize              ! The length of the Elev array.
      INTEGER,          INTENT(IN)     :: NumX                ! The length of the X array.
      INTEGER,          INTENT(IN)     :: NumY                ! The length of the Y array.
      INTEGER,          INTENT(IN)     :: SCsize              ! The length of the SplCoef array.

      REAL(R8Ki),       INTENT(IN)     :: Elev   (ElSize)     ! The "elevations" for the x and y coordinates for which spline coefficients will be computed.
      REAL(R8Ki),       INTENT(OUT)    :: SplCoef(SCsize)     ! The resulting spline coefficients.
      REAL(R8Ki),       INTENT(IN)     :: X      (NumX)       ! X coordinates for the Elev array.
      REAL(R8Ki),       INTENT(IN)     :: Y      (NumY)       ! Y coordinates for the Elev array.
                          
      INTEGER,          INTENT(OUT)    :: ErrStat             ! The error status to be returned to the calling program.
      INTEGER, PARAMETER               :: SplDeg     = 3      ! The degree of the spline polynomials.
                          
      CHARACTER(*),     INTENT(OUT)    :: ErrMsg              ! Error message.
!NOTE:  For the surface elevation, regrid() needs a rank-1 array with the Y values varying most rapidly.  It would be best to make that for AoA and use X for Control, which there will most often be only one value.


!         ! Local declarations.
!
!      REAL(DbKi)                                :: ResidSq                    ! The sum of squared residuals of the spline approximation
!      REAL(DbKi), ALLOCATABLE                   :: ReWorkAry  (:)             ! A working array of type REAL.
!!MLB: Should SmthFact be a user-specified parameter?  It almost seems like it would be dependent on airfoils.  Andrew used 0.1 for lift and 0.2 for drag because of the relative sizes of the curves.
!!     We decided that for the first cut, force it to zero (no smoothing), but maybe eventually making it user-specified.
!      REAL(DbKi), PARAMETER                     :: SmthFact   = 0.0_SiKi      ! The non-negative smoothing factor.  Hard-coded to 0 for no smoothing.
!
!      INTEGER(IntKi)                            :: ErrStatLcl                 ! The local version of the error status.
!      INTEGER(IntKi), ALLOCATABLE               :: InWorkAry  (:)             ! A working array of type INTEGER.
!      INTEGER(IntKi)                            :: LenIWrkAry                 ! The size of InWorkAry.
!      INTEGER(IntKi)                            :: LenRWrkAry                 ! The size of ReWorkAry.
!      INTEGER(IntKi)                            :: LenSpCoAry                 ! The size of SplCoef.
!      INTEGER(IntKi)                            :: MaxXknots                  ! The maximum number of knots in the X direction.
!      INTEGER(IntKi)                            :: MaxYknots                  ! The maximum number of knots in the Y direction.
!      INTEGER(IntKi), PARAMETER                 :: Option     = -1            ! The degree of the spline polynomials.
!
!
         ! This routine is just a dummy.

      SplCoef(:) = 0.0_R8Ki

      CALL ExitThisRoutine ( ErrID_Fatal, ' Error >> FitPack_RegridLite is not available for 8-byte real precision.' )
!
!         ! Check to see if this program was compiled with non-standard REALs (not four bytes).
!         ! We do this because the FitPack routines *may* require the word length of REALs and INTEGERs to be the same.
!                                  
!      IF ( KIND( 1.0 ) /= SiKi )  THEN
!         CALL ExitThisRoutine( ErrID_SEVERE, ' The FitPack routines have not been tested with REALs that are not four bytes.' )
!      END IF         
!
!
!         ! Let's set some array limits.
!
!      MaxXknots  = NumX + SplDeg + 1
!      MaxYknots  = NumY + SplDeg + 1
!      LenIWrkAry = 3 + NumX + NumY + MaxXknots + MaxYknots 
!      LenRWrkAry = 4 + MaxXknots*( NumY + 2*SplDeg + 5 ) + MaxYknots*( 2*SplDeg + 5 ) + NumX*( SplDeg + 1 ) &
!                 + NumY*( SplDeg + 1 ) + MAX( NumY, MaxXknots )
!      LenSpCoAry = ( MaxXknots - SplDeg - 1 )*( MaxYknots - SplDeg - 1 )
!
!
!         ! Allocate the arrays.
!
!      CALL AllocAry ( SplCoef, LenSpCoAry, "array for spline coefficients for airfoil data" &
!                    , ErrStatLcl, ErrMsg )
!      IF ( ErrStatLcl /= 0 )  THEN
!         CALL ExitThisRoutine( ErrID_FATAL, ErrMsg )
!      END IF
!
!      CALL AllocAry ( InWorkAry, LenIWrkAry, 'the INTEGER work array for regrid', ErrStatLcL, ErrMsg )
!      IF ( ErrStatLcl /= 0 )  THEN
!         CALL ExitThisRoutine( ErrID_FATAL, ErrMsg )
!      END IF
!
!      CALL AllocAry ( ReWorkAry, LenRWrkAry, 'the REAL work array for regrid', ErrStatLcL, ErrMsg )
!      IF ( ErrStatLcl /= 0 )  THEN
!         CALL ExitThisRoutine( ErrID_FATAL, ErrMsg )
!      END IF
!
!
!         ! Compute the spline coefficients.
!
!      CALL regrid ( iopt=Option, mx=NumX, x=X, my=NumY, y=Y, z=Elev, xb=X(1), xe=X(NumX), yb=Y(1), ye=Y(NumY), kx=SplDeg &
!                  , ky=SplDeg, s=SmthFact, nxest=MaxXknots, nyest=MaxYknots, nx=NumX, tx=X, ny=NumY, ty=Y, c=SplCoef &
!                  , fp=ResidSq, wrk=ReWorkAry, lwrk=LenRWrkAry, iwrk=InWorkAry, kwrk=LenIWrkAry, ier=ErrStatLcl )
!
!
!         ! Check for errors.
!
!      IF ( ErrStatLcl == 1 )  THEN
!         CALL ExitThisRoutine( ErrID_FATAL, 'regrid: The required storage space exceeds the available storage space, as ' &
!                                          //'specified by the parameters nxest and nyest.' )
!      ELSEIF ( ErrStatLcl == 2 )  THEN
!         CALL ExitThisRoutine( ErrID_FATAL, 'regrid: A theoretically impossible result was found during the iteration process ' &
!                                          //'for finding a smoothing spline.' )
!      ELSEIF ( ErrStatLcl == 3 )  THEN
!         CALL ExitThisRoutine( ErrID_FATAL, 'regrid: The maximal number of iterations maxit (set to 20 by the program) allowed ' &
!                                          //'for finding a smoothing spline with fp=s has been reached.' )
!      ELSEIF ( ErrStatLcl == 10 )  THEN
!         CALL ExitThisRoutine( ErrID_FATAL, 'regrid: One or more input conditions were invalid.' )
!      ELSE
!         CALL ExitThisRoutine( ErrID_None, '' )
!      END IF


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


            ! Deallocate the work arrays.

!         IF ( ALLOCATED( ReWorkAry ) ) DEALLOCATE( ReWorkAry )
!         IF ( ALLOCATED( InWorkAry ) ) DEALLOCATE( InWorkAry )


            ! Set error status/message

         ErrStat = ErrID
         ErrMsg  = Msg


         RETURN

      END SUBROUTINE ExitThisRoutine ! ( ErrID, Msg )

   END SUBROUTINE FitPack_dRegridLite
!=======================================================================
   SUBROUTINE FitPack_sRegridLite ( NumX     , X &
                                  , NumY     , Y &
                                  , ElSize   , Elev &
                                  , NumXknots, Xknots &
                                  , NumYknots, Yknots &
                                  , SCsize   , SplCoef &
                                  , ErrStat  , ErrMsg )


      ! This Fortran routine is a wrapper for Paul Dierckx's FitPack regrid.f routine.
      ! This works only with single-precision values (SiKi).
      ! It *may* assume that the length of REALs and INTEGERs are the same.
      ! Hard-coded assumptions:
      !     The calling routine passes it a rectangular grid of data stored as a vector.
      !     The knots will be at the original X,Y points (Option = -1).  MLB
      !     It will generate cubic splines in both directions unless there are fewer than four points.
      !        If there are less than four points in a given direction, the spline degree will be NumPts-1.
      !     No smoothing will be done.
      ! You can blame Marshall Buhl for this mess.

      ! This routine calls regrid and regrid uses the following routines either directly or indirectly:
      !     fpback, fpbspl, fpchec, fpdisc, fpgivs, fpgrre, fpknot, fpregr, and fprota.
      !
      !  regrid
      !     fpchec
      !     fpregr
      !        fpgrre
      !           fpbspl
      !           fpdisc
      !           fpgivs
      !           fprota
      !           fpback
      !        fpknot

      USE                                          NWTC_Base                  ! We only need the precision and error level constants

      IMPLICIT                                     NONE


         ! Argument declarations:

      INTEGER,      INTENT(IN)                  :: ElSize              ! The length of the Elev array.
      INTEGER,      INTENT(IN)                  :: NumX                ! The length of the X array.
      INTEGER,      INTENT(IN)                  :: NumY                ! The length of the Y array.
      INTEGER,      INTENT(IN)                  :: SCsize              ! The length of the SplCoef array.
                                                
      REAL(ReKi),   INTENT(IN)                  :: Elev   (ElSize)     ! The "elevations" for the x and y coordinates for which spline coefficients will be computed.
      REAL(ReKi),   INTENT(OUT)                 :: SplCoef(SCsize)     ! The resulting spline coefficients.
      REAL(ReKi),   INTENT(IN)                  :: X      (NumX)       ! X coordinates for the Elev array.
      REAL(ReKi),   INTENT(IN)                  :: Y      (NumY)       ! Y coordinates for the Elev array.
      REAL(ReKi),   INTENT(OUT), ALLOCATABLE    :: Xknots (:)          ! The computed knots in the X direction.
      REAL(ReKi),   INTENT(OUT), ALLOCATABLE    :: Yknots (:)          ! The computed knots in the Y direction.
                                   
      INTEGER,      INTENT(OUT)                 :: ErrStat             ! The error status to be returned to the calling program.
      INTEGER,      INTENT(OUT)                 :: NumXknots           ! The number of knots in the X direction.
      INTEGER,      INTENT(OUT)                 :: NumYknots           ! The number of knots in the Y direction.
                                                
      CHARACTER(*), INTENT(OUT)                 :: ErrMsg              ! Error message.


         ! Local declarations.

!NOTE:  For the surface elevation, regrid() needs a rank-1 array with the Y values varying most rapidly.  It would be best to make that for AoA and use X for Control, which there will most often be only one value.

      REAL(ReKi)                              :: ResidSq                   ! The sum of squared residuals of the spline approximation
      REAL(ReKi), ALLOCATABLE                 :: ReWorkAry  (:)            ! A working array of type REAL.
!MLB: Should SmthFact be a user-specified parameter?  It almost seems like it would be dependent on airfoils.  Andrew used 0.1 for lift and 0.2 for drag because of the relative sizes of the curves.
!     We decided that for the first cut, force it to zero (no smoothing), but maybe eventually making it user-specified.
      REAL(ReKi), PARAMETER                   :: SmthFact   = 0.0_ReKi     ! The non-negative smoothing factor.  Hard-coded to 0 for no smoothing.

      INTEGER, PARAMETER                      :: CubicSpl   = 3            ! The degree of the cubic-spline polynomials.
      INTEGER(IntKi)                          :: ErrStatLcl                ! The local version of the error status.
      INTEGER(IntKi)                          :: I                         ! A DO counter.
      INTEGER(IntKi), ALLOCATABLE             :: InWorkAry  (:)            ! A working array of type INTEGER.
      INTEGER(IntKi)                          :: LenIWrkAry                ! The size of InWorkAry.
      INTEGER(IntKi)                          :: LenRWrkAry                ! The size of ReWorkAry.
      INTEGER(IntKi)                          :: MaxXknots                 ! The maximum number of knots in the X direction.
      INTEGER(IntKi)                          :: MaxYknots                 ! The maximum number of knots in the Y direction.
!mlb: testing      INTEGER(IntKi), PARAMETER               :: Option     = -1           ! The degree of the spline polynomials.
      INTEGER(IntKi), PARAMETER               :: Option     = 1            ! The degree of the spline polynomials.
      INTEGER                                 :: SplDegX                   ! The degree of the spline polynomials in the X direction.
      INTEGER                                 :: SplDegY                   ! The degree of the spline polynomials in the Y direction.



         ! Check to see if this program was compiled with non-standard REALs (not four bytes).
         ! We do this because the FitPack routines *may* require the word length of REALs and INTEGERs to be the same.
                                  
      !IF ( KIND( 1.0 ) /= SiKi )  THEN
      !   CALL ExitThisRoutine( ErrID_SEVERE, ' The FitPack routines have not been tested with REALs that are not four bytes.' )
      !   RETURN
      !END IF         


         ! Determine what degree of polynomial to use for each direction.

      SplDegX = MIN( CubicSpl, NumX-1 )
      SplDegY = MIN( CubicSpl, NumY-1 )


         ! Let's set some array limits.

      MaxXknots  = NumX + SplDegX + 1
      MaxYknots  = NumY + SplDegY + 1
!      MaxXknots  = NumX + 2*( SplDegX + 1 )
!      MaxYknots  = NumY + 2*( SplDegY + 1 )
!      NumXknots  = NumX
!      NumYknots  = NumY
      LenIWrkAry = 3 + NumX + NumY + MaxXknots + MaxYknots 
      LenRWrkAry = 4 + MaxXknots*( NumY + 2*SplDegX + 5 ) + MaxYknots*( 2*SplDegY + 5 ) + NumX*( SplDegX + 1 ) &
                 + NumY*( SplDegY + 1 ) + MAX( NumY, MaxXknots )


         ! Allocate the working arrays.

      CALL AllocAry ( Xknots, MaxXknots, 'the array of knots in the X direction for regrid', ErrStatLcL, ErrMsg )
      IF ( ErrStatLcl /= 0 )  THEN
         CALL ExitThisRoutine( ErrID_FATAL, ErrMsg )
         RETURN
      END IF

      CALL AllocAry ( Yknots, MaxYknots, 'the array of knots in the Y direction for regrid', ErrStatLcL, ErrMsg )
      IF ( ErrStatLcl /= 0 )  THEN
         CALL ExitThisRoutine( ErrID_FATAL, ErrMsg )
         RETURN
      END IF

      CALL AllocAry ( InWorkAry, LenIWrkAry, 'the INTEGER work array for regrid', ErrStatLcL, ErrMsg )
      IF ( ErrStatLcl /= 0 )  THEN
         CALL ExitThisRoutine( ErrID_FATAL, ErrMsg )
         RETURN
      END IF

      CALL AllocAry ( ReWorkAry, LenRWrkAry, 'the REAL work array for regrid', ErrStatLcL, ErrMsg )
      IF ( ErrStatLcl /= 0 )  THEN
         CALL ExitThisRoutine( ErrID_FATAL, ErrMsg )
         RETURN
      END IF


         ! Prefill the arrays of knots.

!      DO I=1,NumXknots
!         Xknots(I+SplDegX+1) = X(I)
!      END DO ! I
!
!      DO I=1,NumYknots
!         Yknots(I+SplDegY+1) = Y(I)
!      END DO ! I


         ! Compute the spline coefficients.

      !CALL regrid ( iopt=Option, mx=NumX, x=X, my=NumY, y=Y, z=Elev, xb=X(1), xe=X(NumX), yb=Y(1), ye=Y(NumY), kx=SplDegX &
      !            , ky=SplDegY, s=SmthFact, nxest=MaxXknots, nyest=MaxYknots, nx=NumXknots, tx=Xknots, ny=NumYknots, ty=Yknots &
      !            , c=SplCoef, fp=ResidSq, wrk=ReWorkAry, lwrk=LenRWrkAry, iwrk=InWorkAry, kwrk=LenIWrkAry, ier=ErrStatLcl )
      CALL regrid ( Option, NumX, X, NumY, Y, Elev, X(1), X(NumX), Y(1), Y(NumY), SplDegX &
                  , SplDegY, SmthFact, MaxXknots, MaxYknots, NumXknots, Xknots, NumYknots, Yknots &
                  , SplCoef, ResidSq, ReWorkAry, LenRWrkAry, InWorkAry, LenIWrkAry, ErrStatLcl )


         ! Check for errors.

      IF ( ErrStatLcl == 1 )  THEN
         CALL ExitThisRoutine( ErrID_FATAL &
                             , 'Error >> regrid: The required storage space exceeds the available storage space, as ' &
                                          //'specified by the parameters MaxXknots and MaxYknots.' )
      ELSEIF ( ErrStatLcl == 2 )  THEN
         CALL ExitThisRoutine( ErrID_FATAL &
                             , 'Error >> regrid: A theoretically impossible result was found during the iteration process ' &
                             //'for finding a smoothing spline.' )
      ELSEIF ( ErrStatLcl == 3 )  THEN
         CALL ExitThisRoutine( ErrID_FATAL &
                             , 'Error >> regrid: The maximal number of iterations maxit (set to 20 by the program) allowed ' &
                             //'for finding a smoothing spline with ResidSq=SmthFact has been reached.' )
      ELSEIF ( ErrStatLcl == 10 )  THEN
         CALL ExitThisRoutine( ErrID_FATAL &
                             , 'Error >> regrid: One or more input conditions were invalid.' )
      ELSE
         CALL ExitThisRoutine( ErrID_None, '' )
      END IF


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


            ! Deallocate the work arrays.

         IF ( ALLOCATED( ReWorkAry ) ) DEALLOCATE( ReWorkAry )
         IF ( ALLOCATED( InWorkAry ) ) DEALLOCATE( InWorkAry )


            ! Set error status/message

         ErrStat = ErrID
         ErrMsg  = Msg


         RETURN

      END SUBROUTINE ExitThisRoutine ! ( ErrID, Msg )

   END SUBROUTINE FitPack_sRegridLite
!=======================================================================
!   SUBROUTINE FitPack_dsurfit ( Option, NumPts, X, Y, Z, Wts, Xbeg, Xend, Ybeg, Yend, SmthFact, nxest,nyest, nmax, Eps, Tx, Ty, nxest, nyest, quiet, ErrStat, ErrMsg )
!      subroutine surfit(iopt,m,x,y,z,w,xb,xe,yb,ye,kx,ky,s,nxest,nyest,
!     *  nmax,eps,nx,tx,ny,ty,c,fp,wrk1,lwrk1,wrk2,lwrk2,iwrk,kwrk,ier)
!
!
!      ! This Fortran routine is a wrapper for Paul Dierckx's FitPack surfit.f routine.
!      ! Effort on this routine was abandoned in favor of requiring users to specify rectangular grids of aerodynamic coefficients.
!
!      ! Find a bivariate B-spline representation of a surface.
!
!      ! Given a set of data points (x[i], y[i], z[i]) representing a surface
!      ! z=f(x,y), compute a B-spline representation of the surface. Based on
!      ! the routine SURFIT from FITPACK.
!
!      ! If Task=0, find knots in X and Y and coefficients for a given smoothing factor, S. [default]
!      ! If Task=1, find knots and coefficients for another value of the smoothing factor, S. bisplrep must have been previously called with task=0 or task=1.
!      ! If Task=-1, find coefficients for a given set of knots Tx, Ty.
!
!
!      USE                                          NWTC_Base                  ! We only need the precision and error level constants
!
!      IMPLICIT                                     NONE
!
!
!         ! Argument declarations:
!
!      INTEGER(IntKi),       INTENT(IN)          :: NumPts                     ! The length of the W, X, Y, and Z arrays.
!
!      REAL(DbKi), OPTIONAL, INTENT(IN)          :: Eps                        ! A threshold for determining the effective rank of an over-determined linear system of equations (0 < eps < 1).
!                                                                              !     Eps is not likely to need changing.  Default is 1.0e-16.
!      REAL(DbKi), OPTIONAL, INTENT(IN)          :: SmthFact                   ! A non-negative smoothing factor. If weights correspond to the inverse of the standard-deviation of the errors in Z,
!                                                                              !     then a good value for SmthFact would be found in the range (NumPts-SQRT(2*NumPts),NumPts+SQRT(2*NumPts)).
!      REAL(DbKi), OPTIONAL, INTENT(IN)          :: Tx(NumPts)                 ! The arrays of the X-direction knots of the spline for task=-1.
!      REAL(DbKi), OPTIONAL, INTENT(IN)          :: Ty(NumPts)                 ! The arrays of the Y-direction knots of the spline for task=-1.
!      REAL(DbKi), OPTIONAL, INTENT(IN)          :: Ty(NumPts)                 ! Weights of the data points for fitting a smooth surface.
!      REAL(DbKi), OPTIONAL, INTENT(IN)          :: Wts(NumPts)                ! Weights of the data points for fitting a smooth surface.
!      REAL(DbKi),           INTENT(IN)          :: X (NumPts)                 ! X coordinates for random data points to fit a surface to.
!      REAL(DbKi),           INTENT(IN)          :: Xbeg                       ! Beginning of the X region to be approximated.  Default is the minimum X value.
!      REAL(DbKi),           INTENT(IN)          :: Xend                       ! Ending of the X region to be approximated.  Default is the maximum X value.
!      REAL(DbKi),           INTENT(IN)          :: Y (NumPts)                 ! Y coordinates for random data points to fit a surface to.
!      REAL(DbKi),           INTENT(IN)          :: Ybeg                       ! Beginning of the Y region to be approximated.  Default is the minimum X value.
!      REAL(DbKi),           INTENT(IN)          :: Yend                       ! Ending of the Y region to be approximated.  Default is the maximum Y value.
!      REAL(DbKi),           INTENT(IN)          :: Z (NumPts)                 ! Z is f(X,Y) values of the random data points to fit a surface to.
!
!      INTEGER(IntKi), OPTIONAL, INTENT(IN)      :: Kx                         ! The degree of the splines in the X direction.  Default is 3.
!      INTEGER(IntKi), OPTIONAL, INTENT(IN)      :: Ky                         ! The degree of the splines in the Y direction.  Default is 3.
!      INTEGER(IntKi),           INTENT(IN)      :: Option                     ! The type of task to do.  Possibilities are -1, 0, 1.  Default is 0.
!
!      CHARACTER(*),             INTENT(OUT)     :: ErrMsg                     ! Error message.
!
!
!         ! Local declarations.
!
!      REAL(DbKi), OPTIONAL, INTENT(IN)          :: Epslon                     ! The local version of the threshold for determining the effective rank of an over-determined linear system of equations (0 < eps < 1).
!      REAL(DbKi)                                :: SmthFact                   ! The local version of the non-negative smoothing factor.
!      REAL(DbKi), ALLOCATABLE                   :: Wts(:)                     ! The local version of the weights.
!
!      INTEGER(IntKi)                            :: ErrStatLcl                 ! The local version of the error status.
!      INTEGER(IntKi)                            :: SplDegX                    ! The local version of the degree of the splines in the X direction.
!      INTEGER(IntKi)                            :: SplDegY                    ! The local version of the degree of the splines in the Y direction.
!
!
!
!      ErrStat = ErrID_None
!      ErrMsg  = ""
!
!
!      IF ( .NOT. PRESENT( W ) )  THEN
!         CALL AllocAry ( Wts, NumPts, 'Array of weights in BiSplRep', ErrStatLcl, ErrMsg )
!         IF ( ErrStatLcl /= 0 )  THEN
!            CALL ExitThisRoutine ( ErrID_Fatal, ErrMsg )
!         ENDIF
!         Wts(:) = 1.0
!      ENDIF
!
!      IF ( PRESENT( Kx ) )  THEN
!         SplDegX = Kx
!      ELSE
!         SplDegX = 3
!      ENDIF
!
!      IF ( PRESENT( Ky ) )  THEN
!         SplDegY = Ky
!      ELSE
!         SplDegY = 3
!      ENDIF
!
!      CALL DSurFit ( Opt, NumPts, X, Y, Z, Wts, Xbeg, Xend,Ybeg,Yend,kx,ky,SmthFact, nmax,eps,nx,tx,ny,ty,c,fp,wrk1,lwrk1,wrk2,lwrk2,iwrk,kwrk,ier)
!
!
!
!      IF (INFO /= 0) THEN
!         ErrStat = ErrID_FATAL
!         WRITE( ErrMsg, * ) INFO
!         IF (INFO < 0) THEN
!            ErrMsg  = "LAPACK_DGBSV: illegal value in argument "//TRIM(ErrMsg)//"."
!         ELSE
!            ErrMsg = 'LAPACK_DGBSV: U('//TRIM(ErrMsg)//','//TRIM(ErrMsg)//')=0. Factor U is exactly singular.'
!         END IF
!      END IF
!
!
!   RETURN
!   END SUBROUTINE FitPack_dsurfit
!!=======================================================================
!   SUBROUTINE FitPack_ssurfit( N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, ErrStat, ErrMsg )
!
!
!      ! passed parameters
!
!      INTEGER,         intent(in   ) :: KL                ! The number of subdiagonals within the band of A.  KL >= 0.
!      INTEGER,         intent(in   ) :: KU                ! The number of superdiagonals within the band of A.  KU >= 0.
!      INTEGER,         intent(in   ) :: N                 ! The number of linear equations, i.e., the order of the matrix A.  N >= 0.
!      INTEGER,         intent(in   ) :: NRHS              ! The number of right hand sides, i.e., the number of columns of the matrix B.  NRHS >= 0.
!      INTEGER,         intent(in   ) :: LDAB              ! The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
!      INTEGER,         intent(in   ) :: LDB               ! The leading dimension of the array B.  LDB >= max(1,N).
!
!      !     .. Array Arguments ..
!      REAL(SiKi)      ,intent(inout) :: AB( LDAB, * )      ! On entry, the matrix A in band storage, in rows KL+1 to 2*KL+KU+1; rows 1 to KL of the array need not be set.
!                                                          ! The j-th column of A is stored in the j-th column of the array AB as follows:
!                                                          !    AB(KL+KU+1+i-j,j) = A(i,j) for max(1,j-KU)<=i<=min(N,j+KL)
!                                                          ! On exit, details of the factorization: U is stored as an upper triangular band matrix with KL+KU superdiagonals in
!                                                          ! rows 1 to KL+KU+1, and the multipliers used during the factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
!      REAL(SiKi)      ,intent(inout) :: B( LDB, * )       ! On entry, the N-by-NRHS right hand side matrix B.  On exit, if INFO = 0, the N-by-NRHS solution matrix X.
!      INTEGER,         intent(  out) :: IPIV( * )         ! The pivot indices that define the permutation matrix P; row i of the matrix was interchanged with row IPIV(i).
!
!      INTEGER(IntKi),  intent(  out) :: ErrStat           ! Error level
!      CHARACTER(*),    intent(  out) :: ErrMsg            ! Message describing error
!
!         ! local variables
!      INTEGER                        :: INFO              ! = 0:  successful exit; < 0:  if INFO = -i, the i-th argument had an illegal value; > 0: if INFO = i, U(i,i) is exactly zero.
!                                                          ! > 0: if INFO = i, U(i,i) is exactly zero. The factorization has been completed, but the factor U is exactly singular, so the solution could not be computed.
!
!
!
!      ErrStat = ErrID_None
!      ErrMsg  = ""
!
!      CALL SGBSV (N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO)
!
!      IF (INFO /= 0) THEN
!         ErrStat = ErrID_FATAL
!         WRITE( ErrMsg, * ) INFO
!         IF (INFO < 0) THEN
!            ErrMsg  = "LAPACK_SGBSV: illegal value in argument "//TRIM(ErrMsg)//"."
!         ELSE
!            ErrMsg = 'LAPACK_SGBSV: U('//TRIM(ErrMsg)//','//TRIM(ErrMsg)//')=0. Factor U is exactly singular.'
!         END IF
!      END IF
!
!
!   RETURN
!   END SUBROUTINE FitPack_ssurfit
!=======================================================================
END MODULE NWTC_FitPack
