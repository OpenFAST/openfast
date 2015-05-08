MODULE AirfoilInfo


   ! This module contains airfoil-related routines with non-system-specific logic and references.

! Redo this routing to get rid of some of the phases.  For instance, AFI_Init should be calle directly.

   USE                                             AirfoilInfo_Types
   USE                                             NWTC_FitPack
   USE                                          :: ISO_FORTRAN_ENV , ONLY : IOSTAT_EOR

   IMPLICIT NONE

   PRIVATE

   PUBLIC                                       :: AFI_Init, ReadAFFile, AFI_GetAirfoilParams

   TYPE(ProgDesc), PARAMETER                    :: AFI_Ver = ProgDesc( 'AirfoilInfo', 'v1.00.00c-mlb', '03-Mar-2014')               ! The name, version, and date of AirfoilInfo.


CONTAINS

   ! SUBROUTINE AFI_Init         ( InitInput, p, ErrStat, ErrMsg )                                          ! Initialize.  Read in data.  Compute spline coefficients.
   ! SUBROUTINE AirfoilInfo      ( Phase, FileNum, AoA, Re, Ctrl, InitInput, Cl, Cd, Cm, Cpmin, ErrStat, ErrMsg )    ! Main interface routine.
   ! SUBROUTINE BEMT_GenSplines   ( FileNum, Re, Ctrl, InitInput, ErrStatLcl, ErrMsg )                                ! Phase #2: For BEM, generate splines for Cl and Cd at specific Re/Ctrl setting.
   ! SUBROUTINE BEMT_GetClCd      ( FileNum, AoA, Re, Ctrl, InitInput, Cl, Cd, ErrStatLcl, ErrMsg )                   ! Phase #3: For BEM, compute Cl and Cd for the given AoA.
   ! SUBROUTINE GetAllCoefs      ( FileNum, AoA, Re, Ctrl, InitInput, Cl, Cd, Cm, Cpmin, ErrStatLcl, ErrMsg )        ! Phase #4: Compute all requested airfoil coefficients for the given AoA.
   ! SUBROUTINE ReadAFfile       ( AFInfo, UA_Model, Col_Cm, Col_Cpmin, ErrStat, ErrMsg, InpErrs )             ! Read an airfoil file.

   !=============================================================================
   SUBROUTINE AFI_Init ( InitInput, p, ErrStat, ErrMsg, UnEc )


         ! This routine initializes AirfoilInfo by reading the airfoil files and generating the spline coefficients.


         ! Argument declarations.

      INTEGER(IntKi), INTENT(OUT)               :: ErrStat                    ! Error status.

      INTEGER, INTENT(IN), OPTIONAL             :: UnEc                       ! I/O unit for echo file. If present and > 0, write to UnEc.      CHARACTER(*), INTENT(IN)               :: AFfile                        ! The file to be read.

      CHARACTER(*), INTENT(OUT)                 :: ErrMsg                     ! Error message.

      TYPE (AFI_InitInputType), INTENT(INOUT)   :: InitInput                  ! This structure stores values that are set by the calling routine during the initialization phase.
      TYPE (AFI_ParameterType), INTENT(INOUT)   :: p                          ! This structure stores all the module parameters that are set by AirfoilInfo during the initialization phase.


         ! Local declarations.

      REAL                                   :: AlphaReg      (-180:180)      ! Regularly spaced alphas for interpolated data at the specified Re.
      REAL(SiKi), ALLOCATABLE                :: Coef          (:)             ! The coefficients to send to the regrid routine for 2D splines.
      REAL, ALLOCATABLE                      :: Coefs         (:,:,:,:)       ! Coefficients of cubic polynomials for the airfoil coefficients (maybe be two sets).
      REAL(ReKi)                             :: DelAlpha                      ! The distance between alpha values in AlphaReg.
      REAL, ALLOCATABLE                      :: IntAFCoefs    (:,:)           ! The interpolated airfoil coefficients.
      REAL, ALLOCATABLE                      :: IntAFCoefsHR  (:,:)           ! The interpolated airfoil coefficients in high resolution (0.1 deg).
      REAL                                   :: LogRe       = LOG( 1.75 )     ! LOG of Reynolds Number (in millions) to interpolate to.
      REAL, ALLOCATABLE                      :: LogReAry      (:)             ! Array of LOG( Re ) values to interpolate.
      REAL                                   :: Wt1                           ! Weighting factor for the low Re.
      REAL                                   :: Wt2                           ! Weighting factor for the high Re.

      INTEGER                                :: Co                            ! The index into the coefficients array.
      INTEGER                                :: CoefSiz                       ! The size of the Coef array.
      INTEGER                                :: ErrStatLcl                    ! Local error status.
      INTEGER                                :: File                          ! The file index.
      INTEGER                                :: I                             ! Index into the arrays.
      INTEGER                                :: IA                            ! Index into the alpha array.
      INTEGER                                :: Ind                           ! Index into the Coefs array.  Indicates which Re in the table is below the actual Re.
      INTEGER                                :: Indx                          ! Index into the arrays.
      INTEGER                                :: NumAl                         ! The length of the alpha array being sent to the regrid interface routine.
      INTEGER                                :: NumRe                         ! The length of the Re array being sent to the regrid interface routine.
      INTEGER                                :: Table                         ! Index into the Table array.

      TYPE (InpErrsType)                     :: InpErrs                       ! The derived type for holding input errors.


      ErrStat = ErrID_None
      ErrMsg  = ""

         ! Display the version for this module.

      CALL DispNVD ( AFI_Ver )

        

         ! Process the airfoil files.

      ALLOCATE ( p%AFInfo( InitInput%NumAFfiles ), STAT=ErrStatLcL )
      IF ( ErrStatLcl /= ErrID_None )  THEN
         CALL ExitThisRoutine ( ErrID_Fatal, ' >> Error allocating memory for the p%AFInfo array in AFI_Init.' )
         RETURN
      ENDIF

      InpErrs%NumErrs = 0

      DO File=1,InitInput%NumAFfiles

         IF ( PRESENT(UnEc) ) THEN
            IF ( UnEc > 0 )  THEN
               WRITE (UnEc,'("--",/,A)')  'Contents of "'//TRIM( InitInput%FileNames(File) )//'":'
            END IF
         END IF
         

         CALL ReadAFfile ( InitInput%FileNames(File), InitInput%UA_Model, InitInput%NumCoefs, InitInput%InCol_Alfa &
                         , InitInput%InCol_Cl, InitInput%InCol_Cd, InitInput%InCol_Cm, InitInput%InCol_Cpmin, p%AFInfo(File) &
                         , ErrStatLcl, ErrMsg, UnEc ) !, InpErrs )
         IF ( ErrStatLcl /= ErrID_None )  THEN
            CALL ExitThisRoutine ( ErrID_Fatal, ErrMsg )
            RETURN
         ENDIF


            ! Set the column indices for the variouse airfoil coefficients.

         p%ColCl = 1
         p%ColCd = 2
         IF ( InitInput%InCol_Cm > 0 )  THEN
            p%ColCm = 3
            IF ( InitInput%InCol_Cpmin > 0 )  THEN
               p%ColCpmin = 4
            END IF
         ELSE IF ( InitInput%InCol_Cpmin > 0 )  THEN
               p%ColCpmin = 3
         END IF


            ! Make sure that all the tables meet the current restrictions.

         IF ( p%AFInfo(File)%NumTabs > 1 )  THEN


               ! Do the tables have the same number and set of alphas.

            DO Table=2,p%AFInfo(File)%NumTabs

               IF ( p%AFInfo(File)%Table(Table)%NumAlf /= p%AFInfo(File)%Table(1)%NumAlf )  THEN
          !     IF ( p%AFInfo(File)%NumTabs /= p%AFInfo(File)%Table(1)%NumAlf )  THEN
                  CALL ExitThisRoutine ( ErrID_Fatal, '  >> Fatal Error: For airfoile file "'//TRIM( InitInput%FileNames(File) ) &
                                        //'", Table #'//TRIM( Num2LStr( Table ) ) &
                                        //' does not have the same set of alphas as the first table.' )
                  RETURN
               ENDIF

               DO IA=1,p%AFInfo(File)%Table(1)%NumAlf
                  IF ( p%AFInfo(File)%Table(Table)%Alpha(IA) /= p%AFInfo(File)%Table(1)%Alpha(IA) )  THEN
                     CALL ExitThisRoutine ( ErrID_Fatal, '  >> Fatal Error: For airfoile file "' &
                                           //TRIM( InitInput%FileNames(File) ) &
                                           //'", Table #'//TRIM( Num2LStr( Table ) ) &
                                           //' does not have the same set of alphas as the first table.' )
                     RETURN
                  ENDIF
               END DO ! IA

               IF ( p%AFInfo(File)%Table(Table)%Ctrl /= p%AFInfo(File)%Table(1)%Ctrl )  THEN
                  CALL ExitThisRoutine ( ErrID_Fatal, '  >> Fatal Error: For airfoile file "'//TRIM( InitInput%FileNames(File) ) &
                                        //'", Table #'//TRIM( Num2LStr( Table ) ) &
                                        //' does not have the same value for Ctrl as the first table.' )
                  RETURN
               ENDIF
            ENDDO ! Tab

         ENDIF ! ( p%AFInfo(File)%NumTabs > 1 )


            ! Create an array of Re for interpolation and fill it with the LOG(Re) for each table in the file.
!MLB: Maybe we should just store LogRe when we read in the Re.
         ALLOCATE ( p%AFInfo(File)%LogRe( p%AFInfo(File)%NumTabs ), STAT=ErrStatLcl )
         IF ( ErrStatLcl /= ErrID_None )  THEN
            CALL ExitThisRoutine ( ErrID_Fatal, '  >> Error allocating memory for the LogReAry array in BEMT_GenSplines.' )
            RETURN
         ENDIF

         DO Table=1,p%AFInfo(File)%NumTabs
            p%AFInfo(File)%LogRe(Table) = LOG( p%AFInfo(File)%Table(Table)%Re )
         END DO ! Table

! We need to deal with constant data.

            ! Generate the aifoil coefficients array and spline-coefficients arrays from the airfoil coefficients in a form usable by regrid.

         IF ( p%AFInfo(File)%NumTabs > 1 )  THEN                  ! We use 2D cubic spline interpolation.


               ! Generate the airfoil coefficients array from the airfoil coefficients in a form usable by regrid.
               ! This is a scratch array that will be reused for each airfoil coefficient.

            NumAl   = p%AFInfo(File)%Table(1)%NumAlf
            NumRe   = p%AFInfo(File)%NumTabs
    !        CoefSiz = NumAl*NumRe
            CoefSiz = NumAl*NumRe + 8
            CALL AllocAry ( Coef, CoefSiz, "array for bicubic-spline coefficients for airfoil data", ErrStatLcl, ErrMsg )
            IF ( ErrStatLcl /= 0 )  THEN
               CALL ExitThisRoutine( ErrID_FATAL, ErrMsg )
               RETURN
            END IF


               ! Compute the splines for the drag coefficient.

            CALL AllocAry ( p%AFInfo(File)%CdSpCoef2D, CoefSiz, "array for bicubic-spline coefficients for Cd data" &
                          , ErrStatLcl, ErrMsg )
            IF ( ErrStatLcl /= 0 )  THEN
               CALL ExitThisRoutine( ErrID_FATAL, ErrMsg )
               RETURN
            END IF

            Indx= 0
            DO IA=1,NumAl
               DO Table=1,NumRe
                  Indx       = Indx + 1
                  Coef(Indx) = p%AFInfo(File)%Table(Table)%Coefs(IA,p%ColCd)
               ENDDO ! Table
            ENDDO ! IA

            CALL FitPack_RegridLite ( NumAl                     , p%AFInfo(File)%Table(1)%Alpha &
                                    , NumRe                     , p%AFInfo(File)%LogRe &
                                    , CoefSiz                   , Coef &
                                    , p%AFInfo(File)%NumCdAoAkts, p%AFInfo(File)%CdAoAknots &
                                    , p%AFInfo(File)%NumCdReKts , p%AFInfo(File)%CdReKnots &
                                    , CoefSiz                   , p%AFInfo(File)%CdSpCoef2D &
                                    , ErrStatLcl                , ErrMsg )
            IF ( ErrStatLcl /= ErrID_None )  THEN
               CALL ExitThisRoutine ( ErrStatLcl, ErrMsg )
               RETURN
            ENDIF


               ! Compute the splines for the lift coefficient.

            CALL AllocAry ( p%AFInfo(File)%ClSpCoef2D, CoefSiz, "array for bicubic-spline coefficients for Cl data" &
                          , ErrStatLcl, ErrMsg )
            IF ( ErrStatLcl /= 0 )  THEN
               CALL ExitThisRoutine( ErrID_FATAL, ErrMsg )
               RETURN
            END IF

            Indx= 0
            DO IA=1,NumAl
               DO Table=1,NumRe
                  Indx       = Indx + 1
                  Coef(Indx) = p%AFInfo(File)%Table(Table)%Coefs(IA,p%ColCl)
               ENDDO ! Table
            ENDDO ! IA

            CALL FitPack_RegridLite ( NumAl                     , p%AFInfo(File)%Table(1)%Alpha &
                                    , NumRe                     , p%AFInfo(File)%LogRe &
                                    , CoefSiz                   , Coef &
                                    , p%AFInfo(File)%NumClAoAkts, p%AFInfo(File)%ClAoAknots &
                                    , p%AFInfo(File)%NumClReKts , p%AFInfo(File)%ClReKnots &
                                    , CoefSiz                   , p%AFInfo(File)%ClSpCoef2D &
                                    , ErrStatLcl                , ErrMsg )
            IF ( ErrStatLcl /= ErrID_None )  THEN
               CALL ExitThisRoutine ( ErrStatLcl, ErrMsg )
               RETURN
            ENDIF


               ! Optionally compute the splines for the pitching-moment coefficient.

            IF ( p%ColCm > 0 )  THEN

               CALL AllocAry ( p%AFInfo(File)%CmSpCoef2D, CoefSiz, "array for bicubic-spline coefficients for Cm data" &
                             , ErrStatLcl, ErrMsg )
               IF ( ErrStatLcl /= 0 )  THEN
                  CALL ExitThisRoutine( ErrID_FATAL, ErrMsg )
                  RETURN
               END IF

               Indx= 0
               DO IA=1,NumAl
                  DO Table=1,NumRe
                     Indx       = Indx + 1
                     Coef(Indx) = p%AFInfo(File)%Table(Table)%Coefs(IA,p%ColCm)
                  ENDDO ! Table
               ENDDO ! IA

               CALL FitPack_RegridLite ( NumAl                     , p%AFInfo(File)%Table(1)%Alpha &
                                       , NumRe                     , p%AFInfo(File)%LogRe &
                                       , CoefSiz                   , Coef &
                                       , p%AFInfo(File)%NumCmAoAkts, p%AFInfo(File)%CmAoAknots &
                                       , p%AFInfo(File)%NumCmReKts , p%AFInfo(File)%CmReKnots &
                                       , CoefSiz                   , p%AFInfo(File)%CmSpCoef2D &
                                       , ErrStatLcl                , ErrMsg )
               IF ( ErrStatLcl /= ErrID_None )  THEN
                  CALL ExitThisRoutine ( ErrStatLcl, ErrMsg )
                  RETURN
               ENDIF

            ENDIF ! ( p%ColCm > 0 )


               ! Optionally compute the splines for the minimum pressure coefficient.

            IF ( p%ColCpmin > 0 )  THEN

               CALL AllocAry ( p%AFInfo(File)%CpminSpCoef2D, CoefSiz, "array for bicubic-spline coefficients for Cpmin data" &
                             , ErrStatLcl, ErrMsg )
               IF ( ErrStatLcl /= 0 )  THEN
                  CALL ExitThisRoutine( ErrID_FATAL, ErrMsg )
                  RETURN
               END IF

               Indx= 0
               DO IA=1,NumAl
                  DO Table=1,NumRe
                     Indx       = Indx + 1
                     Coef(Indx) = p%AFInfo(File)%Table(Table)%Coefs(IA,p%ColCpmin)
                  ENDDO ! Table
               ENDDO ! IA

               CALL FitPack_RegridLite ( NumAl                        , p%AFInfo(File)%Table(1)%Alpha &
                                       , NumRe                        , p%AFInfo(File)%LogRe &
                                       , CoefSiz                      , Coef &
                                       , p%AFInfo(File)%NumCpminAoAkts, p%AFInfo(File)%CpminAoAknots &
                                       , p%AFInfo(File)%NumCpminReKts , p%AFInfo(File)%CpminReKnots &
                                       , CoefSiz                      , p%AFInfo(File)%CpminSpCoef2D &
                                       , ErrStatLcl                   , ErrMsg )
               IF ( ErrStatLcl /= ErrID_None )  THEN
                  CALL ExitThisRoutine ( ErrStatLcl, ErrMsg )
                  RETURN
               ENDIF

            ENDIF ! ( p%ColCpmin > 0 )


               ! Deallocate the temporary Coef array.

            IF ( ALLOCATED( Coef ) ) DEALLOCATE( Coef )


            ! Compute the spline coefficients of the piecewise cubic polynomials for the irregularly-spaced airfoil data in each file.
            ! Unless the data are constant.

            DO Co=1,InitInput%NumCoefs

!BUG: What happened to the 2D interpolation?  Maybe I haven't written this yet.

               CALL ExitThisRoutine ( ErrID_FATAL, ' >> The part to compute the bi-cubic splines in AFI_Init is not written yet!' )

            END DO ! Co

         ELSE ! We have only one table.  Use 1-D interpolation.


               ! Allocate the arrays to hold spline coefficients.

            ALLOCATE ( p%AFInfo(File)%Table(1)%SplineCoefs( p%AFInfo(File)%Table(1)%NumAlf-1 &
                     , InitInput%NumCoefs, 0:3 ), STAT=ErrStatLcl )
            IF ( ErrStatLcl /= ErrID_None )  THEN
               CALL ExitThisRoutine ( ErrStatLcl, ' >> Error allocating memory for the SplineCoefs array in AFI_Init.' )
               RETURN
            ENDIF


               ! Compute the one set of coefficients of the piecewise polynomials for the irregularly-spaced data.
               ! Unlike the 2-D interpolation in which we use diffent knots for each airfoil coefficient, we can do
               ! the 1-D stuff all at once.

            CALL CubicSplineInitM ( p%AFInfo(File)%Table(1)%Alpha &
                                  , p%AFInfo(File)%Table(1)%Coefs &
                                  , p%AFInfo(File)%Table(1)%SplineCoefs &
                                  , ErrStatLcl, ErrMsg )
            IF ( ErrStatLcl /= ErrID_None )  THEN
               CALL ExitThisRoutine ( ErrStatLcl, ErrMsg )
               RETURN
            ENDIF

         ENDIF ! ( p%AFInfo(File)%NumTabs > 1 )


            ! Compute the spline coefficients of the piecewise cubic polynomials for the irregularly-spaced airfoil data in each file.
            ! Unless the data are constant.

         DO Co=1,InitInput%NumCoefs


            IF ( p%AFInfo(File)%NumTabs > 1 )  THEN                  ! We use 2D cubic spline interpolation.

!BUG: What happened to the 2D interpolation?  Maybe I haven't written this yet.

               CALL ExitThisRoutine ( ErrID_FATAL, ' >> The part to compute the bi-cubic splines in AFI_Init is not written yet!' )

! We also need to deal with constant data.

            ELSE                                                     ! We use 1D cubic spline interpolation if the data are not constant.

               IF ( p%AFInfo(File)%Table(1)%ConstData )  THEN

                     ! We need to deal with constant data.

                  CALL ExitThisRoutine ( ErrID_FATAL, ' >> The part to deal with constant data in AFI_Init is not written yet!' )

               ENDIF ! p%AFInfo(File)%Table(1)%ConstData

            ENDIF ! ( p%AFInfo(File)%NumTabs > 1 )

         END DO ! Co

      END DO ! File

      CALL ExitThisRoutine ( ErrID_None, ErrMsg )

      RETURN

   !=======================================================================
   CONTAINS
   !=======================================================================
      SUBROUTINE ExitThisRoutine ( ErrID, Msg )

         ! This subroutine cleans up the parent routine before exiting.


            ! Argument declarations.

         INTEGER(IntKi), INTENT(IN)       :: ErrID                            ! The error identifier (ErrLev)

         CHARACTER(*),   INTENT(IN)       :: Msg                              ! The error message (ErrMsg)


            ! Set error status/message

         ErrStat = ErrID
         ErrMsg  = Msg


            ! Deallocate the temporary Coef array.

         IF ( ALLOCATED( Coef ) ) DEALLOCATE( Coef )


         RETURN

      END SUBROUTINE ExitThisRoutine ! ( ErrID, Msg )

   END SUBROUTINE AFI_Init ! ( InitInput, ErrStatLcl, ErrMsg )
   
  
   !=============================================================================
   SUBROUTINE ReadAFfile ( AFfile, UA_Model, NumCoefs, InCol_Alfa, InCol_Cl, InCol_Cd, InCol_Cm, InCol_Cpmin, AFInfo &
                         , ErrStat, ErrMsg, UnEc )!, InpErrs )


         ! This routine reads an airfoil file.


         ! Argument declarations.

      INTEGER(IntKi), INTENT(IN)             :: InCol_Alfa                    ! The airfoil-table input column for angle of attack.
      INTEGER(IntKi), INTENT(IN)             :: InCol_Cd                      ! The airfoil-table input column for drag coefficient.
      INTEGER(IntKi), INTENT(IN)             :: InCol_Cl                      ! The airfoil-table input column for lift coefficient.
      INTEGER(IntKi), INTENT(IN)             :: InCol_Cm                      ! The airfoil-table input column for pitching-moment coefficient.
      INTEGER(IntKi), INTENT(IN)             :: InCol_Cpmin                   ! The airfoil-table input column for minimum pressure coefficient.
      INTEGER(IntKi), INTENT(OUT)            :: ErrStat                       ! Error status.
      INTEGER(IntKi), INTENT(IN)             :: NumCoefs                      ! The number of aerodynamic coefficients to be stored.
      INTEGER(IntKi), INTENT(IN)             :: UA_Model                      ! The type of unsteady-aero model.

      INTEGER, INTENT(IN), OPTIONAL          :: UnEc                          ! I/O unit for echo file. If present and > 0, write to UnEc.      CHARACTER(*), INTENT(IN)               :: AFfile                        ! The file to be read.

      CHARACTER(*), INTENT( IN)              :: AFfile                        ! The file to read in.
      CHARACTER(*), INTENT(OUT)              :: ErrMsg                        ! Error message.

      TYPE (AFInfoType), INTENT(INOUT)       :: AFInfo                        ! The derived type for holding the constant parameters for this airfoil.
!      TYPE (InpErrsType), INTENT(INOUT)      :: InpErrs                       ! The derived type for holding input errors.


         ! Local declarations.

      REAL(SiKi)                             :: Coords   (2)                  ! An array to hold data from the airfoil-shape table.
      REAL(SiKi), ALLOCATABLE                :: SiAry    (:)                  ! A temporary array to hold data from a table.

      INTEGER                                :: Coef                          ! A DO index that points into the coefficient array.
      INTEGER                                :: Cols2Parse                    ! The number of columns that must be read from the coefficient tables.
      INTEGER                                :: CurLine                       ! The current line to be parsed in the FileInfo structure.
      INTEGER(IntKi)                         :: ErrStatLcl                    ! Error status local to this routine.
      INTEGER                                :: NumAlf                        ! The number of rows in the current table.
      INTEGER                                :: Row                           ! The row of a table to be parsed in the FileInfo structure.
      INTEGER                                :: Table                         ! The DO index for the tables.

      LOGICAL                                :: BadVals                       ! A flag that indicates if the values in a table are invalid.

      TYPE (FileInfoType)                    :: FileInfo                      ! The derived type for holding the file information.

      ErrStat = ErrID_None
      ErrMsg  = ""


         ! Process the (possibly) nested set of files.  This copies the decommented contents of
         ! AFI_FileInfo%FileName and the files it includes (both directly and indirectly) into
         ! the FileInfo structure that we can then parse.

      CALL ProcessComFile ( AFfile, FileInfo, ErrStatLcl, ErrMsg )
      IF ( ErrStatLcl > 0 )  THEN
         CALL ExitThisRoutine ( ErrID_Fatal, ErrMsg )
         RETURN
      ENDIF


         ! Process the airfoil shape information if it is included.

      CurLine = 1


         ! These first two parameters are optional, so we don't check for errors.

      CALL ParseVar ( FileInfo, CurLine, 'NonDimArea', AFInfo%NonDimArea, ErrStatLcl, ErrMsg, UnEc )
      CALL ParseVar ( FileInfo, CurLine, 'NumCoords' , AFInfo%NumCoords , ErrStatLcl, ErrMsg, UnEc )

      IF ( AFInfo%NumCoords > 0 )  THEN

         ALLOCATE ( AFInfo%X_Coord( AFInfo%NumCoords ) , STAT=ErrStatLcl )
         IF ( ErrStatLcl /= 0 )  THEN
            CALL ExitThisRoutine ( ErrID_Fatal, '  >> Error allocating memory for the AFInfo%X_Coord array in ReadAFfile.' )
            RETURN
         ENDIF

         ALLOCATE ( AFInfo%Y_Coord( AFInfo%NumCoords ) , STAT=ErrStatLcl )
         IF ( ErrStatLcl /= 0 )  THEN
            CALL ExitThisRoutine ( ErrID_Fatal, '  >> Error allocating memory for the AFInfo%Y_Coord array in ReadAFfile.'  )
            RETURN
         ENDIF

         DO Row=1,AFInfo%NumCoords
            CALL ParseAry ( FileInfo, CurLine, 'X_Coord/Y_Coord', Coords, 2, ErrStatLcl, ErrMsg, UnEc )
            IF ( ErrStatLcl /= 0 )  THEN
               CALL ExitThisRoutine ( ErrID_Fatal, ErrMsg )
               RETURN
            ENDIF
            AFInfo%X_Coord(Row) = Coords(1)
            AFInfo%Y_Coord(Row) = Coords(2)
         ENDDO ! Row

      ENDIF


         ! How many columns do we need to read in the input and how many total coefficients will be used?

      Cols2Parse = MAX( InCol_Alfa, InCol_Cl, InCol_Cd, InCol_Cm, InCol_Cpmin )
      ALLOCATE ( SiAry( Cols2Parse ) , STAT=ErrStatLcl )
      IF ( ErrStatLcl /= 0 )  THEN
         CALL ExitThisRoutine ( ErrID_Fatal, ' >> Error allocating memory for the SiAry array in ReadAFfile.' )
         RETURN
      ENDIF


         ! Work through the multiple tables.

      CALL ParseVar ( FileInfo, CurLine, 'NumTabs' , AFInfo%NumTabs , ErrStatLcl, ErrMsg, UnEc )
      IF ( ErrStatLcl /= 0 )  THEN
         CALL ExitThisRoutine ( ErrID_Fatal, ErrMsg )
         RETURN
      ENDIF

      IF ( AFInfo%NumTabs < 1 )  THEN
         CALL ExitThisRoutine ( ErrID_Fatal, ' >> NumTabs must be > 0 in "'//TRIM( AFfile )//'".' )
         RETURN
      ENDIF

      ALLOCATE ( AFInfo%Table( AFInfo%NumTabs ) , STAT=ErrStatLcl )
      IF ( ErrStatLcl /= 0 )  THEN
         CALL ExitThisRoutine ( ErrID_Fatal, ' >> Error allocating memory for the AFInfo%Table array in ReadAFfile.' )
         RETURN
      ENDIF

      DO Table=1,AFInfo%NumTabs

         CALL ParseVar ( FileInfo, CurLine, 'Re', AFInfo%Table(Table)%Re, ErrStatLcl, ErrMsg, UnEc )
         IF ( ErrStatLcl /= 0 )  THEN
            CALL ExitThisRoutine ( ErrID_Fatal, ErrMsg )
            RETURN
         ENDIF
         IF ( AFInfo%Table(Table)%Re <= 0.0 )  THEN
!            IF ( InpErrs%NumErrs < InpErrs%MaxErrs )  THEN
!               InpErrs%NumErrs                   = InpErrs%NumErrs + 1
!               InpErrs%FileList(InpErrs%NumErrs) = AFfile
!               InpErrs%FileLine(InpErrs%NumErrs) = FileInfo%FileLine(CurLine-1)
!               InpErrs%ErrMsgs (InpErrs%NumErrs) = 'Re must be > 0'
!            ELSE
               CALL ExitThisRoutine ( ErrID_Severe, '  >> Re must be > 0 in "'//TRIM( AFfile ) &
                                                  //'".'//NewLine//'  >> The error occurred on line #' &
                                                  //TRIM( Num2LStr( FileInfo%FileLine(CurLine-1) ) )//'.' )
               RETURN
!            ENDIF ! ( InpErrs%NumErrs < InpErrs%MaxErrs )
         ENDIF ! ( AFInfo%Table(Table)%Re <= 0.0 )

         CALL ParseVar ( FileInfo, CurLine, 'Ctrl', AFInfo%Table(Table)%Ctrl, ErrStatLcl, ErrMsg, UnEc )
         IF ( ErrStatLcl /= 0 )  THEN
            CALL ExitThisRoutine ( ErrID_Fatal, ErrMsg )
            RETURN
         ENDIF
         IF ( AFInfo%Table(Table)%Ctrl /= 0.0 )  THEN
!            IF ( InpErrs%NumErrs < InpErrs%MaxErrs )  THEN
!               InpErrs%NumErrs                   = InpErrs%NumErrs + 1
!               InpErrs%FileList(InpErrs%NumErrs) = AFfile
!               InpErrs%FileLine(InpErrs%NumErrs) = FileInfo%FileLine(CurLine-1)
!               InpErrs%ErrMsgs (InpErrs%NumErrs) = 'Ctrl must equal 0.0'
!            ELSE
               CALL ExitThisRoutine ( ErrID_Severe, '  >> Ctrl must equal 0.0 in "'//TRIM( AFfile ) &
                                                  //'".'//NewLine//'  >> The error occurred on line #' &
                                                  //TRIM( Num2LStr( FileInfo%FileLine(CurLine-1) ) )//'.' )
               RETURN
!            ENDIF ! ( InpErrs%NumErrs < InpErrs%MaxErrs )
         ENDIF ! ( AFInfo%Table(Table)%Ctrl <= 0.0 )

         CALL ParseVar ( FileInfo, CurLine, 'InclUAdata', AFInfo%Table(Table)%InclUAdata, ErrStatLcl, ErrMsg, UnEc )
         IF ( ErrStatLcl /= 0 )  THEN
            CALL ExitThisRoutine ( ErrID_Fatal, ErrMsg )
            RETURN
         ENDIF

         IF ( AFInfo%Table(Table)%InclUAdata )  THEN

            IF ( UA_Model == 1 )  THEN

               CALL ParseVar ( FileInfo, CurLine, 'BL_AOL', AFInfo%Table(Table)%UA_BL%BL_AOL, ErrStatLcl, ErrMsg, UnEc )
               IF ( ErrStatLcl /= 0 )  THEN
                  CALL ExitThisRoutine ( ErrID_Fatal, ErrMsg )
                  RETURN
               ENDIF

               CALL ParseVar ( FileInfo, CurLine, 'BL_CnA', AFInfo%Table(Table)%UA_BL%BL_CnA, ErrStatLcl, ErrMsg, UnEc )
               IF ( ErrStatLcl /= 0 )  THEN
                  CALL ExitThisRoutine ( ErrID_Fatal, ErrMsg )
                  RETURN
               ENDIF

               CALL ParseVar ( FileInfo, CurLine, 'BL_CnS', AFInfo%Table(Table)%UA_BL%BL_CnS, ErrStatLcl, ErrMsg, UnEc )
               IF ( ErrStatLcl /= 0 )  THEN
                  CALL ExitThisRoutine ( ErrID_Fatal, ErrMsg )
                  RETURN
               ENDIF

               CALL ParseVar ( FileInfo, CurLine, 'BL_CnSL', AFInfo%Table(Table)%UA_BL%BL_CnSL, ErrStatLcl, ErrMsg, UnEc )
               IF ( ErrStatLcl /= 0 )  THEN
                  CALL ExitThisRoutine ( ErrID_Fatal, ErrMsg )
                  RETURN
               ENDIF

               CALL ParseVar ( FileInfo, CurLine, 'BL_AOD', AFInfo%Table(Table)%UA_BL%BL_AOD, ErrStatLcl, ErrMsg, UnEc )
               IF ( ErrStatLcl /= 0 )  THEN
                  CALL ExitThisRoutine ( ErrID_Fatal, ErrMsg )
                  RETURN
               ENDIF

               CALL ParseVar ( FileInfo, CurLine, 'BL_Cd0', AFInfo%Table(Table)%UA_BL%BL_Cd0, ErrStatLcl, ErrMsg, UnEc )
               IF ( ErrStatLcl /= 0 )  THEN
                  CALL ExitThisRoutine ( ErrID_Fatal, ErrMsg )
                  RETURN
               ENDIF

            ENDIF ! ( UA_Model == 1 )

         ELSEIF ( UA_Model == 1 )  THEN

            CALL ExitThisRoutine ( ErrID_Fatal &
            , ' >> You must supply Beddoes-Leishman unsteady aerodynamics parameters for all airfoils if you want to use that' &
            //' model.  You did not do so for Table #'//TRIM( Num2LStr( Table ) )//' in the "'//TRIM( AFfile )//'" airfoil file.' )
            RETURN

         ENDIF ! ( AFInfo%Table(Table)%InclUAdata )

         CALL ParseVar ( FileInfo, CurLine, 'NumAlf', AFInfo%Table(Table)%NumAlf, ErrStatLcl, ErrMsg, UnEc )
         IF ( ErrStatLcl /= 0 )  THEN
            CALL ExitThisRoutine ( ErrID_Fatal, ErrMsg )
            RETURN
         ENDIF

         IF ( AFInfo%Table(Table)%NumAlf < 1 )  THEN
            CALL ExitThisRoutine( ErrID_Fatal, '  >> NumAlf must be a positive number on line #' &
                           //TRIM( Num2LStr( FileInfo%FileLine(CurLine-1) ) )//' in "'//TRIM( AFfile )//'".' )
            RETURN
         ELSEIF ( AFInfo%Table(Table)%NumAlf < 3 )  THEN
            AFInfo%Table(Table)%ConstData = .TRUE.
         ELSE
            AFInfo%Table(Table)%ConstData = .FALSE.
         ENDIF ! ( Test for valid values for NumAlf )


            ! Allocate the arrays for the airfoil coefficients.

         ALLOCATE ( AFInfo%Table(Table)%Alpha( AFInfo%Table(Table)%NumAlf ), STAT=ErrStatLcl )
         IF ( ErrStatLcl /= 0 )  THEN
            CALL ExitThisRoutine( ErrID_Fatal, '  >> Error allocating memory for the AFInfo%Table%Alpha array in ReadAFfile.' )
            RETURN
         ENDIF

         ALLOCATE ( AFInfo%Table(Table)%Coefs( AFInfo%Table(Table)%NumAlf, NumCoefs ), STAT=ErrStatLcl )
         IF ( ErrStatLcl /= 0 )  THEN
            CALL ExitThisRoutine( ErrID_Fatal, '  >> Error allocating memory for the AFInfo%Table%Coefs array in ReadAFfile.' )
            RETURN
         ENDIF

         DO Row=1,AFInfo%Table(Table)%NumAlf

            CALL ParseAry ( FileInfo, CurLine, 'CoeffData', SiAry, Cols2Parse, ErrStatLcl, ErrMsg, UnEc )
            IF ( ErrStatLcl /= 0 )  THEN
               CALL ExitThisRoutine( ErrID_Fatal, ErrMsg )
               RETURN
            ENDIF

   !mlb: testing         AFInfo%Table(Table)%Alpha(Row  ) = SiAry(InCol_Alfa)*D2R
            AFInfo%Table(Table)%Alpha(Row  ) = SiAry(InCol_Alfa)
            AFInfo%Table(Table)%Coefs(Row,1) = SiAry(InCol_Cl  )
            AFInfo%Table(Table)%Coefs(Row,2) = SiAry(InCol_Cd  )

            IF ( InCol_Cm > 0 )  THEN
               AFInfo%Table(Table)%Coefs(Row,3) = SiAry(InCol_Cm)
               IF ( InCol_Cpmin > 0 )  AFInfo%Table(Table)%Coefs(Row,4) = SiAry(InCol_Cpmin)
            ELSE
               IF ( InCol_Cpmin > 0 )  AFInfo%Table(Table)%Coefs(Row,3) = SiAry(InCol_Cpmin)
            ENDIF ! IF ( Col_Cm > 0 )  THEN

         ENDDO ! Row


            ! Let's make sure that the data go from -Pi to Pi and that the values are the same for both
            ! unless there is only one point.

         IF ( .NOT. AFInfo%Table(Table)%ConstData )  THEN
            NumAlf  = AFInfo%Table(Table)%NumAlf
            BadVals = .FALSE.
  !mlb: testing          IF ( .NOT. EqualRealNos( AFInfo%Table(Table)%Alpha(1), -Pi ) )  THEN
            IF ( .NOT. EqualRealNos( AFInfo%Table(Table)%Alpha(1), -180.0 ) )  THEN
               BadVals = .TRUE.
            ENDIF
  !mlb: testing          IF ( .NOT. EqualRealNos( AFInfo%Table(Table)%Alpha(NumAlf), Pi ) )  THEN
            IF ( .NOT. EqualRealNos( AFInfo%Table(Table)%Alpha(NumAlf), 180.0 ) )  THEN
               BadVals = .TRUE.
            ENDIF
            DO Coef=1,NumCoefs
               IF ( .NOT. EqualRealNos( AFInfo%Table(Table)%Coefs(1,Coef), AFInfo%Table(Table)%Coefs(NumAlf,Coef) ) )  THEN
                  BadVals = .TRUE.
               ENDIF
            ENDDO ! Coef
            IF ( BadVals )  THEN
               CALL ExitThisRoutine( ErrID_Warn, &
                  ' >> Airfoil data should go from -180 to 180 and the coefficients at the ends should be the same.' )
               RETURN
            ENDIF
         ENDIF ! ( .NOT. AFInfo%Table(Table)%ConstData )

      ENDDO ! Table


      CALL ExitThisRoutine( ErrID_None, '' )

      RETURN

      !=======================================================================
      CONTAINS
      !=======================================================================
         SUBROUTINE ExitThisRoutine ( ErrID, Msg )

            ! This subroutine cleans up all the allocatable arrays, sets the error status/message and closes the binary file

               ! Passed arguments.

            INTEGER(IntKi), INTENT(IN)     :: ErrID        ! The error identifier (ErrLev)

            CHARACTER(*),   INTENT(IN)     :: Msg          ! The error message (ErrMsg)


               ! Set error status/message

            ErrStat = ErrID
            ErrMsg  = Msg


               ! If there was an error, deallocate the arrays in the FileInfo structure.

            IF ( ErrStat /= 0 )  THEN
               IF ( ALLOCATED( FileInfo%FileLine ) ) DEALLOCATE( FileInfo%FileLine )
               IF ( ALLOCATED( FileInfo%FileIndx ) ) DEALLOCATE( FileInfo%FileIndx )
               IF ( ALLOCATED( FileInfo%FileList ) ) DEALLOCATE( FileInfo%FileList )
               IF ( ALLOCATED( FileInfo%Lines    ) ) DEALLOCATE( FileInfo%Lines    )
            END IF ! ( ErrLev /= 0 )

         END SUBROUTINE ExitThisRoutine ! ( ErrID, Msg )


      END SUBROUTINE ReadAFfile ! ( FileInfo, ErrStat, ErrMsg )
      
      
subroutine AFI_GetAirfoilParams( AFInfo, M, Re, alpha, alpha0, alpha1, alpha2, eta_e, C_nalpha, C_nalpha_circ, T_f0, T_V0, T_p, T_VL, St_sh, &
                                 b1, b2, b5, A1, A2, A5, S1, S2, S3, S4, Cn1, Cn2, Cd0, Cm0, k0, k1, k2, k3, k1_hat, x_cp_bar, errMsg, errStat )     

   type(AFInfoType), intent(in   )       :: AFInfo                        ! The derived type for holding the constant parameters for this airfoil.
   real(ReKi),       intent(in   )       :: M                             ! mach number
   real(ReKi),       intent(in   )       :: Re                            ! Reynold's number
   real(ReKi),       intent(in   )       :: alpha                         !
   real(ReKi),       intent(  out)       :: alpha0                        ! zero lift angle of attack (radians)
   real(ReKi),       intent(  out)       :: alpha1                        ! angle of attack at f = 0.7, approximately the stall angle; for alpha >= alpha0 (radians)
   real(ReKi),       intent(  out)       :: alpha2                        ! angle of attack at f = 0.7, approximately the stall angle; for alpha < alpha0 (radians)
   real(ReKi),       intent(  out)       :: eta_e                         !
   real(ReKi),       intent(  out)       :: C_nalpha                      !
   real(ReKi),       intent(  out)       :: C_nalpha_circ                 ! slope of the circulatory normal force coefficient vs alpha curve
   real(ReKi),       intent(  out)       :: T_f0                          ! initial value of T_f, airfoil specific, used to compute D_f and fprimeprime
   real(ReKi),       intent(  out)       :: T_V0                          ! initial value of T_V, airfoil specific, time parameter associated with the vortex lift decay process, used in Cn_v
   real(ReKi),       intent(  out)       :: T_p                           ! boundary-layer, leading edge pressure gradient time parameter; used in D_p; airfoil specific
   real(ReKi),       intent(  out)       :: T_VL                          ! time variable associated with the vortex advection process; it represents the non-dimensional time in semi-chords needed for a vortex to travel from leading edge to trailing edge
   real(ReKi),       intent(  out)       :: St_sh                         !
   real(ReKi),       intent(  out)       :: b1                            ! airfoil constant derived from experimental results, usually 0.14
   real(ReKi),       intent(  out)       :: b2                            ! airfoil constant derived from experimental results, usually 0.53
   real(ReKi),       intent(  out)       :: b5                            ! airfoil constant derived from experimental results, usually 5.0
   real(ReKi),       intent(  out)       :: A1                            ! airfoil constant derived from experimental results, usually 0.3 
   real(ReKi),       intent(  out)       :: A2                            ! airfoil constant derived from experimental results, usually 0.7
   real(ReKi),       intent(  out)       :: A5                            ! airfoil constant derived from experimental results, usually 1.0
   real(ReKi),       intent(  out)       :: S1                            ! constant in the f-curve best-fit, alpha >= alpha0 
   real(ReKi),       intent(  out)       :: S2                            ! constant in the f-curve best-fit, alpha >= alpha0 
   real(ReKi),       intent(  out)       :: S3                            ! constant in the f-curve best-fit, alpha <  alpha0 
   real(ReKi),       intent(  out)       :: S4                            ! constant in the f-curve best-fit, alpha <  alpha0 
   real(ReKi),       intent(  out)       :: Cn1                           ! critical value of Cn_prime at LE separation for alpha >= alpha0
   real(ReKi),       intent(  out)       :: Cn2                           ! critical value of Cn_prime at LE separation for alpha < alpha0
   real(ReKi),       intent(  out)       :: Cd0                           !
   real(ReKi),       intent(  out)       :: Cm0                           ! 2D pitching moment coefficient at zero lift, positive if nose is up
   real(ReKi),       intent(  out)       :: k0                            ! airfoil parameter in the x_cp_hat curve best-fit
   real(ReKi),       intent(  out)       :: k1                            ! airfoil parameter in the x_cp_hat curve best-fit
   real(ReKi),       intent(  out)       :: k2                            ! airfoil parameter in the x_cp_hat curve best-fit
   real(ReKi),       intent(  out)       :: k3                            ! airfoil parameter in the x_cp_hat curve best-fit
   real(ReKi),       intent(  out)       :: k1_hat                        !
   real(ReKi),       intent(  out)       :: x_cp_bar                      ! airfoil parameter for calulating x_cp_v
   integer(IntKi),   intent(  out)       :: errStat                       ! Error status. 
   character(*),     intent(  out)       :: errMsg                        ! Error message.
      
   errMsg         = ''
   errStat        = ErrID_None
   
   alpha0         =  -0.35_ReKi*pi/180.0
   alpha1         =  11.0_ReKi*pi/180.0
   alpha2         =  -11.0_ReKi*pi/180.0
   eta_e          =  0.90               ! Recovery factor in the range [0.85 - 0.95]
   C_nalpha       =  2*pi 
   C_nalpha_circ  =  C_nalpha / sqrt(1.0_ReKi-M**2)
   T_f0           =  3.0_ReKi  ! seconds
   T_V0           =  6.0_ReKi
   T_p            =  1.7_ReKi
   T_VL           =  11.0_ReKi
   b1             =  0.14_ReKi
   b2             =  0.53_ReKi
   b5             =  5.0_ReKi
   A1             =  0.3_ReKi
   A2             =  0.70_ReKi
   A5             =  1.0_ReKi
   S1             =  0.0262_ReKi   !!!!!!!!!!
   S2             =  0.0201_ReKi   !!!!!!!!!!
   S3             =  -0.0262_ReKi   !!!!!!!!!!
   S4             =  -0.0201_ReKi   !!!!!!!!!!
   Cn1            =  1.264_ReKi  ! Stall values of Cn
   Cn2            =  -0.833_ReKi
   St_sh          =  0.19_ReKi
   Cd0            =  0.012_ReKi
   Cm0            =  0.0_ReKi
   k0             =  0.0_ReKi
   k1             =  0.0_ReKi
   k2             =  0.0_ReKi
   k3             =  0.0_ReKi
   k1_hat         =  0.0_ReKi
   x_cp_bar       =  0.2_ReKi
   
  ! Cn1=1.9 Tp=1.7 Tf=3., Tv=6 Tvl=11, Cd0=0.012
   
end subroutine AFI_GetAirfoilParams
                                 
!=============================================================================
END MODULE AirfoilInfo
