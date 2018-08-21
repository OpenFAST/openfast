!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2015-2016  National Renewable Energy Laboratory
!
!    This file is part of AeroDyn.
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
MODULE AirfoilInfo


   ! This module contains airfoil-related routines with non-system-specific logic and references.

! Redo this routing to get rid of some of the phases.  For instance, AFI_Init should be calle directly.

   USE                                             AirfoilInfo_Types
   USE                                             NWTC_FitPack
   USE                                          :: ISO_FORTRAN_ENV , ONLY : IOSTAT_EOR

   IMPLICIT NONE

   PRIVATE

   PUBLIC                                       :: AFI_Init ! routine to initialize AirfoilInfo parameters
   PUBLIC                                       :: AFI_GetAirfoilParams ! routine to calculate Airfoil parameters

   TYPE(ProgDesc), PARAMETER                    :: AFI_Ver = ProgDesc( 'AirfoilInfo', '', '')               ! The name, version, and date of AirfoilInfo.


CONTAINS

   ! SUBROUTINE AFI_Init         ( InitInput, p, ErrStat, ErrMsg )                                          ! Initialize.  Read in data.  Compute spline coefficients.
   ! SUBROUTINE AirfoilInfo      ( Phase, FileNum, AoA, Re, Ctrl, InitInput, Cl, Cd, Cm, Cpmin, ErrStat, ErrMsg )    ! Main interface routine.
   ! SUBROUTINE BEMT_GenSplines   ( FileNum, Re, Ctrl, InitInput, ErrStatLcl, ErrMsg )                                ! Phase #2: For BEM, generate splines for Cl and Cd at specific Re/Ctrl setting.
   ! SUBROUTINE BEMT_GetClCd      ( FileNum, AoA, Re, Ctrl, InitInput, Cl, Cd, ErrStatLcl, ErrMsg )                   ! Phase #3: For BEM, compute Cl and Cd for the given AoA.
   ! SUBROUTINE GetAllCoefs      ( FileNum, AoA, Re, Ctrl, InitInput, Cl, Cd, Cm, Cpmin, ErrStatLcl, ErrMsg )        ! Phase #4: Compute all requested airfoil coefficients for the given AoA.
   ! SUBROUTINE ReadAFfile       ( AFInfo, NumCoefs, Col_Cm, Col_Cpmin, ErrStat, ErrMsg )             ! Read an airfoil file.

   !=============================================================================
   SUBROUTINE AFI_Init ( InitInput, p, ErrStat, ErrMsg, UnEcho )


         ! This routine initializes AirfoilInfo by reading the airfoil files and generating the spline coefficients.


         ! Argument declarations.

      INTEGER(IntKi), INTENT(OUT)               :: ErrStat                    ! Error status.

      INTEGER, INTENT(IN), OPTIONAL             :: UnEcho                     ! I/O unit for echo file. If present and > 0, write to UnEcho.

      CHARACTER(*), INTENT(OUT)                 :: ErrMsg                     ! Error message.

      TYPE (AFI_InitInputType), INTENT(IN   )   :: InitInput                  ! This structure stores values that are set by the calling routine during the initialization phase.
      TYPE (AFI_ParameterType), INTENT(  OUT)   :: p                          ! This structure stores all the module parameters that are set by AirfoilInfo during the initialization phase.


         ! Local declarations.

      REAL(ReKi), ALLOCATABLE                :: Coef          (:)             ! The coefficients to send to the regrid routine for 2D splines.

      INTEGER                                :: Co                            ! The index into the coefficients array.
      INTEGER                                :: CoefSiz                       ! The size of the Coef array.
      INTEGER                                :: File                          ! The file index.
      INTEGER                                :: IA                            ! Index into the alpha array.
      INTEGER                                :: Indx                          ! Index into the arrays.
      INTEGER                                :: NumAl                         ! The length of the alpha array being sent to the regrid interface routine.
      INTEGER                                :: NumRe                         ! The length of the Re array being sent to the regrid interface routine.
      INTEGER                                :: Table                         ! Index into the Table array.
      INTEGER                                :: UnEc                          ! Local echo file unit number
      INTEGER                                :: NumCoefs                      ! The number of aerodynamic coefficients to be stored

      
      INTEGER                                :: ErrStat2                      ! Local error status.
      CHARACTER(ErrMsgLen)                   :: ErrMsg2
      CHARACTER(*), PARAMETER                :: RoutineName = 'AFI_Init'
      
      ErrStat = ErrID_None
      ErrMsg  = ""

         ! Display the version for this module.

      CALL DispNVD ( AFI_Ver )

      IF ( PRESENT(UnEcho) ) THEN
         UnEc = UnEcho
      ELSE
         UnEc = -1
      END IF
             
      
         ! Set the column indices for the various airfoil coefficients.

      p%ColCl    = 1  
      p%ColCd    = 2  
      p%ColCm    = 0 ! These may or may not be used; initialize to zero in case they aren't used
      p%ColCpmin = 0 ! These may or may not be used; initialize to zero in case they aren't used
      
      IF ( InitInput%InCol_Cm > 0 )  THEN
         p%ColCm = 3
         IF ( InitInput%InCol_Cpmin > 0 )  THEN
            p%ColCpmin = 4
         END IF
      ELSE IF ( InitInput%InCol_Cpmin > 0 )  THEN
            p%ColCpmin = 3
      END IF      
      NumCoefs = MAX(p%ColCd, p%ColCm,p%ColCpmin) ! number of non-zero coefficient columns
      
              
         ! Process the airfoil files.

      ALLOCATE ( p%AFInfo( InitInput%NumAFfiles ), STAT=ErrStat2 )
      IF ( ErrStat2 /= 0 )  THEN
         CALL SetErrStat ( ErrID_Fatal, 'Error allocating memory for the p%AFInfo array.', ErrStat, ErrMsg, RoutineName )
         RETURN
      ENDIF

      p%AFInfo( :)%ColCpmin=p%ColCpmin
      p%AFInfo( :)%ColCm=p%ColCm
    

      DO File=1,InitInput%NumAFfiles

         IF ( UnEc > 0 )  THEN
            WRITE (UnEc,'("--",/,A)')  'Contents of "'//TRIM( InitInput%FileNames(File) )//'":'
         END IF
         

         CALL ReadAFfile ( InitInput%FileNames(File), NumCoefs, InitInput%InCol_Alfa &
                         , InitInput%InCol_Cl, InitInput%InCol_Cd, InitInput%InCol_Cm, InitInput%InCol_Cpmin, p%AFInfo(File) &
                         , ErrStat2, ErrMsg2, UnEc ) 
            CALL SetErrStat ( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            IF ( ErrStat >= AbortErrLev )  THEN
               CALL Cleanup ( )
               RETURN
            ENDIF


            ! Make sure that all the tables meet the current restrictions.

         IF ( p%AFInfo(File)%NumTabs > 1 )  THEN


               ! Do the tables have the same number and set of alphas.

            DO Table=2,p%AFInfo(File)%NumTabs

               IF ( p%AFInfo(File)%Table(Table)%NumAlf /= p%AFInfo(File)%Table(1)%NumAlf )  THEN
          !     IF ( p%AFInfo(File)%NumTabs /= p%AFInfo(File)%Table(1)%NumAlf )  THEN
                  CALL SetErrStat ( ErrID_Fatal, 'Fatal Error: Airfoile file "'//TRIM( InitInput%FileNames(File) ) &
                                    //'", Table #'//TRIM( Num2LStr( Table ) ) &
                                    //' does not have the same set of alphas as the first table.', ErrStat, ErrMsg, RoutineName )
                  CALL Cleanup()
                  RETURN
               ENDIF

               DO IA=1,p%AFInfo(File)%Table(1)%NumAlf
                  IF ( p%AFInfo(File)%Table(Table)%Alpha(IA) /= p%AFInfo(File)%Table(1)%Alpha(IA) )  THEN
                     CALL SetErrStat ( ErrID_Fatal, 'Fatal Error: Airfoile file "'//TRIM( InitInput%FileNames(File) ) &
                                    //'", Table #'//TRIM( Num2LStr( Table ) ) &
                                    //' does not have the same set of alphas as the first table.', ErrStat, ErrMsg, RoutineName )
                     CALL Cleanup()
                     RETURN
                  ENDIF
               END DO ! IA

               IF ( p%AFInfo(File)%Table(Table)%Ctrl /= p%AFInfo(File)%Table(1)%Ctrl )  THEN
                  CALL SetErrStat ( ErrID_Fatal, 'Fatal Error: airfoile file "'//TRIM( InitInput%FileNames(File) ) &
                                 //'", Table #'//TRIM( Num2LStr( Table ) ) &
                                 //' does not have the same value for Ctrl as the first table.', ErrStat, ErrMsg, RoutineName )
                  CALL Cleanup()
                  RETURN
               ENDIF
            ENDDO ! Tab

         ENDIF ! ( p%AFInfo(File)%NumTabs > 1 )


            ! Create an array of Re for interpolation and fill it with the LOG(Re) for each table in the file.
!MLB: Maybe we should just store LogRe when we read in the Re.
         ALLOCATE ( p%AFInfo(File)%LogRe( p%AFInfo(File)%NumTabs ), STAT=ErrStat2 )
         IF ( ErrStat2 /= 0 )  THEN
            CALL SetErrStat ( ErrID_Fatal, 'Error allocating memory for the LogReAry array.', ErrStat, ErrMsg, RoutineName )
            CALL Cleanup()
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
            CALL AllocAry ( Coef, CoefSiz, "array for bicubic-spline coefficients for airfoil data", ErrStat2, ErrMsg2 )
               CALL SetErrStat ( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               ! Compute the splines for the drag coefficient.

            CALL AllocAry ( p%AFInfo(File)%CdSpCoef2D, CoefSiz, "array for bicubic-spline coefficients for Cd data" &
                          , ErrStat2, ErrMsg2 )
               CALL SetErrStat ( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               
            IF ( ErrStat >= AbortErrLev )  THEN
               CALL Cleanup ( )
               RETURN
            ENDIF

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
                                    , ErrStat2                  , ErrMsg2 )
               CALL SetErrStat ( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )


               ! Compute the splines for the lift coefficient.

            CALL AllocAry ( p%AFInfo(File)%ClSpCoef2D, CoefSiz, "array for bicubic-spline coefficients for Cl data" &
                          , ErrStat2, ErrMsg2 )
               CALL SetErrStat ( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

            IF ( ErrStat >= AbortErrLev )  THEN
               CALL Cleanup ( )
               RETURN
            ENDIF
               
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
                                    , ErrStat2                  , ErrMsg2 )
               CALL SetErrStat ( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )


               ! Optionally compute the splines for the pitching-moment coefficient.

            IF ( p%ColCm > 0 )  THEN

               CALL AllocAry ( p%AFInfo(File)%CmSpCoef2D, CoefSiz, "array for bicubic-spline coefficients for Cm data" &
                             , ErrStat2, ErrMsg2 )
               CALL SetErrStat ( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               IF ( ErrStat >= AbortErrLev )  THEN
                  CALL Cleanup ( )
                  RETURN
               ENDIF
               
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
                                       , ErrStat2                  , ErrMsg2 )
               CALL SetErrStat ( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

            ENDIF ! ( p%ColCm > 0 )


               ! Optionally compute the splines for the minimum pressure coefficient.

            IF ( p%ColCpmin > 0 )  THEN

               CALL AllocAry ( p%AFInfo(File)%CpminSpCoef2D, CoefSiz, "array for bicubic-spline coefficients for Cpmin data" &
                             , ErrStat2, ErrMsg2 )
               CALL SetErrStat ( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               IF ( ErrStat >= AbortErrLev )  THEN
                  CALL Cleanup ( )
                  RETURN
               ENDIF

               
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
                                       , ErrStat2                     , ErrMsg2 )
               CALL SetErrStat ( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

            ENDIF ! ( p%ColCpmin > 0 )


               ! Deallocate the temporary Coef array.

            IF ( ALLOCATED( Coef ) ) DEALLOCATE( Coef )


            ! Compute the spline coefficients of the piecewise cubic polynomials for the irregularly-spaced airfoil data in each file.
            ! Unless the data are constant.

            DO Co=1,NumCoefs

!BUG: What happened to the 2D interpolation?  Maybe I haven't written this yet.

               CALL SetErrStat ( ErrID_FATAL, 'The part to compute the bi-cubic splines in AFI_Init is not written yet!', ErrStat, ErrMsg, RoutineName )
               CALL Cleanup()
               RETURN

            END DO ! Co

         ELSE ! We have only one table.  Use 1-D interpolation.


               ! Allocate the arrays to hold spline coefficients.

            ALLOCATE ( p%AFInfo(File)%Table(1)%SplineCoefs( p%AFInfo(File)%Table(1)%NumAlf-1 &
                     , NumCoefs, 0:3 ), STAT=ErrStat2 )
            IF ( ErrStat2 /= 0 )  THEN
               CALL SetErrStat ( ErrStat2, 'Error allocating memory for the SplineCoefs array.', ErrStat, ErrMsg, RoutineName )
               CALL Cleanup()
               RETURN
            ENDIF

            
               ! Compute the one set of coefficients of the piecewise polynomials for the irregularly-spaced data.
               ! Unlike the 2-D interpolation in which we use diffent knots for each airfoil coefficient, we can do
               ! the 1-D stuff all at once.

            
            
            if ( p%AFInfo(File)%InterpOrd == 3_IntKi ) then
               
               ! bjj: what happens at the end points (these are periodic, so we should maybe extend the tables to make sure the end point?) 

                  ! use this for cubic splines:
               CALL CubicSplineInitM ( p%AFInfo(File)%Table(1)%Alpha &
                                     , p%AFInfo(File)%Table(1)%Coefs &
                                     , p%AFInfo(File)%Table(1)%SplineCoefs &
                                     , ErrStat2, ErrMsg2 )
               CALL SetErrStat ( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                           
            else if ( p%AFInfo(File)%InterpOrd == 1_IntKi ) then
               
                  ! use this for linear interpolation (sets the higher order coeffs to zero):
               
                  ! This is not the greatest way to get linear interpolation, but then we can use the same cubic spline routine
                  ! later without checking interp order there
               CALL CubicLinSplineInitM ( p%AFInfo(File)%Table(1)%Alpha &
                                     , p%AFInfo(File)%Table(1)%Coefs &
                                     , p%AFInfo(File)%Table(1)%SplineCoefs &
                                     , ErrStat2, ErrMsg2 )
               CALL SetErrStat ( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               
            else
               
               CALL SetErrStat ( ErrID_FATAL, 'Airfoil file "'//TRIM( InitInput%FileNames(File) ) &
                                         //'": InterpOrd must be 1 (linear) or 3 (cubic spline).', ErrStat, ErrMsg, RoutineName )
               CALL Cleanup()
               RETURN
               
            end if
            
                              
         ENDIF ! ( p%AFInfo(File)%NumTabs > 1 )


            ! Compute the spline coefficients of the piecewise cubic polynomials for the irregularly-spaced airfoil data in each file.
            ! Unless the data are constant.

         DO Co=1,NumCoefs


            IF ( p%AFInfo(File)%NumTabs > 1 )  THEN                  ! We use 2D cubic spline interpolation.

!BUG: What happened to the 2D interpolation?  Maybe I haven't written this yet.

               CALL SetErrStat ( ErrID_FATAL, 'The part to compute the bi-cubic splines in AFI_Init is not written yet!', ErrStat, ErrMsg, RoutineName )
               CALL Cleanup()
               RETURN

! We also need to deal with constant data.

            ELSE                                                     ! We use 1D cubic spline interpolation if the data are not constant.

               IF ( p%AFInfo(File)%Table(1)%ConstData )  THEN

                     ! We need to deal with constant data.

                  CALL SetErrStat ( ErrID_FATAL, 'The part to deal with constant data in AFI_Init is not written yet!', ErrStat, ErrMsg, RoutineName )
                  CALL Cleanup()
                  RETURN

               ENDIF ! p%AFInfo(File)%Table(1)%ConstData

            ENDIF ! ( p%AFInfo(File)%NumTabs > 1 )

         END DO ! Co

      END DO ! File

      CALL Cleanup ( )

      RETURN

   !=======================================================================
   CONTAINS
   !=======================================================================
      SUBROUTINE Cleanup ( )

         ! This subroutine cleans up the parent routine before exiting.
            ! Deallocate the temporary Coef array.

         IF ( ALLOCATED( Coef ) ) DEALLOCATE( Coef )

         RETURN

      END SUBROUTINE Cleanup 

   END SUBROUTINE AFI_Init
   
  
   !=============================================================================
   SUBROUTINE ReadAFfile ( AFfile, NumCoefs, InCol_Alfa, InCol_Cl, InCol_Cd, InCol_Cm, InCol_Cpmin, AFInfo &
                         , ErrStat, ErrMsg, UnEc )


         ! This routine reads an airfoil file.


         ! Argument declarations.

      INTEGER(IntKi),    INTENT(IN)           :: InCol_Alfa                    ! The airfoil-table input column for angle of attack.
      INTEGER(IntKi),    INTENT(IN)           :: InCol_Cd                      ! The airfoil-table input column for drag coefficient.
      INTEGER(IntKi),    INTENT(IN)           :: InCol_Cl                      ! The airfoil-table input column for lift coefficient.
      INTEGER(IntKi),    INTENT(IN)           :: InCol_Cm                      ! The airfoil-table input column for pitching-moment coefficient.
      INTEGER(IntKi),    INTENT(IN)           :: InCol_Cpmin                   ! The airfoil-table input column for minimum pressure coefficient.
      INTEGER(IntKi),    INTENT(  OUT)        :: ErrStat                       ! Error status.
      INTEGER(IntKi),    INTENT(IN)           :: NumCoefs                      ! The number of aerodynamic coefficients to be stored.

      INTEGER,           INTENT(IN)           :: UnEc                          ! I/O unit for echo file. If present and > 0, write to UnEc.      CHARACTER(*), INTENT(IN)               :: AFfile                        ! The file to be read.

      CHARACTER(*),      INTENT(IN)           :: AFfile                        ! The file to read in.
      CHARACTER(*),      INTENT(  OUT)        :: ErrMsg                        ! Error message.

      TYPE (AFInfoType), INTENT(INOUT)        :: AFInfo                        ! The derived type for holding the constant parameters for this airfoil.


         ! Local declarations.

      REAL(ReKi)                              :: Coords   (2)                  ! An array to hold data from the airfoil-shape table.
      REAL(ReKi), ALLOCATABLE                 :: SiAry    (:)                  ! A temporary array to hold data from a table.
                                              
      INTEGER                                 :: Coef                          ! A DO index that points into the coefficient array.
      INTEGER                                 :: Cols2Parse                    ! The number of columns that must be read from the coefficient tables.
      INTEGER                                 :: CurLine                       ! The current line to be parsed in the FileInfo structure.
      INTEGER                                 :: NumAlf                        ! The number of rows in the current table.
      INTEGER                                 :: Row                           ! The row of a table to be parsed in the FileInfo structure.
      INTEGER                                 :: Table                         ! The DO index for the tables.
                                              
      LOGICAL                                 :: BadVals                       ! A flag that indicates if the values in a table are invalid.

      TYPE (FileInfoType)                     :: FileInfo                      ! The derived type for holding the file information.

      INTEGER(IntKi)                          :: DefaultInterpOrd              ! value of default interp order
      INTEGER(IntKi)                          :: ErrStat2                      ! Error status local to this routine.
      CHARACTER(ErrMsgLen)                    :: ErrMsg2
      CHARACTER(*), PARAMETER                 :: RoutineName = 'ReadAFfile'
      CHARACTER(10)                           :: defaultStr
      
      ErrStat = ErrID_None
      ErrMsg  = ""
      defaultStr = ""
      
         ! Process the (possibly) nested set of files.  This copies the decommented contents of
         ! AFI_FileInfo%FileName and the files it includes (both directly and indirectly) into
         ! the FileInfo structure that we can then parse.

      CALL ProcessComFile ( AFfile, FileInfo, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF (ErrStat >= AbortErrLev) THEN
            CALL Cleanup()
            RETURN
         END IF
         

         ! Process the airfoil shape information if it is included.

      CurLine = 1

   DefaultInterpOrd = 3      
#ifdef LINEAR_INTERP
   DefaultInterpOrd = 1      
#endif
      
      CALL ParseVarWDefault ( FileInfo, CurLine, 'InterpOrd', AFInfo%InterpOrd, DefaultInterpOrd, ErrStat2, ErrMsg2, UnEc )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

         ! These first two parameters are optional, so we don't check for errors. ! bjj: huh???? 
      
      CALL ParseVar ( FileInfo, CurLine, 'NonDimArea', AFInfo%NonDimArea, ErrStat2, ErrMsg2, UnEc )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL ParseVar ( FileInfo, CurLine, 'NumCoords' , AFInfo%NumCoords , ErrStat2, ErrMsg2, UnEc )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF (ErrStat >= AbortErrLev) THEN
            CALL Cleanup()
            RETURN
         END IF

      IF ( AFInfo%NumCoords > 0 )  THEN

         ALLOCATE ( AFInfo%X_Coord( AFInfo%NumCoords ) , STAT=ErrStat2 )
         IF ( ErrStat2 /= 0 )  THEN
            CALL SetErrStat ( ErrID_Fatal, 'Error allocating memory for AFInfo%X_Coord.', ErrStat, ErrMsg, RoutineName )
            CALL Cleanup()
            RETURN
         ENDIF

         ALLOCATE ( AFInfo%Y_Coord( AFInfo%NumCoords ) , STAT=ErrStat2 )
         IF ( ErrStat2 /= 0 )  THEN
            CALL SetErrStat ( ErrID_Fatal, 'Error allocating memory for AFInfo%Y_Coord.', ErrStat, ErrMsg, RoutineName )
            CALL Cleanup()
            RETURN
         ENDIF

         DO Row=1,AFInfo%NumCoords
            CALL ParseAry ( FileInfo, CurLine, 'X_Coord/Y_Coord', Coords, 2, ErrStat2, ErrMsg2, UnEc )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               IF (ErrStat >= AbortErrLev) THEN
                  CALL Cleanup()
                  RETURN
               END IF
            AFInfo%X_Coord(Row) = Coords(1)
            AFInfo%Y_Coord(Row) = Coords(2)
         ENDDO ! Row

      ENDIF


         ! How many columns do we need to read in the input and how many total coefficients will be used?

      Cols2Parse = MAX( InCol_Alfa, InCol_Cl, InCol_Cd, InCol_Cm, InCol_Cpmin )
      ALLOCATE ( SiAry( Cols2Parse ) , STAT=ErrStat2 )
      IF ( ErrStat2 /= 0 )  THEN
         CALL SetErrStat ( ErrID_Fatal, 'Error allocating memory for SiAry.', ErrStat, ErrMsg, RoutineName )
         CALL Cleanup()
         RETURN
      ENDIF


         ! Work through the multiple tables.

      CALL ParseVar ( FileInfo, CurLine, 'NumTabs' , AFInfo%NumTabs , ErrStat2, ErrMsg2, UnEc )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF (ErrStat >= AbortErrLev) THEN
            CALL Cleanup()
            RETURN
         END IF

      IF ( AFInfo%NumTabs < 1 )  THEN
         CALL SetErrStat ( ErrID_Fatal, 'NumTabs must be > 0.', ErrStat, ErrMsg, RoutineName )
         CALL Cleanup()
         RETURN
      ENDIF

      ALLOCATE ( AFInfo%Table( AFInfo%NumTabs ) , STAT=ErrStat2 )
      IF ( ErrStat2 /= 0 )  THEN
         CALL SetErrStat ( ErrID_Fatal, 'Error allocating memory for AFInfo%Table.', ErrStat, ErrMsg, RoutineName )
         CALL Cleanup()
         RETURN
      ENDIF

      DO Table=1,AFInfo%NumTabs

         CALL ParseVar ( FileInfo, CurLine, 'Re', AFInfo%Table(Table)%Re, ErrStat2, ErrMsg2, UnEc )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            IF (ErrStat >= AbortErrLev) THEN
               CALL Cleanup()
               RETURN
            END IF
         IF ( AFInfo%Table(Table)%Re <= 0.0 )  THEN
               CALL SetErrStat ( ErrID_Severe, 'Re must be > 0 in "'//TRIM( AFfile ) &
                                      //'".'//NewLine//'  >> The error occurred on line #' &
                                      //TRIM( Num2LStr( FileInfo%FileLine(CurLine-1) ) )//'.', ErrStat, ErrMsg, RoutineName )
               CALL Cleanup()
               RETURN
         ENDIF ! ( AFInfo%Table(Table)%Re <= 0.0 )

         CALL ParseVar ( FileInfo, CurLine, 'Ctrl', AFInfo%Table(Table)%Ctrl, ErrStat2, ErrMsg2, UnEc )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            IF (ErrStat >= AbortErrLev) THEN
               CALL Cleanup()
               RETURN
            END IF
         IF ( AFInfo%Table(Table)%Ctrl /= 0.0 )  THEN !bjj: use EqualRealNos?
               CALL SetErrStat ( ErrID_Severe, 'Ctrl must equal 0.0 in "'//TRIM( AFfile ) &
                                       //'".'//NewLine//'  >> The error occurred on line #' &
                                       //TRIM( Num2LStr( FileInfo%FileLine(CurLine-1) ) )//'.', ErrStat, ErrMsg, RoutineName )
               CALL Cleanup()
               RETURN
         ENDIF ! ( AFInfo%Table(Table)%Ctrl <= 0.0 )

         CALL ParseVar ( FileInfo, CurLine, 'InclUAdata', AFInfo%Table(Table)%InclUAdata, ErrStat2, ErrMsg2, UnEc )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            IF (ErrStat >= AbortErrLev) THEN
               CALL Cleanup()
               RETURN
            END IF

         IF ( AFInfo%Table(Table)%InclUAdata )  THEN

               CALL ParseVar ( FileInfo, CurLine, 'alpha0', AFInfo%Table(Table)%UA_BL%alpha0, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVar ( FileInfo, CurLine, 'alpha1', AFInfo%Table(Table)%UA_BL%alpha1, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVar ( FileInfo, CurLine, 'alpha2', AFInfo%Table(Table)%UA_BL%alpha2, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVar ( FileInfo, CurLine, 'eta_e', AFInfo%Table(Table)%UA_BL%eta_e, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVar ( FileInfo, CurLine, 'C_nalpha', AFInfo%Table(Table)%UA_BL%C_nalpha, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            
               CALL ParseVarWDefault ( FileInfo, CurLine, 'T_f0', AFInfo%Table(Table)%UA_BL%T_f0, 3.0_ReKi, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         
               CALL ParseVarWDefault ( FileInfo, CurLine, 'T_V0', AFInfo%Table(Table)%UA_BL%T_V0, 6.0_ReKi, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               
               CALL ParseVarWDefault ( FileInfo, CurLine, 'T_p', AFInfo%Table(Table)%UA_BL%T_p, 1.7_ReKi, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVarWDefault ( FileInfo, CurLine, 'T_VL', AFInfo%Table(Table)%UA_BL%T_VL, 11.0_ReKi, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVarWDefault ( FileInfo, CurLine, 'b1', AFInfo%Table(Table)%UA_BL%b1, .14_ReKi, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVarWDefault ( FileInfo, CurLine, 'b2', AFInfo%Table(Table)%UA_BL%b2, .53_ReKi, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVarWDefault ( FileInfo, CurLine, 'b5', AFInfo%Table(Table)%UA_BL%b5, 5.0_ReKi, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVarWDefault ( FileInfo, CurLine, 'A1', AFInfo%Table(Table)%UA_BL%A1, .3_ReKi, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVarWDefault ( FileInfo, CurLine, 'A2', AFInfo%Table(Table)%UA_BL%A2, .7_ReKi, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVarWDefault ( FileInfo, CurLine, 'A5', AFInfo%Table(Table)%UA_BL%A5, 1.0_ReKi, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVar ( FileInfo, CurLine, 'S1', AFInfo%Table(Table)%UA_BL%S1, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVar ( FileInfo, CurLine, 'S2', AFInfo%Table(Table)%UA_BL%S2, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVar ( FileInfo, CurLine, 'S3', AFInfo%Table(Table)%UA_BL%S3, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVar ( FileInfo, CurLine, 'S4', AFInfo%Table(Table)%UA_BL%S4, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVar ( FileInfo, CurLine, 'Cn1', AFInfo%Table(Table)%UA_BL%Cn1, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVar ( FileInfo, CurLine, 'Cn2', AFInfo%Table(Table)%UA_BL%Cn2, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVarWDefault ( FileInfo, CurLine, 'St_sh', AFInfo%Table(Table)%UA_BL%St_sh, .19_ReKi, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVar ( FileInfo, CurLine, 'Cd0', AFInfo%Table(Table)%UA_BL%Cd0, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVar ( FileInfo, CurLine, 'Cm0', AFInfo%Table(Table)%UA_BL%Cm0, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVar ( FileInfo, CurLine, 'k0', AFInfo%Table(Table)%UA_BL%k0, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVar ( FileInfo, CurLine, 'k1', AFInfo%Table(Table)%UA_BL%k1, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVar ( FileInfo, CurLine, 'k2', AFInfo%Table(Table)%UA_BL%k2, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVar ( FileInfo, CurLine, 'k3', AFInfo%Table(Table)%UA_BL%k3, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVar ( FileInfo, CurLine, 'k1_hat', AFInfo%Table(Table)%UA_BL%k1_hat, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVarWDefault ( FileInfo, CurLine, 'x_cp_bar', AFInfo%Table(Table)%UA_BL%x_cp_bar, .2_ReKi, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                  
               CALL ParseVarWDefault ( FileInfo, CurLine, 'UACutout', AFInfo%Table(Table)%UA_BL%UACutout, 45.0_ReKi, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               AFInfo%Table(Table)%UA_BL%UACutout = AFInfo%Table(Table)%UA_BL%UACutout*D2R

               CALL ParseVarWDefault ( FileInfo, CurLine, 'filtCutOff', AFInfo%Table(Table)%UA_BL%filtCutOff, 20.0_ReKi, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                  
               IF (ErrStat >= AbortErrLev) THEN
                  CALL Cleanup()
                  RETURN
               END IF
                  
         ELSE !ELSEIF ( UA_Model == 1 )  THEN

            !CALL SetErrStat ( ErrID_Fatal &
            !, 'You must supply Beddoes-Leishman unsteady aerodynamics parameters for all airfoils if you want to use that' &
            !//' model. You did not do so for Table #'//TRIM( Num2LStr( Table ) )//' in the "'//TRIM( AFfile )//'" airfoil file.', ErrStat, ErrMsg, RoutineName )
            !CALL Cleanup()
            !RETURN

         ENDIF ! ( AFInfo%Table(Table)%InclUAdata )

         CALL ParseVar ( FileInfo, CurLine, 'NumAlf', AFInfo%Table(Table)%NumAlf, ErrStat2, ErrMsg2, UnEc )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            IF (ErrStat >= AbortErrLev) THEN
               CALL Cleanup()
               RETURN
            END IF

         IF ( AFInfo%Table(Table)%NumAlf < 1 )  THEN
            CALL SetErrStat( ErrID_Fatal, 'NumAlf must be a positive number on line #' &
                           //TRIM( Num2LStr( FileInfo%FileLine(CurLine-1) ) )//' in "'//TRIM( AFfile )//'".', ErrStat, ErrMsg, RoutineName )
            CALL Cleanup()
            RETURN
         ELSEIF ( AFInfo%Table(Table)%NumAlf < 3 )  THEN
            AFInfo%Table(Table)%ConstData = .TRUE.
         ELSE
            AFInfo%Table(Table)%ConstData = .FALSE.
         ENDIF ! ( Test for valid values for NumAlf )


            ! Allocate the arrays for the airfoil coefficients.

         ALLOCATE ( AFInfo%Table(Table)%Alpha( AFInfo%Table(Table)%NumAlf ), STAT=ErrStat2 )
         IF ( ErrStat2 /= 0 )  THEN
            CALL SetErrStat( ErrID_Fatal, 'Error allocating memory for AFInfo%Table%Alpha.', ErrStat, ErrMsg, RoutineName )
            CALL Cleanup()
            RETURN
         ENDIF

         ALLOCATE ( AFInfo%Table(Table)%Coefs( AFInfo%Table(Table)%NumAlf, NumCoefs ), STAT=ErrStat2 )
         IF ( ErrStat2 /= 0 )  THEN
            CALL SetErrStat( ErrID_Fatal, 'Error allocating memory for AFInfo%Table%Coefs.', ErrStat, ErrMsg, RoutineName )
            CALL Cleanup()
            RETURN
         ENDIF

         DO Row=1,AFInfo%Table(Table)%NumAlf

            CALL ParseAry ( FileInfo, CurLine, 'CoeffData', SiAry, Cols2Parse, ErrStat2, ErrMsg2, UnEc )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               IF (ErrStat >= AbortErrLev) THEN
                  CALL Cleanup()
                  RETURN
               END IF

            AFInfo%Table(Table)%Alpha(Row  ) = SiAry(InCol_Alfa)*D2R
           !AFInfo%Table(Table)%Alpha(Row  ) = SiAry(InCol_Alfa)
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
            IF ( .NOT. EqualRealNos( AFInfo%Table(Table)%Alpha(1), -Pi ) )  THEN
            !IF ( .NOT. EqualRealNos( AFInfo%Table(Table)%Alpha(1), -180.0_ReKi ) )  THEN
               BadVals = .TRUE.
            ENDIF
            IF ( .NOT. EqualRealNos( AFInfo%Table(Table)%Alpha(NumAlf), Pi ) )  THEN
            !IF ( .NOT. EqualRealNos( AFInfo%Table(Table)%Alpha(NumAlf), 180.0_ReKi ) )  THEN
               BadVals = .TRUE.
            ENDIF
            DO Coef=1,NumCoefs
               IF ( .NOT. EqualRealNos( AFInfo%Table(Table)%Coefs(1,Coef), AFInfo%Table(Table)%Coefs(NumAlf,Coef) ) )  THEN
                  BadVals = .TRUE.
               ENDIF
            ENDDO ! Coef
            IF ( BadVals )  THEN
               CALL SetErrStat( ErrID_Warn, &
                  'Airfoil data should go from -180 degrees to 180 degrees and the coefficients at the ends should be the same.', ErrStat, ErrMsg, RoutineName )
               CALL Cleanup()
               RETURN
            ENDIF
         ENDIF ! ( .NOT. AFInfo%Table(Table)%ConstData )

      ENDDO ! Table


      CALL Cleanup( )

      RETURN

      !=======================================================================
      CONTAINS
      !=======================================================================
                  
         SUBROUTINE Cleanup ()

            ! This subroutine cleans up all the allocatable arrays, sets the error status/message and closes the binary file
               
            CALL NWTC_Library_DestroyFileInfoType (FileInfo, ErrStat2, ErrMsg2)
            IF ( ALLOCATED(SiAry) ) DEALLOCATE(SiAry)
            
         END SUBROUTINE Cleanup 


   END SUBROUTINE ReadAFfile
      
      
   subroutine AFI_GetAirfoilParams( AFInfo, M, Re, BL_p, C_nalpha_circ, errMsg, errStat )     

   type(AFInfoType),     intent(in   )       :: AFInfo                        ! The derived type for holding the constant parameters for this airfoil.
   real(ReKi),           intent(in   )       :: M                             ! mach number
   real(ReKi),           intent(in   )       :: Re                            ! Reynold's number
   type(AFI_UA_BL_Type), intent(  out)       :: BL_p                          ! airfoil constants (UA BL parameters)

   !real(ReKi),       intent(  out)       :: alpha0                        ! zero lift angle of attack (radians)
   !real(ReKi),       intent(  out)       :: alpha1                        ! angle of attack at f = 0.7, approximately the stall angle; for alpha >= alpha0 (radians)
   !real(ReKi),       intent(  out)       :: alpha2                        ! angle of attack at f = 0.7, approximately the stall angle; for alpha < alpha0 (radians)
   !real(ReKi),       intent(  out)       :: eta_e                         !
   !real(ReKi),       intent(  out)       :: C_nalpha                      !
   real(ReKi),       intent(  out)       :: C_nalpha_circ                 ! slope of the circulatory normal force coefficient vs alpha curve
   !real(ReKi),       intent(  out)       :: T_f0                          ! initial value of T_f, airfoil specific, used to compute D_f and fprimeprime
   !real(ReKi),       intent(  out)       :: T_V0                          ! initial value of T_V, airfoil specific, time parameter associated with the vortex lift decay process, used in Cn_v
   !real(ReKi),       intent(  out)       :: T_p                           ! boundary-layer, leading edge pressure gradient time parameter; used in D_p; airfoil specific
   !real(ReKi),       intent(  out)       :: T_VL                          ! time variable associated with the vortex advection process; it represents the non-dimensional time in semi-chords needed for a vortex to travel from leading edge to trailing edge
   !real(ReKi),       intent(  out)       :: St_sh                         !
   !real(ReKi),       intent(  out)       :: b1                            ! airfoil constant derived from experimental results, usually 0.14
   !real(ReKi),       intent(  out)       :: b2                            ! airfoil constant derived from experimental results, usually 0.53
   !real(ReKi),       intent(  out)       :: b5                            ! airfoil constant derived from experimental results, usually 5.0
   !real(ReKi),       intent(  out)       :: A1                            ! airfoil constant derived from experimental results, usually 0.3 
   !real(ReKi),       intent(  out)       :: A2                            ! airfoil constant derived from experimental results, usually 0.7
   !real(ReKi),       intent(  out)       :: A5                            ! airfoil constant derived from experimental results, usually 1.0
   !real(ReKi),       intent(  out)       :: S1                            ! constant in the f-curve best-fit, alpha >= alpha0 
   !real(ReKi),       intent(  out)       :: S2                            ! constant in the f-curve best-fit, alpha >= alpha0 
   !real(ReKi),       intent(  out)       :: S3                            ! constant in the f-curve best-fit, alpha <  alpha0 
   !real(ReKi),       intent(  out)       :: S4                            ! constant in the f-curve best-fit, alpha <  alpha0 
   !real(ReKi),       intent(  out)       :: Cn1                           ! critical value of Cn_prime at LE separation for alpha >= alpha0
   !real(ReKi),       intent(  out)       :: Cn2                           ! critical value of Cn_prime at LE separation for alpha < alpha0
   !real(ReKi),       intent(  out)       :: Cd0                           !
   !real(ReKi),       intent(  out)       :: Cm0                           ! 2D pitching moment coefficient at zero lift, positive if nose is up
   !real(ReKi),       intent(  out)       :: k0                            ! airfoil parameter in the x_cp_hat curve best-fit
   !real(ReKi),       intent(  out)       :: k1                            ! airfoil parameter in the x_cp_hat curve best-fit
   !real(ReKi),       intent(  out)       :: k2                            ! airfoil parameter in the x_cp_hat curve best-fit
   !real(ReKi),       intent(  out)       :: k3                            ! airfoil parameter in the x_cp_hat curve best-fit
   !real(ReKi),       intent(  out)       :: k1_hat                        !
   !real(ReKi),       intent(  out)       :: x_cp_bar                      ! airfoil parameter for calulating x_cp_v
   !real(ReKi),       intent(  out)       :: filtCutOff                    ! airfoil parameter for the low-pass cut-off frequency for pitching rate and accelerations (Hz)
   integer(IntKi),   intent(  out)       :: errStat                       ! Error status. 
   character(*),     intent(  out)       :: errMsg                        ! Error message.
      
      errMsg         = ''
      errStat        = ErrID_None
           
      
      ! These coefs are stored in the AFInfo data structures based on Re
   BL_p%alpha0         =  AFInfo%Table(1)%UA_BL%alpha0   * D2R   ! Convert to radians
   BL_p%alpha1         =  AFInfo%Table(1)%UA_BL%alpha1   * D2R   ! Convert to radians
   BL_p%alpha2         =  AFInfo%Table(1)%UA_BL%alpha2   * D2R   ! Convert to radians
   BL_p%eta_e          =  AFInfo%Table(1)%UA_BL%eta_e         !0.90               ! Recovery factor in the range [0.85 - 0.95]
   BL_p%C_nalpha       =  AFInfo%Table(1)%UA_BL%C_nalpha      !2*pi  
   BL_p%T_f0           =  AFInfo%Table(1)%UA_BL%T_f0          !3.0_ReKi  ! seconds
   BL_p%T_V0           =  AFInfo%Table(1)%UA_BL%T_V0          !6.0_ReKi
   BL_p%T_p            =  AFInfo%Table(1)%UA_BL%T_p           !1.7_ReKi
   BL_p%T_VL           =  AFInfo%Table(1)%UA_BL%T_VL          !11.0_ReKi
   BL_p%b1             =  AFInfo%Table(1)%UA_BL%b1            !0.14_ReKi
   BL_p%b2             =  AFInfo%Table(1)%UA_BL%b2            !0.53_ReKi
   BL_p%b5             =  AFInfo%Table(1)%UA_BL%b5            !5.0_ReKi
   BL_p%A1             =  AFInfo%Table(1)%UA_BL%A1            !0.3_ReKi
   BL_p%A2             =  AFInfo%Table(1)%UA_BL%A2            !0.70_ReKi
   BL_p%A5             =  AFInfo%Table(1)%UA_BL%A5            !1.0_ReKi
   BL_p%S1             =  AFInfo%Table(1)%UA_BL%S1            !0.0262_ReKi   !!!!!!!!!!
   BL_p%S2             =  AFInfo%Table(1)%UA_BL%S2            !0.0201_ReKi   !!!!!!!!!!
   BL_p%S3             =  AFInfo%Table(1)%UA_BL%S3            !0.0262_ReKi   !!!!!!!!!!
   BL_p%S4             =  AFInfo%Table(1)%UA_BL%S4            !0.0201_ReKi   !!!!!!!!!!
   BL_p%Cn1            =  AFInfo%Table(1)%UA_BL%Cn1           !1.264_ReKi  ! Stall values of Cn
   BL_p%Cn2            =  AFInfo%Table(1)%UA_BL%Cn2           !-0.833_ReKi
   BL_p%St_sh          =  AFInfo%Table(1)%UA_BL%St_sh         !0.19_ReKi
   BL_p%Cd0            =  AFInfo%Table(1)%UA_BL%Cd0           !0.012_ReKi
   BL_p%Cm0            =  AFInfo%Table(1)%UA_BL%Cm0           !0.0_ReKi
   BL_p%k0             =  AFInfo%Table(1)%UA_BL%k0            !0.0_ReKi
   BL_p%k1             =  AFInfo%Table(1)%UA_BL%k1            !0.0_ReKi
   BL_p%k2             =  AFInfo%Table(1)%UA_BL%k2            !0.0_ReKi
   BL_p%k3             =  AFInfo%Table(1)%UA_BL%k3            !0.0_ReKi
   BL_p%k1_hat         =  AFInfo%Table(1)%UA_BL%k1_hat        !0.0_ReKi
   BL_p%x_cp_bar       =  AFInfo%Table(1)%UA_BL%x_cp_bar      !0.2_ReKi
   BL_p%filtCutOff     =  AFInfo%Table(1)%UA_BL%filtCutOff    ! 5.0_ReKi  Hz
   
   C_nalpha_circ  =  BL_p%C_nalpha / sqrt(1.0_ReKi-M**2)
     ! Cn1=1.9 Tp=1.7 Tf=3., Tv=6 Tvl=11, Cd0=0.012
   
   end subroutine AFI_GetAirfoilParams
                                 
!=============================================================================
END MODULE AirfoilInfo
