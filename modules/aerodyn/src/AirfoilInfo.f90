!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2015-2018  National Renewable Energy Laboratory
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
   USE                                          :: ISO_FORTRAN_ENV , ONLY : IOSTAT_EOR
   USE                                          :: NWTC_LAPACK

   IMPLICIT NONE

   PRIVATE

   PUBLIC                                       :: AFI_Init ! routine to initialize AirfoilInfo parameters
   PUBLIC                                       :: AFI_ComputeUACoefs        ! routine to calculate Airfoil BL parameters for UA
   PUBLIC                                       :: AFI_ComputeAirfoilCoefs   ! routine to perform 1D (AOA) or 2D (AOA, Re) or (AOA, UserProp) lookup of the airfoil coefs
   PUBLIC                                       :: AFI_WrHeader
   PUBLIC                                       :: AFI_WrData
   PUBLIC                                       :: AFI_WrTables

   TYPE(ProgDesc), PARAMETER                    :: AFI_Ver = ProgDesc( 'AirfoilInfo', '', '')    ! The name, version, and date of AirfoilInfo.

   integer, parameter                           :: MaxNumAFCoeffs = 7 !cl,cd,cm,cpMin, UA:f_st, FullySeparate, FullyAttached

CONTAINS


   function CheckValuesAreUniqueMonotonicIncreasing(secondVals)
     
      real(ReKi),  intent(in   )  :: secondVals(:)
      logical CheckValuesAreUniqueMonotonicIncreasing
      
      
      integer(IntKi) :: i
      
      CheckValuesAreUniqueMonotonicIncreasing = .true.
      
      do i = 2, size(secondVals)
         if ( EqualRealNos(secondVals(i), secondVals(i-1)) .or. (secondVals(i) < secondVals(i-1))) then
            CheckValuesAreUniqueMonotonicIncreasing = .false.
            exit
         end if
      end do
      
       
   end function CheckValuesAreUniqueMonotonicIncreasing
   
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

      INTEGER                                :: UnEc                          ! Local echo file unit number
      INTEGER                                :: NumCoefs                      ! The number of aerodynamic coefficients to be stored
      integer                                :: iTable                        ! Iterator for airfoil tables
      
      
      INTEGER                                :: ErrStat2                      ! Local error status.
      CHARACTER(ErrMsgLen)                   :: ErrMsg2
      CHARACTER(*), PARAMETER                :: RoutineName = 'AFI_Init'

      ErrStat = ErrID_None
      ErrMsg  = ""

         ! Display the version for this module.

      !CALL DispNVD ( AFI_Ver )
      p%FileName = InitInput%FileName ! store this for error messages later (e.g., in UA)

      
      CALL AFI_ValidateInitInput(InitInput, ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat >= AbortErrLev) return
      
      IF ( PRESENT(UnEcho) ) THEN
         UnEc = UnEcho
      ELSE
         UnEc = -1
      END IF
             
         ! Set the lookup model:  1 = 1D, 2 = 2D based on (AoA,Re), 3 = 2D based on (AoA,UserProp)
      p%AFTabMod   = InitInput%AFTabMod
      
         ! Set the column indices for the various airfoil coefficients as they will be stored in our data structures, 
         !   NOT as they are recorded in the airfoil input file (InitInput%InCol_*) !
      p%ColCl    = 1  
      p%ColCd    = 2  
      p%ColCm    = 0 ! These may or may not be used; initialize to zero in case they aren't used
      p%ColCpmin = 0 ! These may or may not be used; initialize to zero in case they aren't used
      p%ColUAf   = 0 ! These may or may not be used; initialize to zero in case they aren't used
      
      IF ( InitInput%InCol_Cm > 0 )  THEN
         p%ColCm = 3
         IF ( InitInput%InCol_Cpmin > 0 )  THEN
            p%ColCpmin = 4
         END IF
      ELSE IF ( InitInput%InCol_Cpmin > 0 )  THEN
            p%ColCpmin = 3
      END IF      
      NumCoefs = MAX(p%ColCd, p%ColCm,p%ColCpmin) ! number of non-zero coefficient columns
      
              
      ! Process the airfoil file.
      

      IF ( UnEc > 0 )  THEN
         WRITE (UnEc,'("--",/,A)')  'Contents of "'//TRIM( InitInput%FileName )//'":'
      END IF
         

      CALL ReadAFfile ( InitInput, NumCoefs, p, ErrStat2, ErrMsg2, UnEc ) 
         CALL SetErrStat ( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF ( ErrStat >= AbortErrLev )  THEN
            CALL Cleanup ( )
            RETURN
         ENDIF


         ! Make sure that all the tables meet the current restrictions.
         
      IF ( p%NumTabs > 1 )  THEN

         IF ( p%AFTabMod == AFITable_1 ) THEN
               ! 1D interpolation was specified even though multiple tables exist in the file
            
            p%NumTabs = 1
            CALL SetErrStat ( ErrID_Warn, 'DimModel = 1D, therefore using only the first airfoil table in the file: "'//TRIM( InitInput%FileName ), ErrStat, ErrMsg, RoutineName )
         
         ELSE !IF ( p%AFTabMod /= AFITable_1 ) THEN
               
            !--------------
            ! Check that secondary values are unique and monotonically increasing:
            !--------------
            ALLOCATE(p%secondVals(p%NumTabs), STAT=ErrStat2 )
               IF ( ErrStat2 /= 0 )  THEN
                  CALL SetErrStat ( ErrID_Fatal, 'Error allocating memory for the secondVals array.', ErrStat, ErrMsg, RoutineName )
                  RETURN
               ENDIF
                  
               
            IF (p%AFTabMod == AFITable_2Re) THEN
               ! Ctrl Value must be the same, Re must be monotonic increasing without repeats
               
               DO iTable=2,p%NumTabs
                     
                  IF ( p%Table(iTable)%UserProp /= p%Table(1)%UserProp )  THEN
                     CALL SetErrStat ( ErrID_Fatal, 'Fatal Error: airfoil file "'//TRIM( InitInput%FileName ) &
                                    //'", Table #'//TRIM( Num2LStr( iTable ) ) &
                                    //' does not have the same value for Ctrl Property (UserProp) as the first table.', ErrStat, ErrMsg, RoutineName )
                     CALL Cleanup()
                     RETURN
                  ENDIF
                  
               END DO ! iTable

               DO iTable=1,p%NumTabs
                  if (p%Table(iTable)%Re < 0.0_ReKi) then
                     CALL SetErrStat ( ErrID_Fatal, 'Fatal Error: airfoil file "'//TRIM( InitInput%FileName ) &
                                    //'", Table #'//TRIM( Num2LStr( iTable ) ) &
                                    //' has a negative Reynolds Number.', ErrStat, ErrMsg, RoutineName )
                     CALL Cleanup()
                     RETURN
                  end if
                  
                  p%Table(iTable)%Re   = max( p%Table(iTable)%Re, 0.001_ReKi )
                  
#ifndef AFI_USE_LINEAR_RE
                  p%secondVals(iTable) = log( p%Table(iTable)%Re )
#else
                  p%secondVals(iTable) =      p%Table(iTable)%Re
#endif
               END DO ! iTable

            ELSE IF (p%AFTabMod == AFITable_2User) THEN
               ! Re must be the same, Ctrl must be different
                     
               p%secondVals(1) = p%Table(1)%UserProp
               
               DO iTable=2,p%NumTabs
                  IF ( p%Table(iTable)%Re /= p%Table(1)%Re )  THEN
                     CALL SetErrStat ( ErrID_Fatal, 'Fatal Error: airfoil file "'//TRIM( InitInput%FileName ) &
                                    //'", Table #'//TRIM( Num2LStr( iTable ) ) &
                                    //' does not have the same value for Re Property (Re) as the first table.', ErrStat, ErrMsg, RoutineName )
                     CALL Cleanup()
                     RETURN
                  ENDIF
                  p%secondVals(iTable) = p%Table(iTable)%UserProp
               END DO ! iTable
                     
            END IF
                  
               
            IF (.NOT. CheckValuesAreUniqueMonotonicIncreasing(p%secondVals)) THEN
            
               ErrMsg2 = 'Fatal Error: airfoil file "'//TRIM( InitInput%FileName ) &
                           //'", is not monotonic and increasing in the '
               IF (p%AFTabMod == AFITable_2Re) THEN
                  ErrMsg2 = trim(ErrMsg2)//' Re Property (Re).'
               ELSE
                  ErrMsg2 = trim(ErrMsg2)//' Ctrl Property (UserProp).'
               END IF
                     
               CALL SetErrStat ( ErrID_Fatal, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               CALL Cleanup()
               RETURN
               
            END IF
            
                
         END IF
      ELSE
         p%AFTabMod   = AFITable_1
      ENDIF ! ( p%NumTabs > 1 )

      
      do iTable = 1, p%NumTabs
               ! We need to deal with constant data.
         IF ( p%Table(iTable)%ConstData )  CYCLE  ! skip this table; it's just constant

            ! Allocate the arrays to hold spline coefficients.

         allocate ( p%Table(iTable)%SplineCoefs( p%Table(iTable)%NumAlf-1, size(p%Table(iTable)%Coefs,2), 0:3 ), STAT=ErrStat2 )
         if ( ErrStat2 /= 0 )  then
            call SetErrStat ( ErrStat2, 'Error allocating memory for the SplineCoefs array.', ErrStat, ErrMsg, RoutineName )
            call Cleanup()
            return
         end if

            
            ! Compute the one set of coefficients of the piecewise polynomials for the irregularly-spaced data.
            ! Unlike the 2-D interpolation in which we use diffent knots for each airfoil coefficient, we can do
            ! the 1-D stuff all at once.

            
            
         if ( p%InterpOrd == 3_IntKi ) then
               
            ! bjj: what happens at the end points (these are periodic, so we should maybe extend the tables to make sure the end point?) 

               ! use this for cubic splines:
            call CubicSplineInitM ( p%Table(iTable)%Alpha &
                                    , p%Table(iTable)%Coefs &
                                    , p%Table(iTable)%SplineCoefs &
                                    , ErrStat2, ErrMsg2 )
            call SetErrStat ( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                           
         else if ( p%InterpOrd == 1_IntKi ) then
               
               ! use this for linear interpolation (sets the higher order coeffs to zero):
               
               ! This is not the greatest way to get linear interpolation, but then we can use the same cubic spline routine
               ! later without checking interp order there
            call CubicLinSplineInitM ( p%Table(iTable)%Alpha &
                                       , p%Table(iTable)%Coefs &
                                       , p%Table(iTable)%SplineCoefs &
                                       , ErrStat2, ErrMsg2 )
            call SetErrStat ( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               
         else
               
            call SetErrStat ( ErrID_FATAL, 'Airfoil file "'//trim( InitInput%FileName ) &
                                       //'": InterpOrd must be 1 (linear) or 3 (cubic spline).', ErrStat, ErrMsg, RoutineName )
            call Cleanup()
            return
               
         end if
            
      end do

      CALL Cleanup ( )

      RETURN

   !=======================================================================
   CONTAINS
   !=======================================================================
      SUBROUTINE Cleanup ( )

         ! This subroutine cleans up the parent routine before exiting.
            ! Deallocate the temporary Coef array.

         IF ( ALLOCATED( Coef       ) ) DEALLOCATE( Coef )

         RETURN

      END SUBROUTINE Cleanup 

   END SUBROUTINE AFI_Init
   
   !=============================================================================
   !> This routine checks the init input values for AFI and makes sure they are valid
   !! before using them.
   SUBROUTINE AFI_ValidateInitInput(InitInput, ErrStat, ErrMsg)
      TYPE (AFI_InitInputType), INTENT(IN   )   :: InitInput                  ! This structure stores values that are set by the calling routine during the initialization phase.

      INTEGER(IntKi), INTENT(OUT)               :: ErrStat                    ! Error status
      CHARACTER(*), INTENT(OUT)                 :: ErrMsg                     ! Error message
      CHARACTER(*), PARAMETER                   :: RoutineName = 'AFI_Validate'

      ErrStat = ErrID_None
      ErrMsg  = ""
      
      if (InitInput%InCol_Alfa  < 0) call SetErrStat( ErrID_Fatal, 'InCol_Alfa must not be a negative number.', ErrStat, ErrMsg, RoutineName )
      if (InitInput%InCol_Cl    < 0) call SetErrStat( ErrID_Fatal, 'InCol_Cl must not be a negative number.', ErrStat, ErrMsg, RoutineName )
      if (InitInput%InCol_Cd    < 0) call SetErrStat( ErrID_Fatal, 'InCol_Cd must not be a negative number.', ErrStat, ErrMsg, RoutineName )
      if (InitInput%InCol_Cm    < 0) call SetErrStat( ErrID_Fatal, 'InCol_Cm must not be a negative number.', ErrStat, ErrMsg, RoutineName )
      if (InitInput%InCol_Cpmin < 0) call SetErrStat( ErrID_Fatal, 'InCol_Cpmin must not be a negative number.', ErrStat, ErrMsg, RoutineName )
      if (InitInput%AFTabMod /= AFITable_1 .and. InitInput%AFTabMod /= AFITable_2Re .and. InitInput%AFTabMod /= AFITable_2User) then
         call SetErrStat( ErrID_Fatal, 'AFTabMod must be 1, 2, or 3.', ErrStat, ErrMsg, RoutineName )
      end if
      
   
   END SUBROUTINE AFI_ValidateInitInput
  
   !=============================================================================
   SUBROUTINE ReadAFfile ( InitInp, NumCoefsIn, p, ErrStat, ErrMsg, UnEc )


         ! This routine reads an airfoil file.


         ! Argument declarations.

      TYPE (AFI_InitInputType), INTENT(IN)    :: InitInp                       ! This structure stores values that are set by the calling routine during the initialization phase.
      INTEGER(IntKi),    INTENT(  OUT)        :: ErrStat                       ! Error status
      INTEGER(IntKi),    INTENT(IN   )        :: NumCoefsIn                    ! The number of aerodynamic coefficients to be stored
      

      INTEGER,           INTENT(IN)           :: UnEc                          ! I/O unit for echo file. If present and > 0, write to UnEc.

      CHARACTER(*),      INTENT(  OUT)        :: ErrMsg                        ! Error message.

      TYPE (AFI_ParameterType), INTENT(INOUT)   :: p                          ! This structure stores all the module parameters that are set by AirfoilInfo during the initialization phase.


         ! Local declarations.

      REAL(ReKi)                              :: Coords   (2)                  ! An array to hold data from the airfoil-shape table.
      REAL(ReKi), ALLOCATABLE                 :: SiAry    (:)                  ! A temporary array to hold data from a table.
                                              
      INTEGER                                 :: Coef                          ! A DO index that points into the coefficient array.
      INTEGER                                 :: Cols2Parse                    ! The number of columns that must be read from the coefficient tables.
      INTEGER                                 :: CurLine                       ! The current line to be parsed in the FileInfo structure.
      INTEGER                                 :: NumAlf                        ! The number of rows in the current table.
      INTEGER                                 :: Row                           ! The row of a table to be parsed in the FileInfo structure.
      INTEGER                                 :: iTable                        ! The DO index for the tables.
                                              
      LOGICAL                                 :: BadVals                       ! A flag that indicates if the values in a table are invalid.

      TYPE (FileInfoType)                     :: FileInfo                      ! The derived type for holding the file information.
      INTEGER(IntKi)                          :: NumCoefsTab                   ! The number of aerodynamic coefficients to be stored.

      INTEGER(IntKi), parameter               :: DefaultInterpOrd = 1          ! value of default interp order
      INTEGER(IntKi)                          :: ErrStat2                      ! Error status local to this routine.
      CHARACTER(ErrMsgLen)                    :: ErrMsg2
      CHARACTER(*), PARAMETER                 :: RoutineName = 'ReadAFfile'
      CHARACTER(10)                           :: defaultStr
      CHARACTER(1024)                         :: PriPath
      
      TYPE (AFI_UA_BL_Default_Type), ALLOCATABLE :: CalcDefaults(:)            ! Whether to calculate default values for the UA parameters
      
      ErrStat = ErrID_None
      ErrMsg  = ""
      defaultStr = ""
      
      ! Getting parent folder of airfoils data (e.g. "Arifoils/")
      CALL GetPath( InitInp%FileName, PriPath )

         ! Process the (possibly) nested set of files.  This copies the decommented contents of
         ! AFI_FileInfo%FileName and the files it includes (both directly and indirectly) into
         ! the FileInfo structure that we can then parse.

      CALL ProcessComFile ( InitInp%FileName, FileInfo, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF (ErrStat >= AbortErrLev) THEN
            CALL Cleanup()
            RETURN
         END IF
         

         ! Process the airfoil shape information if it is included.

      CurLine = 1
      
      CALL ParseVarWDefault ( FileInfo, CurLine, 'InterpOrd', p%InterpOrd, DefaultInterpOrd, ErrStat2, ErrMsg2, UnEc )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF (ErrStat >= AbortErrLev) THEN
            CALL Cleanup()
            RETURN
         END IF
         
         ! RelThickness, default is 0.2 if user doesn't know it, only used for Boeing-Vertol UA model = 7
      CALL ParseVarWDefault ( FileInfo, CurLine, 'RelThickness', p%RelThickness, 0.2_ReKi, ErrStat2, ErrMsg2, UnEc )
         if (ErrStat2 >= AbortErrLev) then ! if the line is missing, set RelThickness = -1 and move on...
            p%RelThickness=-1 ! To trigger an error
            !call WrScr('Skipping. RelThickness not found on line 7 of Profile file: '//trim(InitInp%FileName) )
         else
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         end if
         
         
         ! NonDimArea is currently unused by AirfoilInfo or codes using AirfoilInfo.  GJH 9/13/2017
      CALL ParseVar ( FileInfo, CurLine, 'NonDimArea', p%NonDimArea, ErrStat2, ErrMsg2, UnEc )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         
         ! NumCoords, with the Coords data, is used for determining the blade shape (currently used 
         !   for visualization only).  This data (blade coordinates) is passed to the caller via 
         !   the InitOut%BladeShape data structure, and stored in p%XCoord, etc.,
         !   but is currently unused by AFI module.  GJH 9/13/2017
      CALL ParseVar ( FileInfo, CurLine, 'NumCoords' , p%NumCoords , ErrStat2, ErrMsg2, UnEc )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF (ErrStat >= AbortErrLev) THEN
            CALL Cleanup()
            RETURN
         END IF

      IF ( p%NumCoords > 0 )  THEN

         ALLOCATE ( p%X_Coord( p%NumCoords ) , STAT=ErrStat2 )
         IF ( ErrStat2 /= 0 )  THEN
            CALL SetErrStat ( ErrID_Fatal, 'Error allocating memory for p%X_Coord.', ErrStat, ErrMsg, RoutineName )
            CALL Cleanup()
            RETURN
         ENDIF

         ALLOCATE ( p%Y_Coord( p%NumCoords ) , STAT=ErrStat2 )
         IF ( ErrStat2 /= 0 )  THEN
            CALL SetErrStat ( ErrID_Fatal, 'Error allocating memory for p%Y_Coord.', ErrStat, ErrMsg, RoutineName )
            CALL Cleanup()
            RETURN
         ENDIF

         DO Row=1,p%NumCoords
            CALL ParseAry ( FileInfo, CurLine, 'X_Coord/Y_Coord', Coords, 2, ErrStat2, ErrMsg2, UnEc )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               IF (ErrStat >= AbortErrLev) THEN
                  CALL Cleanup()
                  RETURN
               END IF
            p%X_Coord(Row) = Coords(1)
            p%Y_Coord(Row) = Coords(2)
         ENDDO ! Row

      ENDIF
      
      ! Reading Boundary layer file for aeroacoustics
      CALL ParseVar ( FileInfo, CurLine, 'BL_file' , p%BL_file , ErrStat2, ErrMsg2, UnEc, IsPath=.true. )
         IF (ErrStat2 >= AbortErrLev) p%BL_file = "NOT_SET_IN_AIRFOIL_FILE"
         !CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( PathIsRelative( p%BL_file ) )  p%BL_file=trim(PriPath)//trim(p%BL_file)


         ! How many columns do we need to read in the input and how many total coefficients will be used?
      Cols2Parse = MAX( InitInp%InCol_Alfa, InitInp%InCol_Cl, InitInp%InCol_Cd, InitInp%InCol_Cm, InitInp%InCol_Cpmin )
      ALLOCATE ( SiAry( Cols2Parse ) , STAT=ErrStat2 )
      IF ( ErrStat2 /= 0 )  THEN
         CALL SetErrStat ( ErrID_Fatal, 'Error allocating memory for SiAry.', ErrStat, ErrMsg, RoutineName )
         CALL Cleanup()
         RETURN
      ENDIF


         ! Work through the multiple tables.

      CALL ParseVar ( FileInfo, CurLine, 'NumTabs' , p%NumTabs , ErrStat2, ErrMsg2, UnEc )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF (ErrStat >= AbortErrLev) THEN
            CALL Cleanup()
            RETURN
         END IF

      IF ( p%NumTabs < 1 )  THEN
         CALL SetErrStat ( ErrID_Fatal, 'NumTabs must be > 0.', ErrStat, ErrMsg, RoutineName )
         CALL Cleanup()
         RETURN
      ENDIF

      ALLOCATE ( p%Table( p%NumTabs ) , CalcDefaults(p%NumTabs), STAT=ErrStat2 )
      IF ( ErrStat2 /= 0 )  THEN
         CALL SetErrStat ( ErrID_Fatal, 'Error allocating memory for p%Table.', ErrStat, ErrMsg, RoutineName )
         CALL Cleanup()
         RETURN
      ENDIF

      DO iTable=1,p%NumTabs
        NumCoefsTab = NumCoefsIn
        
         CALL ParseVar ( FileInfo, CurLine, 'Re', p%Table(iTable)%Re, ErrStat2, ErrMsg2, UnEc )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            IF (ErrStat >= AbortErrLev) THEN
               CALL Cleanup()
               RETURN
            END IF
         p%Table(iTable)%Re = p%Table(iTable)%Re * 1.0e6  ! Entered in millions, so multiply here
         IF ( p%Table(iTable)%Re <= 0.0 )  THEN
               CALL SetErrStat ( ErrID_Severe, 'Re must be > 0 in "'//TRIM( InitInp%FileName ) &
                                      //'".'//NewLine//'  >> The error occurred on line #' &
                                      //TRIM( Num2LStr( FileInfo%FileLine(CurLine-1) ) )//'.', ErrStat, ErrMsg, RoutineName )
               CALL Cleanup()
               RETURN
         ENDIF ! ( p%Table(iTable)%Re <= 0.0 )

         CALL ParseVar ( FileInfo, CurLine, 'UserProp', p%Table(iTable)%UserProp, ErrStat2, ErrMsg2, UnEc )
         if (ErrStat2 >= AbortErrLev) then ! for backward compatibility of input files
            CALL ParseVar ( FileInfo, CurLine, 'Ctrl', p%Table(iTable)%UserProp, ErrStat2, ErrMsg2, UnEc )
         end if
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            IF (ErrStat >= AbortErrLev) THEN
               CALL Cleanup()
               RETURN
            END IF

         CALL ParseVar ( FileInfo, CurLine, 'InclUAdata', p%Table(iTable)%InclUAdata, ErrStat2, ErrMsg2, UnEc )
            if (ErrStat2 >= AbortErrLev) p%Table(iTable)%InclUAdata = .false. ! assume we don't have any UA data included, so we'll calculate it later.

         IF ( p%Table(iTable)%InclUAdata )  THEN

               CALL ParseVar ( FileInfo, CurLine, 'alpha0', p%Table(iTable)%UA_BL%alpha0, ErrStat2, ErrMsg2, UnEc )
               CalcDefaults(iTable)%alpha0 = ErrStat2 >= AbortErrLev
               if (.not. CalcDefaults(iTable)%alpha0) p%Table(iTable)%UA_BL%alpha0 = p%Table(iTable)%UA_BL%alpha0*D2R

               CALL ParseVar ( FileInfo, CurLine, 'alpha1', p%Table(iTable)%UA_BL%alpha1, ErrStat2, ErrMsg2, UnEc )
               CalcDefaults(iTable)%alpha1 = ErrStat2 >= AbortErrLev
               if (.not. CalcDefaults(iTable)%alpha1) p%Table(iTable)%UA_BL%alpha1 = p%Table(iTable)%UA_BL%alpha1*D2R

               CALL ParseVar ( FileInfo, CurLine, 'alpha2', p%Table(iTable)%UA_BL%alpha2, ErrStat2, ErrMsg2, UnEc )
               CalcDefaults(iTable)%alpha2 = ErrStat2 >= AbortErrLev
               if (.not. CalcDefaults(iTable)%alpha2) p%Table(iTable)%UA_BL%alpha2 = p%Table(iTable)%UA_BL%alpha2*D2R
               
               CALL ParseVar ( FileInfo, CurLine, 'alphaUpper', p%Table(iTable)%UA_BL%alphaUpper, ErrStat2, ErrMsg2, UnEc )
               CalcDefaults(iTable)%alphaUpper = ErrStat2 >= AbortErrLev
               if (.not. CalcDefaults(iTable)%alphaUpper) p%Table(iTable)%UA_BL%alphaUpper = p%Table(iTable)%UA_BL%alphaUpper*D2R
               
               CALL ParseVar ( FileInfo, CurLine, 'alphaLower', p%Table(iTable)%UA_BL%alphaLower, ErrStat2, ErrMsg2, UnEc )
               CalcDefaults(iTable)%alphaLower = ErrStat2 >= AbortErrLev
               if (.not. CalcDefaults(iTable)%alphaLower) p%Table(iTable)%UA_BL%alphaLower = p%Table(iTable)%UA_BL%alphaLower*D2R

               CALL ParseVar ( FileInfo, CurLine, 'eta_e', p%Table(iTable)%UA_BL%eta_e, ErrStat2, ErrMsg2, UnEc )
               CalcDefaults(iTable)%eta_e = ErrStat2 >= AbortErrLev

               CALL ParseVar ( FileInfo, CurLine, 'C_nalpha', p%Table(iTable)%UA_BL%C_nalpha, ErrStat2, ErrMsg2, UnEc )
               CalcDefaults(iTable)%C_nalpha = ErrStat2 >= AbortErrLev
            
               CALL ParseVar ( FileInfo, CurLine, 'C_lalpha', p%Table(iTable)%UA_BL%C_lalpha, ErrStat2, ErrMsg2, UnEc )
               CalcDefaults(iTable)%C_lalpha = ErrStat2 >= AbortErrLev
               
               CALL ParseVarWDefault ( FileInfo, CurLine, 'T_f0', p%Table(iTable)%UA_BL%T_f0, 3.0_ReKi, ErrStat2, ErrMsg2, UnEc )
               CalcDefaults(iTable)%T_f0 = ErrStat2 >= AbortErrLev
         
               CALL ParseVarWDefault ( FileInfo, CurLine, 'T_V0', p%Table(iTable)%UA_BL%T_V0, 6.0_ReKi, ErrStat2, ErrMsg2, UnEc )
               CalcDefaults(iTable)%T_V0 = ErrStat2 >= AbortErrLev
               
               CALL ParseVarWDefault ( FileInfo, CurLine, 'T_p', p%Table(iTable)%UA_BL%T_p, 1.7_ReKi, ErrStat2, ErrMsg2, UnEc )
               CalcDefaults(iTable)%T_p = ErrStat2 >= AbortErrLev

               CALL ParseVarWDefault ( FileInfo, CurLine, 'T_VL', p%Table(iTable)%UA_BL%T_VL, 11.0_ReKi, ErrStat2, ErrMsg2, UnEc )
               CalcDefaults(iTable)%T_VL = ErrStat2 >= AbortErrLev

               CALL ParseVarWDefault ( FileInfo, CurLine, 'b1', p%Table(iTable)%UA_BL%b1, .14_ReKi, ErrStat2, ErrMsg2, UnEc )
               CalcDefaults(iTable)%b1 = ErrStat2 >= AbortErrLev

               CALL ParseVarWDefault ( FileInfo, CurLine, 'b2', p%Table(iTable)%UA_BL%b2, .53_ReKi, ErrStat2, ErrMsg2, UnEc )
               CalcDefaults(iTable)%b2 = ErrStat2 >= AbortErrLev

               CALL ParseVarWDefault ( FileInfo, CurLine, 'b5', p%Table(iTable)%UA_BL%b5, 5.0_ReKi, ErrStat2, ErrMsg2, UnEc )
               CalcDefaults(iTable)%b5 = ErrStat2 >= AbortErrLev

               CALL ParseVarWDefault ( FileInfo, CurLine, 'A1', p%Table(iTable)%UA_BL%A1, .3_ReKi, ErrStat2, ErrMsg2, UnEc )
               CalcDefaults(iTable)%A1 = ErrStat2 >= AbortErrLev

               CALL ParseVarWDefault ( FileInfo, CurLine, 'A2', p%Table(iTable)%UA_BL%A2, .7_ReKi, ErrStat2, ErrMsg2, UnEc )
               CalcDefaults(iTable)%A2 = ErrStat2 >= AbortErrLev

               CALL ParseVarWDefault ( FileInfo, CurLine, 'A5', p%Table(iTable)%UA_BL%A5, 1.0_ReKi, ErrStat2, ErrMsg2, UnEc )
               CalcDefaults(iTable)%A5 = ErrStat2 >= AbortErrLev

               CALL ParseVar ( FileInfo, CurLine, 'S1', p%Table(iTable)%UA_BL%S1, ErrStat2, ErrMsg2, UnEc )
               CalcDefaults(iTable)%S1 = ErrStat2 >= AbortErrLev

               CALL ParseVar ( FileInfo, CurLine, 'S2', p%Table(iTable)%UA_BL%S2, ErrStat2, ErrMsg2, UnEc )
               CalcDefaults(iTable)%S2 = ErrStat2 >= AbortErrLev

               CALL ParseVar ( FileInfo, CurLine, 'S3', p%Table(iTable)%UA_BL%S3, ErrStat2, ErrMsg2, UnEc )
               CalcDefaults(iTable)%S3 = ErrStat2 >= AbortErrLev

               CALL ParseVar ( FileInfo, CurLine, 'S4', p%Table(iTable)%UA_BL%S4, ErrStat2, ErrMsg2, UnEc )
               CalcDefaults(iTable)%S4 = ErrStat2 >= AbortErrLev

               CALL ParseVar ( FileInfo, CurLine, 'Cn1', p%Table(iTable)%UA_BL%Cn1, ErrStat2, ErrMsg2, UnEc )
               CalcDefaults(iTable)%Cn1 = ErrStat2 >= AbortErrLev

               CALL ParseVar ( FileInfo, CurLine, 'Cn2', p%Table(iTable)%UA_BL%Cn2, ErrStat2, ErrMsg2, UnEc )
               CalcDefaults(iTable)%Cn2 = ErrStat2 >= AbortErrLev

               CALL ParseVarWDefault ( FileInfo, CurLine, 'St_sh', p%Table(iTable)%UA_BL%St_sh, .19_ReKi, ErrStat2, ErrMsg2, UnEc )
               CalcDefaults(iTable)%St_sh = ErrStat2 >= AbortErrLev

               CALL ParseVar ( FileInfo, CurLine, 'Cd0', p%Table(iTable)%UA_BL%Cd0, ErrStat2, ErrMsg2, UnEc )
               CalcDefaults(iTable)%Cd0 = ErrStat2 >= AbortErrLev

               CALL ParseVar ( FileInfo, CurLine, 'Cm0', p%Table(iTable)%UA_BL%Cm0, ErrStat2, ErrMsg2, UnEc )
               CalcDefaults(iTable)%Cm0 = ErrStat2 >= AbortErrLev

               CALL ParseVar ( FileInfo, CurLine, 'k0', p%Table(iTable)%UA_BL%k0, ErrStat2, ErrMsg2, UnEc )
               CalcDefaults(iTable)%k0 = ErrStat2 >= AbortErrLev

               CALL ParseVar ( FileInfo, CurLine, 'k1', p%Table(iTable)%UA_BL%k1, ErrStat2, ErrMsg2, UnEc )
               CalcDefaults(iTable)%k1 = ErrStat2 >= AbortErrLev

               CALL ParseVar ( FileInfo, CurLine, 'k2', p%Table(iTable)%UA_BL%k2, ErrStat2, ErrMsg2, UnEc )
               CalcDefaults(iTable)%k2 = ErrStat2 >= AbortErrLev

               CALL ParseVar ( FileInfo, CurLine, 'k3', p%Table(iTable)%UA_BL%k3, ErrStat2, ErrMsg2, UnEc )
               CalcDefaults(iTable)%k3 = ErrStat2 >= AbortErrLev

               CALL ParseVar ( FileInfo, CurLine, 'k1_hat', p%Table(iTable)%UA_BL%k1_hat, ErrStat2, ErrMsg2, UnEc )
               CalcDefaults(iTable)%k1_hat = ErrStat2 >= AbortErrLev

               CALL ParseVarWDefault ( FileInfo, CurLine, 'x_cp_bar', p%Table(iTable)%UA_BL%x_cp_bar, .2_ReKi, ErrStat2, ErrMsg2, UnEc )
               CalcDefaults(iTable)%x_cp_bar = ErrStat2 >= AbortErrLev
                  
               CALL ParseVarWDefault ( FileInfo, CurLine, 'UACutout', p%Table(iTable)%UA_BL%UACutout, 45.0_ReKi, ErrStat2, ErrMsg2, UnEc )
               CalcDefaults(iTable)%UACutout = ErrStat2 >= AbortErrLev
               if (.not. CalcDefaults(iTable)%UACutout ) p%Table(iTable)%UA_BL%UACutout = p%Table(iTable)%UA_BL%UACutout*D2R

               CALL ParseVarWDefault ( FileInfo, CurLine, 'UACutout_delta', p%Table(iTable)%UA_BL%UACutout_delta, 5.0_ReKi, ErrStat2, ErrMsg2, UnEc )
               CalcDefaults(iTable)%UACutout_delta = ErrStat2 >= AbortErrLev
               if (.not. CalcDefaults(iTable)%UACutout_delta) p%Table(iTable)%UA_BL%UACutout_delta = p%Table(iTable)%UA_BL%UACutout_delta*D2R

               CALL ParseVarWDefault ( FileInfo, CurLine, 'filtCutOff', p%Table(iTable)%UA_BL%filtCutOff, 0.5_ReKi, ErrStat2, ErrMsg2, UnEc )
               CalcDefaults(iTable)%filtCutOff = ErrStat2 >= AbortErrLev
                  
               IF (ErrStat >= AbortErrLev) THEN
                  CALL Cleanup()
                  RETURN
               END IF
         ELSE
         
            ! everything is default ( we could just attempt to read these from the file, but it will be faster to use the default settings of .true. for each variable)
         
            p%Table(iTable)%InclUAdata = .true. ! make sure we are calculating this for each table

         ENDIF ! ( p%Table(iTable)%InclUAdata )

         if ( p%Table(iTable)%InclUAdata ) then
            p%ColUAf    = NumCoefsIn + 1 ! column for f_st for the UA models
            NumCoefsTab = NumCoefsIn + 3 ! total number of columns if we have UA on (for f_st, cl/cn_fs, cl/cn_fa)
         end if

         
         CALL ParseVar ( FileInfo, CurLine, 'NumAlf', p%Table(iTable)%NumAlf, ErrStat2, ErrMsg2, UnEc )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            IF (ErrStat >= AbortErrLev) THEN
               CALL Cleanup()
               RETURN
            END IF

         IF ( p%Table(iTable)%NumAlf < 1 )  THEN
            CALL SetErrStat( ErrID_Fatal, 'NumAlf must be a positive number on line #' &
                           //TRIM( Num2LStr( FileInfo%FileLine(CurLine-1) ) )//' in "'//TRIM( InitInp%FileName )//'".', ErrStat, ErrMsg, RoutineName )
            CALL Cleanup()
            RETURN
         ENDIF ! ( Test for valid values for NumAlf )


            ! Allocate the arrays for the airfoil coefficients.

         ALLOCATE ( p%Table(iTable)%Alpha( p%Table(iTable)%NumAlf ), STAT=ErrStat2 )
         IF ( ErrStat2 /= 0 )  THEN
            CALL SetErrStat( ErrID_Fatal, 'Error allocating memory for p%Table%Alpha.', ErrStat, ErrMsg, RoutineName )
            CALL Cleanup()
            RETURN
         ENDIF

         ALLOCATE ( p%Table(iTable)%Coefs( p%Table(iTable)%NumAlf, NumCoefsTab ), STAT=ErrStat2 )
         IF ( ErrStat2 /= 0 )  THEN
            CALL SetErrStat( ErrID_Fatal, 'Error allocating memory for p%Table%Coefs.', ErrStat, ErrMsg, RoutineName )
            CALL Cleanup()
            RETURN
         ENDIF
         p%Table(iTable)%Coefs = 0.0_ReKi

         DO Row=1,p%Table(iTable)%NumAlf

            CALL ParseAry ( FileInfo, CurLine, 'CoeffData', SiAry, Cols2Parse, ErrStat2, ErrMsg2, UnEc )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               IF (ErrStat >= AbortErrLev) THEN
                  CALL Cleanup()
                  RETURN
               END IF

            p%Table(iTable)%Alpha(Row        ) = SiAry(InitInp%InCol_Alfa)*D2R
            p%Table(iTable)%Coefs(Row,p%ColCl) = SiAry(InitInp%InCol_Cl  )
            p%Table(iTable)%Coefs(Row,p%ColCd) = SiAry(InitInp%InCol_Cd  )

            IF ( InitInp%InCol_Cm    > 0 ) p%Table(iTable)%Coefs(Row,p%ColCm   ) = SiAry(InitInp%InCol_Cm)
            IF ( InitInp%InCol_Cpmin > 0 ) p%Table(iTable)%Coefs(Row,p%ColCpmin) = SiAry(InitInp%InCol_Cpmin)

         ENDDO ! Row

            ! check that not all the values are constant
         IF ( p%Table(iTable)%NumAlf < 3 )  THEN
            p%Table(iTable)%ConstData = .TRUE. ! we can't do splines with this many points, so it must be constant
         ELSE
               ! check if the columns change with alpha
            p%Table(iTable)%ConstData = .TRUE.
ALPHA_LOOP: DO Row=1,p%Table(iTable)%NumAlf-1
               DO Coef=1,NumCoefsIn ! don't check the additional columns from UA
                  IF ( .NOT. EqualRealNos( p%Table(iTable)%Coefs(Row,Coef), p%Table(iTable)%Coefs(Row+1,Coef) ) )  THEN
                     p%Table(iTable)%ConstData = .FALSE.
                     EXIT ALPHA_LOOP
                  ENDIF
               END DO ! Coef
            END DO ALPHA_LOOP
         END IF
         
         
         if ( p%Table(iTable)%ConstData ) then
            p%Table(iTable)%InclUAdata = .false.
         else
            call CalculateUACoeffs(CalcDefaults(iTable), p%Table(iTable), p%ColCl, p%ColCd, p%ColCm, p%ColUAf, InitInp%UAMod)
         end if

            ! Let's make sure that the data go from -Pi to Pi and that the values are the same for both
            ! unless there is only one point.

         NumAlf  = p%Table(iTable)%NumAlf
         if (NumAlf > 1) then
            BadVals = .FALSE.
            
            if (.not. p%Table(iTable)%ConstData) then
               IF ( .NOT. EqualRealNos( p%Table(iTable)%Alpha(1), -Pi ) )  THEN
                  BadVals = .TRUE.
               ENDIF
               IF ( .NOT. EqualRealNos( p%Table(iTable)%Alpha(NumAlf), Pi ) )  THEN
                  BadVals = .TRUE.
               ENDIF
            end if
               
            DO Coef=1,NumCoefsIn ! don't check the additional columns from UA
               IF ( .NOT. EqualRealNos( p%Table(iTable)%Coefs(1,Coef), p%Table(iTable)%Coefs(NumAlf,Coef) ) )  THEN
                  BadVals = .TRUE.
               ENDIF
            ENDDO ! Coef
            
            IF ( BadVals )  THEN
               if (p%Table(iTable)%ConstData) then
                  ErrStat2 = ErrID_Fatal
               else
                  ErrStat2 = ErrID_Warn
               end if
            
               CALL SetErrStat( ErrStat2, &
                  'Airfoil data should go from -180 degrees to 180 degrees and the coefficients at the ends should be the same.', ErrStat, ErrMsg, RoutineName )
            ENDIF
         end if
      ENDDO ! iTable


      DO iTable=1,p%NumTabs
         if ( .not. p%Table(iTable)%InclUAdata )  then
            p%ColUAf = 0 ! in case some tables have UA data and others don't; this is not set on a per-table basis
            exit ! exit loop
         end if
      ENDDO ! iTable

      CALL Cleanup( )

      RETURN

      !=======================================================================
      CONTAINS
      !=======================================================================
                  
         SUBROUTINE Cleanup ()

            ! This subroutine cleans up all the allocatable arrays, sets the error status/message and closes the binary file
               
            CALL NWTC_Library_DestroyFileInfoType (FileInfo, ErrStat2, ErrMsg2)
            IF ( ALLOCATED(CalcDefaults) ) DEALLOCATE(CalcDefaults) ! this type contains only logicals, so no need to call a destroy routine
            IF ( ALLOCATED(SiAry) ) DEALLOCATE(SiAry)
            
         END SUBROUTINE Cleanup 


   END SUBROUTINE ReadAFfile
!----------------------------------------------------------------------------------------------------------------------------------  
   SUBROUTINE CalculateUACoeffs(CalcDefaults,p,ColCl,ColCd,ColCm,ColUAf,UAMod)
      TYPE (AFI_UA_BL_Default_Type),intent(in):: CalcDefaults
      TYPE (AFI_Table_Type),    intent(inout) :: p                             ! This structure stores all the module parameters that are set by AirfoilInfo during the initialization phase.
      integer(IntKi),           intent(in   ) :: ColCl                         ! column for cl
      integer(IntKi),           intent(in   ) :: ColCd                         ! column for cd
      integer(IntKi),           intent(in   ) :: ColCm                         ! column for cm
      integer(IntKi),           intent(in   ) :: ColUAf                        ! column for UA f_st (based on Cl or cn)
      integer(IntKi),           intent(in   ) :: UAMod                         ! UA model; determines how to compute f_st?
   
      INTEGER(IntKi)                          :: Row                           ! The row of a table to be parsed in the FileInfo structure.
      INTEGER(IntKi)                          :: col_fs                        ! column for UA cn/cl_fs (fully separated cn or cl)
      INTEGER(IntKi)                          :: col_fa                        ! column for UA cn/cl_fa (fully attached cn or cl)
      
      INTEGER(IntKi)                          :: iCdMin
      INTEGER(IntKi)                          :: iHighLimit, iLowLimit
      INTEGER(IntKi)                          :: iHigh,  iLow
      INTEGER(IntKi)                          :: iHigh2, iLow2
      INTEGER(IntKi)                          :: iGuess, iUpper, iLower, i(1)
      INTEGER(IntKi)                          :: nRoots
      
      LOGICAL                                 :: UA_f_cn

      ! note that we don't get here with constant data, so NumAlf>2
      REAL(ReKi)                              :: roots(p%NumAlf)
      REAL(ReKi)                              :: cn(p%NumAlf)
      REAL(ReKi)                              :: ClSlope_raw(p%NumAlf-1)
      REAL(ReKi)                              :: CnSlope_raw(p%NumAlf-1)
      
      REAL(ReKi)                              :: ClSlope_(p%NumAlf-1)
      REAL(ReKi)                              :: CnSlope_(p%NumAlf-1)
      REAL(ReKi)                              :: alpha_(p%NumAlf-1)
      REAL(ReKi)                              :: alphaAtCdMin
      REAL(ReKi)                              :: CnSlopeAtCdMin
      REAL(ReKi)                              :: maxCnSlope
      
      REAL(ReKi) , PARAMETER                  :: CnSlopeThreshold = 0.90;
      REAL(ReKi) , PARAMETER                  :: fAtCriticalCn    = 0.7;
      REAL(ReKi)                              :: LimitAlphaRange
      REAL(ReKi)                              :: Default_Cn_alpha
      REAL(ReKi)                              :: Default_Cl_alpha
      REAL(ReKi)                              :: Default_alpha0
      REAL(ReKi)                              :: alphaMargin
      

      INTEGER(IntKi)                          :: ErrStat2                      ! Error status local to this routine.
      CHARACTER(ErrMsgLen)                    :: ErrMsg2
      CHARACTER(*), PARAMETER                 :: RoutineName = 'CalculateUACoeffs'

      if ( UAMod == UA_HGMV360 ) then
         LimitAlphaRange = TwoPi ! range we're limiting our equations to (in radians)
      else
         LimitAlphaRange =  20.0_ReKi * D2R ! range we're limiting our equations to (in radians)
      end if

      col_fs = ColUAf + 1
      col_fa = col_fs  + 1
      UA_f_cn = UAMod /= UA_HGM .and. UAMod /= UA_Oye ! these models use cl instead of cn
      
      if ( p%InclUAdata )  then

            ! these variables are for the unused UAMod=1 (UA_Baseline) model;
            ! TO DO: fix the calculations if we ever allow this model to be used
         if (CalcDefaults%eta_e          ) p%UA_BL%eta_e          =  1.00_ReKi ! it says it's supposed to be in the range 085 - 0.95, but then when FLookup is true, it uses 1.
         if (CalcDefaults%k0             ) p%UA_BL%k0             =  0.00_ReKi
         if (CalcDefaults%k1             ) p%UA_BL%k1             =  0.00_ReKi
         if (CalcDefaults%k2             ) p%UA_BL%k2             =  0.00_ReKi
         if (CalcDefaults%k3             ) p%UA_BL%k3             =  0.00_ReKi
         if (CalcDefaults%k1_hat         ) p%UA_BL%k1_hat         =  0.00_ReKi
         if (CalcDefaults%S1             ) p%UA_BL%S1             =  0.00_ReKi
         if (CalcDefaults%S2             ) p%UA_BL%S2             =  0.00_ReKi
         if (CalcDefaults%S3             ) p%UA_BL%S3             =  0.00_ReKi
         if (CalcDefaults%S4             ) p%UA_BL%S4             =  0.00_ReKi
         
            ! Set to default values:
         if (CalcDefaults%T_f0           ) p%UA_BL%T_f0           =  3.00_ReKi
         if (CalcDefaults%T_V0           ) p%UA_BL%T_V0           =  6.00_ReKi
         if (CalcDefaults%T_p            ) p%UA_BL%T_p            =  1.70_ReKi
         if (CalcDefaults%T_VL           ) p%UA_BL%T_VL           = 11.00_ReKi
         if (CalcDefaults%b1             ) p%UA_BL%b1             =  0.14_ReKi
         if (CalcDefaults%b2             ) p%UA_BL%b2             =  0.53_ReKi
         if (CalcDefaults%b5             ) p%UA_BL%b5             =  5.00_ReKi
         if (CalcDefaults%A1             ) p%UA_BL%A1             =  0.30_ReKi
         if (CalcDefaults%A2             ) p%UA_BL%A2             =  0.70_ReKi
         if (CalcDefaults%A5             ) p%UA_BL%A5             =  1.00_ReKi
         if (CalcDefaults%x_cp_bar       ) p%UA_BL%x_cp_bar       =  0.20_ReKi
         if (CalcDefaults%filtCutOff     ) p%UA_BL%filtCutOff     =  0.50_ReKi
         
         if (UAMod == UA_HGMV360) then ! set defaults for this model (note: we don't turn off UA)
            if (CalcDefaults%St_sh          ) p%UA_BL%St_sh          =  0.14_ReKi
            if (CalcDefaults%UACutout       ) p%UA_BL%UACutout       =  TwoPi*2              ! don't turn off UA for this model
            if (CalcDefaults%UACutout_delta ) p%UA_BL%UACutout_delta =  D2R                  ! begin turning off 1 degrees before UAcutout (if UACutout is large enough, we don't turn off UAcutout)
         else
            if (CalcDefaults%St_sh          ) p%UA_BL%St_sh          =  0.19_ReKi
            if (CalcDefaults%UACutout       ) p%UA_BL%UACutout       = 45.00_ReKi*D2R        ! turn off UA at 45 degrees
            if (CalcDefaults%UACutout_delta ) p%UA_BL%UACutout_delta =  5.00_ReKi*D2R        ! begin turning off 5 degrees before UAcutout
         end if
         
         p%UA_BL%UACutout_blend  = max(0.0_ReKi, abs(p%UA_BL%UACutout) - abs(p%UA_BL%UACutout_delta))
         
         !-------------------------------------
         ! Calculate based on airfoil polar:
         !-------------------------------------
         iCdMin = minloc(p%Coefs(:,ColCd),DIM=1, MASK=abs(p%alpha)<=LimitAlphaRange)
         
         if ( (maxval(p%Coefs(:,ColCd),DIM=1, MASK=abs(p%alpha)<=LimitAlphaRange) ==        &
               minval(p%Coefs(:,ColCd),DIM=1, MASK=abs(p%alpha)<=LimitAlphaRange)    ) .or. &
               maxval(p%Coefs(:,ColCl),DIM=1, MASK=abs(p%alpha)<=LimitAlphaRange) < 0.01 ) then
                 
               ! Cylinder polar perhaps?
                 
               if (CalcDefaults%Cd0)        p%UA_BL%Cd0 = p%Coefs(iCdMin,ColCd)
               if (CalcDefaults%alpha0)     p%UA_BL%alpha0 = 0;
               if (CalcDefaults%C_nalpha)   p%UA_BL%C_nalpha = 0;
               if (CalcDefaults%C_lalpha)   p%UA_BL%C_lalpha = 0;
               if (CalcDefaults%Cm0)        p%UA_BL%Cm0 = 0;
               if (CalcDefaults%alpha1)     p%UA_BL%alpha1 = 10*D2R;
               if (CalcDefaults%alpha2)     p%UA_BL%alpha2 = -10*D2R;
               if (CalcDefaults%Cn1)        p%UA_BL%Cn1 = 0;
               if (CalcDefaults%Cn2)        p%UA_BL%Cn2 = 0;
               if (CalcDefaults%alphaLower) p%UA_BL%alphaLower = -5*D2R;
               if (CalcDefaults%alphaUpper) p%UA_BL%alphaUpper = 5*D2R;
               
            if (.not. UA_f_cn) then
               call ComputeUASeparationFunction_onCl(p, ColCl, ColUAf, col_fs, col_fa)
            else
            
               if ( UAMod == UA_HGMV360 ) then
                  p%UA_BL%Cd0 = 0.0_ReKi  ! setting this to 0 so that Cn gets calculated properly elsewhere in this code
               end if
               Cn = Calculate_Cn(alpha=p%alpha, cl=p%Coefs(:,ColCl), cd=p%Coefs(:,ColCd), cd0=p%UA_BL%Cd0)
               call ComputeUA360_AttachedFlow(p, ColUAf, Cn, iLower, iUpper)
               call ComputeUA360_updateSeparationF( p, ColUAf, Cn, iLower, iUpper )
               call ComputeUA360_updateCnSeparated( p, ColUAf, Cn, iLower )

            end if
         else
               ! if Cd is constant, does this cause issues???
            alphaAtCdMin = p%alpha(iCdMin)
         
               ! compute cn:
            if (UAMod == UA_HGMV360) then
               !call Compute_iLoweriUpper(p, iLower, iUpper)
               !if (CalcDefaults%Cd0) p%UA_BL%Cd0 = minval( p%Coefs(iLower:iUpper, ColCd) )
               p%UA_BL%Cd0 = 0.0_ReKi  ! setting this to 0 so that Cn gets calculated properly elsewhere in this code
            else
               if (CalcDefaults%Cd0) p%UA_BL%Cd0 = p%Coefs(iCdMin,ColCd)
            end if
            cn = Calculate_Cn(alpha=p%alpha, cl=p%Coefs(:,ColCl), cd=p%Coefs(:,ColCd), cd0=p%UA_BL%Cd0)

               ! compute cn and cl slopes (raw):
            do Row=1,p%NumAlf-1
               CnSlope_raw(Row) = (      cn(Row+1)       - cn(Row)            ) / (p%alpha(Row+1) - p%alpha(Row))
               ClSlope_raw(Row) = ( p%Coefs(Row+1,ColCl) - p%Coefs(Row,ColCl) ) / (p%alpha(Row+1) - p%alpha(Row))
               alpha_(     Row) = 0.5_ReKi * (p%alpha(Row+1) + p%alpha(Row))
            end do
         
               ! smooth cn slope for better calculations later:
            call kernelSmoothing(alpha_, CnSlope_raw, kernelType_TRIWEIGHT, 2.0_ReKi*D2R, CnSlope_)
            call kernelSmoothing(alpha_, ClSlope_raw, kernelType_TRIWEIGHT, 2.0_ReKi*D2R, ClSlope_)
         
            iGuess = iCdMin
            CnSlopeAtCdMin = InterpStp( alphaAtCdMin, alpha_, CnSlope_, iGuess, p%NumAlf-1 )
         
               ! find bounding indices for limitAlphaRange
            iHighLimit = min( maxloc( alpha_ , DIM=1, MASK=alpha_ <  LimitAlphaRange) + 1, size(alpha_) ) ! we can limit this to some range
            iLowLimit  = max( minloc( alpha_ , DIM=1, MASK=alpha_ > -LimitAlphaRange) - 1, 1            ) ! we can limit this to some range
            if (iHighLimit - iLowLimit < 3) iHighLimit = min(iLowLimit+2,size(alpha_)) ! this could still be an issue if we don't have very many points in the airfoil table. If that's the case, this data is not worth anything anyway
            if (iHighLimit - iLowLimit < 3) iLowLimit  = max(iHighLimit-2,1) ! this could still be an issue if we don't have very many points in the airfoil table. If that's the case, this data is not worth anything anyway
            
               ! find alphaUpper (using smoothed Cn values):
            if (CalcDefaults%alphaUpper) then
               iHigh = iHighLimit
               if (iHigh<iCdMin) iHigh = p%NumAlf-1 !this should be an error?
            
               maxCnSlope = CnSlopeAtCdMin
               do Row=iCdMin, iHigh
                  iHigh2 = Row
                  if (CnSlope_(Row) > maxCnSlope) then
                     maxCnSlope = CnSlope_(Row)
                  else if (CnSlope_(Row) < CnSlopeThreshold*maxCnSlope) then
                     exit
                  end if
               end do

               if (iHigh2 == iCdMin) then
                  p%UA_BL%alphaUpper = alphaAtCdMin;
               else
                  iHigh2 = min(max(1, iHigh2-1), p%NumAlf-1 )
                  p%UA_BL%alphaUpper = alpha_(iHigh2);
               end if

            else
               iHigh2 = iHighLimit ! initialize for use in alphaLower if no alphaUpper default is requested
            end if
         
            
               !find alphaLower
            if (CalcDefaults%alphaLower) then
               maxCnSlope = CnSlopeAtCdMin
               
               iLow = iLowLimit
               iHigh = max( 1, min(iHigh2, iCdMin, p%NumAlf-1) )
               if (iHigh < iLow) iLow = 1

               do Row = iHigh,iLow,-1
                  iLow2 = Row
                  if (CnSlope_(Row) > maxCnSlope) then
                     maxCnSlope = CnSlope_(Row);
                  else if ( CnSlope_(Row) < CnSlopeThreshold*maxCnSlope ) then
                     exit
                  end if
               end do
            
               if (iLow2 == iCdMin) then
                  p%UA_BL%alphaLower = alphaAtCdMin;
               else
                  iLow2 = min(max(1, iLow2+1), p%NumAlf-1 )
                  p%UA_BL%alphaLower = alpha_(iLow2);
               end if
            end if
         
            !------------------------------------
            ! Note: C_nalpha, C_lalpha, and alpha0 are not used in HGMV360
            !------------------------------------
            
            if (CalcDefaults%C_nalpha .or. CalcDefaults%C_lalpha .or. CalcDefaults%alpha0) then
               
               alphaMargin = 0.2*( p%UA_BL%alphaUpper - p%UA_BL%alphaLower );
               !mask = p%alpha >= p%UA_BL%alphaLower+alphaMargin & p%alpha <= p%UA_BL%alphaUpper-alphaMargin;
            
               iLow2 = iLowLimit
               do while (iLow2 < iHighLimit-1 .and. p%alpha(iLow2) <  p%UA_BL%alphaLower + alphaMargin) 
                  iLow2 = iLow2 + 1
               end do

               iHigh2 = iHighLimit
               do while (iHigh2 > iLow2+1 .and. p%alpha(iHigh2) >  p%UA_BL%alphaUpper - alphaMargin) 
                  iHigh2 = iHigh2 - 1
               end do

               call Calculate_C_alpha(p%alpha(iLow2:iHigh2), Cn(iLow2:iHigh2), p%Coefs(iLow2:iHigh2,ColCl), Default_Cn_alpha, Default_Cl_alpha, Default_alpha0, ErrStat2, ErrMsg2)
         
               if (CalcDefaults%C_nalpha) p%UA_BL%C_nalpha = Default_Cn_alpha
               if (CalcDefaults%C_lalpha) p%UA_BL%C_lalpha = Default_Cl_alpha
               if (CalcDefaults%alpha0)   p%UA_BL%alpha0   = Default_alpha0
            end if

            if (CalcDefaults%Cm0) then
               if (ColCm > 0) then
                  iGuess = p%NumAlf/2 ! guess: start in the center
                  p%UA_BL%Cm0 = InterpStp( p%UA_BL%alpha0, p%alpha, p%Coefs(:,ColCm), iGuess, p%NumAlf )
               else
                  p%UA_BL%Cm0 = 0.0_ReKi
               end if
            end if
         
            if (.not. UA_f_cn) then !
               call ComputeUASeparationFunction_onCl(p, ColCl, ColUAf, col_fs, col_fa)
               call Compute_iLoweriUpper(p, iLower, iUpper) ! calculating iLower and iUpper here (for alpha1 and alpha2)
            else
               call ComputeUA360_AttachedFlow(p, ColUAf, Cn, iLower, iUpper)
               call ComputeUA360_updateSeparationF( p, ColUAf, Cn, iLower, iUpper )
               call ComputeUA360_updateCnSeparated( p, ColUAf, Cn, iLower )
            end if
            

               ! alpha1
            if (CalcDefaults%alpha1) then
               iGuess = max(1, minloc( p%alpha , DIM=1, MASK=p%alpha >= p%UA_BL%alphaUpper .and. p%Coefs(:,ColUAf) <= fAtCriticalCn))
               call fZeros(p%alpha(iUpper:), fAtCriticalCn - p%Coefs(iUpper:,ColUAf), roots, nRoots)

               if (nRoots==1) then
                  p%UA_BL%alpha1 = roots(1)
               elseif (nRoots>1) then
                  i = minloc( abs(roots(1:nRoots) - p%alpha(iGuess) ), DIM=1 ) ! find root closest to guess
                  p%UA_BL%alpha1 = roots(i(1))
               else
                  p%UA_BL%alpha1 = p%alpha(iGuess)
               end if
            end if
         
               ! alpha2
            if (CalcDefaults%alpha2) then
               iGuess = maxloc( p%alpha , DIM=1, MASK=p%alpha <= p%UA_BL%alphaLower .and. p%Coefs(:,ColUAf) <= fAtCriticalCn)
               call fZeros(p%alpha(:iLower), fAtCriticalCn - p%Coefs(:iLower,ColUAf), roots, nRoots)

               if (nRoots==1) then
                  p%UA_BL%alpha2 = roots(1)
               elseif (nRoots>1) then
                  i = minloc( abs(roots(1:nRoots) - p%alpha(iGuess) ), DIM=1 ) ! find root closest to guess
                  p%UA_BL%alpha2 = roots(i(1))
               else
                  p%UA_BL%alpha2 = p%alpha(iGuess)
               end if
            
            end if
         
               ! Cn1
            if (CalcDefaults%Cn1) then
               iGuess = iHighLimit
               p%UA_BL%Cn1 = InterpStp( p%UA_BL%alpha1, p%alpha, cn, iGuess, p%NumAlf )
            end if
         
               ! Cn2
            if (CalcDefaults%Cn2) then
               iGuess = iLowLimit
               p%UA_BL%Cn2 = InterpStp( p%UA_BL%alpha2, p%alpha, cn, iGuess, p%NumAlf )
            end if
         
         end if ! not a circular polar
         
         if ( UA_f_cn ) then
            iGuess = iLowLimit
            p%UA_BL%c_alphaLower = InterpStp(p%UA_BL%alphaLower, p%alpha, cn, iGuess, p%NumAlf)
            iGuess = iHighLimit
            p%UA_BL%c_alphaUpper = InterpStp(p%UA_BL%alphaUpper, p%alpha, cn, iGuess, p%NumAlf)
         else
            iGuess = iLowLimit
            p%UA_BL%c_alphaLower = InterpStp(p%UA_BL%alphaLower, p%alpha, p%Coefs(:,ColCl), iGuess, p%NumAlf)
            iGuess = iHighLimit
            p%UA_BL%c_alphaUpper = InterpStp(p%UA_BL%alphaUpper, p%alpha, p%Coefs(:,ColCl), iGuess, p%NumAlf)
         end if
         
      end if ! UA is included

   END SUBROUTINE CalculateUACoeffs
!----------------------------------------------------------------------------------------------------------------------------------
   FUNCTION Calculate_Cn (alpha, Cl, Cd, Cd0) RESULT(Cn)
      REAL(ReKi),               intent(in   ) :: alpha(:)                   ! alpha
      REAL(ReKi),               intent(in   ) :: Cl(:)                      ! cl
      REAL(ReKi),               intent(in   ) :: Cd(:)                      ! cd
      REAL(ReKi),               intent(in   ) :: Cd0
      REAL(ReKi)                              :: Cn(size(alpha))            ! cn (result of this function)
      
      integer(IntKi)                          :: NumAlf
      integer(IntKi)                          :: Row
   
      NumAlf = size(alpha)
      
      do Row=1,NumAlf
         cn(Row) = Cl(Row)*cos(alpha(Row)) + (Cd(Row) - Cd0)*sin(alpha(Row))
      end do
      
   END FUNCTION Calculate_Cn
!----------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE Calculate_C_alpha(alpha, Cn, Cl, Default_Cn_alpha, Default_Cl_alpha, Default_alpha0, ErrStat, ErrMsg)
      REAL(ReKi),               intent(in   ) :: alpha(:)                   ! alpha
      REAL(ReKi),               intent(in   ) :: Cn(:)                      ! cn
      REAL(ReKi),               intent(in   ) :: Cl(:)                      ! cl
   
      REAL(ReKi),               intent(  out) :: Default_Cn_alpha
      REAL(ReKi),               intent(  out) :: Default_Cl_alpha
      REAL(ReKi),               intent(  out) :: Default_alpha0
      integer(IntKi),           intent(  out) :: errStat                    ! Error status of the operation
      character(*),             intent(  out) :: errMsg                     ! Error message if ErrStat /= ErrID_None 
      
      REAL(ReKi)                              :: A(      size(alpha), 2)
      REAL(ReKi)                              :: B(max(2,size(alpha)),2)

      if (SIZE(Cn) < 2 .OR. SIZE(Cl) < 2) then
         ErrMsg='Calculate_C_alpha: Not enough data points to compute Cn and Cl slopes.'
         ErrStat=ErrID_Fatal
         Default_Cn_alpha = EPSILON(Default_Cn_alpha)
         Default_Cl_alpha = EPSILON(Default_Cl_alpha)
         Default_alpha0 = 0.0_ReKi
         return
      end if

      A(:,1) = alpha
      A(:,2) = 1.0_ReKi
      
      if (size(Cn) == 1) then
         B(:,1) = Cn(1)
         B(:,2) = Cl(1)
      else
         B(:,1) = Cn
         B(:,2) = Cl
      end if
      
      CALL LAPACK_gels('N', A, B, ErrStat, ErrMsg)
   
      Default_Cn_alpha = B(1,1)
      Default_Cl_alpha = B(1,2)
      
      if (.not. EqualRealNos(B(1,1),0.0_ReKi)) then
         Default_alpha0  = -B(2,1)/B(1,1) ! using the values from Cn_alpha
      else
         Default_alpha0 = 0.0_ReKi
      end if
         
   END SUBROUTINE Calculate_C_alpha
!----------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE ComputeUASeparationFunction_onCl(p, ColCl, ColUAf, col_fs, col_fa)
      TYPE (AFI_Table_Type),    intent(inout) :: p                             ! This structure stores all the module parameters that are set by AirfoilInfo during the initialization phase.
      integer(IntKi),           intent(in   ) :: ColUAf                        ! column for UA f_st (based on Cl or cn)
      integer(IntKi),           intent(in   ) :: ColCl                         ! column for cl
      INTEGER(IntKi),           intent(in   ) :: col_fs                        ! column for UA cn/cl_fs (fully separated cn or cl)
      INTEGER(IntKi),           intent(in   ) :: col_fa                        ! column for UA cn/cl_fa (fully attached cn or cl); NOT USED IN THE MODELS ! note that col_fa is not used in this model, but we set the values to ensure files get written properly

      
      integer                                 :: Row
      REAL(ReKi)                              :: cl_ratio
      REAL(ReKi)                              :: cl_inv
      REAL(ReKi)                              :: f_st
      REAL(ReKi)                              :: fullySeparate
      
         !------------------------------------------------
         ! calculate f_st, cl_fs, and cl_fa for HGM model
         !------------------------------------------------
         if (EqualRealNos(p%UA_BL%c_lalpha,0.0_ReKi)) then
            p%Coefs(:,ColUAf) = 0.0_ReKi                           ! Eq. 59
            p%Coefs(:,col_fs) = p%Coefs(:,ColCl)                   ! Eq. 61
            p%Coefs(:,col_fa) = 0.0_ReKi
            call ComputeUASeparationFunction_zero(p, ColUAf, p%Coefs(:,ColCl)) ! just to initialize these values... UA will turn off without using them
         else
            
               do Row=1,p%NumAlf
            
                  if (EqualRealNos( p%alpha(Row), p%UA_BL%alpha0)) then
                     f_st  = 1.0_ReKi                                         ! Eq. 59
                     fullySeparate = p%Coefs(Row,ColCl) / 2.0_ReKi            ! Eq. 61 (which should be very close to 0 because definition of alpha0 says cl(alpha0) = 0 )
                  else
            
                     cl_ratio = p%Coefs(Row,ColCl) / ( p%UA_BL%c_lalpha*(p%alpha(Row) - p%UA_BL%alpha0))
                     cl_ratio = max(0.0_ReKi, cl_ratio)
            
                     f_st = ( 2.0_ReKi * sqrt(cl_ratio) - 1.0_ReKi )**2
                  
                     if (f_st < 1.0_ReKi) then 
                        ! Region where f_st<1, merge
                        f_st  = max(0.0_ReKi, f_st) ! make sure it is not negative
                        fullySeparate = (p%Coefs(Row,ColCl) - p%UA_BL%c_lalpha* (p%alpha(Row) - p%UA_BL%alpha0)*f_st) / (1.0_ReKi - f_st) ! Eq 61
                     else
                        ! Initialize to linear region (in fact only at singularity, where f_st=1)
                        f_st = 1.0_ReKi
                        fullySeparate = p%Coefs(Row,ColCl) / 2.0_ReKi                      ! Eq. 61
                     end if
                     
                  end if
               
                  p%Coefs(Row,ColUAf) = f_st
                  p%Coefs(Row,col_fs) = fullySeparate
                  p%Coefs(Row,col_fa) = p%UA_BL%c_lalpha * (p%alpha(Row) - p%UA_BL%alpha0) ! not used in the UA model (it's specified directly), but computed here for completeness

               end do

               ! These variables aren't used with the models that use Cl instead of Cn, but it's a way to initialize the values.
               ! They make sure that the separation function is monotonic before p%UA_BL%alphaLower and after p%UA_BL%alphaUpper:
               call ComputeUASeparationFunction_zero(p, ColUAf, p%Coefs(:,ColCl)) ! this was comparing with alpha0, but now we compare with alphaUpper and alphaLower

            
               ! Ensuring everything is in harmony 
               do Row=1,p%NumAlf
                  fullySeparate = p%Coefs(Row,col_fs)
               
                  cl_inv = p%UA_BL%c_lalpha*(p%alpha(Row) - p%UA_BL%alpha0)     ! Eq. 64
                  if (.not. EqualRealNos(cl_inv, fullySeparate)) then
                     f_st=(p%Coefs(Row,ColCl) - fullySeparate) / (cl_inv - fullySeparate);        ! Eq. 60
                     f_st = max(0.0_ReKi, f_st)
                     f_st = min(1.0_ReKi, f_st)
                  
                     p%Coefs(Row,ColUAf) = f_st
                  else 
                     p%Coefs(Row,ColUAf) = 1.0_ReKi
                  end if
               end do
               
               
            end if ! c_lalpha == 0

   END SUBROUTINE ComputeUASeparationFunction_onCl
!----------------------------------------------------------------------------------------------------------------------------------  
   SUBROUTINE Compute_iLoweriUpper(p, iLower, iUpper)
      TYPE (AFI_Table_Type),    intent(in   ) :: p                             ! This structure stores all the module parameters that are set by AirfoilInfo during the initialization phase.
      INTEGER(IntKi)          , intent(  out) :: iLower                        ! The lower index separating the region around 0
      INTEGER(IntKi)          , intent(  out) :: iUpper                        ! The upper index separating the region around 0
      
      !------------------------------------------------
      ! get bounds
      !------------------------------------------------
      iLower = minloc( p%alpha , DIM=1, MASK=p%alpha >= p%UA_BL%alphaLower)
      iUpper = maxloc( p%alpha , DIM=1, MASK=p%alpha <= p%UA_BL%alphaUpper)
      
      iLower = max(1, min(p%NumAlf-1,iLower)) ! 1 <= iLower <= NumAlf-1
      iUpper = max(2, min(p%NumAlf  ,iUpper)) ! 2 <= iUpper <= NumAlf

   END SUBROUTINE Compute_iLoweriUpper
!----------------------------------------------------------------------------------------------------------------------------------  
   SUBROUTINE ComputeUASeparationFunction_zero(p, ColUAf, cn_cl)
      TYPE (AFI_Table_Type),    intent(inout) :: p                             ! This structure stores all the module parameters that are set by AirfoilInfo during the initialization phase.
      integer(IntKi),           intent(in   ) :: ColUAf                        ! column for UA f_st (based on Cl or cn)
      REAL(ReKi),               intent(in   ) :: cn_cl(:)                      ! cn or cl, whichever variable we are computing this on
   
      REAL(ReKi)                              :: c_RateBreak                   ! the slope of the wrap-around region
      INTEGER(IntKi)                          :: Row                           ! The row of a table to be parsed in the FileInfo structure.
      INTEGER(IntKi)                          :: col_fs                        ! column for UA cn/cl_fs (fully separated cn or cl)
      INTEGER(IntKi)                          :: col_fa                        ! column for UA cn/cl_fa (fully attached cn or cl)
      INTEGER(IntKi)                          :: iHigh, iLow
      INTEGER(IntKi)                          :: iTemp
      
      !------------------------------------------------
      ! set column numbers
      !------------------------------------------------
      col_fs = ColUAf + 1
      col_fa = col_fs + 1
      
         ! initialize so that we can find the minimum f on each side of the attached region
     !iLow  = minloc(p%Coefs(:,ColUAf), DIM=1, MASK=p%alpha < p%UA_BL%alphaLower, BACK=.TRUE.) ! because not all compilers allow keyword "BACK" from the F2008 standard, we implement this way:
      iTemp  = minloc(p%Coefs(:,ColUAf), DIM=1, MASK=p%alpha < p%UA_BL%alphaLower) ! because not all compilers (gcc) allow keyword "BACK" from the F2008 standard, we implement this way
      iLow  = maxloc( p%alpha, DIM=1, MASK=p%alpha < p%UA_BL%alphaLower .and. p%Coefs(:,ColUAf) == p%Coefs(iTemp,ColUAf) )

      iHigh = minloc(p%Coefs(:,ColUAf), DIM=1, MASK=p%alpha > p%UA_BL%alphaUpper)

      ! Compute variables to help x3 state with +/-180-degree wrap-around issues
      p%UA_BL%alphaBreakUpper  = p%alpha(iHigh)
      p%UA_BL%alphaBreakLower  = p%alpha(iLow)
      p%UA_BL%CnBreakUpper     = p%Coefs(iHigh,col_fa)
      p%UA_BL%CnBreakLower     = p%Coefs(iLow,col_fa)
      
      c_RateBreak       = (p%UA_BL%CnBreakUpper - p%UA_BL%CnBreakLower) / ( (p%UA_BL%alphaBreakUpper-TwoPi) - p%UA_BL%alphaBreakLower)
      
         ! make sure that the separation function is monotonic before iLow and after iHigh:
      do Row=1,iLow
         p%Coefs(Row,col_fa) = (p%alpha(Row) - p%UA_BL%alphaBreakLower) * c_RateBreak + p%UA_BL%CnBreakLower
         p%Coefs(Row,col_fs) = cn_cl(Row)
         p%Coefs(Row,ColUAf) = 0.0_ReKi
      end do
      do Row=iHigh,p%NumAlf
         p%Coefs(Row,col_fa) = (p%alpha(Row) - p%UA_BL%alphaBreakUpper) * c_RateBreak + p%UA_BL%CnBreakUpper
         p%Coefs(Row,col_fs) = cn_cl(Row)
         p%Coefs(Row,ColUAf) = 0.0_ReKi
      end do
      
   END SUBROUTINE ComputeUASeparationFunction_zero
!----------------------------------------------------------------------------------------------------------------------------------  
   SUBROUTINE ComputeUA360_AttachedFlow(p, ColUAf, cn_cl, iLower, iUpper)
      TYPE (AFI_Table_Type),    intent(inout) :: p                             ! This structure stores all the module parameters that are set by AirfoilInfo during the initialization phase.
      integer(IntKi),           intent(in   ) :: ColUAf                        ! column for UA f_st (based on Cl or cn)
      REAL(ReKi),               intent(in   ) :: cn_cl(:)                      ! cn or cl, whichever variable we are computing this on
      INTEGER(IntKi)          , intent(  out) :: iLower                        ! The lower index separating the region around 0
      INTEGER(IntKi)          , intent(  out) :: iUpper                        ! The upper index separating the region around 0
   
      REAL(ReKi)                              :: roots(p%NumAlf)
      REAL(ReKi)                              :: x_(3), f_(3)
      
      REAL(ReKi)                              :: CnSlopeUpper, alpha0Upper
      REAL(ReKi)                              :: CnSlopeLower, alpha0Lower
      REAL(ReKi)                              :: CnSlopeReverseFlow           ! Cn slope versus angle of attack for reverse flow, 1/rad

      
      INTEGER(IntKi)                          :: Row                           ! The row of a table to be parsed in the FileInfo structure.
      INTEGER(IntKi)                          :: iRoot
      INTEGER(IntKi)                          :: col_fa                        ! column for UA cn/cl_fa (fully attached cn or cl)
      INTEGER(IntKi)                          :: Indx
      INTEGER(IntKi)                          :: nZeros
      
      
      !------------------------------------------------
      ! set column numbers
      !------------------------------------------------
      col_fa = ColUAf + 2

      !------------------------------------------------
      ! get bounds
      !------------------------------------------------
      call Compute_iLoweriUpper(p, iLower, iUpper)

      p%UA_BL%alphaLower = p%alpha(iLower) ! note we are overwritting values here to make them consistent in the linear equation
      p%UA_BL%alphaUpper = p%alpha(iUpper) ! note we are overwritting values here to make them consistent in the linear equation
      
      p%UA_BL%c_alphaLower = cn_cl(iLower) ! for vortex calculations (x5, HGMV model)
      p%UA_BL%c_alphaUpper = cn_cl(iUpper) ! for vortex calculations (x5, HGMV model)
      !------------------------------------------------
      ! From dynamicStallLUT.m/updateCnAttached()
      !------------------------------------------------
      CnSlopeUpper = ( cn_cl(iUpper-1) - cn_cl(iUpper) ) / ( p%alpha(iUpper-1) - p%alpha(iUpper) )
      if (EqualRealNos(CnSlopeUpper, 0.0_ReKi)) then
         alpha0Upper = p%alpha(iUpper)
      else
         alpha0Upper  = p%alpha(iUpper) - cn_cl(iUpper)/CnSlopeUpper;
      end if
      
      CnSlopeLower = ( cn_cl(iLower) - cn_cl(iLower+1) ) / ( p%alpha(iLower) - p%alpha(iLower+1) )
      if (EqualRealNos(CnSlopeLower, 0.0_ReKi)) then
         alpha0Lower = p%alpha(iLower)
      else
         alpha0Lower  = p%alpha(iLower) - cn_cl(iLower)/CnSlopeLower;
      end if
      
      ! Find reverse flow Cn = 0 near positive 180 deg (and not in the range (- 45, 45) degrees)
      call fZeros(p%alpha, cn_cl, roots, nZeros, Period=TwoPi)
      p%UA_BL%alpha0ReverseFlow = p%alpha(1) !  default value, in case there aren't any roots. Maybe this should be an error?
      if (nZeros > 0) then
         iRoot = maxloc( abs(roots(1:nZeros)), DIM=1, MASK=abs(roots(1:nZeros)) >= 45.0_ReKi*D2R )
         if (iRoot > 0) then
            p%UA_BL%alpha0ReverseFlow = roots(iRoot)
            if (p%UA_BL%alpha0ReverseFlow < -PiBy2) p%UA_BL%alpha0ReverseFlow = p%UA_BL%alpha0ReverseFlow + TwoPi !bjj check this value along with alphaBreakLower subtracting the TwoPi
         end if
      end if
      CnSlopeReverseFlow = -TwoPi;

      
      ! Find intersections
      p%UA_BL%alphaBreakUpper = ( CnSlopeReverseFlow *  p%UA_BL%alpha0ReverseFlow          - CnSlopeUpper*alpha0Upper ) / ( CnSlopeReverseFlow - CnSlopeUpper );
      p%UA_BL%CnBreakUpper    = CnSlopeUpper*( p%UA_BL%alphaBreakUpper - alpha0Upper );
            
      p%UA_BL%alphaBreakLower = ( CnSlopeReverseFlow * (p%UA_BL%alpha0ReverseFlow - TwoPi) - CnSlopeLower*alpha0Lower ) / ( CnSlopeReverseFlow - CnSlopeLower );
      p%UA_BL%CnBreakLower    = CnSlopeLower*( p%UA_BL%alphaBreakLower - alpha0Lower );

      ! set fully attached values:
      Indx = 1
      x_ = (/ p%UA_BL%alpha0ReverseFlow-TwoPi, p%UA_BL%alphaBreakLower, p%alpha(iLower) /)
      f_ = (/ 0.0_ReKi,                        p%UA_BL%CnBreakLower,    cn_cl(iLower) /)
      do Row=1,iLower-1
         p%Coefs(Row,col_fa) = InterpExtrapStp(p%alpha(Row), x_, f_, Indx, size(x_))
      end do
      
      do Row=iLower,iUpper
         p%Coefs(Row,col_fa) = cn_cl(Row)
      end do
      
      x_ = (/ p%alpha(iUpper), p%UA_BL%alphaBreakUpper, p%UA_BL%alpha0ReverseFlow /)
      f_ = (/ cn_cl(iUpper)  , p%UA_BL%CnBreakUpper,    0.0_ReKi /)
      do Row=iUpper+1,p%NumAlf
         p%Coefs(Row,col_fa) = InterpExtrapStp(p%alpha(Row), x_, f_, Indx, size(x_))
      end do
      
   END SUBROUTINE ComputeUA360_AttachedFlow
!----------------------------------------------------------------------------------------------------------------------------------  
   SUBROUTINE ComputeUA360_updateSeparationF( p, ColUAf, cn_cl, iLower, iUpper )
      TYPE (AFI_Table_Type),    intent(inout) :: p                             ! This structure stores all the module parameters that are set by AirfoilInfo during the initialization phase.
      integer(IntKi),           intent(in   ) :: ColUAf                        ! column for UA f_st (based on Cl or cn)
      REAL(ReKi),               intent(in   ) :: cn_cl(:)                      ! cn or cl, whichever variable we are computing this on
      INTEGER(IntKi)          , intent(in   ) :: iLower                        ! The lower index separating the region around 0
      INTEGER(IntKi)          , intent(in   ) :: iUpper                        ! The upper index separating the region around 0
   
      REAL(ReKi)                              :: Offset
      REAL(ReKi)                              :: CnRatio
      REAL(ReKi)                              :: alpha_(p%NumAlf)              ! temporary for calculating periodic f_st
      REAL(ReKi)                              :: f_st(  p%NumAlf)              ! temporary for calculating periodic f_st

      INTEGER(IntKi)                          :: Row                           ! The row of a table to be parsed in the FileInfo structure.
      INTEGER(IntKi)                          :: col_fa                        ! column for UA cn/cl_fa (fully attached cn or cl)
      INTEGER(IntKi)                          :: iReverseFlow                  ! The index where f_st is at a local max near +/-180
      INTEGER(IntKi)                          :: iUpperBreak                   ! The upper index separating the region around +/-180
      INTEGER(IntKi)                          :: iLowerBreak                   ! The lower index separating the region around +/-180
      
      
      !------------------------------------------------
      ! set column numbers
      !------------------------------------------------
      col_fa = ColUAf + 2 ! fully attached (column values computed in ComputeUA360_AttachedFlow())
      
      ! compute f_st (separation function, f = p%Coefs(Row,ColUAf))
      do Row=1,p%NumAlf
         offset = ComputeUA360_CnOffset(p, cn_cl, Row, iLower)
         if (EqualRealNos(p%Coefs(Row,col_fa),offset)) then
            CnRatio = 1.0_ReKi
         else
            CnRatio = (cn_cl(Row)-offset) / (p%Coefs(Row,col_fa)-offset);  ! offset needed to ensure numerator and denomonator have same sign since sqrt is used next
         end if
         CnRatio = max( 0.25_ReKi, CnRatio ); ! below 1/4 we assume full separation and f = 0

         p%Coefs(Row,ColUAf) = ( 2.0_ReKi * sqrt( CnRatio ) - 1.0_ReKi )**2
            
         p%Coefs(Row,ColUAf) = min( p%Coefs(Row,ColUAf), 1.0_ReKi )  ! f <= 1
         p%Coefs(Row,ColUAf) = max( 0.0_ReKi, p%Coefs(Row,ColUAf) )  ! f >= 0

         !if (EqualRealNos( p%Coefs(Row,col_fa), cn_cl(Row)) p%Coefs(Row,ColUAf) = 1.0_ReKi ! Set this below without EqualRealNos()
      end do
      
         ! Where p%Coefs(Row,col_fa) == cn_cl(Row), set f = 1
      do Row=iLower,iUpper
         p%Coefs(Row,ColUAf) = 1.0_ReKi 
      end do

      !-----------------------------------------------------------
      ! now fix issues if there is a second peak near 180 degrees:
      !-----------------------------------------------------------
      iLowerBreak = maxloc( p%alpha , DIM=1, MASK=p%alpha <= p%UA_BL%alphaBreakLower)
      alpha_ = cshift(p%alpha,iLowerBreak)
      f_st   = cshift(p%Coefs(:,ColUAf),iLowerBreak)
      do Row = 2,p%NumAlf
         if (alpha_(Row) < alpha_(Row-1)) alpha_(Row) = alpha_(Row)+TwoPi
      end do
      
      iReverseFlow = maxloc( f_st, DIM=1, MASK= alpha_ > p%UA_BL%alphaBreakUpper )
      iUpperBreak = minloc( alpha_ , DIM=1, MASK=alpha_ >= p%UA_BL%alphaBreakUpper)

      ! make sure this is monotonically decreasing from a single peak:
      do Row=iReverseFlow-1,iUpperBreak+1,-1
!        if ( f_st(Row-1) > f_st(Row) )    f_st(Row-1) = max(0.0_ReKi, f_st(Row) - ABS( (f_st(Row+1) - f_st(Row) )/(alpha_(Row+1) - alpha_(Row)) * (alpha_(Row)-alpha_(Row-1))))
         if (EqualRealNos(f_st(Row),0.0_ReKi)) f_st(Row-1) = 0.0_ReKi
         if ( f_st(Row-1) > f_st(Row) )    f_st(Row) = 0.5_ReKi * (f_st(Row+1) + f_st(Row-1))
      end do
      do Row=iReverseFlow+1,p%NumAlf-1
!        if ( f_st(Row+1) > f_st(Row) )    f_st(Row+1) = max(0.0_ReKi, f_st(Row) - ABS( (f_st(Row-1) - f_st(Row) )/(alpha_(Row-1) - alpha_(Row)) * (alpha_(Row+1) - alpha_(Row))))
         if (EqualRealNos(f_st(Row),0.0_ReKi)) f_st(Row+1) = 0.0_ReKi
         if ( f_st(Row+1) > f_st(Row) )    f_st(Row) = 0.5_ReKi * (f_st(Row+1) + f_st(Row-1))
      end do

      p%Coefs(:,ColUAf)   = cshift(f_st,-iLowerBreak)

      
   END SUBROUTINE ComputeUA360_updateSeparationF
!----------------------------------------------------------------------------------------------------------------------------------  
   SUBROUTINE ComputeUA360_updateCnSeparated( p, ColUAf, cn_cl, iLower )
      TYPE (AFI_Table_Type),    intent(inout) :: p                             ! This structure stores all the module parameters that are set by AirfoilInfo during the initialization phase.
      integer(IntKi),           intent(in   ) :: ColUAf                        ! column for UA f_st (based on Cl or cn)
      REAL(ReKi),               intent(in   ) :: cn_cl(:)                      ! cn or cl, whichever variable we are computing this on
      INTEGER(IntKi)          , intent(in   ) :: iLower                        ! The lower index separating the region around 0
   
      REAL(ReKi)                              :: Offset                       
      INTEGER(IntKi)                          :: Row                           ! The row of a table to be parsed in the FileInfo structure.
      INTEGER(IntKi)                          :: col_fa                        ! column for UA cn/cl_fa (fully attached cn or cl)
      INTEGER(IntKi)                          :: col_fs                        ! column for UA cn/cl_fa (fully separated cn or cl)
      
      !------------------------------------------------
      ! set column numbers
      !------------------------------------------------
      col_fa = ColUAf + 2 ! fully attached
      col_fs = ColUAf + 1 ! fully separate

      do Row=1,p%NumAlf
         if (EqualRealNos( p%Coefs(Row,ColUAf), 1.0_ReKi )) then
            offset = ComputeUA360_CnOffset(p, cn_cl, Row, iLower)
            p%Coefs(Row,col_fs) = 0.5_ReKi * (cn_cl(Row) + offset)
         else
            p%Coefs(Row,col_fs) = ( cn_cl(Row) - p%Coefs(Row,col_fa) * p%Coefs(Row,ColUAf) ) / ( 1.0_ReKi - p%Coefs(Row,ColUAf) )
         end if
      end do

   END SUBROUTINE ComputeUA360_updateCnSeparated
!----------------------------------------------------------------------------------------------------------------------------------  
   REAL(ReKi) FUNCTION ComputeUA360_CnOffset(p, cn_cl, Row, iLower) RESULT(offset)
      TYPE (AFI_Table_Type),    intent(in   ) :: p                             ! This structure stores all the module parameters that are set by AirfoilInfo during the initialization phase.
      REAL(ReKi),               intent(in   ) :: cn_cl(:)                      ! cn or cl, whichever variable we are computing this on
      INTEGER(IntKi)          , intent(in   ) :: Row                           ! The row of a table to be parsed in the FileInfo structure.
      INTEGER(IntKi)          , intent(in   ) :: iLower                        ! The lower index separating the region around 0
   
      REAL(ReKi)                              :: CnOffset                     ! Mathematical trick: offset to Cn making formulation of f-separation behave for strange polars with negative stall at positive Cn values (usually soiled polars for thick airfoils)
      REAL(ReKi)                              :: SlopeScale
         
   
      ! compute cnOffset
      if (cn_cl(iLower) > -0.05) then
         CnOffset = cn_cl(iLower) + 0.05
      else
         CnOffset = 0.0_ReKi
      end if
      
      SlopeScale = 0.1_ReKi*R2D
      offset = CnOffset * ( tanh(SlopeScale*(p%alpha(Row)+PiBy2)) - tanh(SlopeScale*(p%alpha(Row)-PiBy2)) ) / 2.0_ReKi; !Only apply Cn offset in vicinity of AoA 0 deg
   END FUNCTION ComputeUA360_CnOffset
!----------------------------------------------------------------------------------------------------------------------------------  
subroutine FindBoundingTables(p, secondaryDepVal, lowerTable, upperTable, xVals)

   TYPE (AFI_ParameterType), intent(in   ) :: p                          ! This structure stores all the module parameters that are set by AirfoilInfo during the initialization phase.
   real(ReKi),               intent(in   ) :: secondaryDepVal            ! Interpolate secondary values (Re or UserProp) based on this value

   integer(IntKi),           intent(  out) :: lowerTable                 ! index of the lower airfoil table p%secondVals(lowerTable) <= secondaryDepVal <= p%secondVals(upperTable)
   integer(IntKi),           intent(  out) :: upperTable                 ! index of the upper airfoil table ( = lowerTable+1 )
   real(ReKi),               intent(  out) :: xVals(2)                   ! this takes the place of time in the extrapInterp routines generated by the Registry
   
   integer(IntKi)                          :: iMid

   ! compare algorithm with InterpBinReal
   lowerTable  = 1
   upperTable  = p%NumTabs

   DO WHILE ( upperTable-lowerTable > 1 )

      iMid = ( upperTable + lowerTable )/2

      IF ( secondaryDepVal >= p%secondVals(iMid) ) THEN
         lowerTable = iMid
      ELSE
         upperTable = iMid
      END IF

   END DO
   
   xVals(1) = p%secondVals(lowerTable)
   xVals(2) = p%secondVals(upperTable)

end subroutine FindBoundingTables
!----------------------------------------------------------------------------------------------------------------------------------  
subroutine AFI_ComputeUACoefs2D( secondaryDepVal, p, UA_BL, errStat, errMsg )
! This routine is calculates the UA parameters for a set of tables which are dependent a 2nd user-defined varible, could be Re or Cntrl, etc.
! If the requested yVar is not associated with a given table, then the two tables which contain yVar are found and, a cubic spline interpolation is performed at the requested AOA.
! for each of those two tables. Then a linear intepolation is performed on the 2nd dimension to find the final Cl,Cd,Cm, and Cpmin values.
! If the requested yVar corresponds to a table, then only a single cubic interpolation based on the requested AOA is performed.
!..................................................................................................................................
   real(ReKi),               intent(in   ) :: secondaryDepVal            ! Interpolate based on this value
   TYPE (AFI_ParameterType), intent(in   ) :: p                          ! This structure stores all the module parameters that are set by AirfoilInfo during the initialization phase.
   type(AFI_UA_BL_Type),     intent(  out) :: UA_BL
   integer(IntKi),           intent(  out) :: errStat                    ! Error status of the operation
   character(*),             intent(  out) :: errMsg                     ! Error message if ErrStat /= ErrID_None 
   
   real(ReKi)                              :: xVals(2)                   ! secondary interpolation values associated with the tables (this takes the place of time in the extrapInterp routines generated by the Registry)
   
   integer                                 :: lowerTable, upperTable
   character(*), parameter                 :: RoutineName = 'AFI_ComputeUACoefs2D'

   
   ErrStat = ErrID_None
   ErrMsg  = ''
   
      ! find boundaries for 
   
         ! Let's check the limits first.

   IF ( secondaryDepVal <= p%secondVals( 1 ) )  THEN
   ! was call SetErrStat (ErrID_Fatal, "Specified Reynold's number, "//trim(num2lstr(secondaryDepVal))//" , is outside the range of Re specified in the airfoil input file tables.", ErrStat, ErrMsg, RoutineName )
   ! or  call SetErrStat (ErrID_Fatal, "Specified User Property's value, "//trim(num2lstr(secondaryDepVal))//" , is outside the range of User Property values specified in the airfoil input file tables.", ErrStat, ErrMsg, RoutineName )
      call AFI_CopyUA_BL_Type( p%Table(1)%UA_BL, UA_BL, MESH_NEWCOPY, errStat, errMsg )  ! this doesn't have a mesh, so the control code is irrelevant
      return
   ELSE IF ( secondaryDepVal >= p%secondVals( p%NumTabs ) ) THEN
   ! was call SetErrStat (ErrID_Fatal, "Specified Reynold's number, "//trim(num2lstr(secondaryDepVal))//" , is outside the range of Re specified in the airfoil input file tables.", ErrStat, ErrMsg, RoutineName )
   ! or  call SetErrStat (ErrID_Fatal, "Specified User Property's value, "//trim(num2lstr(secondaryDepVal))//" , is outside the range of User Property values specified in the airfoil input file tables.", ErrStat, ErrMsg, RoutineName )
      call AFI_CopyUA_BL_Type( p%Table(p%NumTabs)%UA_BL, UA_BL, MESH_NEWCOPY, errStat, errMsg )  ! this doesn't have a mesh, so the control code is irrelevant
      return
   END IF

   call FindBoundingTables(p, secondaryDepVal, lowerTable, upperTable, xVals)

      ! linearly interpolate
   call AFI_UA_BL_Type_ExtrapInterp1(p%Table(lowerTable)%UA_BL, p%Table(upperTable)%UA_BL, xVals, UA_BL, secondaryDepVal, ErrStat, ErrMsg )

   
end subroutine AFI_ComputeUACoefs2D  
!----------------------------------------------------------------------------------------------------------------------------------  
subroutine AFI_ComputeAirfoilCoefs2D( AOA, secondaryDepVal, p, AFI_interp, errStat, errMsg )
! This routine is calculates Cl, Cd, Cm, (and Cpmin) for a set of tables which are dependent on AOA as well as a 2nd user-defined varible, could be Re or Cntrl, etc.
! If the requested yVar is not associated with a given table, then the two tables which contain yVar are found and, a cubic spline interpolation is performed at the requested AOA.
! for each of those two tables. Then a linear intepolation is performed on the 2nd dimension to find the final Cl,Cd,Cm, and Cpmin values.
! If the requested yVar corresponds to a table, then only a single cubic interpolation based on the requested AOA is performed.
!..................................................................................................................................
   real(ReKi),               intent(in   ) :: AOA
   real(ReKi),               intent(in   ) :: secondaryDepVal           ! Unused in the current version!     
   TYPE (AFI_ParameterType), intent(in   ) :: p                          ! This structure stores all the module parameters that are set by AirfoilInfo during the initialization phase.
   type(AFI_OutputType),     intent(  out) :: AFI_interp                 ! contains    real(ReKi),               intent(  out) :: Cl, Cd, Cm, Cpmin
   integer(IntKi),           intent(  out) :: errStat                    ! Error status of the operation
   character(*),             intent(  out) :: errMsg                     ! Error message if ErrStat /= ErrID_None 
   
   
   integer                                 :: lowerTable, upperTable
   real(ReKi)                              :: xVals(2)
   type(AFI_OutputType)                    :: AFI_lower
   type(AFI_OutputType)                    :: AFI_upper
      
   ErrStat = ErrID_None
   ErrMsg  = ''
   
   IF ( secondaryDepVal <= p%secondVals( 1 ) )  THEN
   ! was call SetErrStat (ErrID_Fatal, "Specified Reynold's number, "//trim(num2lstr(secondaryDepVal))//" , is outside the range of Re specified in the airfoil input file tables.", ErrStat, ErrMsg, RoutineName )
   ! or  call SetErrStat (ErrID_Fatal, "Specified User Property's value, "//trim(num2lstr(secondaryDepVal))//" , is outside the range of User Property values specified in the airfoil input file tables.", ErrStat, ErrMsg, RoutineName )
      call AFI_ComputeAirfoilCoefs1D( AOA, p, AFI_interp, errStat, errMsg, 1 )
      return
   ELSE IF ( secondaryDepVal >= p%secondVals( p%NumTabs ) ) THEN
   ! was call SetErrStat (ErrID_Fatal, "Specified Reynold's number, "//trim(num2lstr(secondaryDepVal))//" , is outside the range of Re specified in the airfoil input file tables.", ErrStat, ErrMsg, RoutineName )
   ! or  call SetErrStat (ErrID_Fatal, "Specified User Property's value, "//trim(num2lstr(secondaryDepVal))//" , is outside the range of User Property values specified in the airfoil input file tables.", ErrStat, ErrMsg, RoutineName )
      call AFI_ComputeAirfoilCoefs1D( AOA, p, AFI_interp, errStat, errMsg, p%NumTabs )
      return
   END IF
   
   call FindBoundingTables(p, secondaryDepVal, lowerTable, upperTable, xVals)
   
!fixme ERROR HANDLING!   
   call AFI_ComputeAirfoilCoefs1D( AOA, p, AFI_lower, errStat, errMsg, lowerTable )
      if (ErrStat >= AbortErrLev) return
   call AFI_ComputeAirfoilCoefs1D( AOA, p, AFI_upper, errStat, errMsg, upperTable )
      if (ErrStat >= AbortErrLev) return

       ! linearly interpolate these values
   call AFI_Output_ExtrapInterp1(AFI_lower, AFI_upper, xVals, AFI_interp, secondaryDepVal, ErrStat, ErrMsg )
   
   
end subroutine AFI_ComputeAirfoilCoefs2D  
         
!----------------------------------------------------------------------------------------------------------------------------------  
subroutine AFI_ComputeAirfoilCoefs1D( AOA, p, AFI_interp, errStat, errMsg, TableNum )
! If the requested yVar is not associated with a given table, then the two tables which contain yVar are found and, a cubic spline interpolation is performed at the requested AOA.
! for each of those two tables. Then a linear intepolation is performed on the 2nd dimension to find the final Cl,Cd,Cm, and Cpmin values.
! If the requested yVar corresponds to a table, then only a single cubic interpolation based on the requested AOA is performed.
!..................................................................................................................................
   real(ReKi),               intent(in   ) :: AOA 
   TYPE (AFI_ParameterType), intent(in   ) :: p                          ! This structure stores all the module parameters that are set by AirfoilInfo during the initialization phase.
   type(AFI_OutputType)                    :: AFI_interp                 !  Cl, Cd, Cm, Cpmin
   integer(IntKi),           intent(  out) :: errStat                    ! Error status of the operation
   character(*),             intent(  out) :: errMsg                     ! Error message if ErrStat /= ErrID_None
   integer(IntKi), optional, intent(in   ) :: TableNum
   
   
   real                                    :: IntAFCoefs(MaxNumAFCoeffs)                ! The interpolated airfoil coefficients.
   real(reki)                              :: Alpha
   integer                                 :: s1
   integer                                 :: iTab

      
   ErrStat = ErrID_None
   ErrMsg  = ''

   if (present(TableNum)) then
      iTab = TableNum
   else
      iTab = 1
   end if
   
   IntAFCoefs = 0.0_ReKi ! initialize in case we only don't have MaxNumAFCoeffs columns in the airfoil data (e.g., so cm is zero if not in the file)
 
   s1 = size(p%Table(iTab)%Coefs,2)
   
   if (p%Table(iTab)%ConstData) then
      IntAFCoefs(1:s1) = p%Table(iTab)%Coefs(1,:)   ! all the rows are constant, so we can just return the values at any alpha (e.g., row 1)
   else
      Alpha = AOA
      call MPi2Pi ( Alpha ) ! change AOA into range of -pi to pi
   
   
         ! Spline interpolation of lower table based on requested AOA
   
      IntAFCoefs(1:s1) = CubicSplineInterpM( Alpha  &
                                             , p%Table(iTab)%Alpha &
                                             , p%Table(iTab)%Coefs &
                                             , p%Table(iTab)%SplineCoefs &
                                             , ErrStat, ErrMsg )
   end if
  
   AFI_interp%Cl    = IntAFCoefs(p%ColCl)
   AFI_interp%Cd    = IntAFCoefs(p%ColCd)
     
   if ( p%ColCm > 0 ) then
      AFI_interp%Cm = IntAFCoefs(p%ColCm)
   else
      AFI_interp%Cm    = 0.0_Reki  !Set these to zero unless there is data to be read in
   end if
   
   if ( p%ColCpmin > 0 ) then
      AFI_interp%Cpmin = IntAFCoefs(p%ColCpmin)
   else
      AFI_interp%Cpmin = 0.0_Reki
   end if

   if ( p%ColUAf > 0 ) then
      AFI_interp%f_st          = IntAFCoefs(p%ColUAf)   ! separation function
      AFI_interp%fullySeparate = IntAFCoefs(p%ColUAf+1) ! fully separated cn or cl
      AFI_interp%fullyAttached = IntAFCoefs(p%ColUAf+2) ! fully attached cn or cl
   else
      AFI_interp%f_st          = 0.0_ReKi
      AFI_interp%fullySeparate = 0.0_ReKi
      AFI_interp%fullyAttached = 0.0_ReKi
   end if
   
      ! needed if using UnsteadyAero:
   if (p%Table(iTab)%InclUAdata) then
      AFI_interp%Cd0 = p%Table(iTab)%UA_BL%Cd0
      AFI_interp%Cm0 = p%Table(iTab)%UA_BL%Cm0
   else
      AFI_interp%Cd0 = 0.0_ReKi
      AFI_interp%Cm0 = 0.0_ReKi
   end if
   
   
end subroutine AFI_ComputeAirfoilCoefs1D

!----------------------------------------------------------------------------------------------------------------------------------  
!> This routine calculates Cl, Cd, Cm, (and Cpmin) for a set of tables which are dependent on AOA as well as a 2nd user-defined varible, could be Re or Cntrl, etc.
subroutine AFI_ComputeAirfoilCoefs( AOA, Re, UserProp, p, AFI_interp, errStat, errMsg )

   real(ReKi),               intent(in   ) :: AOA
   real(ReKi),               intent(in   ) :: Re                         ! Reynold's Number
   real(ReKi),               intent(in   ) :: UserProp                   !< User property for interpolating airfoil tables
   TYPE (AFI_ParameterType), intent(in   ) :: p                          ! This structure stores all the module parameters that are set by AirfoilInfo during the initialization phase.
   type(AFI_OutputType),     intent(  out) :: AFI_interp                 ! contains   real(ReKi),               intent(  out) :: Cl, Cd, Cm, Cpmin
   integer(IntKi),           intent(  out) :: errStat                    ! Error status of the operation
   character(*),             intent(  out) :: errMsg                     ! Error message if ErrStat /= ErrID_None 

   real(ReKi)                              :: ReInterp

      ! These coefs are stored in the p data structures based on Re
   
   if ( p%AFTabMod == AFITable_1 ) then 
      call AFI_ComputeAirfoilCoefs1D( AOA, p, AFI_interp, errStat, errMsg, 1 )
   elseif ( p%AFTabMod == AFITable_2Re ) then
#ifndef AFI_USE_LINEAR_RE
      ReInterp = log( Re )
#else
      ReInterp =      Re
#endif
      call AFI_ComputeAirfoilCoefs2D( AOA, ReInterp, p, AFI_interp, errStat, errMsg )
   else !if ( p%AFTabMod == AFITable_2User ) then
      call AFI_ComputeAirfoilCoefs2D( AOA, UserProp, p, AFI_interp, errStat, errMsg )
   end if
   
   ! put some limits on the separation function:
   AFI_interp%f_st = min( max( AFI_interp%f_st, 0.0_ReKi), 1.0_ReKi)  ! separation function

   
end subroutine AFI_ComputeAirfoilCoefs

!----------------------------------------------------------------------------------------------------------------------------------  
!> This routine calculates Cl, Cd, Cm, (and Cpmin) for a set of tables which are dependent on AOA as well as a 2nd user-defined varible, could be Re or Cntrl, etc.
subroutine AFI_ComputeUACoefs( p, Re, UserProp, UA_BL, errMsg, errStat )

   type(AFI_ParameterType), intent(in   ) :: p                             !< This structure stores all the module parameters that are set by AirfoilInfo during the initialization phase.
   real(ReKi),              intent(in   ) :: Re                            !< Reynold's number
   real(ReKi),              intent(in   ) :: UserProp                      !< User property for interpolating airfoil tables
   type(AFI_UA_BL_Type),    intent(  out) :: UA_BL                         !< airfoil constants (UA Beddoes-Leishman parameters )

   integer(IntKi),          intent(  out) :: errStat                       !< Error status
   character(*),            intent(  out) :: errMsg                        !< Error message

   real(ReKi)                             :: ReInterp
   

      ! These coefs are stored in the p data structures based on Re
   
   if ( p%AFTabMod == AFITable_1 ) then 
      call AFI_CopyUA_BL_Type( p%Table(1)%UA_BL, UA_BL, MESH_NEWCOPY, errStat, errMsg )  ! this doesn't have a mesh, so the control code is irrelevant
      return
   elseif ( p%AFTabMod == AFITable_2Re ) then
#ifndef AFI_USE_LINEAR_RE
      ReInterp = log( Re )
#else
      ReInterp =      Re
#endif
      call AFI_ComputeUACoefs2D( ReInterp, p, UA_BL, errStat, errMsg )
   else !if ( p%AFTabMod == AFITable_2User ) then
      call AFI_ComputeUACoefs2D( UserProp, p, UA_BL, errStat, errMsg )
   end if
   
   call MPi2Pi( UA_BL%alpha0 )
   call MPi2Pi( UA_BL%alpha1 )
   call MPi2Pi( UA_BL%alpha2 )
   
   ! Cn1=1.9 Tp=1.7 Tf=3., Tv=6 Tvl=11, Cd0=0.012
   
end subroutine AFI_ComputeUACoefs

!=============================================================================
subroutine AFI_WrHeader(delim, FileName, unOutFile, ErrStat, ErrMsg)

   character(*),                 intent(in   )  :: delim
   character(*),                 intent(in   )  :: FileName
   integer(IntKi),               intent(  out)  :: unOutFile
   integer(IntKi),               intent(  out)  :: ErrStat     ! Error status of the operation
   character(*),                 intent(  out)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   integer(IntKi)                               :: i
   integer(IntKi)                               :: ErrStat2
   character(ErrMsgLen)                         :: ErrMsg2
   character(*), parameter                      :: RoutineName = 'AFI_WrHeader'
   
   integer, parameter                           :: MaxLen = 17
   integer, parameter                           :: NumChans = 46
   character(MaxLen)                            :: ChanName( NumChans)
   character(MaxLen)                            :: ChanUnit( NumChans)
   
   
   i=1
   ChanName(i) = 'AirfoilNumber';     ChanUnit(i) = '(-)';        i = i+1;
   ChanName(i) = 'TableNumber';       ChanUnit(i) = '(-)';        i = i+1;
   ChanName(i) = 'alpha0';            ChanUnit(i) = '(deg)';      i = i+1;
   ChanName(i) = 'alpha1';            ChanUnit(i) = '(deg)';      i = i+1;
   ChanName(i) = 'alpha2';            ChanUnit(i) = '(deg)';      i = i+1;
   ChanName(i) = 'eta_e';             ChanUnit(i) = '(-)';        i = i+1;
   ChanName(i) = 'C_nalpha';          ChanUnit(i) = '(-/rad)';    i = i+1;
   ChanName(i) = 'C_lalpha';          ChanUnit(i) = '(-/rad)';    i = i+1;
   ChanName(i) = 'T_f0';              ChanUnit(i) = '(-)';        i = i+1;
   ChanName(i) = 'T_V0';              ChanUnit(i) = '(-)';        i = i+1;
   ChanName(i) = 'T_p';               ChanUnit(i) = '(-)';        i = i+1;
   ChanName(i) = 'T_VL';              ChanUnit(i) = '(-)';        i = i+1;
   ChanName(i) = 'b1';                ChanUnit(i) = '(-)';        i = i+1;
   ChanName(i) = 'b2';                ChanUnit(i) = '(-)';        i = i+1;
   ChanName(i) = 'b5';                ChanUnit(i) = '(-)';        i = i+1;
   ChanName(i) = 'A1';                ChanUnit(i) = '(-)';        i = i+1;
   ChanName(i) = 'A2';                ChanUnit(i) = '(-)';        i = i+1;
   ChanName(i) = 'A5';                ChanUnit(i) = '(-)';        i = i+1;
   ChanName(i) = 'S1';                ChanUnit(i) = '(-)';        i = i+1;
   ChanName(i) = 'S2';                ChanUnit(i) = '(-)';        i = i+1;
   ChanName(i) = 'S3';                ChanUnit(i) = '(-)';        i = i+1;
   ChanName(i) = 'S4';                ChanUnit(i) = '(-)';        i = i+1;
   ChanName(i) = 'Cn1';               ChanUnit(i) = '(-)';        i = i+1;
   ChanName(i) = 'Cn2';               ChanUnit(i) = '(-)';        i = i+1;
   ChanName(i) = 'St_sh';             ChanUnit(i) = '(-)';        i = i+1;
   ChanName(i) = 'Cd0';               ChanUnit(i) = '(-)';        i = i+1;
   ChanName(i) = 'Cm0';               ChanUnit(i) = '(-)';        i = i+1;
   ChanName(i) = 'k0';                ChanUnit(i) = '(-)';        i = i+1;
   ChanName(i) = 'k1';                ChanUnit(i) = '(-)';        i = i+1;
   ChanName(i) = 'k2';                ChanUnit(i) = '(-)';        i = i+1;
   ChanName(i) = 'k3';                ChanUnit(i) = '(-)';        i = i+1;
   ChanName(i) = 'k1_hat';            ChanUnit(i) = '(-)';        i = i+1;
   ChanName(i) = 'x_cp_bar';          ChanUnit(i) = '(-)';        i = i+1;
   ChanName(i) = 'UACutout';          ChanUnit(i) = '(deg)';      i = i+1;
   ChanName(i) = 'UACutout_delta';    ChanUnit(i) = '(deg)';      i = i+1;
   ChanName(i) = 'UACutout_blend';    ChanUnit(i) = '(deg)';      i = i+1;
   ChanName(i) = 'filtCutOff';        ChanUnit(i) = '(-)';        i = i+1;
   ChanName(i) = 'alphaLower';        ChanUnit(i) = '(deg)';      i = i+1;
   ChanName(i) = 'alphaUpper';        ChanUnit(i) = '(deg)';      i = i+1;
   ChanName(i) = 'c_alphaLower';      ChanUnit(i) = '(-)';        i = i+1;
   ChanName(i) = 'c_alphaUpper';      ChanUnit(i) = '(-)';        i = i+1;
   ChanName(i) = 'alpha0ReverseFlow'; ChanUnit(i) = '(deg)';      i = i+1;
   ChanName(i) = 'alphaBreakUpper';   ChanUnit(i) = '(deg)';      i = i+1;
   ChanName(i) = 'CnBreakUpper';      ChanUnit(i) = '(-)';        i = i+1;
   ChanName(i) = 'alphaBreakLower';   ChanUnit(i) = '(deg)';      i = i+1;
   ChanName(i) = 'CnBreakLower';      ChanUnit(i) = '(-)';        i = i+1;

   !$OMP critical(fileopen_critical)
   CALL GetNewUnit( unOutFile, ErrStat, ErrMsg )
   if (ErrStat < AbortErrLev) then
      CALL OpenFOutFile ( unOutFile, trim(FileName), ErrStat2, ErrMsg2 )
   endif
   !$OMP end critical(fileopen_critical)
         
   ! Generate file outputs

   write (unOutFile,'(/,A)')  'Predictions were generated on '//CurDate()//' at '//CurTime() !//' using '//trim(GetNVD(version))
   write (unOutFile,'(1X,A)') trim(ProgName)
   write (unOutFile,'()' )    !print a blank line
   write (unOutFile,'()' )    !print a blank line
   write (unOutFile,'()' )    !print a blank line
              

      !......................................................
      ! Write the names of the output parameters on one line:
      !......................................................
   call WrFileNR ( unOutFile, ChanName(1) )
   do i=2,size(ChanName)
      call WrFileNR ( unOutFile, delim//ChanName(i) )
   end do
   write (unOutFile,'()')

      !......................................................
      ! Write the units of the output parameters on one line:
      !......................................................
   call WrFileNR ( unOutFile, ChanUnit(1) )
   do i=2,size(ChanName)
      call WrFileNR ( unOutFile, delim//ChanUnit(i) )
   end do
   write (unOutFile,'()')
   
end subroutine AFI_WrHeader
!=============================================================================
subroutine AFI_WrData(k, unOutFile, delim, AFInfo)
   type(AFI_ParameterType),      intent(in   )  :: AFInfo      ! The airfoil parameter data (for all airfoils)
   integer,                      intent(in   )  :: k
   integer(IntKi),               intent(in   )  :: unOutFile
   character(*),                 intent(in   )  :: delim
   
   integer(IntKi)                               :: i
   
   integer, parameter                           :: MaxLen = 17
   integer, parameter                           :: NumChans = 46
   real(ReKi)                                   :: TmpValues(NumChans)
   character(3)                                 :: MaxLenStr
   character(80)                                :: Fmt
   
   MaxLenStr = trim(num2lstr(MaxLen))
   TmpValues = 0.0_ReKi ! initialize in case UAdata is not included in the airfoil table

   Fmt = '(I'//MaxLenStr//',"'//delim//'",I'//MaxLenStr//','//trim(num2lstr(NumChans))//'("'//delim//'",F'//MaxLenStr//'.5))'
   
   do i=1,size(AFInfo%Table)
      IF (AFInfo%Table(i)%InclUAdata) then
         WRITE(unOutFile, Fmt) k, i, &
                                    AFInfo%Table(i)%UA_BL%alpha0*R2D        , &
                                    AFInfo%Table(i)%UA_BL%alpha1*R2D        , &
                                    AFInfo%Table(i)%UA_BL%alpha2*R2D        , &
                                    AFInfo%Table(i)%UA_BL%eta_e             , &
                                    AFInfo%Table(i)%UA_BL%C_nalpha          , &
                                    AFInfo%Table(i)%UA_BL%C_lalpha          , &
                                    AFInfo%Table(i)%UA_BL%T_f0              , &
                                    AFInfo%Table(i)%UA_BL%T_V0              , &
                                    AFInfo%Table(i)%UA_BL%T_p               , &
                                    AFInfo%Table(i)%UA_BL%T_VL              , &
                                    AFInfo%Table(i)%UA_BL%b1                , &
                                    AFInfo%Table(i)%UA_BL%b2                , &
                                    AFInfo%Table(i)%UA_BL%b5                , &
                                    AFInfo%Table(i)%UA_BL%A1                , &
                                    AFInfo%Table(i)%UA_BL%A2                , &
                                    AFInfo%Table(i)%UA_BL%A5                , &
                                    AFInfo%Table(i)%UA_BL%S1                , &
                                    AFInfo%Table(i)%UA_BL%S2                , &
                                    AFInfo%Table(i)%UA_BL%S3                , &
                                    AFInfo%Table(i)%UA_BL%S4                , &
                                    AFInfo%Table(i)%UA_BL%Cn1               , &
                                    AFInfo%Table(i)%UA_BL%Cn2               , &
                                    AFInfo%Table(i)%UA_BL%St_sh             , &
                                    AFInfo%Table(i)%UA_BL%Cd0               , &
                                    AFInfo%Table(i)%UA_BL%Cm0               , &
                                    AFInfo%Table(i)%UA_BL%k0                , &
                                    AFInfo%Table(i)%UA_BL%k1                , &
                                    AFInfo%Table(i)%UA_BL%k2                , &
                                    AFInfo%Table(i)%UA_BL%k3                , &
                                    AFInfo%Table(i)%UA_BL%k1_hat            , &
                                    AFInfo%Table(i)%UA_BL%x_cp_bar          , &
                                    AFInfo%Table(i)%UA_BL%UACutout*R2D      , &
                                    AFInfo%Table(i)%UA_BL%UACutout_delta*R2D, &
                                    AFInfo%Table(i)%UA_BL%UACutout_blend*R2D, &
                                    AFInfo%Table(i)%UA_BL%filtCutOff        , &
                                    AFInfo%Table(i)%UA_BL%alphaLower*R2D    , &
                                    AFInfo%Table(i)%UA_BL%alphaUpper*R2D    , &
                                    AFInfo%Table(i)%UA_BL%c_alphaLower      , &
                                    AFInfo%Table(i)%UA_BL%c_alphaUpper      , &
                                    AFInfo%Table(i)%UA_BL%alpha0ReverseFlow*R2D, &
                                    AFInfo%Table(i)%UA_BL%alphaBreakUpper*R2D, &
                                    AFInfo%Table(i)%UA_BL%CnBreakUpper       , &
                                    AFInfo%Table(i)%UA_BL%alphaBreakLower*R2D, &
                                    AFInfo%Table(i)%UA_BL%CnBreakLower

      ELSE
         WRITE(unOutFile, Fmt) k, i, TmpValues(3:)
      END IF
   end do
      
end subroutine AFI_WrData
!=============================================================================
subroutine AFI_WrTables(AFI_Params,UAMod,OutRootName)
   
   type(AFI_ParameterType), intent(in), target  :: AFI_Params
   integer(IntKi),          intent(in)          :: UAMod
   character(*),            intent(in)          :: OutRootName
      
   integer(IntKi)                               :: unOutFile
   integer(IntKi)                               :: ErrStat
   character(ErrMsgLen)                         :: ErrMsg
   
   Real(ReKi),      allocatable                 :: cl_smooth(:)
   Real(ReKi),      allocatable                 :: cn_smooth(:)
   Real(ReKi),      allocatable                 :: cn(:)
   Real(ReKi),      allocatable                 :: cc(:)
   Real(ReKi),      allocatable                 :: cl_lin(:)
   Real(ReKi),      allocatable                 :: cn_lin(:)
   
   character(len=3)                             :: Prefix
   character(len=11)                            :: sFullyAtt
   character(len=8)                             :: sCm
   integer                                      :: iTab, iRow
   type(AFI_Table_Type), pointer                :: table !< Alias
   
   if (UAMod /= UA_HGM .and. UAMod /= UA_Oye) then
      Prefix='Cn_'
      sFullyAtt='Cn_FullyAtt'
   else
      Prefix='Cl_'
      sFullyAtt='Dummy'
   endif
   if (AFI_Params%ColCm > 0) then
      sCm='Cm'
   else
      sCm='Cm_Dummy'
   endif
   
      
   ! Loop on tables, write a different file for each table.
   do iTab = 1, size(AFI_Params%Table)
      table => AFI_Params%Table(iTab)
      
      ! Compute derived parameters from cl and cd, and UA_BL
      allocate(cl_smooth(table%NumAlf))
      allocate(cn_smooth(table%NumAlf))
      allocate(cn       (table%NumAlf))
      allocate(cc       (table%NumAlf))
      allocate(cl_lin   (table%NumAlf))
      allocate(cn_lin   (table%NumAlf))
      

      cn     = table%Coefs(:,AFI_Params%ColCl) * cos(table%alpha) + (table%Coefs(:,AFI_Params%ColCd) - table%UA_BL%Cd0) * sin(table%alpha);
      cc     = table%Coefs(:,AFI_Params%ColCl) * sin(table%alpha) - (table%Coefs(:,AFI_Params%ColCd) - table%UA_BL%Cd0) * cos(table%alpha);
      cn_lin = table%UA_BL%C_nalpha * (table%alpha - table%UA_BL%alpha0)
      cl_lin = table%UA_BL%C_lalpha * (table%alpha - table%UA_BL%alpha0)

      do iRow = 1, table%NumAlf
         if ((table%alpha(iRow)<table%UA_BL%alphaBreakLower).or. table%alpha(iRow)>table%UA_BL%alphaBreakUpper) then
            cl_lin(iRow) =0.0_ReKi
            cn_lin(iRow) =0.0_ReKi
         endif
      enddo

      ! Smoothing (used priot to compute slope in CalculateUACoeffs)
      call kernelSmoothing(table%alpha, cn                             , kernelType_TRIWEIGHT, 2.0_ReKi*D2R, cn_smooth)
      call kernelSmoothing(table%alpha, table%Coefs(:,AFI_Params%ColCl), kernelType_TRIWEIGHT, 2.0_ReKi*D2R, cl_smooth)

      
      ! Write to file
      !$OMP critical(fileopen_critical)
      CALL GetNewUnit( unOutFile, ErrStat, ErrMsg )
      if (ErrStat < AbortErrLev) then
         CALL OpenFOutFile ( unOutFile, trim(OutRootName)//'.Coefs.'//trim(num2lstr(iTab))//'.out', ErrStat, ErrMsg )
      endif
      !$OMP end critical(fileopen_critical)
         if (ErrStat >= AbortErrLev) then
            call WrScr(Trim(ErrMsg))
            return
         end if
   
      WRITE (unOutFile,'(/,A/)') 'These predictions were generated by AirfoilInfo on '//CurDate()//' at '//CurTime()//'.'
      WRITE (unOutFile,'(/,A/)')  ' '

      if (AFI_Params%ColUAf > 0) then
         WRITE(unOutFile, '(20(A20,1x))') 'Alpha', 'Cl',  'Cd',  sCm,   'Cn',  'Cc', 'f_st', Prefix//'FullySep', sFullyAtt , 'Cl_lin','Cn_lin','Cl_smooth', 'Cn_smooth'
         WRITE(unOutFile, '(20(A20,1x))') '(deg)', '(-)', '(-)', '(-)', '(-)', '(-)','(-)', '(-)'             , '(-)'     ,  '(-)'  , '(-)'  , '(-)'    ,'(-)'

         ! TODO, we could do something with ColCpmim and ColUAf
         if (AFI_Params%ColCm > 0) then
            do iRow=1,size(table%Alpha)
               WRITE(unOutFile, '(20(F20.6,1x))') table%Alpha(iRow)*R2D, table%Coefs(iRow,AFI_Params%ColCl), table%Coefs(iRow,AFI_Params%ColCd), table%Coefs(iRow,AFI_Params%ColCm), &
                                              cn(iRow), cc(iRow),  table%Coefs(iRow,AFI_Params%ColUAf:), cl_lin(iRow), cn_lin(iRow), cl_smooth(iRow), cn_smooth(iRow)
            end do
         else
            do iRow=1,size(table%Alpha)
               WRITE(unOutFile, '(20(F20.6,1x))') table%Alpha(iRow)*R2D, table%Coefs(iRow,AFI_Params%ColCl), table%Coefs(iRow,AFI_Params%ColCd), 0.0_ReKi, &
                                              cn(iRow), cc(iRow), table%Coefs(iRow,AFI_Params%ColUAf:), cl_lin(iRow), cn_lin(iRow), cl_smooth(iRow), cn_smooth(iRow)
            end do
         endif
      else
         WRITE(unOutFile, '(20(A20,1x))') 'Alpha', 'Cl',  'Cd',  sCm,   'Cn',  'Cc',  'Cl_lin','Cn_lin','Cl_smooth', 'Cn_smooth'
         WRITE(unOutFile, '(20(A20,1x))') '(deg)', '(-)', '(-)', '(-)', '(-)', '(-)', '(-)'   ,'(-)'   ,'(-)'      ,'(-)'

         ! TODO, we could do something with ColCpmim and ColUAf
         if (AFI_Params%ColCm > 0) then
            do iRow=1,size(table%Alpha)
               WRITE(unOutFile, '(20(F20.6,1x))') table%Alpha(iRow)*R2D, table%Coefs(iRow,AFI_Params%ColCl), table%Coefs(iRow,AFI_Params%ColCd), table%Coefs(iRow,AFI_Params%ColCm), &
                                              cn(iRow), cn(iRow), cl_lin(iRow), cn_lin(iRow), cl_smooth(iRow), cn_smooth(iRow)
            end do
         else
            do iRow=1,size(table%Alpha)
               WRITE(unOutFile, '(20(F20.6,1x))') table%Alpha(iRow)*R2D, table%Coefs(iRow,AFI_Params%ColCl), table%Coefs(iRow,AFI_Params%ColCd), 0.0_ReKi, &
                                              cn(iRow), cn(iRow), cl_lin(iRow), cn_lin(iRow), cl_smooth(iRow), cn_smooth(iRow)
            end do
         endif
      
      end if
      
      CLOSE(unOutFile)
      
      if(allocated(cl_smooth)) deallocate(cl_smooth)
      if(allocated(cn_smooth)) deallocate(cn_smooth)
      if(allocated(cn       )) deallocate(cn       )
      if(allocated(cc       )) deallocate(cc       )
      if(allocated(cl_lin   )) deallocate(cl_lin   )
      if(allocated(cn_lin   )) deallocate(cn_lin   )
   enddo

end subroutine AFI_WrTables
!=============================================================================
   
END MODULE AirfoilInfo
