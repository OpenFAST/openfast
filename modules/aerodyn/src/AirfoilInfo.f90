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

   IMPLICIT NONE

   PRIVATE

   PUBLIC                                       :: AFI_Init ! routine to initialize AirfoilInfo parameters
   PUBLIC                                       :: AFI_ComputeUACoefs        ! routine to calculate Airfoil BL parameters for UA
   PUBLIC                                       :: AFI_ComputeAirfoilCoefs   ! routine to perform 1D (AOA) or 2D (AOA, Re) or (AOA, UserProp) lookup of the airfoil coefs
   TYPE(ProgDesc), PARAMETER                    :: AFI_Ver = ProgDesc( 'AirfoilInfo', '', '')    ! The name, version, and date of AirfoilInfo.

   integer, parameter                           :: MaxNumAFCoeffs = 6 !cl,cd,cm,cpMin, UA_HGM:f_st, UA_HGM:cl_fs

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

      INTEGER                                :: Co                            ! The index into the coefficients array.
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
         

      CALL ReadAFfile ( InitInput%FileName, NumCoefs, InitInput%InCol_Alfa &
                        , InitInput%InCol_Cl, InitInput%InCol_Cd, InitInput%InCol_Cm, InitInput%InCol_Cpmin, p &
                        , ErrStat2, ErrMsg2, UnEc ) 
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

      
! We need to deal with constant data.


      do iTable = 1, p%NumTabs
               ! We need to deal with constant data.
         IF ( p%Table(iTable)%ConstData )  THEN

            CALL SetErrStat ( ErrID_FATAL, 'The part to deal with constant data in AFI_Init is not written yet!', ErrStat, ErrMsg, RoutineName )
            CALL Cleanup()
            RETURN

         END IF
         
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

         IF ( ALLOCATED( Coef ) ) DEALLOCATE( Coef )

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
   SUBROUTINE ReadAFfile ( AFfile, NumCoefsIn, InCol_Alfa, InCol_Cl, InCol_Cd, InCol_Cm, InCol_Cpmin, p &
                         , ErrStat, ErrMsg, UnEc )


         ! This routine reads an airfoil file.


         ! Argument declarations.

      INTEGER(IntKi),    INTENT(IN)           :: InCol_Alfa                    ! The airfoil-table input column for angle of attack.
      INTEGER(IntKi),    INTENT(IN)           :: InCol_Cd                      ! The airfoil-table input column for drag coefficient.
      INTEGER(IntKi),    INTENT(IN)           :: InCol_Cl                      ! The airfoil-table input column for lift coefficient.
      INTEGER(IntKi),    INTENT(IN)           :: InCol_Cm                      ! The airfoil-table input column for pitching-moment coefficient.
      INTEGER(IntKi),    INTENT(IN)           :: InCol_Cpmin                   ! The airfoil-table input column for minimum pressure coefficient.
      INTEGER(IntKi),    INTENT(  OUT)        :: ErrStat                       ! Error status.
      INTEGER(IntKi),    INTENT(IN   )        :: NumCoefsIn                    ! The number of aerodynamic coefficients to be stored.

      INTEGER,           INTENT(IN)           :: UnEc                          ! I/O unit for echo file. If present and > 0, write to UnEc.      CHARACTER(*), INTENT(IN)               :: AFfile                        ! The file to be read.

      CHARACTER(*),      INTENT(IN)           :: AFfile                        ! The file to read in.
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
      INTEGER(IntKi)                          :: NumCoefsTab                   ! The number of aerodynamic coefficients to be stored for this table.

      INTEGER(IntKi)                          :: DefaultInterpOrd              ! value of default interp order
      INTEGER(IntKi)                          :: ErrStat2                      ! Error status local to this routine.
      CHARACTER(ErrMsgLen)                    :: ErrMsg2
      CHARACTER(*), PARAMETER                 :: RoutineName = 'ReadAFfile'
      CHARACTER(10)                           :: defaultStr
      CHARACTER(1024)                           :: PriPath
    
      ErrStat = ErrID_None
      ErrMsg  = ""
      defaultStr = ""

      ! Getting parent folder of airfoils data (e.g. "Arifoils/")
      CALL GetPath( AFFile, PriPath )
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
      
      ! Default to linear interpolation
   DefaultInterpOrd = 1      
      
      CALL ParseVarWDefault ( FileInfo, CurLine, 'InterpOrd', p%InterpOrd, DefaultInterpOrd, ErrStat2, ErrMsg2, UnEc )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

         
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
      
      ! Reading Boundary layer file  for aeroacoustics
      CALL ParseVar ( FileInfo, CurLine, 'BL_file' , p%BL_file , ErrStat2, ErrMsg2, UnEc )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( PathIsRelative( p%BL_file ) )  p%BL_file=trim(PriPath)//trim(p%BL_file)
         
    ! How many columns do we need to read in the input and how many total coefficients will be used?

      Cols2Parse = MAX( InCol_Alfa, InCol_Cl, InCol_Cd, InCol_Cm, InCol_Cpmin )
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

      ALLOCATE ( p%Table( p%NumTabs ) , STAT=ErrStat2 )
      IF ( ErrStat2 /= 0 )  THEN
         CALL SetErrStat ( ErrID_Fatal, 'Error allocating memory for p%Table.', ErrStat, ErrMsg, RoutineName )
         CALL Cleanup()
         RETURN
      ENDIF

      DO iTable=1,p%NumTabs
         NumCoefsTab = NumCoefsIn  ! Reset this counter for each table

         CALL ParseVar ( FileInfo, CurLine, 'Re', p%Table(iTable)%Re, ErrStat2, ErrMsg2, UnEc )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            IF (ErrStat >= AbortErrLev) THEN
               CALL Cleanup()
               RETURN
            END IF
         p%Table(iTable)%Re = p%Table(iTable)%Re * 1.0e6  ! Entered in millions, so multiply here
         IF ( p%Table(iTable)%Re <= 0.0 )  THEN
               CALL SetErrStat ( ErrID_Severe, 'Re must be > 0 in "'//TRIM( AFfile ) &
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
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            IF (ErrStat >= AbortErrLev) THEN
               CALL Cleanup()
               RETURN
            END IF

         IF ( p%Table(iTable)%InclUAdata )  THEN

               CALL ParseVar ( FileInfo, CurLine, 'alpha0', p%Table(iTable)%UA_BL%alpha0, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               p%Table(iTable)%UA_BL%alpha0 = p%Table(iTable)%UA_BL%alpha0*D2R

               CALL ParseVar ( FileInfo, CurLine, 'alpha1', p%Table(iTable)%UA_BL%alpha1, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               p%Table(iTable)%UA_BL%alpha1 = p%Table(iTable)%UA_BL%alpha1*D2R

               CALL ParseVar ( FileInfo, CurLine, 'alpha2', p%Table(iTable)%UA_BL%alpha2, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               p%Table(iTable)%UA_BL%alpha2 = p%Table(iTable)%UA_BL%alpha2*D2R

               CALL ParseVar ( FileInfo, CurLine, 'eta_e', p%Table(iTable)%UA_BL%eta_e, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVar ( FileInfo, CurLine, 'C_nalpha', p%Table(iTable)%UA_BL%C_nalpha, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
              
!>>> add after this feature gets tested better:
!               CALL ParseVar ( FileInfo, CurLine, 'C_lalpha', p%Table(iTable)%UA_BL%C_lalpha, ErrStat2, ErrMsg2, UnEc )
!                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
!<<<
!>>> remove after this feature gets tested better:
p%Table(iTable)%UA_BL%C_lalpha = p%Table(iTable)%UA_BL%C_nalpha
!<<<
            
               CALL ParseVarWDefault ( FileInfo, CurLine, 'T_f0', p%Table(iTable)%UA_BL%T_f0, 3.0_ReKi, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         
               CALL ParseVarWDefault ( FileInfo, CurLine, 'T_V0', p%Table(iTable)%UA_BL%T_V0, 6.0_ReKi, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               
               CALL ParseVarWDefault ( FileInfo, CurLine, 'T_p', p%Table(iTable)%UA_BL%T_p, 1.7_ReKi, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVarWDefault ( FileInfo, CurLine, 'T_VL', p%Table(iTable)%UA_BL%T_VL, 11.0_ReKi, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVarWDefault ( FileInfo, CurLine, 'b1', p%Table(iTable)%UA_BL%b1, .14_ReKi, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVarWDefault ( FileInfo, CurLine, 'b2', p%Table(iTable)%UA_BL%b2, .53_ReKi, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVarWDefault ( FileInfo, CurLine, 'b5', p%Table(iTable)%UA_BL%b5, 5.0_ReKi, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVarWDefault ( FileInfo, CurLine, 'A1', p%Table(iTable)%UA_BL%A1, .3_ReKi, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVarWDefault ( FileInfo, CurLine, 'A2', p%Table(iTable)%UA_BL%A2, .7_ReKi, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVarWDefault ( FileInfo, CurLine, 'A5', p%Table(iTable)%UA_BL%A5, 1.0_ReKi, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVar ( FileInfo, CurLine, 'S1', p%Table(iTable)%UA_BL%S1, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVar ( FileInfo, CurLine, 'S2', p%Table(iTable)%UA_BL%S2, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVar ( FileInfo, CurLine, 'S3', p%Table(iTable)%UA_BL%S3, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVar ( FileInfo, CurLine, 'S4', p%Table(iTable)%UA_BL%S4, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVar ( FileInfo, CurLine, 'Cn1', p%Table(iTable)%UA_BL%Cn1, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVar ( FileInfo, CurLine, 'Cn2', p%Table(iTable)%UA_BL%Cn2, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVarWDefault ( FileInfo, CurLine, 'St_sh', p%Table(iTable)%UA_BL%St_sh, .19_ReKi, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVar ( FileInfo, CurLine, 'Cd0', p%Table(iTable)%UA_BL%Cd0, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVar ( FileInfo, CurLine, 'Cm0', p%Table(iTable)%UA_BL%Cm0, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVar ( FileInfo, CurLine, 'k0', p%Table(iTable)%UA_BL%k0, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVar ( FileInfo, CurLine, 'k1', p%Table(iTable)%UA_BL%k1, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVar ( FileInfo, CurLine, 'k2', p%Table(iTable)%UA_BL%k2, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVar ( FileInfo, CurLine, 'k3', p%Table(iTable)%UA_BL%k3, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVar ( FileInfo, CurLine, 'k1_hat', p%Table(iTable)%UA_BL%k1_hat, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVarWDefault ( FileInfo, CurLine, 'x_cp_bar', p%Table(iTable)%UA_BL%x_cp_bar, .2_ReKi, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                  
               CALL ParseVarWDefault ( FileInfo, CurLine, 'UACutout', p%Table(iTable)%UA_BL%UACutout, 45.0_ReKi, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               p%Table(iTable)%UA_BL%UACutout = p%Table(iTable)%UA_BL%UACutout*D2R

               CALL ParseVarWDefault ( FileInfo, CurLine, 'filtCutOff', p%Table(iTable)%UA_BL%filtCutOff, 0.5_ReKi, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                  
               p%ColUAf    = NumCoefsIn + 1 ! column for f_st
               NumCoefsTab = p%ColUAf    + 1 ! precompute f_st and cl_fs for the HGM model
               
               IF (ErrStat >= AbortErrLev) THEN
                  CALL Cleanup()
                  RETURN
               END IF
                  
         ELSE !ELSEIF ( UA_Model == 1 )  THEN

            !CALL SetErrStat ( ErrID_Fatal &
            !, 'You must supply Beddoes-Leishman unsteady aerodynamics parameters for all airfoils if you want to use that' &
            !//' model. You did not do so for Table #'//TRIM( Num2LStr( iTable ) )//' in the "'//TRIM( AFfile )//'" airfoil file.', ErrStat, ErrMsg, RoutineName )
            !CALL Cleanup()
            !RETURN

         ENDIF ! ( p%Table(iTable)%InclUAdata )

         CALL ParseVar ( FileInfo, CurLine, 'NumAlf', p%Table(iTable)%NumAlf, ErrStat2, ErrMsg2, UnEc )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            IF (ErrStat >= AbortErrLev) THEN
               CALL Cleanup()
               RETURN
            END IF

         IF ( p%Table(iTable)%NumAlf < 1 )  THEN
            CALL SetErrStat( ErrID_Fatal, 'NumAlf must be a positive number on line #' &
                           //TRIM( Num2LStr( FileInfo%FileLine(CurLine-1) ) )//' in "'//TRIM( AFfile )//'".', ErrStat, ErrMsg, RoutineName )
            CALL Cleanup()
            RETURN
         ELSEIF ( p%Table(iTable)%NumAlf < 3 )  THEN
            p%Table(iTable)%ConstData = .TRUE.
         ELSE
            p%Table(iTable)%ConstData = .FALSE.
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

         DO Row=1,p%Table(iTable)%NumAlf

            CALL ParseAry ( FileInfo, CurLine, 'CoeffData', SiAry, Cols2Parse, ErrStat2, ErrMsg2, UnEc )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               IF (ErrStat >= AbortErrLev) THEN
                  CALL Cleanup()
                  RETURN
               END IF

            p%Table(iTable)%Alpha(Row        ) = SiAry(InCol_Alfa)*D2R
            p%Table(iTable)%Coefs(Row,p%ColCl) = SiAry(InCol_Cl  )
            p%Table(iTable)%Coefs(Row,p%ColCd) = SiAry(InCol_Cd  )

            IF ( InCol_Cm    > 0 ) p%Table(iTable)%Coefs(Row,p%ColCm   ) = SiAry(InCol_Cm)
            IF ( InCol_Cpmin > 0 ) p%Table(iTable)%Coefs(Row,p%ColCpmin) = SiAry(InCol_Cpmin)

         ENDDO ! Row

         call CalculateUACoeffs(p%Table(iTable), p%ColCl, p%ColUAf)

            ! Let's make sure that the data go from -Pi to Pi and that the values are the same for both
            ! unless there is only one point.

         IF ( .NOT. p%Table(iTable)%ConstData )  THEN
            NumAlf  = p%Table(iTable)%NumAlf
            BadVals = .FALSE.
            IF ( .NOT. EqualRealNos( p%Table(iTable)%Alpha(1), -Pi ) )  THEN
               BadVals = .TRUE.
            ENDIF
            IF ( .NOT. EqualRealNos( p%Table(iTable)%Alpha(NumAlf), Pi ) )  THEN
               BadVals = .TRUE.
            ENDIF
            DO Coef=1,NumCoefsTab
               IF ( .NOT. EqualRealNos( p%Table(iTable)%Coefs(1,Coef), p%Table(iTable)%Coefs(NumAlf,Coef) ) )  THEN
                  BadVals = .TRUE.
               ENDIF
            ENDDO ! Coef
            IF ( BadVals )  THEN
!               CALL SetErrStat( ErrID_Fatal &
               CALL SetErrStat( ErrID_Warn, &
                  'Airfoil data should go from -180 degrees to 180 degrees and the coefficients at the ends should be the same.', ErrStat, ErrMsg, RoutineName )
               !CALL Cleanup()
               !RETURN
            ENDIF
         ENDIF ! ( .NOT. p%Table(iTable)%ConstData )

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
            IF ( ALLOCATED(SiAry) ) DEALLOCATE(SiAry)
            
         END SUBROUTINE Cleanup 


   END SUBROUTINE ReadAFfile
!----------------------------------------------------------------------------------------------------------------------------------  
   SUBROUTINE CalculateUACoeffs(p,ColCl,ColUAf)
      TYPE (AFI_Table_Type),    intent(inout) :: p                             ! This structure stores all the module parameters that are set by AirfoilInfo during the initialization phase.
      integer(IntKi),           intent(in   ) :: ColCl                         ! column for cl
      integer(IntKi),           intent(in   ) :: ColUAf                        ! column for UA f_st
   
      INTEGER(IntKi)                          :: Row                           ! The row of a table to be parsed in the FileInfo structure.
      INTEGER(IntKi)                          :: col_clFs                      ! column for UA cl_fs
      
      REAL(ReKi)                              :: cl_ratio, cl_inv
      REAL(ReKi)                              :: f_st, cl_fs
      REAL(ReKi)                              :: f_iHigh, f_iLow
      INTEGER(IntKi)                          :: iHigh, iLow

      col_clFs = ColUAf + 1

      if ( p%InclUAdata )  then
         p%UA_BL%UACutout_blend = max(0.0_ReKi, p%UA_BL%UACutout - 5.0_ReKi*D2R) ! begin turning off 5 degrees before (or at 0 degrees)
      
         if (EqualRealNos(p%UA_BL%c_lalpha,0.0_ReKi)) then
            p%Coefs(:,ColUAf)   = 0.0_ReKi                           ! Eq. 59
            p%Coefs(:,col_clFs) = p%Coefs(:,ColCl)                   ! Eq. 61
         else
         
            f_iHigh = huge(f_iHigh)
            f_iLow  = f_iHigh
            iHigh = 0
            iLow = 0
         
            do Row=1,p%NumAlf
            
               if (EqualRealNos( p%alpha(Row), p%UA_BL%alpha0)) then
                  f_st  = 1.0_ReKi                                         ! Eq. 59
                  p%Coefs(Row,col_clFs) = p%Coefs(Row,ColCl) / 2.0_ReKi    ! Eq. 61 (which should be very close to 0 because definition of alpha0 says cl(alpha0) = 0 )
               else
            
                  cl_ratio = p%Coefs(Row,ColCl) / ( p%UA_BL%c_lalpha*(p%alpha(Row) - p%UA_BL%alpha0))
                  cl_ratio = max(0.0_ReKi, cl_ratio)

                  f_st = ( 2.0_ReKi * sqrt(cl_ratio) - 1.0_ReKi )**2
                  
                  if (f_st < 1.0_ReKi) then 
                     ! Region where f_st<1, merge
                     f_st  = max(0.0_ReKi, f_st) ! make sure it is not negative
                     cl_fs = (p%Coefs(Row,ColCl) - p%UA_BL%c_lalpha* (p%alpha(Row) - p%UA_BL%alpha0)*f_st) / (1.0_ReKi - f_st) ! Eq 61
                  else
                     ! Initialize to linear region (in fact only at singularity, where f_st=1)
                     f_st = 1.0_ReKi
                     cl_fs = p%Coefs(Row,ColCl) / 2.0_ReKi                      ! Eq. 61
                  end if
                  
                  if (p%alpha(Row) < p%UA_BL%alpha0) then
                     if (f_st <= f_iLow) then
                        f_iLow = f_st
                        iLow = Row
                     end if
                  else !p%alpha(Row) > p%UA_BL%alpha0 (note that they can't be equal)
                     if (f_st < f_iHigh) then
                        f_iHigh = f_st
                        iHigh = Row
                     end if
                  end if
               end if
               
               p%Coefs(Row,ColUAf)   = f_st
               p%Coefs(Row,col_clFs) = cl_fs
            end do
            if (iLow >0) p%Coefs(1:iLow,col_clFs) = p%Coefs(1:iLow,ColCl)
            if (iHigh>0) p%Coefs(iHigh:,col_clFs) = p%Coefs(iHigh:,ColCl)
            
            
            ! Ensuring everything is in harmony 
            do Row=1,p%NumAlf
               cl_fs = p%Coefs(Row,col_clFs)
               
               cl_inv = p%UA_BL%c_lalpha*(p%alpha(Row) - p%UA_BL%alpha0)     ! Eq. 64
               if (.not. EqualRealNos(cl_inv, cl_fs)) then
                  f_st=(p%Coefs(Row,ColCl) - cl_fs) / (cl_inv-cl_fs);        ! Eq. 60
                  f_st = max(0.0_ReKi, f_st)
                  f_st = min(1.0_ReKi, f_st)
                  
                  p%Coefs(Row,ColUAf) = f_st
               else 
                  p%Coefs(Row,ColUAf) = 1.0_ReKi
               end if

            end do
         end if ! c_lalpha == 0
         
      end if

            

      
   END SUBROUTINE CalculateUACoeffs
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
   
   Alpha = AOA
   call MPi2Pi ( Alpha ) ! change AOA into range of -pi to pi
   
   
      ! Spline interpolation of lower table based on requested AOA
   
   IntAFCoefs(1:s1) = CubicSplineInterpM( Alpha  &
                                          , p%Table(iTab)%Alpha &
                                          , p%Table(iTab)%Coefs &
                                          , p%Table(iTab)%SplineCoefs &
                                          , ErrStat, ErrMsg )
   
  
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
      AFI_interp%f_st = IntAFCoefs(p%ColUAf)
      AFI_interp%cl_fs = IntAFCoefs(p%ColUAf+1)
   else
      AFI_interp%f_st = 0.0_ReKi
      AFI_interp%cl_fs = 0.0_ReKi
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
END MODULE AirfoilInfo
