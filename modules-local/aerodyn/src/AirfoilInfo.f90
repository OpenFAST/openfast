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
   USE                                             NWTC_FitPack
   USE                                          :: ISO_FORTRAN_ENV , ONLY : IOSTAT_EOR
   
   IMPLICIT NONE

   PRIVATE

   PUBLIC                                       :: AFI_Init ! routine to initialize AirfoilInfo parameters
   PUBLIC                                       :: AFI_GetAirfoilParams ! routine to calculate Airfoil parameters
   PUBLIC                                       :: AFI_ComputeUACoefs2D      ! perform multi-table interpolation of UA parameters
   PUBLIC                                       :: AFI_ComputeAirfoilCoefs2D ! routine to perform a 2D (AOA, Re) lookup of the airfoil coefs
   PUBLIC                                       :: AFI_ComputeAirfoilCoefs1D ! routine to perform a 1D (AOA) lookup of the airfoil coefs
   
   TYPE(ProgDesc), PARAMETER                    :: AFI_Ver = ProgDesc( 'AirfoilInfo', '', '')    ! The name, version, and date of AirfoilInfo.


   CONTAINS


   function CheckValuesAreUniqueMonotonicIncreasing(secondVals)
     
      real(ReKi),  intent(in   )  :: secondVals(:)
      logical CheckValuesAreUniqueMonotonicIncreasing
      
      
      integer(IntKi) :: i
      
      CheckValuesAreUniqueMonotonicIncreasing = .true.
      
      do i = 2, size(secondVals)
         if ( EqualRealNos(secondVals(i)-secondVals(i-1),0.0_ReKi) .or. (secondVals(i) < secondVals(i-1))) then
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
      INTEGER                                :: CoefSiz                       ! The size of the Coef array.
      INTEGER                                :: File                          ! The file index.
      INTEGER                                :: IA                            ! Index into the alpha array.
      INTEGER                                :: Indx                          ! Index into the arrays.
      INTEGER                                :: NumAl                         ! The length of the alpha array being sent to the regrid interface routine.
      INTEGER                                :: NumRe                         ! The length of the Re array being sent to the regrid interface routine.
      INTEGER                                :: Table                         ! Index into the Table array.
      INTEGER                                :: UnEc                          ! Local echo file unit number
      INTEGER                                :: NumCoefs                      ! The number of aerodynamic coefficients to be stored
      integer                                :: iTable                        ! Iterator for airfoil tables
      
      INTEGER                                :: ErrStat2                      ! Local error status.
      CHARACTER(ErrMsgLen)                   :: ErrMsg2
      CHARACTER(*), PARAMETER                :: RoutineName = 'AFI_Init'
      REAL(ReKi), ALLOCATABLE                :: secondVals(:)                 ! The values of the 2nd dependent variable when using multiple airfoil tables
      ErrStat = ErrID_None
      ErrMsg  = ""
 
         ! Display the version for this module.

      CALL DispNVD ( AFI_Ver )

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

            IF ( p%AFTabMod == 2 .or. p%AFTabMod == 3 ) THEN
               
               ALLOCATE(secondVals(p%NumTabs), STAT=ErrStat2 )
                  IF ( ErrStat2 /= 0 )  THEN
                     CALL SetErrStat ( ErrID_Fatal, 'Error allocating memory for the secondVals array.', ErrStat, ErrMsg, RoutineName )
                     RETURN
                  ENDIF
               IF (p%AFTabMod == 2) THEN
                     
              
                     secondVals(1) = p%Table(1)%Re
                     
                  ELSE IF (p%AFTabMod == 3) THEN
                     
               
                     secondVals(1) = p%Table(1)%UserProp
                     
                  END IF   
               ! Do the tables have the same number and set of alphas.

               DO Table=2,p%NumTabs

                  IF ( p%Table(Table)%NumAlf /= p%Table(1)%NumAlf )  THEN
                     CALL SetErrStat ( ErrID_Fatal, 'Fatal Error: Airfoil file "'//TRIM( InitInput%FileName ) &
                                       //'", Table #'//TRIM( Num2LStr( Table ) ) &
                                       //' does not have the same number of alphas as the first table.', ErrStat, ErrMsg, RoutineName )
                     CALL Cleanup()
                     RETURN
                  ENDIF

                  DO IA=1,p%Table(1)%NumAlf
                     IF ( p%Table(Table)%Alpha(IA) /= p%Table(1)%Alpha(IA) )  THEN
                        CALL SetErrStat ( ErrID_Fatal, 'Fatal Error: Airfoil file "'//TRIM( InitInput%FileName ) &
                                       //'", Table #'//TRIM( Num2LStr( Table ) ) &
                                       //' does not have the same set of alphas as the first table.', ErrStat, ErrMsg, RoutineName )
                        CALL Cleanup()
                        RETURN
                     ENDIF
                  END DO ! IA
                  
                  IF (p%AFTabMod == 2) THEN
                     
                        ! Ctrl Value must be the same, Re must be monotonic increasing without repeats
                     IF ( p%Table(Table)%UserProp /= p%Table(1)%UserProp )  THEN
                        CALL SetErrStat ( ErrID_Fatal, 'Fatal Error: airfoil file "'//TRIM( InitInput%FileName ) &
                                       //'", Table #'//TRIM( Num2LStr( Table ) ) &
                                       //' does not have the same value for Ctrl Property (CtrlProp) as the first table.', ErrStat, ErrMsg, RoutineName )
                        CALL Cleanup()
                        RETURN
                     ENDIF
                     secondVals(Table) = p%Table(Table)%Re
                     
                  ELSE IF (p%AFTabMod == 3) THEN
                     
                        ! Re must be the same, Ctrl must be different
                     IF ( p%Table(Table)%Re /= p%Table(1)%Re )  THEN
                        CALL SetErrStat ( ErrID_Fatal, 'Fatal Error: airfoil file "'//TRIM( InitInput%FileName ) &
                                       //'", Table #'//TRIM( Num2LStr( Table ) ) &
                                       //' does not have the same value for Re Property (ReProp) as the first table.', ErrStat, ErrMsg, RoutineName )
                        CALL Cleanup()
                        RETURN
                     ENDIF
                     secondVals(Table) = p%Table(Table)%UserProp
                     
                  END IF
                  
               ENDDO ! Tab
               
                  ! Make sure all Re's are monotonic and increasing with no repeats
               IF (p%AFTabMod > 1) THEN
                  
                  IF (.NOT. CheckValuesAreUniqueMonotonicIncreasing(secondVals)) THEN
                     ErrMsg2 = 'Fatal Error: airfoil file "'//TRIM( InitInput%FileName ) &
                                 //'", is not monotonic and increasing in the '
                     IF (p%AFTabMod == 2) THEN
                        ErrMsg2 = trim(ErrMsg2)//' Re Property (ReProp).'
                     ELSE
                        ErrMsg2 = trim(ErrMsg2)//' Ctrl Property (CtrlProp).'
                     END IF
                     
                     CALL SetErrStat ( ErrID_Fatal, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                     CALL Cleanup()
                     RETURN
                  ENDIF
                  deallocate(secondVals)
               END IF
                
            ELSE
               p%NumTabs = 1
               CALL SetErrStat ( ErrID_Warn, 'DimModel = 1D, therefore only using the first airfoil table in the file: "'//TRIM( InitInput%FileName ), ErrStat, ErrMsg, RoutineName )
            END IF
            
         ENDIF ! ( p%NumTabs > 1 )

! We need to deal with constant data.


         do iTable = 1, p%NumTabs
               ! Allocate the arrays to hold spline coefficients.

            allocate ( p%Table(iTable)%SplineCoefs( p%Table(iTable)%NumAlf-1 &
                     , NumCoefs, 0:3 ), STAT=ErrStat2 )
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


            ! Compute the spline coefficients of the piecewise cubic polynomials for the irregularly-spaced airfoil data in each file.
            ! Unless the data are constant.

         DO Co=1,NumCoefs
            
                  ! We use 1D cubic spline interpolation if the data are not constant.
               IF ( p%Table(1)%ConstData )  THEN

                     ! We need to deal with constant data.

                  CALL SetErrStat ( ErrID_FATAL, 'The part to deal with constant data in AFI_Init is not written yet!', ErrStat, ErrMsg, RoutineName )
                  CALL Cleanup()
                  RETURN

               ENDIF ! p%Table(1)%ConstData

         END DO ! Co


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
   SUBROUTINE ReadAFfile ( AFfile, NumCoefs, InCol_Alfa, InCol_Cl, InCol_Cd, InCol_Cm, InCol_Cpmin, p &
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

      TYPE (AFI_ParameterType), INTENT(INOUT)   :: p                          ! This structure stores all the module parameters that are set by AirfoilInfo during the initialization phase.


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

      DO Table=1,p%NumTabs

         CALL ParseVar ( FileInfo, CurLine, 'Re', p%Table(Table)%Re, ErrStat2, ErrMsg2, UnEc )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            IF (ErrStat >= AbortErrLev) THEN
               CALL Cleanup()
               RETURN
            END IF
         IF ( p%Table(Table)%Re <= 0.0 )  THEN
               CALL SetErrStat ( ErrID_Severe, 'Re must be > 0 in "'//TRIM( AFfile ) &
                                      //'".'//NewLine//'  >> The error occurred on line #' &
                                      //TRIM( Num2LStr( FileInfo%FileLine(CurLine-1) ) )//'.', ErrStat, ErrMsg, RoutineName )
               CALL Cleanup()
               RETURN
         ENDIF ! ( p%Table(Table)%Re <= 0.0 )

         CALL ParseVar ( FileInfo, CurLine, 'Ctrl', p%Table(Table)%UserProp, ErrStat2, ErrMsg2, UnEc )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            IF (ErrStat >= AbortErrLev) THEN
               CALL Cleanup()
               RETURN
            END IF
         !IF ( p%Table(Table)%UserProp /= 0.0 )  THEN !bjj: use EqualRealNos?
         !      CALL SetErrStat ( ErrID_Severe, 'UserProp must equal 0.0 in "'//TRIM( AFfile ) &
         !                              //'".'//NewLine//'  >> The error occurred on line #' &
         !                              //TRIM( Num2LStr( FileInfo%FileLine(CurLine-1) ) )//'.', ErrStat, ErrMsg, RoutineName )
         !      CALL Cleanup()
         !      RETURN
         !ENDIF ! ( p%Table(Table)%UserProp <= 0.0 )

         CALL ParseVar ( FileInfo, CurLine, 'InclUAdata', p%Table(Table)%InclUAdata, ErrStat2, ErrMsg2, UnEc )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            IF (ErrStat >= AbortErrLev) THEN
               CALL Cleanup()
               RETURN
            END IF

         IF ( p%Table(Table)%InclUAdata )  THEN

               CALL ParseVar ( FileInfo, CurLine, 'alpha0', p%Table(Table)%UA_BL%alpha0, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVar ( FileInfo, CurLine, 'alpha1', p%Table(Table)%UA_BL%alpha1, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVar ( FileInfo, CurLine, 'alpha2', p%Table(Table)%UA_BL%alpha2, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVar ( FileInfo, CurLine, 'eta_e', p%Table(Table)%UA_BL%eta_e, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVar ( FileInfo, CurLine, 'C_nalpha', p%Table(Table)%UA_BL%C_nalpha, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            
               CALL ParseVarWDefault ( FileInfo, CurLine, 'T_f0', p%Table(Table)%UA_BL%T_f0, 3.0_ReKi, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         
               CALL ParseVarWDefault ( FileInfo, CurLine, 'T_V0', p%Table(Table)%UA_BL%T_V0, 6.0_ReKi, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               
               CALL ParseVarWDefault ( FileInfo, CurLine, 'T_p', p%Table(Table)%UA_BL%T_p, 1.7_ReKi, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVarWDefault ( FileInfo, CurLine, 'T_VL', p%Table(Table)%UA_BL%T_VL, 11.0_ReKi, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVarWDefault ( FileInfo, CurLine, 'b1', p%Table(Table)%UA_BL%b1, .14_ReKi, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVarWDefault ( FileInfo, CurLine, 'b2', p%Table(Table)%UA_BL%b2, .53_ReKi, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVarWDefault ( FileInfo, CurLine, 'b5', p%Table(Table)%UA_BL%b5, 5.0_ReKi, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVarWDefault ( FileInfo, CurLine, 'A1', p%Table(Table)%UA_BL%A1, .3_ReKi, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVarWDefault ( FileInfo, CurLine, 'A2', p%Table(Table)%UA_BL%A2, .7_ReKi, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVarWDefault ( FileInfo, CurLine, 'A5', p%Table(Table)%UA_BL%A5, 1.0_ReKi, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVar ( FileInfo, CurLine, 'S1', p%Table(Table)%UA_BL%S1, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVar ( FileInfo, CurLine, 'S2', p%Table(Table)%UA_BL%S2, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVar ( FileInfo, CurLine, 'S3', p%Table(Table)%UA_BL%S3, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVar ( FileInfo, CurLine, 'S4', p%Table(Table)%UA_BL%S4, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVar ( FileInfo, CurLine, 'Cn1', p%Table(Table)%UA_BL%Cn1, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVar ( FileInfo, CurLine, 'Cn2', p%Table(Table)%UA_BL%Cn2, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVarWDefault ( FileInfo, CurLine, 'St_sh', p%Table(Table)%UA_BL%St_sh, .19_ReKi, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVar ( FileInfo, CurLine, 'Cd0', p%Table(Table)%UA_BL%Cd0, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVar ( FileInfo, CurLine, 'Cm0', p%Table(Table)%UA_BL%Cm0, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVar ( FileInfo, CurLine, 'k0', p%Table(Table)%UA_BL%k0, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVar ( FileInfo, CurLine, 'k1', p%Table(Table)%UA_BL%k1, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVar ( FileInfo, CurLine, 'k2', p%Table(Table)%UA_BL%k2, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVar ( FileInfo, CurLine, 'k3', p%Table(Table)%UA_BL%k3, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVar ( FileInfo, CurLine, 'k1_hat', p%Table(Table)%UA_BL%k1_hat, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               CALL ParseVarWDefault ( FileInfo, CurLine, 'x_cp_bar', p%Table(Table)%UA_BL%x_cp_bar, .2_ReKi, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                  
               CALL ParseVarWDefault ( FileInfo, CurLine, 'UACutout', p%Table(Table)%UA_BL%UACutout, 45.0_ReKi, ErrStat2, ErrMsg2, UnEc )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               p%Table(Table)%UA_BL%UACutout = p%Table(Table)%UA_BL%UACutout*D2R

               CALL ParseVarWDefault ( FileInfo, CurLine, 'filtCutOff', p%Table(Table)%UA_BL%filtCutOff, 20.0_ReKi, ErrStat2, ErrMsg2, UnEc )
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

         ENDIF ! ( p%Table(Table)%InclUAdata )

         CALL ParseVar ( FileInfo, CurLine, 'NumAlf', p%Table(Table)%NumAlf, ErrStat2, ErrMsg2, UnEc )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            IF (ErrStat >= AbortErrLev) THEN
               CALL Cleanup()
               RETURN
            END IF

         IF ( p%Table(Table)%NumAlf < 1 )  THEN
            CALL SetErrStat( ErrID_Fatal, 'NumAlf must be a positive number on line #' &
                           //TRIM( Num2LStr( FileInfo%FileLine(CurLine-1) ) )//' in "'//TRIM( AFfile )//'".', ErrStat, ErrMsg, RoutineName )
            CALL Cleanup()
            RETURN
         ELSEIF ( p%Table(Table)%NumAlf < 3 )  THEN
            p%Table(Table)%ConstData = .TRUE.
         ELSE
            p%Table(Table)%ConstData = .FALSE.
         ENDIF ! ( Test for valid values for NumAlf )


            ! Allocate the arrays for the airfoil coefficients.

         ALLOCATE ( p%Table(Table)%Alpha( p%Table(Table)%NumAlf ), STAT=ErrStat2 )
         IF ( ErrStat2 /= 0 )  THEN
            CALL SetErrStat( ErrID_Fatal, 'Error allocating memory for p%Table%Alpha.', ErrStat, ErrMsg, RoutineName )
            CALL Cleanup()
            RETURN
         ENDIF

         ALLOCATE ( p%Table(Table)%Coefs( p%Table(Table)%NumAlf, NumCoefs ), STAT=ErrStat2 )
         IF ( ErrStat2 /= 0 )  THEN
            CALL SetErrStat( ErrID_Fatal, 'Error allocating memory for p%Table%Coefs.', ErrStat, ErrMsg, RoutineName )
            CALL Cleanup()
            RETURN
         ENDIF

         DO Row=1,p%Table(Table)%NumAlf

            CALL ParseAry ( FileInfo, CurLine, 'CoeffData', SiAry, Cols2Parse, ErrStat2, ErrMsg2, UnEc )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               IF (ErrStat >= AbortErrLev) THEN
                  CALL Cleanup()
                  RETURN
               END IF

            p%Table(Table)%Alpha(Row  ) = SiAry(InCol_Alfa)*D2R
           !p%Table(Table)%Alpha(Row  ) = SiAry(InCol_Alfa)
            p%Table(Table)%Coefs(Row,1) = SiAry(InCol_Cl  )
            p%Table(Table)%Coefs(Row,2) = SiAry(InCol_Cd  )

            IF ( InCol_Cm > 0 )  THEN
               p%Table(Table)%Coefs(Row,3) = SiAry(InCol_Cm)
               IF ( InCol_Cpmin > 0 )  p%Table(Table)%Coefs(Row,4) = SiAry(InCol_Cpmin)
            ELSE
               IF ( InCol_Cpmin > 0 )  p%Table(Table)%Coefs(Row,3) = SiAry(InCol_Cpmin)
            ENDIF ! IF ( Col_Cm > 0 )  THEN

         ENDDO ! Row


            ! Let's make sure that the data go from -Pi to Pi and that the values are the same for both
            ! unless there is only one point.

         IF ( .NOT. p%Table(Table)%ConstData )  THEN
            NumAlf  = p%Table(Table)%NumAlf
            BadVals = .FALSE.
            IF ( .NOT. EqualRealNos( p%Table(Table)%Alpha(1), -Pi ) )  THEN
            !IF ( .NOT. EqualRealNos( p%Table(Table)%Alpha(1), -180.0_ReKi ) )  THEN
               BadVals = .TRUE.
            ENDIF
            IF ( .NOT. EqualRealNos( p%Table(Table)%Alpha(NumAlf), Pi ) )  THEN
            !IF ( .NOT. EqualRealNos( p%Table(Table)%Alpha(NumAlf), 180.0_ReKi ) )  THEN
               BadVals = .TRUE.
            ENDIF
            DO Coef=1,NumCoefs
               IF ( .NOT. EqualRealNos( p%Table(Table)%Coefs(1,Coef), p%Table(Table)%Coefs(NumAlf,Coef) ) )  THEN
                  BadVals = .TRUE.
               ENDIF
            ENDDO ! Coef
            IF ( BadVals )  THEN
               CALL SetErrStat( ErrID_Warn, &
                  'Airfoil data should go from -180 degrees to 180 degrees and the coefficients at the ends should be the same.', ErrStat, ErrMsg, RoutineName )
               CALL Cleanup()
               RETURN
            ENDIF
         ENDIF ! ( .NOT. p%Table(Table)%ConstData )

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
!----------------------------------------------------------------------------------------------------------------------------------  
subroutine AFI_ComputeUACoefs2D( interpolant, secondaryDepVal, p, &
                      UA_BL, errStat, errMsg )
! This routine is calculates the UA parameters for a set of tables which are dependent a 2nd user-defined varible, could be Re or Cntrl, etc.
! If the requested yVar is not associated with a given table, then the two tables which contain yVar are found and, a cubic spline interpolation is performed at the requested AOA.
! for each of those two tables. Then a linear intepolation is performed on the 2nd dimension to find the final Cl,Cd,Cm, and Cpmin values.
! If the requested yVar corresponds to a table, then only a single cubic interpolation based on the requested AOA is performed.
!..................................................................................................................................
   integer(IntKi),           intent(in   ) :: interpolant
   real(ReKi),               intent(in   ) :: secondaryDepVal            ! Interpolate based on this value
   TYPE (AFI_ParameterType), intent(in   ) :: p                          ! This structure stores all the module parameters that are set by AirfoilInfo during the initialization phase.
   type(AFI_UA_BL_Type),     intent(  out) :: UA_BL
   integer(IntKi),           intent(  out) :: errStat                    ! Error status of the operation
   character(*),             intent(  out) :: errMsg                     ! Error message if ErrStat /= ErrID_None 
   
   
   type(AFI_UA_BL_Type) :: UA_BL1, UA_BL2
   integer                         :: s1, lowerTable, upperTable, i
   real(ReKi)                      :: s, lowerVal, dVal
      
   ErrStat = ErrID_None
   ErrMsg  = ''
   
         ! check that we have tables which contain specified secondaryDepVal.  If not, throw error
   lowerTable = 0
   upperTable = 0
   do i = 1, p%NumTabs
      if (interpolant == 1) then
         if (p%Table(i)%Re <= secondaryDepVal) lowerTable = i
         if (p%Table(i)%Re >= secondaryDepVal) upperTable = i
      else
         if (p%Table(i)%UserProp <= secondaryDepVal) lowerTable = i
         if (p%Table(i)%UserProp >= secondaryDepVal) upperTable = i
      end if
      
   end do
   
   if ( lowerTable == 0 .or. upperTable == 0 ) then
      if (interpolant == 1) then
         ErrMsg  = "Specified Reynold's number, "//trim(num2lstr(secondaryDepVal))//" , is outside the range of Re specified in the airfoil input file tables."
      else
         ErrMsg  = "Specified User Property's value, "//trim(num2lstr(secondaryDepVal))//" , is outside the range of User Property values specified in the airfoil input file tables."
      end if
      
      call SetErrStat ( ErrID_Fatal, ErrMsg, ErrStat, ErrMsg, 'AFI_ComputeAirfoilCoefs2DRe' )
      return
   end if

   UA_BL1 = p%Table(lowerTable)%UA_BL
   UA_BL2 = p%Table(upperTable)%UA_BL
   
      ! Linearly interpolate between tables
   if (interpolant == 1) then
      dVal = (p%Table(upperTable)%Re - p%Table(lowerTable)%Re)
      lowerVal = p%Table(lowerTable)%Re
   else
      dVal = (p%Table(upperTable)%UserProp - p%Table(lowerTable)%UserProp)
      lowerVal = p%Table(lowerTable)%UserProp
   end if
   
   if (EqualRealNos(dVal, 0.0_ReKi)) then
      s = 0.0_ReKi
   else                        
      s = (secondaryDepVal - lowerVal) / dVal
   end if
   
   UA_BL%alpha0 = UA_BL1%alpha0 + s*UA_BL2%alpha0
   
   
end subroutine AFI_ComputeUACoefs2D  
   
!----------------------------------------------------------------------------------------------------------------------------------  
subroutine AFI_ComputeAirfoilCoefs2D( interpolant, AOA, secondaryDepVal, p, &
                      Cl, Cd, Cm, Cpmin, errStat, errMsg )
! This routine is calculates Cl, Cd, Cm, (and Cpmin) for a set of tables which are dependent on AOA as well as a 2nd user-defined varible, could be Re or Cntrl, etc.
! If the requested yVar is not associated with a given table, then the two tables which contain yVar are found and, a cubic spline interpolation is performed at the requested AOA.
! for each of those two tables. Then a linear intepolation is performed on the 2nd dimension to find the final Cl,Cd,Cm, and Cpmin values.
! If the requested yVar corresponds to a table, then only a single cubic interpolation based on the requested AOA is performed.
!..................................................................................................................................
   integer(IntKi),           intent(in   ) :: interpolant
   real(ReKi),               intent(in   ) :: AOA
   real(ReKi),               intent(in   ) :: secondaryDepVal           ! Unused in the current version!     
   TYPE (AFI_ParameterType), intent(in   ) :: p                          ! This structure stores all the module parameters that are set by AirfoilInfo during the initialization phase.
   real(ReKi),               intent(  out) :: Cl, Cd, Cm, Cpmin
   integer(IntKi),           intent(  out) :: errStat                    ! Error status of the operation
   character(*),             intent(  out) :: errMsg                     ! Error message if ErrStat /= ErrID_None 
   
   
   real                            :: IntAFCoefs(4)                ! The interpolated airfoil coefficients.
   real(reki)                      :: Alpha
   integer                         :: s1, lowerTable, upperTable, i
   real(ReKi)                      :: s, lowerVal, dVal, Cl1, Cl2, Cd1, Cd2, Cm1, Cm2, Cpmin1, Cpmin2
      
   ErrStat = ErrID_None
   ErrMsg  = ''
   
   
      ! check that we have tables which contain specified secondaryDepVal.  If not, throw error
   lowerTable = 0
   upperTable = 0
   do i = 1, p%NumTabs
      if (interpolant == 1) then
         if (p%Table(i)%Re <= secondaryDepVal) lowerTable = i
         if (p%Table(i)%Re >= secondaryDepVal) upperTable = i
      else
         if (p%Table(i)%UserProp <= secondaryDepVal) lowerTable = i
         if (p%Table(i)%UserProp >= secondaryDepVal) upperTable = i
      end if
      
   end do
   
   if ( lowerTable == 0 .or. upperTable == 0 ) then
      if (interpolant == 1) then
         ErrMsg  = "Specified Reynold's number, "//trim(num2lstr(secondaryDepVal))//" , is outside the range of Re specified in the airfoil input file tables."
      else
         ErrMsg  = "Specified User Property's value, "//trim(num2lstr(secondaryDepVal))//" , is outside the range of User Property values specified in the airfoil input file tables."
      end if
      
      call SetErrStat ( ErrID_Fatal, ErrMsg, ErrStat, ErrMsg, 'AFI_ComputeAirfoilCoefs2DRe' )
      return
   end if
   
   IntAFCoefs = 0.0_ReKi ! initialize in case we only don't have 4 columns in the airfoil data (i.e., so cm is zero if not in the file)
 
      ! NOTE: we use Table(1) because the right now we can only interpolate with AOA and not Re or other variables.  If we had multiple tables stored
      ! for changes in other variables (Re, Mach #, etc) then then we would need to interpolate across tables.
      !
   s1 = size(p%Table(lowerTable)%Coefs,2)
   
   Alpha = AOA
   call MPi2Pi ( Alpha ) ! change AOA into range of -pi to pi
   
   
      ! Spline interpolation of lower table based on requested AOA
   
   IntAFCoefs(1:s1) = CubicSplineInterpM( Alpha  &
                                          , p%Table(lowerTable)%Alpha &
                                          , p%Table(lowerTable)%Coefs &
                                          , p%Table(lowerTable)%SplineCoefs &
                                          , ErrStat, ErrMsg )
   
  
   Cl1    = IntAFCoefs(1)
   Cd1    = IntAFCoefs(2)
   Cm1    = 0.0_Reki  !Set these to zero unless there is data to be read in
   Cpmin1 = 0.0_Reki
     
   IF ( p%ColCm > 0 ) Cm1 = IntAFCoefs(p%ColCm)
         
   IF ( p%ColCpmin > 0 ) Cpmin1 = IntAFCoefs(p%ColCpmin)
      
   
      ! Spline interpolation of upper table based on requested AOA
   
   IntAFCoefs(1:s1) = CubicSplineInterpM( Alpha  &
                                          , p%Table(upperTable)%Alpha &
                                          , p%Table(upperTable)%Coefs &
                                          , p%Table(upperTable)%SplineCoefs &
                                          , ErrStat, ErrMsg )
   
  
   Cl2    = IntAFCoefs(1)
   Cd2    = IntAFCoefs(2)
   Cm2    = 0.0_Reki  !Set these to zero unless there is data to be read in
   Cpmin2 = 0.0_Reki
     
   IF ( p%ColCm > 0 ) Cm2 = IntAFCoefs(p%ColCm)
         
   IF ( p%ColCpmin > 0 ) Cpmin2 = IntAFCoefs(p%ColCpmin)
   
   
      ! Linearly interpolate between tables
   if (interpolant == 1) then
      dVal = (p%Table(upperTable)%Re - p%Table(lowerTable)%Re)
      lowerVal = p%Table(lowerTable)%Re
   else
      dVal = (p%Table(upperTable)%UserProp - p%Table(lowerTable)%UserProp)
      lowerVal = p%Table(lowerTable)%UserProp
   end if
   
   if (EqualRealNos(dVal, 0.0_ReKi)) then
      s = 0.0_ReKi
   else                        
      s = (secondaryDepVal - lowerVal) / dVal
   end if
   
   Cl = Cl1 + s*Cl2
   Cd = Cd1 + s*Cd2
   Cm = Cm1 + s*Cm2
   Cpmin = Cpmin1 + s*Cpmin2
   
end subroutine AFI_ComputeAirfoilCoefs2D  
         
!----------------------------------------------------------------------------------------------------------------------------------  
subroutine AFI_ComputeAirfoilCoefs1D( AOA, p, &
                      Cl, Cd, Cm, Cpmin, errStat, errMsg )
! This routine is calculates Cl, Cd, Cm, (and Cpmin) for a set of tables which are dependent on AOA as well as a 2nd user-defined varible, could be Re or Cntrl, etc.
! If the requested yVar is not associated with a given table, then the two tables which contain yVar are found and, a cubic spline interpolation is performed at the requested AOA.
! for each of those two tables. Then a linear intepolation is performed on the 2nd dimension to find the final Cl,Cd,Cm, and Cpmin values.
! If the requested yVar corresponds to a table, then only a single cubic interpolation based on the requested AOA is performed.
!..................................................................................................................................
   real(ReKi),               intent(in   ) :: AOA 
   TYPE (AFI_ParameterType), intent(in   ) :: p                          ! This structure stores all the module parameters that are set by AirfoilInfo during the initialization phase.
   real(ReKi),               intent(  out) :: Cl, Cd, Cm, Cpmin
   integer(IntKi),           intent(  out) :: errStat                    ! Error status of the operation
   character(*),             intent(  out) :: errMsg                     ! Error message if ErrStat /= ErrID_None 
   
   
   real                            :: IntAFCoefs(4)                ! The interpolated airfoil coefficients.
   real(reki)                      :: Alpha
   integer                         :: s1
      
   ErrStat = ErrID_None
   ErrMsg  = ''

   IntAFCoefs = 0.0_ReKi ! initialize in case we only don't have 4 columns in the airfoil data (i.e., so cm is zero if not in the file)
 
   s1 = size(p%Table(1)%Coefs,2)
   
   Alpha = AOA
   call MPi2Pi ( Alpha ) ! change AOA into range of -pi to pi
   
   
      ! Spline interpolation of lower table based on requested AOA
   
   IntAFCoefs(1:s1) = CubicSplineInterpM( Alpha  &
                                          , p%Table(1)%Alpha &
                                          , p%Table(1)%Coefs &
                                          , p%Table(1)%SplineCoefs &
                                          , ErrStat, ErrMsg )
   
  
   Cl    = IntAFCoefs(1)
   Cd    = IntAFCoefs(2)
   Cm    = 0.0_Reki  !Set these to zero unless there is data to be read in
   Cpmin = 0.0_Reki
     
   IF ( p%ColCm > 0 ) Cm = IntAFCoefs(p%ColCm)
         
   IF ( p%ColCpmin > 0 ) Cpmin = IntAFCoefs(p%ColCpmin)
  
end subroutine AFI_ComputeAirfoilCoefs1D                        
                      
   subroutine AFI_GetAirfoilParams( p, M, Re, BL_p, C_nalpha_circ, errMsg, errStat )     

   type(AFI_ParameterType), intent(in   ) :: p                          ! This structure stores all the module parameters that are set by AirfoilInfo during the initialization phase.
   real(ReKi),              intent(in   ) :: M                             ! mach number
   real(ReKi),              intent(in   ) :: Re                            ! Reynold's number
   type(AFI_UA_BL_Type), intent(  out)       :: BL_p                          ! airfoil constants (UA BL parameters)
                          
   !real(ReKi),       intent(  out)       :: alpha0                        ! zero lift angle of attack (radians)
   !real(ReKi),       intent(  out)       :: alpha1                        ! angle of attack at f = 0.7, approximately the stall angle; for alpha >= alpha0 (radians)
   !real(ReKi),       intent(  out)       :: alpha2                        ! angle of attack at f = 0.7, approximately the stall angle; for alpha < alpha0 (radians)
   !real(ReKi),       intent(  out)       :: eta_e                         !
   !real(ReKi),       intent(  out)       :: C_nalpha                      !
   real(ReKi),              intent(  out) :: C_nalpha_circ                 ! slope of the circulatory normal force coefficient vs alpha curve
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
   integer(IntKi),          intent(  out) :: errStat                       ! Error status. 
   character(*),            intent(  out) :: errMsg                        ! Error message.
      
      errMsg         = ''
      errStat        = ErrID_None
        
      ! These coefs are stored in the p data structures based on Re
   BL_p%alpha0         =  p%Table(1)%UA_BL%alpha0   * D2R   ! Convert to radians
   BL_p%alpha1         =  p%Table(1)%UA_BL%alpha1   * D2R   ! Convert to radians
   BL_p%alpha2         =  p%Table(1)%UA_BL%alpha2   * D2R   ! Convert to radians
   BL_p%eta_e          =  p%Table(1)%UA_BL%eta_e         !0.90               ! Recovery factor in the range [0.85 - 0.95]
   BL_p%C_nalpha       =  p%Table(1)%UA_BL%C_nalpha      !2*pi  
   BL_p%T_f0           =  p%Table(1)%UA_BL%T_f0          !3.0_ReKi  ! seconds
   BL_p%T_V0           =  p%Table(1)%UA_BL%T_V0          !6.0_ReKi
   BL_p%T_p            =  p%Table(1)%UA_BL%T_p           !1.7_ReKi
   BL_p%T_VL           =  p%Table(1)%UA_BL%T_VL          !11.0_ReKi
   BL_p%b1             =  p%Table(1)%UA_BL%b1            !0.14_ReKi
   BL_p%b2             =  p%Table(1)%UA_BL%b2            !0.53_ReKi
   BL_p%b5             =  p%Table(1)%UA_BL%b5            !5.0_ReKi
   BL_p%A1             =  p%Table(1)%UA_BL%A1            !0.3_ReKi
   BL_p%A2             =  p%Table(1)%UA_BL%A2            !0.70_ReKi
   BL_p%A5             =  p%Table(1)%UA_BL%A5            !1.0_ReKi
   BL_p%S1             =  p%Table(1)%UA_BL%S1            !0.0262_ReKi   !!!!!!!!!!
   BL_p%S2             =  p%Table(1)%UA_BL%S2            !0.0201_ReKi   !!!!!!!!!!
   BL_p%S3             =  p%Table(1)%UA_BL%S3            !0.0262_ReKi   !!!!!!!!!!
   BL_p%S4             =  p%Table(1)%UA_BL%S4            !0.0201_ReKi   !!!!!!!!!!
   BL_p%Cn1            =  p%Table(1)%UA_BL%Cn1           !1.264_ReKi  ! Stall values of Cn
   BL_p%Cn2            =  p%Table(1)%UA_BL%Cn2           !-0.833_ReKi
   BL_p%St_sh          =  p%Table(1)%UA_BL%St_sh         !0.19_ReKi
   BL_p%Cd0            =  p%Table(1)%UA_BL%Cd0           !0.012_ReKi
   BL_p%Cm0            =  p%Table(1)%UA_BL%Cm0           !0.0_ReKi
   BL_p%k0             =  p%Table(1)%UA_BL%k0            !0.0_ReKi
   BL_p%k1             =  p%Table(1)%UA_BL%k1            !0.0_ReKi
   BL_p%k2             =  p%Table(1)%UA_BL%k2            !0.0_ReKi
   BL_p%k3             =  p%Table(1)%UA_BL%k3            !0.0_ReKi
   BL_p%k1_hat         =  p%Table(1)%UA_BL%k1_hat        !0.0_ReKi
   BL_p%x_cp_bar       =  p%Table(1)%UA_BL%x_cp_bar      !0.2_ReKi
   BL_p%filtCutOff     =  p%Table(1)%UA_BL%filtCutOff    ! 5.0_ReKi  Hz
   
   C_nalpha_circ  =  BL_p%C_nalpha / sqrt(1.0_ReKi-M**2)
     ! Cn1=1.9 Tp=1.7 Tf=3., Tv=6 Tvl=11, Cd0=0.012
   
   end subroutine AFI_GetAirfoilParams
                                 
!=============================================================================
END MODULE AirfoilInfo
