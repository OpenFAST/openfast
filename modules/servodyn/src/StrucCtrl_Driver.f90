!*********************************************************************************************************************************
! StrucCtrl_Driver: This code tests the template modules
!..................................................................................................................................
! LICENSING
! Copyright (C) 2014 William La Cava & Matt Lackner, UMass Amherst
! Copyright (C) 2012 National Renewable Energy Laboratory
!
! This file is part of StrucCtrl. 
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
module read_file_module
! this module reads in external nacelle data for testing the module.  
  USE NWTC_Library
    USE StrucCtrl_Types
   implicit none
contains

   SUBROUTINE U_ReadInput(APvec,AVvec,AAvec,LAvec, NumSteps, ErrStat, ErrMsg )
! This subroutine reads the input file and stores all the data in the input vectors.
! It does not perform data validation.
!..................................................................................................................................

      ! Passed variables
   INTEGER(IntKi),           INTENT(IN)       :: NumSteps    ! Number of steps

   Real(ReKi), dimension(9, NumSteps), intent(inout) :: APvec
   Real(ReKi), dimension(3,NumSteps), intent(inout):: AVvec
   Real(ReKi), dimension(3,NumSteps), intent(inout):: AAvec
   Real(ReKi), dimension(3,NumSteps), intent(inout):: LAvec
   
   INTEGER(IntKi),       INTENT(OUT)      :: ErrStat        ! The error status code
   CHARACTER(*),         INTENT(OUT)      :: ErrMsg         ! The error message, if an error occurred

      ! local variables

   INTEGER(IntKi)                         :: UnEcho         ! Unit number for the echo file
   INTEGER(IntKi)                         :: ErrStat2       ! The error status code
   CHARACTER(LEN(ErrMsg))                 :: ErrMsg2        ! The error message, if an error occurred
   CHARACTER(1024) :: AV_file = 'AngVel_NO_Input_Data.inp'
   CHARACTER(1024) :: AA_file = 'AngAccel_NO_Input_Data.inp'
   CHARACTER(1024) :: AP_file = 'AngPos_NO_Input_Data.inp'
   CHARACTER(1024) :: RA_file = 'rddot_NO_Input_Data.inp'
    
      ! initialize values: 
   
   ErrStat = ErrID_None
   ErrMsg  = ""

     
      ! get the primary/platform input-file data
   !DO i = 1,NumSteps
   CALL ReadAngPosFile(AP_file, APvec, NumSteps,UnEcho, ErrStat2, ErrMsg2 )
      CALL CheckError(ErrStat2,ErrMsg2)
      IF ( ErrStat >= AbortErrLev ) RETURN
   CALL ReadAngVelFile( AV_file, AVvec,NumSteps,UnEcho, ErrStat2, ErrMsg2 )
      CALL CheckError(ErrStat2,ErrMsg2)
      IF ( ErrStat >= AbortErrLev ) RETURN
   CALL ReadAngAccelFile( AA_file, AAvec,NumSteps, UnEcho, ErrStat2, ErrMsg2 )
      CALL CheckError(ErrStat2,ErrMsg2)
      IF ( ErrStat >= AbortErrLev ) RETURN
   CALL ReadAccelFile( RA_file, LAvec,NumSteps, UnEcho, ErrStat2, ErrMsg2 )
      CALL CheckError(ErrStat2,ErrMsg2)
      IF ( ErrStat >= AbortErrLev ) RETURN
  ! END DO
      ! we may need to read additional files here (e.g., Bladed Interface)
   
      
      ! close any echo file that was opened
      
   IF ( UnEcho > 0 ) CLOSE( UnEcho )        

CONTAINS
   !...............................................................................................................................
   SUBROUTINE CheckError(ErrID,Msg)
   ! This subroutine sets the error message and level and cleans up if the error is >= AbortErrLev
   !...............................................................................................................................

         ! Passed arguments
      INTEGER(IntKi), INTENT(IN) :: ErrID       ! The error identifier (ErrStat)
      CHARACTER(*),   INTENT(IN) :: Msg         ! The error message (ErrMsg)


      !............................................................................................................................
      ! Set error status/message;
      !............................................................................................................................

      IF ( ErrID /= ErrID_None ) THEN

         IF (ErrStat /= ErrID_None) ErrMsg = TRIM(ErrMsg)//NewLine
         ErrMsg = TRIM(ErrMsg)//'StC_ReadInput:'//TRIM(Msg)
         ErrStat = MAX(ErrStat, ErrID)

         !.........................................................................................................................
         ! Clean up if we're going to return on error: close files, deallocate local arrays
         !.........................................................................................................................
         IF ( ErrStat >= AbortErrLev ) THEN
            IF ( UnEcho > 0 ) CLOSE( UnEcho )
         END IF

      END IF


   END SUBROUTINE CheckError     

END SUBROUTINE U_ReadInput
!----------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE ReadAngPosFile( InputFile, APvec, NumSteps, UnEc, ErrStat, ErrMsg )
! This routine reads in the nacelle angular position.
!   It opens and prints to an echo file if requested.
!..................................................................................................................................


   IMPLICIT                        NONE

      ! Passed variables
    INTEGER(IntKi),           INTENT(IN)       :: NumSteps    ! The default DT (from glue code)
   INTEGER(IntKi),     INTENT(OUT)     :: UnEc                                ! I/O unit for echo file. If > 0, file is open for writing.
   INTEGER(IntKi),     INTENT(OUT)     :: ErrStat                             ! Error status

   CHARACTER(*),       INTENT(IN)      :: InputFile                           ! Name of the file containing the primary input data
   CHARACTER(*),       INTENT(OUT)     :: ErrMsg                              ! Error message

   !TYPE(StC_InputFile), INTENT(INOUT) :: InputFileData                       ! All the data in the StC input file
    Real(ReKi), dimension(9, NumSteps), intent(inout) :: APvec
      ! Local variables:
   REAL(ReKi)                    :: TmpRAry9(9)      ! Temporary variable to read table from file
   INTEGER                :: I                                         ! loop counter
   INTEGER                :: J                                       ! loop counter
   INTEGER(IntKi)                :: UnIn                                      ! Unit number for reading file
     
   INTEGER(IntKi)                :: ErrStat2                                  ! Temporary Error status
   LOGICAL                       :: Echo                                      ! Determines if an echo file should be written
   CHARACTER(LEN(ErrMsg))        :: ErrMsg2                                   ! Temporary Error message
   CHARACTER(1024)               :: PriPath                                   ! Path name of the primary file

   
      ! Initialize some variables:
   ErrStat = ErrID_None
   ErrMsg  = ""
      
   UnEc = -1
   Echo = .FALSE.   
   CALL GetPath( InputFile, PriPath )     ! Input files will be relative to the path where the primary input file is located.
   

   !CALL AllocAry( InputFileData%OutList, MaxOutPts, "ServoDyn Input File's Outlist", ErrStat2, ErrMsg2 )
   !   CALL CheckError( ErrStat2, ErrMsg2 )
   !   IF ( ErrStat >= AbortErrLev ) RETURN   
      
   
      ! Get an available unit number for the file.

   CALL GetNewUnit( UnIn, ErrStat2, ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! Open the Primary input file.

   CALL OpenFInpFile ( UnIn, InputFile, ErrStat2, ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN
                  
      
   ! Read the lines up/including to the "Echo" simulation control variable
   ! If echo is FALSE, don't write these lines to the echo file. 
   ! If Echo is TRUE, rewind and write on the second try.
   
   I = 1 !set the number of times we've read the file
   !
   DO I=1,NumSteps
      CALL ReadAry( UnIn, InputFile, TmpRAry9, 9, 'AngPos_NO', 'Nacelle Rotation Matrix', ErrStat2, ErrMsg2, UnEc )

         DO J = 1,9
           APvec(J,I) = TmpRAry9(J)
         END DO
   END DO
      
   CLOSE ( UnIn )
   RETURN


CONTAINS
   !...............................................................................................................................
   SUBROUTINE CheckError(ErrID,Msg)
   ! This subroutine sets the error message and level
   !...............................................................................................................................

         ! Passed arguments
      INTEGER(IntKi), INTENT(IN) :: ErrID       ! The error identifier (ErrStat)
      CHARACTER(*),   INTENT(IN) :: Msg         ! The error message (ErrMsg)


      !............................................................................................................................
      ! Set error status/message;
      !............................................................................................................................

      IF ( ErrID /= ErrID_None ) THEN

         IF (ErrStat /= ErrID_None) ErrMsg = TRIM(ErrMsg)//NewLine
         ErrMsg = TRIM(ErrMsg)//'ReadPrimaryFile:'//TRIM(Msg)
         ErrStat = MAX(ErrStat, ErrID)

         !.........................................................................................................................
         ! Clean up if we're going to return on error: close file, deallocate local arrays
         !.........................................................................................................................
         IF ( ErrStat >= AbortErrLev ) THEN
            CLOSE( UnIn )
!            IF ( UnEc > 0 ) CLOSE ( UnEc )
         END IF

      END IF


   END SUBROUTINE CheckError
   !...............................................................................................................................
   END SUBROUTINE ReadAngPosFile      
   SUBROUTINE ReadAngVelFile( InputFile, AVvec, NumSteps, UnEc, ErrStat, ErrMsg )
! This routine reads in the nacelle angular velocity.
!   It opens and prints to an echo file if requested.
!..................................................................................................................................


   IMPLICIT                        NONE

      ! Passed variables
    INTEGER(IntKi),           INTENT(IN)       :: NumSteps    ! The default DT (from glue code)
   INTEGER(IntKi),     INTENT(OUT)     :: UnEc                                ! I/O unit for echo file. If > 0, file is open for writing.
   INTEGER(IntKi),     INTENT(OUT)     :: ErrStat                             ! Error status

   CHARACTER(*),       INTENT(IN)      :: InputFile                           ! Name of the file containing the primary input data
   CHARACTER(*),       INTENT(OUT)     :: ErrMsg                              ! Error message

   !TYPE(StC_InputFile), INTENT(INOUT) :: InputFileData                       ! All the data in the StC input file
   Real(ReKi), dimension(3,NumSteps), intent(inout) :: AVvec
      ! Local variables:
   REAL(ReKi)                    :: TmpRAry3(3)      ! Temporary variable to read table from file
   INTEGER(IntKi)                :: I                                         ! loop counter
      INTEGER(IntKi)                :: J                                         ! loop counter
   INTEGER(IntKi)                :: UnIn                                      ! Unit number for reading file
     
   INTEGER(IntKi)                :: ErrStat2                                  ! Temporary Error status
   LOGICAL                       :: Echo                                      ! Determines if an echo file should be written
   CHARACTER(LEN(ErrMsg))        :: ErrMsg2                                   ! Temporary Error message
   CHARACTER(1024)               :: PriPath                                   ! Path name of the primary file

   
      ! Initialize some variables:
   ErrStat = ErrID_None
   ErrMsg  = ""
      
   UnEc = -1
   Echo = .FALSE.   
   CALL GetPath( InputFile, PriPath )     ! Input files will be relative to the path where the primary input file is located.
   

   !CALL AllocAry( InputFileData%OutList, MaxOutPts, "ServoDyn Input File's Outlist", ErrStat2, ErrMsg2 )
   !   CALL CheckError( ErrStat2, ErrMsg2 )
   !   IF ( ErrStat >= AbortErrLev ) RETURN   
      
   
      ! Get an available unit number for the file.

   CALL GetNewUnit( UnIn, ErrStat2, ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! Open the Primary input file.

   CALL OpenFInpFile ( UnIn, InputFile, ErrStat2, ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN
                  
      
   ! Read the lines up/including to the "Echo" simulation control variable
   ! If echo is FALSE, don't write these lines to the echo file. 
   ! If Echo is TRUE, rewind and write on the second try.
   
   I = 1 !set the number of times we've read the file
   !
   DO I=1,NumSteps
      CALL ReadAry( UnIn, InputFile, TmpRAry3, 3, 'AngVel_NO', 'Nacelle Angular Velocity', ErrStat2, ErrMsg2, UnEc )
      DO J = 1,3
         AVvec(J,I) = TmpRAry3(J)
      END DO
   END DO
   
   CLOSE ( UnIn )
   RETURN


CONTAINS
   !...............................................................................................................................
   SUBROUTINE CheckError(ErrID,Msg)
   ! This subroutine sets the error message and level
   !...............................................................................................................................

         ! Passed arguments
      INTEGER(IntKi), INTENT(IN) :: ErrID       ! The error identifier (ErrStat)
      CHARACTER(*),   INTENT(IN) :: Msg         ! The error message (ErrMsg)


      !............................................................................................................................
      ! Set error status/message;
      !............................................................................................................................

      IF ( ErrID /= ErrID_None ) THEN

         IF (ErrStat /= ErrID_None) ErrMsg = TRIM(ErrMsg)//NewLine
         ErrMsg = TRIM(ErrMsg)//'ReadPrimaryFile:'//TRIM(Msg)
         ErrStat = MAX(ErrStat, ErrID)

         !.........................................................................................................................
         ! Clean up if we're going to return on error: close file, deallocate local arrays
         !.........................................................................................................................
         IF ( ErrStat >= AbortErrLev ) THEN
            CLOSE( UnIn )
!            IF ( UnEc > 0 ) CLOSE ( UnEc )
         END IF

      END IF


   END SUBROUTINE CheckError
   !...............................................................................................................................
   END SUBROUTINE ReadAngVelFile     
   SUBROUTINE ReadAngAccelFile( InputFile, AAvec, NumSteps, UnEc, ErrStat, ErrMsg )
! This routine reads in the nacelle angular acceleration.
!   It opens and prints to an echo file if requested.
!..................................................................................................................................


   IMPLICIT                        NONE

      ! Passed variables
    INTEGER(IntKi),           INTENT(IN)       :: NumSteps    ! The default DT (from glue code)
   INTEGER(IntKi),     INTENT(OUT)     :: UnEc                                ! I/O unit for echo file. If > 0, file is open for writing.
   INTEGER(IntKi),     INTENT(OUT)     :: ErrStat                             ! Error status

   CHARACTER(*),       INTENT(IN)      :: InputFile                           ! Name of the file containing the primary input data
   CHARACTER(*),       INTENT(OUT)     :: ErrMsg                              ! Error message

   !TYPE(StC_InputFile), INTENT(INOUT) :: InputFileData                       ! All the data in the StC input file
    Real(ReKi), dimension(3,NumSteps), intent(inout) :: AAvec
      ! Local variables:
   REAL(ReKi)                    :: TmpRAry3(3)      ! Temporary variable to read table from file
   INTEGER(IntKi)                :: I                                         ! loop counter
   INTEGER(IntKi)                :: J                                         ! loop counter
   INTEGER(IntKi)                :: UnIn                                      ! Unit number for reading file
     
   INTEGER(IntKi)                :: ErrStat2                                  ! Temporary Error status
   LOGICAL                       :: Echo                                      ! Determines if an echo file should be written
   CHARACTER(LEN(ErrMsg))        :: ErrMsg2                                   ! Temporary Error message
   CHARACTER(1024)               :: PriPath                                   ! Path name of the primary file

   
      ! Initialize some variables:
   ErrStat = ErrID_None
   ErrMsg  = ""
      
   UnEc = -1
   Echo = .FALSE.   
   CALL GetPath( InputFile, PriPath )     ! Input files will be relative to the path where the primary input file is located.
   

   !CALL AllocAry( InputFileData%OutList, MaxOutPts, "ServoDyn Input File's Outlist", ErrStat2, ErrMsg2 )
   !   CALL CheckError( ErrStat2, ErrMsg2 )
   !   IF ( ErrStat >= AbortErrLev ) RETURN   
      
   
      ! Get an available unit number for the file.

   CALL GetNewUnit( UnIn, ErrStat2, ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! Open the Primary input file.

   CALL OpenFInpFile ( UnIn, InputFile, ErrStat2, ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN
                  
      
   ! Read the lines up/including to the "Echo" simulation control variable
   ! If echo is FALSE, don't write these lines to the echo file. 
   ! If Echo is TRUE, rewind and write on the second try.
   
   I = 1 !set the number of times we've read the file
   !
   DO I = 1,NumSteps
      CALL ReadAry( UnIn, InputFile, TmpRAry3, 3, 'AngAccel_NO', 'Nacelle Angular Acceleration', ErrStat2, ErrMsg2, UnEc )
      DO J = 1,3
         AAvec(J,I) = TmpRAry3(J)
      END DO
   END DO
   CLOSE ( UnIn )
   RETURN


CONTAINS
   !...............................................................................................................................
   SUBROUTINE CheckError(ErrID,Msg)
   ! This subroutine sets the error message and level
   !...............................................................................................................................

         ! Passed arguments
      INTEGER(IntKi), INTENT(IN) :: ErrID       ! The error identifier (ErrStat)
      CHARACTER(*),   INTENT(IN) :: Msg         ! The error message (ErrMsg)


      !............................................................................................................................
      ! Set error status/message;
      !............................................................................................................................

      IF ( ErrID /= ErrID_None ) THEN

         IF (ErrStat /= ErrID_None) ErrMsg = TRIM(ErrMsg)//NewLine
         ErrMsg = TRIM(ErrMsg)//'ReadPrimaryFile:'//TRIM(Msg)
         ErrStat = MAX(ErrStat, ErrID)

         !.........................................................................................................................
         ! Clean up if we're going to return on error: close file, deallocate local arrays
         !.........................................................................................................................
         IF ( ErrStat >= AbortErrLev ) THEN
            CLOSE( UnIn )
!            IF ( UnEc > 0 ) CLOSE ( UnEc )
         END IF

      END IF


   END SUBROUTINE CheckError
   !...............................................................................................................................
END SUBROUTINE ReadAngAccelFile      
   SUBROUTINE ReadAccelFile( InputFile, LAvec, NumSteps, UnEc, ErrStat, ErrMsg )
! This routine reads in the nacelle translational acceleration.
!   It opens and prints to an echo file if requested.
!..................................................................................................................................


   IMPLICIT                        NONE

      ! Passed variables
    INTEGER(IntKi),           INTENT(IN)       :: NumSteps    ! The default DT (from glue code)
   INTEGER(IntKi),     INTENT(OUT)     :: UnEc                                ! I/O unit for echo file. If > 0, file is open for writing.
   INTEGER(IntKi),     INTENT(OUT)     :: ErrStat                             ! Error status

   CHARACTER(*),       INTENT(IN)      :: InputFile                           ! Name of the file containing the primary input data
   CHARACTER(*),       INTENT(OUT)     :: ErrMsg                              ! Error message 

   !TYPE(StC_InputFile), INTENT(INOUT) :: InputFileData                       ! All the data in the StC input file
    Real(ReKi), dimension(3,NumSteps), intent(inout) :: LAvec
      ! Local variables:
   REAL(ReKi)                    :: TmpRAry3(3)      ! Temporary variable to read table from file
   INTEGER(IntKi)                :: I                                         ! loop counter
   INTEGER(IntKi)                :: J                                        ! loop counter
   INTEGER(IntKi)                :: UnIn                                      ! Unit number for reading file
     
   INTEGER(IntKi)                :: ErrStat2                                  ! Temporary Error status
   LOGICAL                       :: Echo                                      ! Determines if an echo file should be written
   CHARACTER(LEN(ErrMsg))        :: ErrMsg2                                   ! Temporary Error message
   CHARACTER(1024)               :: PriPath                                   ! Path name of the primary file

   
      ! Initialize some variables:
   ErrStat = ErrID_None
   ErrMsg  = ""
      
   UnEc = -1
   Echo = .FALSE.   
   CALL GetPath( InputFile, PriPath )     ! Input files will be relative to the path where the primary input file is located.
   

   !CALL AllocAry( InputFileData%OutList, MaxOutPts, "ServoDyn Input File's Outlist", ErrStat2, ErrMsg2 )
   !   CALL CheckError( ErrStat2, ErrMsg2 )
   !   IF ( ErrStat >= AbortErrLev ) RETURN   
      
   
      ! Get an available unit number for the file.

   CALL GetNewUnit( UnIn, ErrStat2, ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! Open the Primary input file.

   CALL OpenFInpFile ( UnIn, InputFile, ErrStat2, ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN
                  
      
   ! Read the lines up/including to the "Echo" simulation control variable
   ! If echo is FALSE, don't write these lines to the echo file. 
   ! If Echo is TRUE, rewind and write on the second try.
   
  ! I = 1 !set the number of times we've read the file
   !
   DO I=1,NumSteps
      CALL ReadAry( UnIn, InputFile, TmpRAry3, 3, 'rddot_NO', 'Nacelle Linear Acceleration', ErrStat2, ErrMsg2, UnEc )

      DO J = 1,3
         LAvec(J,I) = TmpRAry3(J)
      END DO
   END DO
   
   CLOSE ( UnIn )
   RETURN


CONTAINS
   !...............................................................................................................................
   SUBROUTINE CheckError(ErrID,Msg)
   ! This subroutine sets the error message and level
   !...............................................................................................................................

         ! Passed arguments
      INTEGER(IntKi), INTENT(IN) :: ErrID       ! The error identifier (ErrStat)
      CHARACTER(*),   INTENT(IN) :: Msg         ! The error message (ErrMsg)


      !............................................................................................................................
      ! Set error status/message;
      !............................................................................................................................

      IF ( ErrID /= ErrID_None ) THEN

         IF (ErrStat /= ErrID_None) ErrMsg = TRIM(ErrMsg)//NewLine
         ErrMsg = TRIM(ErrMsg)//'ReadPrimaryFile:'//TRIM(Msg)
         ErrStat = MAX(ErrStat, ErrID)

         !.........................................................................................................................
         ! Clean up if we're going to return on error: close file, deallocate local arrays
         !.........................................................................................................................
         IF ( ErrStat >= AbortErrLev ) THEN
            CLOSE( UnIn )
!            IF ( UnEc > 0 ) CLOSE ( UnEc )
         END IF

      END IF


   END SUBROUTINE CheckError
   !...............................................................................................................................
END SUBROUTINE ReadAccelFile     
 

!-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE StC_OpenOutputFile(OutputFile,UnIn,ErrStat,ErrMsg)
! This routine is called by the driver, not this module.
   CHARACTER(1024), Intent(IN)         :: OutputFile    ! Name of the file containing the primary input data
   INTEGER(IntKi), INTENT(OUT)         :: UnIn          ! Unit number for writing file
   INTEGER(IntKi), INTENT(OUT)         :: ErrStat       ! Temporary error ID   
   CHARACTER(*), INTENT(OUT)           :: ErrMsg        ! Temporary message describing error
   CHARACTER(1024)                     :: Header1   
   CHARACTER(1024)                     :: Header2 
   
   ErrStat = ErrID_None
   ErrMsg  = ''
   !OutputFile = 'StC_Output_Data.txt'
   !Fmt = "F10.2))/"
   
   CALL GetNewUnit( UnIn, ErrStat, ErrMsg )
      !CALL CheckError( ErrStat, ErrMsg)
      !IF ( ErrStat >= AbortErrLev ) RETURN


      ! Open the output file.

   CALL OpenFOutFile ( UnIn, OutputFile, ErrStat, ErrMsg )
   Header1 = "-------------- StrucCtrl Output ------------------------------"
   Header2 = "x    dxdt     y     dydt       fx       fy     fz       mx       my       mz"
   
    WRITE( UnIn, *, IOSTAT=ErrStat ) TRIM(Header1)
     WRITE( UnIn, *, IOSTAT=ErrStat ) TRIM(Header2)
     
END SUBROUTINE StC_OpenOutputFile
!-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE StC_CloseOutputFile(Un)
! This routine is called by the driver, not this module.
 
 INTEGER(IntKi), INTENT(IN)                :: Un                                      ! Unit number for writing file
  CLOSE ( Un )
END SUBROUTINE StC_CloseOutputFile
!-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE StC_WriteOutputFile( x, y, UnIn, ErrStat, ErrMsg )
! This routine is called by the driver, not this module.
! write output file with StC states and forces.

   TYPE(StC_ContinuousStateType), INTENT(IN   )    :: x           ! Continuous states at Time
   TYPE(StC_OutputType), INTENT(IN   )             :: y                   ! state outputs
   INTEGER(IntKi), INTENT(IN)                :: UnIn                                      ! Unit number for writing file
   INTEGER(IntKi), INTENT(OUT)                     :: ErrStat       ! Temporary error ID   
   CHARACTER(*), INTENT(OUT)            :: ErrMsg        ! Temporary message describing error
   !REAL(DbKi),                      INTENT(IN   )  :: Time        ! Current simulation time in seconds
   
    CHARACTER(1024)              :: Fmt !text format
    REAL(ReKi), dimension(10)                   :: OutAry
    INTEGER(IntKi)                              :: i
    INTEGER(IntKi)                              :: i_pt  ! index into mesh point
    ErrStat = ErrID_None
   ErrMsg  = ''
 
!FIXME: allow different sizes for StC_x second dimension -- loop over i_pt
!FIXME: allow for different size meshes                  -- loop over i_pt
   i_pt=1

   ! create output array
   DO i=1,4
      OutAry(i) = x%StC_x(i,i_pt)
   END DO
   DO i=5,7
      OutAry(i) = y%Mesh(i_pt)%Force(i-4,1)
   END DO
   DO i=8,10
      OutAry(i) = y%Mesh(i_pt)%Moment(i-7,1)
   END DO
   !Write output 
    Fmt = '(10(1x,F10.2))'
   WRITE( UnIn, Fmt, IOSTAT=ErrStat ) OutAry(:)
   IF (ErrStat /= 0) THEN
      CALL WrScr('Error '//TRIM(Num2LStr(ErrStat))//' writing matrix in WrMatrix1R4().')
      RETURN
   END IF
   !CALL WrMatrix( x%StC_x, UnIn, Fmt )   
   
END SUBROUTINE StC_WriteOutputFile
!-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


end module read_file_module 

PROGRAM StrucCtrl_Driver
    
    USE NWTC_Library
    USE StrucCtrl
    USE StrucCtrl_Types
    USE read_file_module
    
    IMPLICIT NONE
    
    INTEGER(IntKi), PARAMETER :: NumInp = 2 ! Number of inputs sent to StC_UpdateStates
    INTEGER(IntKi), PARAMETER :: NumSteps = 100 ! Number of time steps
        ! Program variables
    REAL(DbKi) :: Time ! Variable for storing time, in seconds
    REAL(DbKi) :: TimeInterval ! Interval between time steps, in seconds
    REAL(DbKi) :: InputTime(NumInp) ! Variable for storing time associated with inputs, in seconds
   
    TYPE(StC_InitInputType) :: InitInData ! Input data for initialization
    TYPE(StC_InitOutputType) :: InitOutData ! Output data from initialization
    
    TYPE(StC_ContinuousStateType) :: x ! Continuous states
    TYPE(StC_DiscreteStateType) :: xd ! Discrete states
    TYPE(StC_ConstraintStateType) :: z ! Constraint states
    TYPE(StC_ConstraintStateType) :: Z_residual ! Residual of the constraint state functions (Z)
    TYPE(StC_OtherStateType) :: OtherState ! Other states
    TYPE(StC_MiscVarType) :: m ! misc variables
    
    TYPE(StC_ParameterType) :: p ! Parameters
    TYPE(StC_InputType) :: u(NumInp) ! System inputs
    TYPE(StC_OutputType) :: y ! System outputs
    
    TYPE(StC_ContinuousStateType) :: dxdt ! First time derivatives of the continuous states
   integer(IntKi)             :: UnOut !output data file number
   

    INTEGER(IntKi) :: n ! Loop counter (for time step)
    INTEGER(IntKi) :: i ! Loop counter (for time step)
    INTEGER(IntKi) :: j ! Loop counter (for time step)
    INTEGER(IntKi) :: count ! Loop counter (for time step)
    INTEGER(IntKi) :: ErrStat ! Status of error message
    CHARACTER(1024) :: ErrMsg ! Error message if ErrStat /= ErrID_None
    
    REAL(ReKi), ALLOCATABLE :: Re_SaveAry (:) ! Array to store reals in packed data structure
    REAL(DbKi), ALLOCATABLE :: Db_SaveAry (:) ! Array to store doubles in packed data structure
    INTEGER(IntKi), ALLOCATABLE :: Int_SaveAry (:) ! Array to store integers in packed data structure
    Real(ReKi), dimension(9,NumSteps) :: APvec
    Real(ReKi), dimension(3,NumSteps) :: AVvec
    Real(ReKi), dimension(3,NumSteps) :: AAvec
    Real(ReKi), dimension(3,NumSteps) :: LAvec
    CHARACTER(1024)              :: OutputName !text file output

    integer(IntKi)                     :: i_pt     ! index counter to points 
!...............................................................................................................................
! Routines called in initialization
!...............................................................................................................................
    ! Populate the InitInData data structure here:
    ! input file with StC settings
    InitInData%InputFile = 'StC_Input_test.dat'
    ! gravity
    InitInData%Gravity = 9.80665
    ! StC origin and orientation
    call AllocAry(InitInData%InitPosition,       3, 1, 'InitPosition',      ErrStat,ErrMsg)
    call AllocAry(InitInData%InitOrientation, 3, 3, 1, 'InitOrientation',   ErrStat,ErrMsg)
    InitInData%InitPosition(1:3,1) = (/ 0.0_ReKi, 0.0_ReKi, 0.0_ReKi /)
    InitInData%InitOrientation = 0.0_R8Ki
    do i=1,3
      InitInData%InitOrientation(i,i,1) = 1.0_R8Ki
    enddo
    
    ! Set the driver's request for time interval here:
    
    TimeInterval = 0.25 ! Glue code's request for delta time (likely based on information from other modules)

    ! Initialize the module
    
    CALL StC_Init( InitInData, u(1), p, x, xd, z, OtherState, y, m, TimeInterval, InitOutData, ErrStat, ErrMsg )
    IF ( ErrStat /= ErrID_None ) THEN ! Check if there was an error and do something about it if necessary
        IF (ErrStat >= AbortErrLev) CALL ProgAbort( ErrMsg )
        CALL WrScr( ErrMsg )
    END IF
    CALL StC_CopyInput( u(1), u(2), MESH_NEWCOPY, ErrStat, ErrMsg )
  
    ! read in nacelle data from file
    CALL U_ReadInput(APvec,AVvec,AAvec,LAvec, NumSteps, ErrStat, ErrMsg )
    
    ! Destroy initialization data
    CALL StC_DestroyInitInput( InitInData, ErrStat, ErrMsg )
    CALL StC_DestroyInitOutput( InitOutData, ErrStat, ErrMsg )
    !...............................................................................................................................
    ! Routines called in loose coupling -- the glue code may implement this in various ways
    !....................................................................................................
    ! setup the output file:
    OutputName = 'StC_Output_Data.txt'
    CALL StC_OpenOutputFile(OutputName,UnOut,ErrStat,ErrMsg)
    
    ! run simulation 

   !FIXME: allow for more than one point?
   i_pt = 1    ! index counter of number of points we are simulating

DO n = 0,NumSteps-1
       count=1
       ! Modify u (likely from the outputs of another module or a set of test conditions) here:
       IF (n>0) THEN
         CALL StC_CopyInput( u(2), u(1), MESH_UPDATECOPY, ErrStat, ErrMsg )
       !  u(1) = u(2) ! save past input as first element in input vector
       END IF
       i=1
       j=1
       ! setup input mesh with data from nacelle positions:
       do i = 1,3
          do j=1,3
               u(2)%Mesh(i_pt)%Orientation(i,j,1) = APvec(count,n+1)
               count = count+1
          end do
          u(2)%Mesh(i_pt)%RotationVel(i,1) = AVvec(i,n+1)
          u(2)%Mesh(i_pt)%RotationAcc(i,1) = AAvec(i,n+1)
          u(2)%Mesh(i_pt)%TranslationAcc(i,1) = LAvec(i,n+1)
       end do
        if (n==0) then
          InputTime(1) = 0
          InputTime(2) = TimeInterval
        else 
           InputTime(1) = Time
           Time = n*TimeInterval
           InputTime(2) = Time
        end if
        
        
        ! Calculate outputs at n
        CALL StC_CalcOutput( Time, u(1), p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )
        IF ( ErrStat /= ErrID_None ) THEN ! Check if there was an error and do something about it if necessary
            CALL WrScr( ErrMsg )
        END IF
        ! Get state variables at next step: INPUT at step n, OUTPUT at step n + 1
        CALL StC_UpdateStates( Time, n, u, InputTime, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
        IF ( ErrStat /= ErrID_None ) THEN ! Check if there was an error and do something about it if necessary
            CALL WrScr( ErrMsg )
        END IF
         
        ! write outputs to file
       CALL StC_WriteOutputFile(x, y, UnOut,ErrStat,ErrMsg)
END DO
! close the output file
CALL StC_CloseOutputFile(UnOut)
    !...............................................................................................................................
    ! Routines called in tight coupling -- time marching only
    !...............................................................................................................................
    !DO n = 0,10
    !    Time = n * TimeInterval ! Note that the discrete states must be updated only at the TimeInterval defined in initialization
    !    ! set inputs (u) here:
    !    ! u =
    !    ! Update constraint states at Time
    !    ! DO
    !    !CALL StC_CalcConstrStateResidual( Time, u(1), p, x, xd, z, OtherState, Z_residual, ErrStat, ErrMsg )
    !    !
    !    !IF ( ErrStat /= ErrID_None ) THEN ! Check if there was an error and do something about it if necessary
    !    !    CALL WrScr( ErrMsg )
    !    !END IF
    !    ! z =
    !    ! END DO
    !    ! Calculate the outputs at Time
    !    CALL StC_CalcOutput( Time, u(1), p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
    !    
    !    IF ( ErrStat /= ErrID_None ) THEN ! Check if there was an error and do something about it if necessary
    !        CALL WrScr( ErrMsg )
    !    END IF
    !    ! Calculate the continuous state derivatives at Time
    !    CALL StC_CalcContStateDeriv( Time, u(1), p, x, xd, z, OtherState, dxdt, ErrStat, ErrMsg )
    !    
    !    IF ( ErrStat /= ErrID_None ) THEN ! Check if there was an error and do something about it if necessary
    !        CALL WrScr( ErrMsg )
    !    END IF
    !    ! Update the discrete state from step n to step n+1
    !    ! Note that the discrete states must be updated only at the TimeInterval defined in initialization
    !    !CALL StC_UpdateDiscState( Time, n, u(1), p, x, xd, z, OtherState, ErrStat, ErrMsg )
    !    !
    !    !IF ( ErrStat /= ErrID_None ) THEN ! Check if there was an error and do something about it if necessary
    !    !    CALL WrScr( ErrMsg )
    !    !END IF
    !    ! Driver should integrate (update) continuous states here:
    !    !x = function of dxdt, x
    !    ! Jacobians required:
    !    !CALL StC_JacobianPInput( Time, u(1), p, x, xd, z, OtherState, dYdu=dYdu, dZdu=dZdu, ErrStat=ErrStat, ErrMsg=ErrMsg )
    !    !
    !    !IF ( ErrStat /= ErrID_None ) THEN ! Check if there was an error and do something about it if necessary
    !    !    CALL WrScr( ErrMsg )
    !    !END IF
    !    !
    !    !CALL StC_JacobianPConstrState( Time, u(1), p, x, xd, z, OtherState, dYdz=dYdz, dZdz=dZdz, &
    !    !ErrStat=ErrStat, ErrMsg=ErrMsg )
    !    !
    !    !IF ( ErrStat /= ErrID_None ) THEN ! Check if there was an error and do something about it if necessary
    !    !    CALL WrScr( ErrMsg )
    !    !END IF
    !END DO
    ! Destroy Z_residual and dxdt because they are not necessary anymore
    CALL StC_DestroyConstrState( Z_residual, ErrStat, ErrMsg )
    
    IF ( ErrStat /= ErrID_None ) THEN ! Check if there was an error and do something about it if necessary
        CALL WrScr( ErrMsg )
    END IF
    
    CALL StC_DestroyContState( dxdt, ErrStat, ErrMsg )
    
    IF ( ErrStat /= ErrID_None ) THEN ! Check if there was an error and do something about it if necessary
        CALL WrScr( ErrMsg )
    END IF
    
    !...............................................................................................................................
    ! Jacobian routines called in tight coupling
    !...............................................................................................................................
    !CALL StC_JacobianPInput( Time, u(1), p, x, xd, z, OtherState, dYdu, dXdu, dXddu, dZdu, ErrStat, ErrMsg )
    !
    !IF ( ErrStat /= ErrID_None ) THEN ! Check if there was an error and do something about it if necessary
    !END IF
    !
    !CALL StC_JacobianPContState( Time, u(1), p, x, xd, z, OtherState, dYdx, dXdx, dXddx, dZdx, ErrStat, ErrMsg )
    !
    !IF ( ErrStat /= ErrID_None ) THEN ! Check if there was an error and do something about it if necessary
    !    CALL WrScr( ErrMsg )
    !END IF
    !
    !CALL StC_JacobianPDiscState( Time, u(1), p, x, xd, z, OtherState, dYdxd, dXdxd, dXddxd, dZdxd, ErrStat, ErrMsg )
    !
    !IF ( ErrStat /= ErrID_None ) THEN ! Check if there was an error and do something about it if necessary
    !    CALL WrScr( ErrMsg )
    !END IF
    !
    !CALL StC_JacobianPConstrState( Time, u(1), p, x, xd, z, OtherState, dYdz, dXdz, dXddz, dZdz, ErrStat, ErrMsg )
    !
    !IF ( ErrStat /= ErrID_None ) THEN ! Check if there was an error and do something about it if necessary
    !    CALL WrScr( ErrMsg )
    !END IF
    !...............................................................................................................................
    ! Routines to pack data (to restart later)
    !...............................................................................................................................
    !CALL StC_Pack(Re_SaveAry, Db_SaveAry, Int_SaveAry, u(1), p, x, xd, z, OtherState, y, ErrStat, ErrMsg)
    
    IF ( ErrStat /= ErrID_None ) THEN
        CALL WrScr( ErrMsg )
    END IF
    !...............................................................................................................................
    ! Routine to terminate program execution
    !...............................................................................................................................
    CALL StC_End( u(1), p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )
    
    IF ( ErrStat /= ErrID_None ) THEN
        CALL WrScr( ErrMsg )
    END IF
    !...............................................................................................................................
    ! Routines to retreive packed data (unpack for restart)
    !...............................................................................................................................
    !CALL StC_Unpack( Re_SaveAry, Db_SaveAry, Int_SaveAry, u(1), p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
    
    IF ( ErrStat /= ErrID_None ) THEN
        CALL WrScr( ErrMsg )
    END IF
    !...............................................................................................................................
    ! Routines to copy data (not already tested)
    !...............................................................................................................................
    !...............................................................................................................................
    ! Routines to destroy data (not already tested)
    !...............................................................................................................................
    IF ( ALLOCATED( Re_SaveAry ) ) DEALLOCATE( Re_SaveAry )
    IF ( ALLOCATED( Db_SaveAry ) ) DEALLOCATE( Db_SaveAry )
    IF ( ALLOCATED( Int_SaveAry ) ) DEALLOCATE( Int_SaveAry )
    ! CALL StC_DestroyPartialOutputPInput ( ) ! Jacobian Routine not yet implemented
    !...............................................................................................................................
    ! Routine to terminate program execution (again)
    !...............................................................................................................................
    CALL StC_End( u(1), p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )
    IF ( ErrStat /= ErrID_None ) THEN
        CALL WrScr( ErrMsg )
    END IF
    
END PROGRAM StrucCtrl_Driver
