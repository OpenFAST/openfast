!*********************************************************************************************************************************
! LICENSING
! Copyright (C) 2015-2025  National Renewable Energy Laboratory
! Copyright (C) 2025-2026  National Laboratory of the Rockies
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
!*********************************************************************************************************************************
PROGRAM Test_ReadComFile


      ! This program tests the ScanComFile and ReadComFile routines in the NWTC Library.
      ! Created by Marshall Buhl.


   USE                                       :: NWTC_Library

   IMPLICIT                                     NONE


      ! Type declarations.

   INTEGER(IntKi)                            :: ErrStat   = 0                 ! Error status.

   INTEGER                                   :: Line                          ! Index into the arrays.

   CHARACTER(256)                            :: ErrMsg    = 'No Error'        ! Error message.
   CHARACTER(512)                            :: TopFileName                   ! The top-level file name.  It may include other files.

   TYPE (FileInfoType)                       :: FileInfo                      ! The derived type for holding the file information.



   CALL NWTC_Init ( 'Test_ReadComFile', 'v1.02.00, 23-Apr-2013' )
   CALL WrScr     ( NewLine//' Running Test_ReadComFile (v1.02.00, 23-Apr-2013)' )


      ! Get the name of a file with comments in it.

   CALL GET_COMMAND_ARGUMENT( 1, TopFileName, STATUS=ErrStat )
   IF ( ErrStat > 0 )  THEN
      CALL ExitThisRoutine ( ErrID_Fatal, ' Syntax: Test_ReadComFile <input_file>' )
      CALL ProgAbort ( ErrMsg )
   ENDIF


      ! Process the (possibly) nested set of files.

   CALL ProcessComFile ( TopFileName, FileInfo, ErrStat, ErrMsg )
   IF ( ErrStat > 0 )  THEN
      CALL ExitThisRoutine ( ErrID_Fatal, ErrMsg )
      CALL ProgAbort ( ErrMsg )
   ENDIF


      ! Print the decommented file(s) to the screen.

   PRINT '(/,A,/)', '  The decommented contents of "'//TRIM( FileInfo%FileList(1) )//'" and its included files are:'
   PRINT '(    A)', '  AggLine      FileName      FileLine         LineText'
   PRINT '(    A)', '  -------      --------      --------         --------'

   DO Line=1,FileInfo%NumLines
      PRINT '(5X,I2.2,4X,A16,5X,I2.2,6X,A)', Line, '"'//TRIM( FileInfo%FileList(FileInfo%FileIndx(Line)) )//'"', &
                                          FileInfo%FileLine(Line), '"'//TRIM( FileInfo%Lines(Line) )//'"'
   ENDDO ! Line

   CALL ExitThisRoutine( ErrID_None, '' )

   STOP

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

      END SUBROUTINE ExitThisRoutine ! ( ErrID, Msg )

END PROGRAM Test_ReadComFile
!=======================================================================
