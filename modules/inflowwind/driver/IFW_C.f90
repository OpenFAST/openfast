!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2021 Nicole Mendoza
!
! This file is part of InflowWind.
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
MODULE InflowWindAPI

USE ISO_C_BINDING
USE pFUnit_mod
USE InflowWind_Subs                !<-- NEED FOR InflowWind_ParseInputFileInfo (Fortran)

IMPLICIT NONE

PUBLIC :: IFW_INIT_C
PUBLIC :: IFW_CALCOUTPUT_C
PUBLIC :: IFW_END_C

CONTAINS

SUBROUTINE IFW_INIT_C(input_file_string,ErrStat_C,ErrMsg_C) BIND (C, NAME='IFW_INIT_C')

    TYPE(C_PTR)                   , INTENT(IN   )      :: input_file_string_c_ptr(*)
    CHARACTER(KIND=C_CHAR),POINTER, INTENT(INOUT)      :: input_file_string_c(:)
    INTEGER(C_INT),               , INTENT(  OUT)      :: ErrStat_C
    CHARACTER(KIND=C_CHAR)        , INTENT(  OUT)      :: ErrMsg_C

    ! Local Variables    
    CHARACTER(1024), DIMENSION(55), INTENT(INOUT)      :: data
    TYPE(FileInfoType)            , INTENT(INOUT)      :: InputFileInfo

! Convert python string array to fortran character array
CALL C_F_POINTER(input_file_string_c_ptr, input_file_string_c, 55)
DO I = 1:55:1
    data(I) = input_file_string_c(I)
END DO

! Convert Fortran character array to FileInfoType
CALL InitFileInfo(data, InputFileInfo, ErrStat, ErrMsg)           ! in, out (FileInfoType), out, out

! Parse the input data - ASK ANDY/RAF WHAT THE FIRST ARGUMENT IS!!
CALL InflowWind_ParseInputFileInfo(InputFileData , InputFileInfo, PriPath, TmpErrStat, TmpErrMsg)           ! inout, in, in , out, out

! Pass FileInfoType to InflowWind_Init
! ?

! Convert the outputs of InflowWind_Init from Fortran to C


END SUBROUTINE IFW_INIT_C

SUBROUTINE IFW_CALCOUTPUT_C()

! CODE GOES HERE

END SUBROUTINE IFW_CALCOUTPUT_C

SUBROUTINE IFW_END_C()

! CODE GOES HERE

END SUBROUTINE IFW_END_C

SUBROUTINE InitFileInfo( StringArray, FileInfo, ErrStat, ErrMsg )

      CHARACTER(*), DIMENSION(:), INTENT(IN   ) :: StringArray
      TYPE(FileInfoType),         INTENT(  OUT) :: FileInfo
      INTEGER(IntKi),             INTENT(  OUT) :: ErrStat
      CHARACTER(*),               INTENT(  OUT) :: ErrMsg

      character(len=len(StringArray))  :: TmpStringArray(size(StringArray))
      character(len=len(StringArray))  :: Line
      integer                          :: TmpFileLine(size(StringArray))

      CHARACTER(*), PARAMETER :: RoutineName = 'InitFileInfo'
      INTEGER :: i, NumLines, IC, NumCommChars, LineLen, FirstComm, CommLoc

      ErrStat = ErrID_None
      ErrMsg  = ""
      NumLines = 0      ! Initialize counter for non-comment populated lines
      TmpFileLine = 0   ! Line number that was passed in 
      NumCommChars = LEN_TRIM( CommChars )   ! Number of globally specified CommChars

         ! Find how many non-comment lines we have
      do i=1,size(StringArray)
         Line=StringArray(i)
         LineLen      = LEN_TRIM( Line )
         IF ( ( NumCommChars == 0 ) .OR. ( LineLen == 0 ) ) CYCLE 
 
         FirstComm = MIN( LEN( Line ), LineLen + 1 )
 
         DO IC=1,NumCommChars
            CommLoc = INDEX( Line, CommChars(IC:IC) )
            IF ( CommLoc > 0 )  THEN
               FirstComm = MIN( CommLoc, FirstComm )
            ENDIF
         END DO

            ! Only keep lines with no comments and some sort of length
         if ( LEN_TRIM( Line(:FirstComm-1) ) > 0 ) then
            NumLines=NumLines+1
            TmpStringArray(NumLines) = Line(:FirstComm-1)   ! Store non-comment line
            TmpFileLine(NumLines) = i                        ! Corresponding line number of passed in info
         endif
      enddo

         ! Now save the FileInfo
      FileInfo%NumLines = NumLines        ! only lines that contained anything 
      FileInfo%NumFiles = 1
      ALLOCATE( FileInfo%Lines(FileInfo%NumLines) )
      ALLOCATE( FileInfo%FileLine(FileInfo%NumLines) )
      ALLOCATE( FileInfo%FileIndx(FileInfo%NumLines) )
      ALLOCATE( FileInfo%FileList(FileInfo%NumFiles) )

      DO i = 1, FileInfo%NumLines
         IF ( LEN(TmpStringArray(i)) > LEN(FileInfo%Lines(i)) ) THEN
            CALL SetErrStat( ErrID_Fatal, 'Input string exceeds the bounds of FileInfoType.' , ErrStat, ErrMsg, RoutineName )
            RETURN
         END IF
         FileInfo%Lines(i)    = TmpStringArray(i)
         FileInfo%FileLine(i) = TmpFileLine(i)
      END DO      
      FileInfo%FileIndx = FileInfo%NumFiles
      FileInfo%FileList = (/ "passed file info" /)

END SUBROUTINE InitFileInfo

END MODULE