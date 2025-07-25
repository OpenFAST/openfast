!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2025  National Renewable Energy Laboratory
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
!**********************************************************************************************************************************

MODULE NWTC_C_Binding

USE ISO_C_BINDING
USE Precision
USE NWTC_Base, ONLY: ErrMsgLen, ErrID_None, ErrID_Info, ErrID_Warn, ErrID_Severe, ErrID_Fatal, SetErrStat

IMPLICIT NONE

!------------------------------------------------------------------------------------
!  Error handling
!     This must exactly match the value in the python-lib. If ErrMsgLen changes at
!     some point in the nwtc-library, this should be updated, but the logic exists
!     to correctly handle different lengths of the strings
INTEGER(IntKi), PARAMETER   :: ErrMsgLen_C = ErrMsgLen + 1  ! Currently, this is 8197
INTEGER(IntKi), PARAMETER   :: IntfStrLen  = 1025           ! length of other strings through the C interface such as file paths

CONTAINS

!> This routine sets the error status in C_CHAR for export to calling code.
SUBROUTINE SetErrStat_F2C(ErrStat_F, ErrMsg_F, ErrStat_C, ErrMsg_C)
    INTEGER,                INTENT(IN   )  :: ErrStat_F             !< aggregated error status  (fortran type)
    CHARACTER(ErrMsgLen),   INTENT(IN   )  :: ErrMsg_F              !< aggregated error message (fortran type)
    INTEGER(C_INT),         INTENT(  OUT)  :: ErrStat_C
    CHARACTER(KIND=C_CHAR), INTENT(  OUT)  :: ErrMsg_C(ErrMsgLen_C)

    ErrStat_C = ErrStat_F     ! We will send back the same error status that is used in OpenFAST
    if (ErrMsgLen > ErrMsgLen_C-1) then   ! If ErrMsgLen is > the space in ErrMsg_C, do not copy everything over
        ErrMsg_C = TRANSFER( TRIM(ErrMsg_F(1:ErrMsgLen_C-1))//C_NULL_CHAR, ErrMsg_C )
    else
        ErrMsg_C = TRANSFER( TRIM(ErrMsg_F)//C_NULL_CHAR, ErrMsg_C )
    endif
END SUBROUTINE SetErrStat_F2C

!> This subroutine incorporates the local error status and error messages into the global error
!! status and message. It expects both local and global error messages to be null-terminated
!! C strings.
!! The routine name must be a Fortran string with assumed length.
SUBROUTINE SetErrStat_C(ErrStatLocal_C, ErrMessLocal_C, ErrStatGlobal_C, ErrMessGlobal_C, RoutineName_F)

    INTEGER(C_INT),                          INTENT(IN   ) :: ErrStatLocal_C                  ! Error status of the operation
    CHARACTER(KIND=C_CHAR),                  INTENT(IN   ) :: ErrMessLocal_C(ErrMsgLen_C)     ! Error message if ErrStat /= ErrID_None
    INTEGER(C_INT),                          INTENT(INOUT) :: ErrStatGlobal_C                 ! Error status of the operation
    CHARACTER(KIND=C_CHAR),                  INTENT(INOUT) :: ErrMessGlobal_C(ErrMsgLen_C)    ! Error message if ErrStat /= ErrID_None
    CHARACTER(*),                            INTENT(IN   ) :: RoutineName_F                   ! Name of the routine error occurred in

    INTEGER(IntKi)           :: ErrStatLocal_F, ErrStatGlobal_F
    character(len=ErrMsgLen) :: ErrMessLocal_F, ErrMessGlobal_F

    ! Convert C strings to Fortran for easier processing within this subroutine
    CALL StringConvert_C2F(ErrMessLocal_C, ErrMessLocal_F)
    CALL StringConvert_C2F(ErrMessGlobal_C, ErrMessGlobal_F)
    ErrStatLocal_F = INT(ErrStatLocal_C, IntKi)
    ErrStatGlobal_F = INT(ErrStatGlobal_C, IntKi)

    ! Return with no-op if the local error status is None
    IF ( ErrStatLocal_F == ErrID_None ) RETURN

    ! Call the standard NWTC Library error handling routine
    CALL SetErrStat(ErrStatLocal_F, ErrMessLocal_F, ErrStatGlobal_F, ErrMessGlobal_F, RoutineName_F)

    ! Convert outputs back to C types
    ErrStatGlobal_C = INT(ErrStatGlobal_F, C_INT)
    CALL StringConvert_F2C(ErrMessGlobal_F, ErrMessGlobal_C)

END SUBROUTINE

SUBROUTINE StringConvert_F2C(String_F, String_C)

    ! Convert a Fortran string into a null-terminated C-string
    ! NOTE this does not check whether String_C is long enough to hold the Fortran string
    ! If not, it will simply overrun String_C, so the calling code must be sure.

    ! This was inspired by https://github.com/vmagnin/gtk-fortran/blob/gtk4/src/gtk-sup.f90.

    CHARACTER(LEN=*), INTENT(IN) :: String_F
    CHARACTER(KIND=C_CHAR), INTENT(OUT) :: String_C(:)

    INTEGER :: i
    INTEGER :: STRING_LEN
    LOGICAL :: ADD_NULL

    ! Determine if the null terminator needs to be added
    STRING_LEN = LEN_TRIM(String_F)

    ! If the string is empty, add a null terminator
    IF (STRING_LEN == 0) THEN
        STRING_LEN = STRING_LEN + 1
        ADD_NULL = .true.

    ! Otherwise, if the last character is not a null terminator, then add it
    ELSE IF (String_F(STRING_LEN:STRING_LEN) /= C_NULL_CHAR) THEN
        STRING_LEN = STRING_LEN + 1
        ADD_NULL = .true.

    ! Otherwise, do not add a null terminator
    ELSE
        ADD_NULL = .false.

    END IF

    DO i = 1, STRING_LEN - 1
        String_C(i) = String_F(i:i)
    END DO

    IF (ADD_NULL) String_C(STRING_LEN) = C_NULL_CHAR

END SUBROUTINE

SUBROUTINE StringConvert_C2F(String_C, String_F)
    ! Convert a null-terminated C-string to a Fortran string
    ! If the C string is longer than the Fortran string, it will be truncated.

    ! This was inspired by https://github.com/vmagnin/gtk-fortran/blob/gtk4/src/gtk-sup.f90.

    CHARACTER(KIND=C_CHAR), INTENT(IN) :: String_C(:)
    CHARACTER(LEN=*), INTENT(OUT) :: String_F

    INTEGER :: i

    DO i = 1, SIZE(String_C)
        IF (String_C(i) == C_NULL_CHAR) EXIT
        IF (i > LEN(String_F)) RETURN
        String_F(i:i) = String_C(i)
    END DO

    String_F(i:) = ''

END SUBROUTINE

FUNCTION RemoveCStringNullChar(StringLength_C, String_C)
    INTEGER(C_INT), INTENT(IN)                              :: StringLength_C
    CHARACTER(KIND=C_CHAR), INTENT(IN)                      :: String_C(StringLength_C)
    CHARACTER(LEN=StringLength_C)                           :: RemoveCStringNullChar

    integer :: i

    CALL StringConvert_C2F(String_C, RemoveCStringNullChar)

    ! if this has a c null character at the end, remove it
    i = INDEX(RemoveCStringNullChar, C_NULL_CHAR) - 1
    IF ( i > 0 ) RemoveCStringNullChar = RemoveCStringNullChar(1:i)
    RETURN

END FUNCTION

FUNCTION FileNameFromCString(String_C, StringLength_C)
    !> This function takes a string from the C interface and returns a file name
    !> that is compatible with the Fortran interface.  The C string may have
    !> trailing null characters that need to be removed.
    !> By convention, the filename must have fewer characters than IntfStrLen.
    INTEGER(C_INT), INTENT(IN)                              :: StringLength_C       !< length of input string from C interface
    CHARACTER(KIND=C_CHAR, LEN=StringLength_C), INTENT(IN)  :: String_C             !< input string from C interface
    CHARACTER(LEN=IntfStrLen)                               :: FileNameFromCString  !< output file name (fortran type)

    INTEGER :: i

    FileNameFromCString = ''
    i = MIN(IntfStrLen, StringLength_C)
    FileNameFromCString(1:i) = String_C(1:i)
    
    FileNameFromCString = RemoveCStringNullChar(IntfStrLen, FileNameFromCString)

    RETURN

END FUNCTION

END MODULE
