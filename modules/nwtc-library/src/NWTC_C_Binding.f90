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
USE NWTC_Base, ONLY: ErrMsgLen

!------------------------------------------------------------------------------------
!  Error handling
!     This must exactly match the value in the python-lib. If ErrMsgLen changes at
!     some point in the nwtc-library, this should be updated, but the logic exists
!     to correctly handle different lengths of the strings
INTEGER(IntKi), PARAMETER   :: ErrMsgLen_C = 8197
INTEGER(IntKi), PARAMETER   :: IntfStrLen  = 1025       ! length of other strings through the C interface

CONTAINS

!> This routine sets the error status in C_CHAR for export to calling code.
!! Make absolutely certain that we do not overrun the end of ErrMsg_C.  That is hard coded to 1025,
!! but ErrMsgLen is set in the nwtc_library, and could change without updates here.  We don't want an
!! inadvertant buffer overrun -- that can lead to bad things.
SUBROUTINE SetErrStat_F2C(ErrStat, ErrMsg, ErrStat_C, ErrMsg_C)
    INTEGER,                INTENT(IN   )  :: ErrStat                 !< aggregated error message (fortran type)
    CHARACTER(ErrMsgLen),   INTENT(IN   )  :: ErrMsg                  !< aggregated error message (fortran type)
    INTEGER(C_INT),         INTENT(  OUT)  :: ErrStat_C
    CHARACTER(KIND=C_CHAR), INTENT(  OUT)  :: ErrMsg_C(ErrMsgLen_C)

    ErrStat_C = ErrStat     ! We will send back the same error status that is used in OpenFAST
    if (ErrMsgLen > ErrMsgLen_C-1) then   ! If ErrMsgLen is > the space in ErrMsg_C, do not copy everything over
        ErrMsg_C = TRANSFER( TRIM(ErrMsg(1:ErrMsgLen_C-1))//C_NULL_CHAR, ErrMsg_C )
    else
        ErrMsg_C = TRANSFER( TRIM(ErrMsg)//C_NULL_CHAR, ErrMsg_C )
    endif
END SUBROUTINE SetErrStat_F2C

!> This subroutine incorporates the local error status and error messages into the global error
!! status and message. It expects both local and global error messages to be null-terminated
!! C strings.
SUBROUTINE SetErrStat_C(ErrStatLocal, ErrMessLocal, ErrStatGlobal, ErrMessGlobal, RoutineName)

    INTEGER(C_INT),                          INTENT(IN   ) :: ErrStatLocal                  ! Error status of the operation
    CHARACTER(KIND=C_CHAR, LEN=ErrMsgLen_C), INTENT(IN   ) :: ErrMessLocal                  ! Error message if ErrStat /= ErrID_None
    INTEGER(C_INT),                          INTENT(INOUT) :: ErrStatGlobal                 ! Error status of the operation
    CHARACTER(KIND=C_CHAR),                  INTENT(INOUT) :: ErrMessGlobal(ErrMsgLen_C)    ! Error message if ErrStat /= ErrID_None
    CHARACTER(*),                            INTENT(IN   ) :: RoutineName                   ! Name of the routine error occurred in

    IF ( ErrStatLocal == ErrID_None ) RETURN

    IF (ErrStatGlobal /= ErrID_None) THEN
        ! print *, "in if", ErrStatGlobal, ErrID_None
        ! ErrMessGlobal = TRIM(ErrMessGlobal)//new_line('a')
        ! print *, "ErrMessGlobal", ErrMessGlobal
    ENDIF
    ! TODO: Does the line below need to remove the null char from ErrMessGlobal prior to appending to it?
    ErrMessGlobal = TRANSFER( ErrMessGlobal//TRIM(RoutineName)//':'//TRIM(ErrMessLocal)//C_NULL_CHAR, ErrMessGlobal )
    ! ErrMessGlobal = TRIM(ErrMessGlobal)//TRIM(RoutineName)//':'//TRIM(ErrMessLocal)
    ErrStatGlobal = MAX(ErrStatGlobal, ErrStatLocal)

END SUBROUTINE

SUBROUTINE StringConvert_F2C(F_String, C_String, Length)
    !> This was taken from https://github.com/vmagnin/gtk-fortran/blob/gtk4/src/gtk-sup.f90#L640.
    !! The original implementation has textptr as an allocatable array, but in OpenFAST
    !! the C string arrays will already be created in the calling code.

    ! Convert a Fortran string into a null-terminated C-string
    !
    ! F_STRING |  F_String |  required |  The Fortran string to convert
    ! TEXTPR |  string |  required |  A C type string, (allocatable).
    ! LENGTH |  c_int |  optional |  The Length of the generated C string.
    !-

    character(len=*), intent(in) :: F_String
    character(kind=c_char), dimension(:), intent(out) :: C_String   ! , allocatable
    integer(c_int), intent(out), optional :: Length

    integer :: lcstr, j
    logical :: add_null

    lcstr = len_trim(F_String)
    if (lcstr == 0) then
        lcstr = lcstr + 1
        add_null = .true.
    else if (F_String(lcstr:lcstr) /= c_null_char) then
        lcstr = lcstr + 1
        add_null = .true.
    else
        add_null = .false.
    end if

    ! allocate(C_String(lcstr))
    if (present(Length)) Length = lcstr

    do j = 1, len_trim(F_String)
        C_String(j) = F_String(j:j)
    end do

    if (add_null) C_String(lcstr) = c_null_char

END SUBROUTINE

FUNCTION RemoveCStringNullChar(String_C, StringLength_C)
    INTEGER(C_INT), INTENT(IN)                              :: StringLength_C
    CHARACTER(KIND=C_CHAR, LEN=StringLength_C), INTENT(IN)  :: String_C
    CHARACTER(LEN=StringLength_C)                           :: RemoveCStringNullChar

    RemoveCStringNullChar = String_C

    ! if this has a c null character at the end, remove it
    i = INDEX(RemoveCStringNullChar, C_NULL_CHAR) - 1
    if ( i > 0 ) RemoveCStringNullChar = RemoveCStringNullChar(1:I)
    RETURN

END FUNCTION

FUNCTION FileNameFromCString(FileString_C, FileStringLength_C)
    !> This function takes a string from the C interface and returns a file name
    !> that is compatible with the Fortran interface.  The C string may have
    !> trailing null characters that need to be removed.
    !> By convention, the filename must have fewer characters than IntfStrLen.
    INTEGER(C_INT), INTENT(IN)                                  :: FileStringLength_C   !< length of input string from C interface
    CHARACTER(KIND=C_CHAR, LEN=FileStringLength_C), INTENT(IN)  :: FileString_C         !< input string from C interface
    CHARACTER(LEN=IntfStrLen)                                   :: FileNameFromCString   !< output file name (fortran type)

    INTEGER :: i

    FileNameFromCString = ''
    i = MIN(IntfStrLen, FileStringLength_C)
    FileNameFromCString(1:i) = FileString_C(1:i)
    
    FileNameFromCString = RemoveCStringNullChar(FileNameFromCString, IntfStrLen)

    RETURN

END FUNCTION

END MODULE