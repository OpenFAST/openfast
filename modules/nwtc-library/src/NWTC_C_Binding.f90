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

!> This module contains I/O-related variables and routines with non-system-specific logic.
MODULE NWTC_C_Binding

USE iso_c_binding
USE Precision
USE NWTC_Base, ONLY: ErrMsgLen

!------------------------------------------------------------------------------------
!  Error handling
!     This must exactly match the value in the python-lib. If ErrMsgLen changes at
!     some point in the nwtc-library, this should be updated, but the logic exists
!     to correctly handle different lengths of the strings
INTEGER(IntKi), PARAMETER            :: ErrMsgLen_C = 1025
INTEGER(IntKi), PARAMETER            :: IntfStrLen  = 1025       ! length of other strings through the C interface

CONTAINS

!> This routine sets the error status in C_CHAR for export to calling code.
!! Make absolutely certain that we do not overrun the end of ErrMsg_C.  That is hard coded to 1025,
!! but ErrMsgLen is set in the nwtc_library, and could change without updates here.  We don't want an
!! inadvertant buffer overrun -- that can lead to bad things.
subroutine SetErr(ErrStat, ErrMsg, ErrStat_C, ErrMsg_C)
    integer,                intent(in   )  :: ErrStat                 !< aggregated error message (fortran type)
    character(ErrMsgLen),   intent(in   )  :: ErrMsg                  !< aggregated error message (fortran type)
    integer(c_int),         intent(  out)  :: ErrStat_C
    character(kind=c_char), intent(  out)  :: ErrMsg_C(ErrMsgLen_C)
    ! integer                                :: i
    ErrStat_C = ErrStat     ! We will send back the same error status that is used in OpenFAST
    if (ErrMsgLen > ErrMsgLen_C-1) then   ! If ErrMsgLen is > the space in ErrMsg_C, do not copy everything over
        ErrMsg_C = TRANSFER( trim(ErrMsg(1:ErrMsgLen_C-1))//C_NULL_CHAR, ErrMsg_C )
    else
        ErrMsg_C = TRANSFER( trim(ErrMsg)//C_NULL_CHAR, ErrMsg_C )
    endif
end subroutine SetErr

function RemoveCStringNullChar(String_c, StringLength_c)
    integer(c_int), intent(in) :: StringLength_c
    character(kind=c_char, len=StringLength_c), intent(in)  :: String_c
    character(len=StringLength_c) :: RemoveCStringNullChar

    RemoveCStringNullChar = String_c

    ! if this has a c null character at the end, remove it
    i = INDEX(RemoveCStringNullChar, C_NULL_CHAR) - 1
    if ( i > 0 ) RemoveCStringNullChar = RemoveCStringNullChar(1:I)
    return

end function

function FileNameFromCString(FileString_c, FileStringLength_c)
    !> This function takes a string from the C interface and returns a file name
    !> that is compatible with the Fortran interface.  The C string may have
    !> trailing null characters that need to be removed.
    !> By convention, the filename must have fewer characters than IntfStrLen.
    integer(c_int), intent(in)                                  :: FileStringLength_c   !< length of input string from C interface
    character(kind=c_char, len=FileStringLength_c), intent(in)  :: FileString_c         !< input string from C interface
    character(len=IntfStrLen)                                   :: FileNameFromCString   !< output file name (fortran type)

    integer :: i

    FileNameFromCString = ''
    i = min(IntfStrLen, FileStringLength_C)
    FileNameFromCString(1:i) = FileString_c(1:i)
    
    FileNameFromCString = RemoveCStringNullChar(FileNameFromCString, IntfStrLen)

    return

end function

END MODULE