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
USE InflowWind_short
USE InflowWind_Types               
USE NWTC_Library

IMPLICIT NONE

PUBLIC :: IFW_INIT_C
PUBLIC :: IFW_CALCOUTPUT_C
PUBLIC :: IFW_END_C

CONTAINS

SUBROUTINE IFW_INIT_C(input_file_string_c_ptr,ErrStat_C,ErrMsg_C) BIND (C, NAME='IFW_INIT_C')

    TYPE(C_PTR)                   , INTENT(IN   )      :: input_file_string_c_ptr(*)
!    CHARACTER(KIND=C_CHAR),POINTER, INTENT(INOUT)      :: input_file_string_c(:)
    INTEGER(C_INT)                , INTENT(  OUT)      :: ErrStat_C
    CHARACTER(KIND=C_CHAR)        , INTENT(  OUT)      :: ErrMsg_C

    ! Local Variables    
    CHARACTER(1024), DIMENSION(55)   :: data
    INTEGER                          :: ErrStat
    CHARACTER(ErrMsgLen)             :: ErrMsg
    TYPE(InflowWind_InitInputType)   :: InitInp           

! Convert python string array to fortran character array
! CALL C_F_POINTER(input_file_string_c_ptr, input_file_string_c, 55)
! DO I = 1:55:1
!     data(I) = input_file_string_c(I)
! END DO

! Convert Fortran character array to FileInfoType within InflowWind_InitInputType
CALL InitFileInfo(data, InitInp%PassedFileData, ErrStat, ErrMsg)           ! in, out (FileInfoType), out, out
PRINT*, ErrStat
PRINT*, "Done calling InitFileInfo ....."

InitInp%UseInputFile = .false.
InitInp%InputFileName = "passed_ifw_file" ! dummy
InitInp%RootName = "ifwRoot" ! used for making echo files
InitInp%NumWindPoints = 100 ! arbirtrary number - how many data points to interrogate for each time step, sets array sizes

! Pass FileInfoType into InflowWind_Init
CALL InflowWind_Init(InitInp, ErrStat, ErrMsg)
PRINT*, ErrStat
PRINT*, "Done calling InflowWind_Init ....."

! Convert the outputs of InflowWind_Init from Fortran to C
PRINT*, ErrMsg
   if (ErrStat /= 0) then
      ErrStat_c = ErrID_Fatal
   else
      ErrStat_c = ErrID_None
   end if
ErrMsg_C = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_c )
PRINT*, "DONE WITH IFW_INIT_C! RETURNING"

END SUBROUTINE IFW_INIT_C

SUBROUTINE IFW_CALCOUTPUT_C()

! CODE GOES HERE

END SUBROUTINE IFW_CALCOUTPUT_C

SUBROUTINE IFW_END_C()

! CODE GOES HERE

END SUBROUTINE IFW_END_C

END MODULE