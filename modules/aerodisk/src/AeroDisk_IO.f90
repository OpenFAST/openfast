!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2020  National Renewable Energy Laboratory
!
!    This file is part of AeroDisk
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
MODULE AeroDisk_IO

   USE AeroDisk_Types
   USE AeroDisk_Output_Params
   USE NWTC_Library

   implicit none


contains

SUBROUTINE ADsk_ValidateInput( InitInp, InputFileData, ErrStat, ErrMsg )
   type(ADsk_InitInputType),  intent(in   )  :: InitInp              !< Input data for initialization
   type(ADsk_InputFile),      intent(inout)  :: InputFileData        !< The data for initialization
   integer(IntKi),            intent(  out)  :: ErrStat              !< Error status  from this subroutine
   character(*),              intent(  out)  :: ErrMsg               !< Error message from this subroutine
   integer(IntKi)                            :: ErrStat2             !< Temporary error status  for subroutine and function calls
   character(ErrMsgLen)                      :: ErrMsg2              !< Temporary error message for subroutine and function calls
   integer(IntKi)                            :: I                    !< Generic counter
   character(*),              parameter      :: RoutineName="ADsk_ValidateInput"
   integer(IntKi)                            :: IOS                  !< Temporary error status

      ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ""

END SUBROUTINE ADsk_ValidateInput


!FIXME: add SetOutParam here


END MODULE AeroDisk_IO

