!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2013-2016  National Renewable Energy Laboratory
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
!
!**********************************************************************************************************************************
!> This module stores basic constants and routines that are not system-specific, but may be used in the system-specific routines.
MODULE NWTC_Base


   USE, INTRINSIC               :: ISO_C_Binding
   USE                             Precision

   IMPLICIT  NONE

!=======================================================================

   !logical :: debug_print = .false.
   
      ! General constants:

   INTEGER, PARAMETER            :: BITS_IN_ADDR  = C_INTPTR_T*8                  !< The number of bits in an address (32-bit or 64-bit).
   INTEGER, PARAMETER            :: ErrMsgLen = 1024                              !< The maximum number of characters in an error message in the FAST framework
   
   INTEGER(IntKi), PARAMETER     :: ChanLen   = 20                                !< The maximum allowable length of channel names (i.e., width of output columns) in the FAST framework
   INTEGER(IntKi), PARAMETER     :: OutStrLenM1 = ChanLen - 1                     !< The maximum allowable length of channel names without optional "-" or "M" at the beginning to indicate the negative of the channel
   
   INTEGER(IntKi), PARAMETER     :: MinChanLen = 10                               !< The min allowable length of channel names (i.e., width of output columns), used because some modules (like Bladed DLL outputs) have excessively long names
   INTEGER(IntKi), PARAMETER     :: LinChanLen = 200                              !< The allowable length of row/column names in linearization files
   INTEGER(IntKi), PARAMETER     :: MaxFileInfoLineLen = 1024                     !< The allowable length of an input line stored in FileInfoType%Lines

   INTEGER(IntKi), PARAMETER     :: NWTC_Verbose = 10                             !< The maximum level of verbosity
   INTEGER(IntKi), PARAMETER     :: NWTC_VerboseLevel = 5                         !< a number in [0, NWTC_Verbose]: 0 = no output; NWTC_Verbose=verbose; 

      ! Global Error-level variables:

   INTEGER(IntKi), PARAMETER     :: ErrID_None   = 0                              !< ErrStat parameter indicating "no error"
   INTEGER(IntKi), PARAMETER     :: ErrID_Info   = 1                              !< ErrStat parameter indicating "informational message"
   INTEGER(IntKi), PARAMETER     :: ErrID_Warn   = 2                              !< ErrStat parameter indicating "warning"
   INTEGER(IntKi), PARAMETER     :: ErrID_Severe = 3                              !< ErrStat parameter indicating "severe error"; 
   INTEGER(IntKi), PARAMETER     :: ErrID_Fatal  = 4                              !< ErrStat parameter indicating "fatal error"; simulation should end

   INTEGER(IntKi)                :: AbortErrLev  = ErrID_Fatal                    !< ErrStat that indicates the error level when program should end; ErrID_Fatal by default. Note that this is not a PARAMETER

   
   INTEGER(IntKi), PARAMETER     :: NWTC_MAX_DLL_PROC  = 5                        !< maximum number of procedures that can be dynamically loaded from a DLL (see DLL_Type nwtc_base::dll_type)
   
#ifdef FLANG_COMPILER
   TYPE(C_FUNPTR), PARAMETER     :: NULL_PROC_ADDR(NWTC_MAX_DLL_PROC) = C_NULL_FUNPTR  !< this is a hack so the Flang compiler will initialize ProcAddr to C_NULL_FUNPTR in DLL_Type (remove if no longer needed)
#endif

      !> Type definition for dynamically loaded libraries:
      !! Note that changes here may need to be reflected in DLLTypePack() (nwtc_io::dlltypepack) DLLTypeUnPack() (nwtc_io::dlltypeunpack), 
      !! and the FAST Registry executable.

   TYPE DLL_Type 

      INTEGER(C_INTPTR_T)       :: FileAddr  = INT(0,C_INTPTR_T)                   !< The address of file FileName.         (RETURN value from LoadLibrary ) [Windows]
      TYPE(C_PTR)               :: FileAddrX = C_NULL_PTR                          !< The address of file FileName.         (RETURN value from dlopen ) [Linux]
#ifdef FLANG_COMPILER
      TYPE(C_FUNPTR)            :: ProcAddr(NWTC_MAX_DLL_PROC) = NULL_PROC_ADDR    !< The address of procedure ProcName.    (RETURN value from GetProcAddress or dlsym) [initialized to Null for pack/unpack]
#else
      TYPE(C_FUNPTR)            :: ProcAddr(NWTC_MAX_DLL_PROC) = C_NULL_FUNPTR     !< The address of procedure ProcName.    (RETURN value from GetProcAddress or dlsym) [initialized to Null for pack/unpack]
#endif

      CHARACTER(1024)           :: FileName                                        !< The name of the DLL file including the full path to the current working directory.
      CHARACTER(1024)           :: ProcName(NWTC_MAX_DLL_PROC)  = ""               !< The name of the procedure in the DLL that will be called.

   END TYPE DLL_Type


END MODULE NWTC_Base
