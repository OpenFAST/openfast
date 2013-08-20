!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2013  National Renewable Energy Laboratory
!
!    This file is part of the NWTC Subroutine Library.
!
!    NWTC Subroutine Library is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
!    published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License along with NWTC Subroutine Library.
!    If not, see <http://www.gnu.org/licenses/>.
!
!**********************************************************************************************************************************
MODULE NWTC_Base


      ! This module stores basic constants and routines that are not system-specific, but may be used in the system-specific routines.
   USE, INTRINSIC               :: ISO_C_Binding
   USE                             Precision

   IMPLICIT  NONE

!=======================================================================

      ! General constants:

   INTEGER, PARAMETER            :: BITS_IN_ADDR  = C_INTPTR_T*8                  ! The number of bits in an address (32-bit or 64-bit).

   INTEGER(IntKi), PARAMETER     :: ChanLen   = 10                                ! The allowable length of channel names (i.e., width of output columns) in the FAST framework


      ! Global Error-level variables:

   INTEGER(IntKi), PARAMETER     :: ErrID_None   = 0
   INTEGER(IntKi), PARAMETER     :: ErrID_Info   = 1
   INTEGER(IntKi), PARAMETER     :: ErrID_Warn   = 2
   INTEGER(IntKi), PARAMETER     :: ErrID_Severe = 3
   INTEGER(IntKi), PARAMETER     :: ErrID_Fatal  = 4

   INTEGER(IntKi)                :: AbortErrLev  = ErrID_Fatal                     ! Note that this is not a PARAMETER


      ! Type definition for dynamically loaded libraries:

   TYPE DLL_Type

      INTEGER(C_INTPTR_T)       :: FileAddr                                        ! The address of file FileName.         (RETURN value from LoadLibrary ) [Windows]
      TYPE(C_PTR)               :: FileAddrX                                       ! The address of file FileName.         (RETURN value dlopen ) [Linux]
      TYPE(C_FUNPTR)            :: ProcAddr                                        ! The address of procedure ProcName.    (RETURN value from GetProcAddress or dlsym)

      CHARACTER(1024)           :: FileName                                        ! The name of the DLL file including the full path to the current working directory.
      CHARACTER(1024)           :: ProcName                                        ! The name of the procedure in the DLL that will be called.

   END TYPE DLL_Type


END MODULE NWTC_Base
