!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2013-2015  National Renewable Energy Laboratory
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
MODULE SysSubs

   ! This module contains routines with system-specific logic and references, including all references to the console unit, CU.
   ! It also contains standard (but not system-specific) routines it uses.
   ! SysGnuWin.f90 is specifically for the GNU Fortran (gfortran) compiler on Windows.
   ! It contains the following routines:
   !     FUNCTION    FileSize( Unit )                                         ! Returns the size (in bytes) of an open file.
   !     FUNCTION    Is_NaN( DblNum )                                         ! Please use IEEE_IS_NAN() instead
   !     FUNCTION    NWTC_ERF( x )
   !     FUNCTION    NWTC_gamma( x )                                          ! Returns the gamma value of its argument.   
   !     SUBROUTINE  FlushOut ( Unit )
   !     SUBROUTINE  GET_CWD( DirName, Status )
   !     SUBROUTINE  MKDIR( new_directory_path )
   !     SUBROUTINE  OpenCon
   !     SUBROUTINE  OpenUnfInpBEFile ( Un, InFile, RecLen, Error )
   !     SUBROUTINE  ProgExit ( StatCode )
   !     SUBROUTINE  Set_IEEE_Constants( NaN_D, Inf_D, NaN, Inf, NaN_S, Inf_S )   
   !     SUBROUTINE  UsrAlarm
   !     SUBROUTINE  WrNR ( Str )
   !     SUBROUTINE  WrOver ( Str )
   !     SUBROUTINE  WriteScr ( Str, Frm )
   !     SUBROUTINE LoadDynamicLib( DLL, ErrStat, ErrMsg )
   !     SUBROUTINE FreeDynamicLib( DLL, ErrStat, ErrMsg )

   USE NWTC_Base

   IMPLICIT NONE

   INTERFACE NWTC_ERF ! Returns the ERF value of its argument
      MODULE PROCEDURE NWTC_ERFR4
      MODULE PROCEDURE NWTC_ERFR8
   END INTERFACE

   INTERFACE NWTC_gamma ! Returns the gamma value of its argument
         ! note: gamma is part of the F08 standard, but may not be implemented everywhere...
      MODULE PROCEDURE NWTC_gammaR4
      MODULE PROCEDURE NWTC_gammaR8
   END INTERFACE

   INTEGER, PARAMETER            :: ConRecL     = 120                               ! The record length for console output.
   INTEGER, PARAMETER            :: CU          = 6                                 ! The I/O unit for the console.
   INTEGER, PARAMETER            :: MaxWrScrLen = 98                                ! The maximum number of characters allowed to be written to a line in WrScr
   LOGICAL, PARAMETER            :: KBInputOK   = .TRUE.                            ! A flag to tell the program that keyboard input is allowed in the environment.
   CHARACTER(*),  PARAMETER      :: NewLine     = ACHAR(10)                         ! The delimiter for New Lines [ Windows is CHAR(13)//CHAR(10); MAC is CHAR(13); Unix is CHAR(10) {CHAR(13)=\r is a line feed, CHAR(10)=\n is a new line}]
   CHARACTER(*),  PARAMETER      :: OS_Desc     = 'GNU Fortran for Windows'         ! Description of the language/OS
   CHARACTER( 1), PARAMETER      :: PathSep     = '\'                               ! The path separator.
   CHARACTER( 1), PARAMETER      :: SwChar      = '/'                               ! The switch character for command-line options.
   CHARACTER(11), PARAMETER      :: UnfForm     = 'UNFORMATTED'                     ! The string to specify unformatted I/O files.

CONTAINS

!=======================================================================
FUNCTION FileSize( Unit )

   ! This function calls the portability routine, FSTAT, to obtain the file size
   ! in bytes corresponding to a file unit number or returns -1 on error.

   INTEGER(B8Ki)                             :: FileSize                      ! The size of the file in bytes to be returned.
   INTEGER, INTENT(IN)                       :: Unit                          ! The I/O unit number of the pre-opened file.
   INTEGER(4)                                :: StatArray(13)                 ! An array returned by FSTAT that includes the file size.
   INTEGER(4)                                :: Status                        ! The status returned by

   Status = FSTAT( INT( Unit, B4Ki ), StatArray )

   IF ( Status /= 0 ) THEN
      FileSize = -1
   ELSE
      FileSize = StatArray(8)
   END IF

   RETURN
END FUNCTION FileSize ! ( Unit )
!=======================================================================
FUNCTION Is_NaN( DblNum )

   ! This routine determines if a REAL(DbKi) variable holds a proper number.
   ! BJJ: this routine is used in CRUNCH.
   ! Note that IsNaN does not exist in earlier versions of gfortran (e.g., 4.2.1),
   ! but does exist in version 4.4. It should be replaced with the standard
   ! IEEE_IS_NAN when gfortran implements it.

   REAL(DbKi), INTENT(IN)       :: DblNum
   LOGICAL                      :: Is_Nan

   Is_NaN = IsNaN( DblNum )

   RETURN
END FUNCTION Is_NaN ! ( DblNum )
!=======================================================================
FUNCTION NWTC_ERFR4( x )

   ! Returns the ERF value of its argument. The result has a value equal  
   ! to the error function: 2/pi * integral_from_0_to_x of e^(-t^2) dt. 

   REAL(SiKi), INTENT(IN)     :: x           ! input 
   REAL(SiKi)                 :: NWTC_ERFR4  ! result
   
   NWTC_ERFR4 = ERF( x )

END FUNCTION NWTC_ERFR4
!=======================================================================
FUNCTION NWTC_ERFR8( x )

   ! Returns the ERF value of its argument. The result has a value equal  
   ! to the error function: 2/pi * integral_from_0_to_x of e^(-t^2) dt. 

   REAL(R8Ki), INTENT(IN)     :: x             ! input 
   REAL(R8Ki)                 :: NWTC_ERFR8    ! result
   
   NWTC_ERFR8 = ERF( x )

END FUNCTION NWTC_ERFR8
!=======================================================================
FUNCTION NWTC_GammaR4( x )

   ! Returns the gamma value of its argument. The result has a value equal  
   ! to a processor-dependent approximation to the gamma function of x. 

   REAL(SiKi), INTENT(IN)     :: x             ! input 
   REAL(SiKi)                 :: NWTC_GammaR4  ! result
   
   NWTC_GammaR4 = gamma( x )

END FUNCTION NWTC_GammaR4
!=======================================================================
FUNCTION NWTC_GammaR8( x )

   ! Returns the gamma value of its argument. The result has a value equal  
   ! to a processor-dependent approximation to the gamma function of x. 

   REAL(R8Ki), INTENT(IN)     :: x             ! input 
   REAL(R8Ki)                 :: NWTC_GammaR8  ! result
   
   NWTC_GammaR8 = gamma( x )

END FUNCTION NWTC_GammaR8
!=======================================================================
SUBROUTINE FlushOut ( Unit )

   ! This subroutine flushes the buffer on the specified Unit.
   ! It is especially useful when printing "running..." type messages.
   
   INTEGER, INTENT(IN)          :: Unit                                         ! The unit number of the file being flushed.

   ! CALL FLUSH ( Unit )

   RETURN
END SUBROUTINE FlushOut ! ( Unit )
!=======================================================================
!bjj note: this subroutine is not tested for this compiler
SUBROUTINE Get_CWD ( DirName, Status )

   ! This routine retrieves the path of the current working directory.

   IMPLICIT NONE

   CHARACTER(*), INTENT(OUT)    :: DirName                                         ! A CHARACTER string containing the path of the current working directory.
   INTEGER,      INTENT(OUT)    :: Status                                          ! Status returned by the call to a portability routine.

   Status = GETCWD ( DirName )

   RETURN
END SUBROUTINE Get_CWD
!=======================================================================
SUBROUTINE MKDIR ( new_directory_path )

   ! This routine creates a given directory if it does not exist.

   implicit none

   character(*), intent(in) :: new_directory_path
   character(1024)          :: make_command
   logical                  :: directory_exists

   ! Check if the directory exists first
   inquire( file=trim(new_directory_path), exist=directory_exists )

   if ( .NOT. directory_exists ) then
      make_command = 'mkdir "'//trim(new_directory_path)//'"'
      call system( make_command )
   endif

END SUBROUTINE MKDIR
!=======================================================================
SUBROUTINE OpenCon

   ! This routine opens the console for standard output.

!bjj: removed for use with CygWin; Because CU = 6 now, this statement is not necessary
!   OPEN ( CU , FILE='/dev/stdout' , STATUS='OLD' )

!   CALL FlushOut ( CU )

   RETURN
END SUBROUTINE OpenCon
!=======================================================================
SUBROUTINE OpenUnfInpBEFile ( Un, InFile, RecLen, Error )

   ! This routine opens a binary input file with data stored in Big Endian format (created on a UNIX machine.)
   ! Data are stored in RecLen-byte records.

   IMPLICIT                        NONE

   INTEGER, INTENT(IN)          :: Un                                           ! Logical unit for the input file.
   CHARACTER(*), INTENT(IN)     :: InFile                                       ! Name of the input file.
   INTEGER, INTENT(IN)          :: RecLen                                       ! Size of records in the input file, in bytes.
   LOGICAL, INTENT(OUT)         :: Error                                        ! Flag to indicate the open failed.
   INTEGER                      :: IOS                                          ! I/O status of OPEN.

   ! Open input file.  Make sure it worked.

   ! The non-standard CONVERT keyword allows us to read UNIX binary files, whose bytes are in reverse order (i.e., stored in BIG ENDIAN format).
   ! NOTE: using RecLen in bytes requires using the /assume:byterecl compiler option!
   OPEN ( Un, FILE=TRIM( InFile ), STATUS='OLD', FORM='UNFORMATTED', ACCESS='DIRECT', RECL=RecLen, IOSTAT=IOS, &
                  ACTION='READ'  )                                              ! Use this for UNIX systems.
   !            ACTION='READ', CONVERT='BIG_ENDIAN' )                         ! Use this for PC systems.

   IF ( IOS /= 0 )  THEN
      Error = .TRUE.
   ELSE
      Error = .FALSE.
   END IF

   RETURN
END SUBROUTINE OpenUnfInpBEFile
!=======================================================================
SUBROUTINE ProgExit ( StatCode )

   ! This routine stops the program.  If the compiler supports the EXIT routine,
   ! pass the program status to it.  Otherwise, do a STOP.

   INTEGER, INTENT(IN)          :: StatCode                                      ! The status code to pass to the OS.

   CALL EXIT ( StatCode )

   ! IF ( StatCode == 0 ) THEN
   !    STOP 0
   ! ELSE
   !    IF ( StatCode < 0 ) THEN
   !       CALL WrScr( 'Invalid STOP code.' )
   !    END IF
   !    STOP 1
   ! END IF

END SUBROUTINE ProgExit ! ( StatCode )
!=======================================================================
SUBROUTINE Set_IEEE_Constants( NaN_D, Inf_D, NaN, Inf, NaN_S, Inf_S )   
      
   ! routine that sets the values of NaN_D, Inf_D, NaN, Inf (IEEE 
   ! values for not-a-number and infinity in sindle and double 
   ! precision) F03 has standard intrinsic routines to do this,  
   ! but Gnu has not yet implemented it. This code will fail if  
   ! the compiler checks for floating-point-error, hence the  
   ! compiler directive FPE_TRAP_ENABLED.

   REAL(DbKi), INTENT(inout)           :: Inf_D          ! IEEE value for NaN (not-a-number) in double precision
   REAL(DbKi), INTENT(inout)           :: NaN_D          ! IEEE value for Inf (infinity) in double precision
   REAL(ReKi), INTENT(inout)           :: Inf            ! IEEE value for NaN (not-a-number)
   REAL(ReKi), INTENT(inout)           :: NaN            ! IEEE value for Inf (infinity)
   REAL(SiKi), INTENT(inout)           :: Inf_S          ! IEEE value for NaN (not-a-number) in single precision
   REAL(SiKi), INTENT(inout)           :: NaN_S          ! IEEE value for Inf (infinity) in single precision

      ! local variables for getting values of NaN and Inf (not necessary when using ieee_arithmetic)
   REAL(DbKi)                          :: Neg_D          ! a negative real(DbKi) number
   REAL(ReKi)                          :: Neg            ! a negative real(ReKi) number
   REAL(SiKi)                          :: Neg_S          ! a negative real(SiKi) number

   
      ! if compiling with floating-point-exception traps, this will not work, so we've added a compiler directive.
      !  note that anything that refers to NaN or Inf will be incorrect in that case.
      
#ifndef FPE_TRAP_ENABLED      
      ! set variables to negative numbers to calculate NaNs (compilers may complain when taking sqrt of negative constants)
   Neg_D = -1.0_DbKi
   Neg   = -1.0_ReKi
   Neg_S = -1.0_SiKi

   NaN_D = SQRT ( Neg_D )
   NaN   = SQRT ( Neg )
   NaN_S = SQRT ( Neg_S )

      ! set variables to zero to calculate Infs (using division by zero)
   Neg_D = 0.0_DbKi
   Neg   = 0.0_ReKi
   Neg_S = 0.0_SiKi
   
   Inf_D = 1.0_DbKi / Neg_D
   Inf   = 1.0_ReKi / Neg
   Inf_S = 1.0_SiKi / Neg_S
#endif 

END SUBROUTINE Set_IEEE_Constants
!=======================================================================
SUBROUTINE UsrAlarm

   ! This routine generates an alarm to warn the user that something went wrong.

   CALL WrNR ( CHAR( 7 ) )

   RETURN
END SUBROUTINE UsrAlarm
!=======================================================================
SUBROUTINE WrNR ( Str )

   ! This routine writes out a string to the screen without following it with a new line.

   CHARACTER(*), INTENT(IN)     :: Str                                          ! The string to write to the screen.

   WRITE (CU,'(1X,A)',ADVANCE='NO')  Str

   RETURN
END SUBROUTINE WrNR ! ( Str )
!=======================================================================
SUBROUTINE WrOver ( Str )

   ! This routine writes out a string that overwrites the previous line.

   CHARACTER(*), INTENT(IN)     :: Str                                          ! The string to write to the screen.
   INTEGER(IntKi)               :: NChars                                       ! Number of characters to write
   CHARACTER(1), PARAMETER      :: CR = ACHAR( 13 )                             ! The carriage return character.
   CHARACTER(25)                :: Fmt = '(2A,   (" "))'                        ! The format specifier for the output.

   NChars = MaxWrScrLen - LEN( Str )

   IF ( NChars > 0 ) THEN
      ! bjj: note that this will produce an error if NChars > 999
      WRITE (Fmt(5:7),'(I3)')  NChars
      WRITE (CU,Fmt,ADVANCE='NO')  CR, Str
   ELSE
      ! bjj: note that this will almost certainly write more than MaxWrScrLen characters on a line
      WRITE (CU,'(2A)',ADVANCE='NO')  CR, Str
   END IF

   RETURN
END SUBROUTINE WrOver ! ( Str )
!=======================================================================
SUBROUTINE WriteScr ( Str, Frm )

   ! This routine writes out a string to the screen.

   IMPLICIT NONE

   CHARACTER(*), INTENT(IN)     :: Str                                         ! The input string to write to the screen.
   CHARACTER(*), INTENT(IN)     :: Frm                                         ! Format specifier for the output.
   INTEGER                      :: ErrStat                                     ! Error status of write operation (so code doesn't crash)

   IF ( LEN_TRIM(Str)  < 1 ) THEN
      WRITE ( CU, '()', IOSTAT=ErrStat )
   ELSE
      WRITE ( CU, Frm, IOSTAT=ErrStat ) TRIM(Str)
   END IF

   IF ( ErrStat /= 0 ) THEN
         print *, ' WriteScr produced an error. Please inform the programmer.'
   ENDIF

END SUBROUTINE WriteScr ! ( Str )
!=======================================================================
SUBROUTINE LoadDynamicLib ( DLL, ErrStat, ErrMsg )

   ! This SUBROUTINE is used to dynamically load a DLL.

   TYPE (DLL_Type),           INTENT(INOUT)  :: DLL         ! The DLL to be loaded.
   INTEGER(IntKi),            INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),              INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   INTERFACE  ! Definitions of Windows API routines

      !...........................
      !bjj: I have been unable to find a solution that works with both IVF and gfortran...
      !bjj: note that "Intel Fortran does not support use of STDCALL with BIND(C) at this time"
      !     See this link: http://software.intel.com/en-us/articles/replacing-intel-fortran-attributes-with-c-interoperability-features
      !bjj: Until this is fixed, Intel uses kernel32.f90 definitions instead of the interface below:
      !...........................

      FUNCTION LoadLibrary(lpFileName) BIND(C,NAME='LoadLibraryA')
         USE, INTRINSIC :: ISO_C_BINDING
         IMPLICIT NONE
         !GCC$ ATTRIBUTES STDCALL :: LoadLibrary
         INTEGER(C_INTPTR_T)        :: LoadLibrary
         CHARACTER(KIND=C_CHAR)     :: lpFileName(*)
      END FUNCTION LoadLibrary

   END INTERFACE

   ErrStat = ErrID_None
   ErrMsg = ''

   ! Load the DLL and get the file address:
   DLL%FileAddr = LoadLibrary( TRIM(DLL%FileName)//C_NULL_CHAR )  !the "C_NULL_CHAR" converts the Fortran string to a C-type string (i.e., adds //CHAR(0) to the end)
   IF ( DLL%FileAddr == INT(0,C_INTPTR_T) ) THEN
      ErrStat = ErrID_Fatal
      WRITE(ErrMsg,'(I2)') BITS_IN_ADDR
      ErrMsg  = 'The dynamic library '//TRIM(DLL%FileName)//' could not be loaded. Check that the file '// &
               'exists in the specified location and that it is compiled for '//TRIM(ErrMsg)//'-bit applications.'
      RETURN
   END IF

   ! Get the procedure address:
   CALL LoadDynamicLibProc ( DLL, ErrStat, ErrMsg )

   RETURN
END SUBROUTINE LoadDynamicLib
!=======================================================================
SUBROUTINE LoadDynamicLibProc ( DLL, ErrStat, ErrMsg )

   ! This SUBROUTINE is used to dynamically load a procedure in a DLL.

   TYPE (DLL_Type),           INTENT(INOUT)  :: DLL         ! The DLL to be loaded.
   INTEGER(IntKi),            INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),              INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   INTEGER(IntKi)                            :: i

   INTERFACE  ! Definitions of Windows API routines

      !...........................
      !bjj: I have been unable to find a solution that works with both IVF and gfortran...
      !bjj: note that "Intel Fortran does not support use of STDCALL with BIND(C) at this time"
      !     See this link: http://software.intel.com/en-us/articles/replacing-intel-fortran-attributes-with-c-interoperability-features
      !bjj: Until this is fixed, Intel uses kernel32.f90 definitions instead of the interface below:
      !...........................

      FUNCTION GetProcAddress(hModule, lpProcName) BIND(C, NAME='GetProcAddress')
         USE, INTRINSIC :: ISO_C_BINDING
         IMPLICIT NONE
         !GCC$ ATTRIBUTES STDCALL :: GetProcAddress
         TYPE(C_FUNPTR)             :: GetProcAddress
         INTEGER(C_INTPTR_T),VALUE  :: hModule
         CHARACTER(KIND=C_CHAR)     :: lpProcName(*)
      END FUNCTION GetProcAddress

   END INTERFACE

   ErrStat = ErrID_None
   ErrMsg = ''

   ! Get the procedure addresses:
   do i=1,NWTC_MAX_DLL_PROC
      if ( len_trim( DLL%ProcName(i) ) > 0 ) then
         DLL%ProcAddr(i) = GetProcAddress( DLL%FileAddr, TRIM(DLL%ProcName(i))//C_NULL_CHAR )  !the "C_NULL_CHAR" converts the Fortran string to a C-type string (i.e., adds //CHAR(0) to the end)
         IF(.NOT. C_ASSOCIATED(DLL%ProcAddr(i))) THEN
            ErrStat = ErrID_Fatal + i - 1
            ErrMsg  = 'The procedure '//TRIM(DLL%ProcName(i))//' in file '//TRIM(DLL%FileName)//' could not be loaded.'
            RETURN
         END IF
      end if
   end do

   RETURN
END SUBROUTINE LoadDynamicLibProc
!=======================================================================
SUBROUTINE FreeDynamicLib ( DLL, ErrStat, ErrMsg )

   ! This SUBROUTINE is used to free a dynamically loaded DLL (loaded in LoadDynamicLib).

   TYPE (DLL_Type),           INTENT(INOUT)  :: DLL         ! The DLL to be freed.
   INTEGER(IntKi),            INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),              INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   ! Local variable:
   INTEGER(C_INT)                            :: Success     ! Whether or not the call to FreeLibrary was successful
   INTEGER(C_INT), PARAMETER                 :: FALSE  = 0

   INTERFACE  ! Definitions of Windows API routines

      FUNCTION FreeLibrary(hLibModule) BIND(C, NAME='FreeLibrary')
         USE, INTRINSIC :: ISO_C_BINDING
         IMPLICIT NONE
         !GCC$ ATTRIBUTES STDCALL :: FreeLibrary
         INTEGER(C_INT)             :: FreeLibrary ! BOOL
         INTEGER(C_INTPTR_T),VALUE  :: hLibModule ! HMODULE hLibModule
      END FUNCTION

   END INTERFACE

   ErrStat = ErrID_None
   ErrMsg = ''
      
   ! Free the DLL:
   IF ( DLL%FileAddr == INT(0,C_INTPTR_T) ) RETURN

   Success = FreeLibrary( DLL%FileAddr ) !If the function succeeds, the return value is nonzero. If the function fails, the return value is zero.

   IF ( Success == FALSE ) THEN !BJJ: note that this isn't the same as the Fortran LOGICAL .FALSE.
      ErrStat = ErrID_Fatal
      ErrMsg  = 'The dynamic library could not be freed.'
      RETURN
   ELSE
      ErrStat = ErrID_None
      ErrMsg = ''
      DLL%FileAddr = INT(0,C_INTPTR_T)
   END IF

   RETURN
END SUBROUTINE FreeDynamicLib
!=======================================================================
END MODULE SysSubs
