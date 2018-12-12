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
!**********************************************************************************************************************************
MODULE SysSubs

   ! This module contains routines with system-specific logic and references, including all references to the console unit, CU.
   ! It also contains standard (but not system-specific) routines it uses.
   ! SysIVF.f90 is specifically for the Intel Visual Fortran for Windows compiler.   
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
   !     SUBROUTINE  Set_IEEE_Constants( NaN_D, Inf_D, NaN, Inf )   
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
      MODULE PROCEDURE NWTC_ERFR16
   END INTERFACE   

   INTERFACE NWTC_gamma ! Returns the gamma value of its argument
         ! note: gamma is part of the F08 standard, but may not be implemented everywhere...
      MODULE PROCEDURE NWTC_gammaR4
      MODULE PROCEDURE NWTC_gammaR8
      MODULE PROCEDURE NWTC_gammaR16
   END INTERFACE

   INTEGER, PARAMETER            :: ConRecL     = 120                               ! The record length for console output.
   INTEGER, PARAMETER            :: CU          = 7                                 ! The I/O unit for the console.  Unit 6 causes ADAMS to crash.
   INTEGER, PARAMETER            :: MaxWrScrLen = 98                                ! The maximum number of characters allowed to be written to a line in WrScr
   LOGICAL, PARAMETER            :: KBInputOK   = .TRUE.                            ! A flag to tell the program that keyboard input is allowed in the environment.
   CHARACTER(*),  PARAMETER      :: NewLine     = ACHAR(10)                         ! The delimiter for New Lines [ Windows is CHAR(13)//CHAR(10); MAC is CHAR(13); Unix is CHAR(10) {CHAR(13)=\r is a line feed, CHAR(10)=\n is a new line}]
   CHARACTER(*),  PARAMETER      :: OS_Desc     = 'Intel Visual Fortran for Windows'! Description of the language/OS
   CHARACTER( 1), PARAMETER      :: PathSep     = '\'                               ! The path separator.
   CHARACTER( 1), PARAMETER      :: SwChar      = '/'                               ! The switch character for command-line options.
   CHARACTER(11), PARAMETER      :: UnfForm     = 'BINARY'                          ! The string to specify unformatted I/O files. (used in OpenUOutFile and OpenUInpFile [see TurbSim's .bin files])

CONTAINS

!=======================================================================
FUNCTION FileSize( Unit )
   ! This function calls the portability routine, FSTAT, to obtain the file size
   ! in bytes corresponding to a file unit number or returns -1 on error.

USE IFPORT

INTEGER(B8Ki)                             :: FileSize                      !< The size of the file in bytes to be returned.
INTEGER, INTENT(IN)                       :: Unit                          !< The I/O unit number of the pre-opened file.
INTEGER                                   :: StatArray(12)                 ! An array returned by FSTAT that includes the file size.
INTEGER                                   :: Status                        ! The status returned by

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

   ! This routine determines if a REAL(DbKi) variable holds a NaN (not-a-number) value.
   ! BJJ: this routine is used in CRUNCH.
   ! It should be replaced with IEEE_IS_NAN, but remains here for
   ! backwards compatibility (some older compilers aren't yet current).

USE, INTRINSIC :: ieee_arithmetic

REAL(DbKi), INTENT(IN)       :: DblNum  !< scalar value in question
LOGICAL                      :: Is_Nan

Is_NaN = IEEE_IS_NAN( DblNum )

RETURN
END FUNCTION Is_NaN ! ( DblNum )
!=======================================================================
!> This function returns the ERF value of its argument. The result has a value equal  
!! to the Gauss error function: 
!! \f{equation}{
!! \mathrm{erf}(x)=\frac{2}{\sqrt{\pi}}\int_0^x e^{-t^2}\,dt
!! \f} \n
!! Use NWTC_ERF (syssubs::nwtc_erf) instead of directly calling a specific routine in the generic interface.
FUNCTION NWTC_ERFR4( x )

   ! Returns the ERF value of its argument. The result has a value equal  
   ! to the error function: 2/pi * integral_from_0_to_x of e^(-t^2) dt. 

   REAL(SiKi), INTENT(IN)     :: x           !< input 
   REAL(SiKi)                 :: NWTC_ERFR4  !< \f$\mathrm{erf}(x)\f$
   
   NWTC_ERFR4 = ERF( x )

END FUNCTION NWTC_ERFR4
!=======================================================================
!> \copydoc syssubs::nwtc_erfr4
FUNCTION NWTC_ERFR8( x )

   REAL(R8Ki), INTENT(IN)     :: x             ! input 
   REAL(R8Ki)                 :: NWTC_ERFR8    ! this function
   
   NWTC_ERFR8 = ERF( x )

END FUNCTION NWTC_ERFR8
!=======================================================================
!> \copydoc syssubs::nwtc_erfr4
FUNCTION NWTC_ERFR16( x )

   ! Returns the ERF value of its argument. The result has a value equal  
   ! to the error function: 2/pi * integral_from_0_to_x of e^(-t^2) dt. 

   REAL(QuKi), INTENT(IN)     :: x             ! input 
   REAL(QuKi)                 :: NWTC_ERFR16   ! result
   
   NWTC_ERFR16 = ERF( x )

END FUNCTION NWTC_ERFR16
!=======================================================================
!> Returns the gamma value of its argument. The result has a value equal  
!! to a processor-dependent approximation to the gamma function of x:
!! \f{equation}{
!! \mathrm{ \Gamma }(x) = \int_0^\infty t^{x-1}e^{-t}\,dt
!! \f} \n
!! Use NWTC_Gamma (syssubs::nwtc_gamma) instead of directly calling a specific routine in the generic interface.
FUNCTION NWTC_GammaR4( x )

! \mathrm{ \Gamma }(x) = \int_0^\infinity t^{x-1}e^{-t}\,dt

   REAL(SiKi), INTENT(IN)     :: x             !< input 
   REAL(SiKi)                 :: NWTC_GammaR4  !< \f$\mathrm{\Gamma}(x)\f$
   
   NWTC_GammaR4 = gamma( x )

END FUNCTION NWTC_GammaR4
!=======================================================================
!> \copydoc syssubs::nwtc_gammar4
FUNCTION NWTC_GammaR8( x )

   REAL(R8Ki), INTENT(IN)     :: x             ! input 
   REAL(R8Ki)                 :: NWTC_GammaR8  ! result
   
   NWTC_GammaR8 = gamma( x )

END FUNCTION NWTC_GammaR8
!=======================================================================
!> \copydoc syssubs::nwtc_gammar4
FUNCTION NWTC_GammaR16( x )

   REAL(QuKi), INTENT(IN)     :: x             ! input 
   REAL(QuKi)                 :: NWTC_GammaR16  ! result
   
   NWTC_GammaR16 = gamma( x )

END FUNCTION NWTC_GammaR16
!=======================================================================
!> This subroutine flushes the buffer on the specified Unit.
!! It is especially useful when printing "running..." type messages.
SUBROUTINE FlushOut ( Unit )

   USE IFPORT

   INTEGER, INTENT(IN)          :: Unit                                         !< The unit number of the file being flushed.

   CALL FLUSH ( INT(Unit, B4Ki) )

   RETURN
END SUBROUTINE FlushOut ! ( Unit )
!=======================================================================
!> This routine retrieves the path of the current working directory.
SUBROUTINE Get_CWD ( DirName, Status )

   USE IFPORT, ONLY: GETCWD

   IMPLICIT NONE

   CHARACTER(*), INTENT(OUT)    :: DirName                                         !< A CHARACTER string containing the path of the current working directory.
   INTEGER,      INTENT(OUT)    :: Status                                          !< Status returned by the call to a portability routine.

   Status = GETCWD ( DirName )

   RETURN
END SUBROUTINE Get_CWD
!=======================================================================
!> This routine creates a given directory if it does not already exist.
SUBROUTINE MKDIR ( new_directory_path )

   implicit none

   character(*), intent(in) :: new_directory_path
   character(1024)          :: make_command
   logical                  :: directory_exists

   ! Check if the directory exists first
   inquire( directory=trim(new_directory_path), exist=directory_exists )

   if ( .NOT. directory_exists ) then
      make_command = 'mkdir "'//trim(new_directory_path)//'"'
      call system( make_command )
   endif

END SUBROUTINE MKDIR
!=======================================================================
!> This routine opens the console for standard output.
SUBROUTINE OpenCon()

   USE IFPORT

   OPEN ( CU , FILE='CON' , STATUS='UNKNOWN' , CARRIAGECONTROL='FORTRAN', RECL=ConRecL )

   CALL FlushOut ( CU )

   RETURN
END SUBROUTINE OpenCon
!=======================================================================
!> This routine opens a binary input file with data stored in Big Endian format (created on a UNIX machine.)
!! Data are stored in RecLen-byte records. (This routine is used for 
SUBROUTINE OpenUnfInpBEFile ( Un, InFile, RecLen, Error )

   IMPLICIT NONE

   INTEGER, INTENT(IN)          :: Un                                           !< Logical unit for the input file.
   CHARACTER(*), INTENT(IN)     :: InFile                                       !< Name of the input file.
   INTEGER, INTENT(IN)          :: RecLen                                       !< Size of records in the input file, in bytes.
   LOGICAL, INTENT(OUT)         :: Error                                        !< Flag to indicate the open failed (true indicates an error occurred).
   INTEGER                      :: IOS                                          ! I/O status of OPEN.

   ! Open input file.  Make sure it worked.

   ! The non-standard CONVERT keyword allows us to read UNIX binary files, whose bytes are in reverse order (i.e., stored in BIG ENDIAN format).
   ! NOTE: using RecLen in bytes requires using the /assume:byterecl compiler option!
   OPEN ( Un, FILE=TRIM( InFile ), STATUS='OLD', FORM='UNFORMATTED', ACCESS='DIRECT', RECL=RecLen, IOSTAT=IOS, &
   !                  ACTION='READ'  )                                              ! Use this for UNIX systems.
                  ACTION='READ', CONVERT='BIG_ENDIAN' )                         ! Use this for PC systems.

   IF ( IOS /= 0 )  THEN
      Error = .TRUE.
   ELSE
      Error = .FALSE.
   END IF

   RETURN
END SUBROUTINE OpenUnfInpBEFile
!=======================================================================
!> This routine stops the program immediately.  If the compiler supports the EXIT routine,
!! pass the program status to it.  Otherwise, do a STOP.
SUBROUTINE ProgExit ( StatCode )

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
!> This routine sets the values of NaN_D, Inf_D, NaN, Inf (IEEE 
!! values for not-a-number and infinity in sindle and double 
!! precision) This uses standard F03 intrinsic routines,  
!! however Gnu has not yet implemented it, so we've placed this
!! routine in the system-specific code.
SUBROUTINE Set_IEEE_Constants( NaN_D, Inf_D, NaN, Inf )   

   USE, INTRINSIC :: ieee_arithmetic  ! use this for compilers that have implemented ieee_arithmetic from F03 standard (otherwise see logic in SysGnu*.f90)

   REAL(DbKi), INTENT(inout)           :: Inf_D          !< IEEE value for NaN (not-a-number) in double precision
   REAL(DbKi), INTENT(inout)           :: NaN_D          !< IEEE value for Inf (infinity) in double precision

   REAL(ReKi), INTENT(inout)           :: Inf            !< IEEE value for NaN (not-a-number)
   REAL(ReKi), INTENT(inout)           :: NaN            !< IEEE value for Inf (infinity)
   
   NaN_D = ieee_value(0.0_DbKi, ieee_quiet_nan)
   Inf_D = ieee_value(0.0_DbKi, ieee_positive_inf)

   NaN   = ieee_value(0.0_ReKi, ieee_quiet_nan)
   Inf   = ieee_value(0.0_ReKi, ieee_positive_inf)   

END SUBROUTINE Set_IEEE_Constants  
!=======================================================================
!> This routine generates an alarm to warn the user that something went wrong.
SUBROUTINE UsrAlarm()

   CALL WrNR ( CHAR( 7 ) )

   RETURN
END SUBROUTINE UsrAlarm
!=======================================================================
!> This routine writes out a string to the screen without following it with a new line.
SUBROUTINE WrNR ( Str )

   CHARACTER(*), INTENT(IN)     :: Str                                          !< The string to write to the screen.

   WRITE (CU,'(1X,A)',ADVANCE='NO')  Str

   RETURN
END SUBROUTINE WrNR ! ( Str )
!=======================================================================
!> This routine writes out a string that overwrites the previous line
SUBROUTINE WrOver ( Str )

   CHARACTER(*), INTENT(IN)     :: Str                                          !< The string to write to the screen.

   WRITE (CU,'(''+'',A)')  Str

   RETURN
END SUBROUTINE WrOver ! ( Str )
!=======================================================================
!> This routine writes out a string to the screen.
SUBROUTINE WriteScr ( Str, Frm )

   IMPLICIT NONE

   CHARACTER(*), INTENT(IN)     :: Str         !< The input string to write to the screen.
   CHARACTER(*), INTENT(IN)     :: Frm         !< Format specifier for the output.
   INTEGER                      :: ErrStat     ! Error status of write operation (so code doesn't crash)

   IF ( LEN_TRIM(Str)  < 1 ) THEN
      WRITE ( CU, '()', IOSTAT=ErrStat )
   ELSE
      WRITE ( CU, Frm, IOSTAT=ErrStat ) TRIM(Str)
   END IF

END SUBROUTINE WriteScr ! ( Str )
!=======================================================================
!> This subroutine is used to dynamically load a DLL, using operating-system API routines to do so.
SUBROUTINE LoadDynamicLib ( DLL, ErrStat, ErrMsg )

   USE IFWINTY,  ONLY : HANDLE
   USE kernel32, ONLY : LoadLibrary

   TYPE (DLL_Type),           INTENT(INOUT)  :: DLL         !< The DLL to be loaded.
   INTEGER(IntKi),            INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),              INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None
   INTEGER(HANDLE)                           :: FileAddr    ! The address of file FileName.         (RETURN value from LoadLibrary in kernel32.f90)

   ErrStat = ErrID_None
   ErrMsg = ''

      ! Load the DLL and get the file address:
   FileAddr = LoadLibrary( TRIM(DLL%FileName)//C_NULL_CHAR )  !the "C_NULL_CHAR" converts the Fortran string to a C-type string (i.e., adds //CHAR(0) to the end)
   DLL%FileAddr = TRANSFER(FileAddr, DLL%FileAddr)             !convert INTEGER(HANDLE) to INTEGER(C_INTPTR_T) [used only for compatibility with gfortran]

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

   ! This subroutine is used to dynamically load a procedure in a DLL, using operating-system API routines to do so.

   USE IFWINTY,  ONLY : LPVOID
   USE kernel32, ONLY : GetProcAddress

   TYPE (DLL_Type),           INTENT(INOUT)  :: DLL         !< The DLL to be loaded.
   INTEGER(IntKi),            INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),              INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None
   INTEGER(LPVOID)                           :: ProcAddr    ! The address of procedure ProcName.    (RETURN value from GetProcAddress in kernel32.f90)
   INTEGER(IntKi)                            :: i

   ErrStat = ErrID_None
   ErrMsg = ''
         
      ! Get the procedure addresses:
   do i=1,NWTC_MAX_DLL_PROC
      if ( len_trim( DLL%ProcName(i) ) > 0 ) then
   
         ProcAddr = GetProcAddress( DLL%FileAddr, TRIM(DLL%ProcName(i))//C_NULL_CHAR )  !the "C_NULL_CHAR" converts the Fortran string to a C-type string (i.e., adds //CHAR(0) to the end)
         DLL%ProcAddr(i) = TRANSFER(ProcAddr, DLL%ProcAddr(i))  !convert INTEGER(LPVOID) to INTEGER(C_FUNPTR) [used only for compatibility with gfortran]

         IF(.NOT. C_ASSOCIATED(DLL%ProcAddr(i))) THEN
            ErrStat = ErrID_Fatal
            ErrMsg  = 'The procedure '//TRIM(DLL%ProcName(i))//' in file '//TRIM(DLL%FileName)//' could not be loaded.'
            RETURN
         END IF
         
      end if
   end do

   RETURN
END SUBROUTINE LoadDynamicLibProc
!=======================================================================
SUBROUTINE FreeDynamicLib ( DLL, ErrStat, ErrMsg )
! This subroutine is used to free a dynamically loaded DLL (loaded in LoadDynamicLib (syssubs::loaddynamiclib)), using operating-system API routines to do so.

   USE IFWINTY,  ONLY : BOOL, HANDLE, FALSE !, LPVOID
   USE kernel32, ONLY : FreeLibrary

   TYPE (DLL_Type),           INTENT(INOUT)  :: DLL         !< The DLL to be freed.
   INTEGER(IntKi),            INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),              INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None
   INTEGER(HANDLE)                           :: FileAddr    ! The address of file FileName.  (RETURN value from LoadLibrary in kernel32.f90)
   INTEGER(BOOL)                             :: Success     ! Whether or not the call to FreeLibrary was successful

   IF ( DLL%FileAddr == INT(0,C_INTPTR_T) ) RETURN
   
   FileAddr = TRANSFER(DLL%FileAddr, FileAddr) !convert INTEGER(C_INTPTR_T) to INTEGER(HANDLE) [used only for compatibility with gfortran]

      ! Free the DLL:
   Success = FreeLibrary( FileAddr ) !If the function succeeds, the return value is nonzero. If the function fails, the return value is zero.

   IF ( Success == FALSE ) THEN !BJJ: note that this is the Windows BOOL type so FALSE isn't the same as the Fortran LOGICAL .FALSE.
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
