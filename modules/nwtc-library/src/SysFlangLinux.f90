!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2013  National Renewable Energy Laboratory
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
    ! SysGnuLinux.f90 is specifically for the GNU Fortran (gfortran) compiler on Linux and macOS.
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
       MODULE PROCEDURE NWTC_gammaR4
       MODULE PROCEDURE NWTC_gammaR8
    END INTERFACE
 
    INTEGER, PARAMETER            :: ConRecL     = 120                               ! The record length for console output.
    INTEGER, PARAMETER            :: CU          = 6                                 ! The I/O unit for the console.  Unit 6 causes ADAMS to crash.
    INTEGER, PARAMETER            :: MaxWrScrLen = 98                                ! The maximum number of characters allowed to be written to a line in WrScr
    LOGICAL, PARAMETER            :: KBInputOK   = .TRUE.                            ! A flag to tell the program that keyboard input is allowed in the environment.
    CHARACTER(*),  PARAMETER      :: NewLine     = ACHAR(10)                         ! The delimiter for New Lines [ Windows is CHAR(13)//CHAR(10); MAC is CHAR(13); Unix is CHAR(10) {CHAR(13)=\r is a line feed, CHAR(10)=\n is a new line}]
    CHARACTER(*),  PARAMETER      :: OS_Desc     = 'GNU Fortran for Linux'           ! Description of the language/OS
    CHARACTER( 1), PARAMETER      :: PathSep     = '/'                               ! The path separator.
    CHARACTER( 1), PARAMETER      :: SwChar      = '-'                               ! The switch character for command-line options.
    CHARACTER(11), PARAMETER      :: UnfForm     = 'UNFORMATTED'                     ! The string to specify unformatted I/O files.

    #ifdef USE_DLL_INTERFACE

    INTEGER(C_INT), PARAMETER :: RTLD_LAZY=1            ! "Perform lazy binding. Only resolve symbols as the code that references them is executed. If the symbol is never referenced, then it is never resolved. (Lazy binding is only performed for function references; references to variables are always immediately bound when the library is loaded.) "
    INTEGER(C_INT), PARAMETER :: RTLD_NOW=2             ! "If this value is specified, or the environment variable LD_BIND_NOW is set to a nonempty string, all undefined symbols in the library are resolved before dlopen() returns. If this cannot be done, an error is returned."
    INTEGER(C_INT), PARAMETER :: RTLD_GLOBAL=256        ! "The symbols defined by this library will be made available for symbol resolution of subsequently loaded libraries"
    INTEGER(C_INT), PARAMETER :: RTLD_LOCAL=0           ! "This is the converse of RTLD_GLOBAL, and the default if neither flag is specified. Symbols defined in this library are not made available to resolve references in subsequently loaded libraries."
 
    INTERFACE !linux API routines
       !bjj see http://linux.die.net/man/3/dlopen
       !    and https://developer.apple.com/library/mac/documentation/Darwin/Reference/ManPages/man3/dlopen.3.html
 
       FUNCTION dlOpen(filename,mode) BIND(C,NAME="dlopen")
       ! void *dlopen(const char *filename, int mode);
          USE ISO_C_BINDING
          IMPLICIT NONE
          TYPE(C_PTR)                   :: dlOpen
          CHARACTER(C_CHAR), INTENT(IN) :: filename(*)
          INTEGER(C_INT), VALUE         :: mode
       END FUNCTION

       FUNCTION dlSym(handle,name) BIND(C,NAME="dlsym")
        ! void *dlsym(void *handle, const char *name);
           USE ISO_C_BINDING
           IMPLICIT NONE
           TYPE(C_FUNPTR)                :: dlSym ! A function pointer
           TYPE(C_PTR), VALUE            :: handle
           CHARACTER(C_CHAR), INTENT(IN) :: name(*)
        END FUNCTION

        FUNCTION dlClose(handle) BIND(C,NAME="dlclose")
        ! int dlclose(void *handle);
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT)       :: dlClose
            TYPE(C_PTR), VALUE   :: handle
        END FUNCTION
 
    END INTERFACE

    #endif
 
 CONTAINS
 
 
 FUNCTION Is_NaN( DblNum )
 
    ! This routine determines if a REAL(DbKi) variable holds a proper number.
    ! BJJ: this routine is used in CRUNCH.
    ! Note that IsNaN does not exist in earlier versions of gfortran (e.g., 4.2.1),
    ! but does exist in version 4.4. It should be replaced with the standard
    ! IEEE_IS_NAN when gfortran implements it.
 
    REAL(DbKi), INTENT(IN)       :: DblNum
    LOGICAL                      :: Is_Nan
 
    Is_NaN = IsNaN( DblNum )
 
 END FUNCTION
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
 !bjj note: this subroutine is not tested for this compiler
 SUBROUTINE Get_CWD ( DirName, Status )
    CHARACTER(1024)              :: pwd
    INTEGER                      :: length
    CHARACTER(*), INTENT(OUT)    :: DirName                                         ! A CHARACTER string containing the path of the current working directory.
    INTEGER,      INTENT(OUT)    :: Status                                          ! Status returned by the call to a portability routine.
    call get_environment_variable('PWD', DirName, length, Status)
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
       make_command = 'mkdir -p '//trim(new_directory_path)
       call system( make_command )
    endif
 
 END SUBROUTINE MKDIR
 !=======================================================================
 SUBROUTINE OpenCon
 END SUBROUTINE OpenCon
 !=======================================================================
 SUBROUTINE OpenUnfInpBEFile ( Un, InFile, RecLen, Error )
 
    ! This routine opens a binary input file. Data is stored in RecLen-byte records.
 
    IMPLICIT                        NONE
 
    INTEGER, INTENT(IN)          :: Un                                           ! Logical unit for the input file.
    CHARACTER(*), INTENT(IN)     :: InFile                                       ! Name of the input file.
    INTEGER, INTENT(IN)          :: RecLen                                       ! Size of records in the input file, in bytes.
    LOGICAL, INTENT(OUT)         :: Error                                        ! Flag to indicate the open failed.
    INTEGER                      :: IOS                                          ! I/O status of OPEN.
 
    ! Open input file.  Make sure it worked.
    OPEN ( Un, FILE=TRIM( InFile ), STATUS='OLD', FORM='UNFORMATTED', &
            ACCESS='DIRECT', RECL=RecLen, IOSTAT=IOS, ACTION='READ'  )
 
    IF ( IOS /= 0 )  THEN
       Error = .TRUE.
    ELSE
       Error = .FALSE.
    END IF
 
 END SUBROUTINE OpenUnfInpBEFile
 !=======================================================================
 SUBROUTINE ProgExit ( StatCode )
 
    ! This routine stops the program.  If the compiler supports the EXIT routine,
    ! pass the program status to it.  Otherwise, do a STOP.
 
    INTEGER, INTENT(IN) :: StatCode ! The status code to pass to the OS.
 
    CALL EXIT ( StatCode )
 
 END SUBROUTINE 
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
 
 END SUBROUTINE UsrAlarm
 !=======================================================================
 SUBROUTINE WrNR ( Str )
 
    ! This routine writes out a string to the screen without following it with a new line.
 
    CHARACTER(*), INTENT(IN)     :: Str                                          ! The string to write to the screen.
 
    WRITE (CU,'(A)',ADVANCE='NO')  Str
 
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
 
 #ifdef USE_DLL_INTERFACE         
 
    ErrStat = ErrID_None
    ErrMsg = ''
 
    ! Load the DLL and get the file address:
 
    DLL%FileAddrX = dlOpen( TRIM(DLL%FileName)//C_NULL_CHAR, RTLD_LAZY )  !the "C_NULL_CHAR" converts the Fortran string to a C-type string (i.e., adds //CHAR(0) to the end)
 
    IF( .NOT. C_ASSOCIATED(DLL%FileAddrX) ) THEN
       ErrStat = ErrID_Fatal
       WRITE(ErrMsg,'(I2)') BITS_IN_ADDR
       ErrMsg  = 'The dynamic library '//TRIM(DLL%FileName)//' could not be loaded. Check that the file '// &
                'exists in the specified location and that it is compiled for '//TRIM(ErrMsg)//'-bit applications.'
       RETURN
    END IF
 
    ! Get the procedure address:
 
    CALL LoadDynamicLibProc ( DLL, ErrStat, ErrMsg )
 #else
 
    ErrStat = ErrID_Fatal
    ErrMsg = ' LoadDynamicLib: Not compiled with -DUSE_DLL_INTERFACE for '//TRIM(OS_Desc)
 
 #endif
 
 END SUBROUTINE LoadDynamicLib
 !=======================================================================
 SUBROUTINE LoadDynamicLibProc ( DLL, ErrStat, ErrMsg )
 
    ! This SUBROUTINE is used to dynamically load a procedure from a DLL.
 
    TYPE (DLL_Type),           INTENT(INOUT)  :: DLL         ! The DLL to be loaded.
    INTEGER(IntKi),            INTENT(  OUT)  :: ErrStat     ! Error status of the operation
    CHARACTER(*),              INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
    INTEGER(IntKi)                            :: i
 
 #ifdef USE_DLL_INTERFACE           
 
    ErrStat = ErrID_None
    ErrMsg = ''
 
    do i=1,NWTC_MAX_DLL_PROC
       if ( len_trim( DLL%ProcName(i) ) > 0 ) then
    
          DLL%ProcAddr(i) = dlSym( DLL%FileAddrX, TRIM(DLL%ProcName(i))//C_NULL_CHAR )  !the "C_NULL_CHAR" converts the Fortran string to a C-type string (i.e., adds //CHAR(0) to the end)
 
          IF(.NOT. C_ASSOCIATED(DLL%ProcAddr(i))) THEN
             ErrStat = ErrID_Fatal + i - 1
             ErrMsg  = 'The procedure '//TRIM(DLL%ProcName(i))//' in file '//TRIM(DLL%FileName)//' could not be loaded.'
             RETURN
          END IF
          
       end if
    end do
    
 #else
 
    ErrStat = ErrID_Fatal
    ErrMsg = ' LoadDynamicLibProc: Not compiled with -DUSE_DLL_INTERFACE for '//TRIM(OS_Desc)
    
 #endif
 
 END SUBROUTINE LoadDynamicLibProc
 !=======================================================================
 SUBROUTINE FreeDynamicLib ( DLL, ErrStat, ErrMsg )
 
    ! This SUBROUTINE is used to free a dynamically loaded DLL (loaded in LoadDynamicLib).
 
    TYPE (DLL_Type),           INTENT(INOUT)  :: DLL         ! The DLL to be freed.
    INTEGER(IntKi),            INTENT(  OUT)  :: ErrStat     ! Error status of the operation
    CHARACTER(*),              INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
    INTEGER(C_INT)                            :: Success     ! Whether or not the call to dlClose was successful
    INTEGER(C_INT), PARAMETER                 :: TRUE  = 0
 
 #ifdef USE_DLL_INTERFACE           
 
    ! Close the library:
 
    IF( .NOT. C_ASSOCIATED(DLL%FileAddrX) ) RETURN
    Success = dlClose( DLL%FileAddrX ) !The function dlclose() returns 0 on success, and nonzero on error.
 
    IF ( Success /= TRUE ) THEN !bjj: note that this is not the same as LOGICAL .TRUE.
       ErrStat = ErrID_Fatal
       ErrMsg  = 'The dynamic library could not be freed.'
       RETURN
    ELSE
       ErrStat = ErrID_None
       ErrMsg = ''
       DLL%FileAddrX = C_NULL_PTR
    END IF
    
 #else
 
    ErrStat = ErrID_Fatal
    ErrMsg = ' FreeDynamicLib: Not compiled with -DUSE_DLL_INTERFACE for '//TRIM(OS_Desc)
       
 #endif
 
 END SUBROUTINE FreeDynamicLib
 !=======================================================================
 END MODULE SysSubs
 