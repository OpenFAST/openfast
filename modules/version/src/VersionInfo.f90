!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2015-2016  National Renewable Energy Laboratory
! Copyright (C) 2016-2017  Envision Energy USA, LTD
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
MODULE VersionInfo

   use NWTC_Library
   implicit none

contains

FUNCTION QueryGitVersion()

   CHARACTER(200) :: QueryGitVersion

! The Visual Studio project sets the path for where to find the header file with version info
#ifdef GIT_INCLUDE_FILE
#include GIT_INCLUDE_FILE
#endif

#ifdef GIT_VERSION_INFO
   QueryGitVersion = GIT_VERSION_INFO
#else
   QueryGitVersion = 'unversioned'
#endif

   RETURN
END FUNCTION QueryGitVersion


!----------------------------------------------------------------------------------------------------------------------------------
!> This function returns a string describing the glue code and some of the compilation options we're using.
FUNCTION GetVersion(ThisProgVer, Cmpl4SFun, Cmpl4LV)

   ! Passed Variables:

   TYPE(ProgDesc), INTENT( IN    ) :: ThisProgVer     !< program name/date/version description
   LOGICAL, INTENT(IN), OPTIONAL   :: Cmpl4SFun
   LOGICAL, INTENT(IN), OPTIONAL   :: Cmpl4LV
   CHARACTER(1024)                 :: GetVersion      !< String containing a description of the compiled precision.

   CHARACTER(200)                  :: git_commit

   GetVersion = TRIM(GetNVD(ThisProgVer))//', compiled on '//__DATE__//' at '//__TIME__

   if (present(Cmpl4SFun)) then
      IF ( Cmpl4SFun )  THEN     ! FAST has been compiled as an S-Function for Simulink
         GetVersion = TRIM(GetVersion)//' as a DLL S-Function for Simulink'
      ENDIF
   endif
   
   if (present(Cmpl4LV)) then
      IF ( Cmpl4LV )  THEN     ! FAST has been compiled as a DLL for Labview
         GetVersion = TRIM(GetVersion)//' as a DLL for LabVIEW'
      ENDIF
   endif
   
   GetVersion = TRIM(GetVersion)//' as a '//TRIM(Num2LStr(BITS_IN_ADDR))//'-bit application using'

   ! determine precision

      IF ( ReKi == SiKi )  THEN     ! Single precision
         GetVersion = TRIM(GetVersion)//' single'
      ELSEIF ( ReKi == R8Ki )  THEN ! Double precision
         GetVersion = TRIM(GetVersion)// ' double'
      ELSE                          ! Unknown precision
         GetVersion = TRIM(GetVersion)//' unknown'
      ENDIF


!   GetVersion = TRIM(GetVersion)//' precision with '//OS_Desc
   GetVersion = TRIM(GetVersion)//' precision'

   ! add git info
   git_commit = QueryGitVersion()
   GetVersion = TRIM(GetVersion)//' at commit '//git_commit

   RETURN
END FUNCTION GetVersion

!=======================================================================
!>
   SUBROUTINE DispCompileRuntimeInfo(name)
#ifdef _OPENMP    
      USE OMP_LIB 
#endif
#ifdef HAS_FORTRAN2008_FEATURES
      USE iso_fortran_env, ONLY: compiler_version
#endif
      CHARACTER(*), intent(in) :: name
      CHARACTER(200) :: compiler_version_str
      CHARACTER(200) :: git_commit, architecture, compiled_precision
      CHARACTER(200) :: execution_date, execution_time, execution_zone

!     name = ProgName
      git_commit = QueryGitVersion()
      architecture = TRIM(Num2LStr(BITS_IN_ADDR))//' bit'
      IF (ReKi == SiKi) THEN
         compiled_precision = 'single'
      ELSE IF (ReKi == R8Ki) THEN
         compiled_precision = 'double'
      ELSE
         compiled_precision = 'unknown'
      END IF

#if defined(HAS_FORTRAN2008_FEATURES)
      compiler_version_str = compiler_version()
#elif defined(__INTEL_COMPILER)
      compiler_version_str = 'Intel(R) Fortran Compiler '//num2lstr(__INTEL_COMPILER)
#else
      compiler_version_str = OS_Desc
#endif

      CALL WrScr(trim(name)//'-'//trim(git_commit))
      CALL WrScr('Compile Info:')
      call wrscr(' - Compiler: '//trim(compiler_version_str))
      CALL WrScr(' - Architecture: '//trim(architecture))
      CALL WrScr(' - Precision: '//trim(compiled_precision))
#ifdef _OPENMP   
      !$OMP PARALLEL default(shared)
      if (omp_get_thread_num()==0) then
         call WrScr(' - OpenMP: Yes, number of threads: '//trim(Num2LStr(omp_get_num_threads()))//'/'//trim(Num2LStr(omp_get_max_threads())))
      endif
      !$OMP END PARALLEL 
#else
   call WrScr(' - OpenMP: No')
#endif
      CALL WrScr(' - Date: '//__DATE__)
      CALL WrScr(' - Time: '//__TIME__)
      ! call wrscr(' - Options: '//trim(compiler_options()))

      CALL DATE_AND_TIME(execution_date, execution_time, execution_zone)

      CALL WrScr('Execution Info:')
      CALL WrScr(' - Date: '//TRIM(execution_date(5:6)//'/'//execution_date(7:8)//'/'//execution_date(1:4)))
      CALL WrScr(' - Time: '//TRIM(execution_time(1:2)//':'//execution_time(3:4)//':'//execution_time(5:6))//TRIM(execution_zone))
      CALL WrScr('')

   END SUBROUTINE DispCompileRuntimeInfo
   
!=======================================================================
!> This subroutine checks for command-line arguments.
   SUBROUTINE CheckArgs ( Arg1, ErrStat, Arg2, Flag, InputArgArray )

      ! Argument declarations:
   CHARACTER(*), INTENT(INOUT)           :: Arg1               !< The first non-flag argument; generally, the name of the input file.
   INTEGER,      INTENT(  OUT), OPTIONAL :: ErrStat            !< An optional argument for catching errors; if present, program does not abort on error.
   CHARACTER(*), INTENT(  OUT), OPTIONAL :: Arg2               !< An optional 2nd non-flag argument.
   CHARACTER(*), INTENT(  OUT), OPTIONAL :: Flag               !< An optional flag argument; the first argument starting with a switch character. 
   CHARACTER(*), INTENT(IN   ), DIMENSION(:), OPTIONAL :: InputArgArray  !< An optional argument containing the arguments to parse; primarily used for unit testing.

      ! Local declarations:
   INTEGER                                    :: I, J          ! Iterator variables
   CHARACTER(1024)                            :: Arg, FlagIter
   CHARACTER(1024), DIMENSION(:), ALLOCATABLE :: ArgArray, TempArray, Flags
   LOGICAL :: FirstArgumentSet, SecondArgumentSet

   FirstArgumentSet = .FALSE.
   SecondArgumentSet = .FALSE.
                                        
   IF ( PRESENT(Arg2) ) Arg2 = ""
   IF ( PRESENT(Flag) ) Flag = ""
   
      ! Save all arguments in a single argument array; this is primarily used to enable unit testing
   IF ( PRESENT(InputArgArray) ) THEN
      ALLOCATE( ArgArray( SIZE(InputArgArray) ) )
      ArgArray = InputArgArray
   ELSE
      ALLOCATE( ArgArray( COMMAND_ARGUMENT_COUNT() ) )
      DO I = 1, SIZE(ArgArray)
         CALL GET_COMMAND_LINE_ARG( I, ArgArray(I) )
      END DO
   END IF

      ! Early return if no arguments and no default input file given
   IF ( SIZE(ArgArray) == 0 .AND. LEN( TRIM(Arg1) ) == 0 ) THEN
      CALL INVALID_SYNTAX( 'no command-line arguments given.' )
      CALL CLEANUP()
      RETURN
   END IF

      ! Split arguments into flags and non-flags
   ALLOCATE( Flags(0) )
   DO I = 1, SIZE(ArgArray)
      Arg = TRIM(ArgArray(I))
      IF ( IsFlag(Arg) ) THEN
            ! This is how we can dynamically resize an array in Fortran...
            ! Dont do this where performance matters.
         ALLOCATE( TempArray( SIZE(Flags) + 1 ) )
         DO J = 1, SIZE(Flags)
            TempArray(J) = Flags(J)
         END DO
         TempArray(SIZE(Flags) + 1) = TRIM(Arg)
         DEALLOCATE(Flags)
         CALL MOVE_ALLOC(TempArray, Flags)
      ELSE IF ( .NOT. FirstArgumentSet ) THEN
         Arg1 = TRIM(Arg)
         FirstArgumentSet = .TRUE.
      ELSE IF ( .NOT. SecondArgumentSet ) THEN
         Arg2 = TRIM(Arg)
         SecondArgumentSet = .True.
      ELSE
         CALL INVALID_SYNTAX( 'too many command-line arguments given.' )
         CALL CLEANUP()
         RETURN
      END IF
   END DO

   DO I = 1, SIZE(Flags)

      FlagIter = Flags(I)(2:) ! This results in the flag without the switch character
      CALL Conv2UC( FlagIter )
      IF ( PRESENT(Flag) ) Flag = FlagIter

      SELECT CASE ( TRIM(FlagIter) )

      CASE ('H')
         CALL DispCopyrightLicense( ProgName )
         CALL DispCompileRuntimeInfo( ProgName )
         CALL NWTC_DisplaySyntax( Arg1, ProgName )
         IF ( PRESENT( ErrStat ) ) ErrStat = ErrID_None
         CALL CLEANUP()
         RETURN

      CASE ('V', 'VERSION')
         CALL DispCopyrightLicense( ProgName )
         CALL DispCompileRuntimeInfo( ProgName )
         IF ( PRESENT( ErrStat ) ) ErrStat = ErrID_None
         CALL CLEANUP()
         RETURN

      CASE ('RESTART')
         IF ( FirstArgumentSet .AND. .NOT. SecondArgumentSet ) THEN
            Arg2 = Arg1
            Arg1 = ""
         END IF
         IF ( .NOT. FirstArgumentSet .AND. .NOT. SecondArgumentSet ) THEN
            CALL INVALID_SYNTAX( 'the restart capability requires at least one argument: <input_file (OPTIONAL)> -restart <checkpoint_file>' )
            CALL CLEANUP()
               RETURN
         END IF

      CASE ('VTKLIN')
         IF ( FirstArgumentSet .AND. .NOT. SecondArgumentSet ) THEN
            Arg2 = Arg1
            Arg1 = ""
         END IF
         IF ( .NOT. FirstArgumentSet .AND. .NOT. SecondArgumentSet ) THEN
            CALL INVALID_SYNTAX( 'the restart capability for vtk mode shapes requires at least one argument: <input_file (OPTIONAL)> -vtklin <visualization_input_file>' )
            CALL CLEANUP()
               RETURN
         END IF

      CASE DEFAULT
         CALL INVALID_SYNTAX( 'unknown command-line argument given: '//TRIM(FlagIter) )
         CALL CLEANUP()
         RETURN
                                                
      END SELECT

   END DO

   IF ( PRESENT( ErrStat ) ) ErrStat = ErrID_None
   CALL CLEANUP()

                  RETURN

   CONTAINS
      SUBROUTINE CLEANUP()
         IF ( ALLOCATED(ArgArray) ) DEALLOCATE(ArgArray)
         IF ( ALLOCATED(Flags) ) DEALLOCATE(Flags)
         IF ( ALLOCATED(TempArray) ) DEALLOCATE(TempArray)
      END SUBROUTINE

      SUBROUTINE INVALID_SYNTAX(ErrorMessage)

         CHARACTER(*), INTENT(IN) :: ErrorMessage

         CALL DispCopyrightLicense( ProgName )
         CALL DispCompileRuntimeInfo( ProgName )
         CALL NWTC_DisplaySyntax( Arg1, ProgName )
         CALL ProgAbort( ' Invalid syntax: '//TRIM(ErrorMessage), PRESENT(ErrStat) )
         IF ( PRESENT(ErrStat) ) ErrStat = ErrID_Fatal

      END SUBROUTINE

      SUBROUTINE GET_COMMAND_LINE_ARG(ArgIndex, ArgGiven)

         INTEGER, INTENT(IN) :: ArgIndex           !< Index location of the argument to get.
         CHARACTER(*), INTENT(OUT) :: ArgGiven  !< The gotten command-line argument.
         INTEGER :: Error                          !< Indicates if there was an error getting an argument.

         CALL GET_COMMAND_ARGUMENT( ArgIndex, ArgGiven, STATUS=Error )
         ArgGiven = TRIM(ArgGiven)
         IF ( Error /= 0 )  THEN
            CALL ProgAbort ( ' Error getting command-line argument #'//TRIM( Int2LStr( ArgIndex ) )//'.', PRESENT(ErrStat) )
            IF ( PRESENT(ErrStat) ) ErrStat = ErrID_Fatal
         END IF

      END SUBROUTINE

      FUNCTION IsFlag(ArgString)

         CHARACTER(*), INTENT(IN) :: ArgString
         LOGICAL :: IsFlag

         IF ( ArgString(1:1) == SwChar .OR. ArgString(1:1) == '-' ) THEN
            IsFlag = .TRUE.
         ELSE
            IsFlag = .FALSE.
         END IF

      END FUNCTION

   END SUBROUTINE CheckArgs
   
END MODULE
