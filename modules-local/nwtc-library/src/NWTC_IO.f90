MODULE NWTC_IO


   ! This module contains I/O-related variables and routines with non-system-specific logic.


   ! It contains the following routines:

   !     SUBROUTINE CheckArgs     ( InputFile [, ErrStat] )
   !     SUBROUTINE CheckIOS      ( IOS, Fil, Variable, VarType [, TrapErrors] )
   !     SUBROUTINE CloseEcho     ( )
   !     SUBROUTINE Conv2UC       ( Str )
   !     FUNCTION   CountWords    ( Line )
   !     FUNCTION   CurDate       ( )
   !     FUNCTION   CurTime       ( )
   !     SUBROUTINE DispNVD       ( )
   !     FUNCTION   Flt2LStr      ( FltNum )
   !     SUBROUTINE GetNewUnit    ( UnIn )
   !     SUBROUTINE GetPath       ( GivenFil, PathName )
   !     SUBROUTINE GetRoot       ( GivenFil, RootName )
   !     SUBROUTINE GetTokens     ( Line, NumTok, Tokens, Error )
   !     SUBROUTINE GetWords      ( Line, Words, NumWords )
   !     FUNCTION   Int2LStr      ( Intgr )
   !     SUBROUTINE NameOFile     ( InArg, OutExten, OutFile )
   !     SUBROUTINE NormStop      ( )
   !     SUBROUTINE OpenBin       ( Un, OutFile, RecLen [, ErrStat] )
   !     SUBROUTINE OpenBInpFile  ( Un, InFile [, ErrStat] )
   !     SUBROUTINE OpenEcho      ( Un, InFile [, ErrStat] )
   !     SUBROUTINE OpenFInpFile  ( Un, InFile [, ErrStat] )
   !     SUBROUTINE OpenFOutFile  ( Un, OutFile [, ErrStat] )
   !     SUBROUTINE OpenFUnkFile  ( Un, OutFile, FailAbt, Failed, Exists [, ErrStat] )
   !     SUBROUTINE OpenUInfile   ( Un, InFile [, ErrStat] )
   !     SUBROUTINE OpenUInBEFile ( Un, InFile, RecLen [, ErrStat] )
   !     SUBROUTINE OpenUOutfile  ( Un, OutFile [, ErrStat] )
   !     FUNCTION   PathIsRelative( GivenFil )
   !     SUBROUTINE PremEOF       ( Fil , Variable [, TrapErrors] )
   !     SUBROUTINE ProgAbort     ( Message [, TrapErrors] )
   !     SUBROUTINE ProgWarn      ( Message )
   !     SUBROUTINE ReadAry       ( UnIn, Fil, Ary, AryLen, AryName, AryDescr [, ErrStat] )         ! Generic interface for ReadCAry, ReadIAry, ReadLAry, and ReadRAry.
   !     SUBROUTINE ReadAryLines  ( UnIn, Fil, Ary, AryLen, AryName, AryDescr [, ErrStat] )         ! Generic interface for ReadCAryLines, and ReadRAryLines.
   !     SUBROUTINE ReadCAry      ( UnIn, Fil, CharAry, AryLen, AryName, AryDescr [, ErrStat] )
   !     SUBROUTINE ReadCAryLines ( UnIn, Fil, CharAry, AryLen, AryName, AryDescr [, ErrStat] )
   !     SUBROUTINE ReadCom       ( UnIn, Fil, ComName [, ErrStat] )
   !     SUBROUTINE ReadCVar      ( UnIn, Fil, CharVar, VarName, VarDescr [, ErrStat] )
   !     SUBROUTINE ReadIAry      ( UnIn, Fil, IntAry, AryLen, AryName, AryDescr [, ErrStat] )
   !     SUBROUTINE ReadIVar      ( UnIn, Fil, IntVar, VarName, VarDescr [, ErrStat] )
   !     SUBROUTINE ReadLAry      ( UnIn, Fil, LogAry, AryLen, AryName, AryDescr [, ErrStat] )
   !     SUBROUTINE ReadLVar      ( UnIn, Fil, LogVar, VarName, VarDescr [, ErrStat] )
   !     SUBROUTINE ReadNum       ( UnIn, Fil, Word, VarName, ErrStat )
   !     SUBROUTINE ReadOutputList( UnIn, Fil, CharAry, AryLenRead, AryName, AryDescr, ErrStat )
   !     SUBROUTINE ReadRAry      ( UnIn, Fil, RealAry, AryLen, AryName, AryDescr [, ErrStat] )
   !     SUBROUTINE ReadRAryLines ( UnIn, Fil, RealAry, AryLen, AryName, AryDescr [, ErrStat] )
   !     SUBROUTINE ReadRVar      ( UnIn, Fil, RealVar, VarName, VarDescr [, ErrStat] )
   !     SUBROUTINE ReadStr       ( UnIn, Fil, CharVar, VarName, VarDescr [, ErrStat] )
   !     SUBROUTINE ReadVar       ( UnIn, Fil, Var, VarName, VarDescr [, ErrStat] )                 ! Generic interface for ReadCVar, ReadIVar, ReadLVar, and ReadRVar.
   !     SUBROUTINE WaitTime      ( WaitSecs )
   !     SUBROUTINE WrPr          ( Str )
   !     SUBROUTINE WrFileNR      ( Unit, Str )
   !     SUBROUTINE WrML          ( Str )
   !     SUBROUTINE WrScr1        ( Str )


   USE                             SysSubs


!=======================================================================


      ! Global I/O-related variables.

   INTEGER, PARAMETER           :: FlgType  = 1                                 ! Switch for telling if a variable is a flag.
   INTEGER, PARAMETER           :: NumType  = 2                                 ! Switch for telling if a variable is a number.
   INTEGER, PARAMETER           :: StrType  = 3                                 ! Switch for telling if a variable is a string.
   INTEGER                      :: UnEc     = 19                                ! I/O unit number for the echo file.

   LOGICAL                      :: Beep     = .TRUE.                            ! Flag that specifies whether or not to beep for error messages and program terminations.
   LOGICAL                      :: Echo     = .FALSE.                           ! Flag that specifies whether or not to produce an echo file.

   CHARACTER(23)                :: NWTCName = 'NWTC Subroutine Library'         ! The name of the NWTC subroutine library.
   CHARACTER(29)                :: NWTCVer  = ' (v1.04.01, 21-Feb-2012)'        ! The version (including date) of the NWTC Subroutine Library.
   CHARACTER(20)                :: ProgName = ' '                               ! The name of the calling program.
   CHARACTER(99)                :: ProgVer                                      ! The version (including date) of the calling program.
   CHARACTER(1), PARAMETER      :: Tab      = CHAR( 9 )                         ! The tab character.


!=======================================================================


      ! Create interface for a generic AllocAry that actually uses specific routines.

   INTERFACE AllocAry
      MODULE PROCEDURE AllCAry1
      MODULE PROCEDURE AllCAry2
      MODULE PROCEDURE AllCAry3
   !   MODULE PROCEDURE AllCAry4                               Not yet coded.
      MODULE PROCEDURE AllIAry1
      MODULE PROCEDURE AllIAry2
      MODULE PROCEDURE AllIAry3
   !   MODULE PROCEDURE AllIAry4                               Not yet coded.
      MODULE PROCEDURE AllLAry1
      MODULE PROCEDURE AllLAry2
      MODULE PROCEDURE AllLAry3
   !   MODULE PROCEDURE AllLAry4                               Not yet coded.
      MODULE PROCEDURE AllRAry1
      MODULE PROCEDURE AllRAry2
      MODULE PROCEDURE AllRAry3
      MODULE PROCEDURE AllRAry4
   END INTERFACE


      ! Create interface for a generic ReadVar that actually uses specific routines.

   INTERFACE ReadVar
      MODULE PROCEDURE ReadCVar
      MODULE PROCEDURE ReadIVar
      MODULE PROCEDURE ReadLVar
      MODULE PROCEDURE ReadRVar
   END INTERFACE


      ! Create interface for a generic ReadAry that actually uses specific routines.

   INTERFACE ReadAry
      MODULE PROCEDURE ReadCAry
      MODULE PROCEDURE ReadIAry
      MODULE PROCEDURE ReadLAry
      MODULE PROCEDURE ReadRAry
   END INTERFACE

   INTERFACE ReadAryLines
      MODULE PROCEDURE ReadCAryLines
!     MODULE PROCEDURE ReadIAryLines         ! Not coded yet
!     MODULE PROCEDURE ReadLAryLines         ! Not coded yet
      MODULE PROCEDURE ReadRAryLines
   END INTERFACE


      ! Create interface for a generic Num2LStr that actually uses specific routines.

   INTERFACE Num2LStr
      MODULE PROCEDURE Int2LStr
      MODULE PROCEDURE Flt2LStr
   END INTERFACE


CONTAINS

!=======================================================================
   SUBROUTINE AllCAry1 ( Ary, AryDim, Descr, ErrStat )

      ! This routine allocates a 1-D CHARACTER array.


      ! Argument declarations.

   INTEGER, INTENT(IN)          :: AryDim                                      ! The size of the array.
   INTEGER, INTENT(OUT),OPTIONAL:: ErrStat                                     ! Error status; if present, program does not abort on error


   CHARACTER(*), ALLOCATABLE    :: Ary    (:)                                 ! Array to be allocated
   CHARACTER(*), INTENT(IN)     :: Descr                                       ! Brief array description.


      ! Local declarations.

   INTEGER                      :: Sttus                                       ! Status of allocation attempt.



   ALLOCATE ( Ary(AryDim) , STAT=Sttus )


   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort ( ' Error allocating memory for the '//TRIM( Descr )//' array.', PRESENT(ErrStat) )
   END IF

   IF ( PRESENT(ErrStat) ) ErrStat = Sttus

   RETURN
   END SUBROUTINE AllCAry1 ! ( Ary, AryDim, Descr )
!=======================================================================
   SUBROUTINE AllCAry2 (  Ary, AryDim1, AryDim2, Descr, ErrStat )

      ! This routine allocates a 2-D CHARACTER array.


      ! Argument declarations.

   INTEGER, INTENT(IN)          :: AryDim1                                     ! The size of the first dimension of the array.
   INTEGER, INTENT(IN)          :: AryDim2                                     ! The size of the second dimension of the array.
   INTEGER, INTENT(OUT),OPTIONAL:: ErrStat                                     ! Error status; if present, program does not abort on error

   CHARACTER(*), ALLOCATABLE    :: Ary    (:,:)                                ! Array to be allocated
   CHARACTER(*), INTENT(IN)     :: Descr                                       ! Brief array description.


      ! Local declarations.

   INTEGER                      :: Sttus                                       ! Status of allocation attempt.



   ALLOCATE ( Ary(AryDim1,AryDim2) , STAT=Sttus )

   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort ( ' Error allocating memory for the '//TRIM( Descr )//' array.', PRESENT(ErrStat) )
   END IF

   IF ( PRESENT(ErrStat) ) ErrStat = Sttus


   RETURN
   END SUBROUTINE AllCAry2 ! (  Ary, AryDim1, AryDim2, Descr )
!=======================================================================
   SUBROUTINE AllCAry3 (  Ary, AryDim1, AryDim2, AryDim3, Descr, ErrStat )


      ! This routine allocates a 3-D CHARACTER array.


      ! Argument declarations.

   INTEGER, INTENT(IN)          :: AryDim1                                     ! The size of the first dimension of the array.
   INTEGER, INTENT(IN)          :: AryDim2                                     ! The size of the second dimension of the array.
   INTEGER, INTENT(IN)          :: AryDim3                                     ! The size of the third dimension of the array.
   INTEGER, INTENT(OUT),OPTIONAL:: ErrStat                                     ! Error status; if present, program does not abort on error


   CHARACTER(*), ALLOCATABLE    :: Ary    (:,:,:)                              ! Array to be allocated
   CHARACTER(*), INTENT(IN)     :: Descr                                       ! Brief array description.


      ! Local declarations.

   INTEGER                      :: Sttus                                       ! Status of allocation attempt.



   ALLOCATE ( Ary(AryDim1,AryDim2,AryDim3) , STAT=Sttus )

   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort ( ' Error allocating memory for the '//TRIM( Descr )//' array.', PRESENT(ErrStat) )
   END IF

   IF ( PRESENT(ErrStat) ) ErrStat = Sttus


   RETURN
   END SUBROUTINE AllCAry3 ! (  Ary, AryDim1, AryDim2, AryDim3, Descr )
!=======================================================================
   SUBROUTINE AllIAry1 ( Ary, AryDim, Descr, ErrStat )


      ! This routine allocates a 1-D INTEGER array.


      ! Argument declarations.

   INTEGER, ALLOCATABLE         :: Ary    (:)                                  ! Array to be allocated
   INTEGER, INTENT(IN)          :: AryDim                                      ! The size of the array.
   INTEGER, INTENT(OUT),OPTIONAL:: ErrStat                                     ! Error status; if present, program does not abort on error

   CHARACTER(*), INTENT(IN)     :: Descr                                       ! Brief array description.


      ! Local declarations.

   INTEGER                      :: Sttus                                       ! Status of allocation attempt.



   ALLOCATE ( Ary(AryDim) , STAT=Sttus )

   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort ( ' Error allocating memory for the '//TRIM( Descr )//' array.', PRESENT(ErrStat) )
   END IF

   IF ( PRESENT(ErrStat) ) ErrStat = Sttus

   RETURN
   END SUBROUTINE AllIAry1 ! ( Ary, AryDim, Descr )
!=======================================================================
   SUBROUTINE AllIAry2 (  Ary, AryDim1, AryDim2, Descr, ErrStat )


      ! This routine allocates a 2-D INTEGER array.


      ! Argument declarations.

   INTEGER, ALLOCATABLE         :: Ary    (:,:)                                ! Array to be allocated
   INTEGER, INTENT(IN)          :: AryDim1                                     ! The size of the first dimension of the array.
   INTEGER, INTENT(IN)          :: AryDim2                                     ! The size of the second dimension of the array.
   INTEGER, INTENT(OUT),OPTIONAL:: ErrStat                                     ! Error status; if present, program does not abort on error

   CHARACTER(*), INTENT(IN)     :: Descr                                       ! Brief array description.


      ! Local declarations.

   INTEGER                      :: Sttus                                       ! Status of allocation attempt.



   ALLOCATE ( Ary(AryDim1,AryDim2) , STAT=Sttus )

   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort ( ' Error allocating memory for the '//TRIM( Descr )//' array.', PRESENT(ErrStat) )
   END IF

   IF ( PRESENT(ErrStat) ) ErrStat = Sttus

   RETURN
   END SUBROUTINE AllIAry2 ! (  Ary, AryDim1, AryDim2, Descr )
!=======================================================================
   SUBROUTINE AllIAry3 (  Ary, AryDim1, AryDim2, AryDim3, Descr, ErrStat )


      ! This routine allocates a 3-D INTEGER array.


      ! Argument declarations.

   INTEGER, ALLOCATABLE         :: Ary    (:,:,:)                              ! Array to be allocated
   INTEGER, INTENT(IN)          :: AryDim1                                     ! The size of the first dimension of the array.
   INTEGER, INTENT(IN)          :: AryDim2                                     ! The size of the second dimension of the array.
   INTEGER, INTENT(IN)          :: AryDim3                                     ! The size of the third dimension of the array.
   INTEGER, INTENT(OUT),OPTIONAL:: ErrStat                                     ! Error status; if present, program does not abort on error

   CHARACTER(*), INTENT(IN)     :: Descr                                       ! Brief array description.


      ! Local declarations.

   INTEGER                      :: Sttus                                       ! Status of allocation attempt.



   ALLOCATE ( Ary(AryDim1,AryDim2,AryDim3) , STAT=Sttus )

   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort ( ' Error allocating memory for the '//TRIM( Descr )//' array.', PRESENT(ErrStat) )
   END IF

   IF ( PRESENT(ErrStat) ) ErrStat = Sttus


   RETURN
   END SUBROUTINE AllIAry3 ! (  Ary, AryDim1, AryDim2, AryDim3, Descr )
!=======================================================================
   SUBROUTINE AllLAry1 ( Ary, AryDim, Descr, ErrStat )


      ! This routine allocates a 1-D LOGICAL array.


      ! Argument declarations.

   INTEGER, INTENT(IN)          :: AryDim                                      ! The size of the array.
   INTEGER, INTENT(OUT),OPTIONAL:: ErrStat                                     ! Error status; if present, program does not abort on error

   LOGICAL, ALLOCATABLE         :: Ary    (:)                                  ! Array to be allocated

   CHARACTER(*), INTENT(IN)     :: Descr                                       ! Brief array description.


      ! Local declarations.

   INTEGER                      :: Sttus                                       ! Status of allocation attempt.



   ALLOCATE ( Ary(AryDim) , STAT=Sttus )

   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort ( ' Error allocating memory for the '//TRIM( Descr )//' array.', PRESENT(ErrStat) )
   END IF

   IF ( PRESENT(ErrStat) ) ErrStat = Sttus


   RETURN
   END SUBROUTINE AllLAry1 ! ( Ary, AryDim, Descr )
!=======================================================================
   SUBROUTINE AllLAry2 (  Ary, AryDim1, AryDim2, Descr, ErrStat )


      ! This routine allocates a 2-D LOGICAL array.


      ! Argument declarations.

   INTEGER, INTENT(IN)          :: AryDim1                                     ! The size of the first dimension of the array.
   INTEGER, INTENT(IN)          :: AryDim2                                     ! The size of the second dimension of the array.
   INTEGER, INTENT(OUT),OPTIONAL:: ErrStat                                     ! Error status; if present, program does not abort on error

   LOGICAL, ALLOCATABLE         :: Ary    (:,:)                                ! Array to be allocated

   CHARACTER(*), INTENT(IN)     :: Descr                                       ! Brief array description.


      ! Local declarations.

   INTEGER                      :: Sttus                                       ! Status of allocation attempt.



   ALLOCATE ( Ary(AryDim1,AryDim2) , STAT=Sttus )

   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort ( ' Error allocating memory for the '//TRIM( Descr )//' array.', PRESENT(ErrStat) )
   END IF

   IF ( PRESENT(ErrStat) ) ErrStat = Sttus


   RETURN
   END SUBROUTINE AllLAry2 ! (  Ary, AryDim1, AryDim2, Descr )
!=======================================================================
   SUBROUTINE AllLAry3 (  Ary, AryDim1, AryDim2, AryDim3, Descr, ErrStat )


      ! This routine allocates a 3-D LOGICAL array.


      ! Argument declarations.

   INTEGER, INTENT(IN)          :: AryDim1                                     ! The size of the first dimension of the array.
   INTEGER, INTENT(IN)          :: AryDim2                                     ! The size of the second dimension of the array.
   INTEGER, INTENT(IN)          :: AryDim3                                     ! The size of the third dimension of the array.
   INTEGER, INTENT(OUT),OPTIONAL:: ErrStat                                     ! Error status; if present, program does not abort on error

   LOGICAL, ALLOCATABLE         :: Ary    (:,:,:)                              ! Array to be allocated

   CHARACTER(*), INTENT(IN)     :: Descr                                       ! Brief array description.


      ! Local declarations.

   INTEGER                      :: Sttus                                       ! Status of allocation attempt.



   ALLOCATE ( Ary(AryDim1,AryDim2,AryDim3) , STAT=Sttus )

   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort ( ' Error allocating memory for the '//TRIM( Descr )//' array.', PRESENT(ErrStat) )
   END IF

   IF ( PRESENT(ErrStat) ) ErrStat = Sttus


   RETURN
   END SUBROUTINE AllLAry3 ! (  Ary, AryDim1, AryDim2, AryDim3, Descr )
!=======================================================================
   SUBROUTINE AllRAry1 ( Ary, AryDim, Descr, ErrStat )


      ! This routine allocates a 1-D REAL array.


      ! Argument declarations.

   REAL(ReKi), ALLOCATABLE      :: Ary    (:)                                  ! Array to be allocated

   INTEGER, INTENT(IN)          :: AryDim                                      ! The size of the array.
   INTEGER, INTENT(OUT),OPTIONAL:: ErrStat                                     ! Error status; if present, program does not abort on error

   CHARACTER(*), INTENT(IN)     :: Descr                                       ! Brief array description.


      ! Local declarations.

   INTEGER                      :: Sttus                                       ! Status of allocation attempt.



   ALLOCATE ( Ary(AryDim) , STAT=Sttus )

   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort ( ' Error allocating memory for the '//TRIM( Descr )//' array.', PRESENT(ErrStat) )
   END IF

   IF ( PRESENT(ErrStat) ) ErrStat = Sttus


   RETURN
   END SUBROUTINE AllRAry1 ! ( Ary, AryDim, Descr )
!=======================================================================
   SUBROUTINE AllRAry2 (  Ary, AryDim1, AryDim2, Descr, ErrStat )


      ! This routine allocates a 2-D REAL array.


      ! Argument declarations.

   REAL(ReKi), ALLOCATABLE      :: Ary    (:,:)                                ! Array to be allocated

   INTEGER, INTENT(IN)          :: AryDim1                                     ! The size of the first dimension of the array.
   INTEGER, INTENT(IN)          :: AryDim2                                     ! The size of the second dimension of the array.
   INTEGER, INTENT(OUT),OPTIONAL:: ErrStat                                     ! Error status; if present, program does not abort on error

   CHARACTER(*), INTENT(IN)     :: Descr                                       ! Brief array description.


      ! Local declarations.

   INTEGER                      :: Sttus                                       ! Status of allocation attempt.



   ALLOCATE ( Ary(AryDim1,AryDim2) , STAT=Sttus )

   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort ( ' Error allocating memory for the '//TRIM( Descr )//' array.', PRESENT(ErrStat) )
   END IF

   IF ( PRESENT(ErrStat) ) ErrStat = Sttus


   RETURN
   END SUBROUTINE AllRAry2 ! (  Ary, AryDim1, AryDim2, Descr )
!=======================================================================
   SUBROUTINE AllRAry3 (  Ary, AryDim1, AryDim2, AryDim3, Descr, ErrStat )


      ! This routine allocates a 3-D REAL array.


      ! Argument declarations.

   REAL(ReKi), ALLOCATABLE      :: Ary    (:,:,:)                              ! Array to be allocated

   INTEGER, INTENT(IN)          :: AryDim1                                     ! The size of the first dimension of the array.
   INTEGER, INTENT(IN)          :: AryDim2                                     ! The size of the second dimension of the array.
   INTEGER, INTENT(IN)          :: AryDim3                                     ! The size of the third dimension of the array.
   INTEGER, INTENT(OUT),OPTIONAL:: ErrStat                                     ! Error status; if present, program does not abort on error

   CHARACTER(*), INTENT(IN)     :: Descr                                       ! Brief array description.


      ! Local declarations.

   INTEGER                      :: Sttus                                       ! Status of allocation attempt.



   ALLOCATE ( Ary(AryDim1,AryDim2,AryDim3) , STAT=Sttus )

   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort ( ' Error allocating memory for the '//TRIM( Descr )//' array.', PRESENT(ErrStat) )
   END IF

   IF ( PRESENT(ErrStat) ) ErrStat = Sttus


   RETURN
  END SUBROUTINE AllRAry3 ! (  Ary, AryDim1, AryDim2, AryDim3, Descr )
!=======================================================================
   SUBROUTINE AllRAry4 (  Ary, AryDim1, AryDim2, AryDim3, AryDim4, Descr, ErrStat )


      ! This routine allocates a 4-D REAL array.


      ! Argument declarations.

   REAL(ReKi), ALLOCATABLE      :: Ary    (:,:,:,:)                            ! Array to be allocated

   INTEGER, INTENT(IN)          :: AryDim1                                     ! The size of the first dimension of the array.
   INTEGER, INTENT(IN)          :: AryDim2                                     ! The size of the second dimension of the array.
   INTEGER, INTENT(IN)          :: AryDim3                                     ! The size of the third dimension of the array.
   INTEGER, INTENT(IN)          :: AryDim4                                     ! The size of the fourth dimension of the array.
   INTEGER, INTENT(OUT),OPTIONAL:: ErrStat                                     ! Error status; if present, program does not abort on error

   CHARACTER(*), INTENT(IN)     :: Descr                                       ! Brief array description.


      ! Local declarations.

   INTEGER                      :: Sttus                                       ! Status of allocation attempt.



   ALLOCATE ( Ary(AryDim1,AryDim2,AryDim3,AryDim4) , STAT=Sttus )

   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort ( ' Error allocating memory for the '//TRIM( Descr )//' array.', PRESENT(ErrStat) )
   END IF

   IF ( PRESENT(ErrStat) ) ErrStat = Sttus


   RETURN
  END SUBROUTINE AllRAry4 ! (  Ary, AryDim1, AryDim2, AryDim3, AryDim4, Descr )
!=======================================================================
!bjj: shouldn't this come after the next subroutine, alphabetically?
   SUBROUTINE CheckIOS ( IOS, Fil, Variable, VarType, TrapErrors )


      ! This routine checks the I/O status and prints either an end-of-file or
      ! an invalid-input message, and then aborts the program.


      ! Argument declarations.

   INTEGER, INTENT(IN)          :: IOS                                          ! I/O status.
   INTEGER, INTENT(IN)          :: VarType                                      ! Type of variable.

   LOGICAL, INTENT(IN), OPTIONAL:: TrapErrors                                   ! Determines if the program should abort or return to calling function
   LOGICAL                      :: TrapThisError                                ! The local version of TrapErrors

   CHARACTER(*), INTENT(IN)     :: Fil                                          ! Name of input file.
   CHARACTER(*), INTENT(IN)     :: Variable                                     ! Variable name.


   IF ( PRESENT( TrapErrors ) ) THEN
      TrapThisError = TrapErrors
   ELSE
      TrapThisError = .FALSE.
   END IF


   IF ( IOS < 0 )  THEN

      CALL PremEOF ( TRIM( Fil ), Variable, TrapThisError )

   ELSE IF ( IOS > 0 )  THEN

      SELECTCASE ( VarType )

      CASE ( NumType )
         CALL WrScr1 ( ' Invalid numerical input for file "'//TRIM( Fil )//'".' )
      CASE ( FlgType )
         CALL WrScr1 ( ' Invalid logical input for file "'//TRIM( Fil )//'".' )
      CASE ( StrType )
         CALL WrScr1 ( ' Invalid character input for file "'//TRIM( Fil )//'".' )
      ENDSELECT

      CALL ProgAbort  ( ' The error occurred while trying to read '//TRIM( Variable )//'.', TrapThisError )

   END IF


   RETURN
   END SUBROUTINE CheckIOS ! ( IOS, Fil, Variable, VarType )
!=======================================================================
   SUBROUTINE CheckArgs ( InputFile, ErrStat )


      ! This subroutine is used to check for command-line arguments.


      ! Argument declarations:
   INTEGER, INTENT(OUT),OPTIONAL:: ErrStat                                      ! Error status; if present, program does not abort on error
   CHARACTER(*), INTENT(OUT)    :: InputFile                                    ! The name of the input file specified on the command line.


      ! Local declarations:

   INTEGER                      :: IArg                                         ! The argument number.
   INTEGER                      :: NumArg                                       ! The number of arguments on the command line.

   LOGICAL                      :: Error                                        ! Flag indicating if there was an error getting an argument.

   CHARACTER(LEN(InputFile))    :: Arg                                          ! A command-line argument.




      ! Find out how many arguments were entered on the command line.

   CALL Get_Arg_Num ( NumArg )


      ! Parse them.

   IF ( NumArg .GT. 0 )  THEN

      DO IArg=1,NumArg

         CALL Get_Arg ( IArg , Arg , Error )

         IF ( Error )  THEN
            CALL ProgAbort ( ' Error getting command-line argument #'//TRIM( Int2LStr( IArg ) )//'.', PRESENT(ErrStat) )
            IF ( PRESENT(ErrStat) ) THEN
               ErrStat = 1
               RETURN
            END IF
         END IF

         IF ( Arg(1:1) == SwChar )  THEN

            CALL WrScr1   ( ' Syntax is:' )
            CALL WrScr1   ( '    '//TRIM( ProgName )//' ['//SwChar//'h] [<InputFile>]' )
            CALL WrScr1   ( ' where:' )
            CALL WrScr1   ( '    '//SwChar//'h generates this help message.' )
            CALL WrScr    ( '    <InputFile> is the name of the primary input file ['//TRIM( InputFile )//'].' )
            CALL WrScr    ( ' ')

            IF ( INDEX( 'Hh?', Arg(2:2)  ) > 0 )  THEN
               IF ( PRESENT(ErrStat) ) THEN
                  ErrStat = -1
                  RETURN
               ELSE
                  CALL ProgExit ( 1 )
               END IF
            ELSE
               CALL ProgAbort ( ' Invalid command-line switch "'//SwChar//TRIM( Arg(2:) )//'".', PRESENT(ErrStat) )
               IF ( PRESENT(ErrStat) ) THEN
                  ErrStat = 1
                  RETURN
               END IF
            END IF

         ELSE
            InputFile = Arg
         END IF

      END DO

   END IF

   IF ( PRESENT( ErrStat ) ) ErrStat = 0

   RETURN
   END SUBROUTINE CheckArgs
!=======================================================================
   SUBROUTINE CloseEcho( )

      ! This subroutine closes the echo file and sets Echo to false

      CLOSE ( UnEc )

      Echo  = .FALSE.

   END SUBROUTINE CloseEcho
!=======================================================================
   SUBROUTINE Conv2UC ( Str )


      ! This routine converts all the text in a string to upper case.


      ! Argument declarations.

   CHARACTER(*), INTENT(INOUT)  :: Str                                          ! The string to be converted to UC.


      ! Local declarations.

   INTEGER                      :: IC                                           ! Character index



   DO IC=1,LEN_TRIM( Str )

      IF ( ( Str(IC:IC) >= 'a' ).AND.( Str(IC:IC) <= 'z' ) )  THEN
         Str(IC:IC) = CHAR( ICHAR( Str(IC:IC) ) - 32 )
      ELSE
         Str(IC:IC) = Str(IC:IC)
      END IF

   END DO ! IC


   RETURN
   END SUBROUTINE Conv2UC !  ( Str )
!=======================================================================
   FUNCTION CountWords ( Line )


      ! This subroutine is used to count the number of "words" in a line of text.
      ! It uses spaces, tabs, commas, semicolons, single quotes, and double quotes ("whitespace")
      !  as word separators.


      ! Function declaration.

   INTEGER                      :: CountWords                                   ! This function.


      ! Argument declarations.

   CHARACTER(*), INTENT(IN)     :: Line                                         ! Count the words in this text string.


      ! Local declarations.

   INTEGER                      :: Ch                                           ! Character position.
   INTEGER                      :: NextWhite                                    ! Position of the next white space.



      ! Let's initialize the number of columns and the character pointer.

   CountWords = 0


      ! Let's make sure we have text on this line.

   IF ( LEN_TRIM( Line ) == 0 )  RETURN


      ! Count words separated by any combination of spaces, tabs, commas,
      ! semicolons, single quotes, and double quotes ("whitespace").

   Ch = 0

   DO

      NextWhite = SCAN( Line(Ch+1:) , ' ,;''"'//Tab )
      Ch        = Ch + NextWhite

      IF ( NextWhite > 1 )  THEN
         CountWords = CountWords + 1
      ELSE IF ( NextWhite == 1 )  THEN
         CYCLE
      ELSE
         EXIT
      END IF

   END DO


   RETURN
   END FUNCTION CountWords ! ( Line )
!=======================================================================
   FUNCTION CurDate( )


      ! This function returns a character string encoded with the date in the form dd-mmm-ccyy.


      ! Function declaration.

   CHARACTER(11)                :: CurDate                                      ! This function


      ! Local declarations.

   CHARACTER(8)                 :: CDate                                        ! String to hold the returned value from the DATE_AND_TIME subroutine call.



   !  Call the system date function.

   CALL DATE_AND_TIME ( CDate )


   !  Parse out the day.

   CurDate(1:3) = CDate(7:8)//'-'


   !  Parse out the month.

   SELECT CASE ( CDate(5:6) )
      CASE ( '01' )
         CurDate(4:6) = 'Jan'
      CASE ( '02' )
         CurDate(4:6) = 'Feb'
      CASE ( '03' )
         CurDate(4:6) = 'Mar'
      CASE ( '04' )
         CurDate(4:6) = 'Apr'
      CASE ( '05' )
         CurDate(4:6) = 'May'
      CASE ( '06' )
         CurDate(4:6) = 'Jun'
      CASE ( '07' )
         CurDate(4:6) = 'Jul'
      CASE ( '08' )
         CurDate(4:6) = 'Aug'
      CASE ( '09' )
         CurDate(4:6) = 'Sep'
      CASE ( '10' )
         CurDate(4:6) = 'Oct'
      CASE ( '11' )
         CurDate(4:6) = 'Nov'
      CASE ( '12' )
         CurDate(4:6) = 'Dec'
   END SELECT


   !  Parse out the year.

   CurDate(7:11) = '-'//CDate(1:4)


   RETURN
   END FUNCTION CurDate ! ()
!=======================================================================
   FUNCTION CurTime( )


      ! This function returns a character string encoded with the time in the form "hh:mm:ss".


      ! Function declaration.

   CHARACTER(8)                 :: CurTime                                      ! This function.


      ! Local declarations.

   CHARACTER(10)                :: CTime                                        ! String to hold the returned value from the DATE_AND_TIME subroutine call.



   CALL DATE_AND_TIME ( TIME=CTime )

   CurTime = CTime(1:2)//':'//CTime(3:4)//':'//CTime(5:6)


   RETURN
   END FUNCTION CurTime ! ()
!=======================================================================
   SUBROUTINE DispNVD


      ! This routine displays the name of the program, its version, and its release date.


      ! Print out program name, version, and date.

   CALL WrScr1 ( ' Running '//TRIM( ProgName )//' '//Trim( ProgVer )//'.' )


   RETURN
   END SUBROUTINE DispNVD
!=======================================================================
   FUNCTION Flt2LStr ( FltNum )


      ! This function converts a floating point number to a left-aligned
      ! string.  It eliminates trailing zeroes and even the decimal
      ! point if it is not a fraction.


      ! Function declaration.

   CHARACTER(15)                :: Flt2LStr                                        ! This function.


      ! Argument declarations.

   REAL(ReKi), INTENT(IN)       :: FltNum                                          ! The floating-point number to convert.


      ! Local declarations.

   INTEGER                      :: IC                                              ! Character index.



      ! Return a 0 if that's what we have.

   IF ( FltNum == 0.0 )  THEN
      Flt2LStr = '0'
      RETURN
   END IF


      ! Write the number into the string using G format and left justify it.

   WRITE (Flt2LStr,'(1PG15.5)')  FltNum

   Flt2LStr = ADJUSTL( Flt2LStr )


      ! Replace trailing zeros and possibly the decimal point with blanks.
      ! Stop trimming once we find the decimal point or a nonzero.

   IF (INDEX( Flt2Lstr, "E" ) > 0 ) RETURN
   IF (INDEX( Flt2Lstr, "e" ) > 0 ) RETURN


   DO IC=LEN_TRIM( Flt2LStr ),1,-1

      IF ( Flt2LStr(IC:IC) == '.' )  THEN
         Flt2LStr(IC:IC) = ' '
         RETURN
      ELSE IF ( Flt2LStr(IC:IC) /= '0' )  THEN
         RETURN
      END IF

      Flt2LStr(IC:IC) = ' '

   END DO ! IC


   RETURN
   END FUNCTION Flt2LStr !  ( FltNum )
!=======================================================================
   SUBROUTINE GetNewUnit ( UnIn )

      ! This routine returns a unit number not currently in use.
      

      ! Argument declarations.

   INTEGER, INTENT(OUT)         :: UnIn                                         ! Logical unit for the file.


      ! Local declarations.

   INTEGER, SAVE                :: Un   = 10                                    ! Unit number; saved between calls (and a GLOBAL) variable
   LOGICAL                      :: Opened                                       ! Flag indicating whether or not a file is opened.
   


      ! See if unit is connected to an open file. Check the next largest number until it is not opened.

   DO

      INQUIRE ( UNIT=Un , OPENED=Opened )
      
      IF ( .NOT. Opened )  EXIT
      Un = Un + 1
      
!      IF ( Un > 99 ) Un = 10                                                     !some systems don't like unit numbers > 99, but we also don't want an infinite loop

   END DO

   UnIn = Un
   
   RETURN   
   END SUBROUTINE GetNewUnit !  ( UnIn )
!=======================================================================

   SUBROUTINE GetPath ( GivenFil, PathName )


      ! Let's parse the path name from the name of the given file.
      ! We'll count everything before (and including) the last "\" or "/".


      ! Argument declarations.

   CHARACTER(*), INTENT(IN)     :: GivenFil                                     ! The name of the given file.
   CHARACTER(*), INTENT(OUT)    :: PathName                                     ! The path name of the given file.


      ! Local declarations.

   INTEGER                      :: I                                            ! DO index for character position.


      ! Look for path separators

   I = INDEX( GivenFil, '\', BACK=.TRUE. )
   I = MAX( I, INDEX( GivenFil, '/', BACK=.TRUE. ) )

   IF ( I == 0 ) THEN
      ! we don't have a path specified, return '.'
      PathName = '.'//PathSep
   ELSE
      PathName = GivenFil(:I)
   END IF


   RETURN
   END SUBROUTINE GetPath ! ( GivenFil, PathName )
!=======================================================================
   SUBROUTINE GetRoot ( GivenFil, RootName )


      ! Let's parse the root file name from the name of the given file.
      ! We'll count everything after the last period as the extension.


      ! Argument declarations.

   CHARACTER(*), INTENT(IN)     :: GivenFil                                     ! The name of the given file.
   CHARACTER(*), INTENT(OUT)    :: RootName                                     ! The parsed root name of the given file.


      ! Local declarations.

   INTEGER                      :: I                                            ! DO index for character position.



      ! Deal with a couple of special cases.

   IF ( ( TRIM( GivenFil ) == "." ) .OR. (  TRIM( GivenFil ) == ".." ) )  THEN
      RootName = TRIM( GivenFil )
      RETURN
   END IF


      ! More-normal cases.

   DO I=LEN_TRIM( GivenFil ),1,-1


      IF ( GivenFil(I:I) == '.' )  THEN


         IF ( I < LEN_TRIM( GivenFil ) ) THEN                   ! Make sure the index I is okay
            IF ( INDEX( '\/', GivenFil(I+1:I+1)) == 0 ) THEN    ! Make sure we don't have the RootName in a different directory
               RootName = GivenFil(:I-1)
            ELSE
               RootName = GivenFil                              ! This does not have a file extension
            END IF
         ELSE
            IF ( I == 1 ) THEN
               RootName = ''
            ELSE
               RootName = GivenFil(:I-1)
            END IF
         END IF

         RETURN

      END IF
   END DO ! I

   RootName =  GivenFil


   RETURN
   END SUBROUTINE GetRoot ! ( GivenFil, RootName )
!=======================================================================
   SUBROUTINE GetTokens ( Line, NumTok, Tokens, Error )


      ! This routine will parse Line for NumTok "tokens" and return them in the Tokens array.
      ! THis routine differs from GetWords() in that it uses only spaces as token separators.


      ! Argument declarations.

   INTEGER, INTENT(IN)          :: NumTok                                       ! The number of "words" to look for.

   LOGICAL, INTENT(OUT)         :: Error                                        ! Error flag to indicate an insuffient number of tokens were found.

   CHARACTER(*), INTENT(INOUT)  :: Line                                         ! The string to search.
   CHARACTER(*), INTENT(OUT)    :: Tokens  (NumTok)                             ! The tokens that were found.


      ! Local declarations.

   INTEGER                      :: IT                                           ! Token index.
   INTEGER                      :: NextBlank                                    ! The location of the next blank character.



   NextBlank = 0

   DO IT=1,NumTok

      Line      = ADJUSTL( Line(NextBlank+1:) )
      NextBlank = INDEX  ( Line , ' ' )

      IF ( NextBlank == 0 )  THEN
        Error = .TRUE.
        RETURN
      END IF

      Tokens(IT) = Line(1:NextBlank-1)

   END DO ! IT

   Error = .FALSE.


   RETURN
   END SUBROUTINE GetTokens ! ( Line, NumTok, Tokens, Error )
!=======================================================================
   SUBROUTINE GetWords ( Line, Words, NumWords )


      ! This subroutine is used to get NumWords "words" from a line of text.


      ! Argument declarations.

   INTEGER, INTENT(IN)          :: NumWords                                     ! The number of words to look for.

   CHARACTER(*), INTENT(IN)     :: Line                                         ! The string to search.
   CHARACTER(*), INTENT(OUT)    :: Words(NumWords)                              ! The array of found words.


      ! Local declarations.

   INTEGER                      :: Ch                                           ! Character position within the string.
   INTEGER                      :: IW                                           ! Word index.
   INTEGER                      :: NextWhite                                    ! The location of the next whitespace in the string.



      ! Let's prefill the array with blanks.

   DO IW=1,NumWords
      Words(IW) = ' '
   END DO ! IW


      ! Let's make sure we have text on this line.

   IF ( LEN_TRIM( Line ) == 0 )  RETURN


      ! Parse words separated by any combination of spaces, tabs, commas,
      ! semicolons, single quotes, and double quotes ("whitespace").

   Ch = 0
   IW = 0

   DO

      NextWhite = SCAN( Line(Ch+1:) , ' ,;''"'//Tab )

      IF ( NextWhite > 1 )  THEN

         IW        = IW + 1
         Words(IW) = Line(Ch+1:Ch+NextWhite-1)

         IF ( IW == NumWords )  EXIT

         Ch = Ch + NextWhite

      ELSE IF ( NextWhite == 1 )  THEN

         Ch = Ch + 1

         CYCLE

      ELSE

         EXIT

      END IF

   END DO


   RETURN
   END SUBROUTINE GetWords ! ( Line, Words, NumWords )
!=======================================================================
   FUNCTION Int2LStr ( Intgr )


      ! This function returns a left-adjusted string representing the passed integer.



   CHARACTER(11)                :: Int2LStr                                     ! This function.


      ! Argument declarations.

   INTEGER, INTENT(IN)          :: Intgr                                        ! The integer to convert to a left-justified string.



   WRITE (Int2LStr,'(I11)')  Intgr

   Int2Lstr = ADJUSTL( Int2LStr )


   RETURN
   END FUNCTION Int2LStr ! ( Intgr )
!=======================================================================
   SUBROUTINE NameOFile ( InArg, OutExten, OutFile, ErrStat )


      ! Get the name of the input file from the InArgth command-line argument.
      ! Remove the extension if there is one, and append OutExten to the end.


      ! Argument declarations.

   INTEGER, INTENT(OUT), OPTIONAL :: ErrStat                                     ! Error status; if present, program does not abort on error
   INTEGER, INTENT(IN)          :: InArg                                        ! The number of the command-line argument that should hold the input file name.

   CHARACTER(*), INTENT(IN)     :: OutExten                                     ! The requested extension for the output file.
   CHARACTER(*), INTENT(OUT)    :: OutFile                                      ! The name of the output file.


      ! Local declarations.

   CHARACTER(100)               :: InFile                                       ! The name of the input file.
   CHARACTER(100)               :: RootName                                     ! The root name of the input file.



      ! See if the command line has enough arguments.

   IF ( InArg > COMMAND_ARGUMENT_COUNT() )  THEN
      CALL ProgAbort ( 'Insufficient arguments on the command line (at least '//&
                         TRIM( Int2LStr( InArg ) )//' were expected).', PRESENT(ErrStat) )
      IF ( PRESENT( ErrStat ) ) ErrStat = 1
      RETURN
   END IF


      ! Get the root of the input file name (strip off the extension).

   CALL GET_COMMAND_ARGUMENT( InArg, InFile )
   CALL GetRoot ( TRIM( InFile ), RootName )

   OutFile = TRIM( RootName )//'.'//OutExten

   IF ( PRESENT( ErrStat ) ) ErrStat = 0

   RETURN
   END SUBROUTINE NameOFile ! ( InArg, OutExten, OutFile [, ErrStat])
!=======================================================================
   SUBROUTINE NormStop


      ! This routine performs a normal termination of the program.


   CALL WrScr1   ( ' '//TRIM( ProgName )//' terminated normally.' )
   CALL WrScr    ( '' )
   CALL ProgExit ( 0 )


   END SUBROUTINE NormStop
!=======================================================================
   SUBROUTINE OpenBin ( Un, OutFile, RecLen, ErrStat )

      ! This routine opens a binary output file.


      ! Argument declarations.

   INTEGER, INTENT(IN)          :: Un                                           ! Logical unit for the output file.
   INTEGER, INTENT(IN)          :: RecLen                                       ! Length of binary record.
   INTEGER, INTENT(OUT), OPTIONAL :: ErrStat                                     ! Error status; if present, program does not abort on error


   CHARACTER(*), INTENT(IN)     :: OutFile                                      ! Name of the output file.


      ! Local declarations.

   LOGICAL                      :: Error                                        ! Flag to indicate the open failed.


      ! Open output file.  Make sure it worked.

   CALL OpenBinFile ( Un, OutFile, RecLen, Error )

   IF ( Error )  THEN
      CALL ProgAbort ( ' Cannot open file "'//TRIM( OutFile )// &
                       '".  Another program may have locked it for writing.', PRESENT(ErrStat) )
      IF ( PRESENT(ErrStat) ) ErrStat = 1
      RETURN
   ELSE
      IF ( PRESENT(ErrStat) ) ErrStat = 0
   ENDIF


   RETURN
   END SUBROUTINE OpenBin ! ( Un, OutFile, RecLen [, ErrStat] )
!=======================================================================
   SUBROUTINE OpenBInpFile ( Un, InFile, ErrStat )


      ! This routine opens a binary input file.

   IMPLICIT                        NONE



      ! Argument declarations.

   INTEGER, INTENT(IN)            :: Un                                          ! Logical unit for the input file.
   INTEGER, INTENT(OUT), OPTIONAL :: ErrStat                                     ! Error status; if present, program does not abort on error

   CHARACTER(*), INTENT(IN)       :: InFile                                      ! Name of the input file.


      ! Local declarations.

      ! NOTE: Do not explicitly declare the precision of this variable [as in
      !       LOGICAL(1)] so that the statements using this variable work with
      !       any compiler:
   LOGICAL                      :: Exists                                       ! Flag indicating whether or not a file Exists.
   LOGICAL                      :: Error                                        ! Flag to indicate the open failed.



      ! See if input file Exists.

   INQUIRE ( FILE=TRIM( InFile ) , EXIST=Exists )

   IF ( .NOT. Exists )  THEN
      CALL ProgAbort ( ' The input file, "'//TRIM( InFile )//'", was not found.', PRESENT( ErrStat) )
      IF ( PRESENT( ErrStat) ) ErrStat = -1
      RETURN
   END IF


      ! Open input file.  Make sure it worked.

   CALL OpenBinInpFile ( Un, InFile, Error )

   IF ( Error )  THEN

      CALL ProgAbort ( ' Cannot open file "'//TRIM( InFile )//'".  Another program may have locked it.', PRESENT( ErrStat ) )
      IF ( PRESENT( ErrStat) ) ErrStat = 1
      RETURN

   ELSE

      IF ( PRESENT( ErrStat) ) ErrStat = 0

   END IF


   RETURN
   END SUBROUTINE OpenBInpFile
!=======================================================================
   SUBROUTINE OpenEcho ( Un, OutFile, ErrStat )


      ! This routine opens a formatted output file for the echo file.


      ! Argument declarations.

   INTEGER, INTENT(OUT), OPTIONAL :: ErrStat                                   ! Error status; if present, program does not abort on error
   INTEGER, INTENT(IN)            :: Un                                        ! Logical unit for the input file.

   CHARACTER(*), INTENT(IN)       :: OutFile                                   ! Name of the input file.




   UnEc = Un

   IF ( PRESENT(ErrStat) ) THEN

      CALL OpenFOutFile( UnEc, OutFile, ErrStat )

   ELSE

      CALL OpenFOutFile( UnEc, OutFile )

   ENDIF

   Echo = .TRUE.

   RETURN
   END SUBROUTINE OpenEcho ! ( Un, OutFile [, ErrStat]  )
!=======================================================================
   SUBROUTINE OpenFInpFile ( Un, InFile, ErrStat )


      ! This routine opens a formatted input file.


      ! Argument declarations.

   INTEGER, INTENT(IN)          :: Un                                           ! Logical unit for the input file.
   INTEGER, INTENT(OUT),OPTIONAL:: ErrStat                                      ! Error status; if present, program does not abort on error
   CHARACTER(*), INTENT(IN)     :: InFile                                       ! Name of the input file.


      ! Local declarations.

   INTEGER                      :: IOS                                          ! I/O status of OPEN.

   LOGICAL                      :: Exists                                       ! Flag indicating whether or not a file Exists.



      ! See if input file Exists.

   INQUIRE ( FILE=TRIM( InFile ) , EXIST=Exists )

   IF ( .NOT. Exists )  THEN
      CALL ProgAbort ( ' The input file, "'//TRIM( InFile )//'", was not found.', PRESENT(ErrStat) )
      IF (PRESENT(ErrStat)) ErrStat = -1
      RETURN
   END IF


      ! Open input file.  Make sure it worked.

   OPEN( Un, FILE=TRIM( InFile ), STATUS='OLD', FORM='FORMATTED', IOSTAT=IOS, ACTION='READ' )

   IF ( IOS /= 0 )  THEN
      CALL ProgAbort ( ' Cannot open file "'//TRIM( InFile ) &
                      //'".  Another program like MS Excel may have locked it for writing.',PRESENT(ErrStat)  )
      IF (PRESENT(ErrStat)) ErrStat = 1
      RETURN
   ELSE
      IF (PRESENT(ErrStat)) ErrStat = 0
   END IF


   RETURN
   END SUBROUTINE OpenFInpFile ! ( Un, InFile [, ErrStat] )
!=======================================================================
   SUBROUTINE OpenFOutFile ( Un, OutFile, ErrStat )


      ! This routine opens a formatted output file.


      ! Argument declarations.

   INTEGER, INTENT(OUT), OPTIONAL :: ErrStat                                    ! Error status; if present, program does not abort on error
   INTEGER, INTENT(IN)          :: Un                                           ! Logical unit for the output file.

   CHARACTER(*), INTENT(IN)     :: OutFile                                      ! Name of the output file.


      ! Local declarations.

   INTEGER                      :: IOS                                          ! I/O status of OPEN.



      ! Open output file.  Make sure it worked.

   OPEN( Un, FILE=TRIM( OutFile ), STATUS='UNKNOWN', FORM='FORMATTED', IOSTAT=IOS, ACTION="WRITE" )



   IF ( PRESENT(ErrStat) ) ErrStat = IOS

   IF ( IOS /= 0 )  CALL ProgAbort( ' Cannot open file "'//TRIM( OutFile )// &
                                    '".  Another program like MS Excel may have locked it for writing.', PRESENT(ErrStat) )


   RETURN
   END SUBROUTINE OpenFOutFile ! ( Un, OutFile [, ErrStat] )
!=======================================================================
   SUBROUTINE OpenFUnkFile ( Un, OutFile, FailAbt, Failed, Exists, ErrStat )


      ! This routine opens a formatted output file and returns a flag telling if it already existed.


      ! Argument declarations.

   INTEGER, INTENT(OUT), OPTIONAL :: ErrStat                                    ! Error status; if present, program does not abort on error
   INTEGER, INTENT(IN)          :: Un                                           ! Logical unit for the output file.

   LOGICAL, INTENT(OUT)         :: Exists                                       ! Flag that indicates if the file already existedo.
   LOGICAL, INTENT(IN)          :: FailAbt                                      ! Flag that tells this routine to abort if the open fails.
   LOGICAL, INTENT(OUT)         :: Failed                                       ! Flag that indicates if the open failed.

   CHARACTER(*), INTENT(IN)     :: OutFile                                      ! Name of the output file.


      ! Local declarations.

   INTEGER                      :: IOS                                          ! I/O status of OPEN.



      ! Check to see if the file already exists.

   INQUIRE ( FILE=TRIM( OutFile ) , EXIST=Exists )

!bjj: should we be checking something here?


      ! Open output file.  Make sure it worked.

   OPEN( Un, FILE=TRIM( OutFile ), STATUS='UNKNOWN', FORM='FORMATTED', IOSTAT=IOS )

   IF ( PRESENT(ErrStat) ) ErrStat = IOS

   IF ( IOS /= 0 )  THEN
      Failed = .TRUE.
      IF ( FailAbt )  CALL ProgAbort ( ' Cannot open file "'//TRIM( OutFile ) &
                                 //'".  Another program like MS Excel may have locked it for writing.', PRESENT(ErrStat) )
   ELSE
      Failed = .FALSE.
   END IF


   RETURN
   END SUBROUTINE OpenFUnkFile ! ( Un, OutFile, FailAbt, Failed, Exists [,ErrStat] )
!=======================================================================
   SUBROUTINE OpenUInfile ( Un, InFile, ErrStat )


      !  This routine opens an unformatted input file.


      ! Argument declarations.

   INTEGER, INTENT(OUT), OPTIONAL :: ErrStat                                    ! Error status; if present, program does not abort on error
   INTEGER, INTENT(IN)         ::  Un                                           ! Logical unit for the input file

   CHARACTER(*), INTENT(IN)    ::  InFile                                       ! Name of the input file


      ! Local declarations.

   INTEGER                     ::  IOS                                          ! Returned input/output status.

   LOGICAL                      :: Exists                                       ! Flag indicating whether or not a file Exists.



      ! See if input file Exists.

   INQUIRE ( FILE=TRIM( InFile ) , EXIST=Exists )

   IF ( .NOT. Exists )  THEN
      CALL ProgAbort ( ' The input file, "'//TRIM( InFile )//'", was not found.', PRESENT(ErrStat) )
      IF ( PRESENT(ErrStat) ) ErrStat = -1
      RETURN
   END IF


      ! Open the file.

   OPEN ( Un, FILE=TRIM( InFile ), STATUS='UNKNOWN', FORM=UnfForm, ACCESS='SEQUENTIAL', IOSTAT=IOS, ACTION='READ' )

   IF ( PRESENT(ErrStat) ) ErrStat = IOS

   IF ( IOS /= 0 )  CALL ProgAbort( ' Cannot open file "'//TRIM( InFile )// &
                                    '".  Another program may have locked it for writing.', PRESENT(ErrStat) )



   RETURN
   END SUBROUTINE OpenUInfile ! ( Un, InFile [,ErrStat] )
!=======================================================================
SUBROUTINE OpenUInBEFile( Un, InFile, RecLen, ErrStat )

      !  This routine opens an unformatted input file of RecLen-byte data records
      !  stored in Big Endian format.


      ! Argument declarations.

   INTEGER, INTENT(IN)           ::  Un                                         ! Logical unit for the input file
   CHARACTER(*), INTENT(IN)      ::  InFile                                     ! Name of the input file
   INTEGER, INTENT(IN)           ::  RecLen                                     ! The input file's record length in bytes
   INTEGER, INTENT(OUT),OPTIONAL ::  ErrStat                                    ! Error status; if present, program does not abort on error



      ! Local declarations.

   LOGICAL                       :: Exists                                       ! Flag to indicate if a file exists
   LOGICAL                       :: Error                                        ! Flag to indicate the open failed



      ! See if input file Exists.

   INQUIRE ( FILE=TRIM( InFile ) , EXIST=Exists )

   IF ( .NOT. Exists )  THEN
      CALL ProgAbort ( ' The input file, "'//TRIM( InFile )//'", was not found.', PRESENT(ErrStat) )
      IF ( PRESENT(ErrStat) ) ErrStat = -1
      RETURN
   END IF


      ! Open the file.

   CALL OpenUnfInpBEFile ( Un, InFile, RecLen, Error )

   IF ( Error )  THEN

      CALL ProgAbort ( ' Cannot open file "'//TRIM( InFile )//'".  Another program may have locked it.', PRESENT( ErrStat ) )
      IF ( PRESENT( ErrStat) ) ErrStat = 1
      RETURN

   ELSE

      IF ( PRESENT( ErrStat) ) ErrStat = 0

   END IF


   RETURN

END SUBROUTINE OpenUInBEFile !( Un, InFile, RecLen [, ErrStat] )
!=======================================================================
   SUBROUTINE OpenUOutfile ( Un, OutFile, ErrStat )


      !  This routine opens an unformatted output file.


      ! Argument declarations.

   INTEGER, INTENT(IN)            ::  Un                                        ! Logical unit for the output file
   INTEGER, INTENT(OUT), OPTIONAL :: ErrStat                                    ! Error status; if present, program does not abort on error

   CHARACTER(*), INTENT(IN)       ::  OutFile                                   ! Name of the output file


      ! Local declarations.

   INTEGER                        ::  IOS                                       ! Returned input/output status.



      ! Open the file.

   OPEN ( Un, FILE=TRIM( OutFile ), STATUS='UNKNOWN', FORM=UnfForm, ACCESS='SEQUENTIAL', IOSTAT=IOS, ACTION='WRITE' )


   IF ( PRESENT( ErrStat ) )   ErrStat = IOS

   IF ( IOS /= 0 )  CALL ProgAbort( ' Cannot open file "'//TRIM( OutFile )// &
                                    '".  Another program may have locked it for writing.', PRESENT( ErrStat ) )


   RETURN
   END SUBROUTINE OpenUOutfile ! ( Un, InFile [,ErrStat] )
!=======================================================================
   FUNCTION PathIsRelative ( GivenFil )


      ! Let's determine in the given file name is absolute or relative.
      !
      ! We'll consider an absolute path one that satisfies one of the 
      ! following four criteria:
      !     1) It contains ":/"
      !     2) It contains ":\"
      !     3) It starts with "/"
      !     4) It starts with "\"
      ! All others are considered relative.



      ! Argument declarations.

   CHARACTER(*), INTENT(IN)     :: GivenFil                                            ! The name of the given file.
   LOGICAL                      :: PathIsRelative                                      ! The function return value


      ! Determine if file name begins with an absolute path name or if it is relative

   PathIsRelative = .FALSE.

   IF ( ( INDEX( GivenFil, ':/') == 0 ) .AND. ( INDEX( GivenFil, ':\') == 0 ) ) THEN   ! No drive is specified (by ':\' or ':/')
   
      IF ( INDEX( '/\', GivenFil(1:1) ) == 0 ) THEN                                    ! The file name doesn't start with '\' or '/'

         PathIsRelative = .TRUE.

      END IF

   END IF

   RETURN
   END FUNCTION PathIsRelative ! ( GivenFil )
!=======================================================================
   SUBROUTINE PremEOF ( Fil , Variable, TrapErrors )


      ! This routine prints out an EOF message and aborts the program.


      ! Argument declarations.

   LOGICAL, INTENT(IN), OPTIONAL:: TrapErrors                                   ! Determines if the program should abort or return to calling function
   LOGICAL                      :: TrapThisError                                ! The local version of TrapErrors

   CHARACTER(*), INTENT(IN)     :: Fil                                          ! The name of the file that ran out of data.
   CHARACTER(*), INTENT(IN)     :: Variable                                     ! The name of the variable we were trying to read at the time.


   IF ( PRESENT( TrapErrors ) ) THEN
      TrapThisError = TrapErrors
   ELSE
      TrapThisError = .FALSE.
   END IF

   CALL WrScr1 ( ' Premature EOF for file "'//TRIM( Fil )//'".' )

   CALL ProgAbort  ( ' The error occurred while trying to read '//TRIM( Variable )//'.', TrapThisError )


   RETURN
   END SUBROUTINE PremEOF ! ( Fil , Variable [, TrapErrors] )
!=======================================================================
   SUBROUTINE ProgAbort ( Message, TrapErrors )


      ! This routine outputs fatal error messages and stops the program.


      ! Argument declarations.

   LOGICAL, INTENT(IN), OPTIONAL:: TrapErrors                                   ! Determines if the program should abort or return to calling function
   CHARACTER(*), INTENT(IN)     :: Message                                      ! Error message.



   IF ( Beep )  CALL UsrAlarm

   CALL WrScr    ( Message )
   IF ( PRESENT(TrapErrors) )  THEN
      IF ( TrapErrors ) RETURN
   END IF

   IF ( LEN_TRIM(ProgName) > 0 ) THEN
      CALL WrScr1   ( ' Aborting '//TRIM( ProgName )//'.' )
   ELSE
      CALL WrScr1   ( ' Aborting program.' )
   END IF
   CALL WrScr    ( ' ' )
   CALL ProgExit ( 1 )


   END SUBROUTINE ProgAbort ! ( Message [, TrapErrors] )
!=======================================================================
   SUBROUTINE ProgWarn ( Message )


      ! This routine outputs non-fatal warning messages and returns to the calling routine.


      ! Argument declarations.

   CHARACTER(*), INTENT(IN)     :: Message                                      ! Warning message.



   IF ( Beep )  CALL UsrAlarm
   CALL WrScr ( ' WARNING:  '//Message )


   RETURN
   END SUBROUTINE ProgWarn ! ( Message )
!=======================================================================
   SUBROUTINE ReadCAry ( UnIn, Fil, CharAry, AryLen, AryName, AryDescr, ErrStat )


      ! This routine reads a AryLen values into a character array from the next AryLen lines of the input file.


      ! Argument declarations:

   INTEGER, INTENT(IN)          :: AryLen                                          ! Length of the array.
   INTEGER, INTENT(IN)          :: UnIn                                            ! I/O unit for input file.
   INTEGER, INTENT(OUT),OPTIONAL:: ErrStat                                         ! Error status; if present, program does not abort on error

   CHARACTER(*), INTENT(OUT)    :: CharAry(AryLen)                                 ! Real variable being read.
   CHARACTER(*), INTENT(IN)     :: AryDescr                                        ! Text string describing the variable.
   CHARACTER(*), INTENT(IN)     :: AryName                                         ! Text string containing the variable name.
   CHARACTER(*), INTENT(IN)     :: Fil                                             ! Name of the input file.


      ! Local declarations:

   INTEGER                      :: Ind                                             ! Index into the string array.  Assumed to be one digit.
   INTEGER                      :: IOS                                             ! I/O status returned from the read statement.

   CHARACTER(44)                :: Frmt = "(15X,A,T30,' - ',A,/,2X,100('""',A,'""',:,1X))"    ! Output format for string parameters.



   READ (UnIn,*,IOSTAT=IOS)  ( CharAry(Ind), Ind=1,AryLen )

   CALL CheckIOS ( IOS, Fil, TRIM( AryName ), StrType, PRESENT(ErrStat) )

   IF (PRESENT(ErrStat) ) THEN
      ErrStat = IOS
      IF (IOS /= 0) RETURN
   ENDIF

   IF ( Echo )  THEN
      WRITE (UnEc,Frmt)  TRIM( AryName ), AryDescr, ( TRIM( CharAry(Ind) ), Ind=1,AryLen )
   END IF


   RETURN
   END SUBROUTINE ReadCAry ! ( UnIn, Fil, CharAry, AryLen, AryName, AryDescr [, ErrStat] )
!=======================================================================
   SUBROUTINE ReadCAryLines ( UnIn, Fil, CharAry, AryLen, AryName, AryDescr, ErrStat )


      ! This routine reads a AryLen values into a real array from the next AryLen lines of the input file.


      ! Argument declarations:

   INTEGER, INTENT(IN)          :: AryLen                                          ! Length of the array.
   INTEGER, INTENT(IN)          :: UnIn                                            ! I/O unit for input file.
   INTEGER, INTENT(OUT),OPTIONAL:: ErrStat                                         ! Error status; if present, program does not abort on error

   CHARACTER(*), INTENT(OUT)    :: CharAry(AryLen)                                 ! Char variable being read.

   CHARACTER(*), INTENT(IN)     :: Fil                                             ! Name of the input file.
   CHARACTER(*), INTENT(IN)     :: AryDescr                                        ! Text string describing the variable.
   CHARACTER(*), INTENT(IN)     :: AryName                                         ! Text string containing the variable name.


      ! Local declarations:

   INTEGER                      :: Ind                                             ! Index into the real array.  Assumed to be one digit.
   INTEGER                      :: IOS                                             ! I/O status returned from the read statement.

   CHARACTER(35)                :: Frmt = "( 15X, A, T30, ' - ', A, /, 2X, A )"    ! Output format for string parameters.


   IF ( PRESENT(ErrStat) ) ErrStat = 0

   DO Ind=1,AryLen
      READ (UnIn,*,IOSTAT=IOS)  CharAry(Ind)

      CALL CheckIOS ( IOS, Fil, TRIM( AryName )//'('//TRIM( Int2LStr( Ind ) )//')', StrType, PRESENT(ErrStat) )

      IF (IOS /= 0) THEN
         IF ( PRESENT(ErrStat) ) ErrStat = IOS
         RETURN
      ENDIF

      IF ( Echo )  THEN
         WRITE (UnEc,Frmt)  TRIM( AryName )//'('//TRIM( Int2LStr( Ind ) )//')', AryDescr, TRIM(CharAry(Ind))
      END IF
   END DO

   RETURN
   END SUBROUTINE ReadCAryLines ! ( UnIn, Fil, RealAry, AryLen, AryName, AryDescr [, ErrStat] )
!=======================================================================
   SUBROUTINE ReadCom ( UnIn, Fil, ComName, ErrStat )

      ! This routine reads a comment from the next line of the input file.


      ! Argument declarations:

   INTEGER, INTENT(IN)          :: UnIn                                            ! I/O unit for input file.
   INTEGER, INTENT(OUT),OPTIONAL:: ErrStat                                         ! Error status; if present, program does not abort on error

   CHARACTER(*), INTENT(IN)   :: Fil                                             ! Name of the input file.
   CHARACTER(*), INTENT(IN)   :: ComName                                         ! Text string containing the comment name.


      ! Local declarations:

   INTEGER                      :: IOS                                             ! I/O status returned from the read statement.

   CHARACTER(200)               :: Comment                                         ! Text string containing the comment.



   READ (UnIn,'(A)',IOSTAT=IOS)  Comment

   CALL CheckIOS ( IOS, Fil, ComName, StrType, PRESENT(ErrStat) )

   IF ( PRESENT(ErrStat) ) THEN
      ErrStat = IOS
      IF (IOS /= 0) RETURN
   ENDIF

   IF ( Echo )  THEN
      WRITE (UnEc,'(A)')  Comment
   END IF


   RETURN
   END SUBROUTINE ReadCom ! ( UnIn, Fil, ComName [, ErrStat] )
!=======================================================================
   SUBROUTINE ReadCVar ( UnIn, Fil, CharVar, VarName, VarDescr, ErrStat )


      ! This routine reads a single character variable from the next line of the input file.


      ! Argument declarations:

   INTEGER, INTENT(IN)          :: UnIn                                            ! I/O unit for input file.
   INTEGER, INTENT(OUT),OPTIONAL:: ErrStat                                         ! Error status; if present, program does not abort on error

   CHARACTER(*), INTENT(OUT)    :: CharVar                                         ! Integer variable being read.
   CHARACTER(*), INTENT(IN)     :: Fil                                             ! Name of the input file.
   CHARACTER(*), INTENT(IN)     :: VarDescr                                        ! Text string describing the variable.
   CHARACTER(*), INTENT(IN)     :: VarName                                         ! Text string containing the variable name.


      ! Local declarations:

   INTEGER                      :: IOS                                             ! I/O status returned from the read statement.

   CHARACTER(35)                :: Frmt = "( 15X, A, T30, ' - ', A, /, 2X, A )"    ! Output format for string parameters.




   READ (UnIn,*,IOSTAT=IOS)  CharVar

   CALL CheckIOS ( IOS, Fil, VarName, StrType, PRESENT(ErrStat) )

   IF ( PRESENT(ErrStat) ) THEN
      ErrStat = IOS
      IF (IOS /= 0) RETURN
   ENDIF

   IF ( Echo )  THEN
      WRITE (UnEc,Frmt)  VarName, VarDescr, '"'//TRIM( CharVar )//'"'
   END IF


   RETURN
   END SUBROUTINE ReadCVar ! ( UnIn, Fil, CharVar, VarName, VarDescr [, ErrStat] )
!=======================================================================
   SUBROUTINE ReadIAry ( UnIn, Fil, IntAry, AryLen, AryName, AryDescr, ErrStat )


      ! This routine reads a AryLen values into an integer array from the next AryLen lines of the input file.


      ! Argument declarations:

   INTEGER, INTENT(IN)          :: AryLen                                          ! Length of the array.
   INTEGER, INTENT(OUT)         :: IntAry(AryLen)                                  ! Integer array being read.
   INTEGER, INTENT(IN)          :: UnIn                                            ! I/O unit for input file.
   INTEGER, INTENT(OUT),OPTIONAL:: ErrStat                                         ! Error status; if present, program does not abort on error

   CHARACTER(*), INTENT(IN)     :: Fil                                             ! Name of the input file.
   CHARACTER(*), INTENT(IN)     :: AryDescr                                        ! Text string describing the variable.
   CHARACTER(*), INTENT(IN)     :: AryName                                         ! Text string containing the variable name.


      ! Local declarations:

   INTEGER                      :: Ind                                             ! Index into the integer array.  Assumed to be one digit.
   INTEGER                      :: IOS                                             ! I/O status returned from the read statement.

   CHARACTER(38)                :: Frmt = "( 2X, I11, 2X, A, T30, ' - ', A )"      ! Output format for integer array parameters



   READ (UnIn,*,IOSTAT=IOS)  ( IntAry(Ind), Ind=1,AryLen )

   CALL CheckIOS ( IOS, Fil, TRIM( AryName ), NumType, PRESENT(ErrStat) )

   IF ( PRESENT(ErrStat) ) THEN
      ErrStat = IOS
      IF (IOS /= 0) RETURN
   ENDIF

   IF ( Echo )  THEN
      DO Ind=1,AryLen
         WRITE (UnEc,Frmt)  IntAry(Ind), TRIM( AryName ), AryDescr
      END DO ! Ind
   END IF


   RETURN
   END SUBROUTINE ReadIAry ! ( UnIn, Fil, IntAry, AryLen, AryName, AryDescr [, ErrStat])
!=======================================================================
   SUBROUTINE ReadIVar ( UnIn, Fil, IntVar, VarName, VarDescr, ErrStat )


      ! This routine reads a single integer variable from the next line of the input file.


      ! Argument declarations:

   INTEGER, INTENT(OUT)         :: IntVar                                          ! Integer variable being read.
   INTEGER, INTENT(IN)          :: UnIn                                            ! I/O unit for input file.
   INTEGER, INTENT(OUT),OPTIONAL:: ErrStat                                         ! Error status; if present, program does not abort on error

   CHARACTER(*), INTENT(IN)     :: Fil                                             ! Name of the input file.
   CHARACTER(*), INTENT(IN)     :: VarDescr                                        ! Text string describing the variable.
   CHARACTER(*), INTENT(IN)     :: VarName                                         ! Text string containing the variable name.


      ! Local declarations:

   INTEGER                      :: IOS                                             ! I/O status returned from the read statement.

   CHARACTER(33)                :: Frmt = "( 2X, I11, 2X, A, T30, ' - ', A )"      ! Output format for integer parameters.
   CHARACTER(30)                :: Word                                            ! String to hold the first word on the line.



   IF ( PRESENT(ErrStat) ) THEN
      CALL ReadNum ( UnIn, Fil, Word, VarName, ErrStat )
      IF (ErrStat /= 0) RETURN
   ELSE
      CALL ReadNum ( UnIn, Fil, Word, VarName )
   END IF

   READ (Word,*,IOSTAT=IOS)  IntVar

   CALL CheckIOS ( IOS, Fil, VarName, NumType, PRESENT(ErrStat) )

   IF ( PRESENT(ErrStat) ) THEN
      ErrStat = IOS
      IF (IOS /= 0) RETURN
   ENDIF

   IF ( Echo )  THEN
      WRITE (UnEc,Frmt)  IntVar, VarName, VarDescr
   END IF


   RETURN
   END SUBROUTINE ReadIVar ! ( UnIn, Fil, IntVar, VarName, VarDescr [, ErrStat])
!=======================================================================
   SUBROUTINE ReadLAry ( UnIn, Fil, LogAry, AryLen, AryName, AryDescr, ErrStat )


      ! This routine reads a AryLen values into an logical array from the next AryLen lines of the input file.


      ! Argument declarations:

   INTEGER, INTENT(IN)          :: AryLen                                          ! Length of the array.
   INTEGER, INTENT(IN)          :: UnIn                                            ! I/O unit for input file.
   INTEGER, INTENT(OUT),OPTIONAL:: ErrStat                                         ! Error status; if present, program does not abort on error

   LOGICAL, INTENT(OUT)         :: LogAry(AryLen)                                  ! Logical array being read.

   CHARACTER(*), INTENT(IN)     :: Fil                                             ! Name of the input file.
   CHARACTER(*), INTENT(IN)     :: AryDescr                                        ! Text string describing the variable.
   CHARACTER(*), INTENT(IN)     :: AryName                                         ! Text string containing the variable name.


      ! Local declarations:

   INTEGER                      :: Ind                                             ! Index into the integer array.  Assumed to be one digit.
   INTEGER                      :: IOS                                             ! I/O status returned from the read statement.

   CHARACTER(38)                :: Frmt = "( 2X, I11, 2X, A, T30, ' - ', A )"      ! Output format for integer array parameters



   READ (UnIn,*,IOSTAT=IOS)  ( LogAry(Ind), Ind=1,AryLen )

   CALL CheckIOS ( IOS, Fil, TRIM( AryName ), NumType, PRESENT(ErrStat) )

   IF ( PRESENT(ErrStat) ) THEN
      ErrStat = IOS
      IF (IOS /= 0) RETURN
   ENDIF

   IF ( Echo )  THEN
      DO Ind=1,AryLen
         WRITE (UnEc,Frmt)  LogAry(Ind), TRIM( AryName ), AryDescr
      END DO ! Ind
   END IF


   RETURN
   END SUBROUTINE ReadLAry ! ( UnIn, Fil, LogAry, AryLen, AryName, AryDescr [, ErrStat])
!=======================================================================
   SUBROUTINE ReadLVar ( UnIn, Fil, LogVar, VarName, VarDescr, ErrStat )


      ! This routine reads a single logical variable from the next line of the input file.


      ! Argument declarations:

   INTEGER, INTENT(IN)          :: UnIn                                            ! I/O unit for input file.
   INTEGER, INTENT(OUT),OPTIONAL:: ErrStat                                         ! Error status; if present, program does not abort on error

   LOGICAL, INTENT(OUT)         :: LogVar                                          ! Logical variable being read.

   CHARACTER(*), INTENT(IN)     :: Fil                                             ! Name of the input file.
   CHARACTER(*), INTENT(IN)     :: VarDescr                                        ! Text string describing the variable.
   CHARACTER(*), INTENT(IN)     :: VarName                                         ! Text string containing the variable name.


      ! Local declarations:

   INTEGER                      :: IOS                                             ! I/O status returned from the read statement.

   CHARACTER(33)                :: Frmt  = "( 2X, L11, 2X, A, T30, ' - ', A )"     ! Output format for logical parameters.
   CHARACTER( 4)                :: VName                                           ! Temporary holder for the variable name.




   READ (UnIn,*,IOSTAT=IOS)  LogVar

   CALL CheckIOS ( IOS, Fil, VarName, FlgType, PRESENT(ErrStat) )

   IF ( PRESENT(ErrStat) ) THEN
      ErrStat = IOS
      IF (IOS /= 0) RETURN
   ENDIF

   VName = VarName

   CALL Conv2UC ( VName )

   IF ( Echo .AND. ( VName /= 'ECHO' ) )  THEN
      WRITE (UnEc,Frmt)  LogVar, VarName, VarDescr
   END IF


   RETURN
   END SUBROUTINE ReadLVar ! ( UnIn, Fil, LogVar, VarName, VarDescr [, ErrStat])
!=======================================================================
   SUBROUTINE ReadNum ( UnIn, Fil, Word, VarName, ErrStat )


      ! This routine reads a single word from a file and tests to see if it's a pure number (no true or false).


      ! Argument declarations:

   INTEGER, INTENT(IN)            :: UnIn                                          ! I/O unit for input file.
   INTEGER, INTENT(OUT), OPTIONAL :: ErrStat                                       ! Error status; if present, program does not abort on error

   CHARACTER(*), INTENT(IN)       :: Fil                                           ! Name of the input file.
   CHARACTER(*), INTENT(IN)       :: VarName                                       ! Text string containing the variable name.
   CHARACTER(*), INTENT(Out)      :: Word                                          ! Text string containing the first word from the input line.


      ! Local declarations:

   INTEGER                        :: IOS                                           ! I/O status returned from the read statement.


   IF ( PRESENT(ErrStat) ) ErrStat = 0


      ! Read in the first word of the input line.  Check I/O status.

   READ (UnIn,*,IOSTAT=IOS)  Word

   CALL CheckIOS ( IOS, Fil, VarName, NumType, PRESENT(ErrStat) )

   IF ( PRESENT(ErrStat) ) THEN
      ErrStat = IOS
      IF (IOS /= 0) RETURN
   ENDIF


      ! See if the word starts with a T or F.  If so, flag it as an invalid number.

   IF ( INDEX( 'FTft', Word(:1) ) > 0 )  THEN
      CALL WrScr ( '' )
      CALL ProgAbort( ' Invalid numeric input.  "'//TRIM( Word )//'" found when trying to read the number, '// &
                      TRIM( VarName )//'.', PRESENT(ErrStat) )

      IF ( PRESENT(ErrStat) ) ErrStat = 1
   END IF



   RETURN
   END SUBROUTINE ReadNum ! ( UnIn, Fil, Word, VarName [, ErrStat] )
!=======================================================================
   SUBROUTINE ReadOutputList ( UnIn, Fil, CharAry, AryLenRead, AryName, AryDescr, ErrStat )


      ! This routine reads a AryLen values into a real array from the next AryLen lines of the input file.


      ! Argument declarations:

   INTEGER, INTENT(OUT)         :: AryLenRead                                      ! Length of the array that was actually read.
   INTEGER, INTENT(IN)          :: UnIn                                            ! I/O unit for input file.
   INTEGER, INTENT(OUT),OPTIONAL:: ErrStat                                         ! Error status; if present, program does not abort on error

   CHARACTER(*), INTENT(OUT)    :: CharAry(:)                                      ! Character array being read (calling routine dimensions it to max allowable size).

   CHARACTER(*), INTENT(IN)     :: Fil                                             ! Name of the input file.
   CHARACTER(*), INTENT(IN)     :: AryDescr                                        ! Text string describing the variable.
   CHARACTER(*), INTENT(IN)     :: AryName                                         ! Text string containing the variable name.


      ! Local declarations:

   INTEGER                      :: MaxAryLen                                       ! Maximum length of the array being read
   INTEGER                      :: NumWords                                        ! Number of words contained on a line


   CHARACTER(1000)              :: OutLine
   CHARACTER(3)                 :: EndOfFile
   
!   CHARACTER(38)                :: Frmt = "( 2X, ES11.4e2, 2X, A, T30, ' - ', A )" ! Output format for real array parameters

      ! Initialize some values

   IF ( PRESENT(ErrStat) ) ErrStat = 0
   MaxAryLen  = SIZE(CharAry)
   AryLenRead = 0
   
   CharAry = ''   


 
      ! Read in all of the lines containing output parameters and store them in CharAry(:).
      ! The end of this list is specified with the line beginning with END.

   DO

      IF ( PRESENT(ErrStat) ) THEN
         CALL ReadVar ( UnIn, Fil, OutLine, AryName, AryDescr, ErrStat )
         IF ( ErrStat /= 0 ) RETURN
      ELSE
         CALL ReadVar ( UnIn, Fil, OutLine, AryName, AryDescr  )
      END IF      

      EndOfFile = OutLine(1:3)            ! EndOfFile is the 1st 3 characters of OutLine
      CALL Conv2UC( EndOfFile )           ! Convert EndOfFile to upper case
      IF ( EndOfFile == 'END' )  EXIT     ! End of OutList has been reached; therefore, exit this DO

      NumWords = CountWords( OutLine )    ! The number of words in OutLine.

      AryLenRead = AryLenRead + NumWords  ! The total number of output channels read in so far.
      
         ! Check to see if the maximum # allowable in the array has been reached.
         
      IF ( AryLenRead > MaxAryLen )  THEN    
      
         CALL ProgAbort ( ' The maximum number of output channels allowed is ' &
                     //TRIM( Int2LStr(MaxAryLen) )//'.', PRESENT(ErrStat)      )
         ErrStat = 1
         RETURN
         
      ELSE
      
         CALL GetWords ( OutLine, CharAry((AryLenRead - NumWords + 1):AryLenRead), NumWords )      
         
      END IF                     

   END DO        

   RETURN
   END SUBROUTINE ReadOutputList ! ( UnIn, Fil, CharAry, AryLenRead, AryName, AryDescr [, ErrStat] )
!=======================================================================
   SUBROUTINE ReadRAry ( UnIn, Fil, RealAry, AryLen, AryName, AryDescr, ErrStat )


      ! This routine reads a AryLen values into a real array separated by white space (possibly on the same line of the input file).


      ! Argument declarations:

   INTEGER, INTENT(IN)          :: AryLen                                          ! Length of the array.
   INTEGER, INTENT(IN)          :: UnIn                                            ! I/O unit for input file.
   INTEGER, INTENT(OUT),OPTIONAL:: ErrStat                                         ! Error status; if present, program does not abort on error

   REAL(ReKi), INTENT(INOUT)    :: RealAry(AryLen)                                 ! Real array being read.

   CHARACTER(*), INTENT(IN)     :: Fil                                             ! Name of the input file.
   CHARACTER(*), INTENT(IN)     :: AryDescr                                        ! Text string describing the variable.
   CHARACTER(*), INTENT(IN)     :: AryName                                         ! Text string containing the variable name.


      ! Local declarations:

   INTEGER                      :: Ind                                             ! Index into the real array.  Assumed to be one digit.
   INTEGER                      :: IOS                                             ! I/O status returned from the read statement.

   CHARACTER(38)                :: Frmt = "( 2X, ES11.4e2, 2X, A, T30, ' - ', A )" ! Output format for real array parameters



   READ (UnIn,*,IOSTAT=IOS)  ( RealAry(Ind), Ind=1,AryLen )

   CALL CheckIOS ( IOS, Fil, TRIM( AryName ), NumType, PRESENT(ErrStat) )

   IF ( PRESENT(ErrStat) ) THEN
      ErrStat = IOS
      IF (IOS /= 0) RETURN
   ENDIF

   IF ( Echo )  THEN
      DO Ind=1,AryLen
         WRITE (UnEc,Frmt)  RealAry(Ind), TRIM( AryName ), AryDescr
      END DO ! Ind
   END IF


   RETURN
   END SUBROUTINE ReadRAry ! ( UnIn, Fil, RealAry, AryLen, AryName, AryDescr [, ErrStat] )
!=======================================================================
   SUBROUTINE ReadRAryLines ( UnIn, Fil, RealAry, AryLen, AryName, AryDescr, ErrStat )


      ! This routine reads a AryLen values into a real array from the next AryLen lines of the input file.


      ! Argument declarations:

   INTEGER, INTENT(IN)          :: AryLen                                          ! Length of the array.
   INTEGER, INTENT(IN)          :: UnIn                                            ! I/O unit for input file.
   INTEGER, INTENT(OUT),OPTIONAL:: ErrStat                                         ! Error status; if present, program does not abort on error

   REAL(ReKi), INTENT(OUT)      :: RealAry(AryLen)                                 ! Real array being read.

   CHARACTER(*), INTENT(IN)     :: Fil                                             ! Name of the input file.
   CHARACTER(*), INTENT(IN)     :: AryDescr                                        ! Text string describing the variable.
   CHARACTER(*), INTENT(IN)     :: AryName                                         ! Text string containing the variable name.


      ! Local declarations:

   INTEGER                      :: Ind                                             ! Index into the real array.  Assumed to be one digit.
   INTEGER                      :: IOS                                             ! I/O status returned from the read statement.

   CHARACTER(38)                :: Frmt = "( 2X, ES11.4e2, 2X, A, T30, ' - ', A )" ! Output format for real array parameters


   IF ( PRESENT(ErrStat) ) ErrStat = 0

   DO Ind=1,AryLen
      READ (UnIn,*,IOSTAT=IOS)  RealAry(Ind)

      CALL CheckIOS ( IOS, Fil, TRIM( AryName )//'('//TRIM( Int2LStr( Ind ) )//')', NumType, PRESENT(ErrStat) )

      IF (IOS /= 0) THEN
         IF ( PRESENT(ErrStat) ) ErrStat = IOS
         RETURN
      ENDIF

      IF ( Echo )  THEN
         WRITE (UnEc,Frmt)  RealAry(Ind), TRIM( AryName )//'('//TRIM( Int2LStr( Ind ) )//')', AryDescr
      END IF
   END DO

   RETURN
   END SUBROUTINE ReadRAryLines ! ( UnIn, Fil, RealAry, AryLen, AryName, AryDescr [, ErrStat] )
!=======================================================================
  SUBROUTINE ReadRVar ( UnIn, Fil, RealVar, VarName, VarDescr, ErrStat )


      ! This routine reads a single real variable from the next line of the input file.


      ! Argument declarations:

   REAL(ReKi), INTENT(OUT)      :: RealVar                                         ! Real variable being read.
   INTEGER, INTENT(OUT),OPTIONAL:: ErrStat                                         ! Error status; if present, program does not abort on error

   INTEGER, INTENT(IN)          :: UnIn                                            ! I/O unit for input file.

   CHARACTER( *), INTENT(IN)    :: Fil                                             ! Name of the input file.
   CHARACTER( *), INTENT(IN)    :: VarDescr                                        ! Text string describing the variable.
   CHARACTER( *), INTENT(IN)    :: VarName                                         ! Text string containing the variable name.


      ! Local declarations:

   INTEGER                      :: IOS                                             ! I/O status returned from the read statement.

   CHARACTER(38)                :: Frmt = "( 2X, ES11.4e2, 2X, A, T30, ' - ', A )" ! Output format for real parameters
   CHARACTER(30)                :: Word                                            ! String to hold the first word on the line.



   IF ( PRESENT(ErrStat) ) THEN
      CALL ReadNum ( UnIn, Fil, Word, VarName, ErrStat )
      IF (ErrStat /= 0) RETURN
   ELSE
      CALL ReadNum ( UnIn, Fil, Word, VarName )
   END IF

   READ (Word,*,IOSTAT=IOS)  RealVar

   CALL CheckIOS ( IOS, Fil, VarName, NumType, PRESENT(ErrStat) )

   IF ( PRESENT(ErrStat) ) THEN
      ErrStat = IOS
      IF (IOS /= 0) RETURN
   END IF

   IF ( Echo )  THEN
      WRITE (UnEc,Frmt)  RealVar, VarName, VarDescr
   END IF


   RETURN
   END SUBROUTINE ReadRVar ! ( UnIn, Fil, RealVar, VarName, VarDescr [, ErrStat] )
!=======================================================================
   SUBROUTINE ReadStr ( UnIn, Fil, CharVar, VarName, VarDescr, ErrStat )


      ! This routine reads a string from the next line of the input file.


      ! Argument declarations:

   INTEGER, INTENT(IN)          :: UnIn                                            ! I/O unit for input file.
   INTEGER, INTENT(OUT),OPTIONAL:: ErrStat                                         ! Error status; if present, program does not abort on error

   CHARACTER(*), INTENT(OUT)    :: CharVar                                         ! Integer variable being read.
   CHARACTER(*), INTENT(IN)     :: Fil                                             ! Name of the input file.
   CHARACTER(*), INTENT(IN)     :: VarDescr                                        ! Text string describing the variable.
   CHARACTER(*), INTENT(IN)     :: VarName                                         ! Text string containing the variable name.


      ! Local declarations:

   INTEGER                      :: IOS                                             ! I/O status returned from the read statement.

   CHARACTER(35)                :: Frmt = "( 15X, A, T30, ' - ', A, /, 2X, A )"    ! Output format for string parameters.




   READ (UnIn,'(A)',IOSTAT=IOS)  CharVar

   CALL CheckIOS ( IOS, Fil, VarName, StrType )

   IF ( PRESENT(ErrStat) ) THEN
      ErrStat = IOS
      IF (IOS /= 0) RETURN
   END IF

   IF ( Echo )  THEN
      WRITE (UnEc,Frmt)  VarName, VarDescr, '"'//TRIM( CharVar )//'"'
   END IF


   RETURN
   END SUBROUTINE ReadStr ! ( UnIn, Fil, CharVar, VarName, VarDescr [, ErrStat] )
!=======================================================================
   SUBROUTINE WaitTime ( WaitSecs )


      ! This routine writes out a prompt to the screen without
      ! following it with a new line, though a new line precedes it.


   IMPLICIT NONE


      ! Argument declarations:

   REAL(ReKi), INTENT(IN)       :: WaitSecs                                        ! The number of seconds to wait.


      ! Local declarations:

   REAL(ReKi)                   :: EndCounts                                       ! The number of counts when wait time is over.

   INTEGER                      :: Counts                                          ! Current number of counts on the system clock.
   INTEGER                      :: CountMax                                        ! Maximum number of counts possible on the system clock.
   INTEGER                      :: CountRate                                       ! Number of counts per second on the system clock.



   CALL SYSTEM_CLOCK ( Counts, CountRate, CountMax )
   EndCounts = Counts + INT( WaitSecs*CountRate )

   DO
      CALL SYSTEM_CLOCK ( Counts, CountRate, CountMax )
      IF ( Counts > EndCounts )  EXIT
   END DO


   RETURN
   END SUBROUTINE WaitTime ! ( Seconds )
!=======================================================================
   SUBROUTINE WrPr ( Str )


      ! This routine writes out a prompt to the screen without
      ! following it with a new line, though a new line precedes it.


      ! Argument declarations:

   CHARACTER(*), INTENT(IN)     :: Str                                          ! The prompt string to print.



   CALL WrScr ( ' ' )
   CALL WrNR  ( TRIM( Str )//' > ' )


   RETURN
   END SUBROUTINE WrPr ! ( Str )
!=======================================================================
   SUBROUTINE WrFileNR ( Unit, Str )


      ! This routine writes out a string to the file connected to Unit without following it with a new line.


      ! Argument declarations.

   INTEGER, INTENT(IN)          :: Unit                                         ! I/O unit for input file.

   CHARACTER(*), INTENT(IN)     :: Str                                          ! String to be written without a newline at the end.



   WRITE (Unit,'(A)',ADVANCE='NO')  Str


   RETURN
   END SUBROUTINE WrFileNR ! ( Unit, Str )
!=======================================================================
   SUBROUTINE WrML ( Str )


      ! This routine writes out a string in the middle of a line.


      ! Argument declarations.

   CHARACTER(*)                 :: Str



   CALL WrNR ( Str )


   RETURN
   END SUBROUTINE WrML ! ( Str )
!=======================================================================
   SUBROUTINE WrScr1 ( Str )


      ! This routine writes out a string to the screen after a blank line.


      ! Argument declarations.

   CHARACTER(*)                 :: Str                                         ! The string to print.



   CALL WrScr ( ' ' )
   CALL WrScr ( TRIM( Str ) )


   RETURN
   END SUBROUTINE WrScr1 ! ( Str )
!=======================================================================

END MODULE NWTC_IO
