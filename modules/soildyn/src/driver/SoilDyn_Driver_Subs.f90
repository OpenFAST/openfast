!**********************************************************************************************************************************
!
!  MODULE: SoilDyn_Driver_Subs  - This module contains subroutines used by the SoilDyn Driver program
!
!**********************************************************************************************************************************
!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2020  National Renewable Energy Laboratory
!
!    This file is part of SoilDyn.
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
MODULE SoilDyn_Driver_Subs

   USE NWTC_Library
   USE SoilDyn_Driver_Types
   IMPLICIT NONE

!  NOTE: This is loosely based on the InflowWind driver code.

CONTAINS
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!> Print out help information
SUBROUTINE DispHelpText()
      !  Statement about usage
   CALL WrScr("")
   CALL WrScr("  Syntax:  SoilDyn_Driver <filename> [options]")
   CALL WrScr("")
   CALL WrScr("       where:     <filename>     -- Name of driver input file to use")
   CALL WrScr("     options:     "//SWChar//"sld           -- treat <filename> as name of SoilDyn input file")
   CALL WrScr("                                    (no driver input file)")
   CALL WrScr("")
   CALL WrScr("              The following options will overwrite values in the driver input file:")
   CALL WrScr("                  "//SwChar//"DT[#]         -- timestep                                        ")
   CALL WrScr("                  "//SwChar//"TStart[#]     -- start time                                      ")
   CALL WrScr("                  "//SwChar//"TSteps[#]     -- number of timesteps                             ")
   CALL WrScr("                  "//SwChar//"v             -- verbose output ")
   CALL WrScr("                  "//SwChar//"vv            -- very verbose output ")
   CALL WrScr("                  "//SwChar//"help          -- print this help menu and exit")
   CALL WrScr("")
   CALL WrScr("   Notes:")
   CALL WrScr("   -- Options are not case sensitive.")
   CALL WrScr("")
END SUBROUTINE DispHelpText


subroutine InitSettingsFlags( ProgInfo, CLSettings, CLFlags )
   implicit none
      ! Storing the arguments
   type( ProgDesc ),                   intent(in   )  :: ProgInfo
   type( SlDDriver_Settings ),         intent(  out)  :: CLSettings        !< Command line arguments passed in
   type( SlDDriver_Flags ),            intent(  out)  :: CLFlags           !< Flags indicating which command line arguments were specified

      ! Set some CLSettings to null/default values
   CLSettings%DvrIptFileName   =  ""             ! No input name name until set
   CLSettings%SlDIptFileName   =  ""             ! No SlD input file name until set
   CLSettings%SummaryFileName  =  ""             ! No summary file name until set
   CLSettings%NumTimeSteps     =  0_IntKi
   CLSettings%DT               =  0.0_DbKi
   CLSettings%TStart           =  0.0_ReKi
   CLSettings%ProgInfo         =  ProgInfo       ! Driver info

      ! Set some CLFlags to null/default values
   CLFlags%DvrIptFile          =  .FALSE.        ! Driver     input filename given as command line argument
   CLFlags%SlDIptFile          =  .FALSE.        ! InflowWind input filename given as command line argument
   CLFlags%TStart              =  .FALSE.        ! specified time to start at
   CLFlags%StiffMatOut         =  .FALSE.        ! stiffness matrix output at start and end
   CLFlags%NumTimeSteps        =  .FALSE.        ! specified a number of timesteps
   CLFlags%NumTimeStepsDefault =  .FALSE.        ! specified 'DEFAULT' for number of timesteps
   CLFlags%DT                  =  .FALSE.        ! specified a resolution in time
   CLFlags%DTDefault           =  .FALSE.        ! specified 'DEFAULT' for resolution in time
   CLFlags%Verbose             =  .FALSE.        ! Turn on verbose error reporting?
   CLFlags%VVerbose            =  .FALSE.        ! Turn on very verbose error reporting?

end subroutine InitSettingsFlags

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!> This subroutine retrieves the command line arguments and passes them to the
!! SoilDyn_driver_subs::parsearg routine for processing.
SUBROUTINE RetrieveArgs( CLSettings, CLFlags, ErrStat, ErrMsg )
      ! Storing the arguments
   type( SlDDriver_Flags ),            intent(  out)  :: CLFlags      !< Flags indicating which command line arguments were specified
   type( SlDDriver_Settings ),         intent(  out)  :: CLSettings           !< Command line arguments passed in
   integer(IntKi),                     intent(  out)  :: ErrStat
   CHARACTER(*),                       intent(  out)  :: ErrMsg

      ! Local variable
   integer(IntKi)                                     :: i                    !< Generic counter
   character(1024)                                    :: Arg                  !< argument given
   character(1024)                                    :: ArgUC                !< Upper case argument to check
   integer(IntKi)                                     :: NumInputArgs         !< Number of argements passed in from command line
   logical                                            :: sldFlag              !< The -sld flag was set
   character(1024)                                    :: FileName             !< Filename from the command line.
   logical                                            :: FileNameGiven        !< Flag indicating if a filename was given.
   integer(IntKi)                                     :: ErrStatTmp           !< Temporary error status (for calls)
   character(1024)                                    :: ErrMsgTmp            !< Temporary error message (for calls)

      ! initialize some things
   CLFlags%DvrIptFile =  .FALSE.
   ErrStat            =  ErrID_None
   ErrStatTmp         =  ErrID_None
   ErrMsg             =  ''
   ErrMsgTmp          =  ''
   sldFlag            =  .FALSE.
   FileNameGiven      =  .FALSE.
   FileName           =  ''

      ! Check how many arguments are passed in
   NumInputArgs = COMMAND_ARGUMENT_COUNT()

      ! exit if we don't have enough
   IF (NumInputArgs == 0) THEN
      CALL SetErrStat(ErrID_Fatal," Insufficient Arguments. Use option "//SwChar//"help for help menu.",  &
         ErrStat,ErrMsg,'RetrieveArgs')
      RETURN
   ENDIF


      ! Loop through all the arguments, and store them
   DO i=1,NumInputArgs
         ! get the ith argument
      CALL get_command_argument(i, Arg)
      ArgUC =  Arg

         ! convert to uppercase
      CALL Conv2UC( ArgUC )

         ! Check to see if it is a control parameter or the filename
      IF ( INDEX( SwChar, ArgUC(1:1) ) > 0 ) THEN

            ! check to see if we asked for help
         IF ( ArgUC(2:5) == "HELP" ) THEN
            CALL DispHelpText()
            CALL ProgExit(0)
         ENDIF


            ! Check the argument and put it where it belongs
            ! chop the SwChar off before passing the argument
         CALL ParseArg( CLSettings, CLFlags, ArgUC(2:), Arg(2:), sldFlag, ErrStatTmp, ErrMsgTmp )
         CALL SetErrStat(ErrStatTmp,ErrMsgTmp,ErrStat,ErrMsg,'RetrieveArgs')
         IF (ErrStat>AbortErrLev) RETURN

      ELSE

            ! since there is no switch character, assume it is the filename, unless we already set one
         IF ( FileNameGiven ) THEN
            CALL SetErrStat(ErrID_Fatal," Multiple driver input filenames given: "//TRIM(FileName)//", "//TRIM(Arg),   &
               ErrStat,ErrMsg,'RetrieveArgs')
            RETURN
         ELSE
            FileName       = TRIM(Arg)
            FileNameGiven  = .TRUE.
         ENDIF

      ENDIF
   END DO


      ! Was a filename given?
   IF ( .NOT. FileNameGiven ) THEN
      CALL SetErrStat( ErrID_Fatal, " No filename given.", ErrStat, ErrMsg, 'RetrieveArgs' )
      RETURN
   ENDIF

      ! Was the -sld flag set?  If so, the filename is the SoilDyn input file.  Otherwise
      ! it is the driver input file.
   IF ( sldFlag ) THEN
      CLSettings%SlDIptFileName  =  TRIM(FileName)
      CLFlags%SlDIptFile =  .TRUE.
   ELSE
      CLSettings%DvrIptFileName  =  TRIM(FileName)
      CLFlags%DvrIptFile =  .TRUE.
   ENDIF



   !-------------------------------------------------------------------------------
   !-------------------------------------------------------------------------------
      CONTAINS


   !-------------------------------------------------------------------------------
   !> Convert a string to a real number
   FUNCTION StringToReal( StringIn, ErrStat )
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat
      CHARACTER(*),                       INTENT(IN   )  :: StringIn

      REAL(ReKi)                                         :: StringToReal
      INTEGER(IntKi)                                     :: ErrStatTmp         ! Temporary variable to hold the error status

         read( StringIn, *, iostat=ErrStatTmp) StringToReal

            ! If that isn't a number, only warn since we can continue by skipping this value
         IF ( ErrStatTmp .ne. 0 ) ErrStat  = ErrID_Warn

   END FUNCTION StringToReal



   !-------------------------------------------------------------------------------
   SUBROUTINE ParseArg( CLSettings, CLFlags, ThisArgUC, ThisArg, sldFlagSet, ErrStat, ErrMsg )
         !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-!
         ! Parse and store the input argument  !
         !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-!

      USE NWTC_Library
      USE SoilDyn_Driver_Types
      USE SoilDyn_Types

      IMPLICIT NONE

         ! Storing the arguments
      TYPE( SlDDriver_Flags ),            INTENT(INOUT)  :: CLFlags      ! Flags indicating which arguments were specified
      TYPE( SlDDriver_Settings ),         INTENT(INOUT)  :: CLSettings           ! Arguments passed in

      CHARACTER(*),                       INTENT(IN   )  :: ThisArgUC            ! The current argument (upper case for testing)
      CHARACTER(*),                       INTENT(IN   )  :: ThisArg              ! The current argument (as passed in for error messages)
      LOGICAL,                            INTENT(INOUT)  :: sldFlagSet           ! Was the -sld flag given?

         ! Error Handling
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg


         ! local variables
      INTEGER(IntKi)                                     :: Delim1               ! where the [ is
      INTEGER(IntKi)                                     :: Delim2               ! where the ] is
      INTEGER(IntKi)                                     :: DelimSep             ! where the : is
      INTEGER(IntKi)                                     :: DelimSep2            ! where the : is
      INTEGER(IntKi)                                     :: DelimSep3            ! where the : is
      REAL(ReKi)                                         :: TempReal             ! temp variable to hold a real

      INTEGER(IntKi)                                     :: ErrStatTmp           ! Temporary error status for calls



         ! Initialize some things
      ErrStat     = ErrID_None
      ErrStatTmp  =  ErrID_None
      ErrMsg      = ''

         ! Get the delimiters -- returns 0 if there isn't one
      Delim1   = INDEX(ThisArgUC,'[')
      Delim2   = INDEX(ThisArgUC,']')
      DelimSep = INDEX(ThisArgUC,':')


         ! check that if there is an opening bracket, then there is a closing one
      IF ( (Delim1 > 0_IntKi ) .and. (Delim2 < Delim1) ) THEN
         CALL SetErrStat(ErrID_Warn," Syntax error in option: '"//SwChar//TRIM(ThisArg)//"'. Ignoring.",   &
            ErrStat,ErrMsg,'ParseArg')
         RETURN
      ENDIF

         ! check that if there is a colon, then there are brackets
      IF ( (DelimSep > 0_IntKi) .and. (Delim1 == 0_IntKi) ) THEN
         CALL SetErrStat(ErrID_Warn," Syntax error in option: '"//SwChar//TRIM(ThisArg)//"'. Ignoring.",   &
            ErrStat,ErrMsg,'ParseArg')
         RETURN
      ENDIF


         ! If no delimeters were given, than this option is simply a flag
      IF ( Delim1 == 0_IntKi ) THEN
            ! check to see if the filename is the name of the SlD input file
         IF       ( ThisArgUC(1:3) == "SLD" )   THEN
            sldFlagSet              = .TRUE.             ! More logic in the routine that calls this one to set things.
            RETURN
         ELSEIF   ( ThisArgUC(1:2) == "VV"  )   THEN
            CLFlags%VVerbose        = .TRUE.
            RETURN
         ELSEIF   ( ThisArgUC(1:1) == "V"   )   THEN
            CLFlags%Verbose         = .TRUE.
            RETURN
         ELSE
            CALL SetErrStat( ErrID_Warn," Unrecognized option '"//SwChar//TRIM(ThisArg)//"'. Ignoring. Use option "//SwChar//"help for list of options.",  &
               ErrStat,ErrMsg,'ParseArg')
         ENDIF

      ENDIF


         ! "DT[#]"
      IF( ThisArgUC(1:Delim1) == "DT["      ) THEN
         TempReal = StringToReal( ThisArgUC(Delim1+1:Delim2-1), ErrStat )
         IF ( ErrStat == ErrID_None ) THEN
            CLFlags%Dt           = .TRUE.
            CLSettings%DT        = abs(TempReal)
         ELSE
            CLFlags%Dt           = .FALSE.
            IF ( ErrStat == ErrID_Warn ) THEN
               CALL SetErrStat(ErrStatTmp," Invalid number in option '"//SwChar//TRIM(ThisArg)//"'. Ignoring.",   &
                  ErrStat,ErrMsg,'ParseArgs')
            ELSE
               CALL SetErrStat( ErrID_Fatal," Something failed in parsing option '"//SwChar//TRIM(ThisArg)//"'.", &
                  ErrStat, ErrMsg, 'ParseArg')
            ENDIF
            RETURN
         ENDIF


         ! "TSTEPS[#]"
      ELSEIF( ThisArgUC(1:Delim1) == "TSTEPS["      ) THEN
         TempReal = StringToReal( ThisArgUC(Delim1+1:Delim2-1), ErrStat )
         IF ( ErrStat == ErrID_None ) THEN
            CLFlags%NumTimeSteps  = .TRUE.
            CLSettings%NumTimeSteps       = nint(abs(TempReal))
         ELSE
            CLFlags%NumTimeSteps  = .FALSE.
            CLSettings%NumTimeSteps = 1_IntKi
            IF ( ErrStat == ErrID_Warn ) THEN
               CALL SetErrStat(ErrStatTmp," Invalid number in option '"//SwChar//TRIM(ThisArg)//"'. Ignoring.",   &
                  ErrStat,ErrMsg,'ParseArgs')
            ELSE
               CALL SetErrStat( ErrID_Fatal," Something failed in parsing option '"//SwChar//TRIM(ThisArg)//"'.", &
                  ErrStat, ErrMsg, 'ParseArg')
            ENDIF
            RETURN
         ENDIF



         ! "TSTART[#]"
      ELSEIF( ThisArgUC(1:Delim1) == "TSTART["      ) THEN
         TempReal = StringToReal( ThisArgUC(Delim1+1:Delim2-1), ErrStat )
         IF ( ErrStat == ErrID_None ) THEN
            CLFlags%TStart          = .TRUE.
            CLSettings%TStart       = abs(TempReal)
         ELSE
            CLFlags%TStart          = .FALSE.
            IF ( ErrStat == ErrID_Warn ) THEN
               CALL SetErrStat(ErrStatTmp," Invalid number in option '"//SwChar//TRIM(ThisArg)//"'. Ignoring.",   &
                  ErrStat,ErrMsg,'ParseArgs')
            ELSE
               CALL SetErrStat( ErrID_Fatal," Something failed in parsing option '"//SwChar//TRIM(ThisArg)//"'.", &
                  ErrStat, ErrMsg, 'ParseArg')
            ENDIF
            RETURN
         ENDIF


      ELSE
         ErrMsg  = " Unrecognized option: '"//SwChar//TRIM(ThisArg)//"'. Ignoring. Use option "//SwChar//"help for list of options."
         ErrStat = ErrID_Warn
      ENDIF

   END SUBROUTINE ParseArg
   !-------------------------------------------------------------------------------

END SUBROUTINE RetrieveArgs


!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!> This subroutine reads the driver input file and sets up the flags and settings
!! for the driver code.  Any settings from the command line options will override
!! this.
SUBROUTINE ReadDvrIptFile( DvrFileName, DvrFlags, DvrSettings, ProgInfo, ErrStat, ErrMsg )

   CHARACTER(1024),                    INTENT(IN   )  :: DvrFileName
   TYPE(SlDDriver_Flags),              INTENT(INOUT)  :: DvrFlags
   TYPE(SlDDriver_Settings),           INTENT(INOUT)  :: DvrSettings
   TYPE(ProgDesc),                     INTENT(IN   )  :: ProgInfo
   INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat              ! returns a non-zero value when an error occurs
   CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg               ! Error message if ErrStat /= ErrID_None

      ! Local variables
   INTEGER(IntKi)                                     :: UnIn                 ! Unit number for the driver input file
   CHARACTER(1024)                                    :: FileName             ! Name of SoilDyn driver input file

      ! Input file echoing
   LOGICAL                                            :: EchoFileContents     ! Do we echo the driver file out or not?
   INTEGER(IntKi)                                     :: UnEchoLocal          ! The local unit number for this module's echo file
   CHARACTER(1024)                                    :: EchoFileName         ! Name of SoilDyn driver echo file

      ! Time steps
   CHARACTER(1024)                                    :: NumTimeStepsChr      ! Character string for number of timesteps (to handle DEFAULT value)
   CHARACTER(1024)                                    :: DTChr                ! Character string for timesteps size (to handle DEFAULT value)

      ! Gridded data
   INTEGER(IntKi)                                     :: TmpIntAr3(3)         ! Temporary array for reading in a pair of integer values from the input file
   REAL(ReKi)                                         :: TmpRealAr3(3)        ! Temporary array for reading in a pair of real values from the input file
   REAL(ReKi)                                         :: GridCtrCoord(3)      ! Center coordinate of the grid read in

      ! Local error handling
   INTEGER(IntKi)                                     :: ios                  !< I/O status
   INTEGER(IntKi)                                     :: ErrStatTmp           !< Temporary error status for calls
   CHARACTER(1024)                                    :: ErrMsgTmp            !< Temporary error messages for calls


      ! Initialize the echo file unit to -1 which is the default to prevent echoing, we will alter this based on user input
   UnEchoLocal = -1

   FileName = TRIM(DvrFileName)

   CALL GetNewUnit( UnIn )
   CALL OpenFInpFile( UnIn, FileName, ErrStatTmp, ErrMsgTmp )
   IF ( ErrStatTmp /= ErrID_None ) THEN
      CALL SetErrStat(ErrID_Fatal,' Failed to open SoilDyn Driver input file: '//FileName,   &
         ErrStat,ErrMsg,'ReadDvrIptFile')
      CLOSE( UnIn )
      RETURN
   ENDIF


   CALL WrScr( 'Opening SoilDyn Driver input file:  '//trim(FileName) )


   !-------------------------------------------------------------------------------------------------
   ! File header
   !-------------------------------------------------------------------------------------------------

   CALL ReadCom( UnIn, FileName,' SoilDyn Driver input file header line 1', ErrStatTmp, ErrMsgTmp )
   if (Failed()) return

   CALL ReadCom( UnIn, FileName, 'SoilDyn Driver input file header line 2', ErrStatTmp, ErrMsgTmp )
   if (Failed()) return

   CALL ReadCom( UnIn, FileName, 'SoilDyn Driver input file seperator line', ErrStatTmp, ErrMsgTmp )
   if (Failed()) return

     ! Echo Input Files.
   CALL ReadVar ( UnIn, FileName, EchoFileContents, 'Echo', 'Echo Input', ErrStatTmp, ErrMsgTmp )
   if (Failed()) return


      ! If we are Echoing the input then we should re-read the first three lines so that we can echo them
      ! using the NWTC_Library routines.  The echoing is done inside those routines via a global variable
      ! which we must store, set, and then replace on error or completion.

   IF ( EchoFileContents ) THEN

      EchoFileName = TRIM(FileName)//'.ech'
      CALL GetNewUnit( UnEchoLocal )
      CALL OpenEcho ( UnEchoLocal, EchoFileName, ErrStatTmp, ErrMsgTmp, ProgInfo )
      if (Failed()) return

      REWIND(UnIn)

         ! Reread and echo
      CALL ReadCom( UnIn, FileName,' SoilDyn Driver input file header line 1', ErrStatTmp, ErrMsgTmp, UnEchoLocal )
      if (Failed()) return

      CALL ReadCom( UnIn, FileName, 'SoilDyn Driver input file header line 2', ErrStatTmp, ErrMsgTmp, UnEchoLocal )
      if (Failed()) return

      CALL ReadCom( UnIn, FileName, 'SoilDyn Driver input file seperator line', ErrStatTmp, ErrMsgTmp, UnEchoLocal )
      if (Failed()) return

        ! Echo Input Files.
      CALL ReadVar ( UnIn, FileName, EchoFileContents, 'Echo', 'Echo Input', ErrStatTmp, ErrMsgTmp, UnEchoLocal )
      if (Failed()) return

   ENDIF


   !-------------------------------------------------------------------------------------------------
   !  Driver setup section
   !-------------------------------------------------------------------------------------------------

      ! Header
   CALL ReadCom( UnIn, FileName,' Driver setup section, comment line', ErrStatTmp, ErrMsgTmp, UnEchoLocal )
   if (Failed()) return

      ! SoilDyn input file
   CALL ReadVar( UnIn, FileName,DvrSettings%SlDIptFileName,'SlDIptFileName',' SoilDyn input filename',   &
      ErrStatTmp,ErrMsgTmp, UnEchoLocal )
   if (Failed()) then
      return
   else
      DvrFlags%SlDIptFile  =  .TRUE.
   endif


      ! TStart    -- start time
   CALL ReadVar( UnIn, FileName,DvrSettings%TStart,'TStart',' Time in wind file to start parsing.',   &
      ErrStatTmp,ErrMsgTmp, UnEchoLocal )
   if (Failed()) then
      return
   else
      DvrFlags%TStart   =  .TRUE.
   endif


      ! DT    -- Timestep size for the driver to take (or DEFAULT for what the file contains)
   CALL ReadVar( UnIn, FileName,DTChr,'DTChr',' Character string for Timestep size for the driver to take (or DEFAULT for what the file contains).',  &
      ErrStatTmp,ErrMsgTmp, UnEchoLocal )
   if (Failed()) return

      ! Check if we asked for the DEFAULT (use what is in the file)
   CALL Conv2UC( DTChr )
   IF ( TRIM(DTChr) == 'DEFAULT' ) THEN     ! we asked for the default value
      DvrFlags%DT         =  .TRUE.
      DvrFlags%DTDefault  =  .TRUE.         ! This flag tells us to use the inflow wind file values
   ELSE
         !  We probably have a number if it isn't 'DEFAULT', so do an internal read and check to
         !  make sure that it was appropriately interpretted.
      READ (DTChr,*,IOSTAT=IOS)   DvrSettings%DT
      IF ( IOS /= 0 )  THEN  ! problem in the read, so parse the error.
         CALL CheckIOS ( IOS, '', 'DT',NumType, ErrStatTmp, ErrMsgTmp )
         if (Failed()) return
      ELSE     ! Was ok, so set the flags
         DvrFlags%DT         =  .TRUE.
         DvrFlags%DTDefault  =  .FALSE.
      ENDIF
   ENDIF


      ! Number of timesteps
   CALL ReadVar( UnIn, FileName,NumTimeStepsChr,'NumTimeStepsChr',' Character string for number of timesteps to read.',   &
      ErrStatTmp,ErrMsgTmp, UnEchoLocal )
   if (Failed()) return

      ! Check if we asked for the DEFAULT (use what is in the file)
   CALL Conv2UC( NumTimeStepsChr )
   IF ( TRIM(NumTimeStepsChr) == 'DEFAULT' ) THEN     ! we asked for the default value
      DvrFlags%NumTimeSteps         =  .TRUE.
      DvrFlags%NumTimeStepsDefault  =  .TRUE.         ! This flag tells us to use the inflow wind file values
   ELSE
         !  We probably have a number if it isn't 'DEFAULT', so do an internal read and check to
         !  make sure that it was appropriately interpretted.
      READ (NumTimeStepsChr,*,IOSTAT=IOS)   DvrSettings%NumTimeSteps
      IF ( IOS /= 0 )  THEN  ! problem in the read, so parse the error.
         CALL CheckIOS ( IOS, '', 'NumTimeSteps',NumType, ErrStatTmp, ErrMsgTmp )
         if (Failed()) return
      ELSE     ! Was ok, so set the flags
         DvrFlags%NumTimeSteps         =  .TRUE.
         DvrFlags%NumTimeStepsDefault  =  .FALSE.
      ENDIF
   ENDIF


      ! Stiffness matrix
   CALL ReadVar( UnIn, FileName,DvrFlags%StiffMatOut,'StiffMatOut',' Output stiffness matrices at start and end',   &
      ErrStatTmp,ErrMsgTmp, UnEchoLocal )
   if (Failed()) return


   !-------------------------------------------------------------------------------------------------
   !  SoilDyn time series input
   !-------------------------------------------------------------------------------------------------
!FIXME: add these
      ! Close the echo and input file
   CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
   CLOSE( UnIn )


CONTAINS

   !> Set error status, close stuff, and return
   logical function Failed()
      CALL SetErrStat(ErrStatTmp,ErrMsgTmp,ErrStat,ErrMsg,'ReadDvrIptFile')
      if (ErrStat >= AbortErrLev) then
         CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
         CLOSE( UnIn )
      endif
      Failed =    ErrStat >= AbortErrLev
   end function Failed

   !> Clean up the module echo file
   subroutine CleanupEchoFile( EchoFlag, UnEcho)
      logical,                      intent(in   )  :: EchoFlag          ! local version of echo flag
      integer(IntKi),               intent(in   )  :: UnEcho            !  echo unit number

         ! Close this module's echo file
      if ( EchoFlag ) then
         close(UnEcho)
      endif
   END SUBROUTINE CleanupEchoFile

END SUBROUTINE ReadDvrIptFile


!> This subroutine copies an command line (CL) settings over to the program settings.  Warnings are
!! issued if anything is changed from what the driver input file requested.
SUBROUTINE UpdateSettingsWithCL( DvrFlags, DvrSettings, CLFlags, CLSettings, DVRIPT, ErrStat, ErrMsg )

   TYPE(SlDDriver_Flags),              INTENT(INOUT)  :: DvrFlags
   TYPE(SlDDriver_Settings),           INTENT(INOUT)  :: DvrSettings
   TYPE(SlDDriver_Flags),              INTENT(IN   )  :: CLFlags
   TYPE(SlDDriver_Settings),           INTENT(IN   )  :: CLSettings
   LOGICAL,                            INTENT(IN   )  :: DVRIPT
   INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat
   CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg


      ! Local variables
   INTEGER(IntKi)                                     :: ErrStatTmp        !< Temporary error status for calls
   CHARACTER(1024)                                    :: ErrMsgTmp         !< Temporary error status for calls
   LOGICAL                                            :: WindGridModify    !< Did we modify any of the WindGrid related settings?

      ! Initialization
   WindGridModify =  .FALSE.

      ! Initialize the error handling
   ErrStat     =  ErrID_None
   ErrMsg      =  ''
   ErrStatTmp  =  ErrID_None
   ErrMsgTmp   =  ''


      !--------------------------------------------
      !  Did we change any time information?
      !--------------------------------------------

      !  Check TStart
   IF ( CLFlags%TStart ) THEN
      IF ( DvrFlags%TStart .AND. ( .NOT. EqualRealNos(DvrSettings%TStart, CLSettings%TStart) ) )   THEN
         CALL SetErrStat( ErrID_Warn, ' Overriding driver input value for TStart with '//TRIM(Num2LStr(CLSettings%TStart))//'.', &
            ErrStat,ErrMsg,'UpdateSettingsWithCL')
      ELSE
         DvrFlags%TStart   =  .TRUE.
      ENDIF
      DvrSettings%TStart   =  CLSettings%TStart
   ENDIF

      ! Check DT
   IF ( CLFlags%DT ) THEN
      IF ( DvrFlags%DT .AND. ( .NOT. EqualRealNos(DvrSettings%DT, CLSettings%DT) ) )   THEN
         CALL SetErrStat( ErrID_Warn, ' Overriding driver input value for DT with '//TRIM(Num2LStr(CLSettings%DT))//'.',  &
            ErrStat,ErrMsg,'UpdateSettingsWithCL')
      ELSE
         DvrFlags%DT       =  .TRUE.
      ENDIF
      DvrSettings%DT       =  CLSettings%DT
      DvrFlags%DTDefault   =  .FALSE.
   ENDIF

      ! Check NumTimeSteps
   IF ( CLFlags%NumTimeSteps ) THEN
      IF ( DvrFlags%NumTimeSteps .AND. ( DvrSettings%NumTimeSteps /= CLSettings%NumTimeSteps ) )   THEN
         CALL SetErrStat( ErrID_Warn, ' Overriding driver input value for NumTimeSteps with '//  &
            TRIM(Num2LStr(CLSettings%NumTimeSteps))//'.',&
            ErrStat,ErrMsg,'UpdateSettingsWithCL')
      ELSE
         DvrFlags%NumTimeSteps      =  .TRUE.
      ENDIF
      DvrSettings%NumTimeSteps      =  CLSettings%NumTimeSteps
      DvrFlags%NumTimeStepsDefault  =  .FALSE.
   ENDIF

      ! Make sure there is at least one timestep
   DvrSettings%NumTimeSteps   =  MAX(DvrSettings%NumTimeSteps,1_IntKi)



      !--------------------------------------------
      ! If there was no driver input file, we need to set a few things.
      !--------------------------------------------

   IF ( .NOT. DVRIPT ) THEN

         ! Do we need to set the NumTimeStepsDefault flag?
      IF ( .NOT. DvrFlags%NumTimeSteps )  THEN
         DvrFlags%NumTimeStepsDefault  =  .TRUE.
         CALL SetErrStat( ErrID_Info,' The number of timesteps is not specified.  Defaulting to what is in the wind file.', &
            ErrStat,ErrMsg,'UpdateSettingsWithCL')
      ENDIF
   ENDIF


!FIXME: remove this after parsing rest of input file.
      ! If no DT value has been set (DEFAULT requested), we need to set a default to pass into SlD
   IF ( .NOT. DvrFlags%DT ) THEN
      DvrSettings%DT =  0.025_DbKi     ! This value gets passed into the SlD_Init routine, so something must be set.
   ENDIF


END SUBROUTINE UpdateSettingsWithCL



!> This routine exists only to support the development of the module.  It will not be needed after the module is complete.
SUBROUTINE  printSettings( DvrFlags, DvrSettings )
      ! The arguments
   TYPE( SlDDriver_Flags ),            INTENT(IN   )  :: DvrFlags           !< Flags indicating which settings were set
   TYPE( SlDDriver_Settings ),         INTENT(IN   )  :: DvrSettings        !< Stored settings

   CALL WrsCr(TRIM(GetNVD(DvrSettings%ProgInfo)))
   CALL WrScr(' DvrIptFile:          '//FLAG(DvrFlags%DvrIptFile)//        '      '//TRIM(DvrSettings%DvrIptFileName))
   CALL WrScr(' SlDIptFile:          '//FLAG(DvrFlags%SlDIptFile)//        '      '//TRIM(DvrSettings%SlDIptFileName))
   CALL WrScr(' TStart:              '//FLAG(DvrFlags%TStart)//            '      '//TRIM(Num2LStr(DvrSettings%TStart)))
   IF ( DvrFlags%DTDefault) THEN
      CALL WrScr(' DT:                  '//FLAG(DvrFlags%DT)//                '      DEFAULT')
   ELSE
      CALL WrScr(' DT:                  '//FLAG(DvrFlags%DT)//                '      '//TRIM(Num2LStr(DvrSettings%DT)))
   ENDIF
   IF ( DvrFlags%NumTimeStepsDefault) THEN
      CALL WrScr(' NumTimeSteps:        '//FLAG(DvrFlags%NumTimeSteps)//      '      DEFAULT')
   ELSE
      CALL WrScr(' NumTimeSteps:        '//FLAG(DvrFlags%NumTimeSteps)//      '      '//TRIM(Num2LStr(DvrSettings%NumTimeSteps)))
   ENDIF
   CALL WrScr(' StiffMatOut:         '//FLAG(DvrFlags%StiffMatOut))
   CALL WrScr(' Verbose:             '//FLAG(DvrFlags%Verbose))
   CALL WrScr(' VVerbose:            '//FLAG(DvrFlags%VVerbose))
   RETURN
END SUBROUTINE printSettings


!> This routine exists only to support the development of the module.  It will not be kept after the module is complete.
!! This routine takes a flag setting (LOGICAL) and exports either 'T' or '-' for T/F (respectively)
FUNCTION FLAG(flagval)
   LOGICAL,       INTENT(IN   )  :: flagval        !< Value of the flag
   CHARACTER(1)                  :: FLAG           !< character interpretation (for prettiness when printing)
   IF ( flagval ) THEN
      FLAG  =  'T'
   ELSE
      FLAG  =  '-'
   ENDIF
   RETURN
END FUNCTION FLAG

END MODULE SoilDyn_Driver_Subs
