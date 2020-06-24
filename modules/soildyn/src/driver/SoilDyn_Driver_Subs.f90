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
   CALL WrScr("                  "//SwChar//"NonLinear     -- only return non-linear portion of reaction force")
   CALL WrScr("                  "//SwChar//"help          -- print this help menu and exit")
   CALL WrScr("")
   CALL WrScr("   Notes:")
   CALL WrScr("   -- Options are not case sensitive.")
   CALL WrScr("")
!FIXME: update this
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
   CLSettings%InputDispFile    =  ""             ! No SlD input displacement timeseries file name until set
   CLSettings%NumTimeSteps     =  0_IntKi
   CLSettings%DT               =  0.0_DbKi
   CLSettings%TStart           =  0.0_ReKi
   CLSettings%ProgInfo         =  ProgInfo       ! Driver info

      ! Set some CLFlags to null/default values
   CLFlags%DvrIptFile          =  .FALSE.        ! Driver     input filename given as command line argument
   CLFlags%SlDIptFile          =  .FALSE.        ! SoilDyn input filename given as command line argument
   CLFlags%InputDispFile       =  .FALSE.        ! No SlD input displacement timeseries file name until set
   CLFlags%TStart              =  .FALSE.        ! specified time to start at
   CLFlags%StiffMatOut         =  .FALSE.        ! stiffness matrix output at start and end
   CLFlags%NumTimeSteps        =  .FALSE.        ! specified a number of timesteps
   CLFlags%NumTimeStepsDefault =  .FALSE.        ! specified 'DEFAULT' for number of timesteps
   CLFlags%DT                  =  .FALSE.        ! specified a resolution in time
   CLFlags%DTDefault           =  .FALSE.        ! specified 'DEFAULT' for resolution in time
   CLFlags%Verbose             =  .FALSE.        ! Turn on verbose error reporting?
   CLFlags%VVerbose            =  .FALSE.        ! Turn on very verbose error reporting?
   CLFlags%SlDNonLinearForcePortionOnly =  .FALSE. ! Report only non-linear portion of forces

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
         IF   ( ThisArgUC(1:9) == "NONLINEAR"   )   THEN
            CLFlags%SlDNonLinearForcePortionOnly = .TRUE.
            RETURN
         ELSEIF   ( ThisArgUC(1:3) == "SLD" )   THEN
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
!FIXME: add in the other inputs here.

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
   CHARACTER(1024)                                    :: InputChr             ! Character string for timesteps and input file names (to handle DEFAULT or NONE value)

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
   CALL ReadVar( UnIn, FileName,InputChr,'InputChr',' Character string for Timestep size for the driver to take (or DEFAULT for what the file contains).',  &
      ErrStatTmp,ErrMsgTmp, UnEchoLocal )
   if (Failed()) return

      ! Check if we asked for the DEFAULT (use what is in the file)
   CALL Conv2UC( InputChr )
   IF ( TRIM(InputChr) == 'DEFAULT' ) THEN     ! we asked for the default value
      DvrFlags%DT         =  .TRUE.
      DvrFlags%DTDefault  =  .TRUE.         ! This flag tells us to use the inflow wind file values
   ELSE
         !  We probably have a number if it isn't 'DEFAULT', so do an internal read and check to
         !  make sure that it was appropriately interpretted.
      READ (InputChr,*,IOSTAT=IOS)   DvrSettings%DT
      IF ( IOS /= 0 )  THEN  ! problem in the read, so parse the error.
         CALL CheckIOS ( IOS, '', 'DT',NumType, ErrStatTmp, ErrMsgTmp )
         if (Failed()) return
      ELSE     ! Was ok, so set the flags
         DvrFlags%DT         =  .TRUE.
         DvrFlags%DTDefault  =  .FALSE.
      ENDIF
   ENDIF


      ! Number of timesteps
   CALL ReadVar( UnIn, FileName,InputChr,'InputChr',' Character string for number of timesteps to read.',   &
      ErrStatTmp,ErrMsgTmp, UnEchoLocal )
   if (Failed()) return

      ! Check if we asked for the DEFAULT (use what is in the file)
   CALL Conv2UC( InputChr )
   IF ( TRIM(InputChr) == 'DEFAULT' ) THEN     ! we asked for the default value
      DvrFlags%NumTimeSteps         =  .FALSE.
      DvrFlags%NumTimeStepsDefault  =  .TRUE.         ! This flag tells us to use the inflow wind file values
   ELSE
         !  We probably have a number if it isn't 'DEFAULT', so do an internal read and check to
         !  make sure that it was appropriately interpretted.
      READ (InputChr,*,IOSTAT=IOS)   DvrSettings%NumTimeSteps
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

      ! Non-linear reaction portion only
   CALL ReadVar( UnIn, FileName,DvrFlags%SlDNonLinearForcePortionOnly,'SlDNonLinearForcePortionOnly',' Only report the non-linear portion of the reaction force.',   &
      ErrStatTmp,ErrMsgTmp, UnEchoLocal )
   if (Failed()) return



   !-------------------------------------------------------------------------------------------------
   !  SoilDyn time series input -- this is read from a file of 7 columns (time and 6 dof)
   !-------------------------------------------------------------------------------------------------

      ! InputDispFile input file
   CALL ReadVar( UnIn, FileName,InputChr,'InputDispFile',' SoilDyn input displacements filename',   &
      ErrStatTmp,ErrMsgTmp, UnEchoLocal )
   if (Failed()) return

   DvrSettings%InputDispFile  = InputChr
   call Conv2UC( InputChr )
   if (trim(InputChr) == 'NONE') then
      DvrSettings%InputDispFile  = ''
      DvrFlags%InputDispFile  =  .FALSE.
   else
      DvrFlags%InputDispFile  =  .TRUE.
   endif


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
   character(*), parameter                            :: RoutineName = 'UpdateSettingsWithCL'

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
            ErrStat,ErrMsg,RoutineName)
      ELSE
         DvrFlags%TStart   =  .TRUE.
      ENDIF
      DvrSettings%TStart   =  CLSettings%TStart
   ENDIF

      ! Check DT
   IF ( CLFlags%DT ) THEN
      IF ( DvrFlags%DT .AND. ( .NOT. EqualRealNos(DvrSettings%DT, CLSettings%DT) ) )   THEN
         CALL SetErrStat( ErrID_Warn, ' Overriding driver input value for DT with '//TRIM(Num2LStr(CLSettings%DT))//'.',  &
            ErrStat,ErrMsg,RoutineName)
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
            ErrStat,ErrMsg,RoutineName)
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
         CALL SetErrStat( ErrID_Info,' The number of timesteps is not specified.  Defaulting to what is in the input series file.', &
            ErrStat,ErrMsg,RoutineName)
      ENDIF
   ENDIF


!FIXME: remove this after parsing rest of input file.
      ! If no DT value has been set (DEFAULT requested), we need to set a default to pass into SlD
   IF ( .NOT. DvrFlags%DT ) THEN
      DvrSettings%DT =  0.025_DbKi     ! This value gets passed into the SlD_Init routine, so something must be set.
   ENDIF


END SUBROUTINE UpdateSettingsWithCL



SUBROUTINE ReadInputDispFile( InputDispFile, DisplacementList, ErrStat, ErrMsg )
   CHARACTER(1024),                    INTENT(IN   )  :: InputDispFile       !< Name of the points file to read
   REAL(R8Ki), ALLOCATABLE,            INTENT(  OUT)  :: DisplacementList(:,:)       !< The coordinates we read in: idx 1 = timestep, idx 2 = values
   INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat              !< The error status
   CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg               !< The message for the status

      ! Local variables
   CHARACTER(1024)                                    :: ErrMsgTmp            !< Temporary error message for calls
   INTEGER(IntKi)                                     :: ErrStatTmp           !< Temporary error status for calls
   INTEGER(IntKi)                                     :: FiUnitPoints         !< Unit number for points file to open

   INTEGER(IntKi)                                     :: NumDataColumns       !< Number of data columns
   INTEGER(IntKi)                                     :: NumDataPoints        !< Number of lines of data (one point per line)
   INTEGER(IntKi)                                     :: NumHeaderLines       !< Number of header lines to ignore

   INTEGER(IntKi)                                     :: I                    !< Generic counter
   character(*), parameter                            :: RoutineName = 'ReadInputDispFile'

      ! Initialization of subroutine
   ErrMsg      =  ''
   ErrMsgTmp   =  ''
   ErrStat     =  ErrID_None
   ErrStatTmp  =  ErrID_None


      ! Now open file
   CALL GetNewUnit( FiUnitPoints, ErrStatTmp, ErrMsgTmp ); if (Failed()) return
   CALL OpenFInpFile(   FiUnitPoints,  TRIM(InputDispFile), ErrStatTmp, ErrMsgTmp )   ! Unformatted input file
      if (Failed()) return

      ! Find out how long the file is
   CALL GetFileLength( FiUnitPoints, InputDispFile, NumDataColumns, NumDataPoints, NumHeaderLines, ErrMsgTmp, ErrStatTmp )
      if (Failed()) return
   IF ( NumDataColumns /= 7 ) THEN
      ErrStatTmp  = ErrID_Fatal
      ErrMsgTmp   = ' Expecting seven columns in '//TRIM(InputDispFile)//' corresponding to '//   &
         'time, dX, dY, dZ, dTheta_X, dTheta_Y, dTheta_Z  coordinates.  Instead found '//TRIM(Num2LStr(NumDataColumns))//' columns.'
      if (Failed()) return
   ENDIF


      ! Allocate the storage for the data
   CALL AllocAry( DisplacementList, NumDataPoints, 7, "Array of Points data", ErrStatTmp, ErrMsgTmp )
      if (Failed()) return


      ! Read in the headers and throw them away
   DO I=1,NumHeaderLines
      CALL ReadCom( FiUnitPoints, InputDispFile,' Points file header line', ErrStatTmp, ErrMsgTmp )
      if (Failed()) return
   ENDDO

      ! Read in the datapoints   -- This is arranged with time in first index for speed in later interpolation operations
   DO I=1,NumDataPoints
      CALL ReadAry ( FiUnitPoints, InputDispFile, DisplacementList(I,:), 7, 'DisplacementList', &
         'Coordinate point from Points file', ErrStatTmp, ErrMsgTmp)
      if (Failed()) return
   ENDDO

   CLOSE( FiUnitPoints )

CONTAINS
   !> Set error status, close stuff, and return
   logical function Failed()
      CALL SetErrStat(ErrStatTmp,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
      if (ErrStat >= AbortErrLev .and. FiUnitPoints >0) close( FiUnitPoints )
      Failed =    ErrStat >= AbortErrLev
   end function Failed


  !-------------------------------------------------------------------------------------------------------------------------------
   !>    This subroutine looks at a file that has been opened and finds out how many header lines there are, how many columns there
   !!    are, and    how many lines of data there are in the file.
   !!
   !!    A few things are assumed about the file:
   !!       1. Any header lines are the first thing in the file.
   !!       2. No text appears anyplace other than in first part of the file
   !!       3. The datalines only contain numbers that can be read in as reals.
   !!
   !!    Limitations:
   !!       1. only handles up to 20 words (columns) on a line
   !!       2. empty lines are considered text lines
   !!       3. All data rows must contain the same number of columns
   !!
   !!
   SUBROUTINE GetFileLength(UnitDataFile, DataFileName, NumDataColumns, NumDataLines, NumHeaderLines, ErrMsg, ErrStat)

      INTEGER(IntKi),                     INTENT(IN   )  :: UnitDataFile      !< Unit number of the file we are looking at.
      CHARACTER(*),                       INTENT(IN   )  :: DataFileName      !< The name of the file we are looking at.
      INTEGER(IntKi),                     INTENT(  OUT)  :: NumDataColumns    !< The number of columns in the data file.
      INTEGER(IntKi),                     INTENT(  OUT)  :: NumDataLines      !< Number of lines containing data
      INTEGER(IntKi),                     INTENT(  OUT)  :: NumHeaderLines    !< Number of header lines at the start of the file
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg            !< Error Message to return (empty if all good)
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat           !< Status flag if there were any problems (ErrID_None if all good)

         ! Local Variables
      CHARACTER(2048)                                    :: ErrMsgTmp         !< Temporary message variable.  Used in calls.
      INTEGER(IntKi)                                     :: ErrStatTmp        !< Temporary error status.  Used in calls.
      INTEGER(IntKi)                                     :: LclErrStat        !< Temporary error status.  Used locally to indicate when we have reached the end of the file.
      INTEGER(IntKi)                                     :: TmpIOErrStat      !< Temporary error status for the internal read of the first word to a real number
      LOGICAL                                            :: IsRealNum         !< Flag indicating if the first word on the line was a real number

      CHARACTER(1024)                                    :: TextLine          !< One line of text read from the file
      INTEGER(IntKi)                                     :: LineLen           !< The length of the line read in
      CHARACTER(1024)                                    :: StrRead           !< String containing the first word read in
      REAL(R8Ki)                                         :: RealRead          !< Returns value of the number (if there was one), or NaN (as set by NWTC_Num) if there wasn't
      CHARACTER(24)                                      :: Words(20)         !< Array of words we extract from a line.  We shouldn't have more than 20.
      INTEGER(IntKi)                                     :: i                 !< simple integer counters
      INTEGER(IntKi)                                     :: LineNumber        !< the line I am on
      LOGICAL                                            :: LineHasText       !< Flag indicating if the line I just read has text.  If so, it is a header line.
      LOGICAL                                            :: HaveReadData      !< Flag indicating if I have started reading data.
      INTEGER(IntKi)                                     :: NumWords          !< Number of words on a line
      INTEGER(IntKi)                                     :: FirstDataLineNum  !< Line number of the first row of data in the file

         ! Initialize the error handling
      ErrStat     = ErrID_None
      ErrStatTmp  = ErrID_None
      LclErrStat  = ErrID_None
      ErrMsg      = ''
      ErrMsgTmp   = ''

         ! Set some of the flags and counters
      HaveReadData   = .FALSE.
      NumDataColumns = 0
      NumHeaderLines = 0
      NumDataLines   = 0
      LineNumber     = 0

         ! Just in case we were handed a file that we are part way through reading (should never be true), rewind to the start

      REWIND( UnitDataFile )

         !------------------------------------
         !> The variable LclErrStat is used to indicate when we have reached the end of the file or had an error from
         !! ReadLine.  Until that occurs, we read each line, and decide if it contained any non-numeric data.  The
         !! first group of lines containing non-numeric data is considered the header.  The first line of all numeric
         !! data is considered the start of the data section.  Any non-numeric containing found within the data section
         !! will be considered as an invalid file format at which point we will return a fatal error from this routine.

      DO WHILE ( LclErrStat == ErrID_None )

            !> Reset the indicator flag for the non-numeric content
         LineHasText = .FALSE.

            !> Read in a single line from the file
         CALL ReadLine( UnitDataFile, '', TextLine, LineLen, LclErrStat )

            !> If there was an error in reading the file, then exit.
            !!    Possible causes: reading beyond end of file in which case we are done so don't process it.
         IF ( LclErrStat /= ErrID_None ) EXIT

            !> Increment the line counter.
         LineNumber  = LineNumber + 1

            !> Read all the words on the line into the array called 'Words'.  Only the first words will be encountered
            !! will be stored.  The others are empty (i.e. only three words on the line, so the remaining 17 are empty).
         CALL GetWords( TextLine, Words, 20 )

            !> Cycle through and count how many are not empty.  Once an empty value is encountered, all the rest should
            !! be empty if GetWords worked correctly.  The index of the last non-empty value is stored.
         DO i=1,20
            IF (TRIM(Words(i)) .ne. '') NumWords=i
         ENDDO


            !> Now cycle through the first 'NumWords' of non-empty values stored in 'Words'.  Words should contain
            !! everything that is one the line.  The subroutine ReadRealNumberFromString will set a flag 'IsRealNum'
            !! when the value in Words(i) can be read as a real(R8Ki).  'StrRead' will contain the string equivalent.
         DO i=1,NumWords
            CALL ReadRealNumberFromString( Words(i), RealRead, StrRead, IsRealNum, ErrStatTmp, ErrMsgTmp, TmpIOErrStat )
            IF ( .NOT. IsRealNum)   LineHasText = .TRUE.
         ENDDO

            !> If all the words on that line had no text in them, then it must have been a line of data.
            !! If not, then we have either a header line, which is ok, or a line containing text in the middle of the
            !! the data section, which is not good (the flag HaveReadData tells us which case this is).
         IF ( LineHasText ) THEN
            IF ( HaveReadData ) THEN      ! Uh oh, we have already read a line of data before now, so there is a problem
               CALL SetErrStat( ErrID_Fatal, ' Found text on line '//TRIM(Num2LStr(LineNumber))//' of '//TRIM(DataFileName)// &
                           ' when real numbers were expected.  There may be a problem with format of the file: '// &
                           TRIM(DataFileName)//'.', ErrStat, ErrMsg, RoutineName)
               IF ( ErrStat >= AbortErrLev )    RETURN
            ELSE
               NumHeaderLines = NumHeaderLines + 1
            ENDIF
         ELSE     ! No text, must be data line
            NumDataLines = NumDataLines + 1
               ! If this is the first row of data, then store the number of words that were on the line
            IF ( .NOT. HaveReadData )  THEN
                  ! If this is the first line of data, keep some relevant info about it and the number of columns in it
               HaveReadData      = .TRUE.
               FirstDataLineNum  = LineNumber         ! Keep the line number of the first row of data (for error reporting)
               NumDataColumns    = NumWords
            ELSE
                  ! Make sure that the number columns on the row matches the number of columnns on the first row of data.
               IF ( NumWords /= NumDataColumns ) THEN
                  CALL SetErrStat( ErrID_Fatal, ' Error in file: '//TRIM(DataFileName)//'.'// &
                           ' The number of data columns on line '//TRIM(Num2LStr(LineNumber))// &
                           '('//TRIM(Num2LStr(NumWords))//' columns) is different than the number of columns on first row of data '// &
                           ' (line: '//TRIM(Num2LStr(FirstDataLineNum))//', '//TRIM(Num2LStr(NumDataColumns))//' columns).', &
                           ErrStat, ErrMsg, RoutineName)
                  IF ( ErrStat >= AbortErrLev )    RETURN
               ENDIF
            ENDIF
         ENDIF

      ENDDO

      REWIND( UnitDataFile )

   END SUBROUTINE GetFileLength

   !-------------------------------------------------------------------------------
   !> This subroutine takes a line of text that is passed in and reads the first
   !! word to see if it is a number.  An internal read is used to do this.  If
   !! it is a number, it is started in ValueRead and returned. The flag IsRealNum
   !! is set to true.  Otherwise, ValueRead is set to NaN (value from the NWTC_Num)
   !! and the flag is set to false.
   !!
   !! The IsRealNum flag is set to indicate if we actually have a real number or
   !! not.  After calling this routine, a simple if statement can be used:
   !!
   !!       @code
   !!    IF (IsRealNum) THEN
   !!       ! do something
   !!    ELSE
   !!       ! do something else
   !!    ENDIF
   !!       @endcode
   !!
   !-------------------------------------------------------------------------------
   SUBROUTINE ReadRealNumberFromString(StringToParse, ValueRead, StrRead, IsRealNum, ErrStat, ErrMsg, IOErrStat)
      CHARACTER(*),        INTENT(IN   )           :: StringToParse  !< The string we were handed.
      REAL(R8Ki),          INTENT(  OUT)           :: ValueRead      !< The variable being read.  Returns as NaN (library defined) if not a Real.
      CHARACTER(*),        INTENT(  OUT)           :: StrRead        !< A string containing what was read from the ReadNum routine.
      LOGICAL,             INTENT(  OUT)           :: IsRealNum      !< Flag indicating if we successfully read a Real
      INTEGER(IntKi),      INTENT(  OUT)           :: ErrStat        !< ErrID level returned from ReadNum
      CHARACTER(*),        INTENT(  OUT)           :: ErrMsg         !< Error message including message from ReadNum
      INTEGER(IntKi),      INTENT(  OUT)           :: IOErrStat      !< Error status from the internal read. Useful for diagnostics.

      ErrStat     = ErrID_None
      ErrMsg      = ''

         ! ReadNum returns a string contained in StrRead.  So, we now try to do an internal read to VarRead and then trap errors.
      read(StringToParse,*,IOSTAT=IOErrStat)   StrRead
      read(StringToParse,*,IOSTAT=IOErrStat)   ValueRead

         ! If IOErrStat==0, then we have a real number, anything else is a problem.
      if (IOErrStat==0) then
         IsRealNum   = .TRUE.
      else
         IsRealNum   = .FALSE.
         ValueRead   = NaN                ! This is NaN as defined in the NWTC_Num.
         ErrMsg      = 'Not a real number. '//TRIM(ErrMsgTmp)//NewLine
         ErrSTat     = ErrID_Severe
      endif

      RETURN
   END SUBROUTINE ReadRealNumberFromString

   !-------------------------------------------------------------------------------
   !> This subroutine works with the ReadNum routine from the library.  ReadNum is
   !! called to read a word from the input file.  An internal read is then done to
   !! convert the string to a number that is stored in VarRead and returned.
   !!
   !! The IsRealNum flag is set to indicate if we actually have a real number or
   !! not.  After calling this routine, a simple if statement can be used:
   !!
   !!       @code
   !!    IF (ISRealNum) THEN
   !!       ! do something
   !!    ELSE
   !!       ! do something else
   !!    ENDIF
   !!       @endcode
   !!
   !-------------------------------------------------------------------------------
   SUBROUTINE ReadRealNumber(UnitNum, FileName, VarName, VarRead, StrRead, IsRealNum, ErrStat, ErrMsg, IOErrStat)
      INTEGER(IntKi),      INTENT(IN   )           :: UnitNum        !< The unit number of the file being read
      CHARACTER(*),        INTENT(IN   )           :: FileName       !< The name of the file being read.  Used in the ErrMsg from ReadNum (Library routine).
      CHARACTER(*),        INTENT(IN   )           :: VarName        !< The variable we are reading.  Used in the ErrMsg from ReadNum (Library routine)'.
      REAL(R8Ki),          INTENT(  OUT)           :: VarRead        !< The variable being read.  Returns as NaN (library defined) if not a Real.
      CHARACTER(*),        INTENT(  OUT)           :: StrRead        !< A string containing what was read from the ReadNum routine.
      LOGICAL,             INTENT(  OUT)           :: IsRealNum      !< Flag indicating if we successfully read a Real
      INTEGER(IntKi),      INTENT(  OUT)           :: ErrStat        !< ErrID level returned from ReadNum
      CHARACTER(*),        INTENT(  OUT)           :: ErrMsg         !< Error message including message from ReadNum
      INTEGER(IntKi),      INTENT(  OUT)           :: IOErrStat      !< Error status from the internal read. Useful for diagnostics.

      INTEGER(IntKi)                      :: ErrStatTmp
      CHARACTER(2048)                     :: ErrMsgTmp

      ErrStat     = ErrID_None
      ErrMsg      = ''

         ! Now call the ReadNum routine to get the number
         ! If it is a word that does not start with T or F, then ReadNum won't give any errors.
      CALL ReadNum( UnitNum, FileName, StrRead, VarName, ErrStatTmp, ErrMsgTmp)

         ! ReadNum returns a string contained in StrRead.  So, we now try to do an internal read to VarRead and then trap errors.
      read(StrRead,*,IOSTAT=IOErrStat)   VarRead

         ! If IOErrStat==0, then we have a real number, anything else is a problem.
      if (IOErrStat==0) then
         IsRealNum   = .TRUE.
      else
         IsRealNum   = .FALSE.
         VarRead     = NaN                ! This is NaN as defined in the NWTC_Num.
         ErrMsg      = 'Not a real number. '//TRIM(ErrMsgTmp)//NewLine
         ErrStat     = ErrStatTmp         ! The ErrStatTmp returned by the ReadNum routine is an ErrID level.
      endif
      RETURN
   END SUBROUTINE ReadRealNumber

END SUBROUTINE ReadInputDispFile




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


SUBROUTINE Dvr_InitializeOutputFile(OutUnit,IntOutput,RootName,ErrStat,ErrMsg)
   integer(IntKi),              intent(  out):: OutUnit
   type(SlD_InitOutputType),    intent(in   ):: IntOutput     ! Output for initialization routine
   integer(IntKi),              intent(  out):: ErrStat     ! Error status of the operation
   character(*),                intent(  out):: ErrMsg      ! Error message if ErrStat /= ErrID_None
   character(*),                intent(in   ):: RootName
   integer(IntKi)                            :: i
   integer(IntKi)                            :: numOuts
   integer(IntKi)                            :: ErrStat2                     ! Temporary Error status
   character(ErrMsgLen)                      :: ErrMsg2                      ! Temporary Error message
   character(*), parameter                   :: RoutineName = 'Dvr_InitializeOutputFile'

   ErrStat = ErrID_none
   ErrMsg  = ""

   CALL GetNewUnit(OutUnit,ErrStat2,ErrMsg2)
   CALL OpenFOutFile ( OutUnit, trim(RootName)//'.out', ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if (ErrStat >= AbortErrLev) return

   write (OutUnit,'(/,A)')  'Predictions were generated on '//CurDate()//' at '//CurTime()//' using '//trim(GetNVD(IntOutput%Ver))
   write (OutUnit,'()' )    !print a blank line

   numOuts = size(IntOutput%WriteOutputHdr)
   !......................................................
   ! Write the names of the output parameters on one line:
   !......................................................

   write (OutUnit,'()')
   write (OutUnit,'()')
   write (OutUnit,'()')

   call WrFileNR ( OutUnit, 'Time' )

   do i=1,NumOuts
      call WrFileNR ( OutUnit, tab//IntOutput%WriteOutputHdr(i) )
   end do ! i

   write (OutUnit,'()')

      !......................................................
      ! Write the units of the output parameters on one line:
      !......................................................

   call WrFileNR ( OutUnit, '(s)' )

   do i=1,NumOuts
      call WrFileNR ( Outunit, tab//trim(IntOutput%WriteOutputUnt(i)) )
   end do ! i

   write (OutUnit,'()')


END SUBROUTINE Dvr_InitializeOutputFile

!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Dvr_WriteOutputLine(t,OutUnit, OutFmt, Output)
   real(DbKi)             ,  intent(in   )   :: t                    ! simulation time (s)
   integer(IntKi)         ,  intent(in   )   :: OutUnit              ! Status of error message
   character(*)           ,  intent(in   )   :: OutFmt
   type(SlD_OutputType),     intent(in   )   :: Output
   integer(IntKi)                            :: errStat              ! Status of error message (we're going to ignore errors in writing to the file)
   character(ErrMsgLen)                      :: errMsg               ! Error message if ErrStat /= ErrID_None
   character(200)                            :: frmt                 ! A string to hold a format specifier
   character(15)                             :: tmpStr               ! temporary string to print the time output as text

   frmt = '"'//tab//'"'//trim(OutFmt)      ! format for array elements from individual modules

      ! time
   write( tmpStr, '(F15.6)' ) t
   call WrFileNR( OutUnit, tmpStr )
   call WrNumAryFileNR ( OutUnit, Output%WriteOutput,  frmt, errStat, errMsg )

     ! write a new line (advance to the next line)
   write (OutUnit,'()')
end subroutine Dvr_WriteOutputLine


END MODULE SoilDyn_Driver_Subs
