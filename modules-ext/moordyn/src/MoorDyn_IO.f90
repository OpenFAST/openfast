!..................................................................................................................................

!    This file is part of MoorDyn.
!

MODULE MoorDyn_IO

  ! This MODULE stores variables used for input and output and provides i/o subs

  USE                              NWTC_Library
  USE                              MoorDyn_Types
  IMPLICIT                         NONE

  INTEGER(IntKi)                 :: UnOutFile      ! unit number of main output file

  !REAL(ReKi),      ALLOCATABLE   :: MDWrOutput(:)  ! one line of output data (duplicate of y%WriteOutput, should fix)

  INTEGER(IntKi),  ALLOCATABLE   :: UnLineOuts(:)  ! Unit numbers of line output files
  REAL(ReKi),      ALLOCATABLE   :: LineWriteOutputs(:) ! one line of output data for Lines


  PRIVATE

  ! --------------------------- Output definitions -----------------------------------------

  ! The following are some definitions for use with the output options in MoorDyn.
  ! These are for the global output quantities specified by OutList, not line-specific outputs.
  ! Output definitions follow the structure described by the MD_OutParmType .
  ! Each output channel is described by the following fields:
  !  Name   - (string) what appears at the top of the output column
  !  Units  - (string) selected from UnitList (see below) based on index QType
  !  OType  - (int) the type of object the output is from. 1=line, 2=connect (0=invalid)
  !  ObjID  - (int) the ID number of the line or connect
  !  QType  - (int) the type of quantity to output.  0=tension, 1=x pos, etc.  see the parameters below
  !  NodeID - (int) the ID number of the node of the output quantity

  ! These are the "OTypes": 0=Connect object, 1=Line Object
  ! (will just use 0 and 1 rather than parameter names)

  ! Indices for computing output channels:  - customized for the MD_OutParmType approach
  ! these are the "QTypes"
  INTEGER, PARAMETER             :: Time      =    0
  INTEGER, PARAMETER             :: PosX      =    1
  INTEGER, PARAMETER             :: PosY      =    2
  INTEGER, PARAMETER             :: PosZ      =    3
  INTEGER, PARAMETER             :: VelX      =    4
  INTEGER, PARAMETER             :: VelY      =    5
  INTEGER, PARAMETER             :: VelZ      =    6
  INTEGER, PARAMETER             :: AccX      =    7
  INTEGER, PARAMETER             :: AccY      =    8
  INTEGER, PARAMETER             :: AccZ      =    9
  INTEGER, PARAMETER             :: Ten      =    10
  INTEGER, PARAMETER             :: FX      =    11
  INTEGER, PARAMETER             :: FY      =    12
  INTEGER, PARAMETER             :: FZ      =    13

  ! List of units corresponding to the quantities parameters for QTypes
  CHARACTER(ChanLen), PARAMETER :: UnitList(0:13) =  (/ &
                               "(s)       ","(m)       ","(m)       ","(m)       ", &
                               "(m/s)     ","(m/s)     ","(m/s)     ", &
                               "(m/s2)    ","(m/s2)    ","(m/s2)    ", &
                               "(N)       ","(N)       ","(N)       ","(N)       " /)

  CHARACTER(28), PARAMETER  :: OutPFmt = "( I4, 3X,A 10,1 X, A10 )"   ! Output format parameter output list.
  CHARACTER(28), PARAMETER  :: OutSFmt = "ES10.3E2"


  ! output naming scheme is as
  ! examples:
  !  FairTen1, AnchTen1
  !  Con1pX
  !  Con3vY (connection 3, y velocity)
  !  L2N4pX (line 2, node 4, x position)

  ! ---------------------------------------------------------------------------------------------------------




   PUBLIC :: MDIO_ReadInput
   PUBLIC :: MDIO_OpenOutput
   PUBLIC :: MDIO_CloseOutput
   PUBLIC :: MDIO_ProcessOutList
   PUBLIC :: MDIO_WriteOutputs


CONTAINS




  !====================================================================================================
  SUBROUTINE MDIO_ReadInput( InitInp, p, other, ErrStat, ErrMsg )

  ! This subroutine reads the input required for MoorDyn from the file whose name is an
  ! input parameter.  It sets the size of p%NTypes, NConnects, and NLines,
  ! allocates LineTypeList, ConnectList, and LineList, and puts all the read contents of
  ! the input file into the respective slots in those lists of types.


    ! Passed variables

    TYPE(MD_InitInputType), INTENT( INOUT )   :: InitInp              ! the MoorDyn data
    TYPE(MD_ParameterType),       INTENT(  OUT)     :: p           ! INTENT( OUT) : Parameters
    TYPE(MD_OtherStateType),      INTENT(  OUT)     :: other       ! INTENT( OUT) : Initial other/optimization states
    INTEGER,                      INTENT(   OUT )   :: ErrStat              ! returns a non-zero value when an error occurs
    CHARACTER(*),                 INTENT(   OUT )   :: ErrMsg               ! Error message if ErrStat /= ErrID_None


    ! Local variables

    INTEGER                      :: I                    ! generic integer for counting
    INTEGER                      :: J                    ! generic integer for counting
    INTEGER                      :: UnIn                 ! Unit number for the input file
    INTEGER                      :: UnEc                 ! The local unit number for this module's echo file
    CHARACTER(1024)              :: EchoFile             ! Name of MoorDyn echo file
    CHARACTER(1024)              :: Line                 ! String to temporarially hold value of read line
    CHARACTER(20)                :: OptString            ! String to temporarially hold name of option variable
    CHARACTER(20)                :: OptValue             ! String to temporarially hold value of options variable input
    CHARACTER(1024)              :: FileName             !

    !
    UnEc = -1

    ! Initialize ErrStat
    ErrStat = ErrID_None
    ErrMsg  = ""

    !-------------------------------------------------------------------------------------------------
    ! Open the file
    !-------------------------------------------------------------------------------------------------
    FileName = TRIM(InitInp%FileName)

    CALL GetNewUnit( UnIn )
    CALL OpenFInpFile( UnIn, FileName, ErrStat )

    IF ( ErrStat /= ErrID_None ) THEN
       ErrMsg  = ' Failed to open MoorDyn input file: '//FileName
       CALL CleanUp()
       RETURN
    END IF


    CALL WrScr( '  MD_Init: Opening MoorDyn input file:  '//FileName )


    !-------------------------------------------------------------------------------------------------
    ! File header
    !-------------------------------------------------------------------------------------------------

    CALL ReadCom( UnIn, FileName, 'MoorDyn input file header line 1', ErrStat )

    IF ( ErrStat /= ErrID_None ) THEN
       ErrMsg  = ' Failed to read MoorDyn input file header line 1.'
       !CALL SetErrStat(ErrID_Fatal, 'Failed to read MoorDyn input file header line 1.',ErrStat,ErrMsg,'MDIO_ReadInput')
       CALL CleanUp()
       RETURN
    END IF


    CALL ReadCom( UnIn, FileName, 'MoorDyn input file header line 2', ErrStat )

    IF ( ErrStat /= ErrID_None ) THEN
       ErrMsg  = ' Failed to read MoorDyn input file header line 2.'
       CALL CleanUp()
       RETURN
    END IF


    ! Echo Input Files.
    CALL ReadVar ( UnIn, FileName, InitInp%Echo, 'Echo', 'Echo Input', ErrStat )
    IF ( ErrStat /= ErrID_None ) THEN
       ErrMsg  = ' Failed to read Echo parameter.'
       CALL CleanUp()
       RETURN
    END IF


    ! If we are Echoing the input then we should re-read the first three lines so that we can echo them
    ! using the NWTC_Library routines.  The echoing is done inside those routines via a global variable
    ! which we must store, set, and then replace on error or completion.

    IF ( InitInp%Echo ) THEN

        EchoFile = TRIM(FileName)//'.ech'                      ! open an echo file for writing
        CALL GetNewUnit( UnEc )
        CALL OpenEcho ( UnEc, EchoFile, ErrStat, ErrMsg )
        IF ( ErrStat /= ErrID_None ) THEN
           ErrMsg  = ' Failed to open Echo file.'
           CALL CleanUp()
           RETURN
        END IF

        REWIND(UnIn)      ! rewind to start of input file to re-read the first few lines

       CALL ReadCom( UnIn, FileName, 'MoorDyn input file header line 1', ErrStat, ErrMsg, UnEc )

       IF ( ErrStat /= ErrID_None ) THEN
          ErrMsg  = ' Failed to read MoorDyn input file header line 1.'
          CALL CleanUp()
          RETURN
       END IF


       CALL ReadCom( UnIn, FileName, 'MoorDyn input file header line 2', ErrStat, ErrMsg, UnEc )

       IF ( ErrStat /= ErrID_None ) THEN
          ErrMsg  = ' Failed to read MoorDyn input file header line 2.'
          CALL CleanUp()
          RETURN
       END IF


       ! Echo Input Files. Note this line is prevented from being echoed by the ReadVar routine.
       CALL ReadVar ( UnIn, FileName, InitInp%Echo, 'Echo', 'Echo the input file data', ErrStat, ErrMsg, UnEc )
       !WRITE (UnEc,Frmt      ) InitInp%Echo, 'Echo', 'Echo input file'
       IF ( ErrStat /= ErrID_None ) THEN
          ErrMsg  = ' Failed to read Echo parameter.'
          CALL CleanUp()
          RETURN
       END IF

    END IF


    !-------------------------------------------------------------------------------------------------
    !  Line Types Properties Section
    !-------------------------------------------------------------------------------------------------

    CALL ReadCom( UnIn, FileName, 'Line types header', ErrStat, ErrMsg, UnEc )

    IF ( ErrStat /= ErrID_None ) THEN
       ErrMsg  = ' Failed to read line types header line.'
       CALL CleanUp()
       RETURN
    END IF


    CALL ReadVar ( UnIn, FileName, p%NTypes, 'NTypes', 'Number of line types', ErrStat, ErrMsg, UnEc )

    IF ( ErrStat /= ErrID_None ) THEN
       ErrMsg  = ' Failed to read NTypes parameter.'
       CALL CleanUp()
       RETURN
    END IF


    ! Table header
    DO I = 1, 2
       CALL ReadCom( UnIn, FileName, 'Line types table header', ErrStat, ErrMsg, UnEc )

       IF ( ErrStat /= ErrID_None ) THEN
          ErrMsg  = ' Failed to read Line types table header line.'
          CALL CleanUp()
          RETURN
       END IF
    END DO

    ! make sure NTypes isn't zero
    IF ( p%NTypes < 1 ) THEN
       ErrMsg  = ' NTypes parameter must be greater than zero.'
       CALL CleanUp()
       RETURN
    END IF

    ! Allocate memory for LineTypeList array to hold line type properties
    ALLOCATE ( other%LineTypeList(p%NTypes), STAT = ErrStat )
    IF ( ErrStat /= ErrID_None ) THEN
       ErrMsg  = ' Error allocating space for LineTypeList array.'
       CALL CleanUp()
       RETURN
    END IF

    ! read each line
    DO I = 1,p%NTypes
          ! read the table entries   Name      Diam    MassDenInAir    EA        cIntDamp     Can     Cat    Cdn     Cdt     in the MoorDyn input file
       READ(UnIn,'(A)',IOSTAT=ErrStat) Line      !read into a line

       IF (ErrStat == 0) THEN
          READ(Line,*,IOSTAT=ErrStat) other%LineTypeList(I)%name, other%LineTypeList(I)%d,  &
             other%LineTypeList(I)%w, other%LineTypeList(I)%EA, other%LineTypeList(I)%BA, &
             other%LineTypeList(I)%Can, other%LineTypeList(I)%Cat, other%LineTypeList(I)%Cdn, other%LineTypeList(I)%Cdt
       END IF

       other%LineTypeList(I)%IdNum = I  ! specify IdNum of line type for error checking


       IF ( ErrStat /= ErrID_None ) THEN
          ErrMsg  = ' Failed to read line type properties for line '//trim(Num2LStr(I))
          CALL CleanUp()
          RETURN
       END IF

       IF ( InitInp%Echo ) THEN
          WRITE( UnEc, '(A)' ) TRIM(Line)
       END IF

    END DO



    !-------------------------------------------------------------------------------------------------
    !  Connections Section
    !-------------------------------------------------------------------------------------------------

    CALL ReadCom( UnIn, FileName, 'Connections header', ErrStat, ErrMsg, UnEc )

    IF ( ErrStat /= ErrID_None ) THEN
       ErrMsg  = ' Failed to read connections header line.'
       CALL CleanUp()
       RETURN
    END IF


    CALL ReadVar ( UnIn, FileName, p%NConnects, 'NConnects', 'Number of Connects', ErrStat, ErrMsg, UnEc )

    IF ( ErrStat /= ErrID_None ) THEN
       ErrMsg  = ' Failed to read NConnects parameter.'
       CALL CleanUp()
       RETURN
    END IF


    ! Table header
    DO I = 1, 2
       CALL ReadCom( UnIn, FileName, 'Connects header', ErrStat, ErrMsg, UnEc )

       IF ( ErrStat /= ErrID_None ) THEN
          ErrMsg  = ' Failed to read Connects table header line.'
          CALL CleanUp()
          RETURN
       END IF
    END DO

    ! make sure NConnects is at least two
    IF ( p%NConnects < 2 ) THEN
       ErrMsg  = ' NConnects parameter must be at least 2.'
       CALL CleanUp()
       RETURN
    END IF

     ! allocate ConnectList
    ALLOCATE ( other%ConnectList(p%NConnects), STAT = ErrStat )
    IF ( ErrStat /= ErrID_None ) THEN
       ErrMsg  = ' Error allocating space for ConnectList array.'
       CALL CleanUp()
       RETURN
    END IF

    ! read each line
    DO I = 1,p%NConnects
          ! read the table entries   Node      Type      X        Y         Z        M        V        FX       FY      FZ
       READ(UnIn,'(A)',IOSTAT=ErrStat) Line      !read into a line

       IF (ErrStat == 0) THEN
          READ(Line,*,IOSTAT=ErrStat) other%ConnectList(I)%IdNum, other%ConnectList(I)%type, other%ConnectList(I)%conX, &
               other%ConnectList(I)%conY, other%ConnectList(I)%conZ, other%ConnectList(I)%conM, &
               other%ConnectList(I)%conV, other%ConnectList(I)%conFX, other%ConnectList(I)%conFY, &
                other%ConnectList(I)%conFZ
       END IF

       ! check for sequential IdNums
       IF ( other%ConnectList(I)%IdNum .NE. I ) THEN
         ErrMsg  = ' Node numbers must be sequential starting from 1.'
         CALL CleanUp()
         RETURN
       END IF


       IF ( ErrStat /= ErrID_None ) THEN
          ErrMsg  = ' Failed to read connects.'  ! would be nice to specify which line
          CALL CleanUp()
          RETURN
       END IF

       IF ( InitInp%Echo ) THEN
          WRITE( UnEc, '(A)' ) TRIM(Line)
       END IF

    END DO




    !-------------------------------------------------------------------------------------------------
    !  Lines Section
    !-------------------------------------------------------------------------------------------------

    CALL ReadCom( UnIn, FileName, 'Lines header', ErrStat, ErrMsg, UnEc )

    IF ( ErrStat /= ErrID_None ) THEN
       ErrMsg  = ' Failed to read lines header line.'
       CALL CleanUp()
       RETURN
    END IF


    CALL ReadVar ( UnIn, FileName, p%NLines, 'NLines', 'Number of Lines', ErrStat, ErrMsg, UnEc )

    IF ( ErrStat /= ErrID_None ) THEN
       ErrMsg  = ' Failed to read NLines parameter.'
       CALL CleanUp()
       RETURN
    END IF


    ! Table header
    DO I = 1, 2
       CALL ReadCom( UnIn, FileName, 'Lines header', ErrStat, ErrMsg, UnEc )

       IF ( ErrStat /= ErrID_None ) THEN
          ErrMsg  = ' Failed to read Lines table header line.'
          CALL CleanUp()
          RETURN
       END IF
    END DO

    ! make sure NLines is at least one
    IF ( p%NLines < 1 ) THEN
       ErrMsg  = ' NLines parameter must be at least 1.'
       CALL CleanUp()
       RETURN
    END IF

     ! allocate LineList
    ALLOCATE ( other%LineList(p%NLines), STAT = ErrStat )
    IF ( ErrStat /= ErrID_None ) THEN
       ErrMsg  = ' Error allocating space for LineList array.'
       CALL CleanUp()
       RETURN
    END IF

    ! read each line
    DO I = 1,p%NLines
          ! read the table entries   Line     LineType  UnstrLen  NumSegs   NodeAnch  NodeFair  Flags/Outputs
       READ(UnIn,'(A)',IOSTAT=ErrStat) Line      !read into a line

       IF (ErrStat == 0) THEN
          READ(Line,*,IOSTAT=ErrStat) other%LineList(I)%IdNum, other%LineList(I)%type, other%LineList(I)%UnstrLen, &
            other%LineList(I)%N, other%LineList(I)%AnchConnect, other%LineList(I)%FairConnect, other%LineList(I)%OutFlags
       END IF

       ! check for sequential IdNums
       IF ( other%LineList(I)%IdNum .NE. I ) THEN
         ErrMsg  = ' Line numbers must be sequential starting from 1.'
         CALL CleanUp()
         RETURN
       END IF

       ! identify index of line type
       DO J = 1,p%NTypes
         IF (trim(other%LineList(I)%type) == trim(other%LineTypeList(J)%name)) THEN
           other%LineList(I)%PropsIdNum = J
           EXIT
           IF (J == p%NTypes) THEN   ! call an error if there is no match
             ErrMsg = 'Unable to find matching line type name for Line '//trim(Num2LStr(I))
             ErrStat = ErrID_Severe
           END IF
         END IF
       END DO

       IF ( ErrStat /= ErrID_None ) THEN
          ErrMsg  = ' Failed to read line data for Line '//trim(Num2LStr(I))  ! would be nice to specify which line
          CALL CleanUp()
          RETURN
       END IF

       IF ( InitInp%Echo ) THEN
          WRITE( UnEc, '(A)' ) TRIM(Line)
       END IF

    END DO  ! I


    !-------------------------------------------------------------------------------------------------
    ! Read any options lines
    !-------------------------------------------------------------------------------------------------

    CALL ReadCom( UnIn, FileName, 'Options header', ErrStat, ErrMsg, UnEc )

    IF ( ErrStat /= ErrID_None ) THEN
       ErrMsg  = ' Failed to read options header line.'
       CALL CleanUp()
       RETURN
    END IF

     ! loop through any remaining input lines, and use them to set options (overwriting default values in many cases).
     ! doing this manually since I'm not sure that there is a built in subroutine for reading any input value on any line number.
     DO

       READ(UnIn,'(A)',IOSTAT=ErrStat) Line      !read into a line

       IF (ErrStat == 0) THEN
         IF (( Line(1:3) == '---' ) .OR. ( Line(1:3) == 'END' ) .OR. ( Line(1:3) == 'end' ))  EXIT  ! check if it's the end line

         READ(Line,*,IOSTAT=ErrStat) OptValue, OptString  ! look at first two entries, ignore remaining words in line, which should be comments
       END IF

       CALL Conv2UC(OptString)

       ! check all possible options types and see if OptString is one of them, in which case set the variable.
       if ( OptString == 'DT') THEN
         read (OptValue,*) InitInp%DTmooring   
       else if ( OptString == 'G') then
         read (OptValue,*) p%G
       else if ( OptString == 'RHOW') then
         read (OptValue,*) p%rhoW
       else if ( OptString == 'WTRDPTH') then
         read (OptValue,*) p%WtrDpth
       else if ( OptString == 'KBOT')  then
         read (OptValue,*) p%kBot
       else if ( OptString == 'CBOT')  then
         read (OptValue,*) p%cBot
       else if ( OptString == 'DTIC')  then
         read (OptValue,*) InitInp%dtIC
       else if ( OptString == 'TMAXIC')  then
         read (OptValue,*) InitInp%TMaxIC
       else if ( OptString == 'CDSCALEIC')  then
         read (OptValue,*) InitInp%CdScaleIC
       else if ( OptString == 'THRESHIC')  then
         read (OptValue,*) InitInp%threshIC
       else
         ErrStat = ErrID_Warn
 print *, 'unable to interpret input ', OptString
       end if

       IF ( ErrStat > ErrID_Warn ) THEN
          ErrMsg  = ' Failed to read options.'  ! would be nice to specify which line had the error
          CALL CleanUp()
          RETURN
       END IF

       IF ( InitInp%Echo ) THEN
          WRITE( UnEc, '(A)' ) TRIM(Line)
       END IF

     END DO


     !-------------------------------------------------------------------------------------------------
     ! Read the FAST-style outputs list in the final section, if there is one
     !-------------------------------------------------------------------------------------------------
!     we don't read in the outputs header line because it's already been read in for detecting the end of the variable-length options section
!     CALL ReadCom( UnIn, FileName, 'Outputs header', ErrStat, ErrMsg, UnEc )

    ! allocate InitInp%Outliest (to a really big number for now...)
    CALL AllocAry( InitInp%OutList, 1000, "MoorDyn Input File's Outlist", ErrStat, ErrMsg )
    !  CALL CheckError( ErrStat2, ErrMsg2 )
    !  IF ( ErrStat >= AbortErrLev ) RETURN

    ! OutList - List of user-requested output channels (-):
    CALL ReadOutputList ( UnIn, FileName, InitInp%OutList, p%NumOuts, 'OutList', "List of user-requested output channels", ErrStat, ErrMsg, UnEc  )     ! Routine in NWTC Subroutine Library
   !   CALL CheckError( ErrStat2, ErrMsg2 )
   !   IF ( ErrStat >= AbortErrLev ) RETURN
    IF ( ErrStat > ErrID_Warn ) THEN
       ErrMsg  = ' Failed to read output list.'
       CALL CleanUp()
       RETURN
    END IF

  !print *, 'NumOuts is ', p%NumOuts
  print *, '  OutList is ', InitInp%OutList(1:p%NumOuts)


     !-------------------------------------------------------------------------------------------------
     ! This is the end of the input file
     !-------------------------------------------------------------------------------------------------

           CLOSE( UnIn )
        IF (InitInp%Echo) CLOSE( UnEc )


  CONTAINS
     ! subroutine to set ErrState and close the files if an error occurs
     SUBROUTINE CleanUp()

        ErrStat = ErrID_Fatal
        CLOSE( UnIn )
        IF (InitInp%Echo) CLOSE( UnEc )

     END SUBROUTINE

  END SUBROUTINE MDIO_ReadInput
  ! ====================================================================================================



  ! ====================================================================================================
  SUBROUTINE MDIO_ProcessOutList(OutList, p, other, y, ErrStat, ErrMsg )

  ! This routine processes the output channels requested by OutList, checking for validity and setting
  ! the p%OutParam structures (of type MD_OutParmType) for each valid output.
  ! It assumes the value p%NumOuts has been set beforehand, and sets the values of p%OutParam.


    IMPLICIT                        NONE

    ! Passed variables
    CHARACTER(ChanLen),        INTENT(IN)     :: OutList(:)                  ! The list out user-requested outputs
    TYPE(MD_ParameterType),    INTENT(INOUT)  :: p                           ! The module parameters
    TYPE(MD_OtherStateType),   INTENT(IN)     :: other
    TYPE(MD_OutputType),       INTENT(INOUT)  :: y                           ! Initial system outputs (outputs are not calculated; only the output mesh is initialized)
    INTEGER(IntKi),            INTENT(OUT)    :: ErrStat                     ! The error status code
    CHARACTER(*),              INTENT(OUT)    :: ErrMsg                      ! The error message, if an error occurred

    ! Local variables
    INTEGER                      :: I                                        ! Generic loop-counting index
    INTEGER                      :: J                                        ! Generic loop-counting index
    INTEGER                      :: INDX                                     ! Index for valid arrays

    CHARACTER(ChanLen)           :: OutListTmp                               ! A string to temporarily hold OutList(I), the name of each output channel
    CHARACTER(ChanLen)           :: qVal                                     ! quantity type string to match to list of valid options

    INTEGER                      :: oID                                      ! ID number of connect or line object
    INTEGER                      :: nID                                      ! ID number of node object
    INTEGER                      :: i1,i2,i3,i4                              ! indices of start of numbers or letters in OutListTmp string, for parsing


    ! see the top of the module for info on the output labelling types

    ! Initialize values
    ErrStat = ErrID_None
    ErrMsg = ""


    ALLOCATE ( p%OutParam(1:p%NumOuts) , STAT=ErrStat )   ! note: I'm skipping the time output entry at index 0 for simplicity
    IF ( ErrStat /= 0_IntKi )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = "Error allocating memory for the MoorDyn OutParam array."
      RETURN
    ELSE
      ErrStat = ErrID_None
    ENDIF


    ! Set index, name, and units for the time output channel: ! note: I'm skipping the time output entry at index 0
    !p%OutParam(0)%Indx  = Time
    !p%OutParam(0)%Name  = "Time"    ! OutParam(0) is the time channel by default.
    !p%OutParam(0)%Units = "(s)"
    !p%OutParam(0)%SignM = 1


    ! Set index, name, and units for all of the output channels.
    ! If a selected output channel is not valid set ErrStat = ErrID_Warn.


    ! go through list of requested output names and process (this is a bit of a mess)

    DO I = 1,p%NumOuts

      OutListTmp          = OutList(I)  ! current requested output name
      !p%OutParam(I)%Name  = OutListTmp
      CALL Conv2UC(OutListTmp)       ! convert to all uppercase for string matching purposes

      ! find indicies of changes in number-vs-letter in characters of OutListTmp
      i1 = scan( OutListTmp , '1234567890' )        ! first number in the string
      i2 = i1+verify( OutListTmp(i1+1:) , '1234567890' ) ! second letter start (assuming first character is a letter, i.e. i1>1)
      i3 = i2+scan( OutListTmp(i2+1:) , '1234567890' )   ! second number start
      i4 = i3+verify( OutListTmp(i3+1:) , '1234567890' ) ! third letter start
      !i5 = scan( OutListTmp(i1:) , '1234567890' )   ! find first letter after first number

      ! error check
      IF (i1 <= 1) THEN
        CALL DenoteInvalidOutput(p%OutParam(I)) ! flag as invalid
        CALL WrScr('Warning: invalid output specifier.')
        CONTINUE
      END IF

        p%OutParam(I)%Name = OutListTmp  ! label channel with whatever name was inputted, for now
!    print *, 'processing output channel request ', OutListTmp

!    print *, 'indices are ', i1, i2, i3, i4

      ! figure out what type of output it is and process accordingly

      ! fairlead tension case
      IF (OutListTmp(1:i1-1) == 'FAIRTEN') THEN
        p%OutParam(I)%OType = 1                ! line object type
        p%OutParam(I)%QType = Ten              ! tension quantity type
        p%OutParam(I)%Units = UnitList(Ten)    ! set units according to QType
        READ (OutListTmp(i1:),*) oID
        p%OutParam(I)%ObjID =  oID
        p%OutParam(I)%NodeID =  other%LineList(oID)%N  ! line type

      ! more general case
      ELSE

        ! what object type?
        ! Line case                                          ... L?N?xxxx
        IF (OutListTmp(1:i1-1) == 'L') THEN
          p%OutParam(I)%OType = 1                ! Line object type
          ! for now we'll just assume the next character(s) are "n" to represent node number:
          READ (OutListTmp(i3:i4-1),*) nID
          p%OutParam%NodeID = nID
          qVal = OutListTmp(i4:)  ! isolate quantity type string
        ! Connect case                                     ... C?xxx or Con?xxx
        ELSE IF (OutListTmp(1:1) == 'C') THEN
          p%OutParam(I)%OType = 2                ! Connect object type
          qVal = OutListTmp(i2:)  ! isolate quantity type string

        ! should do fairlead option also!

        ! error
        ELSE
          CALL DenoteInvalidOutput(p%OutParam(I)) ! flag as invalid
          CALL WrScr('Warning: invalid output specifier - type must be L or C.')
          CONTINUE
        END IF

        ! object number
        READ (OutListTmp(i1:i2-1),*) oID
        p%OutParam(I)%ObjID =  oID             ! line or connect ID number

        ! which kind of quantity?
        IF (qVal == 'PX') THEN
          p%OutParam(I)%QType = PosX
          p%OutParam(I)%Units = UnitList(PosX)
        ELSE IF (qVal == 'PY') THEN
          p%OutParam(I)%QType = PosY
          p%OutParam(I)%Units = UnitList(PosY)
        ELSE IF (qVal == 'PZ') THEN
          p%OutParam(I)%QType = PosZ
          p%OutParam(I)%Units = UnitList(PosZ)
        ELSE IF (qVal == 'VX') THEN
          p%OutParam(I)%QType = VelX
          p%OutParam(I)%Units = UnitList(VelX)
        ELSE IF (qVal == 'VY') THEN
          p%OutParam(I)%QType = VelY
          p%OutParam(I)%Units = UnitList(VelY)
        ELSE IF (qVal == 'VZ') THEN
          p%OutParam(I)%QType = VelZ
          p%OutParam(I)%Units = UnitList(VelZ)
        ELSE IF (qVal == 'AX') THEN
          p%OutParam(I)%QType = AccX
          p%OutParam(I)%Units = UnitList(AccX)
        ELSE IF (qVal == 'Ay') THEN
          p%OutParam(I)%QType = AccY
          p%OutParam(I)%Units = UnitList(AccY)
        ELSE IF (qVal == 'AZ') THEN
          p%OutParam(I)%QType = AccZ
          p%OutParam(I)%Units = UnitList(AccZ)
        ELSE IF ((qVal == 'T') .or. (qval == 'Ten')) THEN
          p%OutParam(I)%QType = Ten
          p%OutParam(I)%Units = UnitList(Ten)
        ELSE IF (qVal == 'FX') THEN
          p%OutParam(I)%QType = FX
          p%OutParam(I)%Units = UnitList(FX)
        ELSE IF (qVal == 'FY') THEN
          p%OutParam(I)%QType = FY
          p%OutParam(I)%Units = UnitList(FY)
        ELSE IF (qVal == 'FZ') THEN
          p%OutParam(I)%QType = FZ
          p%OutParam(I)%Units = UnitList(FZ)
        ELSE
          CALL DenoteInvalidOutput(p%OutParam(I)) ! flag as invalid
          CALL WrScr('Warning: invalid output specifier - quantity type not recognized.')  ! need to figure out how to add numbers/strings to these warning messages...
          CONTINUE
        END IF

      END IF

      ! also check whether each object index and node index (if applicable) is in range
      IF (p%OutParam(I)%OType==2) THEN
        IF (p%OutParam(I)%ObjID > p%NConnects) THEN
          print *, 'warning, output Connect index excedes number of Connects'
          CALL DenoteInvalidOutput(p%OutParam(I)) ! flag as invalid
        END IF
      ELSE IF (p%OutParam(I)%OType==1) THEN
        IF (p%OutParam(I)%ObjID > p%NLines) THEN
          print *, 'warning, output Line index excedes number of Line'
          CALL DenoteInvalidOutput(p%OutParam(I)) ! flag as invalid
        END IF
        IF (p%OutParam(I)%NodeID > other%LineList(p%OutParam(I)%ObjID)%N) THEN
          print *, 'warning, output node index excedes number of nodes'
          CALL DenoteInvalidOutput(p%OutParam(I)) ! flag as invalid
        ELSE IF (p%OutParam(I)%NodeID < 0) THEN
          print *, 'warning, output node index is less than zero'
          CALL DenoteInvalidOutput(p%OutParam(I)) ! flag as invalid
        END IF

      END IF

      ! is the reverse sign functionality necessary?
      !      ! Reverse the sign (+/-) of the output channel if the user prefixed the
      !      !   channel name with a "-", "_", "m", or "M" character indicating "minus".

     END DO  ! I ... looping through OutList


!!   ! Allocate MDWrOutput which is used to store a time step's worth of output channels, prior to writing to a file.
!    ALLOCATE( MDWrOutput( p%NumOuts),  STAT = ErrStat )
!    IF ( ErrStat /= ErrID_None ) THEN
!      ErrMsg  = ' Error allocating space for MDWrOutput array.'
!      ErrStat = ErrID_Fatal
!      RETURN
!    END IF


    ! Allocate MDWrOuput2 which is used to store a time step's worth of output data for each line, just making it really big for now <<<<<<<<<<<<<<
    ALLOCATE( LineWriteOutputs( 200),  STAT = ErrStat )
    IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Error allocating space for LineWriteOutputs array.'
      ErrStat = ErrID_Fatal
      RETURN
    END IF

    !Allocate WriteOuput
    ALLOCATE( y%WriteOutput( p%NumOuts),  STAT = ErrStat )
    IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Error allocating space for y%WriteOutput array.'
      ErrStat = ErrID_Fatal
      RETURN
    END IF

    !print *, "y%WriteOutput allocated to size ", size(y%WriteOutput)

   ! These variables are to help follow the framework template, but the data in them is simply a copy of data
   ! already available in the OutParam data structure
  !  ALLOCATE ( InitOut%WriteOutputHdr(p%NumOuts+p%OutAllint*p%OutAllDims), STAT = ErrStat )
  !  ALLOCATE ( InitOut%WriteOutputUnt(p%NumOuts+p%OutAllint*p%OutAllDims), STAT = ErrStat )

  !   DO I = 1,p%NumOuts+p%OutAllint*p%OutAllDims
  !    InitOut%WriteOutputHdr(I) = TRIM( p%OutParam(I)%Name  )
  !    InitOut%WriteOutputUnt(I) = TRIM( p%OutParam(I)%Units )
  !   END DO


  CONTAINS

    SUBROUTINE DenoteInvalidOutput( OutParm )
      TYPE(MD_OutParmType), INTENT (INOUT)  :: OutParm

      OutParm%OType = 0  ! flag as invalid
      OutParm%Name = 'Invalid'
      OutParm%Units = ' - '

    END SUBROUTINE DenoteInvalidOutput

  END SUBROUTINE MDIO_ProcessOutList
  !---------------------------------------------------------------------------------------------------





   !====================================================================================================
   SUBROUTINE MDIO_OpenOutput( OutRootName,  p, other, InitOut, ErrStat, ErrMsg )
   !----------------------------------------------------------------------------------------------------
   
      CHARACTER(*),                  INTENT( IN    ) :: OutRootName          ! Root name for the output file
      TYPE(MD_ParameterType),        INTENT( INOUT ) :: p
      TYPE(MD_OtherStateType),       INTENT (IN)     :: other
      TYPE(MD_InitOutPutType ),      INTENT( IN    ) :: InitOut              !
      INTEGER,                       INTENT(   OUT ) :: ErrStat              ! a non-zero value indicates an error occurred
      CHARACTER(*),                  INTENT(   OUT ) :: ErrMsg               ! Error message if ErrStat /= ErrID_None
   
      INTEGER                                        :: I                    ! Generic loop counter
      INTEGER                                        :: J                    ! Generic loop counter
      CHARACTER(1024)                                :: OutFileName          ! The name of the output file  including the full path.
      CHARACTER(200)                                 :: Frmt                 ! a string to hold a format statement
      INTEGER                                        :: ErrStat2
   
     
      ErrStat = ErrID_None
      ErrMsg  = ""
   
       p%Delim = ' '  ! for now
   
      !-------------------------------------------------------------------------------------------------
      ! Open the output file, if necessary, and write the header
      !-------------------------------------------------------------------------------------------------
   
      IF ( ALLOCATED( p%OutParam ) .AND. p%NumOuts > 0 ) THEN           ! Output has been requested so let's open an output file
   
            ! Open the file for output
         OutFileName = 'MoorDyn.out'
         CALL GetNewUnit( UnOutFile )
   
         CALL OpenFOutFile ( UnOutFile, OutFileName, ErrStat, ErrMsg )
         IF ( ErrStat >= AbortErrLev ) THEN
            ErrMsg = ' Error opening MoorDyn-level output file: '//TRIM(ErrMsg)
            RETURN
         END IF
   
   
         !Write the names of the output parameters:
   
         Frmt = '(A10,'//TRIM(Int2LStr(p%NumOuts))//'(A1,A10))'
   
         WRITE(UnOutFile,Frmt, IOSTAT=ErrStat2)  TRIM( 'Time' ), ( p%Delim, TRIM( p%OutParam(I)%Name), I=1,p%NumOuts )
   
         WRITE(UnOutFile,Frmt)  TRIM( '(s)' ), ( p%Delim, TRIM( p%OutParam(I)%Units ), I=1,p%NumOuts )
   
      ELSE  ! if no outputs requested
   
         print *, 'note, MDIO_OpenOutput thinks that no outputs have been requested.'
   
      END IF
   
      !--------------------------------------------------------------------------
      ! -------------- now do the same for line output files --------------------
   
      ! allocate UnLineOuts
      ALLOCATE(UnLineOuts(p%NLines))  ! should add error checking
   
      DO I = 1,p%NLines
   
         ! Open the file for output
         OutFileName = 'Line'//TRIM(Int2LStr(I))//'.out'
         CALL GetNewUnit( UnLineOuts(I) )
   
         CALL OpenFOutFile ( UnLineOuts(I), OutFileName, ErrStat, ErrMsg )
         IF ( ErrStat >= AbortErrLev ) THEN
            ErrMsg = ' Error opening Line output file '//TRIM(ErrMsg)
            RETURN
         END IF   
   
         !Write the names of the output parameters:
   
         Frmt = '(A10,'//TRIM(Int2LStr(3+3*other%LineList(I)%N))//'(A1,A10))'
      
         WRITE(UnLineOuts(I),Frmt, IOSTAT=ErrStat2)  TRIM( 'Time' ), ( p%Delim, 'Node'//TRIM(Int2Lstr(J))//'px', p%Delim, 'Node'//TRIM(Int2Lstr(J))//'py', p%Delim, 'Node'//TRIM(Int2Lstr(J))//'pz', J=0,(other%LineList(I)%N) )
   
      END DO ! I
   
   
   END SUBROUTINE MDIO_OpenOutput
   
   !====================================================================================================


   !====================================================================================================
   SUBROUTINE MDIO_CloseOutput ( p, ErrStat, ErrMsg )
      ! This function cleans up after running the MoorDyn output module.
      ! It closes the output file and releases memory.

      TYPE(MD_ParameterType),  INTENT( INOUT )  :: p                    ! data for this instance of the floating platform module
      INTEGER,                       INTENT(   OUT ) :: ErrStat              ! a non-zero value indicates an error occurred
      CHARACTER(*),                  INTENT(   OUT ) :: ErrMsg               ! Error message if ErrStat /= ErrID_None


      INTEGER(IntKi)       :: I  ! generic counter


      ErrStat = 0
      ErrMsg  = ""


      ! close main MoorDyn output file
      CLOSE( UnOutFile, IOSTAT = ErrStat )
         IF ( ErrStat /= 0 ) THEN
            ErrMsg = 'Error closing output file'
            CALL WrScr(ErrMsg)
         END IF

      ! close individual line output files
      DO I=1,p%NLines
         CLOSE( UnLineOuts(I), IOSTAT = ErrStat )
            IF ( ErrStat /= 0 ) THEN
               ErrMsg = 'Error closing line output file'
               CALL WrScr(ErrMsg)
            END IF
      END DO

      ! deallocate output arrays
      IF (ALLOCATED(UnLineOuts)) THEN
         DEALLOCATE(UnLineOuts)
      ENDIF
      IF (ALLOCATED(LineWriteOutputs)) THEN
         DEALLOCATE(LineWriteOutputs)
      ENDIF


   END SUBROUTINE MDIO_CloseOutput
   !====================================================================================================


   !====================================================================================================
   SUBROUTINE MDIO_WriteOutputs( Time, p, other, y, ErrStat, ErrMsg )
      ! This subroutine gathers the output data defined by the OutParams list and
      ! writes it to the output file opened in MDIO_OutInit()
   
      REAL(DbKi),                   INTENT( IN    ) :: Time                 ! Time for this output
      TYPE(MD_ParameterType),       INTENT( IN    ) :: p                    ! MoorDyn module's parameter data
      TYPE(MD_OutputType),          INTENT(INOUT)  :: y           ! INTENT( OUT) : Initial system outputs (outputs are not calculated; only the output mesh is initialized)
      TYPE(MD_OtherStateType),      INTENT( IN    ) :: other                ! MoorDyn module's other data
      INTEGER,                      INTENT(   OUT ) :: ErrStat              ! returns a non-zero value when an error occurs
      CHARACTER(*),                 INTENT(   OUT ) :: ErrMsg               ! Error message if ErrStat /= ErrID_None
    
      INTEGER                                :: I                           ! Generic loop counter
      INTEGER                                :: J                           ! Generic loop counter
      INTEGER                                :: K                           ! Generic loop counter
      CHARACTER(200)                         :: Frmt                        ! a string to hold a format statement
    
    
      IF ( .NOT. ALLOCATED( p%OutParam ) .OR. UnOutFile < 0 )  THEN
         ErrStat = ErrID_Fatal
         ErrMsg  = ' To write outputs for MoorDyn there must be a valid file ID and OutParam must be allocated.'
         RETURN
      ELSE
         ErrStat = ErrID_None
      END IF
    
      ! gather the required output quantities (INCOMPLETE!)
      DO I = 1,p%NumOuts
    
        IF (p%OutParam(I)%OType == 2) THEN  ! if dealing with a Connect output
          SELECT CASE (p%OutParam(I)%QType)
            CASE (PosX)
              y%WriteOutput(I) = other%ConnectList(p%OutParam(I)%ObjID)%r(1)  ! x position
            CASE (PosY)
              y%WriteOutput(I) = other%ConnectList(p%OutParam(I)%ObjID)%r(2) ! y position
            CASE (PosZ)
              y%WriteOutput(I) = other%ConnectList(p%OutParam(I)%ObjID)%r(3) ! z position
            CASE (VelX)
              y%WriteOutput(I) = other%ConnectList(p%OutParam(I)%ObjID)%rd(1) ! x velocity
            CASE (VelY)
              y%WriteOutput(I) = other%ConnectList(p%OutParam(I)%ObjID)%rd(2) ! y velocity
            CASE (VelZ)
              y%WriteOutput(I) = other%ConnectList(p%OutParam(I)%ObjID)%rd(3) ! z velocity
            CASE DEFAULT
              y%WriteOutput(I) = 0.0_ReKi
              ErrStat = ErrID_Warn
              ErrMsg = ' Unsupported output quantity from Connect object requested.'
          END SELECT
    
        ELSE IF (p%OutParam(I)%OType == 1) THEN  ! if dealing with a Line output
    
          SELECT CASE (p%OutParam(I)%QType)
            CASE (PosX)
              y%WriteOutput(I) = other%LineList(p%OutParam(I)%ObjID)%r(1,p%OutParam(I)%NodeID)  ! x position
            CASE (PosY)
              y%WriteOutput(I) = other%LineList(p%OutParam(I)%ObjID)%r(2,p%OutParam(I)%NodeID) ! y position
            CASE (PosZ)
              y%WriteOutput(I) = other%LineList(p%OutParam(I)%ObjID)%r(3,p%OutParam(I)%NodeID) ! z position
            CASE (VelX)
              y%WriteOutput(I) = other%LineList(p%OutParam(I)%ObjID)%rd(1,p%OutParam(I)%NodeID) ! x velocity
            CASE (VelY)
              y%WriteOutput(I) = other%LineList(p%OutParam(I)%ObjID)%rd(2,p%OutParam(I)%NodeID) ! y velocity
            CASE (VelZ)
              y%WriteOutput(I) = other%LineList(p%OutParam(I)%ObjID)%rd(3,p%OutParam(I)%NodeID) ! z velocity
            CASE (Ten)
              y%WriteOutput(I) = TwoNorm(other%LineList(p%OutParam(I)%ObjID)%T(:,p%OutParam(I)%NodeID))  ! tension this isn't quite right, since it's segment tension...
            CASE DEFAULT
              y%WriteOutput(I) = 0.0_ReKi
              ErrStat = ErrID_Warn
              ErrMsg = ' Unsupported output quantity from Line object requested.'
          END SELECT
    
        ELSE  ! it must be an invalid output, so write zero
          y%WriteOutput(I) = 0.0_ReKi
    
        END IF
    
      END DO ! I, loop through OutParam
  
  
      ! Write the output parameters to the file
  
      Frmt = '(F10.4,'//TRIM(Int2LStr(p%NumOuts))//'(A1,e10.4))'   ! should evenutally use user specified format?
    
      WRITE(UnOutFile,Frmt)  Time, ( p%Delim, y%WriteOutput(I), I=1,p%NumOuts )
  
  
  
  
  
      !------------------------------------------------------------------------
      ! now do the outputs for each line!  <<< so far this is just writing node positions without any user options
   
      DO I=1,p%NLines
        Frmt = '(F10.4,'//TRIM(Int2LStr(3+3*other%LineList(I)%N))//'(A1,e10.4))'   ! should evenutally use user specified format?
   
        DO J = 0,other%LineList(I)%N  ! note index starts at zero because these are nodes
          DO K = 1,3
            LineWriteOutputs(3*J+K) = other%LineList(I)%r(K,J)
          END DO
        END DO
   
        WRITE(UnLineOuts(I),Frmt)  Time, ( p%Delim, LineWriteOutputs(J), J=1,(3+3*other%LineList(I)%N) )
   
      END DO ! I
  
   END SUBROUTINE MDIO_WriteOutputs
   !====================================================================================================


END MODULE MoorDyn_IO
