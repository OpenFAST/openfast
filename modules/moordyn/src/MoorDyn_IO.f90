!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2015  Matthew Hall
!
!    This file is part of MoorDyn.
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
MODULE MoorDyn_IO

  ! This MODULE stores variables used for input and output and provides i/o subs

  USE                              NWTC_Library
  USE                              MoorDyn_Types
  IMPLICIT                         NONE


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
   SUBROUTINE MDIO_ReadInput( InitInp, p, m, ErrStat, ErrMsg )

   ! This subroutine reads the input required for MoorDyn from the file whose name is an
   ! input parameter.  It sets the size of p%NTypes, NConnects, and NLines,
   ! allocates LineTypeList, ConnectList, and LineList, and puts all the read contents of
   ! the input file into the respective slots in those lists of types.


    ! Passed variables

    TYPE(MD_InitInputType),       INTENT( INOUT )   :: InitInp              ! the MoorDyn data
    TYPE(MD_ParameterType),       INTENT(INOUT)     :: p                    ! Parameters
    TYPE(MD_MiscVarType),         INTENT(  OUT)     :: m                    ! INTENT( OUT) : Initial misc/optimization vars
    INTEGER,                      INTENT(   OUT )   :: ErrStat              ! returns a non-zero value when an error occurs
    CHARACTER(*),                 INTENT(   OUT )   :: ErrMsg               ! Error message if ErrStat /= ErrID_None


    ! Local variables

    INTEGER                      :: I                    ! generic integer for counting
    INTEGER                      :: J                    ! generic integer for counting
    INTEGER                      :: UnIn                 ! Unit number for the input file
    INTEGER                      :: UnEc                 ! The local unit number for this module's echo file
    CHARACTER(1024)              :: EchoFile             ! Name of MoorDyn echo file
    CHARACTER(1024)              :: Line                 ! String to temporarially hold value of read line
    CHARACTER(20)                :: LineOutString        ! String to temporarially hold characters specifying line output options
    CHARACTER(20)                :: OptString            ! String to temporarially hold name of option variable
    CHARACTER(20)                :: OptValue             ! String to temporarially hold value of options variable input
    CHARACTER(1024)              :: FileName             !

    INTEGER(IntKi)               :: ErrStat2
    CHARACTER(ErrMsgLen)         :: ErrMsg2
    CHARACTER(*), PARAMETER      :: RoutineName = 'MDIO_ReadInput'
    
    
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
    CALL OpenFInpFile( UnIn, FileName, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       IF ( ErrStat >= AbortErrLev ) THEN
          CALL CleanUp()
          RETURN
       END IF


    CALL WrScr( '  MD_Init: Opening MoorDyn input file:  '//FileName )


    !-------------------------------------------------------------------------------------------------
    ! File header
    !-------------------------------------------------------------------------------------------------

    CALL ReadCom( UnIn, FileName, 'MoorDyn input file header line 1', ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       IF ( ErrStat >= AbortErrLev ) THEN
          CALL CleanUp()
          RETURN
       END IF


    CALL ReadCom( UnIn, FileName, 'MoorDyn input file header line 2', ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       IF ( ErrStat >= AbortErrLev ) THEN
          CALL CleanUp()
          RETURN
       END IF


    ! Echo Input Files.
    CALL ReadVar ( UnIn, FileName, InitInp%Echo, 'Echo', 'Echo Input', ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       IF ( ErrStat >= AbortErrLev ) THEN
          CALL CleanUp()
          RETURN
       END IF


    ! If we are Echoing the input then we should re-read the first three lines so that we can echo them
    ! using the NWTC_Library routines.  The echoing is done inside those routines via a global variable
    ! which we must store, set, and then replace on error or completion.

      IF ( InitInp%Echo ) THEN

         !print *, 'gonna try to open echo file'

         EchoFile = TRIM(p%RootName)//'.ech'                      ! open an echo file for writing

         !print *, 'name is ', EchoFile

         CALL GetNewUnit( UnEc )
         CALL OpenEcho ( UnEc, EchoFile, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            IF ( ErrStat >= AbortErrLev ) THEN
               CALL CleanUp()
               RETURN
            END IF

         REWIND(UnIn)      ! rewind to start of input file to re-read the first few lines




       CALL ReadCom( UnIn, FileName, 'MoorDyn input file header line 1', ErrStat2, ErrMsg2, UnEc )
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
          IF ( ErrStat >= AbortErrLev ) THEN
             CALL CleanUp()
             RETURN
          END IF

       CALL ReadCom( UnIn, FileName, 'MoorDyn input file header line 2', ErrStat2, ErrMsg2, UnEc )
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
          IF ( ErrStat >= AbortErrLev ) THEN
             CALL CleanUp()
             RETURN
          END IF


       ! Echo Input Files. Note this line is prevented from being echoed by the ReadVar routine.
       CALL ReadVar ( UnIn, FileName, InitInp%Echo, 'Echo', 'Echo the input file data', ErrStat2, ErrMsg2, UnEc )
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
          IF ( ErrStat >= AbortErrLev ) THEN
             CALL CleanUp()
             RETURN
          END IF

      !print *, 'at end of echo if statement'

    END IF


    !-------------------------------------------------------------------------------------------------
    !  Line Types Properties Section
    !-------------------------------------------------------------------------------------------------

    CALL ReadCom( UnIn, FileName, 'Line types header', ErrStat2, ErrMsg2, UnEc )
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
          IF ( ErrStat >= AbortErrLev ) THEN
             CALL CleanUp()
             RETURN
          END IF


    CALL ReadVar ( UnIn, FileName, p%NTypes, 'NTypes', 'Number of line types', ErrStat2, ErrMsg2, UnEc )
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
          IF ( ErrStat >= AbortErrLev ) THEN
             CALL CleanUp()
             RETURN
          END IF


    ! Table header
    DO I = 1, 2
       CALL ReadCom( UnIn, FileName, 'Line types table header', ErrStat2, ErrMsg2, UnEc )
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
          IF ( ErrStat >= AbortErrLev ) THEN
             CALL CleanUp()
             RETURN
          END IF
    END DO

    ! make sure NTypes isn't zero
    IF ( p%NTypes < 1 ) THEN
       CALL SetErrStat( ErrID_Fatal, 'NTypes parameter must be greater than zero.', ErrStat, ErrMsg, RoutineName )
       CALL CleanUp()
       RETURN
    END IF

    ! Allocate memory for LineTypeList array to hold line type properties
    ALLOCATE ( m%LineTypeList(p%NTypes), STAT = ErrStat2 )
    IF ( ErrStat2 /= 0 ) THEN
       CALL SetErrStat( ErrID_Fatal, 'Error allocating space for LineTypeList array.', ErrStat, ErrMsg, RoutineName )
       CALL CleanUp()
       RETURN
    END IF

    ! read each line
    DO I = 1,p%NTypes
          ! read the table entries   Name      Diam    MassDenInAir    EA        cIntDamp     Can     Cat    Cdn     Cdt     in the MoorDyn input file
       READ(UnIn,'(A)',IOSTAT=ErrStat2) Line      !read into a line

       IF (ErrStat2 == 0) THEN
          READ(Line,*,IOSTAT=ErrStat2) m%LineTypeList(I)%name, m%LineTypeList(I)%d,  &
             m%LineTypeList(I)%w, m%LineTypeList(I)%EA, m%LineTypeList(I)%BA, &
             m%LineTypeList(I)%Can, m%LineTypeList(I)%Cat, m%LineTypeList(I)%Cdn, m%LineTypeList(I)%Cdt
       END IF

       m%LineTypeList(I)%IdNum = I  ! specify IdNum of line type for error checking


       IF ( ErrStat2 /= ErrID_None ) THEN
          CALL SetErrStat( ErrID_Fatal, 'Failed to read line type properties for line '//trim(Num2LStr(I)), ErrStat, ErrMsg, RoutineName )
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

    CALL ReadCom( UnIn, FileName, 'Connections header', ErrStat2, ErrMsg2, UnEc )
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
          IF ( ErrStat >= AbortErrLev ) THEN
             CALL CleanUp()
             RETURN
          END IF


    CALL ReadVar ( UnIn, FileName, p%NConnects, 'NConnects', 'Number of Connects', ErrStat2, ErrMsg2, UnEc )
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
          IF ( ErrStat >= AbortErrLev ) THEN
             CALL CleanUp()
             RETURN
          END IF


    ! Table header
    DO I = 1, 2
       CALL ReadCom( UnIn, FileName, 'Connects header', ErrStat2, ErrMsg2, UnEc )
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
          IF ( ErrStat >= AbortErrLev ) THEN
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
    ALLOCATE ( m%ConnectList(p%NConnects), STAT = ErrStat2 )
    IF ( ErrStat2 /= 0 ) THEN
       CALL SetErrStat( ErrID_Fatal, 'Error allocating space for ConnectList array.', ErrStat, ErrMsg, RoutineName )
       CALL CleanUp()
       RETURN
    END IF
    

    ! read each line
    DO I = 1,p%NConnects
          ! read the table entries   Node      Type      X        Y         Z        M        V        FX       FY      FZ  Cda Ca
       READ(UnIn,'(A)',IOSTAT=ErrStat2) Line      !read into a line

       IF (ErrStat2 == 0) THEN
          READ(Line,*,IOSTAT=ErrStat2) m%ConnectList(I)%IdNum, m%ConnectList(I)%type, m%ConnectList(I)%conX, &
               m%ConnectList(I)%conY, m%ConnectList(I)%conZ, m%ConnectList(I)%conM, &
               m%ConnectList(I)%conV, m%ConnectList(I)%conFX, m%ConnectList(I)%conFY, &
                m%ConnectList(I)%conFZ, m%ConnectList(I)%conCdA, m%ConnectList(I)%conCa
       END IF

       IF ( ErrStat2 /= 0 ) THEN
          CALL WrScr('   Unable to parse Connection '//trim(Num2LStr(I))//' row in input file.')  ! Specific screen output because errors likely
          CALL WrScr('   Ensure row has all 12 columns, including CdA and Ca.')           ! to be caused by non-updated input file formats.
             CALL SetErrStat( ErrID_Fatal, 'Failed to read connects.' , ErrStat, ErrMsg, RoutineName ) ! would be nice to specify which line <<<<<<<<<
          CALL CleanUp()
          RETURN
       END IF
       
       ! check for sequential IdNums
       IF ( m%ConnectList(I)%IdNum .NE. I ) THEN
         CALL SetErrStat( ErrID_Fatal, 'Node numbers must be sequential starting from 1.', ErrStat, ErrMsg, RoutineName )
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

    CALL ReadCom( UnIn, FileName, 'Lines header', ErrStat2, ErrMsg2, UnEc )
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
          IF ( ErrStat >= AbortErrLev ) THEN
             CALL CleanUp()
             RETURN
          END IF


    CALL ReadVar ( UnIn, FileName, p%NLines, 'NLines', 'Number of Lines', ErrStat2, ErrMsg2, UnEc )
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
          IF ( ErrStat >= AbortErrLev ) THEN
             CALL CleanUp()
             RETURN
          END IF


    ! Table header
    DO I = 1, 2
       CALL ReadCom( UnIn, FileName, 'Lines header', ErrStat2, ErrMsg2, UnEc )
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
          IF ( ErrStat >= AbortErrLev ) THEN
             CALL CleanUp()
             RETURN
          END IF
    END DO

    ! make sure NLines is at least one
    IF ( p%NLines < 1 ) THEN
       CALL SetErrStat( ErrID_Fatal, 'NLines parameter must be at least 1.', ErrStat, ErrMsg, RoutineName )
       CALL CleanUp()
       RETURN
    END IF

     ! allocate LineList
    ALLOCATE ( m%LineList(p%NLines), STAT = ErrStat2 )
    IF ( ErrStat2 /= 0 ) THEN
       CALL SetErrStat( ErrID_Fatal, 'Error allocating space for LineList array.', ErrStat, ErrMsg, RoutineName )
       CALL CleanUp()
       RETURN
    END IF

    ! read each line
    DO I = 1,p%NLines
          ! read the table entries   Line     LineType  UnstrLen  NumSegs   NodeAnch  NodeFair  Flags/Outputs
       READ(UnIn,'(A)',IOSTAT=ErrStat2) Line      !read into a line


       IF (ErrStat2 == 0) THEN
          READ(Line,*,IOSTAT=ErrStat2) m%LineList(I)%IdNum, m%LineList(I)%type, m%LineList(I)%UnstrLen, &
            m%LineList(I)%N, m%LineList(I)%AnchConnect, m%LineList(I)%FairConnect, LineOutString, m%LineList(I)%CtrlChan
       END IF

       IF ( ErrStat2 /= 0 ) THEN
          CALL SetErrStat( ErrID_Fatal, 'Failed to read line data for Line '//trim(Num2LStr(I)), ErrStat, ErrMsg, RoutineName )
          CALL CleanUp()
          RETURN
       END IF
       
       
       ! check for sequential IdNums
       IF ( m%LineList(I)%IdNum .NE. I ) THEN
         CALL SetErrStat( ErrID_Fatal, 'Line numbers must be sequential starting from 1.', ErrStat, ErrMsg, RoutineName )
         CALL CleanUp()
         RETURN
       END IF

       ! identify index of line type
       DO J = 1,p%NTypes
         IF (trim(m%LineList(I)%type) == trim(m%LineTypeList(J)%name)) THEN
           m%LineList(I)%PropsIdNum = J
           EXIT
           IF (J == p%NTypes) THEN   ! call an error if there is no match
               CALL SetErrStat( ErrID_Severe, 'Unable to find matching line type name for Line '//trim(Num2LStr(I)), ErrStat, ErrMsg, RoutineName )
           END IF
         END IF
       END DO
       
       ! process output flag characters (LineOutString) and set line output flag array (OutFlagList)
       m%LineList(I)%OutFlagList = 0  ! first set array all to zero
       IF ( scan( LineOutString, 'p') > 0 )  m%LineList(I)%OutFlagList(2) = 1 
       IF ( scan( LineOutString, 'v') > 0 )  m%LineList(I)%OutFlagList(3) = 1
       IF ( scan( LineOutString, 'U') > 0 )  m%LineList(I)%OutFlagList(4) = 1
       IF ( scan( LineOutString, 'D') > 0 )  m%LineList(I)%OutFlagList(5) = 1
       IF ( scan( LineOutString, 't') > 0 )  m%LineList(I)%OutFlagList(6) = 1
       IF ( scan( LineOutString, 'c') > 0 )  m%LineList(I)%OutFlagList(7) = 1
       IF ( scan( LineOutString, 's') > 0 )  m%LineList(I)%OutFlagList(8) = 1
       IF ( scan( LineOutString, 'd') > 0 )  m%LineList(I)%OutFlagList(9) = 1
       IF ( scan( LineOutString, 'l') > 0 )  m%LineList(I)%OutFlagList(10)= 1
       IF (SUM(m%LineList(I)%OutFlagList) > 0)   m%LineList(I)%OutFlagList(1) = 1  ! this first entry signals whether to create any output file at all
       ! the above letter-index combinations define which OutFlagList entry corresponds to which output type
       
   
       ! check errors
       IF ( ErrStat /= ErrID_None ) THEN
          ErrMsg  = ' Failed to read line data for Line '//trim(Num2LStr(I))
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

    CALL ReadCom( UnIn, FileName, 'Options header', ErrStat2, ErrMsg2, UnEc )
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
          IF ( ErrStat >= AbortErrLev ) THEN
             CALL CleanUp()
             RETURN
          END IF

     ! loop through any remaining input lines, and use them to set options (overwriting default values in many cases).
     ! doing this manually since I'm not sure that there is a built in subroutine for reading any input value on any line number.
     DO

       READ(UnIn,'(A)',IOSTAT=ErrStat2) Line      !read into a line

       IF (ErrStat2 == 0) THEN
         IF (( Line(1:3) == '---' ) .OR. ( Line(1:3) == 'END' ) .OR. ( Line(1:3) == 'end' ))  EXIT  ! check if it's the end line

         READ(Line,*,IOSTAT=ErrStat2) OptValue, OptString  ! look at first two entries, ignore remaining words in line, which should be comments
       END IF

       IF ( ErrStat2 /= 0 ) THEN
          CALL SetErrStat( ErrID_Fatal, 'Failed to read options.', ErrStat, ErrMsg, RoutineName ) ! would be nice to specify which line had the error
          CALL CleanUp()
          RETURN
       END IF
                     
       CALL Conv2UC(OptString)

       ! check all possible options types and see if OptString is one of them, in which case set the variable.
       if ( OptString == 'DTM') THEN
         read (OptValue,*) p%dtM0   ! InitInp%DTmooring
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
         CALL SetErrStat( ErrID_Warn, 'unable to interpret input '//trim(OptString), ErrStat, ErrMsg, RoutineName ) 
       end if

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
    CALL AllocAry( InitInp%OutList, 1000, "MoorDyn Input File's Outlist", ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF

    ! OutList - List of user-requested output channels (-):
    CALL ReadOutputList ( UnIn, FileName, InitInp%OutList, p%NumOuts, 'OutList', "List of user-requested output channels", ErrStat2, ErrMsg2, UnEc )
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
          IF ( ErrStat >= AbortErrLev ) THEN
             CALL CleanUp()
             RETURN
          END IF

   !print *, 'NumOuts is ', p%NumOuts
   !print *, '  OutList is ', InitInp%OutList(1:p%NumOuts)


     !-------------------------------------------------------------------------------------------------
     ! This is the end of the input file
     !-------------------------------------------------------------------------------------------------

         CALL CleanUp()

   CONTAINS
     ! subroutine to set ErrState and close the files if an error occurs
     SUBROUTINE CleanUp()

        ! ErrStat = ErrID_Fatal  
        CLOSE( UnIn )
        IF (InitInp%Echo) CLOSE( UnEc )

     END SUBROUTINE

   END SUBROUTINE MDIO_ReadInput
   ! ====================================================================================================



  ! ====================================================================================================
  SUBROUTINE MDIO_ProcessOutList(OutList, p, m, y, InitOut, ErrStat, ErrMsg )

  ! This routine processes the output channels requested by OutList, checking for validity and setting
  ! the p%OutParam structures (of type MD_OutParmType) for each valid output.
  ! It assumes the value p%NumOuts has been set beforehand, and sets the values of p%OutParam.


    IMPLICIT                        NONE

    ! Passed variables
    CHARACTER(ChanLen),        INTENT(IN)     :: OutList(:)                  ! The list of user-requested outputs
    TYPE(MD_ParameterType),    INTENT(INOUT)  :: p                           ! The module parameters
    TYPE(MD_MiscVarType),      INTENT(INOUT)  :: m
    TYPE(MD_OutputType),       INTENT(INOUT)  :: y                           ! Initial system outputs (outputs are not calculated; only the output mesh is initialized)
    TYPE(MD_InitOutputType),   INTENT(INOUT)  :: InitOut                     ! Output for initialization routine
    INTEGER(IntKi),            INTENT(OUT)    :: ErrStat                     ! The error status code
    CHARACTER(*),              INTENT(OUT)    :: ErrMsg                      ! The error message, if an error occurred

    ! Local variables
    INTEGER                      :: I                                        ! Generic loop-counting index
!    INTEGER                      :: J                                        ! Generic loop-counting index
!    INTEGER                      :: INDX                                     ! Index for valid arrays

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
      i1 = scan( OutListTmp , '1234567890' )              ! first number in the string
      i2 = i1+verify( OutListTmp(i1+1:) , '1234567890' )  ! second letter start (assuming first character is a letter, i.e. i1>1)
      i3 = i2+scan( OutListTmp(i2+1:) , '1234567890' )    ! second number start
      i4 = i3+verify( OutListTmp(i3+1:) , '1234567890' )  ! third letter start
      !i5 = scan( OutListTmp(i1:) , '1234567890' )        ! find first letter after first number

      ! error check
      IF (i1 <= 1) THEN
         CALL DenoteInvalidOutput(p%OutParam(I)) ! flag as invalid
         CALL WrScr('Warning: invalid output specifier '//trim(OutListTmp)//'.  Starting character must be C or L.')
         CYCLE    ! <<<<<<<<<<< check correct usage
      END IF

        p%OutParam(I)%Name = OutListTmp  ! label channel with whatever name was inputted, for now


      ! figure out what type of output it is and process accordingly

      ! fairlead tension case (updated)
      IF (OutListTmp(1:i1-1) == 'FAIRTEN') THEN
        p%OutParam(I)%OType = 2                                     ! connection object type
        p%OutParam(I)%QType = Ten                                   ! tension quantity type
        p%OutParam(I)%Units = UnitList(Ten)                         ! set units according to QType
        READ (OutListTmp(i1:),*) oID                                ! this is the line number
        p%OutParam(I)%ObjID = m%LineList(oID)%FairConnect           ! get the connection ID of the fairlead
        p%OutParam(I)%NodeID = -1                                   ! not used.    m%LineList(oID)%N  ! specify node N (fairlead)

      ! achor tension case
      ELSE IF (OutListTmp(1:i1-1) == 'ANCHTEN') THEN
        p%OutParam(I)%OType = 2                                     ! connectoin object type
        p%OutParam(I)%QType = Ten                                   ! tension quantity type
        p%OutParam(I)%Units = UnitList(Ten)                         ! set units according to QType
        READ (OutListTmp(i1:),*) oID                                ! this is the line number
        p%OutParam(I)%ObjID = m%LineList(oID)%AnchConnect           ! get the connection ID of the fairlead
        p%OutParam(I)%NodeID = -1                                   ! not used.    m%LineList(oID)%0  ! specify node 0 (anchor)

      ! more general case
      ELSE

        ! what object type?
        ! Line case                                          ... L?N?xxxx
        IF (OutListTmp(1:i1-1) == 'L') THEN
          p%OutParam(I)%OType = 1                ! Line object type
          ! for now we'll just assume the next character(s) are "n" to represent node number:
          READ (OutListTmp(i3:i4-1),*) nID
          p%OutParam(I)%NodeID = nID
          qVal = OutListTmp(i4:)                 ! isolate quantity type string
        ! Connect case                                     ... C?xxx or Con?xxx
        ELSE IF (OutListTmp(1:1) == 'C') THEN
          p%OutParam(I)%OType = 2                ! Connect object type
          qVal = OutListTmp(i2:)                 ! isolate quantity type string

        ! should do fairlead option also!

        ! error
        ELSE
          CALL DenoteInvalidOutput(p%OutParam(I)) ! flag as invalid
          CALL WrScr('Warning: invalid output specifier '//trim(OutListTmp)//'.  Type must be L or C.')
          CYCLE
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
        ELSE IF (qVal == 'AY') THEN   ! fixed typo Nov 24
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
          CALL WrScr('Warning: invalid output specifier '//trim(OutListTmp)//'.  Quantity type not recognized.')
          CONTINUE
        END IF

      END IF

      ! also check whether each object index and node index (if applicable) is in range
      IF (p%OutParam(I)%OType==2) THEN
        IF (p%OutParam(I)%ObjID > p%NConnects) THEN
          CALL WrScr('Warning: output Connect index excedes number of Connects in requested output '//trim(OutListTmp)//'.')
          CALL DenoteInvalidOutput(p%OutParam(I)) ! flag as invalid
        END IF
      ELSE IF (p%OutParam(I)%OType==1) THEN
        IF (p%OutParam(I)%ObjID > p%NLines) THEN
          CALL WrScr('Warning: output Line index excedes number of Lines in requested output '//trim(OutListTmp)//'.')
          CALL DenoteInvalidOutput(p%OutParam(I)) ! flag as invalid
        END IF
        IF (p%OutParam(I)%NodeID > m%LineList(p%OutParam(I)%ObjID)%N) THEN
          CALL WrScr('Warning: output node index excedes number of nodes in requested output '//trim(OutListTmp)//'.')
          CALL DenoteInvalidOutput(p%OutParam(I)) ! flag as invalid
        ELSE IF (p%OutParam(I)%NodeID < 0) THEN
          CALL WrScr('Warning: output node index is less than zero in requested output '//trim(OutListTmp)//'.')
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
      ! <<<<<<<<<<< should do this for each line instead.
   !   ALLOCATE( LineWriteOutputs( 200),  STAT = ErrStat )
   !   IF ( ErrStat /= ErrID_None ) THEN
   !      ErrMsg  = ' Error allocating space for LineWriteOutputs array.'
   !      ErrStat = ErrID_Fatal
   !      RETURN
   !   END IF

      !Allocate WriteOuput
      ALLOCATE(        y%WriteOutput(  p%NumOuts), &
              InitOut%WriteOutputHdr(p%NumOuts), &
              InitOut%WriteOutputUnt(p%NumOuts),  STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating space for y%WriteOutput array.'
         ErrStat = ErrID_Fatal
         RETURN
      END IF

      ! allocate output array in each Line
      DO I=1,p%NLines
         ALLOCATE(m%LineList(I)%LineWrOutput( 1 + 3*(m%LineList(I)%N + 1)*SUM(m%LineList(I)%OutFlagList(2:5)) + m%LineList(I)%N*SUM(m%LineList(I)%OutFlagList(6:10)) ), STAT = ErrStat)  
         IF ( ErrStat /= ErrID_None ) THEN
            ErrMsg  = ' Error allocating space for a LineWrOutput array'
            ErrStat = ErrID_Fatal
            RETURN
         END IF
      END DO  ! I

      !print *, "y%WriteOutput allocated to size ", size(y%WriteOutput)

      ! These variables are to help follow the framework template, but the data in them is simply a copy of data
      ! already available in the OutParam data structure
      !  ALLOCATE ( InitOut%WriteOutputHdr(p%NumOuts+p%OutAllint*p%OutAllDims), STAT = ErrStat )
      !  ALLOCATE ( InitOut%WriteOutputUnt(p%NumOuts+p%OutAllint*p%OutAllDims), STAT = ErrStat )

      DO I = 1,p%NumOuts
         InitOut%WriteOutputHdr(I) = p%OutParam(I)%Name
         InitOut%WriteOutputUnt(I) = p%OutParam(I)%Units
      END DO


   CONTAINS

      SUBROUTINE DenoteInvalidOutput( OutParm )
         TYPE(MD_OutParmType), INTENT (INOUT)  :: OutParm

         OutParm%OType = 0  ! flag as invalid
         OutParm%Name = 'Invalid'
         OutParm%Units = ' - '

      END SUBROUTINE DenoteInvalidOutput

   END SUBROUTINE MDIO_ProcessOutList
   !====================================================================================================





   !====================================================================================================
   SUBROUTINE MDIO_OpenOutput( OutRootName,  p, m, InitOut, ErrStat, ErrMsg )
   !----------------------------------------------------------------------------------------------------

      CHARACTER(*),                  INTENT( IN    ) :: OutRootName          ! Root name for the output file
      TYPE(MD_ParameterType),        INTENT( INOUT ) :: p
      TYPE(MD_MiscVarType),          INTENT( INOUT ) :: m
      TYPE(MD_InitOutPutType ),      INTENT( IN    ) :: InitOut              !
      INTEGER,                       INTENT(   OUT ) :: ErrStat              ! a non-zero value indicates an error occurred
      CHARACTER(*),                  INTENT(   OUT ) :: ErrMsg               ! Error message if ErrStat /= ErrID_None

      INTEGER                                        :: I                    ! Generic loop counter
      INTEGER                                        :: J                    ! Generic loop counter
      CHARACTER(1024)                                :: OutFileName          ! The name of the output file  including the full path.
!      INTEGER                                        :: L                           ! counter for index in LineWrOutput
      INTEGER                                        :: LineNumOuts                 ! number of entries in LineWrOutput for each line
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
         OutFileName = TRIM(p%RootName)//'.out'
         CALL GetNewUnit( p%MDUnOut )

         CALL OpenFOutFile ( p%MDUnOut, OutFileName, ErrStat, ErrMsg )
         IF ( ErrStat > ErrID_None ) THEN
            ErrMsg = ' Error opening MoorDyn-level output file: '//TRIM(ErrMsg)
            ErrStat = ErrID_Fatal
            RETURN
         END IF


         !Write the names of the output parameters:

         Frmt = '(A10,'//TRIM(Int2LStr(p%NumOuts))//'(A1,A10))'

         WRITE(p%MDUnOut,Frmt, IOSTAT=ErrStat2)  TRIM( 'Time' ), ( p%Delim, TRIM( p%OutParam(I)%Name), I=1,p%NumOuts )

         WRITE(p%MDUnOut,Frmt)  TRIM( '(s)' ), ( p%Delim, TRIM( p%OutParam(I)%Units ), I=1,p%NumOuts )

 !     ELSE  ! if no outputs requested

 !        call wrscr('note, MDIO_OpenOutput thinks that no outputs have been requested.')

      END IF

      !--------------------------------------------------------------------------
      !                    now do the same for line output files
      !--------------------------------------------------------------------------

      !! allocate UnLineOuts
      !ALLOCATE(UnLineOuts(p%NLines))  ! should add error checking

      DO I = 1,p%NLines

         
         IF (m%LineList(I)%OutFlagList(1) == 1) THEN   ! only proceed if the line is flagged to output a file
           
            ! Open the file for output
            OutFileName = TRIM(p%RootName)//'.Line'//TRIM(Int2LStr(I))//'.out'
            CALL GetNewUnit( m%LineList(I)%LineUnOut )

            CALL OpenFOutFile ( m%LineList(I)%LineUnOut, OutFileName, ErrStat, ErrMsg )
            IF ( ErrStat > ErrID_None ) THEN
               ErrMsg = ' Error opening Line output file '//TRIM(ErrMsg)
               ErrStat = ErrID_Fatal
               RETURN
            END IF

                        
            ! calculate number of output entries (including time) to write for this line
            LineNumOuts = 1 + 3*(m%LineList(I)%N + 1)*SUM(m%LineList(I)%OutFlagList(2:5)) + m%LineList(I)%N*SUM(m%LineList(I)%OutFlagList(6:10))

            Frmt = '(A10,'//TRIM(Int2LStr(LineNumOuts))//'(A1,A10))'   ! should evenutally use user specified format?
            !Frmt = '(A10,'//TRIM(Int2LStr(3+3*m%LineList(I)%N))//'(A1,A10))'
            
            ! Write the names of the output parameters:  (these use "implied DO" loops)

            WRITE(m%LineList(I)%LineUnOut,'(A10)', advance='no', IOSTAT=ErrStat2)  TRIM( 'Time' )
            IF (m%LineList(I)%OutFlagList(2) == 1) THEN
               WRITE(m%LineList(I)%LineUnOut,'('//TRIM(Int2LStr((3+3*m%LineList(I)%N)))//'(A1,A10))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, 'Node'//TRIM(Int2Lstr(J))//'px', p%Delim, 'Node'//TRIM(Int2Lstr(J))//'py', p%Delim, 'Node'//TRIM(Int2Lstr(J))//'pz', J=0,(m%LineList(I)%N) )
            END IF
            IF (m%LineList(I)%OutFlagList(3) == 1) THEN
               WRITE(m%LineList(I)%LineUnOut,'('//TRIM(Int2LStr((3+3*m%LineList(I)%N)))//'(A1,A10))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, 'Node'//TRIM(Int2Lstr(J))//'vx', p%Delim, 'Node'//TRIM(Int2Lstr(J))//'vy', p%Delim, 'Node'//TRIM(Int2Lstr(J))//'vz', J=0,(m%LineList(I)%N) )
            END IF
            IF (m%LineList(I)%OutFlagList(4) == 1) THEN
               WRITE(m%LineList(I)%LineUnOut,'('//TRIM(Int2LStr((3+3*m%LineList(I)%N)))//'(A1,A10))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, 'Node'//TRIM(Int2Lstr(J))//'Ux', p%Delim, 'Node'//TRIM(Int2Lstr(J))//'Uy', p%Delim, 'Node'//TRIM(Int2Lstr(J))//'Uz', J=0,(m%LineList(I)%N) )
            END IF
            IF (m%LineList(I)%OutFlagList(5) == 1) THEN
               WRITE(m%LineList(I)%LineUnOut,'('//TRIM(Int2LStr((3+3*m%LineList(I)%N)))//'(A1,A10))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, 'Node'//TRIM(Int2Lstr(J))//'Dx', p%Delim, 'Node'//TRIM(Int2Lstr(J))//'Dy', p%Delim, 'Node'//TRIM(Int2Lstr(J))//'Dz', J=0,(m%LineList(I)%N) )
            END IF
            IF (m%LineList(I)%OutFlagList(6) == 1) THEN
               WRITE(m%LineList(I)%LineUnOut,'('//TRIM(Int2LStr((m%LineList(I)%N)))//'(A1,A10))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, 'Seg'//TRIM(Int2Lstr(J))//'Ten', J=1,(m%LineList(I)%N) )
            END IF
            IF (m%LineList(I)%OutFlagList(7) == 1) THEN
               WRITE(m%LineList(I)%LineUnOut,'('//TRIM(Int2LStr((m%LineList(I)%N)))//'(A1,A10))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, 'Seg'//TRIM(Int2Lstr(J))//'Dmp', J=1,(m%LineList(I)%N) )
            END IF
            IF (m%LineList(I)%OutFlagList(8) == 1) THEN
               WRITE(m%LineList(I)%LineUnOut,'('//TRIM(Int2LStr((m%LineList(I)%N)))//'(A1,A10))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, 'Seg'//TRIM(Int2Lstr(J))//'Str', J=1,(m%LineList(I)%N) )
            END IF
            IF (m%LineList(I)%OutFlagList(9) == 1) THEN
               WRITE(m%LineList(I)%LineUnOut,'('//TRIM(Int2LStr((m%LineList(I)%N)))//'(A1,A10))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, 'Seg'//TRIM(Int2Lstr(J))//'SRt', J=1,(m%LineList(I)%N) )
            END IF
            IF (m%LineList(I)%OutFlagList(10)== 1) THEN
               WRITE(m%LineList(I)%LineUnOut,'('//TRIM(Int2LStr((m%LineList(I)%N)))//'(A1,A10))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, 'Seg'//TRIM(Int2Lstr(J))//'Lst', J=1,(m%LineList(I)%N) )
            END IF
            
            WRITE(m%LineList(I)%LineUnOut,'(A1)', IOSTAT=ErrStat2) ' '  ! make line break at the end
            
            ! Now write the units line

            WRITE(m%LineList(I)%LineUnOut,'(A10)', advance='no', IOSTAT=ErrStat2)  TRIM( '(s)' )
            IF (m%LineList(I)%OutFlagList(2) == 1) THEN
               WRITE(m%LineList(I)%LineUnOut,'('//TRIM(Int2LStr((3+3*m%LineList(I)%N)))//'(A1,A10))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, '(m)', p%Delim, '(m)', p%Delim, '(m)', J=0,(m%LineList(I)%N) )
            END IF
            IF (m%LineList(I)%OutFlagList(3) == 1) THEN
               WRITE(m%LineList(I)%LineUnOut,'('//TRIM(Int2LStr((3+3*m%LineList(I)%N)))//'(A1,A10))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, '(m/s)', p%Delim, '(m/s)', p%Delim, '(m/s)', J=0,(m%LineList(I)%N) )
            END IF
            IF (m%LineList(I)%OutFlagList(4) == 1) THEN
               WRITE(m%LineList(I)%LineUnOut,'('//TRIM(Int2LStr((3+3*m%LineList(I)%N)))//'(A1,A10))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, '(m/s)', p%Delim, '(m/s)', p%Delim, '(m/s)', J=0,(m%LineList(I)%N) )
            END IF
            IF (m%LineList(I)%OutFlagList(5) == 1) THEN
               WRITE(m%LineList(I)%LineUnOut,'('//TRIM(Int2LStr((3+3*m%LineList(I)%N)))//'(A1,A10))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, '(N)', p%Delim, '(N)', p%Delim, '(N)', J=0,(m%LineList(I)%N) )
            END IF
            IF (m%LineList(I)%OutFlagList(6) == 1) THEN
               WRITE(m%LineList(I)%LineUnOut,'('//TRIM(Int2LStr((m%LineList(I)%N)))//'(A1,A10))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, '(N)', J=1,(m%LineList(I)%N) )
            END IF
            IF (m%LineList(I)%OutFlagList(7) == 1) THEN
               WRITE(m%LineList(I)%LineUnOut,'('//TRIM(Int2LStr((m%LineList(I)%N)))//'(A1,A10))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, '(N)', J=1,(m%LineList(I)%N) )
            END IF
            IF (m%LineList(I)%OutFlagList(8) == 1) THEN
               WRITE(m%LineList(I)%LineUnOut,'('//TRIM(Int2LStr((m%LineList(I)%N)))//'(A1,A10))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, '(-)', J=1,(m%LineList(I)%N) )
            END IF
            IF (m%LineList(I)%OutFlagList(9) == 1) THEN
               WRITE(m%LineList(I)%LineUnOut,'('//TRIM(Int2LStr((m%LineList(I)%N)))//'(A1,A10))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, '(1/s)', J=1,(m%LineList(I)%N) )
            END IF
            IF (m%LineList(I)%OutFlagList(10)== 1) THEN
               WRITE(m%LineList(I)%LineUnOut,'('//TRIM(Int2LStr((m%LineList(I)%N)))//'(A1,A10))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, '(m)', J=1,(m%LineList(I)%N) )
            END IF
            
            WRITE(m%LineList(I)%LineUnOut,'(A1)', IOSTAT=ErrStat2) ' '  ! make line break at the end
            
         END IF  ! if line is flagged for output file
         
      END DO ! I - line number


      ! need to fix error handling in this sub

   END SUBROUTINE MDIO_OpenOutput
   !====================================================================================================


   !====================================================================================================
   SUBROUTINE MDIO_CloseOutput ( p, m, ErrStat, ErrMsg )
      ! This function cleans up after running the MoorDyn output module.
      ! It closes the output files and releases memory.

      TYPE(MD_ParameterType),       INTENT( INOUT )  :: p                    ! data for this instance of the floating platform module
      TYPE(MD_MiscVarType),         INTENT( INOUT )  :: m                    ! data for this instance of the floating platform module
      INTEGER,                      INTENT(   OUT )  :: ErrStat              ! a non-zero value indicates an error occurred
      CHARACTER(*),                 INTENT(   OUT )  :: ErrMsg               ! Error message if ErrStat /= ErrID_None

      INTEGER(IntKi)       :: I  ! generic counter


      ErrStat = 0
      ErrMsg  = ""


      ! close main MoorDyn output file
      CLOSE( p%MDUnOut, IOSTAT = ErrStat )
         IF ( ErrStat /= 0 ) THEN
            ErrMsg = 'Error closing output file'
         END IF

      ! close individual line output files
      DO I=1,p%NLines
         CLOSE( m%LineList(I)%LineUnOut, IOSTAT = ErrStat )
            IF ( ErrStat /= 0 ) THEN
               ErrMsg = 'Error closing line output file'
            END IF
      END DO

      ! deallocate output arrays
      IF (ALLOCATED(m%MDWrOutput)) THEN
         DEALLOCATE(m%MDWrOutput)
      ENDIF
      DO I=1,p%NLines
         IF (ALLOCATED(m%LineList(I)%LineWrOutput)) THEN
            DEALLOCATE(m%LineList(I)%LineWrOutput)       ! this may be unnecessary and handled by Line destructor
         ENDIF
      END DO

   END SUBROUTINE MDIO_CloseOutput
   !====================================================================================================


   !====================================================================================================
   SUBROUTINE MDIO_WriteOutputs( Time, p, m, y, ErrStat, ErrMsg )
      ! This subroutine gathers the output data defined by the OutParams list and
      ! writes it to the output file opened in MDIO_OutInit()

      REAL(DbKi),                   INTENT( IN    ) :: Time                 ! Time for this output
      TYPE(MD_ParameterType),       INTENT( IN    ) :: p                    ! MoorDyn module's parameter data
      TYPE(MD_OutputType),          INTENT( INOUT ) :: y                    ! INTENT( OUT) : Initial system outputs (outputs are not calculated; only the output mesh is initialized)
      TYPE(MD_MiscVarType),         INTENT( INOUT ) :: m                    ! MoorDyn module's m data
      INTEGER,                      INTENT(   OUT ) :: ErrStat              ! returns a non-zero value when an error occurs
      CHARACTER(*),                 INTENT(   OUT ) :: ErrMsg               ! Error message if ErrStat /= ErrID_None

      INTEGER                                :: I                           ! Generic loop counter
      INTEGER                                :: J                           ! Generic loop counter
      INTEGER                                :: K                           ! Generic loop counter
      INTEGER                                :: L                           ! counter for index in LineWrOutput
      INTEGER                                :: LineNumOuts                 ! number of entries in LineWrOutput for each line
      CHARACTER(200)                         :: Frmt                        ! a string to hold a format statement


      IF ( .NOT. ALLOCATED( p%OutParam ) .OR. p%MDUnOut < 0 )  THEN
         ErrStat = ErrID_Fatal
         ErrMsg  = ' To write outputs for MoorDyn there must be a valid file ID and OutParam must be allocated.'
         RETURN
      ELSE
         ErrStat = ErrID_None
         ErrMsg  = ''
      END IF

      ! Return if there are no outputs
      if ( p%NumOuts < 1_IntKi ) return


      ! gather the required output quantities (INCOMPLETE!)
      DO I = 1,p%NumOuts

         IF (p%OutParam(I)%OType == 2) THEN  ! if dealing with a Connect output
            SELECT CASE (p%OutParam(I)%QType)
               CASE (PosX)
                  y%WriteOutput(I) = m%ConnectList(p%OutParam(I)%ObjID)%r(1)  ! x position
               CASE (PosY)
                  y%WriteOutput(I) = m%ConnectList(p%OutParam(I)%ObjID)%r(2) ! y position
               CASE (PosZ)
                  y%WriteOutput(I) = m%ConnectList(p%OutParam(I)%ObjID)%r(3) ! z position
               CASE (VelX)
                  y%WriteOutput(I) = m%ConnectList(p%OutParam(I)%ObjID)%rd(1) ! x velocity
               CASE (VelY)
                  y%WriteOutput(I) = m%ConnectList(p%OutParam(I)%ObjID)%rd(2) ! y velocity
               CASE (VelZ)
                  y%WriteOutput(I) = m%ConnectList(p%OutParam(I)%ObjID)%rd(3) ! z velocity
               CASE (Ten)
                  y%WriteOutput(I) = TwoNorm(m%ConnectList(p%OutParam(I)%ObjID)%Ftot)  ! total force magnitude on a connect (used eg. for fairlead and anchor tensions)
               CASE (FX)
                  y%WriteOutput(I) = m%ConnectList(p%OutParam(I)%ObjID)%Ftot(1)  ! total force in x - added Nov 24
               CASE (FY)
                  y%WriteOutput(I) = m%ConnectList(p%OutParam(I)%ObjID)%Ftot(2)  ! total force in y
               CASE (FZ)
                  y%WriteOutput(I) = m%ConnectList(p%OutParam(I)%ObjID)%Ftot(3)  ! total force in z
               CASE DEFAULT
                  y%WriteOutput(I) = 0.0_DbKi
                  ErrStat = ErrID_Warn
                  ErrMsg = ' Unsupported output quantity '//TRIM(Num2Lstr(p%OutParam(I)%QType))//' requested from Connection '//TRIM(Num2Lstr(p%OutParam(I)%ObjID))//'.'
            END SELECT

         ELSE IF (p%OutParam(I)%OType == 1) THEN  ! if dealing with a Line output

            SELECT CASE (p%OutParam(I)%QType)
               CASE (PosX)
                 y%WriteOutput(I) = m%LineList(p%OutParam(I)%ObjID)%r(1,p%OutParam(I)%NodeID)  ! x position
               CASE (PosY)
                 y%WriteOutput(I) = m%LineList(p%OutParam(I)%ObjID)%r(2,p%OutParam(I)%NodeID) ! y position
               CASE (PosZ)
                 y%WriteOutput(I) = m%LineList(p%OutParam(I)%ObjID)%r(3,p%OutParam(I)%NodeID) ! z position
               CASE (VelX)
                 y%WriteOutput(I) = m%LineList(p%OutParam(I)%ObjID)%rd(1,p%OutParam(I)%NodeID) ! x velocity
               CASE (VelY)
                 y%WriteOutput(I) = m%LineList(p%OutParam(I)%ObjID)%rd(2,p%OutParam(I)%NodeID) ! y velocity
               CASE (VelZ)
                 y%WriteOutput(I) = m%LineList(p%OutParam(I)%ObjID)%rd(3,p%OutParam(I)%NodeID) ! z velocity
               CASE (Ten)
                 y%WriteOutput(I) = TwoNorm(m%LineList(p%OutParam(I)%ObjID)%T(:,p%OutParam(I)%NodeID))  ! this is actually the segment tension ( 1 < NodeID < N )  Should deal with properly!
               CASE DEFAULT
                 y%WriteOutput(I) = 0.0_DbKi
                 ErrStat = ErrID_Warn
                 ErrMsg = ' Unsupported output quantity '//TRIM(Num2Lstr(p%OutParam(I)%QType))//' requested from Line '//TRIM(Num2Lstr(p%OutParam(I)%ObjID))//'.'
            END SELECT

         ELSE  ! it must be an invalid output, so write zero
            y%WriteOutput(I) = 0.0_DbKi

         END IF

      END DO ! I, loop through OutParam


      ! Write the output parameters to the file

      Frmt = '(F10.4,'//TRIM(Int2LStr(p%NumOuts))//'(A1,e12.6))'   ! should evenutally use user specified format?

      WRITE(p%MDUnOut,Frmt)  Time, ( p%Delim, y%WriteOutput(I), I=1,p%NumOuts )





      !------------------------------------------------------------------------
      ! now do the outputs for each line!  
      
      DO I=1,p%NLines
        
        IF (m%LineList(I)%OutFlagList(1) == 1) THEN    ! only proceed if the line is flagged to output a file
           
           ! calculate number of output entries to write for this line
           LineNumOuts = 3*(m%LineList(I)%N + 1)*SUM(m%LineList(I)%OutFlagList(2:5)) + m%LineList(I)%N*SUM(m%LineList(I)%OutFlagList(6:10))
           
           
           Frmt = '(F10.4,'//TRIM(Int2LStr(LineNumOuts))//'(A1,e12.6))'   ! should evenutally use user specified format?

           L = 1 ! start of index of line output file at first entry
           
           ! Time
      !     m%LineList(I)%LineWrOutput(L) = Time
      !     L = L+1
           
           ! Node positions
           IF (m%LineList(I)%OutFlagList(2) == 1) THEN
              DO J = 0,m%LineList(I)%N  ! note index starts at zero because these are nodes
                DO K = 1,3
                  m%LineList(I)%LineWrOutput(L) = m%LineList(I)%r(K,J)
                  L = L+1
                END DO
              END DO
           END IF         
           
           ! Node velocities
           IF (m%LineList(I)%OutFlagList(3) == 1) THEN
              DO J = 0,m%LineList(I)%N  ! note index starts at zero because these are nodes
                DO K = 1,3
                  m%LineList(I)%LineWrOutput(L) = m%LineList(I)%rd(K,J)
                  L = L+1
                END DO
              END DO
           END IF
           
           
           ! Node wave velocities (not implemented yet)
           IF (m%LineList(I)%OutFlagList(4) == 1) THEN
              DO J = 0,m%LineList(I)%N  ! note index starts at zero because these are nodes
                DO K = 1,3
                  m%LineList(I)%LineWrOutput(L) = 0.0
                  L = L+1
                END DO
              END DO
           END IF
           
           
           ! Node total hydrodynamic forces (except added mass - just drag for now)
           IF (m%LineList(I)%OutFlagList(5) == 1) THEN
              DO J = 0,m%LineList(I)%N  ! note index starts at zero because these are nodes
                DO K = 1,3
                  m%LineList(I)%LineWrOutput(L) = m%LineList(I)%Dp(K,J) + m%LineList(I)%Dq(K,J)
                  L = L+1
                END DO
              END DO
           END IF
           
           
           ! Segment tension force (excludes damping term, just EA)
           IF (m%LineList(I)%OutFlagList(6) == 1) THEN
              DO J = 1,m%LineList(I)%N  
                m%LineList(I)%LineWrOutput(L) = TwoNorm(m%LineList(I)%T(:,J) )
                L = L+1
              END DO
           END IF
           
           ! Segment internal damping force
           IF (m%LineList(I)%OutFlagList(7) == 1) THEN
              DO J = 1,m%LineList(I)%N  
                 IF (( m%LineList(I)%Td(3,J)*m%LineList(I)%T(3,J) ) > 0)  THEN  ! if statement for handling sign (positive = tension)
                    m%LineList(I)%LineWrOutput(L) = TwoNorm(m%LineList(I)%Td(:,J) )
                 ELSE
                    m%LineList(I)%LineWrOutput(L) = -TwoNorm(m%LineList(I)%Td(:,J) )
                 END IF
                 L = L+1
              END DO
           END IF
           
           ! Segment strain
           IF (m%LineList(I)%OutFlagList(8) == 1) THEN
              DO J = 1,m%LineList(I)%N  
                m%LineList(I)%LineWrOutput(L) = m%LineList(I)%lstr(J)/m%LineList(I)%l(J) - 1.0 
                L = L+1
              END DO
           END IF
           
           ! Segment strain rate
           IF (m%LineList(I)%OutFlagList(9) == 1) THEN
              DO J = 1,m%LineList(I)%N  
                m%LineList(I)%LineWrOutput(L) = m%LineList(I)%lstrd(J)/m%LineList(I)%l(J)
                L = L+1
              END DO
           END IF
           
           ! Segment length
           IF (m%LineList(I)%OutFlagList(10) == 1) THEN
              DO J = 1,m%LineList(I)%N  
                m%LineList(I)%LineWrOutput(L) = m%LineList(I)%lstr(J)
                L = L+1
              END DO
           END IF
                    
           
           WRITE(m%LineList(I)%LineUnOut,Frmt) Time, ( p%Delim, m%LineList(I)%LineWrOutput(J), J=1,(LineNumOuts) )
           !WRITE(m%LineList(I)%LineUnOut,Frmt)  Time, ( p%Delim, m%LineList(I)%LineWrOutput(J), J=1,(3+3*m%LineList(I)%N) )

         END IF  ! if line output file flag is on
           
      END DO ! I

   END SUBROUTINE MDIO_WriteOutputs
   !====================================================================================================


END MODULE MoorDyn_IO
