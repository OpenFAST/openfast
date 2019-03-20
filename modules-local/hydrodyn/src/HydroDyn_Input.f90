!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2013-2015  National Renewable Energy Laboratory
!
!    This file is part of HydroDyn.
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
MODULE HydroDyn_Input

      ! This MODULE stores variables used for input.

   USE                              NWTC_Library
   USE                              HydroDyn_Types
   USE                              HydroDyn_Output
   USE                              Waves
   USE                              Morison
   USE                              WAMIT_Output
   USE                              WAMIT2_Output
   USE                              Waves2_Output
   USE                              Morison_Output
   IMPLICIT                         NONE

   PRIVATE :: CleanupEchoFile
   PRIVATE :: CheckMeshOutput

CONTAINS

!====================================================================================================
FUNCTION CheckMeshOutput( output, numMemberOut, MOutLst, numJointOut )
!     The routine
!----------------------------------------------------------------------------------------------------
!
   CHARACTER(10),             INTENT ( IN    )  :: output
   INTEGER,                   INTENT ( IN    )  :: numMemberOut
   TYPE(Morison_MOutput),     INTENT ( IN    )  :: MOutLst(:)
   INTEGER,                   INTENT ( IN    )  :: numJointOut
   !INTEGER,                   INTENT (   OUT )  :: ErrStat              ! returns a non-zero value when an error occurs
   !CHARACTER(*),              INTENT (   OUT )  :: ErrMsg               ! Error message if ErrStat /= ErrID_None

   LOGICAL                                      :: CheckMeshOutput

   INTEGER                                      :: ErrStat
   CHARACTER(10)                                :: outputTmp
   INTEGER                                      :: indx1, indx2
   CHARACTER(4)                                 :: testStr
   outputTmp         = TRIM(output)

   testStr = outputTmp(1:4)
   CALL Conv2UC( testStr )

      ! Reverse the sign (+/-) of the output channel if the user prefixed the
      !   channel name with a '-', '_', 'm', or 'M' character indicating "minus".

      IF      ( INDEX( '-_', outputTmp(1:1) ) > 0 ) THEN

            ! ex, '-TipDxc1' causes the sign of TipDxc1 to be switched.
         outputTmp                   = outputTmp(2:)
         testStr = outputTmp(1:4)
         CALL Conv2UC( testStr )

      ELSE IF ( INDEX( 'mM', outputTmp(1:1) ) > 0 ) THEN ! We'll assume this is a variable name for now, (if not, we will check later if OutListTmp(2:) is also a variable name)

         IF ( ( INDEX( 'mM', outputTmp(2:2) ) > 0 ) .OR. ( INDEX( 'jJ', outputTmp(2:2) ) > 0 ) )  THEN
            outputTmp                   = outputTmp(2:)

         END IF

      ELSE IF ( INDEX( 'jJ', outputTmp(1:1) ) == 0  .AND. ( testStr /= 'WAVE' )  ) THEN
         ! Invalid output label because the label does not start: -M,-m,-J,-j,_M,_m,_J,_j,MM,mM,Mm,mm,MJ,mJ,Mj,mj, j,J,m,M
         CheckMeshOutput = .FALSE.
         RETURN
      END IF

      IF (( INDEX( 'mM', outputTmp(1:1) ) > 0 ) .OR. ( INDEX( 'jJ', outputTmp(1:1) ) > 0 )) THEN
         ! Read the second character, it should be a number from 1 to 9
      
         READ( outputTmp(2:2), '(i1)', IOSTAT = ErrStat) indx1
         IF ( ErrStat /=0 ) THEN
            ! Not a numerical digit!!!
            CheckMeshOutput = .FALSE.
            RETURN
         END IF
      
            ! Examine members
         IF ( INDEX( 'mM', outputTmp(1:1) ) > 0 ) THEN 
            IF ( indx1 > numMemberOut ) THEN
               CheckMeshOutput = .FALSE.
               RETURN
            END IF
               ! Now make sure the next letter is n or N and then look for the second index
               IF ( INDEX( 'nN', outputTmp(3:3) ) == 0 ) THEN
                     ! Invalid member label
                  CheckMeshOutput = .FALSE.
                  RETURN
               END IF
               READ( outputTmp(4:4), '(i1)', IOSTAT = ErrStat) indx2
               IF ( indx2 > MOutLst(indx1)%NOutLoc ) THEN
                  CheckMeshOutput = .FALSE.
                  RETURN
               END IF
            
         
         END IF 
      
         IF ( INDEX( 'jJ', outputTmp(1:1) ) > 0 ) THEN 
            IF ( indx1 > numJointOut ) THEN
               CheckMeshOutput = .FALSE.
               RETURN
            END IF
         END IF 
   ELSE 
         ! This should be a wave elevation channel
      READ( outputTmp(5:5), '(i1)', IOSTAT = ErrStat) indx1
      IF ( ErrStat /=0 ) THEN
         ! Not a numerical digit!!!
         CheckMeshOutput = .FALSE.
         RETURN
      END IF   
   END IF

      CheckMeshOutput = .TRUE.

END FUNCTION CheckMeshOutput

!====================================================================================================
SUBROUTINE PrintBadChannelWarning(NUserOutputs, UserOutputs , foundMask, ErrStat, ErrMsg )
!     The routine prints out warning messages if the user has requested invalid output channel names
!     The errstat is set to ErrID_Warning if any element in foundMask is .FALSE.
!----------------------------------------------------------------------------------------------------  
   INTEGER,                       INTENT( IN    ) :: NUserOutputs         ! Number of user-specified output channels
   CHARACTER(10),                 INTENT( IN    ) :: UserOutputs (:)      ! An array holding the names of the requested output channels. 
   LOGICAL,                       INTENT( IN    ) :: foundMask (:)        ! A mask indicating whether a user requested channel belongs to a module's output channels.
   INTEGER,                       INTENT(   OUT ) :: ErrStat              ! returns a non-zero value when an error occurs  
   CHARACTER(*),                  INTENT(   OUT ) :: ErrMsg               ! Error message if ErrStat /= ErrID_None
   INTEGER                                        :: I
   
   ErrStat = ErrID_None
   ErrMsg  = ''
   
   DO I = 1, NUserOutputs
      IF (.NOT. foundMask(I)) THEN
         ErrMsg  = ' A requested output channel is invalid'         
         CALL ProgWarn( 'The requested output channel is invalid: ' // UserOutputs(I) )
         ErrStat = ErrID_Warn
      END IF
   END DO
   
   
   
END SUBROUTINE PrintBadChannelWarning


!====================================================================================================
SUBROUTINE CleanupEchoFile( EchoFlag, UnEcho)
!     The routine cleans up the module echo file and resets the NWTC_Library, reattaching it to
!     any existing echo information
!----------------------------------------------------------------------------------------------------
   LOGICAL,                       INTENT( IN    )   :: EchoFlag             ! local version of echo flag
   INTEGER,                       INTENT( IN    )   :: UnEcho               !  echo unit number


      ! Close this module's echo file

   IF ( EchoFlag ) THEN
    CLOSE(UnEcho)
   END IF



END SUBROUTINE CleanupEchoFile




!====================================================================================================
SUBROUTINE HydroDynInput_GetInput( InitInp, ErrStat, ErrMsg )
!     This public subroutine reads the input required for HydroDyn from the file whose name is an
!     input parameter.
!----------------------------------------------------------------------------------------------------


      ! Passed variables

   TYPE(HydroDyn_InitInputType),  INTENT( INOUT )   :: InitInp              ! the hydrodyn data
   INTEGER,                       INTENT(   OUT )   :: ErrStat              ! returns a non-zero value when an error occurs
   CHARACTER(*),                  INTENT(   OUT )   :: ErrMsg               ! Error message if ErrStat /= ErrID_None


      ! Local variables

   INTEGER                                          :: I                    ! generic integer for counting
!   INTEGER                                          :: J                    ! generic integer for counting
   CHARACTER(   2)                                  :: strI                 ! string version of the loop counter

   INTEGER                                          :: UnIn                 ! Unit number for the input file
!   LOGICAL                                          :: EchoStore            ! Stored version of NWTC_Library Echo variable
!   INTEGER                                          :: UnEchoStore          ! Stored unit name for another module's echo file
   INTEGER                                          :: UnEchoLocal          ! The local unit number for this module's echo file
   CHARACTER(1024)                                  :: EchoFile             ! Name of HydroDyn echo file
   CHARACTER(1024)                                  :: Line                 ! String to temporarially hold value of read line
!   CHARACTER(1024)                                  :: TmpPath              ! Temporary storage for relative path name
!   CHARACTER(1024)                                  :: TmpFmt               ! Temporary storage for format statement
   CHARACTER(1024)                                  :: FileName             ! Name of HydroDyn input file
   CHARACTER(  35)                                  :: Frmt                 ! Output format for logical parameters. (matches NWTC Subroutine Library format)
!   INTEGER                                          :: JointID              ! Temporary storage of JointID read from HydroDyn input file
!   INTEGER                                          :: PropSetID            ! Temporary storage of PropSetID read from HydroDyn input file
!   INTEGER                                          :: MemberID             ! Temporary storage of MemberID read from HydroDyn input file
   INTEGER, ALLOCATABLE                             :: tmpArray(:)          ! Temporary array storage of the joint output list


   INTEGER(IntKi)                                   :: ErrStat2
   CHARACTER(ErrMsgLen)                             :: ErrMsg2
   
   
      ! Initialize local data

   UnEchoLocal  = -1
   Frmt         = "( 2X, L11, 2X, A, T30, ' - ', A )"         
   ErrStat      = ErrID_None         
   ErrMsg       = ""   
   InitInp%Echo = .FALSE.  ! initialize for error handling (cleanup() routine)
   
   !-------------------------------------------------------------------------------------------------
   ! Open the file
   !-------------------------------------------------------------------------------------------------
   FileName = TRIM(InitInp%InputFile)

   CALL GetNewUnit( UnIn, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
   CALL OpenFInpFile( UnIn, FileName, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup() 
         RETURN
      END IF


   !CALL WrScr( 'Opening HydroDyn input file:  '//FileName )


   !-------------------------------------------------------------------------------------------------
   ! File header
   !-------------------------------------------------------------------------------------------------

   CALL ReadCom( UnIn, FileName, 'HydroDyn input file header line 1', ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup() 
         RETURN
      END IF
      
      
   CALL ReadCom( UnIn, FileName, 'HydroDyn input file header line 2', ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup() 
         RETURN
      END IF

     ! Echo Input Files.

   CALL ReadVar ( UnIn, FileName, InitInp%Echo, 'Echo', 'Echo Input', ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )      
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup() 
         RETURN
      END IF

      ! If we are Echoing the input then we should re-read the first three lines so that we can echo them
      ! using the NWTC_Library routines.  The echoing is done inside those routines via a global variable
      ! which we must store, set, and then replace on error or completion.

   IF ( InitInp%Echo ) THEN

      EchoFile = TRIM(InitInp%OutRootName)//'.HD.ech'
      CALL OpenEcho ( UnEchoLocal, TRIM(EchoFile), ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
         IF (ErrStat >= AbortErrLev) THEN
            CALL CleanUp()
            RETURN
         END IF
         
      REWIND(UnIn)

      CALL ReadCom( UnIn, FileName, 'HydroDyn input file header line 1', ErrStat2, ErrMsg2, UnEchoLocal )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' ) 

      CALL ReadCom( UnIn, FileName, 'HydroDyn input file header line 2', ErrStat2, ErrMsg2, UnEchoLocal )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )

         ! Echo Input Files. Note this line is prevented from being echoed by the ReadVar routine. (bjj: is that still true?)

      CALL ReadVar ( UnIn, FileName, InitInp%Echo, 'Echo', 'Echo the input file data', ErrStat2, ErrMsg2, UnEchoLocal )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )

   END IF
   !-------------------------------------------------------------------------------------------------
   ! Environmental conditions section
   !-------------------------------------------------------------------------------------------------

      ! Header

   CALL ReadCom( UnIn, FileName, 'Environmental conditions header', ErrStat2, ErrMsg2, UnEchoLocal )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF


      ! WtrDens - Water density.

   CALL ReadVar ( UnIn, FileName, InitInp%Waves%WtrDens, 'WtrDens', 'Water density', ErrStat2, ErrMsg2, UnEchoLocal )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF


      ! WtrDpth - Water depth

   CALL ReadVar ( UnIn, FileName, InitInp%Morison%WtrDpth, 'WtrDpth', 'Water depth', ErrStat2, ErrMsg2, UnEchoLocal )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF


      ! MSL2SWL

   CALL ReadVar ( UnIn, FileName, InitInp%Morison%MSL2SWL, 'MSL2SWL', 'MSL to SWL offset', ErrStat2, ErrMsg2, UnEchoLocal )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF


   !-------------------------------------------------------------------------------------------------
   ! Data section for waves
   !-------------------------------------------------------------------------------------------------

      ! Header

   CALL ReadCom( UnIn, FileName, 'Wave header', ErrStat2, ErrMsg2, UnEchoLocal )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF


      ! WaveMod - Wave kinematics model switch.

   CALL ReadVar ( UnIn, FileName, InitInp%Waves%WaveModChr, 'WaveMod', 'Wave kinematics model switch', ErrStat2, ErrMsg2, UnEchoLocal )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF

   CALL Conv2UC( InitInp%Waves%WaveModChr )    ! Convert Line to upper case.

   InitInp%Waves%WavePhase = 0.0
   InitInp%Waves%WaveNDAmp = .FALSE.


      ! WaveStMod - Model switch for stretching incident wave kinematics to instantaneous free surface.

   CALL ReadVar ( UnIn, FileName, InitInp%Waves%WaveStMod, 'WaveStMod', &
      'Model switch for stretching incident wave kinematics to instantaneous free surface', ErrStat2, ErrMsg2, UnEchoLocal )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF



      ! WaveTMax - Analysis time for incident wave calculations.


   CALL ReadVar ( UnIn, FileName, InitInp%Waves%WaveTMax, 'WaveTMax', &
                              'Analysis time for incident wave calculations', ErrStat2, ErrMsg2, UnEchoLocal )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF


      ! WaveDT - Time step for incident wave calculations


   CALL ReadVar ( UnIn, FileName, InitInp%Waves%WaveDT, 'WaveDT', &
                        'Time step for incident wave calculations', ErrStat2, ErrMsg2, UnEchoLocal )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF




      ! WaveHs - Significant wave height


   CALL ReadVar ( UnIn, FileName, InitInp%Waves%WaveHs, 'WaveHs', 'Significant wave height', ErrStat2, ErrMsg2, UnEchoLocal )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF



      ! WaveTp - Peak spectral period.

   CALL ReadVar ( UnIn, FileName, InitInp%Waves%WaveTp, 'WaveTp', 'Peak spectral period', ErrStat2, ErrMsg2, UnEchoLocal )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF


      ! WavePkShp - Peak shape parameter.

   CALL ReadVar ( UnIn, FileName, InitInp%Waves%WavePkShpChr, 'WavePkShp', 'Peak shape parameter', ErrStat2, ErrMsg2, UnEchoLocal )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF


      ! WvLowCOff - Low Cut-off frequency or lower frequency limit of the wave spectrum beyond which the wave spectrum is zeroed (rad/s).  

   CALL ReadVar ( UnIn, FileName, InitInp%Waves%WvLowCOff, 'WvLowCOff', 'Lower wave cut-off frequency', ErrStat2, ErrMsg2, UnEchoLocal )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF

     ! WvHiCOff - High Cut-off frequency or upper frequency limit of the wave spectrum beyond which the wave spectrum is zeroed (rad/s).  

   CALL ReadVar ( UnIn, FileName, InitInp%Waves%WvHiCOff, 'WvHiCOff', 'Upper wave cut-off frequency', ErrStat2, ErrMsg2, UnEchoLocal )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF
   
      ! WaveDir - Mean wave heading direction.

   CALL ReadVar ( UnIn, FileName, InitInp%Waves%WaveDir, 'WaveDir', 'Mean wave heading direction', ErrStat2, ErrMsg2, UnEchoLocal )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF


      ! WaveDirMod -  Directional spreading function {0: None, 1: COS2S}       (-) [Used only if WaveMod=2]

   CALL ReadVar ( UnIn, FileName, InitInp%Waves%WaveDirMod, 'WaveDirMod', 'Directional spreading function', ErrStat2, ErrMsg2, UnEchoLocal )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF


      ! WaveDirSpread -  Spreading coefficient [only used if WaveMod=2 and WaveDirMod=1]

   CALL ReadVar ( UnIn, FileName, InitInp%Waves%WaveDirSpread, 'WaveDirSpread', 'Wave direction spreading coefficient', ErrStat2, ErrMsg2, UnEchoLocal )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF


      ! WaveNDir -  The number of wave directions to calculate [must be odd; only used if WaveDirMod=1]

   CALL ReadVar (UnIn, FileName, InitInp%Waves%WaveNDir, 'WaveNDir', 'Number of wave directions to calculate', ErrStat2, ErrMsg2, UnEchoLocal )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF


      ! WaveDirRange - Full range of the wave directions from WaveDir - WaveDirRange/2 to WaveDir + WaveDirRange/2 (only used if WaveMod=2 and WaveDirMod=1)

   CALL ReadVar ( UnIn, FileName, InitInp%Waves%WaveDirRange, 'WaveDirRange', 'Maximum wave heading direction', ErrStat2, ErrMsg2, UnEchoLocal )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF

      ! Negative values should be treated as positive.
   InitInp%Waves%WaveDirRange =  ABS( InitInp%Waves%WaveDirRange )


      ! WaveSeed(1), !WaveSeed(2)

   DO I = 1,2

      WRITE(Line,'(I2)') I

      CALL ReadVar ( UnIn, FileName, InitInp%Waves%WaveSeed(I), 'WaveSeed('//TRIM(Line)//')', &
                                    'Random seed #'//TRIM(Line), ErrStat2, ErrMsg2, UnEchoLocal )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
         IF (ErrStat >= AbortErrLev) THEN
            CALL CleanUp()
            RETURN
         END IF

   END DO !I


      ! WaveNDAmp - Flag for normally distributed amplitudes.

   CALL ReadVar ( UnIn, FileName, InitInp%Waves%WaveNDAmp, 'WaveNDAmp', 'Normally distributed amplitudes', ErrStat2, ErrMsg2, UnEchoLocal )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF



      ! WvKinFile

   CALL ReadVar ( UnIn, FileName, InitInp%Waves%WvKinFile, 'WvKinFile', &
                                    'Root name of wave kinematics files', ErrStat2, ErrMsg2, UnEchoLocal )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF




      ! NWaveElev

   CALL ReadVar ( UnIn, FileName, InitInp%Waves%NWaveElev, 'NWaveElev', &
                                  'Number of points where the incident wave elevations can be output', ErrStat2, ErrMsg2, UnEchoLocal )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF


      ! This check is needed here instead of being located in HydroDynInput_ProcessInputData() because
      ! we need to allocate arrays.  If _GetInput() was skipped, then these array would already have
      ! been allocated and populated.

   IF ( InitInp%Waves%NWaveElev < 0 .OR. InitInp%Waves%NWaveElev > 9 ) THEN

      CALL SetErrStat( ErrID_Fatal, 'NWaveElev must be greater than or equal to zero and less than 10.', ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      CALL CleanUp()
      RETURN

   ELSE

         ! allocate space for the output location arrays:
      CALL AllocAry( InitInp%Waves%WaveElevxi, InitInp%Waves%NWaveElev, 'WaveElevxi' , ErrStat2, ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      CALL AllocAry( InitInp%Waves%WaveElevyi, InitInp%Waves%NWaveElev, 'WaveElevyi' , ErrStat2, ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )

      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF      
      
   END IF

      ! WaveElevxi

   CALL ReadAry ( UnIn, FileName, InitInp%Waves%WaveElevxi, InitInp%Waves%NWaveElev, 'WaveElevxi', &
                           'List of xi-coordinates for points where the incident wave elevations can be output', ErrStat2,  ErrMsg2, UnEchoLocal)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF


      ! WaveElevyi

   CALL ReadAry ( UnIn, FileName, InitInp%Waves%WaveElevyi, InitInp%Waves%NWaveElev, 'WaveElevyi', &
                           'List of yi-coordinates for points where the incident wave elevations can be output', ErrStat2,  ErrMsg2, UnEchoLocal )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF



   !-------------------------------------------------------------------------------------------------
   ! Data section for 2nd Order Waves 
   !-------------------------------------------------------------------------------------------------

      ! Header

   CALL ReadCom( UnIn, FileName, 'Waves 2nd order', ErrStat2, ErrMsg2, UnEchoLocal )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF


      ! WvDiffQTFF     - Second order waves -- difference forces

   CALL ReadVar ( UnIn, FileName, InitInp%Waves2%WvDiffQTFF, 'WvDiffQTFF', 'Full difference QTF second order kinematic forces flag', ErrStat2, ErrMsg2, UnEchoLocal )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF


      ! WvSumQTFF      - Second order waves -- sum forces

   CALL ReadVar ( UnIn, FileName, InitInp%Waves2%WvSumQTFF, 'WvSumQTFF', 'Full sum QTF  second order kinematic forces flag', ErrStat2, ErrMsg2, UnEchoLocal )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF



        ! WvLowCOffD   -- Minimum frequency used in the difference methods (rad/s)              [Only used if DiffQTF /= 0]

   CALL ReadVar ( UnIn, FileName, InitInp%Waves2%WvLowCOffD, 'WvLowCOffD', 'Minimum frequency used in second order difference forces', ErrStat2, ErrMsg2, UnEchoLocal )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF


        ! WvHiCOffD   -- Maximum frequency used in the difference methods  (rad/s)              [Only used if DiffQTF /= 0]

   CALL ReadVar ( UnIn, FileName, InitInp%Waves2%WvHiCOffD, 'WvHiCOffD', 'Maximum frequency used in second order difference forces', ErrStat2, ErrMsg2, UnEchoLocal )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF


        ! WvLowCOffS   -- Minimum frequency used in the        sum-QTF     (rad/s)              [Only used if  SumQTF /= 0]

   CALL ReadVar ( UnIn, FileName, InitInp%Waves2%WvLowCOffS, 'WvLowCOffS', 'Minimum frequency used in second order sum forces', ErrStat2, ErrMsg2, UnEchoLocal )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF


        ! WvHiCOffS   -- Maximum frequency used in the        sum-QTF      (rad/s)              [Only used if  SumQTF /= 0]

   CALL ReadVar ( UnIn, FileName, InitInp%Waves2%WvHiCOffS, 'WvHiCOffS', 'Maximum frequency used in second order sum forces', ErrStat2, ErrMsg2, UnEchoLocal )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF



   !-------------------------------------------------------------------------------------------------
   ! Data section for current
   !-------------------------------------------------------------------------------------------------

      ! Header

   CALL ReadCom( UnIn, FileName, 'Current header', ErrStat2, ErrMsg2, UnEchoLocal )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF


      ! CurrMod - Current profile model switch

   CALL ReadVar ( UnIn, FileName, InitInp%Current%CurrMod, 'CurrMod', 'Current profile model switch', ErrStat2, ErrMsg2, UnEchoLocal )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF



      ! CurrSSV0 - Sub-surface current velocity at still water level

   CALL ReadVar ( UnIn, FileName, InitInp%Current%CurrSSV0, 'CurrSSV0', 'Sub-surface current velocity at still water level', &
                         ErrStat2, ErrMsg2, UnEchoLocal )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF



      ! CurrSSDirChr - Sub-surface current heading direction

   CALL ReadVar ( UnIn, FileName, InitInp%Current%CurrSSDirChr, 'CurrSSDirChr', 'Sub-surface current heading direction', ErrStat2, ErrMsg2, UnEchoLocal )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF


   CALL Conv2UC( InitInp%Current%CurrSSDirChr )    ! Convert Line to upper case.


      ! CurrNSRef - Near-surface current reference depth.

   CALL ReadVar ( UnIn, FileName, InitInp%Current%CurrNSRef, 'CurrNSRef', 'Near-surface current reference depth', ErrStat2, ErrMsg2, UnEchoLocal )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF


      ! CurrNSV0 - Near-surface current velocity at still water level.

   CALL ReadVar ( UnIn, FileName, InitInp%Current%CurrNSV0, 'CurrNSV0', 'Near-surface current velocity at still water level', &
                           ErrStat2, ErrMsg2, UnEchoLocal )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF


      ! CurrNSDir - Near-surface current heading direction.

   CALL ReadVar ( UnIn, FileName, InitInp%Current%CurrNSDir, 'CurrNSDir', 'Near-surface current heading direction', ErrStat2, ErrMsg2, UnEchoLocal )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF


      ! CurrDIV - Depth-independent current velocity.

   CALL ReadVar ( UnIn, FileName, InitInp%Current%CurrDIV, 'CurrDIV', 'Depth-independent current velocity', ErrStat2, ErrMsg2, UnEchoLocal )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF


      ! CurrDIDir - Depth-independent current heading direction.

   CALL ReadVar ( UnIn, FileName, InitInp%Current%CurrDIDir, 'CurrDIDir', 'Depth-independent current heading direction', ErrStat2, ErrMsg2, UnEchoLocal )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF



   !-------------------------------------------------------------------------------------------------
   ! Data section for floating platform
   !-------------------------------------------------------------------------------------------------

      ! Header

   CALL ReadCom( UnIn, FileName, 'Floating platform header', ErrStat2, ErrMsg2, UnEchoLocal )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF



      ! PotMod - State indicating potential flow model used in the simulation. 0=none, 1=WAMIT, 2=FIT

   CALL ReadVar ( UnIn, FileName, InitInp%PotMod, 'PotMod', 'Potential flow model', ErrStat2, ErrMsg2, UnEchoLocal )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF


      ! PotFile - Root name of Potential flow data files (Could be WAMIT files or the FIT input file)

   CALL ReadVar ( UnIn, FileName, InitInp%PotFile, 'PotFile', 'Root name of Potential flow model files', ErrStat2, ErrMsg2, UnEchoLocal )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF


      ! WAMITULEN - WAMIT characteristic body length scale

   CALL ReadVar ( UnIn, FileName, InitInp%WAMIT%WAMITULEN, 'WAMITULEN', 'WAMIT characteristic body length scale', ErrStat2, ErrMsg2, UnEchoLocal )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF


      ! PtfmVol0 - Displaced volume of water when the platform is in its undisplaced position

   CALL ReadVar ( UnIn, FileName, InitInp%WAMIT%PtfmVol0, 'PtfmVol0', &
      'Displaced volume of water when the platform is in its undisplaced position', ErrStat2, ErrMsg2, UnEchoLocal )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF


      ! PtfmCOBxt  - The xt offset of the center of buoyancy (COB) from the WAMIT reference point

   CALL ReadVar ( UnIn, FileName, InitInp%WAMIT%PtfmCOBxt, 'PtfmCOBxt', &
      'xt offset of the center of buoyancy (COB) from the WAMIT reference point', ErrStat2, ErrMsg2, UnEchoLocal )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF


      ! PtfmCOByt - The yt offset of the center of buoyancy (COB) from the WAMIT reference point

   CALL ReadVar ( UnIn, FileName, InitInp%WAMIT%PtfmCOByt, 'PtfmCOByt', &
      'yt offset of the center of buoyancy (COB) from the WAMIT reference point', ErrStat2, ErrMsg2, UnEchoLocal )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF

      ! ExctnMod  - Wave Excitation model {0: None, 1: DFT, 2: state-space} (switch)
      ! [STATE-SPACE REQUIRES *.ssexctn INPUT FILE]

   CALL ReadVar ( UnIn, FileName, InitInp%WAMIT%ExctnMod, 'ExctnMod', &
                                 'Wave Excitation model', ErrStat2, ErrMsg2, UnEchoLocal )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF

      ! RdtnMod  - Radiation memory-effect model {1: convolution, 2: state-space} (switch)
      ! [STATE-SPACE REQUIRES *.ss INPUT FILE]

  CALL ReadVar ( UnIn, FileName, InitInp%WAMIT%RdtnMod, 'RdtnMod', &
                                 'Radiation memory-effect model', ErrStat2, ErrMsg2, UnEchoLocal )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF


      ! RdtnTMax - Analysis time for wave radiation kernel calculations
      ! NOTE: Use RdtnTMax = 0.0 to eliminate wave radiation damping

   CALL ReadVar ( UnIn, FileName, InitInp%WAMIT%RdtnTMax, 'RdtnTMax', &
                                 'Analysis time for wave radiation kernel calculations', ErrStat2, ErrMsg2, UnEchoLocal )

      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF


      ! RdtnDT - Time step for wave radiation kernel calculations


   CALL ReadVar ( UnIn, FileName, InitInp%WAMIT%Conv_Rdtn%RdtnDTChr, 'RdtnDT', 'Time step for wave radiation kernel calculations', ErrStat2, ErrMsg2, UnEchoLocal )

      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF

   
!bjj: should we add this?
!test for numerical stability
!      IF ( FP_InitData%RdtnDT <= FP_InitData%RdtnTMax*EPSILON(FP_InitData%RdtnDT) )  THEN  ! Test RdtnDT and RdtnTMax to ensure numerical stability -- HINT: see the use of OnePlusEps."
!         ErrMsg  = ' RdtnDT must be greater than '//TRIM ( Num2LStr( RdtnTMax*EPSILON(RdtnDT) ) )//' seconds.'
!         ErrStat = ErrID_Fatal
!         CLOSE( UnIn )
!         RETURN
!      END IF



   !-------------------------------------------------------------------------------------------------
   ! Data section for 2nd order WAMIT forces
   !-------------------------------------------------------------------------------------------------


     ! Header

   CALL ReadCom( UnIn, FileName, '2nd order forces header (WAMIT2 module)', ErrStat2, ErrMsg2, UnEchoLocal )

      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF


        ! MnDrift    -- Mean drift forces computed from WAMIT file: {0: No mean drift, [7, 8, 9, 10, 11, or 12]: WAMIT file to use}

   CALL ReadVar ( UnIn, FileName, InitInp%WAMIT2%MnDrift, 'MnDrift', 'Mean drift forces computed from WAMIT file: {0: No mean drift, [7, 8, 9, 10, 11, or 12]: WAMIT file to use}', ErrStat2, ErrMsg2, UnEchoLocal )

      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF


        ! NewmanApp  -- Slow drift forces computed with Newman's approximation from  WAMIT file: {0: No mean drift, [7, 8, 9, 10, 11, or 12]: WAMIT file to use}

   CALL ReadVar ( UnIn, FileName, InitInp%WAMIT2%NewmanApp, 'NewmanApp', 'Mean drift forces computed from WAMIT file: {0: No mean drift, [7, 8, 9, 10, 11, or 12]: WAMIT file to use}', ErrStat2, ErrMsg2, UnEchoLocal )

      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF


        ! DiffQTF    -- Full Difference-Frequency forces computed with full QTFs from WAMIT file: {0: No difference-frequency forces, [10, 11, or 12]: WAMIT file to use} -- Only one of MnDrift, NewmanApp, or DiffQYT can be non-zero

   CALL ReadVar ( UnIn, FileName, InitInp%WAMIT2%DiffQTF, 'DiffQTF', 'Full Difference-Frequency forces computed with full QTFs from WAMIT file: '// &
       '{0: No difference-frequency forces, [10, 11, or 12]: WAMIT file to use} -- Only one of MnDrift, NewmanApp, or DiffQYT can be non-zero', ErrStat2, ErrMsg2, UnEchoLocal )

      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF


        ! SumQTF     -- Full        Sum-Frequency forces computed with full QTFs from WAMIT file: {0: No        Sum-frequency forces, [10, 11, or 12]: WAMIT file to use}

   CALL ReadVar ( UnIn, FileName, InitInp%WAMIT2%SumQTF, 'SumQTF', 'Full Sum-Frequency forces computed with full QTFs from WAMIT file: {0: No Sum-frequency forces, [10, 11, or 12]: WAMIT file to use}', ErrStat2, ErrMsg2, UnEchoLocal )

      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF



   !-------------------------------------------------------------------------------------------------
   ! Data section for Floating platform force flags
   !-------------------------------------------------------------------------------------------------

      ! Header

   CALL ReadCom( UnIn, FileName, 'Floating platform force flags header', ErrStat2, ErrMsg2, UnEchoLocal )

      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF


       ! PtfmSgFChr - Platform horizontal surge translation force flag

   CALL ReadVar ( UnIn, FileName, InitInp%PtfmSgFChr, 'PtfmSgFChr', 'Platform horizontal surge translation force flag', ErrStat2, ErrMsg2, UnEchoLocal )

      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF

   CALL Conv2UC( InitInp%PtfmSgFChr )    ! Convert Line to upper case.


      ! PtfmSwFChr - Platform horizontal sway translation force flag

   CALL ReadVar ( UnIn, FileName, InitInp%PtfmSwFChr, 'PtfmSwFChr', 'Platform horizontal sway translation force flag', ErrStat2, ErrMsg2, UnEchoLocal )

      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF

   CALL Conv2UC( InitInp%PtfmSwFChr )    ! Convert Line to upper case.


       ! PtfmHvFChr - Platform vertical heave translation force flag

   CALL ReadVar ( UnIn, FileName, InitInp%PtfmHvFChr, 'PtfmHvFChr', 'Platform vertical heave translation force flag', ErrStat2, ErrMsg2, UnEchoLocal )

      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF

   CALL Conv2UC( InitInp%PtfmHvFChr )    ! Convert Line to upper case.


        ! PtfmRFChr - Platform roll tilt rotation force flag

   CALL ReadVar ( UnIn, FileName, InitInp%PtfmRFChr, 'PtfmRFChr', 'Platform roll tilt rotation force flag', ErrStat2, ErrMsg2, UnEchoLocal )

      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF

   CALL Conv2UC( InitInp%PtfmRFChr )    ! Convert Line to upper case.


        ! PtfmPFChr - Platform pitch tilt rotation force flag

   CALL ReadVar ( UnIn, FileName, InitInp%PtfmPFChr, 'PtfmPFChr', 'Platform pitch tilt rotation force flag', ErrStat2, ErrMsg2, UnEchoLocal )

      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF

   CALL Conv2UC( InitInp%PtfmPFChr )    ! Convert Line to upper case.


        ! PtfmYFChr - Platform yaw rotation force flag

   CALL ReadVar ( UnIn, FileName, InitInp%PtfmYFChr, 'PtfmYFChr', 'Platform yaw rotation force flag', ErrStat2, ErrMsg2, UnEchoLocal )

      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF

   CALL Conv2UC( InitInp%PtfmYFChr )    ! Convert Line to upper case.



   !-------------------------------------------------------------------------------------------------
   ! Floating Platform Additional Stiffness and Damping Section
   !-------------------------------------------------------------------------------------------------


     ! Header

   CALL ReadCom( UnIn, FileName, 'Additional stiffness and damping header', ErrStat2, ErrMsg2, UnEchoLocal )

      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF


      ! AddF0 - Additional preload
      
   CALL ReadAry ( UnIn, FileName, InitInp%AddF0, 6, 'AddF0', &
                           ' Additional preload vector', ErrStat2,  ErrMsg2, UnEchoLocal )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF
   
   
      ! AddCLin

   DO I=1,6

      WRITE(strI,'(I1)') I
      CALL ReadAry ( UnIn, FileName, InitInp%AddCLin(I,:), 6, 'AddCLin', &
                           ' Row '//strI//' of the additional linear stiffness matrix', ErrStat2,  ErrMsg2, UnEchoLocal )

      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF
   END DO


       ! AddBLin

   DO I=1,6

      WRITE(strI,'(I1)') I
      CALL ReadAry ( UnIn, FileName, InitInp%AddBLin(I,:), 6, 'AddBLin', &
                           ' Row '//strI//' of the additional linear damping matrix', ErrStat2,  ErrMsg2, UnEchoLocal )

      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF
   END DO


       ! AddBQuad

   DO I=1,6

      WRITE(strI,'(I1)') I
      CALL ReadAry ( UnIn, FileName, InitInp%AddBQuad(I,:), 6, 'AddBQuad', &
                           ' Row '//strI//' of the additional quadratic damping matrix', ErrStat2,  ErrMsg2, UnEchoLocal )

      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF
   END DO


   !-------------------------------------------------------------------------------------------------
   !  Axial Coefficients Section
   !-------------------------------------------------------------------------------------------------
   
   
       ! Header
      
   CALL ReadCom( UnIn, FileName, 'Axial coefs header', ErrStat2, ErrMsg2, UnEchoLocal )
   
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF
   
   
      ! NAxCoef - Number of axial coefficients
   
   CALL ReadVar ( UnIn, FileName, InitInp%Morison%NAxCoefs, 'NAxCoefs', 'Number of axial coefficients', ErrStat2, ErrMsg2, UnEchoLocal )

      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF
   
         ! Table header
      
      CALL ReadCom( UnIn, FileName, 'Axial coefficient table header', ErrStat2, ErrMsg2, UnEchoLocal )
   
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF
      
         ! Table header
      
      CALL ReadCom( UnIn, FileName, 'Axial coefficient table header', ErrStat2, ErrMsg2, UnEchoLocal )
   
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF
   
   IF ( InitInp%Morison%NAxCoefs > 0 ) THEN
      
      
         ! Allocate memory for Axial Coef-related arrays
         
      ALLOCATE ( InitInp%Morison%AxialCoefs(InitInp%Morison%NAxCoefs), STAT = ErrStat2 )
      IF ( ErrStat2 /= 0 ) THEN      
         CALL SetErrStat( ErrID_Fatal, 'Error allocating space for AxialCoefs array.', ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
         CALL CleanUp()
         RETURN
      END IF
          
      DO I = 1,InitInp%Morison%NAxCoefs
            ! read the table entries   AxCoefID   CdAx  CaAx    in the HydroDyn input file
         READ(UnIn,'(A)',IOSTAT=ErrStat2) Line      !read into a line 
            
         IF (ErrStat2 == 0) THEN
            READ(Line,*,IOSTAT=ErrStat2) InitInp%Morison%AxialCoefs(I)%AxCoefID, InitInp%Morison%AxialCoefs(I)%AxCd, InitInp%Morison%AxialCoefs(I)%AxCa, InitInp%Morison%AxialCoefs(I)%AxCp
         END IF      
       
         IF ( ErrStat2 /= 0 ) THEN
            CALL SetErrStat( ErrID_Fatal, 'Failed to read axial coefficients.', ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
            CALL CleanUp()
            RETURN
         END IF 
         
         IF ( InitInp%Echo ) THEN
            WRITE( UnEchoLocal, '(A)' ) TRIM(Line)
         END IF
         
      END DO
      
   END IF
   
   
   !-------------------------------------------------------------------------------------------------
   ! Member Joints Section
   !-------------------------------------------------------------------------------------------------


      ! Header

   CALL ReadCom( UnIn, FileName, 'Member joints header', ErrStat2, ErrMsg2, UnEchoLocal )

      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF


      ! NJoints - Number of member joints

   CALL ReadVar ( UnIn, FileName, InitInp%Morison%NJoints, 'NJoints', 'Number of member joints', ErrStat2, ErrMsg2, UnEchoLocal )

      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF

         ! Table header

      CALL ReadCom( UnIn, FileName, 'Member joints table header', ErrStat2, ErrMsg2, UnEchoLocal )

      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF

         ! Table header

      CALL ReadCom( UnIn, FileName, 'Member joints table header', ErrStat2, ErrMsg2, UnEchoLocal )

      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF

   IF ( InitInp%Morison%NJoints > 0 ) THEN


         ! Allocate memory for Joint-related arrays

      ALLOCATE ( InitInp%Morison%InpJoints(InitInp%Morison%NJoints), STAT = ErrStat2 )
      IF ( ErrStat2 /= 0 ) THEN
         CALL SetErrStat( ErrID_Fatal, 'Error allocating space for InpJoints array.', ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
         CALL CleanUp()
         RETURN
      END IF

      DO I = 1,InitInp%Morison%NJoints
            ! read the table entries   JointID   Jointxi     Jointyi    Jointzi      JointAxID   JointOvrlp    in the HydroDyn input file
         READ(UnIn,'(A)',IOSTAT=ErrStat2) Line      !read into a line

         IF (ErrStat2 == 0) THEN
            READ(Line,*,IOSTAT=ErrStat2) InitInp%Morison%InpJoints(I)%JointID, InitInp%Morison%InpJoints(I)%JointPos(1), InitInp%Morison%InpJoints(I)%JointPos(2), InitInp%Morison%InpJoints(I)%JointPos(3), InitInp%Morison%InpJoints(I)%JointAxID, InitInp%Morison%InpJoints(I)%JointOvrlp
         END IF

         IF ( ErrStat2 /= 0 ) THEN
            CALL SetErrStat( ErrID_Fatal, 'Failed to read joints.', ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
            CALL CleanUp()
            RETURN
         END IF

         IF ( InitInp%Echo ) THEN
            WRITE( UnEchoLocal, '(A)' ) TRIM(Line)
         END IF

      END DO

   END IF




   !-------------------------------------------------------------------------------------------------
   ! Member Cross-section Properties Section
   !-------------------------------------------------------------------------------------------------


      ! Header

   CALL ReadCom( UnIn, FileName, 'Member cross-section properties header', ErrStat2, ErrMsg2, UnEchoLocal )

      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF


      ! NPropSets - Number of member cross-section property sets

   CALL ReadVar ( UnIn, FileName, InitInp%Morison%NPropSets, 'NPropSets', 'Number of member cross-section property sets', ErrStat2, ErrMsg2, UnEchoLocal )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF

      ! Table header

   CALL ReadCom( UnIn, FileName, 'Member cross-section properties table header', ErrStat2, ErrMsg2, UnEchoLocal )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF

      ! Table header

   CALL ReadCom( UnIn, FileName, 'Member cross-section properties table header', ErrStat2, ErrMsg2, UnEchoLocal )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF

   IF ( InitInp%Morison%NPropSets > 0 ) THEN


         ! Allocate memory for Member cross-section property set-related arrays

      ALLOCATE ( InitInp%Morison%MPropSets(InitInp%Morison%NPropSets), STAT = ErrStat2 )
      IF ( ErrStat2 /= 0 ) THEN
         CALL SetErrStat( ErrID_Fatal, 'Error allocating space for MPropSets array.', ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
         CALL CleanUp()
         RETURN
      END IF


      DO I = 1,InitInp%Morison%NPropSets

         READ(UnIn,'(A)',IOSTAT=ErrStat2) Line      !read into a line

         IF (ErrStat2 == 0) THEN
            READ(Line,*,IOSTAT=ErrStat) InitInp%Morison%MPropSets(I)%PropSetID, InitInp%Morison%MPropSets(I)%PropD, InitInp%Morison%MPropSets(I)%PropThck
         END IF

         IF ( ErrStat2 /= 0 ) THEN
            CALL SetErrStat( ErrID_Fatal, 'Failed to read member cross-section properties.', ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
            CALL CleanUp()
            RETURN
         END IF
            
         IF ( InitInp%Echo ) THEN
            WRITE( UnEchoLocal, '(A)' ) TRIM(Line)
         END IF

      END DO

   END IF





   !-------------------------------------------------------------------------------------------------
   ! Simple hydrodynamic coefficients Section
   !-------------------------------------------------------------------------------------------------


      ! Header

   CALL ReadCom( UnIn, FileName, 'Simple hydrodynamic coefficients header', ErrStat2, ErrMsg2, UnEchoLocal )

      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF

      ! Header

   CALL ReadCom( UnIn, FileName, 'Simple hydrodynamic coefficients table header', ErrStat2, ErrMsg2, UnEchoLocal )

      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF

      ! Header

   CALL ReadCom( UnIn, FileName, 'Simple hydrodynamic coefficients table header', ErrStat2, ErrMsg2, UnEchoLocal )

      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF


   READ(UnIn,'(A)',IOSTAT=ErrStat2) Line      !read into a line

   IF (ErrStat2 == 0) THEN
      READ(Line,*,IOSTAT=ErrStat2) InitInp%Morison%SimplCd, InitInp%Morison%SimplCdMG, InitInp%Morison%SimplCa, InitInp%Morison%SimplCaMG, InitInp%Morison%SimplCp, InitInp%Morison%SimplCpMG, InitInp%Morison%SimplAxCa, InitInp%Morison%SimplAxCaMG, InitInp%Morison%SimplAxCp, InitInp%Morison%SimplAxCpMG
   END IF

   IF ( ErrStat2 /= 0 ) THEN
      CALL SetErrStat( ErrID_Fatal, 'Failed to read simple hydrodynamic coefficients.', ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      CALL CleanUp()
      RETURN
   END IF

   IF ( InitInp%Echo ) THEN
      WRITE( UnEchoLocal, '(A)' ) TRIM(Line)
   END IF




   !-------------------------------------------------------------------------------------------------
   ! Depth-based Hydrodynamic Coefficients Section
   !-------------------------------------------------------------------------------------------------


      ! Header

   CALL ReadCom( UnIn, FileName, 'Depth-based hydrodynamic coefficients header', ErrStat2, ErrMsg2, UnEchoLocal )

      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF


      ! NCoefDpth - Number of depth-based hydrodynamic coefficient property sets

   CALL ReadVar ( UnIn, FileName, InitInp%Morison%NCoefDpth, 'NCoefDpth', 'Number of depth-based hydrodynamic coefficient property sets', ErrStat2, ErrMsg2, UnEchoLocal )

      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF

      ! Table header

   CALL ReadCom( UnIn, FileName, 'Depth-based hydrodynamic coefficients table header', ErrStat2, ErrMsg2, UnEchoLocal )

      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF

      ! Table header

   CALL ReadCom( UnIn, FileName, 'Depth-based hydrodynamic coefficients table header', ErrStat2, ErrMsg2, UnEchoLocal )

      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF


   IF ( InitInp%Morison%NCoefDpth > 0 ) THEN


         ! Allocate memory for depth-based coefficient arrays

      ALLOCATE ( InitInp%Morison%CoefDpths(InitInp%Morison%NCoefDpth), STAT = ErrStat2 )
      IF ( ErrStat2 /= 0 ) THEN
         CALL SetErrStat( ErrID_Fatal, 'Error allocating space for CoefDpths array.', ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
         CALL CleanUp()
         RETURN
      END IF
                  
      DO I = 1,InitInp%Morison%NCoefDpth

         READ(UnIn,'(A)',IOSTAT=ErrStat2) Line      !read into a line

         IF (ErrStat2 == 0) THEN
            READ(Line,*,IOSTAT=ErrStat2) InitInp%Morison%CoefDpths(I)%Dpth, InitInp%Morison%CoefDpths(I)%DpthCd, InitInp%Morison%CoefDpths(I)%DpthCdMG, &
                                         InitInp%Morison%CoefDpths(I)%DpthCa, InitInp%Morison%CoefDpths(I)%DpthCaMG, InitInp%Morison%CoefDpths(I)%DpthCp, InitInp%Morison%CoefDpths(I)%DpthCpMG, &
                                         InitInp%Morison%CoefDpths(I)%DpthAxCa, InitInp%Morison%CoefDpths(I)%DpthAxCaMG, InitInp%Morison%CoefDpths(I)%DpthAxCp, InitInp%Morison%CoefDpths(I)%DpthAxCpMG
         END IF

         IF (ErrStat2 /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, 'Failed to read depth-based coefficient array.', ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
            CALL CleanUp()
            RETURN
         END IF

         IF ( InitInp%Echo ) THEN
            WRITE( UnEchoLocal, '(A)' ) TRIM(Line)
         END IF

      END DO

   END IF


   !-------------------------------------------------------------------------------------------------
   ! Member-based Hydrodynamic Coefficients Section
   !-------------------------------------------------------------------------------------------------


      ! Header

   CALL ReadCom( UnIn, FileName, 'Member-based hydrodynamic coefficients header', ErrStat2, ErrMsg2, UnEchoLocal )

      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF


      ! NCoefMembers - Number of member-based hydrodynamic coefficient property sets

   CALL ReadVar ( UnIn, FileName, InitInp%Morison%NCoefMembers, 'NCoefMembers', 'Number of member-based hydrodynamic coefficient property sets', ErrStat2, ErrMsg2, UnEchoLocal )

      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF


      ! Table header

   CALL ReadCom( UnIn, FileName, 'Member-based hydrodynamic coefficients table header', ErrStat2, ErrMsg2, UnEchoLocal )

      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF

      ! Table header

   CALL ReadCom( UnIn, FileName, 'Member-based hydrodynamic coefficients table header', ErrStat2, ErrMsg2, UnEchoLocal )

      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF


   IF ( InitInp%Morison%NCoefMembers > 0 ) THEN


         ! Allocate memory for Member-based coefficient arrays

      ALLOCATE ( InitInp%Morison%CoefMembers(InitInp%Morison%NCoefMembers), STAT = ErrStat2 )
      IF ( ErrStat2 /= 0 ) THEN
         CALL SetErrStat( ErrID_Fatal, 'Error allocating space for CoefMembers array.', ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
         CALL CleanUp()
         RETURN
      END IF

      DO I = 1,InitInp%Morison%NCoefMembers

         READ(UnIn,'(A)',IOSTAT=ErrStat2) Line      !read into a line

         IF (ErrStat2 == 0) THEN
            READ(Line,*,IOSTAT=ErrStat2) InitInp%Morison%CoefMembers(I)%MemberID,      &
                                         InitInp%Morison%CoefMembers(I)%MemberCd1,     InitInp%Morison%CoefMembers(I)%MemberCd2,     &
                                         InitInp%Morison%CoefMembers(I)%MemberCdMG1,   InitInp%Morison%CoefMembers(I)%MemberCdMG2,   &
                                         InitInp%Morison%CoefMembers(I)%MemberCa1,     InitInp%Morison%CoefMembers(I)%MemberCa2,     &
                                         InitInp%Morison%CoefMembers(I)%MemberCaMG1,   InitInp%Morison%CoefMembers(I)%MemberCaMG2,   &
                                         InitInp%Morison%CoefMembers(I)%MemberCp1,     InitInp%Morison%CoefMembers(I)%MemberCp2,     &
                                         InitInp%Morison%CoefMembers(I)%MemberCpMG1,   InitInp%Morison%CoefMembers(I)%MemberCpMG2,   &
                                         InitInp%Morison%CoefMembers(I)%MemberAxCa1,   InitInp%Morison%CoefMembers(I)%MemberAxCa2,   &
                                         InitInp%Morison%CoefMembers(I)%MemberAxCaMG1, InitInp%Morison%CoefMembers(I)%MemberAxCaMG2, &
                                         InitInp%Morison%CoefMembers(I)%MemberAxCp1,   InitInp%Morison%CoefMembers(I)%MemberAxCp2,   &
                                         InitInp%Morison%CoefMembers(I)%MemberAxCpMG1, InitInp%Morison%CoefMembers(I)%MemberAxCpMG2
         END IF

       
         IF ( ErrStat2 /= 0 ) THEN
            CALL SetErrStat( ErrID_Fatal, 'Failed to read member cross-section properties.', ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
            CALL CleanUp()
            RETURN
         END IF

         IF ( InitInp%Echo ) THEN
            WRITE( UnEchoLocal, '(A)' ) TRIM(Line)
         END IF

      END DO

   END IF



   !-------------------------------------------------------------------------------------------------
   ! Members Section
   !-------------------------------------------------------------------------------------------------


      ! Header

   CALL ReadCom( UnIn, FileName, 'Members header', ErrStat2, ErrMsg2, UnEchoLocal )

      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF


      ! NMembers - Number of members in the input file

   CALL ReadVar ( UnIn, FileName, InitInp%Morison%NMembers, 'NMembers', 'Number of members', ErrStat2, ErrMsg2, UnEchoLocal )

      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF

      ! Table header

   CALL ReadCom( UnIn, FileName, 'Members table header', ErrStat2, ErrMsg2, UnEchoLocal )

      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF

      ! Table header

   CALL ReadCom( UnIn, FileName, 'Members table header', ErrStat2, ErrMsg2, UnEchoLocal )

      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF


   IF ( InitInp%Morison%NMembers > 0 ) THEN


         ! Allocate memory for Members arrays

      ALLOCATE ( InitInp%Morison%InpMembers(InitInp%Morison%NMembers), STAT = ErrStat2 )
      IF ( ErrStat2 /= 0 ) THEN         
         CALL SetErrStat( ErrID_Fatal, 'Error allocating space for InpMembers array.', ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
         CALL CleanUp()
         RETURN
      END IF

      DO I = 1,InitInp%Morison%NMembers
         !ReadStr ( UnIn, Fil, Line, 'Joint table', VarDescr, ErrStat )
         READ(UnIn,'(A)',IOSTAT=ErrStat2) Line      !read into a line

         IF (ErrStat2 == 0) THEN
            READ(Line,*,IOSTAT=ErrStat2) InitInp%Morison%InpMembers(I)%MemberID,   InitInp%Morison%InpMembers(I)%MJointID1,    &
                                        InitInp%Morison%InpMembers(I)%MJointID2,   InitInp%Morison%InpMembers(I)%MPropSetID1,  &
                                        InitInp%Morison%InpMembers(I)%MPropSetID2, InitInp%Morison%InpMembers(I)%MDivSize,     &
                                        InitInp%Morison%InpMembers(I)%MCoefMod,    InitInp%Morison%InpMembers(I)%PropPot
         END IF

         IF ( ErrStat2 /= 0 ) THEN         
            CALL SetErrStat( ErrID_Fatal, 'Failed to read member properties.', ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
            CALL CleanUp()
            RETURN
         END IF

         IF ( InitInp%Echo ) THEN
            WRITE( UnEchoLocal, '(A)' ) TRIM(Line)
         END IF

      END DO

   END IF


   !-------------------------------------------------------------------------------------------------
   ! Filled Members Section
   !-------------------------------------------------------------------------------------------------


      ! Header

   CALL ReadCom( UnIn, FileName, 'Filled members header', ErrStat2, ErrMsg2, UnEchoLocal )

      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF


      ! NFillGroups - Number of fill groups

   CALL ReadVar ( UnIn, FileName, InitInp%Morison%NFillGroups, 'NFillGroups', 'Number of fill groups', ErrStat2, ErrMsg2, UnEchoLocal )

      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF

      ! Table header

   CALL ReadCom( UnIn, FileName, 'Fill groups table header', ErrStat2, ErrMsg2, UnEchoLocal )

      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF

      ! Table header

   CALL ReadCom( UnIn, FileName, 'Fill groups table header', ErrStat2, ErrMsg2, UnEchoLocal )

      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF


   IF ( InitInp%Morison%NFillGroups > 0 ) THEN


         ! Allocate memory for filled group arrays

      ALLOCATE ( InitInp%Morison%FilledGroups(InitInp%Morison%NFillGroups), STAT = ErrStat2 )
      IF ( ErrStat2 /= 0 ) THEN
         CALL SetErrStat( ErrID_Fatal, 'Error allocating space for FilledGroups array.', ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
         CALL CleanUp()
         RETURN
      END IF

      DO I = 1,InitInp%Morison%NFillGroups

         READ(UnIn,'(A)',IOSTAT=ErrStat2) Line      !read into a line

         IF (ErrStat2 == 0) THEN

            READ(Line,*,IOSTAT=ErrStat2) InitInp%Morison%FilledGroups(I)%FillNumM
            IF ( ErrStat2 /= 0 ) THEN
               CALL SetErrStat( ErrID_Fatal, 'Failed to read FillNumM.', ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
               CALL CleanUp()
               RETURN
            END IF


            ALLOCATE ( InitInp%Morison%FilledGroups(I)%FillMList(InitInp%Morison%FilledGroups(I)%FillNumM), STAT = ErrStat2 )
            IF ( ErrStat2 /= 0 ) THEN
               CALL SetErrStat( ErrID_Fatal, 'Error allocating space for FillMList array.', ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
               CALL CleanUp()
               RETURN
            END IF


            READ(Line,*,IOSTAT=ErrStat2) InitInp%Morison%FilledGroups(I)%FillNumM,  InitInp%Morison%FilledGroups(I)%FillMList,   &
                                         InitInp%Morison%FilledGroups(I)%FillFSLoc, InitInp%Morison%FilledGroups(I)%FillDensChr

            IF ( ErrStat2 /= 0 ) THEN
               CALL SetErrStat( ErrID_Fatal, 'Failed to read filled group properties.', ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
               CALL CleanUp()
               RETURN
            END IF
            
            IF ( InitInp%Echo ) THEN
               WRITE( UnEchoLocal, '(A)' ) TRIM(Line)
            END IF

         END IF

      END DO

   END IF


   !-------------------------------------------------------------------------------------------------
   ! Marine Growth by Depth Section
   !-------------------------------------------------------------------------------------------------


      ! Header

   CALL ReadCom( UnIn, FileName, 'Marine growth by depth header', ErrStat2, ErrMsg2, UnEchoLocal )

      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF



      ! NMGDepths - Number marine growth depths

   CALL ReadVar ( UnIn, FileName, InitInp%Morison%NMGDepths, 'NMGDepths', 'Number marine growth depths', ErrStat2, ErrMsg2, UnEchoLocal )

      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF

      ! Table header

   CALL ReadCom( UnIn, FileName, 'Marine growth by depth table header', ErrStat2, ErrMsg2, UnEchoLocal )

      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF

      ! Table header

   CALL ReadCom( UnIn, FileName, 'Marine growth by depth table header', ErrStat2, ErrMsg2, UnEchoLocal )

      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF


   IF ( InitInp%Morison%NMGDepths > 0 ) THEN


         ! Allocate memory for marine growth depths array

      ALLOCATE ( InitInp%Morison%MGDepths(InitInp%Morison%NMGDepths), STAT = ErrStat2 )
      IF ( ErrStat2 /= 0 ) THEN
         CALL SetErrStat( ErrID_Fatal, 'Error allocating space for MGDepths array.', ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
         CALL CleanUp()
         RETURN
      END IF

      DO I = 1,InitInp%Morison%NMGDepths

         READ(UnIn,'(A)',IOSTAT=ErrStat2) Line      !read into a line

         IF (ErrStat2 == 0) THEN
            READ(Line,*,IOSTAT=ErrStat2) InitInp%Morison%MGDepths(I)%MGDpth, InitInp%Morison%MGDepths(I)%MGThck, InitInp%Morison%MGDepths(I)%MGDens
         END IF

         IF ( ErrStat2 /= 0 ) THEN
            CALL SetErrStat( ErrID_Fatal, 'Failed to read marine growth depth properties.', ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
            CALL CleanUp()
            RETURN
         END IF
         
         IF ( InitInp%Echo ) THEN
            WRITE( UnEchoLocal, '(A)' ) TRIM(Line)
         END IF

      END DO

   END IF


   !-------------------------------------------------------------------------------------------------
   ! Member Output List Section
   !-------------------------------------------------------------------------------------------------

      ! Header

   CALL ReadCom( UnIn, FileName, 'Member output list header', ErrStat2, ErrMsg2, UnEchoLocal )

      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF


      ! NMOutputs - Number of members to output

   CALL ReadVar ( UnIn, FileName, InitInp%Morison%NMOutputs, 'NMOutputs', 'Number of members to output', ErrStat2, ErrMsg2, UnEchoLocal )

      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF

      ! Table header

   CALL ReadCom( UnIn, FileName, 'Member output list table header', ErrStat2, ErrMsg2, UnEchoLocal )

      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF

      ! Table header

   CALL ReadCom( UnIn, FileName, 'Member output list table header', ErrStat2, ErrMsg2, UnEchoLocal )

      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF


   IF ( InitInp%Morison%NMOutputs > 0 ) THEN


         ! Allocate memory for filled group arrays

      ALLOCATE ( InitInp%Morison%MOutLst(InitInp%Morison%NMOutputs), STAT = ErrStat2 )
      IF ( ErrStat2 /= 0 ) THEN
         CALL SetErrStat( ErrID_Fatal, 'Error allocating space for MOutLst array.', ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
         CALL CleanUp()
         RETURN
      END IF      


      DO I = 1,InitInp%Morison%NMOutputs

         READ(UnIn,'(A)',IOSTAT=ErrStat2) Line      !read into a line

         IF (ErrStat2 == 0) THEN

            READ(Line,*,IOSTAT=ErrStat2) InitInp%Morison%MOutLst(I)%MemberID, InitInp%Morison%MOutLst(I)%NOutLoc
            IF ( ErrStat2 /= 0 ) THEN
               CALL SetErrStat( ErrID_Fatal, 'Failed to read NOutLoc.', ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
               CALL CleanUp()
               RETURN
            END IF      

            
            ALLOCATE ( InitInp%Morison%MOutLst(I)%NodeLocs(InitInp%Morison%MOutLst(I)%NOutLoc), STAT = ErrStat2 )
            IF ( ErrStat2 /= 0 ) THEN
               CALL SetErrStat( ErrID_Fatal, 'Error allocating space for NodeLocs array.', ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
               CALL CleanUp()
               RETURN
            END IF      


            ALLOCATE ( InitInp%Morison%MOutLst(I)%Marker1(InitInp%Morison%MOutLst(I)%NOutLoc), STAT = ErrStat2 )
            IF ( ErrStat2 /= 0 ) THEN
               CALL SetErrStat( ErrID_Fatal, 'Error allocating space for Marker1 array.', ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
               CALL CleanUp()
               RETURN
            END IF      


            ALLOCATE ( InitInp%Morison%MOutLst(I)%Marker2(InitInp%Morison%MOutLst(I)%NOutLoc), STAT = ErrStat2 )
            IF ( ErrStat2 /= 0 ) THEN
               CALL SetErrStat( ErrID_Fatal, 'Error allocating space for Marker2 array.', ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
               CALL CleanUp()
               RETURN
            END IF      

            ALLOCATE ( InitInp%Morison%MOutLst(I)%s(InitInp%Morison%MOutLst(I)%NOutLoc), STAT = ErrStat2 )

            IF ( ErrStat2 /= 0 ) THEN
               CALL SetErrStat( ErrID_Fatal, 'Error allocating space for s array.', ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
               CALL CleanUp()
               RETURN
            END IF      

            READ(Line,*,IOSTAT=ErrStat2) InitInp%Morison%MOutLst(I)%MemberID,  InitInp%Morison%MOutLst(I)%NOutLoc,  &
                                         InitInp%Morison%MOutLst(I)%NodeLocs

            IF ( ErrStat2 /= 0 ) THEN
               CALL SetErrStat( ErrID_Fatal, 'Failed to read member output list properties.', ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
               CALL CleanUp()
               RETURN
            END IF

            IF ( InitInp%Echo ) THEN
               WRITE( UnEchoLocal, '(A)' ) TRIM(Line)
            END IF

         END IF

      END DO

   END IF

   !-------------------------------------------------------------------------------------------------
   ! Joint Output List Section
   !-------------------------------------------------------------------------------------------------

      ! Header

   CALL ReadCom( UnIn, FileName, 'Joint output list header', ErrStat2, ErrMsg2, UnEchoLocal )

      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF


      ! NJOutputs - Number of joints to output

   CALL ReadVar ( UnIn, FileName, InitInp%Morison%NJOutputs, 'NJOutputs', 'Number of joints to output', ErrStat2, ErrMsg2, UnEchoLocal )

      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF

   IF ( InitInp%Morison%NJOutputs > 0 ) THEN

      ALLOCATE ( InitInp%Morison%JOutLst(InitInp%Morison%NJOutputs), STAT = ErrStat2 )

      IF ( ErrStat2 /= 0 ) THEN
         CALL SetErrStat( ErrID_Fatal, 'Error allocating space for JOutLst data structures.', ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
         CALL CleanUp()
         RETURN
      END IF      

      CALL AllocAry( tmpArray, InitInp%Morison%NJOutputs, 'temporary array for Joint outputs', ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
         IF (ErrStat >= AbortErrLev) THEN
            CALL CleanUp()
            RETURN
         END IF

      CALL ReadAry ( UnIn, FileName, tmpArray, InitInp%Morison%NJOutputs, 'JOutLst', 'Joint output list', ErrStat2,  ErrMsg2, UnEchoLocal )

      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF

      DO I = 1,InitInp%Morison%NJOutputs

         InitInp%Morison%JOutLst(I)%JointID = tmpArray(I)

      END DO

      DEALLOCATE(tmpArray)   
      
   ELSE
      
      ! There are no Joint Outputs, but there is a line to be parsed in the input file!
      
      ALLOCATE ( tmpArray(1), STAT = ErrStat2 )
      IF ( ErrStat2 /= 0 ) THEN
         CALL SetErrStat( ErrID_Fatal, 'Error allocating space for temporary array for Joint outputs.', ErrStat, ErrMsg, 'HydroDynInput_GetInput' )         
         CALL Cleanup()         
         RETURN
      END IF
      
      CALL ReadAry ( UnIn, FileName, tmpArray, 1, 'JOutLst', 'Joint output list', ErrStat2,  ErrMsg2, UnEchoLocal )      
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
         IF (ErrStat >= AbortErrLev) THEN
            CALL CleanUp()
            RETURN
         END IF
      
      DEALLOCATE(tmpArray)
      !we just want to read the line for echoing purposes when there are actually 0 Joint Outputs.
      
   END IF
   
   !-------------------------------------------------------------------------------------------------
   ! Data section for OUTPUT
   !-------------------------------------------------------------------------------------------------

      ! Header

   CALL ReadCom( UnIn, FileName, 'Output header', ErrStat2, ErrMsg2, UnEchoLocal )

      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF


         ! HDSum - Whether or not to generate a summary file

   CALL ReadVar ( UnIn, FileName, InitInp%HDSum, 'HDSum', 'Generate a HydroDyn summary file', ErrStat2, ErrMsg2, UnEchoLocal )

      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF


         ! OutAll - Whether or not to output information for every member and joint

   CALL ReadVar ( UnIn, FileName, InitInp%OutAll, 'OutAll', 'Generate all member and joint outputs', ErrStat2, ErrMsg2, UnEchoLocal )

      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF


         ! OutSwtch - Specify how to write to an output file

   CALL ReadVar ( UnIn, FileName, InitInp%OutSwtch, 'OutSwtch', 'Specify how to write to an output file', ErrStat2, ErrMsg2, UnEchoLocal )

      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF


        ! OutFmt - Format for numerical outputs

   CALL ReadVar ( UnIn, FileName, InitInp%OutFmt, 'OutFmt', 'Format for numerical outputs', ErrStat2, ErrMsg2, UnEchoLocal )

      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF

         ! OutSFmt - Format for output column headers

   CALL ReadVar ( UnIn, FileName, InitInp%OutSFmt, 'OutSFmt', 'Format for output column headers', ErrStat2, ErrMsg2, UnEchoLocal )

      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF


   !-------------------------------------------------------------------------------------------------
   ! Data section for FLOATING PLATFORM OUTPUTS
   !-------------------------------------------------------------------------------------------------

      ! Header

   !CALL ReadCom( UnIn, FileName, 'Floating Platform Outputs header', ErrStat2, ErrMsg2, UnEchoLocal )
   !
   !   CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
   !   IF (ErrStat >= AbortErrLev) THEN
   !      CALL CleanUp()
   !      RETURN
   !   END IF
   CALL ReadCom( UnIn, FileName, 'Outputs header', ErrStat2, ErrMsg2, UnEchoLocal )
   
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF
      
         ! OutList - list of requested parameters to output to a file

   !CALL ReadOutputList ( UnIn, FileName, InitInp%WAMIT%OutList, InitInp%WAMIT%NumOuts, &
   !                                           'OutList', 'List of floating platform outputs requested', ErrStat2, ErrMsg2, UnEchoLocal )
   ALLOCATE( InitInp%UserOutputs(2778), Stat=ErrStat2)  !todo: bjj: what is this 2778? 
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat( ErrID_Fatal, 'Error allocating UserOutputs.', ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
         CALL CleanUp()
         RETURN
      END IF
   
   CALL ReadOutputList ( UnIn, FileName, InitInp%UserOutputs, InitInp%NUserOutputs, &
                                              'OutList', 'List of user requested outputs', ErrStat2, ErrMsg2, UnEchoLocal )
   
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF

   !
   !
   !!-------------------------------------------------------------------------------------------------
   !! Data section for MESH-BASED OUTPUTS
   !!-------------------------------------------------------------------------------------------------
   !
   !   ! Header
   !   
   !CALL ReadCom( UnIn, FileName, 'Mesh-based Outputs header', ErrStat2, ErrMsg2, UnEchoLocal )
   !
   !   CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
   !   IF (ErrStat >= AbortErrLev) THEN
   !      CALL CleanUp()
   !      RETURN
   !   END IF
   !
   !      ! OutList - list of requested parameters to output to a file
   !
   !CALL ReadOutputList ( UnIn, FileName, InitInp%Morison%OutList, InitInp%Morison%NumOuts, &
   !                                           'OutList', 'List of mesh-based outputs requested', ErrStat2, ErrMsg2, UnEchoLocal )
   !
   !CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
   !IF (ErrStat >= AbortErrLev) THEN
   !   CALL CleanUp()
   !   RETURN
   !END IF
   
   !-------------------------------------------------------------------------------------------------
   ! This is the end of the input file
   !-------------------------------------------------------------------------------------------------
   
   CALL Cleanup()

   RETURN

CONTAINS
   !..............................
   SUBROUTINE Cleanup()
   
      IF (ALLOCATED(tmpArray)) DEALLOCATE(tmpArray)         
   
         ! Close input file
      CLOSE ( UnIn )

         ! Cleanup the Echo file and global variables
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      
      
   END SUBROUTINE Cleanup
   

END SUBROUTINE HydroDynInput_GetInput



  

!====================================================================================================
SUBROUTINE HydroDynInput_ProcessInitData( InitInp, ErrStat, ErrMsg )
!     This private subroutine verifies the input required for HydroDyn is correctly specified.
!----------------------------------------------------------------------------------------------------


      ! Passed variables

   TYPE(HydroDyn_InitInputType),  INTENT( INOUT )   :: InitInp              ! the hydrodyn data
   INTEGER,                       INTENT(   OUT )   :: ErrStat              ! returns a non-zero value when an error occurs
   CHARACTER(*),                  INTENT(   OUT )   :: ErrMsg               ! Error message if ErrStat /= ErrID_None

   INTEGER                                          :: I                    ! Generic loop counter index
   INTEGER                                          :: J                    ! Generic loop counter index
   INTEGER                                          :: K                    ! Generic loop counter index
   CHARACTER(1024)                                  :: TmpPath              ! Temporary storage for relative path name
   LOGICAL                                          :: FoundID              ! Boolean flag indicating whether an ID from one tables is found in one of the other input table
   REAL(ReKi)                                       :: MinDepth             ! The minimum depth entry in the Depth-based Hydrodynamic coefficents table
   REAL(ReKi)                                       :: MaxDepth             ! The maximum depth entry in the Depth-based Hydrodynamic coefficents table
   REAL(ReKi)                                       :: z1
   REAL(ReKi)                                       :: z2
   REAL(ReKi)                                       :: MinMembrDpth
   REAL(ReKi)                                       :: MaxMembrDpth
!   CHARACTER(10), ALLOCATABLE                       :: tmpOutLst(:)         !
   CHARACTER(3)                                     :: TmpExtension         ! Temporary variable for holding the file extension for 10d, 11d, 12d, 10s, 11s, 12s WAMIT files
   LOGICAL                                          :: TmpFileExist         ! Temporary variable in checking the existance of an input file.
   LOGICAL                                          :: JointUsed
   REAL(ReKi)                                       :: l
   REAL(ReKi)                                       :: lvec(3)
   LOGICAL, ALLOCATABLE                             :: foundMask(:)
   INTEGER                                          :: WaveModIn
   
   INTEGER(IntKi)                                   :: ErrStat2, IOS
   CHARACTER(1024)                                  :: ErrMsg2
   CHARACTER(*), PARAMETER                          :: RoutineName = 'HydroDynInput_ProcessInitData'
      ! Initialize ErrStat
         
   ErrStat = ErrID_None         
   ErrStat2 = ErrID_None
   ErrMsg  = ""    
   ErrMsg2  = ""
      
     
      !-------------------------------------------------------------------------
      ! Check environmental conditions
      !-------------------------------------------------------------------------

 
      ! WtrDens - Water density.

   IF ( InitInp%Waves%WtrDens < 0.0 )  THEN
      CALL SetErrStat( ErrID_Fatal,'WtrDens must not be negative.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF


      ! WtrDpth - Water depth

   IF ( InitInp%Morison%WtrDpth <= 0.0 )  THEN
      CALL SetErrStat( ErrID_Fatal,'WtrDpth must be greater than zero.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF


      ! MSL2SWL - Mean sea level to still water level

   IF ( InitInp%PotMod == 1 .AND. .NOT. EqualRealNos(InitInp%Morison%MSL2SWL, 0.0_ReKi) ) THEN
      CALL SetErrStat( ErrID_Fatal,'MSL2SWL must be 0 when PotMod = 1 (WAMIT).',ErrStat,ErrMsg,RoutineName)        
      RETURN
   END IF
   
   !IF ( .NOT. EqualRealNos(InitInp%Morison%MSL2SWL, 0.0_ReKi) ) THEN  !TODO  Alter this check when we support MSL2SWL
   !   CALL SetErrStat( ErrID_Fatal,'MSL2SWL must be 0. Future versions of HydroDyn will once again support any value of MSL2SWL.'
   !   RETURN
   !END IF
      
   
      ! WaveMod - Wave kinematics model switch.

   IF ( LEN_TRIM(InitInp%Waves%WaveModChr) > 1 ) THEN

      IF ( InitInp%Waves%WaveModChr(1:2) == '1P' )  THEN                     ! The user wants to specify the phase in place of a random phase

         READ (InitInp%Waves%WaveModChr(3:),*,IOSTAT=IOS )  InitInp%Waves%WavePhase
            CALL CheckIOS ( IOS, "", 'WavePhase', NumType, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg,RoutineName)
            IF ( ErrStat >= AbortErrLev ) RETURN
            
         WaveModIn               = 1
         InitInp%Waves%WaveMod   = 10                                ! Internally define WaveMod = 10 to mean regular waves with a specified (nonrandom) phase
         InitInp%Waves%WavePhase = InitInp%Waves%WavePhase*D2R       ! Convert the phase from degrees to radians

      ELSE                                               ! The user must have specified WaveMod incorrectly.
         CALL SetErrStat( ErrID_Fatal,'WaveMod incorrectly specified',ErrStat,ErrMsg,RoutineName)
         RETURN
      END IF

   ELSE
         ! The line below only works for 1 digit reads
      READ( InitInp%Waves%WaveModChr, *, IOSTAT=IOS ) InitInp%Waves%WaveMod
         CALL CheckIOS ( IOS, "", 'WaveMod', NumType, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) RETURN
         
      WaveModIn               = InitInp%Waves%WaveMod

   END IF ! LEN_TRIM(InitInp%Waves%WaveModChr)

   IF ( (WaveModIn == 6) .AND. .NOT. EqualRealNos(InitInp%Morison%MSL2SWL, 0.0_ReKi) ) THEN
      CALL SetErrStat( ErrID_Fatal,'MSL2SWL must be 0 when WaveMod = 6.',ErrStat,ErrMsg,RoutineName)        
      RETURN
   END IF
   

   IF ( WaveModIn < 0 .OR. WaveModIn > 6 ) THEN
      IF ( InitInp%PotMod == 1  ) THEN
         CALL SetErrStat( ErrID_Fatal,'WaveMod must be 0, 1, 1P#, 2, 3, 4, 5, or 6.',ErrStat,ErrMsg,RoutineName)
         RETURN
!ADP: This seems like a strange test on ErrStat...
      ELSE IF ( ErrStat /= ErrID_None .OR. WaveModIn /= 5)  THEN
         CALL SetErrStat( ErrID_Fatal,'WaveMod must be 0, 1, 1P#, 2, 3, 4, or 5.',ErrStat,ErrMsg,RoutineName)
         RETURN
      END IF
   END IF

      ! Linearization Checks
   ! LIN-TODO:
   !errors if:
   !if (                                                                   &
   !     (WaveModIn /= 0)                                             .or. &
   !     (InitInp%Waves2%WvDiffQTFF /= .false.)                       .or. &
   !     (InitInp%Waves2%WvSumQTFF /= .false.)                        .or. &
   !     (InitInp%PotMod /= 0 .or. InitInp%PotMod /=1)                .or. &
   !     (InitInp%WAMIT%ExctnMod /=0 .or. InitInp%WAMIT%ExctnMod /=2) .or. &
   !     (InitInp%WAMIT%RdtnMod  /=0 .or. InitInp%WAMIT%RdtnMod  /=2) .or. &
   !     (InitInp%WAMIT2%MnDrift /=0)                                 .or. &
   !     (InitInp%WAMIT2%NewmanApp /= 0)                              .or. &
   !     (InitInp%WAMIT2%SumQTF /= 0 )                                     ) then
   !   
   !end if
        
   
         ! WaveStMod - Model switch for stretching incident wave kinematics to instantaneous free surface.

         ! TODO: We are only implementing WaveStMod = 0 (No stretching) at this point in time. 1 Mar 2013 GJH

   IF ( InitInp%Waves%WaveStMod /= 0 ) THEN
      CALL SetErrStat( ErrID_Fatal,'WaveStMod must be 0. Future versions of HydroDyn will once again support other wave stretching models.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF

   IF ( InitInp%Waves%WaveMod /= 6 .AND. InitInp%Morison%NMembers > 0 .AND. InitInp%Waves%WaveMod > 0 ) THEN
      
      IF ( ( InitInp%Waves%WaveStMod /= 0 ) .AND. ( InitInp%Waves%WaveStMod /= 1 ) .AND. &
            ( InitInp%Waves%WaveStMod /= 2 ) ) THEN ! (TODO: future version will support 3) .AND. ( InitInp%Waves%WaveStMod /= 3 ) )  THEN
         ErrMsg  = ' WaveStMod must be 0, 1, or 2.' !, or 3.'
         ErrStat = ErrID_Fatal
   
         RETURN
      END IF
   
      !IF ( ( InitInp%Waves%WaveStMod /= 3 ) .AND. ( InitInp%Waves%WaveMod == 5 ) )  THEN
      !   ErrMsg  = ' WaveStMod must be set to 3 when WaveMod is set to 5.'
      !   ErrStat = ErrID_Fatal
      !
      !   RETURN
      !END IF
      
         
   
   ELSE !don't use this one
   
         ! NOTE: Do not read in WaveStMod for floating platforms since it is
         !       inconsistent to use stretching (which is a nonlinear correction) for
         !       the viscous drag term in Morison's equation while not accounting for
         !       stretching in the diffraction and radiation problems (according to
         !       Paul Sclavounos, there are such corrections).  Instead, the viscous
         !       drag term from Morison's equation is computed by integrating up to
         !       the MSL, regardless of the instantaneous free surface elevation.
   
      InitInp%Waves%WaveStMod = 0
   
   END IF


      ! WaveTMax - Analysis time for incident wave calculations.

   IF ( InitInp%Waves%WaveMod == 0 )  THEN   ! .TRUE if we have incident waves.
      
      ! TODO: Issue warning if WaveTMax was not already 0.0 in this case.
      IF ( .NOT. EqualRealNos(InitInp%Waves%WaveTMax, 0.0_DbKi) ) THEN
         CALL WrScr( '  Setting WaveTMax to 0.0 since WaveMod = 0' )
      InitInp%Waves%WaveTMax = 0.0
      END IF
   ELSEIF ( InitInp%Waves%WaveMod == 5 ) THEN   ! User wave elevation file reading in
      IF (InitInp%TMax > InitInp%Waves%WaveTMax ) THEN
         CALL SetErrstat( ErrID_Fatal, '  WaveTMax must be larger than the simulation time for user wave elevations (WaveMod == 5).',ErrStat,ErrMsg,RoutineName)
         RETURN
      END IF
   ELSE
      IF (InitInp%TMax > InitInp%Waves%WaveTMax ) THEN
         CALL WrScr( '  WaveTMax is less then the simulation time.  Wave data will repeat every WaveTMax seconds.')
      END IF
   END IF   


      ! WaveDT - Time step for incident wave calculations

   IF ( InitInp%Waves%WaveMod > 0 )  THEN   ! .TRUE if we have incident waves.

      IF ( InitInp%Waves%WaveDT <= 0.0 )  THEN
         CALL SetErrStat( ErrID_Fatal,'WaveDT must be greater than zero.',ErrStat,ErrMsg,RoutineName)
         RETURN
      END IF
      
      IF ( (InitInp%Waves%WaveMod == 6) .AND. (.NOT. EqualRealNos(InitInp%Waves%WaveDT, InitInp%DT)) ) THEN
         CALL SetErrStat( ErrID_Fatal,'WaveDT must equal the simulation DT value when WaveMod = 6.',ErrStat,ErrMsg,RoutineName)
         RETURN
      END IF
   ELSE

      InitInp%Waves%WaveDT = 0.0

   END IF


       ! WaveHs - Significant wave height

   IF ( ( InitInp%Waves%WaveMod /= 0 ) .AND. ( InitInp%Waves%WaveMod /= 4 ) .AND. ( InitInp%Waves%WaveMod /= 5 ) ) THEN   ! .TRUE. (when WaveMod = 1, 2, 3, or 10) if we have plane progressive (regular), JONSWAP/Pierson-Moskowitz spectrum (irregular) waves, or white-noise waves, but not user-defined or GH Bladed wave data.

      IF ( InitInp%Waves%WaveHs <= 0.0 )  THEN
         CALL SetErrStat( ErrID_Fatal,'WaveHs must be greater than zero.',ErrStat,ErrMsg,RoutineName)
         RETURN
      END IF

   ELSE

      InitInp%Waves%WaveHs = 0.0

   END IF


      ! WaveTp - Peak spectral period.
   ! We commented out the if else block due to a bug when WaveMod == 3, and then WaveTp is hence set to 0.0.  See line 1092 of Waves.f90 (as of 11/24/2014) GJH
   !IF ( ( InitInp%Waves%WaveMod == 1 ) .OR. ( InitInp%Waves%WaveMod == 2 ) .OR. ( InitInp%Waves%WaveMod == 10 ) ) THEN   ! .TRUE. (when WaveMod = 1, 2, or 10) if we have plane progressive (regular), JONSWAP/Pierson-Moskowitz spectrum (irregular) waves.

      IF ( InitInp%Waves%WaveTp <= 0.0 )  THEN
         CALL SetErrStat( ErrID_Fatal,'WaveTp must be greater than zero.',ErrStat,ErrMsg,RoutineName)
         RETURN
      END IF

  ! ELSE

  !    InitInp%Waves%WaveTp = 0.0

  ! END IF


       ! WavePkShp - Peak shape parameter.

   CALL Conv2UC( InitInp%Waves%WavePkShpChr )    ! Convert Line to upper case.

   IF ( InitInp%Waves%WaveMod == 2 ) THEN   ! .TRUE if we have JONSWAP/Pierson-Moskowitz spectrum (irregular) waves, but not GH Bladed wave data.

      IF ( TRIM(InitInp%Waves%WavePkShpChr) == 'DEFAULT' )  THEN   ! .TRUE. when one wants to use the default value of the peak shape parameter, conditioned on significant wave height and peak spectral period.

         InitInp%Waves%WavePkShp = WavePkShpDefault ( InitInp%Waves%WaveHs, InitInp%Waves%WaveTp )

      ELSE                                   ! The input must have been specified numerically.

         READ (InitInp%Waves%WavePkShpChr,*,IOSTAT=IOS)  InitInp%Waves%WavePkShp
            CALL CheckIOS ( IOS, "", 'WavePkShp', NumType, ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2, ErrMsg2,ErrStat,ErrMsg,RoutineName)
            IF ( ErrStat >= AbortErrLev ) RETURN

         IF ( ( InitInp%Waves%WavePkShp < 1.0 ) .OR. ( InitInp%Waves%WavePkShp > 7.0 ) )  THEN
            CALL SetErrStat( ErrID_Fatal,'WavePkShp must be greater than or equal to 1 and less than or equal to 7.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF

      END IF

   ELSE

      InitInp%Waves%WavePkShp = 1.0

   END IF


      ! WvLowCOff and WvHiCOff - Wave Cut-off frequency
    
   IF ( InitInp%Waves%WvLowCOff < 0 ) THEN
      CALL SetErrStat( ErrID_Fatal,'WvLowCOff must be greater than or equal to zero.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF
   
      ! Threshold upper cut-off based on sampling rate
   IF ( EqualRealNos(InitInp%Waves%WaveDT, 0.0_DbKi) ) THEN
      InitInp%Waves%WvHiCOff = 10000.0;  ! This is not going to be used because WaveDT is zero.
   ELSE
      InitInp%Waves%WvHiCOff =  MIN( REAL( Pi/InitInp%Waves%WaveDT,SiKi), InitInp%Waves%WvHiCOff ) 
   END IF
   
   !TODO Issue warning if we changed WvHiCOff  GJH 7/24/13
   
   IF ( InitInp%Waves%WvLowCOff >= InitInp%Waves%WvHiCOff ) THEN
      CALL SetErrSTat( ErrID_Fatal,'WvLowCOff must be less than WvHiCOff.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF
   
   
        ! Copy over the first order frequency limits to the WAMIT2 module which needs them.
   InitInp%WAMIT2%WvLowCOff  = InitInp%Waves%WvLowCOff
   InitInp%WAMIT2%WvHiCOff   = InitInp%Waves%WvHiCOff


      ! WaveDir - Wave heading direction.

   IF ( ( InitInp%Waves%WaveMod > 0 ) .AND. ( InitInp%Waves%WaveMod /= 6 ) )  THEN   ! .TRUE if we have incident waves, but not user input wave data.

      IF ( ( InitInp%Waves%WaveDir <= -180.0 ) .OR. ( InitInp%Waves%WaveDir > 180.0 ) )  THEN
         CALL SetErrStat( ErrID_Fatal,'WaveDir must be greater than -180 and less than or equal to 180.',ErrStat,ErrMsg,RoutineName)
         RETURN
      END IF

   ELSE

      InitInp%Waves%WaveDir = 0.0

   END IF


      ! Multi-directional waves

      ! Check the WaveDirMod value
   IF ( InitInp%Waves%WaveDirMod < 0 .OR. InitInp%Waves%WaveDirMod > 1 ) THEN
      CALL SetErrStat( ErrID_Fatal,'WaveDirMod must be either 0 (No spreading) or 1 (COS2S spreading function)',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF

      ! Check if we are doing multidirectional waves or not.
      ! We can only use multi directional waves on WaveMod=2,3,4
   InitInp%Waves%WaveMultiDir = .FALSE.         ! Set flag to false to start
   IF ( InitInp%Waves%WaveMod >= 2 .AND. InitInp%Waves%WaveMod <= 4 .AND. InitInp%Waves%WaveDirMod == 1 ) THEN
      InitInp%Waves%WaveMultiDir = .TRUE.
   ELSEIF ( (InitInp%Waves%WaveMod < 2 .OR. InitInp%Waves%WaveMod >4) .AND. InitInp%Waves%WaveDirMod == 1 ) THEN
      CALL SetErrStat( ErrID_Warn,'WaveDirMod unused unless WaveMod == 2, 3, or 4.  Ignoring WaveDirMod.',ErrStat,ErrMsg,RoutineName)
   ENDIF


      !  Check to see if the for some reason the wave direction spreading range is set to zero.  If it is, 
      !  we don't have any spreading, so we will turn off the multidirectional waves.
   IF ( InitInp%Waves%WaveMultiDir .AND. EqualRealNos( InitInp%Waves%WaveDirRange, 0.0_SiKi ) ) THEN
      CALL SetErrStat( ErrID_Warn,' WaveDirRange set to zero, so multidirectional waves are turned off.',ErrStat,ErrMsg,RoutineName)
      InitInp%Waves%WaveMultiDir = .FALSE.
   ENDIF



      ! We check the following only if we set WaveMultiDir to true, otherwise ignore them and set them to zero
   IF ( InitInp%Waves%WaveMultiDir ) THEN

         ! Check WaveDirSpread
      IF ( InitInp%Waves%WaveDirSpread <= 0.0 ) THEN

         CALL SetErrStat( ErrID_Fatal,'WaveDirSpread cannot negative or zero.',ErrStat,ErrMsg,RoutineName)
         RETURN

      ENDIF


         ! Check that the number of wave directions is a positive odd number.
         !     -> If it is less than 0, error out.
         !     -> If it is even, we will increment it by 1.
      IF ( InitInp%Waves%WaveNDir <= 0_IntKi ) THEN
         CALL SetErrStat( ErrID_Fatal,' WaveNDir must be an odd number greater than 0.',ErrStat,ErrMsg,RoutineName)
         RETURN
      ENDIF

         ! Check that the value for WaveNDir is odd
      IF ( MODULO( InitInp%Waves%WaveNDir, 2_IntKi) == 0_IntKi ) THEN
         InitInp%Waves%WaveNDir  = InitInp%Waves%WaveNDir + 1
         CALL SetErrStat( ErrID_Warn,'WaveNDir must be odd.  Changing the value to '//Num2LStr(InitInp%Waves%WaveNDir),ErrStat,ErrMsg,RoutineName)
      ENDIF

         ! Now check that the WaveDirRange is less than 360 degrees (not sure why we would want that)
      IF ( InitInp%Waves%WaveDirRange > 360.0_ReKi ) THEN
         CALL SetErrStat( ErrID_Fatal,' WaveDirRange should be less than a full circle.',ErrStat,ErrMsg,RoutineName)
      ENDIF

   ELSE  ! Set everything to zero if we aren't going to use it

      InitInp%Waves%WaveNDir        = 1         ! Only one direction set -- this shouldn't get used later anyhow
      InitInp%Waves%WaveDirRange    = PiBy2     ! This is so that the constant C=1 in the COS2S function (it shouldn't get called, but in case it does)
      InitInp%Waves%WaveDirSpread   = 0.0

   END IF


       ! WaveSeed(1), !WaveSeed(2)

   IF ( .NOT. ( ( InitInp%Waves%WaveMod > 0 ) .AND. ( InitInp%Waves%WaveMod /= 5 ) .AND. ( InitInp%Waves%WaveMod /= 10 ) ) ) THEN   !.TRUE. for plane progressive (regular) with random phase or irregular wave 

      DO I = 1,2

         InitInp%Waves%WaveSeed(I) = 0

      END DO !I

   END IF


      ! WvKinFile

   IF ( InitInp%Waves%WaveMod == 5 .OR. InitInp%Waves%WaveMod == 6 ) THEN      ! .TRUE if we are to read user-supplied wave elevation or wave kinematics file(s).

      IF ( LEN_TRIM( InitInp%Waves%WvKinFile ) == 0 )  THEN
         CALL SetErrStat( ErrID_Fatal,'WvKinFile must not be an empty string.',ErrStat,ErrMsg,RoutineName)
         RETURN
      END IF

      IF ( PathIsRelative( InitInp%Waves%WvKinFile ) ) THEN
         CALL GetPath( TRIM(InitInp%InputFile), TmpPath )
         InitInp%Waves%WvKinFile    = TRIM(TmpPath)//TRIM(InitInp%Waves%WvKinFile)
      END IF
      InitInp%Waves%WriteWvKin = .FALSE.
   ELSE !don't use this one
      
#ifdef WRITE_WV_KIN
      IF ( LEN_TRIM( InitInp%Waves%WvKinFile ) == 0 )  THEN
         InitInp%Waves%WriteWvKin = .FALSE.
      ELSE
         InitInp%Waves%WriteWvKin = .TRUE.
         IF ( PathIsRelative( InitInp%Waves%WvKinFile ) ) THEN
            CALL GetPath( TRIM(InitInp%InputFile), TmpPath )
            InitInp%Waves%WvKinFile    = TRIM(TmpPath)//TRIM(InitInp%Waves%WvKinFile)
         END IF
      END IF
      
#else
      InitInp%Waves%WvKinFile = ""
      InitInp%Waves%WriteWvKin = .FALSE.
#endif
   END IF


      ! NWaveElev

   IF ( InitInp%Waves%NWaveElev < 0 ) THEN

      CALL SetErrStat( ErrID_Fatal,'NWaveElev must not be negative.',ErrStat,ErrMsg,RoutineName)
      RETURN

   END IF



      !-------------------------------------------------------------------------
      ! Check 2nd Order Waves section
      !-------------------------------------------------------------------------


      ! Difference frequency cutoffs

      ! WvLowCOffD and WvHiCOffD - Wave Cut-off frequency
   IF ( InitInp%Waves2%WvLowCOffD < 0 ) THEN
      CALL SetErrStat( ErrID_Fatal,'WvLowCOffD must be greater than or equal to zero.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF

      ! Check that the order given makes sense. 
   IF ( InitInp%Waves2%WvLowCOffD >= InitInp%Waves2%WvHiCOffD ) THEN
      CALL SetErrStat( ErrID_Fatal,'WvLowCOffD must be less than WvHiCOffD.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF
   

      ! Sum frequency cutoffs

      ! WvLowCOffS and WvHiCOffD - Wave Cut-off frequency
   IF ( InitInp%Waves2%WvLowCOffS < 0 ) THEN
      CALL SetErrStat( ErrID_Fatal,'WvLowCOffS must be greater than or equal to zero.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF

      ! Check that the order given makes sense. 
   IF ( InitInp%Waves2%WvLowCOffS >= InitInp%Waves2%WvHiCOffS ) THEN
      CALL SetErrStat( ErrID_Fatal,'WvLowCOffS must be less than WvHiCOffS.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF


        ! Copy over the 2nd order limits to the WAMIT2 module which needs them.
   InitInp%WAMIT2%WvLowCOffD  = InitInp%Waves2%WvLowCOffD
   InitInp%WAMIT2%WvHiCOffD   = InitInp%Waves2%WvHiCOffD
   InitInp%WAMIT2%WvLowCOffS  = InitInp%Waves2%WvLowCOffS
   InitInp%WAMIT2%WvHiCOffS   = InitInp%Waves2%WvHiCOffS



      !-------------------------------------------------------------------------
      ! Check Current section
      !-------------------------------------------------------------------------
      

      ! CurrMod - Current profile model switch

   IF ( ( InitInp%Current%CurrMod /= 0 ) .AND. ( InitInp%Current%CurrMod /= 1 ) .AND. ( InitInp%Current%CurrMod /= 2 ) )  THEN
      CALL SetErrStat( ErrID_Fatal,'CurrMod must be 0, 1, or 2.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF

   IF ( ( InitInp%Current%CurrMod /= 0 ) .AND. ( InitInp%Waves%WaveMod == 6 ) )  THEN
      CALL SetErrStat( ErrID_Fatal,'CurrMod must be set to 0 when WaveMod is set to 6: user-input wave data.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF


      ! CurrSSV0 - Sub-surface current velocity at still water level

   IF ( InitInp%Current%CurrMod == 1 )  THEN  ! .TRUE if we have standard current.

      IF ( InitInp%Current%CurrSSV0 < 0.0 )  THEN
         CALL SetErrStat( ErrID_Fatal,'CurrSSV0 must not be less than zero.',ErrStat,ErrMsg,RoutineName)
         RETURN
      END IF

   ELSE

      InitInp%Current%CurrSSV0 = 0.0

   END IF


      ! CurrSSDirChr - Sub-surface current heading direction

   IF ( InitInp%Current%CurrMod == 1 )  THEN  ! .TRUE if we have standard current.


      IF ( TRIM(InitInp%Current%CurrSSDirChr) == 'DEFAULT' )  THEN   ! .TRUE. when one wants to use the default value of codirectionality between sub-surface current and incident wave propogation heading directions.

         IF ( InitInp%Waves%WaveMod == 0 ) THEN
            CALL SetErrStat( ErrID_Fatal,'CurrSSDir must not be set to ''DEFAULT'' when WaveMod is set to 0.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF

         InitInp%Current%CurrSSDir = InitInp%Waves%WaveDir

      ELSE                                   ! The input must have been specified numerically.

         READ (InitInp%Current%CurrSSDirChr,*,IOSTAT=IOS)  InitInp%Current%CurrSSDir
            CALL CheckIOS ( IOS, "", 'CurrSSDir', NumType, ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2, ErrMsg2,ErrStat,ErrMsg,RoutineName)
            IF ( ErrStat >= AbortErrLev ) RETURN

         IF ( ( InitInp%Current%CurrSSDir <= -180.0 ) .OR. ( InitInp%Current%CurrSSDir > 180.0 ) )  THEN
            CALL SetErrStat( ErrID_Fatal,'CurrSSDir must be greater than -180 and less than or equal to 180.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF

      END IF


   ELSE

      InitInp%Current%CurrSSDir = 0.0

   END IF


      ! CurrNSRef - Near-surface current reference depth.

   IF ( InitInp%Current%CurrMod == 1 )  THEN  ! .TRUE if we have standard current.

      IF ( InitInp%Current%CurrNSRef <= 0.0 ) THEN
         CALL SetErrStat( ErrID_Fatal,'CurrNSRef must be greater than zero.',ErrStat,ErrMsg,RoutineName)
         RETURN
      END IF

   ELSE

      InitInp%Current%CurrNSRef = 0.0

   END IF



        ! CurrNSV0 - Near-surface current velocity at still water level.

   IF ( InitInp%Current%CurrMod == 1 )  THEN  ! .TRUE if we have standard current.

      IF ( InitInp%Current%CurrNSV0 < 0.0 ) THEN
         CALL SetErrStat( ErrID_Fatal,'CurrNSV0 must not be less than zero.',ErrStat,ErrMsg,RoutineName)
         RETURN
      END IF

   ELSE

      InitInp%Current%CurrNSV0 = 0.0

   END IF


      ! CurrNSDir - Near-surface current heading direction.

   IF ( InitInp%Current%CurrMod == 1 )  THEN  ! .TRUE if we have standard current.

      IF ( ( InitInp%Current%CurrNSDir <= -180.0 ) .OR. ( InitInp%Current%CurrNSDir > 180.0 ) )  THEN
         CALL SetErrStat( ErrID_Fatal,'CurrNSDir must be greater than -180 and less than or equal to 180.',ErrStat,ErrMsg,RoutineName)
         RETURN
      END IF

   ELSE

      InitInp%Current%CurrNSDir = 0.0

   END IF


      ! CurrDIV - Depth-independent current velocity.

   IF ( InitInp%Current%CurrMod == 1 )  THEN  ! .TRUE if we have standard current.

      IF ( InitInp%Current%CurrDIV < 0.0 ) THEN
         CALL SetErrStat( ErrID_Fatal,'CurrDIV must not be less than zero.',ErrStat,ErrMsg,RoutineName)
         RETURN
      END IF

   ELSE

      InitInp%Current%CurrDIV = 0.0

   END IF


      ! CurrDIDir - Depth-independent current heading direction.

   IF ( InitInp%Current%CurrMod == 1 )  THEN  ! .TRUE if we have standard current.

      IF ( ( InitInp%Current%CurrDIDir <= -180.0 ) .OR. ( InitInp%Current%CurrDIDir > 180.0 ) ) THEN
         CALL SetErrStat( ErrID_Fatal,'CurrDIDir must be greater than -180 and less than or equal to 180.',ErrStat,ErrMsg,RoutineName)
         RETURN
      END IF

   ELSE

      InitInp%Current%CurrDIDir = 0.0

   END IF

       ! PotFile - Root name of potential flow files

   IF ( InitInp%PotMod > 0 ) THEN
       IF ( LEN_TRIM( InitInp%PotFile ) == 0 ) THEN
         CALL SetErrStat( ErrID_Fatal,'PotFile must not be an empty string.',ErrStat,ErrMsg,RoutineName)
         RETURN
      END IF

         ! if this is a relative path, let's make it relative to the location of the main input file
         ! tell the WAMIT and WAMIT2 modules what the filename is

      IF ( PathIsRelative( InitInp%PotFile ) ) THEN
         CALL GetPath( TRIM(InitInp%InputFile), TmpPath )
         InitInp%PotFile            = TRIM(TmpPath)//TRIM(InitInp%PotFile)
      END IF
      InitInp%WAMIT%WAMITFile    = InitInp%PotFile
      InitInp%WAMIT2%WAMITFile   = InitInp%PotFile
      
         ! Set the flag for multidirectional waves for WAMIT2 module.  It needs to know since the Newman approximation
         ! can only use uni-directional waves.
      InitInp%WAMIT2%WaveMultiDir = InitInp%Waves%WaveMultiDir

   ELSE
      InitInp%PotFile            = ""
      InitInp%WAMIT%WAMITFile    = ""
      InitInp%WAMIT2%WAMITFile   = ""     
   END IF

      ! Set the WAMIT file name on the Convolution module
   InitInp%WAMIT%Conv_Rdtn%WAMITFile = InitInp%WAMIT%WAMITFile

      ! WAMITULEN - WAMIT characteristic body length scale

   IF ( InitInp%PotMod == 1 ) THEN

      InitInp%WAMIT2%WAMITULEN = InitInp%WAMIT%WAMITULEN    ! Copy to the WAMIT2 module info
      IF ( InitInp%WAMIT%WAMITULEN < 0.0 ) THEN
         CALL SetErrStat( ErrID_Fatal,'WAMITULEN must be positive.',ErrStat,ErrMsg,RoutineName)
         RETURN
      END IF
   ELSE

      InitInp%WAMIT%WAMITULEN = 1.0
      InitInp%WAMIT2%WAMITULEN = 1.0

   END IF


      ! PtfmVol0 - Displaced volume of water when the platform is in its undisplaced position

   IF ( InitInp%PotMod == 1 ) THEN

      IF ( InitInp%WAMIT%PtfmVol0 < 0.0 ) THEN
         CALL SetErrStat( ErrID_Fatal,'PtfmVol0 must not be negative.',ErrStat,ErrMsg,RoutineName)
         RETURN
      END IF

   ELSE

      InitInp%WAMIT%PtfmVol0 = 0.0

   END IF


      ! RdtnTMax - Analysis time for wave radiation kernel calculations
      ! NOTE: Use RdtnTMax = 0.0 to eliminate wave radiation damping

   IF ( InitInp%PotMod == 1 ) THEN

      IF ( InitInp%WAMIT%RdtnTMax < 0.0 ) THEN
         CALL SetErrStat( ErrID_Fatal,'RdtnTMax must not be negative.',ErrStat,ErrMsg,RoutineName)
         RETURN
      END IF

   ELSE

      InitInp%WAMIT%RdtnTMax = 0.0

   END IF

       ! RdtnDT - Time step for wave radiation kernel calculations

   IF ( InitInp%PotMod == 1 ) THEN

      CALL Conv2UC( InitInp%WAMIT%Conv_Rdtn%RdtnDTChr )    ! Convert Line to upper case.
      
      IF ( TRIM(InitInp%WAMIT%Conv_Rdtn%RdtnDTChr) == 'DEFAULT' )  THEN   ! .TRUE. when one wants to use the default value timestep provided by the glue code.

         InitInp%WAMIT%Conv_Rdtn%RdtnDT = InitInp%DT

      ELSE                                   ! The input must have been specified numerically.

         READ (InitInp%WAMIT%Conv_Rdtn%RdtnDTChr,*,IOSTAT=IOS)  InitInp%WAMIT%Conv_Rdtn%RdtnDT
            CALL CheckIOS ( IOS, "", 'RdtnDT', NumType, ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2, ErrMsg2,ErrStat,ErrMsg,RoutineName)
            IF ( ErrStat >= AbortErrLev ) RETURN
         
      END IF
      
      IF ( InitInp%WAMIT%Conv_Rdtn%RdtnDT <= 0.0 ) THEN
         CALL SetErrStat( ErrID_Fatal,'RdtnDT must be greater than zero.',ErrStat,ErrMsg,RoutineName)
         RETURN
      END IF

      if ( (.not. ( EqualRealNos(InitInp%DT, InitInp%WAMIT%Conv_Rdtn%RdtnDT) ) ) .and. ( (InitInp%WAMIT%ExctnMod > 1) .or. (InitInp%WAMIT%RdtnMod > 0) ) ) then
         call SetErrStat( ErrID_Fatal,'RdtnDT must be equal to the glue-code DT if PotMod = 1 and using RdtnMod > 0 or ExctnMod > 1.',ErrStat,ErrMsg,RoutineName)
         return
      end if
      
   ELSE

      InitInp%WAMIT%Conv_Rdtn%RdtnDT = 0.0

   END IF

   !-------------------------------------------------------------------------------------------------
   ! Data section for Floating platform force flags
   !-------------------------------------------------------------------------------------------------

!FIXME: ADP -- the error handling in this section is broken.

   ! If DEFAULT was requested, then the required value has already been set by the calling program
   IF ( TRIM(InitInp%PtfmSgFChr) /= 'DEFAULT' )  THEN

      READ (InitInp%PtfmSgFChr,*,IOSTAT=IOS)  InitInp%PtfmSgF
         CALL CheckIOS ( IOS, "", 'PtfmSgF', NumType, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2, ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) RETURN

   END IF

   IF ( TRIM(InitInp%PtfmSwFChr) /= 'DEFAULT' )  THEN

      READ (InitInp%PtfmSwFChr,*,IOSTAT=IOS)  InitInp%PtfmSwF
         CALL CheckIOS ( IOS, "", 'PtfmSwF', NumType, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2, ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) RETURN

   END IF

   IF ( TRIM(InitInp%PtfmHvFChr) /= 'DEFAULT' )  THEN

      READ (InitInp%PtfmHvFChr,*,IOSTAT=IOS)  InitInp%PtfmHvF
         CALL CheckIOS ( IOS, "", 'PtfmHvF', NumType, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2, ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) RETURN

   END IF

   IF ( TRIM(InitInp%PtfmRFChr) /= 'DEFAULT' )  THEN

      READ (InitInp%PtfmRFChr,*,IOSTAT=IOS)  InitInp%PtfmRF
         CALL CheckIOS ( IOS, "", 'PtfmRF', NumType, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2, ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) RETURN

   END IF

   IF ( TRIM(InitInp%PtfmPFChr) /= 'DEFAULT' )  THEN

      READ (InitInp%PtfmPFChr,*,IOSTAT=IOS)  InitInp%PtfmPF
         CALL CheckIOS ( IOS, "", 'PtfmPF', NumType, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2, ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) RETURN

   END IF

   IF ( TRIM(InitInp%PtfmYFChr) /= 'DEFAULT' )  THEN

      READ (InitInp%PtfmYFChr,*,IOSTAT=IOS)  InitInp%PtfmYF
         CALL CheckIOS ( IOS, "", 'PtfmYF', NumType, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2, ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) RETURN

   END IF


      ! Add checks that all platform DOF flags are true.  TODO:  Allow true or false once these have been implemented

   IF ( ( .NOT. InitInp%PtfmSgF ) .OR.  ( .NOT. InitInp%PtfmSwF ) .OR. ( .NOT. InitInp%PtfmHvF ) .OR. ( .NOT. InitInp%PtfmRF ) .OR. ( .NOT. InitInp%PtfmPF ) .OR. ( .NOT. InitInp%PtfmYF ) )THEN
!      CALL SetErrStat( ErrID_Fatal,'All platform DOF parameters must be set to TRUE.  Future versions of HydroDyn will support values of TRUE,  FALSE, or DEFAULT.',ErrStat,ErrMsg,RoutineName)
      CALL SetErrStat( ErrID_Warn,' Only the second-order floating platform force calculations (WAMIT2 sub-module) allow for selectively dissabling force DOF parameters, the first order (WAMIT sub-module) does not and will calculate all dimensions.  Future versions of HydroDyn will support values of TRUE,  FALSE, or DEFAULT for both modules.',ErrStat,ErrMsg,RoutineName)
   END IF



   !-------------------------------------------------------------------------------------------------
   ! Second order Forces Flags (WAMIT2 Module)
   !-------------------------------------------------------------------------------------------------
   !  We don't have separate inputs for the second order force component flags, rather they are taken
   !  from the platform section in the input file and copied into the InitInp%WAMIT2 derived type.
   !  Within the WAMIT2_Init subroutine, they are reset if necessary (some second order output files
   !  from WAMIT don't support all force components -- i.e. the *.8 files).

   InitInp%WAMIT2%PtfmSgF2    =  InitInp%PtfmSgF
   InitInp%WAMIT2%PtfmSwF2    =  InitInp%PtfmSwF
   InitInp%WAMIT2%PtfmHvF2    =  InitInp%PtfmHvF
   InitInp%WAMIT2%PtfmRF2     =  InitInp%PtfmRF
   InitInp%WAMIT2%PtfmPF2     =  InitInp%PtfmPF
   InitInp%WAMIT2%PtfmYF2     =  InitInp%PtfmYF



   !-------------------------------------------------------------------------------------------------
   ! Second order Forces due to Waves section (WAMIT2 Module)
   !-------------------------------------------------------------------------------------------------

   
      ! Check that we only specified one of MnDrift, NewmanApp, or DiffQTF
      !        (compared pairwise -- if any two are both true, we have a problem)
   IF ( ( InitInp%WAMIT2%MnDrift /= 0 .AND. InitInp%WAMIT2%NewmanApp /= 0 ) .OR. &
        ( InitInp%WAMIT2%DiffQTF /= 0 .AND. InitInp%WAMIT2%NewmanApp /= 0 ) .OR. &
        ( InitInp%WAMIT2%MnDrift /= 0 .AND. InitInp%WAMIT2%DiffQTF   /= 0 ) ) THEN
      CALL SetErrStat( ErrID_Fatal,'Only one of MnDrift, NewmanApp, or DiffQTF can be non-zero.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF


      ! Check MnDrift and set the flag indicating WAMIT2 should perform the mean drift calculation.
      ! Also make sure we have a valid input value for the file extension

   IF ( InitInp%WAMIT2%MnDrift == 0 ) THEN      ! not using MnDrift
      InitInp%WAMIT2%MnDriftF = .FALSE.
   ELSE IF ( InitInp%WAMIT2%MnDrift == 7  .OR. InitInp%WAMIT2%MnDrift == 8  .OR. InitInp%WAMIT2%MnDrift == 9 .OR. &
             InitInp%WAMIT2%MnDrift == 10 .OR. InitInp%WAMIT2%MnDrift == 11 .OR. InitInp%WAMIT2%MnDrift == 12 ) THEN   ! Valid values for MnDrift
      InitInp%WAMIT2%MnDriftF = .TRUE.
   ELSE     ! Must have received an invalid value
      CALL SetErrStat( ErrID_Fatal,'MnDrift can only have values of 0, 7, 8, 9, 10, 11, or 12.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF


      ! Check NewmanApp and set the flag indicating WAMIT2 should perform the mean drift calculation.
      ! Also make sure we have a valid input value for the file extension

   IF ( InitInp%WAMIT2%NewmanApp == 0 ) THEN    ! not using NewmanApp
      InitInp%WAMIT2%NewmanAppF = .FALSE.
   ELSE IF ( InitInp%WAMIT2%NewmanApp == 7  .OR. InitInp%WAMIT2%NewmanApp == 8  .OR. InitInp%WAMIT2%NewmanApp == 9 .OR. &
             InitInp%WAMIT2%NewmanApp == 10 .OR. InitInp%WAMIT2%NewmanApp == 11 .OR. InitInp%WAMIT2%NewmanApp == 12 ) THEN ! Valid values for NewmanApp
      InitInp%WAMIT2%NewmanAppF = .TRUE.
   ELSE     ! Must have received an invalid value
      CALL SetErrStat( ErrID_Fatal,'NewmanApp can only have values of 0, 7, 8, 9, 10, 11, or 12.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF


      ! Check DiffQTF and set the flag indicating WAMIT2 should perform the mean drift calculation.
      ! Also make sure we have a valid input value for the file extension

   IF ( InitInp%WAMIT2%DiffQTF == 0 ) THEN      ! not using DiffQTF method
       InitInp%WAMIT2%DiffQTFF = .FALSE.
   ELSE IF ( InitInp%WAMIT2%DiffQTF == 10 .OR. InitInp%WAMIT2%DiffQTF == 11 .OR. InitInp%WAMIT2%DiffQTF == 12 ) THEN    ! Valid values for DiffQTF
      InitInp%WAMIT2%DiffQTFF = .TRUE.
   ELSE
      CALL SetErrStat( ErrID_Fatal,'DiffQTF can only have values of 0, 10, 11, or 12.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF


      ! Check SumQTF and set the flag indicating WAMIT2 should perform the mean drift calculation.
      ! Also make sure we have a valid input value for the file extension

   IF ( InitInp%WAMIT2%SumQTF == 0 ) THEN       ! not using SumQTF method
      InitInp%WAMIT2%SumQTFF = .FALSE.
   ELSE IF ( InitInp%WAMIT2%SumQTF == 10 .OR. InitInp%WAMIT2%SumQTF == 11 .OR. InitInp%WAMIT2%SumQTF == 12 ) THEN       ! Valid values for SumQTF
      InitInp%WAMIT2%SumQTFF = .TRUE.
   ELSE
      CALL SetErrStat( ErrID_Fatal,'SumQTF can only have values of 0, 10, 11, or 12.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF


      ! Check that the min / max diff frequencies make sense if using any DiffQTF method
   IF ( InitInp%WAMIT2%DiffQTF /= 0 .OR. InitInp%WAMIT2%MnDrift /= 0 .OR. InitInp%WAMIT2%NewmanApp /=0 ) THEN
      IF ( ( InitInp%WAMIT2%WvHiCOffD < InitInp%WAMIT2%WvLowCOffD ) .OR. ( InitInp%WAMIT2%WvLowCOffD < 0.0 ) ) THEN
         CALL SetErrStat( ErrID_Fatal,'WvHiCOffD must be larger than WvLowCOffD. Both must be positive.',ErrStat,ErrMsg,RoutineName)
         RETURN
      END IF
   ELSE  ! set to zero since we don't need them
      InitInp%WAMIT2%WvLowCOffD  = 0.0
      InitInp%WAMIT2%WvHiCOffD  = 0.0
   END IF


      ! Check that the min / max diff frequencies make sense if using SumQTF
   IF ( InitInp%WAMIT2%SumQTF /= 0 ) THEN
      IF ( ( InitInp%WAMIT2%WvHiCOffS < InitInp%WAMIT2%WvLowCOffS ) .OR. ( InitInp%WAMIT2%WvLowCOffS < 0.0 ) ) THEN
         CALL SetErrStat( ErrID_Fatal,'WvHiCOffS must be larger than WvLowCOffS. Both must be positive.',ErrStat,ErrMsg,RoutineName)
         RETURN
      END IF
   ELSE  ! set to zero since we don't need them
      InitInp%WAMIT2%WvLowCOffS  = 0.0
      InitInp%WAMIT2%WvHiCOffS  = 0.0
   END IF


      ! now that it has been established that the input parameters for second order are good, we check to make sure that the WAMIT files actually exist.
      ! Check MnDrift file
   IF ( InitInp%WAMIT2%MnDrift /= 0) THEN
      ! Check if using QTF file types (10d, 11d, 12d) or not (7,8,9)
      IF ( InitInp%WAMIT2%MnDrift <= 9 ) THEN
         TmpExtension = TRIM(Num2LStr(InitInp%WAMIT2%MnDrift))
         INQUIRE( file=TRIM(InitInp%WAMIT2%WAMITFile)//'.'//TRIM(TmpExtension), exist=TmpFileExist )
      ELSE  ! 10, 11, 12
         TmpExtension = TRIM(Num2LStr(InitInp%WAMIT2%MnDrift))//'d'
         INQUIRE( file=TRIM(InitInp%WAMIT2%WAMITFile)//'.'//TRIM(TmpExtension), exist=TmpFileExist )
      ENDIF
      IF ( .not. TmpFileExist ) THEN
         CALL SetErrStat( ErrID_Fatal,'Cannot find the WAMIT file '//TRIM(InitInp%WAMIT2%WAMITFile)//'.'//TRIM(TmpExtension)// &
                    ' required by the MnDrift option.',ErrStat,ErrMsg,RoutineName)
         RETURN
      END IF
   END IF

      ! Check existence of NewmanApp file
   IF ( InitInp%WAMIT2%NewmanApp /= 0) THEN
      ! Check if using QTF file types (10d, 11d, 12d) or not (7,8,9)
      IF ( InitInp%WAMIT2%NewmanApp <= 9 ) THEN
         TmpExtension = TRIM(Num2LStr(InitInp%WAMIT2%NewmanApp))
         INQUIRE( file=TRIM(InitInp%WAMIT2%WAMITFile)//'.'//TRIM(TmpExtension), exist=TmpFileExist )
      ELSE  ! 10, 11, 12
         TmpExtension = TRIM(Num2LStr(InitInp%WAMIT2%NewmanApp))//'d'
         INQUIRE( file=TRIM(InitInp%WAMIT2%WAMITFile)//'.'//TRIM(TmpExtension), exist=TmpFileExist )
      ENDIF
      IF ( .not. TmpFileExist ) THEN
         CALL SetErrStat( ErrID_Fatal,'Cannot find the WAMIT file '//TRIM(InitInp%WAMIT2%WAMITFile)//'.'//TRIM(TmpExtension)// &
                    ' required by the NewmanApp option.',ErrStat,ErrMsg,RoutineName)
         RETURN
      END IF
   END IF

   IF ( InitInp%WAMIT2%DiffQTF /= 0) THEN
      TmpExtension = TRIM(Num2LStr(InitInp%WAMIT2%DiffQTF))//'d'
      INQUIRE( file=TRIM(InitInp%WAMIT2%WAMITFile)//'.'//TRIM(TmpExtension), exist=TmpFileExist )
      IF ( .not. TmpFileExist ) THEN
         CALL SetErrStat( ErrID_Fatal,'Cannot find the WAMIT file '//TRIM(InitInp%WAMIT2%WAMITFile)//'.'//TRIM(TmpExtension)// &
                    ' required by the DiffQTF option.',ErrStat,ErrMsg,RoutineName)
         RETURN
      END IF
   END IF

   IF ( InitInp%WAMIT2%SumQTF /= 0) THEN
      TmpExtension = TRIM(Num2LStr(InitInp%WAMIT2%SumQTF))//'s'
      INQUIRE( file=TRIM(InitInp%WAMIT2%WAMITFile)//'.'//TRIM(TmpExtension), exist=TmpFileExist )
      IF ( .not. TmpFileExist ) THEN
         CALL SetErrStat( ErrID_Fatal,'Cannot find the WAMIT file '//TRIM(InitInp%WAMIT2%WAMITFile)//'.'//TRIM(TmpExtension)// &
                    ' required by the SumQTF option.',ErrStat,ErrMsg,RoutineName)
         RETURN
      END IF
   END IF

   !..................
   ! check for ExctnMod = 2 requirements
   !..................
   if ( (InitInp%WAMIT%ExctnMod == 2) ) then

      if ( InitInp%Waves%WaveMod == 6 ) then
         call SetErrStat( ErrID_Fatal, 'Externally generated full wave-kinematics time series cannot be used with state-space wave excitations. Set WaveMod 0, 1, 1P#, 2, 3, 4, or 5.', ErrStat, ErrMsg, RoutineName )
      end if
      
      if ( InitInp%Waves%WaveDirMod /= 0 ) then
         call SetErrStat( ErrID_Fatal, 'Directional spreading cannot be used with state-space wave excitations. Set WaveDirMod=0.', ErrStat, ErrMsg, RoutineName )
      end if
      
      if ( InitInp%Waves2%WvDiffQTFF ) then
         call SetErrStat( ErrID_Fatal, 'Cannot use full difference-frequency 2nd-order wave kinematics with state-space wave excitations. Set WvDiffQTF=FALSE.', ErrStat, ErrMsg, RoutineName )
      end if
      
      if ( InitInp%Waves2%WvSumQTFF ) then
         call SetErrStat( ErrID_Fatal, 'Cannot use full summation-frequency 2nd-order wave kinematics with state-space wave excitations. Set WvSumQTF=FALSE.', ErrStat, ErrMsg, RoutineName )
      end if

      if ( InitInp%PotMod /= 1 ) then
         call SetErrStat( ErrID_Fatal, 'Potential-flow model via WAMIT must be used with state-space wave excitations. Set PotMod= 1.', ErrStat, ErrMsg, RoutineName )
      end if
      
      if ( InitInp%WAMIT2%MnDrift /= 0 ) then
         call SetErrStat( ErrID_Fatal, 'Mean-drift 2nd-order forces cannot be used with state-space wave excitations. Set MnDrift=0.', ErrStat, ErrMsg, RoutineName )
      end if

      if ( InitInp%WAMIT2%NewmanApp /= 0 ) then
         call SetErrStat( ErrID_Fatal, "Mean- and slow-drift 2nd-order forces computed with Newman's approximation cannot be used with state-space wave excitations. Set NewmanApp=0.", ErrStat, ErrMsg, RoutineName )
      end if
      
      if ( InitInp%WAMIT2%DiffQTF /= 0 ) then
         call SetErrStat( ErrID_Fatal, 'Full difference-frequency 2nd-order forces computed with full QTF cannot be used with state-space wave excitations. Set DiffQTF=0.', ErrStat, ErrMsg, RoutineName )
      end if

      if ( InitInp%WAMIT2%SumQTF /= 0 ) then
         call SetErrStat( ErrID_Fatal, 'Full summation-frequency 2nd-order forces computed with full QTF cannot be used with State-space wave excitations. Set SumQTF=0.', ErrStat, ErrMsg, RoutineName )
      end if

   end if

   !..................
   ! check for linearization
   !..................
   if (InitInp%Linearize) then
      
      if ( InitInp%Waves%WaveMod /= 0 ) then
         call SetErrStat( ErrID_Fatal, 'Still water conditions must be used for linearization. Set WaveMod=0.', ErrStat, ErrMsg, RoutineName )
      end if
      
      if ( InitInp%Waves%WaveDirMod /= 0 ) then
         call SetErrStat( ErrID_Fatal, 'No directional spreading must be used for linearization. Set WaveDirMod=0.', ErrStat, ErrMsg, RoutineName )
      end if
      
      if ( InitInp%Waves2%WvDiffQTFF ) then
         call SetErrStat( ErrID_Fatal, 'Cannot use full difference-frequency 2nd-order wave kinematics for linearization. Set WvDiffQTF=FALSE.', ErrStat, ErrMsg, RoutineName )
      end if
      
      if ( InitInp%Waves2%WvSumQTFF ) then
         call SetErrStat( ErrID_Fatal, 'Cannot use full summation-frequency 2nd-order wave kinematics for linearization. Set WvSumQTF=FALSE.', ErrStat, ErrMsg, RoutineName )
      end if

      if ( InitInp%PotMod > 1 ) then
         call SetErrStat( ErrID_Fatal, 'Potential-flow model cannot be set to FIT for linearization. Set PotMod= 0 or 1.', ErrStat, ErrMsg, RoutineName )
      end if
      
      if ( (InitInp%WAMIT%ExctnMod == 1) ) then
         call SetErrStat( ErrID_Fatal, 'Cannot set wave excitation model to DFT for linearization. Set ExctnMod=0 or 2.', ErrStat, ErrMsg, RoutineName )
      end if

      if ( InitInp%WAMIT%RdtnMod == 1 ) then
         call SetErrStat( ErrID_Fatal, 'Cannot set wave radiation model to convolution for linearization. Set RdtnMod=0 or 2.', ErrStat, ErrMsg, RoutineName )
      end if
      
      if ( InitInp%WAMIT2%MnDrift /= 0 ) then
         call SetErrStat( ErrID_Fatal, 'Mean-drift 2nd-order forces cannot be used for linearization. Set MnDrift=0.', ErrStat, ErrMsg, RoutineName )
      end if

      if ( InitInp%WAMIT2%NewmanApp /= 0 ) then
         call SetErrStat( ErrID_Fatal, "Mean- and slow-drift 2nd-order forces computed with Newman's approximation cannot be used for linearization. Set NewmanApp=0.", ErrStat, ErrMsg, RoutineName )
      end if
      
      if ( InitInp%WAMIT2%DiffQTF /= 0 ) then
         call SetErrStat( ErrID_Fatal, 'Full difference-frequency 2nd-order forces computed with full QTF cannot be used for linearization. Set DiffQTF=0.', ErrStat, ErrMsg, RoutineName )
      end if

      if ( InitInp%WAMIT2%SumQTF /= 0 ) then
         call SetErrStat( ErrID_Fatal, 'Full summation-frequency 2nd-order forces computed with full QTF cannot be used for linearization. Set SumQTF=0.', ErrStat, ErrMsg, RoutineName )
      end if

   end if
 




   !-------------------------------------------------------------------------------------------------
   ! Member Joints Section
   !-------------------------------------------------------------------------------------------------

   IF ( InitInp%Morison%NJoints < 0 ) THEN
      CALL SetErrStat( ErrID_Fatal,'NJoints parameter cannot be negative.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF

   IF ( InitInp%Morison%NJoints == 1 ) THEN
      CALL SetErrStat( ErrID_Fatal,'NJoints parameter cannot be set to 1.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF

     
     
      ! Check the axial coefs are >= 0 and IDs are unique
   IF ( InitInp%Morison%NAxCoefs > 0 ) THEN
   
      DO I = 1,InitInp%Morison%NAxCoefs 
         
         !IF (  .NOT. EqualRealNos(InitInp%Morison%AxialCoefs(I)%AxCd, 0.0) ) THEN
         !   ErrMsg  = ' AxCd must be equal to zero.  Future versions will allow for non-zero axial coefficients.'
         !   ErrStat = ErrID_Fatal
         !   RETURN
         !END IF   
         !IF (  .NOT. EqualRealNos(InitInp%Morison%AxialCoefs(I)%AxCa, 0.0) ) THEN
         !   ErrMsg  = ' AxCa must be equal to zero.  Future versions will allow for non-zero axial coefficients.'
         !   ErrStat = ErrID_Fatal
         !   RETURN
         !END IF   
         
         ! TODO: Once Axial Coefs are working remove the above checks and uncomment the checks below.  GJH 9/29/2013
         IF (  InitInp%Morison%AxialCoefs(I)%AxCd < 0 ) THEN
            CALL SetErrStat( ErrID_Fatal,'AxCd must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF   
         IF (  InitInp%Morison%AxialCoefs(I)%AxCa < 0 ) THEN
            CALL SetErrStat( ErrID_Fatal,'AxCa must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF   
         
            ! Make sure that the current AxCoefID is not used elsewhere in the table.
         DO J = I+1,InitInp%Morison%NAxCoefs
            IF ( InitInp%Morison%AxialCoefs(I)%AxCoefID == InitInp%Morison%AxialCoefs(J)%AxCoefID ) THEN
               CALL SetErrStat( ErrID_Fatal,'Duplicate AxCoefIDs were found in the Axial Coefficients table.',ErrStat,ErrMsg,RoutineName)
               RETURN
            END IF
         END DO
   
      END DO
      
   END IF


      ! Check JointOvrlp values
   InitInp%Morison%TotalPossibleSuperMembers = 0
   
   IF ( InitInp%Morison%NJoints > 1 ) THEN

      ! Initialize Joints
      DO I = 1,InitInp%Morison%NJoints
         InitInp%Morison%InpJoints(I)%NConnections   = 0
      END DO

      
      
      
      DO I = 1,InitInp%Morison%NJoints

            ! Make sure that the current JointID is not used elsewhere in the table.
         DO J = I+1,InitInp%Morison%NJoints
            IF ( InitInp%Morison%InpJoints(I)%JointID == InitInp%Morison%InpJoints(J)%JointID ) THEN
               CALL SetErrStat( ErrID_Fatal,'Duplicate JointIDs were found in the Member Joints table.',ErrStat,ErrMsg,RoutineName)
               RETURN
            END IF
         END DO

            ! Add up total number of joints flagged with JoinOvrlp = 1 option
         IF ( InitInp%Morison%InpJoints(I)%JointOvrlp == 1 ) THEN
            InitInp%Morison%TotalPossibleSuperMembers = InitInp%Morison%TotalPossibleSuperMembers + 1
         END IF

            ! Check that every joint id is used at least once in the members table
         JointUsed = .FALSE.
         DO J = 1, InitInp%Morison%NMembers
         
            IF ( InitInp%Morison%InpMembers(J)%MJointID1 == InitInp%Morison%InpJoints(I)%JointID ) THEN
               JointUsed = .TRUE.
               EXIT
            END IF
            IF ( InitInp%Morison%InpMembers(J)%MJointID2 == InitInp%Morison%InpJoints(I)%JointID ) THEN
               JointUsed = .TRUE.
               EXIT
            END IF
         END DO
         
         IF ( .NOT. JointUsed ) THEN
            CALL SetErrStat( ErrID_Fatal,'Every JointID in the Joints table must appear once in the Members table.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF  
   ! TODO : Implement Super member elements. GJH 7/24/13
   
         IF ( InitInp%Morison%InpJoints(I)%JointOvrlp /= 0  ) THEN
            CALL SetErrStat( ErrID_Fatal,'JointOvrlp parameter must be set to 0.  Future versions of HydroDyn will support vales of 0 or 1.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF
         !IF ( ( InitInp%Morison%InpJoints(I)%JointOvrlp < 0 ) .OR. ( InitInp%Morison%InpJoints(I)%JointOvrlp > 1 ) ) THEN
         !   ErrMsg  = ' JointOvrlp parameter must be set to 0 or 1.'
         !   ErrStat = ErrID_Fatal
         !   RETURN
         !END IF
         
            ! Make sure the axial coef id appears in the Ax table
         IF ( InitInp%Morison%NAxCoefs > 0 ) THEN
            InitInp%Morison%InpJoints(I)%JointAxIDIndx = -1
            DO J = 1,InitInp%Morison%NAxCoefs         
               IF ( InitInp%Morison%InpJoints(I)%JointAxID == InitInp%Morison%AxialCoefs(J)%AxCoefID ) &
                  InitInp%Morison%InpJoints(I)%JointAxIDIndx = J   
            END DO
            IF ( InitInp%Morison%InpJoints(I)%JointAxIDIndx == -1 ) THEN
               CALL SetErrStat( ErrID_Fatal,'The specified JointAxID in the Joints Table does not appear in the Axial Coefficients table.',ErrStat,ErrMsg,RoutineName)
               RETURN
            END IF
         ELSE
            ! TODO: Issue error because we need at least one Axial coef table entry GJH  8/1/31
         END IF
      
      END DO
   END IF


   !-------------------------------------------------------------------------------------------------
   ! Member Cross-section Properties Section
   !-------------------------------------------------------------------------------------------------

   IF ( InitInp%Morison%NPropSets < 0 ) THEN
      CALL SetErrStat( ErrID_Fatal,'Number of member cross-section property sets must be greater than zero.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF

   IF ( InitInp%Morison%NPropSets > 0 ) THEN

      DO I = 1,InitInp%Morison%NPropSets

            ! Make sure that the current JointID is not used elsewhere in the table.
         DO J = I+1,InitInp%Morison%NPropSets
            IF ( InitInp%Morison%MPropSets(I)%PropSetID == InitInp%Morison%MPropSets(J)%PropSetID ) THEN
               CALL SetErrStat( ErrID_Fatal,'Duplicate PropSetIDs were found in the Member Cross-section Properties table.',ErrStat,ErrMsg,RoutineName)
               RETURN
            END IF
         END DO

         IF ( ( InitInp%Morison%MPropSets(I)%PropD < 0 ) .OR.  ( InitInp%Morison%MPropSets(I)%PropThck < 0 ) .OR. ( ( InitInp%Morison%MPropSets(I)%PropD - InitInp%Morison%MPropSets(I)%PropThck / 2.0 ) < 0) ) THEN
            CALL SetErrStat( ErrID_Fatal,'PropD and PropThck must be greater than zero and (PropD - propThck/2 ) must be greater than zero.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF
      END DO

   END IF


   !-------------------------------------------------------------------------------------------------
   ! Simple hydrodynamic coefficients Section
   !-------------------------------------------------------------------------------------------------

   IF ( InitInp%Morison%SimplCd < 0 ) THEN
      CALL SetErrStat( ErrID_Fatal,'SimplCd must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF
   IF ( InitInp%Morison%SimplCdMG < 0 ) THEN
      CALL SetErrStat( ErrID_Fatal,'SimplCdMG must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF
   IF ( InitInp%Morison%SimplCa < 0 ) THEN
      CALL SetErrStat( ErrID_Fatal,'SimplCa must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF
   IF ( InitInp%Morison%SimplCaMG < 0 ) THEN
      CALL SetErrStat( ErrID_Fatal,'SimplCaMG must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF
   IF ( InitInp%Morison%SimplAxCa < 0 ) THEN
      CALL SetErrStat( ErrID_Fatal,'SimplCa must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF
   IF ( InitInp%Morison%SimplAxCaMG < 0 ) THEN
      CALL SetErrStat( ErrID_Fatal,'SimplCaMG must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF


   !-------------------------------------------------------------------------------------------------
   ! Depth-based Hydrodynamic Coefficients Section
   !-------------------------------------------------------------------------------------------------

   IF ( InitInp%Morison%NCoefDpth < 0 ) THEN
      CALL SetErrStat( ErrID_Fatal,'NCoefDpth must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF


   IF ( InitInp%Morison%NCoefDpth > 0 ) THEN
      MinDepth = 99999999.0
      MaxDepth = -99999999.0
      DO I = 1,InitInp%Morison%NCoefDpth

            ! Record the minimum and maximum depths covered by this table.  This will be used as part of a consistency check
            ! in the members table, below.
         IF (  InitInp%Morison%CoefDpths(I)%Dpth < MinDepth ) THEN
            MinDepth = InitInp%Morison%CoefDpths(I)%Dpth
         ELSE
            CALL SetErrStat( ErrID_Fatal,'The rows of the Depth-based Hydrodynamic Coefficients table must be ordered with increasing depth (decreasing Z).',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF
         IF ( InitInp%Morison%CoefDpths(I)%Dpth > MaxDepth ) THEN
            MaxDepth = InitInp%Morison%CoefDpths(I)%Dpth
         END IF

            ! Make sure that the current Dpth is not used elsewhere in the table.
         DO J = I+1,InitInp%Morison%NCoefDpth
            IF ( EqualRealNos( InitInp%Morison%CoefDpths(I)%Dpth, InitInp%Morison%CoefDpths(J)%Dpth ) ) THEN
               CALL SetErrStat( ErrID_Fatal,'Duplicate Dpths were found in the Depth-based Hydrodynamic Coefficients table.',ErrStat,ErrMsg,RoutineName)
               RETURN
            END IF
         END DO

         IF ( InitInp%Morison%CoefDpths(I)%DpthCd < 0 ) THEN
            CALL SetErrStat( ErrID_Fatal,'In the Depth-based hydrodynamic coefficients table, DpthCd must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF
         IF ( InitInp%Morison%CoefDpths(I)%DpthCdMG < 0 ) THEN
            CALL SetErrStat( ErrID_Fatal,'In the Depth-based hydrodynamic coefficients table, DpthCdMG must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF
         IF ( InitInp%Morison%CoefDpths(I)%DpthCa < 0 ) THEN
            CALL SetErrStat( ErrID_Fatal,'In the Depth-based hydrodynamic coefficients table, DpthCa must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF
         IF ( InitInp%Morison%CoefDpths(I)%DpthCaMG < 0 ) THEN
            CALL SetErrStat( ErrID_Fatal,'In the Depth-based hydrodynamic coefficients table, DpthCaMG must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF
         IF ( InitInp%Morison%CoefDpths(I)%DpthAxCa < 0 ) THEN 
            CALL SetErrStat( ErrID_Fatal,'In the Depth-based hydrodynamic coefficients table, DpthAxCa must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF
         IF ( InitInp%Morison%CoefDpths(I)%DpthAxCaMG < 0 ) THEN
            CALL SetErrStat( ErrID_Fatal,'In the Depth-based hydrodynamic coefficients table, DpthAxCaMG must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF
         IF ( InitInp%Morison%CoefDpths(I)%DpthAxCp < 0 ) THEN 
            CALL SetErrStat( ErrID_Fatal,'In the Depth-based hydrodynamic coefficients table, DpthAxCp must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF
         IF ( InitInp%Morison%CoefDpths(I)%DpthAxCpMG < 0 ) THEN
            CALL SetErrStat( ErrID_Fatal,'In the Depth-based hydrodynamic coefficients table, DpthAxCpMG must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF
      END DO

      ! TODO: Sort the table based on depth so that a linear interpolation can be easily performed between entries.

   END IF


   !-------------------------------------------------------------------------------------------------
   ! Member-based Hydrodynamic Coefficients Section
   !-------------------------------------------------------------------------------------------------

   IF ( InitInp%Morison%NCoefMembers < 0 ) THEN
      CALL SetErrStat( ErrID_Fatal,'NCoefMembers must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF

   IF ( InitInp%Morison%NCoefMembers > 0 ) THEN

      DO I = 1,InitInp%Morison%NCoefMembers

            ! Make sure that the current MemberID is not used elsewhere in the table.
         DO J = I+1,InitInp%Morison%NCoefMembers
            IF ( InitInp%Morison%CoefMembers(I)%MemberID == InitInp%Morison%CoefMembers(J)%MemberID ) THEN
               CALL SetErrStat( ErrID_Fatal,'Duplicate MemberIDs were found in the Member-based Hydrodynamic coefficients table.',ErrStat,ErrMsg,RoutineName)
               RETURN
            END IF
         END DO



         IF ( InitInp%Morison%CoefMembers(I)%MemberCd1 < 0 ) THEN
            CALL SetErrStat( ErrID_Fatal,'In the member-based hydrodynamic coefficients table, MemberCd1 must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF
         IF ( InitInp%Morison%CoefMembers(I)%MemberCd2 < 0 ) THEN
            CALL SetErrStat( ErrID_Fatal,'In the member-based hydrodynamic coefficients table, MemberCd2 must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF
         IF ( InitInp%Morison%CoefMembers(I)%MemberCdMG1 < 0 ) THEN
            CALL SetErrStat( ErrID_Fatal,'In the member-based hydrodynamic coefficients table, MemberCdMG1 must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF
         IF ( InitInp%Morison%CoefMembers(I)%MemberCdMG2 < 0 ) THEN
            CALL SetErrStat( ErrID_Fatal,'In the member-based hydrodynamic coefficients table, MemberCdMG2 must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF
         IF ( InitInp%Morison%CoefMembers(I)%MemberCa1 < 0 ) THEN
            CALL SetErrStat( ErrID_Fatal,'In the member-based hydrodynamic coefficients table, MemberCa1 must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF
         IF ( InitInp%Morison%CoefMembers(I)%MemberCa2 < 0 ) THEN
            CALL SetErrStat( ErrID_Fatal,'In the member-based hydrodynamic coefficients table, MemberCa2 must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF
         IF ( InitInp%Morison%CoefMembers(I)%MemberCaMG1 < 0 ) THEN
            CALL SetErrStat( ErrID_Fatal,'In the member-based hydrodynamic coefficients table, MemberCaMG1 must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF
         IF ( InitInp%Morison%CoefMembers(I)%MemberCaMG2 < 0 ) THEN
            CALL SetErrStat( ErrID_Fatal,'In the member-based hydrodynamic coefficients table, MemberCaMG2 must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF
         IF ( InitInp%Morison%CoefMembers(I)%MemberAxCa1 < 0 ) THEN 
            CALL SetErrStat( ErrID_Fatal,'In the member-based hydrodynamic coefficients table, MemberCa1 must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF
         IF ( InitInp%Morison%CoefMembers(I)%MemberAxCa2 < 0 ) THEN 
            CALL SetErrStat( ErrID_Fatal,'In the member-based hydrodynamic coefficients table, MemberCa2 must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF
         IF ( InitInp%Morison%CoefMembers(I)%MemberAxCaMG1 < 0 ) THEN 
            CALL SetErrStat( ErrID_Fatal,'In the member-based hydrodynamic coefficients table, MemberCaMG1 must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF
         IF ( InitInp%Morison%CoefMembers(I)%MemberAxCaMG2 < 0 ) THEN 
            CALL SetErrStat( ErrID_Fatal,'In the member-based hydrodynamic coefficients table, MemberCaMG2 must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF
      END DO

   END IF


   !-------------------------------------------------------------------------------------------------
   ! Members Section
   !-------------------------------------------------------------------------------------------------

   IF ( InitInp%Morison%NMembers < 0 ) THEN
      CALL SetErrStat( ErrID_Fatal,'NMembers in the Members table must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF

   IF ( InitInp%Morison%NMembers > 0 ) THEN

         ! Initialize all member data
      DO I = 1,InitInp%Morison%NMembers
         InitInp%Morison%InpMembers(I)%MJointID1Indx    = -1
         InitInp%Morison%InpMembers(I)%MJointID2Indx    = -1
         InitInp%Morison%InpMembers(I)%MPropSetID1Indx  = -1
         InitInp%Morison%InpMembers(I)%MPropSetID2Indx  = -1
         InitInp%Morison%InpMembers(I)%MmbrFilledIDIndx = -1
         InitInp%Morison%InpMembers(I)%MmbrCoefIDIndx   = -1
         InitInp%Morison%InpMembers(I)%NumSplits        = 0
         InitInp%Morison%InpMembers(I)%Splits           = 0.0_ReKi
      END DO

      DO I = 1,InitInp%Morison%NMembers

            ! Make sure that the current MemberID is not used elsewhere in the table.
         DO J = I+1,InitInp%Morison%NMembers
            IF ( InitInp%Morison%InpMembers(I)%MemberID == InitInp%Morison%InpMembers(J)%MemberID ) THEN
               CALL SetErrStat( ErrID_Fatal,'Duplicate MemberIDs were found in the Members table.',ErrStat,ErrMsg,RoutineName)
               RETURN
            END IF
         END DO

            ! Find JointID1 and JointID2 in the Joint table and then record their index locations in the Joint table
         DO J = 1,InitInp%Morison%NJoints
            IF ( InitInp%Morison%InpMembers(I)%MJointID1 == InitInp%Morison%InpJoints(J)%JointID ) THEN
               InitInp%Morison%InpMembers(I)%MJointID1Indx = J
               InitInp%Morison%InpJoints(J)%NConnections = InitInp%Morison%InpJoints(J)%NConnections + 1
               InitInp%Morison%InpJoints(J)%ConnectionList(InitInp%Morison%InpJoints(J)%NConnections) = I
            END IF
            IF ( InitInp%Morison%InpMembers(I)%MJointID2 == InitInp%Morison%InpJoints(J)%JointID ) THEN
               InitInp%Morison%InpMembers(I)%MJointID2Indx = J
               InitInp%Morison%InpJoints(J)%NConnections = InitInp%Morison%InpJoints(J)%NConnections + 1
               InitInp%Morison%InpJoints(J)%ConnectionList(InitInp%Morison%InpJoints(J)%NConnections) = I
            END IF
         END DO
         
            ! Make sure that a JointID entry in the Joints table was found
         IF ( InitInp%Morison%InpMembers(I)%MJointID1Indx == -1 ) THEN
            CALL SetErrStat( ErrID_Fatal,'JointID1 in the Members table does not appear in the Joints table.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF
         IF ( InitInp%Morison%InpMembers(I)%MJointID2Indx == -1 ) THEN
            CALL SetErrStat( ErrID_Fatal,'JointID2 in the Members table does not appear in the Joints table.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF

            ! Make sure we do not have any zero length members
         lvec = InitInp%Morison%InpJoints(InitInp%Morison%InpMembers(I)%MJointID1Indx)%JointPos - InitInp%Morison%InpJoints(InitInp%Morison%InpMembers(I)%MJointID2Indx)%JointPos
         l = sqrt( lvec(1)*lvec(1) + lvec(2)*lvec(2) + lvec(3)*lvec(3) )
         IF ( EqualRealNos(0.0_ReKi, l) ) THEN
            CALL SetErrStat( ErrID_Fatal,'A member cannot have zero length.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF

            ! Find MPropSetID1 and MPropSetID2 in the Member cross-section properties table and then record their index locations
         DO J = 1,InitInp%Morison%NPropSets



            IF ( InitInp%Morison%InpMembers(I)%MPropSetID1 == InitInp%Morison%MPropSets(J)%PropSetID ) THEN
               InitInp%Morison%InpMembers(I)%MPropSetID1Indx = J
            END IF
            IF ( InitInp%Morison%InpMembers(I)%MPropSetID2 == InitInp%Morison%MPropSets(J)%PropSetID ) THEN
               InitInp%Morison%InpMembers(I)%MPropSetID2Indx = J
            END IF
         END DO

            ! Make sure that a PropSetID entry in the Member cross-section properties table was found
         IF ( InitInp%Morison%InpMembers(I)%MPropSetID1Indx == -1 ) THEN
            CALL SetErrStat( ErrID_Fatal,'MPropSetID1 in the Members table does not appear in the Member cross-section properties table.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF
         IF ( InitInp%Morison%InpMembers(I)%MPropSetID2Indx == -1 ) THEN
            CALL SetErrStat( ErrID_Fatal,'MPropSetID2 in the Members table does not appear in the Member cross-section properties table.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF


         ! NOTE: We cannot test that MDivSize > MemberLength yet because there may be a joint overlap which is going to alter the final length of this member

         IF ( InitInp%Morison%InpMembers(I)%MDivSize <= 0 ) THEN
            CALL SetErrStat( ErrID_Fatal,'MDivSize must be greater than zero.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF


         IF ( ( InitInp%Morison%InpMembers(I)%MCoefMod /= 1 ) .AND. ( InitInp%Morison%InpMembers(I)%MCoefMod /= 2 ) .AND. ( InitInp%Morison%InpMembers(I)%MCoefMod /= 3 ) )  THEN
            CALL SetErrStat( ErrID_Fatal,'MCoefMod must be 1, 2, or 3.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF

         IF ( InitInp%Morison%InpMembers(I)%MCoefMod == 2 ) THEN
            IF ( InitInp%Morison%NCoefDpth == 0 ) THEN
               CALL SetErrStat( ErrID_Fatal,'NCoefDpth must be greater than zero when a member is using a depth-based coefficient model.',ErrStat,ErrMsg,RoutineName)
               RETURN
            END IF
               ! We will not extrapolate depth-based coefficient values, so make sure that the depth-based table has values that are outside the depth range of this member
               ! NOTE: This is actually potentially overly conservative because the final member may be shorter due to joint overlap handling.
            z1 = InitInp%Morison%InpJoints( InitInp%Morison%InpMembers(I)%MJointID1Indx )%JointPos(3)
            z2 = InitInp%Morison%InpJoints( InitInp%Morison%InpMembers(I)%MJointID2Indx )%JointPos(3)
            MinMembrDpth = min( z1, z2 )
            MaxMembrDpth = max( z1, z2 )
            IF ( ( MinMembrDpth < MinDepth ) .OR. ( MaxMembrDpth > MaxDepth ) ) THEN
               CALL SetErrStat( ErrID_Fatal,'This member uses a depth-based coefficient model, but the member depth is outside the range of values provided in the depth-based hydrodynamic coefficients table.',ErrStat,ErrMsg,RoutineName)
               RETURN
            END IF

         END IF


         IF ( InitInp%Morison%InpMembers(I)%MCoefMod == 3 ) THEN
            IF ( InitInp%Morison%NCoefMembers == 0 ) THEN
               CALL SetErrStat( ErrID_Fatal,'NCoefMembers must be greater than zero when a member is using a member-based coefficient model.',ErrStat,ErrMsg,RoutineName)
               RETURN
            END IF
               ! Make sure this id appears in the Members table and mark it's location for future use
            FoundID = .FALSE.
            DO J = 1,InitInp%Morison%NCoefMembers
               IF ( InitInp%Morison%CoefMembers(J)%MemberID == InitInp%Morison%InpMembers(I)%MemberID ) THEN
                  FoundID = .TRUE.
                  InitInp%Morison%InpMembers(I)%MmbrCoefIDIndx = J
               END IF
            END DO

            IF ( .NOT. FoundID ) THEN
               CALL SetErrStat( ErrID_Fatal,'Could not locate the MemberID referenced in the Members table in the associated Member-based Hydrodynamic coefficients table.',ErrStat,ErrMsg,RoutineName)
               RETURN
            END IF
         END IF

         IF ( InitInp%Morison%InpMembers(I)%PropPot .AND. InitInp%PotMod == 0  ) THEN
            CALL SetErrStat( ErrID_Fatal,'A member cannot have PropPot set to TRUE if PotMod = 0 in the FLOATING PLATFORM section.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF


         
      END DO

   END IF


   !-------------------------------------------------------------------------------------------------
   ! Filled Members Section
   !-------------------------------------------------------------------------------------------------

   IF ( InitInp%Morison%NFillGroups < 0 ) THEN
      CALL SetErrStat( ErrID_Fatal,'NFillGroups in the Filled-members table must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF

   IF ( InitInp%Morison%NFillGroups > 0 ) THEN

      DO I = 1,InitInp%Morison%NFillGroups

         IF ( InitInp%Morison%FilledGroups(I)%FillNumM < 1 ) THEN
            CALL SetErrStat( ErrID_Fatal,'FillNumM in the Filled-members table must be greater than zero.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF

         DO J = 1,InitInp%Morison%FilledGroups(I)%FillNumM

            DO K=1,InitInp%Morison%NMembers
               IF ( InitInp%Morison%FilledGroups(I)%FillMList(J) == InitInp%Morison%InpMembers(K)%MemberID ) THEN
                  FoundID = .TRUE.
                     ! Check to make sure this member is not already part of another fill group!
                  IF ( InitInp%Morison%InpMembers(K)%MmbrFilledIDIndx /= -1 ) THEN
                     CALL SetErrStat( ErrID_Fatal,'A member cannot be a part of more than one fill group!',ErrStat,ErrMsg,RoutineName)
                  END IF

                  InitInp%Morison%InpMembers(k)%MmbrFilledIDIndx = I

               END IF
            END DO

         END DO



            ! Make sure that the filled group members are connected
            ! NOTE: This would be easier if the input mesh was already a FAST Framework mesh because then you could use the mesh routines to determine connectivity.

            !InitInp%Morison%FilledGroups(I)%FillMList(J)

            ! Make sure the FillFSLoc is within one of the group members
            !InitInp%Morison%FilledGroups(I)%FillFSLoc


               ! Deal with DEFAULT or create a REAL from the string

         IF ( TRIM(InitInp%Morison%FilledGroups(I)%FillDensChr) /= 'DEFAULT' )  THEN

            READ (InitInp%Morison%FilledGroups(I)%FillDensChr,*,IOSTAT=IOS)  InitInp%Morison%FilledGroups(I)%FillDens
               CALL CheckIOS ( IOS, "", 'FillDens', NumType, ErrStat2, ErrMsg2 )
               CALL SetErrStat(ErrStat2, ErrMsg2,ErrStat,ErrMsg,RoutineName)
               IF ( ErrStat >= AbortErrLev ) RETURN
         ELSE
            InitInp%Morison%FilledGroups(I)%FillDens = InitInp%Waves%WtrDens
         END IF

      END DO

   END IF


   !-------------------------------------------------------------------------------------------------
   ! Marine Growth by Depth Section
   !-------------------------------------------------------------------------------------------------

   IF ( InitInp%Morison%NMGDepths < 0 ) THEN
      CALL SetErrStat( ErrID_Fatal,'NMGDepths in the Marine growth table must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF


   IF ( InitInp%Morison%NMGDepths > 0 ) THEN

      InitInp%Morison%MGTop    = -999999.0
      InitInp%Morison%MGBottom =  999999.0

      DO I = 1,InitInp%Morison%NMGDepths

            ! Store the boundaries of the marine growth zone
         IF ( InitInp%Morison%MGDepths(I)%MGDpth > InitInp%Morison%MGTop ) THEN
            InitInp%Morison%MGTop    = InitInp%Morison%MGDepths(I)%MGDpth
         END IF
         IF ( InitInp%Morison%MGDepths(I)%MGDpth < InitInp%Morison%MGBottom ) THEN
            InitInp%Morison%MGBottom = InitInp%Morison%MGDepths(I)%MGDpth
         ELSE
            CALL SetErrStat( ErrID_Fatal,'The rows of the marine growth table must be ordered with increasing depth (decreasing Z).',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF

            ! Make sure that the current MGDpth is not used elsewhere in the table.
         DO J = I+1,InitInp%Morison%NMGDepths
            IF ( EqualRealNos( InitInp%Morison%MGDepths(I)%MGDpth, InitInp%Morison%MGDepths(J)%MGDpth ) ) THEN
               CALL SetErrStat( ErrID_Fatal,'Duplicate MGDpth were found in the Marine Growth table.',ErrStat,ErrMsg,RoutineName)
               RETURN
            END IF
         END DO

         IF ( InitInp%Morison%MGDepths(I)%MGThck < 0 ) THEN
            CALL SetErrStat( ErrID_Fatal,'MGThck in the Marine growth table must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF
         IF ( InitInp%Morison%MGDepths(I)%MGDens < 0 ) THEN
            CALL SetErrStat( ErrID_Fatal,'MGDens in the Marine growth table must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF

      END DO

   END IF


   !-------------------------------------------------------------------------------------------------
   ! Member Output List Section
   !-------------------------------------------------------------------------------------------------

   IF ( ( InitInp%Morison%NMOutputs < 0 ) .OR. ( InitInp%Morison%NMOutputs > 9 ) ) THEN
      CALL SetErrStat( ErrID_Fatal,'NMOutputs in the Member output list must be greater or equal to zero and less than 10.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF

   IF ( InitInp%Morison%NMOutputs > 0 ) THEN


      DO I = 1,InitInp%Morison%NMOutputs

         InitInp%Morison%MOutLst(I)%MemberIDIndx = -1

            ! Find MemberID in this Member output list table in the Members table
         DO J = 1,InitInp%Morison%NMembers
            IF ( InitInp%Morison%InpMembers(J)%MemberID == InitInp%Morison%MOutLst(I)%MemberID ) THEN
               InitInp%Morison%MOutLst(I)%MemberIDIndx = J
            END IF
         END DO

            ! Make sure that a PropSetID entry in the Member cross-section properties table was found
         IF ( InitInp%Morison%MOutLst(I)%MemberIDIndx == -1 ) THEN
            CALL SetErrStat( ErrID_Fatal,'MemberID in the Member output list table does not appear in the Members table.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF

         IF ( ( InitInp%Morison%MOutLst(I)%NOutLoc < 1 ) .OR. ( InitInp%Morison%MOutLst(I)%NOutLoc > 9) ) THEN
            CALL SetErrStat( ErrID_Fatal,'NOutLoc in the Member output list must be greater than zero and less than 10.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF

         DO J = 1,InitInp%Morison%MOutLst(I)%NOutLoc
            IF ( ( InitInp%Morison%MOutLst(I)%NodeLocs(J) < 0.0 ) .OR. ( InitInp%Morison%MOutLst(I)%NodeLocs(J) > 1.0 ) ) THEN
               CALL SetErrStat( ErrID_Fatal,'NodeLocs in the Member output list must be greater or equal to 0.0 and less than or equal to 1.0.',ErrStat,ErrMsg,RoutineName)
               RETURN
            END IF
         END DO


      END DO

   END IF

   !-------------------------------------------------------------------------------------------------
   ! Joint Output List Section
   !-------------------------------------------------------------------------------------------------

   !IF ( InitInp%Morison%NJOutputs /= 0 ) THEN  ! TODO Remove this check and add back the other checks once Joint Outputs are supported
   !CALL SetErrStat( ErrID_Fatal,'NJOutputs in the Joint output list must be equal to zero.  Future versions of HydroDyn will support values greater or equal to zero and less than 10.'
   !   ErrStat = ErrID_Fatal
   !   RETURN
   !END IF   
   
   
   IF ( ( InitInp%Morison%NJOutputs < 0 ) .OR. ( InitInp%Morison%NMOutputs > 9 ) ) THEN
      CALL SetErrStat( ErrID_Fatal,'NJOutputs in the Joint output list must be greater or equal to zero and less than 10.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF

   IF ( InitInp%Morison%NJOutputs > 0 ) THEN


      DO I=1,InitInp%Morison%NJOutputs
           
           InitInp%Morison%JOutLst(I)%JointIDIndx = -1
         ! Find MemberID in this Member output list table in the Members table
         DO J = 1,InitInp%Morison%NJoints
            IF ( InitInp%Morison%InpJoints(J)%JointID == InitInp%Morison%JOutLst(I)%JointID ) THEN
               InitInp%Morison%JOutLst(I)%JointIDIndx = J
               EXIT
            END IF 
         END DO
         
            ! Make sure that a PropSetID entry in the Member cross-section properties table was found
         IF ( InitInp%Morison%JOutLst(I)%JointIDIndx == -1 ) THEN
            CALL SetErrStat( ErrID_Fatal,'JointID in the Joint output list table does not appear in the Joints table.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF
      END DO
   END IF
   
   !-------------------------------------------------------------------------------------------------
   ! Data section for OUTPUT
   !-------------------------------------------------------------------------------------------------


      ! OutAll - output all member and joint data

   IF ( InitInp%OutAll ) THEN    !TODO: Alter this check once OutAll is supported
         CALL SetErrStat( ErrID_Fatal,'OutAll must be FALSE. Future versions of HydroDyn will once again support values of either TRUE or FALSE.',ErrStat,ErrMsg,RoutineName)
         RETURN
   END IF


      ! OutSwtch - output file switch

   IF ( InitInp%OutSwtch /= 1 .AND. InitInp%OutSwtch /= 2 .AND. InitInp%OutSwtch /= 3 ) THEN
      CALL SetErrStat( ErrID_Fatal,'OutSwitch must be set to 1, 2, or 3.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF

   !InitInp%OutFmt
   !InitInp%OutSFmt


         ! OutList - list of requested parameters to output to a file


   !----------------------------------------------------------
   !  Output List
   !----------------------------------------------------------

      ! First we need to extract module-specific output lists from the user-input list.
      ! Any unidentified channels will be attached to the HydroDyn module's output list.
   IF (  InitInp%NUserOutputs > 0 ) THEN
      ALLOCATE ( foundMask(InitInp%NUserOutputs) , STAT = ErrStat2 )
      IF ( ErrStat2 /= ErrID_None ) THEN
         CALL SetErrStat( ErrID_Fatal,'Error allocating space for temporary array: foundMask in the HydroDynInput_GetInput subroutine.',ErrStat,ErrMsg,RoutineName)
         
         RETURN
      END IF
      foundMask = .FALSE.
         ! Extract Waves2 list
      InitInp%Waves2%NumOuts  = GetWaves2Channels   ( InitInp%NUserOutputs, InitInp%UserOutputs, InitInp%Waves2%OutList, foundMask, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   
         ! Extract WAMIT list
      InitInp%WAMIT%NumOuts   = GetWAMITChannels    ( InitInp%NUserOutputs, InitInp%UserOutputs, InitInp%WAMIT%OutList,  foundMask, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

         ! Extract WAMIT2 list
      InitInp%WAMIT2%NumOuts  = GetWAMIT2Channels   ( InitInp%NUserOutputs, InitInp%UserOutputs, InitInp%WAMIT2%OutList, foundMask, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

         ! Extract Morison list
         !foundMask = .FALSE.
      InitInp%Morison%NumOuts = GetMorisonChannels  ( InitInp%NUserOutputs, InitInp%UserOutputs, InitInp%Morison%OutList, foundMask, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   
         ! Attach remaining items to the HydroDyn list
         !foundMask = .FALSE.
      InitInp%NumOuts       = HDOut_GetChannels ( InitInp%NUserOutputs, InitInp%UserOutputs, InitInp%OutList        , foundMask, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      CALL PrintBadChannelWarning(InitInp%NUserOutputs, InitInp%UserOutputs , foundMask, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      IF (ErrStat >= AbortErrLev ) RETURN

      DEALLOCATE(foundMask)
   END IF
      ! Now that we have the sub-lists organized, lets do some additional validation.
   
   
   
   
   !----------------------------------------------------------
   ! Mesh-related Output List
   !----------------------------------------------------------

   IF ( InitInp%Morison%NumOuts > 0 ) THEN

         ! Create an  output list for validated outputs
      ALLOCATE ( InitInp%Morison%ValidOutList(InitInp%Morison%NumOuts), STAT = ErrStat2 )
      IF ( ErrStat2 /= 0 ) THEN
         CALL SetErrStat( ErrID_Fatal,'Error allocating valid output list array.',ErrStat,ErrMsg,RoutineName)
         RETURN
      END IF

      DO I =1, InitInp%Morison%NumOuts

         InitInp%Morison%ValidOutList(I) = CheckMeshOutput( InitInp%Morison%OutList(I), InitInp%Morison%NMOutputs, InitInp%Morison%MOutLst, InitInp%Morison%NJOutputs )

      END DO

   END IF


   !----------------------------------------------------------
   ! Populate data in sub-types from parent or other module types
   !----------------------------------------------------------

      ! Current
         ! For wave kinematic calculations, the effective water depth is the user input water depth (positive valued) + MSL2SWL (positive when SWL is above MSL).
      InitInp%Current%WtrDpth    = InitInp%Morison%WtrDpth + InitInp%Morison%MSL2SWL ! Adjust for the MSL2SWL.  
                                                       
      
      ! Waves
      InitInp%Waves%Gravity      = InitInp%Gravity
      InitInp%Waves%UnSum        = InitInp%UnSum
         ! For wave kinematic calculations, the effective water depth is the user input water depth (positive valued) + MSL2SWL (positive when SWL is above MSL).
      InitInp%Waves%WtrDpth      = InitInp%Morison%WtrDpth + InitInp%Morison%MSL2SWL ! Adjust for the MSL2SWL
      
      ! Waves2
      IF (InitInp%Waves2%WvDiffQTFF .OR. InitInp%Waves2%WvSumQTFF ) THEN
         InitInp%Waves2%WtrDens     = InitInp%Waves%WtrDens
         InitInp%Waves2%Gravity     = InitInp%Gravity
         InitInp%Waves2%UnSum       = InitInp%UnSum
         InitInp%Waves2%WtrDpth     = InitInp%Waves%WtrDpth
         InitInp%Waves2%WaveStMod   = InitInp%Waves%WaveStMod
         InitInp%Waves2%NWaveElev   = InitInp%Waves%NWaveElev
         CALL AllocAry( InitInp%Waves2%WaveElevxi, InitInp%Waves2%NWaveElev, 'WaveElevxi' , ErrStat2, ErrMsg2)
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
         CALL AllocAry( InitInp%Waves2%WaveElevyi, InitInp%Waves2%NWaveElev, 'WaveElevyi' , ErrStat2, ErrMsg2)
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
         IF ( ErrStat >= AbortErrLev ) RETURN 
         InitInp%Waves2%WaveElevxi  = InitInp%Waves%WaveElevxi
         InitInp%Waves2%WaveElevyi  = InitInp%Waves%WaveElevyi
      ENDIF

      ! WAMIT
      InitInp%WAMIT%WtrDens      = InitInp%Waves%WtrDens
      InitInp%WAMIT%WaveMod      = InitInp%Waves%WaveMod
      InitInp%WAMIT%OutAll       = InitInp%OutAll
      InitInp%WAMIT%HasWAMIT     = InitInp%PotMod == 1
      ! WAMIT2
      InitInp%WAMIT2%WtrDens     = InitInp%Waves%WtrDens
      InitInp%WAMIT2%WaveMod     = InitInp%Waves%WaveMod
      InitInp%WAMIT2%OutAll      = InitInp%OutAll
      InitInp%WAMIT2%HasWAMIT    = InitInp%PotMod == 1
      ! Morison
      InitInp%Morison%UnSum      = InitInp%UnSum
      InitInp%Morison%Gravity    = InitInp%Gravity
      InitInp%Morison%WtrDens    = InitInp%Waves%WtrDens
      InitInp%Morison%OutAll     = InitInp%OutAll

         ! Process the input geometry and generate the simulation mesh representation
      CALL Morison_ProcessMorisonGeometry( InitInp%Morison, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF ( ErrStat >= AbortErrLev ) RETURN

         ! Set the number and global Z locations for the X and Y components of the current velocities
      InitInp%Current%NMorisonNodes = InitInp%Morison%NNodes

      ALLOCATE ( InitInp%Current%MorisonNodezi(InitInp%Morison%NNodes), STAT = ErrStat2 )
      IF ( ErrStat2 /= ErrID_None ) THEN
         CALL SetErrStat( ErrID_Fatal,'Error allocating space for MorisonNodezi array.',ErrStat,ErrMsg,RoutineName)
         RETURN
      END IF



         ! Establish the number and locations where the wave kinematics will be computed
      InitInp%Waves%NWaveKin   = InitInp%Morison%NNodes                          ! Number of points where the incident wave kinematics will be computed (-)
      ALLOCATE ( InitInp%Waves%WaveKinxi(InitInp%Waves%NWaveKin), STAT = ErrStat2 )
      IF ( ErrStat2 /= ErrID_None ) THEN
         CALL SetErrStat( ErrID_Fatal,'Error allocating space for WaveKinxi array.',ErrStat,ErrMsg,RoutineName)

         RETURN
      END IF
      ALLOCATE ( InitInp%Waves%WaveKinyi(InitInp%Waves%NWaveKin), STAT = ErrStat2 )
      IF ( ErrStat2 /= ErrID_None ) THEN
         CALL SetErrStat( ErrID_Fatal,'Error allocating space for WaveKinyi array.',ErrStat,ErrMsg,RoutineName)
         RETURN
      END IF
      ALLOCATE ( InitInp%Waves%WaveKinzi(InitInp%Waves%NWaveKin), STAT = ErrStat2 )
      IF ( ErrStat2 /= ErrID_None ) THEN
         CALL SetErrStat( ErrID_Fatal,'Error allocating space for WaveKinzi array.',ErrStat,ErrMsg,RoutineName)
         RETURN
      END IF
      DO I=1,InitInp%Morison%NNodes
         InitInp%Waves%WaveKinxi(I)      = InitInp%Morison%Nodes(I)%JointPos(1)                          ! xi-coordinates for points where the incident wave kinematics will be computed;
         InitInp%Waves%WaveKinyi(I)      = InitInp%Morison%Nodes(I)%JointPos(2)                          ! yi-coordinates for points where the incident wave kinematics will be computed;
         InitInp%Waves%WaveKinzi(I)      = InitInp%Morison%Nodes(I)%JointPos(3) - InitInp%Morison%MSL2SWL   ! zi-coordinates for points where the incident wave kinematics will be computed, adjusted to the still water level(meters)     
         InitInp%Current%MorisonNodezi(I) = InitInp%Waves%WaveKinzi(I)
      END DO


            ! If we are using the Waves module, the node information must be copied over.
      InitInp%Waves2%NWaveKin   = InitInp%Waves%NWaveKin                          ! Number of points where the incident wave kinematics will be computed (-)
      IF ( InitInp%Waves2%WvDiffQTFF .OR. InitInp%Waves2%WvSumQTFF ) THEN
         ALLOCATE ( InitInp%Waves2%WaveKinxi(InitInp%Waves2%NWaveKin), STAT = ErrStat2 )
         IF ( ErrStat2 /= ErrID_None ) THEN
            CALL SetErrStat( ErrID_Fatal,'Error allocating space for WaveKinxi array for Waves2 module.',ErrStat,ErrMsg,RoutineName)

            RETURN
         END IF
         ALLOCATE ( InitInp%Waves2%WaveKinyi(InitInp%Waves2%NWaveKin), STAT = ErrStat2 )
         IF ( ErrStat2 /= ErrID_None ) THEN
            CALL SetErrStat( ErrID_Fatal,'Error allocating space for WaveKinyi array for Waves2 module.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF
         ALLOCATE ( InitInp%Waves2%WaveKinzi(InitInp%Waves2%NWaveKin), STAT = ErrStat2 )
         IF ( ErrStat2 /= ErrID_None ) THEN
            CALL SetErrStat( ErrID_Fatal,'Error allocating space for WaveKinzi array for Waves2 module.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF

         InitInp%Waves2%WaveKinxi  = InitInp%Waves%WaveKinxi
         InitInp%Waves2%WaveKinyi  = InitInp%Waves%WaveKinyi
         InitInp%Waves2%WaveKinzi  = InitInp%Waves%WaveKinzi

      ENDIF

END SUBROUTINE HydroDynInput_ProcessInitData

END MODULE HydroDyn_Input
