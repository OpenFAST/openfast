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
   USE                              Waves2_Output
   USE                              Morison
   USE                              Morison_Output
   USE                              NWTC_RandomNumber
   IMPLICIT                         NONE

   PRIVATE :: CheckMeshOutput

CONTAINS
   
!====================================================================================================
FUNCTION CheckMeshOutput( output, numMemberOut, MOutLst, numJointOut )
!     The routine
!----------------------------------------------------------------------------------------------------
!
   CHARACTER(ChanLen),        INTENT ( IN    )  :: output
   INTEGER,                   INTENT ( IN    )  :: numMemberOut
   TYPE(Morison_MOutput),     INTENT ( IN    )  :: MOutLst(:)
   INTEGER,                   INTENT ( IN    )  :: numJointOut
   !INTEGER,                   INTENT (   OUT )  :: ErrStat              ! returns a non-zero value when an error occurs
   !CHARACTER(*),              INTENT (   OUT )  :: ErrMsg               ! Error message if ErrStat /= ErrID_None

   LOGICAL                                      :: CheckMeshOutput

   INTEGER                                      :: ErrStat
   CHARACTER(ChanLen)                           :: outputTmp
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
   CHARACTER(ChanLen),            INTENT( IN    ) :: UserOutputs (:)      ! An array holding the names of the requested output channels. 
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
SUBROUTINE HydroDyn_ParseInput( InputFileName, OutRootName, defWtrDens, defWtrDpth, defMSL2SWL, FileInfo_In, InputFileData, ErrStat, ErrMsg )
!     This public subroutine reads the input required for HydroDyn from the file whose name is an
!     input parameter.
!----------------------------------------------------------------------------------------------------

      ! Passed variables
   CHARACTER(*),                  intent(in   ) :: InputFileName        !< The name of the input file, for putting in echo file.
   CHARACTER(*),                  intent(in   ) :: OutRootName          !< The rootname of the echo file, possibly opened in this routine
   real(ReKi),                    intent(in   ) :: defWtrDens           !< default value for water density
   real(ReKi),                    intent(in   ) :: defWtrDpth           !< default value for water depth
   real(ReKi),                    intent(in   ) :: defMSL2SWL           !< default value for mean sea level to still water level
   TYPE(FileInfoType),            INTENT(IN   ) :: FileInfo_In          !< The derived type for holding the file information
   TYPE(HydroDyn_InputFile),      INTENT(INOUT) :: InputFileData        ! the hydrodyn input file data
   INTEGER,                       INTENT(  OUT) :: ErrStat              ! returns a non-zero value when an error occurs
   CHARACTER(*),                  INTENT(  OUT) :: ErrMsg               ! Error message if ErrStat /= ErrID_None

      ! Local variables
   INTEGER                                      :: I, j                 ! generic integer for counting
   CHARACTER(   2)                              :: strI                 ! string version of the loop counter
   INTEGER                                      :: UnEc                 ! The local unit number for this module's echo file
   CHARACTER(1024)                              :: EchoFile             ! Name of HydroDyn echo file
   CHARACTER(1024)                              :: Line                 ! String to temporarially hold value of read line
   real(ReKi), ALLOCATABLE                      :: tmpVec1(:), tmpVec2(:) ! Temporary arrays for WAMIT data
   integer(IntKi)                               :: startIndx, endIndx   ! indices into working arrays
   INTEGER, ALLOCATABLE                         :: tmpArray(:)          ! Temporary array storage of the joint output list
   REAL(ReKi), ALLOCATABLE                      :: tmpReArray(:)        ! Temporary array storage of the joint output list
   CHARACTER(1)                                 :: Line1                ! The first character of an input line
   INTEGER(IntKi)                               :: CurLine              !< Current entry in FileInfo_In%Lines array
   INTEGER(IntKi)                               :: ErrStat2
   CHARACTER(ErrMsgLen)                         :: ErrMsg2
   CHARACTER(*),  PARAMETER                     :: RoutineName = 'HydroDyn_ParaseInput'
   
      ! Initialize local data
   UnEc     = -1
   ErrStat  =  ErrID_None         
   ErrMsg   =  ""   
   InputFileData%Echo = .FALSE.  ! initialize for error handling (cleanup() routine)
   

   !-------------------------------------------------------------------------------------------------
   ! General settings 
   !-------------------------------------------------------------------------------------------------

   CurLine = 3    ! Skip the first three lines as they are known to be header lines and separators
   call ParseVar( FileInfo_In, CurLine, 'Echo', InputFileData%Echo, ErrStat2, ErrMsg2 )
         if (Failed()) return;

   if ( InputFileData%Echo ) then
      EchoFile = TRIM(OutRootName)//'.HD.ech'
      CALL OpenEcho ( UnEc, TRIM(EchoFile), ErrStat2, ErrMsg2 )
         if (Failed())  return;
      WRITE(UnEc, '(A)') 'Echo file for AeroDyn 15 primary input file: '//trim(InputFileName)
      ! Write the first three lines into the echo file
      WRITE(UnEc, '(A)') trim(FileInfo_In%Lines(1))
      WRITE(UnEc, '(A)') trim(FileInfo_In%Lines(2))

      CurLine = 3
      call ParseVar( FileInfo_In, CurLine, 'Echo', InputFileData%Echo, ErrStat2, ErrMsg2, UnEc )
         if (Failed()) return
   endif


   !-------------------------------------------------------------------------------------------------
   ! Environmental conditions section
   !-------------------------------------------------------------------------------------------------
   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') trim(FileInfo_In%Lines(CurLine))    ! Write section break to echo
   CurLine = CurLine + 1

      ! WtrDens - Water density.
   CALL ParseVarWDefault ( FileInfo_In, CurLine, 'WtrDens', InputFileData%Waves%WtrDens, real(defWtrDens,SiKi), ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! WtrDpth - Water depth
   CALL ParseVarWDefault ( FileInfo_In, CurLine, 'WtrDpth', InputFileData%Morison%WtrDpth, defWtrDpth, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! MSL2SWL
   CALL ParseVarWDefault ( FileInfo_In, CurLine, 'MSL2SWL', InputFileData%Morison%MSL2SWL, defMSL2SWL, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;


   !-------------------------------------------------------------------------------------------------
   ! Data section for waves
   !-------------------------------------------------------------------------------------------------
   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') trim(FileInfo_In%Lines(CurLine))    ! Write section break to echo
   CurLine = CurLine + 1

      ! WaveMod - Wave kinematics model switch.
   call ParseVar( FileInfo_In, CurLine, 'WaveMod', InputFileData%Waves%WaveModChr, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

   CALL Conv2UC( InputFileData%Waves%WaveModChr )    ! Convert Line to upper case.

   InputFileData%Waves%WavePhase = 0.0
   InputFileData%Waves%WaveNDAmp = .FALSE.


      ! WaveStMod - Model switch for stretching incident wave kinematics to instantaneous free surface.
   call ParseVar( FileInfo_In, CurLine, 'WaveStMod', InputFileData%Waves%WaveStMod, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! WaveTMax - Analysis time for incident wave calculations.
   call ParseVar( FileInfo_In, CurLine, 'WaveTMax', InputFileData%Waves%WaveTMax, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! WaveDT - Time step for incident wave calculations
   call ParseVar( FileInfo_In, CurLine, 'WaveDT', InputFileData%Waves%WaveDT, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! WaveHs - Significant wave height
   call ParseVar( FileInfo_In, CurLine, 'WaveHs', InputFileData%Waves%WaveHs, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! WaveTp - Peak spectral period.
   call ParseVar( FileInfo_In, CurLine, 'WaveTp', InputFileData%Waves%WaveTp, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! WavePkShp - Peak shape parameter.
   call ParseVar( FileInfo_In, CurLine, 'WavePkShp', InputFileData%Waves%WavePkShpChr, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! WvLowCOff - Low Cut-off frequency or lower frequency limit of the wave spectrum beyond which the wave spectrum is zeroed (rad/s).  
   call ParseVar( FileInfo_In, CurLine, 'WvLowCOff', InputFileData%Waves%WvLowCOff, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

     ! WvHiCOff - High Cut-off frequency or upper frequency limit of the wave spectrum beyond which the wave spectrum is zeroed (rad/s).  
   call ParseVar( FileInfo_In, CurLine, 'WvHiCOff', InputFileData%Waves%WvHiCOff, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;
   
      ! WaveDir - Mean wave heading direction.
   call ParseVar( FileInfo_In, CurLine, 'WaveDir', InputFileData%Waves%WaveDir, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! WaveDirMod -  Directional spreading function {0: None, 1: COS2S}       (-) [Used only if WaveMod=2]
   call ParseVar( FileInfo_In, CurLine, 'WaveDirMod', InputFileData%Waves%WaveDirMod, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! WaveDirSpread -  Spreading coefficient [only used if WaveMod=2 and WaveDirMod=1]
   call ParseVar( FileInfo_In, CurLine, 'WaveDirSpread', InputFileData%Waves%WaveDirSpread, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! WaveNDir -  The number of wave directions to calculate [must be odd; only used if WaveDirMod=1]
   call ParseVar( FileInfo_In, CurLine, 'WaveNDir', InputFileData%Waves%WaveNDir, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! WaveDirRange - Full range of the wave directions from WaveDir - WaveDirRange/2 to WaveDir + WaveDirRange/2 (only used if WaveMod=2 and WaveDirMod=1)
   call ParseVar( FileInfo_In, CurLine, 'WaveDirRange', InputFileData%Waves%WaveDirRange, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! Negative values should be treated as positive.
   InputFileData%Waves%WaveDirRange =  ABS( InputFileData%Waves%WaveDirRange )


      ! WaveSeed(1)
   call ParseVar( FileInfo_In, CurLine, 'WaveSeed(1)', InputFileData%Waves%WaveSeed(1), ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;
   InputFileData%Waves%RNG%RandSeed(1) = InputFileData%Waves%WaveSeed(1)

      !WaveSeed(2)
   call ParseVar( FileInfo_In, CurLine, 'WaveSeed(2)', Line, ErrStat2, ErrMsg2, UnEc )    ! Read into a string and then parse
      if (Failed())  return;

   READ (Line,*,IOSTAT=ErrStat2) Line1  ! check the first character to make sure we don't have T/F, which can be interpreted as 1/-1 or 0 in Fortran
   CALL Conv2UC( Line1 )
   IF ( (Line1 == 'T') .OR. (Line1 == 'F') ) THEN
      ErrStat2 = ErrID_Fatal
      ErrMsg2  = ' WaveSeed(2): Invalid RNG type.'
      if (Failed())  return;
   ENDIF

!FIXME: there is something a little strange here.  RandSeed(2) is an integer, but what if we get an error on the next read?
   READ (Line,*,IOSTAT=ErrStat2) InputFileData%Waves%WaveSeed(2)
   InputFileData%Waves%RNG%RandSeed(2) = InputFileData%Waves%WaveSeed(2)

   IF (ErrStat2 == 0) THEN ! the user entered a number
      InputFileData%Waves%RNG%RNG_type = "NORMAL"
      InputFileData%Waves%RNG%pRNG = pRNG_INTRINSIC

   ELSE

      InputFileData%Waves%RNG%RNG_type = ADJUSTL( Line )
      CALL Conv2UC( InputFileData%Waves%RNG%RNG_type )

      IF ( InputFileData%Waves%RNG%RNG_type == "RANLUX") THEN
         InputFileData%Waves%RNG%pRNG = pRNG_RANLUX
      ELSE
         ErrStat2 = ErrID_Fatal
         ErrMsg2  = ' WaveSeed(2): Invalid alternative random number generator.'
         if (Failed())  return;
      ENDIF

   ENDIF


      ! WaveNDAmp - Flag for normally distributed amplitudes.
   call ParseVar( FileInfo_In, CurLine, 'WaveNDAmp', InputFileData%Waves%WaveNDAmp, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! WvKinFile
   call ParseVar( FileInfo_In, CurLine, 'WvKinFile', InputFileData%Waves%WvKinFile, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! NWaveElev
   call ParseVar( FileInfo_In, CurLine, 'NWaveElev', InputFileData%Waves%NWaveElev, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! This check is needed here instead of being located in HydroDynInput_ProcessInputData() because
      ! we need to allocate arrays.  If _GetInput() was skipped, then these array would already have
      ! been allocated and populated.

   IF ( InputFileData%Waves%NWaveElev < 0 .OR. InputFileData%Waves%NWaveElev > 9 ) THEN
      ErrStat2 = ErrID_Fatal
      ErrMsg2  = 'NWaveElev must be greater than or equal to zero and less than 10.'
      if (Failed())  return;
   END IF

      ! allocate space for the output location arrays:
   CALL AllocAry( InputFileData%Waves%WaveElevxi, InputFileData%Waves%NWaveElev, 'WaveElevxi' , ErrStat2, ErrMsg2);  if (Failed())  return;
   CALL AllocAry( InputFileData%Waves%WaveElevyi, InputFileData%Waves%NWaveElev, 'WaveElevyi' , ErrStat2, ErrMsg2);  if (Failed())  return;
      
      ! WaveElevxi
   call ParseAry ( FileInfo_In, CurLine, 'WaveElevxi.', InputFileData%Waves%WaveElevxi, InputFileData%Waves%NWaveElev, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! WaveElevyi
   call ParseAry ( FileInfo_In, CurLine, 'WaveElevyi.', InputFileData%Waves%WaveElevyi, InputFileData%Waves%NWaveElev, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;


   !-------------------------------------------------------------------------------------------------
   ! Data section for 2nd Order Waves 
   !-------------------------------------------------------------------------------------------------
   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') trim(FileInfo_In%Lines(CurLine))    ! Write section break to echo
   CurLine = CurLine + 1

      ! WvDiffQTFF     - Second order waves -- difference forces
   call ParseVar( FileInfo_In, CurLine, 'WvDiffQTF', InputFileData%Waves2%WvDiffQTFF, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! WvSumQTFF      - Second order waves -- sum forces
   call ParseVar( FileInfo_In, CurLine, 'WvSumQTF', InputFileData%Waves2%WvSumQTFF, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! WvLowCOffD   -- Minimum frequency used in the difference methods (rad/s)              [Only used if DiffQTF /= 0]
   call ParseVar( FileInfo_In, CurLine, 'WvLowCOffD', InputFileData%Waves2%WvLowCOffD, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! WvHiCOffD   -- Maximum frequency used in the difference methods  (rad/s)              [Only used if DiffQTF /= 0]
   call ParseVar( FileInfo_In, CurLine, 'WvHiCOffD', InputFileData%Waves2%WvHiCOffD, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! WvLowCOffS   -- Minimum frequency used in the        sum-QTF     (rad/s)              [Only used if  SumQTF /= 0]
   call ParseVar( FileInfo_In, CurLine, 'WvLowCOffS', InputFileData%Waves2%WvLowCOffS, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! WvHiCOffS   -- Maximum frequency used in the        sum-QTF      (rad/s)              [Only used if  SumQTF /= 0]
   call ParseVar( FileInfo_In, CurLine, 'WvHiCOffS', InputFileData%Waves2%WvHiCOffS, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;


   !-------------------------------------------------------------------------------------------------
   ! Data section for current
   !-------------------------------------------------------------------------------------------------
   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') trim(FileInfo_In%Lines(CurLine))    ! Write section break to echo
   CurLine = CurLine + 1

      ! CurrMod - Current profile model switch
   call ParseVar( FileInfo_In, CurLine, 'CurrMod', InputFileData%Current%CurrMod, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! CurrSSV0 - Sub-surface current velocity at still water level
   call ParseVar( FileInfo_In, CurLine, 'CurrSSV0', InputFileData%Current%CurrSSV0, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;


      ! CurrSSDirChr - Sub-surface current heading direction
   call ParseVar( FileInfo_In, CurLine, 'CurrSSDir', InputFileData%Current%CurrSSDirChr, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

   CALL Conv2UC( InputFileData%Current%CurrSSDirChr )    ! Convert Line to upper case.


      ! CurrNSRef - Near-surface current reference depth.
   call ParseVar( FileInfo_In, CurLine, 'CurrNSRef', InputFileData%Current%CurrNSRef, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! CurrNSV0 - Near-surface current velocity at still water level.
   call ParseVar( FileInfo_In, CurLine, 'CurrNSV0', InputFileData%Current%CurrNSV0, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! CurrNSDir - Near-surface current heading direction.
   call ParseVar( FileInfo_In, CurLine, 'CurrNSDir', InputFileData%Current%CurrNSDir, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! CurrDIV - Depth-independent current velocity.
   call ParseVar( FileInfo_In, CurLine, 'CurrDIV', InputFileData%Current%CurrDIV, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! CurrDIDir - Depth-independent current heading direction.
   call ParseVar( FileInfo_In, CurLine, 'CurrDIDir', InputFileData%Current%CurrDIDir, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;


   !-------------------------------------------------------------------------------------------------
   ! Data section for floating platform
   !-------------------------------------------------------------------------------------------------
   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') trim(FileInfo_In%Lines(CurLine))    ! Write section break to echo
   CurLine = CurLine + 1

      ! PotMod - State indicating potential flow model used in the simulation. 0=none, 1=WAMIT, 2=FIT
   call ParseVar( FileInfo_In, CurLine, 'PotMod', InputFileData%PotMod, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! ExctnMod  - Wave Excitation model {0: None, 1: DFT, 2: state-space} (switch)
      ! [STATE-SPACE REQUIRES *.ssexctn INPUT FILE]
   call ParseVar( FileInfo_In, CurLine, 'ExctnMod', InputFileData%WAMIT%ExctnMod, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! RdtnMod  - Radiation memory-effect model {1: convolution, 2: state-space} (switch)
      ! [STATE-SPACE REQUIRES *.ss INPUT FILE]
   call ParseVar( FileInfo_In, CurLine, 'RdtnMod', InputFileData%WAMIT%RdtnMod, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! RdtnTMax - Analysis time for wave radiation kernel calculations
      ! NOTE: Use RdtnTMax = 0.0 to eliminate wave radiation damping
   call ParseVar( FileInfo_In, CurLine, 'RdtnTMax', InputFileData%WAMIT%RdtnTMax, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! RdtnDT - Time step for wave radiation kernel calculations
   call ParseVar( FileInfo_In, CurLine, 'RdtnDT', InputFileData%WAMIT%Conv_Rdtn%RdtnDTChr, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! NBody - Number of WAMIT bodies to be used (-) [>=1; only used when PotMod=1. If NBodyMod=1, the WAMIT data 
      !         contains a vector of size 6*NBody x 1 and matrices of size 6*NBody x 6*NBody; if NBodyMod>1, there 
      !         are NBody sets of WAMIT data each with a vector of size 6 x 1 and matrices of size 6 x 6]
   call ParseVar( FileInfo_In, CurLine, 'NBody', InputFileData%NBody, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! NBodyMod - Body coupling model {1: include coupling terms between each body and NBody in HydroDyn equals NBODY in WAMIT, 
      !            2: neglect coupling terms between each body and NBODY=1 with XBODY=0 in WAMIT, 3: Neglect coupling terms 
      !            between each body and NBODY=1 with XBODY=/0 in WAMIT} (switch) [only used when PotMod=1]
   call ParseVar( FileInfo_In, CurLine, 'NBodyMod', InputFileData%NBodyMod, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;
      
         ! allocate space for the WAMIT-related data arrays:
   if ( InputFileData%NBodyMod == 1 )  then
      InputFileData%nWAMITObj = 1  ! Special case where all data in a single WAMIT input file as opposed to  InputFileData%NBody number of separate input files.
      InputFileData%vecMultiplier = InputFileData%NBody
   else
      InputFileData%nWAMITObj    = InputFileData%NBody
      InputFileData%vecMultiplier = 1
   end if
   
   CALL AllocAry( InputFileData%PotFile      , InputFileData%nWAMITObj, 'PotFile'      , ErrStat2, ErrMsg2);   if (Failed())  return;
   CALL AllocAry( InputFileData%WAMITULEN    , InputFileData%nWAMITObj, 'WAMITULEN'    , ErrStat2, ErrMsg2);   if (Failed())  return;
   CALL AllocAry( InputFileData%PtfmRefxt    , InputFileData%NBody,     'PtfmRefxt'    , ErrStat2, ErrMsg2);   if (Failed())  return;
   CALL AllocAry( InputFileData%PtfmRefyt    , InputFileData%NBody,     'PtfmRefyt'    , ErrStat2, ErrMsg2);   if (Failed())  return;
   CALL AllocAry( InputFileData%PtfmRefzt    , InputFileData%NBody,     'PtfmRefzt'    , ErrStat2, ErrMsg2);   if (Failed())  return;
   CALL AllocAry( InputFileData%PtfmRefztRot , InputFileData%NBody,     'PtfmRefztRot' , ErrStat2, ErrMsg2);   if (Failed())  return;
   CALL AllocAry( InputFileData%PtfmVol0     , InputFileData%NBody,     'PtfmVol0'     , ErrStat2, ErrMsg2);   if (Failed())  return;
   CALL AllocAry( InputFileData%PtfmCOBxt    , InputFileData%NBody,     'PtfmCOBxt'    , ErrStat2, ErrMsg2);   if (Failed())  return;
   CALL AllocAry( InputFileData%PtfmCOByt    , InputFileData%NBody,     'PtfmCOByt'    , ErrStat2, ErrMsg2);   if (Failed())  return;


      ! PotFile - Root name of Potential flow data files (Could be WAMIT files or the FIT input file)
   call ParseAry( FileInfo_In, CurLine, 'PotFile', InputFileData%PotFile, InputFileData%nWAMITObj, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! WAMITULEN - WAMIT characteristic body length scale
   call ParseAry( FileInfo_In, CurLine, 'WAMITULEN', InputFileData%WAMITULEN, InputFileData%nWAMITObj, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! PtfmRefxt  - The xt offset of the body reference point(s) from (0,0,0) (meters)
   call ParseAry( FileInfo_In, CurLine, 'PtfmRefxt', InputFileData%PtfmRefxt, InputFileData%NBody, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! PtfmRefyt  - The yt offset of the body reference point(s) from (0,0,0) (meters)
   call ParseAry( FileInfo_In, CurLine, 'PtfmRefyt', InputFileData%PtfmRefyt, InputFileData%NBody, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! PtfmRefzt  - The zt offset of the body reference point(s) from (0,0,0) (meters)
   call ParseAry( FileInfo_In, CurLine, 'PtfmRefzt', InputFileData%PtfmRefzt, InputFileData%NBody, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! PtfmRefztRot  - The rotation about zt of the body reference frame(s) from xt/yt (deg)
   call ParseAry( FileInfo_In, CurLine, 'PtfmRefztRot', InputFileData%PtfmRefztRot, InputFileData%NBody, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;
   InputFileData%PtfmRefztRot = InputFileData%PtfmRefztRot*D2R_D ! Convert to radians
   
      ! PtfmVol0 - Displaced volume of water when the platform is in its undisplaced position
   call ParseAry( FileInfo_In, CurLine, 'PtfmVol0', InputFileData%PtfmVol0, InputFileData%NBody, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! PtfmCOBxt  - The xt offset of the center of buoyancy (COB) from the WAMIT reference point
   call ParseAry( FileInfo_In, CurLine, 'PtfmCOBxt', InputFileData%PtfmCOBxt, InputFileData%NBody, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! PtfmCOByt - The yt offset of the center of buoyancy (COB) from the WAMIT reference point
   call ParseAry( FileInfo_In, CurLine, 'PtfmCOByt', InputFileData%PtfmCOByt, InputFileData%NBody, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;
   

!bjj: should we add this?
!test for numerical stability
!      IF ( FP_InitData%RdtnDT <= FP_InitData%RdtnTMax*EPSILON(FP_InitData%RdtnDT) )  THEN  ! Test RdtnDT and RdtnTMax to ensure numerical stability -- HINT: see the use of OnePlusEps."
!         ErrStat = ErrID_Fatal
!         ErrMsg2 = ' RdtnDT must be greater than '//TRIM ( Num2LStr( RdtnTMax*EPSILON(RdtnDT) ) )//' seconds.'
!         if (Failed())  return;
!      END IF



   !-------------------------------------------------------------------------------------------------
   ! Data section for 2nd order WAMIT forces
   !-------------------------------------------------------------------------------------------------
   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') trim(FileInfo_In%Lines(CurLine))    ! Write section break to echo
   CurLine = CurLine + 1

        ! MnDrift    -- Mean drift forces computed from WAMIT file: {0: No mean drift, [7, 8, 9, 10, 11, or 12]: WAMIT file to use}
   call ParseVar( FileInfo_In, CurLine, 'MnDrift', InputFileData%WAMIT2%MnDrift, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

        ! NewmanApp  -- Slow drift forces computed with Newman's approximation from  WAMIT file: {0: No mean drift, [7, 8, 9, 10, 11, or 12]: WAMIT file to use}
   call ParseVar( FileInfo_In, CurLine, 'NewmanApp', InputFileData%WAMIT2%NewmanApp, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

        ! DiffQTF    -- Full Difference-Frequency forces computed with full QTFs from WAMIT file: {0: No difference-frequency forces, [10, 11, or 12]: WAMIT file to use} -- Only one of MnDrift, NewmanApp, or DiffQYT can be non-zero
   call ParseVar( FileInfo_In, CurLine, 'DiffQTF', InputFileData%WAMIT2%DiffQTF, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

        ! SumQTF     -- Full        Sum-Frequency forces computed with full QTFs from WAMIT file: {0: No        Sum-frequency forces, [10, 11, or 12]: WAMIT file to use}
   call ParseVar( FileInfo_In, CurLine, 'SumQTF', InputFileData%WAMIT2%SumQTF, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;


   !-------------------------------------------------------------------------------------------------
   ! Floating Platform Additional Stiffness and Damping Section
   !-------------------------------------------------------------------------------------------------
   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') trim(FileInfo_In%Lines(CurLine))    ! Write section break to echo
   CurLine = CurLine + 1

   ! If NBodyMod = 1 then vecMultiplier = NBody and nWAMITObj = 1
   ! Else                 vecMultiplier = 1     and nWAMITObj = NBody
   CALL AllocAry( InputFileData%AddF0,     InputFileData%vecMultiplier*6, InputFileData%nWAMITObj, 'InputFileData%AddF0'    , ErrStat2, ErrMsg2);                                   if (Failed())  return; 
   CALL AllocAry( InputFileData%AddCLin,   InputFileData%vecMultiplier*6, InputFileData%vecMultiplier*6, InputFileData%nWAMITObj, 'InputFileData%AddCLin'  , ErrStat2, ErrMsg2);    if (Failed())  return; 
   CALL AllocAry( InputFileData%AddBLin,   InputFileData%vecMultiplier*6, InputFileData%vecMultiplier*6, InputFileData%nWAMITObj, 'InputFileData%AddBLin'  , ErrStat2, ErrMsg2);    if (Failed())  return;
   CALL AllocAry( InputFileData%AddBQuad,  InputFileData%vecMultiplier*6, InputFileData%vecMultiplier*6, InputFileData%nWAMITObj, 'InputFileData%AddBQuad' , ErrStat2, ErrMsg2);    if (Failed())  return;
   CALL AllocAry( tmpVec1, InputFileData%nWAMITObj, 'tmpVec1', ErrStat2, ErrMsg2);  if (Failed())  return;
   CALL AllocAry( tmpVec2, 6*InputFileData%NBody,   'tmpVec2', ErrStat2, ErrMsg2);  if (Failed())  return;

      ! AddF0 - Additional preload
   do i = 1,6*InputFileData%vecMultiplier   
      call ParseAry( FileInfo_In, CurLine, 'AddF0', tmpVec1, InputFileData%nWAMITObj, ErrStat2, ErrMsg2, UnEc )
         if (Failed())  return;

      do j = 1, InputFileData%nWAMITObj
         InputFileData%AddF0(i,j) = tmpVec1(j)
      end do
   end do
   
      ! AddCLin
   do i=1,6*InputFileData%vecMultiplier

      write(strI,'(I1)') i
      call ParseAry( FileInfo_In, CurLine, ' Row '//strI//' of the additional linear stiffness matrix', &
                     tmpVec2, 6*InputFileData%NBody, ErrStat2, ErrMsg2, UnEc )
         if (Failed())  return;

      do j = 1, InputFileData%nWAMITObj
         startIndx = 6*InputFileData%vecMultiplier*(j-1) + 1
         endIndx   = startIndx + 6*InputFileData%vecMultiplier - 1
         InputFileData%AddCLin(i,:,j) = tmpVec2(startIndx:endIndx)
      end do
   end do


       ! AddBLin
   DO I=1,6*InputFileData%vecMultiplier

      call ParseAry( FileInfo_In, CurLine, ' Row '//strI//' of the additional linear damping matrix', &
                     tmpVec2, 6*InputFileData%NBody, ErrStat2, ErrMsg2, UnEc )
         if (Failed())  return;

      do j = 1, InputFileData%nWAMITObj
         startIndx = 6*InputFileData%vecMultiplier*(j-1) + 1
         endIndx   = startIndx + 6*InputFileData%vecMultiplier - 1
         InputFileData%AddBLin(I,:,j) = tmpVec2(startIndx:endIndx)
      end do
   END DO


       ! AddBQuad
   DO I=1,6*InputFileData%vecMultiplier

      call ParseAry( FileInfo_In, CurLine, ' Row '//strI//' of the additional quadratic damping matrix', &
                     tmpVec2, 6*InputFileData%NBody, ErrStat2, ErrMsg2, UnEc )
         if (Failed())  return;

      do j = 1, InputFileData%nWAMITObj
         startIndx = 6*InputFileData%vecMultiplier*(j-1) + 1
         endIndx   = startIndx + 6*InputFileData%vecMultiplier - 1
         InputFileData%AddBQuad(I,:,j) = tmpVec2(startIndx:endIndx)
      end do
   END DO


   !-------------------------------------------------------------------------------------------------
   !  Axial Coefficients Section
   !-------------------------------------------------------------------------------------------------
   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') trim(FileInfo_In%Lines(CurLine))    ! Write section break to echo
   CurLine = CurLine + 1
   
      ! NAxCoef - Number of axial coefficients
   call ParseVar( FileInfo_In, CurLine, 'NAxCoef', InputFileData%Morison%NAxCoefs, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;
   
      ! Table header
   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') 'Axial coefficient table header line 1: '//NewLine//trim(FileInfo_In%Lines(CurLine))
   CurLine = CurLine + 1
   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') 'Axial coefficient table header line 2: '//NewLine//trim(FileInfo_In%Lines(CurLine))
   CurLine = CurLine + 1
  
   IF ( InputFileData%Morison%NAxCoefs > 0 ) THEN
      CALL AllocAry( tmpReArray, 4, 'temporary array for AxialCoefs', ErrStat2, ErrMsg2 )
         if (Failed())  return;
      
         ! Allocate memory for Axial Coef-related arrays
      ALLOCATE ( InputFileData%Morison%AxialCoefs(InputFileData%Morison%NAxCoefs), STAT = ErrStat2 )
      IF ( ErrStat2 /= 0 ) THEN
         ErrStat2 = ErrID_Fatal
         ErrMsg2  = 'Error allocating space for AxialCoefs array.'
         if (Failed())  return;
      END IF
          
      DO I = 1,InputFileData%Morison%NAxCoefs
            ! read the table entries   AxCoefID   CdAx  CaAx    in the HydroDyn input file
         call ParseAry( FileInfo_In, CurLine, ' axial coefficients line '//trim( Int2LStr(I)), tmpReArray, size(tmpReArray), ErrStat2, ErrMsg2, UnEc )
            if (Failed())  return;
         InputFileData%Morison%AxialCoefs(I)%AxCoefID = NINT(tmpReArray(1))
         InputFileData%Morison%AxialCoefs(I)%AxCd     =      tmpReArray(2)
         InputFileData%Morison%AxialCoefs(I)%AxCa     =      tmpReArray(3)
         InputFileData%Morison%AxialCoefs(I)%AxCp     =      tmpReArray(4)
      END DO

      if (allocated(tmpReArray))      deallocate(tmpReArray)
   END IF

   
   !-------------------------------------------------------------------------------------------------
   ! Member Joints Section
   !-------------------------------------------------------------------------------------------------
   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') trim(FileInfo_In%Lines(CurLine))    ! Write section break to echo
   CurLine = CurLine + 1

      ! NJoints - Number of member joints
   call ParseVar( FileInfo_In, CurLine, 'NJoints', InputFileData%Morison%NJoints, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! Table header
   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') 'Joints table header line 1: '//NewLine//trim(FileInfo_In%Lines(CurLine))
   CurLine = CurLine + 1
   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') 'Joints table header line 2: '//NewLine//trim(FileInfo_In%Lines(CurLine))
   CurLine = CurLine + 1

   IF ( InputFileData%Morison%NJoints > 0 ) THEN
      CALL AllocAry( tmpReArray, 6, 'temporary array for InpJoints', ErrStat2, ErrMsg2 )
         if (Failed())  return;

         ! Allocate memory for Joint-related arrays
      ALLOCATE ( InputFileData%Morison%InpJoints(InputFileData%Morison%NJoints), STAT = ErrStat2 )
      IF ( ErrStat2 /= 0 ) THEN
         ErrStat2 = ErrID_Fatal
         ErrMsg2  = 'Error allocating space for InpJoints array.'
         if (Failed())  return;
      END IF

      DO I = 1,InputFileData%Morison%NJoints
            ! read the table entries   JointID   Jointxi     Jointyi    Jointzi      JointAxID   JointOvrlp    in the HydroDyn input file
         call ParseAry( FileInfo_In, CurLine, ' joints table line '//trim( Int2LStr(I)), tmpReArray, size(tmpReArray), ErrStat2, ErrMsg2, UnEc )
            if (Failed())  return;
         InputFileData%Morison%InpJoints(I)%JointID      =  NINT(tmpReArray(1))
         InputFileData%Morison%InpJoints(I)%Position(1)  =       tmpReArray(2)
         InputFileData%Morison%InpJoints(I)%Position(2)  =       tmpReArray(3)
         InputFileData%Morison%InpJoints(I)%Position(3)  =       tmpReArray(4)
         InputFileData%Morison%InpJoints(I)%JointAxID    =  NINT(tmpReArray(5))
         InputFileData%Morison%InpJoints(I)%JointOvrlp   =  NINT(tmpReArray(6))
      END DO

      if (allocated(tmpReArray))      deallocate(tmpReArray)
   END IF


   !-------------------------------------------------------------------------------------------------
   ! Member Cross-section Properties Section
   !-------------------------------------------------------------------------------------------------
   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') trim(FileInfo_In%Lines(CurLine))    ! Write section break to echo
   CurLine = CurLine + 1

      ! NPropSets - Number of member cross-section property sets
   call ParseVar( FileInfo_In, CurLine, 'NPropSets', InputFileData%Morison%NPropSets, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! Table header
   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') 'MPropSets table header line 1: '//NewLine//trim(FileInfo_In%Lines(CurLine))
   CurLine = CurLine + 1
   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') 'MPropSets table header line 2: '//NewLine//trim(FileInfo_In%Lines(CurLine))
   CurLine = CurLine + 1

   IF ( InputFileData%Morison%NPropSets > 0 ) THEN

      CALL AllocAry( tmpReArray, 3, 'temporary array for MPropSets', ErrStat2, ErrMsg2 )
         if (Failed())  return;

         ! Allocate memory for Member cross-section property set-related arrays
      ALLOCATE ( InputFileData%Morison%MPropSets(InputFileData%Morison%NPropSets), STAT = ErrStat2 )
      IF ( ErrStat2 /= 0 ) THEN
         ErrStat2 = ErrID_Fatal
         ErrMsg2  = 'Error allocating space for MPropSets array.'
         if (Failed())  return;
      END IF

      DO I = 1,InputFileData%Morison%NPropSets
         call ParseAry( FileInfo_In, CurLine, ' MPropSets line '//trim( Int2LStr(I)), tmpReArray, size(tmpReArray), ErrStat2, ErrMsg2, UnEc )
            if (Failed())  return;
         InputFileData%Morison%MPropSets(I)%PropSetID = NINT(tmpReArray(1))
         InputFileData%Morison%MPropSets(I)%PropD     =      tmpReArray(2)
         InputFileData%Morison%MPropSets(I)%PropThck  =      tmpReArray(3)
      END DO

      if (allocated(tmpReArray))      deallocate(tmpReArray)
   END IF


   !-------------------------------------------------------------------------------------------------
   ! Simple hydrodynamic coefficients Section
   !-------------------------------------------------------------------------------------------------
   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') trim(FileInfo_In%Lines(CurLine))    ! Write section break to echo
   CurLine = CurLine + 1

      ! Table header
   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') 'Simple hydrodynamic coefficients table header line 1: '//NewLine//trim(FileInfo_In%Lines(CurLine))
   CurLine = CurLine + 1
   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') 'Simple hydrodynamic coefficients table header line 2: '//NewLine//trim(FileInfo_In%Lines(CurLine))
   CurLine = CurLine + 1

   CALL AllocAry( tmpReArray, 12, 'temporary array for Simple hydrodynamic coefficients', ErrStat2, ErrMsg2 )
      if (Failed())  return;
   call ParseAry( FileInfo_In, CurLine, 'Simple hydrodynamic coefficients table row '//trim( Int2LStr(I)), tmpReArray, size(tmpReArray), ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

   InputFileData%Morison%SimplCd       = tmpReArray( 1)
   InputFileData%Morison%SimplCdMG     = tmpReArray( 2)
   InputFileData%Morison%SimplCa       = tmpReArray( 3)
   InputFileData%Morison%SimplCaMG     = tmpReArray( 4)
   InputFileData%Morison%SimplCp       = tmpReArray( 5)
   InputFileData%Morison%SimplCpMG     = tmpReArray( 6)
   InputFileData%Morison%SimplAxCd     = tmpReArray( 7)
   InputFileData%Morison%SimplAxCdMG   = tmpReArray( 8)
   InputFileData%Morison%SimplAxCa     = tmpReArray( 9)
   InputFileData%Morison%SimplAxCaMG   = tmpReArray(10)
   InputFileData%Morison%SimplAxCp     = tmpReArray(11)
   InputFileData%Morison%SimplAxCpMG   = tmpReArray(12)

   if (allocated(tmpReArray))      deallocate(tmpReArray)


   !-------------------------------------------------------------------------------------------------
   ! Depth-based Hydrodynamic Coefficients Section
   !-------------------------------------------------------------------------------------------------
   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') trim(FileInfo_In%Lines(CurLine))    ! Write section break to echo
   CurLine = CurLine + 1

      ! NCoefDpth - Number of depth-based hydrodynamic coefficient property sets
   call ParseVar( FileInfo_In, CurLine, 'NCoefDpth', InputFileData%Morison%NCoefDpth, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! Table header
   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') 'Depth-based hydrodynamic coefficients table header line 1: '//NewLine//trim(FileInfo_In%Lines(CurLine))
   CurLine = CurLine + 1
   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') 'Depth-based hydrodynamic coefficients table header line 2: '//NewLine//trim(FileInfo_In%Lines(CurLine))
   CurLine = CurLine + 1

   IF ( InputFileData%Morison%NCoefDpth > 0 ) THEN

      CALL AllocAry( tmpReArray, 13, 'temporary array for CoefDpths', ErrStat2, ErrMsg2 )
         if (Failed())  return;

         ! Allocate memory for depth-based coefficient arrays
      ALLOCATE ( InputFileData%Morison%CoefDpths(InputFileData%Morison%NCoefDpth), STAT = ErrStat2 )
      IF ( ErrStat2 /= 0 ) THEN
         ErrStat2 = ErrID_Fatal
         ErrMsg2  = 'Error allocating space for CoefDpths array.'
         if (Failed())  return;
      END IF
                  
      DO I = 1,InputFileData%Morison%NCoefDpth
         call ParseAry( FileInfo_In, CurLine, ' CoefDpths coefficients table row '//trim( Int2LStr(I)), tmpReArray, size(tmpReArray), ErrStat2, ErrMsg2, UnEc )
            if (Failed())  return;

         InputFileData%Morison%CoefDpths(I)%Dpth         = tmpReArray( 1)
         InputFileData%Morison%CoefDpths(I)%DpthCd       = tmpReArray( 2)
         InputFileData%Morison%CoefDpths(I)%DpthCdMG     = tmpReArray( 3)
         InputFileData%Morison%CoefDpths(I)%DpthCa       = tmpReArray( 4)
         InputFileData%Morison%CoefDpths(I)%DpthCaMG     = tmpReArray( 5)
         InputFileData%Morison%CoefDpths(I)%DpthCp       = tmpReArray( 6)
         InputFileData%Morison%CoefDpths(I)%DpthCpMG     = tmpReArray( 7)
         InputFileData%Morison%CoefDpths(I)%DpthAxCd     = tmpReArray( 8)
         InputFileData%Morison%CoefDpths(I)%DpthAxCdMG   = tmpReArray( 9)
         InputFileData%Morison%CoefDpths(I)%DpthAxCa     = tmpReArray(10)
         InputFileData%Morison%CoefDpths(I)%DpthAxCaMG   = tmpReArray(11)
         InputFileData%Morison%CoefDpths(I)%DpthAxCp     = tmpReArray(12)
         InputFileData%Morison%CoefDpths(I)%DpthAxCpMG   = tmpReArray(13)
      END DO

      if (allocated(tmpReArray))      deallocate(tmpReArray)
   END IF


   !-------------------------------------------------------------------------------------------------
   ! Member-based Hydrodynamic Coefficients Section
   !-------------------------------------------------------------------------------------------------
   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') trim(FileInfo_In%Lines(CurLine))    ! Write section break to echo
   CurLine = CurLine + 1

      ! NCoefMembers - Number of member-based hydrodynamic coefficient property sets
   call ParseVar( FileInfo_In, CurLine, 'NCoefMembers', InputFileData%Morison%NCoefMembers, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! Table header
   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') 'Member-based hydrodynamic coefficients table header line 1: '//NewLine//trim(FileInfo_In%Lines(CurLine))
   CurLine = CurLine + 1
   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') 'Member-based hydrodynamic coefficients table header line 2: '//NewLine//trim(FileInfo_In%Lines(CurLine))
   CurLine = CurLine + 1

   IF ( InputFileData%Morison%NCoefMembers > 0 ) THEN

      CALL AllocAry( tmpReArray, 25, 'temporary array for CoefMembers', ErrStat2, ErrMsg2 )
         if (Failed())  return;

         ! Allocate memory for Member-based coefficient arrays
      ALLOCATE ( InputFileData%Morison%CoefMembers(InputFileData%Morison%NCoefMembers), STAT = ErrStat2 )
      IF ( ErrStat2 /= 0 ) THEN
         ErrStat2 = ErrID_Fatal
         ErrMsg2  = 'Error allocating space for CoefMembers array.'
         if (Failed())  return;
      END IF

      DO I = 1,InputFileData%Morison%NCoefMembers
         call ParseAry( FileInfo_In, CurLine, 'Member-based hydrodynamic coefficients table row '//trim( Int2LStr(I)), tmpReArray, size(tmpReArray), ErrStat2, ErrMsg2, UnEc )
            if (Failed())  return;

         InputFileData%Morison%CoefMembers(I)%MemberID         = NINT(tmpReArray( 1))
         InputFileData%Morison%CoefMembers(I)%MemberCd1        =      tmpReArray( 2)
         InputFileData%Morison%CoefMembers(I)%MemberCd2        =      tmpReArray( 3)
         InputFileData%Morison%CoefMembers(I)%MemberCdMG1      =      tmpReArray( 4)
         InputFileData%Morison%CoefMembers(I)%MemberCdMG2      =      tmpReArray( 5)
         InputFileData%Morison%CoefMembers(I)%MemberCa1        =      tmpReArray( 6)
         InputFileData%Morison%CoefMembers(I)%MemberCa2        =      tmpReArray( 7)
         InputFileData%Morison%CoefMembers(I)%MemberCaMG1      =      tmpReArray( 8)
         InputFileData%Morison%CoefMembers(I)%MemberCaMG2      =      tmpReArray( 9)
         InputFileData%Morison%CoefMembers(I)%MemberCp1        =      tmpReArray(10)
         InputFileData%Morison%CoefMembers(I)%MemberCp2        =      tmpReArray(11)
         InputFileData%Morison%CoefMembers(I)%MemberCpMG1      =      tmpReArray(12)
         InputFileData%Morison%CoefMembers(I)%MemberCpMG2      =      tmpReArray(13)
         InputFileData%Morison%CoefMembers(I)%MemberAxCd1      =      tmpReArray(14)
         InputFileData%Morison%CoefMembers(I)%MemberAxCd2      =      tmpReArray(15)
         InputFileData%Morison%CoefMembers(I)%MemberAxCdMG1    =      tmpReArray(16)
         InputFileData%Morison%CoefMembers(I)%MemberAxCdMG2    =      tmpReArray(17)
         InputFileData%Morison%CoefMembers(I)%MemberAxCa1      =      tmpReArray(18)
         InputFileData%Morison%CoefMembers(I)%MemberAxCa2      =      tmpReArray(19)
         InputFileData%Morison%CoefMembers(I)%MemberAxCaMG1    =      tmpReArray(20)
         InputFileData%Morison%CoefMembers(I)%MemberAxCaMG2    =      tmpReArray(21)
         InputFileData%Morison%CoefMembers(I)%MemberAxCp1      =      tmpReArray(22)
         InputFileData%Morison%CoefMembers(I)%MemberAxCp2      =      tmpReArray(23)
         InputFileData%Morison%CoefMembers(I)%MemberAxCpMG1    =      tmpReArray(24)
         InputFileData%Morison%CoefMembers(I)%MemberAxCpMG2    =      tmpReArray(25)
      END DO

      if (allocated(tmpReArray))      deallocate(tmpReArray)
   END IF


   !-------------------------------------------------------------------------------------------------
   ! Members Section
   !-------------------------------------------------------------------------------------------------
   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') trim(FileInfo_In%Lines(CurLine))   ! Write section break to echo
   CurLine = CurLine + 1

      ! NMembers - Number of members in the input file
   call ParseVar( FileInfo_In, CurLine, 'NMembers', InputFileData%Morison%NMembers, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! Table header
   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') 'Members table header line 1: '//NewLine//trim(FileInfo_In%Lines(CurLine))
   CurLine = CurLine + 1
   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') 'Members table header line 2: '//NewLine//trim(FileInfo_In%Lines(CurLine))
   CurLine = CurLine + 1

   IF ( InputFileData%Morison%NMembers > 0 ) THEN

         ! Allocate memory for Members arrays
      ALLOCATE ( InputFileData%Morison%InpMembers(InputFileData%Morison%NMembers), STAT = ErrStat2 )
      IF ( ErrStat2 /= 0 ) THEN         
         ErrStat2 = ErrID_Fatal
         ErrMsg2  = 'Error allocating space for InpMembers array.'
         if (Failed())  return;
      END IF

      DO I = 1,InputFileData%Morison%NMembers
         ! We can't use the ParseAry here since PropPot is a logical
         Line = FileInfo_In%Lines(CurLine)
         READ(Line,*,IOSTAT=ErrStat2) InputFileData%Morison%InpMembers(I)%MemberID,   InputFileData%Morison%InpMembers(I)%MJointID1,    &
                                     InputFileData%Morison%InpMembers(I)%MJointID2,   InputFileData%Morison%InpMembers(I)%MPropSetID1,  &
                                     InputFileData%Morison%InpMembers(I)%MPropSetID2, InputFileData%Morison%InpMembers(I)%MDivSize,     &
                                     InputFileData%Morison%InpMembers(I)%MCoefMod,    InputFileData%Morison%InpMembers(I)%PropPot
         IF ( ErrStat2 /= 0 ) THEN
            ErrStat2 = ErrID_Fatal
            ErrMsg2  = 'Error reading members table row '//trim( Int2LStr(I))//', line '  &
                        //trim( Int2LStr(FileInfo_In%FileLine(CurLine)))//' of file '//trim(FileInfo_In%FileList(FileInfo_In%FileIndx(CurLine)))
            if (Failed())  return;
         END IF

         if ( InputFileData%Echo )   WRITE(UnEc, '(A)') trim(FileInfo_In%Lines(CurLine))     ! Echo this line
         CurLine = CurLine+1
      END DO

   END IF


   !-------------------------------------------------------------------------------------------------
   ! Filled Members Section
   !-------------------------------------------------------------------------------------------------
   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') trim(FileInfo_In%Lines(CurLine))   ! Write section break to echo
   CurLine = CurLine + 1

      ! NFillGroups - Number of fill groups
   call ParseVar( FileInfo_In, CurLine, 'NFillGroups', InputFileData%Morison%NFillGroups, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! Table header
   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') 'Fill groups table header line 1: '//NewLine//trim(FileInfo_In%Lines(CurLine))
   CurLine = CurLine + 1
   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') 'Fill groups table header line 2: '//NewLine//trim(FileInfo_In%Lines(CurLine))
   CurLine = CurLine + 1

   IF ( InputFileData%Morison%NFillGroups > 0 ) THEN

         ! Allocate memory for filled group arrays
      ALLOCATE ( InputFileData%Morison%FilledGroups(InputFileData%Morison%NFillGroups), STAT = ErrStat2 )
      IF ( ErrStat2 /= 0 ) THEN
         ErrStat2 = ErrID_Fatal
         ErrMsg2  = 'Error allocating space for FilledGroups array.'
         if (Failed())  return;
      END IF

      DO I = 1,InputFileData%Morison%NFillGroups
         ! We can't use the ParseAry here since the number of entries is indicated by the first entry 
         Line = FileInfo_In%Lines(CurLine)

         READ(Line,*,IOSTAT=ErrStat2) InputFileData%Morison%FilledGroups(I)%FillNumM
         IF ( ErrStat2 /= 0 ) THEN
            ErrStat2 = ErrID_Fatal
            ErrMsg2  = 'Failed to read FillNumM.'
            if (Failed())  return;
         END IF

         ALLOCATE ( InputFileData%Morison%FilledGroups(I)%FillMList(InputFileData%Morison%FilledGroups(I)%FillNumM), STAT = ErrStat2 )
         IF ( ErrStat2 /= 0 ) THEN
            ErrStat2 = ErrID_Fatal
            ErrMsg2  = 'Error allocating space for FillMList array.'
            if (Failed())  return;
         END IF

         READ(Line,*,IOSTAT=ErrStat2) InputFileData%Morison%FilledGroups(I)%FillNumM,  InputFileData%Morison%FilledGroups(I)%FillMList,   &
                                      InputFileData%Morison%FilledGroups(I)%FillFSLoc, InputFileData%Morison%FilledGroups(I)%FillDensChr

         IF ( ErrStat2 /= 0 ) THEN
            ErrStat2 = ErrID_Fatal
            ErrMsg2  = 'Failed to read filled group properties.'
            if (Failed())  return;
         END IF

         if ( InputFileData%Echo )   WRITE(UnEc, '(A)') trim(FileInfo_In%Lines(CurLine))     ! Echo this line
         CurLine = CurLine+1
      END DO

   END IF


   !-------------------------------------------------------------------------------------------------
   ! Marine Growth by Depth Section
   !-------------------------------------------------------------------------------------------------
   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') trim(FileInfo_In%Lines(CurLine))    ! Write section break to echo
   CurLine = CurLine + 1

      ! NMGDepths - Number marine growth depths
   call ParseVar( FileInfo_In, CurLine, 'NMGDepths', InputFileData%Morison%NMGDepths, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! Table header
   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') 'Marine growth by depth table header line 1: '//NewLine//trim(FileInfo_In%Lines(CurLine))
   CurLine = CurLine + 1
   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') 'Marine growth by depth table header line 2: '//NewLine//trim(FileInfo_In%Lines(CurLine))
   CurLine = CurLine + 1

   IF ( InputFileData%Morison%NMGDepths > 0 ) THEN
      CALL AllocAry( tmpReArray, 3, 'temporary array for marine growth table', ErrStat2, ErrMsg2 )
         if (Failed())  return;

         ! Allocate memory for marine growth depths array
      ALLOCATE ( InputFileData%Morison%MGDepths(InputFileData%Morison%NMGDepths), STAT = ErrStat2 )
      IF ( ErrStat2 /= 0 ) THEN
         ErrStat2 = ErrID_Fatal
         ErrMsg2  = 'Error allocating space for MGDepths array.'
         if (Failed())  return;
      END IF

      DO I = 1,InputFileData%Morison%NMGDepths
         call ParseAry( FileInfo_In, CurLine, ' Marine growth table row '//trim( Int2LStr(I)), tmpReArray, size(tmpReArray), ErrStat2, ErrMsg2, UnEc )
         InputFileData%Morison%MGDepths(I)%MGDpth  = tmpReArray(1)
         InputFileData%Morison%MGDepths(I)%MGThck  = tmpReArray(2)
         InputFileData%Morison%MGDepths(I)%MGDens  = tmpReArray(3)
      END DO

      if (allocated(tmpReArray))      deallocate(tmpReArray)
   END IF


   !-------------------------------------------------------------------------------------------------
   ! Member Output List Section
   !-------------------------------------------------------------------------------------------------
   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') trim(FileInfo_In%Lines(CurLine))    ! Write section break to echo
   CurLine = CurLine + 1

      ! NMOutputs - Number of members to output
   call ParseVar( FileInfo_In, CurLine, 'NMOutputs', InputFileData%Morison%NMOutputs, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! Table header
   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') 'Member output list table header line 1: '//NewLine//trim(FileInfo_In%Lines(CurLine))
   CurLine = CurLine + 1
   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') 'Member output list table header line 2: '//NewLine//trim(FileInfo_In%Lines(CurLine))
   CurLine = CurLine + 1

   IF ( InputFileData%Morison%NMOutputs > 0 ) THEN

         ! Allocate memory for filled group arrays
      ALLOCATE ( InputFileData%Morison%MOutLst(InputFileData%Morison%NMOutputs), STAT = ErrStat2 )
      IF ( ErrStat2 /= 0 ) THEN
         ErrStat2 = ErrID_Fatal
         ErrMsg2  = 'Error allocating space for MOutLst array.'
         if (Failed())  return;
      END IF      


      DO I = 1,InputFileData%Morison%NMOutputs

         ! We can't use the ParseAry here since the number of entries is indicated by the first entry 
         Line = FileInfo_In%Lines(CurLine)

         READ(Line,*,IOSTAT=ErrStat2) InputFileData%Morison%MOutLst(I)%MemberID, InputFileData%Morison%MOutLst(I)%NOutLoc
         IF ( ErrStat2 /= 0 ) THEN
            ErrStat2 = ErrID_Fatal
            ErrMsg2  = 'Failed to read NOutLoc.'
            if (Failed())  return;
         END IF      
         
         ALLOCATE ( InputFileData%Morison%MOutLst(I)%NodeLocs(InputFileData%Morison%MOutLst(I)%NOutLoc), STAT = ErrStat2 )
         IF ( ErrStat2 /= 0 ) THEN
            ErrStat2 = ErrID_Fatal
            ErrMsg2  = 'Error allocating space for NodeLocs array.'
            if (Failed())  return;
         END IF      

         ALLOCATE ( InputFileData%Morison%MOutLst(I)%MeshIndx1(InputFileData%Morison%MOutLst(I)%NOutLoc), STAT = ErrStat2 )
         IF ( ErrStat2 /= 0 ) THEN
            ErrStat2 = ErrID_Fatal
            ErrMsg2  = 'Error allocating space for %MeshIndx1 array.'
            if (Failed())  return;
         END IF

         ALLOCATE ( InputFileData%Morison%MOutLst(I)%MemberIndx1(InputFileData%Morison%MOutLst(I)%NOutLoc), STAT = ErrStat2 )
         IF ( ErrStat2 /= 0 ) THEN
            ErrStat2 = ErrID_Fatal
            ErrMsg2  = 'Error allocating space for %MemberIndx1 array.'
            if (Failed())  return;
         END IF

         ALLOCATE ( InputFileData%Morison%MOutLst(I)%MeshIndx2(InputFileData%Morison%MOutLst(I)%NOutLoc), STAT = ErrStat2 )
         IF ( ErrStat2 /= 0 ) THEN
            ErrStat2 = ErrID_Fatal
            ErrMsg2  = 'Error allocating space for %MeshIndx2 array.'
            if (Failed())  return;
         END IF      

         ALLOCATE ( InputFileData%Morison%MOutLst(I)%MemberIndx2(InputFileData%Morison%MOutLst(I)%NOutLoc), STAT = ErrStat2 )
         IF ( ErrStat2 /= 0 ) THEN
            ErrStat2 = ErrID_Fatal
            ErrMsg2  = 'Error allocating space for %MemberIndx2 array.'
            if (Failed())  return;
         END IF    

         ALLOCATE ( InputFileData%Morison%MOutLst(I)%s(InputFileData%Morison%MOutLst(I)%NOutLoc), STAT = ErrStat2 )
         IF ( ErrStat2 /= 0 ) THEN
            ErrStat2 = ErrID_Fatal
            ErrMsg2  = 'Error allocating space for s array.'
            if (Failed())  return;
         END IF      

         READ(Line,*,IOSTAT=ErrStat2) InputFileData%Morison%MOutLst(I)%MemberID,  InputFileData%Morison%MOutLst(I)%NOutLoc,  &
                                      InputFileData%Morison%MOutLst(I)%NodeLocs

         IF ( ErrStat2 /= 0 ) THEN
            ErrStat2 = ErrID_Fatal
            ErrMsg2  = 'Failed to read member output list properties.'
            if (Failed())  return;
         END IF

         if ( InputFileData%Echo )   WRITE(UnEc, '(A)') trim(FileInfo_In%Lines(CurLine))     ! Echo this line
         CurLine = CurLine+1
      END DO

   END IF


   !-------------------------------------------------------------------------------------------------
   ! Joint Output List Section
   !-------------------------------------------------------------------------------------------------
   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') trim(FileInfo_In%Lines(CurLine))    ! Write section break to echo
   CurLine = CurLine + 1

      ! NJOutputs - Number of joints to output
   call ParseVar( FileInfo_In, CurLine, 'NJOutputs', InputFileData%Morison%NJOutputs, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

   IF ( InputFileData%Morison%NJOutputs > 0 ) THEN

      ALLOCATE ( InputFileData%Morison%JOutLst(InputFileData%Morison%NJOutputs), STAT = ErrStat2 )
      IF ( ErrStat2 /= 0 ) THEN
         ErrStat2 = ErrID_Fatal
         ErrMsg2  = 'Error allocating space for JOutLst data structures.'
         if (Failed())  return;
      END IF      

      CALL AllocAry( tmpArray, InputFileData%Morison%NJOutputs, 'temporary array for Joint outputs', ErrStat2, ErrMsg2 )
         if (Failed())  return;

      call ParseAry( FileInfo_In, CurLine, 'JOutLst table row '//trim( Int2LStr(I)), tmpArray, InputFileData%Morison%NJOutputs, ErrStat2, ErrMsg2, UnEc )
         if (Failed())  return;

      DO I = 1,InputFileData%Morison%NJOutputs
         InputFileData%Morison%JOutLst(I)%JointID = tmpArray(I)
      END DO

      DEALLOCATE(tmpArray)   
      
   ELSE
      
      ! There are no Joint Outputs, but there is a line here.  We don't parse it since we don't want to error on it
      if ( InputFileData%Echo )   WRITE(UnEc, '(A)') trim(FileInfo_In%Lines(CurLine))    ! Write line to echo
      CurLine = CurLine + 1

   END IF
   

   !-------------------------------------------------------------------------------------------------
   ! Data section for OUTPUT
   !-------------------------------------------------------------------------------------------------
   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') trim(FileInfo_In%Lines(CurLine))    ! Write section break to echo
   CurLine = CurLine + 1

         ! HDSum - Whether or not to generate a summary file
   call ParseVar( FileInfo_In, CurLine, 'HDSum', InputFileData%HDSum, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

         ! OutAll - Whether or not to output information for every member and joint
   call ParseVar( FileInfo_In, CurLine, 'OutAll', InputFileData%OutAll, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

         ! OutSwtch - Specify how to write to an output file
   call ParseVar( FileInfo_In, CurLine, 'OutSwtch', InputFileData%OutSwtch, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

        ! OutFmt - Format for numerical outputs
   call ParseVar( FileInfo_In, CurLine, 'OutFmt', InputFileData%OutFmt, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

         ! OutSFmt - Format for output column headers
   call ParseVar( FileInfo_In, CurLine, 'OutSFmt', InputFileData%OutSFmt, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;


   !-------------------------------------------------------------------------------------------------
   ! Data section for FLOATING PLATFORM OUTPUTS
   !-------------------------------------------------------------------------------------------------
   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') trim(FileInfo_In%Lines(CurLine))    ! Write section break to echo
   CurLine = CurLine + 1
      
      ! OutList - list of requested parameters to output to a file
   call AllocAry( InputFileData%UserOutputs, MaxUserOutputs, 'InputFileData%UserOutputs', ErrStat2, ErrMsg2 )  ! MaxUserOutputs is set in registry 
      if (Failed())  return;
   
   call ReadOutputListFromFileInfo( FileInfo_In, CurLine, InputFileData%UserOutputs, &
            InputFileData%NUserOutputs, 'OutList', "List of user-requested output channels", ErrStat2, ErrMsg2, UnEc )
         if (Failed()) return;

   
   !-------------------------------------------------------------------------------------------------
   ! This is the end of the input file
   !-------------------------------------------------------------------------------------------------
   
   CALL Cleanup()

   RETURN

CONTAINS
   !..............................
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
      if (Failed)    call Cleanup()
   end function Failed
   SUBROUTINE Cleanup()
      IF (ALLOCATED(tmpArray  )) DEALLOCATE(tmpArray  )
      IF (ALLOCATED(tmpReArray)) DEALLOCATE(tmpReArray)
      IF (ALLOCATED(tmpVec1   )) DEALLOCATE(tmpVec1   )
      IF (ALLOCATED(tmpVec2   )) DEALLOCATE(tmpVec2   )
         ! Cleanup the Echo file and global variables
      if (UnEc > 0)  close ( UnEc )
   END SUBROUTINE Cleanup
END SUBROUTINE HydroDyn_ParseInput



  

!====================================================================================================
SUBROUTINE HydroDynInput_ProcessInitData( InitInp, Interval, InputFileData, ErrStat, ErrMsg )
!     This private subroutine verifies the input required for HydroDyn is correctly specified.
!----------------------------------------------------------------------------------------------------


      ! Passed variables

   TYPE(HydroDyn_InitInputType),  INTENT( IN    )   :: InitInp              ! the hydrodyn data
   REAL(DbKi),                    INTENT( IN    )   :: Interval             ! The DT supplied by the glue code/driver
   TYPE(HydroDyn_InputFile),      INTENT( INOUT )   :: InputFileData        ! the hydrodyn input file data
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
!   CHARACTER(ChanLen), ALLOCATABLE                       :: tmpOutLst(:)         !
   CHARACTER(3)                                     :: TmpExtension         ! Temporary variable for holding the file extension for 10d, 11d, 12d, 10s, 11s, 12s WAMIT files
   LOGICAL                                          :: TmpFileExist         ! Temporary variable in checking the existance of an input file.
   LOGICAL                                          :: JointUsed
   REAL(ReKi)                                       :: l
   REAL(ReKi)                                       :: lvec(3)
   LOGICAL, ALLOCATABLE                             :: foundMask(:)
   INTEGER                                          :: WaveModIn
   
   INTEGER(IntKi)                                   :: ErrStat2, IOS
   CHARACTER(ErrMsgLen)                             :: ErrMsg2
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

   IF ( InputFileData%Waves%WtrDens < 0.0 )  THEN
      CALL SetErrStat( ErrID_Fatal,'WtrDens must not be negative.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF


      ! WtrDpth - Water depth
   
   ! First adjust water depth based on MSL2SWL values
   InputFileData%Morison%WtrDpth = InputFileData%Morison%WtrDpth + InputFileData%Morison%MSL2SWL
   
   IF ( InputFileData%Morison%WtrDpth <= 0.0 )  THEN
      CALL SetErrStat( ErrID_Fatal,'WtrDpth must be greater than zero.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF


      ! MSL2SWL - Mean sea level to still water level

   IF ( InputFileData%PotMod == 1 .AND. .NOT. EqualRealNos(InputFileData%Morison%MSL2SWL, 0.0_ReKi) ) THEN
      CALL SetErrStat( ErrID_Fatal,'MSL2SWL must be 0 when PotMod = 1 (WAMIT).',ErrStat,ErrMsg,RoutineName)        
      RETURN
   END IF
     
   
      ! WaveMod - Wave kinematics model switch.

   IF ( LEN_TRIM(InputFileData%Waves%WaveModChr) > 1 ) THEN

      IF ( InputFileData%Waves%WaveModChr(1:2) == '1P' )  THEN                     ! The user wants to specify the phase in place of a random phase

         READ (InputFileData%Waves%WaveModChr(3:),*,IOSTAT=IOS )  InputFileData%Waves%WavePhase
            CALL CheckIOS ( IOS, "", 'WavePhase', NumType, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg,RoutineName)
            IF ( ErrStat >= AbortErrLev ) RETURN
            
         WaveModIn               = 1
         InputFileData%Waves%WaveMod   = 10                                ! Internally define WaveMod = 10 to mean regular waves with a specified (nonrandom) phase
         InputFileData%Waves%WavePhase = InputFileData%Waves%WavePhase*D2R       ! Convert the phase from degrees to radians

      ELSE                                               ! The user must have specified WaveMod incorrectly.
         CALL SetErrStat( ErrID_Fatal,'WaveMod incorrectly specified',ErrStat,ErrMsg,RoutineName)
         RETURN
      END IF

   ELSE
         ! The line below only works for 1 digit reads
      READ( InputFileData%Waves%WaveModChr, *, IOSTAT=IOS ) InputFileData%Waves%WaveMod
         CALL CheckIOS ( IOS, "", 'WaveMod', NumType, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) RETURN
         
      WaveModIn               = InputFileData%Waves%WaveMod

   END IF ! LEN_TRIM(InputFileData%Waves%WaveModChr)

   IF ( (WaveModIn == 6) .AND. .NOT. EqualRealNos(InputFileData%Morison%MSL2SWL, 0.0_ReKi) ) THEN
      CALL SetErrStat( ErrID_Fatal,'MSL2SWL must be 0 when WaveMod = 6.',ErrStat,ErrMsg,RoutineName)        
      RETURN
   END IF
   

   IF ( WaveModIn < 0 .OR. WaveModIn > 6 ) THEN
      IF ( InputFileData%PotMod == 1  ) THEN
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
   !     (InputFileData%Waves2%WvDiffQTFF /= .false.)                       .or. &
   !     (InputFileData%Waves2%WvSumQTFF /= .false.)                        .or. &
   !     (InputFileData%PotMod /= 0 .or. InputFileData%PotMod /=1)                .or. &
   !     (InputFileData%WAMIT%ExctnMod /=0 .or. InputFileData%WAMIT%ExctnMod /=2) .or. &
   !     (InputFileData%WAMIT%RdtnMod  /=0 .or. InputFileData%WAMIT%RdtnMod  /=2) .or. &
   !     (InputFileData%WAMIT2%MnDrift /=0)                                 .or. &
   !     (InputFileData%WAMIT2%NewmanApp /= 0)                              .or. &
   !     (InputFileData%WAMIT2%SumQTF /= 0 )                                     ) then
   !   
   !end if
        
   
         ! WaveStMod - Model switch for stretching incident wave kinematics to instantaneous free surface.

         ! TODO: We are only implementing WaveStMod = 0 (No stretching) at this point in time. 1 Mar 2013 GJH

   IF ( InputFileData%Waves%WaveStMod /= 0 ) THEN
      CALL SetErrStat( ErrID_Fatal,'WaveStMod must be 0. Future versions of HydroDyn will once again support other wave stretching models.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF

   IF ( InputFileData%Waves%WaveMod /= 6 .AND. InputFileData%Morison%NMembers > 0 .AND. InputFileData%Waves%WaveMod > 0 ) THEN
      
      IF ( ( InputFileData%Waves%WaveStMod /= 0 ) .AND. ( InputFileData%Waves%WaveStMod /= 1 ) .AND. &
            ( InputFileData%Waves%WaveStMod /= 2 ) ) THEN ! (TODO: future version will support 3) .AND. ( InputFileData%Waves%WaveStMod /= 3 ) )  THEN
         ErrMsg  = ' WaveStMod must be 0, 1, or 2.' !, or 3.'
         ErrStat = ErrID_Fatal
   
         RETURN
      END IF
   
      !IF ( ( InputFileData%Waves%WaveStMod /= 3 ) .AND. ( InputFileData%Waves%WaveMod == 5 ) )  THEN
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
   
      InputFileData%Waves%WaveStMod = 0
   
   END IF


      ! WaveTMax - Analysis time for incident wave calculations.

   IF ( InputFileData%Waves%WaveMod == 0 )  THEN   ! .TRUE if we have incident waves.
      
      ! TODO: Issue warning if WaveTMax was not already 0.0 in this case.
      IF ( .NOT. EqualRealNos(InputFileData%Waves%WaveTMax, 0.0_DbKi) ) THEN
         CALL WrScr( '  Setting WaveTMax to 0.0 since WaveMod = 0' )
         InputFileData%Waves%WaveTMax = 0.0
      END IF
      IF ( .NOT. EqualRealNos(InputFileData%Waves%WaveDir, 0.0_SiKi) ) THEN
         CALL WrScr( '  Setting WaveDir to 0.0 since WaveMod = 0' )
         InputFileData%Waves%WaveDir = 0.0
      END IF
   ELSEIF ( InputFileData%Waves%WaveMod == 5 ) THEN   ! User wave elevation file reading in
      IF (InitInp%TMax > InputFileData%Waves%WaveTMax ) THEN
         CALL SetErrstat( ErrID_Fatal, '  WaveTMax must be larger than the simulation time for user wave elevations (WaveMod == 5).',ErrStat,ErrMsg,RoutineName)
         RETURN
      END IF
   ELSE
      IF (InitInp%TMax > InputFileData%Waves%WaveTMax ) THEN
         CALL WrScr( '  WaveTMax is less then the simulation time.  Wave data will repeat every WaveTMax seconds.')
      END IF
   END IF   


      ! WaveDT - Time step for incident wave calculations

   IF ( InputFileData%Waves%WaveMod > 0 )  THEN   ! .TRUE if we have incident waves.

      IF ( InputFileData%Waves%WaveDT <= 0.0 )  THEN
         CALL SetErrStat( ErrID_Fatal,'WaveDT must be greater than zero.',ErrStat,ErrMsg,RoutineName)
         RETURN
      END IF
      
      IF ( (InputFileData%Waves%WaveMod == 6) .AND. (.NOT. EqualRealNos(InputFileData%Waves%WaveDT, Interval)) ) THEN
         CALL SetErrStat( ErrID_Fatal,'WaveDT must equal the simulation DT value when WaveMod = 6.',ErrStat,ErrMsg,RoutineName)
         RETURN
      END IF
   ELSE

      InputFileData%Waves%WaveDT = 0.0

   END IF


       ! WaveHs - Significant wave height

   IF ( ( InputFileData%Waves%WaveMod /= 0 ) .AND. ( InputFileData%Waves%WaveMod /= 4 ) .AND. ( InputFileData%Waves%WaveMod /= 5 ) ) THEN   ! .TRUE. (when WaveMod = 1, 2, 3, or 10) if we have plane progressive (regular), JONSWAP/Pierson-Moskowitz spectrum (irregular) waves, or white-noise waves, but not user-defined or GH Bladed wave data.

      IF ( InputFileData%Waves%WaveHs <= 0.0 )  THEN
         CALL SetErrStat( ErrID_Fatal,'WaveHs must be greater than zero.',ErrStat,ErrMsg,RoutineName)
         RETURN
      END IF

   ELSE

      InputFileData%Waves%WaveHs = 0.0

   END IF


      ! WaveTp - Peak spectral period.
   ! We commented out the if else block due to a bug when WaveMod == 3, and then WaveTp is hence set to 0.0.  See line 1092 of Waves.f90 (as of 11/24/2014) GJH
   !IF ( ( InputFileData%Waves%WaveMod == 1 ) .OR. ( InputFileData%Waves%WaveMod == 2 ) .OR. ( InputFileData%Waves%WaveMod == 10 ) ) THEN   ! .TRUE. (when WaveMod = 1, 2, or 10) if we have plane progressive (regular), JONSWAP/Pierson-Moskowitz spectrum (irregular) waves.

      IF ( InputFileData%Waves%WaveTp <= 0.0 )  THEN
         CALL SetErrStat( ErrID_Fatal,'WaveTp must be greater than zero.',ErrStat,ErrMsg,RoutineName)
         RETURN
      END IF

  ! ELSE

  !    InputFileData%Waves%WaveTp = 0.0

  ! END IF


       ! WavePkShp - Peak shape parameter.

   CALL Conv2UC( InputFileData%Waves%WavePkShpChr )    ! Convert Line to upper case.

   IF ( InputFileData%Waves%WaveMod == 2 ) THEN   ! .TRUE if we have JONSWAP/Pierson-Moskowitz spectrum (irregular) waves, but not GH Bladed wave data.

      IF ( TRIM(InputFileData%Waves%WavePkShpChr) == 'DEFAULT' )  THEN   ! .TRUE. when one wants to use the default value of the peak shape parameter, conditioned on significant wave height and peak spectral period.

         InputFileData%Waves%WavePkShp = WavePkShpDefault ( InputFileData%Waves%WaveHs, InputFileData%Waves%WaveTp )

      ELSE                                   ! The input must have been specified numerically.

         READ (InputFileData%Waves%WavePkShpChr,*,IOSTAT=IOS)  InputFileData%Waves%WavePkShp
            CALL CheckIOS ( IOS, "", 'WavePkShp', NumType, ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2, ErrMsg2,ErrStat,ErrMsg,RoutineName)
            IF ( ErrStat >= AbortErrLev ) RETURN

         IF ( ( InputFileData%Waves%WavePkShp < 1.0 ) .OR. ( InputFileData%Waves%WavePkShp > 7.0 ) )  THEN
            CALL SetErrStat( ErrID_Fatal,'WavePkShp must be greater than or equal to 1 and less than or equal to 7.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF

      END IF

   ELSE

      InputFileData%Waves%WavePkShp = 1.0

   END IF


      ! WvLowCOff and WvHiCOff - Wave Cut-off frequency
    
   IF ( InputFileData%Waves%WvLowCOff < 0 ) THEN
      CALL SetErrStat( ErrID_Fatal,'WvLowCOff must be greater than or equal to zero.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF
   
      ! Threshold upper cut-off based on sampling rate
   IF ( EqualRealNos(InputFileData%Waves%WaveDT, 0.0_DbKi) ) THEN
      InputFileData%Waves%WvHiCOff = 10000.0;  ! This is not going to be used because WaveDT is zero.
   ELSE
      InputFileData%Waves%WvHiCOff =  MIN( REAL( Pi/InputFileData%Waves%WaveDT,SiKi), InputFileData%Waves%WvHiCOff ) 
   END IF
   
   !TODO Issue warning if we changed WvHiCOff  GJH 7/24/13
   
   IF ( InputFileData%Waves%WvLowCOff >= InputFileData%Waves%WvHiCOff ) THEN
      CALL SetErrSTat( ErrID_Fatal,'WvLowCOff must be less than WvHiCOff.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF
   
   
        ! Copy over the first order frequency limits to the WAMIT2 module which needs them.
   InputFileData%WAMIT2%WvLowCOff  = InputFileData%Waves%WvLowCOff
   InputFileData%WAMIT2%WvHiCOff   = InputFileData%Waves%WvHiCOff


      ! WaveDir - Wave heading direction.

   IF ( ( InputFileData%Waves%WaveMod > 0 ) .AND. ( InputFileData%Waves%WaveMod /= 6 ) )  THEN   ! .TRUE if we have incident waves, but not user input wave data.

      IF ( ( InputFileData%Waves%WaveDir <= -180.0 ) .OR. ( InputFileData%Waves%WaveDir > 180.0 ) )  THEN
         CALL SetErrStat( ErrID_Fatal,'WaveDir must be greater than -180 and less than or equal to 180.',ErrStat,ErrMsg,RoutineName)
         RETURN
      END IF

   ELSE

      InputFileData%Waves%WaveDir = 0.0

   END IF


      ! Multi-directional waves

      ! Check the WaveDirMod value
   IF ( InputFileData%Waves%WaveDirMod < 0 .OR. InputFileData%Waves%WaveDirMod > 1 ) THEN
      CALL SetErrStat( ErrID_Fatal,'WaveDirMod must be either 0 (No spreading) or 1 (COS2S spreading function)',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF

      ! Check if we are doing multidirectional waves or not.
      ! We can only use multi directional waves on WaveMod=2,3,4
   InputFileData%Waves%WaveMultiDir = .FALSE.         ! Set flag to false to start
   IF ( InputFileData%Waves%WaveMod >= 2 .AND. InputFileData%Waves%WaveMod <= 4 .AND. InputFileData%Waves%WaveDirMod == 1 ) THEN
      InputFileData%Waves%WaveMultiDir = .TRUE.
   ELSEIF ( (InputFileData%Waves%WaveMod < 2 .OR. InputFileData%Waves%WaveMod >4) .AND. InputFileData%Waves%WaveDirMod == 1 ) THEN
      CALL SetErrStat( ErrID_Warn,'WaveDirMod unused unless WaveMod == 2, 3, or 4.  Ignoring WaveDirMod.',ErrStat,ErrMsg,RoutineName)
   ENDIF


      !  Check to see if the for some reason the wave direction spreading range is set to zero.  If it is, 
      !  we don't have any spreading, so we will turn off the multidirectional waves.
   IF ( InputFileData%Waves%WaveMultiDir .AND. EqualRealNos( InputFileData%Waves%WaveDirRange, 0.0_SiKi ) ) THEN
      CALL SetErrStat( ErrID_Warn,' WaveDirRange set to zero, so multidirectional waves are turned off.',ErrStat,ErrMsg,RoutineName)
      InputFileData%Waves%WaveMultiDir = .FALSE.
   ENDIF



      ! We check the following only if we set WaveMultiDir to true, otherwise ignore them and set them to zero
   IF ( InputFileData%Waves%WaveMultiDir ) THEN

         ! Check WaveDirSpread
      IF ( InputFileData%Waves%WaveDirSpread <= 0.0 ) THEN

         CALL SetErrStat( ErrID_Fatal,'WaveDirSpread cannot negative or zero.',ErrStat,ErrMsg,RoutineName)
         RETURN

      ENDIF


         ! Check that the number of wave directions is a positive odd number.
         !     -> If it is less than 0, error out.
         !     -> If it is even, we will increment it by 1.
      IF ( InputFileData%Waves%WaveNDir <= 0_IntKi ) THEN
         CALL SetErrStat( ErrID_Fatal,' WaveNDir must be an odd number greater than 0.',ErrStat,ErrMsg,RoutineName)
         RETURN
      ENDIF

         ! Check that the value for WaveNDir is odd
      IF ( MODULO( InputFileData%Waves%WaveNDir, 2_IntKi) == 0_IntKi ) THEN
         InputFileData%Waves%WaveNDir  = InputFileData%Waves%WaveNDir + 1
         CALL SetErrStat( ErrID_Warn,'WaveNDir must be odd.  Changing the value to '//Num2LStr(InputFileData%Waves%WaveNDir),ErrStat,ErrMsg,RoutineName)
      ENDIF

         ! Now check that the WaveDirRange is less than 360 degrees (not sure why we would want that)
      IF ( InputFileData%Waves%WaveDirRange > 360.0_ReKi ) THEN
         CALL SetErrStat( ErrID_Fatal,' WaveDirRange should be less than a full circle.',ErrStat,ErrMsg,RoutineName)
      ENDIF

   ELSE  ! Set everything to zero if we aren't going to use it

      InputFileData%Waves%WaveNDir        = 1         ! Only one direction set -- this shouldn't get used later anyhow
      InputFileData%Waves%WaveDirRange    = PiBy2     ! This is so that the constant C=1 in the COS2S function (it shouldn't get called, but in case it does)
      InputFileData%Waves%WaveDirSpread   = 0.0

   END IF


       ! WaveSeed(1), !WaveSeed(2)

   IF ( .NOT. ( ( InputFileData%Waves%WaveMod > 0 ) .AND. ( InputFileData%Waves%WaveMod /= 5 ) .AND. ( InputFileData%Waves%WaveMod /= 10 ) ) ) THEN   !.TRUE. for plane progressive (regular) with random phase or irregular wave 

      DO I = 1,2

         InputFileData%Waves%WaveSeed(I) = 0

      END DO !I

   END IF


      ! WvKinFile

   IF ( InputFileData%Waves%WaveMod == 5 .OR. InputFileData%Waves%WaveMod == 6 ) THEN      ! .TRUE if we are to read user-supplied wave elevation or wave kinematics file(s).

      IF ( LEN_TRIM( InputFileData%Waves%WvKinFile ) == 0 )  THEN
         CALL SetErrStat( ErrID_Fatal,'WvKinFile must not be an empty string.',ErrStat,ErrMsg,RoutineName)
         RETURN
      END IF

      IF ( PathIsRelative( InputFileData%Waves%WvKinFile ) ) THEN
         CALL GetPath( TRIM(InitInp%InputFile), TmpPath )
         InputFileData%Waves%WvKinFile    = TRIM(TmpPath)//TRIM(InputFileData%Waves%WvKinFile)
      END IF
      InputFileData%Waves%WriteWvKin = .FALSE.
   ELSE !don't use this one
      
#ifdef WRITE_WV_KIN
      IF ( LEN_TRIM( InputFileData%Waves%WvKinFile ) == 0 )  THEN
         InputFileData%Waves%WriteWvKin = .FALSE.
      ELSE
         InputFileData%Waves%WriteWvKin = .TRUE.
         IF ( PathIsRelative( InputFileData%Waves%WvKinFile ) ) THEN
            CALL GetPath( TRIM(InputFileData%InputFile), TmpPath )
            InputFileData%Waves%WvKinFile    = TRIM(TmpPath)//TRIM(InputFileData%Waves%WvKinFile)
         END IF
      END IF
      
#else
      InputFileData%Waves%WvKinFile = ""
      InputFileData%Waves%WriteWvKin = .FALSE.
#endif
   END IF


      ! NWaveElev

   IF ( InputFileData%Waves%NWaveElev < 0 ) THEN

      CALL SetErrStat( ErrID_Fatal,'NWaveElev must not be negative.',ErrStat,ErrMsg,RoutineName)
      RETURN

   END IF



      !-------------------------------------------------------------------------
      ! Check 2nd Order Waves section
      !-------------------------------------------------------------------------


      ! Difference frequency cutoffs

      ! WvLowCOffD and WvHiCOffD - Wave Cut-off frequency
   IF ( InputFileData%Waves2%WvLowCOffD < 0 ) THEN
      CALL SetErrStat( ErrID_Fatal,'WvLowCOffD must be greater than or equal to zero.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF

      ! Check that the order given makes sense. 
   IF ( InputFileData%Waves2%WvLowCOffD >= InputFileData%Waves2%WvHiCOffD ) THEN
      CALL SetErrStat( ErrID_Fatal,'WvLowCOffD must be less than WvHiCOffD.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF
   

      ! Sum frequency cutoffs

      ! WvLowCOffS and WvHiCOffD - Wave Cut-off frequency
   IF ( InputFileData%Waves2%WvLowCOffS < 0 ) THEN
      CALL SetErrStat( ErrID_Fatal,'WvLowCOffS must be greater than or equal to zero.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF

      ! Check that the order given makes sense. 
   IF ( InputFileData%Waves2%WvLowCOffS >= InputFileData%Waves2%WvHiCOffS ) THEN
      CALL SetErrStat( ErrID_Fatal,'WvLowCOffS must be less than WvHiCOffS.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF


        ! Copy over the 2nd order limits to the WAMIT2 module which needs them.
   InputFileData%WAMIT2%WvLowCOffD  = InputFileData%Waves2%WvLowCOffD
   InputFileData%WAMIT2%WvHiCOffD   = InputFileData%Waves2%WvHiCOffD
   InputFileData%WAMIT2%WvLowCOffS  = InputFileData%Waves2%WvLowCOffS
   InputFileData%WAMIT2%WvHiCOffS   = InputFileData%Waves2%WvHiCOffS



      !-------------------------------------------------------------------------
      ! Check Current section
      !-------------------------------------------------------------------------
      

      ! CurrMod - Current profile model switch

   IF ( ( InputFileData%Current%CurrMod /= 0 ) .AND. ( InputFileData%Current%CurrMod /= 1 ) .AND. ( InputFileData%Current%CurrMod /= 2 ) )  THEN
      CALL SetErrStat( ErrID_Fatal,'CurrMod must be 0, 1, or 2.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF

   IF ( ( InputFileData%Current%CurrMod /= 0 ) .AND. ( InputFileData%Waves%WaveMod == 6 ) )  THEN
      CALL SetErrStat( ErrID_Fatal,'CurrMod must be set to 0 when WaveMod is set to 6: user-input wave data.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF


      ! CurrSSV0 - Sub-surface current velocity at still water level

   IF ( InputFileData%Current%CurrMod == 1 )  THEN  ! .TRUE if we have standard current.

      IF ( InputFileData%Current%CurrSSV0 < 0.0 )  THEN
         CALL SetErrStat( ErrID_Fatal,'CurrSSV0 must not be less than zero.',ErrStat,ErrMsg,RoutineName)
         RETURN
      END IF

   ELSE

      InputFileData%Current%CurrSSV0 = 0.0

   END IF


      ! CurrSSDirChr - Sub-surface current heading direction

   IF ( InputFileData%Current%CurrMod == 1 )  THEN  ! .TRUE if we have standard current.


      IF ( TRIM(InputFileData%Current%CurrSSDirChr) == 'DEFAULT' )  THEN   ! .TRUE. when one wants to use the default value of codirectionality between sub-surface current and incident wave propogation heading directions.

         IF ( InputFileData%Waves%WaveMod == 0 ) THEN
            CALL SetErrStat( ErrID_Fatal,'CurrSSDir must not be set to ''DEFAULT'' when WaveMod is set to 0.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF

         InputFileData%Current%CurrSSDir = InputFileData%Waves%WaveDir

      ELSE                                   ! The input must have been specified numerically.

         READ (InputFileData%Current%CurrSSDirChr,*,IOSTAT=IOS)  InputFileData%Current%CurrSSDir
            CALL CheckIOS ( IOS, "", 'CurrSSDir', NumType, ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2, ErrMsg2,ErrStat,ErrMsg,RoutineName)
            IF ( ErrStat >= AbortErrLev ) RETURN

         IF ( ( InputFileData%Current%CurrSSDir <= -180.0 ) .OR. ( InputFileData%Current%CurrSSDir > 180.0 ) )  THEN
            CALL SetErrStat( ErrID_Fatal,'CurrSSDir must be greater than -180 and less than or equal to 180.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF

      END IF


   ELSE

      InputFileData%Current%CurrSSDir = 0.0

   END IF


      ! CurrNSRef - Near-surface current reference depth.

   IF ( InputFileData%Current%CurrMod == 1 )  THEN  ! .TRUE if we have standard current.

      IF ( InputFileData%Current%CurrNSRef <= 0.0 ) THEN
         CALL SetErrStat( ErrID_Fatal,'CurrNSRef must be greater than zero.',ErrStat,ErrMsg,RoutineName)
         RETURN
      END IF

   ELSE

      InputFileData%Current%CurrNSRef = 0.0

   END IF



        ! CurrNSV0 - Near-surface current velocity at still water level.

   IF ( InputFileData%Current%CurrMod == 1 )  THEN  ! .TRUE if we have standard current.

      IF ( InputFileData%Current%CurrNSV0 < 0.0 ) THEN
         CALL SetErrStat( ErrID_Fatal,'CurrNSV0 must not be less than zero.',ErrStat,ErrMsg,RoutineName)
         RETURN
      END IF

   ELSE

      InputFileData%Current%CurrNSV0 = 0.0

   END IF


      ! CurrNSDir - Near-surface current heading direction.

   IF ( InputFileData%Current%CurrMod == 1 )  THEN  ! .TRUE if we have standard current.

      IF ( ( InputFileData%Current%CurrNSDir <= -180.0 ) .OR. ( InputFileData%Current%CurrNSDir > 180.0 ) )  THEN
         CALL SetErrStat( ErrID_Fatal,'CurrNSDir must be greater than -180 and less than or equal to 180.',ErrStat,ErrMsg,RoutineName)
         RETURN
      END IF

   ELSE

      InputFileData%Current%CurrNSDir = 0.0

   END IF


      ! CurrDIV - Depth-independent current velocity.

   IF ( InputFileData%Current%CurrMod == 1 )  THEN  ! .TRUE if we have standard current.

      IF ( InputFileData%Current%CurrDIV < 0.0 ) THEN
         CALL SetErrStat( ErrID_Fatal,'CurrDIV must not be less than zero.',ErrStat,ErrMsg,RoutineName)
         RETURN
      END IF

   ELSE

      InputFileData%Current%CurrDIV = 0.0

   END IF


      ! CurrDIDir - Depth-independent current heading direction.

   IF ( InputFileData%Current%CurrMod == 1 )  THEN  ! .TRUE if we have standard current.

      IF ( ( InputFileData%Current%CurrDIDir <= -180.0 ) .OR. ( InputFileData%Current%CurrDIDir > 180.0 ) ) THEN
         CALL SetErrStat( ErrID_Fatal,'CurrDIDir must be greater than -180 and less than or equal to 180.',ErrStat,ErrMsg,RoutineName)
         RETURN
      END IF

   ELSE

      InputFileData%Current%CurrDIDir = 0.0

   END IF

       ! PotFile - Root name of potential flow files

   IF ( InputFileData%PotMod > 0 ) THEN
      do i = 1,InputFileData%nWAMITObj
          IF ( LEN_TRIM( InputFileData%PotFile(i) ) == 0 ) THEN
            CALL SetErrStat( ErrID_Fatal,'PotFile must not be an empty string.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF

            ! if this is a relative path, let's make it relative to the location of the main input file
            ! tell the WAMIT and WAMIT2 modules what the filename is

         IF ( PathIsRelative( InputFileData%PotFile(i) ) ) THEN
            CALL GetPath( TRIM(InitInp%InputFile), TmpPath )
            InputFileData%PotFile(i)            = TRIM(TmpPath)//TRIM(InputFileData%PotFile(i))
         END IF
      end do

   !TODO: Move this to where the WAMIT modules are initialized
         InputFileData%WAMIT%WAMITFile    = InputFileData%PotFile(1)
         InputFileData%WAMIT2%WAMITFile   = InputFileData%PotFile(1)
      
         ! Set the flag for multidirectional waves for WAMIT2 module.  It needs to know since the Newman approximation
         ! can only use uni-directional waves.
         InputFileData%WAMIT2%WaveMultiDir = InputFileData%Waves%WaveMultiDir

   ELSE
      InputFileData%PotFile            = ""
      InputFileData%WAMIT%WAMITFile    = ""
      InputFileData%WAMIT2%WAMITFile   = ""     
   END IF

      ! Set the WAMIT file name on the Convolution module
   InputFileData%WAMIT%Conv_Rdtn%WAMITFile = InputFileData%WAMIT%WAMITFile

      ! WAMITULEN - WAMIT characteristic body length scale

   IF ( InputFileData%PotMod == 1 ) THEN
!TODO: Deal with WAMIT2 and check each WAMITULEN not just the first
      InputFileData%WAMIT2%WAMITULEN = InputFileData%WAMITULEN(1)    ! Copy to the WAMIT2 module info
      do i = 1,InputFileData%nWAMITObj
         IF ( InputFileData%WAMITULEN(i) < 0.0 ) THEN
            CALL SetErrStat( ErrID_Fatal,'WAMITULEN must be positive.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF
      end do
   ELSE

      InputFileData%WAMITULEN = 1.0
      InputFileData%WAMIT2%WAMITULEN = 1.0

   END IF


      ! PtfmVol0 - Displaced volume of water when the platform is in its undisplaced position

   IF ( InputFileData%PotMod == 1 ) THEN
      do i = 1,InputFileData%nWAMITObj
         IF ( InputFileData%PtfmVol0(i) < 0.0 ) THEN
            CALL SetErrStat( ErrID_Fatal,'PtfmVol0 must not be negative.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF
      end do
   ELSE

      InputFileData%PtfmVol0 = 0.0

   END IF

      ! PtfmRefzt - The zt offset of the body reference point(s) from (0,0,0) (meters) [1 to NBody] 
      ! NOTE: only used when PotMod=1. If NBodyMod=2,PtfmRefzt=0.0

   IF ( InputFileData%PotMod == 1 .and. InputFileData%NBodyMod == 2) THEN
      do i = 1,InputFileData%NBody
         IF ( .not. EqualRealNos( InputFileData%PtfmRefzt(i), 0.0_ReKi ) )THEN
            CALL SetErrStat( ErrID_Fatal,'PtfmRefzt must be 0.0 for all WAMIT bodies when NBodyMod=2.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF
      end do
   END IF

      ! RdtnTMax - Analysis time for wave radiation kernel calculations
      ! NOTE: Use RdtnTMax = 0.0 to eliminate wave radiation damping

   IF ( InputFileData%PotMod == 1 ) THEN

      IF ( InputFileData%WAMIT%RdtnTMax < 0.0 ) THEN
         CALL SetErrStat( ErrID_Fatal,'RdtnTMax must not be negative.',ErrStat,ErrMsg,RoutineName)
         RETURN
      END IF

   ELSE

      InputFileData%WAMIT%RdtnTMax = 0.0

   END IF

       ! RdtnDT - Time step for wave radiation kernel calculations

   IF ( InputFileData%PotMod == 1 ) THEN

      CALL Conv2UC( InputFileData%WAMIT%Conv_Rdtn%RdtnDTChr )    ! Convert Line to upper case.
      
      IF ( TRIM(InputFileData%WAMIT%Conv_Rdtn%RdtnDTChr) == 'DEFAULT' )  THEN   ! .TRUE. when one wants to use the default value timestep provided by the glue code.

         InputFileData%WAMIT%Conv_Rdtn%RdtnDT = Interval 

      ELSE                                   ! The input must have been specified numerically.

         READ (InputFileData%WAMIT%Conv_Rdtn%RdtnDTChr,*,IOSTAT=IOS)  InputFileData%WAMIT%Conv_Rdtn%RdtnDT
            CALL CheckIOS ( IOS, "", 'RdtnDT', NumType, ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2, ErrMsg2,ErrStat,ErrMsg,RoutineName)
            IF ( ErrStat >= AbortErrLev ) RETURN
         
      END IF

      IF ( InputFileData%WAMIT%Conv_Rdtn%RdtnDT <= 0.0 ) THEN
         CALL SetErrStat( ErrID_Fatal,'RdtnDT must be greater than zero.',ErrStat,ErrMsg,RoutineName)
         RETURN
      END IF

      if ( (.not. ( EqualRealNos(Interval, InputFileData%WAMIT%Conv_Rdtn%RdtnDT) ) ) .and. ( (InputFileData%WAMIT%ExctnMod > 1) .or. (InputFileData%WAMIT%RdtnMod > 0) ) ) then
         call SetErrStat( ErrID_Fatal,'RdtnDT must be equal to the glue-code DT if PotMod = 1 and using RdtnMod > 0 or ExctnMod > 1.',ErrStat,ErrMsg,RoutineName)
         return
      end if
      
   ELSE

      InputFileData%WAMIT%Conv_Rdtn%RdtnDT = 0.0

   END IF


   !-------------------------------------------------------------------------------------------------
   ! Second order Forces due to Waves section (WAMIT2 Module)
   !-------------------------------------------------------------------------------------------------
   
      ! Check that we only specified one of MnDrift, NewmanApp, or DiffQTF
      !        (compared pairwise -- if any two are both true, we have a problem)
   IF ( ( InputFileData%WAMIT2%MnDrift /= 0 .AND. InputFileData%WAMIT2%NewmanApp /= 0 ) .OR. &
        ( InputFileData%WAMIT2%DiffQTF /= 0 .AND. InputFileData%WAMIT2%NewmanApp /= 0 ) .OR. &
        ( InputFileData%WAMIT2%MnDrift /= 0 .AND. InputFileData%WAMIT2%DiffQTF   /= 0 ) ) THEN
      CALL SetErrStat( ErrID_Fatal,'Only one of MnDrift, NewmanApp, or DiffQTF can be non-zero.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF


   if ( InputFileData%NBody > 1 .and. InputFileData%WAMIT2%MnDrift == 8 ) then
      call SetErrStat( ErrID_Fatal,'MnDrift cannot equal 8 when NBody > 1.',ErrStat,ErrMsg,RoutineName)
      return
   end if

   if ( InputFileData%NBody > 1 .and. InputFileData%WAMIT2%NewmanApp == 8 ) then
      call SetErrStat( ErrID_Fatal,'NewmanApp cannot equal 8 when NBody > 1.',ErrStat,ErrMsg,RoutineName)
      return
   end if
   

      ! Check MnDrift and set the flag indicating WAMIT2 should perform the mean drift calculation.
      ! Also make sure we have a valid input value for the file extension

   IF ( InputFileData%WAMIT2%MnDrift == 0 ) THEN      ! not using MnDrift
      InputFileData%WAMIT2%MnDriftF = .FALSE.
   ELSE IF ( InputFileData%WAMIT2%MnDrift == 7  .OR. InputFileData%WAMIT2%MnDrift == 8  .OR. InputFileData%WAMIT2%MnDrift == 9 .OR. &
             InputFileData%WAMIT2%MnDrift == 10 .OR. InputFileData%WAMIT2%MnDrift == 11 .OR. InputFileData%WAMIT2%MnDrift == 12 ) THEN   ! Valid values for MnDrift
      IF ( InputFileData%PotMod /= 1 ) THEN
         CALL SetErrStat( ErrID_warn,'MnDrift can only be used with PotMod==1.  Turning off',ErrStat,ErrMsg,RoutineName)
         InputFileData%WAMIT2%MnDriftF = .FALSE.
      ELSE
         InputFileData%WAMIT2%MnDriftF = .TRUE.
      ENDIF
   ELSE     ! Must have received an invalid value
      CALL SetErrStat( ErrID_Fatal,'MnDrift can only have values of 0, 7, 8, 9, 10, 11, or 12.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF


      ! Check NewmanApp and set the flag indicating WAMIT2 should perform the mean drift calculation.
      ! Also make sure we have a valid input value for the file extension

   IF ( InputFileData%WAMIT2%NewmanApp == 0 ) THEN    ! not using NewmanApp
      InputFileData%WAMIT2%NewmanAppF = .FALSE.
   ELSE IF ( InputFileData%WAMIT2%NewmanApp == 7  .OR. InputFileData%WAMIT2%NewmanApp == 8  .OR. InputFileData%WAMIT2%NewmanApp == 9 .OR. &
             InputFileData%WAMIT2%NewmanApp == 10 .OR. InputFileData%WAMIT2%NewmanApp == 11 .OR. InputFileData%WAMIT2%NewmanApp == 12 ) THEN ! Valid values for NewmanApp
      IF ( InputFileData%PotMod /= 1 ) THEN
         CALL SetErrStat( ErrID_warn,'NewmanApp can only be used with PotMod==1.  Turning off',ErrStat,ErrMsg,RoutineName)
         InputFileData%WAMIT2%NewmanAppF = .FALSE.
      ELSE
         InputFileData%WAMIT2%NewmanAppF = .TRUE.
      ENDIF
   ELSE     ! Must have received an invalid value
      CALL SetErrStat( ErrID_Fatal,'NewmanApp can only have values of 0, 7, 8, 9, 10, 11, or 12.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF


      ! Check DiffQTF and set the flag indicating WAMIT2 should perform the mean drift calculation.
      ! Also make sure we have a valid input value for the file extension

   IF ( InputFileData%WAMIT2%DiffQTF == 0 ) THEN      ! not using DiffQTF method
       InputFileData%WAMIT2%DiffQTFF = .FALSE.
   ELSE IF ( InputFileData%WAMIT2%DiffQTF == 10 .OR. InputFileData%WAMIT2%DiffQTF == 11 .OR. InputFileData%WAMIT2%DiffQTF == 12 ) THEN    ! Valid values for DiffQTF
      IF ( InputFileData%PotMod /= 1 ) THEN
         CALL SetErrStat( ErrID_warn,'DiffQTF can only be used with PotMod==1.  Turning off',ErrStat,ErrMsg,RoutineName)
         InputFileData%WAMIT2%DiffQTFF = .FALSE.
      ELSE
         InputFileData%WAMIT2%DiffQTFF = .TRUE.
      ENDIF
   ELSE
      CALL SetErrStat( ErrID_Fatal,'DiffQTF can only have values of 0, 10, 11, or 12.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF


      ! Check SumQTF and set the flag indicating WAMIT2 should perform the mean drift calculation.
      ! Also make sure we have a valid input value for the file extension

   IF ( InputFileData%WAMIT2%SumQTF == 0 ) THEN       ! not using SumQTF method
      InputFileData%WAMIT2%SumQTFF = .FALSE.
   ELSE IF ( InputFileData%WAMIT2%SumQTF == 10 .OR. InputFileData%WAMIT2%SumQTF == 11 .OR. InputFileData%WAMIT2%SumQTF == 12 ) THEN       ! Valid values for SumQTF
      IF ( InputFileData%PotMod /= 1 ) THEN
         CALL SetErrStat( ErrID_warn,'SumQTF can only be used with PotMod==1.  Turning off',ErrStat,ErrMsg,RoutineName)
         InputFileData%WAMIT2%SumQTFF = .FALSE.
      ELSE
         InputFileData%WAMIT2%SumQTFF = .TRUE.
      ENDIF
   ELSE
      CALL SetErrStat( ErrID_Fatal,'SumQTF can only have values of 0, 10, 11, or 12.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF


      ! Check that the min / max diff frequencies make sense if using any DiffQTF method
   IF ( InputFileData%WAMIT2%DiffQTF /= 0 .OR. InputFileData%WAMIT2%MnDrift /= 0 .OR. InputFileData%WAMIT2%NewmanApp /=0 ) THEN
      IF ( ( InputFileData%WAMIT2%WvHiCOffD < InputFileData%WAMIT2%WvLowCOffD ) .OR. ( InputFileData%WAMIT2%WvLowCOffD < 0.0 ) ) THEN
         CALL SetErrStat( ErrID_Fatal,'WvHiCOffD must be larger than WvLowCOffD. Both must be positive.',ErrStat,ErrMsg,RoutineName)
         RETURN
      END IF
   ELSE  ! set to zero since we don't need them
      InputFileData%WAMIT2%WvLowCOffD  = 0.0
      InputFileData%WAMIT2%WvHiCOffD  = 0.0
   END IF


      ! Check that the min / max diff frequencies make sense if using SumQTF
   IF ( InputFileData%WAMIT2%SumQTF /= 0 ) THEN
      IF ( ( InputFileData%WAMIT2%WvHiCOffS < InputFileData%WAMIT2%WvLowCOffS ) .OR. ( InputFileData%WAMIT2%WvLowCOffS < 0.0 ) ) THEN
         CALL SetErrStat( ErrID_Fatal,'WvHiCOffS must be larger than WvLowCOffS. Both must be positive.',ErrStat,ErrMsg,RoutineName)
         RETURN
      END IF
   ELSE  ! set to zero since we don't need them
      InputFileData%WAMIT2%WvLowCOffS  = 0.0
      InputFileData%WAMIT2%WvHiCOffS  = 0.0
   END IF


      ! now that it has been established that the input parameters for second order are good, we check to make sure that the WAMIT files actually exist.
      ! Check MnDrift file
   IF ( InputFileData%WAMIT2%MnDriftF ) THEN
      ! Check if using QTF file types (10d, 11d, 12d) or not (7,8,9)
      IF ( InputFileData%WAMIT2%MnDrift <= 9 ) THEN
         TmpExtension = TRIM(Num2LStr(InputFileData%WAMIT2%MnDrift))
         INQUIRE( file=TRIM(InputFileData%WAMIT2%WAMITFile)//'.'//TRIM(TmpExtension), exist=TmpFileExist )
      ELSE  ! 10, 11, 12
         TmpExtension = TRIM(Num2LStr(InputFileData%WAMIT2%MnDrift))//'d'
         INQUIRE( file=TRIM(InputFileData%WAMIT2%WAMITFile)//'.'//TRIM(TmpExtension), exist=TmpFileExist )
      ENDIF
      IF ( .not. TmpFileExist ) THEN
         CALL SetErrStat( ErrID_Fatal,'Cannot find the WAMIT file '//TRIM(InputFileData%WAMIT2%WAMITFile)//'.'//TRIM(TmpExtension)// &
                    ' required by the MnDrift option.',ErrStat,ErrMsg,RoutineName)
         RETURN
      END IF
   END IF

      ! Check existence of NewmanApp file
   IF ( InputFileData%WAMIT2%NewmanAppF ) THEN
      ! Check if using QTF file types (10d, 11d, 12d) or not (7,8,9)
      IF ( InputFileData%WAMIT2%NewmanApp <= 9 ) THEN
         TmpExtension = TRIM(Num2LStr(InputFileData%WAMIT2%NewmanApp))
         INQUIRE( file=TRIM(InputFileData%WAMIT2%WAMITFile)//'.'//TRIM(TmpExtension), exist=TmpFileExist )
      ELSE  ! 10, 11, 12
         TmpExtension = TRIM(Num2LStr(InputFileData%WAMIT2%NewmanApp))//'d'
         INQUIRE( file=TRIM(InputFileData%WAMIT2%WAMITFile)//'.'//TRIM(TmpExtension), exist=TmpFileExist )
      ENDIF
      IF ( .not. TmpFileExist ) THEN
         CALL SetErrStat( ErrID_Fatal,'Cannot find the WAMIT file '//TRIM(InputFileData%WAMIT2%WAMITFile)//'.'//TRIM(TmpExtension)// &
                    ' required by the NewmanApp option.',ErrStat,ErrMsg,RoutineName)
         RETURN
      END IF
   END IF

   IF ( InputFileData%WAMIT2%DiffQTFF ) THEN
      TmpExtension = TRIM(Num2LStr(InputFileData%WAMIT2%DiffQTF))//'d'
      INQUIRE( file=TRIM(InputFileData%WAMIT2%WAMITFile)//'.'//TRIM(TmpExtension), exist=TmpFileExist )
      IF ( .not. TmpFileExist ) THEN
         CALL SetErrStat( ErrID_Fatal,'Cannot find the WAMIT file '//TRIM(InputFileData%WAMIT2%WAMITFile)//'.'//TRIM(TmpExtension)// &
                    ' required by the DiffQTF option.',ErrStat,ErrMsg,RoutineName)
         RETURN
      END IF
   END IF

   IF ( InputFileData%WAMIT2%SumQTFF ) THEN
      TmpExtension = TRIM(Num2LStr(InputFileData%WAMIT2%SumQTF))//'s'
      INQUIRE( file=TRIM(InputFileData%WAMIT2%WAMITFile)//'.'//TRIM(TmpExtension), exist=TmpFileExist )
      IF ( .not. TmpFileExist ) THEN
         CALL SetErrStat( ErrID_Fatal,'Cannot find the WAMIT file '//TRIM(InputFileData%WAMIT2%WAMITFile)//'.'//TRIM(TmpExtension)// &
                    ' required by the SumQTF option.',ErrStat,ErrMsg,RoutineName)
         RETURN
      END IF
   END IF

   !..................
   ! check for ExctnMod = 2 requirements
   !..................
   if ( (InputFileData%WAMIT%ExctnMod == 2) ) then

      if ( InputFileData%Waves%WaveMod == 6 ) then
         call SetErrStat( ErrID_Fatal, 'Externally generated full wave-kinematics time series cannot be used with state-space wave excitations. Set WaveMod 0, 1, 1P#, 2, 3, 4, or 5.', ErrStat, ErrMsg, RoutineName )
      end if
      
      if ( InputFileData%Waves%WaveDirMod /= 0 ) then
         call SetErrStat( ErrID_Fatal, 'Directional spreading cannot be used with state-space wave excitations. Set WaveDirMod=0.', ErrStat, ErrMsg, RoutineName )
      end if
      
      if ( InputFileData%Waves2%WvDiffQTFF ) then
         call SetErrStat( ErrID_Fatal, 'Cannot use full difference-frequency 2nd-order wave kinematics with state-space wave excitations. Set WvDiffQTF=FALSE.', ErrStat, ErrMsg, RoutineName )
      end if
      
      if ( InputFileData%Waves2%WvSumQTFF ) then
         call SetErrStat( ErrID_Fatal, 'Cannot use full summation-frequency 2nd-order wave kinematics with state-space wave excitations. Set WvSumQTF=FALSE.', ErrStat, ErrMsg, RoutineName )
      end if

      if ( InputFileData%PotMod /= 1 ) then
         call SetErrStat( ErrID_Fatal, 'Potential-flow model via WAMIT must be used with state-space wave excitations. Set PotMod= 1.', ErrStat, ErrMsg, RoutineName )
      end if
      
      if ( InputFileData%WAMIT2%MnDrift /= 0 ) then
         call SetErrStat( ErrID_Fatal, 'Mean-drift 2nd-order forces cannot be used with state-space wave excitations. Set MnDrift=0.', ErrStat, ErrMsg, RoutineName )
      end if

      if ( InputFileData%WAMIT2%NewmanApp /= 0 ) then
         call SetErrStat( ErrID_Fatal, "Mean- and slow-drift 2nd-order forces computed with Newman's approximation cannot be used with state-space wave excitations. Set NewmanApp=0.", ErrStat, ErrMsg, RoutineName )
      end if
      
      if ( InputFileData%WAMIT2%DiffQTF /= 0 ) then
         call SetErrStat( ErrID_Fatal, 'Full difference-frequency 2nd-order forces computed with full QTF cannot be used with state-space wave excitations. Set DiffQTF=0.', ErrStat, ErrMsg, RoutineName )
      end if

      if ( InputFileData%WAMIT2%SumQTF /= 0 ) then
         call SetErrStat( ErrID_Fatal, 'Full summation-frequency 2nd-order forces computed with full QTF cannot be used with State-space wave excitations. Set SumQTF=0.', ErrStat, ErrMsg, RoutineName )
      end if

   end if

   !..................
   ! check for linearization
   !..................
   if (InitInp%Linearize) then
      
      if ( InputFileData%Waves%WaveMod /= 0 ) then
         call SetErrStat( ErrID_Fatal, 'Still water conditions must be used for linearization. Set WaveMod=0.', ErrStat, ErrMsg, RoutineName )
      end if
      
      if ( InputFileData%Waves%WaveDirMod /= 0 ) then
         call SetErrStat( ErrID_Fatal, 'No directional spreading must be used for linearization. Set WaveDirMod=0.', ErrStat, ErrMsg, RoutineName )
      end if
      
      if ( InputFileData%Waves2%WvDiffQTFF ) then
         call SetErrStat( ErrID_Fatal, 'Cannot use full difference-frequency 2nd-order wave kinematics for linearization. Set WvDiffQTF=FALSE.', ErrStat, ErrMsg, RoutineName )
      end if
      
      if ( InputFileData%Waves2%WvSumQTFF ) then
         call SetErrStat( ErrID_Fatal, 'Cannot use full summation-frequency 2nd-order wave kinematics for linearization. Set WvSumQTF=FALSE.', ErrStat, ErrMsg, RoutineName )
      end if

      if ( InputFileData%PotMod > 1 ) then
         call SetErrStat( ErrID_Fatal, 'Potential-flow model cannot be set to FIT for linearization. Set PotMod= 0 or 1.', ErrStat, ErrMsg, RoutineName )
      end if
      
      if ( (InputFileData%WAMIT%ExctnMod == 1) ) then
         call SetErrStat( ErrID_Fatal, 'Cannot set wave excitation model to DFT for linearization. Set ExctnMod=0 or 2.', ErrStat, ErrMsg, RoutineName )
      end if

      if ( InputFileData%WAMIT%RdtnMod == 1 ) then
         call SetErrStat( ErrID_Fatal, 'Cannot set wave radiation model to convolution for linearization. Set RdtnMod=0 or 2.', ErrStat, ErrMsg, RoutineName )
      end if
      
      if ( InputFileData%WAMIT2%MnDrift /= 0 ) then
         call SetErrStat( ErrID_Fatal, 'Mean-drift 2nd-order forces cannot be used for linearization. Set MnDrift=0.', ErrStat, ErrMsg, RoutineName )
      end if

      if ( InputFileData%WAMIT2%NewmanApp /= 0 ) then
         call SetErrStat( ErrID_Fatal, "Mean- and slow-drift 2nd-order forces computed with Newman's approximation cannot be used for linearization. Set NewmanApp=0.", ErrStat, ErrMsg, RoutineName )
      end if
      
      if ( InputFileData%WAMIT2%DiffQTF /= 0 ) then
         call SetErrStat( ErrID_Fatal, 'Full difference-frequency 2nd-order forces computed with full QTF cannot be used for linearization. Set DiffQTF=0.', ErrStat, ErrMsg, RoutineName )
      end if

      if ( InputFileData%WAMIT2%SumQTF /= 0 ) then
         call SetErrStat( ErrID_Fatal, 'Full summation-frequency 2nd-order forces computed with full QTF cannot be used for linearization. Set SumQTF=0.', ErrStat, ErrMsg, RoutineName )
      end if

   end if
 




   !-------------------------------------------------------------------------------------------------
   ! Member Joints Section
   !-------------------------------------------------------------------------------------------------

   IF ( InputFileData%Morison%NJoints < 0 ) THEN
      CALL SetErrStat( ErrID_Fatal,'NJoints parameter cannot be negative.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF

   IF ( InputFileData%Morison%NJoints == 1 ) THEN
      CALL SetErrStat( ErrID_Fatal,'NJoints parameter cannot be set to 1.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF

     
     
      ! Check the axial coefs are >= 0 and IDs are unique
   IF ( InputFileData%Morison%NAxCoefs > 0 ) THEN
   
      DO I = 1,InputFileData%Morison%NAxCoefs 

         IF (  InputFileData%Morison%AxialCoefs(I)%AxCd < 0 ) THEN
            CALL SetErrStat( ErrID_Fatal,'AxCd must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF   
         IF (  InputFileData%Morison%AxialCoefs(I)%AxCa < 0 ) THEN
            CALL SetErrStat( ErrID_Fatal,'AxCa must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF   
         
            ! Make sure that the current AxCoefID is not used elsewhere in the table.
         DO J = I+1,InputFileData%Morison%NAxCoefs
            IF ( InputFileData%Morison%AxialCoefs(I)%AxCoefID == InputFileData%Morison%AxialCoefs(J)%AxCoefID ) THEN
               CALL SetErrStat( ErrID_Fatal,'Duplicate AxCoefIDs were found in the Axial Coefficients table.',ErrStat,ErrMsg,RoutineName)
               RETURN
            END IF
         END DO
   
      END DO
      
   END IF


      ! Check JointOvrlp values
  !NOTE: This is ignored in the current version of Morison.  3/15/2020 GJH
   
   IF ( InputFileData%Morison%NJoints > 1 ) THEN

      ! Initialize Joints
      DO I = 1,InputFileData%Morison%NJoints
         InputFileData%Morison%InpJoints(I)%NConnections   = 0
      END DO

      
      
      
      DO I = 1,InputFileData%Morison%NJoints

            ! Make sure that the current JointID is not used elsewhere in the table.
         DO J = I+1,InputFileData%Morison%NJoints
            IF ( InputFileData%Morison%InpJoints(I)%JointID == InputFileData%Morison%InpJoints(J)%JointID ) THEN
               CALL SetErrStat( ErrID_Fatal,'Duplicate JointIDs were found in the Member Joints table.',ErrStat,ErrMsg,RoutineName)
               RETURN
            END IF
         END DO

            ! Add up total number of joints flagged with JoinOvrlp = 1 option
         !IF ( InputFileData%Morison%InpJoints(I)%JointOvrlp == 1 ) THEN
         !   InputFileData%Morison%TotalPossibleSuperMembers = InputFileData%Morison%TotalPossibleSuperMembers + 1
         !END IF

            ! Check that every joint id is used at least once in the members table
         JointUsed = .FALSE.
         DO J = 1, InputFileData%Morison%NMembers
         
            IF ( InputFileData%Morison%InpMembers(J)%MJointID1 == InputFileData%Morison%InpJoints(I)%JointID ) THEN
               JointUsed = .TRUE.
               EXIT
            END IF
            IF ( InputFileData%Morison%InpMembers(J)%MJointID2 == InputFileData%Morison%InpJoints(I)%JointID ) THEN
               JointUsed = .TRUE.
               EXIT
            END IF
         END DO
         
         IF ( .NOT. JointUsed ) THEN
            CALL SetErrStat( ErrID_Fatal,'Every JointID in the Joints table must appear once in the Members table.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF  
   ! TODO : Implement Super member elements. GJH 7/24/13
   
         IF ( InputFileData%Morison%InpJoints(I)%JointOvrlp /= 0  ) THEN
            CALL SetErrStat( ErrID_Fatal,'JointOvrlp parameter must be set to 0.  Future versions of HydroDyn will support vales of 0 or 1.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF
         !IF ( ( InputFileData%Morison%InpJoints(I)%JointOvrlp < 0 ) .OR. ( InputFileData%Morison%InpJoints(I)%JointOvrlp > 1 ) ) THEN
         !   ErrMsg  = ' JointOvrlp parameter must be set to 0 or 1.'
         !   ErrStat = ErrID_Fatal
         !   RETURN
         !END IF
         
            ! Make sure the axial coef id appears in the Ax table
         IF ( InputFileData%Morison%NAxCoefs > 0 ) THEN
            InputFileData%Morison%InpJoints(I)%JointAxIDIndx = -1
            DO J = 1,InputFileData%Morison%NAxCoefs         
               IF ( InputFileData%Morison%InpJoints(I)%JointAxID == InputFileData%Morison%AxialCoefs(J)%AxCoefID ) &
                  InputFileData%Morison%InpJoints(I)%JointAxIDIndx = J   
            END DO
            IF ( InputFileData%Morison%InpJoints(I)%JointAxIDIndx == -1 ) THEN
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

   IF ( InputFileData%Morison%NPropSets < 0 ) THEN
      CALL SetErrStat( ErrID_Fatal,'Number of member cross-section property sets must be greater than zero.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF

   IF ( InputFileData%Morison%NPropSets > 0 ) THEN

      DO I = 1,InputFileData%Morison%NPropSets

            ! Make sure that the current JointID is not used elsewhere in the table.
         DO J = I+1,InputFileData%Morison%NPropSets
            IF ( InputFileData%Morison%MPropSets(I)%PropSetID == InputFileData%Morison%MPropSets(J)%PropSetID ) THEN
               CALL SetErrStat( ErrID_Fatal,'Duplicate PropSetIDs were found in the Member Cross-section Properties table.',ErrStat,ErrMsg,RoutineName)
               RETURN
            END IF
         END DO

         IF ( ( InputFileData%Morison%MPropSets(I)%PropD < 0 ) .OR.  ( InputFileData%Morison%MPropSets(I)%PropThck < 0 ) .OR. ( ( InputFileData%Morison%MPropSets(I)%PropD - InputFileData%Morison%MPropSets(I)%PropThck / 2.0 ) < 0) ) THEN
            CALL SetErrStat( ErrID_Fatal,'PropD and PropThck must be greater than zero and (PropD - propThck/2 ) must be greater than zero.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF
      END DO

   END IF


   !-------------------------------------------------------------------------------------------------
   ! Simple hydrodynamic coefficients Section
   !-------------------------------------------------------------------------------------------------

   IF ( InputFileData%Morison%SimplCd < 0 ) THEN
      CALL SetErrStat( ErrID_Fatal,'SimplCd must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF
   IF ( InputFileData%Morison%SimplCdMG < 0 ) THEN
      CALL SetErrStat( ErrID_Fatal,'SimplCdMG must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF
   IF ( InputFileData%Morison%SimplCa < 0 ) THEN
      CALL SetErrStat( ErrID_Fatal,'SimplCa must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF
   IF ( InputFileData%Morison%SimplCaMG < 0 ) THEN
      CALL SetErrStat( ErrID_Fatal,'SimplCaMG must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF
   IF ( InputFileData%Morison%SimplAxCd < 0 ) THEN
      CALL SetErrStat( ErrID_Fatal,'SimplAxCd must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF
   IF ( InputFileData%Morison%SimplAxCdMG < 0 ) THEN
      CALL SetErrStat( ErrID_Fatal,'SimplAxCdMG must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF
   IF ( InputFileData%Morison%SimplAxCa < 0 ) THEN
      CALL SetErrStat( ErrID_Fatal,'SimplAxCa must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF
   IF ( InputFileData%Morison%SimplAxCaMG < 0 ) THEN
      CALL SetErrStat( ErrID_Fatal,'SimplAxCaMG must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF
   !TODO: Do we need a test for AxCp

   !-------------------------------------------------------------------------------------------------
   ! Depth-based Hydrodynamic Coefficients Section
   !-------------------------------------------------------------------------------------------------

   IF ( InputFileData%Morison%NCoefDpth < 0 ) THEN
      CALL SetErrStat( ErrID_Fatal,'NCoefDpth must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF


   IF ( InputFileData%Morison%NCoefDpth > 0 ) THEN
      MinDepth = 99999999.0
      MaxDepth = -99999999.0
      DO I = 1,InputFileData%Morison%NCoefDpth

            ! Record the minimum and maximum depths covered by this table.  This will be used as part of a consistency check
            ! in the members table, below.
         IF (  InputFileData%Morison%CoefDpths(I)%Dpth < MinDepth ) THEN
            MinDepth = InputFileData%Morison%CoefDpths(I)%Dpth
         ELSE
            CALL SetErrStat( ErrID_Fatal,'The rows of the Depth-based Hydrodynamic Coefficients table must be ordered with increasing depth (decreasing Z).',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF
         IF ( InputFileData%Morison%CoefDpths(I)%Dpth > MaxDepth ) THEN
            MaxDepth = InputFileData%Morison%CoefDpths(I)%Dpth
         END IF

            ! Make sure that the current Dpth is not used elsewhere in the table.
         DO J = I+1,InputFileData%Morison%NCoefDpth
            IF ( EqualRealNos( InputFileData%Morison%CoefDpths(I)%Dpth, InputFileData%Morison%CoefDpths(J)%Dpth ) ) THEN
               CALL SetErrStat( ErrID_Fatal,'Duplicate Dpths were found in the Depth-based Hydrodynamic Coefficients table.',ErrStat,ErrMsg,RoutineName)
               RETURN
            END IF
         END DO

         IF ( InputFileData%Morison%CoefDpths(I)%DpthCd < 0 ) THEN
            CALL SetErrStat( ErrID_Fatal,'In the Depth-based hydrodynamic coefficients table, DpthCd must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF
         IF ( InputFileData%Morison%CoefDpths(I)%DpthCdMG < 0 ) THEN
            CALL SetErrStat( ErrID_Fatal,'In the Depth-based hydrodynamic coefficients table, DpthCdMG must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF
         IF ( InputFileData%Morison%CoefDpths(I)%DpthCa < 0 ) THEN
            CALL SetErrStat( ErrID_Fatal,'In the Depth-based hydrodynamic coefficients table, DpthCa must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF
         IF ( InputFileData%Morison%CoefDpths(I)%DpthCaMG < 0 ) THEN
            CALL SetErrStat( ErrID_Fatal,'In the Depth-based hydrodynamic coefficients table, DpthCaMG must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF
         IF ( InputFileData%Morison%CoefDpths(I)%DpthAxCd < 0 ) THEN 
            CALL SetErrStat( ErrID_Fatal,'In the Depth-based hydrodynamic coefficients table, DpthAxCd must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF
         IF ( InputFileData%Morison%CoefDpths(I)%DpthAxCdMG < 0 ) THEN
            CALL SetErrStat( ErrID_Fatal,'In the Depth-based hydrodynamic coefficients table, DpthAxCdMG must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF
         IF ( InputFileData%Morison%CoefDpths(I)%DpthAxCa < 0 ) THEN 
            CALL SetErrStat( ErrID_Fatal,'In the Depth-based hydrodynamic coefficients table, DpthAxCa must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF
         IF ( InputFileData%Morison%CoefDpths(I)%DpthAxCaMG < 0 ) THEN
            CALL SetErrStat( ErrID_Fatal,'In the Depth-based hydrodynamic coefficients table, DpthAxCaMG must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF
         IF ( InputFileData%Morison%CoefDpths(I)%DpthAxCp < 0 ) THEN 
            CALL SetErrStat( ErrID_Fatal,'In the Depth-based hydrodynamic coefficients table, DpthAxCp must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF
         IF ( InputFileData%Morison%CoefDpths(I)%DpthAxCpMG < 0 ) THEN
            CALL SetErrStat( ErrID_Fatal,'In the Depth-based hydrodynamic coefficients table, DpthAxCpMG must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF
      END DO

      ! TODO: Sort the table based on depth so that a linear interpolation can be easily performed between entries.

   END IF


   !-------------------------------------------------------------------------------------------------
   ! Member-based Hydrodynamic Coefficients Section
   !-------------------------------------------------------------------------------------------------

   IF ( InputFileData%Morison%NCoefMembers < 0 ) THEN
      CALL SetErrStat( ErrID_Fatal,'NCoefMembers must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF

   IF ( InputFileData%Morison%NCoefMembers > 0 ) THEN

      DO I = 1,InputFileData%Morison%NCoefMembers

            ! Make sure that the current MemberID is not used elsewhere in the table.
         DO J = I+1,InputFileData%Morison%NCoefMembers
            IF ( InputFileData%Morison%CoefMembers(I)%MemberID == InputFileData%Morison%CoefMembers(J)%MemberID ) THEN
               CALL SetErrStat( ErrID_Fatal,'Duplicate MemberIDs were found in the Member-based Hydrodynamic coefficients table.',ErrStat,ErrMsg,RoutineName)
               RETURN
            END IF
         END DO



         IF ( InputFileData%Morison%CoefMembers(I)%MemberCd1 < 0 ) THEN
            CALL SetErrStat( ErrID_Fatal,'In the member-based hydrodynamic coefficients table, MemberCd1 must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF
         IF ( InputFileData%Morison%CoefMembers(I)%MemberCd2 < 0 ) THEN
            CALL SetErrStat( ErrID_Fatal,'In the member-based hydrodynamic coefficients table, MemberCd2 must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF
         IF ( InputFileData%Morison%CoefMembers(I)%MemberCdMG1 < 0 ) THEN
            CALL SetErrStat( ErrID_Fatal,'In the member-based hydrodynamic coefficients table, MemberCdMG1 must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF
         IF ( InputFileData%Morison%CoefMembers(I)%MemberCdMG2 < 0 ) THEN
            CALL SetErrStat( ErrID_Fatal,'In the member-based hydrodynamic coefficients table, MemberCdMG2 must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF
         IF ( InputFileData%Morison%CoefMembers(I)%MemberCa1 < 0 ) THEN
            CALL SetErrStat( ErrID_Fatal,'In the member-based hydrodynamic coefficients table, MemberCa1 must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF
         IF ( InputFileData%Morison%CoefMembers(I)%MemberCa2 < 0 ) THEN
            CALL SetErrStat( ErrID_Fatal,'In the member-based hydrodynamic coefficients table, MemberCa2 must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF
         IF ( InputFileData%Morison%CoefMembers(I)%MemberCaMG1 < 0 ) THEN
            CALL SetErrStat( ErrID_Fatal,'In the member-based hydrodynamic coefficients table, MemberCaMG1 must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF
         IF ( InputFileData%Morison%CoefMembers(I)%MemberCaMG2 < 0 ) THEN
            CALL SetErrStat( ErrID_Fatal,'In the member-based hydrodynamic coefficients table, MemberCaMG2 must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF
         IF ( InputFileData%Morison%CoefMembers(I)%MemberAxCa1 < 0 ) THEN 
            CALL SetErrStat( ErrID_Fatal,'In the member-based hydrodynamic coefficients table, MemberCa1 must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF
         IF ( InputFileData%Morison%CoefMembers(I)%MemberAxCa2 < 0 ) THEN 
            CALL SetErrStat( ErrID_Fatal,'In the member-based hydrodynamic coefficients table, MemberCa2 must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF
         IF ( InputFileData%Morison%CoefMembers(I)%MemberAxCaMG1 < 0 ) THEN 
            CALL SetErrStat( ErrID_Fatal,'In the member-based hydrodynamic coefficients table, MemberCaMG1 must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF
         IF ( InputFileData%Morison%CoefMembers(I)%MemberAxCaMG2 < 0 ) THEN 
            CALL SetErrStat( ErrID_Fatal,'In the member-based hydrodynamic coefficients table, MemberCaMG2 must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF
      END DO

   END IF


   !-------------------------------------------------------------------------------------------------
   ! Members Section
   !-------------------------------------------------------------------------------------------------

   IF ( InputFileData%Morison%NMembers < 0 ) THEN
      CALL SetErrStat( ErrID_Fatal,'NMembers in the Members table must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF

   IF ( InputFileData%Morison%NMembers > 0 ) THEN

         ! Initialize all member data
      DO I = 1,InputFileData%Morison%NMembers
         InputFileData%Morison%InpMembers(I)%MJointID1Indx    = -1
         InputFileData%Morison%InpMembers(I)%MJointID2Indx    = -1
         InputFileData%Morison%InpMembers(I)%MPropSetID1Indx  = -1
         InputFileData%Morison%InpMembers(I)%MPropSetID2Indx  = -1
         InputFileData%Morison%InpMembers(I)%MmbrFilledIDIndx = -1
         InputFileData%Morison%InpMembers(I)%MmbrCoefIDIndx   = -1
      END DO

      DO I = 1,InputFileData%Morison%NMembers

            ! Make sure that the current MemberID is not used elsewhere in the table.
         DO J = I+1,InputFileData%Morison%NMembers
            IF ( InputFileData%Morison%InpMembers(I)%MemberID == InputFileData%Morison%InpMembers(J)%MemberID ) THEN
               CALL SetErrStat( ErrID_Fatal,'Duplicate MemberIDs were found in the Members table.',ErrStat,ErrMsg,RoutineName)
               RETURN
            END IF
         END DO

            ! Find JointID1 and JointID2 in the Joint table and then record their index locations in the Joint table
         DO J = 1,InputFileData%Morison%NJoints
            IF ( InputFileData%Morison%InpMembers(I)%MJointID1 == InputFileData%Morison%InpJoints(J)%JointID ) THEN
               InputFileData%Morison%InpMembers(I)%MJointID1Indx = J
               InputFileData%Morison%InpJoints(J)%NConnections = InputFileData%Morison%InpJoints(J)%NConnections + 1
               InputFileData%Morison%InpJoints(J)%ConnectionList(InputFileData%Morison%InpJoints(J)%NConnections) = I
            END IF
            IF ( InputFileData%Morison%InpMembers(I)%MJointID2 == InputFileData%Morison%InpJoints(J)%JointID ) THEN
               InputFileData%Morison%InpMembers(I)%MJointID2Indx = J
               InputFileData%Morison%InpJoints(J)%NConnections = InputFileData%Morison%InpJoints(J)%NConnections + 1
               InputFileData%Morison%InpJoints(J)%ConnectionList(InputFileData%Morison%InpJoints(J)%NConnections) = -I !TODO: Come up with a better method for this work GJH 4/6/20
            END IF
         END DO
         
            ! Make sure that a JointID entry in the Joints table was found
         IF ( InputFileData%Morison%InpMembers(I)%MJointID1Indx == -1 ) THEN
            CALL SetErrStat( ErrID_Fatal,'JointID1 in the Members table does not appear in the Joints table.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF
         IF ( InputFileData%Morison%InpMembers(I)%MJointID2Indx == -1 ) THEN
            CALL SetErrStat( ErrID_Fatal,'JointID2 in the Members table does not appear in the Joints table.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF

            ! Make sure we do not have any zero length members
         lvec = InputFileData%Morison%InpJoints(InputFileData%Morison%InpMembers(I)%MJointID1Indx)%Position - InputFileData%Morison%InpJoints(InputFileData%Morison%InpMembers(I)%MJointID2Indx)%Position
         l = sqrt( lvec(1)*lvec(1) + lvec(2)*lvec(2) + lvec(3)*lvec(3) )
         IF ( EqualRealNos(0.0_ReKi, l) ) THEN
            CALL SetErrStat( ErrID_Fatal,'A member cannot have zero length.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF

            ! Find MPropSetID1 and MPropSetID2 in the Member cross-section properties table and then record their index locations
         DO J = 1,InputFileData%Morison%NPropSets



            IF ( InputFileData%Morison%InpMembers(I)%MPropSetID1 == InputFileData%Morison%MPropSets(J)%PropSetID ) THEN
               InputFileData%Morison%InpMembers(I)%MPropSetID1Indx = J
            END IF
            IF ( InputFileData%Morison%InpMembers(I)%MPropSetID2 == InputFileData%Morison%MPropSets(J)%PropSetID ) THEN
               InputFileData%Morison%InpMembers(I)%MPropSetID2Indx = J
            END IF
         END DO

            ! Make sure that a PropSetID entry in the Member cross-section properties table was found
         IF ( InputFileData%Morison%InpMembers(I)%MPropSetID1Indx == -1 ) THEN
            CALL SetErrStat( ErrID_Fatal,'MPropSetID1 in the Members table does not appear in the Member cross-section properties table.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF
         IF ( InputFileData%Morison%InpMembers(I)%MPropSetID2Indx == -1 ) THEN
            CALL SetErrStat( ErrID_Fatal,'MPropSetID2 in the Members table does not appear in the Member cross-section properties table.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF


         ! NOTE: We cannot test that MDivSize > MemberLength yet because there may be a joint overlap which is going to alter the final length of this member

         IF ( InputFileData%Morison%InpMembers(I)%MDivSize <= 0 ) THEN
            CALL SetErrStat( ErrID_Fatal,'MDivSize must be greater than zero.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF


         IF ( ( InputFileData%Morison%InpMembers(I)%MCoefMod /= 1 ) .AND. ( InputFileData%Morison%InpMembers(I)%MCoefMod /= 2 ) .AND. ( InputFileData%Morison%InpMembers(I)%MCoefMod /= 3 ) )  THEN
            CALL SetErrStat( ErrID_Fatal,'MCoefMod must be 1, 2, or 3.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF

         IF ( InputFileData%Morison%InpMembers(I)%MCoefMod == 2 ) THEN
            IF ( InputFileData%Morison%NCoefDpth == 0 ) THEN
               CALL SetErrStat( ErrID_Fatal,'NCoefDpth must be greater than zero when a member is using a depth-based coefficient model.',ErrStat,ErrMsg,RoutineName)
               RETURN
            END IF
               ! We will not extrapolate depth-based coefficient values, so make sure that the depth-based table has values that are outside the depth range of this member
               ! NOTE: This is actually potentially overly conservative because the final member may be shorter due to joint overlap handling.
            z1 = InputFileData%Morison%InpJoints( InputFileData%Morison%InpMembers(I)%MJointID1Indx )%Position(3)
            z2 = InputFileData%Morison%InpJoints( InputFileData%Morison%InpMembers(I)%MJointID2Indx )%Position(3)
            MinMembrDpth = min( z1, z2 )
            MaxMembrDpth = max( z1, z2 )
            IF ( ( MinMembrDpth < MinDepth ) .OR. ( MaxMembrDpth > MaxDepth ) ) THEN
               CALL SetErrStat( ErrID_Fatal,'This member uses a depth-based coefficient model, but the member depth is outside the range of values provided in the depth-based hydrodynamic coefficients table.',ErrStat,ErrMsg,RoutineName)
               RETURN
            END IF

         END IF


         IF ( InputFileData%Morison%InpMembers(I)%MCoefMod == 3 ) THEN
            IF ( InputFileData%Morison%NCoefMembers == 0 ) THEN
               CALL SetErrStat( ErrID_Fatal,'NCoefMembers must be greater than zero when a member is using a member-based coefficient model.',ErrStat,ErrMsg,RoutineName)
               RETURN
            END IF
               ! Make sure this id appears in the Members table and mark it's location for future use
            FoundID = .FALSE.
            DO J = 1,InputFileData%Morison%NCoefMembers
               IF ( InputFileData%Morison%CoefMembers(J)%MemberID == InputFileData%Morison%InpMembers(I)%MemberID ) THEN
                  FoundID = .TRUE.
                  InputFileData%Morison%InpMembers(I)%MmbrCoefIDIndx = J
               END IF
            END DO

            IF ( .NOT. FoundID ) THEN
               CALL SetErrStat( ErrID_Fatal,'Could not locate the MemberID referenced in the Members table in the associated Member-based Hydrodynamic coefficients table.',ErrStat,ErrMsg,RoutineName)
               RETURN
            END IF
         END IF

         IF ( InputFileData%Morison%InpMembers(I)%PropPot .AND. InputFileData%PotMod == 0  ) THEN
            CALL SetErrStat( ErrID_Fatal,'A member cannot have PropPot set to TRUE if PotMod = 0 in the FLOATING PLATFORM section.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF


         
      END DO

   END IF


   !-------------------------------------------------------------------------------------------------
   ! Filled Members Section
   !-------------------------------------------------------------------------------------------------

   IF ( InputFileData%Morison%NFillGroups < 0 ) THEN
      CALL SetErrStat( ErrID_Fatal,'NFillGroups in the Filled-members table must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF

   IF ( InputFileData%Morison%NFillGroups > 0 ) THEN

      DO I = 1,InputFileData%Morison%NFillGroups

         IF ( InputFileData%Morison%FilledGroups(I)%FillNumM < 1 ) THEN
            CALL SetErrStat( ErrID_Fatal,'FillNumM in the Filled-members table must be greater than zero.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF

         DO J = 1,InputFileData%Morison%FilledGroups(I)%FillNumM

            DO K=1,InputFileData%Morison%NMembers
               IF ( InputFileData%Morison%FilledGroups(I)%FillMList(J) == InputFileData%Morison%InpMembers(K)%MemberID ) THEN
                  FoundID = .TRUE.
                     ! Check to make sure this member is not already part of another fill group!
                  IF ( InputFileData%Morison%InpMembers(K)%MmbrFilledIDIndx /= -1 ) THEN
                     CALL SetErrStat( ErrID_Fatal,'A member cannot be a part of more than one fill group!',ErrStat,ErrMsg,RoutineName)
                  END IF

                  InputFileData%Morison%InpMembers(k)%MmbrFilledIDIndx = I

               END IF
            END DO

         END DO



            ! Make sure that the filled group members are connected
            ! NOTE: This would be easier if the input mesh was already a FAST Framework mesh because then you could use the mesh routines to determine connectivity.

            !InputFileData%Morison%FilledGroups(I)%FillMList(J)

            ! Make sure the FillFSLoc is within one of the group members
            !InputFileData%Morison%FilledGroups(I)%FillFSLoc


               ! Deal with DEFAULT or create a REAL from the string

         IF ( TRIM(InputFileData%Morison%FilledGroups(I)%FillDensChr) /= 'DEFAULT' )  THEN

            READ (InputFileData%Morison%FilledGroups(I)%FillDensChr,*,IOSTAT=IOS)  InputFileData%Morison%FilledGroups(I)%FillDens
               CALL CheckIOS ( IOS, "", 'FillDens', NumType, ErrStat2, ErrMsg2 )
               CALL SetErrStat(ErrStat2, ErrMsg2,ErrStat,ErrMsg,RoutineName)
               IF ( ErrStat >= AbortErrLev ) RETURN
         ELSE
            InputFileData%Morison%FilledGroups(I)%FillDens = InputFileData%Waves%WtrDens
         END IF

      END DO

   END IF


   !-------------------------------------------------------------------------------------------------
   ! Marine Growth by Depth Section
   !-------------------------------------------------------------------------------------------------

   IF ( InputFileData%Morison%NMGDepths < 0 ) THEN
      CALL SetErrStat( ErrID_Fatal,'NMGDepths in the Marine growth table must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF


   IF ( InputFileData%Morison%NMGDepths > 0 ) THEN

      InputFileData%Morison%MGTop    = -999999.0
      InputFileData%Morison%MGBottom =  999999.0

      DO I = 1,InputFileData%Morison%NMGDepths
            ! Store the boundaries of the marine growth zone
         IF ( InputFileData%Morison%MGDepths(I)%MGDpth > InputFileData%Morison%MGTop ) THEN
            InputFileData%Morison%MGTop    = InputFileData%Morison%MGDepths(I)%MGDpth
         END IF
         IF ( InputFileData%Morison%MGDepths(I)%MGDpth < InputFileData%Morison%MGBottom ) THEN
            InputFileData%Morison%MGBottom = InputFileData%Morison%MGDepths(I)%MGDpth
         ELSE
            CALL SetErrStat( ErrID_Fatal,'The rows of the marine growth table must be ordered with increasing depth (decreasing Z).',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF

            ! Make sure that the current MGDpth is not used elsewhere in the table.
         DO J = I+1,InputFileData%Morison%NMGDepths
            IF ( EqualRealNos( InputFileData%Morison%MGDepths(I)%MGDpth, InputFileData%Morison%MGDepths(J)%MGDpth ) ) THEN
               CALL SetErrStat( ErrID_Fatal,'Duplicate MGDpth were found in the Marine Growth table.',ErrStat,ErrMsg,RoutineName)
               RETURN
            END IF
         END DO

         IF ( InputFileData%Morison%MGDepths(I)%MGThck < 0 ) THEN
            CALL SetErrStat( ErrID_Fatal,'MGThck in the Marine growth table must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF
         IF ( InputFileData%Morison%MGDepths(I)%MGDens < 0 ) THEN
            CALL SetErrStat( ErrID_Fatal,'MGDens in the Marine growth table must be greater or equal to zero.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF

      END DO

   END IF


   !-------------------------------------------------------------------------------------------------
   ! Member Output List Section
   !-------------------------------------------------------------------------------------------------

   IF ( ( InputFileData%Morison%NMOutputs < 0 ) .OR. ( InputFileData%Morison%NMOutputs > 9 ) ) THEN
      CALL SetErrStat( ErrID_Fatal,'NMOutputs in the Member output list must be greater or equal to zero and less than 10.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF

   IF ( InputFileData%Morison%NMOutputs > 0 ) THEN


      DO I = 1,InputFileData%Morison%NMOutputs

         InputFileData%Morison%MOutLst(I)%MemberIDIndx = -1

            ! Find MemberID in this Member output list table in the Members table
         DO J = 1,InputFileData%Morison%NMembers
            IF ( InputFileData%Morison%InpMembers(J)%MemberID == InputFileData%Morison%MOutLst(I)%MemberID ) THEN
               InputFileData%Morison%MOutLst(I)%MemberIDIndx = J
            END IF
         END DO

            ! Make sure that a PropSetID entry in the Member cross-section properties table was found
         IF ( InputFileData%Morison%MOutLst(I)%MemberIDIndx == -1 ) THEN
            CALL SetErrStat( ErrID_Fatal,'MemberID in the Member output list table does not appear in the Members table.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF

         IF ( ( InputFileData%Morison%MOutLst(I)%NOutLoc < 1 ) .OR. ( InputFileData%Morison%MOutLst(I)%NOutLoc > 9) ) THEN
            CALL SetErrStat( ErrID_Fatal,'NOutLoc in the Member output list must be greater than zero and less than 10.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF

         DO J = 1,InputFileData%Morison%MOutLst(I)%NOutLoc
            IF ( ( InputFileData%Morison%MOutLst(I)%NodeLocs(J) < 0.0 ) .OR. ( InputFileData%Morison%MOutLst(I)%NodeLocs(J) > 1.0 ) ) THEN
               CALL SetErrStat( ErrID_Fatal,'NodeLocs in the Member output list must be greater or equal to 0.0 and less than or equal to 1.0.',ErrStat,ErrMsg,RoutineName)
               RETURN
            END IF
         END DO


      END DO

   END IF

   !-------------------------------------------------------------------------------------------------
   ! Joint Output List Section
   !-------------------------------------------------------------------------------------------------

   IF ( ( InputFileData%Morison%NJOutputs < 0 ) .OR. ( InputFileData%Morison%NMOutputs > 9 ) ) THEN
      CALL SetErrStat( ErrID_Fatal,'NJOutputs in the Joint output list must be greater or equal to zero and less than 10.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF

   IF ( InputFileData%Morison%NJOutputs > 0 ) THEN


      DO I=1,InputFileData%Morison%NJOutputs
           
           InputFileData%Morison%JOutLst(I)%JointIDIndx = -1
         ! Find MemberID in this Member output list table in the Members table
         DO J = 1,InputFileData%Morison%NJoints
            IF ( InputFileData%Morison%InpJoints(J)%JointID == InputFileData%Morison%JOutLst(I)%JointID ) THEN
               InputFileData%Morison%JOutLst(I)%JointIDIndx = J
               EXIT
            END IF 
         END DO
         
            ! Make sure that a Joint Output ID found in the JOutLst is in the Joints table
         IF ( InputFileData%Morison%JOutLst(I)%JointIDIndx == -1 ) THEN
            CALL SetErrStat( ErrID_Fatal,'JointID in the Joint output list table does not appear in the Joints table.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF
      END DO
   END IF
   
   !-------------------------------------------------------------------------------------------------
   ! Data section for OUTPUT
   !-------------------------------------------------------------------------------------------------


      ! OutAll - output all member and joint data

   IF ( InputFileData%OutAll ) THEN    !TODO: Alter this check once OutAll is supported
         CALL SetErrStat( ErrID_Fatal,'OutAll must be FALSE. Future versions of HydroDyn will once again support values of either TRUE or FALSE.',ErrStat,ErrMsg,RoutineName)
         RETURN
   END IF


      ! OutSwtch - output file switch

   IF ( InputFileData%OutSwtch /= 1 .AND. InputFileData%OutSwtch /= 2 .AND. InputFileData%OutSwtch /= 3 ) THEN
      CALL SetErrStat( ErrID_Fatal,'OutSwitch must be set to 1, 2, or 3.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF

   !InputFileData%OutFmt
   !InputFileData%OutSFmt


         ! OutList - list of requested parameters to output to a file


   !----------------------------------------------------------
   !  Output List
   !----------------------------------------------------------

      ! First we need to extract module-specific output lists from the user-input list.
      ! Any unidentified channels will be attached to the HydroDyn module's output list.
   IF (  InputFileData%NUserOutputs > 0 ) THEN
      ALLOCATE ( foundMask(InputFileData%NUserOutputs) , STAT = ErrStat2 )
      IF ( ErrStat2 /= ErrID_None ) THEN
         CALL SetErrStat( ErrID_Fatal,'Error allocating space for temporary array: foundMask in the HydroDynInput_GetInput subroutine.',ErrStat,ErrMsg,RoutineName)
         
         RETURN
      END IF
      foundMask = .FALSE.
         ! Extract Waves2 list
      InputFileData%Waves2%NumOuts  = GetWaves2Channels   ( InputFileData%NUserOutputs, InputFileData%UserOutputs, InputFileData%Waves2%OutList, foundMask, ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   
!         ! Extract WAMIT2 list
!      InputFileData%WAMIT2%NumOuts  = GetWAMIT2Channels   ( InputFileData%NUserOutputs, InputFileData%UserOutputs, InputFileData%WAMIT2%OutList, foundMask, ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
!
         ! Extract Morison list
         !foundMask = .FALSE.
      InputFileData%Morison%NumOuts = GetMorisonChannels  ( InputFileData%NUserOutputs, InputFileData%UserOutputs, InputFileData%Morison%OutList, foundMask, ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   
         ! Attach remaining items to the HydroDyn list
         !foundMask = .FALSE.
      call Allocary(InputFileData%OutList, InputFileData%NUserOutputs, "InputFileData%OutList", ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      InputFileData%NumOuts       = HDOut_GetChannels ( InputFileData%NUserOutputs, InputFileData%UserOutputs, InputFileData%OutList        , foundMask, ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      CALL PrintBadChannelWarning(InputFileData%NUserOutputs, InputFileData%UserOutputs , foundMask, ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
      IF (ErrStat >= AbortErrLev ) RETURN

      DEALLOCATE(foundMask)
      
   END IF
      ! Now that we have the sub-lists organized, lets do some additional validation.
   
   
   
   
   !----------------------------------------------------------
   ! Mesh-related Output List
   !----------------------------------------------------------

   IF ( InputFileData%Morison%NumOuts > 0 ) THEN

         ! Create an  output list for validated outputs
      ALLOCATE ( InputFileData%Morison%ValidOutList(InputFileData%Morison%NumOuts), STAT = ErrStat2 )
      IF ( ErrStat2 /= 0 ) THEN
         CALL SetErrStat( ErrID_Fatal,'Error allocating valid output list array.',ErrStat,ErrMsg,RoutineName)
         RETURN
      END IF

      DO I =1, InputFileData%Morison%NumOuts

         InputFileData%Morison%ValidOutList(I) = CheckMeshOutput( InputFileData%Morison%OutList(I), InputFileData%Morison%NMOutputs, InputFileData%Morison%MOutLst, InputFileData%Morison%NJOutputs )

      END DO

   END IF


   !----------------------------------------------------------
   ! Populate data in sub-types from parent or other module types
   !----------------------------------------------------------

      ! Current
         ! For wave kinematic calculations, the effective water depth is the user input water depth (positive valued) + MSL2SWL (positive when SWL is above MSL).
      InputFileData%Current%WtrDpth    = InputFileData%Morison%WtrDpth ! already adjusted for the MSL2SWL.  
                                                       
      
      ! Waves
      InputFileData%Waves%Gravity      = InitInp%Gravity
      InputFileData%Waves%UnSum        = InputFileData%UnSum
         ! For wave kinematic calculations, the effective water depth is the user input water depth (positive valued) + MSL2SWL (positive when SWL is above MSL).
      InputFileData%Waves%WtrDpth      = InputFileData%Morison%WtrDpth ! already adjusted for the MSL2SWL.
      
      ! Waves2
      IF (InputFileData%Waves2%WvDiffQTFF .OR. InputFileData%Waves2%WvSumQTFF ) THEN
         InputFileData%Waves2%WtrDens     = InputFileData%Waves%WtrDens
         InputFileData%Waves2%Gravity     = InitInp%Gravity
         InputFileData%Waves2%UnSum       = InputFileData%UnSum
         InputFileData%Waves2%WtrDpth     = InputFileData%Waves%WtrDpth
         InputFileData%Waves2%WaveStMod   = InputFileData%Waves%WaveStMod
         InputFileData%Waves2%NWaveElev   = InputFileData%Waves%NWaveElev
         CALL AllocAry( InputFileData%Waves2%WaveElevxi, InputFileData%Waves2%NWaveElev, 'WaveElevxi' , ErrStat2, ErrMsg2)
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
         CALL AllocAry( InputFileData%Waves2%WaveElevyi, InputFileData%Waves2%NWaveElev, 'WaveElevyi' , ErrStat2, ErrMsg2)
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
         IF ( ErrStat >= AbortErrLev ) RETURN 
         InputFileData%Waves2%WaveElevxi  = InputFileData%Waves%WaveElevxi
         InputFileData%Waves2%WaveElevyi  = InputFileData%Waves%WaveElevyi
      ENDIF

      ! WAMIT
      InputFileData%WAMIT%WtrDens      = InputFileData%Waves%WtrDens
      InputFileData%WAMIT%WaveMod      = InputFileData%Waves%WaveMod
      InputFileData%WAMIT%OutAll       = InputFileData%OutAll
      InputFileData%WAMIT%HasWAMIT     = InputFileData%PotMod == 1
      ! WAMIT2
      InputFileData%WAMIT2%WtrDens     = InputFileData%Waves%WtrDens
      InputFileData%WAMIT2%WaveMod     = InputFileData%Waves%WaveMod
      InputFileData%WAMIT2%HasWAMIT    = InputFileData%PotMod == 1
      ! Morison
      InputFileData%Morison%UnSum      = InputFileData%UnSum
      InputFileData%Morison%Gravity    = InitInp%Gravity
      InputFileData%Morison%WtrDens    = InputFileData%Waves%WtrDens
      InputFileData%Morison%OutAll     = InputFileData%OutAll

         ! Process the input geometry and generate the simulation mesh representation
      call Morison_GenerateSimulationNodes( InputFileData%Morison%MSL2SWL, InputFileData%Morison%NJoints, InputFileData%Morison%InpJoints, InputFileData%Morison%NMembers, InputFileData%Morison%InpMembers, InputFileData%Morison%NNodes, InputFileData%Morison%Nodes, errStat2, errMsg2 )
      !CALL Morison_ProcessMorisonGeometry( InputFileData%Morison, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDynInput_GetInput' )
      IF ( ErrStat >= AbortErrLev ) RETURN

         ! Set the number and global Z locations for the X and Y components of the current velocities
      InputFileData%Current%NMorisonNodes = InputFileData%Morison%NNodes

      ALLOCATE ( InputFileData%Current%MorisonNodezi(InputFileData%Morison%NNodes), STAT = ErrStat2 )
      IF ( ErrStat2 /= ErrID_None ) THEN
         CALL SetErrStat( ErrID_Fatal,'Error allocating space for MorisonNodezi array.',ErrStat,ErrMsg,RoutineName)
         RETURN
      END IF



         ! Establish the number and locations where the wave kinematics will be computed
      InputFileData%Waves%NWaveKin   = InputFileData%Morison%NNodes                          ! Number of points where the incident wave kinematics will be computed (-)
      ALLOCATE ( InputFileData%Waves%WaveKinxi(InputFileData%Waves%NWaveKin), STAT = ErrStat2 )
      IF ( ErrStat2 /= ErrID_None ) THEN
         CALL SetErrStat( ErrID_Fatal,'Error allocating space for WaveKinxi array.',ErrStat,ErrMsg,RoutineName)

         RETURN
      END IF
      ALLOCATE ( InputFileData%Waves%WaveKinyi(InputFileData%Waves%NWaveKin), STAT = ErrStat2 )
      IF ( ErrStat2 /= ErrID_None ) THEN
         CALL SetErrStat( ErrID_Fatal,'Error allocating space for WaveKinyi array.',ErrStat,ErrMsg,RoutineName)
         RETURN
      END IF
      ALLOCATE ( InputFileData%Waves%WaveKinzi(InputFileData%Waves%NWaveKin), STAT = ErrStat2 )
      IF ( ErrStat2 /= ErrID_None ) THEN
         CALL SetErrStat( ErrID_Fatal,'Error allocating space for WaveKinzi array.',ErrStat,ErrMsg,RoutineName)
         RETURN
      END IF
      DO I=1,InputFileData%Morison%NNodes
         InputFileData%Waves%WaveKinxi(I)      = InputFileData%Morison%Nodes(I)%Position(1)   ! xi-coordinates for points where the incident wave kinematics will be computed;
         InputFileData%Waves%WaveKinyi(I)      = InputFileData%Morison%Nodes(I)%Position(2)   ! yi-coordinates for points where the incident wave kinematics will be computed;
         InputFileData%Waves%WaveKinzi(I)      = InputFileData%Morison%Nodes(I)%Position(3)   ! zi-coordinates for points where the incident wave kinematics will be computed; 
         InputFileData%Current%MorisonNodezi(I) = InputFileData%Waves%WaveKinzi(I)
      END DO


            ! If we are using the Waves module, the node information must be copied over.
      InputFileData%Waves2%NWaveKin   = InputFileData%Waves%NWaveKin                          ! Number of points where the incident wave kinematics will be computed (-)
      IF ( InputFileData%Waves2%WvDiffQTFF .OR. InputFileData%Waves2%WvSumQTFF ) THEN
         ALLOCATE ( InputFileData%Waves2%WaveKinxi(InputFileData%Waves2%NWaveKin), STAT = ErrStat2 )
         IF ( ErrStat2 /= ErrID_None ) THEN
            CALL SetErrStat( ErrID_Fatal,'Error allocating space for WaveKinxi array for Waves2 module.',ErrStat,ErrMsg,RoutineName)

            RETURN
         END IF
         ALLOCATE ( InputFileData%Waves2%WaveKinyi(InputFileData%Waves2%NWaveKin), STAT = ErrStat2 )
         IF ( ErrStat2 /= ErrID_None ) THEN
            CALL SetErrStat( ErrID_Fatal,'Error allocating space for WaveKinyi array for Waves2 module.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF
         ALLOCATE ( InputFileData%Waves2%WaveKinzi(InputFileData%Waves2%NWaveKin), STAT = ErrStat2 )
         IF ( ErrStat2 /= ErrID_None ) THEN
            CALL SetErrStat( ErrID_Fatal,'Error allocating space for WaveKinzi array for Waves2 module.',ErrStat,ErrMsg,RoutineName)
            RETURN
         END IF

         InputFileData%Waves2%WaveKinxi  = InputFileData%Waves%WaveKinxi
         InputFileData%Waves2%WaveKinyi  = InputFileData%Waves%WaveKinyi
         InputFileData%Waves2%WaveKinzi  = InputFileData%Waves%WaveKinzi

      ENDIF

END SUBROUTINE HydroDynInput_ProcessInitData

END MODULE HydroDyn_Input
