!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2013  National Renewable Energy Laboratory
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
! File last committed: $Date: 2013-10-03 12:00:12 -0600 (Thu, 03 Oct 2013) $
! (File) Revision #: $Rev: 258 $
! URL: $HeadURL: https://windsvn.nrel.gov/HydroDyn/branches/HydroDyn_Modularization/Source/HydroDyn_Input.f90 $
!**********************************************************************************************************************************
MODULE HydroDyn_Input

      ! This MODULE stores variables used for input.

   USE                              NWTC_Library
   USE                              HydroDyn_Types
   USE                              Waves
   USE                              Morison
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
   
   outputTmp         = TRIM(output)
   
   
      ! Reverse the sign (+/-) of the output channel if the user prefixed the
      !   channel name with a '-', '_', 'm', or 'M' character indicating "minus".
      
      IF      ( INDEX( '-_', outputTmp(1:1) ) > 0 ) THEN
         
            ! ex, '-TipDxc1' causes the sign of TipDxc1 to be switched.
         outputTmp                   = outputTmp(2:)
         
      ELSE IF ( INDEX( 'mM', outputTmp(1:1) ) > 0 ) THEN ! We'll assume this is a variable name for now, (if not, we will check later if OutListTmp(2:) is also a variable name)
   
         IF ( ( INDEX( 'mM', outputTmp(2:2) ) > 0 ) .OR. ( INDEX( 'jJ', outputTmp(2:2) ) > 0 ) )  THEN
            outputTmp                   = outputTmp(2:)
         
         END IF       
         
      ELSE IF ( INDEX( 'jJ', outputTmp(1:1) ) == 0 ) THEN
         ! Invalid output label because the label does not start: -M,-m,-J,-j,_M,_m,_J,_j,MM,mM,Mm,mm,MJ,mJ,Mj,mj, j,J,m,M
         CheckMeshOutput = .FALSE.
         RETURN
      END IF
      
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
            
         
      ENDIF 
      
      IF ( INDEX( 'jJ', outputTmp(1:1) ) > 0 ) THEN 
         IF ( indx1 > numJointOut ) THEN
            CheckMeshOutput = .FALSE.
            RETURN
         END IF
      ENDIF 
      
      CheckMeshOutput = .TRUE.
      
END FUNCTION CheckMeshOutput

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
   INTEGER                                          :: J                    ! generic integer for counting
   CHARACTER(   2)                                  :: strI                 ! string version of the loop counter

   INTEGER                                          :: UnIn                 ! Unit number for the input file
   LOGICAL                                          :: EchoStore            ! Stored version of NWTC_Library Echo variable
   INTEGER                                          :: UnEchoStore          ! Stored unit name for another module's echo file
   INTEGER                                          :: UnEchoLocal          ! The local unit number for this module's echo file
   CHARACTER(1024)                                  :: EchoFile             ! Name of HydroDyn echo file  
   CHARACTER(1024)                                  :: Line                 ! String to temporarially hold value of read line   
   CHARACTER(1024)                                  :: TmpPath              ! Temporary storage for relative path name
   CHARACTER(1024)                                  :: TmpFmt               ! Temporary storage for format statement
   CHARACTER(1024)                                  :: FileName             ! Name of HydroDyn input file  
   CHARACTER(  35)                                  :: Frmt                 ! Output format for logical parameters. (matches NWTC Subroutine Library format)
   INTEGER                                          :: JointID              ! Temporary storage of JointID read from HydroDyn input file
   INTEGER                                          :: PropSetID            ! Temporary storage of PropSetID read from HydroDyn input file
   INTEGER                                          :: MemberID             ! Temporary storage of MemberID read from HydroDyn input file
   INTEGER, ALLOCATABLE                             :: tmpArray(:)          ! Temporary array storage of the joint output list
   
   
      ! Initialize local data
      
   UnEchoLocal = -1
   Frmt      = "( 2X, L11, 2X, A, T30, ' - ', A )"
   
   
      ! Initialize ErrStat
         
   ErrStat = ErrID_None         
   ErrMsg  = ""   
   
   
   !-------------------------------------------------------------------------------------------------
   ! Open the file
   !-------------------------------------------------------------------------------------------------   
   FileName = TRIM(InitInp%InputFile)
   
   CALL GetNewUnit( UnIn )   
   CALL OpenFInpFile( UnIn, FileName, ErrStat )
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to open HydroDyn input file: '//FileName
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF

   
   !CALL WrScr( 'Opening HydroDyn input file:  '//FileName )
   
   
   !-------------------------------------------------------------------------------------------------
   ! File header
   !-------------------------------------------------------------------------------------------------
   
   CALL ReadCom( UnIn, FileName, 'HydroDyn input file header line 1', ErrStat )
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read HydroDyn input file header line 1.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF


   CALL ReadCom( UnIn, FileName, 'HydroDyn input file header line 2', ErrStat )
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read HydroDyn input file header line 2.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF

   
     ! Echo Input Files.
      
   CALL ReadVar ( UnIn, FileName, InitInp%Echo, 'Echo', 'Echo Input', ErrStat )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Echo parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
   
      ! If we are Echoing the input then we should re-read the first three lines so that we can echo them
      ! using the NWTC_Library routines.  The echoing is done inside those routines via a global variable
      ! which we must store, set, and then replace on error or completion.
      
   IF ( InitInp%Echo ) THEN
      
      EchoFile = TRIM(FileName)//'.echo'
      CALL GetNewUnit( UnEchoLocal )   
      CALL OpenEcho ( UnEchoLocal, EchoFile, ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) THEN
         !ErrMsg  = ' Failed to open Echo file.'
         ErrStat = ErrID_Fatal
         CLOSE( UnIn )
         RETURN
      END IF
      
      REWIND(UnIn)
      
      CALL ReadCom( UnIn, FileName, 'HydroDyn input file header line 1', ErrStat, ErrMsg, UnEchoLocal )
   
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Failed to read HydroDyn input file header line 1.'
         ErrStat = ErrID_Fatal
         CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      END IF


      CALL ReadCom( UnIn, FileName, 'HydroDyn input file header line 2', ErrStat, ErrMsg, UnEchoLocal )
   
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Failed to read HydroDyn input file header line 2.'
         ErrStat = ErrID_Fatal
         CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      END IF

   
         ! Echo Input Files. Note this line is prevented from being echoed by the ReadVar routine.
      
      CALL ReadVar ( UnIn, FileName, InitInp%Echo, 'Echo', 'Echo the input file data', ErrStat, ErrMsg, UnEchoLocal )
      !WRITE (UnEchoLocal,Frmt      ) InitInp%Echo, 'Echo', 'Echo input file'
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Failed to read Echo parameter.'
         ErrStat = ErrID_Fatal
         CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      END IF
      
   END IF
   !-------------------------------------------------------------------------------------------------
   ! Environmental conditions section
   !-------------------------------------------------------------------------------------------------

      ! Header
      
   CALL ReadCom( UnIn, FileName, 'Environmental conditions header', ErrStat, ErrMsg, UnEchoLocal )
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Comment line.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF


      ! WtrDens - Water density.
      
   CALL ReadVar ( UnIn, FileName, InitInp%Waves%WtrDens, 'WtrDens', 'Water density', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read WtrDens parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF

         
      ! WtrDpth - Water depth   
      
   CALL ReadVar ( UnIn, FileName, InitInp%Morison%WtrDpth, 'WtrDpth', 'Water depth', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read WtrDpth parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF

      
      ! MSL2SWL

   CALL ReadVar ( UnIn, FileName, InitInp%Morison%MSL2SWL, 'MSL2SWL', 'MSL to SWL offset', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read MSL2SWL parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
   

   !-------------------------------------------------------------------------------------------------
   ! Data section for waves
   !-------------------------------------------------------------------------------------------------
      
      ! Header
      
   CALL ReadCom( UnIn, FileName, 'Wave header', ErrStat, ErrMsg, UnEchoLocal )
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Wave header comment line.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF


!BJJ: verify that these tests can be performed (i.e., do we have the correct data to compare here?)
!   IF ( InitInp%StrctType == FixedBtm_Type ) THEN
!
!!!JASON: WHAT LOADING DO WE APPLY TO THE FLEXIBLE PORTION OF THE TOWER EXTENDING BELOW THE SEABED?
! bjj: replace this :
!         IF ( ( TwrDraft - TwrRBHt ) < WtrDpth )  THEN   ! Print out a warning when the flexible portion of the support structure does not extend to the seabed.
! with this:
!         IF ( ( HydroConfig%Substructure%Position(3) ) < -WtrDpth )  THEN   ! Print out a warning when the flexible portion of the support structure does not extend to the seabed.
!            CALL ProgWarn( ' Hydrodynamic loading will only be applied to the flexible portion of the support structure.'// &
!                           ' Make sure that ( TwrDraft - TwrRBHt ) >= WtrDpth if you want hydrodynamic loading applied'// &
!                           ' along the entire submerged portion of the support structure. ')
!         END IF
!      
!   ELSE IF ( InitInp%StrctType == FloatPltfm_Type ) THEN
!   
!      IF ( WtrDpth <= PtfmDraft  )  THEN
!         ErrMsg  = ' WtrDpth must be greater than PtfmDraft.'
!         CLOSE( UnIn )
!         RETURN
!      END IF
!         
!      IF ( FP_InitData%LineMod == 1 )  THEN  ! .TRUE if we have standard quasi-static mooring lines.
!         DO I = 1,FP_InitData%NumLines ! Loop through all mooring lines
!
!            IF ( WtrDpth < -FP_InitData%MooringLine(I)%LAnchzi )  THEN
!               ErrMsg  = ' WtrDpth must not be less than LDpthAnch('//TRIM( Int2LStr( I ) )//').'
!               ErrStat = ErrID_Fatal
!               CLOSE( UnIn )
!               RETURN
!            END IF
!               
!         END DO             ! I - All mooring lines
!      END IF
!      
!   END IF


      
      ! WaveMod - Wave kinematics model switch.

   CALL ReadVar ( UnIn, FileName, InitInp%Waves%WaveModChr, 'WaveMod', 'Wave kinematics model switch', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read WaveMod parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
   
   CALL Conv2UC( InitInp%Waves%WaveModChr )    ! Convert Line to upper case.
   
   InitInp%Waves%WavePhase = 0.0
   InitInp%Waves%WaveNDAmp = .FALSE.   
   
      
      ! WaveStMod - Model switch for stretching incident wave kinematics to instantaneous free surface. 
  
   CALL ReadVar ( UnIn, FileName, InitInp%Waves%WaveStMod, 'WaveStMod', &
      'Model switch for stretching incident wave kinematics to instantaneous free surface', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read WaveStMod parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF        

      
           
      ! WaveTMax - Analysis time for incident wave calculations.  
      
   
   CALL ReadVar ( UnIn, FileName, InitInp%Waves%WaveTMax, 'WaveTMax', &
                              'Analysis time for incident wave calculations', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read WaveTMax parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
         

            
      
      ! WaveDT - Time step for incident wave calculations   
      
   
   CALL ReadVar ( UnIn, FileName, InitInp%Waves%WaveDT, 'WaveDT', &
                        'Time step for incident wave calculations', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read WaveDT parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF

   

      
      ! WaveHs - Significant wave height    
      
   
   CALL ReadVar ( UnIn, FileName, InitInp%Waves%WaveHs, 'WaveHs', 'Significant wave height', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read WaveHs parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF


      
      ! WaveTp - Peak spectral period.

   CALL ReadVar ( UnIn, FileName, InitInp%Waves%WaveTp, 'WaveTp', 'Peak spectral period', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read WaveTp parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF

   
      ! WavePkShp - Peak shape parameter.
      
   CALL ReadVar ( UnIn, FileName, InitInp%Waves%WavePkShpChr, 'WavePkShp', 'Peak shape parameter', ErrStat, ErrMsg, UnEchoLocal ) 
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read WavePkShp parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
   
   
      ! WvLowCOff - Low Cut-off frequency or lower frequency limit of the wave spectrum beyond which the wave spectrum is zeroed (rad/s).  
      
   CALL ReadVar ( UnIn, FileName, InitInp%Waves%WvLowCOff, 'WvLowCOff', 'Lower wave cut-off frequency', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read WvLowCOff parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF   
      
     ! WvHiCOff - High Cut-off frequency or upper frequency limit of the wave spectrum beyond which the wave spectrum is zeroed (rad/s).  
      
   CALL ReadVar ( UnIn, FileName, InitInp%Waves%WvHiCOff, 'WvHiCOff', 'Upper wave cut-off frequency', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read WvHiCOff parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF    
   
      ! WaveDir - Wave heading direction.  
      
   CALL ReadVar ( UnIn, FileName, InitInp%Waves%WaveDir, 'WaveDir', 'Wave heading direction', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read WaveDir parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF   
      
      
      ! WaveSeed(1), !WaveSeed(2)

   DO I = 1,2
      
      WRITE(Line,'(I2)') I
      
      CALL ReadVar ( UnIn, FileName, InitInp%Waves%WaveSeed(I), 'WaveSeed('//TRIM(Line)//')', &
                                    'Random seed #'//TRIM(Line), ErrStat, ErrMsg, UnEchoLocal )
      
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Failed to read WaveSeed parameter.'
         ErrStat = ErrID_Fatal
         CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      END IF
      
   END DO !I


      ! WaveNDAmp - Flag for normally distributed amplitudes.  
      
   CALL ReadVar ( UnIn, FileName, InitInp%Waves%WaveNDAmp, 'WaveNDAmp', 'Normally distributed amplitudes', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read WaveNDAmp parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF   
      
      
      
      ! GHWvFile   
 
   CALL ReadVar ( UnIn, FileName, InitInp%Waves%GHWvFile, 'GHWvFile', &
                                    'Root name of GH Bladed files containing wave data', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read GHWvFile parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
   
   


      ! NWaveElev

   CALL ReadVar ( UnIn, FileName, InitInp%Waves%NWaveElev, 'NWaveElev', &
                                  'Number of points where the incident wave elevations can be output', ErrStat, ErrMsg, UnEchoLocal )
      
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read NWaveElev parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
      
   END IF
   
   
      ! This check is needed here instead of being located in HydroDynInput_ProcessInputData() because
      ! we need to allocate arrays.  If _GetInput() was skipped, then these array would already have
      ! been allocated and populated.
      
   IF ( InitInp%Waves%NWaveElev < 0 .OR. InitInp%Waves%NWaveElev > 9 ) THEN
      
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal ) 
      ErrMsg  = ' NWaveElev must be greater than or equal to zero and less than 10.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
      
   ELSE
      
      
         ! allocate space for the output location arrays: 
      
      ALLOCATE ( InitInp%Waves%WaveElevxi(InitInp%Waves%NWaveElev), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating space for WaveElevxi array.'
         ErrStat = ErrID_Fatal
         CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      END IF
      
      ALLOCATE ( InitInp%Waves%WaveElevyi(InitInp%Waves%NWaveElev), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating space for WaveElevyi array.'
         ErrStat = ErrID_Fatal
         CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      END IF
      
   END IF

      ! WaveElevxi
         
   CALL ReadAry ( UnIn, FileName, InitInp%Waves%WaveElevxi, InitInp%Waves%NWaveElev, 'WaveElevxi', &
                           'List of xi-coordinates for points where the incident wave elevations can be output', ErrStat,  ErrMsg, UnEchoLocal)         
       
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read WaveElevxi parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
      
      
      ! WaveElevyi

   CALL ReadAry ( UnIn, FileName, InitInp%Waves%WaveElevyi, InitInp%Waves%NWaveElev, 'WaveElevyi', &
                           'List of yi-coordinates for points where the incident wave elevations can be output', ErrStat,  ErrMsg, UnEchoLocal )         
      
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read WaveElevyi parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF

  

   !-------------------------------------------------------------------------------------------------
   ! Data section for current
   !-------------------------------------------------------------------------------------------------

      ! Header
      
   CALL ReadCom( UnIn, FileName, 'Current header', ErrStat, ErrMsg, UnEchoLocal )
   
   IF (ErrStat /= ErrID_None) THEN
      ErrMsg  = ' Failed to read Current header line.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF


      ! CurrMod - Current profile model switch
      
   CALL ReadVar ( UnIn, FileName, InitInp%Current%CurrMod, 'CurrMod', 'Current profile model switch', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read CurrMod parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF

   

      ! CurrSSV0 - Sub-surface current velocity at still water level
      
   CALL ReadVar ( UnIn, FileName, InitInp%Current%CurrSSV0, 'CurrSSV0', 'Sub-surface current velocity at still water level', &
                         ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read CurrSSV0 parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
  
   
     
      ! CurrSSDirChr - Sub-surface current heading direction
      
   CALL ReadVar ( UnIn, FileName, InitInp%Current%CurrSSDirChr, 'CurrSSDirChr', 'Sub-surface current heading direction', ErrStat, ErrMsg, UnEchoLocal )
      
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read CurrSSDirChr parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
            
            
   CALL Conv2UC( InitInp%Current%CurrSSDirChr )    ! Convert Line to upper case.
      
   
      ! CurrNSRef - Near-surface current reference depth.
   
   CALL ReadVar ( UnIn, FileName, InitInp%Current%CurrNSRef, 'CurrNSRef', 'Near-surface current reference depth', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read CurrNSRef parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
   
   
      ! CurrNSV0 - Near-surface current velocity at still water level.
   
   CALL ReadVar ( UnIn, FileName, InitInp%Current%CurrNSV0, 'CurrNSV0', 'Near-surface current velocity at still water level', &
                           ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read CurrNSV0 parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF

   
      ! CurrNSDir - Near-surface current heading direction.

   CALL ReadVar ( UnIn, FileName, InitInp%Current%CurrNSDir, 'CurrNSDir', 'Near-surface current heading direction', ErrStat, ErrMsg, UnEchoLocal )
  
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read CurrNSDir parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
   
   
      ! CurrDIV - Depth-independent current velocity.

   CALL ReadVar ( UnIn, FileName, InitInp%Current%CurrDIV, 'CurrDIV', 'Depth-independent current velocity', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read CurrDIV parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
   
   
      ! CurrDIDir - Depth-independent current heading direction.
   
   CALL ReadVar ( UnIn, FileName, InitInp%Current%CurrDIDir, 'CurrDIDir', 'Depth-independent current heading direction', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read CurrDIDir parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF


   
   !-------------------------------------------------------------------------------------------------
   ! Data section for floating platform
   !-------------------------------------------------------------------------------------------------

      ! Header
      
   CALL ReadCom( UnIn, FileName, 'Floating platform header', ErrStat, ErrMsg, UnEchoLocal )
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Floating platform header line.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
   !!!!!!!!!!!!!!!!!!!!!!!!
   
!   IF ( HD_Data%StrctType == FloatPltfm_Type ) THEN ! (should this be done in FAST?)
!
!instead of the following 2 checks, now check that the HydroConfig marker for platform is at 0,0,0
!      IF ( TwrDraft > 0.0 ) THEN
!         ErrMsg  = ' TwrDraft must be less than or equal to zero when PtfmLdMod is set to "'//TRIM(Line)//'".'  ! Do not allow the combination of tower hydrodynamics using Morison's equation and platform hydrodynamics using the true form of the using the true form of the hydrodynamics equations since the true equations require that the shape of the platform does not change above the MSL (platform reference point)--Consider the linear hydrostatic restoring matrix, for example.
!         ErrStat = ErrID_Fatal
!         CLOSE( UnIn )
!         RETURN
!      END IF
!
!      IF ( PtfmRef /= 0.0 ) THEN
!         ErrMsg  = ' PtfmRef must be zero when PtfmLdMod is set to "'//TRIM(Line)//'".'
!         ErrStat = ErrID_Fatal
!         CLOSE( UnIn )
!         RETURN
!      END IF
!bjj: like this:
!      IF ( HydroConfig%Substructure%Position  /= ( 0,0,0 ) .OR.
!           HydroConfig%Substructure%Orientation /= eye(3) ) THEN
!      
!         ErrMsg  = ' HydroConfig%Substructure%Position must be zero and Orientation must be the identity matrix when PtfmLdMod is set to "'//TRIM(Line)//'".'
!           
!      END IF
!
!   END IF

      ! HasWAMIT - Flag indicating whether or not WAMIT is used in the simulation.  
      
   CALL ReadVar ( UnIn, FileName, InitInp%HasWAMIT, 'HasWAMIT', 'Using WAMIT', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read HasWAMIT parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF   

   
      ! WAMITFile - Root name of WAMIT output files
 
   CALL ReadVar ( UnIn, FileName, InitInp%WAMIT%WAMITFile, 'WAMITFile', 'Root name of WAMIT output files', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read WAMITFile parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF

   
      ! WAMITULEN - WAMIT characteristic body length scale
      
   CALL ReadVar ( UnIn, FileName, InitInp%WAMIT%WAMITULEN, 'WAMITULEN', 'WAMIT characteristic body length scale', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read WAMITULEN parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF       


      ! PtfmVol0 - Displaced volume of water when the platform is in its undisplaced position

   CALL ReadVar ( UnIn, FileName, InitInp%WAMIT%PtfmVol0, 'PtfmVol0', &
      'Displaced volume of water when the platform is in its undisplaced position', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read PtfmVol0 parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF

   
      ! PtfmCOBxt  - The xt offset of the center of buoyancy (COB) from the platform reference point

   CALL ReadVar ( UnIn, FileName, InitInp%WAMIT%PtfmCOBxt, 'PtfmCOBxt', &
      'xt offset of the center of buoyancy (COB) from the platform reference point', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read PtfmCOBxt parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
   
   
      ! PtfmCOByt - The yt offset of the center of buoyancy (COB) from the platform reference point

   CALL ReadVar ( UnIn, FileName, InitInp%WAMIT%PtfmCOByt, 'PtfmCOByt', &
      'yt offset of the center of buoyancy (COB) from the platform reference point', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read PtfmCOByt parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF

   
      ! RdtnMod  - Radiation memory-effect model {1: convolution, 2: state-space} (switch) 
      ! [STATE-SPACE REQUIRES *.ss INPUT FILE]
      
  CALL ReadVar ( UnIn, FileName, InitInp%WAMIT%RdtnMod, 'RdtnMod', &
                                 'Radiation memory-effect model', ErrStat, ErrMsg, UnEchoLocal )  
  IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read RdtnMod parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
  END IF
  
  
      ! RdtnTMax - Analysis time for wave radiation kernel calculations
      ! NOTE: Use RdtnTMax = 0.0 to eliminate wave radiation damping
   
   CALL ReadVar ( UnIn, FileName, InitInp%WAMIT%RdtnTMax, 'RdtnTMax', &
                                 'Analysis time for wave radiation kernel calculations', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read RdtnTMax parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
      

      ! RdtnDT - Time step for wave radiation kernel calculations

   
   CALL ReadVar ( UnIn, FileName, InitInp%WAMIT%Conv_Rdtn%RdtnDT, 'RdtnDT', 'Time step for wave radiation kernel calculations', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read RdtnDT parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
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
   ! Data section for Floating platform force flags
   !-------------------------------------------------------------------------------------------------

      ! Header
      
   CALL ReadCom( UnIn, FileName, 'Floating platform force flags header', ErrStat, ErrMsg, UnEchoLocal )
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Floating platform force flags header line.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
   
   
       ! PtfmSgFChr - Platform horizontal surge translation force flag
   
   CALL ReadVar ( UnIn, FileName, InitInp%PtfmSgFChr, 'PtfmSgFChr', 'Platform horizontal surge translation force flag', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read PtfmSgF parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
 
   CALL Conv2UC( InitInp%PtfmSgFChr )    ! Convert Line to upper case.
   
   
      ! PtfmSwFChr - Platform horizontal sway translation force flag
   
   CALL ReadVar ( UnIn, FileName, InitInp%PtfmSwFChr, 'PtfmSwFChr', 'Platform horizontal sway translation force flag', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read PtfmSwF parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
   
   CALL Conv2UC( InitInp%PtfmSwFChr )    ! Convert Line to upper case.
   
   
       ! PtfmHvFChr - Platform vertical heave translation force flag
  
   CALL ReadVar ( UnIn, FileName, InitInp%PtfmHvFChr, 'PtfmHvFChr', 'Platform vertical heave translation force flag', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read PtfmHvF parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
   
   CALL Conv2UC( InitInp%PtfmHvFChr )    ! Convert Line to upper case.
   
   
        ! PtfmRFChr - Platform roll tilt rotation force flag
   
   CALL ReadVar ( UnIn, FileName, InitInp%PtfmRFChr, 'PtfmRFChr', 'Platform roll tilt rotation force flag', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read PtfmRF parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
   
   CALL Conv2UC( InitInp%PtfmRFChr )    ! Convert Line to upper case.
   
   
        ! PtfmPFChr - Platform pitch tilt rotation force flag
   
   CALL ReadVar ( UnIn, FileName, InitInp%PtfmPFChr, 'PtfmPFChr', 'Platform pitch tilt rotation force flag', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read PtfmPF parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
   
   CALL Conv2UC( InitInp%PtfmPFChr )    ! Convert Line to upper case.
   
   
        ! PtfmYFChr - Platform yaw rotation force flag
   
   CALL ReadVar ( UnIn, FileName, InitInp%PtfmYFChr, 'PtfmYFChr', 'Platform yaw rotation force flag', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read PtfmYF parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
   
   CALL Conv2UC( InitInp%PtfmYFChr )    ! Convert Line to upper case.
   
   
   !-------------------------------------------------------------------------------------------------
   ! Floating Platform Additional Stiffness and Damping Section
   !-------------------------------------------------------------------------------------------------
   
   
     ! Header
      
   CALL ReadCom( UnIn, FileName, 'Additional stiffness and damping header', ErrStat, ErrMsg, UnEchoLocal )
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read additional stiffness and damping header line.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
   
   
      ! AddF0 - Additional preload
      
   CALL ReadAry ( UnIn, FileName, InitInp%WAMIT%AddF0, 6, 'AddF0', &
                           ' Additional preload vector', ErrStat,  ErrMsg, UnEchoLocal )
   
      ! AddCLin
      
   DO I=1,6 
        
      WRITE(strI,'(I1)') I
      CALL ReadAry ( UnIn, FileName, InitInp%WAMIT%AddCLin(I,:), 6, 'AddCLin', &
                           ' Row '//strI//' of the additional linear stiffness matrix', ErrStat,  ErrMsg, UnEchoLocal )         
      
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Failed to read AddCLin parameter.'
         ErrStat = ErrID_Fatal
         CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      END IF 
   END DO
   
   
       ! AddBLin
      
   DO I=1,6 
        
      WRITE(strI,'(I1)') I
      CALL ReadAry ( UnIn, FileName, InitInp%WAMIT%AddBLin(I,:), 6, 'AddBLin', &
                           ' Row '//strI//' of the additional linear damping matrix', ErrStat,  ErrMsg, UnEchoLocal )         
      
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Failed to read AddBLin parameter.'
         ErrStat = ErrID_Fatal
         CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      END IF 
   END DO
   
   
       ! AddBQuad
      
   DO I=1,6 
        
      WRITE(strI,'(I1)') I
      CALL ReadAry ( UnIn, FileName, InitInp%WAMIT%AddBQuad(I,:), 6, 'AddBQuad', &
                           ' Row '//strI//' of the additional quadratic damping matrix', ErrStat,  ErrMsg, UnEchoLocal )         
      
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Failed to read AddBQuad parameter.'
         ErrStat = ErrID_Fatal
         CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      END IF 
   END DO
   
   
   !-------------------------------------------------------------------------------------------------
   !  Heave Coefficients Section
   !-------------------------------------------------------------------------------------------------
   
   
       ! Header
      
   CALL ReadCom( UnIn, FileName, 'Heave coefs header', ErrStat, ErrMsg, UnEchoLocal )
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read heave coefs header line.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
   
   
      ! NHvCoef - Number of heave coefficients
   
   CALL ReadVar ( UnIn, FileName, InitInp%Morison%NHvCoefs, 'NHvCoefs', 'Number of heave coefficients', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read NHvCoefs parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
   
         ! Table header
      
      CALL ReadCom( UnIn, FileName, 'Heave coefficient table header', ErrStat, ErrMsg, UnEchoLocal )
   
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Failed to read heave coefficient table header line.'
         ErrStat = ErrID_Fatal
         CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      END IF
      
         ! Table header
      
      CALL ReadCom( UnIn, FileName, 'Heave coefficient table header', ErrStat, ErrMsg, UnEchoLocal )
   
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Failed to read heave coefficient table header line.'
         ErrStat = ErrID_Fatal
         CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      END IF
   
   IF ( InitInp%Morison%NHvCoefs > 0 ) THEN
      
      
         ! Allocate memory for Heave Coef-related arrays
         
      ALLOCATE ( InitInp%Morison%HeaveCoefs(InitInp%Morison%NHvCoefs), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating space for HeaveCoefs array.'
         ErrStat = ErrID_Fatal
         CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      END IF    
          
      DO I = 1,InitInp%Morison%NHvCoefs
            ! read the table entries   HvCoefID   CdHv  CaHv    in the HydroDyn input file
         READ(UnIn,'(A)',IOSTAT=ErrStat) Line      !read into a line 
            
         IF (ErrStat == 0) THEN
            READ(Line,*,IOSTAT=ErrStat) InitInp%Morison%HeaveCoefs(I)%HvCoefID, InitInp%Morison%HeaveCoefs(I)%HvCd, InitInp%Morison%HeaveCoefs(I)%HvCa
         END IF      
       
         IF ( ErrStat /= ErrID_None ) THEN
            ErrMsg  = ' Failed to read heave coefficients.'
            ErrStat = ErrID_Fatal
            CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
            CLOSE( UnIn )
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
      
   CALL ReadCom( UnIn, FileName, 'Member joints header', ErrStat, ErrMsg, UnEchoLocal )
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Member joints header line.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
   
   
      ! NJoints - Number of member joints
   
   CALL ReadVar ( UnIn, FileName, InitInp%Morison%NJoints, 'NJoints', 'Number of member joints', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read NJoints parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
   
         ! Table header
      
      CALL ReadCom( UnIn, FileName, 'Member joints table header', ErrStat, ErrMsg, UnEchoLocal )
   
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Failed to read Member joints table header line.'
         ErrStat = ErrID_Fatal
         CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      END IF
      
         ! Table header
      
      CALL ReadCom( UnIn, FileName, 'Member joints table header', ErrStat, ErrMsg, UnEchoLocal )
   
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Failed to read Member joints table header line.'
         ErrStat = ErrID_Fatal
         CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      END IF
   
   IF ( InitInp%Morison%NJoints > 0 ) THEN
      
      
         ! Allocate memory for Joint-related arrays
         
      ALLOCATE ( InitInp%Morison%InpJoints(InitInp%Morison%NJoints), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating space for InpJoints array.'
         ErrStat = ErrID_Fatal
         CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      END IF    
          
      DO I = 1,InitInp%Morison%NJoints
            ! read the table entries   JointID   Jointxi     Jointyi    Jointzi      JointHvID   JointOvrlp    in the HydroDyn input file
         READ(UnIn,'(A)',IOSTAT=ErrStat) Line      !read into a line 
            
         IF (ErrStat == 0) THEN
            READ(Line,*,IOSTAT=ErrStat) InitInp%Morison%InpJoints(I)%JointID, InitInp%Morison%InpJoints(I)%JointPos(1), InitInp%Morison%InpJoints(I)%JointPos(2), InitInp%Morison%InpJoints(I)%JointPos(3), InitInp%Morison%InpJoints(I)%JointHvID, InitInp%Morison%InpJoints(I)%JointOvrlp
         END IF      
       
         IF ( ErrStat /= ErrID_None ) THEN
            ErrMsg  = ' Failed to read joints.'
            ErrStat = ErrID_Fatal
            CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
            CLOSE( UnIn )
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
      
   CALL ReadCom( UnIn, FileName, 'Member cross-section properties header', ErrStat, ErrMsg, UnEchoLocal )
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Member cross-section properties header line.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
   
   
      ! NPropSets - Number of member cross-section property sets
   
   CALL ReadVar ( UnIn, FileName, InitInp%Morison%NPropSets, 'NPropSets', 'Number of member cross-section property sets', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read NPropSets parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
   
      ! Table header
      
   CALL ReadCom( UnIn, FileName, 'Member cross-section properties table header', ErrStat, ErrMsg, UnEchoLocal )
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Member cross-section properties table header line.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
      
      ! Table header
      
   CALL ReadCom( UnIn, FileName, 'Member cross-section properties table header', ErrStat, ErrMsg, UnEchoLocal )
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Member cross-section properties table header line.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
   
   IF ( InitInp%Morison%NPropSets > 0 ) THEN
      
      
         ! Allocate memory for Member cross-section property set-related arrays
         
      ALLOCATE ( InitInp%Morison%MPropSets(InitInp%Morison%NPropSets), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating space for MPropSets array.'
         ErrStat = ErrID_Fatal
         CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      END IF
          
          
      DO I = 1,InitInp%Morison%NPropSets
         
         READ(UnIn,'(A)',IOSTAT=ErrStat) Line      !read into a line 

         IF (ErrStat == 0) THEN
            READ(Line,*,IOSTAT=ErrStat) InitInp%Morison%MPropSets(I)%PropSetID, InitInp%Morison%MPropSets(I)%PropD, InitInp%Morison%MPropSets(I)%PropThck
         END IF      
       
         IF ( ErrStat /= ErrID_None ) THEN
            ErrMsg  = ' Failed to read member cross-section properties.'
            ErrStat = ErrID_Fatal
            CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
            CLOSE( UnIn )
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
      
   CALL ReadCom( UnIn, FileName, 'Simple hydrodynamic coefficients header', ErrStat, ErrMsg, UnEchoLocal )
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Simple hydrodynamic coefficients header line.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
   
      ! Header
      
   CALL ReadCom( UnIn, FileName, 'Simple hydrodynamic coefficients table header', ErrStat, ErrMsg, UnEchoLocal )
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Simple hydrodynamic coefficients table header line.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
   
      ! Header
      
   CALL ReadCom( UnIn, FileName, 'Simple hydrodynamic coefficients table header', ErrStat, ErrMsg, UnEchoLocal )
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Simple hydrodynamic coefficients table header line.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
   
   
   READ(UnIn,'(A)',IOSTAT=ErrStat) Line      !read into a line 

   IF (ErrStat == 0) THEN
      READ(Line,*,IOSTAT=ErrStat) InitInp%Morison%SimplCd, InitInp%Morison%SimplCdMG, InitInp%Morison%SimplCa, InitInp%Morison%SimplCaMG
   END IF      
       
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read simple hydrodynamic coefficients.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF 
   
   IF ( InitInp%Echo ) THEN
      WRITE( UnEchoLocal, '(A)' ) TRIM(Line)
   END IF
   
   
   
   
   !-------------------------------------------------------------------------------------------------
   ! Depth-based Hydrodynamic Coefficients Section
   !-------------------------------------------------------------------------------------------------
   
   
      ! Header
      
   CALL ReadCom( UnIn, FileName, 'Depth-based hydrodynamic coefficients header', ErrStat, ErrMsg, UnEchoLocal )
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read depth-based hydrodynamic coefficients header line.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
   
   
      ! NCoefDpth - Number of depth-based hydrodynamic coefficient property sets
   
   CALL ReadVar ( UnIn, FileName, InitInp%Morison%NCoefDpth, 'NCoefDpth', 'Number of depth-based hydrodynamic coefficient property sets', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read NCoefDpth parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
   
      ! Table header
      
   CALL ReadCom( UnIn, FileName, 'Depth-based hydrodynamic coefficients table header', ErrStat, ErrMsg, UnEchoLocal )
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read depth-based hydrodynamic coefficients table header line.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
      
      ! Table header
      
   CALL ReadCom( UnIn, FileName, 'Depth-based hydrodynamic coefficients table header', ErrStat, ErrMsg, UnEchoLocal )
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read depth-based hydrodynamic coefficients table header line.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
   
   
   IF ( InitInp%Morison%NCoefDpth > 0 ) THEN
      
      
         ! Allocate memory for depth-based coefficient arrays
         
      ALLOCATE ( InitInp%Morison%CoefDpths(InitInp%Morison%NCoefDpth), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating space for CoefDpths array.'
         ErrStat = ErrID_Fatal
         CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      END IF
            
      DO I = 1,InitInp%Morison%NCoefDpth
         
         READ(UnIn,'(A)',IOSTAT=ErrStat) Line      !read into a line 

         IF (ErrStat == 0) THEN
            READ(Line,*,IOSTAT=ErrStat) InitInp%Morison%CoefDpths(I)%Dpth, InitInp%Morison%CoefDpths(I)%DpthCd, InitInp%Morison%CoefDpths(I)%DpthCdMG, InitInp%Morison%CoefDpths(I)%DpthCa, InitInp%Morison%CoefDpths(I)%DpthCaMG
         END IF      
       
         IF ( ErrStat /= ErrID_None ) THEN
            ErrMsg  = ' Failed to read member cross-section properties.'
            ErrStat = ErrID_Fatal
            CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
            CLOSE( UnIn )
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
      
   CALL ReadCom( UnIn, FileName, 'Member-based hydrodynamic coefficients header', ErrStat, ErrMsg, UnEchoLocal )
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read member-based hydrodynamic coefficients header line.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
   
   
      ! NCoefMembers - Number of member-based hydrodynamic coefficient property sets
   
   CALL ReadVar ( UnIn, FileName, InitInp%Morison%NCoefMembers, 'NCoefMembers', 'Number of member-based hydrodynamic coefficient property sets', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read NCoefMembers parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
   
      ! Table header
      
   CALL ReadCom( UnIn, FileName, 'Member-based hydrodynamic coefficients table header', ErrStat, ErrMsg, UnEchoLocal )
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read member-based hydrodynamic coefficients table header line.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
      
      ! Table header
      
   CALL ReadCom( UnIn, FileName, 'Member-based hydrodynamic coefficients table header', ErrStat, ErrMsg, UnEchoLocal )
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read member-based hydrodynamic coefficients table header line.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
   
   
   IF ( InitInp%Morison%NCoefMembers > 0 ) THEN
      
      
         ! Allocate memory for Member-based coefficient arrays
         
      ALLOCATE ( InitInp%Morison%CoefMembers(InitInp%Morison%NCoefMembers), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating space for CoefMembers array.'
         ErrStat = ErrID_Fatal
         CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      END IF
          
      DO I = 1,InitInp%Morison%NCoefMembers
         
         READ(UnIn,'(A)',IOSTAT=ErrStat) Line      !read into a line 

         IF (ErrStat == 0) THEN
            READ(Line,*,IOSTAT=ErrStat) InitInp%Morison%CoefMembers(I)%MemberID,    &
                                        InitInp%Morison%CoefMembers(I)%MemberCd1,   InitInp%Morison%CoefMembers(I)%MemberCd2,   &
                                        InitInp%Morison%CoefMembers(I)%MemberCdMG1, InitInp%Morison%CoefMembers(I)%MemberCdMG2, &
                                        InitInp%Morison%CoefMembers(I)%MemberCa1,   InitInp%Morison%CoefMembers(I)%MemberCa2,   &
                                        InitInp%Morison%CoefMembers(I)%MemberCaMG1, InitInp%Morison%CoefMembers(I)%MemberCaMG2
         END IF      
       
         IF ( ErrStat /= ErrID_None ) THEN
            ErrMsg  = ' Failed to read member cross-section properties.'
            ErrStat = ErrID_Fatal
            CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
            CLOSE( UnIn )
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
      
   CALL ReadCom( UnIn, FileName, 'Members header', ErrStat, ErrMsg, UnEchoLocal )
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read members header line.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
   
   
      ! NMembers - Number of members in the input file
   
   CALL ReadVar ( UnIn, FileName, InitInp%Morison%NMembers, 'NMembers', 'Number of members', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read NMembers parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
   
      ! Table header
      
   CALL ReadCom( UnIn, FileName, 'Members table header', ErrStat, ErrMsg, UnEchoLocal )
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read members table header line.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
      
      ! Table header
      
   CALL ReadCom( UnIn, FileName, 'Members table header', ErrStat, ErrMsg, UnEchoLocal )
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read members table header line.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
   
   
   IF ( InitInp%Morison%NMembers > 0 ) THEN
      
      
         ! Allocate memory for Members arrays
         
      ALLOCATE ( InitInp%Morison%InpMembers(InitInp%Morison%NMembers), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating space for InpMembers array.'
         ErrStat = ErrID_Fatal
         CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      END IF
          
      DO I = 1,InitInp%Morison%NMembers
         !ReadStr ( UnIn, Fil, Line, 'Joint table', VarDescr, ErrStat )
         READ(UnIn,'(A)',IOSTAT=ErrStat) Line      !read into a line 

         IF (ErrStat == 0) THEN
            READ(Line,*,IOSTAT=ErrStat) InitInp%Morison%InpMembers(I)%MemberID,    InitInp%Morison%InpMembers(I)%MJointID1,    &
                                        InitInp%Morison%InpMembers(I)%MJointID2,   InitInp%Morison%InpMembers(I)%MPropSetID1,  &
                                        InitInp%Morison%InpMembers(I)%MPropSetID2, InitInp%Morison%InpMembers(I)%MDivSize,     &  
                                        InitInp%Morison%InpMembers(I)%MCoefMod,    InitInp%Morison%InpMembers(I)%PropWAMIT
         END IF      
       
         IF ( ErrStat /= ErrID_None ) THEN
            ErrMsg  = ' Failed to read member properties.'
            ErrStat = ErrID_Fatal
            CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
            CLOSE( UnIn )
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
      
   CALL ReadCom( UnIn, FileName, 'Filled members header', ErrStat, ErrMsg, UnEchoLocal )
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read filled members header line.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
   
   
      ! NFillGroups - Number of fill groups
   
   CALL ReadVar ( UnIn, FileName, InitInp%Morison%NFillGroups, 'NFillGroups', 'Number of fill groups', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read NFillGroup parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
   
      ! Table header
      
   CALL ReadCom( UnIn, FileName, 'Fill groups table header', ErrStat, ErrMsg, UnEchoLocal )
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read fill groups table header line.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
      
      ! Table header
      
   CALL ReadCom( UnIn, FileName, 'Fill groups table header', ErrStat, ErrMsg, UnEchoLocal )
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read fill groups table header line.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
   
   
   IF ( InitInp%Morison%NFillGroups > 0 ) THEN
      
      
         ! Allocate memory for filled group arrays
         
      ALLOCATE ( InitInp%Morison%FilledGroups(InitInp%Morison%NFillGroups), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating space for FilledGroups array.'
         ErrStat = ErrID_Fatal
         CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      END IF
          
      DO I = 1,InitInp%Morison%NFillGroups
         
         READ(UnIn,'(A)',IOSTAT=ErrStat) Line      !read into a line 

         IF (ErrStat == 0) THEN
            
            READ(Line,*,IOSTAT=ErrStat) InitInp%Morison%FilledGroups(I)%FillNumM
            
            ALLOCATE ( InitInp%Morison%FilledGroups(I)%FillMList(InitInp%Morison%FilledGroups(I)%FillNumM), STAT = ErrStat )
            
            IF ( ErrStat /= ErrID_None ) THEN
               ErrMsg  = ' Error allocating space for FillMList array.'
               ErrStat = ErrID_Fatal
               CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
               CLOSE( UnIn )
               RETURN
            END IF
            
            READ(Line,*,IOSTAT=ErrStat) InitInp%Morison%FilledGroups(I)%FillNumM,  InitInp%Morison%FilledGroups(I)%FillMList,   &
                                        InitInp%Morison%FilledGroups(I)%FillFSLoc, InitInp%Morison%FilledGroups(I)%FillDensChr
             
            IF ( ErrStat /= ErrID_None ) THEN
               
               ErrMsg  = ' Failed to read filled group properties.'
               ErrStat = ErrID_Fatal
               CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
               CLOSE( UnIn )
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
      
   CALL ReadCom( UnIn, FileName, 'Marine growth by depth header', ErrStat, ErrMsg, UnEchoLocal )
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read marine growth by depth header line.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
   
   
      ! NMGDepths - Number marine growth depths
   
   CALL ReadVar ( UnIn, FileName, InitInp%Morison%NMGDepths, 'NMGDepths', 'Number marine growth depths', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read NMGDepths parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
   
      ! Table header
      
   CALL ReadCom( UnIn, FileName, 'Marine growth by depth table header', ErrStat, ErrMsg, UnEchoLocal )
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read marine growth by depth table header line.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
      
      ! Table header
      
   CALL ReadCom( UnIn, FileName, 'Marine growth by depth table header', ErrStat, ErrMsg, UnEchoLocal )
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read marine growth by depth table header line.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
   
   
   IF ( InitInp%Morison%NMGDepths > 0 ) THEN
      
      
         ! Allocate memory for marine growth depths array
         
      ALLOCATE ( InitInp%Morison%MGDepths(InitInp%Morison%NMGDepths), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating space for MGDepths array.'
         ErrStat = ErrID_Fatal
         CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      END IF
          
      DO I = 1,InitInp%Morison%NMGDepths
         
         READ(UnIn,'(A)',IOSTAT=ErrStat) Line      !read into a line 

         IF (ErrStat == 0) THEN
            READ(Line,*,IOSTAT=ErrStat) InitInp%Morison%MGDepths(I)%MGDpth, InitInp%Morison%MGDepths(I)%MGThck, InitInp%Morison%MGDepths(I)%MGDens     
         END IF      
       
         IF ( ErrStat /= ErrID_None ) THEN
            ErrMsg  = ' Failed to read marine growth depth properties.'
            ErrStat = ErrID_Fatal
            CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
            CLOSE( UnIn )
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
      
   CALL ReadCom( UnIn, FileName, 'Member output list header', ErrStat, ErrMsg, UnEchoLocal )
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read member output list header line.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
   
   
      ! NMOutputs - Number of members to output
   
   CALL ReadVar ( UnIn, FileName, InitInp%Morison%NMOutputs, 'NMOutputs', 'Number of members to output', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read NMOutputs parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
   
      ! Table header
      
   CALL ReadCom( UnIn, FileName, 'Member output list table header', ErrStat, ErrMsg, UnEchoLocal )
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read member output list table header line.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
      
      ! Table header
      
   CALL ReadCom( UnIn, FileName, 'Member output list table header', ErrStat, ErrMsg, UnEchoLocal )
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read member output list table header line.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
   
   
   IF ( InitInp%Morison%NMOutputs > 0 ) THEN
      
      
         ! Allocate memory for filled group arrays
         
      ALLOCATE ( InitInp%Morison%MOutLst(InitInp%Morison%NMOutputs), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating space for MOutLst array.'
         ErrStat = ErrID_Fatal
         CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      END IF
          
      DO I = 1,InitInp%Morison%NMOutputs
         
         READ(UnIn,'(A)',IOSTAT=ErrStat) Line      !read into a line 

         IF (ErrStat == 0) THEN
            
            READ(Line,*,IOSTAT=ErrStat) InitInp%Morison%MOutLst(I)%MemberID, InitInp%Morison%MOutLst(I)%NOutLoc
            
            ALLOCATE ( InitInp%Morison%MOutLst(I)%NodeLocs(InitInp%Morison%MOutLst(I)%NOutLoc), STAT = ErrStat )
            
            IF ( ErrStat /= ErrID_None ) THEN
               ErrMsg  = ' Error allocating space for NodeLocs array.'
               ErrStat = ErrID_Fatal
               CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
               CLOSE( UnIn )
               RETURN
            END IF
            
            ALLOCATE ( InitInp%Morison%MOutLst(I)%Marker1(InitInp%Morison%MOutLst(I)%NOutLoc), STAT = ErrStat )
            
            IF ( ErrStat /= ErrID_None ) THEN
               ErrMsg  = ' Error allocating space for Marker1 array.'
               ErrStat = ErrID_Fatal
               CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
               CLOSE( UnIn )
               RETURN
            END IF
            
            ALLOCATE ( InitInp%Morison%MOutLst(I)%Marker2(InitInp%Morison%MOutLst(I)%NOutLoc), STAT = ErrStat )
            
            IF ( ErrStat /= ErrID_None ) THEN
               ErrMsg  = ' Error allocating space for Marker2 array.'
               ErrStat = ErrID_Fatal
               CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
               CLOSE( UnIn )
               RETURN
            END IF
            
            ALLOCATE ( InitInp%Morison%MOutLst(I)%s(InitInp%Morison%MOutLst(I)%NOutLoc), STAT = ErrStat )
            
            IF ( ErrStat /= ErrID_None ) THEN
               ErrMsg  = ' Error allocating space for s array.'
               ErrStat = ErrID_Fatal
               CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
               CLOSE( UnIn )
               RETURN
            END IF
            
            READ(Line,*,IOSTAT=ErrStat) InitInp%Morison%MOutLst(I)%MemberID,  InitInp%Morison%MOutLst(I)%NOutLoc,  &
                                        InitInp%Morison%MOutLst(I)%NodeLocs
             
            IF ( ErrStat /= ErrID_None ) THEN
               
               ErrMsg  = ' Failed to read member output list properties.'
               ErrStat = ErrID_Fatal
               CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
               CLOSE( UnIn )
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
      
   CALL ReadCom( UnIn, FileName, 'Joint output list header', ErrStat, ErrMsg, UnEchoLocal )
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read joint output list header line.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
   
   
      ! NJOutputs - Number of joints to output
   
   CALL ReadVar ( UnIn, FileName, InitInp%Morison%NJOutputs, 'NJOutputs', 'Number of joints to output', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read NJOutputs parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
   
   IF ( InitInp%Morison%NJOutputs > 0 ) THEN
      
      ALLOCATE ( InitInp%Morison%JOutLst(InitInp%Morison%NJOutputs), STAT = ErrStat )
             
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating space for JOutLst data structures.'
         ErrStat = ErrID_Fatal
         CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      END IF
      
      ALLOCATE ( tmpArray(InitInp%Morison%NJOutputs), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating space for temporary array for Joint outputs.'
         ErrStat = ErrID_Fatal
         CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      END IF
      
      CALL ReadAry ( UnIn, FileName, tmpArray, InitInp%Morison%NJOutputs, 'JOutLst', 'Joint output list', ErrStat,  ErrMsg, UnEchoLocal )      
      
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error reading JOutLst array.'
         ErrStat = ErrID_Fatal
         CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      END IF
   
      DO I = 1,InitInp%Morison%NJOutputs
         
         InitInp%Morison%JOutLst(I)%JointID = tmpArray(I)
      
         
      END DO
         
      
      
      
      
   END IF
   
   !-------------------------------------------------------------------------------------------------
   ! Data section for OUTPUT
   !-------------------------------------------------------------------------------------------------

      ! Header
      
   CALL ReadCom( UnIn, FileName, 'Output header', ErrStat, ErrMsg, UnEchoLocal )
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Output header line.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
   
   
         ! HDSum - Whether or not to generate a summary file
   
   CALL ReadVar ( UnIn, FileName, InitInp%HDSum, 'HDSum', 'Generate a HydroDyn summary file', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read HDSum parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
   
   
         ! OutAll - Whether or not to output information for every member and joint
   
   CALL ReadVar ( UnIn, FileName, InitInp%OutAll, 'OutAll', 'Generate all member and joint outputs', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read OutAll parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
   
   
         ! OutSwtch - Specify how to write to an output file
   
   CALL ReadVar ( UnIn, FileName, InitInp%OutSwtch, 'OutSwtch', 'Specify how to write to an output file', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read OutSwtch parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF  
   
   
        ! OutFmt - Format for numerical outputs
   
   CALL ReadVar ( UnIn, FileName, InitInp%OutFmt, 'OutFmt', 'Format for numerical outputs', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read OutFmt parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
   
         ! OutSFmt - Format for output column headers
   
   CALL ReadVar ( UnIn, FileName, InitInp%OutSFmt, 'OutSFmt', 'Format for output column headers', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read OutSFmt parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
   
   
   !-------------------------------------------------------------------------------------------------
   ! Data section for FLOATING PLATFORM OUTPUTS
   !-------------------------------------------------------------------------------------------------

      ! Header
      
   CALL ReadCom( UnIn, FileName, 'Floating Platform Outputs header', ErrStat, ErrMsg, UnEchoLocal )
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Floating Platform Outputs header line.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
   
         ! OutList - list of requested parameters to output to a file

   CALL ReadOutputList ( UnIn, FileName, InitInp%WAMIT%OutList, InitInp%WAMIT%NumOuts, &
                                              'OutList', 'List of floating platform outputs requested', ErrStat, ErrMsg, UnEchoLocal )
   
   IF ( ErrStat /= 0 ) THEN
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
   
   
   !-------------------------------------------------------------------------------------------------
   ! Data section for MESH-BASED OUTPUTS
   !-------------------------------------------------------------------------------------------------

      ! Header
      
   CALL ReadCom( UnIn, FileName, 'Mesh-based Outputs header', ErrStat, ErrMsg, UnEchoLocal )
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Mesh-based Outputs header line.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
   
         ! OutList - list of requested parameters to output to a file

   CALL ReadOutputList ( UnIn, FileName, InitInp%Morison%OutList, InitInp%Morison%NumOuts, &
                                              'OutList', 'List of mesh-based outputs requested', ErrStat, ErrMsg, UnEchoLocal )
   
   IF ( ErrStat /= 0 ) THEN
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
   
   !-------------------------------------------------------------------------------------------------
   ! This is the end of the input file
   !-------------------------------------------------------------------------------------------------
   CLOSE ( UnIn )
   
   
      ! Cleanup the Echo file and global variables
      
   CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
   
   RETURN    
   
   
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
   CHARACTER(10), ALLOCATABLE                       :: tmpOutLst(:)         !
   LOGICAL                                          :: JointUsed
   REAL(ReKi)                                       :: l
   REAL(ReKi)                                       :: lvec(3)
   
      ! Initialize ErrStat
         
   ErrStat = ErrID_None         
   ErrMsg  = ""    
      
      
      ! WtrDens - Water density.
      
   IF ( InitInp%Waves%WtrDens < 0.0 )  THEN
      ErrMsg  = ' WtrDens must not be negative.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF
      
   
      ! WtrDpth - Water depth 
      
   IF ( InitInp%Morison%WtrDpth <= 0.0 )  THEN
      ErrMsg  = ' WtrDpth must be greater than zero.'
      ErrStat = ErrID_Fatal     
      RETURN
   END IF

   
      ! MSL2SWL - Mean sea level to still water level
   
   IF ( .NOT. EqualRealNos(InitInp%Morison%MSL2SWL, 0.0) ) THEN  !TODO  Alter this check when we support MSL2SWL
      ErrMsg  = ' MSL2SWL must be 0. Future versions of HydroDyn will once again support any value of MSL2SWL.'
      ErrStat = ErrID_Fatal         
      RETURN
   END IF
      
   
      ! WaveMod - Wave kinematics model switch.

   IF ( LEN_TRIM(InitInp%Waves%WaveModChr) > 1 ) THEN
         
      IF ( InitInp%Waves%WaveModChr(1:2) == '1P' )  THEN                     ! The user wants to specify the phase in place of a random phase

         READ (InitInp%Waves%WaveModChr(3:),*,IOSTAT=ErrStat)  InitInp%Waves%WavePhase
         IF ( ErrStat /= ErrID_None ) THEN
            CALL CheckIOS ( ErrStat, "", 'WavePhase', NumType, .TRUE. )         
            RETURN
         END IF
         InitInp%Waves%WaveMod   = 10                               ! Internally define WaveMod = 10 to mean regular waves with a specified (nonrandom) phase
         InitInp%Waves%WavePhase = InitInp%Waves%WavePhase*D2R     ! Convert the phase from degrees to radians

      ELSE                                               ! The user must have specified WaveMod incorrectly.
         ErrStat = ErrID_Fatal
      END IF
   
   ELSE
   
      READ( InitInp%Waves%WaveModChr, *, IOSTAT=ErrStat ) InitInp%Waves%WaveMod
      
      IF ( ErrStat /= ErrID_None ) THEN
         CALL CheckIOS ( ErrStat, "", 'WaveMod', NumType, .TRUE. )
         RETURN
      END IF

   END IF ! LEN_TRIM(InitInp%Waves%WaveModChr)
   

   IF ( InitInp%Waves%WaveMod == 5 ) THEN
      ErrMsg  = 'GH Bladed wave model is currently disabled in this release.  Future versions will support GH Bladed wave models.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF
   
   IF ( InitInp%Waves%WaveMod < 0 .OR. InitInp%Waves%WaveMod > 4 ) THEN
   
      IF ( InitInp%HasWAMIT  ) THEN
      
         ErrMsg  = ' WaveMod must be 0, 1, 1P#, 2, 3, or 4.'
         ErrStat = ErrID_Fatal
         
         RETURN
      
      ELSE IF ( ErrStat /= ErrID_None .OR. InitInp%Waves%WaveMod /= 5)  THEN
            
         ErrMsg  = ' WaveMod must be 0, 1, 1P#, 2, 3, 4, or 5.'
         ErrStat = ErrID_Fatal
         
         RETURN
         
      END IF
      
   END IF
   
         
         ! WaveStMod - Model switch for stretching incident wave kinematics to instantaneous free surface. 
         
         ! TODO: We are only implementing WaveStMod = 0 (No stretching) at this point in time. 1 Mar 2013 GJH
         
   IF ( InitInp%Waves%WaveStMod /= 0 ) THEN 
      ErrMsg  = ' WaveStMod must be 0. Future versions of HydroDyn will once again support other wave stretching models.'
      ErrStat = ErrID_Fatal         
      RETURN
   END IF
   
   !IF ( InitInp%HasWAMIT == .FALSE. .AND. InitInp%Waves%WaveMod > 0 ) THEN 
   !
   !   IF ( ( InitInp%Waves%WaveStMod /= 0 ) .AND. ( InitInp%Waves%WaveStMod /= 1 ) .AND. &
   !        ( InitInp%Waves%WaveStMod /= 2 ) .AND. ( InitInp%Waves%WaveStMod /= 3 ) )  THEN
   !      ErrMsg  = ' WaveStMod must be 0, 1, 2, or 3.'
   !      ErrStat = ErrID_Fatal
   !      
   !      RETURN
   !   END IF
   !
   !   IF ( ( InitInp%Waves%WaveStMod /= 3 ) .AND. ( InitInp%Waves%WaveMod == 5 ) )  THEN
   !      ErrMsg  = ' WaveStMod must be set to 3 when WaveMod is set to 5.'
   !      ErrStat = ErrID_Fatal
   !      
   !      RETURN
   !   END IF
   !   
   !ELSE !don't use this one
   !
   !      ! NOTE: Do not read in WaveStMod for floating platforms since it is
   !      !       inconsistent to use stretching (which is a nonlinear correction) for
   !      !       the viscous drag term in Morison's equation while not accounting for
   !      !       stretching in the diffraction and radiation problems (according to
   !      !       Paul Sclavounos, there are such corrections).  Instead, the viscous
   !      !       drag term from Morison's equation is computed by integrating up to
   !      !       the MSL, regardless of the instantaneous free surface elevation.
   !
   !   InitInp%Waves%WaveStMod = 0
   !
   !END IF      
   
   
      ! WaveTMax - Analysis time for incident wave calculations.
      
   IF ( InitInp%Waves%WaveMod == 0 )  THEN   ! .TRUE if we have incident waves.
      InitInp%Waves%WaveTMax = 0.0
      ! TODO: Issue warning if WaveTMax was not already 0.0 in this case.
   END IF   
   
   
      ! WaveDT - Time step for incident wave calculations   
      
   IF ( InitInp%Waves%WaveMod > 0 )  THEN   ! .TRUE if we have incident waves.
      
      IF ( InitInp%Waves%WaveDT <= 0.0 )  THEN
         ErrMsg  = ' WaveDT must be greater than zero.'
         ErrStat = ErrID_Fatal
         RETURN
      END IF
      
   ELSE
      
      InitInp%Waves%WaveDT = 0.0
      
   END IF
   
   
       ! WaveHs - Significant wave height
       
   IF ( ( InitInp%Waves%WaveMod /= 0 ) .AND. ( InitInp%Waves%WaveMod /= 4 ) .AND. ( InitInp%Waves%WaveMod /= 5 ) ) THEN   ! .TRUE. (when WaveMod = 1, 2, 3, or 10) if we have plane progressive (regular), JONSWAP/Pierson-Moskowitz spectrum (irregular) waves, or white-noise waves, but not user-defined or GH Bladed wave data.
     
      IF ( InitInp%Waves%WaveHs <= 0.0 )  THEN
         ErrMsg  = ' WaveHs must be greater than zero.'
         ErrStat = ErrID_Fatal        
         RETURN
      END IF

   ELSE
   
      InitInp%Waves%WaveHs = 0.0  
   
   END IF      
   
   
      ! WaveTp - Peak spectral period.

   IF ( ( InitInp%Waves%WaveMod /= 0 ) .AND. ( InitInp%Waves%WaveMod /= 4 ) .AND. ( InitInp%Waves%WaveMod /= 5 ) ) THEN   ! .TRUE. (when WaveMod = 1, 2, 3, or 10) if we have plane progressive (regular), JONSWAP/Pierson-Moskowitz spectrum (irregular) waves, or white-noise waves, but not user-defined or GH Bladed wave data.

      IF ( InitInp%Waves%WaveTp <= 0.0 )  THEN
         ErrMsg  = ' WaveTp must be greater than zero.'
         ErrStat = ErrID_Fatal
         
         RETURN
      END IF

   ELSE
   
      InitInp%Waves%WaveTp = 0.0   

   END IF  
   
   
       ! WavePkShp - Peak shape parameter.
       
   CALL Conv2UC( InitInp%Waves%WavePkShpChr )    ! Convert Line to upper case.
   
   IF ( InitInp%Waves%WaveMod == 2 ) THEN   ! .TRUE if we have JONSWAP/Pierson-Moskowitz spectrum (irregular) waves, but not GH Bladed wave data.   

      IF ( TRIM(InitInp%Waves%WavePkShpChr) == 'DEFAULT' )  THEN   ! .TRUE. when one wants to use the default value of the peak shape parameter, conditioned on significant wave height and peak spectral period.

         InitInp%Waves%WavePkShp = WavePkShpDefault ( InitInp%Waves%WaveHs, InitInp%Waves%WaveTp )

      ELSE                                   ! The input must have been specified numerically.

         READ (InitInp%Waves%WavePkShpChr,*,IOSTAT=ErrStat)  InitInp%Waves%WavePkShp
         
         IF ( ErrStat /= ErrID_None ) THEN
            CALL CheckIOS ( ErrStat, "", 'WavePkShp', NumType, .TRUE. )
            
            RETURN
         END IF

         IF ( ( InitInp%Waves%WavePkShp < 1.0 ) .OR. ( InitInp%Waves%WavePkShp > 7.0 ) )  THEN        
            ErrMsg  = ' WavePkShp must be greater than or equal to 1 and less than or equal to 7.'
            ErrStat = ErrID_Fatal
            
            RETURN
         END IF

      END IF

   ELSE

      InitInp%Waves%WavePkShp = 1.0
          
   END IF      
   
   
      ! WvLowCOff and WvHiCOff - Wave Cut-off frequency
    
   IF ( InitInp%Waves%WvLowCOff < 0 ) THEN
      ErrMsg  = ' WvLowCOff must be greater than or equal to zero.'
      ErrStat = ErrID_Fatal        
      RETURN
   END IF
   
      ! Threshold upper cut-off based on sampling rate
!bjj: changed WaveDT to ReKi in equation to avoid error about mixed data types in MIN() function.     
   InitInp%Waves%WvHiCOff =  MIN( Pi/REAL(InitInp%Waves%WaveDT,ReKi), InitInp%Waves%WvHiCOff ) 
   !TODO Issue warning if we changed WvHiCOff  GJH 7/24/13
   
   IF ( InitInp%Waves%WvLowCOff >= InitInp%Waves%WvHiCOff ) THEN
      ErrMsg  = ' WvLowCOff must be less than WvHiCOff.'
      ErrStat = ErrID_Fatal        
      RETURN
   END IF
   
   
      ! WaveDir - Wave heading direction.  
      
   IF ( ( InitInp%Waves%WaveMod > 0 ) .AND. ( InitInp%Waves%WaveMod /= 5 ) )  THEN   ! .TRUE if we have incident waves, but not GH Bladed wave data.
    
      IF ( ( InitInp%Waves%WaveDir <= -180.0 ) .OR. ( InitInp%Waves%WaveDir > 180.0 ) )  THEN
         ErrMsg  = ' WaveDir must be greater than -180 and less than or equal to 180.'
         ErrStat = ErrID_Fatal        
         RETURN
      END IF

   ELSE

      InitInp%Waves%WaveDir = 0.0

   END IF      
   
   
       ! WaveSeed(1), !WaveSeed(2)

   IF ( .NOT. ( ( InitInp%Waves%WaveMod > 0 ) .AND. ( InitInp%Waves%WaveMod /= 5 ) .AND. ( InitInp%Waves%WaveMod /= 10 ) ) ) THEN   !.TRUE. for plane progressive (regular) with random phase or irregular wave 

      DO I = 1,2

         InitInp%Waves%WaveSeed(I) = 0
   
      END DO !I
   
   END IF
   

      ! GHWvFile
      
   IF ( InitInp%Waves%WaveMod == 5 ) THEN      ! .TRUE if we are to use GH Bladed wave data.

      IF ( LEN_TRIM( InitInp%Waves%GHWvFile ) == 0 )  THEN      
         ErrMsg  = ' GHWvFile must not be an empty string.'
         ErrStat = ErrID_Fatal
         RETURN
      END IF

   ELSE !don't use this one

      InitInp%Waves%GHWvFile = ""
      
   END IF
   
   
      ! NWaveElev
      
   IF ( InitInp%Waves%NWaveElev < 0 ) THEN
       
      ErrMsg  = ' NWaveElev must not be negative.'
      ErrStat = ErrID_Fatal
      RETURN
         
   END IF
    
   
      ! CurrMod - Current profile model switch
      
   IF ( ( InitInp%Current%CurrMod /= 0 ) .AND. ( InitInp%Current%CurrMod /= 1 ) .AND. ( InitInp%Current%CurrMod /= 2 ) )  THEN
      ErrMsg  = ' CurrMod must be 0, 1, or 2.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF

   IF ( ( InitInp%Current%CurrMod /= 0 ) .AND. ( InitInp%Waves%WaveMod == 5 ) )  THEN
      ErrMsg  = ' CurrMod must be set to 0 when WaveMod is set to 5: GH Bladed wave data.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF
   
   
      ! CurrSSV0 - Sub-surface current velocity at still water level
         
   IF ( InitInp%Current%CurrMod == 1 )  THEN  ! .TRUE if we have standard current.
    
      IF ( InitInp%Current%CurrSSV0 < 0.0 )  THEN
         ErrMsg  = ' CurrSSV0 must not be less than zero.'
         ErrStat = ErrID_Fatal
         RETURN         
      END IF  
  
   ELSE

      InitInp%Current%CurrSSV0 = 0.0      
        
   END IF
   
   
      ! CurrSSDirChr - Sub-surface current heading direction
      
   IF ( InitInp%Current%CurrMod == 1 )  THEN  ! .TRUE if we have standard current.
      

      IF ( TRIM(InitInp%Current%CurrSSDirChr) == 'DEFAULT' )  THEN   ! .TRUE. when one wants to use the default value of codirectionality between sub-surface current and incident wave propogation heading directions.

         IF ( InitInp%Waves%WaveMod == 0 ) THEN
            ErrMsg  = ' CurrSSDir must not be set to ''DEFAULT'' when WaveMod is set to 0.'
            ErrStat = ErrID_Fatal
            RETURN         
         END IF  

         InitInp%Current%CurrSSDir = InitInp%Waves%WaveDir

      ELSE                                   ! The input must have been specified numerically.

         READ (InitInp%Current%CurrSSDirChr,*,IOSTAT=ErrStat)  InitInp%Current%CurrSSDir
         CALL CheckIOS ( ErrStat, "", 'CurrSSDir', NumType, .TRUE. )

         IF ( ErrStat /= ErrID_None ) THEN
            RETURN
         END IF

         IF ( ( InitInp%Current%CurrSSDir <= -180.0 ) .OR. ( InitInp%Current%CurrSSDir > 180.0 ) )  THEN 
            ErrMsg  = ' CurrSSDir must be greater than -180 and less than or equal to 180.' 
            ErrStat = ErrID_Fatal
            RETURN         
         END IF  

      END IF
  
  
   ELSE

      InitInp%Current%CurrSSDir = 0.0    
   
   END IF
   
   
      ! CurrNSRef - Near-surface current reference depth.
   
   IF ( InitInp%Current%CurrMod == 1 )  THEN  ! .TRUE if we have standard current.

      IF ( InitInp%Current%CurrNSRef <= 0.0 ) THEN
         ErrMsg  = ' CurrNSRef must be greater than zero.'
         ErrStat = ErrID_Fatal
         RETURN
      END IF

   ELSE
   
      InitInp%Current%CurrNSRef = 0.0

   END IF   
   
   
   
        ! CurrNSV0 - Near-surface current velocity at still water level.
   
   IF ( InitInp%Current%CurrMod == 1 )  THEN  ! .TRUE if we have standard current.
  
      IF ( InitInp%Current%CurrNSV0 < 0.0 ) THEN
         ErrMsg  = ' CurrNSV0 must not be less than zero.'
         ErrStat = ErrID_Fatal
         RETURN
      END IF
    
   ELSE
   
      InitInp%Current%CurrNSV0 = 0.0
     
   END IF

   
      ! CurrNSDir - Near-surface current heading direction.

   IF ( InitInp%Current%CurrMod == 1 )  THEN  ! .TRUE if we have standard current. 
  
      IF ( ( InitInp%Current%CurrNSDir <= -180.0 ) .OR. ( InitInp%Current%CurrNSDir > 180.0 ) )  THEN
         ErrMsg  = ' CurrNSDir must be greater than -180 and less than or equal to 180.'
         ErrStat = ErrID_Fatal
         RETURN
      END IF
           
   ELSE
   
      InitInp%Current%CurrNSDir = 0.0

   END IF
   
   
      ! CurrDIV - Depth-independent current velocity.

   IF ( InitInp%Current%CurrMod == 1 )  THEN  ! .TRUE if we have standard current.

      IF ( InitInp%Current%CurrDIV < 0.0 ) THEN
         ErrMsg  = ' CurrDIV must not be less than zero.'
         ErrStat = ErrID_Fatal
         RETURN
      END IF
         
   ELSE

      InitInp%Current%CurrDIV = 0.0
   
   END IF
   
   
      ! CurrDIDir - Depth-independent current heading direction.
   
   IF ( InitInp%Current%CurrMod == 1 )  THEN  ! .TRUE if we have standard current.

      IF ( ( InitInp%Current%CurrDIDir <= -180.0 ) .OR. ( InitInp%Current%CurrDIDir > 180.0 ) ) THEN
         ErrMsg  = ' CurrDIDir must be greater than -180 and less than or equal to 180.'
         ErrStat = ErrID_Fatal
         RETURN
      END IF
    
   ELSE

      InitInp%Current%CurrDIDir = 0.0

   END IF
   
       ! WAMITFile - Root name of WAMIT output files
 
   IF ( InitInp%HasWAMIT ) THEN
       IF ( LEN_TRIM( InitInp%WAMIT%WAMITFile ) == 0 ) THEN
         ErrMsg  = ' WAMITFile must not be an empty string.'
         ErrStat = ErrID_Fatal
         RETURN
      END IF
      
      
         ! if this is a relative path, let's make it relative to the location of the main input file
         
      IF ( PathIsRelative( InitInp%WAMIT%WAMITFile ) ) THEN      
         CALL GetPath( TRIM(InitInp%InputFile), TmpPath ) 
         InitInp%WAMIT%WAMITFile = TRIM(TmpPath)//TRIM(InitInp%WAMIT%WAMITFile)
      END IF
         

   ELSE         
   
      InitInp%WAMIT%WAMITFile = ""
 
      
   END IF


      ! WAMITULEN - WAMIT characteristic body length scale

   IF ( InitInp%HasWAMIT ) THEN
         
      IF ( InitInp%WAMIT%WAMITULEN < 0.0 ) THEN
         ErrMsg  = ' WAMITULEN must be positive.'
         ErrStat = ErrID_Fatal
         RETURN
      END IF
   ELSE
      
      InitInp%WAMIT%WAMITULEN = 1.0
   
   END IF
   
   
      ! PtfmVol0 - Displaced volume of water when the platform is in its undisplaced position

   IF ( InitInp%HasWAMIT ) THEN
   
      IF ( InitInp%WAMIT%PtfmVol0 < 0.0 ) THEN
         ErrMsg  = ' PtfmVol0 must not be negative.'
         ErrStat = ErrID_Fatal
         RETURN
      END IF
      
   ELSE         
   
      InitInp%WAMIT%PtfmVol0 = 0.0

   END IF
   
   
      ! RdtnTMax - Analysis time for wave radiation kernel calculations
      ! NOTE: Use RdtnTMax = 0.0 to eliminate wave radiation damping
      
   IF ( InitInp%HasWAMIT ) THEN
     
      IF ( InitInp%WAMIT%RdtnTMax < 0.0 ) THEN
         ErrMsg  = ' RdtnTMax must not be negative.'
         ErrStat = ErrID_Fatal
         RETURN
      END IF

   ELSE         
   
      InitInp%WAMIT%RdtnTMax = 0.0

   END IF  
   
       ! RdtnDT - Time step for wave radiation kernel calculations

   IF ( InitInp%HasWAMIT ) THEN

      IF ( InitInp%WAMIT%Conv_Rdtn%RdtnDT <= 0.0 ) THEN
         ErrMsg  = ' RdtnDT must be greater than zero.'
         ErrStat = ErrID_Fatal
         RETURN
      END IF

   ELSE         
   
      InitInp%WAMIT%Conv_Rdtn%RdtnDT = 0.0
      
   END IF
   
   !-------------------------------------------------------------------------------------------------
   ! Data section for Floating platform force flags
   !-------------------------------------------------------------------------------------------------
   
   ! If DEFAULT was requested, then the required value has already been set by the calling program
   IF ( TRIM(InitInp%PtfmSgFChr) /= 'DEFAULT' )  THEN
      
      READ (InitInp%PtfmSgFChr,*,IOSTAT=ErrStat)  InitInp%PtfmSgF
      
      CALL CheckIOS ( ErrStat, "", 'PtfmSgF', NumType, .TRUE. )

      IF ( ErrStat /= ErrID_None ) THEN
         RETURN
      END IF
      
   END IF
   
   IF ( TRIM(InitInp%PtfmSwFChr) /= 'DEFAULT' )  THEN
      
      READ (InitInp%PtfmSwFChr,*,IOSTAT=ErrStat)  InitInp%PtfmSwF
      
      CALL CheckIOS ( ErrStat, "", 'PtfmSwF', NumType, .TRUE. )

      IF ( ErrStat /= ErrID_None ) THEN
         RETURN
      END IF
      
   END IF
   
   IF ( TRIM(InitInp%PtfmHvFChr) /= 'DEFAULT' )  THEN
      
      READ (InitInp%PtfmHvFChr,*,IOSTAT=ErrStat)  InitInp%PtfmHvF
      
      CALL CheckIOS ( ErrStat, "", 'PtfmHvF', NumType, .TRUE. )

      IF ( ErrStat /= ErrID_None ) THEN
         RETURN
      END IF
      
   END IF
   
   IF ( TRIM(InitInp%PtfmRFChr) /= 'DEFAULT' )  THEN
      
      READ (InitInp%PtfmRFChr,*,IOSTAT=ErrStat)  InitInp%PtfmRF
      
      CALL CheckIOS ( ErrStat, "", 'PtfmRF', NumType, .TRUE. )

      IF ( ErrStat /= ErrID_None ) THEN
         RETURN
      END IF
      
   END IF
   
   IF ( TRIM(InitInp%PtfmPFChr) /= 'DEFAULT' )  THEN
      
      READ (InitInp%PtfmPFChr,*,IOSTAT=ErrStat)  InitInp%PtfmPF
      
      CALL CheckIOS ( ErrStat, "", 'PtfmPF', NumType, .TRUE. )

      IF ( ErrStat /= ErrID_None ) THEN
         RETURN
      END IF
      
   END IF
   
   IF ( TRIM(InitInp%PtfmYFChr) /= 'DEFAULT' )  THEN
      
      READ (InitInp%PtfmYFChr,*,IOSTAT=ErrStat)  InitInp%PtfmYF
      
      CALL CheckIOS ( ErrStat, "", 'PtfmYF', NumType, .TRUE. )

      IF ( ErrStat /= ErrID_None ) THEN
         RETURN
      END IF
      
   END IF
   
   
      ! Add checks that all platform DOF flags are true.  TODO:  Allow true or false once these have been implemented
      
   IF ( ( .NOT. InitInp%PtfmSgF ) .OR.  ( .NOT. InitInp%PtfmSwF ) .OR. ( .NOT. InitInp%PtfmHvF ) .OR. ( .NOT. InitInp%PtfmRF ) .OR. ( .NOT. InitInp%PtfmPF ) .OR. ( .NOT. InitInp%PtfmYF ) )THEN
      ErrMsg  = ' All platform DOF parameters must be set to TRUE.  Future versions of HydroDyn will support values of TRUE,  FALSE, or DEFAULT.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF 
   
   
   !-------------------------------------------------------------------------------------------------
   ! Member Joints Section
   !-------------------------------------------------------------------------------------------------
   
   IF ( InitInp%Morison%NJoints < 0 ) THEN
      ErrMsg  = ' NJoints parameter cannot be negative.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF
   
   IF ( InitInp%Morison%NJoints == 1 ) THEN
      ErrMsg  = ' NJoints parameter cannot be set to 1.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF 
   
   
     
      ! Check the heave coefs are >= 0 and IDs are unique
   IF ( InitInp%Morison%NHvCoefs > 0 ) THEN
   
      DO I = 1,InitInp%Morison%NHvCoefs 
         
         !IF (  .NOT. EqualRealNos(InitInp%Morison%HeaveCoefs(I)%HvCd, 0.0) ) THEN
         !   ErrMsg  = ' HvCd must be equal to zero.  Future versions will allow for non-zero heave coefficients.'
         !   ErrStat = ErrID_Fatal
         !   RETURN
         !END IF   
         !IF (  .NOT. EqualRealNos(InitInp%Morison%HeaveCoefs(I)%HvCa, 0.0) ) THEN
         !   ErrMsg  = ' HvCa must be equal to zero.  Future versions will allow for non-zero heave coefficients.'
         !   ErrStat = ErrID_Fatal
         !   RETURN
         !END IF   
         
         ! TODO: Once Heave Coefs are working remove the above checks and uncomment the checks below.  GJH 9/29/2013
         IF (  InitInp%Morison%HeaveCoefs(I)%HvCd < 0 ) THEN
            ErrMsg  = ' HvCd must be greater or equal to zero.'
            ErrStat = ErrID_Fatal
            RETURN
         END IF   
         IF (  InitInp%Morison%HeaveCoefs(I)%HvCa < 0 ) THEN
            ErrMsg  = ' HvCa must be greater or equal to zero.'
            ErrStat = ErrID_Fatal
            RETURN
         END IF   
         
            ! Make sure that the current HvCoefID is not used elsewhere in the table.
         DO J = I+1,InitInp%Morison%NHvCoefs
            IF ( InitInp%Morison%HeaveCoefs(I)%HvCoefID == InitInp%Morison%HeaveCoefs(J)%HvCoefID ) THEN
               ErrMsg = ' Duplicate HvCoefIDs were found in the Heave Coefficients table.'
               ErrStat = ErrID_Fatal
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
               ErrMsg = ' Duplicate JointIDs were found in the Member Joints table.'
               ErrStat = ErrID_Fatal
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
            ErrMsg  = ' Every JointID in the Joints table must appear once in the Members table.'
            ErrStat = ErrID_Fatal
            RETURN
         END IF  
   ! TODO : Implement Super member elements. GJH 7/24/13
   
         IF ( InitInp%Morison%InpJoints(I)%JointOvrlp /= 0  ) THEN
            ErrMsg  = ' JointOvrlp parameter must be set to 0.  Future versions of HydroDyn will support vales of 0 or 1.'
            ErrStat = ErrID_Fatal
            RETURN
         END IF   
         !IF ( ( InitInp%Morison%InpJoints(I)%JointOvrlp < 0 ) .OR. ( InitInp%Morison%InpJoints(I)%JointOvrlp > 1 ) ) THEN
         !   ErrMsg  = ' JointOvrlp parameter must be set to 0 or 1.'
         !   ErrStat = ErrID_Fatal
         !   RETURN
         !END IF   
         
            ! Make sure the heave coef id appears in the hv table
         IF ( InitInp%Morison%NHvCoefs > 0 ) THEN
            InitInp%Morison%InpJoints(I)%JointHvIDIndx = -1
            DO J = 1,InitInp%Morison%NHvCoefs         
               IF ( InitInp%Morison%InpJoints(I)%JointHvID == InitInp%Morison%HeaveCoefs(J)%HvCoefID ) &
                  InitInp%Morison%InpJoints(I)%JointHvIDIndx = J   
            END DO
            IF ( InitInp%Morison%InpJoints(I)%JointHvIDIndx == -1 ) THEN
               ErrMsg  = ' The specified JointHvID in the Joints Table does not appear in the Heave Coefficients table.'
               ErrStat = ErrID_Fatal
               RETURN
            END IF
         ELSE
            ! TODO: Issue error because we need at least one Heave coef table entry GJH  8/1/31
         END IF
      
      END DO
   END IF
  

   !-------------------------------------------------------------------------------------------------
   ! Member Cross-section Properties Section
   !-------------------------------------------------------------------------------------------------   
   
   IF ( InitInp%Morison%NPropSets < 0 ) THEN
      ErrMsg  = ' Number of member cross-section property sets must be greater than zero.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF
   
   IF ( InitInp%Morison%NPropSets > 0 ) THEN
               
      DO I = 1,InitInp%Morison%NPropSets
         
            ! Make sure that the current JointID is not used elsewhere in the table.
         DO J = I+1,InitInp%Morison%NPropSets
            IF ( InitInp%Morison%MPropSets(I)%PropSetID == InitInp%Morison%MPropSets(J)%PropSetID ) THEN
               ErrMsg = ' Duplicate PropSetIDs were found in the Member Cross-section Properties table.'
               ErrStat = ErrID_Fatal
               RETURN
            END IF
         END DO
         
         IF ( ( InitInp%Morison%MPropSets(I)%PropD < 0 ) .OR.  ( InitInp%Morison%MPropSets(I)%PropThck < 0 ) .OR. ( ( InitInp%Morison%MPropSets(I)%PropD - InitInp%Morison%MPropSets(I)%PropThck / 2.0 ) < 0) ) THEN
            ErrMsg  = ' PropD and PropThck must be greater than zero and (PropD - propThck/2 ) must be greater than zero.'
            ErrStat = ErrID_Fatal
            RETURN             
         END IF
      END DO
              
   END IF

   
   !-------------------------------------------------------------------------------------------------
   ! Simple hydrodynamic coefficients Section
   !-------------------------------------------------------------------------------------------------
   
   IF ( InitInp%Morison%SimplCd < 0 ) THEN 
      ErrMsg  = ' SimplCd must be greater or equal to zero.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF
   IF ( InitInp%Morison%SimplCdMG < 0 ) THEN 
      ErrMsg  = ' SimplCdMG must be greater or equal to zero.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF
   IF ( InitInp%Morison%SimplCa < 0 ) THEN
      ErrMsg  = ' SimplCa must be greater or equal to zero.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF
   IF ( InitInp%Morison%SimplCaMG < 0 ) THEN
      ErrMsg  = ' SimplCaMG must be greater or equal to zero.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF
   
   
   !-------------------------------------------------------------------------------------------------
   ! Depth-based Hydrodynamic Coefficients Section
   !-------------------------------------------------------------------------------------------------
  
   IF ( InitInp%Morison%NCoefDpth < 0 ) THEN
      ErrMsg  = ' NCoefDpth must be greater or equal to zero.'
      ErrStat = ErrID_Fatal
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
         END IF
         IF ( InitInp%Morison%CoefDpths(I)%Dpth > MaxDepth ) THEN
            MaxDepth = InitInp%Morison%CoefDpths(I)%Dpth
         END IF
            
            ! Make sure that the current Dpth is not used elsewhere in the table.
         DO J = I+1,InitInp%Morison%NCoefDpth
            IF ( EqualRealNos( InitInp%Morison%CoefDpths(I)%Dpth, InitInp%Morison%CoefDpths(J)%Dpth ) ) THEN  
               ErrMsg = ' Duplicate Dpths were found in the Depth-based Hydrodynamic Coefficients table.'
               ErrStat = ErrID_Fatal
               RETURN
            END IF
         END DO
         
         IF ( InitInp%Morison%CoefDpths(I)%DpthCd < 0 ) THEN
            ErrMsg  = ' In the Depth-based hydrodynamic coefficients table, DpthCd must be greater or equal to zero.'
            ErrStat = ErrID_Fatal
            RETURN
         END IF
         IF ( InitInp%Morison%CoefDpths(I)%DpthCdMG < 0 ) THEN 
            ErrMsg  = ' In the Depth-based hydrodynamic coefficients table, DpthCdMG must be greater or equal to zero.'
            ErrStat = ErrID_Fatal
            RETURN
         END IF
         IF ( InitInp%Morison%CoefDpths(I)%DpthCa < 0 ) THEN 
            ErrMsg  = ' In the Depth-based hydrodynamic coefficients table, DpthCa must be greater or equal to zero.'
            ErrStat = ErrID_Fatal
            RETURN
         END IF
         IF ( InitInp%Morison%CoefDpths(I)%DpthCaMG < 0 ) THEN
            ErrMsg  = ' In the Depth-based hydrodynamic coefficients table, DpthCaMG must be greater or equal to zero.'
            ErrStat = ErrID_Fatal
            RETURN
         END IF
         
      END DO
      
      ! TODO: Sort the table based on depth so that a linear interpolation can be easily performed between entries.
      
   END IF
   
   
   !-------------------------------------------------------------------------------------------------
   ! Member-based Hydrodynamic Coefficients Section
   !-------------------------------------------------------------------------------------------------

   IF ( InitInp%Morison%NCoefMembers < 0 ) THEN
      ErrMsg  = ' NCoefMembers must be greater or equal to zero.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF
   
   IF ( InitInp%Morison%NCoefMembers > 0 ) THEN
            
      DO I = 1,InitInp%Morison%NCoefMembers
         
            ! Make sure that the current MemberID is not used elsewhere in the table.
         DO J = I+1,InitInp%Morison%NCoefMembers
            IF ( InitInp%Morison%CoefMembers(I)%MemberID == InitInp%Morison%CoefMembers(J)%MemberID ) THEN
               ErrMsg = ' Duplicate MemberIDs were found in the Member-based Hydrodynamic coefficients table.'
               ErrStat = ErrID_Fatal
               RETURN
            END IF
         END DO        
         
         
         
         IF ( InitInp%Morison%CoefMembers(I)%MemberCd1 < 0 ) THEN 
            ErrMsg  = ' In the member-based hydrodynamic coefficients table, MemberCd1 must be greater or equal to zero.'
            ErrStat = ErrID_Fatal
            RETURN
         END IF
         IF ( InitInp%Morison%CoefMembers(I)%MemberCd2 < 0 ) THEN 
            ErrMsg  = ' In the member-based hydrodynamic coefficients table, MemberCd2 must be greater or equal to zero.'
            ErrStat = ErrID_Fatal
            RETURN
         END IF
         IF ( InitInp%Morison%CoefMembers(I)%MemberCdMG1 < 0 ) THEN 
            ErrMsg  = ' In the member-based hydrodynamic coefficients table, MemberCdMG1 must be greater or equal to zero.'
            ErrStat = ErrID_Fatal
            RETURN
         END IF
         IF ( InitInp%Morison%CoefMembers(I)%MemberCdMG2 < 0 ) THEN 
            ErrMsg  = ' In the member-based hydrodynamic coefficients table, MemberCdMG2 must be greater or equal to zero.'
            ErrStat = ErrID_Fatal
            RETURN
         END IF
         IF ( InitInp%Morison%CoefMembers(I)%MemberCa1 < 0 ) THEN 
            ErrMsg  = ' In the member-based hydrodynamic coefficients table, MemberCa1 must be greater or equal to zero.'
            ErrStat = ErrID_Fatal
            RETURN
         END IF
         IF ( InitInp%Morison%CoefMembers(I)%MemberCa2 < 0 ) THEN 
            ErrMsg  = ' In the member-based hydrodynamic coefficients table, MemberCa2 must be greater or equal to zero.'
            ErrStat = ErrID_Fatal
            RETURN
         END IF
         IF ( InitInp%Morison%CoefMembers(I)%MemberCaMG1 < 0 ) THEN 
            ErrMsg  = ' In the member-based hydrodynamic coefficients table, MemberCaMG1 must be greater or equal to zero.'
            ErrStat = ErrID_Fatal
            RETURN
         END IF
         IF ( InitInp%Morison%CoefMembers(I)%MemberCaMG2 < 0 ) THEN 
            ErrMsg  = ' In the member-based hydrodynamic coefficients table, MemberCaMG2 must be greater or equal to zero.'
            ErrStat = ErrID_Fatal
            RETURN
         END IF
         
      END DO
      
   END IF

   
   !-------------------------------------------------------------------------------------------------
   ! Members Section
   !-------------------------------------------------------------------------------------------------
   
   IF ( InitInp%Morison%NMembers < 0 ) THEN
      ErrMsg  = ' NMembers in the Members table must be greater or equal to zero.'
      ErrStat = ErrID_Fatal
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
      END DO
      
      DO I = 1,InitInp%Morison%NMembers
            
            ! Make sure that the current MemberID is not used elsewhere in the table.
         DO J = I+1,InitInp%Morison%NMembers
            IF ( InitInp%Morison%InpMembers(I)%MemberID == InitInp%Morison%InpMembers(J)%MemberID ) THEN
               ErrMsg = ' Duplicate MemberIDs were found in the Members table.'
               ErrStat = ErrID_Fatal
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
            ErrMsg  = ' JointID1 in the Members table does not appear in the Joints table.'
            ErrStat = ErrID_Fatal
         END IF
         IF ( InitInp%Morison%InpMembers(I)%MJointID2Indx == -1 ) THEN
            ErrMsg  = ' JointID2 in the Members table does not appear in the Joints table.'
            ErrStat = ErrID_Fatal
         END IF
         
            ! Make sure we do not have any zero length members
         lvec = InitInp%Morison%InpJoints(InitInp%Morison%InpMembers(I)%MJointID1Indx)%JointPos - InitInp%Morison%InpJoints(InitInp%Morison%InpMembers(I)%MJointID2Indx)%JointPos
         l = sqrt( lvec(1)*lvec(1) + lvec(2)*lvec(2) + lvec(3)*lvec(3) )
         IF ( EqualRealNos(0.0, l) ) THEN
            ErrMsg  = ' A member cannot have zero length.'
            ErrStat = ErrID_Fatal
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
            ErrMsg  = ' MPropSetID1 in the Members table does not appear in the Member cross-section properties table.'
            ErrStat = ErrID_Fatal
         END IF
         IF ( InitInp%Morison%InpMembers(I)%MPropSetID2Indx == -1 ) THEN
            ErrMsg  = ' MPropSetID2 in the Members table does not appear in the Member cross-section properties table.'
            ErrStat = ErrID_Fatal
         END IF 
     
         
         ! NOTE: We cannot test that MDivSize > MemberLength yet because there may be a joint overlap which is going to alter the final length of this member
         
         IF ( InitInp%Morison%InpMembers(I)%MDivSize <= 0 ) THEN
            ErrMsg  = ' MDivSize must be greater than zero.'
            ErrStat = ErrID_Fatal
         END IF
         
            
         IF ( ( InitInp%Morison%InpMembers(I)%MCoefMod /= 1 ) .AND. ( InitInp%Morison%InpMembers(I)%MCoefMod /= 2 ) .AND. ( InitInp%Morison%InpMembers(I)%MCoefMod /= 3 ) )  THEN
            ErrMsg  = ' MCoefMod must be 1, 2, or 3.'
            ErrStat = ErrID_Fatal
            RETURN
         END IF
         
         IF ( InitInp%Morison%InpMembers(I)%MCoefMod == 2 ) THEN
            IF ( InitInp%Morison%NCoefDpth == 0 ) THEN
               ErrMsg  = ' NCoefDpth must be greater than zero when a member is using a depth-based coefficient model.'
               ErrStat = ErrID_Fatal
               RETURN
            END IF
               ! We will not extrapolate depth-based coefficient values, so make sure that the depth-based table has values that are outside the depth range of this member 
               ! NOTE: This is actually potentially overly conservative because the final member may be shorter due to joint overlap handling.
            z1 = InitInp%Morison%InpJoints( InitInp%Morison%InpMembers(I)%MJointID1Indx )%JointPos(3)
            z2 = InitInp%Morison%InpJoints( InitInp%Morison%InpMembers(I)%MJointID2Indx )%JointPos(3)
            MinMembrDpth = min( z1, z2 )
            MaxMembrDpth = max( z1, z2 )
            IF ( ( MinMembrDpth < MinDepth ) .OR. ( MaxMembrDpth > MaxDepth ) ) THEN
               ErrMsg  = ' This member uses a depth-based coefficient model, but the member depth is outside the range of values provided in the depth-based hydrodynamic coefficients table.'
               ErrStat = ErrID_Fatal
               RETURN
            END IF
            
         END IF
         
         
         IF ( InitInp%Morison%InpMembers(I)%MCoefMod == 3 ) THEN
            IF ( InitInp%Morison%NCoefMembers == 0 ) THEN
               ErrMsg  = ' NCoefMembers must be greater than zero when a member is using a member-based coefficient model.'
               ErrStat = ErrID_Fatal
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
               ErrMsg = ' Could not locate the MemberID referenced in the Members table in the associated Member-based Hydrodynamic coefficients table.'
               ErrStat = ErrID_Fatal
               RETURN
            END IF
         END IF
         
         IF ( InitInp%Morison%InpMembers(I)%PropWAMIT .AND. .NOT. InitInp%HasWAMIT ) THEN
            ErrMsg = ' A member cannot have PropWAMIT set to TRUE if HasWAMIT is set to FALSE in the FLOATING PLATFORM section.'
            ErrStat = ErrID_Fatal
            RETURN
         END IF
         
         
         
      END DO
      
   END IF
   
   
   !-------------------------------------------------------------------------------------------------
   ! Filled Members Section
   !-------------------------------------------------------------------------------------------------
   
   IF ( InitInp%Morison%NFillGroups < 0 ) THEN
      ErrMsg  = ' NFillGroups in the Filled-members table must be greater or equal to zero.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF      

   IF ( InitInp%Morison%NFillGroups > 0 ) THEN
          
      DO I = 1,InitInp%Morison%NFillGroups
         
         IF ( InitInp%Morison%FilledGroups(I)%FillNumM < 1 ) THEN
            ErrMsg  = ' FillNumM in the Filled-members table must be greater than zero.'
            ErrStat = ErrID_Fatal
            RETURN
         END IF
            
         DO J = 1,InitInp%Morison%FilledGroups(I)%FillNumM
            
            DO K=1,InitInp%Morison%NMembers
               IF ( InitInp%Morison%FilledGroups(I)%FillMList(J) == InitInp%Morison%InpMembers(K)%MemberID ) THEN
                  FoundID = .TRUE.
                     ! Check to make sure this member is not already part of another fill group!
                  IF ( InitInp%Morison%InpMembers(K)%MmbrFilledIDIndx /= -1 ) THEN
                     ErrMsg  = ' A member cannot be a part of more than one fill group!'
                     ErrStat = ErrStat
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
      
            READ (InitInp%Morison%FilledGroups(I)%FillDensChr,*,IOSTAT=ErrStat)  InitInp%Morison%FilledGroups(I)%FillDens
      
            CALL CheckIOS ( ErrStat, "", 'FillDens', NumType, .TRUE. )

            IF ( ErrStat /= ErrID_None ) THEN
               RETURN
            END IF
         ELSE
            InitInp%Morison%FilledGroups(I)%FillDens = InitInp%Waves%WtrDens
         END IF

      END DO
      
   END IF
   
   
   !-------------------------------------------------------------------------------------------------
   ! Marine Growth by Depth Section
   !-------------------------------------------------------------------------------------------------
 
   IF ( InitInp%Morison%NMGDepths < 0 ) THEN
      ErrMsg  = ' NMGDepths in the Marine growth table must be greater or equal to zero.'
      ErrStat = ErrID_Fatal
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
         END IF
         
            ! Make sure that the current MGDpth is not used elsewhere in the table.
         DO J = I+1,InitInp%Morison%NMGDepths
            IF ( EqualRealNos( InitInp%Morison%MGDepths(I)%MGDpth, InitInp%Morison%MGDepths(J)%MGDpth ) ) THEN  
               ErrMsg = ' Duplicate MGDpth were found in the Marine Growth table.'
               ErrStat = ErrID_Fatal
               RETURN
            END IF
         END DO
         
         IF ( InitInp%Morison%MGDepths(I)%MGThck < 0 ) THEN
            ErrMsg  = ' MGThck in the Marine growth table must be greater or equal to zero.'
            ErrStat = ErrID_Fatal
            RETURN
         END IF  
         IF ( InitInp%Morison%MGDepths(I)%MGDens < 0 ) THEN
            ErrMsg  = ' MGDens in the Marine growth table must be greater or equal to zero.'
            ErrStat = ErrID_Fatal
            RETURN
         END IF       
         
      END DO
      
   END IF
   
   
   !-------------------------------------------------------------------------------------------------
   ! Member Output List Section
   !-------------------------------------------------------------------------------------------------

   IF ( ( InitInp%Morison%NMOutputs < 0 ) .OR. ( InitInp%Morison%NMOutputs > 9 ) ) THEN
      ErrMsg  = ' NMOutputs in the Member output list must be greater or equal to zero and less than 10.'
      ErrStat = ErrID_Fatal
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
            ErrMsg  = ' MemberID in the Member output list table does not appear in the Members table.'
            ErrStat = ErrID_Fatal
            RETURN
         END IF
         
         IF ( ( InitInp%Morison%MOutLst(I)%NOutLoc < 1 ) .OR. ( InitInp%Morison%MOutLst(I)%NOutLoc > 9) ) THEN
            ErrMsg  = ' NOutLoc in the Member output list must be greater than zero and less than 10.'
            ErrStat = ErrID_Fatal
            RETURN
         END IF 
         
         DO J = 1,InitInp%Morison%MOutLst(I)%NOutLoc
            IF ( ( InitInp%Morison%MOutLst(I)%NodeLocs(J) < 0.0 ) .OR. ( InitInp%Morison%MOutLst(I)%NodeLocs(J) > 1.0 ) ) THEN
               ErrMsg  = ' NodeLocs in the Member output list must be greater or equal to 0.0 and less than or equal to 1.0.'
               ErrStat = ErrID_Fatal
               RETURN
            END IF
         END DO
           
         
      END DO
      
   END IF
   
   !-------------------------------------------------------------------------------------------------
   ! Joint Output List Section
   !-------------------------------------------------------------------------------------------------
   
   !IF ( InitInp%Morison%NJOutputs /= 0 ) THEN  ! TODO Remove this check and add back the other checks once Joint Outputs are supported
   !ErrMsg  = ' NJOutputs in the Joint output list must be equal to zero.  Future versions of HydroDyn will support values greater or equal to zero and less than 10.'
   !   ErrStat = ErrID_Fatal
   !   RETURN
   !END IF   
   
   
   IF ( ( InitInp%Morison%NJOutputs < 0 ) .OR. ( InitInp%Morison%NMOutputs > 9 ) ) THEN
      ErrMsg  = ' NJOutputs in the Joint output list must be greater or equal to zero and less than 10.'
      ErrStat = ErrID_Fatal
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
            ErrMsg  = ' JointID in the Joint output list table does not appear in the Joints table.'
            ErrStat = ErrID_Fatal
            RETURN
         END IF
      END DO
   END IF
   
   !-------------------------------------------------------------------------------------------------
   ! Data section for OUTPUT
   !-------------------------------------------------------------------------------------------------


      ! OutAll - output all member and joint data
      
   IF ( InitInp%OutAll ) THEN    !TODO: Alter this check once OutAll is supported
         ErrMsg  = ' OutAll must be FALSE. Future versions of HydroDyn will once again support values of either TRUE or FALSE.'
         ErrStat = ErrID_Fatal
   END IF
   
   
      ! OutSwtch - output file switch
   
   IF ( InitInp%OutSwtch /= 1 .AND. InitInp%OutSwtch /= 2 .AND. InitInp%OutSwtch /= 3 ) THEN
      ErrMsg  = ' OutSwitch must be set to 1, 2, or 3.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF
   
   !InitInp%OutFmt
   !InitInp%OutSFmt
   
   
         ! OutList - list of requested parameters to output to a file

  
   
   !----------------------------------------------------------
   ! Mesh-related Output List
   !----------------------------------------------------------
   
   IF ( InitInp%Morison%NumOuts > 0 ) THEN
      
         ! Create an  output list for validated outputs
      ALLOCATE ( InitInp%Morison%ValidOutList(InitInp%Morison%NumOuts), STAT = ErrStat )
      IF ( ErrStat /= 0 ) THEN
         ErrMsg  = ' Error allocating valid output list array.'
         ErrStat = ErrID_Fatal
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
      InitInp%Current%WtrDpth    = InitInp%Morison%WtrDpth - InitInp%Morison%MSL2SWL ! Adjust for the MSL2SWL
      
      ! Waves
      InitInp%Waves%Gravity      = InitInp%Gravity
      InitInp%Waves%UnSum        = InitInp%UnSum
      InitInp%Waves%WtrDpth      = InitInp%Morison%WtrDpth - InitInp%Morison%MSL2SWL ! Adjust for the MSL2SWL
      ! WAMIT
      InitInp%WAMIT%WtrDens      = InitInp%Waves%WtrDens
      InitInp%WAMIT%WaveDir      = InitInp%Waves%WaveDir
      InitInp%WAMIT%WaveMod      = InitInp%Waves%WaveMod
      InitInp%WAMIT%OutAll       = InitInp%OutAll
      InitInp%WAMIT%HasWAMIT     = InitInp%HasWAMIT
      ! Morison
      InitInp%Morison%UnSum      = InitInp%UnSum
      InitInp%Morison%Gravity    = InitInp%Gravity
      InitInp%Morison%WtrDens    = InitInp%Waves%WtrDens
      InitInp%Morison%OutAll     = InitInp%OutAll
      
         ! Process the input geometry and generate the simulation mesh representation
      CALL Morison_ProcessMorisonGeometry( InitInp%Morison, ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) THEN  
         RETURN
      END IF
      
         ! Set the number and global Z locations for the X and Y components of the current velocities
      InitInp%Current%NMorisonNodes = InitInp%Morison%NNodes
      
      ALLOCATE ( InitInp%Current%MorisonNodezi(InitInp%Morison%NNodes), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating space for MorisonNodezi array.'
         ErrStat = ErrID_Fatal     
         RETURN
      END IF
      
      
         
         ! Establish the number and locations where the wave kinematics will be computed
      InitInp%Waves%NWaveKin0   = InitInp%Morison%NNodes                          ! Number of points where the incident wave kinematics will be computed (-)
      ALLOCATE ( InitInp%Waves%WaveKinxi0(InitInp%Waves%NWaveKin0), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating space for WaveKinxi0 array.'
         ErrStat = ErrID_Fatal
         
         RETURN
      END IF
      ALLOCATE ( InitInp%Waves%WaveKinyi0(InitInp%Waves%NWaveKin0), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating space for WaveKinyi0 array.'
         ErrStat = ErrID_Fatal
         
         RETURN
      END IF
      ALLOCATE ( InitInp%Waves%WaveKinzi0(InitInp%Waves%NWaveKin0), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating space for WaveKinzi0 array.'
         ErrStat = ErrID_Fatal
         
         RETURN
      END IF
      DO I=1,InitInp%Morison%NNodes
         InitInp%Waves%WaveKinxi0(I)      = InitInp%Morison%Nodes(I)%JointPos(1)                          ! xi-coordinates for points where the incident wave kinematics will be computed; 
         InitInp%Waves%WaveKinyi0(I)      = InitInp%Morison%Nodes(I)%JointPos(2)                          ! yi-coordinates for points where the incident wave kinematics will be computed; 
         InitInp%Waves%WaveKinzi0(I)      = InitInp%Morison%Nodes(I)%JointPos(3) - InitInp%Morison%MSL2SWL   ! zi-coordinates for points where the incident wave kinematics will be computed, adjusted to the mean see level (meters)     
         InitInp%Current%MorisonNodezi(I) = InitInp%Waves%WaveKinzi0(I)
      END DO
      
END SUBROUTINE HydroDynInput_ProcessInitData
   
END MODULE HydroDyn_Input
