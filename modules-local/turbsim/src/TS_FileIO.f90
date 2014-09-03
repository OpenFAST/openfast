!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2014  National Renewable Energy Laboratory
!
!    This file is part of TurbSim.
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
MODULE TS_FileIO

   USE                     NWTC_Library   

   use TS_Profiles
   use TSSubs
   use ts_errors
   use TS_RandNum
   
   IMPLICIT                NONE

CONTAINS

!=======================================================================
!> This subroutine reads parameters from the primary TurbSim input file.
!> It validates most of the meteorology data (because it is used to 
!> calculate default values later in the routine)
SUBROUTINE ReadInputFile(InFile, p, p_cohStr, OtherSt_RandNum, ErrStat, ErrMsg)


  ! This subroutine is used to read parameters from the input file.

   IMPLICIT             NONE

   CHARACTER(*),                 INTENT(IN)    :: InFile             !< name of the primary TurbSim input file
   TYPE(TurbSim_ParameterType),  INTENT(INOUT) ::  p                 !< TurbSim's parameters
   TYPE(CohStr_ParameterType),   INTENT(INOUT) :: p_CohStr           !< parameters for coherent structures
   TYPE(RandNum_OtherStateType), INTENT(INOUT) :: OtherSt_RandNum    !< other states for random numbers (next seed, etc)


   INTEGER(IntKi)  ,             INTENT(OUT)   :: ErrStat            !< allocation status
   CHARACTER(*) ,                INTENT(OUT)   :: ErrMsg             !< error message



      ! Local variables

   REAL(ReKi)                    :: InCVar     (2)                   ! Contains the coherence parameters (used for input)
   REAL(ReKi)                    :: tmp                              ! variable for estimating Ustar
   REAL(ReKi)                    :: TmpUary (3)                      !Temporary vector to store windSpeed(z) values
   REAL(ReKi)                    :: TmpUstar(3)                      !Temporary vector to store ustar(z) values
   REAL(ReKi)                    :: TmpUstarD                        !Temporary ustarD value
   REAL(ReKi)                    :: TmpZary (3)                      !Temporary vector to store height(z) values
   REAL(ReKi)                    :: TmpZLary(3)                      !Temporary vector to store zL(z) values
                              
   INTEGER                       :: TmpIndex                         ! Contains the index number when searching for substrings
   INTEGER                       :: UI                               ! I/O unit for input file
   INTEGER                       :: UnEc                             ! I/O unit for echo file
                              
   LOGICAL                       :: getPLExp                         ! Whether the power law exponent needs to be calculated
   LOGICAL                       :: Randomize                        ! Whether to randomize the coherent turbulence scaling
   LOGICAL                       :: UseDefault                       ! Whether or not to use a default value
   LOGICAL                       :: IsUnusedParameter                ! Whether or not this variable will be ignored
                              
   CHARACTER(200)                :: Line                             ! An input line
   CHARACTER(1)                  :: Line1                            ! The first character of an input line

   INTEGER(IntKi)                :: ErrStat2                         ! Temporary Error status
   INTEGER(IntKi)                :: I                                ! Loop counter (number of times file has been read)

   LOGICAL                       :: Echo                             ! Determines if an echo file should be written
   CHARACTER(MaxMsgLen)          :: ErrMsg2                          ! Temporary Error message
   CHARACTER(1024)               :: PriPath                          ! Path name of the primary file

   CHARACTER(1024)               :: UserFile   

      ! Initialize some variables:
   ErrStat = ErrID_None
   ErrMsg  = ""
      
   UnEc = -1
   Echo = .FALSE.   
   CALL GetPath( InFile, PriPath )     ! Input files will be relative to the path where the primary input file is located.
        
   
   !===============================================================================================================================
   ! Open input file
   !===============================================================================================================================

   CALL GetNewUnit( UI, ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')
      
   CALL OpenFInpFile( UI, InFile, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF

   CALL WrScr1(' Reading the input file "'//TRIM(InFile)//'".' )

   ! Read the lines up/including to the "Echo" simulation control variable
   ! If echo is FALSE, don't write these lines to the echo file. 
   ! If Echo is TRUE, rewind and write on the second try.
   
   I = 1 !set the number of times we've read the file
   DO 
   !-------------------------- HEADER ---------------------------------------------
   
   
      !===============================================================================================================================
      ! Read the runtime options.
      !===============================================================================================================================

      CALL ReadCom( UI, InFile, "File Heading Line 1",    ErrStat2, ErrMsg2, UnEc)
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')
      
      CALL ReadCom( UI, InFile, "File Heading Line 2",    ErrStat2, ErrMsg2, UnEc)
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')
      
      CALL ReadCom( UI, InFile, "Runtime Options Heading",ErrStat2, ErrMsg2, UnEc)
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')
                  
      CALL ReadVar( UI, InFile, Echo, 'Echo', 'Echo switch', ErrStat2, ErrMsg2, UnEc )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')
         IF (ErrStat >= AbortErrLev) THEN
            CALL Cleanup()
            RETURN
         END IF   
         
         
      IF (.NOT. Echo .OR. I > 1) EXIT !exit this loop
   
         ! Otherwise, open the echo file, then rewind the input file and echo everything we've read
      
      I = I + 1         ! make sure we do this only once (increment counter that says how many times we've read this file)
      
      
      CALL OpenEcho ( UnEc, TRIM(p%RootName)//'.ech', ErrStat2, ErrMsg2, TurbSim_Ver )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF   
      
      IF ( UnEc > 0 )  WRITE (UnEc,'(/,A,/)')  'Data from '//TRIM(TurbSim_Ver%Name)//' primary input file "'//TRIM( InFile )//'":'
      
      REWIND( UI, IOSTAT=ErrStat2 )  
         IF (ErrStat2 /= 0_IntKi ) THEN
            CALL SetErrStat( ErrID_Fatal, 'Error rewinding file "'//TRIM(InFile)//'".', ErrStat, ErrMsg, 'ReadInputFile')    
            RETURN
         END IF         
                  
   END DO
   
         
   
      ! RandSeed(1)
   CALL ReadVar( UI, InFile, p%RNG%RandSeed(1), "RandSeed(1)", "Random seed #1",ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')

      ! RandSeed(2)
   CALL ReadVar( UI, InFile, Line, "RandSeed(2)", "Random seed #2",ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF

         !  Check if alternate random number generator is to be used >>>>>>>>>>>>>>>>

      READ (Line,*,IOSTAT=ErrStat2) Line1  ! check the first character to make sure we don't have T/F, which can be interpreted as 1/-1 or 0 in Fortran

      CALL Conv2UC( Line1 )
      IF ( (Line1 == 'T') .OR.  (Line1 == 'F') ) THEN
         CALL SetErrStat( ErrID_Fatal, ' RandSeed(2): Invalid RNG type.', ErrStat, ErrMsg, 'ReadInputFile')
         CALL Cleanup()
         RETURN
      ENDIF

      READ (Line,*,IOSTAT=ErrStat2)  p%RNG%RandSeed(2)

      IF (ErrStat2 == 0) THEN ! the user entered a number
         p%RNG%RNG_type = "NORMAL"
         p%RNG%pRNG = pRNG_INTRINSIC
      ELSE

         p%RNG%RNG_type = ADJUSTL( Line )
         CALL Conv2UC( p%RNG%RNG_type )

         IF ( p%RNG%RNG_type == "RANLUX") THEN
            p%RNG%pRNG = pRNG_RANLUX
         ELSE IF ( p%RNG%RNG_type == "RNSNLW") THEN
            p%RNG%pRNG = pRNG_SNLW3
         ELSE         
            CALL SetErrStat( ErrID_Fatal, ' RandSeed(2): Invalid alternative random number generator.', ErrStat, ErrMsg, 'ReadInputFile')
            CALL Cleanup()
            RETURN
         ENDIF

      ENDIF
      
      !<<<<<<<<<<<<<<<<<<<<<< end rng


      ! --------- Read the flag for writing the binary HH (GenPro) turbulence parameters. -------------
   CALL ReadVar( UI, InFile, p%WrFile(FileExt_BIN), "WrBHHTP", "Output binary HH turbulence parameters? [RootName.bin]",ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')

      ! --------- Read the flag for writing the formatted turbulence parameters. ----------------------
   CALL ReadVar( UI, InFile, p%WrFile(FileExt_DAT), "WrFHHTP", "Output formatted turbulence parameters? [RootName.dat]",ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')

      ! ---------- Read the flag for writing the AeroDyn HH files. -------------------------------------
   CALL ReadVar( UI, InFile, p%WrFile(FileExt_HH), "WrADHH", "Output AeroDyn HH files? [RootName.hh]",ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')

      ! ---------- Read the flag for writing the AeroDyn FF files. ---------------------------------------
   CALL ReadVar( UI, InFile, p%WrFile(FileExt_BTS), "WrADFF", "Output AeroDyn FF files? [RootName.bts]",ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')

      ! ---------- Read the flag for writing the BLADED FF files. -----------------------------------------
   CALL ReadVar( UI, InFile, p%WrFile(FileExt_WND) , "WrBLFF", "Output BLADED FF files? [RootName.wnd]",ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')

      ! ---------- Read the flag for writing the AeroDyn tower files. --------------------------------------
   CALL ReadVar( UI, InFile, p%WrFile(FileExt_TWR), "WrADTWR", "Output tower data? [RootName.twr]",ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')

      ! ---------- Read the flag for writing the formatted FF files. ---------------------------------------
   CALL ReadVar( UI, InFile, p%WrFile(FileExt_UVW), "WrFMTFF", "Output formatted FF files? [RootName.u, .v, .w]",ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')

      ! ---------- Read the flag for writing coherent time series files. --------------------------------------
   CALL ReadVar( UI, InFile, p%WrFile(FileExt_CTS), "WrACT", "Output coherent time series files? [RootName.cts]",ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')

      ! ---------- Read the flag for turbine rotation. -----------------------------------------------------------
   CALL ReadVar( UI, InFile, p%grid%Clockwise, "Clockwise", "Clockwise rotation when looking downwind?",ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')

      ! ---------- Read the flag for determining IEC scaling -----------------------------------------------------
   CALL ReadVar( UI, InFile, p%IEC%ScaleIEC, "ScaleIEC, the switch for scaling IEC turbulence", &
                  "Scale IEC turbulence models to specified standard deviation?",ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')

      
      ! we'll check the errors before going to the next section of the input file
   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF

   !===============================================================================================================================
   ! Read the turbine/model specifications.
   !===============================================================================================================================

   CALL ReadCom( UI, InFile, "Turbine/Model Specifications Heading Line 1",ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')
      
   CALL ReadCom( UI, InFile, "Turbine/Model Specifications Heading Line 2",ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')

      ! ------------ Read in the vertical matrix dimension. ---------------------------------------------
   CALL ReadVar( UI, InFile, p%grid%NumGrid_Z, "NumGrid_Z", "Vertical grid-point matrix dimension [-]",ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')

      ! ------------ Read in the lateral matrix dimension. ---------------------------------------------
   CALL ReadVar( UI, InFile, p%grid%NumGrid_Y, "NumGrid_Y", "Horizontal grid-point matrix dimension [-]",ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')

      ! ------------ Read in the time step. ---------------------------------------------
   CALL ReadVar( UI, InFile, p%grid%TimeStep, "TimeStep", "Time step [seconds]",ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')

      ! ------------ Read in the analysis time. ---------------------------------------------
   CALL ReadVar( UI, InFile, p%grid%AnalysisTime, "AnalysisTime", "Analysis time [seconds]",ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')

      ! ------------ Read in the usable time. ---------------------------------------------
   CALL ReadVar( UI, InFile, Line, "UsableTime", "Usable output time [seconds]",ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! Check if usable time is "ALL" (for periodic files) >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
         READ( Line, *, IOSTAT=ErrStat2) p%grid%UsableTime
   
         IF ( ErrStat2 /= 0 ) THEN ! Line didn't contian a number
            CALL Conv2UC( Line )
            IF ( TRIM(Line) == 'ALL' ) THEN
               p%grid%Periodic   = .TRUE.
               p%grid%UsableTime = p%grid%AnalysisTime
            ELSE
               CALL SetErrStat( ErrID_Fatal, 'The usable output time must be a number greater than zero (or the string "ALL").', &
                                ErrStat, ErrMsg, 'ReadInputFile' )
               CALL Cleanup()
               RETURN
            END IF
         ELSE
            p%grid%Periodic = .FALSE.
         END IF
      ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< end check for UsableTime = "ALL" (periodic)

      ! ------------ Read in the hub height. ---------------------------------------------
   CALL ReadVar( UI, InFile, p%grid%HubHt, "HubHt", "Hub height [m]",ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')

      ! ------------ Read in the grid height. ---------------------------------------------
   CALL ReadVar( UI, InFile, p%grid%GridHeight, "GridHeight", "Grid height [m]",ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')

      ! ------------ Read in the grid width. ---------------------------------------------
   CALL ReadVar( UI, InFile, p%grid%GridWidth, "GridWidth", "Grid width [m]",ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')

      ! ------------ Read in the vertical flow angle. ---------------------------------------------
   CALL ReadVar( UI, InFile, p%met%VFlowAng, "VFlowAng", "Vertical flow angle [degrees]",ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')

      ! ------------ Read in the horizontal flow angle. ---------------------------------------------
   CALL ReadVar( UI, InFile, p%met%HFlowAng, "HFlowAng", "Horizontal flow angle [degrees]",ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')


!..................................................................................................................................
!  Do some error checking on the runtime options and turbine/model specifications before we read the meteorology data
!..................................................................................................................................

   IF ( p%IEC%ScaleIEC > 2 .OR. p%IEC%ScaleIEC < 0 ) CALL SetErrStat( ErrID_Fatal, 'The value for parameter ScaleIEC must be 0, 1, or 2.',    ErrStat, ErrMsg, 'ReadInputFile')
   IF ( p%grid%NumGrid_Z < 2 )                       CALL SetErrStat( ErrID_Fatal, 'The matrix must be >= 2x2.',                              ErrStat, ErrMsg, 'ReadInputFile')
   IF ( p%grid%NumGrid_Y < 2 )                       CALL SetErrStat( ErrID_Fatal, 'The matrix must be >= 2x2.',                              ErrStat, ErrMsg, 'ReadInputFile') 
   IF ( 0.5*p%grid%GridHeight > p%grid%HubHt  )      CALL SetErrStat( ErrID_Fatal, 'The hub must be higher than half of the grid height.',    ErrStat, ErrMsg, 'ReadInputFile')
   IF ( p%grid%GridWidth <=  0.0_ReKi )              CALL SetErrStat( ErrID_Fatal, 'The grid width must be greater than zero.',               ErrStat, ErrMsg, 'ReadInputFile')
   IF ( p%grid%HubHt <=  0.0 )                       CALL SetErrStat( ErrID_Fatal, 'The hub height must be greater than zero.',               ErrStat, ErrMsg, 'ReadInputFile')
   IF ( p%grid%AnalysisTime <=  0.0 )                CALL SetErrStat( ErrID_Fatal, 'The analysis time must be greater than zero.',            ErrStat, ErrMsg, 'ReadInputFile')
   IF ( p%grid%TimeStep <=  0.0 )                    CALL SetErrStat( ErrID_Fatal, 'The time step must be greater than zero.',                ErrStat, ErrMsg, 'ReadInputFile')
   IF ( ABS( p%met%VFlowAng ) > 45.0 )               CALL SetErrStat( ErrID_Fatal, 'The vertical flow angle must not exceed +/- 45 degrees.', ErrStat, ErrMsg, 'ReadInputFile')
   IF ( p%grid%UsableTime <=  0.0 )                  CALL SetErrStat( ErrID_Fatal, 'The usable output time must be a number greater than zero'&
                                                                                                                   //' or the string "ALL".', ErrStat, ErrMsg, 'ReadInputFile')
      
!..................................................................................................................................
!  initialize secondary parameters that will be used to calculate default values in the meteorological boundary conditions section
!..................................................................................................................................
      ! Initialize the RNG (for computing "default" values that contain random variates)
   CALL RandNum_Init(p%RNG, OtherSt_RandNum, ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')
            
   ! ***** Calculate the diameter of the rotor disk *****
   p%grid%RotorDiameter = MIN( p%grid%GridWidth, p%grid%GridHeight )
   
      ! we'll check the errors before going to the next section of the input file
   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF

   !===============================================================================================================================
   ! Read the meteorological boundary conditions.
   !===============================================================================================================================
   
   CALL ReadCom( UI, InFile, "Meteorological Boundary Conditions Heading Line 1",ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')
      
   CALL ReadCom( UI, InFile, "Meteorological Boundary Conditions Heading Line 2",ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')

      ! ------------ Read in the turbulence model. ---------------------------------------------
   CALL ReadVar( UI, InFile, p%met%TurbModel, "TurbModel", "spectral model",ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')
      
      ! ------------ Read in the UserFile------------------- ---------------------------------------------
   CALL ReadVar( UI, InFile, UserFile, "UserFile", "Name of the input file for user-defined spectra or time-series inputs",ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')
   !UserFile = "UsrSpec.inp"
   IF ( PathIsRelative( UserFile ) ) UserFile = TRIM(PriPath)//TRIM(UserFile)

   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF   
   
         ! Verify turbulence model is valid (used for default values later) and read supplemental files 
         !  for user-defined spectra or time-series >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      p%met%TurbModel = ADJUSTL( p%met%TurbModel )
      CALL Conv2UC( p%met%TurbModel )

      p%met%IsIECModel = .FALSE.
      p%usr%nPoints = 0
      SELECT CASE ( TRIM(p%met%TurbModel) )
         CASE ( 'IECKAI' )
            p%met%TMName = 'IEC Kaimal'
            p%met%TurbModel_ID = SpecModel_IECKAI
            p%met%IsIECModel = .TRUE.
         CASE ( 'IECVKM' )
            p%met%TMName = 'IEC von Karman'
            p%met%TurbModel_ID = SpecModel_IECVKM
            p%met%IsIECModel = .TRUE.
         CASE ( 'TIDAL' )
            p%met%TMName = 'Tidal Channel Turbulence'
            p%met%TurbModel_ID = SpecModel_TIDAL
         CASE ( 'RIVER' )
            p%met%TMName = 'River Turbulence'
            p%met%TurbModel_ID = SpecModel_RIVER
         CASE ( 'SMOOTH' )
            p%met%TMName = 'RISO Smooth Terrain'
            p%met%TurbModel_ID = SpecModel_SMOOTH
         CASE ( 'WF_UPW' )
            p%met%TMName = 'NREL Wind Farm Upwind'
            p%met%TurbModel_ID = SpecModel_WF_UPW
         CASE ( 'WF_07D' )
            p%met%TMName = 'NREL 7D Spacing Wind Farm'
            p%met%TurbModel_ID = SpecModel_WF_07D
         CASE ( 'WF_14D' )
            p%met%TMName = 'NREL 14D Spacing Wind Farm'
            p%met%TurbModel_ID = SpecModel_WF_14D
         CASE ( 'NONE'   )
            p%met%TMName = 'Steady wind components'
            p%met%TurbModel_ID = SpecModel_NONE
         CASE ( 'MODVKM' )
            p%met%TMName = 'Modified von Karman'
            p%met%TurbModel_ID = SpecModel_MODVKM
            p%met%IsIECModel = .TRUE.
         CASE ( 'API' )
            p%met%TMName = 'API'
            p%met%TurbModel_ID = SpecModel_API
            p%met%IsIECModel = .TRUE.
         CASE ( 'NWTCUP' )
            p%met%TMName = 'NREL National Wind Technology Center'
            p%met%TurbModel_ID = SpecModel_NWTCUP
         CASE ( 'GP_LLJ' )
            p%met%TMName = 'Great Plains Low-Level Jet'
            p%met%TurbModel_ID = SpecModel_GP_LLJ
         CASE ( 'USRVKM' )
            p%met%TMName = 'von Karman model with user-defined specifications'
            p%met%TurbModel_ID = SpecModel_USRVKM
         CASE ( 'USRINP' )
            p%met%TMName = 'User-input uniform spectra'
            p%met%TurbModel_ID = SpecModel_USER
            
            CALL GetUSRspec(UserFile, p, UnEc, ErrStat2, ErrMsg2)
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')
            
         CASE ( 'TIMESR' )
            p%met%TMName = 'User-input time series'
            p%met%TurbModel_ID = SpecModel_TimeSer 
            
            CALL ReadUSRTimeSeries(UserFile, p, UnEc, ErrStat2, ErrMsg2)
               CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat,ErrMsg, 'ReadInputFile')
               IF (ErrStat >= AbortErrLev) THEN
                  CALL Cleanup()
                  RETURN
               END IF
               
            CALL TimeSeriesToSpectra( p, ErrStat2, ErrMsg2 )
               CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat,ErrMsg, 'ReadInputFile')

print *, 'check the default values of all the rest of the inputs for the TIMESR model!!!'                  
         CASE DEFAULT
   !BONNIE: todo: add the UsrVKM model to this list when the model is complete 
            CALL SetErrStat( ErrID_Fatal, 'The turbulence model must be one of the following: "IECKAI", "IECVKM", "SMOOTH",' &
                       //' "WF_UPW", "WF_07D", "WF_14D", "NWTCUP", "GP_LLJ", "TIDAL", "RIVER", "API", "USRINP", "TIMESR" "NONE".', ErrStat, ErrMsg, 'ReadInputFile')
            CALL Cleanup()
            RETURN

         END SELECT  ! TurbModel
      
         ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< end TurbModel verification
         
!bjj: todo: verify that the API model sets the parameters for IECKAI as well (because it's using IECKAI for the v and w components)

      ! ------------ Read in the IEC standard and edition numbers. ---------------------------------------------
   CALL ReadVar( UI, InFile, Line, "IECstandard", "Number of the IEC standard",ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! Process this line for IEC standard & edition & IECeditionStr >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      CALL ProcessLine_IECstandard( Line, p%met%IsIECModel, p%met%TurbModel_ID, p%IEC%IECstandard, p%IEC%IECedition, p%IEC%IECeditionStr, ErrStat2, ErrMsg2)
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')
         IF (ErrStat >= AbortErrLev) THEN
            CALL Cleanup()
            RETURN
         END IF         
      ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< end processing of IECstandard input variable

      ! ------------ Read in the IEC turbulence characteristic. ---------------------------------------------
   CALL ReadVar( UI, InFile, Line, "IECturbc", "IEC turbulence characteristic",ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')

      ! Process this line for NumTurbInp, IECPerTurbInt, IECTurbC, and KHtest >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      CALL ProcessLine_IECTurbc(Line, p%met%IsIECModel, p%IEC%IECstandard, p%IEC%IECedition, p%IEC%IECeditionStr, &
                                p%IEC%NumTurbInp, p%IEC%IECTurbC, p%IEC%PerTurbInt, p%met%KHtest, ErrStat2, ErrMsg2)
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')
         IF (ErrStat >= AbortErrLev) THEN
            CALL Cleanup()
            RETURN
         END IF
!*** TO DO: fix this (move it somewhere more appropriate)
      IF (p%met%KHtest) THEN
         IF ( p%met%TurbModel_ID /= SpecModel_NWTCUP ) CALL SetErrStat( ErrID_Fatal, 'The KH test can be used with the "NWTCUP" spectral model only.', ErrStat, ErrMsg, 'ReadInputFile')

         IF ( .NOT. p%WrFile(FileExt_CTS) ) THEN
            CALL WRScr( ' Coherent turbulence time step files must be generated when using the "KHTEST" option.' )
            p%WrFile(FileExt_CTS)  = .TRUE.
         ENDIF
      END IF
      ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< end processing of IECturbc input variable
   
   ! ------------ Read in the IEC wind turbulence type ---------------------------------------------
CALL ReadVar( UI, InFile, Line, "IEC_WindType", "IEC turbulence type",ErrStat2, ErrMsg2, UnEc)
   CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')

   IF ( p%met%IsIECModel .AND. p%met%TurbModel_ID /= SpecModel_MODVKM  ) THEN

      CALL Conv2UC( Line )

      p%IEC%IECTurbE   = Line(1:1)

      ! Let's see if the first character is a number (for the ETM case)
      SELECT CASE ( p%IEC%IECTurbE )
         CASE ('1')
            p%IEC%Vref = 50.0
            Line = Line(2:)
         CASE ('2')
            p%IEC%Vref = 42.5
            Line = Line(2:)
         CASE ('3')
            p%IEC%Vref = 37.5
            Line = Line(2:)
         CASE DEFAULT
               ! There's no number at the start of the string so let's move on (it's NTM).
            p%IEC%Vref = -999.9
            p%IEC%IECTurbE = ' '
      END SELECT

      SELECT CASE ( TRIM( Line ) )
         CASE ( 'NTM'   )
            p%IEC%IEC_WindType = IEC_NTM
            p%IEC%IEC_WindDesc = 'Normal Turbulence Model'
         CASE ( 'ETM'   )
            p%IEC%IEC_WindType = IEC_ETM
            p%IEC%IEC_WindDesc = 'Extreme Turbulence Model'
         CASE ( 'EWM1'  )
            p%IEC%IEC_WindType = IEC_EWM1
            p%IEC%IEC_WindDesc = 'Extreme 1-Year Wind Speed Model'
         CASE ( 'EWM50' )
            p%IEC%IEC_WindType = IEC_EWM50
            p%IEC%IEC_WindDesc = 'Extreme 50-Year Wind Speed Model'
         !CASE ( 'EWM100' )
         !   p%IEC%IEC_WindType = IEC_EWM100
         !   p%IEC%IEC_WindDesc = 'Extreme 100-Year Wind Speed Model'
         CASE DEFAULT
            CALL TS_Abort ( ' Valid entries for the IEC wind turbulence are "NTM", "xETM", "xEWM1", or "xEWM50", '// &
                             'where x is the wind turbine class (1, 2, or 3).' )
      END SELECT

      IF ( p%IEC%IEC_WindType /= IEC_NTM ) THEN

         IF (p%IEC%IECedition /= 3 .OR. p%IEC%IECstandard == 2) THEN
            CALL TS_Abort ( ' The extreme turbulence and extreme wind speed models are available with '// &
                         'the IEC 61400-1 Ed. 3 or 61400-3 scaling only.')
         ENDIF

         IF (p%IEC%Vref < 0. ) THEN
            CALL TS_Abort ( ' A wind turbine class (1, 2, or 3) must be specified with the '// &
                         'extreme turbulence and extreme wind types. (i.e. "1ETM")')
         ENDIF

         IF ( p%IEC%NumTurbInp ) THEN
            CALL TS_Abort ( ' When the turbulence intensity is entered as a percent, the IEC wind type must be "NTM".' )
         ENDIF

      ELSE

         p%IEC%IECTurbE = ' '

      ENDIF

   ELSE
      p%IEC%IEC_WindType = IEC_NTM
   ENDIF

   ! ------------ Read in the ETM c parameter (IEC 61400-1, Ed 3: Section 6.3.2.3, Eq. 19) ----------------------
UseDefault = .TRUE.
p%IEC%ETMc = 2.0_ReKi;
CALL ReadRVarDefault( UI, InFile, p%IEC%ETMc, "ETMc", 'IEC Extreme Turbulence Model (ETM) "c" parameter [m/s]', UnEc, &
                      UseDefault, ErrStat2, ErrMsg2, IGNORE=(p%IEC%IEC_WindType /= IEC_ETM ))
   CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')
   
   IF ( p%IEC%ETMc <= 0. ) THEN
      CALL SetErrStat( ErrID_Fatal, 'The ETM "c" parameter must be a positive number', ErrStat, ErrMsg, 'ReadInputFile')
   ENDIF

   
   ! ------------ Read in the wind profile type -----------------------------------------------------------------

UseDefault = .TRUE.         ! Calculate the default value
SELECT CASE ( p%met%TurbModel_ID )
   CASE ( SpecModel_TimeSer )
      IF ( p%usr%NPoints > 1 ) THEN
         p%met%WindProfileType = 'TS'
      ELSE
         p%met%WindProfileType = 'PL'
PRINT *, 'check default wind profile for time-series spectral model.'         
      END IF
   CASE ( SpecModel_GP_LLJ )
      p%met%WindProfileType = 'JET'
   CASE ( SpecModel_IECKAI,SpecModel_IECVKM,SpecModel_MODVKM )
      p%met%WindProfileType = 'IEC'
   CASE ( SpecModel_TIDAL )
      p%met%WindProfileType = 'H2L'
   CASE ( SpecModel_USRVKM )
      p%met%WindProfileType = 'USR'
   CASE ( SpecModel_API )
      p%met%WindProfileType = 'API'  ! ADDED BY YG
   CASE DEFAULT
      p%met%WindProfileType = 'IEC'
END SELECT

CALL ReadCVarDefault( UI, InFile, p%met%WindProfileType, "WindProfileType", "Wind profile type", UnEc, UseDefault, ErrStat2, ErrMsg2) !converts WindProfileType to upper case
   CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')

      ! Make sure the variable is valid for this turbulence model

   SELECT CASE ( TRIM(p%met%WindProfileType) )
      CASE ( 'JET' )
         IF ( p%met%TurbModel_ID /= SpecModel_GP_LLJ ) THEN
            CALL SetErrStat( ErrID_Fatal, 'The jet wind profile is available with the GP_LLJ spectral model only.', ErrStat, ErrMsg, 'ReadInputFile')
         ENDIF
      CASE ( 'LOG')
         IF (p%IEC%IEC_WindType /= IEC_NTM ) THEN
            CALL SetErrStat( ErrID_Fatal, 'The IEC turbulence type must be NTM for the logarithmic wind profile.', ErrStat, ErrMsg, 'ReadInputFile')
!bjj check that IEC_WindType == IEC_NTM for non-IEC
         ENDIF
      CASE ( 'PL'  )
      CASE ( 'H2L' )
         IF ( p%met%TurbModel_ID /= SpecModel_TIDAL ) THEN
            CALL SetErrStat( ErrID_Fatal, 'The "H2L" mean profile type should be used only with the "TIDAL" spectral model.', ErrStat, ErrMsg, 'ReadInputFile')
         ENDIF
      CASE ( 'IEC' )
      CASE ( 'USR' )
      CASE ( 'TS' )
         IF ( p%met%TurbModel_ID /= SpecModel_TimeSer ) THEN
            CALL SetErrStat( ErrID_Fatal, 'The "TS" mean profile type is valid only with the "TIMESR" spectral model.', ErrStat, ErrMsg, 'ReadInputFile')
         ENDIF 
      CASE ( 'API' )   ! ADDED BY Y.GUO
      CASE DEFAULT
         CALL SetErrStat( ErrID_Fatal, 'The wind profile type must be "JET", "LOG", "PL", "IEC", "USR", "H2L", or default.' , ErrStat, ErrMsg, 'ReadInputFile')
   END SELECT

   IF ( p%met%TurbModel_ID == SpecModel_TIDAL .AND. TRIM(p%met%WindProfileType) /= "H2L" ) THEN
      p%met%WindProfileType = 'H2L'
      CALL SetErrStat( ErrID_Warn, 'Overwriting wind profile type to "H2L" for the "TIDAL" spectral model.', ErrStat, ErrMsg, 'ReadInputFile')
   ENDIF

   IF ( p%met%KHtest ) THEN
      IF ( TRIM(p%met%WindProfileType) /= 'IEC' .AND. TRIM(p%met%WindProfileType) /= 'PL' ) THEN
         p%met%WindProfileType = 'IEC'
         CALL SetErrStat( ErrID_Warn, 'Overwriting wind profile type for the KH test.', ErrStat, ErrMsg, 'ReadInputFile')         
      ENDIF
   ENDIF

   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF         
   
   ! ------------ Read in the height for the reference wind speed. ---------------------------------------------

CALL ReadVar( UI, InFile, p%met%RefHt, "RefHt", "Reference height [m]",ErrStat2, ErrMsg2, UnEc)
   CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')

   IF ( p%met%RefHt <=  0.0 .AND. p%met%WindProfileType(1:1) /= 'U' )  THEN
      CALL TS_Abort ( 'The reference height must be greater than zero.' )
   ENDIF


   ! ------------ Read in the reference wind speed. -----------------------------------------------------

UseDefault = .FALSE.
p%met%URef       = -999.9
IsUnusedParameter = p%IEC%IEC_WindType > IEC_ETM  .OR. p%met%WindProfileType(1:1) == 'U' ! p%IEC%IEC_WindType > IEC_ETM == EWM models

! If we specify a Ustar (i.e. if Ustar /= "default") then we can enter "default" here,
! otherwise, we get circular logic...

   CALL ReadRVarDefault( UI, InFile, p%met%URef, "URef", "Reference wind speed [m/s]", UnEc, UseDefault, ErrStat2, ErrMsg2, &
                     IGNORE=IsUnusedParameter ) ! p%IEC%IEC_WindType > IEC_ETM == EWM models
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')

   p%met%NumUSRz = 0  ! initialize the number of points in a user-defined wind profile

   IF ( ( p%met%WindProfileType(1:1) /= 'J' .OR. .NOT. UseDefault) .AND. .NOT. IsUnUsedParameter ) THEN
      IF ( p%met%URef <=  0.0 )  THEN
         CALL TS_Abort ( 'The reference wind speed must be greater than zero.' )
      ENDIF

   ELSEIF ( p%met%WindProfileType(1:1) == 'U' ) THEN ! for user-defined wind profiles, we overwrite RefHt and URef because they don't mean anything otherwise
      p%met%RefHt = p%grid%HubHt
      CALL GetUSR( UI, InFile, 39, p%met, ErrStat, ErrMsg ) !Read the last several lines of the file, then return to line 39
      IF (ErrStat >= AbortErrLev) RETURN
      p%met%URef = getVelocity(p, p%met%URef, p%met%RefHt, p%met%RefHt) !This is UHub
   ENDIF   ! Otherwise, we're using a Jet profile with default wind speed (for now it's -999.9)

   
   ! ------------ Read in the jet height -------------------------------------------------------------

UseDefault = .FALSE.
p%met%ZJetMax    = -999.9
IsUnUsedParameter = p%met%WindProfileType(1:1) /= 'J'
CALL ReadRVarDefault( UI, InFile, p%met%ZJetMax, "ZJetMax", "Jet height [m]", UnEc, UseDefault, ErrStat2, ErrMsg2, IGNORE=IsUnusedParameter)
   CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')

   IF ( .NOT. IsUnusedParameter .AND. .NOT. UseDefault ) THEN
      IF ( p%met%ZJetMax <  70.0 .OR. p%met%ZJetMax > 490.0 )  THEN
         CALL TS_Abort ( 'The height of the maximum jet wind speed must be between 70 and 490 m.' )
      ENDIF
   ENDIF


   ! ------------ Read in the power law exponent, PLExp ---------------------------------------------

SELECT CASE ( p%met%TurbModel_ID )
   CASE (SpecModel_WF_UPW, SpecModel_WF_07D, SpecModel_WF_14D, SpecModel_NWTCUP)
      IF ( p%met%KHtest ) THEN
         UseDefault = .TRUE.
         p%met%PLExp      = 0.3
      ELSE
         UseDefault = .FALSE.             ! This case needs the Richardson number to get a default
         p%met%PLExp      = 0.
      ENDIF

   CASE DEFAULT
      UseDefault = .TRUE.
      p%met%PLExp      = PowerLawExp( p%met%Rich_No )  ! These cases do not use the Richardson number to get a default

END SELECT
getPLExp = .NOT. UseDefault
IsUnusedParameter = INDEX( 'JLUHA', p%met%WindProfileType(1:1) ) > 0 .OR. p%IEC%IEC_WindType > IEC_ETM  !i.e. PL or IEC

CALL ReadRVarDefault( UI, InFile, p%met%PLExp, "PLExp", "Power law exponent", UnEc, UseDefault, ErrStat2, ErrMsg2, IGNORE=IsUnusedParameter)
   CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')

   IF ( .NOT. IsUnusedParameter .AND. .NOT. UseDefault ) THEN  ! We didn't use a default (we entered a number on the line)
      getPLExp = .FALSE.

      IF ( p%met%KHtest ) THEN
         IF ( p%met%PLExp /= 0.3 ) THEN
            p%met%PLExp = 0.3
            CALL TS_Warn  ( 'Overwriting the power law exponent for KH test.', -1 )
         ENDIF
      ENDIF
   ENDIF


   ! ------------ Read in the surface roughness length, Z0 (that's z-zero) ---------------------------------------------

UseDefault = .TRUE.
SELECT CASE ( p%met%TurbModel_ID )
   CASE (SpecModel_SMOOTH)
      p%met%Z0 = 0.010
   CASE (SpecModel_GP_LLJ )
      p%met%Z0 = 0.005
   CASE (SpecModel_WF_UPW )
      p%met%Z0 = 0.018
   CASE (SpecModel_NWTCUP )
      p%met%Z0 = 0.021
   CASE (SpecModel_WF_07D )
      p%met%Z0 = 0.233
   CASE (SpecModel_WF_14D )
      p%met%Z0 = 0.064
   CASE DEFAULT !IEC values
      p%met%Z0 = 0.030 ! Represents smooth, homogenous terrain
END SELECT
IsUnusedParameter = p%met%TurbModel_ID==SpecModel_TIDAL

CALL ReadRVarDefault( UI, InFile, p%met%Z0, "Z0", "Surface roughness length [m]", UnEc, UseDefault, ErrStat2, ErrMsg2, &
                       IGNORE=IsUnusedParameter)
   CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')

   IF ( p%met%Z0 <= 0.0 ) THEN
      CALL TS_Abort ( 'The surface roughness length must be a positive number or "default".')
   ENDIF

      
   
   !===============================================================================================================================
   ! Read the meteorological boundary conditions for non-IEC models. !
   !===============================================================================================================================

IF ( .NOT. p%met%IsIECModel .OR. p%met%TurbModel_ID == SpecModel_MODVKM  ) THEN  

   CALL ReadCom( UI, InFile, "Non-IEC Meteorological Boundary Conditions Heading Line 1", ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')

   CALL ReadCom( UI, InFile, "Non-IEC Meteorological Boundary Conditions Heading Line 2", ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')



      ! ------------ Read in the site latitude, LATITUDE. ---------------------------------------------

   UseDefault = .TRUE.
   p%met%Latitude   = 45.0

   CALL ReadRVarDefault( UI, InFile, p%met%Latitude, "Latitude", "Site latitude [degrees]", UnEc, UseDefault, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')

      IF ( ABS(p%met%Latitude) < 5.0 .OR. ABS(p%met%Latitude) > 90.0 ) THEN
         CALL TS_Abort( 'The latitude must be between -90 and 90 degrees but not between -5 and 5 degrees.' )
      ENDIF

   p%met%Fc = 2.0 * Omega * SIN( ABS(p%met%Latitude*D2R) )  ! Calculate Coriolis parameter from latitude

ELSE

   p%met%Latitude = 0.0                               !Not used in IEC specs
   p%met%Fc = 0.0

ENDIF    ! Not IECKAI and Not IECVKM




IF ( .NOT. p%met%IsIECModel  ) THEN


      ! ------------ Read in the gradient Richardson number, RICH_NO. ---------------------------------------------

   CALL ReadVar( UI, InFile, p%met%Rich_No, "RICH_NO", "Gradient Richardson number",ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')

   IF ( p%met%KHtest ) THEN
      IF ( p%met%Rich_No /= 0.02 ) THEN
         p%met%Rich_No = 0.02
         CALL TS_Warn ( 'Overwriting the Richardson Number for KH test.', -1 )
      ENDIF
   ENDIF

   IF ( p%met%TurbModel_ID == SpecModel_USER .OR. p%met%TurbModel_ID == SpecModel_USRVKM ) THEN
      IF ( p%met%Rich_No /= 0.0 ) THEN
         p%met%Rich_No = 0.0
         CALL TS_Warn ( 'Overwriting the Richardson Number for the '//TRIM(p%met%TurbModel)//' model.', -1 )
      ENDIF
   ELSEIF ( p%met%TurbModel_ID == SpecModel_NWTCUP .or. p%met%TurbModel_ID == SpecModel_GP_LLJ) THEN
      p%met%Rich_No = MIN( MAX( p%met%Rich_No, REAL(-1.0,ReKi) ), REAL(1.0,ReKi) )  ! Ensure that: -1 <= RICH_NO <= 1
   ENDIF

      ! now that we have Rich_No, we can calculate ZL and L
   CALL Calc_MO_zL(p%met%TurbModel_ID, p%met%Rich_No, p%grid%HubHt, p%met%ZL, p%met%L )

   
      ! ***** Calculate power law exponent, if needed *****

   IF ( getPLExp ) THEN
      p%met%PLExp = PowerLawExp( p%met%Rich_No )
   ENDIF

      ! ------------ Read in the shear/friction velocity, Ustar, first calculating UstarDiab ------------------------

         ! Set up the heights for the zl- and ustar-profile averages across the rotor disk
      TmpZary  = (/ p%grid%HubHt-p%grid%RotorDiameter/2., p%grid%HubHt, p%grid%HubHt+p%grid%RotorDiameter/2. /)
      IF (TmpZary(3) .GE. profileZmin .AND. TmpZary(1) .LE. profileZmax ) THEN  !set height limits so we don't extrapolate too far
         DO TmpIndex = 1,3
            TmpZary(TmpIndex) = MAX( MIN(TmpZary(TmpIndex), profileZmax), profileZmin)
         ENDDO
      ENDIF


   UseDefault = .TRUE.

   p%met%UstarDiab  = getUstarDiab(p%met%URef, p%met%RefHt, p%met%z0, p%met%ZL)
   p%met%Ustar      = p%met%UstarDiab

   SELECT CASE ( p%met%TurbModel_ID )

      CASE (SpecModel_WF_UPW)

         IF ( p%met%ZL < 0.0 ) THEN
            p%met%Ustar = 1.162 * p%met%UstarDiab**( 2.0 / 3.0 )
         ELSE ! Include the neutral case to avoid strange discontinuities
            p%met%Ustar = 0.911 * p%met%UstarDiab**( 2.0 / 3.0 )
         ENDIF

      CASE ( SpecModel_WF_07D, SpecModel_WF_14D  )

         IF ( p%met%ZL < 0.0 ) THEN
            p%met%Ustar = 1.484 * p%met%UstarDiab**( 2.0 / 3.0 )
         ELSE ! Include the neutral case with the stable one to avoid strange discontinuities
            p%met%Ustar = 1.370 * p%met%UstarDiab**( 2.0 / 3.0 )
         ENDIF

      CASE (SpecModel_GP_LLJ )
         IF ( p%met%URef < 0 ) THEN ! (1) We can't get a wind speed because Uref was default
            UseDefault = .FALSE. ! We'll calculate the default value later
         ELSE
            p%met%Ustar = 0.17454 + 0.72045*p%met%UstarDiab**1.36242
         ENDIF

      CASE ( SpecModel_NWTCUP  )
         p%met%Ustar = 0.2716 + 0.7573*p%met%UstarDiab**1.2599

      CASE ( SpecModel_TIDAL , SpecModel_RIVER )
         ! Use a constant drag coefficient for the HYDRO spectral models.
         p%met%Ustar = p%met%Uref*0.05 ! This corresponds to a drag coefficient of 0.0025.
         !p%met%Ustar = p%met%Uref*0.04 ! This corresponds to a drag coefficient of 0.0016.

   END SELECT

   CALL ReadRVarDefault( UI, InFile, p%met%Ustar, "UStar", "Friction or shear velocity [m/s]", UnEc, UseDefault, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')


   IF ( p%met%Uref < 0.0 .AND. UseDefault ) THEN  ! This occurs if "default" was entered for both GP_LLJ wind speed and UStar
      CALL TS_Abort( 'The reference wind speed and friction velocity cannot both be "default."')
   ELSEIF (p%met%Ustar <= 0) THEN
      CALL TS_Abort( 'The friction velocity must be a positive number.')
   ENDIF


         ! ***** Calculate wind speed at hub height *****

   IF ( p%met%WindProfileType(1:1) == 'J' ) THEN
      IF ( p%met%ZJetMax < 0 ) THEN ! Calculate a default value
         p%met%ZJetMax = -14.820561*p%met%Rich_No + 56.488123*p%met%ZL + 166.499069*p%met%Ustar + 188.253377
         p%met%ZJetMax =   1.9326  *p%met%ZJetMax - 252.7267  ! Correct with the residual

         CALL RndJetHeight( p%RNG, OtherSt_RandNum, tmp ) ! Add a random amount

         p%met%ZJetMax = MIN( MAX(p%met%ZJetMax + tmp, 70.0_ReKi ), 490.0_ReKi )
      ENDIF

      IF ( p%met%URef < 0 ) THEN ! Calculate a default value

         p%met%UJetMax = MAX( -21.5515_ReKi + 6.6827_ReKi*LOG(p%met%ZJetMax), 5.0_ReKi ) !Jet max must be at least 5 m/s (occurs ~50 m); shouldn't happen, but just in case....

         CALL Rnd3ParmNorm( p%RNG, OtherSt_RandNum, tmp, 0.1076_ReKi, -0.1404_ReKi, 3.6111_ReKi,  -15.0_ReKi, 20.0_ReKi )

         IF (p%met%UJetMax + tmp > 0 ) p%met%UJetMax = p%met%UJetMax + tmp

         CALL GetChebCoefs( p%met%UJetMax, p%met%ZJetMax ) ! These coefficients are a function of UJetMax, ZJetMax, RICH_NO, and p%met%Ustar

         p%met%URef = getVelocity(p, p%met%UJetMax, p%met%ZJetMax, p%met%RefHt)

      ELSE
         CALL GetChebCoefs(p%met%URef, p%met%RefHt)
      ENDIF

   ENDIF !Jet wind profile

   p%UHub = getVelocity(p, p%met%URef, p%met%RefHt, p%grid%HubHt)

         ! ***** Get p%met%Ustar- and zl-profile values, if required, and determine offsets *****
      IF ( p%met%TurbModel_ID /= SpecModel_GP_LLJ ) THEN
         p%met%UstarSlope = 1.0_ReKi         

         TmpUary   = (/ 0.0_ReKi, 0.0_ReKi, 0.0_ReKi /)
         TmpUstar  = (/ 0.0_ReKi, 0.0_ReKi, 0.0_ReKi /)
         
         p%met%UstarOffset= 0.0_ReKi
      ELSE
         p%met%UstarSlope = 1.0_ReKi         
         p%met%UstarDiab = getUstarDiab(p%met%URef, p%met%RefHt, p%met%z0, p%met%ZL) !bjj: is this problematic for anything else?

         TmpUary   = getVelocityProfile(p, p%met%URef, p%met%RefHt, TmpZary)                  
         TmpUstar  = getUstarARY( TmpUary,     TmpZary, 0.0_ReKi, p%met%UstarSlope )
            
         p%met%UstarOffset = p%met%Ustar - SUM(TmpUstar) / SIZE(TmpUstar)    ! Ustar minus the average of those 3 points
         TmpUstar(:) = TmpUstar(:) + p%met%UstarOffset
      ENDIF

      TmpZLary = getZLARY(TmpUary, TmpZary, p%met%Rich_No, p%met%ZL, p%met%L, 0.0_ReKi, p%met%WindProfileType)
      p%met%zlOffset = p%met%ZL - SUM(TmpZLary) / SIZE(TmpZLary)


      ! ------------- Read in the mixing layer depth, ZI ---------------------------------------------

   UseDefault = .TRUE.
   IF ( p%met%ZL >= 0.0 .AND. p%met%TurbModel_ID /= SpecModel_GP_LLJ ) THEN  !We must calculate ZI for stable GP_LLJ. z/L profile can change signs so ZI must be defined for spectra.
      p%met%ZI = 0.0
   ELSE
      IF ( p%met%Ustar < p%met%UstarDiab ) THEN
         p%met%ZI = ( 0.04 * p%met%Uref ) / ( 1.0E-4 * LOG10( p%met%RefHt / p%met%Z0 ) )  !for "very" windy days
      ELSE
         !Should Wind Farm models use the other definition since that was what was used in creating those models?
         p%met%ZI = p%met%Ustar / (6.0 * p%met%Fc)
      ENDIF
   ENDIF

   CALL ReadRVarDefault( UI, InFile, p%met%ZI, "ZI", "Mixing layer depth [m]", UnEc, UseDefault, ErrStat2, ErrMsg2, &
                                                                  IGNORE=(p%met%ZL>=0. .AND. p%met%TurbModel_ID /= SpecModel_GP_LLJ) )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')

      IF ( ( p%met%ZL < 0.0 ) .AND. ( p%met%ZI <= 0.0 ) ) THEN
         CALL TS_Abort ( 'The mixing layer depth must be a positive number for unstable flows.')
      ENDIF


      ! Get the default mean Reynolds stresses

   CALL GetDefaultRS(  p%met%PC_UW, p%met%PC_UV, p%met%PC_VW, p%met%UWskip, p%met%UVskip, p%met%VWskip, TmpUStar(2) )


       ! ----------- Read in the mean hub u'w' Reynolds stress, PC_UW ---------------------------------------------
   UseDefault = .TRUE.
   CALL ReadRVarDefault( UI, InFile, p%met%PC_UW, "PC_UW", &
                                            "Mean hub u'w' Reynolds stress", UnEc, UseDefault, ErrStat2, ErrMsg2, IGNORESTR = p%met%UWskip )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')

   IF (.NOT. p%met%UWskip) THEN
      TmpUstarD = ( TmpUstar(1)- 2.0*TmpUstar(2) + TmpUstar(3) )

      IF ( TmpUstarD /= 0.0 ) THEN
         p%met%UstarSlope  = 3.0*(p%met%Ustar -  SQRT( ABS(p%met%PC_UW) ) ) / TmpUstarD
         p%met%UstarOffset = SQRT( ABS(p%met%PC_UW) ) - p%met%UstarSlope*(TmpUstar(2) - p%met%UstarOffset)
      ELSE
         p%met%UstarSlope  = 0.0
         p%met%UstarOffset = SQRT( ABS(p%met%PC_UW) )
      ENDIF
   ENDIF

      ! ------------ Read in the mean hub u'v' Reynolds stress, PC_UV ---------------------------------------------
   UseDefault = .TRUE.
   CALL ReadRVarDefault( UI, InFile, p%met%PC_UV, "PC_UV", &
                                            "Mean hub u'v' Reynolds stress", UnEc, UseDefault, ErrStat2, ErrMsg2, IGNORESTR = p%met%UVskip )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')

      ! ------------ Read in the mean hub v'w' Reynolds stress, PC_VW ---------------------------------------------
   UseDefault = .TRUE.
   CALL ReadRVarDefault( UI, InFile, p%met%PC_VW, "PC_VW", &
                                            "Mean hub v'w' Reynolds stress", UnEc, UseDefault, ErrStat2, ErrMsg2, IGNORESTR = p%met%VWskip )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')

      ! ------------ Read in the u component coherence decrement (coh-squared def), InCDec(1) = InCDecU ------------
   CALL GetDefaultCoh( p%met%TurbModel_ID, p%met%RICH_NO, p%UHub, p%grid%HubHt, p%met%IncDec, p%met%InCohB )
   UseDefault = .TRUE.
   InCVar(1) = p%met%InCDec(1)
   InCVar(2) = p%met%InCohB(1)

   CALL ReadRAryDefault( UI, InFile, InCVar, "InCDec1", &
                                             "u-component coherence parameters", UnEc, UseDefault, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')
   
   p%met%InCDec(1) = InCVar(1)
   p%met%InCohB(1) = InCVar(2)

      IF ( p%met%InCDec(1) <= 0.0 ) THEN
         CALL SetErrStat( ErrID_Fatal, 'The u-component coherence decrement must be a positive number.', ErrStat, ErrMsg, 'ReadInputFile')
      ENDIF

      ! ------------ Read in the v component coherence decrement (coh-squared def), InCDec(2) = InCDecV ----------

   UseDefault = .TRUE.
   InCVar(1) = p%met%InCDec(2)
   InCVar(2) = p%met%InCohB(2)

   CALL ReadRAryDefault( UI, InFile, InCVar, "InCDec2",  &
                                             "v-component coherence parameters", UnEc, UseDefault, ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')

   p%met%InCDec(2) = InCVar(1)
   p%met%InCohB(2) = InCVar(2)

      IF ( p%met%InCDec(2) <= 0.0 ) THEN
         CALL SetErrStat( ErrID_Fatal, 'The v-component coherence decrement must be a positive number.', ErrStat, ErrMsg, 'ReadInputFile')
      ENDIF

      ! ------------ Read in the w component coherence decrement (coh-squared def), InCDec(3) = InCDecW -------

   UseDefault = .TRUE.
   InCVar(1) = p%met%InCDec(3)
   InCVar(2) = p%met%InCohB(3)

   CALL ReadRAryDefault( UI, InFile, InCVar, "InCDec3", &
                                             "w-component coherence parameters", UnEc, UseDefault, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')

   p%met%InCDec(3) = InCVar(1)
   p%met%InCohB(3) = InCVar(2)

      IF ( p%met%InCDec(3) <= 0.0 ) THEN
         CALL SetErrStat( ErrID_Fatal, 'The w-component coherence decrement must be a positive number.', ErrStat, ErrMsg, 'ReadInputFile')
      ENDIF

         ! ------------ Read in the coherence exponent, COHEXP -----------------------------------

   UseDefault = .TRUE.
   p%met%CohExp     = 0.0    ! was 0.25
   CALL ReadRVarDefault( UI, InFile, p%met%CohExp, "CohExp", "Coherence exponent", UnEc, UseDefault, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')

      IF ( p%met%COHEXP < 0.0 ) THEN
         CALL SetErrStat( ErrID_Fatal, 'The coherence exponent must be non-negative.', ErrStat, ErrMsg, 'ReadInputFile')
      ENDIF


   !===============================================================================================================================
   ! Read the Coherent Turbulence Scaling Parameters, if necessary.  
   !===============================================================================================================================

   IF ( p%WrFile(FileExt_CTS) ) THEN

      CALL ReadCom( UI, InFile, "Coherent Turbulence Scaling Parameters Heading Line 1", ErrStat2, ErrMsg2, UnEc )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')

      CALL ReadCom( UI, InFile, "Coherent Turbulence Scaling Parameters Heading Line 2", ErrStat2, ErrMsg2, UnEc )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')


         ! ------------ Read the name of the path containg event file definitions, CTEventPath --------------------------

      CALL ReadVar( UI, InFile, p_CohStr%CTEventPath, "CTEventPath", "Coherence events path",ErrStat2, ErrMsg2, UnEc)
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')

      CALL ReadVar( UI, InFile, Line, "CTEventFile", "Event file type",ErrStat2, ErrMsg2, UnEc)
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')


      IF ( p%met%KHtest ) THEN

         p_CohStr%CText = 'les'
         p_CohStr%CTEventFile = TRIM(p_CohStr%CTEventPath)//PathSep//'Events.xtm'

         CALL WrScr( ' LES events will be used for the KH test.' )

      ELSE

         p_CohStr%CText = Line  !This will preserve the case formatting, in case it matters.

         CALL Conv2UC( Line )

         IF (Line(1:6) == "RANDOM") THEN
             CALL RndUnif( p%RNG, OtherSt_RandNum, tmp )

             IF ( tmp <= 0.5 ) THEN
                 p_CohStr%CText = 'les'
             ELSE
                 p_CohStr%CText = 'dns'
             ENDIF
         ENDIF

         p_CohStr%CTEventFile = TRIM(p_CohStr%CTEventPath)//PathSep//'Events.'//TRIM(p_CohStr%CText)

      ENDIF


         ! ------------ Read the Randomization Flag, Randomize -----------------------------------
      CALL ReadVar( UI, InFile, Randomize, "Randomize", "Randomize CT Scaling",ErrStat2, ErrMsg2, UnEc)
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')

      IF ( p%met%KHtest ) THEN
         Randomize = .FALSE.
         CALL WrScr( ' Billow will cover rotor disk for KH test. ' )
      ENDIF

         ! ------------ Read the Disturbance Scale, DistScl ---------------------------------------------
      CALL ReadVar( UI, InFile, p_CohStr%DistScl, "DistScl", "Disturbance scale",ErrStat2, ErrMsg2, UnEc)
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')

         ! ------------ Read the Lateral Fractional Location of tower centerline in wave, CTLy ----------
      CALL ReadVar( UI, InFile, p_CohStr%CTLy, "CTLy", "Location of tower centerline",ErrStat2, ErrMsg2, UnEc)
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')

         ! ------------ Read the Vertical Fraction Location of hub in wave, CTLz ------------------------
      CALL ReadVar( UI, InFile, p_CohStr%CTLz, "CTLz", "Location of hub height",ErrStat2, ErrMsg2, UnEc)
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')

      IF ( p%met%KHtest ) THEN
            p_CohStr%DistScl = 1.0
            p_CohStr%CTLy    = 0.5
            p_CohStr%CTLz    = 0.5
      ELSEIF ( Randomize ) THEN

         CALL RndUnif( p%RNG, OtherSt_RandNum, tmp )

            ! Assume a 75% chance of coherent turbulence being the size of the rotor
            ! If the rotor is too small, assume a 100% chance.
            ! If the turbulence is not the size of the rotor, assume it's half the size
            ! of the disk, with equal probablilty of being in the lower or upper half.

         IF ( tmp > 0.25 .OR. p%grid%RotorDiameter <= 30.0 ) THEN

            p_CohStr%DistScl = 1.0
            p_CohStr%CTLy    = 0.5
            p_CohStr%CTLz    = 0.5

         ELSE

            p_CohStr%DistScl = 0.5
            p_CohStr%CTLy    = 0.5

            IF ( tmp < 0.125 ) THEN
               p_CohStr%CTLz = 0.0 ! The hub is on the bottom of the dataset (i.e. the billow is on the top of the disk)
            ELSE
               p_CohStr%CTLz = 1.0 ! The hub is on the top of the dataset
            ENDIF

         ENDIF

      ELSE  !Don't randomize:

         IF ( p_CohStr%DistScl < 0.0 ) THEN
            CALL TS_Abort ('The disturbance scale must be a positive.')
         ELSEIF ( p%grid%RotorDiameter <= 30.0 .AND. p_CohStr%DistScl < 1.0 ) THEN
            CALL TS_Abort ('The disturbance scale must be at least 1.0 for rotor diameters less than 30.')
         ELSEIF ( p%grid%RotorDiameter*p_CohStr%DistScl <= 15.0  ) THEN
            CALL TS_Abort ('The coherent turbulence must be greater than 15 meters in height.  '//&
                        'Increase the rotor diameter or the disturbance scale. ')
         ENDIF

      ENDIF


         ! ---------- Read the Minimum event start time, CTStartTime --------------------------------------------

      CALL ReadVar( UI, InFile, p_CohStr%CTStartTime, "CTStartTime", "CTS Start Time",ErrStat2, ErrMsg2, UnEc)
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadInputFile')

      p_CohStr%CTStartTime = MAX( p_CohStr%CTStartTime, REAL(0.0,ReKi) ) ! A Negative start time doesn't really make sense...

   ENDIF    ! WrFile(FileExt_CTS)


ELSE  ! IECVKM, IECKAI, MODVKM, OR API models

   p%met%Rich_No = 0.0                       ! Richardson Number in neutral conditions
   p%met%COHEXP  = 0.0                       ! Coherence exponent
   
      ! The following variables are not used in the IEC calculations

   p%met%ZL      = 0.0                       ! M-O z/L parameter
   p%met%L       = 0.0                       ! M-O length scale
   p%met%Ustar   = 0.0                       ! Shear or friction velocity
   p%met%ZI      = 0.0                       ! Mixing layer depth
   p%met%PC_UW   = 0.0                       ! u'w' x-axis correlation coefficient
   p%met%PC_UV   = 0.0                       ! u'v' x-axis correlation coefficient
   p%met%PC_VW   = 0.0                       ! v'w' x-axis correlation coefficient

   p%met%UWskip  = .TRUE.
   p%met%UVskip  = .TRUE.
   p%met%VWskip  = .TRUE.
   
   
   IF ( p%IEC%NumTurbInp .AND. p%IEC%PerTurbInt == 0.0 ) THEN    ! This will produce constant winds, instead of an error when the transfer matrix is singular
      p%met%TurbModel = 'NONE'
      p%met%TurbModel_ID = SpecModel_NONE
   ENDIF

      ! Calculate wind speed at hub height

   p%UHub    = getVelocity(p, p%met%URef, p%met%RefHt, p%grid%HubHt)


ENDIF




   ! Done reading the input file.

CALL Cleanup()

RETURN
CONTAINS
!.........................................
SUBROUTINE Cleanup()
   
      IF ( UI   > 0 ) CLOSE( UI)
      IF ( UnEc > 0 ) CLOSE( UnEc )
   
   END SUBROUTINE Cleanup   
!.........................................
END SUBROUTINE ReadInputFile
!=======================================================================
SUBROUTINE OpenSummaryFile(RootName, DescStr)

  ! This subroutine is used to open the summary output file.

USE              TSMods

IMPLICIT         NONE

CHARACTER(*), INTENT(IN)  :: RootName  ! rootname of the primary TurbSim input file
CHARACTER(*), INTENT(OUT) :: DescStr   ! string describing time TurbSim files were generated




   ! Open summary file.

CALL OpenFOutFile( US, TRIM( RootName )//'.sum' ) ! Formatted output file


   ! Write the program name and version, date and time into the summary file.

   ! Let's make sure the binary file and the full-field file have the same date and time.
DescStr = 'generated by '//TRIM( GetNVD(TurbSim_Ver) )//' on '//CurDate()//' at '//CurTime()//'.'

WRITE (US,"( / 'This summary file was ', A / )")  TRIM(DescStr)

   ! Capitalize the first letter of the string.

DescStr = 'This full-field file was '//TRIM(DescStr)


RETURN
END SUBROUTINE OpenSummaryFile
!=======================================================================
SUBROUTINE GetUSR(U_in, FileName, NLines, p_met, ErrStat, ErrMsg)

   USE                                   TSMods, ONLY: p
   IMPLICIT                              NONE

   TYPE(Meteorology_ParameterType), intent(inout) :: p_met
   INTEGER(IntKi),                  intent(  out) :: ErrStat                         ! Error level
   CHARACTER(*),                    intent(  out) :: ErrMsg                          ! Message describing error
   CHARACTER(*),                    INTENT(IN   ) :: FileName                        ! Name of the input file
   
   ! local variables
   
   INTEGER(IntKi)                                 :: ErrStat2                        ! Error level (local)
   CHARACTER(MaxMsgLen)                           :: ErrMsg2                         ! Message describing error (local)
   
   
   CHARACTER(200)                     :: LINE

   REAL(ReKi)                         :: L_Usr_Tmp
   REAL(ReKi)                         :: Sigma_USR_Tmp
   REAL(ReKi)                         :: U_USR_Tmp
   REAL(ReKi)                         :: WindDir_USR_Tmp
   REAL(ReKi)                         :: Z_USR_Tmp

   INTEGER                            :: I
   INTEGER                            :: Indx
   INTEGER                            :: J
   INTEGER                            :: IOAstat                        ! Input/Output/Allocate status
   INTEGER, INTENT(IN), OPTIONAL      :: NLines                         ! Number of lines to be skipped, if the file must be rewound
   INTEGER, INTENT(IN)                :: U_in                           ! Input unit.  This file is assumed to be already open

   LOGICAL                            :: ReadSigL                       ! Whether or not to read the last 2 columns

   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
      ! Find the end of the input file, where the "User-Defined Variables" are located
   READ ( U_in, '(A)', IOSTAT=IOAstat ) LINE

   IF ( IOAstat /= 0 ) THEN
      CALL TS_Abort( 'Could not read entire input file for user-defined variables.' )
   ENDIF

   CALL Conv2UC ( LINE )

   DO WHILE ( INDEX( LINE, 'USER-DEFINED' ) == 0 )

      READ ( U_in, '(A)', IOSTAT=IOAstat ) LINE

      IF ( IOAstat /= 0 ) THEN
         CALL TS_Abort( 'Could not read entire input file for user-defined variables.' )
      ENDIF

      CALL Conv2UC( LINE )

   ENDDO


      ! ---------- Read the size of the arrays --------------------------------------------
   CALL ReadVar( U_in, FileName, p_met%NumUSRz, "NumUSRz", "Number of heights in the user-defined profiles" )

   IF ( p_met%NumUSRz < 1 ) THEN
      CALL TS_Abort( 'The number of heights specified in the user-defined profiles must be at least 1.')
   ENDIF

   DO I=1,3
         ! ---------- Read the size of the arrays --------------------------------------------
      CALL ReadVar( U_in, FileName, p_met%USR_StdScale(I), "USR_StdScale", "Scaling value for user-defined standard deviation profile" )

      IF ( p_met%USR_StdScale(I) <= 0. ) THEN
         CALL TS_Abort( 'The scaling value for the user-defined standard deviation profile must be positive.')
      ENDIF
   ENDDO

      ! Allocate the data arrays
   CALL AllocAry(p_met%USR_Z,       p_met%NumUSRz, 'USR_Z (user-defined height)',               ErrStat2, ErrMsg2); CALL SetErrStat(ErrSTat2, ErrMsg2, ErrStat, ErrMsg, 'GetUSR')
   CALL AllocAry(p_met%USR_U,       p_met%NumUSRz, 'USR_U (user-defined wind speed)',           ErrStat2, ErrMsg2); CALL SetErrStat(ErrSTat2, ErrMsg2, ErrStat, ErrMsg, 'GetUSR')
   CALL AllocAry(p_met%USR_WindDir, p_met%NumUSRz, 'USR_WindDir (user-defined wind direction)', ErrStat2, ErrMsg2); CALL SetErrStat(ErrSTat2, ErrMsg2, ErrStat, ErrMsg, 'GetUSR')


   IF ( p_met%TurbModel_ID == SpecModel_USRVKM ) THEN
      ReadSigL = .TRUE.

      CALL AllocAry(p_met%USR_Sigma, p_met%NumUSRz, 'USR_Sigma (user-defined sigma)',    ErrStat2, ErrMsg2); CALL SetErrStat(ErrSTat2, ErrMsg2, ErrStat, ErrMsg, 'GetUSR')
      CALL AllocAry(p_met%USR_L,     p_met%NumUSRz, 'USR_L (user-defined length scale)', ErrStat2, ErrMsg2); CALL SetErrStat(ErrSTat2, ErrMsg2, ErrStat, ErrMsg, 'GetUSR')
      
   ELSE
      ReadSigL = .FALSE.
   ENDIF

      ! ---------- Skip 4 lines --------------------------------------------
   DO I=1,4
      CALL ReadCom( U_in, FileName, "Headers for user-defined variables" )
   ENDDO

   DO I=1,p_met%NumUSRz

      IF ( ReadSigL ) THEN
         READ( U_in, *, IOSTAT=IOAstat ) p_met%USR_Z(I), p_met%USR_U(I), p_met%USR_WindDir(I), p_met%USR_Sigma(I), p_met%USR_L(I)
      ELSE
         READ( U_in, *, IOSTAT=IOAstat ) p_met%USR_Z(I), p_met%USR_U(I), p_met%USR_WindDir(I)
      ENDIF

      IF ( IOAstat /= 0 ) THEN
         CALL TS_Abort( 'Could not read entire user-defined variable list on line '//Int2LStr(I)//'.' )
      ENDIF

      IF ( ReadSigL ) THEN
         IF ( p_met%USR_Sigma(I) <= REAL( 0., ReKi ) ) THEN
            CALL TS_Abort( 'The standard deviation must be a positive number.' );
         ELSEIF ( p_met%USR_L(I) <= REAL( 0., ReKi ) ) THEN
            CALL TS_Abort( 'The length scale must be a positive number.' );
         ENDIF
      ENDIF

      IF ( p_met%USR_WindDir(I) > 360. ) THEN
         J = INT ( p_met%USR_WindDir(I) / 360. )
         p_met%USR_WindDir(I) = p_met%USR_WindDir(I) - J * 360.
      ELSEIF ( p_met%USR_WindDir(I) < 0. ) THEN
         J = INT ( -p_met%USR_WindDir(I) / 360. ) +1
         p_met%USR_WindDir(I) = p_met%USR_WindDir(I) + J * 360.
      ENDIF
   ENDDO

      ! Sort the arrays
   DO I=2,p_met%NumUSRz
      IF ( p_met%USR_Z(I) < p_met%USR_Z(I-1) ) THEN

         Indx = 1
         DO J=I-2,1,-1
            IF ( p_met%USR_Z(I) > p_met%USR_Z(J) ) THEN
               Indx = J+1
               EXIT
            ELSEIF ( p_met%USR_Z(I) == p_met%USR_Z(J) ) THEN
               CALL TS_Abort( 'Error: user-defined values must contain unique heights.' )
            ENDIF
         ENDDO

         Z_USR_Tmp       = p_met%USR_Z(I)
         U_USR_Tmp       = p_met%USR_U(I)
         WindDir_USR_Tmp = p_met%USR_WindDir(I)

         DO J=I,Indx+1,-1
            p_met%USR_Z(J)       = p_met%USR_Z(J-1)
            p_met%USR_U(J)       = p_met%USR_U(J-1)
            p_met%USR_WindDir(J) = p_met%USR_WindDir(J-1)
         ENDDO

         p_met%USR_Z(Indx)       = Z_USR_Tmp
         p_met%USR_U(Indx)       = U_USR_Tmp
         p_met%USR_WindDir(Indx) = WindDir_USR_Tmp

         IF ( ReadSigL ) THEN
            Sigma_USR_Tmp   = p_met%USR_Sigma(I)
            L_USR_Tmp       = p_met%USR_L(I)

            DO J=I,Indx+1,-1
               p_met%USR_Sigma(J)   = p_met%USR_Sigma(J-1)
               p_met%USR_L(J)       = p_met%USR_L(J-1)
            ENDDO

            p_met%USR_Sigma(Indx)   = Sigma_USR_Tmp
            p_met%USR_L(Indx)       = L_USR_Tmp
         ENDIF ! ReadSigL

      ENDIF
   ENDDO

      ! Rewind the file, if necessary.
   IF ( PRESENT(NLines) ) THEN
      REWIND( U_in , IOSTAT=IOAstat )

      IF ( IOAstat /= 0 ) THEN
         CALL TS_Abort( 'Error rewinding the file '//TRIM(FileName)//'.' )
      ENDIF

      DO I = 1,NLines
         CALL ReadCom( U_in, FileName, "Line "//Int2LStr(I) )
      ENDDO
   ENDIF

END SUBROUTINE GetUSR
!=======================================================================
!> Read the input file for user-defined spectra.
SUBROUTINE GetUSRSpec(FileName, p, UnEc, ErrStat, ErrMsg)

   IMPLICIT NONE

   TYPE(TurbSim_ParameterType),     INTENT(INOUT) :: p                              !< Simulation parameters
   INTEGER(IntKi),                  INTENT(IN   ) :: UnEc                           !< Echo file unit number
   INTEGER(IntKi),                  INTENT(  OUT) :: ErrStat                        !< Error level
   CHARACTER(*),                    INTENT(  OUT) :: ErrMsg                         !< Message describing error
   CHARACTER(*),                    INTENT(IN)    :: FileName                       !< Name of the input file

      ! local variables
   REAL(ReKi)                         :: Freq_USR_Tmp
   REAL(ReKi)                         :: U_USR_Tmp
   REAL(ReKi)                         :: V_USR_Tmp
   REAL(ReKi)                         :: W_USR_Tmp
   REAL(ReKi)                         :: SpecScale (3)

   INTEGER                            :: I
   INTEGER, PARAMETER                 :: iPoint = 1 ! spectra are input for only one point 
   INTEGER                            :: Indx
   INTEGER                            :: J
   INTEGER                            :: USpec                          ! I/O unit for user-defined spectra

   
   INTEGER(IntKi)                                 :: ErrStat2                         ! Error level (local)
   CHARACTER(MaxMsgLen)                           :: ErrMsg2                          ! Message describing error (local)

   ErrStat = ErrID_None
   ErrMSg  = ""
   
      ! --------- Open the file ---------------

   CALL GetNewUnit( USpec, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2 , ErrStat, ErrMsg, 'GetUSR')
      
   CALL OpenFInpFile( USpec, FileName, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2 , ErrStat, ErrMsg, 'GetUSR')
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      ENDIF


   CALL WrScr1(' Reading the user-defined spectra input file "'//TRIM(FileName)//'".' )


      ! --------- Read the comment lines at the beginning of the file ---------------
   DO I=1,3
      CALL ReadCom( USpec, FileName, "user-spectra header line #"//TRIM(Num2LStr(I)), ErrStat2, ErrMsg2, UnEc)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'GetUSR')
   ENDDO


      ! ---------- Read the size of the arrays --------------------------------------------
   CALL ReadVar( USpec, FileName, p%usr%nFreq, "nFreq", "Number of frequencies in the user-defined spectra", ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'GetUSR')
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      ENDIF


   DO I=1,3
         ! ---------- Read the scaling for the arrays --------------------------------------------
      CALL ReadVar( USpec, FileName, SpecScale(I), "SpecScale", "Scaling value for user-defined standard deviation profile", ErrStat2, ErrMsg2, UnEc )
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'GetUSR')

   ENDDO
      
   IF ( p%usr%nFreq < 3      ) CALL SetErrStat(ErrID_Fatal, 'The number of frequencies specified in the user-defined spectra must be at least 3.' , ErrStat, ErrMsg, 'GetUSR')
   IF ( ANY(SpecScale <= 0.) ) CALL SetErrStat(ErrID_Fatal, 'The scaling value for the user-defined spectra must be positive.' , ErrStat, ErrMsg, 'GetUSR')
   
   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   ENDIF   
   
      ! Allocate the data arrays
   CALL AllocAry( p%usr%f,      p%usr%nFreq,    'f (user-defined frequencies)'  ,ErrStat2,ErrMsg2); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'GetUSR')
   CALL AllocAry( p%usr%S,      p%usr%nFreq,1,3,'S (user-defined spectra)'      ,ErrStat2,ErrMsg2); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'GetUSR')
   CALL AllocAry( p%usr%pointzi, iPoint        , 'pointzi (user-defined spectra',ErrStat2,ErrMsg2); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'GetUSR')   
   
   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF
   
   p%usr%pointzi = 0.0_ReKi ! we don't care what this is; it's only potentially used so we can use the same interpolation routine as the user time-series input

      ! ---------- Skip 4 lines --------------------------------------------
   DO I=1,4
      CALL ReadCom( USpec, FileName, "Headers for user-defined variables", ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'GetUSR')
   ENDDO

      ! ---------- Read the data lines --------------------------------------
   DO I=1,p%usr%nFreq

      READ( USpec, *, IOSTAT=ErrStat2 ) p%usr%f(I), p%usr%S(I,iPoint,1), p%usr%S(I,iPoint,2), p%usr%S(I,iPoint,3)

      IF ( ErrStat2 /= 0 ) THEN
         CALL SetErrStat(ErrID_Fatal, 'Could not read entire user-defined spectra on line '//Int2LStr(I)//'.' , ErrStat, ErrMsg, 'GetUSR')
         CALL Cleanup()
         RETURN
      ENDIF

      IF ( ANY( p%usr%S(I,iPoint,:) <=  0._ReKi ) ) THEN

         CALL SetErrStat(ErrID_Fatal, 'The spectra must contain positive numbers.' , ErrStat, ErrMsg, 'GetUSR')
         CALL Cleanup()
         RETURN
         
!      ELSEIF ( p%usr%f(I) <= REAL( 0., ReKi ) ) THEN
!         CALL TS_Abort( 'The frequencies must be a positive number.' );
      ENDIF

         ! Scale by the factors earlier in the input file

      p%usr%S(I,iPoint,1) = p%usr%S(I,iPoint,1)*SpecScale(1)
      p%usr%S(I,iPoint,2) = p%usr%S(I,iPoint,2)*SpecScale(2)
      p%usr%S(I,iPoint,3) = p%usr%S(I,iPoint,3)*SpecScale(3)

   ENDDO

      ! ------- Sort the arrays by frequency -----------------------------------
   DO I=2,p%usr%nFreq
      IF ( p%usr%f(I) < p%usr%f(I-1) ) THEN

         Indx = 1
         DO J=I-2,1,-1
            IF ( p%usr%f(I) > p%usr%f(J) ) THEN
               Indx = J+1
               EXIT
            ELSEIF ( EqualRealNos( p%usr%f(I), p%usr%f(J) ) ) THEN
               CALL SetErrStat(ErrID_Fatal, 'Error: user-defined spectra must contain unique frequencies.' , ErrStat, ErrMsg, 'GetUSR')
               CALL Cleanup()
               RETURN
            ENDIF
         ENDDO

         Freq_USR_Tmp    = p%usr%f(I)
         U_USR_Tmp       = p%usr%S(I,iPoint,1)
         V_USR_Tmp       = p%usr%S(I,iPoint,2)
         W_USR_Tmp       = p%usr%S(I,iPoint,3)

         DO J=I,Indx+1,-1
            p%usr%f(J)          = p%usr%f(J-1)
            p%usr%S(J,iPoint,1) = p%usr%S(J-1,iPoint,1)
            p%usr%S(J,iPoint,2) = p%usr%S(J-1,iPoint,2)
            p%usr%S(J,iPoint,3) = p%usr%S(J-1,iPoint,3)
         ENDDO

         p%usr%f(Indx)       = Freq_USR_Tmp
         p%usr%S(I,iPoint,1) = U_USR_Tmp
         p%usr%S(I,iPoint,2) = V_USR_Tmp
         p%usr%S(I,iPoint,3) = W_USR_Tmp

      ENDIF
   ENDDO

      ! --------- Close the file ---------------------------------------

   CALL Cleanup()
   RETURN
   
CONTAINS
   SUBROUTINE Cleanup()
   
      CLOSE( USpec )
   
   
   END SUBROUTINE Cleanup
   
END SUBROUTINE GetUSRSpec

!=======================================================================
!> Read the input file for user-defined time series
SUBROUTINE ReadUSRTimeSeries(FileName, p, UnEc, ErrStat, ErrMsg)

   IMPLICIT NONE

   TYPE(TurbSim_ParameterType),     INTENT(INOUT) :: p                              !< Simulation parameters
   INTEGER(IntKi),                  INTENT(IN   ) :: UnEc                           !< Echo file unit number
   INTEGER(IntKi),                  INTENT(  OUT) :: ErrStat                        !< Error level
   CHARACTER(*),                    INTENT(  OUT) :: ErrMsg                         !< Message describing error
   CHARACTER(*),                    INTENT(IN)    :: FileName                       !< Name of the input file

      ! local variables
   INTEGER(IntKi), PARAMETER                      :: NumLinesBeforeTS = 10          ! Number of lines in the input file before the time series start (need to add nPoint lines). IMPORTANT:  any changes to the number of lines in the header must be reflected here
      
   INTEGER(IntKi)                                 :: UnIn                           ! unit number for reading input file
   INTEGER(IntKi)                                 :: I, J                           ! loop counters
   INTEGER(IntKi)                                 :: IPoint                         ! loop counter on number of points
   INTEGER(IntKi)                                 :: IVec                           ! loop counter on velocity components being read
   INTEGER(IntKi)                                 :: ErrStat2                       ! Error level (local)
   CHARACTER(MaxMsgLen)                           :: ErrMsg2                        ! Message describing error (local)
   
   CHARACTER(200)                                 :: FormStr          
   CHARACTER(1)                                   :: tmpChar          
   real(reKi)                                     :: tmpAry(3)
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
      ! --------- Open the file ---------------

   CALL GetNewUnit( UnIn, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2 , ErrStat, ErrMsg, 'ReadUSRTimeSeries')
      
   CALL OpenFInpFile( UnIn, FileName, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2 , ErrStat, ErrMsg, 'ReadUSRTimeSeries')
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      ENDIF


   CALL WrScr1(' Reading the user-defined time-series input file "'//TRIM(FileName)//'".' )
   
   IF ( UnEc > 0 )  WRITE (UnEc,'(/,A,/)')  'Data from '//TRIM(TurbSim_Ver%Name)//' user time-series input file "'//TRIM( FileName )//'":'

   
   do i=1,3
      CALL ReadCom( UnIn, FileName, "Header #"//TRIM(Num2Lstr(i))//"for user time-series input", ErrStat2, ErrMsg2, UnEc )   
         CALL SetErrStat(ErrStat2, ErrMsg2 , ErrStat, ErrMsg, 'ReadUSRTimeSeries')
   end do
      
   CALL ReadVar( UnIn, FileName, p%usr%nComp, 'nComp', 'How many velocity components will be input? (1=u component only; 2=u&v components; 3=u,v,and w)', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat(ErrStat2, ErrMsg2 , ErrStat, ErrMsg, 'ReadUSRTimeSeries')
         
   CALL ReadVar( UnIn, FileName, p%usr%nPoints, 'nPoints', 'Number of time series points contained in this file', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat(ErrStat2, ErrMsg2 , ErrStat, ErrMsg, 'ReadUSRTimeSeries')
         
   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF
         
   IF ( p%usr%nPoints < 1 ) THEN
      CALL SetErrStat(ErrID_Fatal, 'user time-series input: nPoints must be at least 1.', ErrStat, ErrMsg, 'ReadUSRTimeSeries')
      CALL Cleanup()
      RETURN
   END IF
      
   CALL AllocAry(p%usr%PointID, p%usr%nPoints, 'PointID', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2 , ErrStat, ErrMsg, 'ReadUSRTimeSeries')
   CALL AllocAry(p%usr%Pointyi, p%usr%nPoints, 'Pointyi', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2 , ErrStat, ErrMsg, 'ReadUSRTimeSeries')
   CALL AllocAry(p%usr%Pointzi, p%usr%nPoints, 'Pointzi', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2 , ErrStat, ErrMsg, 'ReadUSRTimeSeries')

   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF
      
   do i=1,2
      CALL ReadCom( UnIn, FileName, "Point location header #"//TRIM(Num2Lstr(i)), ErrStat2, ErrMsg2, UnEc )   
         CALL SetErrStat(ErrStat2, ErrMsg2 , ErrStat, ErrMsg, 'ReadUSRTimeSeries')
   end do
      
   do i=1,p%usr%nPoints
      CALL ReadAry( UnIn, FileName, TmpAry, 3, "point"//trim(Num2Lstr(i)), "locations of points", ErrStat2, ErrMsg2, UnEc )
         CALL SetErrStat(ErrStat2, ErrMsg2 , ErrStat, ErrMsg, 'ReadUSRTimeSeries')
         
      p%usr%Pointyi = TmpAry(2)
      p%usr%Pointzi = TmpAry(3)         
   end do
      
         
   do i=1,3
      CALL ReadCom( UnIn, FileName, "Time Series header #"//TRIM(Num2Lstr(i)), ErrStat2, ErrMsg2, UnEc )   
         CALL SetErrStat(ErrStat2, ErrMsg2 , ErrStat, ErrMsg, 'ReadUSRTimeSeries')
   end do      
      
   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF      
      
      
   !.......
         ! find out how many rows there are to the end of the file
   p%usr%NTimes   = -1
   ErrStat2 = 0
      
   
   DO WHILE ( ErrStat2 == 0 )
      
     p%usr%NTimes = p%usr%NTimes + 1
     READ(UnIn, *, IOSTAT=ErrStat2) tmpAry(1)
      
   END DO
      
   CALL WrScr( '   Found '//TRIM(Num2LStr(p%usr%NTimes))//' lines of time-series data.' )
   
   IF (p%usr%NTimes < 2) THEN
      CALL SetErrStat(ErrID_Fatal, 'The user time-series input file must contain at least 2 rows of time data.', ErrStat, ErrMsg, 'ReadUSRTimeSeries')
      CALL Cleanup()
      RETURN
   END IF
   
      ! now rewind and skip the first few lines. 
   REWIND( UnIn, IOSTAT=ErrStat2 )  
      IF (ErrStat2 /= 0_IntKi ) THEN
         CALL SetErrStat( ErrID_Fatal, 'Error rewinding file "'//TRIM(FileName)//'".', ErrStat, ErrMsg, 'ReadUSRTimeSeries')   
         CALL Cleanup()
      END IF 

      !IMPORTANT: any changes to the number of lines in the header must be reflected in NumLinesBeforeTS
   DO I=1,NumLinesBeforeTS + p%usr%nPoints        
      READ( UnIn, '(A)', IOSTAT=ErrStat2 ) TmpChar   ! I'm going to ignore this error because we should have caught any issues the first time we read the file.      
   END DO
   
   !.......
   
   if (p%usr%nComp < 1 .OR. p%usr%nComp > 3) then
      CALL SetErrStat( ErrID_Fatal, 'Number of velocity components in file must be 1, 2 or 3.', ErrStat, ErrMsg, 'ReadUSRTimeSeries')   
      CALL Cleanup()
   END IF 
   
   
   
   CALL AllocAry(p%usr%t, p%usr%nTimes,                             't', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2 , ErrStat, ErrMsg, 'ReadUSRTimeSeries')
   CALL AllocAry(p%usr%v, p%usr%nTimes, p%usr%nPoints, p%usr%nComp, 'v', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2 , ErrStat, ErrMsg, 'ReadUSRTimeSeries')
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF
   
      
   DO i=1,p%usr%nTimes
      READ( UnIn, *, IOSTAT=ErrStat2 ) p%usr%t(i), ( (p%usr%v(i,iPoint,iVec), iVec=1,p%usr%nComp), iPoint=1,p%usr%nPoints )
      IF (ErrStat2 /=0) THEN
         CALL SetErrStat( ErrID_Fatal, 'Error reading from time series line '//trim(num2lstr(i))//'.', ErrStat, ErrMsg, 'ReadUSRTimeSeries')
         CALL Cleanup()
         RETURN
      END IF      
   END DO   
   
   IF (UnEc > 0 ) THEN
      FormStr = '('//trim(num2lstr(1+p%usr%nComp*p%usr%nPoints))//'(F13.4," "))'
      DO i=1,p%usr%nTimes
         WRITE( UnEc, FormStr) p%usr%t(i), ( (p%usr%v(i,iPoint,iVec), iVec=1,p%usr%nComp), iPoint=1,p%usr%nPoints )
      END DO
   END IF      

   
   !.........................................................
   ! a little bit of error checking:
   !.........................................................

   !bjj: verify that the locations are okay; for now, we're going to make sure they are unique and that the z values are in increasing order.
   do i=2,p%usr%nPoints
      do j=1,i-1
         IF ( EqualRealNos( p%usr%Pointyi(i), p%usr%Pointyi(j) ) .AND. EqualRealNos( p%usr%Pointzi(i), p%usr%Pointzi(j) ) ) THEN
            CALL SetErrStat(ErrID_Fatal, 'Locations of points specified in the user time-series input file must be unique.', ErrStat, ErrMsg, 'ReadUSRTimeSeries')
            CALL Cleanup()
            RETURN
         END IF
      end do
      
      !bjj: fix this in the future. Currently the interpolation routine won't work if z is not ordered properly. Also, interpolation doesn't take y into account, so we may want to fix that.
      IF ( p%usr%Pointzi(i) < p%usr%Pointzi(i-1) ) THEN
         CALL SetErrStat(ErrID_Fatal, 'The current implementation of user time-series input requires that the points be entered in the order of increasing height.', ErrStat, ErrMsg, 'ReadUSRTimeSeries')
         CALL Cleanup()
         RETURN
      END IF      
   end do
   
   
      ! a little bit of error checking:
   DO i = 2,p%usr%nTimes
      IF (.NOT. EqualRealNos( p%usr%t(i-1) + p%grid%TimeStep, p%usr%t(i) ) ) THEN
         call SetErrStat(ErrID_Fatal, 'the delta time in the file must be constant and must be equal to input file variable TimeStep.', ErrStat, ErrMsg, 'ReadUSRTimeSeries')
         EXIT
      END IF
   END DO
   
   
   
   CALL Cleanup()
   RETURN
   
CONTAINS
!...............................................
   SUBROUTINE Cleanup()
   
      CLOSE( UnIn )
   
   
   END SUBROUTINE Cleanup
!...............................................
END SUBROUTINE ReadUSRTimeSeries 
!=======================================================================
SUBROUTINE ReadCVarDefault ( UnIn, Fil, CharVar, VarName, VarDescr, UnEc, Def, ErrStat, ErrMsg, IGNORE )


      ! This routine reads a single character variable from the next line of the input file.
      ! The input is allowed to be "default"

      ! Argument declarations:


   INTEGER, INTENT(IN)          :: UnIn                                            ! I/O unit for input file.
   INTEGER, INTENT(IN)          :: UnEc                                           ! I/O unit for echo/summary file.
   INTEGER(IntKi), INTENT(OUT)  :: ErrStat                                         ! Error status; if present, program does not abort on error
   CHARACTER(*),   INTENT(OUT)  :: ErrMsg                                          ! Error message

   LOGICAL, INTENT(INOUT)       :: Def                                             ! - on input whether or not to use the default - on output, whether a default was used
   LOGICAL, INTENT(IN), OPTIONAL:: IGNORE                                          ! whether to ignore this input

   CHARACTER(250)               :: CharLine                                        ! Character string being read.
   CHARACTER(*), INTENT(INOUT)  :: CharVar                                         ! Character variable being read.
   CHARACTER( *), INTENT(IN)    :: Fil                                             ! Name of the input file.
   CHARACTER( *), INTENT(IN)    :: VarDescr                                        ! Text string describing the variable.
   CHARACTER( *), INTENT(IN)    :: VarName                                         ! Text string containing the variable name.

      ! Local declarations:


   CALL ReadVar( UnIn, Fil, CharLine, VarName, VarDescr, ErrStat, ErrMsg, UnEc)

   IF ( PRESENT(IGNORE) ) THEN
      IF ( IGNORE ) THEN
         Def = .TRUE.
         RETURN
      ENDIF

   ENDIF

   CALL Conv2UC( CharLine )

   IF ( TRIM(CharLine) == 'DEFAULT' ) THEN

      CALL WrScr ( '    A default value will be used for '//TRIM(VarName)//'.' )
      Def = .TRUE.

   ELSE

      CharVar = CharLine
      Def = .FALSE.

   ENDIF

   RETURN
END SUBROUTINE ReadCVarDefault ! ( UnIn, Fil, RealVar, VarName, VarDescr )
!=======================================================================
SUBROUTINE ReadRAryDefault ( UnIn, Fil, RealAry, VarName, VarDescr, UnEc, Def, ErrStat, ErrMsg, IGNORE )

      ! This routine reads a real array from the next line of the input file.
      ! The input is allowed to be "default"

      ! Argument declarations:

   REAL(ReKi), INTENT(INOUT)    :: RealAry (:)                                     ! Real variable being read.

   INTEGER, INTENT(IN)          :: UnIn                                            ! I/O unit for input file.
   INTEGER, INTENT(IN)          :: UnEc                                           ! I/O unit for echo/summary file.
   INTEGER(IntKi), INTENT(OUT)  :: ErrStat                                         ! Error status; if present, program does not abort on error
   CHARACTER(*),   INTENT(OUT)  :: ErrMsg                                          ! Error message

   LOGICAL, INTENT(INOUT)       :: Def                                             ! - on input whether or not to use the default - on output, whether a default was used
   LOGICAL, INTENT(IN), OPTIONAL:: IGNORE                                          ! whether or not to ignore this input

   CHARACTER(250)               :: CharLine                                        ! Character string being read.
   CHARACTER( *), INTENT(IN)    :: Fil                                             ! Name of the input file.
   CHARACTER( *), INTENT(IN)    :: VarDescr                                        ! Text string describing the variable.
   CHARACTER( *), INTENT(IN)    :: VarName                                         ! Text string containing the variable name.

      ! Local declarations:

   INTEGER                      :: IOS                                             ! I/O status returned from the read statement.


   CALL ReadVar( UnIn, Fil, CharLine, VarName, VarDescr, ErrStat, ErrMsg, UnEc)  !Maybe I should read this in explicitly...

   IF ( PRESENT(IGNORE) ) THEN
      IF ( IGNORE ) THEN
         Def = .TRUE.
         RETURN
      ENDIF

   ENDIF

   CALL Conv2UC( CharLine )

   IF ( TRIM(CharLine) == 'DEFAULT' ) THEN

      CALL WrScr ( '    A default value will be used for '//TRIM(VarName)//'.' )
      Def = .TRUE.

   ELSE

      IF ( INDEX( CharLine(1:1), 'TF') > 0 ) THEN    ! We don't want 'T' or 'F' read as -1 or 0.
         CALL WrScr1 ( ' Invalid numerical input for "'//TRIM( VarName )//'".' )
      ENDIF

      READ (CharLine,*,IOSTAT=IOS)  RealAry

      IF (IOS /=0) THEN
         READ (CharLine,*,IOSTAT=IOS)  RealAry(1)  ! Try reading only the first element
      ENDIF

      CALL CheckIOS ( IOS, Fil, VarName, NumType )
      Def = .FALSE.

   ENDIF


   RETURN

END SUBROUTINE ReadRAryDefault
!=======================================================================
SUBROUTINE ReadRVarDefault ( UnIn, Fil, RealVar, VarName, VarDescr, UnEc, Def, ErrStat, ErrMsg, IGNORE, IGNORESTR )

      ! This routine reads a single real variable from the next line of the input file.
      ! The input is allowed to be "default"

      ! Argument declarations:

   REAL(ReKi), INTENT(INOUT)    :: RealVar                                         ! Real variable being read.

   INTEGER, INTENT(IN)          :: UnIn                                            ! I/O unit for input file.
   INTEGER, INTENT(IN)          :: UnEc                                            ! I/O unit for echo/summary file.

   INTEGER(IntKi), INTENT(OUT)  :: ErrStat                                         ! Error status; if present, program does not abort on error
   CHARACTER(*),   INTENT(OUT)  :: ErrMsg                                          ! Error message
      
   LOGICAL, INTENT(INOUT)       :: Def                                             ! - on input whether or not to use the default - on output, whether a default was used
   LOGICAL, INTENT(IN), OPTIONAL:: IGNORE                                          ! whether or not to ignore this input
   LOGICAL, INTENT(OUT),OPTIONAL:: IGNORESTR                                       ! whether or not user requested to ignore this input

   CHARACTER(250)               :: CharLine                                        ! Character string being read.
   CHARACTER( *), INTENT(IN)    :: Fil                                             ! Name of the input file.
   CHARACTER( *), INTENT(IN)    :: VarDescr                                        ! Text string describing the variable.
   CHARACTER( *), INTENT(IN)    :: VarName                                         ! Text string containing the variable name.

      ! Local declarations:

   INTEGER                      :: IOS                                             ! I/O status returned from the read statement.


   CALL ReadVar( UnIn, Fil, CharLine, VarName, VarDescr, ErrStat, ErrMsg, UnEc )

   IF ( PRESENT(IGNORE) ) THEN

      IF ( IGNORE ) THEN
         Def = .TRUE.
         RETURN
      ENDIF

   ENDIF


   CALL Conv2UC( CharLine )

   IF ( PRESENT(IGNORESTR) ) THEN
      IF ( TRIM( CharLine ) == 'NONE' ) THEN
         IGNORESTR = .TRUE.
         Def = .TRUE.
         RealVar = 0.0  ! This is set for the Reynolds stress inputs, but if IGNORESTR is used for other inputs, it may need to be changed
         RETURN
      ENDIF
   ENDIF

   IF ( TRIM(CharLine) == 'DEFAULT' ) THEN

      CALL WrScr ( '    A default value will be used for '//TRIM(VarName)//'.' )

      IF ( PRESENT(IGNORESTR) ) THEN
         IF ( IGNORESTR ) THEN  !We've told it to ignore this input, by default
            RealVar = 0.0  ! This is set for the Reynolds stress inputs, but if IGNORESTR is used for other inputs, it may need to be changed
            RETURN
         ENDIF
      ENDIF

      Def = .TRUE.

   ELSE

      IF ( INDEX( CharLine(1:1), 'TF') > 0 ) THEN    ! We don't want 'T' or 'F' read as -1 or 0.
         CALL WrScr1 ( ' Invalid numerical input for "'//TRIM( VarName )//'".' )
      ENDIF

      READ (CharLine,*,IOSTAT=IOS)  RealVar

      CALL CheckIOS ( IOS, Fil, VarName, NumType )

      Def = .FALSE.

      IF ( PRESENT(IGNORESTR) ) THEN
         IGNORESTR = .FALSE.
      ENDIF

   ENDIF


   RETURN
END SUBROUTINE ReadRVarDefault
!=======================================================================
!> This routine writes the velocity grid to a binary file.
!! The file has a .wnd extension; scaling information is written in a summary
!! file. A tower file with extension .twr is generated if requested, too.
SUBROUTINE WrBinBLADED(p, V, USig, VSig, WSig, US, ErrStat, ErrMsg)

   IMPLICIT                    NONE

   TYPE(TurbSim_ParameterType),     INTENT(IN)    :: p                             !< TurbSim's parameters
   REAL(ReKi),                      INTENT(IN)    :: V     (:,:,:)                 !< An array containing the summations of the rows of H (NumSteps,NPoints,3).
   REAL(ReKi),                      INTENT(IN)    :: USig                          !< Standard deviation of U
   REAL(ReKi),                      INTENT(IN)    :: VSig                          !< Standard deviation of V
   REAL(ReKi),                      INTENT(IN)    :: WSig                          !< Standard deviation of W
   INTEGER(IntKi)                 , INTENT(IN)    :: US                            !< unit number of file in which to print a summary of the scaling used. If < 1, will not print summary.
   INTEGER(IntKi),                  intent(  out) :: ErrStat                       !< Error level
   CHARACTER(*),                    intent(  out) :: ErrMsg                        !< Message describing error


   REAL(ReKi)                  :: NewUSig                       ! Value of USig that will be used to scale values in the file
   
   REAL(ReKi)                  :: U_C1                          ! Scale for converting BLADED U data
   REAL(ReKi)                  :: U_C2                          ! Offset for converting BLADED U data
   REAL(ReKi)                  :: V_C                           ! Scale for converting BLADED V data
   REAL(ReKi)                  :: W_C                           ! Scale for converting BLADED W data
   REAL(ReKi)                  :: TI(3)                         ! Turbulence intensity for scaling data
   REAL(ReKi)                  :: TmpU                          ! Max value of |V(:,:,1)-UHub|

   INTEGER(B4Ki)               :: CFirst
   INTEGER(B4Ki)               :: CLast
   INTEGER(B4Ki)               :: CStep
   INTEGER(B4Ki)               :: II
   INTEGER(B4Ki)               :: IT
   INTEGER(B4Ki)               :: IY
   INTEGER(B4Ki)               :: IZ

   INTEGER(B4Ki)               :: IP
   INTEGER(B2Ki)               :: TmpVarray(3*p%grid%NumGrid_Y*p%grid%NumGrid_Z) ! This array holds the normalized velocities before being written to the binary file
   INTEGER(B2Ki),ALLOCATABLE   :: TmpTWRarray(:)                   ! This array holds the normalized tower velocities   
   
   INTEGER                     :: AllocStat
   INTEGER                     :: UBFFW                                ! I/O unit for BLADED FF data (*.wnd file).
   INTEGER                     :: UATWR                                ! I/O unit for AeroDyn tower data (*.twr file).

   CHARACTER(200)               :: FormStr                                  ! String used to store format specifiers.


   
      ! We need to take into account the shear across the grid in the sigma calculations for scaling the data, 
      ! and ensure that 32.767*Usig >= |V-UHub| so that we don't get values out of the range of our scaling values
      ! in this BLADED-style binary output.  TmpU is |V-UHub|
   TmpU    = MAX( ABS(MAXVAL(V(:,:,1))-p%UHub), ABS(MINVAL(V(:,:,1))-p%UHub) )  !Get the range of wind speed values for scaling in BLADED-format .wnd files         
   NewUSig = MAX(USig,0.05*TmpU)
   
   
      ! Put normalizing factors into the summary file.  The user can use them to
      ! tell a simulation program how to rescale the data.

   TI(1)  = MAX(100.0*Tolerance, NewUSig) / p%UHub
   TI(2)  = MAX(100.0*Tolerance, VSig)    / p%UHub
   TI(3)  = MAX(100.0*Tolerance, WSig)    / p%UHub

   WRITE (US,"(//,'Normalizing Parameters for Binary Data (approximate statistics):',/)")

   FormStr = "(3X,A,' =',F9.4,A)"
   WRITE (US,FormStr)  'UBar ', p%UHub,      ' m/s'
   WRITE (US,FormStr)  'TI(u)', 100.0*TI(1), ' %'
   WRITE (US,FormStr)  'TI(v)', 100.0*TI(2), ' %'
   WRITE (US,FormStr)  'TI(w)', 100.0*TI(3), ' %'

   WRITE (US,'()')
   WRITE (US,FormStr)  'Height Offset', ( p%grid%HubHt - p%grid%GridHeight / 2.0 - p%grid%Z(1) ),  ' m'
   WRITE (US,FormStr)  'Grid Base    ', p%grid%Z(1),                                               ' m'
   
   WRITE (US,'()'   )
   WRITE (US,'( A)' ) 'Creating a PERIODIC output file.'
 
      ! Calculate some numbers for normalizing the data.

   U_C1 = 1000.0/( p%UHub*TI(1) )
   U_C2 = 1000.0/TI(1)
   V_C  = 1000.0/( p%UHub*TI(2) )
   W_C  = 1000.0/( p%UHub*TI(3) )


   IF ( p%WrFile(FileExt_WND) )  THEN

      CALL GetNewUnit( UBFFW )
      CALL OpenBOutFile ( UBFFW, TRIM(p%RootName)//'.wnd', ErrStat, ErrMsg )
      IF (ErrStat >= AbortErrLev) RETURN
      
      CALL WrScr ( ' Generating BLADED binary time-series file "'//TRIM( p%RootName )//'.wnd"' )

               ! Put header information into the binary data file.

      WRITE (UBFFW)   INT(  -99                                   , B2Ki )     ! -99 = New Bladed format
      WRITE (UBFFW)   INT(    4                                   , B2Ki )     ! 4 = improved von karman (but needed for next 7 inputs)
      WRITE (UBFFW)   INT(    3                                   , B4Ki )     ! 3 = number of wind components
      WRITE (UBFFW)  REAL( p%met%Latitude                         , SiKi )     ! Latitude (degrees)
      WRITE (UBFFW)  REAL( p%met%Z0                               , SiKi )     ! Roughness length (m)
      WRITE (UBFFW)  REAL( p%grid%Zbottom + p%grid%GridHeight/2.0 , SiKi )     ! Reference Height (m) ( Z(1) + GridHeight / 2.0 ) !This is the vertical center of the grid
      WRITE (UBFFW)  REAL( 100.0*TI(1)                            , SiKi )     ! Longitudinal turbulence intensity (%)
      WRITE (UBFFW)  REAL( 100.0*TI(2)                            , SiKi )     ! Lateral turbulence intensity (%)
      WRITE (UBFFW)  REAL( 100.0*TI(3)                            , SiKi )     ! Vertical turbulence intensity (%)

      WRITE (UBFFW)  REAL( p%grid%GridRes_Z                       , SiKi )     ! grid spacing in vertical direction, in m
      WRITE (UBFFW)  REAL( p%grid%GridRes_Y                       , SiKi )     ! grid spacing in lateral direction, in m
      WRITE (UBFFW)  REAL( p%grid%TimeStep*p%UHub                 , SiKi )     ! grid spacing in longitudinal direciton, in m
      WRITE (UBFFW)   INT( p%grid%NumOutSteps/2                   , B4Ki )     ! half the number of points in alongwind direction
      WRITE (UBFFW)  REAL( p%UHub                                 , SiKi )     ! the mean wind speed in m/s
      WRITE (UBFFW)  REAL( 0                                      , SiKi )     ! the vertical length scale of the longitudinal component in m
      WRITE (UBFFW)  REAL( 0                                      , SiKi )     ! the lateral length scale of the longitudinal component in m
      WRITE (UBFFW)  REAL( 0                                      , SiKi )     ! the longitudinal length scale of the longitudinal component in m
      WRITE (UBFFW)   INT( 0                                      , B4Ki )     ! an unused integer
      WRITE (UBFFW)   INT( p%RNG%RandSeed(1)                      , B4Ki )     ! the random number seed
      WRITE (UBFFW)   INT( p%grid%NumGrid_Z                       , B4Ki )     ! the number of grid points vertically
      WRITE (UBFFW)   INT( p%grid%NumGrid_Y                       , B4Ki )     ! the number of grid points laterally
      WRITE (UBFFW)   INT( 0                                      , B4Ki )     ! the vertical length scale of the lateral component, not used
      WRITE (UBFFW)   INT( 0                                      , B4Ki )     ! the lateral length scale of the lateral component, not used
      WRITE (UBFFW)   INT( 0                                      , B4Ki )     ! the longitudinal length scale of the lateral component, not used
      WRITE (UBFFW)   INT( 0                                      , B4Ki )     ! the vertical length scale of the vertical component, not used
      WRITE (UBFFW)   INT( 0                                      , B4Ki )     ! the lateral length scale of the vertical component, not used
      WRITE (UBFFW)   INT( 0                                      , B4Ki )     ! the longitudinal length scale of the vertical component, not used


         ! Compute parameters for ordering output for FF AeroDyn files. (This is for BLADED compatibility.)

      IF ( p%grid%Clockwise )  THEN
         CFirst = p%grid%NumGrid_Y
         CLast  = 1
         CStep  = -1
      ELSE
         CFirst = 1
         CLast  = p%grid%NumGrid_Y
         CStep  = 1
      ENDIF


         ! Loop through time.

      DO IT=1,p%grid%NumOutSteps  !Use only the number of timesteps requested originally

            ! Write out grid data in binary form.
         IP = 1
         DO IZ=1,p%grid%NumGrid_Z
            DO IY=CFirst,CLast,CStep

               II = ( IZ - 1 )*p%grid%NumGrid_Y + IY

               TmpVarray(IP)   = NINT( U_C1*V(IT,p%grid%GridPtIndx(II),1) - U_C2, B2Ki )  ! Put the data into a temp array so that the WRITE() command works faster
               TmpVarray(IP+1) = NINT( V_C *V(IT,p%grid%GridPtIndx(II),2)       , B2Ki )
               TmpVarray(IP+2) = NINT( W_C *V(IT,p%grid%GridPtIndx(II),3)       , B2Ki )

               IP = IP + 3;
            ENDDO ! IY
         ENDDO ! IZ

         WRITE ( UBFFW )  TmpVarray ! bjj: We cannot write the array including time because of stack overflow errors.. otherwise use compile option to put this on the heap instead of the stack?

      ENDDO ! IT

      CLOSE ( UBFFW )

      !.......................................................
      ! Now write tower data file if necessary:
      !.......................................................      
      
      IF ( p%WrFile(FileExt_TWR) ) THEN
                           
         CALL GetNewUnit( UATWR )
         CALL OpenBOutFile ( UATWR, TRIM( p%RootName )//'.twr', ErrStat, ErrMsg )
         IF (ErrStat >= AbortErrLev) RETURN
         
         CALL WrScr ( ' Generating tower binary time-series file "'//TRIM( p%RootName )//'.twr"' )


         WRITE (UATWR)  REAL( p%grid%GridRes_Z       ,   SiKi )         ! grid spacing in vertical direction, in m
         WRITE (UATWR)  REAL( p%grid%TimeStep*p%UHub ,   SiKi )         ! grid spacing in longitudinal direciton, in m
         WRITE (UATWR)  REAL( p%grid%ZBottom         ,   SiKi )         ! The vertical location of the highest tower grid point in m
         WRITE (UATWR)   INT( p%grid%NumOutSteps     ,   B4Ki )         ! The number of points in alongwind direction
         WRITE (UATWR)   INT( SIZE(p%grid%TwrPtIndx) ,   B4Ki )         ! the number of grid points vertically
         WRITE (UATWR)  REAL( p%UHub ,                   SiKi )         ! the mean wind speed in m/s
         WRITE (UATWR)  REAL( 100.0*TI(1),               SiKi )         ! Longitudinal turbulence intensity
         WRITE (UATWR)  REAL( 100.0*TI(2),               SiKi )         ! Lateral turbulence intensity
         WRITE (UATWR)  REAL( 100.0*TI(3),               SiKi )         ! Vertical turbulence intensity


         ALLOCATE ( TmpTWRarray( 3*SIZE(p%grid%TwrPtIndx) ) , STAT=AllocStat ) 

         IF ( AllocStat /= 0 )  THEN
            CALL TS_Abort ( 'Error allocating memory for temporary tower wind speed array.' )
         ENDIF


         DO IT=1,p%grid%NumOutSteps

            IP = 1
            DO II = 1,SIZE(p%grid%TwrPtIndx)
               TmpTWRarray(IP  ) = NINT( U_C1*V(IT,p%grid%TwrPtIndx(II),1) - U_C2 , B2Ki )
               TmpTWRarray(IP+1) = NINT( V_C *V(IT,p%grid%TwrPtIndx(II),2)        , B2Ki )
               TmpTWRarray(IP+2) = NINT( W_C *V(IT,p%grid%TwrPtIndx(II),3)        , B2Ki )

               IP = IP + 3
            ENDDO    ! II

            WRITE (UATWR) TmpTWRarray(:)

         ENDDO ! IT

         CLOSE ( UATWR )

         IF ( ALLOCATED( TmpTWRarray) ) DEALLOCATE( TmpTWRarray )

      ENDIF !WrADWTR

   ENDIF ! p%WrFile(FileExt_WND)
   

END SUBROUTINE WrBinBLADED
!=======================================================================
!> This routine writes the velocity grid to a binary file.
!! The file has a .bts extension.
!=======================================================================
SUBROUTINE WrBinTURBSIM(p, V, ErrStat, ErrMsg)

   IMPLICIT                  NONE

      ! passed variables
   TYPE(TurbSim_ParameterType),     INTENT(IN)     :: p                             !< TurbSim's parameters
   REAL(ReKi),                      INTENT(IN)     :: V     (:,:,:)                 !< An array containing the summations of the rows of H (NumSteps,NPoints,3).
   INTEGER(IntKi),                  intent(  out)  :: ErrStat                       !< Error level
   CHARACTER(*),                    intent(  out)  :: ErrMsg                        !< Message describing error

      ! local variables
   REAL(SiKi), PARAMETER     :: IntMax   =  32767.0
   REAL(SiKi), PARAMETER     :: IntMin   = -32768.0
   REAL(SiKi), PARAMETER     :: IntRng   = IntMax - IntMin ! Max Range of 2-byte integer

   REAL(SiKi)                :: UOff                       ! Offset for the U component
   REAL(SiKi)                :: UScl                       ! Slope  for the U component
   REAL(ReKi)                :: VMax(3)                    ! Maximum value of the 3 wind components
   REAL(ReKi)                :: VMin(3)                    ! Minimum value of the 3 wind components
   REAL(SiKi)                :: VOff                       ! Offset for the V component
   REAL(SiKi)                :: VScl                       ! Slope  for the V component
   REAL(SiKi)                :: WOff                       ! Offset for the W component
   REAL(SiKi)                :: WScl                       ! Slope  for the W component

   INTEGER, PARAMETER        :: DecRound  = 3              ! Number of decimal places to round to
   INTEGER(B2Ki)             :: FileID                     ! File ID, determines specific output format (if periodic or not)
   INTEGER                   :: IC                         ! counter for the velocity component of V
   INTEGER                   :: II                         ! counter for the point on the grid/tower
   INTEGER                   :: IT                         ! counter for the timestep
   INTEGER(B4Ki)             :: LenDesc                    ! Length of the description string
   INTEGER(B4Ki)             :: NumGrid                    ! Number of points on the grid
   INTEGER(B4Ki)             :: NumTower                   ! Number of points on the tower

   INTEGER(B4Ki)             :: IP
   INTEGER(B2Ki),ALLOCATABLE :: TmpVarray(:)               ! This array holds the normalized velocities before being written to the binary file

   INTEGER                   :: AllocStat
   INTEGER                   :: UAFFW                      ! I/O unit for AeroDyn FF data (*.bts file).



      ! Set the file format ID
      
   IF ( p%grid%Periodic ) THEN 
      FileID = 8
   ELSE
      FileID = 7
   END IF      


      ! Find the range of our velocity

   DO IC=1,3

         ! Initialize the Min/Max values

      VMin(IC) = V(1,1,IC)
      VMax(IC) = V(1,1,IC)

      DO II=1,p%grid%NPoints   ! Let's check all of the points
         DO IT=1,p%grid%NumOutSteps  ! Use only the number of timesteps requested originally

            IF ( V(IT,II,IC) > VMax(IC) ) THEN

               VMax(IC) = V(IT,II,IC)

            ELSEIF ( V(IT,II,IC) < VMin(IC) ) THEN

               VMin(IC) = V(IT,II,IC)

            ENDIF

         ENDDO !IT
      ENDDO !II

   ENDDO !IC


      ! Calculate the scaling parameters for each component


   IF ( VMax(1) == VMin(1) ) THEN
      UScl = 1
   ELSE
      UScl = IntRng/REAL( VMax(1) - VMin(1) , SiKi )
   ENDIF

   IF ( VMax(2) == VMin(2) ) THEN
      VScl = 1
   ELSE
      VScl = IntRng/REAL( VMax(2) - VMin(2) , SiKi )
   ENDIF

   IF ( VMax(3) == VMin(3) ) THEN
      WScl = 1
   ELSE
      WScl = IntRng/REAL( VMax(3) - VMin(3) , SiKi )
   ENDIF


   UOff = IntMin - UScl*REAL( VMin(1)    , SiKi )
   VOff = IntMin - VScl*REAL( VMin(2)    , SiKi )
   WOff = IntMin - WScl*REAL( VMin(3)    , SiKi )


      ! Find the first tower point

   NumGrid  = SIZE(p%grid%GridPtIndx) !p%grid%NumGrid_Y*p%grid%NumGrid_Z

   IF ( p%WrFile(FileExt_TWR) ) THEN

      NumTower = SIZE(p%grid%TwrPtIndx)

   ELSE

      NumTower = 0

   ENDIF


   LenDesc = LEN_TRIM( p%DescStr )             ! Length of the string that contains program name, version, date, and time

   CALL WrScr ( ' Generating AeroDyn binary time-series file "'//TRIM( p%RootName )//'.bts"' )
   
   CALL GetNewUnit(UAFFW, ErrStat, ErrMsg)
   CALL OpenBOutFile ( UAFFW, TRIM(p%RootName)//'.bts', ErrStat, ErrMsg )
   IF (ErrStat >= AbortErrLev) RETURN



      ! Write the header

   WRITE (UAFFW)   INT( FileID                    , B2Ki )          ! TurbSim format (7=not PERIODIC, 8=PERIODIC)

   WRITE (UAFFW)   INT( p%grid%NumGrid_Z          , B4Ki )          ! the number of grid points vertically
   WRITE (UAFFW)   INT( p%grid%NumGrid_Y          , B4Ki )          ! the number of grid points laterally
   WRITE (UAFFW)   INT( NumTower                  , B4Ki )          ! the number of tower points
   WRITE (UAFFW)   INT( p%grid%NumOutSteps        , B4Ki )          ! the number of time steps

   WRITE (UAFFW)  REAL( p%grid%GridRes_Z          , SiKi )          ! grid spacing in vertical direction, in m
   WRITE (UAFFW)  REAL( p%grid%GridRes_Y          , SiKi )          ! grid spacing in lateral direction, in m
   WRITE (UAFFW)  REAL( p%grid%TimeStep           , SiKi )          ! grid spacing in delta time, in m/s
   WRITE (UAFFW)  REAL( p%UHub                    , SiKi )          ! the mean wind speed in m/s at hub height
   WRITE (UAFFW)  REAL( p%grid%HubHt              , SiKi )          ! the hub height, in m
   WRITE (UAFFW)  REAL( p%grid%Z(1)               , SiKi )          ! the height of the grid bottom, in m

   WRITE (UAFFW)  REAL( UScl                      , SiKi )          ! the U-component slope for scaling
   WRITE (UAFFW)  REAL( UOff                      , SiKi )          ! the U-component offset for scaling
   WRITE (UAFFW)  REAL( VScl                      , SiKi )          ! the V-component slope for scaling
   WRITE (UAFFW)  REAL( VOff                      , SiKi )          ! the V-component offset for scaling
   WRITE (UAFFW)  REAL( WScl                      , SiKi )          ! the W-component slope for scaling
   WRITE (UAFFW)  REAL( WOff                      , SiKi )          ! the W-component offset for scaling

   WRITE (UAFFW)   INT( LenDesc                   , B4Ki )          ! the number of characters in the string, max 200

   DO II=1,LenDesc

      WRITE (UAFFW)  INT( IACHAR( p%DescStr(II:II) ), B1Ki )   ! converted ASCII characters

   ENDDO

   ALLOCATE ( TmpVarray( 3*(NumGrid + NumTower) ) , STAT=AllocStat )

   IF ( AllocStat /= 0 )  THEN
      CALL TS_Abort ( 'Error allocating memory for temporary wind speed array.' )
   ENDIF

      ! Loop through time.

   DO IT=1,p%grid%NumOutSteps  !Use only the number of timesteps requested originally

         ! Write out grid data in binary form. II = (IZ - 1)*NumGrid_Y + IY, IY varies most rapidly

      IP = 1

      DO II=1,NumGrid

         TmpVarray(IP)   =  NINT( Max( Min( REAL(UScl*V(IT,p%grid%GridPtIndx(II),1) + UOff, SiKi), IntMax ),IntMin) , B2Ki )
         TmpVarray(IP+1) =  NINT( Max( Min( REAL(VScl*V(IT,p%grid%GridPtIndx(II),2) + VOff, SiKi), IntMax ),IntMin) , B2Ki )
         TmpVarray(IP+2) =  NINT( Max( Min( REAL(WScl*V(IT,p%grid%GridPtIndx(II),3) + Woff, SiKi), IntMax ),IntMin) , B2Ki )

         IP = IP + 3
      ENDDO ! II


      IF ( p%WrFile(FileExt_TWR) ) THEN

            ! Write out the tower data in binary form

         DO II=1,SIZE(p%grid%TwrPtIndx)
                ! Values of tower data
            TmpVarray(IP)   =  NINT( Max( Min( REAL(UScl*V(IT,p%grid%TwrPtIndx(II),1) + UOff, SiKi), IntMax ),IntMin) , B2Ki )
            TmpVarray(IP+1) =  NINT( Max( Min( REAL(VScl*V(IT,p%grid%TwrPtIndx(II),2) + VOff, SiKi), IntMax ),IntMin) , B2Ki )
            TmpVarray(IP+2) =  NINT( Max( Min( REAL(WScl*V(IT,p%grid%TwrPtIndx(II),3) + Woff, SiKi), IntMax ),IntMin) , B2Ki )

            IP = IP + 3
         ENDDO ! II

      ENDIF

      WRITE ( UAFFW ) TmpVarray(:)
   ENDDO ! IT

   CLOSE ( UAFFW )

   IF ( ALLOCATED( TmpVarray ) ) DEALLOCATE( TmpVarray )


END SUBROUTINE WrBinTURBSIM
!=======================================================================
SUBROUTINE WrFormattedFF(RootName, p_grid, UHub, V)

USE TurbSim_Types

IMPLICIT                NONE

CHARACTER(*),                    intent(in   ) :: RootName             ! Rootname of output file
TYPE(Grid_ParameterType),        INTENT(IN)    :: p_grid
REAL(ReKi),                      INTENT(IN)    :: UHub                 ! The steady hub-height velocity
REAL(ReKi),                      INTENT(IN)    :: V       (:,:,:)      ! The Velocities to write to a file


REAL(ReKi), ALLOCATABLE      :: ZRow       (:)                           ! The horizontal locations of the grid points (NumGrid_Y) at each height.

INTEGER                      :: AllocStat
INTEGER                      :: II
INTEGER                      :: IT
INTEGER                      :: IVec
INTEGER                      :: IY
INTEGER                      :: IZ

CHARACTER(200)               :: FormStr5                                 ! String used to store format specifiers.

INTEGER                      :: UFFF                                     ! I/O unit for formatted FF data.


   FormStr5 = "(1X,"//trim(num2lstr(max(p_grid%NumGrid_Z,p_grid%NumGrid_Y)))//"(F8.3),:)"

      ! Allocate the array of wind speeds.

   ALLOCATE ( ZRow(p_grid%NumGrid_Y) , STAT=AllocStat )

   IF ( AllocStat /= 0 )  THEN
      CALL TS_Abort ( 'Error allocating memory for array of wind speeds.' )
   ENDIF

   CALL GetNewUnit(UFFF)

   DO IVec=1,3

      CALL WrScr ( ' Generating full-field formatted file "'//TRIM(RootName)//'.'//Comp(IVec)//'".' )
      CALL OpenFOutFile ( UFFF, TRIM( RootName )//'.'//Comp(IVec) )


         ! Create file header.

      WRITE (UFFF,"( / 'This full-field turbulence file was generated by ' , A , ' on ' , A , ' at ' , A , '.' / )" )  TRIM(GetNVD(TurbSim_Ver)), CurDate(), CurTime()

      WRITE (UFFF,"( ' | ', A,'-comp |  Y  x  Z  | Grid Resolution (Y x Z) | Time-step | Hub Elev | Mean U |')")  Comp(IVec)

      WRITE (UFFF,"(I14,I6,F11.3,F11.3,F15.3,F11.2,F10.2)")  p_grid%NumGrid_Y, p_grid%NumGrid_Z, p_grid%GridRes_Y, p_grid%GridRes_Z, p_grid%TimeStep, p_grid%HubHt, UHub
      WRITE (UFFF,"(/,' Z Coordinates (m):')")
      WRITE (UFFF,FormStr5)  ( p_grid%Z(p_grid%GridPtIndx(IZ))-p_grid%HubHt, IZ=1,p_grid%NPoints,p_grid%NumGrid_Y )
      WRITE (UFFF,"(/,' Y Coordinates (m):')")
      WRITE (UFFF,FormStr5)  ( p_grid%Y(p_grid%GridPtIndx(IY)), IY=1,p_grid%NumGrid_Y )

         ! Write out elapsed time & hub-level value before component grid.

      DO IT=1,p_grid%NumOutSteps

         WRITE(UFFF,"(/,1X,2(F8.3))")  p_grid%TimeStep*( IT - 1 ), V(IT,p_grid%HubIndx,IVec)

         DO IZ=1,p_grid%NumGrid_Z  ! From the top to the bottom

            II = ( p_grid%NumGrid_Z - IZ )*p_grid%NumGrid_Y

            DO IY=1,p_grid%NumGrid_Y  ! From the left to the right
               ZRow(IY) = V(IT, p_grid%GridPtIndx(II+IY) ,IVec)
            ENDDO ! IY

            WRITE (UFFF,FormStr5)  ( ZRow(IY), IY=1,p_grid%NumGrid_Y )

         ENDDO ! IZ

      ENDDO ! IT

      CLOSE ( UFFF )

   ENDDO ! IVec

      ! Deallocate the array of wind speeds.

   IF ( ALLOCATED( ZRow ) )  DEALLOCATE( ZRow )

END SUBROUTINE WrFormattedFF
!=======================================================================

SUBROUTINE WrSum_UserInput( p_met, US  )


   TYPE(Meteorology_ParameterType), INTENT(IN)  :: p_met                       ! parameters for TurbSim 
   
   INTEGER, INTENT(IN) :: US
   integer             :: i 


   IF ( p_met%NumUSRz > 0 ) THEN
      WRITE (US,"( // 'User-Defined profiles:' / )")
   
      IF ( ALLOCATED( p_met%USR_L ) ) THEN
         WRITE (US,"(A)") '  Height   Wind Speed   Horizontal Angle   u Std. Dev.   v Std. Dev.   w Std. Dev.   Length Scale'
         WRITE (US,"(A)") '   (m)        (m/s)          (deg)            (m/s)         (m/s)         (m/s)          (m)     '
         WRITE (US,"(A)") '  ------   ----------   ----------------   -----------   -----------   -----------   ------------'
      
         DO I=p_met%NumUSRz,1,-1
            WRITE (US,"( 1X,F7.2, 2X,F9.2,2X, 3X,F10.2,6X, 3(4X,F7.2,3X), 3X,F10.2 )")  &
                                p_met%USR_Z(I), p_met%USR_U(I), p_met%USR_WindDir(I), &
                                p_met%USR_Sigma(I)*p_met%USR_StdScale(1), p_met%USR_Sigma(I)*p_met%USR_StdScale(2), &
                                p_met%USR_Sigma(I)*p_met%USR_StdScale(3), p_met%USR_L(I)      
         ENDDO   
      ELSE
         WRITE (US,"(A)") '  Height   Wind Speed   Horizontal Angle'
         WRITE (US,"(A)") '   (m)        (m/s)          (deg)      '
         WRITE (US,"(A)") '  ------   ----------   ----------------'
     
         DO I=p_met%NumUSRz,1,-1
            WRITE (US,"( 1X,F7.2, 2X,F9.2,2X, 3X,F10.2,6X)") p_met%USR_Z(I), p_met%USR_U(I), p_met%USR_WindDir(I)      
         ENDDO   
      ENDIF
   
   ENDIF

END SUBROUTINE WrSum_UserInput
!=======================================================================
SUBROUTINE WrSum_SpecModel(US, U, HWindDir, p_IEC  )

use TSMods


   INTEGER,                 INTENT(IN) :: US
   TYPE(IEC_ParameterType), INTENT(IN) :: p_IEC                      ! parameters for IEC models   
   REAL(ReKi),              INTENT(IN) :: HWindDir(:)                ! profile of horizontal wind direction
   REAL(ReKi),              INTENT(IN) :: U       (:)                ! profile of steady wind speed

   
   ! local variables:   

   REAL(ReKi)              ::  HalfRotDiam                           ! Half of the rotor diameter
   
   REAL(ReKi)              ::  CVFA                                  ! Cosine of the vertical flow angle
   REAL(ReKi)              ::  SVFA                                  ! Sine of the vertical flow angle
   REAL(ReKi)              ::  CHFA                                  ! Cosine of the horizontal flow angle
   REAL(ReKi)              ::  SHFA                                  ! Sine of the horizontal flow angle
      
   
   REAL(ReKi)              ::  UTmp                                  ! The best fit of observed peak Uh at het height vs jet height
   REAL(ReKi)              ::  U_zb                                  ! The velocity at the bottom of the rotor disk (for estimating log fit)
   REAL(DbKi)              ::  U_zt                                  ! The velocity at the top of the rotor disk (for estimating log fit)
      
   INTEGER                 :: iz, jz                                 ! loop counter/indices of points
   LOGICAL                 ::  HubPr                                 ! Flag to indicate if the hub height is to be printed separately in the summary file

   CHARACTER(200)          :: FormStr                                ! String used to store format specifiers.
   CHARACTER(*),PARAMETER  :: FormStr1 = "('   ',A,' =' ,I9  ,A)"    ! String used to store format specifiers.
   CHARACTER(*),PARAMETER  :: FormStr2 = "('   ',A,' =  ',A)"        ! String used to store format specifiers.
   
   
   ! write to the summary file:
   
   
   WRITE (US,"( // 'Turbulence Simulation Scaling Parameter Summary:' / )")
   WRITE (US,    "('   Turbulence model used                            =  ' , A )")  TRIM(p%met%TMName)

   FormStr  = "('   ',A,' =' ,F9.3,A)"
   
   
   
      ! Write out a parameter summary to the summary file.

IF ( ( p%met%TurbModel_ID  == SpecModel_IECKAI ) .OR. &
     ( p%met%TurbModel_ID  == SpecModel_IECVKM ) .OR. &
     ( p%met%TurbModel_ID  == SpecModel_MODVKM ) .OR. &
     ( p%met%TurbModel_ID  == SpecModel_API    ) )  THEN  ! ADDED BY YGUO on April 192013 snow day!!!
      
      
   IF ( p_IEC%NumTurbInp ) THEN
      WRITE (US,FormStr2)      "Turbulence characteristic                       ", "User-specified"
   ELSE
      WRITE (US,FormStr2)      "Turbulence characteristic                       ", TRIM(p_IEC%IECTurbE)//p_IEC%IECTurbC
      WRITE (US,FormStr2)      "IEC turbulence type                             ", TRIM(p_IEC%IEC_WindDesc)
      
      IF ( p_IEC%IEC_WindType /= IEC_NTM ) THEN       
         WRITE (US,FormStr)    "Reference wind speed average over 10 minutes    ", p_IEC%Vref,                " m/s"
         WRITE (US,FormStr)    "Annual wind speed average at hub height         ", p_IEC%Vave,                " m/s"
      ENDIF
   ENDIF      
   
   WRITE (US,FormStr2)         "IEC standard                                    ", TRIM(p_IEC%IECeditionSTR)
   
   IF ( p%met%TurbModel_ID  /= SpecModel_MODVKM ) THEN
      ! Write out a parameter summary to the summary file.

      WRITE (US,FormStr)       "Mean wind speed at hub height                   ", p%UHub,                      " m/s"

      IF (.NOT. p_IEC%NumTurbInp) THEN ! "A", "B", or "C" turbulence:
         IF ( p_IEC%IECedition == 2 ) THEN
            WRITE (US,FormStr) "Char value of turbulence intensity at 15 m/s    ", 100.0*p_IEC%TurbInt15,     "%"
            WRITE (US,FormStr) "Standard deviation slope                        ", p_IEC%SigmaSlope,          ""
         ELSE                                                                                                 
               ! This is supposed to be the expected value of what is measured at a site.                     
               ! We actually calculate the 90th percentile value to use in the code as the                    
               ! "Characteristic Value".                                                                      
            WRITE (US,FormStr) "Expected value of turbulence intensity at 15 m/s", 100.0*p_IEC%TurbInt15,           "%"
         ENDIF                                                                                                
                                                                                                              
      ENDIF                                                                                                   
                                                                                                              
      WRITE (US,FormStr)       "Characteristic value of standard deviation      ", p_IEC%SigmaIEC(1),          " m/s"
      WRITE (US,FormStr)       "Turbulence scale                                ", p_IEC%Lambda(1),            " m"
                                                                                                              
      IF ( p%met%TurbModel_ID  == SpecModel_IECKAI )  THEN                                                                     
         WRITE (US,FormStr)    "u-component integral scale                      ", p_IEC%IntegralScale(1),     " m"
         WRITE (US,FormStr)    "Coherency scale                                 ", p_IEC%LC,                   " m"
      ELSEIF ( p%met%TurbModel_ID  == SpecModel_IECVKM )  THEN                                                                 
         WRITE (US,FormStr)    "Isotropic integral scale                        ", p_IEC%IntegralScale(1),     " m"
      ENDIF                                                                                                   
                                                                                                              
      WRITE (US,FormStr)       "Characteristic value of hub turbulence intensity", 100.0*p_IEC%TurbInt,              "%"
                                                                                                              
   ELSE   ! ModVKM                                                                                                    
!bjj this is never set in TurbSim:       WRITE (US,FormStr1)      "Boundary layer depth                            ", NINT(h),                   " m"
      WRITE (US,FormStr)       "Site Latitude                                   ", p%met%Latitude,                  " degs"
      WRITE (US,FormStr)       "Hub mean streamwise velocity                    ", p%UHub,                      " m/s"
      WRITE (US,FormStr)       "Hub local u*                                    ", p%met%Ustar,                     " m/s" !BONNIE: is this LOCAL? of Disk-avg
      WRITE (US,FormStr)       "Target IEC Turbulence Intensity                 ", 100.0*p_IEC%TurbInt,             "%"
      WRITE (US,FormStr)       "Target IEC u-component standard deviation       ", p_IEC%SigmaIEC(1),         " m/s"
      WRITE (US,FormStr)       "u-component integral scale                      ", p_IEC%Lambda(1),                 " m"
      WRITE (US,FormStr)       "v-component integral scale                      ", p_IEC%Lambda(2),                 " m"
      WRITE (US,FormStr)       "w-component integral scale                      ", p_IEC%Lambda(3),                 " m"
      WRITE (US,FormStr)       "Isotropic integral scale                        ", p_IEC%LC,                        " m"
   ENDIF                                                                                                      
   WRITE (US,FormStr)          "Gradient Richardson number                      ", 0.0,                       ""

! p%met%Ustar = SigmaIEC/2.15 ! Value based on equating original Kaimal spectrum with IEC formulation

ELSEIF ( p%met%TurbModel_ID == SpecModel_TIDAL ) THEN
   WRITE (US,FormStr2)         "Gradient Richardson number                      ", "N/A"
   WRITE (US,FormStr)          "Mean velocity at hub height                     ", p%UHub,                      " m/s"     
   
ELSE   
 
   WRITE (US,FormStr)          "Gradient Richardson number                      ", p%met%Rich_No,                   ""
   WRITE (US,FormStr)          "Monin-Obukhov (M-O) z/L parameter               ", p%met%ZL,                        ""
                                                                                                              
   IF ( .not. EqualRealNos( p%met%ZL, 0.0_ReKi ) ) THEN                                                                                      
      WRITE (US,FormStr)       "Monin-Obukhov (M-O) length scale                ", p%met%L,                         " m"
   ELSE                                                                                                       
      WRITE (US,FormStr2)      "Monin-Obukhov (M-O) length scale                ", "Infinite"                 
   ENDIF                                                                                                      
   WRITE (US,FormStr)          "Mean wind speed at hub height                   ", p%UHub,                      " m/s"     
    
ENDIF !  TurbModel  == 'IECKAI', 'IECVKM', or 'MODVKM'


HalfRotDiam = 0.5*p%grid%RotorDiameter
U_zt        = getVelocity(p, p%UHub,p%grid%HubHt, p%grid%HubHt+HalfRotDiam)   !Velocity at the top of rotor
U_zb        = getVelocity(p, p%UHub,p%grid%HubHt, p%grid%HubHt-HalfRotDiam)   !Velocity at the bottom of the rotor
      
WRITE(US,'()')   ! A BLANK LINE

SELECT CASE ( TRIM(p%met%WindProfileType) )
   CASE ('JET')
      p%met%PLexp = LOG( U_zt/U_zb ) / LOG( (p%grid%HubHt+HalfRotDiam)/(p%grid%HubHt-HalfRotDiam) )  !TmpReal = p%grid%RotorDiameter/2
      UTmp  = 0.0422*p%met%ZJetMax+10.1979 ! Best fit of observed peak Uh at jet height vs jet height
      
      WRITE (US,FormStr2)      "Wind profile type                               ", "Low-level jet"      
      WRITE (US,FormStr)       "Jet height                                      ",  p%met%ZJetMax,                  " m"
      WRITE (US,FormStr)       "Jet wind speed                                  ",  p%met%UJetMax,                  " m/s"
      WRITE (US,FormStr)       "Upper limit of observed jet wind speed          ",        UTmp,                     " m/s"
      WRITE (US,FormStr)       "Equivalent power law exponent across rotor disk ",  p%met%PLexp,                    ""
      
      IF ( UTmp < p%met%UJetMax ) THEN
         CALL TS_Warn( 'The computed jet wind speed is larger than the ' &
                     //'maximum observed jet wind speed at this height.', -1 )
      ENDIF            
                    
   CASE ('LOG')
      p%met%PLexp = LOG( U_zt/U_zb ) / LOG( (p%grid%HubHt+HalfRotDiam)/(p%grid%HubHt-HalfRotDiam) )  
      
      WRITE (US,FormStr2)      "Wind profile type                               ", "Logarithmic"      
      WRITE (US,FormStr)       "Equivalent power law exponent across rotor disk ",  p%met%PLexp,                    ""

   CASE ('H2L')
      p%met%PLexp = LOG( U_zt/U_zb ) / LOG( (p%grid%HubHt+HalfRotDiam)/(p%grid%HubHt-HalfRotDiam) ) 
      
      WRITE (US,FormStr2)      "Velocity profile type                           ", "Logarithmic (H2L)"      
      WRITE (US,FormStr)       "Equivalent power law exponent across rotor disk ",  p%met%PLexp,                    ""

   CASE ('PL')
      WRITE (US,FormStr2)      "Wind profile type                               ", "Power law"      
      WRITE (US,FormStr)       "Power law exponent                              ",  p%met%PLExp,                    ""
      
   CASE ('USR')
      p%met%PLexp = LOG( U_zt/U_zb ) / LOG( (p%grid%HubHt+HalfRotDiam)/(p%grid%HubHt-HalfRotDiam) )  

      WRITE (US,FormStr2)      "Wind profile type                               ", "Linear interpolation of user-defined profile"
      WRITE (US,FormStr)       "Equivalent power law exponent across rotor disk ",  p%met%PLexp,                    ""
                               
   CASE ('API')
!bjj : fix me:!!! 
      WRITE (US,FormStr2)      "Wind profile type                               ", "API"

      
   CASE DEFAULT                
      WRITE (US,FormStr2)      "Wind profile type                               ", "Power law on rotor disk, logarithmic elsewhere"
      WRITE (US,FormStr)       "Power law exponent                              ",  p%met%PLExp,                    ""
      
END SELECT

WRITE(US,FormStr)              "Mean shear across rotor disk                    ", (U_zt-U_zb)/p%grid%RotorDiameter, " (m/s)/m"
WRITE(US,FormStr)              "Assumed rotor diameter                          ", p%grid%RotorDiameter,             " m"      
WRITE(US,FormStr)              "Surface roughness length                        ", p%met%Z0,                        " m"      
WRITE(US,'()')                                                                                                 ! A BLANK LINE
WRITE(US,FormStr1)             "Number of time steps in the FFT                 ", p%grid%NumSteps,                  ""       
WRITE(US,FormStr1)             "Number of time steps output                     ", p%grid%NumOutSteps,               ""          

IF (p%met%KHtest) THEN
   WRITE(US,"(/'KH Billow Test Parameters:' / )") ! HEADER
   WRITE(US,FormStr)           "Gradient Richardson number                      ", p%met%Rich_No,                   ""
   WRITE(US,FormStr)           "Power law exponent                              ", p%met%PLexp,                     ""
   WRITE(US,FormStr)           "Length of coherent structures                   ", p%grid%UsableTime / 2.0,          " s"
   WRITE(US,FormStr)           "Minimum coherent TKE                            ", 30.0,                      " (m/s)^2"
ENDIF


   ! Write mean flow angles and wind speed profile to the summary file.

WRITE(US,"(//,'Mean Flow Angles:',/)")

FormStr = "(3X,A,F6.1,' degrees')"
WRITE(US,FormStr)  'Vertical   =', p%met%VFlowAng
WRITE(US,FormStr)  'Horizontal =', p%met%HFlowAng


WRITE(US,"(/'Mean Wind Speed Profile:')")

IF ( ALLOCATED( p%met%ZL_profile ) .AND. ALLOCATED( p%met%Ustar_profile ) ) THEN
   WRITE(US,"(/,'   Height    Wind Speed   Horizontal Angle  U-comp (X)   V-comp (Y)   W-comp (Z)   z/L(z)    u*(z)')")
   WRITE(US,"(  '     (m)        (m/s)         (degrees)       (m/s)        (m/s)        (m/s)       (-)      (m/s)')")
   WRITE(US,"(  '   ------    ----------   ----------------  ----------   ----------   ----------   ------   ------')")

   FormStr = '(1X,F8.1,1X,F11.2,5x,F11.2,4x,3(2X,F8.2,3X),2(1X,F8.3))'
ELSE
   WRITE(US,"(/,'   Height    Wind Speed   Horizontal Angle  U-comp (X)   V-comp (Y)   W-comp (Z)')")
   WRITE(US,"(  '     (m)        (m/s)         (degrees)       (m/s)        (m/s)        (m/s)   ')")
   WRITE(US,"(  '   ------    ----------   ----------------  ----------   ----------   ----------')")

   FormStr = '(1X,F8.1,1X,F11.2,5x,F11.2,4x,3(2X,F8.2,3X))'
ENDIF



   ! Get the angles to rotate the wind components from streamwise orientation to the X-Y-Z grid at the Hub
            
CVFA = COS( p%met%VFlowAng*D2R )
SVFA = SIN( p%met%VFlowAng*D2R ) 
CHFA = COS( p%met%HH_HFlowAng*D2R )
SHFA = SIN( p%met%HH_HFlowAng*D2R )

HubPr = .NOT. p%grid%HubOnGrid     !If the hub height is not on the z-grid, print it, too.


   ! Write out the grid points & the hub

DO JZ = p%grid%NumGrid_Z,1, -1
   
   IZ = p%grid%GridPtIndx( (JZ-1)*p%grid%NumGrid_Y+1 )
   
   IF ( HubPr  .AND. ( p%grid%Z(IZ) < p%grid%HubHt ) ) THEN
         
      CHFA = COS( HWindDir(p%grid%HubIndx)*D2R )
      SHFA = SIN( HWindDir(p%grid%HubIndx)*D2R )
         
      IF ( ALLOCATED( p%met%ZL_profile ) ) THEN
            
         WRITE(US,FormStr)  p%grid%Z(p%grid%HubIndx), U(p%grid%HubIndx), HWindDir(p%grid%HubIndx), &
                            U(p%grid%HubIndx)*CHFA*CVFA, U(p%grid%HubIndx)*SHFA*CVFA, U(p%grid%HubIndx)*SVFA, &
                            p%met%ZL_profile(p%grid%HubIndx), p%met%Ustar_profile(p%grid%HubIndx)
      ELSE
         WRITE(US,FormStr)  p%grid%Z(p%grid%HubIndx), U(p%grid%HubIndx), HWindDir(p%grid%HubIndx), &
                            U(p%grid%HubIndx)*CHFA*CVFA, U(p%grid%HubIndx)*SHFA*CVFA, U(p%grid%HubIndx)*SVFA
      ENDIF
   
      HubPr = .FALSE. ! we've printed the hub values, so we don't need to check this anymore
   ENDIF
   
   CHFA = COS( HWindDir(IZ)*D2R )
   SHFA = SIN( HWindDir(IZ)*D2R )

   IF ( ALLOCATED( p%met%ZL_profile ) ) THEN
      WRITE(US,FormStr)  p%grid%Z(IZ), U(IZ), HWindDir(IZ), U(IZ)*CHFA*CVFA, U(IZ)*SHFA*CVFA, U(IZ)*SVFA, &
                           p%met%ZL_profile(IZ), p%met%Ustar_profile(IZ)
   ELSE
      WRITE(US,FormStr)  p%grid%Z(IZ), U(IZ), HWindDir(IZ), U(IZ)*CHFA*CVFA, U(IZ)*SHFA*CVFA, U(IZ)*SVFA
   ENDIF

ENDDO ! JZ
   
   ! Write out the tower points   
DO JZ = 2, SIZE(p%grid%TwrPtIndx)
   IZ =  p%grid%TwrPtIndx(JZ)
   
   CHFA = COS( HWindDir(IZ)*D2R )
   SHFA = SIN( HWindDir(IZ)*D2R )

   IF ( ALLOCATED( p%met%ZL_profile ) ) THEN
      WRITE(US,FormStr)  p%grid%Z(IZ), U(IZ), HWindDir(IZ), U(IZ)*CHFA*CVFA, U(IZ)*SHFA*CVFA, U(IZ)*SVFA, &
                           p%met%ZL_profile(IZ), p%met%Ustar_profile(IZ)
   ELSE
      WRITE(US,FormStr)  p%grid%Z(IZ), U(IZ), HWindDir(IZ), U(IZ)*CHFA*CVFA, U(IZ)*SHFA*CVFA, U(IZ)*SVFA
   ENDIF

ENDDO ! JZ


END SUBROUTINE WrSum_SpecModel
!=======================================================================
SUBROUTINE WrHH_ADtxtfile(V, RootName, TurbInt, ErrStat, ErrMsg)

USE TSMods

   REAL(ReKi),                      INTENT(IN)    :: V     (:,:,:)        !< An array containing the summations of the rows of H (NumSteps,NPoints,3).
   REAL(ReKi),                      INTENT(IN)    :: TurbInt              !< IEC target Turbulence Intensity 
   INTEGER(IntKi),                  intent(  out) :: ErrStat              !< Error level
   CHARACTER(*),                    intent(  out) :: ErrMsg               !< Message describing error
   CHARACTER(*),                    intent(in   ) :: RootName             !< Rootname of output file


   REAL(ReKi)              ::  CVFA                            ! Cosine of the vertical flow angle
   REAL(ReKi)              ::  SVFA                            ! Sine of the vertical flow angle
   REAL(ReKi)              ::  CHFA                            ! Cosine of the horizontal flow angle
   REAL(ReKi)              ::  SHFA                            ! Sine of the horizontal flow angle

   REAL(ReKi)              :: V_Inertial(3)                    ! U,V,W components (inertial)
   REAL(ReKi)              :: UH                               ! horizontal wind speed (U+V components)

   REAL(ReKi)              :: Time                             ! The instantaneous Time (s)
   INTEGER(IntKi)          :: IT                               ! loop counter (time step)
   INTEGER                 :: UAHH                             ! I/O unit for AeroDyn HH data (*.hh  file).

   
   
   CHFA = COS( p%met%HH_HFlowAng*D2R )
   SHFA = SIN( p%met%HH_HFlowAng*D2R )

   CVFA = COS( p%met%VFlowAng*D2R )
   SVFA = SIN( p%met%VFlowAng*D2R )

   CALL GetNewUnit( UAHH, ErrStat, ErrMsg)
   CALL OpenFOutFile ( UAHH, TRIM( RootName)//'.hh', ErrStat, ErrMsg )
   IF (ErrStat >= AbortErrLev) RETURN

   CALL WrScr ( ' Hub-height AeroDyn data were written to "'//TRIM( RootName )//'.hh"' )

   WRITE (UAHH,"( '! This hub-height wind-speed file was generated by ' , A , ' on ' , A , ' at ' , A , '.' )")  TRIM(GetNVD(TurbSim_Ver)), CurDate(), CurTime()
   WRITE (UAHH,"( '!' )")
   WRITE (UAHH,"( '! The requested statistics for this data were:' )")
   WRITE (UAHH,"( '!    Mean Total Wind Speed = ' , F8.3 , ' m/s' )")  p%UHub
IF ( p%met%TurbModel_ID == SpecModel_IECKAI .OR. p%met%TurbModel_ID == SpecModel_IECVKM .OR. p%met%TurbModel_ID == SpecModel_MODVKM ) THEN
   WRITE (UAHH,"( '!    Turbulence Intensity  = ' , F8.3 , '%' )")  100.0*TurbInt
ELSE
   WRITE (UAHH,"( '!' )")
ENDIF
   WRITE (UAHH,"( '!' )")
   WRITE (UAHH,"( '!   Time  HorSpd  WndDir  VerSpd  HorShr  VerShr  LnVShr  GstSpd' )")
   WRITE (UAHH,"( '!  (sec)   (m/s)   (deg)   (m/s)     (-)     (-)     (-)   (m/s)' )")

   DO IT = 1, p%grid%NumOutSteps
      
      Time = p%grid%TimeStep*( IT - 1 )
            
      CALL CalculateWindComponents(V(IT,p%grid%HubIndx,:), p%UHub, p%met%HH_HFlowAng, p%met%VFlowAng, V_Inertial, UH)      
      
      WRITE (UAHH,'(F8.3,3F8.2,3F8.3,F8.2)')  Time, UH, -1.0*R2D*ATAN2( V_Inertial(2) , V_Inertial(1) ), &
                                                    V_Inertial(3), 0.0, p%met%PLExp, 0.0, 0.0 
!bjj: Should we output instantaneous horizontal shear, instead of 0?  
!     Should the power law exponent be an instantaneous value, too?
!                                                   
   END DO

   CLOSE(UAHH)


END SUBROUTINE WrHH_ADtxtfile
!=======================================================================
SUBROUTINE WrHH_binary(V, RootName, ErrStat, ErrMsg)

   ! Output HH binary turbulence parameters for GenPro analysis.
   ! Output order:  Time,U,Uh,Ut,V,W,u',v',w',u'w',u'v',v'w',TKE,CTKE.


USE TSMods

   REAL(ReKi),                      INTENT(IN)    :: V     (:,:,:)                   !< An array containing the summations of the rows of H (NumSteps,NPoints,3).
   CHARACTER(*),                    intent(in   ) :: RootName                        ! Rootname of output file
   INTEGER(IntKi),                  intent(  out) :: ErrStat                         ! Error level
   CHARACTER(*),                    intent(  out) :: ErrMsg                          ! Message describing error


   ! local variables
   
   REAL(ReKi)              :: V_Inertial(3)                    ! U,V,W components (inertial)
   REAL(ReKi)              :: UH                               ! horizontal wind speed (U+V components)
   REAL(ReKi)              :: UT                               ! total wind speed (U+V+W components)
   REAL(ReKi)              :: uv                               ! The instantaneous u'v' Reynolds stress at the hub
   REAL(ReKi)              :: uw                               ! The instantaneous u'w' Reynolds stress at the hub
   REAL(ReKi)              :: vw                               ! The instantaneous v'w' Reynolds stress at the hub
   REAL(ReKi)              :: TKE                              ! The instantaneous TKE at the hub
   REAL(ReKi)              :: CTKE                             ! The instantaneous CTKE the hub
                                                               
   REAL(ReKi)              :: Time                             ! The instantaneous Time (s)
   INTEGER(IntKi)          :: IT                               ! loop counter (time step)
   INTEGER                 :: UGTP                             ! I/O unit for GenPro HH turbulence properties.



   CALL GetNewUnit(UGTP, ErrStat, ErrMsg)
   
   CALL OpenUOutfile ( UGTP , TRIM( RootName)//'.bin', ErrStat, ErrMsg ) 
   IF (ErrStat >= AbortErrLev) RETURN

   CALL WrScr ( ' Hub-height binary turbulence parameters were written to "'//TRIM( RootName )//'.bin"' )

   DO IT = 1, p%grid%NumOutSteps
      
      Time = p%grid%TimeStep*( IT - 1 )
      
      CALL CalculateWindComponents(V(IT,p%grid%HubIndx,:), p%UHub, p%met%HH_HFlowAng, p%met%VFlowAng, V_Inertial, UH, UT)
      CALL CalculateStresses(      V(IT,p%grid%HubIndx,:), uv, uw, vw, TKE, CTKE )
      
      WRITE (UGTP)  REAL(Time,SiKi), REAL(V_Inertial(1),SiKi), REAL(UH,SiKi), REAL(UT,SiKi), &      
                                     REAL(V_Inertial(2),SiKi), REAL(V_Inertial(3),SiKi), &
                                     REAL(V(IT,p%grid%HubIndx,1),SiKi), &
                                     REAL(V(IT,p%grid%HubIndx,2),SiKi), &
                                     REAL(V(IT,p%grid%HubIndx,3),SiKi), &      
                     REAL(uw,SiKi), REAL(uv,SiKi), REAL(vw,SiKi), REAL(TKE,SiKi), REAL(CTKE,SiKi)                        
                                                   
   END DO

   CLOSE(UGTP)
!p%WrFile(FileExt_BIN)

END SUBROUTINE WrHH_binary
!=======================================================================
SUBROUTINE WrHH_text(V, RootName, ErrStat, ErrMsg)

   ! Output HH text turbulence parameters
   ! Output order:  Time,U,Uh,Ut,V,W,u',v',w',u'w',u'v',v'w',TKE,CTKE.


USE TSMods

   REAL(ReKi),                      INTENT(IN)    :: V     (:,:,:)                   !< An array containing the summations of the rows of H (NumSteps,NPoints,3).
   CHARACTER(*),                    intent(in   ) :: RootName                        !< Rootname of output file
   INTEGER(IntKi),                  intent(  out) :: ErrStat                         !< Error level
   CHARACTER(*),                    intent(  out) :: ErrMsg                          !< Message describing error

   ! local variables
   
   REAL(ReKi)              :: V_Inertial(3)                    ! U,V,W components (inertial)
   REAL(ReKi)              :: UH                               ! horizontal wind speed (U+V components)
   REAL(ReKi)              :: UT                               ! total wind speed (U+V+W components)
   REAL(ReKi)              :: uv                               ! The instantaneous u'v' Reynolds stress at the hub
   REAL(ReKi)              :: uw                               ! The instantaneous u'w' Reynolds stress at the hub
   REAL(ReKi)              :: vw                               ! The instantaneous v'w' Reynolds stress at the hub
   REAL(ReKi)              :: TKE                              ! The instantaneous TKE at the hub
   REAL(ReKi)              :: CTKE                             ! The instantaneous CTKE the hub
                                                               
   REAL(ReKi)              :: Time                             ! The instantaneous Time (s)
   INTEGER(IntKi)          :: IT                               ! loop counter (time step)
   INTEGER(IntKi)          :: UFTP                             ! I/O unit for formatted HH turbulence properties

   
   ! p%WrFile(FileExt_DAT)

   CALL GetNewUnit( UFTP, ErrStat, ErrMsg )
   CALL OpenFOutFile ( UFTP, TRIM( RootName)//'.dat', ErrStat, ErrMsg )
   IF (ErrStat >= AbortErrLev) RETURN

   CALL WrScr ( ' Hub-height formatted turbulence parameters were written to "'//TRIM( RootName )//'.dat"' )

   WRITE (UFTP,"( / 'This hub-height turbulence-parameter file was generated by ' , A , ' on ' , A , ' at ' , A , '.' / )") &
                TRIM(GetNVD(TurbSim_Ver)), CurDate(), CurTime()
 
   WRITE (UFTP,"('   Time',6X,'U',7X,'Uh',7X,'Ut',8X,'V',8X,'W',8X,'u''',7X,'v''',7X,'w'''," &
                        //"6X,'u''w''',5X,'u''v''',5X,'v''w''',5X,'TKE',6X,'CTKE')")   
   
   
   DO IT = 1, p%grid%NumOutSteps
      
      Time = p%grid%TimeStep*( IT - 1 )
      
      CALL CalculateWindComponents(V(IT,p%grid%HubIndx,:), p%UHub, p%met%HH_HFlowAng, p%met%VFlowAng, V_Inertial, UH, UT)
      CALL CalculateStresses(      V(IT,p%grid%HubIndx,:), uv, uw, vw, TKE, CTKE )
      
                             
      WRITE(UFTP,'(F7.2,13F9.3)') Time,V_Inertial(1),UH,UT,V_Inertial(2),V_Inertial(3), &
                                    V(IT,p%grid%HubIndx,1), V(IT,p%grid%HubIndx,2), V(IT,p%grid%HubIndx,3), &
                                    uw, uv, vw, TKE, CTKE
      
      
   END DO

   CLOSE(UFTP)

END SUBROUTINE WrHH_text
!=======================================================================
SUBROUTINE WrSum_Stats(V, USig, VSig, WSig, UXBar, UXSig, ErrStat, ErrMsg)


use TSMods

REAL(ReKi),    INTENT(IN   )   :: V           (:,:,:)              ! An array containing the summations of the rows of H (NumSteps,NPoints,3).
REAL(ReKi),    INTENT(  OUT)   ::  USig                            ! Standard deviation of the u-component wind speed at the hub
REAL(ReKi),    INTENT(  OUT)   ::  VSig                            ! Standard deviation of the v-component wind speed at the hub
REAL(ReKi),    INTENT(  OUT)   ::  WSig                            ! Standard deviation of the w-component wind speed at the hub
REAL(ReKi),    INTENT(  OUT)   ::  UXSig                           ! Standard deviation of the U-component wind speed at the hub
REAL(DbKi),    INTENT(  OUT)   ::  UXBar                           ! The mean U-component (u rotated; x-direction) wind speed at the hub

INTEGER(IntKi),intent(  out)   :: ErrStat                          ! Error level
CHARACTER(*),  intent(  out)   :: ErrMsg                           ! Message describing error


REAL(DbKi)                   ::  SumS                            ! Sum of the velocity-squared, used for calculating standard deviations in the summary file
REAL(DbKi)                   ::  UBar                            ! The mean u-component wind speed at the hub
REAL(DbKi)                   ::  UHBar                           ! The mean horizontal wind speed at the hub
REAL(DbKi)                   ::  UHSum2                          ! The sum of the squared horizontal wind speed at the hub
REAL(DbKi)                   ::  UHTmp                           ! The instantaneous horizontal wind speed at the hub
REAL(DbKi)                   ::  UHTmp2                          ! The instantaneous squared horizontal wind speed at the hub
REAL(DbKi)                   ::  USum2                           ! The sum of the squared u-component wind speed at the hub
REAL(DbKi)                   ::  UTBar                           ! The mean total wind speed at the hub
REAL(DbKi)                   ::  UTmp                            ! The instantaneous u-component wind speed at the hub
REAL(DbKi)                   ::  UTmp2                           ! The instantaneous squared u-component wind speed at the hub
REAL(DbKi)                   ::  UTSum2                          ! The sum of the squared total wind speed at the hub
REAL(DbKi)                   ::  UTTmp                           ! The instantaneous total wind speed at the hub
REAL(DbKi)                   ::  UTTmp2                          ! The instantaneous squared total wind speed at the hub
REAL(DbKi)                   ::  UXSum                           ! The sum of the U-component (u rotated) wind speed at the hub
REAL(DbKi)                   ::  UXSum2                          ! The sum of the squared U-component (u rotated) wind speed at the hub
REAL(DbKi)                   ::  UXTmp                           ! The instantaneous U-component (u rotated) wind speed at the hub
REAL(DbKi)                   ::  UXTmp2                          ! The instantaneous squared U-component (u rotated) wind speed at the hub
REAL(DbKi)                   ::  UYBar                           ! The mean V-component (v rotated; y-direction) wind speed at the hub
REAL(DbKi)                   ::  UYSum                           ! The sum of the V-component (v rotated) wind speed at the hub
REAL(DbKi)                   ::  UYSum2                          ! The sum of the squared V-component (v rotated) wind speed at the hub
REAL(DbKi)                   ::  UYTmp                           ! The instantaneous V-component (v rotated) wind speed at the hub
REAL(DbKi)                   ::  UYTmp2                          ! The instantaneous squared V-component (v rotated) wind speed at the hub
REAL(DbKi)                   ::  UZBar                           ! The mean W-component (w rotated; z-direction) wind speed at the hub
REAL(DbKi)                   ::  UZSum                           ! The sum of the W-component (w rotated) wind speed at the hub
REAL(DbKi)                   ::  UZSum2                          ! The sum of the squared W-component (w rotated) wind speed at the hub
REAL(DbKi)                   ::  UZTmp                           ! The instantaneous W-component (w rotated) wind speed at the hub
REAL(DbKi)                   ::  UZTmp2                          ! The instantaneous squared W-component (w rotated) wind speed at the hub
REAL(DbKi)                   ::  VBar                            ! The mean v-component wind speed at the hub
REAL(DbKi)                   ::  VSum2                           ! The sum of the squared v-component wind speed at the hub
REAL(DbKi)                   ::  VTmp                            ! The instantaneous v-component wind speed at the hub
REAL(DbKi)                   ::  VTmp2                           ! The instantaneous squared v-component wind speed at the hub
REAL(DbKi)                   ::  WBar                            ! The mean w-component wind speed at the hub
REAL(DbKi)                   ::  WSum2                           ! The sum of the squared w-component wind speed at the hub
REAL(DbKi)                   ::  WTmp                            ! The instantaneous w-component wind speed at the hub
REAL(DbKi)                   ::  WTmp2                           ! The instantaneous squared w-component wind speed at the hub
                             
REAL(ReKi)                   ::  CHFA                            ! Cosine of the Horizontal Flow Angle
REAL(ReKi)                   ::  CVFA                            ! Cosine of the Vertical Flow Angle
REAL(ReKi)                   ::  SVFA                            ! Sine of the Vertical Flow Angle
REAL(ReKi)                   ::  SHFA                            ! Sine of the Horizontal Flow Angle
REAL(ReKi)                   ::  CTKEmax                         ! Maximum instantaneous Coherent Turbulent Kenetic Energy at the hub
REAL(ReKi)                   ::  TKEmax                          ! Maximum instantaneous Turbulent Kenetic Energy at the hub
                             
                             
REAL(ReKi)                   ::  UHmax                           ! Maximum horizontal wind speed at the hub
REAL(ReKi)                   ::  UHmin                           ! Minimum horizontal wind speed at the hub
REAL(ReKi)                   ::  Umax                            ! Maximum u-component wind speed at the hub
REAL(ReKi)                   ::  Umin                            ! Minimum u-component wind speed at the hub
REAL(ReKi)                   ::  UTSig                           ! Standard deviation of the total wind speed at the hub
REAL(ReKi)                   ::  UT_TI                           ! Turbulent Intensity of the total wind speed at the hub
REAL(ReKi)                   ::  UTmax                           ! Maximum total wind speed at the hub
REAL(ReKi)                   ::  UTmin                           ! Minimum total wind speed at the hub
REAL(ReKi)                   ::  UVMax                           ! Maximum u'v' Reynolds Stress at the hub
REAL(ReKi)                   ::  UVMin                           ! Minimum u'v' Reynolds Stress at the hub
REAL(ReKi)                   ::  UVTmp                           ! The instantaneous u'v' Reynolds stress at the hub
REAL(ReKi)                   ::  UV_RS                           ! The average u'v' Reynolds stress at the hub
REAL(ReKi)                   ::  UVcor                           ! The u-v cross component correlation coefficient at the hub
REAL(ReKi)                   ::  UVsum                           ! The sum of the u'v' Reynolds stress component at the hub
REAL(ReKi)                   ::  UWMax                           ! Maximum u'w' Reynolds Stress at the hub
REAL(ReKi)                   ::  UWMin                           ! Minimum u'w' Reynolds Stress at the hub
REAL(ReKi)                   ::  UWTmp                           ! The instantaneous u'w' Reynolds stress at the hub
REAL(ReKi)                   ::  UW_RS                           ! The average u'w' Reynolds stress at the hub
REAL(ReKi)                   ::  UWcor                           ! The u-w cross component correlation coefficient at the hub
REAL(ReKi)                   ::  UWsum                           ! The sum of the u'w' Reynolds stress component at the hub
REAL(ReKi)                   ::  UXmax                           ! Maximum U-component (X-direction) wind speed at the hub
REAL(ReKi)                   ::  UXmin                           ! Minimum U-component wind speed at the hub
REAL(ReKi)                   ::  UYmax                           ! Maximum V-component (Y-direction) wind speed at the hub
REAL(ReKi)                   ::  UYmin                           ! Minimum V-component wind speed at the hub
REAL(ReKi)                   ::  UYSig                           ! Standard deviation of the V-component wind speed at the hub
REAL(ReKi)                   ::  UZmax                           ! Maximum W-component (Z-direction) wind speed at the hub
REAL(ReKi)                   ::  UZmin                           ! Minimum W-component wind speed at the hub
REAL(ReKi)                   ::  UZSig                           ! Standard deviation of the W-component wind speed at the hub
REAL(ReKi)                   ::  U_TI                            ! The u-component turbulence intensity at the hub
REAL(ReKi)                   ::  Vmax                            ! Maximum v-component wind speed at the hub
REAL(ReKi)                   ::  Vmin                            ! Minimum v-component wind speed at the hub
REAL(ReKi)                   ::  VWMax                           ! Maximum v'w' Reynolds Stress at the hub
REAL(ReKi)                   ::  VWMin                           ! Minimum v'w' Reynolds Stress at the hub
REAL(ReKi)                   ::  VWTmp                           ! The instantaneous v'w' Reynolds stress at the hub
REAL(ReKi)                   ::  VW_RS                           ! The average v'w' Reynolds stress at the hub
REAL(ReKi)                   ::  VWcor                           ! The v-w cross component correlation coefficient at the hub
REAL(ReKi)                   ::  VWsum                           ! The sum of the v'w' Reynolds stress component at the hub
REAL(ReKi)                   ::  V_TI                            ! The v-component turbulence intensity at the hub
REAL(ReKi)                   ::  Wmax                            ! Maximum w-component wind speed at the hub
REAL(ReKi)                   ::  Wmin                            ! Minimum w-component wind speed at the hub
REAL(ReKi)                   ::  W_TI                            ! The w-component turbulence intensity at the hub

REAL(ReKi)                   ::  CTKE                            ! Coherent Turbulent Kenetic Energy at the hub
REAL(ReKi)                   ::  TKE                             ! Turbulent Kenetic Energy at the hub
REAL(ReKi)                   ::  INumSteps                       ! Multiplicative Inverse of the Number of time Steps
REAL(ReKi)                   ::  UHSig                           ! Approximate sigma of the horizontal wind speed at the hub point
REAL(ReKi)                   ::  SUstar                          ! Simulated U-star at the hub
REAL(ReKi), ALLOCATABLE      ::  SDary      (:)                  ! The array of standard deviations (NumGrid_Z,NumGrid_Y).

REAL(ReKi)                   ::  UH_TI                           ! TI of the horizontal wind speed at the hub point


INTEGER(IntKi)               :: IT, IVec, IY, IZ, II

   CHARACTER(200)            :: FormStr                                  ! String used to store format specifiers.
   CHARACTER(*),PARAMETER    :: FormStr2 = "(6X,A,' component: ',F8.3,' m/s')"        ! String used to store format specifiers.


   ! Initialize statistical quantities for hub-height turbulence parameters.

   CALL WrScr ( ' Computing hub-height statistics' )

   INumSteps   = 1.0/p%grid%NumSteps

   CTKEmax = -HUGE( CTKEmax )
   TKEmax  = -HUGE(  TKEmax )
   UBar    =      0.0
   UHBar   =      0.0
   UHmax   = -HUGE( UHmax )
   UHmin   =  HUGE( UHmin )
   UHSum2  =      0.0
   Umax    = -HUGE( Umax )
   Umin    =  HUGE( Umin )
   USum2   =      0.0
   UTBar   =      0.0
   UTmax   = -HUGE( UTmax )
   UTmin   =  HUGE( UTmin )
   UTSum2  =      0.0
   UV_RS   =      0.0
   UVMax   =  V(1,p%grid%HubIndx,1)*V(1,p%grid%HubIndx,2) 
   UVMin   =  HUGE( UVMin )
   UVsum   =      0.0
   UW_RS   =      0.0
   UWMax   =  V(1,p%grid%HubIndx,1)*V(1,p%grid%HubIndx,3) 
   UWMin   =  HUGE( UWMin )
   UWsum   =      0.0
   VBar    =      0.0
   Vmax    = -HUGE( Vmax )
   Vmin    =  HUGE( Vmin )
   VSum2   =      0.0
   VW_RS   =      0.0
   VWMax   =  V(1,p%grid%HubIndx,2)*V(1,p%grid%HubIndx,3)
   VWMin   =  HUGE( VWMin )
   VWsum   =      0.0
   WBar    =      0.0
   Wmax    = -HUGE( Wmax )
   Wmin    =  HUGE( Wmin )
   WSum2   =      0.0
   UXBar   =      0.0
   UXmax   = -HUGE( UXmax )
   UXmin   =  HUGE( UXmin )
   UXSum   =      0.0
   UXSum2  =      0.0
   UXTmp   =      0.0
   UXTmp2  =      0.0
   UYBar   =      0.0
   UYmax   = -HUGE( UYmax )
   UYmin   =  HUGE( UYmin )
   UYSum   =      0.0
   UYSum2  =      0.0
   UYTmp   =      0.0
   UYTmp2  =      0.0
   UZBar   =      0.0
   UZmax   = -HUGE( UZmax )
   UZmin   =  HUGE( UZmin )
   UZSum   =      0.0
   UZSum2  =      0.0
   UZTmp   =      0.0
   UZTmp2  =      0.0

   CHFA = COS( p%met%HH_HFlowAng*D2R )
   SHFA = SIN( p%met%HH_HFlowAng*D2R )

   CVFA = COS( p%met%VFlowAng*D2R )
   SVFA = SIN( p%met%VFlowAng*D2R )

   DO IT=1,p%grid%NumSteps

         ! Calculate longitudinal (UTmp), lateral (VTmp), and upward (WTmp)
         ! values for hub station, as well as rotated (UXTmp, UYTmp, UZTmp) 
         ! components applying specified flow angles.

         ! Add mean wind speed to the streamwise component
      UTmp = V(IT,p%grid%HubIndx,1) + p%UHub
      VTmp = V(IT,p%grid%HubIndx,2)
      WTmp = V(IT,p%grid%HubIndx,3)
   
         ! Rotate the wind components from streamwise orientation to the X-Y-Z grid at the Hub      
      UXTmp = UTmp*CHFA*CVFA - VTmp*SHFA - WTmp*CHFA*SVFA
      UYTmp = UTmp*SHFA*CVFA + VTmp*CHFA - WTmp*SHFA*SVFA  
      UZTmp = UTmp*SVFA                  + WTmp*CVFA

         ! Calculate hub horizontal wind speed (UHTmp) and Total wind speed (UTTmp)
      UTmp2 = UTmp*UTmp          !flow coordinates
      VTmp2 = VTmp*VTmp
      WTmp2 = WTmp*WTmp
   
      UXTmp2 = UXTmp*UXTmp       !inertial frame coordinates
      UYTmp2 = UYTmp*UYTmp
      UZTmp2 = UZTmp*UZTmp
   
      UHTmp2 = UXTmp2 + UYTmp2   !inertial frame coordinates
      UTTmp2 = UHTmp2 + UZTmp2

      UHTmp = SQRT( UHTmp2 )     !inertial frame coordinates
      UTTmp = SQRT( UTTmp2 )

         ! Form running sums for hub standard deviations

      UBar   = UBar   + UTmp     !flow coordinates
      VBar   = VBar   + VTmp     !flow coordinates
      WBar   = WBar   + WTmp     !flow coordinates

      USum2  = USum2  + UTmp2    !flow coordinates
      VSum2  = VSum2  + VTmp2    !flow coordinates
      WSum2  = WSum2  + WTmp2    !flow coordinates

      UXBar   = UXBar + UXTmp
      UYBar   = UYBar + UYTmp
      UZBar   = UZBar + UZTmp

      UXSum2  = UXSum2 + UXTmp2
      UYSum2  = UYSum2 + UYTmp2
      UZSum2  = UZSum2 + UZTmp2

      UHBar  = UHBar  + UHTmp
      UTBar  = UTBar  + UTTmp

      UHSum2 = UHSum2 + UHTmp2
      UTSum2 = UTSum2 + UTTmp2


         ! Determine hub extremes.
      
      IF ( UTmp  > Umax  )  Umax  = UTmp     !flow coordinates,
      IF ( UTmp  < Umin  )  Umin  = UTmp     !flow coordinates,

      IF ( VTmp  > Vmax  )  Vmax  = VTmp     !flow coordinates,
      IF ( VTmp  < Vmin  )  Vmin  = VTmp     !flow coordinates,

      IF ( WTmp  > Wmax  )  Wmax  = WTmp     !flow coordinates,
      IF ( WTmp  < Wmin  )  Wmin  = WTmp     !flow coordinates,

      IF ( UXTmp > UXmax )  UXmax = UXTmp
      IF ( UXTmp < UXmin )  UXmin = UXTmp

      IF ( UYTmp > UYmax )  UYmax = UYTmp
      IF ( UYTmp < UYmin )  UYmin = UYTmp

      IF ( UZTmp > UZmax )  UZmax = UZTmp
      IF ( UZTmp < UZmin )  UZmin = UZTmp

      IF ( UHTmp > UHmax )  UHmax = UHTmp
      IF ( UHTmp < UHmin )  UHmin = UHTmp

      IF ( UTTmp > UTmax )  UTmax = UTTmp
      IF ( UTTmp < UTmin )  UTmin = UTTmp

         ! Find maxes and mins of instantaneous hub Reynolds stresses u'w', u'v', and v'w'

      UVTmp = V(IT,p%grid%HubIndx,1)*V(IT,p%grid%HubIndx,2)
      UWTmp = V(IT,p%grid%HubIndx,1)*V(IT,p%grid%HubIndx,3)
      VWTmp = V(IT,p%grid%HubIndx,2)*V(IT,p%grid%HubIndx,3)

      IF     ( UVTmp < UVMin )  THEN
         UVMin = UVTmp
      ELSEIF ( UVTmp > UVMax )  THEN
         UVMax = UVTmp
      ENDIF

      IF     ( UWTmp < UWMin )  THEN
         UWMin = UWTmp
      ELSEIF ( UWTmp > UWMax )  THEN
         UWMax = UWTmp
      ENDIF

      IF     ( VWTmp < VWMin )  THEN
         VWMin = VWTmp
      ELSEIF ( VWTmp > VWMax )  THEN
         VWMax = VWTmp
      ENDIF

         ! Find maximum of instantaneous TKE and CTKE.

      TKE  = 0.5*(V(IT,p%grid%HubIndx,1)*V(IT,p%grid%HubIndx,1) + V(IT,p%grid%HubIndx,2)*V(IT,p%grid%HubIndx,2) + V(IT,p%grid%HubIndx,3)*V(IT,p%grid%HubIndx,3))
      CTKE = 0.5*SQRT(UVTmp*UVTmp + UWTmp*UWTmp + VWTmp*VWTmp)

      IF (CTKE > CTKEmax) CTKEmax = CTKE
      IF ( TKE >  TKEmax)  TKEmax =  TKE

         ! Find sums for mean and square Reynolds stresses for hub-level simulation.
      UVsum = UVsum + UVTmp
      UWsum = UWsum + UWTmp
      VWsum = VWsum + VWTmp

   ENDDO ! IT



      ! Calculate mean hub-height Reynolds stresses.
   UW_RS =  UWsum*INumSteps
   UV_RS =  UVsum*INumSteps
   VW_RS =  VWsum*INumSteps

      ! Simulated Hub UStar.
   SUstar = SQRT( ABS( UW_RS ) )

      ! Calculate mean values for hub station.

   UBar = UBar*INumSteps
   VBar = VBar*INumSteps
   WBar = WBar*INumSteps

   UXBar = UXBar*INumSteps
   UYBar = UYBar*INumSteps
   UZBar = UZBar*INumSteps

   UHBar = UHBar*INumSteps
   UTBar = UTBar*INumSteps


      ! Calculate the standard deviations for hub station.
      ! (SNWind/SNLwind-3D) NOTE: This algorithm is the approximate algorithm.
      ! bjj: do the algebra and you'll find that it's std() using the 1/n definition   

   USig  = SQRT( MAX( USum2 *INumSteps-UBar *UBar , 0.0_DbKi ) )
   VSig  = SQRT( MAX( VSum2 *INumSteps-VBar *VBar , 0.0_DbKi ) )
   WSig  = SQRT( MAX( WSum2 *INumSteps-WBar *WBar , 0.0_DbKi ) )

   UXSig = SQRT( MAX( UXSum2*INumSteps-UXBar*UXBar, 0.0_DbKi ) )
   UYSig = SQRT( MAX( UYSum2*INumSteps-UYBar*UYBar, 0.0_DbKi ) )
   UZSig = SQRT( MAX( UZSum2*INumSteps-UZBar*UZBar, 0.0_DbKi ) )

   UHSig = SQRT( MAX( UHSum2*INumSteps-UHBar*UHBar, 0.0_DbKi ) )
   UTSig = SQRT( MAX( UTSum2*INumSteps-UTBar*UTBar, 0.0_DbKi ) )


      ! Calculate Cross-component correlation coefficients
   UWcor = ( UW_RS ) / (USig * WSig) ! this definition assumes u' and w' have zero mean
   UVcor = ( UV_RS ) / (USig * VSig)
   VWcor = ( VW_RS ) / (VSig * WSig)


      ! Calculate turbulence intensities.
   U_TI = USig/UBar    
   V_TI = VSig/UBar
   W_TI = WSig/UBar

   UH_TI = UHSig/UHBar
   UT_TI = UTSig/UTBar


      ! Write out the hub-level stats to the summary file.

   CALL WrScr ( ' Writing statistics to summary file' )

   WRITE(US,"(//,'Hub-Height Simulated Turbulence Statistical Summary:')")
   WRITE(US,"(/,3X,'Type of Wind        Min (m/s)   Mean (m/s)    Max (m/s)  Sigma (m/s)       TI (%)')")
   WRITE(US,"(  3X,'----------------    ---------   ----------    ---------  -----------       ------')")

   FormStr = "(3X,A,F13.2,2F13.2,2F13.3)"
   !bjj for analysis, extra precision: FormStr = "(3X,A,F13.2,2F13.2,2F13.6)"

   WRITE (US,FormStr)  'Longitudinal (u)',  Umin,  UBar,  Umax,  USig, 100.0* U_TI
   WRITE (US,FormStr)  'Lateral (v)     ',  Vmin,  VBar,  Vmax,  VSig, 100.0* V_TI
   WRITE (US,FormStr)  'Vertical (w)    ',  Wmin,  WBar,  Wmax,  WSig, 100.0* W_TI
   WRITE (US,FormStr)  'U component     ', UXmin, UXBar, UXmax, UXSig, 100.0*UXSig/UXBar
   WRITE (US,FormStr)  'V component     ', UYmin, UYBar, UYmax, UYSig, 100.0*UYSig/UXBar
   WRITE (US,FormStr)  'W component     ', UZmin, UZBar, UZmax, UZSig, 100.0*UZSig/UXBar
   WRITE (US,FormStr)  'Horizontal (U&V)', UHmin, UHBar, UHmax, UHSig, 100.0*UH_TI
   WRITE (US,FormStr)  'Total           ', UTmin, UTBar, UTmax, UTSig, 100.0*UT_TI

   WRITE(US,"(/,3X,'                    Min Reynolds     Mean Reynolds    Max Reynolds    Correlation')")
   WRITE(US,"(  3X,'Product             Stress (m/s)^2   Stress (m/s)^2   Stress (m/s)^2  Coefficient')")
   WRITE(US,"(  3X,'----------------    --------------   --------------   --------------  -----------')")

   FormStr = "(3X,A,3(3X,F12.3,3X),F11.3)"
   WRITE (US,FormStr)  "u'w'            ",  UWMin,  UW_RS, UWMax, UWcor
   WRITE (US,FormStr)  "u'v'            ",  UVMin,  UV_RS, UVMax, UVcor
   WRITE (US,FormStr)  "v'w'            ",  VWMin,  VW_RS, VWMax, VWcor

   FormStr = "(3X,A,' = ',F10.3,A)"
   WRITE(US,"(/)")   ! blank line
   WRITE(US,FormStr)  "Friction Velocity (Ustar) ", SUstar,  " m/s"
   WRITE(US,FormStr)  "Maximum Instantaneous TKE ", TKEmax,  " (m/s)^2"
   WRITE(US,FormStr)  "Maximum Instantaneous CTKE", CTKEmax, " (m/s)^2"

      !  Allocate the array of standard deviations.

   CALL AllocAry( SDary, p%grid%NumGrid_Y, 'SDary (standard deviations)', ErrStat, ErrMsg)
   IF (ErrStat >= AbortErrLev) RETURN


      ! Calculate standard deviations for each grid point.  Write them to summary file.

   WRITE(US,"(//,'Grid Point Variance Summary:',/)")
   WRITE(US,"(3X,'Y-coord',"//TRIM(Num2LStr(p%grid%NumGrid_Y))//"F8.2)")  p%grid%Y(1:p%grid%NumGrid_Y)


   UTmp = 0
   VTmp = 0
   WTmp = 0

   DO IVec=1,3

      WRITE(US,"(/,3X'Height   Standard deviation at grid points for the ',A,' component:')")  Comp(IVec)

      DO IZ=p%grid%NumGrid_Z,1,-1

         DO IY=1,p%grid%NumGrid_Y

            II   = (IZ-1)*p%grid%NumGrid_Y+IY
            II   = p%grid%GridPtIndx(II)
            
            SumS = 0.0

            DO IT=1,p%grid%NumSteps
               SumS = SumS + V(IT,II,IVec)**2            
            ENDDO ! IT         

            SDary(IY) = SQRT(SumS*INumSteps)   !  Was:  SDary(IZ,IY) = SQRT(SumS*INumSteps)/U(IZ,NumGrid/2)

         ENDDO ! IY

         WRITE(US,"(F9.2,1X,"//TRIM(Num2LStr(p%grid%NumGrid_Y))//"F8.3)") p%grid%Z( p%grid%GridPtIndx((IZ-1)*p%grid%NumGrid_Y+1) ), SDary(1:p%grid%NumGrid_Y)

         IF ( IVec == 1 ) THEN
            UTmp = UTmp + SUM( SDary )
         ELSEIF ( IVec == 2 ) THEN
            VTmp = VTmp + SUM( SDary )
         ELSE
            WTmp = WTmp + SUM( SDary )
         ENDIF
      ENDDO ! IZ

   ENDDO ! Ivec

   
   WRITE(US,"(/'   Mean standard deviation across all grid points:')")
   WRITE(US,FormStr2) Comp(1), UTmp / ( p%grid%NumGrid_Y*p%grid%NumGrid_Z )
   WRITE(US,FormStr2) Comp(2), VTmp / ( p%grid%NumGrid_Y*p%grid%NumGrid_Z ) 
   WRITE(US,FormStr2) Comp(3), WTmp / ( p%grid%NumGrid_Y*p%grid%NumGrid_Z ) 


      !  Deallocate the array of standard deviations.

   IF ( ALLOCATED( SDary ) )  DEALLOCATE( SDary )



END SUBROUTINE WrSum_Stats
!=======================================================================
!> Calculate the mean velocity and turbulence intensity of the U-component
!! of the interpolated hub point for comparison with InflowWind output.
SUBROUTINE WrSum_InterpolatedHubStats(p, V, US)

      ! passed variables:
   TYPE(TurbSim_ParameterType),     INTENT(IN)     ::  p                               !< TurbSim's parameters
   REAL(ReKi),                      INTENT(INOUT)  ::  V(:,:,:)                        !< velocity, aligned along the streamwise direction without mean values added 
   INTEGER(IntKi)                 , INTENT(IN)     ::  US                              !< unit number of file in which to print a summary of the scaling used. If < 1, will not print summary.

      ! local variables:
REAL(DbKi)              ::  CGridSum                        ! The sums of the velocity components at the points surrounding the hub (or at the hub if it's on the grid)
REAL(DbKi)              ::  CGridSum2                       ! The sums of the squared velocity components at the points surrouding the hub 

REAL(ReKi)              ::  TmpV                            ! Temporarily holds the value of the v component
REAL(ReKi)              ::  TmpY                            ! Temp variable for interpolated hub point
REAL(ReKi)              ::  TmpZ                            ! Temp variable for interpolated hub point
REAL(ReKi)              ::  Tmp_YL_Z                        ! Temp variable for interpolated hub point
REAL(ReKi)              ::  Tmp_YH_Z                        ! Temp variable for interpolated hub point

REAL(ReKi)              ::  UGridMean                       ! Average wind speed at the points surrounding the hub
REAL(ReKi)              ::  UGridSig                        ! Standard deviation of the wind speed at the points surrounding the hub
REAL(ReKi)              ::  UGridTI                         ! Turbulent Intensity of the points surrounding the hub

INTEGER                 ::  ZHi_YHi                         ! Index for interpolation of hub point, if necessary
INTEGER                 ::  ZHi_YLo                         ! Index for interpolation of hub point, if necessary
INTEGER                 ::  ZLo_YHi                         ! Index for interpolation of hub point, if necessary
INTEGER                 ::  ZLo_YLo                         ! Index for interpolation of hub point, if necessary
INTEGER                 ::  IT                              ! Index for time step

INTEGER                 ::  IZ_Lo, IY_Lo                    ! Index for lower bound of box surrounding hub point

            
   
      ! Calculate mean value & turb intensity of U-component of the interpolated hub point (for comparison w/ AeroDyn output)

   ! Note that this uses the InflowWind interpolation scheme, which may be updated some day so that it doesn't
   ! depend on which dimension we interpolate first.
      
   IY_Lo = INT(   0.5_ReKi * p%grid%GridWidth     / p%grid%GridRes_Y ) + 1
   IZ_Lo = INT( ( p%grid%HubHt - p%grid%Zbottom ) / p%grid%GridRes_Z ) + 1     
   
   
      ! Get points for bi-linear interpolation  ( indx @ (iy,iz) is (iz-1)*numgrid_y + iy, assuming a full grid (needs to be modified for user-defined spectra)
   ZLo_YLo   = p%grid%GridPtIndx( ( IZ_Lo - 1 )*p%grid%NumGrid_Y + IY_Lo     )
   ZHi_YLo   = p%grid%GridPtIndx( ( IZ_Lo     )*p%grid%NumGrid_Y + IY_Lo     )
   ZLo_YHi   = p%grid%GridPtIndx( ( IZ_Lo - 1 )*p%grid%NumGrid_Y + IY_Lo + 1 )
   ZHi_YHi   = p%grid%GridPtIndx( ( IZ_Lo     )*p%grid%NumGrid_Y + IY_Lo + 1 )
    
   TmpZ      = (p%grid%HubHt - p%grid%Z(p%grid%GridPtIndx((IZ_Lo-1)*p%grid%NumGrid_Y + 1)))/p%grid%GridRes_Z
   TmpY      = ( 0.0_ReKi    - p%grid%Y(p%grid%GridPtIndx( IY_Lo                        )))/p%grid%GridRes_Y
   CGridSum  = 0.0
   CGridSum2 = 0.0

   DO IT=1,p%grid%NumSteps
         
   ! Interpolate within the grid for this time step.

      Tmp_YL_Z  = ( V( IT, ZHi_YLo, 1 ) - V( IT, ZLo_YLo, 1 ) )*TmpZ + V( IT, ZLo_YLo, 1 )
      Tmp_YH_Z  = ( V( IT, ZHi_YHi, 1 ) - V( IT, ZLo_YHi, 1 ) )*TmpZ + V( IT, ZLo_YHi, 1 )
      TmpV      = ( Tmp_YH_Z - Tmp_YL_Z )*TmpY + Tmp_YL_Z

      CGridSum  = CGridSum  + TmpV
      CGridSum2 = CGridSum2 + TmpV*TmpV
   ENDDO ! IT

   UGridMean = CGridSum/p%grid%NumSteps
   UGridSig  = SQRT( ABS( (CGridSum2/p%grid%NumSteps) - UGridMean*UGridMean ) )
   UGridTI   = 100.0*UGridSig/UGridMean


      ! Put the average statistics of the four center points in the summary file.

   WRITE (US,"(//,'U-component (X) statistics from the interpolated hub point:',/)")
   WRITE(US,"(3X,A,' =',F9.4,A)")  'Mean' , UGridMean, ' m/s'
   WRITE(US,"(3X,A,' =',F9.4,A)")  'TI  ' , UGridTI  , ' %'


END SUBROUTINE WrSum_InterpolatedHubStats
!=======================================================================
SUBROUTINE WrSum_EchoInputs()

use TSMods

   INTEGER                      :: I          ! loop counter
   CHARACTER(10)                :: TmpStr     ! temporary string used to write output to summary file



!..................................................................................................................................
   WRITE (US,"( / 'Runtime Options:' / )")
   WRITE (US,"( I10 , 2X , 'Random seed #1' )"                            )  p%RNG%RandSeed(1)
   
   IF (p%RNG%pRNG == pRNG_INTRINSIC) THEN
      WRITE (US,"( I10 , 2X , 'Random seed #2' )"                         )  p%RNG%RandSeed(2)
   ELSE
      WRITE (US,"( 4X, A6, 2X, 'Type of random number generator' )"       )  p%RNG%RNG_type
   ENDIF
      
   WRITE (US,"( L10 , 2X , 'Output binary HH turbulence parameters?' )"   )  p%WrFile(FileExt_BIN)
   WRITE (US,"( L10 , 2X , 'Output formatted turbulence parameters?' )"   )  p%WrFile(FileExt_DAT)
   WRITE (US,"( L10 , 2X , 'Output AeroDyn HH files?' )"                  )  p%WrFile(FileExt_HH)
   WRITE (US,"( L10 , 2X , 'Output AeroDyn FF files?' )"                  )  p%WrFile(FileExt_BTS)
   WRITE (US,"( L10 , 2X , 'Output BLADED FF files?' )"                   )  p%WrFile(FileExt_WND)
   WRITE (US,"( L10 , 2X , 'Output tower data?' )"                        )  p%WrFile(FileExt_TWR)
   WRITE (US,"( L10 , 2X , 'Output formatted FF files?' )"                )  p%WrFile(FileExt_UVW)
   WRITE (US,"( L10 , 2X , 'Output coherent turbulence time step file?' )")  p%WrFile(FileExt_CTS)
   WRITE (US,"( L10 , 2X , 'Clockwise rotation when looking downwind?' )" )  p%grid%Clockwise   
   
   SELECT CASE ( p%IEC%ScaleIEC )
      CASE (0)
         TmpStr= "NONE"
      CASE (1, -1)   ! included the -1 for reading t/f on other systems
         TmpStr = "HUB"
      CASE (2)
         TmpStr = "ALL"
   ENDSELECT   
   
   WRITE (US,"( I2, ' - ', A5, 2X , 'IEC turbulence models scaled to exact specified standard deviation' )")  p%IEC%ScaleIEC, TRIM(TmpStr)
   
   
!..................................................................................................................................
   WRITE (US,"( // 'Turbine/Model Specifications:' / )")
   WRITE (US,"( I10   , 2X , 'Vertical grid-point matrix dimension' )"  )  p%grid%NumGrid_Z
   WRITE (US,"( I10   , 2X , 'Horizontal grid-point matrix dimension' )")  p%grid%NumGrid_Y
   WRITE (US,"( F10.3 , 2X , 'Time step [seconds]' )"                   )  p%grid%TimeStep
   WRITE (US,"( F10.3 , 2X , 'Analysis time [seconds]' )"               )  p%grid%AnalysisTime
   WRITE (US,"( F10.3 , 2X , 'Usable output time [seconds]' )"          )  p%grid%UsableTime
   WRITE (US,"( F10.3 , 2X , 'Hub height [m]' )"                        )  p%grid%HubHt  
   WRITE (US,"( F10.3 , 2X , 'Grid height [m]' )"                       )  p%grid%GridHeight
   WRITE (US,"( F10.3 , 2X , 'Grid width [m]' )"                        )  p%grid%GridWidth
   WRITE (US,"( F10.3 , 2X , 'Vertical flow angle [degrees]' )"         )  p%met%VFlowAng
   WRITE (US,"( F10.3 , 2X , 'Horizontal flow angle [degrees]' )"       )  p%met%HFlowAng
   

!..................................................................................................................................
   WRITE (US,"( // 'Meteorological Boundary Conditions:' / )")
   WRITE (US, "( 4X , A6 , 2X , '"//TRIM( p%met%TMName )//" spectral model' )")  p%met%TurbModel
   IF (p%IEC%IECstandard > 0) then
      WRITE (US,"( 7X, I3, 2X, 'IEC standard: ', A )")  p%IEC%IECstandard, TRIM(p%IEC%IECeditionSTR)
      IF (p%IEC%NumTurbInp) THEN
         WRITE (US,"( F10.3 , 2X , 'Percent turbulence intensity, ', A )")  p%IEC%PerTurbInt, TRIM(p%IEC%IECeditionSTR)
      ELSE
         WRITE (US,"( 9X , A1 , 2X , 'IEC turbulence characteristic' )"  )  p%IEC%IECTurbC   
      END IF
      
      SELECT CASE ( p%IEC%IEC_WindType )
         CASE (IEC_NTM)
            TmpStr= "NTM"
         CASE (IEC_ETM)   
            TmpStr = "ETM"
         CASE (IEC_EWM1)
            TmpStr = "EWM1"
         CASE (IEC_EWM50)
            TmpStr = "EWM50"
         CASE (IEC_EWM100)
            TmpStr = "EWM100"
      ENDSELECT
                        
      WRITE (US,"( 4X, A6 , 2X , 'IEC ', A )")  TRIM(p%IEC%IECTurbE)//TRIM(TmpStr), TRIM(p%IEC%IEC_WindDesc)
      
   ELSE
      WRITE (US,"( 7X, A3, 2X, 'IEC standard' )"                        )  'N/A'
      IF (p%met%KHtest) THEN
         WRITE (US,"( 4X, A6, 2X, 'Kelvin-Helmholtz billow test case' )")  'KHTEST'   
      ELSE   
         WRITE (US,"( A10, 2X, 'IEC turbulence characteristic' )"       )  'N/A'   
      END IF      
      WRITE (US,"( A10 , 2X , 'IEC turbulence type' )"                  )  'N/A'
      
   END IF

   IF ( p%IEC%IEC_WindType == IEC_ETM ) THEN
      WRITE (US,"( F10.3, 2X, 'IEC Extreme Turbulence Model (ETM) ""c"" parameter [m/s]' )") p%IEC%ETMc 
   ELSE
      WRITE (US,"(  A10, 2X, 'IEC Extreme Turbulence Model (ETM) ""c"" parameter [m/s]' )")  'N/A'   
   END IF      

   WRITE (US,"(   A10 , 2X , 'Wind profile type' )"                  )  p%met%WindProfileType   
   WRITE (US,"( F10.3 , 2X , 'Reference height [m]' )"               )  p%met%RefHt  !BJJ: TODO: check if refht makes sense (or is used) for USR profile.
   IF ( p%met%WindProfileType == 'USR' ) THEN
      WRITE (US,"(  A10, 2X, 'Reference wind speed [m/s]' )"         )  'N/A'   
   ELSE
      WRITE (US,"( F10.3 , 2X , 'Reference wind speed [m/s]' )"      )  p%met%URef
   END IF
   
      
   IF ( p%met%WindProfileType == 'JET' ) THEN
      WRITE (US,"( F10.3, 2X, 'Jet height [m]' )"                    ) p%met%ZJetMax 
   ELSE
      WRITE (US,"(  A10, 2X, 'Jet height [m]' )"                     )  'N/A'   
   END IF      
      
   IF ( INDEX( 'JLUHA', p%met%WindProfileType(1:1) ) > 0   ) THEN
      WRITE (US,"(  A10, 2X, 'Power law exponent' )"                 )  'N/A'   
   ELSE
      WRITE (US,"( F10.3 , 2X , 'Power law exponent' )"              )  p%met%PLExp
   END IF      
      
   IF ( p%met%TurbModel_ID==SpecModel_TIDAL  ) THEN
      WRITE (US,"(  A10, 2X, 'Surface roughness length [m]' )"       )  'N/A'   
   ELSE      
      WRITE (US,"( F10.3 , 2X , 'Surface roughness length [m]' )"    )  p%met%Z0
   END IF
        
!..................................................................................................................................
!***
IF ( p%met%TurbModel_ID == SpecModel_IECKAI .OR. p%met%TurbModel_ID == SpecModel_IECVKM .OR. p%met%TurbModel_ID == SpecModel_API ) RETURN  
!***

WRITE (US,"( // 'Non-IEC Meteorological Boundary Conditions:' / )")
   
   IF ( p%met%TurbModel_ID /= SpecModel_IECKAI .AND. p%met%TurbModel_ID /= SpecModel_IECVKM .AND. p%met%TurbModel_ID /= SpecModel_API ) THEN
      WRITE (US,"( F10.3 , 2X , 'Site latitude [degrees]' )"    )  p%met%Latitude      
   ELSE
      WRITE (US,"( A10 , 2X , 'Site latitude [degrees]' )"    )  'N/A'            
   END IF

!***
IF ( p%met%TurbModel_ID == SpecModel_ModVKM ) RETURN  
!***

   IF ( .NOT. p%met%IsIECModel .AND. & 
        p%met%TurbModel_ID /= SpecModel_TIDAL   ) THEN
      WRITE (US,"( F10.3 , 2X , 'Gradient Richardson number' )")  p%met%Rich_No
   ELSE
      WRITE (US,"(   a10 , 2X , 'Gradient Richardson number' )")  'N/A'
   END IF
   
   IF ( .NOT. p%met%IsIECModel  ) THEN
      WRITE (US,"( F10.3 , 2X , 'Friction or shear velocity [m/s]' )")  p%met%Ustar
      
      IF (p%met%ZL>=0. .AND. p%met%TurbModel_ID /= SpecModel_GP_LLJ) THEN
         WRITE (US,'(   A10 , 2X , "Mixing layer depth [m]" )'       )  'N/A'
      ELSE         
         WRITE (US,"( F10.3 , 2X , 'Mixing layer depth [m]' )"       )  p%met%ZI
      END IF

   ELSE
      WRITE (US,'(   A10 , 2X , "Friction or shear velocity [m/s]" )')  'N/A'
      WRITE (US,'(   A10 , 2X , "Mixing layer depth [m]" )'          )  'N/A'
   END IF

   IF (.NOT. p%met%UWskip) THEN
      WRITE (US,'( F10.3 , 2X , "Mean hub u''w'' Reynolds stress" )' )  p%met%PC_UW
   ELSE
      WRITE (US,'(   A10 , 2X , "Mean hub u''w'' Reynolds stress" )' )  'N/A'
   END IF
   
   IF (.NOT. p%met%UVskip) THEN
      WRITE (US,'( F10.3 , 2X , "Mean hub u''v'' Reynolds stress" )' )  p%met%PC_UV
   ELSE
      WRITE (US,'(   A10 , 2X , "Mean hub u''v'' Reynolds stress" )' )  'N/A'
   END IF

   IF (.NOT. p%met%VWskip) THEN
      WRITE (US,'( F10.3 , 2X , "Mean hub v''w'' Reynolds stress" )' )  p%met%PC_VW
   ELSE   
      WRITE (US,'(   A10 , 2X , "Mean hub v''w'' Reynolds stress" )' )  'N/A'
   END IF
   
   do i=1,3
      WRITE (US,"( '(',F9.3,',',G10.3,')',2X , A,'-component coherence parameters' )")  p%met%InCDec(i), p%met%InCohB(i), Comp(i)
   end do
   
   WRITE (US,'( F10.3 , 2X , "Coherence exponent" )' )  p%met%CohExp
   
!..................................................................................................................................
!***
IF ( .NOT. p%WrFile(FileExt_CTS) ) RETURN  
!***
      WRITE (US,"( // 'Coherent Turbulence Scaling Parameters:' / )")
   

      IF ( LEN( TRIM(p_CohStr%CTEventPath) ) <= 10 )  THEN
         WRITE (US,"( A10 , 2X , 'Name of the path containing the coherent turbulence data files' )") TRIM(p_CohStr%CTEventPath)
      ELSE
         WRITE (US,"( A, /, 12X , 'Name of the path containing the coherent turbulence data files' )") TRIM(p_CohStr%CTEventPath)
      ENDIF
      WRITE (US,"( 7X, A3, 2X, 'Type of coherent turbulence data files' )") TRIM(p_CohStr%CText)
!      WRITE (US,"( L10 , 2X , 'Randomize the disturbance scale and location?' )")  Randomize
      WRITE (US,"( F10.3 , 2X , 'Disturbance scale (ratio of wave height to rotor diameter)' )")  p_CohStr%DistScl
      WRITE (US,"( F10.3 , 2X , 'Fractional location of tower centerline from right' )")  p_CohStr%CTLy
      WRITE (US,"( F10.3 , 2X , 'Fractional location of hub height from the bottom of the dataset' )")  p_CohStr%CTLz
      WRITE (US,"( F10.3 , 2X , 'Minimum start time for coherent structures [seconds]' )")  p_CohStr%CTStartTime
            
!..................................................................................................................................
   
END SUBROUTINE WrSum_EchoInputs
!=======================================================================
!> This subroutine processes the user input from the IECstandard line
!! and splits it into the standard and edition being used
SUBROUTINE ProcessLine_IECstandard(Line, IsIECModel, TurbModel_ID, IECstandard, IECedition, IECeditionStr, ErrStat, ErrMsg )

   CHARACTER(*),   INTENT(INOUT) :: Line                !< on entry, the line from the input file. may be modified in this routine
   LOGICAL,        INTENT(IN   ) :: IsIECModel          !< Flag to indicate if this is an IEC model
   INTEGER(IntKi), INTENT(IN   ) :: TurbModel_ID        !< Turbulence model identifier 
   INTEGER(IntKi), INTENT(  OUT) :: IECedition          !< IEC edition
   INTEGER(IntKi), INTENT(  OUT) :: IECstandard         !< IEC standard
   CHARACTER(*),   INTENT(  OUT) :: IECeditionStr       !< string describing the IEC standard/edition being used
   
   INTEGER(IntKi), intent(  out) :: ErrStat             !< Error level
   CHARACTER(*),   intent(  out) :: ErrMsg              !< Message describing error
   
   INTEGER(IntKi)                :: IOS                 ! local error code
   INTEGER(IntKi)                :: TmpIndex            ! index into string
   CHARACTER(*), PARAMETER       :: IECstandardErrMsg = 'The IECstandard input parameter must be either "1", "2", or "3"' &
                                   // ' with an optional IEC 61400-1 edition number ("1-ED2"). If specified, the edition number must be "2" or "3".'
   
   CHARACTER( 23), PARAMETER     :: IECeditionStr_p (3) = &   ! strings for the 
                                   (/'IEC 61400-1 Ed. 1: 1993', &
                                     'IEC 61400-1 Ed. 2: 1999', &
                                     'IEC 61400-1 Ed. 3: 2005'/)            ! The string description of the IEC 61400-1 standard being used

   
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   IF ( .NOT. IsIECModel ) THEN  !bjj: SpecModel==SpecModel_MODVKM is not in the IEC standard
      IECstandard   = 0
      IECedition    = 0
      IECeditionStr = ""
      RETURN
   ENDIF ! IEC   
   
   
      ! Did the line contain "T" or "F", which could be interpreted by Fortran as a number?
   CALL Conv2UC( LINE )
   IF ( (Line(1:1) == 'T') .OR.  (Line(1:1) == 'F') ) THEN
      CALL SetErrStat(ErrID_Fatal, IECstandardErrMsg, ErrStat, ErrMsg, 'ProcessLine_IECstandard')
      RETURN
   ENDIF

   
      ! Did the user enter an edition number?
   TmpIndex = INDEX(Line, "-ED")
   IF ( TmpIndex > 0 ) THEN
      READ ( Line(TmpIndex+3:),*,IOSTAT=IOS ) IECedition

      IF (IOS /= 0) THEN
         CALL SetErrStat(ErrID_Fatal, IECstandardErrMsg, ErrStat, ErrMsg, 'ProcessLine_IECstandard')
         RETURN
      END IF
         
      IF ( IECedition < 1 .OR. IECedition > 3 ) THEN
         CALL SetErrStat(ErrID_Fatal, IECstandardErrMsg, ErrStat, ErrMsg, 'ProcessLine_IECstandard')
         RETURN
      ENDIF

      Line = Line(1:TmpIndex-1)
   ELSE
      IECedition = 0
   ENDIF

      ! What standard did the user enter?   
   READ ( Line,*,IOSTAT=IOS ) IECstandard

   SELECT CASE ( IECstandard )
         
   CASE ( 1 ) ! use the IEC 64100-1 standard, either edition 2 or 3
      IF (IECedition < 1 ) THEN  ! Set up the default
         IF ( TurbModel_ID == SpecModel_IECVKM .OR. TurbModel_ID == SpecModel_USRVKM ) THEN
            IECedition = 2   ! The von Karman model is not specified in edition 3 of the -1 standard
         ELSE
            IECedition = 3
         ENDIF
      ELSE
         IF ( IECedition < 2 ) THEN
            CALL SetErrStat(ErrID_Fatal, IECstandardErrMsg, ErrStat, ErrMsg, 'ProcessLine_IECstandard')
            RETURN                  
         ENDIF
      ENDIF
      IECeditionSTR = IECeditionStr_p(IECedition)

   CASE ( 2 ) ! use the IEC 64100-2 (small turbine) standard, which the same as 64100-1, Ed. 2 with "A" or user-specified turbulence
      IF (IECedition < 1 ) THEN  ! Set up the default
         IECedition = 2    ! This is the edition of the -1 standard
      ELSE
         CALL SetErrStat(ErrID_Fatal, IECstandardErrMsg//' The edition number cannot be specified for the 61400-2 standard.', ErrStat, ErrMsg, 'ProcessLine_IECstandard')
         RETURN                  
      ENDIF
      IECeditionSTR = 'IEC 61400-2 Ed. 2: 2005'

   CASE ( 3 ) ! Use the IEC 61400-3 (Offshore) standard, which is the same as 61400-1 except it has a different power law exponent
      IF (IECedition < 1 ) THEN  ! Set up the default

         IF ( TurbModel_ID /= SpecModel_IECKAI ) THEN
            CALL SetErrStat(ErrID_Fatal, ' The Kaimal model (IECKAI) is the only turbulence model valid for the 61400-3 standard.', ErrStat, ErrMsg, 'ProcessLine_IECstandard')
            RETURN                  
         ENDIF
         IECedition = 3    ! This is the edition of the -1 standard

      ELSE
         CALL SetErrStat(ErrID_Fatal, IECstandardErrMsg//' The edition number cannot be specified for the 61400-3 standard.', ErrStat, ErrMsg, 'ProcessLine_IECstandard')
         RETURN                  
      ENDIF
      IECeditionSTR = 'IEC 61400-3 Ed. 1: 2006'

   CASE DEFAULT
      CALL SetErrStat(ErrID_Fatal, IECstandardErrMsg, ErrStat, ErrMsg, 'ProcessLine_IECstandard')
      RETURN                  

   END SELECT
END SUBROUTINE ProcessLine_IECstandard
!=======================================================================
!> This subroutine processes the user input from the IECturbc line
!! and splits it into the NumTurbInp, IECTurbC, and KHtest variables.
SUBROUTINE ProcessLine_IECturbc(Line, IsIECModel, IECstandard, IECedition, IECeditionStr, NumTurbInp, IECTurbC, PerTurbInt, KHtest, ErrStat, ErrMsg)
   REAL(ReKi)  ,   INTENT(  OUT) :: PerTurbInt          !< Percent Turbulence Intensity
   
   INTEGER(IntKi), INTENT(IN   ) :: IECedition          !< IEC edition
   INTEGER(IntKi), INTENT(IN   ) :: IECstandard         !< IEC standard
   
   LOGICAL,        INTENT(IN   ) :: IsIECModel          !< Flag to indicate if this is an IEC model
   LOGICAL,        INTENT(  OUT) :: NumTurbInp          !< Flag to indicate if turbulence is user-specified (as opposed to IEC standard A, B, or C)
   LOGICAL,        INTENT(  OUT) :: KHtest              !< Flag to indicate that turbulence should be extreme, to demonstrate effect of KH billows
   
   CHARACTER(*),   INTENT(IN   ) :: IECeditionStr       !< string describing the IEC standard/edition being used
   CHARACTER(1),   INTENT(  OUT) :: IECTurbC            !< IEC turbulence characteristic
   CHARACTER(*),   INTENT(INOUT) :: Line                !< on entry, the line from the input file. may be modified in this routine
   
   !INTEGER(IntKi), INTENT(IN   ) :: TurbModel_ID        !< Turbulence model identifier 
   
   INTEGER(IntKi), intent(  out) :: ErrStat             !< Error level
   CHARACTER(*),   intent(  out) :: ErrMsg              !< Message describing error
   
   INTEGER(IntKi)                :: IOS                 ! local error code
   CHARACTER(*), PARAMETER       :: IECstandardErrMsg = 'The IECstandard input parameter must be either "1", "2", or "3"' &
                                   // ' with an optional IEC 61400-1 edition number ("1-ED2"). If specified, the edition number must be "2" or "3".'

   ErrStat = ErrID_None
   ErrMsg  = ""
   
   IF ( IsIECModel ) THEN
      KHtest = .FALSE.
      
      READ (Line,*,IOSTAT=IOS) IECTurbC

      CALL Conv2UC( IECTurbC )

      IF ( (IECTurbC == 'T') .OR.  (IECTurbC == 'F') ) THEN
         CALL SetErrStat( ErrID_Fatal, 'The IEC turbulence characteristic must be either "A", "B", "C", or a real number.', ErrStat, ErrMsg, 'ProcessLine_IECturbc')
         RETURN
      ENDIF

      ! Check to see if the entry was a number.

      READ (Line,*,IOSTAT=IOS)  PerTurbInt

      IF ( IOS == 0 )  THEN

         ! Let's use turbulence value.

         NumTurbInp = .TRUE.         
         IECTurbC = ""

      ELSE

         ! Let's use one of the standard turbulence values (A or B or C).

         NumTurbInp = .FALSE.
         PerTurbInt = 0.0           ! will be set later if necessary
         
         IECTurbC   = ADJUSTL( Line )
         CALL Conv2UC(  IECTurbC )
         SELECT CASE ( IECTurbC )
         CASE ( 'A' )
         CASE ( 'B' )
            IF ( IECstandard == 2 ) THEN
               CALL SetErrStat( ErrID_Fatal, 'The IEC 61400-2 turbulence characteristic must be either "A" or a real number.', ErrStat, ErrMsg, 'ProcessLine_IECturbc')
            ENDIF
         CASE ( 'C' )
            IF ( IECstandard == 2 ) THEN
               CALL SetErrStat( ErrID_Fatal, 'The IEC 61400-2 turbulence characteristic must be either "A" or a real number.', ErrStat, ErrMsg, 'ProcessLine_IECturbc')
            ELSEIF ( IECedition < 3 ) THEN
               CALL SetErrStat( ErrID_Fatal, 'The turbulence characteristic for '//TRIM(IECeditionSTR )// &
                                 ' must be either "A", "B", or a real number.', ErrStat, ErrMsg, 'ProcessLine_IECturbc')
            ENDIF
         CASE DEFAULT
            CALL SetErrStat( ErrID_Fatal, 'The IEC turbulence characteristic must be either "A", "B", "C", or a real number.', ErrStat, ErrMsg, 'ProcessLine_IECturbc')
         END SELECT  ! IECTurbC

      ENDIF
      
      

   ELSE  

      Line = ADJUSTL( Line )
      CALL Conv2UC( Line )

      KHtest = ( Line(1:6) == 'KHTEST' )

         ! These variables are not used for non-IEC turbulence

      NumTurbInp = .FALSE.
      PerTurbInt = 0.0
      IECTurbC   = ""
      
   ENDIF   
   
   
END SUBROUTINE ProcessLine_IECturbc

!=======================================================================
END MODULE TS_FileIO
