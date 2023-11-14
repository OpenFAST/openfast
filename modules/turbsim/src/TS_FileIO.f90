!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2014, 2016  National Renewable Energy Laboratory
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
   use TS_RandNum
   
   IMPLICIT                NONE

CONTAINS

!=======================================================================
!> This subroutine reads parameters from the primary TurbSim input file.
!> It validates most of the meteorology data (because it is used to 
!> calculate default values later in the routine)
SUBROUTINE ReadInputFile(InFile, p, OtherSt_RandNum, ErrStat, ErrMsg)


  ! This subroutine is used to read parameters from the input file.

   IMPLICIT             NONE

   CHARACTER(*),                 INTENT(IN)    :: InFile             !< name of the primary TurbSim input file
   TYPE(TurbSim_ParameterType),  INTENT(INOUT) ::  p                 !< TurbSim's parameters
   TYPE(RandNum_OtherStateType), INTENT(INOUT) :: OtherSt_RandNum    !< other states for random numbers (next seed, etc)


   INTEGER(IntKi)  ,             INTENT(OUT)   :: ErrStat            !< allocation status
   CHARACTER(*) ,                INTENT(OUT)   :: ErrMsg             !< error message



      ! Local variables

   REAL(ReKi)                    :: InCVar     (2)                   ! Contains the coherence parameters (used for input)
   REAL(ReKi)                    :: tmp                              ! variable for estimating Ustar and calculating wind speeds
   REAL(ReKi)                    :: TmpUary (3)                      !Temporary vector to store windSpeed(z) values
   REAL(ReKi)                    :: TmpUstar(3)                      !Temporary vector to store ustar(z) values
   REAL(ReKi)                    :: TmpUstarD                        !Temporary ustarD value
   REAL(ReKi)                    :: RotorDiskHeights (3)             !Temporary vector to store height(z) values
   REAL(ReKi)                    :: TmpZLary(3)                      !Temporary vector to store zL(z) values
                              
   INTEGER                       :: TmpIndex                         ! Contains the index number when searching for substrings
   INTEGER                       :: UI                               ! I/O unit for input file
   INTEGER                       :: UnEc                             ! I/O unit for echo file
                              
   LOGICAL                       :: getDefaultPLExp                  ! Whether a default PLExp needs to be calculated
   LOGICAL                       :: getDefaultURef                   ! Whether a default URef needs to be calculated 
   LOGICAL                       :: getDefaultZJetMax                ! Whether a default ZJetMax needs to be calculated 
   
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

   CHARACTER(1024)               :: UserFile                         ! Name of file containing user-defined spectra or time-series files
   CHARACTER(1024)               :: ProfileFile                      ! Name of the file containing profile data for user-defined velocity profiles and/or USRVKM model
   CHARACTER(*), PARAMETER       :: RoutineName = 'ReadInputFile'
   
      ! Initialize some variables:
   ErrStat = ErrID_None
   ErrMsg  = ""
      
   p%met%NumUSRz = 0  ! initialize the number of points in a user-defined wind profile
   
   
   UnEc = -1
   Echo = .FALSE.   
   CALL GetPath( InFile, PriPath )     ! Input files will be relative to the path where the primary input file is located.
        
   
   !===============================================================================================================================
   ! Open input file
   !===============================================================================================================================

   CALL GetNewUnit( UI, ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      
   CALL OpenFInpFile( UI, InFile, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
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
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      
      CALL ReadCom( UI, InFile, "File Heading Line 2",    ErrStat2, ErrMsg2, UnEc)
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      
      CALL ReadCom( UI, InFile, "Runtime Options Heading",ErrStat2, ErrMsg2, UnEc)
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
                  
      CALL ReadVar( UI, InFile, Echo, 'Echo', 'Echo switch', ErrStat2, ErrMsg2, UnEc )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         IF (ErrStat >= AbortErrLev) THEN
            CALL Cleanup()
            RETURN
         END IF   
         
         
      IF (.NOT. Echo .OR. I > 1) EXIT !exit this loop
   
         ! Otherwise, open the echo file, then rewind the input file and echo everything we've read
      
      I = I + 1         ! make sure we do this only once (increment counter that says how many times we've read this file)
      
      
      CALL OpenEcho ( UnEc, TRIM(p%RootName)//'.ech', ErrStat2, ErrMsg2, TurbSim_Ver )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF   
      
      IF ( UnEc > 0 )  WRITE (UnEc,'(/,A,/)')  'Data from '//TRIM(TurbSim_Ver%Name)//' primary input file "'//TRIM( InFile )//'":'
      
      REWIND( UI, IOSTAT=ErrStat2 )  
         IF (ErrStat2 /= 0_IntKi ) THEN
            CALL SetErrStat( ErrID_Fatal, 'Error rewinding file "'//TRIM(InFile)//'".', ErrStat, ErrMsg, RoutineName)    
            RETURN
         END IF         
                  
   END DO
   
         
   
      ! RandSeed(1)
   CALL ReadVar( UI, InFile, p%RNG%RandSeed(1), "RandSeed(1)", "Random seed #1",ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

      ! RandSeed(2)
   CALL ReadVar( UI, InFile, Line, "RandSeed(2)", "Random seed #2",ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF

         !  Check if alternate random number generator is to be used >>>>>>>>>>>>>>>>

      READ (Line,*,IOSTAT=ErrStat2) Line1  ! check the first character to make sure we don't have T/F, which can be interpreted as 1/-1 or 0 in Fortran

      CALL Conv2UC( Line1 )
      IF ( (Line1 == 'T') .OR.  (Line1 == 'F') ) THEN
         CALL SetErrStat( ErrID_Fatal, ' RandSeed(2): Invalid RNG type.', ErrStat, ErrMsg, RoutineName)
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
            CALL SetErrStat( ErrID_Fatal, ' RandSeed(2): Invalid alternative random number generator.', ErrStat, ErrMsg, RoutineName)
            CALL Cleanup()
            RETURN
         ENDIF

      ENDIF
      
      !<<<<<<<<<<<<<<<<<<<<<< end rng


      ! --------- Read the flag for writing the binary HH (GenPro) turbulence parameters. -------------
   CALL ReadVar( UI, InFile, p%WrFile(FileExt_BIN), "WrBHHTP", "Output binary HH turbulence parameters? [RootName.bin]",ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

      ! --------- Read the flag for writing the formatted turbulence parameters. ----------------------
   CALL ReadVar( UI, InFile, p%WrFile(FileExt_DAT), "WrFHHTP", "Output formatted turbulence parameters? [RootName.dat]",ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

      ! ---------- Read the flag for writing the AeroDyn HH files. -------------------------------------
   CALL ReadVar( UI, InFile, p%WrFile(FileExt_HH), "WrADHH", "Output AeroDyn HH files? [RootName.hh]",ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

      ! ---------- Read the flag for writing the AeroDyn FF files. ---------------------------------------
   CALL ReadVar( UI, InFile, p%WrFile(FileExt_BTS), "WrADFF", "Output AeroDyn FF files? [RootName.bts]",ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

      ! ---------- Read the flag for writing the BLADED FF files. -----------------------------------------
   CALL ReadVar( UI, InFile, p%WrFile(FileExt_WND) , "WrBLFF", "Output BLADED FF files? [RootName.wnd]",ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

      ! ---------- Read the flag for writing the AeroDyn tower files. --------------------------------------
   CALL ReadVar( UI, InFile, p%WrFile(FileExt_TWR), "WrADTWR", "Output tower data? [RootName.twr]",ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

      ! ---------- Read the flag for writing the HAWC FF files. ---------------------------------------
   CALL ReadVar( UI, InFile, p%WrFile(FileExt_HAWC), "WrHAWCFF", "Output HAWC FF files? [RootName-u.bin, -v.bin, -w.bin, .hawc]",ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

      ! ---------- Read the flag for writing the formatted FF files. ---------------------------------------
   CALL ReadVar( UI, InFile, p%WrFile(FileExt_UVW), "WrFMTFF", "Output formatted FF files? [RootName.u, .v, .w]",ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

      ! ---------- Read the flag for writing coherent time series files. --------------------------------------
   CALL ReadVar( UI, InFile, p%WrFile(FileExt_CTS), "WrACT", "Output coherent time series files? [RootName.cts]",ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

      ! ---------- Read the flag for determining IEC scaling -----------------------------------------------------
   CALL ReadVar( UI, InFile, p%IEC%ScaleIEC, "ScaleIEC", "Scale IEC turbulence models to specified standard deviation?",&
                                 ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

      
      ! we'll check the errors before going to the next section of the input file
   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF

   !===============================================================================================================================
   ! Read the turbine/model specifications.
   !===============================================================================================================================

   CALL ReadCom( UI, InFile, "Turbine/Model Specifications Heading Line 1",ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      
   CALL ReadCom( UI, InFile, "Turbine/Model Specifications Heading Line 2",ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

      ! ------------ Read in the vertical matrix dimension. ---------------------------------------------
   CALL ReadVar( UI, InFile, p%grid%NumGrid_Z, "NumGrid_Z", "Vertical grid-point matrix dimension [-]",ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

      ! ------------ Read in the lateral matrix dimension. ---------------------------------------------
   CALL ReadVar( UI, InFile, p%grid%NumGrid_Y, "NumGrid_Y", "Horizontal grid-point matrix dimension [-]",ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

      ! ------------ Read in the time step. ---------------------------------------------
   CALL ReadVar( UI, InFile, p%grid%TimeStep, "TimeStep", "Time step [seconds]",ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

      ! ------------ Read in the analysis time. ---------------------------------------------
   CALL ReadVar( UI, InFile, p%grid%AnalysisTime, "AnalysisTime", "Analysis time [seconds]",ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

      ! ------------ Read in the usable time. ---------------------------------------------
   CALL ReadVar( UI, InFile, Line, "UsableTime", "Usable output time [seconds]",ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! Check if usable time is "ALL" (for periodic files) >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
         READ( Line, *, IOSTAT=ErrStat2) p%grid%UsableTime

         IF ( ErrStat2 /= 0 ) THEN ! Line didn't contain a number
            CALL Conv2UC( Line )
            IF ( TRIM(Line) == 'ALL' ) THEN
               p%grid%Periodic   = .TRUE.
               p%grid%UsableTime = p%grid%AnalysisTime
            ELSE
               CALL SetErrStat( ErrID_Fatal, 'The usable output time must be a number greater than zero (or the string "ALL").', &
                                ErrStat, ErrMsg, RoutineName )
               CALL Cleanup()
               RETURN
            END IF
         ELSE
            p%grid%Periodic = .FALSE.

            CALL CheckRealVar( p%grid%UsableTime, 'UsableTime', ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         END IF
      ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< end check for UsableTime = "ALL" (periodic)

      ! ------------ Read in the hub height. ---------------------------------------------
   CALL ReadVar( UI, InFile, p%grid%HubHt, "HubHt", "Hub height [m]",ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

      ! ------------ Read in the grid height. ---------------------------------------------
   CALL ReadVar( UI, InFile, p%grid%GridHeight, "GridHeight", "Grid height [m]",ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

      ! ------------ Read in the grid width. ---------------------------------------------
   CALL ReadVar( UI, InFile, p%grid%GridWidth, "GridWidth", "Grid width [m]",ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

      ! ------------ Read in the vertical flow angle. ---------------------------------------------
   CALL ReadVar( UI, InFile, p%met%VFlowAng, "VFlowAng", "Vertical flow angle [degrees]",ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

      ! ------------ Read in the horizontal flow angle. ---------------------------------------------
   CALL ReadVar( UI, InFile, p%met%HFlowAng, "HFlowAng", "Horizontal flow angle [degrees]",ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)


!..................................................................................................................................
!  Do some error checking on the runtime options and turbine/model specifications before we read the meteorology data
!..................................................................................................................................

   IF ( p%IEC%ScaleIEC > 2 .OR. p%IEC%ScaleIEC < 0 ) CALL SetErrStat( ErrID_Fatal, 'The value for parameter ScaleIEC must be 0, 1, or 2.',    ErrStat, ErrMsg, RoutineName)
   IF ( p%grid%NumGrid_Z < 2 )                       CALL SetErrStat( ErrID_Fatal, 'The matrix must be >= 2x2.',                              ErrStat, ErrMsg, RoutineName)
   IF ( p%grid%NumGrid_Y < 2 )                       CALL SetErrStat( ErrID_Fatal, 'The matrix must be >= 2x2.',                              ErrStat, ErrMsg, RoutineName) 
   IF ( 0.5*p%grid%GridHeight > p%grid%HubHt  )      CALL SetErrStat( ErrID_Fatal, 'The hub must be higher than half of the grid height.',    ErrStat, ErrMsg, RoutineName)
   IF ( p%grid%GridWidth <=  0.0_ReKi )              CALL SetErrStat( ErrID_Fatal, 'The grid width must be greater than zero.',               ErrStat, ErrMsg, RoutineName)
   IF ( p%grid%HubHt <=  0.0 )                       CALL SetErrStat( ErrID_Fatal, 'The hub height must be greater than zero.',               ErrStat, ErrMsg, RoutineName)
   IF ( p%grid%AnalysisTime <=  0.0 )                CALL SetErrStat( ErrID_Fatal, 'The analysis time must be greater than zero.',            ErrStat, ErrMsg, RoutineName)
   IF ( p%grid%TimeStep <=  0.0 )                    CALL SetErrStat( ErrID_Fatal, 'The time step must be greater than zero.',                ErrStat, ErrMsg, RoutineName)
   IF ( ABS( p%met%VFlowAng ) > 45.0 )               CALL SetErrStat( ErrID_Fatal, 'The vertical flow angle must not exceed +/- 45 degrees.', ErrStat, ErrMsg, RoutineName)
   IF ( p%grid%UsableTime <=  0.0 )                  CALL SetErrStat( ErrID_Fatal, 'The usable output time must be a number greater than zero'&
                                                                                                                   //' or the string "ALL".', ErrStat, ErrMsg, RoutineName)
      
!..................................................................................................................................
!  initialize secondary parameters that will be used to calculate default values in the meteorological boundary conditions section
!..................................................................................................................................
      ! Initialize the RNG (for computing "default" values that contain random variates)
   CALL RandNum_Init(p%RNG, OtherSt_RandNum, ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            
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
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      
   CALL ReadCom( UI, InFile, "Meteorological Boundary Conditions Heading Line 2",ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

      ! ------------ Read in the turbulence model. ---------------------------------------------
   CALL ReadVar( UI, InFile, p%met%TurbModel, "TurbModel", "spectral model",ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      
      ! ------------ Read in the UserFile------------------- ---------------------------------------------
   CALL ReadVar( UI, InFile, UserFile, "UserFile", "Name of the input file for user-defined spectra or time-series inputs",ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
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
            p%met%TMName = 'Uniform user-input'
            p%met%TurbModel_ID = SpecModel_USER
            
            CALL GetUSRspec(UserFile, p, UnEc, ErrStat2, ErrMsg2)
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            
         CASE ( 'TIMESR' )
            p%met%TMName = 'User-input time series'
            p%met%TurbModel_ID = SpecModel_TimeSer 
            
            CALL GetUSRTimeSeries(UserFile, p, UnEc, ErrStat2, ErrMsg2)
               CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat,ErrMsg, RoutineName)
               IF (ErrStat >= AbortErrLev) THEN
                  CALL Cleanup()
                  RETURN
               END IF
               
            CALL TimeSeriesToSpectra( p, ErrStat2, ErrMsg2 )
               CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat,ErrMsg, RoutineName)

         CASE DEFAULT
   !BONNIE: todo: add the UsrVKM model to this list when the model is complete 
            CALL SetErrStat( ErrID_Fatal, 'The turbulence model must be one of the following: "IECKAI", "IECVKM", "SMOOTH",' &
                       //' "WF_UPW", "WF_07D", "WF_14D", "NWTCUP", "GP_LLJ", "TIDAL", "RIVER", "API", "USRINP", "TIMESR" "NONE".', ErrStat, ErrMsg, RoutineName)
            CALL Cleanup()
            RETURN

         END SELECT  ! TurbModel
      
         ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< end TurbModel verification
         
!bjj: todo: verify that the API model sets the parameters for IECKAI as well (because it's using IECKAI for the v and w components)

      ! ------------ Read in the IEC standard and edition numbers. ---------------------------------------------
   CALL ReadVar( UI, InFile, Line, "IECstandard", "Number of the IEC standard",ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! Process this line for IEC standard & edition & IECeditionStr >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      CALL ProcessLine_IECstandard( Line, p%met%IsIECModel, p%met%TurbModel_ID, p%IEC%IECstandard, p%IEC%IECedition, p%IEC%IECeditionStr, ErrStat2, ErrMsg2)
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         IF (ErrStat >= AbortErrLev) THEN
            CALL Cleanup()
            RETURN
         END IF         
      ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< end processing of IECstandard input variable

      ! ------------ Read in the IEC turbulence characteristic. ---------------------------------------------
   CALL ReadVar( UI, InFile, Line, "IECturbc", "IEC turbulence characteristic",ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

      ! Process this line for NumTurbInp, IECPerTurbInt, IECTurbC, and KHtest >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      CALL ProcessLine_IECTurbc(Line, p%met%IsIECModel, p%IEC%IECstandard, p%IEC%IECedition, p%IEC%IECeditionStr, &
                                p%IEC%NumTurbInp, p%IEC%IECTurbC, p%IEC%PerTurbInt, p%met%KHtest, ErrStat2, ErrMsg2)
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         IF (ErrStat >= AbortErrLev) THEN
            CALL Cleanup()
            RETURN
         END IF
      ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< end processing of IECturbc input variable
   
      ! ------------ Read in the IEC wind turbulence type ---------------------------------------------
   CALL ReadVar( UI, InFile, Line, "IEC_WindType", "IEC turbulence type",ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   
         ! Process this line for IECTurbE, Vref, IEC_WindType, and IEC_WindDes >>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      CALL ProcessLine_IEC_WindType(Line, p, ErrStat2, ErrMsg2)
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            IF (ErrStat >= AbortErrLev) THEN
               CALL Cleanup()
               RETURN
            END IF
         ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< end processing of IEC_WindType input variable

     
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>         
   ! set default ETMc, WindProfileType, Z0, CohExp, and Latitude 
   ! for use in ReadRVarDefault, ReadRAryDefault, and ReadCVarDefault routines
CALL DefaultMetBndryCndtns(p)     ! Requires turbModel (some require RICH_NO, which we'll have to redo later)    
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
         

   ! ------------ Read in the ETM c parameter (IEC 61400-1, Ed 3: Section 6.3.2.3, Eq. 19) ----------------------
   CALL ReadRVarDefault( UI, InFile, p%IEC%ETMc, "ETMc", 'IEC Extreme Turbulence Model (ETM) "c" parameter [m/s]', UnEc, &
                         UseDefault, ErrStat2, ErrMsg2, IGNORE=(p%IEC%IEC_WindType /= IEC_ETM ))
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      
   ! ------------ Read in the wind profile type -----------------------------------------------------------------
   CALL ReadCVarDefault( UI, InFile, p%met%WindProfileType, "WindProfileType", "Wind profile type", UnEc, UseDefault, ErrStat2, ErrMsg2) !converts WindProfileType to upper case
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

      ! ------------ Read in the ProfileFile------------------- ---------------------------------------------
   CALL ReadVar( UI, InFile, ProfileFile, "ProfileFile", 'Name of the input file for profiles used with WindProfileType="USR" or TurbModel="USRVKM"',ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   IF ( PathIsRelative( ProfileFile ) ) ProfileFile = TRIM(PriPath)//TRIM(ProfileFile)
            
   ! ------------ Read in the height for the reference wind speed. ---------------------------------------------
   CALL ReadVar( UI, InFile, p%met%RefHt, "RefHt", "Reference height [m]",ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   ! ------------ Read in the reference wind speed. -----------------------------------------------------
   IsUnusedParameter = p%IEC%IEC_WindType > IEC_ETM  .OR. INDEX('TU',p%met%WindProfileType(1:1)) > 0 ! p%IEC%IEC_WindType > IEC_ETM == EWM models   
   CALL ReadRVarDefault( UI, InFile, p%met%URef, "URef", "Reference wind speed [m/s]", UnEc, getDefaultURef, ErrStat2, ErrMsg2, &
                     IGNORE=IsUnusedParameter ) ! p%IEC%IEC_WindType > IEC_ETM == EWM models
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      
   getDefaultURef = getDefaultURef .AND. .NOT. IsUnusedParameter
   
   
   ! ------------ Read in the jet height -------------------------------------------------------------
   IsUnUsedParameter = TRIM(p%met%WindProfileType) /= 'JET'
   CALL ReadRVarDefault( UI, InFile, p%met%ZJetMax, "ZJetMax", "Jet height [m]", UnEc, getDefaultZJetMax, ErrStat2, ErrMsg2, IGNORE=IsUnusedParameter)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   getDefaultZJetMax = getDefaultZJetMax .AND. .NOT. IsUnusedParameter ! Jet height for 

   ! ------------ Read in the power law exponent, PLExp ---------------------------------------------
   IsUnusedParameter = (TRIM(p%met%WindProfileType) /= "PL" .AND. TRIM(p%met%WindProfileType) /= "IEC")
   CALL ReadRVarDefault( UI, InFile, p%met%PLExp, "PLExp", "Power law exponent [-]", UnEc, getDefaultPLExp, ErrStat2, ErrMsg2, IGNORE=IsUnusedParameter)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   getDefaultPLExp = getDefaultPLExp .AND. .NOT. IsUnusedParameter  ! we need RICH_NO before we can calculate a default for this value RICH_NO
                      
   ! ------------ Read in the surface roughness length, Z0 (that's z-zero) ---------------------------------------------
   IsUnusedParameter = p%met%TurbModel_ID==SpecModel_TIDAL
   CALL ReadRVarDefault( UI, InFile, p%met%Z0, "Z0", "Surface roughness length [m]", UnEc, UseDefault, ErrStat2, ErrMsg2, &
                                       IGNORE=IsUnusedParameter)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

!..................................................................................................................................
!  Do some error checking on the meteorology before we read the non-IEC meteorology data
!..................................................................................................................................   
   IF ( p%IEC%IEC_WindType == IEC_ETM .AND. p%IEC%ETMc <= 0. ) CALL SetErrStat( ErrID_Fatal, 'The ETM "c" parameter must be a positive number', ErrStat, ErrMsg, RoutineName)
      
      ! Make sure WindProfileType is valid for this turbulence model
   SELECT CASE ( TRIM(p%met%WindProfileType) )
      CASE ( 'JET' )
         IF ( p%met%TurbModel_ID /= SpecModel_GP_LLJ ) CALL SetErrStat( ErrID_Fatal, 'The jet wind profile is available with the GP_LLJ spectral model only.', ErrStat, ErrMsg, RoutineName)
      CASE ( 'LOG')
         IF (p%IEC%IEC_WindType /= IEC_NTM ) CALL SetErrStat( ErrID_Fatal, 'The IEC turbulence type must be NTM for the logarithmic wind profile.', ErrStat, ErrMsg, RoutineName)
      CASE ( 'PL'  ) !this is a valid WindProfileType
      CASE ( 'H2L' )
         IF ( p%met%TurbModel_ID /= SpecModel_TIDAL ) CALL SetErrStat( ErrID_Fatal, 'The "H2L" mean profile type can be used with only the "TIDAL" spectral model.', ErrStat, ErrMsg, RoutineName)
      CASE ( 'IEC' )
      CASE ( 'USR' )
         !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>         
         ! Get parameters for USR wind profile (so that we can use these parameters to get the wind speed later):
            CALL GetUSRProfiles( ProfileFile, p%met, UnEc, ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<         
         
      CASE ( 'TS' )
         IF ( p%met%TurbModel_ID /= SpecModel_TimeSer ) CALL SetErrStat( ErrID_Fatal, 'The "TS" mean profile type is valid only with the "TIMESR" spectral model.', ErrStat, ErrMsg, RoutineName)
      CASE ( 'API' )   ! ADDED BY Y.GUO
!bjj: I think we need to add some checks here??? MLB has comments about difference between RefHt and HubHt and 10 m         
      CASE DEFAULT
         CALL SetErrStat( ErrID_Fatal, 'The wind profile type must be "JET", "LOG", "PL", "IEC", "USR", "H2L", "TS", or default.' , ErrStat, ErrMsg, RoutineName)
   END SELECT

   IF ( p%met%TurbModel_ID == SpecModel_TIDAL .AND. TRIM(p%met%WindProfileType) /= "H2L" ) THEN
      p%met%WindProfileType = 'H2L'
      CALL SetErrStat( ErrID_Warn, 'Overwriting wind profile type to "H2L" for the "TIDAL" spectral model.', ErrStat, ErrMsg, RoutineName)
   ENDIF
        

   IF (p%met%KHtest) THEN
      IF ( p%met%TurbModel_ID /= SpecModel_NWTCUP ) CALL SetErrStat( ErrID_Fatal, 'The KH test can be used with the "NWTCUP" spectral model only.', ErrStat, ErrMsg, RoutineName)

      IF ( TRIM(p%met%WindProfileType) /= 'IEC' .AND. TRIM(p%met%WindProfileType) /= 'PL' ) THEN
         CALL SetErrStat( ErrID_Warn, 'Overwriting wind profile type for the KH test.', ErrStat, ErrMsg, RoutineName)         
         p%met%WindProfileType = 'IEC'
      ENDIF
      
      IF ( .NOT. p%WrFile(FileExt_CTS) ) THEN
         CALL SetErrStat( ErrID_Warn, 'Coherent turbulence time step files must be generated when using the "KHTEST" option.', ErrStat, ErrMsg, RoutineName)
         p%WrFile(FileExt_CTS)  = .TRUE.
      ENDIF           
      
      IF ( .NOT. EqualRealNos(p%met%PLExp, 0.3_ReKi) ) THEN
         CALL SetErrStat( ErrID_Warn, 'Overwriting the power law exponent for KH test.', ErrStat, ErrMsg, RoutineName)
         p%met%PLExp = 0.3
      ENDIF            
   END IF
   
   
   IF ( getDefaultURef ) THEN
      IF ( p%usr%NPoints > 0 ) THEN
         p%met%RefHt = p%usr%pointzi(p%usr%RefPtID)
         p%met%URef  = p%usr%meanU(p%usr%RefPtID,1)
         getDefaultURef = .FALSE.
      ELSEIF( TRIM(p%met%WindProfileType) /= 'JET' ) THEN
         ! Also note that if we specify a "default for Ustar , we cannot enter "default" for URef. Otherwise, we get circular logic. Will check for that later.
         CALL SetErrStat( ErrID_Fatal, 'URef can be "default" for only the "JET" WindProfileType.', ErrStat, ErrMsg, RoutineName)
      END IF
   END IF
         
   IF ( p%met%Z0 <= 0.0_ReKi ) CALL SetErrStat( ErrID_Fatal, 'The surface roughness length must be a positive number or "default".', ErrStat, ErrMsg, RoutineName)
   
   IF ( TRIM(p%met%WindProfileType) == 'JET' .AND. .NOT. getDefaultZJetMax ) THEN
      IF ( p%met%ZJetMax <  ZJetMax_LB .OR. p%met%ZJetMax > ZJetMax_UB )  THEN
         CALL SetErrStat( ErrID_Fatal, 'The height of the maximum jet wind speed must be between '//TRIM(num2lstr(ZJetMax_LB))//&
                                       ' and '//TRIM(num2lstr(ZJetMax_UB))//' m.', ErrStat, ErrMsg, RoutineName)
      ENDIF
   ENDIF
   
      
   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF   
   
   !.................................................
   ! overwrite RefHt and URef for cases where they are unused [USR wind profiles (or TS)]
   !.................................................
   
      ! NOTE: return on abortErrLev before calling getVelocity
   IF ( TRIM(p%met%WindProfileType) == 'USR' .OR. TRIM(p%met%WindProfileType) == 'TS') THEN ! for user-defined wind profiles, we overwrite RefHt and URef because they don't mean anything otherwise
         ! Calculate URef, which is UHub:
         ! note that the 2 "ref" values in the subroutine arguments aren't used for the USR wind profile type.
         ! (also, we do not necessarially know PLExp, yet, so we can't call this routine when we have "PL" or "IEC" wind profile types.)
      CALL getVelocity(p, p%met%URef, p%met%RefHt, p%met%RefHt, tmp, ErrStat2, ErrMsg2) 
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         
     !p%met%RefHt = p%grid%HubHt         bjj changed this on 23-sep-2014
      p%met%URef  = tmp      
   ELSEIF( p%IEC%IEC_WindType > IEC_ETM ) THEN !i.e., any of the EWM models: IEC_EWM1, IEC_EWM50, IEC_EWM100 
      p%met%RefHt = p%grid%HubHt
      p%met%URef  = p%IEC%VRef
   ENDIF   
      
      ! check that RefHt and URef are appropriate values:
   IF ( p%met%RefHt <=  0.0_ReKi )  CALL SetErrStat( ErrID_Fatal, 'The reference height must be greater than zero.', ErrStat, ErrMsg, RoutineName)         
   IF ( .NOT. getDefaultURef ) THEN
      IF ( p%met%URef <=  0.0_ReKi )  CALL SetErrStat( ErrID_Fatal, 'The reference wind speed must be greater than zero.', ErrStat, ErrMsg, RoutineName)
   ENDIF   ! Otherwise, we're using a Jet profile with default wind speed (for now it's -999.9)
   
   
   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF    
   
   !===============================================================================================================================
   ! Read the meteorological boundary conditions for non-IEC models. 
   !===============================================================================================================================

   CALL ReadCom( UI, InFile, "Non-IEC Meteorological Boundary Conditions Heading Line 1", ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   CALL ReadCom( UI, InFile, "Non-IEC Meteorological Boundary Conditions Heading Line 2", ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

      ! ------------ Read in the site latitude, LATITUDE ---------------------------------------------   
   IsUnusedParameter =  p%met%IsIECModel .AND. p%met%TurbModel_ID /= SpecModel_MODVKM  ! Used to caluculte z0 in ModVKM model; also used for default ZI
   CALL ReadRVarDefault( UI, InFile, p%met%Latitude, "Latitude", "Site latitude [degrees]", UnEc, UseDefault, ErrStat2, ErrMsg2, &
                                       IGNORE=IsUnusedParameter)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            
      ! ------------ Read in the gradient Richardson number, RICH_NO ---------------------------------------------
   CALL ReadVar( UI, InFile, p%met%Rich_No, "RICH_NO", "Gradient Richardson number",ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  
      
         ! Convert RICH_NO input to value that will be used in the code:
         
      IF ( p%met%KHtest ) THEN
         IF ( .NOT. EqualRealNos(p%met%Rich_No, 0.02_ReKi) ) THEN
            p%met%Rich_No = 0.02
            CALL SetErrStat( ErrID_Warn, 'Overwriting the Richardson Number for KH test.', ErrStat, ErrMsg, RoutineName)
         ENDIF
      ELSEIF ( p%met%TurbModel_ID == SpecModel_USRVKM ) THEN
         IF ( .NOT. EqualRealNos(p%met%Rich_No, 0.0_ReKi) ) THEN
            CALL SetErrStat( ErrID_Warn, 'Overwriting the Richardson Number for the '//TRIM(p%met%TurbModel)//' model.', ErrStat, ErrMsg, RoutineName)
            p%met%Rich_No = 0.0
         ENDIF
      ELSEIF ( p%met%TurbModel_ID == SpecModel_NWTCUP .OR. p%met%TurbModel_ID == SpecModel_GP_LLJ ) THEN
         p%met%Rich_No = MIN( MAX( p%met%Rich_No, -1.0_ReKi ), 1.0_ReKi )  ! Ensure that: -1 <= RICH_NO <= 1
      ELSEIF (p%met%IsIECModel) THEN
         p%met%Rich_No = 0.0                       ! Richardson Number in neutral conditions
      ENDIF
                     
      ! now that we have Rich_No, we can calculate ZL and L
      ! necessary for DefaultUStar(p)
   CALL Calc_MO_zL(p%met%TurbModel_ID, p%met%Rich_No, p%grid%HubHt, p%met%ZL, p%met%L )
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
         
      ! ------------ Read in the shear/friction velocity, Ustar ------------------------
   CALL ReadRVarDefault( UI, InFile, p%met%Ustar, "UStar", "Friction or shear velocity [m/s]", UnEc, UseDefault, ErrStat2, ErrMsg2, IGNORE=p%met%IsIECModel )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      
      !IF ( p%met%IsIECModel ) THEN     
      !   p%met%Ustar   = 0.0                       ! Shear or friction velocity
      !ELSE
         IF ( UseDefault ) THEN
            IF ( getDefaultURef ) THEN ! This occurs if "default" was entered for both GP_LLJ wind speed and UStar
               CALL SetErrStat( ErrID_Fatal, 'The reference wind speed and friction velocity cannot both be "default."', ErrStat, ErrMsg, RoutineName)
            ELSE
               CALL DefaultUStar(p)
            END IF
         END IF         
      !END IF
   
      
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>           
      ! Calculate Coriolis parameter from latitude   ( Used for default ZI )      
   p%met%Fc = 2.0 * Omega * SIN( ABS(p%met%Latitude*D2R) )  
   
      
      ! We need the hub-height wind speed to calculate default Reynold's Stresses.
      ! We have a few steps to take before we can get that wind speed:
      
      ! ***** Calculate power law exponent, if needed *****   
   IF ( getDefaultPLExp ) p%met%PLExp = DefaultPowerLawExp( p )            
      
      ! ***** Calculate parameters for Jet profile, if needed *****         
   IF ( TRIM(p%met%WindProfileType) == 'JET' ) THEN      
      IF ( getDefaultZJetMax ) CALL DefaultZJetMax(p, OtherSt_RandNum)       ! requires Rich_No, ZL, Ustar
      CALL getJetCoeffs( p, getDefaultURef, OtherSt_RandNum, ErrStat2, ErrMsg2) ! getDefault
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   END IF
      
      ! now that we know URef (in case getDefaultURef was true), set UstarDiab (used in ustar profile and default ZI):
   p%met%UstarDiab  = getUstarDiab(p%met%URef, p%met%RefHt, p%met%z0, p%met%ZL) 

   IF (ErrStat >= AbortErrLev) THEN  ! just in case we had a fatal error, let's check before calling getVelocity 
      CALL Cleanup()
      RETURN
   END IF  
   
   CALL getVelocity(p, p%met%URef, p%met%RefHt, p%grid%HubHt, tmp, ErrStat2, ErrMsg2)  
      p%UHub = tmp
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      
      
      ! We need the (local) Ustar at the hub-height; while we're at it, we'll 
      !   calculate the UstarOffset and UstarSlope it uses:
      
      
      ! Set up the heights for the zl- and ustar-profile averages across the rotor disk
   RotorDiskHeights  = (/ p%grid%HubHt-0.5*p%grid%RotorDiameter, p%grid%HubHt, p%grid%HubHt+0.5*p%grid%RotorDiameter /)        
   DO TmpIndex = 1,SIZE(RotorDiskHeights)  ! set height limits so we don't extrapolate too far
      RotorDiskHeights(TmpIndex) = MAX( MIN(RotorDiskHeights(TmpIndex), profileZmax), profileZmin)
   ENDDO

   IF (p%met%TurbModel_ID == SpecModel_GP_LLJ ) THEN
      p%met%UstarSlope = 1.0_ReKi         

      CALL getVelocityProfile(p, p%met%URef, p%met%RefHt, RotorDiskHeights, TmpUary, ErrStat2, ErrMsg2)    ! Set TmpUary
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)         
      TmpUstar  = getUStarProfile( P, TmpUary, RotorDiskHeights, 0.0_ReKi, p%met%UstarSlope )  ! set offset to 0 here <-
            
      p%met%UstarOffset = p%met%Ustar - SUM(TmpUstar) / SIZE(TmpUstar)    ! Ustar minus the average of those 3 points
      TmpUstar(:) = TmpUstar(:) + p%met%UstarOffset
   ELSE
      p%met%UstarSlope = 1.0_ReKi         

      TmpUary   = (/ 0.0_ReKi, 0.0_ReKi, 0.0_ReKi /)
      TmpUstar  = (/ 0.0_ReKi, 0.0_ReKi, 0.0_ReKi /)
         
      p%met%UstarOffset= 0.0_ReKi      
   ENDIF
                        
   
      ! Get the default mean spatial coherence models
   CALL GetDefaultSCMod( p%met%TurbModel_ID, p%met%SCMod )   
   
      ! Get the default mean Reynolds stresses
      ! (requires uHub, Ustar, Rich_No, ZL, TmpUStar)
   CALL GetDefaultRS(  p, OtherSt_RandNum, TmpUStar(2), ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      
      ! Default coherence parameters and IEC scaling parameters   
   CALL CalcIECScalingParams(p%IEC, p%grid%HubHt, p%UHub, p%met%InCDec, p%met%InCohB, p%met%TurbModel_ID, p%met%IsIECModel, ErrStat2, ErrMsg2)                  
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)      
      
   IF ( .NOT. p%met%IsIECModel  ) THEN
      CALL GetDefaultCoh( p%met%TurbModel_ID, p%met%RICH_NO, p%UHub, p%grid%HubHt, p%met%IncDec, p%met%InCohB )
   END IF
               
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   
      ! ------------- Read in the mixing layer depth, ZI ---------------------------------------------
   IsUnusedParameter = p%met%ZL >= 0.0_ReKi .AND. p%met%TurbModel_ID /= SpecModel_GP_LLJ ! used for unstable flows; GP_LLJ model may have both stable and unstable flows in its ZL_Profile      
   CALL ReadRVarDefault( UI, InFile, p%met%ZI, "ZI", "Mixing layer depth [m]", UnEc, UseDefault, ErrStat2, ErrMsg2, IGNORE=IsUnusedParameter )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   IF ( IsUnusedParameter ) THEN
      p%met%ZI = 999.9_ReKi ! set to a value > 0 that we don't care about
   ELSE
      IF ( UseDefault ) CALL DefaultMixingLayerDepth(p)      
   ENDIF
   
      
       ! ----------- Read in the mean hub u'w' Reynolds stress, PC_UW ---------------------------------------------
   CALL ReadRVarDefault( UI, InFile, p%met%PC_UW, "PC_UW", "Mean hub u'w' Reynolds stress", UnEc, UseDefault, &
                                            ErrStat2, ErrMsg2, IGNORE=p%met%IsIECModel, IGNORESTR = p%met%UWskip )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

      ! ------------ Read in the mean hub u'v' Reynolds stress, PC_UV ---------------------------------------------
   CALL ReadRVarDefault( UI, InFile, p%met%PC_UV, "PC_UV", "Mean hub u'v' Reynolds stress", UnEc, UseDefault, &
                                            ErrStat2, ErrMsg2, IGNORE=p%met%IsIECModel, IGNORESTR = p%met%UVskip )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

      ! ------------ Read in the mean hub v'w' Reynolds stress, PC_VW ---------------------------------------------
   CALL ReadRVarDefault( UI, InFile, p%met%PC_VW, "PC_VW", "Mean hub v'w' Reynolds stress", UnEc, UseDefault, &
                                            ErrStat2, ErrMsg2, IGNORE=p%met%IsIECModel, IGNORESTR = p%met%VWskip )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

            
   !===============================================================================================================================
   ! Read the spatial coherence model section. 
   !===============================================================================================================================

   CALL ReadCom( UI, InFile, "Spatial Coherence Models Heading Line 1", ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   CALL ReadCom( UI, InFile, "Spatial Coherence Models Heading Line 2", ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

      ! ------------ Read in the spatial coherence models, SCMod(1), SCMod(2), SCMod(3). ---------------------------------------------   

   DO I=1,3
      CALL ReadCVarDefault ( UI, InFile, Line, "SCMod"//TRIM(Num2LStr(I)), Comp(I)//"-component coherence model", UnEc, UseDefault, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         IF ( .NOT. UseDefault ) THEN
            SELECT CASE ( TRIM(Line) )
               CASE("GENERAL")
                  p%met%SCMod(I) = CohMod_GENERAL
               CASE ("IEC")
                  p%met%SCMod(I) = CohMod_IEC
               CASE ("NONE")
                  p%met%SCMod(I) = CohMod_NONE
               CASE ("API")
                  p%met%SCMOD(I) = CohMod_API
                  IF (I /= 1) CALL SetErrStat( ErrID_Fatal, "API coherence model is valid only for the u-component", ErrStat, ErrMsg, RoutineName)
               CASE DEFAULT
                  p%met%SCMod(I) = CohMod_NONE
                  IF (I==1) THEN
                     CALL SetErrStat( ErrID_Fatal, 'Unknown value for SCMod'//TRIM(Num2LStr(I))//'. Valid entries are "GENERAL","IEC","API", or "NONE".', ErrStat, ErrMsg, RoutineName)
                  ELSE               
                     CALL SetErrStat( ErrID_Fatal, 'Unknown value for SCMod'//TRIM(Num2LStr(I))//'. Valid entries are "GENERAL","IEC", or "NONE".', ErrStat, ErrMsg, RoutineName)
                  END IF 
            END SELECT
         END IF
   END DO
         
      ! ------------ Read in the u component coherence parameters, InCDec(1) and InCohB(1) ------------
   CALL ReadRAryDefault( UI, InFile, InCVar, "InCDec1", "u-component coherence parameters", UnEc, UseDefault, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      IF ( .NOT. UseDefault ) THEN         
         p%met%InCDec(1) = InCVar(1)
         p%met%InCohB(1) = InCVar(2)
      END IF

      ! ------------ Read in the v component coherence parameters, InCDec(2) and InCohB(2) ----------
   CALL ReadRAryDefault( UI, InFile, InCVar, "InCDec2", "v-component coherence parameters", UnEc, UseDefault, ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      IF ( .NOT. UseDefault ) THEN       ! these are the values we just read in  
         p%met%InCDec(2) = InCVar(1)
         p%met%InCohB(2) = InCVar(2)
      END IF

      ! ------------ Read in the w component coherence parameters, InCDec(3) and InCohB(3) -------
   CALL ReadRAryDefault( UI, InFile, InCVar, "InCDec3", "w-component coherence parameters", UnEc, UseDefault, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      IF ( .NOT. UseDefault ) THEN         
         p%met%InCDec(3) = InCVar(1)
         p%met%InCohB(3) = InCVar(2)
      END IF

      ! ------------ Read in the coherence exponent, COHEXP -----------------------------------
   CALL ReadRVarDefault( UI, InFile, p%met%CohExp, "CohExp", "Coherence exponent", UnEc, UseDefault, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

      
!..................................................................................................................................
!  Do some error checking on the non-IEC meteorology data and coherence
!..................................................................................................................................         
      
   IF ( .NOT. p%met%IsIECModel ) THEN
      IF ( ABS(p%met%Latitude) < 5.0 .OR. ABS(p%met%Latitude) > 90.0 ) THEN
         CALL SetErrStat( ErrID_Fatal, 'The latitude must be between -90 and 90 degrees but not between -5 and 5 degrees.', ErrStat, ErrMsg, RoutineName)
      ENDIF

      IF (p%met%Ustar   <= 0.0_ReKi) CALL SetErrStat( ErrID_Fatal, 'The friction velocity must be a positive number.', ErrStat, ErrMsg, RoutineName)
      IF ( p%met%ZI     <= 0.0_ReKi) CALL SetErrStat( ErrID_Fatal, 'The mixing layer depth must be a positive number for unstable flows.', ErrStat, ErrMsg, RoutineName)         
   END IF
   
   IF ( p%met%COHEXP <  0.0_ReKi) CALL SetErrStat( ErrID_Fatal, 'The coherence exponent must be non-negative.', ErrStat, ErrMsg, RoutineName)

   DO I = 1,3
      IF ( p%met%InCDec(I) <= 0.0_ReKi ) CALL SetErrStat( ErrID_Fatal, 'The '//Comp(I)//'-component coherence decrement must be a positive number.', ErrStat, ErrMsg, RoutineName)
   END DO
           
!..................................................................................................................................
!  Calculate zlOffset
!  Adjust UstarSlope and UstarOffset based on entered PC_UW
!..................................................................................................................................     
      TmpZLary = getZLProfile(TmpUary, RotorDiskHeights, p%met%Rich_No, p%met%ZL, p%met%L, 0.0_ReKi, p%met%WindProfileType)
      p%met%zlOffset = p%met%ZL - SUM(TmpZLary) / SIZE(TmpZLary)

         ! Modify previously calculated UstarSlope and UstarOffset based on our input (target) Reynolds' stress values
      IF (.NOT. p%met%UWskip) THEN
         TmpUstarD = ( TmpUstar(1)- 2.0*TmpUstar(2) + TmpUstar(3) )

         IF ( .NOT. EqualRealNos( TmpUstarD, 0.0_ReKi ) ) THEN
            p%met%UstarSlope  = 3.0*(p%met%Ustar -  SQRT( ABS(p%met%PC_UW) ) ) / TmpUstarD
            p%met%UstarOffset = SQRT( ABS(p%met%PC_UW) ) - p%met%UstarSlope*(TmpUstar(2) - p%met%UstarOffset)
         ELSE
            p%met%UstarSlope  = 0.0
            p%met%UstarOffset = SQRT( ABS(p%met%PC_UW) )
         ENDIF
      ENDIF   
         
      
   
   !===============================================================================================================================
   ! Read the Coherent Turbulence Scaling Parameters, if necessary.  
   !===============================================================================================================================
IF ( .NOT. p%met%IsIECModel  ) THEN

   IF ( p%WrFile(FileExt_CTS) ) THEN

      CALL ReadCom( UI, InFile, "Coherent Turbulence Scaling Parameters Heading Line 1", ErrStat2, ErrMsg2, UnEc )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

      CALL ReadCom( UI, InFile, "Coherent Turbulence Scaling Parameters Heading Line 2", ErrStat2, ErrMsg2, UnEc )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)


         ! ------------ Read the name of the path containg event file definitions, CTEventPath --------------------------

      CALL ReadVar( UI, InFile, p%CohStr%CTEventPath, "CTEventPath", "Coherence events path",ErrStat2, ErrMsg2, UnEc)
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

      CALL ReadVar( UI, InFile, Line, "CTEventFile", "Event file type",ErrStat2, ErrMsg2, UnEc)
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)


      IF ( p%met%KHtest ) THEN

         p%CohStr%CText = 'les'
         p%CohStr%CTEventFile = TRIM(p%CohStr%CTEventPath)//PathSep//'Events.xtm'

         CALL WrScr( ' LES events will be used for the KH test.' )

      ELSE

         p%CohStr%CText = Line  !This will preserve the case formatting, in case it matters.

         CALL Conv2UC( Line )

         IF (Line(1:6) == "RANDOM") THEN
             CALL RndUnif( p%RNG, OtherSt_RandNum, tmp )

             IF ( tmp <= 0.5 ) THEN
                 p%CohStr%CText = 'les'
             ELSE
                 p%CohStr%CText = 'dns'
             ENDIF
         ENDIF

         p%CohStr%CTEventFile = TRIM(p%CohStr%CTEventPath)//PathSep//'Events.'//TRIM(p%CohStr%CText)

      ENDIF


         ! ------------ Read the Randomization Flag, Randomize -----------------------------------
      CALL ReadVar( UI, InFile, Randomize, "Randomize", "Randomize CT Scaling",ErrStat2, ErrMsg2, UnEc)
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

         ! ------------ Read the Disturbance Scale, DistScl ---------------------------------------------
      CALL ReadVar( UI, InFile, p%CohStr%DistScl, "DistScl", "Disturbance scale",ErrStat2, ErrMsg2, UnEc)
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

         ! ------------ Read the Lateral Fractional Location of tower centerline in wave, CTLy ----------
      CALL ReadVar( UI, InFile, p%CohStr%CTLy, "CTLy", "Location of tower centerline",ErrStat2, ErrMsg2, UnEc)
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

         ! ------------ Read the Vertical Fraction Location of hub in wave, CTLz ------------------------
      CALL ReadVar( UI, InFile, p%CohStr%CTLz, "CTLz", "Location of hub height",ErrStat2, ErrMsg2, UnEc)
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

      IF ( p%met%KHtest ) THEN
         p%CohStr%DistScl = 1.0
         p%CohStr%CTLy    = 0.5
         p%CohStr%CTLz    = 0.5
         Randomize = .FALSE.
         CALL SetErrStat( ErrID_Info, 'Billow will cover rotor disk for KH test.', ErrStat, ErrMsg, RoutineName)
            
      ELSEIF ( Randomize ) THEN

         CALL RndUnif( p%RNG, OtherSt_RandNum, tmp )

            ! Assume a 75% chance of coherent turbulence being the size of the rotor
            ! If the rotor is too small, assume a 100% chance.
            ! If the turbulence is not the size of the rotor, assume it's half the size
            ! of the disk, with equal probablilty of being in the lower or upper half.

         IF ( tmp > 0.25 .OR. p%grid%RotorDiameter <= 30.0 ) THEN

            p%CohStr%DistScl = 1.0
            p%CohStr%CTLy    = 0.5
            p%CohStr%CTLz    = 0.5

         ELSE

            p%CohStr%DistScl = 0.5
            p%CohStr%CTLy    = 0.5

            IF ( tmp < 0.125 ) THEN
               p%CohStr%CTLz = 0.0 ! The hub is on the bottom of the dataset (i.e. the billow is on the top of the disk)
            ELSE
               p%CohStr%CTLz = 1.0 ! The hub is on the top of the dataset
            ENDIF

         ENDIF

      ELSE  !Don't randomize:

         IF ( p%CohStr%DistScl < 0.0 ) THEN
            CALL SetErrStat( ErrID_Fatal, 'The disturbance scale must be a positive.', ErrStat, ErrMsg, RoutineName)         
         ELSEIF ( p%grid%RotorDiameter <= 30.0 .AND. p%CohStr%DistScl < 1.0 ) THEN
            CALL SetErrStat( ErrID_Fatal, 'The disturbance scale must be at least 1.0 for rotor diameters less than 30.', ErrStat, ErrMsg, RoutineName)         
         ELSEIF ( p%grid%RotorDiameter*p%CohStr%DistScl <= 15.0  ) THEN
            CALL SetErrStat( ErrID_Fatal, 'The coherent turbulence must be greater than 15 meters in height.  '//&
                        'Increase the rotor diameter or the disturbance scale. ', ErrStat, ErrMsg, RoutineName)         
         ENDIF

      ENDIF


         ! ---------- Read the Minimum event start time, CTStartTime --------------------------------------------
      CALL ReadVar( UI, InFile, p%CohStr%CTStartTime, "CTStartTime", "CTS Start Time",ErrStat2, ErrMsg2, UnEc)
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

      p%CohStr%CTStartTime = MAX( p%CohStr%CTStartTime, 0.0_ReKi ) ! A Negative start time doesn't really make sense...

   ENDIF    ! WrFile(FileExt_CTS)


ELSE  ! IECVKM, IECKAI, MODVKM, OR API models
   
   IF ( p%IEC%NumTurbInp .AND. EqualRealNos( p%IEC%PerTurbInt, 0.0_ReKi ) ) THEN    ! This will produce constant winds, instead of an error when the transfer matrix is singular
      p%met%TurbModel = 'NONE'
      p%met%TurbModel_ID = SpecModel_NONE
   ENDIF

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
SUBROUTINE OpenSummaryFile(RootName, US, DescStr, ErrStat, ErrMsg)

  ! This subroutine is used to open the summary output file.

IMPLICIT         NONE

   INTEGER(IntKi),                  INTENT(INOUT) :: US                              ! Unit specifier for summary file
   CHARACTER(*),                    INTENT(IN   ) :: RootName                        ! rootname of the primary TurbSim input file
   CHARACTER(*),                    INTENT(  OUT) :: DescStr                         ! string describing time TurbSim files were generated
   INTEGER(IntKi),                  INTENT(  OUT) :: ErrStat                         ! Error level
   CHARACTER(*),                    INTENT(  OUT) :: ErrMsg                          ! Message describing error




   ! Open summary file.
CALL GetNewUnit( US, ErrStat, ErrMsg )
CALL OpenFOutFile( US, TRIM( RootName )//'.sum', ErrStat, ErrMsg ) ! Formatted output file
if (ErrStat >= AbortErrLev) then
   US = -1
   RETURN
end if




   ! Let's make a string so that the binary file and the full-field file have the same date and time:
DescStr = 'generated by '//TRIM( GetNVD(TurbSim_Ver) )//' on '//CurDate()//' at '//CurTime()//'.'

   ! Write the program name and version, date and time into the summary file.
WRITE (US,"( / 'This summary file was ', A / )")  TRIM(DescStr)

   ! Capitalize the first letter of the string and save it for the full-field files.
DescStr = 'This full-field file was '//TRIM(DescStr)


RETURN
END SUBROUTINE OpenSummaryFile
!=======================================================================
SUBROUTINE GetUSRProfiles(FileName, p_met, UnEc, ErrStat, ErrMsg)

   IMPLICIT                              NONE

   TYPE(Meteorology_ParameterType), INTENT(INOUT) :: p_met
   INTEGER(IntKi),                  INTENT(IN   ) :: UnEc                            ! echo file unit number
   INTEGER(IntKi),                  INTENT(  OUT) :: ErrStat                         ! Error level
   CHARACTER(*),                    INTENT(  OUT) :: ErrMsg                          ! Message describing error
   CHARACTER(*),                    INTENT(IN   ) :: FileName                        ! Name of the input file
   
   
   ! local variables
   
   INTEGER                                        :: U_in                           ! Input unit.
   INTEGER(IntKi)                                 :: ErrStat2                        ! Error level (local)
   CHARACTER(MaxMsgLen)                           :: ErrMsg2                         ! Message describing error (local)
   character(*), parameter                        :: RoutineName = 'GetUSRProfiles'
   CHARACTER(200)                                 :: TempWarn
   
!   CHARACTER(200)                                 :: LINE
                                                  
   REAL(ReKi)                                     :: L_Usr_Tmp
   REAL(ReKi)                                     :: Sigma_USR_Tmp
   REAL(ReKi)                                     :: U_USR_Tmp
   REAL(ReKi)                                     :: WindDir_USR_Tmp
   REAL(ReKi)                                     :: Z_USR_Tmp
                                                  
   INTEGER                                        :: I
   INTEGER                                        :: Indx
   INTEGER                                        :: J
                                                  
   LOGICAL                                        :: ReadSigL                       ! Whether or not to read the last 2 columns

   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   U_in = -1
   CALL GetNewUnit( U_in, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   CALL OpenFInpFile( U_in, FileName, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      
   IF (ErrStat >= AbortErrLev) THEN
      CLOSE(U_in)
      RETURN
   END IF
               
   DO I=1,3
      CALL ReadCom( U_in, FileName, "Header line "//trim(num2lstr(I))//" for user-defined profiles", ErrStat2, ErrMsg2, UnEc )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   END DO
   
   
      ! ---------- Read the size of the arrays --------------------------------------------
   CALL ReadVar( U_in, FileName, p_met%NumUSRz, "NumUSRz", "Number of heights in the user-defined profiles", ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   IF ( p_met%NumUSRz < 1 ) THEN
      CALL SetErrStat( ErrID_Fatal, 'The number of heights specified in the user-defined profiles must be at least 1.', ErrStat, ErrMsg, RoutineName)
   ENDIF

   DO I=1,3
         ! ---------- Read the scaling for the standard deviations --------------------------------------------
      CALL ReadVar( U_in, FileName, p_met%USR_StdScale(I), "USR_StdScale", "Scaling value for user-defined standard deviation profile", ErrStat2, ErrMsg2, UnEc )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)


      IF ( p_met%USR_StdScale(I) <= 0. ) THEN
         CALL SetErrStat( ErrID_Fatal, 'The scaling value for the user-defined standard deviation profile must be positive.', ErrStat, ErrMsg, RoutineName)
      ENDIF
   ENDDO

      ! Allocate the data arrays
   CALL AllocAry(p_met%USR_Z,       p_met%NumUSRz, 'USR_Z (user-defined height)',               ErrStat2, ErrMsg2); CALL SetErrStat(ErrSTat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   CALL AllocAry(p_met%USR_U,       p_met%NumUSRz, 'USR_U (user-defined wind speed)',           ErrStat2, ErrMsg2); CALL SetErrStat(ErrSTat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   CALL AllocAry(p_met%USR_WindDir, p_met%NumUSRz, 'USR_WindDir (user-defined wind direction)', ErrStat2, ErrMsg2); CALL SetErrStat(ErrSTat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)


   IF ( p_met%TurbModel_ID == SpecModel_USRVKM ) THEN
      ReadSigL = .TRUE.

      CALL AllocAry(p_met%USR_Sigma, p_met%NumUSRz, 'USR_Sigma (user-defined sigma)',    ErrStat2, ErrMsg2); CALL SetErrStat(ErrSTat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      CALL AllocAry(p_met%USR_L,     p_met%NumUSRz, 'USR_L (user-defined length scale)', ErrStat2, ErrMsg2); CALL SetErrStat(ErrSTat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      
   ELSE
      ReadSigL = .FALSE.
   ENDIF

   IF (ErrStat >= AbortErrLev) THEN
      CLOSE(U_in)
      RETURN
   END IF
   
      ! ---------- Skip 4 lines --------------------------------------------
   DO I=1,4
      CALL ReadCom( U_in, FileName, "Headers for user-defined variables", ErrStat2, ErrMsg2, UnEc )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   ENDDO

   DO I=1,p_met%NumUSRz

      IF ( ReadSigL ) THEN
         READ( U_in, *, IOSTAT=ErrStat2 ) p_met%USR_Z(I), p_met%USR_U(I), p_met%USR_WindDir(I), p_met%USR_Sigma(I), p_met%USR_L(I)
      ELSE
         READ( U_in, *, IOSTAT=ErrStat2 ) p_met%USR_Z(I), p_met%USR_U(I), p_met%USR_WindDir(I)
      ENDIF

      IF ( ErrStat2 /= 0 ) THEN
         CALL SetErrStat( ErrID_Fatal, 'Could not read entire user-defined variable list on line '//Int2LStr(I)//'.', ErrStat, ErrMsg, RoutineName)
         CLOSE(U_in)
         RETURN
      ENDIF

      TempWarn = 'Error reading user-defined variable list on line '//Int2LStr(I)
      CALL CheckRealVar( p_met%USR_Z(I), TempWarn, ErrStat2, ErrMsg2 )
         call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      CALL CheckRealVar( p_met%USR_U(I), TempWarn, ErrStat2, ErrMsg2 )
         call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      CALL CheckRealVar( p_met%USR_WindDir(I), TempWarn, ErrStat2, ErrMsg2 )
         call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      IF ( ReadSigL ) THEN
         IF ( p_met%USR_Sigma(I) <= REAL( 0., ReKi ) ) THEN
            CALL SetErrStat( ErrID_Fatal, 'The standard deviation must be a positive number.', ErrStat, ErrMsg, RoutineName)
         ELSEIF ( p_met%USR_L(I) <= REAL( 0., ReKi ) ) THEN
            CALL SetErrStat( ErrID_Fatal, 'The length scale must be a positive number.', ErrStat, ErrMsg, RoutineName)
         ENDIF
         
         CALL CheckRealVar( p_met%USR_Sigma(I), TempWarn, ErrStat2, ErrMsg2 )
            call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         CALL CheckRealVar( p_met%USR_L(I), TempWarn, ErrStat2, ErrMsg2 )
            call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      ENDIF
      if (ErrStat >= AbortErrLev) then
         CLOSE(U_in)
         RETURN
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
               CALL SetErrStat( ErrID_Fatal, 'User-defined values must contain unique heights.', ErrStat, ErrMsg, RoutineName)
               CLOSE(U_in)
               RETURN
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

   CLOSE(U_in)

END SUBROUTINE GetUSRProfiles
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
   CHARACTER(*), PARAMETER                        :: RoutineName = 'GetUSRSpec'
   CHARACTER(200)                                 :: TempWarn

   ErrStat = ErrID_None
   ErrMSg  = ""
   
      ! --------- Open the file ---------------

   CALL GetNewUnit( USpec, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2 , ErrStat, ErrMsg, RoutineName)
      
   CALL OpenFInpFile( USpec, FileName, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2 , ErrStat, ErrMsg, RoutineName)
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      ENDIF


   CALL WrScr1(' Reading the user-defined spectra input file "'//TRIM(FileName)//'".' )


      ! --------- Read the comment lines at the beginning of the file ---------------
   DO I=1,3
      CALL ReadCom( USpec, FileName, "user-spectra header line #"//TRIM(Num2LStr(I)), ErrStat2, ErrMsg2, UnEc)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   ENDDO


      ! ---------- Read the size of the arrays --------------------------------------------
   CALL ReadVar( USpec, FileName, p%usr%nFreq, "nFreq", "Number of frequencies in the user-defined spectra", ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      ENDIF


   DO I=1,3
         ! ---------- Read the scaling for the arrays --------------------------------------------
      CALL ReadVar( USpec, FileName, SpecScale(I), "SpecScale", "Scaling value for user-defined standard deviation profile", ErrStat2, ErrMsg2, UnEc )
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   ENDDO
      
   IF ( p%usr%nFreq < 3      ) CALL SetErrStat(ErrID_Fatal, 'The number of frequencies specified in the user-defined spectra must be at least 3.' , ErrStat, ErrMsg, RoutineName)
   IF ( ANY(SpecScale <= 0.) ) CALL SetErrStat(ErrID_Fatal, 'The scaling value for the user-defined spectra must be positive.' , ErrStat, ErrMsg, RoutineName)
   
   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   ENDIF   
   
      ! Allocate the data arrays
   CALL AllocAry( p%usr%f,      p%usr%nFreq,    'f (user-defined frequencies)'  ,ErrStat2,ErrMsg2); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   CALL AllocAry( p%usr%S,      p%usr%nFreq,1,3,'S (user-defined spectra)'      ,ErrStat2,ErrMsg2); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   CALL AllocAry( p%usr%pointzi, iPoint        , 'pointzi (user-defined spectra',ErrStat2,ErrMsg2); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)   
   
   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF
   
   p%usr%pointzi = 0.0_ReKi ! we don't care what this is; it's only potentially used so we can use the same interpolation routine as the user time-series input

      ! ---------- Skip 4 lines --------------------------------------------
   DO I=1,4
      CALL ReadCom( USpec, FileName, "Headers for user-defined variables", ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   ENDDO

      ! ---------- Read the data lines --------------------------------------
   DO I=1,p%usr%nFreq
      TempWarn = 'Error reading user-defined spectra line '//Int2LStr(I)
      
      READ( USpec, *, IOSTAT=ErrStat2 ) p%usr%f(I), p%usr%S(I,iPoint,1), p%usr%S(I,iPoint,2), p%usr%S(I,iPoint,3)

      IF ( ErrStat2 /= 0 ) THEN
         CALL SetErrStat(ErrID_Fatal, 'Could not read entire user-defined spectra on line '//Int2LStr(I)//'.' , ErrStat, ErrMsg, RoutineName)
         CALL Cleanup()
         RETURN
      ENDIF

      IF ( ANY( p%usr%S(I,iPoint,:) <=  0._ReKi ) ) THEN

         CALL SetErrStat(ErrID_Fatal, 'The spectra must contain positive numbers.' , ErrStat, ErrMsg, RoutineName)
         CALL Cleanup()
         RETURN
         
!      ELSEIF ( p%usr%f(I) <= 0.0_ReKi ) THEN
!         CALL SetErrStat(ErrID_Fatal, 'The frequencies must be positive numbers.' , ErrStat, ErrMsg, RoutineName)
!         CALL Cleanup()
!         RETURN
      ENDIF

      call CheckRealVar( p%usr%f(I), TempWarn, ErrStat2, ErrMsg2 )
         call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      call CheckRealVar( p%usr%S(I,iPoint,1), TempWarn, ErrStat2, ErrMsg2 )
         call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      call CheckRealVar( p%usr%S(I,iPoint,2), TempWarn, ErrStat2, ErrMsg2 )
         call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      call CheckRealVar( p%usr%S(I,iPoint,3), TempWarn, ErrStat2, ErrMsg2 )
         call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      IF (ErrStat >= AbortErrLev) then
         CALL Cleanup()
         RETURN
      end if
      
      
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
               CALL SetErrStat(ErrID_Fatal, 'Error: user-defined spectra must contain unique frequencies.' , ErrStat, ErrMsg, RoutineName)
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
SUBROUTINE GetUSRTimeSeries(FileName, p, UnEc, ErrStat, ErrMsg)

   IMPLICIT NONE

   TYPE(TurbSim_ParameterType),     INTENT(INOUT) :: p                              !< Simulation parameters
   INTEGER(IntKi),                  INTENT(IN   ) :: UnEc                           !< Echo file unit number
   INTEGER(IntKi),                  INTENT(  OUT) :: ErrStat                        !< Error level
   CHARACTER(*),                    INTENT(  OUT) :: ErrMsg                         !< Message describing error
   CHARACTER(*),                    INTENT(IN)    :: FileName                       !< Name of the input file

      ! local variables
   real(reKi)                                     :: tmpAry(2)
   real(ReKi)                                     :: dt                             ! difference between consecutive times entered in the file (must be constant)
   INTEGER(IntKi), PARAMETER                      :: NumLinesBeforeTS = 11          ! Number of lines in the input file before the time series start (need to add nPoint lines). IMPORTANT:  any changes to the number of lines in the header must be reflected here
      
   INTEGER(IntKi)                                 :: UnIn                           ! unit number for reading input file
   INTEGER(IntKi)                                 :: I, J                           ! loop counters
   INTEGER(IntKi)                                 :: IPoint                         ! loop counter on number of points
   INTEGER(IntKi)                                 :: IVec                           ! loop counter on velocity components being read
   INTEGER(IntKi)                                 :: ErrStat2                       ! Error level (local)
   CHARACTER(MaxMsgLen)                           :: ErrMsg2                        ! Message describing error (local)
   CHARACTER(*), parameter                        :: RoutineName = 'GetUSRTimeSeries'
   
   CHARACTER(200)                                 :: FormStr          
   CHARACTER(200)                                 :: TempWarn
   CHARACTER(1)                                   :: tmpChar          
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
      ! --------- Open the file ---------------

   CALL GetNewUnit( UnIn, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2 , ErrStat, ErrMsg, RoutineName)
      
   CALL OpenFInpFile( UnIn, FileName, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2 , ErrStat, ErrMsg, RoutineName)
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      ENDIF


   CALL WrScr1(' Reading the user-defined time-series input file "'//TRIM(FileName)//'".' )
   
   IF ( UnEc > 0 )  WRITE (UnEc,'(/,A,/)')  'Data from '//TRIM(TurbSim_Ver%Name)//' user time-series input file "'//TRIM( FileName )//'":'

   
   do i=1,3
      CALL ReadCom( UnIn, FileName, "Header #"//TRIM(Num2Lstr(i))//"for user time-series input", ErrStat2, ErrMsg2, UnEc )   
         CALL SetErrStat(ErrStat2, ErrMsg2 , ErrStat, ErrMsg, RoutineName)
   end do
      
   CALL ReadVar( UnIn, FileName, p%usr%nComp, 'nComp', 'How many velocity components will be input? (1=u component only; 2=u&v components; 3=u,v,and w)', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat(ErrStat2, ErrMsg2 , ErrStat, ErrMsg, RoutineName)
         
   CALL ReadVar( UnIn, FileName, p%usr%nPoints, 'nPoints', 'Number of time series points contained in this file', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat(ErrStat2, ErrMsg2 , ErrStat, ErrMsg, RoutineName)
         
   CALL ReadVar( UnIn, FileName, p%usr%RefPtID, 'RefPtID', 'Index of the reference point (1-nPoints)', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat(ErrStat2, ErrMsg2 , ErrStat, ErrMsg, RoutineName)

   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF
         
   IF ( p%usr%RefPtID < 1 .OR. p%usr%RefPtID > p%usr%nPoints ) THEN
      CALL SetErrStat(ErrID_Fatal, 'RefPtID must be between 1 and nPoints (inclusive).', ErrStat, ErrMsg, RoutineName)
      CALL Cleanup()
      RETURN
   END IF
         
   
   CALL AllocAry(p%usr%Pointyi, p%usr%nPoints, 'Pointyi', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2 , ErrStat, ErrMsg, RoutineName)
   CALL AllocAry(p%usr%Pointzi, p%usr%nPoints, 'Pointzi', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2 , ErrStat, ErrMsg, RoutineName)

   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF
      
   do i=1,2
      CALL ReadCom( UnIn, FileName, "Point location header #"//TRIM(Num2Lstr(i)), ErrStat2, ErrMsg2, UnEc )   
         CALL SetErrStat(ErrStat2, ErrMsg2 , ErrStat, ErrMsg, RoutineName)
   end do
      
   do iPoint=1,p%usr%nPoints
      CALL ReadAry( UnIn, FileName, TmpAry, 2, "point"//trim(Num2Lstr(iPoint)), "locations of points", ErrStat2, ErrMsg2, UnEc )
         CALL SetErrStat(ErrStat2, ErrMsg2 , ErrStat, ErrMsg, RoutineName)
         
      p%usr%Pointyi(iPoint) = TmpAry(1)
      p%usr%Pointzi(iPoint) = TmpAry(2)         
   end do
      
         
   do i=1,3
      CALL ReadCom( UnIn, FileName, "Time Series header #"//TRIM(Num2Lstr(i)), ErrStat2, ErrMsg2, UnEc )   
         CALL SetErrStat(ErrStat2, ErrMsg2 , ErrStat, ErrMsg, RoutineName)
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
      CALL SetErrStat(ErrID_Fatal, 'The user time-series input file must contain at least 2 rows of time data.', ErrStat, ErrMsg, RoutineName)
      CALL Cleanup()
      RETURN
   END IF
   
      ! now rewind and skip the first few lines. 
   REWIND( UnIn, IOSTAT=ErrStat2 )  
      IF (ErrStat2 /= 0_IntKi ) THEN
         CALL SetErrStat( ErrID_Fatal, 'Error rewinding file "'//TRIM(FileName)//'".', ErrStat, ErrMsg, RoutineName)   
         CALL Cleanup()
      END IF 

      !IMPORTANT: any changes to the number of lines in the header must be reflected in NumLinesBeforeTS
   DO I=1,NumLinesBeforeTS + p%usr%nPoints        
      READ( UnIn, '(A)', IOSTAT=ErrStat2 ) TmpChar   ! I'm going to ignore this error because we should have caught any issues the first time we read the file.      
   END DO
   
   !.......
   
   if (p%usr%nComp < 1 .OR. p%usr%nComp > 3) then
      CALL SetErrStat( ErrID_Fatal, 'Number of velocity components in file must be 1, 2 or 3.', ErrStat, ErrMsg, RoutineName)   
      CALL Cleanup()
   END IF 
   
   
   
   CALL AllocAry(p%usr%t, p%usr%nTimes,                             't', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2 , ErrStat, ErrMsg, RoutineName)
   CALL AllocAry(p%usr%v, p%usr%nTimes, p%usr%nPoints, p%usr%nComp, 'v', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2 , ErrStat, ErrMsg, RoutineName)
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF
   
      
   DO i=1,p%usr%nTimes
      TempWarn = 'Error reading from time series line '//trim(num2lstr(i))
      
      READ( UnIn, *, IOSTAT=ErrStat2 ) p%usr%t(i), ( (p%usr%v(i,iPoint,iVec), iVec=1,p%usr%nComp), iPoint=1,p%usr%nPoints )
      IF (ErrStat2 /=0) THEN
         CALL SetErrStat( ErrID_Fatal, trim(TempWarn)//'.', ErrStat, ErrMsg, RoutineName)
         CALL Cleanup()
         RETURN
      END IF

      CALL CheckRealVar( p%usr%t(i), TempWarn, ErrStat2, ErrMsg2 )
         call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      
      do iPoint=1,p%usr%nPoints
         do iVec=1,p%usr%nComp
            call CheckRealVar( p%usr%v(i,iPoint,iVec), TempWarn, ErrStat2, ErrMsg2 )
               call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
               if (ErrStat >= AbortErrLev) then
                  call Cleanup()
                  return
               end if
         end do
      end do
            
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
            CALL SetErrStat(ErrID_Fatal, 'Locations of points specified in the user time-series input file must be unique.', ErrStat, ErrMsg, RoutineName)
            CALL Cleanup()
            RETURN
         END IF
      end do
      
      !bjj: fix this in the future. Currently the interpolation routine won't work if z is not ordered properly. Also, interpolation doesn't take y into account, so we may want to fix that.
      IF ( p%usr%Pointzi(i) < p%usr%Pointzi(i-1) ) THEN
         CALL SetErrStat(ErrID_Fatal, 'The current implementation of user time-series input requires that the points be entered in the order of increasing height.', ErrStat, ErrMsg, RoutineName)
         CALL Cleanup()
         RETURN
      END IF      
   end do
   

   
   !DO i = 2,p%usr%nTimes
   !   IF (.NOT. EqualRealNos( p%usr%t(i-1) + p%grid%TimeStep, p%usr%t(i) ) ) THEN
   !      call SetErrStat(ErrID_Fatal, 'the delta time in the file must be constant and must be equal to input file variable TimeStep.', ErrStat, ErrMsg, RoutineName)
   !      EXIT
   !   END IF
   !END DO
   
      ! check for constant delta t:
   
   dt = p%usr%t(2) - p%usr%t(1)
   
   DO i = 3,p%usr%nTimes
      IF (.NOT. EqualRealNos( p%usr%t(i-1) + dt, p%usr%t(i) ) ) THEN
         call SetErrStat(ErrID_Fatal, 'The time between each row in the file must be constant.', ErrStat, ErrMsg, RoutineName)
         EXIT
      END IF
   END DO
   
   
   if ( .NOT. EqualRealNos( dt, p%grid%TimeStep ) ) THEN
      call SetErrStat(ErrID_Fatal, 'In this version of TurbSim, TimeStep must be the same as the delta time in the user-input time series file.', ErrStat, ErrMsg, RoutineName)
   end if
   
   
   
   CALL Cleanup()
   RETURN
   
CONTAINS
!...............................................
   SUBROUTINE Cleanup()
   
      CLOSE( UnIn )
   
   
   END SUBROUTINE Cleanup
!...............................................
END SUBROUTINE GetUSRTimeSeries 
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

!      CALL WrScr ( '    A default value will be used for '//TRIM(VarName)//'.' )
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
   INTEGER                      :: i
   INTEGER                      :: ErrStat2
   CHARACTER(ErrMsgLen)         :: ErrMsg2
   CHARACTER(*), PARAMETER      :: RoutineName = 'ReadRAryDefault'


   CALL ReadVar( UnIn, Fil, CharLine, VarName, VarDescr, ErrStat, ErrMsg, UnEc)  !Maybe I should read this in explicitly...

   IF ( PRESENT(IGNORE) ) THEN
      IF ( IGNORE ) THEN
         Def = .TRUE.
         RETURN
      ENDIF

   ENDIF

   CALL Conv2UC( CharLine )

   IF ( TRIM(CharLine) == 'DEFAULT' ) THEN

!      CALL WrScr ( '    A default value will be used for '//TRIM(VarName)//'.' )
      Def = .TRUE.

   ELSE

      IF ( INDEX( CharLine(1:1), 'TF') > 0 ) THEN    ! We don't want 'T' or 'F' read as -1 or 0.
         CALL WrScr1 ( ' Invalid numerical input for "'//TRIM( VarName )//'".' )
      ENDIF

      READ (CharLine,*,IOSTAT=IOS)  RealAry

      IF (IOS /=0) THEN
         RealAry = 0.0_ReKi                        ! set these all to 0
         READ (CharLine,*,IOSTAT=IOS)  RealAry(1)  ! Try reading only the first element
      ENDIF

      CALL CheckIOS ( IOS, Fil, VarName, NumType, ErrStat2, ErrMsg2 )
         call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         
      do i=1,size(RealAry)
         CALL CheckRealVar( RealAry(i), VarName, ErrStat2, ErrMsg2 )
            call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      end do
            
      Def = .FALSE.

   ENDIF


   RETURN

END SUBROUTINE ReadRAryDefault
!=======================================================================
SUBROUTINE ReadRVarDefault ( UnIn, Fil, RealVar, VarName, VarDescr, UnEc, Def, ErrStat, ErrMsg, IGNORE, IGNORESTR )

      ! This routine reads a single real variable from the next line of the input file.
      ! The input is allowed to be "default"

      ! Argument declarations:

   REAL(ReKi), INTENT(INOUT)      :: RealVar                                         ! Real variable being read.
                                  
   INTEGER, INTENT(IN)            :: UnIn                                            ! I/O unit for input file.
   INTEGER, INTENT(IN)            :: UnEc                                            ! I/O unit for echo/summary file.
                                  
   INTEGER(IntKi), INTENT(OUT)    :: ErrStat                                         ! Error status; if present, program does not abort on error
   CHARACTER(*),   INTENT(OUT)    :: ErrMsg                                          ! Error message
                                  
   LOGICAL, INTENT(  OUT)         :: Def                                             ! - on input whether or not to use the default - on output, whether a default was used
   LOGICAL, INTENT(IN),   OPTIONAL:: IGNORE                                          ! whether or not to ignore this input
   LOGICAL, INTENT(INOUT),OPTIONAL:: IGNORESTR                                       ! whether or not user requested to ignore this input

   CHARACTER(250)                 :: CharLine                                        ! Character string being read.
   CHARACTER( *), INTENT(IN)      :: Fil                                             ! Name of the input file.
   CHARACTER( *), INTENT(IN)      :: VarDescr                                        ! Text string describing the variable.
   CHARACTER( *), INTENT(IN)      :: VarName                                         ! Text string containing the variable name.

      ! Local declarations:

   INTEGER                        :: IOS                                             ! I/O status returned from the read statement.
   INTEGER                        :: ErrStat2
   CHARACTER(ErrMsgLen)           :: ErrMsg2
   CHARACTER(*), PARAMETER        :: RoutineName = 'ReadRVarDefault'

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
         Def       = .TRUE.
         RETURN
      ENDIF
   ENDIF

   IF ( TRIM(CharLine) == 'DEFAULT' ) THEN

!      CALL WrScr ( '    A default value will be used for '//TRIM(VarName)//'.' )
      Def = .TRUE.
      RETURN

   ELSE

      IF ( INDEX( CharLine(1:1), 'TF') > 0 ) THEN    ! We don't want 'T' or 'F' read as -1 or 0.
         CALL WrScr1 ( ' Invalid numerical input for "'//TRIM( VarName )//'".' )
      ENDIF

      READ (CharLine,*,IOSTAT=IOS)  RealVar

      CALL CheckIOS ( IOS, Fil, VarName, NumType, ErrStat2, ErrMsg2 )
         call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      CALL CheckRealVar( RealVar, VarName, ErrStat2, ErrMsg2 )
         call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

      Def = .FALSE.
      
      IF ( PRESENT(IGNORESTR) ) IGNORESTR = .FALSE.
      
      
   ENDIF


   RETURN
END SUBROUTINE ReadRVarDefault
!=======================================================================
!> This routine writes the velocity grid to a binary file.
!! The file has a .wnd extension; scaling information is written in a summary
!! file. A tower file with extension .twr is generated if requested, too.
SUBROUTINE WrBinBLADED(p, V, USig, VSig, WSig, ErrStat, ErrMsg)

   IMPLICIT                    NONE

   TYPE(TurbSim_ParameterType),     INTENT(IN)    :: p                             !< TurbSim's parameters
   REAL(ReKi),                      INTENT(IN)    :: V     (:,:,:)                 !< An array containing the summations of the rows of H (NumSteps,NPoints,3).
   REAL(ReKi),                      INTENT(IN)    :: USig                          !< Standard deviation of U
   REAL(ReKi),                      INTENT(IN)    :: VSig                          !< Standard deviation of V
   REAL(ReKi),                      INTENT(IN)    :: WSig                          !< Standard deviation of W
   INTEGER(IntKi),                  intent(  out) :: ErrStat                       !< Error level
   CHARACTER(*),                    intent(  out) :: ErrMsg                        !< Message describing error


   REAL(ReKi)                  :: NewUSig                       ! Value of USig that will be used to scale values in the file
   
   REAL(ReKi)                  :: U_C1                          ! Scale for converting BLADED U data
   REAL(ReKi)                  :: U_C2                          ! Offset for converting BLADED U data
   REAL(ReKi)                  :: V_C                           ! Scale for converting BLADED V data
   REAL(ReKi)                  :: W_C                           ! Scale for converting BLADED W data
   REAL(ReKi)                  :: TI(3)                         ! Turbulence intensity for scaling data
   REAL(ReKi)                  :: TmpU                          ! Max value of |V(:,:,1)-UHub|

   INTEGER(B4Ki)               :: II
   INTEGER(B4Ki)               :: IT
   INTEGER(B4Ki)               :: IY
   INTEGER(B4Ki)               :: IZ

   INTEGER(B4Ki)               :: IP
   INTEGER(B2Ki)               :: TmpVarray(3*p%grid%NumGrid_Y*p%grid%NumGrid_Z) ! This array holds the normalized velocities before being written to the binary file
   INTEGER(B2Ki), ALLOCATABLE  :: TmpTWRarray(:)          ! This array holds the normalized tower velocities   
   
   INTEGER                     :: AllocStat
   INTEGER                     :: UBFFW                                ! I/O unit for BLADED FF data (*.wnd file).
   INTEGER                     :: UATWR                                ! I/O unit for AeroDyn tower data (*.twr file).

   CHARACTER(200)              :: FormStr                              ! String used to store format specifiers.

   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
      
   
      ! We need to take into account the shear across the grid in the sigma calculations for scaling the data, 
      ! and ensure that 32.767*Usig >= |V-UHub| so that we don't get values out of the range of our scaling values
      ! in this BLADED-style binary output.  TmpU is |V-UHub|
   
   TmpU    = MAX( ABS(MAXVAL(V(:,:,1))-p%UHub), ABS(MINVAL(V(:,:,1))-p%UHub) )  !Get the range of wind speed values for scaling in BLADED-format .wnd files         
   NewUSig = MAX(USig,0.05*TmpU)
   
   
      ! Put normalizing factors into the summary file.  The user can use them to
      ! tell a simulation program how to rescale the data.
   TI(1)  = MAX(100.0*Tolerance, NewUSig) / p%UHub
   TI(2)  = MAX(100.0*Tolerance, VSig, 0.05*ABS(MAXVAL(V(:,:,2))), 0.05*ABS(MINVAL(V(:,:,2))) ) / p%UHub  ! put the abs() after the maxval() and minval() to avoid stack-overflow issues with large wind files 
   TI(3)  = MAX(100.0*Tolerance, WSig, 0.05*ABS(MAXVAL(V(:,:,3))), 0.05*ABS(MINVAL(V(:,:,3))) ) / p%UHub  ! put the abs() after the maxval() and minval() to avoid stack-overflow issues with large wind files 

   WRITE (p%US,"(//,'Normalizing Parameters for Binary Data (approximate statistics):',/)")

   FormStr = "(3X,A,' =',F9.4,A)"
   WRITE (p%US,FormStr)  'UBar ', p%UHub,      ' m/s'
   WRITE (p%US,FormStr)  'TI(u)', 100.0*TI(1), ' %'
   WRITE (p%US,FormStr)  'TI(v)', 100.0*TI(2), ' %'
   WRITE (p%US,FormStr)  'TI(w)', 100.0*TI(3), ' %'

   WRITE (p%US,'()')
   WRITE (p%US,FormStr)  'Height Offset', ( p%grid%HubHt - p%grid%GridHeight / 2.0 - p%grid%Zbottom ),  ' m'
   WRITE (p%US,FormStr)  'Grid Base    ', p%grid%Zbottom,                                               ' m'
   
   WRITE (p%US,'()'   )
IF ( p%grid%Periodic ) THEN 
   WRITE (p%US,'( A)' ) 'Creating a PERIODIC output file.'
END IF
   WRITE (p%US,'( A)' ) 'Creating a BLADED LEFT-HAND RULE output file.'
 
      ! Calculate some numbers for normalizing the data.

   U_C1 =  1000.0/( p%UHub*TI(1) )
   U_C2 =  1000.0/TI(1)
   V_C  = -1000.0/( p%UHub*TI(2) )      ! Bladed convention is positive V is pointed along negative Y (IEC turbine coordinate)
   W_C  =  1000.0/( p%UHub*TI(3) )


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


         ! Loop through time.

      DO IT=1,p%grid%NumOutSteps  !Use only the number of timesteps requested originally

            ! Write out grid data in binary form.
         IP = 1
         DO IZ=1,p%grid%NumGrid_Z
            DO IY=1,p%grid%NumGrid_Y

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
                           
         CALL GetNewUnit( UATWR, ErrStat, ErrMsg )
         CALL OpenBOutFile ( UATWR, TRIM( p%RootName )//'.twr', ErrStat, ErrMsg )
         IF (ErrStat >= AbortErrLev) RETURN
         
         
      !IF ( ALLOCATED(p%grid%TwrPtIndx) ) THEN
         ALLOCATE( TmpTwrarray( 3*SIZE(p%grid%TwrPtIndx) ), STAT=AllocStat )
         IF ( AllocStat /= 0 ) THEN
            ErrStat = ErrID_Fatal
            ErrMsg = "WrBinBLADED: Error allocating space for temporary tower output array."
            RETURN
         END IF
      !END IF
         
            
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
   WRITE (UAFFW)  REAL( p%grid%Zbottom            , SiKi )          ! the height of the grid bottom, in m

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
      ErrStat = ErrID_Fatal
      ErrMsg = 'WrBinTURBSIM:Error allocating memory for temporary wind speed array.'
      RETURN
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
SUBROUTINE WrBinHAWC( p, V, USig, VSig, WSig, ErrStat, ErrMsg)
      ! passed variables
   TYPE(TurbSim_ParameterType),     INTENT(IN)     :: p                             !< TurbSim's parameters
   REAL(ReKi),                      INTENT(IN)     :: V     (:,:,:)                 !< THe wind-speed array, with the mean added and the field rotated. (NumSteps,NPoints,3).
   REAL(ReKi),                      INTENT(IN)     :: USig                          !< Standard deviation of U
   REAL(ReKi),                      INTENT(IN)     :: VSig                          !< Standard deviation of V
   REAL(ReKi),                      INTENT(IN)     :: WSig                          !< Standard deviation of W
   INTEGER(IntKi),                  intent(  out)  :: ErrStat                       !< Error level
   CHARACTER(*),                    intent(  out)  :: ErrMsg                        !< Message describing error
   
   ! local variables
   CHARACTER(*),     PARAMETER   :: Comp(3) = (/'u','v','w'/)
   INTEGER(IntKi)                :: IC, IT, IY, IZ, II, IH, IPoint
   INTEGER(IntKi)                :: UnWind


   REAL(ReKi), ALLOCATABLE       :: TmpV(:,:,:)                ! This array holds velocities with the wind direction added (with zero mean)
   REAL(ReKi), PARAMETER         :: FileSigmas(3) = (/ 1.0, 0.8, 0.5 /) ! These are the ratios we will use for each file generated for HAWC2
   REAL(ReKi)                    :: ScaleFactors(3)
   INTEGER(IntKi)                :: ErrStat2
   CHARACTER(ErrMsgLen)          :: ErrMsg2
   CHARACTER(*), PARAMETER       :: RoutineName = 'WrBinHAWC'
   CHARACTER(1024)               :: RootWithoutPathName
   
   ErrStat = ErrID_None
   ErrMsg = ""
   
   CALL GetNewUnit( UnWind, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName) 
      IF (ErrStat >= AbortErrLev) RETURN
   
      ! calculate factors (for workaround in HAWC2); factor = Target / actual
   ScaleFactors = (/ USig, VSig, WSig /)
   do ic=1,3
      if (.not. EqualRealNos( ScaleFactors(ic), 0.0_ReKi ) ) then
         ScaleFactors(ic) = FileSigmas(ic) / ScaleFactors(ic)
      else
         ScaleFactors(ic) = 1.0_ReKi
      end if
   end do
      
   !............................................
   ! write the summary file
   !............................................
   ic = INDEX( p%RootName, '\', BACK=.TRUE. )
   ic = MAX( ic, INDEX( p%RootName, '/', BACK=.TRUE. ) )
   RootWithoutPathName = p%RootName((ic+1):)

   CALL OpenFOutFile ( UnWind, trim(p%RootName)//'.hawc', ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         IF (ErrStat >= AbortErrLev) RETURN

   WRITE (UnWind,"( / '; These full-field turbulence files were generated by ' , A , ' on ' , A , ' at ' , A , '.' )" )  TRIM(GetNVD(TurbSim_Ver)), CurDate(), CurTime()
   
   WRITE( UnWind, '()' ) 
   WRITE( UnWind, '(A)' )  '; the following factor_scaling is required to obtain the original data:'
   WRITE( UnWind, '(A, 3(1x,F15.5) )' )  '; factor_scaling', 1.0_ReKi / ScaleFactors
   WRITE( UnWind, '()' ) 
   WRITE( UnWind, '(A)' )  'begin mann;'
   
   DO IC = 1,3
      WRITE( UnWind, '(2x,A, T30, A, " ;")' ) 'filename_'//Comp(IC), trim(RootWithoutPathName)//'-'//Comp(IC)//'.bin' 
   END DO
   WRITE( UnWind, '(2x,A, T30, I8, 1x, F15.5, " ;")' ) 'box_dim_u', p%grid%NumSteps, p%UHub*p%grid%TimeStep    ! Note: these files have to be periodic, so I'm going to output all of the steps.
   WRITE( UnWind, '(2x,A, T30, I8, 1x, F15.5, " ;")' ) 'box_dim_v', p%grid%NumGrid_Y, p%grid%GridRes_Y
   WRITE( UnWind, '(2x,A, T30, I8, 1x, F15.5, " ;")' ) 'box_dim_w', p%grid%NumGrid_Z, p%grid%GridRes_Z
   WRITE( UnWind, '(A)' )  'end mann;'
   CLOSE ( UnWind )
   
   !............................................
   ! write the binary files for each component
   !............................................
   
      ! temp array to store rotated wind speeds
   CALL AllocAry( TmpV, p%grid%NumSteps, p%grid%NumGrid_Y*p%grid%NumGrid_Z, 3, 'TmpV', ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName) 

   ! Create temp file for easier writing (originally added vertical wind component here, first):
   DO IT = 1,p%grid%NumSteps
   
      IH = 1
      DO IY=p%grid%NumGrid_Y,1,-1
         DO IZ=1,p%grid%NumGrid_Z

            II = (IZ-1)*p%grid%NumGrid_Y + IY
            IPoint = p%grid%GridPtIndx(II)
            
            TmpV(IT,IH,:) = V(IT,IPoint,:)
               ! subtract the mean wind speed from the u-component
            TmpV(IT,IH,1) = TmpV(IT,IH,1) - p%UHub
            
            IH = IH+1

         END DO
      END DO
      
   END DO



   DO IC = 1,3
      CALL WrScr ( ' Generating HAWC2 binary time-series file "'//trim(p%RootName)//'-'//Comp(ic)//'.bin"' )

      CALL OpenBOutFile ( UnWind, trim(p%RootName)//'-'//Comp(ic)//'.bin', ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName) 
         IF (ErrStat >= AbortErrLev) EXIT ! exit this do loop to deallocate array

      DO IT = 1,p%grid%NumSteps

         WRITE( UnWind, IOSTAT=ErrStat2 ) REAL( TmpV(IT,:,IC) * ScaleFactors(IC), SiKi )

      END DO

      CLOSE ( UnWind )
   END DO
   
   IF (ALLOCATED(TmpV)) DEALLOCATE(TmpV)

   END SUBROUTINE WrBinHAWC
!=======================================================================
SUBROUTINE WrFormattedFF(RootName, p_grid, UHub, V )

   IMPLICIT  NONE

   CHARACTER(*),                    INTENT(IN   ) :: RootName             ! Rootname of output file
   TYPE(Grid_ParameterType),        INTENT(IN)    :: p_grid
   REAL(ReKi),                      INTENT(IN)    :: UHub                 ! The steady hub-height velocity
   REAL(ReKi),                      INTENT(IN)    :: V       (:,:,:)      ! The Velocities to write to a file



REAL(ReKi), ALLOCATABLE      :: ZRow       (:)                           ! The horizontal locations of the grid points (NumGrid_Y) at each height.

INTEGER                      :: II
INTEGER                      :: IT
INTEGER                      :: IVec
INTEGER                      :: IY
INTEGER                      :: IZ

CHARACTER(200)               :: FormStr5                                 ! String used to store format specifiers.

INTEGER                      :: UFFF                                     ! I/O unit for formatted FF data.
INTEGER(IntKi)               :: ErrStat
CHARACTER(ErrMsgLen)         :: ErrMsg

   FormStr5 = "(1X,"//trim(num2lstr(max(p_grid%NumGrid_Z,p_grid%NumGrid_Y)))//"(F8.3),:)"

      ! Allocate the array of wind speeds.

      
   CALL GetNewUnit(UFFF)

   DO IVec=1,3

      CALL WrScr ( ' Generating full-field formatted file "'//TRIM(RootName)//'.'//Comp(IVec)//'".' )
      CALL OpenFOutFile ( UFFF, TRIM( RootName )//'.'//Comp(IVec), ErrStat, ErrMsg )
        IF (ErrStat /= ErrID_None) then
           call WrScr(Trim(ErrMsg))
           if (ErrStat >= AbortErrLev) cycle
        end if
        

         ! Create file header.

      WRITE (UFFF,"( / 'This full-field turbulence file was generated by ' , A , ' on ' , A , ' at ' , A , '.' / )" )  TRIM(GetNVD(TurbSim_Ver)), CurDate(), CurTime()

      WRITE (UFFF,"( ' | ', A,'-comp |  Y  x  Z  | Grid Resolution (Y x Z) | Time-step | Hub Elev | Mean U |')")  Comp(IVec)

      WRITE (UFFF,"(I14,I6,F11.3,F11.3,F15.3,F11.2,F10.2)")  p_grid%NumGrid_Y, p_grid%NumGrid_Z, p_grid%GridRes_Y, p_grid%GridRes_Z, p_grid%TimeStep, p_grid%HubHt, UHub
      WRITE (UFFF,"(/,' Z Coordinates (m):')")
      WRITE (UFFF,FormStr5)  ( p_grid%Z(p_grid%GridPtIndx(IZ))-p_grid%HubHt, IZ=1,p_grid%NumGrid_Y*p_grid%NumGrid_Z,p_grid%NumGrid_Y )
      WRITE (UFFF,"(/,' Y Coordinates (m):')")
      WRITE (UFFF,FormStr5)  ( p_grid%Y(p_grid%GridPtIndx(IY)), IY=1,p_grid%NumGrid_Y )

         ! Write out elapsed time & hub-level value before component grid.

      DO IT=1,p_grid%NumOutSteps

         WRITE(UFFF,"(/,1X,2(F8.3))")  p_grid%TimeStep*( IT - 1 ), V(IT,p_grid%HubIndx,IVec)

         DO IZ=1,p_grid%NumGrid_Z  ! From the top to the bottom

            II = ( p_grid%NumGrid_Z - IZ )*p_grid%NumGrid_Y

            WRITE (UFFF,FormStr5)  ( V(IT, p_grid%GridPtIndx(II+IY) ,IVec), IY=1,p_grid%NumGrid_Y ) ! From the left to the right

         ENDDO ! IZ

      ENDDO ! IT

      CLOSE ( UFFF )

   ENDDO ! IVec

      ! Deallocate the array of wind speeds.

   IF ( ALLOCATED( ZRow ) )  DEALLOCATE( ZRow )

END SUBROUTINE WrFormattedFF
!=======================================================================

SUBROUTINE WrSum_UserInput( p_met, p_usr, US  )


   TYPE(Meteorology_ParameterType), INTENT(IN)  :: p_met                       ! meteorology parameters for TurbSim 
   TYPE(UserTSSpec_ParameterType),  INTENT(IN)  :: p_usr                       ! user-defined parameters for TurbSim 
   
   INTEGER, INTENT(IN) :: US
   integer             :: i 


   IF ( p_met%NumUSRz > 0 ) THEN
      WRITE (US,"( // 'User-Defined Profiles:' / )")
   
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
            WRITE (US,"( 1X,F7.2, 2X,F9.2,2X, 3X,F10.2)") p_met%USR_Z(I), p_met%USR_U(I), p_met%USR_WindDir(I)      
         ENDDO   
      ENDIF
      
   ENDIF
   
   IF ( p_usr%nPoints > 0 ) THEN
      WRITE (US,"( // 'Profiles from User-Defined Time-Series Input:' / )")
   
      WRITE (US,"(A)") '  Height   Wind Speed   Horizontal Angle   Vertical Angle'
      WRITE (US,"(A)") '   (m)        (m/s)          (deg)              (deg)    '
      WRITE (US,"(A)") '  ------   ----------   ----------------   --------------'
     
      DO I=p_usr%nPoints,1,-1
         WRITE (US,"( 1X,F7.2, 2X,F9.2,2X, 3X,F10.2,10x,F10.2)") p_usr%pointzi(I), p_usr%meanU(I,1), p_usr%meanDir(I), p_usr%meanVAng(I)    
      END DO   
      
   END IF
   

END SUBROUTINE WrSum_UserInput
!=======================================================================
SUBROUTINE WrSum_SpecModel(p, U, HWindDir, VWindDir, ErrStat, ErrMsg  )

   TYPE(TurbSim_ParameterType),     INTENT(INout) :: p                    ! parameters for TurbSim  !BJJ: FIX THIS!!! TODO: create equivalent plExp elsewhere
   REAL(ReKi),                      INTENT(IN)    :: HWindDir(:)          ! profile of horizontal wind direction
   REAL(ReKi),                      INTENT(IN)    :: VWindDir(:)          ! profile of vertical wind direction
   REAL(ReKi),                      INTENT(IN)    :: U       (:)          ! profile of steady wind speed

   INTEGER(IntKi),                  intent(  out) :: ErrStat              ! Error level
   CHARACTER(*),                    intent(  out) :: ErrMsg               ! Message describing error
   
   
   ! local variables:   

   REAL(ReKi)              ::  HalfRotDiam                           ! Half of the rotor diameter
   
      
   
   REAL(ReKi)              ::  UTmp                                  ! The best fit of observed peak Uh at het height vs jet height
   REAL(ReKi)              ::  U_zb                                  ! The velocity at the bottom of the rotor disk (for estimating log fit)
   REAL(ReKi)              ::  U_zt                                  ! The velocity at the top of the rotor disk (for estimating log fit)
      
   INTEGER                 :: iz, jz                                 ! loop counter/indices of points
   LOGICAL                 ::  HubPr                                 ! Flag to indicate if the hub height is to be printed separately in the summary file

   CHARACTER(200)          :: FormStr                                ! String used to store format specifiers.
   CHARACTER(*),PARAMETER  :: FormStr1 = "('   ',A,' =' ,I9  ,A)"    ! String used to store format specifiers.
   CHARACTER(*),PARAMETER  :: FormStr2 = "('   ',A,' =  ',A)"        ! String used to store format specifiers.
   
   
   ! write to the summary file:
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   
   WRITE (p%US,"( // 'Turbulence Simulation Scaling Parameter Summary:' / )")
   WRITE (p%US,    "('   Turbulence model used                            =  ' , A )")  TRIM(p%met%TMName)

   FormStr  = "('   ',A,' =' ,F9.3,A)"
   
   
   
      ! Write out a parameter summary to the summary file.

IF ( ( p%met%TurbModel_ID  == SpecModel_IECKAI ) .OR. &
     ( p%met%TurbModel_ID  == SpecModel_IECVKM ) .OR. &
     ( p%met%TurbModel_ID  == SpecModel_MODVKM ) .OR. &
     ( p%met%TurbModel_ID  == SpecModel_API    ) )  THEN  ! ADDED BY YGUO on April 19, 2013 snow day!!!
      
      
   IF ( p%IEC%NumTurbInp ) THEN
      WRITE (p%US,FormStr2)      "Turbulence characteristic                       ", "User-specified"
   ELSE
      WRITE (p%US,FormStr2)      "Turbulence characteristic                       ", TRIM(p%IEC%IECTurbE)//p%IEC%IECTurbC
      WRITE (p%US,FormStr2)      "IEC turbulence type                             ", TRIM(p%IEC%IEC_WindDesc)
      
      IF ( p%IEC%IEC_WindType /= IEC_NTM ) THEN       
         WRITE (p%US,FormStr)    "Reference wind speed average over 10 minutes    ", p%IEC%Vref,               " m/s"
         WRITE (p%US,FormStr)    "Annual wind speed average at hub height         ", p%IEC%Vave,               " m/s"
      ENDIF
   ENDIF      
   
   WRITE (p%US,FormStr2)         "IEC standard                                    ", TRIM(p%IEC%IECeditionSTR)
   
   IF ( p%met%TurbModel_ID  /= SpecModel_MODVKM ) THEN
      ! Write out a parameter summary to the summary file.

      WRITE (p%US,FormStr)       "Mean wind speed at hub height                   ", p%UHub,                   " m/s"

      IF (.NOT. p%IEC%NumTurbInp) THEN ! "A", "B", or "C" turbulence:
         IF ( p%IEC%IECedition == 2 ) THEN
            WRITE (p%US,FormStr) "Char value of turbulence intensity at 15 m/s    ", 100.0*p%IEC%TurbInt15,    "%"
            WRITE (p%US,FormStr) "Standard deviation slope                        ", p%IEC%SigmaSlope,         ""
         ELSE                                                                                                 
               ! This is supposed to be the expected value of what is measured at a site.                     
               ! We actually calculate the 90th percentile value to use in the code as the                    
               ! "Characteristic Value".                                                                      
            WRITE (p%US,FormStr) "Expected value of turbulence intensity at 15 m/s", 100.0*p%IEC%TurbInt15,     "%"
         ENDIF                                                                                                
                                                                                                              
      ENDIF                                                                                                   
                                                                                                              
      WRITE (p%US,FormStr)       "Characteristic value of standard deviation      ", p%IEC%SigmaIEC(1),         " m/s"
      WRITE (p%US,FormStr)       "Turbulence scale                                ", p%IEC%Lambda(1),           " m"
                                                                                                              
      IF ( p%met%TurbModel_ID  == SpecModel_IECKAI )  THEN                                                                     
         WRITE (p%US,FormStr)    "u-component integral scale                      ", p%IEC%IntegralScale(1),    " m"
         WRITE (p%US,FormStr)    "Coherency scale                                 ", p%IEC%LC,                  " m"
      ELSEIF ( p%met%TurbModel_ID  == SpecModel_IECVKM )  THEN                                                                 
         WRITE (p%US,FormStr)    "Isotropic integral scale                        ", p%IEC%IntegralScale(1),    " m"
      ENDIF                                                                                                   
                                                                                                              
      WRITE (p%US,FormStr)       "Characteristic value of hub turbulence intensity", 100.0*p%IEC%TurbInt,       "%"
                                                                                                              
   ELSE   ! ModVKM                                                                                                    
!bjj this is never set in TurbSim:       WRITE (p%US,FormStr1)      "Boundary layer depth                            ", NINT(h),                   " m"
      WRITE (p%US,FormStr)       "Site Latitude                                   ", p%met%Latitude,            " degs"
      WRITE (p%US,FormStr)       "Hub mean streamwise velocity                    ", p%UHub,                    " m/s"
      WRITE (p%US,FormStr)       "Hub local u*                                    ", p%met%Ustar,               " m/s" !BONNIE: is this LOCAL? of Disk-avg
      WRITE (p%US,FormStr)       "Target IEC Turbulence Intensity                 ", 100.0*p%IEC%TurbInt,       "%"
      WRITE (p%US,FormStr)       "Target IEC u-component standard deviation       ", p%IEC%SigmaIEC(1),         " m/s"
      WRITE (p%US,FormStr)       "u-component integral scale                      ", p%IEC%Lambda(1),           " m"
      WRITE (p%US,FormStr)       "v-component integral scale                      ", p%IEC%Lambda(2),           " m"
      WRITE (p%US,FormStr)       "w-component integral scale                      ", p%IEC%Lambda(3),           " m"
      WRITE (p%US,FormStr)       "Isotropic integral scale                        ", p%IEC%LC,                  " m"
   ENDIF                                                                                                      
   WRITE (p%US,FormStr)          "Gradient Richardson number                      ", 0.0,                       ""

! p%met%Ustar = SigmaIEC/2.15 ! Value based on equating original Kaimal spectrum with IEC formulation

ELSEIF ( p%met%TurbModel_ID == SpecModel_TIDAL ) THEN
   WRITE (p%US,FormStr2)         "Gradient Richardson number                      ", "N/A"
   WRITE (p%US,FormStr)          "Mean velocity at hub height                     ", p%UHub,                    " m/s"     
   
ELSE   
 
   WRITE (p%US,FormStr)          "Gradient Richardson number                      ", p%met%Rich_No,             ""
   WRITE (p%US,FormStr)          "Monin-Obukhov (M-O) z/L parameter               ", p%met%ZL,                  ""
                                                                                                              
   IF ( .not. EqualRealNos( p%met%ZL, 0.0_ReKi ) ) THEN                                                                                      
      WRITE (p%US,FormStr)       "Monin-Obukhov (M-O) length scale                ", p%met%L,                   " m"
   ELSE                                                                                                       
      WRITE (p%US,FormStr2)      "Monin-Obukhov (M-O) length scale                ", "Infinite"                 
   ENDIF                                                                                                      
   WRITE (p%US,FormStr)          "Mean wind speed at hub height                   ", p%UHub,                    " m/s"     
    
ENDIF !  TurbModel  == 'IECKAI', 'IECVKM', or 'MODVKM'


HalfRotDiam = 0.5*p%grid%RotorDiameter
CALL getVelocity(p, p%UHub,p%grid%HubHt, p%grid%HubHt+HalfRotDiam, U_zt, ErrStat, ErrMsg)   !Velocity at the top of rotor
      
IF ( TRIM(p%met%WindProfileType) /= 'PL' .AND. TRIM(p%met%WindProfileType) /= 'IEC' ) THEN
   CALL getVelocity(p, p%UHub,p%grid%HubHt, p%grid%HubHt            , U_zb, ErrStat, ErrMsg)   !Velocity at hub/reference height
   p%met%PLexp = LOG( U_zt/U_zb ) / LOG( (p%grid%HubHt+HalfRotDiam)/(p%grid%HubHt) )
END IF
CALL getVelocity(p, p%UHub,p%grid%HubHt, p%grid%HubHt-HalfRotDiam, U_zb, ErrStat, ErrMsg)   !Velocity at bottom of the rotor


WRITE(p%US,'()')   ! A BLANK LINE

SELECT CASE ( TRIM(p%met%WindProfileType) )
   CASE ('JET')
      UTmp  = 0.0422*p%met%ZJetMax+10.1979 ! Best fit of observed peak Uh at jet height vs jet height
      
      WRITE (p%US,FormStr2)      "Wind profile type                               ", "Low-level jet"      
      WRITE (p%US,FormStr)       "Jet height                                      ",  p%met%ZJetMax,            " m"
      WRITE (p%US,FormStr)       "Jet wind speed                                  ",  p%met%UJetMax,            " m/s"
      WRITE (p%US,FormStr)       "Upper limit of observed jet wind speed          ",        UTmp,               " m/s"
      WRITE (p%US,FormStr)       "Equivalent power law exponent (upper rotor disk)",  p%met%PLexp,              ""
      
      IF ( UTmp < p%met%UJetMax ) THEN
         CALL SetErrStat( ErrID_Warn, 'The computed jet wind speed is larger than the ' &
                     //'maximum observed jet wind speed at this height.', ErrStat, ErrMsg, 'WrSum_SpecModel')
      ENDIF            
                    
   CASE ('LOG')      
      WRITE (p%US,FormStr2)      "Wind profile type                               ", "Logarithmic"      
      WRITE (p%US,FormStr)       "Equivalent power law exponent (upper rotor disk)",  p%met%PLexp,              ""

   CASE ('H2L')      
      WRITE (p%US,FormStr2)      "Velocity profile type                           ", "Logarithmic (H2L)"      
      WRITE (p%US,FormStr)       "Equivalent power law exponent (upper rotor disk)",  p%met%PLexp,              ""

   CASE ('PL')
      WRITE (p%US,FormStr2)      "Wind profile type                               ", "Power law"      
      WRITE (p%US,FormStr)       "Power law exponent                              ",  p%met%PLExp,              ""
      
   CASE ('USR')
      WRITE (p%US,FormStr2)      "Wind profile type                               ", "Linear interpolation of user-defined profile"
      WRITE (p%US,FormStr)       "Equivalent power law exponent (upper rotor disk)",  p%met%PLexp,              ""
   
   CASE ('TS')
      WRITE (p%US,FormStr2)      "Wind profile type                               ", "Linear interpolation of user-defined profile generated by time-series data"
      WRITE (p%US,FormStr)       "Equivalent power law exponent (upper rotor disk)",  p%met%PLexp,              ""
      
   CASE ('API')
!bjj : fix me:!!! 
      WRITE (p%US,FormStr2)      "Wind profile type                               ", "API"
      WRITE (p%US,FormStr)       "Equivalent power law exponent (upper rotor disk)",  p%met%PLexp,              ""
      
   CASE DEFAULT                
      WRITE (p%US,FormStr2)      "Wind profile type                               ", "Power law on rotor disk, logarithmic elsewhere"
      WRITE (p%US,FormStr)       "Power law exponent                              ",  p%met%PLExp,              ""
      
END SELECT

WRITE(p%US,FormStr)              "Mean shear across rotor disk                    ", (U_zt-U_zb)/p%grid%RotorDiameter, " (m/s)/m"
WRITE(p%US,FormStr)              "Assumed rotor diameter                          ", p%grid%RotorDiameter,      " m"      
WRITE(p%US,FormStr)              "Surface roughness length                        ", p%met%Z0,                  " m"      
WRITE(p%US,'()')                                                                                                 ! A BLANK LINE
WRITE(p%US,FormStr )             "Nyquist frequency of turbulent wind field       ", 0.5_ReKi / p%grid%TimeStep," Hz"
WRITE(p%US,'()')                                                                                                 ! A BLANK LINE
WRITE(p%US,FormStr1)             "Number of time steps in the FFT                 ", p%grid%NumSteps,           ""       
WRITE(p%US,FormStr1)             "Number of time steps output                     ", p%grid%NumOutSteps,        ""          
WRITE(p%US,FormStr1)             "Number of points simulated                      ", p%grid%NPoints,            ""          


IF (p%met%KHtest) THEN
   WRITE(p%US,"(/'KH Billow Test Parameters:' / )") ! HEADER
   WRITE(p%US,FormStr)           "Gradient Richardson number                      ", p%met%Rich_No,           ""
   WRITE(p%US,FormStr)           "Power law exponent                              ", p%met%PLexp,             ""
   WRITE(p%US,FormStr)           "Length of coherent structures                   ", p%grid%UsableTime / 2.0, " s"
   WRITE(p%US,FormStr)           "Minimum coherent TKE                            ", 30.0,                    " (m/s)^2"
ENDIF


   ! Write mean flow angles and wind speed profile to the summary file.

WRITE(p%US,"(//,'Mean Flow Angles:',/)")

FormStr = "(3X,A,F6.1,' degrees')"
WRITE(p%US,FormStr)  'Vertical   =', p%met%VFlowAng
WRITE(p%US,FormStr)  'Horizontal =', p%met%HFlowAng


WRITE(p%US,"(/'Mean Wind Speed Profile:')")

IF ( ALLOCATED( p%met%ZL_profile ) .AND. ALLOCATED( p%met%Ustar_profile ) ) THEN
   WRITE(p%US,"(/,'   Height    Wind Speed   Horizontal Angle  Vertical Angle  U-comp (X)   V-comp (Y)   W-comp (Z)   z/L(z)    u*(z)')")
   WRITE(p%US,"(  '     (m)        (m/s)         (degrees)       (degrees)       (m/s)        (m/s)        (m/s)       (-)      (m/s)')")
   WRITE(p%US,"(  '   ------    ----------   ----------------  --------------  ----------   ----------   ----------   ------   ------')")

   FormStr = '(1X,F8.1,1X,F11.2,2(5x,F11.2),4x,3(2X,F8.2,3X),2(1X,F8.3))'
ELSE
   WRITE(p%US,"(/,'   Height    Wind Speed   Horizontal Angle  Vertical Angle  U-comp (X)   V-comp (Y)   W-comp (Z)')")
   WRITE(p%US,"(  '     (m)        (m/s)         (degrees)       (degrees)       (m/s)        (m/s)        (m/s)   ')")
   WRITE(p%US,"(  '   ------    ----------   ----------------  --------------  ----------   ----------   ----------')")

   FormStr = '(1X,F8.1,1X,F11.2,2(5x,F11.2),4x,3(2X,F8.2,3X))'
ENDIF



   ! Get the angles to rotate the wind components from streamwise orientation to the X-Y-Z grid at the Hub
            
HubPr = .NOT. p%grid%HubOnGrid     !If the hub height is not on the z-grid, print it, too.


   ! Write out the grid points & the hub

DO JZ = p%grid%NumGrid_Z,1, -1
   
   IZ = p%grid%GridPtIndx( (JZ-1)*p%grid%NumGrid_Y+1 )
   
   IF ( HubPr  .AND. ( p%grid%Z(IZ) < p%grid%HubHt ) ) THEN
         
      CALL writeLine( p%grid%HubIndx )      
         
      HubPr = .FALSE. ! we've printed the hub values, so we don't need to check this anymore
   ENDIF
   
   CALL writeLine( IZ )
   
ENDDO ! JZ
   
   ! Write out the tower points   
DO JZ = 2, SIZE(p%grid%TwrPtIndx)
   CALL writeLine( p%grid%TwrPtIndx(JZ) )
ENDDO ! JZ

!..................................................
CONTAINS
   SUBROUTINE writeLine(Indx)

      INTEGER(IntKi), INTENT(IN) ::  Indx
   
      REAL(ReKi)                 ::  CVFA                                  ! Cosine of the vertical flow angle
      REAL(ReKi)                 ::  SVFA                                  ! Sine of the vertical flow angle
      REAL(ReKi)                 ::  CHFA                                  ! Cosine of the horizontal flow angle
      REAL(ReKi)                 ::  SHFA                                  ! Sine of the horizontal flow angle
   
      CHFA = COS( HWindDir(Indx)*D2R )
      SHFA = SIN( HWindDir(Indx)*D2R )
   
      CVFA = COS( VWindDir(Indx)*D2R )
      SVFA = SIN( VWindDir(Indx)*D2R ) 
   
   
      IF ( ALLOCATED( p%met%ZL_profile ) ) THEN
         WRITE(p%US,FormStr)  p%grid%Z(Indx), U(Indx), HWindDir(Indx), VWindDir(Indx), U(Indx)*CHFA*CVFA, U(Indx)*SHFA*CVFA, U(Indx)*SVFA, &
                              p%met%ZL_profile(Indx), p%met%Ustar_profile(Indx)
      ELSE
         WRITE(p%US,FormStr)  p%grid%Z(Indx), U(Indx), HWindDir(Indx), VWindDir(Indx), U(Indx)*CHFA*CVFA, U(Indx)*SHFA*CVFA, U(Indx)*SVFA
      ENDIF   
   END SUBROUTINE writeLine 
   
END SUBROUTINE WrSum_SpecModel
!=======================================================================
SUBROUTINE WrHH_ADtxtfile(p, V, TurbInt, ErrStat, ErrMsg)

   TYPE(TurbSim_ParameterType), INTENT(IN)    :: p                    !< parameters 
   REAL(ReKi),                  INTENT(IN)    :: V     (:,:,:)        !< An array containing the summations of the rows of H (NumSteps,NPoints,3).
   REAL(ReKi),                  INTENT(IN)    :: TurbInt              !< IEC target Turbulence Intensity 
   INTEGER(IntKi),              intent(  out) :: ErrStat              !< Error level
   CHARACTER(*),                intent(  out) :: ErrMsg               !< Message describing error


   REAL(ReKi)              :: V_Inertial(3)                    ! U,V,W components (inertial)
   REAL(ReKi)              :: UH                               ! horizontal wind speed (U+V components)

   REAL(ReKi)              :: Time                             ! The instantaneous Time (s)
   INTEGER(IntKi)          :: IT                               ! loop counter (time step)
   INTEGER                 :: UAHH                             ! I/O unit for AeroDyn HH data (*.hh  file).

   
   

   CALL GetNewUnit( UAHH, ErrStat, ErrMsg)
   CALL OpenFOutFile ( UAHH, TRIM( p%RootName)//'.hh', ErrStat, ErrMsg )
   IF (ErrStat >= AbortErrLev) RETURN

   CALL WrScr ( ' Hub-height AeroDyn data were written to "'//TRIM( p%RootName )//'.hh"' )

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
            
      CALL CalculateWindComponents(V(IT,p%grid%HubIndx,:), p%UHub, p%met%HH_HFlowAng, p%met%HH_VFlowAng, V_Inertial, UH)      
      
      WRITE (UAHH,'(F8.3,3F8.2,3F8.3,F8.2)')  Time, UH, -1.0*R2D*ATAN2( V_Inertial(2) , V_Inertial(1) ), &
                                                    V_Inertial(3), 0.0, p%met%PLExp, 0.0, 0.0 
!bjj: Should we output instantaneous horizontal shear, instead of 0?  
!     Should the power law exponent be an instantaneous value, too?
!                                                   
   END DO

   CLOSE(UAHH)


END SUBROUTINE WrHH_ADtxtfile
!=======================================================================
SUBROUTINE WrHH_binary(p, V, ErrStat, ErrMsg)

   ! Output HH binary turbulence parameters for GenPro analysis.
   ! Output order:  Time,U,Uh,Ut,V,W,u',v',w',u'w',u'v',v'w',TKE,CTKE.


   TYPE(TurbSim_ParameterType),     INTENT(IN)    :: p                               !< parameters 
   REAL(ReKi),                      INTENT(IN)    :: V     (:,:,:)                   !< An array containing the summations of the rows of H (NumSteps,NPoints,3).
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
   
   CALL OpenUOutfile ( UGTP , TRIM( p%RootName)//'.bin', ErrStat, ErrMsg ) 
   IF (ErrStat >= AbortErrLev) RETURN

   CALL WrScr ( ' Hub-height binary turbulence parameters were written to "'//TRIM( p%RootName )//'.bin"' )

   DO IT = 1, p%grid%NumOutSteps
      
      Time = p%grid%TimeStep*( IT - 1 )
      
      CALL CalculateWindComponents(V(IT,p%grid%HubIndx,:), p%UHub, p%met%HH_HFlowAng, p%met%HH_VFlowAng, V_Inertial, UH, UT)
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
SUBROUTINE WrHH_text(p, V, ErrStat, ErrMsg)

   ! Output HH text turbulence parameters
   ! Output order:  Time,U,Uh,Ut,V,W,u',v',w',u'w',u'v',v'w',TKE,CTKE.

   TYPE(TurbSim_ParameterType), INTENT(IN)    :: p             !< parameters 
   REAL(ReKi),                  INTENT(IN)    :: V     (:,:,:) !< An array containing the summations of the rows of H (NumSteps,NPoints,3).
   INTEGER(IntKi),              intent(  out) :: ErrStat       !< Error level
   CHARACTER(*),                intent(  out) :: ErrMsg        !< Message describing error

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
   CALL OpenFOutFile ( UFTP, TRIM( p%RootName)//'.dat', ErrStat, ErrMsg )
   IF (ErrStat >= AbortErrLev) RETURN

   CALL WrScr ( ' Hub-height formatted turbulence parameters were written to "'//TRIM( p%RootName )//'.dat"' )

   WRITE (UFTP,"( / 'This hub-height turbulence-parameter file was generated by ' , A , ' on ' , A , ' at ' , A , '.' / )") &
                TRIM(GetNVD(TurbSim_Ver)), CurDate(), CurTime()
 
   WRITE (UFTP,"('   Time',6X,'U',7X,'Uh',7X,'Ut',8X,'V',8X,'W',8X,'u''',7X,'v''',7X,'w'''," &
                        //"6X,'u''w''',5X,'u''v''',5X,'v''w''',5X,'TKE',6X,'CTKE')")   
   
   
   DO IT = 1, p%grid%NumOutSteps
      
      Time = p%grid%TimeStep*( IT - 1 )
      
      CALL CalculateWindComponents(V(IT,p%grid%HubIndx,:), p%UHub, p%met%HH_HFlowAng, p%met%HH_VFlowAng, V_Inertial, UH, UT)
      CALL CalculateStresses(      V(IT,p%grid%HubIndx,:), uv, uw, vw, TKE, CTKE )
      
                             
      WRITE(UFTP,'(F7.2,13F9.3)') Time,V_Inertial(1),UH,UT,V_Inertial(2),V_Inertial(3), &
                                    V(IT,p%grid%HubIndx,1), V(IT,p%grid%HubIndx,2), V(IT,p%grid%HubIndx,3), &
                                    uw, uv, vw, TKE, CTKE
      
      
   END DO

   CLOSE(UFTP)

END SUBROUTINE WrHH_text
!=======================================================================
SUBROUTINE WrSum_Stats(p, V, USig, VSig, WSig, ErrStat, ErrMsg)


TYPE(TurbSim_ParameterType), INTENT(IN   )   :: p                                ! parameters
REAL(ReKi),                  INTENT(IN   )   :: V           (:,:,:)              ! An array containing the summations of the rows of H (NumSteps,NPoints,3).
REAL(ReKi),                  INTENT(  OUT)   ::  USig                            ! Standard deviation of the u-component wind speed at the hub
REAL(ReKi),                  INTENT(  OUT)   ::  VSig                            ! Standard deviation of the v-component wind speed at the hub
REAL(ReKi),                  INTENT(  OUT)   ::  WSig                            ! Standard deviation of the w-component wind speed at the hub
                             
INTEGER(IntKi),              intent(  out)   :: ErrStat                          ! Error level
CHARACTER(*),                intent(  out)   :: ErrMsg                           ! Message describing error


REAL(DbKi)                   ::  denom                           ! denominator of equation
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
REAL(DbKi)                   ::  UXBar                           ! The mean U-component (u rotated; x-direction) wind speed at the hub
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
REAL(ReKi)                   ::  UXSig                           ! Standard deviation of the U-component wind speed at the hub
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

   CVFA = COS( p%met%HH_VFlowAng*D2R )
   SVFA = SIN( p%met%HH_VFlowAng*D2R )

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
   denom = USig * WSig
   if ( EqualRealNos( denom, 0.0_DbKi ) ) then
      UWcor = 0.0
   else
      UWcor = UW_RS / denom ! this definition assumes u' and w' have zero mean
   end if

   denom = USig * VSig
   if ( EqualRealNos( denom, 0.0_DbKi ) ) then
      UVcor = 0.0
   else
      UVcor = UV_RS / denom
   end if
   
   denom = VSig * WSig
   if ( EqualRealNos( denom, 0.0_DbKi ) ) then
      VWcor = 0.0
   else
      VWcor = VW_RS / denom
   end if
      

      ! Calculate turbulence intensities.
   U_TI = USig/UBar    
   V_TI = VSig/UBar
   W_TI = WSig/UBar

   UH_TI = UHSig/UHBar
   UT_TI = UTSig/UTBar


      ! Write out the hub-level stats to the summary file.

   CALL WrScr ( ' Writing statistics to summary file' )

   WRITE(p%US,"(//,'Hub-Height Simulated Turbulence Statistical Summary:')")
   WRITE(p%US,"(/,3X,'Type of Wind        Min (m/s)   Mean (m/s)    Max (m/s)  Sigma (m/s)       TI (%)')")
   WRITE(p%US,"(  3X,'----------------    ---------   ----------    ---------  -----------       ------')")

   FormStr = "(3X,A,F13.2,2F13.2,2F13.3)"
   !bjj for analysis, extra precision: FormStr = "(3X,A,F13.2,2F13.2,2F13.6)"

   WRITE (p%US,FormStr)  'Longitudinal (u)',  Umin,  UBar,  Umax,  USig, 100.0* U_TI
   WRITE (p%US,FormStr)  'Lateral (v)     ',  Vmin,  VBar,  Vmax,  VSig, 100.0* V_TI
   WRITE (p%US,FormStr)  'Vertical (w)    ',  Wmin,  WBar,  Wmax,  WSig, 100.0* W_TI
   WRITE (p%US,FormStr)  'U component     ', UXmin, UXBar, UXmax, UXSig, 100.0*UXSig/UXBar
   WRITE (p%US,FormStr)  'V component     ', UYmin, UYBar, UYmax, UYSig, 100.0*UYSig/UXBar
   WRITE (p%US,FormStr)  'W component     ', UZmin, UZBar, UZmax, UZSig, 100.0*UZSig/UXBar
   WRITE (p%US,FormStr)  'Horizontal (U&V)', UHmin, UHBar, UHmax, UHSig, 100.0*UH_TI
   WRITE (p%US,FormStr)  'Total           ', UTmin, UTBar, UTmax, UTSig, 100.0*UT_TI

   WRITE(p%US,"(/,3X,'                    Min Reynolds     Mean Reynolds    Max Reynolds    Correlation')")
   WRITE(p%US,"(  3X,'Product             Stress (m/s)^2   Stress (m/s)^2   Stress (m/s)^2  Coefficient')")
   WRITE(p%US,"(  3X,'----------------    --------------   --------------   --------------  -----------')")

   FormStr = "(3X,A,3(3X,F12.3,3X),F11.3)"
   WRITE (p%US,FormStr)  "u'w'            ",  UWMin,  UW_RS, UWMax, UWcor
   WRITE (p%US,FormStr)  "u'v'            ",  UVMin,  UV_RS, UVMax, UVcor
   WRITE (p%US,FormStr)  "v'w'            ",  VWMin,  VW_RS, VWMax, VWcor

   FormStr = "(3X,A,' = ',F10.3,A)"
   WRITE(p%US,"(/)")   ! blank line
   WRITE(p%US,FormStr)  "Friction Velocity (Ustar) ", SUstar,  " m/s"
   WRITE(p%US,FormStr)  "Maximum Instantaneous TKE ", TKEmax,  " (m/s)^2"
   WRITE(p%US,FormStr)  "Maximum Instantaneous CTKE", CTKEmax, " (m/s)^2"

      !  Allocate the array of standard deviations.

   CALL AllocAry( SDary, p%grid%NumGrid_Y, 'SDary (standard deviations)', ErrStat, ErrMsg)
   IF (ErrStat >= AbortErrLev) RETURN


      ! Calculate standard deviations for each grid point.  Write them to summary file.

   WRITE(p%US,"(//,'Grid Point Variance Summary:',/)")
   WRITE(p%US,"(3X,'Y-coord',"//TRIM(Num2LStr(p%grid%NumGrid_Y))//"F8.2)")  p%grid%Y( p%grid%GridPtIndx(1:p%grid%NumGrid_Y) )


   UTmp = 0
   VTmp = 0
   WTmp = 0

   DO IVec=1,3

      WRITE(p%US,"(/,3X,'Height   Standard deviation at grid points for the ',A,' component:')")  Comp(IVec)

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

         WRITE(p%US,"(F9.2,1X,"//TRIM(Num2LStr(p%grid%NumGrid_Y))//"F8.3)") p%grid%Z( p%grid%GridPtIndx( (IZ-1)*p%grid%NumGrid_Y+1 ) ), SDary(1:p%grid%NumGrid_Y)

         IF ( IVec == 1 ) THEN
            UTmp = UTmp + SUM( SDary )
         ELSEIF ( IVec == 2 ) THEN
            VTmp = VTmp + SUM( SDary )
         ELSE
            WTmp = WTmp + SUM( SDary )
         ENDIF
      ENDDO ! IZ

   ENDDO ! Ivec

   
   WRITE(p%US,"(/'   Mean standard deviation across all grid points:')")
   WRITE(p%US,FormStr2) Comp(1), UTmp / ( p%grid%NumGrid_Y*p%grid%NumGrid_Z )
   WRITE(p%US,FormStr2) Comp(2), VTmp / ( p%grid%NumGrid_Y*p%grid%NumGrid_Z ) 
   WRITE(p%US,FormStr2) Comp(3), WTmp / ( p%grid%NumGrid_Y*p%grid%NumGrid_Z ) 


      !  Deallocate the array of standard deviations.

   IF ( ALLOCATED( SDary ) )  DEALLOCATE( SDary )



END SUBROUTINE WrSum_Stats
!=======================================================================
!> Calculate the mean velocity and turbulence intensity of the U-component
!! of the interpolated hub point for comparison with InflowWind output.
SUBROUTINE WrSum_InterpolatedHubStats(p, V)

      ! passed variables:
   TYPE(TurbSim_ParameterType),     INTENT(IN)     ::  p                               !< TurbSim's parameters
   REAL(ReKi),                      INTENT(INOUT)  ::  V(:,:,:)                        !< velocity, aligned along the streamwise direction without mean values added 

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

   WRITE(p%US,"(//,'U-component (X) statistics from the interpolated hub point:',/)")
   WRITE(p%US,"(3X,A,' =',F9.4,A)")  'Mean' , UGridMean, ' m/s'
   WRITE(p%US,"(3X,A,' =',F9.4,A)")  'TI  ' , UGridTI  , ' %'


END SUBROUTINE WrSum_InterpolatedHubStats
!=======================================================================
SUBROUTINE WrSum_EchoInputs(p )

      ! passed variables:
   TYPE(TurbSim_ParameterType),     INTENT(IN)     ::  p                               !< TurbSim's parameters
   
   INTEGER                                         ::  I          ! loop counter
   CHARACTER(10)                                   ::  TmpStr     ! temporary string used to write output to summary file



!..................................................................................................................................
   WRITE (p%US,"( / 'Runtime Options:' / )")
   WRITE (p%US,"( I10 , 2X , 'Random seed #1' )"                            )  p%RNG%RandSeed(1)
   
   IF (p%RNG%pRNG == pRNG_INTRINSIC) THEN
      WRITE (p%US,"( I10 , 2X , 'Random seed #2' )"                         )  p%RNG%RandSeed(2)
   ELSE
      WRITE (p%US,"( 4X, A6, 2X, 'Type of random number generator' )"       )  p%RNG%RNG_type
   ENDIF
      
   WRITE (p%US,"( L10 , 2X , 'Output binary HH turbulence parameters?' )"   )  p%WrFile(FileExt_BIN)
   WRITE (p%US,"( L10 , 2X , 'Output formatted turbulence parameters?' )"   )  p%WrFile(FileExt_DAT)
   WRITE (p%US,"( L10 , 2X , 'Output AeroDyn HH files?' )"                  )  p%WrFile(FileExt_HH)
   WRITE (p%US,"( L10 , 2X , 'Output AeroDyn FF files?' )"                  )  p%WrFile(FileExt_BTS)
   WRITE (p%US,"( L10 , 2X , 'Output BLADED FF files?' )"                   )  p%WrFile(FileExt_WND)
   WRITE (p%US,"( L10 , 2X , 'Output tower data?' )"                        )  p%WrFile(FileExt_TWR)
   WRITE (p%US,"( L10 , 2X , 'Output HAWC FF files?' )"                     )  p%WrFile(FileExt_HAWC)
   WRITE (p%US,"( L10 , 2X , 'Output formatted FF files?' )"                )  p%WrFile(FileExt_UVW)
   WRITE (p%US,"( L10 , 2X , 'Output coherent turbulence time step file?' )")  p%WrFile(FileExt_CTS)
   
   SELECT CASE ( p%IEC%ScaleIEC )
      CASE (0)
         TmpStr= "NONE"
      CASE (1, -1)   ! included the -1 for reading t/f on other systems
         TmpStr = "HUB"
      CASE (2)
         TmpStr = "ALL"
   ENDSELECT   
   
   WRITE (p%US,"( I2, ' - ', A5, 2X , 'IEC turbulence models scaled to exact specified standard deviation' )")  p%IEC%ScaleIEC, TRIM(TmpStr)
   
   
!..................................................................................................................................
   WRITE (p%US,"( // 'Turbine/Model Specifications:' / )")
   WRITE (p%US,"( I10   , 2X , 'Vertical grid-point matrix dimension' )"  )  p%grid%NumGrid_Z
   WRITE (p%US,"( I10   , 2X , 'Horizontal grid-point matrix dimension' )")  p%grid%NumGrid_Y
   WRITE (p%US,"( F10.3 , 2X , 'Time step [seconds]' )"                   )  p%grid%TimeStep
   WRITE (p%US,"( F10.3 , 2X , 'Analysis time [seconds]' )"               )  p%grid%AnalysisTime
   WRITE (p%US,"( F10.3 , 2X , 'Usable output time [seconds]' )"          )  p%grid%UsableTime
   WRITE (p%US,"( F10.3 , 2X , 'Hub height [m]' )"                        )  p%grid%HubHt  
   WRITE (p%US,"( F10.3 , 2X , 'Grid height [m]' )"                       )  p%grid%GridHeight
   WRITE (p%US,"( F10.3 , 2X , 'Grid width [m]' )"                        )  p%grid%GridWidth
   WRITE (p%US,"( F10.3 , 2X , 'Vertical flow angle [degrees]' )"         )  p%met%VFlowAng
   WRITE (p%US,"( F10.3 , 2X , 'Horizontal flow angle [degrees]' )"       )  p%met%HFlowAng
   

!..................................................................................................................................
   WRITE (p%US,"( // 'Meteorological Boundary Conditions:' / )")
   WRITE (p%US, "( 4X , A6 , 2X , '"//TRIM( p%met%TMName )//" spectral model' )")  p%met%TurbModel
   IF (p%IEC%IECstandard > 0) then
      WRITE (p%US,"( 7X, I3, 2X, 'IEC standard: ', A )")  p%IEC%IECstandard, TRIM(p%IEC%IECeditionSTR)
      IF (p%IEC%NumTurbInp) THEN
         WRITE (p%US,"( F10.3 , 2X , 'Percent turbulence intensity, ', A )")  p%IEC%PerTurbInt, TRIM(p%IEC%IECeditionSTR)
      ELSE
         WRITE (p%US,"( 9X , A1 , 2X , 'IEC turbulence characteristic' )"  )  p%IEC%IECTurbC   
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
                        
      WRITE (p%US,"( 4X, A6 , 2X , 'IEC ', A )")  TRIM(p%IEC%IECTurbE)//TRIM(TmpStr), TRIM(p%IEC%IEC_WindDesc)
      
   ELSE
      WRITE (p%US,"( 7X, A3, 2X, 'IEC standard' )"                        )  'N/A'
      IF (p%met%KHtest) THEN
         WRITE (p%US,"( 4X, A6, 2X, 'Kelvin-Helmholtz billow test case' )")  'KHTEST'   
      ELSE   
         WRITE (p%US,"( A10, 2X, 'IEC turbulence characteristic' )"       )  'N/A'   
      END IF      
      WRITE (p%US,"( A10 , 2X , 'IEC turbulence type' )"                  )  'N/A'
      
   END IF

   IF ( p%IEC%IEC_WindType == IEC_ETM ) THEN
      WRITE (p%US,"( F10.3, 2X, 'IEC Extreme Turbulence Model (ETM) ""c"" parameter [m/s]' )") p%IEC%ETMc 
   ELSE
      WRITE (p%US,"(  A10, 2X, 'IEC Extreme Turbulence Model (ETM) ""c"" parameter [m/s]' )")  'N/A'   
   END IF      

   WRITE (p%US,"(   A10 , 2X , 'Wind profile type' )"                  )  p%met%WindProfileType   
   WRITE (p%US,"( F10.3 , 2X , 'Reference height [m]' )"               )  p%met%RefHt  !BJJ: TODO: check if refht makes sense (or is used) for USR profile.
   IF ( p%met%WindProfileType == 'USR' .OR. p%met%WindProfileType == 'TS' ) THEN
      WRITE (p%US,"(  A10, 2X, 'Reference wind speed [m/s]' )"         )  'N/A'   
   ELSE
      WRITE (p%US,"( F10.3 , 2X , 'Reference wind speed [m/s]' )"      )  p%met%URef
   END IF
   
      
   IF ( p%met%WindProfileType == 'JET' ) THEN
      WRITE (p%US,"( F10.3, 2X, 'Jet height [m]' )"                    ) p%met%ZJetMax 
   ELSE
      WRITE (p%US,"(  A10, 2X, 'Jet height [m]' )"                     )  'N/A'   
   END IF      
      
   IF ( INDEX( 'JLUHAT', p%met%WindProfileType(1:1) ) > 0   ) THEN
      WRITE (p%US,"(  A10, 2X, 'Power law exponent' )"                 )  'N/A'   
   ELSE
      WRITE (p%US,"( F10.3 , 2X , 'Power law exponent' )"              )  p%met%PLExp
   END IF      
      
   IF ( p%met%TurbModel_ID==SpecModel_TIDAL  ) THEN
      WRITE (p%US,"(  A10, 2X, 'Surface roughness length [m]' )"       )  'N/A'   
   ELSE      
      WRITE (p%US,"( F10.3 , 2X , 'Surface roughness length [m]' )"    )  p%met%Z0
   END IF
        
!..................................................................................................................................

WRITE (p%US,"( // 'Non-IEC Meteorological Boundary Conditions:' / )")
   
   IF ( p%met%TurbModel_ID /= SpecModel_IECKAI .AND. p%met%TurbModel_ID /= SpecModel_IECVKM .AND. p%met%TurbModel_ID /= SpecModel_API ) THEN
      WRITE (p%US,"( F10.3 , 2X , 'Site latitude [degrees]' )"    )  p%met%Latitude      
   ELSE
      WRITE (p%US,"( A10 , 2X , 'Site latitude [degrees]' )"    )  'N/A'            
   END IF


   IF ( .NOT. p%met%IsIECModel .AND. p%met%TurbModel_ID /= SpecModel_TIDAL   ) THEN
      WRITE (p%US,"( F10.3 , 2X , 'Gradient Richardson number' )")  p%met%Rich_No
   ELSE
      WRITE (p%US,"(   a10 , 2X , 'Gradient Richardson number' )")  'N/A'
   END IF
   
   IF ( .NOT. p%met%IsIECModel  ) THEN
      WRITE (p%US,"( F10.3 , 2X , 'Friction or shear velocity [m/s]' )")  p%met%Ustar
      
      IF (p%met%ZL>=0. .AND. p%met%TurbModel_ID /= SpecModel_GP_LLJ) THEN
         WRITE (p%US,'(   A10 , 2X , "Mixing layer depth [m]" )'       )  'N/A'
      ELSE         
         WRITE (p%US,"( F10.3 , 2X , 'Mixing layer depth [m]' )"       )  p%met%ZI
      END IF

   ELSE
      WRITE (p%US,'(   A10 , 2X , "Friction or shear velocity [m/s]" )')  'N/A'
      WRITE (p%US,'(   A10 , 2X , "Mixing layer depth [m]" )'          )  'N/A'
   END IF

   IF (.NOT. p%met%UWskip) THEN
      WRITE (p%US,'( F10.3 , 2X , "Mean hub u''w'' Reynolds stress" )' )  p%met%PC_UW
   ELSE
      WRITE (p%US,'(   A10 , 2X , "Mean hub u''w'' Reynolds stress" )' )  'N/A'
   END IF
   
   IF (.NOT. p%met%UVskip) THEN
      WRITE (p%US,'( F10.3 , 2X , "Mean hub u''v'' Reynolds stress" )' )  p%met%PC_UV
   ELSE
      WRITE (p%US,'(   A10 , 2X , "Mean hub u''v'' Reynolds stress" )' )  'N/A'
   END IF

   IF (.NOT. p%met%VWskip) THEN
      WRITE (p%US,'( F10.3 , 2X , "Mean hub v''w'' Reynolds stress" )' )  p%met%PC_VW
   ELSE   
      WRITE (p%US,'(   A10 , 2X , "Mean hub v''w'' Reynolds stress" )' )  'N/A'
   END IF
   
!..................................................................................................................................

WRITE (p%US,"( // 'Spatial Coherence Models:' / )")
   
   do i=1,3
      SELECT CASE (p%met%SCMod(i))
         CASE (CohMod_GENERAL)
            TmpStr = "GENERAL"
         CASE (CohMod_IEC)
            TmpStr = "IEC"
         CASE (CohMod_NONE)
            TmpStr = "NONE"
         CASE (CohMod_API)
            TmpStr = "API"
      END SELECT
      WRITE (p%US,'(   A10 , 2X , A, "-component coherence model" )' )  TRIM(TmpStr), Comp(i)
   end do

   do i=1,3
      IF ( p%met%SCMod(i) == CohMod_General .OR. p%met%SCMod(i) == CohMod_IEC ) THEN
         WRITE (p%US,"( '(',F9.3,',',G10.3,')',2X , A,'-component coherence parameters' )")  p%met%InCDec(i), p%met%InCohB(i), Comp(i)
      ELSE
         WRITE (p%US,"( A22,2X , A,'-component coherence parameters' )")  'N/A', Comp(i)         
      END IF
   end do
   
   IF ( ANY(p%met%SCMod == CohMod_General) ) THEN
      WRITE (p%US,'( F10.3 , 2X , "Coherence exponent" )' )  p%met%CohExp
   ELSE
      WRITE (p%US,'( A10   , 2X , "Coherence exponent" )' )  'N/A'
   END IF
   
   
!..................................................................................................................................
!***
IF ( .NOT. p%WrFile(FileExt_CTS) .OR. p%met%IsIECModel ) RETURN  
!***
      WRITE (p%US,"( // 'Coherent Turbulence Scaling Parameters:' / )")
   

      IF ( LEN( TRIM(p%CohStr%CTEventPath) ) <= 10 )  THEN
         WRITE (p%US,"( A10 , 2X , 'Name of the path containing the coherent turbulence data files' )") TRIM(p%CohStr%CTEventPath)
      ELSE
         WRITE (p%US,"( A, /, 12X , 'Name of the path containing the coherent turbulence data files' )") TRIM(p%CohStr%CTEventPath)
      ENDIF
      WRITE (p%US,"( 7X, A3, 2X, 'Type of coherent turbulence data files' )") TRIM(p%CohStr%CText)
!      WRITE (p%US,"( L10 , 2X , 'Randomize the disturbance scale and location?' )")  Randomize
      WRITE (p%US,"( F10.3 , 2X , 'Disturbance scale (ratio of wave height to rotor diameter)' )")        p%CohStr%DistScl
      WRITE (p%US,"( F10.3 , 2X , 'Fractional location of tower centerline from right' )")                p%CohStr%CTLy
      WRITE (p%US,"( F10.3 , 2X , 'Fractional location of hub height from the bottom of the dataset' )")  p%CohStr%CTLz
      WRITE (p%US,"( F10.3 , 2X , 'Minimum start time for coherent structures [seconds]' )")              p%CohStr%CTStartTime
            
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

         CALL CheckRealVar( PerTurbInt, 'IECTurbC/PerTurbInt', ErrStat, ErrMsg )
         
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
!> This subroutine processes the user input from the IEC_WindType line
!! and initializes the variables p%IEC%IECTurbE, p%IEC%Vref, 
!! p%IEC%IEC_WindType, and p%IEC%IEC_WindDesc.
SUBROUTINE ProcessLine_IEC_WindType(Line, p, ErrStat, ErrMsg)
   
   TYPE(TurbSim_ParameterType), INTENT(INOUT) :: p          !< TurbSim parameters      
   CHARACTER(*),                INTENT(INOUT) :: Line       !< on entry, the line from the input file. may be modified in this routine
                                                                
   INTEGER(IntKi),              intent(  out) :: ErrStat    !< Error level
   CHARACTER(*),                intent(  out) :: ErrMsg     !< Message describing error
   

   ErrStat = ErrID_None
   ErrMsg  = ""
   
   
   IF ( p%met%IsIECModel .AND. p%met%TurbModel_ID /= SpecModel_MODVKM  ) THEN

      CALL Conv2UC( Line )

      p%IEC%IECTurbE   = Line(1:1)

      ! Let's see if the first character is a number (for the ETM case)
      SELECT CASE ( p%IEC%IECTurbE )
         CASE ('1')
            p%IEC%Vref = 50.0_ReKi
            Line = Line(2:)
         CASE ('2')
            p%IEC%Vref = 42.5_ReKi
            Line = Line(2:)
         CASE ('3')
            p%IEC%Vref = 37.5_ReKi
            Line = Line(2:)
         CASE DEFAULT
               ! There's no number at the start of the string so let's move on (it's NTM).
            p%IEC%Vref = -999.9_ReKi
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
            CALL SetErrStat( ErrID_Fatal, 'Valid entries for the IEC wind turbulence are "NTM", "xETM", "xEWM1", or "xEWM50", '// &
                             'where x is the wind turbine class (1, 2, or 3).', ErrStat, ErrMsg, 'ProcessLine_IEC_WindType')
      END SELECT

      IF ( p%IEC%IEC_WindType /= IEC_NTM ) THEN

         IF (p%IEC%IECedition /= 3 .OR. p%IEC%IECstandard == 2) THEN
            CALL SetErrStat( ErrID_Fatal, 'The extreme turbulence and extreme wind speed models are available with '// &
                         'the IEC 61400-1 Ed. 3 or 61400-3 scaling only.', ErrStat, ErrMsg, 'ProcessLine_IEC_WindType')            
         ENDIF

         IF (p%IEC%Vref < 0. ) THEN
            CALL SetErrStat( ErrID_Fatal, 'A wind turbine class (1, 2, or 3) must be specified with the '// &
                         'extreme turbulence and extreme wind types. (i.e. "1ETM")', ErrStat, ErrMsg, 'ProcessLine_IEC_WindType')
         ENDIF

         IF ( p%IEC%NumTurbInp ) THEN
            CALL SetErrStat( ErrID_Fatal, 'When the turbulence intensity is entered as a percent, '//&
                             'the IEC wind type must be "NTM".', ErrStat, ErrMsg, 'ProcessLine_IEC_WindType')
         ENDIF

      ELSE

         p%IEC%IECTurbE = ' '

      ENDIF

   ELSE
      p%IEC%IEC_WindType = IEC_NTM
      p%IEC%IEC_WindDesc = 'Normal turbulence'
      p%IEC%IECTurbE     = ' '                   ! unused for non-IEC models
      p%IEC%Vref         = -999.9_ReKi           ! unused for non-IEC models
   ENDIF   
   
   
   
END SUBROUTINE ProcessLine_IEC_WindType
!=======================================================================
SUBROUTINE GetDefaultSCMod( TurbModel_ID, SCMod )

   INTEGER(IntKi), INTENT(IN   )            :: TurbModel_ID      ! turbulence model Identifier
   INTEGER(IntKi), INTENT(  OUT)            :: SCMod(3)          ! default spatial coherence model

   
   
   SELECT CASE (TurbModel_ID)
      CASE ( SpecModel_IECKAI, SpecModel_IECVKM, SpecModel_MODVKM, SpecModel_USRVKM)
         SCMod(1)   = CohMod_IEC
         SCMod(2:3) = CohMod_NONE
      
      CASE ( SpecModel_API )
         SCMod(1)   = CohMod_API
         SCMod(2:3) = CohMod_NONE
      
      CASE ( SpecModel_USER )
         SCMod(1)   = CohMod_GENERAL
         SCMod(2:3) = CohMod_NONE
      
      CASE ( SpecModel_None )
         SCMod      = CohMod_NONE
      
      CASE DEFAULT
         SCMod      = CohMod_GENERAL
      
   END SELECT


END SUBROUTINE GetDefaultSCMod
!=======================================================================
SUBROUTINE GetDefaultCoh(TurbModel_ID, RICH_NO, WS, Ht, InCDec, InCohB )
! This routine should NOT be called after CalcIECScalingParams() because it will
! incorrectly overwrite the InCDec and InCohB parameters for the IEC models.

   ! These numbers come from Neil's analysis

   INTEGER(IntKi), INTENT(IN)               :: TurbModel_ID      ! turbulence model Identifier
   REAL(ReKi),     INTENT(IN)               :: RICH_NO           ! Richardson Number (stability)
   REAL(ReKi),     INTENT(IN)               :: Ht                !Height, usually hub height
   REAL(ReKi),     INTENT(IN)               :: WS                !Wind speed, usually = UHub
   REAL(ReKi),     INTENT(  OUT)            :: InCDec(3)         ! default coherence decrement
   REAL(ReKi),     INTENT(  OUT)            :: InCohB(3)         ! default coherence parameter B
   
   
!   REAL(ReKi), PARAMETER                :: a =  0.007697495  !coeffs for WF_xxD best-fit equations
!   REAL(ReKi), PARAMETER                :: b =  0.451759656  !coeffs for WF_xxD best-fit equations
!   REAL(ReKi), PARAMETER                :: c =  6.559106387  !coeffs for WF_xxD best-fit equations
!   REAL(ReKi), PARAMETER                :: d = -0.10471942   !coeffs for WF_xxD best-fit equations
!   REAL(ReKi), PARAMETER                :: e = -1.19488521   !coeffs for WF_xxD best-fit equations
!   REAL(ReKi), PARAMETER                :: f =  0.005529328  !coeffs for WF_xxD best-fit equations
!   REAL(ReKi), PARAMETER                :: g =  0.059157163  !coeffs for WF_xxD best-fit equations


   REAL(ReKi)                           :: Coeffs(10,3)      ! coeffs for WS category coherence decrements
   REAL(ReKi)                           :: Ht1               !Height, set to bounds of the individual models
   REAL(ReKi)                           :: Ht2               !Height squared
   REAL(ReKi)                           :: Ht3               !Height cubed
   REAL(ReKi)                           :: WS1               !Wind speed, set to bounds of individual models
   REAL(ReKi)                           :: RI1               !RICH_NO, set to bounds of individual models
   REAL(ReKi)                           :: RI2               !RICH_NO squared
   REAL(ReKi)                           :: RI3               !RICH_NO  cubed

   INTEGER                              :: I
   INTEGER                              :: Ri_Cat


      IF (RICH_NO <= 0.00_ReKi ) THEN
         IF ( RICH_NO <= - 1.0_ReKi ) THEN
            Ri_Cat = 1
         ELSE
            Ri_Cat = 2
         ENDIF
      ELSEIF ( RICH_NO <= 0.25_ReKi ) THEN
         IF ( RICH_NO <= 0.10_ReKi ) THEN
            Ri_Cat = 3
         ELSE
            Ri_Cat = 4
         ENDIF
      ELSE
            Ri_Cat = 5
      ENDIF

      SELECT CASE (  TurbModel_ID )

         CASE ( SpecModel_GP_LLJ )
            HT1 = MAX( 60.0_ReKi, MIN( Ht, 100.0_ReKi ) )
            IF ( WS <= 14.0 ) THEN
               IF ( WS <= 8.0 ) THEN
                  IF     ( WS <= 6.0 ) THEN
                     coeffs(:,3) = (/  3.1322E+00,  2.2819E-03,  2.9214E+00, -5.2203E-04,  1.1877E+00, &
                                      -5.7605E-02,  3.7233E-06, -3.5021E-01, -1.7555E-03,  3.9712E-04 /)    !W  5
                     IF  ( WS <= 4.0 ) THEN !      WS <=  4
                        RI1 = MAX( 0.0_ReKi, MIN(  RICH_NO, 1.0_ReKi ) )
                        coeffs(:,1) = (/  4.8350E+00, -4.0113E-02,  7.8134E+00, -2.0069E-05, -1.9518E-01, &
                                         -1.4009E-01,  2.3195E-06,  8.2029E-02, -7.4979E-04,  6.1186E-04 /) !U  3
                        coeffs(:,2) = (/  3.2587E+00, -5.9086E-02,  9.7426E+00,  5.7360E-04,  2.1274E-01, &
                                         -1.6398E-01, -8.3786E-07,  6.6896E-02, -3.5254E-03,  6.4833E-04 /) !V  3
                     ELSE                   !  4 < WS <=  6
                        RI1 = MAX( -0.5_ReKi, MIN(  RICH_NO, 1.0_ReKi ) )
                        coeffs(:,1) = (/  9.2474E+00, -4.9849E-02,  6.0887E+00, -5.9124E-04,  4.4312E-02, &
                                         -1.1966E-01,  5.2652E-06, -1.0373E-01,  4.0480E-03,  5.5761E-04 /) !U  5
                        coeffs(:,2) = (/  3.6355E+00,  1.7701E-02,  4.2165E+00, -5.8828E-04,  9.5592E-02, &
                                         -6.5313E-02,  3.3875E-06, -1.7981E-02, -1.6375E-03,  3.0423E-04 /) !V  5
                     ENDIF
                  ELSE                      ! 6  < WS <=  8
                     RI1 = MAX( -0.5_ReKi, MIN(  RICH_NO, 1.0_ReKi ) )
                     coeffs(:,1) = (/  1.1795E+01, -7.5393E-02,  9.5279E+00, -3.4922E-04, -5.8973E-01, &
                                      -1.6753E-01,  4.4267E-06,  2.1797E-01,  7.7887E-04,  7.4912E-04 /)    !U  7
                     coeffs(:,2) = (/  1.7730E+00,  9.6577E-02,  8.1310E+00, -1.2028E-03,  3.0145E-02, &
                                      -1.2282E-01,  4.6866E-06,  3.5748E-02, -2.9013E-03,  4.8368E-04 /)    !V  7
                     coeffs(:,3) = (/  9.1695E-01,  9.1488E-02,  6.7163E+00, -1.2938E-03,  1.0315E+00, &
                                      -1.1976E-01,  5.6039E-06, -2.0416E-01, -3.4698E-03,  6.0175E-04 /)    !W  7
                  ENDIF
               ELSE ! 8.0 < WS <= 14.0
                  IF     (WS <= 10.0) THEN  !  8 < WS <= 10
                     RI1 = MAX( -0.5_ReKi, MIN(  RICH_NO, 1.0_ReKi ) )
                     coeffs(:,1) = (/  8.4674E+00,  1.2922E-01,  8.6170E+00, -3.3048E-03, -3.1928E-02, &
                                      -1.2515E-01,  1.8209E-05,  2.9087E-01, -9.3031E-03,  5.0706E-04 /)    !U  9
                     coeffs(:,2) = (/  2.8145E+00,  1.0257E-01,  4.2987E+00, -1.4901E-03,  4.9698E-02, &
                                      -3.9964E-02,  6.7640E-06,  2.2980E-01, -1.0046E-02,  1.3037E-04 /)    !V  9
                     coeffs(:,3) = (/  2.4952E+00,  5.8000E-02,  1.9851E+00, -9.4027E-04, -4.0135E-02, &
                                      -1.8377E-02,  4.3320E-06, -1.0441E-01,  3.6831E-03,  8.6637E-05 /)    !W  9
                  ELSEIF (WS <= 12.0) THEN  ! 10 < WS <= 12
                     RI1 = MAX( -0.5_ReKi, MIN(  RICH_NO, 1.0_ReKi ) )
                     coeffs(:,1) = (/  1.2473E+01,  3.2270E-02,  1.4508E+01, -2.2856E-03, -1.4652E+00, &
                                      -2.4114E-01,  1.4919E-05,  5.5578E-01, -8.5528E-04,  1.0273E-03 /)    !U  11
                     coeffs(:,2) = (/  1.0882E+00,  1.9425E-01,  8.1533E+00, -2.5574E-03,  4.3113E-01, &
                                      -8.0465E-02,  1.0478E-05,  1.1640E-01, -1.1717E-02,  1.6476E-04 /)    !V  11
                     coeffs(:,3) = (/  5.0280E-01,  1.1637E-01,  4.0130E+00, -1.2034E-03, -2.7592E-01, &
                                      -3.8744E-02,  3.4213E-06, -1.5144E-02,  2.4042E-03,  4.7818E-05 /)    !W  11
                  ELSE                      ! 12 < WS <= 14.0
                     RI1 = MAX( -1.0_ReKi, MIN(  RICH_NO, 1.0_ReKi ) )
                     coeffs(:,1) = (/  8.6311E+00,  2.5614E-01,  1.1165E+01, -5.1685E-03,  3.0895E+00, &
                                      -1.9190E-01,  2.7162E-05, -2.6513E-01, -3.6479E-02,  8.8431E-04 /)    !U  13
                     coeffs(:,2) = (/  1.2842E+00,  2.4007E-01,  5.3653E+00, -3.2589E-03,  3.4715E+00, &
                                      -6.8865E-02,  1.3756E-05, -4.8465E-01, -4.0608E-02,  3.8578E-04 /)    !V  13
                     coeffs(:,3) = (/  4.3681E+00,  1.2251E-02,  1.3826E+00, -1.1592E-04,  3.3654E+00, &
                                      -5.2367E-02, -4.4086E-08, -3.5254E-01, -1.6780E-02,  3.9048E-04 /)    !W  13
                  ENDIF
               ENDIF
            ELSE ! WS > 14
               IF (WS <= 20.0 ) THEN
                  IF     (WS <= 16.0) THEN  ! 14 < WS <= 16
                     RI1 = MAX( -1.0_ReKi, MIN(  RICH_NO, 1.0_ReKi ) )
                     coeffs(:,1) = (/  1.3972E-01,  6.3486E-01,  1.7576E+01, -1.0017E-02,  2.8458E+00, &
                                      -2.5233E-01,  4.6539E-05, -1.8899E-01, -2.6717E-02,  9.5173E-04 /)    !U  15
                     coeffs(:,2) = (/ -7.1243E+00,  5.6768E-01,  1.2886E+01, -7.3277E-03,  3.7880E+00, &
                                      -1.4733E-01,  3.0898E-05, -1.5056E-01, -2.9500E-02,  3.6703E-04 /)    !V  15
                     coeffs(:,3) = (/ -1.1004E+01,  5.3470E-01,  5.3118E+00, -5.8999E-03,  1.9009E+00, &
                                      -2.4063E-02,  2.1755E-05, -4.5798E-01,  1.6885E-02, -3.9974E-04 /)    !W  15
                  ELSEIF (WS <= 18.0) THEN  ! 16 < WS <= 18
                     RI1 = MAX( -0.5_ReKi, MIN(  RICH_NO, 1.0_ReKi ) )
                     coeffs(:,1) = (/ -6.9650E+00,  8.8636E-01,  2.3467E+01, -1.1973E-02, -4.3750E+00, &
                                      -3.5519E-01,  5.0414E-05,  9.1789E-01,  9.8340E-03,  1.5885E-03 /)    !U  17
                     coeffs(:,2) = (/  5.5495E-03,  3.2906E-01,  1.4609E+01, -4.1635E-03, -2.1246E+00, &
                                      -1.8887E-01,  1.6964E-05,  3.7805E-01,  1.1880E-03,  8.8265E-04 /)    !V  17
                     coeffs(:,3) = (/ -1.3195E+00,  2.0022E-01,  2.3490E+00, -2.1308E-03,  3.5582E+00, &
                                       1.4379E-02,  7.6830E-06, -7.6155E-01, -2.4660E-02, -2.0199E-04 /)    !W  17
                  ELSE                      ! 18 < WS <= 20
                     RI1 = MAX( -0.5_ReKi, MIN(  RICH_NO, 1.0_ReKi ) )
                     coeffs(:,1) = (/ -1.3985E+01,  1.3161E+00,  3.4773E+01, -1.9237E-02, -1.9845E+00, &
                                      -5.5817E-01,  8.8310E-05,  1.7142E+00, -4.2907E-02,  2.3932E-03 /)    !U  19
                     coeffs(:,2) = (/ -1.2400E+01,  8.6854E-01,  1.9923E+01, -1.1557E-02, -1.0441E+00, &
                                      -2.4593E-01,  4.9813E-05,  2.7861E-01, -8.6189E-03,  9.4314E-04 /)    !V  19
                     coeffs(:,3) = (/ -9.3436E+00,  6.4950E-01,  1.5316E+01, -8.7208E-03,  1.7329E+00, &
                                      -2.2411E-01,  3.6288E-05, -8.0006E-01, -2.6439E-03,  7.9293E-04 /)    !W  19
                  ENDIF
               ELSE ! WS > 20
                  IF     (WS <= 22.0) THEN  ! 20 < WS <= 22
                     RI1 = MAX( -0.5_ReKi, MIN(  RICH_NO, 1.0_ReKi ) )
                     coeffs(:,1) = (/ -2.4317E+01,  1.8176E+00,  5.3359E+01, -2.5973E-02,  6.0349E+00, &
                                      -7.9927E-01,  1.1558E-04,  1.5926E+00, -1.5005E-01,  3.1688E-03 /)    !U  21
                     coeffs(:,2) = (/  8.0459E+00,  1.8058E-01,  1.9426E+01, -3.6730E-03, -9.9717E-01, &
                                      -1.8249E-01,  1.9237E-05,  4.9173E-01, -1.8255E-02,  6.9371E-04 /)    !V  21
                     coeffs(:,3) = (/ -2.3544E+01,  1.1403E+00,  8.3526E+00, -1.4511E-02,  7.2014E+00, &
                                       5.0216E-02,  5.9947E-05, -1.0659E+00, -7.4769E-02, -9.8390E-04 /)    !W  21
                  ELSEIF (WS <= 24.0) THEN  ! 22 < WS <= 24
                     RI1 = MAX( 0.0_ReKi, MIN(  RICH_NO, 1.0_ReKi ) )
                     coeffs(:,1) = (/ -3.5790E+01,  1.5374E+00,  1.1322E+02, -1.6884E-02, -1.7767E+01, &
                                      -1.8122E+00,  6.8247E-05,  7.2101E+00,  3.5536E-02,  7.9269E-03 /)    !U  23
                     coeffs(:,2) = (/ -7.2883E+01,  2.8210E+00,  8.6392E+01, -3.1084E-02, -2.4938E+01, &
                                      -1.5898E+00,  1.0997E-04,  7.1972E+00,  1.2624E-01,  9.3084E-03 /)    !V  23
                     coeffs(:,3) = (/ -3.2844E+01,  1.2683E+00,  3.2032E+01, -1.3197E-02, -1.1129E+01, &
                                      -3.6741E-01,  4.2852E-05,  4.1336E+00,  2.4775E-02,  1.8431E-03 /)    !W  23
                  ELSE                      ! 24 < WS
                     RI1 = MAX( -0.5_ReKi, MIN(  RICH_NO, 1.0_ReKi ) )
                     coeffs(:,1) = (/  2.2906E+01,  9.3209E-02,  1.5448E+01, -5.7421E-03, -8.9114E+00, &
                                      -3.1547E-02,  4.0144E-05,  5.4544E-01,  5.3557E-02, -3.1299E-04 /)    !U  25
                     coeffs(:,2) = (/ -1.1903E+01,  1.1104E+00,  1.7962E+01, -1.6045E-02, -9.2458E+00, &
                                      -4.4526E-02,  6.9880E-05,  2.8017E+00, -2.7211E-02, -8.4099E-04 /)    !V  25
                     coeffs(:,3) = (/  6.1054E-01,  7.1841E-03,  4.2996E+00,  2.9071E-04, -2.0002E+00, &
                                      -7.0403E-02, -2.8931E-06,  2.3943E-02,  1.8395E-02,  5.0406E-04 /)    !W  25
                  ENDIF
               ENDIF
            ENDIF


            HT2 = HT1*HT1
            HT3 = HT1*HT2
            RI2 = RI1*RI1
            RI3 = RI1*RI2

            DO I = 1,3
               InCDec(I) = coeffs( 1,I) + coeffs(2,I)*Ht1 + coeffs(3,I)*RI1 &
                                        + coeffs(4,I)*Ht2 + coeffs(5,I)*RI2 + coeffs( 6,I)*Ht1*RI1 &
                                        + coeffs(7,I)*Ht3 + coeffs(8,I)*RI3 + coeffs( 9,I)*Ht1*RI2 &
                                                                            + coeffs(10,I)*Ht2*RI1
            ENDDO

            WS1 = MAX(  2.0_ReKi, WS )
            SELECT CASE ( Ri_Cat )
               CASE ( 1, 2)
!                 InCDec   = (/            1.744591004*WS1**0.593219225, &
!                              -0.58750092+1.937230512*WS1**0.400548383, &
!                              -0.57833219+1.450654739*WS1**0.443191083 /)
                  InCohB   = (/-0.00014115+0.006826264/WS1, &
                                           0.014025749/WS1, &
                               0.000480386+0.020982336/WS1 /)

               CASE ( 3 )
!                 InCDec   = (/            1.962126171*WS1**0.575523536, &
!                              -2.79495117+3.698342796*WS1**0.305415750, &
!                                          0.887573173*WS1**0.498317195 /)
                  InCohB   = (/-0.00016838+0.009764148/WS1, &
                                           0.018582932/WS1, &
                               0.001865953+0.061952454/WS1 /)

               CASE ( 4 )
!                 InCDec   = (/            0.817085986*WS1**1.045777184, &
!                                          0.599696362*WS1**1.038373995, &
!                                          1.327586050*WS1**0.590370871 /)
                  InCohB   = (/0.000175033+0.004195814/WS1, &
                                           0.008479460/WS1, &
                               0.002318082+0.027820652/WS1 /)

               CASE ( 5 )
!                 InCDec   = (/            0.959999473*WS1**0.972466847, &
!                              0.082701643+0.867230846*WS1**0.925895412, &
!                                          1.524380209*WS1**0.548060899 /)
                  InCohB   = (/0.000241808+0.004267702/WS1, &
                                           0.005408592/WS1, &
                               0.001150319+0.010744459/WS1 /)
               END SELECT


         CASE ( SpecModel_NWTCUP, SpecModel_USRVKM )
            HT1 = MAX( 25.0_ReKi, MIN( Ht, 50.0_ReKi ) )

            IF ( WS <= 14.0 ) THEN
               RI1 = MAX( -1.0_ReKi, MIN(  RICH_NO, 1.0_ReKi ) )
               IF ( WS <= 8.0 ) THEN
                  IF     (WS <= 4.0 ) THEN  !      WS <=  4
                     coeffs(:,1) = (/  8.1767E+00, -3.1018E-01,  3.3055E-01,  4.4232E-03,  4.6550E-01, &
                                      -2.4582E-02, -5.8568E-06, -8.7873E-02,  1.3070E-02,  3.1871E-04 /)   !U  3
                     coeffs(:,2) = (/  5.8003E+00, -2.0838E-01,  2.8727E-01,  2.8669E-03,  6.9669E-01, &
                                      -8.2249E-03, -2.4732E-06, -1.0826E-01,  9.9973E-03,  1.8546E-05 /)   !V  3
                     coeffs(:,3) = (/  5.9625E+00, -2.9247E-01, -9.3269E-01,  4.4089E-03,  1.3779E-01, &
                                       2.6993E-02, -6.1784E-06, -7.2920E-02,  1.7028E-02, -3.3753E-04 /)   !W  3
                  ELSEIF (WS <= 6.0 ) THEN  !  4 < WS <=  6
                     coeffs(:,1) = (/  1.2891E+01, -4.8265E-01,  3.5549E+00,  6.6099E-03,  8.2275E-01, &
                                      -1.5913E-01, -7.9740E-06, -1.2357E-02,  3.2084E-03,  1.7145E-03 /)   !U  5
                     coeffs(:,2) = (/  8.0267E+00, -2.5275E-01,  1.3801E+00,  3.2447E-03,  1.6004E+00, &
                                      -3.2592E-02, -5.1265E-06, -9.8552E-02, -1.3513E-02,  2.8075E-04 /)   !V  5
                     coeffs(:,3) = (/  7.9593E+00, -3.6336E-01,  1.4974E+00,  5.4012E-03,  9.5041E-01, &
                                      -1.0152E-01, -1.0865E-05,  4.3121E-02, -3.2447E-03,  1.3797E-03 /)   !W  5
                  ELSE                      ! 6  < WS <=  8
                     coeffs(:,1) = (/  1.3702E+01, -4.4674E-01,  3.7943E+00,  5.9350E-03,  9.6026E-01, &
                                      -1.7425E-01, -7.2917E-06, -8.8426E-02,  5.1530E-03,  2.0554E-03 /)   !U  7
                     coeffs(:,2) = (/  9.2471E+00, -2.6247E-01,  1.4504E+00,  3.2436E-03,  1.8823E+00, &
                                      -3.2180E-02, -5.9491E-06, -2.0100E-01, -1.7619E-02,  3.8519E-04 /)   !V  7
                     coeffs(:,3) = (/  8.9439E+00, -3.8885E-01,  2.2175E+00,  5.6207E-03,  7.6040E-01, &
                                      -1.3502E-01, -9.2514E-06,  1.9269E-02,  3.8862E-03,  1.7674E-03 /)   !W  7
                  ENDIF
               ELSE ! 8.0 < WS <= 14.0
                  IF     (WS <= 10.0) THEN  !  8 < WS <= 10
                     coeffs(:,1) = (/  1.9061E+01, -4.5354E-01,  7.5961E+00,  5.2422E-03,  1.5158E+00, &
                                      -2.4908E-01, -2.5277E-06, -1.6660E-01,  1.1369E-02,  3.0156E-03 /)   !U  9
                     coeffs(:,2) = (/  1.3362E+01, -3.3806E-01,  7.0401E+00,  4.5349E-03,  2.6798E+00, &
                                      -2.3637E-01, -9.9075E-06, -2.2373E-01, -1.6644E-03,  2.3879E-03 /)   !V  9
                     coeffs(:,3) = (/  8.8401E+00, -2.9945E-01,  3.7883E+00,  4.4581E-03,  2.0417E+00, &
                                      -2.7852E-01, -7.0750E-06, -6.2618E-02,  1.4646E-02,  3.8512E-03 /)   !W  9
                  ELSEIF (WS <= 12.0) THEN  ! 10 < WS <= 12
                     coeffs(:,1) = (/  3.4011E+01, -1.2590E+00,  1.6320E+01,  1.9225E-02,  6.8346E+00, &
                                      -8.8950E-01, -6.2453E-05, -2.4945E-01, -4.3892E-02,  1.2078E-02 /)   !U  11
                     coeffs(:,2) = (/  1.7135E+01, -4.0754E-01,  1.0282E+01,  5.7832E-03,  6.3056E+00, &
                                      -2.8536E-01, -3.0216E-05, -5.3170E-01, -5.7090E-02,  2.8463E-03 /)   !V  11
                     coeffs(:,3) = (/  1.3002E+01, -4.8326E-01,  3.2819E+00,  7.8800E-03,  2.7094E+00, &
                                      -2.5714E-01, -3.0117E-05, -2.1404E-01, -4.2711E-03,  4.1067E-03 /)   !W  11
                  ELSE                      ! 12 < WS <= 14
                     coeffs(:,1) = (/  2.6682E+01, -9.7229E-01,  1.3191E+01,  1.7604E-02, -1.3537E+00, &
                                      -6.4082E-01, -7.8242E-05,  1.7548E-01,  9.7417E-02,  1.0259E-02 /)   !U  13
                     coeffs(:,2) = (/  1.7083E+01, -4.7346E-01,  1.3515E+01,  7.7832E-03,  5.8633E-01, &
                                      -6.1815E-01, -3.3752E-05, -1.7300E-01,  4.3584E-02,  8.9289E-03 /)   !V  13
                     coeffs(:,3) = (/  1.6015E+01, -6.3912E-01,  1.3137E+01,  9.4757E-03,  2.5549E+00, &
                                      -8.1438E-01, -1.5565E-05,  2.9244E-02,  2.2779E-02,  1.1982E-02 /)   !W  13
                  ENDIF
               ENDIF
            ELSE ! WS > 14
               IF (WS <= 20.0 ) THEN
                  IF     (WS <= 16.0) THEN  ! 14 < WS <= 16
                     RI1 = MAX( -1.0_ReKi, MIN(  RICH_NO, 1.0_ReKi ) )
                     coeffs(:,1) = (/  2.9459E+01, -7.3181E-01,  9.4613E+00,  9.2172E-03,  6.1086E+00, &
                                      -4.9990E-01, -2.9994E-05, -6.9606E-01, -8.5076E-03,  8.1330E-03 /)   !U  15
                     coeffs(:,2) = (/  1.7540E+01, -2.6071E-01,  9.3639E+00,  1.3341E-03,  9.4294E+00, &
                                      -4.2565E-01, -2.7836E-06, -6.7708E-01, -6.9127E-02,  6.2290E-03 /)   !V  15
                     coeffs(:,3) = (/  1.2792E+01, -4.6469E-01,  4.6350E+00,  1.0633E-02,  1.8523E+00, &
                                      -3.2417E-01, -8.5038E-05, -2.2253E-01, -7.3351E-04,  5.4781E-03 /)   !W  15
                  ELSEIF (WS <= 18.0) THEN  ! 16 < WS <= 18
                     RI1 = MAX( -1.0_ReKi, MIN(  RICH_NO, 1.0_ReKi ) )
                     coeffs(:,1) = (/  1.7775E+01,  4.5287E-01,  1.6417E+01, -2.3724E-02,  5.8998E+00, &
                                      -5.3502E-01,  2.6202E-04, -9.9466E-02,  4.1386E-02,  4.5663E-03 /)   !U  17
                     coeffs(:,2) = (/  1.2022E+01,  2.4246E-01,  1.3875E+01, -1.1725E-02,  5.1917E+00, &
                                      -5.4329E-01,  1.1893E-04, -2.0308E-01,  6.5256E-02,  5.6597E-03 /)   !V  17
                     coeffs(:,3) = (/  1.2680E+01, -1.4768E-01,  7.1498E+00, -3.0341E-03,  1.9747E+00, &
                                      -3.8374E-01,  7.0412E-05,  2.2297E-01,  5.9943E-02,  5.3514E-03 /)   !W  17
                  ELSE                      ! 18 < WS <= 20
                     RI1 = MAX( -0.5_ReKi, MIN(  RICH_NO, 1.0_ReKi ) )
                     coeffs(:,1) = (/  3.1187E+01, -6.8540E-01,  7.1288E+00,  1.1923E-02,  8.8547E+00, &
                                       6.3133E-02, -9.4673E-05, -2.5710E+00, -5.4077E-02, -1.2797E-04 /)   !U  19
                     coeffs(:,2) = (/  1.2664E+01,  9.1858E-02,  1.9050E+01, -2.8868E-03,  7.2969E+00, &
                                      -4.4573E-01, -6.1033E-06, -2.0960E+00, -1.9913E-02,  4.9023E-03 /)   !V  19
                     coeffs(:,3) = (/  2.2146E+01, -7.6940E-01,  1.1948E+01,  1.0400E-02,  5.0034E+00, &
                                      -4.3958E-01, -2.5936E-05, -3.0848E-01, -6.3381E-02,  5.1204E-03 /)   !W  19
                  ENDIF
               ELSE ! WS > 20
                  RI1 = MAX( -0.5_ReKi, MIN(  RICH_NO, 1.0_ReKi ) )
                  IF     (WS <= 22.0) THEN  ! 20 < WS <= 22
                     coeffs(:,1) = (/  2.5165E+01, -7.7660E-02,  1.9692E+01, -1.1794E-02,  9.8635E+00, &
                                      -2.5520E-01,  2.0573E-04, -4.9850E+00,  1.1272E-01,  1.3267E-03 /)   !U  21
                     coeffs(:,2) = (/  2.1691E+01, -3.1787E-01,  3.2327E+01, -4.5546E-03,  1.1194E+01, &
                                      -8.0823E-01,  1.4306E-04, -4.3418E+00,  7.3163E-02,  6.3637E-03 /)   !V  21
                     coeffs(:,3) = (/  1.4634E+01, -3.9394E-01,  1.1617E+01,  5.6387E-03,  5.4799E+00, &
                                      -3.9011E-01, -1.0420E-05, -2.4279E+00,  6.6452E-02,  4.9504E-03 /)   !W  21
                  ELSEIF (WS <= 24.0) THEN  ! 22 < WS <= 24
                     coeffs(:,1) = (/  7.3816E+00,  1.0538E+00,  2.1578E+01, -3.3487E-02, -6.4986E+00, &
                                      -8.6782E-01,  3.2397E-04,  1.1412E+00,  2.2982E-01,  1.4660E-02 /)   !U  23
                     coeffs(:,2) = (/  6.5302E+00,  1.0524E+00,  2.4596E+01, -4.1648E-02,  4.0584E+00, &
                                      -6.1130E-01,  4.5468E-04, -3.6547E+00,  2.3176E-01,  8.4385E-03 /)   !V  23
                     coeffs(:,3) = (/  1.3424E+01,  2.6104E-02,  7.6014E+00, -1.2744E-02,  1.0735E+01, &
                                       2.2086E-01,  1.9309E-04, -5.9548E+00,  8.6483E-02, -3.9550E-03 /)   !W  23
                  ELSE                      ! 24 < WS
                     coeffs(:,1) = (/ -1.6629E+01,  1.3094E+00, -4.4183E+00, -8.4860E-03, -1.3800E+01, &
                                      -5.5221E-01, -5.6659E-05,  8.1834E+00, -8.2497E-03,  1.8383E-02 /)   !U  25
                     coeffs(:,2) = (/  3.4796E+00,  7.1144E-01,  1.2153E+01, -2.7309E-02,  1.0003E+00, &
                                      -6.3570E-01,  3.4424E-04, -8.5038E-01,  1.2822E-01,  1.3181E-02 /)   !V  25
                     coeffs(:,3) = (/  2.7014E+00,  1.1794E-01,  2.1378E+00,  4.5539E-03,  1.6899E+00, &
                                       1.2254E-01, -9.6940E-05, -2.3430E-01, -2.3826E-02,  5.5964E-05 /)   !W  25
                  ENDIF
               ENDIF
            ENDIF

            HT2 = HT1*HT1
            HT3 = HT1*HT2
            RI2 = RI1*RI1
            RI3 = RI1*RI2

            DO I = 1,3
               InCDec(I) = coeffs( 1,I) + coeffs(2,I)*Ht1 + coeffs(3,I)*RI1 &
                                        + coeffs(4,I)*Ht2 + coeffs(5,I)*RI2 + coeffs( 6,I)*Ht1*RI1 &
                                        + coeffs(7,I)*Ht3 + coeffs(8,I)*RI3 + coeffs( 9,I)*Ht1*RI2 &
                                                                            + coeffs(10,I)*Ht2*RI1
            ENDDO

            WS1 = MAX(  2.0_ReKi, WS )
            SELECT CASE ( Ri_Cat )
               CASE ( 1 )
!                 InCDec   = (/            1.623224368*WS1**1.015099356, &
!                                          0.884720872*WS1**1.192553093, &
!                                          1.338245093*WS1**0.841757461 /)
                  InCohB   = (/ -2.524e-05+0.002122544/WS1, &
                                           0.004367773*WS1**(-1.14945936), &
                                           0.031284497*WS1**(-0.72509517) /)

               CASE ( 2 )
!                 InCDec   = (/            1.478475074*WS1**0.752442176, &
!                                          1.310684825*WS1**0.624122449, &
!                                          0.849106068*WS1**0.627688235 /)
                  InCohB   = (/            0.003320615*WS1**(-1.18592214), &
                                           0.005402681*WS1**(-0.98637053), &
                                           0.091649927*WS1**(-1.48835650) /)

               CASE ( 3 )
!                 InCDec   = (/            1.596175944*WS1**0.674743966, &
!                                          1.114069218*WS1**0.638049141, &
!                                          0.473225245*WS1**0.784331891 /)
                  InCohB   = (/            0.002387997*WS1**(-0.85956868), &
                                           0.009481901*WS1**(-1.02518835), &
                                           0.052147706*WS1**(-0.88949864) /)

               CASE ( 4 )
!                 InCDec   = (/            1.293345620*WS1**0.955639280, &
!                                          1.296399839*WS1**0.838281755, &
!                                          0.333750239*WS1**1.103784094 /)
                  InCohB   = (/            0.002870978*WS1**(-1.07398490), &
                                           0.002435238*WS1**(-0.68685045), &
                                           0.125356016*WS1**(-1.34791890) /)

               CASE ( 5 )
!                 InCDec   = (/            1.325256941*WS1**1.039629269, &
!                                          1.014004299*WS1**1.082810576, &
!                                          0.206383058*WS1**1.435200799 /)
                  InCohB   = (/            0.003545043*WS1**(-1.03669585), &
                                           0.003996215*WS1**(-0.95313438), &
                                           0.125103070*WS1**(-1.02886635) /)
               END SELECT

         CASE ( SpecModel_WF_UPW )
            HT1 = MAX( 5.0_ReKi, MIN( Ht, 35.0_ReKi ) )
            IF ( WS <= 14.0 ) THEN
               IF ( WS <= 10 ) THEN
                  RI1 = MAX( -0.5_ReKi, MIN(  RICH_NO, 0.15_ReKi ) )
                  IF  ( WS <=  8.0 ) THEN   !      WS <= 8
                     coeffs(:,1) = (/  1.6715E+01, -3.8639E-01,  7.1817E+00,  1.5550E-03, -1.4293E+00, &
                                      -2.0350E-01,  8.5532E-06, -3.4710E+00, -1.9743E-02, -3.9949E-04 /) !Upw_U 7
                     coeffs(:,2) = (/  8.4145E+00, -4.7610E-02,  3.9097E+00, -7.1412E-04,  1.8295E+01, &
                                       2.2583E-01, -1.6965E-05,  2.0769E+01, -9.1670E-02, -8.0300E-03 /) !Upw_V 7
                  ELSE                      !  8 < WS <= 10
                     coeffs(:,1) = (/  1.5432E+01, -2.1254E-01,  5.3075E+00, -2.9928E-03,  2.1647E+00, &
                                       1.1787E-02,  6.7458E-05, -9.0445E-01, -7.5941E-02, -4.7053E-03 /) !Upw_U 9
                     coeffs(:,2) = (/  7.5921E+00,  3.3520E-02,  1.2231E+01, -7.0018E-03,  6.0889E+01, &
                                       2.1810E-01,  1.1718E-04,  7.7287E+01, -1.3828E-01, -9.6568E-03 /) !Upw_V 9
                  ENDIF
               ELSE
                  RI1 = MAX( -0.5_ReKi, MIN(  RICH_NO, 0.05_ReKi ) )
                  IF  ( WS <= 12.0 ) THEN   ! 10 < WS <= 12
                     coeffs(:,1) = (/  1.3539E+01, -8.4892E-02, -1.9237E+00, -1.1485E-03, -4.0840E-01, &
                                       3.0956E-01,  2.4048E-05, -1.1523E+00,  9.6877E-03, -4.0606E-03 /) !Upw_U 11
                     coeffs(:,2) = (/  7.7451E+00, -1.3818E-01, -9.5197E-01,  3.9610E-03,  8.3255E-01, &
                                       7.2166E-02, -4.5012E-05, -2.0948E-01, -2.1400E-02, -2.9788E-04 /) !Upw_V 11
                  ELSE                      ! 12 < WS <= 14
                     coeffs(:,1) = (/  1.2857E+01, -7.9408E-03, -1.5310E+00, -4.1077E-03,  1.0496E+00, &
                                       1.9473E-01,  7.2808E-05,  1.8380E-01, -1.6559E-02, -2.0872E-03 /) !Upw_U 13
                     coeffs(:,2) = (/  7.2452E+00, -6.2662E-02, -2.4865E+00,  3.2123E-03, -1.0281E-01, &
                                       1.9698E-01, -7.5745E-05, -1.1637E+00, -4.6458E-02, -2.7037E-03 /) !Upw_V 13
                  ENDIF
               ENDIF
            ELSE
               RI1 = MAX( -0.5_ReKi, MIN(  RICH_NO, 0.05_ReKi ) )
               IF  ( WS  <= 18.0 ) THEN
                  IF ( WS <= 16.0 ) THEN   ! 14 < WS <= 16
                     coeffs(:,1) = (/  1.4646E+01, -1.5023E-01, -9.7543E-01, -3.5607E-03,  4.8663E+00, &
                                      -9.4360E-03,  1.4932E-04,  5.9503E+00,  7.4028E-02,  5.2698E-03 /) !Upw_U 15
                     coeffs(:,2) = (/  1.0133E+01, -3.1417E-01,  2.5400E+00,  6.6777E-03,  3.0790E+00, &
                                      -2.5801E-01, -4.9501E-05,  2.8879E+00, -1.6722E-02,  4.8297E-03 /) !Upw_V 15
                  ELSE                     ! 16 < WS <= 18
                     coeffs(:,1) = (/  1.5282E+01, -2.7642E-01,  2.5903E+00,  9.8716E-03,  5.9314E-01, &
                                      -4.2790E-01, -1.6474E-04, -7.0065E-01, -3.2694E-02,  2.4583E-03 /) !Upw_U 17
                     coeffs(:,2) = (/  1.2464E+01, -3.4306E-01,  3.6261E+00,  5.8254E-03,  2.2592E+00, &
                                      -1.1498E-01, -6.6196E-05,  1.3610E+00, -1.3345E-02,  1.0932E-03 /) !Upw_V 17
                  ENDIF
               ELSE
                  IF ( WS <= 20.0 ) THEN   ! 18 < WS <= 20
                     coeffs(:,1) = (/  1.5059E+01, -8.0478E-02,  8.7088E+00, -1.7854E-03,  3.9922E+00, &
                                      -6.0268E-01,  4.3906E-05,  3.3463E+00, -6.6490E-02,  1.2290E-02 /) !Upw_U 19
                     coeffs(:,2) = (/  1.0672E+01, -2.8104E-01,  7.8021E+00,  6.6360E-03,  2.4345E+00, &
                                      -4.9103E-01, -8.3745E-05,  4.4084E-01, -9.2432E-02,  8.3096E-03 /) !Upw_V 19
                  ELSE                     ! 20 < WS
                     coeffs(:,1) = (/  1.8592E+01,  1.3888E-01,  1.6732E+01, -1.1880E-02,  2.3622E+01, &
                                       6.8199E-01,  7.3664E-05,  4.1289E+00, -3.8604E-01, -3.0381E-02 /) !Upw_U 21
                     coeffs(:,2) = (/  7.7137E+00,  1.2732E-01,  1.3477E+01,  1.9164E-03,  3.7133E+01, &
                                       3.8975E-01, -2.2818E-04,  1.8816E+01, -7.5304E-01, -2.1856E-02 /) !Upw_V 21
                  ENDIF
               ENDIF
            ENDIF

            HT2 = HT1*HT1
            HT3 = HT1*HT2
            RI2 = RI1*RI1
            RI3 = RI1*RI2

            DO I = 1,2
               InCDec(I) = coeffs( 1,I) + coeffs(2,I)*Ht1 + coeffs(3,I)*RI1 &
                                        + coeffs(4,I)*Ht2 + coeffs(5,I)*RI2 + coeffs( 6,I)*Ht1*RI1 &
                                        + coeffs(7,I)*Ht3 + coeffs(8,I)*RI3 + coeffs( 9,I)*Ht1*RI2 &
                                                                            + coeffs(10,I)*Ht2*RI1
            ENDDO

            WS1 = MAX(  3.0_ReKi, WS )
!           InCDec(1:2)   = (/             5.640176786*WS1**0.269850341, &
!                              6.059554513+18.44124731/WS1**1.5 /)
            InCohB(1:2)   = (/ 0.000448295+0.002502915/WS1, &
                               0.001539069+0.005954785/WS1 /)


            InCDec(3)     =  0.4*InCDec(1)  !cohA(w) = cohA(u)/2.5, number derived from histograms of u/w for NWTC and LLLJP data
            InCohB(3)     = 10.0*InCohB(1)  !cohB(w) = cohB(u)*10, number derived from histograms of w/u for NWTC and LLLJP data

         CASE ( SpecModel_WF_07D, SpecModel_WF_14D )
            HT1 = MAX( 5.0_ReKi, MIN( Ht, 35.0_ReKi ) )
            IF ( WS <= 12.0 ) THEN
               IF     ( WS <=  8.0 ) THEN  !      WS <= 8
                  RI1 = MAX( -0.5_ReKi, MIN(  RICH_NO, 0.15_ReKi ) )
                  coeffs(:,1) = (/  1.0310E+01, -6.4824E-03, -1.3258E+00, -2.7238E-03, -6.8515E+00, &
                                    3.1602E-02,  5.5982E-05, -8.4777E+00,  2.1506E-02,  4.9745E-04 /) !Dwn_U 7
                  coeffs(:,2) = (/  6.9491E+00, -1.3378E-01,  1.7961E-01, -4.9439E-04, -1.8140E+00, &
                                   -4.2321E-02,  4.4962E-05, -3.6939E+00, -8.9465E-03,  4.7867E-04 /) !Dwn_V 7
               ELSEIF ( WS <= 10.0 ) THEN  !  8 < WS <= 10
                  RI1 = MAX( -0.5_ReKi, MIN(  RICH_NO, 0.05_ReKi ) )
                  coeffs(:,1) = (/  9.7420E+00,  6.1610E-02,  5.6636E-02, -5.5949E-03, -1.3014E+00, &
                                    2.0655E-01,  8.9989E-05, -1.9837E+00,  5.4957E-03, -3.5496E-03 /) !Dwn_U 9
                  coeffs(:,2) = (/  7.1063E+00, -1.7021E-01,  1.2560E+00, -4.2616E-04,  9.0937E-01, &
                                   -1.3022E-01,  4.7976E-05,  2.1302E-01, -4.3159E-04,  1.5443E-03 /) !Dwn_V 9
               ELSE                        ! 10 < WS <= 12
                  RI1 = MAX( -0.5_ReKi, MIN(  RICH_NO, 0.05_ReKi ) )
                  coeffs(:,1) = (/  1.0869E+01, -9.1393E-03, -1.1695E+00, -3.3725E-03,  3.2199E-01, &
                                    7.2692E-02,  7.0565E-05,  6.9573E-01,  2.5360E-02,  1.0187E-03 /) !Dwn_U 11
                  coeffs(:,2) = (/  6.9882E+00, -1.3517E-01, -3.0492E-01, -4.6775E-04,  4.6897E-01, &
                                   -2.0102E-03,  3.3908E-05,  1.4604E-02,  1.1729E-02, -6.2775E-05 /) !Dwn_V 11
               ENDIF
            ELSE
               RI1 = MAX( -0.5_ReKi, MIN(  RICH_NO, 0.05_ReKi ) )
               IF     ( WS <= 14.0 ) THEN  ! 12 < WS <= 14
                  coeffs(:,1) = (/  1.1105E+01,  5.3789E-02, -9.4253E-02, -5.4203E-03, -1.0114E+00, &
                                    1.1421E-01,  7.6110E-05, -1.2654E+00,  1.5121E-02, -2.9055E-03 /) !Dwn_U 13
                  coeffs(:,2) = (/  7.5741E+00, -8.3945E-02,  3.7020E+00, -6.0317E-03,  3.1339E-01, &
                                   -2.1921E-01,  1.5598E-04,  6.2478E-01,  5.9490E-02,  3.4785E-03 /) !Dwn_V 13
               ELSE                        ! 14 < WS
                  coeffs(:,1) = (/  1.2256E+01,  2.0131E-02,  1.9465E+00, -7.6608E-03,  1.5031E+00, &
                                   -1.0916E-01,  1.3634E-04,  1.3451E+00, -1.6458E-02,  3.8312E-03 /) !Dwn_U 15
                  coeffs(:,2) = (/  7.7749E+00, -2.2712E-01,  1.3675E+00,  6.7944E-03,  4.2033E-02, &
                                   -6.8887E-02, -9.6117E-05, -1.5526E+00, -2.2357E-02, -1.5311E-03 /) !Dwn_V 15
               ENDIF
            ENDIF

            HT2 = HT1*HT1
            HT3 = HT1*HT2
            RI2 = RI1*RI1
            RI3 = RI1*RI2

            DO I = 1,2
               InCDec(I) = coeffs( 1,I) + coeffs(2,I)*Ht1 + coeffs(3,I)*RI1 &
                                        + coeffs(4,I)*Ht2 + coeffs(5,I)*RI2 + coeffs( 6,I)*Ht1*RI1 &
                                        + coeffs(7,I)*Ht3 + coeffs(8,I)*RI3 + coeffs( 9,I)*Ht1*RI2 &
                                                                            + coeffs(10,I)*Ht2*RI1
            ENDDO

            WS1 = MAX(  3.0_ReKi, WS )
!           WS2 = WS1*WS1
!           WS3 = WS2*WS1
!           InCDec(1:2)   = (/ (a+c*WS1+e*WS2+g*WS3)/(1+b*WS1+d*WS2+f*WS3), &
!                                               3.357892649*WS1**0.1198781 /)
            InCohB(1:2)   = (/ 4.49289e-05+0.004933460/WS1, &
                                0.00158053+0.014268899/WS1 /)
            InCDec(3)     =  0.4_ReKi*InCDec(1)  !cohA(w) = cohA(u)/2.5, number derived from histograms of u/w for NWTC and LLLJP data
            InCohB(3)     = 10.0_ReKi*InCohB(1)  !cohB(w) = cohB(u)*10, number derived from histograms of w/u for NWTC and LLLJP data

         CASE ( SpecModel_USER )
            InCDec = (/   WS, HUGE(InCohB(1)), HUGE(InCohB(1)) /)
            InCohB = 0.0_ReKi ! entire array is zero

         CASE DEFAULT   ! includes CASE ( 'SMOOTH' )

            InCDec = (/1.0_ReKi,   0.75_ReKi,  0.75_ReKi /)*WS  ! The davenport exponential parameter indicates that coh(v) ~ coh(w) in NWTC and LLLJP data
            InCohB = 0.0_ReKi ! entire array is zero

      END SELECT

         !note that the IEC models specify their coherence parameters elsewhere... in CalcIECScalingParams()

!      IF ( p%met%IsIECModel ) THEN we'll get the defaults from CalcIECScalingParams
         
         
END SUBROUTINE GetDefaultCoh
!=======================================================================
!> This subroutine is used to get the default values of the Reynolds stresses.
!! sets p%met%PC_UW,  p%met%PC_UV,  p%met%PC_VW and
!!      p%met%UWskip, p%met%UVskip, p%met%VWskip
SUBROUTINE GetDefaultRS( p, OtherSt_RandNum, TmpUstarHub, ErrStat, ErrMsg )

   ! Needs p%grid information; calls getVelocityProfile(), also 
   ! depends on uHub, ZL, Rich_No, UStar
   
   TYPE(TurbSim_ParameterType),     INTENT(INOUT)  :: p                   ! TurbSim parameters
   TYPE(RandNum_OtherStateType),    INTENT(INOUT)  :: OtherSt_RandNum     ! other states for random numbers (next seed, etc)
                                    
   REAL(ReKi),                      INTENT(IN)    :: TmpUstarHub 
   INTEGER(IntKi),                  intent(  out) :: ErrStat              ! Error level
   CHARACTER(*),                    intent(  out) :: ErrMsg               ! Message describing error
   
   
   REAL(ReKi)                                     :: rndSgn
   REAL(ReKi)                                     :: SignProb
   REAL(ReKi)                                     :: Shr
   REAL(ReKi)                                     :: Ustar2
   REAL(ReKi)                                     :: V(2)
   REAL(ReKi)                                     :: Z(2)
   REAL(ReKi)                                     :: ZLtmp
            
   INTEGER(IntKi)                                 :: ErrStat2
   CHARACTER(MaxMsgLen)                           :: ErrMsg2

   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   
   Z(2) = p%grid%HubHt + 0.5*p%grid%RotorDiameter    ! top of the grid
   Z(1) = Z(2) - p%grid%GridHeight                   ! bottom of the grid
   CALL getVelocityProfile(p, p%UHub, p%grid%HubHt, Z, V, ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'GetDefaultRS')

   Shr = ( V(2)-V(1) ) / p%grid%GridHeight    ! dv/dz

!BJJ: check the ranges of our best-fit parameters, using domains of measured values

   SELECT CASE ( p%met%TurbModel_ID )
      CASE ( SpecModel_GP_LLJ )
         ZLtmp  = MIN( MAX( p%met%ZL,    -1.00_ReKi ), 1.0_ReKi )  !Limit the observed values of z/L
         UStar2 = MIN( MAX( p%met%Ustar,  0.15_ReKi ), 1.0_ReKi )  !Limit the observed values of u*
         Ustar2 = Ustar2*Ustar2
      CASE ( SpecModel_NWTCUP )
         ZLtmp  = MIN( MAX( p%met%ZL,    -0.5_ReKi ), 3.5_ReKi )  !Limit the observed values of z/L
         UStar2 = MIN( MAX( p%met%Ustar,  0.2_ReKi ), 1.4_ReKi )  !Limit the observed values of u*
         Ustar2 = Ustar2*Ustar2
!      CASE ( 'WF_UPW' )
!      CASE ( 'WF_07D' )
!      CASE ( 'WF_14D' )

      CASE DEFAULT
         ZLtmp  = p%met%ZL
         Ustar2 = p%met%Ustar*p%met%Ustar
   END SELECT

   !-------------------------------------------------------------------------------------------------
   ! default UW Reynolds stress
   !-------------------------------------------------------------------------------------------------
   p%met%UWskip     = .FALSE.
   
   CALL  RndUnif( p%RNG, OtherSt_RandNum, rndSgn )
   SELECT CASE ( p%met%TurbModel_ID )

      CASE ( SpecModel_GP_LLJ )

         p%met%PC_UW = TmpUstarHub**2        
      
         IF (p%met%PC_UW <= 0) THEN  !We don't have a local u* value to tie it to; otherwise, assume p%met%PC_UW contains magnitude of value we want
            IF ( p%grid%HubHt >= 100.5 ) THEN     ! 116m
               p%met%PC_UW =  0.0399 - 0.00371*p%UHub - 0.00182*p%met%RICH_NO + 0.00251*ZLtmp - 0.402*Shr + 1.033*Ustar2
            ELSEIF ( p%grid%HubHt >= 76.0 ) THEN  ! 85 m
               p%met%PC_UW = 0.00668 - 0.00184*p%UHub + 0.000709*p%met%RICH_NO  + 0.264*Shr + 1.065*Ustar2  !magnitude
            ELSEIF ( p%grid%HubHt >= 60.5 ) THEN  ! 67 m
               p%met%PC_UW = -0.0216 + 0.00319*p%UHub  - 0.00205*ZLtmp + 0.206*Shr + 0.963*Ustar2    !magnitude
            ELSE                           ! 54 m
               p%met%PC_UW = -0.0373 + 0.00675*p%UHub  - 0.00277*ZLtmp + 0.851*Ustar2                !magnitude
            ENDIF
            p%met%PC_UW = MAX(p%met%PC_UW,0.0_ReKi)

         ENDIF

         IF (p%met%PC_UW > 0) THEN
            SignProb = 0.765 + 0.57/PI * ATAN( 0.78511*LOG(p%met%PC_UW)+3.42584)
            IF (rndSgn <= SignProb) p%met%PC_UW = -p%met%PC_UW
         ENDIF

      CASE ( SpecModel_NWTCUP )

         IF ( p%grid%HubHt > 47.0 ) THEN      ! 58m data
            p%met%PC_UW = 0.165 - 0.0232*p%UHub - 0.0129*p%met%RICH_NO + 1.337*Ustar2 - 0.758*SHR
         ELSEIF ( p%grid%HubHt >= 26.0 ) THEN ! 37m data
            p%met%PC_UW = 0.00279 - 0.00139*p%UHub + 1.074*Ustar2 + 0.179*SHR
         ELSE                          ! 15m data
            p%met%PC_UW = -0.1310 + 0.0239*p%UHub + 0.556*Ustar2
         ENDIF
         p%met%PC_UW = MAX(p%met%PC_UW,0.0_ReKi)

         IF (p%met%PC_UW > 0) THEN !i.e. not equal to zero
            SignProb = 0.765 + 0.57/PI * ATAN( 0.88356*LOG(p%met%PC_UW)+2.47668)
            IF (rndSgn <= SignProb) p%met%PC_UW = -p%met%PC_UW
         ENDIF

      CASE ( SpecModel_WF_14D )

         p%met%PC_UW = -Ustar2
         IF ( rndSgn > 0.9937 )  p%met%PC_UW = -p%met%PC_UW

      CASE ( SpecModel_USER, SpecModel_TimeSer )
         p%met%PC_UW = 0.0
         p%met%UWskip = .TRUE.

      CASE ( SpecModel_TIDAL, SpecModel_RIVER ) ! HYDROTURBSIM specific.
         p%met%PC_UW = -Ustar2*(1-p%grid%HubHt/p%met%RefHt) 
         
      CASE DEFAULT

         p%met%PC_UW = -Ustar2

   END SELECT

   IF ( p%met%IsIECModel ) THEN
      p%met%PC_UW = 0.0
      p%met%UWskip = .TRUE.
   END IF
      
   !-------------------------------------------------------------------------------------------------
   ! default UV Reynolds stress
   !-------------------------------------------------------------------------------------------------
   p%met%UVskip     = .FALSE.

   CALL  RndUnif( p%RNG, OtherSt_RandNum, rndSgn )
   SELECT CASE ( p%met%TurbModel_ID )

      CASE ( SpecModel_GP_LLJ )

         IF ( p%grid%HubHt >= 100.5 ) THEN     ! 116m
            p%met%PC_UV = 0.199 - 0.0167*p%UHub + 0.0115*ZLtmp + 1.143*Ustar2
            p%met%PC_UV = MAX(p%met%PC_UV,0.0_ReKi)
            IF ( rndSgn < 0.6527 ) p%met%PC_UV = -p%met%PC_UV
         ELSEIF ( p%grid%HubHt >= 76.0 ) THEN  ! 85 m
            p%met%PC_UV = 0.190 - 0.0156*p%UHub + 0.00931*ZLtmp + 1.101*Ustar2
            p%met%PC_UV = MAX(p%met%PC_UV,0.0_ReKi)
            IF ( rndSgn < 0.6394 ) p%met%PC_UV = -p%met%PC_UV
         ELSEIF ( p%grid%HubHt >= 60.5 ) THEN  ! 67 m
            p%met%PC_UV = 0.178 - 0.0141*p%UHub + 0.00709*ZLtmp + 1.072*Ustar2
            p%met%PC_UV = MAX(p%met%PC_UV,0.0_ReKi)
            IF ( rndSgn < 0.6326 ) p%met%PC_UV = -p%met%PC_UV
         ELSE                           ! 54 m
            p%met%PC_UV = 0.162 - 0.0123*p%UHub + 0.00784*p%met%RICH_NO + 1.024*Ustar2
            p%met%PC_UV = MAX(p%met%PC_UV,0.0_ReKi)
            IF ( rndSgn < 0.6191 ) p%met%PC_UV = -p%met%PC_UV
         ENDIF

      CASE ( SpecModel_NWTCUP )

            ! Get the magnitude and add the sign
         IF ( p%grid%HubHt > 47.0 ) THEN      ! 58m data
            p%met%PC_UV = 0.669 - 0.0300*p%UHub - 0.0911*p%met%RICH_NO + 1.421*Ustar2 - 1.393*SHR
         ELSEIF ( p%grid%HubHt >= 26.0 ) THEN ! 37m data
            p%met%PC_UV = 1.521 - 0.00635*p%UHub - 0.2200*p%met%RICH_NO + 3.214*Ustar2 - 3.858*SHR
         ELSE                          ! 15m data
            p%met%PC_UV = 0.462 - 0.01400*p%UHub + 1.277*Ustar2
         ENDIF
         p%met%PC_UV = MAX(p%met%PC_UV,0.0_ReKi)
         IF (p%met%PC_UV > 0) THEN !i.e. not equal to zero
            SignProb = 0.33 + 0.64/PI * ATAN( -0.374775*LOG(p%met%PC_UV)-0.205681)
            IF (rndSgn <= SignProb) p%met%PC_UV = -p%met%PC_UV
         ENDIF

      CASE ( SpecModel_WF_UPW )

         p%met%PC_UV = 0.0202 + 0.890*Ustar2 - 2.461*Shr
         p%met%PC_UV = MAX(p%met%PC_UV,0.0_ReKi)
         IF ( rndSgn < 0.7315 ) p%met%PC_UV = -p%met%PC_UV

      CASE ( SpecModel_WF_07D )

         p%met%PC_UV = 0.5040 + 0.177*Ustar2
         p%met%PC_UV = MAX(p%met%PC_UV,0.0_ReKi)
         IF ( rndSgn < 0.7355 ) p%met%PC_UV = -p%met%PC_UV

      CASE ( SpecModel_WF_14D )

         p%met%PC_UV = 0.0430 + 0.258*Ustar2
         p%met%PC_UV = MAX(p%met%PC_UV,0.0_ReKi)
         IF ( rndSgn < 0.4423 ) p%met%PC_UV = -p%met%PC_UV

      CASE DEFAULT

         p%met%PC_UV  = 0.0
         p%met%UVskip = .TRUE.  !use whatever comes our way from the random phases

   END SELECT


   !-------------------------------------------------------------------------------------------------
   ! default VW Reynolds stress
   !-------------------------------------------------------------------------------------------------
   p%met%VWskip     = .FALSE.

   CALL  RndUnif( p%RNG, OtherSt_RandNum, rndSgn )
   SELECT CASE ( p%met%TurbModel_ID )

      CASE ( SpecModel_GP_LLJ )

         IF ( p%grid%HubHt >= 100.5 ) THEN     ! 116m
            p%met%PC_VW =  0.0528  - 0.00210*p%UHub - 0.00531*p%met%RICH_NO - 0.519*Shr + 0.283*Ustar2
            p%met%PC_VW = MAX(p%met%PC_VW,0.0_ReKi)
            IF ( rndSgn < 0.2999 ) p%met%PC_VW = -p%met%PC_VW
         ELSEIF ( p%grid%HubHt >= 76.0 ) THEN  ! 85 m
            p%met%PC_VW =  0.0482  - 0.00264*p%UHub - 0.00391*p%met%RICH_NO - 0.240*Shr + 0.265*Ustar2
            p%met%PC_VW = MAX(p%met%PC_VW,0.0_ReKi)
            IF ( rndSgn < 0.3061 ) p%met%PC_VW = -p%met%PC_VW
         ELSEIF ( p%grid%HubHt >= 60.5 ) THEN  ! 67 m
            p%met%PC_VW =  0.0444  - 0.00249*p%UHub - 0.00403*p%met%RICH_NO - 0.141*Shr + 0.250*Ustar2
            p%met%PC_VW = MAX(p%met%PC_VW,0.0_ReKi)
            IF ( rndSgn < 0.3041 ) p%met%PC_VW = -p%met%PC_VW
         ELSE                           ! 54 m
            p%met%PC_VW =  0.0443  - 0.00261*p%UHub - 0.00371*p%met%RICH_NO - 0.107*Shr + 0.226*Ustar2
            p%met%PC_VW = MAX(p%met%PC_VW,0.0_ReKi)
            IF ( rndSgn < 0.3111 ) p%met%PC_VW = -p%met%PC_VW
         ENDIF

      CASE ( SpecModel_NWTCUP )

         IF ( p%grid%HubHt > 47.0 ) THEN      ! 58m data
            p%met%PC_VW = 0.174 + 0.00154*p%UHub - 0.0270*p%met%RICH_NO + 0.380*Ustar2 - 1.131*Shr - 0.00741*ZLtmp
         ELSEIF ( p%grid%HubHt >= 26.0 ) THEN ! 37m data
            p%met%PC_VW = 0.120 + 0.00283*p%UHub - 0.0227*p%met%RICH_NO + 0.306*Ustar2 - 0.825*Shr
         ELSE                          ! 15m data
            p%met%PC_VW = 0.0165 + 0.00833*p%UHub                 + 0.224*Ustar2
         ENDIF
         p%met%PC_VW = MAX(p%met%PC_VW,0.0_ReKi)
         IF (p%met%PC_VW > 0) THEN !i.e. not equal to zero
            SignProb = 0.725 + 0.65/PI * ATAN( 0.654886_ReKi*LOG(p%met%PC_VW)+1.777198_ReKi)
            IF (rndSgn <= SignProb) p%met%PC_VW = -p%met%PC_VW
         ENDIF

      CASE ( SpecModel_WF_UPW )

         p%met%PC_VW = 0.0263 + 0.273*Ustar2 - 0.684*Shr
         p%met%PC_VW = MAX(p%met%PC_VW,0.0_ReKi)
         IF ( rndSgn < 0.3139_ReKi ) p%met%PC_VW = -p%met%PC_VW

      CASE ( SpecModel_WF_07D )

         p%met%PC_VW = 0.241 + 0.118*Ustar2
         p%met%PC_VW = MAX(p%met%PC_VW,0.0_ReKi)
         IF ( rndSgn < 0.0982_ReKi ) p%met%PC_VW = -p%met%PC_VW

      CASE ( SpecModel_WF_14D )

         p%met%PC_VW =-0.0224 + 0.159*Ustar2
         p%met%PC_VW = MAX(p%met%PC_VW,0.0_ReKi)
         IF ( rndSgn < 0.8436_ReKi ) p%met%PC_VW = -p%met%PC_VW

      CASE DEFAULT

         p%met%PC_VW  = 0.0
         p%met%VWskip = .TRUE.  !use whatever comes our way from the random phases

   END SELECT


RETURN
END SUBROUTINE GetDefaultRS
!=======================================================================
!< This function calculates the default power law exponent.
FUNCTION DefaultPowerLawExp( p )

! necessary requirements: 
! Rich_No, KHtest, TurbModel_ID, IEC_WindType, IECstandard


   IMPLICIT                 NONE

   TYPE(TurbSim_ParameterType), INTENT(IN)     ::  p                             !< parameters 
   REAL(ReKi)                                  :: DefaultPowerLawExp             !< Default Power Law exponent for particular model

   IF ( p%met%KHtest ) THEN
      DefaultPowerLawExp = 0.3
      RETURN
   ENDIF

   SELECT CASE ( p%met%TurbModel_ID )

      CASE (SpecModel_WF_UPW, SpecModel_NWTCUP)
         IF ( p%met%RICH_NO > 0.0 ) THEN
            DefaultPowerLawExp = 0.14733
         ELSE
            DefaultPowerLawExp = 0.087687698 + 0.059641545*EXP(p%met%RICH_NO/0.04717783)
         ENDIF

      CASE ( SpecModel_WF_07D, SpecModel_WF_14D )
         IF ( p%met%RICH_NO > 0.04 ) THEN
            DefaultPowerLawExp = 0.17903
         ELSE
            DefaultPowerLawExp = 0.127704032 + 0.031228952*EXP(p%met%RICH_NO/0.0805173)
         ENDIF

      CASE (SpecModel_SMOOTH, SpecModel_GP_LLJ, SpecModel_TIDAL, SpecModel_RIVER, SpecModel_TimeSer )
         ! A 1/7 power law seems to work ok for HYDRO spectral models also...
         DefaultPowerLawExp = 0.143

      CASE DEFAULT
         IF ( p%IEC%IEC_WindType == IEC_EWM1 .OR. p%IEC%IEC_WindType == IEC_EWM50 .OR. p%IEC%IEC_WindType == IEC_EWM100 ) THEN
            DefaultPowerLawExp = 0.11         ! [IEC 61400-1 6.3.2.1 (14)]
         ELSEIF ( p%IEC%IECstandard == 3 ) THEN
            DefaultPowerLawExp = 0.14         ! [IEC 61400-3 Page 22 (3)]
         ELSE
            DefaultPowerLawExp = 0.2          ! [IEC 61400-1 6.3.1.2 (10)]
         ENDIF

   END SELECT

   RETURN
END FUNCTION DefaultPowerLawExp
!=======================================================================
FUNCTION getUstarDiab(u_ref, z_ref, z0, ZL)


   IMPLICIT                              NONE

   REAL(ReKi),   INTENT(IN)           :: u_ref                       ! Wind speed at reference height
   REAL(ReKi),   INTENT(IN)           :: z_ref                       ! Reference height
   REAL(ReKi),   INTENT(IN)           :: z0                          ! Surface roughness length -- It must be > 0 (which we've already checked for)
   REAL(ReKi),   INTENT(IN)           :: ZL                          ! M-O stability parameter

   REAL(ReKi)                         :: tmp                         ! a temporary value
   REAL(ReKi)                         :: psiM
   REAL(ReKi)                         :: getUstarDiab                ! the diabatic u* value (u*0)

   IF ( ZL >= 0 ) THEN !& ZL < 1
      psiM = -5.0*MIN(ZL, REAL(1.0,ReKi) )
   ELSE
      tmp = (1.0 - 15.0*ZL)**0.25

      !psiM = -2.0*LOG( (1.0 + tmp)/2.0 ) - LOG( (1.0 + tmp*tmp)/2.0 ) + 2.0*ATAN( tmp ) - 0.5 * PI
      psiM = -LOG( 0.125 * ( (1.0 + tmp)**2 * (1.0 + tmp*tmp) ) ) + 2.0*ATAN( tmp ) - 0.5 * PI
      
      !bjj 11-may-2016: because of the negative sign in the equation below, I believe psiM needs to switch signs.
      ! if true, this has been implemented incorrectly for at least 15 years.
      psiM = -psiM

   ENDIF

   getUstarDiab = ( 0.4 * u_ref ) / ( LOG( z_ref / z0 ) - psiM )

END FUNCTION getUstarDiab
!=======================================================================
!> this routine calculates the M-O z/L and L parameters using
!!  Rich_No, SpecModel, and HubHt
SUBROUTINE Calc_MO_zL(SpecModel, Rich_No, HubHt, ZL, L )


   IMPLICIT NONE
                
   REAL(ReKi)    , intent(in)               :: HubHt                                    ! Hub height
   REAL(ReKi)    , intent(in)               :: RICH_NO                                  ! Gradient Richardson number                
   REAL(ReKi)    , intent(  out)            :: L                                        ! M-O length
   REAL(ReKi)    , intent(  out)            :: ZL                                       ! A measure of stability
   
   INTEGER(IntKi), intent(in)               :: SpecModel                                ! Integer value of spectral model (see SpecModel enum)
   

      ! ***** Calculate M-O z/L parameter   :  z/L is a number in (-inf, 1] *****

   IF ( SpecModel == SpecModel_NWTCUP ) THEN
         ! Calculate disk averaged Z/L from turbine layer Ri for NWTC/LIST experiment

      IF ( RICH_NO <= -0.1 ) THEN
         ZL = -0.254 + 1.047*RICH_NO
      ELSEIF ( RICH_NO < 0 ) THEN
         ZL = 10.369*RICH_NO/(1.0 - 19.393*RICH_NO)
      ELSE  !( RICH_NO < 0.155 ) THEN
         ZL = 2.535*MIN( RICH_NO, 0.155_ReKi ) / (1.0 - 6.252*MIN( RICH_NO, 0.155_ReKi ))
      ENDIF


   ELSEIF (SpecModel == SpecModel_GP_LLJ) THEN

      IF ( RICH_NO <= -0.1 ) THEN
         ZL = -0.047 + 1.054*RICH_NO
      ELSEIF ( RICH_NO < 0 ) THEN
         ZL = 2.213*RICH_NO/(1.0 - 4.698*RICH_NO)
      ELSE  !( RICH_NO < 0.1367 ) THEN
         ZL = 3.132*MIN( RICH_NO, 0.1367_ReKi ) / (1.0 - 6.762*MIN( RICH_NO, 0.1367_ReKi ))
      ENDIF

   ELSE ! see Businger, J.A.; Wyngaard, J.C.; Izumi, Y.; Bradley, E.F. (1971). "Flux-Profile Relationships in the Atmospheric Surface Layer." Journal of the Atmospheric Sciences (28); pp.181-189.

      IF ( RICH_NO <= 0.0_ReKi ) THEN
         ZL = RICH_NO
         !PhiM = (1.0 - 16.0*ZL)**-0.25
      ELSEIF ( RICH_NO < 0.16667_ReKi ) THEN
         ZL = MIN(RICH_NO / ( 1.0_ReKi - 5.0_ReKi*RICH_NO ), 1.0_ReKi )  ! The MIN() will take care of rounding issues.
         !PhiM = (1.0 + 5.0*ZL)
      ELSE
         ZL = 1.0_ReKi
      ENDIF

   ENDIF !SpecModels

   ZL = MIN( ZL, 1.0_ReKi )

   
      ! ***** Calculate M-O length scale, L [meters] *****
      ! L should be constant in the surface layer

   IF ( .NOT. EqualRealNos(ZL , 0.0_ReKi) ) THEN
      L = HubHt / ZL ! Since ZL is the average ZL over the rotor disk, we should use HubHt to estimate L instead
   ELSE
      L = HUGE( L )
   ENDIF


END SUBROUTINE Calc_MO_zL
!=======================================================================
SUBROUTINE CalcIECScalingParams( p_IEC, HubHt, UHub, InCDec, InCohB, TurbModel_ID, IsIECModel, ErrStat, ErrMsg )
! REQUires these be set prior to calling:NumTurbInp, IECedition, IECTurbC, IEC_WindType, IsIECModel
! calculates SigmaIEC, Lambda, IntegralScale, Lc

   TYPE(IEC_ParameterType), INTENT(INOUT) :: p_IEC                       ! parameters for IEC models
   REAL(ReKi)             , INTENT(IN)    :: HubHt                       ! Hub-height
   REAL(ReKi)             , INTENT(IN)    :: UHub                        ! Hub-height (total) wind speed (m/s)
   
   REAL(ReKi)             , INTENT(OUT)   :: InCDec     (3)              ! Contains the coherence decrements
   REAL(ReKi)             , INTENT(OUT)   :: InCohB     (3)              ! Contains the coherence b/L (offset) parameters
   INTEGER(IntKi)         , INTENT(IN)    :: TurbModel_ID                ! Integer value of spectral model (see SpecModel enum)      
   LOGICAL                , INTENT(IN)    :: IsIECModel                  ! Determines if this is actually an IEC model, or if we just set the values to 0 and return      
   INTEGER(IntKi),          intent(  out) :: ErrStat                       !< Error level
   CHARACTER(*),            intent(  out) :: ErrMsg                        !< Message describing error
   
      
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   IF ( .NOT. IsIECModel )  THEN  
      
       p_IEC%SigmaIEC = 0
       p_IEC%Lambda = 0
       p_IEC%IntegralScale = 0
       p_IEC%LC = 0.0    ! The length scale is not defined for the non-IEC models
    
       RETURN
       
   ENDIF !  TurbModel  == 'IECKAI', 'IECVKM', 'API', or 'MODVKM'
   
   
   
      ! If IECKAI or IECVKM spectral models are specified, determine turb intensity 
      ! and slope of Sigma wrt wind speed from IEC turbulence characteristic, 
      ! IECTurbC = A, B, or C or from user specified quantity.

   
   IF ( p_IEC%NumTurbInp )  THEN
   
         ! user specified a particular percent TI:
         
      p_IEC%TurbInt     = 0.01*p_IEC%PerTurbInt
      p_IEC%SigmaIEC(1) = p_IEC%TurbInt*UHub
      
      ! bjj: note Vave isn't set in this case, but we only print it to the summary file (and use it) if .not. NumTurbInp      

   ELSE

      
      SELECT CASE (p_IEC%IECedition)

         CASE ( 2 )

            IF ( p_IEC%IECTurbC  == 'A' ) THEN
               p_IEC%TurbInt15  = 0.18
               p_IEC%SigmaSlope = 2.0
            ELSEIF ( p_IEC%IECTurbC  == 'B' ) THEN
               p_IEC%TurbInt15  = 0.16
               p_IEC%SigmaSlope = 3.0
            ELSE   ! We should never get here, but just to be complete...
               ErrStat = ErrID_Fatal
               ErrMsg =  'CalcIECScalingParams: Invalid IEC turbulence characteristic.'
               RETURN
            ENDIF

            p_IEC%SigmaIEC(1) = p_IEC%TurbInt15*( ( 15.0 + p_IEC%SigmaSlope*UHub ) / ( p_IEC%SigmaSlope + 1.0 ) )
            p_IEC%TurbInt     = p_IEC%SigmaIEC(1)/UHub
         
         CASE ( 3 )

            IF ( p_IEC%IECTurbC == 'A' ) THEN
               p_IEC%TurbInt15  = 0.16
            ELSEIF ( p_IEC%IECTurbC == 'B' ) THEN
               p_IEC%TurbInt15  = 0.14
            ELSEIF ( p_IEC%IECTurbC == 'C' ) THEN
               p_IEC%TurbInt15  = 0.12
            ELSE   ! We should never get here, but just to be complete...
               ErrStat = ErrID_Fatal
               ErrMsg =  'CalcIECScalingParams: Invalid IEC turbulence characteristic.'
               RETURN
            ENDIF  

                   
            SELECT CASE ( p_IEC%IEC_WindType )
               CASE ( IEC_NTM )
                  p_IEC%SigmaIEC(1) = p_IEC%TurbInt15*( 0.75*UHub + 5.6 )                                      ! [IEC-1 Ed3 6.3.1.3 (11)]
               CASE ( IEC_ETM )
                  p_IEC%Vave        = 0.2*p_IEC%Vref                                                           ! [IEC-1 Ed3 6.3.1.1 ( 9)]
                  p_IEC%SigmaIEC(1) = p_IEC%ETMc * p_IEC%TurbInt15 * ( 0.072 * &
                                     ( p_IEC%Vave / p_IEC%ETMc + 3.0) * (Uhub / p_IEC%ETMc - 4.0)+10.0 )       ! [IEC-1 Ed3 6.3.2.3 (19)]
               CASE ( IEC_EWM1, IEC_EWM50, IEC_EWM100 )
                  p_IEC%Vave        = 0.2*p_IEC%Vref                                                           ! [IEC-1 Ed3 6.3.1.1 ( 9)]
                  p_IEC%SigmaIEC(1) = 0.11*Uhub                                                                ! [IEC-1 Ed3 6.3.2.1 (16)]
               CASE DEFAULT 
                  ErrStat = ErrID_Fatal
                  ErrMsg =  'CalcIECScalingParams: Invalid IEC wind type.'
                  RETURN
            END SELECT           
            p_IEC%TurbInt  = p_IEC%SigmaIEC(1)/UHub     
            
         CASE DEFAULT ! Likewise, this should never happen...

            ErrStat = ErrID_Fatal
            ErrMsg =  'CalcIECScalingParams: Invalid IEC 61400-1 edition number.'
            RETURN
            
         END SELECT                            
      

   ENDIF

   ! note PLExp for IEC is set elsewhere
   
      ! IEC turbulence scale parameter, Lambda(1), and IEC coherency scale parameter, LC

   IF ( p_IEC%IECedition == 2 ) THEN  
      
         ! section 6.3.1.3 Eq. 9
      IF ( HubHt < 30.0_ReKi )  THEN
         p_IEC%Lambda(1) = 0.7*HubHt
      ELSE
         p_IEC%Lambda(1) = 21.0
      ENDIF

      p_IEC%LC = 3.5*p_IEC%Lambda(1)
      InCDec = (/  8.80_ReKi, HUGE(p_IEC%LC), HUGE(p_IEC%LC) /)   ! u-, v-, and w-component coherence decrement

   ELSE !IF (p_IEC%IECedition == 3 ) THEN
      
         ! section 6.3.1.3 Eq. 9      
         
      IF ( HubHt < 60.0_ReKi )  THEN
         p_IEC%Lambda(1) = 0.7*HubHt
      ELSE
         p_IEC%Lambda(1) = 42.0
      ENDIF

      p_IEC%LC = 8.1*p_IEC%Lambda(1)
      InCDec = (/ 12.00_ReKi, HUGE(p_IEC%LC), HUGE(p_IEC%LC) /)   ! u-, v-, and w-component coherence decrement for IEC Ed. 3

   ENDIF
   
   InCohB  =  (/ 0.12_ReKi/p_IEC%LC, 0.0_ReKi, 0.0_ReKi /)
         
   
      ! Set Lambda for Modified von Karman model:
#ifdef MVK
!bjj: this will probably need to be rethought with TurbSim v2.0
   IF ( MVK .AND. TurbModel_ID  == SpecModel_MODVKM ) THEN
      p%met%z0 = FindZ0(HubHt, p_IEC%SigmaIEC(1), UHub, p%met%Fc)
      CALL ScaleMODVKM(HubHt, UHub, p_IEC%Lambda(1), p_IEC%Lambda(2), p_IEC%Lambda(3))
   ENDIF   
#endif
   
      ! Sigma for v and w components and
      ! Integral scales (which depend on lambda)
   
   IF ( TurbModel_ID == SpecModel_IECVKM ) THEN
      
      p_IEC%SigmaIEC(2)      =  1.0*p_IEC%SigmaIEC(1)
      p_IEC%SigmaIEC(3)      =  1.0*p_IEC%SigmaIEC(1)
      
      p_IEC%IntegralScale(:) =  3.5 *p_IEC%Lambda(1)   !L_k
            
   ELSE
      
      p_IEC%SigmaIEC(2)      =  0.8*p_IEC%SigmaIEC(1)
      p_IEC%SigmaIEC(3)      =  0.5*p_IEC%SigmaIEC(1)
            
      p_IEC%IntegralScale(1) =  8.1 *p_IEC%Lambda(1)   !L_k
      p_IEC%IntegralScale(2) =  2.7 *p_IEC%Lambda(1)   !L_k
      p_IEC%IntegralScale(3) =  0.66*p_IEC%Lambda(1)   !L_k
      
   END IF
   
   

END SUBROUTINE CalcIECScalingParams
!=======================================================================
!> Routine sets the default ETMc, WindProfileType, Z0, Latitude, CohExp.
!! These depend on p%met%TurbModel_ID and p%usr%NPoints.
!!
!! URef, ZJetMax, and PLExp are initialized, but cannot calculate their 
!! default values until the richardson number or ustar are set.
SUBROUTINE DefaultMetBndryCndtns(p)


TYPE(TurbSim_ParameterType), INTENT(INOUT)  :: p                            !< TurbSim parameters

      ! default ETMc
   p%IEC%ETMc = 2.0_ReKi

   
      ! default WindProfileType      
   SELECT CASE ( p%met%TurbModel_ID )
      CASE ( SpecModel_TimeSer )
         IF ( p%usr%NPoints > 1 ) THEN
            p%met%WindProfileType = 'TS'
         ELSE
            p%met%WindProfileType = 'PL'
            call WrScr( 'Warning: WindProfileType will default to power-law profile because only one time-series point was entered.')
         END IF
      CASE ( SpecModel_GP_LLJ )
         p%met%WindProfileType = 'JET'
               
      CASE ( SpecModel_TIDAL )
         p%met%WindProfileType = 'H2L'
         
      CASE ( SpecModel_USRVKM )
         p%met%WindProfileType = 'USR'
         
      CASE ( SpecModel_API )
         p%met%WindProfileType = 'API'  ! ADDED BY YG
         
      CASE DEFAULT
!         p%met%WindProfileType = 'IEC'
         p%met%WindProfileType = 'PL'
   END SELECT


      ! Initialize ZJetMax (will need to set default later)
   p%met%ZJetMax    = 0.0_ReKi  
   
      ! Initialize URef (will need to set default later)
   p%met%URef       = 0.0_ReKi  

      ! Initialize PLExp (will need to set default later)
   p%met%PLExp   = 0.0_ReKi ! DefaultPowerLawExp( p )  For some models, this routine requires RICH_NO, which we do not know, yet. We'll call DefaultPowerLawExp later for all cases
      
      ! Default Z0
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
   
      ! Default Latitude   
   p%met%Latitude   = 45.0      
   
      ! Default CohExp
   p%met%CohExp     = 0.0    ! was 0.25
   
   
END SUBROUTINE DefaultMetBndryCndtns
!=======================================================================
!> Calculate the default mixing layer depth, ZI,
!! based on Ustar, UstarDiab, Uref, RefHt, Z0, Fc
SUBROUTINE DefaultMixingLayerDepth(p)
! 
   TYPE(TurbSim_ParameterType), INTENT(INOUT)  :: p  !< TurbSim parameters

      IF ( p%met%Ustar < p%met%UstarDiab ) THEN
         p%met%ZI = ( 0.04 * p%met%Uref ) / ( 1.0E-4 * LOG10( p%met%RefHt / p%met%Z0 ) )  !for "very" windy days
      ELSE
         !Should Wind Farm models use the other definition since that was what was used in creating those models?
         p%met%ZI = p%met%Ustar / (6.0 * p%met%Fc)
      ENDIF

END SUBROUTINE DefaultMixingLayerDepth
!=======================================================================
!> Calculate the default UStar value, based on
!! UstarDiab (URef, RefHt, Z0, ZL),  TurbModel_ID, ZL, URef
SUBROUTINE DefaultUstar(p)

   TYPE(TurbSim_ParameterType), INTENT(INOUT)  :: p  !< TurbSim parameters
     
   
   p%met%UstarDiab  = getUstarDiab(p%met%URef, p%met%RefHt, p%met%z0, p%met%ZL)
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
         p%met%Ustar = 0.17454 + 0.72045*p%met%UstarDiab**1.36242

      CASE ( SpecModel_NWTCUP  )
         p%met%Ustar = 0.2716 + 0.7573*p%met%UstarDiab**1.2599

      CASE ( SpecModel_TIDAL , SpecModel_RIVER )
         ! Use a constant drag coefficient for the HYDRO spectral models.
         p%met%Ustar = p%met%Uref*0.05 ! This corresponds to a drag coefficient of 0.0025.
         !p%met%Ustar = p%met%Uref*0.04 ! This corresponds to a drag coefficient of 0.0016.
         
      CASE DEFAULT
         p%met%Ustar = p%met%UstarDiab

   END SELECT
END SUBROUTINE DefaultUstar
!=======================================================================
!> Calculate the default ZJetMax value, based on
!! Rich_No, ZL, Ustar, plus A random amount
SUBROUTINE DefaultZJetMax(p, OtherSt_RandNum)

   TYPE(TurbSim_ParameterType), INTENT(INOUT)  :: p                  !< TurbSim parameters
   TYPE(RandNum_OtherStateType),INTENT(INOUT)  :: OtherSt_RandNum    !< other states for random numbers (next seed, etc)
   
   REAL(ReKi)                                  :: RandomValue

   
      ! values based on Neil Kelley's analysis
   p%met%ZJetMax = -14.820561*p%met%Rich_No + 56.488123*p%met%ZL + 166.499069*p%met%Ustar + 188.253377
   p%met%ZJetMax =   1.9326  *p%met%ZJetMax - 252.7267  ! Correct with the residual

   CALL RndJetHeight( p%RNG, OtherSt_RandNum, RandomValue ) ! Add a random amount

   p%met%ZJetMax = MIN( MAX(p%met%ZJetMax + RandomValue, ZJetMax_LB ), ZJetMax_UB )

END SUBROUTINE DefaultZJetMax
!=======================================================================
!!> Calculate the default UStar value, based on
!!! UstarDiab, TurbModel_ID, ZL, URef
!SUBROUTINE Default(p)
!
!   TYPE(TurbSim_ParameterType), INTENT(INOUT)  :: p  !< TurbSim parameters
!
!END SUBROUTINE Default

!=======================================================================
!< Routine to set parameters for the JET profile: UJetMax, ChebyCoef_WS, ChebyCoef_WD
!! will also calculate default URef if requested.
!! needs ZJetMax, RefHt, URef (unless asked to calculated here) set prior to calling
SUBROUTINE getJetCoeffs( p, getDefaultURef, OtherSt_RandNum, ErrStat, ErrMsg )

   TYPE(TurbSim_ParameterType),  INTENT(INOUT) :: p                                              !< parameters for TurbSim
   TYPE(RandNum_OtherStateType), INTENT(INOUT) :: OtherSt_RandNum                                !< other states for random number generation
   LOGICAL                     , INTENT(IN   ) :: getDefaultURef                                 !< determines if we also calculate a default URef
   INTEGER(IntKi)  ,             INTENT(  OUT) :: ErrStat                                        !< error level/status
   CHARACTER(*) ,                INTENT(  OUT) :: ErrMsg                                         !< error message

      ! local variables
   REAL(ReKi)                                  :: RandomValue
   REAL(ReKi)                                  :: URef
   INTEGER(IntKi)                              :: ErrStat2
   CHARACTER(MaxMsgLen)                        :: ErrMsg2
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   
   IF ( getDefaultURef ) THEN ! Calculate a default value

      p%met%UJetMax = MAX( -21.5515_ReKi + 6.6827_ReKi*LOG(p%met%ZJetMax), 5.0_ReKi ) !Jet max must be at least 5 m/s (occurs ~50 m); shouldn't happen, but just in case....

      CALL Rnd3ParmNorm( p%RNG, OtherSt_RandNum, RandomValue, 0.1076_ReKi, -0.1404_ReKi, 3.6111_ReKi,  -15.0_ReKi, 20.0_ReKi, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'getJetCoeffs')

      IF (p%met%UJetMax + RandomValue > 0 ) p%met%UJetMax = p%met%UJetMax + RandomValue

      CALL GetChebCoefs( p, .TRUE. , ErrStat2, ErrMsg2 ) ! These coefficients are a function of UJetMax, ZJetMax, RICH_NO, and p%met%Ustar
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'getJetCoeffs')

      CALL getVelocity(p, p%met%UJetMax, p%met%ZJetMax, p%met%RefHt, URef, ErrStat2, ErrMsg2)  
         p%met%URef = URef
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'getJetCoeffs')
         
         
   ELSE !IF ( trim(p%met%WindProfileType) == 'JET' ) then
         IF ( EqualRealNos( p%met%RefHt, p%met%ZJetMax ) ) THEN
            p%met%UJetMax = p%met%URef
            CALL GetChebCoefs( p, .TRUE. , ErrStat2, ErrMsg2 ) ! These coefficients are a function of UJetMax, ZJetMax, RICH_NO, and p%met%Ustar
         ELSE         
            CALL GetChebCoefs(p, .FALSE., ErrStat2, ErrMsg2) ! also calculate p%met%UJetMax
         END IF
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'getJetCoeffs')
         
   ENDIF !Jet wind profile
   
END SUBROUTINE getJetCoeffs
!=======================================================================

END MODULE TS_FileIO
