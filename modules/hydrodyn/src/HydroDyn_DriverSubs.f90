!**********************************************************************************************************************************
! HydroDyn_DriverSubs
!..................................................................................................................................
! LICENSING
! Copyright (C) 2022 Envision Energy USA LTD
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

MODULE HydroDynDriverSubs

   USE NWTC_Library
   use SeaState
   use SeaState_Types
   USE HydroDyn
   USE HydroDyn_Types
   USE HydroDyn_Output
   USE ModMesh_Types
   USE VersionInfo
   
   IMPLICIT NONE
   
   TYPE HD_Drvr_InitInput
      LOGICAL                 :: Echo
      REAL(ReKi)              :: Gravity
      REAL(ReKi)              :: WtrDens
      REAL(ReKi)              :: WtrDpth
      REAL(ReKi)              :: MSL2SWL
      CHARACTER(1024)         :: HDInputFile
      CHARACTER(1024)         :: SeaStateInputFile
      CHARACTER(1024)         :: OutRootName
      LOGICAL                 :: Linearize
      INTEGER                 :: NSteps
      REAL(DbKi)              :: TimeInterval
      INTEGER                 :: PRPInputsMod
      CHARACTER(1024)         :: PRPInputsFile
      REAL(ReKi)              :: uPRPInSteady(6)
      REAL(ReKi)              :: uDotPRPInSteady(6)
      REAL(ReKi)              :: uDotDotPRPInSteady(6)
      LOGICAL                 :: WaveElevSeriesFlag      !< Should we put together a wave elevation series and save it to file?
      REAL(ReKi)              :: WaveElevdX              !< Spacing in the X direction for wave elevation series              (m)
      REAL(ReKi)              :: WaveElevdY              !< Spacing in the Y direction for the wave elevation series          (m)
      INTEGER(IntKi)          :: WaveElevNX              !< Number of points in the X direction for the wave elevation series (-)
      INTEGER(IntKi)          :: WaveElevNY              !< Number of points in the X direction for the wave elevation series (-)
   END TYPE HD_Drvr_InitInput
   
! -----------------------------------------------------------------------------------   
! NOTE:  this module and the ModMesh.f90 modules must use the Fortran compiler flag:  
!        /fpp                  because of they both have preprocessor statements
! ----------------------------------------------------------------------------------- 

   !
   !INTEGER(IntKi), PARAMETER                           :: NumInp = 1           ! Number of inputs sent to HydroDyn_UpdateStates
   !
   !   ! Program variables
   !
   !REAL(DbKi)                                          :: Time                 ! Variable for storing time, in seconds
   !
   !REAL(DbKi)                                          :: InputTime(NumInp)    ! Variable for storing time associated with inputs, in seconds
   !REAL(DbKi)                                          :: Interval             ! HD module requested time interval
   !
   !type(SeaSt_InitInputType)                        :: InitInData_SeaSt           ! Input data for initialization
   !type(SeaSt_InitOutputType)                       :: InitOutData_SeaSt          ! Output data from initialization
   !
   !type(SeaSt_ContinuousStateType)                  :: x_SeaSt                    ! Continuous states
   !type(SeaSt_DiscreteStateType)                    :: xd_SeaSt                   ! Discrete states
   !type(SeaSt_ConstraintStateType)                  :: z_SeaSt                    ! Constraint states
   !type(SeaSt_OtherStateType)                       :: OtherState_SeaSt           ! Other states
   !type(SeaSt_MiscVarType)                          :: m_SeaSt                    ! Misc/optimization variables
   !
   !type(SeaSt_ParameterType)                        :: p_SeaSt                    ! Parameters
   !!type(SeaSt_InputType)                           :: u                    ! System inputs [OLD STYLE]
   !type(SeaSt_InputType)                            :: u_SeaSt(NumInp)            ! System inputs
   !type(SeaSt_OutputType)                           :: y_SeaSt                    ! System outputs
   !
   !
   !
   !TYPE(HydroDyn_InitInputType)                        :: InitInData           ! Input data for initialization
   !TYPE(HydroDyn_InitOutputType)                       :: InitOutData          ! Output data from initialization
   !
   !TYPE(HydroDyn_ContinuousStateType)                  :: x                    ! Continuous states
   !TYPE(HydroDyn_ContinuousStateType)                  :: x_new                ! Continuous states at updated time
   !TYPE(HydroDyn_DiscreteStateType)                    :: xd                   ! Discrete states
   !TYPE(HydroDyn_DiscreteStateType)                    :: xd_new               ! Discrete states at updated time
   !TYPE(HydroDyn_ConstraintStateType)                  :: z                    ! Constraint states
   !TYPE(HydroDyn_OtherStateType)                       :: OtherState           ! Other states
   !TYPE(HydroDyn_MiscVarType)                          :: m                    ! Misc/optimization variables
   !
   !TYPE(HydroDyn_ParameterType)                        :: p                    ! Parameters
   !!TYPE(HydroDyn_InputType)                           :: u                    ! System inputs [OLD STYLE]
   !TYPE(HydroDyn_InputType)                            :: u(NumInp)            ! System inputs
   !TYPE(HydroDyn_OutputType)                           :: y                    ! System outputs
   !
   !
   !
   !INTEGER(IntKi)                                     :: UnPRPInp            ! PRP Inputs file identifier
   !REAL(ReKi), ALLOCATABLE                            :: PRPin(:,:)          ! Variable for storing time, forces, and body velocities, in m/s or rad/s for PRP
   !
   !INTEGER(IntKi)                                     :: NBody                 ! Number of WAMIT bodies to work with if prescribing kinematics on each body (PRPInputsMod<0)
   !
   !INTEGER(IntKi)                                     :: I                    ! Generic loop counter
   !INTEGER(IntKi)                                     :: J                    ! Generic loop counter
   !INTEGER(IntKi)                                     :: n                    ! Loop counter (for time step)
   !INTEGER(IntKi)                                     :: ErrStat,ErrStat2     ! Status of error message
   !CHARACTER(1024)                                    :: ErrMsg,ErrMsg2       ! Error message if ErrStat /= ErrID_None
   !REAL(R8Ki)                                         :: dcm (3,3)            ! The resulting transformation matrix from X to x, (-).
   !CHARACTER(1024)                                    :: drvrFilename         ! Filename and path for the driver input file.  This is passed in as a command line argument when running the Driver exe.
   !TYPE(HD_Drvr_InitInput)                            :: drvrInitInp          ! Initialization data for the driver program
   !
   !integer                                        :: StrtTime (8)                            ! Start time of simulation (including intialization)
   !integer                                        :: SimStrtTime (8)                         ! Start time of simulation (after initialization)
   !real(ReKi)                                     :: PrevClockTime                           ! Clock time at start of simulation in seconds
   !real(ReKi)                                     :: UsrTime1                                ! User CPU time for simulation initialization
   !real(ReKi)                                     :: UsrTime2                                ! User CPU time for simulation (without intialization)
   !real(DbKi)                                     :: TiLstPrn                                ! The simulation time of the last print
   !real(DbKi)                                     :: SttsTime                                ! Amount of time between screen status messages (sec)
   !integer                                        :: n_SttsTime                              ! Number of time steps between screen status messages (-)
   !
   !type(MeshMapType)                              :: HD_Ref_2_WB_P                           ! Mesh mapping between Reference pt mesh and WAMIT body(ies) mesh
   !type(MeshMapType)                              :: HD_Ref_2_M_P                            ! Mesh mapping between Reference pt mesh and Morison mesh
   !
   !! For testing
   !REAL(DbKi)                                         :: maxAngle             ! For debugging, see what the largest rotational angle input is for the simulation
   !
   !CHARACTER(20)                    :: FlagArg       ! Flag argument from command line
   !
   TYPE(ProgDesc), PARAMETER        :: version   = ProgDesc( 'HydroDyn Driver', '', '' )  ! The version number of this program.


CONTAINS

SUBROUTINE ReadDriverInputFile( inputFile, InitInp, ErrStat, ErrMsg )

   CHARACTER(*),                  INTENT( IN    )   :: inputFile
   TYPE(HD_Drvr_InitInput),       INTENT(   OUT )   :: InitInp
   INTEGER,                       INTENT(   OUT )   :: ErrStat              ! returns a non-zero value when an error occurs  
   CHARACTER(*),                  INTENT(   OUT )   :: ErrMsg               ! Error message if ErrStat /= ErrID_None
   
      ! Local variables  

   INTEGER                                          :: UnIn                 ! Unit number for the input file
   INTEGER                                          :: UnEchoLocal          ! The local unit number for this module's echo file
   CHARACTER(1024)                                  :: EchoFile             ! Name of HydroDyn echo file  
   CHARACTER(1024)                                  :: PriPath              ! Temporary storage for relative path name
   CHARACTER(1024)                                  :: FileName             ! Name of HydroDyn input file  

   integer(IntKi)                                   :: errStat2      ! temporary error status of the operation
   character(ErrMsgLen)                             :: errMsg2       ! temporary error message 
   character(*), parameter                          :: RoutineName = 'ReadDriverInputFile'
   
   
      ! Initialize the echo file unit to -1 which is the default to prevent echoing, we will alter this based on user input
   UnEchoLocal = -1
   ErrStat = ErrID_None
   ErrMsg = ""
   
   FileName = TRIM(inputFile)
   
   CALL GetNewUnit( UnIn ) 
   CALL OpenFInpFile ( UnIn, FileName, ErrStat2, ErrMsg2 ) 
   if (Failed()) return


   CALL WrScr( 'Opening HydroDyn Driver input file:  '//FileName )
   call GetPath( TRIM(inputFile), PriPath ) ! store path name in case any of the file names are relative to the primary input file

   
   !-------------------------------------------------------------------------------------------------
   ! File header
   !-------------------------------------------------------------------------------------------------
   
   CALL ReadCom( UnIn, FileName, 'HydroDyn Driver input file header line 1', ErrStat2, ErrMsg2 )
   if (Failed()) return


   CALL ReadCom( UnIn, FileName, 'HydroDyn Driver input file header line 2', ErrStat2, ErrMsg2 )
   if (Failed()) return


     ! Echo Input Files.
   CALL ReadVar ( UnIn, FileName, InitInp%Echo, 'Echo', 'Echo Input', ErrStat2, ErrMsg2 )
   if (Failed()) return
   
   
      ! If we are Echoing the input then we should re-read the first three lines so that we can echo them
      ! using the NWTC_Library routines.  The echoing is done inside those routines via a global variable
      ! which we must store, set, and then replace on error or completion.
      
   IF ( InitInp%Echo ) THEN
      
      EchoFile = TRIM(FileName)//'.ech'
      CALL GetNewUnit( UnEchoLocal )   
      CALL OpenEcho ( UnEchoLocal, EchoFile, ErrStat2, ErrMsg2 )
      if (Failed()) return

      
      REWIND(UnIn)
      
      CALL ReadCom( UnIn, FileName, 'HydroDyn Driver input file header line 1', ErrStat2, ErrMsg2, UnEchoLocal )
      if (Failed()) return

      CALL ReadCom( UnIn, FileName, 'HydroDyn Driver input file header line 2', ErrStat2, ErrMsg2, UnEchoLocal )
      if (Failed()) return

         ! Echo Input Files. Note this line is prevented from being echoed by the ReadVar routine.
      CALL ReadVar ( UnIn, FileName, InitInp%Echo, 'Echo', 'Echo the input file data', ErrStat2, ErrMsg2, UnEchoLocal )
      if (Failed()) return

      
   END IF
   !-------------------------------------------------------------------------------------------------
   ! Environmental conditions section
   !-------------------------------------------------------------------------------------------------

      ! Header
   CALL ReadCom( UnIn, FileName, 'Environmental conditions header', ErrStat2, ErrMsg2, UnEchoLocal )
   if (Failed()) return

      ! Gravity - Gravity.
   CALL ReadVar ( UnIn, FileName, InitInp%Gravity, 'Gravity', 'Gravity', ErrStat2, ErrMsg2, UnEchoLocal )
   if (Failed()) return

      ! WtrDens - Water density.
   CALL ReadVar ( UnIn, FileName, InitInp%WtrDens, 'WtrDens', 'Water density', ErrStat2, ErrMsg2, UnEchoLocal )
   if (Failed()) return

      ! WtrDpth - Water depth.
   CALL ReadVar ( UnIn, FileName, InitInp%WtrDpth, 'WtrDpth', 'Water depth', ErrStat2, ErrMsg2, UnEchoLocal )
   if (Failed()) return

      ! MSL2SWL - Offset between still-water level and mean sea level.
   CALL ReadVar ( UnIn, FileName, InitInp%MSL2SWL, 'MSL2SWL', 'Offset between still-water level and mean sea level', ErrStat2, ErrMsg2, UnEchoLocal )
   if (Failed()) return
   
   !-------------------------------------------------------------------------------------------------
   ! HYDRODYN section
   !-------------------------------------------------------------------------------------------------

      ! Header
   CALL ReadCom( UnIn, FileName, 'HYDRODYN header', ErrStat2, ErrMsg2, UnEchoLocal )
   if (Failed()) return
   
      ! HDInputFile
   CALL ReadVar ( UnIn, FileName, InitInp%HDInputFile, 'HDInputFile', 'HydroDyn input filename', ErrStat2, ErrMsg2, UnEchoLocal )
   if (Failed()) return
   IF ( PathIsRelative( InitInp%HDInputFile ) ) InitInp%HDInputFile = TRIM(PriPath)//TRIM(InitInp%HDInputFile)

       ! SeaStInputFile
   CALL ReadVar ( UnIn, FileName, InitInp%SeaStateInputFile, 'SeaStateInputFile', 'SeaState input filename', ErrStat2, ErrMsg2, UnEchoLocal )
   if (Failed()) return
   IF ( PathIsRelative( InitInp%SeaStateInputFile ) ) InitInp%SeaStateInputFile = TRIM(PriPath)//TRIM(InitInp%SeaStateInputFile)

      ! OutRootName
   CALL ReadVar ( UnIn, FileName, InitInp%OutRootName, 'OutRootName', 'HydroDyn output root filename', ErrStat2, ErrMsg2, UnEchoLocal )
   if (Failed()) return
   IF ( PathIsRelative( InitInp%OutRootName ) ) InitInp%OutRootName = TRIM(PriPath)//TRIM(InitInp%OutRootName)

       ! Linearize
   CALL ReadVar ( UnIn, FileName, InitInp%Linearize, 'Linearize', 'Linearize parameter', ErrStat2, ErrMsg2, UnEchoLocal )
   if (Failed()) return
  
      ! NSteps
   CALL ReadVar ( UnIn, FileName, InitInp%NSteps, 'NSteps', 'Number of time steps in the HydroDyn simulation', ErrStat2, ErrMsg2, UnEchoLocal )
   if (Failed()) return
   
      ! TimeInterval   
   CALL ReadVar ( UnIn, FileName, InitInp%TimeInterval, 'TimeInterval', 'Time interval for any HydroDyn inputs', ErrStat2, ErrMsg2, UnEchoLocal )
   if (Failed()) return
   
   
   !-------------------------------------------------------------------------------------------------
   ! PRP INPUTS section
   !-------------------------------------------------------------------------------------------------

      ! Header
   CALL ReadCom( UnIn, FileName, 'PRP INPUTS header', ErrStat2, ErrMsg2, UnEchoLocal )
   if (Failed()) return
   
      ! PRPInputsMod      
   CALL ReadVar ( UnIn, FileName, InitInp%PRPInputsMod, 'PRPInputsMod', 'Model for the PRP (principal reference point) inputs', ErrStat2, ErrMsg2, UnEchoLocal )
   if (Failed()) return
   
      ! PRPInputsFile      
   CALL ReadVar ( UnIn, FileName, InitInp%PRPInputsFile, 'PRPInputsFile', 'Filename for the PRP HydroDyn inputs', ErrStat2, ErrMsg2, UnEchoLocal )
   if (Failed()) return
   IF ( PathIsRelative( InitInp%PRPInputsFile ) ) InitInp%PRPInputsFile = TRIM(PriPath)//TRIM(InitInp%PRPInputsFile)
   
   
   !-------------------------------------------------------------------------------------------------
   ! PRP STEADY STATE INPUTS section
   !-------------------------------------------------------------------------------------------------

      ! Header
   CALL ReadCom( UnIn, FileName, 'PRP STEADY STATE INPUTS header', ErrStat2, ErrMsg2, UnEchoLocal )
   if (Failed()) return

      ! uPRPInSteady
   CALL ReadAry ( UnIn, FileName, InitInp%uPRPInSteady, 6, 'uPRPInSteady', 'PRP Steady-state displacements and rotations.', ErrStat2,  ErrMsg2, UnEchoLocal)
   if (Failed()) return
   
      ! uDotPRPInSteady
   CALL ReadAry ( UnIn, FileName, InitInp%uDotPRPInSteady, 6, 'uDotPRPInSteady', 'PRP Steady-state translational and rotational velocities.', ErrStat2,  ErrMsg2, UnEchoLocal)
   if (Failed()) return
      
      ! uDotDotPRPInSteady
   CALL ReadAry ( UnIn, FileName, InitInp%uDotDotPRPInSteady, 6, 'uDotDotPRPInSteady', 'PRP Steady-state translational and rotational accelerations.', ErrStat2,  ErrMsg2, UnEchoLocal)
   if (Failed()) return

      
   IF ( InitInp%PRPInputsMod /= 1 ) THEN
      InitInp%uPRPInSteady       = 0.0
      InitInp%uDotPRPInSteady    = 0.0
      InitInp%uDotDotPRPInSteady = 0.0
   END IF


   CALL cleanup()
   
CONTAINS

   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName) 
      Failed =  ErrStat >= AbortErrLev
      if (Failed) call Cleanup()
        
   end function Failed

   subroutine Cleanup()
      CLOSE( UnIn )
      IF ( UnEchoLocal > 0 ) CLOSE( UnEchoLocal )
   end subroutine Cleanup
   
END SUBROUTINE ReadDriverInputFile

!----------------------------------------------------------------------------------------------------------------------------------

END MODULE HydroDynDriverSubs

