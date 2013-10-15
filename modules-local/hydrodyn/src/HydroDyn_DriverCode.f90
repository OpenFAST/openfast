!**********************************************************************************************************************************
! HydroDyn_DriverCode: This code tests the template modules
!..................................................................................................................................
! LICENSING
! Copyright (C) 2012  National Renewable Energy Laboratory
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

PROGRAM HydroDynDriver

   USE NWTC_Library
   USE HydroDyn
   USE HydroDyn_Types
   USE HydroDyn_Output
   USE ModMesh_Types
   
   IMPLICIT NONE
   
   TYPE HD_Drvr_InitInput
      LOGICAL         :: Echo
      REAL(ReKi)      :: Gravity
      CHARACTER(1024) :: HDInputFile
      CHARACTER(1024) :: OutRootName
      INTEGER         :: NSteps
      REAL(DbKi)      :: TimeInterval
      INTEGER         :: WAMITInputsMod
      CHARACTER(1024) :: WAMITInputsFile
      REAL(ReKi)      :: uWAMITInSteady(6)
      REAL(ReKi)      :: uDotWAMITInSteady(6)
      REAL(ReKi)      :: uDotDotWAMITInSteady(6)
   END TYPE HD_Drvr_InitInput
   
! -----------------------------------------------------------------------------------   
! NOTE:  this module and the ModMesh.f90 modules must use the Fortran compiler flag:  
!        /fpp                  because of they both have preprocessor statements
! ----------------------------------------------------------------------------------- 



   INTEGER(IntKi), PARAMETER                           :: NumInp = 1           ! Number of inputs sent to HydroDyn_UpdateStates
   
      ! Program variables

   REAL(DbKi)                                          :: Time                 ! Variable for storing time, in seconds
  
   REAL(DbKi)                                          :: InputTime(NumInp)    ! Variable for storing time associated with inputs, in seconds
   
   INTEGER(B1Ki), ALLOCATABLE                          :: SaveAry(:)           ! Array to store packed data structure

   TYPE(HydroDyn_InitInputType)                        :: InitInData           ! Input data for initialization
   TYPE(HydroDyn_InitOutputType)                       :: InitOutData          ! Output data from initialization

   TYPE(HydroDyn_ContinuousStateType)                  :: x                    ! Continuous states
   TYPE(HydroDyn_ContinuousStateType)                  :: x_new                ! Continuous states at updated time
   TYPE(HydroDyn_DiscreteStateType)                    :: xd                   ! Discrete states
   TYPE(HydroDyn_DiscreteStateType)                    :: xd_new               ! Discrete states at updated time
   TYPE(HydroDyn_ConstraintStateType)                  :: z                    ! Constraint states
   TYPE(HydroDyn_ConstraintStateType)                  :: z_residual           ! Residual of the constraint state equations (Z)
   TYPE(HydroDyn_OtherStateType)                       :: OtherState           ! Other/optimization states

   TYPE(HydroDyn_ParameterType)                        :: p                    ! Parameters
   !TYPE(HydroDyn_InputType)                           :: u                    ! System inputs [OLD STYLE]
   TYPE(HydroDyn_InputType)                            :: u(NumInp)            ! System inputs
   TYPE(HydroDyn_OutputType)                           :: y                    ! System outputs

   TYPE(HydroDyn_ContinuousStateType)                  :: dxdt                 ! First time derivatives of the continuous states


   INTEGER(IntKi)                                     :: UnWAMITInp            ! WAMIT Inputs file identifier
   INTEGER(IntKi)                                     :: UnHD_Out              ! Output file identifier
   REAL(ReKi), ALLOCATABLE                            :: WAMITin(:,:)          ! Variable for storing time, forces, and body velocities, in m/s or rad/s for WAMIT
   INTEGER(IntKi)                                     :: I                    ! Generic loop counter
   INTEGER(IntKi)                                     :: J                    ! Generic loop counter
   INTEGER(IntKi)                                     :: n                    ! Loop counter (for time step)
   INTEGER(IntKi)                                     :: ErrStat              ! Status of error message
   CHARACTER(1024)                                    :: ErrMsg               ! Error message if ErrStat /= ErrID_None
   REAL(ReKi)                                         :: dcm (3,3)            ! The resulting transformation matrix from X to x, (-).
   CHARACTER(1024)                                    :: drvrFilename         ! Filename and path for the driver input file.  This is passed in as a command line argument when running the Driver exe.
   TYPE(HD_Drvr_InitInput)                            :: drvrInitInp          ! Initialization data for the driver program
   
   ! For testing
   LOGICAL                                            :: DoTight = .FALSE.
   REAL(DbKi)                                         :: maxAngle             ! For debugging, see what the largest rotational angle input is for the simulation
   CHARACTER(10)                                      :: AngleMsg             ! For debugging, a string version of the largest rotation input
   
   !...............................................................................................................................
   ! Routines called in initialization
   !...............................................................................................................................

   
   
   ! TODO: Need to think some more about how to pass DRIVER-level initialization data to the HydroDyn module because if UseInputFile = .FALSE.
   !       then the input processing code will still be querying the *Chr input data to look for the use of the 'DEFAULT' string and to set that
   !       data to the driver's version instead of using a module-specific version.  
   !       Currently, these variables are:
   !          InitInp%Waves%WavePkShpChr
   !          InitInp%Current%CurrSSDirChr
   !          InitInp%PtfmSgFChr
   !          InitInp%PtfmSwFChr
   !          InitInp%PtfmHvFChr
   !          InitInp%PtfmRFChr
   !          InitInp%PtfmPFChr
   !          InitInp%PtfmYFChr
   !          InitInp%Morison%InpMembers(k)%FillDensChr
   !          
   !          
      ! Initialize the library which handle file echos and WrScr, for example
   call nwtc_init()
   
   
   IF ( command_argument_count() > 1 ) THEN
      CALL print_help()
      STOP
   END IF
  
   
      ! Parse the driver input file and run the simulation based on that file
      
   IF ( command_argument_count() == 1 ) THEN
      
      CALL get_command_argument(1, drvrFilename)
      CALL ReadDriverInputFile( drvrFilename, drvrInitInp, ErrStat, ErrMsg )
      IF ( ErrStat /= 0 ) THEN
         CALL WrScr( ErrMsg )
         STOP
      END IF
      InitInData%Gravity      = drvrInitInp%Gravity
      InitInData%UseInputFile = .TRUE. 
      InitInData%InputFile    = drvrInitInp%HDInputFile
      InitInData%OutRootName  = drvrInitInp%OutRootName
      
   END IF
  
 



  
  

!-------------------------------------------------------------------------------------
!       Begin Simulation Setup
!-------------------------------------------------------------------------------------
   

   IF ( drvrInitInp%WAMITInputsMod == 2 ) THEN
      
         ! Open the WAMIT inputs data file
      CALL GetNewUnit( UnWAMITInp ) 
      CALL OpenFInpFile ( UnWAMITInp, drvrInitInp%WAMITInputsFile, ErrStat   )  ! Open WAMIT inputs file.
      
      IF ( ErrStat /= 0 ) THEN
         ErrMsg = 'HydroDyn Driver file not found'
         CALL WrScr( ErrMsg )
         STOP
      END IF
      
      ALLOCATE ( WAMITin(drvrInitInp%NSteps, 19), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating space for WAMITin array.'
         CALL WrScr( ErrMsg )
         CLOSE( UnWAMITInp )
         STOP
      END IF 
      
      DO n = 1,drvrInitInp%NSteps
         READ (UnWAMITInp,*,IOSTAT=ErrStat) (WAMITin (n,J), J=1,19)
            
            IF ( ErrStat /= 0 ) THEN
               ErrMsg = 'File not found'
               CALL WrScr( ErrMsg )
               STOP
            END IF 
      END DO  
      
         ! Close the inputs file 
      CLOSE ( UnWAMITInp ) 
   END IF
   
     
   
         ! Initialize the module

   CALL HydroDyn_Init( InitInData, u(1), p,  x, xd, z, OtherState, y, drvrInitInp%TimeInterval, InitOutData, ErrStat, ErrMsg )
   
   IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
      CALL WrScr( ErrMsg )
      IF ( ErrStat > ErrID_Warn ) THEN
         CALL HydroDyn_End( u(1), p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
         IF ( ErrStat /= ErrID_None ) THEN
            CALL WrScr( ErrMsg )     
         END IF
         STOP
      END IF
   END IF

      ! Destroy initialization data

   CALL HydroDyn_DestroyInitInput(  InitInData,  ErrStat, ErrMsg )
   CALL HydroDyn_DestroyInitOutput( InitOutData, ErrStat, ErrMsg )
   
   !...............................................................................................................................
   ! Routines called in loose coupling -- the glue code may implement this in various ways
   !...............................................................................................................................


   DO n = 1,drvrInitInp%NSteps

      Time = (n-1) * drvrInitInp%TimeInterval
      InputTime(1) = Time
      
         ! Modify u (likely from the outputs of another module or a set of test conditions) here:
         
      IF ( u(1)%WAMIT%Mesh%Initialized ) THEN 
         
         IF ( drvrInitInp%WAMITInputsMod == 2 ) THEN
            
            
            
            u(1)%WAMIT%Mesh%TranslationDisp(:,1)   = WAMITin(n,2:4) 
            
            
               ! Compute direction cosine matrix from the rotation angles
               
            IF ( abs(WAMITin(n,5)) > maxAngle ) maxAngle = abs(WAMITin(n,5))
            IF ( abs(WAMITin(n,6)) > maxAngle ) maxAngle = abs(WAMITin(n,6))
            IF ( abs(WAMITin(n,7)) > maxAngle ) maxAngle = abs(WAMITin(n,7))
            
            CALL SmllRotTrans( 'InputRotation', REAL(WAMITin(n,5)), REAL(WAMITin(n,6)), REAL(WAMITin(n,7)), dcm, 'Junk', ErrStat, ErrMsg )            
            u(1)%WAMIT%Mesh%Orientation(:,:,1)     = dcm 
            
            
            u(1)%WAMIT%Mesh%TranslationVel(:,1)    = WAMITin(n,8:10)  
            u(1)%WAMIT%Mesh%RotationVel(:,1)       = WAMITin(n,11:13) 
            u(1)%WAMIT%Mesh%TranslationAcc(:,1)    = WAMITin(n,14:16)  
            u(1)%WAMIT%Mesh%RotationAcc(:,1)       = WAMITin(n,17:19) 
            
         ELSE
            
            u(1)%WAMIT%Mesh%TranslationDisp(:,1)   = drvrInitInp%uWAMITInSteady(1:3) 
            
            
               ! Compute direction cosine matrix from the rotation angles
            CALL SmllRotTrans( 'InputRotation', REAL(drvrInitInp%uWAMITInSteady(4)), REAL(drvrInitInp%uWAMITInSteady(5)), REAL(drvrInitInp%uWAMITInSteady(6)), dcm, 'Junk', ErrStat, ErrMsg )            
            u(1)%WAMIT%Mesh%Orientation(:,:,1)     = dcm
            
            u(1)%WAMIT%Mesh%TranslationVel(:,1)    = drvrInitInp%uDotWAMITInSteady(1:3)  
            u(1)%WAMIT%Mesh%RotationVel(:,1)       = drvrInitInp%uDotWAMITInSteady(4:6) 
            u(1)%WAMIT%Mesh%TranslationAcc(:,1)    = drvrInitInp%uDotDotWAMITInSteady(1:3)  
            u(1)%WAMIT%Mesh%RotationAcc(:,1)       = drvrInitInp%uDotDotWAMITInSteady(4:6) 
            
         END IF
         
      END IF
      
      !IF ( u(1)%Morison%LumpedMesh%Initialized ) THEN   
      !u(1)%Morison%LumpedMesh%TranslationVel(:,J) = l_ud_in(n,J,1:3)
      !u(1)%Morison%LumpedMesh%RotationVel(:,J) = l_ud_in(n,J,4:6)
      !END IF
      !IF ( u(1)%Morison%DistribMesh%Initialized ) THEN  
      !u(1)%Morison%DistribMesh%TranslationVel(:,J) = d_ud_in(n,J,1:3)
      !u(1)%Morison%DistribMesh%RotationVel(:,J) = d_ud_in(n,J,4:6)
      !END IF
      
      
     

         ! Calculate outputs at n

      CALL HydroDyn_CalcOutput( Time, u(1), p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
         CALL WrScr( ErrMsg )
         IF ( ErrStat > ErrID_Warn ) THEN
            CALL HydroDyn_End( u(1), p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
            IF ( ErrStat /= ErrID_None ) THEN
               CALL WrScr( ErrMsg )      
            END IF
            STOP
         END IF
      END IF

      
      
         ! Get state variables at next step: INPUT at step n, OUTPUT at step n + 1

      CALL HydroDyn_UpdateStates( Time, n, u, InputTime, p, x, xd, z, OtherState, ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
         CALL WrScr( ErrMsg )
         IF ( ErrStat > ErrID_Warn ) THEN
            CALL HydroDyn_End( u(1), p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
            IF ( ErrStat /= ErrID_None ) THEN
               CALL WrScr( ErrMsg )      
            END IF
            STOP
         END IF
      END IF     
      
      
      

      ! Write output to a file which is managed by the driver program and not the individual modules
      ! TODO
      
   END DO

   ! For Debug:  TODO  Remove when no longer needed GJH 9/30/2013 
   !Write(AngleMsg, '(F10.4)') ,maxAngle
   !CALL WrScr('The largest input rotation angle was: '//AngleMsg)
   
   CALL HydroDyn_DestroyDiscState( xd_new, ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
      CALL WrScr( ErrMsg )
      IF ( ErrStat > ErrID_Warn ) THEN
         CALL HydroDyn_End( u(1), p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
         IF ( ErrStat /= ErrID_None ) THEN
            CALL WrScr( ErrMsg )      
         END IF
         STOP
      END IF
   END IF

   CALL HydroDyn_DestroyContState( x_new, ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
      CALL WrScr( ErrMsg )
      IF ( ErrStat > ErrID_Warn ) THEN
         CALL HydroDyn_End( u(1), p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
         IF ( ErrStat /= ErrID_None ) THEN
            CALL WrScr( ErrMsg )      
         END IF
         STOP
      END IF
   END IF


! For now, finish here.

 
   !...............................................................................................................................
   ! Routine to terminate program execution (again)
   !...............................................................................................................................

   CALL HydroDyn_End( u(1), p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) THEN
      CALL WrScr( ErrMsg )
      ! TODO Add additional error handling and terminate
   END IF


CONTAINS

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


SUBROUTINE ReadDriverInputFile( inputFile, InitInp, ErrStat, ErrMsg )

   CHARACTER(1024),               INTENT( IN    )   :: inputFile
   TYPE(HD_Drvr_InitInput),       INTENT(   OUT )   :: InitInp
   INTEGER,                       INTENT(   OUT )   :: ErrStat              ! returns a non-zero value when an error occurs  
   CHARACTER(*),                  INTENT(   OUT )   :: ErrMsg               ! Error message if ErrStat /= ErrID_None
   
      ! Local variables  
         
   INTEGER                                          :: I                    ! generic integer for counting
   INTEGER                                          :: J                    ! generic integer for counting
   CHARACTER(   2)                                  :: strI                 ! string version of the loop counter

   INTEGER                                          :: UnIn                 ! Unit number for the input file
   INTEGER                                          :: UnEchoLocal          ! The local unit number for this module's echo file
   CHARACTER(1024)                                  :: EchoFile             ! Name of HydroDyn echo file  
   CHARACTER(1024)                                  :: Line                 ! String to temporarially hold value of read line   
   CHARACTER(1024)                                  :: TmpPath              ! Temporary storage for relative path name
   CHARACTER(1024)                                  :: TmpFmt               ! Temporary storage for format statement
   CHARACTER(1024)                                  :: FileName             ! Name of HydroDyn input file  
   
   
      ! Initialize the echo file unit to -1 which is the default to prevent echoing, we will alter this based on user input
   UnEchoLocal = -1
   
   FileName = TRIM(inputFile)
   
   CALL GetNewUnit( UnIn )   
   CALL OpenFInpFile( UnIn, FileName, ErrStat )
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to open HydroDyn Driver input file: '//FileName
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF

   
   CALL WrScr( 'Opening HydroDyn Driver input file:  '//FileName )
   
   
   !-------------------------------------------------------------------------------------------------
   ! File header
   !-------------------------------------------------------------------------------------------------
   
   CALL ReadCom( UnIn, FileName, 'HydroDyn Driver input file header line 1', ErrStat )
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read HydroDyn Driver input file header line 1.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF


   CALL ReadCom( UnIn, FileName, 'HydroDyn Driver input file header line 2', ErrStat )
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read HydroDyn Driver input file header line 2.'
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
      
      CALL ReadCom( UnIn, FileName, 'HydroDyn Driver input file header line 1', ErrStat, ErrMsg, UnEchoLocal )
   
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Failed to read HydroDyn Driver input file header line 1.'
         ErrStat = ErrID_Fatal
         CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      END IF


      CALL ReadCom( UnIn, FileName, 'HydroDyn Driver input file header line 2', ErrStat, ErrMsg, UnEchoLocal )
   
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Failed to read HydroDyn Driver input file header line 2.'
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


      ! Gravity - Gravity.
      
   CALL ReadVar ( UnIn, FileName, InitInp%Gravity, 'Gravity', 'Gravity', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Gravity parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF

   
   !-------------------------------------------------------------------------------------------------
   ! HYDRODYN section
   !-------------------------------------------------------------------------------------------------

      ! Header
      
   CALL ReadCom( UnIn, FileName, 'HYDRODYN header', ErrStat, ErrMsg, UnEchoLocal )
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Comment line.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
   
   
      ! HDInputFile
      
   CALL ReadVar ( UnIn, FileName, InitInp%HDInputFile, 'HDInputFile', &
                                    'HydroDyn input filename', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read HDInputFile parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF 
   
   
      ! OutRootName
   
   CALL ReadVar ( UnIn, FileName, InitInp%OutRootName, 'OutRootName', &
                                    'HydroDyn output root filename', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read OutRootName parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF   
     
   
      ! NSteps
   
   CALL ReadVar ( UnIn, FileName, InitInp%NSteps, 'NSteps', &
                                    'Number of time steps in the HydroDyn simulation', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read NSteps parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF   
 
   
      ! TimeInterval   
   
   CALL ReadVar ( UnIn, FileName, InitInp%TimeInterval, 'TimeInterval', &
                                    'Time interval for any HydroDyn inputs', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read TimeInterval parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF   
   
   
   !-------------------------------------------------------------------------------------------------
   ! WAMIT INPUTS section
   !-------------------------------------------------------------------------------------------------

      ! Header
      
   CALL ReadCom( UnIn, FileName, 'WAMIT INPUTS header', ErrStat, ErrMsg, UnEchoLocal )
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Comment line.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
 
   
   
      ! WAMITInputsMod      
       
   CALL ReadVar ( UnIn, FileName, InitInp%WAMITInputsMod, 'WAMITInputsMod', &
                                    'Model for the WAMIT inputs', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read WAMITInputsMod parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF   
   
   
      ! WAMITInputsFile      
       
   CALL ReadVar ( UnIn, FileName, InitInp%WAMITInputsFile, 'WAMITInputsFile', &
                                    'Filename for the HydroDyn inputs', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read WAMITInputsFile parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF   
   
   
   !-------------------------------------------------------------------------------------------------
   ! WAMIT STEADY STATE INPUTS section
   !-------------------------------------------------------------------------------------------------

      ! Header
      
   CALL ReadCom( UnIn, FileName, 'WAMIT STEADY STATE INPUTS header', ErrStat, ErrMsg, UnEchoLocal )
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Comment line.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
   
   IF ( InitInp%WAMITInputsMod == 1 ) THEN
   
         ! uWAMITInSteady
         
      CALL ReadAry ( UnIn, FileName, InitInp%uWAMITInSteady, 6, 'uWAMITInSteady', &
                           'WAMIT Steady-state displacements and rotations.', ErrStat,  ErrMsg, UnEchoLocal)         
       
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Failed to read uWAMITInSteady parameter.'
         ErrStat = ErrID_Fatal
         CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      END IF
   
   
         ! uDotWAMITInSteady
         
      CALL ReadAry ( UnIn, FileName, InitInp%uDotWAMITInSteady, 6, 'uDotWAMITInSteady', &
                           'WAMIT Steady-state translational and rotational velocities.', ErrStat,  ErrMsg, UnEchoLocal)         
       
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Failed to read uDotWAMITInSteady parameter.'
         ErrStat = ErrID_Fatal
         CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      END IF
      
      
         ! uDotDotWAMITInSteady
         
      CALL ReadAry ( UnIn, FileName, InitInp%uDotDotWAMITInSteady, 6, 'uDotDotWAMITInSteady', &
                           'WAMIT Steady-state translational and rotational accelerations.', ErrStat,  ErrMsg, UnEchoLocal)         
       
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Failed to read uDotDotWAMITInSteady parameter.'
         ErrStat = ErrID_Fatal
         CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      END IF
      
   ELSE
      InitInp%uWAMITInSteady       = 0.0
      InitInp%uDotWAMITInSteady    = 0.0
      InitInp%uDotDotWAMITInSteady = 0.0
   END IF
   
   
   CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
   CLOSE( UnIn )
   
END SUBROUTINE ReadDriverInputFile

subroutine print_help()
    print '(a)', 'usage: '
    print '(a)', ''
    print '(a)', 'HydroDyn.exe driverfilename'
    print '(a)', ''
    print '(a)', 'Where driverfilename is the name of the HydroDyn driver input file.'
    print '(a)', ''

end subroutine print_help

END PROGRAM HydroDynDriver

