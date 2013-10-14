!**********************************************************************************************************************************
! SubDyn_DriverCode: This code tests the SubDyn modules
!..................................................................................................................................
! LICENSING
! Copyright (C) 2013  National Renewable Energy Laboratory
!
!    This file is part of SubDyn.
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
PROGRAM TestSubDyn

   USE NWTC_Library
   USE SubDyn
   USE SubDyn_Types
   USE SubDyn_Output

   IMPLICIT NONE

   INTEGER(IntKi), PARAMETER                          :: NumInp = 1           ! Number of inputs sent to SD_UpdateStates
   
   
   TYPE SD_Drvr_InitInput
      LOGICAL         :: Echo
      REAL(ReKi)      :: Gravity
      CHARACTER(1024) :: SDInputFile
      CHARACTER(1024) :: OutRootName
      INTEGER         :: NSteps
      REAL(DbKi)      :: TimeInterval
      REAL(ReKi)      :: TP_RefPoint(3)
      REAL(ReKi)      :: SubRotateZ
      INTEGER         :: InputsMod
      CHARACTER(1024) :: InputsFile
      REAL(ReKi)      :: uTPInSteady(6)
      REAL(ReKi)      :: uDotTPInSteady(6)
      REAL(ReKi)      :: uDotDotTPInSteady(6)
   END TYPE SD_Drvr_InitInput
   
   
      ! Program variables

   REAL(DbKi)                                         :: Time                 ! Variable for storing time, in seconds
   REAL(DbKi)                                         :: TimeInterval         ! Interval between time steps, in seconds
   REAL(DbKi)                                         :: InputTime(NumInp)    ! Variable for storing time associated with inputs, in seconds
   
   TYPE(SD_InitInputType)                        :: InitInData           ! Input data for initialization
   TYPE(SD_InitOutputType)                       :: InitOutData          ! Output data from initialization

   TYPE(SD_ContinuousStateType)                  :: x                    ! Continuous states
   TYPE(SD_DiscreteStateType)                    :: xd                   ! Discrete states
   TYPE(SD_ConstraintStateType)                  :: z, z_next                    ! Constraint states
   TYPE(SD_ConstraintStateType)                  :: Z_residual           ! Residual of the constraint state functions (Z)
   TYPE(SD_OtherStateType)                       :: OtherState           ! Other/optimization states

   TYPE(SD_ParameterType)                        :: p                    ! Parameters
   TYPE(SD_InputType)                            :: u(NumInp)            ! System inputs
   TYPE(SD_OutputType)                           :: y                    ! System outputs

   TYPE(SD_ContinuousStateType)                  :: dxdt                 ! First time derivatives of the continuous states

   TYPE(SD_PartialOutputPInputType)              :: dYdu                 ! Partial derivatives of the output functions
                                                                              !  (Y) with respect to the inputs (u)
   TYPE(SD_PartialContStatePInputType)           :: dXdu                 ! Partial derivatives of the continuous state
                                                                              !  functions (X) with respect to the inputs (u)
   TYPE(SD_PartialDiscStatePInputType)           :: dXddu                ! Partial derivatives of the discrete state
                                                                              !  functions (Xd) with respect to the inputs (u)
   TYPE(SD_PartialConstrStatePInputType)         :: dZdu                 ! Partial derivatives of the constraint state
                                                                              !  functions (Z) with respect to the inputs (u)
   TYPE(SD_PartialOutputPContStateType)          :: dYdx                 ! Partial derivatives of the output functions (Y)
                                                                              !  with respect to the continuous states (x)
   TYPE(SD_PartialContStatePContStateType)       :: dXdx                 ! Partial derivatives of the continuous state funct-
                                                                              !  ions (X) with respect to the continuous states (x)
   TYPE(SD_PartialDiscStatePContStateType)       :: dXddx                ! Partial derivatives of the discrete state funct-
                                                                              !  ions (Xd) with respect to continuous states (x)
   TYPE(SD_PartialConstrStatePContStateType)     :: dZdx                 ! Partial derivatives of the constraint state funct-
                                                                              !  ions (Z) with respect to the continuous states (x)
   TYPE(SD_PartialOutputPDiscStateType)          :: dYdxd                ! Partial derivatives of the output functions (Y)
                                                                              !  with respect to the discrete states (xd)
   TYPE(SD_PartialContStatePDiscStateType)       :: dXdxd                ! Partial derivatives of the continuous state funct-
                                                                              !  ions (X) with respect to the discrete states (xd)
   TYPE(SD_PartialDiscStatePDiscStateType)       :: dXddxd               ! Partial derivatives of the discrete state funct-
                                                                              !  ions (Xd) with respect to the discrete states (xd)
   TYPE(SD_PartialConstrStatePDiscStateType)     :: dZdxd                ! Partial derivatives of the constraint state funct-
                                                                              !  ions (Z) with respect to the discrete states (xd)
   TYPE(SD_PartialOutputPConstrStateType)        :: dYdz                 ! Partial derivatives of the output functions (Y)
                                                                              !  with respect to the constraint states (z)
   TYPE(SD_PartialContStatePConstrStateType)     :: dXdz                 ! Partial derivatives of the continuous state funct-
                                                                              !  ions (X) with respect to the constraint states (z)
   TYPE(SD_PartialDiscStatePConstrStateType)     :: dXddz                ! Partial derivatives of the discrete state funct-
                                                                              !  ions (Xd) with respect to constraint states (z)
   TYPE(SD_PartialConstrStatePConstrStateType)   :: dZdz                 ! Partial derivatives of the constraint state funct-
                                                                              !  ions (Z) with respect to the constraint states (z)


   INTEGER(IntKi)                                     :: n                    ! Loop counter (for time step)
   INTEGER(IntKi)                                     :: ErrStat, ErrStat1, ErrStat2, ErrStat3              ! Status of error message
   CHARACTER(1024)                                    :: ErrMsg, ErrMsg1, ErrMsg2, ErrMsg3              ! Error message if ErrStat /= ErrID_None


   REAL(ReKi), ALLOCATABLE                            :: Re_SaveAry  (:)      ! Array to store reals in packed data structure
   REAL(DbKi), ALLOCATABLE                            :: Db_SaveAry  (:)      ! Array to store doubles in packed data structure
   INTEGER(IntKi), ALLOCATABLE                        :: Int_SaveAry (:)      ! Array to store integers in packed data structure
   CHARACTER(1024)                                    :: drvrFilename         ! Filename and path for the driver input file.  This is passed in as a command line argument when running the Driver exe.
   TYPE(SD_Drvr_InitInput)                            :: drvrInitInp          ! Initialization data for the driver program
   INTEGER(IntKi)                                     :: UnInp                !  Inputs file identifier
   INTEGER(IntKi)                                     :: UnSD_Out             ! Output file identifier
   REAL(ReKi), ALLOCATABLE                            :: SDin(:,:)            ! Variable for storing time, forces, and body velocities, in m/s or rad/s for SubDyn inputs
   INTEGER(IntKi)                                     :: J                    ! Generic loop counter
   REAL(ReKi)                                         :: dcm (3,3)            ! The resulting transformation matrix from X to x, (-).
   REAL(DbKi)                                         :: maxAngle             ! For debugging, see what the largest rotational angle input is for the simulation
   CHARACTER(10)                                      :: AngleMsg             ! For debugging, a string version of the largest rotation input
   
      ! Other/Misc variables
REAL(DbKi)                            :: TiLstPrn                                 ! The time of the last print
REAL(DbKi)                            :: TMax
REAL(DbKi)                            :: OutTime                                  ! Used to determine if output should be generated at this simulation time
REAL(ReKi)                            :: PrevClockTime                            ! Clock time at start of simulation in seconds
REAL                                  :: UsrTime1                                 ! User CPU time for simulation initialization
INTEGER                               :: StrtTime (8)                             ! Start time of simulation
   !...............................................................................................................................
   ! Routines called in initialization
   !...............................................................................................................................

        ! Get the current time
   CALL DATE_AND_TIME ( Values=StrtTime )                               ! Let's time the whole simulation
   CALL CPU_TIME ( UsrTime1 )                                           ! Initial time (this zeros the start time when used as a MATLAB function)
   PrevClockTime = TimeValues2Seconds( StrtTime )                       ! We'll use this time for the SimStats routine
   TiLstPrn      = 0.0_DbKi                                             ! The first value of ZTime, used to write simulation stats to screen (s)
  
   
         ! Initialize the NWTC Subroutine Library

   CALL NWTC_Init( )
   
   IF ( command_argument_count() > 1 ) CALL print_help()

  
      ! Parse the driver input file and run the simulation based on that file
      
   IF ( command_argument_count() == 1 ) THEN
      
      CALL get_command_argument(1, drvrFilename)
      CALL ReadDriverInputFile( drvrFilename, drvrInitInp, ErrStat, ErrMsg )
      IF ( ErrStat /= 0 ) THEN
         CALL WrScr( ErrMsg )
         STOP
      END IF
      InitInData%g            = drvrInitInp%Gravity
      !InitInData%UseInputFile = .TRUE. 
      InitInData%SDInputFile  = drvrInitInp%SDInputFile
      InitInData%RootName  = drvrInitInp%OutRootName
      InitInData%TP_RefPoint  = drvrInitInp%TP_RefPoint
      InitInData%SubRotateZ   = drvrInitInp%SubRotateZ
      TimeInterval            = drvrInitInp%TimeInterval
   ELSE
         ! Called without a driver input file!

      ! InitInData%SDInputFile = '..\BeamFEM\IOFiles\TestBeam2.txt'
      InitInData%SDInputFile = '..\MergedSubDyn\IOFiles\TestBeam3.txt'
      ! InitInData%SDInputFile = '..\BeamFEM\IOFiles\TestFrame.txt'
      InitInData%g =  9.80665
      !InitInData%TP_RefPoint = (/0.0, 0.0, 100.0/)  !testbeam2
      InitInData%TP_RefPoint = (/50.0, 0.0, 50.0/)  !testbeam3
      InitInData%SubRotateZ   = 0.0
      !InitInData%TP_RefPoint = (/0.0, 0.0, 40.0/)  !testframe
         ! Set the driver's request for time interval here:
      TimeInterval = 0.001 ! Glue code's request for delta time (likely based on information from other modules)
   END IF
   
   
  TMax = TimeInterval * drvrInitInp%NSteps
   
         ! Initialize the module
   
   CALL SD_Init( InitInData, u(1), p,  x, xd, z, OtherState, y, TimeInterval, InitOutData, ErrStat1, ErrMsg1 )
   


       ! Read Input time series data from a file
      
   IF ( drvrInitInp%InputsMod == 2 ) THEN
      
         ! Open the  inputs data file
      CALL GetNewUnit( UnInp ) 
      CALL OpenFInpFile ( UnInp, drvrInitInp%InputsFile, ErrStat   )  ! Open  inputs file.
      
      IF ( ErrStat /= 0 ) THEN
         ErrMsg = 'SubDyn input timeseries file not found'
         CALL WrScr( ErrMsg )
         STOP
      END IF
      
      ALLOCATE ( SDin(drvrInitInp%NSteps, 13), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating space for SDin array.'
         CALL WrScr( ErrMsg )
         CLOSE( UnInp )
         STOP
      END IF 
      
      DO n = 1,drvrInitInp%NSteps
         READ (UnInp,*,IOSTAT=ErrStat) (SDin (n,J), J=1,13)
            
            IF ( ErrStat /= 0 ) THEN
               ErrMsg = 'File not found'
               CALL WrScr( ErrMsg )
               CLOSE ( UnInp ) 
               STOP
            END IF 
      END DO  
      
         ! Close the inputs file 
      CLOSE ( UnInp ) 
   END IF 
  
         ! Destroy initialization data

   CALL SD_DestroyInitInput(  InitInData,  ErrStat2, ErrMsg2 )
   CALL SD_DestroyInitOutput( InitOutData, ErrStat3, ErrMsg3 )

   
      ! Handle the initialization error after destroying the data structures
   
   IF ( ErrStat1 /= ErrID_None .OR. ErrStat2 /=0 .OR. ErrStat3 /= 0) THEN          ! Check if there was an error and do something about it if necessary
      CALL WrScr( ErrMsg1 )
      STOP
   END IF
   
   IF ( ErrStat2 /=0 .OR. ErrStat3 /= 0) THEN          ! Check if there was an error and do something about it if necessary
      CALL WrScr( 'Error destroying SubDyn intialization data' )
      STOP
   END IF

   !...............................................................................................................................
   ! Routines called in loose coupling -- the glue code may implement this in various ways
   !...............................................................................................................................

   ! Force the displacement of the interface node in the global Z direction to be the sag of the column under it's own weight

   !u(1)%UFL(3) = -0.001821207  !-0.001821235   !This is for testbeam.txt
   ! u(1)%UFL(3)=-12.958  !this is for testbeam3
    
  
   DO n = 0,drvrInitInp%NSteps

      Time = n*TimeInterval
      InputTime(1) = Time

         ! Modify u (likely from the outputs of another module or a set of test conditions) here:

      IF ( u(1)%TPMesh%Initialized ) THEN 
         
         ! For now, set all hydrodynamic load inputs to 0.0
         u(1)%LMesh%Force  (:,:) = 0.0
         u(1)%LMesh%Moment (:,:) = 0.0
         
         IF ( drvrInitInp%InputsMod == 2 ) THEN
            
            
            
            u(1)%TPMesh%TranslationDisp(:,1)   = SDin(n,2:4) 
            
            
               ! Compute direction cosine matrix from the rotation angles
               
            IF ( abs(SDin(n,5)) > maxAngle ) maxAngle = abs(SDin(n,5))
            IF ( abs(SDin(n,6)) > maxAngle ) maxAngle = abs(SDin(n,6))
            IF ( abs(SDin(n,7)) > maxAngle ) maxAngle = abs(SDin(n,7))
            
            CALL SmllRotTrans( 'InputRotation', REAL(SDin(n,5)), REAL(SDin(n,6)), REAL(SDin(n,7)), dcm, 'Junk', ErrStat, ErrMsg )            
            u(1)%TPMesh%Orientation(:,:,1)     = dcm 
            
            
            u(1)%TPMesh%TranslationVel(:,1)    = SDin(n,8:10)  
            u(1)%TPMesh%RotationVel(:,1)       = SDin(n,11:13) 
            
         ELSE
            
            u(1)%TPMesh%TranslationDisp(:,1)   = drvrInitInp%uTPInSteady(1:3) 
            
            
               ! Compute direction cosine matrix from the rotation angles
            CALL SmllRotTrans( 'InputRotation', REAL(drvrInitInp%uTPInSteady(4)), REAL(drvrInitInp%uTPInSteady(5)), REAL(drvrInitInp%uTPInSteady(6)), dcm, 'Junk', ErrStat, ErrMsg )            
            u(1)%TPMesh%Orientation(:,:,1)     = dcm
            
            u(1)%TPMesh%TranslationVel(:,1)    = drvrInitInp%uDotTPInSteady(1:3)  
            u(1)%TPMesh%RotationVel(:,1)       = drvrInitInp%uDotTPInSteady(4:6) 
            
            u(1)%TPMesh%TranslationAcc(:,1)    = drvrInitInp%uDotDotTPInSteady(1:3)  
            u(1)%TPMesh%RotationAcc(:,1)       = drvrInitInp%uDotDotTPInSteady(4:6) 
            
         END IF
         
      END IF   
         ! Calculate outputs at n
      
      CALL SD_CalcOutput( Time, u(1), p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
         CALL WrScr( ErrMsg )
      END IF

         
         ! Get state variables at next step: INPUT at step n, OUTPUT at step n + 1
                                  
      CALL SD_UpdateStates( Time, n, u,      InputTime,  p, x, xd, z, OtherState, ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
         CALL WrScr( ErrMsg )
      END IF     
      
   !.....................................................
      ! Display simulation status every SttsTime-seconds:
      !.....................................................

      IF ( Time - TiLstPrn >= 1 )  THEN

         CALL SimStatus( TiLstPrn, PrevClockTime, Time, TMax )

      ENDIF   
   END DO


   !...............................................................................................................................
   ! Routines called in tight coupling -- time marching only
   !...............................................................................................................................

!   DO n = 0,2
!
!      Time = n * TimeInterval   ! Note that the discrete states must be updated only at the TimeInterval defined in initialization
!
!         ! set inputs (u) here:
!!      u =
!
!         ! Update constraint states at Time
!
!      ! DO 
!
!   !      CALL SD_CalcConstrStateResidual( Time, u(1), p, x, xd, z, OtherState, Z_residual, ErrStat, ErrMsg )
!   !      IF ( ErrStat /= ErrID_None ) THEN      ! Check if there was an error and do something about it if necessary
!   !         CALL WrScr( ErrMsg )
!   !      END IF
!
!         ! z =
!
!      ! END DO
!
!
!
!         ! Calculate the outputs at Time
!
!      CALL SD_CalcOutput( Time, u(1), p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
!      IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
!         CALL WrScr( ErrMsg )
!      END IF
!
!
!         ! Calculate the continuous state derivatives at Time
!
!   !   CALL SD_CalcContStateDeriv( Time, u(1), p, x, xd, z, OtherState, dxdt, ErrStat, ErrMsg )
!   !   IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
!   !      CALL WrScr( ErrMsg )
!   !   END IF
!
!
!         ! Update the discrete state from step n to step n+1
!         ! Note that the discrete states must be updated only at the TimeInterval defined in initialization
!
!!      CALL SD_UpdateDiscState( Time, n, u(1), p, x, xd, z, OtherState, ErrStat, ErrMsg )
!!      IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
!!         CALL WrScr( ErrMsg )
!!      END IF
!
!
!         ! Driver should integrate (update) continuous states here:
!
!      !x = function of dxdt, x
!
!
!         ! Jacobians required:
!                              
!      CALL SD_JacobianPInput( Time, u(1), p, x, xd, z, OtherState, dYdu, dXdu, dXddu, dZdu, ErrStat, ErrMsg )
!      IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
!         CALL WrScr( ErrMsg )
!      END IF
!                                   
!      CALL SD_JacobianPConstrState( Time, u(1), p, x, xd, z, OtherState, dYdz,dXdz, dXddz, dZdz, ErrStat, ErrMsg  )
!      IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
!         CALL WrScr( ErrMsg )
!      END IF
!
!
!   END DO


      ! Destroy Z_residual and dxdt because they are not necessary anymore

   CALL SD_DestroyConstrState( Z_residual, ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) THEN   ! Check if there was an error and do something about it if necessary
      CALL WrScr( ErrMsg )
   END IF

   CALL SD_DestroyContState( dxdt, ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) THEN   ! Check if there was an error and do something about it if necessary
      CALL WrScr( ErrMsg )
   END IF

   !...............................................................................................................................
   ! Jacobian routines called in tight coupling
   !...............................................................................................................................

   CALL SD_JacobianPInput( Time, u(1), p, x, xd, z, OtherState, dYdu, dXdu, dXddu, dZdu, ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
      CALL WrScr( ErrMsg )
   END IF

   CALL SD_JacobianPContState( Time, u(1), p, x, xd, z, OtherState, dYdx, dXdx, dXddx, dZdx, ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
      CALL WrScr( ErrMsg )
   END IF

   CALL SD_JacobianPDiscState( Time, u(1), p, x, xd, z, OtherState, dYdxd, dXdxd, dXddxd, dZdxd, ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
      CALL WrScr( ErrMsg )
   END IF

   CALL SD_JacobianPConstrState( Time, u(1), p, x, xd, z, OtherState, dYdz, dXdz, dXddz, dZdz, ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
      CALL WrScr( ErrMsg )
   END IF


   !...............................................................................................................................
   ! Routines to pack data (to restart later)
   !...............................................................................................................................  
   CALL SD_Pack(Re_SaveAry, Db_SaveAry, Int_SaveAry, u(1), p, x, xd, z, OtherState, y, ErrStat, ErrMsg) 
     
   IF ( ErrStat /= ErrID_None ) THEN
      CALL WrScr( ErrMsg )
   END IF


   !...............................................................................................................................
   ! Routine to terminate program execution
   !...............................................................................................................................
   !CALL SD_End( u(1), p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
   !
   !IF ( ErrStat /= ErrID_None ) THEN
   !   CALL WrScr( ErrMsg )
   !END IF


   !...............................................................................................................................
   ! Routines to retreive packed data (unpack for restart)
   !...............................................................................................................................
   ! TODO:  BUG with Unpack and the added meshes?  GJH 6/12/13
 !  CALL SD_Unpack( Re_SaveAry, Db_SaveAry, Int_SaveAry, u(1), p, x, xd, z, OtherState, y, ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      CALL WrScr( ErrMsg )
   END IF


   !...............................................................................................................................
   ! Routines to copy data (not already tested)
   !...............................................................................................................................



   !...............................................................................................................................
   ! Routines to destroy data (not already tested)
   !...............................................................................................................................

   IF ( ALLOCATED( Re_SaveAry  ) ) DEALLOCATE( Re_SaveAry )
   IF ( ALLOCATED( Db_SaveAry  ) ) DEALLOCATE( Db_SaveAry )
   IF ( ALLOCATED( Int_SaveAry ) ) DEALLOCATE( Int_SaveAry )

!   CALL SD_DestroyPartialOutputPInput ( )  ! Jacobian Routine not yet implemented


   !...............................................................................................................................
   ! Routine to terminate program execution (again)
   !...............................................................................................................................

   CALL SD_End( u(1), p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) THEN
      CALL WrScr( ErrMsg )
   END IF

   !!............................................................................................................................
   !   ! Set exit error code if there was an error;
   !   !............................................................................................................................
   !   IF (Error) CALL ProgAbort( ' Simulation error level: '//TRIM(GetErrStr(ErrLev)), &   !This assumes PRESENT(ErrID) is .TRUE.
   !                               TrapErrors=.FALSE., TimeWait=3._ReKi )  ! wait 3 seconds (in case they double-clicked and got an error)
   !   
   !   !............................................................................................................................
   !   !  Write simulation times and stop
   !   !............................................................................................................................

      CALL RunTimes( StrtTime, UsrTime1, Time )
   
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
   TYPE(SD_Drvr_InitInput),       INTENT(   OUT )   :: InitInp
   INTEGER,                       INTENT(   OUT )   :: ErrStat              ! returns a non-zero value when an error occurs  
   CHARACTER(*),                  INTENT(   OUT )   :: ErrMsg               ! Error message if ErrStat /= ErrID_None
   
      ! Local variables  
         
   INTEGER                                          :: I                    ! generic integer for counting
   INTEGER                                          :: J                    ! generic integer for counting
   CHARACTER(   2)                                  :: strI                 ! string version of the loop counter

   INTEGER                                          :: UnIn                 ! Unit number for the input file
   INTEGER                                          :: UnEchoLocal          ! The local unit number for this module's echo file
   CHARACTER(1024)                                  :: EchoFile             ! Name of SubDyn echo file  
   CHARACTER(1024)                                  :: Line                 ! String to temporarially hold value of read line   
   CHARACTER(1024)                                  :: TmpPath              ! Temporary storage for relative path name
   CHARACTER(1024)                                  :: TmpFmt               ! Temporary storage for format statement
   CHARACTER(1024)                                  :: FileName             ! Name of SubDyn input file  
   
   UnEChoLocal=-1
   
   FileName = TRIM(inputFile)
   
   CALL GetNewUnit( UnIn )   
   CALL OpenFInpFile( UnIn, FileName, ErrStat )
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to open SubDyn Driver input file: '//FileName
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF

   
   CALL WrScr( 'Opening SubDyn Driver input file:  '//FileName )
   
   
   !-------------------------------------------------------------------------------------------------
   ! File header
   !-------------------------------------------------------------------------------------------------
   
   CALL ReadCom( UnIn, FileName, 'SubDyn Driver input file header line 1', ErrStat )
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read SubDyn Driver input file header line 1.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF


   CALL ReadCom( UnIn, FileName, 'SubDyn Driver input file header line 2', ErrStat )
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read SubDyn Driver input file header line 2.'
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
      
      CALL ReadCom( UnIn, FileName, 'SubDyn Driver input file header line 1', ErrStat, ErrMsg, UnEchoLocal )
   
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Failed to read SubDyn Driver input file header line 1.'
         ErrStat = ErrID_Fatal
         CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      END IF


      CALL ReadCom( UnIn, FileName, 'SubDyn Driver input file header line 2', ErrStat, ErrMsg, UnEchoLocal )
   
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Failed to read SubDyn Driver input file header line 2.'
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
   ! SubDyn section
   !-------------------------------------------------------------------------------------------------

      ! Header
      
   CALL ReadCom( UnIn, FileName, 'SubDyn header', ErrStat, ErrMsg, UnEchoLocal )
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Comment line.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
   
   
      ! HDInputFile
      
   CALL ReadVar ( UnIn, FileName, InitInp%SDInputFile, 'HDInputFile', &
                                    'SubDyn input filename', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read SDInputFile parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF 
   
   
      ! OutRootName
   
   CALL ReadVar ( UnIn, FileName, InitInp%OutRootName, 'OutRootName', &
                                    'SubDyn output root filename', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read OutRootName parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF   
     
   
      ! NSteps
   
   CALL ReadVar ( UnIn, FileName, InitInp%NSteps, 'NSteps', &
                                    'Number of time steps in the SubDyn simulation', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read NSteps parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF   
 
   
      ! TimeInterval   
   
   CALL ReadVar ( UnIn, FileName, InitInp%TimeInterval, 'TimeInterval', &
                                    'Time interval for any SubDyn inputs', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read TimeInterval parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF   
   
      ! TP_RefPoint   
   
   CALL ReadAry ( UnIn, FileName, InitInp%TP_RefPoint, 3, 'TP reference point', &
                                    'TP reference point', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read TP_RefPoint parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF   
   
      ! SubRotateZ   
   
   CALL ReadVar ( UnIn, FileName, InitInp%SubRotateZ, 'SubRotateZ', &
                                    'Rotation angle in degrees about Z axis.', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read SubRotateZ parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF   
   !-------------------------------------------------------------------------------------------------
   !  INPUTS section
   !-------------------------------------------------------------------------------------------------

      ! Header
      
   CALL ReadCom( UnIn, FileName, 'INPUTS header', ErrStat, ErrMsg, UnEchoLocal )
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Comment line.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
 
   
   
      ! InputsMod      
       
   CALL ReadVar ( UnIn, FileName, InitInp%InputsMod, 'InputsMod', &
                                    'Model for the inputs', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read InputsMod parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF   
   
   
      ! InputsFile      
       
   CALL ReadVar ( UnIn, FileName, InitInp%InputsFile, 'InputsFile', &
                                    'Filename for the SubDyn inputs', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read InputsFile parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF   
   
   
   !-------------------------------------------------------------------------------------------------
   ! STEADY STATE INPUTS section
   !-------------------------------------------------------------------------------------------------

      ! Header
      
   CALL ReadCom( UnIn, FileName, 'STEADY STATE INPUTS header', ErrStat, ErrMsg, UnEchoLocal )
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Comment line.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
   
   IF ( InitInp%InputsMod == 1 ) THEN
   
         ! uTPInSteady
         
      CALL ReadAry ( UnIn, FileName, InitInp%uTPInSteady, 6, 'uInSteady', &
                           'Steady-state TP displacements and rotations.', ErrStat,  ErrMsg, UnEchoLocal)         
       
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Failed to read uTPInSteady parameter.'
         ErrStat = ErrID_Fatal
         CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      END IF
   
   
         ! uDotTPInSteady
         
      CALL ReadAry ( UnIn, FileName, InitInp%uDotTPInSteady, 6, 'uDotTPInSteady', &
                           ' Steady-state TP translational and rotational velocities.', ErrStat,  ErrMsg, UnEchoLocal)         
       
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Failed to read uDotTPInSteady parameter.'
         ErrStat = ErrID_Fatal
         CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      END IF
      
      
         ! uDotDotTPInSteady
         
      CALL ReadAry ( UnIn, FileName, InitInp%uDotDotTPInSteady, 6, 'uDotDotTPInSteady', &
                           ' Steady-state TP translational and rotational accelerations.', ErrStat,  ErrMsg, UnEchoLocal)         
       
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Failed to read uDotDotTPInSteady parameter.'
         ErrStat = ErrID_Fatal
         CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      END IF
   ELSE
      InitInp%uTPInSteady    = 0.0
      InitInp%uDotTPInSteady = 0.0
      InitInp%uDotDotTPInSteady = 0.0
   END IF
   
   
   CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
   CLOSE( UnIn )
   
END SUBROUTINE ReadDriverInputFile

subroutine print_help()
    print '(a)', 'usage: cmdline [OPTIONS]'
    print '(a)', ''
    print '(a)', 'Without further options, cmdline prints the date and exits.'
    print '(a)', ''
    print '(a)', 'cmdline options:'
    print '(a)', ''
    print '(a)', '  -v, --version     print version information and exit'
    print '(a)', '  -h, --help        print usage information and exit'
    print '(a)', '  -t, --time        print time'
end subroutine print_help
 


SUBROUTINE RunTimes( StrtTime, UsrTime1, ZTime )
! This routine displays a message that gives that status of the simulation and the predicted end time of day.
!..................................................................................................................................

   IMPLICIT                        NONE

      ! Passed variables

   INTEGER                      :: StrtTime (8)                                    ! Start time of simulation
   REAL                         :: UsrTime1                                        ! User CPU time for simulation initialization.
   REAL(DbKi)                   :: ZTime                                           ! The final simulation time (not necessarially TMax)

      ! Local variables

   REAL                         :: ClckTime                                        ! Elapsed clock time for the simulation phase of the run.
   REAL                         :: Factor                                          ! Ratio of seconds to a specified time period.
   REAL                         :: TRatio                                          ! Ration of simulation time to elapsed clock time.
   REAL(ReKi), PARAMETER        :: SecPerDay = 24*60*60.0_ReKi                     ! Number of seconds per day

   REAL                         :: UsrTime                                         ! User CPU time for entire run.
   INTEGER                      :: EndTimes (8)                                    ! An array holding the ending clock time of the simulation.

   CHARACTER( 8)                :: TimePer


      ! Get the end times to compare with start times.

   CALL DATE_AND_TIME ( VALUES=EndTimes )
   CALL CPU_TIME ( UsrTime )


   ! Calculate the elapsed wall-clock time in seconds.

!bjj: I think this calculation will be wrong at certain times (e.g. if it's near midnight on the last day of the month), but to my knowledge, no one has complained...
   ClckTime =  0.001*( EndTimes(8) - StrtTime(8) ) + ( EndTimes(7) - StrtTime(7) ) + 60.0*( EndTimes(6) - StrtTime(6) ) &
            + 3600.0*( EndTimes(5) - StrtTime(5) ) + SecPerDay*( EndTimes(3) - StrtTime(3) )


      ! Calculate CPU times.

   UsrTime  = UsrTime - UsrTime1


   IF ( .NOT. EqualRealNos( UsrTime, 0.0 ) )  THEN

      TRatio = ZTime / UsrTime

      IF     ( UsrTime > SecPerDay )  THEN
         Factor = 1.0/SecPerDay
         TimePer = ' days'
      ELSEIF ( UsrTime >  3600.0 )  THEN
         Factor = 1.0/3600.0
         TimePer = ' hours'
      ELSEIF ( UsrTime >    60.0 )  THEN
         Factor = 1.0/60.0
         TimePer = ' minutes'
      ELSE
         Factor = 1.0
         TimePer = ' seconds'
      ENDIF

      CALL WrOver( '                                                                                   ' )
      CALL WrScr1( ' Total Real Time:       '//TRIM( Num2LStr( Factor*ClckTime      ) )//TRIM( TimePer ) )
      CALL WrScr ( ' Total CPU Time:        '//TRIM( Num2LStr( Factor*UsrTime       ) )//TRIM( TimePer ) )
      CALL WrScr ( ' Simulated Time:        '//TRIM( Num2LStr( Factor*REAL( ZTime ) ) )//TRIM( TimePer ) )
      CALL WrScr ( ' Time Ratio (Sim/CPU):  '//TRIM( Num2LStr( TRatio ) ) )

   ENDIF

   RETURN
END SUBROUTINE RunTimes
   
   !----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SimStatus( PrevSimTime, PrevClockTime, ZTime, TMax )
! This routine displays a message that gives that status of the simulation and the predicted end time of day.
!..................................................................................................................................

   IMPLICIT                        NONE

      ! Passed variables
   REAL(DbKi), INTENT(IN)       :: ZTime                                           ! Current simulation time (s)
   REAL(DbKi), INTENT(IN)       :: TMax                                            ! Expected simulation time (s)
   REAL(DbKi), INTENT(INOUT)    :: PrevSimTime                                     ! Previous time message was written to screen (s > 0)
   REAL(ReKi), INTENT(INOUT)    :: PrevClockTime                                   ! Previous clock time in seconds past midnight


      ! Local variables.

   REAL(ReKi)                   :: CurrClockTime                                   ! Current time in seconds past midnight.
   REAL(ReKi)                   :: DeltTime                                        ! The amount of time elapsed since the last call.
   REAL(ReKi)                   :: EndTime                                         ! Approximate time of day when simulation will complete.
   REAL(ReKi), PARAMETER        :: InSecHr  = 1.0_ReKi/3600.0_ReKi                 ! Inverse of the number of seconds in an hour
   REAL(ReKi), PARAMETER        :: InSecMn  = 1.0_ReKi/  60.0_ReKi                 ! Inverse of the number of seconds in a minute
   REAL(ReKi)                   :: SimTimeLeft                                     ! Approximate clock time remaining before simulation completes

   REAL(ReKi), PARAMETER        :: SecPerDay = 24*60*60.0_ReKi                     ! Number of seconds per day

   INTEGER(4)                   :: EndHour                                         ! The hour when the simulations is expected to complete.
   INTEGER(4)                   :: EndMin                                          ! The minute when the simulations is expected to complete.
   INTEGER(4)                   :: EndSec                                          ! The second when the simulations is expected to complete.
   INTEGER(4)                   :: TimeAry  (8)                                    ! An array containing the elements of the start time.

   CHARACTER( 8)                :: ETimeStr                                        ! String containing the end time.


   IF ( ZTime <= PrevSimTime ) RETURN


      ! How many seconds past midnight?

   CALL DATE_AND_TIME ( Values=TimeAry )
   CurrClockTime = TimeValues2Seconds( TimeAry )

      ! Calculate elapsed clock time

   DeltTime = CurrClockTime - PrevClockTime


      ! We may have passed midnight since the last revoultion.  We will assume that (ZTime - PrevSimTime) of simulation time doesn't take more than a day.

   IF ( CurrClockTime < PrevClockTime )  THEN
      DeltTime = DeltTime + SecPerDay
   ENDIF


      ! Estimate the end time in hours, minutes, and seconds

   SimTimeLeft = ( TMax - ZTime )*DeltTime/( ZTime - PrevSimTime )            ! DeltTime/( ZTime - PrevSimTime ) is the delta_ClockTime divided by the delta_SimulationTime
   EndTime  =  MOD( CurrClockTime+SimTimeLeft, SecPerDay )
   EndHour  =  INT(   EndTime*InSecHr )
   EndMin   =  INT( ( EndTime - REAL( 3600*EndHour ) )*InSecMn )
   EndSec   = NINT(   EndTime - REAL( 3600*EndHour + 60*EndMin ) )

   WRITE (ETimeStr,"(I2.2,2(':',I2.2))")  EndHour, EndMin, EndSec

   CALL WrOver ( ' Timestep: '//TRIM( Num2LStr( NINT( ZTime ) ) )//' of '//TRIM( Num2LStr( TMax ) )// &
                 ' seconds.  Estimated final completion at '//ETimeStr//'.'                             )


      ! Let's save this time as the previous time for the next call to the routine
   PrevClockTime = CurrClockTime
   PrevSimTime   = ZTime

   RETURN
END SUBROUTINE SimStatus
!----------------------------------------------------------------------------------------------------------------------------------
FUNCTION TimeValues2Seconds( TimeAry )
! This routine takes an array of time values such as that returned from
!     CALL DATE_AND_TIME ( Values=TimeAry )
! and converts TimeAry to the number of seconds past midnight.
!..................................................................................................................................

      ! Passed variables:
   INTEGER, INTENT(IN)          :: TimeAry  (8)                                    ! An array containing the elements of the time
   REAL(ReKi)                   :: TimeValues2Seconds                              ! Current time in seconds past midnight


   TimeValues2Seconds = 3600*TimeAry(5) + 60*TimeAry(6) + TimeAry(7) + 0.001_ReKi*TimeAry(8)

END FUNCTION TimeValues2Seconds
!----------------------------------------------------------------------------------------------------------------------------------

END PROGRAM TestSubDyn
