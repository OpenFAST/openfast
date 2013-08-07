!**********************************************************************************************************************************
! SubDyn_DriverCode: This code tests the SubDyn modules
!..................................................................................................................................
! LICENSING
! Copyright (C) 2013  National Renewable Energy Laboratory
!
!    This file is part of SubDyn.
!
!    SubDyn is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
!    published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License along with SubDyn.
!    If not, see <http://www.gnu.org/licenses/>.
!
!**********************************************************************************************************************************
PROGRAM TestSubDyn

   USE NWTC_Library
   USE SubDyn
   USE TempMod
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
   
   !...............................................................................................................................
   ! Routines called in initialization
   !...............................................................................................................................

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
   
   
  
   
         ! Initialize the module
   
   CALL SubDyn_Init( InitInData, u(1), p,  x, xd, z, OtherState, y, TimeInterval, InitOutData, ErrStat1, ErrMsg1 )
   


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
      
      CALL SubDyn_CalcOutput( Time, u(1), p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
         CALL WrScr( ErrMsg )
      END IF

         
         ! Get state variables at next step: INPUT at step n, OUTPUT at step n + 1
                                     
      CALL SubDyn_UpdateStates( Time, n, u, InputTime, p, x, xd, z, z_next, OtherState, ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
         CALL WrScr( ErrMsg )
      END IF     
      
      
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
!      CALL SubDyn_CalcOutput( Time, u(1), p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
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
   !CALL SubDyn_End( u(1), p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
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

   CALL SubDyn_End( u(1), p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) THEN
      CALL WrScr( ErrMsg )
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
   
END PROGRAM TestSubDyn
