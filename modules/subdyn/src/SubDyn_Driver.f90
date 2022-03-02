!**********************************************************************************************************************************
! SubDyn_DriverCode: This code tests the SubDyn modules
!..................................................................................................................................
! LICENSING
! Copyright (C) 2013-2016  National Renewable Energy Laboratory
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
PROGRAM SubDyn_Driver

   USE NWTC_Library
   USE SubDyn
   USE SubDyn_Types
   USE SubDyn_Output
   USE FEM, only: FINDLOCI
   USE VersionInfo

   IMPLICIT NONE

   INTEGER(IntKi), PARAMETER                          :: NumInp = 1           ! Number of inputs sent to SD_UpdateStates
   
   type ALoadType
      integer(IntKi)         :: NodeID            ! Joint and then Node ID where force is applied
      real(ReKi)             :: SteadyLoad(6)     ! Steady Load (Fx, Fy, Fz, Mx, My, Mz)
      real(ReKi),allocatable :: UnsteadyLoad(:,:) ! Unsteady Load nx7 : (Time, Fx, Fy, Fz, Mx, My, Mz)
      integer(IntKi)         :: iTS               ! Optimization index to find closest time stamp in user defined time series of unsteady load
   end type ALoadType
   
   TYPE SD_dvr_InitInput
      LOGICAL         :: Echo
      REAL(ReKi)      :: Gravity
      CHARACTER(1024) :: SDInputFile
      REAL(ReKi)      :: WtrDpth
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
      type(ALoadType), allocatable :: AppliedLoads(:)  ! 7 x nSteadyForces: JointID, Fx, Fy, Fz, Mx, My, Mz
   END TYPE SD_dvr_InitInput
   
   
      ! Program variables

   REAL(DbKi)                      :: Time                 ! Variable for storing time, in seconds
   REAL(DbKi)                      :: TimeInterval         ! Interval between time steps, in seconds
   REAL(DbKi)                      :: InputTime(NumInp)    ! Variable for storing time associated with inputs, in seconds
   
   TYPE(SD_InitInputType)          :: InitInData           ! Input data for initialization
   TYPE(SD_InitOutputType)         :: InitOutData          ! Output data from initialization

   TYPE(SD_ContinuousStateType)    :: x                    ! Continuous states
   TYPE(SD_DiscreteStateType)      :: xd                   ! Discrete states
   TYPE(SD_ConstraintStateType)    :: z                    ! Constraint states
   TYPE(SD_OtherStateType)         :: OtherState           ! Other states
   TYPE(SD_MiscVarType)            :: m                    ! Misc/optimization variables

   TYPE(SD_ParameterType)          :: p                    ! Parameters
   TYPE(SD_InputType)              :: u(NumInp)            ! System inputs
   TYPE(SD_OutputType)             :: y                    ! System outputs
   TYPE(ALoadType), pointer        :: AL                   ! Applied Load (alias to shorten notations)


   INTEGER(IntKi)                  :: n                    ! Loop counter (for time step)
   INTEGER(IntKi)                  :: ErrStat, ErrStat1, ErrStat2, ErrStat3          ! Status of error message
   CHARACTER(1024)                 :: ErrMsg, ErrMsg1, ErrMsg2, ErrMsg3              ! Error message if ErrStat /= ErrID_None


   CHARACTER(1024)                 :: dvrFilename         ! Filename and path for the driver input file.  This is passed in as a command line argument when running the Driver exe.
   TYPE(SD_dvr_InitInput), target  :: drvrInitInp          ! Initialization data for the driver program
   INTEGER                         :: UnIn                 ! Unit number for the input file
   INTEGER                         :: UnEcho          ! The local unit number for this module's echo file
   INTEGER(IntKi)                  :: UnSD_Out             ! Output file identifier
   REAL(ReKi), ALLOCATABLE         :: SDin(:,:)            ! Variable for storing time, forces, and body velocities, in m/s or rad/s for SubDyn inputs
   INTEGER(IntKi)                  :: I,J                  ! Generic loop counter
   INTEGER(IntKi)                  :: iLoad                ! Index on loads             
   INTEGER(IntKi)                  :: iNode                ! Index on nodes
   INTEGER(IntKi)                  :: JointID              ! JointID
   REAL(ReKi)                      :: dcm (3,3)            ! The resulting transformation matrix from X to x, (-).
   CHARACTER(10)                   :: AngleMsg             ! For debugging, a string version of the largest rotation input
   real(ReKi)                      :: UnsteadyLoad(6)      ! Unsteady Load interpolated at a given time
   
      ! Other/Misc variables
   REAL(DbKi)                      :: TiLstPrn             ! The time of the last print
   REAL(DbKi)                      :: TMax
   REAL(DbKi)                      :: OutTime              ! Used to determine if output should be generated at this simulation time
   REAL(ReKi)                      :: PrevClockTime        ! Clock time at start of simulation in seconds
   REAL(ReKi)                      :: UsrTime1             ! User CPU time for simulation initialization
   INTEGER                         :: StrtTime (8)         ! Start time of simulation
   CHARACTER(200)                  :: git_commit           ! String containing the current git commit hash
   TYPE(ProgDesc), PARAMETER       :: version   = ProgDesc( 'SubDyn Driver', '', '' )  ! The version number of this program.
   !...............................................................................................................................
   ! Routines called in initialization
   !...............................................................................................................................
   ErrMsg  = ""
   ErrStat = ErrID_None
   UnEcho=-1
   UnIn  =-1
   
   ! Get the current time
   CALL DATE_AND_TIME ( Values=StrtTime )                               ! Let's time the whole simulation
   CALL CPU_TIME ( UsrTime1 )                                           ! Initial time (this zeros the start time when used as a MATLAB function)
   PrevClockTime = TimeValues2Seconds( StrtTime )                       ! We'll use this time for the SimStats routine
   TiLstPrn      = 0.0_DbKi                                             ! The first value of ZTime, used to write simulation stats to screen (s)
  
   ! Initialize the NWTC Subroutine Library
   CALL NWTC_Init( )
   
   ! Display the copyright notice
   CALL DispCopyrightLicense( version%Name )
   ! Obtain OpenFAST git commit hash
   git_commit = QueryGitVersion()
   ! Tell our users what they're running
   CALL WrScr( ' Running '//TRIM( version%Name )//' a part of OpenFAST - '//TRIM(git_Commit)//NewLine//' linked with '//TRIM( NWTC_Ver%Name )//NewLine )
   
   ! Set the abort error level to a fatal error
   AbortErrLev = ErrID_Fatal
   
   IF ( command_argument_count() /= 1 )  then
      CALL print_help()
      STOP
   endif

   ! Parse the driver input file and run the simulation based on that file
   IF ( command_argument_count() == 1 ) THEN
      CALL get_command_argument(1, dvrFilename)

      CALL ReadDriverInputFile( dvrFilename, drvrInitInp);
      InitInData%g            = drvrInitInp%Gravity
      InitInData%SDInputFile  = drvrInitInp%SDInputFile
      InitInData%RootName     = drvrInitInp%OutRootName
      InitInData%TP_RefPoint  = drvrInitInp%TP_RefPoint
      InitInData%SubRotateZ   = drvrInitInp%SubRotateZ
      TimeInterval            = drvrInitInp%TimeInterval
      InitInData%WtrDpth      = drvrInitInp%WtrDpth
   END IF

   TMax = TimeInterval * drvrInitInp%NSteps
   
   ! Initialize SubDyn module
   CALL SD_Init( InitInData, u(1), p,  x, xd, z, OtherState, y, m, TimeInterval, InitOutData, ErrStat2, ErrMsg2 ); call AbortIfFailed()

   ! Sanity check for outputs
   if (p%NumOuts==0) then
      call WrScr('Warning: No output channels were selected in SubDyn. No output file will be created!')
   endif
   if (p%OutSwtch==2) then
      p%OutSwtch=1
      call WrScr('Warning: Overring `OutSwitch` to 1 to generate outputs with the driver.')
      ! TODO not pretty, it'd be nicer to tell SubDyn it's running with the driver
      drvrInitInp%OutRootName = TRIM(drvrInitInp%OutRootName)//'.SD'
      CALL SDOUT_OpenOutput( SD_ProgDesc, drvrInitInp%OutRootName, p, InitOutData, ErrStat2, ErrMsg2 ); 
   endif


   ! Read Input time series data from a file
   CALL AllocAry(SDin, drvrInitInp%NSteps, 19, 'SDinput array', ErrStat2, ErrMsg2); call AbortIfFailed()
   SDin(:,:)=0.0_ReKi
   IF ( drvrInitInp%InputsMod == 2 ) THEN
      ! Open the  inputs data file
      CALL GetNewUnit( UnIn ) 
      CALL OpenFInpFile ( UnIn, drvrInitInp%InputsFile, ErrStat2, ErrMsg2); Call AbortIfFailed()
      DO n = 1,drvrInitInp%NSteps
         ! TODO Add safety for backward compatibility if only 13 columns
         READ (UnIn,*,IOSTAT=ErrStat2) (SDin (n,J), J=1,19)
         ErrMsg2 = ' Error reading line '//trim(Num2LStr(n))//' of file: '//trim(drvrInitInp%InputsFile)
         call AbortIfFailed()
      END DO  
      CLOSE ( UnIn ) 
   else
      ! We fill an array with constant values
      do n = 0,drvrInitInp%NSteps-1 ! Loop on time steps, starts at 0
         SDin(n+1,1) = n*TimeInterval
         SDin(n+1,2:7 ) = drvrInitInp%uTPInSteady(1:6)     ! Displacements
         SDin(n+1,8:13) = drvrInitInp%uDotTPInSteady(1:6)  ! Velocities
         !SDin(n+1,14:19) = drvrInitInp%uDotDotTPInSteady(1:6)  ! Accelerations
      enddo
   end if 

   ! Setup Applied Loads 
   do iLoad=1,size(drvrInitInp%AppliedLoads)
      AL => drvrInitInp%AppliedLoads(iLoad)
      JointID = AL%NodeID
      AL%NodeID = findloci(p%NodeID2JointID, JointID) ! Replace JointID with Index
      if (AL%NodeID<0) then
         ErrStat2 = ErrID_Fatal
         ErrMsg2  = 'Applied load JointID '//trim(num2lstr(JointID))//' is not found in SubDyn joint list.'
         call AbortIfFailed()
      endif
      AL%iTS=1 ! important init
   enddo

  
   ! Destroy initialization data
   CALL SD_DestroyInitInput(  InitInData,  ErrStat2, ErrMsg2 ); call AbortIfFailed()
   CALL SD_DestroyInitOutput( InitOutData, ErrStat2, ErrMsg2 ); call AbortIfFailed()

   !...............................................................................................................................
   ! Routines called in loose coupling -- the glue code may implement this in various ways
   !...............................................................................................................................
   ! Force the displacement of the interface node in the global Z direction to be the sag of the column under it's own weight
   ! u(1)%UFL(3) =-12.958  !this is for testbeam3

   ! TEMPORARY HACK FOR CONTROLLABLE CABLES
   !allocate(u(1)%CableDeltaL(5))
   !!u(1)%CableDeltaL= 1.0e7_ReKi
   !u(1)%CableDeltaL= 0.0e7_ReKi

   call WrScr('')
   DO n = 0,drvrInitInp%NSteps-1 ! Loop on time steps, starts at 0

      Time = n*TimeInterval
      InputTime(1) = Time

      ! Set module inputs u (likely from the outputs of another module or a set of test conditions) here:
      IF ( u(1)%TPMesh%Initialized ) THEN 
         ! Input displacements, velocities and potentially accelerations
         u(1)%TPMesh%TranslationDisp(:,1)   = SDin(n+1,2:4) 
         CALL SmllRotTrans( 'InputRotation', REAL(SDin(n+1,5),reki), REAL(SDin(n+1,6),reki), REAL(SDin(n+1,7),reki), dcm, 'Junk', ErrStat, ErrMsg )            
         u(1)%TPMesh%Orientation(:,:,1)     = dcm 
         u(1)%TPMesh%TranslationVel(:,1)    = SDin(n+1,8:10)  
         u(1)%TPMesh%RotationVel(:,1)       = SDin(n+1,11:13) 

         IF ( drvrInitInp%InputsMod == 2 ) THEN
            u(1)%TPMesh%TranslationAcc(:,1)    = SDin(n+1,14:16) 
            u(1)%TPMesh%RotationAcc(:,1)       = SDin(n+1,17:19)
         ELSE ! constant inputs
            u(1)%TPMesh%TranslationAcc(:,1)    = drvrInitInp%uDotDotTPInSteady(1:3)  
            u(1)%TPMesh%RotationAcc(:,1)       = drvrInitInp%uDotDotTPInSteady(4:6) 
         END IF
      END IF   
      ! Set LMesh applied loads
      if ( u(1)%LMesh%Initialized ) then 
         ! Default, set all external load to 0.0
         u(1)%LMesh%Force  (:,:) = 0.0
         u(1)%LMesh%Moment (:,:) = 0.0
         ! Add applied loads
         do iLoad=1, size(drvrInitInp%AppliedLoads)
            AL => drvrInitInp%AppliedLoads(iLoad)
            iNode = AL%NodeID
            u(1)%LMesh%Force(:,iNode)  = u(1)%LMesh%Force(:,iNode)  + AL%SteadyLoad(1:3)
            u(1)%LMesh%Moment(:,iNode) = u(1)%LMesh%Moment(:,iNode) + AL%SteadyLoad(4:6)
            if (allocated(AL%UnsteadyLoad)) then
               call interpTimeValue(AL%UnsteadyLoad, Time, AL%iTS, UnsteadyLoad)
               u(1)%LMesh%Force(:,iNode)  = u(1)%LMesh%Force(:,iNode)  + UnsteadyLoad(1:3)
               u(1)%LMesh%Moment(:,iNode) = u(1)%LMesh%Moment(:,iNode) + UnsteadyLoad(4:6)
            endif
         enddo
      endif


      ! Calculate outputs at n
      CALL SD_CalcOutput( Time, u(1), p, x, xd, z, OtherState, y, m, ErrStat2, ErrMsg2); call AbortIfFailed()
      ! Get state variables at next step: INPUT at step n, OUTPUT at step n + 1
      CALL SD_UpdateStates( Time, n, u, InputTime, p, x, xd, z, OtherState, m, ErrStat2, ErrMsg2); call AbortIfFailed()
      ! Display simulation status every SttsTime-seconds:
      IF ( Time - TiLstPrn >= 1 )  THEN
         CALL SimStatus( TiLstPrn, PrevClockTime, Time, TMax )
      ENDIF   

   END DO ! Loop on n, time steps

   ! Routine to terminate program execution
   CALL SD_End( u(1), p, x, xd, z, OtherState, y, m, ErrStat2, ErrMsg2)
   IF ( ErrStat /= ErrID_None ) THEN
      CALL WrScr( ErrMsg )
   END IF

   ! Write simulation times and stop
   CALL RunTimes( StrtTime, UsrTime1, StrtTime, UsrTime1, Time )
   
CONTAINS
   SUBROUTINE AbortIfFailed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'SubDyn_Driver') 
        IF ( ErrStat /= ErrID_None ) THEN
           CALL WrScr( ErrMsg )
           CALL WrScr('')
        END IF
        if (ErrStat >= AbortErrLev) then
           call CleanUp()
           STOP
        endif
   END SUBROUTINE AbortIfFailed

   SUBROUTINE CleanUp()
      integer :: iForce
      if(allocated(drvrInitInp%AppliedLoads)) then
          do iForce=1,size(drvrInitInp%AppliedLoads) 
             if (allocated(drvrInitInp%AppliedLoads(iForce)%UnsteadyLoad)) then
                deallocate(drvrInitInp%AppliedLoads(iForce)%UnsteadyLoad)
             endif
          enddo
       endif
      if(UnEcho>0) CLOSE(UnEcho)
      if(UnEcho>0) CLOSE( UnIn)
      if(allocated(SDin)) deallocate(SDin)
   END SUBROUTINE CleanUp

   !-------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE ReadDriverInputFile( inputFile, InitInp)
      CHARACTER(*),                 INTENT( IN    )   :: inputFile
      TYPE(SD_dvr_InitInput),       INTENT(   OUT )   :: InitInp
      ! Local variables  
      INTEGER                                          :: I                    ! generic integer for counting
      INTEGER                                          :: J                    ! generic integer for counting
      INTEGER                                          :: iDummy               ! dummy integer
      CHARACTER(   2)                                  :: strI                 ! string version of the loop counter

      CHARACTER(1024)                                  :: EchoFile             ! Name of SubDyn echo file  
      CHARACTER(1024)                                  :: Line                 ! String to temporarially hold value of read line   
      CHARACTER(1024)                                  :: TmpPath              ! Temporary storage for relative path name
      CHARACTER(1024)                                  :: TmpFmt               ! Temporary storage for format statement
      CHARACTER(1024)                                  :: FileName             ! Name of SubDyn input file  
      CHARACTER(1024)                                  :: PriPath             ! Path Name of SubDyn input file  
   
      UnEcho=-1
      UnIn  =-1
   
      FileName = TRIM(inputFile)
      ! Primary path, relative files will be based on it
      CALL GetPath( FileName, PriPath )
   
      CALL GetNewUnit( UnIn )   
      CALL OpenFInpFile( UnIn, FileName, ErrStat2, ErrMsg2);
      call AbortIfFailed()
   
      CALL WrScr( 'Opening SubDyn Driver input file:  '//FileName )
      
      ! Read until "echo"
      CALL ReadCom( UnIn, FileName, 'SubDyn Driver input file header line 1', ErrStat2, ErrMsg2); call AbortIfFailed()
      CALL ReadCom( UnIn, FileName, 'SubDyn Driver input file header line 2', ErrStat2, ErrMsg2); call AbortIfFailed()
      CALL ReadVar ( UnIn, FileName, InitInp%Echo, 'Echo', 'Echo Input', ErrStat2, ErrMsg2); call AbortIfFailed()
      ! If we echo, we rewind
      IF ( InitInp%Echo ) THEN
         EchoFile = TRIM(FileName)//'.echo'
         CALL GetNewUnit( UnEcho )   
         CALL OpenEcho ( UnEcho, EchoFile, ErrStat, ErrMsg ); call AbortIfFailed()
         REWIND(UnIn)
         CALL ReadCom( UnIn, FileName, 'SubDyn Driver input file header line 1', ErrStat2, ErrMsg2, UnEcho); call AbortIfFailed()
         CALL ReadCom( UnIn, FileName, 'SubDyn Driver input file header line 2', ErrStat2, ErrMsg2, UnEcho); call AbortIfFailed()
         CALL ReadVar ( UnIn, FileName, InitInp%Echo, 'Echo', 'Echo the input file data', ErrStat2, ErrMsg2, UnEcho); call AbortIfFailed()
      END IF
      !---------------------- ENVIRONMENTAL CONDITIONS -------------------------------------------------
      CALL ReadCom( UnIn, FileName, 'Environmental conditions header', ErrStat2, ErrMsg2, UnEcho); call AbortIfFailed()
      CALL ReadVar( UnIn, FileName, InitInp%Gravity, 'Gravity', 'Gravity', ErrStat2, ErrMsg2, UnEcho); call AbortIfFailed()
      CALL ReadVar( UnIn, FileName, InitInp%WtrDpth, 'WtrDpth', 'WtrDpth', ErrStat2, ErrMsg2, UnEcho); call AbortIfFailed()
      !---------------------- SubDyn -------------------------------------------------------------------
      CALL ReadCom( UnIn, FileName, 'SubDyn header', ErrStat2, ErrMsg2, UnEcho); call AbortIfFailed()
      CALL ReadVar( UnIn, FileName, InitInp%SDInputFile, 'HDInputFile', 'SubDyn input filename', ErrStat2, ErrMsg2, UnEcho); call AbortIfFailed()
      CALL ReadVar( UnIn, FileName, InitInp%OutRootName, 'OutRootName', 'SubDyn output root filename', ErrStat2, ErrMsg2, UnEcho); call AbortIfFailed()
      CALL ReadVar( UnIn, FileName, InitInp%NSteps     , 'NSteps', 'Number of time steps in the SubDyn simulation', ErrStat2, ErrMsg2, UnEcho); call AbortIfFailed()
      CALL ReadVar( UnIn, FileName, InitInp%TimeInterval, 'TimeInterval', 'Time interval for any SubDyn inputs', ErrStat2, ErrMsg2, UnEcho); call AbortIfFailed()
      CALL ReadAry( UnIn, FileName, InitInp%TP_RefPoint, 3, 'TP reference point', 'TP reference point', ErrStat2, ErrMsg2, UnEcho); call AbortIfFailed()
      CALL ReadVar( UnIn, FileName, InitInp%SubRotateZ, 'SubRotateZ', 'Rotation angle in degrees about Z axis.', ErrStat2, ErrMsg2, UnEcho); call AbortIfFailed()
      !---------------------- INPUTS -------------------------------------------------------------------
      CALL ReadCom( UnIn, FileName, 'INPUTS header', ErrStat2, ErrMsg2, UnEcho); call AbortIfFailed()
      CALL ReadVar( UnIn, FileName, InitInp%InputsMod , 'InputsMod', 'Model for the inputs', ErrStat2, ErrMsg2, UnEcho); call AbortIfFailed()
      CALL ReadVar( UnIn, FileName, InitInp%InputsFile, 'InputsFile', 'Filename for the SubDyn inputs', ErrStat2, ErrMsg2, UnEcho); call AbortIfFailed()
      !---------------------- STEADY INPUTS (for InputsMod = 1) ----------------------------------------
      CALL ReadCom( UnIn, FileName, 'STEADY STATE INPUTS header', ErrStat2, ErrMsg2, UnEcho); call AbortIfFailed()
      IF ( InitInp%InputsMod == 1 ) THEN
         CALL ReadAry ( UnIn, FileName, InitInp%uTPInSteady      , 6, 'uInSteady',         'Steady-state TP displacements and rotations.', ErrStat2,  ErrMsg2, UnEcho)         
         CALL ReadAry ( UnIn, FileName, InitInp%uDotTPInSteady   , 6, 'uDotTPInSteady',    'Steady-state TP translational and rotational velocities.', ErrStat2,  ErrMsg2, UnEcho)         
         CALL ReadAry ( UnIn, FileName, InitInp%uDotDotTPInSteady, 6, 'uDotDotTPInSteady', 'Steady-state TP translational and rotational accelerations.', ErrStat2,  ErrMsg2, UnEcho)         
      ELSE
         InitInp%uTPInSteady       = 0.0
         InitInp%uDotTPInSteady    = 0.0
         InitInp%uDotDotTPInSteady = 0.0
         CALL ReadCom( UnIn, FileName, '0.0   0.0   0.0   0.0   0.0   0.0   uTPInSteady     ', ErrStat2, ErrMsg2, UnEcho); call AbortIfFailed()
         CALL ReadCom( UnIn, FileName, '0.0   0.0   0.0   0.0   0.0   0.0   uDotTPInSteady  ', ErrStat2, ErrMsg2, UnEcho); call AbortIfFailed()
         CALL ReadCom( UnIn, FileName, '0.0   0.0   0.0   0.0   0.0   0.0   uDotTPInSteady  ', ErrStat2, ErrMsg2, UnEcho); call AbortIfFailed()
      END IF
      CALL AbortIfFailed()
      !---------------------- FORCES ----------------------------------------
      CALL ReadCom( UnIn, FileName, '--- FORCES INPUTS header', ErrStat2, ErrMsg2, UnEcho); call AbortIfFailed()
      CALL ReadVar ( UnIn, FileName, iDummy,  'nApplied Forces', 'Number of applied forces', ErrStat2,  ErrMsg2, UnEcho); 
      !call AbortIfFailed()
      if (ErrStat2/=ErrID_None) then
         ! TODO Temporary
         call LegacyWarning('Applied loads input missing.')
         allocate(InitInp%AppliedLoads(0), stat=ErrStat2); ErrMsg2='Allocating Forces'; call AbortIfFailed()
      else
         allocate(InitInp%AppliedLoads(iDummy), stat=ErrStat2); ErrMsg2='Allocating Forces'; call AbortIfFailed()
         CALL ReadCom( UnIn, FileName, 'JointID    Fx     Fy    Fz     Mx     My     Mz', ErrStat2, ErrMsg2, UnEcho); call AbortIfFailed()
         CALL ReadCom( UnIn, FileName, ' (-)       (N)   (N)    (N)   (Nm)   (Nm)   (Nm)', ErrStat2, ErrMsg2, UnEcho); call AbortIfFailed()
         do i=1,iDummy
            ! Read line and extract loads
            read(UnIn, fmt='(A)', iostat=ErrStat2) Line ; ErrMsg2='Erro reading force input line'//num2lstr(i); call AbortIfFailed()
            call readAppliedForce(Line, InitInp%AppliedLoads(i), PriPath, Errstat2, ErrMsg2); call AbortIfFailed()
         enddo
      endif

   
      if(UnEcho>0) CLOSE( UnEcho )
      if(UnIn>0)   CLOSE( UnIn   )
      ! --- Perform input checks and triggers
      ! If no root is provided, use the SDInputFile
      IF ( LEN_TRIM(InitInp%OutRootName) == 0 ) THEN
         CALL GetRoot(InitInp%SDInputFile, InitInp%OutRootName)
      END IF
      IF ( PathIsRelative( InitInp%SDInputFile ) ) then
         InitInp%SDInputFile = TRIM(PriPath)//TRIM(InitInp%SDInputFile)
      END IF
      IF ( PathIsRelative( InitInp%OutRootName ) ) then
         InitInp%OutRootName = TRIM(PriPath)//TRIM(InitInp%OutRootName)
      endif
      IF ( PathIsRelative( InitInp%InputsFile ) ) then
         InitInp%InputsFile = TRIM(PriPath)//TRIM(InitInp%InputsFile)
      endif

   END SUBROUTINE ReadDriverInputFile

   subroutine readAppliedForce(Line, AL, PriPath, Errstat, ErrMsg)
      character(*         ), intent(in   ) :: Line    ! Input line from input file
      type     (ALoadType) , intent(inout) :: AL      ! Applied force
      character(*         ), intent(in   ) :: PriPath ! Path to base relative file from
      integer  (IntKi     ), intent(out  ) :: ErrStat ! Error status of the operation
      character(*         ), intent(out  ) :: ErrMsg  ! Error message if ErrStat /    = ErrID_None
      character(255) :: StrArray(8) ! Array of strings extracted from line
      ! Set default err stat
      AL%SteadyLoad=-9999.9_ReKi ! Init
      if (allocated(AL%UnsteadyLoad)) deallocate(AL%UnsteadyLoad)
      ! Convert line to str array
      StrArray(:)='';
      CALL ReadCAryFromStr(Line, StrArray, 8, 'StrArray', 'StrArray', ErrStat, ErrMsg)! NOTE:No Error handling!
      ErrStat=ErrID_Fatal
      ErrMsg='Error reading force inputs.'//char(10)//'Prolematic line: '//trim(Line)
      ! NodeID
      if (.not. is_int(StrArray(1), AL%NodeID) ) then
         ErrMsg=trim(ErrMsg)//achar(13)//achar(10)//'NodeID needs to be an integer.'
         return
      endif
      ! Steady Load
      if (.not. is_numeric(StrArray(2), AL%SteadyLoad(1)) ) return
      if (.not. is_numeric(StrArray(3), AL%SteadyLoad(2)) ) return
      if (.not. is_numeric(StrArray(4), AL%SteadyLoad(3)) ) return
      if (.not. is_numeric(StrArray(5), AL%SteadyLoad(4)) ) return
      if (.not. is_numeric(StrArray(6), AL%SteadyLoad(5)) ) return
      if (.not. is_numeric(StrArray(7), AL%SteadyLoad(6)) ) return
      ! Unsteady Load
      if (len_trim(StrArray(8))>0) then
         ErrMsg='Unsteady load file not yet supported but the input file `'//trim(StrArray(8))//'` was given.'//char(10)//'Prolematic line: '//trim(Line)
         call ReadDelimFile(StrArray(8), 7, AL%UnsteadyLoad, ErrStat, ErrMsg, 1, priPath)
         return
         ! TODO read file here and allocate. See new AeroDyn driver to read csv

      endif
      print'(A,I5,6(E10.2))','    Applied Load: ',AL%NodeID, AL%SteadyLoad
      ! Success
      ErrStat=ErrID_None
      ErrMsg=''
   end subroutine readAppliedForce

   subroutine print_help()
       print '(a)', 'usage: '
       print '(a)', ''
       print '(a)', 'SubDynDriver.exe driverfilename'
       print '(a)', ''
       print '(a)', 'Where driverfilename is the name of the SubDyn driver input file.'
       print '(a)', ''
   end subroutine print_help

   subroutine LegacyWarning(Message)
      character(len=*), intent(in) :: Message
      call WrScr('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
      call WrScr('Warning: the SubDyn driver input file is not at the latest format!' )
      call WrScr('         Visit: https://openfast.readthedocs.io/en/dev/source/user/api_change.html')
      call WrScr('> Issue: '//trim(Message))
      call WrScr('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
   end subroutine LegacyWarning

   ! --------------------------------------------------------------------------------
   ! --- Generic routines (also present in other modules, e.g. OLAF, AD Driver) 
   ! --------------------------------------------------------------------------------
   function is_numeric(string, x)
      implicit none
      character(len=*), intent(in) :: string
      real(reki), intent(out) :: x
      logical :: is_numeric
      integer :: e,n
      character(len=12) :: fmt
      x = 0.0_reki
      n=len_trim(string)
      write(fmt,'("(F",I0,".0)")') n
      read(string,fmt,iostat=e) x
      is_numeric = e == 0
   end function is_numeric

   function is_int(string, x)
      implicit none
      character(len=*), intent(in) :: string
      integer(IntKi), intent(out) :: x
      logical :: is_int
      integer :: e,n
      character(len=12) :: fmt
      x = 0
      n=len_trim(string)
      write(fmt,'("(I",I0,")")') n
      read(string,fmt,iostat=e) x
      is_int = e == 0
   end function is_int

   !> Read a delimited file with one line of header
   subroutine ReadDelimFile(Filename, nCol, Array, errStat, errMsg, nHeaderLines, priPath)
      character(len=*),                        intent(in)  :: Filename
      integer,                                 intent(in)  :: nCol
      real(ReKi), dimension(:,:), allocatable, intent(out) :: Array
      integer(IntKi)         ,                 intent(out) :: errStat ! Status of error message
      character(*)           ,                 intent(out) :: errMsg  ! Error message if ErrStat /= ErrID_None
      integer(IntKi), optional,                intent(in ) :: nHeaderLines
      character(*)  , optional,                intent(in ) :: priPath  ! Primary path, to use if filename is not absolute
      integer              :: UnIn, i, j, nLine, nHead
      character(len= 2048) :: line
      integer(IntKi)       :: errStat2      ! local status of error message
      character(ErrMsgLen) :: errMsg2       ! temporary Error message
      character(len=2048) :: Filename_Loc   ! filename local to this function
      ErrStat = ErrID_None
      ErrMsg  = ""

      Filename_Loc = Filename
      if (present(priPath)) then
         if (PathIsRelative(Filename_Loc)) Filename_Loc = trim(PriPath)//trim(Filename)
      endif


      ! Open file
      call GetNewUnit(UnIn) 
      call OpenFInpFile(UnIn, Filename_Loc, errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'ReadDelimFile')
      if (errStat >= AbortErrLev) return
      ! Count number of lines
      nLine = line_count(UnIn)
      allocate(Array(nLine-1, nCol), stat=errStat2); errMsg2='allocation failed'; call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'ReadDelimFile')
      if (errStat >= AbortErrLev) return
      ! Read header
      nHead=1
      if (present(nHeaderLines)) nHead = nHeaderLines
      do i=1,nHead
         read(UnIn, *, IOSTAT=errStat2) line
         errMsg2 = ' Error reading line '//trim(Num2LStr(1))//' of file: '//trim(Filename_Loc)
         call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'ReadDelimFile')
         if (errStat >= AbortErrLev) return
      enddo
      ! Read data
      do I = 1,nLine-1
         read (UnIn,*,IOSTAT=errStat2) (Array(I,J), J=1,nCol)
         errMsg2 = ' Error reading line '//trim(Num2LStr(I+1))//' of file: '//trim(Filename_Loc)
         call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'ReadDelimFile')
         if (errStat >= AbortErrLev) return
      end do  
      close(UnIn) 
   end subroutine ReadDelimFile

   !> Counts number of lines in a file
   integer function line_count(iunit)
      integer, intent(in) :: iunit
      character(len=2048) :: line
      ! safety for infinite loop..
      integer :: i
      integer, parameter :: nline_max=100000000 ! 100 M
      line_count=0
      do i=1,nline_max 
         line=''
         read(iunit,'(A)',END=100)line
         line_count=line_count+1
      enddo
      if (line_count==nline_max) then
         print*,'Error: maximum number of line exceeded for line_count'
         STOP
      endif
   100 if(len(trim(line))>0) then
         line_count=line_count+1
      endif
      rewind(iunit)
      return
    end function
   !> Perform linear interpolation of an array, where first column is assumed to be ascending time values
   !! First value is used for times before, and last value is used for time beyond
   subroutine interpTimeValue(array, time, iLast, values)
      real(ReKi), dimension(:,:), intent(in)    :: array !< vector of time steps
      real(DbKi),                 intent(in)    :: time  !< time
      integer,                    intent(inout) :: iLast
      real(ReKi), dimension(:),   intent(out)   :: values !< vector of values at given time
      integer :: i
      real(ReKi) :: alpha
      if (array(iLast,1)> time) then 
         values = array(iLast,2:)
      elseif (iLast == size(array,1)) then 
         values = array(iLast,2:)
      else
         ! Look for index
         do i=iLast,size(array,1)
            if (array(i,1)<=time) then
               iLast=i
            else
               exit
            endif
         enddo
         if (iLast==size(array,1)) then
            values = array(iLast,2:)
         else
            ! Linear interpolation
            alpha = (array(iLast+1,1)-time)/(array(iLast+1,1)-array(iLast,1))
            values = array(iLast,2:)*alpha + array(iLast+1,2:)*(1-alpha)
            !print*,'time', array(iLast,1), '<=', time,'<',  array(iLast+1,1), 'fact', alpha
         endif
      endif
   end subroutine interpTimeValue
!----------------------------------------------------------------------------------------------------------------------------------
END PROGRAM 
