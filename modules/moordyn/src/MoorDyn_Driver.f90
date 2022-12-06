!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2020-2021 Alliance for Sustainable Energy, LLC
! Copyright (C) 2015-2019 Matthew Hall
!
!    This file is part of MoorDyn.
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
PROGRAM MoorDyn_Driver

   USE MoorDyn_Types
   USE MoorDyn
   USE NWTC_Library 
   USE VersionInfo

   IMPLICIT NONE 

   TYPE MD_Drvr_InitInput
      LOGICAL                 :: Echo
      REAL(DbKi)              :: Gravity
      REAL(DbKi)              :: rhoW
      REAL(DbKi)              :: WtrDepth
      
      CHARACTER(1024)         :: MDInputFile
      CHARACTER(1024)         :: OutRootName
      REAL(DbKi)              :: TMax
      REAL(DbKi)              :: dtC
      
      INTEGER                 :: FarmSize
      REAL(DbKi)              :: FarmPositions(8,40)
      
      INTEGER                 :: InputsMod
      CHARACTER(1024)         :: InputsFile
      INTEGER                 :: nTurb
   END TYPE MD_Drvr_InitInput


   INTEGER(IntKi)                        :: ErrStat          ! Status of error message   
   CHARACTER(1024)                       :: ErrMsg           ! Error message if ErrStat /= ErrID_None

   INTEGER(IntKi)                        :: ErrStat2          ! Status of error message   
   CHARACTER(1024)                       :: ErrMsg2           ! Error message if ErrStat /= ErrID_None

   CHARACTER(1024)                       :: drvrFilename     ! Filename and path for the driver input file.  This is passed in as a command line argument when running the Driver exe.
   TYPE(MD_Drvr_InitInput)               :: drvrInitInp      ! Initialization data for the driver program
   INTEGER                               :: UnIn             ! Unit number for the input file
   INTEGER                               :: UnEcho           ! The local unit number for this module's echo file
  

   TYPE (MD_InitInputType)               :: MD_InitInp    
   TYPE (MD_ParameterType)               :: MD_p
   TYPE (MD_ContinuousStateType)         :: MD_x             ! continuous states
   TYPE (MD_InitOutputType)              :: MD_InitOut    
   TYPE (MD_DiscreteStateType)           :: MD_xd            ! discrete states
   TYPE (MD_ConstraintStateType)         :: MD_xc            ! constraint states
   TYPE (MD_OtherStateType)              :: MD_xo            ! other states
   TYPE (MD_MiscVarType)                 :: MD_m

   TYPE (MD_InputType),      ALLOCATABLE :: MD_u(:)
   REAL(DbKi), DIMENSION(:), ALLOCATABLE :: MD_uTimes

   TYPE (MD_OutputType)                  :: MD_y        ! Output file identifier

   ! Motion file parsing
   type(FileInfoType)                    :: FileInfo_PrescribeMtn  !< The derived type for holding the prescribed forces input file for parsing -- we may pass this in the future
   integer(IntKi)                        :: CurLine          !< current entry in FileInfo_In%Lines array
   real(ReKi), ALLOCATABLE               :: TmpRe(:)         !< temporary number array for reading values in

   CHARACTER(100)                        :: Line             ! String to temporarially hold value of read line
   REAL(ReKi), ALLOCATABLE               :: PtfmMotIn(:,:)   ! Variable for storing time, and DOF time series from driver input file
   REAL(ReKi), ALLOCATABLE               :: r_in(:,:)        ! Variable for storing interpolated DOF time series from driver input file
   REAL(ReKi), ALLOCATABLE               :: r_in2(:,:)       ! used for filtering
   REAL(ReKi), ALLOCATABLE               :: rd_in(:,:)       ! Variable for storing 1st derivative of interpolate DOF time series
   REAL(ReKi), ALLOCATABLE               :: rd_in2(:,:)      ! used for filtering
   REAL(ReKi), ALLOCATABLE               :: rdd_in(:,:)      ! Variable for storing 2nd derivative of interpolate DOF time series
   REAL(ReKi), ALLOCATABLE               :: rdd_in2(:,:)     ! used for filtering
   INTEGER(IntKi)                        :: ntIn             ! number of time steps read from driver input file
   INTEGER(IntKi)                        :: ncIn             ! number of channels read from driver input file
   INTEGER(IntKi)                        :: nt               ! number of coupling time steps to use in simulation

   REAL(DbKi)                            :: t                ! current time (s)
   REAL(DbKi)                            :: tMax             ! sim end time (s)
   REAL(DbKi)                            :: dtC              ! fixed/constant global time step
   REAL(DbKi)                            :: frac             ! fraction used in interpolation
         
   INTEGER(IntKi)                        :: MD_interp_order     ! order of interpolation/extrapolation

   ! Local variables
   Integer(IntKi)                        :: i,j,k,l              ! counter for various loops
   Integer(IntKi)                        :: iTurb
   Integer(IntKi)                        :: nTurbines
   Integer(IntKi)                        :: iIn
   integer(intKi)                        :: Un
  
   ! data for SimStatus/RunTimes:
   REAL(DbKi)                            :: PrevSimTime        !< Previous time message was written to screen (s > 0)
   REAL(ReKi)                            :: PrevClockTime      !< Previous clock time in seconds past midnight
   INTEGER                               :: SimStrtTime (8)    !< An array containing the elements of the start time (after initialization).
   INTEGER                               :: ProgStrtTime (8)   !< An array containing the elements of the program start time (before initialization).
   REAL(ReKi)                            :: SimStrtCPU         !< User CPU time for simulation (without intialization)
   REAL(ReKi)                            :: ProgStrtCPU        !< User CPU time for program (with intialization)

  
   CHARACTER(20)                         :: FlagArg              ! flag argument from command line
   CHARACTER(200)                        :: git_commit    ! String containing the current git commit hash
   TYPE(ProgDesc), PARAMETER             :: version = ProgDesc( 'MoorDyn Driver', '', '' )

  
  
   ErrMsg  = ""
   ErrStat = ErrID_None
   UnEcho=-1
   UnIn  =-1
  
   ! TODO: Sort out error handling (two sets of flags currently used)
  
   CALL NWTC_Init( ProgNameIn=version%Name )

   MD_InitInp%FileName = "MoorDyn.dat"  ! initialize to empty string to make sure it's input from the command line
   CALL CheckArgs( MD_InitInp%FileName, Arg2=drvrInitInp%InputsFile, Flag=FlagArg )
   IF ( LEN( TRIM(FlagArg) ) > 0 ) CALL NormStop()

      ! Display the copyright notice
   CALL DispCopyrightLicense( version%Name, 'Copyright (C) 2021 NREL, 2019 Matt Hall' )
      ! Obtain OpenFAST git commit hash
   git_commit = QueryGitVersion()
      ! Tell our users what they're running
   CALL WrScr( ' Running '//TRIM( version%Name )//' a part of OpenFAST - '//TRIM(git_commit)//NewLine//' linked with '//TRIM( NWTC_Ver%Name )//NewLine )


   
   CALL DATE_AND_TIME ( Values=ProgStrtTime )                        ! Let's time the whole simulation
   CALL CPU_TIME ( ProgStrtCPU )                                    ! Initial time (this zeros the start time when used as a MATLAB function)
   

   CALL WrScr( ' MD Driver updated 2022-01-12')

   ! Parse the driver input file and run the simulation based on that file
   CALL get_command_argument(1, drvrFilename)
   CALL ReadDriverInputFile( drvrFilename, drvrInitInp);
   
   ! do any initializing and allocating needed in prep for calling MD_Init   

   ! set the input file name and other environment terms
   !MD_InitInp%NStepWave   = 1        ! an arbitrary number > 0 (to set the size of the wave data, which currently contains all zero values)     
   MD_InitInp%Tmax                    = drvrInitInp%TMax   
   MD_InitInp%g                       = drvrInitInp%Gravity
   MD_InitInp%rhoW                    = drvrInitInp%rhoW
   MD_InitInp%WtrDepth                = drvrInitInp%WtrDepth
   MD_InitInp%FileName                = drvrInitInp%MDInputFile
   MD_InitInp%RootName                = drvrInitInp%OutRootName
   MD_InitInp%UsePrimaryInputFile     = .TRUE.
   !MD_InitInp%PassedPrimaryInputData  = 
   MD_InitInp%Echo                    = drvrInitInp%Echo
   !MD_InitInp%OutList                 = <<<< never used?
   MD_InitInp%Linearize               = .FALSE.
   
   TMax = drvrInitInp%TMax  
   dtC  = drvrInitInp%dtC                   ! desired coupling time step size for communicating with MoorDyn
  
   ! do OpenFAST vs FAST.Farm related setup
      
   MD_InitInp%FarmSize                = drvrInitInp%FarmSize
   
   if (drvrInitInp%FarmSize > 0) then   ! Check if this MoorDyn instance is being run from FAST.Farm (indicated by FarmSize > 0)
      nTurbines = drvrInitInp%FarmSize
   else    ! FarmSize==0 indicates normal, FAST module mode
      nTurbines = 1  ! if a regular FAST module mode, we treat it like a nTurbine=1 farm case
   end if
   
   CALL AllocAry(MD_InitInp%PtfmInit,      6, nTurbines, 'PtfmInit array'     , ErrStat2, ErrMsg2); call AbortIfFailed()
   CALL AllocAry(MD_InitInp%TurbineRefPos, 3, nTurbines, 'TurbineRefPos array', ErrStat2, ErrMsg2); call AbortIfFailed()
  
   do J=1,nTurbines
      MD_InitInp%TurbineRefPos(1,J) = drvrInitInp%FarmPositions(1,J)
      MD_InitInp%TurbineRefPos(2,J) = drvrInitInp%FarmPositions(2,J)
      MD_InitInp%TurbineRefPos(3,J) = 0.0_DbKi
      MD_InitInp%PtfmInit(1,J)      = drvrInitInp%FarmPositions(3,J)
      MD_InitInp%PtfmInit(2,J)      = drvrInitInp%FarmPositions(4,J)
      MD_InitInp%PtfmInit(3,J)      = drvrInitInp%FarmPositions(5,J)
      MD_InitInp%PtfmInit(4,J)      = drvrInitInp%FarmPositions(6,J)*D2R   !3.14159265/180.0
      MD_InitInp%PtfmInit(5,J)      = drvrInitInp%FarmPositions(7,J)*D2R   !3.14159265/180.0
      MD_InitInp%PtfmInit(6,J)      = drvrInitInp%FarmPositions(8,J)*D2R   !3.14159265/180.0
   end do
   
   MD_interp_order = 1
  
   ! allocate Input and Output arrays; used for interpolation and extrapolation
   Allocate(MD_uTimes(MD_interp_order + 1)) 
  
   ! @bonnie : This is in the FAST developers glue code example, but it's probably not needed here. 
   Allocate(MD_u(MD_interp_order + 1))
     

   if (drvrInitInp%InputsMod > 1) then
      ErrStat2 = ErrID_Fatal
      ErrMsg2  = ' ERROR: MoorDyn Driver InputsMod must be 0 or 1.'
      CALL AbortIfFailed()
   end if
   
   
   ! -------------------------------- -----------------------------------
   
   ! fill in the hydrodynamics data
   ALLOCATE( MD_InitInp%WaveVel (2,200,3))
   ALLOCATE( MD_InitInp%WaveAcc (2,200,3))
   ALLOCATE( MD_InitInp%WavePDyn(2,200)  )
   ALLOCATE( MD_InitInp%WaveElev(2,200)  )
   ALLOCATE( MD_InitInp%WaveTime(2)      )
   MD_InitInp%WaveVel  = 0.0_ReKi
   MD_InitInp%WaveAcc  = 0.0_ReKi
   MD_InitInp%WavePDyn = 0.0_ReKi
   MD_InitInp%WaveElev = 0.0_ReKi
   MD_InitInp%WaveTime = 0.0_ReKi
   DO I = 1,SIZE(MD_InitInp%WaveTime)
      MD_InitInp%WaveTime(I) = 600.0*I
   END DO
   
   ! open driver output file >>> not yet used <<<
   !CALL GetNewUnit( Un )
   !OPEN(Unit=Un,FILE='MD.out',STATUS='UNKNOWN')
  
   ! call the initialization routine
   CALL MD_Init( MD_InitInp, MD_u(1), MD_p, MD_x , MD_xd, MD_xc, MD_xo, MD_y, MD_m, dtC, MD_InitOut, ErrStat, ErrMsg2 ); call AbortIfFailed()
   
   CALL MD_DestroyInitInput  ( MD_InitInp , ErrStat, ErrMsg ); call AbortIfFailed()
   CALL MD_DestroyInitOutput ( MD_InitOut , ErrStat, ErrMsg ); call AbortIfFailed()
      
   CALL DispNVD( MD_InitOut%Ver ) 
   
   
   ! determine number of input channels expected from driver input file time series (DOFs including active tensioning channels)
   if (allocated(MD_u(1)%DeltaL)) then
      ncIn = size(MD_u(1)%DeltaL)      ! if unallocated, size will return garbage for some compilers
   else
      ncIn = 0
   endif
   
   do iTurb = 1, MD_p%nTurbines
      ncIn = ncIn + MD_p%nCpldBodies(iTurb)*6 + MD_p%nCpldRods(iTurb)*6 + MD_p%nCpldCons(iTurb)*3
   end do

   call WrScr('MoorDyn has '//trim(num2lstr(ncIn))//' coupled DOFs and/or active-tensioned inputs.')

   
   
   if (drvrInitInp%InputsMod == 1 ) then

      if ( LEN( TRIM(drvrInitInp%InputsFile) ) < 1 ) then
         ErrStat = ErrID_Fatal
         ErrMsg  = ' ERROR: MoorDyn Driver InputFile cannot be empty if InputsMode is 1.'
         CALL AbortIfFailed()
      end if
   
      call WrScr('Reading platform motion input data from '//trim(drvrInitInp%InputsFile))
      call WrScr('  MD driver is expecting '//trim(num2lstr(ncIn))//' columns of input data, plus time, in motion input file.')

      ! Parse the motion file and store in the FileInfoType structure.  This will strip out the header
      ! and leave just the table.  PrescribeMtn%NumLines is the number of timesteps.
      call ProcessComFile( drvrInitInp%InputsFile, FileInfo_PrescribeMtn, ErrStat2, ErrMsg2 )
      call AbortIfFailed()

      ! number of lines in table (number of timesteps)
      ntIn = FileInfo_PrescribeMtn%NumLines
 
      ! Allocate the array (include time column)
      call AllocAry( PtfmMotIn, ntIn, ncIn+1, "Array of motion data", ErrStat2, ErrMsg2 ); call AbortIfFailed()
      call AllocAry( TmpRe, ncIn+1, "TempRe", ErrStat2, ErrMsg2 ); call AbortIfFailed()
 
      ! Loop over all table lines.  Expecting ncIn+1 colunns
      CurLine=1
      do i=1,ntIn
         call ParseAry ( FileInfo_PrescribeMtn, CurLine, 'motions', TmpRe, ncIn+1, ErrStat2, ErrMsg2, UnEcho )
         ErrMsg2='Error reading the input time-series file. Expecting '//TRIM(Int2LStr(ncIn))//' channels plus time.'//NewLine//trim(ErrMsg2)
         call AbortIfFailed()
         PtfmMotIn(i,1:ncIn+1) = TmpRe
      enddo

      deallocate(TmpRe)

      call WrScr("Read "//trim(Num2LStr(ntIn))//" time steps from input file.")
      !print *, PtfmMotIn

      ! trim simulation duration to length of input file if needed
      if (PtfmMotIn(ntIn, 1) < TMax) then
         TMax = PtfmMotIn(ntIn, 1)
      end if   
      

  
      ! specify stepping details 
      nt = tMax/dtC - 1            ! number of coupling time steps

      
      ! allocate space for processed motion array
      ALLOCATE ( r_in(nt, ncIn), r_in2(nt, ncIn), rd_in(nt, ncIn), rd_in2(nt, ncIn), rdd_in(nt, ncIn), rdd_in2(nt, ncIn), STAT=ErrStat2)
      IF ( ErrStat2 /= ErrID_None ) THEN
         ErrStat2 = ErrID_Fatal
         ErrMsg  = '  Error allocating space for r_in or rd_in array.'
         call AbortIfFailed()
      END IF 


      ! go through and interpolate inputs to new regular time steps (if nt=0 this array should be left as zeros)
   
      DO i = 1,nt         
         t = dtC*(i-1)
         
         ! interpolation routine 
         DO iIn = 1,ntIn-1      
            IF (PtfmMotIn(iIn+1, 1) > t) THEN   ! find the right two points to interpolate between (remember that the first column of PtfmMotIn is time)
               frac = (t - PtfmMotIn(iIn, 1) )/( PtfmMotIn(iIn+1, 1) - PtfmMotIn(iIn, 1) )  ! interpolation fraction (0-1) between two interpolation points

               DO J=1,ncIn
                  ! get interpolated position of coupling point
                  r_in(i, J) = PtfmMotIn(iIn, J+1) + frac*(PtfmMotIn(iIn+1, J+1) - PtfmMotIn(iIn, J+1))
                  
                  if (iIn==1) then
                     ! use forward different to estimate velocity of coupling point
                     rd_in(i, J) = (PtfmMotIn(iIn+1, J+1) - PtfmMotIn(iIn, J+1)) / (PtfmMotIn(iIn+1, 1) - PtfmMotIn(iIn, 1))
                  else
                     ! use central different to estimate velocity of coupling point
                     rd_in(i, J) = (PtfmMotIn(iIn+1, J+1) - PtfmMotIn(iIn-1, J+1)) / (PtfmMotIn(iIn+1, 1) - PtfmMotIn(iIn-1, 1))

                  end if
               END DO
               
               EXIT   ! break out of the loop for this time step once we've done its interpolation
            END IF
         END DO
      
      END DO
      
      ! ----- filter position -----
      ! now filter forward
      DO i = 1,nt  
         DO J=1,ncIn
            if (i==1) then
               r_in2(i, J) = r_in(i, J)
            else
               r_in2(i, J) = 0.1*r_in(i, J) + 0.9*r_in2(i-1, J)
            end if
         END DO
      END DO
      ! now filter backward and save back to original variable
      DO i = nt,1,-1  
         DO J=1,ncIn
            if (i==nt) then
               r_in(i, J) = r_in2(i, J)
            else
               r_in(i, J) = 0.1*r_in2(i, J) + 0.9*r_in(i+1, J)
            end if
         END DO
      END DO
      
      
      ! now get derivative after filtering has been applied (derivative no longer needs to be calculated earlier) 
      DO i = 1,nt     
         DO J=1,ncIn
            if (i==1) then
               ! use forward different to estimate velocity of coupling point
               rd_in(i, J) = (r_in(i+1, J) - r_in(i, J)) / dtC
            else if (i==nt) then
               ! use forward different to estimate velocity of coupling point
               rd_in(i, J) = (r_in(i, J) - r_in(i-1, J)) / dtC
            else
               ! use central different to estimate velocity of coupling point
               rd_in(i, J) = (r_in(i+1, J) - r_in(i-1, J)) / (2.0*dtC)
            end if
         END DO
      END DO
      
      
      
      ! ----- filter velocity -----
      ! now filter forward
      DO i = 1,nt  
         DO J=1,ncIn
            if (i==1) then
               rd_in2(i, J) = rd_in(i, J)
            else
               rd_in2(i, J) = 0.1*rd_in(i, J) + 0.9*rd_in2(i-1, J)
            end if
         END DO
      END DO
      ! now filter backward and save back to original variable
      DO i = nt,1,-1  
         DO J=1,ncIn
            if (i==nt) then
               rd_in(i, J) = rd_in2(i, J)
            else
               rd_in(i, J) = 0.1*rd_in2(i, J) + 0.9*rd_in(i+1, J)
            end if
         END DO
      END DO
      
      
      ! now get derivative after filtering has been applied (derivative no longer needs to be calculated earlier) 
      DO i = 1,nt     
         DO J=1,ncIn
            if (i==1) then
               ! use forward different to estimate velocity of coupling point
               rdd_in(i, J) = (rd_in(i+1, J) - rd_in(i, J)) / dtC
            else if (i==nt) then
               ! use forward different to estimate velocity of coupling point
               rdd_in(i, J) = (rd_in(i, J) - rd_in(i-1, J)) / dtC
            else
               ! use central different to estimate velocity of coupling point
               rdd_in(i, J) = (rd_in(i+1, J) - rd_in(i-1, J)) / (2.0*dtC)
            end if
         END DO
      END DO
      
      
      ! ----- filter acceleration -----
      ! now filter forward
      DO i = 1,nt  
         DO J=1,ncIn
            if (i==1) then
               rdd_in2(i, J) = rdd_in(i, J)
            else
               rdd_in2(i, J) = 0.2*rdd_in(i, J) + 0.8*rdd_in2(i-1, J)
            end if
         END DO
      END DO
      ! now filter backward and save back to original variable
      DO i = nt,1,-1  
         DO J=1,ncIn
            if (i==nt) then
               rdd_in(i, J) = rdd_in2(i, J)
            else
               rdd_in(i, J) = 0.2*rdd_in2(i, J) + 0.8*rdd_in(i+1, J)
            end if
         END DO
      END DO
      
      
   else   
      nt = tMax/dtC - 1            ! number of coupling time steps
   end if   
   
   CALL WrScr(" ")
   call WrScr("Tmax - "//trim(Num2LStr(tMax))//" and nt="//trim(Num2LStr(nt)))
   CALL WrScr(" ")
   
   
   ! ---------------------------------------------------------------
   ! Set the initial input values
   ! ---------------------------------------------------------------
      
   ! zero the tension commands
   if (allocated(MD_u(1)%DeltaL)) then
      MD_u(1)%DeltaL = 0.0_ReKi
      MD_u(1)%DeltaLdot = 0.0_ReKi
   endif
   
!   ! zero water inputs (if passing wave info in from glue code)
!   MD_u(1)%U    = 0.0  
!   MD_u(1)%Ud   = 0.0  
!   MD_u(1)%zeta = 0.0  
!   MD_u(1)%PDyn = 0.0  
!   ! now add some current in x for testing
!   MD_u(1)%U(1,:) = 1.0
   
   ! copy inputs to initialize input arrays for higher interp orders if applicable
   DO i = 2, MD_interp_order + 1  
      CALL MD_CopyInput( MD_u(1), MD_u(i), MESH_NEWCOPY, ErrStat2, ErrMsg2 ); call AbortIfFailed()
   END DO  
   DO i = 1, MD_interp_order + 1  
       MD_uTimes(i) = -(i - 1) * dtC
   END DO

   ! get output at initialization (before time stepping)
   t = 0
   ! First, set correct inputs for initialization (errors occur otherwise) 
   if (drvrInitInp%InputsMod == 1 ) then
       
      DO iTurb = 1, MD_p%nTurbines
         i = 1  ! read first timestep data 
         K = 1  ! the index of the coupling points in the input mesh CoupledKinematics
         J = 1  ! the starting index of the relevant DOFs in the input array
         ! any coupled bodies (type -1)
         DO l = 1,MD_p%nCpldBodies(iTurb)
            MD_u(1)%CoupledKinematics(iTurb)%TranslationDisp(:,K) = r_in(i, J:J+2) - MD_u(1)%CoupledKinematics(iTurb)%Position(:,K) - MD_p%TurbineRefPos(:,iTurb)
            MD_u(1)%CoupledKinematics(iTurb)%Orientation(  :,:,K) = EulerConstruct( r_in(i, J+3:J+5) ) ! full Euler angle approach
            MD_u(1)%CoupledKinematics(iTurb)%TranslationVel( :,K) = rd_in(i, J:J+2)
            MD_u(1)%CoupledKinematics(iTurb)%RotationVel(    :,K) = rd_in(i, J+3:J+5)
            MD_u(1)%CoupledKinematics(iTurb)%TranslationAcc( :,K) = rdd_in(i, J:J+2)
            MD_u(1)%CoupledKinematics(iTurb)%RotationAcc(    :,K) = rdd_in(i, J+3:J+5)
         
            K = K + 1
            J = J + 6            
         END DO
         
         ! any coupled rods (type -1 or -2)    >>> need to make rotations ignored if it's a pinned rod <<<
         DO l = 1,MD_p%nCpldRods(iTurb)
         
            MD_u(1)%CoupledKinematics(iTurb)%TranslationDisp(:,K) = r_in(i, J:J+2) - MD_u(1)%CoupledKinematics(iTurb)%Position(:,K) - MD_p%TurbineRefPos(:,iTurb)
            MD_u(1)%CoupledKinematics(iTurb)%Orientation(  :,:,K) = EulerConstruct( r_in(i, J+3:J+5) )
            MD_u(1)%CoupledKinematics(iTurb)%TranslationVel( :,K) = rd_in(i, J:J+2)
            MD_u(1)%CoupledKinematics(iTurb)%RotationVel(    :,K) = rd_in(i, J+3:J+5)
            MD_u(1)%CoupledKinematics(iTurb)%TranslationAcc( :,K) = rdd_in(i, J:J+2)
            MD_u(1)%CoupledKinematics(iTurb)%RotationAcc(    :,K) = rdd_in(i, J+3:J+5)
         
            K = K + 1
            J = J + 6            
         END DO
         
         ! any coupled points (type -1)
         DO l = 1, MD_p%nCpldCons(iTurb)
            
            MD_u(1)%CoupledKinematics(iTurb)%TranslationDisp(:,K) = r_in(i, J:J+2) - MD_u(1)%CoupledKinematics(iTurb)%Position(:,K) - MD_p%TurbineRefPos(:,iTurb)
            MD_u(1)%CoupledKinematics(iTurb)%TranslationVel( :,K) = rd_in(i, J:J+2)
            MD_u(1)%CoupledKinematics(iTurb)%TranslationAcc( :,K) = 0.0_DbKi !rdd_in(i, J:J+2)
            
            !print *, u%PtFairleadDisplacement%Position(:,l) + u%PtFairleadDisplacement%TranslationDisp(:,l)
            !print *, u%PtFairleadDisplacement%TranslationVel(:,l)
            
            K = K + 1
            J = J + 3
         END DO
         
      end do  ! iTurb
      
      ! also provide any active tensioning commands
      if (allocated(MD_u(1)%DeltaL)) then
         do l = 1, size(MD_u(1)%DeltaL) 
            MD_u(1)%DeltaL(   l) = 0.0_DbKi ! r_in(i, J)
            MD_u(1)%DeltaLdot(l) = 0.0_DbKi !rd_in(i, J)
            J = J + 1         
         end do
      endif
   
   end if   ! InputsMod == 1 
   CALL MD_CalcOutput(  t, MD_u(1), MD_p, MD_x, MD_xd, MD_xc , MD_xo, MD_y, MD_m, ErrStat2, ErrMsg2 ); call AbortIfFailed()

  
  
  ! -------------------------------------------------------------------------
  ! BEGIN time marching 
  ! -------------------------------------------------------------------------
  
   call WrScr("Doing time marching now...")
   
   CALL SimStatus_FirstTime( PrevSimTime, PrevClockTime, SimStrtTime, SimStrtCPU, t, tMax )

   DO i = 1,nt

      ! --------------------------------- update inputs ---------------------------------

      t = dtC*(i-1)


      if ( MOD( i, 20 ) == 0 ) THEN         
         CALL SimStatus( PrevSimTime, PrevClockTime, t, tMax )
      end if
      
      ! shift older inputs back in the buffer
      CALL MD_CopyInput( MD_u(1), MD_u(2), MESH_NEWCOPY, ErrStat2, ErrMsg2 ); call AbortIfFailed()  ! copy from 1 to 2 before 1 is updated with latest.
      MD_uTimes(1) = t + dtC
      MD_uTimes(2) = MD_uTimes(1) - dtC 
      !MD_uTimes(3) = MD_uTimes(2) - dtC

      ! update coupled object kinematics iff we're reading input time series
      if (drvrInitInp%InputsMod == 1 ) then
         
         DO iTurb = 1, MD_p%nTurbines
            
            K = 1  ! the index of the coupling points in the input mesh CoupledKinematics
            J = 1  ! the starting index of the relevant DOFs in the input array
            ! any coupled bodies (type -1)
            DO l = 1,MD_p%nCpldBodies(iTurb)
               MD_u(1)%CoupledKinematics(iTurb)%TranslationDisp(:,K) = r_in(i, J:J+2) - MD_u(1)%CoupledKinematics(iTurb)%Position(:,K) - MD_p%TurbineRefPos(:,iTurb)
               MD_u(1)%CoupledKinematics(iTurb)%Orientation(  :,:,K) = EulerConstruct( r_in(i, J+3:J+5) ) ! full Euler angle approach
               MD_u(1)%CoupledKinematics(iTurb)%TranslationVel( :,K) = rd_in(i, J:J+2)
               MD_u(1)%CoupledKinematics(iTurb)%RotationVel(    :,K) = rd_in(i, J+3:J+5)
               MD_u(1)%CoupledKinematics(iTurb)%TranslationAcc( :,K) = rdd_in(i, J:J+2)
               MD_u(1)%CoupledKinematics(iTurb)%RotationAcc(    :,K) = rdd_in(i, J+3:J+5)
            
               K = K + 1
               J = J + 6            
            END DO
            
            ! any coupled rods (type -1 or -2)    >>> need to make rotations ignored if it's a pinned rod <<<
            DO l = 1,MD_p%nCpldRods(iTurb)
            
               MD_u(1)%CoupledKinematics(iTurb)%TranslationDisp(:,K) = r_in(i, J:J+2) - MD_u(1)%CoupledKinematics(iTurb)%Position(:,K) - MD_p%TurbineRefPos(:,iTurb)
               MD_u(1)%CoupledKinematics(iTurb)%Orientation(  :,:,K) = EulerConstruct( r_in(i, J+3:J+5) )
               MD_u(1)%CoupledKinematics(iTurb)%TranslationVel( :,K) = rd_in(i, J:J+2)
               MD_u(1)%CoupledKinematics(iTurb)%RotationVel(    :,K) = rd_in(i, J+3:J+5)
               MD_u(1)%CoupledKinematics(iTurb)%TranslationAcc( :,K) = rdd_in(i, J:J+2)
               MD_u(1)%CoupledKinematics(iTurb)%RotationAcc(    :,K) = rdd_in(i, J+3:J+5)
            
               K = K + 1
               J = J + 6            
            END DO
            
            ! any coupled points (type -1)
            DO l = 1, MD_p%nCpldCons(iTurb)
               
               MD_u(1)%CoupledKinematics(iTurb)%TranslationDisp(:,K) = r_in(i, J:J+2) - MD_u(1)%CoupledKinematics(iTurb)%Position(:,K) - MD_p%TurbineRefPos(:,iTurb)
               MD_u(1)%CoupledKinematics(iTurb)%TranslationVel( :,K) = rd_in(i, J:J+2)
               MD_u(1)%CoupledKinematics(iTurb)%TranslationAcc( :,K) = 0.0_DbKi !rdd_in(i, J:J+2)
               
               !print *, u%PtFairleadDisplacement%Position(:,l) + u%PtFairleadDisplacement%TranslationDisp(:,l)
               !print *, u%PtFairleadDisplacement%TranslationVel(:,l)
               
               K = K + 1
               J = J + 3
            END DO
            
         end do  ! iTurb
         
         ! also provide any active tensioning commands
         if (allocated(MD_u(1)%DeltaL)) then
            do l = 1, size(MD_u(1)%DeltaL) 
               MD_u(1)%DeltaL(   l) = 0.0_DbKi ! r_in(i, J)
               MD_u(1)%DeltaLdot(l) = 0.0_DbKi !rd_in(i, J)
               J = J + 1         
            end do
         endif
      
      end if   ! InputsMod == 1 
      
      ! >>> otherwise, mesh kinematics should all still be zero ... maybe worth checking <<<
      
      
      ! --------------------------------- update states ---------------------------------
      CALL MD_UpdateStates( t, nt, MD_u, MD_uTimes, MD_p, MD_x, MD_xd, MD_xc, MD_xo, MD_m, ErrStat2, ErrMsg2 ); call AbortIfFailed()
      
  
      ! update the global time step by one delta t               <<<< ??? why?  ADP: UpdateStates updtes from t -> t+dt.  Need to calculate outputs at this final time.
      t = t + dtC
     
      ! --------------------------------- calculate outputs ---------------------------------
      CALL MD_CalcOutput(  t, MD_u(1), MD_p, MD_x, MD_xd, MD_xc, MD_xo, MD_y, MD_m, ErrStat2, ErrMsg2 ); call AbortIfFailed()
     
     
      ! >>> should make output vector to hold and print outputs <<<
      !WRITE(Un, *) t, MD_u(1)%CoupledKinematics(1)%TranslationDisp(1,1), ((MD_y%CoupledLoads(1)%Force(k,j), k=1,3),j=1,3)
      !WRITE(*,*) t_global
     
      ! FORMAT(2(1X,F8.3),9(1X,E12.5)) 
     
   END DO
   
   
   ! -------------------------------------------------------------------------
   ! END time marching
   ! -------------------------------------------------------------------------
   
   CALL RunTimes( ProgStrtTime, ProgStrtCPU, SimStrtTime, SimStrtCPU, t )   
   
   ! Destroy all objects
   CALL MD_End( MD_u(1), MD_p, MD_x, MD_xd, MD_xc , MD_xo, MD_y, MD_m, ErrStat2, ErrMsg2 ); call AbortIfFailed()
   
   do j = 2,MD_interp_order+1
      call MD_DestroyInput( MD_u(j), ErrStat, ErrMsg)
   end do  
   
   DEALLOCATE(MD_u)
   DEALLOCATE(MD_uTimes)
   
   IF (ALLOCATED(r_in)  ) DEALLOCATE(r_in  )
   IF (ALLOCATED(PtfmMotIn)) DEALLOCATE(PtfmMotIn)
   
   CALL WrScr( "Program has ended" )
   close (un) 
  

CONTAINS

   SUBROUTINE AbortIfFailed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'MoorDyn_Driver') 
        IF ( ErrStat /= ErrID_None ) THEN
           CALL WrScr( "Local error: "//ErrMsg2 )
           CALL WrScr( "Full error messages: "//ErrMsg )
        END IF
        if (ErrStat >= AbortErrLev) then
           call CleanUp()
           STOP
        endif
   END SUBROUTINE AbortIfFailed

   LOGICAL FUNCTION Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'OutSummary') 
        Failed =  ErrStat >= AbortErrLev
        if (Failed) call CleanUp()
   END FUNCTION Failed

   SUBROUTINE CleanUp()
      if(UnEcho>0) CLOSE(UnEcho)
      if(UnEcho>0) CLOSE( UnIn)
      if(allocated(MD_u)) deallocate(MD_u)
   END SUBROUTINE CleanUp

   !-------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE ReadDriverInputFile( inputFile, InitInp)
      CHARACTER(*),                  INTENT( IN    )   :: inputFile
      TYPE(MD_Drvr_InitInput),       INTENT(   OUT )   :: InitInp
      ! Local variables  
      INTEGER                                          :: I                    ! generic integer for counting
      INTEGER                                          :: J                    ! generic integer for counting
      CHARACTER(   2)                                  :: strI                 ! string version of the loop counter

      CHARACTER(1024)                                  :: EchoFile             ! Name of MoorDyn echo file  
      CHARACTER(1024)                                  :: Line                 ! String to temporarially hold value of read line   
      CHARACTER(1024)                                  :: TmpPath              ! Temporary storage for relative path name
      CHARACTER(1024)                                  :: TmpFmt               ! Temporary storage for format statement
      CHARACTER(1024)                                  :: FileName             ! Name of MoorDyn input file  
      CHARACTER(1024)                                  :: FilePath             ! Path Name of MoorDyn input file  
   
      UnEcho=-1
      UnIn  =-1
   
      FileName = TRIM(inputFile)
   
      CALL GetNewUnit( UnIn )   
      CALL OpenFInpFile( UnIn, FileName, ErrStat2, ErrMsg2);
      call AbortIfFailed()
   
      CALL WrScr( 'Opening MoorDyn Driver input file:  '//FileName )
      
      ! Read until "echo"
      CALL ReadCom( UnIn, FileName, 'MoorDyn Driver input file header line 1', ErrStat2, ErrMsg2); call AbortIfFailed()
      CALL ReadCom( UnIn, FileName, 'MoorDyn Driver input file header line 2', ErrStat2, ErrMsg2); call AbortIfFailed()
      CALL ReadVar ( UnIn, FileName, InitInp%Echo, 'Echo', 'Echo Input', ErrStat2, ErrMsg2); call AbortIfFailed()
      ! If we echo, we rewind
      IF ( InitInp%Echo ) THEN
         EchoFile = TRIM(FileName)//'.echo'
         CALL GetNewUnit( UnEcho )   
         CALL OpenEcho ( UnEcho, EchoFile, ErrStat, ErrMsg ); call AbortIfFailed()
         REWIND(UnIn)
         CALL ReadCom( UnIn, FileName, 'MoorDyn Driver input file header line 1', ErrStat2, ErrMsg2, UnEcho); call AbortIfFailed()
         CALL ReadCom( UnIn, FileName, 'MoorDyn Driver input file header line 2', ErrStat2, ErrMsg2, UnEcho); call AbortIfFailed()
         CALL ReadVar ( UnIn, FileName, InitInp%Echo, 'Echo', 'Echo the input file data', ErrStat2, ErrMsg2, UnEcho); call AbortIfFailed()
      END IF
      !---------------------- ENVIRONMENTAL CONDITIONS -------------------------------------------------
      CALL ReadCom( UnIn, FileName, 'Environmental conditions header', ErrStat2, ErrMsg2, UnEcho); call AbortIfFailed()
      CALL ReadVar( UnIn, FileName, InitInp%Gravity, 'Gravity', 'Gravity', ErrStat2, ErrMsg2, UnEcho); call AbortIfFailed()
      CALL ReadVar( UnIn, FileName, InitInp%rhoW   , 'rhoW', 'water density', ErrStat2, ErrMsg2, UnEcho); call AbortIfFailed()
      CALL ReadVar( UnIn, FileName, InitInp%WtrDepth, 'WtrDepth', 'water depth', ErrStat2, ErrMsg2, UnEcho); call AbortIfFailed()
      !---------------------- MoorDyn -------------------------------------------------------------------
      CALL ReadCom( UnIn, FileName, 'MoorDyn header', ErrStat2, ErrMsg2, UnEcho); call AbortIfFailed()
      CALL ReadVar( UnIn, FileName, InitInp%MDInputFile, 'MDInputFile', 'MoorDyn input filename', ErrStat2, ErrMsg2, UnEcho); call AbortIfFailed()
      CALL ReadVar( UnIn, FileName, InitInp%OutRootName, 'OutRootName', 'MoorDyn output root filename', ErrStat2, ErrMsg2, UnEcho); call AbortIfFailed()
      CALL ReadVar( UnIn, FileName, InitInp%TMax       , 'Tmax', 'Simulation time duration', ErrStat2, ErrMsg2, UnEcho); call AbortIfFailed()
      CALL ReadVar( UnIn, FileName, InitInp%dtC        , 'dtC', 'Time step size for calling MoorDyn', ErrStat2, ErrMsg2, UnEcho); call AbortIfFailed()
      CALL ReadVar( UnIn, FileName, InitInp%InputsMod  , 'InputsMode', 'Mode for the inputs - zero/steady/time-series', ErrStat2, ErrMsg2, UnEcho); call AbortIfFailed()
      CALL ReadVar( UnIn, FileName, InitInp%InputsFile , 'InputsFile', 'Filename for the MoorDyn inputs', ErrStat2, ErrMsg2, UnEcho); call AbortIfFailed()
      CALL ReadVar( UnIn, FileName, InitInp%FarmSize   , 'NumTurbines', 'number of turbines in FAST.Farm', ErrStat2, ErrMsg2, UnEcho); call AbortIfFailed()
      CALL ReadCom( UnIn, FileName, 'Initial positions header', ErrStat2, ErrMsg2); call AbortIfFailed()
      CALL ReadCom( UnIn, FileName, 'Initial positions table header line 1', ErrStat2, ErrMsg2); call AbortIfFailed()
      CALL ReadCom( UnIn, FileName, 'Initial positions table header line 2', ErrStat2, ErrMsg2); call AbortIfFailed()
      do J=1,MAX(1,InitInp%FarmSize)
         CALL ReadAry( UnIn, FileName, InitInp%FarmPositions(:,J), 8, "FarmPositions", "FAST.Farm position inputs", ErrStat2, ErrMsg2, UnEcho); call AbortIfFailed()
      end do

      ! done reading
      if(UnEcho>0) CLOSE( UnEcho )
      if(UnIn>0)   CLOSE( UnIn   )
   
      ! Perform input checks and triggers
      !CALL GetPath( FileName, FilePath )
      !IF ( PathIsRelative( InitInp%MDInputFile ) ) then
      !   InitInp%MDInputFile = TRIM(FilePath)//TRIM(InitInp%MDInputFile)
      !END IF
      !IF ( PathIsRelative( InitInp%OutRootName ) ) then
      !   InitInp%OutRootName = TRIM(FilePath)//TRIM(InitInp%OutRootName)
      !endif
      !IF ( PathIsRelative( InitInp%InputsFile ) ) then
      !   InitInp%InputsFile = TRIM(FilePath)//TRIM(InitInp%InputsFile)
      !endif

   END SUBROUTINE ReadDriverInputFile

   subroutine print_help()
       print '(a)', 'usage: '
       print '(a)', ''
       print '(a)', 'MoorDynDriver.exe driverfilename'
       print '(a)', ''
       print '(a)', 'Where driverfilename is the name of the MoorDyn driver input file.'
       print '(a)', ''
   end subroutine print_help
!----------------------------------------------------------------------------------------------------------------------------------

END PROGRAM
