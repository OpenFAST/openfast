!**********************************************************************************************************************************
! UnsteadyAero_DriverCode: This code tests a stand-alone version of the UnsteadyAero module
!..................................................................................................................................
! LICENSING
! Copyright (C) 2015  National Renewable Energy Laboratory
!
!    This file is part of UnsteadyAero.
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


   
   
program UnsteadyAero_Driver

   use NWTC_Library
   use AirfoilInfo
   use AirfoilInfo_Types
   use UnsteadyAero_Types
   use UnsteadyAero
   use UA_Dvr_Subs
   USE VersionInfo

   implicit none

   
   
   
   
    ! Variables
   integer(IntKi), parameter                     :: NumInp = 2           ! Number of inputs sent to UA_UpdateStates (must be at least 2)
   
   real(DbKi)  :: dt, t, uTimes(NumInp)
   integer     :: i, j, n, iu
   type(UA_InitInputType)                        :: InitInData           ! Input data for initialization
   type(UA_InitOutputType)                       :: InitOutData          ! Output data from initialization
   type(UA_ContinuousStateType)                  :: x                    ! Continuous states
   type(UA_DiscreteStateType)                    :: xd                   ! Discrete states
   type(UA_OtherStateType)                       :: OtherState           ! Other/optimization states
   type(UA_MiscVarType)                          :: m                    ! Misc/optimization variables
   type(UA_ParameterType)                        :: p                    ! Parameters
   type(UA_InputType)                            :: u(NumInp)            ! System inputs
   type(UA_OutputType)                           :: y                    ! System outputs
   integer(IntKi)                                :: ErrStat              ! Status of error message
   character(ErrMsgLen)                          :: ErrMsg               ! Error message if ErrStat /= ErrID_None
   
   integer, parameter                            :: NumAFfiles = 1
   character(1024)                               :: afNames(NumAFfiles)
   type(AFI_ParameterType)                       :: AFI_Params(NumAFfiles)
   integer, allocatable                          :: AFIndx(:,:)
   CHARACTER(1024)                               :: dvrFilename          ! Filename and path for the driver input file.  This is passed in as a command line argument when running the Driver exe.
   TYPE(UA_Dvr_InitInput)                        :: dvrInitInp           ! Initialization data for the driver program
   real(DbKi)                                    :: simTime  
   integer                                       :: nSimSteps
   character(*), parameter                       :: RoutineName = 'UnsteadyAero_Driver'
   real(DbKi), allocatable                       :: timeArr(:)
   real(ReKi), allocatable                       :: AOAarr(:)
   real(ReKi), allocatable                       :: Uarr(:)
   real(ReKi), allocatable                       :: OmegaArr(:)
   
   CHARACTER(200)                                :: git_commit
   TYPE(ProgDesc), PARAMETER   :: version   = ProgDesc( 'UnsteadyAero Driver', '', '' )  ! The version number of this program.
      ! Initialize the NWTC library
   call NWTC_Init()
   
      ! Initialize error handling variables
   ErrMsg  = ''
   ErrStat = ErrID_None
   
   
      ! Display the copyright notice
   CALL DispCopyrightLicense( version%Name )
      ! Obtain OpenFAST git commit hash
   git_commit = QueryGitVersion()
      ! Tell our users what they're running
   CALL WrScr( ' Running '//TRIM( version%Name )//' a part of OpenFAST - '//TRIM(git_Commit)//NewLine//' linked with '//TRIM( NWTC_Ver%Name )//NewLine )
   
   
      ! Parse the driver file if one was provided, if not, then set driver parameters using hardcoded values
   if ( command_argument_count() > 1 ) then
      call print_help()
      call checkError()
   end if
  
   
      ! Establish initialization inputs which are fixed for the stand-alone driver, but would be
      ! variable for a coupled simulation
   InitInData%nNodesPerBlade  = 1 
   InitInData%numBlades       = 1
   
      ! Set up initialization data
   allocate(AFIndx(InitInData%nNodesPerBlade,InitInData%numBlades), STAT = ErrStat)
      if ( ErrStat /= 0 ) then
         call SetErrStat( ErrID_Fatal, 'Error trying to allocate InitInData%AFIndx.', ErrStat, ErrMsg, RoutineName)
         call checkError()
      end if
   
   allocate(InitInData%c(InitInData%nNodesPerBlade,InitInData%numBlades), STAT = ErrStat)
      if ( ErrStat /= 0 ) then
         call SetErrStat( ErrID_Fatal, 'Error trying to allocate InitInData%c.', ErrStat, ErrMsg, RoutineName)
         call checkError()
      end if
   
   allocate( InitInData%UAOff_innerNode(InitInData%numBlades), InitInData%UAOff_outerNode(InitInData%numBlades), STAT = ErrStat)
      if ( ErrStat /= 0 ) then
         call SetErrStat( ErrID_Fatal, 'Error trying to allocate UAOff_innerNode and UAOff_outerNode.', ErrStat, ErrMsg, RoutineName)
         call checkError()
      end if
   ! don't turn off UA based on span location:
   InitInData%UAOff_innerNode = 0
   InitInData%UAOff_outerNode = InitInData%nNodesPerBlade + 1
      
      ! Parse the driver input file and run the simulation based on that file
      
   if ( command_argument_count() == 1 ) then
      
      call get_command_argument(1, dvrFilename)
      call ReadDriverInputFile( dvrFilename, dvrInitInp, errStat, errMsg )
         call checkError()
      InitInData%a_s          = dvrInitInp%SpdSound
      InitInData%c(1,1)       = dvrInitInp%Chord
      InitInData%UAMod        = dvrInitInp%UAMod 
      InitInData%Flookup      = dvrInitInp%Flookup
   
   else
      
      dvrInitInp%OutRootName  = './TestingUA_Driver'
      InitInData%UAMod        = 1  
      InitInData%Flookup      = .FALSE.
      InitInData%a_s          = 340.29 ! m/s  
      InitInData%c(1,1)       = 1.0
      
      dvrInitInp%InflowVel    = 30.0 ! m/s
      dvrInitInp%Re           = 75  ! million
      dvrInitInp%AirFoil1     = './OSU075_FAST.txt'
      dvrInitInp%SimMod       = 1
      dvrInitInp%NCycles      = 3.0
      dvrInitInp%Frequency    = 1.2 ! Hz
      dvrInitInp%StepsPerCycle= 180
      dvrInitInp%Amplitude    = 10.0 ! deg
      dvrInitInp%Mean         =  2.0 ! deg
      dvrInitInp%Phase        =  0 ! steps of a cycle
      dvrInitInp%InputsFile   = ''
      
   end if
   InitInData%OutRootName = dvrInitInp%OutRootName
   
   InitInData%WrSum = .true. ! write all the AFI data

   
   if ( dvrInitInp%SimMod == 1 ) then
         ! Using the frequency and NCycles, determine how long the simulation needs to run
      simTime   = dvrInitInp%NCycles/dvrInitInp%Frequency
      nSimSteps = dvrInitInp%StepsPerCycle*dvrInitInp%NCycles  ! we could add 1 here to make this a complete cycle
      dt        = simTime / nSimSteps
      
   else
         ! Read time-series data file with a 1 line header and then each row contains time-step data with 4, white-space-separated columns
         ! time,  Angle-of-attack, Vrel, omega 
      call ReadTimeSeriesData( dvrInitInp%InputsFile, nSimSteps, timeArr, AOAarr, Uarr, OmegaArr, errStat, errMsg )
         call checkError()
      dt = (timeArr(nSimSteps) - timeArr(1)) / (nSimSteps-1)
      nSimSteps = nSimSteps-NumInp + 1
      
   end if
     
      ! Initialize the Airfoil Info Params
   afNames(1)  = dvrInitInp%AirFoil1 ! All nodes/blades are using the same 2D airfoil
   AFIndx(1,1) = 1
   call Init_AFI( InitInData%UAMod, NumAFfiles, afNames, dvrInitInp%UseCm, AFI_Params, errStat, errMsg )
      call checkError()

!   call WriteAFITables(AFI_Params(1), dvrInitInp%OutRootName)
   
   
    ! Initialize UnsteadyAero (after AFI)
   call UA_Init( InitInData, u(1), p, x, xd, OtherState, y, m, dt, AFI_Params, AFIndx, InitOutData, errStat, errMsg ) 
      call checkError()


   if (p%NumOuts <= 0) then
      ErrStat = ErrID_Fatal
      ErrMsg = "No outputs have been selected. Rebuild the executable with -DUA_OUTS"
      call checkError()
   end if

   ! set inputs:
   !u(1) = time at n=1  (t=   0)
   !u(2) = time at n=0  (t= -dt)
   !u(3) = time at n=-1 (t= -2dt) if NumInp > 2

   DO iu = 1, NumInp-1 !u(NumInp) is overwritten in time-sim loop, so no need to init here 
      call setUAinputs(2-iu,  u(iu), uTimes(iu), dt, dvrInitInp, timeArr, AOAarr, Uarr, OmegaArr)
   END DO
   
      ! Set inputs which do not vary with node or time

      ! time marching loop
   do n = 1, nSimSteps

      i = 1 ! nodes per blade
      j = 1 ! number of blades
     
      ! set inputs:
      DO iu = NumInp-1, 1, -1
         u(     iu+1) = u(     iu)
         uTimes(iu+1) = uTimes(iu)
      END DO
  
      ! first value of uTimes/u contain inputs at t+dt
      call setUAinputs(n+1,  u(1), uTimes(1), dt, dvrInitInp, timeArr, AOAarr, Uarr, OmegaArr)

      t = uTimes(2)

         ! Use existing states to compute the outputs
      call UA_CalcOutput(i, j, t, u(2),  p, x, xd, OtherState, AFI_Params(AFIndx(i,j)), y, m, errStat, errMsg )
         call checkError()
            
         ! Generate file outputs
      call UA_WriteOutputToFile(t, p, y)

      
         ! Prepare states for next time step
      call UA_UpdateStates(i, j, t, n, u, uTimes, p, x, xd, OtherState, AFI_Params(AFIndx(i,j)), m, errStat, errMsg )
         call checkError()
               
      
   end do
   
   
   !-------------------------------------------------------------------------------------------------
   ! Close our output file
   !-------------------------------------------------------------------------------------------------
   
   call Cleanup()
   call NormStop()
   
   contains
   
   !====================================================================================================
   subroutine Cleanup()
   !     The routine cleans up the module echo file and resets the NWTC_Library, reattaching it to 
   !     any existing echo information
   !----------------------------------------------------------------------------------------------------  
      call UA_End(p)
      
      ! probably should also deallocate driver variables here
      
   end subroutine Cleanup

   !----------------------------------------------------------------------------------------------------  
   subroutine checkError()
      
      if (ErrStat >= AbortErrLev) then
         
         call Cleanup()
         call ProgAbort(ErrMsg)
            
      elseif ( ErrStat /= ErrID_None ) then
         
         call WrScr( trim(ErrMsg) )
            
      end if
      
   end subroutine checkError
   !----------------------------------------------------------------------------------------------------  
   subroutine setUAinputs(n,u,t,dt,dvrInitInp,timeArr,AOAarr,Uarr,OmegaArr)
   
   integer,                intent(in)           :: n
   type(UA_InputType),     intent(inout)        :: u            ! System inputs
   real(DbKi),             intent(  out)        :: t
   real(DbKi),             intent(in)           :: dt
   TYPE(UA_Dvr_InitInput), intent(in)           :: dvrInitInp           ! Initialization data for the driver program
   real(DbKi),             intent(in)           :: timeArr(:)
   real(ReKi),             intent(in)           :: AOAarr(:)
   real(ReKi),             intent(in)           :: Uarr(:)
   real(ReKi),             intent(in)           :: OmegaArr(:)
   integer                                      :: indx
   real(ReKi)                                   :: phase

      u%UserProp = 0
      u%Re       = dvrInitInp%Re
   
      if ( dvrInitInp%SimMod == 1 ) then
         t       = (n-1)*dt
         phase = (n+dvrInitInp%Phase-1)*2*pi/dvrInitInp%StepsPerCycle
         u%alpha = (dvrInitInp%Amplitude * sin(phase) + dvrInitInp%Mean)*D2R   ! This needs to be in radians
 !        u%omega =  dvrInitInp%Amplitude * cos(phase) * dvrInitInp%Frequency * pi**2 / 90.0   ! This needs to be in radians derivative: d_alpha /d_t
         u%omega =  dvrInitInp%Amplitude * cos(phase) * 2*pi/dvrInitInp%StepsPerCycle / dt * D2R  ! This needs to be in radians derivative: d_alpha /d_t
         
         u%U     = dvrInitInp%InflowVel  ! m/s
      else
         indx = min(n,size(timeArr))
         indx = max(1, indx) ! use constant data at initialization
         
         ! Load timestep data from the time-series inputs which were previous read from input file
         t       = timeArr(indx)
         u%alpha = AOAarr(indx)*pi/180.0   ! This needs to be in radians
         u%omega = OmegaArr(indx)
         u%U     = Uarr(indx)
         if (n> size(timeArr)) then
            t = t + dt*(n - size(timeArr) ) ! update for NumInp>1;
         elseif (n < 1) then
            t = (n-1)*dt
         end if
      end if
      u%v_ac(1) = sin(u%alpha)*u%U
      u%v_ac(2) = cos(u%alpha)*u%U
   
   end subroutine setUAinputs
   !----------------------------------------------------------------------------------------------------  
   
   subroutine print_help()
    print '(a)', 'usage: '
    print '(a)', ''
    print '(a)', 'UnsteadyAero_Driver.exe [driverfilename]'
    print '(a)', ''
    print '(a)', 'Where the optional argument, driverfilename, is the name of the UnsteadyAero driver input file.'
    print '(a)', ''

   end subroutine print_help
   
end program UnsteadyAero_Driver

