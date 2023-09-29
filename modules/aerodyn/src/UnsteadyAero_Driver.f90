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
   use VersionInfo

   use LinDyn

   implicit none
    ! Variables
   
   real(DbKi)  :: t
   integer     :: i, j, n, iu

   ! --- All Data
   type(Dvr_Data)         :: dvr
   TYPE(UA_Dvr_InitInput) :: dvrInitInp           ! Initialization data for the driver program

   integer(IntKi)                                :: ErrStat              ! Status of error message
   character(ErrMsgLen)                          :: ErrMsg               ! Error message if ErrStat /= ErrID_None
   
   CHARACTER(1024)                               :: dvrFilename          ! Filename and path for the driver input file.  This is passed in as a command line argument when running the Driver exe.
   character(*), parameter                       :: RoutineName = 'UnsteadyAero_Driver'
   
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
   
   
   ! --- Parse the driver file if one
   if ( command_argument_count() > 1 ) then
      call print_help()
      call NormStop()
   endif
   call get_command_argument(1, dvrFilename)
   call ReadDriverInputFile( dvrFilename, dvrInitInp, errStat, errMsg ); call checkError()

   ! --- Driver Data 
   dvr%out%Root = dvrInitInp%OutRootName

   ! --- Time simulation control
   if ( dvrInitInp%SimMod == 1 ) then
      ! Using the frequency and NCycles, determine how long the simulation needs to run
      dvr%simTime   = dvrInitInp%NCycles/dvrInitInp%Frequency
      dvr%numSteps = dvrInitInp%StepsPerCycle*dvrInitInp%NCycles  ! we could add 1 here to make this a complete cycle
      dvr%dt        = dvr%simTime / dvr%numSteps
      
   else if ( dvrInitInp%SimMod == 2 ) then
      ! Read time-series data file with columns:( time,  Angle-of-attack, Vrel, omega )
      call ReadTimeSeriesData( dvrInitInp%InputsFile, dvr%numSteps, dvr%timeArr, dvr%AOAarr, dvr%Uarr, dvr%OmegaArr, errStat, errMsg ); call checkError()
      dvr%dt = (dvr%timeArr(dvr%numSteps) - dvr%timeArr(1)) / (dvr%numSteps-1)
      dvr%numSteps = dvr%numSteps-NumInp + 1

   elseif ( dvrInitInp%SimMod == 3 ) then
      dvr%simTime   = dvrInitInp%TMax
      dvr%dt        = dvrInitInp%dt
      dvr%numSteps = int(dvr%simTime/dvr%dt) ! TODO

      ! --- Initialize Elastic Section
      call LD_InitInputData(3, dvr%LD_InitInData, errStat, errMsg); call checkError()
      dvr%LD_InitInData%dt        = dvr%dt
      dvr%LD_InitInData%IntMethod = 1  ! TODO
      dvr%LD_InitInData%prefix    = '' ! TODO for output channel names
      dvr%LD_InitInData%MM         = dvrInitInp%MM
      dvr%LD_InitInData%CC         = dvrInitInp%CC
      dvr%LD_InitInData%KK         = dvrInitInp%KK
      dvr%LD_InitInData%x0         = dvrInitInp%initPos
      dvr%LD_InitInData%xd0        = dvrInitInp%initVel
      dvr%LD_InitInData%activeDOFs = dvrInitInp%activeDOFs
      call LD_Init(dvr%LD_InitInData, dvr%LD_u(1), dvr%LD_p, dvr%LD_x, dvr%LD_xd, dvr%LD_z, dvr%LD_OtherState, dvr%LD_y, dvr%LD_m, dvr%LD_InitOutData, errStat, errMsg); call checkError()

      call Dvr_InitializeDriverOutputs(dvr, dvr%out, errStat, errMsg); call checkError()

   end if

   ! --- Init UA input data based on driver inputs
   call driverInputsToUAInitData(dvrInitInp, dvr%UA_InitInData, dvr%AFI_Params, dvr%AFIndx, errStat, errMsg); call checkError()

   ! --- Initialize UnsteadyAero (need AFI)
   call UA_Init( dvr%UA_InitInData, dvr%UA_u(1), dvr%UA_p, dvr%UA_x, dvr%UA_xd, dvr%UA_OtherState, dvr%UA_y, dvr%UA_m, dvr%dt, dvr%AFI_Params, dvr%AFIndx, dvr%UA_InitOutData, errStat, errMsg ); call checkError()
   if (dvr%UA_p%NumOuts <= 0) then
      ErrStat = ErrID_Fatal
      ErrMsg = "No outputs have been selected. Rebuild the executable with -DUA_OUTS"
      call checkError()
   end if


   ! --- Initialize Inputs
   !u(1) = time at n=1  (t=   0)
   !u(2) = time at n=0  (t= -dt)
   !u(3) = time at n=-1 (t= -2dt) if NumInp > 2
   if ( dvrInitInp%SimMod == 3 ) then
      ! General inputs
      do iu = 1, NumInp !u(NumInp) is overwritten in time-sim loop, so no need to init here 
         dvr%uTimes(iu) = (2-iu-1)*dvr%dt
      enddo
      ! LD Inputs - Allocs
      do iu = 2,NumInp
         call AllocAry(dvr%LD_u(iu)%Fext, dvr%LD_p%nx, 'Fext', errStat, errMsg); call checkError()
      enddo
      ! UA inputs:
      do iu = 1, NumInp-1 !u(NumInp) is overwritten in time-sim loop, so no need to init here 
         ! TODO TODO TODO
         dvr%UA_u(iu)%UserProp = 0
         dvr%UA_u(iu)%Re       = dvrInitInp%Re
         dvr%UA_u(iu)%omega    =  0.0_ReKi
         dvr%UA_u(iu)%v_ac(1)  = 0.0_ReKi
         dvr%UA_u(iu)%v_ac(2)  = 0.0_ReKi
         dvr%UA_u(iu)%alpha    = 0.0_ReKi
         dvr%UA_u(iu)%U        = 0.0_ReKi
      enddo
   else
      ! UA inputs:
      do iu = 1, NumInp-1 !u(NumInp) is overwritten in time-sim loop, so no need to init here 
         call setUAinputs(2-iu,  dvr%UA_u(iu), dvr%uTimes(iu), dvr%dt, dvrInitInp, dvr%timeArr, dvr%AOAarr, dvr%Uarr, dvr%OmegaArr, errStat, errMsg); call checkError()
      end do
   endif

   i = 1 ! nodes per blade
   j = 1 ! number of blades
   ! --- Time marching loop

   if ( dvrInitInp%SimMod == 3 ) then

      call Dvr_InitializeOutputs(dvr%out, dvr%numSteps, errStat, errMsg)

      dvr%LD_u(1)%Fext=0.0_ReKi ! TODO TODO
      dvr%LD_u(2)%Fext=0.0_ReKi ! TODO TODO

      ! --- time marching loop
      print*,'>>> Time simulation', dvr%uTimes(1), dvr%numSteps*dvr%dt
      do n = 1, dvr%numSteps
         ! set inputs:
         do iu = NumInp-1, 1, -1
            dvr%UA_u(     iu+1) = dvr%UA_u(     iu)
            dvr%LD_u(     iu+1) = dvr%LD_u(     iu)
            dvr%uTimes(iu+1) = dvr%uTimes(iu)
         end do
         !          ! first value of uTimes/u contain inputs at t+dt

         ! Basic inputs
         dvr%uTimes(1) = (n+1-1)*dvr%dt
         ! UA-LD Inputs Solve TODO TODO TODO
         !  call setUAinputs(n+1,  u(1), uTimes(1), dt, dvrInitInp, timeArr, AOAarr, Uarr, OmegaArr, errStat, errMsg); call checkError()
         dvr%UA_u(1)%UserProp = 0
         dvr%UA_u(1)%Re       = dvrInitInp%Re
         dvr%UA_u(1)%omega    = dvr%LD_x%q(6)
         dvr%UA_u(1)%v_ac(1)  = dvrInitInp%Mean -dvr%LD_x%q(4)
         dvr%UA_u(1)%v_ac(2)  =                 -dvr%LD_x%q(5)
         dvr%UA_u(1)%alpha    = 0.0_ReKi
         dvr%UA_u(1)%U        = sqrt(  dvr%UA_u(1)%v_ac(1)**2  +  dvr%UA_u(1)%v_ac(2)**2)

         t = dvr%uTimes(2)
         ! Use existing states to compute the outputs
         call LD_CalcOutput(t, dvr%LD_u(2), dvr%LD_p, dvr%LD_x, dvr%LD_xd, dvr%LD_z, dvr%LD_OtherState, dvr%LD_y, dvr%LD_m, errStat, errMsg); call checkError()
         !! Use existing states to compute the outputs
         call UA_CalcOutput(i, j, t, dvr%UA_u(2),  dvr%UA_p, dvr%UA_x, dvr%UA_xd, dvr%UA_OtherState, dvr%AFI_Params(dvr%AFIndx(i,j)), dvr%UA_y, dvr%UA_m, errStat, errMsg ); call checkError()

         dvr%LD_u(1)%Fext(1) = 0.5_ReKi * dvrInitInp%Chord * dvr%UA_u(1)%U**2  * dvr%UA_y%Cl /100   ! TODO TODO
         dvr%LD_u(1)%Fext(2) = 0.5_ReKi * dvrInitInp%Chord * dvr%UA_u(1)%U**2  * dvr%UA_y%Cd /100   ! TODO TODO
         !y%Cn
         !y%Cc
         !y%Cm
         !y%Cl
         !y%Cd



         ! Generate file outputs
         call UA_WriteOutputToFile(t, dvr%UA_p, dvr%UA_y)
         ! Write outputs for all turbines at nt-1
         call Dvr_WriteOutputs(n, t, dvr, dvr%out, errStat, errMsg); call checkError()


         ! Prepare states for next time step
         call LD_UpdateStates(t, n, dvr%LD_u, dvr%uTimes, dvr%LD_p, dvr%LD_x, dvr%LD_xd, dvr%LD_z, dvr%LD_OtherState, dvr%LD_m, errStat, errMsg); call checkError()
         ! Prepare states for next time step
         call UA_UpdateStates(i, j, t, n, dvr%UA_u, dvr%uTimes, dvr%UA_p, dvr%UA_x, dvr%UA_xd, dvr%UA_OtherState, dvr%AFI_Params(dvr%AFIndx(i,j)), dvr%UA_m, errStat, errMsg ); call checkError()
      end do

      call Dvr_EndSim(dvr, errStat, errMsg)
   else
      ! --- time marching loop
      do n = 1, dvr%numSteps
        
         ! set inputs:
         DO iu = NumInp-1, 1, -1
            dvr%UA_u(  iu+1) = dvr%UA_u(     iu)
            dvr%uTimes(iu+1) = dvr%uTimes(iu)
         END DO
     
         ! first value of uTimes/u contain inputs at t+dt
         call setUAinputs(n+1,  dvr%UA_u(1), dvr%uTimes(1), dvr%dt, dvrInitInp, dvr%timeArr, dvr%AOAarr, dvr%Uarr, dvr%OmegaArr, errStat, errMsg); call checkError()
           
         t = dvr%uTimes(2)

         ! Use existing states to compute the outputs
         call UA_CalcOutput(i, j, t, dvr%UA_u(2), dvr%UA_p, dvr%UA_x, dvr%UA_xd, dvr%UA_OtherState, dvr%AFI_Params(dvr%AFIndx(i,j)), dvr%UA_y, dvr%UA_m, errStat, errMsg ); call checkError()
               
         ! Generate file outputs
         call UA_WriteOutputToFile(t, dvr%UA_p, dvr%UA_y)
         
         ! Prepare states for next time step
         call UA_UpdateStates(i, j, t, n, dvr%UA_u, dvr%uTimes, dvr%UA_p, dvr%UA_x, dvr%UA_xd, dvr%UA_OtherState, dvr%AFI_Params(dvr%AFIndx(i,j)), dvr%UA_m, errStat, errMsg ); call checkError()
         
      end do
   endif
   
   ! --- Exit
   call Cleanup()
   call NormStop()

contains
   
   !====================================================================================================
   subroutine Cleanup()
   !     The routine cleans up the module echo file and resets the NWTC_Library, reattaching it to 
   !     any existing echo information
   !----------------------------------------------------------------------------------------------------  
      call UA_End(dvr%UA_p)
      
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

   !----------------------------------------------------------------------------------------------------  
   subroutine setUAinputs(n,u,t,dt,dvrInitInp,timeArr,AOAarr,Uarr,OmegaArr,errStat,errMsg)
   
   integer,                intent(in)              :: n
   type(UA_InputType),     intent(inout)           :: u            ! System inputs
   real(DbKi),             intent(  out)           :: t
   real(DbKi),             intent(in)              :: dt
   TYPE(UA_Dvr_InitInput), intent(in)              :: dvrInitInp           ! Initialization data for the driver program
   real(DbKi),             intent(in), allocatable :: timeArr(:)
   real(ReKi),             intent(in), allocatable :: AOAarr(:)
   real(ReKi),             intent(in), allocatable :: Uarr(:)
   real(ReKi),             intent(in), allocatable :: OmegaArr(:)
   integer,                intent(out)             :: errStat
   character(len=*),       intent(out)             :: errMsg
   integer                                         :: indx
   real(ReKi)                                      :: phase
   real(ReKi)                                      :: d_ref2AC
   real(ReKi)                                      :: alpha_ref
   real(ReKi)                                      :: U_ref
   real(ReKi)                                      :: v_ref(2)
   real(ReKi)                                      :: v_34(2)
   logical, parameter :: OscillationAtMidChord=.true.  ! for legacy, use false
   logical, parameter :: VelocityAt34         =.true.  ! for legacy, use false

      ! Initialize error handling variables
      ErrMsg  = ''
      ErrStat = ErrID_None

      u%UserProp = 0
      u%Re       = dvrInitInp%Re
   
      if ( dvrInitInp%SimMod == 1 ) then
         if (OscillationAtMidChord) then
            d_ref2AC =-0.25_ReKi  ! -0.25: oscillations at mid_chord
         else
            d_ref2AC = 0.0_ReKi   ! 0: oscillations at AC
         endif
         U_ref = dvrInitInp%InflowVel  ! m/s

         t       = (n-1)*dt
         phase = (n+dvrInitInp%Phase-1)*2*pi/dvrInitInp%StepsPerCycle
         alpha_ref = (dvrInitInp%Amplitude * sin(phase) + dvrInitInp%Mean)*D2R   ! This needs to be in radians
         v_ref(1) = sin(alpha_ref)*U_ref
         v_ref(2) = cos(alpha_ref)*U_ref
         u%omega =  dvrInitInp%Amplitude * cos(phase) * 2*pi/dvrInitInp%StepsPerCycle / dt * D2R  ! This needs to be in radians derivative: d_alpha /d_t

         u%v_ac(1) = v_ref(1) + u%omega * d_ref2AC* dvrInitInp%Chord
         u%v_ac(2) = v_ref(2)

         v_34(1) = u%v_ac(1) + u%omega * 0.5* dvrInitInp%Chord
         v_34(2) = u%v_ac(2)


         u%alpha = atan2(u%v_ac(1), u%v_ac(2) )  ! 
         if (VelocityAt34) then
            u%U =  sqrt(v_34(1)**2 + v_34(2)**2) ! Using U at 3/4
         else
            u%U =  sqrt(u%v_ac(1)**2 + u%v_ac(2)**2) ! Using U at 1/4
         endif


      else
         ! check optional variables and allocation status
         if (all( (/ allocated(timeArr),allocated(AOAarr),allocated(OmegaArr),allocated(Uarr) /) )) then
             
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
            u%v_ac(1) = sin(u%alpha)*u%U
            u%v_ac(2) = cos(u%alpha)*u%U
         else
            errStat = ErrID_Fatal
            errMsg = 'mandatory input arrays are not allocated: timeArr,AOAarr,OmegaArr,Uarr'
         end if
             
      end if
   
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

