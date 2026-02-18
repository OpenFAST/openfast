!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2015-2016  National Renewable Energy Laboratory
!
!    This file is part of AeroDyn.
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
 
!there will also be various control flags... this may be updated as needed:
!TBLflag = {'BPM','TNO'}
!bluntnessFlag = {'DTU','BPM'}
!BPMBLflag = {'true','false'}
!useOrigModelAtSepOnset = {'true','false'}



   
   
!**********************************************************************************************************************************
program AeroAcoustics_Driver
   use AeroAcoustics_Driver_Subs
   use VersionInfo
   implicit none
   
   ! Program variables
   REAL(ReKi)                       :: PrevClockTime ! Clock time at start of simulation in seconds [(s)]
   REAL(ReKi)                       :: UsrTime1      ! User CPU time for simulation initialization [(s)]
   REAL(ReKi)                       :: UsrTime2      ! User CPU time for simulation (without initialization) [(s)]
   INTEGER(IntKi) , DIMENSION(1:8)  :: StrtTime      ! Start time of simulation (including initialization) [-]
   INTEGER(IntKi) , DIMENSION(1:8)  :: SimStrtTime   ! Start time of simulation (after initialization) [-]
   REAL(DbKi)                       :: t_global      ! global-loop time marker
   REAL(DbKi)                       :: TiLstPrn      ! The simulation time of the last print (to file) [(s)]
   
   TYPE(Dvr_Data)                   :: DriverData
   
   character(1024)                  :: InputFile
   integer                          :: nt            !< loop counter (for time step)
   character(20)                    :: FlagArg       ! flag argument from command line
   integer(IntKi)                   :: ErrStat       ! status of error message
   character(ErrMsgLen)             :: ErrMsg        !local error message if ErrStat /= ErrID_None

   
   CALL DATE_AND_TIME ( Values=StrtTime )                 ! Let's time the whole simulation
   CALL CPU_TIME ( UsrTime1 )                             ! Initial time (this zeros the start time when used as a MATLAB function)
   UsrTime1 = MAX( 0.0_ReKi, UsrTime1 )                   ! CPU_TIME: If a meaningful time cannot be returned, a processor-dependent negative value is returned
   UsrTime2 = UsrTime1                                    ! CPU_TIME: Initialize in case of error before getting real data
   SimStrtTime = StrtTime                                    ! CPU_TIME: Initialize in case of error before getting real data
   nt = 0
   
   ! --- Driver initialization
   CALL NWTC_Init( ProgNameIN=version%Name )
   
   InputFile = ""  ! initialize to empty string to make sure it's input from the command line
   CALL CheckArgs( InputFile, Flag=FlagArg )
   IF ( LEN( TRIM(FlagArg) ) > 0 ) CALL NormStop()
   
   ! Display the copyright notice and compile info:
   CALL DispCopyrightLicense( version%Name )
   CALL DispCompileRuntimeInfo( version%Name )


   ! Initialize modules
   call ReadDriverInputFile( InputFile, DriverData, ErrStat, ErrMsg ); call CheckError()
   call Init_AFI(DriverData%Airfoil_FileName, DriverData%AFInfo, ErrStat, ErrMsg); call CheckError()
   call Init_AAmodule(DriverData, ErrStat, ErrMsg); call CheckError()

   ! Init of time estimator
   t_global=0.0_DbKi
   call SimStatus_FirstTime( TiLstPrn, PrevClockTime, SimStrtTime, UsrTime2, t_global, DriverData%TMax )

   ! Time loop
   do nt = 1, DriverData%numSteps
      ! Time update to screen
      t_global=nt * DriverData%dt
      
      if (mod( nt + 1, 10 )==0) call SimStatus(TiLstPrn, PrevClockTime, t_global, DriverData%TMax)
      
      ! update states and calculate output
      call SetInputsForAA(DriverData)
      
      call AA_CalcOutput(t_global, DriverData%u, DriverData%p, DriverData%xd, DriverData%OtherState, DriverData%y, DriverData%m, errStat, errMsg); call CheckError()
      call Dvr_WriteOutputs(t_global, nt, DriverData) ! write to file at this step

   ! Get state variables at next step: INPUT at step nt - 1, OUTPUT at step nt
      call AA_UpdateStates(t_global, nt, DriverData%m, DriverData%u, DriverData%p, DriverData%xd, DriverData%OtherState, errStat, errMsg); call CheckError()
      
   end do !nt=1,numSteps


   call Dvr_End()
contains
!................................   
   subroutine CheckError()
      if (ErrStat /= ErrID_None) then
         call WrScr(TRIM(errMsg))
         if (errStat >= AbortErrLev) then
            call Dvr_End()
         end if
         ErrStat = ErrID_None
      end if
   end subroutine CheckError
!................................   
   subroutine Dvr_End()
      integer(IntKi)       :: errStat2      ! local status of error message
      character(ErrMsgLen) :: errMsg2       ! local error message if ErrStat /= ErrID_None

      call Dvr_EndOutput(DriverData, nt, errStat2, errMsg2)
        if (ErrStat2 /= ErrID_None) call WrScr(TRIM(errMsg2))
      
      call RunTimes(StrtTime, UsrTime1, SimStrtTime, UsrTime2, t_global)
        
      if (ErrStat >= AbortErrLev) then
         call WrScr('')
         CALL ProgAbort( 'AeroAcoustics Driver encountered simulation error level: '&
             //TRIM(GetErrStr(ErrStat)), TrapErrors=.FALSE., TimeWait=3._ReKi )  ! wait 3 seconds (in case they double-clicked and got an error)
      else
         call NormStop()
      end if
   end subroutine Dvr_End
!................................
end program AeroAcoustics_Driver


   
!Inputs that will be supplied externally:
!driver%DT
!
!Need to set in InitInput: 
!   rho [InitInputType%airDens]
!   c0 or co [InitInputType%SpdSound]
!   L [InitInputType%BlSpn]
!   chord [InitInputType%BlChord]
!   [driver%DT = Interval]
!   visc [InitInputType%KinVisc]
!   [InitInputType%HubHeight]
!   Airfoil info:
!      BlAFID
!      AFInfo
!   
!Set in AA Input File:
!   Lturb (already in AA input file) [InputFileData%Lturb]
!   dStarS [m%dstarVar(1), dstarVar1, DSTRS -> interpolated from p%dstarall1 = InputFileData%Suct_DispThick using AoA and Re] meters
!   dStarP [m%dstarVar(2), dstarVar2, DSTRP -> interpolated from p%dstarall2 = InputFileData%Pres_DispThick using AoA and Re] meters
!   TI [InputFileData%TI]
!   cfS [m%CfVar(1), Cfall(1) -> interpolated p%Cfall1=InputFileData%Suct_Cf with Re and AoA]
!   cfP [m%CfVar(2), Cfall(2) -> interpolated p%Cfall2=InputFileData%Pres_Cf with Re and AoA]
!   deltaS [m%d99Var(1) -> interpolated p%d99all1=InputFileData%Suct_BLThick with Re and AoA]
!   deltaP [m%d99Var(2), d99Var2 -> interpolated p%d99all2=InputFileData%Pres_BLThick with Re and AoA ] PRESSURE SIDE BOUNDARY LAYER THICKNESS METERS
!   uEdgeS [m%EdgeVelVar(1), EdgeVelAll(2) -> interpolated p%EdgeVelRat1=InputFileData%Suct_EdgeVelRat with Re and AoA]
!   uEdgeP [m%EdgeVelVar(2), EdgeVelAll(2) -> interpolated p%EdgeVelRat2=InputFileData%Pres_EdgeVelRat with Re and AoA]
!
!Inputs caluculated in AeroDyn (now set in driver input file?):
!   meanWindspeed [u%Inflow]
!   AoA [u%AoANoise]
!   [u%vRel]
!   [AeroCent_G] = u%BladeMotion(j)%Position(:,i) + u%BladeMotion(j)%TranslationDisp(:,i) (global position of the blade node) -> fixed value???
!   [RotGtoL] -> set to identity
!
!Inputs calculated
!   Ma [M or Mach] : calculated M = U  / p%SpdSound        ! MACH NUMBER
!   Re [RC] : calculated RC = U  * C/p%KinVisc  ! Reynolds number; C = chord; U=UNoise=sign( max(abs(u%Vrel(J,I)),0.1), u%Vrel(J,I) )
!
!fSep1p0_alpha (new to  BPM)
!fSpe0p7_alpha (new to  BPM)
   