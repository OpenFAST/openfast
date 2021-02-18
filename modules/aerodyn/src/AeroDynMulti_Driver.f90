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
!**********************************************************************************************************************************
program AeroDynMulti_Driver
   use AeroDynMulti_Driver_Subs, only: dat, DvrM_Init, DvrM_TimeStep, DvrM_CleanUp
   use NWTC_IO
   use NWTC_Num, only: RunTimes, SimStatus, SimStatus_FirstTime
   implicit none   
   ! Program variables
   REAL(ReKi)                       :: PrevClockTime ! Clock time at start of simulation in seconds [(s)]
   REAL(ReKi)                       :: UsrTime1      ! User CPU time for simulation initialization [(s)]
   REAL(ReKi)                       :: UsrTime2      ! User CPU time for simulation (without intialization) [(s)]
   INTEGER(IntKi) , DIMENSION(1:8)  :: StrtTime      ! Start time of simulation (including intialization) [-]
   INTEGER(IntKi) , DIMENSION(1:8)  :: SimStrtTime   ! Start time of simulation (after initialization) [-]
   REAL(DbKi)                       :: t_global         ! global-loop time marker
   REAL(DbKi)                       :: t_final         ! global-loop time marker
   REAL(DbKi)                       :: TiLstPrn      ! The simulation time of the last print (to file) [(s)]
   integer :: nt !< loop counter (for time step)
   CALL DATE_AND_TIME ( Values=StrtTime )                 ! Let's time the whole simulation
   CALL CPU_TIME ( UsrTime1 )                             ! Initial time (this zeros the start time when used as a MATLAB function)
   UsrTime1 = MAX( 0.0_ReKi, UsrTime1 )                   ! CPU_TIME: If a meaningful time cannot be returned, a processor-dependent negative value is returned

   call DvrM_Init(dat%DvrData, dat%AD, dat%IW, dat%errStat, dat%errMsg); call CheckError()
   dat%initialized=.true.

   ! Init of time estimator
   t_global=0.0_DbKi
   t_final=dat%DvrData%numSteps*dat%DvrData%dt
   call SimStatus_FirstTime( TiLstPrn, PrevClockTime, SimStrtTime, UsrTime2, t_global, t_final )

   do nt = 1, dat%DvrData%numSteps

      call DvrM_TimeStep(nt, dat%DvrData, dat%AD, dat%IW, dat%errStat, dat%errMsg); call CheckError()

      ! Time update to screen
      t_global=nt*dat%DvrData%dt
      if (mod( nt + 1, 10 )==0) call SimStatus(TiLstPrn, PrevClockTime, t_global, t_final)

   end do !nt=1,numSteps

   ! display runtime to screen
   call RunTimes(StrtTime, UsrTime1, SimStrtTime, UsrTime2, t_global)
   
   call DvrM_End()
contains
!................................   
   subroutine CheckError()
      if (dat%ErrStat /= ErrID_None) then
         call WrScr(TRIM(dat%errMsg))
         if (dat%errStat >= AbortErrLev) then
            call DvrM_End()
         end if
      end if
   end subroutine CheckError
!................................   
   subroutine DvrM_End()

      call DvrM_CleanUp(dat%DvrData, dat%AD, dat%IW, dat%initialized, dat%errStat, dat%errMsg)

      if (dat%errStat >= AbortErrLev) then      
         CALL ProgAbort( 'AeroDyn Driver encountered simulation error level: '&
             //TRIM(GetErrStr(dat%errStat)), TrapErrors=.FALSE., TimeWait=3._ReKi )  ! wait 3 seconds (in case they double-clicked and got an error)
      else
         call NormStop()
      end if
   end subroutine DvrM_End
!................................   
end program AeroDynMulti_Driver
   
