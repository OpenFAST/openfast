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
   implicit none   
   ! Program variables
   integer :: nt !< loop counter (for time step)

   call DvrM_Init(dat%DvrData, dat%AD, dat%IW, dat%errStat, dat%errMsg); call CheckError()
   dat%initialized=.true.

   do nt = 1, dat%DvrData%numSteps
      call DvrM_TimeStep(nt, dat%DvrData, dat%AD, dat%IW, dat%errStat, dat%errMsg); call CheckError()
   end do !nt=1,numSteps

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
   
