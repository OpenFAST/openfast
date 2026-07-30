!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2026  National Renewable Energy Laboratory
!
!    This file is part of Ambient Wind and Array Effects model for FAST.Farm.
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
MODULE AWAE_Timings

   USE NWTC_Library

#ifdef _OPENMP
   USE OMP_LIB
#endif

   IMPLICIT NONE

   integer(IntKi), parameter, public :: AWAE_MaxTimingEntries = 32

   ! Typically we dissalow the `save` attribute for module variables, but we make an exception here since
   ! these variables are used to store timing data that is collected across multiple subroutine calls, but
   ! only outside of OMP loops.  If we were to allow OMP loops to call these subroutines, then we would
   ! need to use threadprivate variables instead of save variables.
   integer(IntKi), save      :: AWAE_TimingCount = 0
   character(128), save      :: AWAE_TimingLabel(AWAE_MaxTimingEntries)
   real(DbKi), save          :: AWAE_TimingSerial(AWAE_MaxTimingEntries) = 0.0_DbKi
   real(DbKi), save          :: AWAE_TimingPar(AWAE_MaxTimingEntries) = 0.0_DbKi
   character(64), save       :: AWAE_TimingPrefix = ''

CONTAINS

subroutine AWAE_SetTimingStage(stage)
   character(*), intent(in) :: stage
   AWAE_TimingPrefix = stage
end subroutine AWAE_SetTimingStage

subroutine AWAE_ResetTiming()
   AWAE_TimingCount   = 0
   AWAE_TimingLabel   = ''
   AWAE_TimingSerial  = 0.0_DbKi
   AWAE_TimingPar     = 0.0_DbKi
   AWAE_TimingPrefix  = ''
end subroutine AWAE_ResetTiming

subroutine AWAE_GetTimingData(numEntries, labels, serialTimes, parTimes)
   integer(IntKi), intent(out) :: numEntries
   character(*), intent(out)   :: labels(:)
   real(DbKi), intent(out)     :: serialTimes(:)
   real(DbKi), intent(out)     :: parTimes(:)
   integer(IntKi)              :: n, i

   numEntries = 0
   labels = ''
   serialTimes = 0.0_DbKi
   parTimes = 0.0_DbKi

   n = min(AWAE_TimingCount, min(size(labels), min(size(serialTimes), size(parTimes))))
   numEntries = n
   do i = 1, n
      labels(i) = AWAE_TimingLabel(i)
      serialTimes(i) = AWAE_TimingSerial(i)
      parTimes(i) = AWAE_TimingPar(i)
   end do
end subroutine AWAE_GetTimingData

subroutine AWAE_DrainTimingData(numEntries, labels, serialTimes, parTimes)
   integer(IntKi), intent(out) :: numEntries
   character(*), intent(out)   :: labels(:)
   real(DbKi), intent(out)     :: serialTimes(:)
   real(DbKi), intent(out)     :: parTimes(:)

   call AWAE_GetTimingData(numEntries, labels, serialTimes, parTimes)
   AWAE_TimingCount = 0
   AWAE_TimingLabel = ''
   AWAE_TimingSerial = 0.0_DbKi
   AWAE_TimingPar = 0.0_DbKi
end subroutine AWAE_DrainTimingData

real(DbKi) function AWAE_WallTime()
#ifdef _OPENMP
   AWAE_WallTime = omp_get_wtime()
#else
   call cpu_time(AWAE_WallTime)
#endif
end function AWAE_WallTime

subroutine AWAE_AddTiming(label, serialDt, parDt)
   character(*), intent(in) :: label
   real(DbKi),   intent(in) :: serialDt
   real(DbKi),   intent(in) :: parDt
   integer(IntKi)           :: i

   do i = 1, AWAE_TimingCount
      if (trim(AWAE_TimingLabel(i)) == trim(label)) then
         AWAE_TimingSerial(i) = AWAE_TimingSerial(i) + serialDt
         AWAE_TimingPar(i) = AWAE_TimingPar(i) + parDt
         return
      end if
   end do

   if (AWAE_TimingCount < AWAE_MaxTimingEntries) then
      AWAE_TimingCount = AWAE_TimingCount + 1
      AWAE_TimingLabel(AWAE_TimingCount) = label
      AWAE_TimingSerial(AWAE_TimingCount) = serialDt
      AWAE_TimingPar(AWAE_TimingCount) = parDt
   end if
end subroutine AWAE_AddTiming

subroutine AWAE_AddStageTiming(action, serialStart, parStart)
   character(*), intent(in) :: action
   real(DbKi),   intent(in) :: serialStart
   real(DbKi),   intent(in) :: parStart
   real(DbKi)               :: serialEnd, parEnd
   character(128)           :: label

   if (len_trim(AWAE_TimingPrefix) == 0) return

   call cpu_time(serialEnd)
   parEnd = AWAE_WallTime()
   label = trim(AWAE_TimingPrefix)//':'//trim(action)
   call AWAE_AddTiming(label, serialEnd-serialStart, parEnd-parStart)
end subroutine AWAE_AddStageTiming

END MODULE AWAE_Timings
