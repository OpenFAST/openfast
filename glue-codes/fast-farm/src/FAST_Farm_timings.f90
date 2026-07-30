!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2026  National Renewable Energy Laboratory
!
!    This file is part of FAST.Farm.
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
MODULE FAST_Farm_Timings

   USE NWTC_Library
   USE AWAE_Timings, only: AWAE_ResetTiming, AWAE_DrainTimingData, AWAE_MaxTimingEntries

#ifdef _OPENMP
   USE OMP_LIB
#endif

   IMPLICIT NONE

   integer(IntKi), parameter :: FFTiming_MaxEntries = 256
   integer(IntKi), save      :: FFTiming_Count = 0
   character(128), save      :: FFTiming_Label(FFTiming_MaxEntries)
   real(DbKi), save          :: FFTiming_Serial(FFTiming_MaxEntries) = 0.0_DbKi
   real(DbKi), save          :: FFTiming_Par(FFTiming_MaxEntries) = 0.0_DbKi

CONTAINS

subroutine FarmTiming_Reset()
   call AWAE_ResetTiming()
   FFTiming_Count = 0
   FFTiming_Label = ''
   FFTiming_Serial = 0.0_DbKi
   FFTiming_Par = 0.0_DbKi
end subroutine FarmTiming_Reset

real(DbKi) function FarmTiming_WallTime()
#ifdef _OPENMP
   FarmTiming_WallTime = omp_get_wtime()
#else
   call cpu_time(FarmTiming_WallTime)
#endif
end function FarmTiming_WallTime

subroutine FarmTiming_Add(label, serialDt, parDt)
   character(*), intent(in) :: label
   real(DbKi),   intent(in) :: serialDt
   real(DbKi),   intent(in) :: parDt
   integer(IntKi)           :: i

   do i = 1, FFTiming_Count
      if (trim(FFTiming_Label(i)) == trim(label)) then
         FFTiming_Serial(i) = FFTiming_Serial(i) + serialDt
         FFTiming_Par(i) = FFTiming_Par(i) + parDt
         return
      end if
   end do

   if (FFTiming_Count < FFTiming_MaxEntries) then
      FFTiming_Count = FFTiming_Count + 1
      FFTiming_Label(FFTiming_Count) = label
      FFTiming_Serial(FFTiming_Count) = serialDt
      FFTiming_Par(FFTiming_Count) = parDt
   end if
end subroutine FarmTiming_Add

subroutine FarmTiming_AddElapsed(label, serialStart, parStart)
   character(*), intent(in) :: label
   real(DbKi),   intent(in) :: serialStart
   real(DbKi),   intent(in) :: parStart
   real(DbKi)               :: serialEnd, parEnd

   call cpu_time(serialEnd)
   parEnd = FarmTiming_WallTime()
   call FarmTiming_Add(label, serialEnd-serialStart, parEnd-parStart)
end subroutine FarmTiming_AddElapsed

subroutine FarmTiming_DrainAWAEDetails()
   integer(IntKi)                                 :: n, i
   character(128)                                 :: labels(AWAE_MaxTimingEntries)
   real(DbKi)                                     :: serial(AWAE_MaxTimingEntries)
   real(DbKi)                                     :: par(AWAE_MaxTimingEntries)

   call AWAE_DrainTimingData(n, labels, serial, par)
   do i = 1, n
      call FarmTiming_Add(trim(labels(i)), serial(i), par(i))
   end do
end subroutine FarmTiming_DrainAWAEDetails

logical function FarmTiming_IsStageChild(label, stage)
   character(*), intent(in) :: label
   character(*), intent(in) :: stage
   character(24)            :: prefix
   integer(IntKi)           :: nColons

   prefix = trim(stage)//':'
   if (index(trim(label), trim(prefix)) /= 1) then
      FarmTiming_IsStageChild = .false.
      return
   end if

   nColons = FarmTiming_CountColons(label)
   FarmTiming_IsStageChild = (nColons == 2)
end function FarmTiming_IsStageChild

integer(IntKi) function FarmTiming_CountColons(label)
   character(*), intent(in) :: label
   integer(IntKi)           :: i

   FarmTiming_CountColons = 0
   do i = 1, len_trim(label)
      if (label(i:i) == ':') FarmTiming_CountColons = FarmTiming_CountColons + 1
   end do
end function FarmTiming_CountColons

logical function FarmTiming_IsParentChild(label, parent)
   character(*), intent(in) :: label
   character(*), intent(in) :: parent
   character(128)           :: prefix
   integer(IntKi)           :: parentColons, labelColons

   prefix = trim(parent)//':'
   if (index(trim(label), trim(prefix)) /= 1) then
      FarmTiming_IsParentChild = .false.
      return
   end if

   parentColons = FarmTiming_CountColons(parent)
   labelColons = FarmTiming_CountColons(label)
   FarmTiming_IsParentChild = (labelColons == parentColons + 1)
end function FarmTiming_IsParentChild

subroutine FarmTiming_AddParentOthers(parent)
   character(*), intent(in) :: parent
   integer(IntKi)           :: i, iParent
   real(DbKi)               :: serialChildren, parChildren
   real(DbKi)               :: serialOthers, parOthers
   character(128)           :: othersLabel

   iParent = 0
   do i = 1, FFTiming_Count
      if (trim(FFTiming_Label(i)) == trim(parent)) then
         iParent = i
         exit
      end if
   end do
   if (iParent == 0) return

   serialChildren = 0.0_DbKi
   parChildren = 0.0_DbKi
   do i = 1, FFTiming_Count
      if (FarmTiming_IsParentChild(FFTiming_Label(i), parent)) then
         serialChildren = serialChildren + FFTiming_Serial(i)
         parChildren = parChildren + FFTiming_Par(i)
      end if
   end do

   serialOthers = FFTiming_Serial(iParent) - serialChildren
   parOthers = FFTiming_Par(iParent) - parChildren
   if (abs(serialOthers) < 1.0e-9_DbKi) serialOthers = 0.0_DbKi
   if (abs(parOthers) < 1.0e-9_DbKi) parOthers = 0.0_DbKi

   othersLabel = trim(parent)//':Others'
   call FarmTiming_Add(trim(othersLabel), serialOthers, parOthers)
end subroutine FarmTiming_AddParentOthers

subroutine FarmTiming_AddStageOthers(stage)
   character(*), intent(in) :: stage
   integer(IntKi)           :: i
   integer(IntKi)           :: iStage
   real(DbKi)               :: serialChildren, parChildren
   real(DbKi)               :: serialOthers, parOthers
   character(32)            :: othersLabel

   iStage = 0
   do i = 1, FFTiming_Count
      if (trim(FFTiming_Label(i)) == trim(stage)) then
         iStage = i
         exit
      end if
   end do
   if (iStage == 0) return

   serialChildren = 0.0_DbKi
   parChildren = 0.0_DbKi
   do i = 1, FFTiming_Count
      if (FarmTiming_IsStageChild(FFTiming_Label(i), stage)) then
         serialChildren = serialChildren + FFTiming_Serial(i)
         parChildren = parChildren + FFTiming_Par(i)
      end if
   end do

   serialOthers = FFTiming_Serial(iStage) - serialChildren
   parOthers = FFTiming_Par(iStage) - parChildren
   if (abs(serialOthers) < 1.0e-9_DbKi) serialOthers = 0.0_DbKi
   if (abs(parOthers) < 1.0e-9_DbKi) parOthers = 0.0_DbKi

   othersLabel = trim(stage)//':Others'
   call FarmTiming_Add(trim(othersLabel), serialOthers, parOthers)
end subroutine FarmTiming_AddStageOthers

subroutine FarmTiming_PrintLiveStart(stage, label, mode)
   character(*), intent(in) :: stage
   character(*), intent(in) :: label
   character(*), intent(in) :: mode
   character(384)           :: line

   line = 'FFTiming|LIVE|Event=Start|Stage='//trim(stage)//'|Label='//trim(label)//'|Mode='//trim(mode)
   call WrScr(trim(line))
end subroutine FarmTiming_PrintLiveStart

subroutine FarmTiming_PrintLiveDone(stage, label, mode, dt)
   character(*), intent(in) :: stage
   character(*), intent(in) :: label
   character(*), intent(in) :: mode
   real(DbKi),   intent(in) :: dt
   character(384)           :: line

   line = 'FFTiming|LIVE|Event=Done|Stage='//trim(stage)//'|Label='//trim(label)//'|Mode='//trim(mode)//'|Seconds='//trim(num2lstr(dt))
   call WrScr(trim(line))
end subroutine FarmTiming_PrintLiveDone

subroutine FarmTiming_PrintLiveThread(stage, label, itemId, threadId, dt)
   character(*), intent(in) :: stage
   character(*), intent(in) :: label
   integer(IntKi), intent(in) :: itemId
   integer(IntKi), intent(in) :: threadId
   real(DbKi),   intent(in) :: dt
   character(384)           :: line

   line = 'FFTiming|LIVE|Event=ThreadDone|Stage='//trim(stage)//'|Label='//trim(label)//'|Item='//trim(num2lstr(itemId))// &
          '|Thread='//trim(num2lstr(threadId))//'|Seconds='//trim(num2lstr(dt))
   call WrScr(trim(line))
end subroutine FarmTiming_PrintLiveThread

subroutine FarmTiming_PrintSummary()
   integer(IntKi)   :: i
   character(384)   :: line

   call FarmTiming_DrainAWAEDetails()
   call FarmTiming_AddParentOthers('ICO:AWAE:CalcOutput1')
   call FarmTiming_AddParentOthers('ICO:AWAE:CalcOutput2')
   call FarmTiming_AddParentOthers('CO:AWAE:UpdateStates')
   call FarmTiming_AddParentOthers('CO:AWAE:CalcOutput')
   call FarmTiming_AddStageOthers('INIT')
   call FarmTiming_AddStageOthers('ICO')
   call FarmTiming_AddStageOthers('US')
   call FarmTiming_AddStageOthers('CO')
   call FarmTiming_AddStageOthers('END')

   call WrScr('')
   call WrScr('FFTiming|SUMMARY|Begin')
   call WrScr('FFTiming|SUMMARY|Columns=Label|CPU|Wall')
   do i = 1, FFTiming_Count
      line = 'FFTiming|SUMMARY|Label='//trim(FFTiming_Label(i))//'|CPU='//trim(num2lstr(FFTiming_Serial(i)))//'|Wall='//trim(num2lstr(FFTiming_Par(i)))
      call WrScr(trim(line))
   end do
   call WrScr('FFTiming|SUMMARY|End')
   call WrScr('')
end subroutine FarmTiming_PrintSummary

END MODULE FAST_Farm_Timings
