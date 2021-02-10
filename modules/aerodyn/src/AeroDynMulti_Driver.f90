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
   use AeroDynMulti_Driver_Subs   
   implicit none   
   ! Program variables
   real(DbKi)           :: time             !< Variable for storing time, in seconds
   type(DvrM_SimData)   :: DvrData              !< The data required for running the AD driver
   type(AeroDyn_Data)   :: AD                   !< AeroDyn data
   type(InflowWind_Data) :: IW                  !< AeroDyn data
   integer(IntKi)       :: nt                   !< loop counter (for time step)
   integer(IntKi)       :: j                    !< 
   integer(IntKi)       :: errStat              !< Status of error message
   character(ErrMsgLen) :: errMsg               !< Error message if ErrStat /= ErrID_None
   logical              :: ADM_Initialized
   errStat        = ErrID_None
   errMsg         = ''
   ADM_Initialized = .false.
   time           = 0.0 ! seconds

   ! --- Initialize driver
   call DvrM_Init(DvrData, ErrStat, ErrMsg); call CheckError()

   ! --- Initialize aerodyn 
   call Init_AeroDyn(DvrData, AD, DvrData%dT, errStat, errMsg); call CheckError()

   ! --- Initialize Inflow Wind 
   call Init_InflowWind(DvrData, IW, AD, DvrData%dt, errStat, errMsg); call CheckError()

   ! --- Initial AD inputs
   AD%InputTime = -999
   DO j = 1-numInp, 0
      call Set_AD_Inputs(j,DvrData,AD,IW,errStat,errMsg); call CheckError()
   END DO              


   ADM_Initialized = .true.
   call DvrM_InitializeOutputFile(DvrData%OutFileData, errStat, errMsg)
   call CheckError()

   do nt = 1, DvrData%numSteps
      !...............................
      ! set AD inputs for nt (and keep values at nt-1 as well)
      !...............................
      ! u(1) is at nt+1, u(2) is at nt
      call Set_AD_Inputs(nt,DvrData,AD,IW,errStat,errMsg); call CheckError()
      time = AD%InputTime(2)
      if (mod(nd,10)==0) then
         print*,'time',time
      endif
      ! Calculate outputs at nt - 1
      call AD_CalcOutput( time, AD%u(2), AD%p, AD%x, AD%xd, AD%z, AD%OtherState, AD%y, AD%m, errStat, errMsg ); call CheckError()

      call Dvr_WriteOutputLine(time, DvrData%OutFileData, AD%y%WriteOutput, IW%y%WriteOutput, errStat, errMsg); call CheckError()
      ! Get state variables at next step: INPUT at step nt - 1, OUTPUT at step nt

      call AD_UpdateStates( time, nt-1, AD%u, AD%InputTime, AD%p, AD%x, AD%xd, AD%z, AD%OtherState, AD%m, errStat, errMsg )
      call CheckError()
   end do !nt=1,numSteps
   call DvrM_End()
contains
!................................   
   subroutine CheckError()
      if (ErrStat /= ErrID_None) then
         call WrScr(TRIM(ErrMsg))
         if (ErrStat >= AbortErrLev) then
            call DvrM_End()
         end if
      end if
   end subroutine CheckError
!................................   
   subroutine DvrM_End()
      character(ErrMsgLen)                          :: errMsg2                 ! temporary Error message if ErrStat /= ErrID_None
      integer(IntKi)                                :: errStat2                ! temporary Error status of the operation
      character(*), parameter                       :: RoutineName = 'DvrM_End'
      ! Close the output file
      if (DvrData%OutFileData%unOutFile > 0) close(DvrData%OutFileData%unOutFile)
            
      if ( ADM_Initialized ) then
         call AD_End( AD%u(1), AD%p, AD%x, AD%xd, AD%z, AD%OtherState, AD%y, AD%m, errStat2, errMsg2 )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
      end if
           
      call ADM_Dvr_DestroyDvrM_SimData( DvrData, ErrStat2, ErrMsg2 )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )

      call ADM_Dvr_DestroyAeroDyn_Data( AD, ErrStat2, ErrMsg2 )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
               
      if (ErrStat >= AbortErrLev) then      
         CALL ProgAbort( 'AeroDyn Driver encountered simulation error level: '&
             //TRIM(GetErrStr(ErrStat)), TrapErrors=.FALSE., TimeWait=3._ReKi )  ! wait 3 seconds (in case they double-clicked and got an error)
      else
         call NormStop()
      end if
   end subroutine DvrM_End
!................................   
end program AeroDynMulti_Driver
   
