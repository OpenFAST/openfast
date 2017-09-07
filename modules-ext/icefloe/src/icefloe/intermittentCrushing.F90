!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2014  DNV KEMA Renewables, Inc.
!
!    This file is part of the IceFloe suite of subroutines
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
!************************************************************************

! Module to initialize and calculate a time series of intermittent
! loads from ice crushing

module iceIntermittentCrushing
       
   use iceCrushingISO

   implicit none

   public
   
contains

!========================================================================================
!  calculate the time series of ice loads for intermittent crushing 
   subroutine initInterCrushing (iceInput, myIceParams, iceLog)

      type(iceInputType), intent(in)            :: iceInput    ! Parameters from input file for initialization
      type(iceFloe_ParameterType), intent(inout)     :: myIceParams ! saved parameters
      type(iceFloe_LoggingType), intent(inout)   :: iceLog      ! structure with message and error logging variables

      type(inputParams)    :: inParams    ! specific input parameter variable list
      real(ReKi)           :: maxLoad  ! Global maximum crushing load
      integer(IntKi)       :: nL !err, 
      
!  initialize the common parmeters for flexural ice failure
      call initIceCrushISO(iceInput, inParams, myIceParams, iceLog)

      call logMessage(iceLog, newLine//' Setting parameters for intermittent crushing loads' )

      call getIceInput(iceInput, 'interPeriod', inParams%interPeriod, iceLog, 1.0_ReKi)
      call logMessage(iceLog, ' Intermittent period = '//TRIM(Num2LStr(inParams%interPeriod))//' seconds')

      call getIceInput(iceInput, 'riseTime', inParams%riseTime, iceLog, 0.1_ReKi, 0.9_ReKi)
      call logMessage(iceLog, ' Saw tooth rise time fraction = '//TRIM(Num2LStr(inParams%riseTime)))

      call getIceInput(iceInput, 'fallTime', inParams%fallTime, iceLog, 0.1_ReKi, 1.0_ReKi-inParams%riseTime)
      call logMessage(iceLog, ' Saw tooth fall time fraction = '//TRIM(Num2LStr(inParams%fallTime)))

   !  get leg load phase
      if (myIceParams%numLegs>1) then
         do nL = 1, myIceParams%numLegs
            call getIceInput(iceInput, 'loadPhase'//TRIM(Num2LStr(nL)), inParams%twr%leg(nL)%phase, iceLog, 0.0_ReKi, 360.0_ReKi)
            call logMessage(iceLog, ' Load phase for leg '//TRIM(Num2LStr(nL))//' is '                             &
                                    //TRIM(Num2LStr(inParams%twr%leg(nL)%phase))//' degrees')
            inParams%twr%leg(nL)%phase = D2R*inParams%twr%leg(nL)%phase
         enddo
      else
         inParams%twr%leg(1)%phase = 0.0
      endif

!  The method precalculates a time series of loads based
!  on the input parameters and a random distributions of components
!  of a spectral model
      maxLoad = globalCrushLoadISO(inParams)
      call logMessage(iceLog, '** Global crushing load is: '//TRIM(Num2LStr(maxLoad))//' Newtons.' )

      call crushLoadTimeSeriesISO(myIceParams, inParams, iceLog, maxLoad, 0.0_ReKi, inParams%interPeriod,  &
                                  inParams%riseTime, inParams%fallTime)

   end subroutine initInterCrushing

!****************************************************************************
!  Continuous crushing uses the standard inerpolation routine
!  of the precalculated time series
   function outputInterCrushLoad (myIceParams, iceLog, time)  result(iceLoads)
      type(iceFloe_ParameterType), intent(in)      :: myIceParams
      type(iceFloe_LoggingType), intent(inout) :: iceLog   ! structure with message and error logging variables
      real(DbKi), intent(in)                  :: time
      real(ReKi)                              :: iceLoads(6,myIceParams%numLegs)

      iceLoads = outputCrushLoadISO (myIceParams, iceLog, time)
   end function outputInterCrushLoad

end module iceIntermittentCrushing