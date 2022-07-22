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

!****************************************************************
!  Loads from ice crushing at the specified structural frequency
!  per the ISO 19906:2010 with a sawtooth wave form
!
module iceLockInCrushingISO

   use iceCrushingISO

   implicit none

   public   

contains

   subroutine initLockInCrushingISO (iceInput, myIceParams, iceLog)

      type(iceInputType), intent(in)            :: iceInput
      type(iceFloe_ParameterType), intent(inout)     :: myIceParams
      type(iceFloe_LoggingType), intent(inout)   :: iceLog   ! structure with message and error logging variables
      type(inputParams)                         :: inParams ! specific input parameter variable list

      real(ReKi)     :: fallTime, maxLoad
      integer(IntKi) :: nL ! err, 

!  initialize the common parmeters
      call initIceCrushISO(iceInput, inParams, myIceParams, iceLog)

      call logMessage(iceLog, newLine//' Setting ice crushing loads with frequency lock-in parameteres per ISO')

      call getIceInput(iceInput, 'towerFrequency', inParams%twr%freq, iceLog, 0.01_ReKi, 10.0_ReKi)
      call logMessage(iceLog, ' Tower fundamental frequency = '//TRIM(Num2LStr(inParams%twr%freq))//' Hz')

      call getIceInput(iceInput, 'riseTime', inParams%riseTime, iceLog, 0.1_ReKi, 0.9_ReKi)
      call logMessage(iceLog, ' Sawtooth rise time fraction = '//TRIM(Num2LStr(inParams%riseTime)))
      fallTime = 1.0 - inParams%riseTime

      call getIceInput(iceInput, 'minLoadFraction', inParams%minLoadFraction, iceLog, 0.0_ReKi, 1.0_ReKi)
      call logMessage(iceLog, ' Minimum load fraction (of max) = '//TRIM(Num2LStr(inParams%minLoadFraction)))

   !  get leg load phase
      if (myIceParams%numLegs>1) then
         do nL = 1, myIceParams%numLegs
            call getIceInput(iceInput, 'loadPhase'//TRIM(Num2LStr(nL)), inParams%twr%leg(nL)%phase, iceLog, 0.0_ReKi, 360.0_ReKi)
            call logMessage(iceLog, ' Load phase for leg '//TRIM(Num2LStr(nL))//' is '                             &
                                    //TRIM(Num2LStr(inParams%twr%leg(nL)%phase))//' degrees')
            inParams%twr%leg(nL)%phase = D2R*inParams%twr%leg(nL)%phase
         enddo
         call getIceInput(iceInput, 'multiLegFactor_kn', inParams%multiLegFactor_kn, iceLog, 0.0_ReKi, 1.0_ReKi)
         call logMessage(iceLog, ' Multi leg factor = '//TRIM(Num2LStr(inParams%multiLegFactor_kn)))
      else
         inParams%twr%leg(1)%phase  = 0.0
         inParams%multiLegFactor_kn = 1.0
      endif

!  Pre-calculate sawtooth time series of loads
      maxLoad = inParams%multiLegFactor_kn*globalCrushLoadISO(inParams)
      call logMessage(iceLog, '** Global crushing load is: '//TRIM(Num2LStr(maxLoad))//' Newtons.' )

      call crushLoadTimeSeriesISO(myIceParams, inParams, iceLog, maxLoad, inParams%minLoadFraction*maxLoad,    &
                                  1.0_ReKi/inParams%twr%freq, inParams%riseTime, fallTime)

   end subroutine initLockInCrushingISO

!****************************************************************************
   function outputLockInLoadISO (myIceParams, iceLog, time)  result(iceLoads)
      type(iceFloe_ParameterType), intent(in)      :: myIceParams
      type(iceFloe_LoggingType), intent(inout) :: iceLog     ! structure with message and error logging variables
      real(DbKi), intent(in)                  :: time
      real(ReKi)                              :: iceLoads(6,myIceParams%numLegs)

      iceLoads = outputCrushLoadISO (myIceParams, iceLog, time)
   end function outputLockInLoadISO

end module iceLockInCrushingISO
