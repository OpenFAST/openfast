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

!  Module for calculations of loads induced by ice crushing against vertical surfaces
!  As suggested in IEC 61400-3, Ed 1 Annex E
!  This is effectively a lock in condition where the ice loads are 
!  sinusoidal at the specified structural frequency
!
module iceCrushingIEC

   use IceFloeBase

   implicit none

   public
   
contains

!  calculate the time series of ice loads for continuous crushing 
!****************************************************************
   subroutine initCrushingIEC (iceInput, myIceParams, iceLog)

      type(iceInputType), intent(in)            :: iceInput    ! Parameter read from input file for initialization
      type(iceFloe_ParameterType), intent(inout)     :: myIceParams ! saved parameters
      type(iceFloe_LoggingType), intent(inout)   :: iceLog      ! structure with message and error logging variables
      type(inputParams)                         :: inParams    ! specific input parameter variable list

!  local variables
      real(ReKi)     :: maxLoad     ! Maximum static ice load, Newtons
      integer(IntKi) :: nL !err, 

!  initialize the common parmeters
      call initIceFloe(iceInput, inParams, myIceParams, iceLog)

!  initialize the parmeters for flexural ice failure via IEC
      call logMessage(iceLog, newLine//' Setting parameters for crushing on vertical surfaces by IEC/Korzhavin method')

      call getIceInput(iceInput, 'refIceStrength', inParams%refIceStrength, iceLog, 0.5E6_ReKi, 5.0E6_ReKi)
      call logMessage(iceLog, ' Reference ice strength = '//TRIM(Num2LStr(inParams%refIceStrength))//' Pascals')

      call getIceInput(iceInput, 'shapeFactor_k1', inParams%shapeFactor_k1, iceLog, 0.1_ReKi, 1.0_ReKi)
      call logMessage(iceLog, ' Static load shape factor = '//TRIM(Num2LStr(inParams%shapeFactor_k1)))

      call getIceInput(iceInput, 'contactFactor_k2', inParams%contactFactor_k2, iceLog, 0.1_ReKi, 1.0_ReKi)
      call logMessage(iceLog, ' Static load contact factor = '//TRIM(Num2LStr(inParams%contactFactor_k2)))

      call getIceInput(iceInput, 'towerFrequency', inParams%twr%freq, iceLog, 0.01_ReKi, 10.0_ReKi)
      call logMessage(iceLog, ' Tower fundamental frequency = '//TRIM(Num2LStr(inParams%twr%freq))//' Hz')

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

      maxLoad = inParams%multiLegFactor_kn*globalCrushLoadIEC(inParams)
      call logMessage(iceLog, '** Global crushing load is: '//TRIM(Num2LStr(maxLoad))//' Newtons.' )

!  Precalculates a time series of loads based on:
!  Sinusoidal per IEC
      call IECLoadTimeSeries(myIceParams, inParams, iceLog, maxLoad, inParams%twr%freq)

   end subroutine initCrushingIEC
   
!****************************************************************************
! Per IEC 61400-3 Design Requirements for Offshore Wind Turbines Ed 1, 2009
! via the Korshavin equation
   function globalCrushLoadIEC(inParams) result(load)
      type(inputParams), intent(in) :: inParams ! specific input parameter variable list
      real(ReKi) :: load
      real(ReKi) :: aspectRatio_k3  ! Per IEC
      
      aspectRatio_k3 = sqrt(1.0+5.0*inParams%iceThickness/inParams%twr%diam)
      load = inParams%shapeFactor_k1*inParams%contactFactor_k2*aspectRatio_k3
      load = load*inParams%iceThickness*inParams%twr%diam*inParams%refIceStrength
      
   end function globalCrushLoadIEC

!****************************************************************************
!  Continuous crushing uses the standard inerpolation routine
!  of the precalculated time series
   function outputCrushLoadIEC (myIceParams, iceLog, time)  result(iceLoads)
      type(iceFloe_ParameterType), intent(in)        :: myIceParams
      type(iceFloe_LoggingType), intent(inout)   :: iceLog   ! structure with message and error logging variables
      real(DbKi), intent(in)                    :: time
      real(ReKi)                                :: iceLoads(6,myIceParams%numLegs)

      iceLoads = outputIceLoads (myIceParams, iceLog, time)
   end function outputCrushLoadIEC

end module iceCrushingIEC