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
 
!  Module for common calculations of loads induced by ice crushing against vertical surfaces
!  Using the method in ISO 19906:2010, section A.8.2.4.3.3
!
module iceCrushingISO

   use iceFloeBase

   implicit none

   public
   
contains

!  calculate the time series of ice loads using the ISO method
   subroutine initIceCrushISO (iceInput, inParams, myIceParams, iceLog)

      type(iceInputType), intent(in)            :: iceInput
      type(iceFloe_ParameterType), intent(inout)     :: myIceParams
      type(iceFloe_LoggingType), intent(inout)   :: iceLog   ! structure with message and error logging variables

      type(inputParams), intent(out)            :: inParams ! specific input parameter variable list
      
!  initialize the common parmeters
      call initIceFloe(iceInput, inParams, myIceParams, iceLog)

      call logMessage(iceLog, newLine//' Setting common ice crushing input parameters ')

      call getIceInput(iceInput, 'refIceStrength', inParams%refIceStrength, iceLog, 0.1E6_ReKi, 50.0E6_ReKi)
      call logMessage(iceLog, ' Reference ice strength = '//TRIM(Num2LStr(inParams%refIceStrength))//' Pascals')

      call getIceInput(iceInput, 'refIceThick', inParams%refIceThick, iceLog, 1.0_ReKi, 1.0_ReKi)
      call logMessage(iceLog, ' Reference ice thickness = '//TRIM(Num2LStr(inParams%refIceThick))//' meters ')

      call getIceInput(iceInput, 'staticExponent', inParams%staticExponent, iceLog, -0.16_ReKi, -0.16_ReKi)
      call logMessage(iceLog, ' Exponent for static load = '//TRIM(Num2LStr(inParams%staticExponent)))

   end subroutine initIceCrushISO


!==================================================================
! Based on ISO 19906, 2010
   function globalCrushLoadISO(inParams) result(load)
      type(inputParams), intent(in) :: inParams ! specific input parameter variable list
      real(ReKi) :: load   ! Newtons
      real(ReKi) :: exp_n  ! per ISO
      
      exp_n = min(-0.3_ReKi, 0.2_ReKi*inParams%iceThickness-0.5_ReKi)
      
      load = inParams%twr%diam*inParams%iceThickness*inParams%refIceStrength
      load = load*(inParams%iceThickness/inParams%refIceThick)**exp_n
      load = load*(inParams%twr%diam/inParams%iceThickness)**inParams%staticExponent

   end function globalCrushLoadISO

!============================================================================================
!  Calculate a sawtooth waveform with fixed period, rise, and fall times
!  Used by intermittent crushing and lock-in models
   subroutine crushLoadTimeSeriesISO (myIceParams, inParams, iceLog, peakLoad, minLoad, period, riseTime, fallTime)

      type(inputParams), intent(in)             :: inParams ! specific input parameter variable list
      type(iceFloe_ParameterType), intent(inout)     :: myIceParams
      type(iceFloe_LoggingType), intent(inout)   :: iceLog   ! structure with message and error logging variables
      real(ReKi), intent(in)                    :: peakLoad    ! Global maximum load, Newtons
      real(ReKi), intent(in)                    :: minLoad     ! Lower limit of dynamic load, Newtons
      real(ReKi), intent(in)                    :: period      ! Loading/unloading cycle period, seconds
      real(ReKi), intent(in)                    :: riseTime    ! fraction of period for which load is rising
      real(ReKi), intent(in)                    :: fallTime    ! fraction of period for which load is falling

!      integer(IntKi) :: err   ! error status from memory allocation
      integer(IntKi) :: n, ns, nL, nSteps, nRiseSteps, nFallSteps, minLoadSteps   ! various counters
      integer(IntKi) :: nStart      ! where to start in the phase of rising/falling/min load period
      
!  Number of steps in time series of loads
      nSteps = size(myIceParams%loadSeries,1)

!  Loop on all of the legs
      do nL = 1, myIceParams%numLegs
         n = 1
         nRiseSteps = nint(riseTime * period / myIceParams%dt)
         nFallSteps = nint( fallTime * period / myIceParams%dt)
         minLoadSteps = nint( period / myIceParams%dt) - nRiseSteps - nFallSteps
         if (nRiseSteps == 0 .or. nFallSteps == 0) then
            call iceErrorHndlr (iceLog, ErrID_Fatal, 'Time step too large or period too short in CrushxLoadTimeSeriesISO', 1)
            return
         endif

         nStart = nint(period*inParams%twr%leg(nL)%phase/2.0/pi/myIceParams%dt)  ! where to start in period

         seriesLoop: do while (n < nSteps+1)      
   !  begin to increase the load from the minimum to the peak
   !  based on phase we may not start on a rise i.e. nStart > nRiseSteps
            do ns = nStart, nRiseSteps
               myIceParams%loadSeries(n,nL) = minLoad + (peakLoad-minLoad)*float(ns)/float(nRiseSteps)
               n = n + 1
               if (n > nSteps) exit seriesLoop
            enddo

   !  If, dependent on phase, we started during the rise then we don't deal with it any more
            if (nStart <= nRiseSteps) then
               nStart = 1
            else
               nStart = nStart - nRiseSteps
            endif
            
   !  begin to decrease the load from the peak down to the minimum
            do ns = nStart, nFallSteps
               myIceParams%loadSeries(n,nL) = peakLoad - (peakLoad-minLoad)*float(ns-1)/float(nFallSteps)
               n = n + 1
               if (n > nSteps) exit seriesLoop
            enddo

   !  If, dependent on phase, we started during the fall then we don't deal with it any more
            if (nStart <= nRiseSteps+nFallSteps) then
               nStart = 1
            else
               nStart = nStart - nRiseSteps - nFallSteps
            endif

   !  Fill the remaining time steps in the period with the minimum load
            do ns = nStart, minLoadSteps
               myIceParams%loadSeries(n,nL) = minLoad
               n = n + 1
               if (n > nSteps) exit seriesLoop
            enddo
            nStart = 1

         enddo seriesLoop

   !  Apply shelter factor to individual legs
         myIceParams%loadSeries(:,nL) = myIceParams%ks(nL)*myIceParams%loadSeries(:,nL)

      enddo  ! end of leg loop

   end subroutine crushLoadTimeSeriesISO

!****************************************************************************
   function outputCrushLoadISO (myIceParams, iceLog, time)  result(iceLoads)
      type(iceFloe_ParameterType), intent(in)      :: myIceParams
      type(iceFloe_LoggingType), intent(inout) :: iceLog   ! structure with message and error logging variables
      real(DbKi), intent(in)                  :: time
      real(ReKi)                              :: iceLoads(6,myIceParams%numLegs)

      iceLoads = outputIceLoads (myIceParams, iceLog, time)
   end function outputCrushLoadISO

end module iceCrushingISO
