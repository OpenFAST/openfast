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
!  Using a coupled model from Maattanen, M. (1998) "Numerical Model for Ice Induced 
!     vibration load lock-in and Synchronization".  Proceedings of the 14th Int. Symp.
!     on Ice, Potsdam, NY USA, vol 2, pp 923-930.

module iceCpldCrushing

   use iceFloeBase

   implicit none

   public

!  local maximum and minima of stress rate to strength curve (fixed since curve is hard coded)
   real(ReKi), parameter   :: stressRateAtMax = 0.2914592
   real(ReKi), parameter   :: stressRateAtMin = 1.3287178
   real(ReKi), parameter   :: strengthAtMin = 2.00 + stressRateAtMin*(7.80 + &
                              stressRateAtMin*(-18.57 + stressRateAtMin*(13.00 - stressRateAtMin*2.91)))

contains

!  calculate the time series of ice loads using the ISO method
   subroutine initCpldCrushing (iceInput, myIceParams, iceLog)

      type(iceInputType), intent(in)               :: iceInput    ! Input parameters read from file
      type(iceFloe_ParameterType), intent(inout)   :: myIceParams
      type(iceFloe_LoggingType), intent(inout)     :: iceLog      ! structure with message and error logging variables
      type(inputParams)                            :: inParams    ! specific input parameter variable list

      real(ReKi)     :: defaultWidth
!      integer(IntKi) :: err

!  initialize the common parmeters
      call initIceFloe(iceInput, inParams, myIceParams, iceLog)

      call logMessage(iceLog, newLine//' Setting coupled ice crushing input parameters ')

      call getIceInput(iceInput, 'refIceStrength', inParams%refIceStrength, iceLog, 0.5E6_ReKi, 5.0E6_ReKi)
      call logMessage(iceLog, ' Reference ice strength = '//TRIM(Num2LStr(inParams%refIceStrength))//' Pascals')

      call getIceInput(iceInput, 'minStrength', inParams%minStrength, iceLog, 0.0_ReKi, 1.0E9_ReKi)
      call logMessage(iceLog, ' Minimum ice strength from stress rate polynomial = '//        &
                              TRIM(Num2LStr(inParams%minStrength))//' Pascals')

      call getIceInput(iceInput, 'minStrengthNegVel', inParams%minStrengthNegVel, iceLog, 0.0_ReKi, 1.0E9_ReKi)
      call logMessage(iceLog, ' Minimum negative velocity ice strength = '//        &
                              TRIM(Num2LStr(inParams%minStrengthNegVel))//' Pascals')

!  Assign some input parameters to saved parameters
      myIceParams%iceVel            = inParams%iceVelocity
      myIceParams%minStrength       = inParams%minStrength
      myIceParams%minStrengthNegVel = inParams%minStrengthNegVel

!  Precalc some things.
      defaultWidth = inParams%twr%diam
      if (inParams%twr%diam > 2*inParams%iceThickness) then
         call logMessage(iceLog, '*** Tower diameter > 2*ice thickness; width for cpld crushing strength '// &
                                 'set to latter value. ***')
         defaultWidth = 2*inParams%iceThickness
      endif
      myIceParams%defaultArea = sqrt(1.0/(defaultWidth*inParams%iceThickness))

!  This is for force from strength, uses actual tower diameter
      myIceParams%crushArea   = inParams%twr%diam*inParams%iceThickness
!  Scale to MPa
      myIceParams%coeffStressRate = 8.0E-6*inParams%refIceStrength/PI/defaultWidth

   end subroutine initCpldCrushing

!****************************************************************************
   function outputCpldCrushLoad (myIceParams, iceLog, inVels, time)  result(iceLoads)
      type(iceFloe_ParameterType), intent(in)      :: myIceParams
      type(iceFloe_LoggingType), intent(inout) :: iceLog   ! structure with message and error logging variables
      real(ReKi), intent(in)                  :: inVels(2,myIceParams%numLegs)
      real(DbKi), intent(in)                  :: time
      real(ReKi)                              :: iceLoads(6,myIceParams%numLegs)

      real(ReKi)     :: twrVelAngle    !  Angle of tower motion vector + ro x-axis towards y-axis
      real(ReKi)     :: velProjected   !  Velocity component projected onto direction of ice movement
      real(ReKi)     :: stressRate     !  MPa/sec
      real(ReKi)     :: strength       !  Instantaneous ice strength, Pascals
      integer(IntKi) :: nL             !  loop index for number of legs

      do nL = 1, myIceParams%numLegs
   !  get tower velocity in direction of ice floe flow
         twrVelAngle  = atan2(inVels(2,nL),inVels(1,nL)) 
         velProjected = sqrt(inVels(1,nL)**2 + inVels(2,nL)**2)*cos(twrVelAngle - myIceParams%iceDirection)
         stressRate = (myIceParams%iceVel - velProjected)*myIceParams%coeffStressRate

   !  This polynomial is specific to the Maattanen method and is thus hard coded
   !  Only applicable from a stress rate of slightly below zero to +1.3287MPa/sec
   !  In this range the failure mode is ductile crushing
         strength = 2.00 + stressRate*(7.80 + stressRate*(-18.57 + stressRate*(13.00 - stressRate*2.91)))
         strength = 1.0E6*strength*myIceParams%defaultArea
         if (stressRate <= stressRateAtMax) then
   !  The ice strength is limited to a user specified minimum for negative stress rate (tower moving away from ice forming a gap)
            strength = max(strength, myIceParams%minStrengthNegVel)
         else
   !  For high positive relative velocities (high stress rate) the ice failure mode is brittle and constant
   !  The ice strength is limited to the minimum point of the strength vs stress rate curve
            if (stressRate > stressRateAtMin) then
               strength = 1.0E6*strengthAtMin*myIceParams%defaultArea
            endif
   !  Unless the user has specified something even higher for the minimum strength
            strength = max(strength, myIceParams%minStrength)
         endif         

   !  apply the load in the proper direction via x,y components
         iceLoads(:,nL) = myIceParams%ks(nL)*iceLoadDirection(strength*myIceParams%crushArea, myIceParams)
      enddo

   !  If the loads from multiple legs are applied as an "effective" load on one leg at the
   !  support center, then compute that effective load vector and put into the first load vector
   !  Note that only Fx, Fy, and Mz terms re-calculated.  All others remain zero
      if (myIceParams%numLegs > 1 .and. myIceParams%singleLoad)    &
         iceLoads(:,1) = iceLoadEquivalent(iceLoads, myIceParams)
      
   end function outputCpldCrushLoad

end module iceCpldCrushing
