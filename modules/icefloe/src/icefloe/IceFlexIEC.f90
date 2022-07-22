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

!*****************************************************************************************
!  Special type of loading due to ice sheet flexural failure using the Ralston method
!  Ralston, T, "Ice Force Design Considerations for Conical Offshore Structures", POAC, 1977
!  As specified in IEC 61400-3 for offshore turbines
module iceFlexIEC

   use IceFlexBase

   implicit none

   public
   
contains

!============================================================================
!  calculate the time series of ice loads using the Ralston method
 
   subroutine initFlexIEC (iceInput, myIceParams, iceLog)

      type(iceInputType), intent(in)            :: iceInput
      type(iceFloe_ParameterType), intent(inout)     :: myIceParams
      type(iceFloe_LoggingType), intent(inout)   :: iceLog   ! structure with message and error logging variables
      type(inputParams)                         :: inParams ! specific input parameter variable list

!  local variables
      real(ReKi)     :: maxLoad  ! maximum load
      real(ReKi)     :: freq     ! frequency of sinusoidal load
      integer(IntKi) :: nL !err, 

!  initialize the common parmeters for flexural ice failure
      call initIceFlex(iceInput, inParams, myIceParams, iceLog)

      call logMessage(iceLog, newLine//' Setting up flexural failure by Ralston method ')
      call logMessage(iceLog, ' Use the IEC 61400-3 suggested equations for max static loads')

      call getIceInput(iceInput, 'rideUpThickness', inParams%rideUpThickness, iceLog, 0.0_ReKi)
      call logMessage(iceLog, ' Ride up thickness = '//TRIM(Num2LStr(inParams%rideUpThickness))//' meters')

      call getIceInput(iceInput, 'twrConeTopDiam', inParams%twr%coneTopDiam, iceLog, 0.0_ReKi)
      call logMessage(iceLog, ' Tower cone top diameter = '//TRIM(Num2LStr(inParams%twr%coneTopDiam))//' meters')

      call getIceInput(iceInput, 'freqParamK', inParams%freqParamK, iceLog, 4.0_ReKi, 7.0_ReKi)
      call logMessage(iceLog, ' Frequency parameter K = '//TRIM(Num2LStr(inParams%freqParamK))//' ')

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

!  The IEC/Ralston method precalculates a time series of loads based
!  as a sinusoid at the fundamental (or a specified) frequency

!  First calculate the absolute maximum load
!  Various terms are user optional
      maxLoad = 0.0
      if(inParams%includeHr) maxLoad = maxLoad + H_r()
      if(inParams%includeHb) maxLoad = (maxLoad + H_b())
      call logMessage(iceLog, '** Maximum static load for flexural failure = '//TRIM(Num2LStr(maxLoad))//'Newtons')

!  Precalculate the random time series
!   with sinusoid per IEC
      freq = inParams%iceVelocity/inParams%iceThickness/inParams%freqParamK
      call IECLoadTimeSeries(myIceParams, inParams, iceLog, maxLoad, freq)

   contains
!==================================================================

! Complete elliptical integral of the first kind
! Approximation from "Handbook of Mathematical Functions", 
!   U.S. Dept of Commerce, Nat. Bureau of Standards
!   Ed. by Milton Abramovitch, June 1963
   function elliptic_1(alpha) result(K)
      real(ReKi), intent(in)  :: alpha
      real(ReKi)              :: m1
      real(ReKi)              :: K

      m1 = 1.0 - sin(alpha)**2  
      if (m1 == 0.0) then  ! not likely, but just in case
         K = huge(K)
         return
      endif
      
      K = 1.3862944 + 0.1119723*m1 + 0.0725296*m1**2
      K = K + log(1.0/m1)*(0.5 + 0.1213478*m1 + 0.0288729*m1**2)
   end function elliptic_1

! Complete elliptical integral of the second kind
! Approximation from "Handbook of Mathematical Functions", 
!   U.S. Dept of Commerce, Nat. Bureau of Standards
!   Ed. by Milton Abramovitch, June 1963
   function elliptic_2(alpha) result(E)
      real(ReKi), intent(in)  :: alpha
      real(ReKi)              :: m1
      real(ReKi)              :: E

      m1 = 1.0 - sin(alpha)**2  
      if (m1 == 0.0) then  ! not likely, but just in case
         E = 1.0
         return
      endif

      E = 1.0 + 0.4630151*m1 + 0.1077812*m1**2
      E = E + log(1.0/m1)*(0.2452727*m1 + 0.0412496*m1**2)
   end function elliptic_2

   real(ReKi) function g_r()
!    No div by zero for towerConeAngle between 20 and 70 degrees as limited on input
       g_r = (sin(inParams%twr%coneAngle) + inParams%twr%coneAngle/cos(inParams%twr%coneAngle)) /     &
             (2.0*inParams%ice2twrFriction*inParams%twr%coneAngle*cos(inParams%twr%coneAngle) +                 &
             0.5*pi*sin(inParams%twr%coneAngle)**2)
   end function g_r

!------------------------------------------------------------------------
!  Ice breaking load
   real(ReKi) function H_b()
      real(ReKi)  :: Y, G, x, term1, term2, term3

      Y = 2.711  ! Based on Tresca yield criteria for the ice sheet
      G = inParams%iceDensity*grav*inParams%twr%diam**2    &      ! ice weight
        / (4.0*inParams%flexStrength*inParams%iceThickness)            ! ice strength

      x = 1.0 + 1.0/sqrt(3.0*G+0.5*Y)   ! x will always be > 1 so no problems in below equations

      term1 = inParams%flexStrength*inParams%iceThickness**2/3.0
      term2 = tan(inParams%twr%coneAngle)/(1.0-g_r()*inParams%ice2twrFriction)  ! 1-g_r*mu can be zero for large angles and friction but input checking should cover this
      term3 = G*(x-1.0)*(x+2.0) + (1.0+Y*x*log(x))/(x-1.0)

      H_b = term1*term2*term3
   end function H_b

!  Hr is the load required to push the ice blocks up the slope of the structure and through the rubble pile
   real(ReKi) function H_r()
      real(ReKi) :: W, f, numer, denom

      W     = inParams%iceDensity*grav*inParams%rideUpThickness*(inParams%twr%diam**2 - inParams%twr%coneTopDiam**2)   &
             /4.0/cos(inParams%twr%coneAngle)
      f     = inParams%ice2TwrFriction*elliptic_1(inParams%twr%coneAngle)*cos(inParams%twr%coneAngle)              &
            + sin(inParams%twr%coneAngle)
      numer = tan(inParams%twr%coneAngle) + inParams%ice2TwrFriction*(elliptic_2(inParams%twr%coneAngle)           &
            - f*g_r()*cos(inParams%twr%coneAngle))
      denom = (1.0-g_r()*inParams%ice2TwrFriction)   ! 1-g_r*mu can be zero for large angles and friction but input checking should cover this

      H_r = W*numer/denom
   end function H_r

   end subroutine initFlexIEC

!--------------------------------------------------------------------
!  get load for requested time
   function outputFlexLoadIEC (myIceParams, iceLog, time) result(iceLoads)
      type(iceFloe_ParameterType), intent(in)      :: myIceParams
      type(iceFloe_LoggingType), intent(inout) :: iceLog   ! structure with message and error logging variables
      real(DbKi), intent(in)                  :: time
      real(ReKi)                              :: iceLoads(6,myIceParams%numLegs)

      iceLoads = outputFlexLoad (myIceParams, iceLog, time)
   end function outputFlexLoadIEC

end module iceFlexIEC

