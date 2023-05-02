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
!  Special type of flexural failure using the Croasdale method
!  Croasdale, K.R. and Cammaert A.B.  "An Improved Method for the Calculation of Ice Loads on Sloping Structures ..."
!  Hydrotechnical Construction Vol. 28 No. 3 ,March, 1994

module IceFlexISO

   use IceFlexBase

   implicit none

   public
   
contains

   subroutine initFlexISO (iceInput, myIceParams, iceLog)

      type(iceInputType), intent(in)            :: iceInput    ! Parameters from input file for initialization
      type(iceFloe_ParameterType), intent(inout)     :: myIceParams ! saved parameters
      type(iceFloe_LoggingType), intent(inout)   :: iceLog      ! structure with message and error logging variables
      type(inputParams)                         :: inParams    ! specific input parameter variable list

!  local variables
      real(ReKi)  :: maxLoad, Hb

!  initialize the common parmeters for flexural ice failure
      call initIceFlex(iceInput, inParams, myIceParams, iceLog)

      call logMessage(iceLog, newLine//' Setting up flexural failure by ISO/Croasdale method ')
      call logMessage(iceLog, ' Based on the ISO 19906 standard equations for max static loads')

      call getIceInput(iceInput, 'peakLoadCOV', inParams%peakLoadCOV, iceLog, 0.1_ReKi, 0.5_ReKi)
      call logMessage(iceLog, ' peakLoadCOV = '//TRIM(Num2LStr(inParams%peakLoadCOV))//' (-)')

      call getIceInput(iceInput, 'coeffLoadPeaks', inParams%coeffLoadPeaks, iceLog, 0.1_ReKi, 1.0_ReKi)
      call logMessage(iceLog, ' coeffLoadPeaks = '//TRIM(Num2LStr(inParams%coeffLoadPeaks)))

      call getIceInput(iceInput, 'coeffLoadMin', inParams%coeffLoadMin, iceLog, 0.0_ReKi, 1.0_ReKi)
      call logMessage(iceLog, ' coeffLoadMin = '//TRIM(Num2LStr(inParams%coeffLoadMin)))

      call getIceInput(iceInput, 'periodCOV', inParams%periodCOV, iceLog, 0.1_ReKi, 0.9_ReKi)
      call logMessage(iceLog, ' periodCOV = '//TRIM(Num2LStr(inParams%periodCOV)))

      call getIceInput(iceInput, 'tauMin', inParams%tauMin, iceLog, 0.1_ReKi, 0.8_ReKi)
      call logMessage(iceLog, ' tauMin = '//TRIM(Num2LStr(inParams%tauMin)))

      call getIceInput(iceInput, 'tauMax', inParams%tauMax, iceLog, inParams%tauMin, 1.0_ReKi)
      call logMessage(iceLog, ' tauMax = '//TRIM(Num2LStr(inParams%tauMax)))

      call getIceInput(iceInput, 'riseTime', inParams%riseTime, iceLog, 0.1_ReKi, 0.9_ReKi)
      call logMessage(iceLog, ' riseTime = '//TRIM(Num2LStr(inParams%riseTime)))
      inParams%fallTime = 1.0 - inParams%riseTime

      call getIceInput(iceInput, 'coeffBreakLength', inParams%coeffBreakLength, iceLog, 3.0_ReKi, 10.0_ReKi)
      call logMessage(iceLog, ' coeffBreakLength = '//TRIM(Num2LStr(inParams%coeffBreakLength)))

      call getIceInput(iceInput, 'rubbleHeight', inParams%rubbleHeight, iceLog, 0.0_ReKi)
      call logMessage(iceLog, ' rubbleHeight = '//TRIM(Num2LStr(inParams%rubbleHeight))//' meters')

      call getIceInput(iceInput, 'rubblePorosity', inParams%rubblePorosity, iceLog, 0.0_ReKi, 1.0_ReKi)
      call logMessage(iceLog, ' rubblePorosity = '//TRIM(Num2LStr(inParams%rubblePorosity)))

      call getIceInput(iceInput, 'rubbleCohesion', inParams%rubbleCohesion, iceLog, 0.0_ReKi)
      call logMessage(iceLog, ' rubbleCohesion strength = '//TRIM(Num2LStr(inParams%rubbleCohesion))//' Pascals')

      call getIceInput(iceInput, 'rubbleAngle', inParams%rubbleAngle, iceLog, 0.0_ReKi, inParams%twr%coneAngle/D2R)
      call logMessage(iceLog, ' rubbleAngle = '//TRIM(Num2LStr(inParams%rubbleAngle))//' degrees')
      inParams%rubbleAngle = D2R*inParams%rubbleAngle   ! Convert to radians

      call getIceInput(iceInput, 'frictionAngle', inParams%frictionAngle, iceLog,      &
                       inParams%rubbleAngle/D2R, inParams%twr%coneAngle/D2R)
      call logMessage(iceLog, ' frictionAngle = '//TRIM(Num2LStr(inParams%frictionAngle))//' degrees')
      inParams%frictionAngle = D2R*inParams%frictionAngle   ! Convert to radians

      call getIceInput(iceInput, 'ice2iceFriction', inParams%ice2iceFriction, iceLog, 0.0_ReKi, 1.0_ReKi)
      call logMessage(iceLog, ' ice2iceFriction = '//TRIM(Num2LStr(inParams%ice2iceFriction)))

      call getIceInput(iceInput, 'waterDensity', inParams%waterDensity, iceLog, 0.0_ReKi)
      call logMessage(iceLog, ' waterDensity = '//TRIM(Num2LStr(inParams%waterDensity))//' kg/m^3')

      call getIceInput(iceInput, 'iceModulus', inParams%iceModulus, iceLog, 0.0_ReKi)
      call logMessage(iceLog, ' iceModulus = '//TRIM(Num2LStr(inParams%iceModulus))//' Pascals')

      call getIceInput(iceInput, 'poissonRatio', inParams%poissonRatio, iceLog, 0.0_ReKi, 0.5_ReKi)
      call logMessage(iceLog, ' poissonRatio = '//TRIM(Num2LStr(inParams%poissonRatio)))

      call getIceInput(iceInput, 'includeHp', inParams%includeHp, iceLog)
      if (inParams%includeHp) then
         call logMessage(iceLog, ' Including pile up term, Hp')
      else
         call logMessage(iceLog, ' NOT including pile up term, Hp')
      endif

      call getIceInput(iceInput, 'includeHl', inParams%includeHl, iceLog)
      if (inParams%includeHl) then
         call logMessage(iceLog, ' Including rubble lifting term, Hl')
      else
         call logMessage(iceLog, ' NOT including rubble lifting term, Hl')
      endif

      call getIceInput(iceInput, 'includeHt', inParams%includeHt, iceLog)
      if (inParams%includeHt) then
         call logMessage(iceLog, ' Including block rotation term, Ht')
      else
         call logMessage(iceLog, ' NOT including block rotation term, Ht')
      endif

      call getIceInput(iceInput, 'includeLc', inParams%includeLc, iceLog)
      if (inParams%includeLc) then
         call logMessage(iceLog, ' Including crack lenght modification term, Lc')
      else
         call logMessage(iceLog, ' NOT including crack lenght modification term, Lc')
      endif

!  The Croasdale method precalculates a time series of loads based
!  on the input parameters and random distributions

!  First calculate the absolute maximum load
!  Various terms are user optional
      maxLoad = 0.0
      if(inParams%includeHp) then
         maxLoad = maxLoad + H_p()
         call logMessage(iceLog, '** Pile up term, Hp = '//TRIM(Num2LStr(H_p()))//' Newtons')
      endif
      if(inParams%includeHl)  then
         maxLoad = maxLoad + H_l()
         call logMessage(iceLog, '** Rubble lifting term, Hl = '//TRIM(Num2LStr(H_l()))//' Newtons')
      endif
      if(inParams%includeHt)  then
         maxLoad = maxLoad + H_t()
         call logMessage(iceLog, '** Block rotation term, Ht = '//TRIM(Num2LStr(H_t()))//' Newtons')
      endif
      if(inParams%includeHr)  then
         maxLoad = maxLoad + H_r()
         call logMessage(iceLog, '** Rubble push through term, Hr = '//TRIM(Num2LStr(H_r()))//' Newtons')
      endif
      if(inParams%includeHb) then
         Hb = H_b()
         maxLoad = (maxLoad + Hb) / (1.0-Hb/(inParams%flexStrength*l_c()*inParams%iceThickness))
         call logMessage(iceLog, '** Breaking load term, Hb = '//TRIM(Num2LStr(Hb))//' Newtons')
      endif
      call logMessage(iceLog, '** Maximum static load for flexural failure = '//TRIM(Num2LStr(maxLoad))//' Newtons')

!  Precalculate the random time series
      call randomFlexLoadTimeSeries(myIceParams, iceLog, maxLoad)

   contains
!--------------------------------------------------------------------
!  Various load terms for max staic load
!  generally DIV by zeros are protected against by input parameter checking
!
!  crack width term
   real(ReKi) function l_c()
      real(ReKi)  :: Lc
      Lc = sqrt(sqrt((inParams%iceModulus*inParams%iceThickness**3)    &
         / (12.0*inParams%waterDensity*grav*(1.0-inParams%poissonRatio**2)) ))
      l_c = inParams%twr%diam
      if (inParams%includeLc) then
         l_c = l_c + 0.25*Lc*pi**2
      endif
   end function l_c

   real(ReKi) function xi()
!  protection against DIV by zero on this function covered by input parameter checking
      xi = (sin(inParams%twr%coneAngle) + inParams%ice2twrFriction * cos(inParams%twr%coneAngle)) /      &
           (cos(inParams%twr%coneAngle) - inParams%ice2twrFriction * sin(inParams%twr%coneAngle))
   end function xi

   real(ReKi) function cot(angle)
      real(ReKi), intent(in) :: angle
      cot = tan(pi/2.0-angle)
   end function cot

!  term with tangents that appears in many of the load terms
   real(ReKi) function angleRatio(theta, alpha)
      real(ReKi), intent(in) :: theta, alpha
      angleRatio = 1.0 - tan(theta)/tan(alpha)
   end function angleRatio

!  Ice breaking load
   real(ReKi) function H_b()
      H_b = 0.68*xi()*inParams%flexStrength*sqrt(sqrt(inParams%waterDensity*grav*inParams%iceThickness**5      &
          / inParams%iceModulus))*l_c()
   end function H_b

!  Hr is the load required to push the ice blocks up the slope of the structure and through the rubble pile
   real(ReKi) function H_r()
      real(ReKi) :: P  ! effective weight of the rubble pile

      P = 0.5*inParams%ice2iceFriction*(inParams%ice2iceFriction+inParams%ice2twrFriction)         &
         *inParams%iceDensity*grav*(1.0-inParams%rubblePorosity)*(inParams%rubbleHeight**2)        &
         *sin(inParams%twr%coneAngle)*(cot(inParams%rubbleAngle)-cot(inParams%twr%coneAngle))    &
         *angleRatio(inParams%rubbleAngle,inParams%twr%coneAngle)

      P = P + 0.5*(inParams%ice2iceFriction+inParams%ice2twrFriction)*inParams%iceDensity*grav           &
                 *(1.0-inParams%rubblePorosity)*(inParams%rubbleHeight**2)*cos(inParams%twr%coneAngle)  &
                 /tan(inParams%twr%coneAngle)*angleRatio(inParams%rubbleAngle,inParams%twr%coneAngle)

      P = P + inParams%rubbleHeight*inParams%iceThickness*inParams%iceDensity*grav                    &
              *(sin(inParams%twr%coneAngle)+inParams%ice2twrFriction*cos(inParams%twr%coneAngle))   &
              /sin(inParams%twr%coneAngle)

      H_r = inParams%twr%diam*P /(cos(inParams%twr%coneAngle)-inParams%ice2twrFriction*sin(inParams%twr%coneAngle))
   end function H_r

!  Load required to push the ice sheet through the rubble
   real(ReKi) function H_p()
      H_p = inParams%twr%diam*(inParams%rubbleHeight**2)*inParams%ice2iceFriction*inParams%iceDensity     &
            *grav*(1.0-inParams%rubblePorosity)*(angleRatio(inParams%rubbleAngle,inParams%twr%coneAngle))**2  &
            /(2.0*tan(inParams%rubbleAngle))
   end function H_p

!  Force required to rotate the ice blocks at the cone to cylinder transition of the tower
   real(ReKi) function H_t()
      H_t = 1.5*inParams%twr%diam*(inParams%iceThickness**2)*inParams%iceDensity*grav*cos(inParams%twr%coneAngle)
      H_t = H_t/(sin(inParams%twr%coneAngle)-inParams%ice2twrFriction*cos(inParams%twr%coneAngle))
   end function H_t

!  Fore required to lift the rubble pile that is above the ice sheet
   real(ReKi) function H_l()
      H_l = cot(inParams%rubbleAngle)-cot(inParams%twr%coneAngle)
      H_l = H_l + tan(inParams%frictionAngle)*angleRatio(inParams%rubbleAngle,inParams%twr%coneAngle)
      H_l = 0.5*H_l*inParams%rubbleHeight*inParams%iceDensity*grav*(1.0-inParams%rubblePorosity) + inParams%rubbleCohesion
      H_l = H_l*xi()*inParams%twr%diam*inParams%rubbleHeight*angleRatio(inParams%rubbleAngle,inParams%twr%coneAngle)
   end function H_l

!****************************************************************************
!  Calculate a sawtooth loading time series with random periods and peaks
   subroutine randomFlexLoadTimeSeries (myIceParams, iceLog, maxLoad)
      real(ReKi), intent(in)                    :: maxLoad
      type(iceFloe_ParameterType), intent(inout)     :: myIceParams
      type(iceFloe_LoggingType), intent(inout)   :: iceLog   ! structure with message and error logging variables

      integer(IntKi), parameter                 :: LuxLevel = 3
      integer(IntKi) :: err
      integer(IntKi) :: n, ns, nL, nSteps, nRiseSteps, nFallSteps, minLoadSteps
      real(ReKi)     :: meanLoad    ! Mean of the instantaneous peaks of the cycles
      real(ReKi)     :: peakLoad    ! Instantaneous peak for a single cycle
      real(ReKi)     :: minLoad     ! Minimum possible load
      real(ReKi)     :: meanPeriod  ! Mean loading/unloading/dwell cycle period, seconds
      real(ReKi)     :: period      ! randomly distributed period for a single cycle

      real(ReKi)     :: tau(1)  ! for passing into random number routine
      
!  Initialize the random number generator
      CALL RLuxGo ( LuxLevel, abs(inParams%randomSeed), 0, 0 )

!  Calc overall number of steps
      nSteps = size(myIceParams%loadSeries,1)

!  Calculate the mean period and loads
      meanPeriod = inParams%coeffBreakLength*inParams%iceThickness/inParams%iceVelocity
      minLoad    = inParams%coeffLoadMin*maxLoad
      meanLoad   = minLoad*(1.0-inParams%coeffLoadPeaks) + inParams%coeffLoadPeaks*maxLoad

!  Loop on all legs
      do nL = 1, myIceParams%numLegs

   !  Number of whole periods required not known (since period lenght is random)
   !  So iterate until we have enough points in the time series
         n = 1
         seriesLoop: do while (n < nSteps+1)      
   !  Period is time from no load up to peak, down, then dwell at minimum (normal distribution)
            CALL RndNorm( period, meanPeriod, inParams%periodCOV*meanPeriod )
   !  Period has to be limited to +/- 50% of the mean period
            period = min(1.5_ReKi*meanPeriod, max(0.5_ReKi*meanPeriod, period))

   !  sub period is the fraction of a period: time for load to go up then down (uniform distribution)
            CALL RanLux ( tau )
   !  adjust the sub period to be on the requested interval
            tau(1) = inParams%tauMin + (inParams%tauMax - inParams%tauMin)*tau(1)

   !  Peak instantaneous load for the period (normal distribution)
            CALL RndNorm( peakLoad, meanLoad, inParams%peakLoadCOV*meanLoad )
   !  adjust the peak load to the absolute min and instantaneous max
            peakLoad = min(maxLoad, max(minLoad, peakLoad))

   !  begin to increase the load from the minimum to the instantaneous peak
            nRiseSteps = nint( inParams%riseTime * tau(1) * period / myIceParams%dt)
            if (nRiseSteps == 0) then   ! if dt is too large compared to the rise time then get div by zero below
               call iceErrorHndlr (iceLog, ErrID_Fatal, 'Time step too large or period too short in RandomFlexLoadTimeSeries', 1)
               return
            endif
            do ns = 1, nRiseSteps
               myIceParams%loadSeries(n,nL) = minLoad + (peakLoad-minLoad)*float(ns-1)/float(nRiseSteps)
               n = n + 1
               if (n > nSteps) exit seriesLoop
            enddo
            
   !  begin to decrease the load from the peak down to the minimum
            nFallSteps = nint( inParams%fallTime * tau(1) * period / myIceParams%dt)
            if (nFallSteps == 0) then   ! if dt is too large compared to the fall time then get div by zero below
               call iceErrorHndlr (iceLog, ErrID_Fatal, 'Time step too large or period too short in RandomFlexLoadTimeSeries', 1)
               return
            endif
            do ns = 1, nFallSteps
               myIceParams%loadSeries(n,nL) = peakLoad - (peakLoad-minLoad)*float(ns-1)/float(nFallSteps)
               n = n + 1
               if (n > nSteps) exit seriesLoop
            enddo

   !  Fill the remaining time steps in the period with the minimum load
            minLoadSteps = nint( period / myIceParams%dt) - nRiseSteps - nFallSteps
            do ns = 1, minLoadSteps
               myIceParams%loadSeries(n,nL) = minLoad
               n = n + 1
               if (n > nSteps) exit seriesLoop
            enddo

         enddo seriesLoop

   !  Apply shelter factor to individual legs
         myIceParams%loadSeries(:,nL) = myIceParams%ks(nL)*myIceParams%loadSeries(:,nL)

      enddo ! leg loop

   end subroutine randomFlexLoadTimeSeries 

   end subroutine initFlexISO  ! containing subroutine

!--------------------------------------------------------------------
!  get load for requested time
   function outputFlexLoadISO (myIceParams, iceLog, time) result(iceLoads)
      type(iceFloe_ParameterType), intent(in)      :: myIceParams
      type(iceFloe_LoggingType), intent(inout) :: iceLog   ! structure with message and error logging variables
      real(DbKi), intent(in)                  :: time
      real(ReKi)                              :: iceLoads(6,myIceParams%numLegs)

      iceLoads = outputFlexLoad (myIceParams, iceLog, time)
   end function outputFlexLoadISO

end module IceFlexISO
