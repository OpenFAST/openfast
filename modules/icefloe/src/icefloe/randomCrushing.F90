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

!  Module to calculate a time series of random loads due to continuous ice crushing
!  
!  Based on a spectral method from:
!  Karna, T, et. al., "A Spectral Model for Forces due to Ice Crushing", 
!  Journal of Offshore Mechanics and Arctic Engineering, May, 2006
!
module randomCrushing

   use iceCrushingISO

   implicit none

   public
   
contains

!================================================================================
!  Initialize the time series of ice loads for random continuous crushing
   subroutine initRandomCrushing (iceInput, myIceParams, iceLog)

      type(iceInputType), intent(in)            :: iceInput    ! Parameters read from file for initialization
      type(iceFloe_ParameterType), intent(inout)     :: myIceParams ! saved parameters
      type(iceFloe_LoggingType), intent(inout)   :: iceLog      ! structure with message and error logging variables
      type(inputParams)                         :: inParams    ! specific input parameter variable list

!  local variables
      real(ReKi)  :: maxLoad

!  initialize the common parmeters for flexural ice failure
      call initIceCrushISO(iceInput, inParams, myIceParams, iceLog)

      call logMessage(iceLog, newLine//' Setting parameters for random crushing loads time series')

      call getIceInput(iceInput, 'coeffPSD_b', inParams%coeffPSD_b, iceLog, 0.1_ReKi, 3.0_ReKi)
      call logMessage(iceLog, ' coeffPSD_b = '//TRIM(Num2LStr(inParams%coeffPSD_b)))

      call getIceInput(iceInput, 'coeffPSD_Ks', inParams%coeffPSD_Ks, iceLog, 1.0_ReKi, 5.0_ReKi)
      call logMessage(iceLog, ' coeffPSD_Ks = '//TRIM(Num2LStr(inParams%coeffPSD_Ks)))

      call getIceInput(iceInput, 'crushLoadCOV', inParams%crushLoadCOV, iceLog, 0.1_ReKi, 1.0_ReKi)
      call logMessage(iceLog, ' crushLoadCOV = '//TRIM(Num2LStr(inParams%crushLoadCOV)))

      call getIceInput(iceInput, 'stdLoadMult', inParams%stdLoadMult, iceLog, 1.0_ReKi, 6.0_ReKi)
      call logMessage(iceLog, ' stdLoadMult (number of std devs) = '//TRIM(Num2LStr(inParams%stdLoadMult)))

      call getIceInput(iceInput, 'freqStep', inParams%freqStep, iceLog, .001_ReKi, 0.1_ReKi)
      call logMessage(iceLog, ' Frequency step = '//TRIM(Num2LStr(inParams%freqStep))//' Hz')

!  Precalculates a time series of loads
      maxLoad = globalCrushLoadISO(inParams)
      call logMessage(iceLog, '** Global crushing load is: '//TRIM(Num2LStr(maxLoad))//' Newtons.' )
      call randomCrushLoadTimeSeries(myIceParams, iceLog, maxLoad)
   
   contains

!==================================================================
! Generate a random time series by randomizing the phase of sinusoids
! based on a spectrum
   subroutine randomCrushLoadTimeSeries (myIceParams, iceLog, maxLoad)
      type(iceFloe_ParameterType), intent(inout)     :: myIceParams
      type(iceFloe_LoggingType), intent(inout)   :: iceLog   ! structure with message and error logging variables
      real(ReKi), intent(in)                    :: maxLoad  ! global crushing maximum load

      integer(IntKi), parameter  :: LuxLevel = 3            ! for ranlux
      integer(IntKi)             :: err
      integer(IntKi)             :: nf, ns, nL, nSteps, nfSteps
      real(ReKi)                 :: timeStep, meanLoad, stdLoad
      real(ReKi), allocatable    :: frequency(:), amplitude(:), randPhase(:)   ! temporary variables
      real(ReKi)                 :: nyqFreq, a, stdSum

!  get the mean and standard deviation of the loading
!  these will be used to scale the generated time series
      meanLoad = maxLoad/(1.0+inParams%stdLoadMult*inParams%crushLoadCOV)
      call logMessage(iceLog, '** Mean crushing load is: '//TRIM(Num2LStr(meanLoad))//' Newtons.' )
      stdLoad  = meanLoad*inParams%crushLoadCOV
      call logMessage(iceLog, '** Stdev of crushing load is: '//TRIM(Num2LStr(stdLoad))//' Newtons.' )
      
!  Number of steps in time series of loads
      nSteps = size(myIceParams%loadSeries,1)

!  allocate some temporary arrays
      nyqFreq = 0.5/myIceParams%dt
      nfSteps = floor(nyqFreq/inParams%freqStep)
      allocate(randPhase(nfSteps), stat=err)
      if (err /= 0) then
         call iceErrorHndlr (iceLog, ErrID_Fatal, 'Error in phase array allocation in RandomCrushLoadTimeSeries', 1)
         return  !  error in array allocation
      endif   
      allocate(frequency(nfSteps), stat=err)
      allocate(amplitude(nfSteps), stat=err)
      if (err /= 0) then
         call iceErrorHndlr (iceLog, ErrID_Fatal, 'Error in frequency or amplitude array allocation in '//  &
                                                  'RandomCrushLoadTimeSeries', 1)
         return  !  error in array allocation
      endif   

!  Initialize the random number generator
      call RLuxGo ( LuxLevel, abs(inParams%randomSeed), 0, 0 )

!  Loop on number of legs
      do nL = 1, myIceParams%numLegs
   !  get a random phase on the interval +/- pi
         call RanLux ( randPhase )
         randPhase = 2.0*pi*randPhase - pi

   !  calculate the amplitude for each frequency from the PSD
         a = inParams%coeffPSD_b/(inParams%iceVelocity**0.6)
         do nf = 1, nfSteps
            frequency(nf) = float(nf)*inParams%freqStep
   !  Note PSD is amplitude squared so take square root and mult by freq resolution to get discrete
   !  Also times 2 since we are using a single sided PSD
            amplitude(nf) = sqrt( (2.0*a*(stdLoad**2)*inParams%freqStep)            &
                             / (1.0+inParams%coeffPSD_Ks*(a**1.5)*frequency(nf)**2) )
         enddo

   !  calculate time series as a sum of frequencies at each time step 
   !  with amplitudes and random phases at each frequency
   !  get the standard deviation for later scaling
         stdSum = 0.0
         do ns = 1, nSteps
            timeStep = myIceParams%dt*real(ns,ReKi)
            myIceParams%loadSeries(ns,nL) = 0.0
            do nf = 1, nfSteps
               myIceParams%loadSeries(ns,nL) = myIceParams%loadSeries(ns,nL) +         &
                       amplitude(nf)*cos(2.0*pi*frequency(nf)*timeStep+randPhase(nf))
            enddo
!            stdSum = stdSum + myIceParams%loadSeries(ns,nL)**2
         enddo
   !  scale the variations to get the desired standard deviation
   !  Testing suggests that this doesn't give the desired PSD, so remove
!         stdSum = sqrt(stdSum/real(nSteps,ReKi))
!         myIceParams%loadSeries(:,nL) = (stdLoad/stdSum)*myIceParams%loadSeries(:,nL)
   !  add in the mean
         myIceParams%loadSeries(:,nL) = myIceParams%ks(nL)*(meanLoad + myIceParams%loadSeries(:,nL))

      enddo   ! leg loop

      deallocate(frequency)
      deallocate(amplitude)
      deallocate(randPhase)

   end subroutine randomCrushLoadTimeSeries

   end subroutine initRandomCrushing   ! containing subroutine

!****************************************************************************
!  Continuous crushing uses the standard inerpolation routine
!  of the precalculated time series
   function outputRandomCrushLoad (myIceParams, iceLog, time)  result(iceLoads)
      type(iceFloe_ParameterType), intent(in)      :: myIceParams
      type(iceFloe_LoggingType), intent(inout) :: iceLog   ! structure with message and error logging variables
      real(DbKi), intent(in)                  :: time
      real(ReKi)                              :: iceLoads(6,myIceParams%numLegs)

      iceLoads = outputCrushLoadISO (myIceParams, iceLog, time)
   end function outputRandomCrushLoad

end module randomCrushing