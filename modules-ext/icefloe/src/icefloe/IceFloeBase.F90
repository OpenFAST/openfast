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
! base type for ice floes, contains generic
! initialization, update, and utility routines


module IceFloeBase

   use precision
   use IceFloe_Types
   use IceInputParams
   use Ran_Lux_Mod
   use NWTC_IO, only : Num2LStr

   implicit none

   public

! This is compared to gravity in FAST, warning if different
   real(ReKi), parameter :: grav = 9.81

!  ice type parameters
   integer(IntKi), parameter :: randomCrush    = 1
   integer(IntKi), parameter :: interCrush     = 2
   integer(IntKi), parameter :: lockInISO      = 3
   integer(IntKi), parameter :: crushIEC       = 4
   integer(IntKi), parameter :: cpldCrush      = 5
   integer(IntKi), parameter :: flexFailISO    = 6
   integer(IntKi), parameter :: flexFailIEC    = 7

!  for input checking
   integer(IntKi), parameter :: lowTypeLimit   = min(randomCrush, interCrush, lockInISO, crushIEC,    &
                                                     cpldCrush, flexFailISO, flexFailIEC)
   integer(IntKi), parameter :: hiTypeLimit    = max(randomCrush, interCrush, lockInISO, crushIEC,    &
                                                     cpldCrush, flexFailISO, flexFailIEC)

!  structure to hold some tower data



!----------------------------------------------------------------------------------------------------
contains

!  Generic initialization of the most common parameters
!  Only function is to assign input parameters to convenient variables in module
   subroutine initIceFloe (iceInput, inParams, myIceParams, iceLog)
      type(iceInputType), intent(in)               :: iceInput    ! Parameters from input file for initialization
      type(iceFloe_ParameterType), intent(inout)   :: myIceParams ! saved parameters
      type(iceFloe_LoggingType), intent(inout)     :: iceLog      ! structure with message and error logging variables
      type(inputParams), intent(out)            :: inParams    ! specific input parameter variable list
      integer(IntKi)                            :: err
      integer(IntKi)                            :: nL
      real(ReKi)                                :: spacing     ! get leg spacing from x,y coords

      call logMessage(iceLog, newLine//' Setting common ice flow input parameters ')

      call getIceInput(iceInput, 'iceType', inParams%iceType, iceLog, lowTypeLimit, hiTypeLimit)
      call logMessage(iceLog, ' Ice type = '//TRIM(Num2LStr(inParams%iceType)))

      call getIceInput(iceInput, 'iceThickness', inParams%iceThickness, iceLog, 0.001_ReKi, 100.0_ReKi)
      call logMessage(iceLog, ' Ice thickness = '//TRIM(Num2LStr(inParams%iceThickness))//' meters')

      call getIceInput(iceInput, 'iceVelocity', inParams%iceVelocity, iceLog, 0.001_ReKi, 10.0_ReKi)
      call logMessage(iceLog, ' Ice velocity = '//TRIM(Num2LStr(inParams%iceVelocity))//' m/s')

      call getIceInput(iceInput, 'iceDirection', inParams%iceDirection, iceLog, -360.0_ReKi, 360.0_ReKi)
      call logMessage(iceLog, ' Ice direction = '//TRIM(Num2LStr(inParams%iceDirection))//' degrees')
      inParams%iceDirection = D2R*inParams%iceDirection
      myIceParams%iceDirection = inParams%iceDirection  ! must save this parameter for later output calcs
      
      call getIceInput(iceInput, 'timeStep', inParams%timeStep, iceLog, 0.0_ReKi)
      call logMessage(iceLog, ' Time step = '//TRIM(Num2LStr(inParams%timeStep))//' sec')

      call getIceInput(iceInput, 'towerDiameter', inParams%twr%diam, iceLog, 0.1_ReKi, 100.0_ReKi)
      call logMessage(iceLog, ' Tower diameter = '//TRIM(Num2LStr(inParams%twr%diam))//' meters')

      call getIceInput(iceInput, 'randomSeed', inParams%randomSeed, iceLog, 0)
      call logMessage(iceLog, ' Random seed = '//TRIM(Num2LStr(inParams%randomSeed)))

      allocate (inParams%twr%leg(myIceParams%numLegs), stat=err)
      if (err /= 0) then
         call iceErrorHndlr (iceLog, ErrID_Fatal, 'Error in allocation of support leg array', 1)
         return  !  error in array allocation
      endif   

      if (myIceParams%numLegs > 1) then
         call logMessage(iceLog, ' Number of legs in support structure = '//TRIM(Num2LStr(myIceParams%numLegs)))

         call getIceInput(iceInput, 'singleLoad', inParams%singleLoad, iceLog)
         myIceParams%singleLoad = inParams%singleLoad  ! must save this parameter for later use
         if (inParams%singleLoad) then
            call logMessage(iceLog, ' Apply multitple leg loads as a single equivalent load')
         else
            call logMessage(iceLog, ' Apply multitple leg loads individually')
         endif

   !  get leg locations 
         do nL = 1, myIceParams%numLegs
            call getIceInput(iceInput, 'legX'//TRIM(Num2LStr(nL)), inParams%twr%leg(nL)%X, iceLog)
            call getIceInput(iceInput, 'legY'//TRIM(Num2LStr(nL)), inParams%twr%leg(nL)%Y, iceLog)
            call logMessage(iceLog, ' Leg '//TRIM(Num2LStr(nL))//' position is ('                                  &
              //TRIM(Num2LStr(inParams%twr%leg(nL)%X))//', '//TRIM(Num2LStr(inParams%twr%leg(nL)%Y))//') meters')
            myIceParams%legX(nL) = inParams%twr%leg(nL)%X
            myIceParams%legY(nL) = inParams%twr%leg(nL)%Y
         enddo

   !  get minimum leg spacing
         spacing = sqrt((inParams%twr%leg(1)%X-inParams%twr%leg(2)%X)**2 + (inParams%twr%leg(1)%Y-inParams%twr%leg(2)%Y)**2)
         if (myIceParams%numLegs > 2) spacing = min(spacing, sqrt((inParams%twr%leg(2)%X-inParams%twr%leg(3)%X)**2 +   &
                                                                  (inParams%twr%leg(2)%Y-inParams%twr%leg(3)%Y)**2))
         if (myIceParams%numLegs > 3) spacing = min(spacing, sqrt((inParams%twr%leg(4)%X-inParams%twr%leg(3)%X)**2 +   &
                                                                  (inParams%twr%leg(4)%Y-inParams%twr%leg(3)%Y)**2))
         spacing = min(spacing, sqrt((inParams%twr%leg(1)%X-inParams%twr%leg(myIceParams%numLegs)%X)**2 +   &
                                     (inParams%twr%leg(1)%Y-inParams%twr%leg(myIceParams%numLegs)%Y)**2))
                                     
         if (spacing/inParams%twr%diam < 4.0) then
               call iceErrorHndlr (iceLog, ErrID_Warn,                     &
               'Leg spacing/diameter is < 4, high potential for jamming forces which are not considered by this routine', 1)
         endif

!   Auto calculate or assign from user input leg load shelter factors
         call getIceInput(iceInput, 'legAutoFactor', inParams%legAutoFactor, iceLog)
         if (inParams%legAutoFactor) then
            call logMessage(iceLog, ' Leg load shelter factors determined automatically')
         else
            call logMessage(iceLog, ' Leg load shelter factors assigned by user')
         endif
         do nL = 1, myIceParams%numLegs
            myIceParams%ks(nL) = shelterFactor(myIceParams%numLegs, spacing,           &
                                 inParams%twr%leg(nL)%X, inParams%twr%leg(nL)%Y, inParams%iceDirection)
            call logMessage(iceLog, ' Auto calculated shelter factor for leg #'//TRIM(Num2LStr(nL))// &
                                    ' = '//TRIM(Num2LStr(myIceParams%ks(nL))))
            if (.not. inParams%legAutoFactor) then
               call getIceInput(iceInput, 'shelterFactor_ks'//TRIM(Num2LStr(nL)), inParams%twr%leg(nL)%ks, iceLog, 0.0_ReKi, 1.0_ReKi)
               call logMessage(iceLog, ' User specified shelter factor for leg #'//TRIM(Num2LStr(nL))//  &
                                       ' = '//TRIM(Num2LStr(inParams%twr%leg(nL)%ks)))
!     Compare auto generated to user assigned factors, warn if big difference
               if (abs(inParams%twr%leg(nL)%ks - myIceParams%ks(nL)) > 0.1) then
                  call iceErrorHndlr (iceLog, ErrID_Warn,                     &
                  'User assigned shelter factor is more then 0.1 different than auto calculated factor - please check', 1)
               endif
               myIceParams%ks(nL) = inParams%twr%leg(nL)%ks
            endif
         enddo

      else
         myIceParams%ks(1) = 1.0
         myIceParams%singleLoad = .true.
      endif

   end subroutine initIceFloe


!=======================================================================
! Use some trigonometry to assign shelter factors to the individual legs
! Based on iceFloe directon and leg coordinates
! This algorithm assumes the legs are spaced evenly around the tower centroid at (0,0)
   function shelterFactor (nLegs, spacing, X, Y, iceDir) result(factor) 
      
      integer(IntKi), intent(in) :: nLegs
      real(ReKi), intent(in)     :: X, Y, spacing, iceDir   ! ice floe direction in radians
      real(ReKi)                 :: factor
      real(ReKi)                 :: angle, rotX, rotY, loAngle, hiAngle

!  For large spacing legs are independent and factor = 1
      if (spacing > 5.0) then
         factor = 1.0
         return
      endif

! calculate the angle between the ice Floe direction vector and 
! the vector from the tower centroid to the specified leg
      rotX =  X*cos(iceDir) + Y*sin(iceDir)
      rotY = -X*sin(iceDir) + Y*cos(iceDir)
      angle = mod(R2D*atan2(rotY,rotX)+360.0_ReKi,360.0_ReKi)

! assign angle limits dependent on number of legs
      if (nLegs < 4) then
         loAngle = 60.0
         hiAngle = 300.0
      else
         loAngle = 67.5
         hiAngle = 292.5
      endif
      
! assign the factors
      factor = 1.0
      if (angle < loAngle .or. angle >= hiAngle) factor = 0.0

   end function shelterFactor   
   
!=======================================================================
! Calculate a sinusoidal load time series per IEC 61400-3 Ed 1 Annex E section E.4.6
   subroutine IECLoadTimeSeries (myIceParams, inParams, iceLog, maxLoad, freq)
      real(ReKi), intent(in)                    :: maxLoad, freq
      type(iceFloe_ParameterType), intent(inout)     :: myIceParams
      type(inputParams), intent(in)             :: inParams ! specific input parameter variable list
      type(iceFloe_LoggingType), intent(inout)   :: iceLog   ! structure with message and error logging variables
!      integer(IntKi) :: err
      integer(IntKi) :: ns, nL, nSteps
      real(ReKi)     :: timeStep
   
!  Calc overall number of steps
      nSteps = size(myIceParams%loadSeries,1)

      do ns = 1, nSteps
         timeStep = myIceParams%dt*real(ns,ReKi)
         do nL = 1, myIceParams%numLegs
            myIceParams%loadSeries(ns,nL) = myIceParams%ks(nL)*maxLoad*(0.75 +         &
                      0.25*sin(2.0*pi*freq*timeStep+inParams%twr%leg(nL)%phase))
         enddo
      enddo
   
   end subroutine IECLoadTimeSeries 

!=======================================================================
!  Generic update function - really just an interpolation
!  routine for the precalculated load time series
   function outputIceLoads (myIceParams, iceLog, time)  result(iceLoads)
      type(iceFloe_ParameterType), intent(in)        :: myIceParams
      type(iceFloe_LoggingType),  intent(inout)  :: iceLog   ! structure with message and error logging variables
      real(DbKi), intent(in)  :: time
      real(ReKi)              :: iceLoads(6,myIceParams%numLegs)
      integer(IntKi)          :: nStepLow, nStepHigh
      real(ReKi)              :: iceForceMag    ! load in Newtons in direction of ice flow
      integer(IntKi)          :: nL

      do nL = 1, myIceParams%numLegs
   !  find the time steps before and after time
         nStepLow  = floor(time/myIceParams%dt) + 1  ! +1 so that time = 0 is point #1
         if (nStepLow < 1) then
            nStepLow = 1
            call iceErrorHndlr (iceLog, ErrID_Warn, 'Calling program is attempting to obtain an '//   &
                                                    'ice load value for time < 0, '//                 &
                                                    't=0 value will be used', 0 )
         endif
         if (nStepLow > size(myIceParams%loadSeries,1)) nStepLow = size(myIceParams%loadSeries,1)
         nStepHigh = min(nStepLow+1, size(myIceParams%loadSeries,1))
         
         if(nStepLow > size(myIceParams%loadSeries,1)) then
            nStepLow = size(myIceParams%loadSeries,1)
            nStepHigh = nStepLow
            call iceErrorHndlr (iceLog, ErrID_Warn, 'Calling program is attempting to obtain an '//   &
                                                    'ice load value past the end time, '//            &
                                                    't=endTime value will be used', 0 )
         endif
         
   !  calculate the magnitude of the load
         iceForceMag = (myIceParams%loadSeries(nStepHigh,nL) - myIceParams%loadSeries(nStepLow,nL))       &
                     * (time/myIceParams%dt - float(nStepLow-1)) + myIceParams%loadSeries(nStepLow,nL)

   !  apply the load in the proper direction via x,y components
         iceLoads(:,nL) = iceLoadDirection(iceForceMag, myIceParams)
      enddo

   !  If the loads from multiple legs are applied as an "effective" load on one leg at the
   !  support center, then compute that effective load vector and put into the first load vector
   !  Note that only Fx, Fy, and Mz terms re-calculated.  All others remain zero
      if (myIceParams%numLegs > 1 .and. myIceParams%singleLoad)    &
         iceLoads(:,1) = iceLoadEquivalent(iceLoads, myIceParams)

   end function outputIceLoads

!===================================================================================================
!  Function to handle ice forces from a specified direction relative to ground coordinates (in horizontal plane)
   function iceLoadDirection(iceForceMag, myIceParams)   result(loadVect)
      real(ReKi), intent(in)             :: iceForceMag
      type(iceFloe_ParameterType), intent(in) :: myIceParams
      real(ReKi)                         :: loadVect(6)

      loadVect(1) = iceForceMag * cos(myIceParams%iceDirection)
      loadVect(2) = iceForceMag * sin(myIceParams%iceDirection)
      loadVect(3) = 0.0
      loadVect(4) = 0.0
      loadVect(5) = 0.0
      loadVect(6) = 0.0

   end function iceLoadDirection

!===================================================================================================
!  Function to calculate a single equivalent load vector from multiple leg locad vectors
   function iceLoadEquivalent(inLoads, myIceParams)   result(equiv)
      type(iceFloe_ParameterType), intent(in) :: myIceParams
      real(ReKi)                         :: inLoads(6,myIceParams%numLegs)
      real(ReKi)                         :: equiv(6)
      integer(IntKi)                     :: nL

      equiv = 0.0
      do nL = 1, myIceParams%numLegs
         equiv(6) = equiv(6) + inLoads(2,nL)*myIceParams%legX(nL) - inLoads(1,nL)*myIceParams%legY(nL)
         equiv(1) = equiv(1) + inLoads(1,nL)
         equiv(2) = equiv(2) + inLoads(2,nL)
      enddo

   end function iceLoadEquivalent

!  Some utility routines

!=======================================================================
!  Generate random numbers on a normal distribution
!  Borrowed from Turbsim v1.50

   SUBROUTINE RndNorm( RandNormNum, mu, sigma )

   IMPLICIT                            NONE

   REAL(ReKi), INTENT(OUT)          :: RandNormNum    ! Normally distributed numbers
   REAL(ReKi), INTENT(IN), OPTIONAL :: mu             ! mean of the distributed numbers - DEFAULT IS 0.0
   REAL(ReKi), INTENT(IN), OPTIONAL :: sigma          ! standard deviation of the distributed numbers - DEFAULT IS 1.0

 ! Internal variable
   REAL(ReKi)                       :: RN(2)          ! Pairs of random numbers
!   integer(IntKi)                   :: n

 ! Generate a normal distribution on (0,1) from a uniform distribution ( ACTUALLY [0,1) )
      CALL RanLux(RN)

      RandNormNum = SQRT(PI/8.0)*LOG((1.0+RN(1))/(1.0-RN(1)))
      IF (RN(2) < 0.5) RandNormNum = -RandNormNum

 ! Give the correct mean and standard deviation, if specified
      IF ( PRESENT( sigma ) ) THEN
         RandNormNum = RandNormNum * sigma
      ENDIF

      IF ( PRESENT( mu ) ) THEN
         RandNormNum = RandNormNum + mu
      ENDIF

   END SUBROUTINE RndNorm
   
end module IceFloeBase
