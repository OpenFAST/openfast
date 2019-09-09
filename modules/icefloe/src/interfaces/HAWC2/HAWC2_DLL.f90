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

!**********************************************************************************************************************************
! File last committed: $Date: 2014-06-24 10:28:50 -0700 (Tue, 24 Jun 2014) $
! (File) Revision #: $Rev: 150 $
! URL: $HeadURL: http://sel1004.verit.dnv.com:8080/svn/LoadSimCtl_SurfaceIce/trunk/IceDyn_IntelFortran/HAWC2_DLL/HAWC2_DLL.f90 $
!**********************************************************************************************************************************

! DLL to be linked with the HAWC2 Aeroelastic wind turbine simulation
! As a type 2 DLL per the user manual
module HAWC2DLL

   use IceLog
   use IceInputParams
   use iceFloeBase
   use IceCrushingISO
   use IceCrushingIEC
   use IceIntermittentCrushing
   use IceFlexISO
   use IceFlexIEC
   use IceLockInCrushingISO
   use randomCrushing
   use IceCpldCrushing
   use NWTC_IO, only : DispNVD, progdesc, curdate, curtime

   TYPE(iceFloe_ParameterType), save  :: icep           ! Parameters, including precalculated time series
   type(iceFloe_LoggingType),save     :: theIceLog   ! structure with message and error logging variables 

end module HAWC2DLL


subroutine icefloe_init (array1, array2)

! Expose subroutine DLL to users of this DLL
!DEC$ ATTRIBUTES DLLEXPORT, C, ALIAS:'icefloe_init' :: icefloe_init
   
   use HAWC2DLL

   implicit none
   
   ! input
   real(8), intent(IN)        :: array1(1)
   ! output
   real(8), intent(OUT)       :: array2(1)

   ! locals
   Real(ReKi)                 :: duration
   INTEGER(IntKi)             :: nSteps
   INTEGER(IntKi)             :: Err     ! Error status of the operation
   LOGICAL, SAVE              :: bInit = .FALSE.      ! Initialization flag
   TYPE(ProgDesc), PARAMETER  :: IceFloe_Ver = ProgDesc( 'IceFloe', 'v1.00.00', 'May-2014' )

! More IceFloe types
   type(iceInputType)         :: iceInput    ! hold list of input names and values from file
   character(132)             :: logFile  ! File name for message logging
   
   ! Initialise on first call
   IF (.NOT.bInit) THEN
      bInit = .TRUE.
   ENDIF

!   Set up error logging
   theIceLog%warnFlag = .false.
   theIceLog%ErrID    = ErrID_None
   theIceLog%ErrMsg   = ""

   ! Initialize the NWTC Subroutine Library

!   CALL NWTC_Init( )

   ! Display the module information

   CALL DispNVD( IceFloe_Ver )

   call openIceLog (theIceLog, 'IceFloe.log')
   call logMessage(theIceLog, ' Running: '//trim(IceFloe_Ver%Name)//trim(IceFloe_Ver%Ver)//trim(IceFloe_Ver%Date))
   call logMessage(theIceLog, ' This run started on: '//curdate()//' at '//curtime()//newLine)

   call countIceInputs('testIce.inp', theIceLog, iceInput)
   call readIceInputs(theIceLog, iceInput)
   call logMessage(theIceLog, ' Input file read complete.'//newLine)

   ! call IceFloe initialization routine and get parameters back
   ! not all parameters used by all ice floe types
   call getIceInput(iceInput, 'iceType',icep%iceType, theIceLog, lowTypeLimit, hiTypeLimit)
   if (theIceLog%ErrID >= AbortErrLev) then
      return   
   endif

!  Set the time step as the minimum of the suggested p%dt and the time step from the ice input file
   call getIceInput(iceInput, 'timeStep',icep%dt, theIceLog, 0.0)
   if (theIceLog%ErrID >= AbortErrLev) then
      return   
   endif

   ! get the duration of the simulation
      call getIceInput(iceInput, 'duration', duration, theIceLog, 0.0)
      call logMessage(theIceLog, ' Load time series duration = '//TRIM(Num2LStr(duration))//' sec')

   ! get the load ramp up time
      call getIceInput(iceInput, 'rampTime', icep%rampTime, theIceLog, 0.0)
      call logMessage(theIceLog, ' Load ramp up time = '//TRIM(Num2LStr(icep%rampTime))//' sec')

   ! get the number of legs on the support structure
      call getIceInput(iceInput, 'numLegs', icep%numLegs, theIceLog, 1, 4)
      if (theIceLog%ErrID >= AbortErrLev) then
         return   
      endif

   ! allocate storage for load series
      nSteps = ceiling(duration/icep%dt) + 1  ! + 1 for zero point
      allocate(icep%loadSeries(nSteps, icep%numLegs), stat=err)
      if (err /= 0) then
         call iceErrorHndlr (theIceLog, ErrID_Fatal, 'Error in time series array allocation in IceFloe_init', 1)
         return
      endif   
   ! point internal iceFloe array to the saved load series
   !   icep%loadSeries => loadseries

   ! allocate storage for the leg positions
      allocate (icep%legX(icep%numLegs), stat=err)
      allocate (icep%legY(icep%numLegs), stat=err)
      allocate (icep%ks(icep%numLegs), stat=err)
      if (err /= 0) then
         call iceErrorHndlr (theIceLog, ErrID_Fatal, 'Error in allocation of leg data in parameters', 1)
         return  !  error in array allocation
      endif   
      icep%legX = 0.0
      icep%legY = 0.0
      icep%ks   = 1.0

   iceType: select case (icep%iceType)
     case (randomCrush)
        call initRandomCrushing(iceInput, icep, theIceLog)
        if (theIceLog%ErrID <= ErrID_Warn)  &
           call logMessage(theIceLog, newLine//' Random continuous ice crushing via ISO/Karna initialized'//newLine)
     case (interCrush)
        call initInterCrushing(iceInput, icep, theIceLog)
        if (theIceLog%ErrID <= ErrID_Warn)  &
           call logMessage(theIceLog, newLine//' Intermittent ice crushing loads initialized'//newLine)
     case (crushIEC)
        call initCrushingIEC(iceInput, icep, theIceLog)
        if (theIceLog%ErrID <= ErrID_Warn)  &
           call logMessage(theIceLog, newLine//' Ice crushing loads IEC/Korzhavin initialized'//newLine)
     case (flexFailISO)
        call initFlexISO(iceInput, icep, theIceLog)
        if (theIceLog%ErrID <= ErrID_Warn)  &
           call logMessage(theIceLog, newLine//' ISO/Croasdale ice flexural failure loads initialized'//newLine)
     case (flexFailIEC)
        call initFlexIEC(iceInput, icep, theIceLog)
        if (theIceLog%ErrID <= ErrID_Warn)  &
           call logMessage(theIceLog, newLine//' IEC/Ralston ice flexural failure loads initialized'//newLine)
     case (lockInISO)
        call initLockInCrushingISO(iceInput, icep, theIceLog)
        if (theIceLog%ErrID <= ErrID_Warn)  &
           call logMessage(theIceLog, newLine//' Frequency lock-in ice crushing loads via ISO initialized'//newLine)
     case (cpldCrush)
        call initCpldCrushing(iceInput, icep, theIceLog)
        if (theIceLog%ErrID <= ErrID_Warn)  &
           call logMessage(theIceLog, newLine//' Coupled crushing ice loads (Maattanen) initialized'//newLine)
     case default
        call iceErrorHndlr (theIceLog, ErrID_Fatal, 'Invalid ice floe type, '                &
                           //TRIM(Num2LStr(icep%IceType))//' is not a valid selection', 1)
   end select iceType

   array2 = 0.0d0

end subroutine icefloe_init


SUBROUTINE icefloe_update(array1, array2)

! Expose subroutine DLL to users of this DLL
!DEC$ ATTRIBUTES DLLEXPORT, C, ALIAS:'icefloe_update' :: icefloe_update

   use HAWC2DLL
   
   implicit none
   ! input
   real(8), intent(IN)     :: array1(1+3*icep%numLegs)  ! 1) Time, 2) Vx, 3) Vy, 4) Vz  all in HAWC coords
   ! output
   real(8), intent(OUT)    :: array2(2*icep%numLegs)  ! Fx and Fy (in HAWC coords)

   REAL(DbKi)        :: t           ! Current simulation time in seconds
   real(ReKi)        :: loadVect(6,icep%numLegs)
   real(ReKi)        :: inVels(2,icep%numLegs)
   integer(IntKi)    :: nL

!  reset up error logging
   theIceLog%warnFlag = .false.
   theIceLog%ErrID    = ErrID_None
   theIceLog%ErrMsg   = ""

   t = array1(1)
   inVels = 0.0
! HAWC has +y axis downwind, x lateral, z + down
! IceFloe has +x axis downwind, y lateral, z + up
   do nL = 1, icep%numLegs
      inVels(1,nL) = array1(3*nL)
      inVels(2,nL) = array1(3*nL-1)
   enddo

! get loads from IceFloe
   iceType: select case (icep%iceType)
     case (randomCrush)
        loadVect = outputRandomCrushLoad(icep, theIceLog, t)
     case (crushIEC)
        loadVect = outputCrushLoadIEC(icep, theIceLog, t)
     case (interCrush)
        loadVect = outputInterCrushLoad(icep, theIceLog, t)
     case (flexFailISO)
        loadVect = outputFlexLoadISO(icep, theIceLog, t)
     case (flexFailIEC)
        loadVect = outputFlexLoadIEC(icep, theIceLog, t)
     case (lockInISO)
        loadVect = outputLockInLoadISO(icep, theIceLog, t)
     case (cpldCrush)
        loadVect = outputCpldCrushLoad(icep, theIceLog, inVels, t)
     case default
        call logFatal (theIceLog, 'Invalid Ice Floe Type Selection', 1)
   end select iceType

!TODO  deal w/ single point application
!Apply ramp for first 10 seconds
   loadVect = loadVect*min(1.0, array1(1)/icep%rampTime)
   do nL = 1, icep%numLegs
      array2(2*nL-1) = dble(loadVect(2,nL))*min(1.0, 0.1*t)
      array2(2*nL)   = dble(loadVect(1,nL))*min(1.0, 0.1*t)
   enddo

end SUBROUTINE icefloe_update
