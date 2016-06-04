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
! File last committed: $Date: 2014-01-31 09:14:28 -0800 (Fri, 31 Jan 2014) $
! (File) Revision #: $Rev: 131 $
! URL: $HeadURL: http://sel1004.verit.dnv.com:8080/svn/LoadSimCtl_SurfaceIce/trunk/IceDyn_IntelFortran/HAWC2_DLL/HAWC2_DLL.f90 $
!**********************************************************************************************************************************

module IceADAMS

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

   type(iceFloe_ParameterType), save   :: icep          ! Parameters, including precalculated time series
   type(iceFloe_LoggingType), save     :: theIceLog     ! structure with message and error logging variables 

contains
!------------------------------------------------------------------------------------
   subroutine icefloe_init()

      implicit none

      ! locals
      Real(ReKi)                 :: duration
      INTEGER(IntKi)             :: nSteps
      INTEGER(IntKi)             :: Err     ! Error status of the operation
      LOGICAL, SAVE              :: bInit = .FALSE.      ! Initialization flag
      TYPE(ProgDesc), PARAMETER  :: IceFloe_Ver = ProgDesc( 'IceFloe', 'v1.00.00', 'May-2014' )

   ! More IceFloe types
      type(iceInputType)         :: iceInput    ! hold list of input names and values from file
      character(132)             :: logFile  ! File name for message logging
      
   !   Set up error logging
      theIceLog%warnFlag = .false.
      theIceLog%ErrID    = ErrID_None
      theIceLog%ErrMsg   = ""

   ! Display the module information

      CALL DispNVD( IceFloe_Ver )

      call openIceLog(theIceLog, 'IceFloe.log')
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

   end subroutine icefloe_init

!**********************************************************************************************************************************
   SUBROUTINE icefloe_update(Time, Vels, Forces)

      implicit none
      ! input
      real(DbKi), intent(IN)     :: Time
      real(ReKi), intent(IN)     :: Vels(3)
      ! output
      real(ReKi), intent(OUT)    :: Forces(3)

      real(ReKi)        :: loadVect(6,icep%numLegs)
      real(ReKi)        :: inVels(2,icep%numLegs)
      integer(IntKi)    :: nL

   !  reset up error logging
      theIceLog%warnFlag = .false.
      theIceLog%ErrID    = ErrID_None
      theIceLog%ErrMsg   = ""

   !  Assign support structure velocities to local variables
   !  CURRENTLY only monopile supported so NL should = 1
      inVels = 0.0
      do nL = 1, icep%numLegs
         inVels(1,nL) = Vels(1)
         inVels(2,nL) = Vels(2)
      enddo

   ! get loads from IceFloe
      iceType: select case (icep%iceType)
        case (randomCrush)
           loadVect = outputRandomCrushLoad(icep, theIceLog, time)
        case (crushIEC)
           loadVect = outputCrushLoadIEC(icep, theIceLog, time)
        case (interCrush)
           loadVect = outputInterCrushLoad(icep, theIceLog, time)
        case (flexFailISO)
           loadVect = outputFlexLoadISO(icep, theIceLog, time)
        case (flexFailIEC)
           loadVect = outputFlexLoadIEC(icep, theIceLog, time)
        case (lockInISO)
           loadVect = outputLockInLoadISO(icep, theIceLog, time)
        case (cpldCrush)
           loadVect = outputCpldCrushLoad(icep, theIceLog, inVels, time)
        case default
           call logFatal (theIceLog, 'Invalid Ice Floe Type Selection', 1)
      end select iceType

   !Apply ramp for first 10 seconds
      do nL = 1, icep%numLegs
         Forces(1) = loadVect(1,nL)*min(1.0, 0.1*time)
         Forces(2) = loadVect(2,nL)*min(1.0, 0.1*time)
         Forces(3) = 0.0
      enddo

   end SUBROUTINE icefloe_update

end module IceADAMS

!**********************************************************************************************************************************

SUBROUTINE VFOSUB ( ID, ATIME, PAR, NPAR, DFLAG, IFLAG, Forces )

! This routine is used to calculate the functional values for a VFORCE statement.
! Specifically to call the IceFloe routines for use in ADAMS

   use IceADAMS

   IMPLICIT    NONE

! Passed variables:

   INTEGER(4), INTENT(IN )       :: ID                                              ! ID of the REQUEST statement calling this routine.
   INTEGER(4), INTENT(IN )       :: NPAR                                            ! Number of PARameters passed through PAR.

   REAL(8),     INTENT(IN )      :: ATIME                                           ! The current simimulation time.
   REAL(8),     INTENT(IN )      :: PAR(NPAR)                                       ! The passeed PARameters.
   REAL(8),     INTENT(OUT)      :: Forces(3)                                       ! Values returned to ADAMS.

   LOGICAL, INTENT(IN)           :: DFLAG                                           ! Partial Derivative flag.
   LOGICAL, INTENT(IN)           :: IFLAG                                           ! Initialization pass flag.

! Local variables
   LOGICAL, save                 :: initFlag = .true.
   logical                       :: errFlg
   integer(4)                    :: IPar3(3), NStates
   REAL(4)                       :: IceForces(3)                                     
   REAL(8)                       :: States(3)                                       

!  get the local velocity at the marker specified in PAR
   IPar3(1) = int(PAR(1))
   IPar3(2) = 1   ! ground marker
   IPar3(3) = 1   ! ground marker
   call SYSARY('TVEL', IPar3, 3, States, NStates, errFlg )

   IceForces = 0.0
   if (IFLAG .and. initFlag) then
      call icefloe_init
      initFlag = .false.
   endif

   if (.not. initFlag) then
!  Get the ice loads
      call icefloe_update(atime, real(States, 4), IceForces)
   endif

   Forces = real(IceForces, 8)/1000.0

END SUBROUTINE VFOSUB

