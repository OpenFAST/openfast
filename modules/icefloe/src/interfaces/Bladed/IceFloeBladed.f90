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
! File last committed: $Date: 2014-01-07 14:11:45 -0800 (Tue, 07 Jan 2014) $
! (File) Revision #: $Rev: 120 $
! URL: $HeadURL: http://sel1004.verit.dnv.com:8080/svn/LoadSimCtl_SurfaceIce/trunk/IceDyn_IntelFortran/IceFloeConsole/source/IceFloe.f90 $
!**********************************************************************************************************************************

!****************************************************************
!  Console program for IceFloe package specific to the Bladed Point load format
!  Used to generate time series of ice loads to be read in by Bladed

program IceFloeBladed

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
   use NWTC_IO, only : DispNVD, progdesc, curdate, curtime, nameofile

   implicit none

   TYPE(iceFloe_ParameterType), save   :: icep           ! Parameters, including precalculated time series
   type(iceFloe_LoggingType),save      :: theIceLog      ! structure with message and error logging variables 
   type(iceInputType)                  :: iceInput       ! hold list of input names and values from file

   Real(ReKi)                 :: duration
   INTEGER(IntKi)             :: nSteps
   INTEGER(IntKi)             :: Err                  ! Error status of the operation
   LOGICAL, SAVE              :: bInit = .FALSE.      ! Initialization flag
   TYPE(ProgDesc), PARAMETER  :: IceFloe_Ver = ProgDesc( 'IceFloe', 'v1.00.00', 'May-2014' )

   CHARACTER(msgLen)          :: Msg      ! Error message 
   character(132)             :: logFile  ! File name for message logging
   character(132)             :: inputFile
   character(132)             :: outFile
   character(4)               :: fileLegNum
   integer(IntKi)             :: outUnitNum, numArg, n, nL

! get parameter input file from command line
   numArg = COMMAND_ARGUMENT_COUNT()
   if (numArg > 0) then
      CALL GET_COMMAND_ARGUMENT( 1, inputFile )
   else
      write(*,*) ' No input file specified on command line.  Exiting program.'
      stop
   endif

! Set up error logging
   theIceLog%warnFlag = .false.
   theIceLog%ErrID    = ErrID_None
   theIceLog%ErrMsg   = ""

! Display the module information

   CALL DispNVD( IceFloe_Ver )

   call NameOFile ( 1, 'log', logFile )   ! The filename comes from the command line
   call openIceLog (theIceLog, logFile)
   call logMessage(theIceLog, ' Running: '//trim(IceFloe_Ver%Name)//trim(IceFloe_Ver%Ver)//trim(IceFloe_Ver%Date))
   call logMessage(theIceLog, ' This run started on: '//curdate()//' at '//curtime()//newLine)

! read the input file twice: once to count the other to read
   call countIceInputs(inputFile, theIceLog, iceInput)
   call readIceInputs(theIceLog, iceInput)
   call logMessage(theIceLog, ' Input file read complete.'//newLine)

! get some top level input parameters
! the type of ice loading
   call getIceInput(iceInput, 'iceType',icep%iceType, theIceLog, lowTypeLimit, hiTypeLimit)
   if (theIceLog%ErrID >= AbortErrLev) then
      stop   
   endif
! the time step
   call getIceInput(iceInput, 'timeStep',icep%dt, theIceLog, 0.0)
   if (theIceLog%ErrID >= AbortErrLev) then
      stop   
   endif
   call getIceInput(iceInput, 'rampTime',icep%rampTime, theIceLog, 0.0)
   if (theIceLog%ErrID >= AbortErrLev) then
      stop   
   endif
! the duration of the simulation
   call getIceInput(iceInput, 'duration', duration, theIceLog, 0.0)
   call logMessage(theIceLog, ' Load time series duration = '//TRIM(Num2LStr(duration))//' sec')
! the number of legs on the support structure
   call getIceInput(iceInput, 'numLegs', icep%numLegs, theIceLog, 1, 4)
   if (theIceLog%ErrID >= AbortErrLev) then
      stop  
   endif

! allocate storage for load time series
   nSteps = ceiling(duration/icep%dt) + 1  ! + 1 for zero point
   allocate(icep%loadSeries(nSteps, icep%numLegs), stat=err)
   if (err /= 0) then
      call iceErrorHndlr (theIceLog, ErrID_Fatal, 'Error in time series array allocation in IceFloe_init', 1)
      stop
   endif   

! allocate storage for the leg positions
   allocate (icep%legX(icep%numLegs), stat=err)
   allocate (icep%legY(icep%numLegs), stat=err)
   allocate (icep%ks(icep%numLegs), stat=err)
   if (err /= 0) then
      call iceErrorHndlr (theIceLog, ErrID_Fatal, 'Error in allocation of leg data in parameters', 1)
      stop  !  error in array allocation
   endif   
   icep%legX = 0.0
   icep%legY = 0.0
   icep%ks   = 1.0

! initialize and calculate load time series
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
     case default
        call iceErrorHndlr (theIceLog, ErrID_Fatal, 'Invalid ice floe type, '                &
                           //TRIM(Num2LStr(icep%IceType))//' is not a valid selection', 1)
   end select iceType

   if (theIceLog%ErrID >= AbortErrLev) stop

!  Ouput time series of loads to ASCII file
!  If only one leg or applying as a single load for multiple legs
!  then only writing one file
   write(*,*) 'Writing load series to file.'
   if (icep%singleLoad .or. icep%numLegs == 1) then
      call NameOFile ( 1, 'dat', outFile )   ! The filename comes from the command line
      call GetNewUnit ( outUnitNum, Err, Msg )
      call OpenFOutFile ( outUnitNum, outFile, Err, Msg )
      if (Err >= ErrID_Severe) then
         call iceErrorHndlr (theIceLog, ErrID_Fatal, 'Error opening output file: '//newLine//trim(Msg), 1)
         stop
      endif
      call writeTimeSeries (outUnitNum, nSteps, 1, icep)
      close (outUnitNum)
   else
!  Need one file for each load point, i.e. each leg
      do nL = 1, icep%numLegs   
         write(fileLegNum,'("dat",I1)') nL
         call NameOFile ( 1, fileLegNum, outFile )   ! The filename comes from the command line
         call GetNewUnit ( outUnitNum, Err, Msg )
         call OpenFOutFile ( outUnitNum, outFile, Err, Msg )
         call writeTimeSeries (outUnitNum, nSteps, nL, icep)
         close (outUnitNum)
      enddo
   endif

   write(*,*) 'Load series write complete.'

end program IceFloeBladed


subroutine writeTimeSeries (outUnitNum, nSteps, nL, icep)

   use iceFloebase
   implicit none

   TYPE(iceFloe_ParameterType), intent(IN) :: icep           ! Parameters, including precalculated time series
   integer(IntKi), intent(IN)              :: nSteps, outUnitNum, nL

   real(ReKi)                          :: loadVect(6,icep%numLegs)        ! Forces and moments for a given leg at a given time step
   integer(IntKi)                      :: n, nn

!  Start with header info and zeros for the first time step as suggested in Bladed User Manual
   write(outUnitNum,'("NUMPOINTS  ", I10)') nSteps
   write(outUnitNum,'(7(1pe12.5,1x))') (0.0, nn=1,7)
   do n = 1, nSteps-2
   !  If the loads from multiple legs are applied as an "effective" load on one leg at the
   !  support center, then compute that effective load vector and put into the first load vector
   !  Note that only Fx, Fy, and Mz terms re-calculated.  All others remain zero
      if (icep%numLegs > 1 .and. icep%singleLoad) then
         do nn = 1, icep%numLegs
            loadVect(:,nn) = iceLoadDirection(icep%loadSeries(n+1,nn), icep)
         enddo
         loadVect(:,1) = iceLoadEquivalent(loadVect, icep)
      else
         loadVect(:,1) = iceLoadDirection(icep%loadSeries(n+1,nL), icep)
      endif      
      loadVect = loadVect*min(1.0, real(n,ReKi)/(icep%rampTime/icep%dt))  ! apply ramp as loads can be large impacts
      write(outUnitNum,'(7(1pe12.5,1x))') real(n,ReKi)*icep%dt, (loadVect(nn,1), nn = 1,6)
   enddo

!  end with zeros for the last time step
   write(outUnitNum,'(7(1pe12.5,1x))') real(nSteps-1,ReKi)*icep%dt, (0.0, nn=1,6)

end subroutine writeTimeSeries
!**********************************************************************************************************************************
