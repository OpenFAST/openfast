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
! Module to handle parameter input from user's file

module iceInputParams

   use , INTRINSIC :: iso_fortran_env , only : iostat_end  ! for end of file detection
   use precision
   use NWTC_IO, only : GetNewUnit, OpenFInpFile
   use iceLog

   implicit none

!  This is used for initial read in of variables by name
!  Extraction of values into specific named variables done by getIceInput functions
   type paramType
      character(132) :: name
      real(ReKi)     :: value
   end type paramType

!  Hold some meta data for the input parameters
   type iceInputType
      integer(IntKi) :: unitNum
      integer(IntKi) :: count
      type(paramType), allocatable :: params(:)
   end type iceInputType

!  generic functions for retrieving input parameters
   interface getIceInput
      module procedure getRealInput
      module procedure getIntInput
      module procedure getLogicalInput
   end interface getIceInput

!  data structures for tower - much of this could come from FAST/ElastoDyn/SubDyn 
!  but need for other codes
   type twrLegs
      real(ReKi)     :: X
      real(ReKi)     :: Y
      real(ReKi)     :: phase
      real(ReKi)     :: ks
   end type twrLegs
   type towerData
      real(ReKi)     :: diam           ! at water line - this could come from FAST but not other codes
      type(twrLegs), allocatable :: leg(:)
      real(ReKi)     :: freq           ! leg or tower fundamental frequency, Hz
      real(ReKi)     :: coneAngle      ! for flexural failure
      real(ReKi)     :: coneTopDiam    ! for flexural failure
   end type towerData

!  This is the list of input parameters
   type inputParams
      integer(IntKi) :: iceType
      real(ReKi)     :: iceThickness
      real(ReKi)     :: iceVelocity
      real(ReKi)     :: iceDirection
      real(ReKi)     :: timeStep
      real(ReKi)     :: rampTime
      integer(IntKi) :: randomSeed

      type(towerData):: twr
      real(ReKi)     :: multiLegFactor_kn ! Applied to peak load lock-in for all legs when numLegs>1
      logical        :: singleLoad        ! true if multiple leg loads will be combined into an equivalent load
      logical        :: legAutoFactor     ! code decides on leg load shelter factors if true

      real(ReKi)     :: refIceThick
      real(ReKi)     :: refIceStrength
      real(ReKi)     :: staticExponent

      real(ReKi)     :: shapeFactor_k1
      real(ReKi)     :: contactFactor_k2

      real(ReKi)     :: coeffPSD_b
      real(ReKi)     :: coeffPSD_ks
      real(ReKi)     :: crushLoadCOV
      real(ReKi)     :: stdLoadMult
      real(ReKi)     :: freqStep

      real(ReKi)     :: interPeriod
      real(ReKi)     :: riseTime
      real(ReKi)     :: fallTime

      real(ReKi)     :: minLoadFraction

      real(ReKi)     :: flexStrength
      real(ReKi)     :: ice2twrFriction
      real(ReKi)     :: iceDensity
      logical        :: includeHb
      logical        :: includeHr

      real(ReKi)     :: peakLoadCOV
      real(ReKi)     :: coeffLoadPeaks
      real(ReKi)     :: coeffLoadMin
      real(ReKi)     :: periodCOV
      real(ReKi)     :: tauMin
      real(ReKi)     :: tauMax
      real(ReKi)     :: coeffBreakLength
      real(ReKi)     :: rubbleHeight
      real(ReKi)     :: rubblePorosity
      real(ReKi)     :: rubbleCohesion
      real(ReKi)     :: rubbleAngle
      real(ReKi)     :: frictionAngle
      real(ReKi)     :: ice2iceFriction
      real(ReKi)     :: waterDensity
      real(ReKi)     :: iceModulus
      real(ReKi)     :: poissonRatio
      logical        :: includeHp
      logical        :: includeHl
      logical        :: includeHt
      logical        :: includeLc

      real(ReKi)     :: rideUpThickness
      real(ReKi)     :: freqParamK

      real(ReKi)     :: minStrength
      real(ReKi)     :: minStrengthNegVel 
    end type inputParams

contains

!========================================================================
   subroutine countIceInputs (fname, iceLog, input)
      character(*), intent(in)                  :: fname    ! name of file with input parameters
      type(iceFloe_LoggingType), intent(inout)   :: iceLog   ! structure with message and error logging variables
      type(iceInputType), intent(out)              :: input

      character(1024)   :: lineIn
      INTEGER(IntKi)    :: ErrStat     ! Error status of the operation
      CHARACTER(msgLen) :: ErrMsg      ! Error message if ErrStat /= ErrID_None
      integer(IntKi)    :: ioStatus

      call GetNewUnit ( input%UnitNum, ErrStat, ErrMsg )
      if (ErrStat >= AbortErrLev) then
         call iceErrorHndlr (iceLog, ErrStat, 'CountIceInputs => GetNewUnit: '//newLine//ErrMsg, 1)
         return
      endif

      call OpenFInpFile ( input%UnitNum, trim(fname), ErrStat, ErrMsg  )
      if (ErrStat >= AbortErrLev) then
         call iceErrorHndlr (iceLog, ErrStat, 'CountIceInputs => OpenFInpFile: '//newLine//ErrMsg, 1)
         return
      endif

      input%count = 0
      readLoop: do while (.true.)
         read(input%UnitNum,'(A)',iostat=ioStatus) lineIn
         if (ioStatus == iostat_end) exit readLoop
         if(lineIn(1:1) == '!' .or.  &   ! check for a comment
            lineIn(1:1) == '#' .or.  &
            lineIn(1:1) == '$' .or.  &
            lineIn(1:1) == '%' ) then
            cycle readLoop
         else
            input%count = input%count + 1
         endif
      enddo readLoop

      allocate(input%params(input%count), stat=errStat)
      if (errStat /= 0) then
         call iceErrorHndlr (iceLog, ErrID_Fatal, 'Error in input parameter array allocation in CountIceInputs', 1)
         return  !  error in array allocation
      endif   

!   Go back to the beginning of the file
      rewind (unit=input%UnitNum, iostat = ErrStat)
      if (ErrStat /= 0) then
         call iceErrorHndlr (iceLog, ErrID_Fatal, 'CountIceInputs => OpenFInpFile: '//newLine//    &
                            ' Error rewinding parameter input file', 1)
         return
      endif

   end subroutine countIceInputs
   
!========================================================================
   subroutine readIceInputs (iceLog, input)
      type(iceFloe_LoggingType), intent(inout)   :: iceLog   ! structure with message and error logging variables
      type(iceInputType), intent(inout)            :: input

      character(1024)   :: lineIn
      character(132)    :: varName
      real(ReKi)        :: varValue
      integer(IntKi)    :: ioStatus, nCount, n
      INTEGER(IntKi)    :: ErrStat     ! Error status of the operation
      CHARACTER(msgLen) :: ErrMsg      ! Error message if ErrStat /= ErrID_None


!  read all of the records and put into the proper variable
      nCount = 0
      readLoop: do while (nCount < input%count)
         read(input%UnitNum,'(A)',iostat=ioStatus) lineIn
         if (ioStatus == iostat_end) exit readLoop
         if(lineIn(1:1) == '!' .or.  &   ! check for a comment
            lineIn(1:1) == '$' .or.  &
            lineIn(1:1) == '%' ) cycle readLoop

!   read in the variable name and its value, convert name to upper case
         read(lineIn,*) varName, varValue
         call Conv2UC (varName)
         nCount = nCount + 1
         input%params(nCount)%name = varName
         input%params(nCount)%value = varValue
         
!   loop through the list of possible input variables check for duplication
         searchLoop: do n = 1, nCount-1
            if(index(varname, input%params(n)%name) > 0) then
               call iceErrorHndlr (iceLog, ErrID_Warn, 'Input parameter '//trim(varName)//  &
                                                       ' has been specified twice.', 0)
               exit searchLoop
            endif
         enddo searchloop
      enddo readLoop

      close (input%UnitNum)

   end subroutine readIceInputs

!========================================================================
!  extract the input parameter value from the list
!  and send back to calling routine
   subroutine getRealInput (input, varName, outVal, iceLog, min, max)
      type(iceInputType), intent(in)             :: input
      type(iceFloe_LoggingType), intent(inout) :: iceLog   ! structure with message and error logging variables
      character(*), intent(in)         :: varName     ! This is the variable requested
      real(ReKi), intent(out)          :: outVal
      real(ReKi), intent(in), optional :: min, max

      integer(IntKi)                :: n
      character(len_trim(varName))  :: tmpName
      logical                       :: foundParam

      tmpName = varName  ! can't change varName since it's an intent(in)
      call Conv2UC (tmpName)

      foundParam = .false.
      do n = 1, input%count
         if(index(input%params(n)%name, tmpName) > 0) then
            outVal = input%params(n)%value
            foundParam = .true.
            exit
         endif
      enddo

      if (.not. foundParam)  then
         call iceErrorHndlr (iceLog, ErrID_Fatal, 'Input parameter '//varName//' was not in the parameter input file.', 1 )
         return
      endif

!  If requeseted, check that the value is withing reasonable bounds
      if (present(min)) then
         if (outVal < min) then
            call iceErrorHndlr (iceLog, ErrID_Warn, 'Input parameter '//varName//' is below the minimum bound', 0 )
         endif
      endif
      if (present(max)) then
         if (outVal > max) then
            call iceErrorHndlr (iceLog, ErrID_Warn, 'Input parameter '//varName//' is above the maximum bound', 0 )
         endif
      endif

   end subroutine getRealInput

!------------------------------------------------------------------
   subroutine getIntInput (input, varName, outVal, iceLog, min, max)
      type(iceInputType), intent(in)             :: input
      type(iceFloe_LoggingType), intent(inout) :: iceLog   ! structure with message and error logging variables
      character(*), intent(in)                :: varName     ! This is the variable requested
      integer(IntKi), intent(out)             :: outVal
      integer(IntKi), intent(in), optional    :: min, max

      integer(IntKi)                :: n
      character(len_trim(varName))  :: tmpInName
      logical                       :: foundParam

      tmpInName = varName  ! can't use varname since Conv2UC changes it
      call Conv2UC (tmpInName)

      foundParam = .false.
      do n = 1, input%count
         if(index(input%params(n)%name, tmpInName) > 0) then
            outVal = input%params(n)%value
            foundParam = .true.
            exit
         endif
      enddo

      if (.not. foundParam)  then
         call iceErrorHndlr (iceLog, ErrID_Fatal, 'Input parameter '//varName//' was not in the parameter input file.', 1 )
         return
      endif

!  If requeseted, check that the value is withing reasonable bounds
      if (present(min)) then
         if (outVal < min) then
            call iceErrorHndlr (iceLog, ErrID_Warn, 'Input parameter '//varName//' is below the minimum bound', 0 )
         endif
      endif
      if (present(max)) then
         if (outVal > max) then
            call iceErrorHndlr (iceLog, ErrID_Warn, 'Input parameter '//varName//' is above the maximum bound', 0 )
         endif
      endif
      
   end subroutine getIntInput

!------------------------------------------------------------------
   subroutine getLogicalInput (input, varName, outVal, iceLog)
      type(iceInputType), intent(in)             :: input
      type(iceFloe_LoggingType), intent(inout)   :: iceLog   ! structure with message and error logging variables
      character(*), intent(in)      :: varName     ! This is the variable requested
      logical, intent(out)          :: outVal
      character(len_trim(varName))  :: tmpName
      integer(IntKi)                :: n
      logical                       :: foundParam

      tmpName = varName
      call Conv2UC (tmpName)

      foundParam = .false.
      do n = 1, input%count
         if(index(input%params(n)%name, tmpName) > 0) then
            !outVal = input%params(n)%value ! warning #6192: Fortran 2003 does not allow this data type conversion. <- real to logical
            outVal = ( NINT(input%params(n)%value) /= 0 )  !Intel: The numeric value of .TRUE. and .FALSE. can be -1 and 0 or 1 and 0 depending on compiler option fpscomp 
            foundParam = .true.
            exit
         endif
      enddo

      if (.not. foundParam)  then
         call iceErrorHndlr (iceLog, ErrID_Fatal, 'Input parameter '//varName//' was not in the parameter input file.', 1 )
         return
      endif
   end subroutine getLogicalInput

end module iceInputParams
