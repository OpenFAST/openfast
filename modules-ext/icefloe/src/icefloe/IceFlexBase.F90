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
! Routines for ice floes that apply loads to the structure
! via flexural failure of the ice sheet
! This module is a base module for flexural failure where a 
! time series is precalculated

module IceFlexBase

   use IceFloeBase

   implicit none

   public

!--------------------------------------------------
contains

!  calculate the time series of ice loads using the Croasdale method
   subroutine initIceFlex (iceInput, inParams, myIceParams, iceLog)

      type(iceInputType), intent(in)            :: iceInput
      type(iceFloe_ParameterType), intent(inout)    :: myIceParams
      type(iceFloe_LoggingType), intent(inout)  :: iceLog   ! structure with message and error logging variables
      type(inputParams), intent(out)            :: inParams ! specific input parameter variable list
      real(ReKi)                                :: frictionLimit

!  initialize the common parmeters
      call initIceFloe(iceInput, inParams, myIceParams, iceLog)

      call logMessage(iceLog, newLine//' Setting common flexural failure input parameters ')

      call getIceInput(iceInput, 'flexStrength', inParams%flexStrength, iceLog, 0.0_ReKi, 1.0E9_ReKi)
      call logMessage(iceLog, ' flexStrength = '//TRIM(Num2LStr(inParams%flexStrength))//' Pascals')

      call getIceInput(iceInput, 'towerConeAngle', inParams%twr%coneAngle, iceLog, 20.0_ReKi, 70.0_ReKi)
      call logMessage(iceLog, ' towerConeAngle = '//TRIM(Num2LStr(inParams%twr%coneAngle))//' degrees')
      inParams%twr%coneAngle = D2R*inParams%twr%coneAngle   ! convert to radians

!   Check on friction bounds to insure against incorrect zero or negative calculation of static load
      if (inParams%twr%coneAngle > 45) then
         frictionLimit = cos(inParams%twr%coneAngle)/sin(inParams%twr%coneAngle) - 0.01
      else
         frictionLimit = sin(inParams%twr%coneAngle)/cos(inParams%twr%coneAngle) - 0.01
      endif
      call getIceInput(iceInput, 'ice2twrFriction', inParams%ice2twrFriction, iceLog, 0.0_ReKi, frictionLimit)
      call logMessage(iceLog, ' ice2twrFriction = '//TRIM(Num2LStr(inParams%ice2twrFriction)))

      call getIceInput(iceInput, 'iceDensity', inParams%iceDensity, iceLog, 0.0_ReKi)
      call logMessage(iceLog, ' iceDensity = '//TRIM(Num2LStr(inParams%iceDensity))//' kg/m^3')

      call getIceInput(iceInput, 'includeHb', inParams%includeHb, iceLog)
      if (inParams%includeHb) call logMessage(iceLog, ' Breaking term, Hb term is included')
      if (.not. inParams%includeHb) call logMessage(iceLog, ' Breaking term, Hb term is NOT included')

      call getIceInput(iceInput, 'includeHr', inParams%includeHr, iceLog)
      if (inParams%includeHr) call logMessage(iceLog, ' Rubble pushing term, Hr term is included')
      if (.not. inParams%includeHr) call logMessage(iceLog, ' Rubble pushing term, term is NOT included')
      
   end subroutine initIceFlex

!****************************************************************************
   function outputFlexLoad (myIceParams, iceLog, time)  result(iceLoads)
      type(iceFloe_ParameterType), intent(in)      :: myIceParams
      type(iceFloe_LoggingType), intent(inout) :: iceLog   ! structure with message and error logging variables
      real(DbKi), intent(in)                  :: time
      real(ReKi)                              :: iceLoads(6,myIceParams%numLegs)

      iceLoads = outputIceLoads (myIceParams, iceLog, time)
   end function outputFlexLoad

end module IceFlexBase
