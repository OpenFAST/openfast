!**********************************************************************************************************************************
!
!  MODULE: ADsk_Driver_Types  - This module contains types used by the AeroDisk Driver program to store arguments passed in
!
!  The types listed here are used within the AeroDisk Driver program to store the settings. These settings are read in as
!  command line arguments, then stored within these types.
!
!**********************************************************************************************************************************
!
!..................................................................................................................................
! LICENSING
! Copyright (C) 2024  National Renewable Energy Laboratory
!
!    This file is part of AeroDisk.
!
!    AeroDisk is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
!    published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License along with AeroDisk.
!    If not, see <http://www.gnu.org/licenses/>.
!
!**********************************************************************************************************************************
module AeroDisk_Driver_Types

   use NWTC_Library
   use AeroDisk_Types

   implicit none

      !> This contains flags to note if the settings were made.  This same data structure is
      !! used both during the driver input file and the command line options.
      !!
      !! NOTE: The WindFileType is only set if it is given as a command line option.  Otherwise
      !!       it is handled internally by InflowWInd.
      !!
      !! NOTE: The wind direction is specified by the AeroDisk input file.
   type     :: ADskDriver_Flags
      logical                 :: DvrIptFile           = .FALSE.      !< Was an input file name given on the command line?
      logical                 :: ADskIptFile          = .FALSE.      !< Was an AeroDisk input file requested?
      logical                 :: OutRootName          = .FALSE.      !< Was an AeroDisk output rootname 
      logical                 :: AirDens              = .FALSE.      !< Air density
      logical                 :: RotorRad             = .FALSE.      !< rotor radius
      logical                 :: RotorHeight          = .FALSE.      !< rotor height
      logical                 :: ShftTilt             = .FALSE.      !< shaft tilt
      logical                 :: TStart               = .FALSE.      !< specified a start time
      logical                 :: NumTimeSteps         = .FALSE.      !< specified a number of timesteps to process
      logical                 :: NumTimeStepsDefault  = .FALSE.      !< specified a 'DEFAULT' for number of timesteps to process
      logical                 :: DT                   = .FALSE.      !< specified a resolution in time
      logical                 :: DTDefault            = .FALSE.      !< specified a 'DEFAULT' for the time resolution
      logical                 :: Verbose              = .FALSE.      !< Verbose error reporting
      logical                 :: VVerbose             = .FALSE.      !< Very Verbose error reporting
   end type    ADskDriver_Flags


      ! This contains all the settings (possible passed in arguments).
   type     :: ADskDriver_Settings
      character(1024)         :: DvrIptFileName                !< Driver input file name
      character(1024)         :: ADskIptFileName               !< Filename of AeroDisk input file to read (if no driver input file)
      character(1024)         :: OutRootName                   !< Output root name

      real(ReKi)              :: AirDens                       !< Air density (kg/m^3)
      real(ReKi)              :: RotorRad                      !< rotor radius (m)
      real(ReKi)              :: RotorHeight                   !< rotor height (m)
      real(ReKi)              :: ShftTilt                      !< shaft tilt (deg)

      real(DbKi)              :: TStart                        !< Start time
      integer(IntKi)          :: NumTimeSteps                  !< Number of timesteps
      real(DbKi)              :: DT                            !< resolution of time

      type(ProgDesc)          :: ProgInfo                      !< Program info
      type(ProgDesc)          :: ADskProgInfo                  !< Program info for AeroDisk
   end type    ADskDriver_Settings


end module AeroDisk_Driver_Types
