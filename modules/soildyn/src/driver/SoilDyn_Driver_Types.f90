!**********************************************************************************************************************************
!
!  MODULE: SlD_Driver_Types  - This module contains types used by the SoilDyn Driver program to store arguments passed in
!
!  The types listed here are used within the SoilDyn Driver program to store the settings. These settings are read in as
!  command line arguments, then stored within these types.
!
!**********************************************************************************************************************************
!
!..................................................................................................................................
! LICENSING
! Copyright (C) 2015  National Renewable Energy Laboratory
!
!    This file is part of SoilDyn.
!
!    SoilDyn is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
!    published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License along with SoilDyn.
!    If not, see <http://www.gnu.org/licenses/>.
!
!**********************************************************************************************************************************

MODULE SoilDyn_Driver_Types

   USE NWTC_Library
   USE SoilDyn_Types

   IMPLICIT NONE

      !> This contains flags to note if the settings were made.  This same data structure is
      !! used both during the driver input file and the command line options.
   TYPE     :: SlDDriver_Flags
      LOGICAL                 :: DvrIptFile           = .FALSE.      !< Was an input file name given on the command line?
      LOGICAL                 :: SlDIptFile           = .FALSE.      !< Was an SoilDyn input file requested?
      LOGICAL                 :: InputDispFile        = .FALSE.      !< Input displacement time series
      LOGICAL                 :: TStart               = .FALSE.      !< specified a start time
      LOGICAL                 :: StiffMatOut          = .FALSE.      !< output stiffness matrices at start and finish
      LOGICAL                 :: NumTimeSteps         = .FALSE.      !< specified a number of timesteps to process
      LOGICAL                 :: NumTimeStepsDefault  = .FALSE.      !< specified a 'DEFAULT' for number of timesteps to process
      LOGICAL                 :: DT                   = .FALSE.      !< specified a resolution in time
      LOGICAL                 :: DTDefault            = .FALSE.      !< specified a 'DEFAULT' for the time resolution
      LOGICAL                 :: Verbose              = .FALSE.      !< Verbose error reporting
      LOGICAL                 :: VVerbose             = .FALSE.      !< Very Verbose error reporting
      LOGICAL                 :: SlDNonLinearForcePortionOnly = .FALSE. !< To only return the non-linear portion of the reaction force
   END TYPE    SlDDriver_Flags


      ! This contains all the settings (possible passed in arguments).
   TYPE     :: SlDDriver_Settings
      CHARACTER(1024)         :: DvrIptFileName                !< Driver input file name
      CHARACTER(1024)         :: SlDIptFileName                !< Filename of SoilDyn input file to read (if no driver input file)
      CHARACTER(1024)         :: InputDispFile                 !< Filename of SoilDyn time series displacements

      INTEGER(IntKi)          :: NumTimeSteps                  !< Number of timesteps
      REAL(DbKi)              :: DT                            !< resolution of time
      REAL(DbKi)              :: TStart                        !< Start time

      TYPE(ProgDesc)          :: ProgInfo                      !< Program info
      TYPE(ProgDesc)          :: SlDProgInfo                   !< Program info for SoilDyn

   END TYPE    SlDDriver_Settings


END MODULE SoilDyn_Driver_Types
