!**********************************************************************************************************************************
!
!  MODULE: IfW_Driver_Types  - This module contains types used by the InflowWind Driver program to store arguments passed in
!
!  The types listed here are used within the InflowWind Driver program to store the settings. These settings are read in as
!  command line arguments, then stored within these types.
!
!**********************************************************************************************************************************
!
!..................................................................................................................................
! LICENSING
! Copyright (C) 2012  National Renewable Energy Laboratory
!
!    This file is part of InflowWind.
!
!    InflowWind is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
!    published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License along with InflowWind.
!    If not, see <http://www.gnu.org/licenses/>.
!
!**********************************************************************************************************************************

MODULE InflowWind_Driver_Types

   USE NWTC_Library
   USE InflowWind_Types

   IMPLICIT NONE

      !> This contains flags to note if the settings were made.  This same data structure is
      !! used both during the driver input file and the command line options.
      !!
      !! NOTE: The WindFileType is only set if it is given as a command line option.  Otherwise
      !!       it is handled internally by InflowWInd.
      !!
      !! NOTE: The wind direction is specified by the InflowWind input file.
   TYPE     :: IfWDriver_Flags
      LOGICAL                 :: DvrIptFile           = .FALSE.      !< Was an input file name given on the command line?
      LOGICAL                 :: IfWIptFile           = .FALSE.      !< Was an InflowWind input file requested?
      LOGICAL                 :: Summary              = .FALSE.      !< create a summary at command line? (data extents in the wind file)
      LOGICAL                 :: SummaryFile          = .FALSE.      !< create a summary file of the output?
      LOGICAL                 :: TStart               = .FALSE.      !< specified a start time
      LOGICAL                 :: NumTimeSteps         = .FALSE.      !< specified a number of timesteps to process
      LOGICAL                 :: NumTimeStepsDefault  = .FALSE.      !< specified a 'DEFAULT' for number of timesteps to process
      LOGICAL                 :: DT                   = .FALSE.      !< specified a resolution in time
      LOGICAL                 :: DTDefault            = .FALSE.      !< specified a 'DEFAULT' for the time resolution

      LOGICAL                 :: TurbineHeight        = .FALSE.      !< specified a turbine height
      LOGICAL                 :: Width                = .FALSE.      !< specified a width

      LOGICAL                 :: FFTcalc              = .FALSE.      !< do an FFT

      LOGICAL                 :: WindGrid             = .FALSE.      !< Requested output of wind data on a grid -- input file option only
      LOGICAL                 :: XRange               = .FALSE.      !< specified a range of x      -- command line option only -- stored as GridCtrCoord and GridDelta
      LOGICAL                 :: YRange               = .FALSE.      !< specified a range of y      -- command line option only -- stored as GridCtrCoord and GridDelta
      LOGICAL                 :: ZRange               = .FALSE.      !< specified a range of z      -- command line option only -- stored as GridCtrCoord and GridDelta
      LOGICAL                 :: Dx                   = .FALSE.      !< specified a resolution in x -- command line option only, 0.0 otherwise
      LOGICAL                 :: Dy                   = .FALSE.      !< speficied a resolution in y
      LOGICAL                 :: Dz                   = .FALSE.      !< specified a resolution in z

      LOGICAL                 :: PointsFile           = .FALSE.      !< points filename to read in  -- command line option only
   END TYPE    IfWDriver_Flags


      ! This contains all the settings (possible passed in arguments).
   TYPE     :: IfWDriver_Settings
      CHARACTER(1024)         :: DvrIptFileName                !< Driver input file name
      CHARACTER(1024)         :: IfWIptFileName                !< Filename of InflowWind input file to read (if no driver input file)
      CHARACTER(1024)         :: SummaryFileName               !< Filename for the summary information output
      INTEGER(IntKi)          :: NumTimeSteps                  !< Number of timesteps
      REAL( DbKi )            :: DT                            !< resolution of time
      REAL( DbKi )            :: TStart                        !< range of time -- end time converted from TRange (command line option only)

      REAL( ReKi )            :: TurbineHeight                 !< Height of the turbine
      REAL( ReKi )            :: Width                         !< Reference width

      REAL( ReKi )            :: FFTcoord(1:3)                 !< (x,y,z) coordinate to do an FFT at

!      REAL( ReKi )            :: GridCtrCoord(1:3)             !< (x,y,z) coordinate for center of the wind grid to output
      REAL( ReKi )            :: GridDelta(1:3)                !< (GridDx,GridDy,GridDz) -- grid point spacing
      INTEGER(IntKi)          :: GridN(1:3)                    !< (GridNx,GridNy,GridNz) -- number of grid points

      REAL( ReKi )            :: XRange(1:2)                   !< Range in the x-direction for the gridded data
      REAL( ReKi )            :: YRange(1:2)                   !< Range in the y-direction for the gridded data
      REAL( ReKi )            :: ZRange(1:2)                   !< Range in the z-direction for the gridded data

      CHARACTER(1024)         :: PointsFileName                !< Filename of points file to read in
   END TYPE    IfWDriver_Settings


END MODULE InflowWind_Driver_Types
