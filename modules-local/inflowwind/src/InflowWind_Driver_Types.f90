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

      ! This contains flags to note if the settings were made.
   TYPE     :: InflowWind_Driver_ArgFlags
      LOGICAL                 :: WindFileType   = .FALSE.      ! specified a windfiletype
      LOGICAL                 :: Height         = .FALSE.      ! specified a height
      LOGICAL                 :: Width          = .FALSE.      ! specified a width
      LOGICAL                 :: XRange         = .FALSE.      ! specified a range of x
      LOGICAL                 :: YRange         = .FALSE.      ! specified a range of y
      LOGICAL                 :: ZRange         = .FALSE.      ! specified a range of z
      LOGICAL                 :: TRange         = .FALSE.      ! specified a range of time
      LOGICAL                 :: XRes           = .FALSE.      ! specified a resolution in x
      LOGICAL                 :: YRes           = .FALSE.      ! speficied a resolution in y
      LOGICAL                 :: ZRes           = .FALSE.      ! specified a resolution in z
      LOGICAL                 :: TRes           = .FALSE.      ! specified a resolution in time
      LOGICAL                 :: ParaPrint      = .FALSE.      ! create a ParaView file?
      LOGICAL                 :: Summary        = .FALSE.      ! create a summary file?
      LOGICAL                 :: fft            = .FALSE.      ! do an FFT
      LOGICAL                 :: PointsFile     = .FALSE.      ! points file specified
   END TYPE    InflowWind_Driver_ArgFlags


      ! This contains all the settings (possible passed in arguments).
   TYPE     :: InflowWind_Driver_Args
      INTEGER                 :: WindFileType   = DEFAULT_WINDNumber ! the kind of windfile     -- set default to simplify things later
      REAL( ReKi )            :: Height                        ! Reference height
      REAL( ReKi )            :: Width                         ! Reference width
      REAL( ReKi )            :: XRange(1:2)                   ! range of x
      REAL( ReKi )            :: YRange(1:2)                   ! range of y
      REAL( ReKi )            :: ZRange(1:2)                   ! range of z
      REAL( DbKi )            :: TRange(1:2)                   ! range of time
      REAL( ReKi )            :: XRes                          ! resolution of x
      REAL( ReKi )            :: YRes                          ! resolution of y
      REAL( ReKi )            :: ZRes                          ! resolution of z
      REAL( DbKi )            :: TRes                          ! resolution of time
      REAL( ReKi )            :: fft(1:3)                      ! Coords to do an FFT
      CHARACTER(1024)         :: PointsFileName                ! Filename of points file
      CHARACTER(1024)         :: InputFileName                  ! Filename of file to process
   END TYPE    InflowWind_Driver_Args


END MODULE InflowWind_Driver_Types
