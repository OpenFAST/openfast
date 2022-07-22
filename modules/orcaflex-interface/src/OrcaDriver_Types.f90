!**********************************************************************************************************************************
!
!  MODULE: Orca_Driver_Types  - This module contains types used by the OrcaFlexInterface Driver program to store arguments passed in
!
!  The types listed here are used within the OrcaFlexInterface Driver program to store the settings. These settings are read in as
!  command line arguments, then stored within these types.
!
!**********************************************************************************************************************************
!
!..................................................................................................................................
! LICENSING
! Copyright (C) 2015  National Renewable Energy Laboratory
!
!    This file is part of OrcaFlexInterface.
!
!    OrcaFlexInterface is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
!    published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License along with OrcaFlexInterface.
!    If not, see <http://www.gnu.org/licenses/>.
!
!**********************************************************************************************************************************
! File last committed: $Date: 2014-07-29 13:30:04 -0600 (Tue, 29 Jul 2014) $
! (File) Revision #: $Rev: 169 $
! URL: $HeadURL: https://windsvn.nrel.gov/OrcaFlexInterface/branches/modularization2/Source/Driver/OrcaDriver_Types.f90 $
!**********************************************************************************************************************************

MODULE OrcaDriver_Types

   USE NWTC_Library
   USE OrcaFlexInterface_Types

   IMPLICIT NONE

      !> This contains flags to note if the settings were made.  This same data structure is
      !! used both during the driver input file and the command line options.
      !!
      !! NOTE: The WindFileType is only set if it is given as a command line option.  Otherwise
      !!       it is handled internally by InflowWInd.
      !!
      !! NOTE: The wind direction is specified by the OrcaFlexInterface input file.
   TYPE     :: OrcaDriver_Flags
      LOGICAL                 :: DvrIptFile           = .FALSE.      !< Was an input file name given on the command line?
      LOGICAL                 :: OrcaIptFile          = .FALSE.      !< Was an OrcaFlexInterface input file requested?
      LOGICAL                 :: AddedMass            = .FALSE.      !< create an added mass table at command line?
      LOGICAL                 :: AddedMassFile        = .FALSE.      !< create an added mass file?
      LOGICAL                 :: DT                   = .FALSE.      !< specified a resolution in time
      LOGICAL                 :: DTDefault            = .FALSE.      !< specified a 'DEFAULT' for the time resolution


      LOGICAL                 :: Degrees              = .FALSE.      !< angles are specified in degrees

      LOGICAL                 :: PtfmCoord            = .FALSE.      !< (x,y,z,R1,R2,R3) coordinate specified
      LOGICAL                 :: PtfmVeloc            = .FALSE.      !< (x,y,z,R1,R2,R3) coordinate specified
      LOGICAL                 :: PtfmAccel            = .FALSE.      !< (x,y,z,R1,R2,R3) coordinate specified


      LOGICAL                 :: PointsFile           = .FALSE.      !< points filename to read in
      LOGICAL                 :: PointsDegrees        = .FALSE.      !< points in the pointsfile are specified in degrees

      LOGICAL                 :: AddedMassOutputInit  = .FALSE.      !< Is the WindGridOut file initialized
      LOGICAL                 :: PointsOutputInit     = .FALSE.      !< Is the Points output file initialized
      LOGICAL                 :: Verbose              = .FALSE.      !< Verbose error reporting
      LOGICAL                 :: VVerbose             = .FALSE.      !< Very Verbose error reporting
   END TYPE    OrcaDriver_Flags


      ! This contains all the settings (possible passed in arguments).
   TYPE     :: OrcaDriver_Settings
      CHARACTER(1024)         :: DvrIptFileName                !< Driver input file name
      CHARACTER(1024)         :: OrcaIptFileName               !< Filename of OrcaFlexInterface input file to read (if no driver input file)
      CHARACTER(1024)         :: AddedMassFileName             !< Filename for the added mass matrix output

      CHARACTER(1024)         :: PointsFileName                !< Filename of points file to read in
      CHARACTER(1024)         :: PointsOutputName              !< Filename for output from points read in from points file
      INTEGER(IntKi)          :: AddedMassOutputUnit           !< Unit number for the output file for the AddedMass matrix
      INTEGER(IntKi)          :: PointsOutputUnit              !< Unit number for the output file for the Points file output
      REAL(DbKi)              :: DT                            !< resolution of time
      REAL(ReKi)              :: TMax                          !< Maximum time (we calculate this based on the number of points and timestep)

      REAL(ReKi)              :: PtfmCoord(1:6)                !< (x,y,z,R1,R2,R3) coordinate and rotations to calculate at
      REAL(ReKi)              :: PtfmVeloc(1:6)                !< instantaneous velocities corresponding to the PtfmCoord
      REAL(ReKi)              :: PtfmAccel(1:6)                !< instantaneous velocities corresponding to the PtfmCoord

      TYPE(ProgDesc)          :: ProgInfo                      !< Program info
      TYPE(ProgDesc)          :: OrcaProgInfo                  !< Program info for OrcaFlexInterface

   END TYPE    OrcaDriver_Settings


END MODULE OrcaDriver_Types
