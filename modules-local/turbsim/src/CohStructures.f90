!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2014  National Renewable Energy Laboratory
!
!    This file is part of TurbSim.
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
!
!**********************************************************************************************************************************
MODULE TS_CohStructures

   USE                     NWTC_Library   
use ts_errors
use TS_Profiles

   IMPLICIT                NONE

   REAL(ReKi), PARAMETER   ::  KHT_LES_dT = 0.036335           ! The average time step in the LES test file, used here for the KH test
   REAL(ReKi), PARAMETER   ::  KHT_LES_Zm = 6.35475            ! The non-dimensional z dimension defined in LES test file, used here for the KH test
   

CONTAINS

SUBROUTINE CohStr_Open()
use TSMods

   

   IF ( WrACT ) THEN
      p_CohStr%ScaleWid = RotorDiameter * DistScl           !  This is the scaled height of the coherent event data set
      p_grid%Zbottom  = HubHt - CTLz*p_CohStr%ScaleWid             !  This is the height of the bottom of the wave in the scaled/shifted coherent event data set

      IF ( KHtest ) THEN      
            ! for LES test case....
         p_CohStr%ScaleVel = p_CohStr%ScaleWid * KHT_LES_dT /  KHT_LES_Zm    
         p_CohStr%ScaleVel = 50 * p_CohStr%ScaleVel                  ! We want 25 hz bandwidth so multiply by 50
      ELSE
      !   TmpPLExp = PLExp 
      !   PLExp    = MIN( 2.0, 1.35*PLExp )        ! Increase the shear of the background (?)

         p_CohStr%ScaleVel =                     getWindSpeed(UHub,HubHt,p_grid%Zbottom+p_CohStr%ScaleWid,RotorDiameter,PROFILE=WindProfileType)   ! Velocity at the top of the wave
         p_CohStr%ScaleVel = p_CohStr%ScaleVel - getWindSpeed(UHub,HubHt,p_grid%Zbottom,                  RotorDiameter,PROFILE=WindProfileType)   ! Shear across the wave
         p_CohStr%ScaleVel = 0.5 * p_CohStr%ScaleVel                                                                               ! U0 is half the difference between the top and bottom of the billow
      
      !   PLExp = TmpPLExp
      ENDIF

      p_CohStr%Uwave = getWindSpeed(UHub,HubHt,p_grid%Zbottom+0.5*p_CohStr%ScaleWid,RotorDiameter,PROFILE=WindProfileType)                 ! WindSpeed at center of wave

   !BONNIE: MAYBE WE SHOULDN'T OPEN THIS FILE UNTIL WE NEED TO WRITE TO IT
      IF (p_CohStr%ScaleVel < 0. ) THEN
         CALL TS_Warn( ' A coherent turbulence time step file cannot be generated with negative shear.', .TRUE. )
         WrACT = .FALSE.
      ENDIF
   ENDIF

END SUBROUTINE CohStr_Open



!=======================================================================

END MODULE TS_CohStructures
