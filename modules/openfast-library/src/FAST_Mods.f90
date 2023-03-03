!**********************************************************************************************************************************
! FAST_Prog.f90, FAST_Subs.f90, FAST_Solver.f90, FAST_Lin.f90, FAST_Types.f90, and FAST_Mods.f90 make up the FAST glue code in 
! the FAST Modularization Framework.
!..................................................................................................................................
! LICENSING
! Copyright (C) 2013-2016  National Renewable Energy Laboratory
!
!    This file is part of FAST.
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
!**********************************************************************************************************************************
MODULE FAST_ModTypes

   USE NWTC_Library
   USE FAST_Types

   TYPE(ProgDesc), PARAMETER :: FAST_Ver    = &
                                ProgDesc( 'OpenFAST', '', '' ) !< The version number of this module
         
   !..................................................................
   
   INTEGER(IntKi), PARAMETER :: Type_LandBased          = 1
   INTEGER(IntKi), PARAMETER :: Type_Offshore_Fixed     = 2
   INTEGER(IntKi), PARAMETER :: Type_Offshore_Floating  = 3
   INTEGER(IntKi), PARAMETER :: Type_MHK_Fixed          = 4
   INTEGER(IntKi), PARAMETER :: Type_MHK_Floating       = 5
   
   ! state array indexes
   INTEGER(IntKi), PARAMETER :: STATE_CURR              = 1          !< index for "current" (t_global) states
   INTEGER(IntKi), PARAMETER :: STATE_PRED              = 2          !< index for "predicted" (t_global_next) states
   
   ! VTK visualization
   INTEGER(IntKi), PARAMETER :: VTK_Unknown             = -1         !< unknown option (will produce error)
   INTEGER(IntKi), PARAMETER :: VTK_None                =  0         !< none (no VTK output)
   INTEGER(IntKi), PARAMETER :: VTK_InitOnly            =  1         !< VTK output only at initialization
   INTEGER(IntKi), PARAMETER :: VTK_Animate             =  2         !< VTK animation output
   INTEGER(IntKi), PARAMETER :: VTK_ModeShapes          =  3         !< VTK output after linearization analysis
      
   INTEGER(IntKi), PARAMETER :: VTK_Surf                =  1         !< output surfaces
   INTEGER(IntKi), PARAMETER :: VTK_Basic               =  2         !< output minimal number of point/line meshes
   INTEGER(IntKi), PARAMETER :: VTK_All                 =  3         !< output all point/line meshes
   INTEGER(IntKi), PARAMETER :: VTK_Old                 =  4         !< output in old binary format (for Matlab viewing)
   REAL(SiKi),     PARAMETER :: VTK_GroundFactor        =  4.0_SiKi  !< factor for number of rotor radii -- sets width of seabed, waves, and still water in VTK surface visualization
         
   ! linearization values
   INTEGER(IntKi), PARAMETER :: LIN_NONE                = 0          !< no inputs/outputs in linearization
   INTEGER(IntKi), PARAMETER :: LIN_STANDARD            = 1          !< use standard inputs/outputs in linearization
   INTEGER(IntKi), PARAMETER :: LIN_ALL                 = 2          !< use all inputs/outputs in linearization
   
   INTEGER(IntKi), PARAMETER :: LIN_INPUT_COL           = 1          !< index for inputs
   INTEGER(IntKi), PARAMETER :: LIN_OUTPUT_COL          = 2          !< index for outputs
   INTEGER(IntKi), PARAMETER :: LIN_ContSTATE_COL       = 3          !< index for continuous states
   
   
   INTEGER(IntKi), PARAMETER :: SizeJac_ED_HD  = 12

   LOGICAL,        PARAMETER :: BD_Solve_Option1 = .TRUE.


END MODULE FAST_ModTypes
!=======================================================================

