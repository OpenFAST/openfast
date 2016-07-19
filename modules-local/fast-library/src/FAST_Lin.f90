!**********************************************************************************************************************************
! FAST_Solver.f90, FAST_Subs.f90, FAST_Lin.f90, and FAST_Mods.f90 make up the FAST glue code in the FAST Modularization Framework.
! FAST_Prog.f90, FAST_Library.f90, FAST_Prog.c are different drivers for this code.
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
!
!**********************************************************************************************************************************
! File last committed: $Date: 2016-05-16 11:25:44 -0600 (Mon, 16 May 2016) $
! (File) Revision #: $Rev: 1280 $
! URL: $HeadURL: https://windsvn2.nrel.gov/FAST/branches/BJonkman/Source/FAST_Subs.f90 $
!**********************************************************************************************************************************
MODULE FAST_Linear

   USE FAST_Solver  ! I mostly just want the modules that are inherited from this module, not the routines in it
   
   IMPLICIT NONE

   CONTAINS
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine that initializes some variables for linearization.
SUBROUTINE Init_Lin(p_FAST, y_FAST, m_FAST, AD, NumBl, ErrStat, ErrMsg)

   TYPE(FAST_ParameterType), INTENT(INOUT) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST              !< Output variables for the glue code
   TYPE(FAST_MiscVarType),   INTENT(INOUT) :: m_FAST              !< Miscellaneous variables
   TYPE(AeroDyn_Data),       INTENT(IN   ) :: AD                  !< AeroDyn data
   INTEGER(IntKi),           INTENT(IN)    :: NumBl               !< Number of blades (for index into ED input array)
   
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None

   INTEGER(IntKi)                          :: i, j                ! loop/temp variables
   INTEGER(IntKi)                          :: ThisModule          ! Module ID # 
   
   INTEGER(IntKi)                          :: i_u, i_u_ext        ! loop/temp variables
   INTEGER(IntKi)                          :: i_y, i_x            ! loop/temp variables

   INTEGER(IntKi)                          :: NextStart(4)        ! allocated to be size(p_FAST%LinStartIndx,2); helps compute the next starting index for the module components
   INTEGER(IntKi)                          :: ErrStat2            ! local error status
   CHARACTER(1024)                         :: ErrMsg2             ! local error message
   CHARACTER(*),             PARAMETER     :: RoutineName = 'Init_Lin' 
   
   
   ErrStat = ErrID_None
   ErrMsg = ""

   !.....................
   ! determine the number of modules that will be linearized:
   !.....................
   p_FAST%Lin_NumMods = 0 
   
      ! InflowWind is first, if activated:
   if ( p_FAST%CompInflow  == Module_IfW ) then 
      p_FAST%Lin_NumMods = p_FAST%Lin_NumMods + 1
      p_FAST%Lin_ModOrder( p_FAST%Lin_NumMods ) = Module_IfW
      
      call Init_Lin_IfW( p_FAST, y_FAST, AD%Input(1) ) ! overwrite some variables based on knowledge from glue code
                  
   end if
   
      ! ServoDyn is next, if activated:
   if ( p_FAST%CompServo  == Module_SrvD ) then 
      p_FAST%Lin_NumMods = p_FAST%Lin_NumMods + 1
      p_FAST%Lin_ModOrder( p_FAST%Lin_NumMods ) = Module_SrvD
   end if
   
   
      ! ElastoDyn is next; it is always activated:
   p_FAST%Lin_NumMods = p_FAST%Lin_NumMods + 1
   p_FAST%Lin_ModOrder( p_FAST%Lin_NumMods ) = Module_ED
   
      ! ED extended inputs
   if (p_FAST%LinInputs /= LIN_NONE) then
      call AllocAry(y_FAST%Lin%Modules(Module_ED)%Names_u_ext,1,'ED Names_u_ext', ErrStat2, ErrMsg2) 
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      call AllocAry(y_FAST%Lin%Modules(Module_ED)%op_u_ext,1,'ED op_u_ext', ErrStat2, ErrMsg2) 
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)                  
         if (ErrStat >= AbortErrLev) return
      y_FAST%Lin%Modules(Module_ED)%Names_u_ext(1) = 'BlPitchComC: collective blade-pitch command, rad' 
   end if
   
   
      ! AeroDyn is next, if activated:
   if ( p_FAST%CompAero  == Module_AD ) then 
      p_FAST%Lin_NumMods = p_FAST%Lin_NumMods + 1
      p_FAST%Lin_ModOrder( p_FAST%Lin_NumMods ) = Module_AD
      
         ! AD extended inputs
      if (p_FAST%LinInputs /= LIN_NONE) then
         call AllocAry(y_FAST%Lin%Modules(Module_AD)%Names_u_ext,3,'AD Names_u_ext', ErrStat2, ErrMsg2) 
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         call AllocAry(y_FAST%Lin%Modules(Module_AD)%op_u_ext,3,'AD op_u_ext', ErrStat2, ErrMsg2) 
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)                  
            if (ErrStat >= AbortErrLev) return
         y_FAST%Lin%Modules(Module_AD)%Names_u_ext(1) = 'HWindSpeed: horizontal wind speed, m/s' 
         y_FAST%Lin%Modules(Module_AD)%Names_u_ext(2) = 'PLexp: power-law shear exponent, -' 
         y_FAST%Lin%Modules(Module_AD)%Names_u_ext(3) = 'PropagationDir: wind propagation direction, rad'       
      end if
   end if
   
   
   !.....................
   ! determine total number of inputs/outputs/contStates:
   !.....................
   p_FAST%SizeLin = 0
   do i = 1,p_FAST%Lin_NumMods
      ThisModule = p_FAST%Lin_ModOrder( i )
      
      if (allocated(y_FAST%Lin%Modules(ThisModule)%Names_u))     p_FAST%SizeLin(ThisModule,LIN_INPUT_COL)     = size(y_FAST%Lin%Modules(ThisModule)%Names_u)
      if (allocated(y_FAST%Lin%Modules(ThisModule)%Names_u_ext)) p_FAST%SizeLin(ThisModule,LIN_ExtINPUT_COL)  = size(y_FAST%Lin%Modules(ThisModule)%Names_u_ext)
      if (allocated(y_FAST%Lin%Modules(ThisModule)%Names_y))     p_FAST%SizeLin(ThisModule,LIN_OUTPUT_COL)    = size(y_FAST%Lin%Modules(ThisModule)%Names_y)  
      if (allocated(y_FAST%Lin%Modules(ThisModule)%Names_x))     p_FAST%SizeLin(ThisModule,LIN_ContSTATE_COL) = size(y_FAST%Lin%Modules(ThisModule)%Names_x)
   end do
   if (p_FAST%CompAero == MODULE_AD) p_FAST%SizeLin(MODULE_AD,LIN_ExtINPUT_COL) = 3 ! horizontal wind speed, power-law exponent, and propagation direction (from InflowWind's steady/uniform wind types)                    
   
   do i=1,size(p_FAST%SizeLin,2)
      p_FAST%SizeLin(NumModules+1,i) = sum( p_FAST%SizeLin(1:NumModules,i) )  ! total number of inputs, extended inputs, outputs, and continuous states
   end do
                               
   !.....................
   ! compute the starting index in the combined (full) matrices:
   !.....................
   p_FAST%LinStartIndx = -1   
   NextStart = 1 ! whole array
   do i = 1,p_FAST%Lin_NumMods
      ThisModule = p_FAST%Lin_ModOrder( i )
      
      do j=1,4
         p_FAST%LinStartIndx(ThisModule,j) = NextStart(j)
         NextStart(j) = NextStart(j) + p_FAST%SizeLin(ThisModule,j)
      end do
   end do
   
   
      ! ...................................
      ! determine which of the module inputs/outputs are written to file
      ! ...................................
   !NumBl = size(u_ED%BlPitchCom)   
   call Init_Lin_InputOutput(p_FAST, y_FAST, m_FAST, NumBl, ErrStat2, ErrMsg2)
      call SetErrStat(errStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   
      ! ...................................
      ! get names of inputs, outputs, and continuous states
      ! ...................................
   call AllocAry( y_FAST%Lin%Glue%names_u, p_FAST%SizeLin(NumModules+1,LIN_INPUT_COL), 'names_u', ErrStat2, ErrMsg2)
      call SetErrStat(errStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call AllocAry( y_FAST%Lin%Glue%names_u_ext, p_FAST%SizeLin(NumModules+1,LIN_ExtINPUT_COL), 'names_u_ext', ErrStat2, ErrMsg2)
      call SetErrStat(errStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call AllocAry( y_FAST%Lin%Glue%names_y, p_FAST%SizeLin(NumModules+1,LIN_OUTPUT_COL), 'names_y', ErrStat2, ErrMsg2)
      call SetErrStat(errStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)      
   call AllocAry( y_FAST%Lin%Glue%names_x, p_FAST%SizeLin(NumModules+1,LIN_ContSTATE_COL), 'names_x', ErrStat2, ErrMsg2)
      call SetErrStat(errStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)       
   call AllocAry( y_FAST%Lin%Glue%Use_u, p_FAST%SizeLin(NumModules+1,LIN_INPUT_COL), 'use_u', ErrStat2, ErrMsg2)
      call SetErrStat(errStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)      
   call AllocAry( y_FAST%Lin%Glue%Use_y, p_FAST%SizeLin(NumModules+1,LIN_OUTPUT_COL), 'use_y', ErrStat2, ErrMsg2)
      call SetErrStat(errStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)      
   call AllocAry( y_FAST%Lin%Glue%RotFrame_u, p_FAST%SizeLin(NumModules+1,LIN_INPUT_COL), 'RotFrame_u', ErrStat2, ErrMsg2)
      call SetErrStat(errStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call AllocAry( y_FAST%Lin%Glue%RotFrame_y, p_FAST%SizeLin(NumModules+1,LIN_OUTPUT_COL), 'RotFrame_y', ErrStat2, ErrMsg2)
      call SetErrStat(errStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)      
   call AllocAry( y_FAST%Lin%Glue%RotFrame_x, p_FAST%SizeLin(NumModules+1,LIN_ContSTATE_COL), 'RotFrame_x', ErrStat2, ErrMsg2)
      call SetErrStat(errStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)       
      
   if (ErrStat >= AbortErrLev) return
               
   
   i_u = 1
   i_u_ext = 1
   i_y = 1      
   i_x = 1      
   do i = 1,p_FAST%Lin_NumMods
      ThisModule = p_FAST%Lin_ModOrder( i )
      
      do j=1,p_FAST%SizeLin(ThisModule,LIN_INPUT_COL)
         y_FAST%Lin%Glue%names_u(i_u) = TRIM(y_FAST%Module_Abrev(ThisModule))//' '//y_FAST%Lin%Modules(ThisModule)%Names_u(j)         
         y_FAST%Lin%Glue%use_u(  i_u) = y_FAST%Lin%Modules(ThisModule)%use_u(j) 
         
         if (allocated(y_FAST%Lin%Modules(ThisModule)%RotFrame_u)) then
            y_FAST%Lin%Glue%RotFrame_u(i_u) = y_FAST%Lin%Modules(ThisModule)%RotFrame_u(j) 
         else 
            y_FAST%Lin%Glue%RotFrame_u(i_u) = .false.
         end if         
         
         i_u = i_u + 1;
      end do
      
      
      do j=1,p_FAST%SizeLin(ThisModule,LIN_ExtINPUT_COL)
         y_FAST%Lin%Glue%names_u_ext(i_u_ext) = TRIM(y_FAST%Module_Abrev(ThisModule))//' '//y_FAST%Lin%Modules(ThisModule)%Names_u_ext(j)         
         i_u_ext = i_u_ext + 1;
      end do
      
      
      do j=1,p_FAST%SizeLin(ThisModule,LIN_OUTPUT_COL)
         y_FAST%Lin%Glue%names_y(i_y) = TRIM(y_FAST%Module_Abrev(ThisModule))//' '//y_FAST%Lin%Modules(ThisModule)%Names_y(j)
         y_FAST%Lin%Glue%use_y(  i_y) = y_FAST%Lin%Modules(ThisModule)%use_y(j)
         if (allocated(y_FAST%Lin%Modules(ThisModule)%RotFrame_y)) then
            y_FAST%Lin%Glue%RotFrame_y(i_y) = y_FAST%Lin%Modules(ThisModule)%RotFrame_y(j)
         else 
            y_FAST%Lin%Glue%RotFrame_y(i_y) = .false.
         end if                  
         i_y = i_y + 1;
      end do      

      do j=1,p_FAST%SizeLin(ThisModule,LIN_ContSTATE_COL)
         y_FAST%Lin%Glue%names_x( i_x) = TRIM(y_FAST%Module_Abrev(ThisModule))//' '//y_FAST%Lin%Modules(ThisModule)%Names_x( j)
         if (allocated(y_FAST%Lin%Modules(ThisModule)%RotFrame_x)) then
            y_FAST%Lin%Glue%RotFrame_x(i_x) = y_FAST%Lin%Modules(ThisModule)%RotFrame_x(j) 
         else 
            y_FAST%Lin%Glue%RotFrame_x(i_x) = .false.
         end if                  
         i_x = i_x + 1;
      end do      
      
   end do
         
   
END SUBROUTINE Init_Lin
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine that initializes the names and rotating frame portion of IfW.
SUBROUTINE Init_Lin_IfW( p_FAST, y_FAST, u_AD )

   TYPE(FAST_ParameterType),       INTENT(IN   )   :: p_FAST      !< FAST parameter data 
   TYPE(FAST_OutputFileType),      INTENT(INOUT)   :: y_FAST      !< Output variables for the glue code
   TYPE(AD_InputType),             INTENT(IN)      :: u_AD        !< The input meshes (already calculated) from AeroDyn   
   
   INTEGER(IntKi)                          :: i, j, k             ! loop counters
   INTEGER(IntKi)                          :: i2, j2              ! loop counters
   INTEGER(IntKi)                          :: Node                ! InflowWind node number
   CHARACTER(25)                           :: NodeDesc            ! Node description
   INTEGER(IntKi)                          :: position            ! position in string
   
            
         ! compare with IfW_InputSolve():
   
      Node = 0 !InflowWind node
            
      ! I'm going to overwrite some of the input/output descriptions 
      if (p_FAST%CompServo == MODULE_SrvD) then
         Node = Node + 1
         NodeDesc = ' (hub)'
         
         do i=1,3
            position = index(y_FAST%Lin%Modules(Module_IfW)%Names_u(i), ',') - 1
            y_FAST%Lin%Modules(Module_IfW)%Names_u(i) = y_FAST%Lin%Modules(Module_IfW)%Names_u(i)(1:position)//trim(NodeDesc)//&
                                                        y_FAST%Lin%Modules(Module_IfW)%Names_u(i)(position+1:)
         end do    
         do i=1,3
            position = index(y_FAST%Lin%Modules(Module_IfW)%Names_y(i), ',') - 1
            y_FAST%Lin%Modules(Module_IfW)%Names_y(i) = y_FAST%Lin%Modules(Module_IfW)%Names_y(i)(1:position)//trim(NodeDesc)//&
                                                        y_FAST%Lin%Modules(Module_IfW)%Names_y(i)(position+1:)
         end do    
      end if
                  
      IF (p_FAST%CompAero == MODULE_AD) THEN 
                           
         DO K = 1,SIZE(u_AD%BladeMotion)
            DO J = 1,u_AD%BladeMotion(k)%Nnodes
               Node = Node + 1 ! InflowWind node
               NodeDesc = ' (blade '//trim(num2lstr(k))//', node '//trim(num2lstr(j))//')'
               
               do i=1,3 !XYZ components of this node
                  i2 = (Node-1)*3 + i
                                    
                  position = index(y_FAST%Lin%Modules(Module_IfW)%Names_u(i2), ',') - 1
                  y_FAST%Lin%Modules(Module_IfW)%Names_u(i2) = y_FAST%Lin%Modules(Module_IfW)%Names_u(i2)(1:position)//trim(NodeDesc)//&
                                                               y_FAST%Lin%Modules(Module_IfW)%Names_u(i2)(position+1:)
                                                       
                  position = index(y_FAST%Lin%Modules(Module_IfW)%Names_y(i2), ',') - 1
                  y_FAST%Lin%Modules(Module_IfW)%Names_y(i2) = y_FAST%Lin%Modules(Module_IfW)%Names_y(i2)(1:position)//trim(NodeDesc)//&
                                                               y_FAST%Lin%Modules(Module_IfW)%Names_y(i2)(position+1:)
                  
                  y_FAST%Lin%Modules(Module_IfW)%RotFrame_u(i2) = .true.
                  y_FAST%Lin%Modules(Module_IfW)%RotFrame_y(i2) = .true.
                  
               end do            
            END DO !J = 1,p%BldNodes ! Loop through the blade nodes / elements
         END DO !K = 1,p%NumBl     
         
            ! tower:
         DO J=1,u_AD%TowerMotion%nnodes
            Node = Node + 1  
            NodeDesc = ' (Tower node '//trim(num2lstr(j))//')'

            do i=1,3 !XYZ components of this node
               i2 = (Node-1)*3 + i
                                    
               position = index(y_FAST%Lin%Modules(Module_IfW)%Names_u(i2), ',') - 1
               y_FAST%Lin%Modules(Module_IfW)%Names_u(i2) = y_FAST%Lin%Modules(Module_IfW)%Names_u(i2)(1:position)//trim(NodeDesc)//&
                                                            y_FAST%Lin%Modules(Module_IfW)%Names_u(i2)(position+1:)
                                     
               position = index(y_FAST%Lin%Modules(Module_IfW)%Names_y(i2), ',') - 1
               y_FAST%Lin%Modules(Module_IfW)%Names_y(i2) = y_FAST%Lin%Modules(Module_IfW)%Names_y(i2)(1:position)//trim(NodeDesc)//&
                                                            y_FAST%Lin%Modules(Module_IfW)%Names_y(i2)(position+1:)                                    
            end do            
         END DO              
         
      END IF     
   
END SUBROUTINE Init_Lin_IfW
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine that initializes some use_u and use_y, which determine which, if any, inputs and outputs are output in the linearization file.
SUBROUTINE Init_Lin_InputOutput(p_FAST, y_FAST, m_FAST, NumBl, ErrStat, ErrMsg)

   TYPE(FAST_ParameterType), INTENT(INOUT) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST              !< Output variables for the glue code
   TYPE(FAST_MiscVarType),   INTENT(INOUT) :: m_FAST              !< Miscellaneous variables
   INTEGER(IntKi),           INTENT(IN   ) :: NumBl               !< Number of blades (for index into ED input array)
   
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None

   INTEGER(IntKi)                          :: i, j, col           ! loop/temp variables
   INTEGER(IntKi)                          :: ThisModule          ! Module ID # 
   integer(intKi)                          :: ED_Start_tmp        ! starting index of ED inputs where standard inputs are located
   
   INTEGER(IntKi)                          :: ErrStat2            ! local error status
   CHARACTER(1024)                         :: ErrMsg2             ! local error message
   CHARACTER(*),             PARAMETER     :: RoutineName = 'Init_Lin_InputOutput' 
   
   
   ErrStat = ErrID_None
   ErrMsg = ""                               
   
      ! ...................................
      ! allocate module arrays
      ! ...................................
   do i = 1,p_FAST%Lin_NumMods
      ThisModule = p_FAST%Lin_ModOrder( i )
      
      call AllocAry ( y_FAST%Lin%Modules(ThisModule)%Use_u, size(y_FAST%Lin%Modules(ThisModule)%Names_u), TRIM(y_FAST%Module_Abrev(ThisModule))//'_'//'Use_u', ErrStat2, ErrMsg2)
         call SetErrStat(errStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      call AllocAry ( y_FAST%Lin%Modules(ThisModule)%Use_y, size(y_FAST%Lin%Modules(ThisModule)%Names_y), TRIM(y_FAST%Module_Abrev(ThisModule))//'_'//'Use_y', ErrStat2, ErrMsg2)
         call SetErrStat(errStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)                     
   end do
   if (ErrStat >= AbortErrLev) return
   
   
      ! ...................................
      ! set true/false flags for inputs:
      ! ...................................
   
   if (p_FAST%LinInputs == LIN_NONE) then
      do i = 1,p_FAST%Lin_NumMods
         ThisModule = p_FAST%Lin_ModOrder( i )
         y_FAST%Lin%Modules(ThisModule)%use_u = .false.
      end do      
   elseif(p_FAST%LinInputs == LIN_STANDARD) then
      do i = 1,p_FAST%Lin_NumMods
         ThisModule = p_FAST%Lin_ModOrder( i )
         y_FAST%Lin%Modules(ThisModule)%use_u = .false.
      end do      
      
      ! ED standard inputs: BlPitchCom, YawMom, GenTrq
      ED_Start_tmp = p_FAST%SizeLin(MODULE_ED,LIN_INPUT_COL) - NumBl - 2 + 1 
      do j=ED_Start_tmp, p_FAST%SizeLin(MODULE_ED,LIN_INPUT_COL) 
         y_FAST%Lin%Modules(MODULE_ED)%use_u(j) = .true.
      end do
            
      ! extended inputs allocated only when needed, so we don't need to look at use_u
      
   elseif(p_FAST%LinInputs == LIN_ALL) then
      do i = 1,p_FAST%Lin_NumMods
         ThisModule = p_FAST%Lin_ModOrder( i )
         y_FAST%Lin%Modules(ThisModule)%use_u = .true.
      end do      
   end if
            
        
      ! ...................................
      ! set true/false flags for outputs:
      ! ...................................
   
   if (p_FAST%LinOutputs == LIN_NONE) then
      do i = 1,p_FAST%Lin_NumMods
         ThisModule = p_FAST%Lin_ModOrder( i )
         y_FAST%Lin%Modules(ThisModule)%use_y = .false.
      end do      
   elseif(p_FAST%LinOutputs == LIN_STANDARD) then
      
      ! WriteOutput values are the last entries of the modules      
      do i = 1,p_FAST%Lin_NumMods         
         ThisModule = p_FAST%Lin_ModOrder( i )
         
         col = p_FAST%SizeLin(ThisModule,LIN_OUTPUT_COL) - y_FAST%NumOuts(ThisModule) !first column where WriteOutput occurs
         do j=1,col
            y_FAST%Lin%Modules(ThisModule)%use_y(j) = .false.
         end do
         do j=col+1,p_FAST%SizeLin(ThisModule,LIN_OUTPUT_COL)
            y_FAST%Lin%Modules(ThisModule)%use_y(j) = .true.
         end do
         
      end do      
      
   elseif(p_FAST%LinOutputs == LIN_ALL) then
      do i = 1,p_FAST%Lin_NumMods
         ThisModule = p_FAST%Lin_ModOrder( i )
         y_FAST%Lin%Modules(ThisModule)%use_y = .true.
      end do      
   end if
   
   
END SUBROUTINE Init_Lin_InputOutput
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine that performs lineaization at current operating point for a turbine. 
SUBROUTINE FAST_Linearize_OP(t_global, p_FAST, y_FAST, m_FAST, ED, BD, SrvD, AD14, AD, IfW, OpFM, HD, SD, MAPp, FEAM, MD, Orca, &
                         IceF, IceD, MeshMapData, ErrStat, ErrMsg )

   REAL(DbKi),               INTENT(IN   ) :: t_global            !< current (global) simulation time

   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST              !< Output variables for the glue code
   TYPE(FAST_MiscVarType),   INTENT(INOUT) :: m_FAST              !< Miscellaneous variables
     
   TYPE(ElastoDyn_Data),     INTENT(INOUT) :: ED                  !< ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(INOUT) :: BD                  !< BeamDyn data
   TYPE(ServoDyn_Data),      INTENT(INOUT) :: SrvD                !< ServoDyn data
   TYPE(AeroDyn14_Data),     INTENT(INOUT) :: AD14                !< AeroDyn14 data
   TYPE(AeroDyn_Data),       INTENT(INOUT) :: AD                  !< AeroDyn data
   TYPE(InflowWind_Data),    INTENT(INOUT) :: IfW                 !< InflowWind data
   TYPE(OpenFOAM_Data),      INTENT(INOUT) :: OpFM                !< OpenFOAM data
   TYPE(HydroDyn_Data),      INTENT(INOUT) :: HD                  !< HydroDyn data
   TYPE(SubDyn_Data),        INTENT(INOUT) :: SD                  !< SubDyn data
   TYPE(MAP_Data),           INTENT(INOUT) :: MAPp                !< MAP data
   TYPE(FEAMooring_Data),    INTENT(INOUT) :: FEAM                !< FEAMooring data
   TYPE(MoorDyn_Data),       INTENT(INOUT) :: MD                  !< Data for the MoorDyn module
   TYPE(OrcaFlex_Data),      INTENT(INOUT) :: Orca                !< OrcaFlex interface data
   TYPE(IceFloe_Data),       INTENT(INOUT) :: IceF                !< IceFloe data
   TYPE(IceDyn_Data),        INTENT(INOUT) :: IceD                !< All the IceDyn data used in time-step loop

   TYPE(FAST_ModuleMapType), INTENT(INOUT) :: MeshMapData         !< Data for mapping between modules
      
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None

      ! local variables
   INTEGER(IntKi)                          :: Un                  ! unit number for linearization output file (written in two parts)
   INTEGER(IntKi)                          :: ErrStat2            ! local error status
   CHARACTER(1024)                         :: ErrMsg2             ! local error message
   CHARACTER(*),             PARAMETER     :: RoutineName = 'FAST_Linearize_OP' 
   
   REAL(ReKi), ALLOCATABLE                 :: dYdz(:,:), dZdz(:,:), dZdu(:,:)
   REAL(ReKi), ALLOCATABLE                 :: ext(:,:)            ! extra columns of specific matrix necessary for the module's extended inputs
   REAL(ReKi), ALLOCATABLE                 :: dUdu(:,:), dUdy(:,:) ! variables for glue-code linearization
   INTEGER(IntKi), ALLOCATABLE             :: ipiv(:)
   integer(intki)                          :: nu, ny, NumBl
   CHARACTER(1024)                         :: LinRootName
   CHARACTER(1024)                         :: OutFileName
   
   
   
   ErrStat = ErrID_None
   ErrMsg = ""

   
   LinRootName = TRIM(p_FAST%OutFileRoot)//'.'//trim(num2lstr(m_FAST%NextLinTimeIndx))
   
   NumBl = size(ED%Input(1)%BlPitchCom) 
   y_FAST%Lin%RotSpeed = ED%Output(1)%RotSpeed
   y_FAST%Lin%Azimuth  = ED%Output(1)%LSSTipPxa
   !.....................
   ! ElastoDyn
   !.....................
      ! get the jacobians
   call ED_JacobianPInput( t_global, ED%Input(1), ED%p, ED%x(STATE_CURR), ED%xd(STATE_CURR), ED%z(STATE_CURR), ED%OtherSt(STATE_CURR), &
                              ED%Output(1), ED%m, ErrStat2, ErrMsg2, dYdu=y_FAST%Lin%Modules(Module_ED)%D, dXdu=y_FAST%Lin%Modules(Module_ED)%B )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
   call ED_JacobianPContState( t_global, ED%Input(1), ED%p, ED%x(STATE_CURR), ED%xd(STATE_CURR), ED%z(STATE_CURR), ED%OtherSt(STATE_CURR), &
                                  ED%Output(1), ED%m, ErrStat2, ErrMsg2, dYdx=y_FAST%Lin%Modules(Module_ED)%C, dXdx=y_FAST%Lin%Modules(Module_ED)%A )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   
      ! get the operating point
   call ED_GetOP( t_global, ED%Input(1), ED%p, ED%x(STATE_CURR), ED%xd(STATE_CURR), ED%z(STATE_CURR), ED%OtherSt(STATE_CURR), &
                     ED%Output(1), ED%m, ErrStat2, ErrMsg2, u_op=y_FAST%Lin%Modules(Module_ED)%op_u, y_op=y_FAST%Lin%Modules(Module_ED)%op_y, &
                    x_op=y_FAST%Lin%Modules(Module_ED)%op_x, dx_op=y_FAST%Lin%Modules(Module_ED)%op_dx )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (ErrStat >=AbortErrLev) then
         call cleanup()
         return
      end if
      
   if ( allocated(y_FAST%Lin%Modules(Module_ED)%op_u_ext) ) then
      call ED_OP_ext( p_FAST, y_FAST, ED%Input(1), ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         
      if (.not. allocated(y_FAST%Lin%Modules(Module_ED)%D_ext) ) then
         call AllocAry(y_FAST%Lin%Modules(Module_ED)%B_ext, size(y_FAST%Lin%Modules(Module_ED)%B,1), &
                        p_FAST%SizeLin(Module_ED, LIN_ExtINPUT_COL), 'ED ext dXdu', ErrStat2,ErrMsg2)
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         call AllocAry(y_FAST%Lin%Modules(Module_ED)%D_ext, size(y_FAST%Lin%Modules(Module_ED)%D,1), &
                        p_FAST%SizeLin(Module_ED, LIN_ExtINPUT_COL), 'ED ext dYdu', ErrStat2,ErrMsg2)
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            
         if (ErrStat >=AbortErrLev) then
            call cleanup()
            return
         end if            
      end if
         
   end if
   
      
      ! write the module matrices:
   if (p_FAST%LinOutMod) then
      
      OutFileName = trim(LinRootName)//'.'//TRIM(y_FAST%Module_Abrev(Module_ED))      
      call WrLinFile_txt_Head(t_global, p_FAST, y_FAST, y_FAST%Lin%Modules(Module_ED), OutFileName, Un, ErrStat2, ErrMsg2 )       
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat >=AbortErrLev) then
            call cleanup()
            return
         end if
         
      if (p_FAST%LinOutJac) then
         ! Jacobians
         !dXdx:
         call WrPartialMatrix( y_FAST%Lin%Modules(Module_ED)%A, Un, p_FAST%OutFmt, 'dXdx' )    
         
         !dXdu, possibly with extended inputs:
         if ( p_FAST%SizeLin(Module_ED, LIN_ExtINPUT_COL) > 0 ) then
            call ED_SetExtInputs( p_FAST, y_FAST%Lin%Modules(Module_ED)%B, 1, 1, NumBl, y_FAST%Lin%Modules(Module_ED)%B_ext)
            call WrPartialMatrix( y_FAST%Lin%Modules(Module_ED)%B, Un, p_FAST%OutFmt, 'dXdu',&
                                    UseCol=y_FAST%Lin%Modules(Module_ED)%use_u, ExtCol = y_FAST%Lin%Modules(Module_ED)%B_ext )
         else
            call WrPartialMatrix( y_FAST%Lin%Modules(Module_ED)%B, Un, p_FAST%OutFmt, 'dXdu', UseCol=y_FAST%Lin%Modules(Module_ED)%use_u )
         end if
         
         ! dYdx:
         call WrPartialMatrix( y_FAST%Lin%Modules(Module_ED)%C, Un, p_FAST%OutFmt, 'dYdx', UseRow=y_FAST%Lin%Modules(Module_ED)%use_y )
         
         !dYdu, possibly with extended inputs:
         if ( p_FAST%SizeLin(Module_ED, LIN_ExtINPUT_COL) > 0 ) then
            call ED_SetExtInputs(p_FAST, y_FAST%Lin%Modules(Module_ED)%D, 1, 1, NumBl, y_FAST%Lin%Modules(Module_ED)%D_ext)
            call WrPartialMatrix( y_FAST%Lin%Modules(Module_ED)%D, Un, p_FAST%OutFmt, 'dYdu', UseRow=y_FAST%Lin%Modules(Module_ED)%use_y, &
                                             UseCol=y_FAST%Lin%Modules(Module_ED)%use_u, ExtCol = y_FAST%Lin%Modules(Module_ED)%D_ext )
         else         
            call WrPartialMatrix( y_FAST%Lin%Modules(Module_ED)%D, Un, p_FAST%OutFmt, 'dYdu', UseRow=y_FAST%Lin%Modules(Module_ED)%use_y, &
                                                                                              UseCol=y_FAST%Lin%Modules(Module_ED)%use_u )
         end if
         
      end if
      
         ! finish writing the file
      call WrLinFile_txt_End(Un, p_FAST, y_FAST%Lin%Modules(Module_ED) )       
               
   end if
   
   !.....................
   ! InflowWind
   !.....................      
   if ( p_FAST%CompInflow  == Module_IfW ) then 
      
         ! get the jacobians
      call InflowWind_JacobianPInput( t_global, IfW%Input(1), IfW%p, IfW%x(STATE_CURR), IfW%xd(STATE_CURR), IfW%z(STATE_CURR), &
                                   IfW%OtherSt(STATE_CURR), IfW%y, IfW%m, ErrStat2, ErrMsg2, dYdu=y_FAST%Lin%Modules(Module_IfW)%D )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
      ! get the operating point
      call InflowWind_GetOP( t_global, IfW%Input(1), IfW%p, IfW%x(STATE_CURR), IfW%xd(STATE_CURR), IfW%z(STATE_CURR), &
                             IfW%OtherSt(STATE_CURR), IfW%y, IfW%m, ErrStat2, ErrMsg2, u_op=y_FAST%Lin%Modules(Module_IfW)%op_u, &
                       y_op=y_FAST%Lin%Modules(Module_IfW)%op_y )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat >=AbortErrLev) then
            call cleanup()
            return
         end if
                      
      
         ! write the module matrices:
      if (p_FAST%LinOutMod) then
               
         OutFileName = trim(LinRootName)//'.'//TRIM(y_FAST%Module_Abrev(Module_IfW))      
         call WrLinFile_txt_Head(t_global, p_FAST, y_FAST, y_FAST%Lin%Modules(Module_IfW), OutFileName, Un, ErrStat2, ErrMsg2 )       
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            if (ErrStat >=AbortErrLev) then
               call cleanup()
               return
            end if
         
         if (p_FAST%LinOutJac) then
            ! Jacobians
            !dYdu, which does not have extended inputs:
            call WrPartialMatrix( y_FAST%Lin%Modules(Module_IfW)%D, Un, p_FAST%OutFmt, 'dYdu', &
               UseRow=y_FAST%Lin%Modules(Module_IfW)%use_y, UseCol=y_FAST%Lin%Modules(Module_IfW)%use_u )
         end if
      
            ! finish writing the file
         call WrLinFile_txt_End(Un, p_FAST, y_FAST%Lin%Modules(Module_IfW) )       
               
      end if      
            
   end if
   
   !.....................
   ! ServoDyn
   !.....................   
   if ( p_FAST%CompServo  == Module_SrvD ) then 
         ! get the jacobians
      call SrvD_JacobianPInput( t_global, SrvD%Input(1), SrvD%p, SrvD%x(STATE_CURR), SrvD%xd(STATE_CURR), SrvD%z(STATE_CURR), &
                                   SrvD%OtherSt(STATE_CURR), SrvD%y, SrvD%m, ErrStat2, ErrMsg2, dYdu=y_FAST%Lin%Modules(Module_SrvD)%D )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
      ! get the operating point
      call SrvD_GetOP( t_global, SrvD%Input(1), SrvD%p, SrvD%x(STATE_CURR), SrvD%xd(STATE_CURR), SrvD%z(STATE_CURR), &
                       SrvD%OtherSt(STATE_CURR), SrvD%y, SrvD%m, ErrStat2, ErrMsg2, u_op=y_FAST%Lin%Modules(Module_SrvD)%op_u, &
                       y_op=y_FAST%Lin%Modules(Module_SrvD)%op_y )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat >=AbortErrLev) then
            call cleanup()
            return
         end if
            
         ! write the module matrices:
      if (p_FAST%LinOutMod) then
      
         OutFileName = trim(LinRootName)//'.'//TRIM(y_FAST%Module_Abrev(Module_SrvD))      
         call WrLinFile_txt_Head(t_global, p_FAST, y_FAST, y_FAST%Lin%Modules(Module_SrvD), OutFileName, Un, ErrStat2, ErrMsg2 )       
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            if (ErrStat >=AbortErrLev) then
               call cleanup()
               return
            end if
         
            ! Jacobians
         if (p_FAST%LinOutJac) then
            !dYdu, which does not have extended inputs:
            call WrPartialMatrix( y_FAST%Lin%Modules(Module_SrvD)%D, Un, p_FAST%OutFmt, 'dYdu', &
               UseRow=y_FAST%Lin%Modules(Module_SrvD)%use_y, UseCol=y_FAST%Lin%Modules(Module_SrvD)%use_u )
         end if
      
            ! finish writing the file
         call WrLinFile_txt_End(Un, p_FAST, y_FAST%Lin%Modules(Module_SrvD) )       
               
      end if      
   end if

   !.....................
   ! AeroDyn
   !.....................
   if ( p_FAST%CompAero  == Module_AD ) then 
         ! get the jacobians
      call AD_JacobianPInput( t_global, AD%Input(1), AD%p, AD%x(STATE_CURR), AD%xd(STATE_CURR), AD%z(STATE_CURR), &
                                   AD%OtherSt(STATE_CURR), AD%y, AD%m, ErrStat2, ErrMsg2, dYdu=y_FAST%Lin%Modules(Module_AD)%D, dZdu=dZdu )      
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
      call AD_JacobianPConstrState( t_global, AD%Input(1), AD%p, AD%x(STATE_CURR), AD%xd(STATE_CURR), AD%z(STATE_CURR), &
                                   AD%OtherSt(STATE_CURR), AD%y, AD%m, ErrStat2, ErrMsg2, dYdz=dYdz, dZdz=dZdz )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      ! get the operating point
      call AD_GetOP( t_global, AD%Input(1), AD%p, AD%x(STATE_CURR), AD%xd(STATE_CURR), AD%z(STATE_CURR), &
                       AD%OtherSt(STATE_CURR), AD%y, AD%m, ErrStat2, ErrMsg2, u_op=y_FAST%Lin%Modules(Module_AD)%op_u, &
                       y_op=y_FAST%Lin%Modules(Module_AD)%op_y, z_op=y_FAST%Lin%Modules(Module_AD)%op_z )               
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat >=AbortErrLev) then
            call cleanup()
            return
         end if
         
      if ( allocated(y_FAST%Lin%Modules(Module_AD)%op_u_ext) ) then
         call AD_OP_ext( p_FAST, y_FAST, IfW%y )
         
         if (.not. allocated(y_FAST%Lin%Modules(Module_AD)%D_ext)) then
            call AllocAry(y_FAST%Lin%Modules(Module_AD)%D_ext, size(y_FAST%Lin%Modules(Module_AD)%D,1),  &
                           p_FAST%SizeLin(Module_AD, LIN_ExtINPUT_COL), 'AD ext dYdu', ErrStat2,ErrMsg2)
               call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
               if (ErrStat >=AbortErrLev) then
                  call cleanup()
                  return
               end if
         end if
         
      end if
         
      
         ! write the module matrices:
      if (p_FAST%LinOutMod) then
      
         OutFileName = trim(LinRootName)//'.'//TRIM(y_FAST%Module_Abrev(Module_AD))      
         call WrLinFile_txt_Head(t_global, p_FAST, y_FAST, y_FAST%Lin%Modules(Module_AD), OutFileName, Un, ErrStat2, ErrMsg2 )       
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            if (ErrStat >=AbortErrLev) then
               call cleanup()
               return
            end if
         
         if (p_FAST%LinOutJac) then
            ! Jacobians
            ! dZdz:
            call WrPartialMatrix( dZdz, Un, p_FAST%OutFmt, 'dZdz' )
                                    
            ! dZdu:
            if ( p_FAST%SizeLin(Module_AD, LIN_ExtINPUT_COL) > 0 ) then
               call AllocAry(ext, size(dZdu,1), p_FAST%SizeLin(Module_AD, LIN_ExtINPUT_COL), 'AD ext dZdu', ErrStat2,ErrMsg2)
                  call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
                  ! i don't want to return on error while this file is open:
               if (ErrStat<AbortErrLev) then
                  call AD_SetExtInputs(p_FAST, y_FAST, dZdu, 1, 1, AD%Input(1), IfW%y, ext)
                  call WrPartialMatrix( dZdu, Un, p_FAST%OutFmt, 'dZdu', UseCol=y_FAST%Lin%Modules(Module_AD)%use_u, ExtCol = ext )
               end if            
               if (allocated(ext)) deallocate(ext)
            else         
               call WrPartialMatrix( dZdu, Un, p_FAST%OutFmt, 'dZdu', UseCol=y_FAST%Lin%Modules(Module_AD)%use_u )
            end if   
            
            ! dYdz:
            call WrPartialMatrix( dYdz, Un, p_FAST%OutFmt, 'dYdz', UseRow=y_FAST%Lin%Modules(Module_AD)%use_y )
            
            !dYdu, possibly with extended inputs:
            if ( p_FAST%SizeLin(Module_AD, LIN_ExtINPUT_COL) > 0 ) then               
               call AD_SetExtInputs(p_FAST, y_FAST, y_FAST%Lin%Modules(Module_AD)%D, 1, 1, AD%Input(1), IfW%y, y_FAST%Lin%Modules(Module_AD)%D_ext)                  
               call WrPartialMatrix( y_FAST%Lin%Modules(Module_AD)%D, Un, p_FAST%OutFmt, 'dYdu', &
                  UseRow=y_FAST%Lin%Modules(Module_AD)%use_y, UseCol=y_FAST%Lin%Modules(Module_AD)%use_u, ExtCol=y_FAST%Lin%Modules(Module_AD)%D_ext )
            else         
               call WrPartialMatrix( y_FAST%Lin%Modules(Module_AD)%D, Un, p_FAST%OutFmt, 'dYdu', &
                  UseRow=y_FAST%Lin%Modules(Module_AD)%use_y, UseCol=y_FAST%Lin%Modules(Module_AD)%use_u )
            end if   
            
         end if
         
      end if
      
         
      call allocAry( ipiv, size(dZdz,1), 'ipiv', ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat >= AbortErrLev) then
            call cleanup() 
            return
         end if
                  
      CALL LAPACK_getrf( M=size(dZdz,1), N=size(dZdz,2), A=dZdz, IPIV=ipiv, ErrStat=ErrStat2, ErrMsg=ErrMsg2 )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat >= AbortErrLev) then
            call cleanup() 
            return
         end if
         
      CALL LAPACK_getrs( trans='N', N=size(dZdz,2), A=dZdz, IPIV=ipiv, B=dZdu, ErrStat=ErrStat2, ErrMsg=ErrMsg2 )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      ! note that after the above solve, dZdu is now matmul(dZdz^-1, dZdu)
      y_FAST%Lin%Modules(Module_AD)%D = y_FAST%Lin%Modules(Module_AD)%D - matmul(dYdz, dZdu )
      
      ! extended inputs:
      if ( p_FAST%SizeLin(Module_AD, LIN_ExtINPUT_COL) > 0 ) then
         call AD_SetExtInputs(p_FAST, y_FAST, y_FAST%Lin%Modules(Module_AD)%D, 1, 1, AD%Input(1), IfW%y, y_FAST%Lin%Modules(Module_AD)%D_ext)                  
      end if
      
      
      if (p_FAST%LinOutMod) then
            ! finish writing the file
         call WrLinFile_txt_End(Un, p_FAST, y_FAST%Lin%Modules(Module_AD) )            
      end if
      
         ! AD doesn't need these any more, and we may need them for other modules
      if (allocated(dYdz)) deallocate(dYdz)
      if (allocated(dZdz)) deallocate(dZdz)
      if (allocated(dZdu)) deallocate(dZdu)
      if (allocated(ipiv)) deallocate(ipiv)     
      
   end if
   
   !.....................
   ! Linearization of glue code Input/Output solve:
   !.....................
   
   !.....................
   ! Glue code (currently a linearization of SolveOption2):
   ! Make sure we avoid any case where the operating point values change earlier in this routine (e.g., by calling the module Jacobian routines).
   !.....................

   call Glue_GetOP(p_FAST, y_FAST, ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (ErrStat >=AbortErrLev) then
         call cleanup()
         return
      end if
   
      ! get the dUdu and dUdy matrices, which linearize SolveOption2 for the modules we've included in linearization
   call Glue_Jacobians( t_global, p_FAST, y_FAST, m_FAST, ED, BD, SrvD, AD14, AD, IfW, OpFM, HD, SD, MAPp, FEAM, MD, Orca, &
                         IceF, IceD, MeshMapData, dUdu, dUdy, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (ErrStat >=AbortErrLev) then
         call cleanup()
         return
      end if
      
      
         
   call WrLinFile_txt_Head(t_global, p_FAST, y_FAST, y_FAST%Lin%Glue, LinRootName, Un, ErrStat2, ErrMsg2 )       
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (ErrStat >=AbortErrLev) then
         call cleanup()
         return
      end if
      
   
   if (p_FAST%LinOutJac) then
      ! Jacobians
      if ( p_FAST%SizeLin(NumModules+1, LIN_ExtINPUT_COL) > 0 ) then
         call AllocAry(ext, size(dUdu,1), p_FAST%SizeLin(NumModules+1, LIN_ExtINPUT_COL), 'ext dUdu', ErrStat2,ErrMsg2)
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            ! i don't want to return on error while this file is open:
         if (ErrStat<AbortErrLev) then
            
            call ED_SetExtInputs(p_FAST, dUdu, p_FAST%LinStartIndx(MODULE_ED, LIN_INPUT_COL), &
                          p_FAST%LinStartIndx(MODULE_ED, LIN_ExtINPUT_COL), NumBl, ext)
            if ( p_FAST%SizeLin(Module_AD, LIN_ExtINPUT_COL) > 0 ) then         
               call AD_SetExtInputs(p_FAST, y_FAST, dUdu, p_FAST%LinStartIndx(MODULE_AD, LIN_INPUT_COL), &
                                    p_FAST%LinStartIndx(MODULE_AD, LIN_ExtINPUT_COL), AD%Input(1), IfW%y, ext)                  
            end if      
                        
            call WrPartialMatrix( dUdu, Un, p_FAST%OutFmt, 'dUdu', UseRow=y_FAST%Lin%Glue%use_u, UseCol=y_FAST%Lin%Glue%use_u, ExtCol = ext )
         end if            
         if (allocated(ext)) deallocate(ext)
      else         
         call WrPartialMatrix( dUdu, Un, p_FAST%OutFmt, 'dUdu', UseRow=y_FAST%Lin%Glue%use_u, UseCol=y_FAST%Lin%Glue%use_u )
      end if   
      call WrPartialMatrix( dUdy, Un, p_FAST%OutFmt, 'dUdy', UseRow=y_FAST%Lin%Glue%use_u, UseCol=y_FAST%Lin%Glue%use_y )
   end if
   
   
      ! calculate the glue-code state matrices
   call Glue_StateMatrices( p_FAST, y_FAST, NumBl, AD, IfW, dUdu, dUdy, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (ErrStat >=AbortErrLev) then
         call cleanup()
         return
      end if
      
      ! Write the results to the file:
   call WrLinFile_txt_End(Un, p_FAST, y_FAST%Lin%Glue )            

contains
   subroutine cleanup()
      if (allocated(dYdz)) deallocate(dYdz)
      if (allocated(dZdz)) deallocate(dZdz)
      if (allocated(dZdu)) deallocate(dZdu)
      if (allocated(ipiv)) deallocate(ipiv)     
      if (allocated(ext))  deallocate(ext)
      
      if (allocated(dUdu)) deallocate(dUdu)
      if (allocated(dUdy)) deallocate(dUdy)
   end subroutine cleanup
END SUBROUTINE FAST_Linearize_OP   
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine calculates the operating point of the ED extended inputs from the operating point of the ED inputs
SUBROUTINE ED_OP_ext( p_FAST, y_FAST, u_ED, ErrStat, ErrMsg )
   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< parameters
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST              !< Output variables for the glue code
   TYPE(ED_InputType),       INTENT(IN   ) :: u_ED                !< ED Inputs at t
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None
   
   integer(intKi)                          :: i                   ! loop counter
   
   ErrStat = ErrID_None
   ErrMsg = ""
   
   y_FAST%Lin%Modules(Module_ED)%names_u_ext(1) = 'BlPitchComC: collective blade-pitch command, rad'   
   y_FAST%Lin%Modules(Module_ED)%op_u_ext(1)    = u_ED%BlPitchCom(1)
   
   do i = 2,size(u_ED%BlPitchCom)
      if (.not. EqualRealNos( u_ED%BlPitchCom(1), u_ED%BlPitchCom(i) ) ) then
         ErrStat = ErrID_Info
         ErrMsg = "Operating point of ElastoDyn's collective pitch extended input is invalid because "//&
                  "the commanded blade pitch angles are not the same for each blade."         
         y_FAST%Lin%Modules(Module_ED)%names_u_ext(1) = 'BlPitchComC: collective blade-pitch command, rad [NOTE INVALID OPERATING POINT]'
         return
      end if      
   end do
      
   
END SUBROUTINE ED_OP_ext 
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine calculates the extended-input columns of inMat, stored in ext.
SUBROUTINE ED_SetExtInputs( p_FAST, inMat, ED_Start, ED_Ext_Start, NumBl, ext)
   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< parameters
   REAL(ReKi),               INTENT(IN   ) :: inMat(:,:)          !< Matrix that contains the ED matrix of du (delta inputs)
   INTEGER(IntKi),           INTENT(IN   ) :: ED_Start            !< Index where ED inputs begin in input matrix, inMat
   INTEGER(IntKi),           INTENT(IN   ) :: ED_Ext_Start        !< Index where ED external inputs begin in output matrix, ext
   INTEGER(IntKi),           INTENT(IN   ) :: NumBl               !< Number of blades (to help get the correct index into ED inputs)
   REAL(ReKi),               INTENT(INOUT) :: ext(:,:)            !< Matrix that contains the extended inputs to be computed
   
   integer(intKi)                          :: ED_Start_BLPitchCom ! start of blade pitch column in ED
   integer(intKi)                          :: i                   ! loop counter
   
   
   ED_Start_BLPitchCom = ED_Start + p_FAST%SizeLin(Module_ED,LIN_INPUT_COL) - NumBl - 2 + 1 ! last NumBl+2 columns are: GenTrq, YawMom, and BlPitchCom
   
      ! collective pitch
   ext(:,1) = inMat(:,ED_Start_BLPitchCom)   
   do i=2,NumBl
      ED_Start_BLPitchCom = ED_Start_BLPitchCom + 1
      ext(:,ED_Ext_Start) = ext(:,ED_Ext_Start) + inMat(:,ED_Start_BLPitchCom)
   end do
      
END SUBROUTINE ED_SetExtInputs 
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine calculates the operating point of the AD extended inputs.
SUBROUTINE AD_OP_ext( p_FAST, y_FAST, y_IfW )
   TYPE(FAST_ParameterType),     INTENT(IN   ) :: p_FAST              !< parameters
   TYPE(FAST_OutputFileType),    INTENT(INOUT) :: y_FAST              !< Output variables for the glue code
   TYPE(InflowWind_OutputType),  INTENT(IN   ) :: y_IfW               !< ED Inputs at t
   
   integer(intKi)                              :: i                   ! loop counter
   
   
   if (p_FAST%CompInflow /= Module_IfW) then         
      y_FAST%Lin%Modules(Module_AD)%op_u_ext = 0.0_ReKi
   else   
      y_FAST%Lin%Modules(Module_AD)%op_u_ext(1) = y_IfW%UniformWind%V
      y_FAST%Lin%Modules(Module_AD)%op_u_ext(2) = y_IfW%UniformWind%VSHR
      y_FAST%Lin%Modules(Module_AD)%op_u_ext(3) = y_FAST%Lin%PropagationDir
   end if   
      
   
END SUBROUTINE AD_OP_ext 
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine calculates the AD extended-input columns of inMat, stored in ext. (HWindSpeed, PlExp, and PropagationDir)
SUBROUTINE AD_SetExtInputs( p_FAST, y_FAST, inMat, AD_Start, AD_Ext_Start, u_AD, y_IfW, ext)

   TYPE(FAST_ParameterType),     INTENT(IN   ) :: p_FAST              !< parameters
   TYPE(FAST_OutputFileType),    INTENT(IN   ) :: y_FAST              !< Output variables for the glue code
   REAL(ReKi),                   INTENT(IN   ) :: inMat(:,:)          !< Matrix that contains the AD matrix of du (delta inputs)
   INTEGER(IntKi),               INTENT(IN   ) :: AD_Start            !< Index where AD inputs begin in input matrix, inMat
   INTEGER(IntKi),               INTENT(IN   ) :: AD_Ext_Start        !< Index where AD external inputs begin in output matrix, ext
   TYPE(AD_InputType),           INTENT(IN   ) :: u_AD                !< AD inputs at t
   TYPE(InflowWind_OutputType),  INTENT(IN   ) :: y_IfW               !< IfW outputs (operating point) at t
   REAL(ReKi),                   INTENT(INOUT) :: ext(:,:)            !< Matrix that contains the extended inputs to be computed
   
   integer(intKi)                              :: i, j, k             ! loop counter
   integer(intKi)                              :: node_start, node    ! IfW node counter
   INTEGER(IntKi)                              :: i2, j2              ! loop counters
   INTEGER(IntKi)                              :: AD_Start_Bl         ! starting index of dUdu (column) where AD blade motion inputs are located
   INTEGER(IntKi)                              :: AD_Start_IfW        ! starting index of dUdu (column) where AD (blade) wind velocity inputs are located
   
   integer(intKi)                              :: i_HWindSpeed , i_PLexp, i_PropDir
   
   real(reKi)                                  :: Position(3)         ! absolute position of node
   real(reKi)                                  :: HWindSpeed(3)       ! partial derivatives with respect to HWindSpeed
   real(reKi)                                  :: PLexp(3)            ! partial derivatives with respect to PLexp
   real(reKi)                                  :: PropDir(3)          ! partial derivatives with respect to PropDir
      
   
   
      ! columns of the AD extended inputs in ext(:,:)
   i_HWindSpeed = AD_Ext_Start
   i_PLexp      = AD_Ext_Start + 1
   i_PropDir    = AD_Ext_Start + 2
   
   ext(:,i_HWindSpeed) = 0.0_ReKi
   ext(:,i_PLexp     ) = 0.0_ReKi
   ext(:,i_PropDir   ) = 0.0_ReKi
   
   
   if (p_FAST%CompInflow == MODULE_IfW) then
      !-------------------------------------------------------------------------------------------------
      ! Set the inputs from inflow wind:
      !-------------------------------------------------------------------------------------------------
                     
         ! blades:
      AD_Start_IfW = AD_Start &
                  + u_AD%TowerMotion%NNodes * 9  & ! 3 fields (MASKID_TRANSLATIONDISP,MASKID_Orientation,MASKID_TRANSLATIONVel) with 3 components
                  + u_AD%HubMotion%NNodes   * 9    ! 3 fields (MASKID_TRANSLATIONDISP,MASKID_Orientation,MASKID_RotationVel) with 3 components
   
      do k = 1,size(u_AD%BladeRootMotion)         
         AD_Start_IfW = AD_Start_IfW + u_AD%BladeRootMotion(k)%NNodes * 3 ! 1 field (MASKID_Orientation) with 3 components
      end do
                 
         ! get the index into the UVW fields (u_AD%InflowOnBlade)
      DO k=1,size(u_AD%BladeMotion)         
         AD_Start_IfW = AD_Start_IfW + u_AD%BladeMotion(k)%NNodes * 6 ! 3 fields (TranslationDisp, Orientation, TranslationVel) with 3 components            
      END DO
      
      do k=1,size(u_AD%InflowOnBlade,3) ! blades (SIZE(u_AD%BladeMotion))
         do j=1,size(u_AD%InflowOnBlade,2) ! nodes (u_AD%BladeMotion(k)%Nnodes)
            Position = u_AD%BladeMotion(k)%TranslationDisp(:,j) + u_AD%BladeMotion(k)%Position(:,j) 

            ! get Jacobians from this 
            call AD_ExtInputs_JacPHWindSpeed(p_FAST, y_FAST, Position, y_IfW, HWindSpeed )
            call AD_ExtInputs_JacPPLexp(p_FAST, y_FAST, Position, y_IfW, PLexp )
            call AD_ExtInputs_JacPPropDir(p_FAST, y_FAST, Position, y_IfW, PropDir )   

            do i=1,3 !velocity component               
               j2 = AD_Start_IfW + i - 1
               
               ext(:,i_HWindSpeed) = ext(:,i_HWindSpeed) + inMat(:,j2) * HWindSpeed(i)
               ext(:,i_PLexp     ) = ext(:,i_PLexp     ) + inMat(:,j2) * PLexp(i)
               ext(:,i_PropDir   ) = ext(:,i_PropDir   ) + inMat(:,j2) * PropDir(i)
                              
            end do
            AD_Start_IfW = AD_Start_IfW + 3
         end do         
      end do
      
               
         ! tower:
      if ( allocated(u_AD%InflowOnTower) ) then         
         do j=1,size(u_AD%InflowOnTower,2) !nodes (u_AD%TowerMotion%nnodes)
            
            Position = u_AD%TowerMotion%TranslationDisp(:,j) + u_AD%TowerMotion%Position(:,j) 
            ! get jacobians
            call AD_ExtInputs_JacPHWindSpeed(p_FAST, y_FAST, Position, y_IfW, HWindSpeed )
            call AD_ExtInputs_JacPPLexp(p_FAST, y_FAST, Position, y_IfW, PLexp )
            call AD_ExtInputs_JacPPropDir(p_FAST, y_FAST, Position, y_IfW, PropDir )            
            
            do i=1,3 !velocity component               
               j2 = AD_Start_IfW + i - 1
               
               ext(:,i_HWindSpeed) = ext(:,i_HWindSpeed) + inMat(:,j2) * HWindSpeed(i)
               ext(:,i_PLexp     ) = ext(:,i_PLexp     ) + inMat(:,j2) * PLexp(i)
               ext(:,i_PropDir   ) = ext(:,i_PropDir   ) + inMat(:,j2) * PropDir(i)               
            end do
            AD_Start_IfW = AD_Start_IfW + 3
         end do      
      end if
      
   end if
   
      
END SUBROUTINE AD_SetExtInputs 
!----------------------------------------------------------------------------------------------------------------------------------
!> this routine computes partial derivatives for the AD extended inputs from steady/uniform inflowWind.
subroutine AD_ExtInputs_JacPHWindSpeed(p_FAST, y_FAST, InputPosition, Ifw_OP, dVdH )
   TYPE(FAST_ParameterType),    INTENT(IN   ) :: p_FAST              !< parameters
   TYPE(FAST_OutputFileType),   INTENT(IN   ) :: y_FAST              !< Output variables for the glue code
   REAL(ReKi),                  INTENT(IN   ) :: InputPosition(3)    !< input position where this jacobian is evaluated
   TYPE(InflowWind_OutputType), INTENT(IN   ) :: Ifw_OP              !< IfW outputs (operating point) at t
   REAL(ReKi),                  INTENT(  OUT) :: dVdH(3)             !< partial derivative of 3 wind velocity components (U-V-W) w.r.t. horizontal wind speed (uniform/steady wind)
   
   REAL(ReKi)                                 :: Position(3)         !< position (rotated by PropagationDir) where this jacobian is evaluated
   REAL(ReKi)                                 :: dV1dH               ! partial derivative of V1 (in ifw_uniformwind.f90::getwindspeed) w.r.t. horizontal wind speed         
   REAL(ReKi)                                 :: CosDelta            ! cosine of Delta
   REAL(ReKi)                                 :: SinDelta            ! sine of Delta
   
   
   if ( InputPosition(3) < 0.0_ReKi .or. EqualRealNos(InputPosition(3), 0.0_ReKi)) then
         dVdH(3) = 0.0
   else
      CosDelta = COS( Ifw_OP%UniformWind%Delta )
      SinDelta = SIN( Ifw_OP%UniformWind%Delta )
      
      Position(1) = InputPosition(1)*y_FAST%Lin%cosPropDir - InputPosition(2)*y_FAST%Lin%sinPropDir
      Position(2) = InputPosition(1)*y_FAST%Lin%sinPropDir + InputPosition(2)*y_FAST%Lin%cosPropDir
      Position(3) = InputPosition(3)
   
      dV1dH =  ( ( Position(3)/y_FAST%Lin%RefHt ) ** Ifw_OP%UniformWind%VShr &                          ! power-law wind shear
             + ( Ifw_OP%UniformWind%HShr   * ( Position(2) * CosDelta + Position(1) * SinDelta ) &      ! horizontal linear shear
             +  Ifw_OP%UniformWind%VLinShr * ( Position(3)-y_FAST%Lin%RefHt ) )/y_FAST%Lin%RefLength  ) ! vertical linear shear
               
   
      dVdH(1) =  dV1dH*( y_FAST%Lin%cosPropDir*CosDelta - y_FAST%Lin%sinPropDir*SinDelta )
      dVdH(2) = -dV1dH*( y_FAST%Lin%sinPropDir*CosDelta + y_FAST%Lin%cosPropDir*SinDelta )
      dVdH(3) = 0.0_ReKi
   end if
   
end subroutine AD_ExtInputs_JacPHWindSpeed
!----------------------------------------------------------------------------------------------------------------------------------
!> this routine computes partial derivatives for the AD extended inputs from steady/uniform inflowWind.
subroutine AD_ExtInputs_JacPPLexp(p_FAST, y_FAST, InputPosition, Ifw_OP, dVdPLexp )
   TYPE(FAST_ParameterType),    INTENT(IN   ) :: p_FAST              !< parameters
   TYPE(FAST_OutputFileType),   INTENT(IN   ) :: y_FAST              !< Output variables for the glue code
   REAL(ReKi),                  INTENT(IN   ) :: InputPosition(3)    !< input position where this jacobian is evaluated
   TYPE(InflowWind_OutputType), INTENT(IN   ) :: Ifw_OP              !< IfW outputs (operating point) at t
   REAL(ReKi),                  INTENT(  OUT) :: dVdPLexp(3)         !< partial derivative of 3 wind velocity components (U-V-W) w.r.t. power-law exponent (uniform/steady wind)
   
   REAL(ReKi)                                 :: Position(3)         !< position (rotated by PropagationDir) where this jacobian is evaluated
   REAL(ReKi)                                 :: dV1dPlExp           ! partial derivative of V1 (in ifw_uniformwind.f90::getwindspeed) w.r.t. power-law exponent          
   REAL(ReKi)                                 :: CosDelta            ! cosine of Delta
   REAL(ReKi)                                 :: SinDelta            ! sine of Delta
   
      
   if ( InputPosition(3) < 0.0_ReKi .or. EqualRealNos(InputPosition(3), 0.0_ReKi)) then
         dVdPLexp(3) = 0.0
   else
      CosDelta = COS( Ifw_OP%UniformWind%Delta )
      SinDelta = SIN( Ifw_OP%UniformWind%Delta )
      
      !Position(1) = InputPosition(1)*y_FAST%Lin%cosPropDir - InputPosition(2)*y_FAST%Lin%sinPropDir
      !Position(2) = InputPosition(1)*y_FAST%Lin%sinPropDir + InputPosition(2)*y_FAST%Lin%cosPropDir
      Position(3) = InputPosition(3)
   
      dV1dPlExp =  Ifw_OP%UniformWind%V * ( Position(3)/y_FAST%Lin%RefHt ) ** Ifw_OP%UniformWind%VShr * log(Position(3)/y_FAST%Lin%RefHt)
               
      dVdPLexp(1) =  dV1dPlExp*( y_FAST%Lin%cosPropDir*CosDelta - y_FAST%Lin%sinPropDir*SinDelta )
      dVdPLexp(2) = -dV1dPlExp*( y_FAST%Lin%sinPropDir*CosDelta + y_FAST%Lin%cosPropDir*SinDelta )
      dVdPLexp(3) = 0.0_ReKi
   end if
   
end subroutine AD_ExtInputs_JacPPLexp
!----------------------------------------------------------------------------------------------------------------------------------
!> this routine computes partial derivatives for the AD extended inputs from steady/uniform inflowWind.
subroutine AD_ExtInputs_JacPPropDir(p_FAST, y_FAST, InputPosition, Ifw_OP, dVdPropDir )
   TYPE(FAST_ParameterType),    INTENT(IN   ) :: p_FAST              !< parameters
   TYPE(FAST_OutputFileType),   INTENT(IN   ) :: y_FAST              !< Output variables for the glue code
   REAL(ReKi),                  INTENT(IN   ) :: InputPosition(3)    !< input position where this jacobian is evaluated
   TYPE(InflowWind_OutputType), INTENT(IN   ) :: Ifw_OP              !< IfW outputs (operating point) at t
   REAL(ReKi),                  INTENT(  OUT) :: dVdPropDir(3)       !< partial derivative of 3 wind velocity components (U-V-W) w.r.t. propagation direction (uniform/steady wind)
   
   REAL(ReKi)                                 :: Position(3)         !< position (rotated by PropagationDir) where this jacobian is evaluated
   REAL(ReKi)                                 :: CosDelta            ! cosine of Delta
   REAL(ReKi)                                 :: SinDelta            ! sine of Delta   
   REAL(ReKi)                                 :: V1                           
   REAL(ReKi)                                 :: dV1dPD                           
   REAL(ReKi)                                 :: dP1, dP2                           
   REAL(ReKi)                                 :: f1, f2                           
   
   
   if ( InputPosition(3) < 0.0_ReKi .or. EqualRealNos(InputPosition(3), 0.0_ReKi)) then
      dVdPropDir(3) = 0.0
   else
      CosDelta = COS( Ifw_OP%UniformWind%Delta )
      SinDelta = SIN( Ifw_OP%UniformWind%Delta )
      
      Position(1) = InputPosition(1)*y_FAST%Lin%cosPropDir - InputPosition(2)*y_FAST%Lin%sinPropDir
      Position(2) = InputPosition(1)*y_FAST%Lin%sinPropDir + InputPosition(2)*y_FAST%Lin%cosPropDir
      Position(3) = InputPosition(3)
      
      dP1 = -Position(2) !-InputPosition(1)*y_FAST%Lin%sinPropDir - InputPosition(2)*y_FAST%Lin%cosPropDir
      dP2 =  Position(1) ! InputPosition(1)*y_FAST%Lin%cosPropDir - InputPosition(2)*y_FAST%Lin%sinPropDir
   
      V1 = Ifw_OP%UniformWind%V * ( ( Position(3)/y_FAST%Lin%RefHt ) ** Ifw_OP%UniformWind%VShr &         ! power-law wind shear
            + ( Ifw_OP%UniformWind%HShr   * ( Position(2) * CosDelta + Position(1) * SinDelta ) &         ! horizontal linear shear
            +  Ifw_OP%UniformWind%VLinShr * ( Position(3)-y_FAST%Lin%RefHt ) )/y_FAST%Lin%RefLength  ) &  ! vertical linear shear
            +  Ifw_OP%UniformWind%VGust                                                                   ! gust speed
         
      dV1dPD = Ifw_OP%UniformWind%V * Ifw_OP%UniformWind%HShr/y_FAST%Lin%RefLength * ( dP2*CosDelta + dP1*SinDelta )
            
      f1 = y_FAST%Lin%cosPropDir*CosDelta - y_FAST%Lin%sinPropDir*SinDelta 
      f2 = y_FAST%Lin%sinPropDir*CosDelta + y_FAST%Lin%cosPropDir*SinDelta 
      dVdPropDir(1) =  dV1dPD * f1 - V1 * f2
      dVdPropDir(2) = -dV1dPD * f2 - V1 * f1
      dVdPropDir(3) = 0.0_ReKi
      
   end if
   
end subroutine AD_ExtInputs_JacPPropDir
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine that writes the A,B,C,D matrices from linearization to a text file. 
SUBROUTINE WrLinFile_txt_Head(t_global, p_FAST, y_FAST, LinData, FileName, Un, ErrStat, ErrMsg)

   INTEGER(IntKi),           INTENT(  OUT) :: Un                  !< unit number
   REAL(DbKi),               INTENT(IN   ) :: t_global            !< current (global) simulation time
   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< parameters
   TYPE(FAST_OutputFileType),INTENT(IN   ) :: y_FAST              !< Output variables for the glue code
   TYPE(FAST_LinType),       INTENT(IN   ) :: LinData             !< Linearization data for individual module or glue (coupled system)
   CHARACTER(*),             INTENT(IN   ) :: FileName            !< root name of the linearization file to open for writing
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None

      ! local variables
   INTEGER(IntKi)                          :: i                   ! loop counter
   
   INTEGER(IntKi)                          :: ErrStat2            ! local error status
   CHARACTER(1024)                         :: ErrMsg2             ! local error message
   CHARACTER(*),             PARAMETER     :: RoutineName = 'WrLinFile_txt_Head'    
   INTEGER(IntKi)                          :: n(6)                ! sizes of arrays to print
   CHARACTER(*),             PARAMETER     :: TypeNames(6) = (/ 'continuous states', &
                                                                'discrete states  ', &
                                                                'constraint states', &
                                                                'inputs           ', &
                                                                'extended inputs  ', &
                                                                'outputs          '  /)
   CHARACTER(50)                           :: Fmt
   CHARACTER(32)                           :: Desc
   
   integer, parameter :: Indx_x      = 1
   integer, parameter :: Indx_xd     = 2
   integer, parameter :: Indx_z      = 3
   integer, parameter :: Indx_u      = 4
   integer, parameter :: Indx_u_ext  = 5
   integer, parameter :: Indx_y      = 6
   
                  
   ErrStat = ErrID_None
   ErrMsg = ""
         
   n = 0;
   if (allocated(LinData%names_x )) n(Indx_x) = size(LinData%names_x )
   if (allocated(LinData%names_xd)) n(Indx_xd) = size(LinData%names_xd)
   if (allocated(LinData%names_z )) n(Indx_z) = size(LinData%names_z )
   !if (allocated(LinData%names_u )) n(Indx_u) = size(LinData%names_u )
   if (allocated(LinData%names_u_ext)) n(Indx_u_ext) = size(LinData%names_u_ext )
   !if (allocated(LinData%names_y )) n(Indx_y) = size(LinData%names_y )
   
   if (allocated(LinData%names_u )) then
      do i=1,size(LinData%use_u)
         if (LinData%use_u(i)) n(Indx_u) = n(Indx_u)+1
      end do
   end if
   
   if (allocated(LinData%names_y )) then
      do i=1,size(LinData%use_y)
         if (LinData%use_y(i)) n(Indx_y) = n(Indx_y)+1
      end do
   end if
   
   
   CALL GetNewUnit( Un, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (ErrStat >= AbortErrLev) then
         Un = -1
         return
      end if
      
   CALL OpenFOutFile ( Un, TRIM(FileName)//'.lin', ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (ErrStat >= AbortErrLev) then
         close( Un ) 
         Un = -1
         return
      end if
      
   ! heading
      
      ! Add some file information:

   WRITE (Un,'(/,A)')  'Linearized model: '//TRIM( y_FAST%FileDescLines(1) )
   WRITE (Un,'(1X,A)') TRIM( y_FAST%FileDescLines(2) )
   WRITE (Un,'()' )    !print a blank line
   WRITE (Un,'(A)'   ) TRIM( y_FAST%FileDescLines(3) )
   WRITE (Un,'()' )    !print a blank line
   
   WRITE(Un, '(A)') 'Simulation information:'      
   fmt = '(3x,A,1x,'//trim(p_FAST%OutFmt_t)//',1x,A)' 
   Desc = 'Simulation time:'; WRITE (Un, fmt) Desc, t_global, 's'
   Desc = 'Rotor Speed:';     WRITE (Un, fmt) Desc, y_FAST%Lin%RotSpeed, 'rad/s'
   Desc = 'Azimuth:';         WRITE (Un, fmt) Desc, y_FAST%Lin%Azimuth,  'rad'
   
   fmt = '(3x,A,1x,I5)'
   do i=1,size(n)
      Desc = 'Number of '//trim(TypeNames(i))//':'
      WRITE(Un, fmt) Desc, n(i)
   end do
   
   Desc = 'Jacobians included in this file?'
   fmt  = '(3x,A,1x,A5)'
   if (p_FAST%LinOutJac) then
      write (Un, fmt) Desc, 'Yes'
   else      
      write (Un, fmt) Desc, 'No'
   end if   
      
   WRITE (Un,'()' )    !print a blank line

   
         !......................................................
         ! Write the names of the output parameters on one line:
         !......................................................
   if (n(Indx_x) > 0) then
      WRITE(Un, '(A)') 'Order of continuous states:'      
      call WrLinFile_txt_Table(p_FAST, Un, "Row/Column", LinData%op_x, LinData%names_x, rotFrame=LinData%RotFrame_x  )      
      
      WRITE(Un, '(A)') 'Order of continuous state derivatives:'      
      call WrLinFile_txt_Table(p_FAST, Un, "Row/Column", LinData%op_dx, LinData%names_x, rotFrame=LinData%RotFrame_x, deriv=.true.  )      
   end if
   
   if (n(Indx_xd) > 0) then
      WRITE(Un, '(A)') 'Order of discrete states:'      
      call WrLinFile_txt_Table(p_FAST, Un, "Row/Column", LinData%op_xd, LinData%names_xd  )      
   end if

   if (n(Indx_z) > 0) then
      WRITE(Un, '(A)') 'Order of constraint states:'      
      call WrLinFile_txt_Table(p_FAST, Un, "Row/Column", LinData%op_z, LinData%names_z  )      
   end if
         
   if (n(Indx_u) > 0) then
      WRITE(Un, '(A)') 'Order of inputs:'   
      call WrLinFile_txt_Table(p_FAST, Un, "Column  ", LinData%op_u, LinData%names_u, rotFrame=LinData%RotFrame_u, UseCol=LinData%use_u  )
   end if
   
   if (n(Indx_u_ext) > 0) then
      WRITE(Un, '(A)') 'Order of extended inputs:'   
      call WrLinFile_txt_Table(p_FAST, Un, "Column  ", LinData%op_u_ext, LinData%names_u_ext )
   end if
   
   if (n(Indx_y) > 0) then
      WRITE(Un, '(A)') 'Order of outputs:'      
      call WrLinFile_txt_Table(p_FAST, Un, "Row  ", LinData%op_y, LinData%names_y, rotFrame=LinData%RotFrame_y, UseCol=LinData%use_y  )      
   end if
      
   !.............
   if (p_FAST%LinOutJac) then
      WRITE (Un,'(/,A,/)' ) 'Jacobian matrices:'    !print a blank line
      ! we'll have the modules write their own Jacobians outside this routine
   end if

         
END SUBROUTINE WrLinFile_txt_Head   
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine that writes the A,B,C,D matrices from linearization to a text file. 
SUBROUTINE WrLinFile_txt_End(Un, p_FAST, LinData)

   INTEGER(IntKi),           INTENT(IN   ) :: Un                  !< unit number
   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< parameters
   TYPE(FAST_LinType),       INTENT(IN   ) :: LinData             !< Linearization data for individual module or glue (coupled system)
   
      ! local variables
   
   WRITE (Un,'(/,A,/)' ) 'Linearized state matrices:'    !print a blank line
   
   ! A matrix
   if (allocated(LinData%A)) call WrPartialMatrix( LinData%A, Un, p_FAST%OutFmt, 'A' )
   ! B matrix   
   if (allocated(LinData%B)) then
      if (allocated(LinData%B_ext)) then
         call WrPartialMatrix( LinData%B, Un, p_FAST%OutFmt, 'B', UseCol=LinData%use_u, ExtCol=LinData%B_ext  )
      else
         call WrPartialMatrix( LinData%B, Un, p_FAST%OutFmt, 'B', UseCol=LinData%use_u )
      end if
   end if
   
   ! C matrix
   if (allocated(LinData%C)) call WrPartialMatrix( LinData%C, Un, p_FAST%OutFmt, 'C', UseRow=LinData%use_y )
   ! D matrix
   if (allocated(LinData%D)) then
      if (allocated(LinData%D_ext)) then
         call WrPartialMatrix( LinData%D, Un, p_FAST%OutFmt, 'D', UseRow=LinData%use_y, UseCol=LinData%use_u, ExtCol=LinData%D_ext )
      else
         call WrPartialMatrix( LinData%D, Un, p_FAST%OutFmt, 'D', UseRow=LinData%use_y, UseCol=LinData%use_u )
      end if
   end if
            

   close(un)
   
END SUBROUTINE WrLinFile_txt_End   
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE WrLinFile_txt_Table(p_FAST, Un, RowCol, op, names, rotFrame, deriv, UseCol)

   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< parameters
   INTEGER(IntKi),           INTENT(IN   ) :: Un                  !< unit number
   CHARACTER(*),             INTENT(IN   ) :: RowCol              !< Row/Column description
   REAL(ReKi),               INTENT(IN   ) :: op(:)               !< operating point values (possibly different size that Desc because of orientations)
   CHARACTER(LinChanLen),    INTENT(IN   ) :: names(:)            !< Descriptions of the channels (names and units)
   logical, optional,        INTENT(IN   ) :: rotFrame(:)         !< determines if this parameter is in the rotating frame
   logical, optional,        intent(in   ) :: deriv               !< flag that tells us if we need to modify the channel names for derivatives (xdot)
   logical, optional,        intent(in   ) :: UseCol(:)           !< flags that tell us if we should use each column or skip it
   
      ! local variables
   INTEGER(IntKi)                          :: TS                  ! tab stop column
   INTEGER(IntKi)                          :: i, i_print          ! loop counter
   INTEGER(IntKi)                          :: i_op                ! loop counter
   
   INTEGER(IntKi)                          :: ErrStat2            ! local error status
   CHARACTER(1024)                         :: ErrMsg2             ! local error message
   logical                                 :: UseDerivNames       !< flag that tells us if we need to modify the channel names for derivatives (xdot)
   logical                                 :: UseThisCol          !< flag that tells us if we should use this particular column or skip it
   logical                                 :: RotatingCol         !< flag that tells us if this column is in the rotating frame
   CHARACTER(*),             PARAMETER     :: RoutineName = 'WrLinFile_txt_Table'    
   CHARACTER(100)                          :: Fmt
   CHARACTER(100)                          :: Fmt_Str
   CHARACTER(100)                          :: FmtOrient
   

   
   if (present(deriv) ) then
      UseDerivNames = deriv
   else
      UseDerivNames = .false.
   end if
   
   
   TS = 14 + 3*p_FAST%FmtWidth+7 ! tab stop after operating point
   
   Fmt       = '(3x,I8,3x,'//trim(p_FAST%OutFmt)//',T'//trim(num2lstr(TS))//',L8,8x,A)'
   FmtOrient = '(3x,I8,3x,'//trim(p_FAST%OutFmt)//',2(", ",'//trim(p_FAST%OutFmt)//'),T'//trim(num2lstr(TS))//',L8,8x,A)'
   Fmt_Str   = '(3x,A10,1x,A,T'//trim(num2lstr(TS))//',A15,1x,A)'
   
   WRITE(Un, Fmt_Str) RowCol,      'Operating Point', 'Rotating Frame?','Description'
   WRITE(Un, Fmt_Str) '----------','---------------', '---------------','-----------'
   
   i_op = 1
   i_print = 1
   do i=1,size(names)
      
      UseThisCol = .true.
      if (present(UseCol)) then
         UseThisCol = useCol(i)
      end if     
      
      RotatingCol = .false.
      if (present(rotFrame)) RotatingCol = rotFrame(i)
                  
      if (index(names(i), ' orientation angle, node ') > 0 ) then  ! make sure this matches what is written in PackMotionMesh_Names()
         if (UseThisCol) then
            WRITE(Un, FmtOrient) i_print, op(i_op), op(i_op+1), op(i_op+2), RotatingCol, trim(names(i))  !//' [OP is a row of the DCM]
            i_print = i_print + 1
         end if
         
         i_op = i_op + 3
      else
         if (UseThisCol) then
            if (UseDerivNames) then
               WRITE(Un, Fmt) i_print, op(i_op), RotatingCol, 'First time derivative of '//trim(names(i))//'/s'
            else
               WRITE(Un, Fmt) i_print, op(i_op), RotatingCol, trim(names(i))
            end if         
            i_print = i_print + 1
         end if         
         
         i_op = i_op + 1
      end if      
   end do
         
   WRITE (Un,'()' )    !print a blank line
   
   
   
END SUBROUTINE WrLinFile_txt_Table   
!----------------------------------------------------------------------------------------------------------------------------------


!> This routine returns the operating points for the entire glue code.
SUBROUTINE Glue_GetOP(p_FAST, y_FAST, ErrStat, ErrMsg)

   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< parameters
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST              !< Output variables for the glue code
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None

   
   INTEGER(IntKi)                          :: i, j                ! loop/temp variables
   INTEGER(IntKi)                          :: ThisModule          ! Module ID # 
   INTEGER(IntKi)                          :: i_u, i_u_ext        ! loop/temp variables
   INTEGER(IntKi)                          :: i_y, i_x            ! loop/temp variables
   INTEGER(IntKi)                          :: ErrStat2            ! local error status
   CHARACTER(1024)                         :: ErrMsg2             ! local error message
   CHARACTER(*),             PARAMETER     :: RoutineName = 'Glue_GetOP'    
   
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   if (.not. allocated(y_FAST%Lin%Glue%op_u)) then  ! assume none of them are allocated
       
         ! calculate the size of the input and output operating points
         ! this size isn't very straightforward since it may contain orientations
      i_u = 0
      i_u_ext = 0
      i_y = 0
      do i = 1,p_FAST%Lin_NumMods
         ThisModule = p_FAST%Lin_ModOrder( i )
         if (allocated(y_FAST%Lin%Modules(ThisModule)%op_u)) then
            i_u = i_u + size(y_FAST%Lin%Modules(ThisModule)%op_u)
         end if
         
         if (allocated(y_FAST%Lin%Modules(ThisModule)%op_u_ext)) then
            i_u_ext = i_u_ext + size(y_FAST%Lin%Modules(ThisModule)%op_u_ext)
         end if
         
         if (allocated(y_FAST%Lin%Modules(ThisModule)%op_y)) then
            i_y = i_y + size(y_FAST%Lin%Modules(ThisModule)%op_y)
         end if
      end do      
      
      call AllocAry( y_FAST%Lin%Glue%op_u, i_u, 'op_u', ErrStat2, ErrMsg2)
         call SetErrStat(errStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      call AllocAry( y_FAST%Lin%Glue%op_u_ext, i_u_ext, 'op_u_ext', ErrStat2, ErrMsg2)
         call SetErrStat(errStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      call AllocAry( y_FAST%Lin%Glue%op_y, i_y, 'op_y', ErrStat2, ErrMsg2)
         call SetErrStat(errStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      call AllocAry( y_FAST%Lin%Glue%op_x, p_FAST%SizeLin(NumModules+1,LIN_ContSTATE_COL), 'op_x', ErrStat2, ErrMsg2)
         call SetErrStat(errStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      call AllocAry( y_FAST%Lin%Glue%op_dx, p_FAST%SizeLin(NumModules+1,LIN_ContSTATE_COL), 'op_dx', ErrStat2, ErrMsg2)
         call SetErrStat(errStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (ErrStat >=AbortErrLev) return
   end if
   
   
   i_u = 1
   i_u_ext = 1
   i_y = 1      
   i_x = 1
   do i = 1,p_FAST%Lin_NumMods
      ThisModule = p_FAST%Lin_ModOrder( i )

      if (allocated(y_FAST%Lin%Modules(ThisModule)%op_u)) then
         do j=1,size(y_FAST%Lin%Modules(ThisModule)%op_u)
            y_FAST%Lin%Glue%op_u(i_u) = y_FAST%Lin%Modules(ThisModule)%op_u(j)
            i_u = i_u + 1;
         end do
      end if
      
      if (allocated(y_FAST%Lin%Modules(ThisModule)%op_u_ext)) then
         do j=1,size(y_FAST%Lin%Modules(ThisModule)%op_u_ext)
            y_FAST%Lin%Glue%op_u_ext(i_u_ext) = y_FAST%Lin%Modules(ThisModule)%op_u_ext(j)
            i_u_ext = i_u_ext + 1;
         end do
      end if      
         
      if (allocated(y_FAST%Lin%Modules(ThisModule)%op_y)) then
         do j=1,size(y_FAST%Lin%Modules(ThisModule)%op_y)
            y_FAST%Lin%Glue%op_y(i_y) = y_FAST%Lin%Modules(ThisModule)%op_y(j)
            i_y = i_y + 1;
         end do      
      end if

      if (allocated(y_FAST%Lin%Modules(ThisModule)%op_x)) then
         do j=1,size(y_FAST%Lin%Modules(ThisModule)%op_x)
            y_FAST%Lin%Glue%op_x(i_x) = y_FAST%Lin%Modules(ThisModule)%op_x(j)
            
            y_FAST%Lin%Glue%op_dx(i_x) = y_FAST%Lin%Modules(ThisModule)%op_dx(j)            
            i_x = i_x + 1;
         end do      
      end if
      
   end do
         
END SUBROUTINE Glue_GetOP
!----------------------------------------------------------------------------------------------------------------------------------

!> This routine forms the Jacobian for the glue-code input-output solves.
SUBROUTINE Glue_Jacobians( t_global, p_FAST, y_FAST, m_FAST, ED, BD, SrvD, AD14, AD, IfW, OpFM, HD, SD, MAPp, FEAM, MD, Orca, &
                         IceF, IceD, MeshMapData, dUdu, dUdy, ErrStat, ErrMsg )

   REAL(DbKi),               INTENT(IN   ) :: t_global            !< current (global) simulation time

   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST              !< Output variables for the glue code
   TYPE(FAST_MiscVarType),   INTENT(INOUT) :: m_FAST              !< Miscellaneous variables
     
   TYPE(ElastoDyn_Data),     INTENT(INOUT) :: ED                  !< ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(INOUT) :: BD                  !< BeamDyn data
   TYPE(ServoDyn_Data),      INTENT(INOUT) :: SrvD                !< ServoDyn data
   TYPE(AeroDyn14_Data),     INTENT(INOUT) :: AD14                !< AeroDyn14 data
   TYPE(AeroDyn_Data),       INTENT(INOUT) :: AD                  !< AeroDyn data
   TYPE(InflowWind_Data),    INTENT(INOUT) :: IfW                 !< InflowWind data
   TYPE(OpenFOAM_Data),      INTENT(INOUT) :: OpFM                !< OpenFOAM data
   TYPE(HydroDyn_Data),      INTENT(INOUT) :: HD                  !< HydroDyn data
   TYPE(SubDyn_Data),        INTENT(INOUT) :: SD                  !< SubDyn data
   TYPE(MAP_Data),           INTENT(INOUT) :: MAPp                !< MAP data
   TYPE(FEAMooring_Data),    INTENT(INOUT) :: FEAM                !< FEAMooring data
   TYPE(MoorDyn_Data),       INTENT(INOUT) :: MD                  !< Data for the MoorDyn module
   TYPE(OrcaFlex_Data),      INTENT(INOUT) :: Orca                !< OrcaFlex interface data
   TYPE(IceFloe_Data),       INTENT(INOUT) :: IceF                !< IceFloe data
   TYPE(IceDyn_Data),        INTENT(INOUT) :: IceD                !< All the IceDyn data used in time-step loop

   TYPE(FAST_ModuleMapType), INTENT(INOUT) :: MeshMapData         !< Data for mapping between modules
   REAL(ReKi), ALLOCATABLE,  INTENT(INOUT) :: dUdu(:,:)           !< Partial derivatives of input-output equations (U(y,u)=0) with respect
                                                                  !!   to the inputs (u)
   REAL(ReKi), ALLOCATABLE,  INTENT(INOUT) :: dUdy(:,:)           !< Partial derivatives of input-output equations (U(y,u)=0) with respect
                                                                  !!   to the outputs (y)
      
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None

   
   
      ! local variables
   INTEGER(IntKi)                          :: ThisModule          ! Module ID
   INTEGER(IntKi)                          :: i                   ! loop counter
   INTEGER(IntKi)                          :: r_start, r_end      ! row start/end of glue matrix
   
   INTEGER(IntKi)                          :: ErrStat2            ! local error status
   CHARACTER(1024)                         :: ErrMsg2             ! local error message
   CHARACTER(*),             PARAMETER     :: RoutineName = 'Glue_Jacobians' 
   
   ErrStat = ErrID_None
   ErrMsg = ""
   
   
   
   
   !.....................................
   ! dUdu 
   !> \f$ \frac{\partial U_\Lambda}{\partial u} =  
   !!  \begin{bmatrix} \frac{\partial U_\Lambda^{IfW}}{\partial u^{IfW}} & \frac{\partial U_\Lambda^{IfW}}{\partial u^{SrvD}} & 
   !!                  \frac{\partial U_\Lambda^{IfW}}{\partial u^{ED}}  & \frac{\partial U_\Lambda^{IfW}}{\partial u^{AD}} \\
   !!                  \frac{\partial U_\Lambda^{SrvD}}{\partial u^{IfW}} & \frac{\partial U_\Lambda^{SrvD}}{\partial u^{SrvD}} & 
   !!                  \frac{\partial U_\Lambda^{SrvD}}{\partial u^{ED}}  & \frac{\partial U_\Lambda^{SrvD}}{\partial u^{AD}} \\
   !!                  \frac{\partial U_\Lambda^{ED}}{\partial u^{IfW}} & \frac{\partial U_\Lambda^{ED}}{\partial u^{SrvD}} & 
   !!                  \frac{\partial U_\Lambda^{ED}}{\partial u^{ED}}  & \frac{\partial U_\Lambda^{ED}}{\partial u^{AD}} \\
   !!                  \frac{\partial U_\Lambda^{AD}}{\partial u^{IfW}} & \frac{\partial U_\Lambda^{AD}}{\partial u^{SrvD}} & 
   !!                  \frac{\partial U_\Lambda^{AD}}{\partial u^{ED}}  & \frac{\partial U_\Lambda^{AD}}{\partial u^{AD}} \\
   !!  \end{bmatrix} = 
   !!  \begin{bmatrix} I & 0 & 0  & \frac{\partial U_\Lambda^{IfW}}{\partial u^{AD}} \\
   !!                  0 & I & 0 & 0 \\
   !!                  0 & 0 & I & \frac{\partial U_\Lambda^{ED}}{\partial u^{AD}} \\
   !!                  0 & 0 & 0 & \frac{\partial U_\Lambda^{AD}}{\partial u^{AD}} \\
   !!  \end{bmatrix} \f$
   !.....................................

   if (.not. allocated(dUdu)) then
      call AllocAry(dUdu, p_FAST%SizeLin(NumModules+1,LIN_INPUT_COL), p_FAST%SizeLin(NumModules+1,LIN_INPUT_COL), 'dUdu', ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat>=AbortErrLev) return
   end if
   
   dUdu = 0.0_ReKi      ! most of this matrix is zero, so we'll just initialize everything and set only the non-zero parts below
   
      !............
      ! \f$ \frac{\partial U_\Lambda^{IfW}}{\partial u^{IfW}} = I \f$ 
      !............
   ThisModule = Module_Ifw
   r_start = p_FAST%LinStartIndx(ThisModule,LIN_INPUT_COL)
   r_end   = r_start + p_FAST%SizeLin(ThisModule,LIN_INPUT_COL) - 1
   do i = r_start,r_end
      dUdu(i,i) = 1.0_ReKi
   end do

      !............
      ! \f$ \frac{\partial U_\Lambda^{SrvD}}{\partial u^{SrvD}} = I \f$ 
      !............
   ThisModule = Module_SrvD
   r_start = p_FAST%LinStartIndx(ThisModule,LIN_INPUT_COL)
   r_end   = r_start + p_FAST%SizeLin(ThisModule,LIN_INPUT_COL) - 1
   do i = r_start,r_end
      dUdu(i,i) = 1.0_ReKi
   end do   
   
      !............
      ! \f$ \frac{\partial U_\Lambda^{ED}}{\partial u^{ED}} = I \f$ 
      !............
   ThisModule = Module_ED
   r_start = p_FAST%LinStartIndx(ThisModule,LIN_INPUT_COL)
   r_end   = r_start + p_FAST%SizeLin(ThisModule,LIN_INPUT_COL) - 1
   do i = r_start,r_end
      dUdu(i,i) = 1.0_ReKi
   end do
   
      !............
      ! \f$ \frac{\partial U_\Lambda^{IfW}}{\partial u^{AD}} \end{bmatrix} = \f$   
      !............
   IF (p_FAST%CompInflow == MODULE_IfW .and. p_FAST%CompAero == MODULE_AD) THEN  
      call Linear_IfW_InputSolve_du_AD( p_FAST, AD%Input(1), dUdu )
   end if ! we're using the InflowWind module
   
      !............
      ! \f$ \frac{\partial U_\Lambda^{ED}}{\partial u^{AD}} \end{bmatrix} = \f$   
      !............   
   IF (p_FAST%CompAero == MODULE_AD) THEN   ! we need to do this regardless of CompElast
      call Linear_ED_InputSolve_du_AD( p_FAST, ED%Input(1), ED%Output(1), AD%y, AD%Input(1), MeshMapData, dUdu, ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   end if ! we're using the InflowWind module
   
      !............
      ! \f$ \frac{\partial U_\Lambda^{AD}}{\partial u^{AD}} \end{bmatrix} = \f$   
      !............
   IF (p_FAST%CompAero == MODULE_AD) THEN 
      call Linear_AD_InputSolve_du_AD( p_FAST, AD%Input(1), ED%Output(1), MeshMapData, dUdu, ErrStat2, ErrMsg2 )         
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   end if ! we're using the InflowWind module
   
   !.....................................
   ! dUdy
   !> \f$ \frac{\partial U_\Lambda}{\partial y} =  
   !!  \begin{bmatrix} \frac{\partial U_\Lambda^{IfW}}{\partial y^{IfW}} & \frac{\partial U_\Lambda^{IfW}}{\partial y^{SrvD}} & 
   !!                  \frac{\partial U_\Lambda^{IfW}}{\partial y^{ED}}  & \frac{\partial U_\Lambda^{IfW}}{\partial y^{AD}} \\
   !!                  \frac{\partial U_\Lambda^{SrvD}}{\partial y^{IfW}} & \frac{\partial U_\Lambda^{SrvD}}{\partial y^{SrvD}} & 
   !!                  \frac{\partial U_\Lambda^{SrvD}}{\partial y^{ED}}  & \frac{\partial U_\Lambda^{SrvD}}{\partial y^{AD}} \\
   !!                  \frac{\partial U_\Lambda^{ED}}{\partial y^{IfW}} & \frac{\partial U_\Lambda^{ED}}{\partial y^{SrvD}} & 
   !!                  \frac{\partial U_\Lambda^{ED}}{\partial y^{ED}}  & \frac{\partial U_\Lambda^{ED}}{\partial y^{AD}} \\
   !!                  \frac{\partial U_\Lambda^{AD}}{\partial y^{IfW}} & \frac{\partial U_\Lambda^{AD}}{\partial y^{SrvD}} & 
   !!                  \frac{\partial U_\Lambda^{AD}}{\partial y^{ED}}  & \frac{\partial U_\Lambda^{AD}}{\partial y^{AD}} \\
   !!  \end{bmatrix} = 
   !!  \begin{bmatrix} 0 & 0 & 0 & 0 \\
   !!                  0 & 0 & \frac{\partial U_\Lambda^{SrvD}}{\partial y^{ED}}  & 0 \\
   !!                  0 & \frac{\partial U_\Lambda^{ED}}{\partial y^{SrvD}} & \frac{\partial U_\Lambda^{ED}}{\partial y^{ED}}  & \frac{\partial U_\Lambda^{ED}}{\partial y^{AD}} \\
   !!                  \frac{\partial U_\Lambda^{AD}}{\partial y^{IfW}} & 0 \frac{\partial U_\Lambda^{AD}}{\partial y^{ED}}  & 0 \\
   !!  \end{bmatrix} \f$
   !.....................................
   if (.not. allocated(dUdy)) then
      call AllocAry(dUdy, p_FAST%SizeLin(NumModules+1,LIN_INPUT_COL), p_FAST%SizeLin(NumModules+1,LIN_OUTPUT_COL), 'dUdy', ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat>=AbortErrLev) return
   end if
         
   dUdy = 0.0_ReKi      ! most of this matrix is zero, so we'll just initialize everything and set only the non-zero parts below
   
      !............
      ! \f$ \frac{\partial U_\Lambda^{AD}}{\partial y^{IfW}} \end{bmatrix} = \f$   
      !............
   if (p_FAST%CompInflow == MODULE_IfW .and. p_FAST%CompAero == MODULE_AD) then   
      call Linear_AD_InputSolve_IfW_dy( p_FAST, AD%Input(1), dUdy )      
   end if
   
   
      !............
      ! \f$ \frac{\partial U_\Lambda^{SrvD}}{\partial y^{ED}} \end{bmatrix} = \f$   
      !............
   if (p_FAST%CompServo == MODULE_SrvD) then   ! need to do this regardless of CompElast
      call Linear_SrvD_InputSolve_dy_ED( p_FAST, y_FAST, dUdy )      
   end if
   
      !............
      ! \f$ \frac{\partial U_\Lambda^{ED}}{\partial y^{SrvD}} \end{bmatrix} = \f$   
      ! \f$ \frac{\partial U_\Lambda^{ED}}{\partial y^{ED}} \end{bmatrix} = \f$   
      ! \f$ \frac{\partial U_\Lambda^{ED}}{\partial y^{AD}} \end{bmatrix} = \f$   
      !............
   call Linear_ED_InputSolve_dy( p_FAST, ED%Input(1), ED%Output(1), AD%y, AD%Input(1), MeshMapData, dUdy, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   
      !............
      ! \f$ \frac{\partial U_\Lambda^{AD}}{\partial y^{ED}} \end{bmatrix} = \f$   
      !............
   if (p_FAST%CompAero == MODULE_AD) then   ! need to do this regardless of CompElast
      call Linear_AD_InputSolve_NoIfW_dy( p_FAST, AD%Input(1), ED%Output(1), MeshMapData, dUdy, ErrStat2, ErrMsg2 )      
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   end if
      
   
   
END SUBROUTINE Glue_Jacobians      
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine forms the dU^{IfW}/du^{AD} block of dUdu.
SUBROUTINE Linear_IfW_InputSolve_du_AD( p_FAST, u_AD, dUdu )

   TYPE(FAST_ParameterType),       INTENT(IN   )   :: p_FAST      !< FAST parameter data 
   TYPE(AD_InputType),             INTENT(IN)      :: u_AD        !< The input meshes (already calculated) from AeroDyn
   real(reki),                     INTENT(INOUT)   :: dUdu(:,:)   !< Jacobian matrix of which we are computing the dU^(IfW)/du^(AD) block
   
   
   INTEGER(IntKi)                          :: i, j, k             ! loop counters
   INTEGER(IntKi)                          :: i2, j2              ! loop counters
   INTEGER(IntKi)                          :: AD_Start_Bl         ! starting index of dUdu (column) where AD blade motion inputs are located
   INTEGER(IntKi)                          :: Node                ! InflowWind node number
   
            
         ! compare with IfW_InputSolve():
   
      Node = 0 !InflowWind node
      if (p_FAST%CompServo == MODULE_SrvD) Node = Node + 1
            
      IF (p_FAST%CompAero == MODULE_AD) THEN 
         
            ! blades:
         AD_Start_Bl = p_FAST%LinStartIndx(MODULE_AD,LIN_INPUT_COL) &
                     + u_AD%TowerMotion%NNodes * 9  & ! 3 fields (MASKID_TRANSLATIONDISP,MASKID_Orientation,MASKID_TRANSLATIONVel) with 3 components
                     + u_AD%HubMotion%NNodes   * 9    ! 3 fields (MASKID_TRANSLATIONDISP,MASKID_Orientation,MASKID_RotationVel) with 3 components
   
         do k = 1,size(u_AD%BladeRootMotion)         
            AD_Start_Bl = AD_Start_Bl + u_AD%BladeRootMotion(k)%NNodes * 3 ! 1 field (MASKID_Orientation) with 3 components
         end do
         ! next is u_AD%BladeMotion(k); note that it has 3 fields and we only need 1
                  
         DO K = 1,SIZE(u_AD%BladeMotion)
            DO J = 1,u_AD%BladeMotion(k)%Nnodes
               Node = Node + 1 ! InflowWind node
               do i=1,3 !XYZ components of this node
                  i2 = p_FAST%LinStartIndx(MODULE_IfW,LIN_INPUT_COL) + (Node-1)*3 + i - 1
                  j2 = AD_Start_Bl + (j-1)*3 + i - 1
                  dUdu( i2, j2 ) = -1.0_ReKi
               end do            
            END DO !J = 1,p%BldNodes ! Loop through the blade nodes / elements
                     
               ! get starting AD index of BladeMotion for next blade
            AD_Start_Bl = AD_Start_Bl + u_AD%BladeMotion(k)%Nnodes * 9  ! 3 fields (MASKID_TRANSLATIONDISP,MASKID_Orientation,MASKID_TRANSLATIONVel) with 3 components
         END DO !K = 1,p%NumBl     
         
            ! tower:
         DO J=1,u_AD%TowerMotion%nnodes
            Node = Node + 1   
            do i=1,3 !XYZ components of this node
               i2 = p_FAST%LinStartIndx(MODULE_IfW,LIN_INPUT_COL) + (Node-1)*3 + i - 1
               j2 = p_FAST%LinStartIndx(MODULE_AD, LIN_INPUT_COL) +    (j-1)*3 + i - 1
               dUdu( i2, j2 ) = -1.0_ReKi
            end do            
         END DO              
         
      END IF     
   
END SUBROUTINE Linear_IfW_InputSolve_du_AD
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine forms the dU^{ED}/du^{AD} block of dUdu.
SUBROUTINE Linear_ED_InputSolve_du_AD( p_FAST, u_ED, y_ED, y_AD, u_AD, MeshMapData, dUdu, ErrStat, ErrMsg )

   TYPE(FAST_ParameterType),       INTENT(IN   )  :: p_FAST         !< Glue-code simulation parameters
   TYPE(ED_InputType),             INTENT(INOUT)  :: u_ED           !< ED Inputs at t
   TYPE(ED_OutputType),            INTENT(IN   )  :: y_ED           !< ElastoDyn outputs (need translation displacement on meshes for loads mapping)
   TYPE(AD_OutputType),            INTENT(IN   )  :: y_AD           !< AeroDyn outputs
   TYPE(AD_InputType),             INTENT(INOUT)  :: u_AD           !< AD inputs (for AD-ED load linerization)
   
   TYPE(FAST_ModuleMapType),       INTENT(INOUT)  :: MeshMapData    !< Data for mapping between modules
   real(reki),                     INTENT(INOUT)  :: dUdu(:,:)      !< Jacobian matrix of which we are computing the dU^(ED)/du^(AD) block
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat        !< Error status
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg         !< Error message
   
      ! local variables
   INTEGER(IntKi)                                 :: i              ! rows/columns
   INTEGER(IntKi)                                 :: J              ! Loops through nodes / elements
   INTEGER(IntKi)                                 :: K              ! Loops through blades
   INTEGER(IntKi)                                 :: AD_Start_Bl    ! starting index of dUdu (column) where AD blade motion inputs are located
   INTEGER(IntKi)                                 :: ED_Start_mt    ! starting index of dUdu (column) where ED blade/tower moment inputs are located
   INTEGER(IntKi)                                 :: ErrStat2       ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                           :: ErrMsg2        ! temporary Error message if ErrStat /= ErrID_None
   
   CHARACTER(*), PARAMETER                        :: RoutineName = 'Linear_ED_InputSolve_du_AD' 
   
   
      ! Initialize error status
      
   ErrStat = ErrID_None
   ErrMsg = ""

      ! ED inputs on blade from AeroDyn
   IF (p_FAST%CompElast == Module_ED) THEN       
      !IF ( p_FAST%CompAero == Module_AD ) THEN !already checked before calling this routine
         
            ! blades:
         AD_Start_Bl = p_FAST%LinStartIndx(Module_AD, LIN_INPUT_COL) &
                     + u_AD%TowerMotion%NNodes * 9  & ! 3 fields (MASKID_TRANSLATIONDISP,MASKID_Orientation,MASKID_TRANSLATIONVel) with 3 components
                     + u_AD%HubMotion%NNodes   * 9    ! 3 fields (MASKID_TRANSLATIONDISP,MASKID_Orientation,MASKID_RotationVel) with 3 components
   
         do k = 1,size(u_AD%BladeRootMotion)         
            AD_Start_Bl = AD_Start_Bl + u_AD%BladeRootMotion(k)%NNodes * 3 ! 1 field (MASKID_Orientation) with 3 components
         end do
         ! next is u_AD%BladeMotion(k); note that it has 3 fields and we only need 1
      
         ED_Start_mt = p_FAST%LinStartIndx(MODULE_ED,LIN_INPUT_COL)
         DO K = 1,SIZE(u_ED%BladePtLoads,1) ! Loop through all blades (p_ED%NumBl)
            ED_Start_mt = ED_Start_mt + u_ED%BladePtLoads(k)%NNodes*3 ! skip the forces on this blade
            
            CALL Linearize_Line2_to_Point( y_AD%BladeLoad(k), u_ED%BladePtLoads(k), MeshMapData%AD_L_2_BDED_B(k), ErrStat2, ErrMsg2, u_AD%BladeMotion(k), y_ED%BladeLn2Mesh(k) )
               CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            
               ! AD is source in the mapping, so we want M_{uSm}               
            if (allocated(MeshMapData%AD_L_2_BDED_B(k)%dM%m_us )) then
               
               do i=1,size(MeshMapData%AD_L_2_BDED_B(k)%dM%m_us,2)
                  do j=1,size(MeshMapData%AD_L_2_BDED_B(k)%dM%m_us,1)
                     dUdu( ED_Start_mt + j - 1, AD_Start_Bl + i - 1 ) = - MeshMapData%AD_L_2_BDED_B(k)%dM%m_us(j,i)
                  end do
               end do
               
            end if
            
               ! get starting index of next blade
            AD_Start_Bl = AD_Start_Bl + u_AD%BladeMotion(k)%Nnodes * 9  ! 3 fields (MASKID_TRANSLATIONDISP,MASKID_Orientation,MASKID_TRANSLATIONVel) with 3 components                         
            ED_Start_mt = ED_Start_mt + u_ED%BladePtLoads(k)%NNodes* 3  ! skip the moments on this blade
               
         END DO
                           
      !END IF      
   END IF
      
   !IF ( p_FAST%CompAero == Module_AD ) THEN    !already checked before calling this routine  
      IF ( y_AD%TowerLoad%Committed ) THEN
         ED_Start_mt = p_FAST%LinStartIndx(MODULE_ED,LIN_INPUT_COL)           
         if (allocated(u_ED%BladePtLoads)) then
            do i=1,size(u_ED%BladePtLoads)
               ED_Start_mt = ED_Start_mt + u_ED%BladePtLoads(i)%NNodes * 6  ! 3 forces + 3 moments at each node on each blade
            end do      
         end if
         ED_Start_mt = ED_Start_mt + u_ED%PlatformPtMesh%NNodes * 6      &      ! 3 forces + 3 moments at each node
                                   + u_ED%TowerPtLoads%NNodes   * 3             ! 3 forces at each node (we're going to start at the moments)
         
         
         CALL Linearize_Line2_to_Point( y_AD%TowerLoad, u_ED%TowerPtLoads, MeshMapData%AD_L_2_ED_P_T, ErrStat2, ErrMsg2, u_AD%TowerMotion, y_ED%TowerLn2Mesh )
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)  
            
            ! AD is source in the mapping, so we want M_{uSm}               
         if (allocated(MeshMapData%AD_L_2_ED_P_T%dM%m_us )) then
               
            do i=1,size(MeshMapData%AD_L_2_ED_P_T%dM%m_us,2)
               do j=1,size(MeshMapData%AD_L_2_ED_P_T%dM%m_us,1)                     
                  dUdu( ED_Start_mt + j - 1, p_FAST%LinStartIndx(MODULE_AD,LIN_INPUT_COL) + i - 1 ) = - MeshMapData%AD_L_2_ED_P_T%dM%m_us(j,i)
               end do
            end do
               
         end if                                    
      END IF
            
   !END IF
               
END SUBROUTINE Linear_ED_InputSolve_du_AD
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine forms the dU^{AD}/du^{AD} block of dUdu.
SUBROUTINE Linear_AD_InputSolve_du_AD( p_FAST, u_AD, y_ED, MeshMapData, dUdu, ErrStat, ErrMsg )

      ! Passed variables
   TYPE(FAST_ParameterType),    INTENT(IN   )   :: p_FAST      !< FAST parameter data    
   TYPE(AD_InputType),          INTENT(INOUT)   :: u_AD        !< The inputs to AeroDyn14
   TYPE(ED_OutputType),         INTENT(IN)      :: y_ED        !< The outputs from the structural dynamics module
   TYPE(FAST_ModuleMapType),    INTENT(INOUT)   :: MeshMapData !< Data for mapping between modules
   real(reki),                  INTENT(INOUT)   :: dUdu(:,:)   !< Jacobian matrix of which we are computing the dU^(ED)/du^(AD) block
   
   INTEGER(IntKi)                               :: ErrStat     !< Error status of the operation
   CHARACTER(*)                                 :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      ! Local variables:

   INTEGER(IntKi)                               :: I              ! Loops through rows/columns
   INTEGER(IntKi)                               :: J              ! Loops through nodes / elements
   INTEGER(IntKi)                               :: K              ! Loops through blades   
   INTEGER(IntKi)                               :: AD_Start_td    ! starting index of dUdu (column) where AD translation displacements are located   
   INTEGER(IntKi)                               :: AD_Start_tv    ! starting index of dUdu (column) where AD translation velocities are located   
   INTEGER(IntKi)                               :: ErrStat2
   CHARACTER(ErrMsgLen)                         :: ErrMsg2 
   CHARACTER(*), PARAMETER                      :: RoutineName = 'Linear_AD_InputSolve_du_AD'

   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   ! note that we assume this block matrix has been initialized to zero before calling this routine
   
   ! first we set this block to have 1s on the diagonal:   
   do j=p_FAST%LinStartIndx(MODULE_AD, LIN_INPUT_COL), &
        p_FAST%LinStartIndx(MODULE_AD, LIN_INPUT_COL)+p_FAST%SizeLin(MODULE_AD, LIN_INPUT_COL) - 1
      dUdu(j,j) = 1.0_ReKi
   end do
      
   ! then we look at how the translational displacement gets transfered to the translational velocity: 
      
   !-------------------------------------------------------------------------------------------------
   ! Set the inputs from ElastoDyn and/or BeamDyn:
   !-------------------------------------------------------------------------------------------------
      
      ! tower
   IF (u_AD%TowerMotion%Committed) THEN
            
      CALL Linearize_Line2_to_Line2( y_ED%TowerLn2Mesh, u_AD%TowerMotion, MeshMapData%ED_L_2_AD_L_T, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//':u_AD%TowerMotion' )     

      AD_Start_td = p_FAST%LinStartIndx(MODULE_AD, LIN_INPUT_COL)   
      AD_Start_tv = AD_Start_td + u_AD%TowerMotion%NNodes * 6 ! 2 fields (TranslationDisp and Orientation) with 3 components before translational velocity field      

      !AD is the destination here, so we need tv_ud
      if (allocated( MeshMapData%ED_L_2_AD_L_T%dM%tv_ud)) then
         
         do i=1,size(MeshMapData%ED_L_2_AD_L_T%dM%tv_ud,2)
            do j=1,size(MeshMapData%ED_L_2_AD_L_T%dM%tv_ud,1)
               dUdu( AD_Start_tv + j - 1, AD_Start_td + i - 1 ) = - MeshMapData%ED_L_2_AD_L_T%dM%tv_ud(j,i)
            end do
         end do
                  
      end if
               
            
   END IF
   
      
   
      ! blades
   IF (p_FAST%CompElast == Module_ED ) THEN
      
      AD_Start_td = p_FAST%LinStartIndx(MODULE_AD, LIN_INPUT_COL) &
                  + u_AD%TowerMotion%NNodes * 9  & ! 3 fields (MASKID_TRANSLATIONDISP,MASKID_Orientation,MASKID_TRANSLATIONVel) with 3 components
                  + u_AD%HubMotion%NNodes   * 9    ! 3 fields (MASKID_TRANSLATIONDISP,MASKID_Orientation,MASKID_RotationVel) with 3 components
   
      do k = 1,size(u_AD%BladeRootMotion)         
         AD_Start_td = AD_Start_td + u_AD%BladeRootMotion(k)%NNodes * 3 ! 1 field (MASKID_Orientation) with 3 components
      end do
      
      
      
      DO k=1,size(y_ED%BladeLn2Mesh)
         CALL Linearize_Line2_to_Line2( y_ED%BladeLn2Mesh(k), u_AD%BladeMotion(k), MeshMapData%BDED_L_2_AD_L_B(k), ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//':u_AD%BladeMotion('//trim(num2lstr(k))//')' )   
            
         AD_Start_tv = AD_Start_td + u_AD%BladeMotion(k)%NNodes * 6 ! 2 fields (TranslationDisp and Orientation) with 3 components before translational velocity field      
            
         !AD is the destination here, so we need tv_ud
         if (allocated( MeshMapData%BDED_L_2_AD_L_B(k)%dM%tv_ud)) then
         
            do i=1,size(MeshMapData%BDED_L_2_AD_L_B(k)%dM%tv_ud,2)
               do j=1,size(MeshMapData%BDED_L_2_AD_L_B(k)%dM%tv_ud,1)
                  dUdu( AD_Start_tv + j - 1, AD_Start_td + i - 1 ) = - MeshMapData%BDED_L_2_AD_L_B(k)%dM%tv_ud(j,i)
               end do
            end do
                  
         end if        
         
         AD_Start_td = AD_Start_td + u_AD%BladeMotion(k)%NNodes * 6 ! 3 fields (TranslationDisp, Orientation, TranslationVel) with 3 components
            
      END DO
      
   !ELSEIF (p_FAST%CompElast == Module_BD ) THEN
                  
   END IF
   
   
END SUBROUTINE Linear_AD_InputSolve_du_AD
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine forms the dU^{SrvD}/dy^{ED} block of dUdy.
SUBROUTINE Linear_SrvD_InputSolve_dy_ED( p_FAST, y_FAST, dUdy  )
!..................................................................................................................................

   TYPE(FAST_ParameterType),       INTENT(IN)     :: p_FAST         !< Glue-code simulation parameters
   TYPE(FAST_OutputFileType),      INTENT(IN)     :: y_FAST         !< Output variables for the glue code
   real(reki),                     INTENT(INOUT)  :: dUdy(:,:)      !< Jacobian matrix of which we are computing the dU^{SrvD}/dy^{ED} block
   
   integer(intKi)                                 :: ED_Start_Yaw   !< starting index of dUdy (column) where ED Yaw/YawRate/HSS_Spd outputs are located (just before WriteOutput)

   
   
   INTEGER(IntKi)                                   :: i            ! loop counter
   
   INTEGER(IntKi)                                   :: ErrStat2                 ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                             :: ErrMsg2                  ! temporary Error message if ErrStat /= ErrID_None
   CHARACTER(*), PARAMETER                          :: RoutineName = 'Linear_SrvD_InputSolve_dy_ED' 
   
   ED_Start_Yaw = p_FAST%LinStartIndx(MODULE_ED,LIN_OUTPUT_COL) + p_FAST%SizeLin(Module_ED,LIN_OUTPUT_COL) & !end of ED outputs (+1)
                  - y_FAST%numOuts(Module_ED) - 3 ! start of ED where Yaw, YawRate, HSS_Spd occur (right before WriteOutputs)
   do i=1,3
      dUdy(p_FAST%LinStartIndx(MODULE_SrvD,LIN_INPUT_COL) + i - 1, ED_Start_Yaw + i - 1) = -1.0_ReKi
   end do
      
   !IF (u_SrvD%NTMD%Mesh%Committed) THEN
   !   
   !   CALL Linearize_Point_to_Point( y_ED%NacelleMotion, u_SrvD%NTMD%Mesh, MeshMapData%ED_P_2_SrvD_P_N, ErrStat2, ErrMsg2 )
   !      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   !         
   !END IF
   !
   !IF (u_SrvD%TTMD%Mesh%Committed) THEN
   !   
   !   CALL Linearize_Line2_to_Point( y_ED%TowerLn2Mesh, u_SrvD%TTMD%Mesh, MeshMapData%ED_L_2_SrvD_P_T, ErrStat2, ErrMsg2 )
   !      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   !         
   !END IF
                                    
END SUBROUTINE Linear_SrvD_InputSolve_dy_ED
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine forms the dU^{ED}/dy^{SrvD}, dU^{ED}/dy^{ED}, and dU^{ED}/dy^{AD} blocks of dUdy.
SUBROUTINE Linear_ED_InputSolve_dy( p_FAST, u_ED, y_ED, y_AD, u_AD, MeshMapData, dUdy, ErrStat, ErrMsg )

   TYPE(FAST_ParameterType),       INTENT(IN   )  :: p_FAST           !< Glue-code simulation parameters
   TYPE(ED_InputType),             INTENT(INOUT)  :: u_ED             !< ED Inputs at t
   TYPE(ED_OutputType),            INTENT(IN   )  :: y_ED             !< ElastoDyn outputs (need translation displacement on meshes for loads mapping)
   TYPE(AD_OutputType),            INTENT(IN   )  :: y_AD             !< AeroDyn outputs
   TYPE(AD_InputType),             INTENT(INOUT)  :: u_AD             !< AD inputs (for AD-ED load linerization)
                                                                      
   TYPE(FAST_ModuleMapType),       INTENT(INOUT)  :: MeshMapData      !< Data for mapping between modules
   real(reki),                     INTENT(INOUT)  :: dUdy(:,:)        !< Jacobian matrix of which we are computing the dU^(ED)/du^(AD) block
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat          !< Error status
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg           !< Error message
   
      ! local variables
   INTEGER(IntKi)                                 :: i                ! rows/columns
   INTEGER(IntKi)                                 :: J                ! Loops through nodes / elements
   INTEGER(IntKi)                                 :: K                ! Loops through blades
   INTEGER(IntKi)                                 :: AD_Start_tmp     ! starting index of dUdy (column) where particular AD fields are located
   INTEGER(IntKi)                                 :: ED_Start_tmp     ! starting index of dUdy (row/column) where particular ED fields are located
   INTEGER(IntKi)                                 :: ED_Out_Start_tmp ! starting index of dUdy (row/column) where ED translation displacement fields are located
   INTEGER(IntKi)                                 :: ErrStat2         ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                           :: ErrMsg2          ! temporary Error message if ErrStat /= ErrID_None
   
   CHARACTER(*), PARAMETER                        :: RoutineName = 'Linear_ED_InputSolve_du_AD' 
   
   
      ! Initialize error status
      
   ErrStat = ErrID_None
   ErrMsg = ""
   
      ! ED inputs from ServoDyn outputs
   IF ( p_FAST%CompServo == Module_SrvD ) THEN

         ! BlPitchCom, YawMom, GenTrq
      ED_Start_tmp = p_FAST%LinStartIndx(MODULE_ED,LIN_INPUT_COL) + p_FAST%SizeLin(MODULE_ED,LIN_INPUT_COL) - size(u_ED%BlPitchCom) - 2
      do i=1,size(u_ED%BlPitchCom)+2
         dUdy(ED_Start_tmp + i - 1, p_FAST%LinStartIndx(Module_SrvD,LIN_OUTPUT_COL) + i - 1) = -1.0_ReKi
      end do
      
      
      !IF (y_SrvD%NTMD%Mesh%Committed) THEN      
      !   CALL Linearize_Point_to_Point( y_SrvD%NTMD%Mesh, u_ED%NacelleLoads, MeshMapData%SrvD_P_2_ED_P_N, ErrStat2, ErrMsg2, u_SrvD%NTMD%Mesh, y_ED%NacelleMotion )
      !      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg,RoutineName//':u_ED%NacelleLoads' )                  
      !END IF
      !
      !IF (y_SrvD%TTMD%Mesh%Committed) THEN      
      !   CALL Linearize_Point_to_Point( y_SrvD%TTMD%Mesh, u_ED%TowerPtLoads, MeshMapData%SrvD_P_2_ED_P_T, ErrStat2, ErrMsg2, u_SrvD%TTMD%Mesh, y_ED%TowerLn2Mesh )
      !      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg,RoutineName//':u_ED%TowerPtLoads' )                  
      !END IF
            
   END IF

            
      ! ED inputs on blade from AeroDyn
   IF (p_FAST%CompElast == Module_ED) THEN 
      IF ( p_FAST%CompAero == Module_AD ) THEN
         
         AD_Start_tmp = p_FAST%LinStartIndx(MODULE_AD,LIN_OUTPUT_COL) + y_AD%TowerLoad%NNodes * 6    ! 2 fields (force, moment) with 3 components
         ED_Start_tmp = p_FAST%LinStartIndx(MODULE_ED,LIN_INPUT_COL) ! blade loads in u_ED
         ED_Out_Start_tmp  = p_FAST%LinStartIndx(MODULE_ED,LIN_OUTPUT_COL) ! blade motions in y_ED
         
         DO K = 1,SIZE(u_ED%BladePtLoads,1) ! Loop through all blades (p_ED%NumBl)
            CALL Linearize_Line2_to_Point( y_AD%BladeLoad(k), u_ED%BladePtLoads(k), MeshMapData%AD_L_2_BDED_B(k), ErrStat2, ErrMsg2, u_AD%BladeMotion(k), y_ED%BladeLn2Mesh(k) )
               CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
               if (errStat >=AbortErrLev) return
               
            ! force equation:
               
               ! force-to-force transfer:
            do i=1,size(MeshMapData%AD_L_2_BDED_B(k)%dM%li,2)
               do j=1,size(MeshMapData%AD_L_2_BDED_B(k)%dM%li,1)
                  dUdy( ED_Start_tmp + j - 1, AD_Start_tmp + i - 1 ) = - MeshMapData%AD_L_2_BDED_B(k)%dM%li(j,i)
               end do
            end do                              

            ! moment equation:
            
            ED_Start_tmp = ED_Start_tmp + u_ED%BladePtLoads(k)%NNodes*3 ! skip the ED forces
                         
               ! translation displacement-to-moment transfer:
            do i=1,size(MeshMapData%AD_L_2_BDED_B(k)%dM%m_uD,2)
               do j=1,size(MeshMapData%AD_L_2_BDED_B(k)%dM%m_uD,1)
                  dUdy( ED_Start_tmp + j - 1, ED_Out_Start_tmp + i - 1 ) = - MeshMapData%AD_L_2_BDED_B(k)%dM%m_uD(j,i)
               end do
            end do            
                        
               ! AD force-to-ED moment transfer:
            do i=1,size(MeshMapData%AD_L_2_BDED_B(k)%dM%m_f,2)
               do j=1,size(MeshMapData%AD_L_2_BDED_B(k)%dM%m_f,1)
                  dUdy( ED_Start_tmp + j - 1, AD_Start_tmp + i - 1 ) = - MeshMapData%AD_L_2_BDED_B(k)%dM%m_f(j,i)
               end do
            end do            
            
            AD_Start_tmp = AD_Start_tmp + y_AD%BladeLoad(k)%NNodes*3 ! skip the AD forces to get to the moments
               
               ! moment-to-moment transfer:
            do i=1,size(MeshMapData%AD_L_2_BDED_B(k)%dM%li,2)
               do j=1,size(MeshMapData%AD_L_2_BDED_B(k)%dM%li,1)
                  dUdy( ED_Start_tmp + j - 1, AD_Start_tmp + i - 1 ) = - MeshMapData%AD_L_2_BDED_B(k)%dM%li(j,i)
               end do
            end do                              
                           
            
            AD_Start_tmp = AD_Start_tmp + y_AD%BladeLoad(k)%NNodes*3       ! skip the moments to get to forces on next blade
            ED_Start_tmp = ED_Start_tmp + u_ED%BladePtLoads(k)%NNodes*3    ! skip the moments to get to forces on next blade
            ED_Out_Start_tmp  = ED_Out_Start_tmp  + y_ED%BladeLn2Mesh(k)%NNodes*18   ! skip 6 fields (with 3 components) to get to next blade  
                                    
         END DO
      END IF
      
   END IF
                        
   IF ( p_FAST%CompAero == Module_AD ) THEN
      
      IF ( y_AD%TowerLoad%Committed ) THEN
         CALL Linearize_Line2_to_Point( y_AD%TowerLoad, u_ED%TowerPtLoads, MeshMapData%AD_L_2_ED_P_T, ErrStat2, ErrMsg2, u_AD%TowerMotion, y_ED%TowerLn2Mesh )
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)         
            
            
         AD_Start_tmp = p_FAST%LinStartIndx(MODULE_AD,LIN_OUTPUT_COL) ! tower is first in y_AD
         
         ! find tower in u_ED:
         ED_Start_tmp = p_FAST%LinStartIndx(MODULE_ED,LIN_INPUT_COL) 
         if (allocated(u_ED%BladePtLoads)) then
            do K = 1,SIZE(u_ED%BladePtLoads,1) ! Loop through all blades (p_ED%NumBl)
               ED_Start_tmp = ED_Start_tmp + u_ED%BladePtLoads(k)%NNodes*6 ! 2 fields (force/moment) with 3 components on each blade
            end do
         end if
         ED_Start_tmp = ED_Start_tmp + u_ED%PlatformPtMesh%NNodes*6 ! 2 fields (force/moment) with 3 components
         
         ! find tower in y_ED:
         ED_Out_Start_tmp  = p_FAST%LinStartIndx(MODULE_ED,LIN_OUTPUT_COL) 
         if (allocated(y_ED%BladeLn2Mesh)) then
            do K = 1,SIZE(y_ED%BladeLn2Mesh,1) ! Loop through all blades (p_ED%NumBl)
               ED_Out_Start_tmp = ED_Out_Start_tmp + y_ED%BladeLn2Mesh(k)%NNodes*18 ! 6 fields with 3 components on each blade
            end do
         end if
         ED_Out_Start_tmp = ED_Out_Start_tmp + y_ED%PlatformPtMesh%NNodes*6 ! 6 fields with 3 components
            
            
            ! force equation:
               
            ! force-to-force transfer:
         do i=1,size(MeshMapData%AD_L_2_ED_P_T%dM%li,2)
            do j=1,size(MeshMapData%AD_L_2_ED_P_T%dM%li,1)
               dUdy( ED_Start_tmp + j - 1, AD_Start_tmp + i - 1 ) = - MeshMapData%AD_L_2_ED_P_T%dM%li(j,i)
            end do
         end do                              

         ! moment equation:
            
         ED_Start_tmp = ED_Start_tmp + u_ED%TowerPtLoads%NNodes*3 ! skip the ED forces to get to the moments
                         
            ! translation displacement-to-moment transfer:
         do i=1,size(MeshMapData%AD_L_2_ED_P_T%dM%m_uD,2)
            do j=1,size(MeshMapData%AD_L_2_ED_P_T%dM%m_uD,1)
               dUdy( ED_Start_tmp + j - 1, ED_Out_Start_tmp + i - 1 ) = - MeshMapData%AD_L_2_ED_P_T%dM%m_uD(j,i)
            end do
         end do            
                        
            ! AD force-to-ED moment transfer:
         do i=1,size(MeshMapData%AD_L_2_ED_P_T%dM%m_f,2)
            do j=1,size(MeshMapData%AD_L_2_ED_P_T%dM%m_f,1)
               dUdy( ED_Start_tmp + j - 1, AD_Start_tmp + i - 1 ) = - MeshMapData%AD_L_2_ED_P_T%dM%m_f(j,i)
            end do
         end do            
            
         AD_Start_tmp = AD_Start_tmp + y_AD%TowerLoad%NNodes*3 ! skip the AD forces to get to the moments
               
            ! moment-to-moment transfer:
         do i=1,size(MeshMapData%AD_L_2_ED_P_T%dM%li,2)
            do j=1,size(MeshMapData%AD_L_2_ED_P_T%dM%li,1)
               dUdy( ED_Start_tmp + j - 1, AD_Start_tmp + i - 1 ) = - MeshMapData%AD_L_2_ED_P_T%dM%li(j,i)
            end do
         end do                   
            
      END IF
            
   END IF
                           
END SUBROUTINE Linear_ED_InputSolve_dy
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine forms the dU^{AD}/dy^{IfW} block of dUdy.
SUBROUTINE Linear_AD_InputSolve_IfW_dy( p_FAST, u_AD, dUdy )

      ! Passed variables
   TYPE(FAST_ParameterType),       INTENT(IN   )  :: p_FAST         !< FAST parameter data    
   TYPE(AD_InputType),             INTENT(INOUT)  :: u_AD           !< The inputs to AeroDyn
   real(reki),                     INTENT(INOUT)  :: dUdy(:,:)      !< Jacobian matrix of which we are computing the dU^{AD}/dy^{IfW} block

      ! Local variables:

   INTEGER(IntKi)                               :: I           ! Loops through components
   INTEGER(IntKi)                               :: J           ! Loops through nodes / elements
   INTEGER(IntKi)                               :: K           ! Loops through blades
   INTEGER(IntKi)                               :: NodeNum     ! Node number for blade/node on mesh
   INTEGER(IntKi)                               :: NNodes
   INTEGER(IntKi)                               :: node
   INTEGER(IntKi)                               :: AD_Start_tmp   ! starting index of dUdy (row) where AD input equations (for specific fields) are located   

                  
   !-------------------------------------------------------------------------------------------------
   ! Set the inputs from inflow wind:
   !-------------------------------------------------------------------------------------------------
   !IF (p_FAST%CompInflow == MODULE_IfW) THEN !already checked in calling routine

      if (p_FAST%CompServo == MODULE_SrvD) then
         node = 2
      else
         node = 1
      end if
      
      
      AD_Start_tmp = p_FAST%LinStartIndx(MODULE_AD, LIN_INPUT_COL) &
                   + u_AD%TowerMotion%NNodes * 9  & ! 3 fields (MASKID_TRANSLATIONDISP,MASKID_Orientation,MASKID_TRANSLATIONVel) with 3 components
                   + u_AD%HubMotion%NNodes   * 9    ! 3 fields (MASKID_TRANSLATIONDISP,MASKID_Orientation,MASKID_RotationVel) with 3 components
   
      do k = 1,size(u_AD%BladeRootMotion)         
         AD_Start_tmp = AD_Start_tmp + u_AD%BladeRootMotion(k)%NNodes * 3 ! 1 field (MASKID_Orientation) with 3 components
      end do                  
      DO k=1,size(u_AD%BladeMotion)         
         AD_Start_tmp = AD_Start_tmp + u_AD%BladeMotion(k)%NNodes * 6 ! 3 fields (TranslationDisp, Orientation, TranslationVel) with 3 components            
      END DO
                  
      
      do k=1,size(u_AD%InflowOnBlade,3) ! blades
         do j=1,size(u_AD%InflowOnBlade,2) ! nodes
            do i=1,3 !velocity component
               dUdy( AD_Start_tmp + i - 1, p_FAST%LinStartIndx(MODULE_IfW,LIN_OUTPUT_COL) + (node-1)*3 + i - 1 ) = -1.0_ReKi                
            end do
            node = node + 1
            AD_Start_tmp = AD_Start_tmp + 3
         end do         
      end do
                  
      if ( allocated(u_AD%InflowOnTower) ) then         
         do j=1,size(u_AD%InflowOnTower,2) !nodes
            do i=1,3 !velocity component
               dUdy( AD_Start_tmp + i - 1, p_FAST%LinStartIndx(MODULE_IfW,LIN_OUTPUT_COL) + (node-1)*3 + i - 1 ) = -1.0_ReKi                
            end do
            node = node + 1
            AD_Start_tmp = AD_Start_tmp + 3
         end do      
      end if
                     
   !END IF
   
   
END SUBROUTINE Linear_AD_InputSolve_IfW_dy
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine forms the dU^{AD}/dy^{ED} block of dUdy.
SUBROUTINE Linear_AD_InputSolve_NoIfW_dy( p_FAST, u_AD, y_ED, MeshMapData, dUdy, ErrStat, ErrMsg )

      ! Passed variables
   TYPE(FAST_ParameterType),    INTENT(IN   )   :: p_FAST      !< FAST parameter data    
   TYPE(AD_InputType),          INTENT(INOUT)   :: u_AD        !< The inputs to AeroDyn14
   TYPE(ED_OutputType),         INTENT(IN)      :: y_ED        !< The outputs from the structural dynamics module
   TYPE(FAST_ModuleMapType),    INTENT(INOUT)   :: MeshMapData !< Data for mapping between modules
   real(reki),                  INTENT(INOUT)   :: dUdy(:,:)   !< Jacobian matrix of which we are computing the dU^{AD}/dy^{ED} block
   
   INTEGER(IntKi)                               :: ErrStat     !< Error status of the operation
   CHARACTER(*)                                 :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      ! Local variables:

   INTEGER(IntKi)                               :: i,j         ! Loops through rows/columns of mesh-mapping linearization matrices
   INTEGER(IntKi)                               :: K           ! Loops through blades
   INTEGER(IntKi)                               :: NodeNum     ! Node number for blade/node on mesh
   INTEGER(IntKi)                               :: node
   INTEGER(IntKi)                               :: AD_Start_tmp     ! starting index of dUdy (column) where particular AD fields are located
   INTEGER(IntKi)                               :: ED_Start_tmp     ! starting index of dUdy (row/column) where particular ED fields are located
   INTEGER(IntKi)                               :: ErrStat2
   CHARACTER(ErrMsgLen)                         :: ErrMsg2 
   CHARACTER(*), PARAMETER                      :: RoutineName = 'Linear_AD_InputSolve_NoIfW_dy'

   
   ErrStat = ErrID_None
   ErrMsg  = ""
               
   
   !-------------------------------------------------------------------------------------------------
   ! Set the inputs from ElastoDyn and/or BeamDyn:
   !-------------------------------------------------------------------------------------------------
      !...................................
      ! tower
      !...................................
   AD_Start_tmp = p_FAST%LinStartIndx(MODULE_AD, LIN_INPUT_COL)
   
      ! tower comes after ED blades and platform:
   ED_Start_tmp = p_FAST%LinStartIndx(Module_ED, LIN_OUTPUT_COL)
   if (allocated(y_ED%BladeLn2Mesh)) then
      DO k=1,size(y_ED%BladeLn2Mesh)
         ED_Start_tmp = ED_Start_tmp + y_ED%BladeLn2Mesh(k)%NNodes * 18 ! 6 fields (translation disp, orientation, translation vel, rotation vel, translation acc, rotation acc) with 3 components            
      END DO      
   end if
   ED_Start_tmp = ED_Start_tmp + y_ED%PlatformPtMesh%NNodes * 18 ! 6 fields (translation disp, orientation, translation vel, rotation vel, translation acc, rotation acc) with 3 components
   
   IF (u_AD%TowerMotion%Committed) THEN
            
      CALL Linearize_Line2_to_Line2( y_ED%TowerLn2Mesh, u_AD%TowerMotion, MeshMapData%ED_L_2_AD_L_T, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//':u_AD%TowerMotion' )      
         if (errStat>=AbortErrLev) return
         
      ! *** AD translational displacement: from ED translational displacement (MeshMapData%ED_L_2_AD_L_T%dM%mi) and orientation (MeshMapData%ED_L_2_AD_L_T%dM%fx_p)                              
      do i=1,size(MeshMapData%ED_L_2_AD_L_T%dM%mi,2)
         do j=1,size(MeshMapData%ED_L_2_AD_L_T%dM%mi,1)
            dUdy( AD_Start_tmp + j - 1, ED_Start_tmp + i - 1 ) = - MeshMapData%ED_L_2_AD_L_T%dM%mi(j,i)
         end do
      end do                              
      ED_Start_tmp = ED_Start_tmp + y_ED%TowerLn2Mesh%NNodes * 3 ! move past the ED translation disp field to orientation field
            
      do i=1,size(MeshMapData%ED_L_2_AD_L_T%dM%fx_p,2)
         do j=1,size(MeshMapData%ED_L_2_AD_L_T%dM%fx_p,1)
            dUdy( AD_Start_tmp + j - 1, ED_Start_tmp + i - 1 ) = - MeshMapData%ED_L_2_AD_L_T%dM%fx_p(j,i)
         end do
      end do
         
      ! *** AD orientation: from ED orientation
      AD_Start_tmp = AD_Start_tmp + u_AD%TowerMotion%NNodes * 3 ! move past the AD translation disp field to orientation field         
      do i=1,size(MeshMapData%ED_L_2_AD_L_T%dM%mi,2)
         do j=1,size(MeshMapData%ED_L_2_AD_L_T%dM%mi,1)
            dUdy( AD_Start_tmp + j - 1, ED_Start_tmp + i - 1 ) = - MeshMapData%ED_L_2_AD_L_T%dM%mi(j,i)
         end do
      end do 
         
      ! *** AD translational velocity: from ED translation displacement, translation velocity, and rotational velocity         
      AD_Start_tmp = AD_Start_tmp + u_AD%TowerMotion%NNodes * 3 ! move past the AD orientation field to translation velocity field              
      ED_Start_tmp = ED_Start_tmp - y_ED%TowerLn2Mesh%NNodes * 3 ! move past the ED translation disp field back to beginning of ED tower         
      do i=1,size(MeshMapData%ED_L_2_AD_L_T%dM%tv_uS,2)
         do j=1,size(MeshMapData%ED_L_2_AD_L_T%dM%tv_uS,1)
            dUdy( AD_Start_tmp + j - 1, ED_Start_tmp + i - 1 ) = - MeshMapData%ED_L_2_AD_L_T%dM%tv_uS(j,i)
         end do
      end do                              

      ED_Start_tmp = ED_Start_tmp + y_ED%TowerLn2Mesh%NNodes * 6 ! move past the ED translation disp and orientation fields to translation vel field
      do i=1,size(MeshMapData%ED_L_2_AD_L_T%dM%mi,2)
         do j=1,size(MeshMapData%ED_L_2_AD_L_T%dM%mi,1)
            dUdy( AD_Start_tmp + j - 1, ED_Start_tmp + i - 1 ) = - MeshMapData%ED_L_2_AD_L_T%dM%mi(j,i)
         end do
      end do 
         
      ED_Start_tmp = ED_Start_tmp + y_ED%TowerLn2Mesh%NNodes * 3 ! move past the ED translational velocity field to rotational velocity field
      do i=1,size(MeshMapData%ED_L_2_AD_L_T%dM%fx_p,2)
         do j=1,size(MeshMapData%ED_L_2_AD_L_T%dM%fx_p,1)
            dUdy( AD_Start_tmp + j - 1, ED_Start_tmp + i - 1 ) = - MeshMapData%ED_L_2_AD_L_T%dM%fx_p(j,i)
         end do
      end do
         
         
      AD_Start_tmp = AD_Start_tmp + u_AD%TowerMotion%NNodes * 3 ! 1 [of 3] fields ([already took care of translation disp, orientation], translation vel) with 3 components
      ED_Start_tmp = ED_Start_tmp + y_ED%TowerLn2Mesh%NNodes * 9 ! 3 [of 6] fields ([already took care of translation disp, orientation, translation vel], rotation vel, translation acc, rotation acc) with 3 components
   ELSE                                          
      ! This should be zero, but we don't have conditional statements elsewhere, so just in case there are nodes on an uncommitted mesh....
      AD_Start_tmp = AD_Start_tmp + u_AD%TowerMotion%NNodes  * 9  ! 3 fields (translation disp, orientation, translation vel) with 3 components
      ED_Start_tmp = ED_Start_tmp + y_ED%TowerLn2Mesh%NNodes * 18 ! 6 fields (translation disp, orientation, translation vel, rotation vel, translation acc, rotation acc) with 3 components            
   END IF
      
      !...................................
      ! hub
      !...................................
   CALL Linearize_Point_to_Point( y_ED%HubPtMotion, u_AD%HubMotion, MeshMapData%ED_P_2_AD_P_H, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//':u_AD%HubMotion' )      
      if (errStat>=AbortErrLev) return

      
   ! *** AD translational displacement: from ED translational displacement (MeshMapData%ED_P_2_AD_P_H%dM%mi) and orientation (MeshMapData%ED_P_2_AD_P_H%dM%fx_p)                              
   do i=1,size(MeshMapData%ED_P_2_AD_P_H%dM%mi,2)
      do j=1,size(MeshMapData%ED_P_2_AD_P_H%dM%mi,1)
         dUdy( AD_Start_tmp + j - 1, ED_Start_tmp + i - 1 ) = - MeshMapData%ED_P_2_AD_P_H%dM%mi(j,i)
      end do
   end do
   
   ED_Start_tmp = ED_Start_tmp + y_ED%HubPtMotion%NNodes * 3 ! move past the ED translation disp field to orientation field            
   do i=1,size(MeshMapData%ED_P_2_AD_P_H%dM%fx_p,2)
      do j=1,size(MeshMapData%ED_P_2_AD_P_H%dM%fx_p,1)
         dUdy( AD_Start_tmp + j - 1, ED_Start_tmp + i - 1 ) = - MeshMapData%ED_P_2_AD_P_H%dM%fx_p(j,i)
      end do
   end do
         
   ! *** AD orientation: from ED orientation
   AD_Start_tmp = AD_Start_tmp + u_AD%HubMotion%NNodes * 3 ! move past the AD translation disp field to orientation field         
   do i=1,size(MeshMapData%ED_P_2_AD_P_H%dM%mi,2)
      do j=1,size(MeshMapData%ED_P_2_AD_P_H%dM%mi,1)
         dUdy( AD_Start_tmp + j - 1, ED_Start_tmp + i - 1 ) = - MeshMapData%ED_P_2_AD_P_H%dM%mi(j,i)
      end do
   end do       
      
      
   ! *** AD rotational velocity: from ED rotational velocity
   AD_Start_tmp = AD_Start_tmp + u_AD%HubMotion%NNodes * 3 ! move past the AD orientation field to rotational velocity field          
   ED_Start_tmp = ED_Start_tmp + y_ED%HubPtMotion%NNodes * 3 ! move past the ED orientation field to rotational velocity field            
   do i=1,size(MeshMapData%ED_P_2_AD_P_H%dM%mi,2)
      do j=1,size(MeshMapData%ED_P_2_AD_P_H%dM%mi,1)
         dUdy( AD_Start_tmp + j - 1, ED_Start_tmp + i - 1 ) = - MeshMapData%ED_P_2_AD_P_H%dM%mi(j,i)
      end do
   end do       
            
   AD_Start_tmp = AD_Start_tmp + u_AD%HubMotion%NNodes * 3 ! 1 [of 3] fields ([already did translation disp, orientation,] rotation vel) with 3 components
   ED_Start_tmp = ED_Start_tmp + y_ED%HubPtMotion%NNodes * 3 ! 1 [of 3] fields ([already did translation disp, orientation,] rotational velocity) with 3 components            

   
      !...................................
      ! blade root   
      !...................................
   DO k=1,size(y_ED%BladeRootMotion)
      CALL Linearize_Point_to_Point( y_ED%BladeRootMotion(k), u_AD%BladeRootMotion(k), MeshMapData%ED_P_2_AD_P_R(k), ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//':u_AD%BladeRootMotion('//trim(num2lstr(k))//')' )      
         if (errStat>=AbortErrLev) return
               
      ! *** AD orientation: from ED orientation
      ED_Start_tmp = ED_Start_tmp + y_ED%BladeRootMotion(k)%NNodes * 3 ! move past ED translation disp field to orientation field
      do i=1,size(MeshMapData%ED_P_2_AD_P_R(k)%dM%mi,2)
         do j=1,size(MeshMapData%ED_P_2_AD_P_R(k)%dM%mi,1)
            dUdy( AD_Start_tmp + j - 1, ED_Start_tmp + i - 1 ) = - MeshMapData%ED_P_2_AD_P_R(k)%dM%mi(j,i)
         end do
      end do
                  
      AD_Start_tmp = AD_Start_tmp + u_AD%BladeRootMotion(k)%NNodes * 3 ! 1 field (orientation) with 3 components               
      ED_Start_tmp = ED_Start_tmp + y_ED%BladeRootMotion(k)%NNodes * 15 ! 5 [of 6] fields ([already took care of translation disp,] orientation, translation vel, rotation vel, translation acc, rotation acc) with 3 components
      
   END DO
      
   
      !...................................
      ! blades
      !...................................
   IF (p_FAST%CompElast == Module_ED ) THEN
            
      ED_Start_tmp = p_FAST%LinStartIndx(Module_ED, LIN_OUTPUT_COL)
      
      DO k=1,size(y_ED%BladeLn2Mesh)
         CALL Linearize_Line2_to_Line2( y_ED%BladeLn2Mesh(k), u_AD%BladeMotion(k), MeshMapData%BDED_L_2_AD_L_B(k), ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//':u_AD%BladeMotion('//trim(num2lstr(k))//')' )   
            if (errStat>=AbortErrLev) return
            
         ! *** AD translational displacement: from ED translational displacement (MeshMapData%BDED_L_2_AD_L_B(k)%dM%mi) and orientation (MeshMapData%BDED_L_2_AD_L_B(k)%dM%fx_p)                              
         do i=1,size(MeshMapData%BDED_L_2_AD_L_B(k)%dM%mi,2)
            do j=1,size(MeshMapData%BDED_L_2_AD_L_B(k)%dM%mi,1)
               dUdy( AD_Start_tmp + j - 1, ED_Start_tmp + i - 1 ) = - MeshMapData%BDED_L_2_AD_L_B(k)%dM%mi(j,i)
            end do
         end do
         
         ED_Start_tmp = ED_Start_tmp + y_ED%BladeLn2Mesh(k)%NNodes * 3 ! move past the ED translation disp field to orientation field            
         do i=1,size(MeshMapData%BDED_L_2_AD_L_B(k)%dM%fx_p,2)
            do j=1,size(MeshMapData%BDED_L_2_AD_L_B(k)%dM%fx_p,1)
               dUdy( AD_Start_tmp + j - 1, ED_Start_tmp + i - 1 ) = - MeshMapData%BDED_L_2_AD_L_B(k)%dM%fx_p(j,i)
            end do
         end do
         
         ! *** AD orientation: from ED orientation
         AD_Start_tmp = AD_Start_tmp + u_AD%BladeMotion(k)%NNodes * 3 ! move past the AD translation disp field to orientation field         
         do i=1,size(MeshMapData%BDED_L_2_AD_L_B(k)%dM%mi,2)
            do j=1,size(MeshMapData%BDED_L_2_AD_L_B(k)%dM%mi,1)
               dUdy( AD_Start_tmp + j - 1, ED_Start_tmp + i - 1 ) = - MeshMapData%BDED_L_2_AD_L_B(k)%dM%mi(j,i)
            end do
         end do 
         
         ! *** AD translational velocity: from ED translation displacement, translation velocity, and rotational velocity         
         AD_Start_tmp = AD_Start_tmp + u_AD%BladeMotion(k)%NNodes * 3 ! move past the AD orientation field to translation velocity field              
         ED_Start_tmp = ED_Start_tmp - y_ED%BladeLn2Mesh(k)%NNodes * 3 ! move past the ED translation disp field back to beginning of ED blades         
         do i=1,size(MeshMapData%BDED_L_2_AD_L_B(k)%dM%tv_uS,2)
            do j=1,size(MeshMapData%BDED_L_2_AD_L_B(k)%dM%tv_uS,1)
               dUdy( AD_Start_tmp + j - 1, ED_Start_tmp + i - 1 ) = - MeshMapData%BDED_L_2_AD_L_B(k)%dM%tv_uS(j,i)
            end do
         end do                              

         ED_Start_tmp = ED_Start_tmp + y_ED%BladeLn2Mesh(k)%NNodes * 6 ! move past the ED translation disp and orientation fields to translation vel field
         do i=1,size(MeshMapData%BDED_L_2_AD_L_B(k)%dM%mi,2)
            do j=1,size(MeshMapData%BDED_L_2_AD_L_B(k)%dM%mi,1)
               dUdy( AD_Start_tmp + j - 1, ED_Start_tmp + i - 1 ) = - MeshMapData%BDED_L_2_AD_L_B(k)%dM%mi(j,i)
            end do
         end do 
         
         ED_Start_tmp = ED_Start_tmp + y_ED%BladeLn2Mesh(k)%NNodes * 3 ! move past the ED translational vel field to rotational vel field
         do i=1,size(MeshMapData%BDED_L_2_AD_L_B(k)%dM%fx_p,2)
            do j=1,size(MeshMapData%BDED_L_2_AD_L_B(k)%dM%fx_p,1)
               dUdy( AD_Start_tmp + j - 1, ED_Start_tmp + i - 1 ) = - MeshMapData%BDED_L_2_AD_L_B(k)%dM%fx_p(j,i)
            end do
         end do
         
         
         AD_Start_tmp = AD_Start_tmp + u_AD%BladeMotion(k)%NNodes * 3 ! 1 [of 3] fields ([already took care of translation disp, orientation], translation vel) with 3 components
         ED_Start_tmp = ED_Start_tmp + y_ED%BladeLn2Mesh(k)%NNodes * 9 ! 3 [of 6] fields ([already took care of translation disp, orientation, translation vel], rotation vel, translation acc, rotation acc) with 3 components
            
      END DO
      
   !ELSEIF (p_FAST%CompElast == Module_BD ) THEN
   !   
   !      ! get them from BeamDyn
   !   DO k=1,size(u_AD%BladeMotion)
   !      CALL Linearize_Line2_to_Line2( BD%y(k)%BldMotion, u_AD%BladeMotion(k), MeshMapData%BDED_L_2_AD_L_B(k), ErrStat2, ErrMsg2 )
   !         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//':u_AD%BladeMotion('//trim(num2lstr(k))//')' )   
   !   END DO
      
            
   END IF
   
   
END SUBROUTINE Linear_AD_InputSolve_NoIfW_dy
!----------------------------------------------------------------------------------------------------------------------------------

!> This routine allocates the state matrices for the glue code and concatenates the module-level state matrices into
!! the first step of computing the full system state matrices. This routine returns
!! \f$ A = A^{ED} \f$, \f$ B = \begin{bmatrix} 0 & 0 & B^{ED} & 0 \end{bmatrix} \f$,
!! \f$ C = \begin{bmatrix} 0 \\ 0 \\ C^{ED} \\ 0 \end{bmatrix} \f$, and 
!! \f$ D = \begin{bmatrix} D^{IfW} & 0 & 0 & 0 \\ 0 &  D^{SrvD} & 0 & 0 \\ 0 & 0 &  D^{ED} & 0 \\ 0 & 0 & 0 &  D^{AD}\end{bmatrix}\f$.
SUBROUTINE Glue_FormDiag( p_FAST, y_FAST, ErrStat, ErrMsg )

   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST              !< Output variables for the glue code
           
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None

   
   
      ! local variables
   INTEGER(IntKi)                          :: ThisModule          ! Module ID # 
   
   INTEGER(IntKi)                          :: i                   ! module loop counter
   INTEGER(IntKi)                          :: r                   ! row loop counter
   INTEGER(IntKi)                          :: c                   ! column loop counter
   INTEGER(IntKi)                          :: r_start             ! row in glue matrix where module block matrix starts
   INTEGER(IntKi)                          :: c_start             ! column in glue matrix where module block matrix starts
   
   INTEGER(IntKi)                          :: ErrStat2            ! local error status
   CHARACTER(1024)                         :: ErrMsg2             ! local error message
   CHARACTER(*),             PARAMETER     :: RoutineName = 'Glue_FormDiag' 
   
   ErrStat = ErrID_None
   ErrMsg = ""
   
   
   
   !.....................................
   ! Allocate the state matrices if necessary:
   !.....................................

   if (.not. allocated(y_FAST%Lin%Glue%A)) then ! assume none of them are allocated
      ! A: rows = x; columns = x
      call AllocAry(y_FAST%Lin%Glue%A, p_FAST%SizeLin(NumModules+1,LIN_ContSTATE_COL), &
                                       p_FAST%SizeLin(NumModules+1,LIN_ContSTATE_COL), 'A', ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         
      !B: rows = x; columns = u
      call AllocAry(y_FAST%Lin%Glue%B, p_FAST%SizeLin(NumModules+1,LIN_ContSTATE_COL), &
                                       p_FAST%SizeLin(NumModules+1,LIN_INPUT_COL), 'B', ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         
      !C: rows = y; columns = x
      call AllocAry(y_FAST%Lin%Glue%C, p_FAST%SizeLin(NumModules+1,LIN_OUTPUT_COL), &
                                       p_FAST%SizeLin(NumModules+1,LIN_ContSTATE_COL), 'C', ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         
      !D: rows = y; columns = u
      call AllocAry(y_FAST%Lin%Glue%D, p_FAST%SizeLin(NumModules+1,LIN_OUTPUT_COL), &
                                       p_FAST%SizeLin(NumModules+1,LIN_INPUT_COL), 'D', ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         
      if ( p_FAST%SizeLin(NumModules+1, LIN_ExtINPUT_COL) > 0 ) then
         !B_ext: rows = x; columns = u extended
         call AllocAry(y_FAST%Lin%Glue%B_ext, p_FAST%SizeLin(NumModules+1,LIN_ContSTATE_COL), &
                                              p_FAST%SizeLin(NumModules+1,LIN_ExtINPUT_COL), 'B_ext', ErrStat2, ErrMsg2)
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         !D_ext: rows = y; columns = u extended
         call AllocAry(y_FAST%Lin%Glue%D_ext, p_FAST%SizeLin(NumModules+1,LIN_OUTPUT_COL), &
                                          p_FAST%SizeLin(NumModules+1,LIN_ExtINPUT_COL), 'D_ext', ErrStat2, ErrMsg2)
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      end if      
         
      if (ErrStat>=AbortErrLev) return
   end if
   
   
   ! The equations of the matrices returned from this routine are really just a general form with the null matrices removed:      
   
   ! A
   y_FAST%Lin%Glue%A = 0.0_ReKi
   r_start = 1
   c_start = 1
   do i = 1,p_FAST%Lin_NumMods
      ThisModule = p_FAST%Lin_ModOrder( i )
      
      if (allocated( y_FAST%Lin%Modules(ThisModule)%A) ) then
         do c=1,size( y_FAST%Lin%Modules(ThisModule)%A, 2)
            do r=1,size( y_FAST%Lin%Modules(ThisModule)%A, 1)
               y_FAST%Lin%Glue%A(r_start + r - 1, c_start + c - 1) = y_FAST%Lin%Modules(ThisModule)%A(r,c)
            end do
         end do
      end if
      
      r_start = r_start + p_FAST%SizeLin(ThisModule,LIN_ContSTATE_COL)
      c_start = c_start + p_FAST%SizeLin(ThisModule,LIN_ContSTATE_COL)      
   end do
   
    
   ! B
   y_FAST%Lin%Glue%B = 0.0_ReKi
   r_start = 1
   c_start = 1
   do i = 1,p_FAST%Lin_NumMods
      ThisModule = p_FAST%Lin_ModOrder( i )
      
      if (allocated( y_FAST%Lin%Modules(ThisModule)%B) ) then
         do c=1,size( y_FAST%Lin%Modules(ThisModule)%B, 2)
            do r=1,size( y_FAST%Lin%Modules(ThisModule)%B, 1)
               y_FAST%Lin%Glue%B(r_start + r - 1, c_start + c - 1) = y_FAST%Lin%Modules(ThisModule)%B(r,c)
            end do
         end do
      end if
      
      r_start = r_start + p_FAST%SizeLin(ThisModule,LIN_ContSTATE_COL)
      c_start = c_start + p_FAST%SizeLin(ThisModule,LIN_INPUT_COL)      
   end do
   
   ! C
   y_FAST%Lin%Glue%C = 0.0_ReKi
   r_start = 1
   c_start = 1
   do i = 1,p_FAST%Lin_NumMods
      ThisModule = p_FAST%Lin_ModOrder( i )
      
      if (allocated( y_FAST%Lin%Modules(ThisModule)%C) ) then
         do c=1,size( y_FAST%Lin%Modules(ThisModule)%C, 2)
            do r=1,size( y_FAST%Lin%Modules(ThisModule)%C, 1)
               y_FAST%Lin%Glue%C(r_start + r - 1, c_start + c - 1) = y_FAST%Lin%Modules(ThisModule)%C(r,c)
            end do
         end do
      end if
      
      r_start = r_start + p_FAST%SizeLin(ThisModule,LIN_OUTPUT_COL)
      c_start = c_start + p_FAST%SizeLin(ThisModule,LIN_ContSTATE_COL)      
   end do   
   
   ! D
   y_FAST%Lin%Glue%D = 0.0_ReKi
   r_start = 1
   c_start = 1
   do i = 1,p_FAST%Lin_NumMods
      ThisModule = p_FAST%Lin_ModOrder( i )
      
      if (allocated( y_FAST%Lin%Modules(ThisModule)%D) ) then
         do c=1,size( y_FAST%Lin%Modules(ThisModule)%D, 2)
            do r=1,size( y_FAST%Lin%Modules(ThisModule)%D, 1)
               y_FAST%Lin%Glue%D(r_start + r - 1, c_start + c - 1) = y_FAST%Lin%Modules(ThisModule)%D(r,c)
            end do
         end do
      end if
      
      r_start = r_start + p_FAST%SizeLin(ThisModule,LIN_OUTPUT_COL)
      c_start = c_start + p_FAST%SizeLin(ThisModule,LIN_INPUT_COL)      
   end do      
   
   
END SUBROUTINE Glue_FormDiag      
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine forms the full-system state matrices for linearization: A, B, C, and D.
!! Note that it uses LAPACK_GEMM instead of MATMUL for matrix multiplications because of stack-space issues (these
!! matrices get large quickly).
SUBROUTINE Glue_StateMatrices( p_FAST, y_FAST, NumBl, AD, IfW, dUdu, dUdy, ErrStat, ErrMsg )

   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST              !< Output variables for the glue code
   INTEGER(IntKi),           INTENT(IN   ) :: NumBl               !< Number of blades (for index into ED input array)
   TYPE(AeroDyn_Data),       INTENT(IN   ) :: AD                  !< AeroDyn data
   TYPE(InflowWind_Data),    INTENT(IN   ) :: IfW                 !< InflowWind data
   REAL(ReKi),               INTENT(INOUT) :: dUdu(:,:)           !< glue-code Jacobian: \f$ \frac{\partial U}{\partial u} \f$; on exit will hold G^{-1}*dUdu
   REAL(ReKi),               INTENT(INOUT) :: dUdy(:,:)           !< glue-code Jacobian: \f$ \frac{\partial U}{\partial y} \f$; on exit will hold G^{-1}*dUdy
           
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None

   
   
      ! local variables
   REAL(ReKi), ALLOCATABLE                 :: G(:,:), tmp(:,:) ! variables for glue-code linearization
   INTEGER(IntKi), ALLOCATABLE             :: ipiv(:)
            
   INTEGER(IntKi)                          :: ErrStat2            ! local error status
   CHARACTER(1024)                         :: ErrMsg2             ! local error message
   CHARACTER(*),             PARAMETER     :: RoutineName = 'Glue_StateMatrices' 
   
   ErrStat = ErrID_None
   ErrMsg = ""
   
   
   
   !.....................................   
   ! allocate the glue-code state matrices; after this call they will contain the state matrices from the 
   ! modules (without glue-code influence) on their diagonals
   !.....................................   
   call Glue_FormDiag( p_FAST, y_FAST, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (ErrStat >=AbortErrLev) then
         call cleanup()
         return
      end if
   
      
   !if (p_FAST%LinInputs == LIN_NONE .or. p_FAST%LinOutputs == LIN_NONE) then
   !       the glue-code input-output solve doesn't affect the rest of the equations, so we'll just return early
   !   call cleanup()
   !   return
   !end if
      
   
   !..................................... 
   ! solve for state matrices:
   !..................................... 
   
   ! *** get G matrix ****
   !----------------------
   if (.not. allocated(G)) then
      call AllocAry(G, size(dUdu,1), size(dUdu,2), 'G', ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
      call AllocAry( ipiv, size(dUdu,1), 'ipiv', ErrStat2, ErrMsg2 ) ! size(G,1)
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            
      if (ErrStat >=AbortErrLev) then
         call cleanup()
         return
      end if            
   end if      
   
   !G = dUdu + matmul( dUdy, y_FAST%Lin%Glue%D )            
   G = dUdu
   call LAPACK_GEMM( 'N', 'N', 1.0_ReKi, dUdy, y_FAST%Lin%Glue%D, 1.0_ReKi, G, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
   ! now we need to form G^-1 * dUdy and G^-1 * dUdu   
      ! factor G for the two solves:
   CALL LAPACK_getrf( M=size(G,1), N=size(G,2), A=G, IPIV=ipiv, ErrStat=ErrStat2, ErrMsg=ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (ErrStat >= AbortErrLev) then
         call cleanup() 
         return
      end if
      
    ! after the this solve, dUdy holds G^-1 * dUdy:
   CALL LAPACK_getrs( trans='N', N=size(G,2), A=G, IPIV=ipiv, B=dUdy, ErrStat=ErrStat2, ErrMsg=ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)   
    ! after the this solve, dUdu holds G^-1 * dUdu:
   CALL LAPACK_getrs( trans='N', N=size(G,2), A=G, IPIV=ipiv, B=dUdu, ErrStat=ErrStat2, ErrMsg=ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)   
         
   deallocate(G)    ! we're finished with the solves, so let's get rid of them
   deallocate(ipiv) ! we're finished with the solves, so let's get rid of them
      
                    
   ! *** get tmp matrix  for A and C calculations ****
   !----------------------         
   call AllocAry(tmp, p_FAST%SizeLin(NumModules+1,LIN_INPUT_COL), p_FAST%SizeLin(NumModules+1,LIN_ContSTATE_COL), 'G^-1*dUdy*C', ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (errStat>=AbortErrLev) then
         call cleanup()
         return
      end if
      
   !tmp = G^-1 * dUdy * diag(C)
   call LAPACK_GEMM( 'N', 'N', 1.0_ReKi, dUdy, y_FAST%Lin%Glue%C, 0.0_ReKi, tmp, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
                        
   !  A
   !----------------------         
   !> \f{equation}{ A = A^{ED} - \begin{bmatrix} 0 & 0 & B^{ED} & 0 \end{bmatrix} \,
   !! \begin{bmatrix} G \end{bmatrix}^{-1} \, \frac{\partial U}{\partial y} \, \begin{bmatrix} 0 \\ 0 \\ C^{ED} \\ 0 \end{bmatrix}
   !! }\f
   !y_FAST%Lin%Glue%A = y_FAST%Lin%Glue%A - matmul( y_FAST%Lin%Glue%B, tmp )  
   call LAPACK_GEMM( 'N', 'N', -1.0_ReKi, y_FAST%Lin%Glue%B, tmp, 1.0_ReKi, y_FAST%Lin%Glue%A, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

   !  C
   !----------------------         
   !> \f{equation}{ C = \begin{bmatrix} 0 \\ 0 \\ C^{ED} \\ 0 \end{bmatrix} - 
   !! \begin{bmatrix} D^{IfW} & 0 & 0 & 0 \\ 0 &  D^{SrvD} & 0 & 0 \\ 0 & 0 &  D^{ED} & 0 \\ 0 & 0 & 0 &  D^{AD}\end{bmatrix} \,
   !! \begin{bmatrix} G \end{bmatrix}^{-1} \, \frac{\partial U}{\partial y} \, \begin{bmatrix} 0 \\ 0 \\ C^{ED} \\ 0 \end{bmatrix}
   !! }\f
   !y_FAST%Lin%Glue%C = y_FAST%Lin%Glue%C - matmul( y_FAST%Lin%Glue%D, tmp ) 
   call LAPACK_GEMM( 'N', 'N', -1.0_ReKi, y_FAST%Lin%Glue%D, tmp, 1.0_ReKi, y_FAST%Lin%Glue%C, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
   deallocate(tmp)
   
         
   !  B
   !----------------------         
   !> \f{equation}{ B = \begin{bmatrix} 0 & 0 & B^{ED} & 0 \end{bmatrix} \,
   !! \begin{bmatrix} G \end{bmatrix}^{-1} \, \frac{\partial U}{\partial u}
   !! }\f
   call AllocAry(tmp,size(y_FAST%Lin%Glue%B,1),size(y_FAST%Lin%Glue%B,2),'tmp',ErrStat2,ErrMsg2)
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (errStat>=AbortErrLev) then
         call cleanup()
         return
      end if
   tmp = y_FAST%Lin%Glue%B   
         
   !y_FAST%Lin%Glue%B = matmul( y_FAST%Lin%Glue%B, dUdu ) 
   call LAPACK_GEMM( 'N', 'N', 1.0_ReKi, tmp, dUdu, 0.0_ReKi, y_FAST%Lin%Glue%B, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   deallocate(tmp)      
      
   if (allocated(y_FAST%Lin%Glue%B_ext)) then
      call ED_SetExtInputs(p_FAST, y_FAST%Lin%Glue%B, p_FAST%LinStartIndx(MODULE_ED, LIN_INPUT_COL), &
                    p_FAST%LinStartIndx(MODULE_ED, LIN_ExtINPUT_COL), NumBl, y_FAST%Lin%Glue%B_ext)
      if ( p_FAST%SizeLin(Module_AD, LIN_ExtINPUT_COL) > 0 ) then         
         call AD_SetExtInputs(p_FAST, y_FAST, y_FAST%Lin%Glue%B, p_FAST%LinStartIndx(MODULE_AD, LIN_INPUT_COL), &
                              p_FAST%LinStartIndx(MODULE_AD, LIN_ExtINPUT_COL), AD%Input(1), IfW%y, y_FAST%Lin%Glue%B_ext)                  
      end if      
   end if
   
    
   !  D
   !----------------------         
   !> \f{equation}{ D = \begin{bmatrix} D^{IfW} & 0 & 0 & 0 \\ 0 &  D^{SrvD} & 0 & 0 \\ 0 & 0 &  D^{ED} & 0 \\ 0 & 0 & 0 &  D^{AD}\end{bmatrix} \,
   !! \begin{bmatrix} G \end{bmatrix}^{-1} \, \frac{\partial U}{\partial u}
   !! }\f
   call AllocAry(tmp,size(y_FAST%Lin%Glue%D,1),size(y_FAST%Lin%Glue%D,2),'tmp',ErrStat2,ErrMsg2)
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (errStat>=AbortErrLev) then
         call cleanup()
         return
      end if
   tmp = y_FAST%Lin%Glue%D
         
   !y_FAST%Lin%Glue%D = matmul( y_FAST%Lin%Glue%D, dUdu )
   call LAPACK_GEMM( 'N', 'N', 1.0_ReKi, tmp, dUdu, 0.0_ReKi, y_FAST%Lin%Glue%D, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)   
   deallocate(tmp)    
   
   if (allocated(y_FAST%Lin%Glue%D_ext)) then
      call ED_SetExtInputs(p_FAST, y_FAST%Lin%Glue%D, p_FAST%LinStartIndx(MODULE_ED, LIN_INPUT_COL), &
                              p_FAST%LinStartIndx(MODULE_ED, LIN_ExtINPUT_COL), NumBl, y_FAST%Lin%Glue%D_ext)
      if ( p_FAST%SizeLin(Module_AD, LIN_ExtINPUT_COL) > 0 ) then         
         call AD_SetExtInputs(p_FAST, y_FAST, y_FAST%Lin%Glue%D, p_FAST%LinStartIndx(MODULE_AD, LIN_INPUT_COL), &
                              p_FAST%LinStartIndx(MODULE_AD, LIN_ExtINPUT_COL), AD%Input(1), IfW%y, y_FAST%Lin%Glue%D_ext)                  
      end if      
   end if

   
   call cleanup()
      
contains
   subroutine cleanup()
      if (allocated(ipiv)) deallocate(ipiv)     
      if (allocated(G)) deallocate(G)
      if (allocated(tmp)) deallocate(tmp)
   end subroutine cleanup   
END SUBROUTINE Glue_StateMatrices      
!----------------------------------------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------------------------------------
END MODULE FAST_Linear
