!**********************************************************************************************************************************
! FAST_Solver.f90, FAST_Subs.f90, FAST_Lin.f90, and FAST_Mods.f90 make up the FAST glue code in the FAST Modularization Framework.
! FAST_Prog.f90, FAST_Library.f90, FAST_Prog.c are different drivers for this code.
!..................................................................................................................................
! LICENSING
! Copyright (C) 2013-2016  National Renewable Energy Laboratory
! Copyright (C) 2018 Envision Energy USA, LTD
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

   INTEGER(IntKi)                          :: i, j, k             ! loop/temp variables
   INTEGER(IntKi)                          :: ThisModule          ! Module ID # 
   INTEGER(IntKi)                          :: NumInstances        ! Number of instances of each module
   
   INTEGER(IntKi)                          :: i_u                 ! loop/temp variables
   INTEGER(IntKi)                          :: i_y, i_x            ! loop/temp variables

   INTEGER(IntKi)                          :: NextStart(3)        ! allocated to be size(LinStartIndx)=size(SizeLin); helps compute the next starting index for the module components
   INTEGER(IntKi)                          :: ErrStat2            ! local error status
   CHARACTER(1024)                         :: ErrMsg2             ! local error message
   CHARACTER(*),             PARAMETER     :: RoutineName = 'Init_Lin' 
   CHARACTER(200)                          :: ModAbrev
   
   
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
      
      ! BeamDyn is next, if activated:
   if (p_FAST%CompElast == Module_BD) then
      p_FAST%Lin_NumMods = p_FAST%Lin_NumMods + 1
      p_FAST%Lin_ModOrder( p_FAST%Lin_NumMods ) = Module_BD
   end if
   
      ! AeroDyn is next, if activated:
   if ( p_FAST%CompAero  == Module_AD ) then 
      p_FAST%Lin_NumMods = p_FAST%Lin_NumMods + 1
      p_FAST%Lin_ModOrder( p_FAST%Lin_NumMods ) = Module_AD
   end if
   
   
   !.....................
   ! determine total number of inputs/outputs/contStates:
   !.....................
   y_FAST%Lin%Glue%SizeLin = 0
   
   do i = 1,p_FAST%Lin_NumMods
      ThisModule = p_FAST%Lin_ModOrder( i )
      
      do k=1,size(y_FAST%Lin%Modules(ThisModule)%Instance)
         y_FAST%Lin%Modules(ThisModule)%Instance(k)%SizeLin = 0
         if (allocated(y_FAST%Lin%Modules(ThisModule)%Instance(k)%Names_u)) y_FAST%Lin%Modules(ThisModule)%Instance(k)%SizeLin(LIN_INPUT_COL)     = size(y_FAST%Lin%Modules(ThisModule)%Instance(k)%Names_u)
         if (allocated(y_FAST%Lin%Modules(ThisModule)%Instance(k)%Names_y)) y_FAST%Lin%Modules(ThisModule)%Instance(k)%SizeLin(LIN_OUTPUT_COL)    = size(y_FAST%Lin%Modules(ThisModule)%Instance(k)%Names_y)
         if (allocated(y_FAST%Lin%Modules(ThisModule)%Instance(k)%Names_x)) y_FAST%Lin%Modules(ThisModule)%Instance(k)%SizeLin(LIN_ContSTATE_COL) = size(y_FAST%Lin%Modules(ThisModule)%Instance(k)%Names_x)
         
         y_FAST%Lin%Glue%SizeLin = y_FAST%Lin%Glue%SizeLin + y_FAST%Lin%Modules(ThisModule)%Instance(k)%SizeLin ! total number of inputs, outputs, and continuous states
      end do
   end do
   
   !.....................
   ! compute the starting index in the combined (full) matrices:
   !.....................
   NextStart = 1 ! whole array
   do i = 1,p_FAST%Lin_NumMods
      ThisModule = p_FAST%Lin_ModOrder( i )
      
      do k=1,size(y_FAST%Lin%Modules(ThisModule)%Instance)
         y_FAST%Lin%Modules(ThisModule)%Instance(k)%LinStartIndx = NextStart
         NextStart = NextStart + y_FAST%Lin%Modules(ThisModule)%Instance(k)%SizeLin
      end do
   end do
   
   
      ! ...................................
      ! determine which of the module inputs/outputs are written to file
      ! ...................................
   !NumBl = size(u_ED%BlPitchCom)   
   call Init_Lin_InputOutput(p_FAST, y_FAST, NumBl, ErrStat2, ErrMsg2)
      call SetErrStat(errStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   
      ! ...................................
      ! get names of inputs, outputs, and continuous states
      ! ...................................
   call AllocAry( y_FAST%Lin%Glue%names_u, y_FAST%Lin%Glue%SizeLin(LIN_INPUT_COL), 'names_u', ErrStat2, ErrMsg2)
      call SetErrStat(errStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call AllocAry( y_FAST%Lin%Glue%names_y, y_FAST%Lin%Glue%SizeLin(LIN_OUTPUT_COL), 'names_y', ErrStat2, ErrMsg2)
      call SetErrStat(errStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)      
   call AllocAry( y_FAST%Lin%Glue%names_x, y_FAST%Lin%Glue%SizeLin(LIN_ContSTATE_COL), 'names_x', ErrStat2, ErrMsg2)
      call SetErrStat(errStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)       
   call AllocAry( y_FAST%Lin%Glue%Use_u, y_FAST%Lin%Glue%SizeLin(LIN_INPUT_COL), 'use_u', ErrStat2, ErrMsg2)
      call SetErrStat(errStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)      
   call AllocAry( y_FAST%Lin%Glue%Use_y, y_FAST%Lin%Glue%SizeLin(LIN_OUTPUT_COL), 'use_y', ErrStat2, ErrMsg2)
      call SetErrStat(errStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)      
   call AllocAry( y_FAST%Lin%Glue%RotFrame_u, y_FAST%Lin%Glue%SizeLin(LIN_INPUT_COL), 'RotFrame_u', ErrStat2, ErrMsg2)
      call SetErrStat(errStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call AllocAry( y_FAST%Lin%Glue%RotFrame_y, y_FAST%Lin%Glue%SizeLin(LIN_OUTPUT_COL), 'RotFrame_y', ErrStat2, ErrMsg2)
      call SetErrStat(errStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)      
   call AllocAry( y_FAST%Lin%Glue%RotFrame_x, y_FAST%Lin%Glue%SizeLin(LIN_ContSTATE_COL), 'RotFrame_x', ErrStat2, ErrMsg2)
      call SetErrStat(errStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)       
   call AllocAry( y_FAST%Lin%Glue%IsLoad_u, y_FAST%Lin%Glue%SizeLin(LIN_INPUT_COL), 'IsLoad_u', ErrStat2, ErrMsg2)
      call SetErrStat(errStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
   if (ErrStat >= AbortErrLev) return
               
   
   i_u = 1
   i_y = 1      
   i_x = 1      
   do i = 1,p_FAST%Lin_NumMods
      ThisModule = p_FAST%Lin_ModOrder( i )
      ModAbrev = y_FAST%Module_Abrev(ThisModule)
      NumInstances = size(y_FAST%Lin%Modules(ThisModule)%Instance)

         ! inputs
      do k=1,NumInstances
         if (NumInstances > 1 .or. trim(y_FAST%Module_Abrev(ThisModule)) == "BD") then
            ModAbrev = TRIM(y_FAST%Module_Abrev(ThisModule))//'_'//trim(num2lstr(k))
         end if
      
         do j=1,y_FAST%Lin%Modules(ThisModule)%Instance(k)%SizeLin(LIN_INPUT_COL)
            y_FAST%Lin%Glue%names_u(i_u) = TRIM(ModAbrev)//' '//y_FAST%Lin%Modules(ThisModule)%Instance(k)%Names_u(j)
            y_FAST%Lin%Glue%use_u(  i_u) = y_FAST%Lin%Modules(ThisModule)%Instance(k)%use_u(j)
            y_FAST%Lin%Glue%IsLoad_u(i_u) = y_FAST%Lin%Modules(ThisModule)%Instance(k)%IsLoad_u(j) 
         
            if (allocated(y_FAST%Lin%Modules(ThisModule)%Instance(k)%RotFrame_u)) then
               y_FAST%Lin%Glue%RotFrame_u(i_u) = y_FAST%Lin%Modules(ThisModule)%Instance(k)%RotFrame_u(j) 
            else 
               y_FAST%Lin%Glue%RotFrame_u(i_u) = .false.
            end if
            i_u = i_u + 1;
         end do

      end do
      
         ! outputs
      do k=1,size(y_FAST%Lin%Modules(ThisModule)%Instance)
         if (NumInstances > 1 .or. trim(y_FAST%Module_Abrev(ThisModule)) == "BD") then
            ModAbrev = TRIM(y_FAST%Module_Abrev(ThisModule))//'_'//trim(num2lstr(k))
         end if
            
         do j=1,y_FAST%Lin%Modules(ThisModule)%Instance(k)%SizeLin(LIN_OUTPUT_COL)
            y_FAST%Lin%Glue%names_y(i_y) = TRIM(ModAbrev)//' '//y_FAST%Lin%Modules(ThisModule)%Instance(k)%Names_y(j)
            y_FAST%Lin%Glue%use_y(  i_y) = y_FAST%Lin%Modules(ThisModule)%Instance(k)%use_y(j)
            if (allocated(y_FAST%Lin%Modules(ThisModule)%Instance(k)%RotFrame_y)) then
               y_FAST%Lin%Glue%RotFrame_y(i_y) = y_FAST%Lin%Modules(ThisModule)%Instance(k)%RotFrame_y(j)
            else 
               y_FAST%Lin%Glue%RotFrame_y(i_y) = .false.
            end if
            i_y = i_y + 1;
         end do
      end do
 
         ! continuous states
      do k=1,size(y_FAST%Lin%Modules(ThisModule)%Instance)
         if (NumInstances > 1 .or. trim(y_FAST%Module_Abrev(ThisModule)) == "BD") then
            ModAbrev = TRIM(y_FAST%Module_Abrev(ThisModule))//'_'//trim(num2lstr(k))
         end if

         do j=1,y_FAST%Lin%Modules(ThisModule)%Instance(k)%SizeLin(LIN_ContSTATE_COL)
            y_FAST%Lin%Glue%names_x( i_x) = TRIM(ModAbrev)//' '//y_FAST%Lin%Modules(ThisModule)%Instance(k)%Names_x( j)
            if (allocated(y_FAST%Lin%Modules(ThisModule)%Instance(k)%RotFrame_x)) then
               y_FAST%Lin%Glue%RotFrame_x(i_x) = y_FAST%Lin%Modules(ThisModule)%Instance(k)%RotFrame_x(j) 
            else 
               y_FAST%Lin%Glue%RotFrame_x(i_x) = .false.
            end if
            i_x = i_x + 1;
         end do
      end do
      
   end do ! each module

END SUBROUTINE Init_Lin
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine that initializes the names and rotating frame portion of IfW.
SUBROUTINE Init_Lin_IfW( p_FAST, y_FAST, u_AD )

   TYPE(FAST_ParameterType),       INTENT(IN   )   :: p_FAST      !< FAST parameter data 
   TYPE(FAST_OutputFileType),      INTENT(INOUT)   :: y_FAST      !< Output variables for the glue code
   TYPE(AD_InputType),             INTENT(IN)      :: u_AD        !< The input meshes (already calculated) from AeroDyn   
   
   INTEGER(IntKi)                          :: i, j, k             ! loop counters
   INTEGER(IntKi)                          :: i2                  ! loop counters
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
            position = index(y_FAST%Lin%Modules(Module_IfW)%Instance(1)%Names_u(i), ',') - 1
            y_FAST%Lin%Modules(Module_IfW)%Instance(1)%Names_u(i) = y_FAST%Lin%Modules(Module_IfW)%Instance(1)%Names_u(i)(1:position)//trim(NodeDesc)//&
                                                                    y_FAST%Lin%Modules(Module_IfW)%Instance(1)%Names_u(i)(position+1:)
         end do    
         do i=1,3
            position = index(y_FAST%Lin%Modules(Module_IfW)%Instance(1)%Names_y(i), ',') - 1
            y_FAST%Lin%Modules(Module_IfW)%Instance(1)%Names_y(i) = y_FAST%Lin%Modules(Module_IfW)%Instance(1)%Names_y(i)(1:position)//trim(NodeDesc)//&
                                                                    y_FAST%Lin%Modules(Module_IfW)%Instance(1)%Names_y(i)(position+1:)
         end do    
      end if
                  
      IF (p_FAST%CompAero == MODULE_AD) THEN 
                           
         DO K = 1,SIZE(u_AD%BladeMotion)
            DO J = 1,u_AD%BladeMotion(k)%Nnodes
               Node = Node + 1 ! InflowWind node
               NodeDesc = ' (blade '//trim(num2lstr(k))//', node '//trim(num2lstr(j))//')'
               
               do i=1,3 !XYZ components of this node
                  i2 = (Node-1)*3 + i
                                    
                  position = index(y_FAST%Lin%Modules(Module_IfW)%Instance(1)%Names_u(i2), ',') - 1
                  y_FAST%Lin%Modules(Module_IfW)%Instance(1)%Names_u(i2) = y_FAST%Lin%Modules(Module_IfW)%Instance(1)%Names_u(i2)(1:position)//trim(NodeDesc)//&
                                                                           y_FAST%Lin%Modules(Module_IfW)%Instance(1)%Names_u(i2)(position+1:)
                                                       
                  position = index(y_FAST%Lin%Modules(Module_IfW)%Instance(1)%Names_y(i2), ',') - 1
                  y_FAST%Lin%Modules(Module_IfW)%Instance(1)%Names_y(i2) = y_FAST%Lin%Modules(Module_IfW)%Instance(1)%Names_y(i2)(1:position)//trim(NodeDesc)//&
                                                                           y_FAST%Lin%Modules(Module_IfW)%Instance(1)%Names_y(i2)(position+1:)
                  
                  ! IfW has inputs and outputs in the global frame
                  !y_FAST%Lin%Modules(Module_IfW)%Instance(1)%RotFrame_u(i2) = .true.
                  !y_FAST%Lin%Modules(Module_IfW)%Instance(1)%RotFrame_y(i2) = .true.
                  
               end do            
            END DO !J = 1,p%BldNodes ! Loop through the blade nodes / elements
         END DO !K = 1,p%NumBl     
         
            ! tower:
         DO J=1,u_AD%TowerMotion%nnodes
            Node = Node + 1  
            NodeDesc = ' (Tower node '//trim(num2lstr(j))//')'

            do i=1,3 !XYZ components of this node
               i2 = (Node-1)*3 + i
                                    
               position = index(y_FAST%Lin%Modules(Module_IfW)%Instance(1)%Names_u(i2), ',') - 1
               y_FAST%Lin%Modules(Module_IfW)%Instance(1)%Names_u(i2) = y_FAST%Lin%Modules(Module_IfW)%Instance(1)%Names_u(i2)(1:position)//trim(NodeDesc)//&
                                                                        y_FAST%Lin%Modules(Module_IfW)%Instance(1)%Names_u(i2)(position+1:)
                                     
               position = index(y_FAST%Lin%Modules(Module_IfW)%Instance(1)%Names_y(i2), ',') - 1
               y_FAST%Lin%Modules(Module_IfW)%Instance(1)%Names_y(i2) = y_FAST%Lin%Modules(Module_IfW)%Instance(1)%Names_y(i2)(1:position)//trim(NodeDesc)//&
                                                                        y_FAST%Lin%Modules(Module_IfW)%Instance(1)%Names_y(i2)(position+1:)
            end do            
         END DO              
         
      END IF     
   
END SUBROUTINE Init_Lin_IfW
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine that initializes some use_u and use_y, which determine which, if any, inputs and outputs are output in the linearization file.
SUBROUTINE Init_Lin_InputOutput(p_FAST, y_FAST, NumBl, ErrStat, ErrMsg)

   TYPE(FAST_ParameterType), INTENT(INOUT) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST              !< Output variables for the glue code
   INTEGER(IntKi),           INTENT(IN   ) :: NumBl               !< Number of blades (for index into ED input array)
   
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None

   INTEGER(IntKi)                          :: i, j, col           ! loop/temp variables
   INTEGER(IntKi)                          :: k                   ! loop/temp variables
   INTEGER(IntKi)                          :: ThisModule          ! Module ID # 
   
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
      
      do k=1,size(y_FAST%Lin%Modules(ThisModule)%Instance)
         call AllocAry ( y_FAST%Lin%Modules(ThisModule)%Instance(k)%Use_u, size(y_FAST%Lin%Modules(ThisModule)%Instance(k)%Names_u), TRIM(y_FAST%Module_Abrev(ThisModule))//'_Use_u', ErrStat2, ErrMsg2)
            call SetErrStat(errStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         call AllocAry ( y_FAST%Lin%Modules(ThisModule)%Instance(k)%Use_y, size(y_FAST%Lin%Modules(ThisModule)%Instance(k)%Names_y), TRIM(y_FAST%Module_Abrev(ThisModule))//'_Use_y', ErrStat2, ErrMsg2)
            call SetErrStat(errStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      end do

   end do
   if (ErrStat >= AbortErrLev) return
   
   
      ! ...................................
      ! set true/false flags for inputs:
      ! ...................................
   
   if (p_FAST%LinInputs == LIN_NONE) then
      do i = 1,p_FAST%Lin_NumMods
         ThisModule = p_FAST%Lin_ModOrder( i )
         do k=1,size(y_FAST%Lin%Modules(ThisModule)%Instance)
            y_FAST%Lin%Modules(ThisModule)%Instance(k)%use_u = .false.
         end do
      end do
   elseif(p_FAST%LinInputs == LIN_STANDARD) then
      do i = 1,p_FAST%Lin_NumMods
         ThisModule = p_FAST%Lin_ModOrder( i )
         do k=1,size(y_FAST%Lin%Modules(ThisModule)%Instance)
            y_FAST%Lin%Modules(ThisModule)%Instance(k)%use_u = .false.
         end do
      end do
      
      ! ED standard inputs: BlPitchCom, YawMom, GenTrq, extended input (collective pitch)
      do j=1,NumBl+3
         y_FAST%Lin%Modules(MODULE_ED)%Instance(1)%use_u(y_FAST%Lin%Modules(MODULE_ED)%Instance(1)%SizeLin(LIN_INPUT_COL)+1-j) = .true.
      end do
      
      ! IfW standard inputs: HWindSpeed, PLexp, PropagationDir
      if (p_FAST%CompInflow == MODULE_IfW) then
         do j = 1,3
            y_FAST%Lin%Modules(MODULE_IfW)%Instance(1)%use_u(y_FAST%Lin%Modules(MODULE_IfW)%Instance(1)%SizeLin(LIN_INPUT_COL)+1-j) = .true.
         end do
      end if
                  
   elseif(p_FAST%LinInputs == LIN_ALL) then
      do i = 1,p_FAST%Lin_NumMods
         ThisModule = p_FAST%Lin_ModOrder( i )
         do k=1,size(y_FAST%Lin%Modules(ThisModule)%Instance)
            y_FAST%Lin%Modules(ThisModule)%Instance(k)%use_u = .true.
         end do
      end do      
   end if
            
        
      ! ...................................
      ! set true/false flags for outputs:
      ! ...................................
   
   if (p_FAST%LinOutputs == LIN_NONE) then
      do i = 1,p_FAST%Lin_NumMods
         ThisModule = p_FAST%Lin_ModOrder( i )
         do k=1,size(y_FAST%Lin%Modules(ThisModule)%Instance)
            y_FAST%Lin%Modules(ThisModule)%Instance(k)%use_y = .false.
         end do
      end do      
   elseif(p_FAST%LinOutputs == LIN_STANDARD) then

      ! WriteOutput values are the last entries of the modules      
      do i = 1,p_FAST%Lin_NumMods         
         ThisModule = p_FAST%Lin_ModOrder( i )
         
         do k=1,size(y_FAST%Lin%Modules(ThisModule)%Instance)
            col = y_FAST%Lin%Modules(ThisModule)%Instance(k)%SizeLin(LIN_OUTPUT_COL) - y_FAST%Lin%Modules(ThisModule)%Instance(k)%NumOutputs !first column where WriteOutput occurs
            do j=1,col
               y_FAST%Lin%Modules(ThisModule)%Instance(k)%use_y(j) = .false.
            end do
            do j=col+1,y_FAST%Lin%Modules(ThisModule)%Instance(k)%SizeLin(LIN_OUTPUT_COL)
               y_FAST%Lin%Modules(ThisModule)%Instance(k)%use_y(j) = .true.
            end do
         end do
      end do      
      
   elseif(p_FAST%LinOutputs == LIN_ALL) then
      do i = 1,p_FAST%Lin_NumMods
         ThisModule = p_FAST%Lin_ModOrder( i )
         do k=1,size(y_FAST%Lin%Modules(ThisModule)%Instance)
            y_FAST%Lin%Modules(ThisModule)%Instance(k)%use_y = .true.
         end do
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
   
   REAL(R8Ki), ALLOCATABLE                 :: dYdz(:,:), dZdz(:,:), dZdu(:,:)
   REAL(R8Ki), ALLOCATABLE                 :: dUdu(:,:), dUdy(:,:) ! variables for glue-code linearization
   INTEGER(IntKi), ALLOCATABLE             :: ipiv(:)
   integer(intki)                          :: NumBl
   integer(intki)                          :: k
   CHARACTER(1024)                         :: LinRootName
   CHARACTER(1024)                         :: OutFileName
   
   
   
   ErrStat = ErrID_None
   ErrMsg = ""
   Un = -1

   
   LinRootName = TRIM(p_FAST%OutFileRoot)//'.'//trim(num2lstr(m_FAST%NextLinTimeIndx))
   
   NumBl = size(ED%Input(1)%BlPitchCom) 
   y_FAST%Lin%RotSpeed = ED%Output(1)%RotSpeed
   y_FAST%Lin%Azimuth  = ED%Output(1)%LSSTipPxa
   !.....................
   ! ElastoDyn
   !.....................
      ! get the jacobians
   call ED_JacobianPInput( t_global, ED%Input(1), ED%p, ED%x(STATE_CURR), ED%xd(STATE_CURR), ED%z(STATE_CURR), ED%OtherSt(STATE_CURR), &
                              ED%Output(1), ED%m, ErrStat2, ErrMsg2, dYdu=y_FAST%Lin%Modules(Module_ED)%Instance(1)%D, dXdu=y_FAST%Lin%Modules(Module_ED)%Instance(1)%B )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
   call ED_JacobianPContState( t_global, ED%Input(1), ED%p, ED%x(STATE_CURR), ED%xd(STATE_CURR), ED%z(STATE_CURR), ED%OtherSt(STATE_CURR), &
                                  ED%Output(1), ED%m, ErrStat2, ErrMsg2, dYdx=y_FAST%Lin%Modules(Module_ED)%Instance(1)%C, dXdx=y_FAST%Lin%Modules(Module_ED)%Instance(1)%A )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   
      ! get the operating point
   call ED_GetOP( t_global, ED%Input(1), ED%p, ED%x(STATE_CURR), ED%xd(STATE_CURR), ED%z(STATE_CURR), ED%OtherSt(STATE_CURR), &
                     ED%Output(1), ED%m, ErrStat2, ErrMsg2, u_op=y_FAST%Lin%Modules(Module_ED)%Instance(1)%op_u, y_op=y_FAST%Lin%Modules(Module_ED)%Instance(1)%op_y, &
                    x_op=y_FAST%Lin%Modules(Module_ED)%Instance(1)%op_x, dx_op=y_FAST%Lin%Modules(Module_ED)%Instance(1)%op_dx )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (ErrStat >=AbortErrLev) then
         call cleanup()
         return
      end if
      
      
      ! write the module matrices:
   if (p_FAST%LinOutMod) then
            
      OutFileName = trim(LinRootName)//'.'//TRIM(y_FAST%Module_Abrev(Module_ED))      
      call WrLinFile_txt_Head(t_global, p_FAST, y_FAST, y_FAST%Lin%Modules(Module_ED)%Instance(1), OutFileName, Un, ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat >=AbortErrLev) then
            call cleanup()
            return
         end if
         
      if (p_FAST%LinOutJac) then
         ! Jacobians
         !dXdx:
         call WrPartialMatrix( y_FAST%Lin%Modules(Module_ED)%Instance(1)%A, Un, p_FAST%OutFmt, 'dXdx' )
         
         !dXdu:
         call WrPartialMatrix( y_FAST%Lin%Modules(Module_ED)%Instance(1)%B, Un, p_FAST%OutFmt, 'dXdu', UseCol=y_FAST%Lin%Modules(Module_ED)%Instance(1)%use_u )
         
         ! dYdx:
         call WrPartialMatrix( y_FAST%Lin%Modules(Module_ED)%Instance(1)%C, Un, p_FAST%OutFmt, 'dYdx', UseRow=y_FAST%Lin%Modules(Module_ED)%Instance(1)%use_y )
         
         !dYdu:
         call WrPartialMatrix( y_FAST%Lin%Modules(Module_ED)%Instance(1)%D, Un, p_FAST%OutFmt, 'dYdu', UseRow=y_FAST%Lin%Modules(Module_ED)%Instance(1)%use_y, &
                                                                                                       UseCol=y_FAST%Lin%Modules(Module_ED)%Instance(1)%use_u )
         
      end if
      
         ! finish writing the file
      call WrLinFile_txt_End(Un, p_FAST, y_FAST%Lin%Modules(Module_ED)%Instance(1) )
               
   end if
   
   !.....................
   ! BeamDyn
   !.....................
   if ( p_FAST%CompElast  == Module_BD ) then
      do k=1,p_FAST%nBeams

         ! get the jacobians
         call BD_JacobianPInput( t_global, BD%Input(1,k), BD%p(k), BD%x(k,STATE_CURR), BD%xd(k,STATE_CURR), BD%z(k,STATE_CURR), BD%OtherSt(k,STATE_CURR), &
                                    BD%y(k), BD%m(k), ErrStat2, ErrMsg2, dYdu=y_FAST%Lin%Modules(Module_BD)%Instance(k)%D, &
                                    dXdu=y_FAST%Lin%Modules(Module_BD)%Instance(k)%B, &
                                    StateRel_x   =y_FAST%Lin%Modules(Module_BD)%Instance(k)%StateRel_x, &
                                    StateRel_xdot=y_FAST%Lin%Modules(Module_BD)%Instance(k)%StateRel_xdot )
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
         call BD_JacobianPContState( t_global, BD%Input(1,k), BD%p(k), BD%x(k,STATE_CURR), BD%xd(k,STATE_CURR), BD%z(k,STATE_CURR), BD%OtherSt(k,STATE_CURR), &
                                    BD%y(k), BD%m(k), ErrStat2, ErrMsg2, dYdx=y_FAST%Lin%Modules(Module_BD)%Instance(k)%C, dXdx=y_FAST%Lin%Modules(Module_BD)%Instance(k)%A, &
                                    StateRotation=y_FAST%Lin%Modules(Module_BD)%Instance(k)%StateRotation)
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   
            ! get the operating point
         call BD_GetOP( t_global, BD%Input(1,k), BD%p(k), BD%x(k,STATE_CURR), BD%xd(k,STATE_CURR), BD%z(k,STATE_CURR), BD%OtherSt(k,STATE_CURR), &
                        BD%y(k), BD%m(k), ErrStat2, ErrMsg2, u_op=y_FAST%Lin%Modules(Module_BD)%Instance(k)%op_u,  y_op=y_FAST%Lin%Modules(Module_BD)%Instance(k)%op_y, &
                                                             x_op=y_FAST%Lin%Modules(Module_BD)%Instance(k)%op_x, dx_op=y_FAST%Lin%Modules(Module_BD)%Instance(k)%op_dx )
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            if (ErrStat >=AbortErrLev) then
               call cleanup()
               return
            end if
      
      
            ! write the module matrices:
         if (p_FAST%LinOutMod) then
            
            OutFileName = trim(LinRootName)//'.'//TRIM(y_FAST%Module_Abrev(Module_BD))//TRIM(num2lstr(k))
            call WrLinFile_txt_Head(t_global, p_FAST, y_FAST, y_FAST%Lin%Modules(Module_BD)%Instance(k), OutFileName, Un, ErrStat2, ErrMsg2 )
               call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
               if (ErrStat >=AbortErrLev) then
                  call cleanup()
                  return
               end if
         
            if (p_FAST%LinOutJac) then
               ! Jacobians
               !dXdx:
               call WrPartialMatrix( y_FAST%Lin%Modules(Module_BD)%Instance(k)%A, Un, p_FAST%OutFmt, 'dXdx' )
         
               !dXdu:
               call WrPartialMatrix( y_FAST%Lin%Modules(Module_BD)%Instance(k)%B, Un, p_FAST%OutFmt, 'dXdu', UseCol=y_FAST%Lin%Modules(Module_BD)%Instance(k)%use_u )
         
               !dYdx:
               call WrPartialMatrix( y_FAST%Lin%Modules(Module_BD)%Instance(k)%C, Un, p_FAST%OutFmt, 'dYdx', UseRow=y_FAST%Lin%Modules(Module_BD)%Instance(k)%use_y )
         
               !dYdu:
               call WrPartialMatrix( y_FAST%Lin%Modules(Module_BD)%Instance(k)%D, Un, p_FAST%OutFmt, 'dYdu', UseRow=y_FAST%Lin%Modules(Module_BD)%Instance(k)%use_y, &
                                                                                                             UseCol=y_FAST%Lin%Modules(Module_BD)%Instance(k)%use_u )
            end if
      
               ! finish writing the file
            call WrLinFile_txt_End(Un, p_FAST, y_FAST%Lin%Modules(Module_BD)%Instance(k) )
         end if

      end do
   end if
   !.....................
   ! InflowWind
   !.....................      
   if ( p_FAST%CompInflow  == Module_IfW ) then 
      
         ! get the jacobians
      call InflowWind_JacobianPInput( t_global, IfW%Input(1), IfW%p, IfW%x(STATE_CURR), IfW%xd(STATE_CURR), IfW%z(STATE_CURR), &
                                   IfW%OtherSt(STATE_CURR), IfW%y, IfW%m, ErrStat2, ErrMsg2, dYdu=y_FAST%Lin%Modules(Module_IfW)%Instance(1)%D )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
      ! get the operating point
      call InflowWind_GetOP( t_global, IfW%Input(1), IfW%p, IfW%x(STATE_CURR), IfW%xd(STATE_CURR), IfW%z(STATE_CURR), &
                             IfW%OtherSt(STATE_CURR), IfW%y, IfW%m, ErrStat2, ErrMsg2, u_op=y_FAST%Lin%Modules(Module_IfW)%Instance(1)%op_u, &
                       y_op=y_FAST%Lin%Modules(Module_IfW)%Instance(1)%op_y )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat >=AbortErrLev) then
            call cleanup()
            return
         end if
                      
      
         ! write the module matrices:
      if (p_FAST%LinOutMod) then
               
         OutFileName = trim(LinRootName)//'.'//TRIM(y_FAST%Module_Abrev(Module_IfW))      
         call WrLinFile_txt_Head(t_global, p_FAST, y_FAST, y_FAST%Lin%Modules(Module_IfW)%Instance(1), OutFileName, Un, ErrStat2, ErrMsg2 )
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            if (ErrStat >=AbortErrLev) then
               call cleanup()
               return
            end if
         
         if (p_FAST%LinOutJac) then
            ! Jacobians
            !dYdu:
            call WrPartialMatrix( y_FAST%Lin%Modules(Module_IfW)%Instance(1)%D, Un, p_FAST%OutFmt, 'dYdu', &
               UseRow=y_FAST%Lin%Modules(Module_IfW)%Instance(1)%use_y, UseCol=y_FAST%Lin%Modules(Module_IfW)%Instance(1)%use_u )
         end if
      
            ! finish writing the file
         call WrLinFile_txt_End(Un, p_FAST, y_FAST%Lin%Modules(Module_IfW)%Instance(1) )
               
      end if      
            
   end if
   
   !.....................
   ! ServoDyn
   !.....................   
   if ( p_FAST%CompServo  == Module_SrvD ) then 
         ! get the jacobians
      call SrvD_JacobianPInput( t_global, SrvD%Input(1), SrvD%p, SrvD%x(STATE_CURR), SrvD%xd(STATE_CURR), SrvD%z(STATE_CURR), &
                                   SrvD%OtherSt(STATE_CURR), SrvD%y, SrvD%m, ErrStat2, ErrMsg2, dYdu=y_FAST%Lin%Modules(Module_SrvD)%Instance(1)%D )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
      ! get the operating point
      call SrvD_GetOP( t_global, SrvD%Input(1), SrvD%p, SrvD%x(STATE_CURR), SrvD%xd(STATE_CURR), SrvD%z(STATE_CURR), &
                       SrvD%OtherSt(STATE_CURR), SrvD%y, SrvD%m, ErrStat2, ErrMsg2, u_op=y_FAST%Lin%Modules(Module_SrvD)%Instance(1)%op_u, &
                       y_op=y_FAST%Lin%Modules(Module_SrvD)%Instance(1)%op_y )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat >=AbortErrLev) then
            call cleanup()
            return
         end if
            
         ! write the module matrices:
      if (p_FAST%LinOutMod) then
      
         OutFileName = trim(LinRootName)//'.'//TRIM(y_FAST%Module_Abrev(Module_SrvD))      
         call WrLinFile_txt_Head(t_global, p_FAST, y_FAST, y_FAST%Lin%Modules(Module_SrvD)%Instance(1), OutFileName, Un, ErrStat2, ErrMsg2 )
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            if (ErrStat >=AbortErrLev) then
               call cleanup()
               return
            end if
         
            ! Jacobians
         if (p_FAST%LinOutJac) then
            !dYdu:
            call WrPartialMatrix( y_FAST%Lin%Modules(Module_SrvD)%Instance(1)%D, Un, p_FAST%OutFmt, 'dYdu', &
               UseRow=y_FAST%Lin%Modules(Module_SrvD)%Instance(1)%use_y, UseCol=y_FAST%Lin%Modules(Module_SrvD)%Instance(1)%use_u )
         end if
      
            ! finish writing the file
         call WrLinFile_txt_End(Un, p_FAST, y_FAST%Lin%Modules(Module_SrvD)%Instance(1) )
               
      end if      
   end if

   !.....................
   ! AeroDyn
   !.....................
   if ( p_FAST%CompAero  == Module_AD ) then 
         ! get the jacobians
      call AD_JacobianPInput( t_global, AD%Input(1), AD%p, AD%x(STATE_CURR), AD%xd(STATE_CURR), AD%z(STATE_CURR), &
                                   AD%OtherSt(STATE_CURR), AD%y, AD%m, ErrStat2, ErrMsg2, dYdu=y_FAST%Lin%Modules(Module_AD)%Instance(1)%D, dZdu=dZdu )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
      call AD_JacobianPConstrState( t_global, AD%Input(1), AD%p, AD%x(STATE_CURR), AD%xd(STATE_CURR), AD%z(STATE_CURR), &
                                   AD%OtherSt(STATE_CURR), AD%y, AD%m, ErrStat2, ErrMsg2, dYdz=dYdz, dZdz=dZdz )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      ! get the operating point
      call AD_GetOP( t_global, AD%Input(1), AD%p, AD%x(STATE_CURR), AD%xd(STATE_CURR), AD%z(STATE_CURR), &
                       AD%OtherSt(STATE_CURR), AD%y, AD%m, ErrStat2, ErrMsg2, u_op=y_FAST%Lin%Modules(Module_AD)%Instance(1)%op_u, &
                       y_op=y_FAST%Lin%Modules(Module_AD)%Instance(1)%op_y, z_op=y_FAST%Lin%Modules(Module_AD)%Instance(1)%op_z )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat >=AbortErrLev) then
            call cleanup()
            return
         end if               
      
         ! write the module matrices:
      if (p_FAST%LinOutMod) then
      
         OutFileName = trim(LinRootName)//'.'//TRIM(y_FAST%Module_Abrev(Module_AD))      
         call WrLinFile_txt_Head(t_global, p_FAST, y_FAST, y_FAST%Lin%Modules(Module_AD)%Instance(1), OutFileName, Un, ErrStat2, ErrMsg2 )
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
            call WrPartialMatrix( dZdu, Un, p_FAST%OutFmt, 'dZdu', UseCol=y_FAST%Lin%Modules(Module_AD)%Instance(1)%use_u )
            
            ! dYdz:
            call WrPartialMatrix( dYdz, Un, p_FAST%OutFmt, 'dYdz', UseRow=y_FAST%Lin%Modules(Module_AD)%Instance(1)%use_y )
            
            !dYdu:
            call WrPartialMatrix( y_FAST%Lin%Modules(Module_AD)%Instance(1)%D, Un, p_FAST%OutFmt, 'dYdu', &
                  UseRow=y_FAST%Lin%Modules(Module_AD)%Instance(1)%use_y, UseCol=y_FAST%Lin%Modules(Module_AD)%Instance(1)%use_u )
            
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
      !y_FAST%Lin%Modules(Module_AD)%D = y_FAST%Lin%Modules(Module_AD)%D - matmul(dYdz, dZdu )
      call LAPACK_GEMM( 'N', 'N', -1.0_R8Ki, dYdz, dZdu, 1.0_R8Ki, y_FAST%Lin%Modules(Module_AD)%Instance(1)%D, ErrStat2, ErrMsg2 )
      
      
      if (p_FAST%LinOutMod) then
            ! finish writing the file
         call WrLinFile_txt_End(Un, p_FAST, y_FAST%Lin%Modules(Module_AD)%Instance(1) )
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
      call WrPartialMatrix( dUdu, Un, p_FAST%OutFmt, 'dUdu', UseRow=y_FAST%Lin%Glue%use_u, UseCol=y_FAST%Lin%Glue%use_u )
      call WrPartialMatrix( dUdy, Un, p_FAST%OutFmt, 'dUdy', UseRow=y_FAST%Lin%Glue%use_u, UseCol=y_FAST%Lin%Glue%use_y )
   end if
   
   
      ! calculate the glue-code state matrices
   call Glue_StateMatrices( p_FAST, y_FAST, dUdu, dUdy, ErrStat2, ErrMsg2 )
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
      
      if (allocated(dUdu)) deallocate(dUdu)
      if (allocated(dUdy)) deallocate(dUdy)
      
      if (Un > 0) close(Un)

   end subroutine cleanup
END SUBROUTINE FAST_Linearize_OP   
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
   INTEGER(IntKi)                          :: n(5)                ! sizes of arrays to print
   CHARACTER(*),             PARAMETER     :: TypeNames(5) = (/ 'continuous states', &
                                                                'discrete states  ', &
                                                                'constraint states', &
                                                                'inputs           ', &
                                                                'outputs          '  /)
   CHARACTER(50)                           :: Fmt
   CHARACTER(32)                           :: Desc
   
   integer, parameter :: Indx_x      = 1
   integer, parameter :: Indx_xd     = 2
   integer, parameter :: Indx_z      = 3
   integer, parameter :: Indx_u      = 4
   integer, parameter :: Indx_y      = 5
   
                  
   ErrStat = ErrID_None
   ErrMsg = ""
         
   n = 0;
   if (allocated(LinData%names_x )) n(Indx_x) = size(LinData%names_x )
   if (allocated(LinData%names_xd)) n(Indx_xd) = size(LinData%names_xd)
   if (allocated(LinData%names_z )) n(Indx_z) = size(LinData%names_z )
   !if (allocated(LinData%names_u )) n(Indx_u) = size(LinData%names_u )
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
      call WrLinFile_txt_Table(p_FAST, Un, "Row/Column", LinData%op_z, LinData%names_z, rotFrame=LinData%RotFrame_z  )      
   end if
         
   if (n(Indx_u) > 0) then
      WRITE(Un, '(A)') 'Order of inputs:'   
      call WrLinFile_txt_Table(p_FAST, Un, "Column  ", LinData%op_u, LinData%names_u, rotFrame=LinData%RotFrame_u, UseCol=LinData%use_u  )
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

   INTEGER(IntKi),           INTENT(INOUT) :: Un                  !< unit number
   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< parameters
   TYPE(FAST_LinType),       INTENT(IN   ) :: LinData             !< Linearization data for individual module or glue (coupled system)
   
      ! local variables
   
   WRITE (Un,'(/,A,/)' ) 'Linearized state matrices:'    !print a blank line
   
   ! A matrix
   if (allocated(LinData%A)) call WrPartialMatrix( LinData%A, Un, p_FAST%OutFmt, 'A' )
   ! B matrix   
   if (allocated(LinData%B)) call WrPartialMatrix( LinData%B, Un, p_FAST%OutFmt, 'B', UseCol=LinData%use_u )
   
   ! C matrix
   if (allocated(LinData%C)) call WrPartialMatrix( LinData%C, Un, p_FAST%OutFmt, 'C', UseRow=LinData%use_y )
   ! D matrix
   if (allocated(LinData%D)) call WrPartialMatrix( LinData%D, Un, p_FAST%OutFmt, 'D', UseRow=LinData%use_y, UseCol=LinData%use_u )

   ! StateRotation matrix
   if (allocated(LinData%StateRotation)) call WrPartialMatrix( LinData%StateRotation, Un, p_FAST%OutFmt, 'StateRotation' )

   ! RelState matrices
   if (allocated(LinData%StateRel_x))    call WrPartialMatrix( LinData%StateRel_x,    Un, p_FAST%OutFmt, 'State_Rel_x' )
   if (allocated(LinData%StateRel_xdot)) call WrPartialMatrix( LinData%StateRel_xdot, Un, p_FAST%OutFmt, 'State_Rel_xdot' )

   close(Un)
   Un = -1
   
END SUBROUTINE WrLinFile_txt_End   
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE WrLinFile_txt_Table(p_FAST, Un, RowCol, op, names, rotFrame, deriv, UseCol,start_indx)

   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< parameters
   INTEGER(IntKi),           INTENT(IN   ) :: Un                  !< unit number
   CHARACTER(*),             INTENT(IN   ) :: RowCol              !< Row/Column description
   REAL(ReKi),               INTENT(IN   ) :: op(:)               !< operating point values (possibly different size that Desc because of orientations)
   CHARACTER(LinChanLen),    INTENT(IN   ) :: names(:)            !< Descriptions of the channels (names and units)
   logical, optional,        INTENT(IN   ) :: rotFrame(:)         !< determines if this parameter is in the rotating frame
   logical, optional,        intent(in   ) :: deriv               !< flag that tells us if we need to modify the channel names for derivatives (xdot)
   logical, optional,        intent(in   ) :: UseCol(:)           !< flags that tell us if we should use each column or skip it
   INTEGER(IntKi),optional,  INTENT(IN   ) :: start_indx          !< starting index (so extended inputs can be numbered starting after the # of inputs)
   
      ! local variables
   INTEGER(IntKi)                          :: TS                  ! tab stop column
   INTEGER(IntKi)                          :: i, i_print          ! loop counter
   INTEGER(IntKi)                          :: i_op                ! loop counter
   
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
   if (present(start_indx)) then
      i_print = start_indx + 1
   else      
      i_print = 1
   end if
   
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

   
   INTEGER(IntKi)                          :: i, j, k             ! loop/temp variables
   INTEGER(IntKi)                          :: ThisModule          ! Module ID # 
   INTEGER(IntKi)                          :: i_u                 ! loop/temp variables
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
      i_y = 0
      do i = 1,p_FAST%Lin_NumMods
         ThisModule = p_FAST%Lin_ModOrder( i )
         
         do k=1,size(y_FAST%Lin%Modules(ThisModule)%Instance)
            if (allocated(y_FAST%Lin%Modules(ThisModule)%Instance(k)%op_u)) then
               i_u = i_u + size(y_FAST%Lin%Modules(ThisModule)%Instance(k)%op_u)
            end if
                  
            if (allocated(y_FAST%Lin%Modules(ThisModule)%Instance(k)%op_y)) then
               i_y = i_y + size(y_FAST%Lin%Modules(ThisModule)%Instance(k)%op_y)
            end if
         end do
      end do      
      
      call AllocAry( y_FAST%Lin%Glue%op_u, i_u, 'op_u', ErrStat2, ErrMsg2)
         call SetErrStat(errStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      call AllocAry( y_FAST%Lin%Glue%op_y, i_y, 'op_y', ErrStat2, ErrMsg2)
         call SetErrStat(errStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      call AllocAry( y_FAST%Lin%Glue%op_x, y_FAST%Lin%Glue%SizeLin(LIN_ContSTATE_COL), 'op_x', ErrStat2, ErrMsg2)
         call SetErrStat(errStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      call AllocAry( y_FAST%Lin%Glue%op_dx, y_FAST%Lin%Glue%SizeLin(LIN_ContSTATE_COL), 'op_dx', ErrStat2, ErrMsg2)
         call SetErrStat(errStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (ErrStat >=AbortErrLev) return
   end if
   
   
   i_u = 1
   i_y = 1      
   i_x = 1
   do i = 1,p_FAST%Lin_NumMods
      ThisModule = p_FAST%Lin_ModOrder( i )

      do k=1,size(y_FAST%Lin%Modules(ThisModule)%Instance)
         if (allocated(y_FAST%Lin%Modules(ThisModule)%Instance(k)%op_u)) then
            do j=1,size(y_FAST%Lin%Modules(ThisModule)%Instance(k)%op_u)
               y_FAST%Lin%Glue%op_u(i_u) = y_FAST%Lin%Modules(ThisModule)%Instance(k)%op_u(j)
               i_u = i_u + 1;
            end do
         end if
               
         if (allocated(y_FAST%Lin%Modules(ThisModule)%Instance(k)%op_y)) then
            do j=1,size(y_FAST%Lin%Modules(ThisModule)%Instance(k)%op_y)
               y_FAST%Lin%Glue%op_y(i_y) = y_FAST%Lin%Modules(ThisModule)%Instance(k)%op_y(j)
               i_y = i_y + 1;
            end do
         end if

         if (allocated(y_FAST%Lin%Modules(ThisModule)%Instance(k)%op_x)) then
            do j=1,size(y_FAST%Lin%Modules(ThisModule)%Instance(k)%op_x)
               y_FAST%Lin%Glue%op_x(i_x) = y_FAST%Lin%Modules(ThisModule)%Instance(k)%op_x(j)
            
               y_FAST%Lin%Glue%op_dx(i_x) = y_FAST%Lin%Modules(ThisModule)%Instance(k)%op_dx(j)
               i_x = i_x + 1;
            end do
         end if
      end do
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
   REAL(R8Ki), ALLOCATABLE,  INTENT(INOUT) :: dUdu(:,:)           !< Partial derivatives of input-output equations (U(y,u)=0) with respect
                                                                  !!   to the inputs (u)
   REAL(R8Ki), ALLOCATABLE,  INTENT(INOUT) :: dUdy(:,:)           !< Partial derivatives of input-output equations (U(y,u)=0) with respect
                                                                  !!   to the outputs (y)
      
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None

   
   
      ! local variables
   INTEGER(IntKi)                          :: ThisModule          ! Module ID
   INTEGER(IntKi)                          :: i, j, k             ! loop counter
   INTEGER(IntKi)                          :: r_start, r_end      ! row start/end of glue matrix
   
   INTEGER(IntKi)                          :: ErrStat2            ! local error status
   CHARACTER(1024)                         :: ErrMsg2             ! local error message
   CHARACTER(*),             PARAMETER     :: RoutineName = 'Glue_Jacobians' 
   
   ErrStat = ErrID_None
   ErrMsg = ""
   
   
   ! Note: Where the same Linearize_*_to_*() routines for mesh mapping are used in both dUdu and dUdy, the dUdy routines assume dUdu 
   ! has already called the routine (and so avoids calling the routines a second time). This means the dUdu routines must be called first.
   
   !.....................................
   ! dUdu 
   !> \f$ \frac{\partial U_\Lambda}{\partial u} =  
   !!  \begin{bmatrix} \frac{\partial U_\Lambda^{IfW}}{\partial u^{IfW}} & \frac{\partial U_\Lambda^{IfW}}{\partial u^{SrvD}} & 
   !!                  \frac{\partial U_\Lambda^{IfW}}{\partial u^{ED}}  & \frac{\partial U_\Lambda^{IfW}}{\partial u^{BD}} & \frac{\partial U_\Lambda^{IfW}}{\partial u^{AD}} \\
   !!                  \frac{\partial U_\Lambda^{SrvD}}{\partial u^{IfW}} & \frac{\partial U_\Lambda^{SrvD}}{\partial u^{SrvD}} & 
   !!                  \frac{\partial U_\Lambda^{SrvD}}{\partial u^{ED}}  & \frac{\partial U_\Lambda^{SrvD}}{\partial u^{BD}} & \frac{\partial U_\Lambda^{SrvD}}{\partial u^{AD}} \\
   !!                  \frac{\partial U_\Lambda^{ED}}{\partial u^{IfW}} & \frac{\partial U_\Lambda^{ED}}{\partial u^{SrvD}} & 
   !!                  \frac{\partial U_\Lambda^{ED}}{\partial u^{ED}}  & \frac{\partial U_\Lambda^{ED}}{\partial u^{BD}} & \frac{\partial U_\Lambda^{ED}}{\partial u^{AD}} \\
   !!                  \frac{\partial U_\Lambda^{BD}}{\partial u^{IfW}} & \frac{\partial U_\Lambda^{BD}}{\partial u^{SrvD}} & 
   !!                  \frac{\partial U_\Lambda^{BD}}{\partial u^{ED}}  & \frac{\partial U_\Lambda^{BD}}{\partial u^{BD}} & \frac{\partial U_\Lambda^{BD}}{\partial u^{AD}} \\
   !!                  \frac{\partial U_\Lambda^{AD}}{\partial u^{IfW}} & \frac{\partial U_\Lambda^{AD}}{\partial u^{SrvD}} & 
   !!                  \frac{\partial U_\Lambda^{AD}}{\partial u^{ED}}  & \frac{\partial U_\Lambda^{AD}}{\partial u^{BD}} & \frac{\partial U_\Lambda^{AD}}{\partial u^{AD}} \\
   !!  \end{bmatrix} = 
   !!  \begin{bmatrix} I & 0 & 0 & 0                                               & \frac{\partial U_\Lambda^{IfW}}{\partial u^{AD}} \\
   !!                  0 & I & 0 & 0                                               & 0 \\
   !!                  0 & 0 & I & \frac{\partial U_\Lambda^{ED}}{\partial u^{BD}} & \frac{\partial U_\Lambda^{ED}}{\partial u^{AD}} \\
   !!                  0 & 0 & 0 & \frac{\partial U_\Lambda^{BD}}{\partial u^{BD}} & \frac{\partial U_\Lambda^{BD}}{\partial u^{AD}} \\
   !!                  0 & 0 & 0 & 0                                               & \frac{\partial U_\Lambda^{AD}}{\partial u^{AD}} \\
   !!  \end{bmatrix} \f$
   !.....................................

   if (.not. allocated(dUdu)) then
      call AllocAry(dUdu, y_FAST%Lin%Glue%SizeLin(LIN_INPUT_COL), y_FAST%Lin%Glue%SizeLin(LIN_INPUT_COL), 'dUdu', ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat>=AbortErrLev) return
   end if
   
   dUdu = 0.0_R8Ki      ! most of this matrix is zero, so we'll just initialize everything and set only the non-zero parts below
   
   
      !............
      !  \f$ \frac{\partial U_\Lambda^{IfW}}{\partial u^{IfW}} = I \f$ \n
      !  \f$ \frac{\partial U_\Lambda^{SrvD}}{\partial u^{SrvD}} = I \f$ \n
      !  \f$ \frac{\partial U_\Lambda^{ED}}{\partial u^{ED}} = I \f$ \n
      !  Note that we're also doing \f$ \frac{\partial U_\Lambda^{BD}}{\partial u^{BD}} = I \f$ and 
      !  \f$ \frac{\partial U_\Lambda^{AD}}{\partial u^{AD}} = I \f$ here; We will add values to the off=diagonal terms of those block matrices later.
      !............
   do j = 1,p_FAST%Lin_NumMods
      ThisModule = p_FAST%Lin_ModOrder(j)
      do k=1,size(y_FAST%Lin%Modules(ThisModule)%Instance)
         r_start =           y_FAST%Lin%Modules(ThisModule)%Instance(k)%LinStartIndx(LIN_INPUT_COL)
         r_end   = r_start + y_FAST%Lin%Modules(ThisModule)%Instance(k)%SizeLin(     LIN_INPUT_COL) - 1
         do i = r_start,r_end
            dUdu(i,i) = 1.0_R8Ki
         end do
      end do
   end do
   
   
      !............
      ! \f$ \frac{\partial U_\Lambda^{IfW}}{\partial u^{AD}} \end{bmatrix} = \f$   (dUdu block row 1=IfW)
      !............
   IF (p_FAST%CompInflow == MODULE_IfW .and. p_FAST%CompAero == MODULE_AD) THEN  
      call Linear_IfW_InputSolve_du_AD( p_FAST, y_FAST, AD%Input(1), dUdu )
   end if ! we're using the InflowWind module
   
      !............ 
      ! \f$ \frac{\partial U_\Lambda^{ED}}{\partial u^{AD}} \end{bmatrix} = \f$ and 
      ! \f$ \frac{\partial U_\Lambda^{ED}}{\partial u^{BD}} \end{bmatrix} = \f$ (dUdu block row 3=ED)
      !............   
   ! we need to do this for CompElast=ED and CompElast=BD
   call Linear_ED_InputSolve_du( p_FAST, y_FAST, ED%Input(1), ED%Output(1), AD%y, AD%Input(1), BD, MeshMapData, dUdu, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   
      !............
      ! \f$ \frac{\partial U_\Lambda^{BD}}{\partial u^{AD}} \end{bmatrix} = \f$ and 
      ! \f$ \frac{\partial U_\Lambda^{BD}}{\partial u^{BD}} \end{bmatrix} = \f$ (dUdu block row 4=BD)
      !............   
   IF (p_FAST%CompElast == Module_BD) THEN
      call Linear_BD_InputSolve_du( p_FAST, y_FAST, ED%Output(1), AD%y, AD%Input(1), BD, MeshMapData, dUdu, ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   END IF

      !............
      ! \f$ \frac{\partial U_\Lambda^{AD}}{\partial u^{AD}} \end{bmatrix} = \f$ (dUdu block row 5=AD)
      !............
   IF (p_FAST%CompAero == MODULE_AD) THEN 
      call Linear_AD_InputSolve_du( p_FAST, y_FAST, AD%Input(1), ED%Output(1), BD, MeshMapData, dUdu, ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   end if ! we're using the InflowWind module
   
   !.....................................
   ! dUdy
   !> \f$ \frac{\partial U_\Lambda}{\partial y} =  
   !!  \begin{bmatrix} \frac{\partial U_\Lambda^{IfW} }{\partial y^{IfW}} & \frac{\partial U_\Lambda^{IfW} }{\partial y^{SrvD}} & 
   !!                  \frac{\partial U_\Lambda^{IfW} }{\partial y^{ED} } & \frac{\partial U_\Lambda^{IfW} }{\partial y^{BD}  } & \frac{\partial U_\Lambda^{IfW} }{\partial y^{AD}} \\
   !!                  \frac{\partial U_\Lambda^{SrvD}}{\partial y^{IfW}} & \frac{\partial U_\Lambda^{SrvD}}{\partial y^{SrvD}} & 
   !!                  \frac{\partial U_\Lambda^{SrvD}}{\partial y^{ED} } & \frac{\partial U_\Lambda^{SrvD}}{\partial y^{BD}  } & \frac{\partial U_\Lambda^{SrvD}}{\partial y^{AD}} \\
   !!                  \frac{\partial U_\Lambda^{ED}  }{\partial y^{IfW}} & \frac{\partial U_\Lambda^{ED}  }{\partial y^{SrvD}} & 
   !!                  \frac{\partial U_\Lambda^{ED}  }{\partial y^{ED} } & \frac{\partial U_\Lambda^{ED}  }{\partial y^{BD}  } & \frac{\partial U_\Lambda^{ED}  }{\partial y^{AD}} \\
   !!                  \frac{\partial U_\Lambda^{BD}  }{\partial y^{IfW}} & \frac{\partial U_\Lambda^{BD}  }{\partial y^{SrvD}} & 
   !!                  \frac{\partial U_\Lambda^{BD}  }{\partial y^{ED} } & \frac{\partial U_\Lambda^{BD}  }{\partial y^{BD}  } & \frac{\partial U_\Lambda^{BD}  }{\partial y^{AD}} \\
   !!                  \frac{\partial U_\Lambda^{AD}  }{\partial y^{IfW}} & \frac{\partial U_\Lambda^{AD}  }{\partial y^{SrvD}} & 
   !!                  \frac{\partial U_\Lambda^{AD}  }{\partial y^{ED} } & \frac{\partial U_\Lambda^{AD}  }{\partial y^{BD}  } & \frac{\partial U_\Lambda^{AD}  }{\partial y^{AD}} \\
   !!  \end{bmatrix} =
   !!  \begin{bmatrix} 0                                                & 0                                                 & 0                                                 & 0                                               & 0 \\
   !!                  0                                                & 0                                                 & \frac{\partial U_\Lambda^{SrvD}}{\partial y^{ED}} & 0                                               & 0 \\
   !!                  0                                                & \frac{\partial U_\Lambda^{ED}}{\partial y^{SrvD}} & \frac{\partial U_\Lambda^{ED}}{\partial y^{ED}}   & \frac{\partial U_\Lambda^{ED}}{\partial y^{BD}} & \frac{\partial U_\Lambda^{ED}}{\partial y^{AD}} \\
   !!                  0                                                & 0                                                 & \frac{\partial U_\Lambda^{BD}}{\partial y^{ED}}   & \frac{\partial U_\Lambda^{BD}}{\partial y^{BD}} & \frac{\partial U_\Lambda^{BD}}{\partial y^{AD}} \\
   !!                  \frac{\partial U_\Lambda^{AD}}{\partial y^{IfW}} & 0                                                 & \frac{\partial U_\Lambda^{AD}}{\partial y^{ED}}   & \frac{\partial U_\Lambda^{AD}}{\partial y^{BD}} & 0 \\
   !!  \end{bmatrix} \f$
   !.....................................
   if (.not. allocated(dUdy)) then
      call AllocAry(dUdy, y_FAST%Lin%Glue%SizeLin(LIN_INPUT_COL), y_FAST%Lin%Glue%SizeLin(LIN_OUTPUT_COL), 'dUdy', ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat>=AbortErrLev) return
   end if
         
   dUdy = 0.0_R8Ki      ! most of this matrix is zero, so we'll just initialize everything and set only the non-zero parts below
   
      !............
      ! \f$ \frac{\partial U_\Lambda^{SrvD}}{\partial y^{ED}} \end{bmatrix} = \f$ (dUdy block row 2=SrvD)
      !............
   if (p_FAST%CompServo == MODULE_SrvD) then   ! need to do this regardless of CompElast
      call Linear_SrvD_InputSolve_dy( p_FAST, y_FAST, dUdy )
   end if
   
      !............
      ! \f$ \frac{\partial U_\Lambda^{ED}}{\partial y^{SrvD}} \end{bmatrix} = \f$
      ! \f$ \frac{\partial U_\Lambda^{ED}}{\partial y^{ED}} \end{bmatrix} = \f$
      ! \f$ \frac{\partial U_\Lambda^{ED}}{\partial y^{BD}} \end{bmatrix} = \f$
      ! \f$ \frac{\partial U_\Lambda^{ED}}{\partial y^{AD}} \end{bmatrix} = \f$ (dUdy block row 3=ED)
      !............
   call Linear_ED_InputSolve_dy( p_FAST, y_FAST, ED%Input(1), ED%Output(1), AD%y, AD%Input(1), BD, MeshMapData, dUdy, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   
      !............
      ! \f$ \frac{\partial U_\Lambda^{BD}}{\partial y^{ED}} \end{bmatrix} = \f$
      ! \f$ \frac{\partial U_\Lambda^{BD}}{\partial y^{BD}} \end{bmatrix} = \f$
      ! \f$ \frac{\partial U_\Lambda^{BD}}{\partial y^{AD}} \end{bmatrix} = \f$ (dUdy block row 4=BD)
      !............
   if (p_FAST%CompElast == MODULE_BD) then
      call Linear_BD_InputSolve_dy( p_FAST, y_FAST, ED%Input(1), ED%Output(1), AD%y, AD%Input(1), BD, MeshMapData, dUdy, ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   end if
   
      !............
      ! \f$ \frac{\partial U_\Lambda^{AD}}{\partial y^{IfW}} \end{bmatrix} = \f$
      ! \f$ \frac{\partial U_\Lambda^{AD}}{\partial y^{ED}} \end{bmatrix} = \f$
      ! \f$ \frac{\partial U_\Lambda^{AD}}{\partial y^{BD}} \end{bmatrix} = \f$ (dUdy block row 5=AD)
      !............
   if (p_FAST%CompAero == MODULE_AD) then   ! need to do this regardless of CompElast
   
      if (p_FAST%CompInflow == MODULE_IfW) then
         call Linear_AD_InputSolve_IfW_dy( p_FAST, y_FAST, AD%Input(1), dUdy )
      end if

      call Linear_AD_InputSolve_NoIfW_dy( p_FAST, y_FAST, AD%Input(1), ED%Output(1), BD, MeshMapData, dUdy, ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

   end if
      
   
   
END SUBROUTINE Glue_Jacobians
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine forms the dU^{IfW}/du^{AD} block of dUdu. (i.e., how do changes in the AD inputs affect IfW inputs?)
SUBROUTINE Linear_IfW_InputSolve_du_AD( p_FAST, y_FAST, u_AD, dUdu )

   TYPE(FAST_ParameterType),       INTENT(IN   )   :: p_FAST      !< FAST parameter data
   TYPE(FAST_OutputFileType),      INTENT(IN   )   :: y_FAST      !< FAST output data (for linearization)
   TYPE(AD_InputType),             INTENT(IN)      :: u_AD        !< The input meshes (already calculated) from AeroDyn
   REAL(R8Ki),                     INTENT(INOUT)   :: dUdu(:,:)   !< Jacobian matrix of which we are computing the dU^(IfW)/du^(AD) block
   
   
   INTEGER(IntKi)                          :: i, j, k             ! loop counters
   INTEGER(IntKi)                          :: i2, j2              ! loop counters
   INTEGER(IntKi)                          :: AD_Start_Bl         ! starting index of dUdu (column) where AD blade motion inputs are located
   INTEGER(IntKi)                          :: Node                ! InflowWind node number
   
            
         ! compare with IfW_InputSolve():
   
      Node = 0 !InflowWind node
      if (p_FAST%CompServo == MODULE_SrvD) Node = Node + 1
            
      IF (p_FAST%CompAero == MODULE_AD) THEN 
         
            ! blades:
         AD_Start_Bl = y_FAST%Lin%Modules(MODULE_AD)%Instance(1)%LinStartIndx(LIN_INPUT_COL) &
                     + u_AD%TowerMotion%NNodes * 9  & ! 3 fields (MASKID_TRANSLATIONDISP,MASKID_Orientation,MASKID_TRANSLATIONVel) with 3 components
                     + u_AD%HubMotion%NNodes   * 9    ! 3 fields (MASKID_TRANSLATIONDISP,MASKID_Orientation,MASKID_RotationVel) with 3 components
   
         do k = 1,size(u_AD%BladeRootMotion)         
            AD_Start_Bl = AD_Start_Bl + u_AD%BladeRootMotion(k)%NNodes * 3 ! 1 field (MASKID_Orientation) with 3 components
         end do
         ! next is u_AD%BladeMotion(k):
                  
         DO K = 1,SIZE(u_AD%BladeMotion)
            DO J = 1,u_AD%BladeMotion(k)%Nnodes
               Node = Node + 1 ! InflowWind node
               do i=1,3 !XYZ components of this node
                  i2 = y_FAST%Lin%Modules(MODULE_IfW)%Instance(1)%LinStartIndx(LIN_INPUT_COL) + (Node-1)*3 + i - 1
                  j2 = AD_Start_Bl + (j-1)*3 + i - 1
                  dUdu( i2, j2 ) = -1.0_R8Ki
               end do            
            END DO !J = 1,p%BldNodes ! Loop through the blade nodes / elements
                     
               ! get starting AD index of BladeMotion for next blade
            AD_Start_Bl = AD_Start_Bl + u_AD%BladeMotion(k)%Nnodes * 9  ! 3 fields (MASKID_TRANSLATIONDISP,MASKID_Orientation,MASKID_TRANSLATIONVel) with 3 components
         END DO !K = 1,p%NumBl     
         
            ! tower:
         DO J=1,u_AD%TowerMotion%nnodes
            Node = Node + 1   
            do i=1,3 !XYZ components of this node
               i2 = y_FAST%Lin%Modules(MODULE_IfW)%Instance(1)%LinStartIndx(LIN_INPUT_COL) + (Node-1)*3 + i - 1
               j2 = y_FAST%Lin%Modules(MODULE_AD )%Instance(1)%LinStartIndx(LIN_INPUT_COL) +    (j-1)*3 + i - 1
               dUdu( i2, j2 ) = -1.0_R8Ki
            end do            
         END DO              
         
      END IF     
   
END SUBROUTINE Linear_IfW_InputSolve_du_AD
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine forms the dU^{ED}/du^{BD} and dU^{ED}/du^{AD} blocks (ED row) of dUdu. (i.e., how do changes in the AD and BD inputs affect the ED inputs?)
SUBROUTINE Linear_ED_InputSolve_du( p_FAST, y_FAST, u_ED, y_ED, y_AD, u_AD, BD, MeshMapData, dUdu, ErrStat, ErrMsg )

   TYPE(FAST_ParameterType),       INTENT(IN   )  :: p_FAST         !< Glue-code simulation parameters
   TYPE(FAST_OutputFileType),      INTENT(IN   )  :: y_FAST         !< Glue-code output parameters (for linearization)
   TYPE(ED_InputType),             INTENT(INOUT)  :: u_ED           !< ED Inputs at t
   TYPE(ED_OutputType),            INTENT(IN   )  :: y_ED           !< ElastoDyn outputs (need translation displacement on meshes for loads mapping)
   TYPE(AD_OutputType),            INTENT(IN   )  :: y_AD           !< AeroDyn outputs
   TYPE(AD_InputType),             INTENT(INOUT)  :: u_AD           !< AD inputs (for AD-ED load linerization)
   TYPE(BeamDyn_Data),             INTENT(INOUT)  :: BD             !< BD data at t

   TYPE(FAST_ModuleMapType),       INTENT(INOUT)  :: MeshMapData    !< Data for mapping between modules
   REAL(R8Ki),                     INTENT(INOUT)  :: dUdu(:,:)      !< Jacobian matrix of which we are computing the dU^(ED)/du^(AD) block
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat        !< Error status
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg         !< Error message
   
      ! local variables
   INTEGER(IntKi)                                 :: i              ! rows/columns
   INTEGER(IntKi)                                 :: K              ! Loops through blades
   INTEGER(IntKi)                                 :: BD_Start       ! starting index of dUdu (column) where BD root motion inputs are located
   INTEGER(IntKi)                                 :: AD_Start_Bl    ! starting index of dUdu (column) where AD blade motion inputs are located
   INTEGER(IntKi)                                 :: ED_Start_mt    ! starting index of dUdu (row) where ED blade/tower or hub moment inputs are located
   INTEGER(IntKi)                                 :: ErrStat2       ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                           :: ErrMsg2        ! temporary Error message if ErrStat /= ErrID_None
   
   CHARACTER(*), PARAMETER                        :: RoutineName = 'Linear_ED_InputSolve_du'
   
   
      ! Initialize error status
      
   ErrStat = ErrID_None
   ErrMsg = ""

   !..........
   ! dU^{ED}/du^{AD}
   !..........
   IF ( p_FAST%CompAero == Module_AD ) THEN
   
         ! ED inputs on blade from AeroDyn
      IF (p_FAST%CompElast == Module_ED) THEN
         
            ! blades:
         AD_Start_Bl = y_FAST%Lin%Modules(MODULE_AD)%Instance(1)%LinStartIndx(LIN_INPUT_COL) &
                     + u_AD%TowerMotion%NNodes * 9  & ! 3 fields (MASKID_TRANSLATIONDISP,MASKID_Orientation,MASKID_TRANSLATIONVel) with 3 components
                     + u_AD%HubMotion%NNodes   * 9    ! 3 fields (MASKID_TRANSLATIONDISP,MASKID_Orientation,MASKID_RotationVel) with 3 components
   
         do k = 1,size(u_AD%BladeRootMotion)         
            AD_Start_Bl = AD_Start_Bl + u_AD%BladeRootMotion(k)%NNodes * 3 ! 1 field (MASKID_Orientation) with 3 components
         end do
         ! next is u_AD%BladeMotion(k); note that it has 3 fields and we only need 1
      
         ED_Start_mt = y_FAST%Lin%Modules(MODULE_ED)%Instance(1)%LinStartIndx(LIN_INPUT_COL)
         DO K = 1,SIZE(u_ED%BladePtLoads,1) ! Loop through all blades (p_ED%NumBl)
            ED_Start_mt = ED_Start_mt + u_ED%BladePtLoads(k)%NNodes*3 ! skip the forces on this blade
            
            CALL Linearize_Line2_to_Point( y_AD%BladeLoad(k), u_ED%BladePtLoads(k), MeshMapData%AD_L_2_BDED_B(k), ErrStat2, ErrMsg2, u_AD%BladeMotion(k), y_ED%BladeLn2Mesh(k) )
               CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            
               ! AD is source in the mapping, so we want M_{uSm}               
            if (allocated(MeshMapData%AD_L_2_BDED_B(k)%dM%m_us )) then
               call SetBlockMatrix( dUdu, MeshMapData%AD_L_2_BDED_B(k)%dM%m_us, ED_Start_mt, AD_Start_Bl )
            end if
            
               ! get starting index of next blade
            AD_Start_Bl = AD_Start_Bl + u_AD%BladeMotion(k)%Nnodes * 9  ! 3 fields (MASKID_TRANSLATIONDISP,MASKID_Orientation,MASKID_TRANSLATIONVel) with 3 components                         
            ED_Start_mt = ED_Start_mt + u_ED%BladePtLoads(k)%NNodes* 3  ! skip the moments on this blade
               
         END DO

      END IF
      
      ! ED inputs on tower from AD:
      
      IF ( y_AD%TowerLoad%Committed ) THEN
         ED_Start_mt = y_FAST%Lin%Modules(MODULE_ED)%Instance(1)%LinStartIndx(LIN_INPUT_COL)
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
            call SetBlockMatrix( dUdu, MeshMapData%AD_L_2_ED_P_T%dM%m_us, ED_Start_mt, y_FAST%Lin%Modules(MODULE_AD)%Instance(1)%LinStartIndx(LIN_INPUT_COL) )
         end if
      END IF
      
   END IF
   
   
   !..........
   ! dU^{ED}/du^{BD}
   !..........
   
   IF ( p_FAST%CompElast == Module_BD ) THEN ! see routine U_ED_SD_HD_BD_Orca_Residual() in SolveOption1

         ED_Start_mt = y_FAST%Lin%Modules(MODULE_ED)%Instance(1)%LinStartIndx(LIN_INPUT_COL)
         ! u_ED%BladePtLoads(i)%NNodes = 0 here
         ED_Start_mt = ED_Start_mt + u_ED%PlatformPtMesh%NNodes * 6      &      ! 3 forces + 3 moments at each node
                                   + u_ED%TowerPtLoads%NNodes   * 6      &      ! 3 forces + 3 moments at each node
                                   + u_ED%HubPtLoad%NNodes      * 3             ! 3 forces at the hub (so we start at the moments)
   
         ! Transfer BD loads to ED hub input:
         ! we're mapping loads, so we also need the sibling meshes' displacements:
         do k=1,p_FAST%nBeams
            BD_Start = y_FAST%Lin%Modules(MODULE_BD)%Instance(k)%LinStartIndx(LIN_INPUT_COL)
            
            CALL Linearize_Point_to_Point( BD%y(k)%ReactionForce, u_ED%HubPtLoad, MeshMapData%BD_P_2_ED_P(k), ErrStat2, ErrMsg2, BD%Input(1,k)%RootMotion, y_ED%HubPtMotion) !u_BD%RootMotion and y_ED%HubPtMotion contain the displaced positions for load calculations
               CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)
               
               ! BD is source in the mapping, so we want M_{uSm}
            if (allocated(MeshMapData%BD_P_2_ED_P(k)%dM%m_us )) then
               call SetBlockMatrix( dUdu, MeshMapData%BD_P_2_ED_P(k)%dM%m_us, ED_Start_mt, BD_Start )
            end if
               
         end do ! k
   
   END IF
      
END SUBROUTINE Linear_ED_InputSolve_du
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine forms the dU^{BD}/du^{BD} and dU^{BD}/du^{AD} blocks (BD row) of dUdu. (i.e., how do changes in the AD and BD inputs 
!! affect the BD inputs?) This should be called only when p_FAST%CompElast == Module_BD.
SUBROUTINE Linear_BD_InputSolve_du( p_FAST, y_FAST, y_ED, y_AD, u_AD, BD, MeshMapData, dUdu, ErrStat, ErrMsg )

   TYPE(FAST_ParameterType),       INTENT(IN   )  :: p_FAST         !< Glue-code simulation parameters
   TYPE(FAST_OutputFileType),      INTENT(IN   )  :: y_FAST         !< Glue-code output parameters (for linearization)
   TYPE(ED_OutputType),            INTENT(IN   )  :: y_ED           !< ElastoDyn outputs (need translation displacement on meshes for loads mapping)
   TYPE(AD_OutputType),            INTENT(IN   )  :: y_AD           !< AeroDyn outputs
   TYPE(AD_InputType),             INTENT(INOUT)  :: u_AD           !< AD inputs (for AD-ED load linerization)
   TYPE(BeamDyn_Data),             INTENT(INOUT)  :: BD             !< BD data at t

   TYPE(FAST_ModuleMapType),       INTENT(INOUT)  :: MeshMapData    !< Data for mapping between modules
   REAL(R8Ki),                     INTENT(INOUT)  :: dUdu(:,:)      !< Jacobian matrix of which we are computing the dU^(ED)/du^(AD) block
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat        !< Error status
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg         !< Error message
   
      ! local variables
   INTEGER(IntKi)                                 :: k              ! Loops through blades
   INTEGER(IntKi)                                 :: BD_Start       ! starting index of dUdu (row) where BD inputs are located
   INTEGER(IntKi)                                 :: AD_Start       ! starting index of dUdu (column) where AD inputs are located
   INTEGER(IntKi)                                 :: ErrStat2       ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                           :: ErrMsg2        ! temporary Error message if ErrStat /= ErrID_None
   
   CHARACTER(*), PARAMETER                        :: RoutineName = 'Linear_BD_InputSolve_du'
   
   
      ! Initialize error status
      
   ErrStat = ErrID_None
   ErrMsg = ""

   !..........
   ! dU^{BD}/du^{AD}
   !..........
   IF ( p_FAST%CompAero == Module_AD ) THEN
   
      ! BD inputs on blade from AeroDyn
   
         
      if (p_FAST%BD_OutputSibling) then
            
         DO K = 1,p_FAST%nBeams ! Loop through all blades
            CALL Linearize_Line2_to_Line2( y_AD%BladeLoad(k), BD%Input(1,k)%DistrLoad, MeshMapData%AD_L_2_BDED_B(k), ErrStat2, ErrMsg2, u_AD%BladeMotion(k), BD%y(k)%BldMotion )
               CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         END DO
         
      else
      
         DO K = 1,p_FAST%nBeams ! Loop through all blades
            !linearization for dUdy will need some matrix multiplies because of the transfer (chain rule!), but we will perform individual linearization calculations here
            !!! need to transfer the BD output blade motions to nodes on a sibling of the BD blade motion mesh:
            CALL Linearize_Line2_to_Line2( BD%y(k)%BldMotion, MeshMapData%y_BD_BldMotion_4Loads(k), MeshMapData%BD_L_2_BD_L(k), ErrStat2, ErrMsg2 )
               CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

            CALL Linearize_Line2_to_Line2( y_AD%BladeLoad(k), BD%Input(1,k)%DistrLoad, MeshMapData%AD_L_2_BDED_B(k), ErrStat2, ErrMsg2, u_AD%BladeMotion(k), MeshMapData%y_BD_BldMotion_4Loads(k) )
               CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         END DO
         
      end if

      
      DO K = 1,p_FAST%nBeams ! Loop through all blades
         
            ! AD is source in the mapping, so we want M_{uSm}
         if (allocated(MeshMapData%AD_L_2_BDED_B(k)%dM%m_us )) then
            AD_Start = Indx_u_AD_Blade_Start(u_AD, y_FAST, k) ! index for the start of u_AD%BladeMotion(k)%translationDisp field
         
            BD_Start = y_FAST%Lin%Modules(MODULE_BD)%Instance(k)%LinStartIndx(LIN_INPUT_COL) &
                     + BD%Input(1,k)%RootMotion%NNodes *18  & ! displacement, rotation, & acceleration fields for each node
                     + BD%Input(1,k)%PointLoad%NNodes  * 6  & ! force + moment fields for each node
                     + BD%Input(1,k)%DistrLoad%NNodes  * 3    ! force field for each node (start with moment field)
                        
            call SetBlockMatrix( dUdu, MeshMapData%AD_L_2_BDED_B(k)%dM%m_us, BD_Start, AD_Start )
         end if
               
      END DO

   END IF
   
   !..........
   ! dU^{BD}/du^{BD}
   ! note that the 1s on the diagonal have already been set, so we will fill in the off diagonal terms.
   !..........
   
   !IF ( p_FAST%CompElast == Module_BD ) THEN ! see routine U_ED_SD_HD_BD_Orca_Residual() in SolveOption1

      ! Transfer ED motions to BD motion input (BD inputs depend on previously calculated BD inputs from ED):
      do k=1,p_FAST%nBeams
            
         call Linearize_Point_to_Point( y_ED%BladeRootMotion(k), BD%Input(1,k)%RootMotion, MeshMapData%ED_P_2_BD_P(k), ErrStat2, ErrMsg2 )
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

            ! BD is destination in the mapping, so we want M_{tv_uD} and M_{ta_uD}
            ! translational velocity:
         if (allocated(MeshMapData%ED_P_2_BD_P(k)%dM%tv_uD )) then
            BD_Start = y_FAST%Lin%Modules(MODULE_BD)%Instance(k)%LinStartIndx(LIN_INPUT_COL) + BD%Input(1,k)%RootMotion%NNodes * 6 ! skip root translational displacement and orientation fields
            call SetBlockMatrix( dUdu, MeshMapData%ED_P_2_BD_P(k)%dM%tv_uD, BD_Start, y_FAST%Lin%Modules(MODULE_BD)%Instance(k)%LinStartIndx(LIN_INPUT_COL) )
         end if

            ! translational acceleration:
         if (allocated(MeshMapData%ED_P_2_BD_P(k)%dM%ta_uD )) then
            BD_Start = y_FAST%Lin%Modules(MODULE_BD)%Instance(k)%LinStartIndx(LIN_INPUT_COL) + BD%Input(1,k)%RootMotion%NNodes * 12 ! skip root translational displacement, orientation, and velocity (translation and rotation) fields

            call SetBlockMatrix( dUdu, MeshMapData%ED_P_2_BD_P(k)%dM%ta_uD, BD_Start, y_FAST%Lin%Modules(MODULE_BD)%Instance(k)%LinStartIndx(LIN_INPUT_COL) )
         end if

      end do ! k
   
   !END IF
      
END SUBROUTINE Linear_BD_InputSolve_du
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine forms the dU^{AD}/du^{AD} block of dUdu. (i.e., how do changes in the AD inputs affect the AD inputs?)
SUBROUTINE Linear_AD_InputSolve_du( p_FAST, y_FAST, u_AD, y_ED, BD, MeshMapData, dUdu, ErrStat, ErrMsg )

      ! Passed variables
   TYPE(FAST_ParameterType),    INTENT(IN   )   :: p_FAST      !< FAST parameter data    
   TYPE(FAST_OutputFileType),   INTENT(IN   )   :: y_FAST      !< FAST output file data (for linearization)
   TYPE(AD_InputType),          INTENT(INOUT)   :: u_AD        !< The inputs to AeroDyn14
   TYPE(ED_OutputType),         INTENT(IN)      :: y_ED        !< The outputs from the structural dynamics module
   TYPE(BeamDyn_Data),          INTENT(INOUT)   :: BD          !< BD data at t
   TYPE(FAST_ModuleMapType),    INTENT(INOUT)   :: MeshMapData !< Data for mapping between modules
   REAL(R8Ki),                  INTENT(INOUT)   :: dUdu(:,:)   !< Jacobian matrix of which we are computing the dU^(ED)/du^(AD) block
   
   INTEGER(IntKi)                               :: ErrStat     !< Error status of the operation
   CHARACTER(*)                                 :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      ! Local variables:

   INTEGER(IntKi)                               :: K              ! Loops through blades   
   INTEGER(IntKi)                               :: AD_Start_td    ! starting index of dUdu (column) where AD translation displacements are located   
   INTEGER(IntKi)                               :: AD_Start_tv    ! starting index of dUdu (column) where AD translation velocities are located   
   INTEGER(IntKi)                               :: ErrStat2
   CHARACTER(ErrMsgLen)                         :: ErrMsg2 
   CHARACTER(*), PARAMETER                      :: RoutineName = 'Linear_AD_InputSolve_du'

   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   ! note that we assume this block matrix has been initialized to the identity matrix before calling this routine
   
   ! look at how the translational displacement gets transfered to the translational velocity:
   !-------------------------------------------------------------------------------------------------
   ! Set the inputs from ElastoDyn and/or BeamDyn:
   !-------------------------------------------------------------------------------------------------
      
      ! tower
   IF (u_AD%TowerMotion%Committed) THEN
            
      CALL Linearize_Line2_to_Line2( y_ED%TowerLn2Mesh, u_AD%TowerMotion, MeshMapData%ED_L_2_AD_L_T, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//':u_AD%TowerMotion' )     

      AD_Start_td = y_FAST%Lin%Modules(MODULE_AD)%Instance(1)%LinStartIndx(LIN_INPUT_COL)
      AD_Start_tv = AD_Start_td + u_AD%TowerMotion%NNodes * 6 ! 2 fields (TranslationDisp and Orientation) with 3 components before translational velocity field      

      !AD is the destination here, so we need tv_ud
      if (allocated( MeshMapData%ED_L_2_AD_L_T%dM%tv_ud)) then
         call SetBlockMatrix( dUdu, MeshMapData%ED_L_2_AD_L_T%dM%tv_ud, AD_Start_tv, AD_Start_td )
      end if
               
            
   END IF
   
      
      ! blades
   IF (p_FAST%CompElast == Module_ED ) THEN
      
      DO k=1,size(u_AD%BladeMotion)
         CALL Linearize_Line2_to_Line2( y_ED%BladeLn2Mesh(k), u_AD%BladeMotion(k), MeshMapData%BDED_L_2_AD_L_B(k), ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//':u_AD%BladeMotion('//trim(num2lstr(k))//')' )
      END DO
      
   ELSEIF (p_FAST%CompElast == Module_BD ) THEN
   
      DO k=1,size(u_AD%BladeMotion)
         CALL Linearize_Line2_to_Line2( BD%y(k)%BldMotion, u_AD%BladeMotion(k), MeshMapData%BDED_L_2_AD_L_B(k), ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//':u_AD%BladeMotion('//trim(num2lstr(k))//')' )
      END DO
         
   END IF
   
   
      ! index for u_AD%BladeMotion(1)%translationDisp field
   AD_Start_td = y_FAST%Lin%Modules(MODULE_AD)%Instance(1)%LinStartIndx(LIN_INPUT_COL) &
               + u_AD%TowerMotion%NNodes * 9  & ! 3 fields (MASKID_TRANSLATIONDISP,MASKID_Orientation,MASKID_TRANSLATIONVel) with 3 components
               + u_AD%HubMotion%NNodes   * 9    ! 3 fields (MASKID_TRANSLATIONDISP,MASKID_Orientation,MASKID_RotationVel) with 3 components
   
   do k = 1,size(u_AD%BladeRootMotion)
      AD_Start_td = AD_Start_td + u_AD%BladeRootMotion(k)%NNodes * 3 ! 1 field (MASKID_Orientation) with 3 components
   end do
   
   DO k=1,size(u_AD%BladeMotion)

         !AD is the destination here, so we need tv_ud
      if (allocated( MeshMapData%BDED_L_2_AD_L_B(k)%dM%tv_ud)) then
            ! index for u_AD%BladeMotion(k+1)%translationVel field
         AD_Start_tv = AD_Start_td + u_AD%BladeMotion(k)%NNodes * 6 ! 2 fields (TranslationDisp and Orientation) with 3 components before translational velocity field      

         call SetBlockMatrix( dUdu, MeshMapData%BDED_L_2_AD_L_B(k)%dM%tv_ud, AD_Start_tv, AD_Start_td )
      end if

         ! index for u_AD%BladeMotion(k+1)%translationDisp field
      AD_Start_td = AD_Start_td + u_AD%BladeMotion(k)%NNodes * 9 ! 3 fields (TranslationDisp, Orientation, TranslationVel) with 3 components

   END DO
      
   
   
END SUBROUTINE Linear_AD_InputSolve_du
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine forms the dU^{SrvD}/dy^{ED} block of dUdy. (i.e., how do changes in the ED outputs affect the SrvD inputs?)
SUBROUTINE Linear_SrvD_InputSolve_dy( p_FAST, y_FAST, dUdy  )
!..................................................................................................................................

   TYPE(FAST_ParameterType),       INTENT(IN)     :: p_FAST         !< Glue-code simulation parameters
   TYPE(FAST_OutputFileType),      INTENT(IN)     :: y_FAST         !< Output variables for the glue code
   REAL(R8Ki),                     INTENT(INOUT)  :: dUdy(:,:)      !< Jacobian matrix of which we are computing the dU^{SrvD}/dy^{ED} block
   
   integer(intKi)                                 :: ED_Start_Yaw   !< starting index of dUdy (column) where ED Yaw/YawRate/HSS_Spd outputs are located (just before WriteOutput)

   
   
   INTEGER(IntKi)                                   :: i            ! loop counter
   
   CHARACTER(*), PARAMETER                          :: RoutineName = 'Linear_SrvD_InputSolve_dy'
   
   ED_Start_Yaw = y_FAST%Lin%Modules(Module_ED)%Instance(1)%LinStartIndx(LIN_OUTPUT_COL) + y_FAST%Lin%Modules(Module_ED)%Instance(1)%SizeLin(LIN_OUTPUT_COL) & !end of ED outputs (+1)
                - y_FAST%Lin%Modules(Module_ED)%Instance(1)%NumOutputs - 3 ! start of ED where Yaw, YawRate, HSS_Spd occur (right before WriteOutputs)
   do i=1,3
      dUdy(y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%LinStartIndx(LIN_INPUT_COL) + i - 1, ED_Start_Yaw + i - 1) = -1.0_ReKi
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
                                    
END SUBROUTINE Linear_SrvD_InputSolve_dy
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine forms the dU^{ED}/dy^{SrvD}, dU^{ED}/dy^{ED}, dU^{ED}/dy^{BD}, and dU^{ED}/dy^{AD} blocks of dUdy. (i.e., how do 
!! changes in the SrvD, ED, BD, and AD outputs effect the ED inputs?)
SUBROUTINE Linear_ED_InputSolve_dy( p_FAST, y_FAST, u_ED, y_ED, y_AD, u_AD, BD, MeshMapData, dUdy, ErrStat, ErrMsg )

   TYPE(FAST_ParameterType),       INTENT(IN   )  :: p_FAST           !< Glue-code simulation parameters
   TYPE(FAST_OutputFileType),      INTENT(IN   )  :: y_FAST           !< FAST output file data (for linearization)
   TYPE(ED_InputType),             INTENT(INOUT)  :: u_ED             !< ED Inputs at t
   TYPE(ED_OutputType),            INTENT(IN   )  :: y_ED             !< ElastoDyn outputs (need translation displacement on meshes for loads mapping)
   TYPE(AD_OutputType),            INTENT(IN   )  :: y_AD             !< AeroDyn outputs
   TYPE(AD_InputType),             INTENT(INOUT)  :: u_AD             !< AD inputs (for AD-ED load linerization)
   TYPE(BeamDyn_Data),             INTENT(INOUT)  :: BD               !< BD data at t
                                                                      
   TYPE(FAST_ModuleMapType),       INTENT(INOUT)  :: MeshMapData      !< Data for mapping between modules
   REAL(R8Ki),                     INTENT(INOUT)  :: dUdy(:,:)        !< Jacobian matrix of which we are computing the dU^(ED)/du^(AD) block
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat          !< Error status
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg           !< Error message
   
      ! local variables
   INTEGER(IntKi)                                 :: i                ! rows/columns
   INTEGER(IntKi)                                 :: K                ! Loops through blades
   INTEGER(IntKi)                                 :: AD_Out_Start     ! starting index of dUdy (column) where particular AD fields are located
   INTEGER(IntKi)                                 :: BD_Out_Start     ! starting index of dUdy (column) where particular BD fields are located
   INTEGER(IntKi)                                 :: ED_Start         ! starting index of dUdy (row) where ED input fields are located
   INTEGER(IntKi)                                 :: ED_Out_Start     ! starting index of dUdy (column) where ED output fields are located
   
   CHARACTER(*), PARAMETER                        :: RoutineName = 'Linear_ED_InputSolve_dy' 
   
   
      ! Initialize error status
      
   ErrStat = ErrID_None
   ErrMsg = ""
   
      ! ED inputs from ServoDyn outputs
   IF ( p_FAST%CompServo == Module_SrvD ) THEN

         ! BlPitchCom, YawMom, GenTrq
      ED_Start = Indx_u_ED_BlPitchCom_Start(u_ED, y_FAST)
      do i=1,size(u_ED%BlPitchCom)+2 ! BlPitchCom, YawMom, GenTrq (NOT collective pitch)
         dUdy(ED_Start + i - 1, y_FAST%Lin%Modules(Module_SrvD)%Instance(1)%LinStartIndx(LIN_OUTPUT_COL) + i - 1) = -1.0_ReKi
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

   ! parts of dU^{ED}/dy^{AD} and dU^{ED}/dy^{ED}:
   
      ! ElastoDyn inputs on blade from AeroDyn and ElastoDyn
   IF ( p_FAST%CompAero == Module_AD ) THEN

      IF (p_FAST%CompElast == Module_ED) THEN 
         AD_Out_Start = y_FAST%Lin%Modules(MODULE_AD)%Instance(1)%LinStartIndx(LIN_OUTPUT_COL) + y_AD%TowerLoad%NNodes * 6    ! start of y_AD%BladeLoad(1)%Force field [2 fields (force, moment) with 3 components]
         
         DO K = 1,SIZE(u_ED%BladePtLoads,1) ! Loop through all blades (p_ED%NumBl)
            !!! ! This linearization was done in forming dUdu (see Linear_ED_InputSolve_du()), so we don't need to re-calculate these matrices 
            !!! ! while forming dUdy, too.
            !CALL Linearize_Line2_to_Point( y_AD%BladeLoad(k), u_ED%BladePtLoads(k), MeshMapData%AD_L_2_BDED_B(k), ErrStat2, ErrMsg2, u_AD%BladeMotion(k), y_ED%BladeLn2Mesh(k) )
               
               ! AD loads-to-ED loads transfer (dU^{ED}/dy^{AD}):
            ED_Start = Indx_u_ED_Blade_Start(u_ED, y_FAST, k) ! start of u_ED%BladePtLoads(k)%Force field
            call Assemble_dUdy_Loads(y_AD%BladeLoad(k), u_ED%BladePtLoads(k), MeshMapData%AD_L_2_BDED_B(k), ED_Start, AD_Out_Start, dUdy)

               ! ED translation displacement-to-ED moment transfer (dU^{ED}/dy^{ED}):
            ED_Start = Indx_u_ED_Blade_Start(u_ED, y_FAST, k) + u_ED%BladePtLoads(k)%NNodes*3   ! start of u_ED%BladePtLoads(k)%Moment field (skip the ED forces)
            ED_Out_Start = Indx_y_ED_Blade_Start(y_ED, y_FAST, k) ! start of y_ED%BladeLn2Mesh(1)%TranslationDisp field
            call SetBlockMatrix( dUdy, MeshMapData%AD_L_2_BDED_B(k)%dM%m_uD, ED_Start, ED_Out_Start )

            AD_Out_Start = AD_Out_Start + y_AD%BladeLoad(k)%NNodes*6        ! start of y_AD%BladeLoad(k+1)%Force field [skip 2 fields to forces on next blade]
         END DO
      END IF ! ED
      
      
      IF ( y_AD%TowerLoad%Committed ) THEN
         !!! ! This linearization was done in forming dUdu (see Linear_ED_InputSolve_du()), so we don't need to re-calculate these matrices 
         !!! ! while forming dUdy, too.
         !CALL Linearize_Line2_to_Point( y_AD%TowerLoad, u_ED%TowerPtLoads, MeshMapData%AD_L_2_ED_P_T, ErrStat2, ErrMsg2, u_AD%TowerMotion, y_ED%TowerLn2Mesh )
            
            ! AD loads-to-ED loads transfer (dU^{ED}/dy^{AD}):
         ED_Start = Indx_u_ED_Tower_Start(u_ED, y_FAST) ! u_ED%TowerPtLoads%Force field
         AD_Out_Start = y_FAST%Lin%Modules(MODULE_AD)%Instance(1)%LinStartIndx(LIN_OUTPUT_COL) ! start of y_AD%Tower%Force
         call Assemble_dUdy_Loads(y_AD%TowerLoad, u_ED%TowerPtLoads, MeshMapData%AD_L_2_ED_P_T, ED_Start, AD_Out_Start, dUdy)

            ! ED translation displacement-to-ED moment transfer (dU^{ED}/dy^{ED}):
         ED_Start = ED_Start + u_ED%TowerPtLoads%NNodes*3 ! start of u_ED%TowerPtLoads%Moment field  [skip the ED forces to get to the moments]
         ED_Out_Start  = Indx_y_ED_Tower_Start(y_ED, y_FAST) ! start of y_ED%TowerLn2Mesh%TranslationDisp field
         call SetBlockMatrix( dUdy, MeshMapData%AD_L_2_ED_P_T%dM%m_uD, ED_Start, ED_Out_Start )
            
      END IF ! tower
      
   END IF ! aero loads
      
      ! U_ED_SD_HD_BD_Orca_Residual() in InputSolve Option 1
   IF (p_FAST%CompElast == Module_BD) THEN
   
      !!! ! This linearization was done in forming dUdu (see Linear_ED_InputSolve_du()), so we don't need to re-calculate these matrices 
      !!! ! while forming dUdy, too.
      !!!DO k=1,p_FAST%nBeams
      !!!   CALL Linearize_Point_to_Point( BD%y(k)%ReactionForce, u_ED%HubPtLoad, MeshMapData%BD_P_2_ED_P(k), ErrStat2, ErrMsg2, BD%Input(1,k)%RootMotion, y_ED%HubPtMotion)
      !!!END DO

         ! BD Reaction force-to-ED force transfer (dU^{ED}/dy^{BD}) from BD root-to-ED hub load transfer:
      ED_Start = Indx_u_ED_Hub_Start(u_ED, y_FAST) ! start of u_ED%HubPtLoad%Force field
      DO k=1,p_FAST%nBeams
         BD_Out_Start = y_FAST%Lin%Modules(MODULE_BD)%Instance(k)%LinStartIndx(LIN_OUTPUT_COL) ! BD%y(k)%ReactionForce%Force field
         call Assemble_dUdy_Loads(BD%y(k)%ReactionForce, u_ED%HubPtLoad, MeshMapData%BD_P_2_ED_P(k), ED_Start, BD_Out_Start, dUdy)
      END DO

         ! ED translation displacement-to-ED moment transfer (dU^{ED}/dy^{ED}) from BD root-to-ED hub load transfer:
      ED_Start = Indx_u_ED_Hub_Start(u_ED, y_FAST) + u_ED%HubPtLoad%NNodes*3 ! start of u_ED%HubPtLoad%Moment field (skip forces)
      DO k=1,p_FAST%nBeams
         ED_Out_Start = Indx_y_ED_BladeRoot_Start(y_ED, y_FAST, k) ! start of y_ED%BladeRootMotion(k)%TranslationDisp field
         CALL SumBlockMatrix( dUdy, MeshMapData%BD_P_2_ED_P(k)%dM%m_ud, ED_Start, ED_Out_Start)
      END DO

   END IF

END SUBROUTINE Linear_ED_InputSolve_dy
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine forms the dU^{BD}/dy^{ED}, dU^{BD}/dy^{BD}, and dU^{BD}/dy^{AD} blocks of dUdy. (i.e., how do 
!! changes in the ED, BD, and AD outputs effect the BD inputs?)
SUBROUTINE Linear_BD_InputSolve_dy( p_FAST, y_FAST, u_ED, y_ED, y_AD, u_AD, BD, MeshMapData, dUdy, ErrStat, ErrMsg )

   TYPE(FAST_ParameterType),       INTENT(IN   )  :: p_FAST           !< Glue-code simulation parameters
   TYPE(FAST_OutputFileType),      INTENT(IN   )  :: y_FAST           !< FAST output file data (for linearization)
   TYPE(ED_InputType),             INTENT(INOUT)  :: u_ED             !< ED Inputs at t
   TYPE(ED_OutputType),            INTENT(IN   )  :: y_ED             !< ElastoDyn outputs (need translation displacement on meshes for loads mapping)
   TYPE(AD_OutputType),            INTENT(IN   )  :: y_AD             !< AeroDyn outputs
   TYPE(AD_InputType),             INTENT(INOUT)  :: u_AD             !< AD inputs (for AD-ED load linerization)
   TYPE(BeamDyn_Data),             INTENT(IN   )  :: BD               !< BD data at t
                                                                      
   TYPE(FAST_ModuleMapType),       INTENT(INOUT)  :: MeshMapData      !< Data for mapping between modules
   REAL(R8Ki),                     INTENT(INOUT)  :: dUdy(:,:)        !< Jacobian matrix of which we are computing the dU^(ED)/du^(AD) block
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat          !< Error status
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg           !< Error message
   
      ! local variables
   INTEGER(IntKi)                                 :: K                ! Loops through blades
   INTEGER(IntKi)                                 :: AD_Out_Start     ! starting index of dUdy (column) where particular AD fields are located
   INTEGER(IntKi)                                 :: BD_Start         ! starting index of dUdy (column) where particular BD fields are located
   INTEGER(IntKi)                                 :: BD_Out_Start     ! starting index of dUdy (column) where BD output fields are located
   INTEGER(IntKi)                                 :: ED_Out_Start     ! starting index of dUdy (column) where particular AD fields are located
   INTEGER(IntKi)                                 :: ErrStat2
   CHARACTER(ErrMsgLen)                           :: ErrMsg2
   REAL(R8Ki), ALLOCATABLE                        :: TempMat(:,:)     ! temporary matrix for getting linearization matrices when BD input and output meshes are not siblings
   CHARACTER(*), PARAMETER                        :: RoutineName = 'Linear_BD_InputSolve_dy' 
   
   
      ! Initialize error status
      
   ErrStat = ErrID_None
   ErrMsg = ""

   ! parts of dU^{BD}/dy^{AD} and dU^{BD}/dy^{BD}:
   
      ! BeamDyn inputs on blade from AeroDyn and BeamDyn
   IF ( p_FAST%CompAero == Module_AD ) THEN

      !!! ! This linearization was done in forming dUdu (see Linear_BD_InputSolve_du()), so we don't need to re-calculate these matrices 
      !!! ! while forming dUdy, too.
      !!!if (p_FAST%BD_OutputSibling) then
      !!!   CALL Linearize_Line2_to_Line2( y_AD%BladeLoad(k), BD%Input(1,k)%DistrLoad, MeshMapData%AD_L_2_BDED_B(k), ErrStat2, ErrMsg2, u_AD%BladeMotion(k), BD%y(k)%BldMotion )
      !!!else
      !!!   CALL Linearize_Line2_to_Line2( BD%y(k)%BldMotion, MeshMapData%y_BD_BldMotion_4Loads(k), MeshMapData%BD_L_2_BD_L(k), ErrStat2, ErrMsg2 )
      !!!   CALL Linearize_Line2_to_Line2( y_AD%BladeLoad(k), BD%Input(1,k)%DistrLoad, MeshMapData%AD_L_2_BDED_B(k), ErrStat2, ErrMsg2, u_AD%BladeMotion(k), MeshMapData%y_BD_BldMotion_4Loads(k) )
      !!!end if

      AD_Out_Start = y_FAST%Lin%Modules(MODULE_AD)%Instance(1)%LinStartIndx(LIN_OUTPUT_COL) + y_AD%TowerLoad%NNodes * 6    ! start of y_AD%BladeLoad(1)%Force field [2 fields (force, moment) with 3 components]
      DO K = 1,p_FAST%nBeams ! Loop through all blades

         BD_Start = y_FAST%Lin%Modules(MODULE_BD)%Instance(k)%LinStartIndx(LIN_INPUT_COL) & ! start of BD%Input(1,k)%DistrLoad%Force field
                  + BD%Input(1,k)%RootMotion%NNodes *18  & ! displacement, rotation, & acceleration fields for each node
                  + BD%Input(1,k)%PointLoad%NNodes  * 6    ! force + moment fields for each node
            
            ! AD loads-to-BD loads transfer (dU^{BD}/dy^{AD}):
         call Assemble_dUdy_Loads(y_AD%BladeLoad(k), BD%Input(1,k)%DistrLoad, MeshMapData%AD_L_2_BDED_B(k), BD_Start, AD_Out_Start, dUdy)
         AD_Out_Start = AD_Out_Start + y_AD%BladeLoad(k)%NNodes*6  ! start of y_AD%BladeLoad(k+1)%Force field [skip the moments to get to forces on next blade]
         
         
            ! BD translation displacement-to-BD moment transfer (dU^{BD}/dy^{BD}):
         BD_Start = BD_Start + BD%Input(1,k)%DistrLoad%NNodes  * 3    ! start of BD%Input(1,k)%DistrLoad%Moment field (start with moment field)
         BD_Out_Start  = y_FAST%Lin%Modules(MODULE_BD)%Instance(k)%LinStartIndx(LIN_OUTPUT_COL) & ! start of BD%y(k)%BldMotion%TranslationDisp field
                       + BD%y(k)%ReactionForce%NNodes * 6 ! 2 fields with 3 components

         if (p_FAST%BD_OutputSibling) then
            call SetBlockMatrix( dUdy, MeshMapData%AD_L_2_BDED_B(k)%dM%m_uD, BD_Start, BD_Out_Start )
         else
            call AllocAry(TempMat, size(MeshMapData%AD_L_2_BDED_B(k)%dM%m_uD,1), size(MeshMapData%BD_L_2_BD_L(k)%dM%mi,2), 'TempMat', ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
               if (ErrStat>=AbortErrLev) return
            
                  ! these blocks should be small enough that we can use matmul instead of calling a LAPACK routine to do it.
            TempMat = matmul(MeshMapData%AD_L_2_BDED_B(k)%dM%m_uD,MeshMapData%BD_L_2_BD_L(k)%dM%mi)
            call SetBlockMatrix( dUdy, TempMat, BD_Start, BD_Out_Start )
            
            BD_Out_Start = BD_Out_Start + BD%y(k)%BldMotion%NNodes*3 ! start of BD%y(k)%BldMotion%Orientation field
            TempMat = matmul(MeshMapData%AD_L_2_BDED_B(k)%dM%m_uD,MeshMapData%BD_L_2_BD_L(k)%dM%fx_p)
            call SetBlockMatrix( dUdy, TempMat, BD_Start, BD_Out_Start )

            deallocate(TempMat) ! the next blade may have a different number of nodes
         end if

      END DO
            
   END IF ! aero loads
      
   ! U_ED_SD_HD_BD_Orca_Residual() in InputSolve Option 1; call to Transfer_ED_to_BD_tmp()
   !IF ( p_FAST%CompElast == Module_BD .and. BD_Solve_Option1) THEN
         
   ! Transfer ED motions to BD inputs (dU^{BD}/dy^{ED}):
   do k = 1,size(y_ED%BladeRootMotion)
      !!! ! This linearization was done in forming dUdu (see Linear_BD_InputSolve_du()), so we don't need to re-calculate these matrices 
      !!! ! while forming dUdy, too.
      !!!CALL Linearize_Point_to_Point( y_ED%BladeRootMotion(k), BD%Input(1,k)%RootMotion, MeshMapData%ED_P_2_BD_P(k), ErrStat2, ErrMsg2 )
   
      BD_Start = y_FAST%Lin%Modules(MODULE_BD)%Instance(k)%LinStartIndx(LIN_INPUT_COL) ! ! start of BD%Input(1,k)%RootMotion%TranslationDisp field
      ED_Out_Start = Indx_y_ED_BladeRoot_Start(y_ED, y_FAST, k) ! start of y_ED%BladeRootMotion(k)%TranslationDisp field
      
      call Assemble_dUdy_Motions(y_ED%BladeRootMotion(k), BD%Input(1,k)%RootMotion, MeshMapData%ED_P_2_BD_P(k), BD_Start, ED_Out_Start, dUdy)
   end do

END SUBROUTINE Linear_BD_InputSolve_dy
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine forms the dU^{AD}/dy^{IfW} block of dUdy. (i.e., how do changes in the IfW outputs affect the AD inputs?)
SUBROUTINE Linear_AD_InputSolve_IfW_dy( p_FAST, y_FAST, u_AD, dUdy )

      ! Passed variables
   TYPE(FAST_ParameterType),       INTENT(IN   )  :: p_FAST         !< FAST parameter data    
   TYPE(FAST_OutputFileType),      INTENT(IN   )  :: y_FAST         !< FAST output file data (for linearization)
   TYPE(AD_InputType),             INTENT(INOUT)  :: u_AD           !< The inputs to AeroDyn
   REAL(R8Ki),                     INTENT(INOUT)  :: dUdy(:,:)      !< Jacobian matrix of which we are computing the dU^{AD}/dy^{IfW} block

      ! Local variables:

   INTEGER(IntKi)                               :: I           ! Loops through components
   INTEGER(IntKi)                               :: J           ! Loops through nodes / elements
   INTEGER(IntKi)                               :: K           ! Loops through blades
   INTEGER(IntKi)                               :: node
   INTEGER(IntKi)                               :: AD_Start   ! starting index of dUdy (row) where AD input equations (for specific fields) are located   

                  
   !-------------------------------------------------------------------------------------------------
   ! Set the inputs from inflow wind:
   !-------------------------------------------------------------------------------------------------
   !IF (p_FAST%CompInflow == MODULE_IfW) THEN !already checked in calling routine

      if (p_FAST%CompServo == MODULE_SrvD) then
         node = 2
      else
         node = 1
      end if
      
      
      AD_Start = Indx_u_AD_BladeInflow_Start(u_AD, y_FAST) ! start of u_AD%InflowOnBlade array
      
      do k=1,size(u_AD%InflowOnBlade,3) ! blades
         do j=1,size(u_AD%InflowOnBlade,2) ! nodes
            do i=1,3 !velocity component
               dUdy( AD_Start + i - 1, y_FAST%Lin%Modules(MODULE_IfW)%Instance(1)%LinStartIndx(LIN_OUTPUT_COL) + (node-1)*3 + i - 1 ) = -1.0_R8Ki
            end do
            node = node + 1
            AD_Start = AD_Start + 3
         end do         
      end do
                  
      if ( allocated(u_AD%InflowOnTower) ) then         
         do j=1,size(u_AD%InflowOnTower,2) !nodes
            do i=1,3 !velocity component
               dUdy( AD_Start + i - 1, y_FAST%Lin%Modules(MODULE_IfW)%Instance(1)%LinStartIndx(LIN_OUTPUT_COL) + (node-1)*3 + i - 1 ) = -1.0_R8Ki
            end do
            node = node + 1
            AD_Start = AD_Start + 3
         end do      
      end if
                     
   !END IF
   
   
END SUBROUTINE Linear_AD_InputSolve_IfW_dy
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine forms the dU^{AD}/dy^{ED} and dU^{AD}/dy^{BD} blocks of dUdy. (i.e., how do changes in the ED and BD outputs affect 
!! the AD inputs?)
SUBROUTINE Linear_AD_InputSolve_NoIfW_dy( p_FAST, y_FAST, u_AD, y_ED, BD, MeshMapData, dUdy, ErrStat, ErrMsg )

      ! Passed variables
   TYPE(FAST_ParameterType),    INTENT(IN   )   :: p_FAST      !< FAST parameter data    
   TYPE(FAST_OutputFileType),   INTENT(IN   )   :: y_FAST      !< FAST output file data (for linearization)
   TYPE(AD_InputType),          INTENT(INOUT)   :: u_AD        !< The inputs to AeroDyn14
   TYPE(ED_OutputType),         INTENT(IN)      :: y_ED        !< The outputs from the structural dynamics module
   TYPE(BeamDyn_Data),          INTENT(IN   )   :: BD          !< BD data at t
   TYPE(FAST_ModuleMapType),    INTENT(INOUT)   :: MeshMapData !< Data for mapping between modules
   REAL(R8Ki),                  INTENT(INOUT)   :: dUdy(:,:)   !< Jacobian matrix of which we are computing the dU^{AD}/dy^{ED} block
   
   INTEGER(IntKi)                               :: ErrStat     !< Error status of the operation
   CHARACTER(*)                                 :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      ! Local variables:

   INTEGER(IntKi)                               :: K           ! Loops through blades
   INTEGER(IntKi)                               :: AD_Start    ! starting index of dUdy (column) where particular AD fields are located
   INTEGER(IntKi)                               :: ED_Out_Start! starting index of dUdy (row) where particular ED fields are located
   INTEGER(IntKi)                               :: BD_Out_Start! starting index of dUdy (row) where particular BD fields are located
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
   IF (u_AD%TowerMotion%Committed) THEN
            
      !!! ! This linearization was done in forming dUdu (see Linear_AD_InputSolve_du()), so we don't need to re-calculate these matrices 
      !!! ! while forming dUdy, too.
      !!!CALL Linearize_Line2_to_Line2( y_ED%TowerLn2Mesh, u_AD%TowerMotion, MeshMapData%ED_L_2_AD_L_T, ErrStat2, ErrMsg2 )
      
      AD_Start = Indx_u_AD_Tower_Start(u_AD, y_FAST) ! start of u_AD%TowerMotion%TranslationDisp field
      
      ED_Out_Start = Indx_y_ED_Tower_Start(y_ED, y_FAST) ! start of y_ED%TowerLn2Mesh%TranslationDisp field
      call Assemble_dUdy_Motions(y_ED%TowerLn2Mesh, u_AD%TowerMotion, MeshMapData%ED_L_2_AD_L_T, AD_Start, ED_Out_Start, dUdy, .true.)
   END IF
      
      !...................................
      ! hub
      !...................................
   CALL Linearize_Point_to_Point( y_ED%HubPtMotion, u_AD%HubMotion, MeshMapData%ED_P_2_AD_P_H, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//':u_AD%HubMotion' )      
      if (errStat>=AbortErrLev) return

   ! *** AD translational displacement: from ED translational displacement (MeshMapData%ED_P_2_AD_P_H%dM%mi) and orientation (MeshMapData%ED_P_2_AD_P_H%dM%fx_p)                              
   AD_Start = Indx_u_AD_Hub_Start(u_AD, y_FAST) ! start of u_AD%HubMotion%TranslationDisp field   
   ED_Out_Start = Indx_y_ED_Hub_Start(y_ED, y_FAST) ! start of y_ED%HubPtMotion%TranslationDisp field
   call SetBlockMatrix( dUdy, MeshMapData%ED_P_2_AD_P_H%dM%mi, AD_Start, ED_Out_Start )
   
   ED_Out_Start = Indx_y_ED_Hub_Start(y_ED, y_FAST) + y_ED%HubPtMotion%NNodes * 3 ! start of y_ED%HubPtMotion%Orientation field
   call SetBlockMatrix( dUdy, MeshMapData%ED_P_2_AD_P_H%dM%fx_p, AD_Start, ED_Out_Start )
         
   ! *** AD orientation: from ED orientation
   AD_Start = AD_Start + u_AD%HubMotion%NNodes * 3 ! move past the AD translation disp field to orientation field         
   call SetBlockMatrix( dUdy, MeshMapData%ED_P_2_AD_P_H%dM%mi, AD_Start, ED_Out_Start )
      
   ! *** AD rotational velocity: from ED rotational velocity
   AD_Start = AD_Start + u_AD%HubMotion%NNodes * 3 ! move past the AD orientation field to rotational velocity field          
   ED_Out_Start = Indx_y_ED_Hub_Start(y_ED, y_FAST) + y_ED%HubPtMotion%NNodes * 6 ! ! start of y_ED%HubPtMotion%RotationVel field
   call SetBlockMatrix( dUdy, MeshMapData%ED_P_2_AD_P_H%dM%mi, AD_Start, ED_Out_Start )
            

   
      !...................................
      ! blade root   
      !...................................
   DO k=1,size(y_ED%BladeRootMotion)
      CALL Linearize_Point_to_Point( y_ED%BladeRootMotion(k), u_AD%BladeRootMotion(k), MeshMapData%ED_P_2_AD_P_R(k), ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//':u_AD%BladeRootMotion('//trim(num2lstr(k))//')' )      
         if (errStat>=AbortErrLev) return
               
      ! *** AD orientation: from ED orientation
      AD_Start = Indx_u_AD_BladeRoot_Start(u_AD, y_FAST, k)       ! start of u_AD%BladeRootMotion(k)%Orientation field
      
      ED_Out_Start = Indx_y_ED_BladeRoot_Start(y_ED, y_FAST, k) & ! start of y_ED%BladeRootMotion(k)%TranslationDisp field
                   + y_ED%BladeRootMotion(k)%NNodes * 3           ! start of y_ED%BladeRootMotion(k)%Orientation field
      call SetBlockMatrix( dUdy, MeshMapData%ED_P_2_AD_P_R(k)%dM%mi, AD_Start, ED_Out_Start )
                  
   END DO
      
   
      !...................................
      ! blades
      !...................................
   IF (p_FAST%CompElast == Module_ED ) THEN
            
      
      DO k=1,size(y_ED%BladeLn2Mesh)
         !!! ! This linearization was done in forming dUdu (see Linear_AD_InputSolve_du()), so we don't need to re-calculate these matrices 
         !!! ! while forming dUdy, too.
         !!!CALL Linearize_Line2_to_Line2( y_ED%BladeLn2Mesh(k), u_AD%BladeMotion(k), MeshMapData%BDED_L_2_AD_L_B(k), ErrStat2, ErrMsg2 )
         
         AD_Start = Indx_u_AD_Blade_Start(u_AD, y_FAST, k)     ! start of u_AD%BladeMotion(k)%TranslationDisp field
         ED_Out_Start = Indx_y_ED_Blade_Start(y_ED, y_FAST, k) ! start of y_ED%BladeLn2Mesh(k)%TranslationDisp field
         CALL Assemble_dUdy_Motions(y_ED%BladeLn2Mesh(k), u_AD%BladeMotion(k), MeshMapData%BDED_L_2_AD_L_B(k), AD_Start, ED_Out_Start, dUdy, .true.)
         
      END DO
      
   ELSEIF (p_FAST%CompElast == Module_BD ) THEN
      !!! ! This linearization was done in forming dUdu (see Linear_AD_InputSolve_du()), so we don't need to re-calculate these matrices 
      !!! ! while forming dUdy, too.
      !!!CALL Linearize_Line2_to_Line2( BD%y(k)%BldMotion, u_AD%BladeMotion(k), MeshMapData%BDED_L_2_AD_L_B(k), ErrStat2, ErrMsg2 )
      
      DO k=1,p_FAST%nBeams
         AD_Start     = Indx_u_AD_Blade_Start(u_AD, y_FAST, k)     ! start of u_AD%BladeMotion(k)%TranslationDisp field
         BD_Out_Start = y_FAST%Lin%Modules(Module_BD)%Instance(k)%LinStartIndx(LIN_OUTPUT_COL)
         
         CALL Assemble_dUdy_Motions(BD%y(k)%BldMotion, u_AD%BladeMotion(k), MeshMapData%BDED_L_2_AD_L_B(k), AD_Start, BD_Out_Start, dUdy, .true.)
      END DO
      
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
   INTEGER(IntKi)                          :: k                   ! module instance loop counter
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
      call AllocAry(y_FAST%Lin%Glue%A, y_FAST%Lin%Glue%SizeLin(LIN_ContSTATE_COL), &
                                       y_FAST%Lin%Glue%SizeLin(LIN_ContSTATE_COL), 'A', ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         
      !B: rows = x; columns = u
      call AllocAry(y_FAST%Lin%Glue%B, y_FAST%Lin%Glue%SizeLin(LIN_ContSTATE_COL), &
                                       y_FAST%Lin%Glue%SizeLin(LIN_INPUT_COL), 'B', ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         
      !C: rows = y; columns = x
      call AllocAry(y_FAST%Lin%Glue%C, y_FAST%Lin%Glue%SizeLin(LIN_OUTPUT_COL), &
                                       y_FAST%Lin%Glue%SizeLin(LIN_ContSTATE_COL), 'C', ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         
      !D: rows = y; columns = u
      call AllocAry(y_FAST%Lin%Glue%D, y_FAST%Lin%Glue%SizeLin(LIN_OUTPUT_COL), &
                                       y_FAST%Lin%Glue%SizeLin(LIN_INPUT_COL), 'D', ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
                  
      if (ErrStat>=AbortErrLev) return
   end if
   
   
   ! The equations of the matrices returned from this routine are really just a general form with the null matrices removed:      
   
   ! A
   y_FAST%Lin%Glue%A = 0.0_R8Ki
   r_start = 1
   c_start = 1
   do i = 1,p_FAST%Lin_NumMods
      ThisModule = p_FAST%Lin_ModOrder( i )
      
      do k=1,size(y_FAST%Lin%Modules(ThisModule)%Instance)
         if (allocated( y_FAST%Lin%Modules(ThisModule)%Instance(k)%A) ) then
            do c=1,size( y_FAST%Lin%Modules(ThisModule)%Instance(k)%A, 2)
               do r=1,size( y_FAST%Lin%Modules(ThisModule)%Instance(k)%A, 1)
                  y_FAST%Lin%Glue%A(r_start + r - 1, c_start + c - 1) = y_FAST%Lin%Modules(ThisModule)%Instance(k)%A(r,c)
               end do
            end do
         end if
      
         r_start = r_start + y_FAST%Lin%Modules(ThisModule)%Instance(k)%SizeLin(LIN_ContSTATE_COL)
         c_start = c_start + y_FAST%Lin%Modules(ThisModule)%Instance(k)%SizeLin(LIN_ContSTATE_COL)
      end do
   end do
   
    
   ! B
   y_FAST%Lin%Glue%B = 0.0_R8Ki
   r_start = 1
   c_start = 1
   do i = 1,p_FAST%Lin_NumMods
      ThisModule = p_FAST%Lin_ModOrder( i )
      
      do k=1,size(y_FAST%Lin%Modules(ThisModule)%Instance)
         if (allocated( y_FAST%Lin%Modules(ThisModule)%Instance(k)%B) ) then
            do c=1,size( y_FAST%Lin%Modules(ThisModule)%Instance(k)%B, 2)
               do r=1,size( y_FAST%Lin%Modules(ThisModule)%Instance(k)%B, 1)
                  y_FAST%Lin%Glue%B(r_start + r - 1, c_start + c - 1) = y_FAST%Lin%Modules(ThisModule)%Instance(k)%B(r,c)
               end do
            end do
         end if
      
         r_start = r_start + y_FAST%Lin%Modules(ThisModule)%Instance(k)%SizeLin(LIN_ContSTATE_COL)
         c_start = c_start + y_FAST%Lin%Modules(ThisModule)%Instance(k)%SizeLin(LIN_INPUT_COL)
      end do
   end do
   
   ! C
   y_FAST%Lin%Glue%C = 0.0_R8Ki
   r_start = 1
   c_start = 1
   do i = 1,p_FAST%Lin_NumMods
      ThisModule = p_FAST%Lin_ModOrder( i )
      
      do k=1,size(y_FAST%Lin%Modules(ThisModule)%Instance)
         if (allocated( y_FAST%Lin%Modules(ThisModule)%Instance(k)%C) ) then
            do c=1,size( y_FAST%Lin%Modules(ThisModule)%Instance(k)%C, 2)
               do r=1,size( y_FAST%Lin%Modules(ThisModule)%Instance(k)%C, 1)
                  y_FAST%Lin%Glue%C(r_start + r - 1, c_start + c - 1) = y_FAST%Lin%Modules(ThisModule)%Instance(k)%C(r,c)
               end do
            end do
         end if
      
         r_start = r_start + y_FAST%Lin%Modules(ThisModule)%Instance(k)%SizeLin(LIN_OUTPUT_COL)
         c_start = c_start + y_FAST%Lin%Modules(ThisModule)%Instance(k)%SizeLin(LIN_ContSTATE_COL)
      end do
   end do   
   
   ! D
   y_FAST%Lin%Glue%D = 0.0_R8Ki
   r_start = 1
   c_start = 1
   do i = 1,p_FAST%Lin_NumMods
      ThisModule = p_FAST%Lin_ModOrder( i )
      
      do k=1,size(y_FAST%Lin%Modules(ThisModule)%Instance)
         if (allocated( y_FAST%Lin%Modules(ThisModule)%Instance(k)%D) ) then
            do c=1,size( y_FAST%Lin%Modules(ThisModule)%Instance(k)%D, 2)
               do r=1,size( y_FAST%Lin%Modules(ThisModule)%Instance(k)%D, 1)
                  y_FAST%Lin%Glue%D(r_start + r - 1, c_start + c - 1) = y_FAST%Lin%Modules(ThisModule)%Instance(k)%D(r,c)
               end do
            end do
         end if
      
         r_start = r_start + y_FAST%Lin%Modules(ThisModule)%Instance(k)%SizeLin(LIN_OUTPUT_COL)
         c_start = c_start + y_FAST%Lin%Modules(ThisModule)%Instance(k)%SizeLin(LIN_INPUT_COL)
      end do
   end do      
   
   
END SUBROUTINE Glue_FormDiag      
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine forms the full-system state matrices for linearization: A, B, C, and D.
!! Note that it uses LAPACK_GEMM instead of MATMUL for matrix multiplications because of stack-space issues (these
!! matrices get large quickly).
SUBROUTINE Glue_StateMatrices( p_FAST, y_FAST, dUdu, dUdy, ErrStat, ErrMsg )

   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST              !< Output variables for the glue code
   REAL(R8Ki),               INTENT(INOUT) :: dUdu(:,:)           !< glue-code Jacobian: \f$ \frac{\partial U}{\partial u} \f$; on exit will hold G^{-1}*dUdu
   REAL(R8Ki),               INTENT(INOUT) :: dUdy(:,:)           !< glue-code Jacobian: \f$ \frac{\partial U}{\partial y} \f$; on exit will hold G^{-1}*dUdy
           
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None

   
   
      ! local variables
   REAL(R8Ki), ALLOCATABLE                 :: G(:,:), tmp(:,:) ! variables for glue-code linearization
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
   call LAPACK_GEMM( 'N', 'N', 1.0_R8Ki, dUdy, y_FAST%Lin%Glue%D, 1.0_R8Ki, G, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
   ! because G can be ill-conditioned, we are going to precondition with G_hat = S^(-1) * G * S
   ! we will also multiply the right-hand-side of the equations that need G inverse so that 
   ! dUdy_hat = S^(-1)*dUdy and dUdu_hat = S^(-1)*dUdu
   call Precondition(p_FAST, y_FAST, G, dUdu, dUdy)

   
   ! now we need to form G_hat^(-1) * (S^-1*dUdy) and G^(-1) * (S^-1*dUdu)   
      ! factor G for the two solves:
   CALL LAPACK_getrf( M=size(G,1), N=size(G,2), A=G, IPIV=ipiv, ErrStat=ErrStat2, ErrMsg=ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (ErrStat >= AbortErrLev) then
         call cleanup() 
         return
      end if
      
    ! after the this solve, dUdy holds G_hat^(-1) * dUdy_hat:
   CALL LAPACK_getrs( trans='N', N=size(G,2), A=G, IPIV=ipiv, B=dUdy, ErrStat=ErrStat2, ErrMsg=ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

    ! after the this solve, dUdu holds G_hat^(-1) * dUdu_hat:
   CALL LAPACK_getrs( trans='N', N=size(G,2), A=G, IPIV=ipiv, B=dUdu, ErrStat=ErrStat2, ErrMsg=ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)   
         
   deallocate(G)    ! we're finished with the solves, so let's get rid of them
   deallocate(ipiv) ! we're finished with the solves, so let's get rid of them

   ! after this call, dUdu holds G^(-1)*dUdu and dUdy holds G^(-1)*dUdy:
   call Postcondition(p_FAST, y_FAST, dUdu, dUdy)

                    
   ! *** get tmp matrix  for A and C calculations ****
   !----------------------         
   call AllocAry(tmp, y_FAST%Lin%Glue%SizeLin(LIN_INPUT_COL), y_FAST%Lin%Glue%SizeLin(LIN_ContSTATE_COL), 'G^-1*dUdy*C', ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (errStat>=AbortErrLev) then
         call cleanup()
         return
      end if
      
   !tmp = G^(-1) * dUdy * diag(C)
   call LAPACK_GEMM( 'N', 'N', 1.0_R8Ki, dUdy, y_FAST%Lin%Glue%C, 0.0_R8Ki, tmp, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
                        
   !  A
   !----------------------         
   !> \f{equation}{ A = 
   !! \begin{bmatrix} A^{ED} & 0 \\ 0 & A^{BD} \end{bmatrix} - 
   !! \begin{bmatrix} 0 & 0 & B^{ED} & 0 & 0 \\ 0 & 0 & 0 & B^{BD} & 0\end{bmatrix} \,
   !! \begin{bmatrix} G \end{bmatrix}^{-1} \, \frac{\partial U}{\partial y} \, \begin{bmatrix} 0 & 0 \\ 0 & 0 \\ C^{ED} & 0 \\ 0 & C^{BD} \\ 0 & 0 \end{bmatrix}
   !! \f}
   !y_FAST%Lin%Glue%A = y_FAST%Lin%Glue%A - matmul( y_FAST%Lin%Glue%B, tmp )  
   call LAPACK_GEMM( 'N', 'N', -1.0_R8Ki, y_FAST%Lin%Glue%B, tmp, 1.0_R8Ki, y_FAST%Lin%Glue%A, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

   !  C
   !----------------------         
   !> \f{equation}{ C = \begin{bmatrix} 0 & 0 \\ 0 & 0 \\ C^{ED} & 0 \\ 0 & C^{BD} \\ 0 & 0 \end{bmatrix} - 
   !! \begin{bmatrix} D^{IfW} & 0 & 0 & 0 & 0 \\ 0 &  D^{SrvD} & 0 & 0 & 0 \\ 0 & 0 &  D^{ED} & 0 & 0 \\ 0 & 0 & 0 & D^{BD} & 0\\ 0 & 0 & 0 & 0 & D^{AD}\end{bmatrix} \,
   !! \begin{bmatrix} G \end{bmatrix}^{-1} \, \frac{\partial U}{\partial y} \, \begin{bmatrix} 0 & 0 \\ 0 & 0 \\ C^{ED} & 0 \\ 0 & C^{BD} \\ 0 & 0 \end{bmatrix}
   !! \f}
   !y_FAST%Lin%Glue%C = y_FAST%Lin%Glue%C - matmul( y_FAST%Lin%Glue%D, tmp ) 
   call LAPACK_GEMM( 'N', 'N', -1.0_R8Ki, y_FAST%Lin%Glue%D, tmp, 1.0_R8Ki, y_FAST%Lin%Glue%C, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
   deallocate(tmp)
   
         
   !  B
   !----------------------         
   !> \f{equation}{ B = \begin{bmatrix} 0 & 0 \\ 0 & 0 \\ B^{ED} & 0 \\ 0 & B^{BD} \\ 0 & 0 \end{bmatrix} \,
   !! \begin{bmatrix} G \end{bmatrix}^{-1} \, \frac{\partial U}{\partial u}
   !! \f}
   call AllocAry(tmp,size(y_FAST%Lin%Glue%B,1),size(y_FAST%Lin%Glue%B,2),'tmp',ErrStat2,ErrMsg2)
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (errStat>=AbortErrLev) then
         call cleanup()
         return
      end if
   tmp = y_FAST%Lin%Glue%B   
         
   !y_FAST%Lin%Glue%B = matmul( y_FAST%Lin%Glue%B, dUdu ) 
   call LAPACK_GEMM( 'N', 'N', 1.0_R8Ki, tmp, dUdu, 0.0_R8Ki, y_FAST%Lin%Glue%B, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   deallocate(tmp)      
      
   !  D
   !----------------------         
   !> \f{equation}{ D = \begin{bmatrix} D^{IfW} & 0 & 0 & 0 & 0 \\ 0 &  D^{SrvD} & 0 & 0 & 0 \\ 0 & 0 &  D^{ED} & 0 & 0 \\ 0 & 0 & 0 & D^{BD} & 0\\ 0 & 0 & 0 & 0 & D^{AD}\end{bmatrix} \,
   !! \begin{bmatrix} G \end{bmatrix}^{-1} \, \frac{\partial U}{\partial u}
   !! \f}
   call AllocAry(tmp,size(y_FAST%Lin%Glue%D,1),size(y_FAST%Lin%Glue%D,2),'tmp',ErrStat2,ErrMsg2)
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (errStat>=AbortErrLev) then
         call cleanup()
         return
      end if
   tmp = y_FAST%Lin%Glue%D
         
   !y_FAST%Lin%Glue%D = matmul( y_FAST%Lin%Glue%D, dUdu )
   call LAPACK_GEMM( 'N', 'N', 1.0_R8Ki, tmp, dUdu, 0.0_R8Ki, y_FAST%Lin%Glue%D, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)   
   deallocate(tmp)    
   
   call cleanup()
      
contains
   subroutine cleanup()
      if (allocated(ipiv)) deallocate(ipiv)     
      if (allocated(G)) deallocate(G)
      if (allocated(tmp)) deallocate(tmp)
   end subroutine cleanup   
END SUBROUTINE Glue_StateMatrices
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine returns the preconditioned matrix, \f$ \hat{G} \f$, such that \f$ \hat{G} = S^(-1) G S \f$ with \f$S^(-1)\f$ defined
!! such that loads are scaled by p_FAST\%UJacSclFact. It also returns the preconditioned matrices \f$ \hat{dUdu} \f$ and 
!! \f$ \hat{dUdy} \f$ such that \f$ \hat{dUdu} = S^(-1) dUdu \f$ and
!! \f$ \hat{dUdy} = S^(-1) dUdy \f$ for the right-hand sides of the equations to be solved.
SUBROUTINE Precondition(p_FAST, y_FAST, G, dUdu, dUdy)


   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST              !< Output variables for the glue code
   REAL(R8Ki),               INTENT(INOUT) :: G(:,:)              !< variable for glue-code linearization (in is G; out is G_hat)
   REAL(R8Ki),               INTENT(INOUT) :: dUdu(:,:)           !< jacobian in FAST linearization from right-hand-side of equation
   REAL(R8Ki),               INTENT(INOUT) :: dUdy(:,:)           !< jacobian in FAST linearization from right-hand-side of equation

   integer :: r, c
   
   !! Change G to G_hat:
   do c = 1,size(y_FAST%Lin%Glue%IsLoad_u)
   
      if ( y_FAST%Lin%Glue%IsLoad_u(c) ) then

         do r = 1,size(y_FAST%Lin%Glue%IsLoad_u)
            if ( .not. y_FAST%Lin%Glue%IsLoad_u(r) ) then
               ! column is load, but row is a motion:
               G(r,c) = G(r,c) * p_FAST%UJacSclFact
            end if
         end do
         
      else
      
         do r = 1,size(y_FAST%Lin%Glue%IsLoad_u)
            if ( y_FAST%Lin%Glue%IsLoad_u(r) ) then
               ! column is motion, but row is a load:
               G(r,c) = G(r,c) / p_FAST%UJacSclFact
            end if
         end do
         
      end if
         
   end do
   
   
   !! Change dUdu to dUdu_hat (note that multiplying on the left multiplies the entire row):
   do r = 1,size(y_FAST%Lin%Glue%IsLoad_u)
   
      if ( y_FAST%Lin%Glue%IsLoad_u(r) ) then
         dUdu(r,:) = dUdu(r,:) / p_FAST%UJacSclFact
      end if
   
   end do

   !! Change dUdy to dUdy_hat:
   do r = 1,size(y_FAST%Lin%Glue%IsLoad_u)
   
      if ( y_FAST%Lin%Glue%IsLoad_u(r) ) then
         dUdy(r,:) = dUdy(r,:) / p_FAST%UJacSclFact
      end if
   
   end do
   
   
END SUBROUTINE Precondition
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine returns the matrices \f$ \tilde{dUdu} \f$ and \f$ \tilde{dUdy} \f$ such that 
!! \f$ \tilde{dUdu} = G^(-1) dUdu \f$ and
!! \f$ \tilde{dUdy} = G^(-1) dUdy \f$, which have been solved using the preconditioned system defined in fast_lin::precondition.
SUBROUTINE Postcondition(p_FAST, y_FAST, dUdu, dUdy)


   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST              !< Output variables for the glue code
   REAL(R8Ki),               INTENT(INOUT) :: dUdu(:,:)           !< jacobian in FAST linearization from right-hand-side of equation
   REAL(R8Ki),               INTENT(INOUT) :: dUdy(:,:)           !< jacobian in FAST linearization from right-hand-side of equation

   integer :: r
   
   !! Change S^(-1) * G_hat^(-1) * dUdu_hat to G^(-1) * dUdu (note that multiplying on the left multiplies the entire row):
   do r = 1,size(y_FAST%Lin%Glue%IsLoad_u)
   
      if ( y_FAST%Lin%Glue%IsLoad_u(r) ) then
         dUdu(r,:) = dUdu(r,:) * p_FAST%UJacSclFact
      end if
   
   end do

   !! Change S^(-1) * G_hat^(-1) * dUdy_hat to G^(-1) * dUdy (note that multiplying on the left multiplies the entire row):
   do r = 1,size(y_FAST%Lin%Glue%IsLoad_u)
   
      if ( y_FAST%Lin%Glue%IsLoad_u(r) ) then
         dUdy(r,:) = dUdy(r,:) * p_FAST%UJacSclFact
      end if
   
   end do
   
   
END SUBROUTINE Postcondition
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SetBlockMatrix( matrix, submatrix, RowStart, ColStart )
   REAL(R8Ki),     INTENT(INOUT)  :: matrix(:,:)      !< matrix that will have the negative of the submatrix block added to it
   REAL(R8Ki),     INTENT(IN )    :: submatrix(:,:)   !< block matrix that needs to be added to matrix
   INTEGER(IntKi), INTENT(IN )    :: RowStart         !< first row in matrix where submatrix should start
   INTEGER(IntKi), INTENT(IN )    :: ColStart         !< first column in matrix where submatrix should start

   INTEGER(IntKi)                 :: col
   INTEGER(IntKi)                 :: row

   
   do col=1,size( submatrix, 2)
      do row=1,size( submatrix, 1)
         matrix(RowStart + row - 1, ColStart + col - 1) = - submatrix(row,col) ! note the negative sign here!!!!
      end do
   end do
   
   
END SUBROUTINE SetBlockMatrix
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SumBlockMatrix( matrix, submatrix, RowStart, ColStart )
   REAL(R8Ki),     INTENT(INOUT)  :: matrix(:,:)      !< matrix that will have the negative of the submatrix block added to it
   REAL(R8Ki),     INTENT(IN )    :: submatrix(:,:)   !< block matrix that needs to be added to matrix
   INTEGER(IntKi), INTENT(IN )    :: RowStart         !< first row in matrix where submatrix should start
   INTEGER(IntKi), INTENT(IN )    :: ColStart         !< first column in matrix where submatrix should start

   INTEGER(IntKi)                 :: col
   INTEGER(IntKi)                 :: row

   
   do col=1,size( submatrix, 2)
      do row=1,size( submatrix, 1)
         matrix(RowStart + row - 1, ColStart + col - 1) = matrix(RowStart + row - 1, ColStart + col - 1) &
                                                        - submatrix(row,col) ! note the negative sign here!!!!
      end do
   end do
   
   
END SUBROUTINE SumBlockMatrix
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine assembles the linearization matrices for transfer of motion fields between two meshes.
!> It set the following block matrix, which is the dUdy block for transfering output (source) mesh \f$y\f$ to the
!! input (destination) mesh \f$u\f$:\n
!! \f$ M = - \begin{bmatrix} M_{mi}      & M_{f_{\times p}} & 0      & 0                & 0      & 0                \\
!!                           0           & M_{mi}           & 0      & 0                & 0      & 0                \\
!!                           M_{tv\_uS}  & 0                & M_{mi} & M_{f_{\times p}} & 0      & 0                \\
!!                           0           & 0                & 0      & M_{mi}           & 0      & 0                \\
!!                           M_{ta\_uS}  & 0                & 0      & M_{ta\_rv}       & M_{mi} & M_{f_{\times p}} \\
!!                           0           & 0                & 0      & 0                & 0      & M_{mi}           \\
!! \end{bmatrix} \f$
!! where the matrices correspond to 
!! \f$ \left\{ \begin{matrix}
!!      \vec{u}^S \\
!!      \vec{\theta}^S \\
!!      \vec{v}^S \\
!!      \vec{\omega}^S \\
!!      \vec{a}^S \\
!!      \vec{\alpha}^S \\
!! \end{matrix} \right\} \f$
SUBROUTINE Assemble_dUdy_Motions(y, u, MeshMap, BlockRowStart, BlockColStart, dUdy, skipRotVel)
   TYPE(MeshType),    INTENT(IN)     :: y             !< the output (source) mesh that is transfering motions
   TYPE(MeshType),    INTENT(IN)     :: u             !< the input (destination) mesh that is receiving motions
   TYPE(MeshMapType), INTENT(IN)     :: MeshMap       !< the mesh mapping from y to u
   INTEGER(IntKi),    INTENT(IN)     :: BlockRowStart !< the index of the row defining the block of dUdy to be set
   INTEGER(IntKi),    INTENT(IN)     :: BlockColStart !< the index of the column defining the block of dUdy to be set
   REAL(R8Ki),        INTENT(INOUT)  :: dUdy(:,:)     !< full Jacobian matrix
   LOGICAL, OPTIONAL, INTENT(IN)     :: skipRotVel    !< if present and true, we skip the rotational velocity and acceleration fields and return early
   
   INTEGER(IntKi)                    :: row
   INTEGER(IntKi)                    :: col
   
!! \f$M_{mi}\f$ is modmesh_mapping::meshmaplinearizationtype::mi (motion identity)\n
!! \f$M_{f_{\times p}}\f$ is modmesh_mapping::meshmaplinearizationtype::fx_p \n
!! \f$M_{tv\_uD}\f$ is modmesh_mapping::meshmaplinearizationtype::tv_uD \n
!! \f$M_{tv\_uS}\f$ is modmesh_mapping::meshmaplinearizationtype::tv_uS \n
!! \f$M_{ta\_uD}\f$ is modmesh_mapping::meshmaplinearizationtype::ta_uD \n
!! \f$M_{ta\_uS}\f$ is modmesh_mapping::meshmaplinearizationtype::ta_uS \n
!! \f$M_{ta\_rv}\f$ is modmesh_mapping::meshmaplinearizationtype::ta_rv \n

      !*** row for translational displacement ***
         ! source translational displacement to destination translational displacement:
      row = BlockRowStart                    ! start of u%TranslationDisp field
      col = BlockColStart                    ! start of y%TranslationDisp field
      call SetBlockMatrix( dUdy, MeshMap%dM%mi, row, col )

         ! source orientation to destination translational displacement:
      row = BlockRowStart                    ! start of u%TranslationDisp field
      col = BlockColStart + y%NNodes*3       ! start of y%Orientation field [skip 1 field with 3 components]
      call SetBlockMatrix( dUdy, MeshMap%dM%fx_p, row, col )


      !*** row for orientation ***
         ! source orientation to destination orientation:
      row = BlockRowStart + u%NNodes*3       ! start of u%Orientation field [skip 1 field with 3 components]
      col = BlockColStart + y%NNodes*3       ! start of y%Orientation field [skip 1 field with 3 components]
      call SetBlockMatrix( dUdy, MeshMap%dM%mi, row, col )


      !*** row for translational velocity ***
         ! source translational displacement to destination translational velocity:
      row = BlockRowStart + u%NNodes*6       ! start of u%TranslationVel field [skip 2 fields with 3 components]
      col = BlockColStart                    ! start of y%TranslationDisp field
      call SetBlockMatrix( dUdy, MeshMap%dM%tv_us, row, col )

         ! source translational velocity to destination translational velocity:
      row = BlockRowStart + u%NNodes*6       ! start of u%TranslationVel field [skip 2 fields with 3 components]
      col = BlockColStart + y%NNodes*6       ! start of y%TranslationVel field [skip 2 fields with 3 components]
      call SetBlockMatrix( dUdy, MeshMap%dM%mi, row, col )

         ! source rotational velocity to destination translational velocity:
      row = BlockRowStart + u%NNodes*6       ! start of u%TranslationVel field [skip 2 fields with 3 components]
      col = BlockColStart + y%NNodes*9       ! start of y%RotationVel field [skip 3 fields with 3 components]
      call SetBlockMatrix( dUdy, MeshMap%dM%fx_p, row, col )


      if (PRESENT(skipRotVel)) then
         if (skipRotVel) return ! destination does not include rotational velocities or accelerations, so we'll just return
      end if


      !*** row for rotational velocity ***
         ! source rotational velocity to destination rotational velocity:
      row = BlockRowStart + u%NNodes*9       ! start of u%RotationVel field [skip 3 fields with 3 components]
      col = BlockColStart + y%NNodes*9       ! start of y%RotationVel field [skip 3 fields with 3 components]
      call SetBlockMatrix( dUdy, MeshMap%dM%mi, row, col )


      !*** row for translational acceleration ***
         ! source translational displacement to destination translational acceleration:
      row = BlockRowStart + u%NNodes*12      ! start of u%TranslationAcc field [skip 4 fields with 3 components]
      col = BlockColStart                    ! start of y%TranslationDisp field
      call SetBlockMatrix( dUdy, MeshMap%dM%ta_us, row, col )

         ! source rotational velocity to destination translational acceleration:
      row = BlockRowStart + u%NNodes*12      ! start of u%TranslationAcc field [skip 4 fields with 3 components]
      col = BlockColStart + y%NNodes*9       ! start of y%RotationVel field [skip 3 fields with 3 components]
      call SetBlockMatrix( dUdy, MeshMap%dM%ta_rv, row, col )

         ! source translational acceleration to destination translational acceleration:
      row = BlockRowStart + u%NNodes*12      ! start of u%TranslationAcc field [skip 4 fields with 3 components]
      col = BlockColStart + y%NNodes*12      ! start of y%TranslationAcc field [skip 4 fields with 3 components]
      call SetBlockMatrix( dUdy, MeshMap%dM%mi, row, col )

         ! source rotational acceleration to destination translational acceleration:
      row = BlockRowStart + u%NNodes*12      ! start of u%TranslationAcc field [skip 4 fields with 3 components]
      col = BlockColStart + y%NNodes*15      ! start of y%RotationAcc field [skip 5 fields with 3 components]
      call SetBlockMatrix( dUdy, MeshMap%dM%fx_p, row, col )


      !*** row for rotational acceleration ***
         ! source rotational acceleration to destination rotational acceleration
      row = BlockRowStart + u%NNodes*15      ! start of u%RotationAcc field [skip 5 fields with 3 components]
      col = BlockColStart + y%NNodes*15      ! start of y%RotationAcc field [skip 5 fields with 3 components]
      call SetBlockMatrix( dUdy, MeshMap%dM%mi, row, col )


END SUBROUTINE Assemble_dUdy_Motions
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine assembles the linearization matrices for transfer of load fields between two meshes.
!> It set the following block matrix, which is the dUdy block for transfering output (source) mesh \f$y\f$ to the
!! input (destination) mesh \f$u\f$:\n
!! \f$ M = - \begin{bmatrix} M_{li}      & 0       \\
!!                           M_{fm}      & M_{li}  \\
!! \end{bmatrix} \f$
!! & M_{mi} & } 
!! \f$ \left\{ \begin{matrix}
!!      \vec{F}^S \\
!!      \vec{M}^S
!! \end{matrix} \right\} \f$
SUBROUTINE Assemble_dUdy_Loads(y, u, MeshMap, BlockRowStart, BlockColStart, dUdy)
   TYPE(MeshType),    INTENT(IN)     :: y             !< the output (source) mesh that is transfering loads
   TYPE(MeshType),    INTENT(IN)     :: u             !< the input (destination) mesh that is receiving loads
   TYPE(MeshMapType), INTENT(IN)     :: MeshMap       !< the mesh mapping from y to u
   INTEGER(IntKi),    INTENT(IN)     :: BlockRowStart !< the index of the row defining the block of dUdy to be set
   INTEGER(IntKi),    INTENT(IN)     :: BlockColStart !< the index of the column defining the block of dUdy to be set
   REAL(R8Ki),        INTENT(INOUT)  :: dUdy(:,:)     !< full Jacobian matrix
   
   INTEGER(IntKi)                    :: row
   INTEGER(IntKi)                    :: col
   
      !*** row for force ***
         ! source force to destination force:
      row = BlockRowStart                    ! start of u%Force field
      col = BlockColStart                    ! start of y%Force field
      call SetBlockMatrix( dUdy, MeshMap%dM%li, row, col )

      !*** row for moment ***
         ! source force to destination moment:
      row = BlockRowStart + u%NNodes*3       ! start of u%Moment field [skip 1 field with 3 components]
      col = BlockColStart                    ! start of y%Force field
      call SetBlockMatrix( dUdy, MeshMap%dM%m_f, row, col )

         ! source moment to moment:
      row = BlockRowStart + u%NNodes*3       ! start of u%Moment field [skip 1 field with 3 components]
      col = BlockColStart + y%NNodes*3       ! start of y%Moment field [skip 1 field with 3 components]
      call SetBlockMatrix( dUdy, MeshMap%dM%li, row, col )

END SUBROUTINE Assemble_dUdy_Loads


!----------------------------------------------------------------------------------------------------------------------------------
!> This routine returns the starting index for the u_ED%BladePtLoads(BladeNum) mesh in the FAST linearization inputs.
FUNCTION Indx_u_ED_Blade_Start(u_ED, y_FAST, BladeNum) RESULT(ED_Start)
   TYPE(FAST_OutputFileType),      INTENT(IN )  :: y_FAST           !< FAST output file data (for linearization)
   TYPE(ED_InputType),             INTENT(IN )  :: u_ED             !< ED Inputs at t
   INTEGER,                        INTENT(IN )  :: BladeNum         !< blade number to find index for
   INTEGER                                      :: k                !< blade number loop

   INTEGER                                      :: ED_Start         !< starting index of this blade mesh in ElastoDyn inputs

   ED_Start = y_FAST%Lin%Modules(Module_ED)%Instance(1)%LinStartIndx(LIN_INPUT_COL) 
   if (allocated(u_ED%BladePtLoads)) then
      do k = 1,min(BladeNum-1, size(u_ED%BladePtLoads))
         ED_Start = ED_Start + u_ED%BladePtLoads(k)%NNodes * 6  ! 3 forces + 3 moments at each node on each blade
      end do
   end if

END FUNCTION Indx_u_ED_Blade_Start
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine returns the starting index for the u_ED%PlatformPtMesh mesh in the FAST linearization inputs.
FUNCTION Indx_u_ED_Platform_Start(u_ED, y_FAST) RESULT(ED_Start)
   TYPE(FAST_OutputFileType),      INTENT(IN )  :: y_FAST           !< FAST output file data (for linearization)
   TYPE(ED_InputType),             INTENT(IN )  :: u_ED             !< ED Inputs at t
   INTEGER                                      :: ED_Start         !< starting index of this mesh

   ED_Start = Indx_u_ED_Blade_Start(u_ED, y_FAST, MaxNBlades+1) ! skip all of the blades to get to start of platform
END FUNCTION Indx_u_ED_Platform_Start
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine returns the starting index for the u_ED%TowerPtLoads mesh in the FAST linearization inputs.
FUNCTION Indx_u_ED_Tower_Start(u_ED, y_FAST) RESULT(ED_Start)
   TYPE(FAST_OutputFileType),      INTENT(IN )  :: y_FAST           !< FAST output file data (for linearization)
   TYPE(ED_InputType),             INTENT(IN )  :: u_ED             !< ED Inputs at t

   INTEGER                                      :: ED_Start         !< starting index of this mesh

   ED_Start = Indx_u_ED_Platform_Start(u_ED, y_FAST)
   ED_Start = ED_Start + u_ED%PlatformPtMesh%NNodes * 6            ! 3 forces + 3 moments at each node
END FUNCTION Indx_u_ED_Tower_Start
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine returns the starting index for the u_ED%HubPtLoad mesh in the FAST linearization inputs.
FUNCTION Indx_u_ED_Hub_Start(u_ED, y_FAST) RESULT(ED_Start)
   TYPE(FAST_OutputFileType),      INTENT(IN )  :: y_FAST           !< FAST output file data (for linearization)
   TYPE(ED_InputType),             INTENT(IN )  :: u_ED             !< ED Inputs at t

   INTEGER                                      :: ED_Start         !< starting index of this mesh

   ED_Start = Indx_u_ED_Tower_Start(u_ED, y_FAST)
   ED_Start = ED_Start + u_ED%TowerPtLoads%NNodes * 6            ! 3 forces + 3 moments at each node
END FUNCTION Indx_u_ED_Hub_Start
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine returns the starting index for the u_ED%NacelleLoads mesh in the FAST linearization inputs.
FUNCTION Indx_u_ED_Nacelle_Start(u_ED, y_FAST) RESULT(ED_Start)
   TYPE(FAST_OutputFileType),      INTENT(IN )  :: y_FAST           !< FAST output file data (for linearization)
   TYPE(ED_InputType),             INTENT(IN )  :: u_ED             !< ED Inputs at t

   INTEGER                                      :: ED_Start         !< starting index of this mesh

   ED_Start = Indx_u_ED_Hub_Start(u_ED, y_FAST)
   ED_Start = ED_Start + u_ED%HubPtLoad%NNodes * 6            ! 3 forces + 3 moments at each node
END FUNCTION Indx_u_ED_Nacelle_Start
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine returns the starting index for the u_ED%BladePitchCom array in the FAST linearization inputs.
FUNCTION Indx_u_ED_BlPitchCom_Start(u_ED, y_FAST) RESULT(ED_Start)
   TYPE(FAST_OutputFileType),      INTENT(IN )  :: y_FAST           !< FAST output file data (for linearization)
   TYPE(ED_InputType),             INTENT(IN )  :: u_ED             !< ED Inputs at t

   INTEGER                                      :: ED_Start         !< starting index of this mesh

   ED_Start = Indx_u_ED_Nacelle_Start(u_ED, y_FAST)
   ED_Start = ED_Start + u_ED%NacelleLoads%NNodes * 6            ! 3 forces + 3 moments at each node
END FUNCTION Indx_u_ED_BlPitchCom_Start
!----------------------------------------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------------------------------------
!> This routine returns the starting index for the y_ED%BladeLn2Mesh(BladeNum) mesh in the FAST linearization outputs.
FUNCTION Indx_y_ED_Blade_Start(y_ED, y_FAST, BladeNum) RESULT(ED_Out_Start)
   TYPE(FAST_OutputFileType),      INTENT(IN )  :: y_FAST           !< FAST output file data (for linearization)
   TYPE(ED_OutputType),            INTENT(IN )  :: y_ED             !< ED outputs at t
   INTEGER,                        INTENT(IN )  :: BladeNum         !< blade number to find index for
   INTEGER                                      :: k                !< blade number loop

   INTEGER                                      :: ED_Out_Start     !< starting index of this blade mesh in ElastoDyn outputs

   ED_Out_Start = y_FAST%Lin%Modules(MODULE_ED)%Instance(1)%LinStartIndx(LIN_OUTPUT_COL) ! start of y_ED%BladeLn2Mesh(1)%TranslationDisp field (blade motions in y_ED)
   if (allocated(y_ED%BladeLn2Mesh)) then
      do k = 1,min(BladeNum-1,SIZE(y_ED%BladeLn2Mesh,1)) ! Loop through all blades (p_ED%NumBl)
         ED_Out_Start = ED_Out_Start + y_ED%BladeLn2Mesh(k)%NNodes*18 ! 6 fields with 3 components on each blade
      end do
   end if

END FUNCTION Indx_y_ED_Blade_Start
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine returns the starting index for the y_ED%PlatformPtMesh mesh in the FAST linearization outputs.
FUNCTION Indx_y_ED_Platform_Start(y_ED, y_FAST) RESULT(ED_Out_Start)
   TYPE(FAST_OutputFileType),      INTENT(IN )  :: y_FAST           !< FAST output file data (for linearization)
   TYPE(ED_OutputType),            INTENT(IN )  :: y_ED             !< ED outputs at t

   INTEGER                                      :: ED_Out_Start     !< starting index of this mesh in ElastoDyn outputs

   ED_Out_Start = Indx_y_ED_Blade_Start(y_ED, y_FAST, MaxNBlades+1) ! skip all of the blades to get to start of platform
END FUNCTION Indx_y_ED_Platform_Start
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine returns the starting index for the y_ED%TowerLn2Mesh mesh in the FAST linearization outputs.
FUNCTION Indx_y_ED_Tower_Start(y_ED, y_FAST) RESULT(ED_Out_Start)
   TYPE(FAST_OutputFileType),      INTENT(IN )  :: y_FAST           !< FAST output file data (for linearization)
   TYPE(ED_OutputType),            INTENT(IN )  :: y_ED             !< ED outputs at t

   INTEGER                                      :: ED_Out_Start     !< starting index of this mesh in ElastoDyn outputs

   ED_Out_Start = Indx_y_ED_Platform_Start(y_ED, y_FAST)
   ED_Out_Start = ED_Out_Start + y_ED%PlatformPtMesh%NNodes*18 ! 6 fields with 3 components
END FUNCTION Indx_y_ED_Tower_Start
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine returns the starting index for the y_ED%HubPtMesh mesh in the FAST linearization outputs.
FUNCTION Indx_y_ED_Hub_Start(y_ED, y_FAST) RESULT(ED_Out_Start)
   TYPE(FAST_OutputFileType),      INTENT(IN )  :: y_FAST           !< FAST output file data (for linearization)
   TYPE(ED_OutputType),            INTENT(IN )  :: y_ED             !< ED outputs at t

   INTEGER                                      :: ED_Out_Start     !< starting index of this mesh in ElastoDyn outputs

   ED_Out_Start = Indx_y_ED_Tower_Start(y_ED, y_FAST)
   ED_Out_Start = ED_Out_Start + y_ED%TowerLn2Mesh%NNodes*18 ! 6 fields with 3 components
END FUNCTION Indx_y_ED_Hub_Start
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine returns the starting index for the y_ED%BladeRootMotion(BladeNum) mesh in the FAST linearization outputs.
FUNCTION Indx_y_ED_BladeRoot_Start(y_ED, y_FAST, BladeNum) RESULT(ED_Out_Start)
   TYPE(FAST_OutputFileType),      INTENT(IN )  :: y_FAST           !< FAST output file data (for linearization)
   TYPE(ED_OutputType),            INTENT(IN )  :: y_ED             !< ED outputs at t
   INTEGER,                        INTENT(IN )  :: BladeNum         !< blade number to find index for
   INTEGER                                      :: k                !< blade number loop

   INTEGER                                      :: ED_Out_Start     !< starting index of this blade mesh in ElastoDyn outputs

   ED_Out_Start = Indx_y_ED_Hub_Start(y_ED, y_FAST)
   ED_Out_Start = ED_Out_Start + y_ED%HubPtMotion%NNodes*9 ! 3 fields with 3 components
   
   do k = 1,min(BladeNum-1,size(y_ED%BladeRootMotion))
      ED_Out_Start = ED_Out_Start + y_ED%BladeRootMotion(k)%NNodes*18
   end do
END FUNCTION Indx_y_ED_BladeRoot_Start
!----------------------------------------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------------------------------------
!> This routine returns the starting index for the u_AD%TowerMotion mesh in the FAST linearization inputs.
FUNCTION Indx_u_AD_Tower_Start(u_AD, y_FAST) RESULT(AD_Start)
   TYPE(FAST_OutputFileType),      INTENT(IN )  :: y_FAST           !< FAST output file data (for linearization)
   TYPE(AD_InputType),             INTENT(IN )  :: u_AD             !< AD Inputs at t

   INTEGER                                      :: AD_Start         !< starting index of this mesh in AeroDyn inputs

   AD_Start = y_FAST%Lin%Modules(Module_AD)%Instance(1)%LinStartIndx(LIN_INPUT_COL) 

END FUNCTION Indx_u_AD_Tower_Start
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine returns the starting index for the u_AD%HubMotion mesh in the FAST linearization inputs.
FUNCTION Indx_u_AD_Hub_Start(u_AD, y_FAST) RESULT(AD_Start)
   TYPE(FAST_OutputFileType),      INTENT(IN )  :: y_FAST           !< FAST output file data (for linearization)
   TYPE(AD_InputType),             INTENT(IN )  :: u_AD             !< AD Inputs at t

   INTEGER                                      :: AD_Start         !< starting index of this mesh in AeroDyn inputs

   AD_Start = Indx_u_AD_Tower_Start(u_AD, y_FAST) + u_AD%TowerMotion%NNodes * 9  ! 3 fields (MASKID_TRANSLATIONDISP,MASKID_Orientation,MASKID_TRANSLATIONVel) with 3 components

END FUNCTION Indx_u_AD_Hub_Start
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine returns the starting index for the u_AD%BladeRootMotion(k) mesh in the FAST linearization inputs.
FUNCTION Indx_u_AD_BladeRoot_Start(u_AD, y_FAST, BladeNum) RESULT(AD_Start)
   TYPE(FAST_OutputFileType),      INTENT(IN )  :: y_FAST           !< FAST output file data (for linearization)
   TYPE(AD_InputType),             INTENT(IN )  :: u_AD             !< AD Inputs at t
   INTEGER,                        INTENT(IN )  :: BladeNum         !< blade number to find index for
   INTEGER                                      :: k                !< blade number loop

   INTEGER                                      :: AD_Start         !< starting index of this mesh in AeroDyn inputs

   AD_Start = Indx_u_AD_Hub_Start(u_AD, y_FAST) + u_AD%HubMotion%NNodes * 9  ! 3 fields (MASKID_TRANSLATIONDISP,MASKID_Orientation,MASKID_RotationVel) with 3 components
   
   do k = 1,min(BladeNum-1,size(u_AD%BladeRootMotion))
      AD_Start = AD_Start + u_AD%BladeRootMotion(k)%NNodes * 3 ! 1 field (MASKID_Orientation) with 3 components
   end do
END FUNCTION Indx_u_AD_BladeRoot_Start
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine returns the starting index for the u_AD%BladeMotion(k) mesh in the FAST linearization inputs.
FUNCTION Indx_u_AD_Blade_Start(u_AD, y_FAST, BladeNum) RESULT(AD_Start)
   TYPE(FAST_OutputFileType),      INTENT(IN )  :: y_FAST           !< FAST output file data (for linearization)
   TYPE(AD_InputType),             INTENT(IN )  :: u_AD             !< AD Inputs at t
   INTEGER,                        INTENT(IN )  :: BladeNum         !< blade number to find index for
   INTEGER                                      :: k                !< blade number loop

   INTEGER                                      :: AD_Start         !< starting index of this mesh in AeroDyn inputs

   AD_Start = Indx_u_AD_BladeRoot_Start(u_AD, y_FAST, MaxNBlades+1)
   
   do k = 1,min(BladeNum-1,size(u_AD%BladeMotion))
      AD_Start = AD_Start + u_AD%BladeMotion(k)%NNodes * 9 ! 3 fields (TranslationDisp, MASKID_Orientation, TranslationVel) with 3 components
   end do
END FUNCTION Indx_u_AD_Blade_Start
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine returns the starting index for the u_AD%InflowOnBlade array in the FAST linearization inputs.
FUNCTION Indx_u_AD_BladeInflow_Start(u_AD, y_FAST) RESULT(AD_Start)
   TYPE(FAST_OutputFileType),      INTENT(IN )  :: y_FAST           !< FAST output file data (for linearization)
   TYPE(AD_InputType),             INTENT(IN )  :: u_AD             !< AD Inputs at t

   INTEGER                                      :: AD_Start         !< starting index of this array in AeroDyn inputs

   AD_Start = Indx_u_AD_Blade_Start(u_AD, y_FAST, MaxNBlades+1)

END FUNCTION Indx_u_AD_BladeInflow_Start
!----------------------------------------------------------------------------------------------------------------------------------

END MODULE FAST_Linear