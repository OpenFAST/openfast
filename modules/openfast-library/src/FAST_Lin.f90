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
SUBROUTINE Init_Lin(p_FAST, y_FAST, m_FAST, AD, ED, NumBl, NumBlNodes, ErrStat, ErrMsg)

   TYPE(FAST_ParameterType), INTENT(INOUT) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST              !< Output variables for the glue code
   TYPE(FAST_MiscVarType),   INTENT(INOUT) :: m_FAST              !< Miscellaneous variables
   TYPE(AeroDyn_Data),       INTENT(IN   ) :: AD                  !< AeroDyn data
   TYPE(ElastoDyn_Data),     INTENT(IN   ) :: ED                  !< ElastoDyn data
   INTEGER(IntKi),           INTENT(IN   ) :: NumBl               !< Number of blades (for index into ED,AD input array)
   INTEGER(IntKi),           INTENT(IN   ) :: NumBlNodes          !< Number of blade nodes (for index into AD input array)
   
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None

   INTEGER(IntKi)                          :: i, j, k             ! loop/temp variables
   INTEGER(IntKi)                          :: ThisModule          ! Module ID # 
   INTEGER(IntKi)                          :: NumInstances        ! Number of instances of each module
   INTEGER(IntKi)                          :: NumStates           ! Number of states required for the x_eig arrays
   
   INTEGER(IntKi)                          :: i_u                 ! loop/temp variables
   INTEGER(IntKi)                          :: i_y, i_x            ! loop/temp variables

   INTEGER(IntKi)                          :: NextStart(3)        ! allocated to be size(LinStartIndx)=size(SizeLin); helps compute the next starting index for the module components
   INTEGER(IntKi)                          :: ErrStat2            ! local error status
   CHARACTER(ErrMsgLen)                    :: ErrMsg2             ! local error message
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

      ! HydroDyn is next, if activated:
   if ( p_FAST%CompHydro  == Module_HD ) then 
      p_FAST%Lin_NumMods = p_FAST%Lin_NumMods + 1
      p_FAST%Lin_ModOrder( p_FAST%Lin_NumMods ) = Module_HD
   end if

  
      ! SD or ExtPtfm is next, if activated:
   if ( p_FAST%CompSub  == Module_SD ) then 
      p_FAST%Lin_NumMods = p_FAST%Lin_NumMods + 1
      p_FAST%Lin_ModOrder( p_FAST%Lin_NumMods ) = Module_SD
   else if ( p_FAST%CompSub  == Module_ExtPtfm ) then  
      p_FAST%Lin_NumMods = p_FAST%Lin_NumMods + 1
      p_FAST%Lin_ModOrder( p_FAST%Lin_NumMods ) = Module_ExtPtfm
   end if
   
      ! MAP is next, if activated:
   if ( p_FAST%CompMooring  == Module_MAP ) then 
      p_FAST%Lin_NumMods = p_FAST%Lin_NumMods + 1
      p_FAST%Lin_ModOrder( p_FAST%Lin_NumMods ) = Module_MAP
   else if ( p_FAST%CompMooring  == Module_MD ) then 
      p_FAST%Lin_NumMods = p_FAST%Lin_NumMods + 1
      p_FAST%Lin_ModOrder( p_FAST%Lin_NumMods ) = Module_MD
   end if


   !.....................
   ! determine total number of inputs/outputs/contStates:
   !.....................
   y_FAST%Lin%Glue%SizeLin = 0
   y_FAST%Lin%Glue%NumOutputs = 0
   
   do i = 1,p_FAST%Lin_NumMods
      ThisModule = p_FAST%Lin_ModOrder( i )
      
      do k=1,size(y_FAST%Lin%Modules(ThisModule)%Instance)
         y_FAST%Lin%Modules(ThisModule)%Instance(k)%SizeLin = 0
         if (allocated(y_FAST%Lin%Modules(ThisModule)%Instance(k)%Names_u)) y_FAST%Lin%Modules(ThisModule)%Instance(k)%SizeLin(LIN_INPUT_COL)     = size(y_FAST%Lin%Modules(ThisModule)%Instance(k)%Names_u)
         if (allocated(y_FAST%Lin%Modules(ThisModule)%Instance(k)%Names_y)) y_FAST%Lin%Modules(ThisModule)%Instance(k)%SizeLin(LIN_OUTPUT_COL)    = size(y_FAST%Lin%Modules(ThisModule)%Instance(k)%Names_y)
         if (allocated(y_FAST%Lin%Modules(ThisModule)%Instance(k)%Names_x)) y_FAST%Lin%Modules(ThisModule)%Instance(k)%SizeLin(LIN_ContSTATE_COL) = size(y_FAST%Lin%Modules(ThisModule)%Instance(k)%Names_x)
         
         y_FAST%Lin%Glue%SizeLin = y_FAST%Lin%Glue%SizeLin + y_FAST%Lin%Modules(ThisModule)%Instance(k)%SizeLin ! total number of inputs, outputs, and continuous states
         
         y_FAST%Lin%Glue%NumOutputs = y_FAST%Lin%Glue%NumOutputs + y_FAST%Lin%Modules(ThisModule)%Instance(k)%NumOutputs ! total number of WriteOutputs
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
   call Init_Lin_InputOutput(p_FAST, y_FAST, NumBl, NumBlNodes, ErrStat2, ErrMsg2)
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
   call AllocAry( y_FAST%Lin%Glue%DerivOrder_x, y_FAST%Lin%Glue%SizeLin(LIN_ContSTATE_COL), 'DerivOrder_x', ErrStat2, ErrMsg2)
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
      do k=1,NumInstances
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
      do k=1,NumInstances
         if (NumInstances > 1 .or. trim(y_FAST%Module_Abrev(ThisModule)) == "BD") then
            ModAbrev = TRIM(y_FAST%Module_Abrev(ThisModule))//'_'//trim(num2lstr(k))
         end if

         if (y_FAST%Lin%Modules(ThisModule)%Instance(k)%SizeLin(LIN_ContSTATE_COL) > 0) then
            if (p_FAST%WrVTK == VTK_ModeShapes) then ! allocate these for restart later
               if (ThisModule == Module_ED) then
                  ! ED has only the active DOFs as the continuous states, but to perturb the OP [Perterb_OP()], we need all of the DOFs
                  NumStates = ED%p%NDOF*2
               else
                  NumStates = y_FAST%Lin%Modules(ThisModule)%Instance(k)%SizeLin(LIN_ContSTATE_COL)
               end if

               call AllocAry( y_FAST%Lin%Modules(ThisModule)%Instance(k)%op_x_eig_mag, NumStates, 'op_x_eig_mag', ErrStat2, ErrMsg2)
                  call SetErrStat(errStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
               call AllocAry( y_FAST%Lin%Modules(ThisModule)%Instance(k)%op_x_eig_phase, NumStates, 'op_x_eig_phase', ErrStat2, ErrMsg2)
                  call SetErrStat(errStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
                  
               if (ErrStat >= AbortErrLev) return

               y_FAST%Lin%Modules(ThisModule)%Instance(k)%op_x_eig_mag = 0.0_R8Ki
               y_FAST%Lin%Modules(ThisModule)%Instance(k)%op_x_eig_phase = 0.0_R8Ki
            end if
         end if
         
            
         do j=1,y_FAST%Lin%Modules(ThisModule)%Instance(k)%SizeLin(LIN_ContSTATE_COL)
            y_FAST%Lin%Glue%names_x( i_x) = TRIM(ModAbrev)//' '//y_FAST%Lin%Modules(ThisModule)%Instance(k)%Names_x( j)
            if (allocated(y_FAST%Lin%Modules(ThisModule)%Instance(k)%RotFrame_x)) then
               y_FAST%Lin%Glue%RotFrame_x(i_x) = y_FAST%Lin%Modules(ThisModule)%Instance(k)%RotFrame_x(j) 
            else 
               y_FAST%Lin%Glue%RotFrame_x(i_x) = .false.
            end if
            
            if (allocated(y_FAST%Lin%Modules(ThisModule)%Instance(k)%DerivOrder_x)) then
               y_FAST%Lin%Glue%DerivOrder_x(i_x) = y_FAST%Lin%Modules(ThisModule)%Instance(k)%DerivOrder_x(j) 
            else 
               y_FAST%Lin%Glue%DerivOrder_x(i_x) = 0
            end if
            i_x = i_x + 1;
         end do
      end do
      
   end do ! each module


   !.....................
   ! initialize variables for periodic steady state solution
   !.....................
   
   m_FAST%Lin%NextLinTimeIndx = 1 
   m_FAST%Lin%CopyOP_CtrlCode = MESH_NEWCOPY
   m_FAST%Lin%n_rot           = 0
   m_FAST%Lin%IsConverged     = .false.
   m_FAST%Lin%FoundSteady     = .false.
   m_FAST%Lin%ForceLin        = .false.
   m_FAST%Lin%AzimIndx        = 1
   
   p_FAST%AzimDelta   = TwoPi / p_FAST%NLinTimes
   
      ! allocate space to save operating points
   if (p_FAST%CalcSteady .or. p_FAST%WrVTK==VTK_ModeShapes) then
      
      call AllocateOP(p_FAST, y_FAST, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         
         ! allocate spaces for variables needed to determine 
      if (p_FAST%CalcSteady) then
      
        !call AllocAry(m_FAST%Lin%AzimTarget, p_FAST%NLinTimes,'AzimTarget', ErrStat2, ErrMsg2)
         allocate( m_FAST%Lin%AzimTarget(0 : p_FAST%NLinTimes+1), stat=ErrStat2 )
         if (ErrStat2 /= 0) then
            call SetErrStat(ErrID_Fatal,"Unable to allocate space for AzimTarget.",ErrStat,ErrMsg,RoutineName)
         end if
         
         call AllocAry( m_FAST%Lin%LinTimes, p_FAST%NLinTimes, 'LinTimes', ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

         call AllocAry( m_FAST%Lin%Psi, p_FAST%LinInterpOrder+1, 'Psi', ErrStat2, ErrMsg2)
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

            ! these flattened output arrays will contain spaces for %WriteOutputs, which are being ignored for purposes of CalcSteady computations
         call AllocAry( m_FAST%Lin%y_interp, y_FAST%Lin%Glue%SizeLin(LIN_OUTPUT_COL), 'y_interp', ErrStat2, ErrMsg2)
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

         call AllocAry( m_FAST%Lin%Y_prevRot, y_FAST%Lin%Glue%SizeLin(LIN_OUTPUT_COL), p_FAST%NLinTimes, 'Y_prevRot', ErrStat2, ErrMsg2)
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

         call AllocAry( m_FAST%Lin%y_ref, y_FAST%Lin%Glue%SizeLin(LIN_OUTPUT_COL), 'y_ref', ErrStat2, ErrMsg2)
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

         if (ErrStat < AbortErrLev) then
            m_FAST%Lin%y_interp = 0.0_R8Ki
            m_FAST%Lin%Y_prevRot = 0.0_R8Ki
            m_FAST%Lin%y_ref = 1.0_R8Ki
         end if
         
      end if
      
   end if
   
   
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
                           
         DO K = 1,SIZE(u_AD%rotors(1)%BladeMotion)
            DO J = 1,u_AD%rotors(1)%BladeMotion(k)%Nnodes
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
         DO J=1,u_AD%rotors(1)%TowerMotion%nnodes
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
SUBROUTINE Init_Lin_InputOutput(p_FAST, y_FAST, NumBl, NumBlNodes, ErrStat, ErrMsg)

   TYPE(FAST_ParameterType), INTENT(INOUT) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST              !< Output variables for the glue code
   INTEGER(IntKi),           INTENT(IN   ) :: NumBl               !< Number of blades (for index into ED,AD input array)
   INTEGER(IntKi),           INTENT(IN   ) :: NumBlNodes          !< Number of blades nodes (for index into AD input array)
   
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None

   INTEGER(IntKi)                          :: i, j, col           ! loop/temp variables
   INTEGER(IntKi)                          :: k                   ! loop/temp variables
   INTEGER(IntKi)                          :: ThisModule          ! Module ID # 
   
   INTEGER(IntKi)                          :: ErrStat2            ! local error status
   CHARACTER(ErrMsgLen)                    :: ErrMsg2             ! local error message
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
      
      ! AD standard inputs: UserProp(NumBlNodes,NumBl)
      if (p_FAST%CompAero == MODULE_AD) then
         do j=1,NumBl*NumBlNodes
            y_FAST%Lin%Modules(MODULE_AD)%Instance(1)%use_u(y_FAST%Lin%Modules(MODULE_AD)%Instance(1)%SizeLin(LIN_INPUT_COL)+1-j) = .true.
         end do
      end if
      
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

      ! HD standard inputs: WaveElev0
      if (p_FAST%CompHydro == MODULE_HD) then
            y_FAST%Lin%Modules(MODULE_HD)%Instance(1)%use_u(y_FAST%Lin%Modules(MODULE_HD)%Instance(1)%SizeLin(LIN_INPUT_COL)) = .true.
      end if
      
     ! SD has no standard inputs

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
            col = y_FAST%Lin%Modules(ThisModule)%Instance(k)%SizeLin(LIN_OUTPUT_COL) - y_FAST%Lin%Modules(ThisModule)%Instance(k)%NumOutputs !last column before WriteOutput occurs
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
SUBROUTINE FAST_Linearize_OP(t_global, p_FAST, y_FAST, m_FAST, ED, BD, SrvD, AD, IfW, OpFM, HD, SD, ExtPtfm, MAPp, FEAM, MD, Orca, &
                         IceF, IceD, MeshMapData, ErrStat, ErrMsg )

   REAL(DbKi),               INTENT(IN   ) :: t_global            !< current (global) simulation time

   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST              !< Output variables for the glue code
   TYPE(FAST_MiscVarType),   INTENT(INOUT) :: m_FAST              !< Miscellaneous variables
     
   TYPE(ElastoDyn_Data),     INTENT(INOUT) :: ED                  !< ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(INOUT) :: BD                  !< BeamDyn data
   TYPE(ServoDyn_Data),      INTENT(INOUT) :: SrvD                !< ServoDyn data
   TYPE(AeroDyn_Data),       INTENT(INOUT) :: AD                  !< AeroDyn data
   TYPE(InflowWind_Data),    INTENT(INOUT) :: IfW                 !< InflowWind data
   TYPE(OpenFOAM_Data),      INTENT(INOUT) :: OpFM                !< OpenFOAM data
   TYPE(HydroDyn_Data),      INTENT(INOUT) :: HD                  !< HydroDyn data
   TYPE(SubDyn_Data),        INTENT(INOUT) :: SD                  !< SubDyn data
   TYPE(ExtPtfm_Data),       INTENT(INOUT) :: ExtPtfm             !< ExtPtfm_MCKF data
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
   CHARACTER(ErrMsgLen)                    :: ErrMsg2             ! local error message
   CHARACTER(*),             PARAMETER     :: RoutineName = 'FAST_Linearize_OP' 
   
   REAL(R8Ki), ALLOCATABLE                 :: dUdu(:,:), dUdy(:,:) ! variables for glue-code linearization
   integer(intki)                          :: NumBl
   integer(intki)                          :: k
   CHARACTER(1024)                         :: LinRootName
   CHARACTER(1024)                         :: OutFileName
   CHARACTER(200)                          :: SimStr
   CHARACTER(MaxWrScrLen)                  :: BlankLine
   CHARACTER(*), PARAMETER                 :: Fmt = 'F10.2'
   
   
   
   ErrStat = ErrID_None
   ErrMsg = ""
   Un = -1
   
   !.....................
      SimStr = '(RotSpeed='//trim(num2lstr(ED%y%RotSpeed*RPS2RPM,Fmt))//' rpm, BldPitch1='//trim(num2lstr(ED%y%BlPitch(1)*R2D,Fmt))//' deg)'

   BlankLine = ""
   CALL WrOver( BlankLine )  ! BlankLine contains MaxWrScrLen spaces
   CALL WrOver ( ' Performing linearization '//trim(num2lstr(m_FAST%Lin%NextLinTimeIndx))//' at simulation time '//TRIM( Num2LStr(t_global) )//' s. '//trim(SimStr) )
   CALL WrScr('')
   
   !.....................
   
   LinRootName = TRIM(p_FAST%OutFileRoot)//'.'//trim(num2lstr(m_FAST%Lin%NextLinTimeIndx))
   
   if (p_FAST%WrVTK == VTK_ModeShapes .and. .not. p_FAST%CalcSteady) then ! we already saved these for the CalcSteady case
      call SaveOP(m_FAST%Lin%NextLinTimeIndx, p_FAST, y_FAST, ED, BD, SrvD, AD, IfW, OpFM, HD, SD, ExtPtfm, MAPp, FEAM, MD, Orca, &
                            IceF, IceD, ErrStat, ErrMsg, m_FAST%Lin%CopyOP_CtrlCode )
      !m_FAST%Lin%CopyOP_CtrlCode = MESH_UPDATECOPY ! we need a new copy for each LinTime
   end if

   
      NumBl = size(ED%Input(1)%BlPitchCom) 
      y_FAST%Lin%RotSpeed = ED%y%RotSpeed
      y_FAST%Lin%Azimuth  = ED%y%LSSTipPxa
      !.....................
      ! ElastoDyn
      !.....................
         ! get the jacobians
      call ED_JacobianPInput( t_global, ED%Input(1), ED%p, ED%x(STATE_CURR), ED%xd(STATE_CURR), ED%z(STATE_CURR), ED%OtherSt(STATE_CURR), &
                                 ED%y, ED%m, ErrStat2, ErrMsg2, dYdu=y_FAST%Lin%Modules(Module_ED)%Instance(1)%D, dXdu=y_FAST%Lin%Modules(Module_ED)%Instance(1)%B )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
      call ED_JacobianPContState( t_global, ED%Input(1), ED%p, ED%x(STATE_CURR), ED%xd(STATE_CURR), ED%z(STATE_CURR), ED%OtherSt(STATE_CURR), &
                                     ED%y, ED%m, ErrStat2, ErrMsg2, dYdx=y_FAST%Lin%Modules(Module_ED)%Instance(1)%C, dXdx=y_FAST%Lin%Modules(Module_ED)%Instance(1)%A )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   
         ! get the operating point
      call ED_GetOP( t_global, ED%Input(1), ED%p, ED%x(STATE_CURR), ED%xd(STATE_CURR), ED%z(STATE_CURR), ED%OtherSt(STATE_CURR), &
                        ED%y, ED%m, ErrStat2, ErrMsg2, u_op=y_FAST%Lin%Modules(Module_ED)%Instance(1)%op_u, &
                                                       y_op=y_FAST%Lin%Modules(Module_ED)%Instance(1)%op_y, &
                                                       x_op=y_FAST%Lin%Modules(Module_ED)%Instance(1)%op_x, &
                                                      dx_op=y_FAST%Lin%Modules(Module_ED)%Instance(1)%op_dx )
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
      end if !BeamDyn
      

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
                                   SrvD%OtherSt(STATE_CURR), SrvD%y, SrvD%m, ErrStat2, ErrMsg2,    &
                                   dXdu=y_FAST%Lin%Modules(Module_SrvD)%Instance(1)%B, &
                                   dYdu=y_FAST%Lin%Modules(Module_SrvD)%Instance(1)%D )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
      call SrvD_JacobianPContState( t_global, SrvD%Input(1), SrvD%p, SrvD%x(STATE_CURR), SrvD%xd(STATE_CURR), SrvD%z(STATE_CURR), &
                                       SrvD%OtherSt(STATE_CURR), SrvD%y, SrvD%m, ErrStat2, ErrMsg2, &
                                       dYdx=y_FAST%Lin%Modules(Module_SrvD)%Instance(1)%C, &
                                       dXdx=y_FAST%Lin%Modules(Module_SrvD)%Instance(1)%A )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      ! get the operating point
      call SrvD_GetOP( t_global, SrvD%Input(1), SrvD%p, SrvD%x(STATE_CURR), SrvD%xd(STATE_CURR), SrvD%z(STATE_CURR), &
                       SrvD%OtherSt(STATE_CURR), SrvD%y, SrvD%m, ErrStat2, ErrMsg2, &
                                                       u_op=y_FAST%Lin%Modules(Module_SrvD)%Instance(1)%op_u,  &
                                                      dx_op=y_FAST%Lin%Modules(Module_SrvD)%Instance(1)%op_dx, &
                                                       x_op=y_FAST%Lin%Modules(Module_SrvD)%Instance(1)%op_x,  &
                                                       y_op=y_FAST%Lin%Modules(Module_SrvD)%Instance(1)%op_y   )
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
            ! Jacobians
            call WrPartialMatrix(y_FAST%Lin%Modules(Module_SrvD)%Instance(1)%A, Un, p_FAST%OutFmt, 'dXdx')
            call WrPartialMatrix(y_FAST%Lin%Modules(Module_SrvD)%Instance(1)%B, Un, p_FAST%OutFmt, 'dXdu', UseCol=y_FAST%Lin%Modules(Module_SrvD)%Instance(1)%use_u)
            call WrPartialMatrix(y_FAST%Lin%Modules(Module_SrvD)%Instance(1)%C, Un, p_FAST%OutFmt, 'dYdx', UseRow=y_FAST%Lin%Modules(Module_SrvD)%Instance(1)%use_y)
            call WrPartialMatrix(y_FAST%Lin%Modules(Module_SrvD)%Instance(1)%D, Un, p_FAST%OutFmt, 'dYdu', UseRow=y_FAST%Lin%Modules(Module_SrvD)%Instance(1)%use_y, &
                                                                                                           UseCol=y_FAST%Lin%Modules(Module_SrvD)%Instance(1)%use_u)
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
                                   AD%OtherSt(STATE_CURR), AD%y, AD%m, ErrStat2, ErrMsg2, &
                                   dXdu=y_FAST%Lin%Modules(Module_AD)%Instance(1)%B, &
                                   dYdu=y_FAST%Lin%Modules(Module_AD)%Instance(1)%D )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      call AD_JacobianPContState( t_global, AD%Input(1), AD%p, AD%x(STATE_CURR), AD%xd(STATE_CURR), AD%z(STATE_CURR), &
                                   AD%OtherSt(STATE_CURR), AD%y, AD%m, ErrStat2, ErrMsg2, &
                                   dXdx=y_FAST%Lin%Modules(Module_AD)%Instance(1)%A, &
                                   dYdx=y_FAST%Lin%Modules(Module_AD)%Instance(1)%C )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         
      ! get the operating point
      call AD_GetOP( t_global, AD%Input(1), AD%p, AD%x(STATE_CURR), AD%xd(STATE_CURR), AD%z(STATE_CURR), &
                       AD%OtherSt(STATE_CURR), AD%y, AD%m, ErrStat2, ErrMsg2, &
                       u_op=y_FAST%Lin%Modules(Module_AD)%Instance(1)%op_u, &
                       y_op=y_FAST%Lin%Modules(Module_AD)%Instance(1)%op_y, &
                       x_op=y_FAST%Lin%Modules(Module_AD)%Instance(1)%op_x, &
                      dx_op=y_FAST%Lin%Modules(Module_AD)%Instance(1)%op_dx )
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
            call WrPartialMatrix( y_FAST%Lin%Modules(Module_AD)%Instance(1)%A, Un, p_FAST%OutFmt, 'dXdx' )
                           
            call WrPartialMatrix( y_FAST%Lin%Modules(Module_AD)%Instance(1)%B, Un, p_FAST%OutFmt, 'dXdu', &
                           UseCol=y_FAST%Lin%Modules(Module_AD)%Instance(1)%use_u )
                           
            call WrPartialMatrix( y_FAST%Lin%Modules(Module_AD)%Instance(1)%C, Un, p_FAST%OutFmt, 'dYdx', &
                           UseRow=y_FAST%Lin%Modules(Module_AD)%Instance(1)%use_y )
                           
            call WrPartialMatrix( y_FAST%Lin%Modules(Module_AD)%Instance(1)%D, Un, p_FAST%OutFmt, 'dYdu', &
                           UseRow=y_FAST%Lin%Modules(Module_AD)%Instance(1)%use_y, &
                           UseCol=y_FAST%Lin%Modules(Module_AD)%Instance(1)%use_u )
         end if

            ! finish writing the file
         call WrLinFile_txt_End(Un, p_FAST, y_FAST%Lin%Modules(Module_AD)%Instance(1) )
      end if

   end if
   
   
   !.....................
   ! HydroDyn
   !.....................
   if ( p_FAST%CompHydro  == Module_HD ) then 
         ! get the jacobians
      call HD_JacobianPInput( t_global, HD%Input(1), HD%p, HD%x(STATE_CURR), HD%xd(STATE_CURR), HD%z(STATE_CURR), HD%OtherSt(STATE_CURR), &
                                 HD%y, HD%m, ErrStat2, ErrMsg2, dYdu=y_FAST%Lin%Modules(Module_HD)%Instance(1)%D, dXdu=y_FAST%Lin%Modules(Module_HD)%Instance(1)%B )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
      call HD_JacobianPContState( t_global, HD%Input(1), HD%p, HD%x(STATE_CURR), HD%xd(STATE_CURR), HD%z(STATE_CURR), HD%OtherSt(STATE_CURR), &
                                     HD%y, HD%m, ErrStat2, ErrMsg2, dYdx=y_FAST%Lin%Modules(Module_HD)%Instance(1)%C, dXdx=y_FAST%Lin%Modules(Module_HD)%Instance(1)%A )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   
         ! get the operating point
      call HD_GetOP( t_global, HD%Input(1), HD%p, HD%x(STATE_CURR), HD%xd(STATE_CURR), HD%z(STATE_CURR), HD%OtherSt(STATE_CURR), &
                        HD%y, HD%m, ErrStat2, ErrMsg2, u_op=y_FAST%Lin%Modules(Module_HD)%Instance(1)%op_u, y_op=y_FAST%Lin%Modules(Module_HD)%Instance(1)%op_y, &
                       x_op=y_FAST%Lin%Modules(Module_HD)%Instance(1)%op_x, dx_op=y_FAST%Lin%Modules(Module_HD)%Instance(1)%op_dx )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat >=AbortErrLev) then
            call cleanup()
            return
         end if
         
      
      
         ! write the module matrices:
      if (p_FAST%LinOutMod) then
            
         OutFileName = trim(LinRootName)//'.'//TRIM(y_FAST%Module_Abrev(Module_HD))      
         call WrLinFile_txt_Head(t_global, p_FAST, y_FAST, y_FAST%Lin%Modules(Module_HD)%Instance(1), OutFileName, Un, ErrStat2, ErrMsg2 )
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            if (ErrStat >=AbortErrLev) then
               call cleanup()
               return
            end if
         
         if (p_FAST%LinOutJac) then
            ! Jacobians
            !dXdx:
            call WrPartialMatrix( y_FAST%Lin%Modules(Module_HD)%Instance(1)%A, Un, p_FAST%OutFmt, 'dXdx' )
         
            !dXdu:
            call WrPartialMatrix( y_FAST%Lin%Modules(Module_HD)%Instance(1)%B, Un, p_FAST%OutFmt, 'dXdu', UseCol=y_FAST%Lin%Modules(Module_HD)%Instance(1)%use_u )
         
            !dYdx:
            call WrPartialMatrix( y_FAST%Lin%Modules(Module_HD)%Instance(1)%C, Un, p_FAST%OutFmt, 'dYdx', UseRow=y_FAST%Lin%Modules(Module_HD)%Instance(1)%use_y )
         
            !dYdu:
            call WrPartialMatrix( y_FAST%Lin%Modules(Module_HD)%Instance(1)%D, Un, p_FAST%OutFmt, 'dYdu', UseRow=y_FAST%Lin%Modules(Module_HD)%Instance(1)%use_y, &
                                                                                                          UseCol=y_FAST%Lin%Modules(Module_HD)%Instance(1)%use_u )
         
         end if 
      
             ! finish writing the file
         call WrLinFile_txt_End(Un, p_FAST, y_FAST%Lin%Modules(Module_HD)%Instance(1) )
               
      end if  
   end if
   
   !.....................
   ! SubDyn / ExtPtfm
   !.....................

   if ( p_FAST%CompSub == Module_SD ) then 
      ! get the jacobians
      call SD_JacobianPInput( t_global, SD%Input(1), SD%p, SD%x(STATE_CURR), SD%xd(STATE_CURR), &
              SD%z(STATE_CURR), SD%OtherSt(STATE_CURR),  SD%y, SD%m, ErrStat2, ErrMsg2, &
              dYdu=y_FAST%Lin%Modules(Module_SD)%Instance(1)%D, dXdu=y_FAST%Lin%Modules(Module_SD)%Instance(1)%B )
      if(Failed()) return;

      call SD_JacobianPContState( t_global, SD%Input(1), SD%p, SD%x(STATE_CURR), SD%xd(STATE_CURR), &
          SD%z(STATE_CURR), SD%OtherSt(STATE_CURR), SD%y, SD%m, ErrStat2, ErrMsg2,&
          dYdx=y_FAST%Lin%Modules(Module_SD)%Instance(1)%C, dXdx=y_FAST%Lin%Modules(Module_SD)%Instance(1)%A )
      if(Failed()) return;

      ! get the operating point
      call SD_GetOP(t_global, SD%Input(1), SD%p, SD%x(STATE_CURR), SD%xd(STATE_CURR), SD%z(STATE_CURR),&
          SD%OtherSt(STATE_CURR), SD%y, SD%m, ErrStat2, ErrMsg2, u_op=y_FAST%Lin%Modules(Module_SD)%Instance(1)%op_u,&
          y_op=y_FAST%Lin%Modules(Module_SD)%Instance(1)%op_y, &
          x_op=y_FAST%Lin%Modules(Module_SD)%Instance(1)%op_x, dx_op=y_FAST%Lin%Modules(Module_SD)%Instance(1)%op_dx)
      if(Failed()) return;

      ! write the module matrices:
      if (p_FAST%LinOutMod) then
         OutFileName = trim(LinRootName)//'.'//TRIM(y_FAST%Module_Abrev(Module_SD))      
         call WrLinFile_txt_Head(t_global, p_FAST, y_FAST, y_FAST%Lin%Modules(Module_SD)%Instance(1), OutFileName, Un, ErrStat2, ErrMsg2)
         if(Failed()) return;
            
         if (p_FAST%LinOutJac) then
            ! Jacobians
            call WrPartialMatrix(y_FAST%Lin%Modules(Module_SD)%Instance(1)%A, Un, p_FAST%OutFmt, 'dXdx')
            call WrPartialMatrix(y_FAST%Lin%Modules(Module_SD)%Instance(1)%B, Un, p_FAST%OutFmt, 'dXdu', UseCol=y_FAST%Lin%Modules(Module_SD)%Instance(1)%use_u)
            call WrPartialMatrix(y_FAST%Lin%Modules(Module_SD)%Instance(1)%C, Un, p_FAST%OutFmt, 'dYdx', UseRow=y_FAST%Lin%Modules(Module_SD)%Instance(1)%use_y)
            call WrPartialMatrix(y_FAST%Lin%Modules(Module_SD)%Instance(1)%D, Un, p_FAST%OutFmt, 'dYdu', UseRow=y_FAST%Lin%Modules(Module_SD)%Instance(1)%use_y, &
                                                                                                               UseCol=y_FAST%Lin%Modules(Module_SD)%Instance(1)%use_u)
         end if
         ! finish writing the file
         call WrLinFile_txt_End(Un, p_FAST, y_FAST%Lin%Modules(Module_SD)%Instance(1) )
      end if
   elseif ( p_FAST%CompSub == Module_ExtPtfm ) then 
      ! get the jacobians
      call ExtPtfm_JacobianPInput( t_global, ExtPtfm%Input(1), ExtPtfm%p, ExtPtfm%x(STATE_CURR), ExtPtfm%xd(STATE_CURR), &
              ExtPtfm%z(STATE_CURR), ExtPtfm%OtherSt(STATE_CURR),  ExtPtfm%y, ExtPtfm%m, ErrStat2, ErrMsg2, &
              dYdu=y_FAST%Lin%Modules(Module_ExtPtfm)%Instance(1)%D, dXdu=y_FAST%Lin%Modules(Module_ExtPtfm)%Instance(1)%B )
      if(Failed()) return;

      call ExtPtfm_JacobianPContState( t_global, ExtPtfm%Input(1), ExtPtfm%p, ExtPtfm%x(STATE_CURR), ExtPtfm%xd(STATE_CURR), &
          ExtPtfm%z(STATE_CURR), ExtPtfm%OtherSt(STATE_CURR), ExtPtfm%y, ExtPtfm%m, ErrStat2, ErrMsg2,&
          dYdx=y_FAST%Lin%Modules(Module_ExtPtfm)%Instance(1)%C, dXdx=y_FAST%Lin%Modules(Module_ExtPtfm)%Instance(1)%A )
      if(Failed()) return;

      ! get the operating point
      call ExtPtfm_GetOP(t_global, ExtPtfm%Input(1), ExtPtfm%p, ExtPtfm%x(STATE_CURR), ExtPtfm%xd(STATE_CURR), ExtPtfm%z(STATE_CURR),&
          ExtPtfm%OtherSt(STATE_CURR), ExtPtfm%y, ExtPtfm%m, ErrStat2, ErrMsg2, u_op=y_FAST%Lin%Modules(Module_ExtPtfm)%Instance(1)%op_u,&
          y_op=y_FAST%Lin%Modules(Module_ExtPtfm)%Instance(1)%op_y, &
          x_op=y_FAST%Lin%Modules(Module_ExtPtfm)%Instance(1)%op_x, dx_op=y_FAST%Lin%Modules(Module_ExtPtfm)%Instance(1)%op_dx)
      if(Failed()) return;

      ! write the module matrices:
      if (p_FAST%LinOutMod) then
         OutFileName = trim(LinRootName)//'.'//TRIM(y_FAST%Module_Abrev(Module_ExtPtfm))      
         call WrLinFile_txt_Head(t_global, p_FAST, y_FAST, y_FAST%Lin%Modules(Module_ExtPtfm)%Instance(1), OutFileName, Un, ErrStat2, ErrMsg2)
         if(Failed()) return;
            
         if (p_FAST%LinOutJac) then
            ! Jacobians
            call WrPartialMatrix(y_FAST%Lin%Modules(Module_ExtPtfm)%Instance(1)%A, Un, p_FAST%OutFmt, 'dXdx')
            call WrPartialMatrix(y_FAST%Lin%Modules(Module_ExtPtfm)%Instance(1)%B, Un, p_FAST%OutFmt, 'dXdu', UseCol=y_FAST%Lin%Modules(Module_ExtPtfm)%Instance(1)%use_u)
            call WrPartialMatrix(y_FAST%Lin%Modules(Module_ExtPtfm)%Instance(1)%C, Un, p_FAST%OutFmt, 'dYdx', UseRow=y_FAST%Lin%Modules(Module_ExtPtfm)%Instance(1)%use_y)
            call WrPartialMatrix(y_FAST%Lin%Modules(Module_ExtPtfm)%Instance(1)%D, Un, p_FAST%OutFmt, 'dYdu', UseRow=y_FAST%Lin%Modules(Module_ExtPtfm)%Instance(1)%use_y, &
                                                                                                               UseCol=y_FAST%Lin%Modules(Module_ExtPtfm)%Instance(1)%use_u)
         end if
         ! finish writing the file
         call WrLinFile_txt_End(Un, p_FAST, y_FAST%Lin%Modules(Module_ExtPtfm)%Instance(1) )
     end if
   end if ! SubDyn/ExtPtfm
   
   
   !.....................
   ! MAP
   !.....................
   if ( p_FAST%CompMooring  == Module_MAP ) then
      ! LIN-TODO: We need this to compute the dYdu total derivative which is D for MAP, and the template uses OtherSt(STATE_CURR), but the FAST MAP DATA has OtherSt as a scalar
      call MAP_JacobianPInput( t_global, MAPp%Input(1), MAPp%p, MAPp%x(STATE_CURR), MAPp%xd(STATE_CURR), MAPp%z(STATE_CURR), &
                                   MAPp%OtherSt, MAPp%y, ErrStat2, ErrMsg2, y_FAST%Lin%Modules(Module_MAP)%Instance(1)%D )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
      ! get the operating point
      !LIN-TODO: template uses OtherSt(STATE_CURR), but the FAST MAP DATA has OtherSt as a scalar
      !          email bonnie for a discussion on this.
      call MAP_GetOP( t_global, MAPp%Input(1), MAPp%p, MAPp%x(STATE_CURR), MAPp%xd(STATE_CURR), MAPp%z(STATE_CURR), &
                             MAPp%OtherSt, MAPp%y, ErrStat2, ErrMsg2,  &
                       y_FAST%Lin%Modules(Module_MAP)%Instance(1)%op_u, y_FAST%Lin%Modules(Module_MAP)%Instance(1)%op_y )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (ErrStat >=AbortErrLev) then
         call cleanup()
         return
      end if
      
      ! write the module matrices:
      if (p_FAST%LinOutMod) then
            
         OutFileName = trim(LinRootName)//'.'//TRIM(y_FAST%Module_Abrev(Module_MAP))      
         call WrLinFile_txt_Head(t_global, p_FAST, y_FAST, y_FAST%Lin%Modules(Module_MAP)%Instance(1), OutFileName, Un, ErrStat2, ErrMsg2 )
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            if (ErrStat >=AbortErrLev) then
               call cleanup()
               return
            end if
         
         if (p_FAST%LinOutJac) then
            ! Jacobians
            !dYdu:
            call WrPartialMatrix( y_FAST%Lin%Modules(Module_MAP)%Instance(1)%D, Un, p_FAST%OutFmt, 'dYdu', &
                                         UseRow=y_FAST%Lin%Modules(Module_MAP)%Instance(1)%use_y, &
                                         UseCol=y_FAST%Lin%Modules(Module_MAP)%Instance(1)%use_u )         
         end if 
      
             ! finish writing the file
         call WrLinFile_txt_End(Un, p_FAST, y_FAST%Lin%Modules(Module_MAP)%Instance(1) )
               
      end if  ! if ( p_FAST%LinOutMod )
   end if     ! if ( p_FAST%CompMooring  == Module_MAP )
   
   
   !.....................
   ! MoorDyn
   !.....................
   if ( p_FAST%CompMooring  == Module_MD ) then
      
      call MD_JacobianPInput( t_global, MD%Input(1), MD%p, MD%x(STATE_CURR), MD%xd(STATE_CURR), MD%z(STATE_CURR), &
                                   MD%OtherSt(STATE_CURR), MD%y, MD%m, ErrStat2, ErrMsg2, &
                                   dXdu=y_FAST%Lin%Modules(Module_MD)%Instance(1)%B, &
                                   dYdu=y_FAST%Lin%Modules(Module_MD)%Instance(1)%D )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
      call MD_JacobianPContState( t_global, MD%Input(1), MD%p, MD%x(STATE_CURR), MD%xd(STATE_CURR), MD%z(STATE_CURR), MD%OtherSt(STATE_CURR), &
                                     MD%y, MD%m, ErrStat2, ErrMsg2, dYdx=y_FAST%Lin%Modules(Module_MD)%Instance(1)%C, &
                                                                    dXdx=y_FAST%Lin%Modules(Module_MD)%Instance(1)%A )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         
      ! get the operating point
      call MD_GetOP( t_global, MD%Input(1), MD%p, MD%x(STATE_CURR), MD%xd(STATE_CURR), MD%z(STATE_CURR), &
                             MD%OtherSt(STATE_CURR), MD%y, MD%m, ErrStat2, ErrMsg2,  &
                              u_op=y_FAST%Lin%Modules(Module_MD)%Instance(1)%op_u, &
                              y_op=y_FAST%Lin%Modules(Module_MD)%Instance(1)%op_y, &
                              x_op=y_FAST%Lin%Modules(Module_MD)%Instance(1)%op_x, &
                             dx_op=y_FAST%Lin%Modules(Module_MD)%Instance(1)%op_dx )                       
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (ErrStat >=AbortErrLev) then
         call cleanup()
         return
      end if
      
      ! write the module matrices:
      if (p_FAST%LinOutMod) then
            
         OutFileName = trim(LinRootName)//'.'//TRIM(y_FAST%Module_Abrev(Module_MD))      
         call WrLinFile_txt_Head(t_global, p_FAST, y_FAST, y_FAST%Lin%Modules(Module_MD)%Instance(1), OutFileName, Un, ErrStat2, ErrMsg2 )
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            if (ErrStat >=AbortErrLev) then
               call cleanup()
               return
            end if
         
         if (p_FAST%LinOutJac) then
            ! Jacobians
            ! dXdx, dXdu, dYdx, dYdu:
            call WrPartialMatrix( y_FAST%Lin%Modules(Module_MD)%Instance(1)%A, Un, p_FAST%OutFmt, 'dXdx' )
            call WrPartialMatrix( y_FAST%Lin%Modules(Module_MD)%Instance(1)%B, Un, p_FAST%OutFmt, 'dXdu', UseCol=y_FAST%Lin%Modules(Module_MD)%Instance(1)%use_u )
            call WrPartialMatrix( y_FAST%Lin%Modules(Module_MD)%Instance(1)%C, Un, p_FAST%OutFmt, 'dYdx', UseRow=y_FAST%Lin%Modules(Module_MD)%Instance(1)%use_y )
            call WrPartialMatrix( y_FAST%Lin%Modules(Module_MD)%Instance(1)%D, Un, p_FAST%OutFmt, 'dYdu', UseRow=y_FAST%Lin%Modules(Module_MD)%Instance(1)%use_y, &
                                                                                                          UseCol=y_FAST%Lin%Modules(Module_MD)%Instance(1)%use_u )         
         end if 
      
             ! finish writing the file
         call WrLinFile_txt_End(Un, p_FAST, y_FAST%Lin%Modules(Module_MD)%Instance(1) )
               
      end if  ! if ( p_FAST%LinOutMod )
   end if     ! if ( p_FAST%CompMooring  == Module_MD )
   
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
   call Glue_Jacobians( p_FAST, y_FAST, m_FAST, ED, BD, SrvD, AD, IfW, OpFM, HD, SD, MAPp, FEAM, MD, Orca, &
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

   m_FAST%Lin%NextLinTimeIndx = m_FAST%Lin%NextLinTimeIndx + 1

contains
    logical function Failed()
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      Failed =  ErrStat >= AbortErrLev
      if(Failed) call cleanup()
   end function Failed
   subroutine cleanup()
      
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
   CHARACTER(ErrMsgLen)                    :: ErrMsg2             ! local error message
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
   Desc = 'Wind Speed:';      WRITE (Un, fmt) Desc, y_FAST%Lin%WindSpeed,  'm/s'
   
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
      call WrLinFile_txt_Table(p_FAST, Un, "Row/Column", LinData%op_x, LinData%names_x, rotFrame=LinData%RotFrame_x, derivOrder=LinData%DerivOrder_x  )      
      
      WRITE(Un, '(A)') 'Order of continuous state derivatives:'      
      call WrLinFile_txt_Table(p_FAST, Un, "Row/Column", LinData%op_dx, LinData%names_x, rotFrame=LinData%RotFrame_x, deriv=.true., derivOrder=LinData%DerivOrder_x  )      
   end if
   
   if (n(Indx_xd) > 0) then
      WRITE(Un, '(A)') 'Order of discrete states:'      
      call WrLinFile_txt_Table(p_FAST, Un, "Row/Column", LinData%op_xd, LinData%names_xd )      
   end if

   if (n(Indx_z) > 0) then
      WRITE(Un, '(A)') 'Order of constraint states:' 
      call WrLinFile_txt_Table(p_FAST, Un, "Row/Column", LinData%op_z, LinData%names_z, rotFrame=LinData%RotFrame_z )      
   end if
         
   if (n(Indx_u) > 0) then
      WRITE(Un, '(A)') 'Order of inputs:'   
      call WrLinFile_txt_Table(p_FAST, Un, "Column  ", LinData%op_u, LinData%names_u, rotFrame=LinData%RotFrame_u, UseCol=LinData%use_u )
   end if
   
   if (n(Indx_y) > 0) then
      WRITE(Un, '(A)') 'Order of outputs:'      
      call WrLinFile_txt_Table(p_FAST, Un, "Row  ", LinData%op_y, LinData%names_y, rotFrame=LinData%RotFrame_y, UseCol=LinData%use_y )      
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
SUBROUTINE WrLinFile_txt_Table(p_FAST, Un, RowCol, op, names, rotFrame, deriv, derivOrder, UseCol,start_indx)

   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< parameters
   INTEGER(IntKi),           INTENT(IN   ) :: Un                  !< unit number
   CHARACTER(*),             INTENT(IN   ) :: RowCol              !< Row/Column description
   REAL(ReKi),               INTENT(IN   ) :: op(:)               !< operating point values (possibly different size that Desc because of orientations)
   CHARACTER(LinChanLen),    INTENT(IN   ) :: names(:)            !< Descriptions of the channels (names and units)
   logical, optional,        INTENT(IN   ) :: rotFrame(:)         !< determines if this parameter is in the rotating frame
   logical, optional,        intent(in   ) :: deriv               !< flag that tells us if we need to modify the channel names for derivatives (xdot)
   integer(IntKi), optional, intent(in   ) :: derivOrder(:)       !< Order of the time derivatives associated with the channel
   logical, optional,        intent(in   ) :: UseCol(:)           !< flags that tell us if we should use each column or skip it
   INTEGER(IntKi),optional,  INTENT(IN   ) :: start_indx          !< starting index (so extended inputs can be numbered starting after the # of inputs)
   
      ! local variables
   INTEGER(IntKi)                          :: TS                  ! tab stop column
   INTEGER(IntKi)                          :: i, i_print          ! loop counter
   INTEGER(IntKi)                          :: i_op                ! loop counter
   
   logical                                 :: UseDerivNames       !< flag that tells us if we need to modify the channel names for derivatives (xdot)
   logical                                 :: UseThisCol          !< flag that tells us if we should use this particular column or skip it
   logical                                 :: RotatingCol         !< flag that tells us if this column is in the rotating frame
   integer(IntKi)                          :: DerivOrdCol         !< integer indicating the maximum time-derivative order of a channel (this will be 0 for anything that is not a continuous state)
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
   
   Fmt       = '(3x,I8,3x,'//trim(p_FAST%OutFmt)//',T'//trim(num2lstr(TS))//',L8,8x,I8,9x,A)'
   FmtOrient = '(3x,I8,3x,'//trim(p_FAST%OutFmt)//',2(", ",'//trim(p_FAST%OutFmt)//'),T'//trim(num2lstr(TS))//',L8,8x,I8,9x,A)'
   Fmt_Str   = '(3x,A10,1x,A,T'//trim(num2lstr(TS))//',A15,1x,A16,1x,A)'
   
   WRITE(Un, Fmt_Str) RowCol,      'Operating Point', 'Rotating Frame?', 'Derivative Order', 'Description'
   WRITE(Un, Fmt_Str) '----------','---------------', '---------------', '----------------', '-----------'
   
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
      
      DerivOrdCol = 0
      if (present(derivOrder)) DerivOrdCol = derivOrder(i)
      
      RotatingCol = .false.
      if (present(rotFrame)) RotatingCol = rotFrame(i)
                  
      if (index(names(i), ' orientation angle, node ') > 0 ) then  ! make sure this matches what is written in PackMotionMesh_Names()
         if (UseThisCol) then
            WRITE(Un, FmtOrient) i_print, op(i_op), op(i_op+1), op(i_op+2), RotatingCol, DerivOrdCol, trim(names(i))  !//' [OP is a row of the DCM]
            i_print = i_print + 1
         end if
         
         i_op = i_op + 3
      else
         if (UseThisCol) then
            if (UseDerivNames) then
               WRITE(Un, Fmt) i_print, op(i_op), RotatingCol, DerivOrdCol, 'First time derivative of '//trim(names(i))//'/s'
            else
               WRITE(Un, Fmt) i_print, op(i_op), RotatingCol, DerivOrdCol, trim(names(i))
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
   CHARACTER(ErrMsgLen)                    :: ErrMsg2             ! local error message
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
SUBROUTINE Glue_Jacobians( p_FAST, y_FAST, m_FAST, ED, BD, SrvD, AD, IfW, OpFM, HD, SD, MAPp, FEAM, MD, Orca, &
                         IceF, IceD, MeshMapData, dUdu, dUdy, ErrStat, ErrMsg )

   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST              !< Output variables for the glue code
   TYPE(FAST_MiscVarType),   INTENT(INOUT) :: m_FAST              !< Miscellaneous variables
     
   TYPE(ElastoDyn_Data),     INTENT(INOUT) :: ED                  !< ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(INOUT) :: BD                  !< BeamDyn data
   TYPE(ServoDyn_Data),      INTENT(INOUT) :: SrvD                !< ServoDyn data
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
   CHARACTER(ErrMsgLen)                    :: ErrMsg2             ! local error message
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
! LIN-TODO: Add doc strings for new modules: SrvD & HD & MAP & SD
   
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
      ! \f$ \frac{\partial U_\Lambda^{SrvD}}{\partial u^{SrvD}} \end{bmatrix} = \f$ (dUdu block row 2=SrvD)
      !............
   if (p_FAST%CompServo == MODULE_SrvD) then
      call Linear_SrvD_InputSolve_du( p_FAST, y_FAST, SrvD%p, SrvD%Input(1), ED%y, BD, SD, MeshMapData, dUdu, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   end if
   
      !............ 
      ! \f$ \frac{\partial U_\Lambda^{ED}}{\partial u^{SrvD}} \end{bmatrix} = \f$ and 
      ! \f$ \frac{\partial U_\Lambda^{ED}}{\partial u^{BD}}   \end{bmatrix} = \f$ and 
      ! \f$ \frac{\partial U_\Lambda^{ED}}{\partial u^{AD}}   \end{bmatrix} = \f$ and 
      ! \f$ \frac{\partial U_\Lambda^{ED}}{\partial u^{HD}}   \end{bmatrix} = \f$ and 
      ! \f$ \frac{\partial U_\Lambda^{ED}}{\partial u^{SD}}   \end{bmatrix} = \f$ and 
      ! \f$ \frac{\partial U_\Lambda^{ED}}{\partial u^{MAP}}  \end{bmatrix} = \f$ (dUdu block row 3=ED)
      !............   
   ! we need to do this for CompElast=ED and CompElast=BD

   call Linear_ED_InputSolve_du( p_FAST, y_FAST, SrvD, ED%Input(1), ED%y, AD%y, AD%Input(1), BD, HD, SD, MAPp, MD, MeshMapData, dUdu, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   
      !............
      ! \f$ \frac{\partial U_\Lambda^{BD}}{\partial u^{BD}} \end{bmatrix} = \f$ and 
      ! \f$ \frac{\partial U_\Lambda^{BD}}{\partial u^{AD}} \end{bmatrix} = \f$ (dUdu block row 4=BD)
      !............   
   IF (p_FAST%CompElast == Module_BD) THEN
      call Linear_BD_InputSolve_du( p_FAST, y_FAST, SrvD, ED%y, AD%y, AD%Input(1), BD, MeshMapData, dUdu, ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   END IF
      
      !............
      ! \f$ \frac{\partial U_\Lambda^{AD}}{\partial u^{AD}} \end{bmatrix} = \f$ (dUdu block row 5=AD)
      !............
   IF (p_FAST%CompAero == MODULE_AD) THEN 
      call Linear_AD_InputSolve_du( p_FAST, y_FAST, AD%Input(1), ED%y, BD, MeshMapData, dUdu, ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   end if 

  

      !............
      ! \f$ \frac{\partial U_\Lambda^{HD}}{\partial u^{HD}} \end{bmatrix} = \f$ (dUdu block row 6=HD)
      !............
   IF (p_FAST%CompHydro == MODULE_HD) THEN 
      call Linear_HD_InputSolve_du( p_FAST, y_FAST, HD%Input(1), ED%y, SD%y, MeshMapData, dUdu, ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   end if 

      !............
      ! \f$ \frac{\partial U_\Lambda^{SD}}{\partial u^{HD}} \end{bmatrix} = \f$ and  
      ! \f$ \frac{\partial U_\Lambda^{SD}}{\partial u^{SD}} \end{bmatrix} = \f$ and  
      ! \f$ \frac{\partial U_\Lambda^{SD}}{\partial u^{MAP}} \end{bmatrix} = \f$ (dUdu block row 7=SD)
      !............
   IF (p_FAST%CompSub == MODULE_SD) THEN 
      call Linear_SD_InputSolve_du( p_FAST, y_FAST, SrvD, SD%Input(1), SD%y, ED%y, HD, MAPp, MD, MeshMapData, dUdu, ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   ELSE IF (p_FAST%CompSub == Module_ExtPtfm) THEN
       CALL WrScr('>>> FAST_LIN: Linear_ExtPtfm_InputSolve_du, TODO')
   ENDIF 
   
      !............
      ! \f$ \frac{\partial U_\Lambda^{MD}}{\partial u^{MD}} \end{bmatrix} = \f$ (dUdu block row 9=MD) <<<<
      !............
   if (p_FAST%CompMooring == MODULE_MD) then
      call Linear_MD_InputSolve_du( p_FAST, y_FAST, MD%Input(1), ED%y, SD%y, MeshMapData, dUdu, ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   end if

   ! LIN-TODO: Update the doc lines below to include SrvD, HD, SD, and MAP
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
      call Linear_SrvD_InputSolve_dy( p_FAST, y_FAST, SrvD%p, SrvD%Input(1), ED%y, BD, SD%y, MeshMapData, dUdy, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   end if
   

      !............
      ! \f$ \frac{\partial U_\Lambda^{ED}}{\partial y^{SrvD}} \end{bmatrix} = \f$
      ! \f$ \frac{\partial U_\Lambda^{ED}}{\partial y^{ED}}   \end{bmatrix} = \f$
      ! \f$ \frac{\partial U_\Lambda^{ED}}{\partial y^{BD}}   \end{bmatrix} = \f$
      ! \f$ \frac{\partial U_\Lambda^{ED}}{\partial y^{AD}}   \end{bmatrix} = \f$
      ! \f$ \frac{\partial U_\Lambda^{ED}}{\partial y^{HD}}   \end{bmatrix} = \f$
      ! \f$ \frac{\partial U_\Lambda^{ED}}{\partial y^{SD}}   \end{bmatrix} = \f$
      ! \f$ \frac{\partial U_\Lambda^{ED}}{\partial y^{MAP}}  \end{bmatrix} = \f$ (dUdy block row 3=ED)
      !............

   call Linear_ED_InputSolve_dy( p_FAST, y_FAST, SrvD, ED%Input(1), ED%y, AD%y, AD%Input(1), BD, HD, SD, MAPp, MD, MeshMapData, dUdy, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   
      !............
      ! \f$ \frac{\partial U_\Lambda^{BD}}{\partial y^{ED}} \end{bmatrix} = \f$
      ! \f$ \frac{\partial U_\Lambda^{BD}}{\partial y^{BD}} \end{bmatrix} = \f$
      ! \f$ \frac{\partial U_\Lambda^{BD}}{\partial y^{AD}} \end{bmatrix} = \f$ (dUdy block row 4=BD)
      !............
   if (p_FAST%CompElast == MODULE_BD) then
      call Linear_BD_InputSolve_dy( p_FAST, y_FAST, SrvD, ED%Input(1), ED%y, AD%y, AD%Input(1), BD, MeshMapData, dUdy, ErrStat2, ErrMsg2 )
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

      call Linear_AD_InputSolve_NoIfW_dy( p_FAST, y_FAST, AD%Input(1), ED%y, BD, MeshMapData, dUdy, ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

   end if


      !............
      ! \f$ \frac{\partial U_\Lambda^{SD}}{\partial y^{ED}} \end{bmatrix} = \f$
      ! \f$ \frac{\partial U_\Lambda^{SD}}{\partial y^{HD}} \end{bmatrix} = \f$
      ! \f$ \frac{\partial U_\Lambda^{SD}}{\partial y^{SD}} \end{bmatrix} = \f$
      ! \f$ \frac{\partial U_\Lambda^{SD}}{\partial y^{MAP}} \end{bmatrix} = \f$ (dUdy block row 7=SD)
      !............
   if (p_FAST%CompHydro == MODULE_HD) then
      call Linear_HD_InputSolve_dy( p_FAST, y_FAST, HD%Input(1), ED%y, SD%y, MeshMapData, dUdy, ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   end if
   
   !LIN-TODO: Add doc strings and look at above doc string
   IF (p_FAST%CompSub == Module_SD) THEN
       call Linear_SD_InputSolve_dy( p_FAST, y_FAST, SrvD, SD%Input(1), SD%y, ED%y, HD, MAPp, MD, MeshMapData, dUdy, ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   ELSE IF (p_FAST%CompSub == Module_ExtPtfm) THEN
       write(*,*)'>>> FAST_LIN: Linear_ExtPtfm_InputSolve_dy, TODO'
   ENDIF
   
      !............
      ! \f$ \frac{\partial U_\Lambda^{MAP}}{\partial y^{ED}} \end{bmatrix} = \f$
      ! \f$ \frac{\partial U_\Lambda^{MAP}}{\partial y^{SD}} \end{bmatrix} = \f$ (dUdy block row 8=MAP)
      !............
   if (p_FAST%CompMooring == MODULE_MAP) then
      call Linear_MAP_InputSolve_dy( p_FAST, y_FAST, MAPp%Input(1), ED%y, SD%y, MeshMapData, dUdy, ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   end if
      !............
      ! \f$ \frac{\partial U_\Lambda^{MD}}{\partial y^{ED}} \end{bmatrix} = \f$
      ! \f$ \frac{\partial U_\Lambda^{MD}}{\partial y^{SD}} \end{bmatrix} = \f$ (dUdy block row 9=MD) <<<<
      !............
   if (p_FAST%CompMooring == MODULE_MD) then
      call Linear_MD_InputSolve_dy( p_FAST, y_FAST, MD%Input(1), ED%y, SD%y, MeshMapData, dUdy, ErrStat2, ErrMsg2 )
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
                     + u_AD%rotors(1)%TowerMotion%NNodes * 9  & ! 3 fields (MASKID_TRANSLATIONDISP,MASKID_Orientation,MASKID_TRANSLATIONVel) with 3 components
                     + u_AD%rotors(1)%HubMotion%NNodes   * 9    ! 3 fields (MASKID_TRANSLATIONDISP,MASKID_Orientation,MASKID_RotationVel) with 3 components
   
         do k = 1,size(u_AD%rotors(1)%BladeRootMotion)         
            AD_Start_Bl = AD_Start_Bl + u_AD%rotors(1)%BladeRootMotion(k)%NNodes * 3 ! 1 field (MASKID_Orientation) with 3 components
         end do
         ! next is u_AD%BladeMotion(k):
                  
         DO K = 1,SIZE(u_AD%rotors(1)%BladeMotion)
            DO J = 1,u_AD%rotors(1)%BladeMotion(k)%Nnodes
               Node = Node + 1 ! InflowWind node
               do i=1,3 !XYZ components of this node
                  i2 = y_FAST%Lin%Modules(MODULE_IfW)%Instance(1)%LinStartIndx(LIN_INPUT_COL) + (Node-1)*3 + i - 1
                  j2 = AD_Start_Bl + (j-1)*3 + i - 1
                  dUdu( i2, j2 ) = -1.0_R8Ki
               end do            
            END DO !J = 1,p%BldNodes ! Loop through the blade nodes / elements
                     
               ! get starting AD index of BladeMotion for next blade
            AD_Start_Bl = AD_Start_Bl + u_AD%rotors(1)%BladeMotion(k)%Nnodes * 9  ! 3 fields (MASKID_TRANSLATIONDISP,MASKID_Orientation,MASKID_TRANSLATIONVel) with 3 components
         END DO !K = 1,p%NumBl     
         
            ! tower:
         DO J=1,u_AD%rotors(1)%TowerMotion%nnodes
            Node = Node + 1   
            do i=1,3 !XYZ components of this node
               i2 = y_FAST%Lin%Modules(MODULE_IfW)%Instance(1)%LinStartIndx(LIN_INPUT_COL) + (Node-1)*3 + i - 1
               j2 = y_FAST%Lin%Modules(MODULE_AD )%Instance(1)%LinStartIndx(LIN_INPUT_COL) +    (j-1)*3 + i - 1
               dUdu( i2, j2 ) = -1.0_R8Ki
            end do            
         END DO              
         
         ! HubPosition and HubOrientation from ElastoDyn are missing from this
      END IF     
END SUBROUTINE Linear_IfW_InputSolve_du_AD


!----------------------------------------------------------------------------------------------------------------------------------
!> This routine forms the dU^{ED}/du^{BD} and dU^{ED}/du^{AD} blocks (ED row) of dUdu. (i.e., how do changes in the AD and BD inputs affect the ED inputs?)
SUBROUTINE Linear_ED_InputSolve_du( p_FAST, y_FAST, SrvD, u_ED, y_ED, y_AD, u_AD, BD, HD, SD, MAPp, MD, MeshMapData, dUdu, ErrStat, ErrMsg )

   TYPE(FAST_ParameterType),       INTENT(IN   )  :: p_FAST         !< Glue-code simulation parameters
   TYPE(FAST_OutputFileType),      INTENT(IN   )  :: y_FAST         !< Glue-code output parameters (for linearization)
   type(ServoDyn_Data),            intent(in   )  :: SrvD           !< SrvD parameters
   TYPE(ED_InputType),             INTENT(INOUT)  :: u_ED           !< ED Inputs at t
   TYPE(ED_OutputType),            INTENT(IN   )  :: y_ED           !< ElastoDyn outputs (need translation displacement on meshes for loads mapping)
   TYPE(AD_OutputType),            INTENT(IN   )  :: y_AD           !< AeroDyn outputs
   TYPE(AD_InputType),             INTENT(INOUT)  :: u_AD           !< AD inputs (for AD-ED load linerization)
   TYPE(BeamDyn_Data),             INTENT(INOUT)  :: BD             !< BD data at t
   TYPE(HydroDyn_Data),            INTENT(INOUT)  :: HD             !< HD data at t
   TYPE(SubDyn_Data),              INTENT(INOUT)  :: SD             !< SD data at t
   TYPE(MAP_Data),                 INTENT(INOUT)  :: MAPp           !< MAP data at t
   TYPE(MoorDyn_Data),             INTENT(INOUT)  :: MD             !< MD data at t
   TYPE(FAST_ModuleMapType),       INTENT(INOUT)  :: MeshMapData    !< Data for mapping between modules
   REAL(R8Ki),                     INTENT(INOUT)  :: dUdu(:,:)      !< Jacobian matrix of which we are computing the dU^(ED)/du^(AD) block
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat        !< Error status
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg         !< Error message
   
      ! local variables
   INTEGER(IntKi)                                 :: j              ! Loops through StC instances
   INTEGER(IntKi)                                 :: K              ! Loops through blades
   INTEGER(IntKi)                                 :: SrvD_Start     ! starting index of dUdu (column) where SrvD StC load is
   INTEGER(IntKi)                                 :: BD_Start       ! starting index of dUdu (column) where BD root motion inputs are located
   INTEGER(IntKi)                                 :: AD_Start_Bl    ! starting index of dUdu (column) where AD blade motion inputs are located
   INTEGER(IntKi)                                 :: ED_Start_mt    ! starting index of dUdu (row) where ED blade/tower or hub moment inputs are located
   INTEGER(IntKi)                                 :: HD_Start       ! starting index of dUdu (column) where HD motion inputs are located
   INTEGER(IntKi)                                 :: SD_Start       ! starting index of dUdu (column) where SD TP motion inputs are located
   INTEGER(IntKi)                                 :: MAP_Start      ! starting index of dUdu (column) where MAP fairlead motion inputs are located
   INTEGER(IntKi)                                 :: MD_Start       ! starting index of dUdu (column) where MD fairlead motion inputs are located
   INTEGER(IntKi)                                 :: ErrStat2       ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                           :: ErrMsg2        ! temporary Error message if ErrStat /= ErrID_None
   
   CHARACTER(*), PARAMETER                        :: RoutineName = 'Linear_ED_InputSolve_du'
   
   
      ! Initialize error status
      
   ErrStat = ErrID_None
   ErrMsg = ""

   !..........
   ! dU^{ED}/du^{SrvD}
   !..........
   if (p_FAST%CompServo == MODULE_SrvD) then
      ED_Start_mt = y_FAST%Lin%Modules(MODULE_ED)%Instance(1)%LinStartIndx(LIN_INPUT_COL)
      !--------------------
      ! Blade (BD or ED)
      if ( p_FAST%CompElast == Module_ED ) then
         if ( allocated(SrvD%y%BStCLoadMesh) ) then
            do j=1,size(SrvD%y%BStCLoadMesh,2)
               do K = 1,SIZE(u_ED%BladePtLoads,1) ! Loop through all blades (p_ED%NumBl)
                  if (SrvD%y%BStCLoadMesh(K,j)%Committed) then
                     CALL Linearize_Point_to_Point( SrvD%y%BStCLoadMesh(k,j), u_ED%BladePtLoads(k), MeshMapData%BStC_P_2_ED_P_B(k,j), ErrStat2, ErrMsg2, SrvD%Input(1)%BStCMotionMesh(k,j), y_ED%BladeLn2Mesh(k) )
                        call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

                     ED_Start_mt = ED_Start_mt + u_ED%BladePtLoads(k)%NNodes*3         ! 3 forces at each node (we're going to start at the moments since the M_us matrix is for moments...)
                     SrvD_Start = y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%LinStartIndx(LIN_INPUT_COL) - 1 + SrvD%p%Jac_Idx_BStC_u(1,k,j)
                     ! SrvD is source in the mapping, so we want M_{uSm} (moments)
                     if (allocated(MeshMapData%BStC_P_2_ED_P_B(k,j)%dM%m_us )) then
                        call SetBlockMatrix( dUdu, MeshMapData%BStC_P_2_ED_P_B(k,j)%dM%m_us, ED_Start_mt, SrvD_Start )
                     end if
                     ! get starting index of next blade
                     ED_Start_mt = ED_Start_mt + u_ED%BladePtLoads(k)%NNodes* 3
                  endif
               enddo
            enddo
         endif
      endif
      !--------------------
      ! Nacelle (ED only)
      if ( allocated(SrvD%y%NStCLoadMesh) ) then
         do j = 1,size(SrvD%y%NStCLoadMesh)
            if (SrvD%y%NStCLoadMesh(j)%Committed) then
               call Linearize_Point_to_Point( SrvD%y%NStCLoadMesh(j), u_ED%NacelleLoads, MeshMapData%NStC_P_2_ED_P_N(j), ErrStat2, ErrMsg2, SrvD%Input(1)%NStCMotionMesh(j), y_ED%NacelleMotion )
                  call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

               ED_Start_mt = Indx_u_ED_Nacelle_Start(u_ED, y_FAST) &
                             + u_ED%NacelleLoads%NNodes      * 3             ! 3 forces at the nacelle (so we start at the moments)
               SrvD_Start = y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%LinStartIndx(LIN_INPUT_COL) - 1 + SrvD%p%Jac_Idx_NStC_u(1,j)
               ! SrvD is source in the mapping, so we want M_{uSm} (moments)
               if (allocated(MeshMapData%NStC_P_2_ED_P_N(j)%dM%m_us )) then
                  call SetBlockMatrix( dUdu, MeshMapData%NStC_P_2_ED_P_N(j)%dM%m_us, ED_Start_mt, SrvD_Start )
               end if
            endif
         enddo
      endif
      !--------------------
      ! Tower (ED only)
      if ( allocated(SrvD%y%TStCLoadMesh) ) then
         do j = 1,size(SrvD%y%TStCLoadMesh)
            if (SrvD%y%TStCLoadMesh(j)%Committed) then
               call Linearize_Point_to_Point( SrvD%y%TStCLoadMesh(j), u_ED%TowerPtLoads, MeshMapData%TStC_P_2_ED_P_T(j), ErrStat2, ErrMsg2, SrvD%Input(1)%TStCMotionMesh(j), y_ED%TowerLn2Mesh )
                  call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

               ED_Start_mt = Indx_u_ED_Tower_Start(u_ED, y_FAST) &
                             + u_ED%TowerPtLoads%NNodes      * 3             ! 3 forces at the nacelle (so we start at the moments)
               SrvD_Start = y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%LinStartIndx(LIN_INPUT_COL) - 1 + SrvD%p%Jac_Idx_TStC_u(1,j)
               ! SrvD is source in the mapping, so we want M_{uSm} (moments)
               if (allocated(MeshMapData%TStC_P_2_ED_P_T(j)%dM%m_us )) then
                  call SetBlockMatrix( dUdu, MeshMapData%TStC_P_2_ED_P_T(j)%dM%m_us, ED_Start_mt, SrvD_Start )
               endif
            endif
         enddo
      endif
      !--------------------
      ! Substructure (SD or ED)
      if (p_FAST%CompSub /= MODULE_SD) then
         if ( allocated(SrvD%y%SStCLoadMesh) ) then
            do j=1,size(SrvD%y%SStCLoadMesh)
               if (SrvD%y%SStCLoadMesh(j)%Committed) then
                  call Linearize_Point_to_Point( SrvD%y%SStCLoadMesh(j), u_ED%PlatformPtMesh, MeshMapData%SStC_P_P_2_SubStructure(j), ErrStat2, ErrMsg2, SrvD%Input(1)%SStCMotionMesh(j), y_ED%PlatformPtMesh )
                     call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

                  ED_Start_mt = Indx_u_ED_Platform_Start(u_ED, y_FAST) &
                                + u_ED%PlatformPtMesh%NNodes      * 3             ! 3 forces at the nacelle (so we start at the moments)
                  SrvD_Start = y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%LinStartIndx(LIN_INPUT_COL) - 1 + SrvD%p%Jac_Idx_SStC_u(1,j)
                  ! SrvD is source in the mapping, so we want M_{uSm} (moments)
                  if (allocated(MeshMapData%SStC_P_P_2_SubStructure(j)%dM%m_us )) then
                     call SetBlockMatrix( dUdu, MeshMapData%SStC_P_P_2_SubStructure(j)%dM%m_us, ED_Start_mt, SrvD_Start )
                  endif
               endif
            enddo
         endif
      endif
   endif


   !..........
   ! dU^{ED}/du^{AD}
   !..........
   IF ( p_FAST%CompAero == Module_AD ) THEN
   
         ! ED inputs on blade from AeroDyn
      IF (p_FAST%CompElast == Module_ED) THEN
         
         ED_Start_mt = y_FAST%Lin%Modules(MODULE_ED)%Instance(1)%LinStartIndx(LIN_INPUT_COL)
         
         DO K = 1,SIZE(u_ED%BladePtLoads,1) ! Loop through all blades (p_ED%NumBl)
            ED_Start_mt = ED_Start_mt + u_ED%BladePtLoads(k)%NNodes*3 ! skip the forces on this blade
            AD_Start_Bl = Indx_u_AD_Blade_Start(u_AD, y_FAST, k) 
            
            CALL Linearize_Line2_to_Point( y_AD%rotors(1)%BladeLoad(k), u_ED%BladePtLoads(k), MeshMapData%AD_L_2_BDED_B(k), ErrStat2, ErrMsg2, u_AD%rotors(1)%BladeMotion(k), y_ED%BladeLn2Mesh(k) )
               CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            
               ! AD is source in the mapping, so we want M_{uSm}               
            if (allocated(MeshMapData%AD_L_2_BDED_B(k)%dM%m_us )) then
               call SetBlockMatrix( dUdu, MeshMapData%AD_L_2_BDED_B(k)%dM%m_us, ED_Start_mt, AD_Start_Bl )
            end if
            
               ! get starting index of next blade
            ED_Start_mt = ED_Start_mt + u_ED%BladePtLoads(k)%NNodes* 3  ! skip the moments on this blade
               
         END DO

      END IF
      
      ! ED inputs on tower from AD:
      
      IF ( y_AD%rotors(1)%TowerLoad%Committed ) THEN
         ED_Start_mt = Indx_u_ED_Tower_Start(u_ED, y_FAST) &
                       + u_ED%TowerPtLoads%NNodes   * 3             ! 3 forces at each node (we're going to start at the moments)
         
         CALL Linearize_Line2_to_Point( y_AD%rotors(1)%TowerLoad, u_ED%TowerPtLoads, MeshMapData%AD_L_2_ED_P_T, ErrStat2, ErrMsg2, u_AD%rotors(1)%TowerMotion, y_ED%TowerLn2Mesh )
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
         ED_Start_mt = Indx_u_ED_Hub_Start(u_ED, y_FAST) &
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
   
   ED_Start_mt = Indx_u_ED_Platform_Start(u_ED, y_FAST) &
                       + u_ED%PlatformPtMesh%NNodes * 3         ! 3 forces at each node (we're going to start at the moments)
   
   if ( p_FAST%CompSub == Module_SD ) then
      !..........
      ! dU^{ED}/du^{SD}
      !..........
      ! Transfer SD load outputs to ED PlatformPtMesh input:
         ! we're mapping loads, so we also need the sibling meshes' displacements:
      SD_Start = Indx_u_SD_TPMesh_Start(SD%Input(1), y_FAST)
            
      call Linearize_Point_to_Point( SD%y%Y1Mesh, u_ED%PlatformPtMesh, MeshMapData%SD_TP_2_ED_P, ErrStat2, ErrMsg2, SD%Input(1)%TPMesh, y_ED%PlatformPtMesh) !SD%Input(1)%TPMesh and y_ED%PlatformPtMesh contain the displaced positions for load calculations
         call SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)
               
         ! SD is source in the mapping, so we want M_{uSm}
      if (allocated(MeshMapData%SD_TP_2_ED_P%dM%m_us )) then
         call SetBlockMatrix( dUdu, MeshMapData%SD_TP_2_ED_P%dM%m_us, ED_Start_mt, SD_Start )
      end if
      
   else if ( p_FAST%CompSub == Module_None ) then
      !..........
      ! dU^{ED}/du^{HD}
      !..........
   
      if ( p_FAST%CompHydro == Module_HD ) then ! HydroDyn-{ElastoDyn} 
         
            ! we're just going to assume u_ED%PlatformPtMesh is committed
         if ( HD%y%Morison%Mesh%Committed  ) then ! meshes for floating
               
            
               ! Transfer HD load outputs to ED PlatformPtMesh input:
               ! we're mapping loads, so we also need the sibling meshes' displacements:
            HD_Start = Indx_u_HD_Morison_Start(HD%Input(1), y_FAST)
            
            call Linearize_Point_to_Point( HD%y%Morison%Mesh, u_ED%PlatformPtMesh, MeshMapData%HD_M_P_2_SubStructure, ErrStat2, ErrMsg2, HD%Input(1)%Morison%Mesh, y_ED%PlatformPtMesh) !HD%Input(1)%Morison and y_ED%PlatformPtMesh contain the displaced positions for load calculations
               call SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)
               
               ! HD is source in the mapping, so we want M_{uSm}
            if (allocated(MeshMapData%HD_M_P_2_SubStructure%dM%m_us )) then
               call SetBlockMatrix( dUdu, MeshMapData%HD_M_P_2_SubStructure%dM%m_us, ED_Start_mt, HD_Start )
            end if
            
         end if   
         if ( HD%y%WAMITMesh%Committed  ) then ! meshes for floating
            
               ! Transfer HD load outputs to ED PlatformPtMesh input:
               ! we're mapping loads, so we also need the sibling meshes' displacements:
            HD_Start = Indx_u_HD_WAMIT_Start(HD%Input(1), y_FAST)
            
            call Linearize_Point_to_Point( HD%y%WAMITMesh, u_ED%PlatformPtMesh, MeshMapData%HD_W_P_2_SubStructure, ErrStat2, ErrMsg2, HD%Input(1)%WAMITMesh, y_ED%PlatformPtMesh) !HD%Input(1)%WAMITMesh and y_ED%PlatformPtMesh contain the displaced positions for load calculations
               call SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)
               
               ! HD is source in the mapping, so we want M_{uSm}
            if (allocated(MeshMapData%HD_W_P_2_SubStructure%dM%m_us )) then
               call SetBlockMatrix( dUdu, MeshMapData%HD_W_P_2_SubStructure%dM%m_us, ED_Start_mt, HD_Start )
            end if
            
         end if
      end if
   
      !..........
      ! dU^{ED}/du^{MAP}
      !..........
      
      if ( p_FAST%CompMooring == Module_MAP ) then 

         ED_Start_mt = Indx_u_ED_Platform_Start(u_ED, y_FAST) &
                        + u_ED%PlatformPtMesh%NNodes * 3         ! 3 forces at each node (we're going to start at the moments)
      
            ! Transfer MAP loads to ED PlatformPtmesh input:
            ! we're mapping loads, so we also need the sibling meshes' displacements:

         MAP_Start = y_FAST%Lin%Modules(MODULE_MAP)%Instance(1)%LinStartIndx(LIN_INPUT_COL)
      
         ! NOTE: Assumes at least one MAP Fairlead point    
      
         CALL Linearize_Point_to_Point( MAPp%y%ptFairleadLoad, u_ED%PlatformPtMesh, MeshMapData%Mooring_2_Structure, ErrStat2, ErrMsg2, MAPp%Input(1)%PtFairDisplacement, y_ED%PlatformPtMesh) !MAPp%Input(1)%ptFairleadLoad and y_ED%PlatformPtMesh contain the displaced positions for load calculations
            CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)
               
            ! MAP is source in the mapping, so we want M_{uSm}
         if (allocated(MeshMapData%Mooring_2_Structure%dM%m_us )) then
            call SetBlockMatrix( dUdu, MeshMapData%Mooring_2_Structure%dM%m_us, ED_Start_mt, MAP_Start )
         end if
   
      !..........
      ! dU^{ED}/du^{MD}
      !..........
      else if ( p_FAST%CompMooring == Module_MD ) then 

         ED_Start_mt = Indx_u_ED_Platform_Start(u_ED, y_FAST) &
                        + u_ED%PlatformPtMesh%NNodes * 3         ! 3 forces at each node (we're going to start at the moments)
      
            ! Transfer MD loads to ED PlatformPtmesh input:
            ! we're mapping loads, so we also need the sibling meshes' displacements:

         MD_Start = y_FAST%Lin%Modules(Module_MD)%Instance(1)%LinStartIndx(LIN_INPUT_COL)
      
         ! NOTE: Assumes at least one coupled MD object
      
         CALL Linearize_Point_to_Point( MD%y%CoupledLoads(1), u_ED%PlatformPtMesh, MeshMapData%Mooring_2_Structure, ErrStat2, ErrMsg2, MD%Input(1)%CoupledKinematics(1), y_ED%PlatformPtMesh)
            CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)
               
            ! HD is source in the mapping, so we want M_{uSm}
         if (allocated(MeshMapData%Mooring_2_Structure%dM%m_us )) then
            call SetBlockMatrix( dUdu, MeshMapData%Mooring_2_Structure%dM%m_us, ED_Start_mt, MD_Start )
         end if
   
      end if
      
   end if    
END SUBROUTINE Linear_ED_InputSolve_du


!----------------------------------------------------------------------------------------------------------------------------------
!> This routine forms the dU^{SD}/du^{SrvD}, dU^{SD}/du^{HD}, dU^{SD}/du^{SD}, and dU^{SD}/du^{MAP} blocks (SD row) of dUdu. (i.e., how do changes in SrvD, HD, SD, and MAP inputs affect the SD inputs?)
SUBROUTINE Linear_SD_InputSolve_du( p_FAST, y_FAST, SrvD, u_SD, y_SD, y_ED, HD, MAPp, MD, MeshMapData, dUdu, ErrStat, ErrMsg )

   TYPE(FAST_ParameterType),       INTENT(IN   )  :: p_FAST         !< Glue-code simulation parameters
   TYPE(FAST_OutputFileType),      INTENT(IN   )  :: y_FAST         !< Glue-code output parameters (for linearization)
   type(ServoDyn_Data),            intent(in   )  :: SrvD           !< SrvD parameters
   TYPE(SD_InputType),             INTENT(INOUT)  :: u_SD           !< SD Inputs at t
   TYPE(SD_OutputType),            INTENT(IN   )  :: y_SD           !< SubDyn outputs (need translation displacement on meshes for loads mapping)
   TYPE(ED_OutputType),            INTENT(IN   )  :: y_ED           !< ElastoDyn outputs
   TYPE(HydroDyn_Data),            INTENT(INOUT)  :: HD             !< HD data at t
   TYPE(MAP_Data),                 INTENT(INOUT)  :: MAPp           !< MAP data at t
   TYPE(MoorDyn_Data),             INTENT(INOUT)  :: MD             !< MD data at t
   TYPE(FAST_ModuleMapType),       INTENT(INOUT)  :: MeshMapData    !< Data for mapping between modules
   REAL(R8Ki),                     INTENT(INOUT)  :: dUdu(:,:)      !< Jacobian matrix of which we are computing the dU^(SD)/du^(AD) block
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat        !< Error status
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg         !< Error message
   
      ! local variables
   INTEGER(IntKi)                                 :: j, SrvD_Start
   INTEGER(IntKi)                                 :: HD_Start
   INTEGER(IntKi)                                 :: MAP_Start, MD_Start
   INTEGER(IntKi)                                 :: SD_Start, SD_Start_td, SD_Start_tr
   INTEGER(IntKi)                                 :: ErrStat2       ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                           :: ErrMsg2        ! temporary Error message if ErrStat /= ErrID_None
   
   CHARACTER(*), PARAMETER                        :: RoutineName = 'Linear_SD_InputSolve_du'
   
   
      ! Initialize error status
      
   ErrStat = ErrID_None
   ErrMsg = ""

   IF ( p_FAST%CompSub == Module_SD ) THEN ! see routine U_ED_SD_HD_BD_Orca_Residual() in SolveOption1
   

   !..........
   ! dU^{SD}/du^{SrvD}
   !..........
   if (p_FAST%CompServo == MODULE_SrvD) then
      !--------------------
      ! Substructure (SD or ED)
      if ( allocated(SrvD%y%SStCLoadMesh) ) then
         SD_Start = Indx_u_SD_LMesh_Start(u_SD, y_FAST) &
                              + u_SD%LMesh%NNodes * 3         ! 3 forces at each node (we're going to start at the moments)
         do j=1,size(SrvD%y%SStCLoadMesh)
            if (SrvD%y%SStCLoadMesh(j)%Committed) then
               call Linearize_Point_to_Point( SrvD%y%SStCLoadMesh(j), u_SD%LMesh, MeshMapData%SStC_P_P_2_SubStructure(j), ErrStat2, ErrMsg2, SrvD%Input(1)%SStCMotionMesh(j), y_SD%Y3Mesh )
                  call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

               SrvD_Start = y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%LinStartIndx(LIN_INPUT_COL) - 1 + SrvD%p%Jac_Idx_SStC_u(1,j)
               ! SrvD is source in the mapping, so we want M_{uSm} (moments)
               if (allocated(MeshMapData%SStC_P_P_2_SubStructure(j)%dM%m_us )) then
                  call SetBlockMatrix( dUdu, MeshMapData%SStC_P_P_2_SubStructure(j)%dM%m_us, SD_Start, SrvD_Start )
               endif
            endif
         enddo
      endif
   endif


   !..........
   ! dU^{SD}/du^{SD}
   !..........
   
   call Linearize_Point_to_Point( y_ED%PlatformPtMesh, u_SD%TPMesh, MeshMapData%ED_P_2_SD_TP, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         
   ! SD is destination in the mapping, so we want M_{tv_uD} and M_{ta_uD}
            
   SD_Start_td = y_FAST%Lin%Modules(MODULE_SD)%Instance(1)%LinStartIndx(LIN_INPUT_COL)  
   SD_Start_tr = SD_Start_td + u_SD%TPMesh%NNodes * 6 ! skip 2 fields (TranslationDisp and Orientation) with 3 components before translational velocity field      
         
      ! translational velocity:
   if (allocated(MeshMapData%ED_P_2_SD_TP%dM%tv_uD )) then             
      call SetBlockMatrix( dUdu, MeshMapData%ED_P_2_SD_TP%dM%tv_ud, SD_Start_tr, SD_Start_td )
   end if
         
      ! translational acceleration:
   SD_Start_tr = SD_Start_tr + u_SD%TPMesh%NNodes * 6 ! skip 2 fields ( TranslationVel and RotationVel)
   if (allocated(MeshMapData%ED_P_2_SD_TP%dM%ta_uD )) then            
      call SetBlockMatrix( dUdu, MeshMapData%ED_P_2_SD_TP%dM%ta_ud, SD_Start_tr, SD_Start_td )
   end if

   
  
   !..........
   ! dU^{SD}/du^{HD}
   !..........
   
    ! we're just going to assume u_SD%LMesh is committed
   SD_Start = Indx_u_SD_LMesh_Start(u_SD, y_FAST) &
                        + u_SD%LMesh%NNodes * 3         ! 3 forces at each node (we're going to start at the moments)
   
   if ( p_FAST%CompHydro == Module_HD ) then ! HydroDyn-{ElastoDyn or SubDyn}
           
         
         ! Transfer HD load outputs to SD LMesh input:
         if ( HD%y%Morison%Mesh%Committed  ) then ! meshes for floating
             
            ! dU^{SD}/du^{HD}  
            
               ! we're mapping loads, so we also need the sibling meshes' displacements:
            HD_Start = Indx_u_HD_Morison_Start(HD%Input(1), y_FAST)
            
            call Linearize_Point_to_Point( HD%y%Morison%Mesh, u_SD%LMesh, MeshMapData%HD_M_P_2_SubStructure, ErrStat2, ErrMsg2, HD%Input(1)%Morison%Mesh, y_SD%Y2Mesh) !HD%Input(1)%Mesh and y_ED%PlatformPtMesh contain the displaced positions for load calculations
               call SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)
               
               ! HD is source in the mapping, so we want M_{uSm}
            if (allocated(MeshMapData%HD_M_P_2_SubStructure%dM%m_us )) then
               call SetBlockMatrix( dUdu, MeshMapData%HD_M_P_2_SubStructure%dM%m_us, SD_Start, HD_Start )
            end if
            
        
            
            
         end if      
         if ( HD%y%WAMITMesh%Committed  ) then ! meshes for floating
             
               ! Transfer HD load outputs to ED PlatformPtMesh input:
               ! we're mapping loads, so we also need the sibling meshes' displacements:
            HD_Start = Indx_u_HD_WAMIT_Start(HD%Input(1), y_FAST)
            
            call Linearize_Point_to_Point( HD%y%WAMITMesh, u_SD%LMesh, MeshMapData%HD_W_P_2_SubStructure, ErrStat2, ErrMsg2, HD%Input(1)%WAMITMesh, y_SD%Y2Mesh) !HD%Input(1)%Mesh and y_ED%PlatformPtMesh contain the displaced positions for load calculations
               call SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)
               
               ! HD is source in the mapping, so we want M_{uSm}
            if (allocated(MeshMapData%HD_W_P_2_SubStructure%dM%m_us )) then
               call SetBlockMatrix( dUdu, MeshMapData%HD_W_P_2_SubStructure%dM%m_us, SD_Start, HD_Start )
            end if
            
           
         end if
   end if
   
   !..........
   ! dU^{SD}/du^{MAP}
   !..........

   if ( p_FAST%CompMooring == Module_MAP ) then 

         ! Transfer MAP loads to ED PlatformPtmesh input:
         ! we're mapping loads, so we also need the sibling meshes' displacements:

      MAP_Start = y_FAST%Lin%Modules(MODULE_MAP)%Instance(1)%LinStartIndx(LIN_INPUT_COL)
   
         MAP_Start = y_FAST%Lin%Modules(MODULE_MAP)%Instance(1)%LinStartIndx(LIN_INPUT_COL)
      
         ! NOTE: Assumes at least one MAP Fairlead point    
     
         CALL Linearize_Point_to_Point( MAPp%y%ptFairleadLoad, u_SD%LMesh, MeshMapData%Mooring_2_Structure, ErrStat2, ErrMsg2, MAPp%Input(1)%PtFairDisplacement, y_SD%Y3Mesh) !MAPp%Input(1)%ptFairleadLoad and y_SD%Y3Mesh contain the displaced positions for load calculations
            CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)
               
            ! SD is source in the mapping, so we want M_{uSm}
         if (allocated(MeshMapData%Mooring_2_Structure%dM%m_us )) then
            call SetBlockMatrix( dUdu, MeshMapData%Mooring_2_Structure%dM%m_us, SD_Start, MAP_Start )
         end if
   
   !..........
   ! dU^{SD}/du^{MD}
   !..........
   else if ( p_FAST%CompMooring == Module_MD ) then 

         ! Transfer MD loads to ED PlatformPtmesh input:
         ! we're mapping loads, so we also need the sibling meshes' displacements:

      MD_Start = y_FAST%Lin%Modules(Module_MD)%Instance(1)%LinStartIndx(LIN_INPUT_COL)
   
      ! NOTE: Assumes at least one coupled MD object
  
      CALL Linearize_Point_to_Point( MD%y%CoupledLoads(1), u_SD%LMesh, MeshMapData%Mooring_2_Structure, ErrStat2, ErrMsg2, MD%Input(1)%CoupledKinematics(1), y_SD%Y3Mesh)
         CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)
            
         ! SD is source in the mapping, so we want M_{uSm}
      if (allocated(MeshMapData%Mooring_2_Structure%dM%m_us )) then
         call SetBlockMatrix( dUdu, MeshMapData%Mooring_2_Structure%dM%m_us, SD_Start, MD_Start )
      end if

   end if
      
   END IF   
END SUBROUTINE Linear_SD_InputSolve_du


!----------------------------------------------------------------------------------------------------------------------------------
!> This routine forms the dU^{SD}/dy^{SrvD}, dU^{SD}/dy^{HD} and dU^{SD}/dy^{SD} blocks (SD row) of dUdu. (i.e., how do changes in SrvD, HD, and SD inputs affect the SD inputs?)
SUBROUTINE Linear_SD_InputSolve_dy( p_FAST, y_FAST, SrvD, u_SD, y_SD, y_ED, HD, MAPp, MD, MeshMapData, dUdy, ErrStat, ErrMsg )

   TYPE(FAST_ParameterType),       INTENT(IN   )  :: p_FAST         !< Glue-code simulation parameters
   TYPE(FAST_OutputFileType),      INTENT(IN   )  :: y_FAST         !< Glue-code output parameters (for linearization)
   type(ServoDyn_Data),            intent(in   )  :: SrvD           !< SrvD parameters
   TYPE(SD_InputType),             INTENT(INOUT)  :: u_SD           !< SD Inputs at t
   TYPE(SD_OutputType),            INTENT(IN   )  :: y_SD           !< SubDyn outputs (need translation displacement on meshes for loads mapping)
   TYPE(ED_OutputType),            INTENT(IN   )  :: y_ED           !< ElastoDyn outputs
   TYPE(HydroDyn_Data),            INTENT(INOUT)  :: HD             !< HD data at t
   TYPE(MAP_Data),                 INTENT(INOUT)  :: MAPp           !< MAP data at t
   TYPE(MoorDyn_Data),             INTENT(INOUT)  :: MD             !< MD data at t
   TYPE(FAST_ModuleMapType),       INTENT(INOUT)  :: MeshMapData    !< Data for mapping between modules
   REAL(R8Ki),                     INTENT(INOUT)  :: dUdy(:,:)      !< Jacobian matrix of which we are computing the dU^(SD)/dy^(SD) block
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat        !< Error status
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg         !< Error message
   
      ! local variables
   INTEGER(IntKi)                                 :: j, SrvD_Out_Start, SD_Start, SD_Out_Start, HD_Start, HD_Out_Start, ED_Out_Start, MAP_Out_Start, MD_Out_Start
   INTEGER(IntKi)                                 :: MAP_Start, MD_Start
!   INTEGER(IntKi)                                 :: ErrStat2       ! temporary Error status of the operation
!   CHARACTER(ErrMsgLen)                           :: ErrMsg2        ! temporary Error message if ErrStat /= ErrID_None
   
   CHARACTER(*), PARAMETER                        :: RoutineName = 'Linear_SD_InputSolve_du'
   
   
      ! Initialize error status
      
   ErrStat = ErrID_None
   ErrMsg = ""
   if ( p_FAST%CompSub /= Module_SD ) return
   
   !..........
   ! dU^{SD}/dy^{SrvD}
   !..........
   if (p_FAST%CompServo == MODULE_SrvD) then
      !--------------------
      ! Substructure (SD or ED)
      if ( allocated(SrvD%y%SStCLoadMesh) ) then
         SD_Start   = Indx_u_SD_LMesh_Start(u_SD, y_FAST)   ! start of u_SD%LMesh%Force field
         do j=1,size(SrvD%y%SStCLoadMesh)
            if (SrvD%y%SStCLoadMesh(j)%Committed) then
               SrvD_Out_Start = y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%LinStartIndx(LIN_INPUT_COL) - 1 + SrvD%p%Jac_Idx_SStC_y(1,j)
               call Assemble_dUdy_Loads(SrvD%y%SStCLoadMesh(j), u_SD%LMesh, MeshMapData%SStC_P_P_2_SubStructure(j), SD_Start, SrvD_Out_Start, dUdy)
            endif
         enddo
      endif
   endif

   !..........
   ! dU^{SD}/dy^{ED}
   !..........
  
   !!! ! This linearization was done in forming dUdu (see Linear_SD_InputSolve_du()), so we don't need to re-calculate these matrices 
   !!! ! while forming dUdy, too.
   !!!call Linearize_Point_to_Line2( y_ED%PlatformPtMesh, u_SD%TPMesh, MeshMapData%ED_P_2_SD_TP, ErrStat2, ErrMsg2 )
      
   SD_Start     = Indx_u_SD_TPMesh_Start(u_SD, y_FAST)  ! start of u_SD%MTPMesh%TranslationDisp field     
   ED_Out_Start = Indx_y_ED_Platform_Start(y_ED, y_FAST) ! start of y_ED%PlatformPtMesh%TranslationDisp field
   call Assemble_dUdy_Motions(y_ED%PlatformPtMesh, u_SD%TPMesh, MeshMapData%ED_P_2_SD_TP, SD_Start, ED_Out_Start, dUdy, .false.)
   
   !..........
   ! dU^{SD}/dy^{HD}
   !..........
   ! HD
   ! parts of dU^{SD}/dy^{HD} and dU^{SD}/dy^{SD}:
   if ( p_FAST%CompHydro == Module_HD ) then ! HydroDyn-SubDyn
      SD_Out_Start = Indx_y_SD_Y2Mesh_Start(y_SD, y_FAST) ! start of y_SD%Y2Mesh%TranslationDisp field
         ! we're just going to assume u_SD%LMesh is committed
      if ( HD%y%Morison%Mesh%Committed  ) then ! meshes for floating
         !!! ! This linearization was done in forming dUdu (see Linear_SD_InputSolve_du()), so we don't need to re-calculate these matrices 
         !!! ! while forming dUdy, too.
         ! call Linearize_Point_to_Point( HD%y%Morison%Mesh, u_SD%LMesh, MeshMapData%HD_M_P_2_SubStructure, ErrStat2, ErrMsg2, HD%Input(1)%Morison%Mesh, y_SD%Y2Mesh)
         HD_Out_Start = Indx_y_HD_Morison_Start(HD%y, y_FAST)
         SD_Start     = Indx_u_SD_LMesh_Start(u_SD, y_FAST) ! start of u_SD%LMesh%Force field
         call Assemble_dUdy_Loads(HD%y%Morison%Mesh, u_SD%LMesh, MeshMapData%HD_M_P_2_SubStructure, SD_Start, HD_Out_Start, dUdy)
         
            ! SD translation displacement-to-SD moment transfer (dU^{SD}/dy^{SD}):
         SD_Start = Indx_u_SD_LMesh_Start(u_SD, y_FAST) + u_SD%LMesh%NNodes*3   ! start of u_SD%LMesh%Moment field (skip the SD forces) 
         call SetBlockMatrix( dUdy, MeshMapData%HD_M_P_2_SubStructure%dM%m_uD, SD_Start, SD_Out_Start )
! maybe this should be SumBlockMatrix with future changes to linearized modules???            
      end if      
      if ( HD%y%WAMITMesh%Committed  ) then ! meshes for floating
         !!! ! This linearization was done in forming dUdu (see Linear_SD_InputSolve_du()), so we don't need to re-calculate these matrices 
         !!! ! while forming dUdy, too.
         ! call Linearize_Point_to_Point( HD%y%WAMITMesh, u_SD%LMesh, MeshMapData%HD_W_P_2_SubStructure, ErrStat2, ErrMsg2, HD%Input(1)%WAMITMesh, y_SD%Y2Mesh)
         HD_Out_Start = Indx_y_HD_WAMIT_Start(HD%y, y_FAST)
         SD_Start     = Indx_u_SD_LMesh_Start(u_SD, y_FAST) ! start of u_SD%LMesh%Force field
         call Assemble_dUdy_Loads(HD%y%WAMITMesh, u_SD%LMesh, MeshMapData%HD_W_P_2_SubStructure, SD_Start, HD_Out_Start, dUdy)
         
            ! SD translation displacement-to-SD moment transfer (dU^{SD}/dy^{SD}):
         SD_Start = Indx_u_SD_LMesh_Start(u_SD, y_FAST) + u_SD%LMesh%NNodes*3   ! start of u_SD%LMesh%Moment field (skip the SD forces)  
         call SumBlockMatrix( dUdy, MeshMapData%HD_W_P_2_SubStructure%dM%m_uD, SD_Start, SD_Out_Start )
! maybe this should be SumBlockMatrix with future changes to linearized modules???            
      end if
   end if
   
   !..........
   ! dU^{SD}/dy^{MAP}
   !..........
   if ( p_FAST%CompMooring == Module_MAP ) then
      if ( MAPp%y%ptFairleadLoad%Committed  ) then ! meshes for floating
         !!! ! This linearization was done in forming dUdu (see Linear_SD_InputSolve_du()), so we don't need to re-calculate these matrices 
         !!! ! while forming dUdy, too.
         !       CALL Linearize_Point_to_Point( MAPp%y%ptFairleadLoad, u_SD%LMesh, MeshMapData%Mooring_2_Structure, ErrStat2, ErrMsg2, MAPp%Input(1)%PtFairDisplacement, y_SD%Y3Mesh) !MAPp%Input(1)%ptFairleadLoad and y_ED%Y3Mesh contain the displaced positions for load calculations
         MAP_Out_Start = y_FAST%Lin%Modules(MODULE_MAP)%Instance(1)%LinStartIndx(LIN_OUTPUT_COL)
         SD_Start      = Indx_u_SD_LMesh_Start(u_SD, y_FAST) ! start of u_SD%LMesh%TranslationDisp field
         call Assemble_dUdy_Loads(MAPp%y%ptFairLeadLoad, u_SD%LMesh, MeshMapData%Mooring_2_Structure, SD_Start, MAP_Out_Start, dUdy)
      
         ! SD translation displacement-to-SD moment transfer (dU^{SD}/dy^{SD}):
         SD_Start = Indx_u_SD_LMesh_Start(u_SD, y_FAST) + u_SD%LMesh%NNodes*3   ! start of u_ED%LMesh%Moment field (skip the SD forces)
         SD_Out_Start = Indx_y_SD_Y3Mesh_Start(y_SD, y_FAST) ! start of y_SD%Y3Mesh%TranslationDisp field
         call SumBlockMatrix( dUdy, MeshMapData%Mooring_2_Structure%dM%m_uD, SD_Start, SD_Out_Start )
      end if     
      
   !..........
   ! dU^{SD}/dy^{MD}
   !..........
   else if ( p_FAST%CompMooring == Module_MD ) then
      if ( MD%y%CoupledLoads(1)%Committed  ) then ! meshes for floating
         !!! ! This linearization was done in forming dUdu (see Linear_SD_InputSolve_du()), so we don't need to re-calculate these matrices 
         !!! ! while forming dUdy, too.
         MD_Out_Start = y_FAST%Lin%Modules(Module_MD)%Instance(1)%LinStartIndx(LIN_OUTPUT_COL)
         SD_Start      = Indx_u_SD_LMesh_Start(u_SD, y_FAST) ! start of u_SD%LMesh%TranslationDisp field
         call Assemble_dUdy_Loads(MD%y%CoupledLoads(1), u_SD%LMesh, MeshMapData%Mooring_2_Structure, SD_Start, MD_Out_Start, dUdy)
      
         ! SD translation displacement-to-SD moment transfer (dU^{SD}/dy^{SD}):
         SD_Start = Indx_u_SD_LMesh_Start(u_SD, y_FAST) + u_SD%LMesh%NNodes*3   ! start of u_ED%LMesh%Moment field (skip the SD forces)
         SD_Out_Start = Indx_y_SD_Y3Mesh_Start(y_SD, y_FAST) ! start of y_SD%Y3Mesh%TranslationDisp field
         call SumBlockMatrix( dUdy, MeshMapData%Mooring_2_Structure%dM%m_uD, SD_Start, SD_Out_Start )
      end if     
   end if
END SUBROUTINE Linear_SD_InputSolve_dy  
   

!----------------------------------------------------------------------------------------------------------------------------------
!> This routine forms the dU^{BD}/du^{BD} and dU^{BD}/du^{AD} blocks (BD row) of dUdu. (i.e., how do changes in the AD and BD inputs 
!! affect the BD inputs?) This should be called only when p_FAST%CompElast == Module_BD.
SUBROUTINE Linear_BD_InputSolve_du( p_FAST, y_FAST, SrvD, y_ED, y_AD, u_AD, BD, MeshMapData, dUdu, ErrStat, ErrMsg )

   TYPE(FAST_ParameterType),       INTENT(IN   )  :: p_FAST         !< Glue-code simulation parameters
   TYPE(FAST_OutputFileType),      INTENT(IN   )  :: y_FAST         !< Glue-code output parameters (for linearization)
   type(ServoDyn_Data),            intent(in   )  :: SrvD           !< SrvD parameters
   TYPE(ED_OutputType),            INTENT(IN   )  :: y_ED           !< ElastoDyn outputs (need translation displacement on meshes for loads mapping)
   TYPE(AD_OutputType),            INTENT(IN   )  :: y_AD           !< AeroDyn outputs
   TYPE(AD_InputType),             INTENT(INOUT)  :: u_AD           !< AD inputs (for AD-ED load linerization)
   TYPE(BeamDyn_Data),             INTENT(INOUT)  :: BD             !< BD data at t

   TYPE(FAST_ModuleMapType),       INTENT(INOUT)  :: MeshMapData    !< Data for mapping between modules
   REAL(R8Ki),                     INTENT(INOUT)  :: dUdu(:,:)      !< Jacobian matrix of which we are computing the dU^(ED)/du^(AD) block
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat        !< Error status
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg         !< Error message
   
      ! local variables
   INTEGER(IntKi)                                 :: j              ! Loops through StC instances
   INTEGER(IntKi)                                 :: k              ! Loops through blades
   INTEGER(IntKi)                                 :: SrvD_Start     ! starting index of dUdu (row) where BD inputs are located
   INTEGER(IntKi)                                 :: BD_Start       ! starting index of dUdu (row) where BD inputs are located
   INTEGER(IntKi)                                 :: AD_Start       ! starting index of dUdu (column) where AD inputs are located
   INTEGER(IntKi)                                 :: ErrStat2       ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                           :: ErrMsg2        ! temporary Error message if ErrStat /= ErrID_None
   
   CHARACTER(*), PARAMETER                        :: RoutineName = 'Linear_BD_InputSolve_du'
   
   
      ! Initialize error status
      
   ErrStat = ErrID_None
   ErrMsg = ""

   !..........
   ! dU^{BD}/du^{SrvD}
   !..........
   if (p_FAST%CompServo == MODULE_SrvD) then
      !--------------------
      ! Blade (BD or ED)
      if ( allocated(SrvD%y%BStCLoadMesh) ) then
         do j=1,size(SrvD%y%BStCLoadMesh,2)
            do K = 1,p_FAST%nBeams ! Loop through all blades
               if (SrvD%y%BStCLoadMesh(K,j)%Committed) then
                  CALL Linearize_Point_to_Line2( SrvD%y%BStCLoadMesh(k,j), BD%Input(1,k)%DistrLoad, MeshMapData%BStC_P_2_BD_P_B(k,j), ErrStat2, ErrMsg2, SrvD%Input(1)%BStCMotionMesh(k,j), BD%y(k)%BldMotion )
                     call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

                  BD_Start = y_FAST%Lin%Modules(MODULE_BD)%Instance(k)%LinStartIndx(LIN_INPUT_COL) &
                           + BD%Input(1,k)%RootMotion%NNodes *18  & ! displacement, rotation, & acceleration fields for each node
                           + BD%Input(1,k)%PointLoad%NNodes  * 6  & ! force + moment fields for each node
                           + BD%Input(1,k)%DistrLoad%NNodes  * 3    ! force field for each node (start with moment field)
                  SrvD_Start = y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%LinStartIndx(LIN_INPUT_COL) - 1 + SrvD%p%Jac_Idx_BStC_u(1,k,j)
                  ! SrvD is source in the mapping, so we want M_{uSm} (moments)
                  if (allocated(MeshMapData%BStC_P_2_BD_P_B(k,j)%dM%m_us )) then
                     call SetBlockMatrix( dUdu, MeshMapData%BStC_P_2_BD_P_B(k,j)%dM%m_us, BD_Start, SrvD_Start )
                  end if
               endif
            enddo
         enddo
      endif
   endif

   !..........
   ! dU^{BD}/du^{AD}
   !..........
   IF ( p_FAST%CompAero == Module_AD ) THEN
   
      ! BD inputs on blade from AeroDyn
   
         
      if (p_FAST%BD_OutputSibling) then
            
         DO K = 1,p_FAST%nBeams ! Loop through all blades
            CALL Linearize_Line2_to_Line2( y_AD%rotors(1)%BladeLoad(k), BD%Input(1,k)%DistrLoad, MeshMapData%AD_L_2_BDED_B(k), ErrStat2, ErrMsg2, u_AD%rotors(1)%BladeMotion(k), BD%y(k)%BldMotion )
               CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         END DO
         
      else
      
         DO K = 1,p_FAST%nBeams ! Loop through all blades
            !linearization for dUdy will need some matrix multiplies because of the transfer (chain rule!), but we will perform individual linearization calculations here
            !!! need to transfer the BD output blade motions to nodes on a sibling of the BD blade motion mesh:
            CALL Linearize_Line2_to_Line2( BD%y(k)%BldMotion, MeshMapData%y_BD_BldMotion_4Loads(k), MeshMapData%BD_L_2_BD_L(k), ErrStat2, ErrMsg2 )
               CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

            CALL Linearize_Line2_to_Line2( y_AD%rotors(1)%BladeLoad(k), BD%Input(1,k)%DistrLoad, MeshMapData%AD_L_2_BDED_B(k), ErrStat2, ErrMsg2, u_AD%rotors(1)%BladeMotion(k), MeshMapData%y_BD_BldMotion_4Loads(k) )
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
   INTEGER(IntKi)                               :: AD_Start_ta    ! starting index of dUdu (column) where AD translation accelerations are located
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
   IF (u_AD%rotors(1)%TowerMotion%Committed) THEN

         CALL Linearize_Line2_to_Line2( y_ED%TowerLn2Mesh, u_AD%rotors(1)%TowerMotion, MeshMapData%ED_L_2_AD_L_T, ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//':u_AD%TowerMotion' )     

      
      !AD is the destination here, so we need tv_ud
      if (allocated( MeshMapData%ED_L_2_AD_L_T%dM%tv_ud)) then
         AD_Start_td = y_FAST%Lin%Modules(MODULE_AD)%Instance(1)%LinStartIndx(LIN_INPUT_COL)
         AD_Start_tv = AD_Start_td + u_AD%rotors(1)%TowerMotion%NNodes * 6 ! 2 fields (TranslationDisp and Orientation) with 3 components before translational velocity field      

         call SetBlockMatrix( dUdu, MeshMapData%ED_L_2_AD_L_T%dM%tv_ud, AD_Start_tv, AD_Start_td )
      end if
               
            
   END IF
   
      
      ! blades
   IF (p_FAST%CompElast == Module_ED ) THEN
      
      DO k=1,size(u_AD%rotors(1)%BladeMotion)
         CALL Linearize_Line2_to_Line2( y_ED%BladeLn2Mesh(k), u_AD%rotors(1)%BladeMotion(k), MeshMapData%BDED_L_2_AD_L_B(k), ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//':u_AD%BladeMotion('//trim(num2lstr(k))//')' )
      END DO
      
   ELSEIF (p_FAST%CompElast == Module_BD ) THEN
   
      DO k=1,size(u_AD%rotors(1)%BladeMotion)
         CALL Linearize_Line2_to_Line2( BD%y(k)%BldMotion, u_AD%rotors(1)%BladeMotion(k), MeshMapData%BDED_L_2_AD_L_B(k), ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//':u_AD%BladeMotion('//trim(num2lstr(k))//')' )
      END DO
         
   END IF
   
   
   
   DO k=1,size(u_AD%rotors(1)%BladeMotion)
      AD_Start_td = Indx_u_AD_Blade_Start(u_AD, y_FAST, k) ! index for u_AD%BladeMotion(k)%translationDisp field

         !AD is the destination here, so we need tv_ud
      if (allocated( MeshMapData%BDED_L_2_AD_L_B(k)%dM%tv_ud)) then
            ! index for u_AD%BladeMotion(k+1)%translationVel field
         AD_Start_tv = AD_Start_td + u_AD%rotors(1)%BladeMotion(k)%NNodes * 6 ! 2 fields (TranslationDisp and Orientation) with 3 components before translational velocity field

         call SetBlockMatrix( dUdu, MeshMapData%BDED_L_2_AD_L_B(k)%dM%tv_ud, AD_Start_tv, AD_Start_td )
      end if
         
      if (allocated( MeshMapData%BDED_L_2_AD_L_B(k)%dM%ta_ud)) then
         AD_Start_ta = AD_Start_td + u_AD%rotors(1)%BladeMotion(k)%NNodes * 12 ! 4 fields (TranslationDisp, Orientation, TranslationVel, and RotationVel) with 3 components before translational velocity field
         
         call SetBlockMatrix( dUdu, MeshMapData%BDED_L_2_AD_L_B(k)%dM%ta_ud, AD_Start_ta, AD_Start_td )
      end if

   END DO
END SUBROUTINE Linear_AD_InputSolve_du



!----------------------------------------------------------------------------------------------------------------------------------
!> This routine forms the dU^{SrvD}/du^{SrvD} block (SrvD row) of dUdu.
!! (i.e., how do changes in the SrvD inputs affect the SrvD inputs?)
SUBROUTINE Linear_SrvD_InputSolve_du( p_FAST, y_FAST, p_SrvD, u_SrvD, y_ED, BD, SD, MeshMapData, dUdu, ErrStat, ErrMsg )
   type(FAST_ParameterType),     intent(in   )  :: p_FAST         !< Glue-code simulation parameters
   type(FAST_OutputFileType),    intent(in   )  :: y_FAST         !< Glue-code output parameters (for linearization)
   type(SrvD_ParameterType),     intent(in   )  :: p_SrvD         !< SrvD parameters
   type(SrvD_InputType),         intent(inout)  :: u_SrvD         !< SrvD Inputs at t
   type(ED_OutputType),          intent(in   )  :: y_ED           !< ElastoDyn outputs (need translation displacement on meshes for loads mapping)
   type(BeamDyn_Data),           intent(in   )  :: BD             !< BD data at t
   type(SubDyn_Data),            intent(in   )  :: SD             !< SD data at t
   type(FAST_ModuleMapType),     intent(inout)  :: MeshMapData    !< Data for mapping between modules
   real(R8Ki),                   intent(inout)  :: dUdu(:,:)      !< Jacobian matrix of which we are computing the dU^(ED)/du^(AD) block
   integer(IntKi),               intent(  out)  :: ErrStat        !< Error status
   character(*),                 intent(  out)  :: ErrMsg         !< Error message
   
      ! local variables
   integer(IntKi)                               :: i,j,k          ! Generic counters 
   INTEGER(IntKi)                               :: SrvD_Start     ! starting index of dUdu (column) where the StC motion inputs are located
   integer(IntKi)                               :: ErrStat2       ! temporary Error status of the operation
   character(ErrMsgLen)                         :: ErrMsg2        ! temporary Error message if ErrStat /= ErrID_None
   character(*), parameter                      :: RoutineName = 'Linear_SrvD_InputSolve_du'
   
      ! Initialize error status
   ErrStat = ErrID_None
   ErrMsg = ""

   !--------------------
   ! dU^{SrvD}/du^{SrvD}
   !--------------------
   ! Blade StrucCtrl
   if ( p_FAST%CompElast == Module_ED ) then
      if ( ALLOCATED(u_SrvD%BStCMotionMesh) ) then
         do j=1,size(u_SrvD%BStCMotionMesh,2)
            do K = 1,size(y_ED%BladeLn2Mesh)
               if (u_SrvD%BStCMotionMesh(K,j)%Committed) then
                  CALL Linearize_Line2_to_Point( y_ED%BladeLn2Mesh(K), u_SrvD%BStCMotionMesh(K,j), MeshMapData%ED_L_2_BStC_P_B(K,j), ErrStat2, ErrMsg2 )
                     call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

                  ! SrvD is destination in the mapping, so we want M_{tv_uD} and M_{ta_uD}
                  ! translational velocity:
                  if (allocated(MeshMapData%ED_L_2_BStC_P_B(K,j)%dM%tv_uD )) then
                     SrvD_Start = y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%LinStartIndx(LIN_INPUT_COL) - 1 + (p_SrvD%Jac_Idx_BStC_u(1,k,j) + 6) ! skip translational displacement and orientation fields
                     call SetBlockMatrix( dUdu, MeshMapData%ED_L_2_BStC_P_B(K,j)%dM%tv_uD, SrvD_Start, y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%LinStartIndx(LIN_INPUT_COL) )
                  end if

                  ! translational acceleration:
                  if (allocated(MeshMapData%ED_L_2_BStC_P_B(K,j)%dM%ta_uD )) then
                     SrvD_Start = y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%LinStartIndx(LIN_INPUT_COL) - 1 + (p_SrvD%Jac_Idx_BStC_u(1,k,j) + 12) ! skip translational displacement and orientation fields
                     call SetBlockMatrix( dUdu, MeshMapData%ED_L_2_BStC_P_B(K,j)%dM%ta_uD, SrvD_Start, y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%LinStartIndx(LIN_INPUT_COL) )
                  end if
                endif
             enddo
          enddo
       endif
    elseif ( p_FAST%CompElast == Module_BD ) then
       if ( ALLOCATED(u_SrvD%BStCMotionMesh) ) then
          do j=1,size(u_SrvD%BStCMotionMesh,2)
             do K = 1,p_FAST%nBeams
                if (u_SrvD%BStCMotionMesh(K,j)%Committed) then
                   CALL Linearize_Line2_to_Point( BD%y(k)%BldMotion, u_SrvD%BStCMotionMesh(K,j), MeshMapData%BD_L_2_BStC_P_B(K,j), ErrStat2, ErrMsg2 )
                      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

                  ! SrvD is destination in the mapping, so we want M_{tv_uD} and M_{ta_uD}
                  ! translational velocity:
                  if (allocated(MeshMapData%BD_L_2_BStC_P_B(K,j)%dM%tv_uD )) then
                     SrvD_Start = y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%LinStartIndx(LIN_INPUT_COL) - 1 + (p_SrvD%Jac_Idx_BStC_u(1,k,j) + 6) ! skip translational displacement and orientation fields
                     call SetBlockMatrix( dUdu, MeshMapData%BD_L_2_BStC_P_B(K,j)%dM%tv_uD, SrvD_Start, y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%LinStartIndx(LIN_INPUT_COL) )
                  end if

                  ! translational acceleration:
                  if (allocated(MeshMapData%BD_L_2_BStC_P_B(K,j)%dM%ta_uD )) then
                     SrvD_Start = y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%LinStartIndx(LIN_INPUT_COL) - 1 + (p_SrvD%Jac_Idx_BStC_u(1,k,j) + 12) ! skip translational displacement and orientation fields
                     call SetBlockMatrix( dUdu, MeshMapData%BD_L_2_BStC_P_B(K,j)%dM%ta_uD, SrvD_Start, y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%LinStartIndx(LIN_INPUT_COL) )
                  end if
                endif
             enddo
          enddo
       endif
    endif
   !--------------------
   ! Nacelle (ED only)
   if ( ALLOCATED(u_SrvD%NStCMotionMesh) ) then
      do j = 1,size(u_SrvD%NStCMotionMesh)
         if (u_SrvD%NStCMotionMesh(j)%Committed) then
            call Linearize_Point_to_Point( y_ED%NacelleMotion, u_SrvD%NStCMotionMesh(j), MeshMapData%ED_P_2_NStC_P_N(j), ErrStat2, ErrMsg2 )
               call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

            ! SrvD is destination in the mapping, so we want M_{tv_uD} and M_{ta_uD}
            ! translational velocity:
            if (allocated(MeshMapData%ED_P_2_NStC_P_N(j)%dM%tv_uD )) then
               SrvD_Start = y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%LinStartIndx(LIN_INPUT_COL) - 1 + (p_SrvD%Jac_Idx_NStC_u(1,j) + 6) ! skip translational displacement and orientation fields
               call SetBlockMatrix( dUdu, MeshMapData%ED_P_2_NStC_P_N(j)%dM%tv_uD, SrvD_Start, y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%LinStartIndx(LIN_INPUT_COL) )
            end if

            ! translational acceleration:
            if (allocated(MeshMapData%ED_P_2_NStC_P_N(j)%dM%ta_uD )) then
               SrvD_Start = y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%LinStartIndx(LIN_INPUT_COL) - 1 + (p_SrvD%Jac_Idx_NStC_u(1,j) + 12) ! skip translational displacement and orientation fields
               call SetBlockMatrix( dUdu, MeshMapData%ED_P_2_NStC_P_N(j)%dM%ta_uD, SrvD_Start, y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%LinStartIndx(LIN_INPUT_COL) )
            end if
         endif
      enddo
   endif
   !--------------------
   ! Tower
   if ( ALLOCATED(u_SrvD%TStCMotionMesh) ) then
      do j = 1,size(u_SrvD%TStCMotionMesh)
         if (u_SrvD%TStCMotionMesh(j)%Committed) then
            call Linearize_Line2_to_Point( y_ED%TowerLn2Mesh, u_SrvD%TStCMotionMesh(j), MeshMapData%ED_L_2_TStC_P_T(j), ErrStat2, ErrMsg2 )
               call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

            ! SrvD is destination in the mapping, so we want M_{tv_uD} and M_{ta_uD}
            ! translational velocity:
            if (allocated(MeshMapData%ED_L_2_TStC_P_T(j)%dM%tv_uD )) then
               SrvD_Start = y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%LinStartIndx(LIN_INPUT_COL) - 1 + (p_SrvD%Jac_Idx_TStC_u(1,j) + 6) ! skip translational displacement and orientation fields
               call SetBlockMatrix( dUdu, MeshMapData%ED_L_2_TStC_P_T(j)%dM%tv_uD, SrvD_Start, y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%LinStartIndx(LIN_INPUT_COL) )
            end if

            ! translational acceleration:
            if (allocated(MeshMapData%ED_L_2_TStC_P_T(j)%dM%ta_uD )) then
               SrvD_Start = y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%LinStartIndx(LIN_INPUT_COL) - 1 + (p_SrvD%Jac_Idx_TStC_u(1,j) + 12) ! skip translational displacement and orientation fields
               call SetBlockMatrix( dUdu, MeshMapData%ED_L_2_TStC_P_T(j)%dM%ta_uD, SrvD_Start, y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%LinStartIndx(LIN_INPUT_COL) )
            end if
         endif
      enddo
   endif
   !--------------------
   ! Substructure (SD or ED)
   if (p_FAST%CompSub /= MODULE_SD) then
      if ( ALLOCATED(u_SrvD%SStCMotionMesh) ) then
         do j=1,size(u_SrvD%SStCMotionMesh)
            if (u_SrvD%SStCMotionMesh(j)%Committed) then
               CALL Linearize_Point_to_Point( y_ED%PlatformPtMesh, u_SrvD%SStCMotionMesh(j), MeshMapData%Substructure_2_SStC_P_P(j), ErrStat2, ErrMsg2 )
                  call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

               ! SrvD is destination in the mapping, so we want M_{tv_uD} and M_{ta_uD}
               ! translational velocity:
               if (allocated(MeshMapData%Substructure_2_SStC_P_P(j)%dM%tv_uD )) then
                  SrvD_Start = y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%LinStartIndx(LIN_INPUT_COL) - 1 + (p_SrvD%Jac_Idx_SStC_u(1,j) + 6) ! skip translational displacement and orientation fields
                  call SetBlockMatrix( dUdu, MeshMapData%Substructure_2_SStC_P_P(j)%dM%tv_uD, SrvD_Start, y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%LinStartIndx(LIN_INPUT_COL) )
               end if
   
               ! translational acceleration:
               if (allocated(MeshMapData%Substructure_2_SStC_P_P(j)%dM%ta_uD )) then
                  SrvD_Start = y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%LinStartIndx(LIN_INPUT_COL) - 1 + (p_SrvD%Jac_Idx_SStC_u(1,j) + 12) ! skip translational displacement and orientation fields
                  call SetBlockMatrix( dUdu, MeshMapData%Substructure_2_SStC_P_P(j)%dM%ta_uD, SrvD_Start, y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%LinStartIndx(LIN_INPUT_COL) )
               end if
            endif
         enddo
      endif
   else
      if ( ALLOCATED(u_SrvD%SStCMotionMesh) ) then
         do j=1,size(u_SrvD%SStCMotionMesh)
            IF (u_SrvD%SStCMotionMesh(j)%Committed) then
               CALL Linearize_Point_to_Point( SD%y%y3Mesh, u_SrvD%SStCMotionMesh(j), MeshMapData%SubStructure_2_SStC_P_P(j), ErrStat2, ErrMsg2 )
                  call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

               ! SrvD is destination in the mapping, so we want M_{tv_uD} and M_{ta_uD}
               ! translational velocity:
               if (allocated(MeshMapData%SubStructure_2_SStC_P_P(j)%dM%tv_uD )) then
                  SrvD_Start = y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%LinStartIndx(LIN_INPUT_COL) - 1 + (p_SrvD%Jac_Idx_SStC_u(1,j) + 6) ! skip translational displacement and orientation fields
                  call SetBlockMatrix( dUdu, MeshMapData%SubStructure_2_SStC_P_P(j)%dM%tv_uD, SrvD_Start, y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%LinStartIndx(LIN_INPUT_COL) )
               end if
   
               ! translational acceleration:
               if (allocated(MeshMapData%SubStructure_2_SStC_P_P(j)%dM%ta_uD )) then
                  SrvD_Start = y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%LinStartIndx(LIN_INPUT_COL) - 1 + (p_SrvD%Jac_Idx_SStC_u(1,j) + 12) ! skip translational displacement and orientation fields
                  call SetBlockMatrix( dUdu, MeshMapData%SubStructure_2_SStC_P_P(j)%dM%ta_uD, SrvD_Start, y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%LinStartIndx(LIN_INPUT_COL) )
               end if
            endif
         enddo
      endif
   endif
END SUBROUTINE Linear_SrvD_InputSolve_du



!----------------------------------------------------------------------------------------------------------------------------------
!> This routine forms the dU^{SrvD}/dy^{ED}, dU^{SrvD}/dy^{BD}, dU^{SrvD}/dy^{SD} block of dUdy.
!! (i.e., how do changes in the ED, SD, BD outputs affect the SrvD inputs?)
!! NOTE: Linearze_Point_to_Point routines done in Linear_SrvD_InputSolve_du
SUBROUTINE Linear_SrvD_InputSolve_dy( p_FAST, y_FAST, p_SrvD, u_SrvD, y_ED, BD, y_SD, MeshMapData, dUdy, ErrStat, ErrMsg )
!..................................................................................................................................
   type(FAST_ParameterType),        intent(in   )  :: p_FAST           !< Glue-code simulation parameters
   type(FAST_OutputFileType),       intent(in   )  :: y_FAST           !< Output variables for the glue code
   type(SrvD_ParameterType),        intent(in   )  :: p_SrvD           !< SrvD parameters (holds indices for jacobian entries for each StC)
   type(SrvD_InputType),            intent(inout)  :: u_SrvD           !< SrvD Inputs at t
   type(ED_OutputType),             intent(in   )  :: y_ED             !< ElastoDyn outputs (need translation displacement on meshes for loads mapping)
   type(BeamDyn_Data),              intent(in   )  :: BD               !< BeamDyn   data
   type(SD_OutputType),             intent(in   )  :: y_SD             !< SubDyn    outputs (need translation displacement on meshes for loads mapping)
                                    
   type(FAST_ModuleMapType),        intent(inout)  :: MeshMapData      !< Data for mapping between modules
   real(R8Ki),                      intent(inout)  :: dUdy(:,:)        !< Jacobian matrix of which we are computing the dU^(ED)/du^(AD) block
   integer(IntKi),                  intent(  out)  :: ErrStat          !< Error status
   character(*),                    intent(  out)  :: ErrMsg           !< Error message
                                    
   integer(IntKi)                                  :: i,j,k            ! loop counters
   integer(intKi)                                  :: ED_Start_Yaw     !< starting index of dUdy (column) where ED Yaw/YawRate/HSS_Spd outputs are located (just before WriteOutput)
   integer(IntKi)                                  :: SrvD_Start, ED_Out_Start, BD_Out_Start, SD_Out_Start
   character(*), parameter                         :: RoutineName = 'Linear_SrvD_InputSolve_dy' 
   
      ! Initialize error status
   ErrStat = ErrID_None
   ErrMsg = ""

   !--------------------
   ! dU^{SrvD}/dy^{ED}
   !--------------------
   ED_Start_Yaw = Indx_y_Yaw_Start(y_FAST, Module_ED) ! start of ED where Yaw, YawRate, HSS_Spd occur (right before WriteOutputs)
   do i=1,size(SrvD_Indx_Y_BlPitchCom)    ! first 3 columns
      dUdy(y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%LinStartIndx(LIN_INPUT_COL) + SrvD_Indx_Y_BlPitchCom(i) - 1, ED_Start_Yaw + i - 1) = -1.0_ReKi
   end do

   !----------------------------------------
   ! Structural controls
   !----------------------------------------
   ! Blade
   if ( p_FAST%CompElast == Module_ED ) then
      !--------------------
      ! dU^{SrvD}/dy^{ED}
      !--------------------
      if ( ALLOCATED(u_SrvD%BStCMotionMesh) ) then
         do j=1,size(u_SrvD%BStCMotionMesh,2)
            do K = 1,size(y_ED%BladeLn2Mesh)
               if (u_SrvD%BStCMotionMesh(K,j)%Committed) then
                  SrvD_Start   = y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%LinStartIndx(LIN_INPUT_COL) - 1 + (p_SrvD%Jac_Idx_BStC_u(1,k,j))
                  ED_Out_Start = Indx_y_ED_Blade_Start(y_ED, y_FAST, k)    ! start of %TranslationDisp field
                  call Assemble_dUdy_Motions( y_ED%BladeLn2Mesh(K), u_SrvD%BStCMotionMesh(K,j), MeshMapData%ED_L_2_BStC_P_B(K,j), SrvD_Start, ED_Out_Start, dUdy, .false.)
               endif
            enddo
         enddo
      endif
   elseif ( p_FAST%CompElast == Module_BD ) then
      !--------------------
      ! dU^{SrvD}/dy^{BD}
      !--------------------
      if ( ALLOCATED(u_SrvD%BStCMotionMesh) ) then
         do j=1,size(u_SrvD%BStCMotionMesh,2)
            do K = 1,p_FAST%nBeams
               if (u_SrvD%BStCMotionMesh(K,j)%Committed) then
                  SrvD_Start   = y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%LinStartIndx(LIN_INPUT_COL) - 1 + (p_SrvD%Jac_Idx_BStC_u(1,k,j))
                  BD_Out_Start = y_FAST%Lin%Modules(MODULE_BD)%Instance(k)%LinStartIndx(LIN_OUTPUT_COL)    ! start of %TranslationDisp field
                  call Assemble_dUdy_Motions( BD%y(k)%BldMotion, u_SrvD%BStCMotionMesh(K,j), MeshMapData%BD_L_2_BStC_P_B(K,j), SrvD_Start, BD_Out_Start, dUdy, .false.)
               endif
            enddo
         enddo
      endif
   endif

   !--------------------
   ! Nacelle -- dU^{SrvD}/dy^{ED}
   !--------------------
   if ( ALLOCATED(u_SrvD%NStCMotionMesh) ) then
      do j = 1,size(u_SrvD%NStCMotionMesh)
         if (u_SrvD%NStCMotionMesh(j)%Committed) then
            SrvD_Start   = y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%LinStartIndx(LIN_INPUT_COL) - 1 + (p_SrvD%Jac_Idx_NStC_u(1,j))
            ED_Out_Start = Indx_y_ED_Nacelle_Start(y_ED, y_FAST)    ! start of %TranslationDisp field
            call Assemble_dUdy_Motions( y_ED%NacelleMotion, u_SrvD%NStCMotionMesh(j), MeshMapData%ED_P_2_NStC_P_N(j), SrvD_Start, ED_Out_Start, dUdy, .false.)
         endif
      enddo
   endif

   !--------------------
   ! Tower -- dU^{SrvD}/dy^{ED}
   !--------------------
   if ( ALLOCATED(u_SrvD%TStCMotionMesh) ) then
      do j = 1,size(u_SrvD%TStCMotionMesh)
         if (u_SrvD%TStCMotionMesh(j)%Committed) then
            SrvD_Start   = y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%LinStartIndx(LIN_INPUT_COL) - 1 + (p_SrvD%Jac_Idx_TStC_u(1,j))
            ED_Out_Start = Indx_y_ED_Tower_Start(y_ED, y_FAST)    ! start of %TranslationDisp field
            call Assemble_dUdy_Motions( y_ED%TowerLn2Mesh, u_SrvD%TStCMotionMesh(j), MeshMapData%ED_L_2_TStC_P_T(j), SrvD_Start, ED_Out_Start, dUdy, .false.)
         endif
      enddo
   endif

   !--------------------
   ! Substructure (SD or ED)
   !--------------------
   if (p_FAST%CompSub /= MODULE_SD) then
      !--------------------
      ! dU^{SrvD}/dy^{ED}
      !--------------------
      if ( ALLOCATED(u_SrvD%SStCMotionMesh) ) then
         do j=1,size(u_SrvD%SStCMotionMesh)
            if (u_SrvD%SStCMotionMesh(j)%Committed) then
               SrvD_Start   = y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%LinStartIndx(LIN_INPUT_COL) - 1 + (p_SrvD%Jac_Idx_SStC_u(1,j))
               ED_Out_Start = Indx_y_ED_Platform_Start(y_ED, y_FAST) ! start of %TranslationDisp field
               call Assemble_dUdy_Motions( y_ED%PlatformPtMesh, u_SrvD%SStCMotionMesh(j), MeshMapData%Substructure_2_SStC_P_P(j), SrvD_Start, ED_Out_Start, dUdy, .false.)
            endif
         enddo
      endif
   else
      !--------------------
      ! dU^{SrvD}/dy^{SD}
      !--------------------
      if ( ALLOCATED(u_SrvD%SStCMotionMesh) ) then
         do j=1,size(u_SrvD%SStCMotionMesh)
            if (u_SrvD%SStCMotionMesh(j)%Committed) then
               SrvD_Start   = y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%LinStartIndx(LIN_INPUT_COL) - 1 + (p_SrvD%Jac_Idx_SStC_u(1,j))
               SD_Out_Start = Indx_y_SD_Y3Mesh_Start(y_SD, y_FAST)   ! start of %TranslationDisp field
               call Assemble_dUdy_Motions( y_SD%y3Mesh, u_SrvD%SStCMotionMesh(j), MeshMapData%SubStructure_2_SStC_P_P(j), SrvD_Start, SD_Out_Start, dUdy, .false.)
            endif
         enddo
      endif
   endif
END SUBROUTINE Linear_SrvD_InputSolve_dy



!----------------------------------------------------------------------------------------------------------------------------------
!> This routine forms the dU^{ED}/dy^{SrvD}, dU^{ED}/dy^{ED}, dU^{ED}/dy^{BD},  dU^{ED}/dy^{AD}, dU^{ED}/dy^{HD}, and dU^{ED}/dy^{MAP}
!! blocks of dUdy. (i.e., how do changes in the SrvD, ED, BD, AD, HD, and MAP outputs effect the ED inputs?)
SUBROUTINE Linear_ED_InputSolve_dy( p_FAST, y_FAST, SrvD, u_ED, y_ED, y_AD, u_AD, BD, HD, SD, MAPp, MD, MeshMapData, dUdy, ErrStat, ErrMsg )

   TYPE(FAST_ParameterType),       INTENT(IN   )  :: p_FAST           !< Glue-code simulation parameters
   TYPE(FAST_OutputFileType),      INTENT(IN   )  :: y_FAST           !< FAST output file data (for linearization)
   type(ServoDyn_Data),            intent(in   )  :: SrvD             !< SrvD parameters
   TYPE(ED_InputType),             INTENT(INOUT)  :: u_ED             !< ED Inputs at t
   TYPE(ED_OutputType),            INTENT(IN   )  :: y_ED             !< ElastoDyn outputs (need translation displacement on meshes for loads mapping)
   TYPE(AD_OutputType),            INTENT(IN   )  :: y_AD             !< AeroDyn outputs
   TYPE(AD_InputType),             INTENT(INOUT)  :: u_AD             !< AD inputs (for AD-ED load linerization)
   TYPE(BeamDyn_Data),             INTENT(INOUT)  :: BD               !< BD data at t
   TYPE(HydroDyn_Data),            INTENT(INOUT)  :: HD               !< HD data at t
   TYPE(SubDyn_Data),              INTENT(INOUT)  :: SD               !< SD data at t
   TYPE(MAP_Data),                 INTENT(INOUT)  :: MAPp             !< MAP data at t
   TYPE(MoorDyn_Data),             INTENT(INOUT)  :: MD               !< MD data at t
                                                                      
   TYPE(FAST_ModuleMapType),       INTENT(INOUT)  :: MeshMapData      !< Data for mapping between modules
   REAL(R8Ki),                     INTENT(INOUT)  :: dUdy(:,:)        !< Jacobian matrix of which we are computing the dU^(ED)/du^(AD) block
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat          !< Error status
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg           !< Error message
   
      ! local variables
   INTEGER(IntKi)                                 :: i                ! rows/columns
   INTEGER(IntKi)                                 :: j                ! Loops through StC instance 
   INTEGER(IntKi)                                 :: K                ! Loops through blades
   INTEGER(IntKi)                                 :: SrvD_Out_Start   ! starting index of dUdy (column) where the StC motion inputs are located
   INTEGER(IntKi)                                 :: AD_Out_Start     ! starting index of dUdy (column) where particular AD fields are located
   INTEGER(IntKi)                                 :: BD_Out_Start     ! starting index of dUdy (column) where particular BD fields are located
   INTEGER(IntKi)                                 :: ED_Start         ! starting index of dUdy (row) where ED input fields are located
   INTEGER(IntKi)                                 :: ED_Out_Start     ! starting index of dUdy (column) where ED output fields are located
   INTEGER(IntKi)                                 :: HD_Out_Start     ! starting index of dUdy (column) where HD output fields are located
   INTEGER(IntKi)                                 :: SD_Out_Start     ! starting index of dUdy (column) where SD output fields are located
   INTEGER(IntKi)                                 :: MAP_Out_Start    ! starting index of dUdy (column) where MAP output fields are located
   INTEGER(IntKi)                                 :: MD_Out_Start     ! starting index of dUdy (column) where MoorDyn output fields are located
   CHARACTER(*), PARAMETER                        :: RoutineName = 'Linear_ED_InputSolve_dy' 
   
   
      ! Initialize error status
      
   ErrStat = ErrID_None
   ErrMsg = ""

   !..........
   ! dU^{ED}/dy^{SrvD}
   !  ED inputs from ServoDyn outputs
   !..........
   IF ( p_FAST%CompServo == Module_SrvD ) THEN

         ! BlPitchCom, YawMom, GenTrq
      ED_Start = Indx_u_ED_BlPitchCom_Start(u_ED, y_FAST)
      do i=1,size(u_ED%BlPitchCom)+2 ! BlPitchCom, YawMom, GenTrq (NOT collective pitch)
         dUdy(ED_Start + i - 1, y_FAST%Lin%Modules(Module_SrvD)%Instance(1)%LinStartIndx(LIN_OUTPUT_COL) + i - 1) = -1.0_ReKi !SrvD_Indx_Y_BlPitchCom
      end do
      !--------------------
      ! Blade (BD or ED)
      if ( p_FAST%CompElast == Module_ED ) then
         if ( allocated(SrvD%y%BStCLoadMesh) ) then
            do j=1,size(SrvD%y%BStCLoadMesh,2)
               do K = 1,SIZE(u_ED%BladePtLoads,1) ! Loop through all blades (p_ED%NumBl)
                  if (SrvD%y%BStCLoadMesh(K,j)%Committed) then
                     ED_Start       = Indx_u_ED_Blade_Start(u_ED, y_FAST, k) ! start of u_ED%BladePtLoads(k)%Force field
                     SrvD_Out_Start = y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%LinStartIndx(LIN_INPUT_COL) - 1 + SrvD%p%Jac_Idx_BStC_y(1,k,j)
                     call Assemble_dUdy_Loads(SrvD%y%BStCLoadMesh(k,j), u_ED%BladePtLoads(k), MeshMapData%BStC_P_2_ED_P_B(k,j), ED_Start, SrvD_Out_Start, dUdy)

                     ! ED translation displacement-to-ED moment transfer (dU^{ED}/dy^{ED}):
                     ED_Start     = Indx_u_ED_Blade_Start(u_ED, y_FAST, k) + u_ED%BladePtLoads(k)%NNodes*3  ! start of u_ED%BladePtLoads(k)%Moment field (skip the ED forces)
                     ED_Out_Start = Indx_y_ED_Blade_Start(y_ED, y_FAST, k)                                  ! start of y_ED%BladeLn2Mesh(1)%TranslationDisp field
                     call SumBlockMatrix( dUdy, MeshMapData%BStC_P_2_ED_P_B(k,j)%dM%m_uD, ED_Start, ED_Out_Start )
                  endif
               enddo
            enddo
         endif
      endif
      !--------------------
      ! Nacelle (ED only)
      if ( allocated(SrvD%y%NStCLoadMesh) ) then
         do j = 1,size(SrvD%y%NStCLoadMesh)
            if (SrvD%y%NStCLoadMesh(j)%Committed) then
               ED_Start       = Indx_u_ED_Nacelle_Start(u_ED, y_FAST)
               SrvD_Out_Start = y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%LinStartIndx(LIN_INPUT_COL) - 1 + SrvD%p%Jac_Idx_NStC_y(1,j)
               call Assemble_dUdy_Loads(SrvD%y%NStCLoadMesh(j), u_ED%NacelleLoads, MeshMapData%NStC_P_2_ED_P_N(j), ED_Start, SrvD_Out_Start, dUdy)

               ! ED translation displacement-to-ED moment transfer (dU^{ED}/dy^{ED}):
               ED_Start     = Indx_u_ED_Nacelle_Start(u_ED, y_FAST) + u_ED%NacelleLoads%NNodes*3   ! start of u_ED%NacelleLoads%Moment field (skip the ED forces)
               ED_Out_Start = Indx_y_ED_Nacelle_Start(y_ED, y_FAST)                                ! start of y_ED%NacelleMotion%TranslationDisp field
               call SumBlockMatrix( dUdy, MeshMapData%NStC_P_2_ED_P_N(j)%dM%m_uD, ED_Start, ED_Out_Start )
            endif
         enddo
      endif
      !--------------------
      ! Tower (ED only)
      if ( allocated(SrvD%y%TStCLoadMesh) ) then
         do j = 1,size(SrvD%y%TStCLoadMesh)
            if (SrvD%y%TStCLoadMesh(j)%Committed) then
               ED_Start       = Indx_u_ED_Tower_Start(u_ED, y_FAST)
               SrvD_Out_Start = y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%LinStartIndx(LIN_INPUT_COL) - 1 + SrvD%p%Jac_Idx_TStC_y(1,j)
               call Assemble_dUdy_Loads(SrvD%y%TStCLoadMesh(j), u_ED%TowerPtLoads, MeshMapData%TStC_P_2_ED_P_T(j), ED_Start, SrvD_Out_Start, dUdy)

                  ! ED translation displacement-to-ED moment transfer (dU^{ED}/dy^{ED}):
               ED_Start     = Indx_u_ED_Tower_Start(u_ED, y_FAST) + u_ED%TowerPtLoads%NNodes*3  ! start of u_ED%TowerPtLoads%Moment field  [skip the ED forces to get to the moments]
               ED_Out_Start = Indx_y_ED_Tower_Start(y_ED, y_FAST)                               ! start of y_ED%TowerLn2Mesh%TranslationDisp field
               call SumBlockMatrix( dUdy, MeshMapData%TStC_P_2_ED_P_T(j)%dM%m_uD, ED_Start, ED_Out_Start )
            endif
         enddo
      endif
      !--------------------
      ! Substructure (SD or ED)
      if (p_FAST%CompSub /= MODULE_SD) then
         if ( allocated(SrvD%y%SStCLoadMesh) ) then
            do j=1,size(SrvD%y%SStCLoadMesh)
               if (SrvD%y%SStCLoadMesh(j)%Committed) then
                  ED_Start       = Indx_u_ED_Platform_Start(u_ED, y_FAST)
                  SrvD_Out_Start = y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%LinStartIndx(LIN_INPUT_COL) - 1 + SrvD%p%Jac_Idx_SStC_y(1,j)
                  call Assemble_dUdy_Loads(SrvD%y%SStCLoadMesh(j), u_ED%PlatformPtMesh, MeshMapData%SStC_P_P_2_SubStructure(j), ED_Start, SrvD_Out_Start, dUdy)

                     ! ED translation displacement-to-ED moment transfer (dU^{ED}/dy^{ED}):
                  ED_Start     = Indx_u_ED_Platform_Start(u_ED, y_FAST) + u_ED%PlatformPtMesh%NNodes*3   ! start of u_ED%PlatformPtMesh%Moment field (skip the ED forces)
                  ED_Out_Start = Indx_y_ED_Platform_Start(y_ED, y_FAST)                                  ! start of y_ED%PlatformPtMesh%TranslationDisp field
                  call SumBlockMatrix( dUdy, MeshMapData%HD_M_P_2_SubStructure%dM%m_uD, ED_Start, ED_Out_Start )
               endif
            enddo
         endif
      endif
   END IF

   ! parts of dU^{ED}/dy^{AD} and dU^{ED}/dy^{ED}:
   
      ! ElastoDyn inputs on blade from AeroDyn and ElastoDyn
   IF ( p_FAST%CompAero == Module_AD ) THEN

      IF (p_FAST%CompElast == Module_ED) THEN 
         AD_Out_Start = y_FAST%Lin%Modules(MODULE_AD)%Instance(1)%LinStartIndx(LIN_OUTPUT_COL) + y_AD%rotors(1)%TowerLoad%NNodes * 6    ! start of y_AD%BladeLoad(1)%Force field [2 fields (force, moment) with 3 components]
         
         DO K = 1,SIZE(u_ED%BladePtLoads,1) ! Loop through all blades (p_ED%NumBl)
            !!! ! This linearization was done in forming dUdu (see Linear_ED_InputSolve_du()), so we don't need to re-calculate these matrices 
            !!! ! while forming dUdy, too.
            !CALL Linearize_Line2_to_Point( y_AD%BladeLoad(k), u_ED%BladePtLoads(k), MeshMapData%AD_L_2_BDED_B(k), ErrStat2, ErrMsg2, u_AD%BladeMotion(k), y_ED%BladeLn2Mesh(k) )
               
               ! AD loads-to-ED loads transfer (dU^{ED}/dy^{AD}):
            ED_Start = Indx_u_ED_Blade_Start(u_ED, y_FAST, k) ! start of u_ED%BladePtLoads(k)%Force field
            call Assemble_dUdy_Loads(y_AD%rotors(1)%BladeLoad(k), u_ED%BladePtLoads(k), MeshMapData%AD_L_2_BDED_B(k), ED_Start, AD_Out_Start, dUdy)

               ! ED translation displacement-to-ED moment transfer (dU^{ED}/dy^{ED}):
            ED_Start = Indx_u_ED_Blade_Start(u_ED, y_FAST, k) + u_ED%BladePtLoads(k)%NNodes*3   ! start of u_ED%BladePtLoads(k)%Moment field (skip the ED forces)
            ED_Out_Start = Indx_y_ED_Blade_Start(y_ED, y_FAST, k) ! start of y_ED%BladeLn2Mesh(1)%TranslationDisp field
            call SumBlockMatrix( dUdy, MeshMapData%AD_L_2_BDED_B(k)%dM%m_uD, ED_Start, ED_Out_Start )

            AD_Out_Start = AD_Out_Start + y_AD%rotors(1)%BladeLoad(k)%NNodes*6        ! start of y_AD%BladeLoad(k+1)%Force field [skip 2 fields to forces on next blade]
         END DO
      END IF ! ED
      
      
      IF ( y_AD%rotors(1)%TowerLoad%Committed ) THEN
         !!! ! This linearization was done in forming dUdu (see Linear_ED_InputSolve_du()), so we don't need to re-calculate these matrices 
         !!! ! while forming dUdy, too.
         !CALL Linearize_Line2_to_Point( y_AD%TowerLoad, u_ED%TowerPtLoads, MeshMapData%AD_L_2_ED_P_T, ErrStat2, ErrMsg2, u_AD%TowerMotion, y_ED%TowerLn2Mesh )
            
            ! AD loads-to-ED loads transfer (dU^{ED}/dy^{AD}):
         ED_Start = Indx_u_ED_Tower_Start(u_ED, y_FAST) ! u_ED%TowerPtLoads%Force field
         AD_Out_Start = y_FAST%Lin%Modules(MODULE_AD)%Instance(1)%LinStartIndx(LIN_OUTPUT_COL) ! start of y_AD%Tower%Force
         call Assemble_dUdy_Loads(y_AD%rotors(1)%TowerLoad, u_ED%TowerPtLoads, MeshMapData%AD_L_2_ED_P_T, ED_Start, AD_Out_Start, dUdy)

            ! ED translation displacement-to-ED moment transfer (dU^{ED}/dy^{ED}):
         ED_Start = ED_Start + u_ED%TowerPtLoads%NNodes*3 ! start of u_ED%TowerPtLoads%Moment field  [skip the ED forces to get to the moments]
         ED_Out_Start  = Indx_y_ED_Tower_Start(y_ED, y_FAST) ! start of y_ED%TowerLn2Mesh%TranslationDisp field
         call SumBlockMatrix( dUdy, MeshMapData%AD_L_2_ED_P_T%dM%m_uD, ED_Start, ED_Out_Start )
            
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
         call SumBlockMatrix( dUdy, MeshMapData%BD_P_2_ED_P(k)%dM%m_ud, ED_Start, ED_Out_Start)
      END DO

   END IF

   if ( p_FAST%CompSub == Module_None ) then !This also occurs with ExtPtfm (though that's not linearized, yet)
      ! HD
      ! parts of dU^{ED}/dy^{HD} and dU^{ED}/dy^{ED}:
      if ( p_FAST%CompHydro == Module_HD ) then ! HydroDyn-{ElastoDyn or SubDyn}
            ED_Out_Start = Indx_y_ED_Platform_Start(y_ED, y_FAST) ! start of y_ED%PlatformPtMesh%TranslationDisp field
               ! we're just going to assume u_ED%PlatformPtMesh is committed
            if ( HD%y%Morison%Mesh%Committed  ) then ! meshes for floating
               !!! ! This linearization was done in forming dUdu (see Linear_ED_InputSolve_du()), so we don't need to re-calculate these matrices 
               !!! ! while forming dUdy, too.
               ! call Linearize_Point_to_Point( HD%y%Morison, u_ED%PlatformPtMesh, MeshMapData%HD_M_P_2_SubStructure, ErrStat2, ErrMsg2, HD%Input(1)%Morison, y_ED%PlatformPtMesh) !HD%Input(1)%Morison and y_ED%PlatformPtMesh contain the displaced positions for load calculations
               HD_Out_Start = Indx_y_HD_Morison_Start(HD%y, y_FAST)
               ED_Start     = Indx_u_ED_Platform_Start(u_ED, y_FAST) ! start of u_ED%PlatformPtMesh%Force field
               call Assemble_dUdy_Loads(HD%y%Morison%Mesh, u_ED%PlatformPtMesh, MeshMapData%HD_M_P_2_SubStructure, ED_Start, HD_Out_Start, dUdy)
         
                  ! ED translation displacement-to-ED moment transfer (dU^{ED}/dy^{ED}):
               ED_Start = Indx_u_ED_Platform_Start(u_ED, y_FAST) + u_ED%PlatformPtMesh%NNodes*3   ! start of u_ED%PlatformPtMesh%Moment field (skip the ED forces) 
               call SumBlockMatrix( dUdy, MeshMapData%HD_M_P_2_SubStructure%dM%m_uD, ED_Start, ED_Out_Start )
              
            end if      
            if ( HD%y%WAMITMesh%Committed  ) then ! meshes for floating
               !!! ! This linearization was done in forming dUdu (see Linear_ED_InputSolve_du()), so we don't need to re-calculate these matrices 
               !!! ! while forming dUdy, too.
               ! call Linearize_Point_to_Point( HD%y%WAMITMesh, u_ED%PlatformPtMesh, MeshMapData%HD_W_P_2_SubStructure, ErrStat2, ErrMsg2, HD%Input(1)%WAMITMesh, y_ED%PlatformPtMesh) !HD%Input(1)%WAMITMesh and y_ED%PlatformPtMesh contain the displaced positions for load calculations
               HD_Out_Start = Indx_y_HD_WAMIT_Start(HD%y, y_FAST)
               ED_Start     = Indx_u_ED_Platform_Start(u_ED, y_FAST) ! start of u_ED%PlatformPtMesh%Force field
               call Assemble_dUdy_Loads(HD%y%WAMITMesh, u_ED%PlatformPtMesh, MeshMapData%HD_W_P_2_SubStructure, ED_Start, HD_Out_Start, dUdy)
         
                  ! ED translation displacement-to-ED moment transfer (dU^{ED}/dy^{ED}):
               ED_Start = Indx_u_ED_Platform_Start(u_ED, y_FAST) + u_ED%PlatformPtMesh%NNodes*3   ! start of u_ED%PlatformPtMesh%Moment field (skip the ED forces)
               call SumBlockMatrix( dUdy, MeshMapData%HD_W_P_2_SubStructure%dM%m_uD, ED_Start, ED_Out_Start )
       
            end if

            
      end if
   
      ! MAP
      ! parts of dU^{ED}/dy^{MAP} and dU^{ED}/dy^{ED}:
      if ( p_FAST%CompMooring == Module_MAP ) then
         if ( MAPp%y%ptFairleadLoad%Committed  ) then ! meshes for floating
            !!! ! This linearization was done in forming dUdu (see Linear_ED_InputSolve_du()), so we don't need to re-calculate these matrices 
            !!! ! while forming dUdy, too.
            !       CALL Linearize_Point_to_Point( MAPp%y%ptFairleadLoad, u_ED%PlatformPtMesh, MeshMapData%Mooring_2_Structure, ErrStat2, ErrMsg2, MAPp%Input(1)%PtFairDisplacement, y_ED%PlatformPtMesh) !MAPp%Input(1)%ptFairleadLoad and y_ED%PlatformPtMesh contain the displaced positions for load calculations
            MAP_Out_Start = y_FAST%Lin%Modules(MODULE_MAP)%Instance(1)%LinStartIndx(LIN_OUTPUT_COL)
            ED_Start      = Indx_u_ED_Platform_Start(u_ED, y_FAST) ! start of u_ED%PlatformPtMesh%Force field
            call Assemble_dUdy_Loads(MAPp%y%ptFairLeadLoad, u_ED%PlatformPtMesh, MeshMapData%Mooring_2_Structure, ED_Start, MAP_Out_Start, dUdy)
      
            ! ED translation displacement-to-ED moment transfer (dU^{ED}/dy^{ED}):
            ED_Start = Indx_u_ED_Platform_Start(u_ED, y_FAST) + u_ED%PlatformPtMesh%NNodes*3   ! start of u_ED%PlatformPtMesh%Moment field (skip the ED forces)
            ED_Out_Start = Indx_y_ED_Platform_Start(y_ED, y_FAST) ! start of y_ED%PlatformPtMesh%TranslationDisp field
            call SumBlockMatrix( dUdy, MeshMapData%Mooring_2_Structure%dM%m_uD, ED_Start, ED_Out_Start )
         end if
      ! MoorDyn
      ! parts of dU^{ED}/dy^{MD} and dU^{ED}/dy^{ED}:
      else if ( p_FAST%CompMooring == Module_MD ) then
         if ( MD%y%CoupledLoads(1)%Committed  ) then ! meshes for floating
            !!! ! This linearization was done in forming dUdu (see Linear_ED_InputSolve_du()), so we don't need to re-calculate these matrices 
            !!! ! while forming dUdy, too.
            MD_Out_Start = y_FAST%Lin%Modules(Module_MD)%Instance(1)%LinStartIndx(LIN_OUTPUT_COL)
            ED_Start      = Indx_u_ED_Platform_Start(u_ED, y_FAST) ! start of u_ED%PlatformPtMesh%TranslationDisp field
            call Assemble_dUdy_Loads(MD%y%CoupledLoads(1), u_ED%PlatformPtMesh, MeshMapData%Mooring_2_Structure, ED_Start, MD_Out_Start, dUdy)
      
            ! ED translation displacement-to-ED moment transfer (dU^{ED}/dy^{ED}):
            ED_Start = Indx_u_ED_Platform_Start(u_ED, y_FAST) + u_ED%PlatformPtMesh%NNodes*3   ! start of u_ED%PlatformPtMesh%Moment field (skip the ED forces)
            ED_Out_Start = Indx_y_ED_Platform_Start(y_ED, y_FAST) ! start of y_ED%PlatformPtMesh%TranslationDisp field
            call SumBlockMatrix( dUdy, MeshMapData%Mooring_2_Structure%dM%m_uD, ED_Start, ED_Out_Start )
         end if
      end if
   else if ( p_FAST%CompSub == Module_SD ) then
      ! SubDyn
      ! parts of dU^{ED}/dy^{SD} and dU^{ED}/dy^{ED}:
      !!! ! This linearization was done in forming dUdu (see Linear_ED_InputSolve_du()), so we don't need to re-calculate these matrices 
      !!! ! while forming dUdy, too.
      !       CALL Linearize_Point_to_Point( SD%y%Y1Mesh, u_ED%PlatformPtMesh, MeshMapData%SD_TP_2_ED_P, ErrStat2, ErrMsg2, SD%Input(1)%TPMesh, y_ED%PlatformPtMesh) !SD%Input(1)%TPMesh and y_ED%PlatformPtMesh contain the displaced positions for load calculations
      SD_Out_Start = y_FAST%Lin%Modules(MODULE_SD)%Instance(1)%LinStartIndx(LIN_OUTPUT_COL)
      ED_Start      = Indx_u_ED_Platform_Start(u_ED, y_FAST) ! start of u_ED%PlatformPtMesh%TranslationDisp field
      call Assemble_dUdy_Loads(SD%y%Y1Mesh, u_ED%PlatformPtMesh, MeshMapData%SD_TP_2_ED_P, ED_Start, SD_Out_Start, dUdy)
           
      ! ED translation displacement-to-ED moment transfer (dU^{ED}/dy^{ED}):
      ED_Start = Indx_u_ED_Platform_Start(u_ED, y_FAST) + u_ED%PlatformPtMesh%NNodes*3   ! start of u_ED%PlatformPtMesh%Moment field (skip the ED forces)
      ED_Out_Start = Indx_y_ED_Platform_Start(y_ED, y_FAST) ! start of y_ED%PlatformPtMesh%TranslationDisp field
      call SetBlockMatrix( dUdy, MeshMapData%SD_TP_2_ED_P%dM%m_uD, ED_Start, ED_Out_Start )
      
      !Mooring gets set in the Linear_SD_InputSolve_ routines
   end if
END SUBROUTINE Linear_ED_InputSolve_dy
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine forms the dU^{BD}/dy^{ED}, dU^{BD}/dy^{BD}, and dU^{BD}/dy^{AD} blocks of dUdy. (i.e., how do 
!! changes in the ED, BD, and AD outputs effect the BD inputs?)
SUBROUTINE Linear_BD_InputSolve_dy( p_FAST, y_FAST, SrvD, u_ED, y_ED, y_AD, u_AD, BD, MeshMapData, dUdy, ErrStat, ErrMsg )

   TYPE(FAST_ParameterType),       INTENT(IN   )  :: p_FAST           !< Glue-code simulation parameters
   TYPE(FAST_OutputFileType),      INTENT(IN   )  :: y_FAST           !< FAST output file data (for linearization)
   type(ServoDyn_Data),            intent(in   )  :: SrvD             !< SrvD parameters
   TYPE(ED_InputType),             INTENT(INOUT)  :: u_ED             !< ED Inputs at t
   TYPE(ED_OutputType),            INTENT(IN   )  :: y_ED             !< ElastoDyn outputs (need translation displacement on meshes for loads mapping)
   TYPE(AD_OutputType),            INTENT(IN   )  :: y_AD             !< AeroDyn outputs
   TYPE(AD_InputType),             INTENT(INOUT)  :: u_AD             !< AD inputs (for AD-ED load linearization)
   TYPE(BeamDyn_Data),             INTENT(IN   )  :: BD               !< BD data at t
                                                                      
   TYPE(FAST_ModuleMapType),       INTENT(INOUT)  :: MeshMapData      !< Data for mapping between modules
   REAL(R8Ki),                     INTENT(INOUT)  :: dUdy(:,:)        !< Jacobian matrix of which we are computing the dU^(ED)/du^(AD) block
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat          !< Error status
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg           !< Error message
   
      ! local variables
   INTEGER(IntKi)                                 :: j                ! Loops through StC instance
   INTEGER(IntKi)                                 :: K                ! Loops through blades
   INTEGER(IntKi)                                 :: SrvD_Out_Start   ! starting index of dUdy (column) where particular SrvD fields are located
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

   !..........
   ! dU^{ED}/du^{SrvD}
   !..........
   if (p_FAST%CompServo == MODULE_SrvD) then
      !--------------------
      ! Blade (BD or ED)
      if ( allocated(SrvD%y%BStCLoadMesh) ) then
         do j=1,size(SrvD%y%BStCLoadMesh,2)
            do K = 1,p_FAST%nBeams ! Loop through all blades
               if (SrvD%y%BStCLoadMesh(K,j)%Committed) then
                  ! Start of DistrLoad nodes
                  BD_Start = y_FAST%Lin%Modules(MODULE_BD)%Instance(k)%LinStartIndx(LIN_INPUT_COL) &
                           + BD%Input(1,k)%RootMotion%NNodes *18  & ! displacement, rotation, & acceleration fields for each node
                           + BD%Input(1,k)%PointLoad%NNodes  * 6    ! force + moment fields for each node
                  SrvD_Out_Start = y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%LinStartIndx(LIN_INPUT_COL) - 1 + SrvD%p%Jac_Idx_BStC_y(1,k,j)
                  call Assemble_dUdy_Loads(SrvD%y%BStCLoadMesh(k,j), BD%Input(1,k)%DistrLoad, MeshMapData%BStC_P_2_BD_P_B(k,j), BD_Start, SrvD_Out_Start, dUdy)
               endif
            enddo
         enddo
      endif
   endif


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

      AD_Out_Start = y_FAST%Lin%Modules(MODULE_AD)%Instance(1)%LinStartIndx(LIN_OUTPUT_COL) + y_AD%rotors(1)%TowerLoad%NNodes * 6    ! start of y_AD%BladeLoad(1)%Force field [2 fields (force, moment) with 3 components]
      DO K = 1,p_FAST%nBeams ! Loop through all blades

         BD_Start = y_FAST%Lin%Modules(MODULE_BD)%Instance(k)%LinStartIndx(LIN_INPUT_COL) & ! start of BD%Input(1,k)%DistrLoad%Force field
                  + BD%Input(1,k)%RootMotion%NNodes *18  & ! displacement, rotation, & acceleration fields for each node
                  + BD%Input(1,k)%PointLoad%NNodes  * 6    ! force + moment fields for each node
            
            ! AD loads-to-BD loads transfer (dU^{BD}/dy^{AD}):
         call Assemble_dUdy_Loads(y_AD%rotors(1)%BladeLoad(k), BD%Input(1,k)%DistrLoad, MeshMapData%AD_L_2_BDED_B(k), BD_Start, AD_Out_Start, dUdy)
         AD_Out_Start = AD_Out_Start + y_AD%rotors(1)%BladeLoad(k)%NNodes*6  ! start of y_AD%BladeLoad(k+1)%Force field [skip the moments to get to forces on next blade]
         
         
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
      
      do k=1,size(u_AD%rotors(1)%InflowOnBlade,3) ! blades
         do j=1,size(u_AD%rotors(1)%InflowOnBlade,2) ! nodes
            do i=1,3 !velocity component
               dUdy( AD_Start + i - 1, y_FAST%Lin%Modules(MODULE_IfW)%Instance(1)%LinStartIndx(LIN_OUTPUT_COL) + (node-1)*3 + i - 1 ) = -1.0_R8Ki
            end do
            node = node + 1
            AD_Start = AD_Start + 3
         end do         
      end do
                  
      if ( allocated(u_AD%rotors(1)%InflowOnTower) ) then         
         do j=1,size(u_AD%rotors(1)%InflowOnTower,2) !nodes
            do i=1,3 !velocity component
               dUdy( AD_Start + i - 1, y_FAST%Lin%Modules(MODULE_IfW)%Instance(1)%LinStartIndx(LIN_OUTPUT_COL) + (node-1)*3 + i - 1 ) = -1.0_R8Ki
            end do
            node = node + 1
            AD_Start = AD_Start + 3
         end do      
      end if

      do i=1,3 !rotor-disk velocity component (DiskVel)
         dUdy( AD_Start + i - 1, y_FAST%Lin%Modules(MODULE_IfW)%Instance(1)%LinStartIndx(LIN_OUTPUT_COL) + (node-1)*3 + i - 1 ) = -1.0_R8Ki
      end do
      node = node + 1
      AD_Start = AD_Start + 3

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
   IF (u_AD%rotors(1)%TowerMotion%Committed) THEN
            
      !!! ! This linearization was done in forming dUdu (see Linear_AD_InputSolve_du()), so we don't need to re-calculate these matrices 
      !!! ! while forming dUdy, too.
      !!!CALL Linearize_Line2_to_Line2( y_ED%TowerLn2Mesh, u_AD%TowerMotion, MeshMapData%ED_L_2_AD_L_T, ErrStat2, ErrMsg2 )
      
      AD_Start = Indx_u_AD_Tower_Start(u_AD, y_FAST) ! start of u_AD%TowerMotion%TranslationDisp field
      
         ED_Out_Start = Indx_y_ED_Tower_Start(y_ED, y_FAST) ! start of y_ED%TowerLn2Mesh%TranslationDisp field
         call Assemble_dUdy_Motions(y_ED%TowerLn2Mesh, u_AD%rotors(1)%TowerMotion, MeshMapData%ED_L_2_AD_L_T, AD_Start, ED_Out_Start, dUdy, skipRotVel=.true.)
      
   END IF
      
      !...................................
      ! hub
      !...................................
      CALL Linearize_Point_to_Point( y_ED%HubPtMotion, u_AD%rotors(1)%HubMotion, MeshMapData%ED_P_2_AD_P_H, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//':u_AD%HubMotion' )
         if (errStat>=AbortErrLev) return
         
      ! *** AD translational displacement: from ED translational displacement (MeshMapData%ED_P_2_AD_P_H%dM%mi) and orientation (MeshMapData%ED_P_2_AD_P_H%dM%fx_p)
      AD_Start = Indx_u_AD_Hub_Start(u_AD, y_FAST) ! start of u_AD%HubMotion%TranslationDisp field   
      ED_Out_Start = Indx_y_ED_Hub_Start(y_ED, y_FAST) ! start of y_ED%HubPtMotion%TranslationDisp field
      call SetBlockMatrix( dUdy, MeshMapData%ED_P_2_AD_P_H%dM%mi, AD_Start, ED_Out_Start )
   
      ED_Out_Start = Indx_y_ED_Hub_Start(y_ED, y_FAST) + y_ED%HubPtMotion%NNodes * 3 ! start of y_ED%HubPtMotion%Orientation field
      call SetBlockMatrix( dUdy, MeshMapData%ED_P_2_AD_P_H%dM%fx_p, AD_Start, ED_Out_Start )
         
      ! *** AD orientation: from ED orientation
      AD_Start = AD_Start + u_AD%rotors(1)%HubMotion%NNodes * 3 ! move past the AD translation disp field to orientation field         
      call SetBlockMatrix( dUdy, MeshMapData%ED_P_2_AD_P_H%dM%mi, AD_Start, ED_Out_Start )
      
      ! *** AD rotational velocity: from ED rotational velocity
      AD_Start = AD_Start + u_AD%rotors(1)%HubMotion%NNodes * 3 ! move past the AD orientation field to rotational velocity field          
      ED_Out_Start = Indx_y_ED_Hub_Start(y_ED, y_FAST) + y_ED%HubPtMotion%NNodes * 6 ! ! start of y_ED%HubPtMotion%RotationVel field
      call SetBlockMatrix( dUdy, MeshMapData%ED_P_2_AD_P_H%dM%mi, AD_Start, ED_Out_Start )
         

   
      !...................................
      ! blade root   
      !...................................
      DO k=1,size(y_ED%BladeRootMotion)
         CALL Linearize_Point_to_Point( y_ED%BladeRootMotion(k), u_AD%rotors(1)%BladeRootMotion(k), MeshMapData%ED_P_2_AD_P_R(k), ErrStat2, ErrMsg2 )
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
         CALL Assemble_dUdy_Motions(y_ED%BladeLn2Mesh(k), u_AD%rotors(1)%BladeMotion(k), MeshMapData%BDED_L_2_AD_L_B(k), AD_Start, ED_Out_Start, dUdy, skipRotAcc=.true.)
         
      END DO
      
   ELSEIF (p_FAST%CompElast == Module_BD ) THEN
      !!! ! This linearization was done in forming dUdu (see Linear_AD_InputSolve_du()), so we don't need to re-calculate these matrices 
      !!! ! while forming dUdy, too.
      !!!CALL Linearize_Line2_to_Line2( BD%y(k)%BldMotion, u_AD%BladeMotion(k), MeshMapData%BDED_L_2_AD_L_B(k), ErrStat2, ErrMsg2 )
      
      DO k=1,p_FAST%nBeams
         AD_Start     = Indx_u_AD_Blade_Start(u_AD, y_FAST, k)     ! start of u_AD%BladeMotion(k)%TranslationDisp field
         BD_Out_Start = y_FAST%Lin%Modules(Module_BD)%Instance(k)%LinStartIndx(LIN_OUTPUT_COL)
         
         CALL Assemble_dUdy_Motions(BD%y(k)%BldMotion, u_AD%rotors(1)%BladeMotion(k), MeshMapData%BDED_L_2_AD_L_B(k), AD_Start, BD_Out_Start, dUdy, skipRotAcc=.true.)
      END DO
   
   END IF
   
   
END SUBROUTINE Linear_AD_InputSolve_NoIfW_dy
!----------------------------------------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------------------------------------
!> This routine forms the dU^{HD}/du^{HD} blocks of dUdu.
SUBROUTINE Linear_HD_InputSolve_du( p_FAST, y_FAST, u_HD, y_ED, y_SD, MeshMapData, dUdu, ErrStat, ErrMsg )

      ! Passed variables
   TYPE(FAST_ParameterType),    INTENT(IN   )   :: p_FAST      !< FAST parameter data    
   TYPE(FAST_OutputFileType),   INTENT(IN   )   :: y_FAST      !< FAST output file data (for linearization)
   TYPE(HydroDyn_InputType),    INTENT(INOUT)   :: u_HD        !< The inputs to HydroDyn
   TYPE(ED_OutputType),TARGET,  INTENT(IN)      :: y_ED        !< The outputs from the ElastoDyn structural dynamics module
   TYPE(SD_OutputType),TARGET,  INTENT(IN)      :: y_SD        !< The outputs from the SubDyn structural dynamics module
   TYPE(FAST_ModuleMapType),    INTENT(INOUT)   :: MeshMapData !< Data for mapping between modules
   REAL(R8Ki),                  INTENT(INOUT)   :: dUdu(:,:)   !< Jacobian matrix of which we are computing the dU^{HD}/du^{HD} block
   
   INTEGER(IntKi)                               :: ErrStat     !< Error status of the operation
   CHARACTER(*)                                 :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      ! Local variables:

   INTEGER(IntKi)                               :: HD_Start_td ! starting index of dUdu (column) where particular HD fields are located
   INTEGER(IntKi)                               :: HD_Start_tr ! starting index of dUdu (row) where particular HD fields are located
   INTEGER(IntKi)                               :: ErrStat2
   CHARACTER(ErrMsgLen)                         :: ErrMsg2 
   CHARACTER(*), PARAMETER                      :: RoutineName = 'Linear_HD_InputSolve_du'
   TYPE(MeshType), POINTER                      :: PlatformMotion
   TYPE(MeshType), POINTER                      :: SubstructureMotion2HD

   
   ErrStat = ErrID_None
   ErrMsg  = ""
               
   

      PlatformMotion => y_ED%PlatformPtMesh

   IF (p_FAST%CompSub == Module_SD) THEN
      SubstructureMotion2HD => y_SD%Y2Mesh
   ELSE
      SubstructureMotion2HD => PlatformMotion
   END IF   
   ! look at how the translational displacement gets transfered to the translational velocity and translational acceleration:
   !-------------------------------------------------------------------------------------------------
   ! Set the inputs from ElastoDyn:
   !-------------------------------------------------------------------------------------------------
      
   !..........
   ! dU^{HD}/du^{HD}
   ! note that the 1s on the diagonal have already been set, so we will fill in the off diagonal terms.
   !..........
   
   if ( p_FAST%CompHydro == Module_HD ) then ! HydroDyn-{ElastoDyn or SubDyn}
   
      !===================================================   
      !  y_ED%PlatformPtMesh and u_HD%PRPMesh ! this is always done with ED, even if using SD
      !=================================================== 
      
      ! Transfer ED motions to HD motion input (HD inputs depend on previously calculated HD inputs from ED):
      if ( u_HD%PRPMesh%Committed ) then
         call Linearize_Point_to_Point( PlatformMotion, u_HD%PRPMesh, MeshMapData%ED_P_2_HD_PRP_P, ErrStat2, ErrMsg2 )
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         
         ! HD is destination in the mapping, so we want M_{tv_uD} and M_{ta_uD}
            
         HD_Start_td = Indx_u_HD_PRP_Start(u_HD, y_FAST)  
         HD_Start_tr = HD_Start_td + u_HD%PRPMesh%NNodes * 6 ! skip 2 fields (TranslationDisp and Orientation) with 3 components before translational velocity field      
         
            ! translational velocity:
         if (allocated(MeshMapData%ED_P_2_HD_PRP_P%dM%tv_uD )) then             
            call SetBlockMatrix( dUdu, MeshMapData%ED_P_2_HD_PRP_P%dM%tv_ud, HD_Start_tr, HD_Start_td )
         end if
         
            ! translational acceleration:
         HD_Start_tr = HD_Start_tr + u_HD%PRPMesh%NNodes * 6 ! skip 2 fields ( TranslationVel and RotationVel)
         if (allocated(MeshMapData%ED_P_2_HD_PRP_P%dM%ta_uD )) then            
            call SetBlockMatrix( dUdu, MeshMapData%ED_P_2_HD_PRP_P%dM%ta_ud, HD_Start_tr, HD_Start_td )
         end if
      end if
      
         
      !===================================================
      !  y_ED%PlatformPtMesh or SD%y2Mesh and u_HD%Morison%Mesh
      !===================================================
         ! Transfer ED motions to HD motion input (HD inputs depend on previously calculated HD inputs from ED):
      if ( u_HD%Morison%Mesh%Committed ) then
         call Linearize_Point_to_Point( SubstructureMotion2HD, u_HD%Morison%Mesh, MeshMapData%SubStructure_2_HD_M_P, ErrStat2, ErrMsg2 )
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         
         ! HD is destination in the mapping, so we want M_{tv_uD} and M_{ta_uD}
            
         HD_Start_td = Indx_u_HD_Morison_Start(u_HD, y_FAST)
         HD_Start_tr = HD_Start_td + u_HD%Morison%Mesh%NNodes * 6 ! skip 2 fields (TranslationDisp and Orientation) with 3 components before translational velocity field
         
            ! translational velocity:
         if (allocated(MeshMapData%SubStructure_2_HD_M_P%dM%tv_uD )) then             
            call SetBlockMatrix( dUdu, MeshMapData%SubStructure_2_HD_M_P%dM%tv_ud, HD_Start_tr, HD_Start_td )
         end if
         
            ! translational acceleration:
         HD_Start_tr = HD_Start_tr + u_HD%Morison%Mesh%NNodes * 6 ! skip 2 fields ( TranslationVel and RotationVel)
         if (allocated(MeshMapData%SubStructure_2_HD_M_P%dM%ta_uD )) then            
            call SetBlockMatrix( dUdu, MeshMapData%SubStructure_2_HD_M_P%dM%ta_ud, HD_Start_tr, HD_Start_td )
         end if
      end if
      
      !===================================================   
      !  y_ED%PlatformPtMesh or SD%y2Mesh and u_HD%WAMITMesh
      !=================================================== 
      if ( u_HD%WAMITMesh%Committed ) then

         call Linearize_Point_to_Point( SubstructureMotion2HD, u_HD%WAMITMesh, MeshMapData%SubStructure_2_HD_W_P, ErrStat2, ErrMsg2 )
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
               
         HD_Start_td = Indx_u_HD_WAMIT_Start(u_HD, y_FAST)   
         HD_Start_tr = HD_Start_td + u_HD%WAMITMesh%NNodes  * 6 ! skip 2 fields (TranslationDisp and Orientation) with 3 components before translational velocity field
            
            ! translational velocity:
         if (allocated(MeshMapData%SubStructure_2_HD_W_P%dM%tv_uD )) then             
            call SetBlockMatrix( dUdu, MeshMapData%SubStructure_2_HD_W_P%dM%tv_ud, HD_Start_tr, HD_Start_td )
         end if

            ! translational acceleration:
         HD_Start_tr = HD_Start_tr + u_HD%WAMITMesh%NNodes * 6 ! skip 2 fields ( TranslationVel and RotationVel)
         if (allocated(MeshMapData%SubStructure_2_HD_W_P%dM%ta_uD )) then
            call SetBlockMatrix( dUdu, MeshMapData%SubStructure_2_HD_W_P%dM%ta_ud, HD_Start_tr, HD_Start_td )
         end if
      end if
      
   end if
   
   
END SUBROUTINE Linear_HD_InputSolve_du
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine forms the dU^{HD}/dy^{ED} block of dUdy. (i.e., how do changes in the ED outputs affect 
!! the HD inputs?)
SUBROUTINE Linear_HD_InputSolve_dy( p_FAST, y_FAST, u_HD, y_ED, y_SD, MeshMapData, dUdy, ErrStat, ErrMsg )

      ! Passed variables
   TYPE(FAST_ParameterType),    INTENT(IN   )   :: p_FAST      !< FAST parameter data    
   TYPE(FAST_OutputFileType),   INTENT(IN   )   :: y_FAST      !< FAST output file data (for linearization)
   TYPE(HydroDyn_InputType),    INTENT(INOUT)   :: u_HD        !< The inputs to HydroDyn
   TYPE(ED_OutputType), TARGET, INTENT(IN)      :: y_ED        !< The outputs from the ElastoDyn structural dynamics module
   TYPE(SD_OutputType), TARGET, INTENT(IN)      :: y_SD        !< The outputs from the SubDyn structural dynamics module
   TYPE(FAST_ModuleMapType),    INTENT(INOUT)   :: MeshMapData !< Data for mapping between modules
   REAL(R8Ki),                  INTENT(INOUT)   :: dUdy(:,:)   !< Jacobian matrix of which we are computing the dU^{HD}/dy^{ED} block
   
   INTEGER(IntKi)                               :: ErrStat     !< Error status of the operation
   CHARACTER(*)                                 :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      ! Local variables:

   INTEGER(IntKi)                               :: HD_Start    ! starting index of dUdy (column) where particular HD fields are located
   INTEGER(IntKi)                               :: Platform_Out_Start! starting index of dUdy (row) where particular ED fields are located
   INTEGER(IntKi)                               :: SubStructure_Out_Start! starting index of dUdy (row) where particular SD/ED fields are located
   TYPE(MeshType), POINTER                      :: PlatformMotion
   TYPE(MeshType), POINTER                      :: SubstructureMotion2HD

   CHARACTER(*), PARAMETER                      :: RoutineName = 'Linear_HD_InputSolve_dy'

   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
           
      PlatformMotion => y_ED%PlatformPtMesh
      Platform_Out_Start = Indx_y_ED_Platform_Start(y_ED, y_FAST) ! start of y_ED%PlatformPtMesh%TranslationDisp field
      
   IF (p_FAST%CompSub == Module_SD) THEN
      SubstructureMotion2HD => y_SD%y2Mesh
      SubStructure_Out_Start = Indx_y_SD_Y2Mesh_Start(y_SD, y_FAST)   ! start of y_SD%Y2Mesh%TranslationDisp field
   ELSE
      SubstructureMotion2HD => PlatformMotion
      SubStructure_Out_Start = Platform_Out_Start
   END IF
   
   
   !...................................
   ! HD PRP Mesh
   !...................................
   !  use Indx_u_HD_PRP_Start
   HD_Start     = Indx_u_HD_PRP_Start(u_HD, y_FAST)  ! start of u_HD%Morison%Mesh%TranslationDisp field     
   call Assemble_dUdy_Motions(PlatformMotion, u_HD%PRPMesh, MeshMapData%ED_P_2_HD_PRP_P, HD_Start, Platform_Out_Start, dUdy, .false.)
   
   
   ! dU^{HD}/dy^{ED} or ! dU^{HD}/dy^{SD}
      !...................................
      ! Morison Mesh
      !...................................
   IF (u_HD%Morison%Mesh%Committed) THEN
            
      !!! ! This linearization was done in forming dUdu (see Linear_HD_InputSolve_du()), so we don't need to re-calculate these matrices 
      !!! ! while forming dUdy, too.
      !!!call Linearize_Point_to_Line2( SubstructureMotion2HD, u_HD%Morison%Mesh, MeshMapData%SubStructure_2_HD_M_P, ErrStat2, ErrMsg2 )
      
      HD_Start     = Indx_u_HD_Morison_Start(u_HD, y_FAST)  ! start of u_HD%Morison%Mesh%TranslationDisp field
      call Assemble_dUdy_Motions(SubstructureMotion2HD, u_HD%Morison%Mesh, MeshMapData%SubStructure_2_HD_M_P, HD_Start, SubStructure_Out_Start, dUdy, .false.)
   END IF

      !...................................
      ! Lumped Platform Reference Pt Mesh
      !...................................
   IF (u_HD%WAMITMesh%Committed) THEN
            
      !!! ! This linearization was done in forming dUdu (see Linear_HD_InputSolve_du()), so we don't need to re-calculate these matrices 
      !!! ! while forming dUdy, too.
      !!!call Linearize_Point_to_Point( SubstructureMotion2HD, u_HD%Mesh, MeshMapData%SubStructure_2_HD_W_P, ErrStat2, ErrMsg2 )
      
      HD_Start     = Indx_u_HD_WAMIT_Start(u_HD, y_FAST)  ! start of u_HD%Mesh%TranslationDisp field
      call Assemble_dUdy_Motions(SubstructureMotion2HD, u_HD%WAMITMesh, MeshMapData%SubStructure_2_HD_W_P, HD_Start, SubStructure_Out_Start, dUdy, .false.)
   END IF

   
END SUBROUTINE Linear_HD_InputSolve_dy

!----------------------------------------------------------------------------------------------------------------------------------
!> This routine forms the dU^{MAP}/dy^{ED} block of dUdy. (i.e., how do changes in the ED outputs affect 
!! the MAP inputs?)
SUBROUTINE Linear_MAP_InputSolve_dy( p_FAST, y_FAST, u_MAP, y_ED, y_SD, MeshMapData, dUdy, ErrStat, ErrMsg )

      ! Passed variables
   TYPE(FAST_ParameterType),    INTENT(IN   )   :: p_FAST      !< FAST parameter data    
   TYPE(FAST_OutputFileType),   INTENT(IN   )   :: y_FAST      !< FAST output file data (for linearization)
   TYPE(MAP_InputType),         INTENT(INOUT)   :: u_MAP       !< The inputs to MAP
   TYPE(ED_OutputType), TARGET, INTENT(IN)      :: y_ED        !< The outputs from the ElastoDyn structural dynamics module
   TYPE(SD_OutputType), TARGET, INTENT(IN)      :: y_SD        !< The outputs from the SubDyn structural dynamics module
   TYPE(FAST_ModuleMapType),    INTENT(INOUT)   :: MeshMapData !< Data for mapping between modules
   REAL(R8Ki),                  INTENT(INOUT)   :: dUdy(:,:)   !< Jacobian matrix of which we are computing the dU^{MAP}/dy^{ED} block
   
   INTEGER(IntKi)                               :: ErrStat     !< Error status of the operation
   CHARACTER(*)                                 :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      ! Local variables:

   INTEGER(IntKi)                               :: MAP_Start   ! starting index of dUdy (column) where particular MAP fields are located

   INTEGER(IntKi)                               :: SubStructure_Out_Start! starting index of dUdy (row) where particular SD/ED fields are located
   TYPE(MeshType), POINTER                      :: SubstructureMotion
   
   INTEGER(IntKi)                               :: ErrStat2
   CHARACTER(ErrMsgLen)                         :: ErrMsg2 
   CHARACTER(*), PARAMETER                      :: RoutineName = 'Linear_MAP_InputSolve_dy'

   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   IF (p_FAST%CompSub == Module_SD) THEN
      SubstructureMotion => y_SD%y3Mesh
      SubStructure_Out_Start = Indx_y_SD_Y3Mesh_Start(y_SD, y_FAST) ! start of y_SD%Y3Mesh%TranslationDisp field
   ELSE
      SubstructureMotion => y_ED%PlatformPtMesh
      SubStructure_Out_Start = Indx_y_ED_Platform_Start(y_ED, y_FAST) ! start of y_ED%PlatformPtMesh%TranslationDisp field
   END IF
   
   IF (u_MAP%PtFairDisplacement%Committed) THEN 
         !...................................
         ! FairLead Mesh
         !...................................

      MAP_Start    = y_FAST%Lin%Modules(MODULE_MAP)%Instance(1)%LinStartIndx(LIN_INPUT_COL)
      
      ! dU^{MAP}/dy^{SD} or ! dU^{MAP}/dy^{ED}
      call Linearize_Point_to_Point( SubstructureMotion, u_MAP%PtFairDisplacement, MeshMapData%Structure_2_Mooring, ErrStat2, ErrMsg2 )
      call Assemble_dUdy_Motions(y_ED%PlatformPtMesh, u_MAP%PtFairDisplacement, MeshMapData%Structure_2_Mooring, MAP_Start, SubStructure_Out_Start, dUdy, OnlyTranslationDisp=.true.)

   END IF
END SUBROUTINE Linear_MAP_InputSolve_dy

!----------------------------------------------------------------------------------------------------------------------------------
!> This routine forms the dU^{MD}/du^{MD} block of dUdu. (i.e., how do changes in the MD outputs affect 
!! the MD inputs?)
SUBROUTINE Linear_MD_InputSolve_du( p_FAST, y_FAST, u_MD, y_ED, y_SD, MeshMapData, dUdu, ErrStat, ErrMsg )

      ! Passed variables
   TYPE(FAST_ParameterType),    INTENT(IN   )   :: p_FAST      !< FAST parameter data    
   TYPE(FAST_OutputFileType),   INTENT(IN   )   :: y_FAST      !< FAST output file data (for linearization)
   TYPE(MD_InputType),          INTENT(INOUT)   :: u_MD        !< The inputs to MoorDyn
   TYPE(ED_OutputType), TARGET, INTENT(IN)      :: y_ED        !< The outputs from the ElastoDyn structural dynamics module
   TYPE(SD_OutputType), TARGET, INTENT(IN)      :: y_SD        !< The outputs from the SubDyn structural dynamics module
   TYPE(FAST_ModuleMapType),    INTENT(INOUT)   :: MeshMapData !< Data for mapping between modules
   REAL(R8Ki),                  INTENT(INOUT)   :: dUdu(:,:)   !< Jacobian matrix of which we are computing the dU^{MD}/dy^{ED} block
   
   INTEGER(IntKi)                               :: ErrStat     !< Error status of the operation
   CHARACTER(*)                                 :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      ! Local variables:

   INTEGER(IntKi)                               :: MD_Start_td ! starting index of dUdu (column) where particular MD fields are located
   INTEGER(IntKi)                               :: MD_Start_tr ! starting index of dUdu (row) where particular MD fields are located
   
   TYPE(MeshType), POINTER                      :: SubstructureMotion
   
   INTEGER(IntKi)                               :: ErrStat2
   CHARACTER(ErrMsgLen)                         :: ErrMsg2 
   CHARACTER(*), PARAMETER                      :: RoutineName = 'Linear_MD_InputSolve_du'

   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   IF (p_FAST%CompSub == Module_SD) THEN
      SubstructureMotion => y_SD%y3Mesh
   ELSE
      SubstructureMotion => y_ED%PlatformPtMesh
   END IF
   
   IF (u_MD%CoupledKinematics(1)%Committed) THEN 
      !...................................
      ! FairLead Mesh
      !...................................

         ! dU^{MD}/du^{MD}
         call Linearize_Point_to_Point( SubstructureMotion, u_MD%CoupledKinematics(1), MeshMapData%Structure_2_Mooring, ErrStat2, ErrMsg2 )

      ! MD is destination in the mapping, so we want M_{tv_uD} and M_{ta_uD}
      MD_Start_td = y_FAST%Lin%Modules(MODULE_MD)%Instance(1)%LinStartIndx(LIN_INPUT_COL)
      MD_Start_tr = MD_Start_td + u_MD%CoupledKinematics(1)%NNodes * 6 ! skip 2 fields (TranslationDisp and Orientation) with 3 components before translational velocity field      
         
            ! translational velocity:
         if (allocated(MeshMapData%Structure_2_Mooring%dM%tv_uD )) then             
            call SetBlockMatrix( dUdu, MeshMapData%Structure_2_Mooring%dM%tv_ud, MD_Start_tr, MD_Start_td )
         end if
         
            ! translational acceleration:
         MD_Start_tr = MD_Start_tr + u_MD%CoupledKinematics(1)%NNodes * 6 ! skip 2 fields ( TranslationVel and RotationVel)
         if (allocated(MeshMapData%Structure_2_Mooring%dM%ta_uD )) then            
            call SetBlockMatrix( dUdu, MeshMapData%Structure_2_Mooring%dM%ta_ud, MD_Start_tr, MD_Start_td )
         end if

   END IF
END SUBROUTINE Linear_MD_InputSolve_du

!----------------------------------------------------------------------------------------------------------------------------------
!> This routine forms the dU^{MD}/dy^{ED} block of dUdy. (i.e., how do changes in the ED outputs affect 
!! the MD inputs?)
SUBROUTINE Linear_MD_InputSolve_dy( p_FAST, y_FAST, u_MD, y_ED, y_SD, MeshMapData, dUdy, ErrStat, ErrMsg )

      ! Passed variables
   TYPE(FAST_ParameterType),    INTENT(IN   )   :: p_FAST      !< FAST parameter data    
   TYPE(FAST_OutputFileType),   INTENT(IN   )   :: y_FAST      !< FAST output file data (for linearization)
   TYPE(MD_InputType),          INTENT(INOUT)   :: u_MD        !< The inputs to MoorDyn
   TYPE(ED_OutputType), TARGET, INTENT(IN)      :: y_ED        !< The outputs from the ElastoDyn structural dynamics module
   TYPE(SD_OutputType), TARGET, INTENT(IN)      :: y_SD        !< The outputs from the SubDyn structural dynamics module
   TYPE(FAST_ModuleMapType),    INTENT(INOUT)   :: MeshMapData !< Data for mapping between modules
   REAL(R8Ki),                  INTENT(INOUT)   :: dUdy(:,:)   !< Jacobian matrix of which we are computing the dU^{MD}/dy^{ED} block
   
   INTEGER(IntKi)                               :: ErrStat     !< Error status of the operation
   CHARACTER(*)                                 :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      ! Local variables:

   INTEGER(IntKi)                               :: MD_Start    ! starting index of dUdy (column) where particular MD fields are located
   INTEGER(IntKi)                               :: SubStructure_Out_Start! starting index of dUdy (row) where particular SD/ED fields are located
   TYPE(MeshType), POINTER                      :: SubstructureMotion
   
   CHARACTER(*), PARAMETER                      :: RoutineName = 'Linear_MD_InputSolve_dy'

   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   IF (u_MD%CoupledKinematics(1)%Committed) THEN

      IF (p_FAST%CompSub == Module_SD) THEN
         SubstructureMotion => y_SD%y3Mesh
         SubStructure_Out_Start = Indx_y_SD_Y3Mesh_Start(y_SD, y_FAST) ! start of y_SD%Y3Mesh%TranslationDisp field
      ELSE
         SubstructureMotion => y_ED%PlatformPtMesh
         SubStructure_Out_Start = Indx_y_ED_Platform_Start(y_ED, y_FAST) ! start of y_ED%PlatformPtMesh%TranslationDisp field
      END IF
   
         !...................................
         ! FairLead Mesh
         !...................................

      MD_Start    = y_FAST%Lin%Modules(MODULE_MD)%Instance(1)%LinStartIndx(LIN_INPUT_COL)
      
      ! dU^{MD}/dy^{SD} or   dU^{MD}/dy^{ED}

      !!! ! This linearization was done in forming dUdu (see Linear_MD_InputSolve_du()), so we don't need to re-calculate these matrices 
      !!! ! while forming dUdy, too.
      !!!call Linearize_Point_to_Point( SubstructureMotion, u_MD%CoupledKinematics(1), MeshMapData%Structure_2_Mooring, ErrStat2, ErrMsg2 )

      call Assemble_dUdy_Motions(    SubstructureMotion, u_MD%CoupledKinematics(1), MeshMapData%Structure_2_Mooring, MD_Start, SubStructure_Out_Start, dUdy, OnlyTranslationDisp=.false.)
      
   END IF
END SUBROUTINE Linear_MD_InputSolve_dy

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
   CHARACTER(ErrMsgLen)                    :: ErrMsg2             ! local error message
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
   CHARACTER(ErrMsgLen)                    :: ErrMsg2             ! local error message
   CHARACTER(*),             PARAMETER     :: RoutineName = 'Glue_StateMatrices' 
   
   ErrStat = ErrID_None
   ErrMsg = ""
   
   !LIN-TODO: Update doc string comments below for SrvD, HD, MAP, SD
   
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
   !! \begin{bmatrix} A^{ED} & 0 & 0 \\ 0 & A^{BD} & 0 \\ 0 & 0 & A^{HD}\end{bmatrix} - 
   !! \begin{bmatrix} 0 & 0 & B^{ED} & 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & B^{BD} & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 & 0 & B^{HD}\end{bmatrix} \,
   !! \begin{bmatrix} G \end{bmatrix}^{-1} \, \frac{\partial U}{\partial y} \, \begin{bmatrix} 0 & 0 & 0 \\ 0 & 0 & 0 \\ C^{ED} & 0 & 0 \\ 0 & C^{BD} & 0 \\ 0 & 0 & 0 \\ 0 & 0 & C^{HD} \\ 0 & 0 & 0\end{bmatrix}
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
SUBROUTINE Assemble_dUdy_Motions(y, u, MeshMap, BlockRowStart, BlockColStart, dUdy, skipRotVel, skipRotAcc, onlyTranslationDisp)
   TYPE(MeshType),    INTENT(IN)     :: y                      !< the output (source) mesh that is transfering motions
   TYPE(MeshType),    INTENT(IN)     :: u                      !< the input (destination) mesh that is receiving motions
   TYPE(MeshMapType), INTENT(IN)     :: MeshMap                !< the mesh mapping from y to u
   INTEGER(IntKi),    INTENT(IN)     :: BlockRowStart          !< the index of the row defining the block of dUdy to be set
   INTEGER(IntKi),    INTENT(IN)     :: BlockColStart          !< the index of the column defining the block of dUdy to be set
   REAL(R8Ki),        INTENT(INOUT)  :: dUdy(:,:)              !< full Jacobian matrix
   LOGICAL, OPTIONAL, INTENT(IN)     :: skipRotVel             !< if present and true, we skip the rotational velocity and both acceleration fields and return early
   LOGICAL, OPTIONAL, INTENT(IN)     :: onlyTranslationDisp    !< if present and true, we set only the destination translationDisp fields and return early
   LOGICAL, OPTIONAL, INTENT(IN)     :: skipRotAcc             !< if present and true, we skip the rotational acceleration field
   
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


      if (PRESENT(onlyTranslationDisp)) then
         if (onlyTranslationDisp) return ! destination includes only the translational displacement field, so we'll just return
      end if


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


      if (PRESENT(skipRotAcc)) then
         if (skipRotAcc) return ! destination does not include rotational accelerations, so we'll just return
      end if


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
      
      if (allocated(y%Moment)) then
         ! source moment to destination moment:
         row = BlockRowStart + u%NNodes*3       ! start of u%Moment field [skip 1 field with 3 components]
         col = BlockColStart + y%NNodes*3       ! start of y%Moment field [skip 1 field with 3 components]
         call SetBlockMatrix( dUdy, MeshMap%dM%li, row, col )
      end if
      
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
!> This routine returns the starting index for the y_ED%NacelleMotion mesh in the FAST linearization outputs.
FUNCTION Indx_y_ED_Nacelle_Start(y_ED, y_FAST) RESULT(ED_Out_Start)
   TYPE(FAST_OutputFileType),      INTENT(IN )  :: y_FAST           !< FAST output file data (for linearization)
   TYPE(ED_OutputType),            INTENT(IN )  :: y_ED             !< ED outputs at t
   INTEGER                                      :: k                !< blade number loop

   INTEGER                                      :: ED_Out_Start     !< starting index of this blade mesh in ElastoDyn outputs

   ED_Out_Start = Indx_y_ED_BladeRoot_Start(y_ED, y_FAST, size(y_ED%BladeRootMotion))        ! start of last blade root
   ED_Out_Start = ED_Out_Start + y_ED%BladeRootMotion(size(y_ED%BladeRootMotion))%NNodes*18  ! N blade roots, 6 fields with 3 components per blade.
END FUNCTION Indx_y_ED_Nacelle_Start
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine returns the starting index for y_ED%Yaw in the FAST linearization outputs.
FUNCTION Indx_y_Yaw_Start(y_FAST, ThisModule) RESULT(ED_Out_Start)
   TYPE(FAST_OutputFileType),      INTENT(IN )  :: y_FAST           !< FAST output file data (for linearization)
   INTEGER,                        INTENT(IN )  :: ThisModule       !< which structural module this is for

   INTEGER                                      :: ED_Out_Start     !< starting index of this blade mesh in ElastoDyn outputs

   
   ED_Out_Start = y_FAST%Lin%Modules(thisModule)%Instance(1)%LinStartIndx(LIN_OUTPUT_COL) + y_FAST%Lin%Modules(thisModule)%Instance(1)%SizeLin(LIN_OUTPUT_COL) & !end of ED outputs (+1)
                - y_FAST%Lin%Modules(thisModule)%Instance(1)%NumOutputs - 3 ! start of ED where Yaw, YawRate, HSS_Spd occur (right before WriteOutputs)
                
END FUNCTION Indx_y_Yaw_Start
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

   AD_Start = Indx_u_AD_Tower_Start(u_AD, y_FAST) + u_AD%rotors(1)%TowerMotion%NNodes * 9  ! 3 fields (MASKID_TRANSLATIONDISP,MASKID_Orientation,MASKID_TRANSLATIONVel) with 3 components

END FUNCTION Indx_u_AD_Hub_Start
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine returns the starting index for the u_AD%BladeRootMotion(k) mesh in the FAST linearization inputs.
FUNCTION Indx_u_AD_BladeRoot_Start(u_AD, y_FAST, BladeNum) RESULT(AD_Start)
   TYPE(FAST_OutputFileType),      INTENT(IN )  :: y_FAST           !< FAST output file data (for linearization)
   TYPE(AD_InputType),             INTENT(IN )  :: u_AD             !< AD Inputs at t
   INTEGER,                        INTENT(IN )  :: BladeNum         !< blade number to find index for
   INTEGER                                      :: k                !< blade number loop

   INTEGER                                      :: AD_Start         !< starting index of this mesh in AeroDyn inputs

   AD_Start = Indx_u_AD_Hub_Start(u_AD, y_FAST) + u_AD%rotors(1)%HubMotion%NNodes * 9  ! 3 fields (MASKID_TRANSLATIONDISP,MASKID_Orientation,MASKID_RotationVel) with 3 components
   
   do k = 1,min(BladeNum-1,size(u_AD%rotors(1)%BladeRootMotion))
      AD_Start = AD_Start + u_AD%rotors(1)%BladeRootMotion(k)%NNodes * 3 ! 1 field (MASKID_Orientation) with 3 components
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
   
   do k = 1,min(BladeNum-1,size(u_AD%rotors(1)%BladeMotion))
      AD_Start = AD_Start + u_AD%rotors(1)%BladeMotion(k)%NNodes * 15 ! 5 fields (TranslationDisp, MASKID_Orientation, TranslationVel, RotationVel, TranslationAcc) with 3 components
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

!----------------------------------------------------------------------------------------------------------------------------------
!> This routine returns the starting index for the u_SD%TPMesh mesh in the FAST linearization inputs.
FUNCTION Indx_u_SD_TPMesh_Start(u_SD, y_FAST) RESULT(SD_Start)
   TYPE(FAST_OutputFileType),      INTENT(IN )  :: y_FAST           !< FAST output file data (for linearization)
   TYPE(SD_InputType),         INTENT(IN )  :: u_SD             !< SD Inputs at t

   INTEGER                                      :: SD_Start         !< starting index of this mesh in SubDyn inputs

   SD_Start = y_FAST%Lin%Modules(Module_SD)%Instance(1)%LinStartIndx(LIN_INPUT_COL) 

END FUNCTION Indx_u_SD_TPMesh_Start
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine returns the starting index for the y_SD%Y1Mesh mesh in the FAST linearization outputs.
FUNCTION Indx_y_SD_Y1Mesh_Start(y_SD, y_FAST) RESULT(SD_Out_Start)
   TYPE(FAST_OutputFileType),      INTENT(IN )  :: y_FAST           !< FAST output file data (for linearization)
   TYPE(SD_OutputType),            INTENT(IN )  :: y_SD             !< SD outputs at t

   INTEGER                                      :: SD_Out_Start     !< starting index of this mesh in ElastoDyn outputs

   SD_Out_Start = y_FAST%Lin%Modules(MODULE_SD)%Instance(1)%LinStartIndx(LIN_OUTPUT_COL)
END FUNCTION Indx_y_SD_Y1Mesh_Start
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine returns the starting index for the y_SD%Y2Mesh mesh in the FAST linearization outputs.
FUNCTION Indx_y_SD_Y2Mesh_Start(y_SD, y_FAST) RESULT(SD_Out_Start)
   TYPE(FAST_OutputFileType),      INTENT(IN )  :: y_FAST           !< FAST output file data (for linearization)
   TYPE(SD_OutputType),            INTENT(IN )  :: y_SD             !< SD outputs at t

   INTEGER                                      :: SD_Out_Start     !< starting index of this mesh in ElastoDyn outputs

   SD_Out_Start = Indx_y_SD_Y1Mesh_Start(y_SD, y_FAST) + y_SD%Y1Mesh%NNodes * 6            ! 3 forces + 3 moments at each node! skip all of the Y1Mesh data and get to the beginning of 
END FUNCTION Indx_y_SD_Y2Mesh_Start
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine returns the starting index for the y_SD%Y3Mesh mesh in the FAST linearization outputs.
FUNCTION Indx_y_SD_Y3Mesh_Start(y_SD, y_FAST) RESULT(SD_Out_Start)
   TYPE(FAST_OutputFileType),      INTENT(IN )  :: y_FAST           !< FAST output file data (for linearization)
   TYPE(SD_OutputType),            INTENT(IN )  :: y_SD             !< SD outputs at t

   INTEGER                                      :: SD_Out_Start     !< starting index of this mesh in ElastoDyn outputs

   SD_Out_Start = Indx_y_SD_Y2Mesh_Start(y_SD, y_FAST) + y_SD%Y2Mesh%NNodes * 6            ! 3 forces + 3 moments at each node! skip all of the Y1Mesh data and get to the beginning of 
END FUNCTION Indx_y_SD_Y3Mesh_Start
!----------------------------------------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------------------------------------
!> This routine returns the starting index for the u_SD%TPMesh mesh in the FAST linearization inputs.
FUNCTION Indx_u_SD_LMesh_Start(u_SD, y_FAST) RESULT(SD_Start)
   TYPE(FAST_OutputFileType),      INTENT(IN )  :: y_FAST           !< FAST output file data (for linearization)
   TYPE(SD_InputType),         INTENT(IN )  :: u_SD             !< SD Inputs at t

   INTEGER                                      :: SD_Start         !< starting index of this mesh in SubDyn inputs

   SD_Start = Indx_u_SD_TPMesh_Start(u_SD, y_FAST) + u_SD%TPMesh%NNodes*18 ! 6 fields with 3 components 

END FUNCTION Indx_u_SD_LMesh_Start
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine returns the starting index for the u_HD%Morison%Mesh mesh in the FAST linearization inputs.
FUNCTION Indx_u_HD_Morison_Start(u_HD, y_FAST) RESULT(HD_Start)
   TYPE(FAST_OutputFileType),      INTENT(IN )  :: y_FAST           !< FAST output file data (for linearization)
   TYPE(HydroDyn_InputType),       INTENT(IN )  :: u_HD             !< HD Inputs at t

   INTEGER                                      :: HD_Start         !< starting index of this mesh in HydroDyn inputs

   HD_Start = y_FAST%Lin%Modules(Module_HD)%Instance(1)%LinStartIndx(LIN_INPUT_COL) 

END FUNCTION Indx_u_HD_Morison_Start
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine returns the starting index for the u_HD%WAMITMesh mesh in the FAST linearization inputs.
FUNCTION Indx_u_HD_WAMIT_Start(u_HD, y_FAST) RESULT(HD_Start)
   TYPE(FAST_OutputFileType),      INTENT(IN )  :: y_FAST           !< FAST output file data (for linearization)
   TYPE(HydroDyn_InputType),       INTENT(IN )  :: u_HD             !< HD Inputs at t

   INTEGER                                      :: HD_Start         !< starting index of this mesh in HydroDyn inputs

   HD_Start = Indx_u_HD_Morison_Start(u_HD, y_FAST) 
   if (u_HD%Morison%Mesh%committed)  HD_Start =  HD_Start + u_HD%Morison%Mesh%NNodes * 18  ! 6 fields (MASKID_TRANSLATIONDISP,MASKID_Orientation,MASKID_TRANSLATIONVel,MASKID_ROTATIONVel,MASKID_TRANSLATIONAcc,MASKID_ROTATIONAcc) with 3 components

END FUNCTION Indx_u_HD_WAMIT_Start
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine returns the starting index for the u_HD%PRPMesh mesh in the FAST linearization inputs.
FUNCTION Indx_u_HD_PRP_Start(u_HD, y_FAST) RESULT(HD_Start)
   TYPE(FAST_OutputFileType),      INTENT(IN )  :: y_FAST           !< FAST output file data (for linearization)
   TYPE(HydroDyn_InputType),       INTENT(IN )  :: u_HD             !< HD Inputs at t

   INTEGER                                      :: HD_Start         !< starting index of this mesh in HydroDyn inputs

   HD_Start = Indx_u_HD_WAMIT_Start(u_HD, y_FAST) 
   if (u_HD%WAMITMesh%committed)  HD_Start =  HD_Start + u_HD%WAMITMesh%NNodes * 18  ! 6 fields (MASKID_TRANSLATIONDISP,MASKID_Orientation,MASKID_TRANSLATIONVel,MASKID_ROTATIONVel,MASKID_TRANSLATIONAcc,MASKID_ROTATIONAcc) with 3 components

   END FUNCTION Indx_u_HD_PRP_Start
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine returns the starting index for the y_HD%Morison%DistribMesh mesh in the FAST linearization outputs.
FUNCTION Indx_y_HD_Morison_Start(y_HD, y_FAST) RESULT(HD_Start)
   TYPE(FAST_OutputFileType),      INTENT(IN )  :: y_FAST           !< FAST output file data (for linearization)
   TYPE(HydroDyn_OutputType),      INTENT(IN )  :: y_HD             !< HD Outputs at t

   INTEGER                                      :: HD_Start         !< starting index of this mesh in HydroDyn Outputs

   HD_Start = y_FAST%Lin%Modules(Module_HD)%Instance(1)%LinStartIndx(LIN_OUTPUT_COL) 

END FUNCTION Indx_y_HD_Morison_Start
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine returns the starting index for the y_HD%Mesh mesh in the FAST linearization outputs.
FUNCTION Indx_y_HD_WAMIT_Start(y_HD, y_FAST) RESULT(HD_Start)
   TYPE(FAST_OutputFileType),      INTENT(IN )  :: y_FAST           !< FAST output file data (for linearization)
   TYPE(HydroDyn_OutputType),      INTENT(IN )  :: y_HD             !< HD Outputs at t

   INTEGER                                      :: HD_Start
   
      !< starting index of this mesh in HydroDyn Outputs

   HD_Start = Indx_y_HD_Morison_Start(y_HD, y_FAST) 
   if (y_HD%Morison%Mesh%committed)  HD_Start =  HD_Start + y_HD%Morison%Mesh%NNodes * 6  ! 2 fields (MASKID_FORCE,MASKID_MOMENT) with 3 components

   END FUNCTION Indx_y_HD_WAMIT_Start
!----------------------------------------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine allocates the arrays that store the operating point at each linearization time for later producing VTK
!! files of the mode shapes.
SUBROUTINE AllocateOP(p_FAST, y_FAST, ErrStat, ErrMsg )

   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST              !< Output variables for the glue code
     
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None

   ! local variables
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(*), PARAMETER                 :: RoutineName = 'AllocateOP'


   ErrStat = ErrID_None
   ErrMsg  = ""
   
   !----------------------------------------------------------------------------------------
   !! copy the operating point of the states and inputs at LinTimes(i)
   !----------------------------------------------------------------------------------------
   

      ALLOCATE(       y_FAST%op%x_ED(p_FAST%NLinTimes), STAT=ErrStat2 )
         if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
      ALLOCATE(      y_FAST%op%xd_ED(p_FAST%NLinTimes), STAT=ErrStat2 )
         if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
      ALLOCATE(       y_FAST%op%z_ED(p_FAST%NLinTimes), STAT=ErrStat2 )
         if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
      ALLOCATE( y_FAST%op%OtherSt_ED(p_FAST%NLinTimes), STAT=ErrStat2 )
         if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
      ALLOCATE(       y_FAST%op%u_ED(p_FAST%NLinTimes), STAT=ErrStat2 )
         if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
      
      IF ( p_FAST%CompElast == Module_BD ) THEN
         ALLOCATE(       y_FAST%op%x_BD(p_FAST%nBeams, p_FAST%NLinTimes), STAT=ErrStat2 )
            if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
         ALLOCATE(      y_FAST%op%xd_BD(p_FAST%nBeams, p_FAST%NLinTimes), STAT=ErrStat2 )
            if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
         ALLOCATE(       y_FAST%op%z_BD(p_FAST%nBeams, p_FAST%NLinTimes), STAT=ErrStat2 )
            if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
         ALLOCATE( y_FAST%op%OtherSt_BD(p_FAST%nBeams, p_FAST%NLinTimes), STAT=ErrStat2 )
            if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
         ALLOCATE(       y_FAST%op%u_BD(p_FAST%nBeams, p_FAST%NLinTimes), STAT=ErrStat2 )
            if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
      END IF

   
      
   !IF ( p_FAST%CompAero == Module_AD14 ) THEN
   !ELSE
   IF ( p_FAST%CompAero == Module_AD ) THEN
      ALLOCATE(       y_FAST%op%x_AD(p_FAST%NLinTimes), STAT=ErrStat2 )
         if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
      ALLOCATE(      y_FAST%op%xd_AD(p_FAST%NLinTimes), STAT=ErrStat2 )
         if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
      ALLOCATE(       y_FAST%op%z_AD(p_FAST%NLinTimes), STAT=ErrStat2 )
         if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
      ALLOCATE( y_FAST%op%OtherSt_AD(p_FAST%NLinTimes), STAT=ErrStat2 )
         if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
      ALLOCATE(       y_FAST%op%u_AD(p_FAST%NLinTimes), STAT=ErrStat2 )
         if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
   END IF
         
   IF ( p_FAST%CompInflow == Module_IfW ) THEN
      ALLOCATE(       y_FAST%op%x_IfW(p_FAST%NLinTimes), STAT=ErrStat2 )
         if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
      ALLOCATE(      y_FAST%op%xd_IfW(p_FAST%NLinTimes), STAT=ErrStat2 )
         if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
      ALLOCATE(       y_FAST%op%z_IfW(p_FAST%NLinTimes), STAT=ErrStat2 )
         if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
      ALLOCATE( y_FAST%op%OtherSt_IfW(p_FAST%NLinTimes), STAT=ErrStat2 )
         if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
      ALLOCATE(       y_FAST%op%u_IfW(p_FAST%NLinTimes), STAT=ErrStat2 )
         if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
   END IF
            
      
   IF ( p_FAST%CompServo == Module_SrvD ) THEN
      ALLOCATE(       y_FAST%op%x_SrvD(p_FAST%NLinTimes), STAT=ErrStat2 )
         if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
      ALLOCATE(      y_FAST%op%xd_SrvD(p_FAST%NLinTimes), STAT=ErrStat2 )
         if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
      ALLOCATE(       y_FAST%op%z_SrvD(p_FAST%NLinTimes), STAT=ErrStat2 )
         if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
      ALLOCATE( y_FAST%op%OtherSt_SrvD(p_FAST%NLinTimes), STAT=ErrStat2 )
         if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
      ALLOCATE(       y_FAST%op%u_SrvD(p_FAST%NLinTimes), STAT=ErrStat2 )
         if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
   END IF
      
      
   IF ( p_FAST%CompHydro == Module_HD ) THEN
      ALLOCATE(       y_FAST%op%x_HD(p_FAST%NLinTimes), STAT=ErrStat2 )
         if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
      ALLOCATE(      y_FAST%op%xd_HD(p_FAST%NLinTimes), STAT=ErrStat2 )
         if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
      ALLOCATE(       y_FAST%op%z_HD(p_FAST%NLinTimes), STAT=ErrStat2 )
         if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
      ALLOCATE( y_FAST%op%OtherSt_HD(p_FAST%NLinTimes), STAT=ErrStat2 )
         if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
      ALLOCATE(       y_FAST%op%u_HD(p_FAST%NLinTimes), STAT=ErrStat2 )
         if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
   END IF
            
            
   ! SubDyn: copy final predictions to actual states
   IF ( p_FAST%CompSub == Module_SD ) THEN
      ALLOCATE(       y_FAST%op%x_SD(p_FAST%NLinTimes), STAT=ErrStat2 )
         if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
      ALLOCATE(      y_FAST%op%xd_SD(p_FAST%NLinTimes), STAT=ErrStat2 )
         if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
      ALLOCATE(       y_FAST%op%z_SD(p_FAST%NLinTimes), STAT=ErrStat2 )
         if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
      ALLOCATE( y_FAST%op%OtherSt_SD(p_FAST%NLinTimes), STAT=ErrStat2 )
         if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
      ALLOCATE(       y_FAST%op%u_SD(p_FAST%NLinTimes), STAT=ErrStat2 )
         if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
   ELSE IF ( p_FAST%CompSub == Module_ExtPtfm ) THEN
      ALLOCATE(       y_FAST%op%x_ExtPtfm(p_FAST%NLinTimes), STAT=ErrStat2 )
         if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
      ALLOCATE(      y_FAST%op%xd_ExtPtfm(p_FAST%NLinTimes), STAT=ErrStat2 )
         if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
      ALLOCATE(       y_FAST%op%z_ExtPtfm(p_FAST%NLinTimes), STAT=ErrStat2 )
         if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
      ALLOCATE( y_FAST%op%OtherSt_ExtPtfm(p_FAST%NLinTimes), STAT=ErrStat2 )
         if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
      ALLOCATE(       y_FAST%op%u_ExtPtfm(p_FAST%NLinTimes), STAT=ErrStat2 )
         if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
   END IF
         
      
   ! MAP/MoorDyn/FEAM: copy states and inputs to OP array
   IF (p_FAST%CompMooring == Module_MAP) THEN
      ALLOCATE(       y_FAST%op%x_MAP(p_FAST%NLinTimes), STAT=ErrStat2 )
         if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
      ALLOCATE(      y_FAST%op%xd_MAP(p_FAST%NLinTimes), STAT=ErrStat2 )
         if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
      ALLOCATE(       y_FAST%op%z_MAP(p_FAST%NLinTimes), STAT=ErrStat2 )
         if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
      !ALLOCATE( y_FAST%op%OtherSt_MAP(p_FAST%NLinTimes), STAT=ErrStat2 )
      !   if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
      ALLOCATE(       y_FAST%op%u_MAP(p_FAST%NLinTimes), STAT=ErrStat2 )
         if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
   ELSEIF (p_FAST%CompMooring == Module_MD) THEN
      ALLOCATE(       y_FAST%op%x_MD(p_FAST%NLinTimes), STAT=ErrStat2 )
         if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
      ALLOCATE(      y_FAST%op%xd_MD(p_FAST%NLinTimes), STAT=ErrStat2 )
         if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
      ALLOCATE(       y_FAST%op%z_MD(p_FAST%NLinTimes), STAT=ErrStat2 )
         if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
      ALLOCATE( y_FAST%op%OtherSt_MD(p_FAST%NLinTimes), STAT=ErrStat2 )
         if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
      ALLOCATE(       y_FAST%op%u_MD(p_FAST%NLinTimes), STAT=ErrStat2 )
         if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
   ELSEIF (p_FAST%CompMooring == Module_FEAM) THEN
      ALLOCATE(       y_FAST%op%x_FEAM(p_FAST%NLinTimes), STAT=ErrStat2 )
         if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
      ALLOCATE(      y_FAST%op%xd_FEAM(p_FAST%NLinTimes), STAT=ErrStat2 )
         if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
      ALLOCATE(       y_FAST%op%z_FEAM(p_FAST%NLinTimes), STAT=ErrStat2 )
         if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
      ALLOCATE( y_FAST%op%OtherSt_FEAM(p_FAST%NLinTimes), STAT=ErrStat2 )
         if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
      ALLOCATE(       y_FAST%op%u_FEAM(p_FAST%NLinTimes), STAT=ErrStat2 )
         if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
   !ELSEIF (p_FAST%CompMooring == Module_Orca) THEN
   END IF
             
         ! IceFloe/IceDyn: copy states and inputs to OP array
   IF ( p_FAST%CompIce == Module_IceF ) THEN
      ALLOCATE(       y_FAST%op%x_IceF(p_FAST%NLinTimes), STAT=ErrStat2 )
         if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
      ALLOCATE(      y_FAST%op%xd_IceF(p_FAST%NLinTimes), STAT=ErrStat2 )
         if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
      ALLOCATE(       y_FAST%op%z_IceF(p_FAST%NLinTimes), STAT=ErrStat2 )
         if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
      ALLOCATE( y_FAST%op%OtherSt_IceF(p_FAST%NLinTimes), STAT=ErrStat2 )
         if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
      ALLOCATE(       y_FAST%op%u_IceF(p_FAST%NLinTimes), STAT=ErrStat2 )
         if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
   ELSEIF ( p_FAST%CompIce == Module_IceD ) THEN
      ALLOCATE(       y_FAST%op%x_IceD(p_FAST%numIceLegs, p_FAST%NLinTimes), STAT=ErrStat2 )
         if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
      ALLOCATE(      y_FAST%op%xd_IceD(p_FAST%numIceLegs, p_FAST%NLinTimes), STAT=ErrStat2 )
         if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
      ALLOCATE(       y_FAST%op%z_IceD(p_FAST%numIceLegs, p_FAST%NLinTimes), STAT=ErrStat2 )
         if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
      ALLOCATE( y_FAST%op%OtherSt_IceD(p_FAST%numIceLegs, p_FAST%NLinTimes), STAT=ErrStat2 )
         if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
      ALLOCATE(       y_FAST%op%u_IceD(p_FAST%numIceLegs, p_FAST%NLinTimes), STAT=ErrStat2 )
         if (ErrStat2 /= 0) call SetErrStat( ErrID_Fatal, 'Error allocating arrays for VTK operating points.', ErrStat, ErrMsg, RoutineName)
   END IF
   
END SUBROUTINE AllocateOP
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine is the inverse of SetOperatingPoint(). It saves the current operating points so they can be retrieved 
!> when visualizing mode shapes.
SUBROUTINE SaveOP(i, p_FAST, y_FAST, ED, BD, SrvD, AD, IfW, OpFM, HD, SD, ExtPtfm, MAPp, FEAM, MD, Orca, &
                         IceF, IceD, ErrStat, ErrMsg, CtrlCode )

   INTEGER(IntKi)          , INTENT(IN   ) :: i                   !< current index into LinTimes
   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST              !< Output variables for the glue code
     
   TYPE(ElastoDyn_Data),     INTENT(INOUT) :: ED                  !< ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(INOUT) :: BD                  !< BeamDyn data
   TYPE(ServoDyn_Data),      INTENT(INOUT) :: SrvD                !< ServoDyn data
   TYPE(AeroDyn_Data),       INTENT(INOUT) :: AD                  !< AeroDyn data
   TYPE(InflowWind_Data),    INTENT(INOUT) :: IfW                 !< InflowWind data
   TYPE(OpenFOAM_Data),      INTENT(INOUT) :: OpFM                !< OpenFOAM data
   TYPE(HydroDyn_Data),      INTENT(INOUT) :: HD                  !< HydroDyn data
   TYPE(SubDyn_Data),        INTENT(INOUT) :: SD                  !< SubDyn data
   TYPE(ExtPtfm_Data),       INTENT(INOUT) :: ExtPtfm             !< ExtPtfm_MCKF data
   TYPE(MAP_Data),           INTENT(INOUT) :: MAPp                !< MAP data
   TYPE(FEAMooring_Data),    INTENT(INOUT) :: FEAM                !< FEAMooring data
   TYPE(MoorDyn_Data),       INTENT(INOUT) :: MD                  !< Data for the MoorDyn module
   TYPE(OrcaFlex_Data),      INTENT(INOUT) :: Orca                !< OrcaFlex interface data
   TYPE(IceFloe_Data),       INTENT(INOUT) :: IceF                !< IceFloe data
   TYPE(IceDyn_Data),        INTENT(INOUT) :: IceD                !< All the IceDyn data used in time-step loop

   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None
   INTEGER(IntKi),           INTENT(IN   ) :: CtrlCode            !< mesh copy control code (new, vs update)

   ! local variables
   INTEGER(IntKi)                          :: k                   ! generic loop counters
   
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMsg2
   CHARACTER(*), PARAMETER                 :: RoutineName = 'SaveOP'


   ErrStat = ErrID_None
   ErrMsg  = ""
   
   
   !----------------------------------------------------------------------------------------
   !! copy the operating point of the states and inputs at LinTimes(i)
   !----------------------------------------------------------------------------------------
      
         ! ElastoDyn: copy states and inputs to OP array
      CALL ED_CopyContState   (ED%x( STATE_CURR), y_FAST%op%x_ED( i), CtrlCode, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL ED_CopyDiscState   (ED%xd(STATE_CURR), y_FAST%op%xd_ED( i), CtrlCode, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL ED_CopyConstrState (ED%z( STATE_CURR), y_FAST%op%z_ED( i), CtrlCode, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL ED_CopyOtherState (ED%OtherSt( STATE_CURR), y_FAST%op%OtherSt_ED( i), CtrlCode, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
      CALL ED_CopyInput (ED%Input(1), y_FAST%op%u_ED( i), CtrlCode, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
         ! BeamDyn: copy states and inputs to OP array
      IF ( p_FAST%CompElast == Module_BD ) THEN
         DO k=1,p_FAST%nBeams
            CALL BD_CopyContState   (BD%x( k,STATE_CURR), y_FAST%op%x_BD(k, i), CtrlCode, Errstat2, ErrMsg2)
               CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            CALL BD_CopyDiscState   (BD%xd(k,STATE_CURR), y_FAST%op%xd_BD(k, i), CtrlCode, Errstat2, ErrMsg2)
               CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            CALL BD_CopyConstrState (BD%z( k,STATE_CURR), y_FAST%op%z_BD(k, i), CtrlCode, Errstat2, ErrMsg2)
               CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            CALL BD_CopyOtherState (BD%OtherSt( k,STATE_CURR), y_FAST%op%OtherSt_BD(k, i), CtrlCode, Errstat2, ErrMsg2)
               CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                     
            CALL BD_CopyInput (BD%Input(1,k), y_FAST%op%u_BD(k, i), CtrlCode, Errstat2, ErrMsg2)
               CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                     
         END DO
      END IF

   
      
      ! AeroDyn: copy states and inputs to OP array
   !IF ( p_FAST%CompAero == Module_AD14 ) THEN
   !ELSE
   IF ( p_FAST%CompAero == Module_AD ) THEN
      CALL AD_CopyContState   (AD%x( STATE_CURR), y_FAST%op%x_AD( i), CtrlCode, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AD_CopyDiscState   (AD%xd(STATE_CURR), y_FAST%op%xd_AD( i), CtrlCode, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AD_CopyConstrState (AD%z( STATE_CURR), y_FAST%op%z_AD( i), CtrlCode, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AD_CopyOtherState (AD%OtherSt(STATE_CURR), y_FAST%op%OtherSt_AD( i), CtrlCode, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      CALL AD_CopyInput (AD%Input(1), y_FAST%op%u_AD(i), CtrlCode, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   END IF
         
   ! InflowWind: copy states and inputs to OP array
   IF ( p_FAST%CompInflow == Module_IfW ) THEN
      CALL InflowWind_CopyContState   (IfW%x( STATE_CURR), y_FAST%op%x_IfW( i), CtrlCode, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL InflowWind_CopyDiscState   (IfW%xd(STATE_CURR), y_FAST%op%xd_IfW( i), CtrlCode, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL InflowWind_CopyConstrState (IfW%z( STATE_CURR), y_FAST%op%z_IfW( i), CtrlCode, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL InflowWind_CopyOtherState( IfW%OtherSt( STATE_CURR), y_FAST%op%OtherSt_IfW( i), CtrlCode, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               
      CALL InflowWind_CopyInput (IfW%Input(1), y_FAST%op%u_IfW(i), CtrlCode, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               
   END IF
            
      
   ! ServoDyn: copy states and inputs to OP array
   IF ( p_FAST%CompServo == Module_SrvD ) THEN
      CALL SrvD_CopyContState   (SrvD%x( STATE_CURR), y_FAST%op%x_SrvD( i), CtrlCode, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL SrvD_CopyDiscState   (SrvD%xd(STATE_CURR), y_FAST%op%xd_SrvD( i), CtrlCode, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL SrvD_CopyConstrState (SrvD%z( STATE_CURR), y_FAST%op%z_SrvD( i), CtrlCode, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL SrvD_CopyOtherState (SrvD%OtherSt( STATE_CURR), y_FAST%op%OtherSt_SrvD( i), CtrlCode, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               
      CALL SrvD_CopyInput (SrvD%Input(1), y_FAST%op%u_SrvD(i), CtrlCode, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   END IF
      
      
   ! HydroDyn: copy states and inputs to OP array
   IF ( p_FAST%CompHydro == Module_HD ) THEN
      CALL HydroDyn_CopyContState   (HD%x( STATE_CURR), y_FAST%op%x_HD( i), CtrlCode, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL HydroDyn_CopyDiscState   (HD%xd(STATE_CURR), y_FAST%op%xd_HD( i), CtrlCode, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL HydroDyn_CopyConstrState (HD%z( STATE_CURR), y_FAST%op%z_HD( i), CtrlCode, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL HydroDyn_CopyOtherState (HD%OtherSt(STATE_CURR), y_FAST%op%OtherSt_HD( i), CtrlCode, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               
      CALL HydroDyn_CopyInput (HD%Input(1), y_FAST%op%u_HD(i), CtrlCode, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   END IF
            
            
   ! SubDyn: copy final predictions to actual states
   IF ( p_FAST%CompSub == Module_SD ) THEN
      CALL SD_CopyContState   (y_FAST%op%x_SD(i), SD%x( STATE_CURR), CtrlCode, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL SD_CopyDiscState   (y_FAST%op%xd_SD(i), SD%xd(STATE_CURR), CtrlCode, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL SD_CopyConstrState( y_FAST%op%z_SD(i), SD%z( STATE_CURR), CtrlCode, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL SD_CopyOtherState (y_FAST%op%OtherSt_SD(i), SD%OtherSt(STATE_CURR), CtrlCode, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      CALL SD_CopyInput (y_FAST%op%u_SD(i), SD%Input(1), CtrlCode, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   ELSE IF ( p_FAST%CompSub == Module_ExtPtfm ) THEN
      CALL ExtPtfm_CopyContState   (ExtPtfm%x( STATE_CURR), y_FAST%op%x_ExtPtfm(i), CtrlCode, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL ExtPtfm_CopyDiscState   (ExtPtfm%xd(STATE_CURR), y_FAST%op%xd_ExtPtfm(i), CtrlCode, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL ExtPtfm_CopyConstrState (ExtPtfm%z( STATE_CURR), y_FAST%op%z_ExtPtfm(i), CtrlCode, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL ExtPtfm_CopyOtherState (ExtPtfm%OtherSt(STATE_CURR), y_FAST%op%OtherSt_ExtPtfm(i), CtrlCode, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      CALL ExtPtfm_CopyInput (ExtPtfm%Input(1), y_FAST%op%u_ExtPtfm(i), CtrlCode, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   END IF
         
      
   ! MAP/MoorDyn/FEAM: copy states and inputs to OP array
   IF (p_FAST%CompMooring == Module_MAP) THEN
      CALL MAP_CopyContState   (MAPp%x( STATE_CURR), y_FAST%op%x_MAP(i), CtrlCode, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL MAP_CopyDiscState   (MAPp%xd(STATE_CURR), y_FAST%op%xd_MAP(i), CtrlCode, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL MAP_CopyConstrState (MAPp%z( STATE_CURR), y_FAST%op%z_MAP(i), CtrlCode, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      !CALL MAP_CopyOtherState (MAPp%OtherSt(STATE_CURR), y_FAST%op%OtherSt_MAP(i), CtrlCode, Errstat2, ErrMsg2)
      !   CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      CALL MAP_CopyInput (MAPp%Input(1), y_FAST%op%u_MAP(i), CtrlCode, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         
   ELSEIF (p_FAST%CompMooring == Module_MD) THEN
      CALL MD_CopyContState   (MD%x( STATE_CURR), y_FAST%op%x_MD(i), CtrlCode, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL MD_CopyDiscState   (MD%xd(STATE_CURR), y_FAST%op%xd_MD(i), CtrlCode, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL MD_CopyConstrState (MD%z( STATE_CURR), y_FAST%op%z_MD(i), CtrlCode, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL MD_CopyOtherState (MD%OtherSt(STATE_CURR), y_FAST%op%OtherSt_MD(i), CtrlCode, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               
      CALL MD_CopyInput (MD%Input(1), y_FAST%op%u_MD(i), CtrlCode, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         
   ELSEIF (p_FAST%CompMooring == Module_FEAM) THEN
      CALL FEAM_CopyContState   (FEAM%x( STATE_CURR), y_FAST%op%x_FEAM(i), CtrlCode, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL FEAM_CopyDiscState   (FEAM%xd(STATE_CURR), y_FAST%op%xd_FEAM(i), CtrlCode, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL FEAM_CopyConstrState (FEAM%z( STATE_CURR), y_FAST%op%z_FEAM(i), CtrlCode, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL FEAM_CopyOtherState (FEAM%OtherSt( STATE_CURR), y_FAST%op%OtherSt_FEAM(i), CtrlCode, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               
      CALL FEAM_CopyInput (FEAM%Input(1), y_FAST%op%u_FEAM(i), CtrlCode, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   !ELSEIF (p_FAST%CompMooring == Module_Orca) THEN
   END IF
             
         ! IceFloe/IceDyn: copy states and inputs to OP array
   IF ( p_FAST%CompIce == Module_IceF ) THEN
      CALL IceFloe_CopyContState   (IceF%x( STATE_CURR), y_FAST%op%x_IceF(i), CtrlCode, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL IceFloe_CopyDiscState   (IceF%xd(STATE_CURR), y_FAST%op%xd_IceF(i), CtrlCode, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL IceFloe_CopyConstrState (IceF%z( STATE_CURR), y_FAST%op%z_IceF(i), CtrlCode, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL IceFloe_CopyOtherState (IceF%OtherSt(STATE_CURR), y_FAST%op%OtherSt_IceF(i), CtrlCode, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               
      CALL IceFloe_CopyInput (IceF%Input(1), y_FAST%op%u_IceF(i), CtrlCode, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   ELSEIF ( p_FAST%CompIce == Module_IceD ) THEN
      DO k=1,p_FAST%numIceLegs
         CALL IceD_CopyContState   (IceD%x( k,STATE_CURR), y_FAST%op%x_IceD(k, i), CtrlCode, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         CALL IceD_CopyDiscState   (IceD%xd(k,STATE_CURR), y_FAST%op%xd_IceD(k, i), CtrlCode, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         CALL IceD_CopyConstrState (IceD%z( k,STATE_CURR), y_FAST%op%z_IceD(k, i), CtrlCode, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         CALL IceD_CopyOtherState (IceD%OtherSt( k,STATE_CURR), y_FAST%op%OtherSt_IceD(k, i), CtrlCode, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                  
         CALL IceD_CopyInput (IceD%Input(1,k), y_FAST%op%u_IceD(k, i), CtrlCode, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END DO
   END IF   
   
   
END SUBROUTINE SaveOP
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine takes arrays representing the eigenvector of the states and uses it to modify the operating points for 
!! continuous states. It is highly tied to the module organizaton.
SUBROUTINE PerturbOP(t, iLinTime, iMode, p_FAST, y_FAST, ED, BD, SrvD, AD, IfW, OpFM, HD, SD, ExtPtfm, MAPp, FEAM, MD, Orca, &
                         IceF, IceD, ErrStat, ErrMsg )

   REAL(DbKi),               INTENT(IN   ) :: t
   INTEGER(IntKi),           INTENT(IN   ) :: iLinTime            !< index into LinTimes dimension of arrays (azimuth)
   INTEGER(IntKi),           INTENT(IN   ) :: iMode               !< index into Mode dimension of arrays

   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST              !< Output variables for the glue code
     
   TYPE(ElastoDyn_Data),     INTENT(INOUT) :: ED                  !< ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(INOUT) :: BD                  !< BeamDyn data
   TYPE(ServoDyn_Data),      INTENT(INOUT) :: SrvD                !< ServoDyn data
   TYPE(AeroDyn_Data),       INTENT(INOUT) :: AD                  !< AeroDyn data
   TYPE(InflowWind_Data),    INTENT(INOUT) :: IfW                 !< InflowWind data
   TYPE(OpenFOAM_Data),      INTENT(INOUT) :: OpFM                !< OpenFOAM data
   TYPE(HydroDyn_Data),      INTENT(INOUT) :: HD                  !< HydroDyn data
   TYPE(SubDyn_Data),        INTENT(INOUT) :: SD                  !< SubDyn data
   TYPE(ExtPtfm_Data),       INTENT(INOUT) :: ExtPtfm             !< ExtPtfm_MCKF data
   TYPE(MAP_Data),           INTENT(INOUT) :: MAPp                !< MAP data
   TYPE(FEAMooring_Data),    INTENT(INOUT) :: FEAM                !< FEAMooring data
   TYPE(MoorDyn_Data),       INTENT(INOUT) :: MD                  !< Data for the MoorDyn module
   TYPE(OrcaFlex_Data),      INTENT(INOUT) :: Orca                !< OrcaFlex interface data
   TYPE(IceFloe_Data),       INTENT(INOUT) :: IceF                !< IceFloe data
   TYPE(IceDyn_Data),        INTENT(INOUT) :: IceD                !< All the IceDyn data used in time-step loop

     
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None

   ! local variables
   INTEGER(IntKi)                          :: k                   ! generic loop counters
   INTEGER(IntKi)                          :: i, iStart           ! generic loop counters
   INTEGER(IntKi)                          :: iBody               ! WAMIT body loop counter
   INTEGER(IntKi)                          :: j                   ! generic loop counters
   INTEGER(IntKi)                          :: indx                ! generic loop counters
   INTEGER(IntKi)                          :: indx_last           ! generic loop counters
   INTEGER(IntKi)                          :: i_x                 ! index into packed array
   INTEGER(IntKi)                          :: nStates             ! number of second-order states
   INTEGER(IntKi)                          :: ThisModule          ! identifier of current module
   
   CHARACTER(*), PARAMETER                 :: RoutineName = 'PerturbOP'


   ErrStat = ErrID_None
   ErrMsg  = ""


   i_x = 1

   do i = 1,p_FAST%Lin_NumMods
      ThisModule = p_FAST%Lin_ModOrder( i )

      do k=1,size(y_FAST%Lin%Modules(ThisModule)%Instance)

         if (allocated(y_FAST%Lin%Modules(ThisModule)%Instance(k)%op_x_eig_mag)) then
            do j=1,size(y_FAST%Lin%Modules(ThisModule)%Instance(k)%op_x)  ! use this for the loop because ED may have a larger op_x_eig_mag array than op_x
         
            ! this is a hack because not all modules pack the continuous states in the same way:
               if (ThisModule == Module_ED) then
                  if (j<= ED%p%DOFs%NActvDOF) then
                     indx = ED%p%DOFs%PS(j)
                  else
                     indx = ED%p%DOFs%PS(j-ED%p%DOFs%NActvDOF) + ED%p%NDOF
                  end if
               else
                  indx = j
               end if
               y_FAST%Lin%Modules(ThisModule)%Instance(k)%op_x_eig_mag(  indx) = p_FAST%VTK_modes%x_eig_magnitude(i_x, iLinTime, iMode)  ! this is going to hold the magnitude of the eigenvector
               y_FAST%Lin%Modules(ThisModule)%Instance(k)%op_x_eig_phase(indx) = p_FAST%VTK_modes%x_eig_phase(    i_x, iLinTime, iMode) ! this is going to hold the phase of the eigenvector
               i_x = i_x + 1;
            end do
         end if
         
      end do
   end do
   
   
   
         ! ElastoDyn:
      ThisModule = Module_ED
      if (allocated(y_FAST%Lin%Modules(ThisModule)%Instance(1)%op_x_eig_mag)) then
         nStates = size(y_FAST%Lin%Modules(ThisModule)%Instance(1)%op_x_eig_mag)/2
         
         call GetStateAry(p_FAST, iMode, t, ED%x( STATE_CURR)%QT,  y_FAST%Lin%Modules(ThisModule)%Instance(1)%op_x_eig_mag(         :nStates), y_FAST%Lin%Modules(ThisModule)%Instance(1)%op_x_eig_phase(         :nStates))
         call GetStateAry(p_FAST, iMode, t, ED%x( STATE_CURR)%QDT, y_FAST%Lin%Modules(ThisModule)%Instance(1)%op_x_eig_mag(1+nStates:       ), y_FAST%Lin%Modules(ThisModule)%Instance(1)%op_x_eig_phase(1+nStates:       ))
      end if

         ! BeamDyn:
      IF ( p_FAST%CompElast == Module_BD ) THEN
         ThisModule = Module_BD
         DO k=1,p_FAST%nBeams
            if (allocated(y_FAST%Lin%Modules(ThisModule)%Instance(k)%op_x_eig_mag)) then
               nStates = size(y_FAST%Lin%Modules(ThisModule)%Instance(k)%op_x_eig_mag)/2

               indx = 1
               do i=2,BD%p(k)%node_total
                  indx_last = indx + BD%p(k)%dof_node - 1
                  call GetStateAry(p_FAST, iMode, t, BD%x(k, STATE_CURR)%q(   :,i), y_FAST%Lin%Modules(ThisModule)%Instance(k)%op_x_eig_mag(        indx:indx_last        ), y_FAST%Lin%Modules(ThisModule)%Instance(k)%op_x_eig_phase(        indx:indx_last        ))
                  call GetStateAry(p_FAST, iMode, t, BD%x(k, STATE_CURR)%dqdt(:,i), y_FAST%Lin%Modules(ThisModule)%Instance(k)%op_x_eig_mag(nStates+indx:indx_last+nStates), y_FAST%Lin%Modules(ThisModule)%Instance(k)%op_x_eig_phase(nStates+indx:indx_last+nStates))
                  indx = indx_last+1
               end do
               
            end if

         END DO
      END IF
   
      
   !!!   ! AeroDyn: copy final predictions to actual states; copy current outputs to next 
   !!!!IF ( p_FAST%CompAero == Module_AD14 ) THEN
   !!!!ELSE
   IF ( p_FAST%CompAero == Module_AD ) THEN
      ThisModule = Module_AD
      if (allocated(y_FAST%Lin%Modules(ThisModule)%Instance(1)%op_x_eig_mag)) then
      
         indx = 1
            ! set linearization operating points:
         if (AD%p%rotors(1)%BEMT%DBEMT%lin_nx>0) then
            do j=1,size(AD%x(STATE_CURR)%rotors(1)%BEMT%DBEMT%element,2)
               do i=1,size(AD%x(STATE_CURR)%rotors(1)%BEMT%DBEMT%element,1)
                  indx_last = indx + size(AD%x(STATE_CURR)%rotors(1)%BEMT%DBEMT%element(i,j)%vind) - 1
                  call GetStateAry(p_FAST, iMode, t, AD%x(STATE_CURR)%rotors(1)%BEMT%DBEMT%element(i,j)%vind, y_FAST%Lin%Modules(ThisModule)%Instance(1)%op_x_eig_mag(  indx : indx_last), &
                                                                                                              y_FAST%Lin%Modules(ThisModule)%Instance(1)%op_x_eig_phase(indx : indx_last) )
                  indx = indx_last + 1
               end do
            end do
   
            do j=1,size(AD%x(STATE_CURR)%rotors(1)%BEMT%DBEMT%element,2)
               do i=1,size(AD%x(STATE_CURR)%rotors(1)%BEMT%DBEMT%element,1)
                  indx_last = indx + size(AD%x(STATE_CURR)%rotors(1)%BEMT%DBEMT%element(i,j)%vind_1) - 1
                  call GetStateAry(p_FAST, iMode, t, AD%x(STATE_CURR)%rotors(1)%BEMT%DBEMT%element(i,j)%vind_1, y_FAST%Lin%Modules(ThisModule)%Instance(1)%op_x_eig_mag(  indx : indx_last), &
                                                                                                                  y_FAST%Lin%Modules(ThisModule)%Instance(1)%op_x_eig_phase(indx : indx_last) )
                  indx = indx_last + 1
               end do
            end do
      
         end if
   
         if (AD%p%rotors(1)%BEMT%UA%lin_nx>0) then
            do j=1,size(AD%x(STATE_CURR)%rotors(1)%BEMT%UA%element,2)
               do i=1,size(AD%x(STATE_CURR)%rotors(1)%BEMT%UA%element,1)
                  indx_last = indx + size(AD%x(STATE_CURR)%rotors(1)%BEMT%UA%element(i,j)%x) - 1
                  call GetStateAry(p_FAST, iMode, t, AD%x(STATE_CURR)%rotors(1)%BEMT%UA%element(i,j)%x,  y_FAST%Lin%Modules(ThisModule)%Instance(1)%op_x_eig_mag(  indx : indx_last), &
                                                                                                         y_FAST%Lin%Modules(ThisModule)%Instance(1)%op_x_eig_phase(indx : indx_last) )
                  indx = indx_last + 1
               end do
            end do
         end if
      
      end if
   END IF
   !!!      
   !!!! InflowWind: copy op to actual states and inputs
   !!!IF ( p_FAST%CompInflow == Module_IfW ) THEN
   !!!END IF
   !!!         
   !!!   
   !!!! ServoDyn: copy op to actual states and inputs
   !!!IF ( p_FAST%CompServo == Module_SrvD ) THEN
   !!!END IF
      
   ! HydroDyn: copy op to actual states and inputs
   IF ( p_FAST%CompHydro == Module_HD ) THEN
      ThisModule = Module_HD
      if (allocated(y_FAST%Lin%Modules(ThisModule)%Instance(1)%op_x_eig_mag)) then
         ! WAMIT parameter and continuous states are now an arrays of length NBody
         ! All Excitation states are stored first and then Radiation states.
         ! We will try to loop over each of the NBody(s) and add each body's states to the overall state array
         iStart = 1
         do iBody = 1, HD%p%NBody
            nStates = HD%p%WAMIT(iBody)%SS_Exctn%numStates
            if (nStates > 0) then
               call GetStateAry(p_FAST, iMode, t, HD%x( STATE_CURR)%WAMIT(iBody)%SS_Exctn%x, y_FAST%Lin%Modules(ThisModule)%Instance(1)%op_x_eig_mag(iStart:iStart+nStates-1), y_FAST%Lin%Modules(ThisModule)%Instance(1)%op_x_eig_phase(iStart:iStart+nStates-1))
               iStart = iStart + nStates
            end if
         end do
         do iBody = 1, HD%p%NBody
            nStates = HD%p%WAMIT(iBody)%SS_Rdtn%numStates
            if (nStates > 0) then
               call GetStateAry(p_FAST, iMode, t, HD%x( STATE_CURR)%WAMIT(iBody)%SS_Rdtn%x,  y_FAST%Lin%Modules(ThisModule)%Instance(1)%op_x_eig_mag(iStart:iStart+nStates-1), y_FAST%Lin%Modules(ThisModule)%Instance(1)%op_x_eig_phase(iStart:iStart+nStates-1))
               iStart = iStart + nStates
            end if
         end do
      end if
   END IF
                 
   
   ! SubDyn: copy final predictions to actual states
   IF ( p_FAST%CompSub == Module_SD ) THEN
      ThisModule = Module_SD
      if (allocated(y_FAST%Lin%Modules(ThisModule)%Instance(1)%op_x_eig_mag)) then
         nStates = size(y_FAST%Lin%Modules(ThisModule)%Instance(1)%op_x_eig_mag)/2
         call GetStateAry(p_FAST, iMode, t, SD%x( STATE_CURR)%qm,    y_FAST%Lin%Modules(ThisModule)%Instance(1)%op_x_eig_mag(         :nStates), y_FAST%Lin%Modules(ThisModule)%Instance(1)%op_x_eig_phase(         :nStates))
         call GetStateAry(p_FAST, iMode, t, SD%x( STATE_CURR)%qmdot, y_FAST%Lin%Modules(ThisModule)%Instance(1)%op_x_eig_mag(1+nStates:       ), y_FAST%Lin%Modules(ThisModule)%Instance(1)%op_x_eig_phase(1+nStates:       ))
      end if
   ELSE IF ( p_FAST%CompSub == Module_ExtPtfm ) THEN
   END IF   

END SUBROUTINE PerturbOP
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SetOperatingPoint(i, p_FAST, y_FAST, m_FAST, ED, BD, SrvD, AD, IfW, OpFM, HD, SD, ExtPtfm, &
                         MAPp, FEAM, MD, Orca, IceF, IceD, ErrStat, ErrMsg )

   INTEGER(IntKi),           INTENT(IN   ) :: i                   !< Index into LinTimes (to determine which operating point to copy)
   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST              !< Output variables for the glue code
   TYPE(FAST_MiscVarType),   INTENT(INOUT) :: m_FAST              !< Miscellaneous variables
     
   TYPE(ElastoDyn_Data),     INTENT(INOUT) :: ED                  !< ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(INOUT) :: BD                  !< BeamDyn data
   TYPE(ServoDyn_Data),      INTENT(INOUT) :: SrvD                !< ServoDyn data
   TYPE(AeroDyn_Data),       INTENT(INOUT) :: AD                  !< AeroDyn data
   TYPE(InflowWind_Data),    INTENT(INOUT) :: IfW                 !< InflowWind data
   TYPE(OpenFOAM_Data),      INTENT(INOUT) :: OpFM                !< OpenFOAM data
   TYPE(HydroDyn_Data),      INTENT(INOUT) :: HD                  !< HydroDyn data
   TYPE(SubDyn_Data),        INTENT(INOUT) :: SD                  !< SubDyn data
   TYPE(ExtPtfm_Data),       INTENT(INOUT) :: ExtPtfm             !< ExtPtfm_MCKF data
   TYPE(MAP_Data),           INTENT(INOUT) :: MAPp                !< MAP data
   TYPE(FEAMooring_Data),    INTENT(INOUT) :: FEAM                !< FEAMooring data
   TYPE(MoorDyn_Data),       INTENT(INOUT) :: MD                  !< Data for the MoorDyn module
   TYPE(OrcaFlex_Data),      INTENT(INOUT) :: Orca                !< OrcaFlex interface data
   TYPE(IceFloe_Data),       INTENT(INOUT) :: IceF                !< IceFloe data
   TYPE(IceDyn_Data),        INTENT(INOUT) :: IceD                !< All the IceDyn data used in time-step loop

   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None
   
   ! local variables
   INTEGER(IntKi)                          :: k                   ! generic loop counters
   
   
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMsg2
   CHARACTER(*), PARAMETER                 :: RoutineName = 'SetOperatingPoint'


   ErrStat = ErrID_None
   ErrMsg  = ""
   
   
   !----------------------------------------------------------------------------------------
   !! copy the operating point of the states and inputs at LinTimes(i)
   !----------------------------------------------------------------------------------------
         ! ElastoDyn: copy op to actual states and inputs
      CALL ED_CopyContState   (y_FAST%op%x_ED( i), ED%x( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL ED_CopyDiscState   (y_FAST%op%xd_ED( i), ED%xd(STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL ED_CopyConstrState (y_FAST%op%z_ED( i), ED%z( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL ED_CopyOtherState (y_FAST%op%OtherSt_ED( i), ED%OtherSt( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
      CALL ED_CopyInput (y_FAST%op%u_ED( i), ED%Input(1), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
         ! BeamDyn: copy op to actual states and inputs
      IF ( p_FAST%CompElast == Module_BD ) THEN
         DO k=1,p_FAST%nBeams
            CALL BD_CopyContState   (y_FAST%op%x_BD(k, i), BD%x( k,STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
               CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            CALL BD_CopyDiscState   (y_FAST%op%xd_BD(k, i), BD%xd(k,STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
               CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            CALL BD_CopyConstrState (y_FAST%op%z_BD(k, i), BD%z( k,STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
               CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            CALL BD_CopyOtherState (y_FAST%op%OtherSt_BD(k, i), BD%OtherSt( k,STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
               CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                     
            CALL BD_CopyInput (y_FAST%op%u_BD(k, i), BD%Input(1,k), MESH_UPDATECOPY, Errstat2, ErrMsg2)
               CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                     
         END DO
      END IF
      
      ! AeroDyn: copy final predictions to actual states; copy current outputs to next 
   !IF ( p_FAST%CompAero == Module_AD14 ) THEN
   !ELSE
   IF ( p_FAST%CompAero == Module_AD ) THEN
      CALL AD_CopyContState   (y_FAST%op%x_AD( i), AD%x( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AD_CopyDiscState   (y_FAST%op%xd_AD( i), AD%xd(STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AD_CopyConstrState (y_FAST%op%z_AD( i), AD%z( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AD_CopyOtherState (y_FAST%op%OtherSt_AD( i), AD%OtherSt(STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      CALL AD_CopyInput (y_FAST%op%u_AD(i), AD%Input(1), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   END IF
         
   ! InflowWind: copy op to actual states and inputs
   IF ( p_FAST%CompInflow == Module_IfW ) THEN
      CALL InflowWind_CopyContState   (y_FAST%op%x_IfW( i), IfW%x( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL InflowWind_CopyDiscState   (y_FAST%op%xd_IfW( i), IfW%xd(STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL InflowWind_CopyConstrState (y_FAST%op%z_IfW( i), IfW%z( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL InflowWind_CopyOtherState (y_FAST%op%OtherSt_IfW( i), IfW%OtherSt( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               
      CALL InflowWind_CopyInput (y_FAST%op%u_IfW(i), IfW%Input(1), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               
   END IF
            
      
   ! ServoDyn: copy op to actual states and inputs
   IF ( p_FAST%CompServo == Module_SrvD ) THEN
      CALL SrvD_CopyContState   (y_FAST%op%x_SrvD( i), SrvD%x( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL SrvD_CopyDiscState   (y_FAST%op%xd_SrvD( i), SrvD%xd(STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL SrvD_CopyConstrState (y_FAST%op%z_SrvD( i), SrvD%z( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL SrvD_CopyOtherState (y_FAST%op%OtherSt_SrvD( i), SrvD%OtherSt( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               
      CALL SrvD_CopyInput (y_FAST%op%u_SrvD(i), SrvD%Input(1), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   END IF
      
      
   ! HydroDyn: copy op to actual states and inputs
   IF ( p_FAST%CompHydro == Module_HD ) THEN         
      CALL HydroDyn_CopyContState   (y_FAST%op%x_HD( i), HD%x( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL HydroDyn_CopyDiscState   (y_FAST%op%xd_HD( i), HD%xd(STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL HydroDyn_CopyConstrState (y_FAST%op%z_HD( i), HD%z( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL HydroDyn_CopyOtherState (y_FAST%op%OtherSt_HD( i), HD%OtherSt(STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               
      CALL HydroDyn_CopyInput (y_FAST%op%u_HD(i), HD%Input(1), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   END IF
            
            
   ! SubDyn: copy final predictions to actual states
   IF ( p_FAST%CompSub == Module_SD ) THEN
      CALL SD_CopyContState   (y_FAST%op%x_SD(i), SD%x( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL SD_CopyDiscState   (y_FAST%op%xd_SD(i), SD%xd(STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL SD_CopyConstrState( y_FAST%op%z_SD(i), SD%z( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL SD_CopyOtherState (y_FAST%op%OtherSt_SD(i), SD%OtherSt(STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      CALL SD_CopyInput (y_FAST%op%u_SD(i), SD%Input(1), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   ELSE IF ( p_FAST%CompSub == Module_ExtPtfm ) THEN
      CALL ExtPtfm_CopyContState   (y_FAST%op%x_ExtPtfm(i), ExtPtfm%x( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL ExtPtfm_CopyDiscState   (y_FAST%op%xd_ExtPtfm(i), ExtPtfm%xd(STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL ExtPtfm_CopyConstrState (y_FAST%op%z_ExtPtfm(i), ExtPtfm%z( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL ExtPtfm_CopyOtherState (y_FAST%op%OtherSt_ExtPtfm(i), ExtPtfm%OtherSt(STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      CALL ExtPtfm_CopyInput (y_FAST%op%u_ExtPtfm(i), ExtPtfm%Input(1), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   END IF
         
      
   ! MAP/MoorDyn/FEAM: copy op to actual states and inputs
   IF (p_FAST%CompMooring == Module_MAP) THEN
      CALL MAP_CopyContState   (y_FAST%op%x_MAP(i), MAPp%x( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL MAP_CopyDiscState   (y_FAST%op%xd_MAP(i), MAPp%xd(STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL MAP_CopyConstrState (y_FAST%op%z_MAP(i), MAPp%z( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      !CALL MAP_CopyOtherState (y_FAST%op%OtherSt_MAP(i), MAPp%OtherSt(STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
      !   CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      CALL MAP_CopyInput (y_FAST%op%u_MAP(i), MAPp%Input(1), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   ELSEIF (p_FAST%CompMooring == Module_MD) THEN
      CALL MD_CopyContState   (y_FAST%op%x_MD(i), MD%x( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL MD_CopyDiscState   (y_FAST%op%xd_MD(i), MD%xd(STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL MD_CopyConstrState (y_FAST%op%z_MD(i), MD%z( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL MD_CopyOtherState (y_FAST%op%OtherSt_MD(i), MD%OtherSt(STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               
      CALL MD_CopyInput (y_FAST%op%u_MD(i), MD%Input(1), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   ELSEIF (p_FAST%CompMooring == Module_FEAM) THEN
      CALL FEAM_CopyContState   (y_FAST%op%x_FEAM(i), FEAM%x( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL FEAM_CopyDiscState   (y_FAST%op%xd_FEAM(i), FEAM%xd(STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL FEAM_CopyConstrState (y_FAST%op%z_FEAM(i), FEAM%z( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL FEAM_CopyOtherState (y_FAST%op%OtherSt_FEAM(i), FEAM%OtherSt( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               
      CALL FEAM_CopyInput (y_FAST%op%u_FEAM(i), FEAM%Input(1), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   !ELSEIF (p_FAST%CompMooring == Module_Orca) THEN
   END IF
             
         ! IceFloe/IceDyn: copy op to actual states and inputs
   IF ( p_FAST%CompIce == Module_IceF ) THEN
      CALL IceFloe_CopyContState   (y_FAST%op%x_IceF(i), IceF%x( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL IceFloe_CopyDiscState   (y_FAST%op%xd_IceF(i), IceF%xd(STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL IceFloe_CopyConstrState (y_FAST%op%z_IceF(i), IceF%z( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL IceFloe_CopyOtherState (y_FAST%op%OtherSt_IceF(i), IceF%OtherSt(STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               
      CALL IceFloe_CopyInput (y_FAST%op%u_IceF(i), IceF%Input(1), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   ELSEIF ( p_FAST%CompIce == Module_IceD ) THEN
      DO k=1,p_FAST%numIceLegs
         CALL IceD_CopyContState   (y_FAST%op%x_IceD(k, i), IceD%x( k,STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         CALL IceD_CopyDiscState   (y_FAST%op%xd_IceD(k, i), IceD%xd(k,STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         CALL IceD_CopyConstrState (y_FAST%op%z_IceD(k, i), IceD%z( k,STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         CALL IceD_CopyOtherState (y_FAST%op%OtherSt_IceD(k, i), IceD%OtherSt( k,STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                  
         CALL IceD_CopyInput (y_FAST%op%u_IceD(k, i), IceD%Input(1,k), MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END DO
   END IF

END SUBROUTINE SetOperatingPoint
!----------------------------------------------------------------------------------------------------------------------------------
subroutine GetStateAry(p_FAST, iMode, t, x, x_eig_magnitude, x_eig_phase)
   INTEGER(IntKi),           INTENT(IN   ) :: iMode               !< index into Mode dimension of arrays
   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
   REAL(DbKi)              , INTENT(IN   ) :: t                   !< time
   REAL(R8Ki),               INTENT(INOUT) :: x(:)                !< in: state at its operating point; out: added perturbation
   REAL(R8Ki),               INTENT(IN)    :: x_eig_magnitude(:)  !< magnitude of the eigenvector
   REAL(R8Ki),               INTENT(IN)    :: x_eig_phase(:)      !< phase of the eigenvector
      
   ! note that this assumes p_FAST%VTK_modes%VTKLinPhase is zero for VTKLinTim=2
   x = x + x_eig_magnitude * p_FAST%VTK_modes%VTKLinScale * cos( TwoPi_D * p_FAST%VTK_modes%DampedFreq_Hz(iMode)*t + x_eig_phase + p_FAST%VTK_modes%VTKLinPhase )
end subroutine GetStateAry



!----------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine performs the algorithm for computing a periodic steady-state solution.
SUBROUTINE FAST_CalcSteady( n_t_global, t_global, p_FAST, y_FAST, m_FAST, ED, BD, SrvD, AD, IfW, OpFM, HD, SD, ExtPtfm, MAPp, FEAM, MD, Orca, &
                   IceF, IceD, ErrStat, ErrMsg )

   INTEGER(IntKi),           INTENT(IN   ) :: n_t_global          !< integer time step
   REAL(DbKi),               INTENT(IN   ) :: t_global            ! current simulation time
   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST              !< Output variables for the glue code
   TYPE(FAST_MiscVarType),   INTENT(INOUT) :: m_FAST              !< Miscellaneous variables
     
   TYPE(ElastoDyn_Data),     INTENT(INOUT) :: ED                  !< ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(INOUT) :: BD                  !< BeamDyn data
   TYPE(ServoDyn_Data),      INTENT(INOUT) :: SrvD                !< ServoDyn data
   TYPE(AeroDyn_Data),       INTENT(INOUT) :: AD                  !< AeroDyn data
   TYPE(InflowWind_Data),    INTENT(INOUT) :: IfW                 !< InflowWind data
   TYPE(OpenFOAM_Data),      INTENT(INOUT) :: OpFM                !< OpenFOAM data
   TYPE(HydroDyn_Data),      INTENT(INOUT) :: HD                  !< HydroDyn data
   TYPE(SubDyn_Data),        INTENT(INOUT) :: SD                  !< SubDyn data
   TYPE(ExtPtfm_Data),       INTENT(INOUT) :: ExtPtfm             !< ExtPtfm data
   TYPE(MAP_Data),           INTENT(INOUT) :: MAPp                !< MAP data
   TYPE(FEAMooring_Data),    INTENT(INOUT) :: FEAM                !< FEAMooring data
   TYPE(MoorDyn_Data),       INTENT(INOUT) :: MD                  !< Data for the MoorDyn module
   TYPE(OrcaFlex_Data),      INTENT(INOUT) :: Orca                !< OrcaFlex interface data
   TYPE(IceFloe_Data),       INTENT(INOUT) :: IceF                !< IceFloe data
   TYPE(IceDyn_Data),        INTENT(INOUT) :: IceD                !< All the IceDyn data used in time-step loop

   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None

   ! local variables
   REAL(DbKi)                              :: DeltaAzim
   REAL(DbKi)                              :: psi                 !< psi (rotor azimuth) at which the outputs are defined
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMsg2
   LOGICAL                                 :: NextAzimuth
   
   CHARACTER(*), PARAMETER                 :: RoutineName = 'FAST_CalcSteady'
   
    
      ErrStat = ErrID_None
      ErrMsg  = ""


         ! get azimuth angle

      psi = ED%y%LSSTipPxa
      call Zero2TwoPi( psi )

      if (n_t_global == 0) then
            ! initialize a few things on the first call:
         call FAST_InitSteadyOutputs( psi, p_FAST, m_FAST, ED, BD, SrvD, AD, IfW, HD, SD, ExtPtfm, MAPp, FEAM, MD, Orca, &
                   IceF, IceD, ErrStat2, ErrMsg2 )
            call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      else
         DeltaAzim =  psi - m_FAST%Lin%Psi(1)
         call Zero2TwoPi(DeltaAzim)
      
         if (DeltaAzim > p_FAST%AzimDelta) then
            call SetErrStat(ErrID_Fatal, "The rotor is spinning too fast. The time step or NLinTimes is too large when CalcSteady=true.", ErrStat, ErrMsg, RoutineName)
            return
         end if
         
            ! save the outputs and azimuth angle for possible interpolation later
         call FAST_SaveOutputs( psi, p_FAST, m_FAST, ED, BD, SrvD, AD, IfW, HD, SD, ExtPtfm, MAPp, FEAM, MD, Orca, &
                      IceF, IceD, ErrStat2, ErrMsg2 )
            call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      end if
      if (ErrStat >= AbortErrLev) return

      
      
      if ( m_FAST%Lin%AzimTarget(m_FAST%Lin%AzimIndx-1) <= m_FAST%Lin%AzimTarget(m_FAST%Lin%AzimIndx) ) then ! the equal sign takes care of the zero-rpm case
         NextAzimuth = psi >= m_FAST%Lin%AzimTarget(m_FAST%Lin%AzimIndx)
      else
         ! this is the 2pi boundary, so we are either larger than the last target azimuth or less than the next one
         NextAzimuth = psi >= m_FAST%Lin%AzimTarget(m_FAST%Lin%AzimIndx) .and. psi < m_FAST%Lin%AzimTarget(m_FAST%Lin%AzimIndx-1)
      end if

      ! Forcing linearization if it's the last step
      if (t_global >= p_FAST%TMax - 0.5_DbKi*p_FAST%DT) then
         call WrScr('')
         call WrScr('[WARNING] Steady state not found before end of simulation. Forcing linearization.')
         m_FAST%Lin%ForceLin = .True.
         m_FAST%Lin%AzimIndx = 1
         NextAzimuth         = .True.
      endif
      
      if (NextAzimuth) then
      
            ! interpolate to find y at the target azimuth
         call FAST_DiffInterpOutputs( m_FAST%Lin%AzimTarget(m_FAST%Lin%AzimIndx), p_FAST, y_FAST, m_FAST, ED, BD, SrvD, AD, IfW, HD, SD, ExtPtfm, MAPp, FEAM, MD, Orca, &
                   IceF, IceD, ErrStat, ErrMsg )
         ! If linearization is forced
         if (m_FAST%Lin%ForceLin) then
            m_FAST%Lin%IsConverged = .True.
         endif
                   
         if (m_FAST%Lin%IsConverged .or. m_FAST%Lin%n_rot == 0) then ! save this operating point for linearization later
            m_FAST%Lin%LinTimes(m_FAST%Lin%AzimIndx) = t_global  
            call SaveOP(m_FAST%Lin%AzimIndx, p_FAST, y_FAST, ED, BD, SrvD, AD, IfW, OpFM, HD, SD, ExtPtfm, MAPp, FEAM, MD, Orca, &
                                  IceF, IceD, ErrStat, ErrMsg, m_FAST%Lin%CopyOP_CtrlCode )
         end if
         
             ! increment the counter to check the next azimuth:
         m_FAST%Lin%AzimIndx = m_FAST%Lin%AzimIndx + 1
         
             ! check if we've completed one rotor revolution
         if (m_FAST%Lin%AzimIndx > p_FAST%NLinTimes) then
            m_FAST%Lin%n_rot = m_FAST%Lin%n_rot + 1
         
            m_FAST%Lin%FoundSteady = m_FAST%Lin%IsConverged
            
            if (.not. m_FAST%Lin%FoundSteady) then
               ! compute the reference values for this rotor revolution
               call ComputeOutputRanges(p_FAST, y_FAST, m_FAST, SrvD%y)
               m_FAST%Lin%IsConverged = .true. ! check errors next rotor revolution
               m_FAST%Lin%AzimIndx = 1
               m_FAST%Lin%CopyOP_CtrlCode = MESH_UPDATECOPY
            end if
            ! Forcing linearization if time is close to tmax (with sufficient margin) 
            if (.not.m_FAST%Lin%FoundSteady) then
               if (ED%p%RotSpeed>0) then
                  ! If simulation is at least 10 revolutions, and error in rotor speed less than 0.1%
                  if ((p_FAST%TMax>10*(TwoPi_D)/ED%p%RotSpeed) .and. ( t_global >= p_FAST%TMax - 2._DbKi*(TwoPi_D)/ED%p%RotSpeed)) then
                     if (abs(ED%y%RotSpeed-ED%p%RotSpeed)/ED%p%RotSpeed<0.001) then
                        call WrScr('')
                        call WrScr('[WARNING] Steady state not found before end of simulation. Forcing linearization.')
                        m_FAST%Lin%ForceLin = .True. 
                     endif
                  endif
               else
                  if (t_global >= p_FAST%TMax - 1.5_DbKi*p_FAST%DT) then
                     call WrScr('')
                     call WrScr('[WARNING] Steady state not found before end of simulation. Forcing linearization.')
                     m_FAST%Lin%ForceLin = .True. 
                  endif
               endif
            endif
         end if
         
      end if
      if (m_FAST%Lin%ForceLin) then
         m_FAST%Lin%IsConverged=.true.
         m_FAST%Lin%FoundSteady=.true.
      endif
         

END SUBROUTINE FAST_CalcSteady
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine initializes variables for calculating periodic steady-state solution.
SUBROUTINE FAST_InitSteadyOutputs( psi, p_FAST, m_FAST, ED, BD, SrvD, AD, IfW, HD, SD, ExtPtfm, MAPp, FEAM, MD, Orca, &
                   IceF, IceD, ErrStat, ErrMsg )

   REAL(DbKi),               INTENT(IN   ) :: psi                 !< psi (rotor azimuth) at which the outputs are defined
   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_MiscVarType),   INTENT(INOUT) :: m_FAST              !< Miscellaneous variables
     
   TYPE(ElastoDyn_Data),     INTENT(INOUT) :: ED                  !< ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(INOUT) :: BD                  !< BeamDyn data
   TYPE(ServoDyn_Data),      INTENT(INOUT) :: SrvD                !< ServoDyn data
   TYPE(AeroDyn_Data),       INTENT(INOUT) :: AD                  !< AeroDyn data
   TYPE(InflowWind_Data),    INTENT(INOUT) :: IfW                 !< InflowWind data
   TYPE(HydroDyn_Data),      INTENT(INOUT) :: HD                  !< HydroDyn data
   TYPE(SubDyn_Data),        INTENT(INOUT) :: SD                  !< SubDyn data
   TYPE(ExtPtfm_Data),       INTENT(INOUT) :: ExtPtfm             !< ExtPtfm data
   TYPE(MAP_Data),           INTENT(INOUT) :: MAPp                !< MAP data
   TYPE(FEAMooring_Data),    INTENT(INOUT) :: FEAM                !< FEAMooring data
   TYPE(MoorDyn_Data),       INTENT(INOUT) :: MD                  !< Data for the MoorDyn module
   TYPE(OrcaFlex_Data),      INTENT(INOUT) :: Orca                !< OrcaFlex interface data
   TYPE(IceFloe_Data),       INTENT(INOUT) :: IceF                !< IceFloe data
   TYPE(IceDyn_Data),        INTENT(INOUT) :: IceD                !< All the IceDyn data used in time-step loop

   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None

   ! local variables
   INTEGER(IntKi)                          :: j, k                ! loop counters
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMsg2
   
   CHARACTER(*), PARAMETER                 :: RoutineName = 'FAST_InitSteadyOutputs'
   
    
      ErrStat = ErrID_None
      ErrMsg  = ""
      
      do j=1,p_FAST%NLinTimes
         m_FAST%Lin%AzimTarget(j) = (j-1) * p_FAST%AzimDelta + psi
         call Zero2TwoPi( m_FAST%Lin%AzimTarget(j) )
      end do
      ! this is circular, so I am going to add points at the beginning and end to avoid 
      ! more IF statements later
      m_FAST%Lin%AzimTarget(0) = m_FAST%Lin%AzimTarget(p_FAST%NLinTimes)
      m_FAST%Lin%AzimTarget(p_FAST%NLinTimes+1) = m_FAST%Lin%AzimTarget(1)
      

         ! Azimuth angles that correspond to Output arrays for interpolation:
      !m_FAST%Lin%Psi  = psi ! initialize entire array (note that we won't be able to interpolate with a constant array
      DO j = 1, p_FAST%LinInterpOrder + 1
         m_FAST%Lin%Psi(j) = psi - (j - 1) * D2R_D  ! arbitrarily say azimuth is one degree different
      END DO
      
      
         ! ElastoDyn
         allocate( ED%Output( p_FAST%LinInterpOrder+1 ), STAT = ErrStat2 )
         if (ErrStat2 /= 0) then
            call SetErrStat(ErrID_Fatal, "Error allocating ED%Output.", ErrStat, ErrMsg, RoutineName )
         else
            do j = 1, p_FAST%LinInterpOrder + 1
               call ED_CopyOutput(ED%y, ED%Output(j), MESH_NEWCOPY, ErrStat2, ErrMsg2)
                  call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            end do
            
            call ED_CopyOutput(ED%y, ED%y_interp, MESH_NEWCOPY, ErrStat2, ErrMsg2)
               call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
         end if
      
         ! BeamDyn
         IF (p_FAST%CompElast == Module_BD) THEN
         
            allocate( BD%Output( p_FAST%LinInterpOrder+1, p_FAST%nBeams ), STAT = ErrStat2 )
            if (ErrStat2 /= 0) then
               call SetErrStat(ErrID_Fatal, "Error allocating BD%Output.", ErrStat, ErrMsg, RoutineName )
            else
               do k=1,p_FAST%nBeams
                  do j = 1, p_FAST%LinInterpOrder + 1
                     call BD_CopyOutput(BD%y(k), BD%Output(j,k), MESH_NEWCOPY, ErrStat2, ErrMsg2)
                        call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
                  end do
               end do

               allocate( BD%y_interp( p_FAST%nBeams ), STAT = ErrStat2 )
               if (ErrStat2 /= 0) then
                  call SetErrStat(ErrID_Fatal, "Error allocating BD%Output.", ErrStat, ErrMsg, RoutineName )
               else
                  do k=1,p_FAST%nBeams
                     call BD_CopyOutput(BD%y(k), BD%y_interp(k), MESH_NEWCOPY, ErrStat2, ErrMsg2)
                        call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
                  end do
               end if
               
            end if
            
         END IF  ! BeamDyn 
         
      ! AeroDyn
      IF ( p_FAST%CompAero == Module_AD ) THEN
         
         allocate( AD%Output( p_FAST%LinInterpOrder+1 ), STAT = ErrStat2 )
         if (ErrStat2 /= 0) then
            call SetErrStat(ErrID_Fatal, "Error allocating AD%Output.", ErrStat, ErrMsg, RoutineName )
         else
            do j = 1, p_FAST%LinInterpOrder + 1
               call AD_CopyOutput(AD%y, AD%Output(j), MESH_NEWCOPY, ErrStat2, ErrMsg2)
                  call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            end do
            
            call AD_CopyOutput(AD%y, AD%y_interp, MESH_NEWCOPY, ErrStat2, ErrMsg2)
               call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
         end if
         
      END IF  ! CompAero
      
         
      ! InflowWind
      IF ( p_FAST%CompInflow == Module_IfW ) THEN
         
         allocate( IfW%Output( p_FAST%LinInterpOrder+1 ), STAT = ErrStat2 )
         if (ErrStat2 /= 0) then
            call SetErrStat(ErrID_Fatal, "Error allocating IfW%Output.", ErrStat, ErrMsg, RoutineName )
         else
            do j = 1, p_FAST%LinInterpOrder + 1
               call InflowWind_CopyOutput(IfW%y, IfW%Output(j), MESH_NEWCOPY, ErrStat2, ErrMsg2)
                  call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            end do
            
            call InflowWind_CopyOutput(IfW%y, IfW%y_interp, MESH_NEWCOPY, ErrStat2, ErrMsg2)
               call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
         end if
            
      END IF  ! CompInflow
      
      
      ! ServoDyn
      IF ( p_FAST%CompServo == Module_SrvD ) THEN
         
         allocate( SrvD%Output( p_FAST%LinInterpOrder+1 ), STAT = ErrStat2 )
         if (ErrStat2 /= 0) then
            call SetErrStat(ErrID_Fatal, "Error allocating SrvD%Output.", ErrStat, ErrMsg, RoutineName )
         else
            do j = 1, p_FAST%LinInterpOrder + 1
               call SrvD_CopyOutput(SrvD%y, SrvD%Output(j), MESH_NEWCOPY, ErrStat2, ErrMsg2)
                  call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            end do
            
            call SrvD_CopyOutput(SrvD%y, SrvD%y_interp, MESH_NEWCOPY, ErrStat2, ErrMsg2)
               call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
         end if
            
      END IF  ! ServoDyn
      
      ! HydroDyn
      IF ( p_FAST%CompHydro == Module_HD ) THEN

         allocate( HD%Output( p_FAST%LinInterpOrder+1 ), STAT = ErrStat2 )
         if (ErrStat2 /= 0) then
            call SetErrStat(ErrID_Fatal, "Error allocating HD%Output.", ErrStat, ErrMsg, RoutineName )
         else
            do j = 1, p_FAST%LinInterpOrder + 1
               call HydroDyn_CopyOutput(HD%y, HD%Output(j), MESH_NEWCOPY, ErrStat2, ErrMsg2)
                  call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            end do
            
            call HydroDyn_CopyOutput(HD%y, HD%y_interp, MESH_NEWCOPY, ErrStat2, ErrMsg2)
               call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
         end if
            
      END IF  ! HydroDyn
      
      !! SubDyn/ExtPtfm_MCKF
      IF ( p_FAST%CompSub == Module_SD ) THEN
         allocate( SD%Output( p_FAST%LinInterpOrder+1 ), STAT = ErrStat2 )
         if (ErrStat2 /= 0) then
            call SetErrStat(ErrID_Fatal, "Error allocating SD%Output.", ErrStat, ErrMsg, RoutineName )
         else
            do j = 1, p_FAST%LinInterpOrder + 1
               call SD_CopyOutput(SD%y, SD%Output(j), MESH_NEWCOPY, ErrStat2, ErrMsg2)
                  call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            end do
            
            call SD_CopyOutput(SD%y, SD%y_interp, MESH_NEWCOPY, ErrStat2, ErrMsg2)
               call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
         end if   
      ELSE IF ( p_FAST%CompSub == Module_ExtPtfm ) THEN
      END IF  ! SubDyn/ExtPtfm_MCKF
      
      
      ! Mooring (MAP , FEAM , MoorDyn)
      ! MAP
      IF ( p_FAST%CompMooring == Module_MAP ) THEN
         
         allocate( MAPp%Output( p_FAST%LinInterpOrder+1 ), STAT = ErrStat2 )
         if (ErrStat2 /= 0) then
            call SetErrStat(ErrID_Fatal, "Error allocating MAPp%Output.", ErrStat, ErrMsg, RoutineName )
         else
            do j = 1, p_FAST%LinInterpOrder + 1
               call MAP_CopyOutput(MAPp%y, MAPp%Output(j), MESH_NEWCOPY, ErrStat2, ErrMsg2)
                  call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            end do
            
            call MAP_CopyOutput(MAPp%y, MAPp%y_interp, MESH_NEWCOPY, ErrStat2, ErrMsg2)
               call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
         end if
            
      ! MoorDyn
      ELSEIF ( p_FAST%CompMooring == Module_MD ) THEN
      
         allocate( MD%Output( p_FAST%LinInterpOrder+1 ), STAT = ErrStat2 )
         if (ErrStat2 /= 0) then
            call SetErrStat(ErrID_Fatal, "Error allocating MD%Output.", ErrStat, ErrMsg, RoutineName )
         else
            do j = 1, p_FAST%LinInterpOrder + 1
               call MD_CopyOutput(MD%y, MD%Output(j), MESH_NEWCOPY, ErrStat2, ErrMsg2)
                  call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            end do
            
            call MD_CopyOutput(MD%y, MD%y_interp, MESH_NEWCOPY, ErrStat2, ErrMsg2)
               call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
         end if
      
      
      
      
      !! FEAM
      !ELSEIF ( p_FAST%CompMooring == Module_FEAM ) THEN
      !! OrcaFlex
      !ELSEIF ( p_FAST%CompMooring == Module_Orca ) THEN
         
      END IF  ! MAP/FEAM/MoorDyn/OrcaFlex
      
           
            
      !! Ice (IceFloe or IceDyn)
      !! IceFloe
      !IF ( p_FAST%CompIce == Module_IceF ) THEN
      !   
      !! IceDyn
      !ELSEIF ( p_FAST%CompIce == Module_IceD ) THEN
      !
      !END IF  ! IceFloe/IceDyn


END SUBROUTINE FAST_InitSteadyOutputs
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine saves outputs for future interpolation at a desired azimuth.
SUBROUTINE FAST_SaveOutputs( psi, p_FAST, m_FAST, ED, BD, SrvD, AD, IfW, HD, SD, ExtPtfm, MAPp, FEAM, MD, Orca, &
                   IceF, IceD, ErrStat, ErrMsg )

   REAL(DbKi),               INTENT(IN   ) :: psi                 !< psi (rotor azimuth) at which the outputs are defined
   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_MiscVarType),   INTENT(INOUT) :: m_FAST              !< Miscellaneous variables
     
   TYPE(ElastoDyn_Data),     INTENT(INOUT) :: ED                  !< ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(INOUT) :: BD                  !< BeamDyn data
   TYPE(ServoDyn_Data),      INTENT(INOUT) :: SrvD                !< ServoDyn data
   TYPE(AeroDyn_Data),       INTENT(INOUT) :: AD                  !< AeroDyn data
   TYPE(InflowWind_Data),    INTENT(INOUT) :: IfW                 !< InflowWind data
   TYPE(HydroDyn_Data),      INTENT(INOUT) :: HD                  !< HydroDyn data
   TYPE(SubDyn_Data),        INTENT(INOUT) :: SD                  !< SubDyn data
   TYPE(ExtPtfm_Data),       INTENT(INOUT) :: ExtPtfm             !< ExtPtfm data
   TYPE(MAP_Data),           INTENT(INOUT) :: MAPp                !< MAP data
   TYPE(FEAMooring_Data),    INTENT(INOUT) :: FEAM                !< FEAMooring data
   TYPE(MoorDyn_Data),       INTENT(INOUT) :: MD                  !< Data for the MoorDyn module
   TYPE(OrcaFlex_Data),      INTENT(INOUT) :: Orca                !< OrcaFlex interface data
   TYPE(IceFloe_Data),       INTENT(INOUT) :: IceF                !< IceFloe data
   TYPE(IceDyn_Data),        INTENT(INOUT) :: IceD                !< All the IceDyn data used in time-step loop

   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None

   ! local variables
   INTEGER(IntKi)                          :: j, k                ! loop counters
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMsg2
   
   CHARACTER(*), PARAMETER                 :: RoutineName = 'FAST_SaveOutputs'
   
    
      ErrStat = ErrID_None
      ErrMsg  = ""
      
      DO j = p_FAST%LinInterpOrder, 1, -1
         m_FAST%Lin%Psi(j+1) = m_FAST%Lin%Psi(j)
      END DO
      
      if (psi < m_FAST%Lin%Psi(1)) then
         ! if we go around a 2pi boundary, we will subtract 2pi from the saved values so that interpolation works as expected
         m_FAST%Lin%Psi = m_FAST%Lin%Psi - TwoPi_D
      end if
      m_FAST%Lin%Psi(1)  = psi

         ! ElastoDyn
         DO j = p_FAST%LinInterpOrder, 1, -1
            CALL ED_CopyOutput(ED%Output(j), ED%Output(j+1), MESH_UPDATECOPY, Errstat2, ErrMsg2)
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
         END DO
  
         CALL ED_CopyOutput (ED%y,  ED%Output(1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
      
         ! BeamDyn
         IF (p_FAST%CompElast == Module_BD) THEN
         
            DO k = 1,p_FAST%nBeams
         
               DO j = p_FAST%LinInterpOrder, 1, -1
                  CALL BD_CopyOutput (BD%Output(j,k),  BD%Output(j+1,k),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
                     CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
               END DO
  
               CALL BD_CopyOutput (BD%y(k),  BD%Output(1,k),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
                  CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            
            END DO ! k=p_FAST%nBeams
         
         END IF  ! BeamDyn 
         
      
      ! AeroDyn
      IF ( p_FAST%CompAero == Module_AD ) THEN
         
         DO j = p_FAST%LinInterpOrder, 1, -1
            CALL AD_CopyOutput (AD%Output(j),  AD%Output(j+1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
         END DO
  
         CALL AD_CopyOutput (AD%y,  AD%Output(1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
         
      END IF  ! CompAero
      
         
      ! InflowWind
      IF ( p_FAST%CompInflow == Module_IfW ) THEN
         
         DO j = p_FAST%LinInterpOrder, 1, -1
            CALL InflowWind_CopyOutput (IfW%Output(j),  IfW%Output(j+1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
         END DO
  
         CALL InflowWind_CopyOutput (IfW%y,  IfW%Output(1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            
      END IF  ! CompInflow
      
      
      ! ServoDyn
      IF ( p_FAST%CompServo == Module_SrvD ) THEN
         
         DO j = p_FAST%LinInterpOrder, 1, -1
            CALL SrvD_CopyOutput (SrvD%Output(j),  SrvD%Output(j+1),  MESH_UPDATECOPY, ErrStat2, ErrMsg2)
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
         END DO
  
         CALL SrvD_CopyOutput (SrvD%y,  SrvD%Output(1),  MESH_UPDATECOPY, ErrStat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            
      END IF  ! ServoDyn       
      
      ! HydroDyn
      IF ( p_FAST%CompHydro == Module_HD ) THEN

         DO j = p_FAST%LinInterpOrder, 1, -1

            CALL HydroDyn_CopyOutput (HD%Output(j),  HD%Output(j+1),  MESH_UPDATECOPY, ErrStat2, ErrMsg2)
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
         END DO

         CALL HydroDyn_CopyOutput (HD%y,  HD%Output(1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            
      END IF  ! HydroDyn

      !! SubDyn/ExtPtfm_MCKF
      IF ( p_FAST%CompSub == Module_SD ) THEN
         DO j = p_FAST%LinInterpOrder, 1, -1
            CALL SD_CopyOutput (SD%Output(j),  SD%Output(j+1),  MESH_UPDATECOPY, ErrStat2, ErrMsg2)
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
         END DO

         CALL SD_CopyOutput (SD%y,  SD%Output(1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
      ELSE IF ( p_FAST%CompSub == Module_ExtPtfm ) THEN
      END IF  ! SubDyn/ExtPtfm_MCKF
      
      
      ! Mooring (MAP , FEAM , MoorDyn)
      ! MAP
      IF ( p_FAST%CompMooring == Module_MAP ) THEN
         
         DO j = p_FAST%LinInterpOrder, 1, -1
            CALL MAP_CopyOutput (MAPp%Output(j),  MAPp%Output(j+1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
         END DO
  
         CALL MAP_CopyOutput (MAPp%y,  MAPp%Output(1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            
      ! MoorDyn
      ELSEIF ( p_FAST%CompMooring == Module_MD ) THEN
      
         DO j = p_FAST%LinInterpOrder, 1, -1
            CALL MD_CopyOutput (MD%Output(j),  MD%Output(j+1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
         END DO
  
         CALL MD_CopyOutput (MD%y,  MD%Output(1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            
      !! FEAM
      !ELSEIF ( p_FAST%CompMooring == Module_FEAM ) THEN
      !! OrcaFlex
      !ELSEIF ( p_FAST%CompMooring == Module_Orca ) THEN
         
      END IF  ! MAP/FEAM/MoorDyn/OrcaFlex
      
           
            
      !! Ice (IceFloe or IceDyn)
      !! IceFloe
      !IF ( p_FAST%CompIce == Module_IceF ) THEN
      !   
      !! IceDyn
      !ELSEIF ( p_FAST%CompIce == Module_IceD ) THEN
      !
      !END IF  ! IceFloe/IceDyn


END SUBROUTINE FAST_SaveOutputs
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine interpolates the outputs at the target azimuths, computes the compared to the previous rotation, and stores 
!! them for future rotation .
SUBROUTINE FAST_DiffInterpOutputs( psi_target, p_FAST, y_FAST, m_FAST, ED, BD, SrvD, AD, IfW, HD, SD, ExtPtfm, MAPp, FEAM, MD, Orca, &
                   IceF, IceD, ErrStat, ErrMsg )

   REAL(DbKi),               INTENT(IN   ) :: psi_target          !< psi (rotor azimuth) at which the outputs are requested
   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST              !< Output variables for the glue code
   TYPE(FAST_MiscVarType),   INTENT(INOUT) :: m_FAST              !< Miscellaneous variables
     
   TYPE(ElastoDyn_Data),     INTENT(INOUT) :: ED                  !< ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(INOUT) :: BD                  !< BeamDyn data
   TYPE(ServoDyn_Data),      INTENT(INOUT) :: SrvD                !< ServoDyn data
   TYPE(AeroDyn_Data),       INTENT(INOUT) :: AD                  !< AeroDyn data
   TYPE(InflowWind_Data),    INTENT(INOUT) :: IfW                 !< InflowWind data
   TYPE(HydroDyn_Data),      INTENT(INOUT) :: HD                  !< HydroDyn data
   TYPE(SubDyn_Data),        INTENT(INOUT) :: SD                  !< SubDyn data
   TYPE(ExtPtfm_Data),       INTENT(INOUT) :: ExtPtfm             !< ExtPtfm data
   TYPE(MAP_Data),           INTENT(INOUT) :: MAPp                !< MAP data
   TYPE(FEAMooring_Data),    INTENT(INOUT) :: FEAM                !< FEAMooring data
   TYPE(MoorDyn_Data),       INTENT(INOUT) :: MD                  !< Data for the MoorDyn module
   TYPE(OrcaFlex_Data),      INTENT(INOUT) :: Orca                !< OrcaFlex interface data
   TYPE(IceFloe_Data),       INTENT(INOUT) :: IceF                !< IceFloe data
   TYPE(IceDyn_Data),        INTENT(INOUT) :: IceD                !< All the IceDyn data used in time-step loop

   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None

   ! local variables
   INTEGER(IntKi)                          :: k                   ! loop counters
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMsg2
   REAL(DbKi)                              :: t_global
   REAL(ReKi)                              :: eps_squared
   
   CHARACTER(*), PARAMETER                 :: RoutineName = 'FAST_DiffInterpOutputs'
    
      ErrStat = ErrID_None
      ErrMsg  = ""
      t_global = 0.0_DbKi ! we don't really need this to get the output OPs

      !................................................................................................
      ! Extrapolate outputs to the target azimuth and pack into OP arrays
      !................................................................................................

         ! ElastoDyn
         CALL ED_Output_ExtrapInterp (ED%Output, m_FAST%Lin%Psi,  ED%y_interp, psi_target, ErrStat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
      
         call ED_GetOP( t_global, ED%Input(1), ED%p, ED%x(STATE_CURR), ED%xd(STATE_CURR), ED%z(STATE_CURR), ED%OtherSt(STATE_CURR), &
                           ED%y_interp, ED%m, ErrStat2, ErrMsg2, y_op=y_FAST%Lin%Modules(Module_ED)%Instance(1)%op_y, NeedTrimOP=.true.)
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            
         ! BeamDyn
         IF (p_FAST%CompElast == Module_BD) THEN
         
            DO k = 1,p_FAST%nBeams
         
               CALL BD_Output_ExtrapInterp (BD%Output(:,k), m_FAST%Lin%Psi,  BD%y_interp(k), psi_target, ErrStat2, ErrMsg2)
                  CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            
               call BD_GetOP( t_global, BD%Input(1,k), BD%p(k), BD%x(k,STATE_CURR), BD%xd(k,STATE_CURR), BD%z(k,STATE_CURR), BD%OtherSt(k,STATE_CURR), &
                           BD%y_interp(k), BD%m(k), ErrStat2, ErrMsg2, y_op=y_FAST%Lin%Modules(Module_BD)%Instance(k)%op_y, NeedTrimOP=.true.)
                  call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            END DO ! k=p_FAST%nBeams
         
         END IF  ! BeamDyn 
         
      
      ! AeroDyn
      IF ( p_FAST%CompAero == Module_AD ) THEN
         
         CALL AD_Output_ExtrapInterp (AD%Output, m_FAST%Lin%Psi,  AD%y_interp, psi_target, ErrStat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
         
         call AD_GetOP( t_global, AD%Input(1), AD%p, AD%x(STATE_CURR), AD%xd(STATE_CURR), AD%z(STATE_CURR), AD%OtherSt(STATE_CURR), &
                           AD%y_interp, AD%m, ErrStat2, ErrMsg2, y_op=y_FAST%Lin%Modules(Module_AD)%Instance(1)%op_y)
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      END IF  ! CompAero
      
         
      ! InflowWind
      IF ( p_FAST%CompInflow == Module_IfW ) THEN
         
         CALL InflowWind_Output_ExtrapInterp (IfW%Output, m_FAST%Lin%Psi,  IfW%y_interp, psi_target, ErrStat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            
         call InflowWind_GetOP( t_global, IfW%Input(1), IfW%p, IfW%x(STATE_CURR), IfW%xd(STATE_CURR), IfW%z(STATE_CURR), IfW%OtherSt(STATE_CURR), &
                           IfW%y_interp, IfW%m, ErrStat2, ErrMsg2, y_op=y_FAST%Lin%Modules(Module_IfW)%Instance(1)%op_y)
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      END IF  ! CompInflow
      
      
      ! ServoDyn
      IF ( p_FAST%CompServo == Module_SrvD ) THEN
         
         CALL SrvD_Output_ExtrapInterp (SrvD%Output, m_FAST%Lin%Psi,  SrvD%y_interp, psi_target, ErrStat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            
         call SrvD_GetOP( t_global, SrvD%Input(1), SrvD%p, SrvD%x(STATE_CURR), SrvD%xd(STATE_CURR), SrvD%z(STATE_CURR), SrvD%OtherSt(STATE_CURR), &
                           SrvD%y_interp, SrvD%m, ErrStat2, ErrMsg2, y_op=y_FAST%Lin%Modules(Module_SrvD)%Instance(1)%op_y)
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      END IF  ! ServoDyn
      
      ! HydroDyn
      IF ( p_FAST%CompHydro == Module_HD ) THEN

         CALL HydroDyn_Output_ExtrapInterp (HD%Output, m_FAST%Lin%Psi,  HD%y_interp, psi_target, ErrStat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            
         call HD_GetOP( t_global, HD%Input(1), HD%p, HD%x(STATE_CURR), HD%xd(STATE_CURR), HD%z(STATE_CURR), HD%OtherSt(STATE_CURR), &
                           HD%y_interp, HD%m, ErrStat2, ErrMsg2, y_op=y_FAST%Lin%Modules(Module_HD)%Instance(1)%op_y)
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      END IF  ! HydroDyn

  
      !! SubDyn/ExtPtfm_MCKF
      IF ( p_FAST%CompSub == Module_SD ) THEN
      
         CALL SD_Output_ExtrapInterp (SD%Output, m_FAST%Lin%Psi,  SD%y_interp, psi_target, ErrStat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            
         call SD_GetOP( t_global, SD%Input(1), SD%p, SD%x(STATE_CURR), SD%xd(STATE_CURR), SD%z(STATE_CURR), SD%OtherSt(STATE_CURR), &
                           SD%y_interp, SD%m, ErrStat2, ErrMsg2, y_op=y_FAST%Lin%Modules(Module_SD)%Instance(1)%op_y, NeedTrimOP=.true.)
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      ELSE IF ( p_FAST%CompSub == Module_ExtPtfm ) THEN
      END IF  ! SubDyn/ExtPtfm_MCKF
      
      
      ! Mooring (MAP , FEAM , MoorDyn)
      ! MAP
      IF ( p_FAST%CompMooring == Module_MAP ) THEN
         
         CALL MAP_Output_ExtrapInterp (MAPp%Output, m_FAST%Lin%Psi,  MAPp%y_interp, psi_target, ErrStat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            
         call MAP_GetOP( t_global, MAPp%Input(1), MAPp%p, MAPp%x(STATE_CURR), MAPp%xd(STATE_CURR), MAPp%z(STATE_CURR), MAPp%OtherSt, &
                           MAPp%y_interp, ErrStat2, ErrMsg2, y_op=y_FAST%Lin%Modules(Module_MAP)%Instance(1)%op_y)
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      ! MoorDyn
      ELSEIF ( p_FAST%CompMooring == Module_MD ) THEN
      
         CALL MD_Output_ExtrapInterp (MD%Output, m_FAST%Lin%Psi,  MD%y_interp, psi_target, ErrStat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            
         call MD_GetOP( t_global, MD%Input(1), MD%p, MD%x(STATE_CURR), MD%xd(STATE_CURR), MD%z(STATE_CURR), MD%OtherSt(STATE_CURR), &
                           MD%y_interp, MD%m, ErrStat2, ErrMsg2, y_op=y_FAST%Lin%Modules(Module_MD)%Instance(1)%op_y)
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
      !! FEAM
      !ELSEIF ( p_FAST%CompMooring == Module_FEAM ) THEN
      !! OrcaFlex
      !ELSEIF ( p_FAST%CompMooring == Module_Orca ) THEN
         
      END IF  ! MAP/FEAM/MoorDyn/OrcaFlex
      
           
            
      !! Ice (IceFloe or IceDyn)
      !! IceFloe
      !IF ( p_FAST%CompIce == Module_IceF ) THEN
      !   
      !! IceDyn
      !ELSEIF ( p_FAST%CompIce == Module_IceD ) THEN
      !
      !END IF  ! IceFloe/IceDyn

      
      call pack_in_array(p_FAST, y_FAST, m_FAST)
      
      if (m_FAST%Lin%IsConverged) then
         ! check that error equation is less than TrimTol !!!call 
         call calc_error(p_FAST, y_FAST, m_FAST, SrvD%y, eps_squared)
         m_FAST%Lin%IsConverged = eps_squared < p_FAST%TrimTol
      end if
      
      
      m_FAST%Lin%Y_prevRot(:,m_FAST%Lin%AzimIndx) = m_FAST%Lin%y_interp
      
END SUBROUTINE FAST_DiffInterpOutputs
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE pack_in_array(p_FAST, y_FAST, m_FAST)

   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_MiscVarType),   INTENT(INOUT) :: m_FAST              !< Miscellaneous variables
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST              !< Output variables for the glue code
   
   INTEGER(IntKi)                          :: ThisModule          !< module identifier
   INTEGER(IntKi)                          :: ThisInstance        !< index of the module instance

   integer                                 :: i, j
   integer                                 :: ny
   integer                                 :: indx
   
   ! note that op_y may be larger than SizeLin if there are orientations; also, we are NOT including the WriteOutputs

   do i = 1,p_FAST%Lin_NumMods
      ThisModule = p_FAST%Lin_ModOrder( i )
         
      do ThisInstance=1,size(y_FAST%Lin%Modules(ThisModule)%Instance)
      
         ny = y_FAST%Lin%Modules(ThisModule)%Instance(ThisInstance)%SizeLin(LIN_OUTPUT_COL) - y_FAST%Lin%Modules(ThisModule)%Instance(ThisInstance)%NumOutputs !last column before WriteOutput occurs
         do j=1,ny
            indx = y_FAST%Lin%Modules(ThisModule)%Instance(ThisInstance)%LinStartIndx(LIN_OUTPUT_COL) + j - 1
            
            m_FAST%Lin%y_interp( indx ) = y_FAST%Lin%Modules(ThisModule)%Instance(ThisInstance)%op_y(j)
         end do
         
      end do
   end do
   
END SUBROUTINE pack_in_array
!----------------------------------------------------------------------------------------------------------------------------------
!> This function computes the error function between this rotor revolution and the previous one.
!! Angles represented in m_FAST%Lin%y_interp may have 2pi added or subtracted to allow the angles to be closer to the previous
!! rotor revolution.
SUBROUTINE calc_error(p_FAST, y_FAST, m_FAST, y_SrvD, eps_squared)

   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_MiscVarType),   INTENT(INOUT) :: m_FAST              !< Miscellaneous variables
   TYPE(FAST_OutputFileType),INTENT(IN   ) :: y_FAST              !< Output variables for the glue code
   TYPE(SrvD_OutputType),    INTENT(IN   ) :: y_SrvD              !< Output variables for the glue code
   REAL(ReKi)               ,INTENT(  OUT) :: eps_squared         !< epsilon squared
   
   INTEGER(IntKi)                          :: ThisModule          !< module identifier
   INTEGER(IntKi)                          :: ThisInstance        !< index of the module instance

   integer                                 :: i, j
   integer                                 :: ny
   integer                                 :: indx
   real(ReKi)                              :: diff
   
   
   ! special cases for angles:
      indx = Indx_y_Yaw_Start(y_FAST, Module_ED)  ! start of ED where Yaw, YawRate, HSS_Spd occur (right before WriteOutputs)
   call AddOrSub2Pi(m_FAST%Lin%Y_prevRot( indx, m_FAST%Lin%AzimIndx ), m_FAST%Lin%y_interp( indx ))

   if (p_FAST%CompServo == Module_SrvD) then
      do i = 1, size( y_SrvD%BlPitchCom )
         indx = y_FAST%Lin%Modules(Module_SrvD)%Instance(1)%LinStartIndx(LIN_OUTPUT_COL) + i - 1
         call AddOrSub2Pi(m_FAST%Lin%Y_prevRot( indx, m_FAST%Lin%AzimIndx ), m_FAST%Lin%y_interp( indx ))
      end do
   end if
   

   ! compute the error:
   eps_squared = 0.0_ReKi
   
   do i = 1,p_FAST%Lin_NumMods
      ThisModule = p_FAST%Lin_ModOrder( i )
         
      do ThisInstance=1,size(y_FAST%Lin%Modules(ThisModule)%Instance)
      
         ny = y_FAST%Lin%Modules(ThisModule)%Instance(ThisInstance)%SizeLin(LIN_OUTPUT_COL) - y_FAST%Lin%Modules(ThisModule)%Instance(ThisInstance)%NumOutputs !last column before WriteOutput occurs
         
         do j=1,ny
            indx = y_FAST%Lin%Modules(ThisModule)%Instance(ThisInstance)%LinStartIndx(LIN_OUTPUT_COL) + j - 1
            
            if (EqualRealNos(m_FAST%Lin%y_interp( indx ), m_FAST%Lin%Y_prevRot( indx, m_FAST%Lin%AzimIndx ))) then
               diff = 0.0_ReKi ! take care of some potential numerical issues
            else
               diff = m_FAST%Lin%y_interp( indx ) - m_FAST%Lin%Y_prevRot( indx, m_FAST%Lin%AzimIndx )
            end if
            
            eps_squared = eps_squared + ( diff / m_FAST%Lin%y_ref( indx ) ) ** 2
         end do
         
      end do
   end do
   

   !.................................
   ! Normalize:
   !.................................
   eps_squared = eps_squared / ( y_FAST%Lin%Glue%SizeLin(LIN_OUTPUT_COL) - y_FAST%Lin%Glue%NumOutputs )
   
!   write(50+m_FAST%Lin%AzimIndx,'(3000(F15.7,1x))') m_FAST%Lin%y_interp, eps_squared
END SUBROUTINE calc_error
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ComputeOutputRanges(p_FAST, y_FAST, m_FAST, y_SrvD)

   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_MiscVarType),   INTENT(INOUT) :: m_FAST              !< Miscellaneous variables
   TYPE(FAST_OutputFileType),INTENT(IN   ) :: y_FAST              !< Output variables for the glue code
   TYPE(SrvD_OutputType),    INTENT(IN   ) :: y_SrvD              !< Output variables for the glue code
   
   integer                                 :: indx
   integer                                 :: i
   
   ! note that op_y may be larger than SizeLin if there are orientations; also, we are NOT including the WriteOutputs

   do indx = 1,y_FAST%Lin%Glue%SizeLin(LIN_OUTPUT_COL)
      m_FAST%Lin%y_ref(indx) = maxval( m_FAST%Lin%Y_prevRot( indx, : ) ) - minval( m_FAST%Lin%Y_prevRot( indx, : ) )
      m_FAST%Lin%y_ref(indx) = max( m_FAST%Lin%y_ref(indx), 0.01_ReKi )
!      if (m_FAST%Lin%y_ref(indx) < 1.0e-4) m_FAST%Lin%y_ref(indx) = 1.0_ReKi ! not sure why we wouldn't just do m_FAST%Lin%y_ref(indx) = max(1.0_ReKi, m_FAST%Lin%y_ref(indx)) or max(1e-4, y_ref(indx))
   end do
   
   ! special case for angles:
      indx = Indx_y_Yaw_Start(y_FAST, Module_ED)  ! start of ED where Yaw, YawRate, HSS_Spd occur (right before WriteOutputs)
   m_FAST%Lin%y_ref(indx) = min( m_FAST%Lin%y_ref(indx), Pi )

   if (p_FAST%CompServo == Module_SrvD) then
      do i = 1, size( y_SrvD%BlPitchCom )
         indx = y_FAST%Lin%Modules(Module_SrvD)%Instance(1)%LinStartIndx(LIN_OUTPUT_COL) + i - 1
         m_FAST%Lin%y_ref(indx) = min( m_FAST%Lin%y_ref(indx), Pi )
      end do
   end if
   
   ! Note: I'm ignoring the periodicity of the log maps that represent orientations
   
END SUBROUTINE ComputeOutputRanges
!----------------------------------------------------------------------------------------------------------------------------------

END MODULE FAST_Linear
