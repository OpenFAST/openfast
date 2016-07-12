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
SUBROUTINE Init_Lin(p_FAST, y_FAST, m_FAST, ErrStat, ErrMsg)

   TYPE(FAST_ParameterType), INTENT(INOUT) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST              !< Output variables for the glue code
   TYPE(FAST_MiscVarType),   INTENT(INOUT) :: m_FAST              !< Miscellaneous variables
   
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None

   INTEGER(IntKi)                          :: position            ! position in string
   INTEGER(IntKi)                          :: i, j                ! loop/temp variables
   INTEGER(IntKi)                          :: ThisModule          ! Module ID # 
   
   INTEGER(IntKi)                          :: i_u, i_y, i_x       ! loop/temp variables

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
      
      
      ! I'm going to overwrite some of the input/output descriptions 
      if (p_FAST%CompServo == MODULE_SrvD) then
         do i=1,3
            position = index(y_FAST%Lin%Modules(Module_IfW)%Names_u(i), ',') - 1
            y_FAST%Lin%Modules(Module_IfW)%Names_u(i) = y_FAST%Lin%Modules(Module_IfW)%Names_u(i)(1:position)//' (hub)'//&
                                                         y_FAST%Lin%Modules(Module_IfW)%Names_u(i)(position+1:)
         end do                         
      end if
                  
   end if
   
      ! ServoDyn is next, if activated:
   if ( p_FAST%CompServo  == Module_SrvD ) then 
      p_FAST%Lin_NumMods = p_FAST%Lin_NumMods + 1
      p_FAST%Lin_ModOrder( p_FAST%Lin_NumMods ) = Module_SrvD
   end if
   
      ! ElastoDyn is next; it is always activated:
   p_FAST%Lin_NumMods = p_FAST%Lin_NumMods + 1
   p_FAST%Lin_ModOrder( p_FAST%Lin_NumMods ) = Module_ED
   
      ! AeroDyn is next, if activated:
   if ( p_FAST%CompAero  == Module_AD ) then 
      p_FAST%Lin_NumMods = p_FAST%Lin_NumMods + 1
      p_FAST%Lin_ModOrder( p_FAST%Lin_NumMods ) = Module_AD
   end if
   
   
   !.....................
   ! determine total number of inputs/outputs/contStates:
   !.....................
   p_FAST%SizeLin = 0
   do i = 1,p_FAST%Lin_NumMods
      ThisModule = p_FAST%Lin_ModOrder( i )
      p_FAST%SizeLin(ThisModule,LIN_INPUT_COL) = size(y_FAST%Lin%Modules(ThisModule)%Names_u)
      p_FAST%SizeLin(ThisModule,LIN_OUTPUT_COL) = size(y_FAST%Lin%Modules(ThisModule)%Names_y) !- y_FAST%numOuts(ThisModule)  
      p_FAST%SizeLin(ThisModule,LIN_ContSTATE_COL) = size(y_FAST%Lin%Modules(ThisModule)%Names_x)      
   end do
               
   !if (p_FAST%LinInputs == LIN_NONE) p_FAST%SizeLin(:,LIN_INPUT_COL) = 0
   !if (p_FAST%LinOutputs == LIN_NONE) p_FAST%SizeLin(:,LIN_OUTPUT_COL) = 0
         
   p_FAST%SizeLin(NumModules+1,LIN_INPUT_COL) = sum( p_FAST%SizeLin(1:NumModules,LIN_INPUT_COL) )  ! total number of inputs
   p_FAST%SizeLin(NumModules+1,LIN_OUTPUT_COL) = sum( p_FAST%SizeLin(1:NumModules,LIN_OUTPUT_COL) )  ! total number of outputs
   p_FAST%SizeLin(NumModules+1,LIN_ContSTATE_COL) = sum( p_FAST%SizeLin(1:NumModules,LIN_ContSTATE_COL) )  ! total number of continuous states
    
      ! ...................................
      ! get names of inputs, outputs, and continuous states
      ! ...................................
   call AllocAry( y_FAST%Lin%Glue%names_u, p_FAST%SizeLin(NumModules+1,LIN_INPUT_COL), 'names_u', ErrStat2, ErrMsg2)
      call SetErrStat(errStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call AllocAry( y_FAST%Lin%Glue%names_y, p_FAST%SizeLin(NumModules+1,LIN_OUTPUT_COL), 'names_y', ErrStat2, ErrMsg2)
      call SetErrStat(errStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)      
   call AllocAry( y_FAST%Lin%Glue%names_x, p_FAST%SizeLin(NumModules+1,LIN_ContSTATE_COL), 'names_x', ErrStat2, ErrMsg2)
      call SetErrStat(errStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)      
   if (ErrStat >= AbortErrLev) return
   
   
   i_u = 1
   i_y = 1      
   i_x = 1      
   do i = 1,p_FAST%Lin_NumMods
      ThisModule = p_FAST%Lin_ModOrder( i )

      do j=1,p_FAST%SizeLin(ThisModule,LIN_INPUT_COL)
         y_FAST%Lin%Glue%names_u(i_u) = TRIM(y_FAST%Module_Abrev(ThisModule))//' '//y_FAST%Lin%Modules(ThisModule)%Names_u(j)
         i_u = i_u + 1;
      end do

      do j=1,p_FAST%SizeLin(ThisModule,LIN_OUTPUT_COL)
         y_FAST%Lin%Glue%names_y(i_y) = TRIM(y_FAST%Module_Abrev(ThisModule))//' '//y_FAST%Lin%Modules(ThisModule)%Names_y(j)
         i_y = i_y + 1;
      end do      

      do j=1,p_FAST%SizeLin(ThisModule,LIN_ContSTATE_COL)
         y_FAST%Lin%Glue%names_x( i_x) = TRIM(y_FAST%Module_Abrev(ThisModule))//' '//y_FAST%Lin%Modules(ThisModule)%Names_x( j)
         i_x = i_x + 1;
      end do      
   end do
         
   
END SUBROUTINE Init_Lin
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
   REAL(ReKi), ALLOCATABLE                 :: dUdu(:,:), dUdy(:,:), G(:,:), tmp(:,:) ! variables for glue-code linearization
   INTEGER(IntKi), ALLOCATABLE             :: ipiv(:)
   integer(intki)                          :: nu, ny
   CHARACTER(1024)                         :: LinRootName
   CHARACTER(1024)                         :: OutFileName
   
   
   
   ErrStat = ErrID_None
   ErrMsg = ""

   
   LinRootName = TRIM(p_FAST%OutFileRoot)//'.'//trim(num2lstr(m_FAST%NextLinTimeIndx))
                     
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
      
      
      ! write the module matrices:
   if (p_FAST%LinOutMod) then
      
      OutFileName = trim(LinRootName)//'.'//TRIM(y_FAST%Module_Abrev(Module_ED))      
      call WrLinFile_txt_Head(t_global, p_FAST, y_FAST, y_FAST%Lin%Modules(Module_ED), OutFileName, Un, ErrStat2, ErrMsg2 )       
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat >=AbortErrLev) return
         
      if (p_FAST%LinOutJac) then
         ! Jacobians
         call WrMatrix( y_FAST%Lin%Modules(Module_ED)%D, Un, p_FAST%OutFmt, 'dYdu' )
         call WrMatrix( y_FAST%Lin%Modules(Module_ED)%B, Un, p_FAST%OutFmt, 'dXdu' )
         call WrMatrix( y_FAST%Lin%Modules(Module_ED)%C, Un, p_FAST%OutFmt, 'dYdx' )
         call WrMatrix( y_FAST%Lin%Modules(Module_ED)%A, Un, p_FAST%OutFmt, 'dXdx' )         
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
            call WrMatrix( y_FAST%Lin%Modules(Module_IfW)%D, Un, p_FAST%OutFmt, 'dYdu' )
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
            call WrMatrix( y_FAST%Lin%Modules(Module_SrvD)%D, Un, p_FAST%OutFmt, 'dYdu' )
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
            call WrMatrix( y_FAST%Lin%Modules(Module_AD)%D, Un, p_FAST%OutFmt, 'dYdu' )
            call WrMatrix( dZdu, Un, p_FAST%OutFmt, 'dZdu' )
            call WrMatrix( dYdz, Un, p_FAST%OutFmt, 'dYdz' )
            call WrMatrix( dZdz, Un, p_FAST%OutFmt, 'dZdz' )
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
      
      
      ! allocate the glue-code state matrices
   call Glue_FormDiag( p_FAST, y_FAST, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (ErrStat >=AbortErrLev) then
         call cleanup()
         return
      end if
   
      
      ! allocate and form matrix G:
   if (.not. allocated(G)) then
      call AllocAry(G, size(dUdu,1), size(dUdu,2), 'G', ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      call AllocAry(tmp, p_FAST%SizeLin(NumModules+1,LIN_INPUT_COL), p_FAST%SizeLin(NumModules+1,LIN_ContSTATE_COL), 'GyC', ErrStat2, ErrMsg2) ! product of G^-1 * dUdy * C
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
      if (ErrStat >=AbortErrLev) then
         call cleanup()
         return
      end if            
   end if
   
   
   call WrLinFile_txt_Head(t_global, p_FAST, y_FAST, y_FAST%Lin%Glue, LinRootName, Un, ErrStat2, ErrMsg2 )       
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (ErrStat >=AbortErrLev) then
         call cleanup()
         return
      end if
      
   
   if (p_FAST%LinOutJac) then
      ! Jacobians
      call WrMatrix( dUdu, Un, p_FAST%OutFmt, 'dUdu' )
      call WrMatrix( dUdy, Un, p_FAST%OutFmt, 'dUdy' )
   end if
   
   G = dUdu
   call LAPACK_GEMM( 'N', 'N', 1.0_ReKi, dUdy, y_FAST%Lin%Glue%D, 1.0_ReKi, G, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   !G = dUdu + matmul( dUdy, y_FAST%Lin%Glue%D )
   
   if (p_FAST%LinOutJac) call WrMatrix( G, Un, p_FAST%OutFmt, 'G' )
   
   
   ! now we need to form G^-1 * dUdy and G^-1 * dUdu 
   call allocAry( ipiv, size(G,1), 'ipiv', ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (ErrStat >= AbortErrLev) then
         call cleanup() 
         return
      end if
                  
   CALL LAPACK_getrf( M=size(G,1), N=size(G,2), A=G, IPIV=ipiv, ErrStat=ErrStat2, ErrMsg=ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (ErrStat >= AbortErrLev) then
         call cleanup() 
         return
      end if
               
   CALL LAPACK_getrs( trans='N', N=size(G,2), A=G, IPIV=ipiv, B=dUdy, ErrStat=ErrStat2, ErrMsg=ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)   
   CALL LAPACK_getrs( trans='N', N=size(G,2), A=G, IPIV=ipiv, B=dUdu, ErrStat=ErrStat2, ErrMsg=ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)   
   
   ! after the above solves, dUdy holds G^-1 * dUdy and dUdu holds G^-1 * dUdu
   ! now we can find the complete state matrices (note that I am not using matmul() because it easily gives stack overflow errors):
   
      ! temp variable:
   !tmp = G^-1 * dUdy * C = matmul( dUdy, y_FAST%Lin%Glue%C)
   call LAPACK_GEMM( 'N', 'N', 1.0_ReKi, dUdy, y_FAST%Lin%Glue%C, 0.0_ReKi, tmp, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
   !  A
   !y_FAST%Lin%Glue%A = y_FAST%Lin%Glue%A - matmul( y_FAST%Lin%Glue%B, tmp )  
   call LAPACK_GEMM( 'N', 'N', -1.0_ReKi, y_FAST%Lin%Glue%B, tmp, 1.0_ReKi, y_FAST%Lin%Glue%A, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

   !  C
   !y_FAST%Lin%Glue%C = y_FAST%Lin%Glue%C - matmul( y_FAST%Lin%Glue%D, tmp ) 
   call LAPACK_GEMM( 'N', 'N', -1.0_ReKi, y_FAST%Lin%Glue%D, tmp, 1.0_ReKi, y_FAST%Lin%Glue%C, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
   !  B
   deallocate(tmp)
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
      
    
   !  D
   deallocate(tmp)
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
   
      ! Write the results to the file:
   call WrLinFile_txt_End(Un, p_FAST, y_FAST%Lin%Glue )            

contains
   subroutine cleanup()
         ! these variables are currently for AeroDyn only:
      if (allocated(dYdz)) deallocate(dYdz)
      if (allocated(dZdz)) deallocate(dZdz)
      if (allocated(dZdu)) deallocate(dZdu)
      if (allocated(ipiv)) deallocate(ipiv)     
      
      if (allocated(dUdu)) deallocate(dUdu)
      if (allocated(dUdy)) deallocate(dUdy)
      if (allocated(G)) deallocate(G)
      if (allocated(tmp)) deallocate(tmp)
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
   INTEGER(IntKi)                          :: n(5)                ! sizes of arrays to print
   
   INTEGER(IntKi)                          :: ErrStat2            ! local error status
   CHARACTER(1024)                         :: ErrMsg2             ! local error message
   CHARACTER(*),             PARAMETER     :: RoutineName = 'WrLinFile_txt_Head'    
   CHARACTER(*),             PARAMETER     :: TypeNames(5) = (/ 'continuous states', &
                                                                'discrete states  ', &
                                                                'constraint states', &
                                                                'inputs           ', &
                                                                'outputs          '  /)
   CHARACTER(50)                           :: Fmt
   CHARACTER(32)                           :: Desc
   
                  
   ErrStat = ErrID_None
   ErrMsg = ""
         
   n = 0;
   if (allocated(LinData%names_x )) n(1) = size(LinData%names_x )
   if (allocated(LinData%names_xd)) n(2) = size(LinData%names_xd)
   if (allocated(LinData%names_z )) n(3) = size(LinData%names_z )
   if (allocated(LinData%names_u )) n(4) = size(LinData%names_u )
   if (allocated(LinData%names_y )) n(5) = size(LinData%names_y )
   
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
   
   WRITE (Un, '(A,'//trim(p_FAST%OutFmt_t)//',A,)') 'Linearized model at simulation time = ', t_global, ' s'
   
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
   if (allocated(LinData%names_x)) then
      WRITE(Un, '(A)') 'Order of continuous states:'      
      call WrLinFile_txt_Table(p_FAST, Un, "Row/Column", LinData%op_x, LinData%names_x  )      
      
      WRITE(Un, '(A)') 'Order of continuous state derivatives:'      
      call WrLinFile_txt_Table(p_FAST, Un, "Row/Column", LinData%op_dx, LinData%names_x, .true.  )      
   end if
   
   if (allocated(LinData%names_xd)) then
      WRITE(Un, '(A)') 'Order of discrete states:'      
      call WrLinFile_txt_Table(p_FAST, Un, "Row/Column", LinData%op_xd, LinData%names_xd  )      
   end if

   if (allocated(LinData%names_z)) then
      WRITE(Un, '(A)') 'Order of constraint states:'      
      call WrLinFile_txt_Table(p_FAST, Un, "Row/Column", LinData%op_z, LinData%names_z  )      
   end if
         
   if (allocated(LinData%names_u)) then
      WRITE(Un, '(A)') 'Order of inputs:'   
      call WrLinFile_txt_Table(p_FAST, Un, "Column  ", LinData%op_u, LinData%names_u  )      
   end if
   
   if (allocated(LinData%names_y)) then
      WRITE(Un, '(A)') 'Order of outputs:'      
      call WrLinFile_txt_Table(p_FAST, Un, "Row  ", LinData%op_y, LinData%names_y  )      
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
   if (allocated(LinData%A)) call WrMatrix( LinData%A, Un, p_FAST%OutFmt, 'A' )
   ! B matrix
   if (allocated(LinData%B)) call WrMatrix( LinData%B, Un, p_FAST%OutFmt, 'B' )
   ! C matrix
   if (allocated(LinData%C)) call WrMatrix( LinData%C, Un, p_FAST%OutFmt, 'C' )
   ! D matrix
   if (allocated(LinData%D)) call WrMatrix( LinData%D, Un, p_FAST%OutFmt, 'D' )
            

   close(un)
   
END SUBROUTINE WrLinFile_txt_End   
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE WrLinFile_txt_Table(p_FAST, Un, RowCol, op, names, deriv  )

   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< parameters
   INTEGER(IntKi),           INTENT(IN   ) :: Un                  !< unit number
   CHARACTER(*),             INTENT(IN   ) :: RowCol              !< Row/Column description
   REAL(ReKi),               INTENT(IN   ) :: op(:)               !< operating point values (possibly different size that Desc because of orientations)
   CHARACTER(LinChanLen),    INTENT(IN   ) :: names(:)            !< Descriptions of the channels (names and units)
   logical, optional,        intent(in   ) :: deriv               !< flag that tells us if we need to modify the channel names for derivatives (xdot)
   
      ! local variables
   INTEGER(IntKi)                          :: TS                  ! tab stop column
   INTEGER(IntKi)                          :: i                   ! loop counter
   INTEGER(IntKi)                          :: i_op                ! loop counter
   
   INTEGER(IntKi)                          :: ErrStat2            ! local error status
   CHARACTER(1024)                         :: ErrMsg2             ! local error message
   logical                                 :: UseDerivNames       !< flag that tells us if we need to modify the channel names for derivatives (xdot)
   CHARACTER(*),             PARAMETER     :: RoutineName = 'WrLinFile_txt_Table'    
   CHARACTER(100)                          :: Fmt
   CHARACTER(100)                          :: Fmt_Str
   CHARACTER(100)                          :: FmtOrient
   

   
   if (present(deriv) ) then
      UseDerivNames = deriv
   else
      UseDerivNames = .false.
   end if
   
   
   TS = 14 + 3*p_FAST%FmtWidth+7 
   
   Fmt       = '(3x,I8,3x,'//trim(p_FAST%OutFmt)//',T'//trim(num2lstr(TS))//',A)'
   FmtOrient = '(3x,I8,3x,'//trim(p_FAST%OutFmt)//',2(", ",'//trim(p_FAST%OutFmt)//'),2x,A)'
   Fmt_Str   = '(3x,A10,1x,A,T'//trim(num2lstr(TS))//',A)'
   
   WRITE(Un, Fmt_Str) RowCol,      'Operating Point', 'Description'
   WRITE(Un, Fmt_Str) '----------','---------------', '-----------'
   
   i_op = 1
   do i=1,size(names)
      if (index(names(i), ' orientation angle, node ') > 0 ) then  ! make sure this matches what is written in PackMotionMesh_Names()
         WRITE(Un, FmtOrient) i, op(i_op), op(i_op+1), op(i_op+2), trim(names(i))  !//' [OP is a row of the DCM]
         i_op = i_op + 3
      else
         if (UseDerivNames) then
            WRITE(Un, Fmt) i, op(i_op), 'First time derivative of '//trim(names(i))//'/s'
         else
            WRITE(Un, Fmt) i, op(i_op), trim(names(i))
         end if         
         i_op = i_op + 1
      end if      
   end do
         
   WRITE (Un,'()' )    !print a blank line
   
   
   
END SUBROUTINE WrLinFile_txt_Table   
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Glue_GetOP(p_FAST, y_FAST, ErrStat, ErrMsg)

   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< parameters
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST              !< Output variables for the glue code
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None

   
   INTEGER(IntKi)                          :: i, j                ! loop/temp variables
   INTEGER(IntKi)                          :: ThisModule          ! Module ID # 
   INTEGER(IntKi)                          :: i_u, i_y, i_x       ! loop/temp variables
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
         if (allocated(y_FAST%Lin%Modules(ThisModule)%op_u)) then
            i_u = i_u + size(y_FAST%Lin%Modules(ThisModule)%op_u)
         end if
         
         if (allocated(y_FAST%Lin%Modules(ThisModule)%op_y)) then
            i_y = i_y + size(y_FAST%Lin%Modules(ThisModule)%op_y)
         end if
      end do      
      
      call AllocAry( y_FAST%Lin%Glue%op_u, i_u, 'op_u', ErrStat2, ErrMsg2)
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
   INTEGER(IntKi)                          :: i, j, k             ! loop counters
   INTEGER(IntKi)                          :: i2, j2, k2          ! loop counters
   INTEGER(IntKi)                          :: r_start, r_end      ! row start/end of glue matrix
   INTEGER(IntKi)                          :: c_start             ! column start of glue matrix
   INTEGER(IntKi)                          :: Node                ! loop counters
   
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
   r_start = 1
   r_end   = p_FAST%SizeLin(Module_IfW,LIN_INPUT_COL)
   do i = r_start,r_end
      dUdu(i,i) = 1.0_ReKi
   end do

      !............
      ! \f$ \frac{\partial U_\Lambda^{SrvD}}{\partial u^{SrvD}} = I \f$ 
      !............
   r_start = r_end + 1
   r_end   = r_end + p_FAST%SizeLin(Module_SrvD,LIN_INPUT_COL)
   do i = r_start,r_end
      dUdu(i,i) = 1.0_ReKi
   end do
   
      !............
      ! \f$ \frac{\partial U_\Lambda^{ED}}{\partial u^{ED}} = I \f$ 
      !............
   r_start = r_end + 1
   r_end   = r_end + p_FAST%SizeLin(Module_ED,LIN_INPUT_COL)
   do i = r_start,r_end
      dUdu(i,i) = 1.0_ReKi
   end do
   
      !............
      ! \f$ \frac{\partial U_\Lambda^{IfW}}{\partial u^{AD}} \end{bmatrix} = \f$   
      !............
   IF (p_FAST%CompInflow == MODULE_IfW .and. p_FAST%CompAero == MODULE_AD) THEN  
      r_start = 1 ! start of IfW input equations
      c_start = p_FAST%SizeLin(Module_IfW,LIN_INPUT_COL) + p_FAST%SizeLin(Module_SrvD,LIN_INPUT_COL) + p_FAST%SizeLin(Module_ED,LIN_INPUT_COL) + 1 ! start of AD inputs
      call Linear_IfW_InputSolve_du_AD( p_FAST, AD%Input(1), dUdu, IfW_start=r_start, AD_start=c_start )
   end if ! we're using the InflowWind module
   
      !............
      ! \f$ \frac{\partial U_\Lambda^{ED}}{\partial u^{AD}} \end{bmatrix} = \f$   
      !............   
   IF (p_FAST%CompAero == MODULE_AD) THEN   ! we need to do this regardless of CompElast
      r_start = p_FAST%SizeLin(Module_IfW,LIN_INPUT_COL) + p_FAST%SizeLin(Module_SrvD,LIN_INPUT_COL) + 1 ! start of ED input equations
      c_start = p_FAST%SizeLin(Module_IfW,LIN_INPUT_COL) + p_FAST%SizeLin(Module_SrvD,LIN_INPUT_COL) + p_FAST%SizeLin(Module_ED,LIN_INPUT_COL) + 1 ! start of AD inputs
      call Linear_ED_InputSolve_du_AD( p_FAST, ED%Input(1), ED%Output(1), AD%y, AD%Input(1), MeshMapData, dUdu, ED_start=r_start, AD_start=c_start, ErrStat=ErrStat2, ErrMsg=ErrMsg2 )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   end if ! we're using the InflowWind module
   
      !............
      ! \f$ \frac{\partial U_\Lambda^{AD}}{\partial u^{AD}} \end{bmatrix} = \f$   
      !............
   IF (p_FAST%CompAero == MODULE_AD) THEN 
      r_start = p_FAST%SizeLin(Module_IfW,LIN_INPUT_COL) + p_FAST%SizeLin(Module_SrvD,LIN_INPUT_COL) + p_FAST%SizeLin(Module_ED,LIN_INPUT_COL) + 1 ! start of AD input equations
      r_end   = p_FAST%SizeLin(NumModules+1,LIN_INPUT_COL) ! end of AD inputs (& input equations)
      call Linear_AD_InputSolve_du_AD( p_FAST, AD%Input(1), ED%Output(1), MeshMapData, dUdu, AD_start=r_start, AD_end=r_end, ErrStat=ErrStat2, ErrMsg=ErrMsg2 )         
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
   
   !call Linear_IfW_InputSolve_dy( p_FAST, m_FAST, IfW%Input(1), IfW%p, AD14%Input(1), AD%Input(1), ED%Output(1), ErrStat2, ErrMsg2 )
      !............
      ! \f$ \frac{\partial U_\Lambda^{AD}}{\partial y^{IfW}} \end{bmatrix} = \f$   
      !............
   if (p_FAST%CompInflow == MODULE_IfW .and. p_FAST%CompAero == MODULE_AD) then   
      r_start = p_FAST%SizeLin(Module_IfW,LIN_INPUT_COL) + p_FAST%SizeLin(Module_SrvD,LIN_INPUT_COL) + p_FAST%SizeLin(Module_ED,LIN_INPUT_COL) + 1 !row=u_AD; start of AD input equations
      c_start = 1 !start of IfW outputs
      call Linear_AD_InputSolve_IfW_dy( p_FAST, AD%Input(1), dUdy, AD_Start=r_start, IfW_Start=c_start )      
   end if
   
   
      !............
      ! \f$ \frac{\partial U_\Lambda^{SrvD}}{\partial y^{ED}} \end{bmatrix} = \f$   
      !............
   if (p_FAST%CompServo == MODULE_SrvD) then   ! need to do this regardless of CompElast
      r_start = p_FAST%SizeLin(Module_IfW,LIN_INPUT_COL) + 1 !row = start of SrvD input equations

      c_start = p_FAST%SizeLin(Module_IfW,LIN_OUTPUT_COL) + p_FAST%SizeLin(Module_SrvD,LIN_OUTPUT_COL) + 1 !column=start of ED outputs
      c_start = c_start + p_FAST%SizeLin(Module_ED,LIN_OUTPUT_COL) - y_FAST%numOuts(Module_ED) - 3 ! start of ED where Yaw, YawRate, HSS_Spd occur (right before WriteOutputs)
            
      call Linear_SrvD_InputSolve_dy_ED( p_FAST, dUdy, SrvD_Start=r_start, ED_Start_Yaw=c_start )      
   end if
   
      !............
      ! \f$ \frac{\partial U_\Lambda^{ED}}{\partial y^{SrvD}} \end{bmatrix} = \f$   
      ! \f$ \frac{\partial U_\Lambda^{ED}}{\partial y^{ED}} \end{bmatrix} = \f$   
      ! \f$ \frac{\partial U_\Lambda^{ED}}{\partial y^{AD}} \end{bmatrix} = \f$   
      !............
   i2 = p_FAST%SizeLin(Module_IfW,LIN_OUTPUT_COL) + 1 ! y_SrvD start
   j2 = p_FAST%SizeLin(Module_IfW,LIN_OUTPUT_COL) + p_FAST%SizeLin(Module_SrvD,LIN_OUTPUT_COL) + 1 ! y_ED start
   k2 = p_FAST%SizeLin(Module_IfW,LIN_OUTPUT_COL) + p_FAST%SizeLin(Module_SrvD,LIN_OUTPUT_COL) + p_FAST%SizeLin(Module_ED,LIN_OUTPUT_COL) + 1 ! y_AD start
   
   r_start = p_FAST%SizeLin(Module_IfW,LIN_INPUT_COL) + p_FAST%SizeLin(Module_SrvD,LIN_INPUT_COL) + 1 ! start of u_ED
   r_end   = p_FAST%SizeLin(Module_IfW,LIN_INPUT_COL) + p_FAST%SizeLin(Module_SrvD,LIN_INPUT_COL) + p_FAST%SizeLin(Module_ED,LIN_INPUT_COL) ! end of u_ED
   
   call Linear_ED_InputSolve_dy( p_FAST, ED%Input(1), ED%Output(1), AD%y, AD%Input(1), MeshMapData, dUdy, &
              SrvD_start=i2, ED_start=r_start, ED_end=r_end, ED_Out_Start=j2, AD_start=k2, ErrStat=ErrStat2, ErrMsg=ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   
      !............
      ! \f$ \frac{\partial U_\Lambda^{AD}}{\partial y^{ED}} \end{bmatrix} = \f$   
      !............
   if (p_FAST%CompAero == MODULE_AD) then   ! need to do this regardless of CompElast
      r_start = p_FAST%SizeLin(Module_IfW,LIN_INPUT_COL) + p_FAST%SizeLin(Module_SrvD,LIN_INPUT_COL) + 1 ! start of ED input equations
      c_start = p_FAST%SizeLin(Module_IfW,LIN_OUTPUT_COL) + p_FAST%SizeLin(Module_SrvD,LIN_OUTPUT_COL) + 1 ! y_ED start
      
      call Linear_AD_InputSolve_NoIfW_dy( p_FAST, AD%Input(1), ED%Output(1), MeshMapData, dUdy, AD_Start=r_start, ED_Start=c_start, ErrStat=ErrStat2, ErrMsg=ErrMsg2 )      
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   end if
      
   
   
END SUBROUTINE Glue_Jacobians      
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine forms the dU^{IfW}/du^{AD} block of dUdu.
SUBROUTINE Linear_IfW_InputSolve_du_AD( p_FAST, u_AD, dUdu, ifW_Start, AD_Start )

   TYPE(FAST_ParameterType),       INTENT(IN   )   :: p_FAST      !< FAST parameter data 
   TYPE(AD_InputType),             INTENT(IN)      :: u_AD        !< The input meshes (already calculated) from AeroDyn
   real(reki),                     INTENT(INOUT)   :: dUdu(:,:)   !< Jacobian matrix of which we are computing the dU^(IfW)/du^(AD) block
   integer(intKi),                 INTENT(IN)      :: ifW_Start   !< starting index of dUdu (row) where ifW input equations are located 
   integer(intKi),                 INTENT(IN)      :: AD_Start    !< starting index of dUdu (column) where AD inputs are located 
   
   
   INTEGER(IntKi)                          :: i, j, k             ! loop counters
   INTEGER(IntKi)                          :: i2, j2              ! loop counters
   INTEGER(IntKi)                          :: AD_Start_Bl         ! starting index of dUdu (column) where AD blade motion inputs are located
   INTEGER(IntKi)                          :: Node                ! InflowWind node number
   
            
         ! compare with IfW_InputSolve():
   
      Node = 0 !InflowWind node
      if (p_FAST%CompServo == MODULE_SrvD) Node = Node + 1
            
      IF (p_FAST%CompAero == MODULE_AD) THEN 
         
            ! blades:
         AD_Start_Bl = AD_Start + u_AD%TowerMotion%NNodes * 9  & ! 3 fields (MASKID_TRANSLATIONDISP,MASKID_Orientation,MASKID_TRANSLATIONVel) with 3 components
                                + u_AD%HubMotion%NNodes   * 9    ! 3 fields (MASKID_TRANSLATIONDISP,MASKID_Orientation,MASKID_RotationVel) with 3 components
   
         do k = 1,size(u_AD%BladeRootMotion)         
            AD_Start_Bl = AD_Start_Bl + u_AD%BladeRootMotion(k)%NNodes * 3 ! 1 field (MASKID_Orientation) with 3 components
         end do
         ! next is u_AD%BladeMotion(k); note that it has 3 fields and we only need 1
                  
         DO K = 1,SIZE(u_AD%BladeMotion)
            DO J = 1,u_AD%BladeMotion(k)%Nnodes
               Node = Node + 1 ! InflowWind node
               do i=1,3 !XYZ components of this node
                  i2 = ifW_Start + (Node-1)*3 + i - 1
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
               i2 = (Node-1)*3 + i ! = i_start + (Node-1)*3 + i - 1
               j2 = AD_Start + (j-1)*3 + i - 1
               dUdu( i2, j2 ) = -1.0_ReKi
            end do            
         END DO              
         
      END IF     
   
END SUBROUTINE Linear_IfW_InputSolve_du_AD
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine forms the dU^{ED}/du^{AD} block of dUdu.
SUBROUTINE Linear_ED_InputSolve_du_AD( p_FAST, u_ED, y_ED, y_AD, u_AD, MeshMapData, dUdu, ED_start, AD_start, ErrStat, ErrMsg )

   TYPE(FAST_ParameterType),       INTENT(IN   )  :: p_FAST         !< Glue-code simulation parameters
   TYPE(ED_InputType),             INTENT(INOUT)  :: u_ED           !< ED Inputs at t
   TYPE(ED_OutputType),            INTENT(IN   )  :: y_ED           !< ElastoDyn outputs (need translation displacement on meshes for loads mapping)
   TYPE(AD_OutputType),            INTENT(IN   )  :: y_AD           !< AeroDyn outputs
   TYPE(AD_InputType),             INTENT(INOUT)  :: u_AD           !< AD inputs (for AD-ED load linerization)
   
   TYPE(FAST_ModuleMapType),       INTENT(INOUT)  :: MeshMapData    !< Data for mapping between modules
   real(reki),                     INTENT(INOUT)  :: dUdu(:,:)      !< Jacobian matrix of which we are computing the dU^(ED)/du^(AD) block
   integer(intKi),                 INTENT(IN)     :: ED_Start       !< starting index of dUdu (row) where ED input equations are located 
   integer(intKi),                 INTENT(IN)     :: AD_Start       !< starting index of dUdu (column) where AD inputs are located 
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
         AD_Start_Bl = AD_Start + u_AD%TowerMotion%NNodes * 9  & ! 3 fields (MASKID_TRANSLATIONDISP,MASKID_Orientation,MASKID_TRANSLATIONVel) with 3 components
                                + u_AD%HubMotion%NNodes   * 9    ! 3 fields (MASKID_TRANSLATIONDISP,MASKID_Orientation,MASKID_RotationVel) with 3 components
   
         do k = 1,size(u_AD%BladeRootMotion)         
            AD_Start_Bl = AD_Start_Bl + u_AD%BladeRootMotion(k)%NNodes * 3 ! 1 field (MASKID_Orientation) with 3 components
         end do
         ! next is u_AD%BladeMotion(k); note that it has 3 fields and we only need 1
      
         ED_Start_mt = ED_Start
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
         ED_Start_mt = ED_Start           
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
                  dUdu( ED_Start_mt + j - 1, AD_Start + i - 1 ) = - MeshMapData%AD_L_2_ED_P_T%dM%m_us(j,i)
               end do
            end do
               
         end if                                    
      END IF
            
   !END IF
               
END SUBROUTINE Linear_ED_InputSolve_du_AD
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine forms the dU^{AD}/du^{AD} block of dUdu.
SUBROUTINE Linear_AD_InputSolve_du_AD( p_FAST, u_AD, y_ED, MeshMapData, dUdu, AD_start, AD_end, ErrStat, ErrMsg )

      ! Passed variables
   TYPE(FAST_ParameterType),    INTENT(IN   )   :: p_FAST      !< FAST parameter data    
   TYPE(AD_InputType),          INTENT(INOUT)   :: u_AD        !< The inputs to AeroDyn14
   TYPE(ED_OutputType),         INTENT(IN)      :: y_ED        !< The outputs from the structural dynamics module
   TYPE(FAST_ModuleMapType),    INTENT(INOUT)   :: MeshMapData !< Data for mapping between modules
   real(reki),                  INTENT(INOUT)   :: dUdu(:,:)   !< Jacobian matrix of which we are computing the dU^(ED)/du^(AD) block
   integer(intKi),              INTENT(IN)      :: AD_Start    !< starting index of dUdu (column) where AD inputs are located 
   integer(intKi),              INTENT(IN)      :: AD_End      !< last index of dUdu (column) where AD inputs are located 
   
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
   do j=AD_Start, AD_End
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

      AD_Start_td = AD_Start   
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
      
      AD_Start_td = AD_Start + u_AD%TowerMotion%NNodes * 9  & ! 3 fields (MASKID_TRANSLATIONDISP,MASKID_Orientation,MASKID_TRANSLATIONVel) with 3 components
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
SUBROUTINE Linear_SrvD_InputSolve_dy_ED( p_FAST, dUdy, SrvD_Start, ED_Start_Yaw )
!..................................................................................................................................

   TYPE(FAST_ParameterType),       INTENT(IN)     :: p_FAST         !< Glue-code simulation parameters
   real(reki),                     INTENT(INOUT)  :: dUdy(:,:)      !< Jacobian matrix of which we are computing the dU^{SrvD}/dy^{ED} block
   integer(intKi),                 INTENT(IN)     :: SrvD_Start     !< starting index of dUdy (row) where SrvD input equations are located 
   integer(intKi),                 INTENT(IN)     :: ED_Start_Yaw   !< starting index of dUdy (column) where ED Yaw/YawRate/HSS_Spd outputs are located (just before WriteOutput)

   
   
   INTEGER(IntKi)                                   :: i            ! loop counter
   
   INTEGER(IntKi)                                   :: ErrStat2                 ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                             :: ErrMsg2                  ! temporary Error message if ErrStat /= ErrID_None
   CHARACTER(*), PARAMETER                          :: RoutineName = 'Linear_SrvD_InputSolve_dy_ED' 
   
   
   do i=1,3
      dUdy(SrvD_Start + i - 1, ED_Start_Yaw + i - 1) = -1.0_ReKi
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
SUBROUTINE Linear_ED_InputSolve_dy( p_FAST, u_ED, y_ED, y_AD, u_AD, MeshMapData, dUdy, SrvD_start, ED_start, ED_end, ED_Out_Start, AD_start, ErrStat, ErrMsg )

   TYPE(FAST_ParameterType),       INTENT(IN   )  :: p_FAST           !< Glue-code simulation parameters
   TYPE(ED_InputType),             INTENT(INOUT)  :: u_ED             !< ED Inputs at t
   TYPE(ED_OutputType),            INTENT(IN   )  :: y_ED             !< ElastoDyn outputs (need translation displacement on meshes for loads mapping)
   TYPE(AD_OutputType),            INTENT(IN   )  :: y_AD             !< AeroDyn outputs
   TYPE(AD_InputType),             INTENT(INOUT)  :: u_AD             !< AD inputs (for AD-ED load linerization)
                                                                      
   TYPE(FAST_ModuleMapType),       INTENT(INOUT)  :: MeshMapData      !< Data for mapping between modules
   real(reki),                     INTENT(INOUT)  :: dUdy(:,:)        !< Jacobian matrix of which we are computing the dU^(ED)/du^(AD) block
   integer(intKi),                 INTENT(IN)     :: ED_Start         !< starting index of dUdy (row) where ED input equations are located 
   integer(intKi),                 INTENT(IN)     :: ED_End           !< last index of dUdy (row) where ED input equations are located 
   integer(intKi),                 INTENT(IN)     :: ED_Out_Start     !< starting index of dUdy (column) where ED outputs are located 
   integer(intKi),                 INTENT(IN)     :: SrvD_Start       !< starting index of dUdy (column) where SrvD outputs are located 
   integer(intKi),                 INTENT(IN)     :: AD_Start         !< starting index of dUdy (column) where AD outputs are located 
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
   
      ! ED inputs from ServoDyn
   IF ( p_FAST%CompServo == Module_SrvD ) THEN

         ! BlPitchCom, YawMom, GenTrq
      ED_Start_tmp = ED_End - size(u_ED%BlPitchCom) - 2
      do i=1,size(u_ED%BlPitchCom)+2
         dUdy(ED_Start_tmp + i - 1, SrvD_Start + i - 1) = -1.0_ReKi
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
         
         AD_Start_tmp = AD_Start + u_AD%TowerMotion%NNodes * 6    ! 2 fields (force, moment) with 3 components
         ED_Start_tmp = ED_Start     ! blades in u_ED
         ED_Out_Start_tmp  = ED_Out_Start ! blades in y_ED
         
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
            
            
         AD_Start_tmp = AD_Start ! tower is first in y_AD
         
         ! find tower in u_ED:
         ED_Start_tmp = ED_Start 
         if (allocated(u_ED%BladePtLoads)) then
            do K = 1,SIZE(u_ED%BladePtLoads,1) ! Loop through all blades (p_ED%NumBl)
               ED_Start_tmp = ED_Start_tmp + u_ED%BladePtLoads(k)%NNodes*6 ! 2 fields (force/moment) with 3 components on each blade
            end do
         end if
         ED_Start_tmp = ED_Start_tmp + u_ED%PlatformPtMesh%NNodes*6 ! 2 fields (force/moment) with 3 components
         
         ! find tower in y_ED:
         ED_Out_Start_tmp  = ED_Out_Start 
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
SUBROUTINE Linear_AD_InputSolve_IfW_dy( p_FAST, u_AD, dUdy, AD_Start, IfW_Start )

      ! Passed variables
   TYPE(FAST_ParameterType),       INTENT(IN   )  :: p_FAST         !< FAST parameter data    
   TYPE(AD_InputType),             INTENT(INOUT)  :: u_AD           !< The inputs to AeroDyn
   real(reki),                     INTENT(INOUT)  :: dUdy(:,:)      !< Jacobian matrix of which we are computing the dU^{AD}/dy^{IfW} block
   integer(intKi),                 INTENT(IN)     :: AD_Start       !< starting index of dUdy (row) where AD input equations are located 
   integer(intKi),                 INTENT(IN)     :: IfW_Start      !< starting index of dUdy (column) where IfW outputs are located

      ! Local variables:

   INTEGER(IntKi)                               :: I           ! Loops through components
   INTEGER(IntKi)                               :: J           ! Loops through nodes / elements
   INTEGER(IntKi)                               :: K           ! Loops through blades
   INTEGER(IntKi)                               :: NodeNum     ! Node number for blade/node on mesh
   INTEGER(IntKi)                               :: NumBl
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
      
      
      AD_Start_tmp = AD_Start + u_AD%TowerMotion%NNodes * 9  & ! 3 fields (MASKID_TRANSLATIONDISP,MASKID_Orientation,MASKID_TRANSLATIONVel) with 3 components
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
               dUdy( AD_Start_tmp + i - 1, IfW_Start + (node-1)*3 + i - 1 ) = -1.0_ReKi                
            end do
            node = node + 1
            AD_Start_tmp = AD_Start_tmp + 3
         end do         
      end do
                  
      if ( allocated(u_AD%InflowOnTower) ) then         
         do j=1,size(u_AD%InflowOnTower,2) !nodes
            do i=1,3 !velocity component
               dUdy( AD_Start_tmp + i - 1, IfW_Start + (node-1)*3 + i - 1 ) = -1.0_ReKi                
            end do
            node = node + 1
            AD_Start_tmp = AD_Start_tmp + 3
         end do      
      end if
                     
   !END IF
   
   
END SUBROUTINE Linear_AD_InputSolve_IfW_dy
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine forms the dU^{AD}/dy^{ED} block of dUdy.
SUBROUTINE Linear_AD_InputSolve_NoIfW_dy( p_FAST, u_AD, y_ED, MeshMapData, dUdy, AD_Start, ED_Start, ErrStat, ErrMsg )

      ! Passed variables
   TYPE(FAST_ParameterType),    INTENT(IN   )   :: p_FAST      !< FAST parameter data    
   TYPE(AD_InputType),          INTENT(INOUT)   :: u_AD        !< The inputs to AeroDyn14
   TYPE(ED_OutputType),         INTENT(IN)      :: y_ED        !< The outputs from the structural dynamics module
   TYPE(FAST_ModuleMapType),    INTENT(INOUT)   :: MeshMapData !< Data for mapping between modules
   real(reki),                  INTENT(INOUT)   :: dUdy(:,:)   !< Jacobian matrix of which we are computing the dU^{AD}/dy^{ED} block
   integer(intKi),              INTENT(IN)      :: AD_Start    !< starting index of dUdy (row) where AD input equations are located 
   integer(intKi),              INTENT(IN)      :: ED_Start    !< starting index of dUdy (column) where ED outputs are located
   
   INTEGER(IntKi)                               :: ErrStat     !< Error status of the operation
   CHARACTER(*)                                 :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      ! Local variables:

   INTEGER(IntKi)                               :: i,j         ! Loops through rows/columns of mesh-mapping linearization matrices
   INTEGER(IntKi)                               :: K           ! Loops through blades
   INTEGER(IntKi)                               :: NodeNum     ! Node number for blade/node on mesh
   INTEGER(IntKi)                               :: NumBl
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
   AD_Start_tmp = AD_Start
   
      ! tower comes after ED blades and platform:
   ED_Start_tmp = ED_Start
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
            
      ED_Start_tmp = ED_Start
      
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


!----------------------------------------------------------------------------------------------------------------------------------
END MODULE FAST_Linear
