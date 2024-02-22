!**********************************************************************************************************************************
! FAST_ModLin.f90 performs linearization using the ModVars module.
!..................................................................................................................................
! LICENSING
! Copyright (C) 2024  National Renewable Energy Laboratory
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
module FAST_ModLin

use NWTC_Library
use NWTC_LAPACK

use FAST_Types
use FAST_Funcs
use FAST_Mapping

implicit none

private
public :: ModLin_Init, ModLin_Linearize_OP

contains

subroutine ModLin_Init(ModGlue, Mods, p, m, p_FAST, m_FAST, Turbine, ErrStat, ErrMsg)

   type(ModDataType), intent(inout)                :: ModGlue  !< Module data for glue code
   type(ModDataType), allocatable, intent(inout)   :: Mods(:)  !< Data for all modules
   type(ML_ParameterType), intent(inout)           :: p        !< ModLin parameters
   type(ML_MiscVarType), intent(inout)             :: m        !< ModLin miscvars
   type(FAST_ParameterType), intent(inout)         :: p_FAST
   type(FAST_MiscVarType), intent(inout)           :: m_FAST
   type(FAST_TurbineType), intent(inout)           :: Turbine
   integer(IntKi), intent(out)                     :: ErrStat
   character(*), intent(out)                       :: ErrMsg

   character(*), parameter                         :: RoutineName = 'ModLin_Init'
   integer(IntKi)                                  :: ErrStat2
   character(ErrMsgLen)                            :: ErrMsg2
   integer(IntKi), allocatable                     :: modIDs(:), modIdx(:)
   integer(IntKi)                                  :: i, j, k
   integer(IntKi)                                  :: FlagFilters
   character(LinChanLen), allocatable              :: xLinNames(:), uLinNames(:), yLinNames(:)
   character(20)                                   :: NamePrefix

   ! Initialize error return
   ErrStat = ErrID_None
   ErrMsg = ""

   !----------------------------------------------------------------------------
   ! FAST Lin Settings
   !----------------------------------------------------------------------------

   m_FAST%Lin%NextLinTimeIndx = 1
   m_FAST%Lin%CopyOP_CtrlCode = MESH_NEWCOPY
   m_FAST%Lin%n_rot = 0
   m_FAST%Lin%IsConverged = .false.
   m_FAST%Lin%FoundSteady = .false.
   m_FAST%Lin%ForceLin = .false.
   m_FAST%Lin%AzimIndx = 1

   p_FAST%AzimDelta = TwoPi/p_FAST%NLinTimes

   !----------------------------------------------------------------------------
   ! Module order and indexing
   !----------------------------------------------------------------------------

   ! If no modules were added, return error
   if (.not. allocated(Mods)) then
      call SetErrStat(ErrID_Fatal, "No modules were used", ErrStat, ErrMsg, RoutineName)
      return
   end if

   ! Create array of indices for Mods array
   modIdx = [(i, i=1, size(Mods))]

   ! Get array of module IDs
   modIDs = [(Mods(i)%ID, i=1, size(Mods))]

   ! Establish module index order for linearization
   p%iMod = [pack(modIdx, ModIDs == Module_IfW), &   ! InflowWind
             pack(modIdx, ModIDs == Module_SrvD), &  ! ServoDyn
             pack(modIdx, ModIDs == Module_ED), &    ! ElastoDyn
             pack(modIdx, ModIDs == Module_BD), &    ! BeamDyn
             pack(modIdx, ModIDs == Module_AD), &    ! AeroDyn
             pack(modIdx, ModIDs == Module_SeaSt), & ! SeaState
             pack(modIdx, ModIDs == Module_HD), &    ! HydroDyn
             pack(modIdx, ModIDs == Module_SD), &    ! SubDyn
             pack(modIdx, ModIDs == Module_MAP), &   ! MAP++
             pack(modIdx, ModIDs == Module_MD)]      ! MoorDyn

   ! Loop through modules, if module is not in index, return with error
   do i = 1, size(Mods)
      if (.not. any(i == p%iMod)) then
         call SetErrStat(ErrID_Fatal, "Module "//trim(Mods(i)%Abbr)//" not supported in linearization", &
                         ErrStat, ErrMsg, RoutineName)
         return
      end if
   end do

   !----------------------------------------------------------------------------
   ! Glue Module Variables
   !----------------------------------------------------------------------------

   ! Allocate variable structure for glue
   allocate (ModGlue%Vars)

   ! Initialize number of values in each variable group
   ModGlue%Vars%Nx = 0
   ModGlue%Vars%Nxd = 0
   ModGlue%Vars%Nz = 0
   ModGlue%Vars%Nu = 0
   ModGlue%Vars%Ny = 0

   ! Allocate arrays for glue variables
   allocate (ModGlue%Vars%x(0), ModGlue%Vars%xd(0), ModGlue%Vars%z(0), ModGlue%Vars%u(0), ModGlue%Vars%y(0))

   ! Loop through each module by index
   do i = 1, size(p%iMod)
      associate (ModData => Mods(p%iMod(i)))

         ! Create variable name prefix for linearization names. Add instance
         ! number to module abbreviation if more than 1 instance or the module is BeamDyn
         NamePrefix = ModData%Abbr
         if ((ModData%ID == Module_BD) .or. (count(modIDs == ModData%ID) > 1)) then
            NamePrefix = trim(NamePrefix)//"_"//Num2LStr(ModData%Ins)
         end if

         !----------------------------------------------------------------------
         ! Module continuous state variables
         !----------------------------------------------------------------------

         ! Set linearize flag on all variables
         do j = 1, size(ModData%Vars%x)
            call MV_SetFlags(ModData%Vars%x(j), VF_Linearize)
         end do

         ! Set module data start index in global arrays, increment data size
         ModData%ixg = ModGlue%Vars%Nx + 1
         ModGlue%Vars%Nx = ModGlue%Vars%Nx + ModData%Vars%Nx

         ! Save start index of module variables and append to glue code variables
         k = size(ModGlue%Vars%x) + 1
         ModGlue%Vars%x = [ModGlue%Vars%x, ModData%Vars%x]

         ! Loop through added variables and add name prefix to linearization names
         call AddLinNamePrefix(ModGlue%Vars%x(k:), NamePrefix)

         !----------------------------------------------------------------------
         ! Module discrete state variables
         !----------------------------------------------------------------------

         ! Set module data start index in global arrays, increment data size
         ModData%ixdg = ModGlue%Vars%Nxd + 1
         ModGlue%Vars%Nxd = ModGlue%Vars%Nxd + ModData%Vars%Nxd

         ! Save start index of module variables and append to glue code variables
         k = size(ModGlue%Vars%xd) + 1
         ModGlue%Vars%xd = [ModGlue%Vars%xd, ModData%Vars%xd]

         ! Loop through added variables and add name prefix to linearization names
         call AddLinNamePrefix(ModGlue%Vars%xd(k:), NamePrefix)

         !----------------------------------------------------------------------
         ! Module constraint state variables
         !----------------------------------------------------------------------

         ! Set module data start index in global arrays, increment data size
         ModData%izg = ModGlue%Vars%Nz + 1
         ModGlue%Vars%Nz = ModGlue%Vars%Nz + ModData%Vars%Nz

         ! Save start index of module variables and append to glue code variables
         k = size(ModGlue%Vars%z) + 1
         ModGlue%Vars%z = [ModGlue%Vars%z, ModData%Vars%z]

         ! Loop through added variables and add name prefix to linearization names
         call AddLinNamePrefix(ModGlue%Vars%z(k:), NamePrefix)

         !----------------------------------------------------------------------
         ! Module input variables
         !----------------------------------------------------------------------

         ! Add or remove linearize flag based on requested output
         select case (p_FAST%LinInputs)
         case (LIN_NONE)
            do j = 1, size(ModData%Vars%u)
               call MV_UnsetFlags(ModData%Vars%u(j), VF_Linearize)
            end do
         case (LIN_STANDARD)
            ! For standard inputs, use VF_Linearize flag set in the module
         case (LIN_ALL)
            do j = 1, size(ModData%Vars%u)
               call MV_SetFlags(ModData%Vars%u(j), VF_Linearize)
            end do
         end select

         ! Set module data start index in global arrays, increment data size
         ModData%iug = ModGlue%Vars%Nu + 1
         ModGlue%Vars%Nu = ModGlue%Vars%Nu + ModData%Vars%Nu

         ! Save start index of module variables and append to glue code variables
         k = size(ModGlue%Vars%u) + 1
         ModGlue%Vars%u = [ModGlue%Vars%u, ModData%Vars%u]

         ! Loop through added variables and add name prefix to linearization names
         call AddLinNamePrefix(ModGlue%Vars%u(k:), NamePrefix)

         !----------------------------------------------------------------------
         ! Module output variables
         !----------------------------------------------------------------------

         ! Add or remove linearize flag based on requested output
         select case (p_FAST%LinOutputs)
         case (LIN_NONE)
            do j = 1, size(ModData%Vars%y)
               call MV_UnsetFlags(ModData%Vars%y(j), VF_Linearize)
            end do
         case (LIN_STANDARD)  ! Set linearize flag for write output variables
            do j = 1, size(ModData%Vars%y)
               if (MV_HasFlags(ModData%Vars%y(j), VF_WriteOut)) then
                  call MV_SetFlags(ModData%Vars%y(j), VF_Linearize)
               else
                  call MV_UnsetFlags(ModData%Vars%y(j), VF_Linearize)
               end if
            end do
         case (LIN_ALL)
            do j = 1, size(ModData%Vars%y)
               call MV_SetFlags(ModData%Vars%y(j), VF_Linearize)
            end do
         end select

         ! Set module data start index in global arrays, increment data size
         ModData%iyg = ModGlue%Vars%Ny + 1
         ModGlue%Vars%Ny = ModGlue%Vars%Ny + ModData%Vars%Ny

         ! Save start index of module variables and append to glue code variables
         k = size(ModGlue%Vars%y) + 1
         ModGlue%Vars%y = [ModGlue%Vars%y, ModData%Vars%y]

         ! Loop through added variables and add name prefix to linearization names
         call AddLinNamePrefix(ModGlue%Vars%y(k:), NamePrefix)

         ! Initialize module linearization variable indexing
         call MV_InitVarIdx(ModData%Vars, ModData%Vars%IdxLin, VF_Linearize, ErrStat2, ErrMsg2); if (Failed()) return

      end associate
   end do

   ! Calculate number of values in each group and set data location index
   call CalcVarDataLoc(ModGlue%Vars%x, ModGlue%Vars%Nx)
   call CalcVarDataLoc(ModGlue%Vars%xd, ModGlue%Vars%Nxd)
   call CalcVarDataLoc(ModGlue%Vars%z, ModGlue%Vars%Nz)
   call CalcVarDataLoc(ModGlue%Vars%u, ModGlue%Vars%Nu)
   call CalcVarDataLoc(ModGlue%Vars%y, ModGlue%Vars%Ny)

   ! Initialize linearization index filtering
   call MV_InitVarIdx(ModGlue%Vars, ModGlue%Vars%IdxLin, VF_Linearize, ErrStat2, ErrMsg2); if (Failed()) return

   !----------------------------------------------------------------------------
   ! Allocate linearization arrays and matrices
   !----------------------------------------------------------------------------

   ! Allocate linearization arrays
   call AllocAry(ModGlue%Lin%x, ModGlue%Vars%Nx, "x", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(ModGlue%Lin%dx, ModGlue%Vars%Nx, "dx", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(ModGlue%Lin%xd, ModGlue%Vars%Nxd, "xd", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(ModGlue%Lin%z, ModGlue%Vars%Nz, "z", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(ModGlue%Lin%u, ModGlue%Vars%Nu, "u", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(ModGlue%Lin%y, ModGlue%Vars%Ny, "y", ErrStat2, ErrMsg2); if (Failed()) return

   ! Allocate full Jacobian matrices
   call AllocAry(ModGlue%Lin%dYdu, ModGlue%Vars%Ny, ModGlue%Vars%Nu, "dYdu", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(ModGlue%Lin%dXdu, ModGlue%Vars%Nx, ModGlue%Vars%Nu, "dXdu", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(ModGlue%Lin%dYdx, ModGlue%Vars%Ny, ModGlue%Vars%Nx, "dYdx", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(ModGlue%Lin%dXdx, ModGlue%Vars%Nx, ModGlue%Vars%Nx, "dXdx", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(ModGlue%Lin%dUdu, ModGlue%Vars%Nu, ModGlue%Vars%Nu, "dUdu", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(ModGlue%Lin%dUdy, ModGlue%Vars%Nu, ModGlue%Vars%Ny, "dUdy", ErrStat2, ErrMsg2); if (Failed()) return

   !----------------------------------------------------------------------------
   ! Mesh Mapping
   !----------------------------------------------------------------------------

   call FAST_InitMappings(Mods, m%Mappings, Turbine, ErrStat2, ErrMsg2); if (Failed()) return

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function Failed
   subroutine CalcVarDataLoc(VarAry, DataSize)
      type(ModVarType), intent(inout)  :: VarAry(:)
      integer(IntKi), intent(out)      :: DataSize
      DataSize = 0
      do i = 1, size(VarAry)
         VarAry(i)%iLoc = [DataSize + 1, DataSize + VarAry(i)%Num]
         DataSize = DataSize + VarAry(i)%Num
      end do
   end subroutine
end subroutine

subroutine AddLinNamePrefix(VarAry, Prefix)
   type(ModVarType), intent(inout)  :: VarAry(:)
   character(*), intent(in)         :: Prefix
   integer(IntKi)                   :: i, j
   do i = 1, size(VarAry)
      if (allocated(VarAry(i)%LinNames)) then
         do j = 1, size(VarAry(i)%LinNames)
            VarAry(i)%LinNames(j) = trim(Prefix)//" "//VarAry(i)%LinNames(j)
         end do
      end if
   end do
end subroutine

subroutine ModLin_Linearize_OP(Turbine, ModGlue, Mods, p, m, p_FAST, m_FAST, y_FAST, t_global, ErrStat, ErrMsg)

   type(ModDataType), intent(inout)          :: ModGlue  !< Module data for glue code
   type(ModDataType), intent(inout)          :: Mods(:)  !< Data for all modules
   type(ML_ParameterType), intent(inout)     :: p        !< ModLin parameters
   type(ML_MiscVarType), intent(inout)       :: m        !< ModLin MiscVars
   type(FAST_ParameterType), intent(in)      :: p_FAST
   type(FAST_MiscVarType), intent(inout)     :: m_FAST
   type(FAST_OutputFileType), intent(inout)  :: y_FAST
   type(FAST_TurbineType), intent(inout)     :: Turbine  !< Turbine type
   real(DbKi), intent(IN)                    :: t_global !< current (global) simulation time
   integer(IntKi), intent(out)               :: ErrStat
   character(*), intent(out)                 :: ErrMsg

   character(*), parameter                   :: RoutineName = 'ModLin_Init'
   integer(IntKi)                            :: ErrStat2
   character(ErrMsgLen)                      :: ErrMsg2
   integer(IntKi)                            :: i, j, k
   integer(IntKi)                            :: ix, ixd, iz, iu, iy
   integer(IntKi)                            :: Un
   character(200)                            :: SimStr
   character(MaxWrScrLen)                    :: BlankLine
   character(1024)                           :: LinRootName
   character(1024)                           :: OutFileName
   character(*), parameter                   :: Fmt = 'F10.2'

   ! Initialize error return
   ErrStat = ErrID_None
   ErrMsg = ""

   ! Write message to screen
   BlankLine = ""
   call WrOver(BlankLine)  ! BlankLine contains MaxWrScrLen spaces
   SimStr = '(RotSpeed='//trim(Num2LStr(Turbine%ED%y%RotSpeed*RPS2RPM, Fmt))//' rpm, BldPitch1='//trim(Num2LStr(Turbine%ED%y%BlPitch(1)*R2D, Fmt))//' deg)'
   call WrOver(' Performing linearization '//trim(Num2LStr(Turbine%m_FAST%Lin%NextLinTimeIndx))//' at simulation time '//TRIM(Num2LStr(t_global))//' s. '//trim(SimStr))
   call WrScr('')

   ! Get parameters
   ! NumBl = size(T%ED%Input(1)%BlPitchCom)
   y_FAST%Lin%RotSpeed = Turbine%ED%y%RotSpeed
   y_FAST%Lin%Azimuth = Turbine%ED%y%LSSTipPxa

   ! Assemble linearization root file name
   LinRootName = trim(p_FAST%OutFileRoot)//'.'//trim(Num2LStr(m_FAST%Lin%NextLinTimeIndx))

   ! Get unit number for writing files
   call GetNewUnit(Un, ErrStat2, ErrMsg2); if (Failed()) return

   ! Initialize the index numbers
   ix = 1
   ixd = 1
   iz = 1
   iu = 1
   iy = 1

   ! Initialize data in Jacobian matrices to zero
   ModGlue%Lin%dYdu = 0.0_R8Ki
   ModGlue%Lin%dXdu = 0.0_R8Ki
   ModGlue%Lin%dYdx = 0.0_R8Ki
   ModGlue%Lin%dXdx = 0.0_R8Ki

   ! Loop through modules by index
   do i = 1, size(p%iMod)
      associate (ModData => Mods(p%iMod(i)))

         ! Operating point values
         call FAST_GetOP(ModData, t_global, STATE_CURR, Turbine, ErrStat2, ErrMsg2, &
                         u_op=ModData%Lin%u, y_op=ModData%Lin%y, &
                         x_op=ModData%Lin%x, dx_op=ModData%Lin%dx)
         if (Failed()) return

         ! Derivatives wrt input
         call FAST_JacobianPInput(ModData, t_global, STATE_CURR, Turbine, ErrStat2, ErrMsg2, &
                                  dYdu=ModData%Lin%dYdu, dXdu=ModData%Lin%dXdu)
         if (Failed()) return


         ! Derivatives wrt continuous state
         call FAST_JacobianPContState(ModData, t_global, STATE_CURR, Turbine, ErrStat2, ErrMsg2, &
                                      dYdx=ModData%Lin%dYdx, dXdx=ModData%Lin%dXdx, &
                                      StateRotation=ModData%Lin%StateRotation)
         if (Failed()) return

         ! Copy module linearization arrays into glue linearization arrays
         if ((size(ModGlue%Lin%x) > 0) .and. allocated(ModData%Lin%x)) ModGlue%Lin%x(ix:ix + ModData%Vars%Nx - 1) = ModData%Lin%x
         if ((size(ModGlue%Lin%dx) > 0) .and. allocated(ModData%Lin%dx)) ModGlue%Lin%dx(ix:ix + ModData%Vars%Nx - 1) = ModData%Lin%dx
         if ((size(ModGlue%Lin%xd) > 0) .and. allocated(ModData%Lin%xd)) ModGlue%Lin%xd(ixd:ixd + ModData%Vars%Nxd - 1) = ModData%Lin%xd
         if ((size(ModGlue%Lin%z) > 0) .and. allocated(ModData%Lin%z)) ModGlue%Lin%z(iz:iz + ModData%Vars%Nz - 1) = ModData%Lin%z
         if ((size(ModGlue%Lin%u) > 0) .and. allocated(ModData%Lin%u)) ModGlue%Lin%u(iu:iu + ModData%Vars%Nu - 1) = ModData%Lin%u
         if ((size(ModGlue%Lin%y) > 0) .and. allocated(ModData%Lin%y)) ModGlue%Lin%y(iy:iy + ModData%Vars%Ny - 1) = ModData%Lin%y

         ! Copy module Jacobians into glue code Jacobians
         if ((size(ModGlue%Lin%dYdu) > 0) .and. allocated(ModData%Lin%dYdu)) ModGlue%Lin%dYdu(iy:iy + ModData%Vars%Ny - 1, iu:iu + ModData%Vars%Nu - 1) = ModData%Lin%dYdu
         if ((size(ModGlue%Lin%dXdu) > 0) .and. allocated(ModData%Lin%dXdu)) ModGlue%Lin%dXdu(ix:ix + ModData%Vars%Nx - 1, iu:iu + ModData%Vars%Nu - 1) = ModData%Lin%dXdu
         if ((size(ModGlue%Lin%dYdx) > 0) .and. allocated(ModData%Lin%dYdx)) ModGlue%Lin%dYdx(iy:iy + ModData%Vars%Ny - 1, ix:ix + ModData%Vars%Nx - 1) = ModData%Lin%dYdx
         if ((size(ModGlue%Lin%dXdx) > 0) .and. allocated(ModData%Lin%dXdx)) ModGlue%Lin%dXdx(ix:ix + ModData%Vars%Nx - 1, ix:ix + ModData%Vars%Nx - 1) = ModData%Lin%dXdx

         ! Increment starting index for next module
         ix = ix + ModData%Vars%Nx
         ixd = ixd + ModData%Vars%Nxd
         iz = iz + ModData%Vars%Nz
         iu = iu + ModData%Vars%Nu
         iy = iy + ModData%Vars%Ny

         ! If writing the module matrices was requested
         if (p_FAST%LinOutMod) then

            ! Assemble output file name based on module abbreviation
            ! If module is BeamDyn or more than one instance, include instance
            OutFileName = trim(LinRootName)//'.'//trim(ModData%Abbr)//".lin"
            if ((ModData%ID == Module_BD) .or. (count(Mods%ID == ModData%ID) > 1)) then
               OutFileName = trim(LinRootName)//'.'//trim(ModData%Abbr)//trim(Num2LStr(ModData%Ins))//".lin"
            end if

            ! Write linearization matrices
            call WriteModuleLinearMatrices(ModData, ModData%Vars%IdxLin, p_FAST, y_FAST, t_global, Un, OutFileName, ErrStat2, ErrMsg2)
            if (Failed()) return

         end if
         
         ! Check for NaNs or infinity in module Jacobian matrices
         if (allocated(ModData%Lin%dYdu)) then
            if (any(isnan(ModData%Lin%dYdu))) then
               ErrStat = ErrID_Fatal
               ErrMsg = 'NaNs detected in dYdu for module '//ModData%Abbr
               return
            end if
         end if
         if (allocated(ModData%Lin%dXdu)) then
            if (any(isnan(ModData%Lin%dXdu))) then
               ErrStat = ErrID_Fatal
               ErrMsg = 'NaNs detected in dXdu for module '//ModData%Abbr
               return
            end if
         end if
         if (allocated(ModData%Lin%dYdx)) then
            if (any(isnan(ModData%Lin%dYdx))) then
               ErrStat = ErrID_Fatal
               ErrMsg = 'NaNs detected in dYdx for module '//ModData%Abbr
               return
            end if
         end if
         if (allocated(ModData%Lin%dXdx)) then
            if (any(isnan(ModData%Lin%dXdx))) then
               ErrStat = ErrID_Fatal
               ErrMsg = 'NaNs detected in dXdx for module '//ModData%Abbr
               return
            end if
         end if

      end associate
   end do

   ! Linearize mesh mappings to popoulate dUdy and dUdu
   ModGlue%Lin%dUdy = 0.0_R8Ki
   call Eye2D(ModGlue%Lin%dUdu, ErrStat2, ErrMsg2); if (Failed()) return
   call FAST_LinearizeMappings(Turbine, Mods, m%Mappings, p%iMod, ErrStat2, ErrMsg2, ModGlue%Lin%dUdu, ModGlue%Lin%dUdy)
   if (Failed()) return

   ! Calculate the glue code state matrices (A, B, C, D)
   call ModLin_StateMatrices(ModGlue, real(p_FAST%UJacSclFact, R8Ki), ErrStat2, ErrMsg2)
   if (Failed()) return

   ! Write glue code data
   OutFileName = trim(LinRootName)//".lin"
   call WriteModuleLinearMatrices(ModGlue, ModGlue%Vars%IdxLin, p_FAST, y_FAST, t_global, Un, OutFileName, ErrStat2, ErrMsg2, IsGlue=.true.)
   if (Failed()) return

   ! Update index for next linearization time
   m_FAST%Lin%NextLinTimeIndx = m_FAST%Lin%NextLinTimeIndx + 1

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function Failed
end subroutine

!> ModLin_StateMatrices forms the full-system state matrices for linearization: A, B, C, and D.
!! Note that it uses LAPACK_GEMM instead of MATMUL for matrix multiplications because of stack-space issues (these
!! matrices get large quickly).
subroutine ModLin_StateMatrices(ModGlue, JacScaleFactor, ErrStat, ErrMsg)
   type(ModDataType), intent(inout) :: ModGlue        !< Glue module data
   real(R8Ki), intent(in)           :: JacScaleFactor !< Scale factor for conditioning the Jacobians
   integer(IntKi), intent(out)      :: ErrStat        !< Error status of the operation
   character(*), intent(out)        :: ErrMsg         !< Error message if ErrStat /= ErrID_None

   character(*), parameter          :: RoutineName = 'ModLin_StateMatrices'
   integer(IntKi)                   :: ErrStat2
   character(ErrMsgLen)             :: ErrMsg2
   real(R8Ki), allocatable          :: G(:, :), tmp(:, :), dUdu(:, :), dUdy(:, :)
   integer(IntKi), allocatable      :: ipiv(:)

   ! A = dXdx
   ! B = dXdu
   ! C = dYdx
   ! D = dYdu

   ! Create copies of dUdu and dUdy for calculating matrices
   call AllocAry(dUdu, size(ModGlue%Lin%dUdu, 1), size(ModGlue%Lin%dUdu, 2), 'dUdu', ErrStat2, ErrMsg2)
   call AllocAry(dUdy, size(ModGlue%Lin%dUdy, 1), size(ModGlue%Lin%dUdy, 2), 'dUdy', ErrStat2, ErrMsg2)
   dUdu = ModGlue%Lin%dUdu
   dUdy = ModGlue%Lin%dUdy

   ! *** get G matrix ****
   !----------------------
   if (.not. allocated(G)) then
      call AllocAry(G, size(dUdu, 1), size(dUdu, 2), 'G', ErrStat2, ErrMsg2)
      if (Failed()) return

      call AllocAry(ipiv, ModGlue%Vars%Nu, 'ipiv', ErrStat2, ErrMsg2)
      if (Failed()) return
   end if

   !G = dUdu + matmul( dUdy, y_FAST%Lin%Glue%D )
   G = dUdu
   call LAPACK_GEMM('N', 'N', 1.0_R8Ki, dUdy, ModGlue%Lin%dYdu, 1.0_R8Ki, G, ErrStat2, ErrMsg2)
   if (Failed()) return

   ! G can be ill-conditioned, so we are going to precondition with G_hat = S^(-1) * G * S
   ! we will also multiply the right-hand-side of the equations that need G inverse so that
   ! dUdy_hat = S^(-1)*dUdy and dUdu_hat = S^(-1)*dUdu
   call Precondition(ModGlue%Vars%u, G, dUdu, dUdy, JacScaleFactor)

   ! Form G_hat^(-1) * (S^-1*dUdy) and G^(-1) * (S^-1*dUdu)
   ! factor G for the two solves:
   call LAPACK_getrf(M=size(G, 1), N=size(G, 2), A=G, IPIV=ipiv, ErrStat=ErrStat2, ErrMsg=ErrMsg2)
   if (Failed()) return

   ! after the this solve, dUdy holds G_hat^(-1) * dUdy_hat:
   call LAPACK_getrs(trans='N', N=size(G, 2), A=G, IPIV=ipiv, B=dUdy, ErrStat=ErrStat2, ErrMsg=ErrMsg2)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   ! after the this solve, dUdu holds G_hat^(-1) * dUdu_hat:
   call LAPACK_getrs(trans='N', N=size(G, 2), A=G, IPIV=ipiv, B=dUdu, ErrStat=ErrStat2, ErrMsg=ErrMsg2)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   ! Deallocate G and ipiv because the solves are complete
   deallocate (G)
   deallocate (ipiv)

   ! after this call, dUdu holds G^(-1)*dUdu and dUdy holds G^(-1)*dUdy:
   call Postcondition(ModGlue%Vars%u, dUdu, dUdy, JacScaleFactor)

   ! Allocate tmp matrix for A and C calculations
   call AllocAry(tmp, ModGlue%Vars%Nu, ModGlue%Vars%Nx, 'G^-1*dUdy*C', ErrStat2, ErrMsg2)
   if (Failed()) return

   ! tmp = G^(-1) * dUdy * diag(C)
   call LAPACK_GEMM('N', 'N', 1.0_R8Ki, dUdy, ModGlue%Lin%dYdx, 0.0_R8Ki, tmp, ErrStat2, ErrMsg2)
   if (Failed()) return

   ! A
   ! dXdx = dXdx - matmul( dXdu, tmp )
   call LAPACK_GEMM('N', 'N', -1.0_R8Ki, ModGlue%Lin%dXdu, tmp, 1.0_R8Ki, ModGlue%Lin%dXdx, ErrStat2, ErrMsg2)
   if (Failed()) return

   ! C
   ! dYdx = dYdx - matmul( dYdu, tmp )
   call LAPACK_GEMM('N', 'N', -1.0_R8Ki, ModGlue%Lin%dYdu, tmp, 1.0_R8Ki, ModGlue%Lin%dYdx, ErrStat2, ErrMsg2)
   if (Failed()) return

   ! B
   if (Failed()) return
   tmp = ModGlue%Lin%dXdu
   ! dXdu = matmul( dXdu, dUdu )
   call LAPACK_GEMM('N', 'N', 1.0_R8Ki, tmp, dUdu, 0.0_R8Ki, ModGlue%Lin%dXdu, ErrStat2, ErrMsg2)
   if (Failed()) return

   ! D
   if (Failed()) return
   tmp = ModGlue%Lin%dYdu
   ! D = matmul( dYdu, dUdu )
   call LAPACK_GEMM('N', 'N', 1.0_R8Ki, tmp, dUdu, 0.0_R8Ki, ModGlue%Lin%dYdu, ErrStat2, ErrMsg2)
   if (Failed()) return

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function Failed
end subroutine

!> Precondition returns the preconditioned matrix, hat{G}, such that hat{G} = S^(-1) G S withS^(-1 defined
!! such that loads are scaled by p_FAST%UJacSclFact. It also returns the preconditioned matrices hat{dUdu} and
!! hat{dUdy} such that hat{dUdu} = S^(-1) dUdu and
!! hat{dUdy} = S^(-1) dUdy for the right-hand sides of the equations to be solved.
subroutine Precondition(uVars, G, dUdu, dUdy, JacScaleFactor)
   type(ModVarType), intent(in)  :: uVars(:)       !< Input variables from glue code
   real(R8Ki), intent(inout)     :: G(:, :)        !< variable for glue-code linearization (in is G; out is G_hat)
   real(R8Ki), intent(inout)     :: dUdu(:, :)     !< jacobian in FAST linearization from right-hand-side of equation
   real(R8Ki), intent(inout)     :: dUdy(:, :)     !< jacobian in FAST linearization from right-hand-side of equation
   real(R8Ki), intent(in)        :: JacScaleFactor !< jacobian scale factor
   real(R8Ki), allocatable       :: diag(:)        !< diagonal elements of G
   integer(IntKi)                :: i

   ! Copy diagonal of G into temporary array, to be restored after conditioning,
   ! this is done to avoid loss of precision in the diagonal terms
   allocate (diag(size(G, 1)))
   do i = 1, size(diag)
      diag(i) = G(i, i)
   end do

   ! Loop through glue code input varies
   do i = 1, size(uVars)

      ! If variable is not a load (force or moment), continue
      if (.not. MV_HasFlags(uVars(i), ior(VF_Force, VF_Moment))) cycle

      ! Otherwise get variable start and end indices in matrix
      associate (iLoc => uVars(i)%iLoc)

         ! Multiply columns of G
         G(:, iLoc(1):iLoc(2)) = G(:, iLoc(1):iLoc(2))*JacScaleFactor

         ! Divide rows of G
         G(iLoc(1):iLoc(2), :) = G(iLoc(1):iLoc(2), :)/JacScaleFactor

         ! Divide rows of dUdu
         dUdu(iLoc(1):iLoc(2), :) = dUdu(iLoc(1):iLoc(2), :)/JacScaleFactor

         ! Divide rows of dUdy
         dUdy(iLoc(1):iLoc(2), :) = dUdy(iLoc(1):iLoc(2), :)/JacScaleFactor

      end associate
   end do

   ! Restore diagonal of G from temporary array
   do i = 1, size(diag)
      G(i, i) = diag(i)
   end do

end subroutine

!> This routine returns the matrices tilde{dUdu} and tilde{dUdy} such that
!! tilde{dUdu} = G^(-1) dUdu and
!! tilde{dUdy} = G^(-1) dUdy, which have been solved using the preconditioned system defined in fast_lin::precondition.
subroutine Postcondition(uVars, dUdu, dUdy, JacScaleFactor)
   type(ModVarType), intent(in)  :: uVars(:)       !< Input variables from glue code
   real(R8Ki), intent(in)        :: JacScaleFactor !< jacobian scale factor
   real(R8Ki), intent(inout)     :: dUdu(:, :)     !< jacobian in FAST linearization from right-hand-side of equation
   real(R8Ki), intent(inout)     :: dUdy(:, :)     !< jacobian in FAST linearization from right-hand-side of equation
   integer(IntKi)                :: i

   ! Loop through glue code input varies
   do i = 1, size(uVars)

      ! If variable is not a load (force or moment), continue
      if (.not. MV_HasFlags(uVars(i), ior(VF_Force, VF_Moment))) cycle

      ! Otherwise get variable start and end indices in matrix
      associate (iLoc => uVars(i)%iLoc)

         ! Multiply rows of dUdu
         dUdu(iLoc(1):iLoc(2), :) = dUdu(iLoc(1):iLoc(2), :)*JacScaleFactor

         ! Multiply rows of dUdy
         dUdy(iLoc(1):iLoc(2), :) = dUdy(iLoc(1):iLoc(2), :)*JacScaleFactor

      end associate
   end do

end subroutine

subroutine WriteModuleLinearMatrices(ModData, VarIdx, p_FAST, y_FAST, t_global, Un, OutFileName, ErrStat, ErrMsg, IsGlue)

   type(ModDataType), intent(in)    :: ModData        !< Module data
   type(VarsIdxType), intent(in)    :: VarIdx         !< Variable index
   type(FAST_ParameterType)         :: p_FAST         !< Parameters
   type(FAST_OutputFileType)        :: y_FAST         !< Output variables
   real(DbKi), intent(in)           :: t_global       !< current time step (written in file)
   integer(IntKi), intent(out)      :: Un             !< Unit number for file
   character(*), intent(in)         :: OutFileName    !< output file name
   integer(IntKi), intent(out)      :: ErrStat        !< Error status of the operation
   character(*), intent(out)        :: ErrMsg         !< Error message if ErrStat /= ErrID_None
   logical, optional                :: IsGlue         !< Flag indicating this is writing glue code matrices

   character(*), parameter          :: RoutineName = 'WriteModuleLinearMatrices'
   integer(IntKi)                   :: ErrStat2
   character(ErrMsgLen)             :: ErrMsg2
   character(32)                    :: Desc
   integer(IntKi)                   :: i
   character(50)                    :: Fmt
   logical, allocatable             :: uUse(:), yUse(:)
   logical                          :: IsGlueLoc

   ErrStat = ErrID_None
   ErrMsg = ""

   ! Set local flag for if glue code matrices are being written
   IsGlueLoc = .false.
   if (present(IsGlue)) IsGlueLoc = IsGlue

   ! Open linearization file
   call OpenFOutFile(Un, OutFileName, ErrStat2, ErrMsg2); if (Failed()) return

   !----------------------------------------------------------------------------
   ! Header
   !----------------------------------------------------------------------------

   write (Un, '(/,A)') 'Linearized model: '//trim(y_FAST%FileDescLines(1))
   write (Un, '(1X,A,/)') trim(y_FAST%FileDescLines(2))
   write (Un, '(A,/)') trim(y_FAST%FileDescLines(3))

   write (Un, '(A)') 'Simulation information:'

   fmt = '(3x,A,1x,'//trim(p_FAST%OutFmt_t)//',1x,A)'
   Desc = 'Simulation time:'; write (Un, fmt) Desc, t_global, 's'
   Desc = 'Rotor Speed:    '; write (Un, fmt) Desc, y_FAST%Lin%RotSpeed, 'rad/s'
   Desc = 'Azimuth:        '; write (Un, fmt) Desc, y_FAST%Lin%Azimuth, 'rad'
   Desc = 'Wind Speed:     '; write (Un, fmt) Desc, y_FAST%Lin%WindSpeed, 'm/s'

   fmt = '(3x,A,1x,I5)'
   Desc = 'Number of continuous states: '; write (Un, fmt) Desc, VarIdx%Nx
   Desc = 'Number of discrete states:   '; write (Un, fmt) Desc, VarIdx%Nxd
   Desc = 'Number of constraint states: '; write (Un, fmt) Desc, VarIdx%Nz
   Desc = 'Number of inputs:            '; write (Un, fmt) Desc, VarIdx%Nu
   Desc = 'Number of outputs:           '; write (Un, fmt) Desc, VarIdx%Ny

   Desc = 'Jacobians included in this file?'
   fmt = '(3x,A,1x,A5)'
   if (p_FAST%LinOutJac) then
      write (Un, fmt) Desc, 'Yes'
   else
      write (Un, fmt) Desc, 'No'
   end if

   write (Un, '()')    !print a blank line

   if (VarIdx%Nx > 0) then
      write (Un, '(A)') 'Order of continuous states:'
      call WrLinFile_txt_Table(ModData%Vars%x, VarIdx%FlagFilter, p_FAST, Un, "Row/Column", ModData%Lin%x)

      write (Un, '(A)') 'Order of continuous state derivatives:'
      call WrLinFile_txt_Table(ModData%Vars%x, VarIdx%FlagFilter, p_FAST, Un, "Row/Column", ModData%Lin%dx, IsDeriv=.true.)
   end if

   if (VarIdx%Nxd > 0) then
      write (Un, '(A)') 'Order of discrete states:'
      call WrLinFile_txt_Table(ModData%Vars%xd, VarIdx%FlagFilter, p_FAST, Un, "Row/Column", ModData%Lin%xd)
   end if

   if (VarIdx%Nz > 0) then
      write (Un, '(A)') 'Order of constraint states:'
      call WrLinFile_txt_Table(ModData%Vars%z, VarIdx%FlagFilter, p_FAST, Un, "Row/Column", ModData%Lin%z)
   end if

   if (VarIdx%Nu > 0) then
      write (Un, '(A)') 'Order of inputs:'
      call WrLinFile_txt_Table(ModData%Vars%u, VarIdx%FlagFilter, p_FAST, Un, "Column  ", ModData%Lin%u, ShowRot=.true.)
   end if

   if (VarIdx%Ny > 0) then
      write (Un, '(A)') 'Order of outputs:'
      call WrLinFile_txt_Table(ModData%Vars%y, VarIdx%FlagFilter, p_FAST, Un, "Row  ", ModData%Lin%y, ShowRot=.true.)
   end if

   allocate (uUse(ModData%Vars%Nu))
   uUse = .false.
   uUse(VarIdx%iu) = .true.

   allocate (yUse(ModData%Vars%Ny))
   yUse = .false.
   yUse(VarIdx%iy) = .true.

   if (p_FAST%LinOutJac) then
      write (Un, '(/,A,/)') 'Jacobian matrices:'
      if (IsGlueLoc) then
         call WrPartialMatrix(ModData%Lin%dUdu, Un, p_FAST%OutFmt, 'dUdu', UseRow=uUse, UseCol=uUse)
         call WrPartialMatrix(ModData%Lin%dUdy, Un, p_FAST%OutFmt, 'dUdy', UseRow=uUse, UseCol=yUse)
      else
         if (allocated(ModData%Lin%dXdx)) call WrPartialMatrix(ModData%Lin%dXdx, Un, p_FAST%OutFmt, 'dXdx')
         if (allocated(ModData%Lin%dXdu)) call WrPartialMatrix(ModData%Lin%dXdu, Un, p_FAST%OutFmt, 'dXdu', UseCol=uUse)
         if (allocated(ModData%Lin%dYdx)) call WrPartialMatrix(ModData%Lin%dYdx, Un, p_FAST%OutFmt, 'dYdx', UseRow=yUse)
         if (allocated(ModData%Lin%dYdu)) call WrPartialMatrix(ModData%Lin%dYdu, Un, p_FAST%OutFmt, 'dYdu', UseRow=yUse, UseCol=uUse)
      end if
   end if

   write (Un, '(/,A,/)') 'Linearized state matrices:'
   if (allocated(ModData%Lin%dXdx)) call WrPartialMatrix(ModData%Lin%dXdx, Un, p_FAST%OutFmt, 'A')
   if (allocated(ModData%Lin%dXdu)) call WrPartialMatrix(ModData%Lin%dXdu, Un, p_FAST%OutFmt, 'B', UseCol=uUse)
   if (allocated(ModData%Lin%dYdx)) call WrPartialMatrix(ModData%Lin%dYdx, Un, p_FAST%OutFmt, 'C', UseRow=yUse)
   if (allocated(ModData%Lin%dYdu)) call WrPartialMatrix(ModData%Lin%dYdu, Un, p_FAST%OutFmt, 'D', UseRow=yUse, UseCol=uUse)
   if (allocated(ModData%Lin%StateRotation)) call WrPartialMatrix(ModData%Lin%StateRotation, Un, p_FAST%OutFmt, 'StateRotation')

   ! Close file
   close (Un)

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
      if (Failed) close (Un)
   end function Failed
end subroutine WriteModuleLinearMatrices

subroutine WrLinFile_txt_Table(VarAry, FlagFilter, p_FAST, Un, RowCol, op, IsDeriv, ShowRot)

   type(ModVarType), intent(in)  :: VarAry(:)   !< variable array
   integer(IntKi), intent(in)    :: FlagFilter  !< unit number
   type(FAST_ParameterType)      :: p_FAST      !< Parameters
   integer(IntKi), intent(in)    :: Un          !< unit number
   character(*), intent(in)      :: RowCol      !< Row/Column description
   real(R8Ki), intent(in)        :: op(:)       !< operating point values (possibly different size that Desc because of orientations)
   logical, optional, intent(in) :: IsDeriv     !< flag that tells us if we need to modify the channel names for derivatives (xdot)
   logical, optional, intent(in) :: ShowRot     !< flag to show rotation matrix if field is orientation

   character(*), parameter       :: RoutineName = 'WrLinFile_txt_Table'
   integer(IntKi)                :: TS             ! Tab stop column
   integer(IntKi)                :: i_op           ! Index of value in operating piont
   logical                       :: IsDerivLoc     ! flag that tells us if we need to modify the channel names for derivatives (xdot)
   logical                       :: VarRotFrame    ! flag that tells us if this column is in the rotating frame
   integer(IntKi)                :: VarDerivOrder  ! integer indicating the maximum time-derivative order of a channel (this will be 0 for anything that is not a continuous state)
   character(100)                :: Fmt, FmtStr, FmtRot
   character(25)                 :: DerivStr, DerivUnitStr
   logical                       :: ShowRotLoc
   real(R8Ki)                    :: DCM(3, 3)
   integer(IntKi)                :: i, j, RowColIdx

   ShowRotLoc = .false.
   if (present(ShowRot)) ShowRotLoc = ShowRot

   IsDerivLoc = .false.
   if (present(IsDeriv)) IsDerivLoc = IsDeriv

   if (IsDerivLoc) then
      if (p_FAST%CompAeroMaps .and. p_FAST%CompElast /= MODULE_BD) then ! this might not work if we are using some other (not BD, ED) module with states
         DerivStr = 'Second time derivative of'
         DerivUnitStr = '/s^2'
      else
         DerivStr = 'First time derivative of'
         DerivUnitStr = '/s'
      end if
   else
      DerivStr = ''
      DerivUnitStr = ''
   end if

   ! tab stop after operating point
   TS = 14 + 3*p_FAST%FmtWidth + 7

   ! Construct write formats
   Fmt = '(3x,I8,3x,'//trim(p_FAST%OutFmt)//',T'//trim(Num2LStr(TS))//',L8,8x,I8,9x,A)'
   FmtRot = '(3x,I8,3x,'//trim(p_FAST%OutFmt)//',2(", ",'//trim(p_FAST%OutFmt)//'),T'//trim(Num2LStr(TS))//',L8,8x,I8,9x,A)'
   FmtStr = '(3x,A10,1x,A,T'//trim(Num2LStr(TS))//',A15,1x,A16,1x,A)'

   ! Write header
   write (Un, FmtStr) RowCol, 'Operating Point', 'Rotating Frame?', 'Derivative Order', 'Description'
   write (Un, FmtStr) '----------', '---------------', '---------------', '----------------', '-----------'

   ! Loop through variables in array
   RowColIdx = 0
   do i = 1, size(VarAry)
      associate (Var => VarAry(i))

         ! If variable does not have the filter flag, continue
         if (.not. MV_HasFlags(Var, FlagFilter)) cycle

         ! Is variable in the rotating frame?
         VarRotFrame = MV_HasFlags(Var, VF_RotFrame)

         ! Get variable derivative order
         if (MV_HasFlags(Var, VF_DerivOrder2)) then
            VarDerivOrder = 2
         else if (MV_HasFlags(Var, VF_DerivOrder1)) then
            VarDerivOrder = 1
         else
            VarDerivOrder = 0
         end if

         ! Loop through values in variable
         do j = 1, Var%Num

            ! Increment value counter
            RowColIdx = RowColIdx + 1

            ! Index in operating point array
            i_op = Var%iLoc(1) + j - 1

            ! If variable is orientation and show rotation matrix flag is true
            if (ShowRotLoc .and. (Var%Field == VF_Orientation)) then

               ! Skip writing if not the first value in orientation (3 values)
               if (mod(j - 1, 3) /= 0) cycle

               ! Convert WM parameters to DCM
               DCM = wm_to_dcm(real(op(i_op:i_op + 2), R8Ki))

               ! Write 3 rows of data (full dcm)
               write (Un, FmtRot) RowColIdx + 0, dcm(1, 1), dcm(1, 2), dcm(1, 3), VarRotFrame, VarDerivOrder, trim(Var%LinNames(j + 0))
               write (Un, FmtRot) RowColIdx + 1, dcm(2, 1), dcm(2, 2), dcm(2, 3), VarRotFrame, VarDerivOrder, trim(Var%LinNames(j + 1))
               write (Un, FmtRot) RowColIdx + 2, dcm(3, 1), dcm(3, 2), dcm(3, 3), VarRotFrame, VarDerivOrder, trim(Var%LinNames(j + 2))

            else if (IsDerivLoc) then
               write (Un, Fmt) RowColIdx, op(i_op), VarRotFrame, VarDerivOrder, trim(DerivStr)//' '//trim(Var%LinNames(j))//trim(DerivUnitStr)
            else
               write (Un, Fmt) RowColIdx, op(i_op), VarRotFrame, VarDerivOrder, trim(Var%LinNames(j))
            end if

         end do
      end associate
   end do

   write (Un, '()')    !print a blank line

end subroutine WrLinFile_txt_Table

end module
