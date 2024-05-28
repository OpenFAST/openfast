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
module FAST_ModGlue

use NWTC_Library
use NWTC_LAPACK

use FAST_ModTypes
use FAST_Types
use FAST_Funcs
use FAST_Mapping

implicit none

private
public :: ModGlue_Init, ModGlue_Linearize_OP, MV_AddModule

contains

subroutine ModGlue_Init(p, m, y, p_FAST, m_FAST, Turbine, ErrStat, ErrMsg)

   type(Glue_ParameterType), intent(inout)         :: p        !< Glue Parameters
   type(Glue_MiscVarType), intent(inout)           :: m        !< Glue MiscVars
   type(Glue_OutputFileType), intent(inout)        :: y        !< Glue Output
   type(FAST_ParameterType), intent(inout)         :: p_FAST   !< FAST Parameters
   type(FAST_MiscVarType), intent(inout)           :: m_FAST   !< FAST MiscVars
   type(FAST_TurbineType), intent(inout)           :: Turbine
   integer(IntKi), intent(out)                     :: ErrStat
   character(*), intent(out)                       :: ErrMsg

   character(*), parameter                         :: RoutineName = 'ModLin_Init'
   integer(IntKi)                                  :: ErrStat2
   character(ErrMsgLen)                            :: ErrMsg2
   integer(IntKi), allocatable                     :: modIDs(:), modIdx(:)
   integer(IntKi)                                  :: i, j, k
   integer(IntKi)                                  :: FlagFilters
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
   if (.not. allocated(m%ModData)) then
      call SetErrStat(ErrID_Fatal, "No modules were used", ErrStat, ErrMsg, RoutineName)
      return
   end if

   ! Create array of indices for Mods array
   modIdx = [(i, i=1, size(m%ModData))]

   ! Get array of module IDs
   modIDs = [(m%ModData(i)%ID, i=1, size(m%ModData))]

   ! Establish module index order for linearization
   allocate (p%iMod(0))
   do i = 1, size(LinMods)
      p%iMod = [p%iMod, pack(modIdx, ModIDs == LinMods(i))]
   end do

   ! Loop through modules, if module is not in index, return with error
   do i = 1, size(m%ModData)
      if (.not. any(i == p%iMod)) then
         call SetErrStat(ErrID_Fatal, "Module "//trim(m%ModData(i)%Abbr)//" not supported in linearization", &
                         ErrStat, ErrMsg, RoutineName)
         return
      end if
   end do

   !----------------------------------------------------------------------------
   ! Glue Module Variables
   !----------------------------------------------------------------------------

   ! Allocate variable structure for glue
   allocate (y%ModGlue%Vars)

   ! Initialize number of values in each variable group
   y%ModGlue%Vars%Nx = 0
   y%ModGlue%Vars%Nxd = 0
   y%ModGlue%Vars%Nz = 0
   y%ModGlue%Vars%Nu = 0
   y%ModGlue%Vars%Ny = 0

   ! Allocate arrays for glue variables
   allocate (y%ModGlue%Vars%x(0), y%ModGlue%Vars%xd(0), y%ModGlue%Vars%z(0), y%ModGlue%Vars%u(0), y%ModGlue%Vars%y(0))

   ! Loop through each module by index
   do i = 1, size(p%iMod)
      associate (ModData => m%ModData(p%iMod(i)))

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
         y%ModGlue%Vars%Nx = y%ModGlue%Vars%Nx + ModData%Vars%Nx

         ! Save start index of module variables and append to glue code variables
         k = size(y%ModGlue%Vars%x) + 1
         y%ModGlue%Vars%x = [y%ModGlue%Vars%x, ModData%Vars%x]

         ! Loop through added variables and add name prefix to linearization names
         call AddLinNamePrefix(y%ModGlue%Vars%x(k:), NamePrefix)

         !----------------------------------------------------------------------
         ! Module discrete state variables
         !----------------------------------------------------------------------

         ! Set module data start index in global arrays, increment data size
         y%ModGlue%Vars%Nxd = y%ModGlue%Vars%Nxd + ModData%Vars%Nxd

         ! Save start index of module variables and append to glue code variables
         k = size(y%ModGlue%Vars%xd) + 1
         y%ModGlue%Vars%xd = [y%ModGlue%Vars%xd, ModData%Vars%xd]

         ! Loop through added variables and add name prefix to linearization names
         call AddLinNamePrefix(y%ModGlue%Vars%xd(k:), NamePrefix)

         !----------------------------------------------------------------------
         ! Module constraint state variables
         !----------------------------------------------------------------------

         ! Set module data start index in global arrays, increment data size
         y%ModGlue%Vars%Nz = y%ModGlue%Vars%Nz + ModData%Vars%Nz

         ! Save start index of module variables and append to glue code variables
         k = size(y%ModGlue%Vars%z) + 1
         y%ModGlue%Vars%z = [y%ModGlue%Vars%z, ModData%Vars%z]

         ! Loop through added variables and add name prefix to linearization names
         call AddLinNamePrefix(y%ModGlue%Vars%z(k:), NamePrefix)

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
            ! For standard inputs, use VF_Linearize flag as set in the module
         case (LIN_ALL)
            do j = 1, size(ModData%Vars%u)
               call MV_SetFlags(ModData%Vars%u(j), VF_Linearize)
            end do
         end select

         ! Set module data start index in global arrays, increment data size
         y%ModGlue%Vars%Nu = y%ModGlue%Vars%Nu + ModData%Vars%Nu

         ! Save start index of module variables and append to glue code variables
         k = size(y%ModGlue%Vars%u) + 1
         y%ModGlue%Vars%u = [y%ModGlue%Vars%u, ModData%Vars%u]

         ! Loop through added variables and add name prefix to linearization names
         call AddLinNamePrefix(y%ModGlue%Vars%u(k:), NamePrefix)

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
         y%ModGlue%Vars%Ny = y%ModGlue%Vars%Ny + ModData%Vars%Ny

         ! Save start index of module variables and append to glue code variables
         k = size(y%ModGlue%Vars%y) + 1
         y%ModGlue%Vars%y = [y%ModGlue%Vars%y, ModData%Vars%y]

         ! Loop through added variables and add name prefix to linearization names
         call AddLinNamePrefix(y%ModGlue%Vars%y(k:), NamePrefix)

      end associate
   end do

   ! Calculate number of values in each group and set data location index
   call CalcVarDataLoc(y%ModGlue%Vars%x, y%ModGlue%Vars%Nx)
   call CalcVarDataLoc(y%ModGlue%Vars%xd, y%ModGlue%Vars%Nxd)
   call CalcVarDataLoc(y%ModGlue%Vars%z, y%ModGlue%Vars%Nz)
   call CalcVarDataLoc(y%ModGlue%Vars%u, y%ModGlue%Vars%Nu)
   call CalcVarDataLoc(y%ModGlue%Vars%y, y%ModGlue%Vars%Ny)

   !----------------------------------------------------------------------------
   ! Mesh Mapping
   !----------------------------------------------------------------------------

   call FAST_InitMappings(m%ModData, m%Mappings, Turbine, ErrStat2, ErrMsg2); if (Failed()) return

   !----------------------------------------------------------------------------
   ! Allocate linearization arrays and matrices
   !----------------------------------------------------------------------------

   ! If linearization is enabled
   if (p_FAST%Linearize) then

      ! Initialize linearization index
      call Idx_Init(m%ModData, p%iMod, p%IdxLin, VF_None, ErrStat2, ErrMsg2); if (Failed()) return

      ! Allocate linearization arrays
      call AllocAry(y%ModGlue%Lin%x, y%ModGlue%Vars%Nx, "x", ErrStat2, ErrMsg2); if (Failed()) return
      call AllocAry(y%ModGlue%Lin%dx, y%ModGlue%Vars%Nx, "dx", ErrStat2, ErrMsg2); if (Failed()) return
      call AllocAry(y%ModGlue%Lin%xd, y%ModGlue%Vars%Nxd, "xd", ErrStat2, ErrMsg2); if (Failed()) return
      call AllocAry(y%ModGlue%Lin%z, y%ModGlue%Vars%Nz, "z", ErrStat2, ErrMsg2); if (Failed()) return
      call AllocAry(y%ModGlue%Lin%u, y%ModGlue%Vars%Nu, "u", ErrStat2, ErrMsg2); if (Failed()) return
      call AllocAry(y%ModGlue%Lin%y, y%ModGlue%Vars%Ny, "y", ErrStat2, ErrMsg2); if (Failed()) return

      ! Allocate full Jacobian matrices
      call AllocAry(y%ModGlue%Lin%dYdu, y%ModGlue%Vars%Ny, y%ModGlue%Vars%Nu, "dYdu", ErrStat2, ErrMsg2); if (Failed()) return
      call AllocAry(y%ModGlue%Lin%dXdu, y%ModGlue%Vars%Nx, y%ModGlue%Vars%Nu, "dXdu", ErrStat2, ErrMsg2); if (Failed()) return
      call AllocAry(y%ModGlue%Lin%dYdx, y%ModGlue%Vars%Ny, y%ModGlue%Vars%Nx, "dYdx", ErrStat2, ErrMsg2); if (Failed()) return
      call AllocAry(y%ModGlue%Lin%dXdx, y%ModGlue%Vars%Nx, y%ModGlue%Vars%Nx, "dXdx", ErrStat2, ErrMsg2); if (Failed()) return
      call AllocAry(y%ModGlue%Lin%dUdu, y%ModGlue%Vars%Nu, y%ModGlue%Vars%Nu, "dUdu", ErrStat2, ErrMsg2); if (Failed()) return
      call AllocAry(y%ModGlue%Lin%dUdy, y%ModGlue%Vars%Nu, y%ModGlue%Vars%Ny, "dUdy", ErrStat2, ErrMsg2); if (Failed()) return

      ! Initialize arrays to store operating point states and input
      call AllocAry(y%OP%x, y%ModGlue%Vars%Nx, p_FAST%NLinTimes, "x", ErrStat2, ErrMsg2); if (Failed()) return
      call AllocAry(y%OP%xd, y%ModGlue%Vars%Nxd, p_FAST%NLinTimes, "xd", ErrStat2, ErrMsg2); if (Failed()) return
      call AllocAry(y%OP%z, y%ModGlue%Vars%Nz, p_FAST%NLinTimes, "z", ErrStat2, ErrMsg2); if (Failed()) return
      call AllocAry(y%OP%u, y%ModGlue%Vars%Nu, p_FAST%NLinTimes, "u", ErrStat2, ErrMsg2); if (Failed()) return
   end if

contains

   subroutine CalcVarDataLoc(VarAry, DataSize)
      type(ModVarType), intent(inout)  :: VarAry(:)
      integer(IntKi), intent(out)      :: DataSize
      DataSize = 0
      do i = 1, size(VarAry)
         VarAry(i)%iLoc = [DataSize + 1, DataSize + VarAry(i)%Num]
         DataSize = DataSize + VarAry(i)%Num
      end do
   end subroutine

   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function Failed

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

subroutine ModGlue_Linearize_OP(Turbine, p, m, y, p_FAST, m_FAST, y_FAST, t_global, ErrStat, ErrMsg)

   type(Glue_ParameterType), intent(inout)   :: p        !< Glue parameters
   type(Glue_MiscVarType), intent(inout)     :: m        !< Glue MiscVars
   type(Glue_OutputFileType), intent(inout)  :: y        !< Glue Output
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
   y%ModGlue%Lin%dYdu = 0.0_R8Ki
   y%ModGlue%Lin%dXdu = 0.0_R8Ki
   y%ModGlue%Lin%dYdx = 0.0_R8Ki
   y%ModGlue%Lin%dXdx = 0.0_R8Ki

   ! Loop through modules by index
   do i = 1, size(p%iMod)
      associate (ModData => m%ModData(p%iMod(i)))

         ! Derivatives wrt input
         call FAST_JacobianPInput(ModData, t_global, STATE_CURR, Turbine, ErrStat2, ErrMsg2, &
                                  dYdu=ModData%Lin%dYdu, dXdu=ModData%Lin%dXdu)
         if (Failed()) return

         ! Derivatives wrt continuous state
         call FAST_JacobianPContState(ModData, t_global, STATE_CURR, Turbine, ErrStat2, ErrMsg2, &
                                      dYdx=ModData%Lin%dYdx, dXdx=ModData%Lin%dXdx, &
                                      StateRotation=ModData%Lin%StateRotation)
         if (Failed()) return

         ! Operating point values (must come after Jacobian routines because
         ! some modules calculate OP in those routines [MD])
         call FAST_GetOP(ModData, t_global, STATE_CURR, Turbine, ErrStat2, ErrMsg2, &
                         u_op=ModData%Lin%u, y_op=ModData%Lin%y, &
                         x_op=ModData%Lin%x, dx_op=ModData%Lin%dx)
         if (Failed()) return

         ! Copy module linearization arrays into glue linearization arrays
         if ((size(y%ModGlue%Lin%x) > 0) .and. allocated(ModData%Lin%x)) y%ModGlue%Lin%x(ix:ix + ModData%Vars%Nx - 1) = ModData%Lin%x
         if ((size(y%ModGlue%Lin%dx) > 0) .and. allocated(ModData%Lin%dx)) y%ModGlue%Lin%dx(ix:ix + ModData%Vars%Nx - 1) = ModData%Lin%dx
         if ((size(y%ModGlue%Lin%xd) > 0) .and. allocated(ModData%Lin%xd)) y%ModGlue%Lin%xd(ixd:ixd + ModData%Vars%Nxd - 1) = ModData%Lin%xd
         if ((size(y%ModGlue%Lin%z) > 0) .and. allocated(ModData%Lin%z)) y%ModGlue%Lin%z(iz:iz + ModData%Vars%Nz - 1) = ModData%Lin%z
         if ((size(y%ModGlue%Lin%u) > 0) .and. allocated(ModData%Lin%u)) y%ModGlue%Lin%u(iu:iu + ModData%Vars%Nu - 1) = ModData%Lin%u
         if ((size(y%ModGlue%Lin%y) > 0) .and. allocated(ModData%Lin%y)) y%ModGlue%Lin%y(iy:iy + ModData%Vars%Ny - 1) = ModData%Lin%y

         ! Copy module Jacobians into glue code Jacobian
         if ((size(y%ModGlue%Lin%dYdu) > 0) .and. allocated(ModData%Lin%dYdu)) y%ModGlue%Lin%dYdu(iy:iy + ModData%Vars%Ny - 1, iu:iu + ModData%Vars%Nu - 1) = ModData%Lin%dYdu
         if ((size(y%ModGlue%Lin%dXdu) > 0) .and. allocated(ModData%Lin%dXdu)) y%ModGlue%Lin%dXdu(ix:ix + ModData%Vars%Nx - 1, iu:iu + ModData%Vars%Nu - 1) = ModData%Lin%dXdu
         if ((size(y%ModGlue%Lin%dYdx) > 0) .and. allocated(ModData%Lin%dYdx)) y%ModGlue%Lin%dYdx(iy:iy + ModData%Vars%Ny - 1, ix:ix + ModData%Vars%Nx - 1) = ModData%Lin%dYdx
         if ((size(y%ModGlue%Lin%dXdx) > 0) .and. allocated(ModData%Lin%dXdx)) y%ModGlue%Lin%dXdx(ix:ix + ModData%Vars%Nx - 1, ix:ix + ModData%Vars%Nx - 1) = ModData%Lin%dXdx

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
            if ((ModData%ID == Module_BD) .or. (count(m%ModData%ID == ModData%ID) > 1)) then
               OutFileName = trim(LinRootName)//'.'//trim(ModData%Abbr)//trim(Num2LStr(ModData%Ins))//".lin"
            end if

            ! Write linearization matrices
            call CalcWriteLinearMatrices(ModData, p_FAST, y_FAST, t_global, Un, OutFileName, .false., ErrStat2, ErrMsg2)
            if (Failed()) return

         end if

         ! Check for NaNs or infinity in module Jacobian matrices
         if (JacobianHasNaNs(ModData%Lin%dYdu, "dYdu", ModData%Abbr)) return
         if (JacobianHasNaNs(ModData%Lin%dXdu, "dXdu", ModData%Abbr)) return
         if (JacobianHasNaNs(ModData%Lin%dYdx, "dYdx", ModData%Abbr)) return
         if (JacobianHasNaNs(ModData%Lin%dXdx, "dXdx", ModData%Abbr)) return

         ! Copy arrays into linearization operating points
         if (size(y%ModGlue%Lin%x) > 0) y%OP%x(:,m_FAST%Lin%NextLinTimeIndx) = y%ModGlue%Lin%x
         if (size(y%ModGlue%Lin%xd) > 0) y%OP%xd(:,m_FAST%Lin%NextLinTimeIndx) = y%ModGlue%Lin%xd
         if (size(y%ModGlue%Lin%z) > 0) y%OP%z(:,m_FAST%Lin%NextLinTimeIndx) = y%ModGlue%Lin%z
         if (size(y%ModGlue%Lin%u) > 0) y%OP%u(:,m_FAST%Lin%NextLinTimeIndx) = y%ModGlue%Lin%u

      end associate
   end do

   ! Linearize mesh mappings to populate dUdy and dUdu
   y%ModGlue%Lin%dUdy = 0.0_R8Ki
   call Eye2D(y%ModGlue%Lin%dUdu, ErrStat2, ErrMsg2); if (Failed()) return
   call FAST_LinearizeMappings(Turbine, m%ModData, m%Mappings, p%iMod, p%IdxLin, ErrStat2, ErrMsg2, y%ModGlue%Lin%dUdu, y%ModGlue%Lin%dUdy)
   if (Failed()) return

   ! Write glue code matrices to file
   OutFileName = trim(LinRootName)//".lin"
   call CalcWriteLinearMatrices(y%ModGlue, p_FAST, y_FAST, t_global, Un, OutFileName, .true., ErrStat2, ErrMsg2)
   if (Failed()) return

   ! Update index for next linearization time
   m_FAST%Lin%NextLinTimeIndx = m_FAST%Lin%NextLinTimeIndx + 1

contains
   logical function JacobianHasNaNs(Jac, label, abbr)
      real(R8Ki), allocatable, intent(in) :: Jac(:,:)
      character(*), intent(in)            :: label, abbr
      JacobianHasNaNs = .false.
      if (.not. allocated(Jac)) return
      if (size(Jac) == 0) return
      if (.not. any(isnan(Jac))) return
      ErrStat = ErrID_Fatal
      ErrMsg = 'NaNs detected in dXdx for module '//abbr
      JacobianHasNaNs = .true.
   end function
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function Failed
end subroutine

!> ModLin_StateMatrices forms the full-system state matrices for linearization: A, B, C, and D.
!! Note that it uses LAPACK_GEMM instead of MATMUL for matrix multiplications because of stack-space issues (these
!! matrices get large quickly).
subroutine CalcGlueStateMatrices(ModGlue, JacScaleFactor, ErrStat, ErrMsg)
   type(ModDataType), intent(inout) :: ModGlue        !< Glue module data
   real(R8Ki), intent(in)           :: JacScaleFactor !< Scale factor for conditioning the Jacobians
   integer(IntKi), intent(out)      :: ErrStat        !< Error status of the operation
   character(*), intent(out)        :: ErrMsg         !< Error message if ErrStat /= ErrID_None

   character(*), parameter          :: RoutineName = 'ModLin_StateMatrices'
   integer(IntKi)                   :: ErrStat2
   character(ErrMsgLen)             :: ErrMsg2
   real(R8Ki), allocatable          :: G(:, :), tmp(:, :)
   integer(IntKi), allocatable      :: ipiv(:)

   ! A = dXdx
   ! B = dXdu
   ! C = dYdx
   ! D = dYdu

   ! call DumpMatrix(1000, "dUdu.bin", ModGlue%Lin%dUdu, ErrStat2, ErrMsg2); if (Failed()) return
   ! call DumpMatrix(1000, "dUdy.bin", ModGlue%Lin%dUdy, ErrStat2, ErrMsg2); if (Failed()) return
   ! call DumpMatrix(1000, "A.bin", ModGlue%Lin%dXdx, ErrStat2, ErrMsg2); if (Failed()) return
   ! call DumpMatrix(1000, "B.bin", ModGlue%Lin%dXdu, ErrStat2, ErrMsg2); if (Failed()) return
   ! call DumpMatrix(1000, "C.bin", ModGlue%Lin%dYdx, ErrStat2, ErrMsg2); if (Failed()) return
   ! call DumpMatrix(1000, "D.bin", ModGlue%Lin%dYdu, ErrStat2, ErrMsg2); if (Failed()) return

   ! *** get G matrix ****
   !----------------------
   call AllocAry(G, size(ModGlue%Lin%dUdu, 1), size(ModGlue%Lin%dUdu, 2), 'G', ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(ipiv, ModGlue%Vars%Nu, 'ipiv', ErrStat2, ErrMsg2); if (Failed()) return

   !G = dUdu + matmul(dUdy, y_FAST%Lin%Glue%D)
   G = ModGlue%Lin%dUdu
   call LAPACK_GEMM('N', 'N', 1.0_R8Ki, ModGlue%Lin%dUdy, ModGlue%Lin%dYdu, 1.0_R8Ki, G, ErrStat2, ErrMsg2); if (Failed()) return

   ! G can be ill-conditioned, so we are going to precondition with G_hat = S^(-1) * G * S
   ! we will also multiply the right-hand-side of the equations that need G inverse so that
   ! dUdy_hat = S^(-1)*dUdy and dUdu_hat = S^(-1)*dUdu
   call Precondition(ModGlue%Vars%u, G, ModGlue%Lin%dUdu, ModGlue%Lin%dUdy, JacScaleFactor)

   ! Form G_hat^(-1) * (S^-1*dUdy) and G^(-1) * (S^-1*dUdu)
   ! factor G for the two solves:
   call LAPACK_getrf(M=size(G, 1), N=size(G, 2), A=G, IPIV=ipiv, ErrStat=ErrStat2, ErrMsg=ErrMsg2); if (Failed()) return

   ! after the this solve, dUdy holds G_hat^(-1) * dUdy_hat:
   call LAPACK_getrs(trans='N', N=size(G, 2), A=G, IPIV=ipiv, B=ModGlue%Lin%dUdy, ErrStat=ErrStat2, ErrMsg=ErrMsg2); if (Failed()) return

   ! after the this solve, dUdu holds G_hat^(-1) * dUdu_hat:
   call LAPACK_getrs(trans='N', N=size(G, 2), A=G, IPIV=ipiv, B=ModGlue%Lin%dUdu, ErrStat=ErrStat2, ErrMsg=ErrMsg2); if (Failed()) return

   ! Deallocate G and ipiv because the solves are complete
   deallocate (G)
   deallocate (ipiv)

   ! After this call, dUdu holds G^(-1)*dUdu and dUdy holds G^(-1)*dUdy
   call Postcondition(ModGlue%Vars%u, ModGlue%Lin%dUdu, ModGlue%Lin%dUdy, JacScaleFactor)

   ! Allocate tmp matrix for A and C calculations
   call AllocAry(tmp, ModGlue%Vars%Nu, ModGlue%Vars%Nx, 'G^-1*dUdy*C', ErrStat2, ErrMsg2); if (Failed()) return

   ! tmp = G^(-1) * dUdy * diag(C)
   call LAPACK_GEMM('N', 'N', 1.0_R8Ki, ModGlue%Lin%dUdy, ModGlue%Lin%dYdx, 0.0_R8Ki, tmp, ErrStat2, ErrMsg2); if (Failed()) return

   ! A
   ! dXdx = dXdx - matmul(dXdu, tmp)
   call LAPACK_GEMM('N', 'N', -1.0_R8Ki, ModGlue%Lin%dXdu, tmp, 1.0_R8Ki, ModGlue%Lin%dXdx, ErrStat2, ErrMsg2); if (Failed()) return

   ! C
   ! dYdx = dYdx - matmul(dYdu, tmp)
   call LAPACK_GEMM('N', 'N', -1.0_R8Ki, ModGlue%Lin%dYdu, tmp, 1.0_R8Ki, ModGlue%Lin%dYdx, ErrStat2, ErrMsg2); if (Failed()) return

   ! B
   tmp = ModGlue%Lin%dXdu
   ! dXdu = matmul(dXdu, dUdu)
   call LAPACK_GEMM('N', 'N', 1.0_R8Ki, tmp, ModGlue%Lin%dUdu, 0.0_R8Ki, ModGlue%Lin%dXdu, ErrStat2, ErrMsg2); if (Failed()) return

   ! D
   tmp = ModGlue%Lin%dYdu
   ! D = matmul(dYdu, dUdu)
   call LAPACK_GEMM('N', 'N', 1.0_R8Ki, tmp, ModGlue%Lin%dUdu, 0.0_R8Ki, ModGlue%Lin%dYdu, ErrStat2, ErrMsg2); if (Failed()) return

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
   integer(IntKi)                :: LoadFlags
   integer(IntKi)                :: i, j, k
   logical                       :: isRowLoad, isColLoad
   logical, allocatable          :: isLoad(:)

   allocate(isLoad(size(dUdu,1)))
   isLoad=.false.

   ! Loop through glue code input variables (cols)
   do i = 1, size(uVars)

      ! Get if col variable is a load
      isColLoad = uVars(i)%Field == VF_Force .or. uVars(i)%Field == VF_Moment

      ! Get col variable start and end indices in matrix
      associate (iLoc => uVars(i)%iLoc)

         isLoad(iLoc(1):iLoc(2)) = isColLoad

         ! Loop through glue code input variables (rows)
         do j = 1, size(uVars)

            ! Get if row variable is a load
            isRowLoad = uVars(j)%Field == VF_Force .or. uVars(j)%Field == VF_Moment

            ! Get row variable start and end indices in matrix
            associate (jLoc => uVars(j)%iLoc)

               if (isColLoad .and. (.not. isRowLoad)) then

                  ! Multiply columns of G
                  G(jLoc(1):jLoc(2), iLoc(1):iLoc(2)) = G(jLoc(1):jLoc(2), iLoc(1):iLoc(2))*JacScaleFactor

               else if (isRowLoad .and. (.not. isColLoad)) then

                  ! Divide rows of G
                  G(jLoc(1):jLoc(2), iLoc(1):iLoc(2)) = G(jLoc(1):jLoc(2), iLoc(1):iLoc(2))/JacScaleFactor

               end if

            end associate

         end do

         ! Divide rows of dUdu and dUdy by scale factor
         if (isColLoad) then
            dUdu(iLoc(1):iLoc(2), :) = dUdu(iLoc(1):iLoc(2), :)/JacScaleFactor
            dUdy(iLoc(1):iLoc(2), :) = dUdy(iLoc(1):iLoc(2), :)/JacScaleFactor
         end if

      end associate

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

      ! If variable is a (force or moment), apply post-conditioner
      if (uVars(i)%Field == VF_Force .or. uVars(i)%Field == VF_Moment) then

         ! Otherwise get variable start and end indices in matrix
         associate (iLoc => uVars(i)%iLoc)

            ! Multiply rows of dUdu
            dUdu(iLoc(1):iLoc(2), :) = dUdu(iLoc(1):iLoc(2), :)*JacScaleFactor

            ! Multiply rows of dUdy
            dUdy(iLoc(1):iLoc(2), :) = dUdy(iLoc(1):iLoc(2), :)*JacScaleFactor

         end associate

      end if
   end do

end subroutine

subroutine CalcWriteLinearMatrices(ModData, p_FAST, y_FAST, t_global, Un, OutFileName, IsGlue, ErrStat, ErrMsg)

   type(ModDataType), intent(inout) :: ModData        !< Module data
   type(FAST_ParameterType)         :: p_FAST         !< Parameters
   type(FAST_OutputFileType)        :: y_FAST         !< Output variables
   real(DbKi), intent(in)           :: t_global       !< current time step (written in file)
   integer(IntKi), intent(out)      :: Un             !< Unit number for file
   character(*), intent(in)         :: OutFileName    !< output file name
   logical                          :: IsGlue         !< Flag indicating this is writing glue code matrices
   integer(IntKi), intent(out)      :: ErrStat        !< Error status of the operation
   character(*), intent(out)        :: ErrMsg         !< Error message if ErrStat /= ErrID_None

   character(*), parameter          :: RoutineName = 'WriteModuleLinearMatrices'
   integer(IntKi)                   :: ErrStat2
   character(ErrMsgLen)             :: ErrMsg2
   character(32)                    :: Desc
   integer(IntKi)                   :: i
   integer(IntKi)                   :: Nx, Nxd, Nz, Nu, Ny
   character(50)                    :: Fmt
   logical, allocatable             :: uUse(:), yUse(:)

   ErrStat = ErrID_None
   ErrMsg = ""

   ! Open linearization file
   call OpenFOutFile(Un, OutFileName, ErrStat2, ErrMsg2); if (Failed()) return

   ! Calculate number of values in variable after applying filter
   Nx = MV_NumVars(ModData%Vars%x, VF_Linearize)
   Nxd = MV_NumVars(ModData%Vars%xd, VF_Linearize)
   Nz = MV_NumVars(ModData%Vars%z, VF_Linearize)
   Nu = MV_NumVars(ModData%Vars%u, VF_Linearize)
   Ny = MV_NumVars(ModData%Vars%y, VF_Linearize)

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
   Desc = 'Number of continuous states: '; write (Un, fmt) Desc, Nx
   Desc = 'Number of discrete states:   '; write (Un, fmt) Desc, Nxd
   Desc = 'Number of constraint states: '; write (Un, fmt) Desc, Nz
   Desc = 'Number of inputs:            '; write (Un, fmt) Desc, Nu
   Desc = 'Number of outputs:           '; write (Un, fmt) Desc, Ny

   Desc = 'Jacobians included in this file?'
   fmt = '(3x,A,1x,A5)'
   if (p_FAST%LinOutJac) then
      write (Un, fmt) Desc, 'Yes'
   else
      write (Un, fmt) Desc, 'No'
   end if

   write (Un, '()')    !print a blank line

   if (Nx > 0) then
      write (Un, '(A)') 'Order of continuous states:'
      call WrLinFile_txt_Table(ModData%Vars%x, VF_Linearize, p_FAST, Un, "Row/Column", ModData%Lin%x)

      write (Un, '(A)') 'Order of continuous state derivatives:'
      call WrLinFile_txt_Table(ModData%Vars%x, VF_Linearize, p_FAST, Un, "Row/Column", ModData%Lin%dx, IsDeriv=.true.)
   end if

   if (Nxd > 0) then
      write (Un, '(A)') 'Order of discrete states:'
      call WrLinFile_txt_Table(ModData%Vars%xd, VF_Linearize, p_FAST, Un, "Row/Column", ModData%Lin%xd)
   end if

   if (Nz > 0) then
      write (Un, '(A)') 'Order of constraint states:'
      call WrLinFile_txt_Table(ModData%Vars%z, VF_Linearize, p_FAST, Un, "Row/Column", ModData%Lin%z)
   end if

   if (Nu > 0) then
      write (Un, '(A)') 'Order of inputs:'
      call WrLinFile_txt_Table(ModData%Vars%u, VF_Linearize, p_FAST, Un, "Column  ", ModData%Lin%u, ShowRot=.true.)
   end if

   if (Ny > 0) then
      write (Un, '(A)') 'Order of outputs:'
      call WrLinFile_txt_Table(ModData%Vars%y, VF_Linearize, p_FAST, Un, "Row  ", ModData%Lin%y, ShowRot=.true.)
   end if

   ! Create boolean array indicating which input values to write
   allocate (uUse(ModData%Vars%Nu))
   uUse = .false.
   do i = 1, size(ModData%Vars%u)
      associate (Var => ModData%Vars%u(i))
         if (MV_HasFlags(Var, VF_Linearize)) uUse(Var%iLoc(1):Var%iLoc(2)) = .true.
      end associate
   end do

   ! Create boolean array indicating which output values to write
   allocate (yUse(ModData%Vars%Ny))
   yUse = .false.
   do i = 1, size(ModData%Vars%y)
      associate (Var => ModData%Vars%y(i))
         if (MV_HasFlags(Var, VF_Linearize)) yUse(Var%iLoc(1):Var%iLoc(2)) = .true.
      end associate
   end do

   ! If Jacobian matrix output is requested
   if (p_FAST%LinOutJac) then
      write (Un, '(/,A,/)') 'Jacobian matrices:'
      if (IsGlue) then
         call WrPartialMatrix(ModData%Lin%dUdu, Un, p_FAST%OutFmt, 'dUdu', UseRow=uUse, UseCol=uUse)
         call WrPartialMatrix(ModData%Lin%dUdy, Un, p_FAST%OutFmt, 'dUdy', UseRow=uUse, UseCol=yUse)
      else
         if (allocated(ModData%Lin%dXdx)) call WrPartialMatrix(ModData%Lin%dXdx, Un, p_FAST%OutFmt, 'dXdx')
         if (allocated(ModData%Lin%dXdu)) call WrPartialMatrix(ModData%Lin%dXdu, Un, p_FAST%OutFmt, 'dXdu', UseCol=uUse)
         if (allocated(ModData%Lin%dYdx)) call WrPartialMatrix(ModData%Lin%dYdx, Un, p_FAST%OutFmt, 'dYdx', UseRow=yUse)
         if (allocated(ModData%Lin%dYdu)) call WrPartialMatrix(ModData%Lin%dYdu, Un, p_FAST%OutFmt, 'dYdu', UseRow=yUse, UseCol=uUse)
      end if
   end if

   ! If this is glue code module, calculate the glue code state matrices (A, B, C, D)
   ! Called here, after writing dUdu and dUdy, because those matrices are overwritten
   ! in the process of calculating the other state matrices
   if (IsGlue) then
      call CalcGlueStateMatrices(ModData, real(p_FAST%UJacSclFact, R8Ki), ErrStat2, ErrMsg2)
      if (Failed()) return
   end if

   ! Write the linearized state matrices
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
end subroutine CalcWriteLinearMatrices

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

               ! Convert quaternion parameters to DCM
               DCM = quat_to_dcm(real(op(i_op:i_op + 2), R8Ki))

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

subroutine Idx_Init(Mods, ModOrder, Idx, FlagFilter, ErrStat, ErrMsg)
   type(ModDataType), intent(in)       :: Mods(:)
   integer(IntKi), intent(in)          :: ModOrder(:)
   type(VarsIdxType), intent(out)      :: Idx
   integer(IntKi), intent(in)          :: FlagFilter
   integer(IntKi), intent(out)         :: ErrStat
   character(ErrMsgLen), intent(out)   :: ErrMsg

   character(*), parameter             :: RoutineName = 'Idx_Init'
   integer(IntKi)                      :: ErrStat2
   character(ErrMsgLen)                :: ErrMsg2
   integer(IntKi)                      :: NumVars
   integer(IntKi)                      :: iGbl(2)
   integer(IntKi)                      :: i, j

   ! Initialize error return
   ErrStat = ErrID_None
   ErrMsg = ""

   ! Destroy VarIdx in case it has been previously used
   call Glue_DestroyVarsIdxType(Idx, ErrStat2, ErrMsg2); if (Failed()) return

   ! Save filter in index
   Idx%FlagFilter = FlagFilter

   !----------------------------------------------------------------------------
   ! Indexing Data Description
   !----------------------------------------------------------------------------

   ! For each variable (x, u, y, etc.) there are two arrays:
   !     1) Variable local and global value indices (ValLocGbl)
   !     2) Module variable start index (ModVarStart)
   ! ValLocGbl has 4 rows and N columns where N is the total number of variables
   ! for all modules in Mods. The columns are as follows:
   !     1) Values start index inside module arrays/matrices (iLoc(1))
   !     2) Values end index inside module arrays/matrices (iLoc(2))
   !     3) Values start index in global arrays/matrices (iGbl(1))
   !     4) Values end index in global arrays/matrices (iLoc(2))
   ! ModVarStart contains N rows where N is the total number of modules in Mods.
   ! The values in this array contain the variable start index offset for each
   ! module into ValLocGbl so value indices can be looked up given module index
   ! and variable index. Keeping all value indices in one matrix makes data
   ! storage much simpler at the cost of of having to maintain the array of
   ! module offsets.

   !----------------------------------------------------------------------------
   ! Build index for continuous state variables
   !----------------------------------------------------------------------------

   ! Allocate array of module variable start indices for each module, init to 0
   call AllocAry(Idx%x%ModVarStart, size(Mods) + 1, "VarIdx%x%ModVarStart", ErrStat2, ErrMsg2); if (Failed()) return
   Idx%x%ModVarStart(1) = 0

   ! Populate ModVarStart with variable offsets and calculate total number of variables
   NumVars = 0
   do i = 1, size(Mods)
      NumVars = NumVars + size(Mods(i)%Vars%x)
      Idx%x%ModVarStart(i + 1) = NumVars
   end do

   ! Allocate variable value index matrix and initialize to zero
   call AllocAry(Idx%x%ValLocGbl, 4, NumVars, "VarIdx%x%ValLocGbl", ErrStat2, ErrMsg2); if (Failed()) return
   Idx%x%ValLocGbl = 0

   ! Initialize global index to zero
   iGbl = 0

   ! Loop through modules and variables, add value indices to index if variable has filter flags, increment global indices
   do i = 1, size(ModOrder)
      associate (ModData => Mods(ModOrder(i)))
         do j = 1, size(ModData%Vars%x)
            if (MV_HasFlags(ModData%Vars%x(j), FlagFilter)) then
               iGbl(1) = iGbl(2) + 1
               iGbl(2) = iGbl(1) + ModData%Vars%x(j)%Num - 1
               Idx%x%ValLocGbl(:, Idx%x%ModVarStart(ModData%Idx) + j) = [ModData%Vars%x(j)%iLoc, iGbl]
            end if
         end do
      end associate
   end do

   ! Save total number of values
   Idx%Nx = iGbl(2)

   !----------------------------------------------------------------------------
   ! Build index for discrete state variables
   !----------------------------------------------------------------------------

   ! Allocate array of module variable start indices for each module, init to 0
   call AllocAry(Idx%xd%ModVarStart, size(Mods) + 1, "VarIdx%xd%ModVarStart", ErrStat2, ErrMsg2); if (Failed()) return
   Idx%xd%ModVarStart(1) = 0

   ! Populate ModVarStart with variable offsets and calculate total number of variables and values
   NumVars = 0
   do i = 1, size(Mods)
      NumVars = NumVars + size(Mods(i)%Vars%xd)
      Idx%xd%ModVarStart(i + 1) = NumVars
   end do

   ! Allocate variable value index matrix and initialize to zero
   call AllocAry(Idx%xd%ValLocGbl, 4, NumVars, "VarIdx%xd%ValLocGbl", ErrStat2, ErrMsg2); if (Failed()) return
   Idx%xd%ValLocGbl = 0

   ! Initialize global index and number of values to zero
   iGbl = 0

   ! Loop through modules and variables, add value indices to index if variable has filter flags, increment global indices
   do i = 1, size(ModOrder)
      associate (ModData => Mods(ModOrder(i)))
         do j = 1, size(ModData%Vars%xd)
            if (MV_HasFlags(ModData%Vars%xd(j), FlagFilter)) then
               iGbl(1) = iGbl(2) + 1
               iGbl(2) = iGbl(1) + ModData%Vars%xd(j)%Num - 1
               Idx%xd%ValLocGbl(:, Idx%xd%ModVarStart(ModData%Idx) + j) = [ModData%Vars%xd(j)%iLoc, iGbl]
            end if
         end do
      end associate
   end do

   ! Save total number of values
   Idx%Nxd = iGbl(2)

   !----------------------------------------------------------------------------
   ! Build index for constraint state variables
   !----------------------------------------------------------------------------

   ! Allocate array of module variable start indices for each module, init to 0
   call AllocAry(Idx%z%ModVarStart, size(Mods) + 1, "VarIdx%z%ModVarStart", ErrStat2, ErrMsg2); if (Failed()) return
   Idx%z%ModVarStart(1) = 0

   ! Populate ModVarStart with variable offsets and calculate total number of variables
   NumVars = 0
   do i = 1, size(Mods)
      NumVars = NumVars + size(Mods(i)%Vars%z)
      Idx%z%ModVarStart(i + 1) = NumVars
   end do

   ! Allocate variable value index matrix and initialize to zero
   call AllocAry(Idx%z%ValLocGbl, 4, NumVars, "VarIdx%z%ValLocGbl", ErrStat2, ErrMsg2); if (Failed()) return
   Idx%z%ValLocGbl = 0

   ! Initialize global index to zero
   iGbl = 0

   ! Loop through modules and variables, add value indices to index if variable has filter flags, increment global indices
   do i = 1, size(ModOrder)
      associate (ModData => Mods(ModOrder(i)))
         do j = 1, size(ModData%Vars%z)
            if (MV_HasFlags(ModData%Vars%z(j), FlagFilter)) then
               iGbl(1) = iGbl(2) + 1
               iGbl(2) = iGbl(1) + ModData%Vars%z(j)%Num - 1
               Idx%z%ValLocGbl(:, Idx%z%ModVarStart(ModData%Idx) + j) = [ModData%Vars%z(j)%iLoc, iGbl]
            end if
         end do
      end associate
   end do

   ! Save total number of values
   Idx%Nz = iGbl(2)

   !----------------------------------------------------------------------------
   ! Build index for input variables
   !----------------------------------------------------------------------------

   ! Allocate array of module variable start indices for each module, init to 0
   call AllocAry(Idx%u%ModVarStart, size(Mods) + 1, "VarIdx%u%ModVarStart", ErrStat2, ErrMsg2); if (Failed()) return
   Idx%u%ModVarStart(1) = 0

   ! Populate ModVarStart with variable offsets and calculate total number of variables
   NumVars = 0
   do i = 1, size(Mods)
      NumVars = NumVars + size(Mods(i)%Vars%u)
      Idx%u%ModVarStart(i + 1) = NumVars
   end do

   ! Allocate variable value index matrix and initialize to zero
   call AllocAry(Idx%u%ValLocGbl, 4, NumVars, "VarIdx%u%ValLocGbl", ErrStat2, ErrMsg2); if (Failed()) return
   Idx%u%ValLocGbl = 0

   ! Initialize global index to zero
   iGbl = 0

   ! Loop through modules and variables, add value indices to index if variable has filter flags, increment global indices
   do i = 1, size(ModOrder)
      associate (ModData => Mods(ModOrder(i)))
         do j = 1, size(ModData%Vars%u)
            if (MV_HasFlags(ModData%Vars%u(j), FlagFilter)) then
               iGbl(1) = iGbl(2) + 1
               iGbl(2) = iGbl(1) + ModData%Vars%u(j)%Num - 1
               Idx%u%ValLocGbl(:, Idx%u%ModVarStart(ModData%Idx) + j) = [ModData%Vars%u(j)%iLoc, iGbl]
            end if
         end do
      end associate
   end do

   ! Save total number of values
   Idx%Nu = iGbl(2)

   !----------------------------------------------------------------------------
   ! Build index for output variables
   !----------------------------------------------------------------------------

   ! Allocate array of module variable start indices for each module, init to 0
   call AllocAry(Idx%y%ModVarStart, size(Mods) + 1, "VarIdx%y%ModVarStart", ErrStat2, ErrMsg2); if (Failed()) return
   Idx%y%ModVarStart(1) = 0

   ! Populate ModVarStart with variable offsets and calculate total number of variables
   NumVars = 0
   do i = 1, size(Mods)
      NumVars = NumVars + size(Mods(i)%Vars%y)
      Idx%y%ModVarStart(i + 1) = NumVars
   end do

   ! Allocate variable value index matrix and initialize to zero
   call AllocAry(Idx%y%ValLocGbl, 4, NumVars, "VarIdx%y%ValLocGbl", ErrStat2, ErrMsg2); if (Failed()) return
   Idx%y%ValLocGbl = 0

   ! Initialize global index to zero
   iGbl = 0

   ! Loop through modules and variables, add value indices to index if variable has filter flags, increment global indices
   do i = 1, size(ModOrder)
      associate (ModData => Mods(ModOrder(i)))
         do j = 1, size(ModData%Vars%y)
            if (MV_HasFlags(ModData%Vars%y(j), FlagFilter)) then
               iGbl(1) = iGbl(2) + 1
               iGbl(2) = iGbl(1) + ModData%Vars%y(j)%Num - 1
               Idx%y%ValLocGbl(:, Idx%y%ModVarStart(ModData%Idx) + j) = [ModData%Vars%y(j)%iLoc, iGbl]
            end if
         end do
      end associate
   end do

   ! Save total number of values
   Idx%Ny = iGbl(2)

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

logical function Idx_GetLocGbl(Idx, ModIdx, VarIdx, iLoc, iGbl)
   type(VarIdxType), intent(in)  :: Idx
   integer(IntKi), intent(in)    :: ModIdx, VarIdx
   integer(IntKi), intent(out)   :: iLoc(2), iGbl(2)
   integer(IntKi)                :: iLocGbl(4)
   iLocGbl = Idx%ValLocGbl(:, Idx%ModVarStart(ModIdx) + VarIdx)
   iLoc = iLocGbl(1:2)
   iGbl = iLocGbl(3:4)
   Idx_GetLocGbl = iLocGbl(3) /= 0  ! Variable has global index
end function

subroutine MV_AddModule(Mods, ModID, ModAbbr, Instance, ModDT, SolverDT, Vars, ErrStat, ErrMsg)
   type(ModDataType), allocatable, intent(inout)   :: Mods(:)
   integer(IntKi), intent(in)                      :: ModID
   character(*), intent(in)                        :: ModAbbr
   integer(IntKi), intent(in)                      :: Instance
   real(R8Ki), intent(in)                          :: ModDT
   real(R8Ki), intent(in)                          :: SolverDT
   type(ModVarsType), pointer, intent(in)          :: Vars
   integer(IntKi), intent(out)                     :: ErrStat
   character(ErrMsgLen), intent(out)               :: ErrMsg

   character(*), parameter                         :: RoutineName = 'MV_AddModule'
   integer(IntKi)                                  :: ErrStat2
   character(ErrMsgLen)                            :: ErrMsg2
   type(ModDataType)                               :: ModData

   ErrStat = ErrID_None
   ErrMsg = ''

   ! If module array hasn't been allocated, allocate with zero size
   if (.not. allocated(Mods)) allocate (Mods(0))

   ! Populate ModuleDataType derived type
   ModData = ModDataType(Idx=size(Mods) + 1, ID=ModID, Abbr=ModAbbr, &
                         Ins=Instance, DT=ModDT, Vars=Vars)

   ! Allocate source and destination mapping arrays
   call AllocAry(ModData%SrcMaps, 0, "ModData%SrcMaps", ErrStat2, ErrMsg2)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return
   call AllocAry(ModData%DstMaps, 0, "ModData%DstMaps", ErrStat2, ErrMsg2)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return

   !----------------------------------------------------------------------------
   ! Calculate Module Substepping
   !----------------------------------------------------------------------------

   ! If module time step is same as global time step, set substeps to 1
   if (EqualRealNos(ModData%DT, SolverDT)) then
      ModData%SubSteps = 1
   else
      ! If the module time step is greater than the global time step, set error
      if (ModData%DT > SolverDT) then
         call SetErrStat(ErrID_Fatal, "The "//trim(ModData%Abbr)// &
                         " module time step ("//trim(Num2LStr(ModData%DT))//" s) "// &
                         "cannot be larger than FAST time step ("//trim(Num2LStr(SolverDT))//" s).", &
                         ErrStat, ErrMsg, RoutineName)
         return
      end if

      ! Calculate the number of substeps
      ModData%SubSteps = nint(SolverDT/ModData%DT)

      ! If the module DT is not an exact integer divisor of the global time step, set error
      if (.not. EqualRealNos(SolverDT, ModData%DT*ModData%SubSteps)) then
         call SetErrStat(ErrID_Fatal, "The "//trim(ModData%Abbr)// &
                         " module time step ("//trim(Num2LStr(ModData%DT))//" s) "// &
                         "must be an integer divisor of the FAST time step ("//trim(Num2LStr(SolverDT))//" s).", &
                         ErrStat, ErrMsg, RoutineName)
         return
      end if
   end if

   !----------------------------------------------------------------------------
   ! Add module data to array
   !----------------------------------------------------------------------------

   Mods = [Mods, ModData]

end subroutine

end module
