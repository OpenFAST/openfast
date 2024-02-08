!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2023  National Renewable Energy Laboratory
!
!    This file is part of the NWTC Subroutine Library.
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
!> The modules ModVar and ModVar_Types provide data structures and subroutines for representing and manipulating meshes
!! and meshed data in the FAST modular framework.
!!
!! Module variables provide a structured way for documenting, locating, and orchestrating the interdependencies between modules.
!!

module ModVar
use NWTC_Library_Types
use NWTC_IO
use ModMesh
implicit none

private
public :: MV_InitVarsJac, MV_Pack, MV_Unpack
public :: MV_ComputeCentralDiff, MV_Perturb, MV_ComputeDiff
public :: MV_AddVar, MV_AddMeshVar, MV_AddModule
public :: MV_HasFlags, MV_SetFlags, MV_UnsetFlags, MV_NumVars
public :: LoadFields, MotionFields, TransFields, AngularFields
public :: wm_to_dcm, wm_compose, wm_from_dcm, wm_inv, wm_to_rvec, wm_from_rvec
public :: MV_FieldString, IdxStr
public :: MV_InitVarIdx

integer(IntKi), parameter :: &
   LoadFields(*) = [VF_Force, VF_Moment], &
   TransFields(*) = [VF_TransDisp, VF_TransVel, VF_TransAcc], &
   AngularFields(*) = [VF_Orientation, VF_AngularDisp, VF_AngularVel, VF_AngularAcc], &
   MotionFields(*) = [VF_TransDisp, VF_Orientation, VF_TransVel, VF_AngularVel, VF_TransAcc, VF_AngularAcc]

interface MV_Pack
   module procedure MV_PackVarRank0R4, MV_PackVarRank1R4, MV_PackVarRank2R4
   module procedure MV_PackVarRank0R8, MV_PackVarRank1R8, MV_PackVarRank2R8
   module procedure MV_PackMesh
end interface

interface MV_Unpack
   module procedure MV_UnpackVarRank0R4, MV_UnpackVarRank1R4, MV_UnpackVarRank2R4
   module procedure MV_UnpackVarRank0R8, MV_UnpackVarRank1R8, MV_UnpackVarRank2R8
   module procedure MV_UnpackMesh
end interface

contains

function MV_FieldString(Field) result(str)
   integer(IntKi), intent(in) :: Field
   character(16)              :: str
   select case (Field)
   case (VF_AngularAcc)
      str = "VF_AngularAcc"
   case (VF_AngularDisp)
      str = "VF_AngularDisp"
   case (VF_AngularVel)
      str = "VF_AngularVel"
   case (VF_Force)
      str = "VF_Force"
   case (VF_Moment)
      str = "VF_Moment"
   case (VF_Orientation)
      str = "VF_Orientation"
   case (VF_TransAcc)
      str = "VF_TransAcc"
   case (VF_TransDisp)
      str = "VF_TransDisp"
   case (VF_TransVel)
      str = "VF_TransVel"
   case default
      str = "Unknown"
   end select
end function

subroutine MV_InitVarsJac(Vars, Jac, Linearize, ErrStat, ErrMsg)
   type(ModVarsType), pointer, intent(inout)    :: Vars
   type(ModJacType), intent(inout)              :: Jac
   logical, intent(in)                          :: Linearize
   integer(IntKi), intent(out)                  :: ErrStat
   character(ErrMsgLen), intent(out)            :: ErrMsg

   character(*), parameter       :: RoutineName = 'MV_InitVarsLin'
   integer(IntKi)                :: ErrStat2
   character(ErrMsgLen)          :: ErrMsg2
   integer(IntKi)                :: i, StartIndex
   type(VarsIdxType), pointer    :: VarIdx

   ! Initialize error outputs
   ErrStat = ErrID_None
   ErrMsg = ''

   ! Initialize number of variables in each group
   Vars%Nx = 0
   Vars%Nxd = 0
   Vars%Nz = 0
   Vars%Nu = 0
   Vars%Ny = 0

   ! Initialize continuous state variables
   if (.not. allocated(Vars%x)) allocate (Vars%x(0))
   StartIndex = 1
   do i = 1, size(Vars%x)
      call ModVarType_Init(Vars%x(i), StartIndex, Linearize, ErrStat2, ErrMsg2)
      if (Failed()) return
   end do
   Vars%Nx = sum(Vars%x%Num)

   ! Initialize discrete state variables
   if (.not. allocated(Vars%xd)) allocate (Vars%xd(0))
   StartIndex = 1
   do i = 1, size(Vars%xd)
      call ModVarType_Init(Vars%xd(i), StartIndex, Linearize, ErrStat2, ErrMsg2)
      if (Failed()) return
   end do
   Vars%Nxd = sum(Vars%xd%Num)

   ! Initialize constraint state variables
   if (.not. allocated(Vars%z)) allocate (Vars%z(0))
   StartIndex = 1
   do i = 1, size(Vars%z)
      call ModVarType_Init(Vars%z(i), StartIndex, Linearize, ErrStat2, ErrMsg2)
      if (Failed()) return
   end do
   Vars%Nz = sum(Vars%z%Num)

   ! Initialize input variables
   if (.not. allocated(Vars%u)) allocate (Vars%u(0))
   StartIndex = 1
   do i = 1, size(Vars%u)
      call ModVarType_Init(Vars%u(i), StartIndex, Linearize, ErrStat2, ErrMsg2)
      if (Failed()) return
   end do
   Vars%Nu = sum(Vars%u%Num)

   ! Initialize output variables
   if (.not. allocated(Vars%y)) allocate (Vars%y(0))
   StartIndex = 1
   do i = 1, size(Vars%y)
      call ModVarType_Init(Vars%y(i), StartIndex, Linearize, ErrStat2, ErrMsg2)
      if (Failed()) return
   end do
   Vars%Ny = sum(Vars%y%Num)

   ! Initialize full linearization variable indexing (all variables)
   call MV_InitVarIdx(Vars, Vars%IdxLin, VF_None, ErrStat2, ErrMsg2); if (Failed()) return

   ! Allocate arrays
   if (Vars%Nx > 0) then
      call AllocAry(Jac%x, Vars%Nx, "Lin%x", ErrStat2, ErrMsg2); if (Failed()) return
      call AllocAry(Jac%dx, Vars%Nx, "Lin%dx", ErrStat2, ErrMsg2); if (Failed()) return
      call AllocAry(Jac%x_perturb, Vars%Nx, "Lin%x_perturb", ErrStat2, ErrMsg2); if (Failed()) return
      call AllocAry(Jac%x_pos, Vars%Nx, "Lin%x_pos", ErrStat2, ErrMsg2); if (Failed()) return
      call AllocAry(Jac%x_neg, Vars%Nx, "Lin%x_neg", ErrStat2, ErrMsg2); if (Failed()) return
   end if
   if (Vars%Nxd > 0) then
      call AllocAry(Jac%xd, Vars%Nxd, "Lin%xd", ErrStat2, ErrMsg2); if (Failed()) return
   end if
   if (Vars%Nz > 0) then
      call AllocAry(Jac%z, Vars%Nz, "Lin%z", ErrStat2, ErrMsg2); if (Failed()) return
   end if
   if (Vars%Nu > 0) then
      call AllocAry(Jac%u, Vars%Nu, "Lin%u", ErrStat2, ErrMsg2); if (Failed()) return
      call AllocAry(Jac%u_perturb, Vars%Nu, "Lin%u_perturb", ErrStat2, ErrMsg2); if (Failed()) return
   end if
   if (Vars%Ny > 0) then
      call AllocAry(Jac%y, Vars%Ny, "Lin%y", ErrStat2, ErrMsg2); if (Failed()) return
      call AllocAry(Jac%y_pos, Vars%Ny, "Lin%y_pos", ErrStat2, ErrMsg2); if (Failed()) return
      call AllocAry(Jac%y_neg, Vars%Ny, "Lin%y_neg", ErrStat2, ErrMsg2); if (Failed()) return
   end if

contains

   function Failed()
      logical Failed
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function

   function FailedAlloc()
      logical FailedAlloc
      FailedAlloc = ErrStat2 /= 0
      if (FailedAlloc) call SetErrStat(ErrID_Fatal, 'error allocating Vals', ErrStat, ErrMsg, RoutineName)
   end function

end subroutine

elemental function IsMesh(Var) result(r)
   type(ModVarType), intent(in)  :: Var
   logical                 :: r
   r = iand(Var%Flags, VF_Mesh) > 0
end function

subroutine ModVarType_Init(Var, Index, Linearize, ErrStat, ErrMsg)
   type(ModVarType), intent(inout)     :: Var
   integer(IntKi), intent(inout)       :: Index
   logical, intent(in)                 :: Linearize
   integer(IntKi), intent(out)         :: ErrStat
   character(ErrMsgLen), intent(out)   :: ErrMsg

   character(*), parameter       :: RoutineName = 'ModVarsType_Init'
   integer(IntKi)                :: ErrStat2
   character(ErrMsgLen)          :: ErrMsg2
   integer(IntKi)                :: i, j
   character(1), parameter       :: Comp(3) = ['X', 'Y', 'Z']
   character(*), parameter       :: Fmt = '(A," ",A,", node",I0,", ",A)'
   character(2)                  :: UnitDesc

   ! Initialize error outputs
   ErrStat = ErrID_None
   ErrMsg = ''

   !----------------------------------------------------------------------------
   ! Mesh
   !----------------------------------------------------------------------------

   ! If this variable belongs to a mesh
   if (MV_HasFlags(Var, VF_Mesh)) then

      ! Size is the number of nodes in a mesh
      Var%Nodes = Var%Num

      ! Number of values
      Var%Num = Var%Nodes*3

      ! If linearization enabled
      if (Linearize) then

         ! Set unit description for line mesh
         UnitDesc = ''
         if (MV_HasFlags(Var, VF_Line)) UnitDesc = "/m"

         ! Switch based on field number
         select case (Var%Field)
         case (VF_Force)
            Var%LinNames = [character(LinChanLen) ::((trim(Var%Name)//" "//Comp(j)//" force, node "//trim(num2lstr(i))//', N'//UnitDesc, j=1, 3), i=1, Var%Nodes)]
         case (VF_Moment)
            Var%LinNames = [character(LinChanLen) ::((trim(Var%Name)//" "//Comp(j)//" moment, node "//trim(num2lstr(i))//', Nm'//UnitDesc, j=1, 3), i=1, Var%Nodes)]
         case (VF_TransDisp)
            Var%LinNames = [character(LinChanLen) ::((trim(Var%Name)//" "//Comp(j)//" translation displacement, node "//trim(num2lstr(i))//', m', j=1, 3), i=1, Var%Nodes)]
         case (VF_Orientation)
            Var%LinNames = [character(LinChanLen) ::((trim(Var%Name)//" "//Comp(j)//" orientation angle, node "//trim(num2lstr(i))//', rad', j=1, 3), i=1, Var%Nodes)]
         case (VF_TransVel)
            Var%LinNames = [character(LinChanLen) ::((trim(Var%Name)//" "//Comp(j)//" translation velocity, node "//trim(num2lstr(i))//', m/s', j=1, 3), i=1, Var%Nodes)]
         case (VF_AngularVel)
            Var%LinNames = [character(LinChanLen) ::((trim(Var%Name)//" "//Comp(j)//" rotation velocity, node "//trim(num2lstr(i))//', rad/s', j=1, 3), i=1, Var%Nodes)]
         case (VF_TransAcc)
            Var%LinNames = [character(LinChanLen) ::((trim(Var%Name)//" "//Comp(j)//" translation acceleration, node "//trim(num2lstr(i))//', m/s^2', j=1, 3), i=1, Var%Nodes)]
         case (VF_AngularAcc)
            Var%LinNames = [character(LinChanLen) ::((trim(Var%Name)//" "//Comp(j)//" rotation acceleration, node "//trim(num2lstr(i))//', rad/s^2', j=1, 3), i=1, Var%Nodes)]
         case default
            call SetErrStat(ErrID_Fatal, "Invalid mesh field type", ErrStat, ErrMsg, RoutineName)
            return
         end select

      end if
   end if

   !----------------------------------------------------------------------------
   ! Linearization
   !----------------------------------------------------------------------------

   if (Linearize) then
      if (.not. allocated(Var%LinNames)) then
         call SetErrStat(ErrID_Fatal, "LinNames not allocated for "//Var%Name, ErrStat, ErrMsg, RoutineName)
         return
      else if (size(Var%LinNames) < Var%Num) then
         call SetErrStat(ErrID_Fatal, "insufficient LinNames given for "//Var%Name, ErrStat, ErrMsg, RoutineName)
         return
      else if (size(Var%LinNames) > Var%Num) then
         call SetErrStat(ErrID_Fatal, "excessive LinNames given for "//Var%Name, ErrStat, ErrMsg, RoutineName)
         return
      end if
   else
      ! Deallocate linearization names if linearization is not enabled
      if (allocated(Var%LinNames)) deallocate (Var%LinNames)
   end if

   !----------------------------------------------------------------------------
   ! Indices
   !----------------------------------------------------------------------------

   ! Set start and end indices for local matrices
   Var%iLoc = [index, index + Var%Num - 1]

   ! Update index based on variable size
   index = index + Var%Num

contains
   function Failed()
      logical :: Failed
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

!-------------------------------------------------------------------------------
! Functions for packing and unpacking data by variable
!-------------------------------------------------------------------------------

subroutine MV_PackMatrix(RowVarAry, ColVarAry, FlagFilter, M, SubM)
   type(ModVarType), intent(in)  :: RowVarAry(:), ColVarAry(:)
   real(R8Ki), intent(in)        :: M(:, :)
   real(R8Ki), intent(inout)     :: SubM(:, :)
   integer(IntKi), intent(in)    :: FlagFilter
   integer(IntKi)                :: i, j
   integer(IntKi)                :: row, col
   col = 1
   row = 1
   do i = 1, size(ColVarAry)
      if (iand(ColVarAry(i)%Flags, FlagFilter) == 0) cycle
      do j = 1, size(RowVarAry)
         if (iand(RowVarAry(j)%Flags, FlagFilter) == 0) cycle
         associate (rVar => RowVarAry(i), cVar => ColVarAry(i))
            SubM(row:row + rVar%Num - 1, col:col + cVar%Num - 1) = M(rVar%iLoc(1):rVar%iLoc(2), cVar%iLoc(1):cVar%iLoc(2))
         end associate
         row = row + RowVarAry(j)%Num - 1
      end do
      col = col + ColVarAry(i)%Num - 1
   end do
end subroutine

subroutine MV_PackVarRank0R4(VarAry, iVar, Val, Ary)
   type(ModVarType), intent(in)  :: VarAry(:)
   integer(IntKi), intent(in)    :: iVar
   real(R4Ki), intent(in)        :: Val
   real(R8Ki), intent(inout)     :: Ary(:)
   if (iVar == 0) return
   Ary(VarAry(iVar)%iLoc(1)) = real(Val, R8Ki)
end subroutine

subroutine MV_PackVarRank0R8(VarAry, iVar, Val, Ary)
   type(ModVarType), intent(in)  :: VarAry(:)
   integer(IntKi), intent(in)    :: iVar
   real(R8Ki), intent(in)        :: Val
   real(R8Ki), intent(inout)     :: Ary(:)
   if (iVar == 0) return
   Ary(VarAry(iVar)%iLoc(1)) = Val
end subroutine

subroutine MV_PackVarRank1R4(VarAry, iVar, Vals, Ary)
   type(ModVarType), intent(in)  :: VarAry(:)
   integer(IntKi), intent(in)    :: iVar
   real(R4Ki), intent(in)        :: Vals(:)
   real(R8Ki), intent(inout)     :: Ary(:)
   associate (iLoc => VarAry(iVar)%iLoc)
      Ary(iLoc(1):iLoc(2)) = real(Vals, R8Ki)
   end associate
end subroutine

subroutine MV_PackVarRank1R8(VarAry, iVar, Vals, Ary)
   type(ModVarType), intent(in)  :: VarAry(:)
   integer(IntKi), intent(in)    :: iVar
   real(R8Ki), intent(in)        :: Vals(:)
   real(R8Ki), intent(inout)     :: Ary(:)
   if (iVar == 0) return
   associate (iLoc => VarAry(iVar)%iLoc)
      Ary(iLoc(1):iLoc(2)) = Vals
   end associate
end subroutine

subroutine MV_PackVarRank2R4(VarAry, iVar, Vals, Ary)
   type(ModVarType), intent(in)  :: VarAry(:)
   integer(IntKi), intent(in)    :: iVar
   real(R4Ki), intent(in)        :: Vals(:, :)
   real(R8Ki), intent(inout)     :: Ary(:)
   if (iVar == 0) return
   associate (iLoc => VarAry(iVar)%iLoc)
      Ary(iLoc(1):iLoc(2)) = pack(real(Vals, R8Ki), .true.)
   end associate
end subroutine

subroutine MV_PackVarRank2R8(VarAry, iVar, Vals, Ary)
   type(ModVarType), intent(in)  :: VarAry(:)
   integer(IntKi), intent(in)    :: iVar
   real(R8Ki), intent(in)        :: Vals(:, :)
   real(R8Ki), intent(inout)     :: Ary(:)
   if (iVar == 0) return
   associate (iLoc => VarAry(iVar)%iLoc)
      Ary(iLoc(1):iLoc(2)) = pack(Vals, .true.)
   end associate
end subroutine

subroutine MV_UnpackVarRank0R4(VarAry, iVar, Ary, Val)
   type(ModVarType), intent(in)  :: VarAry(:)
   integer(IntKi), intent(in)    :: iVar
   real(R8Ki), intent(in)        :: Ary(:)
   real(R4Ki), intent(inout)     :: Val
   Val = Ary(VarAry(iVar)%iLoc(1))
end subroutine

subroutine MV_UnpackVarRank0R8(VarAry, iVar, Ary, Vals)
   type(ModVarType), intent(in)  :: VarAry(:)
   integer(IntKi), intent(in)    :: iVar
   real(R8Ki), intent(in)        :: Ary(:)
   real(R8Ki), intent(inout)     :: Vals
   Vals = Ary(VarAry(iVar)%iLoc(1))
end subroutine

subroutine MV_UnpackVarRank1R4(VarAry, iVar, Ary, Vals)
   type(ModVarType), intent(in)  :: VarAry(:)
   integer(IntKi), intent(in)    :: iVar
   real(R8Ki), intent(in)        :: Ary(:)
   real(R4Ki), intent(inout)     :: Vals(:)
   if (iVar == 0) return
   associate (iLoc => VarAry(iVar)%iLoc)
      Vals = real(Ary(iLoc(1):iLoc(2)), R4Ki)
   end associate
end subroutine

subroutine MV_UnpackVarRank1R8(VarAry, iVar, Ary, Vals)
   type(ModVarType), intent(in)  :: VarAry(:)
   integer(IntKi), intent(in)    :: iVar
   real(R8Ki), intent(in)        :: Ary(:)
   real(R8Ki), intent(inout)     :: Vals(:)
   if (iVar == 0) return
   associate (iLoc => VarAry(iVar)%iLoc)
      Vals = Ary(iLoc(1):iLoc(2))
   end associate
end subroutine

subroutine MV_UnpackVarRank2R4(VarAry, iVar, Ary, Vals)
   type(ModVarType), intent(in)  :: VarAry(:)
   integer(IntKi), intent(in)    :: iVar
   real(R8Ki), intent(in)        :: Ary(:)
   real(R4Ki), intent(inout)     :: Vals(:, :)
   if (iVar == 0) return
   associate (iLoc => VarAry(iVar)%iLoc)
      Vals = reshape(real(Ary(iLoc(1):iLoc(2)), R4Ki), shape(Vals))
   end associate
end subroutine

subroutine MV_UnpackVarRank2R8(VarAry, iVar, Ary, Vals)
   type(ModVarType), intent(in)  :: VarAry(:)
   integer(IntKi), intent(in)    :: iVar
   real(R8Ki), intent(in)        :: Ary(:)
   real(R8Ki), intent(inout)     :: Vals(:, :)
   if (iVar == 0) return
   associate (iLoc => VarAry(iVar)%iLoc)
      Vals = reshape(Ary(iLoc(1):iLoc(2)), shape(Vals))
   end associate
end subroutine

subroutine MV_PackMesh(VarAry, iVar, Mesh, Values)
   type(ModVarType), intent(in)  :: VarAry(:)
   integer(IntKi), intent(in)    :: iVar
   type(MeshType), intent(in)    :: Mesh
   real(R8Ki), intent(inout)     :: Values(:)
   integer(IntKi)                :: MeshID, i, j
   if (iVar == 0) return
   MeshID = VarAry(iVar)%MeshID
   do i = iVar, size(VarAry)
      if (VarAry(i)%MeshID /= MeshID) exit
      associate (iLoc => VarAry(i)%iLoc)
         select case (VarAry(i)%Field)
         case (VF_Force)
            Values(iLoc(1):iLoc(2)) = pack(Mesh%Force, .true.)
         case (VF_Moment)
            Values(iLoc(1):iLoc(2)) = pack(Mesh%Moment, .true.)
         case (VF_TransDisp)
            Values(iLoc(1):iLoc(2)) = pack(Mesh%TranslationDisp, .true.)
         case (VF_Orientation)
            do j = 1, VarAry(i)%Nodes
               Values(iLoc(1) + 3*(j - 1):iLoc(1) + 3*j) = wm_from_dcm(Mesh%Orientation(:, :, j))
            end do
         case (VF_TransVel)
            Values(iLoc(1):iLoc(2)) = pack(Mesh%TranslationVel, .true.)
         case (VF_AngularVel)
            Values(iLoc(1):iLoc(2)) = pack(Mesh%RotationVel, .true.)
         case (VF_TransAcc)
            Values(iLoc(1):iLoc(2)) = pack(Mesh%TranslationAcc, .true.)
         case (VF_AngularAcc)
            Values(iLoc(1):iLoc(2)) = pack(Mesh%RotationAcc, .true.)
         case (VF_Scalar)
            Values(iLoc(1):iLoc(2)) = pack(Mesh%Scalars, .true.)
         end select
      end associate
   end do
end subroutine

subroutine MV_UnpackMesh(VarAry, iVar, Values, Mesh)
   type(ModVarType), intent(in)  :: VarAry(:)
   integer(IntKi), intent(in) :: iVar
   real(R8Ki), intent(in)        :: Values(:)
   type(MeshType), intent(inout) :: Mesh
   integer(IntKi)                :: MeshID, i, j
   MeshID = VarAry(iVar)%MeshID
   do i = iVar, size(VarAry)
      if (VarAry(i)%MeshID /= MeshID) exit
      associate (iLoc => VarAry(i)%iLoc)
         select case (VarAry(i)%Field)
         case (VF_Force)
            Mesh%Force = reshape(Values(iLoc(1):iLoc(2)), shape(Mesh%Force))
         case (VF_Moment)
            Mesh%Moment = reshape(Values(iLoc(1):iLoc(2)), shape(Mesh%Moment))
         case (VF_TransDisp)
            Mesh%TranslationDisp = reshape(Values(iLoc(1):iLoc(2)), shape(Mesh%TranslationDisp))
         case (VF_Orientation)
            do j = 1, VarAry(i)%Nodes
               Mesh%Orientation(:, :, j) = wm_to_dcm(Values(iLoc(1) + 3*(j - 1):iLoc(1) + 3*j))
            end do
         case (VF_TransVel)
            Mesh%TranslationVel = reshape(Values(iLoc(1):iLoc(2)), shape(Mesh%TranslationVel))
         case (VF_AngularVel)
            Mesh%RotationVel = reshape(Values(iLoc(1):iLoc(2)), shape(Mesh%RotationVel))
         case (VF_TransAcc)
            Mesh%TranslationAcc = reshape(Values(iLoc(1):iLoc(2)), shape(Mesh%TranslationAcc))
         case (VF_AngularAcc)
            Mesh%RotationAcc = reshape(Values(iLoc(1):iLoc(2)), shape(Mesh%RotationAcc))
         case (VF_Scalar)
            Mesh%Scalars = reshape(Values(iLoc(1):iLoc(2)), shape(Mesh%Scalars))
         end select
      end associate
   end do
end subroutine

subroutine MV_Perturb(Var, iLin, PerturbSign, BaseAry, PerturbAry)
   type(ModVarType), intent(in)     :: Var
   integer(IntKi), intent(in)       :: iLin
   integer(IntKi), intent(in)       :: PerturbSign
   real(R8Ki), intent(in)           :: BaseAry(:)
   real(R8Ki), intent(inout)        :: PerturbAry(:)

   real(R8Ki)                       :: Perturb
   real(R8Ki)                       :: WM(3), rotvec(3)
   integer(IntKi)                   :: i, j

   ! Copy base array to perturbed array
   PerturbAry = BaseAry

   ! Get variable perturbation and combine with sign
   Perturb = Var%Perturb*real(PerturbSign, R8Ki)

   ! Index of perturbation value in array
   i = Var%iLoc(1) + iLin - 1

   ! If variable field is orientation, perturbation is in WM parameters
   if (Var%Field == VF_Orientation) then
      j = mod(iLin - 1, 3)                                        ! component being modified (0, 1, 2)
      rotvec = 0.0_R8Ki                                           ! Init WM perturbation to zero
      rotvec(j + 1) = Perturb                                     ! WM perturbation around X,Y,Z axis
      i = i - j                                                   ! index of start of WM parameters (3)
      WM = PerturbAry(i:i + 2)                                    ! Current WM parameters value
      PerturbAry(i:i + 2) = wm_compose(wm_from_rvec(rotvec), WM)  ! Compose value and perturbation
   else
      PerturbAry(i) = PerturbAry(i) + Perturb                     ! Add perturbation
   end if

end subroutine

subroutine MV_ComputeDiff(VarAry, PosAry, NegAry, DiffAry)
   type(ModVarType), intent(in)  :: VarAry(:)      ! Array of variables
   real(R8Ki), intent(in)        :: PosAry(:)      ! Positive result array
   real(R8Ki), intent(in)        :: NegAry(:)      ! Negative result array
   real(R8Ki), intent(inout)     :: DiffAry(:)     ! Array containing difference
   integer(IntKi)                :: i, j, k
   real(R8Ki)                    :: DeltaWM(3), R(3, 3), C1(3), C2(3)

   ! Loop through variables
   do i = 1, size(VarAry)

      ! If variable field is orientation
      if (VarAry(i)%Field == VF_Orientation) then

         ! Loop through nodes
         do j = 1, VarAry(i)%Nodes

            ! Get vector of indicies of WM rotation parameters in array
            k = VarAry(i)%iLoc(1) + 3*(j - 1)

            ! Compose WM parameters to go from negative to positive array
            DeltaWM = wm_compose((PosAry(k:k + 2)), wm_inv(NegAry(k:k + 2)))

            ! Calculate change in rotation in XYZ in radians
            DiffAry(k:k + 2) = wm_to_rvec(DeltaWM)
         end do

      else

         ! Subtract negative array from positive array
         associate (iLoc => VarAry(i)%iLoc)
            DiffAry(iLoc(1):iLoc(2)) = PosAry(iLoc(1):iLoc(2)) - NegAry(iLoc(1):iLoc(2))
         end associate
      end if
   end do
end subroutine

subroutine MV_ComputeCentralDiff(VarAry, Delta, PosAry, NegAry, DerivAry)
   type(ModVarType), intent(in)  :: VarAry(:)      ! Array of variables
   real(R8Ki), intent(in)        :: Delta          ! Positive perturbation value
   real(R8Ki), intent(in)        :: PosAry(:)      ! Positive perturbation result array
   real(R8Ki), intent(in)        :: NegAry(:)      ! Negative perturbation result array
   real(R8Ki), intent(inout)     :: DerivAry(:)    ! Array containing derivative

   ! Compute difference between all values
   call MV_ComputeDiff(VarAry, PosAry, NegAry, DerivAry)

   ! Divide derivative array by twice delta
   DerivAry = DerivAry/(2.0_R8Ki*Delta)

end subroutine

!-------------------------------------------------------------------------------
! Functions for adding Variables an Modules
!-------------------------------------------------------------------------------

subroutine MV_AddModule(ModAry, ModID, ModAbbr, Instance, ModDT, SolverDT, Vars, ErrStat, ErrMsg)
   type(ModDataType), allocatable, intent(inout)   :: ModAry(:)
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
   if (.not. allocated(ModAry)) allocate (ModAry(0))

   ! Populate ModuleDataType derived type
   ModData = ModDataType(Idx=size(ModAry) + 1, ID=ModID, Abbr=ModAbbr, &
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

   ModAry = [ModAry, ModData]

end subroutine

subroutine MV_AddMeshVar(VarAry, Name, Fields, Mesh, Flags, Perturbs, VarIdx, Active)
   type(ModVarType), allocatable, intent(inout) :: VarAry(:)
   character(*), intent(in)                     :: Name
   integer(IntKi), intent(in)                   :: Fields(:)
   type(MeshType), intent(inout)                :: Mesh
   integer(IntKi), optional, intent(in)         :: Flags
   real(R8Ki), optional, intent(in)             :: Perturbs(:)
   integer(IntKi), optional, intent(out)        :: VarIdx
   logical, optional, intent(in)                :: Active
   integer(IntKi)                               :: FlagsLocal
   logical                                      :: ActiveLocal
   real(R8Ki), allocatable                      :: PerturbsLocal(:)
   integer(IntKi)                               :: i, idx

   ! Initialize variable index (variable is not active or mesh is not commited)
   if (present(VarIdx)) VarIdx = 0

   ! If active argument specified and not active, return
   if (present(Active)) then
      if (.not. Active) return
   end if

   ! If mesh has not been committed, return
   if (.not. Mesh%committed) return

   ! Set mesh ID
   if (allocated(VarAry)) then
      Mesh%ID = size(VarAry) + 1
   else
      Mesh%ID = 1
   end if

   ! If present, set variable index from mesh ID
   if (present(VarIdx)) VarIdx = Mesh%ID

   ! Apply flags if specified
   FlagsLocal = VF_Mesh
   if (present(Flags)) FlagsLocal = ior(FlagsLocal, Flags)

   ! Set perturbations if specified
   PerturbsLocal = [(0.0_R8Ki, i=1, size(Fields))]
   if (present(Perturbs)) PerturbsLocal = Perturbs

   ! Loop through fields in mesh
   do i = 1, size(Fields)

      ! Add variable
      call MV_AddVar(VarAry, Name, Fields(i), VarIdx=idx, &
                     Num=Mesh%Nnodes, &
                     Flags=FlagsLocal, &
                     Perturb=PerturbsLocal(i))

      ! Save mesh ID
      VarAry(size(VarAry))%MeshID = Mesh%ID
   end do
end subroutine

subroutine MV_AddVar(VarAry, Name, Field, Num, Flags, iUsr, jUsr, DerivOrder, Perturb, LinNames, VarIdx, Active)
   type(ModVarType), allocatable, intent(inout) :: VarAry(:)
   character(*), intent(in)                     :: Name
   integer(IntKi), intent(in)                   :: Field
   integer(IntKi), optional, intent(in)         :: Num, Flags, iUsr, jUsr
   real(R8Ki), optional, intent(in)             :: Perturb
   integer(IntKi), optional, intent(in)         :: DerivOrder
   character(*), optional, intent(in)           :: LinNames(:)
   integer(IntKi), optional, intent(out)        :: VarIdx
   logical, optional, intent(in)                :: Active
   integer(IntKi)                               :: i
   type(ModVarType)                             :: Var

   ! If active argument specified and not active, return
   if (present(Active)) then
      if (.not. Active) then
         ! Set variable index to zero if present
         if (present(VarIdx)) VarIdx = 0
         return
      end if
   end if

   ! Initialize var with default values
   Var = ModVarType(Name=Name, Field=Field)

   ! Set optional values
   if (present(Num)) Var%Num = Num
   if (present(Flags)) Var%Flags = Flags
   if (present(iUsr)) Var%iUsr = [iUsr, iUsr + Var%Num - 1]
   if (present(jUsr)) Var%jUsr = jUsr
   if (present(Perturb)) Var%Perturb = Perturb
   if (present(LinNames)) then
      allocate (Var%LinNames(size(LinNames)))
      do i = 1, size(LinNames)
         Var%LinNames(i) = LinNames(i)
      end do
   end if

   ! Set Derivative Order
   if (present(DerivOrder)) then
      Var%DerivOrder = DerivOrder
   else
      select case (Var%Field)
      case (VF_Orientation, VF_TransDisp, VF_AngularDisp)   ! Position/displacement
         Var%DerivOrder = 0
      case (VF_TransVel, VF_AngularVel)                     ! Velocity
         Var%DerivOrder = 1
      case (VF_TransAcc, VF_AngularAcc)                     ! Acceleration
         Var%DerivOrder = 2
      case default
         Var%DerivOrder = -1
      end select
   end if

   ! Append Var to VarArray
   if (allocated(VarAry)) then
      VarAry = [VarAry, Var]
   else
      VarAry = [Var]
   end if

   ! Set variable index if present
   if (present(VarIdx)) VarIdx = size(VarAry)
end subroutine

subroutine MV_InitVarIdx(Vars, VarIdx, FlagFilter, ErrStat, ErrMsg)
   type(ModVarsType), intent(in)       :: Vars
   type(VarsIdxType), intent(out)      :: VarIdx
   integer(IntKi), intent(in)          :: FlagFilter
   integer(IntKi), intent(out)         :: ErrStat
   character(ErrMsgLen), intent(out)   :: ErrMsg

   character(*), parameter             :: RoutineName = 'MV_InitVarIdx'
   integer(IntKi)                      :: ErrStat2
   character(ErrMsgLen)                :: ErrMsg2
   type(ModDataType)                   :: ModData
   integer(IntKi)                      :: ivar, inum

   ! Initialize error return
   ErrStat = ErrID_None
   ErrMsg = ""

   ! Save filter in index
   VarIdx%FlagFilter = FlagFilter

   ! Get number of filtered variables
   VarIdx%Nx = MV_NumVars(Vars%x, FlagFilter)
   VarIdx%Nu = MV_NumVars(Vars%u, FlagFilter)
   VarIdx%Ny = MV_NumVars(Vars%y, FlagFilter)

   ! Allocate index arrays
   call AllocAry(VarIdx%ix, VarIdx%Nx, "ix", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(VarIdx%idx, VarIdx%Nx, "idx", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(VarIdx%ixd, VarIdx%Nxd, "ixd", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(VarIdx%iz, VarIdx%Nz, "iz", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(VarIdx%iu, VarIdx%Nu, "iu", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(VarIdx%iy, VarIdx%Ny, "iy", ErrStat2, ErrMsg2); if (Failed()) return

   ! Get filtered value indices
   call GetIndices(Vars%x, VarIdx%ix, FlagFilter)
   call GetIndices(Vars%xd, VarIdx%ixd, FlagFilter)
   call GetIndices(Vars%z, VarIdx%iz, FlagFilter)
   call GetIndices(Vars%u, VarIdx%iu, FlagFilter)
   call GetIndices(Vars%y, VarIdx%iy, FlagFilter)

   ! Copy state variable indices to state variable derivative indices
   VarIdx%idx = VarIdx%ix

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine GetIndices(VarAry, Indices, Mask)
   type(ModVarType), intent(in)  :: VarAry(:)
   integer(IntKi), intent(in)    :: Mask
   integer(IntKi), intent(inout) :: Indices(:)
   integer(IntKi)                :: i, j, k
   k = 1
   do i = 1, size(VarAry)
      if (.not. MV_HasFlags(VarAry(i), Mask)) cycle
      do j = 0, VarAry(i)%Num - 1
         Indices(k) = VarAry(i)%iLoc(1) + j
         k = k + 1
      end do
   end do
end subroutine

function MV_NumVars(VarAry, FlagFilter) result(Num)
   type(ModVarType), intent(in)           :: VarAry(:)
   integer(IntKi), optional, intent(in)   :: FlagFilter
   integer(IntKi)                         :: Num, i
   if (present(FlagFilter)) then
      Num = 0
      do i = 1, size(VarAry)
         if ((FlagFilter == VF_None) .or. (iand(VarAry(i)%Flags, FlagFilter) /= 0)) Num = Num + VarAry(i)%Num
      end do
   else
      Num = sum(VarAry%Num)
   end if
end function

!-------------------------------------------------------------------------------
! Flag Utilities
!-------------------------------------------------------------------------------

pure logical function MV_HasFlags(Var, Flags)
   type(ModVarType), intent(in)  :: Var
   integer(IntKi), intent(in)    :: Flags
   MV_HasFlags = iand(Var%Flags, Flags) == Flags
end function

subroutine MV_SetFlags(Var, Flags)
   type(ModVarType), intent(inout)  :: Var
   integer(IntKi), intent(in)       :: Flags
   integer(IntKi)                   :: i
   Var%Flags = ior(Var%Flags, Flags)
end subroutine

subroutine MV_UnsetFlags(Var, Flags)
   type(ModVarType), intent(inout)  :: Var
   integer(IntKi), intent(in)       :: Flags
   integer(IntKi)                   :: i
   Var%Flags = iand(Var%Flags, not(Flags))
end subroutine

!-------------------------------------------------------------------------------
! String Utilities
!-------------------------------------------------------------------------------

function IdxStr(i1, i2, i3, i4, i5) result(s)
   integer(IntKi), intent(in)             :: i1
   integer(IntKi), optional, intent(in)   :: i2, i3, i4, i5
   character(100)                         :: s
   if (present(i5)) then
      s = '('//trim(Num2LStr(i1))//','//trim(Num2LStr(i2))//','//trim(Num2LStr(i3))//','//trim(Num2LStr(i4))//','//trim(Num2LStr(i5))//')'
   else if (present(i4)) then
      s = '('//trim(Num2LStr(i1))//','//trim(Num2LStr(i2))//','//trim(Num2LStr(i3))//','//trim(Num2LStr(i4))//')'
   else if (present(i3)) then
      s = '('//trim(Num2LStr(i1))//','//trim(Num2LStr(i2))//','//trim(Num2LStr(i3))//')'
   else if (present(i2)) then
      s = '('//trim(Num2LStr(i1))//','//trim(Num2LStr(i2))//')'
   else
      s = '('//trim(Num2LStr(i1))//')'
   end if
end function

!-------------------------------------------------------------------------------
! Rotation Utilities
!-------------------------------------------------------------------------------

pure function quat_from_dcm(R) result(q)
   real(R8Ki), intent(in)  :: R(3, 3)
   real(R8Ki)              :: q(4), C
   integer(IntKi)          :: j

   q = [(1.0_R8Ki + R(1, 1) - R(2, 2) - R(3, 3)), &
        (1.0_R8Ki - R(1, 1) + R(2, 2) - R(3, 3)), &
        (1.0_R8Ki - R(1, 1) - R(2, 2) + R(3, 3)), &
        (1.0_R8Ki + R(1, 1) + R(2, 2) + R(3, 3))]

   ! Get index of max value in q
   j = maxloc(q, dim=1)

   ! Calculate quaternion from direction cosine matrix
   C = q(j)
   select case (j)
   case (1)
      q = [C, (R(1, 2) + R(2, 1)), (R(3, 1) + R(1, 3)), (R(2, 3) - R(3, 2))]
   case (2)
      q = [(R(1, 2) + R(2, 1)), C, (R(2, 3) + R(3, 2)), (R(3, 1) - R(1, 3))]
   case (3)
      q = [(R(3, 1) + R(1, 3)), (R(2, 3) + R(3, 2)), C, (R(1, 2) - R(2, 1))]
   case (4)
      q = [(R(2, 3) - R(3, 2)), (R(3, 1) - R(1, 3)), (R(1, 2) - R(2, 1)), C]
   end select
   q = q/(2.0_R8Ki*sqrt(C))
   if (q(4) < 0.0_R8Ki) q = -q
end function

pure function quat_to_dcm(q) result(R)
   real(R8Ki), intent(in)  :: q(4)
   real(R8Ki)              :: R(3, 3)
   real(R8Ki)              :: q1, q2, q3, q4
   q1 = q(1); q2 = q(2); q3 = q(3); q4 = q(4)
   R(1, :) = [q4*q4 + q1*q1 - q2*q2 - q3*q3, 2*(q1*q2 + q3*q4), 2*(q1*q3 - q2*q4)]
   R(2, :) = [2*(q1*q2 - q3*q4), q4*q4 - q1*q1 + q2*q2 - q3*q3, 2*(q2*q3 + q1*q4)]
   R(3, :) = [2*(q1*q3 + q2*q4), 2*(q2*q3 - q1*q4), q4*q4 - q1*q1 - q2*q2 + q3*q3]
end function

pure function wm_to_quat(c) result(q)
   real(R8Ki), intent(in)  :: c(3)
   real(R8Ki)              :: q(4)
   real(R8Ki)              :: c0, e0, e(3)
   c0 = 2.0_R8Ki - dot_product(c, c)/8.0_R8Ki
   e0 = c0/(4.0_R8Ki - c0)
   e = c/(4.0_R8Ki - c0)
   q = [e, e0]
end function

pure function wm_from_quat(q) result(c)
   real(R8Ki), intent(in)  :: q(4)
   real(R8Ki)              :: c(3)
   real(R8Ki)              :: e0, e(3)
   e0 = q(4)
   e = q(1:3)
   c = 4.0_R8Ki*e/(1.0_R8Ki + e0)
end function

pure function wm_to_dcm(c) result(R)
   real(R8Ki), intent(in)  :: c(3)
   real(R8Ki)              :: R(3, 3), c0, c1, c2, c3
   integer(IntKi)          :: i, j
   c1 = c(1)
   c2 = c(2)
   c3 = c(3)
   c0 = 2.0_R8Ki - dot_product(c, c)/8.0_R8Ki
   R(:, 1) = [c0*c0 + c1*c1 - c2*c2 - c3*c3, &
              2.0_R8Ki*(c1*c2 - c0*c3), &
              2.0_R8Ki*(c1*c3 + c0*c2)]
   R(:, 2) = [2.0_R8Ki*(c1*c2 + c0*c3), &
              c0*c0 - c1*c1 + c2*c2 - c3*c3, &
              2.0_R8Ki*(c2*c3 - c0*c1)]
   R(:, 3) = [2.0_R8Ki*(c1*c3 - c0*c2), &
              2.0_R8Ki*(c2*c3 + c0*c1), &
              c0*c0 - c1*c1 - c2*c2 + c3*c3]
   R = R/(4.0_R8Ki - c0)**2
end function

pure function wm_from_dcm(dcm) result(c)
   real(R8Ki), intent(in)  :: dcm(3, 3)
   real(R8Ki)              :: pivot(4) ! Trace of the rotation matrix and diagonal elements
   real(R8Ki)              :: sm(0:3)
   real(R8Ki)              :: em
   real(R8Ki)              :: Rr(3, 3), c(3)
   integer                 :: i        ! case indicator

   Rr = transpose(dcm)

   ! mjs--find max value of T := Tr(Rr) and diagonal elements of Rr
   ! This tells us which denominator is largest (and less likely to produce numerical noise)
   pivot = [Rr(1, 1) + Rr(2, 2) + Rr(3, 3), Rr(1, 1), Rr(2, 2), Rr(3, 3)]
   i = maxloc(pivot, 1) - 1 ! our sm array starts at 0, so we need to subtract 1 here to get the correct index

   select case (i)
   case (3)
      sm(0) = Rr(2, 1) - Rr(1, 2)                           !  4 c_0 c_3 t_{r0}
      sm(1) = Rr(1, 3) + Rr(3, 1)                           !  4 c_1 c_3 t_{r0}
      sm(2) = Rr(2, 3) + Rr(3, 2)                           !  4 c_2 c_3 t_{r0}
      sm(3) = 1.0_R8Ki - Rr(1, 1) - Rr(2, 2) + Rr(3, 3)      !  4 c_3 c_3 t_{r0}
   case (2)
      sm(0) = Rr(1, 3) - Rr(3, 1)                           !  4 c_0 c_2 t_{r0}
      sm(1) = Rr(1, 2) + Rr(2, 1)                           !  4 c_1 c_2 t_{r0}
      sm(2) = 1.0_R8Ki - Rr(1, 1) + Rr(2, 2) - Rr(3, 3)      !  4 c_2 c_2 t_{r0}
      sm(3) = Rr(2, 3) + Rr(3, 2)                           !  4 c_3 c_2 t_{r0}
   case (1)
      sm(0) = Rr(3, 2) - Rr(2, 3)                           !  4 c_0 c_1 t_{r0}
      sm(1) = 1.0_R8Ki + Rr(1, 1) - Rr(2, 2) - Rr(3, 3)      !  4 c_1 c_1 t_{r0}
      sm(2) = Rr(1, 2) + Rr(2, 1)                           !  4 c_2 c_1 t_{r0}
      sm(3) = Rr(1, 3) + Rr(3, 1)                           !  4 c_3 c_1 t_{r0}
   case (0)
      sm(0) = 1.0_R8Ki + Rr(1, 1) + Rr(2, 2) + Rr(3, 3)      !  4 c_0 c_0 t_{r0}
      sm(1) = Rr(3, 2) - Rr(2, 3)                           !  4 c_1 c_0 t_{r0}
      sm(2) = Rr(1, 3) - Rr(3, 1)                           !  4 c_2 c_0 t_{r0}
      sm(3) = Rr(2, 1) - Rr(1, 2)                           !  4 c_3 c_0 t_{r0}
   end select

   em = sm(0) + SIGN(2.0_R8Ki*SQRT(sm(i)), sm(0))
   em = 4.0_R8Ki/em                                        ! 1 / ( 4 t_{r0} c_{i} ), assuming 0 <= c_0 < 4 and c_{i} > 0
   c = em*sm(1:3)
end function

pure function wm_to_rvec(c) result(rvec)
   real(R8Ki), intent(in) :: c(3)
   real(R8Ki)             :: phi, m, rvec(3)
   m = sqrt(dot_product(c, c))
   if (m == 0.0_R8Ki) then
      rvec = 0.0_R8Ki
      return
   end if
   phi = 4.0_R8Ki*atan(m/4.0_R8Ki)
   rvec = phi*c/m
end function

pure function wm_from_rvec(rvec) result(c)
   real(R8Ki), intent(in) :: rvec(3)
   real(R8Ki)             :: phi, c(3)
   phi = sqrt(dot_product(rvec, rvec))
   if (phi == 0.0_R8Ki) then
      c = 0.0_R8Ki
      return
   end if
   c = 4.0_R8Ki*tan(phi/4.0_R8Ki)*rvec/phi
end function

pure function wm_compose(p, q) result(r)
   real(R8Ki), intent(in)  :: p(3), q(3)
   real(R8Ki)              :: r(3)
   real(R8Ki)              :: p0, q0, D1, D2
   p0 = 2.0_R8Ki - dot_product(p, p)/8.0_R8Ki
   q0 = 2.0_R8Ki - dot_product(q, q)/8.0_R8Ki
   D1 = (4.0_R8Ki - p0)*(4.0_R8Ki - q0)
   D2 = p0*q0 - dot_product(p, q)
   r = 4.0_R8Ki*(q0*p + p0*q + cross(p, q))
   if (D2 >= 0.0_R8Ki) then
      r = r/(D1 + D2)
   else
      r = -r/(D1 - D2)
   end if
end function

pure function wm_inv(c) result(cinv)
   real(R8Ki), intent(in)  :: c(3)
   real(R8Ki)              :: cinv(3)
   cinv = -c
end function

pure function cross(a, b) result(c)
   real(R8Ki), intent(in) :: a(3), b(3)
   real(R8Ki)             :: c(3)
   c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - b(1)*a(2)]
end function

end module
