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

integer(IntKi), parameter :: &
   LoadFields(*) = [VF_Force, VF_Moment], &
   TransFields(*) = [VF_TransDisp, VF_TransVel, VF_TransAcc], &
   AngularFields(*) = [VF_Orientation, VF_AngularDisp, VF_AngularVel, VF_AngularAcc], &
   MotionFields(*) = [VF_TransDisp, VF_Orientation, VF_TransVel, VF_AngularVel, VF_TransAcc, VF_AngularAcc], &
   MeshFields(*) = [LoadFields, MotionFields]

interface MV_PackVar
   module procedure MV_PackVarR4, MV_PackVarR4Ary
   module procedure MV_PackVarR8, MV_PackVarR8Ary
end interface

interface MV_UnpackVar
   module procedure MV_UnpackVarR4, MV_UnpackVarR4Ary
   module procedure MV_UnpackVarR8, MV_UnpackVarR8Ary
end interface

public :: MV_InitVarsVals, MV_LinkOutputInput, MV_VarIndex, MV_PackMesh, MV_UnpackMesh, MV_PackVar, MV_UnpackVar
public :: MV_ComputeCentralDiff, MV_Perturb, MV_ComputeDiff
public :: MV_AddVar, MV_AddMeshVar, MV_AddModule, SetFlags
public :: LoadFields, MotionFields, TransFields, AngularFields, MeshFields
public :: wm_to_dcm, wm_compose, wm_from_dcm, wm_from_xyz, wm_inv
public :: MV_FieldString

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

subroutine MV_InitVarsVals(Vars, Vals, Linearize, ErrStat, ErrMsg)
   type(ModVarsType), intent(inout)    :: Vars
   type(ModValsType), intent(inout)    :: Vals
   logical, intent(in)                 :: Linearize
   integer(IntKi), intent(out)         :: ErrStat
   character(ErrMsgLen), intent(out)   :: ErrMsg

   character(*), parameter  :: RoutineName = 'MV_InitMod'
   integer(IntKi)           :: ErrStat2
   character(ErrMsgLen)     :: ErrMsg2
   integer(IntKi)           :: i, StartIndex

   ! Initialize error outputs
   ErrStat = ErrID_None
   ErrMsg = ''

   ! Initialize state variables
   StartIndex = 1
   do i = 1, size(Vars%x)
      call ModVarType_Init(Vars%x(i), StartIndex, Linearize, ErrStat2, ErrMsg2); if (Failed()) return
   end do

   ! Initialize input variables
   StartIndex = 1
   do i = 1, size(Vars%u)
      call ModVarType_Init(Vars%u(i), StartIndex, Linearize, ErrStat2, ErrMsg2); if (Failed()) return
   end do

   ! Initialize output variables
   StartIndex = 1
   do i = 1, size(Vars%y)
      call ModVarType_Init(Vars%y(i), StartIndex, Linearize, ErrStat2, ErrMsg2); if (Failed()) return
   end do

   ! Calculate number of variables in group (exclude non linearization vars)
   Vars%Nx = sum(Vars%x%Size, iand(Vars%x%Flags, VF_NoLin) == 0)
   Vars%Nu = sum(Vars%u%Size, iand(Vars%u%Flags, VF_NoLin) == 0)
   Vars%Ny = sum(Vars%y%Size, iand(Vars%y%Flags, VF_NoLin) == 0)

   ! Allocate state, input, and output values
   call AllocAry(Vals%x, Vars%Nx, "Vals%x", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(Vals%dxdt, Vars%Nx, "Vals%dxdt", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(Vals%u, Vars%Nu, "Vals%u", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(Vals%y, Vars%Ny, "Vals%y", ErrStat2, ErrMsg2); if (Failed()) return

   ! Allocate perturbation input and output values
   call AllocAry(Vals%u_perturb, Vars%Nu, "Vals%u_perturb", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(Vals%x_perturb, Vars%Nx, "Vals%x_perturb", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(Vals%xp, Vars%Nx, "Vals%xp", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(Vals%xn, Vars%Nx, "Vals%xn", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(Vals%yp, Vars%Ny, "Vals%yp", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(Vals%yn, Vars%Ny, "Vals%yn", ErrStat2, ErrMsg2); if (Failed()) return

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
   integer(IntKi)                :: nNodes
   character(1), parameter       :: Comp(3) = ['X', 'Y', 'Z']
   character(*), parameter       :: Fmt = '(A," ",A,", node",I0,", ",A)'
   character(2)                  :: UnitDesc

   ! Initialize error outputs
   ErrStat = ErrID_None
   ErrMsg = ''

   !----------------------------------------------------------------------------
   ! Basic Variable
   !----------------------------------------------------------------------------

   Var%NumLin = Var%Size

   !----------------------------------------------------------------------------
   ! Mesh
   !----------------------------------------------------------------------------

   ! If this variable belongs to a mesh
   if (iand(Var%Flags, VF_Mesh) > 0) then

      ! Size is the number of nodes in a mesh
      Var%Nodes = Var%Size

      ! Number of linearization values
      Var%NumLin = Var%Nodes*3
      Var%Size = Var%Nodes*3

      ! If linearization requested
      if (Linearize) then

         ! Set unit description for line mesh
         UnitDesc = ''
         if (iand(Var%Flags, VF_Line) > 0) UnitDesc = "/m"

         ! Switch based on field number
         select case (Var%Field)
         case (VF_Force)
            Var%LinNames = [character(LinChanLen) ::((trim(Var%Name)//" "//Comp(j)//" force, node "//trim(num2lstr(i))//', N'//UnitDesc, j=1, 3), i=1, nNodes)]
         case (VF_Moment)
            Var%LinNames = [character(LinChanLen) ::((trim(Var%Name)//" "//Comp(j)//" moment, node "//trim(num2lstr(i))//', Nm'//UnitDesc, j=1, 3), i=1, nNodes)]
         case (VF_TransDisp)
            Var%LinNames = [character(LinChanLen) ::((trim(Var%Name)//" "//Comp(j)//" translation displacement, node "//trim(num2lstr(i))//', m', j=1, 3), i=1, nNodes)]
         case (VF_Orientation)
            Var%LinNames = [character(LinChanLen) ::((trim(Var%Name)//" "//Comp(j)//" orientation angle, node "//trim(num2lstr(i))//', rad', j=1, 3), i=1, nNodes)]
         case (VF_TransVel)
            Var%LinNames = [character(LinChanLen) ::((trim(Var%Name)//" "//Comp(j)//" translation velocity, node "//trim(num2lstr(i))//', m/s', j=1, 3), i=1, nNodes)]
         case (VF_AngularVel)
            Var%LinNames = [character(LinChanLen) ::((trim(Var%Name)//" "//Comp(j)//" rotation velocity, node "//trim(num2lstr(i))//', rad/s', j=1, 3), i=1, nNodes)]
         case (VF_TransAcc)
            Var%LinNames = [character(LinChanLen) ::((trim(Var%Name)//" "//Comp(j)//" translation acceleration, node "//trim(num2lstr(i))//', m/s^2', j=1, 3), i=1, nNodes)]
         case (VF_AngularAcc)
            Var%LinNames = [character(LinChanLen) ::((trim(Var%Name)//" "//Comp(j)//" rotation acceleration, node "//trim(num2lstr(i))//', rad/s^2', j=1, 3), i=1, nNodes)]
         case default
            call SetErrStat(ErrID_Fatal, "Invalid mesh field type", ErrStat, ErrMsg, RoutineName)
            return
         end select

      end if
   end if

   !----------------------------------------------------------------------------
   ! No Linearization
   !----------------------------------------------------------------------------

   if (iand(Var%Flags, VF_NoLin) > 0) then   ! No Linearization

      ! Number of linearization values is zero if NoLin flag is set
      Var%NumLin = 0

   else if (Linearize) then

      ! If insufficient linearization names, return error
      if (size(Var%LinNames) < Var%NumLin) then
         call SetErrStat(ErrID_Fatal, "insufficient LinNames given for "//Var%Name, ErrStat, ErrMsg, RoutineName)
         return
      end if
   end if

   !----------------------------------------------------------------------------
   ! Indices
   !----------------------------------------------------------------------------

   ! Initialize local index
   call AllocAry(Var%iLoc, Var%Size, "Var%iLoc", ErrStat2, ErrMsg2); if (Failed()) return
   Var%iLoc = [(index + i, i=0, Var%Size - 1)]

   ! Update index based on variable size
   index = index + Var%Size

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

subroutine MV_PackVarR4(VarAry, iVar, Val, Ary)
   type(ModVarType), intent(in)  :: VarAry(:)
   integer(IntKi), intent(inout) :: iVar
   real(R4Ki), intent(in)        :: Val
   real(R8Ki), intent(inout)     :: Ary(:)
   Ary(VarAry(iVar)%iLoc(1)) = real(Val, R8Ki)
   iVar = iVar + 1
end subroutine

subroutine MV_PackVarR8(VarAry, iVar, Val, Ary)
   type(ModVarType), intent(in)  :: VarAry(:)
   integer(IntKi), intent(inout) :: iVar
   real(R8Ki), intent(in)        :: Val
   real(R8Ki), intent(inout)     :: Ary(:)
   Ary(VarAry(iVar)%iLoc(1)) = Val
   iVar = iVar + 1
end subroutine

subroutine MV_PackVarR4Ary(VarAry, iVar, Val, Ary)
   type(ModVarType), intent(in)  :: VarAry(:)
   integer(IntKi), intent(inout) :: iVar
   real(R4Ki), intent(in)        :: Val(:)
   real(R8Ki), intent(inout)     :: Ary(:)
   Ary(VarAry(iVar)%iLoc) = real(pack(Val, .true.), R4Ki)
   iVar = iVar + 1
end subroutine

subroutine MV_PackVarR8Ary(VarAry, iVar, Vals, Ary)
   type(ModVarType), intent(in)  :: VarAry(:)
   integer(IntKi), intent(inout) :: iVar
   real(R8Ki), intent(in)        :: Vals(:)
   real(R8Ki), intent(inout)     :: Ary(:)
   Ary(VarAry(iVar)%iLoc) = pack(Vals, .true.)
   iVar = iVar + 1
end subroutine

subroutine MV_UnpackVarR4(VarAry, iVar, Ary, Val)
   type(ModVarType), intent(in)  :: VarAry(:)
   integer(IntKi), intent(inout) :: iVar
   real(R4Ki), intent(in)        :: Ary(:)
   real(R8Ki), intent(inout)     :: Val
   Val = Ary(VarAry(iVar)%iLoc(1))
   iVar = iVar + 1
end subroutine

subroutine MV_UnpackVarR4Ary(VarAry, iVar, Ary, Vals)
   type(ModVarType), intent(in)  :: VarAry(:)
   integer(IntKi), intent(inout) :: iVar
   real(R4Ki), intent(in)        :: Ary(:)
   real(R8Ki), intent(inout)     :: Vals(:)
   Vals = Ary(VarAry(iVar)%iLoc)
   iVar = iVar + 1
end subroutine

subroutine MV_UnpackVarR8(VarAry, iVar, Ary, Vals)
   type(ModVarType), intent(in)  :: VarAry(:)
   integer(IntKi), intent(inout) :: iVar
   real(R8Ki), intent(in)        :: Ary(:)
   real(R8Ki), intent(inout)     :: Vals
   Vals = Ary(VarAry(iVar)%iLoc(1))
   iVar = iVar + 1
end subroutine

subroutine MV_UnpackVarR8Ary(VarAry, iVar, Ary, Vals)
   type(ModVarType), intent(in)  :: VarAry(:)
   integer(IntKi), intent(inout) :: iVar
   real(R8Ki), intent(in)        :: Ary(:)
   real(R8Ki), intent(inout)     :: Vals(:)
   Vals = Ary(VarAry(iVar)%iLoc)
   iVar = iVar + 1
end subroutine

subroutine MV_PackMesh(VarAry, iVar, Mesh, Values)
   type(ModVarType), intent(in)  :: VarAry(:)
   integer(IntKi), intent(inout) :: iVar
   type(MeshType), intent(in)    :: Mesh
   real(R8Ki), intent(inout)     :: Values(:)
   character(VarNameLen)         :: MeshName
   integer(IntKi)                :: j
   MeshName = VarAry(iVar)%Name
   do while (VarAry(iVar)%Name == MeshName)
      select case (VarAry(iVar)%Field)
      case (VF_Force)
         Values(VarAry(iVar)%iLoc) = pack(Mesh%Force, .true.)
      case (VF_Moment)
         Values(VarAry(iVar)%iLoc) = pack(Mesh%Moment, .true.)
      case (VF_TransDisp)
         Values(VarAry(iVar)%iLoc) = pack(Mesh%TranslationDisp, .true.)
      case (VF_Orientation)
         do j = 1, VarAry(iVar)%Nodes
            Values(VarAry(iVar)%iLoc(3*(j - 1) + 1:3*j)) = wm_from_dcm(Mesh%Orientation(:, :, j))
         end do
      case (VF_TransVel)
         Values(VarAry(iVar)%iLoc) = pack(Mesh%TranslationVel, .true.)
      case (VF_AngularVel)
         Values(VarAry(iVar)%iLoc) = pack(Mesh%RotationVel, .true.)
      case (VF_TransAcc)
         Values(VarAry(iVar)%iLoc) = pack(Mesh%TranslationAcc, .true.)
      case (VF_AngularAcc)
         Values(VarAry(iVar)%iLoc) = pack(Mesh%RotationAcc, .true.)
      case (VF_Scalar)
         Values(VarAry(iVar)%iLoc) = pack(Mesh%Scalars, .true.)
      end select
      iVar = iVar + 1
      if (iVar > size(VarAry)) exit
   end do
end subroutine

subroutine MV_UnpackMesh(VarAry, iVar, Values, Mesh)
   type(ModVarType), intent(in)  :: VarAry(:)
   integer(IntKi), intent(inout) :: iVar
   real(R8Ki), intent(in)        :: Values(:)
   type(MeshType), intent(inout) :: Mesh
   character(VarNameLen)         :: MeshName
   integer(IntKi)                :: j
   MeshName = VarAry(iVar)%Name
   do while (VarAry(iVar)%Name == MeshName)
      select case (VarAry(iVar)%Field)
      case (VF_Force)
         Mesh%Force = reshape(Values(VarAry(iVar)%iLoc), shape(Mesh%Force))
      case (VF_Moment)
         Mesh%Moment = reshape(Values(VarAry(iVar)%iLoc), shape(Mesh%Moment))
      case (VF_TransDisp)
         Mesh%TranslationDisp = reshape(Values(VarAry(iVar)%iLoc), shape(Mesh%TranslationDisp))
      case (VF_Orientation)
         do j = 1, VarAry(iVar)%Nodes
            Mesh%Orientation(:, :, j) = wm_to_dcm(Values(VarAry(iVar)%iLoc(3*(j - 1) + 1:3*j)))
         end do
      case (VF_TransVel)
         Mesh%TranslationVel = reshape(Values(VarAry(iVar)%iLoc), shape(Mesh%TranslationVel))
      case (VF_AngularVel)
         Mesh%RotationVel = reshape(Values(VarAry(iVar)%iLoc), shape(Mesh%RotationVel))
      case (VF_TransAcc)
         Mesh%TranslationAcc = reshape(Values(VarAry(iVar)%iLoc), shape(Mesh%TranslationAcc))
      case (VF_AngularAcc)
         Mesh%RotationAcc = reshape(Values(VarAry(iVar)%iLoc), shape(Mesh%RotationAcc))
      case (VF_Scalar)
         Mesh%Scalars = reshape(Values(VarAry(iVar)%iLoc), shape(Mesh%Scalars))
      end select
      iVar = iVar + 1
      if (iVar > size(VarAry)) exit
   end do
end subroutine

subroutine MV_Perturb(Var, iLin, PerturbSign, BaseAry, PerturbAry, iPerturb)
   type(ModVarType), intent(in)     :: Var
   integer(IntKi), intent(in)       :: iLin
   integer(IntKi), intent(in)       :: PerturbSign
   real(R8Ki), intent(in)           :: BaseAry(:)
   real(R8Ki), intent(inout)        :: PerturbAry(:)
   integer(IntKi), intent(out)      :: iPerturb
   real(R8Ki)                       :: Perturb
   real(R8Ki)                       :: WM(3), WMp(3)
   integer(IntKi)                   :: i, j, iLoc(3)

   ! Copy base array to perturbed array
   PerturbAry = BaseAry

   ! Get variable perturbation and combine with sign
   Perturb = Var%Perturb*real(PerturbSign, R8Ki)

   ! Perturbation index within array
   iPerturb = Var%iLoc(iLin)

   ! If variable field is orientation, perturbation is in WM parameters
   if (Var%Field == VF_Orientation) then
      j = mod(iLin - 1, 3)                         ! component being modified (0, 1, 2)
      i = iLin - j                                 ! index of start of WM parameters (3)
      iLoc = Var%iLoc(i:i + 2)                     ! array index vector
      WMp = 0.0_R8Ki                               ! Init WM perturbation to zero
      WMp(j + 1) = 4.0_R8Ki*tan(Perturb/4.0_R8Ki)  ! WM perturbation around X,Y,Z axis
      ! WMp(j + 1) = Perturb                       ! WM perturbation around X,Y,Z axis
      WM = PerturbAry(iLoc)                        ! Current WM parameters value
      PerturbAry(iLoc) = wm_compose(WM, WMp)       ! Compose value and perturbation
   else
      PerturbAry(Var%iLoc(iLin)) = PerturbAry(Var%iLoc(iLin)) + Perturb
   end if

end subroutine

subroutine MV_ComputeDiff(VarAry, PosAry, NegAry, DiffAry)
   type(ModVarType), intent(in)  :: VarAry(:)      ! Array of variables
   real(R8Ki), intent(in)        :: PosAry(:)      ! Positive result array
   real(R8Ki), intent(in)        :: NegAry(:)      ! Negative result array
   real(R8Ki), intent(inout)     :: DiffAry(:)     ! Array containing difference
   integer(IntKi)                :: i, j, ind(3)
   real(R8Ki)                    :: DeltaWM(3)
   
   ! Loop through variables
   do i = 1, size(VarAry)

      ! If variable field is orientation
      if (VarAry(i)%Field == VF_Orientation) then

         ! Loop through nodes
         do j = 1, VarAry(i)%Nodes

            ! Get vector of indicies of WM rotation parameters in array
            ind = VarAry(i)%iLoc(3*(j - 1) + 1:3*j)

            ! Compose WM parameters to go from negative to positive array
            DeltaWM = wm_compose(wm_inv(NegAry(ind)), PosAry(ind))

            ! Calculate change in rotation in XYZ in radians
            DiffAry(ind) = 4.0_R8Ki*atan(DeltaWM/4.0_R8Ki)
         end do

      else

         ! Subtract negative array from positive array
         DiffAry(VarAry(i)%iLoc) = PosAry(VarAry(i)%iLoc) - NegAry(VarAry(i)%iLoc)
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

   ! Populate ModuleDataType derived type
   ModData = ModDataType(ID=ModID, Abbr=ModAbbr, Ins=Instance, DT=ModDT, Vars=Vars)

   ! Allocate mapping index with zero length
   call AllocAry(ModData%iMapsOpt1, 0, "ModData%iMapsOpt1", ErrStat2, ErrMsg2)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return
   call AllocAry(ModData%iMapsOpt2, 0, "ModData%iMapsOpt2", ErrStat2, ErrMsg2)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return
   call AllocAry(ModData%iMapsAll, 0, "ModData%iMapsAll", ErrStat2, ErrMsg2)
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

   if (allocated(ModAry)) then
      ModAry = [ModAry, ModData]
   else
      ModAry = [ModData]
   end if
end subroutine

subroutine MV_AddMeshVar(VarAry, Name, Fields, Nodes, Flags, Perturbs, Active)
   type(ModVarType), allocatable, intent(inout) :: VarAry(:)
   character(*), intent(in)                     :: Name
   integer(IntKi), intent(in)                   :: Fields(:)
   integer(IntKi), intent(in)                   :: Nodes
   integer(IntKi), optional, intent(in)         :: Flags
   real(R8Ki), optional, intent(in)             :: Perturbs(:)
   logical, optional, intent(in)                :: Active
   integer(IntKi)                               :: i, FlagsLocal
   logical                                      :: ActiveLocal
   real(R8Ki), allocatable                      :: PerturbsLocal(:)
   FlagsLocal = 0
   if (present(Flags)) FlagsLocal = Flags
   FlagsLocal = ior(FlagsLocal, VF_Mesh)
   PerturbsLocal = [(0.0_R8Ki, i=1, size(Fields))]
   if (present(Perturbs)) PerturbsLocal = Perturbs
   ActiveLocal = .true.
   if (present(Active)) ActiveLocal = Active
   do i = 1, size(Fields)
      call MV_AddVar(VarAry, Name, Fields(i), Num=Nodes, Flags=FlagsLocal, &
                     Perturb=PerturbsLocal(i), Active=ActiveLocal)
   end do
end subroutine

subroutine MV_AddVar(VarAry, Name, Field, Num, Flags, iUsr, Perturb, LinNames, Active)
   type(ModVarType), allocatable, intent(inout) :: VarAry(:)
   character(*), intent(in)                     :: Name
   integer(IntKi), intent(in)                   :: Field
   integer(IntKi), optional, intent(in)         :: Num, Flags, iUsr(:)
   real(R8Ki), optional, intent(in)             :: Perturb
   logical, optional, intent(in)                :: Active
   character(*), optional, intent(in)           :: LinNames(:)
   integer(IntKi)                               :: i
   type(ModVarType)                             :: Var

   ! If active argument specified and not active, return
   if (present(Active)) then
      if (.not. Active) return
   end if

   ! Initialize var with default values
   Var = ModVarType(Name=Name, Field=Field)

   ! Set optional values
   if (present(Num)) Var%Size = Num
   if (present(Flags)) Var%Flags = Flags
   if (present(iUsr)) Var%iUsr = iUsr
   if (present(Perturb)) Var%Perturb = Perturb
   if (present(LinNames)) then
      allocate (Var%LinNames(size(LinNames)))
      do i = 1, size(LinNames)
         Var%LinNames(i) = LinNames(i)
      end do
   end if

   ! Set Derivative Order
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

   ! Append Var to VarArray
   if (allocated(VarAry)) then
      VarAry = [VarAry, Var]
   else
      VarAry = [Var]
   end if
end subroutine

! Get index of variable in array matching name and field
function MV_VarIndex(VarAry, Name, Field) result(Indx)
   type(ModVarType), intent(in)  :: VarAry(:)
   character(*), intent(in)      :: Name
   integer(IntKi), intent(in)    :: Field
   integer(IntKi)                :: Indx
   do Indx = 1, size(VarAry)
      if (string_equal_ci(VarAry(Indx)%Name, Name) .and. &
          VarAry(Indx)%Field == Field) exit
   end do
   if (Indx > size(VarAry)) Indx = 0
end function

!-------------------------------------------------------------------------------
! Functions for linking variables (Output and Input)
!-------------------------------------------------------------------------------

subroutine MV_LinkOutputInput(OutVars, InpVars, OutName, InpName, Field, ErrStat, ErrMsg)
   type(ModVarsType), intent(inout)    :: OutVars, InpVars
   character(*), intent(in)            :: OutName, InpName
   integer(IntKi), intent(in)          :: Field
   integer(IntKi), intent(out)         :: ErrStat
   character(ErrMsgLen), intent(out)   :: ErrMsg

   character(*), parameter  :: RoutineName = 'MV_LinkOutputInput'
   ! integer(IntKi)           :: ErrStat2
   ! character(ErrMsgLen)     :: ErrMsg2
   ! integer(IntKi)           :: i
   integer(IntKi)           :: InpVarIndex, OutVarIndex

   ! Initialize error outputs
   ErrStat = ErrID_None
   ErrMsg = ''

   ! Find name/field in input vars
   InpVarIndex = MV_VarIndex(InpVars%u, InpName, Field)
   if (InpVarIndex == 0) then
      call SetErrStat(ErrID_Fatal, 'Input variable "'//InpName//'" with field '// &
                      trim(num2lstr(Field))//' not found', ErrStat, ErrMsg, RoutineName)
      return
   end if

   ! Find name/field in output vars
   OutVarIndex = MV_VarIndex(OutVars%u, OutName, Field)
   if (OutVarIndex == 0) then
      call SetErrStat(ErrID_Fatal, 'Output variable "'//OutName//'" with field '// &
                      trim(num2lstr(Field))//' not found', ErrStat, ErrMsg, RoutineName)
      return
   end if

   ! If error finding input or output variable, return
   if (ErrStat >= AbortErrLev) return

   ! TODO: figure out what to do here

end subroutine

!-------------------------------------------------------------------------------
! Flag Utilities
!-------------------------------------------------------------------------------

subroutine SetFlags(Var, Mask)
   type(ModVarType), intent(inout)  :: Var
   integer(IntKi), intent(in)       :: Mask
   integer(IntKi)                   :: i
   Var%Flags = ior(Var%Flags, Mask)
end subroutine

!-------------------------------------------------------------------------------
! String Utilities
!-------------------------------------------------------------------------------

! Compare strings s1 and s2 while ignoring case
function string_equal_ci(s1, s2) result(is_equal)
   character(*), intent(in)   :: s1, s2
   logical                    :: is_equal
   integer(IntKi), parameter  :: ca = iachar("a")
   integer(IntKi), parameter  :: cz = iachar("z")
   integer(IntKi)             :: i, j
   integer(IntKi)             :: c1, c2
   is_equal = .false.
   i = len_trim(s1)
   j = len_trim(s2)
   if (i /= j) return
   do i = 1, j
      c1 = iachar(s1(i:i))
      c2 = iachar(s2(i:i))
      if (c1 == c2) cycle
      if (c1 >= ca .and. c1 <= cz) c1 = c1 - 32
      if (c2 >= ca .and. c2 <= cz) c2 = c2 - 32
      if (c1 /= c2) return
   end do
   is_equal = .true.
end function

!-------------------------------------------------------------------------------
! Rotation Utilities
!-------------------------------------------------------------------------------

pure function dcm_to_quat(R) result(q)
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

! pure function wm_to_dcm(c) result(R)
!    real(R8Ki), intent(in)  :: c(3)
!    real(R8Ki)              :: R(3, 3)
!    R = quat_to_dcm(wm_to_quat(c))
! end function

pure function wm_to_dcm(c) result(R)
   real(R8Ki), intent(in)  :: c(3)
   real(R8Ki)              :: R(3, 3), cct, F(3, 3)
   integer(IntKi)          :: i, j
   cct = dot_product(c, c)
   F = reshape([0.0_R8Ki, -c(3), c(2), c(3), 0.0_R8Ki, -c(1), -c(2), c(1), 0.0_R8Ki], [3, 3])/2.0_R8Ki
   do i = 1, 3
      F(i, i) = F(i, i) + 1.0_R8Ki - cct/16.0_R8Ki
      do j = 1, 3
         F(i, j) = F(i, j) + c(i)*c(j)/8.0_R8Ki
      end do
   end do
   F = F/(1.0_R8Ki + cct/16.0_R8Ki)
   R = matmul(F, F)
end function

pure function wm_from_dcm(R) result(c)
   real(R8Ki), intent(in)  :: R(3, 3)
   real(R8Ki)              :: c(3), cct
   c = wm_from_quat(dcm_to_quat(R))
   cct = dot_product(c, c)
   if (cct > 16.0_R8Ki) c = 16.0_R8Ki*c/cct
end function

! pure function wm_from_dcm(R) result(c)
!    real(R8Ki), intent(in)  :: R(3, 3)
!    real(R8Ki)              :: c(3), t1, t2, cct
!    t1 = 1.0_R8Ki + R(1,1) + R(2,2) + R(3,3)
!    t2 = 2.0_R8Ki*sqrt(t1)
!    c(1) = (R(3,2) - R(2,3))
!    c(2) = (R(1,3) - R(3,1))
!    c(3) = (R(2,1) - R(1,2))
!    c = 4.0_R8Ki * c / (t1 + t2)
!    cct = dot_product(c,c)
!    if (cct > 16.0_R8Ki) c = 16.0_R8Ki*c / cct
! end function

pure function wm_compose(p, q) result(r)
   real(R8Ki), intent(in)  :: p(3), q(3)
   real(R8Ki)              :: r(3)
   real(R8Ki)              :: p0, q0, D1, D2
   p0 = 2.0_R8Ki - dot_product(p, p)/8.0_R8Ki
   q0 = 2.0_R8Ki - dot_product(q, q)/8.0_R8Ki
   D1 = (4.0_R8Ki - p0)*(4.0_R8Ki - q0)
   D2 = p0*q0 - dot_product(p, q)
   if (D2 >= 0.0_R8Ki) then
      r = 4*(q0*p + p0*q + cross(p, q))/(D1 + D2)
   else
      r = -4*(q0*p + p0*q + cross(p, q))/(D1 - D2)
   end if
end function

pure function wm_to_zyx(c) result(zyx)
   real(R8Ki), intent(in)  :: c(3)
   real(R8Ki)              :: zyx(3)
   real(R8Ki)              :: q(4), qx, qy, qz, qw
   q = wm_to_quat(c)
   qx = q(1); qy = q(2); qz = q(3); qw = q(4)
   zyx(1) = atan2(2*(qw*qx + qy*qz), 1.0_R8Ki - 2.0_R8Ki*(qx*qx + qy*qy))
   zyx(2) = -PiBy2_D + 2.0_R8Ki*atan2(sqrt(1.0_R8Ki + 2.0_R8Ki*(qw*qy - qx*qz)), &
                                      sqrt(1.0_R8Ki - 2.0_R8Ki*(qw*qy - qx*qz)))
   zyx(3) = atan2(2.0_R8Ki*(qw*qz + qx*qy), 1.0_R8Ki - 2.0_R8Ki*(qy*qy + qz*qz))
end function

function wm_to_xyz(c) result(xyz)
   real(R8Ki), intent(in)  :: c(3)
   real(R8Ki)              :: xyz(3)
   xyz = EulerExtract(wm_to_dcm(c))
end function

pure function wm_inv(c) result(cinv)
   real(R8Ki), intent(in)  :: c(3)
   real(R8Ki)              :: cinv(3)
   cinv = -c
end function

pure function wm_from_xyz(xyz) result(c)
   real(R8Ki), intent(in)  :: xyz(3)
   real(R8Ki)              :: c(3)
   real(R8Ki)              :: n(3)
   c = 0.0_R8Ki
   c = wm_compose([4, 0, 0]*tan(xyz(1)/4.0_R8Ki), c)  ! X
   c = wm_compose([0, 4, 0]*tan(xyz(2)/4.0_R8Ki), c)  ! Y
   c = wm_compose([0, 0, 4]*tan(xyz(3)/4.0_R8Ki), c)  ! Z
end function

pure function cross(a, b) result(c)
   real(R8Ki), intent(in) :: a(3), b(3)
   real(R8Ki)             :: c(3)
   c = [a(2)*b(3) - a(3)*b(2), -a(3)*b(1) + a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
end function

end module
