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
use NWTC_Num
use ModMesh
implicit none

private
public :: MV_InitVarsJac, MV_Pack, MV_Unpack
public :: MV_ComputeCentralDiff, MV_Perturb, MV_ComputeDiff, MV_ExtrapInterp, MV_AddDelta
public :: MV_AddVar, MV_AddMeshVar
public :: MV_HasFlags, MV_SetFlags, MV_ClearFlags, MV_NumVars, MV_NumVals, MV_FindVarDatLoc
public :: LoadFields, MotionFields, TransFields, AngularFields
public :: quat_to_dcm, dcm_to_quat, quat_inv, quat_to_rvec, rvec_to_quat, wm_to_quat, quat_to_wm, wm_inv
public :: MV_FieldString, MV_IsLoad, IdxStr
public :: DumpMatrix, MV_AddModule
public :: MV_EqualDL

integer(IntKi), parameter :: &
   LoadFields(*) = [FieldForce, FieldMoment], &
   TransFields(*) = [FieldTransDisp, FieldTransVel, FieldTransAcc], &
   AngularFields(*) = [FieldOrientation, FieldAngularVel, FieldAngularAcc, FieldAngularDisp], &
   MotionFields(*) = [FieldTransDisp, FieldOrientation, FieldTransVel, FieldAngularVel, FieldTransAcc, FieldAngularAcc]

interface MV_Pack
   module procedure MV_PackVarRank0R4, MV_PackVarRank1R4
   module procedure MV_PackVarRank0R8, MV_PackVarRank1R8
   module procedure MV_PackMesh
end interface

interface MV_Unpack
   module procedure MV_UnpackVarRank0R4, MV_UnpackVarRank1R4
   module procedure MV_UnpackVarRank0R8, MV_UnpackVarRank1R8
   module procedure MV_UnpackMesh
end interface

logical, parameter   :: UseSmallRotAngles = .true.

contains

!-------------------------------------------------------------------------------
! MV_Pack
!-------------------------------------------------------------------------------

subroutine MV_PackVarRank0R4(Var, SrcVal, DstAry)
   type(ModVarType), intent(in)  :: Var
   real(R4Ki), intent(in)        :: SrcVal
   real(R8Ki), intent(inout)     :: DstAry(:)
   DstAry(Var%iLoc(1)) = real(SrcVal, R8Ki)
end subroutine

subroutine MV_PackVarRank0R8(Var, SrcVal, DstAry)
   type(ModVarType), intent(in)  :: Var
   real(R8Ki), intent(in)        :: SrcVal
   real(R8Ki), intent(inout)     :: DstAry(:)
   DstAry(Var%iLoc(1)) = SrcVal
end subroutine

subroutine MV_PackVarRank1R4(Var, SrcAry, DstAry)
   type(ModVarType), intent(in)  :: Var
   real(R4Ki), intent(in)        :: SrcAry(:)
   real(R8Ki), intent(inout)     :: DstAry(:)
   DstAry(Var%iLoc(1):Var%iLoc(2)) = real(SrcAry, R8Ki)
end subroutine

subroutine MV_PackVarRank1R8(Var, SrcAry, Ary)
   type(ModVarType), intent(in)  :: Var
   real(R8Ki), intent(in)        :: SrcAry(:)
   real(R8Ki), intent(inout)     :: Ary(:)
   Ary(Var%iLoc(1):Var%iLoc(2)) = SrcAry
end subroutine

subroutine MV_PackMesh(Var, Mesh, DstAry)
   type(ModVarType), intent(in)  :: Var
   type(MeshType), intent(in)    :: Mesh
   real(R8Ki), intent(inout)     :: DstAry(:)
   integer(IntKi)                :: i, j, k
   select case (Var%Field)
   case (FieldForce)
      DstAry(Var%iLoc(1):Var%iLoc(2)) = pack(real(Mesh%Force, R8Ki), .true.)
   case (FieldMoment)
      DstAry(Var%iLoc(1):Var%iLoc(2)) = pack(real(Mesh%Moment, R8Ki), .true.)
   case (FieldTransDisp)
      DstAry(Var%iLoc(1):Var%iLoc(2)) = pack(real(Mesh%TranslationDisp, R8Ki), .true.)
   case (FieldOrientation)
      k = Var%iLoc(1)
      do j = 1, Var%Nodes
         DstAry(k:k + 2) = dcm_to_quat(Mesh%Orientation(:, :, j))
         k = k + 3
      end do
   case (FieldTransVel)
      DstAry(Var%iLoc(1):Var%iLoc(2)) = pack(real(Mesh%TranslationVel, R8Ki), .true.)
   case (FieldAngularVel)
      DstAry(Var%iLoc(1):Var%iLoc(2)) = pack(real(Mesh%RotationVel, R8Ki), .true.)
   case (FieldTransAcc)
      DstAry(Var%iLoc(1):Var%iLoc(2)) = pack(real(Mesh%TranslationAcc, R8Ki), .true.)
   case (FieldAngularAcc)
      DstAry(Var%iLoc(1):Var%iLoc(2)) = pack(real(Mesh%RotationAcc, R8Ki), .true.)
   case (FieldScalar)
      DstAry(Var%iLoc(1):Var%iLoc(2)) = pack(real(Mesh%Scalars, R8Ki), .true.)
   end select
end subroutine

!-------------------------------------------------------------------------------
! MV_Unpack
!-------------------------------------------------------------------------------

subroutine MV_UnpackVarRank0R4(Var, SrcAry, DstVal)
   type(ModVarType), intent(in)  :: Var
   real(R8Ki), intent(in)        :: SrcAry(:)
   real(R4Ki), intent(inout)     :: DstVal
   DstVal = real(SrcAry(Var%iLoc(1)), R4Ki)
end subroutine

subroutine MV_UnpackVarRank0R8(Var, SrcAry, DstVal)
   type(ModVarType), intent(in)  :: Var
   real(R8Ki), intent(in)        :: SrcAry(:)
   real(R8Ki), intent(inout)     :: DstVal
   DstVal = SrcAry(Var%iLoc(1))
end subroutine

subroutine MV_UnpackVarRank1R4(Var, SrcAry, DstAry)
   type(ModVarType), intent(in)  :: Var
   real(R8Ki), intent(in)        :: SrcAry(:)
   real(R4Ki), intent(inout)     :: DstAry(:)
   DstAry = real(SrcAry(Var%iLoc(1):Var%iLoc(2)), R4Ki)
end subroutine

subroutine MV_UnpackVarRank1R8(Var, SrcAry, DstAry)
   type(ModVarType), intent(in)  :: Var
   real(R8Ki), intent(in)        :: SrcAry(:)
   real(R8Ki), intent(inout)     :: DstAry(:)
   DstAry = SrcAry(Var%iLoc(1):Var%iLoc(2))
end subroutine

subroutine MV_UnpackMesh(Var, SrcAry, Mesh)
   type(ModVarType), intent(in)  :: Var
   real(R8Ki), intent(in)        :: SrcAry(:)
   type(MeshType), intent(inout) :: Mesh
   integer(IntKi)                :: i, j, k
   select case (Var%Field)
   case (FieldForce)
      Mesh%Force = reshape(SrcAry(Var%iLoc(1):Var%iLoc(2)), shape(Mesh%Force))
   case (FieldMoment)
      Mesh%Moment = reshape(SrcAry(Var%iLoc(1):Var%iLoc(2)), shape(Mesh%Moment))
   case (FieldTransDisp)
      Mesh%TranslationDisp = reshape(SrcAry(Var%iLoc(1):Var%iLoc(2)), shape(Mesh%TranslationDisp))
   case (FieldOrientation)
      k = Var%iLoc(1)
      do j = 1, Var%Nodes
         Mesh%Orientation(:, :, j) = quat_to_dcm(SrcAry(k:k + 2))
         k = k + 3
      end do
   case (FieldTransVel)
      Mesh%TranslationVel = reshape(SrcAry(Var%iLoc(1):Var%iLoc(2)), shape(Mesh%TranslationVel))
   case (FieldAngularVel)
      Mesh%RotationVel = reshape(SrcAry(Var%iLoc(1):Var%iLoc(2)), shape(Mesh%RotationVel))
   case (FieldTransAcc)
      Mesh%TranslationAcc = reshape(SrcAry(Var%iLoc(1):Var%iLoc(2)), shape(Mesh%TranslationAcc))
   case (FieldAngularAcc)
      Mesh%RotationAcc = reshape(SrcAry(Var%iLoc(1):Var%iLoc(2)), shape(Mesh%RotationAcc))
   case (FieldScalar)
      Mesh%Scalars = reshape(SrcAry(Var%iLoc(1):Var%iLoc(2)), shape(Mesh%Scalars))
   end select
end subroutine

!-------------------------------------------------------------------------------
! Field Names
!-------------------------------------------------------------------------------

function MV_FieldString(Field) result(str)
   integer(IntKi), intent(in) :: Field
   character(16)              :: str
   select case (Field)
   case (FieldAngularAcc)
      str = "FieldAngularAcc"
   case (FieldAngularDisp)
      str = "FieldAngularDisp"
   case (FieldAngularVel)
      str = "FieldAngularVel"
   case (FieldForce)
      str = "FieldForce"
   case (FieldMoment)
      str = "FieldMoment"
   case (FieldOrientation)
      str = "FieldOrientation"
   case (FieldTransAcc)
      str = "FieldTransAcc"
   case (FieldTransDisp)
      str = "FieldTransDisp"
   case (FieldTransVel)
      str = "FieldTransVel"
   case default
      str = "Unknown"
   end select
end function

subroutine MV_InitVarsJac(Vars, Jac, Linearize, ErrStat, ErrMsg)
   type(ModVarsType), intent(inout)    :: Vars
   type(ModJacType), intent(inout)     :: Jac
   logical, intent(in)                 :: Linearize
   integer(IntKi), intent(out)         :: ErrStat
   character(ErrMsgLen), intent(out)   :: ErrMsg

   character(*), parameter       :: RoutineName = 'MV_InitVarsJac'
   integer(IntKi)                :: ErrStat2
   character(ErrMsgLen)          :: ErrMsg2
   integer(IntKi)                :: i, StartIndex

   ! Initialize error outputs
   ErrStat = ErrID_None
   ErrMsg = ''

   ! Initialize number of variables in each group
   Vars%Nx = 0
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
   Jac%Nx = Vars%Nx

   ! Initialize constraint state variables
   if (.not. allocated(Vars%z)) allocate (Vars%z(0))
   StartIndex = 1
   do i = 1, size(Vars%z)
      call ModVarType_Init(Vars%z(i), StartIndex, Linearize, ErrStat2, ErrMsg2)
      if (Failed()) return
   end do
   Vars%Nz = sum(Vars%z%Num)
   Jac%Nz = Vars%Nz

   ! Initialize input variables
   if (.not. allocated(Vars%u)) allocate (Vars%u(0))
   StartIndex = 1
   do i = 1, size(Vars%u)
      call ModVarType_Init(Vars%u(i), StartIndex, Linearize, ErrStat2, ErrMsg2)
      if (Failed()) return
   end do
   Vars%Nu = sum(Vars%u%Num)
   Jac%Nu = Vars%Nu

   ! Initialize output variables
   if (.not. allocated(Vars%y)) allocate (Vars%y(0))
   StartIndex = 1
   do i = 1, size(Vars%y)
      call ModVarType_Init(Vars%y(i), StartIndex, Linearize, ErrStat2, ErrMsg2)
      if (Failed()) return
   end do
   Vars%Ny = sum(Vars%y%Num)
   Jac%Ny = Vars%Ny

   ! Allocate Jacobian data arrays
   if (Linearize) then
      if (Jac%Nx > 0) then
         call AllocAry(Jac%x, Jac%Nx, "Lin%x", ErrStat2, ErrMsg2); if (Failed()) return
         call AllocAry(Jac%x_perturb, Jac%Nx, "Lin%x_perturb", ErrStat2, ErrMsg2); if (Failed()) return
         call AllocAry(Jac%x_pos, Jac%Nx, "Lin%x_pos", ErrStat2, ErrMsg2); if (Failed()) return
         call AllocAry(Jac%x_neg, Jac%Nx, "Lin%x_neg", ErrStat2, ErrMsg2); if (Failed()) return
      end if
      if (Jac%Nz > 0) then
         call AllocAry(Jac%z, Jac%Nz, "Lin%z", ErrStat2, ErrMsg2); if (Failed()) return
      end if
      if (Jac%Nu > 0) then
         call AllocAry(Jac%u, Jac%Nu, "Lin%u", ErrStat2, ErrMsg2); if (Failed()) return
         call AllocAry(Jac%u_perturb, Jac%Nu, "Lin%u_perturb", ErrStat2, ErrMsg2); if (Failed()) return
      end if
      if (Jac%Ny > 0) then
         call AllocAry(Jac%y, Jac%Ny, "Lin%y", ErrStat2, ErrMsg2); if (Failed()) return
         call AllocAry(Jac%y_pos, Jac%Ny, "Lin%y_pos", ErrStat2, ErrMsg2); if (Failed()) return
         call AllocAry(Jac%y_neg, Jac%Ny, "Lin%y_neg", ErrStat2, ErrMsg2); if (Failed()) return
      end if
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
         case (FieldForce)
            Var%LinNames = [character(LinChanLen) ::((trim(Var%Name)//" "//Comp(j)//" force, node "//trim(num2lstr(i))//', N'//UnitDesc, j=1, 3), i=1, Var%Nodes)]
         case (FieldMoment)
            Var%LinNames = [character(LinChanLen) ::((trim(Var%Name)//" "//Comp(j)//" moment, node "//trim(num2lstr(i))//', Nm'//UnitDesc, j=1, 3), i=1, Var%Nodes)]
         case (FieldTransDisp)
            Var%LinNames = [character(LinChanLen) ::((trim(Var%Name)//" "//Comp(j)//" translation displacement, node "//trim(num2lstr(i))//', m', j=1, 3), i=1, Var%Nodes)]
         case (FieldOrientation)
            Var%LinNames = [character(LinChanLen) ::((trim(Var%Name)//" "//Comp(j)//" orientation angle, node "//trim(num2lstr(i))//', rad', j=1, 3), i=1, Var%Nodes)]
         case (FieldTransVel)
            Var%LinNames = [character(LinChanLen) ::((trim(Var%Name)//" "//Comp(j)//" translation velocity, node "//trim(num2lstr(i))//', m/s', j=1, 3), i=1, Var%Nodes)]
         case (FieldAngularVel)
            Var%LinNames = [character(LinChanLen) ::((trim(Var%Name)//" "//Comp(j)//" rotation velocity, node "//trim(num2lstr(i))//', rad/s', j=1, 3), i=1, Var%Nodes)]
         case (FieldTransAcc)
            Var%LinNames = [character(LinChanLen) ::((trim(Var%Name)//" "//Comp(j)//" translation acceleration, node "//trim(num2lstr(i))//', m/s^2', j=1, 3), i=1, Var%Nodes)]
         case (FieldAngularAcc)
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

subroutine MV_AddModule(ModDataAry, ModID, ModAbbr, Instance, ModDT, SolverDT, Vars, Linearize, ErrStat, ErrMsg)
   type(ModDataType), allocatable, intent(inout)   :: ModDataAry(:)
   integer(IntKi), intent(in)                      :: ModID
   character(*), intent(in)                        :: ModAbbr
   integer(IntKi), intent(in)                      :: Instance
   real(R8Ki), intent(in)                          :: ModDT
   real(R8Ki), intent(in)                          :: SolverDT
   type(ModVarsType), intent(in)                   :: Vars
   logical, intent(in)                             :: Linearize
   integer(IntKi), intent(out)                     :: ErrStat
   character(ErrMsgLen), intent(out)               :: ErrMsg

   character(*), parameter                         :: RoutineName = 'MV_AddModule'
   integer(IntKi)                                  :: ErrStat2
   character(ErrMsgLen)                            :: ErrMsg2
   type(ModDataType)                               :: ModData
   integer(IntKi)                                  :: i, StartIndex

   ErrStat = ErrID_None
   ErrMsg = ''

   ! Populate module information
   if (allocated(ModDataAry)) then
      ModData%iMod = size(ModDataAry) + 1
   else
      ModData%iMod = 1
   end if
   ModData%ID = ModID
   ModData%Abbr = ModAbbr
   ModData%Ins = Instance
   ModData%DT = ModDT
   call NWTC_Library_CopyModVarsType(Vars, ModData%Vars, MESH_NEWCOPY, ErrStat2, ErrMsg2); if (Failed()) return

   !----------------------------------------------------------------------------
   ! Calculate Module Sub-stepping
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
   ! Add module info to array
   !----------------------------------------------------------------------------

   if (.not. allocated(ModDataAry)) then
      ModDataAry = [ModData]
   else
      ModDataAry = [ModDataAry, ModData]
   end if

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine GetModuleOrder(ModDataAry, ModIDs, ModOrder)
   type(ModDataType), intent(in)             :: ModDataAry(:)     !< Array of module data structures
   integer(IntKi), intent(in)                :: ModIDs(:)   !< List of module IDs to keep in order
   integer(IntKi), allocatable, intent(out)  :: ModOrder(:) !< Module data indices in order of ModIDs
   integer(IntKi), allocatable               :: ModIDAry(:), indices(:)
   integer(IntKi)                            :: i

   ! Create array 1 to size(Mod) representing the index of each module data
   indices = [(i, i=1, size(ModDataAry))]

   ! Get array of module IDs from array of module data
   ModIDAry = [(ModDataAry(i)%ID, i=1, size(ModDataAry))]

   ! Initialize module order array with no size
   allocate (ModOrder(0))

   ! Loop through module IDs to keep, add module data indices that match module ID to order array
   do i = 1, size(ModIDs)
      ModOrder = [ModOrder, pack(indices, ModIDAry == ModIDs(i))]
   end do

end subroutine

!-------------------------------------------------------------------------------
! Functions for packing and unpacking data by variable
!-------------------------------------------------------------------------------


subroutine MV_Perturb(Var, iLin, PerturbSign, BaseAry, PerturbAry)
   type(ModVarType), intent(in)     :: Var
   integer(IntKi), intent(in)       :: iLin
   integer(IntKi), intent(in)       :: PerturbSign
   real(R8Ki), intent(in)           :: BaseAry(:)
   real(R8Ki), intent(inout)        :: PerturbAry(:)

   real(R8Ki)                       :: Perturb
   real(R8Ki)                       :: quat(3), quat_p(3)
   real(R8Ki)                       :: rv(3), dcm(3,3)
   integer(IntKi)                   :: i, j
   integer(IntKi)                   :: ErrStat
   character(ErrMsgLen)             :: ErrMsg

   ! Copy base array to perturbed array
   PerturbAry = BaseAry

   ! Get variable perturbation and combine with sign
   Perturb = Var%Perturb*real(PerturbSign, R8Ki)

   ! Index of perturbation value in array
   i = Var%iLoc(1) + iLin - 1

   ! If variable field is orientation, perturbation is in radians
   if (Var%Field == FieldOrientation) then
      j = mod(iLin - 1, 3)                      ! component being modified (0, 1, 2)
      i = i - j                                 ! index of start of quaternion parameters (3)
      quat = BaseAry(i:i + 2)                   ! Current quat parameters value
      if (MV_HasFlags(Var, VF_SmallAngle)) then
         dcm = quat_to_dcm(quat)
         rv = GetSmllRotAngs(dcm, ErrStat, ErrMsg)
         rv(j+1) = rv(j+1) + Perturb
         call SmllRotTrans('linearization perturbation', rv(1), rv(2), rv(3), dcm, ErrStat=ErrStat, ErrMsg=ErrMsg)
         quat = dcm_to_quat(dcm)
      else
         quat_p = perturb_quat(Perturb, j + 1)     ! Quaternion of perturbed angle
         quat = quat_compose(quat, quat_p)         ! Compose perturbation and current rotation
      end if
      PerturbAry(i:i + 2) = quat                ! Save perturbed quaternion in array
   else
      PerturbAry(i) = PerturbAry(i) + Perturb   ! Add perturbation directly
   end if

end subroutine

subroutine MV_ComputeDiff(VarAry, PosAry, NegAry, DiffAry)
   type(ModVarType), intent(in)  :: VarAry(:)      ! Array of variables
   real(R8Ki), intent(in)        :: PosAry(:)      ! Positive result array
   real(R8Ki), intent(in)        :: NegAry(:)      ! Negative result array
   real(R8Ki), intent(inout)     :: DiffAry(:)     ! Array containing difference
   integer(IntKi)                :: i, j, k
   real(R8Ki)                    :: delta(3), R(3, 3), quat_pos(3), quat_neg(3)
   real(R8Ki)                    :: ang_pos(3), ang_neg(3)
   integer(IntKi)                :: ErrStat
   character(ErrMsgLen)          :: ErrMsg

   ! Loop through variables
   do i = 1, size(VarAry)

      ! If variable field is orientation
      if (VarAry(i)%Field == FieldOrientation) then

         ! Starting index into arrays
         k = VarAry(i)%iLoc(1)

         ! Loop through nodes
         do j = 1, VarAry(i)%Nodes

            ! Quaternions from negative and positive perturbations
            quat_neg = NegAry(k:k + 2)
            quat_pos = PosAry(k:k + 2)

            ! If flag set to use small angle rotations
            if (UseSmallRotAngles) then

               ! If variable has flag to use small angles when computing difference
               if (MV_HasFlags(VarAry(i), VF_SmallAngle)) then

                  ang_pos = GetSmllRotAngs(quat_to_dcm(quat_pos), ErrStat, ErrMsg)
                  ang_neg = GetSmllRotAngs(quat_to_dcm(quat_neg), ErrStat, ErrMsg)

                  DiffAry(k:k + 2) = ang_pos - ang_neg
               else

                  ! Calculate relative rotation from negative to positive perturbation
                  delta = quat_compose(-quat_neg, quat_pos)

                  ! Convert relative rotation from quaternion to rotation vector
                  DiffAry(k:k + 2) = GetSmllRotAngs(quat_to_dcm(delta), ErrStat, ErrMsg)
               end if

            else

               ! Calculate relative rotation from negative to positive perturbation
               delta = quat_compose(-quat_neg, quat_pos)

               ! Convert delta quaternion to rotation vector and store in diff array
               DiffAry(k:k + 2) = quat_to_rvec(delta)

            end if

            ! Increment starting index
            k = k + 3

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

!> MV_ExtrapInterp interpolates arrays of variable data to the target x value from
!! the array of x values. Supports constant, linear, and quadratic interpolation
!! similar to the ExtrapInterp routines created by the registry.
subroutine MV_ExtrapInterp(VarAry, y, tin, y_out, tin_out, ErrStat, ErrMsg)
   type(ModVarType), intent(in)  :: VarAry(:)      ! Array of variables
   real(R8Ki), intent(in)        :: y(:, :)
   real(R8Ki), intent(in)        :: tin(:)
   real(R8Ki), intent(inout)     :: y_out(:)
   real(R8Ki), intent(in)        :: tin_out
   integer(IntKi), intent(out)   :: ErrStat
   character(*), intent(out)     :: ErrMsg

   character(*), parameter       :: RoutineName = 'MV_ExtrapInterp'
   integer(IntKi)                :: ErrStat2
   character(ErrMsgLen)          :: ErrMsg2
   integer(IntKi)                :: InterpOrder
   real(R8Ki)                    :: t(3), t_out, a1, a2, a3
   real(R8Ki)                    :: q1(4), q2(4), q3(4), q(4)
   integer(IntKi)                :: i, j, k

   ErrStat = ErrID_None
   ErrMsg = ''

   ! Check that array sizes match
   if (size(t) /= size(y, 2)) then
      call SetErrStat(ErrID_Fatal, 'size(t) must equal size(y)', ErrStat, ErrMsg, RoutineName)
      return
   end if

   ! Calculate interpolation order
   InterpOrder = size(t) - 1

   ! Switch based on interpolation order
   select case (InterpOrder)

   case (0) ! Constant interpolation (copy)

      y_out = y(:, 1)

   case (1) ! Linear Interpolation

      t(1:2) = tin - tin(1)
      t_out = tin_out - tin(1)
      a1 = -(t_out - t(2))/t(2)
      a2 = t_out/t(2)
      y_out = a1*y(:, 1) + a2*y(:, 2)

      ! Loop through glue output variables
      do i = 1, size(VarAry)

         ! Switch based on variable field type
         select case (VarAry(i)%Field)

         case (FieldOrientation)   ! SLERP for orientation quaternions

            k = VarAry(i)%iLoc(1)
            do j = 1, VarAry(i)%Nodes

               ! Get quaternion 1 from array, calculate scalar
               q1(2:4) = y(k:k + 2, 1)
               q1(1) = quat_scalar(q1(2:4))

               ! Get quaternion 2 from array, calculate scalar
               q2(2:4) = y(k:k + 2, 2)
               q2(1) = quat_scalar(q2(2:4))

               ! Calculate dot product of two quaternions
               ! Make quaternion 2 consistent with quaternion 1 for interp
               if (dot_product(q1, q2) < 0.0_R8Ki) q2 = -q2

               ! Interpolate quaternion components
               q = a1*q1 + a2*q2

               ! Store canonical quaternion in output array
               y_out(k:k + 2) = quat_canonical(q(1), q(2:4))

               ! Increment quaternion index
               k = k + 3
            end do

         case (FieldScalar) ! Scalar field

            ! If field is on the range [0,2PI], perform angular interp
            if (MV_HasFlags(VarAry(i), VF_2PI)) then

               k = VarAry(i)%iLoc(1)
               do j = 1, VarAry(i)%Num
                  call Angles_ExtrapInterp(y(k, 1), y(k, 2), t(1:2), y_out(k), t_out)
                  k = k + 1
               end do

            end if

         end select

      end do

   case (2) ! Quadratic Interpolation

      t = tin - tin(1)
      t_out = tin_out - tin(1)
      a1 = (t_out - t(2))*(t_out - t(3))/((t(1) - t(2))*(t(1) - t(3)))
      a2 = (t_out - t(1))*(t_out - t(3))/((t(2) - t(1))*(t(2) - t(3)))
      a3 = (t_out - t(1))*(t_out - t(2))/((t(3) - t(1))*(t(3) - t(2)))
      y_out = a1*y(:, 1) + a2*y(:, 2) + a3*y(:, 3)

      ! Loop through glue output variables
      do i = 1, size(VarAry)

         ! Switch based on variable field type
         select case (VarAry(i)%Field)

         case (FieldOrientation)   ! SLERP for orientation quaternions

            k = VarAry(i)%iLoc(1)
            do j = 1, VarAry(i)%Nodes

               ! Get quaternion 1 from array, calculate scalar
               q1(2:4) = y(k:k + 2, 1)
               q1(1) = quat_scalar(q1(2:4))

               ! Get quaternion 2 from array, calculate scalar
               q2(2:4) = y(k:k + 2, 2)
               q2(1) = quat_scalar(q2(2:4))

               ! Get quaternion 3 from array, calculate scalar
               q3(2:4) = y(k:k + 2, 3)
               q3(1) = quat_scalar(q2(2:4))

               ! Make quaternions 2 and 3 consistent with quaternion 1
               if (dot_product(q1, q2) < 0.0_R8Ki) q2 = -q2
               if (dot_product(q1, q3) < 0.0_R8Ki) q3 = -q3

               ! Interpolate quaternion components
               q = a1*q1 + a2*q2 + a3*q3

               ! Store canonical quaternion in output array
               y_out(k:k + 2) = quat_canonical(q(1), q(2:4))

               ! Increment quaternion index
               k = k + 3
            end do

         case (FieldScalar) ! Scalar field

            ! If field is on the range [0,2PI], perform angular interp
            if (MV_HasFlags(VarAry(i), VF_2PI)) then

               k = VarAry(i)%iLoc(1)
               do j = 1, VarAry(i)%Num
                  call Angles_ExtrapInterp(y(k, 1), y(k, 2), y(k, 3), t, y_out(k), t_out)
                  k = k + 1
               end do

            end if

         end select

      end do

   case default

      ! Unsupported Interpolation
      call SetErrStat(ErrID_Fatal, 'size(t) must be less than 4 (order must be less than 3).', ErrStat, ErrMsg, RoutineName)
      return
   end select

end subroutine

subroutine MV_AddDelta(VarAry, DeltaAry, DataAry)
   type(ModVarType), intent(in)  :: VarAry(:)      ! Array of variables
   real(R8Ki), intent(in)        :: DeltaAry(:)    ! Array of delta values
   real(R8Ki), intent(inout)     :: DataAry(:)     ! Array to be modified
   integer(IntKi)                :: i, j, k
   real(R8Ki)                    :: quat_base(3), quat_delta(3)

   ! Loop through variables
   do i = 1, size(VarAry)
      associate (iLoc => VarAry(i)%iLoc)
         select case (VarAry(i)%Field)
         case (FieldOrientation)

            ! Starting index into arrays
            k = iLoc(1)

            ! Loop through nodes
            do j = 1, VarAry(i)%Nodes

               ! Quaternions from negative and positive perturbations
               quat_base = DataAry(k:k + 2)
               quat_delta = rvec_to_quat(DeltaAry(k:k + 2))

               ! Calculate composition of base quaternion and delta quaternion
               DataAry(k:k + 2) = quat_compose(quat_base, quat_delta)

               ! Increment starting index
               k = k + 3
            end do

         case default
            DataAry(iLoc(1):iLoc(2)) = DataAry(iLoc(1):iLoc(2)) + DeltaAry(iLoc(1):iLoc(2))
         end select
      end associate
   end do
end subroutine

!-------------------------------------------------------------------------------
! Functions for adding Variables
!-------------------------------------------------------------------------------

subroutine MV_AddMeshVar(VarAry, Name, Fields, DL, Mesh, Flags, Perturbs, Active, iVar)
   type(ModVarType), allocatable, intent(inout) :: VarAry(:)
   character(*), intent(in)                     :: Name
   integer(IntKi), intent(in)                   :: Fields(:)
   type(DatLoc), intent(in)                     :: DL
   type(MeshType), intent(inout)                :: Mesh
   integer(IntKi), optional, intent(in)         :: Flags
   real(R8Ki), optional, intent(in)             :: Perturbs(:)
   logical, optional, intent(in)                :: Active
   integer(IntKi)                               :: FlagsLocal
   logical                                      :: ActiveLocal
   real(R8Ki), allocatable                      :: PerturbsLocal(:)
   integer(IntKi), optional, intent(out)        :: iVar
   integer(IntKi)                               :: i

   ! If variable index is present, initialize to zero in case variable is inactive
   if (present(iVar)) iVar = 0

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

   ! Save variable index
   if (present(iVar)) iVar = Mesh%ID

   ! Apply flags if specified
   FlagsLocal = VF_Mesh
   if (present(Flags)) FlagsLocal = ior(FlagsLocal, Flags)

   ! Set perturbations if specified
   PerturbsLocal = [(0.0_R8Ki, i=1, size(Fields))]
   if (present(Perturbs)) PerturbsLocal = Perturbs

   ! Loop through fields in mesh
   do i = 1, size(Fields)

      ! Add variable
      call MV_AddVar(VarAry, Name, Fields(i), &
                     DL=DL, &
                     Num=Mesh%Nnodes, &
                     Flags=FlagsLocal, &
                     Perturb=PerturbsLocal(i))

      ! Save mesh ID
      VarAry(size(VarAry))%MeshID = Mesh%ID
   end do
end subroutine

subroutine MV_AddVar(VarAry, Name, Field, DL, Num, iAry, jAry, kAry, Flags, DerivOrder, Perturb, LinNames, Active)
   type(ModVarType), allocatable, intent(inout) :: VarAry(:)
   character(*), intent(in)                     :: Name
   integer(IntKi), intent(in)                   :: Field
   type(DatLoc), intent(in)                     :: DL
   integer(IntKi), optional, intent(in)         :: iAry, jAry, kAry
   integer(IntKi), optional, intent(in)         :: Num, Flags
   real(R8Ki), optional, intent(in)             :: Perturb
   integer(IntKi), optional, intent(in)         :: DerivOrder
   character(*), optional, intent(in)           :: LinNames(:)
   logical, optional, intent(in)                :: Active
   integer(IntKi)                               :: i
   type(ModVarType)                             :: Var

   ! If active argument specified and not active, return
   if (present(Active)) then
      if (.not. Active) return
   end if

   ! Initialize var with default values
   Var = ModVarType(Name=Name, Field=Field, DL=DL, Num=1)

   ! If number of values is zero, return
   if (present(Num)) then
      if (Num == 0) return
      Var%Num = Num
   end if

   ! Set optional values
   if (present(Flags)) Var%Flags = Flags
   if (present(iAry)) then
      Var%iAry = [iAry, iAry + Var%Num - 1]
   else
      Var%iAry = [1, Var%Num]
   end if
   if (present(jAry)) Var%jAry = jAry
   if (present(kAry)) Var%kAry = kAry
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
      case (FieldOrientation, FieldTransDisp, FieldAngularDisp)   ! Position/displacement
         Var%DerivOrder = 0
      case (FieldTransVel, FieldAngularVel)                     ! Velocity
         Var%DerivOrder = 1
      case (FieldTransAcc, FieldAngularAcc)                     ! Acceleration
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

end subroutine

function MV_NumVals(VarAry, FlagFilter) result(Num)
   type(ModVarType), intent(in)           :: VarAry(:)
   integer(IntKi), optional, intent(in)   :: FlagFilter
   integer(IntKi)                         :: Num, i
   if (present(FlagFilter)) then
      Num = 0
      do i = 1, size(VarAry)
         if (MV_HasFlags(VarAry(i), FlagFilter)) Num = Num + VarAry(i)%Num
      end do
   else
      Num = sum(VarAry%Num)
   end if
end function

function MV_NumVars(VarAry, FlagFilter) result(Num)
   type(ModVarType), intent(in)           :: VarAry(:)
   integer(IntKi), optional, intent(in)   :: FlagFilter
   integer(IntKi)                         :: Num, i
   if (present(FlagFilter)) then
      Num = 0
      do i = 1, size(VarAry)
         if (MV_HasFlags(VarAry(i), FlagFilter)) Num = Num + 1
      end do
   else
      Num = size(VarAry)
   end if
end function

! MV_IsLoad returns true if the variable field is FieldForce or FieldMoment
pure logical function MV_IsLoad(Var)
   type(ModVarType), intent(in)  :: Var
   MV_IsLoad = Var%Field == FieldForce .or. Var%Field == FieldMoment
end function

! MV_EqualDL returns true if data location numbers are greater than zero and
! all components of the data location are the same.
pure logical function MV_EqualDL(DL1, DL2)
   type(DatLoc), intent(in)   :: DL1, DL2
   MV_EqualDL = DL1%Num /= 0 .and. DL2%Num /= 0 .and. &
                DL1%Num == DL2%Num .and. &
                DL1%i1 == DL2%i1 .and. &
                DL1%i2 == DL2%i2 .and. &
                DL1%i3 == DL2%i3
end function

! Find variable index in array based on DatLoc number
pure function MV_FindVarDatLoc(VarAry, DL) result(iVar)
   type(ModVarType), intent(in)  :: VarAry(:)
   type(DatLoc), intent(in)      :: DL
   integer(IntKi)                :: iVar
   do iVar = 1, size(VarAry)
      if (VarAry(iVar)%DL%Num /= DL%Num) cycle
      if (VarAry(iVar)%DL%i1 /= DL%i1) cycle
      if (VarAry(iVar)%DL%i2 /= DL%i2) cycle
      if (VarAry(iVar)%DL%i3 /= DL%i3) cycle
      return
   end do
   iVar = 0
end function

!-------------------------------------------------------------------------------
! Flag Utilities
!-------------------------------------------------------------------------------

!> MV_HasFlags returns true if Flags is VF_None or if variable contains all
!> flags in Flags.
pure logical function MV_HasFlags(Var, Flags)
   type(ModVarType), intent(in)  :: Var
   integer(IntKi), intent(in)    :: Flags
   MV_HasFlags = iand(Var%Flags, Flags) == Flags
end function

!> MV_SetFlags adds the given flags to the variable.
subroutine MV_SetFlags(Var, Flags)
   type(ModVarType), intent(inout)  :: Var
   integer(IntKi), intent(in)       :: Flags
   integer(IntKi)                   :: i
   Var%Flags = ior(Var%Flags, Flags)
end subroutine

!> MV_ClearFlags removes the given flags from the variable.
subroutine MV_ClearFlags(Var, Flags)
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

function perturb_quat(theta, idir) result(q)
   real(R8Ki), intent(in)     :: theta
   integer(IntKi), intent(in) :: idir
   real(R8Ki)                 :: rvec(3), q(3), dcm(3, 3)
   integer(IntKi)             :: ErrStat
   character(ErrMsgLen)       :: ErrMsg

   select case (idir)
   case (1)
      rvec = [theta, 0.0_R8Ki, 0.0_R8Ki]
   case (2)
      rvec = [0.0_R8Ki, theta, 0.0_R8Ki]
   case (3)
      rvec = [0.0_R8Ki, 0.0_R8Ki, theta]
   end select

   if (UseSmallRotAngles) then
      call SmllRotTrans('linearization perturbation', rvec(1), rvec(2), rvec(3), dcm, ErrStat=ErrStat, ErrMsg=ErrMsg)
      q = dcm_to_quat(dcm)
   else
      q = rvec_to_quat(rvec)
   end if
end function

pure function quat_scalar(q) result(w)
   real(R8Ki), intent(in)  :: q(3)
   real(R8Ki)              :: im, w
   ! Calculate magnitude of imaginary part of quaternion
   im = dot_product(q, q)
   if (im < 1.0_R8Ki) then
      w = sqrt(1.0_R8Ki - im)
   else if (im > 1.0_R8Ki) then
      w = 0.0_R8Ki
   else
      w = 0.0_R8Ki
   end if
end function

pure function quat_canonical(q0, q) result(qc)
   real(R8Ki), intent(in)  :: q0, q(3)
   real(R8Ki)              :: qc(3), m
   integer(IntKi)          :: i
   m = q0*q0 + dot_product(q, q)
   qc = q/m
   if (q0 < 0.0_R8Ki) then
      qc = -q/m
   else
      qc = q/m
   end if
   ! if (q0 > 0.0_R8Ki) return
   ! if (q0 < 0.0_R8Ki) then
   !    qc = -qc
   !    return
   ! end if
   ! do i = 1, 3
   !    if (q(i) > 0.0_R8Ki) return
   !    if (q(i) < 0.0_R8Ki) then
   !       qc = -qc
   !       return
   !    end if
   ! end do
end function

function dcm_to_quat(dcm) result(q)
   real(R8Ki), intent(in)  :: dcm(3, 3)
   real(R8Ki)              :: q(3)
   real(R8Ki)              :: t, s, qw

   ! Trace of matrix
   t = dcm(1, 1) + dcm(2, 2) + dcm(3, 3)

   if (t > 0.0_R8Ki) then
      S = sqrt(t + 1.0_R8Ki)*2.0_R8Ki  ! S=4*qw
      qw = 0.25_R8Ki*S
      q(1) = (dcm(3, 2) - dcm(2, 3))/S
      q(2) = (dcm(1, 3) - dcm(3, 1))/S
      q(3) = (dcm(2, 1) - dcm(1, 2))/S
   elseif ((dcm(1, 1) > dcm(2, 2)) .and. (dcm(1, 1) > dcm(3, 3))) then
      S = sqrt(1.0_R8Ki + dcm(1, 1) - dcm(2, 2) - dcm(3, 3))*2.0_R8Ki  ! S=4*qx
      qw = (dcm(3, 2) - dcm(2, 3))/S
      q(1) = 0.25_R8Ki*S
      q(2) = (dcm(1, 2) + dcm(2, 1))/S
      q(3) = (dcm(1, 3) + dcm(3, 1))/S
   elseif (dcm(2, 2) > dcm(3, 3)) then
      S = sqrt(1.0_R8Ki + dcm(2, 2) - dcm(1, 1) - dcm(3, 3))*2.0_R8Ki  ! S=4*qy
      qw = (dcm(1, 3) - dcm(3, 1))/S
      q(1) = (dcm(1, 2) + dcm(2, 1))/S
      q(2) = 0.25_R8Ki*S
      q(3) = (dcm(2, 3) + dcm(3, 2))/S
   else
      S = sqrt(1.0_R8Ki + dcm(3, 3) - dcm(1, 1) - dcm(2, 2))*2.0_R8Ki  ! S=4*qz
      qw = (dcm(2, 1) - dcm(1, 2))/S
      q(1) = (dcm(1, 3) + dcm(3, 1))/S
      q(2) = (dcm(2, 3) + dcm(3, 2))/S
      q(3) = 0.25_R8Ki*S
   end if

   q = quat_canonical(qw, q)
end function

! dcm_to_quat2 returns a quaternion from a DCM based on eigenanalysis
! https://en.wikipedia.org/wiki/Rotation_matrix#Quaternion
function dcm_to_quat2(dcm) result(q)
   real(R8Ki), intent(in)     :: dcm(3, 3)
   real(R8Ki)                 :: q(3)
   integer(IntKi), parameter  :: n = 4
   real(R8Ki)                 :: Qxx, Qxy, Qxz, Qyx, Qyy, Qyz, Qzx, Qzy, Qzz
   real(R8Ki)                 :: A(n, n), wr(n), wi(n), vl(n, n), vr(n, n), work(4*n)
   integer(IntKi)             :: info, lwork, i

   Qxx = dcm(1, 1)
   Qyx = dcm(2, 1)
   Qzx = dcm(3, 1)
   Qxy = dcm(1, 2)
   Qyy = dcm(2, 2)
   Qzy = dcm(3, 2)
   Qxz = dcm(1, 3)
   Qyz = dcm(2, 3)
   Qzz = dcm(3, 3)

   A(:, 1) = [Qxx - Qyy - Qzz, Qyx + Qxy, Qzx + Qxz, Qzy - Qyz]/3.0_R8Ki
   A(:, 2) = [Qyx + Qxy, Qyy - Qxx - Qzz, Qzy + Qyz, Qxz - Qzx]/3.0_R8Ki
   A(:, 3) = [Qzx + Qxz, Qzy + Qyz, Qzz - Qxx - Qyy, Qyx - Qxy]/3.0_R8Ki
   A(:, 4) = [Qzy - Qyz, Qxz - Qzx, Qyx - Qxy, Qxx + Qyy + Qzz]/3.0_R8Ki

   lwork = 4*n

   call dgeev('N', 'V', n, A, n, wr, wi, vl, n, vr, n, work, lwork, info)

   ! If error calculating eigenvector/eigenvalues
   if (info /= 0) then
      q = 0.0_R8Ki
      return
   end if

   ! Get index of maximum real eigenvalue
   i = maxloc(wr, dim=1)

   ! Canonical form of quaternion
   q = quat_canonical(vr(4, i), vr(1:3, i))
end function

! quat_to_dcm returns a dcm based on the quaternion where q is a unit quaternion with a positive scalar component
! https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation#Quaternion-derived_rotation_matrix
pure function quat_to_dcm(q) result(dcm)
   real(R8Ki), intent(in)  :: q(3)
   real(R8Ki)              :: dcm(3, 3)
   real(R8Ki)              :: w, ww, xx, yy, zz, n, s
   real(R8Ki)              :: xy, yz, xz, wx, wy, wz

   ! Calculate scalar component
   w = quat_scalar(q)

   ww = w*w
   xx = q(1)*q(1)
   yy = q(2)*q(2)
   zz = q(3)*q(3)

   xy = q(1)*q(2)
   yz = q(2)*q(3)
   xz = q(1)*q(3)

   wx = q(1)*w
   wy = q(2)*w
   wz = q(3)*w

   n = ww + xx + yy + zz
   if (n < epsilon(n)) then
      s = 0.0_R8Ki
   else
      s = 2.0_R8Ki/n
   end if

   dcm(1, 1) = 1.0_R8Ki - s*(yy + zz)
   dcm(2, 1) = s*(xy + wz)
   dcm(3, 1) = s*(xz - wy)

   dcm(1, 2) = s*(xy - wz)
   dcm(2, 2) = 1.0_R8Ki - s*(xx + zz)
   dcm(3, 2) = s*(yz + wx)

   dcm(1, 3) = s*(xz + wy)
   dcm(2, 3) = s*(yz - wx)
   dcm(3, 3) = 1.0_R8Ki - s*(xx + yy)

end function

pure function quat_compose(q1, q2) result(q)
   real(R8Ki), intent(in)  :: q1(3), q2(3)
   real(R8Ki)              :: q(3), q0
   real(R8Ki)              :: w1, x1, y1, z1
   real(R8Ki)              :: w2, x2, y2, z2
   w1 = quat_scalar(q1)
   x1 = q1(1); y1 = q1(2); z1 = q1(3)
   w2 = quat_scalar(q2)
   x2 = q2(1); y2 = q2(2); z2 = q2(3)
   q0 = w1*w2 - x1*x2 - y1*y2 - z1*z2
   q(1) = w1*x2 + x1*w2 + y1*z2 - z1*y2
   q(2) = w1*y2 - x1*z2 + y1*w2 + z1*x2
   q(3) = w1*z2 + x1*y2 - y1*x2 + z1*w2
   q = quat_canonical(q0, q)
end function

pure function quat_inv(q) result(qi)
   real(R8Ki), intent(in)  :: q(3)
   real(R8Ki)              :: qi(3)
   qi = -q
end function

! https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation#Recovering_the_axis-angle_representation
pure function quat_to_rvec(q) result(rvec)
   real(R8Ki), intent(in)  :: q(3)
   real(R8Ki)              :: qr, theta, tmp, rvec(3), m

   ! Magnitude of imaginary part
   m = sqrt(dot_product(q, q))

   ! If this is an identity quaternion, qr == 1, rotation vector is zero
   if (m < epsilon(m)) then
      rvec = 0.0_R8Ki
   else
      qr = sqrt(1.0_R8Ki - m*m)        ! Scalar part
      theta = 2.0_R8Ki*atan2(m, qr)  ! Angle
      rvec = -theta*q/m             ! Negative sign doesn't make sense, but needed for quaternions
   end if
end function

pure function rvec_to_quat(rvec) result(q)
   real(R8Ki), intent(in)  :: rvec(3)
   real(R8Ki)              :: theta, half_theta, q0, q(3)
   theta = sqrt(dot_product(rvec, rvec))
   if (theta < epsilon(theta)) then
      ! Angle is zero, quaternion is identity
      q = 0.0_R8Ki
   else
      half_theta = theta/2.0_R8Ki
      q0 = cos(half_theta)
      q = rvec/theta*sin(half_theta)
      q = -quat_canonical(q0, q) ! Negative sign doesn't make sense, but needed for quaternions
   end if
end function

pure function wm_to_quat(c) result(q)
   real(R8Ki), intent(in)  :: c(3)
   real(R8Ki)              :: c0, q0, q(3)
   c0 = 2.0_R8Ki - dot_product(c, c)/8.0_R8Ki
   q0 = c0/(4.0_R8Ki - c0)
   q = c/(4.0_R8Ki - c0)
   q = quat_canonical(q0, q)
end function

pure function quat_to_wm(q) result(c)
   real(R8Ki), intent(in)  :: q(3)
   real(R8Ki)              :: c(3)
   real(R8Ki)              :: q0
   q0 = quat_scalar(q)
   c = 4.0_R8Ki*q/(1.0_R8Ki + q0)
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

!-------------------------------------------------------------------------------
! Debugging
!-------------------------------------------------------------------------------

subroutine DumpMatrix(unit, filename, A, ErrStat, ErrMsg)
   integer(IntKi), intent(in)             :: unit
   character(*), intent(in)               :: filename
   real(R8Ki), intent(in)                 :: A(:, :)
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg

   character(*), parameter                :: RoutineName = 'DumpMatrix'
   integer(IntKi)                         :: ErrStat2
   character(ErrMsgLen)                   :: ErrMsg2

   ErrStat = ErrID_None
   ErrMsg = ''

   call OpenBOutFile(unit, filename, ErrStat2, ErrMsg2)
   write (unit) int(shape(A), B4Ki)
   write (unit) pack(A, .true.)
   close (unit)
end subroutine

end module
