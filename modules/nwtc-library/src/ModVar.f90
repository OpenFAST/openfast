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
   AngularFields(*) = [VF_AngularDisp, VF_AngularVel, VF_AngularAcc], &
   MotionFields(*) = [VF_TransDisp, VF_Orientation, VF_TransVel, VF_AngularVel, VF_TransAcc, VF_AngularAcc]

public :: MV_InitVars, MV_LinkOutputInput, MV_VarIndex, MV_PackMesh, MV_UnpackMesh
public :: MV_ComputeCentralDiff, MV_Perturb
public :: MV_AddVar, MV_AddMeshVar, MV_AddModule
public :: LoadFields, MotionFields, TransFields, AngularFields
public :: wm_to_dcm, MV_CollectGlobalIndices, wm_compose

contains

subroutine MV_Perturb(Var, iLin, PerturbSign, BaseArr, PerturbArr, iPerturb)
   type(ModVarType), intent(in)           :: Var
   integer(IntKi), intent(in)       :: iLin
   integer(IntKi), intent(in)       :: PerturbSign
   real(R8Ki), intent(in)           :: BaseArr(:)
   real(R8Ki), intent(inout)        :: PerturbArr(:)
   integer(IntKi), intent(out)      :: iPerturb
   real(R8Ki)                       :: Perturb
   real(R8Ki)                       :: WM(3), WMp(3)
   integer(IntKi)                   :: i, j, iLocArr(3)

   ! Copy base array to perturbed array
   PerturbArr = BaseArr

   ! Get variable perturbation and combine with sign
   Perturb = Var%Perturb*real(PerturbSign, R8Ki)

   ! Perturbation index within array
   iPerturb = Var%iLoc(iLin)

   ! If variable is in a mesh and field is orientation
   if (Var%Field == VF_Orientation) then
      j = mod(iLin - 1, 3)                         ! component being modified (0, 1, 2)
      i = iLin - j                                 ! index of start of WM parameters (3)
      iLocArr = Var%iLoc(i:i + 2)                  ! array index vector
      WMp = 0.0_R8Ki                               ! Init WM perturbation to zero
      WMp(j + 1) = 4*tan(Perturb/4.0_R8Ki)         ! WM perturbation around X,Y,Z axis
      WM = PerturbArr(iLocArr)                     ! Current WM parameters value
      PerturbArr(iLocArr) = wm_compose(WM, WMp)    ! Compose value and perturbation
   else
      PerturbArr(Var%iLoc(iLin)) = PerturbArr(Var%iLoc(iLin)) + Perturb
   end if

end subroutine

subroutine MV_InitVars(Vars, Vals, ErrStat, ErrMsg)
   type(ModVarsType), intent(inout)          :: Vars
   type(ModValsType), pointer, intent(inout) :: Vals
   integer(IntKi), intent(out)               :: ErrStat
   character(ErrMsgLen), intent(out)         :: ErrMsg

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
      call ModVarType_Init(Vars%x(i), StartIndex, ErrStat2, ErrMsg2); if (Failed()) return
   end do

   ! Initialize input variables
   StartIndex = 1
   do i = 1, size(Vars%u)
      call ModVarType_Init(Vars%u(i), StartIndex, ErrStat2, ErrMsg2); if (Failed()) return
   end do

   ! Initialize output variables
   StartIndex = 1
   do i = 1, size(Vars%y)
      call ModVarType_Init(Vars%y(i), StartIndex, ErrStat2, ErrMsg2); if (Failed()) return
   end do

   ! Calculate number of variables in group (exclude non linearization vars)
   Vars%Nx = sum(Vars%x%Size, iand(Vars%x%Flags, VF_NoLin) == 0)
   Vars%Nu = sum(Vars%u%Size, iand(Vars%u%Flags, VF_NoLin) == 0)
   Vars%Ny = sum(Vars%y%Size, iand(Vars%y%Flags, VF_NoLin) == 0)

   call AllocAry(Vars%ixg, Vars%Nx, "Vars%ixg", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(Vars%iug, Vars%Nu, "Vars%iug", ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(Vars%iyg, Vars%Ny, "Vars%iyg", ErrStat2, ErrMsg2); if (Failed()) return

   ! Allocate Vals derived type
   allocate (Vals, stat=ErrStat2); if (FailedAlloc()) return

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

subroutine MV_CollectGlobalIndices(VarArr, iGbl)
   type(ModVarType), intent(in)  :: VarArr(:)
   integer(IntKi), intent(out)   :: iGbl(:)
   integer(IntKi)                :: i
   do i = 1, size(VarArr)
      iGbl(VarArr(i)%iLoc) = VarArr(i)%iGbl
   end do
end subroutine

elemental function IsMesh(Var) result(r)
   type(ModVarType), intent(in)  :: Var
   logical                 :: r
   r = iand(Var%Flags, VF_Mesh) > 0
end function

subroutine ModVarType_Init(Var, Index, ErrStat, ErrMsg)
   type(ModVarType), intent(inout)           :: Var
   integer(IntKi), intent(inout)       :: Index
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

   !----------------------------------------------------------------------------
   ! No Linearization
   !----------------------------------------------------------------------------

   if (iand(Var%Flags, VF_NoLin) > 0) then   ! No Linearization

      ! Number of linearization values is zero if NoLin flag is set
      Var%NumLin = 0

   else

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

   ! Initialize global index
   call AllocAry(Var%iGbl, Var%Size, "Var%iGbl", ErrStat2, ErrMsg2); if (Failed()) return
   Var%iGbl = 0

   ! Update index based on variable size
   index = index + Var%Size

   !----------------------------------------------------------------------------
   ! Derivative Order
   !----------------------------------------------------------------------------

   select case (Var%Field)
   case (VF_Orientation, VF_TransDisp, VF_AngularDisp)   ! Position
      Var%DerivOrder = 0
   case (VF_TransVel, VF_AngularVel)                     ! Velocity
      Var%DerivOrder = 1
   case (VF_TransAcc, VF_AngularAcc)                     ! Acceleration
      Var%DerivOrder = 2
   case default
      Var%DerivOrder = -1
   end select

   !----------------------------------------------------------------------------
   ! Module Index
   !----------------------------------------------------------------------------

   ! If module index has been allocated and size does not mach variable size, return error
   if (allocated(Var%iUsr)) then
      if (size(Var%iUsr) < Var%Size) then
         call SetErrStat(ErrID_Fatal, "insufficient iMod given for "//Var%Name, ErrStat, ErrMsg, RoutineName)
         return
      end if
   end if

contains
   function Failed()
      logical :: Failed
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine MV_PackMesh(VarArr, Mesh, arr)
   type(ModVarType), intent(in)        :: VarArr(:)
   type(MeshType), intent(in)    :: Mesh
   real(R8Ki), intent(inout)     :: arr(:)
   integer(IntKi)                :: i, j
   do i = 1, size(VarArr)
      select case (VarArr(i)%Field)
      case (VF_Force)
         arr(VarArr(i)%iLoc) = pack(Mesh%Force, .true.)
      case (VF_Moment)
         arr(VarArr(i)%iLoc) = pack(Mesh%Moment, .true.)
      case (VF_TransDisp)
         arr(VarArr(i)%iLoc) = pack(Mesh%TranslationDisp, .true.)
      case (VF_Orientation)
         do j = 1, VarArr(i)%Nodes
            arr(VarArr(i)%iLoc(3*(j - 1) + 1:3*j)) = wm_from_dcm(Mesh%Orientation(:, :, j))
         end do
      case (VF_TransVel)
         arr(VarArr(i)%iLoc) = pack(Mesh%TranslationVel, .true.)
      case (VF_AngularVel)
         arr(VarArr(i)%iLoc) = pack(Mesh%RotationVel, .true.)
      case (VF_TransAcc)
         arr(VarArr(i)%iLoc) = pack(Mesh%TranslationAcc, .true.)
      case (VF_AngularAcc)
         arr(VarArr(i)%iLoc) = pack(Mesh%RotationAcc, .true.)
      case (VF_Scalar)
         arr(VarArr(i)%iLoc) = pack(Mesh%Scalars, .true.)
      end select
   end do
end subroutine

subroutine MV_UnpackMesh(VarArr, arr, Mesh)
   type(ModVarType), intent(in)        :: VarArr(:)
   real(R8Ki), intent(in)        :: arr(:)
   type(MeshType), intent(inout) :: Mesh
   integer(IntKi)                :: i, j
   do i = 1, size(VarArr)
      select case (VarArr(i)%Field)
      case (VF_Force)
         Mesh%Force = reshape(arr(VarArr(i)%iLoc), shape(Mesh%Force))
      case (VF_Moment)
         Mesh%Moment = reshape(arr(VarArr(i)%iLoc), shape(Mesh%Moment))
      case (VF_TransDisp)
         Mesh%TranslationDisp = reshape(arr(VarArr(i)%iLoc), shape(Mesh%TranslationDisp))
      case (VF_Orientation)
         do j = 1, VarArr(i)%Nodes
            Mesh%Orientation(:, :, j) = wm_to_dcm(arr(VarArr(i)%iLoc(3*(j - 1) + 1:3*j)))
         end do
      case (VF_TransVel)
         Mesh%TranslationVel = reshape(arr(VarArr(i)%iLoc), shape(Mesh%TranslationVel))
      case (VF_AngularVel)
         Mesh%RotationVel = reshape(arr(VarArr(i)%iLoc), shape(Mesh%RotationVel))
      case (VF_TransAcc)
         Mesh%TranslationAcc = reshape(arr(VarArr(i)%iLoc), shape(Mesh%TranslationAcc))
      case (VF_AngularAcc)
         Mesh%RotationAcc = reshape(arr(VarArr(i)%iLoc), shape(Mesh%RotationAcc))
      case (VF_Scalar)
         Mesh%Scalars = reshape(arr(VarArr(i)%iLoc), shape(Mesh%Scalars))
      end select
   end do
end subroutine

subroutine MV_ComputeCentralDiff(VarArr, Delta, PosArr, NegArr, DerivArr)
   type(ModVarType), intent(in)     :: VarArr(:)   ! Array of variables
   real(R8Ki), intent(in)     :: Delta       ! Positive perturbation value
   real(R8Ki), intent(in)     :: PosArr(:)   ! Positive perturbation result array
   real(R8Ki), intent(in)     :: NegArr(:)   ! Negative perturbation result array
   real(R8Ki), intent(inout)  :: DerivArr(:) ! Array containing derivative
   integer(IntKi)             :: i, j, rloc(3)
   real(R8Ki)                 :: WMp(3), WMn(3)

   ! Compute difference between all values
   DerivArr = PosArr - NegArr

   ! Loop through variables
   do i = 1, size(VarArr)

      ! If variable is mesh rotation
      if (VarArr(i)%Field == VF_Orientation) then

         ! Loop through nodes
         do j = 1, VarArr(i)%Nodes

            ! Get vector of indicies of WM rotation parameters in array
            rloc = VarArr(i)%iLoc(3*(j - 1) + 1:3*j)

            ! Get rotation from positive and negative perturbation
            WMp = PosArr(rloc)
            WMn = NegArr(rloc)

            ! Calculate change in rotation and add to array
            DerivArr(rloc) = wm_compose(wm_inv(WMn), WMp)
            ! arrDelta(rloc) = wm_to_zyx(wm_compose(wm_inv(WMn), WMp))
         end do
      end if
   end do

   ! Divide array by 2*delta
   DerivArr = DerivArr/(2.0_R8Ki*Delta)

end subroutine

!-------------------------------------------------------------------------------
! Functions for adding Variables an Modules
!-------------------------------------------------------------------------------

subroutine MV_AddModule(ModArr, ModID, ModAbbr, Instance, DT, Vars, Vals)
   type(ModDataType), allocatable, intent(inout)    :: ModArr(:)
   integer(IntKi), intent(in)                       :: ModID
   character(*), intent(in)                         :: ModAbbr
   integer(IntKi), intent(in)                       :: Instance
   real(R8Ki), intent(in)                           :: DT
   type(ModVarsType), pointer, intent(in)           :: Vars
   type(ModValsType), pointer, intent(in)           :: Vals
   type(ModDataType)                                :: MData
   MData = ModDataType(ID=ModID, Abbr=ModAbbr, Instance=Instance, DT=DT, Vars=Vars, Vals=Vals)
   if (allocated(ModArr)) then
      ModArr = [ModArr, MData]
   else
      ModArr = [MData]
   end if
end subroutine

subroutine MV_AddMeshVar(VarArr, Name, Fields, Nodes, Flags, Perturbs, Active)
   type(ModVarType), allocatable, intent(inout)       :: VarArr(:)
   character(*), intent(in)                     :: Name
   integer(IntKi), intent(in)                   :: Fields(:)
   integer(IntKi), optional, intent(in)         :: Nodes, Flags
   real(R8Ki), optional, intent(in)             :: Perturbs(:)
   logical, optional, intent(in)                :: Active
   integer(IntKi)                               :: i, NodesLocal, FlagsLocal
   logical                                      :: ActiveLocal
   real(R8Ki), allocatable                      :: PerturbsLocal(:)
   NodesLocal = 1
   if (present(Nodes)) NodesLocal = Nodes
   FlagsLocal = 0
   if (present(Flags)) FlagsLocal = Flags
   FlagsLocal = ior(FlagsLocal, VF_Mesh)
   PerturbsLocal = [(0.0_R8Ki, i=1, size(Fields))]
   if (present(Perturbs)) PerturbsLocal = Perturbs
   ActiveLocal = .true.
   if (present(Active)) ActiveLocal = Active
   do i = 1, size(Fields)
      call MV_AddVar(VarArr, Name, Fields(i), Num=NodesLocal, Flags=FlagsLocal, &
                     Perturb=PerturbsLocal(i), Active=ActiveLocal)
   end do
end subroutine

subroutine MV_AddVar(VarArr, Name, Field, Num, Flags, iUsr, Perturb, LinNames, Active)
   type(ModVarType), allocatable, intent(inout)       :: VarArr(:)
   character(*), intent(in)                     :: Name
   integer(IntKi), intent(in)                   :: Field
   integer(IntKi), optional, intent(in)         :: Num, Flags, iUsr(:)
   real(R8Ki), optional, intent(in)             :: Perturb
   logical, optional, intent(in)                :: Active
   character(*), optional, intent(in)           :: LinNames(:)
   integer(IntKi)                               :: i
   type(ModVarType)                                   :: Var

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

   if (allocated(VarArr)) then
      VarArr = [VarArr, Var]
   else
      VarArr = [Var]
   end if
end subroutine

function MV_VarIndex(VarArr, Name, Field) result(Indx)
   type(ModVarType), intent(in)        :: VarArr(:)
   character(*), intent(in)      :: Name
   integer(IntKi), intent(in)    :: Field
   integer(IntKi)                :: Indx
   do Indx = 1, size(VarArr)
      if (string_equal_ci(VarArr(Indx)%Name, Name) .and. &
          VarArr(Indx)%Field == Field) exit
   end do
   if (Indx > size(VarArr)) Indx = 0
end function

!-------------------------------------------------------------------------------
! Functions for linking variables (Output and Input)
!-------------------------------------------------------------------------------

subroutine MV_LinkOutputInput(OutVars, InpVars, OutName, InpName, Field, ErrStat, ErrMsg)
   type(ModVarsType), intent(inout)        :: OutVars, InpVars
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

pure function wm_to_dcm(c) result(R)
   real(R8Ki), intent(in)  :: c(3)
   real(R8Ki)              :: R(3, 3)
   R = quat_to_dcm(wm_to_quat(c))
end function

pure function wm_from_dcm(R) result(c)
   real(R8Ki), intent(in)  :: R(3, 3)
   real(R8Ki)              :: c(3)
   c = wm_from_quat(dcm_to_quat(R))
end function

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

pure function wm_inv(c) result(cinv)
   real(R8Ki), intent(in)  :: c(3)
   real(R8Ki)              :: cinv(3)
   cinv = -c
end function

pure function cross(a, b) result(c)
   real(R8Ki), intent(in) :: a(3), b(3)
   real(R8Ki)             :: c(3)
   c = [a(2)*b(3) - a(3)*b(2), &
        -a(3)*b(1) + a(1)*b(3), &
        a(1)*b(2) - a(2)*b(1)]
end function

end module
