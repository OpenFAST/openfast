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
module FAST_ModData

use NWTC_Library
use NWTC_LAPACK

use FAST_Types

implicit none

private
public :: ModD_AddModule
public :: ModD_GetValLoc, GetModuleOrder
public :: ModD_PackAry, ModD_PackMatrix, ModD_CombineModules

contains

subroutine GetModuleOrder(Mods, ModIDs, ModOrder)
   type(ModDataType), intent(in)             :: Mods(:)     !< Array of module data structures
   integer(IntKi), intent(in)                :: ModIDs(:)   !< List of module IDs to keep in order
   integer(IntKi), allocatable, intent(out)  :: ModOrder(:) !< Module data indices in order of ModIDs
   integer(IntKi), allocatable               :: ModIDAry(:), indices(:)
   integer(IntKi)                            :: i

   ! Create array 1 to size(Mod) representing the index of each module data
   indices = [(i, i=1, size(Mods))]

   ! Get array of module IDs from array of module data
   ModIDAry = [(Mods(i)%ID, i=1, size(Mods))]

   ! Initialize module order array with no size
   allocate (ModOrder(0))

   ! Loop through module IDs to keep, add module data indices that match module ID to order array
   do i = 1, size(ModIDs)
      ModOrder = [ModOrder, pack(indices, ModIDAry == ModIDs(i))]
   end do

end subroutine

subroutine ModD_CombineModules(ModAry, iModOrder, FlagFilter, Linearize, ModOut, ErrStat, ErrMsg)
   type(ModDataType), intent(inout)    :: ModAry(:)
   integer(IntKi), intent(in)          :: iModOrder(:)
   integer(IntKi), intent(in)          :: FlagFilter
   type(ModDataType), intent(out)      :: ModOut
   logical, intent(in)                 :: Linearize
   integer(IntKi), intent(out)         :: ErrStat
   character(ErrMsgLen), intent(out)   :: ErrMsg

   character(*), parameter             :: RoutineName = 'ModD_Build'
   integer(IntKi)                      :: ErrStat2
   character(ErrMsgLen)                :: ErrMsg2
   integer(IntKi)                      :: NumVars
   integer(IntKi)                      :: iGbl(2)
   integer(IntKi)                      :: i, j, k
   integer(IntKi)                      :: iMod
   integer(IntKi)                      :: xNum, xdNum, zNum, uNum, yNum
   integer(IntKi)                      :: ix, ixd, iz, iu, iy
   character(20)                       :: NamePrefix

   ! Initialize error return
   ErrStat = ErrID_None
   ErrMsg = ""

   ! If no modules or order is empty, return error
   if ((size(ModAry) == 0) .or. (size(iModOrder) == 0)) then
      call SetErrStat(ErrID_Fatal, "No modules were used", ErrStat, ErrMsg, RoutineName)
      return
   end if

   !----------------------------------------------------------------------------
   ! Construct index to lookup variables
   !----------------------------------------------------------------------------

   ! Allocate variable index array with size equal to number of modules
   allocate (ModOut%Xfr(size(ModAry)), stat=ErrStat2)
   if (FailedAlloc("ModOut%Xfr")) return

   !----------------------------------------------------------------------------
   ! Combine modules into output module
   !----------------------------------------------------------------------------

   ! Clear module linearization abbreviation
   ModOut%Lin%Abbr = ""

   ! Allocate variable structure for glue
   allocate (ModOut%Vars)

   ! Initialize number of variables in each group
   xNum = 0; xdNum = 0; zNum = 0; uNum = 0; yNum = 0

   ! Loop through each module and sum the number of variables that will be in
   ! the combined module
   do i = 1, size(iModOrder)
      iMod = iModOrder(i)
      associate (ModData => ModAry(iMod))

         ! Continuous state
         call CountVariablesFiltered(ModData%Vars%x, NumVars)
         allocate (ModOut%Xfr(iMod)%x(NumVars), stat=ErrStat2)
         if (FailedAlloc("ModOut%Xfr(iMod)%x")) return
         xNum = xNum + NumVars

         ! Discrete state
         call CountVariablesFiltered(ModData%Vars%xd, NumVars)
         allocate (ModOut%Xfr(iMod)%xd(NumVars), stat=ErrStat2)
         if (FailedAlloc("ModOut%Xfr(iMod)%xd")) return
         xdNum = xdNum + NumVars

         ! Constraint state
         call CountVariablesFiltered(ModData%Vars%z, NumVars)
         allocate (ModOut%Xfr(iMod)%z(NumVars), stat=ErrStat2)
         if (FailedAlloc("ModOut%Xfr(iMod)%z")) return
         zNum = zNum + NumVars

         ! Input
         call CountVariablesFiltered(ModData%Vars%u, NumVars)
         allocate (ModOut%Xfr(iMod)%u(NumVars), stat=ErrStat2)
         if (FailedAlloc("ModOut%Xfr(iMod)%u")) return
         uNum = uNum + NumVars

         ! Output
         call CountVariablesFiltered(ModData%Vars%y, NumVars)
         allocate (ModOut%Xfr(iMod)%y(NumVars), stat=ErrStat2)
         if (FailedAlloc("ModOut%Xfr(iMod)%y")) return
         yNum = yNum + NumVars

      end associate
   end do

   ! Allocate arrays for to hold combined variables
   allocate (ModOut%Vars%x(xNum), stat=ErrStat2); if (FailedAlloc("ModOut%Vars%x")) return
   allocate (ModOut%Vars%xd(xdNum), stat=ErrStat2); if (FailedAlloc("ModOut%Vars%xd")) return
   allocate (ModOut%Vars%z(zNum), stat=ErrStat2); if (FailedAlloc("ModOut%Vars%z")) return
   allocate (ModOut%Vars%u(uNum), stat=ErrStat2); if (FailedAlloc("ModOut%Vars%u")) return
   allocate (ModOut%Vars%y(yNum), stat=ErrStat2); if (FailedAlloc("ModOut%Vars%y")) return

   ! Initialize variable index counters
   ix = 1; ixd = 1; iz = 1; iu = 1; iy = 1

   ! Initialize number of values in each group variable group
   ModOut%Vars%Nx = 0
   ModOut%Vars%Nxd = 0
   ModOut%Vars%Nz = 0
   ModOut%Vars%Nu = 0
   ModOut%Vars%Ny = 0

   ! Loop through each module by index and add variables
   do i = 1, size(iModOrder)
      iMod = iModOrder(i)
      associate (ModData => ModAry(iMod))

         ! Create variable name prefix for linearization names. Add instance
         ! number to module abbreviation if more than 1 instance or the module is BeamDyn
         if ((ModData%ID == Module_BD) .or. (count(ModAry%ID == ModData%ID) > 1)) then
            NamePrefix = trim(ModData%Abbr)//"_"//Num2LStr(ModData%Ins)
            ModData%Lin%Abbr = "."//trim(ModData%Abbr)//Num2LStr(ModData%Ins)
         else
            NamePrefix = ModData%Abbr
            ModData%Lin%Abbr = "."//ModData%Abbr
         end if

         if (size(ModData%Vars%x) > 0) call AddVariables(ModData%Vars%x, ModOut%Vars%x, ModOut%Xfr(iMod)%x, ix, ModOut%Vars%Nx)       ! Continuous state
         if (size(ModData%Vars%xd) > 0) call AddVariables(ModData%Vars%xd, ModOut%Vars%xd, ModOut%Xfr(iMod)%xd, ixd, ModOut%Vars%Nxd)  ! Discrete state
         if (size(ModData%Vars%z) > 0) call AddVariables(ModData%Vars%z, ModOut%Vars%z, ModOut%Xfr(iMod)%z, iz, ModOut%Vars%Nz)       ! Constraint state
         if (size(ModData%Vars%u) > 0) call AddVariables(ModData%Vars%u, ModOut%Vars%u, ModOut%Xfr(iMod)%u, iu, ModOut%Vars%Nu)       ! Input
         if (size(ModData%Vars%y) > 0) call AddVariables(ModData%Vars%y, ModOut%Vars%y, ModOut%Xfr(iMod)%y, iy, ModOut%Vars%Ny)       ! Output
      end associate
   end do

contains

   subroutine AddVariables(VarAryIn, VarAryOut, VarXfr, iVar, iVal)
      type(ModVarType), intent(in)     :: VarAryIn(:)
      type(ModVarType), intent(inout)  :: VarAryOut(:)
      type(VarXfrType), intent(inout)  :: VarXfr(:)
      integer(IntKi), intent(inout)    :: iVar
      integer(IntKi), intent(inout)    :: iVal

      integer(IntKi)                   :: NumVals, iXfr

      iXfr = 1

      ! Loop through variables in original module
      do k = 1, size(VarAryIn)

         ! If filter flag is not none and variable doesn't have flag, cycle
         if (.not. MV_HasFlags(VarAryIn(k), FlagFilter) .and. FlagFilter /= VF_None) cycle

         associate (Var => VarAryOut(iVar))

            ! Add variable to module
            VarAryOut(iVar) = VarAryIn(k)

            ! Get number of values in variable
            NumVals = VarAryIn(k)%Num

            ! Set value indices in combined module
            Var%iLoc = [iVal + 1, iVal + NumVals]

            ! Increment global value index
            iVal = iVal + NumVals

            ! Set transfer index
            VarXfr(iXfr)%iVar = k                ! Variable number in source module
            VarXfr(iXfr)%NumVals = NumVals       ! Number of values in variable
            VarXfr(iXfr)%iSrc = VarAryIn(k)%iLoc ! value start-end indices in source module
            VarXfr(iXfr)%iDst = Var%iLoc         ! Value start-end indices in destination module

            ! Increment transfer index
            iXfr = iXfr + 1

            ! Prepend module names
            call AddLinNamePrefix(Var, NamePrefix)

            ! Increment variable index
            iVar = iVar + 1

         end associate

      end do

   end subroutine

   subroutine CountVariablesFiltered(VarAry, nVars)
      type(ModVarType), intent(in)  :: VarAry(:)
      integer(IntKi), intent(out)   :: nVars
      nVars = 0
      ! If no filter
      if (FlagFilter == VF_None) then
         ! Count all variables in array
         nVars = size(VarAry)
      else
         ! Loop through filters and increment nVars if they have the flag
         do k = 1, size(VarAry)
            if (MV_HasFlags(VarAry(k), FlagFilter)) nVars = nVars + 1
         end do
      end if
   end subroutine

   subroutine AddLinNamePrefix(Var, Prefix)
      type(ModVarType), intent(inout)  :: Var
      character(*), intent(in)         :: Prefix
      integer(IntKi)                   :: j
      if (allocated(Var%LinNames)) then
         do j = 1, size(Var%LinNames)
            Var%LinNames(j) = trim(Prefix)//" "//Var%LinNames(j)
         end do
      end if
   end subroutine

   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function

   logical function FailedAlloc(name)
      character(*), intent(in)   :: name
      if (ErrStat2 == 0) then
         FailedAlloc = .false.
      else
         call SetErrStat(ErrID_Fatal, "Failed to allocate "//name, ErrStat, ErrMsg, RoutineName)
         FailedAlloc = .true.
      end if
   end function

end subroutine

! ModD_GetValLoc is used to get the global or module value indices based on module index and variable index.
! iMod is the start and end indices of the values in the module data
! iGbl is teh start and end indices of the values in the global data
logical function ModD_GetValLoc(VarXfrAry, iVar, iSrc, iDst) result(Active)
   type(VarXfrType), intent(in)           :: VarXfrAry(:)
   integer(IntKi), intent(in)             :: iVar
   integer(IntKi), optional, intent(out)  :: iSrc(2), iDst(2)
   integer(IntKi)                         :: i
   do i = 1, size(VarXfrAry)
      if (VarXfrAry(i)%iVar /= iVar) cycle
      if (present(iSrc)) iSrc = VarXfrAry(i)%iSrc
      if (present(iDst)) iDst = VarXfrAry(i)%iDst
      Active = .true.
      return
   end do
   Active = .false.
end function

subroutine ModD_PackAry(VarXfrAry, SrcAry, DstAry)
   type(VarXfrType), intent(in)     :: VarXfrAry(:)
   real(R8Ki), intent(in)           :: SrcAry(:)
   real(R8Ki), intent(inout)        :: DstAry(:)
   integer(IntKi)                   :: i
   do i = 1, size(VarXfrAry)
      DstAry(VarXfrAry(i)%iDst(1):VarXfrAry(i)%iDst(2)) = &
         SrcAry(VarXfrAry(i)%iSrc(1):VarXfrAry(i)%iSrc(2))
   end do
end subroutine

subroutine ModD_PackMatrix(RowXfrAry, ColXfrAry, SrcMat, DstMat)
   type(VarXfrType), intent(in)           :: RowXfrAry(:), ColXfrAry(:)
   real(R8Ki), intent(in)                 :: SrcMat(:, :)
   real(R8Ki), intent(inout)              :: DstMat(:, :)
   integer(IntKi)                         :: i, j
   do i = 1, size(RowXfrAry)
      do j = 1, size(ColXfrAry)
         DstMat(RowXfrAry(i)%iDst(1):RowXfrAry(i)%iDst(2), ColXfrAry(j)%iDst(1):ColXfrAry(j)%iDst(2)) = &
            SrcMat(RowXfrAry(i)%iSrc(1):RowXfrAry(i)%iSrc(2), ColXfrAry(j)%iSrc(1):ColXfrAry(j)%iSrc(2))
      end do
   end do
end subroutine

subroutine ModD_AddModule(Mods, ModID, ModAbbr, Instance, ModDT, SolverDT, Vars, ErrStat, ErrMsg)
   type(ModDataType), allocatable, intent(inout)   :: Mods(:)
   integer(IntKi), intent(in)                      :: ModID
   character(*), intent(in)                        :: ModAbbr
   integer(IntKi), intent(in)                      :: Instance
   real(R8Ki), intent(in)                          :: ModDT
   real(R8Ki), intent(in)                          :: SolverDT
   type(ModVarsType), pointer, intent(in)          :: Vars
   integer(IntKi), intent(out)                     :: ErrStat
   character(ErrMsgLen), intent(out)               :: ErrMsg

   character(*), parameter                         :: RoutineName = 'ModD_AddModule'
   integer(IntKi)                                  :: ErrStat2
   character(ErrMsgLen)                            :: ErrMsg2
   type(ModDataType)                               :: ModData
   integer(IntKi)                                  :: iMod

   ErrStat = ErrID_None
   ErrMsg = ''

   ! If module array hasn't been allocated, allocate with zero size
   if (.not. allocated(Mods)) allocate (Mods(0))

   ! Populate ModuleDataType derived type
   ModData = ModDataType(iMod=size(Mods) + 1, ID=ModID, Abbr=ModAbbr, &
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

   ! Get index of newly added module
   iMod = size(Mods)

   ! Set module index in each variable
   ModData%Vars%x%iMod = iMod
   ModData%Vars%xd%iMod = iMod
   ModData%Vars%z%iMod = iMod
   ModData%Vars%u%iMod = iMod
   ModData%Vars%y%iMod = iMod

end subroutine

end module
