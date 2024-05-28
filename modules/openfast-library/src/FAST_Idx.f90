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
module FAST_Idx

use NWTC_Library
use NWTC_LAPACK

use FAST_Types

implicit none

private
public :: Idx_Init, Idx_GetValLoc, GetModuleOrder

contains

subroutine GetModuleOrder(Mods, ModIDs, ModOrder)
   type(ModDataType), intent(in)             :: Mods(:)     !< Array of module data structures
   integer(IntKi), intent(in)                :: ModIDs(:)   !< List of module IDs to keep in order
   integer(IntKi), allocatable, intent(out)  :: ModOrder(:) !< Module data indices in order of ModIDs
   integer(IntKi), allocatable               :: ModIDAry(:), indices(:)
   integer(IntKi)                            :: i

   ! Create array 1 to size(Mod) representing the index of each module data
   indices = [(i, i = 1, size(Mods))]

   ! Get array of module IDs from array of module data
   ModIDAry = [(Mods(i)%ID, i = 1, size(Mods))]

   ! Initialize module order array with no size
   allocate (ModOrder(0))

   ! Loop through module IDs to keep, add module data indices that match module ID to order array
   do i = 1, size(ModIDs)
      ModOrder = [ModOrder, pack(indices, ModIDAry == ModIDs(i))]
   end do

end subroutine

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

! Idx_GetValLoc is used to get the global or module value indices based on module index and variable index.
! iMod is the start and end indices of the values in the module data
! iGbl is teh start and end indices of the values in the global data
subroutine Idx_GetValLoc(Idx, ModIdx, VarIdx, iMod, iGbl)
   type(VarIdxType), intent(in)           :: Idx
   integer(IntKi), intent(in)             :: ModIdx, VarIdx
   integer(IntKi), optional, intent(out)  :: iMod(2), iGbl(2)
   integer(IntKi)                         :: col
   col = Idx%ModVarStart(ModIdx) + VarIdx
   if (present(iMod)) iMod = Idx%ValLocGbl(1:2, col)
   if (present(iGbl)) iGbl = Idx%ValLocGbl(3:4, col)
end subroutine

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
