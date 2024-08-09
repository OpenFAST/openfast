!**********************************************************************************************************************************
! FAST_ModGlue.f90 performs linearization using the ModVars module.
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
public :: ModGlue_Init
public :: ModGlue_Linearize_OP, ModGlue_CalcSteady
public :: ModGlue_SaveOperatingPoint, ModGlue_RestoreOperatingPoint
public :: CalcWriteLinearMatrices, Glue_CombineModules

contains

subroutine Glue_CombineModules(ModGlue, ModDataAry, Mappings, iModAry, FlagFilter, Linearize, ErrStat, ErrMsg)
   type(ModGlueType), intent(out)      :: ModGlue
   type(ModDataType), intent(in)       :: ModDataAry(:)
   integer(IntKi), intent(in)          :: iModAry(:)
   integer(IntKi), intent(in)          :: FlagFilter
   logical, intent(in)                 :: Linearize
   type(MappingType), intent(in)       :: Mappings(:)    !< Mesh and variable mappings
   integer(IntKi), intent(out)         :: ErrStat
   character(ErrMsgLen), intent(out)   :: ErrMsg

   character(*), parameter             :: RoutineName = 'Glue_CombineModules'
   integer(IntKi)                      :: ErrStat2
   character(ErrMsgLen)                :: ErrMsg2
   integer(IntKi)                      :: iGbl(2)
   integer(IntKi)                      :: i, j, k
   integer(IntKi)                      :: iMod, iVarGlue
   integer(IntKi)                      :: xNumVals, zNumVals, uNumVals, yNumVals
   integer(IntKi)                      :: xNumVars, zNumVars, uNumVars, yNumVars
   integer(IntKi)                      :: ix, iz, iu, iy
   character(20)                       :: NamePrefix
   type(VarMapType)                    :: ModMap

   ! Initialize error return
   ErrStat = ErrID_None
   ErrMsg = ""

   ! If no modules or order is empty, return error
   if ((size(ModDataAry) == 0) .or. (size(iModAry) == 0)) then
      call SetErrStat(ErrID_Fatal, "No modules were used", ErrStat, ErrMsg, RoutineName)
      return
   end if

   !----------------------------------------------------------------------------
   ! Allocate module data array
   !----------------------------------------------------------------------------

   ! Allocate module info array based on number of modules in iMod
   allocate (ModGlue%ModData(size(iModAry)), stat=ErrStat2)
   if (FailedAlloc("ModOut%VarsAry")) return

   !----------------------------------------------------------------------------
   ! Combine modules into glue module
   !----------------------------------------------------------------------------

   ! Initialize number of variables and values in each group
   xNumVars = 0; zNumVars = 0; uNumVars = 0; yNumVars = 0
   xNumVals = 0; zNumVals = 0; uNumVals = 0; yNumVals = 0

   ! Loop through each module and sum the number of variables that will be in
   ! the combined module
   do i = 1, size(iModAry)
      associate (ModData => ModDataAry(iModAry(i)), GlueModData => ModGlue%ModData(i))

         ! Copy values from source module info
         GlueModData%Abbr = ModData%Abbr
         GlueModData%ID = ModData%ID
         GlueModData%iMod = ModData%iMod  ! Keep original module index for input solve
         GlueModData%Ins = ModData%Ins
         GlueModData%DT = ModData%DT
         GlueModData%SubSteps = ModData%SubSteps

         ! Continuous state
         call CopyVariables(ModData%Vars%x, GlueModData%Vars%x, xNumVals); if (Failed()) return
         GlueModData%Vars%Nx = ModData%Vars%Nx ! Same as original module
         xNumVars = xNumVars + size(GlueModData%Vars%x)

         ! Constraint state
         call CopyVariables(ModData%Vars%z, GlueModData%Vars%z, zNumVals); if (Failed()) return
         GlueModData%Vars%Nz = ModData%Vars%Nz ! Same as original module
         zNumVars = zNumVars + size(GlueModData%Vars%z)

         ! Input
         call CopyVariables(ModData%Vars%u, GlueModData%Vars%u, uNumVals); if (Failed()) return
         GlueModData%Vars%Nu = ModData%Vars%Nu ! Same as original module
         uNumVars = uNumVars + size(GlueModData%Vars%u)

         ! Output
         call CopyVariables(ModData%Vars%y, GlueModData%Vars%y, yNumVals); if (Failed()) return
         GlueModData%Vars%Ny = ModData%Vars%Ny ! Same as original module
         yNumVars = yNumVars + size(GlueModData%Vars%y)

      end associate
   end do

   ! Set total number of values in glue module
   ModGlue%Vars%Nx = xNumVals
   ModGlue%Vars%Nz = zNumVals
   ModGlue%Vars%Nu = uNumVals
   ModGlue%Vars%Ny = yNumVals

   ! Allocate arrays for to hold combined variables
   allocate (ModGlue%Vars%x(xNumVars), stat=ErrStat2); if (FailedAlloc("ModOut%Vars%x")) return
   allocate (ModGlue%Vars%z(zNumVars), stat=ErrStat2); if (FailedAlloc("ModOut%Vars%z")) return
   allocate (ModGlue%Vars%u(uNumVars), stat=ErrStat2); if (FailedAlloc("ModOut%Vars%u")) return
   allocate (ModGlue%Vars%y(yNumVars), stat=ErrStat2); if (FailedAlloc("ModOut%Vars%y")) return

   ! Loop through module info in glue module
   ix = 0; iz = 0; iu = 0; iy = 0
   do i = 1, size(ModGlue%ModData)

      associate (GlueModData => ModGlue%ModData(i))

         ! Determine module name prefix for linearization
         if ((GlueModData%ID == Module_BD) .or. (count(ModDataAry%ID == GlueModData%ID) > 1)) then
            NamePrefix = trim(GlueModData%Abbr)//"_"//Num2LStr(GlueModData%Ins)
            GlueModData%Abbr = trim(GlueModData%Abbr)//Num2LStr(GlueModData%Ins)
         else
            NamePrefix = GlueModData%Abbr
            GlueModData%Abbr = GlueModData%Abbr
         end if

         ! Continuous state
         do j = 1, size(GlueModData%Vars%x)
            ix = ix + 1
            ModGlue%Vars%x(ix) = GlueModData%Vars%x(j)
            ModGlue%Vars%x(ix)%iLoc = ModGlue%Vars%x(ix)%iGlu  ! Set local indices to glue indices
            ModGlue%Vars%x(ix)%iGlu = 0                        ! Set glue indices to 0
            call AddLinNamePrefix(ModGlue%Vars%x(ix), NamePrefix)
         end do

         ! Constraint state
         do j = 1, size(GlueModData%Vars%z)
            iz = iz + 1
            ModGlue%Vars%z(iz) = GlueModData%Vars%z(j)
            ModGlue%Vars%z(iz)%iLoc = ModGlue%Vars%z(iz)%iGlu  ! Set local indices to glue indices
            ModGlue%Vars%z(iz)%iGlu = 0                        ! Set glue indices to 0
            call AddLinNamePrefix(ModGlue%Vars%z(iz), NamePrefix)
         end do

         ! Input
         do j = 1, size(GlueModData%Vars%u)
            iu = iu + 1
            ModGlue%Vars%u(iu) = GlueModData%Vars%u(j)
            ModGlue%Vars%u(iu)%iLoc = ModGlue%Vars%u(iu)%iGlu  ! Set local indices to glue indices
            ModGlue%Vars%u(iu)%iGlu = 0                        ! Set glue indices to 0
            call AddLinNamePrefix(ModGlue%Vars%u(iu), NamePrefix)
         end do

         ! Output
         do j = 1, size(GlueModData%Vars%y)
            iy = iy + 1
            ModGlue%Vars%y(iy) = GlueModData%Vars%y(j)
            ModGlue%Vars%y(iy)%iLoc = ModGlue%Vars%y(iy)%iGlu  ! Set local indices to glue indices
            ModGlue%Vars%y(iy)%iGlu = 0                        ! Set glue indices to 0
            call AddLinNamePrefix(ModGlue%Vars%y(iy), NamePrefix)
         end do

      end associate
   end do

   !----------------------------------------------------------------------------
   ! Determine mappings which apply to the modules in this glue module
   !----------------------------------------------------------------------------

   allocate (ModGlue%VarMaps(0))

   ! Loop through mappings
   do i = 1, size(Mappings)

      ! Find index of source module in glue module, cycle if not found
      ModMap%iModSrc = 0
      do j = 1, size(iModAry)
         if (iModAry(j) == Mappings(i)%iModSrc) then
            ModMap%iModSrc = j
            exit
         end if
      end do
      if (ModMap%iModSrc == 0) cycle

      ! Find index of destination module in glue module, cycle if not found
      ModMap%iModDst = 0
      do j = 1, size(iModAry)
         if (iModAry(j) == Mappings(i)%iModDst) then
            ModMap%iModDst = j
            exit
         end if
      end do
      if (ModMap%iModDst == 0) cycle

      ! Get source and destination modules from glue module data array
      associate (Mapping => Mappings(i), &
                 ModSrc => ModGlue%ModData(ModMap%iModSrc), &
                 ModDst => ModGlue%ModData(ModMap%iModDst))

         ! Set mapping index and clear variable indices
         ModMap%iMapping = i
         ModMap%iVarSrc = 0
         ModMap%iVarSrcDisp = 0
         ModMap%iVarDst = 0
         ModMap%iVarDstDisp = 0

         ! Init variable indices and find indices that apply to the source data location
         select case (Mapping%MapType)
         case (Map_Variable)

            do j = 1, size(ModSrc%Vars%y)
               if (MV_EqualDL(ModSrc%Vars%y(j)%DL, Mapping%SrcDL)) ModMap%iVarSrc(1) = j
            end do

         case (Map_LoadMesh, Map_MotionMesh)

            do j = 1, size(ModSrc%Vars%y)
               if (MV_EqualDL(ModSrc%Vars%y(j)%DL, Mapping%SrcDL)) ModMap%iVarSrc(ModSrc%Vars%y(j)%Field) = j
            end do

            if (Mapping%MapType == Map_LoadMesh) then
               do j = 1, size(ModSrc%Vars%u)
                  if (MV_EqualDL(ModSrc%Vars%u(j)%DL, Mapping%SrcDispDL)) ModMap%iVarSrcDisp(ModSrc%Vars%u(j)%Field) = j
               end do
            end if

         end select

         ! If no source variable indices found, cycle
         if (all(ModMap%iVarSrc == 0)) cycle
         if (Mapping%MapType == Map_LoadMesh .and. all(ModMap%iVarSrcDisp == 0)) cycle

         ! Init variable indices and find indices that apply to the destination data location
         select case (Mapping%MapType)
         case (Map_Variable)

            do j = 1, size(ModDst%Vars%u)
               if (MV_EqualDL(ModDst%Vars%u(j)%DL, Mapping%DstDL)) ModMap%iVarDst(1) = j
            end do

         case (Map_LoadMesh, Map_MotionMesh)

            do j = 1, size(ModDst%Vars%u)
               if (MV_EqualDL(ModDst%Vars%u(j)%DL, Mapping%DstDL)) ModMap%iVarDst(ModDst%Vars%u(j)%Field) = j
            end do

            if (Mapping%MapType == Map_LoadMesh) then
               do j = 1, size(ModDst%Vars%y)
                  if (MV_EqualDL(ModDst%Vars%y(j)%DL, Mapping%DstDispDL)) ModMap%iVarDstDisp(ModDst%Vars%y(j)%Field) = j
               end do
            end if

         end select

         ! If no destination variable indices found, cycle
         if (all(ModMap%iVarDst == 0)) cycle
         if (Mapping%MapType == Map_LoadMesh .and. all(ModMap%iVarDstDisp == 0)) cycle

         ! Add new module mapping to array
         ModGlue%VarMaps = [ModGlue%VarMaps, ModMap]

      end associate
   end do

   !----------------------------------------------------------------------------
   ! Linearization
   !----------------------------------------------------------------------------

   if (.not. Linearize) return

   ! Allocate linearization arrays
   if (ModGlue%Vars%Nx > 0) then
      call AllocAry(ModGlue%Lin%x, ModGlue%Vars%Nx, "x", ErrStat2, ErrMsg2)
      if (Failed()) return
   end if
   if (ModGlue%Vars%Nx > 0) then
      call AllocAry(ModGlue%Lin%dx, ModGlue%Vars%Nx, "dx", ErrStat2, ErrMsg2)
      if (Failed()) return
   end if
   if (ModGlue%Vars%Nz > 0) then
      call AllocAry(ModGlue%Lin%z, ModGlue%Vars%Nz, "z", ErrStat2, ErrMsg2)
      if (Failed()) return
   end if
   if (ModGlue%Vars%Nu > 0) then
      call AllocAry(ModGlue%Lin%u, ModGlue%Vars%Nu, "u", ErrStat2, ErrMsg2)
      if (Failed()) return
   end if
   if (ModGlue%Vars%Ny > 0) then
      call AllocAry(ModGlue%Lin%y, ModGlue%Vars%Ny, "y", ErrStat2, ErrMsg2)
      if (Failed()) return
   end if

   ! Allocate full Jacobian matrices
   if (ModGlue%Vars%Ny > 0 .and. ModGlue%Vars%Nu > 0) then
      call AllocAry(ModGlue%Lin%dYdu, ModGlue%Vars%Ny, ModGlue%Vars%Nu, "dYdu", ErrStat2, ErrMsg2)
      if (Failed()) return
   end if
   if (ModGlue%Vars%Nx > 0 .and. ModGlue%Vars%Nu > 0) then
      call AllocAry(ModGlue%Lin%dXdu, ModGlue%Vars%Nx, ModGlue%Vars%Nu, "dXdu", ErrStat2, ErrMsg2)
      if (Failed()) return
   end if
   if (ModGlue%Vars%Ny > 0 .and. ModGlue%Vars%Nx > 0) then
      call AllocAry(ModGlue%Lin%dYdx, ModGlue%Vars%Ny, ModGlue%Vars%Nx, "dYdx", ErrStat2, ErrMsg2)
      if (Failed()) return
   end if
   if (ModGlue%Vars%Nx > 0 .and. ModGlue%Vars%Nx > 0) then
      call AllocAry(ModGlue%Lin%dXdx, ModGlue%Vars%Nx, ModGlue%Vars%Nx, "dXdx", ErrStat2, ErrMsg2)
      if (Failed()) return
   end if
   if (ModGlue%Vars%Nu > 0 .and. ModGlue%Vars%Nu > 0) then
      call AllocAry(ModGlue%Lin%dUdu, ModGlue%Vars%Nu, ModGlue%Vars%Nu, "dUdu", ErrStat2, ErrMsg2)
      if (Failed()) return
   end if
   if (ModGlue%Vars%Nu > 0 .and. ModGlue%Vars%Ny > 0) then
      call AllocAry(ModGlue%Lin%dUdy, ModGlue%Vars%Nu, ModGlue%Vars%Ny, "dUdy", ErrStat2, ErrMsg2)
      if (Failed()) return
   end if

contains

   subroutine CopyVariables(VarAryIn, VarAryOut, iVal)
      type(ModVarType), intent(in)                 :: VarAryIn(:)
      type(ModVarType), allocatable, intent(inout) :: VarAryOut(:)
      integer(IntKi), intent(inout)                :: iVal

      integer(IntKi)                   :: NumVars, NumVals, iVar

      ! Get number of variables that have flag
      NumVars = 0
      do k = 1, size(VarAryIn)
         if (MV_HasFlagsAny(VarAryIn(k), FlagFilter)) NumVars = NumVars + 1
      end do

      ! Allocate output array of variables
      allocate (VarAryOut(NumVars), stat=ErrStat2)
      if (ErrStat2 /= 0) then
         ErrStat2 = ErrID_Fatal
         ErrMsg2 = "Failed to allocate vars"
         return
      end if

      iVar = 1

      ! Loop through variables in original module
      do k = 1, size(VarAryIn)

         ! If variable doesn't have flag, cycle
         if (.not. MV_HasFlagsAny(VarAryIn(k), FlagFilter)) cycle

         associate (Var => VarAryOut(iVar))

            ! Copy variable
            Var = VarAryIn(k)

            ! Get number of values in variable
            NumVals = VarAryIn(k)%Num

            ! Set value indices in combined module
            Var%iGlu = [iVal + 1, iVal + NumVals]

            ! Increment global value index
            iVal = iVal + NumVals

            ! Increment variable index in module info variable array
            iVar = iVar + 1

            ! Deallocate linearization names if not doing linearization
            if (.not. Linearize .and. allocated(Var%LinNames)) deallocate (Var%LinNames)

         end associate

      end do

   end subroutine

   subroutine AddLinNamePrefix(Var, Prefix)
      type(ModVarType), intent(inout)  :: Var
      character(*), intent(in)         :: Prefix
      integer(IntKi)                   :: m
      if (allocated(Var%LinNames)) then
         do m = 1, size(Var%LinNames)
            Var%LinNames(m) = trim(Prefix)//" "//Var%LinNames(m)
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

subroutine ModGlue_Init(p, m, y, p_FAST, m_FAST, Turbine, ErrStat, ErrMsg)
   type(Glue_ParameterType), intent(inout)         :: p        !< Glue Parameters
   type(Glue_MiscVarType), intent(inout)           :: m        !< Glue MiscVars
   type(Glue_OutputFileType), intent(inout)        :: y        !< Glue Output
   type(FAST_ParameterType), intent(inout)         :: p_FAST   !< FAST Parameters
   type(FAST_MiscVarType), intent(inout)           :: m_FAST   !< FAST MiscVars
   type(FAST_TurbineType), intent(inout)           :: Turbine
   integer(IntKi), intent(out)                     :: ErrStat
   character(*), intent(out)                       :: ErrMsg

   character(*), parameter                         :: RoutineName = 'ModGlue_Init'
   integer(IntKi)                                  :: ErrStat2
   character(ErrMsgLen)                            :: ErrMsg2
   integer(IntKi), allocatable                     :: modIDs(:), modIdx(:)
   integer(IntKi)                                  :: i, j, k
   integer(IntKi)                                  :: FlagFilters

   ! Initialize error return
   ErrStat = ErrID_None
   ErrMsg = ""

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
   allocate (p%Lin%iMod(0))
   do i = 1, size(LinMods)
      p%Lin%iMod = [p%Lin%iMod, pack(modIdx, ModIDs == LinMods(i))]
   end do

   ! Loop through modules, if module is not in index, return with error
   if (p_FAST%Linearize) then
      do i = 1, size(m%ModData)
         if (.not. any(i == p%Lin%iMod)) then
            call SetErrStat(ErrID_Fatal, "Module "//trim(m%ModData(i)%Abbr)// &
                           " not supported in linearization", ErrStat, ErrMsg, RoutineName)
            return
         end if
      end do
   end if

   !----------------------------------------------------------------------------
   ! Set Variable Flags for linearization
   !----------------------------------------------------------------------------

   ! Loop through each module by index
   do i = 1, size(p%Lin%iMod)
      associate (ModData => m%ModData(p%Lin%iMod(i)))

         ! Set linearize flag on all continuous state variables
         do j = 1, size(ModData%Vars%x)
            call MV_SetFlags(ModData%Vars%x(j), VF_Linearize)
         end do

         ! Add or remove linearize flag based on requested input
         select case (p_FAST%LinInputs)
         case (LIN_NONE)
            do j = 1, size(ModData%Vars%u)
               call MV_ClearFlags(ModData%Vars%u(j), VF_Linearize)
            end do
         case (LIN_STANDARD)
            ! For standard inputs, use VF_Linearize flag as set in the module
         case (LIN_ALL)
            do j = 1, size(ModData%Vars%u)
               call MV_SetFlags(ModData%Vars%u(j), VF_Linearize)
            end do
         end select

         ! Add or remove linearize flag based on requested output
         select case (p_FAST%LinOutputs)
         case (LIN_NONE)
            do j = 1, size(ModData%Vars%y)
               call MV_ClearFlags(ModData%Vars%y(j), VF_Linearize)
            end do
         case (LIN_STANDARD)  ! Set linearize flag for write output variables
            do j = 1, size(ModData%Vars%y)
               if (MV_HasFlagsAll(ModData%Vars%y(j), VF_WriteOut)) then
                  call MV_SetFlags(ModData%Vars%y(j), VF_Linearize)
               else
                  call MV_ClearFlags(ModData%Vars%y(j), VF_Linearize)
               end if
            end do
         case (LIN_ALL)
            do j = 1, size(ModData%Vars%y)
               call MV_SetFlags(ModData%Vars%y(j), VF_Linearize)
            end do
         end select

      end associate
   end do

   !----------------------------------------------------------------------------
   ! Glue Module
   !----------------------------------------------------------------------------

   call Glue_CombineModules(m%ModGlue, m%ModData, m%Mappings, p%Lin%iMod, VF_None, p_FAST%Linearize, ErrStat2, ErrMsg2); if (Failed()) return

   !----------------------------------------------------------------------------
   ! Allocate linearization arrays and matrices
   !----------------------------------------------------------------------------

   ! If linearization is enabled
   if (p_FAST%Linearize) then

      ! Copy linearization parameters
      p%Lin%NumTimes = p_FAST%NLinTimes
      p%Lin%InterpOrder = p_FAST%InterpOrder
      if (allocated(m_FAST%Lin%LinTimes)) then
         y%Lin%Times = m_FAST%Lin%LinTimes
      end if

      ! Initialize indices
      m%Lin%TimeIndex = 1
      m%Lin%AzimuthIndex = 1

      ! Set flag to save operating points during linearization if mode shapes requested
      p%Lin%SaveOPs = p_FAST%WrVTK == VTK_ModeShapes

      ! Initialize arrays to store operating point states and input
      call AllocAry(y%Lin%x, m%ModGlue%Vars%Nx, p%Lin%NumTimes, "Lin%x", ErrStat2, ErrMsg2); if (Failed()) return
      call AllocAry(y%Lin%z, m%ModGlue%Vars%Nz, p%Lin%NumTimes, "Lin%z", ErrStat2, ErrMsg2); if (Failed()) return
      call AllocAry(y%Lin%u, m%ModGlue%Vars%Nu, p%Lin%NumTimes, "Lin%u", ErrStat2, ErrMsg2); if (Failed()) return

   end if

   ! If linearization and steady state calculation is enabled
   if (p_FAST%Linearize .and. p_FAST%CalcSteady) then

      ! Disable saving of OPs during linearization as ModGlue_CalcSteady saves them automatically
      p%Lin%SaveOPs = .false.

      ! Initialize variables
      m%CS%AzimuthDelta = TwoPi_D/p%Lin%NumTimes
      m%CS%NumRotations = 0
      m%CS%IsConverged = .false.
      m%CS%FoundSteady = .false.
      m%CS%ForceLin = .false.

      ! Calculate number of output values (ignoring write outputs)
      m%CS%NumOutputs = 0
      do i = 1, size(m%ModGlue%Vars%y)
         associate (Var => m%ModGlue%Vars%y(i))
            if (.not. MV_HasFlagsAll(Var, VF_WriteOut)) m%CS%NumOutputs = m%CS%NumOutputs + Var%Num
         end associate
      end do

      ! Allocate arrays
      call AllocAry(y%Lin%Times, p%Lin%NumTimes, "Lin%Times", ErrStat2, ErrMsg2); if (Failed()) return
      call AllocAry(m%CS%AzimuthTarget, p%Lin%NumTimes, "CS%AzimuthTarget", ErrStat2, ErrMsg2); if (Failed()) return
      call AllocAry(m%CS%psi_buffer, p_FAST%LinInterpOrder + 1, "CS%psi_buffer", ErrStat2, ErrMsg2); if (Failed()) return
      call AllocAry(m%CS%y_buffer, m%ModGlue%Vars%Ny, p_FAST%LinInterpOrder + 1, "CS%y_buffer", ErrStat2, ErrMsg2); if (Failed()) return
      call AllocAry(m%CS%y_interp, m%ModGlue%Vars%Ny, "CS%y_interp", ErrStat2, ErrMsg2); if (Failed()) return
      call AllocAry(m%CS%y_diff, m%ModGlue%Vars%Ny, "CS%y_diff", ErrStat2, ErrMsg2); if (Failed()) return
      call AllocAry(m%CS%y_azimuth, m%ModGlue%Vars%Ny, p%Lin%NumTimes, "CS%y_azimuth", ErrStat2, ErrMsg2); if (Failed()) return
      call AllocAry(m%CS%y_ref, m%ModGlue%Vars%Ny, "CS%y_ref", ErrStat2, ErrMsg2); if (Failed()) return

      ! Initialize arrays to zero
      m%CS%psi_buffer = 0.0_R8Ki
      m%CS%y_buffer = 0.0_R8Ki
      m%CS%y_interp = 0.0_R8Ki
      m%CS%y_diff = 0.0_R8Ki
      m%CS%y_azimuth = 0.0_R8Ki
      m%CS%y_ref = 1.0_R8Ki

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

subroutine ModGlue_CalcSteady(n_t_global, t_global, p, m, y, p_FAST, m_FAST, T, ErrStat, ErrMsg)

   integer(IntKi), intent(IN)                :: n_t_global     !< integer time step
   real(DbKi), intent(IN)                    :: t_global       !< current simulation time
   type(Glue_ParameterType), intent(inout)   :: p              !< Glue Parameters
   type(Glue_MiscVarType), intent(inout)     :: m              !< Glue MiscVars
   type(Glue_OutputFileType), intent(inout)  :: y              !< Glue Output
   type(FAST_ParameterType), intent(inout)   :: p_FAST         !< FAST Parameters
   type(FAST_MiscVarType), intent(inout)     :: m_FAST         !< FAST MiscVars
   type(FAST_TurbineType), intent(inout)     :: T              !< Turbine Type
   integer(IntKi), intent(OUT)               :: ErrStat        !< Error status of the operation
   character(*), intent(OUT)                 :: ErrMsg         !< Error message if ErrStat /= ErrID_None

   character(*), parameter :: RoutineName = 'ModGlue_CalcSteady'
   integer(IntKi)          :: ErrStat2
   character(ErrMsgLen)    :: ErrMsg2
   real(DbKi)              :: DeltaAzimuth, AzimuthTargetDelta, AzimuthTarget
   real(DbKi)              :: psi                              !< psi (rotor azimuth) at which the outputs are defined
   real(DbKi)              :: error
   logical                 :: ProcessAzimuth
   integer(IntKi)          :: i, j, iy

   ErrStat = ErrID_None
   ErrMsg = ""

   ! Get current azimuth angle from ElastoDyn output
   psi = real(T%ED%y%LSSTipPxa, R8Ki)
   call Zero2TwoPi(psi)

   ! Cyclic shift psi buffer and set first index to new psi
   do i = size(m%CS%psi_buffer) - 1, 1, -1
      m%CS%psi_buffer(i + 1) = m%CS%psi_buffer(i)
   end do
   ! If passing the 2PI boundary, subtract 2PI from saved values so interpolation works correctly
   if (psi < m%CS%psi_buffer(1)) m%CS%psi_buffer = m%CS%psi_buffer - TwoPi_D
   m%CS%psi_buffer(1) = psi

   ! Cyclic shift output buffer and collect outputs from all modules
   do i = size(m%CS%psi_buffer) - 1, 1, -1
      m%CS%y_buffer(:, i + 1) = m%CS%y_buffer(:, i)
   end do

   ! Loop through modules and collect output

   do j = 1, size(m%ModGlue%ModData)
      associate (ModData => m%ModGlue%ModData(j))

         ! Skip of module has no outputs
         if (size(ModData%Vars%y) == 0) cycle

         ! Get outputs
         call FAST_GetOP(ModData, t_global, INPUT_CURR, STATE_CURR, T, ErrStat2, ErrMsg2, y_op=m%ModGlue%Lin%y, y_glue=m%ModGlue%Lin%y)
         if (Failed()) return

      end associate
   end do

   ! Copy outputs to buffer (can't be used directly since it's not allocatable)
   m%CS%y_buffer(:, 1) = m%ModGlue%Lin%y

   ! If first call
   if (n_t_global == 0) then

      ! Initialize azimuth targets
      do i = 1, p%Lin%NumTimes
         m%CS%AzimuthTarget(i) = (i - 1)*m%CS%AzimuthDelta + psi
         call Zero2TwoPi(m%CS%AzimuthTarget(i))
      end do

      ! Initialize psi buffer for interpolation based on time step and rotor speed
      do i = 1, size(m%CS%psi_buffer)
         m%CS%psi_buffer(i) = psi - (i - 1)*p_FAST%DT*T%ED%y%LSS_Spd
      end do

      ! Initialize output buffer by copying outputs from first buffer location
      do i = 2, size(m%CS%y_buffer, 2)
         m%CS%y_buffer(:, i) = m%CS%y_buffer(:, 1)
      end do

   end if

   ! Calculate change in azimuth from last call, if change is too great, return error
   DeltaAzimuth = psi - m%CS%psi_buffer(1)
   call Zero2TwoPi(DeltaAzimuth)
   if (DeltaAzimuth > m%CS%AzimuthDelta) then
      call SetErrStat(ErrID_Fatal, "The rotor is spinning too fast. The time step or NLinTimes is too large when CalcSteady=true.", ErrStat, ErrMsg, RoutineName)
      return
   end if

   ! Get the current azimuth target
   AzimuthTarget = m%CS%AzimuthTarget(m%Lin%AzimuthIndex)

   ! Difference between current azimuth and the target
   AzimuthTargetDelta = psi - AzimuthTarget

   ! Set flag to process next azimuth if psi is greater than the next azimuth target
   ! and the difference between psi and the target is less than the AzimuthDelta (difference between targets)
   ProcessAzimuth = (AzimuthTargetDelta >= 0.0_R8Ki) .and. (AzimuthTargetDelta < m%CS%AzimuthDelta)

   ! If this is the last step, force linearization
   if (t_global >= p_FAST%TMax - 0.5_DbKi*p_FAST%DT) then
      m%CS%ForceLin = .true.
      m%Lin%AzimuthIndex = 1
      ProcessAzimuth = .true.
   end if

   ! If flag is set to process azimuth
   if (ProcessAzimuth) then

      ! Interpolate outputs to target azimuth
      call MV_ExtrapInterp(m%ModGlue%Vars%y, m%CS%y_buffer, m%CS%psi_buffer, &
                           m%CS%y_interp, AzimuthTarget, ErrStat2, ErrMsg2)

      ! If converged
      if (m%CS%IsConverged) then

         ! Calculate error between interpolated outputs and outputs at this
         ! azimuth from the previous rotation
         error = CalcOutputErrorAtAzimuth()

         ! Update converged flag based on error and tolerance
         m%CS%IsConverged = error < p_FAST%TrimTol
      end if

      ! Save interpolated outputs for this azimuth
      m%CS%y_azimuth(:, m%Lin%AzimuthIndex) = m%CS%y_interp

      ! If linearization is forced
      if (m%CS%ForceLin) m%CS%IsConverged = .true.

      ! If converged or in first rotation, save this operating point for linearization later
      if (m%CS%IsConverged .or. m%CS%NumRotations == 0) then !
         y%Lin%Times(m%Lin%AzimuthIndex) = t_global
         call ModGlue_SaveOperatingPoint(p, m, m%Lin%AzimuthIndex, m%CS%NumRotations == 0, T, ErrStat2, ErrMsg2)
         if (Failed()) return
      end if

      ! Increment the azimuth index counter
      m%Lin%AzimuthIndex = m%Lin%AzimuthIndex + 1

      ! If we've completed one rotor revolution
      if (m%Lin%AzimuthIndex > p%Lin%NumTimes) then

         ! Increment number of rotations
         m%CS%NumRotations = m%CS%NumRotations + 1

         ! Save if steady state has been found
         m%CS%FoundSteady = m%CS%IsConverged

         ! If steady state has been found, return
         if (m%CS%FoundSteady) return

         ! Compute the reference values for this rotor revolution
         m%CS%y_ref = max(maxval(m%CS%y_azimuth, dim=2) - minval(m%CS%y_azimuth, dim=2), 0.01_R8Ki)

         ! Check errors next rotor revolution
         m%CS%IsConverged = .true.

         ! Reset the azimuth index
         m%Lin%AzimuthIndex = 1

         ! Forcing linearization if time is close to tmax (with sufficient margin)

         ! If rotor has nonzero speed
         if (T%ED%p%RotSpeed > 0) then

            ! If simulation is at least 10 revolutions, and error in rotor speed less than 0.1%
            if ((p_FAST%TMax > 10*(TwoPi_D)/T%ED%p%RotSpeed) .and. &
                (t_global >= p_FAST%TMax - 2._DbKi*(TwoPi_D)/T%ED%p%RotSpeed)) then
               if (abs(T%ED%y%RotSpeed - T%ED%p%RotSpeed)/T%ED%p%RotSpeed < 0.001) then
                  m%CS%ForceLin = .true.
               end if
            end if
         else
            if (t_global >= p_FAST%TMax - 1.5_DbKi*p_FAST%DT) then
               m%CS%ForceLin = .true.
            end if
         end if

      end if
   end if

   ! If linearization is being forced, set flags and display message
   if (m%CS%ForceLin) then
      m%CS%IsConverged = .true.
      m%CS%FoundSteady = .true.
      call WrScr('')
      call WrScr('[WARNING] Steady state not found before end of simulation. Forcing linearization.')
   end if

contains

   function CalcOutputErrorAtAzimuth() result(eps_squared)
      real(R8Ki)  :: eps_squared_sum, eps_squared

      ! Calculate difference between interpolated outputs for this rotation and
      ! interpolated outputs from previous rotation
      call MV_ComputeDiff(m%ModGlue%Vars%y, m%CS%y_interp, m%CS%y_azimuth(:, m%Lin%AzimuthIndex), m%CS%y_diff)

      ! Initialize epsilon squared sum
      eps_squared_sum = 0

      ! Loop through glue output variables
      do i = 1, size(m%ModGlue%Vars%y)
         associate (Var => m%ModGlue%Vars%y(i))

            ! Skip write outputs
            if (MV_HasFlagsAll(Var, VF_WriteOut)) cycle

            ! Loop through values in variable
            do j = Var%iLoc(1), Var%iLoc(2)

               ! If difference is not essentially zero, sum difference
               if (.not. EqualRealNos(m%CS%y_diff(j), 0.0_R8Ki)) then
                  eps_squared_sum = eps_squared_sum + (m%CS%y_diff(j)/m%CS%y_ref(j))**2
               end if
            end do
         end associate
      end do

      ! Normalize error by number of outputs
      eps_squared = eps_squared_sum/m%CS%NumOutputs
   end function

   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function Failed
end subroutine

subroutine ModGlue_Linearize_OP(p, m, y, p_FAST, m_FAST, y_FAST, t_global, Turbine, ErrStat, ErrMsg)

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

   character(*), parameter                   :: RoutineName = 'ModGlue_Linearize_OP'
   integer(IntKi)                            :: ErrStat2
   character(ErrMsgLen)                      :: ErrMsg2
   integer(IntKi)                            :: i, j, k
   integer(IntKi)                            :: ix, iz, iu, iy
   integer(IntKi)                            :: Un
   integer(IntKi)                            :: StateLinIndex, InputLinIndex
   character(200)                            :: SimStr
   character(MaxWrScrLen)                    :: BlankLine
   character(1024)                           :: LinRootName
   character(*), parameter                   :: Fmt = 'F10.2'

   ! Initialize error return
   ErrStat = ErrID_None
   ErrMsg = ""

   ! Write message to screen
   BlankLine = ""
   call WrOver(BlankLine)  ! BlankLine contains MaxWrScrLen spaces
   SimStr = '(RotSpeed='//trim(Num2LStr(Turbine%ED%y%RotSpeed*RPS2RPM, Fmt))//' rpm, BldPitch1='//trim(Num2LStr(Turbine%ED%y%BlPitch(1)*R2D, Fmt))//' deg)'
   call WrOver(' Performing linearization '//trim(Num2LStr(m%Lin%TimeIndex))//' at simulation time '//TRIM(Num2LStr(t_global))//' s. '//trim(SimStr))
   call WrScr('')

   !----------------------------------------------------------------------------
   ! Save operating point
   !----------------------------------------------------------------------------

   ! If flag set to save operating points during linearization
   if (p%Lin%SaveOPs) then
      call ModGlue_SaveOperatingPoint(p, m, m%Lin%TimeIndex, .true., Turbine, ErrStat2, ErrMsg2)
      if (Failed()) return
   end if

   !----------------------------------------------------------------------------
   ! Initialization
   !----------------------------------------------------------------------------

   ! Get parameters
   y_FAST%Lin%RotSpeed = Turbine%ED%y%RotSpeed
   y_FAST%Lin%Azimuth = Turbine%ED%y%LSSTipPxa

   ! Assemble linearization root file name
   LinRootName = trim(p_FAST%OutFileRoot)//'.'//trim(Num2LStr(m%Lin%TimeIndex))

   ! Get unit number for writing files
   call GetNewUnit(Un, ErrStat2, ErrMsg2); if (Failed()) return

   ! Initialize the index numbers
   ix = 1
   iz = 1
   iu = 1
   iy = 1

   ! Initialize data in Jacobian matrices to zero
   m%ModGlue%Lin%dYdu = 0.0_R8Ki
   m%ModGlue%Lin%dXdu = 0.0_R8Ki
   m%ModGlue%Lin%dYdx = 0.0_R8Ki
   m%ModGlue%Lin%dXdx = 0.0_R8Ki

   ! Loop through linearization modules by index
   do i = 1, size(m%ModGlue%ModData)
      associate (ModData => m%ModGlue%ModData(i))

         ! Derivatives with respect to input
         call FAST_JacobianPInput(ModData, t_global, INPUT_CURR, STATE_CURR, Turbine, ErrStat2, ErrMsg2, &
                                  dYdu=ModData%Lin%dYdu, dYdu_glue=m%ModGlue%Lin%dYdu, &
                                  dXdu=ModData%Lin%dXdu, dXdu_glue=m%ModGlue%Lin%dXdu)
         if (Failed()) return

         ! Derivatives with respect to continuous state
         call FAST_JacobianPContState(ModData, t_global, INPUT_CURR, STATE_CURR, Turbine, ErrStat2, ErrMsg2, &
                                      dYdx=ModData%Lin%dYdx, dYdx_glue=m%ModGlue%Lin%dYdx, &
                                      dXdx=ModData%Lin%dXdx, dXdx_glue=m%ModGlue%Lin%dXdx)
         if (Failed()) return

         ! Operating point values (must come after Jacobian routines because
         ! some modules calculate OP in those routines [MD])
         call FAST_GetOP(ModData, t_global, INPUT_CURR, STATE_CURR, Turbine, ErrStat2, ErrMsg2, &
                         u_op=ModData%Lin%u, u_glue=m%ModGlue%Lin%u, &
                         y_op=ModData%Lin%y, y_glue=m%ModGlue%Lin%y, &
                         x_op=ModData%Lin%x, x_glue=m%ModGlue%Lin%x, &
                         dx_op=ModData%Lin%dx, dx_glue=m%ModGlue%Lin%dx)
         if (Failed()) return

         ! If requested, write the module linearization matrices was requested
         if (p_FAST%LinOutMod) then
            call CalcWriteLinearMatrices(ModData%Vars, ModData%Lin, p_FAST, y_FAST, t_global, Un, &
                                         LinRootName, VF_Linearize, ErrStat2, ErrMsg2, ModSuffix=ModData%Abbr)
            if (Failed()) return
         end if

      end associate
   end do

   ! Copy arrays into linearization operating points
   if (allocated(m%ModGlue%Lin%x)) y%Lin%x(:, m%Lin%TimeIndex) = m%ModGlue%Lin%x
   if (allocated(m%ModGlue%Lin%z)) y%Lin%z(:, m%Lin%TimeIndex) = m%ModGlue%Lin%z
   if (allocated(m%ModGlue%Lin%u)) y%Lin%u(:, m%Lin%TimeIndex) = m%ModGlue%Lin%u

   ! Linearize mesh mappings to populate dUdy and dUdu
   call FAST_LinearizeMappings(m%ModGlue, m%Mappings, Turbine, ErrStat2, ErrMsg2)
   if (Failed()) return

   ! Write glue code matrices to file
   call CalcWriteLinearMatrices(m%ModGlue%Vars, m%ModGlue%Lin, p_FAST, y_FAST, t_global, Un, LinRootName, VF_Linearize, ErrStat2, ErrMsg2)
   if (Failed()) return

   ! Update index for next linearization time
   m%Lin%TimeIndex = m%Lin%TimeIndex + 1

contains
   logical function JacobianHasNaNs(Jac, label, abbr)
      real(R8Ki), allocatable, intent(in) :: Jac(:, :)
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

subroutine ModGlue_SaveOperatingPoint(p, m, OPIndex, NewCopy, Turbine, ErrStat, ErrMsg)
   type(Glue_ParameterType), intent(in)   :: p
   type(Glue_MiscVarType), intent(inout)  :: m
   integer(IntKi), intent(in)             :: OPIndex
   logical, intent(in)                    :: NewCopy
   type(FAST_TurbineType), intent(inout)  :: Turbine
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg

   character(*), parameter  :: RoutineName = 'ModGlue_SaveOperatingPoint'
   integer(IntKi)           :: ErrStat2
   character(ErrMsgLen)     :: ErrMsg2
   integer(IntKi)           :: StateIndex, InputIndex, CtrlCode, i

   ErrStat = ErrID_None
   ErrMsg = ''

   ! Set CtrlCode based on NewCopy flag
   if (NewCopy) then
      CtrlCode = MESH_NEWCOPY
   else
      CtrlCode = MESH_UPDATECOPY
   end if

   ! Index into state array where linearization data will be stored for this OP
   StateIndex = NumStateTimes + OPIndex

   ! Index into input save array where linearization data will be stored for OP
   InputIndex = Turbine%p_FAST%InterpOrder + 1 + OPIndex

   ! Loop through modules by index
   do i = 1, size(p%Lin%iMod)
      associate (ModData => m%ModData(p%Lin%iMod(i)))

         ! Copy current module state to linearization save location
         call FAST_CopyStates(ModData, Turbine, STATE_CURR, StateIndex, CtrlCode, ErrStat2, ErrMsg2)
         if (Failed()) return

         ! Copy current module output to linearization save location
         call FAST_CopyInput(ModData, Turbine, INPUT_CURR, -InputIndex, CtrlCode, ErrStat2, ErrMsg2)
         if (Failed()) return

      end associate
   end do

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function Failed
end subroutine

subroutine ModGlue_RestoreOperatingPoint(p, m, OPIndex, Turbine, ErrStat, ErrMsg)
   type(Glue_ParameterType), intent(in)   :: p
   type(Glue_MiscVarType), intent(inout)  :: m
   integer(IntKi), intent(in)             :: OPIndex
   type(FAST_TurbineType), intent(inout)  :: Turbine
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg

   character(*), parameter  :: RoutineName = 'ModGlue_RestoreOperatingPoint'
   integer(IntKi)           :: ErrStat2
   character(ErrMsgLen)     :: ErrMsg2
   integer(IntKi)           :: StateIndex, InputIndex, i

   ErrStat = ErrID_None
   ErrMsg = ''

   ! Index into state array where linearization data will be stored for this OP
   StateIndex = NumStateTimes + OPIndex

   ! Index into input save array where linearization data will be stored for OP
   InputIndex = Turbine%p_FAST%InterpOrder + 1 + OPIndex

   ! Loop through modules by index
   do i = 1, size(p%Lin%iMod)
      associate (ModData => m%ModData(p%Lin%iMod(i)))

         ! Copy current module state to linearization save location
         call FAST_CopyStates(ModData, Turbine, StateIndex, STATE_CURR, MESH_UPDATECOPY, ErrStat2, ErrMsg2)
         if (Failed()) return

         ! Copy current module input to linearization save location
         call FAST_CopyInput(ModData, Turbine, -InputIndex, INPUT_CURR, MESH_UPDATECOPY, ErrStat2, ErrMsg2)
         if (Failed()) return

      end associate
   end do

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function Failed
end subroutine

!> CalcGlueStateMatrices forms the full-system state matrices for linearization: A, B, C, and D.
!! Note that it uses LAPACK_GEMM instead of MATMUL for matrix multiplications because of stack-space issues (these
!! matrices get large quickly).
subroutine CalcGlueStateMatrices(Vars, Lin, JacScaleFactor, ErrStat, ErrMsg)
   type(ModVarsType), intent(in)    :: Vars           !< Glue variable data
   type(ModLinType), intent(inout)  :: Lin            !< Glue linearization data
   real(R8Ki), intent(in)           :: JacScaleFactor !< Scale factor for conditioning the Jacobians
   integer(IntKi), intent(out)      :: ErrStat        !< Error status of the operation
   character(*), intent(out)        :: ErrMsg         !< Error message if ErrStat /= ErrID_None

   character(*), parameter          :: RoutineName = 'CalcGlueStateMatrices'
   integer(IntKi)                   :: ErrStat2
   character(ErrMsgLen)             :: ErrMsg2
   real(R8Ki), allocatable          :: G(:, :), tmp(:, :)
   integer(IntKi), allocatable      :: ipiv(:)

   ! A = dXdx
   ! B = dXdu
   ! C = dYdx
   ! D = dYdu

   ! call DumpMatrix(1000, "dUdu.bin", Lin%dUdu, ErrStat2, ErrMsg2); if (Failed()) return
   ! call DumpMatrix(1000, "dUdy.bin", Lin%dUdy, ErrStat2, ErrMsg2); if (Failed()) return
   ! call DumpMatrix(1000, "A.bin", Lin%dXdx, ErrStat2, ErrMsg2); if (Failed()) return
   ! call DumpMatrix(1000, "B.bin", Lin%dXdu, ErrStat2, ErrMsg2); if (Failed()) return
   ! call DumpMatrix(1000, "C.bin", Lin%dYdx, ErrStat2, ErrMsg2); if (Failed()) return
   ! call DumpMatrix(1000, "D.bin", Lin%dYdu, ErrStat2, ErrMsg2); if (Failed()) return

   ! *** get G matrix ****
   !----------------------
   call AllocAry(G, size(Lin%dUdu, 1), size(Lin%dUdu, 2), 'G', ErrStat2, ErrMsg2); if (Failed()) return
   call AllocAry(ipiv, Vars%Nu, 'ipiv', ErrStat2, ErrMsg2); if (Failed()) return

   ! G = dUdu + matmul(dUdy, y_FAST%Lin%Glue%D)
   G = Lin%dUdu
   call LAPACK_GEMM('N', 'N', 1.0_R8Ki, Lin%dUdy, Lin%dYdu, 1.0_R8Ki, G, ErrStat2, ErrMsg2); if (Failed()) return

   ! G can be ill-conditioned, so we are going to precondition with G_hat = S^(-1) * G * S
   ! we will also multiply the right-hand-side of the equations that need G inverse so that
   ! dUdy_hat = S^(-1)*dUdy and dUdu_hat = S^(-1)*dUdu
   call Precondition(Vars%u, G, Lin%dUdu, Lin%dUdy, JacScaleFactor)

   ! Form G_hat^(-1) * (S^-1*dUdy) and G^(-1) * (S^-1*dUdu)
   ! factor G for the two solves:
   call LAPACK_getrf(M=size(G, 1), N=size(G, 2), A=G, IPIV=ipiv, ErrStat=ErrStat2, ErrMsg=ErrMsg2); if (Failed()) return

   ! after the this solve, dUdy holds G_hat^(-1) * dUdy_hat:
   call LAPACK_getrs(trans='N', N=size(G, 2), A=G, IPIV=ipiv, B=Lin%dUdy, ErrStat=ErrStat2, ErrMsg=ErrMsg2); if (Failed()) return

   ! after the this solve, dUdu holds G_hat^(-1) * dUdu_hat:
   call LAPACK_getrs(trans='N', N=size(G, 2), A=G, IPIV=ipiv, B=Lin%dUdu, ErrStat=ErrStat2, ErrMsg=ErrMsg2); if (Failed()) return

   ! Deallocate G and ipiv because the solves are complete
   deallocate (G)
   deallocate (ipiv)

   ! After this call, dUdu holds G^(-1)*dUdu and dUdy holds G^(-1)*dUdy
   call Postcondition(Vars%u, Lin%dUdu, Lin%dUdy, JacScaleFactor)

   ! Allocate tmp matrix for A and C calculations
   call AllocAry(tmp, Vars%Nu, Vars%Nx, 'G^-1*dUdy*C', ErrStat2, ErrMsg2); if (Failed()) return

   ! tmp = G^(-1) * dUdy * diag(C)
   call LAPACK_GEMM('N', 'N', 1.0_R8Ki, Lin%dUdy, Lin%dYdx, 0.0_R8Ki, tmp, ErrStat2, ErrMsg2); if (Failed()) return

   ! A
   ! dXdx = dXdx - matmul(dXdu, tmp)
   call LAPACK_GEMM('N', 'N', -1.0_R8Ki, Lin%dXdu, tmp, 1.0_R8Ki, Lin%dXdx, ErrStat2, ErrMsg2); if (Failed()) return

   ! C
   ! dYdx = dYdx - matmul(dYdu, tmp)
   call LAPACK_GEMM('N', 'N', -1.0_R8Ki, Lin%dYdu, tmp, 1.0_R8Ki, Lin%dYdx, ErrStat2, ErrMsg2); if (Failed()) return

   ! B
   tmp = Lin%dXdu
   ! dXdu = matmul(dXdu, dUdu)
   call LAPACK_GEMM('N', 'N', 1.0_R8Ki, tmp, Lin%dUdu, 0.0_R8Ki, Lin%dXdu, ErrStat2, ErrMsg2); if (Failed()) return

   ! D
   tmp = Lin%dYdu
   ! D = matmul(dYdu, dUdu)
   call LAPACK_GEMM('N', 'N', 1.0_R8Ki, tmp, Lin%dUdu, 0.0_R8Ki, Lin%dYdu, ErrStat2, ErrMsg2); if (Failed()) return

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

   allocate (isLoad(size(dUdu, 1)))
   isLoad = .false.

   ! Loop through glue code input variables (cols)
   do i = 1, size(uVars)

      ! Get if col variable is a load
      isColLoad = MV_IsLoad(uVars(i))

      ! Get col variable start and end indices in matrix
      associate (iLoc => uVars(i)%iLoc)

         isLoad(iLoc(1):iLoc(2)) = isColLoad

         ! Loop through glue code input variables (rows)
         do j = 1, size(uVars)

            ! Get if row variable is a load
            isRowLoad = MV_IsLoad(uVars(j))

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
      if (uVars(i)%Field == FieldForce .or. uVars(i)%Field == FieldMoment) then

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

subroutine CalcWriteLinearMatrices(Vars, Lin, p_FAST, y_FAST, t_global, Un, LinRootName, FilterFlag, ErrStat, ErrMsg, ModSuffix, CalcGlue, FullOutput)
   type(ModVarsType), intent(in)             :: Vars           !< Variable data
   type(ModLinType), intent(inout)           :: Lin            !< Linearization data
   type(FAST_ParameterType), intent(in)      :: p_FAST         !< Parameters
   type(FAST_OutputFileType), intent(in)     :: y_FAST         !< Output variables
   real(DbKi), intent(in)                    :: t_global       !< current time step (written in file)
   integer(IntKi), intent(in)                :: Un             !< Unit number for file
   character(*), intent(in)                  :: LinRootName    !< output file name
   integer(IntKi), intent(in)                :: FilterFlag     !< Variable flag for filtering
   integer(IntKi), intent(out)               :: ErrStat        !< Error status of the operation
   character(*), intent(out)                 :: ErrMsg         !< Error message if ErrStat /= ErrID_None
   character(*), optional, intent(in)        :: ModSuffix      !< Module suffix for file name
   logical, optional, intent(in)             :: CalcGlue       !< Flag to calculate glue state matrices
   logical, optional, intent(in)             :: FullOutput     !< Flag to output all Jacobians

   character(*), parameter          :: RoutineName = 'WriteModuleLinearMatrices'
   integer(IntKi)                   :: ErrStat2
   character(ErrMsgLen)             :: ErrMsg2
   character(32)                    :: Desc
   character(1024)                  :: OutFileName
   integer(IntKi)                   :: i
   integer(IntKi)                   :: Nx, Nxd, Nz, Nu, Ny
   character(50)                    :: Fmt
   logical, allocatable             :: uUse(:), yUse(:), xUse(:)
   logical                          :: CalcGlueLoc, FullOutputLoc

   ErrStat = ErrID_None
   ErrMsg = ""

   ! Assemble output file name based on glue linearization abbreviation
   if (present(ModSuffix)) then
      OutFileName = trim(LinRootName)//"."//trim(ModSuffix)//".lin"
      CalcGlueLoc = .false.
   else
      OutFileName = trim(LinRootName)//".lin"
      CalcGlueLoc = .true.
   end if

   if (present(FullOutput)) then
      FullOutputLoc = FullOutput
   else
      FullOutputLoc = p_FAST%LinOutJac
   end if

   ! Set flag to calculate glue matrices based on optional parameter
   if (present(CalcGlue)) CalcGlueLoc = CalcGlue

   ! Open linearization file
   call OpenFOutFile(Un, OutFileName, ErrStat2, ErrMsg2); if (Failed()) return

   ! Calculate number of values in variable after applying filter
   Nx = MV_NumVals(Vars%x, FilterFlag)
   Nxd = 0
   Nz = MV_NumVals(Vars%z, FilterFlag)
   Nu = MV_NumVals(Vars%u, FilterFlag)
   Ny = MV_NumVals(Vars%y, FilterFlag)

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

   if (Nx > 0 .and. allocated(Lin%x)) then
      write (Un, '(A)') 'Order of continuous states:'
      call WrLinFile_txt_Table(Vars%x, FilterFlag, p_FAST, Un, "Row/Column", Lin%x)
   end if

   if (Nx > 0 .and. allocated(Lin%dx)) then
      write (Un, '(A)') 'Order of continuous state derivatives:'
      call WrLinFile_txt_Table(Vars%x, FilterFlag, p_FAST, Un, "Row/Column", Lin%dx, IsDeriv=.true.)
   end if

   if (Nz > 0 .and. allocated(Lin%z)) then
      write (Un, '(A)') 'Order of constraint states:'
      call WrLinFile_txt_Table(Vars%z, FilterFlag, p_FAST, Un, "Row/Column", Lin%z)
   end if

   if (Nu > 0 .and. allocated(Lin%u)) then
      write (Un, '(A)') 'Order of inputs:'
      call WrLinFile_txt_Table(Vars%u, FilterFlag, p_FAST, Un, "Column  ", Lin%u, ShowRot=.true.)
   end if

   if (Ny > 0 .and. allocated(Lin%y)) then
      write (Un, '(A)') 'Order of outputs:'
      call WrLinFile_txt_Table(Vars%y, FilterFlag, p_FAST, Un, "Row  ", Lin%y, ShowRot=.true.)
   end if

   ! Create boolean array indicating which continuous state values to write
   allocate (xUse(Vars%Nx))
   xUse = .false.
   do i = 1, size(Vars%x)
      associate (Var => Vars%x(i))
         if (MV_HasFlagsAll(Var, FilterFlag)) xUse(Var%iLoc(1):Var%iLoc(2)) = .true.
      end associate
   end do

   ! Create boolean array indicating which input values to write
   allocate (uUse(Vars%Nu))
   uUse = .false.
   do i = 1, size(Vars%u)
      associate (Var => Vars%u(i))
         if (MV_HasFlagsAll(Var, FilterFlag)) uUse(Var%iLoc(1):Var%iLoc(2)) = .true.
      end associate
   end do

   ! Create boolean array indicating which output values to write
   allocate (yUse(Vars%Ny))
   yUse = .false.
   do i = 1, size(Vars%y)
      associate (Var => Vars%y(i))
         if (MV_HasFlagsAll(Var, FilterFlag)) yUse(Var%iLoc(1):Var%iLoc(2)) = .true.
      end associate
   end do

   ! If Jacobian matrix output is requested
   if (FullOutputLoc) then
      write (Un, '(/,A,/)') 'Jacobian matrices:'
      if (allocated(Lin%dUdu)) call WrPartialMatrix(Lin%dUdu, Un, p_FAST%OutFmt, 'dUdu', UseRow=uUse, UseCol=uUse)
      if (allocated(Lin%dUdy)) call WrPartialMatrix(Lin%dUdy, Un, p_FAST%OutFmt, 'dUdy', UseRow=uUse, UseCol=yUse)
      if (allocated(Lin%dXdy)) call WrPartialMatrix(Lin%dXdy, Un, p_FAST%OutFmt, 'dXdy', UseRow=xUse, UseCol=yUse)
      if (allocated(Lin%J)) call WrPartialMatrix(Lin%J, Un, p_FAST%OutFmt, 'J')
      if (present(ModSuffix)) then
         if (allocated(Lin%dXdx)) call WrPartialMatrix(Lin%dXdx, Un, p_FAST%OutFmt, 'dXdx', UseRow=xUse, UseCol=xUse)
         if (allocated(Lin%dXdu)) call WrPartialMatrix(Lin%dXdu, Un, p_FAST%OutFmt, 'dXdu', UseRow=xUse, UseCol=uUse)
         if (allocated(Lin%dYdx)) call WrPartialMatrix(Lin%dYdx, Un, p_FAST%OutFmt, 'dYdx', UseRow=yUse, UseCol=xUse)
         if (allocated(Lin%dYdu)) call WrPartialMatrix(Lin%dYdu, Un, p_FAST%OutFmt, 'dYdu', UseRow=yUse, UseCol=uUse)
      end if
   end if

   ! If this is glue code module, calculate the glue code state matrices (A, B, C, D)
   ! Called here, after writing dUdu and dUdy, because those matrices are overwritten
   ! in the process of calculating the other state matrices
   if (CalcGlueLoc) then
      call CalcGlueStateMatrices(Vars, Lin, real(p_FAST%UJacSclFact, R8Ki), ErrStat2, ErrMsg2)
      if (Failed()) return
   end if

   ! Write the linearized state matrices
   write (Un, '(/,A,/)') 'Linearized state matrices:'
   if (allocated(Lin%dXdx)) call WrPartialMatrix(Lin%dXdx, Un, p_FAST%OutFmt, 'A', UseRow=xUse, UseCol=xUse)
   if (allocated(Lin%dXdu)) call WrPartialMatrix(Lin%dXdu, Un, p_FAST%OutFmt, 'B', UseRow=xUse, UseCol=uUse)
   if (allocated(Lin%dYdx)) call WrPartialMatrix(Lin%dYdx, Un, p_FAST%OutFmt, 'C', UseRow=yUse, UseCol=xUse)
   if (allocated(Lin%dYdu)) call WrPartialMatrix(Lin%dYdu, Un, p_FAST%OutFmt, 'D', UseRow=yUse, UseCol=uUse)
   if (allocated(Lin%StateRotation)) call WrPartialMatrix(Lin%StateRotation, Un, p_FAST%OutFmt, 'StateRotation')

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
   real(R8Ki)                    :: DCM(3, 3), wm(3)
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
         if (.not. MV_HasFlagsAll(Var, FlagFilter)) cycle

         ! Is variable in the rotating frame?
         VarRotFrame = MV_HasFlagsAll(Var, VF_RotFrame)

         ! Get variable derivative order
         if (MV_HasFlagsAll(Var, VF_DerivOrder2)) then
            VarDerivOrder = 2
         else if (MV_HasFlagsAll(Var, VF_DerivOrder1)) then
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
            if (ShowRotLoc .and. (Var%Field == FieldOrientation)) then

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

            else if (MV_HasFlagsAll(Var, VF_WM_Rot)) then ! BeamDyn Wiener-Milenkovic orientation

               ! Skip writing if not the first value in orientation (3 values)
               if (mod(j - 1, 3) /= 0) cycle

               ! Convert from quaternion in operating point to BeamDyn WM parameter
               wm = -quat_to_wm(op(i_op:i_op + 2))

               ! Write all components of WM parameters
               write (Un, Fmt) RowColIdx, wm(1), VarRotFrame, VarDerivOrder, trim(Var%LinNames(j))
               write (Un, Fmt) RowColIdx, wm(2), VarRotFrame, VarDerivOrder, trim(Var%LinNames(j))
               write (Un, Fmt) RowColIdx, wm(3), VarRotFrame, VarDerivOrder, trim(Var%LinNames(j))

            else

               write (Un, Fmt) RowColIdx, op(i_op), VarRotFrame, VarDerivOrder, trim(Var%LinNames(j))

            end if

         end do
      end associate
   end do

   write (Un, '()')    !print a blank line

end subroutine WrLinFile_txt_Table

end module
