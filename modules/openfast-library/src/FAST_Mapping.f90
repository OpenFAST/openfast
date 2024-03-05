!**********************************************************************************************************************************
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
!> This module contains routines used by FAST to map meshes and values between modules for transfering data and doing linearization.

module FAST_Mapping

use FAST_Types
use FAST_ModTypes

implicit none

private
public :: FAST_InitMappings, FAST_LinearizeMappings, FAST_ResetRemapFlags, FAST_InputSolve

integer(IntKi), parameter  :: AD_rotor = 1

integer(IntKi), parameter  :: Xfr_Invalid = 0, &
                              Xfr_Point_to_Point = 1, &
                              Xfr_Line2_to_Point = 2, &
                              Xfr_Point_to_Line2 = 3, &
                              Xfr_Line2_to_Line2 = 4

character(24), parameter   :: Custom_ED_to_ExtLd = 'ED -> ExtLd', &
                              Custom_AD_to_ExtLd = 'AD -> ExtLd', &
                              Custom_SrvD_to_AD = 'SrvD -> AD', &
                              Custom_ED_to_IfW = 'ED -> IfW', &
                              Custom_SrvD_to_IfW = 'SrvD -> IfW', &
                              Custom_BD_to_SrvD = 'BD -> SrvD', &
                              Custom_ED_to_SrvD = 'ED -> SrvD', &
                              Custom_IfW_to_SrvD = 'IfW -> SrvD', &
                              Custom_ExtInfw_to_SrvD = 'ExtInfw -> SrvD', &
                              Custom_SrvD_to_SD = 'SrvD -> SD', &
                              Custom_SrvD_to_MD = 'SrvD -> MD', &
                              Custom_ExtLd_to_SrvD = 'ExtLd -> SrvD'

contains

subroutine FAST_InputMeshPointer(ModData, Turbine, MeshLoc, UseU, Mesh, ErrStat, ErrMsg)
   type(ModDataType), intent(in)                :: ModData
   type(FAST_TurbineType), target, intent(in)   :: Turbine
   type(MeshLocType), intent(in)                :: MeshLoc
   logical, intent(in)                          :: UseU
   type(MeshType), pointer, intent(out)         :: Mesh
   integer(IntKi), intent(out)                  :: ErrStat
   character(*), intent(out)                    :: ErrMsg

   ErrStat = ErrID_None
   ErrMsg = ""

   nullify (Mesh)

   select case (ModData%ID)
   case (Module_AD)
      if (UseU) then
         Mesh => AD_InputMeshPointer(Turbine%AD%u, MeshLoc)
      else
         Mesh => AD_InputMeshPointer(Turbine%AD%Input(1), MeshLoc)
      end if
   case (Module_BD)
      if (UseU) then
         Mesh => BD_InputMeshPointer(Turbine%BD%u(ModData%Ins), MeshLoc)
      else
         Mesh => BD_InputMeshPointer(Turbine%BD%Input(1, ModData%Ins), MeshLoc)
      end if
   case (Module_ED)
      if (UseU) then
         Mesh => ED_InputMeshPointer(Turbine%ED%u, MeshLoc)
      else
         Mesh => ED_InputMeshPointer(Turbine%ED%Input(1), MeshLoc)
      end if
   case (Module_ExtInfw)
      if (UseU) then
         Mesh => ExtInfw_InputMeshPointer(Turbine%ExtInfw%u, MeshLoc)
      else
         ! ExtInfw doesn't have the typical input structure, using u for both
         ! Mesh => ExtInfw_InputMeshPointer(Turbine%ExtInfw%Input(1), MeshLoc)
         Mesh => ExtInfw_InputMeshPointer(Turbine%ExtInfw%u, MeshLoc)
      end if
   case (Module_ExtPtfm)
      if (UseU) then
         Mesh => ExtPtfm_InputMeshPointer(Turbine%ExtPtfm%u, MeshLoc)
      else
         Mesh => ExtPtfm_InputMeshPointer(Turbine%ExtPtfm%Input(1), MeshLoc)
      end if
   case (Module_FEAM)
      if (UseU) then
         Mesh => FEAM_InputMeshPointer(Turbine%FEAM%u, MeshLoc)
      else
         Mesh => FEAM_InputMeshPointer(Turbine%FEAM%Input(1), MeshLoc)
      end if
   case (Module_HD)
      if (UseU) then
         Mesh => HydroDyn_InputMeshPointer(Turbine%HD%u, MeshLoc)
      else
         Mesh => HydroDyn_InputMeshPointer(Turbine%HD%Input(1), MeshLoc)
      end if
   case (Module_IceD)
      if (UseU) then
         Mesh => IceD_InputMeshPointer(Turbine%IceD%u(ModData%Ins), MeshLoc)
      else
         Mesh => IceD_InputMeshPointer(Turbine%IceD%Input(1, ModData%Ins), MeshLoc)
      end if
   case (Module_IceF)
      if (UseU) then
         Mesh => IceFloe_InputMeshPointer(Turbine%IceF%u, MeshLoc)
      else
         Mesh => IceFloe_InputMeshPointer(Turbine%IceF%Input(1), MeshLoc)
      end if
   case (Module_IfW)
      if (UseU) then
         Mesh => InflowWind_InputMeshPointer(Turbine%IfW%u, MeshLoc)
      else
         Mesh => InflowWind_InputMeshPointer(Turbine%IfW%Input(1), MeshLoc)
      end if
   case (Module_MAP)
      if (UseU) then
         Mesh => MAP_InputMeshPointer(Turbine%MAP%u, MeshLoc)
      else
         Mesh => MAP_InputMeshPointer(Turbine%MAP%Input(1), MeshLoc)
      end if
   case (Module_MD)
      if (UseU) then
         Mesh => MD_InputMeshPointer(Turbine%MD%u, MeshLoc)
      else
         Mesh => MD_InputMeshPointer(Turbine%MD%Input(1), MeshLoc)
      end if
   case (Module_Orca)
      if (UseU) then
         Mesh => Orca_InputMeshPointer(Turbine%Orca%u, MeshLoc)
      else
         Mesh => Orca_InputMeshPointer(Turbine%Orca%Input(1), MeshLoc)
      end if
   case (Module_SD)
      if (UseU) then
         Mesh => SD_InputMeshPointer(Turbine%SD%u, MeshLoc)
      else
         Mesh => SD_InputMeshPointer(Turbine%SD%Input(1), MeshLoc)
      end if
   case (Module_SeaSt)
      if (UseU) then
         Mesh => SeaSt_InputMeshPointer(Turbine%SeaSt%u, MeshLoc)
      else
         Mesh => SeaSt_InputMeshPointer(Turbine%SeaSt%Input(1), MeshLoc)
      end if
   case (Module_SrvD)
      if (UseU) then
         Mesh => SrvD_InputMeshPointer(Turbine%SrvD%u, MeshLoc)
      else
         Mesh => SrvD_InputMeshPointer(Turbine%SrvD%Input(1), MeshLoc)
      end if
   case default
      ErrStat = ErrID_Fatal
      ErrMsg = "Unsupported module: "//ModData%Abbr
      return
   end select

   if (.not. associated(Mesh)) then
      ErrStat = ErrID_Fatal
      ErrMsg = "Mesh not found in module "//ModData%Abbr// &
               ", Num="//trim(Num2LStr(MeshLoc%Num))// &
               ", i1="//trim(Num2LStr(MeshLoc%i1))// &
               ", i2="//trim(Num2LStr(MeshLoc%i2))// &
               ", i3="//trim(Num2LStr(MeshLoc%i3))
      return
   end if
end subroutine

subroutine FAST_OutputMeshPointer(ModData, Turbine, MeshLoc, Mesh, ErrStat, ErrMsg)
   type(ModDataType), intent(in)                :: ModData
   type(FAST_TurbineType), target, intent(in)   :: Turbine
   type(MeshLocType), intent(in)                :: MeshLoc
   type(MeshType), pointer, intent(out)         :: Mesh
   integer(IntKi), intent(out)                  :: ErrStat
   character(*), intent(out)                    :: ErrMsg

   ErrStat = ErrID_None
   ErrMsg = ""

   nullify (Mesh)

   select case (ModData%ID)
   case (Module_AD)
      Mesh => AD_OutputMeshPointer(Turbine%AD%y, MeshLoc)
   case (Module_BD)
      Mesh => BD_OutputMeshPointer(Turbine%BD%y(ModData%Ins), MeshLoc)
   case (Module_ED)
      Mesh => ED_OutputMeshPointer(Turbine%ED%y, MeshLoc)
   case (Module_ExtInfw)
      Mesh => ExtInfw_OutputMeshPointer(Turbine%ExtInfw%y, MeshLoc)
   case (Module_ExtPtfm)
      Mesh => ExtPtfm_OutputMeshPointer(Turbine%ExtPtfm%y, MeshLoc)
   case (Module_FEAM)
      Mesh => FEAM_OutputMeshPointer(Turbine%FEAM%y, MeshLoc)
   case (Module_HD)
      Mesh => HydroDyn_OutputMeshPointer(Turbine%HD%y, MeshLoc)
   case (Module_IceD)
      Mesh => IceD_OutputMeshPointer(Turbine%IceD%y(ModData%Ins), MeshLoc)
   case (Module_IceF)
      Mesh => IceFloe_OutputMeshPointer(Turbine%IceF%y, MeshLoc)
   case (Module_IfW)
      Mesh => InflowWind_OutputMeshPointer(Turbine%IfW%y, MeshLoc)
   case (Module_MAP)
      Mesh => MAP_OutputMeshPointer(Turbine%MAP%y, MeshLoc)
   case (Module_MD)
      Mesh => MD_OutputMeshPointer(Turbine%MD%y, MeshLoc)
   case (Module_Orca)
      Mesh => Orca_OutputMeshPointer(Turbine%Orca%y, MeshLoc)
   case (Module_SD)
      Mesh => SD_OutputMeshPointer(Turbine%SD%y, MeshLoc)
   case (Module_SeaSt)
      Mesh => SeaSt_OutputMeshPointer(Turbine%SeaSt%y, MeshLoc)
   case (Module_SrvD)
      Mesh => SrvD_OutputMeshPointer(Turbine%SrvD%y, MeshLoc)
   case default
      ErrStat = ErrID_Fatal
      ErrMsg = "Unsupported module: "//ModData%Abbr
      return
   end select

   if (.not. associated(Mesh)) then
      ErrStat = ErrID_Fatal
      ErrMsg = "Mesh not found in module "//ModData%Abbr// &
               ", Num="//trim(Num2LStr(MeshLoc%Num))// &
               ", i1="//trim(Num2LStr(MeshLoc%i1))// &
               ", i2="//trim(Num2LStr(MeshLoc%i2))// &
               ", i3="//trim(Num2LStr(MeshLoc%i3))
      return
   end if
end subroutine

function FAST_InputMeshName(ModData, MeshLoc) result(Name)
   type(ModDataType), intent(in)                :: ModData
   type(MeshLocType), intent(in)                :: MeshLoc
   character(32)                                :: Name
   Name = "Unknown mesh in "//ModData%Abbr
   select case (ModData%ID)
   case (Module_AD)
      Name = trim(ModData%Abbr)//"%"//AD_InputMeshName(MeshLoc)
   case (Module_BD)
      Name = trim(ModData%Abbr)//"("//trim(Num2LStr(ModData%Ins))//")%"//BD_InputMeshName(MeshLoc)
   case (Module_ED)
      Name = trim(ModData%Abbr)//"%"//ED_InputMeshName(MeshLoc)
   case (Module_ExtInfw)
      Name = trim(ModData%Abbr)//"%"//ExtInfw_InputMeshName(MeshLoc)
   case (Module_ExtPtfm)
      Name = trim(ModData%Abbr)//"%"//ExtPtfm_InputMeshName(MeshLoc)
   case (Module_FEAM)
      Name = trim(ModData%Abbr)//"%"//FEAM_InputMeshName(MeshLoc)
   case (Module_HD)
      Name = trim(ModData%Abbr)//"%"//HydroDyn_InputMeshName(MeshLoc)
   case (Module_IceD)
      Name = trim(ModData%Abbr)//"("//trim(Num2LStr(ModData%Ins))//")%"//IceD_InputMeshName(MeshLoc)
   case (Module_IceF)
      Name = trim(ModData%Abbr)//"%"//IceFloe_InputMeshName(MeshLoc)
   case (Module_IfW)
      Name = trim(ModData%Abbr)//"%"//InflowWind_InputMeshName(MeshLoc)
   case (Module_MAP)
      Name = trim(ModData%Abbr)//"%"//MAP_InputMeshName(MeshLoc)
   case (Module_MD)
      Name = trim(ModData%Abbr)//"%"//MD_InputMeshName(MeshLoc)
   case (Module_Orca)
      Name = trim(ModData%Abbr)//"%"//Orca_InputMeshName(MeshLoc)
   case (Module_SD)
      Name = trim(ModData%Abbr)//"%"//SD_InputMeshName(MeshLoc)
   case (Module_SeaSt)
      Name = trim(ModData%Abbr)//"%"//SeaSt_InputMeshName(MeshLoc)
   case (Module_SrvD)
      Name = trim(ModData%Abbr)//"%"//SrvD_InputMeshName(MeshLoc)
   end select
end function

function FAST_OutputMeshName(ModData, MeshLoc) result(Name)
   type(ModDataType), intent(in)                :: ModData
   type(MeshLocType), intent(in)                :: MeshLoc
   character(32)                                :: Name
   Name = "Unknown mesh in "//ModData%Abbr
   select case (ModData%ID)
   case (Module_AD)
      Name = trim(ModData%Abbr)//"%"//AD_OutputMeshName(MeshLoc)
   case (Module_BD)
      Name = trim(ModData%Abbr)//"("//trim(Num2LStr(ModData%Ins))//")%"//BD_OutputMeshName(MeshLoc)
   case (Module_ED)
      Name = trim(ModData%Abbr)//"%"//ED_OutputMeshName(MeshLoc)
   case (Module_ExtInfw)
      Name = trim(ModData%Abbr)//"%"//ExtInfw_OutputMeshName(MeshLoc)
   case (Module_ExtPtfm)
      Name = trim(ModData%Abbr)//"%"//ExtPtfm_OutputMeshName(MeshLoc)
   case (Module_FEAM)
      Name = trim(ModData%Abbr)//"%"//FEAM_OutputMeshName(MeshLoc)
   case (Module_HD)
      Name = trim(ModData%Abbr)//"%"//HydroDyn_OutputMeshName(MeshLoc)
   case (Module_IceD)
      Name = trim(ModData%Abbr)//"("//trim(Num2LStr(ModData%Ins))//")%"//IceD_OutputMeshName(MeshLoc)
   case (Module_IceF)
      Name = trim(ModData%Abbr)//"%"//IceFloe_OutputMeshName(MeshLoc)
   case (Module_IfW)
      Name = trim(ModData%Abbr)//"%"//InflowWind_OutputMeshName(MeshLoc)
   case (Module_MAP)
      Name = trim(ModData%Abbr)//"%"//MAP_OutputMeshName(MeshLoc)
   case (Module_MD)
      Name = trim(ModData%Abbr)//"%"//MD_OutputMeshName(MeshLoc)
   case (Module_Orca)
      Name = trim(ModData%Abbr)//"%"//Orca_OutputMeshName(MeshLoc)
   case (Module_SD)
      Name = trim(ModData%Abbr)//"%"//SD_OutputMeshName(MeshLoc)
   case (Module_SeaSt)
      Name = trim(ModData%Abbr)//"%"//SeaSt_OutputMeshName(MeshLoc)
   case (Module_SrvD)
      Name = trim(ModData%Abbr)//"%"//SrvD_OutputMeshName(MeshLoc)
   end select
end function

subroutine FAST_InitMappings(Mods, Mappings, Turbine, ErrStat, ErrMsg)
   type(ModDataType), intent(inout)                   :: Mods(:)     !< Module data
   type(TC_MappingType), allocatable, intent(inout)   :: Mappings(:)
   type(FAST_TurbineType), intent(inout)              :: Turbine     !< Turbine type
   integer(IntKi), intent(out)                        :: ErrStat
   character(*), intent(out)                          :: ErrMsg

   character(*), parameter    :: RoutineName = 'FAST_InitMappings'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2
   integer(IntKi)             :: i, j, k
   integer(IntKi)             :: iMap, ModIns, iModIn, iModSrc, iModDst

   ErrStat = ErrID_None
   ErrMsg = ''

   !----------------------------------------------------------------------------
   ! Define mesh mappings between modules
   !----------------------------------------------------------------------------

   ! Define a list of all possible module mesh mappings between modules
   allocate (Mappings(0), stat=ErrStat2)
   if (ErrStat2 /= 0) then
      call SetErrStat(ErrID_Fatal, "Error allocating mappings", ErrStat, ErrMsg, RoutineName)
      return
   end if

   ! Loop through module pairings
   do iModSrc = 1, size(Mods)
      do iModDst = 1, size(Mods)

         ! Switch by destination module (inputs)
         select case (Mods(IModDst)%ID)
         case (Module_AD)
            call InitMappings_AD(Mappings, Mods(iModSrc), Mods(iModDst), Turbine, ErrStat2, ErrMsg2)
         case (Module_BD)
            call InitMappings_BD(Mappings, Mods(iModSrc), Mods(iModDst), Turbine, ErrStat2, ErrMsg2)
         case (Module_ED)
            call InitMappings_ED(Mappings, Mods(iModSrc), Mods(iModDst), Turbine, ErrStat2, ErrMsg2)
         case (Module_ExtInfw)
            ! call InitMappings_ExtInfw(Mappings, Mods(iModSrc), Mods(iModDst), Turbine, ErrStat2, ErrMsg2)
         case (Module_ExtLd)
            call InitMappings_ExtLd(Mappings, Mods(iModSrc), Mods(iModDst), Turbine, ErrStat2, ErrMsg2)
         case (Module_ExtPtfm)
            call InitMappings_ExtPtfm(Mappings, Mods(iModSrc), Mods(iModDst), Turbine, ErrStat2, ErrMsg2)
         case (Module_FEAM)
            call InitMappings_FEAM(Mappings, Mods(iModSrc), Mods(iModDst), Turbine, ErrStat2, ErrMsg2)
         case (Module_HD)
            call InitMappings_HD(Mappings, Mods(iModSrc), Mods(iModDst), Turbine, ErrStat2, ErrMsg2)
         case (Module_IceD)
            call InitMappings_IceD(Mappings, Mods(iModSrc), Mods(iModDst), Turbine, ErrStat2, ErrMsg2)
         case (Module_IceF)
            call InitMappings_IceF(Mappings, Mods(iModSrc), Mods(iModDst), Turbine, ErrStat2, ErrMsg2)
         case (Module_IfW)
            call InitMappings_IfW(Mappings, Mods(iModSrc), Mods(iModDst), Turbine, ErrStat2, ErrMsg2)
         case (Module_MAP)
            call InitMappings_MAP(Mappings, Mods(iModSrc), Mods(iModDst), Turbine, ErrStat2, ErrMsg2)
         case (Module_MD)
            call InitMappings_MD(Mappings, Mods(iModSrc), Mods(iModDst), Turbine, ErrStat2, ErrMsg2)
         case (Module_Orca)
            call InitMappings_Orca(Mappings, Mods(iModSrc), Mods(iModDst), Turbine, ErrStat2, ErrMsg2)
         case (Module_SD)
            call InitMappings_SD(Mappings, Mods(iModSrc), Mods(iModDst), Turbine, ErrStat2, ErrMsg2)
         case (Module_SeaSt)
            ! call InitMappings_SeaSt(Mappings, Mods(iModSrc), Mods(iModDst), Turbine, ErrStat2, ErrMsg2)
         case (Module_SrvD)
            call InitMappings_SrvD(Mappings, Mods(iModSrc), Mods(iModDst), Turbine, ErrStat2, ErrMsg2)
         end select
         if (Failed()) return
      end do
   end do

   !----------------------------------------------------------------------------
   ! Get module indices in ModData and determine which mappings are active
   !----------------------------------------------------------------------------

   ! Reorder the mappings so that motion maps come before the load maps
   Mappings = [pack(Mappings, Mappings%MapType == Map_MotionMesh), &
               pack(Mappings, Mappings%MapType == Map_LoadMesh), &
               pack(Mappings, Mappings%MapType == Map_Variable)]

   ! Loop through mappings
   do iMap = 1, size(Mappings)
      associate (SrcMod => Mods(Mappings(iMap)%SrcModIdx), &
                 DstMod => Mods(Mappings(iMap)%DstModIdx))

         ! Add mapping index to sorce and destination module mapping arrays
         SrcMod%SrcMaps = [SrcMod%SrcMaps, iMap]
         DstMod%DstMaps = [DstMod%DstMaps, iMap]

         write (*, *) "Mapping: ", Mappings(iMap)%Desc
         ! write (*, *) " Src: ", trim(SrcMod%Abbr), SrcMod%Ins
         ! if (Mappings(iMap)%iLocSrcTransDisp(1) > 0) write (*, *) "   iLocTransDisp       ", Mappings(iMap)%iLocSrcTransDisp + SrcMod%iyg - 1
         ! if (Mappings(iMap)%iLocSrcTransVel(1) > 0) write (*, *) "   iLocTransVel        ", Mappings(iMap)%iLocSrcTransVel + SrcMod%iyg - 1
         ! if (Mappings(iMap)%iLocSrcTransAcc(1) > 0) write (*, *) "   iLocTransAcc        ", Mappings(iMap)%iLocSrcTransAcc + SrcMod%iyg - 1
         ! if (Mappings(iMap)%iLocSrcOrientation(1) > 0) write (*, *) "   iLocOrientation     ", Mappings(iMap)%iLocSrcOrientation + SrcMod%iyg - 1
         ! if (Mappings(iMap)%iLocSrcAngularVel(1) > 0) write (*, *) "   iLocAngularVel      ", Mappings(iMap)%iLocSrcAngularVel + SrcMod%iyg - 1
         ! if (Mappings(iMap)%iLocSrcAngularAcc(1) > 0) write (*, *) "   iLocAngularAcc      ", Mappings(iMap)%iLocSrcAngularAcc + SrcMod%iyg - 1
         ! if (Mappings(iMap)%iLocSrcForce(1) > 0) write (*, *) "   iLocForce           ", Mappings(iMap)%iLocSrcForce + SrcMod%iyg - 1
         ! if (Mappings(iMap)%iLocSrcMoment(1) > 0) write (*, *) "   iLocMoment          ", Mappings(iMap)%iLocSrcMoment + SrcMod%iyg - 1
         ! if (Mappings(iMap)%iLocSrcDispTransDisp(1) > 0) write (*, *) "   iLocDispTransDisp   ", Mappings(iMap)%iLocSrcDispTransDisp + SrcMod%iyg - 1
         ! write (*, *) " Dst: ", trim(DstMod%Abbr), DstMod%Ins
         ! if (Mappings(iMap)%iLocDstTransDisp(1) > 0) write (*, *) "   iLocTransDisp       ", Mappings(iMap)%iLocDstTransDisp
         ! if (Mappings(iMap)%iLocDstTransVel(1) > 0) write (*, *) "   iLocTransVel        ", Mappings(iMap)%iLocDstTransVel
         ! if (Mappings(iMap)%iLocDstTransAcc(1) > 0) write (*, *) "   iLocTransAcc        ", Mappings(iMap)%iLocDstTransAcc
         ! if (Mappings(iMap)%iLocDstOrientation(1) > 0) write (*, *) "   iLocOrientation     ", Mappings(iMap)%iLocDstOrientation
         ! if (Mappings(iMap)%iLocDstAngularVel(1) > 0) write (*, *) "   iLocAngularVel      ", Mappings(iMap)%iLocDstAngularVel
         ! if (Mappings(iMap)%iLocDstAngularAcc(1) > 0) write (*, *) "   iLocAngularAcc      ", Mappings(iMap)%iLocDstAngularAcc
         ! if (Mappings(iMap)%iLocDstForce(1) > 0) write (*, *) "   iLocForce           ", Mappings(iMap)%iLocDstForce
         ! if (Mappings(iMap)%iLocDstMoment(1) > 0) write (*, *) "   iLocMoment          ", Mappings(iMap)%iLocDstMoment
         ! if (Mappings(iMap)%iLocDstDispTransDisp(1) > 0) write (*, *) "   iLocDispTransDisp   ", Mappings(iMap)%iLocDstDispTransDisp
         ! if (Mappings(iMap)%iLocDstDispOrientation(1) > 0) write (*, *) "   iLocDispOrientation ", Mappings(iMap)%iLocDstDispOrientation

      end associate
   end do

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine InitMappings_AD(Mappings, SrcMod, DstMod, Turbine, ErrStat, ErrMsg)
   type(TC_MappingType), allocatable      :: Mappings(:)
   type(ModDataType), intent(in)          :: SrcMod, DstMod
   type(FAST_TurbineType), intent(inout)  :: Turbine
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg

   character(*), parameter    :: RoutineName = 'InitMappings_AD'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2
   integer(IntKi)             :: i

   ErrStat = ErrID_None
   ErrMsg = ''

   select case (SrcMod%ID)

   case (Module_BD)

      call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                         SrcMeshLoc=MeshLocType(BD_y_BldMotion), &                                 ! BD%y(SrcMod%Ins)%BldMotion
                         DstMeshLoc=MeshLocType(AD_u_rotors_BladeMotion, AD_rotor, SrcMod%Ins), &  ! AD%u%rotors(1)%BladeMotion(SrcMod%Ins)
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

   case (Module_ED)

      if (Turbine%p_FAST%CompElast == Module_ED) then
         do i = 1, size(Turbine%ED%y%BladeLn2Mesh)
            call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                               SrcMeshLoc=MeshLocType(ED_y_BladeLn2Mesh, i), &                  ! ED%y%BladeLn2Mesh(i)
                               DstMeshLoc=MeshLocType(AD_u_rotors_BladeMotion, AD_rotor, i), &  ! AD%u%rotors(1)%BladeMotion(i)
                               ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return
         end do
      end if

      call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                         SrcMeshLoc=MeshLocType(ED_y_TowerLn2Mesh), &                     ! ED%y%TowerLn2Mesh
                         DstMeshLoc=MeshLocType(AD_u_rotors_TowerMotion, AD_rotor), &     ! AD%u%rotors(1)%TowerMotion
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

      do i = 1, size(Turbine%ED%y%BladeRootMotion)
         call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                            SrcMeshLoc=MeshLocType(ED_y_BladeRootMotion, i), &                  ! ED%y%BladeRootMotion(i)
                            DstMeshLoc=MeshLocType(AD_u_rotors_BladeRootMotion, AD_rotor, i), & ! AD%u%rotors(1)%BladeRootMotion(i)
                            ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return
      end do

      call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                         SrcMeshLoc=MeshLocType(ED_y_HubPtMotion), &                      ! ED%y%HubPtMotion
                         DstMeshLoc=MeshLocType(AD_u_rotors_HubMotion, AD_rotor), &       ! AD%u%rotors(1)%HubMotion
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

      call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                         SrcMeshLoc=MeshLocType(ED_y_NacelleMotion), &                    ! ED%y%NacelleMotion
                         DstMeshLoc=MeshLocType(AD_u_rotors_NacelleMotion, AD_rotor), &   ! AD%u%rotors(1)%NacelleMotion
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

      call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                         SrcMeshLoc=MeshLocType(ED_y_TFinCMMotion), &                     ! ED%y%TFinCMMotion
                         DstMeshLoc=MeshLocType(AD_u_rotors_TFinMotion, AD_rotor), &      ! AD%u%rotors(1)%TFinMotion
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

   case (Module_IfW)

      call MapVariable(Mappings, "IfW HWindSpeed -> AD HWindSpeed", &
                       SrcMod=SrcMod, DstMod=DstMod, &
                       iVarSrc=Turbine%IfW%p%iVarHWindSpeed, &
                       iVarDst=Turbine%AD%p%rotors(DstMod%Ins)%iVarHWindSpeed, &
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2); if (Failed()) return

      call MapVariable(Mappings, "IfW PLExp -> AD PLExp", &
                       SrcMod=SrcMod, DstMod=DstMod, &
                       iVarSrc=Turbine%IfW%p%iVarPLExp, &
                       iVarDst=Turbine%AD%p%rotors(DstMod%Ins)%iVarPLExp, &
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2); if (Failed()) return

      call MapVariable(Mappings, "IfW PropagationDir -> AD PropagationDir", &
                       SrcMod=SrcMod, DstMod=DstMod, &
                       iVarSrc=Turbine%IfW%p%iVarPropagationDir, &
                       iVarDst=Turbine%AD%p%rotors(DstMod%Ins)%iVarPropagationDir, &
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2); if (Failed()) return

   case (Module_SrvD)

      call MapCustom(Mappings, Custom_SrvD_to_AD, SrcMod, DstMod)

   end select

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine InitMappings_BD(Mappings, SrcMod, DstMod, Turbine, ErrStat, ErrMsg)
   type(TC_MappingType), allocatable      :: Mappings(:)
   type(ModDataType), intent(in)          :: SrcMod, DstMod
   type(FAST_TurbineType), intent(inout)  :: Turbine           !< Turbine type
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg

   character(*), parameter    :: RoutineName = 'InitMappings_BD'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2
   integer(IntKi)             :: i

   ErrStat = ErrID_None
   ErrMsg = ''

   select case (SrcMod%ID)

   case (Module_AD)

      call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                       SrcMeshLoc=MeshLocType(AD_y_rotors_BladeLoad, AD_rotor, DstMod%Ins), &         ! AD%y%rotors(1)%BladeLoad(DstMod%Ins)
                       SrcDispMeshLoc=MeshLocType(AD_u_rotors_BladeMotion, AD_rotor, DstMod%Ins), &   ! AD%u%rotors(1)%BladeMotion(DstMod%Ins)
                       DstMeshLoc=MeshLocType(BD_u_DistrLoad), &                                      ! BD%u(DstMod%Ins)%DistrLoad
                       DstDispMeshLoc=MeshLocType(BD_y_BldMotion), &                                  ! BD%y(DstMod%Ins)%BldMotion
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

   case (Module_ED)

      call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                         SrcMeshLoc=MeshLocType(ED_y_BladeRootMotion, DstMod%Ins), &   ! ED%y%BladeRootMotion(DstMod%Ins)
                         DstMeshLoc=MeshLocType(BD_u_RootMotion), &                    ! BD%u(DstMod%Ins)%RootMotion
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

      ! Hub motion not used
      ! call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
      !                    SrcMeshLoc=MeshLocType(ED_y_HubPtMotion), &                   ! ED%y%HubED_y_HubPtMotion
      !                    DstMeshLoc=MeshLocType(BD_u_HubMotion), &                     ! BD%Input(1, DstMod%Ins)%HubMotion
      !                    ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

   case (Module_ExtLd)

      ! TODO
      ! CALL MeshMapCreate( ExtLd%y%BladeLoad(K), BD%Input(1,k)%DistrLoad,  MeshMapData%ExtLd_P_2_BDED_B(K), ErrStat2, ErrMsg2 )

   case (Module_SrvD)

      do i = 1, Turbine%SrvD%p%NumBStC
         call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                          SrcMeshLoc=MeshLocType(SrvD_y_BStCLoadMesh, DstMod%Ins, i), &        ! SrvD%y%BStCLoadMesh(DstMod%Ins, i), &
                          SrcDispMeshLoc=MeshLocType(SrvD_u_BStCMotionMesh, DstMod%Ins, i), &  ! SrvD%u%BStCMotionMesh(DstMod%Ins, i)
                          DstMeshLoc=MeshLocType(BD_u_DistrLoad), &                            ! BD%Input(1, DstMod%Ins)%DistrLoad
                          DstDispMeshLoc=MeshLocType(BD_y_BldMotion), &                        ! BD%y(DstMod%Ins)%BldMotion
                          ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return
      end do

   end select

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine InitMappings_ED(Mappings, SrcMod, DstMod, Turbine, ErrStat, ErrMsg)
   type(TC_MappingType), allocatable      :: Mappings(:)
   type(ModDataType), intent(in)          :: SrcMod, DstMod
   type(FAST_TurbineType), intent(inout)  :: Turbine           !< Turbine type
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg

   character(*), parameter    :: RoutineName = 'InitMappings_ED'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2
   integer(IntKi)             :: i, j

   ErrStat = ErrID_None
   ErrMsg = ''

   select case (SrcMod%ID)

   case (Module_AD)

      do i = 1, Turbine%ED%p%NumBl
         call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                          SrcMeshLoc=MeshLocType(AD_y_rotors_BladeLoad, AD_rotor, i), &         ! AD%y%rotors(iR)%BladeLoad(i)
                          SrcDispMeshLoc=MeshLocType(AD_u_rotors_BladeMotion, AD_rotor, i), &   ! AD%u%rotors(iR)%BladeMotion(i)
                          DstMeshLoc=MeshLocType(ED_u_BladePtLoads, i), &                       ! ED%u%BladePtLoads(i)
                          DstDispMeshLoc=MeshLocType(ED_y_BladeLn2Mesh, i), &                   ! ED%y%BladeLn2Mesh(i)
                          Active=Turbine%p_FAST%CompElast == Module_ED, &
                          ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return
      end do

      call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                       SrcMeshLoc=MeshLocType(AD_y_rotors_HubLoad, AD_Rotor), &           ! AD%y%rotors(1)%HubLoad
                       SrcDispMeshLoc=MeshLocType(AD_u_rotors_HubMotion, AD_rotor), &     ! AD%u%rotors(1)%HubMotion
                       DstMeshLoc=MeshLocType(ED_u_HubPtLoad), &                          ! ED%u%HubPtLoad
                       DstDispMeshLoc=MeshLocType(ED_y_HubPtMotion), &                    ! ED%y%HubPtMotion
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

      call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                       SrcMeshLoc=MeshLocType(AD_y_rotors_NacelleLoad, AD_Rotor), &       ! AD%y%rotors(1)%NacelleLoad
                       SrcDispMeshLoc=MeshLocType(AD_u_rotors_NacelleMotion, AD_rotor), & ! AD%u%rotors(1)%NacelleMotion
                       DstMeshLoc=MeshLocType(ED_u_NacelleLoads), &                       ! ED%u%NacelleLoads
                       DstDispMeshLoc=MeshLocType(ED_y_NacelleMotion), &                  ! ED%y%NacelleMotion
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

      call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                       SrcMeshLoc=MeshLocType(AD_y_rotors_TFinLoad, AD_Rotor), &          ! AD%y%rotors(1)%TFinLoad
                       SrcDispMeshLoc=MeshLocType(AD_u_rotors_TFinMotion, AD_rotor), &    ! AD%u%rotors(1)%TFinMotion
                       DstMeshLoc=MeshLocType(ED_u_TFinCMLoads), &                        ! ED%u%TFinCMLoads
                       DstDispMeshLoc=MeshLocType(ED_y_TFinCMMotion), &                   ! ED%y%TFinCMMotion
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

      call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                       SrcMeshLoc=MeshLocType(AD_y_rotors_TowerLoad, AD_Rotor), &         ! AD%y%rotors(1)%TowerLoad
                       SrcDispMeshLoc=MeshLocType(AD_u_rotors_TowerMotion, AD_rotor), &   ! AD%u%rotors(1)%TowerMotion
                       DstMeshLoc=MeshLocType(ED_u_TowerPtLoads), &                       ! ED%Input(1)%TowerPtLoads
                       DstDispMeshLoc=MeshLocType(ED_y_TowerLn2Mesh), &                   ! ED%y%TowerLn2Mesh
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

   case (Module_BD)

      call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                       SrcMeshLoc=MeshLocType(BD_y_ReactionForce), &                      ! BD%y(SrcMod%Ins)%ReactionForce
                       SrcDispMeshLoc=MeshLocType(BD_u_RootMotion), &                     ! BD%Input(1, SrcMod%Ins)%RootMotion
                       DstMeshLoc=MeshLocType(ED_u_HubPtLoad), &                          ! ED%Input(1)%HubPtLoad
                       DstDispMeshLoc=MeshLocType(ED_y_HubPtMotion), &                    ! ED%y%HubPtMotion
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

   case (Module_ExtLd)

      ! TODO
      ! CALL MeshMapCreate( ExtLd%y%TowerLoad, ED%Input(1)%TowerPtLoads,  MeshMapData%ExtLd_P_2_ED_P_T, ErrStat2, ErrMsg2 )

   case (Module_ExtPtfm)

      call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                       SrcMeshLoc=MeshLocType(ExtPtfm_y_PtfmMesh), &             ! ExtPtfm%y%PtfmMesh
                       SrcDispMeshLoc=MeshLocType(ExtPtfm_u_PtfmMesh), &         ! ExtPtfm%u%PtfmMesh
                       DstMeshLoc=MeshLocType(ED_u_PlatformPtMesh), &            ! ED%Input(1)%PlatformPtMesh
                       DstDispMeshLoc=MeshLocType(ED_y_PlatformPtMesh), &        ! ED%y%PlatformPtMesh
                       Active=Turbine%p_FAST%CompSub /= Module_SD, &
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

   case (Module_FEAM)

      call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                       SrcMeshLoc=MeshLocType(FEAM_y_PtFairleadLoad), &             ! FEAM%y%PtFairleadLoad, &
                       SrcDispMeshLoc=MeshLocType(FEAM_u_PtFairleadDisplacement), & ! FEAM%u%PtFairleadDisplacement
                       DstMeshLoc=MeshLocType(ED_u_PlatformPtMesh), &               ! ED%u%PlatformPtMesh
                       DstDispMeshLoc=MeshLocType(ED_y_PlatformPtMesh), &           ! ED%y%PlatformPtMesh
                       Active=Turbine%p_FAST%CompSub /= Module_SD, &
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

   case (Module_HD)

      ! Platform loads (SubDyn not active)
      call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                       SrcMeshLoc=MeshLocType(HydroDyn_y_Morison_Mesh), &        ! HD%y%Morison%Mesh
                       SrcDispMeshLoc=MeshLocType(HydroDyn_u_Morison_Mesh), &    ! HD%u%Morison%Mesh
                       DstMeshLoc=MeshLocType(ED_u_PlatformPtMesh), &            ! ED%u%PlatformPtMesh
                       DstDispMeshLoc=MeshLocType(ED_y_PlatformPtMesh), &        ! ED%y%PlatformPtMesh
                       Active=Turbine%p_FAST%CompSub == Module_None, &
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

      ! Platform loads (SubDyn not active)
      call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                       SrcMeshLoc=MeshLocType(HydroDyn_y_WAMITMesh), &           ! HD%y%WAMITMesh
                       SrcDispMeshLoc=MeshLocType(HydroDyn_u_WAMITMesh), &       ! HD%u%WAMITMesh
                       DstMeshLoc=MeshLocType(ED_u_PlatformPtMesh), &            ! ED%u%PlatformPtMesh
                       DstDispMeshLoc=MeshLocType(ED_y_PlatformPtMesh), &        ! ED%y%PlatformPtMesh
                       Active=Turbine%p_FAST%CompSub == Module_None, &
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

   case (Module_IceD)

      ! Platform loads (SubDyn not active)
      call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                       SrcMeshLoc=MeshLocType(IceD_y_PointMesh), &               ! IceD%y%PointMesh
                       SrcDispMeshLoc=MeshLocType(IceD_u_PointMesh), &           ! IceD%u%PointMesh
                       DstMeshLoc=MeshLocType(ED_u_PlatformPtMesh), &            ! ED%Input(1)%PlatformPtMesh
                       DstDispMeshLoc=MeshLocType(ED_y_PlatformPtMesh), &        ! ED%y%PlatformPtMesh
                       Active=Turbine%p_FAST%CompSub /= Module_SD, &
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

   case (Module_IceF)

      ! Platform loads (SubDyn not active)
      call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                       SrcMeshLoc=MeshLocType(IceFloe_y_iceMesh), &              ! IceFloe%y%iceMesh
                       SrcDispMeshLoc=MeshLocType(IceFloe_u_iceMesh), &          ! IceFloe%u%iceMesh
                       DstMeshLoc=MeshLocType(ED_u_PlatformPtMesh), &            ! ED%Input(1)%PlatformPtMesh
                       DstDispMeshLoc=MeshLocType(ED_y_PlatformPtMesh), &        ! ED%y%PlatformPtMesh
                       Active=Turbine%p_FAST%CompSub /= Module_SD, &
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

   case (Module_MAP)

      ! Platform loads (SubDyn not active)
      call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                       SrcMeshLoc=MeshLocType(MAP_y_ptFairleadLoad), &           ! MAP%y%PtFairleadLoad
                       SrcDispMeshLoc=MeshLocType(MAP_u_PtFairDisplacement), &   ! MAP%u%PtFairDisplacement
                       DstMeshLoc=MeshLocType(ED_u_PlatformPtMesh), &            ! ED%Input(1)%PlatformPtMesh
                       DstDispMeshLoc=MeshLocType(ED_y_PlatformPtMesh), &        ! ED%y%PlatformPtMesh
                       Active=Turbine%p_FAST%CompSub /= Module_SD, &
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

   case (Module_MD)

      ! Platform loads (SubDyn not active)
      call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                       SrcMeshLoc=MeshLocType(MD_y_CoupledLoads, 1), &           ! MD%y%CoupledLoads(1)
                       SrcDispMeshLoc=MeshLocType(MD_u_CoupledKinematics, 1), &  ! MD%u%CoupledKinematics(1)
                       DstMeshLoc=MeshLocType(ED_u_PlatformPtMesh), &            ! ED%Input(1)%PlatformPtMesh
                       DstDispMeshLoc=MeshLocType(ED_y_PlatformPtMesh), &        ! ED%y%PlatformPtMesh
                       Active=Turbine%p_FAST%CompSub /= Module_SD, &
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

   case (Module_Orca)

      ! Platform loads (SubDyn not active)
      call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                       SrcMeshLoc=MeshLocType(Orca_y_PtfmMesh), &                ! Orca%y%PtfmMesh
                       SrcDispMeshLoc=MeshLocType(Orca_u_PtfmMesh), &            ! Orca%u%PtfmMesh
                       DstMeshLoc=MeshLocType(ED_u_PlatformPtMesh), &            ! ED%Input(1)%PlatformPtMesh
                       DstDispMeshLoc=MeshLocType(ED_y_PlatformPtMesh), &        ! ED%y%PlatformPtMesh
                       Active=Turbine%p_FAST%CompSub /= Module_SD, &
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

   case (Module_SD)

      call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                       SrcMeshLoc=MeshLocType(SD_y_Y1Mesh), &                    ! SD%y%Y1mesh, &
                       SrcDispMeshLoc=MeshLocType(SD_u_TPMesh), &                ! SD%Input(1)%TPMesh
                       DstMeshLoc=MeshLocType(ED_u_PlatformPtMesh), &            ! ED%Input(1)%PlatformPtMesh
                       DstDispMeshLoc=MeshLocType(ED_y_PlatformPtMesh), &        ! ED%y%PlatformPtMesh
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

   case (Module_SrvD)

      call MapVariable(Mappings, "SrvD BlPitchCom -> ED BlPitchCom", &
                       SrcMod=SrcMod, iVarSrc=Turbine%SrvD%p%iVarBlPitchCom, &
                       DstMod=DstMod, iVarDst=Turbine%ED%p%iVarBlPitchCom, &
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2); if (Failed()) return

      call MapVariable(Mappings, "SrvD YawMom -> ED YawMom", &
                       SrcMod=SrcMod, iVarSrc=Turbine%SrvD%p%iVarYawMom, &
                       DstMod=DstMod, iVarDst=Turbine%ED%p%iVarYawMom, &
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2); if (Failed()) return

      call MapVariable(Mappings, "SrvD GenTrq -> ED GenTrq", &
                       SrcMod=SrcMod, iVarSrc=Turbine%SrvD%p%iVarGenTrq, &
                       DstMod=DstMod, iVarDst=Turbine%ED%p%iVarGenTrq, &
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2); if (Failed()) return

      ! Blade Structural Controller (if ElastoDyn is used for blades)
      do j = 1, Turbine%SrvD%p%NumBStC
         do i = 1, Turbine%ED%p%NumBl
            call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                             SrcMeshLoc=MeshLocType(SrvD_y_BStCLoadMesh, i, j), &        ! SrvD%y%BStCLoadMesh(i, j), &
                             SrcDispMeshLoc=MeshLocType(SrvD_u_BStCMotionMesh, i, j), &  ! SrvD%u%BStCMotionMesh(i, j)
                             DstMeshLoc=MeshLocType(ED_u_BladePtLoads, j), &             ! ED%Input(1)%BladePtLoads(j)
                             DstDispMeshLoc=MeshLocType(ED_y_BladeLn2Mesh, j), &         ! ED%y%BladeLn2Mesh(j)
                             Active=Turbine%p_FAST%CompElast == Module_ED, &
                             ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return
         end do
      end do

      ! Nacelle Structural Controller
      do j = 1, Turbine%SrvD%p%NumNStC
         call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                          SrcMeshLoc=MeshLocType(SrvD_y_NStCLoadMesh, j), &         ! SrvD%y%NStCLoadMesh(j), &
                          SrcDispMeshLoc=MeshLocType(SrvD_u_NStCMotionMesh, j), &   ! SrvD%u%NStCMotionMesh(j)
                          DstMeshLoc=MeshLocType(ED_u_NacelleLoads), &              ! ED%Input(1)%NacelleLoads
                          DstDispMeshLoc=MeshLocType(ED_y_NacelleMotion), &         ! ED%y%NacelleMotion
                          ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return
      end do

      ! Tower Structural Controller
      do j = 1, Turbine%SrvD%p%NumTStC
         call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                          SrcMeshLoc=MeshLocType(SrvD_y_TStCLoadMesh, j), &         ! SrvD%y%TStCLoadMesh(j), &
                          SrcDispMeshLoc=MeshLocType(SrvD_u_TStCMotionMesh, j), &   ! SrvD%u%TStCMotionMesh(j)
                          DstMeshLoc=MeshLocType(ED_u_TowerPtLoads), &              ! ED%Input(1)%TowerLoads
                          DstDispMeshLoc=MeshLocType(ED_y_TowerLn2Mesh), &          ! ED%y%TowerLn2Mesh
                          ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return
      end do

      ! Substructure Structural Controller
      do j = 1, Turbine%SrvD%p%NumSStC
         call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                          SrcMeshLoc=MeshLocType(SrvD_y_SStCLoadMesh, j), &         ! SrvD%y%SStCLoadMesh(j), &
                          SrcDispMeshLoc=MeshLocType(SrvD_u_SStCMotionMesh, j), &   ! SrvD%u%SStCMotionMesh(j)
                          DstMeshLoc=MeshLocType(ED_u_PlatformPtMesh), &            ! ED%Input(1)%PlatformPtMesh
                          DstDispMeshLoc=MeshLocType(ED_y_PlatformPtMesh), &        ! ED%y%PlatformPtMesh
                          Active=Turbine%p_FAST%CompSub /= Module_SD, &
                          ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return
      end do

   end select

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine InitMappings_ExtLd(Mappings, SrcMod, DstMod, Turbine, ErrStat, ErrMsg)
   type(TC_MappingType), allocatable      :: Mappings(:)
   type(ModDataType), intent(in)          :: SrcMod, DstMod
   type(FAST_TurbineType), intent(inout)  :: Turbine           !< Turbine type
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg

   character(*), parameter    :: RoutineName = 'InitMappings_ExtLd'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2
   integer(IntKi)             :: i

   ErrStat = ErrID_None
   ErrMsg = ''

   select case (SrcMod%ID)
   case (Module_AD)

      ! call MapCustom(Mappings, Custom_AD_to_ExtLd, SrcMod, DstMod, ErrStat2, ErrMsg2); if (Failed()) return

      ! TODO Add mapping from aerodyn blade and tower to new input meshes
      ! MeshMapCreate( AD%y%rotors(1)%BladeLoad(k), ExtLd%y%BladeLoadAD(k), MeshMapData%AD_L_2_ExtLd_B(k), ErrStat2, ErrMsg2)
      ! MeshMapCreate( AD%y%rotors(1)%TowerLoad, ExtLd%y%TowerLoadAD,  MeshMapData%AD_L_2_ExtLd_T, ErrStat2, ErrMsg2 )

   case (Module_BD)

      call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                         SrcMeshLoc=MeshLocType(BD_y_BldMotion), &                   ! BD%y(SrcMod%Ins)%BldMotion
                         DstMeshLoc=MeshLocType(ExtLd_u_BladeMotion, SrcMod%Ins), &  ! ExtLd%u%BladeMotion(SrcMod%Ins)
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

   case (Module_ED)

      call MapCustom(Mappings, Custom_ED_to_ExtLd, SrcMod, DstMod)

      do i = 1, Turbine%ED%p%NumBl
         call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                            SrcMeshLoc=MeshLocType(ED_y_BladeLn2Mesh, i), &        ! ED%y%BladeLn2Mesh(i)
                            DstMeshLoc=MeshLocType(ExtLd_u_BladeMotion, i), &      ! ExtLd%u%BladeMotion(i)
                            Active=Turbine%p_FAST%CompElast == Module_ED, &
                            ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return
      end do

      do i = 1, Turbine%ED%p%NumBl
         call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                            SrcMeshLoc=MeshLocType(ED_y_BladeRootMotion, i), &     ! ED%y%BladeRootMotion(i)
                            DstMeshLoc=MeshLocType(ExtLd_u_BladeRootMotion, i), &  ! ExtLd%u%BladeRootMotion(i)
                            ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return
      end do

      call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                         SrcMeshLoc=MeshLocType(ED_y_TowerLn2Mesh), &              ! ED%y%TowerLn2Mesh
                         DstMeshLoc=MeshLocType(ExtLd_u_TowerMotion), &            ! ExtLd%u%TowerMotion
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

      call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                         SrcMeshLoc=MeshLocType(ED_y_HubPtMotion), &               ! ED%y%HubPtMotion
                         DstMeshLoc=MeshLocType(ExtLd_u_HubMotion), &              ! ExtLd%u%HubMotion
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

      call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                         SrcMeshLoc=MeshLocType(ED_y_NacelleMotion), &             ! ED%y%NacelleMotion
                         DstMeshLoc=MeshLocType(ExtLd_u_NacelleMotion), &          ! ExtLd%u%NacelleMotion
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

   end select

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine InitMappings_ExtPtfm(Mappings, SrcMod, DstMod, Turbine, ErrStat, ErrMsg)
   type(TC_MappingType), allocatable      :: Mappings(:)
   type(ModDataType), intent(in)          :: SrcMod, DstMod
   type(FAST_TurbineType), intent(inout)  :: Turbine           !< Turbine type
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg

   character(*), parameter    :: RoutineName = 'InitMappings_ExtPtfm'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2
   integer(IntKi)             :: i

   ErrStat = ErrID_None
   ErrMsg = ''

   select case (SrcMod%ID)

   case (Module_ED)

      if (Turbine%p_FAST%CompSub /= Module_SD) then
         call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                            SrcMeshLoc=MeshLocType(ED_y_PlatformPtMesh), &       ! ED%y%PlatformPtMesh
                            DstMeshLoc=MeshLocType(ExtPtfm_u_PtfmMesh), &        ! ExtPtfm%u%PtfmMesh
                            ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return
      end if

   case (Module_SD)

      call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                         SrcMeshLoc=MeshLocType(ED_y_PlatformPtMesh), &           ! ED%y%PlatformPtMesh
                         DstMeshLoc=MeshLocType(SD_u_TPMesh), &                   ! SD%u%TPMesh
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

   end select

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine InitMappings_FEAM(Mappings, SrcMod, DstMod, Turbine, ErrStat, ErrMsg)
   type(TC_MappingType), allocatable      :: Mappings(:)
   type(ModDataType), intent(in)          :: SrcMod, DstMod
   type(FAST_TurbineType), intent(inout)  :: Turbine           !< Turbine type
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg

   character(*), parameter    :: RoutineName = 'InitMappings_FEAM'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2
   integer(IntKi)             :: i

   ErrStat = ErrID_None
   ErrMsg = ''

   select case (SrcMod%ID)

   case (Module_ED)

      call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                         SrcMeshLoc=MeshLocType(ED_y_PlatformPtMesh), &             ! ED%y%PlatformPtMesh
                         DstMeshLoc=MeshLocType(FEAM_u_PtFairleadDisplacement), &   ! FEAM%u%PtFairleadDisplacement
                         Active=Turbine%p_FAST%CompSub /= Module_SD, &
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

   case (Module_SD)

      call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                         SrcMeshLoc=MeshLocType(SD_y_Y3Mesh), &                     ! SD%y%y3Mesh
                         DstMeshLoc=MeshLocType(FEAM_u_PtFairleadDisplacement), &   ! FEAM%u%PtFairleadDisplacement
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

   end select

   select case (SrcMod%ID)
   case (Module_ED)

      if (Turbine%p_FAST%CompSub /= Module_SD) then
         ! CALL MeshMapCreate( SubstructureMotion, FEAM%Input(1)%PtFairleadDisplacement,  MeshMapData%Structure_2_Mooring, ErrStat2, ErrMsg2 )
      end if

   case (Module_SD)

      ! CALL MeshMapCreate( SubstructureMotion, FEAM%Input(1)%PtFairleadDisplacement,  MeshMapData%Structure_2_Mooring, ErrStat2, ErrMsg2 )

   end select

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine InitMappings_HD(Mappings, SrcMod, DstMod, Turbine, ErrStat, ErrMsg)
   type(TC_MappingType), allocatable      :: Mappings(:)
   type(ModDataType), intent(in)          :: SrcMod, DstMod
   type(FAST_TurbineType), intent(inout)  :: Turbine           !< Turbine type
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg

   character(*), parameter    :: RoutineName = 'InitMappings_HD'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2
   integer(IntKi)             :: i

   ErrStat = ErrID_None
   ErrMsg = ''

   select case (SrcMod%ID)

   case (Module_ED)

      call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                         SrcMeshLoc=MeshLocType(ED_y_PlatformPtMesh), &          ! ED%y%PlatformPtMesh
                         DstMeshLoc=MeshLocType(HydroDyn_u_PRPMesh), &           ! HD%Input(1)%PRPMesh
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

      call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                         SrcMeshLoc=MeshLocType(ED_y_PlatformPtMesh), &       ! ED%y%PlatformPtMesh
                         DstMeshLoc=MeshLocType(HydroDyn_u_Morison_Mesh), &   ! HD%Input(1)%Morison%Mesh
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2, &
                         Active=Turbine%p_FAST%CompSub /= Module_SD); if(Failed()) return

      call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                         SrcMeshLoc=MeshLocType(ED_y_PlatformPtMesh), &       ! ED%y%PlatformPtMesh
                         DstMeshLoc=MeshLocType(HydroDyn_u_WAMITMesh), &      ! HD%Input(1)%WAMITMesh
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2, &
                         Active=Turbine%p_FAST%CompSub /= Module_SD); if(Failed()) return

   case (Module_SD)

      call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                         SrcMeshLoc=MeshLocType(SD_y_Y2Mesh), &                  ! SD%y%Y2Mesh
                         DstMeshLoc=MeshLocType(HydroDyn_u_Morison_Mesh), &      ! HD%Input(1)%Morison%Mesh
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

      call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                         SrcMeshLoc=MeshLocType(SD_y_Y2Mesh), &                  ! SD%y%Y2Mesh
                         DstMeshLoc=MeshLocType(HydroDyn_u_WAMITMesh), &         ! HD%Input(1)%WAMITMesh
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return
   end select

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine InitMappings_IceD(Mappings, SrcMod, DstMod, Turbine, ErrStat, ErrMsg)
   type(TC_MappingType), allocatable      :: Mappings(:)
   type(ModDataType), intent(in)          :: SrcMod, DstMod
   type(FAST_TurbineType), intent(inout)  :: Turbine
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg

   character(*), parameter    :: RoutineName = 'InitMappings_IceD'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2

   ErrStat = ErrID_None
   ErrMsg = ''

   select case (SrcMod%ID)
   case (Module_ED)

      call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                         SrcMeshLoc=MeshLocType(ED_y_PlatformPtMesh), &          ! ED%y%PlatformPtMesh
                         DstMeshLoc=MeshLocType(IceD_u_PointMesh), &             ! IceD%u%PointMesh
                         Active=Turbine%p_FAST%CompSub /= Module_SD, &
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

   case (Module_SD)

      call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                         SrcMeshLoc=MeshLocType(SD_y_Y3Mesh), &                  ! SD%y%y3Mesh
                         DstMeshLoc=MeshLocType(IceD_u_PointMesh), &             ! IceD%u%PointMesh
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

   end select

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine InitMappings_IceF(Mappings, SrcMod, DstMod, Turbine, ErrStat, ErrMsg)
   type(TC_MappingType), allocatable      :: Mappings(:)
   type(ModDataType), intent(in)          :: SrcMod, DstMod
   type(FAST_TurbineType), intent(inout)  :: Turbine
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg

   character(*), parameter    :: RoutineName = 'InitMappings_IceF'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2

   ErrStat = ErrID_None
   ErrMsg = ''

   select case (SrcMod%ID)
   case (Module_ED)

      call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                         SrcMeshLoc=MeshLocType(ED_y_PlatformPtMesh), &          ! ED%y%PlatformPtMesh
                         DstMeshLoc=MeshLocType(IceFloe_u_iceMesh), &            ! IceFloe%u%iceMesh
                         Active=Turbine%p_FAST%CompSub /= Module_SD, &
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

   case (Module_SD)

      call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                         SrcMeshLoc=MeshLocType(SD_y_Y3Mesh), &                  ! SD%y%y3Mesh
                         DstMeshLoc=MeshLocType(IceFloe_u_iceMesh), &            ! IceFloe%u%iceMesh
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

   end select

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine InitMappings_IfW(Mappings, SrcMod, DstMod, T, ErrStat, ErrMsg)
   type(TC_MappingType), allocatable      :: Mappings(:)
   type(ModDataType), intent(in)          :: SrcMod, DstMod
   type(FAST_TurbineType), intent(inout)  :: T           !< Turbine type
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg

   character(*), parameter    :: RoutineName = 'InitMappings_IfW'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2

   ErrStat = ErrID_None
   ErrMsg = ''

   select case (SrcMod%ID)
   case (Module_ED)
      call MapCustom(Mappings, Custom_ED_to_IfW, SrcMod, DstMod)
   case (Module_SrvD)
      call MapCustom(Mappings, Custom_SrvD_to_IfW, SrcMod=SrcMod, DstMod=DstMod)
   end select

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine InitMappings_MAP(Mappings, SrcMod, DstMod, Turbine, ErrStat, ErrMsg)
   type(TC_MappingType), allocatable      :: Mappings(:)
   type(ModDataType), intent(in)          :: SrcMod, DstMod
   type(FAST_TurbineType), intent(inout)  :: Turbine           !< Turbine type
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg

   character(*), parameter    :: RoutineName = 'InitMappings_MAP'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2

   ErrStat = ErrID_None
   ErrMsg = ''

   select case (SrcMod%ID)
   case (Module_ED)

      call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                         SrcMeshLoc=MeshLocType(ED_y_PlatformPtMesh), &       ! ED%y%PlatformPtMesh
                         DstMeshLoc=MeshLocType(MAP_u_PtFairDisplacement), &  ! MAPp%u%PtFairDisplacement
                         Active=Turbine%p_FAST%CompSub /= Module_SD, &
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

   case (Module_SD)

      call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                         SrcMeshLoc=MeshLocType(SD_y_Y3Mesh), &                  ! SD%y%y3Mesh
                         DstMeshLoc=MeshLocType(MAP_u_PtFairDisplacement), &     ! MAPp%u%PtFairDisplacement
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

   end select

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine InitMappings_MD(Mappings, SrcMod, DstMod, Turbine, ErrStat, ErrMsg)
   type(TC_MappingType), allocatable      :: Mappings(:)
   type(ModDataType), intent(in)          :: SrcMod, DstMod
   type(FAST_TurbineType), intent(inout)  :: Turbine           !< Turbine type
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg

   character(*), parameter    :: RoutineName = 'InitMappings_MD'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2

   ErrStat = ErrID_None
   ErrMsg = ''

   select case (SrcMod%ID)

   case (Module_ED)

      call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                         SrcMeshLoc=MeshLocType(ED_y_PlatformPtMesh), &          ! ED%y%PlatformPtMesh
                         DstMeshLoc=MeshLocType(MD_u_CoupledKinematics, 1), &    ! MD%u%CoupledKinematics(1)
                         Active=Turbine%p_FAST%CompSub /= Module_SD, &
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

   case (Module_SD)

      call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                         SrcMeshLoc=MeshLocType(SD_y_Y3Mesh), &                  ! SD%y%y3Mesh
                         DstMeshLoc=MeshLocType(MD_u_CoupledKinematics, 1), &    ! MD%u%CoupledKinematics(1)
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

   case (Module_SrvD)

      call MapCustom(Mappings, Custom_SrvD_to_MD, SrcMod, DstMod)

   end select

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine InitMappings_Orca(Mappings, SrcMod, DstMod, Turbine, ErrStat, ErrMsg)
   type(TC_MappingType), allocatable      :: Mappings(:)
   type(ModDataType), intent(in)          :: SrcMod, DstMod
   type(FAST_TurbineType), intent(inout)  :: Turbine           !< Turbine type
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg

   character(*), parameter    :: RoutineName = 'InitMappings_Orca'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2

   ErrStat = ErrID_None
   ErrMsg = ''

   select case (SrcMod%ID)
   case (Module_ED)

      call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                         SrcMeshLoc=MeshLocType(ED_y_PlatformPtMesh), &          ! ED%y%PlatformPtMesh
                         DstMeshLoc=MeshLocType(Orca_u_PtfmMesh), &              ! Orca%u%PtfmMesh
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

   end select

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine InitMappings_SD(Mappings, SrcMod, DstMod, Turbine, ErrStat, ErrMsg)
   type(TC_MappingType), allocatable      :: Mappings(:)
   type(ModDataType), intent(in)          :: SrcMod, DstMod
   type(FAST_TurbineType), intent(inout)  :: Turbine           !< Turbine type
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg

   character(*), parameter    :: RoutineName = 'InitMappings_BD'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2
   integer(IntKi)             :: i, j

   ErrStat = ErrID_None
   ErrMsg = ''

   select case (SrcMod%ID)

   case (Module_ED)

      call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                         SrcMeshLoc=MeshLocType(ED_y_PlatformPtMesh), &          ! ED%y%PlatformPtMesh
                         DstMeshLoc=MeshLocType(SD_u_TPMesh), &                  ! SD%u%TPMesh
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

   case (Module_FEAM)

      call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                       SrcMeshLoc=MeshLocType(FEAM_y_PtFairleadLoad), &             ! FEAM%y%PtFairleadLoad, &
                       SrcDispMeshLoc=MeshLocType(FEAM_u_PtFairleadDisplacement), & ! FEAM%u%PtFairleadDisplacement
                       DstMeshLoc=MeshLocType(SD_u_LMesh), &                        ! SD%u%LMesh
                       DstDispMeshLoc=MeshLocType(SD_y_y3Mesh), &                   ! SD%y%y3Mesh
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

   case (Module_HD)

      call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                       SrcMeshLoc=MeshLocType(HydroDyn_y_Morison_Mesh), &          ! HD%y%Morison%Mesh
                       SrcDispMeshLoc=MeshLocType(HydroDyn_u_Morison_Mesh), &      ! HD%u%Morison%Mesh
                       DstMeshLoc=MeshLocType(SD_u_LMesh), &                       ! SD%u%LMesh
                       DstDispMeshLoc=MeshLocType(SD_y_y3Mesh), &                  ! SD%y%y3Mesh
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

      call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                       SrcMeshLoc=MeshLocType(HydroDyn_y_WAMITMesh), &             ! HD%y%WAMITMesh
                       SrcDispMeshLoc=MeshLocType(HydroDyn_u_WAMITMesh), &         ! HD%u%WAMITMesh
                       DstMeshLoc=MeshLocType(SD_u_LMesh), &                       ! SD%u%LMesh
                       DstDispMeshLoc=MeshLocType(SD_y_y3Mesh), &                  ! SD%y%y3Mesh
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

   case (Module_IceD)

      call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                       SrcMeshLoc=MeshLocType(IceD_y_PointMesh), &                 ! IceD%y%PointMesh
                       SrcDispMeshLoc=MeshLocType(IceD_u_PointMesh), &             ! IceD%u%PointMesh
                       DstMeshLoc=MeshLocType(SD_u_LMesh), &                       ! SD%u%LMesh
                       DstDispMeshLoc=MeshLocType(SD_y_y3Mesh), &                  ! SD%y%y3Mesh
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

   case (Module_IceF)

      call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                       SrcMeshLoc=MeshLocType(IceFloe_y_iceMesh), &                ! IceFloe%y%iceMesh
                       SrcDispMeshLoc=MeshLocType(IceFloe_u_iceMesh), &            ! IceFloe%u%iceMesh
                       DstMeshLoc=MeshLocType(SD_u_LMesh), &                       ! SD%u%LMesh
                       DstDispMeshLoc=MeshLocType(SD_y_y3Mesh), &                  ! SD%y%y3Mesh
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

   case (Module_MAP)

      call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                       SrcMeshLoc=MeshLocType(MAP_y_ptFairleadLoad), &             ! MAP%y%PtFairleadLoad
                       SrcDispMeshLoc=MeshLocType(MAP_u_PtFairDisplacement), &     ! MAP%u%PtFairDisplacement
                       DstMeshLoc=MeshLocType(SD_u_LMesh), &                       ! SD%u%LMesh
                       DstDispMeshLoc=MeshLocType(SD_y_y3Mesh), &                  ! SD%y%y3Mesh
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

   case (Module_MD)

      call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                       SrcMeshLoc=MeshLocType(MD_y_CoupledLoads, 1), &             ! MD%y%CoupledLoads(1)
                       SrcDispMeshLoc=MeshLocType(MD_u_CoupledKinematics, 1), &    ! MD%u%CoupledKinematics(1)
                       DstMeshLoc=MeshLocType(SD_u_LMesh), &                       ! SD%u%LMesh
                       DstDispMeshLoc=MeshLocType(SD_y_y3Mesh), &                  ! SD%y%y3Mesh
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

   case (Module_SrvD)

      call MapCustom(Mappings, Custom_SrvD_to_SD, SrcMod, DstMod)

      ! Substructure Structural Controller
      do j = 1, Turbine%SrvD%p%NumSStC
         call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                          SrcMeshLoc=MeshLocType(SrvD_y_SStCLoadMesh, j), &         ! SrvD%y%SStCLoadMesh(j), &
                          SrcDispMeshLoc=MeshLocType(SrvD_u_SStCMotionMesh, j), &   ! SrvD%u%SStCMotionMesh(j)
                          DstMeshLoc=MeshLocType(SD_u_LMesh), &                     ! SD%u%LMesh
                          DstDispMeshLoc=MeshLocType(SD_y_y3Mesh), &                ! SD%y%y3Mesh
                          ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return
      end do

   end select

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine InitMappings_SrvD(Mappings, SrcMod, DstMod, Turbine, ErrStat, ErrMsg)
   type(TC_MappingType), allocatable      :: Mappings(:)
   type(ModDataType), intent(in)          :: SrcMod, DstMod
   type(FAST_TurbineType), intent(inout)  :: Turbine           !< Turbine type
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg

   character(*), parameter    :: RoutineName = 'InitMappings_BD'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2
   integer(IntKi)             :: i, j

   ErrStat = ErrID_None
   ErrMsg = ''

   ! MeshMapCreate( PlatformMotion, SrvD%u%PtfmMotionMesh, MeshMapData%ED_P_2_SrvD_P_P, ErrStat2, ErrMsg2 )


   select case (SrcMod%ID)

   case (Module_BD)

      call MapCustom(Mappings, Custom_BD_to_SrvD, SrcMod, DstMod)

      ! Blade Structural Controller
      do i = 1, Turbine%SrvD%p%NumBStC
         call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                            SrcMeshLoc=MeshLocType(BD_y_BldMotion), &                         ! BD%y%BldMotion
                            DstMeshLoc=MeshLocType(SrvD_u_BStCMotionMesh, DstMod%Ins, i), &   ! SrvD%u%BStCMotionMesh(i, j)
                            ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return
      end do

   case (Module_ED)

      call MapCustom(Mappings, Custom_ED_to_SrvD, SrcMod, DstMod)

      call MapVariable(Mappings, "ED Yaw -> SrvD Yaw", &
                       SrcMod=SrcMod, iVarSrc=Turbine%ED%p%iVarYaw, &
                       DstMod=DstMod, iVarDst=Turbine%SrvD%p%iVarYaw, &
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2); if (Failed()) return

      call MapVariable(Mappings, "ED YawRate -> SrvD YawRate", &
                       SrcMod=SrcMod, iVarSrc=Turbine%ED%p%iVarYawRate, &
                       DstMod=DstMod, iVarDst=Turbine%SrvD%p%iVarYawRate, &
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2); if (Failed()) return

      call MapVariable(Mappings, "ED HSS_Spd -> SrvD HSS_Spd", &
                       SrcMod=SrcMod, iVarSrc=Turbine%ED%p%iVarHSS_Spd, &
                       DstMod=DstMod, iVarDst=Turbine%SrvD%p%iVarHSS_Spd, &
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2); if (Failed()) return

      call MapVariable(Mappings, "ED HSS_Spd -> SrvD HSS_Spd", &
                       SrcMod=SrcMod, iVarSrc=Turbine%ED%p%iVarHSS_Spd, &
                       DstMod=DstMod, iVarDst=Turbine%SrvD%p%iVarHSS_Spd, &
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2); if (Failed()) return

      ! Nacelle Structural Controller
      do j = 1, Turbine%SrvD%p%NumNStC
         call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                            SrcMeshLoc=MeshLocType(ED_y_NacelleMotion), &          ! ED%y%NacelleMotion
                            DstMeshLoc=MeshLocType(SrvD_u_NStCMotionMesh, j), &    ! SrvD%u%NStCMotionMesh(j)
                            ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return
      end do

      ! Tower Structural Controller
      do j = 1, Turbine%SrvD%p%NumTStC
         call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                            SrcMeshLoc=MeshLocType(ED_y_TowerLn2Mesh), &           ! ED%y%TowerMotion
                            DstMeshLoc=MeshLocType(SrvD_u_TStCMotionMesh, j), &    ! SrvD%u%TStCMotionMesh(j)
                            ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return
      end do

      ! Blade Structural Controller (if ElastoDyn blades)
      do j = 1, Turbine%SrvD%p%NumBStC
         do i = 1, Turbine%ED%p%NumBl
            call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                               SrcMeshLoc=MeshLocType(ED_y_BladeLn2Mesh, i), &         ! ED%y%BladeLn2Mesh(i)
                               DstMeshLoc=MeshLocType(SrvD_u_BStCMotionMesh, i, j), &  ! SrvD%u%BStCMotionMesh(i, j)
                               Active=Turbine%p_FAST%CompElast == Module_ED, &
                               ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return
         end do
      end do

      ! Substructure Structural Controller (if not using SubDyn)
      do j = 1, Turbine%SrvD%p%NumSStC
         call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                            SrcMeshLoc=MeshLocType(ED_y_PlatformPtMesh), &         ! ED%y%PlatformPtMesh
                            DstMeshLoc=MeshLocType(SrvD_u_SStCMotionMesh, j), &    ! SrvD%u%SStCMotionMesh(j)
                            Active=Turbine%p_FAST%CompSub /= Module_SD, &
                            ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return
      end do

   case (Module_IfW)

      call MapCustom(Mappings, Custom_IfW_to_SrvD, SrcMod=SrcMod, DstMod=DstMod)

   case (Module_SD)

      ! Substructure Structural Controller
      do j = 1, Turbine%SrvD%p%NumSStC
         call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                            SrcMeshLoc=MeshLocType(SD_y_y3Mesh), &               ! SD%y%y3Mesh
                            DstMeshLoc=MeshLocType(SrvD_u_SStCMotionMesh, j), &  ! SrvD%u%SStCMotionMesh(j)
                            ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return
      end do

   end select

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine MapLoadMesh(Turbine, Mappings, SrcMod, SrcMeshLoc, SrcDispMeshLoc, &
                       DstMod, DstMeshLoc, DstDispMeshLoc, ErrStat, ErrMsg, Active)
   type(FAST_TurbineType), target         :: Turbine
   type(TC_MappingType), allocatable      :: Mappings(:)
   type(ModDataType), intent(in)          :: SrcMod, DstMod
   type(MeshLocType), intent(in)          :: SrcMeshLoc, DstMeshLoc
   type(MeshLocType), intent(in)          :: SrcDispMeshLoc, DstDispMeshLoc
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg
   logical, optional, intent(in)          :: Active

   character(*), parameter                :: RoutineName = 'MapLoadMesh'
   integer(IntKi)                         :: ErrStat2
   character(ErrMsgLen)                   :: ErrMsg2
   type(TC_MappingType)                   :: Mapping
   type(MeshType), pointer                :: SrcMesh, SrcDispMesh
   type(MeshType), pointer                :: DstMesh, DstDispMesh
   type(MeshType)                         :: DstMotionMesh

   ErrStat = ErrID_None
   ErrMsg = ''

   ! If active argument is set to false, return
   if (present(Active)) then
      if (.not. Active) return
   end if

   ! Get mesh pointers
   call FAST_OutputMeshPointer(SrcMod, Turbine, SrcMeshLoc, SrcMesh, ErrStat2, ErrMsg2); if (Failed()) return
   call FAST_InputMeshPointer(SrcMod, Turbine, SrcDispMeshLoc, .false., SrcDispMesh, ErrStat2, ErrMsg2); if (Failed()) return
   call FAST_InputMeshPointer(DstMod, Turbine, DstMeshLoc, .false., DstMesh, ErrStat2, ErrMsg2); if (Failed()) return
   call FAST_OutputMeshPointer(DstMod, Turbine, DstDispMeshLoc, DstDispMesh, ErrStat2, ErrMsg2); if (Failed()) return

   ! If any meshes aren't commited, return
   if (.not. (SrcMesh%committed .and. DstMesh%committed .and. SrcDispMesh%committed .and. DstDispMesh%committed)) return

   ! Check that all meshes in mapping have nonzero identifiers
   if (SrcMesh%ID == 0) then
      call SetErrStat(ErrID_Fatal, 'SrcMesh "'//trim(FAST_OutputMeshName(SrcMod, SrcMeshLoc))//'" not in module variables', &
                      ErrStat, ErrMsg, RoutineName)
      return
   else if (SrcDispMesh%ID == 0) then
      call SetErrStat(ErrID_Fatal, 'SrcDispMesh "'//trim(FAST_InputMeshName(SrcMod, SrcDispMeshLoc))//'" not in module variables', &
                      ErrStat, ErrMsg, RoutineName)
      return
   else if (DstMesh%ID == 0) then
      call SetErrStat(ErrID_Fatal, 'DstMesh "'//trim(FAST_InputMeshName(DstMod, DstMeshLoc))//'" not in module variables', &
                      ErrStat, ErrMsg, RoutineName)
      return
   else if (DstDispMesh%ID == 0) then
      call SetErrStat(ErrID_Fatal, 'DstDispMesh "'//trim(FAST_OutputMeshName(DstMod, DstDispMeshLoc))//'" not in module variables', &
                      ErrStat, ErrMsg, RoutineName)
      return
   end if

   call FAST_InputMeshPointer(DstMod, Turbine, DstMeshLoc, .false., DstMesh, ErrStat2, ErrMsg2); if (Failed()) return
   call FAST_OutputMeshPointer(DstMod, Turbine, DstDispMeshLoc, DstDispMesh, ErrStat2, ErrMsg2); if (Failed()) return

   ! Create mapping description
   Mapping%Desc = trim(FAST_OutputMeshName(SrcMod, SrcMeshLoc))//" -> "// &
                  trim(FAST_InputMeshName(DstMod, DstMeshLoc))// &
                  " ["//trim(FAST_InputMeshName(SrcMod, SrcDispMeshLoc))// &
                  " -> "//trim(FAST_OutputMeshName(DstMod, DstDispMeshLoc))//"]"

   ! Initialize mapping structure
   Mapping%MapType = Map_LoadMesh
   Mapping%SrcModIdx = SrcMod%Idx
   Mapping%SrcModID = SrcMod%ID
   Mapping%SrcIns = SrcMod%Ins
   Mapping%DstModIdx = DstMod%Idx
   Mapping%DstModID = DstMod%ID
   Mapping%DstIns = DstMod%Ins
   Mapping%SrcMeshLoc = SrcMeshLoc
   Mapping%SrcDispMeshLoc = SrcDispMeshLoc
   Mapping%DstMeshLoc = DstMeshLoc
   Mapping%DstDispMeshLoc = DstDispMeshLoc

   ! Create mesh mapping
   call MeshMapCreate(SrcMesh, DstMesh, Mapping%MeshMap, ErrStat2, ErrMsg2); if (Failed()) return

   ! Create a copy of destination mesh in mapping for load summation
   call MeshCopy(DstMesh, Mapping%TmpLoadMesh, MESH_NEWCOPY, ErrStat2, ErrMsg2); if (Failed()) return

   ! Get mapping indices for linearized mesh mapping
   call InitMeshVarLocs(Mapping, SrcMod, DstMod, SrcMesh, DstMesh, SrcDispMesh, DstDispMesh)

   ! If the destination displacement mesh is not a sibling of the load mesh
   Mapping%DstUsesSibling = IsSiblingMesh(DstMesh, DstDispMesh)
   if (.not. Mapping%DstUsesSibling) then

      ! Print warning
      call WrScr('Warning: load mesh transfer "'//trim(Mapping%Desc)//'" does not use sibling mesh')

      ! Create temporary motion mesh as cousin of load mesh, this will be used for an intermediate transfer
      ! of the destination motion to the destination load locations
      call MeshCopy(SrcMesh=DstMesh, &
                    DestMesh=Mapping%TmpMotionMesh, &
                    CtrlCode=MESH_COUSIN, &
                    IOS=COMPONENT_OUTPUT, &
                    TranslationDisp=.true., &
                    Orientation=.true., &
                    RotationVel=.true., &
                    TranslationVel=.true., &
                    RotationAcc=.true., &
                    TranslationAcc=.true., &
                    ErrStat=ErrStat2, &
                    ErrMess=ErrMsg2)
      if (Failed()) return

      ! Determine transfer/linearization type for this auxiliary transfer
      Mapping%XfrTypeAux = MeshTransferType(DstDispMesh, Mapping%TmpMotionMesh)

      ! Create motion mapping from destination displacement to temporary motion mesh
      call MeshMapCreate(DstDispMesh, Mapping%TmpMotionMesh, Mapping%MeshMapAux, ErrStat2, ErrMsg2); if (Failed()) return

   end if

   ! Add mapping to array of mappings
   Mappings = [Mappings, Mapping]

contains
   logical function Failed()
      Failed = ErrStat2 >= AbortErrLev
      if (Failed) call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   end function

   ! IsSiblingMesh returns true if MeshB is a sibling of MeshA
   ! (can't just check pointers as they won't match after restart,
   ! also there can only be one sibling mesh so doesn't work for cousins)
   logical function IsSiblingMesh(MeshA, MeshB)
      type(MeshType), intent(in) :: MeshA, MeshB
      integer(IntKi)             :: i, j
      IsSiblingMesh = (MeshA%Nnodes == MeshB%Nnodes)
      if (.not. IsSiblingMesh) return
      IsSiblingMesh = IsSiblingMesh .and. &
                      all(MeshA%Position == MeshB%Position) .and. &
                      all(MeshA%RefOrientation == MeshB%RefOrientation)
      do i = 1, NELEMKINDS
         IsSiblingMesh = IsSiblingMesh .and. &
                         (MeshA%ElemTable(i)%nelem == MeshB%ElemTable(i)%nelem) .and. &
                         (MeshA%ElemTable(i)%XElement == MeshB%ElemTable(i)%XElement)
         if (.not. IsSiblingMesh) return
         do j = 1, MeshA%ElemTable(i)%nelem
            IsSiblingMesh = IsSiblingMesh .and. all(MeshA%ElemTable(i)%Elements(j)%ElemNodes == &
                                                    MeshB%ElemTable(i)%Elements(j)%ElemNodes)
         end do
      end do
   end function
end subroutine

subroutine MapMotionMesh(Turbine, Mappings, SrcMod, SrcMeshLoc, DstMod, DstMeshLoc, ErrStat, ErrMsg, Active)
   type(FAST_TurbineType), target         :: Turbine
   type(TC_MappingType), allocatable      :: Mappings(:)
   type(ModDataType), intent(in)          :: SrcMod, DstMod
   type(MeshLocType), intent(in)          :: SrcMeshLoc, DstMeshLoc
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg
   logical, optional, intent(in)          :: Active

   character(*), parameter                :: RoutineName = 'MapMotionMesh'
   integer(IntKi)                         :: ErrStat2
   character(ErrMsgLen)                   :: ErrMsg2
   type(TC_MappingType)                   :: Mapping
   type(MeshType), pointer                :: SrcMesh, DstMesh

   ErrStat = ErrID_None
   ErrMsg = ''

   ! If active argument is set to false, return
   if (present(Active)) then
      if (.not. Active) return
   end if

   ! Get mesh pointers
   call FAST_OutputMeshPointer(SrcMod, Turbine, SrcMeshLoc, SrcMesh, ErrStat2, ErrMsg2); if (Failed()) return
   call FAST_InputMeshPointer(DstMod, Turbine, DstMeshLoc, .false., DstMesh, ErrStat2, ErrMsg2); if (Failed()) return

   ! If source or destination meshes aren't commited, return
   if (.not. (SrcMesh%committed .and. DstMesh%committed)) return

   ! Check that all meshes in mapping have nonzero identifiers
   if (SrcMesh%ID == 0) then
      call SetErrStat(ErrID_Fatal, 'SrcMesh "'//trim(FAST_OutputMeshName(SrcMod, SrcMeshLoc))//'" not in module variables', &
                      ErrStat, ErrMsg, RoutineName)
      return
   else if (DstMesh%ID == 0) then
      call SetErrStat(ErrID_Fatal, 'DstMesh "'//trim(FAST_InputMeshName(DstMod, DstMeshLoc))//'" not in module variables', &
                      ErrStat, ErrMsg, RoutineName)
      return
   end if

   ! Create mapping description
   Mapping%Desc = trim(FAST_OutputMeshName(SrcMod, SrcMeshLoc))//" -> "// &
                  trim(FAST_InputMeshName(DstMod, DstMeshLoc))

   ! Initialize mapping structure
   Mapping%MapType = Map_MotionMesh
   Mapping%SrcModIdx = SrcMod%Idx
   Mapping%SrcModID = SrcMod%ID
   Mapping%SrcIns = SrcMod%Ins
   Mapping%DstModIdx = DstMod%Idx
   Mapping%DstModID = DstMod%ID
   Mapping%DstIns = DstMod%Ins
   Mapping%SrcMeshLoc = SrcMeshLoc
   Mapping%DstMeshLoc = DstMeshLoc

   ! Create mesh mapping
   call MeshMapCreate(SrcMesh, DstMesh, Mapping%MeshMap, ErrStat2, ErrMsg2); if (Failed()) return

   ! Get mapping indices for linearized mesh mapping
   call InitMeshVarLocs(Mapping, SrcMod, DstMod, SrcMesh, DstMesh)

   ! Add mapping to array of mappings
   Mappings = [Mappings, Mapping]

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine MapVariable(Maps, Key, SrcMod, DstMod, iVarSrc, iVarDst, ErrStat, ErrMsg, Active)
   type(TC_MappingType), allocatable      :: Maps(:)
   character(*), intent(in)               :: Key
   type(ModDataType), intent(in)          :: SrcMod, DstMod
   integer(IntKi), intent(in)             :: iVarSrc, iVarDst
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg
   logical, optional, intent(in)          :: Active
   type(TC_MappingType)                   :: Mapping

   ErrStat = ErrID_None
   ErrMsg = ''

   if (present(Active)) then
      if (.not. Active) return
   end if

   ! If either variable index is zero, return error
   if (iVarSrc == 0) then
      ErrStat = ErrID_Fatal
      ErrMsg = "Source variable in mapping '"//Key//"' is not active"
      return
   else if (iVarDst == 0) then
      ErrStat = ErrID_Fatal
      ErrMsg = "Destination variable in mapping '"//Key//"' is not active"
      return
   end if

   ! Verify that variables have compatible sizes
   ! If source variable has size 1, it can be mapped to multiple destination variables
   if ((SrcMod%Vars%y(iVarSrc)%Num > 1) .and. &
       (SrcMod%Vars%y(iVarSrc)%Num /= DstMod%Vars%u(iVarDst)%Num)) then
      ErrStat = ErrID_Fatal
      ErrMsg = "Variables in mapping '"//Key//"' have different sizes"
      return
   end if

   ! Initialize mapping structure
   Mapping%Desc = Key
   Mapping%MapType = Map_Variable
   Mapping%SrcModIdx = SrcMod%Idx
   Mapping%DstModIdx = DstMod%Idx
   Mapping%SrcModID = SrcMod%ID
   Mapping%DstModID = DstMod%ID
   Mapping%SrcIns = SrcMod%Ins
   Mapping%DstIns = DstMod%Ins
   Mapping%iVarSrc = iVarSrc
   Mapping%iVarDst = iVarDst

   Maps = [Maps, Mapping]
end subroutine

!> MapCustom creates a custom mapping that is not included in linearization.
!! Each custom mapping needs an entry in FAST_InputSolve to actually perform the transfer.
subroutine MapCustom(Maps, Desc, SrcMod, DstMod, Active)
   type(TC_MappingType), allocatable      :: Maps(:)
   character(*), intent(in)               :: Desc
   type(ModDataType), intent(in)          :: SrcMod, DstMod
   logical, optional, intent(in)          :: Active
   type(TC_MappingType)                   :: Mapping

   if (present(Active)) then
      if (.not. Active) return
   end if

   ! Initialize mapping structure
   Mapping%Desc = Desc
   Mapping%MapType = Map_Custom
   Mapping%SrcModIdx = SrcMod%Idx
   Mapping%DstModIdx = DstMod%Idx
   Mapping%SrcModID = SrcMod%ID
   Mapping%DstModID = DstMod%ID
   Mapping%SrcIns = SrcMod%Ins
   Mapping%DstIns = DstMod%Ins

   Maps = [Maps, Mapping]
end subroutine

subroutine InitMeshVarLocs(Mapping, SrcMod, DstMod, SrcMesh, DstMesh, SrcDispMesh, DstDispMesh)
   type(TC_MappingType), intent(inout)    :: Mapping
   type(ModDataType), intent(in)          :: SrcMod, DstMod
   type(MeshType), intent(in)             :: SrcMesh, DstMesh
   type(MeshType), optional, intent(in)   :: SrcDispMesh, DstDispMesh

   ! Save source and destination mesh ID
   Mapping%SrcMeshID = SrcMesh%ID
   Mapping%DstMeshID = DstMesh%ID

   ! Determine transfer type
   Mapping%XfrType = MeshTransferType(SrcMesh, DstMesh)

   ! Get data locations for variables of source mesh fields
   call FindVarByMeshAndField(SrcMod%Vars%y, SrcMesh%ID, VF_TransDisp, SrcMod%iyg, Mapping%iLocSrcTransDisp)
   call FindVarByMeshAndField(SrcMod%Vars%y, SrcMesh%ID, VF_TransVel, SrcMod%iyg, Mapping%iLocSrcTransVel)
   call FindVarByMeshAndField(SrcMod%Vars%y, SrcMesh%ID, VF_TransAcc, SrcMod%iyg, Mapping%iLocSrcTransAcc)
   call FindVarByMeshAndField(SrcMod%Vars%y, SrcMesh%ID, VF_Orientation, SrcMod%iyg, Mapping%iLocSrcOrientation)
   call FindVarByMeshAndField(SrcMod%Vars%y, SrcMesh%ID, VF_AngularVel, SrcMod%iyg, Mapping%iLocSrcAngularVel)
   call FindVarByMeshAndField(SrcMod%Vars%y, SrcMesh%ID, VF_AngularAcc, SrcMod%iyg, Mapping%iLocSrcAngularAcc)
   call FindVarByMeshAndField(SrcMod%Vars%y, SrcMesh%ID, VF_Force, SrcMod%iyg, Mapping%iLocSrcForce)
   call FindVarByMeshAndField(SrcMod%Vars%y, SrcMesh%ID, VF_Moment, SrcMod%iyg, Mapping%iLocSrcMoment)

   ! Get data locations for variables of destination mesh fields
   call FindVarByMeshAndField(DstMod%Vars%u, DstMesh%ID, VF_TransDisp, DstMod%iug, Mapping%iLocDstTransDisp)
   call FindVarByMeshAndField(DstMod%Vars%u, DstMesh%ID, VF_TransVel, DstMod%iug, Mapping%iLocDstTransVel)
   call FindVarByMeshAndField(DstMod%Vars%u, DstMesh%ID, VF_TransAcc, DstMod%iug, Mapping%iLocDstTransAcc)
   call FindVarByMeshAndField(DstMod%Vars%u, DstMesh%ID, VF_Orientation, DstMod%iug, Mapping%iLocDstOrientation)
   call FindVarByMeshAndField(DstMod%Vars%u, DstMesh%ID, VF_AngularVel, DstMod%iug, Mapping%iLocDstAngularVel)
   call FindVarByMeshAndField(DstMod%Vars%u, DstMesh%ID, VF_AngularAcc, DstMod%iug, Mapping%iLocDstAngularAcc)
   call FindVarByMeshAndField(DstMod%Vars%u, DstMesh%ID, VF_Force, DstMod%iug, Mapping%iLocDstForce)
   call FindVarByMeshAndField(DstMod%Vars%u, DstMesh%ID, VF_Moment, DstMod%iug, Mapping%iLocDstMoment)

   if (present(SrcDispMesh)) then
      Mapping%SrcDispMeshID = SrcDispMesh%ID
      call FindVarByMeshAndField(SrcMod%Vars%u, SrcDispMesh%ID, VF_TransDisp, SrcMod%iug, Mapping%iLocSrcDispTransDisp)
   end if

   if (present(DstDispMesh)) then
      Mapping%DstDispMeshID = DstDispMesh%ID
      call FindVarByMeshAndField(DstMod%Vars%y, DstDispMesh%ID, VF_TransDisp, DstMod%iyg, Mapping%iLocDstDispTransDisp)
      call FindVarByMeshAndField(DstMod%Vars%y, DstDispMesh%ID, VF_Orientation, DstMod%iyg, Mapping%iLocDstDispOrientation)
   end if

contains
   subroutine FindVarByMeshAndField(VarAry, MeshID, Field, iGbl, iLoc)
      type(ModVarType), intent(in)  :: VarAry(:)
      integer(IntKi), intent(in)    :: MeshID, Field, iGbl
      integer(IntKi), intent(out)   :: iLoc(2)
      integer(IntKi)                :: i

      ! Initialize locations
      iLoc = 0

      ! Loop through variables, if variable's mesh ID and field matches given values, return
      do i = 1, size(VarAry)
         if ((VarAry(i)%MeshID == MeshID) .and. (VarAry(i)%Field == Field)) then
            iLoc = VarAry(i)%iLoc + iGbl - 1
            return
         end if
      end do
   end subroutine
end subroutine

function MeshTransferType(SrcMesh, DstMesh) result(XfrType)
   type(MeshType), intent(in) :: SrcMesh, DstMesh
   integer(IntKi)             :: XfrType
   if ((SrcMesh%ElemTable(ELEMENT_POINT)%nelem > 0) .and. (DstMesh%ElemTable(ELEMENT_POINT)%nelem > 0)) then
      XfrType = Xfr_Point_to_Point
   else if ((SrcMesh%ElemTable(ELEMENT_POINT)%nelem > 0) .and. (DstMesh%ElemTable(ELEMENT_LINE2)%nelem > 0)) then
      XfrType = Xfr_Point_to_Line2
   else if ((SrcMesh%ElemTable(ELEMENT_LINE2)%nelem > 0) .and. (DstMesh%ElemTable(ELEMENT_POINT)%nelem > 0)) then
      XfrType = Xfr_Line2_to_Point
   else if ((SrcMesh%ElemTable(ELEMENT_LINE2)%nelem > 0) .and. (DstMesh%ElemTable(ELEMENT_LINE2)%nelem > 0)) then
      XfrType = Xfr_Line2_to_Line2
   else
      XfrType = Xfr_Invalid
   end if
end function

subroutine FAST_LinearizeMappings(Turbine, Mods, Mappings, ModOrder, ErrStat, ErrMsg, dUdu, dUdy)
   type(FAST_TurbineType), target, intent(inout)   :: Turbine     !< Turbine type
   type(ModDataType), intent(in)                   :: Mods(:)     !< Module data
   type(TC_MappingType), intent(inout)             :: Mappings(:)
   integer(IntKi), intent(in)                      :: ModOrder(:)
   integer(IntKi), intent(out)                     :: ErrStat
   character(*), intent(out)                       :: ErrMsg
   real(R8Ki), intent(inout)                       :: dUdu(:, :), dUdy(:, :)

   character(*), parameter       :: RoutineName = 'FAST_LinearizeMappings'
   integer(IntKi)                :: ErrStat2
   character(ErrMsgLen)          :: ErrMsg2
   integer(IntKi)                :: i, j, k
   type(MeshType), pointer       :: SrcMesh, DstMesh
   type(MeshType), pointer       :: SrcDispMesh, DstDispMesh

   ErrStat = ErrID_None
   ErrMsg = ''

   ! Loop through modules in specified order
   do i = 1, size(ModOrder)

      ! Loop through mappings where this module is the destination
      do j = 1, size(Mods((ModOrder(i)))%DstMaps)
         associate (Mapping => Mappings(Mods((ModOrder(i)))%DstMaps(j)))

            ! Select based on type of mapping
            select case (Mapping%MapType)

            case (Map_Variable)

               associate (SrcMod => Mods(Mapping%SrcModIdx), &
                          DstMod => Mods(Mapping%DstModIdx), &
                          SrcVar => Mods(Mapping%SrcModIdx)%Vars%y(Mapping%iVarSrc), &
                          DstVar => Mods(Mapping%DstModIdx)%Vars%u(Mapping%iVarDst))
                  if (SrcVar%Num == 1) then
                     ! Map rank 0 source var to rank 1 destination var
                     do k = 0, DstVar%Num - 1
                        dUdy(DstMod%iug + DstVar%iLoc(1) + k - 1, SrcMod%iyg + SrcVar%iLoc(1) - 1) = -1.0_R8Ki
                     end do
                  else
                     ! Map rank 1 source var to rank 1 destination var
                     do k = 0, SrcVar%Num - 1
                        dUdy(DstMod%iug + DstVar%iLoc(1) + k - 1, SrcMod%iyg + SrcVar%iLoc(1) + k - 1) = -1.0_R8Ki
                     end do
                  end if
               end associate

            case (Map_MotionMesh)

               ! Get source and destination meshes
               call FAST_OutputMeshPointer(Mods(Mapping%SrcModIdx), Turbine, Mapping%SrcMeshLoc, SrcMesh, ErrStat2, ErrMsg2); if (Failed()) return
               call FAST_InputMeshPointer(Mods(Mapping%DstModIdx), Turbine, Mapping%DstMeshLoc, .false., DstMesh, ErrStat2, ErrMsg2); if (Failed()) return

               ! Perform linearization based on transfer type
               call LinearizeMeshTransfer(Mapping%XfrType, SrcMesh, DstMesh, Mapping%MeshMap); if (Failed()) return

               ! Copy linearization matrices to global dUdy matrix
               call Assemble_dUdy_Motions(Mapping, dUdy)

               ! Copy linearization matrices to global dUdu matrix
               call Assemble_dUdu(Mapping, dUdu)

            case (Map_LoadMesh)

               ! Get source and destination meshes
               call FAST_OutputMeshPointer(Mods(Mapping%SrcModIdx), Turbine, Mapping%SrcMeshLoc, SrcMesh, ErrStat2, ErrMsg2); if (Failed()) return
               call FAST_InputMeshPointer(Mods(Mapping%DstModIdx), Turbine, Mapping%DstMeshLoc, .false., DstMesh, ErrStat2, ErrMsg2); if (Failed()) return

               ! Get source and destination displacement meshes
               call FAST_InputMeshPointer(Mods(Mapping%SrcModIdx), Turbine, Mapping%SrcDispMeshLoc, .false., SrcDispMesh, ErrStat2, ErrMsg2); if (Failed()) return
               call FAST_OutputMeshPointer(Mods(Mapping%DstModIdx), Turbine, Mapping%DstDispMeshLoc, DstDispMesh, ErrStat2, ErrMsg2); if (Failed()) return

               ! If DstDispMesh is a sibling of DstMesh
               if (Mapping%DstUsesSibling) then

                  ! Linearize the load mesh transfer
                  call LinearizeMeshTransfer(Mapping%XfrType, SrcMesh, DstMesh, Mapping%MeshMap, SrcDispMesh, DstDispMesh); if (Failed()) return

               else

                  ! Transfer destination displacement mesh to temporary motion mesh (cousin of destionation load mesh)
                  call TransferMesh(Mapping%XfrTypeAux, DstDispMesh, Mapping%TmpMotionMesh, Mapping%MeshMapAux); if (Failed()) return

                  ! Linearize the motion mesh transfer
                  call LinearizeMeshTransfer(Mapping%XfrTypeAux, DstDispMesh, Mapping%TmpMotionMesh, Mapping%MeshMapAux); if (Failed()) return

                  ! Linearize the load mesh transfer
                  call LinearizeMeshTransfer(Mapping%XfrType, SrcMesh, DstMesh, Mapping%MeshMap, SrcDispMesh, Mapping%TmpMotionMesh); if (Failed()) return
               end if

               ! Copy linearization matrices to global dUdy matrix
               call Assemble_dUdy_Loads(Mapping, dUdy)

               ! Copy linearization matrices to global dUdu matrix
               call Assemble_dUdu(Mapping, dUdu)

            end select

         end associate
      end do
   end do

contains

   ! LinearizeMeshTransfer calls the specific linearization function based on
   ! transfer type (Point_to_Point, Point_to_Line2, etc.)
   subroutine LinearizeMeshTransfer(Typ, Src, Dst, MMap, SrcDisp, DstDisp)
      integer(IntKi), intent(in)             :: Typ
      type(MeshType), intent(in)             :: Src, Dst
      type(MeshMapType), intent(inout)       :: MMap
      type(MeshType), optional, intent(in)   :: SrcDisp, DstDisp
      select case (Typ)
      case (Xfr_Point_to_Point)
         call Linearize_Point_to_Point(Src, Dst, MMap, ErrStat2, ErrMsg2, SrcDisp, DstDisp)
      case (Xfr_Point_to_Line2)
         call Linearize_Point_to_Line2(Src, Dst, MMap, ErrStat2, ErrMsg2, SrcDisp, DstDisp)
      case (Xfr_Line2_to_Point)
         call Linearize_Line2_to_Point(Src, Dst, MMap, ErrStat2, ErrMsg2, SrcDisp, DstDisp)
      case (Xfr_Line2_to_Line2)
         call Linearize_Line2_to_Line2(Src, Dst, MMap, ErrStat2, ErrMsg2, SrcDisp, DstDisp)
      case default
         ErrStat2 = ErrID_Fatal
         ErrMsg2 = "LinearizeMeshTransfer: unknown transfer type: "//Num2LStr(Typ)
      end select
   end subroutine

   ! MeshTransfer calls the specific transfer function based on
   ! transfer type (Point_to_Point, Point_to_Line2, etc.)
   subroutine TransferMesh(Typ, Src, Dst, MMap, SrcDisp, DstDisp)
      integer(IntKi), intent(in)             :: Typ
      type(MeshType), intent(in)             :: Src
      type(MeshType), intent(inout)          :: Dst
      type(MeshMapType), intent(inout)       :: MMap
      type(MeshType), optional, intent(in)   :: SrcDisp, DstDisp
      select case (Typ)
      case (Xfr_Point_to_Point)
         call Transfer_Point_to_Point(Src, Dst, MMap, ErrStat2, ErrMsg2, SrcDisp, DstDisp)
      case (Xfr_Point_to_Line2)
         call Transfer_Point_to_Line2(Src, Dst, MMap, ErrStat2, ErrMsg2, SrcDisp, DstDisp)
      case (Xfr_Line2_to_Point)
         call Transfer_Line2_to_Point(Src, Dst, MMap, ErrStat2, ErrMsg2, SrcDisp, DstDisp)
      case (Xfr_Line2_to_Line2)
         call Transfer_Line2_to_Line2(Src, Dst, MMap, ErrStat2, ErrMsg2, SrcDisp, DstDisp)
      case default
         ErrStat2 = ErrID_Fatal
         ErrMsg2 = "TransferMeshTransfer: unknown transfer type: "//Num2LStr(Typ)
      end select
   end subroutine

   logical function Failed()
      Failed = ErrStat2 >= AbortErrLev
      if (Failed) call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   end function
end subroutine

subroutine Assemble_dUdu(Mapping, dUdu)
   type(TC_MappingType), intent(in) :: Mapping
   real(R8Ki), intent(inout)        :: dUdu(:, :)

   ! Effect of input Translation Displacement on input Translation Velocity
   if (allocated(Mapping%MeshMap%dM%tv_uD)) then
      call SumBlock(Mapping%iLocDstTransVel, Mapping%iLocDstTransDisp, Mapping%MeshMap%dM%tv_uD, dUdu)
   end if

   ! Effect of input Translation Displacement on input Translation Acceleration
   if (allocated(Mapping%MeshMap%dM%ta_uD)) then
      call SumBlock(Mapping%iLocDstTransAcc, Mapping%iLocDstTransDisp, Mapping%MeshMap%dM%ta_uD, dUdu)
   end if

   ! Effect of input Translation Displacement on input Moments
   if (allocated(Mapping%MeshMap%dM%M_uS)) then
      call SumBlock(Mapping%iLocDstMoment, Mapping%iLocSrcDispTransDisp, Mapping%MeshMap%dM%M_uS, dUdu)
   end if
end subroutine

!> Assemble_dUdy_Loads assembles the linearization matrices for transfer of
!! load fields between two meshes. It sets the following block matrix, which
!! is the dUdy block for transfering output (source) mesh  to the input
!! (destination) mesh :
!! M = -| M_li   0    | * M_mi | F^S |
!!      | M_fm   M_li |        | M^S |
subroutine Assemble_dUdy_Loads(Mapping, dUdy)
   type(TC_MappingType), intent(inout) :: Mapping
   real(R8Ki), intent(inout)           :: dUdy(:, :)

   ! Load identity
   if (allocated(Mapping%MeshMap%dM%li)) then
      call SumBlock(Mapping%iLocDstForce, Mapping%iLocSrcForce, Mapping%MeshMap%dM%li, dUdy)
      call SumBlock(Mapping%iLocDstMoment, Mapping%iLocSrcMoment, Mapping%MeshMap%dM%li, dUdy)
   end if

   ! Force to Moment
   if (allocated(Mapping%MeshMap%dM%m_f)) then
      call SumBlock(Mapping%iLocDstMoment, Mapping%iLocSrcForce, Mapping%MeshMap%dM%m_f, dUdy)
   end if

   ! Destination Translation Displacement to Moment
   if (allocated(Mapping%MeshMap%dM%m_uD)) then
      if (Mapping%DstUsesSibling) then
         ! Direct transfer
         call SumBlock(Mapping%iLocDstMoment, Mapping%iLocDstDispTransDisp, Mapping%MeshMap%dM%m_uD, dUdy)
      else
         ! call SumBlock(Mapping%iLocDstMoment, [Mapping%iLocDstDispTransDisp(1), Mapping%iLocDstDispTransDisp(1)  + size(Mapping%MeshMap%dM%m_uD,2) - 1], Mapping%MeshMap%dM%m_uD, dUdy)
         ! Compose linearization of motion and loads
         Mapping%TmpMatrix = matmul(Mapping%MeshMap%dM%m_uD, Mapping%MeshMapAux%dM%mi)
         call SumBlock(Mapping%iLocDstMoment, Mapping%iLocDstDispTransDisp, Mapping%TmpMatrix, dUdy)
         Mapping%TmpMatrix = matmul(Mapping%MeshMap%dM%m_uD, Mapping%MeshMapAux%dM%fx_p)
         call SumBlock(Mapping%iLocDstMoment, Mapping%iLocDstDispOrientation, Mapping%TmpMatrix, dUdy)
      end if
   end if
end subroutine

!> Assemble_dUdy_Motions assembles the linearization matrices for transfer of
!! motion fields between two meshes. It set the following block matrix, which
!! is the dUdy block for transfering output (source) mesh  to the input
!! (destination) mesh :
!! M = -| M_mi       M_f_p    0      0           0      0     |
!!      | 0          M_mi     0      0           0      0     |
!!      | M_tv_uS    0        M_mi   M_f_p       0      0     |
!!      | 0          0        0      M_mi        0      0     |
!!      | M_ta_uS    0        0      M_ta_rv     M_mi   M_f_p |
!!      | 0          0        0      0           0      M_mi  |
!! where the matrices correspond to
!! u^S, theta^S, v^S, omega^S, a^S, alpha^S
subroutine Assemble_dUdy_Motions(Mapping, dUdy)
   type(TC_MappingType), intent(inout) :: Mapping
   real(R8Ki), intent(inout)           :: dUdy(:, :)

   ! Motion identity
   if (allocated(Mapping%MeshMap%dM%mi)) then
      call SumBlock(Mapping%iLocDstTransDisp, Mapping%iLocSrcTransDisp, Mapping%MeshMap%dM%mi, dUdy)
      call SumBlock(Mapping%iLocDstOrientation, Mapping%iLocSrcOrientation, Mapping%MeshMap%dM%mi, dUdy)
      call SumBlock(Mapping%iLocDstTransVel, Mapping%iLocSrcTransVel, Mapping%MeshMap%dM%mi, dUdy)
      call SumBlock(Mapping%iLocDstAngularVel, Mapping%iLocSrcAngularVel, Mapping%MeshMap%dM%mi, dUdy)
      call SumBlock(Mapping%iLocDstTransAcc, Mapping%iLocSrcTransAcc, Mapping%MeshMap%dM%mi, dUdy)
      call SumBlock(Mapping%iLocDstAngularAcc, Mapping%iLocSrcAngularAcc, Mapping%MeshMap%dM%mi, dUdy)
   end if

   ! Rotation to Translation
   if (allocated(Mapping%MeshMap%dM%fx_p)) then
      call SumBlock(Mapping%iLocDstTransDisp, Mapping%iLocSrcOrientation, Mapping%MeshMap%dM%fx_p, dUdy)
      call SumBlock(Mapping%iLocDstTransVel, Mapping%iLocSrcAngularVel, Mapping%MeshMap%dM%fx_p, dUdy)
      call SumBlock(Mapping%iLocDstTransAcc, Mapping%iLocSrcAngularAcc, Mapping%MeshMap%dM%fx_p, dUdy)
   end if

   ! Translation displacement to Translation velocity
   if (allocated(Mapping%MeshMap%dM%tv_us)) then
      call SumBlock(Mapping%iLocDstTransVel, Mapping%iLocSrcTransDisp, Mapping%MeshMap%dM%tv_us, dUdy)
   end if

   ! Translation displacement to Translation acceleration
   if (allocated(Mapping%MeshMap%dM%ta_us)) then
      call SumBlock(Mapping%iLocDstTransAcc, Mapping%iLocSrcTransDisp, Mapping%MeshMap%dM%ta_us, dUdy)
   end if

   ! Angular velocity to Translation acceleration
   if (allocated(Mapping%MeshMap%dM%ta_rv)) then
      call SumBlock(Mapping%iLocDstTransAcc, Mapping%iLocSrcAngularVel, Mapping%MeshMap%dM%ta_rv, dUdy)
   end if
end subroutine

subroutine SumBlock(iLocRow, iLocCol, SrcM, DstM)
   integer(IntKi), intent(in)    :: iLocRow(2), iLocCol(2)
   real(R8Ki), intent(in)        :: SrcM(:, :)
   real(R8Ki), intent(inout)     :: DstM(:, :)
   if (iLocRow(1) > 0 .and. iLocCol(1) > 0) then
      ! Subtracts the source matrix from the destination sub-matrix
      associate (DstSubM => DstM(iLocRow(1):iLocRow(2), iLocCol(1):iLocCol(2)))
         DstSubM = DstSubM - SrcM
      end associate
   end if
end subroutine

subroutine FAST_InputSolve(Turbine, Mods, Mappings, iMod, ErrStat, ErrMsg, UseU)
   type(FAST_TurbineType), intent(inout)  :: Turbine     !< Turbine type
   type(ModDataType), intent(in)          :: Mods(:)     !< Module data
   type(TC_MappingType), intent(inout)    :: Mappings(:)
   integer(IntKi), intent(in)             :: iMod        !< Index of module in Mods to do input solve
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg
   logical, intent(in)                    :: UseU        ! Flag to transfer to u instead of Input

   character(*), parameter       :: RoutineName = 'FAST_InputSolve'
   integer(IntKi)                :: ErrStat2
   character(ErrMsgLen)          :: ErrMsg2
   type(MeshType), pointer       :: SrcMesh, DstMesh
   type(MeshType), pointer       :: SrcDispMesh, DstDispMesh
   integer(IntKi)                :: i, j, k

   ErrStat = ErrID_None
   ErrMsg = ''

   ! Loop through mappings where this module is the destination
   do i = 1, size(Mods(iMod)%DstMaps)
      associate (Mapping => Mappings(Mods(iMod)%DstMaps(i)))

         ! Select based on type of mapping
         select case (Mapping%MapType)

         case (Map_Custom)

            call Custom_InputSolve(Turbine, Mapping, ErrStat2, ErrMsg2, UseU)
            if (Failed()) return

         case (Map_Variable)

         case (Map_MotionMesh)

            ! Get source and destination meshes
            call FAST_OutputMeshPointer(Mods(Mapping%SrcModIdx), Turbine, Mapping%SrcMeshLoc, SrcMesh, ErrStat2, ErrMsg2); if (Failed()) return
            call FAST_InputMeshPointer(Mods(Mapping%DstModIdx), Turbine, Mapping%DstMeshLoc, UseU, DstMesh, ErrStat2, ErrMsg2); if (Failed()) return

            ! Perform linearization based on transfer type
            select case (Mapping%XfrType)
            case (Xfr_Point_to_Point)
               call Linearize_Point_to_Point(SrcMesh, DstMesh, Mapping%MeshMap, ErrStat2, ErrMsg2)
            case (Xfr_Point_to_Line2)
               call Linearize_Point_to_Line2(SrcMesh, DstMesh, Mapping%MeshMap, ErrStat2, ErrMsg2)
            case (Xfr_Line2_to_Point)
               call Linearize_Line2_to_Point(SrcMesh, DstMesh, Mapping%MeshMap, ErrStat2, ErrMsg2)
            case (Xfr_Line2_to_Line2)
               call Linearize_Line2_to_Line2(SrcMesh, DstMesh, Mapping%MeshMap, ErrStat2, ErrMsg2)
            end select
            if (Failed()) return

         case (Map_LoadMesh)

            ! Get source and destination meshes
            call FAST_OutputMeshPointer(Mods(Mapping%SrcModIdx), Turbine, Mapping%SrcMeshLoc, SrcMesh, ErrStat2, ErrMsg2); if (Failed()) return
            call FAST_InputMeshPointer(Mods(Mapping%DstModIdx), Turbine, Mapping%DstMeshLoc, UseU, DstMesh, ErrStat2, ErrMsg2); if (Failed()) return

            ! Get source and destination displacement meshes
            call FAST_InputMeshPointer(Mods(Mapping%SrcModIdx), Turbine, Mapping%SrcDispMeshLoc, UseU, SrcDispMesh, ErrStat2, ErrMsg2); if (Failed()) return
            call FAST_OutputMeshPointer(Mods(Mapping%DstModIdx), Turbine, Mapping%DstDispMeshLoc, DstDispMesh, ErrStat2, ErrMsg2); if (Failed()) return

            ! Perform linearization based on transfer type
            select case (Mapping%XfrType)
            case (Xfr_Point_to_Point)
               call Linearize_Point_to_Point(SrcMesh, DstMesh, Mapping%MeshMap, ErrStat2, ErrMsg2, SrcDispMesh, DstDispMesh)
            case (Xfr_Point_to_Line2)
               call Linearize_Point_to_Line2(SrcMesh, DstMesh, Mapping%MeshMap, ErrStat2, ErrMsg2, SrcDispMesh, DstDispMesh)
            case (Xfr_Line2_to_Point)
               call Linearize_Line2_to_Point(SrcMesh, DstMesh, Mapping%MeshMap, ErrStat2, ErrMsg2, SrcDispMesh, DstDispMesh)
            case (Xfr_Line2_to_Line2)
               call Linearize_Line2_to_Line2(SrcMesh, DstMesh, Mapping%MeshMap, ErrStat2, ErrMsg2, SrcDispMesh, DstDispMesh)
            end select
            if (Failed()) return

         end select

      end associate
   end do

contains
   logical function Failed()
      Failed = ErrStat2 /= ErrID_None
      if (Failed) call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, &
                                  RoutineName//':Module='//trim(Mods(iMod)%Abbr)//', Instance='//Num2LStr(Mods(iMod)%Ins))
   end function
end subroutine

subroutine Custom_InputSolve(T, Mapping, ErrStat, ErrMsg, UseU)
   type(FAST_TurbineType), target, intent(inout)   :: T     !< Turbine type
   type(TC_MappingType), intent(in)                :: Mapping
   integer(IntKi), intent(out)                     :: ErrStat
   character(*), intent(out)                       :: ErrMsg
   logical, intent(in)                             :: UseU        ! Flag to transfer to u instead of Input

   character(*), parameter                :: RoutineName = 'Custom_InputSolve'
   integer(IntKi)                         :: ErrStat2
   character(ErrMsgLen)                   :: ErrMsg2
   integer(IntKi)                         :: i, j, k
   real(ReKi)                             :: z, u, v, mean_vel
   type(AD_InputType), pointer            :: u_AD
   type(ED_InputType), pointer            :: u_ED
   type(ExtLd_InputType), pointer         :: u_ExtLd
   type(InflowWind_InputType), pointer    :: u_IfW
   type(MD_InputType), pointer            :: u_MD
   type(SD_InputType), pointer            :: u_SD
   type(SrvD_InputType), pointer          :: u_SrvD

   ErrStat = ErrID_None
   ErrMsg = ''

   select case (Mapping%DstModID)

   case (Module_AD)
      if (UseU) then
         u_AD => T%AD%u
      else
         u_AD => T%AD%Input(1)
      end if
   case (Module_ED)
      if (UseU) then
         u_ED => T%ED%u
      else
         u_ED => T%ED%Input(1)
      end if
   case (Module_ExtLd)
      u_ExtLd => T%ExtLd%u
   case (Module_IfW)
      if (UseU) then
         u_IfW => T%IfW%u
      else
         u_IfW => T%IfW%Input(1)
      end if
   case (Module_SD)
      if (UseU) then
         u_SD => T%SD%u
      else
         u_SD => T%SD%Input(1)
      end if
   case (Module_SrvD)
      if (UseU) then
         u_SrvD => T%SrvD%u
      else
         u_SrvD => T%SrvD%Input(1)
      end if
   end select

   ! Select based on mapping description
   select case (Mapping%Desc)

!-------------------------------------------------------------------------------
! AeroDyn Inputs
!-------------------------------------------------------------------------------

   case (Custom_SrvD_to_AD)

      ! Set Conrol parameter (i.e. flaps) if using ServoDyn bem:
      ! This takes in flap deflection for each blade (only one flap deflection angle per blade),
      ! from ServoDyn (which comes from Bladed style DLL controller)
      ! Commanded Airfoil UserProp for blade (must be same units as given in AD15 airfoil tables)
      ! This is passed to AD15 to be interpolated with the airfoil table userprop column
      ! (might be used for airfoil flap angles for example)
      ! Must be same units as given in airfoil (no unit conversions handled in code)ß
      do i = 1, size(T%AD%u%rotors(1)%UserProp, dim=2)   ! Blade
         u_AD%rotors(1)%UserProp(:, i) = T%SrvD%y%BlAirfoilCom(i)
      end do

!-------------------------------------------------------------------------------
! ExtLoads Inputs
!-------------------------------------------------------------------------------

   case (Custom_ED_to_ExtLd)

      u_ExtLd%az = T%ED%y%LSSTipPxa
      u_ExtLd%DX_u%bldPitch(:) = T%ED%y%BlPitch

!-------------------------------------------------------------------------------
! InflowWind Inputs
!-------------------------------------------------------------------------------

   case (Custom_ED_to_IfW)

      ! This section should be refactored so that IfW uses a hub point mesh
      u_IfW%HubPosition = T%ED%y%HubPtMotion%Position(:, 1) + &
                          T%ED%y%HubPtMotion%TranslationDisp(:, 1)
      u_IfW%HubOrientation = T%ED%y%HubPtMotion%Orientation(:, :, 1)

      ! Set Lidar position directly from hub motion mesh
      u_IfW%lidar%HubDisplacementX = T%ED%y%HubPtMotion%TranslationDisp(1, 1)
      u_IfW%lidar%HubDisplacementY = T%ED%y%HubPtMotion%TranslationDisp(2, 1)
      u_IfW%lidar%HubDisplacementZ = T%ED%y%HubPtMotion%TranslationDisp(3, 1)

   case (Custom_SrvD_to_IfW)

      ! Set hub position so ServoDyn can get hub wind speed
      u_IfW%PositionXYZ(:,1) = T%ED%y%HubPtMotion%Position(:, 1)

!-------------------------------------------------------------------------------
! MoorDyn Inputs
!-------------------------------------------------------------------------------

   case (Custom_SrvD_to_MD)

      if (allocated(u_MD%DeltaL) .and. allocated(T%SrvD%y%CableDeltaL)) then
         u_MD%DeltaL = T%SrvD%y%CableDeltaL      ! these should be sized identically during init
      end if

      if (allocated(u_MD%DeltaLdot) .and. allocated(T%SrvD%y%CableDeltaLdot)) then
         u_MD%DeltaLdot = T%SrvD%y%CableDeltaLdot   ! these should be sized identically during init
      end if

!-------------------------------------------------------------------------------
! SubDyn Inputs
!-------------------------------------------------------------------------------

   case (Custom_SrvD_to_SD)

      if (allocated(u_SD%CableDeltaL) .and. allocated(T%SrvD%y%CableDeltaL)) then
         u_SD%CableDeltaL = T%SrvD%y%CableDeltaL   ! these should be sized identically during init
      end if

!-------------------------------------------------------------------------------
! ServoDyn Inputs
!-------------------------------------------------------------------------------

   case (Custom_BD_to_SrvD)

      u_SrvD%RootMxc(Mapping%SrcIns) = T%BD%y(Mapping%SrcIns)%RootMxr*cos(T%ED%y%BlPitch(Mapping%SrcIns)) + &
                                       T%BD%y(Mapping%SrcIns)%RootMyr*sin(T%ED%y%BlPitch(Mapping%SrcIns))
      u_SrvD%RootMyc(Mapping%SrcIns) = -T%BD%y(Mapping%SrcIns)%RootMxr*sin(T%ED%y%BlPitch(Mapping%SrcIns)) + &
                                       T%BD%y(Mapping%SrcIns)%RootMyr*cos(T%ED%y%BlPitch(Mapping%SrcIns))

   case (Custom_ED_to_SrvD)

      ! Blade root moment if not using BeamDyn
      if (T%p_FAST%CompElast /= Module_BD) then
         u_SrvD%RootMxc = T%ED%y%RootMxc ! fixed-size arrays: always size 3
         u_SrvD%RootMyc = T%ED%y%RootMyc ! fixed-size arrays: always size 3
      end if

      u_SrvD%YawAngle = T%ED%y%YawAngle ! nacelle yaw plus platform yaw
      u_SrvD%YawErr = u_SrvD%WindDir - u_SrvD%YawAngle ! the nacelle yaw error estimate (positive about zi-axis)

      u_SrvD%Yaw = T%ED%y%Yaw  ! nacelle yaw
      u_SrvD%YawRate = T%ED%y%YawRate
      u_SrvD%BlPitch = T%ED%y%BlPitch
      u_SrvD%LSS_Spd = T%ED%y%LSS_Spd
      u_SrvD%HSS_Spd = T%ED%y%HSS_Spd
      u_SrvD%RotSpeed = T%ED%y%RotSpeed

      u_SrvD%YawBrTAxp = T%ED%y%YawBrTAxp
      u_SrvD%YawBrTAyp = T%ED%y%YawBrTAyp
      u_SrvD%LSSTipPxa = T%ED%y%LSSTipPxa

      u_SrvD%LSSTipMxa = T%ED%y%LSSTipMxa
      u_SrvD%LSSTipMya = T%ED%y%LSSTipMya
      u_SrvD%LSSTipMza = T%ED%y%LSSTipMza
      u_SrvD%LSSTipMys = T%ED%y%LSSTipMys
      u_SrvD%LSSTipMzs = T%ED%y%LSSTipMzs

      u_SrvD%YawBrMyn = T%ED%y%YawBrMyn
      u_SrvD%YawBrMzn = T%ED%y%YawBrMzn
      u_SrvD%NcIMURAxs = T%ED%y%NcIMURAxs
      u_SrvD%NcIMURAys = T%ED%y%NcIMURAys
      u_SrvD%NcIMURAzs = T%ED%y%NcIMURAzs

      u_SrvD%RotPwr = T%ED%y%RotPwr

      u_SrvD%LSShftFxa = T%ED%y%LSShftFxa
      u_SrvD%LSShftFys = T%ED%y%LSShftFys
      u_SrvD%LSShftFzs = T%ED%y%LSShftFzs

   case (Custom_IfW_to_SrvD)

      u_SrvD%WindDir = atan2(T%IfW%y%VelocityUVW(2, 1), T%IfW%y%VelocityUVW(1, 1))
      u_SrvD%HorWindV = sqrt(T%IfW%y%VelocityUVW(1, 1)**2 + T%IfW%y%VelocityUVW(2, 1)**2)
      if (allocated(T%IfW%y%lidar%LidSpeed)) u_SrvD%LidSpeed = T%IfW%y%lidar%LidSpeed
      if (allocated(T%IfW%y%lidar%MsrPositionsX)) u_SrvD%MsrPositionsX = T%IfW%y%lidar%MsrPositionsX
      if (allocated(T%IfW%y%lidar%MsrPositionsY)) u_SrvD%MsrPositionsY = T%IfW%y%lidar%MsrPositionsY
      if (allocated(T%IfW%y%lidar%MsrPositionsZ)) u_SrvD%MsrPositionsZ = T%IfW%y%lidar%MsrPositionsZ
      u_SrvD%YawErr = u_SrvD%WindDir - u_SrvD%YawAngle ! the nacelle yaw error estimate (positive about zi-axis)

   case (Custom_ExtInfw_to_SrvD)

      u_SrvD%WindDir = ATAN2(T%ExtInfw%y%v(1), T%ExtInfw%y%u(1))
      u_SrvD%HorWindV = SQRT(T%ExtInfw%y%u(1)**2 + T%ExtInfw%y%v(1)**2)
      if (allocated(u_SrvD%LidSpeed)) u_SrvD%LidSpeed = 0.0
      if (allocated(u_SrvD%MsrPositionsX)) u_SrvD%MsrPositionsX = 0.0
      if (allocated(u_SrvD%MsrPositionsY)) u_SrvD%MsrPositionsY = 0.0
      if (allocated(u_SrvD%MsrPositionsz)) u_SrvD%MsrPositionsz = 0.0
      u_SrvD%YawErr = u_SrvD%WindDir - u_SrvD%YawAngle ! the nacelle yaw error estimate (positive about zi-axis)

   case (Custom_ExtLd_to_SrvD)

      pi = acos(-1.0)
      z = T%ED%y%HubPtMotion%Position(3, 1)
      mean_vel = T%ExtLd%p%vel_mean*((z/T%ExtLd%p%z_ref)**T%ExtLd%p%shear_exp)
      u = -mean_vel*sin(T%ExtLd%p%wind_dir*pi/180.0)
      v = -mean_vel*cos(T%ExtLd%p%wind_dir*pi/180.0)
      u_SrvD%HorWindV = mean_vel
      u_SrvD%WindDir = atan2(v, u)
      if (allocated(u_SrvD%LidSpeed)) u_SrvD%LidSpeed = 0.0
      if (allocated(u_SrvD%MsrPositionsX)) u_SrvD%MsrPositionsX = 0.0
      if (allocated(u_SrvD%MsrPositionsY)) u_SrvD%MsrPositionsY = 0.0
      if (allocated(u_SrvD%MsrPositionsz)) u_SrvD%MsrPositionsz = 0.0
      u_SrvD%YawErr = u_SrvD%WindDir - u_SrvD%YawAngle ! the nacelle yaw error estimate (positive about zi-axis)

!-------------------------------------------------------------------------------
! Unknown Mapping
!-------------------------------------------------------------------------------

   case default

      ErrStat = ErrID_Fatal
      ErrMsg = "Custom_InputSolve: unknown mapping '"//trim(Mapping%Desc)//"'"

   end select

end subroutine

subroutine SumMeshLoads(SrcMesh, DstMesh, DstResetFlag)
   type(MeshType), intent(in)    :: SrcMesh
   type(MeshType), intent(inout) :: DstMesh
   logical, intent(inout)        :: DstResetFlag
   if (DstResetFlag) then
      DstResetFlag = .false.
      if (DstMesh%fieldmask(MASKID_FORCE)) DstMesh%Force = 0.0_ReKi
      if (DstMesh%fieldmask(MASKID_MOMENT)) DstMesh%Moment = 0.0_ReKi
   end if
   if (DstMesh%fieldmask(MASKID_FORCE)) DstMesh%Force = DstMesh%Force + SrcMesh%Force
   if (DstMesh%fieldmask(MASKID_MOMENT)) DstMesh%Moment = DstMesh%Moment + SrcMesh%Moment
end subroutine

subroutine FAST_ResetRemapFlags(Mods, Maps, T, ErrStat, ErrMsg)
   type(ModDataType), intent(in)           :: Mods(:) !< Module data
   type(TC_MappingType), intent(inout)     :: Maps(:)
   type(FAST_TurbineType), intent(inout)   :: T       !< Turbine type
   integer(IntKi), intent(out)             :: ErrStat
   character(*), intent(out)               :: ErrMsg

   character(*), parameter    :: RoutineName = 'FAST_ResetRemapFlags'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2
   integer(IntKi)             :: i, k

   ErrStat = ErrID_None
   ErrMsg = ''

   ! Reset remap flags in mapping temporary meshes
   do i = 1, size(Maps)
      if (associated(Maps(i)%TmpLoadMesh%RemapFlag)) Maps(i)%TmpLoadMesh%RemapFlag = .false.
   end do

   do i = 1, size(Mods)

      ! Select based on module ID
      select case (Mods(i)%ID)

      case (Module_AD)

         if (T%AD%Input(1)%rotors(1)%HubMotion%Committed) then
            T%AD%Input(1)%rotors(1)%HubMotion%RemapFlag = .false.
            T%AD%y%rotors(1)%HubLoad%RemapFlag = .false.
         end if

         if (T%AD%Input(1)%rotors(1)%TowerMotion%Committed) then
            T%AD%Input(1)%rotors(1)%TowerMotion%RemapFlag = .false.

            if (T%AD%y%rotors(1)%TowerLoad%Committed) then
               T%AD%y%rotors(1)%TowerLoad%RemapFlag = .false.
            end if
         end if

         if (T%AD%Input(1)%rotors(1)%NacelleMotion%Committed) then
            T%AD%Input(1)%rotors(1)%NacelleMotion%RemapFlag = .false.
            T%AD%y%rotors(1)%NacelleLoad%RemapFlag = .false.
         end if

         if (T%AD%Input(1)%rotors(1)%TFinMotion%Committed) then
            T%AD%Input(1)%rotors(1)%TFinMotion%RemapFlag = .false.
            T%AD%y%rotors(1)%TFinLoad%RemapFlag = .false.
         end if

         do k = 1, size(T%AD%Input(1)%rotors(1)%BladeMotion)
            T%AD%Input(1)%rotors(1)%BladeRootMotion(k)%RemapFlag = .false.
            T%AD%Input(1)%rotors(1)%BladeMotion(k)%RemapFlag = .false.
            T%AD%y%rotors(1)%BladeLoad(k)%RemapFlag = .false.
         end do

      case (Module_BD)

         T%BD%Input(1, Mods(i)%Ins)%RootMotion%RemapFlag = .false.
         T%BD%Input(1, Mods(i)%Ins)%PointLoad%RemapFlag = .false.
         T%BD%Input(1, Mods(i)%Ins)%DistrLoad%RemapFlag = .false.
         T%BD%Input(1, Mods(i)%Ins)%HubMotion%RemapFlag = .false.

         T%BD%y(Mods(i)%Ins)%ReactionForce%RemapFlag = .false.
         T%BD%y(Mods(i)%Ins)%BldMotion%RemapFlag = .false.

      case (Module_ED)

         T%ED%Input(1)%PlatformPtMesh%RemapFlag = .false.
         T%ED%y%PlatformPtMesh%RemapFlag = .false.
         T%ED%Input(1)%TowerPtLoads%RemapFlag = .false.
         T%ED%y%TowerLn2Mesh%RemapFlag = .false.
         do K = 1, size(T%ED%y%BladeRootMotion)
            T%ED%y%BladeRootMotion(K)%RemapFlag = .false.
         end do
         if (allocated(T%ED%Input(1)%BladePtLoads)) then
            do K = 1, size(T%ED%Input(1)%BladePtLoads)
               T%ED%Input(1)%BladePtLoads(K)%RemapFlag = .false.
               T%ED%y%BladeLn2Mesh(K)%RemapFlag = .false.
            end do
         end if
         T%ED%Input(1)%NacelleLoads%RemapFlag = .false.
         T%ED%y%NacelleMotion%RemapFlag = .false.
         T%ED%Input(1)%TFinCMLoads%RemapFlag = .false.
         T%ED%y%TFinCMMotion%RemapFlag = .false.
         T%ED%Input(1)%HubPtLoad%RemapFlag = .false.
         T%ED%y%HubPtMotion%RemapFlag = .false.

      case (Module_ExtPtfm)

         if (T%ExtPtfm%Input(1)%PtfmMesh%Committed) then
            T%ExtPtfm%Input(1)%PtfmMesh%RemapFlag = .false.
            T%ExtPtfm%y%PtfmMesh%RemapFlag = .false.
         end if

      case (Module_FEAM)

         T%FEAM%Input(1)%PtFairleadDisplacement%RemapFlag = .false.
         T%FEAM%y%PtFairleadLoad%RemapFlag = .false.

      case (Module_HD)

         T%HD%Input(1)%PRPMesh%RemapFlag = .false.
         if (T%HD%Input(1)%WAMITMesh%Committed) then
            T%HD%Input(1)%WAMITMesh%RemapFlag = .false.
            T%HD%y%WAMITMesh%RemapFlag = .false.
         end if
         if (T%HD%Input(1)%Morison%Mesh%Committed) then
            T%HD%Input(1)%Morison%Mesh%RemapFlag = .false.
            T%HD%y%Morison%Mesh%RemapFlag = .false.
         end if

      case (Module_IceD)

         if (T%IceD%Input(1, Mods(i)%Ins)%PointMesh%Committed) then
            T%IceD%Input(1, Mods(i)%Ins)%PointMesh%RemapFlag = .false.
            T%IceD%y(Mods(i)%Ins)%PointMesh%RemapFlag = .false.
         end if

      case (Module_IceF)

         if (T%IceF%Input(1)%iceMesh%Committed) then
            T%IceF%Input(1)%iceMesh%RemapFlag = .false.
            T%IceF%y%iceMesh%RemapFlag = .false.
         end if

      case (Module_MAP)

         T%MAP%Input(1)%PtFairDisplacement%RemapFlag = .false.
         T%MAP%y%PtFairleadLoad%RemapFlag = .false.

      case (Module_MD)

         T%MD%Input(1)%CoupledKinematics(1)%RemapFlag = .false.
         T%MD%y%CoupledLoads(1)%RemapFlag = .false.

      case (Module_Orca)

         T%Orca%Input(1)%PtfmMesh%RemapFlag = .false.
         T%Orca%y%PtfmMesh%RemapFlag = .false.

      case (Module_SD)

         if (T%SD%Input(1)%TPMesh%Committed) then
            T%SD%Input(1)%TPMesh%RemapFlag = .false.
            T%SD%y%Y1Mesh%RemapFlag = .false.
         end if

         if (T%SD%Input(1)%LMesh%Committed) then
            T%SD%Input(1)%LMesh%RemapFlag = .false.
            T%SD%y%Y2Mesh%RemapFlag = .false.
            T%SD%y%Y3Mesh%RemapFlag = .false.
         end if

      end select

   end do

end subroutine

end module
