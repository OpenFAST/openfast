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
use ExtLoads

implicit none

private
public :: FAST_InitMappings
public :: FAST_LinearizeMappings
public :: FAST_ResetRemapFlags
public :: FAST_InputSolve
public :: FAST_ResetMappingReady
public :: FAST_InputFieldName, FAST_OutputFieldName

integer(IntKi), parameter  :: Xfr_Invalid = 0, &
                              Xfr_Point_to_Point = 1, &
                              Xfr_Line2_to_Point = 2, &
                              Xfr_Point_to_Line2 = 3, &
                              Xfr_Line2_to_Line2 = 4

character(24), parameter   :: Custom_ED_to_ExtLd = 'ED -> ExtLd', &
                              Custom_SrvD_to_AD = 'SrvD -> AD', &
                              Custom_ED_to_ADsk = 'ED -> ADsk', &
                              Custom_SED_to_ADsk = 'SED -> ADsk', &
                              Custom_ED_to_IfW = 'ED -> IfW', &
                              Custom_SrvD_to_IfW = 'SrvD -> IfW', &
                              Custom_BD_to_SrvD = 'BD -> SrvD', &
                              Custom_ED_to_SrvD = 'ED -> SrvD', &
                              Custom_SrvD_to_ED = 'SrvD -> ED', &
                              Custom_SrvD_to_SED = 'SrvD -> SED', &
                              Custom_IfW_to_SrvD = 'IfW -> SrvD', &
                              Custom_ExtInfw_to_SrvD = 'ExtInfw -> SrvD', &
                              Custom_SrvD_to_SD = 'SrvD -> SD', &
                              Custom_SrvD_to_MD = 'SrvD -> MD', &
                              Custom_ExtLd_to_SrvD = 'ExtLd -> SrvD'

contains

subroutine FAST_InputMeshPointer(ModData, Turbine, MeshLoc, Mesh, iInput, ErrStat, ErrMsg)
   type(ModDataType), intent(in)                :: ModData
   type(DatLoc), intent(in)                     :: MeshLoc
   type(FAST_TurbineType), target, intent(in)   :: Turbine
   type(MeshType), pointer, intent(out)         :: Mesh
   integer(IntKi), intent(in)                   :: iInput
   integer(IntKi), intent(out)                  :: ErrStat
   character(*), intent(out)                    :: ErrMsg

   ErrStat = ErrID_None
   ErrMsg = ""

   nullify (Mesh)

   select case (ModData%ID)
   case (Module_AD)
      Mesh => AD_InputMeshPointer(Turbine%AD%Input(iInput)%rotors(ModData%Ins), MeshLoc)
   case (Module_ADsk)
      Mesh => ADsk_InputMeshPointer(Turbine%ADsk%Input(iInput), MeshLoc)
   case (Module_BD)
      Mesh => BD_InputMeshPointer(Turbine%BD%Input(iInput, ModData%Ins), MeshLoc)
   case (Module_ED)
      Mesh => ED_InputMeshPointer(Turbine%ED%Input(iInput), MeshLoc)
   case (Module_SED)
      Mesh => SED_InputMeshPointer(Turbine%SED%Input(iInput), MeshLoc)
   case (Module_ExtInfw)
      ! ExtInfw doesn't have the typical input structure, using u
      Mesh => ExtInfw_InputMeshPointer(Turbine%ExtInfw%u, MeshLoc)
   case (Module_ExtLd)
      ! ExtLd doesn't have the typical input structure, using u
      Mesh => ExtLd_InputMeshPointer(Turbine%ExtLd%u, MeshLoc)
   case (Module_ExtPtfm)
      Mesh => ExtPtfm_InputMeshPointer(Turbine%ExtPtfm%Input(iInput), MeshLoc)
   case (Module_FEAM)
      Mesh => FEAM_InputMeshPointer(Turbine%FEAM%Input(iInput), MeshLoc)
   case (Module_HD)
      Mesh => HydroDyn_InputMeshPointer(Turbine%HD%Input(iInput), MeshLoc)
   case (Module_IceD)
      Mesh => IceD_InputMeshPointer(Turbine%IceD%Input(iInput, ModData%Ins), MeshLoc)
   case (Module_IceF)
      Mesh => IceFloe_InputMeshPointer(Turbine%IceF%Input(iInput), MeshLoc)
   case (Module_IfW)
      Mesh => InflowWind_InputMeshPointer(Turbine%IfW%Input(iInput), MeshLoc)
   case (Module_MAP)
      Mesh => MAP_InputMeshPointer(Turbine%MAP%Input(iInput), MeshLoc)
   case (Module_MD)
      Mesh => MD_InputMeshPointer(Turbine%MD%Input(iInput), MeshLoc)
   case (Module_Orca)
      Mesh => Orca_InputMeshPointer(Turbine%Orca%Input(iInput), MeshLoc)
   case (Module_SD)
      Mesh => SD_InputMeshPointer(Turbine%SD%Input(iInput), MeshLoc)
   case (Module_SeaSt)
      Mesh => SeaSt_InputMeshPointer(Turbine%SeaSt%Input(iInput), MeshLoc)
   case (Module_SrvD)
      Mesh => SrvD_InputMeshPointer(Turbine%SrvD%Input(iInput), MeshLoc)
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
   type(ModDataType), intent(in)                   :: ModData
   type(DatLoc), intent(in)                        :: MeshLoc
   type(FAST_TurbineType), target, intent(inout)   :: Turbine
   type(MeshType), pointer, intent(out)            :: Mesh
   integer(IntKi), intent(out)                     :: ErrStat
   character(*), intent(out)                       :: ErrMsg

   ErrStat = ErrID_None
   ErrMsg = ""

   nullify (Mesh)

   select case (ModData%ID)
   case (Module_AD)
      Mesh => AD_OutputMeshPointer(Turbine%AD%y%rotors(ModData%Ins), MeshLoc)
   case (Module_ADsk)
      Mesh => ADsk_OutputMeshPointer(Turbine%ADsk%y, MeshLoc)
   case (Module_BD)
      Mesh => BD_OutputMeshPointer(Turbine%BD%y(ModData%Ins), MeshLoc)
   case (Module_ED)
      Mesh => ED_OutputMeshPointer(Turbine%ED%y, MeshLoc)
   case (Module_SED)
      Mesh => SED_OutputMeshPointer(Turbine%SED%y, MeshLoc)
   case (Module_ExtInfw)
      Mesh => ExtInfw_OutputMeshPointer(Turbine%ExtInfw%y, MeshLoc)
   case (Module_ExtLd)
      Mesh => ExtLd_OutputMeshPointer(Turbine%ExtLd%y, MeshLoc)
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

function FAST_InputFieldName(ModData, DL) result(Name)
   type(ModDataType), intent(in) :: ModData
   type(DatLoc), intent(in)      :: DL
   character(42)                 :: Name, tmp
   select case (ModData%ID)
   case (Module_AD)
      Name = trim(ModData%Abbr)//"%u%rotors("//trim(Num2LStr(ModData%Ins))//")"
      select case (DL%Num)
      case (1:)
         tmp = AD_InputFieldName(DL)
         Name = trim(Name)//tmp(2:)
      case (AD_u_HWindSpeed)
         Name = 'AD%u%HWindSpeed (Ext)'
      case (AD_u_PLExp)
         Name = 'AD%u%PLExp (Ext)'
      case (AD_u_PropagationDir)
         Name = 'AD%u%PropagationDir (Ext)'
      end select
   case (Module_ADsk)
      Name = trim(ModData%Abbr)//"%"//ADsk_InputFieldName(DL)
   case (Module_BD)
      Name = trim(ModData%Abbr)//"("//trim(Num2LStr(ModData%Ins))//")%"//BD_InputFieldName(DL)
   case (Module_ED)
      select case (DL%Num)
      case (1:)
         Name = trim(ModData%Abbr)//"%"//ED_InputFieldName(DL)
      case (ED_u_BlPitchComC)
         Name = 'ED%u%BlPitchComC (Ext)'
      end select
   case (Module_SED)
      Name = trim(ModData%Abbr)//"%"//SED_InputFieldName(DL)
   case (Module_ExtInfw)
      Name = trim(ModData%Abbr)//"%"//ExtInfw_InputFieldName(DL)
   case (Module_ExtLd)
      Name = trim(ModData%Abbr)//"%"//ExtLd_InputFieldName(DL)
   case (Module_ExtPtfm)
      Name = trim(ModData%Abbr)//"%"//ExtPtfm_InputFieldName(DL)
   case (Module_FEAM)
      Name = trim(ModData%Abbr)//"%"//FEAM_InputFieldName(DL)
   case (Module_HD)
      select case (DL%Num)
      case (1:)
         Name = trim(ModData%Abbr)//"%"//HydroDyn_InputFieldName(DL)
      case (HydroDyn_u_WaveElev0)
         Name = 'HD%u%WaveElev0 (Ext)'
      case (HydroDyn_u_HWindSpeed)
         Name = 'HD%u%HWindSpeed (Ext)'
      case (HydroDyn_u_PLexp)
         Name = 'HD%u%PLexp (Ext)'
      case (HydroDyn_u_PropagationDir)
         Name = 'HD%u%PropagationDir (Ext)'
      end select
   case (Module_IceD)
      Name = trim(ModData%Abbr)//"("//trim(Num2LStr(ModData%Ins))//")%"//IceD_InputFieldName(DL)
   case (Module_IceF)
      Name = trim(ModData%Abbr)//"%"//IceFloe_InputFieldName(DL)
   case (Module_IfW)
      select case (DL%Num)
      case (1:)
         Name = trim(ModData%Abbr)//"%"//InflowWind_InputFieldName(DL)
      case (InflowWind_u_HWindSpeed)
         Name = 'IfW%u%HWindSpeed (Ext)'
      case (InflowWind_u_PLexp)
         Name = 'IfW%u%PLexp (Ext)'
      case (InflowWind_u_PropagationDir)
         Name = 'IfW%u%PropagationDir (Ext)'
      end select
   case (Module_MAP)
      Name = trim(ModData%Abbr)//"%"//MAP_InputFieldName(DL)
   case (Module_MD)
      Name = trim(ModData%Abbr)//"%"//MD_InputFieldName(DL)
   case (Module_Orca)
      Name = trim(ModData%Abbr)//"%"//Orca_InputFieldName(DL)
   case (Module_SD)
      Name = trim(ModData%Abbr)//"%"//SD_InputFieldName(DL)
   case (Module_SeaSt)
      select case (DL%Num)
      case (1:)
         Name = trim(ModData%Abbr)//"%"//SeaSt_InputFieldName(DL)
      case (SeaSt_u_WaveElev0)
         Name = 'SeaSt%u%WaveElev0 (Ext)'
      end select
   case (Module_SrvD)
      Name = trim(ModData%Abbr)//"%"//SrvD_InputFieldName(DL)
   case default
      Name = "Unknown field "//Num2LStr(DL%Num)//" in "//ModData%Abbr
   end select
end function

function FAST_OutputFieldName(ModData, DL) result(Name)
   type(ModDataType), intent(in)    :: ModData
   type(DatLoc), intent(in)         :: DL
   character(42)                    :: Name, tmp
   select case (ModData%ID)
   case (Module_AD)
      tmp = AD_OutputFieldName(DL)
      Name = trim(ModData%Abbr)//"%y%rotors("//trim(Num2LStr(ModData%Ins))//")"//tmp(2:)
   case (Module_ADsk)
      Name = trim(ModData%Abbr)//"%"//ADsk_OutputFieldName(DL)
   case (Module_BD)
      Name = trim(ModData%Abbr)//"("//trim(Num2LStr(ModData%Ins))//")%"//BD_OutputFieldName(DL)
   case (Module_ED)
      Name = trim(ModData%Abbr)//"%"//ED_OutputFieldName(DL)
   case (Module_SED)
      Name = trim(ModData%Abbr)//"%"//SED_OutputFieldName(DL)
   case (Module_ExtInfw)
      Name = trim(ModData%Abbr)//"%"//ExtInfw_OutputFieldName(DL)
   case (Module_ExtLd)
      Name = trim(ModData%Abbr)//"%"//ExtLd_OutputFieldName(DL)
   case (Module_ExtPtfm)
      Name = trim(ModData%Abbr)//"%"//ExtPtfm_OutputFieldName(DL)
   case (Module_FEAM)
      Name = trim(ModData%Abbr)//"%"//FEAM_OutputFieldName(DL)
   case (Module_HD)
      Name = trim(ModData%Abbr)//"%"//HydroDyn_OutputFieldName(DL)
   case (Module_IceD)
      Name = trim(ModData%Abbr)//"("//trim(Num2LStr(ModData%Ins))//")%"//IceD_OutputFieldName(DL)
   case (Module_IceF)
      Name = trim(ModData%Abbr)//"%"//IceFloe_OutputFieldName(DL)
   case (Module_IfW)
      select case (DL%Num)
      case (1:)
         Name = trim(ModData%Abbr)//"%"//InflowWind_OutputFieldName(DL)
      case (InflowWind_y_HWindSpeed)
         Name = 'IfW%y%HWindSpeed (Ext)'
      case (InflowWind_y_PLexp)
         Name = 'IfW%y%PLexp (Ext)'
      case (InflowWind_y_PropagationDir)
         Name = 'IfW%y%PropagationDir (Ext)'
      end select
   case (Module_MAP)
      Name = trim(ModData%Abbr)//"%"//MAP_OutputFieldName(DL)
   case (Module_MD)
      Name = trim(ModData%Abbr)//"%"//MD_OutputFieldName(DL)
   case (Module_Orca)
      Name = trim(ModData%Abbr)//"%"//Orca_OutputFieldName(DL)
   case (Module_SD)
      Name = trim(ModData%Abbr)//"%"//SD_OutputFieldName(DL)
   case (Module_SeaSt)
      select case (DL%Num)
      case (1:)
         Name = trim(ModData%Abbr)//"%"//SeaSt_OutputFieldName(DL)
      case (SeaSt_y_WaveElev0)
         Name = 'SeaSt%y%WaveElev0 (Ext)'
      end select
   case (Module_SrvD)
      Name = trim(ModData%Abbr)//"%"//SrvD_OutputFieldName(DL)
   case default
      Name = "Unknown field "//Num2LStr(DL%Num)//" in "//ModData%Abbr
   end select
end function

subroutine FAST_InitMappings(Mappings, Mods, Turbine, ErrStat, ErrMsg)
   type(MappingType), allocatable, intent(out)        :: Mappings(:)
   type(ModDataType), intent(inout)                   :: Mods(:)     !< Module data
   type(FAST_TurbineType), intent(inout)              :: Turbine     !< Turbine type
   integer(IntKi), intent(out)                        :: ErrStat
   character(*), intent(out)                          :: ErrMsg

   character(*), parameter          :: RoutineName = 'FAST_InitMappings'
   integer(IntKi)                   :: ErrStat2
   character(ErrMsgLen)             :: ErrMsg2
   integer(IntKi)                   :: i, j, k
   integer(IntKi)                   :: iMap, ModIns, iModIn, iModSrc, iModDst
   type(MappingType), allocatable   :: MappingsTmp(:)
   integer(IntKi), parameter        :: MappingTypeOrder(*) = [Map_MotionMesh, Map_LoadMesh, Map_Variable, Map_Custom]

   ErrStat = ErrID_None
   ErrMsg = ''

   !----------------------------------------------------------------------------
   ! Define mesh mappings between modules
   !----------------------------------------------------------------------------

   ! Define a list of all possible module mesh mappings between modules
   allocate (MappingsTmp(0), stat=ErrStat2)
   if (ErrStat2 /= 0) then
      call SetErrStat(ErrID_Fatal, "Error allocating temporary mappings", ErrStat, ErrMsg, RoutineName)
      return
   end if

   ! Loop through module pairings
   do iModSrc = 1, size(Mods)
      do iModDst = 1, size(Mods)

         ! Switch by destination module (inputs)
         select case (Mods(IModDst)%ID)
         case (Module_AD)
            call InitMappings_AD(MappingsTmp, Mods(iModSrc), Mods(iModDst), Turbine, ErrStat2, ErrMsg2)
         case (Module_ADsk)
            call InitMappings_ADsk(MappingsTmp, Mods(iModSrc), Mods(iModDst), Turbine, ErrStat2, ErrMsg2)
         case (Module_BD)
            call InitMappings_BD(MappingsTmp, Mods(iModSrc), Mods(iModDst), Turbine, ErrStat2, ErrMsg2)
         case (Module_ED)
            call InitMappings_ED(MappingsTmp, Mods(iModSrc), Mods(iModDst), Turbine, ErrStat2, ErrMsg2)
         case (Module_SED)
            call InitMappings_SED(MappingsTmp, Mods(iModSrc), Mods(iModDst), Turbine, ErrStat2, ErrMsg2)
         case (Module_ExtInfw)
            call InitMappings_ExtInfw(MappingsTmp, Mods(iModSrc), Mods(iModDst), Turbine, ErrStat2, ErrMsg2)
         case (Module_ExtLd)
            call InitMappings_ExtLd(MappingsTmp, Mods(iModSrc), Mods(iModDst), Turbine, ErrStat2, ErrMsg2)
         case (Module_ExtPtfm)
            call InitMappings_ExtPtfm(MappingsTmp, Mods(iModSrc), Mods(iModDst), Turbine, ErrStat2, ErrMsg2)
         case (Module_FEAM)
            call InitMappings_FEAM(MappingsTmp, Mods(iModSrc), Mods(iModDst), Turbine, ErrStat2, ErrMsg2)
         case (Module_HD)
            call InitMappings_HD(MappingsTmp, Mods(iModSrc), Mods(iModDst), Turbine, ErrStat2, ErrMsg2)
         case (Module_IceD)
            call InitMappings_IceD(MappingsTmp, Mods(iModSrc), Mods(iModDst), Turbine, ErrStat2, ErrMsg2)
         case (Module_IceF)
            call InitMappings_IceF(MappingsTmp, Mods(iModSrc), Mods(iModDst), Turbine, ErrStat2, ErrMsg2)
         case (Module_IfW)
            call InitMappings_IfW(MappingsTmp, Mods(iModSrc), Mods(iModDst), Turbine, ErrStat2, ErrMsg2)
         case (Module_MAP)
            call InitMappings_MAP(MappingsTmp, Mods(iModSrc), Mods(iModDst), Turbine, ErrStat2, ErrMsg2)
         case (Module_MD)
            call InitMappings_MD(MappingsTmp, Mods(iModSrc), Mods(iModDst), Turbine, ErrStat2, ErrMsg2)
         case (Module_Orca)
            call InitMappings_Orca(MappingsTmp, Mods(iModSrc), Mods(iModDst), Turbine, ErrStat2, ErrMsg2)
         case (Module_SD)
            call InitMappings_SD(MappingsTmp, Mods(iModSrc), Mods(iModDst), Turbine, ErrStat2, ErrMsg2)
         case (Module_SeaSt)
            call InitMappings_SeaSt(MappingsTmp, Mods(iModSrc), Mods(iModDst), Turbine, ErrStat2, ErrMsg2)
         case (Module_SrvD)
            call InitMappings_SrvD(MappingsTmp, Mods(iModSrc), Mods(iModDst), Turbine, ErrStat2, ErrMsg2)
         end select
         if (Failed()) return
      end do
   end do

   !----------------------------------------------------------------------------
   ! Reorder mappings to be Motion, Load, Variable, Custom
   !----------------------------------------------------------------------------

   allocate(Mappings(size(MappingsTmp)), stat=ErrStat2)
   if (ErrStat2 /= 0) then
      call SetErrStat(ErrID_Fatal, "Error allocating mappings", ErrStat, ErrMsg, RoutineName)
      return
   end if

   ! Loop through MappingTypeOrder and copy mesh to Mappings array if it matches the type
   k = 0
   do i = 1, size(MappingTypeOrder)
      do j = 1, size(MappingsTmp)
         if (MappingsTmp(j)%MapType == MappingTypeOrder(i)) then
            k = k + 1
            call Glue_CopyMappingType(MappingsTmp(j), Mappings(k), MESH_NEWCOPY, ErrStat2, ErrMsg2)
            if (Failed()) return
         end if
      end do
   end do

   ! Destroy temporary mappings
   do i = 1, size(MappingsTmp)
      call Glue_DestroyMappingType(MappingsTmp(i), ErrStat2, ErrMsg2)
      if (Failed()) return
   end do

   ! Loop through mappings
   do iMap = 1, size(Mappings)
      associate (SrcMod => Mods(Mappings(iMap)%iModSrc), &
                 DstMod => Mods(Mappings(iMap)%iModDst))

         write (*, *) "Mapping: ", Mappings(iMap)%Desc

      end associate
   end do

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine InitMappings_AD(Mappings, SrcMod, DstMod, Turbine, ErrStat, ErrMsg)
   type(MappingType), allocatable         :: Mappings(:)
   type(ModDataType), intent(inout)       :: SrcMod, DstMod
   type(FAST_TurbineType), intent(inout)  :: Turbine
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg

   character(*), parameter    :: RoutineName = 'InitMappings_AD'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2
   integer(IntKi)             :: i
   logical                    :: NotCompAeroMaps, CompElastED

   ErrStat = ErrID_None
   ErrMsg = ''

   ! Flag is true if not computing AeroMaps
   NotCompAeroMaps = .not. Turbine%p_FAST%CompAeroMaps

   ! Flag is true if CompElast == Module_ED
   CompElastED = Turbine%p_FAST%CompElast == Module_ED

   ! Select based on source module identifier
   select case (SrcMod%ID)

   case (Module_BD)

      call MapMotionMesh(Turbine, Mappings, &
                         SrcMod=SrcMod, SrcDL=DatLoc(BD_y_BldMotion), &                                    ! BD%y(SrcMod%Ins)%BldMotion
                         DstMod=DstMod, DstDL=DatLoc(AD_u_BladeMotion, SrcMod%Ins), &   ! AD%u%rotors(DstMod%Ins)%BladeMotion(SrcMod%Ins)
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2, &
                         Active=NotCompAeroMaps .or. (SrcMod%Ins == 1))
      if (Failed()) return

   case (Module_ED)

      ! Blade motion
      if (Turbine%p_FAST%CompElast == Module_ED) then
         do i = 1, size(Turbine%ED%y%BladeLn2Mesh)
            call MapMotionMesh(Turbine, Mappings, &
                               SrcMod=SrcMod, SrcDL=DatLoc(ED_y_BladeLn2Mesh, i), &                     ! ED%y%BladeLn2Mesh(i)
                               DstMod=DstMod, DstDL=DatLoc(AD_u_BladeMotion, i), &   ! AD%u%rotors(DstMod%Ins)%BladeMotion(i)
                               ErrStat=ErrStat2, ErrMsg=ErrMsg2, &
                               Active=CompElastED .and. (NotCompAeroMaps .or. (i == 1)))
            if (Failed()) return
         end do
      end if

      ! Blade root motion
      do i = 1, size(Turbine%ED%y%BladeRootMotion)
         call MapMotionMesh(Turbine, Mappings, &
                            SrcMod=SrcMod, SrcDL=DatLoc(ED_y_BladeRootMotion, i), &        ! ED%y%BladeRootMotion(i)
                            DstMod=DstMod, DstDL=DatLoc(AD_u_BladeRootMotion, i), &        ! AD%u%rotors(DstMod%Ins)%BladeRootMotion(i)
                            ErrStat=ErrStat2, ErrMsg=ErrMsg2, &
                            Active=NotCompAeroMaps)
         if (Failed()) return
      end do

      ! Tower motion
      call MapMotionMesh(Turbine, Mappings, &
                         SrcMod=SrcMod, SrcDL=DatLoc(ED_y_TowerLn2Mesh), &                 ! ED%y%TowerLn2Mesh
                         DstMod=DstMod, DstDL=DatLoc(AD_u_TowerMotion), &                  ! AD%u%rotors(DstMod%Ins)%TowerMotion
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2, &
                         Active=NotCompAeroMaps)
      if (Failed()) return

      ! Hub motion
      call MapMotionMesh(Turbine, Mappings, &
                         SrcMod=SrcMod, SrcDL=DatLoc(ED_y_HubPtMotion), &                  ! ED%y%HubPtMotion
                         DstMod=DstMod, DstDL=DatLoc(AD_u_HubMotion), &                    ! AD%u%rotors(DstMod%Ins)%HubMotion
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2, &
                         Active=NotCompAeroMaps)
      if (Failed()) return

      ! Nacelle motion
      call MapMotionMesh(Turbine, Mappings, &
                         SrcMod=SrcMod, SrcDL=DatLoc(ED_y_NacelleMotion), &                ! ED%y%NacelleMotion
                         DstMod=DstMod, DstDL=DatLoc(AD_u_NacelleMotion), &                ! AD%u%rotors(DstMod%Ins)%NacelleMotion
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2, &
                         Active=NotCompAeroMaps)
      if (Failed()) return

      ! TailFin motion
      call MapMotionMesh(Turbine, Mappings, &
                         SrcMod=SrcMod, SrcDL=DatLoc(ED_y_TFinCMMotion), &                 ! ED%y%TFinCMMotion
                         DstMod=DstMod, DstDL=DatLoc(AD_u_TFinMotion), &                   ! AD%u%rotors(DstMod%Ins)%TFinMotion
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2, &
                         Active=NotCompAeroMaps)
      if (Failed()) return

   case (Module_IfW)

      call MapVariable(Mappings, &
                       SrcMod=SrcMod, SrcDL=DatLoc(InflowWind_y_HWindSpeed), &
                       DstMod=DstMod, DstDL=DatLoc(AD_u_HWindSpeed), &
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2, &
                       Active=Turbine%p_FAST%Linearize)
      if (Failed()) return

      call MapVariable(Mappings, &
                       SrcMod=SrcMod, SrcDL=DatLoc(InflowWind_y_PLExp), &
                       DstMod=DstMod, DstDL=DatLoc(AD_u_PLExp), &
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2, &
                       Active=Turbine%p_FAST%Linearize)
      if (Failed()) return

      call MapVariable(Mappings, &
                       SrcMod=SrcMod, SrcDL=DatLoc(InflowWind_y_PropagationDir), &
                       DstMod=DstMod, DstDL=DatLoc(AD_u_PropagationDir), &
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2, &
                       Active=Turbine%p_FAST%Linearize)
      if (Failed()) return

   case (Module_SrvD)

      call MapCustom(Mappings, Custom_SrvD_to_AD, SrcMod, DstMod)

   end select

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine InitMappings_ADsk(Mappings, SrcMod, DstMod, Turbine, ErrStat, ErrMsg)
   type(MappingType), allocatable         :: Mappings(:)
   type(ModDataType), intent(inout)       :: SrcMod, DstMod
   type(FAST_TurbineType), intent(inout)  :: Turbine
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg

   character(*), parameter    :: RoutineName = 'InitMappings_ADsk'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2

   ErrStat = ErrID_None
   ErrMsg = ''

   ! Select based on source module identifier
   select case (SrcMod%ID)

   case (Module_ED)

      ! Hub motion
      call MapMotionMesh(Turbine, Mappings, &
                         SrcMod=SrcMod, SrcDL=DatLoc(ED_y_HubPtMotion), &        ! ED%y%HubPtMotion
                         DstMod=DstMod, DstDL=DatLoc(ADsk_u_HubMotion), &        ! ADsk%u%HubMotion
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2)
      if (Failed()) return

      call MapCustom(Mappings, Custom_ED_to_ADsk, SrcMod, DstMod)

   case (Module_SED)

      ! Hub motion
      call MapMotionMesh(Turbine, Mappings, &
                         SrcMod=SrcMod, SrcDL=DatLoc(SED_y_HubPtMotion), &       ! ED%y%HubPtMotion
                         DstMod=DstMod, DstDL=DatLoc(ADsk_u_HubMotion), &        ! ADsk%u%HubMotion
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2)
      if (Failed()) return

      call MapCustom(Mappings, Custom_SED_to_ADsk, SrcMod, DstMod)

   end select

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine InitMappings_BD(Mappings, SrcMod, DstMod, Turbine, ErrStat, ErrMsg)
   type(MappingType), allocatable      :: Mappings(:)
   type(ModDataType), intent(inout)       :: SrcMod, DstMod
   type(FAST_TurbineType), intent(inout)  :: Turbine           !< Turbine type
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg

   character(*), parameter    :: RoutineName = 'InitMappings_BD'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2
   integer(IntKi)             :: i
   logical                    :: NotCompAeroMaps, CompAeroAD

   ErrStat = ErrID_None
   ErrMsg = ''

   ! Flag is true if not computing AeroMaps
   NotCompAeroMaps = .not. Turbine%p_FAST%CompAeroMaps

   ! Flag is true of CompAero == Module_AD
   CompAeroAD = Turbine%p_FAST%CompAero == Module_AD

   ! Select based on source module identifier
   select case (SrcMod%ID)

   case (Module_AD)

      call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                       SrcDL=DatLoc(AD_y_BladeLoad, DstMod%Ins), &          ! AD%y%rotors(SrcMod%Ins)%BladeLoad(DstMod%Ins)
                       SrcDispDL=DatLoc(AD_u_BladeMotion, DstMod%Ins), &    ! AD%u%rotors(SrcMod%Ins)%BladeMotion(DstMod%Ins)
                       DstDL=DatLoc(BD_u_DistrLoad), &                      ! BD%u(DstMod%Ins)%DistrLoad
                       DstDispDL=DatLoc(BD_y_BldMotion), &                  ! BD%y(DstMod%Ins)%BldMotion
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2, &
                       Active=CompAeroAD .and. (NotCompAeroMaps .or. (DstMod%Ins == 1)))
      if (Failed()) return

   case (Module_ED)

      call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                         SrcDL=DatLoc(ED_y_BladeRootMotion, DstMod%Ins), &   ! ED%y%BladeRootMotion(DstMod%Ins)
                         DstDL=DatLoc(BD_u_RootMotion), &                    ! BD%u(DstMod%Ins)%RootMotion
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2, &
                         Active=NotCompAeroMaps)
      if (Failed()) return

      ! Hub motion not used
      ! call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
      !                    SrcDL=DatLoc(ED_y_HubPtMotion), &                   ! ED%y%HubED_y_HubPtMotion
      !                    DstDL=DatLoc(BD_u_HubMotion), &                     ! BD%Input(1, DstMod%Ins)%HubMotion
      !                    ErrStat=ErrStat2, ErrMsg=ErrMsg2, &
      !                    Active=NotCompAeroMaps)
      ! if (Failed()) return

   case (Module_ExtLd)

      call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                        SrcDL=DatLoc(ExtLd_y_BladeLoad, DstMod%Ins), &           ! ExtLd%y%BladeLoad(DstMod%Ins), &
                        SrcDispDL=DatLoc(ExtLd_u_BladeMotion, DstMod%Ins), &     ! ExtLd%u%BladeMotion(DstMod%Ins)
                        DstDL=DatLoc(BD_u_DistrLoad), &                          ! BD%Input(1, DstMod%Ins)%DistrLoad
                        DstDispDL=DatLoc(BD_y_BldMotion), &                      ! BD%y(DstMod%Ins)%BldMotion
                        ErrStat=ErrStat2, ErrMsg=ErrMsg2)
      if (Failed()) return

   case (Module_SrvD)

      do i = 1, Turbine%SrvD%p%NumBStC
         call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                          SrcDL=DatLoc(SrvD_y_BStCLoadMesh, DstMod%Ins, i), &        ! SrvD%y%BStCLoadMesh(DstMod%Ins, i), &
                          SrcDispDL=DatLoc(SrvD_u_BStCMotionMesh, DstMod%Ins, i), &  ! SrvD%u%BStCMotionMesh(DstMod%Ins, i)
                          DstDL=DatLoc(BD_u_DistrLoad), &                            ! BD%Input(1, DstMod%Ins)%DistrLoad
                          DstDispDL=DatLoc(BD_y_BldMotion), &                        ! BD%y(DstMod%Ins)%BldMotion
                          ErrStat=ErrStat2, ErrMsg=ErrMsg2)
         if (Failed()) return
      end do

   end select

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine InitMappings_ED(Mappings, SrcMod, DstMod, Turbine, ErrStat, ErrMsg)
   type(MappingType), allocatable      :: Mappings(:)
   type(ModDataType), intent(inout)       :: SrcMod, DstMod
   type(FAST_TurbineType), intent(inout)  :: Turbine           !< Turbine type
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg

   character(*), parameter    :: RoutineName = 'InitMappings_ED'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2
   integer(IntKi)             :: i, j
   logical                    :: NotCompAeroMaps, CompAeroAD, CompElastED

   ErrStat = ErrID_None
   ErrMsg = ''

   ! Flag is true if not computing AeroMaps
   NotCompAeroMaps = .not. Turbine%p_FAST%CompAeroMaps

   ! Flag is true of CompAero == Module_AD
   CompAeroAD = Turbine%p_FAST%CompAero == Module_AD

   ! Flag is true of CompElast == Module_ED
   CompElastED = Turbine%p_FAST%CompElast == Module_ED

   ! Select based on source module identifier
   select case (SrcMod%ID)

   case (Module_AD)

      ! Blade Loads
      do i = 1, Turbine%ED%p%NumBl
         call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                          SrcDL=DatLoc(AD_y_BladeLoad, i), &            ! AD%y%rotors(SrcMod%Ins)%BladeLoad(i)
                          SrcDispDL=DatLoc(AD_u_BladeMotion, i), &      ! AD%u%rotors(SrcMod%Ins)%BladeMotion(i)
                          DstDL=DatLoc(ED_u_BladePtLoads, i), &         ! ED%u%BladePtLoads(i)
                          DstDispDL=DatLoc(ED_y_BladeLn2Mesh, i), &     ! ED%y%BladeLn2Mesh(i)
                          ErrStat=ErrStat2, ErrMsg=ErrMsg2, &
                          Active=CompAeroAD .and. CompElastED .and. &
                                 (NotCompAeroMaps .or. (i == 1)))
         if (Failed()) return
      end do

      ! Hub Loads
      call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                       SrcDL=DatLoc(AD_y_HubLoad), &           ! AD%y%rotors(SrcMod%Ins)%HubLoad
                       SrcDispDL=DatLoc(AD_u_HubMotion), &     ! AD%u%rotors(SrcMod%Ins)%HubMotion
                       DstDL=DatLoc(ED_u_HubPtLoad), &         ! ED%u%HubPtLoad
                       DstDispDL=DatLoc(ED_y_HubPtMotion), &   ! ED%y%HubPtMotion
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2, &
                       Active=NotCompAeroMaps)
      if (Failed()) return

      ! Nacelle Loads
      call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                       SrcDL=DatLoc(AD_y_NacelleLoad), &       ! AD%y%rotors(SrcMod%Ins)%NacelleLoad
                       SrcDispDL=DatLoc(AD_u_NacelleMotion), & ! AD%u%rotors(SrcMod%Ins)%NacelleMotion
                       DstDL=DatLoc(ED_u_NacelleLoads), &      ! ED%u%NacelleLoads
                       DstDispDL=DatLoc(ED_y_NacelleMotion), & ! ED%y%NacelleMotion
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2, &
                       Active=NotCompAeroMaps)
      if (Failed()) return

      ! Tail Fin Loads
      call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                       SrcDL=DatLoc(AD_y_TFinLoad), &          ! AD%y%rotors(SrcMod%Ins)%TFinLoad
                       SrcDispDL=DatLoc(AD_u_TFinMotion), &    ! AD%u%rotors(SrcMod%Ins)%TFinMotion
                       DstDL=DatLoc(ED_u_TFinCMLoads), &       ! ED%u%TFinCMLoads
                       DstDispDL=DatLoc(ED_y_TFinCMMotion), &  ! ED%y%TFinCMMotion
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2, &
                       Active=NotCompAeroMaps)
      if (Failed()) return

      ! Tower Loads
      call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                       SrcDL=DatLoc(AD_y_TowerLoad), &         ! AD%y%rotors(SrcMod%Ins)%TowerLoad
                       SrcDispDL=DatLoc(AD_u_TowerMotion), &   ! AD%u%rotors(SrcMod%Ins)%TowerMotion
                       DstDL=DatLoc(ED_u_TowerPtLoads), &      ! ED%u%TowerPtLoads
                       DstDispDL=DatLoc(ED_y_TowerLn2Mesh), &  ! ED%y%TowerLn2Mesh
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2, &
                       Active=CompAeroAD .and. NotCompAeroMaps)
      if (Failed()) return

   case (Module_ADsk)

      ! Hub Loads
      call MapLoadMesh(Turbine, Mappings, &
                       SrcMod=SrcMod, &
                       SrcDL=DatLoc(ADsk_y_AeroLoads), &        ! ADsk%y%AeroLoads
                       SrcDispDL=DatLoc(ADsk_u_HubMotion), &    ! ADsk%u%HubMotion
                       DstMod=DstMod, &
                       DstDL=DatLoc(ED_u_HubPtLoad), &          ! ED%u%HubPtLoad
                       DstDispDL=DatLoc(ED_y_HubPtMotion), &    ! ED%y%HubPtMotion
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2)
      if (Failed()) return

   case (Module_BD)

      ! Hub Loads
      call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                       SrcDL=DatLoc(BD_y_ReactionForce), &      ! BD%y(SrcMod%Ins)%ReactionForce
                       SrcDispDL=DatLoc(BD_u_RootMotion), &     ! BD%u(SrcMod%Ins)%RootMotion
                       DstDL=DatLoc(ED_u_HubPtLoad), &          ! ED%u%HubPtLoad
                       DstDispDL=DatLoc(ED_y_HubPtMotion), &    ! ED%y%HubPtMotion
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2, &
                       Active=NotCompAeroMaps)
      if (Failed()) return

   case (Module_ExtLd)

      ! Blade loads
      do i = 1, Turbine%ED%p%NumBl
         call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                          SrcDL=DatLoc(ExtLd_y_BladeLoad, i), &            ! ExtLd%y%BladeLoad(i)
                          SrcDispDL=DatLoc(ExtLd_u_BladeMotion, i), &      ! ExtLd%u%BladeMotion(i)
                          DstDL=DatLoc(ED_u_BladePtLoads, i), &            ! ED%u%BladePtLoads(i)
                          DstDispDL=DatLoc(ED_y_BladeLn2Mesh, i), &        ! ED%y%BladeLn2Mesh(i)
                          ErrStat=ErrStat2, ErrMsg=ErrMsg2)
         if (Failed()) return
      end do

      ! Tower load
      call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                       SrcDL=DatLoc(ExtLd_y_TowerLoad), &                  ! ExtLd%y%TowerLoad
                       SrcDispDL=DatLoc(ExtLd_u_TowerMotion), &            ! ExtLd%u%TowerMotion
                       DstDL=DatLoc(ED_u_TowerPtLoads), &                  ! ED%u%TowerPtLoads
                       DstDispDL=DatLoc(ED_y_TowerLn2Mesh), &              ! ED%y%TowerLn2Mesh
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2)
      if (Failed()) return

   case (Module_ExtPtfm)

      call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                       SrcDL=DatLoc(ExtPtfm_y_PtfmMesh), &             ! ExtPtfm%y%PtfmMesh
                       SrcDispDL=DatLoc(ExtPtfm_u_PtfmMesh), &         ! ExtPtfm%u%PtfmMesh
                       DstDL=DatLoc(ED_u_PlatformPtMesh), &            ! ED%u%PlatformPtMesh
                       DstDispDL=DatLoc(ED_y_PlatformPtMesh), &        ! ED%y%PlatformPtMesh
                       Active=Turbine%p_FAST%CompSub /= Module_SD, &
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2)
      if (Failed()) return

   case (Module_FEAM)

      call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                       SrcDL=DatLoc(FEAM_y_PtFairleadLoad), &             ! FEAM%y%PtFairleadLoad, &
                       SrcDispDL=DatLoc(FEAM_u_PtFairleadDisplacement), & ! FEAM%u%PtFairleadDisplacement
                       DstDL=DatLoc(ED_u_PlatformPtMesh), &               ! ED%u%PlatformPtMesh
                       DstDispDL=DatLoc(ED_y_PlatformPtMesh), &           ! ED%y%PlatformPtMesh
                       Active=Turbine%p_FAST%CompSub /= Module_SD, &
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2)
      if (Failed()) return

   case (Module_HD)

      ! Platform loads (SubDyn not active)
      call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                       SrcDL=DatLoc(HydroDyn_y_Morison_Mesh), &        ! HD%y%Morison%Mesh
                       SrcDispDL=DatLoc(HydroDyn_u_Morison_Mesh), &    ! HD%u%Morison%Mesh
                       DstDL=DatLoc(ED_u_PlatformPtMesh), &            ! ED%u%PlatformPtMesh
                       DstDispDL=DatLoc(ED_y_PlatformPtMesh), &        ! ED%y%PlatformPtMesh
                       Active=Turbine%p_FAST%CompSub == Module_None, &
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2)
      if (Failed()) return

      ! Platform loads (SubDyn not active)
      call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                       SrcDL=DatLoc(HydroDyn_y_WAMITMesh), &           ! HD%y%WAMITMesh
                       SrcDispDL=DatLoc(HydroDyn_u_WAMITMesh), &       ! HD%u%WAMITMesh
                       DstDL=DatLoc(ED_u_PlatformPtMesh), &            ! ED%u%PlatformPtMesh
                       DstDispDL=DatLoc(ED_y_PlatformPtMesh), &        ! ED%y%PlatformPtMesh
                       Active=Turbine%p_FAST%CompSub == Module_None, &
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2)
      if (Failed()) return

   case (Module_IceD)

      ! Platform loads (SubDyn not active)
      call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                       SrcDL=DatLoc(IceD_y_PointMesh), &               ! IceD%y%PointMesh
                       SrcDispDL=DatLoc(IceD_u_PointMesh), &           ! IceD%u%PointMesh
                       DstDL=DatLoc(ED_u_PlatformPtMesh), &            ! ED%u%PlatformPtMesh
                       DstDispDL=DatLoc(ED_y_PlatformPtMesh), &        ! ED%y%PlatformPtMesh
                       Active=Turbine%p_FAST%CompSub /= Module_SD, &
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2)
      if (Failed()) return

   case (Module_IceF)

      ! Platform loads (SubDyn not active)
      call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                       SrcDL=DatLoc(IceFloe_y_iceMesh), &              ! IceFloe%y%iceMesh
                       SrcDispDL=DatLoc(IceFloe_u_iceMesh), &          ! IceFloe%u%iceMesh
                       DstDL=DatLoc(ED_u_PlatformPtMesh), &            ! ED%u%PlatformPtMesh
                       DstDispDL=DatLoc(ED_y_PlatformPtMesh), &        ! ED%y%PlatformPtMesh
                       Active=Turbine%p_FAST%CompSub /= Module_SD, &
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2)
      if (Failed()) return

   case (Module_MAP)

      ! Platform loads (SubDyn not active)
      call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                       SrcDL=DatLoc(MAP_y_ptFairleadLoad), &           ! MAP%y%PtFairleadLoad
                       SrcDispDL=DatLoc(MAP_u_PtFairDisplacement), &   ! MAP%u%PtFairDisplacement
                       DstDL=DatLoc(ED_u_PlatformPtMesh), &            ! ED%u%PlatformPtMesh
                       DstDispDL=DatLoc(ED_y_PlatformPtMesh), &        ! ED%y%PlatformPtMesh
                       Active=Turbine%p_FAST%CompSub /= Module_SD, &
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2)
      if (Failed()) return

   case (Module_MD)

      ! Platform loads (SubDyn not active)
      call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                       SrcDL=DatLoc(MD_y_CoupledLoads, 1), &           ! MD%y%CoupledLoads(1)
                       SrcDispDL=DatLoc(MD_u_CoupledKinematics, 1), &  ! MD%u%CoupledKinematics(1)
                       DstDL=DatLoc(ED_u_PlatformPtMesh), &            ! ED%u%PlatformPtMesh
                       DstDispDL=DatLoc(ED_y_PlatformPtMesh), &        ! ED%y%PlatformPtMesh
                       Active=Turbine%p_FAST%CompSub /= Module_SD, &
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2)
      if (Failed()) return

   case (Module_Orca)

      ! Platform loads (SubDyn not active)
      call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                       SrcDL=DatLoc(Orca_y_PtfmMesh), &                ! Orca%y%PtfmMesh
                       SrcDispDL=DatLoc(Orca_u_PtfmMesh), &            ! Orca%u%PtfmMesh
                       DstDL=DatLoc(ED_u_PlatformPtMesh), &            ! ED%u%PlatformPtMesh
                       DstDispDL=DatLoc(ED_y_PlatformPtMesh), &        ! ED%y%PlatformPtMesh
                       Active=Turbine%p_FAST%CompSub /= Module_SD, &
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2)
      if (Failed()) return

   case (Module_SD)

      call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                       SrcDL=DatLoc(SD_y_Y1Mesh), &                    ! SD%y%Y1mesh, &
                       SrcDispDL=DatLoc(SD_u_TPMesh), &                ! SD%u%TPMesh
                       DstDL=DatLoc(ED_u_PlatformPtMesh), &            ! ED%u%PlatformPtMesh
                       DstDispDL=DatLoc(ED_y_PlatformPtMesh), &        ! ED%y%PlatformPtMesh
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2)
      if (Failed()) return

   case (Module_SrvD)

      call MapVariable(Mappings, &
                       SrcMod=SrcMod, SrcDL=DatLoc(SrvD_y_BlPitchCom), &
                       DstMod=DstMod, DstDL=DatLoc(ED_u_BlPitchCom), &
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2); if (Failed()) return

      call MapVariable(Mappings, &
                       SrcMod=SrcMod, SrcDL=DatLoc(SrvD_y_YawMom), &
                       DstMod=DstMod, DstDL=DatLoc(ED_u_YawMom), &
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2); if (Failed()) return

      call MapVariable(Mappings, &
                       SrcMod=SrcMod, SrcDL=DatLoc(SrvD_y_GenTrq), &
                       DstMod=DstMod, DstDL=DatLoc(ED_u_GenTrq), &
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2); if (Failed()) return

      call MapCustom(Mappings, Custom_SrvD_to_ED, SrcMod, DstMod)

      ! Blade Structural Controller (if ElastoDyn is used for blades)
      do j = 1, Turbine%SrvD%p%NumBStC
         do i = 1, Turbine%ED%p%NumBl
            call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                             SrcDL=DatLoc(SrvD_y_BStCLoadMesh, i, j), &        ! SrvD%y%BStCLoadMesh(i, j), &
                             SrcDispDL=DatLoc(SrvD_u_BStCMotionMesh, i, j), &  ! SrvD%u%BStCMotionMesh(i, j)
                             DstDL=DatLoc(ED_u_BladePtLoads, i), &             ! ED%u%BladePtLoads(i)
                             DstDispDL=DatLoc(ED_y_BladeLn2Mesh, i), &         ! ED%y%BladeLn2Mesh(i)
                             Active=Turbine%p_FAST%CompElast == Module_ED, &
                             ErrStat=ErrStat2, ErrMsg=ErrMsg2)
            if (Failed()) return
         end do
      end do

      ! Nacelle Structural Controller
      do j = 1, Turbine%SrvD%p%NumNStC
         call MapLoadMesh(Turbine, Mappings, &
                          SrcMod=SrcMod, &
                          SrcDL=DatLoc(SrvD_y_NStCLoadMesh, j), &         ! SrvD%y%NStCLoadMesh(j), &
                          SrcDispDL=DatLoc(SrvD_u_NStCMotionMesh, j), &   ! SrvD%u%NStCMotionMesh(j)
                          DstMod=DstMod, &
                          DstDL=DatLoc(ED_u_NacelleLoads), &              ! ED%u%NacelleLoads
                          DstDispDL=DatLoc(ED_y_NacelleMotion), &         ! ED%y%NacelleMotion
                          ErrStat=ErrStat2, ErrMsg=ErrMsg2)
         if (Failed()) return
      end do

      ! Tower Structural Controller
      do j = 1, Turbine%SrvD%p%NumTStC
         call MapLoadMesh(Turbine, Mappings, &
                          SrcMod=SrcMod, &
                          SrcDL=DatLoc(SrvD_y_TStCLoadMesh, j), &         ! SrvD%y%TStCLoadMesh(j), &
                          SrcDispDL=DatLoc(SrvD_u_TStCMotionMesh, j), &   ! SrvD%u%TStCMotionMesh(j)
                          DstMod=DstMod, &
                          DstDL=DatLoc(ED_u_TowerPtLoads), &              ! ED%u%TowerLoads
                          DstDispDL=DatLoc(ED_y_TowerLn2Mesh), &          ! ED%y%TowerLn2Mesh
                          ErrStat=ErrStat2, ErrMsg=ErrMsg2)
         if (Failed()) return
      end do

      ! Substructure Structural Controller
      do j = 1, Turbine%SrvD%p%NumSStC
         call MapLoadMesh(Turbine, Mappings, &
                          SrcMod=SrcMod, &
                          SrcDL=DatLoc(SrvD_y_SStCLoadMesh, j), &         ! SrvD%y%SStCLoadMesh(j), &
                          SrcDispDL=DatLoc(SrvD_u_SStCMotionMesh, j), &   ! SrvD%u%SStCMotionMesh(j)
                          DstMod=DstMod, &
                          DstDL=DatLoc(ED_u_PlatformPtMesh), &            ! ED%u%PlatformPtMesh
                          DstDispDL=DatLoc(ED_y_PlatformPtMesh), &        ! ED%y%PlatformPtMesh
                          Active=Turbine%p_FAST%CompSub /= Module_SD, &
                          ErrStat=ErrStat2, ErrMsg=ErrMsg2)
         if (Failed()) return
      end do

   end select

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine InitMappings_SED(Mappings, SrcMod, DstMod, Turbine, ErrStat, ErrMsg)
   type(MappingType), allocatable      :: Mappings(:)
   type(ModDataType), intent(inout)       :: SrcMod, DstMod
   type(FAST_TurbineType), intent(inout)  :: Turbine           !< Turbine type
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg

   character(*), parameter    :: RoutineName = 'InitMappings_SED'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2
   integer(IntKi)             :: i, j

   ErrStat = ErrID_None
   ErrMsg = ''

   ! Select based on source module identifier
   select case (SrcMod%ID)

   case (Module_AD)

      ! Blade Loads
      do i = 1, size(Turbine%AD%y%rotors(SrcMod%Ins)%BladeLoad)
         call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                          SrcDL=DatLoc(AD_y_BladeLoad, i), &            ! AD%y%rotors(SrcMod%Ins)%BladeLoad(i)
                          SrcDispDL=DatLoc(AD_u_BladeMotion, i), &      ! AD%u%rotors(SrcMod%Ins)%BladeMotion(i)
                          DstDL=DatLoc(SED_u_HubPtLoad), &              ! SED%u%HubPtLoad
                          DstDispDL=DatLoc(SED_y_HubPtMotion), &        ! SED%y%HubPtMotion
                          ErrStat=ErrStat2, ErrMsg=ErrMsg2)
         if (Failed()) return
      end do

   case (Module_ADsk)

      ! Hub Loads
      call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                       SrcDL=DatLoc(ADsk_y_AeroLoads), &                ! ADsk%y%AeroLoads
                       SrcDispDL=DatLoc(ADsk_u_HubMotion), &            ! ADsk%u%HubMotion
                       DstDL=DatLoc(ED_u_HubPtLoad), &                  ! ED%u%HubPtLoad
                       DstDispDL=DatLoc(ED_y_HubPtMotion), &            ! ED%y%HubPtMotion
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2)
      if (Failed()) return

   case (Module_SrvD)

      call MapCustom(Mappings, Custom_SrvD_to_SED, SrcMod, DstMod)

   end select

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine InitMappings_ExtInfw(Mappings, SrcMod, DstMod, Turbine, ErrStat, ErrMsg)
   type(MappingType), allocatable      :: Mappings(:)
   type(ModDataType), intent(inout)       :: SrcMod, DstMod
   type(FAST_TurbineType), intent(inout)  :: Turbine           !< Turbine type
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg

   character(*), parameter    :: RoutineName = 'InitMappings_ExtInfw'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2
   integer(IntKi)             :: i

   ErrStat = ErrID_None
   ErrMsg = ''

   select case (SrcMod%ID)
   end select

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine InitMappings_ExtLd(Mappings, SrcMod, DstMod, Turbine, ErrStat, ErrMsg)
   type(MappingType), allocatable      :: Mappings(:)
   type(ModDataType), intent(inout)       :: SrcMod, DstMod
   type(FAST_TurbineType), intent(inout)  :: Turbine           !< Turbine type
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg

   character(*), parameter    :: RoutineName = 'InitMappings_ExtLd'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2
   integer(IntKi)             :: i, k
   logical                    :: CompElastED

   ErrStat = ErrID_None
   ErrMsg = ''

   ! Flag is true if CompElast == Module_ED
   CompElastED = Turbine%p_FAST%CompElast == Module_ED

   select case (SrcMod%ID)
   case (Module_AD)

      ! Blade Loads
      do i = 1, Turbine%ED%p%NumBl
         call MapLoadMesh(Turbine, Mappings, &
                          SrcMod=SrcMod, &
                          SrcDL=DatLoc(AD_y_BladeLoad, i), &            ! AD%y%rotors(SrcMod%Ins)%BladeLoad(i)
                          SrcDispDL=DatLoc(AD_u_BladeMotion, i), &      ! AD%u%rotors(SrcMod%Ins)%BladeMotion(i)
                          DstMod=DstMod, &
                          DstDL=DatLoc(ExtLd_u_BladeLoadAD, i), &       ! ExtLd%u%BladeLoadAD(i)
                          DstDispDL=DatLoc(ExtLd_u_BladeMotion, i), &   ! ExtLd%u%BladeMotion(i)
                          ErrStat=ErrStat2, ErrMsg=ErrMsg2)
         if (Failed()) return
      end do

      ! Tower Loads
      call MapLoadMesh(Turbine, Mappings, &
                       SrcMod=SrcMod, &
                       SrcDL=DatLoc(AD_y_TowerLoad), &               ! AD%y%rotors(SrcMod%Ins)%TowerLoad
                       SrcDispDL=DatLoc(AD_u_TowerMotion), &         ! AD%u%rotors(SrcMod%Ins)%TowerMotion
                       DstMod=DstMod, &
                       DstDL=DatLoc(ExtLd_u_TowerLoadAD), &          ! ExtLd%u%TowerLoadAD
                       DstDispDL=DatLoc(ExtLd_u_TowerMotion), &      ! ExtLd%u%TowerMotion
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2)
      if (Failed()) return
      

   case (Module_BD)

      ! Blade motion
      call MapMotionMesh(Turbine, Mappings, &
                         SrcMod=SrcMod, SrcDL=DatLoc(BD_y_BldMotion), &                   ! BD%y(SrcMod%Ins)%BldMotion
                         DstMod=DstMod, DstDL=DatLoc(ExtLd_u_BladeMotion, SrcMod%Ins), &  ! ExtLd%u%BladeMotion(SrcMod%Ins)
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2)
      if(Failed()) return

   case (Module_ED)

      call MapCustom(Mappings, Custom_ED_to_ExtLd, SrcMod, DstMod)

      ! Blade motion
      do i = 1, Turbine%ED%p%NumBl
         call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                            SrcDL=DatLoc(ED_y_BladeLn2Mesh, i), &        ! ED%y%BladeLn2Mesh(i)
                            DstDL=DatLoc(ExtLd_u_BladeMotion, i), &      ! ExtLd%u%BladeMotion(i)
                            Active=CompElastED, &
                            ErrStat=ErrStat2, ErrMsg=ErrMsg2)
         if(Failed()) return
      end do

      ! Blade root motion
      do i = 1, Turbine%ED%p%NumBl
         call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                            SrcDL=DatLoc(ED_y_BladeRootMotion, i), &     ! ED%y%BladeRootMotion(i)
                            DstDL=DatLoc(ExtLd_u_BladeRootMotion, i), &  ! ExtLd%u%BladeRootMotion(i)
                            ErrStat=ErrStat2, ErrMsg=ErrMsg2)
         if(Failed()) return
      end do

      ! Tower motion
      call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                         SrcDL=DatLoc(ED_y_TowerLn2Mesh), &              ! ED%y%TowerLn2Mesh
                         DstDL=DatLoc(ExtLd_u_TowerMotion), &            ! ExtLd%u%TowerMotion
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2)
      if(Failed()) return

      ! Hub motion
      call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                         SrcDL=DatLoc(ED_y_HubPtMotion), &               ! ED%y%HubPtMotion
                         DstDL=DatLoc(ExtLd_u_HubMotion), &              ! ExtLd%u%HubMotion
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2)
      if(Failed()) return

      ! Nacelle motion
      call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                         SrcDL=DatLoc(ED_y_NacelleMotion), &             ! ED%y%NacelleMotion
                         DstDL=DatLoc(ExtLd_u_NacelleMotion), &          ! ExtLd%u%NacelleMotion
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2)
      if(Failed()) return

   end select

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine InitMappings_ExtPtfm(Mappings, SrcMod, DstMod, Turbine, ErrStat, ErrMsg)
   type(MappingType), allocatable      :: Mappings(:)
   type(ModDataType), intent(inout)       :: SrcMod, DstMod
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

      call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                         SrcDL=DatLoc(ED_y_PlatformPtMesh), &                    ! ED%y%PlatformPtMesh
                         DstDL=DatLoc(ExtPtfm_u_PtfmMesh), &                     ! ExtPtfm%u%PtfmMesh
                         Active=Turbine%p_FAST%CompSub /= Module_SD, &
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

   end select

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine InitMappings_FEAM(Mappings, SrcMod, DstMod, Turbine, ErrStat, ErrMsg)
   type(MappingType), allocatable      :: Mappings(:)
   type(ModDataType), intent(inout)       :: SrcMod, DstMod
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
                         SrcDL=DatLoc(ED_y_PlatformPtMesh), &             ! ED%y%PlatformPtMesh
                         DstDL=DatLoc(FEAM_u_PtFairleadDisplacement), &   ! FEAM%u%PtFairleadDisplacement
                         Active=Turbine%p_FAST%CompSub /= Module_SD, &
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

   case (Module_SD)

      call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                         SrcDL=DatLoc(SD_y_Y3Mesh), &                     ! SD%y%y3Mesh
                         DstDL=DatLoc(FEAM_u_PtFairleadDisplacement), &   ! FEAM%u%PtFairleadDisplacement
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

   end select

   select case (SrcMod%ID)
   case (Module_ED)

      if (Turbine%p_FAST%CompSub /= Module_SD) then
         ! CALL MeshMapCreate( SubstructureMotion, FEAM%u%PtFairleadDisplacement,  MeshMapData%Structure_2_Mooring, ErrStat2, ErrMsg2 )
      end if

   case (Module_SD)

      ! CALL MeshMapCreate( SubstructureMotion, FEAM%u%PtFairleadDisplacement,  MeshMapData%Structure_2_Mooring, ErrStat2, ErrMsg2 )

   end select

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine InitMappings_HD(Mappings, SrcMod, DstMod, Turbine, ErrStat, ErrMsg)
   type(MappingType), allocatable      :: Mappings(:)
   type(ModDataType), intent(inout)       :: SrcMod, DstMod
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
                         SrcDL=DatLoc(ED_y_PlatformPtMesh), &          ! ED%y%PlatformPtMesh
                         DstDL=DatLoc(HydroDyn_u_PRPMesh), &           ! HD%u%PRPMesh
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

      call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                         SrcDL=DatLoc(ED_y_PlatformPtMesh), &       ! ED%y%PlatformPtMesh
                         DstDL=DatLoc(HydroDyn_u_Morison_Mesh), &   ! HD%u%Morison%Mesh
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2, &
                         Active=Turbine%p_FAST%CompSub /= Module_SD); if(Failed()) return

      call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                         SrcDL=DatLoc(ED_y_PlatformPtMesh), &       ! ED%y%PlatformPtMesh
                         DstDL=DatLoc(HydroDyn_u_WAMITMesh), &      ! HD%u%WAMITMesh
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2, &
                         Active=Turbine%p_FAST%CompSub /= Module_SD); if(Failed()) return

   case (Module_SeaSt)

      call MapVariable(Mappings, &
                       SrcMod=SrcMod, SrcDL=DatLoc(SeaSt_y_WaveElev0), &
                       DstMod=DstMod, DstDL=DatLoc(HydroDyn_u_WaveElev0), &
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2, &
                       Active=Turbine%p_FAST%Linearize); if (Failed()) return

   case (Module_SD)

      call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                         SrcDL=DatLoc(SD_y_Y2Mesh), &                  ! SD%y%Y2Mesh
                         DstDL=DatLoc(HydroDyn_u_Morison_Mesh), &      ! HD%u%Morison%Mesh
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

      call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                         SrcDL=DatLoc(SD_y_Y2Mesh), &                  ! SD%y%Y2Mesh
                         DstDL=DatLoc(HydroDyn_u_WAMITMesh), &         ! HD%u%WAMITMesh
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return
   end select

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine InitMappings_IceD(Mappings, SrcMod, DstMod, Turbine, ErrStat, ErrMsg)
   type(MappingType), allocatable      :: Mappings(:)
   type(ModDataType), intent(inout)       :: SrcMod, DstMod
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
                         SrcDL=DatLoc(ED_y_PlatformPtMesh), &          ! ED%y%PlatformPtMesh
                         DstDL=DatLoc(IceD_u_PointMesh), &             ! IceD%u%PointMesh
                         Active=Turbine%p_FAST%CompSub /= Module_SD, &
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

   case (Module_SD)

      call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                         SrcDL=DatLoc(SD_y_Y3Mesh), &                  ! SD%y%y3Mesh
                         DstDL=DatLoc(IceD_u_PointMesh), &             ! IceD%u%PointMesh
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

   end select

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine InitMappings_IceF(Mappings, SrcMod, DstMod, Turbine, ErrStat, ErrMsg)
   type(MappingType), allocatable      :: Mappings(:)
   type(ModDataType), intent(inout)       :: SrcMod, DstMod
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
                         SrcDL=DatLoc(ED_y_PlatformPtMesh), &          ! ED%y%PlatformPtMesh
                         DstDL=DatLoc(IceFloe_u_iceMesh), &            ! IceFloe%u%iceMesh
                         Active=Turbine%p_FAST%CompSub /= Module_SD, &
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

   case (Module_SD)

      call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                         SrcDL=DatLoc(SD_y_Y3Mesh), &                  ! SD%y%y3Mesh
                         DstDL=DatLoc(IceFloe_u_iceMesh), &            ! IceFloe%u%iceMesh
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

   end select

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine InitMappings_IfW(Mappings, SrcMod, DstMod, T, ErrStat, ErrMsg)
   type(MappingType), allocatable      :: Mappings(:)
   type(ModDataType), intent(inout)       :: SrcMod, DstMod
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
      call MapCustom(Mappings, Custom_SrvD_to_IfW, SrcMod, DstMod)
   end select

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine InitMappings_MAP(Mappings, SrcMod, DstMod, Turbine, ErrStat, ErrMsg)
   type(MappingType), allocatable      :: Mappings(:)
   type(ModDataType), intent(inout)       :: SrcMod, DstMod
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
                         SrcDL=DatLoc(ED_y_PlatformPtMesh), &       ! ED%y%PlatformPtMesh
                         DstDL=DatLoc(MAP_u_PtFairDisplacement), &  ! MAPp%u%PtFairDisplacement
                         Active=Turbine%p_FAST%CompSub /= Module_SD, &
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

   case (Module_SD)

      call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                         SrcDL=DatLoc(SD_y_Y3Mesh), &                  ! SD%y%y3Mesh
                         DstDL=DatLoc(MAP_u_PtFairDisplacement), &     ! MAPp%u%PtFairDisplacement
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

   end select

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine InitMappings_MD(Mappings, SrcMod, DstMod, Turbine, ErrStat, ErrMsg)
   type(MappingType), allocatable         :: Mappings(:)
   type(ModDataType), intent(inout)       :: SrcMod, DstMod
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
                         SrcDL=DatLoc(ED_y_PlatformPtMesh), &          ! ED%y%PlatformPtMesh
                         DstDL=DatLoc(MD_u_CoupledKinematics, 1), &    ! MD%u%CoupledKinematics(1)
                         Active=Turbine%p_FAST%CompSub /= Module_SD, &
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

   case (Module_SD)

      call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                         SrcDL=DatLoc(SD_y_Y3Mesh), &                  ! SD%y%y3Mesh
                         DstDL=DatLoc(MD_u_CoupledKinematics, 1), &    ! MD%u%CoupledKinematics(1)
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
   type(MappingType), allocatable      :: Mappings(:)
   type(ModDataType), intent(inout)       :: SrcMod, DstMod
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
                         SrcDL=DatLoc(ED_y_PlatformPtMesh), &          ! ED%y%PlatformPtMesh
                         DstDL=DatLoc(Orca_u_PtfmMesh), &              ! Orca%u%PtfmMesh
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

   end select

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine InitMappings_SD(Mappings, SrcMod, DstMod, Turbine, ErrStat, ErrMsg)
   type(MappingType), allocatable      :: Mappings(:)
   type(ModDataType), intent(inout)       :: SrcMod, DstMod
   type(FAST_TurbineType), intent(inout)  :: Turbine           !< Turbine type
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg

   character(*), parameter    :: RoutineName = 'InitMappings_SD'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2
   integer(IntKi)             :: i, j

   ErrStat = ErrID_None
   ErrMsg = ''

   select case (SrcMod%ID)

   case (Module_ED)

      call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                         SrcDL=DatLoc(ED_y_PlatformPtMesh), &          ! ED%y%PlatformPtMesh
                         DstDL=DatLoc(SD_u_TPMesh), &                  ! SD%u%TPMesh
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

   case (Module_FEAM)

      call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                       SrcDL=DatLoc(FEAM_y_PtFairleadLoad), &             ! FEAM%y%PtFairleadLoad, &
                       SrcDispDL=DatLoc(FEAM_u_PtFairleadDisplacement), & ! FEAM%u%PtFairleadDisplacement
                       DstDL=DatLoc(SD_u_LMesh), &                        ! SD%u%LMesh
                       DstDispDL=DatLoc(SD_y_y3Mesh), &                   ! SD%y%y3Mesh
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

   case (Module_HD)

      call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                       SrcDL=DatLoc(HydroDyn_y_Morison_Mesh), &          ! HD%y%Morison%Mesh
                       SrcDispDL=DatLoc(HydroDyn_u_Morison_Mesh), &      ! HD%u%Morison%Mesh
                       DstDL=DatLoc(SD_u_LMesh), &                       ! SD%u%LMesh
                       DstDispDL=DatLoc(SD_y_y2Mesh), &                  ! SD%y%y2Mesh
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

      call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                       SrcDL=DatLoc(HydroDyn_y_WAMITMesh), &             ! HD%y%WAMITMesh
                       SrcDispDL=DatLoc(HydroDyn_u_WAMITMesh), &         ! HD%u%WAMITMesh
                       DstDL=DatLoc(SD_u_LMesh), &                       ! SD%u%LMesh
                       DstDispDL=DatLoc(SD_y_y2Mesh), &                  ! SD%y%y2Mesh
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

   case (Module_IceD)

      call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                       SrcDL=DatLoc(IceD_y_PointMesh), &                 ! IceD%y%PointMesh
                       SrcDispDL=DatLoc(IceD_u_PointMesh), &             ! IceD%u%PointMesh
                       DstDL=DatLoc(SD_u_LMesh), &                       ! SD%u%LMesh
                       DstDispDL=DatLoc(SD_y_y3Mesh), &                  ! SD%y%y3Mesh
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

   case (Module_IceF)

      call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                       SrcDL=DatLoc(IceFloe_y_iceMesh), &                ! IceFloe%y%iceMesh
                       SrcDispDL=DatLoc(IceFloe_u_iceMesh), &            ! IceFloe%u%iceMesh
                       DstDL=DatLoc(SD_u_LMesh), &                       ! SD%u%LMesh
                       DstDispDL=DatLoc(SD_y_y3Mesh), &                  ! SD%y%y3Mesh
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

   case (Module_MAP)

      call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                       SrcDL=DatLoc(MAP_y_ptFairleadLoad), &             ! MAP%y%PtFairleadLoad
                       SrcDispDL=DatLoc(MAP_u_PtFairDisplacement), &     ! MAP%u%PtFairDisplacement
                       DstDL=DatLoc(SD_u_LMesh), &                       ! SD%u%LMesh
                       DstDispDL=DatLoc(SD_y_y3Mesh), &                  ! SD%y%y3Mesh
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

   case (Module_MD)

      call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                       SrcDL=DatLoc(MD_y_CoupledLoads, 1), &             ! MD%y%CoupledLoads(1)
                       SrcDispDL=DatLoc(MD_u_CoupledKinematics, 1), &    ! MD%u%CoupledKinematics(1)
                       DstDL=DatLoc(SD_u_LMesh), &                       ! SD%u%LMesh
                       DstDispDL=DatLoc(SD_y_y3Mesh), &                  ! SD%y%y3Mesh
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return

   case (Module_SrvD)

      call MapCustom(Mappings, Custom_SrvD_to_SD, SrcMod, DstMod)

      ! Substructure Structural Controller
      do j = 1, Turbine%SrvD%p%NumSStC
         call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                          SrcDL=DatLoc(SrvD_y_SStCLoadMesh, j), &         ! SrvD%y%SStCLoadMesh(j), &
                          SrcDispDL=DatLoc(SrvD_u_SStCMotionMesh, j), &   ! SrvD%u%SStCMotionMesh(j)
                          DstDL=DatLoc(SD_u_LMesh), &                     ! SD%u%LMesh
                          DstDispDL=DatLoc(SD_y_y3Mesh), &                ! SD%y%y3Mesh
                          ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return
      end do

   end select

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine InitMappings_SeaSt(Mappings, SrcMod, DstMod, Turbine, ErrStat, ErrMsg)
   type(MappingType), allocatable      :: Mappings(:)
   type(ModDataType), intent(inout)       :: SrcMod, DstMod
   type(FAST_TurbineType), intent(inout)  :: Turbine           !< Turbine type
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg

   character(*), parameter    :: RoutineName = 'InitMappings_SeaSt'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2
   integer(IntKi)             :: i, j

   ErrStat = ErrID_None
   ErrMsg = ''

   select case (SrcMod%ID)

      ! No inputs to SeaState currently

   end select

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine InitMappings_SrvD(Mappings, SrcMod, DstMod, Turbine, ErrStat, ErrMsg)
   type(MappingType), allocatable      :: Mappings(:)
   type(ModDataType), intent(inout)       :: SrcMod, DstMod
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

   case (Module_BD)

      call MapCustom(Mappings, Custom_BD_to_SrvD, SrcMod, DstMod)

      ! Blade Structural Controller
      do i = 1, Turbine%SrvD%p%NumBStC
         call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                            SrcDL=DatLoc(BD_y_BldMotion), &                         ! BD%y%BldMotion
                            DstDL=DatLoc(SrvD_u_BStCMotionMesh, DstMod%Ins, i), &   ! SrvD%u%BStCMotionMesh(i, j)
                            ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return
      end do

   case (Module_ED)

      call MapCustom(Mappings, Custom_ED_to_SrvD, SrcMod, DstMod)

      call MapVariable(Mappings, &
                       SrcMod=SrcMod, SrcDL=DatLoc(ED_y_Yaw), &
                       DstMod=DstMod, DstDL=DatLoc(SrvD_u_Yaw), &
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2); if (Failed()) return

      call MapVariable(Mappings, &
                       SrcMod=SrcMod, SrcDL=DatLoc(ED_y_YawRate), &
                       DstMod=DstMod, DstDL=DatLoc(SrvD_u_YawRate), &
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2); if (Failed()) return

      call MapVariable(Mappings, &
                       SrcMod=SrcMod, SrcDL=DatLoc(ED_y_HSS_Spd), &
                       DstMod=DstMod, DstDL=DatLoc(SrvD_u_HSS_Spd), &
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2); if (Failed()) return

      ! Nacelle Structural Controller
      do j = 1, Turbine%SrvD%p%NumNStC
         call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                            SrcDL=DatLoc(ED_y_NacelleMotion), &          ! ED%y%NacelleMotion
                            DstDL=DatLoc(SrvD_u_NStCMotionMesh, j), &    ! SrvD%u%NStCMotionMesh(j)
                            ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return
      end do

      ! Tower Structural Controller
      do j = 1, Turbine%SrvD%p%NumTStC
         call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                            SrcDL=DatLoc(ED_y_TowerLn2Mesh), &           ! ED%y%TowerMotion
                            DstDL=DatLoc(SrvD_u_TStCMotionMesh, j), &    ! SrvD%u%TStCMotionMesh(j)
                            ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return
      end do

      ! Blade Structural Controller (if ElastoDyn blades)
      do j = 1, Turbine%SrvD%p%NumBStC
         do i = 1, Turbine%ED%p%NumBl
            call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                               SrcDL=DatLoc(ED_y_BladeLn2Mesh, i), &         ! ED%y%BladeLn2Mesh(i)
                               DstDL=DatLoc(SrvD_u_BStCMotionMesh, i, j), &  ! SrvD%u%BStCMotionMesh(i, j)
                               Active=Turbine%p_FAST%CompElast == Module_ED, &
                               ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return
         end do
      end do

      ! Substructure Structural Controller (if not using SubDyn)
      do j = 1, Turbine%SrvD%p%NumSStC
         call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                            SrcDL=DatLoc(ED_y_PlatformPtMesh), &         ! ED%y%PlatformPtMesh
                            DstDL=DatLoc(SrvD_u_SStCMotionMesh, j), &    ! SrvD%u%SStCMotionMesh(j)
                            Active=Turbine%p_FAST%CompSub /= Module_SD, &
                            ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return
      end do

   case (Module_IfW)

      call MapCustom(Mappings, Custom_IfW_to_SrvD, SrcMod=SrcMod, DstMod=DstMod)

   case (Module_SD)

      ! Substructure Structural Controller
      do j = 1, Turbine%SrvD%p%NumSStC
         call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                            SrcDL=DatLoc(SD_y_y3Mesh), &               ! SD%y%y3Mesh
                            DstDL=DatLoc(SrvD_u_SStCMotionMesh, j), &  ! SrvD%u%SStCMotionMesh(j)
                            ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return
      end do

   end select

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine MapLoadMesh(Turbine, Mappings, SrcMod, SrcDL, SrcDispDL, &
                       DstMod, DstDL, DstDispDL, ErrStat, ErrMsg, Active)
   type(FAST_TurbineType), target         :: Turbine
   type(MappingType), allocatable         :: Mappings(:)
   type(ModDataType), intent(inout)       :: SrcMod, DstMod
   type(DatLoc), intent(in)               :: SrcDL, DstDL
   type(DatLoc), intent(in)               :: SrcDispDL, DstDispDL
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg
   logical, optional, intent(in)          :: Active

   character(*), parameter                :: RoutineName = 'MapLoadMesh'
   integer(IntKi)                         :: ErrStat2
   character(ErrMsgLen)                   :: ErrMsg2
   type(MappingType)                      :: Mapping
   type(MeshType), pointer                :: SrcMesh, SrcDispMesh
   type(MeshType), pointer                :: DstMesh, DstDispMesh
   type(MeshType)                         :: DstMotionMesh

   ErrStat = ErrID_None
   ErrMsg = ''

   ! If active argument is set to false, return
   if (present(Active)) then
      if (.not. Active) return
   end if

   ! Get mesh pointers (DstDispMesh may be found in Input for some modules: ExtLd)
   call FAST_OutputMeshPointer(SrcMod, Turbine, SrcDL, SrcMesh, ErrStat2, ErrMsg2); if (Failed()) return
   call FAST_InputMeshPointer(SrcMod, Turbine, SrcDispDL, SrcDispMesh, INPUT_CURR, ErrStat2, ErrMsg2); if (Failed()) return
   call FAST_InputMeshPointer(DstMod, Turbine, DstDL, DstMesh, INPUT_CURR, ErrStat2, ErrMsg2); if (Failed()) return
   call FAST_OutputMeshPointer(DstMod, Turbine, DstDispDL, DstDispMesh, ErrStat2, ErrMsg2)
   if (ErrStat2 == ErrID_Fatal) call FAST_InputMeshPointer(DstMod, Turbine, DstDispDL, DstDispMesh, INPUT_CURR, ErrStat2, ErrMsg2)
   if (Failed()) return

   ! If any meshes aren't committed, return
   if (.not. (SrcMesh%committed .and. DstMesh%committed .and. SrcDispMesh%committed .and. DstDispMesh%committed)) return

   ! Check that all meshes in mapping have nonzero identifiers
   if (SrcMesh%ID == 0) then
      call SetErrStat(ErrID_Fatal, 'SrcMesh "'//trim(FAST_OutputFieldName(SrcMod, SrcDL))//'" not in module variables', &
                      ErrStat, ErrMsg, RoutineName)
      return
   else if (SrcDispMesh%ID == 0) then
      call SetErrStat(ErrID_Fatal, 'SrcDispMesh "'//trim(FAST_InputFieldName(SrcMod, SrcDispDL))//'" not in module variables', &
                      ErrStat, ErrMsg, RoutineName)
      return
   else if (DstMesh%ID == 0) then
      call SetErrStat(ErrID_Fatal, 'DstMesh "'//trim(FAST_InputFieldName(DstMod, DstDL))//'" not in module variables', &
                      ErrStat, ErrMsg, RoutineName)
      return
   else if (DstDispMesh%ID == 0) then
      call SetErrStat(ErrID_Fatal, 'DstDispMesh "'//trim(FAST_OutputFieldName(DstMod, DstDispDL))//'" not in module variables', &
                      ErrStat, ErrMsg, RoutineName)
      return
   end if

   ! Create mapping description
   Mapping%Desc = trim(FAST_OutputFieldName(SrcMod, SrcDL))//" -> "// &
                  trim(FAST_InputFieldName(DstMod, DstDL))// &
                  " ["//trim(FAST_InputFieldName(SrcMod, SrcDispDL))// &
                  " @ "//trim(FAST_OutputFieldName(DstMod, DstDispDL))//"]"

   ! Initialize mapping structure
   Mapping%MapType = Map_LoadMesh
   Mapping%iModSrc = SrcMod%iMod
   Mapping%SrcModID = SrcMod%ID
   Mapping%SrcIns = SrcMod%Ins
   Mapping%iModDst = DstMod%iMod
   Mapping%DstModID = DstMod%ID
   Mapping%DstIns = DstMod%Ins
   Mapping%SrcDL = SrcDL
   Mapping%SrcDispDL = SrcDispDL
   Mapping%DstDL = DstDL
   Mapping%DstDispDL = DstDispDL
   Mapping%XfrType = MeshTransferType(SrcMesh, DstMesh)

   ! Create mesh mapping
   call MeshMapCreate(SrcMesh, DstMesh, Mapping%MeshMap, ErrStat2, ErrMsg2); if (Failed()) return

   ! Create a copy of destination mesh in mapping for load summation
   call MeshCopy(DstMesh, Mapping%TmpLoadMesh, MESH_NEWCOPY, ErrStat2, ErrMsg2); if (Failed()) return

   ! Set VF_Mapping on variables in this mapping
   call SetMapVarFlags(Mapping, SrcMod, DstMod)

   ! If the destination displacement mesh is not a sibling of the load mesh
   Mapping%DstUsesSibling = IsSiblingMesh(DstMesh, DstDispMesh)
   if (.not. Mapping%DstUsesSibling) then

      ! Indicate non-sibling destination displacement mesh in description
      Mapping%Desc = trim(Mapping%Desc)//'*'

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
   pure logical function IsSiblingMesh(MeshA, MeshB)
      type(MeshType), intent(in) :: MeshA, MeshB
      integer(IntKi)             :: i, j
      IsSiblingMesh = .false.
      if (MeshA%Nnodes /= MeshB%Nnodes) return
      if (any(MeshA%Position /= MeshB%Position)) return
      if (any(MeshA%RefOrientation /= MeshB%RefOrientation)) return
      do i = 1, NELEMKINDS
         if (MeshA%ElemTable(i)%nelem /= MeshB%ElemTable(i)%nelem) return
         if (MeshA%ElemTable(i)%XElement /= MeshB%ElemTable(i)%XElement) return
         do j = 1, MeshA%ElemTable(i)%nelem
            if (any(MeshA%ElemTable(i)%Elements(j)%ElemNodes /= MeshB%ElemTable(i)%Elements(j)%ElemNodes)) return
         end do
      end do
      IsSiblingMesh = .true.
   end function
end subroutine

subroutine MapMotionMesh(Turbine, Mappings, SrcMod, SrcDL, DstMod, DstDL, ErrStat, ErrMsg, Active)
   type(FAST_TurbineType), target         :: Turbine
   type(MappingType), allocatable         :: Mappings(:)
   type(ModDataType), intent(inout)       :: SrcMod, DstMod
   type(DatLoc), intent(in)               :: SrcDL, DstDL
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg
   logical, optional, intent(in)          :: Active

   character(*), parameter                :: RoutineName = 'MapMotionMesh'
   integer(IntKi)                         :: ErrStat2
   character(ErrMsgLen)                   :: ErrMsg2
   type(MappingType)                      :: Mapping
   type(MeshType), pointer                :: SrcMesh, DstMesh

   ErrStat = ErrID_None
   ErrMsg = ''

   ! If active argument is set to false, return
   if (present(Active)) then
      if (.not. Active) return
   end if

   ! Get mesh pointers
   call FAST_OutputMeshPointer(SrcMod, Turbine, SrcDL, SrcMesh, ErrStat2, ErrMsg2); if (Failed()) return
   call FAST_InputMeshPointer(DstMod, Turbine, DstDL, DstMesh, INPUT_CURR, ErrStat2, ErrMsg2); if (Failed()) return

   ! If source or destination meshes aren't commited, return
   if (.not. (SrcMesh%committed .and. DstMesh%committed)) return

   ! Check that all meshes in mapping have nonzero identifiers
   if (SrcMesh%ID == 0) then
      call SetErrStat(ErrID_Fatal, 'SrcMesh "'//trim(FAST_OutputFieldName(SrcMod, SrcDL))//'" not in module variables', &
                      ErrStat, ErrMsg, RoutineName)
      return
   else if (DstMesh%ID == 0) then
      call SetErrStat(ErrID_Fatal, 'DstMesh "'//trim(FAST_InputFieldName(DstMod, DstDL))//'" not in module variables', &
                      ErrStat, ErrMsg, RoutineName)
      return
   end if

   ! Create mapping description
   Mapping%Desc = trim(FAST_OutputFieldName(SrcMod, SrcDL))//" -> "// &
                  trim(FAST_InputFieldName(DstMod, DstDL))

   ! Initialize mapping structure
   Mapping%MapType = Map_MotionMesh
   Mapping%iModSrc = SrcMod%iMod
   Mapping%SrcModID = SrcMod%ID
   Mapping%SrcIns = SrcMod%Ins
   Mapping%iModDst = DstMod%iMod
   Mapping%DstModID = DstMod%ID
   Mapping%DstIns = DstMod%Ins
   Mapping%SrcDL = SrcDL
   Mapping%DstDL = DstDL
   Mapping%XfrType = MeshTransferType(SrcMesh, DstMesh)

   ! Set VF_Mapping on variables in this mapping
   call SetMapVarFlags(Mapping, SrcMod, DstMod)

   ! Create mesh mapping
   call MeshMapCreate(SrcMesh, DstMesh, Mapping%MeshMap, ErrStat2, ErrMsg2); if (Failed()) return

   ! Add mapping to array of mappings
   Mappings = [Mappings, Mapping]

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine MapVariable(Maps, SrcMod, SrcDL, DstMod, DstDL, ErrStat, ErrMsg, Active)
   type(MappingType), allocatable      :: Maps(:)
   type(ModDataType), intent(inout)    :: SrcMod, DstMod
   type(DatLoc), intent(in)            :: SrcDL, DstDL
   integer(IntKi), intent(out)         :: ErrStat
   character(*), intent(out)           :: ErrMsg
   logical, optional, intent(in)       :: Active
   type(MappingType)                   :: Mapping
   integer(IntKi)                      :: iVarSrc, iVarDst

   ErrStat = ErrID_None
   ErrMsg = ''

   if (present(Active)) then
      if (.not. Active) return
   end if

   ! Get source and destination variable indices
   iVarSrc = MV_FindVarDatLoc(SrcMod%Vars%y, SrcDL)
   iVarDst = MV_FindVarDatLoc(DstMod%Vars%u, DstDL)

   ! If either variable index is zero, return error
   if (iVarSrc == 0) then
      ErrStat = ErrID_Fatal
      ErrMsg = "Source variable "//trim(Num2LStr(SrcDL%Num))//" in module '"//trim(SrcMod%Abbr)//"' is not active"
      return
   else if (iVarDst == 0) then
      ErrStat = ErrID_Fatal
      ErrMsg = "Destination variable "//trim(Num2LStr(DstDL%Num))//" in module '"//trim(DstMod%Abbr)//"' is not active"
      return
   end if

   ! Create mapping description
   Mapping%Desc = trim(FAST_OutputFieldName(SrcMod, SrcDL))//" -> "// &
                  trim(FAST_InputFieldName(DstMod, DstDL))

   ! Verify that variables have compatible sizes
   ! If source variable has size 1, it can be mapped to multiple destination variables
   if ((SrcMod%Vars%y(iVarSrc)%Num > 1) .and. &
       (SrcMod%Vars%y(iVarSrc)%Num /= DstMod%Vars%u(iVarDst)%Num)) then
      ErrStat = ErrID_Fatal
      ErrMsg = "Variables in mapping '"//trim(Mapping%Desc)//"' have incompatible sizes"
      return
   end if

   ! Initialize mapping structure
   Mapping%MapType = Map_Variable
   Mapping%iModSrc = SrcMod%iMod
   Mapping%iModDst = DstMod%iMod
   Mapping%SrcModID = SrcMod%ID
   Mapping%DstModID = DstMod%ID
   Mapping%SrcIns = SrcMod%Ins
   Mapping%DstIns = DstMod%Ins
   Mapping%SrcDL = SrcDL
   Mapping%DstDL = DstDL

   ! Set VF_Mapping on variables in this mapping
   call SetMapVarFlags(Mapping, SrcMod, DstMod)

   ! Copy source and destination variables and modify for packing/unpacking
   Mapping%SrcVar = SrcMod%Vars%y(iVarSrc)
   Mapping%DstVar = DstMod%Vars%u(iVarDst)
   Mapping%SrcVar%iLoc = [1, Mapping%SrcVar%Num]
   Mapping%DstVar%iLoc = [1, Mapping%DstVar%Num]

   ! Allocate variable data storage
   call AllocAry(Mapping%VarData, max(Mapping%SrcVar%Num, Mapping%DstVar%Num), "VarData", ErrStat, ErrMsg)

   Maps = [Maps, Mapping]
end subroutine

!> MapCustom creates a custom mapping that is not included in linearization.
!! Each custom mapping needs an entry in FAST_InputSolve to actually perform the transfer.
subroutine MapCustom(Maps, Desc, SrcMod, DstMod, Active)
   type(MappingType), allocatable   :: Maps(:)
   character(*), intent(in)         :: Desc
   type(ModDataType), intent(inout) :: SrcMod, DstMod
   logical, optional, intent(in)    :: Active
   type(MappingType)                :: Mapping

   if (present(Active)) then
      if (.not. Active) return
   end if

   ! Initialize mapping structure
   Mapping%Desc = Desc
   Mapping%MapType = Map_Custom
   Mapping%iModSrc = SrcMod%iMod
   Mapping%iModDst = DstMod%iMod
   Mapping%SrcModID = SrcMod%ID
   Mapping%DstModID = DstMod%ID
   Mapping%SrcIns = SrcMod%Ins
   Mapping%DstIns = DstMod%Ins

   Maps = [Maps, Mapping]
end subroutine

subroutine SetMapVarFlags(Mapping, SrcMod, DstMod)
   type(MappingType), intent(in)          :: Mapping
   type(ModDataType), intent(inout)       :: SrcMod, DstMod
   integer(IntKi)                         :: i

   ! Set mapping flag on source variables
   do i = 1, size(SrcMod%Vars%y)
      associate (Var => SrcMod%Vars%y(i))
         if (MV_EqualDL(Mapping%SrcDL, Var%DL)) call MV_SetFlags(Var, VF_Mapping)
      end associate
   end do

   ! Set mapping flag on destination variables
   do i = 1, size(DstMod%Vars%u)
      associate (Var => DstMod%Vars%u(i))
         if (MV_EqualDL(Mapping%DstDL, Var%DL)) call MV_SetFlags(Var, VF_Mapping)
      end associate
   end do

   ! If this a load mesh mapping
   if (Mapping%MapType == Map_LoadMesh) then

      ! Set mapping flag on source displacement mesh variables
      do i = 1, size(SrcMod%Vars%u)
         associate (Var => SrcMod%Vars%u(i))
            if (MV_EqualDL(Mapping%SrcDispDL, Var%DL)) then
               select case (Var%Field)
               case (FieldTransDisp)
                  call MV_SetFlags(Var, VF_Mapping)
               end select
            end if
         end associate
      end do

      ! Set mapping flag on destination displacement mesh variables
      do i = 1, size(DstMod%Vars%y)
         associate (Var => DstMod%Vars%y(i))
            if (MV_EqualDL(Mapping%DstDispDL, Var%DL)) then
               select case (Var%Field)
               case (FieldTransDisp, FieldOrientation)
                  call MV_SetFlags(Var, VF_Mapping)
               end select
            end if
         end associate
      end do
   end if

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

subroutine FAST_LinearizeMappings(ModGlue, Mappings, Turbine, ErrStat, ErrMsg)
   type(ModGlueType), intent(inout)                :: ModGlue     !< Glue module data
   type(MappingType), intent(inout)                :: Mappings(:) !< Variable mappings
   type(FAST_TurbineType), target, intent(inout)   :: Turbine     !< Turbine type
   integer(IntKi), intent(out)                     :: ErrStat
   character(*), intent(out)                       :: ErrMsg

   character(*), parameter       :: RoutineName = 'FAST_LinearizeMappings'
   integer(IntKi)                :: ErrStat2
   character(ErrMsgLen)          :: ErrMsg2
   integer(IntKi)                :: iGluSrc(2), iGluDst(2), nLocSrc, nLocDst
   integer(IntKi)                :: i, j, k
   type(MeshType), pointer       :: SrcMesh, DstMesh
   type(MeshType), pointer       :: SrcDispMesh, DstDispMesh

   ErrStat = ErrID_None
   ErrMsg = ''

   ! Exit function if dUdy and dUdu aren't allocated
   if (.not. (allocated(ModGlue%Lin%dUdy) .and. allocated(ModGlue%Lin%dUdu))) return

   ! Initialize dUdy to zero
   ModGlue%Lin%dUdy = 0.0_R8Ki

   ! Initialize dUdu to identity matrix
   call Eye2D(ModGlue%Lin%dUdu, ErrStat2, ErrMsg2); if (Failed()) return

   ! Loop through variable maps
   do i = 1, size(ModGlue%VarMaps)

      associate (ModMap => ModGlue%VarMaps(i), &
                 Mapping => Mappings(ModGlue%VarMaps(i)%iMapping), &
                 ModSrc => ModGlue%ModData(ModGlue%VarMaps(i)%iModSrc), &
                 ModDst => ModGlue%ModData(ModGlue%VarMaps(i)%iModDst))

         ! Select based on type of mapping
         select case (Mapping%MapType)

         case (Map_Variable)

            ! Get source and destination indices, skip if no variable index for either
            if (ModMap%iVarSrc(1) == 0 .or. ModMap%iVarDst(1) == 0) cycle
            iGluSrc = ModSrc%Vars%y(ModMap%iVarSrc(1))%iGlu
            iGluDst = ModDst%Vars%u(ModMap%iVarDst(1))%iGlu

            ! Get number of source and destination locations
            nLocSrc = iGluSrc(2) - iGluSrc(1) + 1
            nLocDst = iGluDst(2) - iGluDst(1) + 1

            ! If source has multiple locations, destination must have same number, connect 1-to-1
            ! MapVariable checks that variables have same number if nLocSrc > 1
            if (nLocSrc > 1) then
               do k = 0, nLocDst - 1
                  ModGlue%Lin%dUdy(iGluDst(1) + k, iGluSrc(1) + k) = -1.0_R8Ki
               end do
            else if (nLocDst == 1) then
               ! Source and destination have one location
               ModGlue%Lin%dUdy(iGluDst(1), iGluSrc(1)) = -1.0_R8Ki
            else
               ! One source location to many destination locations
               ModGlue%Lin%dUdy(iGluDst(1):iGluDst(2), iGluSrc(1)) = -1.0_R8Ki
            end if

         case (Map_MotionMesh)

            ! Get source and destination meshes
            call FAST_OutputMeshPointer(ModSrc, Turbine, Mapping%SrcDL, SrcMesh, ErrStat2, ErrMsg2); if (Failed()) return
            call FAST_InputMeshPointer(ModDst, Turbine, Mapping%DstDL, DstMesh, INPUT_CURR, ErrStat2, ErrMsg2); if (Failed()) return

            ! Perform linearization based on transfer type
            call LinearizeMeshTransfer(Mapping%XfrType, SrcMesh, DstMesh, Mapping%MeshMap); if (Failed()) return

            ! Copy linearization matrices to global dUdy matrix
            call Assemble_dUdy_Motions(Mapping, ModMap, ModSrc%Vars, ModDst%Vars, ModGlue%Lin%dUdy)

            ! Copy linearization matrices to global dUdu matrix
            call Assemble_dUdu(Mapping, ModMap, ModSrc%Vars, ModDst%Vars, ModGlue%Lin%dUdu)

         case (Map_LoadMesh)

            ! Get source and destination meshes
            call FAST_OutputMeshPointer(ModSrc, Turbine, Mapping%SrcDL, SrcMesh, ErrStat2, ErrMsg2); if (Failed()) return
            call FAST_InputMeshPointer(ModDst, Turbine, Mapping%DstDL, DstMesh, INPUT_CURR, ErrStat2, ErrMsg2); if (Failed()) return

            ! Get source and destination displacement meshes (DstDispMesh must be in output)
            call FAST_InputMeshPointer(ModSrc, Turbine, Mapping%SrcDispDL, SrcDispMesh, INPUT_CURR, ErrStat2, ErrMsg2); if (Failed()) return
            call FAST_OutputMeshPointer(ModDst, Turbine, Mapping%DstDispDL, DstDispMesh, ErrStat2, ErrMsg2); if (Failed()) return

            ! If DstDispMesh is a sibling of DstMesh
            if (Mapping%DstUsesSibling) then

               ! Linearize the load mesh transfer
               call LinearizeMeshTransfer(Mapping%XfrType, SrcMesh, DstMesh, Mapping%MeshMap, SrcDispMesh, DstDispMesh); if (Failed()) return

            else

               ! Transfer destination displacement mesh to temporary motion mesh (cousin of destination load mesh)
               call TransferMesh(Mapping%XfrTypeAux, DstDispMesh, Mapping%TmpMotionMesh, Mapping%MeshMapAux, ErrStat=ErrStat2, ErrMsg=ErrMsg2); if (Failed()) return

               ! Linearize the motion mesh transfer
               call LinearizeMeshTransfer(Mapping%XfrTypeAux, DstDispMesh, Mapping%TmpMotionMesh, Mapping%MeshMapAux); if (Failed()) return

               ! Linearize the load mesh transfer
               call LinearizeMeshTransfer(Mapping%XfrType, SrcMesh, DstMesh, Mapping%MeshMap, SrcDispMesh, Mapping%TmpMotionMesh); if (Failed()) return

            end if

            ! Copy linearization matrices to global dUdy matrix
            call Assemble_dUdy_Loads(Mapping, ModMap, ModSrc%Vars, ModDst%Vars, ModGlue%Lin%dUdy)

            ! Copy linearization matrices to global dUdu matrix
            call Assemble_dUdu(Mapping, ModMap, ModSrc%Vars, ModDst%Vars, ModGlue%Lin%dUdu)

         end select

      end associate

   end do

contains

   ! LinearizeMeshTransfer calls the specific linearization function based on
   ! transfer type (Point_to_Point, Point_to_Line2, etc.)
   subroutine LinearizeMeshTransfer(Typ, Src, Dst, MeshMap, SrcDisp, DstDisp)
      integer(IntKi), intent(in)             :: Typ
      type(MeshType), intent(in)             :: Src, Dst
      type(MeshMapType), intent(inout)       :: MeshMap
      type(MeshType), optional, intent(in)   :: SrcDisp, DstDisp
      select case (Typ)
      case (Xfr_Point_to_Point)
         call Linearize_Point_to_Point(Src, Dst, MeshMap, ErrStat2, ErrMsg2, SrcDisp, DstDisp)
      case (Xfr_Point_to_Line2)
         call Linearize_Point_to_Line2(Src, Dst, MeshMap, ErrStat2, ErrMsg2, SrcDisp, DstDisp)
      case (Xfr_Line2_to_Point)
         call Linearize_Line2_to_Point(Src, Dst, MeshMap, ErrStat2, ErrMsg2, SrcDisp, DstDisp)
      case (Xfr_Line2_to_Line2)
         call Linearize_Line2_to_Line2(Src, Dst, MeshMap, ErrStat2, ErrMsg2, SrcDisp, DstDisp)
      case default
         ErrStat2 = ErrID_Fatal
         ErrMsg2 = "LinearizeMeshTransfer: unknown transfer type: "//Num2LStr(Typ)
      end select
   end subroutine

   subroutine Assemble_dUdu(Mapping, ModMap, VarsSrc, VarsDst, dUdu)
      type(MappingType), intent(in) :: Mapping
      type(VarMapType), intent(in)  :: ModMap
      type(ModVarsType), intent(in) :: VarsSrc, VarsDst
      real(R8Ki), intent(inout)     :: dUdu(:, :)

      ! Effect of input Translation Displacement on input Translation Velocity
      if (allocated(Mapping%MeshMap%dM%tv_uD)) then
         call SumBlock(VarsDst%u, ModMap%iVarDst(FieldTransDisp), VarsDst%u, ModMap%iVarDst(FieldTransVel), Mapping%MeshMap%dM%tv_uD, dUdu)
      end if

      ! Effect of input Translation Displacement on input Translation Acceleration
      if (allocated(Mapping%MeshMap%dM%ta_uD)) then
         call SumBlock(VarsDst%u, ModMap%iVarDst(FieldTransDisp), VarsDst%u, ModMap%iVarDst(FieldTransAcc), Mapping%MeshMap%dM%ta_uD, dUdu)
      end if

      ! Effect of input Translation Displacement on input Moments
      if (allocated(Mapping%MeshMap%dM%M_uS)) then
         call SumBlock(VarsSrc%u, ModMap%iVarSrcDisp(FieldTransDisp), VarsDst%u, ModMap%iVarDst(FieldMoment), Mapping%MeshMap%dM%M_uS, dUdu)
      end if
   end subroutine

   !> Assemble_dUdy_Loads assembles the linearization matrices for transfer of
   !! load fields between two meshes. It sets the following block matrix, which
   !! is the dUdy block for transfering output (source) mesh  to the input
   !! (destination) mesh :
   !! M = -| M_li   0    | * M_mi | F^S |
   !!      | M_fm   M_li |        | M^S |
   subroutine Assemble_dUdy_Loads(Mapping, ModMap, VarsSrc, VarsDst, dUdy)
      type(MappingType), intent(inout) :: Mapping
      type(VarMapType), intent(in)     :: ModMap
      type(ModVarsType), intent(in)    :: VarsSrc, VarsDst
      real(R8Ki), intent(inout)        :: dUdy(:, :)

      ! Load identity
      if (allocated(Mapping%MeshMap%dM%li)) then
         call SumBlock(VarsSrc%y, ModMap%iVarSrc(FieldForce), VarsDst%u, ModMap%iVarDst(FieldForce), Mapping%MeshMap%dM%li, dUdy)
         call SumBlock(VarsSrc%y, ModMap%iVarSrc(FieldMoment), VarsDst%u, ModMap%iVarDst(FieldMoment), Mapping%MeshMap%dM%li, dUdy)
      end if

      ! Force to Moment
      if (allocated(Mapping%MeshMap%dM%m_f)) then
         call SumBlock(VarsSrc%y, ModMap%iVarSrc(FieldForce), VarsDst%u, ModMap%iVarDst(FieldMoment), Mapping%MeshMap%dM%m_f, dUdy)
      end if

      ! Destination Translation Displacement to Moment
      if (allocated(Mapping%MeshMap%dM%m_uD)) then
         if (Mapping%DstUsesSibling) then
            ! Direct transfer
            call SumBlock(VarsDst%y, ModMap%iVarDstDisp(FieldTransDisp), VarsDst%u, ModMap%iVarDst(FieldMoment), Mapping%MeshMap%dM%m_uD, dUdy)
         else
            ! Compose linearization of motion and loads
            Mapping%TmpMatrix = matmul(Mapping%MeshMap%dM%m_uD, Mapping%MeshMapAux%dM%mi)
            call SumBlock(VarsDst%y, ModMap%iVarDstDisp(FieldTransDisp), VarsDst%u, ModMap%iVarDst(FieldMoment), Mapping%TmpMatrix, dUdy)
            Mapping%TmpMatrix = matmul(Mapping%MeshMap%dM%m_uD, Mapping%MeshMapAux%dM%fx_p)
            call SumBlock(VarsDst%y, ModMap%iVarDstDisp(FieldOrientation), VarsDst%u, ModMap%iVarDst(FieldMoment), Mapping%TmpMatrix, dUdy)
         end if
      end if
   end subroutine

   !> Assemble_dUdy_Motions assembles the linearization matrices for transfer of
   !! motion fields between two meshes. It set the following block matrix, which
   !! is the dUdy block for transfering output (source) mesh  to the input
   !! (destination) mesh :
   !! M = -| M_mi      M_f_p   0      0         0      0     |
   !!      | 0         M_mi    0      0         0      0     |
   !!      | M_tv_uS   0       M_mi   M_f_p     0      0     |
   !!      | 0         0       0      M_mi      0      0     |
   !!      | M_ta_uS   0       0      M_ta_rv   M_mi   M_f_p |
   !!      | 0         0       0      0         0      M_mi  |
   !! where the matrices correspond to
   !! u^S, theta^S, v^S, omega^S, a^S, alpha^S
   subroutine Assemble_dUdy_Motions(Mapping, ModMap, VarsSrc, VarsDst, dUdy)
      type(MappingType), intent(in) :: Mapping
      type(VarMapType), intent(in)  :: ModMap
      type(ModVarsType), intent(in) :: VarsSrc, VarsDst
      real(R8Ki), intent(inout)     :: dUdy(:, :)

      ! Motion identity
      if (allocated(Mapping%MeshMap%dM%mi)) then
         call SumBlock(VarsSrc%y, ModMap%iVarSrc(FieldTransDisp), VarsDst%u, ModMap%iVarDst(FieldTransDisp), Mapping%MeshMap%dM%mi, dUdy)
         call SumBlock(VarsSrc%y, ModMap%iVarSrc(FieldOrientation), VarsDst%u, ModMap%iVarDst(FieldOrientation), Mapping%MeshMap%dM%mi, dUdy)
         call SumBlock(VarsSrc%y, ModMap%iVarSrc(FieldTransVel), VarsDst%u, ModMap%iVarDst(FieldTransVel), Mapping%MeshMap%dM%mi, dUdy)
         call SumBlock(VarsSrc%y, ModMap%iVarSrc(FieldAngularVel), VarsDst%u, ModMap%iVarDst(FieldAngularVel), Mapping%MeshMap%dM%mi, dUdy)
         call SumBlock(VarsSrc%y, ModMap%iVarSrc(FieldTransAcc), VarsDst%u, ModMap%iVarDst(FieldTransAcc), Mapping%MeshMap%dM%mi, dUdy)
         call SumBlock(VarsSrc%y, ModMap%iVarSrc(FieldAngularAcc), VarsDst%u, ModMap%iVarDst(FieldAngularAcc), Mapping%MeshMap%dM%mi, dUdy)
      end if

      ! Rotation to Translation
      if (allocated(Mapping%MeshMap%dM%fx_p)) then
         call SumBlock(VarsSrc%y, ModMap%iVarSrc(FieldOrientation), VarsDst%u, ModMap%iVarDst(FieldTransDisp), Mapping%MeshMap%dM%fx_p, dUdy)
         call SumBlock(VarsSrc%y, ModMap%iVarSrc(FieldAngularVel), VarsDst%u, ModMap%iVarDst(FieldTransVel), Mapping%MeshMap%dM%fx_p, dUdy)
         call SumBlock(VarsSrc%y, ModMap%iVarSrc(FieldAngularAcc), VarsDst%u, ModMap%iVarDst(FieldTransAcc), Mapping%MeshMap%dM%fx_p, dUdy)
      end if

      ! Translation displacement to Translation velocity
      if (allocated(Mapping%MeshMap%dM%tv_us)) then
         call SumBlock(VarsSrc%y, ModMap%iVarSrc(FieldTransDisp), VarsDst%u, ModMap%iVarDst(FieldTransVel), Mapping%MeshMap%dM%tv_us, dUdy)
      end if

      ! Translation displacement to Translation acceleration
      if (allocated(Mapping%MeshMap%dM%ta_us)) then
         call SumBlock(VarsSrc%y, ModMap%iVarSrc(FieldTransDisp), VarsDst%u, ModMap%iVarDst(FieldTransAcc), Mapping%MeshMap%dM%ta_us, dUdy)
      end if

      ! Angular velocity to Translation acceleration
      if (allocated(Mapping%MeshMap%dM%ta_rv)) then
         call SumBlock(VarsSrc%y, ModMap%iVarSrc(FieldAngularVel), VarsDst%u, ModMap%iVarDst(FieldTransAcc), Mapping%MeshMap%dM%ta_rv, dUdy)
      end if
   end subroutine

   subroutine SumBlock(VarArySrc, iVarSrc, VarAryDst, iVarDst, SrcM, DstM)
      type(ModVarType), intent(in)     :: VarArySrc(:), VarAryDst(:)
      integer(IntKi), intent(in)       :: iVarSrc, iVarDst
      real(R8Ki), intent(in)           :: SrcM(:, :)
      real(R8Ki), intent(inout)        :: DstM(:, :)

      ! If no variable index for source or destination, return
      if (iVarDst == 0 .or. iVarSrc == 0) return

      ! Get pointers to source and destination locations
      associate (iGluSrc => VarArySrc(iVarSrc)%iGlu, iGluDst => VarAryDst(iVarDst)%iGlu)

         ! Subtracts the source matrix from the destination sub-matrix
         associate (DstSubM => DstM(iGluDst(1):iGluDst(2), iGluSrc(1):iGluSrc(2)))
            DstSubM = DstSubM - SrcM
         end associate

      end associate
   end subroutine

   logical function Failed()
      Failed = ErrStat2 >= AbortErrLev
      if (Failed) call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   end function
end subroutine

subroutine VarUnpackInput(ModData, Var, ValAry, T, iInput, ErrStat, ErrMsg)
   type(ModDataType), intent(in)          :: ModData
   type(ModVarType), intent(in)           :: Var
   real(R8Ki), intent(in)                 :: ValAry(:)
   type(FAST_TurbineType), intent(inout)  :: T           !< Turbine data
   integer(IntKi), intent(in)             :: iInput
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg
   ErrStat = ErrID_None
   ErrMsg = ''
   select case (ModData%ID)
   case (Module_AD)
      call AD_VarUnpackInput(Var, ValAry, T%AD%Input(iInput)%rotors(ModData%Ins))
   case (Module_ADsk)
      call ADsk_VarUnpackInput(Var, ValAry, T%ADsk%Input(iInput))
   case (Module_BD)
      call BD_VarUnpackInput(Var, ValAry, T%BD%Input(iInput, ModData%Ins))
   case (Module_ED)
      call ED_VarUnpackInput(Var, ValAry, T%ED%Input(iInput))
   case (Module_SED)
      call SED_VarUnpackInput(Var, ValAry, T%SED%Input(iInput))
   case (Module_ExtLd)
      call ExtLd_VarUnpackInput(Var, ValAry, T%ExtLd%u)
   case (Module_ExtInfw)
      call ExtInfw_VarUnpackInput(Var, ValAry, T%ExtInfw%u)
   case (Module_ExtPtfm)
      call ExtPtfm_VarUnpackInput(Var, ValAry, T%ExtPtfm%Input(iInput))
   case (Module_FEAM)
      call FEAM_VarUnpackInput(Var, ValAry, T%FEAM%Input(iInput))
   case (Module_HD)
      call HydroDyn_VarUnpackInput(Var, ValAry, T%HD%Input(iInput))
   case (Module_IceD)
      call IceD_VarUnpackInput(Var, ValAry, T%IceD%Input(iInput, ModData%Ins))
   case (Module_IceF)
      call IceFloe_VarUnpackInput(Var, ValAry, T%IceF%Input(iInput))
   case (Module_IfW)
      call InflowWind_VarUnpackInput(Var, ValAry, T%IfW%Input(iInput))
   case (Module_MAP)
      call MAP_VarUnpackInput(Var, ValAry, T%MAP%Input(iInput))
   case (Module_MD)
      call MD_VarUnpackInput(Var, ValAry, T%MD%Input(iInput))
   case (Module_Orca)
      call Orca_VarUnpackInput(Var, ValAry, T%Orca%Input(iInput))
   case (Module_SD)
      call SD_VarUnpackInput(Var, ValAry, T%SD%Input(iInput))
   case (Module_SeaSt)
      call SeaSt_VarUnpackInput(Var, ValAry, T%SeaSt%Input(iInput))
   case (Module_SrvD)
      call SrvD_VarUnpackInput(Var, ValAry, T%SrvD%Input(iInput))
   case default
      call SetErrStat(ErrID_Fatal, "Unsupported module: "//ModData%Abbr, ErrStat, ErrMsg, "VarPackInput")
   end select
end subroutine

subroutine VarPackOutput(ModData, Var, ValAry, T, ErrStat, ErrMsg)
   type(ModDataType), intent(in)       :: ModData
   type(ModVarType), intent(in)        :: Var
   real(R8Ki), intent(inout)           :: ValAry(:)
   type(FAST_TurbineType), intent(in)  :: T           !< Turbine data
   integer(IntKi), intent(out)         :: ErrStat
   character(*), intent(out)           :: ErrMsg
   ErrStat = ErrID_None
   ErrMsg = ''
   select case (ModData%ID)
   case (Module_AD)
      call AD_VarPackOutput(Var, T%AD%y%rotors(ModData%Ins), ValAry)
   case (Module_ADsk)
      call ADsk_VarPackOutput(Var, T%ADsk%y, ValAry)
   case (Module_BD)
      call BD_VarPackOutput(Var, T%BD%y(ModData%Ins), ValAry)
   case (Module_ED)
      call ED_VarPackOutput(Var, T%ED%y, ValAry)
   case (Module_SED)
      call SED_VarPackOutput(Var, T%SED%y, ValAry)
   case (Module_ExtLd)
      call ExtLd_VarPackOutput(Var, T%ExtLd%y, ValAry)
   case (Module_ExtInfw)
      call ExtInfw_VarPackOutput(Var, T%ExtInfw%y, ValAry)
   case (Module_ExtPtfm)
      call ExtPtfm_VarPackOutput(Var, T%ExtPtfm%y, ValAry)
   case (Module_FEAM)
      call FEAM_VarPackOutput(Var, T%FEAM%y, ValAry)
   case (Module_HD)
      call HydroDyn_VarPackOutput(Var, T%HD%y, ValAry)
   case (Module_IceD)
      call IceD_VarPackOutput(Var, T%IceD%y(ModData%Ins), ValAry)
   case (Module_IceF)
      call IceFloe_VarPackOutput(Var, T%IceF%y, ValAry)
   case (Module_IfW)
      call InflowWind_VarPackOutput(Var, T%IfW%y, ValAry)
   case (Module_MAP)
      call MAP_VarPackOutput(Var, T%MAP%y, ValAry)
   case (Module_MD)
      call MD_VarPackOutput(Var, T%MD%y, ValAry)
   case (Module_Orca)
      call Orca_VarPackOutput(Var, T%Orca%y, ValAry)
   case (Module_SD)
      call SD_VarPackOutput(Var, T%SD%y, ValAry)
   case (Module_SeaSt)
      call SeaSt_VarPackOutput(Var, T%SeaSt%y, ValAry)
   case (Module_SrvD)
      call SrvD_VarPackOutput(Var, T%SrvD%y, ValAry)
   case default
      call SetErrStat(ErrID_Fatal, "Unsupported module: "//ModData%Abbr, ErrStat, ErrMsg, "VarPackOutput")
   end select
end subroutine

subroutine FAST_InputSolve(iModDst, ModAry, MapAry, iInput, Turbine, ErrStat, ErrMsg, VarMapAry)
   integer(IntKi), intent(in)                      :: iModDst     !< Destination module index in module data array
   type(ModDataType), intent(in)                   :: ModAry(:)   !< Module data
   type(MappingType), intent(inout)                :: MapAry(:)   !< Mesh and variable mappings
   integer(IntKi), intent(in)                      :: iInput      !< Input index to store data
   type(FAST_TurbineType), target, intent(inout)   :: Turbine     !< Turbine type
   integer(IntKi), intent(out)                     :: ErrStat
   character(*), intent(out)                       :: ErrMsg
   type(VarMapType), optional, intent(in)          :: VarMapAry(:)

   character(*), parameter                         :: RoutineName = 'FAST_InputSolve'
   integer(IntKi)                                  :: ErrStat2
   character(ErrMsgLen)                            :: ErrMsg2

   integer(IntKi)                                  :: i

   ErrStat = ErrID_None
   ErrMsg = ''

   if (present(VarMapAry)) then

      ! Loop through mappings and zero load meshes before transfer
      do i = 1, size(VarMapAry)
         associate (Mapping => MapAry(VarMapAry(i)%iMapping))

            ! Skip mappings where this isn't the destination module
            if (iModDst /= Mapping%iModDst) cycle

            ! Skip mappings that are not ready
            if (.not. Mapping%Ready) cycle
            
            ! If this is a load mesh mapping, clear the loads
            if (Mapping%MapType == Map_LoadMesh) call ZeroDstLoadMesh(Mapping, ModAry(VarMapAry(i)%iModDst))
         end associate
      end do

      ! Loop through mappings and perform input solve
      do i = 1, size(VarMapAry)
         associate (Mapping => MapAry(VarMapAry(i)%iMapping))

            ! Skip mappings where this isn't the destination module
            if (iModDst /= VarMapAry(i)%iModDst) cycle

            ! Skip mappings that are not ready
            if (.not. Mapping%Ready) cycle

            ! Perform input solve
            call InputSolveMapping(MapAry(VarMapAry(i)%iMapping), ModAry(VarMapAry(i)%iModSrc), ModAry(VarMapAry(i)%iModDst))
            if (ErrStat >= AbortErrLev) return
         end associate
      end do
   
   else

      ! Loop through mappings and zero load meshes before transfer
      do i = 1, size(MapAry)

         ! Skip mappings that are not ready
         if (.not. MapAry(i)%Ready) cycle

         ! Skip mappings where this isn't the destination module
         if (iModDst /= MapAry(i)%iModDst) cycle
         
         ! If this is a load mesh mapping, clear the loads
         if (MapAry(i)%MapType == Map_LoadMesh) call ZeroDstLoadMesh(MapAry(i), ModAry(MapAry(i)%iModDst))
      end do

      ! Loop through mappings and perform input solve
      do i = 1, size(MapAry)

         ! Skip mappings where this isn't the destination module
         if (iModDst /= MapAry(i)%iModDst) cycle

         ! Skip mappings that are not ready
         if (.not. MapAry(i)%Ready) cycle

         ! Perform input solve
         call InputSolveMapping(MapAry(i), ModAry(MapAry(i)%iModSrc), ModAry(MapAry(i)%iModDst))
         if (ErrStat >= AbortErrLev) return
      end do
   end if

contains

   subroutine ZeroDstLoadMesh(Mapping, ModDst)
      type(MappingType), intent(inout) :: Mapping
      type(ModDataType), intent(in)    :: ModDst
      type(MeshType), pointer          :: DstMesh

      ! Get pointer to destination load mesh
      call FAST_InputMeshPointer(ModDst, Turbine, Mapping%DstDL, DstMesh, iInput, ErrStat2, ErrMsg2)
      if (Failed()) return

      ! If mesh has force, set it to zero
      if (DstMesh%fieldmask(MASKID_FORCE)) DstMesh%Force = 0.0_ReKi

      ! If mesh has moment, set it to zero
      if (DstMesh%fieldmask(MASKID_MOMENT)) DstMesh%Moment = 0.0_ReKi

   end subroutine

   subroutine InputSolveMapping(Mapping, ModSrc, ModDst)
      type(MappingType), intent(inout) :: Mapping
      type(ModDataType), intent(in)    :: ModSrc, ModDst
      type(MeshType), pointer          :: SrcMesh, DstMesh
      type(MeshType), pointer          :: SrcDispMesh, DstDispMesh

      ! Select based on type of mapping
      select case (Mapping%MapType)

      case (Map_Custom)

         call Custom_InputSolve(Mapping, ModSrc, ModDst, iInput, Turbine, ErrStat2, ErrMsg2)
         if (Failed()) return

      case (Map_Variable)

         ! Pack module output value into array
         call VarPackOutput(ModSrc, Mapping%SrcVar, Mapping%VarData, Turbine, ErrStat2, ErrMsg2)
         if (Failed()) return

         ! If fewer source values than destination values, copy first value to all values
         if (Mapping%SrcVar%Num < Mapping%DstVar%Num) then
            Mapping%VarData = Mapping%VarData(1)
         end if

         ! Unpack array into module input
         call VarUnpackInput(ModDst, Mapping%DstVar, Mapping%VarData, Turbine, iInput, ErrStat2, ErrMsg2)
         if (Failed()) return

      case (Map_MotionMesh)

         ! Get source and destination meshes
         call FAST_OutputMeshPointer(ModSrc, Turbine, Mapping%SrcDL, SrcMesh, ErrStat2, ErrMsg2)
         if (Failed()) return
         call FAST_InputMeshPointer(ModDst, Turbine, Mapping%DstDL, DstMesh, iInput, ErrStat2, ErrMsg2)
         if (Failed()) return

         ! Perform transfer based on type
         call TransferMesh(Mapping%XfrType, SrcMesh, DstMesh, Mapping%MeshMap, ErrStat=ErrStat2, ErrMsg=ErrMsg2)
         if (Failed()) return

      case (Map_LoadMesh)

         ! Get source and destination meshes
         call FAST_OutputMeshPointer(ModSrc, Turbine, Mapping%SrcDL, SrcMesh, ErrStat2, ErrMsg2)
         if (Failed()) return
         call FAST_InputMeshPointer(ModDst, Turbine, Mapping%DstDL, DstMesh, iInput, ErrStat2, ErrMsg2)
         if (Failed()) return

         ! Get source and destination displacement meshes
         ! Note: Displacement meshes always references current input index when in input
         call FAST_InputMeshPointer(ModSrc, Turbine, Mapping%SrcDispDL, SrcDispMesh, INPUT_CURR, ErrStat2, ErrMsg2)
         if (Failed()) return
         call FAST_OutputMeshPointer(ModDst, Turbine, Mapping%DstDispDL, DstDispMesh, ErrStat2, ErrMsg2)
         if (ErrStat2 == ErrID_Fatal) call FAST_InputMeshPointer(ModDst, Turbine, Mapping%DstDispDL, DstDispMesh, INPUT_CURR, ErrStat2, ErrMsg2)
         if (Failed()) return

         ! If DstDispMesh is a sibling of DstMesh
         if (Mapping%DstUsesSibling) then

            ! Transfer the load mesh to the temporary load mesh to be summed below
            call TransferMesh(Mapping%XfrType, SrcMesh, Mapping%TmpLoadMesh, Mapping%MeshMap, SrcDispMesh, DstDispMesh, ErrStat2, ErrMsg2)
            if (Failed()) return

         else

            ! Transfer destination displacement mesh to temporary motion mesh (cousin of destination load mesh)
            call TransferMesh(Mapping%XfrTypeAux, DstDispMesh, Mapping%TmpMotionMesh, Mapping%MeshMapAux, ErrStat=ErrStat2, ErrMsg=ErrMsg2)
            if (Failed()) return
            
            ! Transfer to temporary load mesh using the temporary motion mesh as the destination displacement mesh
            call TransferMesh(Mapping%XfrType, SrcMesh, Mapping%TmpLoadMesh, Mapping%MeshMap, SrcDispMesh, Mapping%TmpMotionMesh, ErrStat2, ErrMsg2)
            if (Failed()) return

         end if

         ! Sum loads from temporary mesh to destination mesh
         if (DstMesh%fieldmask(MASKID_FORCE)) DstMesh%Force = DstMesh%Force + Mapping%TmpLoadMesh%Force
         if (DstMesh%fieldmask(MASKID_MOMENT)) DstMesh%Moment = DstMesh%Moment + Mapping%TmpLoadMesh%Moment

      end select

   end subroutine

   logical function Failed()
      Failed = ErrStat2 /= ErrID_None
      if (Failed) call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, &
                                  RoutineName//':Module='//trim(ModAry(iModDst)%Abbr)// &
                                  ', Instance='//Num2LStr(ModAry(iModDst)%Ins))
   end function
end subroutine

! Reset mapping read flags
subroutine FAST_ResetMappingReady(MapAry)
   type(MappingType), intent(inout) :: MapAry(:)   !< Mesh and variable mappings
   integer(IntKi)                   :: i
   do i = 1, size(MapAry)
      select case (MapAry(i)%SrcModID)
      case default         ! Default to transfer is not ready
         MapAry(i)%Ready = .false.
      end select
   end do
end subroutine

! TransferMesh calls the specific transfer function based on
! transfer type (Point_to_Point, Point_to_Line2, etc.)
subroutine TransferMesh(Typ, Src, Dst, MeshMap, SrcDisp, DstDisp, ErrStat, ErrMsg)
   integer(IntKi), intent(in)             :: Typ
   type(MeshType), intent(in)             :: Src
   type(MeshType), intent(inout)          :: Dst
   type(MeshMapType), intent(inout)       :: MeshMap
   type(MeshType), optional, intent(in)   :: SrcDisp, DstDisp
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg
   select case (Typ)
   case (Xfr_Point_to_Point)
      call Transfer_Point_to_Point(Src, Dst, MeshMap, ErrStat, ErrMsg, SrcDisp, DstDisp)
   case (Xfr_Point_to_Line2)
      call Transfer_Point_to_Line2(Src, Dst, MeshMap, ErrStat, ErrMsg, SrcDisp, DstDisp)
   case (Xfr_Line2_to_Point)
      call Transfer_Line2_to_Point(Src, Dst, MeshMap, ErrStat, ErrMsg, SrcDisp, DstDisp)
   case (Xfr_Line2_to_Line2)
      call Transfer_Line2_to_Line2(Src, Dst, MeshMap, ErrStat, ErrMsg, SrcDisp, DstDisp)
   case default
      ErrStat = ErrID_Fatal
      ErrMsg = "TransferMesh: unknown transfer type: "//Num2LStr(Typ)
   end select
end subroutine

subroutine Custom_InputSolve(Mapping, ModSrc, ModDst, iInput, T, ErrStat, ErrMsg)
   type(MappingType), intent(in)          :: Mapping
   type(ModDataType), intent(in)          :: ModSrc, ModDst
   integer(IntKi), intent(in)             :: iInput
   type(FAST_TurbineType), intent(inout)  :: T           !< Turbine type
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg

   character(*), parameter                :: RoutineName = 'Custom_InputSolve'
   integer(IntKi)                         :: ErrStat2
   character(ErrMsgLen)                   :: ErrMsg2
   integer(IntKi)                         :: i, j, k
   real(ReKi)                             :: z, u, v, mean_vel

   ErrStat = ErrID_None
   ErrMsg = ''

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
      do i = 1, size(T%AD%Input(iInput)%rotors(ModDst%Ins)%UserProp, dim=2)   ! Blade
         T%AD%Input(iInput)%rotors(ModDst%Ins)%UserProp(:, i) = T%SrvD%y%BlAirfoilCom(i)
      end do

!-------------------------------------------------------------------------------
! ADsk Inputs
!-------------------------------------------------------------------------------

   case (Custom_ED_to_ADsk)

      T%ADsk%Input(iInput)%RotSpeed   = T%ED%y%RotSpeed
      T%ADsk%Input(iInput)%BlPitch    = T%ED%y%BlPitch(1)   ! ADsk only uses collective blade pitch

   case (Custom_SED_to_ADsk)

      T%ADsk%Input(iInput)%RotSpeed   = T%SED%y%RotSpeed
      T%ADsk%Input(iInput)%BlPitch    = T%SED%y%BlPitch(1)   ! ADsk only uses collective blade pitch

!-------------------------------------------------------------------------------
! ElastoDyn Inputs
!-------------------------------------------------------------------------------

   case (Custom_SrvD_to_ED)

      T%ED%Input(iInput)%GenTrq = T%SrvD%y%GenTrq
      T%ED%Input(iInput)%HSSBrTrqC = T%SrvD%y%HSSBrTrqC
      T%ED%Input(iInput)%BlPitchCom = T%SrvD%y%BlPitchCom
      T%ED%Input(iInput)%YawMom = T%SrvD%y%YawMom

!-------------------------------------------------------------------------------
! SED Inputs
!-------------------------------------------------------------------------------

   case (Custom_SrvD_to_SED)

      T%SED%Input(iInput)%GenTrq = T%SrvD%y%GenTrq
      T%SED%Input(iInput)%HSSBrTrqC = T%SrvD%y%HSSBrTrqC
      T%SED%Input(iInput)%BlPitchCom = T%SrvD%y%BlPitchCom
      T%SED%Input(iInput)%YawPosCom = T%SrvD%y%YawPosCom
      T%SED%Input(iInput)%YawRateCom = T%SrvD%y%YawRateCom

!-------------------------------------------------------------------------------
! ExtLoads Inputs
!-------------------------------------------------------------------------------

   case (Custom_ED_to_ExtLd)

      T%ExtLd%u%az = T%ED%y%LSSTipPxa
      T%ExtLd%u%DX_u%bldPitch(:) = T%ED%y%BlPitch

      ! Note: this may be better inside CalcOutput
      call ExtLd_ConvertInpDataForExtProg(T%ExtLd%u, T%ExtLd%p, ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) return

!-------------------------------------------------------------------------------
! InflowWind Inputs
!-------------------------------------------------------------------------------

   case (Custom_ED_to_IfW)

      ! This section should be refactored so that IfW uses a hub point mesh
      T%IfW%Input(iInput)%HubPosition = T%ED%y%HubPtMotion%Position(:, 1) + &
                                        T%ED%y%HubPtMotion%TranslationDisp(:, 1)
      T%IfW%Input(iInput)%HubOrientation = T%ED%y%HubPtMotion%Orientation(:, :, 1)

      ! Set Lidar position directly from hub motion mesh
      T%IfW%Input(iInput)%lidar%HubDisplacementX = T%ED%y%HubPtMotion%TranslationDisp(1, 1)
      T%IfW%Input(iInput)%lidar%HubDisplacementY = T%ED%y%HubPtMotion%TranslationDisp(2, 1)
      T%IfW%Input(iInput)%lidar%HubDisplacementZ = T%ED%y%HubPtMotion%TranslationDisp(3, 1)

   case (Custom_SrvD_to_IfW)

!-------------------------------------------------------------------------------
! MoorDyn Inputs
!-------------------------------------------------------------------------------

   case (Custom_SrvD_to_MD)

      if (allocated(T%MD%Input(iInput)%DeltaL) .and. allocated(T%SrvD%y%CableDeltaL)) then
         T%MD%Input(iInput)%DeltaL = T%SrvD%y%CableDeltaL      ! these should be sized identically during init
      end if

      if (allocated(T%MD%Input(iInput)%DeltaLdot) .and. allocated(T%SrvD%y%CableDeltaLdot)) then
         T%MD%Input(iInput)%DeltaLdot = T%SrvD%y%CableDeltaLdot   ! these should be sized identically during init
      end if

!-------------------------------------------------------------------------------
! SubDyn Inputs
!-------------------------------------------------------------------------------

   case (Custom_SrvD_to_SD)

      if (allocated(T%SD%Input(iInput)%CableDeltaL) .and. allocated(T%SrvD%y%CableDeltaL)) then
         T%SD%Input(iInput)%CableDeltaL = T%SrvD%y%CableDeltaL   ! these should be sized identically during init
      end if

!-------------------------------------------------------------------------------
! ServoDyn Inputs
!-------------------------------------------------------------------------------

   case (Custom_BD_to_SrvD)

      T%SrvD%Input(iInput)%RootMxc(Mapping%SrcIns) = T%BD%y(Mapping%SrcIns)%RootMxr*cos(T%ED%y%BlPitch(Mapping%SrcIns)) + &
                                                     T%BD%y(Mapping%SrcIns)%RootMyr*sin(T%ED%y%BlPitch(Mapping%SrcIns))
      T%SrvD%Input(iInput)%RootMyc(Mapping%SrcIns) = -T%BD%y(Mapping%SrcIns)%RootMxr*sin(T%ED%y%BlPitch(Mapping%SrcIns)) + &
                                                     T%BD%y(Mapping%SrcIns)%RootMyr*cos(T%ED%y%BlPitch(Mapping%SrcIns))

   case (Custom_ED_to_SrvD)

      ! Blade root moment if not using BeamDyn
      if (T%p_FAST%CompElast /= Module_BD) then
         T%SrvD%Input(iInput)%RootMxc = T%ED%y%RootMxc ! fixed-size arrays: always size 3
         T%SrvD%Input(iInput)%RootMyc = T%ED%y%RootMyc ! fixed-size arrays: always size 3
      end if

      T%SrvD%Input(iInput)%YawAngle = T%ED%y%YawAngle ! nacelle yaw plus platform yaw
      T%SrvD%Input(iInput)%YawErr = T%SrvD%Input(iInput)%WindDir - T%SrvD%Input(iInput)%YawAngle ! the nacelle yaw error estimate (positive about zi-axis)

      T%SrvD%Input(iInput)%BlPitch = T%ED%y%BlPitch
      T%SrvD%Input(iInput)%LSS_Spd = T%ED%y%LSS_Spd
      T%SrvD%Input(iInput)%RotSpeed = T%ED%y%RotSpeed

      T%SrvD%Input(iInput)%YawBrTAxp = T%ED%y%YawBrTAxp
      T%SrvD%Input(iInput)%YawBrTAyp = T%ED%y%YawBrTAyp
      T%SrvD%Input(iInput)%LSSTipPxa = T%ED%y%LSSTipPxa

      T%SrvD%Input(iInput)%LSSTipMxa = T%ED%y%LSSTipMxa
      T%SrvD%Input(iInput)%LSSTipMya = T%ED%y%LSSTipMya
      T%SrvD%Input(iInput)%LSSTipMza = T%ED%y%LSSTipMza
      T%SrvD%Input(iInput)%LSSTipMys = T%ED%y%LSSTipMys
      T%SrvD%Input(iInput)%LSSTipMzs = T%ED%y%LSSTipMzs

      T%SrvD%Input(iInput)%YawBrMyn = T%ED%y%YawBrMyn
      T%SrvD%Input(iInput)%YawBrMzn = T%ED%y%YawBrMzn
      T%SrvD%Input(iInput)%NcIMURAxs = T%ED%y%NcIMURAxs
      T%SrvD%Input(iInput)%NcIMURAys = T%ED%y%NcIMURAys
      T%SrvD%Input(iInput)%NcIMURAzs = T%ED%y%NcIMURAzs

      T%SrvD%Input(iInput)%RotPwr = T%ED%y%RotPwr

      T%SrvD%Input(iInput)%LSShftFxa = T%ED%y%LSShftFxa
      T%SrvD%Input(iInput)%LSShftFys = T%ED%y%LSShftFys
      T%SrvD%Input(iInput)%LSShftFzs = T%ED%y%LSShftFzs

   case (Custom_IfW_to_SrvD)

      T%SrvD%Input(iInput)%WindDir = atan2(T%IfW%y%HubVel(2), T%IfW%y%HubVel(1))
      T%SrvD%Input(iInput)%HorWindV = sqrt(T%IfW%y%HubVel(1)**2 + T%IfW%y%HubVel(2)**2)
      if (allocated(T%IfW%y%lidar%LidSpeed)) T%SrvD%Input(iInput)%LidSpeed = T%IfW%y%lidar%LidSpeed
      if (allocated(T%IfW%y%lidar%MsrPositionsX)) T%SrvD%Input(iInput)%MsrPositionsX = T%IfW%y%lidar%MsrPositionsX
      if (allocated(T%IfW%y%lidar%MsrPositionsY)) T%SrvD%Input(iInput)%MsrPositionsY = T%IfW%y%lidar%MsrPositionsY
      if (allocated(T%IfW%y%lidar%MsrPositionsZ)) T%SrvD%Input(iInput)%MsrPositionsZ = T%IfW%y%lidar%MsrPositionsZ
      T%SrvD%Input(iInput)%YawErr = T%SrvD%Input(iInput)%WindDir - T%SrvD%Input(iInput)%YawAngle ! the nacelle yaw error estimate (positive about zi-axis)

   case (Custom_ExtInfw_to_SrvD)

      T%SrvD%Input(iInput)%WindDir = ATAN2(T%ExtInfw%y%v(1), T%ExtInfw%y%u(1))
      T%SrvD%Input(iInput)%HorWindV = SQRT(T%ExtInfw%y%u(1)**2 + T%ExtInfw%y%v(1)**2)
      if (allocated(T%SrvD%Input(iInput)%LidSpeed)) T%SrvD%Input(iInput)%LidSpeed = 0.0
      if (allocated(T%SrvD%Input(iInput)%MsrPositionsX)) T%SrvD%Input(iInput)%MsrPositionsX = 0.0
      if (allocated(T%SrvD%Input(iInput)%MsrPositionsY)) T%SrvD%Input(iInput)%MsrPositionsY = 0.0
      if (allocated(T%SrvD%Input(iInput)%MsrPositionsz)) T%SrvD%Input(iInput)%MsrPositionsz = 0.0
      
      ! the nacelle yaw error estimate (positive about zi-axis)
      T%SrvD%Input(iInput)%YawErr = T%SrvD%Input(iInput)%WindDir - T%SrvD%Input(iInput)%YawAngle 

   case (Custom_ExtLd_to_SrvD)

      pi = acos(-1.0)
      z = T%ED%y%HubPtMotion%Position(3, 1)
      mean_vel = T%ExtLd%p%vel_mean*((z/T%ExtLd%p%z_ref)**T%ExtLd%p%shear_exp)
      u = -mean_vel*sin(T%ExtLd%p%wind_dir*pi/180.0)
      v = -mean_vel*cos(T%ExtLd%p%wind_dir*pi/180.0)
      T%SrvD%Input(iInput)%HorWindV = mean_vel
      T%SrvD%Input(iInput)%WindDir = atan2(v, u)
      if (allocated(T%SrvD%Input(iInput)%LidSpeed)) T%SrvD%Input(iInput)%LidSpeed = 0.0
      if (allocated(T%SrvD%Input(iInput)%MsrPositionsX)) T%SrvD%Input(iInput)%MsrPositionsX = 0.0
      if (allocated(T%SrvD%Input(iInput)%MsrPositionsY)) T%SrvD%Input(iInput)%MsrPositionsY = 0.0
      if (allocated(T%SrvD%Input(iInput)%MsrPositionsz)) T%SrvD%Input(iInput)%MsrPositionsz = 0.0
      
      ! the nacelle yaw error estimate (positive about zi-axis)
      T%SrvD%Input(iInput)%YawErr = T%SrvD%Input(iInput)%WindDir - T%SrvD%Input(iInput)%YawAngle 

!-------------------------------------------------------------------------------
! Unknown Mapping
!-------------------------------------------------------------------------------

   case default

      ErrStat = ErrID_Fatal
      ErrMsg = "Custom_InputSolve: unknown mapping '"//trim(Mapping%Desc)//"'"

   end select

end subroutine

subroutine FAST_ResetRemapFlags(Mods, Maps, T, ErrStat, ErrMsg)
   type(ModDataType), intent(in)                   :: Mods(:) !< Module data
   type(MappingType), intent(inout)                :: Maps(:)
   type(FAST_TurbineType), target, intent(inout)   :: T       !< Turbine type
   integer(IntKi), intent(out)                     :: ErrStat
   character(*), intent(out)                       :: ErrMsg

   character(*), parameter    :: RoutineName = 'FAST_ResetRemapFlags'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2
   integer(IntKi)             :: i, k
   type(MeshType), pointer    :: SrcMesh, DstMesh

   ErrStat = ErrID_None
   ErrMsg = ''

   ! Reset remap flags in mapping meshes
   do i = 1, size(Maps)
      select case (Maps(i)%MapType)
      case (Map_LoadMesh, Map_MotionMesh)

         if (associated(Maps(i)%TmpLoadMesh%RemapFlag)) Maps(i)%TmpLoadMesh%RemapFlag = .false.
         if (associated(Maps(i)%TmpMotionMesh%RemapFlag)) Maps(i)%TmpMotionMesh%RemapFlag = .false.

         call FAST_OutputMeshPointer(Mods(Maps(i)%iModSrc), T, Maps(i)%SrcDL, SrcMesh, ErrStat2, ErrMsg2)
         if (Failed()) return
         SrcMesh%RemapFlag = .false.

         call FAST_InputMeshPointer(Mods(Maps(i)%iModDst), T, Maps(i)%DstDL, DstMesh, INPUT_CURR, ErrStat2, ErrMsg2)
         if (Failed()) return
         DstMesh%RemapFlag = .false.

      end select
   end do

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

end module
