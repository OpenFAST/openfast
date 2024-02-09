module FAST_Mesh

use FAST_ModTypes

implicit none

private
public :: FAST_InitMappings, FAST_LinearizeMappings, FAST_ResetRemapFlags, FAST_InputSolve

! Input Solve destinations
integer(IntKi), parameter  :: IS_DstInput = 1, IS_Dstu = 2, IS_Linearize = 3

integer(IntKi), parameter  :: AD_rotor = 1

integer(IntKi), parameter  :: Xfr_Point_to_Point = 1, &
                              Xfr_Line2_to_Point = 2, &
                              Xfr_Point_to_Line2 = 3, &
                              Xfr_Line2_to_Line2 = 4

contains

function FAST_InputMeshPointer(ModData, Turbine, MeshLoc, UseU) result(Mesh)
   type(ModDataType), intent(in)                :: ModData
   type(FAST_TurbineType), target, intent(in)   :: Turbine
   type(MeshLocType), intent(in)                :: MeshLoc
   logical, intent(in)                          :: UseU
   type(MeshType), pointer                      :: Mesh

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
   case (Module_SD)
      if (UseU) then
         Mesh => SD_InputMeshPointer(Turbine%SD%u, MeshLoc)
      else
         Mesh => SD_InputMeshPointer(Turbine%SD%Input(1), MeshLoc)
      end if
   case (Module_SrvD)
      if (UseU) then
         Mesh => SrvD_InputMeshPointer(Turbine%SrvD%u, MeshLoc)
      else
         Mesh => SrvD_InputMeshPointer(Turbine%SrvD%Input(1), MeshLoc)
      end if
   end select
end function

function FAST_InputMeshName(ModData, Turbine, MeshLoc) result(Name)
   type(ModDataType), intent(in)                :: ModData
   type(FAST_TurbineType), target, intent(in)   :: Turbine
   type(MeshLocType), intent(in)                :: MeshLoc
   character(32)                                :: Name
   select case (ModData%ID)
   case (Module_AD)
      Name = trim(ModData%Abbr)//"%"//AD_InputMeshName(Turbine%AD%u, MeshLoc)
   case (Module_BD)
      Name = trim(ModData%Abbr)//"("//trim(Num2LStr(ModData%Ins))//")%"//BD_InputMeshName(Turbine%BD%u(ModData%Ins), MeshLoc)
   case (Module_ED)
      Name = trim(ModData%Abbr)//"%"//ED_InputMeshName(Turbine%ED%u, MeshLoc)
   case (Module_SD)
      Name = trim(ModData%Abbr)//"%"//SD_InputMeshName(Turbine%SD%u, MeshLoc)
   case (Module_SrvD)
      Name = trim(ModData%Abbr)//"%"//SrvD_InputMeshName(Turbine%SrvD%u, MeshLoc)
   end select
end function

function FAST_OutputMeshPointer(ModData, Turbine, MeshLoc) result(Mesh)
   type(ModDataType), intent(in)                :: ModData
   type(FAST_TurbineType), target, intent(in)   :: Turbine
   type(MeshLocType), intent(in)                :: MeshLoc
   type(MeshType), pointer                      :: Mesh
   select case (ModData%ID)
   case (Module_AD)
      Mesh => AD_OutputMeshPointer(Turbine%AD%y, MeshLoc)
   case (Module_BD)
      Mesh => BD_OutputMeshPointer(Turbine%BD%y(ModData%Ins), MeshLoc)
   case (Module_ED)
      Mesh => ED_OutputMeshPointer(Turbine%ED%y, MeshLoc)
   case (Module_SD)
      Mesh => SD_OutputMeshPointer(Turbine%SD%y, MeshLoc)
   case (Module_SrvD)
      Mesh => SrvD_OutputMeshPointer(Turbine%SrvD%y, MeshLoc)
   end select
end function

function FAST_OutputMeshName(ModData, Turbine, MeshLoc) result(Name)
   type(ModDataType), intent(in)                :: ModData
   type(FAST_TurbineType), target, intent(in)   :: Turbine
   type(MeshLocType), intent(in)                :: MeshLoc
   character(32)                                :: Name
   select case (ModData%ID)
   case (Module_AD)
      Name = trim(ModData%Abbr)//"%"//AD_OutputMeshName(Turbine%AD%y, MeshLoc)
   case (Module_BD)
      Name = trim(ModData%Abbr)//"("//trim(Num2LStr(ModData%Ins))//")%"//BD_OutputMeshName(Turbine%BD%y(ModData%Ins), MeshLoc)
   case (Module_ED)
      Name = trim(ModData%Abbr)//"%"//ED_OutputMeshName(Turbine%ED%y, MeshLoc)
   case (Module_SD)
      Name = trim(ModData%Abbr)//"%"//SD_OutputMeshName(Turbine%SD%y, MeshLoc)
   case (Module_SrvD)
      Name = trim(ModData%Abbr)//"%"//SrvD_OutputMeshName(Turbine%SrvD%y, MeshLoc)
   end select
end function

subroutine FAST_InitMappings(Mods, Mappings, T, ErrStat, ErrMsg)
   type(ModDataType), intent(inout)                   :: Mods(:)     !< Module data
   type(TC_MappingType), allocatable, intent(inout)   :: Mappings(:)
   type(FAST_TurbineType), intent(inout)              :: T           !< Turbine type
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
            call AD_InitInputMappings(Mappings, Mods(iModSrc), Mods(iModDst), T, ErrStat, ErrMsg)
         case (Module_BD)
            call BD_InitInputMappings(Mappings, Mods(iModSrc), Mods(iModDst), T, ErrStat, ErrMsg)
         case (Module_ED)
            call ED_InitInputMappings(Mappings, Mods(iModSrc), Mods(iModDst), T, ErrStat, ErrMsg)
         case (Module_HD)
            call HD_InitInputMappings(Mappings, Mods(iModSrc), Mods(iModDst), T, ErrStat, ErrMsg)
         case (Module_IfW)
            call IfW_InitInputMappings(Mappings, Mods(iModSrc), Mods(iModDst), T, ErrStat, ErrMsg)
         case (Module_SD)
            call IfW_InitInputMappings(Mappings, Mods(iModSrc), Mods(iModDst), T, ErrStat, ErrMsg)
         case (Module_SrvD)
            call SrvD_InitInputMappings(Mappings, Mods(iModSrc), Mods(iModDst), T, ErrStat, ErrMsg)
         end select
      end do
   end do

   !----------------------------------------------------------------------------
   ! Get module indices in ModData and determine which mappings are active
   !----------------------------------------------------------------------------

   ! Reorder the mappings so that motion maps come before the load maps
   Mappings = [pack(Mappings, Mappings%MapType == Map_MotionMesh), &
               pack(Mappings, Mappings%MapType == Map_LoadMesh), &
               pack(Mappings, Mappings%MapType == Map_NonMesh)]

   ! Loop through mappings
   do iMap = 1, size(Mappings)
      associate (SrcMod => Mods(Mappings(iMap)%SrcModIdx), &
                 DstMod => Mods(Mappings(iMap)%DstModIdx))

         ! Add mapping index to sorce and destination module mapping arrays
         SrcMod%SrcMaps = [SrcMod%SrcMaps, iMap]
         DstMod%DstMaps = [DstMod%DstMaps, iMap]

      end associate
   end do

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine AD_InitInputMappings(Mappings, SrcMod, DstMod, Turbine, ErrStat, ErrMsg)
   type(TC_MappingType), allocatable      :: Mappings(:)
   type(ModDataType), intent(in)          :: SrcMod, DstMod
   type(FAST_TurbineType), intent(inout)  :: Turbine
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg

   character(*), parameter    :: RoutineName = 'AD_InitInputMappings'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2
   integer(IntKi)             :: i

   ErrStat = ErrID_None
   ErrMsg = ''

   select case (SrcMod%ID)
   case (Module_ED)

      call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                         SrcMeshLoc=MeshLocType(ED_y_TowerLn2Mesh), &                  ! T%ED%y%TowerLn2Mesh
                         DstMeshLoc=MeshLocType(AD_u_rotors_TowerMotion, AD_rotor), &  ! T%AD%Input(1)%rotors(1)%TowerMotion
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2, &
                         Active=Turbine%AD%Input(1)%rotors(1)%TowerMotion%Committed)
      if (Failed()) return

      call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                         SrcMeshLoc=MeshLocType(ED_y_HubPtMotion), &                   ! T%ED%y%HubPtMotion
                         DstMeshLoc=MeshLocType(AD_u_rotors_HubMotion, AD_rotor), &    ! T%AD%Input(1)%rotors(1)%HubMotion
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2)
      if (Failed()) return

      call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                         SrcMeshLoc=MeshLocType(ED_y_NacelleMotion), &                    ! T%ED%y%NacelleMotion
                         DstMeshLoc=MeshLocType(AD_u_rotors_NacelleMotion, AD_rotor), &   ! T%AD%Input(1)%rotors(1)%NacelleMotion
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2, &
                         Active=Turbine%AD%Input(1)%rotors(1)%NacelleMotion%Committed)
      if (Failed()) return

      call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                         SrcMeshLoc=MeshLocType(ED_y_TFinCMMotion), &                  ! T%ED%y%TFinCMMotion
                         DstMeshLoc=MeshLocType(AD_u_rotors_TFinMotion, AD_rotor), &   ! T%AD%Input(1)%rotors(1)%TFinMotion
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2, &
                         Active=Turbine%AD%Input(1)%rotors(1)%TFinMotion%Committed)
      if (Failed()) return

      do i = 1, size(Turbine%ED%y%BladeRootMotion)
         call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                            SrcMeshLoc=MeshLocType(ED_y_BladeRootMotion, i), &                  ! T%ED%y%BladeRootMotion(i)
                            DstMeshLoc=MeshLocType(AD_u_rotors_BladeRootMotion, AD_rotor, i), & ! T%AD%Input(1)%rotors(1)%BladeRootMotion(i)
                            ErrStat=ErrStat2, ErrMsg=ErrMsg2)
         if (Failed()) return
      end do

      if (Turbine%p_FAST%CompElast == Module_ED) then
         do i = 1, size(Turbine%ED%y%BladeLn2Mesh)
            call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                               SrcMeshLoc=MeshLocType(ED_y_BladeLn2Mesh, i), &                  ! T%ED%y%BladeLn2Mesh(i)
                               DstMeshLoc=MeshLocType(AD_u_rotors_BladeMotion, AD_rotor, i), &  ! T%AD%Input(1)%rotors(1)%BladeMotion(i)
                               ErrStat=ErrStat2, ErrMsg=ErrMsg2)
            if (Failed()) return
         end do
      end if

   case (Module_BD)

      call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                         SrcMeshLoc=MeshLocType(BD_y_BldMotion), &                                 ! BD%y(SrcMod%Ins)%BldMotion
                         DstMeshLoc=MeshLocType(AD_u_rotors_BladeMotion, AD_rotor, SrcMod%Ins), &  ! AD%Input(1)%rotors(1)%BladeMotion(SrcMod%Ins)
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2)
      if (Failed()) return

   case (Module_SrvD)

      call NonMeshMap(Mappings, Key='SrvD BlAirfoilCom -> AD UserProp', SrcMod=SrcMod, DstMod=DstMod)

   end select

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine BD_InitInputMappings(Mappings, SrcMod, DstMod, Turbine, ErrStat, ErrMsg)
   type(TC_MappingType), allocatable      :: Mappings(:)
   type(ModDataType), intent(in)          :: SrcMod, DstMod
   type(FAST_TurbineType), intent(inout)  :: Turbine           !< Turbine type
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg

   character(*), parameter    :: RoutineName = 'BD_InitInputMappings'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2
   integer(IntKi)             :: i

   ErrStat = ErrID_None
   ErrMsg = ''

   select case (SrcMod%ID)
   case (Module_ED)

      call MapMotionMesh(Turbine, Mappings, &
                         SrcMod=SrcMod, SrcMeshLoc=MeshLocType(ED_y_BladeRootMotion, DstMod%Ins), &   ! T%ED%y%BladeRootMotion(DstMod%Ins)
                         DstMod=DstMod, DstMeshLoc=MeshLocType(BD_u_RootMotion), &                    ! T%BD%Input(1, DstMod%Ins)%RootMotion
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2)
      if (Failed()) return

   case (Module_AD)

      call MapLoadMesh(Turbine, Mappings, &
                       SrcMod=SrcMod, SrcMeshLoc=MeshLocType(AD_y_rotors_BladeLoad, AD_rotor, DstMod%Ins), &   ! T%AD%y%rotors(1)%BladeLoad(DstMod%Ins)
                       SrcDispMeshLoc=MeshLocType(AD_u_rotors_BladeMotion, AD_rotor, DstMod%Ins), &            ! AD%Input(1)%rotors(1)%BladeMotion(DstMod%Ins)
                       DstMod=DstMod, DstMeshLoc=MeshLocType(BD_u_DistrLoad), &                                ! BD%Input(1, DstMod%Ins)%DistrLoad
                       DstDispMeshLoc=MeshLocType(BD_y_BldMotion), &                                           ! BD%y(DstMod%Ins)%BldMotion
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2)
      if (Failed()) return

   case (Module_SrvD)

   end select

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine ED_InitInputMappings(Mappings, SrcMod, DstMod, Turbine, ErrStat, ErrMsg)
   type(TC_MappingType), allocatable      :: Mappings(:)
   type(ModDataType), intent(in)          :: SrcMod, DstMod
   type(FAST_TurbineType), intent(inout)  :: Turbine           !< Turbine type
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg

   character(*), parameter    :: RoutineName = 'ED_InitInputMappings'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2
   integer(IntKi)             :: i, j

   ErrStat = ErrID_None
   ErrMsg = ''

   select case (SrcMod%ID)
   case (Module_BD)

      call MapLoadMesh(Turbine, Mappings, &
                       SrcMod=SrcMod, &
                       SrcMeshLoc=MeshLocType(BD_y_ReactionForce), &    ! BD%y(SrcMod%Ins)%ReactionForce
                       SrcDispMeshLoc=MeshLocType(BD_u_RootMotion), &   ! BD%Input(1, SrcMod%Ins)%RootMotion
                       DstMod=DstMod, &
                       DstMeshLoc=MeshLocType(ED_u_HubPtLoad), &        ! ED%Input(1)%HubPtLoad
                       DstDispMeshLoc=MeshLocType(ED_y_HubPtMotion), &  ! ED%y%HubPtMotion
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2)
      if (Failed()) return

   case (Module_AD)

      if (Turbine%p_FAST%CompElast == Module_ED) then
         do i = 1, size(Turbine%ED%Input(1)%BladePtLoads)
            call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                             SrcMeshLoc=MeshLocType(), &       ! AD%y%rotors(1)%BladeLoad(i)
                             SrcDispMeshLoc=MeshLocType(), &   ! AD%Input(1)%rotors(1)%BladeMotion(i)
                             DstMeshLoc=MeshLocType(), &       ! ED%Input(1)%BladePtLoads(i)
                             DstDispMeshLoc=MeshLocType(), &   ! ED%y%BladeLn2Mesh(i)
                             ErrStat=ErrStat2, ErrMsg=ErrMsg2)
            if (Failed()) return
         end do
      end if

      call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                       SrcMeshLoc=MeshLocType(AD_y_rotors_TowerLoad, AD_Rotor), &         ! AD%y%rotors(1)%TowerLoad
                       SrcDispMeshLoc=MeshLocType(AD_u_rotors_TowerMotion, AD_rotor), &   ! AD%u%rotors(1)%TowerMotion
                       DstMeshLoc=MeshLocType(ED_u_TowerPtLoads), &                       ! ED%Input(1)%TowerPtLoads
                       DstDispMeshLoc=MeshLocType(ED_y_TowerLn2Mesh), &                   ! ED%y%TowerLn2Mesh
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2, &
                       Active=Turbine%AD%y%rotors(1)%TowerLoad%committed)
      if (Failed()) return

      call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                       SrcMeshLoc=MeshLocType(AD_y_rotors_NacelleLoad, AD_Rotor), &          ! AD%y%rotors(1)%NacelleLoad
                       SrcDispMeshLoc=MeshLocType(AD_u_rotors_NacelleMotion, AD_rotor), &    ! AD%Input(1)%rotors(1)%NacelleMotion
                       DstMeshLoc=MeshLocType(ED_u_NacelleLoads), &                          ! ED%Input(1)%NacelleLoads
                       DstDispMeshLoc=MeshLocType(ED_y_NacelleMotion), &                     ! ED%y%NacelleMotion
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2, &
                       Active=Turbine%AD%Input(1)%rotors(1)%NacelleMotion%committed)
      if (Failed()) return

      call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                       SrcMeshLoc=MeshLocType(AD_y_rotors_HubLoad, AD_Rotor), &              ! AD%y%rotors(1)%HubLoad
                       SrcDispMeshLoc=MeshLocType(AD_u_rotors_HubMotion, AD_rotor), &        ! AD%Input(1)%rotors(1)%HubMotion
                       DstMeshLoc=MeshLocType(ED_u_HubPtLoad), &                             ! ED%Input(1)%HubPtLoad
                       DstDispMeshLoc=MeshLocType(ED_y_HubPtMotion), &                       ! ED%y%HubPtMotion
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2)
      if (Failed()) return

      call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                       SrcMeshLoc=MeshLocType(AD_y_rotors_TFinLoad, AD_Rotor), &             ! AD%y%rotors(1)%TFinLoad
                       SrcDispMeshLoc=MeshLocType(AD_u_rotors_TFinMotion, AD_rotor), &       ! AD%Input(1)%rotors(1)%TFinMotion
                       DstMeshLoc=MeshLocType(ED_u_TFinCMLoads), &                           ! ED%Input(1)%TFinCMLoads
                       DstDispMeshLoc=MeshLocType(ED_y_TFinCMMotion), &                      ! ED%y%TFinCMMotion
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2, &
                       Active=Turbine%AD%Input(1)%rotors(1)%TFinMotion%committed)
      if (Failed()) return

   case (Module_SD)

      call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                       SrcMeshLoc=MeshLocType(), &       ! SD%y%Y1mesh, &
                       SrcDispMeshLoc=MeshLocType(), &   ! SD%Input(1)%TPMesh
                       DstMeshLoc=MeshLocType(), &       ! ED%Input(1)%PlatformPtMesh
                       DstDispMeshLoc=MeshLocType(), &   ! ED%y%PlatformPtMesh
                       ErrStat=ErrStat2, ErrMsg=ErrMsg2)
      if (Failed()) return

   case (Module_SrvD)

      ! Nacelle Structural Controller
      do j = 1, Turbine%SrvD%p%NumNStC
         call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                          SrcMeshLoc=MeshLocType(SrvD_y_NStCLoadMesh, j), &         ! SrvD%y%NStCLoadMesh(j), &
                          SrcDispMeshLoc=MeshLocType(SrvD_u_NStCMotionMesh, j), &   ! SrvD%u%NStCMotionMesh(j)
                          DstMeshLoc=MeshLocType(ED_u_NacelleLoads), &              ! ED%Input(1)%NacelleLoads
                          DstDispMeshLoc=MeshLocType(ED_y_NacelleMotion), &         ! ED%y%NacelleMotion
                          ErrStat=ErrStat2, ErrMsg=ErrMsg2)
         if (Failed()) return
      end do

      ! Tower Structural Controller
      do j = 1, Turbine%SrvD%p%NumTStC
         call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                          SrcMeshLoc=MeshLocType(SrvD_y_TStCLoadMesh, j), &         ! SrvD%y%TStCLoadMesh(j), &
                          SrcDispMeshLoc=MeshLocType(SrvD_u_TStCMotionMesh, j), &   ! SrvD%u%TStCMotionMesh(j)
                          DstMeshLoc=MeshLocType(ED_u_TowerPtLoads), &              ! ED%Input(1)%TowerLoads
                          DstDispMeshLoc=MeshLocType(ED_y_TowerLn2Mesh), &          ! ED%y%TowerLn2Mesh
                          ErrStat=ErrStat2, ErrMsg=ErrMsg2)
         if (Failed()) return
      end do

      ! Blade Structural Controller (if ElastoDyn is used for blades)
      if (Turbine%p_FAST%CompElast == Module_ED) then
         do i = 1, Turbine%SrvD%p%NumBStC
            do j = 1, Turbine%ED%p%NumBl
               call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                                SrcMeshLoc=MeshLocType(SrvD_y_BStCLoadMesh, i, j), &        ! SrvD%y%BStCLoadMesh(i, j), &
                                SrcDispMeshLoc=MeshLocType(SrvD_u_BStCMotionMesh, i, j), &  ! SrvD%u%BStCMotionMesh(i, j)
                                DstMeshLoc=MeshLocType(ED_u_BladePtLoads, j), &             ! ED%Input(1)%BladePtLoads(j)
                                DstDispMeshLoc=MeshLocType(ED_y_BladeLn2Mesh, j), &         ! ED%y%BladeLn2Mesh(j)
                                ErrStat=ErrStat2, ErrMsg=ErrMsg2)
               if (Failed()) return
            end do
         end do
      end if

      ! Substructure Structural Controller (if not using SubDyn)
      if (Turbine%p_FAST%CompSub /= Module_SD) then
         do j = 1, Turbine%SrvD%p%NumSStC
            call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                             SrcMeshLoc=MeshLocType(SrvD_y_SStCLoadMesh, j), &         ! SrvD%y%SStCLoadMesh(j), &
                             SrcDispMeshLoc=MeshLocType(SrvD_u_SStCMotionMesh, j), &   ! SrvD%u%SStCMotionMesh(j)
                             DstMeshLoc=MeshLocType(ED_u_PlatformPtMesh), &            ! ED%Input(1)%PlatformPtMesh
                             DstDispMeshLoc=MeshLocType(ED_y_PlatformPtMesh), &        ! ED%y%PlatformPtMesh
                             ErrStat=ErrStat2, ErrMsg=ErrMsg2)
            if (Failed()) return
         end do
      end if

      call NonMeshMap(Mappings, "SrvD Data -> ED Data", SrcMod=SrcMod, DstMod=DstMod)  ! TODO

   end select

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine HD_InitInputMappings(Mappings, SrcMod, DstMod, Turbine, ErrStat, ErrMsg)
   type(TC_MappingType), allocatable      :: Mappings(:)
   type(ModDataType), intent(in)          :: SrcMod, DstMod
   type(FAST_TurbineType), intent(inout)  :: Turbine           !< Turbine type
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg

   character(*), parameter    :: RoutineName = 'HD_InitInputMappings'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2
   integer(IntKi)             :: i

   ErrStat = ErrID_None
   ErrMsg = ''

   select case (SrcMod%ID)
   case (Module_ED)

      call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                         SrcMeshLoc=MeshLocType(ED_y_PlatformPtMesh), & ! T%ED%y%PlatformPtMesh
                         DstMeshLoc=MeshLocType(HydroDyn_u_PRPMesh), &  ! T%HD%Input(1)%PRPMesh
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2)
      if (Failed()) return

      ! If SubDyn is not active substructure motion/loads come from ElastoDyn
      if (Turbine%p_FAST%CompSub /= Module_SD) then

         call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                            SrcMeshLoc=MeshLocType(ED_y_PlatformPtMesh), &    ! ED%y%PlatformPtMesh
                            DstMeshLoc=MeshLocType(HydroDyn_u_WAMITMesh), &   ! HD%Input(1)%WAMITMesh
                            ErrStat=ErrStat2, ErrMsg=ErrMsg2, &
                            Active=Turbine%HD%y%WAMITMesh%Committed)
         if (Failed()) return

         call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                            SrcMeshLoc=MeshLocType(ED_y_PlatformPtMesh), &       ! ED%y%PlatformPtMesh
                            DstMeshLoc=MeshLocType(HydroDyn_u_Morison_Mesh), &   ! HD%Input(1)%Morison%Mesh
                            ErrStat=ErrStat2, ErrMsg=ErrMsg2, &
                            Active=Turbine%HD%Input(1)%Morison%Mesh%Committed)
         if (Failed()) return
      end if

   case (Module_SD)
      call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                         SrcMeshLoc=MeshLocType(SD_y_Y2Mesh), &            ! SD%y%Y2Mesh
                         DstMeshLoc=MeshLocType(HydroDyn_u_WAMITMesh), &   ! HD%Input(1)%WAMITMesh
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2, &
                         Active=Turbine%HD%y%WAMITMesh%Committed)
      if (Failed()) return

      call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                         SrcMeshLoc=MeshLocType(SD_y_Y2Mesh), &               ! SD%y%Y2Mesh
                         DstMeshLoc=MeshLocType(HydroDyn_u_Morison_Mesh), &   ! HD%Input(1)%Morison%Mesh
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2, &
                         Active=Turbine%HD%Input(1)%Morison%Mesh%Committed)
      if (Failed()) return

   end select

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine IfW_InitInputMappings(Mappings, SrcMod, DstMod, T, ErrStat, ErrMsg)
   type(TC_MappingType), allocatable      :: Mappings(:)
   type(ModDataType), intent(in)          :: SrcMod, DstMod
   type(FAST_TurbineType), intent(inout)  :: T           !< Turbine type
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg

   character(*), parameter    :: RoutineName = 'BD_InitInputMappings'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2
   integer(IntKi)             :: i

   ErrStat = ErrID_None
   ErrMsg = ''

   select case (SrcMod%ID)
   case (Module_ED)
      call NonMeshMap(Mappings, "ED HubMotion -> IfW HubMotion", SrcMod=SrcMod, DstMod=DstMod)
   end select

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine SD_InitInputMappings(Mappings, SrcMod, DstMod, Turbine, ErrStat, ErrMsg)
   type(TC_MappingType), allocatable      :: Mappings(:)
   type(ModDataType), intent(in)          :: SrcMod, DstMod
   type(FAST_TurbineType), intent(inout)  :: Turbine           !< Turbine type
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg

   character(*), parameter    :: RoutineName = 'BD_InitInputMappings'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2
   integer(IntKi)             :: i, j

   ErrStat = ErrID_None
   ErrMsg = ''

   select case (SrcMod%ID)

   case (Module_ED)
      call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                         SrcMeshLoc=MeshLocType(ED_y_PlatformPtMesh), &    ! T%ED%y%PlatformPtMesh
                         DstMeshLoc=MeshLocType(SD_u_TPMesh), &            ! T%SD%Input(1)%TPMesh
                         ErrStat=ErrStat2, ErrMsg=ErrMsg2)
      if (Failed()) return

   case (Module_SrvD)

      ! Substructure Structural Controller
      do j = 1, Turbine%SrvD%p%NumSStC
         call MapLoadMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                          SrcMeshLoc=MeshLocType(SrvD_y_SStCLoadMesh, j), &         ! SrvD%y%SStCLoadMesh(j), &
                          SrcDispMeshLoc=MeshLocType(SrvD_u_SStCMotionMesh, j), &   ! SrvD%u%SStCMotionMesh(j)
                          DstMeshLoc=MeshLocType(SD_u_LMesh), &                     ! SD%Input(1)%LMesh
                          DstDispMeshLoc=MeshLocType(SD_y_y3Mesh), &                ! SD%y%y3Mesh
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

subroutine SrvD_InitInputMappings(Mappings, SrcMod, DstMod, Turbine, ErrStat, ErrMsg)
   type(TC_MappingType), allocatable      :: Mappings(:)
   type(ModDataType), intent(in)          :: SrcMod, DstMod
   type(FAST_TurbineType), intent(inout)  :: Turbine           !< Turbine type
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg

   character(*), parameter    :: RoutineName = 'BD_InitInputMappings'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2
   integer(IntKi)             :: i, j

   ErrStat = ErrID_None
   ErrMsg = ''

   select case (SrcMod%ID)
   case (Module_BD)
      call NonMeshMap(Mappings, "BD Data -> SrvD Data", SrcMod=SrcMod, DstMod=DstMod)
      call NonMeshMap(Mappings, "BD RootM -> SrvD RootM", SrcMod=SrcMod, DstMod=DstMod)

   case (Module_ED)

      call NonMeshMap(Mappings, "ED Data -> SrvD Data", SrcMod=SrcMod, DstMod=DstMod)
      if (Turbine%p_FAST%CompElast == Module_ED) then
         call NonMeshMap(Mappings, "ED RootM -> SrvD RootM", SrcMod=SrcMod, DstMod=DstMod)
      end if

      ! Nacelle Structural Controller
      do j = 1, Turbine%SrvD%p%NumNStC
         call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                            SrcMeshLoc=MeshLocType(ED_y_NacelleMotion), &        ! ED%y%NacelleMotion
                            DstMeshLoc=MeshLocType(SrvD_u_NStCMotionMesh, j), &  ! SrvD%u%NStCMotionMesh(j)
                            ErrStat=ErrStat2, ErrMsg=ErrMsg2)
         if (Failed()) return
      end do

      ! Tower Structural Controller
      do j = 1, Turbine%SrvD%p%NumTStC
         call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                            SrcMeshLoc=MeshLocType(ED_y_TowerLn2Mesh), &         ! ED%y%TowerMotion
                            DstMeshLoc=MeshLocType(SrvD_u_TStCMotionMesh, j), &  ! SrvD%u%TStCMotionMesh(j)
                            ErrStat=ErrStat2, ErrMsg=ErrMsg2)
         if (Failed()) return
      end do

      ! Blade Structural Controller (if ElastoDyn is used for blades)
      if (Turbine%p_FAST%CompElast == Module_ED) then
         do i = 1, Turbine%SrvD%p%NumBStC
            do j = 1, Turbine%ED%p%NumBl
               call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                                  SrcMeshLoc=MeshLocType(ED_y_BladeLn2Mesh, j), &         ! ED%y%BladeLn2Mesh(j)
                                  DstMeshLoc=MeshLocType(SrvD_u_BStCMotionMesh, i, j), &  ! SrvD%u%BStCMotionMesh(i, j)
                                  ErrStat=ErrStat2, ErrMsg=ErrMsg2)
               if (Failed()) return
            end do
         end do
      end if

      ! Substructure Structural Controller (if not using SubDyn)
      if (Turbine%p_FAST%CompSub /= Module_SD) then
         do j = 1, Turbine%SrvD%p%NumSStC
            call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                               SrcMeshLoc=MeshLocType(ED_y_PlatformPtMesh), &       ! ED%y%PlatformPtMesh
                               DstMeshLoc=MeshLocType(SrvD_u_SStCMotionMesh, j), &  ! SrvD%u%SStCMotionMesh(j)
                               ErrStat=ErrStat2, ErrMsg=ErrMsg2)
            if (Failed()) return
         end do
      end if

   case (Module_IfW)
      call NonMeshMap(Mappings, "IfW Data -> SrvD Data", SrcMod=SrcMod, DstMod=DstMod)

   case (Module_SD)

      ! Substructure Structural Controller
      do j = 1, Turbine%SrvD%p%NumSStC
         call MapMotionMesh(Turbine, Mappings, SrcMod=SrcMod, DstMod=DstMod, &
                            SrcMeshLoc=MeshLocType(SD_y_y3Mesh), &               ! SD%y%y3Mesh
                            DstMeshLoc=MeshLocType(SrvD_u_SStCMotionMesh, j), &  ! SrvD%u%SStCMotionMesh(j)
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

   ErrStat = ErrID_None
   ErrMsg = ''

   ! If active argument is set to false, return
   if (present(Active)) then
      if (.not. Active) return
   end if

   ! Get mesh pointers
   SrcMesh => FAST_OutputMeshPointer(SrcMod, Turbine, SrcMeshLoc)
   SrcDispMesh => FAST_InputMeshPointer(SrcMod, Turbine, SrcDispMeshLoc, UseU=.false.)
   DstMesh => FAST_InputMeshPointer(DstMod, Turbine, DstMeshLoc, UseU=.false.)
   DstDispMesh => FAST_OutputMeshPointer(DstMod, Turbine, DstDispMeshLoc)

   ! Check that all meshes in mapping have nonzero identifiers
   if (SrcMesh%ID == 0) then
      call SetErrStat(ErrID_Fatal, 'SrcMesh not in module variable', ErrStat, ErrMsg, RoutineName)
      return
   else if (SrcDispMesh%ID == 0) then
      call SetErrStat(ErrID_Fatal, 'SrcDispMesh not in module variable', ErrStat, ErrMsg, RoutineName)
      return
   else if (DstMesh%ID == 0) then
      call SetErrStat(ErrID_Fatal, 'DstMesh not in module variable', ErrStat, ErrMsg, RoutineName)
      return
   else if (DstDispMesh%ID == 0) then
      call SetErrStat(ErrID_Fatal, 'DstDispMesh not in module variable', ErrStat, ErrMsg, RoutineName)
      return
   end if

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
   call MeshCopy(DstMesh, Mapping%MeshTmp, MESH_NEWCOPY, ErrStat2, ErrMsg2); if (Failed()) return

   ! Get mapping indices for linearized mesh mapping
   call InitMeshLinearization(Mapping, SrcMod, DstMod, SrcMesh, DstMesh, SrcDispMesh, DstDispMesh)

   ! Add mapping to array of mappings
   Mappings = [Mappings, Mapping]

contains
   logical function Failed()
      Failed = ErrStat2 >= AbortErrLev
      if (Failed) call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
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
   SrcMesh => FAST_OutputMeshPointer(SrcMod, Turbine, SrcMeshLoc)
   DstMesh => FAST_InputMeshPointer(DstMod, Turbine, DstMeshLoc, UseU=.false.)

   ! Check that all meshes in mapping have nonzero identifiers
   if (SrcMesh%ID == 0) then
      call SetErrStat(ErrID_Fatal, 'SrcMesh not in module variable', ErrStat, ErrMsg, RoutineName)
      return
   else if (DstMesh%ID == 0) then
      call SetErrStat(ErrID_Fatal, 'DstMesh not in module variable', ErrStat, ErrMsg, RoutineName)
      return
   end if

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
   call InitMeshLinearization(Mapping, SrcMod, DstMod, SrcMesh, DstMesh)

   ! Add mapping to array of mappings
   Mappings = [Mappings, Mapping]

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine NonMeshMap(Maps, Key, SrcMod, DstMod, i1, i2, Active)
   type(TC_MappingType), allocatable      :: Maps(:)
   character(*), intent(in)               :: Key
   type(ModDataType), intent(in)          :: SrcMod, DstMod
   integer(IntKi), optional, intent(in)   :: i1, i2
   logical, optional, intent(in)          :: Active
   type(TC_MappingType)                   :: Mapping

   if (present(Active)) then
      if (.not. Active) return
   end if

   ! Initialize mapping structure
   Mapping%MapType = Map_NonMesh
   Mapping%SrcModIdx = SrcMod%Idx
   Mapping%SrcModID = SrcMod%ID
   Mapping%SrcIns = SrcMod%Ins
   Mapping%DstModIdx = DstMod%Idx
   Mapping%DstModID = DstMod%ID
   Mapping%DstIns = DstMod%Ins

   ! Get optional mapping indicies
   if (present(i1)) Mapping%i1 = i1
   if (present(i2)) Mapping%i2 = i2

   Maps = [Maps, Mapping]
end subroutine

subroutine InitMeshLinearization(Mapping, SrcMod, DstMod, SrcMesh, DstMesh, SrcDispMesh, DstDispMesh)
   type(TC_MappingType), intent(inout)    :: Mapping
   type(ModDataType), intent(in)          :: SrcMod, DstMod
   type(MeshType), intent(in)             :: SrcMesh, DstMesh
   type(MeshType), optional, intent(in)   :: SrcDispMesh, DstDispMesh

   ! Save source and destination mesh ID
   Mapping%SrcMeshID = SrcMesh%ID
   Mapping%DstMeshID = DstMesh%ID

   ! Determine transfer type
   if ((SrcMesh%ElemTable(ELEMENT_POINT)%nelem > 0) .and. (DstMesh%ElemTable(ELEMENT_POINT)%nelem > 0)) then
      Mapping%XfrType = Xfr_Point_to_Point
   else if ((SrcMesh%ElemTable(ELEMENT_POINT)%nelem > 0) .and. (DstMesh%ElemTable(ELEMENT_LINE2)%nelem > 0)) then
      Mapping%XfrType = Xfr_Point_to_Line2
   else if ((SrcMesh%ElemTable(ELEMENT_LINE2)%nelem > 0) .and. (DstMesh%ElemTable(ELEMENT_POINT)%nelem > 0)) then
      Mapping%XfrType = Xfr_Line2_to_Point
   else if ((SrcMesh%ElemTable(ELEMENT_LINE2)%nelem > 0) .and. (DstMesh%ElemTable(ELEMENT_LINE2)%nelem > 0)) then
      Mapping%XfrType = Xfr_Line2_to_Line2
   end if

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

subroutine FAST_LinearizeMappings(Turbine, Mods, Mappings, ModOrder, ErrStat, ErrMsg, dUdu, dUdy)
   type(FAST_TurbineType), target, intent(inout)   :: Turbine     !< Turbine type
   type(ModDataType), intent(in)                   :: Mods(:)     !< Module data
   type(TC_MappingType), intent(inout)             :: Mappings(:)
   integer(IntKi), intent(in)                      :: ModOrder(:)
   integer(IntKi), intent(out)                     :: ErrStat
   character(*), intent(out)                       :: ErrMsg
   real(R8Ki), optional, intent(inout)             :: dUdu(:, :), dUdy(:, :)

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
            case (Map_NonMesh)
               cycle

            case (Map_MotionMesh)

               ! Get source and destination meshes
               SrcMesh => FAST_OutputMeshPointer(Mods(Mapping%SrcModIdx), Turbine, Mapping%SrcMeshLoc)
               DstMesh => FAST_InputMeshPointer(Mods(Mapping%DstModIdx), Turbine, Mapping%DstMeshLoc, UseU=.false.)

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

            case (Map_LoadMesh)

               ! Get source and destination meshes
               SrcMesh => FAST_OutputMeshPointer(Mods(Mapping%SrcModIdx), Turbine, Mapping%SrcMeshLoc)
               DstMesh => FAST_InputMeshPointer(Mods(Mapping%DstModIdx), Turbine, Mapping%DstMeshLoc, UseU=.false.)

               ! Get source and destination displacement meshes
               SrcDispMesh => FAST_InputMeshPointer(Mods(Mapping%SrcModIdx), Turbine, Mapping%SrcDispMeshLoc, UseU=.false.)
               DstDispMesh => FAST_OutputMeshPointer(Mods(Mapping%DstModIdx), Turbine, Mapping%DstDispMeshLoc)

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

            end select

            write (*, *) trim(FAST_OutputMeshName(Mods(Mapping%SrcModIdx), Turbine, Mapping%SrcMeshLoc)), " -> ", &
               FAST_InputMeshName(Mods(Mapping%DstModIdx), Turbine, Mapping%DstMeshLoc)

            ! Check for errors
            call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            if (ErrStat >= AbortErrLev) return

            ! Copy linearization matrices to global dUdu matrix
            if (present(dUdu)) then
               call dUduSetBlocks(Mapping, Mapping%MeshMap%dM)
            end if

            ! Copy linearization matrices to global dUdy matrix
            if (present(dUdy)) then
               call dUdySetBlocks(Mapping, Mapping%MeshMap%dM)
            end if

         end associate
      end do
   end do

contains
   subroutine dUduSetBlocks(Mapping, dM)
      type(TC_MappingType), intent(inout)          :: Mapping  !< Mapping
      type(MeshMapLinearizationType), intent(in)   :: dM       !< Mesh Map Linearization data

      ! Effect of input Translation Velocity on input Translation Displacement
      if (allocated(dM%tv_uD)) then
         ! call SetBlock(DstMod%Vars%u(M%DstVarIdx), VF_TransVel, DstMod%Vars%u(M%DstVarIdx), VF_TransDisp, -MML%tv_uD, dUdu)
         call SetBlock(Mapping%iLocDstTransVel, Mapping%iLocDstTransDisp, -dM%tv_uD, dUdU)
      end if

      ! Effect of input Translation Acceleration on input Translation Displacement
      if (allocated(dM%ta_uD)) then
         ! call SetBlock(DstMod%Vars%u(M%DstVarIdx), VF_TransAcc, DstMod%Vars%u(M%DstVarIdx), VF_TransDisp, -MML%ta_uD, dUdu)
         call SetBlock(Mapping%iLocDstTransAcc, Mapping%iLocDstTransDisp, -dM%ta_uD, dUdU)
      end if

      ! Effect of input Moments on input Translation Displacement
      if (allocated(dM%M_uS)) then
         ! call SetBlock(DstMod%Vars%u(M%DstVarIdx), VF_Moment, SrcMod%Vars%u([M%SrcDispVarIdx]), VF_TransDisp, -MML%M_uS, dUdu)
         call SetBlock(Mapping%iLocDstMoment, Mapping%iLocSrcDispTransDisp, -dM%M_uS, dUdU)
      end if
   end subroutine

   subroutine dUdySetBlocks(Mapping, dM)
      type(TC_MappingType), intent(inout)          :: Mapping     !< Mapping
      type(MeshMapLinearizationType), intent(in)   :: dM          !< Mesh Map Linearization data

      ! Load identity
      if (allocated(dM%li)) then
         ! call SetBlock(DstMod%Vars%u(M%DstVarIdx), VF_Force, SrcMod%Vars%y(M%SrcVarIdx), VF_Force, -MML%li, dUdy)
         call SetBlock(Mapping%iLocDstForce, Mapping%iLocSrcForce, -dM%li, dUdy)
         ! call SetBlock(DstMod%Vars%u(M%DstVarIdx), VF_Moment, SrcMod%Vars%y(M%SrcVarIdx), VF_Moment, -MML%li, dUdy)
         call SetBlock(Mapping%iLocDstMoment, Mapping%iLocSrcMoment, -dM%li, dUdy)
      end if

      ! Moment to Force
      if (allocated(dM%m_f)) then
         ! call SetBlock(DstMod%Vars%u(M%DstVarIdx), VF_Moment, SrcMod%Vars%y(M%SrcVarIdx), VF_Force, -MML%m_f, dUdy)
         call SetBlock(Mapping%iLocDstMoment, Mapping%iLocSrcForce, -dM%m_f, dUdy)
      end if

      ! Moment to destination translation displacement
      if (allocated(dM%m_uD)) then
         ! call SetBlock(DstMod%Vars%u(M%DstVarIdx), VF_Moment, DstMod%Vars%y([M%DstDispVarIdx]), VF_TransDisp, -MML%m_uD, dUdy)
         call SetBlock(Mapping%iLocDstMoment, Mapping%iLocDstDispTransDisp, -dM%m_uD, dUdy)
      end if

      ! Motion identity
      if (allocated(dM%mi)) then
         ! call SetBlock(DstMod%Vars%u(M%DstVarIdx), VF_TransDisp, SrcMod%Vars%y(M%SrcVarIdx), VF_TransDisp, -MML%mi, dUdy)
         call SetBlock(Mapping%iLocDstTransDisp, Mapping%iLocSrcTransDisp, -dM%mi, dUdy)
         ! call SetBlock(DstMod%Vars%u(M%DstVarIdx), VF_Orientation, SrcMod%Vars%y(M%SrcVarIdx), VF_Orientation, -MML%mi, dUdy)
         call SetBlock(Mapping%iLocDstOrientation, Mapping%iLocSrcOrientation, -dM%mi, dUdy)
         ! call SetBlock(DstMod%Vars%u(M%DstVarIdx), VF_TransVel, SrcMod%Vars%y(M%SrcVarIdx), VF_TransVel, -MML%mi, dUdy)
         call SetBlock(Mapping%iLocDstTransVel, Mapping%iLocSrcTransVel, -dM%mi, dUdy)
         ! call SetBlock(DstMod%Vars%u(M%DstVarIdx), VF_AngularVel, SrcMod%Vars%y(M%SrcVarIdx), VF_AngularVel, -MML%mi, dUdy)
         call SetBlock(Mapping%iLocDstAngularVel, Mapping%iLocSrcAngularVel, -dM%mi, dUdy)
         ! call SetBlock(DstMod%Vars%u(M%DstVarIdx), VF_TransAcc, SrcMod%Vars%y(M%SrcVarIdx), VF_TransAcc, -MML%mi, dUdy)
         call SetBlock(Mapping%iLocDstTransAcc, Mapping%iLocSrcTransAcc, -dM%mi, dUdy)
         ! call SetBlock(DstMod%Vars%u(M%DstVarIdx), VF_AngularAcc, SrcMod%Vars%y(M%SrcVarIdx), VF_AngularAcc, -MML%mi, dUdy)
         call SetBlock(Mapping%iLocDstAngularAcc, Mapping%iLocSrcAngularAcc, -dM%mi, dUdy)
      end if

      ! Translation to Rotation
      if (allocated(dM%fx_p)) then
         ! call SetBlock(DstMod%Vars%u(M%DstVarIdx), VF_TransDisp, SrcMod%Vars%y(M%SrcVarIdx), VF_Orientation, -MML%fx_p, dUdy)
         call SetBlock(Mapping%iLocDstTransDisp, Mapping%iLocSrcOrientation, -dM%fx_p, dUdy)
         ! call SetBlock(DstMod%Vars%u(M%DstVarIdx), VF_TransVel, SrcMod%Vars%y(M%SrcVarIdx), VF_AngularVel, -MML%fx_p, dUdy)
         call SetBlock(Mapping%iLocDstTransVel, Mapping%iLocSrcAngularVel, -dM%fx_p, dUdy)
         ! call SetBlock(DstMod%Vars%u(M%DstVarIdx), VF_TransAcc, SrcMod%Vars%y(M%SrcVarIdx), VF_AngularAcc, -MML%fx_p, dUdy)
         call SetBlock(Mapping%iLocDstTransAcc, Mapping%iLocSrcAngularAcc, -dM%fx_p, dUdy)
      end if

      ! Translation velocity to translation displacement
      if (allocated(dM%tv_us)) then
         ! call SetBlock(DstMod%Vars%u(M%DstVarIdx), VF_TransVel, SrcMod%Vars%y(M%SrcVarIdx), VF_TransDisp, -MML%tv_us, dUdy)
         call SetBlock(Mapping%iLocDstTransVel, Mapping%iLocDstDispTransDisp, -dM%tv_us, dUdy)
      end if

      ! Translation acceleration to translation displacement
      if (allocated(dM%ta_us)) then
         ! call SetBlock(DstMod%Vars%u(M%DstVarIdx), VF_TransAcc, SrcMod%Vars%y(M%SrcVarIdx), VF_TransDisp, -MML%ta_us, dUdy)
         call SetBlock(Mapping%iLocDstTransAcc, Mapping%iLocSrcTransDisp, -dM%ta_us, dUdy)
      end if

      ! Translation acceleration to angular velocity
      if (allocated(dM%ta_rv)) then
         ! call SetBlock(DstMod%Vars%u(M%DstVarIdx), VF_TransAcc, SrcMod%Vars%y(M%SrcVarIdx), VF_AngularVel, -MML%ta_rv, dUdy)
         call SetBlock(Mapping%iLocDstTransAcc, Mapping%iLocSrcAngularVel, -dM%ta_rv, dUdy)
      end if
   end subroutine

   subroutine SetBlock(iLocRow, iLocCol, SrcM, DstM)
      integer(IntKi), intent(in)    :: iLocRow(2), iLocCol(2)
      real(R8Ki), intent(in)        :: SrcM(:, :)
      real(R8Ki), intent(inout)     :: DstM(:, :)
      if (iLocRow(1) > 0 .and. iLocCol(1) > 0) then
         associate (DstSubM => DstM(iLocRow(1):iLocRow(1)+size(SrcM,1)-1, iLocCol(1):iLocCol(1)+size(SrcM,2)-1))
         ! associate (DstSubM => DstM(iLocRow(1):iLocRow(2), iLocCol(1):iLocCol(2)))
            ! if ((size(SrcM, 1) /= (iLocRow(2) - iLocRow(1) + 1)) .or. (size(SrcM, 2) /= (iLocCol(2) - iLocCol(1)) + 1)) then
            !    print *, "hello"
            ! end if
            DstSubM = DstSubM + SrcM
         end associate
      end if
   end subroutine
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
         case (Map_NonMesh)
            call NonMesh_InputSolve(Turbine, Mapping, ErrStat2, ErrMsg2, UseU)
            if (Failed()) return

         case (Map_MotionMesh)

            ! Get source and destination meshes
            SrcMesh => FAST_OutputMeshPointer(Mods(Mapping%SrcModIdx), Turbine, Mapping%SrcMeshLoc)
            DstMesh => FAST_InputMeshPointer(Mods(Mapping%DstModIdx), Turbine, Mapping%DstMeshLoc, UseU)

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
            SrcMesh => FAST_OutputMeshPointer(Mods(Mapping%SrcModIdx), Turbine, Mapping%SrcMeshLoc)
            DstMesh => FAST_InputMeshPointer(Mods(Mapping%DstModIdx), Turbine, Mapping%DstMeshLoc, UseU)

            ! Get source and destination displacement meshes
            SrcDispMesh => FAST_InputMeshPointer(Mods(Mapping%SrcModIdx), Turbine, Mapping%SrcDispMeshLoc, UseU)
            DstDispMesh => FAST_OutputMeshPointer(Mods(Mapping%DstModIdx), Turbine, Mapping%DstDispMeshLoc)

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

subroutine NonMesh_InputSolve(Turbine, Mapping, ErrStat, ErrMsg, UseU)
   type(FAST_TurbineType), intent(inout)  :: Turbine     !< Turbine type
   type(TC_MappingType), intent(in)       :: Mapping
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg
   logical, intent(in)                    :: UseU        ! Flag to transfer to u instead of Input

   character(*), parameter       :: RoutineName = 'NonMesh_InputSolve'
   integer(IntKi)                :: ErrStat2
   character(ErrMsgLen)          :: ErrMsg2
   integer(IntKi)                :: i, j, k

   ErrStat = ErrID_None
   ErrMsg = ''

! case ("BD RootM -> SrvD RootM")

   !          u_SrvD%RootMxc(Maps(i)%SrcIns) = T%BD%y(Maps(i)%SrcIns)%RootMxr*cos(T%ED%y%BlPitch(Maps(i)%SrcIns)) + &
   !                                           T%BD%y(Maps(i)%SrcIns)%RootMyr*sin(T%ED%y%BlPitch(Maps(i)%SrcIns))
   !          u_SrvD%RootMyc(Maps(i)%SrcIns) = -T%BD%y(Maps(i)%SrcIns)%RootMxr*sin(T%ED%y%BlPitch(Maps(i)%SrcIns)) + &
   !                                           T%BD%y(Maps(i)%SrcIns)%RootMyr*cos(T%ED%y%BlPitch(Maps(i)%SrcIns))

   !       case ("ED RootM -> SrvD RootM")

   !          u_SrvD%RootMxc = T%ED%y%RootMxc ! fixed-size arrays: always size 3
   !          u_SrvD%RootMyc = T%ED%y%RootMyc ! fixed-size arrays: always size 3

   !       case ("ED Data -> SrvD Data")

   !          u_SrvD%YawAngle = T%ED%y%YawAngle ! nacelle yaw plus platform yaw

   !          u_SrvD%Yaw = T%ED%y%Yaw  ! nacelle yaw
   !          u_SrvD%YawRate = T%ED%y%YawRate
   !          u_SrvD%BlPitch = T%ED%y%BlPitch
   !          u_SrvD%LSS_Spd = T%ED%y%LSS_Spd
   !          u_SrvD%HSS_Spd = T%ED%y%HSS_Spd
   !          u_SrvD%RotSpeed = T%ED%y%RotSpeed

   !          u_SrvD%YawBrTAxp = T%ED%y%YawBrTAxp
   !          u_SrvD%YawBrTAyp = T%ED%y%YawBrTAyp
   !          u_SrvD%LSSTipPxa = T%ED%y%LSSTipPxa

   !          u_SrvD%LSSTipMxa = T%ED%y%LSSTipMxa
   !          u_SrvD%LSSTipMya = T%ED%y%LSSTipMya
   !          u_SrvD%LSSTipMza = T%ED%y%LSSTipMza
   !          u_SrvD%LSSTipMys = T%ED%y%LSSTipMys
   !          u_SrvD%LSSTipMzs = T%ED%y%LSSTipMzs

   !          u_SrvD%YawBrMyn = T%ED%y%YawBrMyn
   !          u_SrvD%YawBrMzn = T%ED%y%YawBrMzn
   !          u_SrvD%NcIMURAxs = T%ED%y%NcIMURAxs
   !          u_SrvD%NcIMURAys = T%ED%y%NcIMURAys
   !          u_SrvD%NcIMURAzs = T%ED%y%NcIMURAzs

   !          u_SrvD%RotPwr = T%ED%y%RotPwr

   !          u_SrvD%LSShftFxa = T%ED%y%LSShftFxa
   !          u_SrvD%LSShftFys = T%ED%y%LSShftFys
   !          u_SrvD%LSShftFzs = T%ED%y%LSShftFzs

   !       case ('ED PlatformMotion -> SrvD PlatformMotion')
   !       case ('ED TowerMotion -> SrvD TowerMotion')
   !       case ('ED NacelleMotion -> SrvD NacelleMotion')
   !       case ('ED BladeMotion -> SrvD BladeMotion')

   !       case ("IfW Data -> SrvD Data")

   !          u_SrvD%WindDir = atan2(T%IfW%y%VelocityUVW(2, 1), T%IfW%y%VelocityUVW(1, 1))
   !          u_SrvD%HorWindV = sqrt(T%IfW%y%VelocityUVW(1, 1)**2 + T%IfW%y%VelocityUVW(2, 1)**2)
   !          if (allocated(T%IfW%y%lidar%LidSpeed)) u_SrvD%LidSpeed = T%IfW%y%lidar%LidSpeed
   !          if (allocated(T%IfW%y%lidar%MsrPositionsX)) u_SrvD%MsrPositionsX = T%IfW%y%lidar%MsrPositionsX
   !          if (allocated(T%IfW%y%lidar%MsrPositionsY)) u_SrvD%MsrPositionsY = T%IfW%y%lidar%MsrPositionsY
   !          if (allocated(T%IfW%y%lidar%MsrPositionsZ)) u_SrvD%MsrPositionsZ = T%IfW%y%lidar%MsrPositionsZ

!    ! Zero tower and platform added mass
!    ! u_ED%TwrAddedMass = 0.0_ReKi
!    ! u_ED%PtfmAddedMass = 0.0_ReKi

!       case ("SrvD Data -> ED Data")
!          if (Linearize) then
!          else
!             u_ED%GenTrq = T%SrvD%y%GenTrq
!             u_ED%HSSBrTrqC = T%SrvD%y%HSSBrTrqC
!             u_ED%BlPitchCom = T%SrvD%y%BlPitchCom
!             u_ED%YawMom = T%SrvD%y%YawMom
!          end if

!       case ('SrvD BlAirfoilCom -> AD UserProp')
!          ! Set Conrol parameter (i.e. flaps) if using ServoDyn bem:
!          ! This takes in flap deflection for each blade (only one flap deflection angle per blade),
!          ! from ServoDyn (which comes from Bladed style DLL controller)
!          ! Commanded Airfoil UserProp for blade (must be same units as given in AD15 airfoil tables)
!          ! This is passed to AD15 to be interpolated with the airfoil table userprop column
!          ! (might be used for airfoil flap angles for example)
!          ! Must be same units as given in airfoil (no unit conversions handled in code)ß
!          ! do k_bl = 1, size(u_AD%rotors(1)%UserProp, dim=2)
!          !    do k_bn = 1, size(u_AD%rotors(1)%UserProp, dim=1)
!          !       u_AD%rotors(1)%UserProp(k_bn, k_bl) = T%SrvD%y%BlAirfoilCom(k_bl)
!          !    end do
!          ! end do

!       case ('ED HubMotion -> IfW HubMotion')

!          u_IfW%PositionXYZ(:, 1) = T%ED%y%HubPtMotion%Position(:, 1)
!          u_IfW%HubPosition = T%ED%y%HubPtMotion%Position(:, 1) + &
!                              T%ED%y%HubPtMotion%TranslationDisp(:, 1)
!          u_IfW%HubOrientation = T%ED%y%HubPtMotion%Orientation(:, :, 1)

!          ! Set Lidar position directly from hub motion mesh
!          u_IfW%lidar%HubDisplacementX = T%ED%y%HubPtMotion%TranslationDisp(1, 1)
!          u_IfW%lidar%HubDisplacementY = T%ED%y%HubPtMotion%TranslationDisp(2, 1)
!          u_IfW%lidar%HubDisplacementZ = T%ED%y%HubPtMotion%TranslationDisp(3, 1)

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
      if (associated(Maps(i)%MeshTmp%RemapFlag)) Maps(i)%MeshTmp%RemapFlag = .false.
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
