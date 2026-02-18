!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2025  National Renewable Energy Laboratory
!
!    This file is a module specific to an experimental wave tank at NREL.
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
!
!  This module provides structural model for the wavetank interface
!
!**********************************************************************************************************************************
module WaveTank_Struct
   use ISO_C_BINDING
   use NWTC_Library
   use WaveTank_Types

   implicit none
   private

   save

   public :: StructCreate
   public :: StructCreateMeshMaps
   public :: StructDestroy
   public :: StructMotionUpdate
   public :: StructLoadsMeshTransfer
   public :: WrVTK_Struct_Ref
   public :: WrVTK_Struct
   public :: FroudeScaleM2F_Disp
   public :: FroudeScaleM2F_TVel
   public :: FroudeScaleM2F_RVel
   public :: FroudeScaleM2F_TAcc
   public :: FroudeScaleM2F_RAcc
   public :: FroudeScaleM2F_Time
   public :: FroudeScaleF2M_Frc
   public :: FroudeScaleF2M_Mom

contains


!> create the structural model, allocate temp data storage, setup mesh mappings
subroutine StructCreate(SimSettings, MeshMotions, MeshLoads, MeshMaps, StructTmp, ErrStat, ErrMsg)
   type(SimSettingsType),  target,  intent(in   )  :: SimSettings
   type(MeshesMotionType), target,  intent(inout)  :: MeshMotions
   type(MeshesLoadsType ), target,  intent(inout)  :: MeshLoads
   type(MeshesMapsType  ),          intent(inout)  :: MeshMaps
   type(StructTmpType   ),          intent(inout)  :: StructTmp
   integer(IntKi),                  intent(  out)  :: ErrStat
   character(ErrMsgLen),            intent(  out)  :: ErrMsg
   integer(IntKi)                                  :: ErrStat2
   character(ErrMsgLen)                            :: ErrMsg2
   character(*),           parameter               :: RoutineName = 'WaveTank::StructCreate'
   real(ReKi)                                      :: TmpPos(3)
   real(DbKi)                                      :: AzBlade        ! temporary var for calculating blade mounting azimuth
   real(DbKi)                                      :: TmpAng(3)      ! temporary euler angle
   real(DbKi)                                      :: Orient(3,3)    ! temporary orientation
   type(TurbConfigType),   pointer                 :: TrbCfg         ! to shorten notation
   type(TurbInitCondType), pointer                 :: TrbInit        ! to shorten notation
   type(MeshType),         pointer                 :: Ptfm, PtfmLd   ! to shorten notation
   type(MeshType),         pointer                 :: Twr,  TwrLd    ! to shorten notation
   type(MeshType),         pointer                 :: Hub,  HubLd    ! to shorten notation
   type(MeshType),         pointer                 :: Root, RootLd   ! to shorten notation
   type(MeshType),         pointer                 :: MoorLd         ! to shorten notation
   integer(IntKi)                                  :: k              ! blade counter
   ErrStat = ErrID_None
   ErrMsg  = ''

   TrbCfg    => SimSettings%TrbCfg
   TrbInit   => SimSettings%TrbInit

   ! Set some state information
   StructTmp%RotSpeed = TrbInit%RotSpeed
   StructTmp%BldPitch = TrbInit%BldPitch
   StructTmp%NacYaw   = TrbInit%NacYaw
   StructTmp%Azimuth  = TrbInit%Azimuth


   !-------------------------------
   ! Wave measurement buoy
   !-------------------------------
   TmpPos = 0.0_ReKi
   TmpPos(1:2) = SimSettings%WaveBuoy%XYLoc(1:2)
   call Eye(Orient, ErrStat2, ErrMsg2);   if (Failed()) return
   call CreateInputPointMesh(MeshMotions%WaveBuoyMotion, TmpPos, Orient, ErrStat2, ErrMsg2, hasMotion=.true., hasLoads=.false.); if (Failed()) return


   !-------------------------------
   ! create PRP platform mesh point
   !-------------------------------
   Ptfm => MeshMotions%PtfmPtMotion
   TmpPos = real(TrbCfg%PtfmRefPos, ReKi)
   Orient=WT_EulerToDCM_fromInput(TrbCfg%PtfmRefOrient)
   call CreateInputPointMesh(Ptfm, TmpPos, Orient, ErrStat2, ErrMsg2, hasMotion=.true., hasLoads=.false.); if (Failed()) return
   Ptfm%RemapFlag = .false.

   ! create platform load mesh
   PtfmLd => MeshLoads%PtfmPtLoads
   call MeshCopy( SrcMesh=Ptfm, DestMesh=PtfmLd, CtrlCode=MESH_SIBLING, IOS=COMPONENT_OUTPUT, ErrStat=ErrStat2, ErrMess=ErrMsg2, Force=.true., Moment=.true. )
   if (Failed()) return
   PtfmLd%RemapFlag = .false.

   ! create a temporary load mesh
   call MeshCopy( PtfmLd, MeshLoads%PtfmPtLoadsTmp, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
   if (Failed()) return
   PtfmLd%RemapFlag = .false.

   ! create a mooring mesh point
   MoorLd => MeshLoads%MooringLoads
   call MeshCopy( SrcMesh=Ptfm, DestMesh=MoorLd, CtrlCode=MESH_COUSIN, IOS=COMPONENT_OUTPUT, ErrStat=ErrStat2, ErrMess=ErrMsg2, Force=.true., Moment=.true. )
   if (Failed()) return
   PtfmLd%RemapFlag = .false.


   !-------------------------------
   ! create 2 point tower mesh
   !-------------------------------
   Twr => MeshMotions%TowerMotion
   call MeshCreate ( BlankMesh = Twr, IOS=COMPONENT_INPUT, Nnodes=2, ErrStat=ErrStat2, ErrMess=ErrMsg2,  &
               Orientation = .true., TranslationDisp = .true., TranslationVel = .true., RotationVel = .true., TranslationAcc  = .TRUE., RotationAcc = .true.)
   if (Failed()) return

   ! Tower bottom
   TmpPos(1:2) = real(TrbCfg%TowerBsPt(1:2) + TrbCfg%PtfmRefPos(1:2), ReKi)      ! relative to PtfmRefPos in (x,y)
   TmpPos(3)   = real(TrbCfg%TowerBsPt(3), ReKi)                                 ! relative to MSL in (z)
   call MeshPositionNode(Twr, 1, TmpPos, errStat2, errMsg2)  ! orientation is identity by default
   if (Failed()) return

   ! Tower top -- assumes vertical tower
   TmpPos(3) = real(TrbCfg%TowerHt,ReKi)   ! c_float to ReKi
   call MeshPositionNode(Twr, 2, TmpPos, errStat2, errMsg2)  ! orientation is identity by default
   if (Failed()) return

   ! create line element
   call MeshConstructElement( Twr, ELEMENT_LINE2, errStat2, errMsg2, p1=1, p2=2 )
   if (Failed()) return
  
   ! commit mesh          
   call MeshCommit(Twr, errStat2, errMsg2 )

   ! initialize location
   Twr%Orientation     = Twr%RefOrientation
   Twr%TranslationDisp = 0.0_R8Ki
   Twr%TranslationVel  = 0.0_ReKi
   Twr%RemapFlag = .false.

   ! create tower load mesh
   TwrLd => MeshLoads%TowerLoads
   call MeshCopy( SrcMesh=Twr, DestMesh=TwrLd, CtrlCode=MESH_SIBLING, IOS=COMPONENT_OUTPUT, ErrStat=ErrStat2, ErrMess=ErrMsg2, Force=.true., Moment=.true. )
   if (Failed()) return
   TwrLd%RemapFlag = .false.
   TwrLd%Force     = 0.0_ReKi
   TwrLd%Moment    = 0.0_ReKi


   !-------------------------------
   ! create hub mesh
   !-------------------------------
   !  NOTE: for a reference mesh position, nacelle yaw should be zero.  since NacYaw is static in this
   !        we are setting it once here.  If it needs to be dynamic, zero it here and update the yaw
   !        in the StructMotionUpdate routine below
   Hub => MeshMotions%HubMotion
   TmpPos(1:3) = Twr%Position(1:3,2)                                                ! Tower top
   TmpPos(1)   = TmpPos(1) + cos(TrbInit%NacYaw) * TrbCfg%OverHang                              ! X, nacelle yaw, and overhang
   TmpPos(2)   = TmpPos(2) + sin(TrbInit%NacYaw) * TrbCfg%OverHang                              ! Y, nacelle yaw, and overhang
   TmpPos(3)   = TmpPos(3) + TrbCfg%Twr2Shft     - abs(TrbCfg%OverHang) * tan(TrbCfg%ShftTilt)  ! Z, shaft height above tower top, and shaft tilt

   TmpAng = (/ 0.0_DbKi, -real(TrbCfg%ShftTilt,DbKi), real(TrbInit%NacYaw,DbKi)  /)             ! Hub/rotor azimuth is zero for reference. Hub axis on upwind points towards nacelle.
   Orient = EulerConstruct(TmpAng)
   call CreateInputPointMesh(Hub, TmpPos, Orient, ErrStat2, ErrMsg2, hasMotion=.true., hasLoads=.false.); if (Failed()) return
   Hub%RemapFlag = .false.

   ! create tower load mesh
   HubLd => MeshLoads%HubLoads
   call MeshCopy( SrcMesh=Hub, DestMesh=HubLd, CtrlCode=MESH_SIBLING, IOS=COMPONENT_OUTPUT, ErrStat=ErrStat2, ErrMess=ErrMsg2, Force=.true., Moment=.true. )
   if (Failed()) return
   HubLd%RemapFlag = .false.
   HubLd%Force     = 0.0_ReKi
   HubLd%Moment    = 0.0_ReKi


   !-------------------------------
   ! create blade root mesh
   !-------------------------------
   !  NOTE: for a reference mesh position, blade pitch should be zero.  since BldPitch is static in this
   !        we are setting it once here.  If it needs to be dynamic, zero it here and update the yaw
   !        in the StructMotionUpdate routine below
   allocate(MeshMotions%BladeRootMotion(TrbCfg%NumBl),STAT=ErrStat2)
   if (ErrStat2 /= 0) then
      ErrStat2 = ErrID_Fatal
      ErrMsg2  = "Could not allocate BladeRootMotion mesh"
      if (Failed()) return
   endif
   do k=1,TrbCfg%NumBl
      Root => MeshMotions%BladeRootMotion(k)
      AzBlade = TwoPi_D * real((k-1),DbKi)/real(TrbCfg%NumBl,DbKi)
      TmpAng = (/ AzBlade, real(TrbCfg%PreCone,DbKi), real(-TrbInit%BldPitch,DbKi) /)     ! Blade pitch does not follow RHR
      Orient = EulerConstruct(TmpAng)
      Orient = matmul(Orient,Hub%Orientation(1:3,1:3,1))
      TmpPos = Hub%Position(1:3,1) + TrbCfg%HubRad * real(Orient(3,1:3),ReKi)
      call CreateInputPointMesh(Root, TmpPos, Orient, ErrStat2, ErrMsg2, hasMotion=.true., hasLoads=.false.); if (Failed()) return
      Root%RemapFlag = .false.
   enddo

   ! create blade root load mesh
   allocate(MeshLoads%BladeRootLoads(TrbCfg%NumBl),STAT=ErrStat2)
   if (ErrStat2 /= 0) then
      ErrStat2 = ErrID_Fatal
      ErrMsg2  = "Could not allocate BladeRootLoads mesh"
      if (Failed()) return
   endif
   do k=1,TrbCfg%NumBl
      Root   => MeshMotions%BladeRootMotion(k)
      RootLd => MeshLoads%BladeRootLoads(k)
      call MeshCopy( SrcMesh=Root, DestMesh=RootLd, CtrlCode=MESH_SIBLING, IOS=COMPONENT_OUTPUT, ErrStat=ErrStat2, ErrMess=ErrMsg2, Force=.true., Moment=.true. )
      if (Failed()) return
      RootLd%RemapFlag = .false.
      RootLd%Force     = 0.0_ReKi
      RootLd%Moment    = 0.0_ReKi
   enddo

contains
   logical function Failed()
      call SetErrStat(errStat2, errMsg2, errStat, errMsg, RoutineName)
      Failed = errStat >= AbortErrLev
   end function Failed
end subroutine



!> create mesh mappings
subroutine StructCreateMeshMaps(SimSettings, MeshMotions, MeshLoads, MeshMaps, ErrStat, ErrMsg)
   type(SimSettingsType),           intent(in   )  :: SimSettings
   type(MeshesMotionType),          intent(inout)  :: MeshMotions
   type(MeshesLoadsType ),          intent(inout)  :: MeshLoads
   type(MeshesMapsType  ),          intent(inout)  :: MeshMaps
   integer(IntKi),                  intent(  out)  :: ErrStat
   character(ErrMsgLen),            intent(  out)  :: ErrMsg
   integer(IntKi)                                  :: ErrStat2
   character(ErrMsgLen)                            :: ErrMsg2
   character(*),           parameter               :: RoutineName = 'WaveTank::StructCreateMeshMaps'
   integer(IntKi)                                  :: k

   ErrStat = ErrID_None
   ErrMsg  = ''

   !-------------------------------
   ! Mapping arrays
   allocate(MeshMaps%Motion_Hub_2_BldRoot(SimSettings%TrbCfg%NumBl),STAT=ErrStat2)
   if (ErrStat2 /= 0) then
      ErrStat = ErrID_Fatal
      ErrMsg  = "Could not allocate Motion_Hub_2_BldRoot mesh mapping"
      return
   endif
   allocate(MeshMaps%Load_BldRoot_2_Hub(SimSettings%TrbCfg%NumBl),STAT=ErrStat2)
   if (ErrStat2 /= 0) then
      ErrStat = ErrID_Fatal
      ErrMsg  = "Could not allocate Load_BldRoot_2_Hub mesh mapping"
      return
   endif

   !-------------------------------
   ! Mesh motion mappings
   call MeshMapCreate(MeshMotions%PtfmPtMotion, MeshMotions%TowerMotion, MeshMaps%Motion_PRP_2_Twr, errStat2, errMsg2); if(Failed())return
   call MeshMapCreate(MeshMotions%PtfmPtMotion, MeshMotions%HubMotion,   MeshMaps%Motion_PRP_2_Hub, errStat2, errMsg2); if(Failed())return
   do k=1,SimSettings%TrbCfg%NumBl
      call MeshMapCreate(MeshMotions%HubMotion, MeshMotions%BladeRootMotion(k), MeshMaps%Motion_Hub_2_BldRoot(k), errStat2, errMsg2); if(Failed())return
   enddo
 
   !-------------------------------
   ! Mesh load mappings
   call MeshMapCreate(MeshLoads%TowerLoads, MeshLoads%PtfmPtLoads, MeshMaps%Load_Twr_2_PRP, errStat2, ErrMsg2); if(Failed()) return
   call MeshMapCreate(MeshLoads%HubLoads,   MeshLoads%PtfmPtLoads, MeshMaps%Load_Hub_2_PRP, errStat2, ErrMsg2); if(Failed()) return
   do k=1,SimSettings%TrbCfg%NumBl
      call MeshMapCreate(MeshLoads%BladeRootLoads(k), MeshLoads%HubLoads, MeshMaps%Load_BldRoot_2_Hub(k), errStat2, errMsg2); if(Failed())return
   enddo
   call MeshMapCreate(MeshLoads%MooringLoads, MeshLoads%PtfmPtLoads, MeshMaps%Load_Moor_2_PRP, errStat2, ErrMsg2); if(Failed()) return
  
contains
   logical function Failed()
      call SetErrStat(errStat2, errMsg2, errStat, errMsg, RoutineName)
      Failed = errStat >= AbortErrLev
   end function Failed
end subroutine



!> updates the structural meshes
subroutine StructMotionUpdate(SimSettings, CalcStepIO, MeshMotions, MeshMaps, StructTmp, ErrStat, ErrMsg)
   type(SimSettingsType),  target,  intent(in   )  :: SimSettings
   type(CalcStepIOdataType),        intent(in   )  :: CalcStepIO
   type(MeshesMotionType), target,  intent(inout)  :: MeshMotions
   type(MeshesMapsType  ),          intent(inout)  :: MeshMaps
   type(StructTmpType   ),          intent(inout)  :: StructTmp
   integer(IntKi),                  intent(  out)  :: ErrStat
   character(ErrMsgLen),            intent(  out)  :: ErrMsg
   integer(IntKi)                                  :: ErrStat2
   character(ErrMsgLen)                            :: ErrMsg2
   character(*),           parameter               :: RoutineName = 'WaveTank::StructMotionUpdate'
   real(R8Ki)                                      :: TmpTransDisp(3)
   real(DbKi)                                      :: TmpAng(3)   ! temporary euler angle
   real(R8Ki)                                      :: Orient(3,3) ! temporary orientation
   type(TurbConfigType),   pointer                 :: TrbCfg      ! to shorten notation
   type(TurbInitCondType), pointer                 :: TrbInit     ! to shorten notation
   type(MeshType),         pointer                 :: Ptfm        ! to shorten notation
   type(MeshType),         pointer                 :: Twr         ! to shorten notation
   type(MeshType),         pointer                 :: Hub         ! to shorten notation
   type(MeshType),         pointer                 :: Root        ! to shorten notation
   real(c_float)                                   :: ScaleFact   ! to shorten notation
   integer(IntKi)                                  :: k

   ErrStat = ErrID_None
   ErrMsg  = ''

   TrbCfg    => SimSettings%TrbCfg
   TrbInit   => SimSettings%TrbInit

   ! scaling factor
   ScaleFact = SimSettings%Sim%ScaleFact

   ! update PtfmPtMotion
   Ptfm => MeshMotions%PtfmPtMotion
   Ptfm%TranslationDisp(1:3,1) = FroudeScaleM2F_Disp(ScaleFact, CalcStepIO%PosAng_c(1:3), Ptfm%Position(1:3,1))
   Ptfm%Orientation(1:3,1:3,1) = WT_EulerToDCM_fromInput(CalcStepIO%PosAng_c(4:6))        ! angles don't scale
   Ptfm%TranslationVel(1:3,1)  = FroudeScaleM2F_TVel(ScaleFact, CalcStepIO%Vel_c(1:3))
   Ptfm%RotationVel(1:3,1)     = FroudeScaleM2F_RVel(ScaleFact, CalcStepIO%Vel_c(4:6))
   Ptfm%TranslationAcc(1:3,1)  = FroudeScaleM2F_TAcc(ScaleFact, CalcStepIO%Acc_c(1:3))
   Ptfm%RotationAcc(1:3,1)     = FroudeScaleM2F_RAcc(ScaleFact, CalcStepIO%Acc_c(4:6))

   !--------------------------------------
   ! transfer Ptfm to Tower
   Twr => MeshMotions%TowerMotion
   call Transfer_Point_to_Line2( Ptfm, Twr, MeshMaps%Motion_PRP_2_Twr, ErrStat2, ErrMsg2 ); if (Failed()) return;

   !--------------------------------------
   ! transfer Ptfm to hub (tower is rigid)
   Hub => MeshMotions%HubMotion
   call Transfer_Point_to_Point( Ptfm, Hub, MeshMaps%Motion_PRP_2_Hub, ErrStat2, ErrMsg2 ); if (Failed()) return;

   ! rotor azimuth
   StructTmp%Azimuth = modulo(real(CalcStepIO%Time_c,ReKi)*StructTmp%RotSpeed + TrbInit%Azimuth, TwoPi )
  
   ! update hub azimuth -- include initial azimuth
   TmpAng = (/ real(StructTmp%Azimuth,DbKi), 0.0_DbKi, 0.0_DbKi /)
   Orient = EulerConstruct(TmpAng)
   Hub%Orientation(1:3,1:3,1) = matmul(Orient,Hub%Orientation(1:3,1:3,1))

   !--------------------------------------
   ! hub to blades
   do k=1,TrbCfg%NumBl
      Root => MeshMotions%BladeRootMotion(k)
      call Transfer_Point_to_Point( Hub, Root, MeshMaps%Motion_Hub_2_BldRoot(k), ErrStat2, ErrMsg2 ); if (Failed()) return;
   enddo

contains
   logical function Failed()
      call SetErrStat(errStat2, errMsg2, errStat, errMsg, RoutineName)
      Failed = errStat >= AbortErrLev
   end function Failed
end subroutine



!> updates the structural load meshes and populates aggregated loads
subroutine StructLoadsMeshTransfer(SimSettings, CalcStepIO, MeshMotions, MeshLoads, MeshMaps, StructTmp, ErrStat, ErrMsg)
   type(SimSettingsType),  target,  intent(in   )  :: SimSettings
   type(CalcStepIOdataType),        intent(in   )  :: CalcStepIO
   type(MeshesMotionType), target,  intent(inout)  :: MeshMotions
   type(MeshesLoadsType),  target,  intent(inout)  :: MeshLoads
   type(MeshesMapsType  ),          intent(inout)  :: MeshMaps
   type(StructTmpType   ),          intent(inout)  :: StructTmp
   integer(IntKi),                  intent(  out)  :: ErrStat
   character(ErrMsgLen),            intent(  out)  :: ErrMsg
   integer(IntKi)                                  :: ErrStat2
   character(ErrMsgLen)                            :: ErrMsg2
   character(*),           parameter               :: RoutineName = 'WaveTank::StructMotionUpdate'
   real(R8Ki)                                      :: TmpTransDisp(3)
   real(DbKi)                                      :: TmpAng(3)               ! temporary euler angle
   real(R8Ki)                                      :: Orient(3,3)             ! temporary orientation
   type(TurbConfigType),   pointer                 :: TrbCfg                  ! to shorten notation
   type(TurbInitCondType), pointer                 :: TrbInit                 ! to shorten notation
   type(MeshType),         pointer                 :: Ptfm, PtfmLd, PtfmLdTmp ! to shorten notation
   type(MeshType),         pointer                 :: Twr,  TwrLd             ! to shorten notation
   type(MeshType),         pointer                 :: Hub,  HubLd             ! to shorten notation
   type(MeshType),         pointer                 :: Root, RootLd            ! to shorten notation
   type(MeshType),         pointer                 :: MoorLd                  ! to shorten notation
   real(c_float)                                   :: ScaleFact               ! to shorten notation
   integer(IntKi)                                  :: k

   ErrStat = ErrID_None
   ErrMsg  = ''

   ! shorthand pointers
   TrbCfg   => SimSettings%TrbCfg
   TrbInit  => SimSettings%TrbInit
   Hub      => MeshMotions%HubMotion
   HubLd    => MeshLoads%HubLoads
   Twr      => MeshMotions%TowerMotion
   TwrLd    => MeshLoads%TowerLoads
   Ptfm     => MeshMotions%PtfmPtMotion
   PtfmLd   => MeshLoads%PtfmPtLoads
   PtfmLdTmp=> MeshLoads%PtfmPtLoadsTmp
   MoorLd   => MeshLoads%MooringLoads

   !-----------------------------------------
   ! Aero loading
   !-----------------------------------------
   ! Transfer blade root loads to hub
   do k=1,TrbCfg%NumBl
      Root   => MeshMotions%BladeRootMotion(k)
      RootLd => MeshLoads%BladeRootLoads(k)
      call Transfer_Point_To_Point( RootLd, HubLd, MeshMaps%Load_BldRoot_2_Hub(k), ErrStat2, ErrMsg2, Root, Hub )
      if (Failed()) return
   enddo

   ! Transfer hub to platform
   call Transfer_Point_To_Point( HubLd, PtfmLd, MeshMaps%Load_Hub_2_PRP, ErrStat2, ErrMsg2, Hub, Ptfm )
   if (Failed()) return


   !-----------------------------------------
   ! Transfer tower to platform
   !  NOTE: no tower loads at present from ADI
   !FIXME: add tower loads output transfer to mesh here
   !call Transfer_Line2_To_Point( TwrLd, PtfmLdTmp, MeshMaps%Load_Twr_2_PRP, ErrStat2, ErrMsg2, Twr, Ptfm )
   !PtfmLd%Force(1:3,1)  = PtfmLd%Force(1:3,1)  + PtfmLdTmp%Force(1:3,1)
   !PtfmLd%Moment(1:3,1) = PtfmLd%Moment(1:3,1) + PtfmLdTmp%Moment(1:3,1)

   ! Store the ADI summed foreces and moments for output
   StructTmp%FrcMom_ADI_at_Ptfm(1:3) = PtfmLd%Force(1:3,1)
   StructTmp%FrcMom_ADI_at_Ptfm(4:6) = PtfmLd%Moment(1:3,1)


   !-----------------------------------------
   ! Mooring loading
   !-----------------------------------------
   ! Transfer mooring load
   call Transfer_Point_To_Point( MoorLd, PtfmLdTmp, MeshMaps%Load_Moor_2_PRP, ErrStat2, ErrMsg2, Ptfm, Ptfm )
   if (Failed()) return
   PtfmLd%Force(1:3,1)  = PtfmLd%Force(1:3,1)  + PtfmLdTmp%Force(1:3,1)
   PtfmLd%Moment(1:3,1) = PtfmLd%Moment(1:3,1) + PtfmLdTmp%Moment(1:3,1)

   ! Store the MD summed foreces and moments for output
   StructTmp%FrcMom_MD_at_Ptfm(1:3) = PtfmLdTmp%Force(1:3,1)
   StructTmp%FrcMom_MD_at_Ptfm(4:6) = PtfmLdTmp%Moment(1:3,1)

contains
   logical function Failed()
      call SetErrStat(errStat2, errMsg2, errStat, errMsg, RoutineName)
      Failed = errStat >= AbortErrLev
   end function Failed
end subroutine



!> destroy all structural model related info
subroutine StructDestroy(MeshMotions, MeshLoads, MeshMaps, StructTmp, ErrStat,ErrMsg)  ! We are actually ignoring all errors from here
   type(MeshesMotionType), intent(inout)  :: MeshMotions
   type(MeshesLoadsType ), intent(inout)  :: MeshLoads
   type(MeshesMapsType  ), intent(inout)  :: MeshMaps
   type(StructTmpType   ), intent(inout)  :: StructTmp
   integer(IntKi),         intent(  out)  :: ErrStat
   character(ErrMsgLen),   intent(  out)  :: ErrMsg
   integer(IntKi)                         :: ErrStat2
   character(ErrMsgLen)                   :: ErrMsg2
   character(*),           parameter      :: RoutineName = 'WaveTank::StructDestroy'
   ErrStat = ErrID_None
   ErrMsg  = ''
   call WT_DestroyMeshesMotionType(MeshMotions, ErrStat2, ErrMsg2)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   call WT_DestroyMeshesLoadsType(MeshLoads, ErrStat2, ErrMsg2)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   call WT_DestroyMeshesMapsType(MeshMaps, ErrStat2, ErrMsg2)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   call WT_DestroyStructTmpType(StructTmp, ErrStat2, ErrMsg2)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
end subroutine



!> Convert an Euler angle set of Roll, Pitch, Yaw ordering to a DCM.
!! this routine exists for two reasons
!!    1. ordering may be different
!!    2. incoming Euler angle is c_float instead of R8Ki
!! NOTE: no Euler angles are exported, so we stick with the OF convention
!!    for all internal conversions
function WT_EulerToDCM_fromInput(Ang) result(DCM)
   real(c_float), intent(in   )  :: Ang(3)
   real(R8Ki)                    :: DCM(3,3)
   !>>> Select one of the two following orders
   ! 3-2-1 intrinsic rotation sequence of the 3 Tait-Bryan angles (1-2-3 extrinsic rotation)
   DCM = EulerConstruct(real(Ang, DbKi))
   !! 1-2-3 intrinsic rotation sequence of the 3 Tait-Bryan angles (3-2-1 extrinsic rotation)
   !DCM = EulerConstructZYX(real(Ang, DbKi))
end function



subroutine WrVTK_Struct_Ref(SimSettings, MeshMotions, MeshLoads, ErrStat, ErrMsg)
   type(SimSettingsType),  target,  intent(in   )  :: SimSettings
   type(MeshesMotionType),          intent(in   )  :: MeshMotions
   type(MeshesLoadsType ),          intent(in   )  :: MeshLoads
   integer(IntKi),                  intent(  out)  :: ErrStat
   character(ErrMsgLen),            intent(  out)  :: ErrMsg
   integer(IntKi)                                  :: ErrStat2
   character(ErrMsgLen)                            :: ErrMsg2
   character(*),           parameter               :: RoutineName = 'WaveTank::WrVTK_Struct_Ref'
   character(1024)                                 :: DirRootName
   real(SiKi)                                      :: RefPt(3)
   integer(IntKi)                                  :: k
   ErrStat = ErrID_None
   ErrMsg  = ''
   RefPt   = (/ 0.0_SiKi, 0.0_SiKi, 0.0_SiKi /)
   DirRootName = trim(SimSettings%Viz%WrVTK_dir)//PathSep//trim(SimSettings%Sim%OutRootName)
   ! Wave elevation measurement buoy
   call MeshWrVTKreference(RefPt, MeshMotions%WaveBuoyMotion, trim(DirRootName)//'.WaveBuoyMotion', ErrStat2, ErrMsg2); if (Failed()) return
   ! Platform point
   call MeshWrVTKreference(RefPt, MeshMotions%PtfmPtMotion, trim(DirRootName)//'.Struct'//'.PtfmPtMotion', ErrStat2, ErrMsg2); if (Failed()) return
   ! Tower
   call MeshWrVTKreference(RefPt, MeshMotions%TowerMotion,  trim(DirRootName)//'.Struct'//'.TowerMotion',  ErrStat2, ErrMsg2); if (Failed()) return
   ! hub point
   call MeshWrVTKreference(RefPt, MeshMotions%HubMotion,    trim(DirRootName)//'.Struct'//'.HubMotion',    ErrStat2, ErrMsg2); if (Failed()) return
   ! RootMotion points
   do k=1,SimSettings%TrbCfg%NumBl
      call MeshWrVTKreference(RefPt, MeshMotions%BladeRootMotion(k), trim(DirRootName)//'.Struct'//'.RootMotion'//trim(Num2LStr(k)), ErrStat2, ErrMsg2); if (Failed()) return
   enddo
contains
   logical function Failed()
      call SetErrStat(errStat2, errMsg2, errStat, errMsg, RoutineName)
      Failed = errStat >= AbortErrLev
   end function Failed
end subroutine

subroutine WrVTK_Struct(n_Global, SimSettings, MeshMotions, MeshLoads, ErrStat, ErrMsg)
   integer(IntKi),                  intent(in   )  :: n_Global
   type(SimSettingsType),  target,  intent(in   )  :: SimSettings
   type(MeshesMotionType),          intent(in   )  :: MeshMotions
   type(MeshesLoadsType ),          intent(in   )  :: MeshLoads
   integer(IntKi),                  intent(  out)  :: ErrStat
   character(ErrMsgLen),            intent(  out)  :: ErrMsg
   integer(IntKi)                                  :: ErrStat2
   character(ErrMsgLen)                            :: ErrMsg2
   character(*),           parameter               :: RoutineName = 'WaveTank::WrVTK_Struct'
   character(1024)                                 :: DirRootName
   real(SiKi)                                      :: RefPt(3)
   integer(IntKi)                                  :: k
   ErrStat = ErrID_None
   ErrMsg  = ''
   RefPt   = (/ 0.0_SiKi, 0.0_SiKi, 0.0_SiKi /)
   DirRootName = trim(SimSettings%Viz%WrVTK_dir)//PathSep//trim(SimSettings%Sim%OutRootName)
   ! Wave elevation measurement buoy
   call MeshWrVTK(RefPt, MeshMotions%WaveBuoyMotion, trim(DirRootName)//'.WaveBuoyMotion', n_Global, .true., ErrStat2, ErrMsg2, Twidth=SimSettings%Viz%Twidth); if (Failed()) return
   ! Platform point
   call MeshWrVTK(RefPt, MeshMotions%PtfmPtMotion, trim(DirRootName)//'.Struct'//'.PtfmPtMotion', n_Global, .true., ErrStat2, ErrMsg2, Twidth=SimSettings%Viz%Twidth); if (Failed()) return
   ! Tower
   call MeshWrVTK(RefPt, MeshMotions%TowerMotion,  trim(DirRootName)//'.Struct'//'.TowerMotion',  n_Global, .true., ErrStat2, ErrMsg2, Twidth=SimSettings%Viz%Twidth); if (Failed()) return
   ! Hub point
   call MeshWrVTK(RefPt, MeshMotions%HubMotion,    trim(DirRootName)//'.Struct'//'.HubMotion',    n_Global, .true., ErrStat2, ErrMsg2, Twidth=SimSettings%Viz%Twidth); if (Failed()) return
   ! RootMotion points
   do k=1,SimSettings%TrbCfg%NumBl
      call MeshWrVTK(RefPt, MeshMotions%BladeRootMotion(k), trim(DirRootName)//'.Struct'//'.RootMotion'//trim(Num2LStr(k)), n_Global, .true., ErrStat2, ErrMsg2, Twidth=SimSettings%Viz%Twidth); if (Failed()) return
   enddo
contains
   logical function Failed()
      call SetErrStat(errStat2, errMsg2, errStat, errMsg, RoutineName)
      Failed = errStat >= AbortErrLev
   end function Failed
end subroutine



!-----------------------------------------------
! Froude scaling from here: https://home.hvl.no/ansatte/gste/ftp/MarinLab_files/Litteratur/NTNU_Scaling_Laws.pdf, page 21
! notation below:
!     model scale: _m
!     full scale:  _f
!     ScaleFact:   length_f/length_m = lambda
!     DensFact:    rho_f/rho_m

!> scale model displacements to full scale 
!! length_full = length_model * lambda
function FroudeScaleM2F_Disp(ScaleFact, Pos_m, refPos_f) result(transDisp_f)
   real(c_float), intent(in   ) :: ScaleFact
   real(c_float), intent(in   ) :: Pos_m(3)
   real(ReKi),    intent(in   ) :: refPos_f(3)
   real(R8Ki)                   :: transdisp_f(3)
   transDisp_f = real(ScaleFact*Pos_m,R8Ki) - real(refPos_f,R8Ki)
end function

!> scale model translational velocity to full scale 
!! TVel_full = TVel_model * sqrt( lambda )         TODO: check this!!!! 
function FroudeScaleM2F_TVel(ScaleFact, TVel_m) result(TVel_f)
   real(c_float), intent(in   ) :: ScaleFact
   real(c_float), intent(in   ) :: TVel_m(3)
   real(ReKi)                   :: TVel_f(3)
   TVel_f = sqrt(real(ScaleFact,ReKi)) * real(TVel_m,ReKi)
end function

!> scale model rotational velocity to full scale 
!! RVel_full = RVel_model * sqrt(lambda)           TODO: check this!!!! 
function FroudeScaleM2F_RVel(ScaleFact, RVel_m) result(RVel_f)
   real(c_float), intent(in   ) :: ScaleFact
   real(c_float), intent(in   ) :: RVel_m(3)
   real(ReKi)                   :: RVel_f(3)
   RVel_f =  real(RVel_m,ReKi) / sqrt(real(ScaleFact,ReKi))
end function

!> scale model translational acceleration to full scale 
!! TAcc_full = TAcc_model ---> no scaling applied
function FroudeScaleM2F_TAcc(ScaleFact, TAcc_m) result(TAcc_f)
   real(c_float), intent(in   ) :: ScaleFact
   real(c_float), intent(in   ) :: TAcc_m(3)
   real(ReKi)                   :: TAcc_f(3)
   TAcc_f = real(TAcc_m,ReKi)
end function

!> scale model rotational acceleration to full scale 
!! RAcc_full = RAcc_model / lambda                 TODO: check this!!!! 
function FroudeScaleM2F_RAcc(ScaleFact, RAcc_m) result(RAcc_f)
   real(c_float), intent(in   ) :: ScaleFact
   real(c_float), intent(in   ) :: RAcc_m(3)
   real(ReKi)                   :: RAcc_f(3)
   RAcc_f = real(RAcc_m,ReKi) / real(ScaleFact,ReKi)
end function

!> scale model time to full scale
!! sqrt(lambda) = sqrt(Length_full/ Length_model)
function FroudeScaleM2F_Time(ScaleFact, Time_m) result(Time_f)
   real(c_float), intent(in   ) :: ScaleFact
   real(c_double),intent(in   ) :: Time_m
   real(R8Ki)                   :: Time_f
   Time_f = sqrt(real(ScaleFact,R8Ki)) * real(Time_m,R8Ki)
end function

!> scale full scale force to model
!! lambda^3 * DensFact * Frc_model = Frc_full
function FroudeScaleF2M_Frc(ScaleFact, DensFact, Frc_f) result(Frc_m)
   real(c_float), intent(in   ) :: ScaleFact
   real(c_float), intent(in   ) :: DensFact
   real(ReKi),    intent(in   ) :: Frc_f(3)
   real(c_float)                :: Frc_m(3)
   Frc_m = real(Frc_f, c_float) / (ScaleFact**3 * DensFact)
end function

!> scale full scale moment to model
!! lambda^4 * DensFact * Mom_model = Mom_full
function FroudeScaleF2M_Mom(ScaleFact, DensFact, Mom_f) result(Mom_m)
   real(c_float), intent(in   ) :: ScaleFact
   real(c_float), intent(in   ) :: DensFact
   real(ReKi),    intent(in   ) :: Mom_f(3)
   real(c_float)                :: Mom_m(3)
   Mom_m = real(Mom_f, c_float) / (ScaleFact**4 * DensFact)
end function


end module
