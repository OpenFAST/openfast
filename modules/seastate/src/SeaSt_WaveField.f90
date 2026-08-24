MODULE SeaSt_WaveField

USE GridInterp
USE SeaSt_WaveField_Types
USE IfW_FlowField, only: IfW_FlowField_GetVelAcc
USE GridInterp_Types
USE Waves,  ONLY: WaveKinKernel_ComputeColumns          ! shared first-order per-column generation kernel (on-demand block population)
USE Waves2, ONLY: WaveKinKernel_AddSecondOrderColumns   ! shared second-order per-column kernel (on-demand block population)
USE VTK, ONLY: VTK_Misc, vtk_misc_init, vtk_new_ascii_file, vtk_dataset_rectilinear, &
               vtk_cell_data_init, vtk_cell_data_scalar, vtk_close_file   ! block VTK export

IMPLICIT NONE

PRIVATE

! Public functions and subroutines
PUBLIC WaveField_GetNodeTotalWaveElev
PUBLIC WaveField_GetMinMaxWaveElevEstimate
PUBLIC WaveField_GetNodeWaveNormal
PUBLIC WaveField_GetNodeWaveKin
PUBLIC WaveField_GetNodeWaveVel
PUBLIC WaveField_GetNodeWaveVelAcc
PUBLIC WaveField_GetWaveKin
PUBLIC WaveField_GetWaveVelAcc_AD
PUBLIC WaveField_GetMeanDynSurfCurr
PUBLIC WaveField_GetDynP
PUBLIC WaveField_BlockStore_Init
PUBLIC WaveField_WriteBlockVTK

CONTAINS

!-------------------- Subroutine for wave elevation ------------------!

FUNCTION WaveField_GetNodeTotalWaveElev( WaveField, WaveField_m, Time, pos, ErrStat, ErrMsg, Elev1, Elev2 )
   type(SeaSt_WaveFieldType),          intent(in   ) :: WaveField
   type(GridInterp_MiscVarType),       intent(inout) :: WaveField_m
   real(DbKi),                         intent(in   ) :: Time
   real(ReKi),                         intent(in   ) :: pos(*)  ! Position at which free-surface elevation is to be calculated. Third entry ignored if present.
   integer(IntKi),                     intent(  out) :: ErrStat ! Error status of the operation
   character(*),                       intent(  out) :: ErrMsg  ! Error message if errStat /= ErrID_None
   real(SiKi), optional,               intent(  out) :: Elev1, Elev2 ! Elev1 and Elev2 components

   real(SiKi)                                        :: WaveField_GetNodeTotalWaveElev
   real(SiKi)                                        :: Zeta1, Zeta2
   character(*),                       parameter     :: RoutineName = 'WaveField_GetNodeTotalWaveElev'
   integer(IntKi)                                    :: errStat2
   character(ErrMsgLen)                              :: errMsg2

   ErrStat   = ErrID_None

   IF (ALLOCATED(WaveField%WaveElev1) .or. ALLOCATED(WaveField%WaveElev2)) then
      CALL WaveField_Interp_Setup3D(Time, pos, WaveField%SrfGridParams, WaveField_m, ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) return
   end if

   IF (ALLOCATED(WaveField%WaveElev1)) THEN
      Zeta1 = GridInterp3D(WaveField%WaveElev1, WaveField_m)
   ELSE
      Zeta1 = 0.0_SiKi
   END IF

   IF (ALLOCATED(WaveField%WaveElev2)) THEN
      Zeta2 = GridInterp3D(WaveField%WaveElev2, WaveField_m)
   ELSE
      Zeta2 = 0.0_SiKi
   END IF

   if (present(Elev1)) Elev1 = Zeta1
   if (present(Elev2)) Elev2 = Zeta2

   WaveField_GetNodeTotalWaveElev = Zeta1 + Zeta2

END FUNCTION WaveField_GetNodeTotalWaveElev


!> Gives an estimate of the min and max wave elevation.  It will overshoot for second order
subroutine WaveField_GetMinMaxWaveElevEstimate( WaveField, MinElev, MaxElev, ErrStat, ErrMsg )
   type(SeaSt_WaveFieldType), pointer, intent(in   ) :: WaveField
   real(SiKi),                         intent(  out) :: MinElev
   real(SiKi),                         intent(  out) :: MaxElev
   integer(IntKi),                     intent(  out) :: ErrStat ! Error status of the operation
   character(*),                       intent(  out) :: ErrMsg  ! Error message if errStat /= ErrID_None
   character(*),                       parameter     :: RoutineName = 'WaveField_GetMinMaxWaveElevEstimate'

   ErrStat   = ErrID_None
   ErrMsg    = ""
   MinElev = 0.0_SiKi
   MaxElev = 0.0_SiKi

   ! Check that data exists
   if (.not. associated(WaveField)) then
      ErrStat = ErrID_Fatal
      ErrMsg  = trim(RoutineName)//": WaveField data does not exist."
      return
   endif

   if (allocated(WaveField%WaveElev1)) then
      MinElev = minval(WaveField%WaveElev1)
      MaxElev = maxval(WaveField%WaveElev1)
   endif
   if (allocated(WaveField%WaveElev2)) then
      MinElev = MinElev + minval(WaveField%WaveElev2)
      MaxElev = MaxElev + maxval(WaveField%WaveElev2)
   endif
end subroutine WaveField_GetMinMaxWaveElevEstimate

SUBROUTINE WaveField_GetNodeWaveNormal( WaveField, WaveField_m, Time, pos, n, ErrStat, ErrMsg )
   type(SeaSt_WaveFieldType),          intent(in   ) :: WaveField
   type(GridInterp_MiscVarType),       intent(inout) :: WaveField_m
   real(DbKi),                         intent(in   ) :: Time
   real(ReKi),                         intent(in   ) :: pos(:)  ! Position at which free-surface normal is to be calculated. Third entry ignored if present.
   real(ReKi),                         intent(  out) :: n(3)    ! Free-surface normal vector
   integer(IntKi),                     intent(  out) :: ErrStat ! Error status of the operation
   character(*),                       intent(  out) :: ErrMsg  ! Error message if errStat /= ErrID_None

   real(SiKi)                                        :: slope(2)
   character(*),                       parameter     :: RoutineName = 'WaveField_GetNodeWaveNormal'
   integer(IntKi)                                    :: errStat2
   character(ErrMsgLen)                              :: errMsg2

   ErrStat   = ErrID_None

   call GridInterpSetupN( (/Real(Time+WaveField%WaveTimeShift,ReKi),pos(1),pos(2)/), WaveField%SrfGridParams, WaveField_m, ErrStat2, ErrMsg2 )
   slope = GridInterpS( WaveField%WaveElev1, WaveField%SrfGridParams, WaveField_m )
   if (ALLOCATED(WaveField%WaveElev2)) then
      slope = slope + GridInterpS( WaveField%WaveElev2, WaveField%SrfGridParams, WaveField_m )
   end if

   n = Real( (/-slope(1),-slope(2),1.0_SiKi/), ReKi)
   n = n / TwoNorm(n)

contains
   logical function Failed()
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
   end function
END SUBROUTINE WaveField_GetNodeWaveNormal


!-------------------- Subroutine for full wave field kinematics --------------------!
SUBROUTINE WaveField_GetNodeWaveKin( WaveField, WaveField_m, Time, pos, forceNodeInWater, fetchDynCurrent, nodeInWater, WaveElev1, WaveElev2, WaveElev, FDynP, FV, FA, FAMCF, ErrStat, ErrMsg )
   type(SeaSt_WaveFieldType),          intent(in   ) :: WaveField
   type(GridInterp_MiscVarType),       intent(inout) :: WaveField_m
   real(DbKi),                         intent(in   ) :: Time
   real(ReKi),                         intent(in   ) :: pos(3)
   logical,                            intent(in   ) :: forceNodeInWater
   logical,                            intent(in   ) :: fetchDynCurrent
   real(SiKi),                         intent(  out) :: WaveElev1
   real(SiKi),                         intent(  out) :: WaveElev2
   real(SiKi),                         intent(  out) :: WaveElev
   real(SiKi),                         intent(  out) :: FV(3)
   real(SiKi),                         intent(  out) :: FA(3)
   real(SiKi),                         intent(  out) :: FAMCF(3)
   real(SiKi),                         intent(  out) :: FDynP
   integer(IntKi),                     intent(  out) :: nodeInWater
   integer(IntKi),                     intent(  out) :: ErrStat ! Error status of the operation
   character(*),                       intent(  out) :: ErrMsg  ! Error message if errStat /= ErrID_None

   real(ReKi)                                        :: posXY(2), posPrime(3), posXY0(3), PosOffset(3), posDummy(3,1)
   integer(IntKi)                                    :: startNode
   real(ReKi), allocatable                           :: FV_DC(:,:), FA_DC(:,:)
   character(*),                       parameter     :: RoutineName = 'WaveField_GetNodeWaveKin'
   integer(IntKi)                                    :: errStat2
   character(ErrMsgLen)                              :: errMsg2

   ErrStat   = ErrID_None

   posXY    = pos(1:2)
   posXY0   = (/pos(1),pos(2),0.0_ReKi/)
   FAMCF(:) = 0.0

   ! Wave elevation (Calls WaveField_Interp_Setup3D internally so WaveField_Interp_3D can be used below)
   WaveElev = WaveField_GetNodeTotalWaveElev(WaveField, WaveField_m, Time, pos, ErrStat2, ErrMsg2, Elev1=WaveElev1, Elev2=WaveElev2)
   if (Failed()) return

   ! Check if point is below the seabed
   if (pos(3)<-WaveField%EffWtrDpth) then
      nodeInWater = 1_IntKi  ! Prevent problems with HydroDyn logic
      FV(:)       = 0.0_SiKi
      FA(:)       = 0.0_SiKi
      FDynP       = 0.0_SiKi
      FAMCF(:)    = 0.0_SiKi
      return
   end if

   IF (WaveField%WaveStMod == 0) THEN ! No wave stretching

      IF ( pos(3) <= 0.0_ReKi) THEN ! Node is at or below the SWL
         nodeInWater = 1_IntKi
         ! Use location to obtain interpolated values of kinematics
         CALL WaveField_Interp_Setup4D( Time+WaveField%WaveTimeShift, pos, WaveField%GridDepth, WaveField%VolGridParams, WaveField_m, ErrStat2, ErrMsg2 ); if (Failed()) return;
         CALL WaveField_InterpVol( WaveField, WaveField_m, Time, ErrStat2, ErrMsg2, FDynP=FDynP, FV=FV, FA=FA, FAMCF=FAMCF ); if (Failed()) return;
      ELSE ! Node is above the SWL
         nodeInWater = 0_IntKi
         FV(:)       = 0.0_SiKi
         FA(:)       = 0.0_SiKi
         FDynP       = 0.0_SiKi
         FAMCF(:)    = 0.0_SiKi
      END IF

   ELSE ! Wave stretching enabled

      IF ( (pos(3) <= WaveElev) .OR. forceNodeInWater ) THEN ! Node is submerged

         nodeInWater = 1_IntKi

         IF ( WaveField%WaveStMod < 3 ) THEN ! Vertical or extrapolated wave stretching

            IF ( pos(3) <= 0.0_SiKi) THEN ! Node is below the SWL - evaluate wave dynamics as usual

               ! Use location to obtain interpolated values of kinematics
               CALL WaveField_Interp_Setup4D( Time+WaveField%WaveTimeShift, pos, WaveField%GridDepth, WaveField%VolGridParams, WaveField_m, ErrStat2, ErrMsg2 ); if (Failed()) return;
               CALL WaveField_InterpVol( WaveField, WaveField_m, Time, ErrStat2, ErrMsg2, FDynP=FDynP, FV=FV, FA=FA, FAMCF=FAMCF ); if (Failed()) return;

            ELSE ! Node is above SWL - need wave stretching

               ! Vertical wave stretching
               CALL WaveField_Interp_Setup4D( Time+WaveField%WaveTimeShift, posXY0, WaveField%GridDepth, WaveField%VolGridParams, WaveField_m, ErrStat2, ErrMsg2 ); if (Failed()) return;
               CALL WaveField_InterpVol( WaveField, WaveField_m, Time, ErrStat2, ErrMsg2, FDynP=FDynP, FV=FV, FA=FA, FAMCF=FAMCF ); if (Failed()) return;

               ! Extrapolated wave stretching
               IF (WaveField%WaveStMod == 2) THEN
                  CALL WaveField_Interp_Setup3D( Time+WaveField%WaveTimeShift, posXY, WaveField%SrfGridParams, WaveField_m, ErrStat2, ErrMsg2 ); if (Failed()) return;
                  FV(:) = FV(:) + GridInterp3DVec( WaveField%PWaveVel0,  WaveField_m ) * pos(3)
                  FA(:) = FA(:) + GridInterp3DVec( WaveField%PWaveAcc0,  WaveField_m ) * pos(3)
                  FDynP = FDynP + GridInterp3D   ( WaveField%PWaveDynP0, WaveField_m ) * pos(3)
                  IF ( WaveField_HasMCF(WaveField) ) THEN
                     FAMCF(:) = FAMCF(:) + GridInterp3DVec( WaveField%PWaveAccMCF0, WaveField_m ) * pos(3)
                  END IF
               END IF

            END IF ! Node is submerged

         ELSE ! Wheeler stretching - no need to check whether the node is above or below SWL

            ! Map the node z-position linearly from [-EffWtrDpth,m%WaveElev(j)] to [-EffWtrDpth,0]
            posPrime    = pos
            posPrime(3) = WaveField%EffWtrDpth*(WaveField%EffWtrDpth+pos(3))/(WaveField%EffWtrDpth+WaveElev)-WaveField%EffWtrDpth
            posPrime(3) = MIN( posPrime(3), 0.0_ReKi) ! Clamp z-position to zero. Needed when forceNodeInWater=.TRUE.

            ! Obtain the wave-field variables by interpolation with the mapped position.
            CALL WaveField_Interp_Setup4D( Time+WaveField%WaveTimeShift, posPrime, WaveField%GridDepth, WaveField%VolGridParams, WaveField_m, ErrStat2, ErrMsg2 ); if (Failed()) return;
            CALL WaveField_InterpVol( WaveField, WaveField_m, Time, ErrStat2, ErrMsg2, FDynP=FDynP, FV=FV, FA=FA, FAMCF=FAMCF ); if (Failed()) return;
         END IF

      ELSE ! Node is out of water - zero-out all wave dynamics

         nodeInWater = 0_IntKi
         FV(:)       = 0.0_SiKi
         FA(:)       = 0.0_SiKi
         FDynP       = 0.0_SiKi
         FAMCF(:)    = 0.0_SiKi

      END IF ! If node is in or out of water

   END IF ! If wave stretching is on or off

   ! Get dynamic current velocity and acceleration
   IF (fetchDynCurrent .AND. WaveField%hasCurrField) THEN
      startNode = -1
      PosOffset = (/0.0_ReKi,0.0_ReKi,WaveField%EffWtrDpth/)
      posDummy(:,1) = pos
      ALLOCATE(FV_DC(3,1), STAT=ErrStat2); if (FailedMsg('Error allocating FV_DC')) return;
      ALLOCATE(FA_DC(3,1), STAT=ErrStat2); if (FailedMsg('Error allocating FA_DC')) return;    
      CALL IfW_FlowField_GetVelAcc(WaveField%CurrField, startNode, Time, posDummy, FV_DC, FA_DC, ErrStat2, ErrMsg2, PosOffset=PosOffset); if (Failed()) return;
      FV = FV + nodeInWater * FV_DC(:,1)
      FA = FA + nodeInWater * FA_DC(:,1)
   END IF

contains
   logical function Failed()
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
   end function
   logical function FailedMsg(ErrMsg2)
      character(*), intent(in   ) :: ErrMsg2
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      FailedMsg = ErrStat >= AbortErrLev
   end function
END SUBROUTINE WaveField_GetNodeWaveKin


!-------------------- Subroutine for dynamic pressure --------------------!
!NOTE: this is a subset of WaveField_GetNodeWaveKin
SUBROUTINE WaveField_GetDynP( WaveField, WaveField_m, Time, pos, forceNodeInWater, nodeInWater, FDynP, ErrStat, ErrMsg )
   type(SeaSt_WaveFieldType),          intent(in   ) :: WaveField
   type(GridInterp_MiscVarType),       intent(inout) :: WaveField_m
   real(DbKi),                         intent(in   ) :: Time
   real(ReKi),                         intent(in   ) :: pos(3)
   logical,                            intent(in   ) :: forceNodeInWater
   real(SiKi),                         intent(  out) :: FDynP
   integer(IntKi),                     intent(  out) :: nodeInWater
   integer(IntKi),                     intent(  out) :: ErrStat ! Error status of the operation
   character(*),                       intent(  out) :: ErrMsg  ! Error message if errStat /= ErrID_None

   real(ReKi)                                        :: posXY(2), posPrime(3), posXY0(3)
   character(*),                       parameter     :: RoutineName = 'WaveField_GetDynP'
   integer(IntKi)                                    :: errStat2
   character(ErrMsgLen)                              :: errMsg2

   ! Temporary vars not kept
   real(SiKi)                                        :: WaveElev

   ErrStat   = ErrID_None
   ErrMsg    = ""

   posXY    = pos(1:2)
   posXY0   = (/pos(1),pos(2),0.0_ReKi/)

   ! Wave elevation (Calls WaveField_Interp_Setup3D internally so WaveField_Interp_3D_vec can be used below)
   WaveElev  = WaveField_GetNodeTotalWaveElev( WaveField, WaveField_m, Time, pos, ErrStat2, ErrMsg2 ); if (Failed()) return;

   ! Check if point is below the seabed
   if (pos(3)<-WaveField%EffWtrDpth) then
      nodeInWater = 1_IntKi  ! Prevent problems with HydroDyn logic
      FDynP       = 0.0_SiKi
      return
   end if

   IF (WaveField%WaveStMod == 0) THEN ! No wave stretching

      IF ( pos(3) <= 0.0_ReKi) THEN ! Node is at or below the SWL
         nodeInWater = 1_IntKi
         ! Use location to obtain interpolated values of kinematics
         CALL WaveField_Interp_Setup4D( Time+WaveField%WaveTimeShift, pos, WaveField%GridDepth, WaveField%VolGridParams, WaveField_m, ErrStat2, ErrMsg2 ); if (Failed()) return;
         CALL WaveField_InterpVol( WaveField, WaveField_m, Time, ErrStat2, ErrMsg2, FDynP=FDynP ); if (Failed()) return;
      ELSE ! Node is above the SWL
         nodeInWater = 0_IntKi
         FDynP       = 0.0_SiKi
      END IF

   ELSE ! Wave stretching enabled
      IF ( (pos(3) <= WaveElev) .OR. forceNodeInWater ) THEN ! Node is submerged
         nodeInWater = 1_IntKi
         IF ( WaveField%WaveStMod < 3 ) THEN ! Vertical or extrapolated wave stretching

            IF ( pos(3) <= 0.0_SiKi) THEN ! Node is below the SWL - evaluate wave dynamics as usual
               ! Use location to obtain interpolated values of kinematics
               CALL WaveField_Interp_Setup4D( Time+WaveField%WaveTimeShift, pos, WaveField%GridDepth, WaveField%VolGridParams, WaveField_m, ErrStat2, ErrMsg2 ); if (Failed()) return;
               CALL WaveField_InterpVol( WaveField, WaveField_m, Time, ErrStat2, ErrMsg2, FDynP=FDynP ); if (Failed()) return;
            ELSE ! Node is above SWL - need wave stretching

               ! Vertical wave stretching
               CALL WaveField_Interp_Setup4D( Time+WaveField%WaveTimeShift, posXY0, WaveField%GridDepth, WaveField%VolGridParams, WaveField_m, ErrStat2, ErrMsg2 ); if (Failed()) return;
               CALL WaveField_InterpVol( WaveField, WaveField_m, Time, ErrStat2, ErrMsg2, FDynP=FDynP ); if (Failed()) return;

               ! Extrapoled wave stretching
               IF (WaveField%WaveStMod == 2) THEN
                  CALL WaveField_Interp_Setup3D( Time+WaveField%WaveTimeShift, posXY, WaveField%SrfGridParams, WaveField_m, ErrStat2, ErrMsg2 ); if (Failed()) return;
                  FDynP = FDynP + GridInterp3D   ( WaveField%PWaveDynP0, WaveField_m ) * pos(3)
               END IF

            END IF ! Node is submerged

         ELSE ! Wheeler stretching - no need to check whether the node is above or below SWL

            ! Map the node z-position linearly from [-EffWtrDpth,m%WaveElev(j)] to [-EffWtrDpth,0]
            posPrime    = pos
            posPrime(3) = WaveField%EffWtrDpth*(WaveField%EffWtrDpth+pos(3))/(WaveField%EffWtrDpth+WaveElev)-WaveField%EffWtrDpth
            posPrime(3) = MIN( posPrime(3), 0.0_ReKi) ! Clamp z-position to zero. Needed when forceNodeInWater=.TRUE.

            ! Obtain the wave-field variables by interpolation with the mapped position.
            CALL WaveField_Interp_Setup4D( Time+WaveField%WaveTimeShift, posPrime, WaveField%GridDepth, WaveField%VolGridParams, WaveField_m, ErrStat2, ErrMsg2 ); if (Failed()) return;
            CALL WaveField_InterpVol( WaveField, WaveField_m, Time, ErrStat2, ErrMsg2, FDynP=FDynP ); if (Failed()) return;
         END IF

      ELSE ! Node is out of water - zero-out all wave dynamics

         nodeInWater = 0_IntKi
         FDynP       = 0.0_SiKi

      END IF ! If node is in or out of water
   END IF ! If wave stretching is on or off

contains
   logical function Failed()
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
   end function
END SUBROUTINE WaveField_GetDynP


!-------------------- Subroutine for wave field velocity only --------------------!
SUBROUTINE WaveField_GetNodeWaveVel( WaveField, WaveField_m, Time, pos, forceNodeInWater, fetchDynCurrent, nodeInWater, FV, ErrStat, ErrMsg )
   type(SeaSt_WaveFieldType),          intent(in   ) :: WaveField
   type(GridInterp_MiscVarType),       intent(inout) :: WaveField_m
   real(DbKi),                         intent(in   ) :: Time
   real(ReKi),                         intent(in   ) :: pos(3)
   logical,                            intent(in   ) :: forceNodeInWater
   logical,                            intent(in   ) :: fetchDynCurrent
   integer(IntKi),                     intent(  out) :: nodeInWater
   real(SiKi),                         intent(  out) :: FV(3)
   integer(IntKi),                     intent(  out) :: ErrStat ! Error status of the operation
   character(*),                       intent(  out) :: ErrMsg  ! Error message if errStat /= ErrID_None

   real(SiKi)                                        :: WaveElev
   real(ReKi)                                        :: posXY(2), posPrime(3), posXY0(3), PosOffset(3), posDummy(3,1)
   real(ReKi), allocatable                           :: FV_DC(:,:), FA_DC(:,:)
   integer(IntKi)                                    :: startNode
   character(*),                       parameter     :: RoutineName = 'WaveField_GetNodeWaveVel'
   integer(IntKi)                                    :: errStat2
   character(ErrMsgLen)                              :: errMsg2

   ErrStat   = ErrID_None

   posXY    = pos(1:2)
   posXY0   = (/pos(1),pos(2),0.0_ReKi/)

   ! Wave elevation (Calls WaveField_Interp_Setup3D internally so WaveField_Interp_3D_vec can be used below)
   WaveElev  = WaveField_GetNodeTotalWaveElev( WaveField, WaveField_m, Time, pos, ErrStat2, ErrMsg2 ); if (Failed()) return;

   ! Check if point is below the seabed
   if (pos(3)<-WaveField%EffWtrDpth) then
      nodeInWater = 1_IntKi  ! Prevent problems with HydroDyn logic
      FV(:)       = 0.0_SiKi
      return
   end if

   IF (WaveField%WaveStMod == 0) THEN ! No wave stretching

      IF ( pos(3) <= 0.0_ReKi) THEN ! Node is at or below the SWL
         nodeInWater = 1_IntKi
         ! Use location to obtain interpolated values of kinematics
         CALL WaveField_Interp_Setup4D( Time+WaveField%WaveTimeShift, pos, WaveField%GridDepth, WaveField%VolGridParams, WaveField_m, ErrStat2, ErrMsg2 ); if (Failed()) return;
         CALL WaveField_InterpVol( WaveField, WaveField_m, Time, ErrStat2, ErrMsg2, FV=FV ); if (Failed()) return;
      ELSE ! Node is above the SWL
         nodeInWater = 0_IntKi
         FV(:)       = 0.0_SiKi
      END IF

   ELSE ! Wave stretching enabled

      IF ( (pos(3) <= WaveElev) .OR. forceNodeInWater ) THEN ! Node is submerged

         nodeInWater = 1_IntKi

         IF ( WaveField%WaveStMod < 3 ) THEN ! Vertical or extrapolated wave stretching

            IF ( pos(3) <= 0.0_SiKi) THEN ! Node is below the SWL - evaluate wave dynamics as usual

               ! Use location to obtain interpolated values of kinematics
               CALL WaveField_Interp_Setup4D( Time+WaveField%WaveTimeShift, pos, WaveField%GridDepth, WaveField%VolGridParams, WaveField_m, ErrStat2, ErrMsg2 ); if (Failed()) return;
               CALL WaveField_InterpVol( WaveField, WaveField_m, Time, ErrStat2, ErrMsg2, FV=FV ); if (Failed()) return;

            ELSE ! Node is above SWL - need wave stretching

               ! Vertical wave stretching
               CALL WaveField_Interp_Setup4D( Time, posXY0, WaveField%GridDepth, WaveField%VolGridParams, WaveField_m, ErrStat2, ErrMsg2 ); if (Failed()) return;
               CALL WaveField_InterpVol( WaveField, WaveField_m, Time, ErrStat2, ErrMsg2, FV=FV ); if (Failed()) return;

               ! Extrapolated wave stretching
               IF (WaveField%WaveStMod == 2) THEN
                  CALL WaveField_Interp_Setup3D( Time+WaveField%WaveTimeShift, posXY, WaveField%SrfGridParams, WaveField_m, ErrStat2, ErrMsg2 ); if (Failed()) return;
                  FV(:) = FV(:) + GridInterp3DVec( WaveField%PWaveVel0, WaveField_m ) * pos(3)
               END IF

            END IF ! Node is submerged

         ELSE ! Wheeler stretching - no need to check whether the node is above or below SWL

            ! Map the node z-position linearly from [-EffWtrDpth,m%WaveElev(j)] to [-EffWtrDpth,0]
            posPrime    = pos
            posPrime(3) = WaveField%EffWtrDpth*(WaveField%EffWtrDpth+pos(3))/(WaveField%EffWtrDpth+WaveElev)-WaveField%EffWtrDpth
            posPrime(3) = MIN( posPrime(3), 0.0_ReKi) ! Clamp z-position to zero. Needed when forceNodeInWater=.TRUE.

            ! Obtain the wave-field variables by interpolation with the mapped position.
            CALL WaveField_Interp_Setup4D( Time+WaveField%WaveTimeShift, posPrime, WaveField%GridDepth, WaveField%VolGridParams, WaveField_m, ErrStat2, ErrMsg2 ); if (Failed()) return;
            CALL WaveField_InterpVol( WaveField, WaveField_m, Time, ErrStat2, ErrMsg2, FV=FV ); if (Failed()) return;

         END IF

      ELSE ! Node is out of water - zero-out all wave dynamics

         nodeInWater = 0_IntKi
         FV(:)       = 0.0_SiKi

      END IF ! If node is in or out of water

   END IF ! If wave stretching is on or off

   ! Get dynamic current velocity
   IF (fetchDynCurrent .AND. WaveField%hasCurrField) THEN
      startNode = -1
      PosOffset = (/0.0_ReKi,0.0_ReKi,WaveField%EffWtrDpth/)
      posDummy(:,1) = pos
      ALLOCATE(FV_DC(3,1), STAT=ErrStat2); if (FailedMsg('Error allocating FV_DC')) return; 
      CALL IfW_FlowField_GetVelAcc(WaveField%CurrField, startNode, Time, posDummy, FV_DC, FA_DC, ErrStat2, ErrMsg2, PosOffset=PosOffset); if (Failed()) return;
      FV = FV + nodeInWater * FV_DC(:,1)
   END IF

contains
   logical function Failed()
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
   end function
   logical function FailedMsg(ErrMsg2)
      character(*), intent(in   ) :: ErrMsg2
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      FailedMsg = ErrStat >= AbortErrLev
   end function
END SUBROUTINE WaveField_GetNodeWaveVel


SUBROUTINE WaveField_GetNodeWaveVelAcc( WaveField, WaveField_m, Time, pos, forceNodeInWater, fetchDynCurrent, nodeInWater, FV, FA, ErrStat, ErrMsg )
   type(SeaSt_WaveFieldType),          intent(in   ) :: WaveField
   type(GridInterp_MiscVarType),       intent(inout) :: WaveField_m
   real(DbKi),                         intent(in   ) :: Time
   real(ReKi),                         intent(in   ) :: pos(3)
   logical,                            intent(in   ) :: forceNodeInWater
   logical,                            intent(in   ) :: fetchDynCurrent
   real(SiKi),                         intent(  out) :: FV(3)
   real(SiKi),                         intent(  out) :: FA(3)
   integer(IntKi),                     intent(  out) :: nodeInWater
   integer(IntKi),                     intent(  out) :: ErrStat ! Error status of the operation
   character(*),                       intent(  out) :: ErrMsg  ! Error message if errStat /= ErrID_None

   real(SiKi)                                        :: WaveElev
   real(ReKi)                                        :: posXY(2), posPrime(3), posXY0(3), PosOffset(3), posDummy(3,1)
   integer(IntKi)                                    :: startNode
   real(ReKi), allocatable                           :: FV_DC(:,:), FA_DC(:,:)
   character(*),                       parameter     :: RoutineName = 'WaveField_GetNodeWaveVelAcc'
   integer(IntKi)                                    :: errStat2
   character(ErrMsgLen)                              :: errMsg2

   ErrStat   = ErrID_None

   posXY    = pos(1:2)
   posXY0   = (/pos(1),pos(2),0.0_ReKi/)
   
   ! Wave elevation
   WaveElev  = WaveField_GetNodeTotalWaveElev( WaveField, WaveField_m, Time, pos, ErrStat2, ErrMsg2 ); if (Failed()) return;

   ! Check if point is below the seabed
   if (pos(3)<-WaveField%EffWtrDpth) then
      nodeInWater = 1_IntKi  ! Prevent problems with HydroDyn logic
      FV(:)       = 0.0_SiKi
      FA(:)       = 0.0_SiKi
      return
   end if

   IF (WaveField%WaveStMod == 0) THEN ! No wave stretching

      IF ( pos(3) <= 0.0_ReKi) THEN ! Node is at or below the SWL
         nodeInWater = 1_IntKi
         ! Use location to obtain interpolated values of kinematics
         CALL WaveField_Interp_Setup4D( Time+WaveField%WaveTimeShift, pos, WaveField%GridDepth, WaveField%VolGridParams, WaveField_m, ErrStat2, ErrMsg2 ); if (Failed()) return;
         CALL WaveField_InterpVol( WaveField, WaveField_m, Time, ErrStat2, ErrMsg2, FV=FV, FA=FA ); if (Failed()) return;
      ELSE ! Node is above the SWL
         nodeInWater = 0_IntKi
         FV(:)       = 0.0_SiKi
         FA(:)       = 0.0_SiKi
      END IF

   ELSE ! Wave stretching enabled

      IF ( (pos(3) <= WaveElev) .OR. forceNodeInWater ) THEN ! Node is submerged

         nodeInWater = 1_IntKi

         IF ( WaveField%WaveStMod < 3 ) THEN ! Vertical or extrapolated wave stretching

            IF ( pos(3) <= 0.0_SiKi) THEN ! Node is below the SWL - evaluate wave dynamics as usual

               ! Use location to obtain interpolated values of kinematics
               CALL WaveField_Interp_Setup4D( Time+WaveField%WaveTimeShift, pos, WaveField%GridDepth, WaveField%VolGridParams, WaveField_m, ErrStat2, ErrMsg2 ); if (Failed()) return;
               CALL WaveField_InterpVol( WaveField, WaveField_m, Time, ErrStat2, ErrMsg2, FV=FV, FA=FA ); if (Failed()) return;

            ELSE ! Node is above SWL - need wave stretching

               ! Vertical wave stretching
               CALL WaveField_Interp_Setup4D( Time, posXY0, WaveField%GridDepth, WaveField%VolGridParams, WaveField_m, ErrStat2, ErrMsg2 ); if (Failed()) return;
               CALL WaveField_InterpVol( WaveField, WaveField_m, Time, ErrStat2, ErrMsg2, FV=FV, FA=FA ); if (Failed()) return;

               ! Extrapolated wave stretching
               IF (WaveField%WaveStMod == 2) THEN
                  CALL WaveField_Interp_Setup3D( Time+WaveField%WaveTimeShift, posXY, WaveField%SrfGridParams, WaveField_m, ErrStat2, ErrMsg2 ); if (Failed()) return;
                  FV(:) = FV(:) + GridInterp3DVec( WaveField%PWaveVel0, WaveField_m ) * pos(3)
                  FA(:) = FA(:) + GridInterp3DVec( WaveField%PWaveAcc0, WaveField_m ) * pos(3)
               END IF

            END IF ! Node is submerged

         ELSE ! Wheeler stretching - no need to check whether the node is above or below SWL

            ! Map the node z-position linearly from [-EffWtrDpth,m%WaveElev(j)] to [-EffWtrDpth,0]
            posPrime    = pos
            posPrime(3) = WaveField%EffWtrDpth*(WaveField%EffWtrDpth+pos(3))/(WaveField%EffWtrDpth+WaveElev)-WaveField%EffWtrDpth
            posPrime(3) = MIN( posPrime(3), 0.0_ReKi) ! Clamp z-position to zero. Needed when forceNodeInWater=.TRUE.

            ! Obtain the wave-field variables by interpolation with the mapped position.
            CALL WaveField_Interp_Setup4D( Time+WaveField%WaveTimeShift, posPrime, WaveField%GridDepth, WaveField%VolGridParams, WaveField_m, ErrStat2, ErrMsg2 ); if (Failed()) return;
            CALL WaveField_InterpVol( WaveField, WaveField_m, Time, ErrStat2, ErrMsg2, FV=FV, FA=FA ); if (Failed()) return;
         END IF

      ELSE ! Node is out of water - zero-out all wave dynamics

         nodeInWater = 0_IntKi
         FV(:)       = 0.0_SiKi
         FA(:)       = 0.0_SiKi
      END IF ! If node is in or out of water

   END IF ! If wave stretching is on or off
   
   ! Get dynamic current velocity and acceleration
   IF (fetchDynCurrent .AND. WaveField%hasCurrField) THEN
      startNode = -1
      PosOffset = (/0.0_ReKi,0.0_ReKi,WaveField%EffWtrDpth/)
      posDummy(:,1) = pos
      ALLOCATE(FV_DC(3,1), STAT=ErrStat2); if (FailedMsg('Error allocating FV_DC')) return;    
      ALLOCATE(FA_DC(3,1), STAT=ErrStat2); if (FailedMsg('Error allocating FA_DC')) return;
      CALL IfW_FlowField_GetVelAcc(WaveField%CurrField, startNode, Time, posDummy, FV_DC, FA_DC, ErrStat2, ErrMsg2, PosOffset=PosOffset); if (Failed()) return;
      FV = FV + nodeInWater * FV_DC(:,1)
      FA = FA + nodeInWater * FA_DC(:,1)
   END IF

contains
   logical function Failed()
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
   end function
   logical function FailedMsg(ErrMsg2)
      character(*), intent(in   ) :: ErrMsg2
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      FailedMsg = ErrStat >= AbortErrLev
   end function
END SUBROUTINE WaveField_GetNodeWaveVelAcc


SUBROUTINE WaveField_GetWaveKin( WaveField, WaveField_m, Time, pos, forceNodeInWater, fetchDynCurrent, nodeInWater, WaveElev1, WaveElev2, WaveElev, FDynP, FV, FA, FAMCF, ErrStat, ErrMsg )
   type(SeaSt_WaveFieldType),          intent(in   ) :: WaveField
   type(GridInterp_MiscVarType),       intent(inout) :: WaveField_m
   real(DbKi),                         intent(in   ) :: Time
   real(ReKi),                         intent(in   ) :: pos(:,:)
   logical,                            intent(in   ) :: forceNodeInWater
   logical,                            intent(in   ) :: fetchDynCurrent
   real(SiKi),                         intent(  out) :: WaveElev1(:)
   real(SiKi),                         intent(  out) :: WaveElev2(:)
   real(SiKi),                         intent(  out) :: WaveElev(:)
   real(ReKi),                         intent(  out) :: FV(:,:)
   real(ReKi),                         intent(  out) :: FA(:,:)
   real(ReKi),                         intent(  out) :: FAMCF(:,:)
   real(ReKi),                         intent(  out) :: FDynP(:)
   integer(IntKi),                     intent(  out) :: nodeInWater(:)
   integer(IntKi),                     intent(  out) :: ErrStat ! Error status of the operation
   character(*),                       intent(  out) :: ErrMsg  ! Error message if errStat /= ErrID_None

   character(*),                       parameter     :: RoutineName = 'WaveField_GetWaveKin'
   integer(IntKi)                                    :: errStat2
   character(ErrMsgLen)                              :: errMsg2

   integer(IntKi)                                    :: NumPoints, i, startNode
   real(SiKi)                                        :: FDynP_node, FV_node(3), FA_node(3), FAMCF_node(3)
   real(ReKi)                                        :: PosOffset(3)
   real(ReKi),                         allocatable   :: FV_DC(:,:), FA_DC(:,:)

   ErrStat   = ErrID_None

   NumPoints = size(pos, dim=2)
   DO i = 1, NumPoints
      CALL WaveField_GetNodeWaveKin( WaveField, WaveField_m, Time, pos(:,i), forceNodeInWater, .FALSE., nodeInWater(i), WaveElev1(i), WaveElev2(i), WaveElev(i), FDynP_node, FV_node, FA_node, FAMCF_node, ErrStat2, ErrMsg2 )
      if (Failed()) return;
      FDynP(i) = REAL(FDynP_node,ReKi)
      FV(:, i) = REAL(FV_node,   ReKi)
      FA(:, i) = REAL(FA_node,   ReKi)
      IF ( WaveField_HasMCF(WaveField) ) THEN
         FAMCF(:,i) = REAL(FAMCF_node,ReKi)
      END IF
   END DO

   ! If dynamic current field from IfW is present, get velocity and acceleration contributions
   IF (fetchDynCurrent .AND. WaveField%hasCurrField) THEN
      startNode = -1
      PosOffset = (/0.0_ReKi,0.0_ReKi,WaveField%EffWtrDpth/)
      ALLOCATE(FV_DC( 3, NumPoints ), STAT=ErrStat2); if (FailedMsg('Error allocating FV_DC')) return;  
      ALLOCATE(FA_DC( 3, NumPoints ), STAT=ErrStat2); if (FailedMsg('Error allocating FA_DC')) return;
      CALL IfW_FlowField_GetVelAcc(WaveField%CurrField, startNode, Time, pos, FV_DC, FA_DC, ErrStat2, ErrMsg2, PosOffset=PosOffset); if (Failed()) return;

      ! Add contributions from IfW current field if node is in water
      DO i = 1, NumPoints
         FV(:,i) = FV(:,i) + nodeInWater(i) * FV_DC(:,i)
         FA(:,i) = FA(:,i) + nodeInWater(i) * FA_DC(:,i)
      END DO

   END IF

contains
   logical function Failed()
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
   end function
   logical function FailedMsg(ErrMsgTmp)
      character(*), intent(in   ) :: ErrMsgTmp
      call SetErrStat( ErrStat2, ErrMsgTmp, ErrStat, ErrMsg, RoutineName )
      FailedMsg = ErrStat >= AbortErrLev
   end function
end subroutine WaveField_GetWaveKin


! This subroutine is intended for AeroDyn when modeling MHK turbines
SUBROUTINE WaveField_GetWaveVelAcc_AD( WaveField, WaveField_m, StartNode, Time, pos, FV, FA, ErrStat, ErrMsg, BoxExceedAllow )
   type(SeaSt_WaveFieldType),          intent(in   ) :: WaveField
   type(GridInterp_MiscVarType),       intent(inout) :: WaveField_m
   integer(IntKi),                     intent(in   ) :: StartNode
   real(DbKi),                         intent(in   ) :: Time
   real(ReKi),                         intent(in   ) :: pos(:,:) ! z=0 at MSL
   real(ReKi),                         intent(  out) :: FV(:,:)
   real(ReKi),       allocatable,      intent(inout) :: FA(:,:)
   integer(IntKi),                     intent(  out) :: ErrStat  ! Error status of the operation
   character(*),                       intent(  out) :: ErrMsg   ! Error message if errStat /= ErrID_None
   logical,          optional,         intent(in   ) :: BoxExceedAllow
   integer(IntKi),   allocatable                     :: nodeInWater(:)
   integer(IntKi)                                    :: NumPoints, i
   real(SiKi)                                        :: FV_node(3), FA_node(3)
   real(ReKi)                                        :: PosOffset(3), MSL2SWL, WtrDpth
   real(ReKi),       allocatable                     :: FV_DC(:,:), FA_DC(:,:)
   logical                                           :: getAcc
   character(*),     parameter                       :: RoutineName = 'WaveField_GetWaveVelAcc_AD'
   integer(IntKi)                                    :: errStat2
   character(ErrMsgLen)                              :: errMsg2

   ErrStat   = ErrID_None

   MSL2SWL   = WaveField%MSL2SWL
   WtrDpth   = WaveField%EffWtrDpth - MSL2SWL
   getAcc    = ALLOCATED(FA)
   NumPoints = size(pos, dim=2)

   ALLOCATE( nodeInWater(NumPoints), STAT=ErrStat2); if (FailedMsg('Error allocating nodeInWater')) return;

   ! Note: SeaState wavefield grid has z=0 on the SWL
   IF (getAcc) THEN
      DO i = 1, NumPoints
         CALL WaveField_GetNodeWaveVelAcc( WaveField, WaveField_m, Time, pos(:,i)-[0.0_ReKi,0.0_ReKi,real(MSL2SWL, ReKi)], .FALSE., .FALSE., nodeInWater(i), FV_node, FA_node, ErrStat2, ErrMsg2 ); if (Failed()) return;
         FV(:, i) = REAL(FV_node,   ReKi)
         FA(:, i) = REAL(FA_node,   ReKi)
      END DO
   ELSE
      DO i = 1, NumPoints
         CALL WaveField_GetNodeWaveVel( WaveField, WaveField_m, Time, pos(:,i)-[0.0_ReKi,0.0_ReKi,real(MSL2SWL, ReKi)], .FALSE., .FALSE., nodeInWater(i), FV_node, ErrStat2, ErrMsg2 ); if (Failed()) return;
         FV(:, i) = REAL(FV_node,   ReKi)
      END DO
   END IF

   ! If dynamic current field from IfW is present, get velocity and acceleration contributions
   IF (WaveField%hasCurrField) THEN
      PosOffset = (/0.0_ReKi,0.0_ReKi,WtrDpth/) ! IfW FlowField grid effectively has z=0 on the seabed
      ALLOCATE(FV_DC( 3, NumPoints ), STAT=ErrStat2); if (FailedMsg('Error allocating FV_DC')) return;
      IF (getAcc) THEN
         ALLOCATE(FA_DC( 3, NumPoints ), STAT=ErrStat2); if (FailedMsg('Error allocating FA_DC')) return;
      END IF

      IF (PRESENT(BoxExceedAllow)) THEN
         CALL IfW_FlowField_GetVelAcc(WaveField%CurrField, StartNode, Time, pos, FV_DC, FA_DC, ErrStat2, ErrMsg2, BoxExceedAllow=BoxExceedAllow, PosOffset=PosOffset); if (Failed()) return
      ELSE
         CALL IfW_FlowField_GetVelAcc(WaveField%CurrField, StartNode, Time, pos, FV_DC, FA_DC, ErrStat2, ErrMsg2, PosOffset=PosOffset); if (Failed()) return
      END IF

      ! Add contributions from IfW current field if node is in water
      DO i = 1, NumPoints
         FV(:,i) = FV(:,i) + nodeInWater(i) * FV_DC(:,i)
      END DO
      IF (getAcc) THEN
         DO i = 1, NumPoints
            FA(:,i) = FA(:,i) + nodeInWater(i) * FA_DC(:,i)
         END DO
      END IF
   END IF

contains
   logical function Failed()
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
   end function
   logical function FailedMsg(ErrMsg2)
      character(*), intent(in   ) :: ErrMsg2
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      FailedMsg = ErrStat >= AbortErrLev
   end function
END SUBROUTINE WaveField_GetWaveVelAcc_AD

!-------------------- Subroutine for wave field velocity only --------------------!
SUBROUTINE WaveField_GetMeanDynSurfCurr( WaveField, WaveTMax, WaveDT, CurrVxi0, CurrVyi0, ErrStat, ErrMsg )
   type(SeaSt_WaveFieldType),          intent(in   ) :: WaveField

   real(DbKi),                         intent(in   ) :: WaveTMax
   real(DbKi),                         intent(in   ) :: WaveDT

   real(SiKi),                         intent(  out) :: CurrVxi0
   real(SiKi),                         intent(  out) :: CurrVyi0
   integer(IntKi),                     intent(  out) :: ErrStat ! Error status of the operation
   character(*),                       intent(  out) :: ErrMsg  ! Error message if errStat /= ErrID_None

   real(ReKi)                                        :: pos(3,1), PosOffset(3)
   real(ReKi), allocatable                           :: FV_DC(:,:), FA_DC(:,:)
   integer(IntKi)                                    :: startNode
   integer(IntKi)                                    :: step
   real(DbKi)                                        :: time
   character(*),                       parameter     :: RoutineName = 'WaveField_GetMeanDynSurfCurr'
   integer(IntKi)                                    :: errStat2
   character(ErrMsgLen)                              :: errMsg2

   ErrStat   = ErrID_None

   CurrVxi0 = 0.0_SiKi
   CurrVyi0 = 0.0_SiKi

   ! Get dynamic current velocity
   IF ( WaveField%hasCurrField ) THEN

      pos       = 0.0_ReKi
      step      = 0_IntKi
      time      = 0.0_DbKi
      startNode = -1
      PosOffset = (/0.0_ReKi,0.0_ReKi,WaveField%EffWtrDpth/)
      ALLOCATE(FV_DC(3,1), STAT=ErrStat2); if (FailedMsg('Error allocating FV_DC')) return;

      DO WHILE ( time <= WaveTMax)
         CALL IfW_FlowField_GetVelAcc(WaveField%CurrField, startNode, Time, pos, FV_DC, FA_DC, ErrStat2, ErrMsg2, PosOffset=PosOffset); if (Failed()) return;
         CurrVxi0 = CurrVxi0 + FV_DC(1,1)
         CurrVyi0 = CurrVyi0 + FV_DC(2,1)
         step = step + 1
         time = time + WaveDT
      END DO
      CurrVxi0 = CurrVxi0 / REAL(step,SiKi)
      CurrVyi0 = CurrVyi0 / REAL(step,SiKi)

   END IF

contains
   logical function Failed()
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
   end function
   logical function FailedMsg(ErrMsg2)
      character(*), intent(in   ) :: ErrMsg2
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      FailedMsg = ErrStat >= AbortErrLev
   end function
END SUBROUTINE WaveField_GetMeanDynSurfCurr

!----------------------------------------------------------------------------------------------------
! On-demand wave-kinematics block partitioning (WvKinBlockMod=True)
!----------------------------------------------------------------------------------------------------

!> True when a MacCamy-Fuchs scaled acceleration field exists — either the eager full-domain array
!! (WvKinBlockMod=False) or per-block MCF data (WvKinBlockMod=True with MCFD>0, where the full-domain array
!! is never allocated).
LOGICAL FUNCTION WaveField_HasMCF( WaveField )
   type(SeaSt_WaveFieldType), intent(in   ) :: WaveField
   WaveField_HasMCF = ALLOCATED(WaveField%WaveAccMCF) .OR. &
                      ( WaveField%WvKinBlockMod .AND. WaveField%MCFD > 0.0_SiKi )
END FUNCTION WaveField_HasMCF


!> Set up the XY block layout of the wave-kinematics volume grid (WvKinBlockMod=True). Blocks partition the
!! grid CELLS; each block additionally stores the grid-point planes its interpolation stencils reach into
!! (the volume interpolation uses a 4-point stencil per dimension, base cell c touching points c-1..c+2),
!! so any query whose base cell lies in a block is served entirely from that block. No block data is
!! allocated here — population happens on first access.
SUBROUTINE WaveField_BlockStore_Init( WaveField, ErrStat, ErrMsg )
   type(SeaSt_WaveFieldType),    intent(in   ) :: WaveField   ! block store reached through the pointer component (mutable)
   integer(IntKi),               intent(  out) :: ErrStat
   character(*),                 intent(  out) :: ErrMsg

   integer(IntKi)                              :: NX, NY, NZ, ib, jb, kb, k, c0, c1, iP0, iP1
   real(ReKi)                                  :: mean_dz
   integer(IntKi)                              :: ErrStat2
   character(*), parameter                     :: RoutineName = 'WaveField_BlockStore_Init'

   ErrStat = ErrID_None
   ErrMsg  = ""

   IF ( .NOT. ASSOCIATED(WaveField%BlockStore) ) THEN
      CALL SetErrStat( ErrID_Fatal, 'The wave-kinematics block store has not been wired up.', ErrStat, ErrMsg, RoutineName )
      RETURN
   END IF

   ASSOCIATE ( Store => WaveField%BlockStore )

      NX = WaveField%VolGridParams%n(2)
      NY = WaveField%VolGridParams%n(3)
      NZ = WaveField%VolGridParams%n(4)

      ! Target edge length snapped to whole grid cells, at least 8 cells (overlap economy), at most the whole grid
      IF ( NX > 1_IntKi ) THEN
         Store%BlkCellsX = MIN( MAX( 8_IntKi, NINT( WaveField%WvKinBlockSize / WaveField%VolGridParams%delta(2) ) ), NX-1_IntKi )
         Store%nBlkX     = ( (NX-1_IntKi) + Store%BlkCellsX - 1_IntKi ) / Store%BlkCellsX
      ELSE
         Store%BlkCellsX = 1_IntKi
         Store%nBlkX     = 1_IntKi
      END IF
      IF ( NY > 1_IntKi ) THEN
         Store%BlkCellsY = MIN( MAX( 8_IntKi, NINT( WaveField%WvKinBlockSize / WaveField%VolGridParams%delta(3) ) ), NY-1_IntKi )
         Store%nBlkY     = ( (NY-1_IntKi) + Store%BlkCellsY - 1_IntKi ) / Store%BlkCellsY
      ELSE
         Store%BlkCellsY = 1_IntKi
         Store%nBlkY     = 1_IntKi
      END IF
      ! z: the grid is cosine-distributed (non-uniform), so snap using the mean spacing; this keeps
      ! z blocks as cubes only in the mean-spacing sense (exact metre-cubes where spacing is uniform).
      IF ( NZ > 1_IntKi ) THEN
         mean_dz = ABS( Store%zGrid(SIZE(Store%zGrid)) - Store%zGrid(1) ) / REAL( MAX(1_IntKi,NZ-1_IntKi), ReKi )
         IF ( mean_dz > 0.0_ReKi ) THEN
            Store%BlkCellsZ = MIN( MAX( 8_IntKi, NINT( WaveField%WvKinBlockSize / mean_dz ) ), NZ-1_IntKi )
            Store%nBlkZ     = ( (NZ-1_IntKi) + Store%BlkCellsZ - 1_IntKi ) / Store%BlkCellsZ
         ELSE
            Store%BlkCellsZ = NZ-1_IntKi
            Store%nBlkZ     = 1_IntKi
         END IF
      ELSE
         Store%BlkCellsZ = 1_IntKi
         Store%nBlkZ     = 1_IntKi
      END IF

      ALLOCATE ( Store%Blocks( Store%nBlkX * Store%nBlkY * Store%nBlkZ ), STAT=ErrStat2 )
      IF ( ErrStat2 /= 0 ) THEN
         CALL SetErrStat( ErrID_Fatal, 'Error allocating the wave-kinematics block array.', ErrStat, ErrMsg, RoutineName )
         RETURN
      END IF

      DO kb = 1, Store%nBlkZ
         DO jb = 1, Store%nBlkY
            DO ib = 1, Store%nBlkX
               k = ((kb-1)*Store%nBlkY + (jb-1))*Store%nBlkX + ib
               ! x extent: cells [c0..c1] (0-based), stored points [c0-1 .. c1+2] clamped to the grid
               IF ( NX > 1_IntKi ) THEN
                  c0  = (ib-1_IntKi)*Store%BlkCellsX
                  c1  = MIN( ib*Store%BlkCellsX - 1_IntKi, NX - 2_IntKi )
                  iP0 = MAX( 0_IntKi, c0 - 1_IntKi )
                  iP1 = MIN( NX - 1_IntKi, c1 + 2_IntKi )
               ELSE
                  iP0 = 0_IntKi
                  iP1 = 0_IntKi
               END IF
               Store%Blocks(k)%iPtX0 = iP0 + 1_IntKi   ! stored 1-based
               Store%Blocks(k)%nPtX  = iP1 - iP0 + 1_IntKi
               ! y extent
               IF ( NY > 1_IntKi ) THEN
                  c0  = (jb-1_IntKi)*Store%BlkCellsY
                  c1  = MIN( jb*Store%BlkCellsY - 1_IntKi, NY - 2_IntKi )
                  iP0 = MAX( 0_IntKi, c0 - 1_IntKi )
                  iP1 = MIN( NY - 1_IntKi, c1 + 2_IntKi )
               ELSE
                  iP0 = 0_IntKi
                  iP1 = 0_IntKi
               END IF
               Store%Blocks(k)%iPtY0 = iP0 + 1_IntKi
               Store%Blocks(k)%nPtY  = iP1 - iP0 + 1_IntKi
               ! z extent
               IF ( NZ > 1_IntKi ) THEN
                  c0  = (kb-1_IntKi)*Store%BlkCellsZ
                  c1  = MIN( kb*Store%BlkCellsZ - 1_IntKi, NZ - 2_IntKi )
                  iP0 = MAX( 0_IntKi, c0 - 1_IntKi )
                  iP1 = MIN( NZ - 1_IntKi, c1 + 2_IntKi )
               ELSE
                  iP0 = 0_IntKi
                  iP1 = 0_IntKi
               END IF
               Store%Blocks(k)%iPtZ0 = iP0 + 1_IntKi
               Store%Blocks(k)%nPtZ  = iP1 - iP0 + 1_IntKi
            END DO
         END DO
      END DO

      CALL WrScr ( ' SeaState wave-kinematics on-demand blocks: '//TRIM(Num2LStr(Store%nBlkX))//' x '// &
                   TRIM(Num2LStr(Store%nBlkY))//' x '//TRIM(Num2LStr(Store%nBlkZ))//' blocks of '// &
                   TRIM(Num2LStr(Store%BlkCellsX))//' x '//TRIM(Num2LStr(Store%BlkCellsY))//' x '// &
                   TRIM(Num2LStr(Store%BlkCellsZ))//' grid cells (+stencil overlap); populated on first access.' )

   END ASSOCIATE

END SUBROUTINE WaveField_BlockStore_Init


!> Locate the wave block containing the interpolation stencil last set up in WaveField_m, populate it on
!! first access (first- and, if enabled, second-order kinematics via the shared column kernels — the same
!! code path as the full-domain fill, so block contents are bit-identical to the mode-0 arrays), and stamp
!! its last-access time. The check-allocate-fill-publish sequence is serialized for thread safety
!! (FAST.Farm shares turbine 1's wave field with farm-level MoorDyn).
SUBROUTINE WaveField_EnsureBlock( WaveField, WaveField_m, Time, iBlk, ErrStat, ErrMsg )
   type(SeaSt_WaveFieldType),    intent(in   ) :: WaveField
   type(GridInterp_MiscVarType), intent(in   ) :: WaveField_m
   real(DbKi),                   intent(in   ) :: Time
   integer(IntKi),               intent(  out) :: iBlk
   integer(IntKi),               intent(  out) :: ErrStat
   character(*),                 intent(  out) :: ErrMsg

   integer(IntKi)                              :: ib, jb, kb, nResident
   real(ReKi)                                  :: MBytes
   integer(IntKi)                              :: ErrStat2
   character(ErrMsgLen)                        :: ErrMsg2
   character(*), parameter                     :: RoutineName = 'WaveField_EnsureBlock'

   ErrStat = ErrID_None
   ErrMsg  = ""

   ASSOCIATE ( Store => WaveField%BlockStore )

      ! Locate the block from the stencil's base cell: Indx(2,dim) is the 0-based lower cell index of the
      ! query, already clamped in range by GridInterpSetup4D. All four stencil points per dimension are
      ! then guaranteed to lie within the block's stored point range (see WaveField_BlockStore_Init).
      ib   = MIN( WaveField_m%Indx(2,2) / Store%BlkCellsX, Store%nBlkX - 1_IntKi ) + 1_IntKi
      jb   = MIN( WaveField_m%Indx(2,3) / Store%BlkCellsY, Store%nBlkY - 1_IntKi ) + 1_IntKi
      kb   = MIN( WaveField_m%Indx(2,4) / Store%BlkCellsZ, Store%nBlkZ - 1_IntKi ) + 1_IntKi
      iBlk = ((kb-1_IntKi)*Store%nBlkY + (jb-1_IntKi))*Store%nBlkX + ib

      IF ( .NOT. Store%Blocks(iBlk)%Populated ) THEN
         !$OMP CRITICAL(SeaSt_BlockPop)
         IF ( .NOT. Store%Blocks(iBlk)%Populated ) THEN

            ASSOCIATE ( Blk => Store%Blocks(iBlk) )

               ALLOCATE ( Blk%WaveDynP(0:WaveField%NStepWave, Blk%nPtX, Blk%nPtY, Blk%nPtZ   ), &
                          Blk%WaveVel (0:WaveField%NStepWave, Blk%nPtX, Blk%nPtY, Blk%nPtZ, 3), &
                          Blk%WaveAcc (0:WaveField%NStepWave, Blk%nPtX, Blk%nPtY, Blk%nPtZ, 3), STAT=ErrStat2 )
               IF ( ErrStat2 == 0 .AND. WaveField%MCFD > 0.0_SiKi ) &
                  ALLOCATE ( Blk%WaveAccMCF(0:WaveField%NStepWave, Blk%nPtX, Blk%nPtY, Blk%nPtZ, 3), STAT=ErrStat2 )
               IF ( ErrStat2 /= 0 ) THEN
                  CALL SetErrStat( ErrID_Fatal, 'Error allocating the arrays of wave block ('// &
                                   TRIM(Num2LStr(ib))//','//TRIM(Num2LStr(jb))//').', ErrStat, ErrMsg, RoutineName )
               END IF

               IF ( ErrStat < AbortErrLev ) THEN
                  IF ( WaveField%MCFD > 0.0_SiKi ) THEN
                     CALL WaveKinKernel_ComputeColumns ( WaveField, Store, Blk%iPtX0, Blk%nPtX, Blk%iPtY0, Blk%nPtY, &
                                                         Blk%iPtZ0, Blk%nPtZ, &
                                                         Blk%WaveDynP, Blk%WaveVel, Blk%WaveAcc, ErrStat2, ErrMsg2, &
                                                         WaveAccMCF=Blk%WaveAccMCF )
                  ELSE
                     CALL WaveKinKernel_ComputeColumns ( WaveField, Store, Blk%iPtX0, Blk%nPtX, Blk%iPtY0, Blk%nPtY, &
                                                         Blk%iPtZ0, Blk%nPtZ, &
                                                         Blk%WaveDynP, Blk%WaveVel, Blk%WaveAcc, ErrStat2, ErrMsg2 )
                  END IF
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               END IF

               IF ( ErrStat < AbortErrLev .AND. ( Store%SecondOrderDiff .OR. Store%SecondOrderSum ) ) THEN
                  CALL WaveKinKernel_AddSecondOrderColumns ( WaveField, Store, Blk%iPtX0, Blk%nPtX, Blk%iPtY0, Blk%nPtY, &
                                                             Blk%iPtZ0, Blk%nPtZ, &
                                                             Blk%WaveDynP, Blk%WaveVel, Blk%WaveAcc, ErrStat2, ErrMsg2 )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               END IF

               IF ( ErrStat < AbortErrLev ) THEN
                  Store%nPopulated    = Store%nPopulated + 1_IntKi
                  nResident           = Store%nPopulated - Store%nEvicted
                  Store%nPeakResident = MAX( Store%nPeakResident, nResident )
                  MBytes = REAL(WaveField%NStepWave+1,ReKi) * Blk%nPtX * Blk%nPtY * Blk%nPtZ * 4.0_ReKi * &
                           MERGE( 10.0_ReKi, 7.0_ReKi, WaveField%MCFD > 0.0_SiKi ) / 1.0E6_ReKi
                  CALL WrScr ( ' SeaState: populated wave block ('//TRIM(Num2LStr(ib))//','//TRIM(Num2LStr(jb))//','// &
                               TRIM(Num2LStr(kb))// &
                               ') covering x('//TRIM(Num2LStr(Blk%iPtX0))//':'//TRIM(Num2LStr(Blk%iPtX0+Blk%nPtX-1))// &
                               '), y('//TRIM(Num2LStr(Blk%iPtY0))//':'//TRIM(Num2LStr(Blk%iPtY0+Blk%nPtY-1))// &
                               '), z('//TRIM(Num2LStr(Blk%iPtZ0))//':'//TRIM(Num2LStr(Blk%iPtZ0+Blk%nPtZ-1))// &
                               ') of the volume grid at t='//TRIM(Num2LStr(REAL(Time,ReKi)))//' s ('// &
                               TRIM(Num2LStr(MBytes))//' MB; '//TRIM(Num2LStr(nResident))//' block(s) resident).' )
                  Blk%LastAccess = Time    ! stamp before any sweep so the fresh block is never a victim
                  Blk%EverPopulated = .TRUE.
                  Blk%Populated = .TRUE.   ! publish last
                  ! This population raised resident memory — free any now-idle blocks while the lock is held.
                  IF ( WaveField%WvKinBlockFreeT > 0.0_DbKi ) CALL WaveField_SweepBlocks( WaveField, Time )
               END IF

            END ASSOCIATE

         END IF
         !$OMP END CRITICAL(SeaSt_BlockPop)
         IF ( ErrStat >= AbortErrLev ) RETURN
      END IF

      Store%Blocks(iBlk)%LastAccess = MAX( Store%Blocks(iBlk)%LastAccess, Time )

      ! Throttled idle-block eviction on any access (the population path above already sweeps under the
      ! lock). The unlocked read of LastSweep is a benign race — re-checked inside the critical section.
      IF ( WaveField%WvKinBlockFreeT > 0.0_DbKi .AND. &
           Time - Store%LastSweep > MAX( WaveField%WvKinBlockFreeT, 60.0_DbKi ) ) THEN
         !$OMP CRITICAL(SeaSt_BlockPop)
         IF ( Time - Store%LastSweep > MAX( WaveField%WvKinBlockFreeT, 60.0_DbKi ) ) THEN
            CALL WaveField_SweepBlocks( WaveField, Time )
         END IF
         !$OMP END CRITICAL(SeaSt_BlockPop)
      END IF

   END ASSOCIATE

END SUBROUTINE WaveField_EnsureBlock


!> Free the arrays of every populated block that has gone longer than WvKinBlockFreeT seconds of
!! simulation time without an access, so its memory is reclaimed until the query returns. Blocks are
!! deterministically regenerated by the shared kernels on the next access, so eviction never changes
!! results. The caller MUST hold the SeaSt_BlockPop critical section — this routine mutates the shared
!! block store and must not race population; it therefore contains no critical region of its own.
!! Stamps Store%LastSweep, counts evictions, and logs one line per freed block. Callers pass only when
!! WvKinBlockFreeT > 0 (a non-positive value disables freeing entirely).
!!
!! Serial-access assumption: today every wave-field query is serial (single-turbine runs are
!! single-threaded and FAST.Farm steps the shared field serially — there is no !$OMP PARALLEL region
!! around the accessors), and WaveField_EnsureBlock never evicts the block it just ensured (its own
!! LastAccess is stamped to the current Time first), so the unlocked block-array reads in
!! WaveField_InterpVol that follow an EnsureBlock call always see live arrays. If accessors are ever
!! driven from a genuinely concurrent (OMP-parallel) region over a shared field, those reads would
!! need to pin their block against eviction (e.g. a per-block reader count checked here) before this
!! sweep could run concurrently with them.
SUBROUTINE WaveField_SweepBlocks( WaveField, Time )
   type(SeaSt_WaveFieldType),    intent(in   ) :: WaveField
   real(DbKi),                   intent(in   ) :: Time

   integer(IntKi)                              :: iBlk, nResident

   ASSOCIATE ( Store => WaveField%BlockStore )

      DO iBlk = 1, SIZE(Store%Blocks)
         ASSOCIATE ( Blk => Store%Blocks(iBlk) )
            IF ( Blk%Populated .AND. ( Time - Blk%LastAccess > WaveField%WvKinBlockFreeT ) ) THEN
               Blk%Populated = .FALSE.   ! unpublish first, before freeing the arrays
               IF ( ALLOCATED(Blk%WaveDynP)   ) DEALLOCATE( Blk%WaveDynP )
               IF ( ALLOCATED(Blk%WaveVel)    ) DEALLOCATE( Blk%WaveVel )
               IF ( ALLOCATED(Blk%WaveAcc)    ) DEALLOCATE( Blk%WaveAcc )
               IF ( ALLOCATED(Blk%WaveAccMCF) ) DEALLOCATE( Blk%WaveAccMCF )
               Store%nEvicted = Store%nEvicted + 1_IntKi
               nResident      = Store%nPopulated - Store%nEvicted
               CALL WrScr ( ' SeaState: freed idle wave block '//TRIM(Num2LStr(iBlk))//' (idle '// &
                            TRIM(Num2LStr(REAL(Time - Blk%LastAccess,ReKi)))//' s > WvKinBlockFreeT='// &
                            TRIM(Num2LStr(REAL(WaveField%WvKinBlockFreeT,ReKi)))//' s) at t='// &
                            TRIM(Num2LStr(REAL(Time,ReKi)))//' s ('//TRIM(Num2LStr(nResident))//' block(s) resident).' )
            END IF
         END ASSOCIATE
      END DO

      Store%LastSweep = Time

   END ASSOCIATE

END SUBROUTINE WaveField_SweepBlocks


!----------------------------------------------------------------------------------------------------------------------------------
!> Write the on-demand block partition as a legacy-VTK rectilinear grid: one cell per block,
!> one CELL_DATA scalar "BlockLife" (-1 = never populated, 0 = evicted, (0,1] = active,
!> normalized time remaining before eviction; pinned at 1 when eviction is disabled).
!> Silent no-op unless WvKinBlockMod=True with an allocated block store. Errors are warnings only.
SUBROUTINE WaveField_WriteBlockVTK ( Time, WaveField, OutRootName, FrameNo, TWidth, ErrStat, ErrMsg )
   REAL(DbKi),                INTENT(IN   ) :: Time
   TYPE(SeaSt_WaveFieldType), INTENT(IN   ) :: WaveField
   CHARACTER(*),              INTENT(IN   ) :: OutRootName
   INTEGER(IntKi),            INTENT(IN   ) :: FrameNo
   INTEGER(IntKi),            INTENT(IN   ) :: TWidth
   INTEGER(IntKi),            INTENT(  OUT) :: ErrStat
   CHARACTER(*),              INTENT(  OUT) :: ErrMsg

   TYPE(SeaSt_WaveBlockStoreType), POINTER :: Store
   TYPE(VTK_Misc)                          :: mvtk
   REAL(ReKi), ALLOCATABLE                 :: xE(:), yE(:), zE(:), Life(:)
   INTEGER(IntKi)                          :: k, nBlk, ErrStatTmp
   CHARACTER(64)                           :: Tstr
   CHARACTER(1024)                         :: FileName
   CHARACTER(64)                           :: Descr
   CHARACTER(*), PARAMETER                 :: RoutineName = 'WaveField_WriteBlockVTK'

   ErrStat = ErrID_None
   ErrMsg  = ''

   IF ( .NOT. WaveField%WvKinBlockMod ) RETURN
   IF ( .NOT. ASSOCIATED(WaveField%BlockStore) ) RETURN
   Store => WaveField%BlockStore
   IF ( .NOT. ALLOCATED(Store%Blocks) ) RETURN

   nBlk = Store%nBlkX * Store%nBlkY * Store%nBlkZ

   ALLOCATE ( xE(Store%nBlkX+1), yE(Store%nBlkY+1), zE(Store%nBlkZ+1), Life(nBlk), STAT=ErrStatTmp )
   IF ( ErrStatTmp /= 0 ) THEN
      CALL SetErrStat( ErrID_Warn, 'Could not allocate block VTK work arrays.', ErrStat, ErrMsg, RoutineName )
      RETURN
   END IF

   CALL BlockEdgeCoords( Store%xGrid, Store%nBlkX, Store%BlkCellsX, xE )
   CALL BlockEdgeCoords( Store%yGrid, Store%nBlkY, Store%BlkCellsY, yE )
   CALL BlockEdgeCoords( Store%zGrid, Store%nBlkZ, Store%BlkCellsZ, zE )

   ! Blocks(k) ordering (x fastest, then y, then z) matches VTK rectilinear cell ordering.
   DO k = 1, nBlk
      IF ( .NOT. Store%Blocks(k)%EverPopulated ) THEN
         Life(k) = -1.0_ReKi
      ELSE IF ( .NOT. Store%Blocks(k)%Populated ) THEN
         Life(k) = 0.0_ReKi
      ELSE IF ( WaveField%WvKinBlockFreeT > 0.0_DbKi ) THEN
         ! floor keeps an overdue-but-unswept active block from displaying as evicted
         Life(k) = REAL( MAX( ( WaveField%WvKinBlockFreeT - (Time - Store%Blocks(k)%LastAccess) ) &
                              / WaveField%WvKinBlockFreeT, 0.001_DbKi ), ReKi )
      ELSE
         Life(k) = 1.0_ReKi   ! eviction disabled: alive forever
      END IF
   END DO

   WRITE (Tstr,'(I'//TRIM(Num2LStr(TWidth))//'.'//TRIM(Num2LStr(TWidth))//')') FrameNo
   FileName = TRIM(OutRootName)//'.SeaSt.WaveBlocks.'//TRIM(Tstr)//'.vtk'
   WRITE (Descr,'(A,F0.4,A)') 'SeaState wave-kinematics blocks, t=', Time, ' s'

   CALL vtk_misc_init( mvtk )
   IF ( .NOT. vtk_new_ascii_file( TRIM(FileName), TRIM(Descr), mvtk ) ) THEN
      CALL SetErrStat( ErrID_Warn, 'Could not open block VTK file "'//TRIM(FileName)//'".', ErrStat, ErrMsg, RoutineName )
      RETURN
   END IF
   CALL vtk_dataset_rectilinear( xE, yE, zE, mvtk )
   mvtk%nData = nBlk   ! the VTK module never sets nData for rectilinear datasets; CELL_DATA needs the cell count
   CALL vtk_cell_data_init( mvtk )
   CALL vtk_cell_data_scalar( Life, 'BlockLife', mvtk )
   CALL vtk_close_file( mvtk )

CONTAINS

   !> Physical coordinates of the nBlkAxis+1 block-boundary planes along one axis.
   !> Boundary i (0-based) sits at grid point i*BlkCells+1 (1-based), clamped to the last point.
   !> A degenerate axis (a single grid point) gets an arbitrary 1 m slab so cells stay renderable.
   SUBROUTINE BlockEdgeCoords ( Grid, nBlkAxis, BlkCells, Edges )
      REAL(SiKi),     INTENT(IN   ) :: Grid(:)
      INTEGER(IntKi), INTENT(IN   ) :: nBlkAxis
      INTEGER(IntKi), INTENT(IN   ) :: BlkCells
      REAL(ReKi),     INTENT(  OUT) :: Edges(:)
      INTEGER(IntKi)                :: i, nPt
      nPt = SIZE(Grid)
      IF ( nPt <= 1_IntKi ) THEN
         Edges(1) = REAL(Grid(1),ReKi) - 0.5_ReKi
         Edges(2) = REAL(Grid(1),ReKi) + 0.5_ReKi
         RETURN
      END IF
      DO i = 0, nBlkAxis
         Edges(i+1) = REAL( Grid( MIN( i*BlkCells + 1_IntKi, nPt ) ), ReKi )
      END DO
   END SUBROUTINE BlockEdgeCoords

END SUBROUTINE WaveField_WriteBlockVTK


!> Interpolate the wave-kinematics volume quantities for the point/time previously set up through
!! WaveField_Interp_Setup4D. This is the only place the volume data is read: in full-domain mode
!! (WvKinBlockMod=False) it reads the eager WaveField arrays; in on-demand mode (WvKinBlockMod=True) it
!! locates (and if needed populates) the wave block containing the interpolation stencil and reads
!! the block-local arrays — same values by construction, since blocks are filled by the same kernels
!! as the full-domain arrays. FAMCF is only written when a MacCamy-Fuchs field exists.
SUBROUTINE WaveField_InterpVol( WaveField, WaveField_m, Time, ErrStat, ErrMsg, FDynP, FV, FA, FAMCF )
   type(SeaSt_WaveFieldType),    intent(in   ) :: WaveField
   type(GridInterp_MiscVarType), intent(inout) :: WaveField_m
   real(DbKi),                   intent(in   ) :: Time      !< Simulation time, for block bookkeeping
   integer(IntKi),               intent(  out) :: ErrStat
   character(*),                 intent(  out) :: ErrMsg
   real(SiKi), optional,         intent(  out) :: FDynP
   real(SiKi), optional,         intent(  out) :: FV(3)
   real(SiKi), optional,         intent(  out) :: FA(3)
   real(SiKi), optional,         intent(  out) :: FAMCF(3)

   type(GridInterp_MiscVarType)                :: m_blk
   integer(IntKi)                              :: iBlk
   character(*), parameter                     :: RoutineName = 'WaveField_InterpVol'

   ErrStat = ErrID_None
   ErrMsg  = ""

   IF ( WaveField%WvKinBlockMod ) THEN

      CALL WaveField_EnsureBlock( WaveField, WaveField_m, Time, iBlk, ErrStat, ErrMsg )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ASSOCIATE ( Blk => WaveField%BlockStore%Blocks(iBlk) )
         ! Shift the stencil indices from global to block-local (both are 0-based inside GridInterp)
         m_blk = WaveField_m
         m_blk%Indx(:,2) = WaveField_m%Indx(:,2) - ( Blk%iPtX0 - 1_IntKi )
         m_blk%Indx(:,3) = WaveField_m%Indx(:,3) - ( Blk%iPtY0 - 1_IntKi )
         m_blk%Indx(:,4) = WaveField_m%Indx(:,4) - ( Blk%iPtZ0 - 1_IntKi )
         IF ( PRESENT(FV)    ) FV(:) = GridInterp4DVec( Blk%WaveVel,  m_blk )
         IF ( PRESENT(FA)    ) FA(:) = GridInterp4DVec( Blk%WaveAcc,  m_blk )
         IF ( PRESENT(FDynP) ) FDynP = GridInterp4D   ( Blk%WaveDynP, m_blk )
         IF ( PRESENT(FAMCF) ) THEN
            IF ( ALLOCATED(Blk%WaveAccMCF) ) THEN
               FAMCF(:) = GridInterp4DVec( Blk%WaveAccMCF, m_blk )
            END IF
         END IF
      END ASSOCIATE

   ELSE

      IF ( PRESENT(FV)    ) FV(:) = GridInterp4DVec( WaveField%WaveVel,  WaveField_m )
      IF ( PRESENT(FA)    ) FA(:) = GridInterp4DVec( WaveField%WaveAcc,  WaveField_m )
      IF ( PRESENT(FDynP) ) FDynP = GridInterp4D   ( WaveField%WaveDynP, WaveField_m )
      IF ( PRESENT(FAMCF) ) THEN
         IF ( ALLOCATED(WaveField%WaveAccMCF) ) THEN
            FAMCF(:) = GridInterp4DVec( WaveField%WaveAccMCF, WaveField_m )
         END IF
      END IF

   END IF

END SUBROUTINE WaveField_InterpVol


!----------------------------------------------------------------------------------------------------
! Interpolation related functions
!----------------------------------------------------------------------------------------------------

!====================================================================================================
!> This routine sets up interpolation of a 3-d or 4-d dataset.
!! This method is described here: http://rjwagner49.com/Mathematics/Interpolation.pdf
subroutine WaveField_Interp_Setup3D( Time, Position, p, m, ErrStat, ErrMsg )
   real(DbKi),                          intent(in   )  :: Time              !< Time from the start of the simulation
   real(ReKi),                          intent(in   )  :: Position(2)       !< Array of XY coordinates, 2
   type(GridInterp_ParameterType),      intent(in   )  :: p                 !< Parameters
   type(GridInterp_MiscVarType),        intent(inout)  :: m                 !< MiscVars
   integer(IntKi),                      intent(  out)  :: ErrStat           !< Error status
   character(*),                        intent(  out)  :: ErrMsg            !< Error message if ErrStat /= ErrID_None

   character(*), parameter              :: RoutineName = 'WaveField_Interp_Setup3D'
   integer(IntKi)                       :: ErrStat2
   character(ErrMsgLen)                 :: ErrMsg2

   ErrStat = ErrID_None

   CALL GridInterpSetup3D((/Real(Time,ReKi),Position(1),Position(2)/), p, m, ErrStat2, ErrMsg2 )
     if (Failed()) return;

contains
   logical function Failed()
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
   end function
END Subroutine WaveField_Interp_Setup3D

subroutine WaveField_Interp_Setup4D( Time, Position, GridDepth, p, m, ErrStat, ErrMsg )
   real(DbKi),                          intent(in   )  :: Time              !< Time from the start of the simulation
   real(ReKi),                          intent(in   )  :: Position(3)       !< Array of XYZ coordinates, 3
   real(SiKi),                          intent(in   )  :: GridDepth         !< Depth (>0) of the wave grid below SWL
   type(GridInterp_ParameterType),      intent(in   )  :: p                 !< Parameters
   type(GridInterp_MiscVarType),        intent(inout)  :: m                 !< MiscVars
   integer(IntKi),                      intent(  out)  :: ErrStat           !< Error status
   character(*),                        intent(  out)  :: ErrMsg            !< Error message if ErrStat /= ErrID_None

   real(ReKi)                           :: kz

   character(*), parameter              :: RoutineName = 'WaveField_Interp_Setup4D'
   integer(IntKi)                       :: ErrStat2
   character(ErrMsgLen)                 :: ErrMsg2

   ErrStat = ErrID_None
   ErrMsg  = ""

   ! Map physical z-coordinate to grid index space
   kz = 0.5_ReKi*Pi - acos( max( -1.0_ReKi, min( 1.0_ReKi, 1.0_ReKi + (Position(3) / GridDepth) ) ) )
   call GridInterpSetup4D( (/Real(Time,ReKi),Position(1),Position(2),kz/), p, m, ErrStat, ErrMsg )

contains
   logical function Failed()
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
   end function
END Subroutine WaveField_Interp_Setup4D


END MODULE SeaSt_WaveField
