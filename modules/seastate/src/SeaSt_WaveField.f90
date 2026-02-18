MODULE SeaSt_WaveField

USE GridInterp
USE SeaSt_WaveField_Types
USE IfW_FlowField, only: IfW_FlowField_GetVelAcc
USE GridInterp_Types

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

   IF (WaveField%WaveStMod == 0) THEN ! No wave stretching

      IF ( pos(3) <= 0.0_ReKi) THEN ! Node is at or below the SWL
         nodeInWater = 1_IntKi
         ! Use location to obtain interpolated values of kinematics
         CALL WaveField_Interp_Setup4D( Time+WaveField%WaveTimeShift, pos, WaveField%GridDepth, WaveField%VolGridParams, WaveField_m, ErrStat2, ErrMsg2 ); if (Failed()) return;
         FV(:) = GridInterp4DVec( WaveField%WaveVel,  WaveField_m )
         FA(:) = GridInterp4DVec( WaveField%WaveAcc,  WaveField_m )
         FDynP = GridInterp4D   ( WaveField%WaveDynP, WaveField_m )
         IF ( ALLOCATED(WaveField%WaveAccMCF) ) THEN
            FAMCF(:) = GridInterp4DVec( WaveField%WaveAccMCF, WaveField_m )
         END IF
      ELSE ! Node is above the SWL
         nodeInWater = 0_IntKi
         FV(:)       = 0.0
         FA(:)       = 0.0
         FDynP       = 0.0
         FAMCF(:)    = 0.0
      END IF

   ELSE ! Wave stretching enabled

      IF ( (pos(3) <= WaveElev) .OR. forceNodeInWater ) THEN ! Node is submerged

         nodeInWater = 1_IntKi

         IF ( WaveField%WaveStMod < 3 ) THEN ! Vertical or extrapolated wave stretching

            IF ( pos(3) <= 0.0_SiKi) THEN ! Node is below the SWL - evaluate wave dynamics as usual

               ! Use location to obtain interpolated values of kinematics
               CALL WaveField_Interp_Setup4D( Time+WaveField%WaveTimeShift, pos, WaveField%GridDepth, WaveField%VolGridParams, WaveField_m, ErrStat2, ErrMsg2 ); if (Failed()) return;
               FV(:) = GridInterp4DVec( WaveField%WaveVel,  WaveField_m )
               FA(:) = GridInterp4DVec( WaveField%WaveAcc,  WaveField_m )
               FDynP = GridInterp4D   ( WaveField%WaveDynP, WaveField_m )
               IF ( ALLOCATED(WaveField%WaveAccMCF) ) THEN
                  FAMCF(:) = GridInterp4DVec( WaveField%WaveAccMCF, WaveField_m )
               END IF

            ELSE ! Node is above SWL - need wave stretching

               ! Vertical wave stretching
               CALL WaveField_Interp_Setup4D( Time+WaveField%WaveTimeShift, posXY0, WaveField%GridDepth, WaveField%VolGridParams, WaveField_m, ErrStat2, ErrMsg2 ); if (Failed()) return;
               FV(:) = GridInterp4DVec( WaveField%WaveVel,  WaveField_m )
               FA(:) = GridInterp4DVec( WaveField%WaveAcc,  WaveField_m )
               FDynP = GridInterp4D   ( WaveField%WaveDynP, WaveField_m )
               IF ( ALLOCATED(WaveField%WaveAccMCF) ) THEN
                  FAMCF(:) = GridInterp4DVec( WaveField%WaveAccMCF, WaveField_m )
               END IF

               ! Extrapolated wave stretching
               IF (WaveField%WaveStMod == 2) THEN
                  CALL WaveField_Interp_Setup3D( Time+WaveField%WaveTimeShift, posXY, WaveField%SrfGridParams, WaveField_m, ErrStat2, ErrMsg2 ); if (Failed()) return;
                  FV(:) = FV(:) + GridInterp3DVec( WaveField%PWaveVel0,  WaveField_m ) * pos(3)
                  FA(:) = FA(:) + GridInterp3DVec( WaveField%PWaveAcc0,  WaveField_m ) * pos(3)
                  FDynP = FDynP + GridInterp3D   ( WaveField%PWaveDynP0, WaveField_m ) * pos(3)
                  IF ( ALLOCATED(WaveField%WaveAccMCF) ) THEN
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
            FV(:) = GridInterp4DVec( WaveField%WaveVel,  WaveField_m )
            FA(:) = GridInterp4DVec( WaveField%WaveAcc,  WaveField_m )
            FDynP = GridInterp4D   ( WaveField%WaveDynP, WaveField_m )
            IF ( ALLOCATED(WaveField%WaveAccMCF) ) THEN
               FAMCF(:) = GridInterp4DVec( WaveField%WaveAccMCF, WaveField_m )
            END IF
         END IF

      ELSE ! Node is out of water - zero-out all wave dynamics

         nodeInWater = 0_IntKi
         FV(:)       = 0.0
         FA(:)       = 0.0
         FDynP       = 0.0
         FAMCF(:)    = 0.0

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

   IF (WaveField%WaveStMod == 0) THEN ! No wave stretching

      IF ( pos(3) <= 0.0_ReKi) THEN ! Node is at or below the SWL
         nodeInWater = 1_IntKi
         ! Use location to obtain interpolated values of kinematics
         CALL WaveField_Interp_Setup4D( Time+WaveField%WaveTimeShift, pos, WaveField%GridDepth, WaveField%VolGridParams, WaveField_m, ErrStat2, ErrMsg2 ); if (Failed()) return;
         FV(:) = GridInterp4DVec( WaveField%WaveVel,  WaveField_m )
      ELSE ! Node is above the SWL
         nodeInWater = 0_IntKi
         FV(:)       = 0.0
      END IF

   ELSE ! Wave stretching enabled

      IF ( (pos(3) <= WaveElev) .OR. forceNodeInWater ) THEN ! Node is submerged

         nodeInWater = 1_IntKi

         IF ( WaveField%WaveStMod < 3 ) THEN ! Vertical or extrapolated wave stretching

            IF ( pos(3) <= 0.0_SiKi) THEN ! Node is below the SWL - evaluate wave dynamics as usual

               ! Use location to obtain interpolated values of kinematics
               CALL WaveField_Interp_Setup4D( Time+WaveField%WaveTimeShift, pos, WaveField%GridDepth, WaveField%VolGridParams, WaveField_m, ErrStat2, ErrMsg2 ); if (Failed()) return;
               FV(:) = GridInterp4DVec( WaveField%WaveVel,  WaveField_m )

            ELSE ! Node is above SWL - need wave stretching

               ! Vertical wave stretching
               CALL WaveField_Interp_Setup4D( Time, posXY0, WaveField%GridDepth, WaveField%VolGridParams, WaveField_m, ErrStat2, ErrMsg2 ); if (Failed()) return;
               FV(:) = GridInterp4DVec( WaveField%WaveVel,  WaveField_m )

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
            FV(:) = GridInterp4DVec( WaveField%WaveVel,  WaveField_m )

         END IF

      ELSE ! Node is out of water - zero-out all wave dynamics

         nodeInWater = 0_IntKi
         FV(:)       = 0.0

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
    
   IF (WaveField%WaveStMod == 0) THEN ! No wave stretching

      IF ( pos(3) <= 0.0_ReKi) THEN ! Node is at or below the SWL
         nodeInWater = 1_IntKi
         ! Use location to obtain interpolated values of kinematics
         CALL WaveField_Interp_Setup4D( Time+WaveField%WaveTimeShift, pos, WaveField%GridDepth, WaveField%VolGridParams, WaveField_m, ErrStat2, ErrMsg2 ); if (Failed()) return;
         FV(:) = GridInterp4DVec( WaveField%WaveVel,  WaveField_m )
         FA(:) = GridInterp4DVec( WaveField%WaveAcc,  WaveField_m )
      ELSE ! Node is above the SWL
         nodeInWater = 0_IntKi
         FV(:)       = 0.0
         FA(:)       = 0.0
      END IF

   ELSE ! Wave stretching enabled

      IF ( (pos(3) <= WaveElev) .OR. forceNodeInWater ) THEN ! Node is submerged

         nodeInWater = 1_IntKi

         IF ( WaveField%WaveStMod < 3 ) THEN ! Vertical or extrapolated wave stretching

            IF ( pos(3) <= 0.0_SiKi) THEN ! Node is below the SWL - evaluate wave dynamics as usual

               ! Use location to obtain interpolated values of kinematics
               CALL WaveField_Interp_Setup4D( Time+WaveField%WaveTimeShift, pos, WaveField%GridDepth, WaveField%VolGridParams, WaveField_m, ErrStat2, ErrMsg2 ); if (Failed()) return;
               FV(:) = GridInterp4DVec( WaveField%WaveVel,  WaveField_m )
               FA(:) = GridInterp4DVec( WaveField%WaveAcc,  WaveField_m )

            ELSE ! Node is above SWL - need wave stretching

               ! Vertical wave stretching
               CALL WaveField_Interp_Setup4D( Time, posXY0, WaveField%GridDepth, WaveField%VolGridParams, WaveField_m, ErrStat2, ErrMsg2 ); if (Failed()) return;
               FV(:) = GridInterp4DVec( WaveField%WaveVel,  WaveField_m )
               FA(:) = GridInterp4DVec( WaveField%WaveAcc,  WaveField_m )

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
            FV(:) = GridInterp4DVec( WaveField%WaveVel,  WaveField_m )
            FA(:) = GridInterp4DVec( WaveField%WaveAcc,  WaveField_m )
         END IF

      ELSE ! Node is out of water - zero-out all wave dynamics

         nodeInWater = 0_IntKi
         FV(:)       = 0.0
         FA(:)       = 0.0
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
      IF (ALLOCATED(WaveField%WaveAccMCF)) THEN
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
SUBROUTINE WaveField_GetWaveVelAcc_AD( WaveField, WaveField_m, StartNode, Time, pos, FV, FA, ErrStat, ErrMsg )
   type(SeaSt_WaveFieldType),          intent(in   ) :: WaveField
   type(GridInterp_MiscVarType),       intent(inout) :: WaveField_m
   integer(IntKi),                     intent(in   ) :: StartNode
   real(DbKi),                         intent(in   ) :: Time
   real(ReKi),                         intent(in   ) :: pos(:,:) ! z=0 at MSL
   real(ReKi),                         intent(  out) :: FV(:,:)
   real(ReKi),       allocatable,      intent(inout) :: FA(:,:)
   integer(IntKi),                     intent(  out) :: ErrStat  ! Error status of the operation
   character(*),                       intent(  out) :: ErrMsg   ! Error message if errStat /= ErrID_None
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
      CALL IfW_FlowField_GetVelAcc(WaveField%CurrField, StartNode, Time, pos, FV_DC, FA_DC, ErrStat2, ErrMsg2, PosOffset=PosOffset); if (Failed()) return;

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
