MODULE SeaSt_WaveField

USE SeaSt_WaveField_Types

IMPLICIT NONE

PRIVATE

! Public functions and subroutines
PUBLIC WaveField_GetNodeWaveElev1
PUBLIC WaveField_GetNodeWaveElev2
PUBLIC WaveField_GetNodeTotalWaveElev
PUBLIC WaveField_GetNodeWaveNormal
PUBLIC WaveField_GetNodeWaveKin
PUBLIC WaveField_GetNodeWaveVel

PUBLIC WaveField_GetWaveKin

public WaveField_Interp_Setup3D, WaveField_Interp_Setup4D

CONTAINS

!-------------------- Subroutine for wave elevation ------------------!
function WaveField_GetNodeWaveElev1( WaveField, WaveField_m, Time, pos, ErrStat, ErrMsg )
   type(SeaSt_WaveFieldType),          intent(in   ) :: WaveField
   type(SeaSt_WaveField_MiscVarType),  intent(inout) :: WaveField_m
   real(DbKi),                         intent(in   ) :: Time
   real(ReKi),                         intent(in   ) :: pos(*)  ! Position at which free-surface elevation is to be calculated. Third entry ignored if present.
   integer(IntKi),                     intent(  out) :: ErrStat ! Error status of the operation
   character(*),                       intent(  out) :: ErrMsg  ! Error message if errStat /= ErrID_None

   real(SiKi)                                        :: WaveField_GetNodeWaveElev1
   real(SiKi)                                        :: Zeta
   character(*),                       parameter     :: RoutineName = 'WaveField_GetNodeWaveElev1'
   integer(IntKi)                                    :: errStat2
   character(ErrMsgLen)                              :: errMsg2

   ErrStat   = ErrID_None
   ErrMsg    = ""

   IF (ALLOCATED(WaveField%WaveElev1)) THEN
      CALL WaveField_Interp_Setup3D( Time, pos, WaveField%GridParams, WaveField_m, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Zeta = WaveField_Interp_3D( WaveField%WaveElev1, WaveField_m )
   ELSE
      Zeta = 0.0_SiKi
   END IF

   WaveField_GetNodeWaveElev1 = Zeta

end function WaveField_GetNodeWaveElev1


function WaveField_GetNodeWaveElev2( WaveField, WaveField_m, Time, pos, ErrStat, ErrMsg )
   type(SeaSt_WaveFieldType),          intent(in   ) :: WaveField
   type(SeaSt_WaveField_MiscVarType),  intent(inout) :: WaveField_m
   real(DbKi),                         intent(in   ) :: Time
   real(ReKi),                         intent(in   ) :: pos(*)  ! Position at which free-surface elevation is to be calculated. Third entry ignored if present.
   integer(IntKi),                     intent(  out) :: ErrStat ! Error status of the operation
   character(*),                       intent(  out) :: ErrMsg  ! Error message if errStat /= ErrID_None

   real(SiKi)                                        :: WaveField_GetNodeWaveElev2
   real(SiKi)                                        :: Zeta
   character(*),                       parameter     :: RoutineName = 'WaveField_GetNodeWaveElev2'
   integer(IntKi)                                    :: errStat2
   character(ErrMsgLen)                              :: errMsg2

   ErrStat   = ErrID_None
   ErrMsg    = ""

   IF (ALLOCATED(WaveField%WaveElev2)) THEN
      CALL WaveField_Interp_Setup3D( Time, pos, WaveField%GridParams, WaveField_m, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Zeta = WaveField_Interp_3D( WaveField%WaveElev2, WaveField_m )
   ELSE
      Zeta = 0.0_SiKi
   END IF

   WaveField_GetNodeWaveElev2 = Zeta

end function WaveField_GetNodeWaveElev2


FUNCTION WaveField_GetNodeTotalWaveElev( WaveField, WaveField_m, Time, pos, ErrStat, ErrMsg )
   type(SeaSt_WaveFieldType),          intent(in   ) :: WaveField
   type(SeaSt_WaveField_MiscVarType),  intent(inout) :: WaveField_m
   real(DbKi),                         intent(in   ) :: Time
   real(ReKi),                         intent(in   ) :: pos(*)  ! Position at which free-surface elevation is to be calculated. Third entry ignored if present.
   integer(IntKi),                     intent(  out) :: ErrStat ! Error status of the operation
   character(*),                       intent(  out) :: ErrMsg  ! Error message if errStat /= ErrID_None

   real(SiKi)                                        :: WaveField_GetNodeTotalWaveElev
   real(SiKi)                                        :: Zeta1, Zeta2
   character(*),                       parameter     :: RoutineName = 'WaveField_GetNodeTotalWaveElev'
   integer(IntKi)                                    :: errStat2
   character(ErrMsgLen)                              :: errMsg2

   ErrStat   = ErrID_None
   ErrMsg    = ""

   Zeta1 = WaveField_GetNodeWaveElev1( WaveField, WaveField_m, Time, pos, ErrStat2, ErrMsg2 ); if (Failed()) return;
   Zeta2 = WaveField_GetNodeWaveElev2( WaveField, WaveField_m, Time, pos, ErrStat2, ErrMsg2 ); if (Failed()) return;
   WaveField_GetNodeTotalWaveElev = Zeta1 + Zeta2

contains
   logical function Failed()
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
   end function
END FUNCTION WaveField_GetNodeTotalWaveElev


SUBROUTINE WaveField_GetNodeWaveNormal( WaveField, WaveField_m, Time, pos, r, n, ErrStat, ErrMsg )
   type(SeaSt_WaveFieldType),          intent(in   ) :: WaveField
   type(SeaSt_WaveField_MiscVarType),  intent(inout) :: WaveField_m
   real(DbKi),                         intent(in   ) :: Time
   real(ReKi),                         intent(in   ) :: pos(*)  ! Position at which free-surface normal is to be calculated. Third entry ignored if present.
   real(ReKi),                         intent(in   ) :: r       ! Distance for central differencing
   real(ReKi),                         intent(  out) :: n(3)    ! Free-surface normal vector
   integer(IntKi),                     intent(  out) :: ErrStat ! Error status of the operation
   character(*),                       intent(  out) :: ErrMsg  ! Error message if errStat /= ErrID_None

   real(SiKi)                                        :: ZetaP,ZetaM
   real(ReKi)                                        :: r1,dZetadx,dZetady
   character(*),                       parameter     :: RoutineName = 'WaveField_GetNodeWaveNormal'
   integer(IntKi)                                    :: errStat2
   character(ErrMsgLen)                              :: errMsg2

   ErrStat   = ErrID_None
   ErrMsg    = ""

   r1 = MAX(r,real(1.0e-6,ReKi)) ! In case r is zero

   ZetaP = WaveField_GetNodeTotalWaveElev( WaveField, WaveField_m, Time, (/pos(1)+r1,pos(2)/), ErrStat2, ErrMsg2 ); if (Failed()) return;
   ZetaM = WaveField_GetNodeTotalWaveElev( WaveField, WaveField_m, Time, (/pos(1)-r1,pos(2)/), ErrStat2, ErrMsg2 ); if (Failed()) return;
   dZetadx = REAL(ZetaP-ZetaM,ReKi)/(2.0_ReKi*r1)

   ZetaP = WaveField_GetNodeTotalWaveElev( WaveField, WaveField_m, Time, (/pos(1),pos(2)+r1/), ErrStat2, ErrMsg2 ); if (Failed()) return;
   ZetaM = WaveField_GetNodeTotalWaveElev( WaveField, WaveField_m, Time, (/pos(1),pos(2)-r1/), ErrStat2, ErrMsg2 ); if (Failed()) return;
   dZetady = REAL(ZetaP-ZetaM,ReKi)/(2.0_ReKi*r1)

   n = (/-dZetadx,-dZetady,1.0_ReKi/)
   n = n / SQRT(Dot_Product(n,n))

contains
   logical function Failed()
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
   end function
END SUBROUTINE WaveField_GetNodeWaveNormal


!-------------------- Subroutine for full wave field kinematics --------------------!
SUBROUTINE WaveField_GetNodeWaveKin( WaveField, WaveField_m, Time, pos, forceNodeInWater, nodeInWater, WaveElev1, WaveElev2, WaveElev, FDynP, FV, FA, FAMCF, ErrStat, ErrMsg )
   type(SeaSt_WaveFieldType),          intent(in   ) :: WaveField
   type(SeaSt_WaveField_MiscVarType),  intent(inout) :: WaveField_m
   real(DbKi),                         intent(in   ) :: Time
   real(ReKi),                         intent(in   ) :: pos(3)
   logical,                            intent(in   ) :: forceNodeInWater
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

   real(ReKi)                                        :: posXY(2), posPrime(3), posXY0(3)
   character(*),                       parameter     :: RoutineName = 'WaveField_GetNodeWaveKin'
   integer(IntKi)                                    :: errStat2
   character(ErrMsgLen)                              :: errMsg2

   ErrStat   = ErrID_None
   ErrMsg    = ""

   posXY    = pos(1:2)
   posXY0   = (/pos(1),pos(2),0.0_ReKi/)
   FAMCF(:) = 0.0

   ! Wave elevation
   WaveElev1 = WaveField_GetNodeWaveElev1( WaveField, WaveField_m, Time, pos, ErrStat2, ErrMsg2 ); if (Failed()) return;
   WaveElev2 = WaveField_GetNodeWaveElev2( WaveField, WaveField_m, Time, pos, ErrStat2, ErrMsg2 ); if (Failed()) return;
   WaveElev  = WaveElev1 + WaveElev2

   IF (WaveField%WaveStMod == 0) THEN ! No wave stretching

      IF ( pos(3) <= 0.0_ReKi) THEN ! Node is at or below the SWL
         nodeInWater = 1_IntKi
         ! Use location to obtain interpolated values of kinematics
         CALL WaveField_Interp_Setup4D( Time, pos, WaveField%GridParams, WaveField_m, ErrStat2, ErrMsg2 ); if (Failed()) return;
         FV(:) = WaveField_Interp_4D_Vec( WaveField%WaveVel,  WaveField_m )
         FA(:) = WaveField_Interp_4D_Vec( WaveField%WaveAcc,  WaveField_m )
         FDynP = WaveField_Interp_4D    ( WaveField%WaveDynP, WaveField_m )
         IF ( ALLOCATED(WaveField%WaveAccMCF) ) THEN
            FAMCF(:) = WaveField_Interp_4D_Vec( WaveField%WaveAccMCF, WaveField_m )
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
               CALL WaveField_Interp_Setup4D( Time, pos, WaveField%GridParams, WaveField_m, ErrStat2, ErrMsg2 ); if (Failed()) return;
               FV(:) = WaveField_Interp_4D_Vec( WaveField%WaveVel,  WaveField_m )
               FA(:) = WaveField_Interp_4D_Vec( WaveField%WaveAcc,  WaveField_m )
               FDynP = WaveField_Interp_4D    ( WaveField%WaveDynP, WaveField_m )
               IF ( ALLOCATED(WaveField%WaveAccMCF) ) THEN
                  FAMCF(:) = WaveField_Interp_4D_Vec( WaveField%WaveAccMCF, WaveField_m )
               END IF

            ELSE ! Node is above SWL - need wave stretching

               ! Vertical wave stretching
               CALL WaveField_Interp_Setup4D( Time, posXY0, WaveField%GridParams, WaveField_m, ErrStat2, ErrMsg2 ); if (Failed()) return;
               FV(:) = WaveField_Interp_4D_vec( WaveField%WaveVel,  WaveField_m )
               FA(:) = WaveField_Interp_4D_vec( WaveField%WaveAcc,  WaveField_m )
               FDynP = WaveField_Interp_4D    ( WaveField%WaveDynP, WaveField_m )
               IF ( ALLOCATED(WaveField%WaveAccMCF) ) THEN
                  FAMCF(:) = WaveField_Interp_4D_vec( WaveField%WaveAccMCF, WaveField_m )
               END IF

               ! Extrapoled wave stretching
               IF (WaveField%WaveStMod == 2) THEN
                  CALL WaveField_Interp_Setup3D( Time, posXY, WaveField%GridParams, WaveField_m, ErrStat2, ErrMsg2 ); if (Failed()) return;
                  FV(:) = FV(:) + WaveField_Interp_3D_vec( WaveField%PWaveVel0,  WaveField_m ) * pos(3)
                  FA(:) = FA(:) + WaveField_Interp_3D_vec( WaveField%PWaveAcc0,  WaveField_m ) * pos(3)
                  FDynP = FDynP + WaveField_Interp_3D    ( WaveField%PWaveDynP0, WaveField_m ) * pos(3)
                  IF ( ALLOCATED(WaveField%WaveAccMCF) ) THEN
                     FAMCF(:) = FAMCF(:) + WaveField_Interp_3D_vec( WaveField%PWaveAccMCF0, WaveField_m ) * pos(3)
                  END IF
               END IF

            END IF ! Node is submerged

         ELSE ! Wheeler stretching - no need to check whether the node is above or below SWL

            ! Map the node z-position linearly from [-EffWtrDpth,m%WaveElev(j)] to [-EffWtrDpth,0]
            posPrime    = pos
            posPrime(3) = WaveField%EffWtrDpth*(WaveField%EffWtrDpth+pos(3))/(WaveField%EffWtrDpth+WaveElev)-WaveField%EffWtrDpth
            posPrime(3) = MIN( posPrime(3), 0.0_ReKi) ! Clamp z-position to zero. Needed when forceNodeInWater=.TRUE.

            ! Obtain the wave-field variables by interpolation with the mapped position.
            CALL WaveField_Interp_Setup4D( Time, posPrime, WaveField%GridParams, WaveField_m, ErrStat2, ErrMsg2 ); if (Failed()) return;
            FV(:) = WaveField_Interp_4D_Vec( WaveField%WaveVel,  WaveField_m )
            FA(:) = WaveField_Interp_4D_Vec( WaveField%WaveAcc,  WaveField_m )
            FDynP = WaveField_Interp_4D    ( WaveField%WaveDynP, WaveField_m )
            IF ( ALLOCATED(WaveField%WaveAccMCF) ) THEN
               FAMCF(:) = WaveField_Interp_4D_Vec( WaveField%WaveAccMCF, WaveField_m )
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

contains
   logical function Failed()
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
   end function
END SUBROUTINE WaveField_GetNodeWaveKin


!-------------------- Subroutine for wave field velocity only --------------------!
SUBROUTINE WaveField_GetNodeWaveVel( WaveField, WaveField_m, Time, pos, forceNodeInWater, nodeInWater, FV, ErrStat, ErrMsg )
   type(SeaSt_WaveFieldType),          intent(in   ) :: WaveField
   type(SeaSt_WaveField_MiscVarType),  intent(inout) :: WaveField_m
   real(DbKi),                         intent(in   ) :: Time
   real(ReKi),                         intent(in   ) :: pos(3)
   logical,                            intent(in   ) :: forceNodeInWater
   integer(IntKi),                     intent(  out) :: nodeInWater
   real(SiKi),                         intent(  out) :: FV(3)
   integer(IntKi),                     intent(  out) :: ErrStat ! Error status of the operation
   character(*),                       intent(  out) :: ErrMsg  ! Error message if errStat /= ErrID_None

   real(SiKi)                                        :: WaveElev
   real(ReKi)                                        :: posXY(2), posPrime(3), posXY0(3)
   character(*),                       parameter     :: RoutineName = 'WaveField_GetNodeWaveVel'
   integer(IntKi)                                    :: errStat2
   character(ErrMsgLen)                              :: errMsg2

   ErrStat   = ErrID_None
   ErrMsg    = ""

   posXY    = pos(1:2)
   posXY0   = (/pos(1),pos(2),0.0_ReKi/)

   ! Wave elevation
   WaveElev  = WaveField_GetNodeTotalWaveElev( WaveField, WaveField_m, Time, pos, ErrStat2, ErrMsg2 ); if (Failed()) return;

   IF (WaveField%WaveStMod == 0) THEN ! No wave stretching

      IF ( pos(3) <= 0.0_ReKi) THEN ! Node is at or below the SWL
         nodeInWater = 1_IntKi
         ! Use location to obtain interpolated values of kinematics
         CALL WaveField_Interp_Setup4D( Time, pos, WaveField%GridParams, WaveField_m, ErrStat2, ErrMsg2 ); if (Failed()) return;
         FV(:) = WaveField_Interp_4D_Vec( WaveField%WaveVel,  WaveField_m )
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
               CALL WaveField_Interp_Setup4D( Time, pos, WaveField%GridParams, WaveField_m, ErrStat2, ErrMsg2 ); if (Failed()) return;
               FV(:) = WaveField_Interp_4D_Vec( WaveField%WaveVel,  WaveField_m )

            ELSE ! Node is above SWL - need wave stretching

               ! Vertical wave stretching
               CALL WaveField_Interp_Setup4D( Time, posXY0, WaveField%GridParams, WaveField_m, ErrStat2, ErrMsg2 ); if (Failed()) return;
               FV(:) = WaveField_Interp_4D_vec( WaveField%WaveVel,  WaveField_m )

               ! Extrapoled wave stretching
               IF (WaveField%WaveStMod == 2) THEN
                  CALL WaveField_Interp_Setup3D( Time, posXY, WaveField%GridParams, WaveField_m, ErrStat2, ErrMsg2 ); if (Failed()) return;
                  FV(:) = FV(:) + WaveField_Interp_3D_vec( WaveField%PWaveVel0, WaveField_m ) * pos(3)
               END IF

            END IF ! Node is submerged

         ELSE ! Wheeler stretching - no need to check whether the node is above or below SWL

            ! Map the node z-position linearly from [-EffWtrDpth,m%WaveElev(j)] to [-EffWtrDpth,0]
            posPrime    = pos
            posPrime(3) = WaveField%EffWtrDpth*(WaveField%EffWtrDpth+pos(3))/(WaveField%EffWtrDpth+WaveElev)-WaveField%EffWtrDpth
            posPrime(3) = MIN( posPrime(3), 0.0_ReKi) ! Clamp z-position to zero. Needed when forceNodeInWater=.TRUE.

            ! Obtain the wave-field variables by interpolation with the mapped position.
            CALL WaveField_Interp_Setup4D( Time, posPrime, WaveField%GridParams, WaveField_m, ErrStat2, ErrMsg2 ); if (Failed()) return;
            FV(:) = WaveField_Interp_4D_Vec( WaveField%WaveVel,  WaveField_m )

         END IF

      ELSE ! Node is out of water - zero-out all wave dynamics

         nodeInWater = 0_IntKi
         FV(:)       = 0.0

      END IF ! If node is in or out of water

   END IF ! If wave stretching is on or off

contains
   logical function Failed()
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
   end function
END SUBROUTINE WaveField_GetNodeWaveVel


SUBROUTINE WaveField_GetWaveKin( WaveField, WaveField_m, Time, pos, forceNodeInWater, nodeInWater, WaveElev1, WaveElev2, WaveElev, FDynP, FV, FA, FAMCF, ErrStat, ErrMsg )
   type(SeaSt_WaveFieldType),          intent(in   ) :: WaveField
   type(SeaSt_WaveField_MiscVarType),  intent(inout) :: WaveField_m
   real(DbKi),                         intent(in   ) :: Time
   real(ReKi),                         intent(in   ) :: pos(:,:)
   logical,                            intent(in   ) :: forceNodeInWater
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

   integer(IntKi)                                    :: NumPoints, i
   real(SiKi)                                        :: FDynP_node, FV_node(3), FA_node(3), FAMCF_node(3)

   ErrStat   = ErrID_None
   ErrMsg    = ""

   NumPoints = size(pos, dim=2)
   DO i = 1, NumPoints
      CALL WaveField_GetNodeWaveKin( WaveField, WaveField_m, Time, pos(:,i), forceNodeInWater, nodeInWater(i), WaveElev1(i), WaveElev2(i), WaveElev(i), FDynP_node, FV_node, FA_node, FAMCF_node, ErrStat2, ErrMsg2 )
      if (Failed()) return;
      FDynP(i) = REAL(FDynP_node,ReKi)
      FV(:, i) = REAL(FV_node,   ReKi)
      FA(:, i) = REAL(FA_node,   ReKi)
      IF (ALLOCATED(WaveField%WaveAccMCF)) THEN
         FAMCF(:,i) = REAL(FAMCF_node,ReKi)
      END IF
   END DO

contains
   logical function Failed()
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine WaveField_GetWaveKin


!----------------------------------------------------------------------------------------------------
! Interpolation related functions
!----------------------------------------------------------------------------------------------------

subroutine SetCartesianXYIndex(p, pZero, delta, nMax, Indx_Lo, Indx_Hi, isopc, FirstWarn, ErrStat, ErrMsg)
   REAL(ReKi),       intent(in   )  :: p
   REAL(ReKi),       intent(in   )  :: pZero
   REAL(ReKi),       intent(in   )  :: delta
   INTEGER(IntKi),   intent(in   )  :: nMax
   INTEGER(IntKi),   intent(inout)  :: Indx_Lo
   INTEGER(IntKi),   intent(inout)  :: Indx_Hi
   real(SiKi),       intent(inout)  :: isopc
   logical,          intent(inout)  :: FirstWarn
   INTEGER(IntKi),   intent(  out)  :: ErrStat
   CHARACTER(*),     intent(  out)  :: ErrMsg

   real(ReKi)                       :: Tmp

   ErrStat = ErrID_None
   ErrMsg  = ""

   isopc   = -1.0
   Indx_Lo = 0
   Indx_Hi = 0

   if ( nMax .EQ. 1_IntKi ) then ! Only one grid point
      Indx_Lo = 1_IntKi
      Indx_Hi = 1_IntKi
      isopc   = 0_SiKi
      return
   end if

   Tmp =  (p-pZero) / delta
   Indx_Lo = INT( Tmp ) + 1    ! convert REAL to INTEGER, then add one since our grid indices start at 1, not 0
   isopc = 2.0_ReKi * (Tmp - REAL(Indx_Lo - 1, ReKi)) - 1.0_ReKi  ! convert to value between -1 and 1

   if ( Indx_Lo < 1 ) then
      Indx_Lo = 1
      isopc = -1.0
      if (FirstWarn) then
         call SetErrStat(ErrID_Warn,'Position has been clamped to the grid boundary. Warning will not be repeated though condition may persist.',ErrStat,ErrMsg,'SetCartesianXYIndex') !error out if time is outside the lower bounds
         FirstWarn = .false.
      end if
   end if

   Indx_Hi = min( Indx_Lo + 1, nMax )     ! make sure it's a valid index, zero-based

   if ( Indx_Lo >= Indx_Hi ) then
      ! Need to clamp to grid boundary
      if (FirstWarn .and. Indx_Lo /= Indx_Hi) then ! don't warn if we are exactly at the boundary
         call SetErrStat(ErrID_Warn,'Position has been clamped to the grid boundary. Warning will not be repeated though condition may persist.',ErrStat,ErrMsg,'SetCartesianXYIndex') !error out if time is outside the lower bounds
         FirstWarn = .false.
      end if
      Indx_Lo = max(Indx_Hi - 1, 1)
      isopc = 1.0
   end if

   !-------------------------------------------------------------------------------------------------
   ! to verify that we don't extrapolate, make sure isopc is bound between -1 and 1 (effectively nearest neighbor)
   !-------------------------------------------------------------------------------------------------
   isopc = min( 1.0_SiKi, isopc )
   isopc = max(-1.0_SiKi, isopc )

end subroutine SetCartesianXYIndex


subroutine SetCartesianZIndex(p, z_depth, delta, nMax, Indx_Lo, Indx_Hi, isopc, FirstWarn, ErrStat, ErrMsg)
   real(ReKi),        intent(in   )  :: p
   real(ReKi),        intent(in   )  :: z_depth
   real(ReKi),        intent(in   )  :: delta
   integer(IntKi),    intent(in   )  :: nMax
   integer(IntKi),    intent(inout)  :: Indx_Lo
   integer(IntKi),    intent(inout)  :: Indx_Hi
   real(SiKi),        intent(inout)  :: isopc
   logical,           intent(inout)  :: FirstWarn
   integer(IntKi),    intent(  out)  :: ErrStat
   character(*),      intent(  out)  :: ErrMsg

   real(ReKi)                        :: Tmp

   ErrStat = ErrID_None
   ErrMsg  = ""

   isopc   = -1.0
   Indx_Lo = 0
   Indx_Hi = 0


   !Tmp =  acos(-p / z_depth) / delta
   Tmp = acos( max(-1.0_ReKi, min(1.0_ReKi, 1+(p / z_depth)) ) ) / delta
   Tmp =  nmax - 1 - Tmp
   Indx_Lo = INT( Tmp ) + 1    ! convert REAL to INTEGER, then add one since our grid indices start at 1, not 0
   isopc = 2.0_ReKi * (Tmp - REAL(Indx_Lo - 1, ReKi)) - 1.0_ReKi  ! convert to value between -1 and 1

   if ( Indx_Lo < 1 ) then
      Indx_Lo = 1
      isopc = -1.0
      if (FirstWarn) then
         call SetErrStat(ErrID_Warn,'Position has been clamped to the grid boundary. Warning will not be repeated though condition may persist.',ErrStat,ErrMsg,'SetCartesianZIndex') !error out if z is outside the lower bounds
         FirstWarn = .false.
      end if
   end if

   Indx_Hi = min( Indx_Lo + 1, nMax )     ! make sure it's a valid index, one-based

   if ( Indx_Lo >= Indx_Hi ) then
      ! Need to clamp to grid boundary
      if (FirstWarn .and. Indx_Lo /= Indx_Hi) then ! don't warn if we are exactly at the boundary
         call SetErrStat(ErrID_Warn,'Position has been clamped to the grid boundary. Warning will not be repeated though condition may persist.',ErrStat,ErrMsg,'SetCartesianZIndex') !error out if z is outside the upper bounds
         FirstWarn = .false.
      end if
      Indx_Lo = max(Indx_Hi - 1, 1)
      isopc = 1.0
   end if

   !-------------------------------------------------------------------------------------------------
   ! to verify that we don't extrapolate, make sure isopc is bound between -1 and 1 (effectively nearest neighbor)
   !-------------------------------------------------------------------------------------------------
   isopc = min( 1.0_SiKi, isopc )
   isopc = max(-1.0_SiKi, isopc )

end subroutine SetCartesianZIndex


subroutine SetTimeIndex(Time, deltaT, nMax, Indx_Lo, Indx_Hi, isopc, ErrStat, ErrMsg)
   real(DbKi),        intent(in   )  :: Time     !< time from the start of the simulation
   real(ReKi),        intent(in   )  :: deltaT
   integer(IntKi),    intent(in   )  :: nMax
   integer(IntKi),    intent(inout)  :: Indx_Lo
   integer(IntKi),    intent(inout)  :: Indx_Hi
   real(SiKi),        intent(inout)  :: isopc
   integer(IntKi),    intent(  out)  :: ErrStat
   character(*),      intent(  out)  :: ErrMsg

   real(ReKi)                        :: Tmp

   ErrStat = ErrID_None
   ErrMsg  = ""

   isopc   = -1.0
   Indx_Lo = 0
   Indx_Hi = 0
   if ( Time < 0.0_DbKi ) then
      CALL SetErrStat(ErrID_Fatal,'Time value must be greater than or equal to zero!',ErrStat,ErrMsg,'SetTimeIndex') !error out if time is outside the lower bounds
      RETURN
   end if

   ! if there are no timesteps, don't proceed
   if (EqualRealNos(deltaT,0.0_ReKi) .or. deltaT < 0.0_ReKi)  return;

! NOTE: nMax is the total number of time values in the grid, since this is zero-based indexing, the max index is nMax-1
!       for example: in a time grid with 11 grid points, the indices run from 0,1,2,3,4,5,6,7,8,9,10
!                    for the repeating waves feature, index 10 is the same as index 0, so if Indx_Lo = 10 then we want to
!                    wrap it back to index 0, if Indx_Lo = 11 we want to wrap back to index 1.

   Tmp =  real( (Time/ real(deltaT,DbKi)) ,ReKi)
   Tmp =  MOD(Tmp,real((nMax), ReKi))
   Indx_Lo = INT( Tmp )     ! convert REAL to INTEGER

   isopc = 2.0_ReKi * (Tmp - REAL(Indx_Lo , ReKi)) - 1.0_ReKi  ! convert to value between -1 and 1

   !-------------------------------------------------------------------------------------------------
   ! to verify that we don't extrapolate, make sure isopc is bound between -1 and 1 (effectively nearest neighbor)
   !-------------------------------------------------------------------------------------------------
   isopc = min( 1.0_SiKi, isopc )
   isopc = max(-1.0_SiKi, isopc )

   Indx_Hi = min( Indx_Lo + 1, nMax  )     ! make sure it's a valid index, zero-based

end subroutine SetTimeIndex


!====================================================================================================
!> This routine sets up interpolation of a 3-d or 4-d dataset.
!! This method is described here: http://rjwagner49.com/Mathematics/Interpolation.pdf
subroutine WaveField_Interp_Setup4D( Time, Position, p, m, ErrStat, ErrMsg )
   real(DbKi),                          intent(in   )  :: Time              !< time from the start of the simulation
   real(ReKi),                          intent(in   )  :: Position(3)       !< Array of XYZ coordinates, 3
   type(SeaSt_WaveField_ParameterType), intent(in   )  :: p                 !< Parameters
   type(SeaSt_WaveField_MiscVarType),   intent(inout)  :: m                 !< MiscVars
   integer(IntKi),                      intent(  out)  :: ErrStat           !< Error status
   character(*),                        intent(  out)  :: ErrMsg            !< Error message if ErrStat /= ErrID_None


   character(*), parameter              :: RoutineName = 'WaveField_Interp_Setup4D'
   integer(IntKi)                       :: i
   real(SiKi)                           :: isopc(4)           ! isoparametric coordinates
   integer(IntKi)                       :: ErrStat2
   character(ErrMsgLen)                 :: ErrMsg2

   ErrStat = ErrID_None
   ErrMsg  = ""

   ! Find the bounding indices for time
   call SetTimeIndex(Time, p%delta(1), p%n(1), m%Indx_Lo(1), m%Indx_Hi(1), isopc(1), ErrStat2, ErrMsg2)
   if (Failed()) return;

   ! Find the bounding indices for XY position
   do i=2,3  ! x and y components
      call SetCartesianXYIndex(Position(i-1), p%pZero(i), p%delta(i), p%n(i), m%Indx_Lo(i), m%Indx_Hi(i), isopc(i), m%FirstWarn_Clamp, ErrStat2, ErrMsg2)
      if (Failed()) return;
   enddo

   ! Find the bounding indices for Z position
   i=4 ! z component
   if (p%Z_Depth>0) then
      call SetCartesianZIndex(Position(i-1), p%Z_Depth, p%delta(i), p%n(i), m%Indx_Lo(i), m%Indx_Hi(i), isopc(i), m%FirstWarn_Clamp, ErrStat2, ErrMsg2)
      if (Failed()) return;
   else ! Regular z-grid
      call SetCartesianXYIndex(Position(i-1), p%pZero(i), p%delta(i), p%n(i), m%Indx_Lo(i), m%Indx_Hi(i), isopc(i), m%FirstWarn_Clamp, ErrStat2, ErrMsg2)
      if (Failed()) return;
   end if

   ! compute weighting factors
   m%N4D( 1) = ( 1.0_SiKi - isopc(1) ) * ( 1.0_SiKi - isopc(2) ) * ( 1.0_SiKi - isopc(3) ) * ( 1.0_SiKi - isopc(4) )
   m%N4D( 2) = ( 1.0_SiKi - isopc(1) ) * ( 1.0_SiKi - isopc(2) ) * ( 1.0_SiKi - isopc(3) ) * ( 1.0_SiKi + isopc(4) )
   m%N4D( 3) = ( 1.0_SiKi - isopc(1) ) * ( 1.0_SiKi - isopc(2) ) * ( 1.0_SiKi + isopc(3) ) * ( 1.0_SiKi - isopc(4) )
   m%N4D( 4) = ( 1.0_SiKi - isopc(1) ) * ( 1.0_SiKi - isopc(2) ) * ( 1.0_SiKi + isopc(3) ) * ( 1.0_SiKi + isopc(4) )
   m%N4D( 5) = ( 1.0_SiKi - isopc(1) ) * ( 1.0_SiKi + isopc(2) ) * ( 1.0_SiKi - isopc(3) ) * ( 1.0_SiKi - isopc(4) )
   m%N4D( 6) = ( 1.0_SiKi - isopc(1) ) * ( 1.0_SiKi + isopc(2) ) * ( 1.0_SiKi - isopc(3) ) * ( 1.0_SiKi + isopc(4) )
   m%N4D( 7) = ( 1.0_SiKi - isopc(1) ) * ( 1.0_SiKi + isopc(2) ) * ( 1.0_SiKi + isopc(3) ) * ( 1.0_SiKi - isopc(4) )
   m%N4D( 8) = ( 1.0_SiKi - isopc(1) ) * ( 1.0_SiKi + isopc(2) ) * ( 1.0_SiKi + isopc(3) ) * ( 1.0_SiKi + isopc(4) )
   m%N4D( 9) = ( 1.0_SiKi + isopc(1) ) * ( 1.0_SiKi - isopc(2) ) * ( 1.0_SiKi - isopc(3) ) * ( 1.0_SiKi - isopc(4) )
   m%N4D(10) = ( 1.0_SiKi + isopc(1) ) * ( 1.0_SiKi - isopc(2) ) * ( 1.0_SiKi - isopc(3) ) * ( 1.0_SiKi + isopc(4) )
   m%N4D(11) = ( 1.0_SiKi + isopc(1) ) * ( 1.0_SiKi - isopc(2) ) * ( 1.0_SiKi + isopc(3) ) * ( 1.0_SiKi - isopc(4) )
   m%N4D(12) = ( 1.0_SiKi + isopc(1) ) * ( 1.0_SiKi - isopc(2) ) * ( 1.0_SiKi + isopc(3) ) * ( 1.0_SiKi + isopc(4) )
   m%N4D(13) = ( 1.0_SiKi + isopc(1) ) * ( 1.0_SiKi + isopc(2) ) * ( 1.0_SiKi - isopc(3) ) * ( 1.0_SiKi - isopc(4) )
   m%N4D(14) = ( 1.0_SiKi + isopc(1) ) * ( 1.0_SiKi + isopc(2) ) * ( 1.0_SiKi - isopc(3) ) * ( 1.0_SiKi + isopc(4) )
   m%N4D(15) = ( 1.0_SiKi + isopc(1) ) * ( 1.0_SiKi + isopc(2) ) * ( 1.0_SiKi + isopc(3) ) * ( 1.0_SiKi - isopc(4) )
   m%N4D(16) = ( 1.0_SiKi + isopc(1) ) * ( 1.0_SiKi + isopc(2) ) * ( 1.0_SiKi + isopc(3) ) * ( 1.0_SiKi + isopc(4) )
   m%N4D     = m%N4D / REAL( SIZE(m%N4D), SiKi )  ! normalize

contains
   logical function Failed()
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
   end function
END Subroutine WaveField_Interp_Setup4D


subroutine WaveField_Interp_Setup3D( Time, Position, p, m, ErrStat, ErrMsg )
   real(DbKi),                          intent(in   )  :: Time              !< time from the start of the simulation
   real(ReKi),                          intent(in   )  :: Position(2)       !< Array of XYZ coordinates, 3
   type(SeaSt_WaveField_ParameterType), intent(in   )  :: p                 !< Parameters
   type(SeaSt_WaveField_MiscVarType),   intent(inout)  :: m                 !< MiscVars
   integer(IntKi),                      intent(  out)  :: ErrStat           !< Error status
   character(*),                        intent(  out)  :: ErrMsg            !< Error message if ErrStat /= ErrID_None

   character(*), parameter              :: RoutineName = 'WaveField_Interp_Setup3D'
   integer(IntKi)                       :: i
   real(SiKi)                           :: isopc(4)           ! isoparametric coordinates
   integer(IntKi)                       :: ErrStat2
   character(ErrMsgLen)                 :: ErrMsg2

   ErrStat = ErrID_None
   ErrMsg  = ""

   ! Find the bounding indices for time
   call SetTimeIndex(Time, p%delta(1), p%n(1), m%Indx_Lo(1), m%Indx_Hi(1), isopc(1), ErrStat2, ErrMsg2)
   if (Failed()) return;

   ! Find the bounding indices for XY position
   do i=2,3  ! x and y components
      call SetCartesianXYIndex(Position(i-1), p%pZero(i), p%delta(i), p%n(i), m%Indx_Lo(i), m%Indx_Hi(i), isopc(i), m%FirstWarn_Clamp, ErrStat2, ErrMsg2)
      if (Failed()) return;
   enddo

   ! compute weighting factors
   m%N3D(1)  = ( 1.0_ReKi + isopc(1) )*( 1.0_ReKi - isopc(2) )*( 1.0_ReKi - isopc(3) )
   m%N3D(2)  = ( 1.0_ReKi + isopc(1) )*( 1.0_ReKi + isopc(2) )*( 1.0_ReKi - isopc(3) )
   m%N3D(3)  = ( 1.0_ReKi - isopc(1) )*( 1.0_ReKi + isopc(2) )*( 1.0_ReKi - isopc(3) )
   m%N3D(4)  = ( 1.0_ReKi - isopc(1) )*( 1.0_ReKi - isopc(2) )*( 1.0_ReKi - isopc(3) )
   m%N3D(5)  = ( 1.0_ReKi + isopc(1) )*( 1.0_ReKi - isopc(2) )*( 1.0_ReKi + isopc(3) )
   m%N3D(6)  = ( 1.0_ReKi + isopc(1) )*( 1.0_ReKi + isopc(2) )*( 1.0_ReKi + isopc(3) )
   m%N3D(7)  = ( 1.0_ReKi - isopc(1) )*( 1.0_ReKi + isopc(2) )*( 1.0_ReKi + isopc(3) )
   m%N3D(8)  = ( 1.0_ReKi - isopc(1) )*( 1.0_ReKi - isopc(2) )*( 1.0_ReKi + isopc(3) )
   m%N3D     = m%N3D / REAL( SIZE(m%N3D), ReKi )  ! normalize

contains
   logical function Failed()
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
   end function
END Subroutine WaveField_Interp_Setup3D


!====================================================================================================
!> This routine interpolates a 4-d dataset.
!! This method is described here: http://rjwagner49.com/Mathematics/WaveFieldolation.pdf
function WaveField_Interp_4D( pKinXX, m )
   real(SiKi),                         intent(in   )  :: pKinXX(0:,:,:,:)
   type(SeaSt_WaveField_MiscVarType),  intent(in   )  :: m

   real(SiKi)                          :: WaveField_Interp_4D
   real(SiKi)                          :: u(16)    ! size 2^n

   ! interpolate
   u( 1) = pKinXX( m%Indx_Lo(1), m%Indx_Lo(2), m%Indx_Lo(3), m%Indx_Lo(4) )
   u( 2) = pKinXX( m%Indx_Lo(1), m%Indx_Lo(2), m%Indx_Lo(3), m%Indx_Hi(4) )
   u( 3) = pKinXX( m%Indx_Lo(1), m%Indx_Lo(2), m%Indx_Hi(3), m%Indx_Lo(4) )
   u( 4) = pKinXX( m%Indx_Lo(1), m%Indx_Lo(2), m%Indx_Hi(3), m%Indx_Hi(4) )
   u( 5) = pKinXX( m%Indx_Lo(1), m%Indx_Hi(2), m%Indx_Lo(3), m%Indx_Lo(4) )
   u( 6) = pKinXX( m%Indx_Lo(1), m%Indx_Hi(2), m%Indx_Lo(3), m%Indx_Hi(4) )
   u( 7) = pKinXX( m%Indx_Lo(1), m%Indx_Hi(2), m%Indx_Hi(3), m%Indx_Lo(4) )
   u( 8) = pKinXX( m%Indx_Lo(1), m%Indx_Hi(2), m%Indx_Hi(3), m%Indx_Hi(4) )
   u( 9) = pKinXX( m%Indx_Hi(1), m%Indx_Lo(2), m%Indx_Lo(3), m%Indx_Lo(4) )
   u(10) = pKinXX( m%Indx_Hi(1), m%Indx_Lo(2), m%Indx_Lo(3), m%Indx_Hi(4) )
   u(11) = pKinXX( m%Indx_Hi(1), m%Indx_Lo(2), m%Indx_Hi(3), m%Indx_Lo(4) )
   u(12) = pKinXX( m%Indx_Hi(1), m%Indx_Lo(2), m%Indx_Hi(3), m%Indx_Hi(4) )
   u(13) = pKinXX( m%Indx_Hi(1), m%Indx_Hi(2), m%Indx_Lo(3), m%Indx_Lo(4) )
   u(14) = pKinXX( m%Indx_Hi(1), m%Indx_Hi(2), m%Indx_Lo(3), m%Indx_Hi(4) )
   u(15) = pKinXX( m%Indx_Hi(1), m%Indx_Hi(2), m%Indx_Hi(3), m%Indx_Lo(4) )
   u(16) = pKinXX( m%Indx_Hi(1), m%Indx_Hi(2), m%Indx_Hi(3), m%Indx_Hi(4) )
   WaveField_Interp_4D = SUM ( m%N4D * u )
end function WaveField_Interp_4D


!====================================================================================================
!> This routine interpolates a 4-d dataset.
!! This method is described here: http://rjwagner49.com/Mathematics/Interpolation.pdf
function WaveField_Interp_4D_Vec( pKinXX, m)
   real(SiKi),                         intent(in   )  :: pKinXX(0:,:,:,:,:)
   type(SeaSt_WaveField_MiscVarType),  intent(in   )  :: m                    !< misc vars for interpolation

   real(SiKi)                                         :: WaveField_Interp_4D_Vec(3)
   real(SiKi)                                         :: u(16)   ! size 2^n
   integer(IntKi)                                     :: iDir

   ! interpolate
   do iDir = 1,3
      u( 1) = pKinXX( m%Indx_Lo(1), m%Indx_Lo(2), m%Indx_Lo(3), m%Indx_Lo(4), iDir )
      u( 2) = pKinXX( m%Indx_Lo(1), m%Indx_Lo(2), m%Indx_Lo(3), m%Indx_Hi(4), iDir )
      u( 3) = pKinXX( m%Indx_Lo(1), m%Indx_Lo(2), m%Indx_Hi(3), m%Indx_Lo(4), iDir )
      u( 4) = pKinXX( m%Indx_Lo(1), m%Indx_Lo(2), m%Indx_Hi(3), m%Indx_Hi(4), iDir )
      u( 5) = pKinXX( m%Indx_Lo(1), m%Indx_Hi(2), m%Indx_Lo(3), m%Indx_Lo(4), iDir )
      u( 6) = pKinXX( m%Indx_Lo(1), m%Indx_Hi(2), m%Indx_Lo(3), m%Indx_Hi(4), iDir )
      u( 7) = pKinXX( m%Indx_Lo(1), m%Indx_Hi(2), m%Indx_Hi(3), m%Indx_Lo(4), iDir )
      u( 8) = pKinXX( m%Indx_Lo(1), m%Indx_Hi(2), m%Indx_Hi(3), m%Indx_Hi(4), iDir )
      u( 9) = pKinXX( m%Indx_Hi(1), m%Indx_Lo(2), m%Indx_Lo(3), m%Indx_Lo(4), iDir )
      u(10) = pKinXX( m%Indx_Hi(1), m%Indx_Lo(2), m%Indx_Lo(3), m%Indx_Hi(4), iDir )
      u(11) = pKinXX( m%Indx_Hi(1), m%Indx_Lo(2), m%Indx_Hi(3), m%Indx_Lo(4), iDir )
      u(12) = pKinXX( m%Indx_Hi(1), m%Indx_Lo(2), m%Indx_Hi(3), m%Indx_Hi(4), iDir )
      u(13) = pKinXX( m%Indx_Hi(1), m%Indx_Hi(2), m%Indx_Lo(3), m%Indx_Lo(4), iDir )
      u(14) = pKinXX( m%Indx_Hi(1), m%Indx_Hi(2), m%Indx_Lo(3), m%Indx_Hi(4), iDir )
      u(15) = pKinXX( m%Indx_Hi(1), m%Indx_Hi(2), m%Indx_Hi(3), m%Indx_Lo(4), iDir )
      u(16) = pKinXX( m%Indx_Hi(1), m%Indx_Hi(2), m%Indx_Hi(3), m%Indx_Hi(4), iDir )
      WaveField_Interp_4D_Vec(iDir) = SUM ( m%N4D * u )
   end do
END FUNCTION WaveField_Interp_4D_Vec


!====================================================================================================
!> This routine interpolates a 4-d dataset.
!! This method is described here: http://rjwagner49.com/Mathematics/Interpolation.pdf
function WaveField_Interp_4D_Vec6( pKinXX, m)
   real(SiKi),                         intent(in   )  :: pKinXX(0:,:,:,:,:)
   type(SeaSt_WaveField_MiscVarType),  intent(in   )  :: m                    !< misc vars for interpolation

   real(SiKi)                                         :: WaveField_Interp_4D_Vec6(6)
   real(SiKi)                                         :: u(16)   ! size 2^n
   integer(IntKi)                                     :: iDir

   ! interpolate
   do iDir = 1,6
      u( 1) = pKinXX( m%Indx_Lo(1), m%Indx_Lo(2), m%Indx_Lo(3), m%Indx_Lo(4), iDir )
      u( 2) = pKinXX( m%Indx_Lo(1), m%Indx_Lo(2), m%Indx_Lo(3), m%Indx_Hi(4), iDir )
      u( 3) = pKinXX( m%Indx_Lo(1), m%Indx_Lo(2), m%Indx_Hi(3), m%Indx_Lo(4), iDir )
      u( 4) = pKinXX( m%Indx_Lo(1), m%Indx_Lo(2), m%Indx_Hi(3), m%Indx_Hi(4), iDir )
      u( 5) = pKinXX( m%Indx_Lo(1), m%Indx_Hi(2), m%Indx_Lo(3), m%Indx_Lo(4), iDir )
      u( 6) = pKinXX( m%Indx_Lo(1), m%Indx_Hi(2), m%Indx_Lo(3), m%Indx_Hi(4), iDir )
      u( 7) = pKinXX( m%Indx_Lo(1), m%Indx_Hi(2), m%Indx_Hi(3), m%Indx_Lo(4), iDir )
      u( 8) = pKinXX( m%Indx_Lo(1), m%Indx_Hi(2), m%Indx_Hi(3), m%Indx_Hi(4), iDir )
      u( 9) = pKinXX( m%Indx_Hi(1), m%Indx_Lo(2), m%Indx_Lo(3), m%Indx_Lo(4), iDir )
      u(10) = pKinXX( m%Indx_Hi(1), m%Indx_Lo(2), m%Indx_Lo(3), m%Indx_Hi(4), iDir )
      u(11) = pKinXX( m%Indx_Hi(1), m%Indx_Lo(2), m%Indx_Hi(3), m%Indx_Lo(4), iDir )
      u(12) = pKinXX( m%Indx_Hi(1), m%Indx_Lo(2), m%Indx_Hi(3), m%Indx_Hi(4), iDir )
      u(13) = pKinXX( m%Indx_Hi(1), m%Indx_Hi(2), m%Indx_Lo(3), m%Indx_Lo(4), iDir )
      u(14) = pKinXX( m%Indx_Hi(1), m%Indx_Hi(2), m%Indx_Lo(3), m%Indx_Hi(4), iDir )
      u(15) = pKinXX( m%Indx_Hi(1), m%Indx_Hi(2), m%Indx_Hi(3), m%Indx_Lo(4), iDir )
      u(16) = pKinXX( m%Indx_Hi(1), m%Indx_Hi(2), m%Indx_Hi(3), m%Indx_Hi(4), iDir )
      WaveField_Interp_4D_Vec6(iDir) = SUM ( m%N4D * u )
   end do
END FUNCTION WaveField_Interp_4D_Vec6


!====================================================================================================
!> This routine interpolates a 3-d dataset with index 1 = time (zero-based indexing), 2 = x-coordinate (1-based indexing), 3 = y-coordinate (1-based indexing)
!! This method is described here: http://rjwagner49.com/Mathematics/Interpolation.pdf
!FIXME: do like the above and call the WaveField_Interp_Setup3D routine ahead
function WaveField_Interp_3D( pKinXX, m )
   real(SiKi),                            intent(in   )  :: pKinXX(0:,:,:)    !< 3D Wave elevation data (SiKi for storage space reasons)
   type(SeaSt_WaveField_MiscVarType),     intent(inout)  :: m                 !< MiscVars

   character(*), parameter                :: RoutineName = 'WaveField_Interp_3D'
   real(SiKi)                             :: WaveField_Interp_3D
   real(SiKi)                             :: u(8)
   integer(IntKi)                         :: i

   ! interpolate
   u(1)  = pKinXX( m%Indx_Hi(1), m%Indx_Lo(2), m%Indx_Lo(3) )
   u(2)  = pKinXX( m%Indx_Hi(1), m%Indx_Hi(2), m%Indx_Lo(3) )
   u(3)  = pKinXX( m%Indx_Lo(1), m%Indx_Hi(2), m%Indx_Lo(3) )
   u(4)  = pKinXX( m%Indx_Lo(1), m%Indx_Lo(2), m%Indx_Lo(3) )
   u(5)  = pKinXX( m%Indx_Hi(1), m%Indx_Lo(2), m%Indx_Hi(3) )
   u(6)  = pKinXX( m%Indx_Hi(1), m%Indx_Hi(2), m%Indx_Hi(3) )
   u(7)  = pKinXX( m%Indx_Lo(1), m%Indx_Hi(2), m%Indx_Hi(3) )
   u(8)  = pKinXX( m%Indx_Lo(1), m%Indx_Lo(2), m%Indx_Hi(3) )
   WaveField_Interp_3D = SUM ( m%N3D * u )
end function WaveField_Interp_3D


FUNCTION WaveField_Interp_3D_VEC( pKinXX, m )
   real(SiKi),                            intent(in   )  :: pKinXX(0:,:,:,:)  !< 3D Wave excitation data (SiKi for storage space reasons)
   type(SeaSt_WaveField_MiscVarType),     intent(inout)  :: m                 !< MiscVars

   character(*), parameter                :: RoutineName = 'WaveField_Interp_3D_VEC'
   real(SiKi)                             :: WaveField_Interp_3D_VEC(3)
   real(SiKi)                             :: u(8)
   integer(IntKi)                         :: i

   ! interpolate
   do i = 1,3
      u(1)  = pKinXX( m%Indx_Hi(1), m%Indx_Lo(2), m%Indx_Lo(3), i )
      u(2)  = pKinXX( m%Indx_Hi(1), m%Indx_Hi(2), m%Indx_Lo(3), i )
      u(3)  = pKinXX( m%Indx_Lo(1), m%Indx_Hi(2), m%Indx_Lo(3), i )
      u(4)  = pKinXX( m%Indx_Lo(1), m%Indx_Lo(2), m%Indx_Lo(3), i )
      u(5)  = pKinXX( m%Indx_Hi(1), m%Indx_Lo(2), m%Indx_Hi(3), i )
      u(6)  = pKinXX( m%Indx_Hi(1), m%Indx_Hi(2), m%Indx_Hi(3), i )
      u(7)  = pKinXX( m%Indx_Lo(1), m%Indx_Hi(2), m%Indx_Hi(3), i )
      u(8)  = pKinXX( m%Indx_Lo(1), m%Indx_Lo(2), m%Indx_Hi(3), i )
      WaveField_Interp_3D_VEC(i) = SUM ( m%N3D * u )
   end do
end function WaveField_Interp_3D_VEC


function Wavefield_Interp_3D_VEC6( pKinXX, m )
   real(SiKi),                            intent(in   )  :: pKinXX(0:,:,:,:)  !< 3D Wave excitation data (SiKi for storage space reasons)
   type(SeaSt_WaveField_MiscVarType),     intent(inout)  :: m                 !< Miscvars

   character(*), parameter                :: RoutineName = 'Wavefield_Interp_3D_VEC6'
   real(SiKi)                             :: Wavefield_Interp_3D_VEC6(6)
   real(SiKi)                             :: u(8)
   integer(IntKi)                         :: i

   ! interpolate
   do i = 1,6
      u(1)  = pKinXX( m%Indx_Hi(1), m%Indx_Lo(2), m%Indx_Lo(3), i )
      u(2)  = pKinXX( m%Indx_Hi(1), m%Indx_Hi(2), m%Indx_Lo(3), i )
      u(3)  = pKinXX( m%Indx_Lo(1), m%Indx_Hi(2), m%Indx_Lo(3), i )
      u(4)  = pKinXX( m%Indx_Lo(1), m%Indx_Lo(2), m%Indx_Lo(3), i )
      u(5)  = pKinXX( m%Indx_Hi(1), m%Indx_Lo(2), m%Indx_Hi(3), i )
      u(6)  = pKinXX( m%Indx_Hi(1), m%Indx_Hi(2), m%Indx_Hi(3), i )
      u(7)  = pKinXX( m%Indx_Lo(1), m%Indx_Hi(2), m%Indx_Hi(3), i )
      u(8)  = pKinXX( m%Indx_Lo(1), m%Indx_Lo(2), m%Indx_Hi(3), i )
      Wavefield_Interp_3D_VEC6(i) = SUM ( m%N3D * u )
   end do
end function Wavefield_Interp_3D_VEC6



END MODULE SeaSt_WaveField
