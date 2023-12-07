MODULE SeaSt_WaveField

USE SeaState_Interp
USE SeaSt_WaveField_Types
USE IfW_FlowField, only: IfW_FlowField_GetVelAcc

IMPLICIT NONE
   
PRIVATE

! Public functions and subroutines
PUBLIC WaveField_GetNodeWaveElev1
PUBLIC WaveField_GetNodeWaveElev2
PUBLIC WaveField_GetNodeTotalWaveElev
PUBLIC WaveField_GetNodeWaveNormal
PUBLIC WaveField_GetNodeWaveKin
PUBLIC WaveField_GetNodeWaveVel
PUBLIC WaveField_GetNodeWaveVelAcc
PUBLIC WaveField_GetWaveKin
PUBLIC WaveField_GetWaveVelAcc_AD

CONTAINS

!-------------------- Subroutine for wave elevation ------------------!
FUNCTION WaveField_GetNodeWaveElev1( WaveField, SeaSt_Interp_m, Time, pos, ErrStat, ErrMsg )
   TYPE(SeaSt_WaveFieldType), INTENT( IN    ) :: WaveField
   TYPE(SeaSt_Interp_MiscVarType), INTENT(INOUT) :: SeaSt_Interp_m
   REAL(DbKi),      INTENT( IN    ) :: Time
   REAL(ReKi),      INTENT( IN    ) :: pos(*)  ! Position at which free-surface elevation is to be calculated. Third entry ignored if present.
   INTEGER(IntKi),  INTENT(   OUT ) :: ErrStat ! Error status of the operation
   CHARACTER(*),    INTENT(   OUT ) :: ErrMsg  ! Error message if errStat /= ErrID_None
  
   REAL(SiKi)                       :: WaveField_GetNodeWaveElev1
   REAL(SiKi)                       :: Zeta
   CHARACTER(*),    PARAMETER       :: RoutineName = 'WaveField_GetNodeWaveElev1'
   INTEGER(IntKi)                   :: errStat2
   CHARACTER(ErrMsgLen)             :: errMsg2
   
   ErrStat   = ErrID_None
   ErrMsg    = ""
   
   IF (ALLOCATED(WaveField%WaveElev1)) THEN
      Zeta = SeaSt_Interp_3D( Time, pos(1:2), WaveField%WaveElev1, WaveField%seast_interp_p, SeaSt_Interp_m%FirstWarn_Clamp, ErrStat2, ErrMsg2 )
        CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   ELSE
      Zeta = 0.0_SiKi
   END IF
   
   WaveField_GetNodeWaveElev1 = Zeta

END FUNCTION WaveField_GetNodeWaveElev1

FUNCTION WaveField_GetNodeWaveElev2( WaveField, SeaSt_Interp_m, Time, pos, ErrStat, ErrMsg )
   TYPE(SeaSt_WaveFieldType), INTENT( IN    ) :: WaveField
   TYPE(SeaSt_Interp_MiscVarType), INTENT(INOUT) :: SeaSt_Interp_m
   REAL(DbKi),      INTENT( IN    ) :: Time
   REAL(ReKi),      INTENT( IN    ) :: pos(*)  ! Position at which free-surface elevation is to be calculated. Third entry ignored if present.
   INTEGER(IntKi),  INTENT(   OUT ) :: ErrStat ! Error status of the operation
   CHARACTER(*),    INTENT(   OUT ) :: ErrMsg  ! Error message if errStat /= ErrID_None
   
   REAL(SiKi)                       :: WaveField_GetNodeWaveElev2
   REAL(SiKi)                       :: Zeta
   CHARACTER(*),    PARAMETER       :: RoutineName = 'WaveField_GetNodeWaveElev2'
   INTEGER(IntKi)                   :: errStat2
   CHARACTER(ErrMsgLen)             :: errMsg2
   
   ErrStat   = ErrID_None
   ErrMsg    = ""
   
   IF (ALLOCATED(WaveField%WaveElev2)) THEN
      Zeta = SeaSt_Interp_3D( Time, pos(1:2), WaveField%WaveElev2, WaveField%seast_interp_p, SeaSt_Interp_m%FirstWarn_Clamp, ErrStat2, ErrMsg2 )
        CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   ELSE
      Zeta = 0.0_SiKi
   END IF

   WaveField_GetNodeWaveElev2 = Zeta

END FUNCTION WaveField_GetNodeWaveElev2

FUNCTION WaveField_GetNodeTotalWaveElev( WaveField, SeaSt_Interp_m, Time, pos, ErrStat, ErrMsg )
   TYPE(SeaSt_WaveFieldType), INTENT( IN    ) :: WaveField
   TYPE(SeaSt_Interp_MiscVarType), INTENT(INOUT) :: SeaSt_Interp_m
   REAL(DbKi),      INTENT( IN    ) :: Time
   REAL(ReKi),      INTENT( IN    ) :: pos(*)  ! Position at which free-surface elevation is to be calculated. Third entry ignored if present.
   INTEGER(IntKi),  INTENT(   OUT ) :: ErrStat ! Error status of the operation
   CHARACTER(*),    INTENT(   OUT ) :: ErrMsg  ! Error message if errStat /= ErrID_None

   REAL(SiKi)                       :: WaveField_GetNodeTotalWaveElev
   REAL(SiKi)                       :: Zeta1, Zeta2
   CHARACTER(*),    PARAMETER       :: RoutineName = 'WaveField_GetNodeTotalWaveElev'
   INTEGER(IntKi)                   :: errStat2
   CHARACTER(ErrMsgLen)             :: errMsg2

   ErrStat   = ErrID_None
   ErrMsg    = ""
   
   Zeta1 = WaveField_GetNodeWaveElev1( WaveField, SeaSt_Interp_m, Time, pos, ErrStat2, ErrMsg2 )
     CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   Zeta2 = WaveField_GetNodeWaveElev2( WaveField, SeaSt_Interp_m, Time, pos, ErrStat2, ErrMsg2 )
     CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   WaveField_GetNodeTotalWaveElev = Zeta1 + Zeta2
   
END FUNCTION WaveField_GetNodeTotalWaveElev

SUBROUTINE WaveField_GetNodeWaveNormal( WaveField, SeaSt_Interp_m, Time, pos, r, n, ErrStat, ErrMsg )

   TYPE(SeaSt_WaveFieldType), INTENT( IN    ) :: WaveField
   TYPE(SeaSt_Interp_MiscVarType), INTENT(INOUT) :: SeaSt_Interp_m
   REAL(DbKi),                INTENT( IN    ) :: Time
   REAL(ReKi),                INTENT( IN    ) :: pos(*)  ! Position at which free-surface normal is to be calculated. Third entry ignored if present.
   REAL(ReKi),                INTENT( IN    ) :: r       ! Distance for central differencing
   REAL(ReKi),                INTENT(   OUT ) :: n(3)    ! Free-surface normal vector
   INTEGER(IntKi),            INTENT(   OUT ) :: ErrStat ! Error status of the operation
   CHARACTER(*),              INTENT(   OUT ) :: ErrMsg  ! Error message if errStat /= ErrID_None
   REAL(SiKi)                                 :: ZetaP,ZetaM
   REAL(ReKi)                                 :: r1,dZetadx,dZetady
   CHARACTER(*),              PARAMETER       :: RoutineName = 'WaveField_GetNodeWaveNormal'
   INTEGER(IntKi)                             :: errStat2
   CHARACTER(ErrMsgLen)                       :: errMsg2
   ErrStat   = ErrID_None
   ErrMsg    = ""

   r1 = MAX(r,real(1.0e-6,ReKi)) ! In case r is zero

   ZetaP = WaveField_GetNodeTotalWaveElev( WaveField, SeaSt_Interp_m, Time, (/pos(1)+r1,pos(2)/), ErrStat2, ErrMsg2 )
     CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   ZetaM = WaveField_GetNodeTotalWaveElev( WaveField, SeaSt_Interp_m, Time, (/pos(1)-r1,pos(2)/), ErrStat2, ErrMsg2 )
     CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   dZetadx = REAL(ZetaP-ZetaM,ReKi)/(2.0_ReKi*r1)
      
   ZetaP = WaveField_GetNodeTotalWaveElev( WaveField, SeaSt_Interp_m, Time, (/pos(1),pos(2)+r1/), ErrStat2, ErrMsg2 )
     CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   ZetaM = WaveField_GetNodeTotalWaveElev( WaveField, SeaSt_Interp_m, Time, (/pos(1),pos(2)-r1/), ErrStat2, ErrMsg2 )
     CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   dZetady = REAL(ZetaP-ZetaM,ReKi)/(2.0_ReKi*r1)
      
   n = (/-dZetadx,-dZetady,1.0_ReKi/)
   n = n / SQRT(Dot_Product(n,n))

END SUBROUTINE WaveField_GetNodeWaveNormal

!-------------------- Subroutine for full wave field kinematics --------------------!
SUBROUTINE WaveField_GetNodeWaveKin( WaveField, SeaSt_Interp_m, Time, pos, forceNodeInWater, fetchDynCurrent, nodeInWater, WaveElev1, WaveElev2, WaveElev, FDynP, FV, FA, FAMCF, ErrStat, ErrMsg )
   TYPE(SeaSt_WaveFieldType), INTENT( IN    ) :: WaveField
   TYPE(SeaSt_Interp_MiscVarType), INTENT( INOUT ) :: SeaSt_Interp_m
   REAL(DbKi),                INTENT( IN    ) :: Time
   REAL(ReKi),                INTENT( IN    ) :: pos(3)
   LOGICAL,                   INTENT( IN    ) :: forceNodeInWater
   LOGICAL,                   INTENT( IN    ) :: fetchDynCurrent
   REAL(SiKi),                INTENT(   OUT ) :: WaveElev1
   REAL(SiKi),                INTENT(   OUT ) :: WaveElev2
   REAL(SiKi),                INTENT(   OUT ) :: WaveElev
   REAL(SiKi),                INTENT(   OUT ) :: FV(3)
   REAL(SiKi),                INTENT(   OUT ) :: FA(3)
   REAL(SiKi),                INTENT(   OUT ) :: FAMCF(3)
   REAL(SiKi),                INTENT(   OUT ) :: FDynP
   INTEGER(IntKi),            INTENT(   OUT ) :: nodeInWater

   INTEGER(IntKi),            INTENT(   OUT ) :: ErrStat ! Error status of the operation
   CHARACTER(*),              INTENT(   OUT ) :: ErrMsg  ! Error message if errStat /= ErrID_None

   REAL(ReKi)                                 :: posXY(2), posPrime(3), posXY0(3), PosOffset(3), posDummy(3,1)
   INTEGER(IntKi)                             :: startNode
   REAL(ReKi), allocatable                    :: FV_DC(:,:), FA_DC(:,:)
   CHARACTER(*),              PARAMETER       :: RoutineName = 'WaveField_GetNodeWaveKin'
   INTEGER(IntKi)                             :: errStat2
   CHARACTER(ErrMsgLen)                       :: errMsg2

   ErrStat   = ErrID_None
   ErrMsg    = ""

   posXY    = pos(1:2)
   posXY0   = (/pos(1),pos(2),0.0_ReKi/)
   FAMCF(:) = 0.0

   ! Wave elevation
   WaveElev1 = WaveField_GetNodeWaveElev1( WaveField, SeaSt_Interp_m, Time, pos, ErrStat2, ErrMsg2 )
     CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   WaveElev2 = WaveField_GetNodeWaveElev2( WaveField, SeaSt_Interp_m, Time, pos, ErrStat2, ErrMsg2 )
     CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   WaveElev  = WaveElev1 + WaveElev2
    
   IF (WaveField%WaveStMod == 0) THEN ! No wave stretching
    
      IF ( pos(3) <= 0.0_ReKi) THEN ! Node is at or below the SWL
         nodeInWater = 1_IntKi
         ! Use location to obtain interpolated values of kinematics         
         CALL SeaSt_Interp_Setup( Time, pos, WaveField%seast_interp_p, SeaSt_Interp_m, ErrStat2, ErrMsg2 ) 
           CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         FV(:) = SeaSt_Interp_4D_Vec( WaveField%WaveVel,  SeaSt_Interp_m, ErrStat2, ErrMsg2 )
           CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         FA(:) = SeaSt_Interp_4D_Vec( WaveField%WaveAcc,  SeaSt_Interp_m, ErrStat2, ErrMsg2 )
           CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         FDynP = SeaSt_Interp_4D    ( WaveField%WaveDynP, SeaSt_Interp_m, ErrStat2, ErrMsg2 )
           CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF ( ALLOCATED(WaveField%WaveAccMCF) ) THEN
            FAMCF(:) = SeaSt_Interp_4D_Vec( WaveField%WaveAccMCF, SeaSt_Interp_m, ErrStat2, ErrMsg2 )
              CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
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
               CALL SeaSt_Interp_Setup( Time, pos, WaveField%seast_interp_p, SeaSt_Interp_m, ErrStat2, ErrMsg2 ) 
                 CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               FV(:) = SeaSt_Interp_4D_Vec( WaveField%WaveVel,  SeaSt_Interp_m, ErrStat2, ErrMsg2 )
                 CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               FA(:) = SeaSt_Interp_4D_Vec( WaveField%WaveAcc,  SeaSt_Interp_m, ErrStat2, ErrMsg2 )
                 CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               FDynP = SeaSt_Interp_4D    ( WaveField%WaveDynP, SeaSt_Interp_m, ErrStat2, ErrMsg2 )
                 CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               IF ( ALLOCATED(WaveField%WaveAccMCF) ) THEN
                  FAMCF(:) = SeaSt_Interp_4D_Vec( WaveField%WaveAccMCF, SeaSt_Interp_m, ErrStat2, ErrMsg2 )
                    CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               END IF

            ELSE ! Node is above SWL - need wave stretching
          
               ! Vertical wave stretching
               CALL SeaSt_Interp_Setup( Time, posXY0, WaveField%seast_interp_p, SeaSt_Interp_m, ErrStat2, ErrMsg2 ) 
                 CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               FV(:) = SeaSt_Interp_4D_vec( WaveField%WaveVel,  seast_interp_m, ErrStat2, ErrMsg2 )
                 CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               FA(:) = SeaSt_Interp_4D_vec( WaveField%WaveAcc,  seast_interp_m, ErrStat2, ErrMsg2 )
                 CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               FDynP = SeaSt_Interp_4D    ( WaveField%WaveDynP, seast_interp_m, ErrStat2, ErrMsg2 )
                 CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               IF ( ALLOCATED(WaveField%WaveAccMCF) ) THEN
                  FAMCF(:) = SeaSt_Interp_4D_vec( WaveField%WaveAccMCF, seast_interp_m, ErrStat2, ErrMsg2 )
                    CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               END IF
                      
               ! Extrapoled wave stretching
               IF (WaveField%WaveStMod == 2) THEN 
                  FV(:) = FV(:) + SeaSt_Interp_3D_vec( Time, posXY, WaveField%PWaveVel0,  WaveField%seast_interp_p, SeaSt_Interp_m%FirstWarn_Clamp, ErrStat2, ErrMsg2 ) * pos(3)
                    CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                  FA(:) = FA(:) + SeaSt_Interp_3D_vec( Time, posXY, WaveField%PWaveAcc0,  WaveField%seast_interp_p, SeaSt_Interp_m%FirstWarn_Clamp, ErrStat2, ErrMsg2 ) * pos(3)
                    CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                  FDynP = FDynP + SeaSt_Interp_3D    ( Time, posXY, WaveField%PWaveDynP0, WaveField%seast_interp_p, SeaSt_Interp_m%FirstWarn_Clamp, ErrStat2, ErrMsg2 ) * pos(3)
                    CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                  IF ( ALLOCATED(WaveField%WaveAccMCF) ) THEN
                     FAMCF(:) = FAMCF(:) + SeaSt_Interp_3D_vec( Time, posXY, WaveField%PWaveAccMCF0, WaveField%seast_interp_p, SeaSt_Interp_m%FirstWarn_Clamp, ErrStat2, ErrMsg2 ) * pos(3)
                       CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                  END IF
               END IF
          
            END IF ! Node is above or below SWL
 
         ELSE ! Wheeler stretching - no need to check whether the node is above or below SWL
                  
            ! Map the node z-position linearly from [-EffWtrDpth,m%WaveElev(j)] to [-EffWtrDpth,0] 
            posPrime    = pos
            posPrime(3) = WaveField%EffWtrDpth*(WaveField%EffWtrDpth+pos(3))/(WaveField%EffWtrDpth+WaveElev)-WaveField%EffWtrDpth
            posPrime(3) = MIN( posPrime(3), 0.0_ReKi) ! Clamp z-position to zero. Needed when forceNodeInWater=.TRUE. 
                  
            ! Obtain the wave-field variables by interpolation with the mapped position.
            CALL SeaSt_Interp_Setup( Time, posPrime, WaveField%seast_interp_p, seast_interp_m, ErrStat2, ErrMsg2 ) 
              CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            FV(:) = SeaSt_Interp_4D_Vec( WaveField%WaveVel,  seast_interp_m, ErrStat2, ErrMsg2 )
              CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            FA(:) = SeaSt_Interp_4D_Vec( WaveField%WaveAcc,  seast_interp_m, ErrStat2, ErrMsg2 )
              CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            FDynP = SeaSt_Interp_4D    ( WaveField%WaveDynP, seast_interp_m, ErrStat2, ErrMsg2 )
              CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            IF ( ALLOCATED(WaveField%WaveAccMCF) ) THEN
               FAMCF(:) = SeaSt_Interp_4D_Vec( WaveField%WaveAccMCF, seast_interp_m, ErrStat2, ErrMsg2 )
                 CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
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
   
   IF (fetchDynCurrent .AND. WaveField%hasCurrField) THEN
      startNode = -1
      PosOffset = (/0.0_ReKi,0.0_ReKi,WaveField%EffWtrDpth/)
      posDummy(:,1) = pos
      ALLOCATE(FV_DC(3,1), STAT=ErrStat2)
      if (ErrStat2 /= 0) then
         call SetErrStat( ErrID_Info, 'Error allocating FV_DC', ErrStat, ErrMsg, RoutineName )
         return
      end if     
      ALLOCATE(FA_DC(3,1), STAT=ErrStat2)
      if (ErrStat2 /= 0) then
         call SetErrStat( ErrID_Info, 'Error allocating FA_DC', ErrStat, ErrMsg, RoutineName )
         return
      end if     
      CALL IfW_FlowField_GetVelAcc(WaveField%CurrField, startNode, Time, posDummy, FV_DC, FA_DC, ErrStat2, ErrMsg2, PosOffset=PosOffset)
        CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      FV = FV + nodeInWater * FV_DC(:,1)
      FA = FA + nodeInWater * FA_DC(:,1)
   END IF

END SUBROUTINE WaveField_GetNodeWaveKin

!-------------------- Subroutine for wave field velocity only --------------------!
SUBROUTINE WaveField_GetNodeWaveVel( WaveField, SeaSt_Interp_m, Time, pos, forceNodeInWater, fetchDynCurrent, nodeInWater, FV, ErrStat, ErrMsg )
   TYPE(SeaSt_WaveFieldType), INTENT( IN    ) :: WaveField
   TYPE(SeaSt_Interp_MiscVarType), INTENT(INOUT) :: SeaSt_Interp_m
   REAL(DbKi),                INTENT( IN    ) :: Time
   REAL(ReKi),                INTENT( IN    ) :: pos(3)
   LOGICAL,                   INTENT( IN    ) :: forceNodeInWater
   LOGICAL,                   INTENT( IN    ) :: fetchDynCurrent
   INTEGER(IntKi),            INTENT(   OUT ) :: nodeInWater
   REAL(SiKi),                INTENT(   OUT ) :: FV(3)
   INTEGER(IntKi),            INTENT(   OUT ) :: ErrStat ! Error status of the operation
   CHARACTER(*),              INTENT(   OUT ) :: ErrMsg  ! Error message if errStat /= ErrID_None

   REAL(SiKi)                                 :: WaveElev
   REAL(ReKi)                                 :: posXY(2), posPrime(3), posXY0(3), PosOffset(3), posDummy(3,1)
   REAL(ReKi), allocatable                    :: FV_DC(:,:), FA_DC(:,:)
   INTEGER(IntKi)                             :: startNode
   CHARACTER(*),              PARAMETER       :: RoutineName = 'WaveField_GetNodeWaveVel'
   INTEGER(IntKi)                             :: errStat2
   CHARACTER(ErrMsgLen)                       :: errMsg2

   ErrStat   = ErrID_None
   ErrMsg    = ""

   posXY    = pos(1:2)
   posXY0   = (/pos(1),pos(2),0.0_ReKi/)

   ! Wave elevation
   WaveElev  = WaveField_GetNodeTotalWaveElev( WaveField, SeaSt_Interp_m, Time, pos, ErrStat2, ErrMsg2 )
     CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
    
   IF (WaveField%WaveStMod == 0) THEN ! No wave stretching
    
      IF ( pos(3) <= 0.0_ReKi) THEN ! Node is at or below the SWL
         nodeInWater = 1_IntKi
         ! Use location to obtain interpolated values of kinematics         
         CALL SeaSt_Interp_Setup( Time, pos, WaveField%seast_interp_p, seast_interp_m, ErrStat2, ErrMsg2 ) 
           CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         FV(:) = SeaSt_Interp_4D_Vec( WaveField%WaveVel,  seast_interp_m, ErrStat2, ErrMsg2 )
           CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
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
               CALL SeaSt_Interp_Setup( Time, pos, WaveField%seast_interp_p, seast_interp_m, ErrStat2, ErrMsg2 ) 
                 CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               FV(:) = SeaSt_Interp_4D_Vec( WaveField%WaveVel,  seast_interp_m, ErrStat2, ErrMsg2 )
                 CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

            ELSE ! Node is above SWL - need wave stretching
          
               ! Vertical wave stretching
               CALL SeaSt_Interp_Setup( Time, posXY0, WaveField%seast_interp_p, seast_interp_m, ErrStat2, ErrMsg2 ) 
                 CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               FV(:) = SeaSt_Interp_4D_vec( WaveField%WaveVel,  seast_interp_m, ErrStat2, ErrMsg2 )
                 CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                      
               ! Extrapoled wave stretching
               IF (WaveField%WaveStMod == 2) THEN 
                  FV(:) = FV(:) + SeaSt_Interp_3D_vec( Time, posXY, WaveField%PWaveVel0,  WaveField%seast_interp_p, SeaSt_Interp_m%FirstWarn_Clamp, ErrStat2, ErrMsg2 ) * pos(3)
                    CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               END IF
          
            END IF ! Node is above or below SWL
 
         ELSE ! Wheeler stretching - no need to check whether the node is above or below SWL
                  
            ! Map the node z-position linearly from [-EffWtrDpth,m%WaveElev(j)] to [-EffWtrDpth,0] 
            posPrime    = pos
            posPrime(3) = WaveField%EffWtrDpth*(WaveField%EffWtrDpth+pos(3))/(WaveField%EffWtrDpth+WaveElev)-WaveField%EffWtrDpth
            posPrime(3) = MIN( posPrime(3), 0.0_ReKi) ! Clamp z-position to zero. Needed when forceNodeInWater=.TRUE. 
                  
            ! Obtain the wave-field variables by interpolation with the mapped position.
            CALL SeaSt_Interp_Setup( Time, posPrime, WaveField%seast_interp_p, seast_interp_m, ErrStat2, ErrMsg2 ) 
              CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            FV(:) = SeaSt_Interp_4D_Vec( WaveField%WaveVel,  seast_interp_m, ErrStat2, ErrMsg2 )
              CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            
         END IF
        
      ELSE ! Node is out of water - zero-out all wave dynamics
          
         nodeInWater = 0_IntKi  
         FV(:)       = 0.0
          
      END IF ! If node is in or out of water
      
   END IF ! If wave stretching is on or off
   
   IF (fetchDynCurrent .AND. WaveField%hasCurrField) THEN
      startNode = -1
      PosOffset = (/0.0_ReKi,0.0_ReKi,WaveField%EffWtrDpth/)
      posDummy(:,1) = pos
      ALLOCATE(FV_DC(3,1), STAT=ErrStat2)
      if (ErrStat2 /= 0) then
         call SetErrStat( ErrID_Info, 'Error allocating FV_DC', ErrStat, ErrMsg, RoutineName )
         return
      end if   
      CALL IfW_FlowField_GetVelAcc(WaveField%CurrField, startNode, Time, posDummy, FV_DC, FA_DC, ErrStat2, ErrMsg2, PosOffset=PosOffset)
        CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      FV = FV + nodeInWater * FV_DC(:,1)
   END IF

END SUBROUTINE WaveField_GetNodeWaveVel

SUBROUTINE WaveField_GetNodeWaveVelAcc( WaveField, SeaSt_Interp_m, Time, pos, forceNodeInWater, fetchDynCurrent, nodeInWater, FV, FA, ErrStat, ErrMsg )
   TYPE(SeaSt_WaveFieldType), INTENT( IN    ) :: WaveField
   TYPE(SeaSt_Interp_MiscVarType), INTENT(INOUT) :: SeaSt_Interp_m
   REAL(DbKi),                INTENT( IN    ) :: Time
   REAL(ReKi),                INTENT( IN    ) :: pos(3)
   LOGICAL,                   INTENT( IN    ) :: forceNodeInWater
   LOGICAL,                   INTENT( IN    ) :: fetchDynCurrent
   REAL(SiKi),                INTENT(   OUT ) :: FV(3)
   REAL(SiKi),                INTENT(   OUT ) :: FA(3)
   INTEGER(IntKi),            INTENT(   OUT ) :: nodeInWater
   INTEGER(IntKi),            INTENT(   OUT ) :: ErrStat ! Error status of the operation
   CHARACTER(*),              INTENT(   OUT ) :: ErrMsg  ! Error message if errStat /= ErrID_None

   REAL(SiKi)                                 :: WaveElev
   REAL(ReKi)                                 :: posXY(2), posPrime(3), posXY0(3), PosOffset(3), posDummy(3,1)
   INTEGER(IntKi)                             :: startNode
   REAL(ReKi), allocatable                    :: FV_DC(:,:), FA_DC(:,:)
   CHARACTER(*),              PARAMETER       :: RoutineName = 'WaveField_GetNodeWaveVelAcc'
   INTEGER(IntKi)                             :: errStat2
   CHARACTER(ErrMsgLen)                       :: errMsg2

   ErrStat   = ErrID_None
   ErrMsg    = ""

   posXY    = pos(1:2)
   posXY0   = (/pos(1),pos(2),0.0_ReKi/)
   
   ! Wave elevation
   WaveElev  = WaveField_GetNodeTotalWaveElev( WaveField, SeaSt_Interp_m, Time, pos, ErrStat2, ErrMsg2 )
     CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
    
   IF (WaveField%WaveStMod == 0) THEN ! No wave stretching
    
      IF ( pos(3) <= 0.0_ReKi) THEN ! Node is at or below the SWL
         nodeInWater = 1_IntKi
         ! Use location to obtain interpolated values of kinematics         
         CALL SeaSt_Interp_Setup( Time, pos, WaveField%seast_interp_p, seast_interp_m, ErrStat2, ErrMsg2 ) 
           CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         FV(:) = SeaSt_Interp_4D_Vec( WaveField%WaveVel,  seast_interp_m, ErrStat2, ErrMsg2 )
           CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         FA(:) = SeaSt_Interp_4D_Vec( WaveField%WaveAcc,  seast_interp_m, ErrStat2, ErrMsg2 )
           CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
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
               CALL SeaSt_Interp_Setup( Time, pos, WaveField%seast_interp_p, seast_interp_m, ErrStat2, ErrMsg2 ) 
                 CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               FV(:) = SeaSt_Interp_4D_Vec( WaveField%WaveVel,  seast_interp_m, ErrStat2, ErrMsg2 )
                 CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               FA(:) = SeaSt_Interp_4D_Vec( WaveField%WaveAcc,  seast_interp_m, ErrStat2, ErrMsg2 )
                 CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

            ELSE ! Node is above SWL - need wave stretching
          
               ! Vertical wave stretching
               CALL SeaSt_Interp_Setup( Time, posXY0, WaveField%seast_interp_p, seast_interp_m, ErrStat2, ErrMsg2 ) 
                 CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               FV(:) = SeaSt_Interp_4D_vec( WaveField%WaveVel,  seast_interp_m, ErrStat2, ErrMsg2 )
                 CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               FA(:) = SeaSt_Interp_4D_vec( WaveField%WaveAcc,  seast_interp_m, ErrStat2, ErrMsg2 )
                 CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                      
               ! Extrapoled wave stretching
               IF (WaveField%WaveStMod == 2) THEN 
                  FV(:) = FV(:) + SeaSt_Interp_3D_vec( Time, posXY, WaveField%PWaveVel0,  WaveField%seast_interp_p, SeaSt_Interp_m%FirstWarn_Clamp, ErrStat2, ErrMsg2 ) * pos(3)
                    CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                  FA(:) = FA(:) + SeaSt_Interp_3D_vec( Time, posXY, WaveField%PWaveAcc0,  WaveField%seast_interp_p, SeaSt_Interp_m%FirstWarn_Clamp, ErrStat2, ErrMsg2 ) * pos(3)
                    CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               END IF
          
            END IF ! Node is above or below SWL
 
         ELSE ! Wheeler stretching - no need to check whether the node is above or below SWL
                  
            ! Map the node z-position linearly from [-EffWtrDpth,m%WaveElev(j)] to [-EffWtrDpth,0] 
            posPrime    = pos
            posPrime(3) = WaveField%EffWtrDpth*(WaveField%EffWtrDpth+pos(3))/(WaveField%EffWtrDpth+WaveElev)-WaveField%EffWtrDpth
            posPrime(3) = MIN( posPrime(3), 0.0_ReKi) ! Clamp z-position to zero. Needed when forceNodeInWater=.TRUE. 
                  
            ! Obtain the wave-field variables by interpolation with the mapped position.
            CALL SeaSt_Interp_Setup( Time, posPrime, WaveField%seast_interp_p, seast_interp_m, ErrStat2, ErrMsg2 ) 
              CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            FV(:) = SeaSt_Interp_4D_Vec( WaveField%WaveVel,  seast_interp_m, ErrStat2, ErrMsg2 )
              CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            FA(:) = SeaSt_Interp_4D_Vec( WaveField%WaveAcc,  seast_interp_m, ErrStat2, ErrMsg2 )
              CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         END IF ! Wave stretching method
        
      ELSE ! Node is out of water - zero-out all wave dynamics
          
         nodeInWater = 0_IntKi  
         FV(:)       = 0.0
         FA(:)       = 0.0
          
      END IF ! If node is in or out of water
      
   END IF ! If wave stretching is on or off
   
   IF (fetchDynCurrent .AND. WaveField%hasCurrField) THEN
      startNode = -1
      PosOffset = (/0.0_ReKi,0.0_ReKi,WaveField%EffWtrDpth/)
      posDummy(:,1) = pos
      ALLOCATE(FV_DC(3,1), STAT=ErrStat2)
      if (ErrStat2 /= 0) then
         call SetErrStat( ErrID_Info, 'Error allocating FV_DC', ErrStat, ErrMsg, RoutineName )
         return
      end if     
      ALLOCATE(FA_DC(3,1), STAT=ErrStat2)
      if (ErrStat2 /= 0) then
         call SetErrStat( ErrID_Info, 'Error allocating FA_DC', ErrStat, ErrMsg, RoutineName )
         return
      end if     
      CALL IfW_FlowField_GetVelAcc(WaveField%CurrField, startNode, Time, posDummy, FV_DC, FA_DC, ErrStat2, ErrMsg2, PosOffset=PosOffset)
        CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      FV = FV + nodeInWater * FV_DC(:,1)
      FA = FA + nodeInWater * FA_DC(:,1)
   END IF

END SUBROUTINE WaveField_GetNodeWaveVelAcc


SUBROUTINE WaveField_GetWaveKin( WaveField, SeaSt_Interp_m, Time, pos, forceNodeInWater, fetchDynCurrent, nodeInWater, WaveElev1, WaveElev2, WaveElev, FDynP, FV, FA, FAMCF, ErrStat, ErrMsg )
   TYPE(SeaSt_WaveFieldType), INTENT( IN    ) :: WaveField
   TYPE(SeaSt_Interp_MiscVarType), INTENT(INOUT) :: SeaSt_Interp_m
   REAL(DbKi),                INTENT( IN    ) :: Time
   REAL(ReKi),                INTENT( IN    ) :: pos(:,:)
   LOGICAL,                   INTENT( IN    ) :: forceNodeInWater
   LOGICAL,                   INTENT( IN    ) :: fetchDynCurrent
   REAL(SiKi),                INTENT(   OUT ) :: WaveElev1(:)
   REAL(SiKi),                INTENT(   OUT ) :: WaveElev2(:)
   REAL(SiKi),                INTENT(   OUT ) :: WaveElev(:)
   REAL(ReKi),                INTENT(   OUT ) :: FV(:,:)
   REAL(ReKi),                INTENT(   OUT ) :: FA(:,:)
   REAL(ReKi),                INTENT(   OUT ) :: FAMCF(:,:)
   REAL(ReKi),                INTENT(   OUT ) :: FDynP(:)
   INTEGER(IntKi),            INTENT(   OUT ) :: nodeInWater(:)
   INTEGER(IntKi),            INTENT(   OUT ) :: ErrStat ! Error status of the operation
   CHARACTER(*),              INTENT(   OUT ) :: ErrMsg  ! Error message if errStat /= ErrID_None

   CHARACTER(*),              PARAMETER       :: RoutineName = 'WaveField_GetWaveKin'
   INTEGER(IntKi)                             :: errStat2
   CHARACTER(ErrMsgLen)                       :: errMsg2

   INTEGER(IntKi)                             :: NumPoints, i, startNode
   REAL(SiKi)                                 :: FDynP_node, FV_node(3), FA_node(3), FAMCF_node(3)
   REAL(ReKi)                                 :: PosOffset(3)

   REAL(ReKi), allocatable                    :: FV_DC(:,:), FA_DC(:,:)

   ErrStat   = ErrID_None
   ErrMsg    = ""

   NumPoints = size(pos, dim=2)
   DO i = 1, NumPoints
      CALL WaveField_GetNodeWaveKin( WaveField, SeaSt_Interp_m, Time, pos(:,i), forceNodeInWater, .FALSE., nodeInWater(i), WaveElev1(i), WaveElev2(i), WaveElev(i), FDynP_node, FV_node, FA_node, FAMCF_node, ErrStat2, ErrMsg2 )
        CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
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
      ALLOCATE(FV_DC( 3, NumPoints ), STAT=ErrStat2)
      if (ErrStat2 /= 0) then
         call SetErrStat( ErrID_Info, 'Error allocating FV_DC', ErrStat, ErrMsg, RoutineName )
         return
      end if     
      ALLOCATE(FA_DC( 3, NumPoints ), STAT=ErrStat2)
      if (ErrStat2 /= 0) then
         call SetErrStat( ErrID_Info, 'Error allocating FA_DC', ErrStat, ErrMsg, RoutineName )
         return
      end if     
      CALL IfW_FlowField_GetVelAcc(WaveField%CurrField, startNode, Time, pos, FV_DC, FA_DC, ErrStat2, ErrMsg2, PosOffset=PosOffset)
        CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      ! Add contributions from IfW current field if node is in water
      DO i = 1, NumPoints
         FV(:,i) = FV(:,i) + nodeInWater(i) * FV_DC(:,i)
         FA(:,i) = FA(:,i) + nodeInWater(i) * FA_DC(:,i)
      END DO

   END IF
   
END SUBROUTINE WaveField_GetWaveKin

! This subroutine is intended for AeroDyn when modeling MHK turbines
SUBROUTINE WaveField_GetWaveVelAcc_AD( WaveField, SeaSt_Interp_m, StartNode, Time, pos, FV, FA, ErrStat, ErrMsg )
   TYPE(SeaSt_WaveFieldType), INTENT( IN    ) :: WaveField
   TYPE(SeaSt_Interp_MiscVarType), INTENT(INOUT) :: SeaSt_Interp_m
   INTEGER(IntKi),            INTENT( IN    ) :: StartNode
   REAL(DbKi),                INTENT( IN    ) :: Time
   REAL(ReKi),                INTENT( IN    ) :: pos(:,:) ! z=0 at MSL
   REAL(ReKi),                INTENT(   OUT ) :: FV(:,:)
   REAL(ReKi), ALLOCATABLE,   INTENT(   OUT ) :: FA(:,:)
   INTEGER(IntKi),            INTENT(   OUT ) :: ErrStat ! Error status of the operation
   CHARACTER(*),              INTENT(   OUT ) :: ErrMsg  ! Error message if errStat /= ErrID_None
   INTEGER(IntKi), ALLOCATABLE                :: nodeInWater(:)
   INTEGER(IntKi)                             :: NumPoints, i
   REAL(SiKi)                                 :: FV_node(3), FA_node(3)
   REAL(ReKi)                                 :: PosOffset(3), MSL2SWL, WtrDpth
   REAL(ReKi), ALLOCATABLE                    :: FV_DC(:,:), FA_DC(:,:)
   LOGICAL                                    :: getAcc

   CHARACTER(*),              PARAMETER       :: RoutineName = 'WaveField_GetWaveVelAcc_AD'
   INTEGER(IntKi)                             :: errStat2
   CHARACTER(ErrMsgLen)                       :: errMsg2

   ErrStat   = ErrID_None
   ErrMsg    = ""

   MSL2SWL   = WaveField%MSL2SWL
   WtrDpth   = WaveField%EffWtrDpth - MSL2SWL
   getAcc    = ALLOCATED(FA)
   NumPoints = size(pos, dim=2)

   ALLOCATE( nodeInWater(NumPoints), STAT=ErrStat2)
   IF (ErrStat2 /= 0) then
      CALL SetErrStat( ErrID_Info, 'Error allocating FA_DC', ErrStat, ErrMsg, RoutineName )
      RETURN
   END IF  

   ! Note: SeaState wavefield grid has z=0 on the SWL
   IF (getAcc) THEN
      DO i = 1, NumPoints
         CALL WaveField_GetNodeWaveVelAcc( WaveField, SeaSt_Interp_m, Time, pos(:,i)-(/0.0,0.0,MSL2SWL/), .FALSE., .FALSE., nodeInWater(i), FV_node, FA_node, ErrStat2, ErrMsg2 )
           CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         FV(:, i) = REAL(FV_node,   ReKi)
         FA(:, i) = REAL(FA_node,   ReKi)
      END DO
   ELSE
     DO i = 1, NumPoints
         CALL WaveField_GetNodeWaveVel( WaveField, SeaSt_Interp_m, Time, pos(:,i)-(/0.0,0.0,MSL2SWL/), .FALSE., .FALSE., nodeInWater(i), FV_node, ErrStat2, ErrMsg2 )
           CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         FV(:, i) = REAL(FV_node,   ReKi)
      END DO
   END IF

   ! If dynamic current field from IfW is present, get velocity and acceleration contributions
   IF (WaveField%hasCurrField) THEN
      PosOffset = (/0.0_ReKi,0.0_ReKi,WtrDpth/) ! IfW FlowField grid effectively has z=0 on the seabed
      ALLOCATE(FV_DC( 3, NumPoints ), STAT=ErrStat2)
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat( ErrID_Info, 'Error allocating FV_DC', ErrStat, ErrMsg, RoutineName )
         RETURN
      END IF
      IF (getAcc) THEN    
         ALLOCATE(FA_DC( 3, NumPoints ), STAT=ErrStat2)
         IF (ErrStat2 /= 0) THEN
            CALL SetErrStat( ErrID_Info, 'Error allocating FA_DC', ErrStat, ErrMsg, RoutineName )
            RETURN
         END IF
      END IF
      CALL IfW_FlowField_GetVelAcc(WaveField%CurrField, StartNode, Time, pos, FV_DC, FA_DC, ErrStat2, ErrMsg2, PosOffset=PosOffset)
        CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

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
   
END SUBROUTINE WaveField_GetWaveVelAcc_AD

END MODULE SeaSt_WaveField
