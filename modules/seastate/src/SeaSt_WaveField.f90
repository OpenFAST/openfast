MODULE SeaSt_WaveField

USE SeaState_Interp
USE SeaSt_WaveField_Types

IMPLICIT NONE
   
PRIVATE

! Public functions and subroutines
PUBLIC WaveField_GetWaveElev1
PUBLIC WaveField_GetWaveElev2
PUBLIC WaveField_GetTotalWaveElev
PUBLIC WaveField_GetWaveNormal
PUBLIC WaveField_GetWaveKin
PUBLIC WaveField_GetWaveVel

CONTAINS

!-------------------- Subroutine for wave elevation ------------------!
FUNCTION WaveField_GetWaveElev1( WaveField, Time, pos, ErrStat, ErrMsg )
   TYPE(SeaSt_WaveFieldType), INTENT( IN    ) :: WaveField
   REAL(DbKi),      INTENT( IN    ) :: Time
   REAL(ReKi),      INTENT( IN    ) :: pos(*)  ! Position at which free-surface elevation is to be calculated. Third entry ignored if present.
   INTEGER(IntKi),  INTENT(   OUT ) :: ErrStat ! Error status of the operation
   CHARACTER(*),    INTENT(   OUT ) :: ErrMsg  ! Error message if errStat /= ErrID_None
  
   REAL(SiKi)                       :: WaveField_GetWaveElev1
   REAL(SiKi)                       :: Zeta
   LOGICAL                          :: FirstWarn_Clamp
   CHARACTER(*),    PARAMETER       :: RoutineName = 'WaveField_GetWaveElev1'
   INTEGER(IntKi)                   :: errStat2
   CHARACTER(ErrMsgLen)             :: errMsg2
   
   ErrStat   = ErrID_None
   ErrMsg    = ""
   
   IF (ALLOCATED(WaveField%WaveElev1)) THEN
      Zeta = SeaSt_Interp_3D( Time, pos(1:2), WaveField%WaveElev1, WaveField%seast_interp_p, FirstWarn_Clamp, ErrStat2, ErrMsg2 )
        CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   ELSE
      Zeta = 0.0_SiKi
   END IF
   
   WaveField_GetWaveElev1 = Zeta

END FUNCTION WaveField_GetWaveElev1

FUNCTION WaveField_GetWaveElev2( WaveField, Time, pos, ErrStat, ErrMsg )
   TYPE(SeaSt_WaveFieldType), INTENT( IN    ) :: WaveField
   REAL(DbKi),      INTENT( IN    ) :: Time
   REAL(ReKi),      INTENT( IN    ) :: pos(*)  ! Position at which free-surface elevation is to be calculated. Third entry ignored if present.
   INTEGER(IntKi),  INTENT(   OUT ) :: ErrStat ! Error status of the operation
   CHARACTER(*),    INTENT(   OUT ) :: ErrMsg  ! Error message if errStat /= ErrID_None
   
   REAL(SiKi)                       :: WaveField_GetWaveElev2
   REAL(SiKi)                       :: Zeta
   LOGICAL                          :: FirstWarn_Clamp
   CHARACTER(*),    PARAMETER       :: RoutineName = 'WaveField_GetWaveElev2'
   INTEGER(IntKi)                   :: errStat2
   CHARACTER(ErrMsgLen)             :: errMsg2
   
   ErrStat   = ErrID_None
   ErrMsg    = ""
   
   IF (ALLOCATED(WaveField%WaveElev2)) THEN
      Zeta = SeaSt_Interp_3D( Time, pos(1:2), WaveField%WaveElev2, WaveField%seast_interp_p, FirstWarn_Clamp, ErrStat2, ErrMsg2 )
        CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   ELSE
      Zeta = 0.0_SiKi
   END IF

   WaveField_GetWaveElev2 = Zeta

END FUNCTION WaveField_GetWaveElev2

FUNCTION WaveField_GetTotalWaveElev( WaveField, Time, pos, ErrStat, ErrMsg )
   TYPE(SeaSt_WaveFieldType), INTENT( IN    ) :: WaveField
   REAL(DbKi),      INTENT( IN    ) :: Time
   REAL(ReKi),      INTENT( IN    ) :: pos(*)  ! Position at which free-surface elevation is to be calculated. Third entry ignored if present.
   INTEGER(IntKi),  INTENT(   OUT ) :: ErrStat ! Error status of the operation
   CHARACTER(*),    INTENT(   OUT ) :: ErrMsg  ! Error message if errStat /= ErrID_None

   REAL(SiKi)                       :: WaveField_GetTotalWaveElev
   REAL(SiKi)                       :: Zeta1, Zeta2
   LOGICAL                          :: FirstWarn_Clamp
   CHARACTER(*),    PARAMETER       :: RoutineName = 'WaveField_GetTotalWaveElev'
   INTEGER(IntKi)                   :: errStat2
   CHARACTER(ErrMsgLen)             :: errMsg2

   ErrStat   = ErrID_None
   ErrMsg    = ""
   
   Zeta1 = WaveField_GetWaveElev1( WaveField, Time, pos, ErrStat2, ErrMsg2 )
     CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   Zeta2 = WaveField_GetWaveElev2( WaveField, Time, pos, ErrStat2, ErrMsg2 )
     CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   WaveField_GetTotalWaveElev = Zeta1 + Zeta2
   
END FUNCTION WaveField_GetTotalWaveElev

SUBROUTINE WaveField_GetWaveNormal( WaveField, Time, pos, r, n, ErrStat, ErrMsg )

   TYPE(SeaSt_WaveFieldType), INTENT( IN    ) :: WaveField
   REAL(DbKi),                INTENT( IN    ) :: Time
   REAL(ReKi),                INTENT( IN    ) :: pos(*)  ! Position at which free-surface normal is to be calculated. Third entry ignored if present.
   REAL(ReKi),                INTENT( IN    ) :: r       ! Distance for central differencing
   REAL(ReKi),                INTENT(   OUT ) :: n(3)    ! Free-surface normal vector
   INTEGER(IntKi),            INTENT(   OUT ) :: ErrStat ! Error status of the operation
   CHARACTER(*),              INTENT(   OUT ) :: ErrMsg  ! Error message if errStat /= ErrID_None
   REAL(SiKi)                                 :: ZetaP,ZetaM
   REAL(ReKi)                                 :: r1,dZetadx,dZetady
   CHARACTER(*),              PARAMETER       :: RoutineName = 'WaveField_GetFreeSurfaceNormal'
   INTEGER(IntKi)                             :: errStat2
   CHARACTER(ErrMsgLen)                       :: errMsg2
   ErrStat   = ErrID_None
   ErrMsg    = ""

   r1 = MAX(r,1.0e-6) ! In case r is zero

   ZetaP = WaveField_GetTotalWaveElev( WaveField, Time, (/pos(1)+r1,pos(2)/), ErrStat2, ErrMsg2 )
     CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   ZetaM = WaveField_GetTotalWaveElev( WaveField, Time, (/pos(1)-r1,pos(2)/), ErrStat2, ErrMsg2 )
     CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   dZetadx = REAL(ZetaP-ZetaM,ReKi)/(2.0_ReKi*r1)
      
   ZetaP = WaveField_GetTotalWaveElev( WaveField, Time, (/pos(1),pos(2)+r1/), ErrStat2, ErrMsg2 )
     CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   ZetaM = WaveField_GetTotalWaveElev( WaveField, Time, (/pos(1),pos(2)-r1/), ErrStat2, ErrMsg2 )
     CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   dZetady = REAL(ZetaP-ZetaM,ReKi)/(2.0_ReKi*r1)
      
   n = (/-dZetadx,-dZetady,1.0_ReKi/)
   n = n / SQRT(Dot_Product(n,n))

END SUBROUTINE WaveField_GetWaveNormal

!-------------------- Subroutine for full wave field kinematics --------------------!
SUBROUTINE WaveField_GetWaveKin( WaveField, Time, pos, nodeInWater, WaveElev1, WaveElev2, WaveElev, FDynP, FV, FA, FAMCF, ErrStat, ErrMsg )
   TYPE(SeaSt_WaveFieldType), INTENT( IN    ) :: WaveField
   REAL(DbKi),                INTENT( IN    ) :: Time
   REAL(ReKi),                INTENT( IN    ) :: pos(3)
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

   REAL(ReKi)                                 :: posXY(2), posPrime(3), posXY0(3)
   TYPE(SeaSt_Interp_MiscVarType)             :: SeaSt_Interp_m
   LOGICAL                                    :: FirstWarn_Clamp
   CHARACTER(*),              PARAMETER       :: RoutineName = 'WaveField_GetWaveKin'
   INTEGER(IntKi)                             :: errStat2
   CHARACTER(ErrMsgLen)                       :: errMsg2

   ErrStat   = ErrID_None
   ErrMsg    = ""

   posXY    = pos(1:2)
   posXY0   = (/pos(1),pos(2),0.0_ReKi/)
   FAMCF(:) = 0.0

   ! Wave elevation
   WaveElev1 = WaveField_GetWaveElev1( WaveField, Time, pos, ErrStat2, ErrMsg2 )
     CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   WaveElev2 = WaveField_GetWaveElev2( WaveField, Time, pos, ErrStat2, ErrMsg2 )
     CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   WaveElev  = WaveElev1 + WaveElev2
    
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
         FDynP = SeaSt_Interp_4D    ( WaveField%WaveDynP, seast_interp_m, ErrStat2, ErrMsg2 )
           CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF ( ALLOCATED(WaveField%WaveAccMCF) ) THEN
            FAMCF(:) = SeaSt_Interp_4D_Vec( WaveField%WaveAccMCF, seast_interp_m, ErrStat2, ErrMsg2 )
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
      
      IF ( pos(3) <= WaveElev ) THEN ! Node is submerged
          
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
               FDynP = SeaSt_Interp_4D    ( WaveField%WaveDynP, seast_interp_m, ErrStat2, ErrMsg2 )
                 CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               IF ( ALLOCATED(WaveField%WaveAccMCF) ) THEN
                  FAMCF(:) = SeaSt_Interp_4D_Vec( WaveField%WaveAccMCF, seast_interp_m, ErrStat2, ErrMsg2 )
                    CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               END IF

            ELSE ! Node is above SWL - need wave stretching
          
               ! Vertical wave stretching
               CALL SeaSt_Interp_Setup( Time, posXY0, WaveField%seast_interp_p, seast_interp_m, ErrStat2, ErrMsg2 ) 
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
                  FV(:) = FV(:) + SeaSt_Interp_3D_vec( Time, posXY, WaveField%PWaveVel0,  WaveField%seast_interp_p, FirstWarn_Clamp, ErrStat2, ErrMsg2 ) * pos(3)
                    CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                  FA(:) = FA(:) + SeaSt_Interp_3D_vec( Time, posXY, WaveField%PWaveAcc0,  WaveField%seast_interp_p, FirstWarn_Clamp, ErrStat2, ErrMsg2 ) * pos(3)
                    CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                  FDynP = FDynP + SeaSt_Interp_3D    ( Time, posXY, WaveField%PWaveDynP0, WaveField%seast_interp_p, FirstWarn_Clamp, ErrStat2, ErrMsg2 ) * pos(3)
                    CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                  IF ( ALLOCATED(WaveField%WaveAccMCF) ) THEN
                     FAMCF(:) = FAMCF(:) + SeaSt_Interp_3D_vec( Time, posXY, WaveField%PWaveAccMCF0, WaveField%seast_interp_p, FirstWarn_Clamp, ErrStat2, ErrMsg2 ) * pos(3)
                       CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                  END IF
               END IF
          
            END IF ! Node is submerged
 
         ELSE ! Wheeler stretching - no need to check whether the node is above or below SWL
                  
            ! Map the node z-position linearly from [-EffWtrDpth,m%WaveElev(j)] to [-EffWtrDpth,0] 
            posPrime    = pos
            posPrime(3) = WaveField%EffWtrDpth*(WaveField%EffWtrDpth+pos(3))/(WaveField%EffWtrDpth+WaveElev)-WaveField%EffWtrDpth
                  
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
   
END SUBROUTINE WaveField_GetWaveKin

!-------------------- Subroutine for wave field velocity only --------------------!
SUBROUTINE WaveField_GetWaveVel( WaveField, Time, pos, nodeInWater, FV, ErrStat, ErrMsg )
   TYPE(SeaSt_WaveFieldType), INTENT( IN    ) :: WaveField
   REAL(DbKi),                INTENT( IN    ) :: Time
   REAL(ReKi),                INTENT( IN    ) :: pos(3)
   INTEGER(IntKi),            INTENT(   OUT ) :: nodeInWater
   REAL(SiKi),                INTENT(   OUT ) :: FV(3)
   INTEGER(IntKi),            INTENT(   OUT ) :: ErrStat ! Error status of the operation
   CHARACTER(*),              INTENT(   OUT ) :: ErrMsg  ! Error message if errStat /= ErrID_None

   REAL(SiKi)                                 :: WaveElev
   REAL(ReKi)                                 :: posXY(2), posPrime(3), posXY0(3)
   TYPE(SeaSt_Interp_MiscVarType)             :: SeaSt_Interp_m
   LOGICAL                                    :: FirstWarn_Clamp
   CHARACTER(*),              PARAMETER       :: RoutineName = 'WaveField_GetWaveVel'
   INTEGER(IntKi)                             :: errStat2
   CHARACTER(ErrMsgLen)                       :: errMsg2

   ErrStat   = ErrID_None
   ErrMsg    = ""

   posXY    = pos(1:2)
   posXY0   = (/pos(1),pos(2),0.0_ReKi/)

   ! Wave elevation
   WaveElev  = WaveField_GetTotalWaveElev( WaveField, Time, pos, ErrStat2, ErrMsg2 )
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
      
      IF ( pos(3) <= WaveElev ) THEN ! Node is submerged
          
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
                  FV(:) = FV(:) + SeaSt_Interp_3D_vec( Time, posXY, WaveField%PWaveVel0,  WaveField%seast_interp_p, FirstWarn_Clamp, ErrStat2, ErrMsg2 ) * pos(3)
                    CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               END IF
          
            END IF ! Node is submerged
 
         ELSE ! Wheeler stretching - no need to check whether the node is above or below SWL
                  
            ! Map the node z-position linearly from [-EffWtrDpth,m%WaveElev(j)] to [-EffWtrDpth,0] 
            posPrime    = pos
            posPrime(3) = WaveField%EffWtrDpth*(WaveField%EffWtrDpth+pos(3))/(WaveField%EffWtrDpth+WaveElev)-WaveField%EffWtrDpth
                  
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
   
END SUBROUTINE WaveField_GetWaveVel

END MODULE SeaSt_WaveField
