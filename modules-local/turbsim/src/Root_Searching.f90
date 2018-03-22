MODULE NEWTON_MOD
IMPLICIT NONE
PUBLIC :: F, FP
CONTAINS
!==============================================================================
FUNCTION F(X,URef,RefHt,HubHt) RESULT (Y)
USE NWTC_Library
REAL(ReKi), INTENT(IN) :: X
!REAL(ReKi), INTENT(IN) :: V_10
REAL(ReKi), INTENT(IN) :: URef
REAL(ReKi), INTENT(IN) :: RefHt
REAL(ReKi), INTENT(IN) :: HubHt
REAL(ReKi)             :: Y
REAL(ReKi)             :: C_1
REAL(ReKi)             :: C_2
REAL(ReKi)             :: D_1
REAL(ReKi)             :: D_2
REAL(ReKi)             :: D_3
REAL(ReKi)             :: D_4

IF ( EqualRealNos( RefHt, 10.0_ReKi ) ) THEN
   C_1=1.0-0.41*0.06*LOG(600.0/3600.0)
   C_2=0.41*0.06*0.043*LOG(600.0/3600.0)
   Y=C_1*X-C_2*X*X-URef
ELSEIF ( EqualRealNos( RefHt, HubHt) )  THEN
   D_1=1.0-0.41*0.06*LOG(600.0/3600.0)*(RefHt/10.0)**(-0.22)
   D_2=0.41*0.06*0.043*LOG(600.0/3600.0)*(RefHt/10.0)**(-0.22)
   D_3=D_1*0.0573*LOG(RefHt/10.0)
   D_4=D_2*0.0573*LOG(RefHt/10.0)
   Y=D_1*X-D_2*X*X+D_3*X*(1.0+0.15*X)**0.5-D_4*X*X*(1.0+0.15*X)**0.5-URef    
ELSE
   Y = -9999999   
ENDIF

END FUNCTION F
!==============================================================================
FUNCTION FP(X,URef,RefHt,HubHt) RESULT (Y)
USE NWTC_Library
REAL(ReKi), INTENT(IN) :: X
REAL(ReKi), INTENT(IN) :: URef
REAL(ReKi), INTENT(IN) :: RefHt
REAL(ReKi), INTENT(IN) :: HubHt
REAL(ReKi)             :: Y
REAL(ReKi)             :: C_1
REAL(ReKi)             :: C_2
REAL(ReKi)             :: D_1
REAL(ReKi)             :: D_2
REAL(ReKi)             :: D_3
REAL(ReKi)             :: D_4

IF ( EqualRealNos( RefHt, 10.0_ReKi ) ) THEN
   C_1 = 1.0-0.41*0.06*LOG(600.0/3600.0)
   C_2 = 0.41*0.06*0.043*LOG(600.0/3600.0)
   Y   = C_1-2.0*C_2*X
ELSEIF ( EqualRealNos( RefHt, HubHt) )  THEN
   C_1 = 1.0 - 0.41  * 0.06 * LOG(600.0/3600.0)
   C_2 = 0.41 * 0.06*0.043  * LOG(600.0/3600.0)
   D_1 = 1.0-0.41*0.06*LOG(600.0/3600.0)*(RefHt/10.0)**(-0.22)
   D_2 = 0.41*0.06*0.043*LOG(600.0/3600.0)*(RefHt/10.0)**(-0.22)
   D_3 = D_1*0.0573*LOG(RefHt/10.0)
   D_4 = D_2*0.0573*LOG(RefHt/10.0)
   Y   = D_1-2.0*D_2*X+D_3*((1.0+0.15*X)**0.5+0.15*X/2.0/(1.0+0.15*X)**0.5)  &
         -D_4*(2.0*X*(1.0+0.15*X)**0.5+0.15*X*X/2.0/(1.0+0.15*X)**0.5)
ELSE
   Y = -9999999
ENDIF

END FUNCTION FP
!==============================================================================
END MODULE NEWTON_MOD
    
    
SUBROUTINE ROOT_SEARCHING(X0,X,URef,RefHt,HubHt)
USE NEWTON_MOD
USE NWTC_Library
IMPLICIT NONE

REAL(ReKi), PARAMETER     :: TOL=1.0E-5

REAL(ReKi), INTENT(IN   ) :: URef
REAL(ReKi), INTENT(IN   ) :: RefHt
REAL(ReKi), INTENT(IN   ) :: HubHt
REAL(ReKi), INTENT(IN   ) :: X0
REAL(ReKi), INTENT(  OUT) :: X

X=X0
DO  ! bjj: I'd like this better if there were absolutely no way to have an infinite loop here...
   X = X - F(X,URef,RefHt,HubHt) / FP(X,URef,RefHt,HubHt)
   IF (ABS(F(X,URef,RefHt,HubHt)) < TOL) THEN
      EXIT
   END IF
   
END DO
END SUBROUTINE ROOT_SEARCHING
    
    

    