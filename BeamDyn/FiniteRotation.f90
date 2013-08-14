      MODULE FiniteRotatonLib
      
      USE GolbalDataFun
      IMPLICIT NONE
      
      PRIVATE

      PUBLIC CrvMatrixR,CrvCompose      

      CONTAINS

      SUBROUTINE CrvMatrixR(cc,Rr) 

      REAL(ReKi),INTENT(IN)::cc(:)
      REAL(ReKi),INTENT(OUT)::Rr(:,:)

      INTEGER:: i, j
      REAL(ReKi):: c1,c2,c3,c0,tr0

      Rr = ZERO
      
      c1 = cc(1)/FOUR
      c2 = cc(2)/FOUR
      c3 = cc(3)/FOUR

      c0 = HALF*(ONE-c1*c1-c2*c2-c3*c3)
      tr0 = ONE - c0
      tr0 = TWO/(tr0*tr0)

      Rr(1,1) = tr0*(c1*c1 + c0*c0) - ONE 
      Rr(2,1) = tr0*(c1*c2 + c0*c3)
      Rr(3,1) = tr0*(c1*c3 - c0*c2)
      Rr(1,2) = tr0*(c1*c2 - c0*c3)
      Rr(2,2) = tr0*(c2*c2 + c0*c0) - ONE 
      Rr(3,2) = tr0*(c2*c3 + c0*c1)
      Rr(1,3) = tr0*(c1*c3 + c0*c2)
      Rr(2,3) = tr0*(c2*c3 - c0*c1)
      Rr(3,3) = tr0*(c3*c3 + c0*c0) - ONE

      END SUBROUTINE CrvMatrixR

      SUBROUTINE CrvCompose(rr,pp,qq,flag)

      REAL(ReKi),INTENT(IN):: pp(:), qq(:)
      INTEGER,INTENT(IN):: flag
      REAL(ReKi),INTENT(OUT):: rr(:)

      REAL(ReKi):: pp0,pp1,pp2,pp3,qq0,qq1,qq2,qq3,tr1,tr2,dd1,dd2

      IF(flag==1 .OR. flag==3) THEN
          pp1 = -pp(1)
          pp2 = -pp(2)
          pp3 = -pp(3)
      ELSE
          pp1 = pp(1)
          pp2 = pp(2)
          pp3 = pp(3)
      ENDIF
      pp0 = TWO - (pp1*pp1 + pp2*pp2 + pp3*pp3)/EIGHT

      IF(flag==2 .OR. flag==3) THEN
          qq1 = -qq(1) 
          qq2 = -qq(2)
          qq3 = -qq(3)
      ELSE
          qq1 = qq(1)
          qq2 = qq(2)
          qq3 = qq(3)
      ENDIF
      qq0 = TWO - (qq1*qq1 + qq2*qq2 + qq3*qq3)/EIGHT

      tr1 = (FOUR - pp0) * (FOUR - qq0)
      tr2 = pp0*qq0 - pp1*qq1 - pp2*qq2 - pp3*qq3
      dd1 = tr1 + tr2
      dd2 = tr1 - tr2

      IF(dd1>dd2) THEN
          tr1 = FOUR/dd1
      ELSE
          tr1 = -FOUR/dd2
      ENDIF

      rr(1) = tr1 * (pp1*qq0 + pp0*qq1 - pp3*qq2 + pp2*qq3)
      rr(2) = tr1 * (pp2*qq0 + pp3*qq1 + pp0*qq2 - pp1*qq3)
      rr(3) = tr1 * (pp3*qq0 - pp2*qq1 + pp1*qq2 + pp0*qq3)

      END SUBROUTINE CrvCompose 
