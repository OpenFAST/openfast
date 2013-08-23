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