   SUBROUTINE CrvCompose_Check(pp,qq,flag,flag2)


   REAL(ReKi),INTENT(IN):: pp(:) 
   REAL(ReKi),INTENT(IN):: qq(:) 
   INTEGER(IntKi),INTENT(IN):: flag 
   INTEGER(IntKi),INTENT(OUT):: flag2 

   REAL(ReKi):: pp0,pp1,pp2,pp3,qq0,qq1,qq2,qq3,tr1,tr2,dd1,dd2

   flag2 = 0

   IF(flag==1 .OR. flag==3) THEN
       pp1 = -pp(1)
       pp2 = -pp(2)
       pp3 = -pp(3)
   ELSE
       pp1 = pp(1)
       pp2 = pp(2)
       pp3 = pp(3)
   ENDIF
   pp0 = 2.0D0 - (pp1 * pp1 + pp2 * pp2 + pp3 * pp3) / 8.0D0

   IF(flag==2 .OR. flag==3) THEN
       qq1 = -qq(1) 
       qq2 = -qq(2)
       qq3 = -qq(3)
   ELSE
       qq1 = qq(1)
       qq2 = qq(2)
       qq3 = qq(3)
   ENDIF
   qq0 = 2.0D0 - (qq1 * qq1 + qq2 * qq2 + qq3 * qq3)/8.0D0

   tr1 = (4.0D0 - pp0) * (4.0D0 - qq0)
   tr2 = pp0 * qq0 - pp1 * qq1 - pp2 * qq2 - pp3 * qq3
   dd1 = tr1 + tr2
   dd2 = tr1 - tr2

   IF(dd1>dd2) THEN
       tr1 = 4.0D0 / dd1
   ELSE
       flag2 = 1       
   ENDIF

   END SUBROUTINE CrvCompose_Check
