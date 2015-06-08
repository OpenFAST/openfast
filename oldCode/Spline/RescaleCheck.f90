   SUBROUTINE RescaleCheck(x,node_total,counter)

   TYPE(BD_ContinuousStateType),      INTENT(INOUT)  :: x
   INTEGER(IntKi),                    INTENT(IN   )  :: node_total
   INTEGER(IntKi),                    INTENT(INOUT)  :: counter

   REAL(ReKi)                                        :: temp_pp(3)
   REAL(ReKi)                                        :: temp_qq(3)
   REAL(ReKi)                                        :: temp_rr(3)
   REAL(ReKi)                                        :: temp_norm(node_total)
   REAL(ReKi)                                        :: temp_min
   INTEGER(IntKi)                                    :: i
   INTEGER(IntKi)                                    :: j
   INTEGER(IntKi)                                    :: k
   INTEGER(IntKi)                                    :: temp_id

   temp_norm(:) = 0.0D0
   DO i=1,node_total
       temp_id = (i-1)*6
       temp_rr(1:3) = x%q(temp_id+4:temp_id+6)
       temp_norm(i) = Norm(temp_rr)
   ENDDO
!WRITE(*,*) temp_norm
   temp_min = MINVAL(temp_norm)
!   temp_min = temp_norm(1)
   IF(temp_min .GE. 4.0D0) THEN
       DO i=1,node_total
           temp_id = (i-1)*6
           temp_pp(:) = 0.0D0
           temp_qq(:) = 0.0D0
           DO j=1,3
               temp_pp(j) = x%q(temp_id+3+j)
           ENDDO
           CALL CrvCompose_temp2(temp_rr,temp_pp,temp_qq,0)
           DO j=1,3
               x%q(temp_id+3+j) = temp_rr(j)
           ENDDO
       ENDDO
       counter = counter + 1
   ENDIF

   END SUBROUTINE RescaleCheck
