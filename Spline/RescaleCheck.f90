   SUBROUTINE RescaleCheck(x,node_total,flag_scale)

   TYPE(BD_ContinuousStateType),      INTENT(INOUT)  :: x
   INTEGER(IntKi),                    INTENT(IN   )  :: node_total
   INTEGER(IntKi),                    INTENT(IN   )  :: flag_scale

   REAL(ReKi)                                        :: temp_pp(3)
   REAL(ReKi)                                        :: temp_qq(3)
   REAL(ReKi)                                        :: temp_rr(3)
   INTEGER(IntKi)                                    :: i
   INTEGER(IntKi)                                    :: j
   INTEGER(IntKi)                                    :: k
   INTEGER(IntKi)                                    :: temp_id
   INTEGER(IntKi)                                    :: flag

   WRITE(*,*) "flag_scale",flag_scale

   IF(flag_scale .EQ. 1) THEN
       DO i=2,node_total
           temp_id = (i-1)*6
           temp_pp(:) = 0.0D0
           temp_qq(:) = 0.0D0
           DO k=1,3
               temp_pp(k) = x%q(temp_id+3+k)
           ENDDO
           CALL CrvCompose_temp2(temp_rr,temp_pp,temp_qq,0)
WRITE(*,*) temp_rr
           DO k=1,3
               x%q(temp_id+3+k) = temp_rr(k)
           ENDDO
       ENDDO
   ELSE
       DO i=2,node_total
           temp_id = (i-1)*6
           temp_pp(:) = 0.0D0
           temp_qq(:) = 0.0D0
           DO j=1,3
               temp_pp(j) = x%q(temp_id+3+j)
           ENDDO
           CALL CrvCompose_Check(temp_pp,temp_qq,0,flag)
           IF(flag .EQ. 1) THEN
               WRITE(*,*) "Rescaled"
               DO j=1,node_total
                   temp_id = (j-1)*6
                   temp_pp(:) = 0.0D0
                   temp_qq(:) = 0.0D0
                   DO k=1,3
                       temp_pp(k) = x%q(temp_id+3+k)
                   ENDDO
                   CALL CrvCompose_temp2(temp_rr,temp_pp,temp_qq,0)
WRITE(*,*) temp_rr
                   DO k=1,3
                       x%q(temp_id+3+k) = temp_rr(k)
                   ENDDO
               ENDDO
               EXIT
           ENDIF
       ENDDO
   ENDIF

   END SUBROUTINE RescaleCheck
