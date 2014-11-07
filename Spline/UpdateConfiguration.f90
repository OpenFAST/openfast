   SUBROUTINE UpdateConfiguration(uinc,uf,node_total,dof_node)

   REAL(ReKi),    INTENT(IN   ):: uinc(:)
   INTEGER(IntKi),INTENT(IN   ):: node_total
   INTEGER(IntKi),INTENT(IN   ):: dof_node
   REAL(ReKi),    INTENT(INOUT):: uf(:)

   REAL(ReKi)                  :: rotf_temp(3)
   REAL(ReKi)                  :: roti_temp(3)
   REAL(ReKi)                  :: rot_temp(3)
   INTEGER(IntKi)              :: i
   INTEGER(IntKi)              :: j
   INTEGER(IntKi)              :: temp_id

   DO i=1, node_total
       temp_id   = (i - 1) * dof_node
       rotf_temp = 0.0D0
       roti_temp = 0.0D0
       rot_temp  = 0.0D0
       DO j=1,3
           uf(temp_id+j) = uf(temp_id+j) + uinc(temp_id+j)
           uf(temp_id+j+3) = uf(temp_id+j+3) + uinc(temp_id+j+3)
           rotf_temp(j)  = uf(temp_id+3+j)
           roti_temp(j)  = uinc(temp_id+3+j)
       ENDDO
       CALL CrvCompose_temp(rot_temp,roti_temp,rotf_temp,0)
       DO j=1,3
!           uf(temp_id+3+j) = rot_temp(j)
       ENDDO
   ENDDO

   END SUBROUTINE UpdateConfiguration
