   SUBROUTINE UpdateConfiguration_AM2(uinc,uf,node_total,dof_node)

   REAL(ReKi),    INTENT(IN   ):: uinc(:)
   INTEGER(IntKi),INTENT(IN   ):: node_total
   INTEGER(IntKi),INTENT(IN   ):: dof_node
   REAL(ReKi),    INTENT(  OUT):: uf(:)

   INTEGER(IntKi)              :: i
   INTEGER(IntKi)              :: j
   INTEGER(IntKi)              :: temp_id

   DO i=1, node_total
       temp_id   = (i - 1) * dof_node
       DO j=1,6
           uf(temp_id+j) = uf(temp_id+j) + uinc(temp_id+j)
       ENDDO
   ENDDO

   END SUBROUTINE UpdateConfiguration_AM2
