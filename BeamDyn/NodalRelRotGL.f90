   SUBROUTINE NodalRelRotGL(Nu,node_elem,dof_node,Nr)

   REAL(ReKi),INTENT(IN)::Nu(:)
   INTEGER(IntKi),INTENT(IN)::node_elem,dof_node

   REAL(ReKi),INTENT(INOUT)::Nr(:)

   INTEGER(IntKi)::i,k,temp_id
   REAL(ReKi)::Nu_temp1(3),Nu_temp(3),Nr_temp(3)

   Nr = 0.0D0
   Nu_temp1 = 0.0D0
   DO i=1,node_elem
       temp_id = (i - 1) * dof_node
       Nu_temp = 0.0D0
       DO k=1,3
           IF(i==1) Nu_temp1(k) = Nu(temp_id+k+3)
           Nu_temp(k) = Nu(temp_id+k+3)
       ENDDO
       Nr_temp = 0.0D0
       CALL CrvCompose(Nr_temp,Nu_temp1,Nu_temp,1)
       DO k=1,3
           temp_id = (i-1)*3+k
           Nr(temp_id) = Nr_temp(k)
       ENDDO
   ENDDO


   END SUBROUTINE NodalRelRotGL
