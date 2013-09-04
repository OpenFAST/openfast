   SUBROUTINE NodalRelRot(Nuuu,node_elem,nelem,norder,dof_node,Nrrr)

   REAL(ReKi),INTENT(IN)::Nuuu(:)
   INTEGER(IntKi),INTENT(IN)::node_elem,nelem,norder,dof_node

   REAL(ReKi),INTENT(INOUT)::Nrrr(:)

   INTEGER(IntKi)::i,k,temp_id
   REAL(ReKi)::Nuuu_temp1(3),Nuuu_temp(3),Nrrr_temp(3)
      
   DO i=1,node_elem
       temp_id = (i - 1) * dof_node
       DO k=1,3
           IF(i==1) Nuuu_temp1(k) = Nuuu(temp_id+k+3)
           Nuuu_temp(k) = Nuuu(temp_id+k+3)
       ENDDO
       CALL CrvCompose(Nrrr_temp,Nuuu_temp1,Nuuu_temp,1)
       DO k=1,3
           temp_id = (i-1)*3+k
           Nrrr(temp_id) = Nrrr_temp(k)
       ENDDO
   ENDDO

   END SUBROUTINE NodalRelRot
