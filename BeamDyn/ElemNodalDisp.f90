 SUBROUTINE ElemNodalDisp(uu,node_elem,dof_node,nelem,norder,Nu)

   REAL(ReKi),INTENT(IN)::uu(:)
   INTEGER,INTENT(IN)::node_elem,dof_node,nelem,norder
   REAL(ReKi),INTENT(INOUT)::Nu(:)

   INTEGER::i,j,temp_id1,temp_id2

   DO i=1,node_elem
       DO j=1,dof_node
           temp_id1 = (i-1)*dof_node+j
           temp_id2 = ((nelem-1)*norder+i-1)*dof_node+j
           Nu(temp_id1) = uu(temp_id2)
       ENDDO
   ENDDO

   END SUBROUTINE ElemNodalDisp