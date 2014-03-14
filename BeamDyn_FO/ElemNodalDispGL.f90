   SUBROUTINE ElemNodalDispGL(uu,node_elem,dof_node,nelem,Nu)

   REAL(ReKi),INTENT(IN):: uu(:)
   INTEGER(IntKi),INTENT(IN):: node_elem,dof_node,nelem
   REAL(ReKi),INTENT(INOUT):: Nu(:)

   INTEGER(IntKi):: i,j,temp_id1,temp_id2

   DO i=1,node_elem
       DO j=1,dof_node
           temp_id1 = (i-1)*dof_node+j
           temp_id2 = ((nelem - 1)*(node_elem-1)+i-1)*dof_node+j
           Nu(temp_id1) = uu(temp_id2)
       ENDDO
   ENDDO

   END SUBROUTINE ElemNodalDispGL
