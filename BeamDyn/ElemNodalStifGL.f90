   SUBROUTINE ElemNodalStifGL(Stif0,node_elem,dof_node,nelem,NStif0)

   REAL(ReKi),INTENT(IN):: Stif0(:,:,:)
   INTEGER(IntKi),INTENT(IN):: node_elem,dof_node,nelem
   REAL(ReKi),INTENT(INOUT):: NStif0(:,:,:)

   INTEGER(IntKi):: i,j,k,temp_id

   DO i=1,node_elem
       DO j=1,dof_node
           DO k=1,dof_node
               temp_id = (nelem - 1)*(node_elem-1)+i
               NStif0(j,k,i) = Stif0(j,k,temp_id)
           ENDDO
       ENDDO
   ENDDO

   END SUBROUTINE ElemNodalStifGL

