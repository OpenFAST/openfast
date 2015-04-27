   SUBROUTINE BD_AssembleStiffK(nelem,node_elem,dof_elem,dof_node,ElemK,GlobalK)
   !-------------------------------------------------------------------------------
   ! This subroutine assembles total stiffness matrix.
   !-------------------------------------------------------------------------------
   REAL(ReKi),INTENT(IN)::ElemK(:,:) ! Element mass matrix
   INTEGER(IntKi),INTENT(IN)::nelem ! Number of elements
   INTEGER(IntKi),INTENT(IN)::node_elem ! Nodes per element
   INTEGER(IntKi),INTENT(IN)::dof_elem ! Degrees of freedom per element
   INTEGER(IntKi),INTENT(IN)::dof_node ! Degrees of freedom per node
   REAL(ReKi),INTENT(INOUT)::GlobalK(:,:) ! Global stiffness matrix

   INTEGER(IntKi)::i,j,temp_id1,temp_id2

   DO i=1,dof_elem
       temp_id1 = (nelem-1)*(node_elem-1)*dof_node+i
       DO j=1,dof_elem
           temp_id2 = (nelem-1)*(node_elem-1)*dof_node+j
           GlobalK(temp_id1,temp_id2) = GlobalK(temp_id1,temp_id2) + ElemK(i,j)
       ENDDO
   ENDDO

   END SUBROUTINE BD_AssembleStiffK
