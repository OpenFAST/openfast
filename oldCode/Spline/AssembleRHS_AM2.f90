   SUBROUTINE AssembleRHS_AM2(nelem,dof_elem,node_elem,dof_node,ElemRHS1,ElemRHS2,GlobalRHS)
   !-------------------------------------------------------------------------------
   ! This subroutine assembles global force vector.
   !-------------------------------------------------------------------------------

   REAL(ReKi),    INTENT(IN   ):: ElemRHS1(:) ! Total element force (Fc, Fd, Fb)
   REAL(ReKi),    INTENT(IN   ):: ElemRHS2(:) ! Total element force (Fc, Fd, Fb)
   INTEGER(IntKi),INTENT(IN   ):: nelem       ! Number of elements
   INTEGER(IntKi),INTENT(IN   ):: dof_elem ! Degrees of freedom per element
   INTEGER(IntKi),INTENT(IN   ):: node_elem ! Nodes per element
   INTEGER(IntKi),INTENT(IN   ):: dof_node ! Degrees of freedom per node
   REAL(ReKi),    INTENT(  OUT):: GlobalRHS(:) ! Global force 

   INTEGER(IntKi)              ::i
   INTEGER(IntKi)              ::temp_id

   DO i=1,dof_elem
       temp_id = (nelem-1)*(node_elem-1)*dof_node+i
       GlobalRHS(temp_id) = GlobalRHS(temp_id)+ElemRHS1(i)
   ENDDO

   DO i=1,dof_elem
       temp_id = (nelem-1)*(node_elem-1)*dof_node+i+dof_elem
       GlobalRHS(temp_id) = GlobalRHS(temp_id)+ElemRHS2(i)
   ENDDO

   END SUBROUTINE AssembleRHS_AM2
