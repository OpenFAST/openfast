   SUBROUTINE AssembleRHSGL(nelem,dof_elem,node_elem,dof_node,ElemRHS,GlobalRHS)
   !-------------------------------------------------------------------------------
   ! This subroutine assembles global force vector.
   !-------------------------------------------------------------------------------

   REAL(ReKi),INTENT(IN)::ElemRHS(:) ! Total element force (Fc, Fd, Fb)
   INTEGER(IntKi),INTENT(IN)::nelem ! Number of elements
   INTEGER(IntKi),INTENT(IN)::dof_elem ! Degrees of freedom per element
   INTEGER(IntKi),INTENT(IN)::node_elem ! Nodes per element
   INTEGER(IntKi),INTENT(IN)::dof_node ! Degrees of freedom per node
   REAL(ReKi),INTENT(INOUT)::GlobalRHS(:) ! Global force 

   INTEGER(IntKi)::i,temp_id

   DO i=1,dof_elem
       temp_id = (nelem-1)*(node_elem-1)*dof_node+i
       GlobalRHS(temp_id) = GlobalRHS(temp_id)+ElemRHS(i)
   ENDDO 

   END SUBROUTINE AssembleRHSGL
