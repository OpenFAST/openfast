   SUBROUTINE AssembleRHSGL(nelem,dof_elem,node_elem,dof_node,ElemRHS,GlobalRHS)

   REAL(ReKi),INTENT(IN)::ElemRHS(:)
   INTEGER(IntKi),INTENT(IN)::nelem,dof_elem,node_elem,dof_node
   REAL(ReKi),INTENT(INOUT)::GlobalRHS(:)

   INTEGER(IntKi)::i,temp_id

   DO i=1,dof_elem
       temp_id = (nelem-1)*(node_elem-1)*dof_node+i
       GlobalRHS(temp_id) = GlobalRHS(temp_id)+ElemRHS(i)
   ENDDO 

   END SUBROUTINE AssembleRHSGL
