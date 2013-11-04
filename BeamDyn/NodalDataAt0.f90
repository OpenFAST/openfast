   SUBROUTINE NodalDataAt0(node_elem,nelem,norder,dof_node,nnode,hhp,Nuu0,Jac,E10)

   INTEGER(IntKi),INTENT(IN)::node_elem,nelem,norder,dof_node,nnode
   REAL(ReKi),INTENT(IN)::hhp(:,:),Nuu0(:),Jac
   REAL(ReKi),INTENT(INOUT)::E10(:)

   INTEGER(IntKi)::i,temp_id


   DO i=1,node_elem
       temp_id = (i-1)*dof_node  
       E10(1) = E10(1) + hhp(i,nnode)*Nuu0(temp_id+1)
       E10(2) = E10(2) + hhp(i,nnode)*Nuu0(temp_id+2)
       E10(3) = E10(3) + hhp(i,nnode)*Nuu0(temp_id+3)
   ENDDO
   E10 = E10 / Jac


   END SUBROUTINE NodalDataAt0
