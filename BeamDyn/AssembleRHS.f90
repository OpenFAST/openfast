      SUBROUTINE AssembleRHS(nelem,dof_elem,norder,dof_node,ElemRHS,GlobalRHS)

      REAL(ReKi),INTENT(IN)::ElemRHS(:)
      INTEGER,INTENT(IN)::nelem,dof_elem,norder,dof_node
      REAL(ReKi),INTENT(INOUT)::GlobalRHS(:)

      INTEGER::i,temp_id
      
      DO i=1,dof_elem
          temp_id = (nelem-1)*norder*dof_node+i
          GlobalRHS(temp_id) = GlobalRHS(temp_id)+ElemRHS(i)
      ENDDO      

      END SUBROUTINE AssembleRHS