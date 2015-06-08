   SUBROUTINE AssembleStiffK_AM2(nelem,node_elem,dof_elem,dof_node,&
                                &ElemK11,ElemK12,ElemK21,ElemK22,GlobalK)
   !-------------------------------------------------------------------------------
   ! This subroutine assembles total stiffness matrix.
   !-------------------------------------------------------------------------------
   REAL(ReKi),    INTENT(IN   ):: ElemK11(:,:) ! Element mass matrix
   REAL(ReKi),    INTENT(IN   ):: ElemK12(:,:) ! Element mass matrix
   REAL(ReKi),    INTENT(IN   ):: ElemK21(:,:) ! Element mass matrix
   REAL(ReKi),    INTENT(IN   ):: ElemK22(:,:) ! Element mass matrix
   INTEGER(IntKi),INTENT(IN   ):: nelem ! Number of elements
   INTEGER(IntKi),INTENT(IN   ):: node_elem ! Nodes per element
   INTEGER(IntKi),INTENT(IN   ):: dof_elem ! Degrees of freedom per element
   INTEGER(IntKi),INTENT(IN   ):: dof_node ! Degrees of freedom per node
   REAL(ReKi),    INTENT(  OUT):: GlobalK(:,:) ! Global stiffness matrix

   INTEGER(IntKi)              :: i
   INTEGER(IntKi)              :: j
   INTEGER(IntKi)              :: temp_id1
   INTEGER(IntKi)              :: temp_id2

   DO i=1,dof_elem
       temp_id1 = (nelem-1)*(node_elem-1)*dof_node+i
       DO j=1,dof_elem
           temp_id2 = (nelem-1)*(node_elem-1)*dof_node+j
           GlobalK(temp_id1,temp_id2) = GlobalK(temp_id1,temp_id2) + ElemK11(i,j)
       ENDDO
   ENDDO

   DO i=1,dof_elem
       temp_id1 = (nelem-1)*(node_elem-1)*dof_node+i
       DO j=1,dof_elem
           temp_id2 = (nelem-1)*(node_elem-1)*dof_node+j+dof_elem
           GlobalK(temp_id1,temp_id2) = GlobalK(temp_id1,temp_id2) + ElemK12(i,j)
       ENDDO
   ENDDO

   DO i=1,dof_elem
       temp_id1 = (nelem-1)*(node_elem-1)*dof_node+i+dof_elem
       DO j=1,dof_elem
           temp_id2 = (nelem-1)*(node_elem-1)*dof_node+j
           GlobalK(temp_id1,temp_id2) = GlobalK(temp_id1,temp_id2) + ElemK21(i,j)
       ENDDO
   ENDDO

   DO i=1,dof_elem
       temp_id1 = (nelem-1)*(node_elem-1)*dof_node+i+dof_elem
       DO j=1,dof_elem
           temp_id2 = (nelem-1)*(node_elem-1)*dof_node+j+dof_elem
           GlobalK(temp_id1,temp_id2) = GlobalK(temp_id1,temp_id2) + ElemK22(i,j)
       ENDDO
   ENDDO

   END SUBROUTINE AssembleStiffK_AM2


