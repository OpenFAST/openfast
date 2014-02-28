   SUBROUTINE ElemNodalMassGL(m00,mEta0,rho0,node_elem,dof_node,nelem,Nm00,NmEta0,Nrho0)

   REAL(ReKi),INTENT(IN):: m00(:),mEta0(:,:),rho0(:,:,:)
   INTEGER(IntKi),INTENT(IN):: node_elem,dof_node,nelem
   REAL(ReKi),INTENT(INOUT):: Nm00(:),NmEta0(:,:),Nrho0(:,:,:)

   INTEGER(IntKi):: i,j,k,temp_id

   DO i=1,node_elem
       temp_id = (nelem - 1)*(node_elem-1)+i
       Nm00(i) = m00(temp_id)
       DO j=1,dof_node/2
           NmEta0(j,i) = mEta0(j,temp_id)
           DO k=1,dof_node/2
               Nrho0(j,k,i) = rho0(j,k,temp_id)
           ENDDO
       ENDDO
   ENDDO

   END SUBROUTINE ElemNodalMassGL


