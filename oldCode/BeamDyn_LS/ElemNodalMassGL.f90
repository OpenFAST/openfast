!**********************************************************************************************************************************
!  
!   This subroutine allocates mass, sectional m\eta, and Sectional tensor of inertia per unit span for each element. 
!
!********************************************************************************************************************************** 
   
   SUBROUTINE ElemNodalMassGL(m00,mEta0,rho0,node_elem,dof_node,nelem,Nm00,NmEta0,Nrho0)

   REAL(ReKi),INTENT(IN):: m00(:) ! Mass of beam per unit span at each node
   REAL(ReKi),INTENT(IN):: mEta0(:,:) ! Sectional m\Eta_0 at each node
   REAL(ReKi),INTENT(IN):: rho0(:,:,:) ! Sectional tensor of inertia per unit span
   INTEGER(IntKi),INTENT(IN):: node_elem ! Nodes per element
   INTEGER(IntKi),INTENT(IN):: dof_node ! Degrees of freedom per node
   INTEGER(IntKi),INTENT(IN):: nelem ! Number of elements in beam
   REAL(ReKi),INTENT(INOUT):: Nm00(:) ! Nodal mass of beam per unit span for each element 
   REAL(ReKi),INTENT(INOUT):: NmEta0(:,:) ! Nodal sectional m\Eta_0 for each element
   REAL(ReKi),INTENT(INOUT):: Nrho0(:,:,:) ! Nodal sectional tensor of inertia per unit span for each element

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


