!*********************************************************************************************************************************
!   This subroutine is called 3 times by GenerateDynamicElement. It creates an array "Nu" which are nodal values for each element
!   for uuN0, uuN, and vvn (initial nodal configuration, nodal displacements, and velocity of mass) these values 
!   are then passed back to GenerateDynamicElement.  
!**********************************************************************************************************************************
   SUBROUTINE BD_ElemNodalDisp(uu,node_elem,dof_node,nelem,Nu)

   REAL(ReKi),INTENT(IN):: uu(:) ! Initial position vector, nodal displacements, and velocity of mass
   INTEGER(IntKi),INTENT(IN):: node_elem ! Node per element
   INTEGER(IntKi),INTENT(IN):: dof_node ! Degrees of freedom per node
   INTEGER(IntKi),INTENT(IN):: nelem ! Number of elements, being looped from 1, elem_total
   REAL(ReKi),INTENT(INOUT):: Nu(:) ! Initial nodal configuration, nodal displacements, and velocity of mass at the nodes of each element

   INTEGER(IntKi):: i ! Index counter
   INTEGER(IntKi):: j ! Index counter
   INTEGER(IntKi):: temp_id1 ! Counter
   INTEGER(IntKi):: temp_id2 ! Counter

   DO i=1,node_elem
       DO j=1,dof_node
           temp_id1 = (i-1)*dof_node+j
           temp_id2 = ((nelem - 1)*(node_elem-1)+i-1)*dof_node+j
           Nu(temp_id1) = uu(temp_id2)
       ENDDO
   ENDDO

   END SUBROUTINE BD_ElemNodalDisp
