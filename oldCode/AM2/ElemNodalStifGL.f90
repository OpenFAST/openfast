!**********************************************************************************************************************************
!  
!   This subroutine allocates the stiffness matrix for each element. 
!
!**********************************************************************************************************************************  
   
   SUBROUTINE ElemNodalStifGL(Stif0,node_elem,dof_node,nelem,NStif0)

   REAL(ReKi),INTENT(IN):: Stif0(:,:,:) ! Sectional material properties at each node
   INTEGER(IntKi),INTENT(IN):: node_elem ! Nodes per element
   INTEGER(IntKi),INTENT(IN):: dof_node ! Degrees of freedom per node
   INTEGER(IntKi),INTENT(IN):: nelem ! Number of elements in beam
   REAL(ReKi),INTENT(INOUT):: NStif0(:,:,:) ! Output Sectional material properties for each element

   INTEGER(IntKi):: i,j,k,temp_id

   DO i=1,node_elem
       DO j=1,dof_node
           DO k=1,dof_node
               temp_id = (nelem - 1)*(node_elem-1)+i
               NStif0(j,k,i) = Stif0(j,k,temp_id)
           ENDDO
       ENDDO
   ENDDO

   END SUBROUTINE ElemNodalStifGL

