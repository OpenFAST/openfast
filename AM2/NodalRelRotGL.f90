!**********************************************************************************************************************************
!   This subroutine is called 2 times by GenerateDynamicElement. It uses the element arrays "Nu" which represents Nuu0, 
!   and Nuuu (Nodal initial position for each element, and Nodal displacement of Mass 1 for each element respectively) 
!   created by ElemNodalDispGL and sends the values to CrvCompose to calculate the rotation parameters. 
!   The rotation parameters are returned to BldGenerateStaticElement and then sent to ElementMatrixLSGL.
!**********************************************************************************************************************************  
  SUBROUTINE NodalRelRotGL(Nu,node_elem,dof_node,Nr)

   REAL(ReKi),INTENT(IN):: Nu(:) ! Nodal initial position for each element, and Nodal displacement of Mass 1 for each element respectively
   INTEGER(IntKi),INTENT(IN):: node_elem ! Node per element
   INTEGER(IntKi),INTENT(IN):: dof_node ! Degrees of freedom per node
   REAL(ReKi),INTENT(INOUT):: Nr(:) ! Output, Nodal rotation parameters.

   INTEGER(IntKi)::i ! Index counter
   INTEGER(IntKi)::k ! Index counter
   INTEGER(IntKi)::temp_id ! Counter to get 4,5,6 DOF from Nu
   REAL(ReKi)::Nu_temp1(3) ! 4th, 5th, 6th DOF of Nu for each node, sent to CrvCompose
   REAL(ReKi)::Nu_temp(3) ! 4th, 5th, 6th DOF of Nu for each node, sent to CrvCompose
   REAL(ReKi)::Nr_temp(3) ! Rotation parameter returned from CrvCompose

   Nr = 0.0D0
   Nu_temp1 = 0.0D0
   DO i=1,node_elem
       temp_id = (i - 1) * dof_node
       Nu_temp = 0.0D0
       DO k=1,3
           IF(i==1) Nu_temp1(k) = Nu(temp_id+k+3)
           Nu_temp(k) = Nu(temp_id+k+3)
       ENDDO
       Nr_temp = 0.0D0
!       CALL CrvCompose_temp(Nr_temp,Nu_temp1,Nu_temp,1)
       CALL CrvCompose(Nr_temp,Nu_temp1,Nu_temp,1)
       DO k=1,3
           temp_id = (i-1)*3+k
           Nr(temp_id) = Nr_temp(k)
       ENDDO
   ENDDO


   END SUBROUTINE NodalRelRotGL
