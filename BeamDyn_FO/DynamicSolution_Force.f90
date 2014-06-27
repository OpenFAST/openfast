   SUBROUTINE DynamicSolution_Force(uuN0,uuN,vvN,Stif0,Mass0,gravity,u,time,&
                                   &node_elem,dof_node,elem_total,dof_total,node_total,ngp,&
                                   &qddot,Force)
   !***************************************************************************************
   ! This subroutine calls other subroutines to apply the force, build the beam element 
   ! stiffness and mass matrices, build nodal force vector.  The output of this subroutine
   ! is the second time derivative of state "q".   
   !***************************************************************************************
   REAL(ReKi),INTENT(IN):: uuN0(:,:) ! Initial position vector
   REAL(ReKi),INTENT(IN):: Stif0(:,:,:) ! Element stiffness matrix
   REAL(ReKi),INTENT(IN):: Mass0(:,:,:) ! Element stiffness matrix
   REAL(ReKi),INTENT(IN):: gravity(:) ! 
   TYPE(BD_InputType),           INTENT(IN   )  :: u           ! Inputs at t
   REAL(ReKi),INTENT(IN):: uuN(:) ! Displacement of Mass 1: m
   REAL(ReKi),INTENT(IN):: vvN(:) ! Velocity of Mass 1: m/s
   REAL(DbKi),INTENT(IN):: time ! Current time
   INTEGER(IntKi),INTENT(IN):: node_elem ! Node per element
   INTEGER(IntKi),INTENT(IN):: dof_node ! Degrees of freedom per element
   INTEGER(IntKi),INTENT(IN):: elem_total ! Total number of elements
   INTEGER(IntKi),INTENT(IN):: dof_total ! Total number of degrees of freedom
   INTEGER(IntKi),INTENT(IN):: node_total ! Total number of nodes
   INTEGER(IntKi),INTENT(IN):: ngp ! Number of Gauss points
   REAL(ReKi),INTENT(IN):: qddot(:) ! Second time derivative of state "q"
   
   REAL(ReKi),INTENT(OUT):: Force(:)

   ! Local variables
   
   REAL(ReKi):: MassM(dof_total,dof_total) 
   REAL(ReKi):: RHS(dof_total) 
   REAL(ReKi):: F_PointLoad(dof_total)  
   REAL(ReKi):: d 
   INTEGER(IntKi):: j 
   INTEGER(IntKi):: k 
   INTEGER(IntKi):: temp_id


   RHS = 0.0D0
   MassM = 0.0D0

   CALL GenerateDynamicElement(uuN0,uuN,vvN,Stif0,Mass0,gravity,u,&
                              &elem_total,node_elem,dof_node,ngp,RHS,MassM)
   DO j=1,node_total
       temp_id = (j-1)*dof_node
       F_PointLoad(temp_id+1:temp_id+3) = u%PointLoad%Force(1:3,j)
       F_PointLoad(temp_id+4:temp_id+6) = u%PointLoad%Moment(1:3,j)
   ENDDO
   RHS(:) = RHS(:) + F_PointLoad(:)
   
   Force(:) = 0.0D0
   Force(:) = RHS(:) - MATMUL(MassM,qddot)
   
   END SUBROUTINE DynamicSolution_Force
