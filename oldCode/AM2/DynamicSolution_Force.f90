   SUBROUTINE DynamicSolution_Force(uuN0,uuN,vvN,aaN,                                      &
                                    Stif0,Mass0,gravity,u,                                 &
                                    damp_flag,beta,                                        &
                                    node_elem,dof_node,elem_total,dof_total,node_total,ngp,&
                                    analysis_type,Force,ReactionForce)
   !***************************************************************************************
   ! This subroutine calls other subroutines to apply the force, build the beam element 
   ! stiffness and mass matrices, build nodal force vector.  The output of this subroutine
   ! is the second time derivative of state "q".   
   !***************************************************************************************
   REAL(ReKi),         INTENT(IN   ):: uuN0(:,:) ! Initial position vector
   REAL(ReKi),         INTENT(IN   ):: uuN(:) ! Displacement of Mass 1: m
   REAL(ReKi),         INTENT(IN   ):: vvN(:) ! Velocity of Mass 1: m/s
   REAL(ReKi),         INTENT(IN   ):: aaN(:) ! Velocity of Mass 1: m/s
   REAL(ReKi),         INTENT(IN   ):: Stif0(:,:,:) ! Element stiffness matrix
   REAL(ReKi),         INTENT(IN   ):: Mass0(:,:,:) ! Element stiffness matrix
   REAL(ReKi),         INTENT(IN   ):: gravity(:) ! 
   INTEGER(IntKi),     INTENT(IN   ):: damp_flag ! Number of Gauss points
   REAL(ReKi),         INTENT(IN   ):: beta(:)
   TYPE(BD_InputType), INTENT(IN   ):: u           ! Inputs at t
   INTEGER(IntKi),     INTENT(IN   ):: node_elem ! Node per element
   INTEGER(IntKi),     INTENT(IN   ):: dof_node ! Degrees of freedom per element
   INTEGER(IntKi),     INTENT(IN   ):: elem_total ! Total number of elements
   INTEGER(IntKi),     INTENT(IN   ):: dof_total ! Total number of degrees of freedom
   INTEGER(IntKi),     INTENT(IN   ):: node_total ! Total number of nodes
   INTEGER(IntKi),     INTENT(IN   ):: ngp ! Number of Gauss points
   INTEGER(IntKi),     INTENT(IN   ):: analysis_type ! Number of Gauss points
   REAL(ReKi),         INTENT(  OUT):: Force(:)
   REAL(ReKi),         INTENT(  OUT):: ReactionForce(:)

   ! Local variables
   
   REAL(ReKi):: RHS(dof_total) 
   REAL(ReKi):: Reaction(6) 
   REAL(ReKi):: d 
   INTEGER(IntKi):: j 
   INTEGER(IntKi):: k 
   INTEGER(IntKi):: temp_id


   RHS(:) = 0.0D0
   Reaction(:) = 0.0D0

   CALL GenerateDynamicElement_Force(uuN0,uuN,vvN,aaN,     &
                                     Stif0,Mass0,gravity,u,&
                                     damp_flag,beta,&
                                     elem_total,node_elem,dof_node,ngp,RHS,Reaction)
   
   Force(:) = 0.0D0
   Force(:) = RHS(:)
   ReactionForce(:) = 0.0D0
   ReactionForce(:) = Reaction(:)
   
   END SUBROUTINE DynamicSolution_Force
