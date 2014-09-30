   SUBROUTINE DynamicSolution(uuN0,uuN,vvN,Stif0,Mass0,gravity,u,time,&
                             &node_elem,dof_node,elem_total,dof_total,node_total,ngp,&
                             &qddot)
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
   REAL(ReKi),INTENT(OUT):: qddot(:) ! Second time derivative of state "q"

   ! Local variables
   
   REAL(ReKi):: MassM(dof_total,dof_total) 
   REAL(ReKi):: RHS(dof_total) 
   REAL(ReKi):: MassM_LU(dof_total-6,dof_total-6) 
   REAL(ReKi):: RHS_LU(dof_total-6) 
   REAL(ReKi):: F_PointLoad(dof_total) 
   REAL(ReKi):: qdd_temp(dof_total-6) 
   REAL(ReKi):: d 
   INTEGER(IntKi):: indx(dof_total-6) 
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
   RHS = RHS + F_PointLoad
   DO j=1,dof_total-6
       RHS_LU(j) = RHS(j+6)
       DO k=1,dof_total-6
           MassM_LU(j,k) = MassM(j+6,k+6)
       ENDDO
   ENDDO
   
   CALL ludcmp(MassM_LU,dof_total-6,indx,d)
   CALL lubksb(MassM_LU,dof_total-6,indx,RHS_LU,qdd_temp)

   qddot = 0.0D0
   DO j=1,dof_total-6
       qddot(j+6) = qdd_temp(j)
   ENDDO

!   DO j=1,18
!       WRITE(*,*) elm(i,j+1), elm(i,j+2),elm(i,j+3),elm(i,j+4),elm(i,j+5),elm(i,j+6)
!       WRITE(*,*) qddot(j)
!   ENDDO
!   STOP
   END SUBROUTINE DynamicSolution
