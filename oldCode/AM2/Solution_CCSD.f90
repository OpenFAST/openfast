   SUBROUTINE Solution_CCSD(uuN0,uuN,vvN,Stif0,Mass0,gravity,u,&
                           &damp_flag,beta,&
                           &node_elem,dof_node,elem_total,dof_total,node_total,ngp,&
                           &xdot)
   !***************************************************************************************
   ! This subroutine calls other subroutines to apply the force, build the beam element 
   ! stiffness and mass matrices, build nodal force vector.  The output of this subroutine
   ! is the second time derivative of state "q".   
   !***************************************************************************************
   REAL(ReKi),                   INTENT(IN   ):: uuN0(:,:) ! Initial position vector
   REAL(ReKi),                   INTENT(IN   ):: Stif0(:,:,:) ! Element stiffness matrix
   REAL(ReKi),                   INTENT(IN   ):: Mass0(:,:,:) ! Element stiffness matrix
   REAL(ReKi),                   INTENT(IN   ):: gravity(:) ! 
   TYPE(BD_InputType),           INTENT(IN   ):: u           ! Inputs at t
   REAL(ReKi),                   INTENT(IN   ):: uuN(:) ! Displacement of Mass 1: m
   REAL(ReKi),                   INTENT(IN   ):: vvN(:) ! Velocity of Mass 1: m/s
   INTEGER(IntKi),               INTENT(IN   ):: damp_flag ! Total number of elements
   REAL(ReKi),                   INTENT(IN   ):: beta(:)
   INTEGER(IntKi),               INTENT(IN   ):: node_elem ! Node per element
   INTEGER(IntKi),               INTENT(IN   ):: dof_node ! Degrees of freedom per element
   INTEGER(IntKi),               INTENT(IN   ):: elem_total ! Total number of elements
   INTEGER(IntKi),               INTENT(IN   ):: dof_total ! Total number of degrees of freedom
   INTEGER(IntKi),               INTENT(IN   ):: node_total ! Total number of nodes
   INTEGER(IntKi),               INTENT(IN   ):: ngp ! Number of Gauss points
   TYPE(BD_ContinuousStateType), INTENT(INOUT):: xdot

   ! Local variables
   
   REAL(ReKi):: MassM(dof_total*2,dof_total*2) 
   REAL(ReKi):: MassM_LU(dof_total*2-12,dof_total*2-12) 
   REAL(ReKi):: RHS(dof_total*2) 
   REAL(ReKi):: RHS_LU(dof_total*2-12) 
   REAL(ReKi):: F_PointLoad(dof_total) 
   REAL(ReKi):: sol_temp(dof_total*2-12) 
   REAL(ReKi):: d 
   INTEGER(IntKi):: indx(dof_total*2-12) 
   INTEGER(IntKi):: j 
   INTEGER(IntKi):: k 
   INTEGER(IntKi):: temp_id

   RHS(:)     = 0.0D0
   MassM(:,:) = 0.0D0

   CALL GenerateDynamicElement_CCSD(uuN0,uuN,vvN,Stif0,Mass0,gravity,u,&
                                   &damp_flag,beta,&
                                   &elem_total,node_elem,dof_total,dof_node,ngp,RHS,MassM)
   DO j=1,node_total
       temp_id = (j-1)*dof_node
       F_PointLoad(temp_id+1:temp_id+3) = u%PointLoad%Force(1:3,j)
       F_PointLoad(temp_id+4:temp_id+6) = u%PointLoad%Moment(1:3,j)
   ENDDO
!WRITE(*,*) F_PointLoad
   RHS(dof_total+1:dof_total*2) = RHS(dof_total+1:dof_total*2) + F_PointLoad(1:dof_total) 
   DO j=1,dof_total-6
       RHS_LU(j) = RHS(j+6)
       RHS_LU(j+dof_total-6) = RHS(j+6+dof_total)
       DO k=1,dof_total-6
           MassM_LU(j,k) = MassM(j+6,k+6)
           MassM_LU(j,k+dof_total-6) = MassM(j+6,k+dof_total+6)
           MassM_LU(j+dof_total-6,k) = MassM(j+dof_total+6,k+6)
           MassM_LU(j+dof_total-6,k+dof_total-6) = MassM(j+dof_total+6,k+dof_total+6)
       ENDDO
   ENDDO

!DO j=1,60
!WRITE(*,*) j,MassM_LU(j,j)
!ENDDO
   sol_temp(:) = 0.0D0
   CALL ludcmp(MassM_LU,dof_total*2-12,indx,d)
   CALL lubksb(MassM_LU,dof_total*2-12,indx,RHS_LU,sol_temp)

!   xdot%q(:) = 0.0D0
!   xdot%dqdt(:) = 0.0D0
   DO j=1,dof_total-6
       xdot%q(j+6)    = sol_temp(j)
       xdot%dqdt(j+6) = sol_temp(dof_total-6+j)
   ENDDO

   END SUBROUTINE Solution_CCSD
