   SUBROUTINE DynamicSolution(uuN0,uuN,vvN,Stif0,m00,mEta0,rho0,time,&
                             &node_elem,dof_node,elem_total,dof_total,node_total,ngp,&
                             &qddot)

   REAL(ReKi),INTENT(IN):: uuN0(:),Stif0(:,:,:),m00(:),mEta0(:,:),rho0(:,:,:)
   REAL(ReKi),INTENT(IN):: uuN(:),vvN(:)
   REAL(DbKi),INTENT(IN):: time
   INTEGER(IntKi),INTENT(IN)::node_elem,dof_node,elem_total,dof_total,node_total,ngp
   REAL(ReKi),INTENT(OUT):: qddot(:)

   REAL(ReKi):: MassM(dof_total,dof_total),RHS(dof_total)
   REAL(ReKi):: MassM_LU(dof_total-6,dof_total-6),RHS_LU(dof_total-6)
   REAL(ReKi):: F_ext(dof_total),qdd_temp(dof_total-6)
   
   REAL(ReKi):: d
   INTEGER(IntKi):: indx(dof_total-6),j,k

   CALL AppliedNodalLoad(F_ext,time,dof_total)

   RHS = 0.0D0
   MassM = 0.0D0

   CALL GenerateDynamicElement(uuN0,uuN,vvN,Stif0,m00,mEta0,rho0,&
                              &elem_total,node_elem,dof_node,ngp,RHS,MassM)

   RHS = RHS + F_ext

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

   END SUBROUTINE DynamicSolution
