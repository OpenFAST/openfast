   SUBROUTINE BDyn_CalcContStateDeriv(uuN0,uuN,vvN,Stif0,m00,mEta0,rho0,&
                                     &node_elem,dof_node,elem_total,dof_total,node_total,ngp,&
                                     &time,qdot,qddot)

   REAL(ReKi),INTENT(IN):: uuN0(:),Stif0(:,:,:),m00(:),mEta0(:,:),rho0(:,:,:)
   REAL(ReKi),INTENT(IN):: uuN(:),vvN(:)
   REAL(ReKi),INTENT(IN):: time
   INTEGER(IntKi),INTENT(IN)::node_elem,dof_node,elem_total,dof_total,node_total,ngp

   REAL(ReKi),INTENT(OUT):: qdot(:),qddot(:)

   REAL(ReKi),ALLOCATABLE:: udN(:)
   INTEGER(IntKi):: i,allo_stat

   ALLOCATE(udN(dof_total), STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   udN = 0.0D0

   CALL ComputeUDN(node_total,dof_node,vvN,uuN,udN)

   DO i=1,dof_total
       qdot(i) = udN(i)
   ENDDO

   CALL DynamicSolution(uuN0,uuN,vvN,Stif0,m00,mEta0,rho0,time,&
                       &node_elem,dof_node,elem_total,dof_total,node_total,ngp,&
                       &qddot)
    
   DEALLOCATE(udN)

   9999 IF(allo_stat/=0) THEN
            IF(ALLOCATED(udN)) DEALLOCATE(udN)
        ENDIF


   END SUBROUTINE BDyn_CalcContStateDeriv
