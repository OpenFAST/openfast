SUBROUTINE BeamGenerateStaticMatrix(uuN0,uuNf,hhp,w,Jacobian,Stif0,&
                                      &node_elem,dof_node,norder,elem_total,dof_total,node_total,&
                                      &StifK,RHS)

   REAL(ReKi),INTENT(IN)::uuN0(:),uuNf(:),hhp(:,:),w(:),Jacobian
   REAL(ReKi),INTENT(IN)::Stif0(:,:,:)
   INTEGER,INTENT(IN)::node_elem,dof_node,norder,elem_total,dof_total,node_total
   REAL(ReKi),INTENT(OUT),ALLOCATABLE::StifK(:,:),RHS(:)

   REAL(ReKi),ALLOCATABLE::Nuu0(:),Nuuu(:),Nrr0(:),Nrrr(:)
   REAL(ReKi),ALLOCATABLE::elk(:,:),elf(:)
   REAL(ReKi)::Stif(6,6)
   INTEGER::i
   INTEGER::dof_elem,rot_elem

   dof_elem = dof_node * node_elem
   rot_elem = 3 * node_elem
   
   ALLOCATE(StifK(dof_total,dof_total),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   StifK = ZERO
   
   ALLOCATE(RHS(dof_total),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   RHS = ZERO
   

   ALLOCATE(Nuu0(dof_elem),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   Nuu0 = ZERO
   
   ALLOCATE(Nuuu(dof_elem),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   Nuuu = ZERO

   ALLOCATE(Nrr0(rot_elem),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   Nrr0 = ZERO

   ALLOCATE(Nrrr(rot_elem),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   Nrrr = ZERO
   
   ALLOCATE(elf(dof_elem),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   elf = ZERO

   ALLOCATE(elk(dof_elem,dof_elem),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   elk = ZERO
 
   DO i=1,elem_total
       ! Get Nodal Displacement Vector Nuu0 from Global Vector uuN0 at t=0    
       CALL ElemNodalDisp(uuN0,node_elem,dof_node,i,norder,Nuu0)
       !Get Nodal Displacement Vector Nuuu from Global Vector uuNf
       CALL ElemNodalDisp(uuNf,node_elem,dof_node,i,norder,Nuuu)
       !Compute Nodal Relative Rotation Vector Nrr0 
       CALL NodalRelRot(Nuu0,node_elem,i,norder,dof_node,Nrr0)
       !Compute Nodal Relative Rotation Vector Nrrr
       CALL NodalRelRot(Nuuu,node_elem,i,norder,dof_node,Nrrr)
       !Compute Elemental Matrices: elf, elk
       elf = ZERO
       elk = ZERO    
       CALL ElementMatrix(Nuu0,Nuuu,Nrr0,Nrrr,hhp,Stif0,Jacobian,&
                          &w,node_elem,i,norder,dof_node,elk,elf)
       !Assemble Elemental Matrices into Global Matrices
       CALL AssembleStiffK(i,dof_elem,norder,dof_node,elk,StifK)
       CALL AssembleRHS(i,dof_elem,norder,dof_node,elf,RHS)
   ENDDO

   9999 IF(allo_stat/=0) THEN
            IF(ALLOCATED(StifK)) DEALLOCATE(StifK)
            IF(ALLOCATED(RHS)) DEALLOCATE(RHS)
            IF(ALLOCATED(Nuu0)) DEALLOCATE(Nuu0)
            IF(ALLOCATED(Nuuu)) DEALLOCATE(Nuuu)
            IF(ALLOCATED(Nrr0)) DEALLOCATE(Nrr0)
            IF(ALLOCATED(Nrrr)) DEALLOCATE(Nrrr)
            IF(ALLOCATED(elf)) DEALLOCATE(elf)
            IF(ALLOCATED(elk)) DEALLOCATE(elk)
        ENDIF

   END SUBROUTINE BeamGenerateStaticMatrix