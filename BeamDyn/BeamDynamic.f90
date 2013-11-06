   SUBROUTINE BeamDynamic(uuN0,uuNf,uuNv,uuNa,hhp,w,Jacobian,Stif0,m00,mEta0,rho0,&
                        &node_elem,dof_node,norder,elem_total,dof_total,node_total,dof_elem,&
                        &StifK,MassM,DampG,RHS)

   REAL(ReKi),INTENT(IN)::uuN0(:),uuNf(:),hhp(:,:),w(:),Jacobian
   REAL(ReKi),INTENT(IN)::uuNv(:),uuNa(:)
   REAL(ReKi),INTENT(IN)::Stif0(:,:,:)
   REAL(ReKi),INTENT(IN)::m00,mEta0(:),rho0(:,:) 
   INTEGER(IntKi),INTENT(IN)::node_elem,dof_node,norder,elem_total,dof_total,node_total
   INTEGER(IntKi),INTENT(IN)::dof_elem
   REAL(ReKi),INTENT(INOUT)::StifK(:,:),RHS(:)
   REAL(ReKi),INTENT(INOUT)::MassM(:,:),DampG(:,:)

   REAL(ReKi),ALLOCATABLE::Nuu0(:),Nuuu(:),Nrr0(:),Nrrr(:),Next(:)
   REAL(ReKi),ALLOCATABLE::Nvvv(:),Naaa(:)
   REAL(ReKi),ALLOCATABLE::elk(:,:),elf(:)
   REAL(ReKi),ALLOCATABLE::elm(:,:),elg(:,:)


   INTEGER(IntKi)::i,allo_stat,j
   INTEGER(IntKi)::rot_elem

   rot_elem = dof_elem / 2
   

   ALLOCATE(Nuu0(dof_elem),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   Nuu0 = 0.0D0
   
   ALLOCATE(Nuuu(dof_elem),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   Nuuu = 0.0D0

   ALLOCATE(Nrr0(rot_elem),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   Nrr0 = 0.0D0

   ALLOCATE(Nrrr(rot_elem),STAT = allo_stat) 
   IF(allo_stat/=0) GOTO 9999
   Nrrr = 0.0D0

   ALLOCATE(Nvvv(dof_elem),STAT = allo_stat) 
   IF(allo_stat/=0) GOTO 9999
   Nvvv = 0.0D0

   ALLOCATE(Naaa(dof_elem),STAT = allo_stat) 
   IF(allo_stat/=0) GOTO 9999
   Naaa = 0.0D0

   ALLOCATE(Next(dof_elem),STAT = allo_stat)
   IF(allo_stat /=0) GOTO 9999
   Next = 0.0D0
   
   ALLOCATE(elf(dof_elem),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   elf = 0.0D0

   ALLOCATE(elk(dof_elem,dof_elem),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   elk = 0.0D0

   ALLOCATE(elm(dof_elem,dof_elem),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   elm = 0.0D0

   ALLOCATE(elg(dof_elem,dof_elem),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   elg = 0.0D0

   DO i=1,elem_total
       ! Get Nodal Displacement Vector Nuu0 from Global Vector uuN0 at t=0    
       CALL ElemNodalDisp(uuN0,node_elem,dof_node,i,norder,Nuu0)
       !Get Nodal Displacement Vector Nuuu from Global Vector uuNf
       CALL ElemNodalDisp(uuNf,node_elem,dof_node,i,norder,Nuuu)

       CALL ElemNodalDisp(uuNv,node_elem,dof_node,i,norder,Nvvv)

       CALL ElemNodalDisp(uuNa,node_elem,dof_node,i,norder,Naaa)

!       CALL ElemNodalDisp(Fext,node_elem,dof_node,i,norder,Next)
       !Compute Nodal Relative Rotation Vector Nrr0 
       CALL NodalRelRot(Nuu0,node_elem,i,norder,dof_node,Nrr0)
       !Compute Nodal Relative Rotation Vector Nrrr
       CALL NodalRelRot(Nuuu,node_elem,i,norder,dof_node,Nrrr)
       !Compute Elemental Matrices: elf, elk
       elf = 0.0D0
       elk = 0.0D0    
       elm = 0.0D0    
       elg = 0.0D0    
!       DO j=1,dof_elem
!           WRITE(*,*) "j = ",j
!           WRITE(*,*) "Naaa = ",Naaa(j)
!       ENDDO       
       CALL ElementMatrixDyn(Nuu0,Nuuu,Nrr0,Nrrr,Next,Nvvv,Naaa,hhp,Stif0,m00,mEta0,rho0,Jacobian,&
                          &w,node_elem,i,norder,dof_node,elk,elf,elm,elg)
       !Assemble Elemental Matrices into Global Matrices
       CALL AssembleStiffK(i,dof_elem,norder,dof_node,elk,StifK)
       CALL AssembleStiffK(i,dof_elem,norder,dof_node,elm,MassM)
       CALL AssembleStiffK(i,dof_elem,norder,dof_node,elg,DampG)
!       DO j=1,dof_elem
!           WRITE(*,*) "j=",j,elf(j)
!       ENDDO
       CALL AssembleRHS(i,dof_elem,norder,dof_node,elf,RHS)
!      WRITE(*,*) "StiffK"
!      WRITE(*,*) StifK
!      STOP
   ENDDO

   DEALLOCATE(Nuu0)
   DEALLOCATE(Nuuu)
   DEALLOCATE(Nrr0)
   DEALLOCATE(Nrrr)
   DEALLOCATE(Next)
   DEALLOCATE(Nvvv)
   DEALLOCATE(Naaa)
   DEALLOCATE(elf)
   DEALLOCATE(elk)
   DEALLOCATE(elm)
   DEALLOCATE(elg)

   9999 IF(allo_stat/=0) THEN
!            IF(ALLOCATED(StifK)) DEALLOCATE(StifK)
!            IF(ALLOCATED(RHS)) DEALLOCATE(RHS)
            IF(ALLOCATED(Nuu0)) DEALLOCATE(Nuu0)
            IF(ALLOCATED(Nuuu)) DEALLOCATE(Nuuu)
            IF(ALLOCATED(Nrr0)) DEALLOCATE(Nrr0)
            IF(ALLOCATED(Nrrr)) DEALLOCATE(Nrrr)
            IF(ALLOCATED(Next)) DEALLOCATE(Next)
            IF(ALLOCATED(elf)) DEALLOCATE(elf)
            IF(ALLOCATED(elk)) DEALLOCATE(elk)
            IF(ALLOCATED(elm)) DEALLOCATE(elm)
            IF(ALLOCATED(elg)) DEALLOCATE(elg)
        ENDIF

   END SUBROUTINE BeamDynamic

