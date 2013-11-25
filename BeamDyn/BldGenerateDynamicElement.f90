   SUBROUTINE BldGenerateDynamicElement(uuN0,uuNf,vvNf,aaNf,Fext,Stif0,m00,mEta0,rho0,&
                                       &elem_total,node_elem,dof_node,ngp,StifK,RHS,MassM,DampG)

   REAL(ReKi),INTENT(IN):: uuN0(:),uuNf(:),Fext(:),Stif0(:,:)
   REAL(ReKi),INTENT(IN):: vvNf(:),aaNf(:),m00,mEta0(:),rho0(:,:)
   INTEGER(IntKi),INTENT(IN):: elem_total,node_elem,dof_node,ngp
   REAL(ReKi),INTENT(OUT):: StifK(:,:),RHS(:) 
   REAL(ReKi),INTENT(OUT):: MassM(:,:),DampG(:,:)

   REAL(ReKi),ALLOCATABLE:: Nuu0(:),Nuuu(:),Next(:),Nrr0(:),Nrrr(:)
   REAL(ReKi),ALLOCATABLE:: Nvvv(:),Naaa(:)
   REAL(ReKi),ALLOCATABLE:: elk(:,:),elf(:)
   REAL(ReKi),ALLOCATABLE:: elm(:,:),elg(:,:)

   INTEGER(IntKi):: dof_elem,rot_elem
   INTEGER(IntKi):: nelem,j
   INTEGER(IntKi):: allo_stat

   dof_elem = dof_node * node_elem
   rot_elem = (dof_node/2) * node_elem

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

   DO nelem=1,elem_total
       CALL ElemNodalDispGL(uuN0,node_elem,dof_node,nelem,Nuu0)
       CALL ElemNodalDispGL(uuNf,node_elem,dof_node,nelem,Nuuu)

       CALL ElemNodalDispGL(Fext,node_elem,dof_node,nelem,Next)

       CALL NodalRelRotGL(Nuu0,node_elem,dof_node,Nrr0)
       CALL NodalRelRotGL(Nuuu,node_elem,dof_node,Nrrr)
       
       CALL ElemNodalDispGL(vvNf,node_elem,dof_node,nelem,Nvvv)
       CALL ElemNodalDispGL(aaNf,node_elem,dof_node,nelem,Naaa)
       
!       CALL ElementMatrixDynGL(Nuu0,Nuuu,Nrr0,Nrrr,Next,Nvvv,Naaa,Stif0,m00,mEta0,rho0,&
!                              &ngp,node_elem,dof_node,elk,elf,elm,elg)

       CALL ElementMatrixDynLSGL(Nuu0,Nuuu,Nrr0,Nrrr,Next,Nvvv,Naaa,Stif0,m00,mEta0,rho0,&
                                &ngp,node_elem,dof_node,elk,elf,elm,elg)

       CALL AssembleStiffKGL(nelem,node_elem,dof_elem,dof_node,elk,StifK)
       CALL AssembleStiffKGL(nelem,node_elem,dof_elem,dof_node,elm,MassM)
       CALL AssembleStiffKGL(nelem,node_elem,dof_elem,dof_node,elg,DampG)
       CALL AssembleRHSGL(nelem,dof_elem,node_elem,dof_node,elf,RHS)
   ENDDO

   DEALLOCATE(Nuu0)
   DEALLOCATE(Nuuu)
   DEALLOCATE(Nrr0)
   DEALLOCATE(Nrrr)
   DEALLOCATE(Nvvv)
   DEALLOCATE(Naaa)
   DEALLOCATE(Next)
   DEALLOCATE(elf)
   DEALLOCATE(elk)
   DEALLOCATE(elm)
   DEALLOCATE(elg)

   9999 IF(allo_stat/=0) THEN
            IF(ALLOCATED(Nuu0)) DEALLOCATE(Nuu0)
            IF(ALLOCATED(Nuuu)) DEALLOCATE(Nuuu)
            IF(ALLOCATED(Nrr0)) DEALLOCATE(Nrr0)
            IF(ALLOCATED(Nrrr)) DEALLOCATE(Nrrr)
            IF(ALLOCATED(Nvvv)) DEALLOCATE(Nvvv)
            IF(ALLOCATED(Naaa)) DEALLOCATE(Naaa)
            IF(ALLOCATED(Next)) DEALLOCATE(Next)
            IF(ALLOCATED(elf)) DEALLOCATE(elf)
            IF(ALLOCATED(elk)) DEALLOCATE(elk)
            IF(ALLOCATED(elm)) DEALLOCATE(elm)
            IF(ALLOCATED(elg)) DEALLOCATE(elg)
        ENDIF


   END SUBROUTINE BldGenerateDynamicElement
