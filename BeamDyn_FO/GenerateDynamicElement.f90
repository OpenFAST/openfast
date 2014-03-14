   SUBROUTINE GenerateDynamicElement(uuN0,uuN,vvN,Stif0,m00,mEta0,rho0,&
                                    &elem_total,node_elem,dof_node,ngp,RHS,MassM)

   REAL(ReKi),INTENT(IN):: uuN0(:),uuN(:),vvN(:),Stif0(:,:,:)
   REAL(ReKi),INTENT(IN):: m00(:),mEta0(:,:),rho0(:,:,:)
   INTEGER(IntKi),INTENT(IN):: elem_total,node_elem,dof_node,ngp
   REAL(ReKi),INTENT(OUT):: MassM(:,:),RHS(:)   

   REAL(ReKi),ALLOCATABLE:: Nuu0(:),Nuuu(:),Nrr0(:),Nrrr(:)
   REAL(ReKi),ALLOCATABLE:: Nvvv(:)
   REAL(ReKi),ALLOCATABLE:: NStif0(:,:,:)
   REAL(ReKi),ALLOCATABLE:: Nm00(:),NmEta0(:,:),Nrho0(:,:,:)
   REAL(ReKi),ALLOCATABLE:: elf(:),elm(:,:)

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

   ALLOCATE(NStif0(dof_node,dof_node,node_elem),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   NStif0 = 0.0D0

   ALLOCATE(Nm00(node_elem),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   Nm00 = 0.0D0
   
   ALLOCATE(NmEta0(3,node_elem),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   NmEta0 = 0.0D0
   
   ALLOCATE(Nrho0(3,3,node_elem),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   Nrho0 = 0.0D0

   ALLOCATE(elf(dof_elem),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   elf = 0.0D0

   ALLOCATE(elm(dof_elem,dof_elem),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   elm = 0.0D0

   DO nelem=1,elem_total
       CALL ElemNodalDispGL(uuN0,node_elem,dof_node,nelem,Nuu0)
       CALL ElemNodalDispGL(uuN,node_elem,dof_node,nelem,Nuuu)
       CALL ElemNodalStifGL(Stif0,node_elem,dof_node,nelem,NStif0) 
       CALL ElemNodalMassGL(m00,mEta0,rho0,node_elem,dof_node,nelem,Nm00,NmEta0,Nrho0)
       
       CALL NodalRelRotGL(Nuu0,node_elem,dof_node,Nrr0)
       CALL NodalRelRotGL(Nuuu,node_elem,dof_node,Nrrr)

       CALL ElemNodalDispGL(vvN,node_elem,dof_node,nelem,Nvvv)

       CALL ElementMatrix(Nuu0,Nuuu,Nrr0,Nrrr,Nvvv,NStif0,Nm00,NmEta0,Nrho0,&
                         &ngp,node_elem,dof_node,elf,elm)

       CALL AssembleStiffKGL(nelem,node_elem,dof_elem,dof_node,elm,MassM)
       CALL AssembleRHSGL(nelem,dof_elem,node_elem,dof_node,elf,RHS)

   ENDDO

   DEALLOCATE(Nuu0)
   DEALLOCATE(Nuuu)
   DEALLOCATE(Nrr0)
   DEALLOCATE(Nrrr)
   DEALLOCATE(Nvvv)
   DEALLOCATE(NStif0)
   DEALLOCATE(Nm00)
   DEALLOCATE(NmEta0)
   DEALLOCATE(Nrho0)
   DEALLOCATE(elf)
   DEALLOCATE(elm)

   9999 IF(allo_stat/=0) THEN
            IF(ALLOCATED(Nuu0)) DEALLOCATE(Nuu0)
            IF(ALLOCATED(Nuuu)) DEALLOCATE(Nuuu)
            IF(ALLOCATED(Nrr0)) DEALLOCATE(Nrr0)
            IF(ALLOCATED(Nrrr)) DEALLOCATE(Nrrr)
            IF(ALLOCATED(Nvvv)) DEALLOCATE(Nvvv)
            IF(ALLOCATED(NStif0)) DEALLOCATE(NStif0)
            IF(ALLOCATED(Nm00)) DEALLOCATE(Nm00)
            IF(ALLOCATED(NmEta0)) DEALLOCATE(NmEta0)
            IF(ALLOCATED(Nrho0)) DEALLOCATE(Nrho0)
            IF(ALLOCATED(elf)) DEALLOCATE(elf)
            IF(ALLOCATED(elm)) DEALLOCATE(elm)
        ENDIF


   END SUBROUTINE GenerateDynamicElement
