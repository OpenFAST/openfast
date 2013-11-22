   SUBROUTINE BldGenerateStaticElement(uuN0,uuNf,Fext,Stif0,elem_total,node_elem,dof_node,ngp,StifK,RHS)

   REAL(ReKi),INTENT(IN):: uuN0(:),uuNf(:),Fext(:),Stif0(:,:)
   INTEGER(IntKi),INTENT(IN):: elem_total,node_elem,dof_node,ngp
   REAL(ReKi),INTENT(INOUT):: StifK(:,:),RHS(:) 

   REAL(ReKi),ALLOCATABLE:: Nuu0(:),Nuuu(:),Next(:),Nrr0(:),Nrrr(:)
   REAL(ReKi),ALLOCATABLE:: elk(:,:), elf(:)

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

   ALLOCATE(Next(dof_elem),STAT = allo_stat)
   IF(allo_stat /=0) GOTO 9999

   ALLOCATE(elf(dof_elem),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   elf = 0.0D0

   ALLOCATE(elk(dof_elem,dof_elem),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   elk = 0.0D0

   DO nelem=1,elem_total
       CALL ElemNodalDispGL(uuN0,node_elem,dof_node,nelem,Nuu0)
       CALL ElemNodalDispGL(uuNf,node_elem,dof_node,nelem,Nuuu)

       CALL ElemNodalDispGL(Fext,node_elem,dof_node,nelem,Next)

       CALL NodalRelRotGL(Nuu0,node_elem,dof_node,Nrr0)
       CALL NodalRelRotGL(Nuuu,node_elem,dof_node,Nrrr)
       
       elk = 0.0D0
       elf = 0.0D0
!       CALL ElementMatrixGL(Nuu0,Nuuu,Nrr0,Nrrr,Next,Stif0,ngp,node_elem,dof_node,elk,elf)
       CALL ElementMatrixLSGL(Nuu0,Nuuu,Nrr0,Nrrr,Next,Stif0,ngp,node_elem,dof_node,elk,elf)

       CALL AssembleStiffKGL(nelem,node_elem,dof_elem,dof_node,elk,StifK)
       CALL AssembleRHSGL(nelem,dof_elem,node_elem,dof_node,elf,RHS)
   ENDDO

   DEALLOCATE(Nuu0)
   DEALLOCATE(Nuuu)
   DEALLOCATE(Nrr0)
   DEALLOCATE(Nrrr)
   DEALLOCATE(Next)
   DEALLOCATE(elf)
   DEALLOCATE(elk)

   9999 IF(allo_stat/=0) THEN
            IF(ALLOCATED(Nuu0)) DEALLOCATE(Nuu0)
            IF(ALLOCATED(Nuuu)) DEALLOCATE(Nuuu)
            IF(ALLOCATED(Nrr0)) DEALLOCATE(Nrr0)
            IF(ALLOCATED(Nrrr)) DEALLOCATE(Nrrr)
            IF(ALLOCATED(Next)) DEALLOCATE(Next)
            IF(ALLOCATED(elf)) DEALLOCATE(elf)
            IF(ALLOCATED(elk)) DEALLOCATE(elk)
        ENDIF


   END SUBROUTINE BldGenerateStaticElement
