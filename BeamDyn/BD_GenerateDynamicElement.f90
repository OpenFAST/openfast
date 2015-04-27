   SUBROUTINE BD_GenerateDynamicElement(uuN0,uuNf,vvNf,aaNf,                 &
                                        Stif0,Mass0,gravity,u,damp_flag,beta,&
                                        elem_total,node_elem,dof_node,ngp,   &
                                        StifK,RHS,MassM,DampG)

   REAL(ReKi),        INTENT(IN   ):: uuN0(:,:)
   REAL(ReKi),        INTENT(IN   ):: uuNf(:)
   REAL(ReKi),        INTENT(IN   ):: vvNf(:)
   REAL(ReKi),        INTENT(IN   ):: aaNf(:)
   REAL(ReKi),        INTENT(IN   ):: Stif0(:,:,:)
   REAL(ReKi),        INTENT(IN   ):: Mass0(:,:,:)
   REAL(ReKi),        INTENT(IN   ):: gravity(:)
   TYPE(BD_InputType),INTENT(IN   ):: u
   INTEGER(IntKi),    INTENT(IN   ):: damp_flag
   REAL(ReKi),        INTENT(IN   ):: beta(:)
   INTEGER(IntKi),    INTENT(IN   ):: elem_total
   INTEGER(IntKi),    INTENT(IN   ):: node_elem
   INTEGER(IntKi),    INTENT(IN   ):: dof_node
   INTEGER(IntKi),    INTENT(IN   ):: ngp

   REAL(ReKi),        INTENT(  OUT):: StifK(:,:)
   REAL(ReKi),        INTENT(  OUT):: RHS(:) 
   REAL(ReKi),        INTENT(  OUT):: MassM(:,:)
   REAL(ReKi),        INTENT(  OUT):: DampG(:,:)

   REAL(ReKi),          ALLOCATABLE:: Nuu0(:)
   REAL(ReKi),          ALLOCATABLE:: Nuuu(:)
   REAL(ReKi),          ALLOCATABLE:: Nrr0(:)
   REAL(ReKi),          ALLOCATABLE:: Nrrr(:)
   REAL(ReKi),          ALLOCATABLE:: Nvvv(:)
   REAL(ReKi),          ALLOCATABLE:: Naaa(:)
   REAL(ReKi),          ALLOCATABLE:: elk(:,:)
   REAL(ReKi),          ALLOCATABLE:: elf(:)
   REAL(ReKi),          ALLOCATABLE:: elm(:,:)
   REAL(ReKi),          ALLOCATABLE:: elg(:,:)
   REAL(ReKi)                      :: EStif0_GL(6,6,node_elem-1)
   REAL(ReKi)                      :: EMass0_GL(6,6,node_elem-1)
   REAL(ReKi)                      :: DistrLoad_GL(6,node_elem-1)
   INTEGER(IntKi)                  :: dof_elem
   INTEGER(IntKi)                  :: rot_elem
   INTEGER(IntKi)                  :: nelem
   INTEGER(IntKi)                  :: j
   INTEGER(IntKi)                  :: temp_id
   INTEGER(IntKi)                  :: allo_stat

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
   
   Nrrr = 0.0D0

   ALLOCATE(Nvvv(dof_elem),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   Nvvv = 0.0D0
   
   ALLOCATE(Naaa(dof_elem),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   Naaa = 0.0D0

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
       Nuu0(:) = uuN0(:,nelem)
       CALL BD_ElemNodalDisp(uuNf,node_elem,dof_node,nelem,Nuuu)
       temp_id = (nelem-1)*ngp
       DO j=1,ngp
           EStif0_GL(1:6,1:6,j) = Stif0(1:6,1:6,temp_id+j)
           EMass0_GL(1:6,1:6,j) = Mass0(1:6,1:6,temp_id+j)
           DistrLoad_GL(1:3,j)  = u%DistrLoad%Force(1:3,temp_id+j+1)
           DistrLoad_GL(4:6,j)  = u%DistrLoad%Moment(1:3,temp_id+j+1)
       ENDDO
       CALL BD_NodalRelRot(Nuu0,node_elem,dof_node,Nrr0)
       CALL BD_NodalRelRot(Nuuu,node_elem,dof_node,Nrrr)
       CALL BD_ElemNodalDisp(vvNf,node_elem,dof_node,nelem,Nvvv)
       CALL BD_ElemNodalDisp(aaNf,node_elem,dof_node,nelem,Naaa)

       CALL BD_ElementMatrixGA2(Nuu0,Nuuu,Nrr0,Nrrr,Nvvv,Naaa,           &
                              EStif0_GL,EMass0_GL,gravity,DistrLoad_GL,&
                              damp_flag,beta,                          &
                              ngp,node_elem,dof_node,elk,elf,elm,elg)

       CALL BD_AssembleStiffK(nelem,node_elem,dof_elem,dof_node,elk,StifK)
       CALL BD_AssembleStiffK(nelem,node_elem,dof_elem,dof_node,elm,MassM)
       CALL BD_AssembleStiffK(nelem,node_elem,dof_elem,dof_node,elg,DampG)
       CALL BD_AssembleRHS(nelem,dof_elem,node_elem,dof_node,elf,RHS)
   ENDDO

   DEALLOCATE(Nuu0)
   DEALLOCATE(Nuuu)
   DEALLOCATE(Nrr0)
   DEALLOCATE(Nrrr)
   DEALLOCATE(Nvvv)
   DEALLOCATE(Naaa)
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
            IF(ALLOCATED(elf)) DEALLOCATE(elf)
            IF(ALLOCATED(elk)) DEALLOCATE(elk)
            IF(ALLOCATED(elm)) DEALLOCATE(elm)
            IF(ALLOCATED(elg)) DEALLOCATE(elg)
        ENDIF


   END SUBROUTINE BD_GenerateDynamicElement
