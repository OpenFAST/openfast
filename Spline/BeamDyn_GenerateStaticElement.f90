   SUBROUTINE BeamDyn_GenerateStaticElement(uuN0,uuNf,Mass0,Stif0,gravity,u,&
                                            elem_total,node_elem,dof_node,ngp,StifK,RHS)

   REAL(ReKi),        INTENT(IN   ):: uuN0(:,:)
   REAL(ReKi),        INTENT(IN   ):: uuNf(:)
   REAL(ReKi),        INTENT(IN   ):: Mass0(:,:,:)
   REAL(ReKi),        INTENT(IN   ):: Stif0(:,:,:)
   REAL(ReKi),        INTENT(IN   ):: gravity(:)
   TYPE(BD_InputType),INTENT(IN   ):: u 
   INTEGER(IntKi),    INTENT(IN   ):: elem_total
   INTEGER(IntKi),    INTENT(IN   ):: node_elem
   INTEGER(IntKi),    INTENT(IN   ):: dof_node
   INTEGER(IntKi),    INTENT(IN   ):: ngp
   REAL(ReKi),        INTENT(  OUT):: StifK(:,:)
   REAL(ReKi),        INTENT(  OUT):: RHS(:) 

   REAL(ReKi),          ALLOCATABLE:: Nuuu(:)
   REAL(ReKi),          ALLOCATABLE:: Nrr0(:)
   REAL(ReKi),          ALLOCATABLE:: Nrrr(:)
   REAL(ReKi),          ALLOCATABLE:: elk(:,:)
   REAL(ReKi),          ALLOCATABLE:: elf(:)
   REAL(ReKi)                      :: DistrLoad_GL(dof_node,ngp)
   INTEGER(IntKi)                  :: dof_elem
   INTEGER(IntKi)                  :: rot_elem
   INTEGER(IntKi)                  :: nelem
   INTEGER(IntKi)                  :: j
   INTEGER(IntKi)                  :: temp_id
   INTEGER(IntKi)                  :: allo_stat

   dof_elem = dof_node * node_elem
   rot_elem = (dof_node/2) * node_elem


   ALLOCATE(Nuuu(dof_elem),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   Nuuu(:) = 0.0D0

   ALLOCATE(Nrr0(rot_elem),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   Nrr0(:) = 0.0D0

   ALLOCATE(Nrrr(rot_elem),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   Nrrr(:) = 0.0D0


   ALLOCATE(elf(dof_elem),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   elf(:) = 0.0D0

   ALLOCATE(elk(dof_elem,dof_elem),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   elk(:,:) = 0.0D0

   DO nelem=1,elem_total
       CALL ElemNodalDispGL(uuNf,node_elem,dof_node,nelem,Nuuu)
       CALL NodalRelRotGL(uuN0(:,nelem),node_elem,dof_node,Nrr0)
       CALL NodalRelRotGL(Nuuu,node_elem,dof_node,Nrrr)
       
       elk = 0.0D0
       elf = 0.0D0
       temp_id = (nelem-1)*ngp
       DO j=1,ngp
           DistrLoad_GL(1:3,j)  = u%DistrLoad%Force(1:3,temp_id+j+1)
           DistrLoad_GL(4:6,j)  = u%DistrLoad%Moment(1:3,temp_id+j+1)
       ENDDO
       CALL BeamDyn_StaticElementMatrix(uuN0(:,nelem),Nuuu,Nrr0,Nrrr,DistrLoad_GL,gravity,&
                                        Mass0(:,:,temp_id+1:temp_id+ngp),Stif0(:,:,temp_id+1:temp_id+ngp),&
                                        ngp,node_elem,dof_node,elk,elf)

       CALL AssembleStiffKGL(nelem,node_elem,dof_elem,dof_node,elk,StifK)
       CALL AssembleRHSGL(nelem,dof_elem,node_elem,dof_node,elf,RHS)
   ENDDO

   DEALLOCATE(Nuuu)
   DEALLOCATE(Nrr0)
   DEALLOCATE(Nrrr)
   DEALLOCATE(elf)
   DEALLOCATE(elk)

   9999 IF(allo_stat/=0) THEN
            IF(ALLOCATED(Nuuu)) DEALLOCATE(Nuuu)
            IF(ALLOCATED(Nrr0)) DEALLOCATE(Nrr0)
            IF(ALLOCATED(Nrrr)) DEALLOCATE(Nrrr)
            IF(ALLOCATED(elf)) DEALLOCATE(elf)
            IF(ALLOCATED(elk)) DEALLOCATE(elk)
        ENDIF


   END SUBROUTINE BeamDyn_GenerateStaticElement
