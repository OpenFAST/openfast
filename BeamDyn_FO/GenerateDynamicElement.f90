   SUBROUTINE GenerateDynamicElement(uuN0,uuN,vvN,Stif0,m00,mEta0,rho0,&
                                    &elem_total,node_elem,dof_node,ngp,RHS,MassM)

   REAL(ReKi),INTENT(IN):: uuN0(:) ! Initial position vector
   REAL(ReKi),INTENT(IN):: uuN(:) ! Displacement of Mass 1: m
   REAL(ReKi),INTENT(IN):: vvN(:) ! Velocity of Mass 1: m/s
   REAL(ReKi),INTENT(IN):: Stif0(:,:,:) ! Element stiffness matrix
   REAL(ReKi),INTENT(IN):: m00(:) ! Mass of beam per unit span at each node
   REAL(ReKi),INTENT(IN):: mEta0(:,:) ! Sectional m\Eta_0 at each node
   REAL(ReKi),INTENT(IN):: rho0(:,:,:) ! Sectional tensor of inertia per unit span
   INTEGER(IntKi),INTENT(IN):: elem_total ! Total number of elements
   INTEGER(IntKi),INTENT(IN):: node_elem ! Node per element
   INTEGER(IntKi),INTENT(IN):: dof_node ! Degrees of freedom per node
   INTEGER(IntKi),INTENT(IN):: ngp ! Number of Gauss points
   REAL(ReKi),INTENT(OUT):: MassM(:,:) ! Mass matrix 
   REAL(ReKi),INTENT(OUT):: RHS(:) ! Right hand side of the equation Ax=B  

   REAL(ReKi),ALLOCATABLE:: Nuu0(:) ! Nodal initial position for each element
   REAL(ReKi),ALLOCATABLE:: Nuuu(:) ! Nodal displacement of Mass 1 for each element
   REAL(ReKi),ALLOCATABLE:: Nrr0(:) ! Nodal rotation parameters for initial position 
   REAL(ReKi),ALLOCATABLE:: Nrrr(:) ! Nodal rotation parameters for displacement of Mass 1
   REAL(ReKi),ALLOCATABLE:: Nvvv(:) ! Nodal velocity of Mass 1: m/s for each element
   REAL(ReKi),ALLOCATABLE:: NStif0(:,:,:) ! Nodal material properties for each element
   REAL(ReKi),ALLOCATABLE:: Nm00(:) ! Nodal mass of beam per unit span for each element
   REAL(ReKi),ALLOCATABLE:: NmEta0(:,:) ! Nodal sectional m\Eta_0 for each element
   REAL(ReKi),ALLOCATABLE:: Nrho0(:,:,:) ! Nodal sectional tensor of inertia per unit span for each element
   REAL(ReKi),ALLOCATABLE:: elf(:) !
   REAL(ReKi),ALLOCATABLE:: elm(:,:) !

   INTEGER(IntKi):: dof_elem ! Degree of freedom per node
   INTEGER(IntKi):: rot_elem ! Rotational degrees of freedom
   INTEGER(IntKi):: nelem ! number of elements
   INTEGER(IntKi):: j ! Index counter
   INTEGER(IntKi):: allo_stat ! Allows for an error code return

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
