   SUBROUTINE GenerateDynamicElement_AM2(uuN0,uuN,vvN,uuN00,vvN00,Stif0,Mass0,gravity,u,u0,&
                                        &damp_flag,beta,&
                                        &elem_total,node_elem,dof_node,ngp,dt,RHS,MassM)
   !----------------------------------------------------------------------------------------
   ! This subroutine computes Global mass matrix and force vector for the beam.
   !----------------------------------------------------------------------------------------
   REAL(ReKi),        INTENT(IN   ):: uuN0(:,:) ! Initial position vector
   REAL(ReKi),        INTENT(IN   ):: uuN(:) ! Displacement of Mass 1: m
   REAL(ReKi),        INTENT(IN   ):: vvN(:) ! Velocity of Mass 1: m/s
   REAL(ReKi),        INTENT(IN   ):: uuN00(:) ! Displacement of Mass 1: m
   REAL(ReKi),        INTENT(IN   ):: vvN00(:) ! Velocity of Mass 1: m/s
   REAL(ReKi),        INTENT(IN   ):: Stif0(:,:,:) ! Element stiffness matrix
   REAL(ReKi),        INTENT(IN   ):: Mass0(:,:,:) ! Element stiffness matrix
   REAL(ReKi),        INTENT(IN   ):: gravity(:) ! Velocity of Mass 1: m/s
   REAL(DbKi),        INTENT(IN   ):: dt       ! Velocity of Mass 1: m/s
   TYPE(BD_InputType),INTENT(IN   ):: u           ! Inputs at t
   TYPE(BD_InputType),INTENT(IN   ):: u0          ! Inputs at t
   INTEGER(IntKi),    INTENT(IN   ):: damp_flag ! Total number of elements
   REAL(ReKi),        INTENT(IN   ):: beta
   INTEGER(IntKi),    INTENT(IN   ):: elem_total ! Total number of elements
   INTEGER(IntKi),    INTENT(IN   ):: node_elem ! Node per element
   INTEGER(IntKi),    INTENT(IN   ):: dof_node ! Degrees of freedom per node
   INTEGER(IntKi),    INTENT(IN   ):: ngp ! Number of Gauss points
   REAL(ReKi),        INTENT(INOUT):: MassM(:,:) ! Mass matrix 
   REAL(ReKi),        INTENT(INOUT):: RHS(:) ! Right hand side of the equation Ax=B  

!   REAL(ReKi),        ALLOCATABLE  :: Nuu0(:) ! Nodal initial position for each element
!   REAL(ReKi),        ALLOCATABLE  :: Nuuu(:) ! Nodal displacement of Mass 1 for each element
!   REAL(ReKi),        ALLOCATABLE  :: Nuuu0(:) ! Nodal displacement of Mass 1 for each element
!   REAL(ReKi),        ALLOCATABLE  :: Nrr0(:) ! Nodal rotation parameters for initial position 
!   REAL(ReKi),        ALLOCATABLE  :: Nrrr(:) ! Nodal rotation parameters for displacement of Mass 1
!   REAL(ReKi),        ALLOCATABLE  :: Nrrr0(:) ! Nodal rotation parameters for displacement of Mass 1
!   REAL(ReKi),        ALLOCATABLE  :: Nvvv(:) ! Nodal velocity of Mass 1: m/s for each element
!   REAL(ReKi),        ALLOCATABLE  :: Nvvv0(:) ! Nodal velocity of Mass 1: m/s for each element
!   REAL(ReKi),        ALLOCATABLE  :: EStif0_GL(:,:,:) ! Nodal material properties for each element
!   REAL(ReKi),        ALLOCATABLE  :: EMass0_GL(:,:,:) ! Nodal material properties for each element
!   REAL(ReKi),        ALLOCATABLE  :: DistrLoad_GL(:,:) ! Nodal material properties for each element
!   REAL(ReKi),        ALLOCATABLE  :: DistrLoad_GL0(:,:) ! Nodal material properties for each element
!   REAL(ReKi),        ALLOCATABLE  :: elf1(:) ! Total element force (Fc, Fd, Fb)
!   REAL(ReKi),        ALLOCATABLE  :: elf2(:) ! Total element force (Fc, Fd, Fb)
!   REAL(ReKi),        ALLOCATABLE  :: elm11(:,:) ! Element mass matrix
!   REAL(ReKi),        ALLOCATABLE  :: elm12(:,:) ! Element mass matrix
!   REAL(ReKi),        ALLOCATABLE  :: elm21(:,:) ! Element mass matrix
!   REAL(ReKi),        ALLOCATABLE  :: elm22(:,:) ! Element mass matrix

   REAL(ReKi) :: Nuu0(dof_node*node_elem) ! Nodal initial position for each element
   REAL(ReKi) :: Nuuu(dof_node*node_elem) ! Nodal displacement of Mass 1 for each element
   REAL(ReKi) :: Nuuu0(dof_node*node_elem) ! Nodal displacement of Mass 1 for each element
   REAL(ReKi) :: Nrr0(3*node_elem) ! Nodal rotation parameters for initial position 
   REAL(ReKi) :: Nrrr(3*node_elem) ! Nodal rotation parameters for displacement of Mass 1
   REAL(ReKi) :: Nrrr0(3*node_elem) ! Nodal rotation parameters for displacement of Mass 1
   REAL(ReKi) :: Nvvv(dof_node*node_elem) ! Nodal velocity of Mass 1: m/s for each element
   REAL(ReKi) :: Nvvv0(dof_node*node_elem) ! Nodal velocity of Mass 1: m/s for each element
   REAL(ReKi) :: EStif0_GL(6,6,node_elem-1) ! Nodal material properties for each element
   REAL(ReKi) :: EMass0_GL(6,6,node_elem-1) ! Nodal material properties for each element
   REAL(ReKi) :: DistrLoad_GL(6,node_elem-1) ! Nodal material properties for each element
   REAL(ReKi) :: DistrLoad_GL0(6,node_elem-1) ! Nodal material properties for each element
   REAL(ReKi) :: elf1(dof_node*node_elem) ! Total element force (Fc, Fd, Fb)
   REAL(ReKi) :: elf2(dof_node*node_elem) ! Total element force (Fc, Fd, Fb)
   REAL(ReKi) :: elm11(dof_node*node_elem,dof_node*node_elem) ! Element mass matrix
   REAL(ReKi) :: elm12(dof_node*node_elem,dof_node*node_elem) ! Element mass matrix
   REAL(ReKi) :: elm21(dof_node*node_elem,dof_node*node_elem) ! Element mass matrix
   REAL(ReKi) :: elm22(dof_node*node_elem,dof_node*node_elem) ! Element mass matrix

   INTEGER(IntKi)                  :: dof_elem ! Degree of freedom per node
   INTEGER(IntKi)                  :: rot_elem ! Rotational degrees of freedom
   INTEGER(IntKi)                  :: nelem ! number of elements
   INTEGER(IntKi)                  :: j ! Index counter
   INTEGER(IntKi)                  :: temp_id ! Index counter
   INTEGER(IntKi)                  :: allo_stat ! Allows for an error code return

   dof_elem = dof_node * node_elem
   rot_elem = (dof_node/2) * node_elem

!   ALLOCATE(Nuu0(dof_elem),STAT = allo_stat)
!   IF(allo_stat/=0) GOTO 9999
!   Nuu0 = 0.0D0
!
!   ALLOCATE(Nuuu(dof_elem),STAT = allo_stat)
!   IF(allo_stat/=0) GOTO 9999
!   Nuuu = 0.0D0
!
!   ALLOCATE(Nuuu0(dof_elem),STAT = allo_stat)
!   IF(allo_stat/=0) GOTO 9999
!   Nuuu0 = 0.0D0
!
!   ALLOCATE(Nrr0(rot_elem),STAT = allo_stat)
!   IF(allo_stat/=0) GOTO 9999
!   Nrr0 = 0.0D0
!
!   ALLOCATE(Nrrr(rot_elem),STAT = allo_stat)
!   IF(allo_stat/=0) GOTO 9999
!   Nrrr = 0.0D0
!
!   ALLOCATE(Nrrr0(rot_elem),STAT = allo_stat)
!   IF(allo_stat/=0) GOTO 9999
!   Nrrr0 = 0.0D0
!
!   ALLOCATE(Nvvv(dof_elem),STAT = allo_stat)
!   IF(allo_stat/=0) GOTO 9999
!   Nvvv = 0.0D0
!
!   ALLOCATE(Nvvv0(dof_elem),STAT = allo_stat)
!   IF(allo_stat/=0) GOTO 9999
!   Nvvv0 = 0.0D0
!
!   ALLOCATE(EStif0_GL(dof_node,dof_node,node_elem-1),STAT = allo_stat)
!   IF(allo_stat/=0) GOTO 9999
!   EStif0_GL = 0.0D0
!
!   ALLOCATE(EMass0_GL(dof_node,dof_node,node_elem-1),STAT = allo_stat)
!   IF(allo_stat/=0) GOTO 9999
!   EMass0_GL = 0.0D0
!
!   ALLOCATE(DistrLoad_GL(dof_node,node_elem-1),STAT = allo_stat)
!   IF(allo_stat/=0) GOTO 9999
!   DistrLoad_GL = 0.0D0
!
!   ALLOCATE(DistrLoad_GL0(dof_node,node_elem-1),STAT = allo_stat)
!   IF(allo_stat/=0) GOTO 9999
!   DistrLoad_GL0 = 0.0D0
!
!   ALLOCATE(elf1(dof_elem),STAT = allo_stat)
!   IF(allo_stat/=0) GOTO 9999
!   elf1 = 0.0D0
!
!   ALLOCATE(elf2(dof_elem),STAT = allo_stat)
!   IF(allo_stat/=0) GOTO 9999
!   elf2 = 0.0D0
!
!   ALLOCATE(elm11(dof_elem,dof_elem),STAT = allo_stat)
!   IF(allo_stat/=0) GOTO 9999
!   elm11 = 0.0D0
!
!   ALLOCATE(elm12(dof_elem,dof_elem),STAT = allo_stat)
!   IF(allo_stat/=0) GOTO 9999
!   elm12 = 0.0D0
!
!   ALLOCATE(elm21(dof_elem,dof_elem),STAT = allo_stat)
!   IF(allo_stat/=0) GOTO 9999
!   elm21 = 0.0D0
!
!   ALLOCATE(elm22(dof_elem,dof_elem),STAT = allo_stat)
!   IF(allo_stat/=0) GOTO 9999
!   elm22 = 0.0D0

   DO nelem=1,elem_total
       Nuu0(:) = uuN0(:,nelem)
       CALL ElemNodalDispGL(uuN,node_elem,dof_node,nelem,Nuuu)
       CALL ElemNodalDispGL(uuN00,node_elem,dof_node,nelem,Nuuu0)

       temp_id = (nelem-1)*ngp
       DO j=1,ngp
           EStif0_GL(1:6,1:6,j) = Stif0(1:6,1:6,temp_id+j)
           EMass0_GL(1:6,1:6,j) = Mass0(1:6,1:6,temp_id+j)
           DistrLoad_GL(1:3,j)  = u%DistrLoad%Force(1:3,temp_id+j+1)
           DistrLoad_GL(4:6,j)  = u%DistrLoad%Moment(1:3,temp_id+j+1)
           DistrLoad_GL0(1:3,j)  = u0%DistrLoad%Force(1:3,temp_id+j+1)
           DistrLoad_GL0(4:6,j)  = u0%DistrLoad%Moment(1:3,temp_id+j+1)
       ENDDO
       
       CALL NodalRelRotGL(Nuu0,node_elem,dof_node,Nrr0)
       CALL NodalRelRotGL(Nuuu,node_elem,dof_node,Nrrr)
       CALL NodalRelRotGL(Nuuu0,node_elem,dof_node,Nrrr0)

       CALL ElemNodalDispGL(vvN,node_elem,dof_node,nelem,Nvvv)
       CALL ElemNodalDispGL(vvN00,node_elem,dof_node,nelem,Nvvv0)

       CALL ElementMatrix_AM2(Nuu0,Nuuu,Nuuu0,Nrr0,Nrrr,Nrrr0,Nvvv,Nvvv0,&
                             &EStif0_GL,EMass0_GL,gravity,DistrLoad_GL,DistrLoad_GL0,&
                             &ngp,dt,node_elem,dof_node,damp_flag,beta,&
                             &elf1,elf2,elm11,elm12,elm21,elm22)

!WRITE(*,*) "TEST3"

       CALL AssembleStiffK_AM2(nelem,node_elem,dof_elem,dof_node,&
                              &elm11,elm12,elm21,elm22,MassM)
       CALL AssembleRHS_AM2(nelem,dof_elem,node_elem,dof_node,elf1,elf2,RHS)

   ENDDO

!   DEALLOCATE(Nuu0)
!   DEALLOCATE(Nuuu)
!   DEALLOCATE(Nuuu0)
!   DEALLOCATE(Nrr0)
!   DEALLOCATE(Nrrr)
!   DEALLOCATE(Nrrr0)
!   DEALLOCATE(Nvvv)
!   DEALLOCATE(Nvvv0)
!   DEALLOCATE(EStif0_GL)
!   DEALLOCATE(EMass0_GL)
!   DEALLOCATE(DistrLoad_GL)
!   DEALLOCATE(DistrLoad_GL0)
!   DEALLOCATE(elf1)
!   DEALLOCATE(elf2)
!   DEALLOCATE(elm11)
!   DEALLOCATE(elm12)
!   DEALLOCATE(elm21)
!   DEALLOCATE(elm22)
!
!   9999 IF(allo_stat/=0) THEN
!            IF(ALLOCATED(Nuu0)) DEALLOCATE(Nuu0)
!            IF(ALLOCATED(Nuuu)) DEALLOCATE(Nuuu)
!            IF(ALLOCATED(Nuuu0)) DEALLOCATE(Nuuu0)
!            IF(ALLOCATED(Nrr0)) DEALLOCATE(Nrr0)
!            IF(ALLOCATED(Nrrr)) DEALLOCATE(Nrrr)
!            IF(ALLOCATED(Nrrr0)) DEALLOCATE(Nrrr0)
!            IF(ALLOCATED(Nvvv)) DEALLOCATE(Nvvv)
!            IF(ALLOCATED(Nvvv0)) DEALLOCATE(Nvvv0)
!            IF(ALLOCATED(EStif0_GL)) DEALLOCATE(EStif0_GL)
!            IF(ALLOCATED(EMass0_GL)) DEALLOCATE(EMass0_GL)
!            IF(ALLOCATED(DistrLoad_GL)) DEALLOCATE(DistrLoad_GL)
!            IF(ALLOCATED(DistrLoad_GL0)) DEALLOCATE(DistrLoad_GL0)
!            IF(ALLOCATED(elf1)) DEALLOCATE(elf1)
!            IF(ALLOCATED(elf2)) DEALLOCATE(elf2)
!            IF(ALLOCATED(elm11)) DEALLOCATE(elm11)
!            IF(ALLOCATED(elm12)) DEALLOCATE(elm12)
!            IF(ALLOCATED(elm21)) DEALLOCATE(elm21)
!            IF(ALLOCATED(elm22)) DEALLOCATE(elm22)
!        ENDIF


   END SUBROUTINE GenerateDynamicElement_AM2
