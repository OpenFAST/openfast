   SUBROUTINE ComputeRootForceNodal(uuN0,uuNf,Stif0,node_elem,dof_node,Fc)

   REAL(ReKi),    INTENT(IN   ):: uuN0(:)
   REAL(ReKi),    INTENT(IN   ):: uuNf(:)
   REAL(ReKi),    INTENT(IN   ):: Stif0(:,:,:)
   INTEGER(IntKi),INTENT(IN   ):: node_elem
   INTEGER(IntKi),INTENT(IN   ):: dof_node
   REAL(ReKi),    INTENT(  OUT):: Fc(:)
   ! Local variables 
   REAL(ReKi),ALLOCATABLE:: Nuu0(:)
   REAL(ReKi),ALLOCATABLE:: Nuuu(:)
   REAL(ReKi),ALLOCATABLE:: Nrr0(:)
   REAL(ReKi),ALLOCATABLE:: Nrrr(:)
   REAL(ReKi),ALLOCATABLE:: GLL_temp(:)
   REAL(ReKi),ALLOCATABLE:: w_temp(:)
   REAL(ReKi),ALLOCATABLE:: hhp(:,:)
   REAL(ReKi)::             E10(3)
   REAL(ReKi)::             E1(3)
   REAL(ReKi)::             RR0(3,3)
   REAL(ReKi)::             kapa(3)
   REAL(ReKi)::             Fd(6)
   REAL(ReKi)::             Oe(6,6)
   REAL(ReKi)::             Pe(6,6)
   REAL(ReKi)::             Qe(6,6)
   REAL(ReKi)::             Stif(6,6)
   REAL(ReKi)::             cet
   REAL(ReKi)::             Jac
   INTEGER(IntKi)::         i
   INTEGER(IntKi)::         dof_elem
   INTEGER(IntKi)::         rot_elem 
   INTEGER(IntKi)::         allo_stat


   dof_elem = dof_node*node_elem
   rot_elem = dof_elem/2

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

   ALLOCATE(GLL_temp(node_elem),STAT = allo_stat) 
   IF(allo_stat/=0) GOTO 9999
   GLL_temp = 0.0D0

   ALLOCATE(w_temp(node_elem),STAT = allo_stat) 
   IF(allo_stat/=0) GOTO 9999
   w_temp = 0.0D0

   ALLOCATE(hhp(node_elem,node_elem),STAT = allo_stat) 
   IF(allo_stat/=0) GOTO 9999
   hhp = 0.0D0

   i = 1
   ! Get Nodal Displacement Vector Nuu0 from Global Vector uuN0 at t=0    
   CALL ElemNodalDispGL(uuN0,node_elem,dof_node,i,Nuu0)
   ! Get Nodal Displacement Vector Nuuu from Global Vector uuNf
   CALL ElemNodalDispGL(uuNf,node_elem,dof_node,i,Nuuu)
   ! Compute Nodal Relative Rotation Vector Nrr0 
   CALL NodalRelRotGL(Nuu0,node_elem,dof_node,Nrr0)
   ! Compute Nodal Relative Rotation Vector Nrrr
   CALL NodalRelRotGL(Nuuu,node_elem,dof_node,Nrrr) 
   E10 = 0.0D0
   E1 = 0.0D0
   RR0 = 0.0D0
   kapa = 0.0D0
   Fc = 0.0D0
   Fd = 0.0D0
   Oe = 0.0D0
   Pe = 0.0D0
   Qe = 0.0D0
   Stif = 0.0D0
   cet = 0.0D0
   CALL BDyn_gen_gll_LSGL(node_elem-1,GLL_temp,w_temp)
   CALL BDyn_gen_deriv_LSGL(node_elem-1,GLL_temp,hhp)
   CALL NewNodalDataAt0(node_elem,dof_node,i,hhp,Nuu0,Jac,E10)
   CALL NodalData(Nuuu,Nrrr,Nuu0,Nrr0,E10,hhp,Stif0,Jac,&
                  &node_elem,1,dof_node,&
                  &E1,RR0,kapa,Stif,cet)
   CALL ElasticForce(E1,RR0,kapa,Stif,cet,Fc,Fd,Oe,Pe,Qe)

   DEALLOCATE(Nuu0)
   DEALLOCATE(Nuuu)
   DEALLOCATE(Nrr0)
   DEALLOCATE(Nrrr)
   DEALLOCATE(GLL_temp)
   DEALLOCATE(w_temp)
   DEALLOCATE(hhp)

   9999 IF(allo_stat/=0) THEN
            IF(ALLOCATED(Nuu0)) DEALLOCATE(Nuu0)
            IF(ALLOCATED(Nuuu)) DEALLOCATE(Nuuu)
            IF(ALLOCATED(Nrr0)) DEALLOCATE(Nrr0)
            IF(ALLOCATED(Nrrr)) DEALLOCATE(Nrrr)
            IF(ALLOCATED(GLL_temp)) DEALLOCATE(GLL_temp)
            IF(ALLOCATED(w_temp)) DEALLOCATE(w_temp)
            IF(ALLOCATED(hhp)) DEALLOCATE(hhp)
        ENDIF


   END SUBROUTINE ComputeRootForceNodal
