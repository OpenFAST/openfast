   SUBROUTINE ComputeRootForce(uuN0,uuNf,Stif0,node_elem,dof_node,ngp,RootForce)
   
   REAL(ReKi),INTENT(IN):: uuN0(:,:),uuNf(:),Stif0(:,:,:)
   INTEGER(IntKi),INTENT(IN):: node_elem,dof_node,ngp
   REAL(ReKi),INTENT(OUT):: RootForce(:)
   
   REAL(ReKi),ALLOCATABLE:: Nuu0(:),Nuuu(:),Nrr0(:),Nrrr(:)
   REAL(ReKi),ALLOCATABLE:: NStif0(:,:,:)
   REAL(ReKi),ALLOCATABLE:: eeeGP(:),fffGP(:)
   REAL(ReKi),ALLOCATABLE:: GLL_temp(:),w_temp(:),gp(:),gw(:)
   REAL(ReKi),ALLOCATABLE:: hhx(:),hpx(:)
   
   REAL(ReKi):: rr(1)
   INTEGER(IntKi):: dof_elem,rot_elem
   INTEGER(IntKi):: i,j,temp_id
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
   
   ALLOCATE(NStif0(dof_node,dof_node,ngp),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   NStif0 = 0.0D0
   
   ALLOCATE(gp(ngp), STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   gp = 0.0D0

   ALLOCATE(gw(ngp), STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   gw = 0.0D0

   ALLOCATE(hhx(ngp), STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   hhx = 0.0D0

   ALLOCATE(hpx(ngp), STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   hpx = 0.0D0

   ALLOCATE(GLL_temp(ngp), STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   GLL_temp = 0.0D0

   ALLOCATE(w_temp(ngp), STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   w_temp = 0.0D0
   
   ALLOCATE(eeeGP(ngp*dof_node), STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   eeeGP = 0.0D0
   
   ALLOCATE(fffGP(ngp*dof_node), STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   fffGP = 0.0D0
   
   i = 1
   Nuu0(:) = uuN0(:,i)
!   CALL ElemNodalDispGL(uuN0,node_elem,dof_node,1,Nuu0)
   CALL ElemNodalDispGL(uuNf,node_elem,dof_node,i,Nuuu)
!   CALL ElemNodalStifGL(Stif0,node_elem,dof_node,1,NStif0)
   DO i=1,ngp
       NStif0(:,:,i) = Stif0(:,:,i)
   ENDDO
   CALL NodalRelRotGL(Nuu0,node_elem,dof_node,Nrr0)
   CALL NodalRelRotGL(Nuuu,node_elem,dof_node,Nrrr)
   
   CALL ComputeStrainForceGP(Nuu0,Nuuu,Nrr0,Nrrr,NStif0,ngp,node_elem,dof_node,eeeGP,fffGP)
   
   CALL BeamDyn_gen_gll_LSGL(ngp-1,GLL_temp,w_temp)
   CALL BldGaussPointWeight(ngp,gp,gw)
   
   rr(1) = -1.0D0/gp(ngp)
   
   CALL diffmtc(ngp-1,1,rr,GLL_temp,1,hhx,hpx)

   RootForce(:) = 0.0D0
   DO i=1,ngp
       temp_id = (i-1)*dof_node
       DO j=1,dof_node
           RootForce(j) = RootForce(j) + hhx(i)*fffGP(temp_id+j)
       ENDDO
   ENDDO
   
   DEALLOCATE(gp)
   DEALLOCATE(gw)
   DEALLOCATE(hhx)
   DEALLOCATE(hpx)
   DEALLOCATE(GLL_temp)
   DEALLOCATE(w_temp)
   DEALLOCATE(Nuu0)
   DEALLOCATE(Nuuu)
   DEALLOCATE(Nrr0)
   DEALLOCATE(Nrrr)
   DEALLOCATE(NStif0)
   DEALLOCATE(eeeGP)
   DEALLOCATE(fffGP)

   9999 IF(allo_stat/=0) THEN
            IF(ALLOCATED(gp))  DEALLOCATE(gp)
            IF(ALLOCATED(gw))  DEALLOCATE(gw)
            IF(ALLOCATED(hhx)) DEALLOCATE(hhx)
            IF(ALLOCATED(hpx)) DEALLOCATE(hpx)
            IF(ALLOCATED(GLL_temp)) DEALLOCATE(GLL_temp)
            IF(ALLOCATED(w_temp)) DEALLOCATE(w_temp)
            IF(ALLOCATED(Nuu0)) DEALLOCATE(Nuu0)
            IF(ALLOCATED(Nuuu)) DEALLOCATE(Nuuu)
            IF(ALLOCATED(Nrr0)) DEALLOCATE(Nrr0)
            IF(ALLOCATED(Nrrr)) DEALLOCATE(Nrrr)
            IF(ALLOCATED(NStif0)) DEALLOCATE(NStif0)
            IF(ALLOCATED(eeeGP)) DEALLOCATE(eeeGP)
            IF(ALLOCATED(fffGP)) DEALLOCATE(fffGP)
        ENDIF
   
   END SUBROUTINE ComputeRootForce
