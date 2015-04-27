   SUBROUTINE BD_ElementMatrixGA2(Nuu0,Nuuu,Nrr0,Nrrr,Nvvv,Naaa,           &
                                  EStif0_GL,EMass0_GL,gravity,DistrLoad_GL,&
                                  damp_flag,beta,                          &
                                  ngp,node_elem,dof_node,elk,elf,elm,elg)

   REAL(ReKi),     INTENT(IN   ):: Nuu0(:)
   REAL(ReKi),     INTENT(IN   ):: Nuuu(:)
   REAL(ReKi),     INTENT(IN   ):: Nrr0(:)
   REAL(ReKi),     INTENT(IN   ):: Nrrr(:)
   REAL(ReKi),     INTENT(IN   ):: Nvvv(:)
   REAL(ReKi),     INTENT(IN   ):: Naaa(:)
   REAL(ReKi),     INTENT(IN   ):: EStif0_GL(:,:,:)
   REAL(ReKi),     INTENT(IN   ):: EMass0_GL(:,:,:)
   REAL(ReKi),     INTENT(IN   ):: gravity(:)
   REAL(ReKi),     INTENT(IN   ):: DistrLoad_GL(:,:)
   INTEGER(IntKi), INTENT(IN   ):: damp_flag
   REAL(ReKi),     INTENT(IN   ):: beta(:)
   INTEGER(IntKi), INTENT(IN   ):: ngp
   INTEGER(IntKi), INTENT(IN   ):: node_elem
   INTEGER(IntKi), INTENT(IN   ):: dof_node
   REAL(ReKi),     INTENT(  OUT):: elk(:,:)
   REAL(ReKi),     INTENT(  OUT):: elf(:)
   REAL(ReKi),     INTENT(  OUT):: elm(:,:)
   REAL(ReKi),     INTENT(  OUT):: elg(:,:)

   REAL(ReKi),       ALLOCATABLE:: gp(:)
   REAL(ReKi),       ALLOCATABLE:: gw(:)
   REAL(ReKi),       ALLOCATABLE:: hhx(:)
   REAL(ReKi),       ALLOCATABLE:: hpx(:)
   REAL(ReKi),       ALLOCATABLE:: GLL_temp(:)
   REAL(ReKi),       ALLOCATABLE:: w_temp(:)
   REAL(ReKi)                   :: uu0(6)
   REAL(ReKi)                   :: E10(3)
   REAL(ReKi)                   :: RR0(3,3)
   REAL(ReKi)                   :: kapa(3)
   REAL(ReKi)                   :: E1(3)
   REAL(ReKi)                   :: Stif(6,6)
   REAL(ReKi)                   :: cet
   REAL(ReKi)                   :: uuu(6)
   REAL(ReKi)                   :: uup(3)
   REAL(ReKi)                   :: Jacobian
   REAL(ReKi)                   :: gpr
   REAL(ReKi)                   :: Fc(6)
   REAL(ReKi)                   :: Fd(6)
   REAL(ReKi)                   :: Fg(6)
   REAL(ReKi)                   :: Oe(6,6)
   REAL(ReKi)                   :: Pe(6,6)
   REAL(ReKi)                   :: Qe(6,6)
   REAL(ReKi)                   :: Sd(6,6)
   REAL(ReKi)                   :: Od(6,6)
   REAL(ReKi)                   :: Pd(6,6)
   REAL(ReKi)                   :: Qd(6,6)
   REAL(ReKi)                   :: betaC(6,6)
   REAL(ReKi)                   :: Gd(6,6)
   REAL(ReKi)                   :: Xd(6,6)
   REAL(ReKi)                   :: Yd(6,6)
   REAL(ReKi)                   :: vvv(6)
   REAL(ReKi)                   :: vvp(6)
   REAL(ReKi)                   :: aaa(6)
   REAL(ReKi)                   :: mmm
   REAL(ReKi)                   :: mEta(3)
   REAL(ReKi)                   :: rho(3,3)
   REAL(ReKi)                   :: Fi(6)
   REAL(ReKi)                   :: Mi(6,6)
   REAL(ReKi)                   :: Ki(6,6)
   REAL(ReKi)                   :: Gi(6,6)
   INTEGER(IntKi)               :: igp
   INTEGER(IntKi)               :: i
   INTEGER(IntKi)               :: j
   INTEGER(IntKi)               :: m
   INTEGER(IntKi)               :: n
   INTEGER(IntKi)               :: temp_id1
   INTEGER(IntKi)               :: temp_id2
   INTEGER(IntKi)               :: allo_stat


   ALLOCATE(gp(ngp), STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   gp = 0.0D0

   ALLOCATE(gw(ngp), STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   gw = 0.0D0

   ALLOCATE(hhx(node_elem), STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   hhx = 0.0D0

   ALLOCATE(hpx(node_elem), STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   hpx = 0.0D0
   
   ALLOCATE(GLL_temp(node_elem), STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   GLL_temp = 0.0D0
   
   ALLOCATE(w_temp(node_elem), STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   w_temp = 0.0D0

   elk(:,:) = 0.0D0
   elf(:)   = 0.0D0
   elg(:,:) = 0.0D0
   elm(:,:) = 0.0D0
   
   CALL BD_GenerateGLL(node_elem-1,GLL_temp,w_temp)
   CALL BD_GaussPointWeight(ngp,gp,gw)
   DO igp=1,ngp
       gpr = gp(igp)
       CALL BD_ComputeJacobian(gpr,Nuu0,node_elem,dof_node,gp,GLL_temp,ngp,igp,hhx,hpx,Jacobian)
       CALL BD_GaussPointDataAt0(hhx,hpx,Nuu0,Nrr0,node_elem,dof_node,uu0,E10)
       Stif(:,:) = 0.0D0
       Stif(1:6,1:6) = EStif0_GL(1:6,1:6,igp)
       CALL BD_GaussPointData(hhx,hpx,Nuuu,Nrrr,uu0,E10,node_elem,dof_node,uuu,uup,E1,RR0,kapa,Stif,cet)
       CALL BD_ElasticForceGA2(E1,RR0,kapa,Stif,cet,Fc,Fd,Oe,Pe,Qe)

       mmm  = 0.0D0
       mEta = 0.0D0
       rho  = 0.0D0
       mmm          =  EMass0_GL(1,1,igp)
       mEta(2)      = -EMass0_GL(1,6,igp)
       mEta(3)      =  EMass0_GL(1,5,igp)
       rho(1:3,1:3) =  EMass0_GL(4:6,4:6,igp)
       CALL BD_GaussPointDataMass(hhx,hpx,Nvvv,Naaa,RR0,node_elem,dof_node,vvv,aaa,vvp,mmm,mEta,rho)
       CALL BD_InertialForce(mmm,mEta,rho,vvv,aaa,Fi,Mi,Gi,Ki)
       IF(damp_flag .NE. 0) THEN
           CALL BD_DissipativeForce(beta,Stif,vvv,vvp,E1,Fc,Fd,Sd,Od,Pd,Qd,betaC,Gd,Xd,Yd)
       ENDIF
       CALL BD_GravityForce(mmm,mEta,gravity,Fg)
       Fd(:) = Fd(:) - Fg(:)

       DO i=1,node_elem
           DO j=1,node_elem
               DO m=1,dof_node
                   temp_id1 = (i-1)*dof_node+m
                   DO n=1,dof_node
                       temp_id2 = (j-1)*dof_node+n
                       elk(temp_id1,temp_id2) = elk(temp_id1,temp_id2) + hhx(i)*Qe(m,n)*hhx(j)*Jacobian*gw(igp)
                       elk(temp_id1,temp_id2) = elk(temp_id1,temp_id2) + hhx(i)*Pe(m,n)*hpx(j)*Jacobian*gw(igp)
                       elk(temp_id1,temp_id2) = elk(temp_id1,temp_id2) + hpx(i)*Oe(m,n)*hhx(j)*Jacobian*gw(igp)
                       elk(temp_id1,temp_id2) = elk(temp_id1,temp_id2) + hpx(i)*Stif(m,n)*hpx(j)*Jacobian*gw(igp)
                       elk(temp_id1,temp_id2) = elk(temp_id1,temp_id2) + hhx(i)*Ki(m,n)*hhx(j)*Jacobian*gw(igp)
                       elm(temp_id1,temp_id2) = elm(temp_id1,temp_id2) + hhx(i)*Mi(m,n)*hhx(j)*Jacobian*gw(igp)
                       elg(temp_id1,temp_id2) = elg(temp_id1,temp_id2) + hhx(i)*Gi(m,n)*hhx(j)*Jacobian*gw(igp)
                       IF(damp_flag .NE. 0) THEN
                           elk(temp_id1,temp_id2) = elk(temp_id1,temp_id2) + hhx(i)*Qd(m,n)*hhx(j)*Jacobian*gw(igp)
                           elk(temp_id1,temp_id2) = elk(temp_id1,temp_id2) + hhx(i)*Pd(m,n)*hpx(j)*Jacobian*gw(igp)
                           elk(temp_id1,temp_id2) = elk(temp_id1,temp_id2) + hpx(i)*Od(m,n)*hhx(j)*Jacobian*gw(igp)
                           elk(temp_id1,temp_id2) = elk(temp_id1,temp_id2) + hpx(i)*Sd(m,n)*hpx(j)*Jacobian*gw(igp)
                           elg(temp_id1,temp_id2) = elg(temp_id1,temp_id2) + hhx(i)*Xd(m,n)*hhx(j)*Jacobian*gw(igp)
                           elg(temp_id1,temp_id2) = elg(temp_id1,temp_id2) + hhx(i)*Yd(m,n)*hpx(j)*Jacobian*gw(igp)
                           elg(temp_id1,temp_id2) = elg(temp_id1,temp_id2) + hpx(i)*Gd(m,n)*hhx(j)*Jacobian*gw(igp)
                           elg(temp_id1,temp_id2) = elg(temp_id1,temp_id2) + hpx(i)*betaC(m,n)*hpx(j)*Jacobian*gw(igp)
                       ENDIF
                   ENDDO
               ENDDO
           ENDDO
       ENDDO 

       DO i=1,node_elem
           DO j=1,dof_node
               temp_id1 = (i-1) * dof_node+j
               elf(temp_id1) = elf(temp_id1) - hhx(i)*Fd(j)*Jacobian*gw(igp)
               elf(temp_id1) = elf(temp_id1) - hpx(i)*Fc(j)*Jacobian*gw(igp)
               elf(temp_id1) = elf(temp_id1) - hhx(i)*Fi(j)*Jacobian*gw(igp)
           ENDDO
       ENDDO

   ENDDO

   DEALLOCATE(gp)
   DEALLOCATE(gw)
   DEALLOCATE(hhx)
   DEALLOCATE(hpx)
   DEALLOCATE(GLL_temp)
   DEALLOCATE(w_temp)

   9999 IF(allo_stat/=0) THEN
            IF(ALLOCATED(gp))  DEALLOCATE(gp)
            IF(ALLOCATED(gw))  DEALLOCATE(gw)
            IF(ALLOCATED(hhx)) DEALLOCATE(hhx)
            IF(ALLOCATED(hpx)) DEALLOCATE(hpx)
            IF(ALLOCATED(GLL_temp)) DEALLOCATE(GLL_temp)
            IF(ALLOCATED(w_temp)) DEALLOCATE(w_temp)
        ENDIF

   END SUBROUTINE BD_ElementMatrixGA2
