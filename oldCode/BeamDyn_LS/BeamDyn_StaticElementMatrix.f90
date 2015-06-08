   SUBROUTINE BeamDyn_StaticElementMatrix(Nuu0,Nuuu,Nrr0,Nrrr,Distr_GL,gravity,&
                                          EMass0_GL,EStif0_GL,ngp,node_elem,dof_node,elk,elf)

   REAL(ReKi),    INTENT(IN   ):: Nuu0(:)
   REAL(ReKi),    INTENT(IN   ):: Nuuu(:)
   REAL(ReKi),    INTENT(IN   ):: Nrr0(:)
   REAL(ReKi),    INTENT(IN   ):: Nrrr(:)
   REAL(ReKi),    INTENT(IN   ):: Distr_GL(:,:)
   REAL(ReKi),    INTENT(IN   ):: gravity(:)
   REAL(ReKi),    INTENT(IN   ):: EMass0_GL(:,:,:)
   REAL(ReKi),    INTENT(IN   ):: EStif0_GL(:,:,:)
   INTEGER(IntKi),INTENT(IN   ):: ngp
   INTEGER(IntKi),INTENT(IN   ):: node_elem
   INTEGER(IntKi),INTENT(IN   ):: dof_node
   REAL(ReKi),    INTENT(  OUT):: elk(:,:)
   REAL(ReKi),    INTENT(  OUT):: elf(:)      

   REAL(ReKi)                  :: gp(ngp)
   REAL(ReKi)                  :: gw(ngp)
   REAL(ReKi)                  :: hhx(node_elem)
   REAL(ReKi)                  :: hpx(node_elem)
   REAL(ReKi)                  :: GLL_temp(node_elem)
   REAL(ReKi)                  :: w_temp(node_elem)
   REAL(ReKi)                  :: uu0(6)
   REAL(ReKi)                  :: E10(3)
   REAL(ReKi)                  :: RR0(3,3)
   REAL(ReKi)                  :: kapa(3)
   REAL(ReKi)                  :: E1(3)
   REAL(ReKi)                  :: Stif(6,6)
   REAL(ReKi)                  :: cet
   REAL(ReKi)                  :: mmm
   REAL(ReKi)                  :: mEta(3)
   REAL(ReKi)                  :: rho(3,3)
   REAL(ReKi)                  :: uuu(6)
   REAL(ReKi)                  :: uup(3)
   REAL(ReKi)                  :: vvv(6)
   REAL(ReKi)                  :: Jacobian
   REAL(ReKi)                  :: gpr
   REAL(ReKi)                  :: Fc(6)
   REAL(ReKi)                  :: Fd(6)
   REAL(ReKi)                  :: Fg(6)
   REAL(ReKi)                  :: Oe(6,6)
   REAL(ReKi)                  :: Pe(6,6)
   REAL(ReKi)                  :: Qe(6,6)
   REAL(ReKi)                  :: Nvvv(dof_node*node_elem)
   INTEGER(IntKi)              :: igp
   INTEGER(IntKi)              :: i
   INTEGER(IntKi)              :: j
   INTEGER(IntKi)              :: m
   INTEGER(IntKi)              :: n
   INTEGER(IntKi)              :: temp_id1
   INTEGER(IntKi)              :: temp_id2
   INTEGER(IntKi)              :: allo_stat

   vvv(:)   = 0.0D0
   Nvvv(:)  = 0.0D0
   elk(:,:) = 0.0D0
   elf(:)   = 0.0D0
   CALL BD_gen_gll_LSGL(node_elem-1,GLL_temp,w_temp)
   CALL BldGaussPointWeight(ngp,gp,gw)
   DO igp=1,ngp
       gpr = gp(igp)
       CALL BldComputeJacobianLSGL(gpr,Nuu0,node_elem,dof_node,gp,GLL_temp,ngp,igp,hhx,hpx,Jacobian)
       CALL BldGaussPointDataAt0(hhx,hpx,Nuu0,Nrr0,node_elem,dof_node,uu0,E10)
       Stif(1:6,1:6) = EStif0_GL(1:6,1:6,igp)
       CALL BldGaussPointData(hhx,hpx,Nuuu,Nrrr,uu0,E10,node_elem,dof_node,uuu,uup,E1,RR0,kapa,Stif,cet)
       CALL BeamDyn_StaticElasticForce(E1,RR0,kapa,Stif,cet,Fc,Fd,Oe,Pe,Qe)
       mmm          = EMass0_GL(1,1,igp)
       mEta(2)      =-EMass0_GL(1,6,igp)
       mEta(3)      = EMass0_GL(1,5,igp)
       rho(1:3,1:3) = EMass0_GL(4:6,4:6,igp)
       CALL BldGaussPointDataMass(hhx,hpx,Nvvv,RR0,node_elem,dof_node,vvv,mmm,mEta,rho)
       CALL GravityLoads(mmm,mEta,gravity,Fg)
       Fd(:) = Fd(:) - Fg(:) - Distr_GL(:,igp)

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
                   ENDDO
               ENDDO
           ENDDO
       ENDDO 

   
       DO i=1,node_elem
           DO j=1,dof_node
               temp_id1 = (i-1) * dof_node+j
               elf(temp_id1) = elf(temp_id1) - hhx(i)*Fd(j)*Jacobian*gw(igp)
               elf(temp_id1) = elf(temp_id1) - hpx(i)*Fc(j)*Jacobian*gw(igp)
           ENDDO
       ENDDO

   ENDDO



   END SUBROUTINE BeamDyn_StaticElementMatrix
